import math
import numpy as np
from matplotlib import pyplot
from scipy.optimize import curve_fit
import skrf as rf
from pylab import *
import subprocess
from scipy.optimize import minimize
import os
from reportlab.pdfgen.canvas import Canvas
from reportlab.pdfbase import pdfdoc
from scipy.interpolate import interp1d
import datetime



# Modify the circuit files to adjust the number of generated frequency points
def modifyCircuitFile(path):
    with open(path, 'r') as file:
        lines = file.readlines()

    index = 0
    for line in lines:
        if line.startswith('.sp'):
            lines[index] = '.sp dec 51 100 6G\n'
        
        index += 1

    with open(path, 'w') as file:
        file.writelines(lines)

    return None


# Generates the error based on the difference between the absolute value of the target and generated Z params
# The function is minimized later in the code as part of the optimization of the generated Z params
def s_param_error_RLC(x):
    R = x[0]
    L = x[1]
    C = x[2]

    params_path = ngspice + 'param.txt'
    with open(params_path, 'w') as paramf:
        paramf.write('.param C=%.20f L=%.20f R=%.20f' % (C, L, R))
    
    subprocess.call('ngspice.exe C.cir', shell=True)

    test_z21 = circuitCalculations()

    test_z21 = np.abs(test_z21)
    target_z = np.abs(z)
    error = np.abs(target_z - test_z21)

    return np.sum(error)


# Generates the error for the complex circuit model
def s_param_error_complex(x):
    Cnom = x[0]
    Rlfloss = x[1]
    Lmount = x[2]
    Rhfskin = x[3]
    Lhfskin = x[4]
    Rmfskin = x[5]
    Llfskin = x[6]
    Rlfskin = x[7]

    params_path = ngspice + 'paramcomplex.txt'
    with open(params_path, 'w') as paramf:
        paramf.write('.param Cnom=%.20f Rlfloss=%.20f Lmount=%.20f Rhfskin=%.20f Lhfskin=%.20f Rmfskin=%.20f Llfskin=%.20f Rlfskin=%.20f' % (Cnom, Rlfloss, Lmount, Rhfskin, Lhfskin, Rmfskin, Llfskin, Rlfskin))

    subprocess.call('ngspice.exe C_complex.cir', shell=True)

    test_z21 = circuitCalculations()

    test_z21 = np.abs(test_z21)
    target_z = np.abs(z)
    error = np.abs(target_z - test_z21)

    return np.sum(error)


# Given the parameters, generates graphs to be used in the PDF report as PNG images
# Defined as a function to avoid boilerplate due to large number of images
def plotGraph(noPlots, title, x, y, label, name, oPlot):
    pyplot.title(title)

    for i in range(noPlots):
        plotLines = '-'
        if oPlot == 1 and i % 2 == 0:
            plotLines = 'o'
        pyplot.plot(x[i], y[i], plotLines, label=label[i])

    pyplot.legend(loc='best')

    pyplot.xlabel('Frequency[Hz]')
    pyplot.ylabel('Impedance[Ohm]')
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.grid()

    pyplot.savefig(name, dpi=300)

    pyplot.clf()
    return None


# Given a path to the S params file, interpolates the values based on the frequencies and converts to absolute Z value, if the model is series
# Again, defined as a function to avoid boilerplate and clean the code up a bit
def getZSeries(path):
    ntwk = rf.Network(path)
    f = ntwk.f
    s = ntwk.s

    s11_interp = interp1d(f, s[:, 0, 0], kind='cubic')
    s12_interp = interp1d(f, s[:, 0, 1], kind='cubic')
    s21_interp = interp1d(f, s[:, 1, 0], kind='cubic')
    s22_interp = interp1d(f, s[:, 1, 1], kind='cubic')

    start_index = math.log10(f[0])
    stop_index = math.log10(f[-1])
    f = np.logspace(start_index, stop_index, 397, base=10)
    f[-1] -= 1

    s_params_interp = np.zeros((len(f), 2, 2), dtype=np.complex128)
    s_params_interp[:, 0, 0] = s11_interp(f) 
    s_params_interp[:, 0, 1] = s12_interp(f)
    s_params_interp[:, 1, 0] = s21_interp(f)
    s_params_interp[:, 1, 1] = s22_interp(f)

    y = rf.s2y(s=s_params_interp, z0=ntwk.z0[0,0])
    y21 = y[:,1,0]
    z = np.array([])

    for i in range(y21.size):
        z = np.append(z, 1 / y21[i])
    
    f_z = [f, z]
    return f_z


# Given a path to the S params file, interpolates the values based on the frequencies and converts to absolute Z value, if the model is shunt
# Again, defined as a function to avoid boilerplate and clean the code up a bit
def getZShunt(path):
    ntwk = rf.Network(path)
    f = ntwk.f
    s = ntwk.s

    s11_interp = interp1d(f, s[:, 0, 0], kind='cubic')
    s12_interp = interp1d(f, s[:, 0, 1], kind='cubic')
    s21_interp = interp1d(f, s[:, 1, 0], kind='cubic')
    s22_interp = interp1d(f, s[:, 1, 1], kind='cubic')

    start_index = math.log10(f[0])
    stop_index = math.log10(f[-1])
    f = np.logspace(start_index, stop_index, 397, base=10)
    f[-1] -= 1

    s_params_interp = np.zeros((len(f), 2, 2), dtype=np.complex128)
    s_params_interp[:, 0, 0] = s11_interp(f) 
    s_params_interp[:, 0, 1] = s12_interp(f)
    s_params_interp[:, 1, 0] = s21_interp(f)
    s_params_interp[:, 1, 1] = s22_interp(f)

    z = rf.s2z(s=s_params_interp, z0=ntwk.z0[0,0])
    z = z[:,1,0]
    
    f_z = [f, z]
    return f_z


# Defines the RLC values to be used as heuristics based on simple impedance formulas for series RLC circuit
# Not necessarily accurate, but the optimization takes care of that
def heuristics(f, z):
    z = np.abs(z)
    r = np.amin(z)

    low_f = float(f[0])
    low_f_z = float(z[0])
    z_c = math.sqrt(pow(low_f_z, 2) - pow(r, 2))

    c = 1 / (z_c * 2 * math.pi * low_f)

    resonant_index = np.where(z == r)
    res_f = float(f[resonant_index])
    resonant_w = res_f * 2 * math.pi

    l = 1 / (c * pow(resonant_w, 2))
    x = [r, l, c]
    return x


# Defines the values to be used as heuristics for the complex model
# Calculation of values is based on Figure 9 of document present in directory
def heuristics_complex(f, z):
    re_z = np.real(z)
    im_z = np.imag(z)

    low_f = float(f[0])
    z_c = math.sqrt(pow(float(im_z[0]), 2))
    c = 1 / (z_c * 2 * math.pi * low_f)

    r_lfloss = pow(float(im_z[0]), 2) / abs(re_z[0])
    r_hfskin = abs(float(re_z[-1]))
    r_mfskin = abs(r_hfskin) / 8

    high_f = float(f[-1])
    z_l = math.sqrt(pow(float(im_z[-1]), 2))
    l = z_l / (2 * math.pi * high_f)

    x = [l, c, r_lfloss, r_hfskin, r_mfskin, r_mfskin]
    return x


# Specifies the type of model, the path to ngspice, the path to the target S params and the complexity of the circuit
# Additional info in readme.txt
def configuration():
    file = open('config.txt')
    lines = file.readlines()
    model = ''
    ngspice = ''
    spath = ''
    circuit = ''
    for line in lines:
        line = line.lower()
        line = line.rstrip()
        if line.startswith('#'):
            continue
        else:
            splits = line.split('=')
            if splits[0] == 'model':
                model = splits[1]
            elif splits[0] == 'ngspice':
                ngspice = splits[1]
            elif splits[0] == 'spath':
                spath = splits[1]
            elif splits[0] == 'circuit':
                circuit = splits[1]

    x = [model, ngspice, spath, circuit]
    return x


# Calculates the Z params after the appropriate model (either series RLC or complex) has been loaded and filled with parameters
def circuitCalculations():
    ntwk = rf.Network('c.s2p')
    test_s = ntwk.s
    test_y = rf.s2y(s=test_s, z0=ntwk.z0[0,0])
    test_z = rf.y2z(y=test_y)
    test_z21 = test_z[:,1,0]

    y = np.array([])

    for i in range(test_z21.size):
        y = np.append(y, test_z21[i])

    return y



x = configuration()
model = x[0]
ngspice = x[1]
spath = x[2]
circuit = x[3]

if model == 'series':
    f_z = getZSeries(spath)
if model == 'shunt':
    f_z = getZShunt(spath)

f = f_z[0]
z = f_z[1]

if circuit == 'rlc':
    x = heuristics(f, z)
    r_heur = x[0]
    l_heur = x[1]
    c_heur = x[2]

if circuit == 'complex':
    x = heuristics_complex(f, z)
    Lmount_heur = x[0]
    Cnom_heur = x[1]
    Rlfloss_heur = x[2]
    Rhfskin_heur = x[3]
    Rmfskin_heur = x[4]
    Rlfskin_heur = x[5]


original_dir = os.getcwd()

os.chdir(ngspice)

# Calculate equivalent values for series RLC circuit
if circuit == 'rlc':
    modifyCircuitFile('C.cir')
    rlc_test = [r_heur, l_heur, c_heur]
    bounds = [(0.01 * r_heur, 100 * r_heur), (0.01 * l_heur, 100 * l_heur), (0.01 * c_heur, 100 * c_heur)]
    result = minimize(s_param_error_RLC, rlc_test, bounds=bounds, method='Nelder-Mead', options={'maxiter': 100})

    r_final = result.x[0]
    l_final = result.x[1]
    c_final = result.x[2]

    params_path = 'param.txt'
    with open(params_path, 'w') as paramf:
        paramf.write('.param C=%.20f L=%.20f R=%.20f' % (c_final, l_final, r_final))
        
    subprocess.call('ngspice.exe C.cir', shell=True)

    y = circuitCalculations()   

    with open(params_path, 'w') as paramf:
        paramf.write('.param C=%.20f L=%.20f R=%.20f' % (c_heur, l_heur, r_heur))

    subprocess.call('ngspice.exe C.cir', shell=True)

    y_heur = circuitCalculations()


# Calculate equivalent values for complex circuit
if circuit == 'complex':
    modifyCircuitFile('C_complex.cir')
    rlc_test = [Cnom_heur, Rlfloss_heur, Lmount_heur, Rhfskin_heur, Lmount_heur * 0.1, Rmfskin_heur, Lmount_heur * 0.1, Rlfskin_heur]
    bounds = [(Cnom_heur, Cnom_heur), (0.01 * Rlfloss_heur, 100 * Rlfloss_heur), (0.01 * Lmount_heur, 100 * Lmount_heur), (0.01 * Rhfskin_heur, 100 * Rhfskin_heur), (0.001 * Lmount_heur, 10 * Lmount_heur), (0.01 * Rmfskin_heur, 100 * Rmfskin_heur), (0.001 * Lmount_heur, 10 * Lmount_heur), (0.01 * Rlfskin_heur, 100 * Rlfskin_heur)]
    result = minimize(s_param_error_complex, rlc_test, bounds=bounds, method='Nelder-Mead', options={'maxiter': 100})

    Cnom = result.x[0]
    Rlfloss = result.x[1]
    Lmount = result.x[2]
    Rhfskin = result.x[3]
    Lhfskin = result.x[4]
    Rmfskin = result.x[5]
    Llfskin = result.x[6]
    Rlfskin = result.x[7]

    params_path = ngspice + 'paramcomplex.txt'
    with open(params_path, 'w') as paramf:
        paramf.write('.param Cnom=%.20f Rlfloss=%.20f Lmount=%.20f Rhfskin=%.20f Lhfskin=%.20f Rmfskin=%.20f Llfskin=%.20f Rlfskin=%.20f' % (Cnom, Rlfloss, Lmount, Rhfskin, Lhfskin, Rmfskin, Llfskin, Rlfskin))

    subprocess.call('ngspice.exe C_complex.cir', shell=True)

    y = circuitCalculations()

    with open(params_path, 'w') as paramf:
        paramf.write('.param Cnom=%.20f Rlfloss=%.20f Lmount=%.20f Rhfskin=%.20f Lhfskin=%.20f Rmfskin=%.20f Llfskin=%.20f Rlfskin=%.20f' % (Cnom_heur, Rlfloss_heur, Lmount_heur, Rhfskin_heur, Lmount_heur * 0.1, Rmfskin_heur, Lmount_heur * 0.1, Rlfskin_heur))

    subprocess.call('ngspice.exe C_complex.cir', shell=True)

    y_heur = circuitCalculations()


os.chdir(original_dir)

target = np.array([])

for i in range(z.size):
    target = np.append(target, z[i])


for i in range(z.size-5, z.size):
    print('Frequency: ' + str(f[i]))
    print('Generated real: ' + str(real(y[i])))
    print('Target real: ' + str(real(z[i])))

plotGraph(2, 'Target and generated absolute impedance values', [f, f], [np.abs(target), np.abs(y)], ['Target data for component', 'Fitted data(absolute value)'], 'abs_target_fitted.png', 1)
plotGraph(2, 'Target and generated equivalent series resistance', [f, f], [np.abs(np.real(target)), np.abs(np.real(y))], ['Target data for component', 'Fitted data(real value)'], 'real_target_fitted.png', 0)

# PDF report generation
now = datetime.datetime.now()

spath = str(spath)
header = spath.split('/')[-1].split('_')[0]
canvas = Canvas('report.pdf')
title = 'Component name: ' + header
canvas.drawString(100, canvas._pagesize[1] - 25, title)
canvas.drawString(100, canvas._pagesize[1] - 50, 'Input file: c.s2p')
canvas.drawString(100, canvas._pagesize[1] - 75, 'Date of generation: ' + now.strftime("%Y-%m-%d"))

if circuit == 'rlc':
    canvas.drawString(100, canvas._pagesize[1] - 100, 'Circuit type: series RLC')

    plotGraph(2, 'Values after applying heuristics', [f, f], [np.abs(target), np.abs(y_heur)], ['Target data for component', 'Data after applying heuristics'], 'heuristic.png', 1)
    canvas.drawString(100, canvas._pagesize[1] - 125, 'Heuristic values: R = %.2e Ohm, L = %.2e H, C = %.2e F' % (r_heur, l_heur, c_heur))
    canvas.drawString(100, canvas._pagesize[1] - 150, 'Values after optimization: R = %.2e Ohm, L = %.2e H, C = %.2e F' % (r_final, l_final, c_final))

if circuit == 'complex':
    canvas.drawString(100, canvas._pagesize[1] - 100, 'Circuit type: complex')

    plotGraph(2, 'Values after applying heuristics', [f, f], [np.abs(target), np.abs(y_heur)], ['Target data for component', 'Data after applying heuristics'], 'heuristic.png', 1)

    canvas.drawString(100, canvas._pagesize[1] - 125, 'Heuristic values:')
    canvas.drawString(100, canvas._pagesize[1] - 150, 'Cnom=%.2e F, Rlfloss=%.2e Ohm, Lmount=%.2e H, Rhfskin=%.2e Ohm,' % (Cnom_heur, Rlfloss_heur, Lmount_heur, Rhfskin_heur))
    canvas.drawString(100, canvas._pagesize[1] - 175, 'Lhfskin=%.2e H, Rmfskin=%.2e Ohm, Llfskin=%.2e H, Rlfskin=%.2e Ohm' % (Lmount_heur * 0.1, Rmfskin_heur, Lmount_heur * 0.1, Rlfskin_heur))
    canvas.drawString(100, canvas._pagesize[1] - 200, 'Values after optimization:')
    canvas.drawString(100, canvas._pagesize[1] - 225, 'Cnom=%.2e F, Rlfloss=%.2e Ohm, Lmount=%.2e H, Rhfskin=%.2e Ohm,' % (Cnom, Rlfloss, Lmount, Rhfskin))
    canvas.drawString(100, canvas._pagesize[1] - 250, 'Lhfskin=%.2e H, Rmfskin=%.2e Ohm, Llfskin=%.2e H, Rlfskin=%.2e Ohm' % (Lhfskin, Rmfskin, Llfskin, Rlfskin))

    canvas.showPage()
    canvas.drawString(100, canvas._pagesize[1] - 25, 'Image representation of complex circuit model:')
    canvas.drawImage('Complex_circuit.png', 50, canvas._pagesize[1] - 525, width=495, height=495)


canvas.showPage()
canvas.drawImage('heuristic.png', 50, canvas._pagesize[1] - 525, width=495, height=495)
canvas.showPage()
canvas.drawImage('abs_target_fitted.png', 50, canvas._pagesize[1] - 525, width=495, height=495)
canvas.showPage()
canvas.drawImage('real_target_fitted.png', 50, canvas._pagesize[1] - 525, width=495, height=495)
canvas.showPage()


canvas.save()