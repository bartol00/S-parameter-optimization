# Generating equivalent impedance parameters of a circuit based on its scatter parameters 

The program is designed to get the impedance data of a capacitor based on its scatter parameters. It works by converting the
supplied scatter parameters into equivalent impedance (Z) parameters and using the Nelder Mead optimization algorithm to get a series 
RLC or complex* circuit with an equivalent absolute impedance for all frequency points. The number of frequency points for the target
and generated parameters is adjusted automatically via interpolation.

The optimization works using heuristics which are calculated automatically and are based on the target data to which we want to match the
generated data. 

Depending on the chosen circuit model (either rlc or complex), the generated data may vary. The absolute impedance of the generated data 
will be matched as closely as possible to the target in both circuit models. The equivalent real/resistive data will be generated much more
accurately using the complex circuit model, with the drawback of a higher time complexity. The decision of which circuit model to use 
should be made in regards to this information.

The graphs generated for an example circuit are stored in the 'images/' folder, while the example S-parameters for circuits are stored in the 
's_param_samples/' folder.



*The complex circuit model refers to the model pictured in the generated report if the chosen circuit model was complex.



Information about working with config.txt:

The following variables must be specified in order for the program to working

Model: the kind of model the S param data references, it can be series or shunt (default is shunt)
Ngspice: the ngspice bin directory, where ngspice.exe is located
Spath: the path to the S params to be analyzed
Circuit: the type of circuit to be analyzed, can be rlc or complex

Additional line comments can be added by putting '#' before a line
