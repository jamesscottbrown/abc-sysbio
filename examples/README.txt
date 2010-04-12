This folder contains the models and user input files required to reproduce the examples in Chapter 2 of the documentation.

##############

Example 1 - SIR model selection using ODEs. There should be three SBML models and a user input file.
	This example requires both scipy and libsbml installed in addition to the minmimum packages (numpy, ScientificPython, matplotlib)

Example 2 - Linear growth model using SDEs. There should be three input *.txt files and a python model file.
	No import of SBML is performed and as such this example requires only the minmimum packages (numpy, ScientificPython, matplotlib)

Example 3 - Dimerisation reaction using Markov jump process and Gillespie algorithm. There should be one SBML model and a user input file.
	This example requires libsbml installed in addition to the minmimum packages (numpy, ScientificPython, matplotlib)

Example 4 - Immigration death model using Markov jump process and Gillespie algorithm on CUDA device. There should be one cuda file and a user input file.
	This example requires pycuda installed in addition to the minmimum packages (numpy, ScientificPython, matplotlib)

