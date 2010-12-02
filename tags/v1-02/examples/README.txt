This folder contains the models and user input files required to reproduce the examples in Chapter 2 of the documentation.

##############

Example 1 - SIR model selection using ODEs. There should be three SBML models and a user input file.
	Requires scipy, libsbml, numpy, matplotlib
	To run:
	   cd Example1
	   run-abc-sysbio --infile input_file_SIR.txt -f

Example 2 - Linear growth model using SDEs. There should be three input *.txt files and a python model file.
	Requires numpy, matplotlib
	To run:
	   cd Example2
	   run-abc-sysbio --infile input_file_lingrow.txt -pm -f

Example 3 - Dimerisation reaction using Markov jump process and Gillespie algorithm. There should be one SBML model and a user input file.
	Requires libsbml, numpy, matplotlib
	To run:
	   cd Example3
	   run-abc-sysbio -f -i input_file_abcRejection.txt 

Example 4 - Immigration death model using Markov jump process and Gillespie algorithm on CUDA device. There should be one cuda file and a user input file.
	Requires pycuda, numpy, matplotlib
	To run this on Tesla card
	   cd Example4
	   run-abc-sysbio -f -i input_file_pop.txt -sd=10 --runmode=2 -of=results -ct=512 -cb=16 --cudacode=immdeath_logistic.cu -pm -t

	NOTE: depending on your system, the number of threadblocks (-cb) and the number of threads (-ct) may have to be adjusted

Example 5 - Repressilator model using ODE solver. There should be one user input file and the python module repressilator.py. 
	Requires scipy, numpy, matplotlib
	To run: 
	   cd Example5
	   run-abc-sysbio -i input_file_repressilator.txt -f -pm


