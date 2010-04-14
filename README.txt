This is the distribution of the Python package abc-sysbio.
abc-sysbio implements likelihood free parameter inference 
and model selection in dynamical systems. It is designed 
to work with both stochastic and deterministic models 
written in Systems Biology Markup Language (SBML). 
abc-sysbio is a Python package that combines three algorithms: 
ABC rejection sampler, ABC SMC for parameter inference and 
ABC SMC for model selection.

#################################
# PREREQUISITES
#################################

Before trying to install and use abc-sysbio you MUST 
first install the following packages

	numpy
	matplotlib

If you wish to use the SBML interface in abc-sysbio 
you must also install

	libSBML

If you wish to use abc-sysbio for inference in ODE 
models you must also install
	
	scipy

If you wish to use abc-sysbio on Nvidia GPUs you 
must also install
	
	pycuda

For instructions and help on how to install these 
packages see DEPENDENCIES in the doc/ folder.


#################################
# LINUX/MAC INSTALLATION 
#################################

Once the required dependent packages have been installed, 
abc-sysbio is installed via the standard distutils interface. 
The command:
	
	python setup.py install 

will place the module abcsysbio into the standard folder, 
usually something /usr/lib/python2.6/site-packages/ on Linux, 
depending on the version of python you are using. It also 
places two executable Python scripts called run-abc-sysbio 
and abc-sysbio-sbml-sum into the standard location, 
usually /usr/bin on Linux.

If you wish to install to a custom location you can use 
the command:

	python setup.py install \
		--prefix=/path/to/custom/location/

this places the scripts run-abc-sysbio and abc-sysbio-sbml-sum 
into /path/to/custom/location/bin/ and the module abcsysbio 
in /path/to/custom/location/lib/python2.6/site-packages/
depending on which python version you are using.


The current directory should also be placed into the python 
path. For bash type shells (default on MAC) do:

	export PYTHONPATH=$PYTHONPATH:.

For C type shells do:

	setenv PYTHONPATH $PYTHONPATH:.

Note that this can be added to your .bashrc or .cshrc.

If everything is installed correctly it should then suffice 
to do:

	run-abc-sysbio -h

to get a list of options printed. 
In addition the abcsysbio module can be used like any other 
Python library.



#################################
# CONTENTS OF PACKAGE		 
#################################

abcsysbio/ is the python module containing all the source 
code for the approximate Bayesian computation algorithms

doc/ contains the package documentation Manual.pdf and 
DEPENDENCIES which contains instructions on how to install 
dependent packages.

examples/ contains the SBML models and input files reproduce 
Examples 1, 2, 3 and 4 in the package documentation.

scripts/ contains two scripts run-abc-sysbio and abc-sysbio-sbml-sum, 
that enable easy use of the package for the most common applications.

MersenneTwister.dat is an initialization file for the parallel 
MersenneTwister obtained from http://www.jcornwall.me.uk

For more details, refer to the package documentation.

##############
