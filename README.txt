This is the distribution of the Python package abc-sysbio

#################################
# INSTALLATION
python setup.py install ( --prefix=/path/to/custom/location/ )

This command will place the module abc into the default	(custom) installation location.
It also places two executable Python scripts (abc-sysbio and abc-sysbio-sbml-sum) into the
default (custom) location.		   

#################################
# CONTENTS			 

##############

abc is a module containing all the Approximate Bayesian Computation code
The abc module requires the following Python libraries. With the exception of libSBML, these are included in the Enthought Python Distribution
libSBML
math
matplotlib
Numpy
os
pickle
re
sys
Scientific
Scipy
For information on how to use the package, refer to the package documentation.

##############

doc contains the package documentation

##############

Examples contains the SBML models, input files and (for Example 2) modified scripts required to reproduce Examples 1, 2, 3 and 4 in the package documentation.

##############

scripts contains two scripts, abc-sysbio and abc-sysbio-sbml-sum, that enable easy use of the package for the most common applications.

##############

data contains MersenneTwister.dat which is an initialization file for the parallel Mersenne Twister obtained from (http://www.jcornwall.me.uk)

##############

For more details, refer to the package documentation.

##############
