
This is the distribution of the Python package abc-sysbio.
abc-sysbio implements likelihood free parameter inference 
and model selection in dynamical systems. It is designed 
to work with both stochastic and deterministic models 
written in python or Systems Biology Markup Language (SBML). 
abc-sysbio combines three algorithms: ABC rejection sampler, 
ABC SMC for parameter inference and ABC SMC for model selection.


#################################
# CONTENTS 
#################################
1) Prerequisites
2) Linux installation
3) Mac OSX installation
4) pyCUDA on Linux
5) Package contents

#################################
# 1) Prerequisites
#################################

Before trying to install and use abc-sysbio you must
first install the following packages

	numpy
	matplotlib
	libSBML (SBML interface)
	scipy (ODE models)
	pycuda (Nvidia GPU)

The first two are essential the latter three need only
be installed if full use of abc-sysbio is required.

#################################
# 2) LINUX installation
#################################

The --prefix=<dir> option can be ignored if the installations are to the default locations.
If custom installation is required then replace <dir> with the full path to a location. This 
will be the location contaning lib and bin directories (usually /usr by default).

# python2.X http://www.python.org/download/
./configure --prefix=<dir> --enable-shared
make
make install

# numpy http://numpy.scipy.org/
python setup.py install --prefix=<dir>

# matplotlib http://matplotlib.sourceforge.net/
python setup.py build
python setup.py install --prefix=<dir>

*** For SBML interface

# swig http://www.swig.org/ (required by libsbml)
./configure --prefix=<dir>
make
make install

# libsbml http://sbml.org/Software/libSBML
./configure --with-python=<dir> --prefix=<dir> --with-swig=<dir>
make 
make install


*** For ODE modelling

# ATLAS http://math-atlas.sourceforge.net/ (required by scipy but maybe included with Linux/MAC distribution)

# scipy
export ATLAS=/full/path/to/atlas
python setup.py build
python setup.py install --prefix=<dir>



Once the required dependent packages have been installed, 
abc-sysbio is installed via the standard distutils interface. 
The command:
	
	python setup.py install 

will place the module abcsysbio into the standard folder, 
usually something /usr/lib/python2.X/site-packages/ on Linux, 
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
path. For bash type shells do:

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
# 3) Mac OSX 10.5 installation
#################################
The following assumes that you have admin rights.

1) Download and install python 2.6
http://www.python.org/ftp/python/2.6.5/python-2.6.5-macosx10.3-2010-03-24.dmg
Note that the default python that comes with OSX will not work. 

2) Download and install numpy
http://downloads.sourceforge.net/project/numpy/NumPy/1.4.1rc2/numpy-1.4.1rc2-py2.6-python.org.dmg

3) Download and install matplotlib
http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-0.99.1/matplotlib-0.99.1.1-py2.6-macosx10.5.dmg

4) Download and install scipy
http://downloads.sourceforge.net/project/scipy/scipy/0.7.2rc2/scipy-0.7.2rc2-py2.6-python.org.dmg

5) Download and install libsbml source code
http://downloads.sourceforge.net/project/sbml/libsbml/4.0.1/libsbml-4.0.1-src.zip
Do following in an xterm

   unzip libsbml-4.0.1-src.zip
   cd libsbml-4.0.1
   ./configure --with-python=/usr/local
   make
   sudo make install
   cp -r /usr/local/lib/python2.6/site-packages/libsbml* /Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/

6) Install abc-sysbio
In the unzipped abc-sysbio package directory do
python2.6 setup.py install

This places the abcsysbio package into 
/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/
and writes the scripts
/Library/Frameworks/Python.framework/Versions/2.6/bin/abc-sysbio-sbml-sum
/Library/Frameworks/Python.framework/Versions/2.6/bin/run-abc-sysbio

Add the script directory to the path (must be done in each session or added to .bashrc file)
export PATH=/Library/Frameworks/Python.framework/Versions/2.6/bin:$PATH

Then  the command

      run-abc-sysbio -h

should display a list of options and you are ready to run the examples.

#################################
# 4) pyCUDA on Linux	 
#################################

# Boost (required for pycuda) see http://wiki.tiker.net/BoostInstallationHowto
  
  ./bootstrap.sh --prefix=<dir> --libdir=<dir>/lib
  ./bjam variant=release link=shared install

# distribute-0.6.8 (required for pycuda) 
# see http://wiki.tiker.net/DistributeVsSetuptools

# Pytools (required for pycuda)
  python setup.py install

# PyCUDA
# see http://wiki.tiker.net/PyCuda/Installation/Linux

#################################
# 5) Package contents	 
#################################

abcsysbio/ is the python module containing all the source 
code for the approximate Bayesian computation algorithms

doc/ contains the package documentation Manual.pdf 

examples/ contains the SBML models and input files reproduce 
Examples 1-5 in the package documentation.

scripts/ contains two scripts run-abc-sysbio and abc-sysbio-sbml-sum, 
that enable easy use of the package for the most common applications.

MersenneTwister.dat is an initialization file for the parallel 
MersenneTwister obtained from http://www.jcornwall.me.uk

For more details, refer to the package documentation.

##############
