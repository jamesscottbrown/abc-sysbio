
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
3) Mac OS X installation

#################################
# 1) Prerequisites
#################################

Before trying to install and use abc-sysbio you must
first install the following packages

	numpy
	matplotlib
	libSBML (SBML interface)
	scipy (ODE models)
	cuda-sim (Nvidia GPU)

While the first two are essential, the latter three need only
be installed if full use of abc-sysbio is required.

#################################
# 2) LINUX installation
#################################
If custom installation is required then replace <dir> 
with the full path to a location. This will be the 
location containing lib and bin directories (usually 
/usr/local by default).

The --prefix=<dir> option is recommended since it will 
guarantee that each package picks up the correct dependency.

These instructions are for old versions of packages and the 
newest versions should be used in their place.

1) Download and install python
http://www.python.org/ftp/python/2.6.5/Python-2.6.5.tgz
	
	tar xzf Python-2.6.5.tgz
	cd Python-2.6.5
	./configure --prefix=<dir> --enable-shared
	make
	make install

Note that the user may have to do 

	sudo make install

if the target <dir> is write protected

Make sure this new version of python is picked up by default
	
	export PATH=<dir>/bin:$PATH  (for BASH)
	setenv PATH <dir>/bin:$PATH

2) Download and install numpy
http://downloads.sourceforge.net/project/numpy/NumPy/1.4.1rc2/numpy-1.4.1rc2.tar.gz
	
	tar xzf numpy-1.4.1rc2.tar.gz	
	cd numpy-1.4.1rc2
	python setup.py install --prefix=<dir>

3) Download and install matplotlib
http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-0.99.1/matplotlib-0.99.1.2.tar.gz

	tar xzf matplotlib-0.99.1.2.tar.gz
	cd matplotlib-0.99.1.1/	
	python setup.py build
	python setup.py install --prefix=<dir>

4) Download and install swig 
Note that this is required by libsbml and it must be at least version 1.3.39 
http://downloads.sourceforge.net/project/swig/swig/swig-1.3.40/swig-1.3.40.tar.gz

	tar -xzf swig-1.3.40.tar.gz
	cd swig-1.3.40
	./configure --prefix=<dir>
	make
	make install

Note that the user may have to do

	sudo make install

if the target <dir> is write protected

5) Download and install libSBML
http://downloads.sourceforge.net/project/sbml/libsbml/4.0.1/libsbml-4.0.1-src.zip

	unzip libsbml-4.0.1-src.zip
	cd libsbml-4.0.1	
	./configure --with-python=<dir> --prefix=<dir> --with-swig=<dir>
	make 
	make install

Note that the user may have to do 

	sudo make install

if the target <dir> is write protected

6) Download and install scipy
http://downloads.sourceforge.net/project/scipy/scipy/0.7.2rc2/scipy-0.7.2rc2.tar.gz

Note that scipy requires ATLAS (http://math-atlas.sourceforge.net/)
but should be included in most linux distributions. You need to locate
the ATLAS libraries

	tar xzf scipy-0.7.2rc2.tar.gz
	cd scipy-0.7.2rc2
	export ATLAS=/full/path/to/atlas/libraries
	python setup.py build
	python setup.py install --prefix=<dir>

7) Download and install cuda-sim (optional)

8) Install abc-sysbio
In the unzipped abc-sysbio package directory do
   
	python setup.py install --prefix=<dir>

This places the abcsysbio package into 
     
	<dir>/lib/python2.6/site-packages/

and writes the scripts
    
	<dir>/bin/abc-sysbio-sbml-sum
	<dir>/bin/run-abc-sysbio

Add the script directory to the path (must be done in each session or added to .bashrc or .cshrc file)

	export PATH=<dir>/bin:$PATH  (bash shells)
	setenv PATH <dir>/bin:$PATH  (c shells)

Then  the command

	run-abc-sysbio -h

should display a list of options and you are ready to run the examples.


#################################
# 3) Mac OSX 10.8 installation
#################################

What follows is one way to install abc-sysbio on Mac OSX 10.8 (though it should work on 10.7 as well). It assumes that you have admin rights. 

Before installing ABC-SysBio we need to install Python the following packages:
Numpy, Scipy, Matplotlib, libSBML 

Luckily the first three can be obtained easily by installing the Scipy Superpack 
http://fonnesbeck.github.io/ScipySuperpack/

To obtain a compatible version of libSBML it is best to install from source. Download libSBML from http://sourceforge.net/projects/sbml/files/libsbml/ and unzip it. In the terminal go to the libSBML directory and type:

   	  ./configure --with-python=/usr/ --prefix=/usr/ --enable-m64
   	  make
   	  sudo make install
   	  sudo cp -r /usr/lib/python2.7/site-packages/* /Library/Python/2.7/site-packages/

Assuming that the previous steps were completed successfully, we can now install ABC-SysBio. Download the ABC-SysBio package from http://sourceforge.net/projects/abc-sysbio/files/ and unzip it. Open a terminal and type:

	 cd abc-sysbio-2.07
	 sudo python setup.py install 

This places the ABC-SysBio package into 

     	 /Library/Python/2.7/site-packages/

and writes the scripts

    	 /usr/bin/abc-sysbio-sbml-sum
	 /usr/bin/run-abc-sysbio

Then the command

     	 run-abc-sysbio -h

will display a list of options and you are ready to run the examples.
