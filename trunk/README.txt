
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
4) Windows Vista installation
5) pycuda on Linuxa
6) Package contents

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

7) Install abc-sysbio
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
# 3) Mac OSX 10.5 installation
#################################

What follows is one way to install abc-sysbio on Mac OSX 10.5. It assumes that you have admin rights. 

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
# 4) Windows Vista installation	 
#################################

Windows is not currently supported but we have succesfully installed
matplotlib, numpy, scipy, libsbml and abc-sysbio on Windows Vista
using WinPython.


