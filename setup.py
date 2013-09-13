import os, sys
from distutils.core import setup
from distutils.command.install import install as _install

# extend the install class
class install(_install):
    def run(self):
        _install.run(self)
        v = sys.version_info
        src_dir = self.install_lib + "abcsysbio/src"
        print "src_dir", src_dir 
        comm = "cd " + src_dir +"; chmod +x install_nm.sh; ./install_nm.sh"
        os.system(comm)

setup(name='abc-sysbio',
      version='2.08',
      description='Approximate Bayesian Computation for systems biology',

      author='Chris Barnes',

      author_email='christopher.barnes@imperial.ac.uk',

      url='http://abc-sysbio.sourceforge.net/',

      packages=['abcsysbio','abcsysbio_parser'],

      scripts=['scripts/run-abc-sysbio', 
               'scripts/abc-sysbio-sbml-sum'],

      package_data={ 'abcsysbio': ['src/*'] },

      cmdclass={"install": install},
      
      requires=['libSBML', 
                'matplotlib',
                'Numpy',
                'Scipy']
                )
