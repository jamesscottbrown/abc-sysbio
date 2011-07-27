
from distutils.core import setup
from distutils.command.install_data import install_data

class post_install(install_data):
    def run(self):
        # Call parent 
        install_data.run(self)
        # Execute commands
        print "Running"



setup(name='abc-sysbio',
      version='2.01',
      description='Approximate Bayesian Computation for systems biology',

      author='Chris Barnes, Juliane Liepe, Erika Cule',

      author_email='christopher.barnes@imperial.ac.uk,juliane.liepe08@imperial.ac.uk',

      url='http://abc-sysbio.sourceforge.net/',

      packages=['abcsysbio','abcsysbio_parser'],

      scripts=['scripts/run-abc-sysbio', 
               'scripts/abc-sysbio-sbml-sum'],

      package_data={ 'abcsysbio': ['MersenneTwister.dat','MersenneTwister.cu','cuLsoda_all.cu'] },

      cmdclass={"install_data": post_install},
      
      requires=['libSBML', 
                'matplotlib',
                'Numpy',
                'Scipy']
                )
