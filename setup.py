
from distutils.core import setup

setup(name='abc-sysbio',
      version='0.2',
      description='Approximate Bayesian Computation for systems biology',
      author='Juliane Liepe, Erika Cule, Chris Barnes',
      author_email='juliane.liepe08@imperial.ac.uk, christopher.barnes@imperial.ac.uk',
      url='http://abc-sysbio.sourceforge.net/',
      packages=['abcsysbio'],
      scripts=['scripts/run-abc-sysbio', 
               'scripts/abc-sysbio-sbml-sum'],
      requires=['libSBML', 
                'math',
                'matplotlib',
                'Numpy',
                'os',
                'pickle',
                're',
                'sys',
                'Scientific',
                'Scipy']
                )
