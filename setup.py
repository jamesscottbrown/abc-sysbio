
from distutils.core import setup

setup(name='abc-sysbio',
      version='0.1',
      description='Approximate Bayesian Computation for systems biology',
      author='Juliane Liepe and Erika Cule',
      author_email='juliane.liepe08@imperial.ac.uk',
      url='http://abc-sysbio.sourceforge.net/',
      packages=['abc'],
      scripts=['scripts/abc-sysbio', 
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
