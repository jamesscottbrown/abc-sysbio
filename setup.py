from distutils.core import setup


setup(name='abc-sysbio',
      version='2.08',
      description='Approximate Bayesian Computation for systems biology',

      author='Chris Barnes',

      author_email='christopher.barnes@imperial.ac.uk',

      url='http://abc-sysbio.sourceforge.net/',

      packages=['abcsysbio'],

      scripts=['scripts/run-abc-sysbio',
               'scripts/abc-sysbio-sbml-sum'],

      requires=['libSBML',
                'matplotlib',
                'Numpy',
                'Scipy']
      )
