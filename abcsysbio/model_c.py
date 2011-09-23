import numpy
import os
import re	
from ctypes import * 

# compilation step to create share object for a correct solver and model
def compile(name, integration):
	integ = integration + "Solver"
	libname = "./lib"+name+".so.1.0"

	ABC_GSL_LIB=os.getenv("GSL_LIB")
	ABC_GSL_INC=os.getenv("GSL_INC")

	if ABC_GSL_LIB == None:
		ABC_GSL_LIB="/usr/local/lib"
	if ABC_GSL_INC == None:
		ABC_GSL_INC="/usr/local/include"

	ABC_NM_LIB=os.path.join(os.path.split(os.path.realpath(__file__))[0],'src/newmat11/')
	ABC_NM_INC=os.path.join(os.path.split(os.path.realpath(__file__))[0],'src/newmat11/')
	ABC_SRC_DIR=os.path.join(os.path.split(os.path.realpath(__file__))[0],'src/')

	command = "make -f " + ABC_SRC_DIR + "makefile --quiet "
	command = command + " MODEL=" + name + " SOLVER=" + integ + " LIBNAME=" + libname + " "
	command = command + "ABC_GSL_LIB=" + ABC_GSL_LIB + " "
	command = command + "ABC_GSL_INC=" + ABC_GSL_INC + " "
	command = command + "ABC_NM_LIB=" + ABC_NM_LIB + " "
	command = command + "ABC_NM_INC=" + ABC_NM_INC + " "
	command = command + "ABC_SRC_DIR=" + ABC_SRC_DIR + " "

	##print "COMPILE:", command
	os.system(command)
	return CDLL(libname)

	
class model:
	# instantiation
 	def __init__(self, name, nspecies, nparameters, prior, x0prior, source, integration, fit, dt, beta, initstep, relativeError, absoluteError, logp):
		gil=re.compile('Gillespie', re.I)
        	ode=re.compile('ODE', re.I)
        	euler=re.compile('Euler', re.I)
        	heun=re.compile('Heun', re.I)
        	milstein=re.compile('Milstein', re.I)		

		solverName = ""
		if gil.search(integration):
			solverName = 'Gillespie'
		elif ode.search(integration):
			solverName = 'ODE'
		elif euler.search(integration):
			solverName = 'EulerSDE'
		elif heun.search(integration):
			solverName = 'HeunSDE'
		elif milstein.search(integration):
			solverName = 'MilsteinSDE'		
		else:
			print "C model : unrecognised integrator : ", solverName

		self.name = name
		self.nspecies = nspecies

		##self.seed = numpy.random.randint(low=1, high=4294967296)
		##print "C model:", solverName, self.seed
		
		## combine the parameters with the species
		self.kparameters = nparameters
		self.nparameters = nparameters + nspecies
        
		self.prior = [x[:] for x in prior]
		for x in x0prior:
			self.prior.append( x[:] )

		self.source = source
		self.integration = solverName
		self.fit = fit 	
 		self.dt = dt	
		self.beta = beta 
		self.initstep = initstep
		self.relativeError = relativeError
		self.absoluteError = absoluteError
		self.logp = logp

		self.lib = compile(self.name, self.integration)

  	   	if self.integration=='ODE':	
  	   		self.simulate = self.simulate_ode
  	   	if self.integration=='EulerSDE':
  	   		self.simulate = self.simulate_sde
   	   	if self.integration=='MilsteinSDE':
   		  	self.simulate = self.simulate_sde
        	if self.integration=='HeunSDE':
   		  	self.simulate = self.simulate_sde
   		elif self.integration=='Gillespie':
  	 	 	self.simulate = self.simulate_gillespie

		 	
   	def simulate_ode(self, p, t, n, beta):
   	# must return a structure of the form [n][ntimepoints][nspecies]
  	# n is the number of models
 		ntimepoints = len(t)
  	 	ret = numpy.zeros( [n,beta,ntimepoints,self.nspecies]) 
		#set output from cpp file to python
		output_arr_type = beta*(self.nspecies+1)*ntimepoints*c_double
		output = output_arr_type()
  		
		#set timepoints; ctypes
		tim_arr_type = ntimepoints * c_double
		ctime = tim_arr_type() 
		for i in range (ntimepoints):
			ctime[i]=t[i]


		#set other ctypes
 		cbeta = c_int(beta)
		cntimepoints = c_int(ntimepoints)
		CNPARAMETERS = c_int(self.nparameters)
		CNSPECIES = c_int(self.nspecies)
		cinitstep = c_double(self.initstep)
		cabsoluteError = c_double(self.absoluteError)
		crelativeError = c_double(self.relativeError)
   		for ni in range(n):
			#set parameters; ctypes 
			par_arr_type = self.nparameters * c_double				
			cparam = par_arr_type() 
			for i in range (self.kparameters):
				if self.logp==False: cparam[i]=p[ni][i]
				else: cparam[i]=numpy.power(10,p[ni][i])

			# set initial values; ctypes
			init_arr_type = self.nspecies * c_double
			cinit = init_arr_type()
			j=0
			for i in range (self.kparameters,self.nparameters):      	
				cinit[j]=p[ni][i]
				j += 1

  			dat = self.lib.MainC(byref(cinit), byref(cparam), cbeta, byref(ctime), cntimepoints, CNPARAMETERS, CNSPECIES, cinitstep, cabsoluteError, crelativeError, byref(output))

			count = 0
			for j in range(beta):
				for h in range(self.nspecies):
					for k in range(ntimepoints):
						ret[ni][j][k][h] = output[count]
						count += 1
  		return ret
	

  	def simulate_sde(self, p, t, n, beta) :
  	#must return a structure of the form [n][ntimepoints][nspecies]
	# where p is an 2D array of parameters, t is a list of timepoints, n is the number of models and beta is number of model 
 		ntimepoints = len(t)
  		ret = numpy.zeros([n, beta, ntimepoints, self.nspecies])
	        #set output from cpp file to python
		output_arr_type = beta*(self.nspecies+1)*ntimepoints*c_double
		output = output_arr_type()
  		
		#set timepoints; ctypes
		tim_arr_type = ntimepoints * c_double
		ctime = tim_arr_type() 
		for i in range (ntimepoints):
			ctime[i]=t[i]
		#set other ctypes
 		cbeta = c_int(beta)
		cntimepoints = c_int(ntimepoints)
		cdt = c_double(self.dt)
		CNPARAMETERS = c_int(self.nparameters)
		CNSPECIES = c_int(self.nspecies)
		
		for ni in range(n):	
			#set parameters; ctypes 
			par_arr_type = self.nparameters * c_double				
			cparam = par_arr_type() 
			for i in range (self.kparameters):
				if self.logp==False: cparam[i]=p[ni][i]
				else: cparam[i]=numpy.power(10,p[ni][i])

			# set initial values; ctypes
			init_arr_type = self.nspecies * c_double
			cinit = init_arr_type()
			j=0
			for i in range (self.kparameters,self.nparameters):      	
				cinit[j]=p[ni][i]
				j += 1
				
			## double* initialValues, double* parameters, int beta, double* timepoints, int ntimepoints, double dt, NPARAMETERS, NSPECIES
  			dat = self.lib.MainC(byref(cinit), byref(cparam), cbeta, byref(ctime),cdt,cntimepoints, CNPARAMETERS, CNSPECIES, byref(output))

			count = 0
			for j in range(beta):
				for h in range(self.nspecies):
					for k in range(ntimepoints):
						ret[ni][j][k][h] = output[count]
						count += 1
  		return ret
       	
  	def simulate_gillespie(self, p, t, n, beta):
##  	# must return a structure of the form [n][ntimepoints][nspecies]
		ntimepoints = len(t)
  		ret = numpy.zeros([n, beta, ntimepoints, self.nspecies])
	        #set output from cpp file to python
		output_arr_type = beta*(self.nspecies+1)*ntimepoints*c_double
		output = output_arr_type()
  		
		#set timepoints; ctypes
		tim_arr_type = ntimepoints * c_double
		ctime = tim_arr_type() 
		for i in range (ntimepoints):
			ctime[i]=t[i]
		#set other ctypes
 		cbeta = c_int(beta)
		cntimepoints = c_int(ntimepoints)
		CNPARAMETERS = c_int(self.nparameters)
		CNSPECIES = c_int(self.nspecies)
		##print "info:", self.name, ntimepoints, beta, self.nparameters, self.nspecies

		for ni in range(n):	
			#set parameters; ctypes 
			par_arr_type = self.nparameters * c_double				
			cparam = par_arr_type() 
			for i in range (self.kparameters):
				if self.logp==False: cparam[i]=p[ni][i]
				else: cparam[i]=numpy.power(10,p[ni][i])

			# set initial values; ctypes
			init_arr_type = self.nspecies * c_double
			cinit = init_arr_type()
			j=0
			for i in range (self.kparameters,self.nparameters):      	
				cinit[j]=p[ni][i]
				j += 1

			self.lib.MainC(byref(cinit), byref(cparam), cbeta, byref(ctime),cntimepoints, CNPARAMETERS, CNSPECIES, byref(output))

			count = 0
			for j in range(beta):
				for h in range(self.nspecies):
					for k in range(ntimepoints):
						ret[ni][j][k][h] = output[count]
						count += 1
  		return ret



class data:
    def __init__(self, timepoints, values):
        self.timepoints = timepoints
        self.values = values
