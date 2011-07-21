import numpy
import os
import re	
from ctypes import * 

# compilation step to create share object for a correct solver and model
def compile(name, integration):
	integ = integration + "Solver"
	libname = "lib"+name+"Model.so.1.0"
	command = "cd /cluster/home/cbarnes/dev/abc-sysbio-area/abc-sysbio/trunk/src; make MODEL=" + name+"Model" + " SOLVER=" + integ + " LIBNAME=" + libname
	
	os.system(command)
	return CDLL("/cluster/home/cbarnes/dev/abc-sysbio-area/abc-sysbio/trunk/src"+libname)

	
class model:
	# instantiation
 	def __init__(self, name, nspecies, nparameters, prior, x0prior, source, integration, fit, dt, beta, initstep, relativeError, absoluteError):
		gil=re.compile('Gillespie', re.I)
        	ode=re.compile('ODE', re.I)
        	euler=re.compile('Euler', re.I)
        	heun=re.compile('Heun', re.I)
        	milstein=re.compile('Milstein', re.I)		

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

		self.name = name
		self.nspecies = nspecies

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
  		
		# set initial values; ctypes
		init_arr_type = self.nspecies * c_double
		cinit = init_arr_type()
		for i in range (self.nspecies):      	
			cinit[i]=self.init[i]
	
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
			for i in range (self.nparameters):
				cparam[i]=p[ni][i]

  			dat = self.lib.MainC(byref(cinit), byref(cparam), cbeta, byref(ctime), cntimepoints, CNPARAMETERS, CNSPECIES, cinitstep, cabsoluteError, crelativeError, byref(output))
				
			iterationNumber = 0
			specieNb = 0
			timepoint = 0
			for i in range(beta*(self.nspecies+1)*ntimepoints):
					if(i<ntimepoints):
						pass				
					else:
						if(i>ntimepoints and (i-ntimepoints)%(ntimepoints*self.nspecies)==0):
							iterationNumber = iterationNumber + 1
							specieNb = 0
						if(i>ntimepoints and i%ntimepoints==0):
							specieNb = specieNb + 1
							timepoint = 0
	
						ret[ni][iterationNumber][timepoint][specieNb] = output[i]
						timepoint = timepoint + 1
  		return ret
	

  	def simulate_sde(self, p, t, n, beta) :
  	#must return a structure of the form [n][ntimepoints][nspecies]
	# where p is an 2D array of parameters, t is a list of timepoints, n is the number of models and beta is number of model 
 		ntimepoints = len(t)
  		ret = numpy.zeros([n, beta, ntimepoints, self.nspecies])
	        #set output from cpp file to python
		output_arr_type = beta*(self.nspecies+1)*ntimepoints*c_double
		output = output_arr_type()
  		
		# set initial values; ctypes
		init_arr_type = self.nspecies * c_double
		cinit = init_arr_type()
		for i in range (self.nspecies):      	
			cinit[i]=self.init[i]
		
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
			for i in range (self.nparameters):
				cparam[i]=p[ni][i]
			## double* initialValues, double* parameters, int beta, double* timepoints, int ntimepoints, double dt, NPARAMETERS, NSPECIES
  			dat = self.lib.MainC(byref(cinit), byref(cparam), cbeta, byref(ctime),cdt,cntimepoints, CNPARAMETERS, CNSPECIES, byref(output))
			iterationNumber = 0
			specieNb = 0
			timepoint = 0
			for i in range(beta*(self.nspecies+1)*ntimepoints):
				if(i<ntimepoints):
					pass				
				else:
					if(i>ntimepoints and (i-ntimepoints)%(ntimepoints*self.nspecies)==0):
						iterationNumber = iterationNumber + 1
						specieNb = 0
					if(i>ntimepoints and i%ntimepoints==0):
						specieNb = specieNb + 1
						timepoint = 0

					ret[ni][iterationNumber][timepoint][specieNb] = output[i]
					timepoint = timepoint + 1
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
		for ni in range(n):	
			#set parameters; ctypes 
			par_arr_type = self.nparameters * c_double				
			cparam = par_arr_type() 
			for i in range (self.kparameters):
				cparam[i]=p[ni][i]

			# set initial values; ctypes
			init_arr_type = self.nspecies * c_double
			cinit = init_arr_type()
			for i in range (self.kparameters,self.nparameters):      	
				cinit[i]=p[ni][i]

			self.lib.MainC(byref(cinit), byref(cparam), cbeta, byref(ctime),cntimepoints, CNPARAMETERS, CNSPECIES, byref(output))
			iterationNumber = 0
			specieNb = 0
			timepoint = 0
			for i in range(beta*(self.nspecies+1)*ntimepoints):
				if(i<ntimepoints):
					pass				
				else:
					if(i>ntimepoints and (i-ntimepoints)%(ntimepoints*self.nspecies)==0):
						iterationNumber = iterationNumber + 1
						specieNb = 0
					if(i>ntimepoints and i%ntimepoints==0):
						specieNb = specieNb + 1
						timepoint = 0

					ret[ni][iterationNumber][timepoint][specieNb] = output[i]
					timepoint = timepoint + 1

  		return ret



class data:
    def __init__(self, timepoints, values):
        self.timepoints = timepoints
        self.values = values
