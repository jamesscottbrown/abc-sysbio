import numpy

from abcsysbio import abcodesolve
from abcsysbio import sdeint
from abcsysbio import GillespieAlgorithm

class model:
    
    # instantiation
    def __init__(self, name, nspecies, nparameters, prior, source, integration, fit, init, dt, atol, rtol):
        self.nspecies = nspecies
        self.nparameters = nparameters
        self.name = name
        self.prior = [x[:] for x in prior]
        self.source = source
        self.integration = integration
        self.fit = fit

        self.module = __import__(source)
        self.init = init
        self.dt = dt
        self.atol = atol
        self.rtol = rtol

        if self.integration=='ODE':
            self.simulate = self.simulate_ode
        elif self.integration=='SDE':
            self.simulate = self.simulate_sde
        elif self.integration=='Gillespie':
            self.simulate = self.simulate_mjp

    def simulate_ode(self, p, t, n, beta):
        # must return a structure of the form [n][nbeta][ntimepoints][nspecies]
        ret = numpy.zeros( [n,beta,len(t),self.nspecies] )

        for i in range(n):
            for j in range(beta):
                dat = abcodesolve.abcodeint(func=self.module, InitValues=self.init, timepoints=t, parameters=p[i], dt=self.dt, atol=self.atol, rtol=self.rtol )
                for k in range(len(t)):
                    for l in range(self.nspecies):
                        ret[i,j,k,l] = dat[k,l]
                        
        return ret
    
    def simulate_sde(self, p, t, n, beta):
        # must return a structure of the form [n][nbeta][ntimepoints][nspecies]
        ret = numpy.zeros( [n,beta,len(t),self.nspecies] )

        for i in range(n):
            for j in range(beta):
                dat = sdeint.sdeint(func=self.module, InitValues=self.init, parameter=p[i], timepoints=t, dt=self.dt )
                for k in range(len(t)):
                    for l in range(self.nspecies):
                        ret[i,j,k,l] = dat[k,l]
                        
        return ret
    
    def simulate_mjp(self, p, t, n, beta):
        # must return a structure of the form [n][nbeta][ntimepoints][nspecies]
        ret = numpy.zeros( [n,beta,len(t),self.nspecies] )

        for i in range(n):
            for j in range(beta):
                dat = GillespieAlgorithm.GillespieInt(func=self.module, initValues=self.init, parameters=p[i], outputtimes=t )
                for k in range(len(t)):
                    for l in range(self.nspecies):
                        ret[i,j,k,l] = dat[k,l]
                        
        return ret

class cuda_model:
    
    # instantiation
    def __init__(self, name, nspecies, nparameters, prior, source, integration, fit, init, dt, beta, timepoints):
        self.nspecies = nspecies
        self.nparameters = nparameters
        self.name = name
        self.prior = [x[:] for x in prior]
        self.source = source
        self.integration = integration
        self.fit = fit
        self.cudaCode = self.name +  '.cu' 
        self.init = init
        self.dt = dt
        self.beta = beta
        self.timepoints = timepoints
       
        if self.integration=='ODE':
            import cudasim.Lsoda as Lsoda
            self.modelInstance = Lsoda.Lsoda(self.timepoints, self.cudaCode, dt=self.dt)
        elif self.integration=='SDE':
            import cudasim.EulerMaruyama as EulerMaruyama
            self.modelInstance = EulerMaruyama.EulerMaruyama(self.timepoints, self.cudaCode, beta=self.beta, dt=self.dt)
        elif self.integration=='Gillespie':
            import cudasim.Gillespie as Gillespie
            self.modelInstance = Gillespie.Gillespie(self.timepoints, self.cudaCode, beta=self.beta, dt=self.dt)

    def simulate(self, p, t, n, beta):
        # note that in this function t and beta are not used as they are specified at compile time
        
        species = []
        for i in range(n):
            species.append( self.init )
        
        result = self.modelInstance.run(p, species)
        return result
    

class data:
    def __init__(self, timepoints, values):
        self.timepoints = timepoints
        self.values = values

        #print self.timepoints
        #print self.values
