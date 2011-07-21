import numpy

from abcsysbio import abcodesolve
from abcsysbio import sdeint
from abcsysbio import GillespieAlgorithm

class model:
    
    # instantiation
    def __init__(self, name, nspecies, nparameters, prior, x0prior, source, integration, fit, dt, atol, rtol):

        self.name = name
        self.nspecies = nspecies

        ## combine the parameters with the species
        self.kparameters = nparameters
        self.nparameters = nparameters + nspecies
        
        self.prior = [x[:] for x in prior]
        for x in x0prior:
            self.prior.append( x[:] )

        self.source = source
        self.integration = integration
        self.fit = fit

        self.module = __import__(source)
        
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
                dat = abcodesolve.abcodeint(func=self.module, InitValues=p[i][self.kparameters:self.nparameters], timepoints=t, parameters=p[i][0:self.kparameters], dt=self.dt, atol=self.atol, rtol=self.rtol )
                for k in range(len(t)):
                    for l in range(self.nspecies):
                        ret[i,j,k,l] = dat[k,l]
                        
        return ret
    
    def simulate_sde(self, p, t, n, beta):
        # must return a structure of the form [n][nbeta][ntimepoints][nspecies]
        ret = numpy.zeros( [n,beta,len(t),self.nspecies] )

        for i in range(n):
            for j in range(beta):
                dat = sdeint.sdeint(func=self.module, InitValues=p[i][self.kparameters:self.nparameters], parameter=p[i][0:self.kparameters], timepoints=t, dt=self.dt )
                for k in range(len(t)):
                    for l in range(self.nspecies):
                        ret[i,j,k,l] = dat[k,l]
                        
        return ret
    
    def simulate_mjp(self, p, t, n, beta):
        # must return a structure of the form [n][nbeta][ntimepoints][nspecies]
        ret = numpy.zeros( [n,beta,len(t),self.nspecies] )

        for i in range(n):
            for j in range(beta):
                dat = GillespieAlgorithm.GillespieInt(func=self.module, initValues=p[i][self.kparameters:self.nparameters], parameters=p[i][0:self.kparameters], outputtimes=t )
                for k in range(len(t)):
                    for l in range(self.nspecies):
                        ret[i,j,k,l] = dat[k,l]
                        
        return ret



