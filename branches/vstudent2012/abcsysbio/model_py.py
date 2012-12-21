import numpy
import copy 

from abcsysbio import abcodesolve
from abcsysbio import sdeint
from abcsysbio import GillespieAlgorithm

class model:
    
    # instantiation
    def __init__(self, name, submodelname, nspecies, prior, pindex, source, integration, fit, dt, atol, rtol, logp):

        self.name = name
        self.submodelname = submodelname
        self.nspecies = nspecies

        self.prior = [x[:] for x in prior]
        self.pindex = pindex
        
        self.nparameters = len(self.prior) 

        self.source = source
        self.integration = integration
        self.fit = fit
        
        self.module = [__import__(source[i]) for i in range(0, len(self.submodelname))] 
                
        self.dt = dt
        self.atol = atol
        self.rtol = rtol
        self.logp = logp

        if self.integration=='ODE':
            self.simulate = self.simulate_ode
        elif self.integration=='SDE':
            self.simulate = self.simulate_sde
        elif self.integration=='Gillespie':
            self.simulate = self.simulate_mjp

        
    def simulate_ode(self, p, t, n, beta):
        # must return a structure of the form [nsubmodel][n][nbeta][ntimepoints][nspecies]

        nsubmodels = len(t) 
        ret = [numpy.zeros( [n,beta,len(t[m]),self.nspecies[m]] ) for m in range(nsubmodels)]
        dat = [[] for m in range(nsubmodels)]
        
        for m in range(nsubmodels): 
            
            this_sub_params = [[] for i in range(n)]
            
            # append parameters
            for i in range(n):
                for s in self.pindex[m][0]:
                    this_sub_params[i].append(p[i][s]) 
                    
            # append initial conditions
            for i in range(n):
                for s in self.pindex[m][1]:
                    this_sub_params[i].append(p[i][s])  

            nparameters = len(this_sub_params[0]) 
            
            for i in range(n):
                for j in range(beta):
                    if self.logp == False: par = this_sub_params[i][0:len(self.pindex[m][0])]
                    else: par = numpy.power( 10, this_sub_params[i][0:len(self.pindex[m][0])] )
                    dat[m] = abcodesolve.abcodeint(func=self.module[m], InitValues=this_sub_params[i][len(self.pindex[m][0]):nparameters], timepoints=t[m], parameters=par, dt=self.dt, atol=self.atol, rtol=self.rtol )  
                    for k in range(len(t[m])):
                        for l in range(self.nspecies[m]):
                            ret[m][i][j][k][l] = dat[m][k,l]


                        
        return ret  
 
    def simulate_sde(self, p, t, n, beta):  
        # must return a structure of the form [nsubmodels][n][nbeta][ntimepoints][nspecies]

        nsubmodels = len(t) 
        ret = [numpy.zeros( [n,beta,len(t[m]),self.nspecies[m]] ) for m in range(nsubmodels)]

        dat = [[] for m in range(nsubmodels)]

        for m in range(nsubmodels):  
            this_sub_params = [[] for i in range(n)]
            
            for i in range(n):
                for s in self.pindex[m][0]:
                    this_sub_params[i].append(p[i][s])
            
            for i in range(n):
                for s in self.pindex[m][1]:
                    this_sub_params[i].append(p[i][s])

            nparameters = len(this_sub_params[0])
            
            for i in range(n):
                for j in range(beta):
                    if self.logp == False: par = this_sub_params[i][0:len(self.pindex[m][0])]
                    else: par = numpy.power( 10, this_sub_params[i][0:len(self.pindex[m][0])] )

                    dat[m] = sdeint.sdeint(func=self.module[m], InitValues=this_sub_params[i][len(self.pindex[m][0]):nparameters], parameter=par, timepoints=t[m], dt=self.dt )

                    for k in range(len(t[m])):
                        for l in range(self.nspecies[m]):
                            ret[m][i][j][k][l] = dat[m][k,l]

        return ret 
    
   
   
   
    def simulate_mjp(self, p, t, n, beta):
        # must return a structure of the form [nsubmodels][n][nbeta][ntimepoints][nspecies]

        nsubmodels = len(t)  
        ret = [numpy.zeros( [n,beta,len(t[m]),self.nspecies[m]] ) for m in range(nsubmodels)]
        dat = [[] for m in range(nsubmodels)]
        
        for m in range(0, len(self.submodelname)): 
            this_sub_params = [[] for i in range(n)]
            
            # append parameters
            for i in range(n):
                for s in self.pindex[m][0]:
                    this_sub_params[i].append(p[i][s])
            
            # append initial conditions
            for i in range(n):
                for s in self.pindex[m][1]:
                    this_sub_params[i].append(p[i][s])  
                        
            nparameters = len(this_sub_params[i])
            
            for i in range(n):
                for j in range(beta):
                    if self.logp == False: par = this_sub_params[i][0:len(self.pindex[m][0])]
                    else: par = numpy.power( 10, this_sub_params[i][0:len(self.pindex[m][0])] )
                    dat[m] = GillespieAlgorithm.GillespieInt(func=self.module[m], initValues=this_sub_params[i][len(self.pindex[m][0]):nparameters], parameters=par, outputtimes=t[m] )

                    for k in range(len(t[m])):
                            for l in range(self.nspecies[m]):
                                ret[m][i,j,k,l] = dat[m][k,l]
                        
        return ret



