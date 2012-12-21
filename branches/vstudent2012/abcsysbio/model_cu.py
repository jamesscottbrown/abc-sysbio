import numpy
import copy 

class cuda_model:
    
    # instantiation
    def __init__(self, name, submodelname, nspecies, prior, pindex, source, integration, fit, dt, beta, timepoints, logp):
        self.nspecies = nspecies
        self.name = name
        self.submodelname = submodelname
        
        self.prior = [x[:] for x in prior]
        self.pindex = pindex
        
        self.nparameters = len(self.prior)
        
        self.source = source
        self.integration = integration
        self.fit = fit 
        self.cudaCode = [self.submodelname[i] +  '.cu' for i in range(len(self.submodelname))]

        self.dt = dt
        self.beta = beta
        self.timepoints = timepoints
        self.logp = logp
       
        self.modelInstance = [0 for i in range(len(self.submodelname))]
        if self.integration=='ODE':
            import cudasim.Lsoda as Lsoda
            for i in range(len(self.submodelname)):
                self.modelInstance[i] = Lsoda.Lsoda(self.timepoints[i], self.cudaCode[i], dt=self.dt) 
        elif self.integration=='SDE':
            import cudasim.EulerMaruyama as EulerMaruyama
            for i in range(len(self.submodelname)):
                self.modelInstance[i] = EulerMaruyama.EulerMaruyama(self.timepoints[i], self.cudaCode[i], beta=self.beta, dt=self.dt)
        elif self.integration=='Gillespie':
            import cudasim.Gillespie as Gillespie
            for i in range(len(self.submodelname)):
                self.modelInstance[i] = Gillespie.Gillespie(self.timepoints[i], self.cudaCode[i], beta=self.beta, dt=self.dt)
        
        #import cudasim.Gillespie as Gillespie
        #self.modelInstance = Gillespie.Gillespie(self.timepoints[0], self.cudaCode[0], beta=self.beta, dt=self.dt)


    def simulate(self, p, t, n, beta):
        # note that in this function t and beta are not used as they are specified at compile time

        result = [[] for m in range(0, len(t))]
        species = [numpy.zeros([n,self.nspecies[m]]) for m in range(0, len(t))]
        for m in range(0, len(t)):
            pp = numpy.zeros([n,len(self.pindex[m][0])])
            species[m] = numpy.zeros([n,self.nspecies[m]])      
            this_sub_params = [[] for i in range(n)]
            
            # append parameters
            for i in range(n):
                for s in self.pindex[m][0]:
                    this_sub_params[i].append(p[i][s])  
            
            # append initial conditions 
            for i in range(n):
                for s in self.pindex[m][1]:
                    this_sub_params[i].append(p[i][s])  
            
            for i in range(n):
                species[m][i,:] = this_sub_params[i][len(self.pindex[m][0]):len(this_sub_params[i])]
            
                if self.logp == False: pp[i,:] = this_sub_params[i][0:len(self.pindex[m][0])]
                else: pp[i,:] = numpy.power(10,this_sub_params[i][0:len(self.pindex[m][0])])

            result[m] = self.modelInstance[m].run(pp, species[m])
        return result

