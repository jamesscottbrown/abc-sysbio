import numpy
from numpy import random as rnd

import copy, time

from abcsysbio import euclidian
from abcsysbio import kernels
from abcsysbio import statistics

class abcsmc_results:
    def __init__(self, 
                 naccepted, 
                 sampled,
                 rate, 
                 trajectories, 
                 distances, 
                 margins, 
                 models, 
                 weights, 
                 parameters,
                 epsilon):
        self.naccepted = naccepted
        self.sampled = sampled
        self.rate = rate
        self.trajectories = trajectories
        self.distances = distances
        self.margins = margins
        self.models = models
        self.weights = weights
        self.parameters = parameters
        self.epsilon = epsilon

class abcsmc:

    # instantiation
    def __init__(self, 
                 models, 
                 nparticles, 
                 modelprior, 
                 data,
                 beta,
                 nbatch,
                 modelKernel,
                 debug,
                 distancefn = euclidian.euclidianDistance,
                 kernelfn = kernels.getKernel,
                 kernelpdffn = kernels.getPdfParameterKernel,
                 perturbfn = kernels.perturbParticle):
        
        self.nmodel = len(models)
        self.models = copy.copy( models )
        
        self.nparticles = nparticles

        self.model_prev      = [0  for i in range(nparticles)]
        self.weights_prev    = [0  for i in range(nparticles)]
        self.parameters_prev = [[] for i in range(nparticles)]
        self.margins_prev    = [0  for i in range(self.nmodel)] 
        
        self.model_curr      = [0  for i in range(nparticles)]
        self.weights_curr    = [0  for i in range(nparticles)]
        self.parameters_curr = [[] for i in range(nparticles)]
        self.margins_curr    = [0  for i in range(self.nmodel)] 
        
        self.b = [0  for i in range(0,nparticles)]

        self.distancefn = distancefn
        self.kernelfn = kernels.getKernel
        self.kernelpdffn = kernels.getPdfParameterKernel
        self.perturbfn = kernels.perturbParticle

        self.beta = beta
        self.dead_models = []
        self.nbatch = nbatch
        self.debug = debug
        self.data = copy.deepcopy(data)
    
        self.modelprior = modelprior[:]
        self.modelKernel = modelKernel
        self.kernels = []
        for i in range(self.nmodel):
            self.kernels.append([])
            for j in range(self.models[i].nparameters):
                self.kernels[i].append( [1, -1, 1 ] ) # uniform kernels
            
        self.hits = []
        self.sampled = []
        self.rate = []
        self.dead_models = []
        self.sample_from_prior = True

    def run_fixed_schedule(self, epsilon, io):
        for pop in range(len(epsilon)):
            start_time = time.time()
            if(pop==0 and self.sample_from_prior==True): 
                results = self.iterate_one_population(epsilon[pop], prior=True)
            else:
                results = self.iterate_one_population(epsilon[pop], prior=False)
            end_time = time.time()
             
            io.write_pickled(self.nmodel, self.model_prev, self.weights_prev, self.parameters_prev, self.margins_prev, self.kernels)
            io.write_data(pop, results, end_time-start_time, self.models, self.data)

            if self.debug == 1:
                print "### population ", pop+1
                print "\t sampling steps / acceptance rate :", self.sampled[pop], "/", self.rate[pop]
                print "\t model marginals:", self.margins_prev
                print "\t dead models    :", self.dead_models
        return
            

    def iterate_one_population(self, next_epsilon, prior):
        if self.debug == 2:print "\n\n****iterate_one_population: next_epsilon, prior", next_epsilon, prior

        naccepted = 0
        sampled = 0

        while(naccepted < self.nparticles):
            if self.debug == 2:print "\t****batch"
            if( prior == False):
                sampled_models = self.sampleTheModel()
                sampled_params = self.sampleTheParameter(sampled_models)
            else:
                sampled_models = self.sampleTheModelFromPrior()
                sampled_params = self.sampleTheParameterFromPrior(sampled_models)
                
            accepted_index, trajectories, distances = self.simulate_and_compare_to_data(sampled_models,sampled_params,next_epsilon)

            for i in range(self.nbatch):
                if naccepted < self.nparticles:
                    sampled = sampled + 1

                if naccepted < self.nparticles and accepted_index[i] > 0 : 
                    
                    self.model_curr[naccepted] = sampled_models[i]
                    if self.debug == 2:print "\t****accepted", i, accepted_index[i], sampled_models[i]
                    
                    for p in range( self.models[ sampled_models[i] ].nparameters ):
                        self.parameters_curr[naccepted].append(sampled_params[i][p])

                    self.b[naccepted] = accepted_index[i]
                    naccepted = naccepted + 1

            
            if self.debug == 2:print "\t****end  batch naccepted/sampled:", naccepted,  sampled

        # Finished loop over particles
        if self.debug == 2:print "**** end of population naccepted/sampled:", naccepted,  sampled

        if( prior == False):
            self.computeParticleWeights()
        else:
            for i in range(self.nparticles):
                self.weights_curr[i] = self.b[i]
        
        self.normalizeWeights()
        self.modelMarginals()
  
        if self.debug == 2:
            print "**** end of population: particles"
            for i in range(self.nparticles):
                print i, self.weights_curr[i], self.model_curr[i], self.parameters_curr[i]
            print self.margins_curr

        #
        # Prepare for next population
        # 
        self.margins_prev = self.margins_curr[ : ]
        self.weights_prev = self.weights_curr[ : ]
        self.model_prev = self.model_curr[ : ]
        self.parameters_prev = []
        for i in range(self.nparticles):
            self.parameters_prev.append(self.parameters_curr[i][ : ])

        self.model_curr      = [0  for j in range(0,self.nparticles)]
        self.weights_curr    = [0  for j in range(0,self.nparticles)]
        self.parameters_curr = [[] for j in range(0,self.nparticles)]
        self.margins_curr    = [0  for j in range(0,self.nmodel)] 

        self.b = [0  for i in range(0,self.nparticles)]

        #
        # Check for dead models
        #
        self.dead_models = []
        for j in range(self.nmodel):
            if self.margins_prev[j] < 1e-6:
                self.dead_models.append(j)
        
        #
        # Compute kernels
        #
        for mod in range(self.nmodel):
            this_model_index = numpy.arange(self.nparticles)[ numpy.array(self.model_prev) == mod ]
            this_population = numpy.zeros([ len(this_model_index), self.models[mod].nparameters ])
           
            if len(this_model_index) > 0:
                for it in range(len(this_model_index)):
                    this_population[it,:] = self.parameters_prev[ this_model_index[it] ][:]
            
                self.kernels[mod] = self.kernelfn( self.kernels[mod], this_population )[:]

        self.hits.append( naccepted )
        self.sampled.append( sampled )
        self.rate.append( naccepted/float(sampled) )

        results = abcsmc_results(naccepted, 
                                 sampled,
                                 naccepted/float(sampled), 
                                 trajectories, 
                                 distances, 
                                 self.margins_prev, 
                                 self.model_prev, 
                                 self.weights_prev, 
                                 self.parameters_prev,
                                 next_epsilon)

        return results

    # end of iterate_one_population
       
    def fill_values(self, particle_data ):
        # particle data comprises of [model_pickled, weights_pickled, parameters_pickled, margins_pickled, kernel]
        self.model_prev = particle_data[0][:]
        self.weights_prev = particle_data[1][:]
        self.margins_prev = particle_data[3][:]

        self.parameters_prev = []
        for it in range(len(particle_data[2])):
            self.parameters_prev.append( particle_data[2][it][:] )
        
        self.kernels = []
        for i in range(self.nmodel):
            self.kernels.append([])
            for j in range(self.models[i].nparameters):
                self.kernels[i].append( particle_data[4][i][j][:] )
        
        self.sample_from_prior = False

    def simulate_and_compare_to_data(self, sampled_models, sampled_params, this_epsilon):
        # Here do the simulations for each model together
        if self.debug == 2:print '\t\t\t***simulate_and_compare_to_data'

        ret = [0 for it in range(self.nbatch)]
        traj = []
        distances = []

        n_to_simulate = [0 for it in range(self.nmodel)]
        mods = numpy.array( sampled_models )
        
        for m in range(self.nmodel):
            # Get the indices for this model and extract
            # the relevant parameters 

            # create a mapping to the original position
            mapping = numpy.arange(self.nbatch)[mods == m]
            if self.debug == 2:print "\t\t\tmodel / mapping:", m , mapping
            n_to_simulate = len( mapping )
                    
            if( n_to_simulate > 0 ):
                # simulate this chunk
                this_model_parameters = []
                for i in range(n_to_simulate):
                    this_model_parameters.append( sampled_params[ mapping[i] ])

                #print "this_model_parameters", this_model_parameters
                sims = self.models[ m ].simulate( this_model_parameters, self.data.timepoints, n_to_simulate, self.beta )
                if self.debug == 2:print '\t\t\tsimulation dimensions:', sims.shape

                for i in range(n_to_simulate):
                    this_dist = []
                    for k in range(self.beta):
                        samplePoints = sims[i,k,:,:]
                        #print 'samplePoints', samplePoints
                        points = howToFitData( self.models[ m ].fit, samplePoints )
                        traj.append( points )
                        distance=self.distancefn(points, self.data.values, this_model_parameters[i], m)
                        this_dist.append( distance )

                        dist = evaluateDistance(distance, this_epsilon )
                        #if(distance <= this_epsilon and distance >= 0):
                        #    ret[mapping[i]] = ret[mapping[i]] + 1

                        if dist == True:
                            ret[mapping[i]] = ret[mapping[i]] + 1

                        if self.debug == 2:print '\t\t\tdistance/this_epsilon/mapping/b:', distance, this_epsilon, mapping[i], ret[mapping[i]]
                    
                    distances.append( this_dist )
        return ret[:], traj, distances
        

    def sampleTheModelFromPrior(self):
        ret = [0 for it in range(self.nbatch)]
        if self.nmodel > 1:
            for i in range(self.nbatch):
                ret[i] = statistics.w_choice( range(self.nmodel), self.modelprior )
        
        return ret[:]

    def sampleTheParameterFromPrior(self, sampled_models):
        ret = []
 
        for i in range(self.nbatch):
            #print "sampleTheParameterFromPrior", i, sampled_models[i], self.models[ sampled_models[i] ].name, self.models[ sampled_models[i] ].nparameters

            reti = [ 0 for it in range(self.models[ sampled_models[i] ].nparameters) ]

            for n in range(self.models[ sampled_models[i] ].nparameters):
                if self.models[ sampled_models[i] ].prior[n][0] == 0: 
                    reti[n]=self.models[ sampled_models[i] ].prior[n][1]
                
                if self.models[ sampled_models[i] ].prior[n][0] == 2: 
                    reti[n]=rnd.uniform( low=self.models[ sampled_models[i] ].prior[n][1],
                                         high=self.models[ sampled_models[i] ].prior[n][2])

                if self.models[ sampled_models[i] ].prior[n][0] == 3: 
                    reti[n]=rnd.lognormal(self.models[ sampled_models[i] ].prior[n][1],
                                          self.models[ sampled_models[i] ].prior[n][2])
            
            ret.append( reti[:] )
            
        return [x[:] for x in ret]


    def sampleTheModel(self):
        ret = [0 for it in range(self.nbatch)]
    
        if self.nmodel > 1:
            for i in range(self.nbatch):
                ret[i] = statistics.w_choice( range(self.nmodel), self.margins_prev )

            if( len(self.dead_models) > self.nmodel-1 ):
            #if(0): 
                # perturb models
                for i in range(self.nbatch):
                    u = rnd.uniform(low=0, high=1)
                    if u > self.modelKernel:
                        # sample randomly from other (non dead) models
                        not_available = self.dead_models[:]
                        not_available.append(ret[i])
                        ss = set( not_available )
                        s = set( range(0,self.nmodel) )
                        # print "perturbing model ", set( range(self.nmodel) ), ss, s-ss, list(s-ss), numpy.array( list(s-ss) )
                        ar = numpy.array( list(s-ss) )
                        rnd.shuffle( ar )
                        perturbed_model = ar[0]
                        # print "perturbation ", ret[i], "->", perturbed_model
                        ret[i] = perturbed_model
        return ret[:]
    
    def sampleTheParameter(self, sampled_models):
        if self.debug == 2:print "\t\t\t***sampleTheParameter"
        ret = []
        
        for i in range(self.nbatch):
            np = self.models[ sampled_models[i] ].nparameters
            reti = [ 0 for it in range(np) ]

            #print '\n\t\t\tsampleTheParameter, model np prior:', sampled_models[i], self.models[ sampled_models[i] ].name, np, self.models[ sampled_models[i] ].prior
            
            prior_prob = -1
            while prior_prob <= 0 :

                # sample putative particle from previous population
                p = sample_particle(self.nparticles, sampled_models[i], self.margins_prev, self.model_prev, self.weights_prev )
                
                for nn in range(np):
                    #print reti[nn], self.parameters_prev[ p ][nn]
                    reti[nn] = self.parameters_prev[ p ][nn]
                
                prior_prob = self.perturbfn( reti, self.models[ sampled_models[i] ].prior, self.kernels[sampled_models[i]] )

                if self.debug == 2:print "\t\t\tsampled p prob:", prior_prob
                if self.debug == 2:print "\t\t\tnew:", reti
                if self.debug == 2:print "\t\t\told:", self.parameters_prev[p]

            ret.append( reti )

        return [x[:] for x in ret]
    
    def computeParticleWeights(self):
        if self.debug == 2:print "\t***computeParticleWeights"

        for k in range(self.nparticles):
            this_model = self.model_curr[k]
            this_param = self.parameters_curr[k]
            
            # calculate model prior probility 
            mprob = self.modelprior[ this_model ]

            np = len(self.parameters_curr[k])
            # particle prior probability
            pprob = 1
            for n in range(0,np):
                x = 1.0
                if self.models[ this_model ].prior[n][0]==0:
                    x=1
                
                if self.models[ this_model ].prior[n][0]==2: 
                    x=statistics.getPdfUniform(self.models[ this_model ].prior[n][1],
                                               self.models[ this_model ].prior[n][2],
                                               this_param[n])
                
                if self.models[ this_model ].prior[n][0]==3: 
                    x=statistics.getPdfLognormal(self.models[ this_model ].prior[n][1],
                                                 self.models[ this_model ].prior[n][2],
                                                 this_param[n])
                pprob = pprob*x

            numer = self.b[k] * mprob * pprob
        
            denom_m = 0
            for i in range(self.nmodel):
                denom_m = denom_m + self.margins_prev[i]*getPdfModelKernel(this_model, i, self.modelKernel, self.nmodel, self.dead_models)
            denom = 0
            # print "Calculating denom\t", selected_model, sampleParameters
            for j in range(self.nparticles):

                if(int(this_model) == int(self.model_prev[j]) ):
                    # print "\t", j, model_prev[j], weights_prev[j], parameters_prev[j]
                    if self.debug == 2: print "\tj, weights_prev, kernelpdf", j, self.weights_prev[j], self.kernelpdffn(this_param, self.parameters_prev[j], self.models[this_model].prior, self.kernels[this_model] )
                    denom = denom + self.weights_prev[j] * self.kernelpdffn(this_param, self.parameters_prev[j], self.models[this_model].prior, self.kernels[this_model] )

                if self.debug == 2: print "\tnumer/denom_m/denom/m(t-1) : ", numer,denom_m, denom, self.margins_prev[this_model]

            self.weights_curr[k] = numer/(denom_m*denom/self.margins_prev[this_model])
        
    def normalizeWeights(self):
        n = sum( self.weights_curr )
        for i in range(self.nparticles):
            self.weights_curr[i] = self.weights_curr[i]/float(n)

    def modelMarginals(self):
        for i in range(self.nmodel):
            self.margins_curr[i] = 0
            for j in range(self.nparticles):
                if int(self.model_curr[j]) == int(i):
                    self.margins_curr[i] = self.margins_curr[i] + self.weights_curr[j]
            
    


# additional work functions

def sample_particle(nparticle, selected_model, margins_prev, model_prev, weights_prev ):
    u = rnd.uniform(low=0, high=margins_prev[selected_model])
    F = 0

    for i in range(0,nparticle):
        if int(model_prev[i]) == int(selected_model) :
            F = F + weights_prev[i]

            if(F > u):
                break

    return i

########### how to fit variable to the given data (if fit is None, data for all variables are available in order of the model)
def howToFitData(fitting_instruction,samplePoints):
    
    if fitting_instruction!=None:
        #print fitting_instruction, fitting_instruction[2]
        points=numpy.zeros([len(samplePoints), len(fitting_instruction)])
        #print "points", points.shape
        #print "len", len(fitting_instruction)
        #print "samples", samplePoints.shape

        for i in range(len(fitting_instruction)):
            #print 'eval', i, eval(fitting_instruction[i])
            points[:,i] = eval(fitting_instruction[i])
           
    else:
        points=samplePoints

    return points[:]

def getPdfModelKernel(m, m0, modelK, nmodel, dead_models):
    ndead = len(dead_models)
   
    if(ndead == nmodel-1):
        return 1.0
    else:

        if(m == m0):
            return modelK
        else:
            return (1-modelK)/(nmodel-ndead)

def evaluateDistance(distance,epsilon):
    #print 'd:', distance
    #print 'e:', epsilon

    accepted = False
    for i in range(len(epsilon)):
        #print "d:", distance[i], epsilon[i][t]
        if(distance[i]<epsilon[i] and distance[i]>=0 ):
            accepted = True
        else: 
            accepted = False
            break

    return accepted
