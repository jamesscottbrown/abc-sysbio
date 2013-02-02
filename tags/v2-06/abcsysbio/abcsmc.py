import numpy
from numpy import random as rnd

import copy, time

from abcsysbio import euclidian
from abcsysbio import kernels
from abcsysbio import statistics

"""
priors: 
              a 3D list.
              The first dimension represents the number of models, the second dimension
              represents the number of parameter for a particular model and the third dimension
              represents the distribution for this parameter and has a constant length of 3.
              The first entry for each parameter is an integer number that stands for a
              specific distribution.
              
              Implemented distributions are:
              0   ---   constant parameter.
                        Example: constant parameter with value 12.3
                        [0 , 12.3 , x] , where x can be any number

              1   ---   normal distribution. 
                        Example: normal distribution in mean 10 and var 1
                        [1 , 10 , 1]

              2   ---   uniform distribution. 
                        Example: uniform distribution in the range 0.1 to 50
                        [2 , 0.1 , 50]

              3   ---   lognormal distribution.
                        Example: lognormal distribution with mean 3 and var 1.5
                        [3 , 3 , 1.5]

              Example:
              1 model with 3 parameter, the first two parameter have uniform prior between 0 and
              5, the third parameter has lognormal prior with mean 1 and varianz 3.
              [ [ [1,0,5],[1,0,5],[2,1,3] ], ]
              2 models where the first model has 2 parameter (the first is constant 3 and the 
              second is uniform between 0 and 1) and the second model has 1 lognormal parameter
              with 0 mean and varianz 0.5
              [ [ [0,3,0],[1,0,1] ] , [ [2,0,0.5],] ]

fit:      
              a 2D list of strings.
              This list contains the fitting instructions and therefore defines how to fit the
              experimental data to the systems variables. The first dimension represents the 
              models to be investigated, the second dimension represents the amount of data.
              
              Example:
              1 model with 7 species, 3 data series
              fit = [ ['species1+species2' , 'species5' , '3*species4-species7'], ]
              2 models with each 7 species, 2 data series
              fit = [ ['species1' , 'species2*species3'] , ['species1' , 'species2*species4'] ]
              
              The default value is 'None', i.e. the order of data series corresponds exactly 
              to the order of species.
"""

## distances are stored as [nparticle][nbeta][d1, d2, d3 .... ]
## trajectories are stored as [nparticle][nbeta][ species ][ times ]

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
        self.trajectories = copy.deepcopy(trajectories)
        self.distances = copy.deepcopy(distances)
        self.margins = copy.deepcopy(margins)
        self.models = copy.deepcopy(models)
        self.weights = copy.deepcopy(weights)
        self.parameters = copy.deepcopy(parameters)
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
                 timing,
                 distancefn = euclidian.euclidianDistance,
                 kernel_type = 1,
                 kernelfn = kernels.getKernel,
                 kernelpdffn = kernels.getPdfParameterKernel,
                 perturbfn = kernels.perturbParticle):
        
        self.nmodel = len(models)
        self.models = copy.copy( models )
        self.data = copy.deepcopy(data)
        
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
        self.distances = []
        self.trajectories = []

        self.distancefn = distancefn
        self.kernel_type = kernel_type
        self.kernelfn = kernels.getKernel
        self.kernelpdffn = kernels.getPdfParameterKernel
        self.perturbfn = kernels.perturbParticle

        self.beta = beta
        self.dead_models = []
        self.nbatch = nbatch
        self.debug = debug
        self.timing = timing
    
        self.modelprior = modelprior[:]
        self.modelKernel = modelKernel
        self.kernel_aux = [0 for i in range(0,nparticles)]

        self.kernels = list()
        # self.kernels is a list of length the number of models
        # self.kernels[i] is a list of length 2 such that :
        # self.kernels[i][0] contains the index of the non constant parameters for the model i
        # self.kernels[i][1] contains the information required to build the kernel and given by the input_file
        # self.kernels[i][2] is filled in during the kernelfn step and contains values/matrix etc depending ont he kernel
        kernel_option=list()
        for i in range(self.nmodel):
            if self.kernel_type==4:
                # Option for K nearest neigbours - user should be able to specify
                kernel_option.append(int(nparticles/4))
            else:
                kernel_option.append(0)

        # get the list of parameters with non constant prior
        for i in range(self.nmodel):
            ind=list()
            for j in range(self.models[i].nparameters):
                if not(self.models[i].prior[j][0]==0):
                    ind.append(j)
            # kernel info will get set after first population
            self.kernels.append( [ind, kernel_option[i], 0 ] ) 

        # get
        self.special_cases = [0 for m in range(self.nmodel)]
        if self.kernel_type == 1:
            
            for m in range(self.nmodel):
                all_uniform = True
                for j in range(self.models[m].nparameters):
                    if not (self.models[m].prior[j][0]==0 or self.models[m].prior[j][0]==2):
                        all_uniform = False
                if all_uniform == True:
                    self.special_cases[m] = 1
                    print "### Found special kernel case 1 for model ", m, "###"
        
        self.hits = []
        self.sampled = []
        self.rate = []
        self.dead_models = []
        self.sample_from_prior = True

    def run_fixed_schedule(self, epsilon, io):
        all_start_time = time.time()
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
                print "\t model marginals                  :", self.margins_prev
                
                if(len(self.dead_models) > 0):
                    print "\t dead models                      :", self.dead_models
                if self.timing == True:
                    print "\t timing:                          :", end_time - start_time

        if self.timing == True:
            print "#### final time:", time.time() - all_start_time

        return
            
    def run_automated_schedule(self, final_epsilon, alpha, io):
        all_start_time = time.time()
        
        done = False
        final = False
        pop = 0
        epsilon = [1e10 for i in final_epsilon]

        while done == False:
            if final==True: done = True

            start_time = time.time()
            if(pop==0 and self.sample_from_prior==True): 
                results = self.iterate_one_population(epsilon, prior=True)
            else:
                results = self.iterate_one_population(epsilon, prior=False)
            end_time = time.time()

            io.write_pickled(self.nmodel, self.model_prev, self.weights_prev, self.parameters_prev, self.margins_prev, self.kernels)
            io.write_data(pop, results, end_time-start_time, self.models, self.data)

            final, epsilon = self.compute_next_epsilon(results, epsilon, final_epsilon, alpha)

            if self.debug == 1:
                print "### population ", pop+1
                print "\t sampling steps / acceptance rate :", self.sampled[pop], "/", self.rate[pop]
                print "\t model marginals                  :", self.margins_prev
                print "\t next epsilon                     :", epsilon
                
                if(len(self.dead_models) > 0):
                    print "\t dead models                      :", self.dead_models
                if self.timing == True:
                    print "\t timing:                          :", end_time - start_time

            pop += 1

        if self.timing == True:
            print "#### final time:", time.time() - all_start_time

        return


    def compute_next_epsilon(self, results, this_epsilon, target_epsilon, alpha):
        nepsilon = len(target_epsilon)
        #print "compute_next_epsilon : this_epsilon, target_epsilon, nepsilon:", this_epsilon, target_epsilon, nepsilon 

        distance_values = []
        for i in range(self.nparticles):
            for j in range(self.beta):
                distance_values.append( results.distances[i][j] )

        # Important to remember that the initial sort on distance is done on the first distance value
        distance_values = numpy.sort(distance_values, axis=0)
        ntar = int( alpha * self.nparticles )
        #print distance_values

        new_epsilon = [ round(distance_values[ntar,ne],4) for ne in range(nepsilon) ]
        ret_epsilon = [0 for i in range(nepsilon)]

        # Set the next epsilon
        for ne in range(nepsilon):
            if new_epsilon[ne] < this_epsilon[ne]:
                ret_epsilon[ne] = new_epsilon[ne]
            else :
                # This is an attempt to reduce epsilon even if the new and previous epsilon are equal
                ret_epsilon[ne] = 0.95*new_epsilon[ne]

        # print "new/ret epsilon:", new_epsilon, ret_epsilon

        # See if we are finished
        finished = True
        for ne in range(nepsilon):
            if ret_epsilon[ne] < target_epsilon[ne] or numpy.fabs(ret_epsilon[ne]-target_epsilon[ne]) < 1e-6:
                ret_epsilon[ne] = target_epsilon[ne]
            else:
                finished = False 

        return finished, ret_epsilon 

    def run_simulations(self, io):

        naccepted = 0
        sampled = 0

        while(naccepted < self.nparticles):
            if self.debug == 2:print "\t****batch"
            sampled_models = self.sampleTheModelFromPrior()
            sampled_params = self.sampleTheParameterFromPrior(sampled_models)
                
            accepted_index, distances, traj = self.simulate_and_compare_to_data(sampled_models,sampled_params,this_epsilon=0,do_comp=False)

            for i in range(self.nbatch):
                if naccepted < self.nparticles:
                    sampled = sampled + 1

                if naccepted < self.nparticles and accepted_index[i] > 0 : 
                    
                    self.model_curr[naccepted] = sampled_models[i]
                    if self.debug == 2:print "\t****accepted", i, accepted_index[i], sampled_models[i]
                    
                    for p in range( self.models[ sampled_models[i] ].nparameters ):
                        self.parameters_curr[naccepted].append(sampled_params[i][p])

                    self.b[naccepted] = accepted_index[i]
                    self.trajectories.append( copy.deepcopy(traj[i]) )
                    self.distances.append( copy.deepcopy(distances[i]) )
                    
                    naccepted = naccepted + 1

            if self.debug == 2:print "\t****end  batch naccepted/sampled:", naccepted,  sampled

        # Finished loop over particles
        if self.debug == 2:print "**** end of population naccepted/sampled:", naccepted,  sampled

        results = abcsmc_results(naccepted, sampled, naccepted/float(sampled), self.trajectories, self.distances, 
                                 0, self.model_curr, 0, self.parameters_curr, 0)

        io.write_data_simulation(0, results, 0, self.models, self.data)
        


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
                
            accepted_index, distances, traj = self.simulate_and_compare_to_data(sampled_models,sampled_params,next_epsilon)

            for i in range(self.nbatch):
                if naccepted < self.nparticles:
                    sampled = sampled + 1

                if naccepted < self.nparticles and accepted_index[i] > 0 : 
                    
                    self.model_curr[naccepted] = sampled_models[i]
                    if self.debug == 2:print "\t****accepted", i, accepted_index[i], sampled_models[i]
                    
                    for p in range( self.models[ sampled_models[i] ].nparameters ):
                        self.parameters_curr[naccepted].append(sampled_params[i][p])

                    self.b[naccepted] = accepted_index[i]
                    self.trajectories.append( copy.deepcopy(traj[i]) )
                    self.distances.append( copy.deepcopy(distances[i]) )
                    
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
            this_weights = numpy.zeros( len(this_model_index) ) 

            # if we have just sampled from the prior we shall initialise the kernels using all available particles
            if prior == True:
                for it in range(len(this_model_index)):
                    this_population[it,:] = self.parameters_prev[ this_model_index[it] ][:]
                    this_weights[it] = self.weights_prev[ this_model_index[it] ]
                tmp_kernel = self.kernelfn( self.kernel_type, self.kernels[mod], this_population, this_weights )
                self.kernels[mod] = tmp_kernel[:]

            else:
                # only update the kernels if there are > 5 particles
                if len(this_model_index) > 5:
                    for it in range(len(this_model_index)):
                        this_population[it,:] = self.parameters_prev[ this_model_index[it] ][:]
                        this_weights[it] = self.weights_prev[ this_model_index[it] ]
                    tmp_kernel = self.kernelfn( self.kernel_type, self.kernels[mod], this_population, this_weights )
                    self.kernels[mod] = tmp_kernel[:]

        #
        # Kernel auxilliary information
        #
        self.kernel_aux = kernels.getAuxilliaryInfo(self.kernel_type, self.model_prev, self.parameters_prev, self.models, self.kernels )[:]
        
        self.hits.append( naccepted )
        self.sampled.append( sampled )
        self.rate.append( naccepted/float(sampled) )

        results = abcsmc_results(naccepted, 
                                 sampled,
                                 naccepted/float(sampled), 
                                 self.trajectories, 
                                 self.distances, 
                                 self.margins_prev, 
                                 self.model_prev, 
                                 self.weights_prev, 
                                 self.parameters_prev,
                                 next_epsilon)

        self.trajectories = []
        self.distances = []

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
            for j in range( len(particle_data[4][i])):
                self.kernels[i].append( particle_data[4][i][j] )
        
        self.sample_from_prior = False

    def simulate_and_compare_to_data(self, sampled_models, sampled_params, this_epsilon, do_comp=True):
        # Here do the simulations for each model together
        if self.debug == 2:print '\t\t\t***simulate_and_compare_to_data'

        ret = [0 for it in range(self.nbatch)]
        traj = [[] for it in range(self.nbatch) ]
        distances = [[] for it in range(self.nbatch)]
    
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
                    # store the trajectories and distances in a list of length beta
                    this_dist = []
                    this_traj = []
                    
                    for k in range(self.beta):
                        samplePoints = sims[i,k,:,:]
                        points = howToFitData( self.models[ m ].fit, samplePoints )
                        if do_comp == True:
                            distance = self.distancefn(points, self.data.values, this_model_parameters[i], m)
                            dist = evaluateDistance(distance, this_epsilon )
                        else:
                            distance = 0
                            dist = True

                        this_dist.append( distance )
                        this_traj.append( points )
                       
                        if dist == True:
                            ret[mapping[i]] = ret[mapping[i]] + 1
                            
                        if self.debug == 2:print '\t\t\tdistance/this_epsilon/mapping/b:', distance, this_epsilon, mapping[i], ret[mapping[i]]

                    traj[ mapping[i] ] = copy.deepcopy( this_traj )
                    distances[ mapping[i] ] = copy.deepcopy( this_dist )
                
        return ret[:], distances, traj    

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

                if self.models[ sampled_models[i] ].prior[n][0] == 1: 
                    reti[n]=rnd.normal( loc=self.models[ sampled_models[i] ].prior[n][1],
                                        scale=numpy.sqrt(self.models[ sampled_models[i] ].prior[n][2]) )
                
                if self.models[ sampled_models[i] ].prior[n][0] == 2: 
                    reti[n]=rnd.uniform( low=self.models[ sampled_models[i] ].prior[n][1],
                                         high=self.models[ sampled_models[i] ].prior[n][2])

                if self.models[ sampled_models[i] ].prior[n][0] == 3: 
                    reti[n]=rnd.lognormal(mean=self.models[ sampled_models[i] ].prior[n][1],
                                          sigma=numpy.sqrt(self.models[ sampled_models[i] ].prior[n][2]) )
            
            ret.append( reti[:] )
            
        return [x[:] for x in ret]


    def sampleTheModel(self):
        ret = [0 for it in range(self.nbatch)]
    
        if self.nmodel > 1:
            for i in range(self.nbatch):
                ret[i] = statistics.w_choice( range(self.nmodel), self.margins_prev )

            if( len(self.dead_models) < self.nmodel-1 ):
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
                
                prior_prob = self.perturbfn( reti, self.models[ sampled_models[i] ].prior, self.kernels[sampled_models[i]], self.kernel_type, self.special_cases[sampled_models[i]] )

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

                if self.models[ this_model ].prior[n][0]==1: 
                    x=statistics.getPdfGauss(self.models[ this_model ].prior[n][1],
                                             numpy.sqrt(self.models[ this_model ].prior[n][2]),
                                             this_param[n])
                
                if self.models[ this_model ].prior[n][0]==2: 
                    x=statistics.getPdfUniform(self.models[ this_model ].prior[n][1],
                                               self.models[ this_model ].prior[n][2],
                                               this_param[n])
    
                if self.models[ this_model ].prior[n][0]==3: 
                    x=statistics.getPdfLognormal(self.models[ this_model ].prior[n][1],
                                                 numpy.sqrt(self.models[ this_model ].prior[n][2]),
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
                    if self.debug == 2:
                        print "\tj, weights_prev, kernelpdf", j, self.weights_prev[j],
                        self.kernelpdffn(this_param, self.parameters_prev[j], self.models[this_model].prior, self.kernels[this_model], self.kernel_aux[j], self.kernel_type )
                    denom = denom + self.weights_prev[j] * self.kernelpdffn(this_param, self.parameters_prev[j], self.models[this_model].prior, self.kernels[this_model], self.kernel_aux[j], self.kernel_type )

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
    
    accepted = False
    for i in range(len(epsilon)):
        #print "d:", distance[i], epsilon[i][t]
        if(distance[i]<epsilon[i] and distance[i]>=0 ):
            accepted = True
        else: 
            accepted = False
            break

    #if accepted == True:
    #    print '\teval:', distance, epsilon
        
    return accepted
