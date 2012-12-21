### ABCSMC.py ###
# This script contains the abcsmc class constructor and associated functions. 

import numpy  
from numpy import random as rnd

import copy, time
import sys 

from abcsysbio import euclidian
from abcsysbio import kernels
from abcsysbio import statistics 
cstats = statistics.c_statistics()  # instantiate statistics class, which contains statistics functions created with ctypes


"""
PREAMBLE 
How Key Information is Stored

priors: 
              a 3D list.
              The first dimension represents the number of models, the second dimension
              represents the number of parameters for a particular model, and the third dimension, which has length=3,
              represents the distribution of this parameter.
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

              4   ---   truncated normal distribution
                        Example: truncated normal with mean 10 and var 1, with an upper bound of 2 (infinity is set to +/- 10,000,001 - change in parse_info.py)
                        [4, 10, 1, neg_infinity, 2]
                        Example: truncated normal with mean 10 and var 1, bounded at 2 and 5
                        [4, 10, 1, 2, 5]

              Examples:
              1 model with 3 parameters, the first two parameters have uniform priors between 0 and
              5, and the third parameter has a lognormal prior with mean 1 and variance 3.
              [ [ [2,0,5],[2,0,5],[3,1,3] ], ]
              
              2 models where the first model has 2 parameters (the first is constant 3 and the 
              second is uniform between 0 and 1) and the second model has 1 lognormal parameter
              with 0 mean and variance 0.5
              [ [ [0,3,0],[2,0,1] ] , [ [3,0,0.5],] ]

fit:      
              a 2D list of strings.
              This list contains the fitting instructions and therefore defines how to fit the
              experimental data to the systems' variables. The first dimension represents the 
              models to be investigated, and the second dimension represents the amount of data.
              
              Example:
              1 model with 7 species, 3 data series
              fit = [ ['species1+species2' , 'species5' , '3*species4-species7'], ]
              2 models with each 7 species, 2 data series
              fit = [ ['species1' , 'species2*species3'] , ['species1' , 'species2*species4'] ]
              
              The default value is 'None', i.e. the order of data series corresponds exactly 
              to the order of species.

kernel_info:
                a list with length = no of models.  Contains the appropriate information for parameter sampling.
                kernel_info[i] is a list of length 4:
                    kernel_info[i][0] is a list of indices of non-constant parameters.
                    kernel_info[i][1] contains a list of upper and lower bounds for each parameter, which are defined by the parameter's prior
                                        and are used for perturbing particles with truncated distributions. 
                    kernel_info[i][2] kernel info defined in the input file.
                    kernel_info[i][3] contains the kernel once it has been built (using getXKinfo).

kernel_functions:
                a list with length = no of models.  Contains the appropriate kernel functions for each model/parameter
                kernel_functions[i] is a list of length 3:
                    kernel_functions[i][0] contains just the function to update the kernel after every population - getKernelInfo
                    kernel_functions[i][1] is a list (length = no of nonconstant parameters for the model) that contains a sampling 
                                        function for each nonconstant parameter (eg. sample from truncated/regular uniform distribution)
                    kernel_functions[i][1] is a list (length = no of nonconstant parameters for the model) that contains a PDF for     
                                        each nonconstant parameter (eg. PDF of a truncated/regular uniform distribution)  

distances are stored as [nparticle][nbeta][d1, d2, d3 .... ]
trajectories are stored as [nparticle][nbeta][ species ][ times ]
"""


##############################################################################################################################
                                          ### ABCSMC_RESULTS CLASS   ###
##############################################################################################################################
class abcsmc_results:
    def __init__(self, 
                 naccepted, 
                 sampled,
                 rate, 
                 trajectories, 
                 distances,
                 alldistances, 
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
        self.alldistances = copy.deepcopy(alldistances)
        self.margins = copy.deepcopy(margins)
        self.models = copy.deepcopy(models)
        self.weights = copy.deepcopy(weights)
        self.parameters = copy.deepcopy(parameters)
        self.epsilon = epsilon




##############################################################################################################################
                                          ### ABCSMC CLASS   ###
##############################################################################################################################
class abcsmc:

    # instantiation
    def __init__(self, 
                 models,                                            # models defined in input
                 nparticles,                                        # no of (accepted) populations per population
                 modelprior,                                        # prior probabilities of each model
                 data,                                              # dataset(s) associated with models
                 beta,                                              # no of times to simulate each sampled particle (1 for ODEs)
                 nbatch,                                            # no of particles to sample at a time; ie. sample nbatch, simulate for nbatch, check distances
                 modelKernel,                                       # kernel for perturbing models (NOT parameters)
                 debug,                                             # debugging option; ie. printing information and errors to standard output
                 timing,                                            # print times to standard output
                 distancefn,                                        # functions to calculate distance between simulated and data
                 kernel_type):                                      # kernel for perturbing parameters

        self.nmodel = len(models)
        self.nsubmodel = len(models[0].submodelname)    
        
        self.models = copy.copy( models )
        self.data = copy.deepcopy(data)
        
        self.nparticles = nparticles

        self.model_prev      = [0  for i in range(nparticles)]      # list to store models from which particles were sampled and accepted in previous population
        self.weights_prev    = [0  for i in range(nparticles)]      # list to store weights of particles that were sampled and accepted in previous population
        self.parameters_prev = [[] for i in range(nparticles)]      # nested list for parameter values of particles in previous population
        self.margins_prev    = [0  for i in range(self.nmodel)]     # marginal probabilities of models in previous population
        
        self.model_curr      = [0  for i in range(nparticles)]
        self.weights_curr    = [0  for i in range(nparticles)]
        self.parameters_curr = [[] for i in range(nparticles)]
        self.margins_curr    = [0  for i in range(self.nmodel)] 
        
        self.b = [0  for i in range(0,nparticles)]                  # list to store accepted indices 
        self.alldistances = []                                      # distances between simulated results and actual data
        self.distances = []
        self.trajectories = []

        self.distancefn = distancefn
        self.submodelweights = [1.0/self.nsubmodel for j in range(self.nsubmodel)] 
        self.kernel_type = kernel_type

        self.beta = beta
        self.dead_models = []                                       # list to store models with very low probabilities
        self.nbatch = nbatch
        self.debug = debug
        self.timing = timing
    
        self.modelprior = modelprior[:]
        self.modelKernel = modelKernel
        self.kernel_aux = [0 for i in range(0,nparticles)]

        self.kernel_functions = list()  
        self.kernel_info = list()  
        self.perturbfn=list() 
        self.getPdfParameter = list() 
        
        self.hits = []
        self.sampled = []
        self.rate = []
        self.sample_from_prior = True


        # If custom kernel is to be used, use a parameter sampling function that checks whether perturbed particle values are within prior bounds
        # Else, sampling function does not perform this check, as truncated distributions will definitely be used
        if self.kernel_type[0] == 6:
            self.sampleTheParameter = self.sampleTheParameterCustom
        else:
            self.sampleTheParameter = self.sampleTheParameterNonCustom

        
        #Infinity values are used to represent infinity in practice.  Required for MVN
        self.positive_infinity = 10000001
        self.negative_infinity = -10000001


        # Go through each parameter of each model.  Extract the relevant parameter information (indices of nonconstand parameters, prior bounds for each parameter)
        # and assign the appropriate kernel functions depending on kernel type and prior distribution.
            
        if self.kernel_type[0] == 1:  # if Uniform CW kernel
            for i in range(self.nmodel):
                self.perturbfn.append(kernels.perturbParticleCW)
                self.getPdfParameter.append(kernels.getPdfParameterCW)
                nc_indices=[]
                sample_function=[]
                pdf_function=[]
                bounds=[]
                for j in range(self.models[i].nparameters):
                    if self.models[i].prior[j][0]==1:     # normal prior
                        nc_indices.append(j)
                        bounds.append( [self.negative_infinity, self.positive_infinity] )
                        sample_function.append(cstats.sampleUniformReg)
                        pdf_function.append(cstats.pdfUniformReg)
                       
                    elif self.models[i].prior[j][0]==2:    # uniform prior   
                        nc_indices.append(j)
                        bounds.append([ self.models[i].prior[j][1], self.models[i].prior[j][2] ])
                        sample_function.append(cstats.sampleUniformTrunc)
                        pdf_function.append(cstats.pdfUniformTrunc) 
                        
                    elif self.models[i].prior[j][0]==3:     # lognormal prior
                        nc_indices.append(j)
                        bounds.append( [0, self.positive_infinity] )
                        sample_function.append(cstats.sampleUniformTrunc)
                        pdf_function.append(cstats.pdfUniformTrunc) 

                    elif self.models[i].prior[j][0]==4:     # truncated normal prior   
                        nc_indices.append(j)
                        bounds.append([ self.models[i].prior[j][3][0], self.models[i].prior[j][3][1] ])
                        sample_function.append(cstats.sampleUniformTrunc)
                        pdf_function.append(cstats.pdfUniformTrunc)
                             
                self.kernel_functions.append( [kernels.getUniformKinfo, sample_function, pdf_function] )
                self.kernel_info.append( [nc_indices, bounds, 0, 0] )


        elif self.kernel_type[0] == 2:  # if Normal CW kernel
            for i in range(self.nmodel):
                self.getPdfParameter.append(kernels.getPdfParameterCW)
                self.perturbfn.append(kernels.perturbParticleCW)
                nc_indices=[]
                sample_function=[]
                pdf_function=[]
                bounds=[]
                for j in range(self.models[i].nparameters):
                    if self.models[i].prior[j][0]==1:    # normal prior
                        nc_indices.append(j)
                        bounds.append( [self.negative_infinity, self.positive_infinity] )
                        sample_function.append(cstats.sampleGaussReg)
                        pdf_function.append(cstats.pdfGaussReg)
                       
                    elif self.models[i].prior[j][0]==2:    # uniform prior      
                        nc_indices.append(j)
                        bounds.append( [self.models[i].prior[j][1], self.models[i].prior[j][2]] )
                        sample_function.append(cstats.sampleGaussTrunc)
                        pdf_function.append(cstats.pdfGaussTrunc)
                        
                    elif self.models[i].prior[j][0]==3:     # lognormal prior
                        nc_indices.append(j)
                        bounds.append( [0, self.positive_infinity] )
                        sample_function.append(cstats.sampleGaussTrunc)
                        pdf_function.append(cstats.pdfGaussTrunc)

                    elif self.models[i].prior[j][0]==4:     # truncated normal prior
                        nc_indices.append(j)
                        bounds.append([ self.models[i].prior[j][3][0], self.models[i].prior[j][3][1] ])
                        sample_function.append(cstats.sampleGaussTrunc)
                        pdf_function.append(cstats.pdfGaussTrunc)
                             
                self.kernel_functions.append( [kernels.getGaussKinfo, sample_function, pdf_function] )
                self.kernel_info.append( [nc_indices, bounds, 0, 0] )


        elif self.kernel_type[0] > 2 and self.kernel_type[0] < 6:   # if Multivariate Normal kernels
            for i in range(self.nmodel):
                nc_indices = []
                bounds=[]
                truncated_present = 0
                for j in range(self.models[i].nparameters):
                    if self.models[i].prior[j][0] == 1: # normal prior
                        nc_indices.append(j)
                        bounds.append([self.negative_infinity, self.positive_infinity])
                        truncated_present = 0
                    elif self.models[i].prior[j][0] == 2: # uniform prior
                        nc_indices.append(j)
                        bounds.append([self.models[i].prior[j][1], self.models[i].prior[j][2]])
                        truncated_present = 1
                    elif self.models[i].prior[j][0] == 3: # lognormal prior
                        nc_indices.append(j)
                        bounds.append([0, self.positive_infinity])
                        truncated_present = 1
                    elif self.models[i].prior[j][0] == 4: # truncated normal prior
                        nc_indices.append(j)
                        bounds.append([ self.models[i].prior[j][3][0], self.models[i].prior[j][3][1] ])
                        truncated_present = 1
                
                self.kernel_info.append([nc_indices, bounds, nparticles, 0])
                               
                if self.kernel_type[0] == 3:                # specific information for MVN kernel variations
                    self.kernel_functions.append([kernels.getMVNkernel,0,0])
                    if truncated_present == 0:
                        self.perturbfn.append(kernels.perturbParticleMVN)
                        self.getPdfParameter.append(kernels.getPdfParameterKernelMVN)
                    else:
                        self.perturbfn.append(kernels.perturbParticleMVNT)
                        self.getPdfParameter.append(kernels.getPdfParameterKernelMVNT)
                elif self.kernel_type[0] == 4:
                    self.kernel_functions.append([kernels.getMVNkernelNN,0,0])
                    if truncated_present == 0:
                        self.perturbfn.append(kernels.perturbParticleMVN_local)
                        self.getPdfParameter.append(kernels.getPdfParameterKernelMVN_local)
                    else:
                        self.perturbfn.append(kernels.perturbParticleMVNT_local)
                        self.getPdfParameter.append(kernels.getPdfParameterKernelMVNT_local)
                elif self.kernel_type[0] == 5:
                    self.kernel_functions.append([kernels.getMVNkernelOpt,0,0])
                    if truncated_present == 0:
                        self.perturbfn.append(kernels.perturbParticleMVN_local)
                        self.getPdfParameter.append(kernels.getPdfParameterKernelMVN_local)
                    else:
                        self.perturbfn.append(kernels.perturbParticleMVNT_local)
                        self.getPdfParameter.append(kernels.getPdfParameterKernelMVNT_local)  

        if self.kernel_type[0] == 4:        # MVN - K-nearest neighbours
            for i in range(self.nmodel):
                k_choice = int(nparticles/4)     # option for K nearest neigbours - user should be able to specify
                self.kernel_info[i][2] == int(k_choice)
                self.kernel_info[i][3] == numpy.identity(k_choice).tolist()


#To use custom kernels, customABC must contain the following functions:
#   1. Get Kinfo
#   2. PerturbParticle - Takes in a particle, returns particle with altered values and prior probability
#   3. GetPDF - Takes in a particle, returns particle probability
        elif self.kernel_type[0] == 6:     #Custom kernel  
            import customABC
            for i in range(self.nmodel):
                self.getPdfParameter.append(customABC.getPdf)
                self.perturbfn.append(customABC.perturbParticle)
                nc_indices=[]
                bounds=[]
                for j in range(self.models[i].nparameters):
                    if self.models[i].prior[j][0]==1:    # normal prior
                        nc_indices.append(j)
                        bounds.append( [self.negative_infinity, self.positive_infinity] )
                       
                    elif self.models[i].prior[j][0]==2:    # uniform prior      
                        nc_indices.append(j)
                        bounds.append( [self.models[i].prior[j][1], self.models[i].prior[j][2]] )
                        
                    elif self.models[i].prior[j][0]==3:     # lognormal prior
                        nc_indices.append(j)
                        bounds.append( [0, self.positive_infinity] )

                    elif self.models[i].prior[j][0]==4:     # truncated normal prior
                        nc_indices.append(j)
                        bounds.append([ self.models[i].prior[j][3][0], self.models[i].prior[j][3][1] ])

                self.kernel_functions.append( [customABC.getKinfo, 0,0] )
                self.kernel_info.append( [nc_indices, bounds, 0, 0] )


        elif self.kernel_type[0] == 7:  # if non-adaptive Uniform CW kernel
            for i in range(self.nmodel):
                self.getPdfParameter.append(kernels.getPdfParameterCW)
                self.perturbfn.append(kernels.perturbParticleCW)
                nc_indices=[]
                sample_function=[]
                pdf_function=[]
                bounds=[]
                kernel_params=[]
                for j in range(self.models[i].nparameters):
                    if self.models[i].prior[j][0]==1:    # normal prior
                        nc_indices.append(j)
                        bounds.append( [self.negative_infinity, self.positive_infinity] )
                        sample_function.append(cstats.sampleUniformReg)
                        pdf_function.append(cstats.pdfUniformReg)
                        kernel_params.append([self.kernel_type[1][0], self.kernel_type[1][1]])
                       
                    elif self.models[i].prior[j][0]==2:    # uniform prior      
                        nc_indices.append(j)
                        bounds.append( [self.models[i].prior[j][1], self.models[i].prior[j][2]] )
                        sample_function.append(cstats.sampleUniformTrunc)
                        pdf_function.append(cstats.pdfUniformTrunc)
                        kernel_params.append([self.kernel_type[1][0], self.kernel_type[1][1]])
                        
                    elif self.models[i].prior[j][0]==3:     # lognormal prior
                        nc_indices.append(j)
                        bounds.append( [0, self.positive_infinity] )
                        sample_function.append(cstats.sampleUniformTrunc)
                        pdf_function.append(cstats.pdfUniformTrunc)
                        kernel_params.append([self.kernel_type[1][0], self.kernel_type[1][1]])

                    elif self.models[i].prior[j][0]==4:     # truncated normal prior
                        nc_indices.append(j)
                        bounds.append([ self.models[i].prior[j][3][0], self.models[i].prior[j][3][1] ])
                        sample_function.append(cstats.sampleUniformTrunc)
                        pdf_function.append(cstats.pdfUniformTrunc)
                        kernel_params.append([self.kernel_type[1][0], self.kernel_type[1][1]])
                             
                self.kernel_functions.append( [0, sample_function, pdf_function] )      
                self.kernel_info.append( [nc_indices, bounds, [0 for i in range(len(nc_indices))], kernel_params] )     # hardcode perturbation parameters


        elif self.kernel_type[0] == 8:  # if non-adaptive Normal CW kernel
            for i in range(self.nmodel):
                self.getPdfParameter.append(kernels.getPdfParameterCW)
                self.perturbfn.append(kernels.perturbParticleCW)
                nc_indices=[]
                sample_function=[]
                pdf_function=[]
                bounds=[]
                kernel_params=[]
                for j in range(self.models[i].nparameters):
                    if self.models[i].prior[j][0]==1:    # normal prior
                        nc_indices.append(j)
                        bounds.append( [self.negative_infinity, self.positive_infinity] )
                        sample_function.append(cstats.sampleGaussReg)
                        pdf_function.append(cstats.pdfGaussReg)
                        kernel_params.append(self.kernel_type[1])
                       
                    elif self.models[i].prior[j][0]==2:    # uniform prior      
                        nc_indices.append(j)
                        bounds.append( [self.models[i].prior[j][1], self.models[i].prior[j][2]] )
                        sample_function.append(cstats.sampleGaussTrunc)
                        pdf_function.append(cstats.pdfGaussTrunc)
                        kernel_params.append(self.kernel_type[1])
                        
                    elif self.models[i].prior[j][0]==3:     # lognormal prior
                        nc_indices.append(j)
                        bounds.append( [0, self.positive_infinity] )
                        sample_function.append(cstats.sampleGaussTrunc)
                        pdf_function.append(cstats.pdfGaussTrunc)
                        kernel_params.append(self.kernel_type[1])

                    elif self.models[i].prior[j][0]==4:     # truncated Normal Prior
                        nc_indices.append(j)
                        bounds.append([ self.models[i].prior[j][3][0], self.models[i].prior[j][3][1] ])
                        sample_function.append(cstats.sampleGaussTrunc)
                        pdf_function.append(cstats.pdfGaussTrunc)
                        kernel_params.append(self.kernel_type[1])
                             
                self.kernel_functions.append( [0, sample_function, pdf_function] )
                self.kernel_info.append( [nc_indices, bounds, kernel_params,  [0 for i in range(len(nc_indices))]] )    # hardcode perturbation parameters

                        
### CLASS FUNCTION: run_fixed_schedule
### Run the ABC SMC algorithm for a user-defined series of epsilons. 
##############################################################################################################################
    def run_fixed_schedule(self, epsilon, io):
        all_start_time = time.time()
        for pop in range(len(epsilon)):
            start_time = time.time()
            if(pop==0 and self.sample_from_prior==True): 
                results = self.iterate_one_population(epsilon[pop], prior=True)
            else:
                results = self.iterate_one_population(epsilon[pop], prior=False)
            end_time = time.time()
   
            io.write_pickled(self.nmodel, self.model_prev, self.weights_prev, self.parameters_prev, self.margins_prev, self.kernel_info) 
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

            
### CLASS FUNCTION: run_automated_schedule
### Run the ABC SMC algorithm for an epsilon series that is generated automatically (using alpha)        
##############################################################################################################################     
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

            io.write_pickled(self.nmodel, self.model_prev, self.weights_prev, self.parameters_prev, self.margins_prev, self.kernel_info)
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


### CLASS FUNCTION: compute_next_epsilon
### Compute the next epsilon for run_automated_schedule
##############################################################################################################################
    def compute_next_epsilon(self, results, this_epsilon, target_epsilon, alpha):
        nepsilon = len(target_epsilon)
 
        distance_values = []
        for i in range(self.nparticles):
            for j in range(self.beta):
                distance_values.append( results.distances[i][j] )

        # important to remember that the initial sort on distance is done on the first distance value
        distance_values = numpy.sort(distance_values, axis=0)
        ntar = int( alpha * self.nparticles )


        new_epsilon = [ round(distance_values[ntar,ne],4) for ne in range(nepsilon) ]
        ret_epsilon = [0 for i in range(nepsilon)]

        # set the next epsilon
        for ne in range(nepsilon):
            if new_epsilon[ne] < this_epsilon[ne]:
                ret_epsilon[ne] = new_epsilon[ne]
            else :
                # this is an attempt to reduce epsilon even if the new and previous epsilon are equal
                ret_epsilon[ne] = 0.95*new_epsilon[ne]



        # see if we are finished
        finished = True
        for ne in range(nepsilon):
            if ret_epsilon[ne] < target_epsilon[ne] or numpy.fabs(ret_epsilon[ne]-target_epsilon[ne]) < 1e-6:
                ret_epsilon[ne] = target_epsilon[ne]
            else:
                finished = False 

        return finished, ret_epsilon 


### CLASS FUNCTION: run_simulations
### Run simulations for simulation mode        
##############################################################################################################################   
    def run_simulations(self, io):

        naccepted = 0
        sampled = 0

        while(naccepted < self.nparticles):
            if self.debug == 2:print "\t****batch"
            sampled_models = self.sampleTheModelFromPrior()
            sampled_params = self.sampleTheParameterFromPrior(sampled_models)

            # distances and trajectories from here have an extra dimension for submodels                
            accepted_index, distances, distanceavgs, traj = self.simulate_and_compare_to_data(sampled_models,sampled_params,this_epsilon=0,do_comp=False) 

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
                    self.alldistances.append( copy.deepcopy(distances[i]) )
                    self.distances.append( copy.deepcopy(distanceavgs[i]) ) 
                    
                    naccepted = naccepted + 1

            if self.debug == 2:print "\t****end  batch naccepted/sampled:", naccepted,  sampled
        
        # finished loop over particles
        if self.debug == 2:print "**** end of population naccepted/sampled:", naccepted,  sampled

        results = abcsmc_results(naccepted, sampled, naccepted/float(sampled), self.trajectories, self.distances, self.alldistances, 
                                 0, self.model_curr, 0, self.parameters_curr, 0)

        io.write_data_simulation(0, results, 0, self.models, self.data)
        
         
### CLASS FUNCTION: iterate_one_population
### Sample models, parameters, simulate and check distances, compute weights, kernels and model marginals.       
##############################################################################################################################       
    def iterate_one_population(self, next_epsilon, prior):
        print "beginning population with", next_epsilon
        if self.debug == 2:print "\n\n****iterate_one_population: next_epsilon, prior", next_epsilon, prior

        naccepted = 0
        sampled = 0

        while(naccepted < self.nparticles):
            print "hits:", naccepted
            if self.debug == 2:print "\t****batch"
            if( prior == False):
                sampled_models = self.sampleTheModel()
                sampled_params = self.sampleTheParameter(sampled_models) 
            else:
                sampled_models = self.sampleTheModelFromPrior()
                sampled_params = self.sampleTheParameterFromPrior(sampled_models) 

            # traj and alldistances have dimensions that account for submodels
            # accepted_index counts the number of accepted submodels.              
            accepted_index, distances, distanceavgs, traj = self.simulate_and_compare_to_data(sampled_models,sampled_params,next_epsilon) 

            for i in range(self.nbatch):
                if naccepted < self.nparticles:
                    sampled = sampled + 1

                if naccepted < self.nparticles and accepted_index[i] > 0 : # accepted only if all submodels are, on average, accepted. 
                    
                    self.model_curr[naccepted] = sampled_models[i]
                    if self.debug == 2:print "\t****accepted", i, accepted_index[i], sampled_models[i]
                   
                    for p in range( self.models[ sampled_models[i] ].nparameters ): 
                        self.parameters_curr[naccepted].append(sampled_params[i][p])

                    self.b[naccepted] = accepted_index[i] 
                    self.trajectories.append(copy.deepcopy(traj[i])  )
                    self.alldistances.append(copy.deepcopy(distances[i])  ) 
                    self.distances.append(copy.deepcopy(distanceavgs[i]) ) 
                    
                    naccepted = naccepted + 1

            if self.debug == 2:print "\t****end  batch naccepted/sampled:", naccepted,  sampled

        # finished loop over particles
        if self.debug == 2:print "**** end of population naccepted/sampled:", naccepted,  sampled


        if( prior == False):
            self.computeParticleWeights()
        else:
            for i in range(self.nparticles):
                self.weights_curr[i] = self.b[i]
        
        self.normaliseWeights()
        self.modelMarginals()
  
        if self.debug == 2:
            print "**** end of population: particles"
            for i in range(self.nparticles):
                print i, self.weights_curr[i], self.model_curr[i], self.parameters_curr[i]
            print self.margins_curr

        # prepare for next population
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


        # check for dead models
        self.dead_models = []
        for j in range(self.nmodel):
            if self.margins_prev[j] < 1e-6:
                self.dead_models.append(j)

        
        # compute kernels
        if self.kernel_type[0] < 7:                # ie. only if adaptive kernels (or if custom kernel)
            for mod in range(self.nmodel):
                this_model_index = numpy.arange(self.nparticles)[ numpy.array(self.model_prev) == mod ]
                this_population = numpy.zeros([ len(this_model_index), self.models[mod].nparameters ])
                this_weights = numpy.zeros( len(this_model_index) ) 
            
                # if we have just sampled from the prior we shall initialise the kernels using all available particles
                if prior == True:
                    for it in range(len(this_model_index)):
                        this_population[it,:] = self.parameters_prev[ this_model_index[it] ][:] 
                        this_weights[it] = self.weights_prev[ this_model_index[it] ]
                    tmp_kernel = self.kernel_functions[mod][0]( self.kernel_info[mod], this_population, this_weights ) 
                    self.kernel_info[mod] = tmp_kernel[:]   

                else:
                    # only update the kernels if there are > 100 particles
                    if len(this_model_index) > 100:
                        for it in range(len(this_model_index)):
                            this_population[it,:] = self.parameters_prev[ this_model_index[it] ][:] 
                            this_weights[it] = self.weights_prev[ this_model_index[it] ]
                        tmp_kernel = self.kernel_functions[mod][0]( self.kernel_info[mod], this_population, this_weights )  
                        self.kernel_info[mod] = tmp_kernel[:]   

        
        self.hits.append( naccepted )
        self.sampled.append( sampled )
        self.rate.append( naccepted/float(sampled) )
        
        results = abcsmc_results(naccepted, 
                                 sampled,
                                 naccepted/float(sampled), 
                                 self.trajectories, 
                                 self.distances, 
                                 self.alldistances, 
                                 self.margins_prev, 
                                 self.model_prev, 
                                 self.weights_prev, 
                                 self.parameters_prev,
                                 next_epsilon)

        self.trajectories = []
        self.distances = []
        self.alldistances = []

        return results

    # end of iterate_one_population
       
  
### CLASS FUNCTION: fill_values
### Retrieves partial run data from a specified directory, and continues algorithm using this information (only for fixed epsilon schedule)  
##############################################################################################################################   
    def fill_values(self, particle_data ):
        # particle data comprises of [model_pickled, weights_pickled, parameters_pickled, margins_pickled, kernel]
        self.model_prev = particle_data[0][:]
        self.weights_prev = particle_data[1][:]
        self.margins_prev = particle_data[3][:]

        self.parameters_prev = []
        for it in range(len(particle_data[2])):
            self.parameters_prev.append( particle_data[2][it][:] )

        for i in range(self.nmodel):
            self.kernel_info[i][2] = particle_data[4][i]
            self.kernel_info[i][3] = particle_data[5][i]

        self.sample_from_prior = False

        
### CLASS FUNCTION: simulate_and_compare_to_data
### Use accepted particles to run model simulations, calculate distances and return accepted indices.    
##############################################################################################################################        
    def simulate_and_compare_to_data(self, sampled_models, sampled_params, this_epsilon, do_comp=True): 
        # Here do the simulations for each model together
        if self.debug == 2:print '\t\t\t***simulate_and_compare_to_data'

        ret = [0 for it in range(self.nbatch)]
        traj = [[[] for n in range(self.nsubmodel)] for it in range(self.nbatch) ] 
        distances = [[0 for n in range(self.nsubmodel)] for it in range(self.nbatch) ] 
        distanceavgs = [[] for it in range(self.nbatch) ] 
            
        n_to_simulate = [0 for it in range(self.nmodel)]
        mods = numpy.array( sampled_models )
        
        submodelweights = [1.0/self.nsubmodel for j in range(self.nsubmodel)]
        contribution = [0.0 for j in range(self.nsubmodel)] 

        # get the indices for this model and extract the relevant parameters
        for m in range(self.nmodel):
            mapping = numpy.arange(self.nbatch)[mods == m]      # create a mapping to the original position
            if self.debug == 2:print "\t\t\tmodel / mapping:", m , mapping
            n_to_simulate = len( mapping )
                    
            if( n_to_simulate > 0 ):
                # simulate this chunk
                this_model_parameters = []
                for i in range(n_to_simulate):
                    this_model_parameters.append( sampled_params[ mapping[i] ]) 

                # all timepoints (for all submodels) passed in one go
                t = [self.data[i].timepoints for i in range(0,len(self.data))]

                # simulate loops over submodels.  Can be via ODE, SDE or Gillespie
                sims = self.models[ m ].simulate( this_model_parameters, t, n_to_simulate, self.beta )
                if self.debug == 2:print '\t\t\tsimulation dimensions:', numpy.shape(sims) 

                # stores the trajectories and distances in a list of length nsubmodel*beta
                for i in range(n_to_simulate):
                    this_dist = [[[] for k in range(self.beta)] for n in range(self.nsubmodel)]   
                    this_traj = [[[] for k in range(self.beta)] for n in range(self.nsubmodel)]    
                    distanceaverage = [[] for k in range(self.beta)]
                    for k in range(self.beta):
                        for n in range(self.nsubmodel):     
                            samplePoints = sims[n][i,k,:,:] 
                            points = howToFitData( self.models[ m ].fit[n], samplePoints ) 
                            if do_comp == True:
                                if self.debug == 2:print '\t\t\t which distance:', self.distancefn[m][n]
                                distance = self.distancefn[m][n](points, self.data[n].values, this_model_parameters[i], m)
                                
                            else:
                                distance = 0
                                dist = True
  
                            this_dist[n][k]= distance  # group of submodels
                            this_traj[n][k]= points
                            #if self.debug == 2:print '\t\t\t distance:', this_dist

                        # distance for each model is weighted average of distances for each submodel
                        if do_comp == True:
                            for j in range(self.nsubmodel):
                                contribution[j] = self.submodelweights[j] * this_dist[j][k][0]  # submodels equally weighted at present - MAR2012
                            distanceaverage[k] = [numpy.sum(contribution)]  
                            dist = evaluateDistance(distanceaverage[k], this_epsilon ) # evaluateDistance called model by model
                       

                        if dist == True:
                            ret[mapping[i]] = ret[mapping[i]] + 1
                              
                        if self.debug == 2:print '\t\t\tdistance/this_epsilon/mapping/b:', distanceaverage[k], this_epsilon, mapping[i], ret[mapping[i]] 

                    traj[ mapping[i] ] = copy.deepcopy( this_traj ) 
                    distances[ mapping[i] ] = copy.deepcopy( this_dist ) 
                    distanceavgs[ mapping[i] ] = copy.deepcopy( distanceaverage )

        return ret[:], distances, distanceavgs, traj    
  

### CLASS FUNCTION: sampleTheModelFromPrior
### Sample a model for first iteration (multinomial sampling)
##############################################################################################################################
    def sampleTheModelFromPrior(self):
        ret = [0 for it in range(self.nbatch)]
        if self.nmodel > 1:
            for i in range(self.nbatch):
                ret[i] = statistics.w_choice( range(self.nmodel), self.modelprior )
        
        return ret[:]

          
### CLASS FUNCTION: sampleTheParameterFromPrior
### Sample particles for each model from the parameter priors.     
##############################################################################################################################        
    def sampleTheParameterFromPrior(self, sampled_models):
        ret = []
        kl=0
        for i in range(self.nbatch):         
            reti = [ 0 for it in range(self.models[ sampled_models[i] ].nparameters) ]  # nbatch particles samplet at once

            for n in range(self.models[ sampled_models[i] ].nparameters):

                if self.models[ sampled_models[i] ].prior[n][0] == 0:       # constant prior
                    reti[n]=self.models[ sampled_models[i] ].prior[n][1]
        
                if self.models[ sampled_models[i] ].prior[n][0] == 1:       # normal prior
                    reti[n]=rnd.normal( loc=self.models[ sampled_models[i] ].prior[n][1],
                                        scale=numpy.sqrt(self.models[ sampled_models[i] ].prior[n][2]) )
                
                if self.models[ sampled_models[i] ].prior[n][0] == 2:       # uniform prior
                    reti[n]=rnd.uniform( low=self.models[ sampled_models[i] ].prior[n][1],
                                         high=self.models[ sampled_models[i] ].prior[n][2])

                if self.models[ sampled_models[i] ].prior[n][0] == 3:       # lognormal prior
                    reti[n]=rnd.lognormal(mean=self.models[ sampled_models[i] ].prior[n][1],
                                          sigma=numpy.sqrt(self.models[ sampled_models[i] ].prior[n][2]) )

                if self.models[ sampled_models[i] ].prior[n][0] == 4:       # truncated normal prior
                    reti[n] = cstats.sampleGaussTrunc( mean_p=self.models[ sampled_models[i] ].prior[n][1], sd_p=numpy.sqrt(self.models[ sampled_models[i] ].prior[n][2]), lower_p = self.models[ sampled_models[i] ].prior[n][3][0], upper_p=self.models[ sampled_models[i] ].prior[n][3][1]  )
               
            ret.append( reti[:] )

        return [x[:] for x in ret]

            
### CLASS FUNCTION: sampleTheModel
### Sample (nbatch) model(s) for a new population     
##############################################################################################################################
    def sampleTheModel(self):
        ret = [0 for it in range(self.nbatch)]
    
        if self.nmodel > 1:
            for i in range(self.nbatch):
                ret[i] = statistics.w_choice( range(self.nmodel), self.margins_prev )  # multinomial sampling

            if( len(self.dead_models) > self.nmodel-1 ):
                for i in range(self.nbatch):
                    u = rnd.uniform(low=0, high=1)
                    if u > self.modelKernel:
                        # sample randomly from other (non dead) models
                        not_available = self.dead_models[:]
                        not_available.append(ret[i])
                        ss = set( not_available )
                        s = set( range(0,self.nmodel) )
                        ar = numpy.array( list(s-ss) )
                        rnd.shuffle( ar )
                        perturbed_model = ar[0]
                        ret[i] = perturbed_model
        return ret[:]
        

### CLASS FUNCTION: sampleTheParameterNonCustom
### Select nbatch particles from the previous population, perturb parameters and return  
### No need to check that prior probabilities of returned particles >0, as sampling from truncated distributions available
##############################################################################################################################         
    
    def sampleTheParameterNonCustom(self, sampled_models):
        if self.debug == 2:print "\t\t\t***sampleTheParameterNonCustom"
        ret = []
        
        for i in range(self.nbatch):
            np = self.models[ sampled_models[i] ].nparameters
            reti = [ 0 for it in range(np) ]
                        
            # sample putative particle from previous population
            p = sample_particle(self.nparticles, sampled_models[i], self.margins_prev, self.model_prev, self.weights_prev )
                
            for nn in range(np):
                reti[nn] = self.parameters_prev[ p ][nn]

            self.perturbfn[sampled_models[i]](reti, self.kernel_functions[sampled_models[i]], self.kernel_info[sampled_models[i]])

            if self.debug == 2:print "\t\t\tnew:", reti
            if self.debug == 2:print "\t\t\told:", self.parameters_prev[p]

            ret.append( reti )

        return [x[:] for x in ret]
        

### CLASS FUNCTION: sampleTheParameterCustom
### Select nbatch particles from the previous population, perturb parameters and return  
### Includes a prior probability check, as cannot be sure that custom perturbation uses truncated distributions
##############################################################################################################################  
    def sampleTheParameterCustom(self, sampled_models):
        if self.debug == 2:print "\t\t\t***sampleTheParameterCustom"
        ret = []
        
        for i in range(self.nbatch):
            np = self.models[ sampled_models[i] ].nparameters
            reti = [ 0 for it in range(np) ]

            prior_prob = -1
            while prior_prob <= 0:      #ie. prior probability must be >0
            
            # sample putative particle from previous population
                p = sample_particle(self.nparticles, sampled_models[i], self.margins_prev, self.model_prev, self.weights_prev )
                
                for nn in range(np):
                    reti[nn] = self.parameters_prev[ p ][nn] 

                prior_prob = self.perturbfn[sampled_models[i]](reti, self.kernel_info[sampled_models[i]], self.models[sampled_models[i]].prior)       # returns prior probability

                if self.debug == 2:print "\t\t\tnew:", reti
                if self.debug == 2:print "\t\t\told:", self.parameters_prev[p]

            ret.append( reti )

        return [x[:] for x in ret]
        

### CLASS FUNCTION: computeParticleWeights
### Calculate particle weights   
##############################################################################################################################         
    
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
                    x=statistics.getPdfGauss(parameter=this_param[n],
                                         mean=self.models[ this_model ].prior[n][1],
                                         scale=numpy.sqrt(self.models[ this_model ].prior[n][2]))
                
                if self.models[ this_model ].prior[n][0]==2: 
                    x=statistics.getPdfUniform(parameter=this_param[n],
                                           scale1=self.models[ this_model ].prior[n][1],
                                           scale2=self.models[ this_model ].prior[n][2])
                                               
    
                if self.models[ this_model ].prior[n][0]==3: 
                    x=statistics.getPdfLognormal(mean=self.models[ this_model ].prior[n][1],
                                                 sigma=numpy.sqrt(self.models[ this_model ].prior[n][2]),
                                                 parameter=this_param[n])

                if self.models[ this_model ].prior[n][0]==4:
                    x = cstats.pdfGaussTrunc(mean_p=self.models[this_model].prior[n][1], sd_p=numpy.sqrt(self.models[this_model].prior[n][2]), lower_p = self.models[this_model].prior[n][3][0], upper_p=self.models[this_model].prior[n][3][1],x_p=this_param[n]  )
                    
                pprob = pprob*x

            numer = self.b[k] * mprob * pprob
            denom_m = 0
            for i in range(self.nmodel):
                denom_m = denom_m + self.margins_prev[i]*getPdfModelKernel(this_model, i, self.modelKernel, self.nmodel, self.dead_models)
            denom = 0
            for j in range(self.nparticles):
                if(int(this_model) == int(self.model_prev[j]) ):
                    if self.debug == 2:
                        print "\tj, weights_prev, kernelpdf",j,self.weights_prev[j], self.getPdfParameter[this_model](this_param, self.parameters_prev[j], self.kernel_functions[this_model], self.kernel_info[this_model], self.positive_infinity, self.negative_infinity)
                    denom = denom + self.weights_prev[j] * self.getPdfParameter[this_model]( this_param, self.parameters_prev[j], self.kernel_functions[this_model], self.kernel_info[this_model], self.positive_infinity, self.negative_infinity )
                if self.debug == 2: print "\tnumer/denom_m/denom/m(t-1) : ", numer,denom_m, denom, self.margins_prev[this_model]
            self.weights_curr[k] = numer/(denom_m*denom/self.margins_prev[this_model])

            
### CLASS FUNCTION: normaliseWeights
### Divide individual particle weights by the sum of all weights; ie. normalise     
##############################################################################################################################   
    def normaliseWeights(self):
        n = sum( self.weights_curr )
        for i in range(self.nparticles):
            self.weights_curr[i] = self.weights_curr[i]/float(n)


### CLASS FUNCTION: modelMarginals
### Caculate marginal probabilities for each model; ie. sum the weights for each particle of a model   
############################################################################################################################## 
    def modelMarginals(self):
        for i in range(self.nmodel):
            self.margins_curr[i] = 0
            for j in range(self.nparticles):
                if int(self.model_curr[j]) == int(i):
                    self.margins_curr[i] = self.margins_curr[i] + self.weights_curr[j]
            
##############################################################################################################################
                                          ### END OF ABCSMC CLASS   ###
##############################################################################################################################    






##############################################################################################################################
                                          ### FUNCTION: SAMPLE PARTICLE   ###
##############################################################################################################################
def sample_particle(nparticle, selected_model, margins_prev, model_prev, weights_prev ):
    u = rnd.uniform(low=0, high=margins_prev[selected_model])
    F = 0

    for i in range(0,nparticle):
        if int(model_prev[i]) == int(selected_model) :
            F = F + weights_prev[i]

            if(F > u):
                break

    return i


##############################################################################################################################
                                          ### FUNCTION: HOW TO FIT DATA   ###
########### how to fit variable to the given data (if fit is None, data for all variables are available in order of the model)
##############################################################################################################################
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


##############################################################################################################################
                                          ### FUNCTION: GET PDF MODEL KERNEL   ###
##############################################################################################################################
def getPdfModelKernel(m, m0, modelK, nmodel, dead_models):
    ndead = len(dead_models)
   
    if(ndead == nmodel-1):
        return 1.0
    else:

        if(m == m0):
            return modelK
        else:
            return (1-modelK)/(nmodel-ndead)


##############################################################################################################################
                                          ### FUNCTION: EVALUATE DISTANCE   ###
##############################################################################################################################
def evaluateDistance(distance,epsilon):
    
    accepted = False
    for i in range(len(epsilon)):
        #print "d:", distance[i], epsilon[i]
        if(distance[i]<epsilon[i] and distance[i]>=0 ):
            accepted = True
        else: 
            accepted = False
            break
        
    return accepted


    
