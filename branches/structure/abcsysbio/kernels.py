# module defining all things kernel

# kernel_type = 1 : component-wise uniform kernels
# kernel_type = 2 : component-wise normal kernels
# kernel_type = 3 : multi-variate normal kernels
# kernel_type = 4 : multi-variate normal kernels using the nearest neighbours of the particles
# kernel_type = 4 : multi-variate normal kernels using an "optimal covariance matrix" 

import numpy
from numpy import random as rnd
from scipy.stats import norm
from abcsysbio import statistics
import link_stats

# kernel is a list of length 3 such that :
# kernel[0] contains the index of the non-constant paramameters
# kernel[1] contains the informations required to build the kernels in function getKernels - those informations are defined in the input file
# kernel[2] contains the kernel (list, matrix or dictionnary) once it has been built

# Compute the kernels WITHIN a model by examining the previous population of particles
# populations, weights refers to particles and weights from previous population for one model
def getKernel(kernel_type, kernel, population, weights, priors, link_info):
    pop_size = population.shape[0]
    npar = population.shape[1]

    # kernel type is always 1 - component-wise uniform kernels
    tmp=list()
    for n in kernel[0]:
        
        if priors[n][0] == 4:
        # if False
            # if switch parameter
            #xx = numpy.array( population[:,n] )
            #pm1 = len( numpy.where(xx == -1)[0] )/float(len(xx))
            #p0  = len( numpy.where(xx ==  0)[0] )/float(len(xx))
            #pp1 = len( numpy.where(xx == +1)[0] )/float(len(xx))
            #tmp.append( [pm1, p0, pp1] )
            #print "switch kernel:", pm1, p0, pp1

            tmp.append( [0, 0] )
        else:
            # if continuous parameter
            minimum=min(population[:,n])
            maximum=max(population[:,n])
            scale=(maximum-minimum)
            tmp.append([-scale/2.0,scale/2.0])

    kernel[2]=tmp
    # kernel[2] is a list of length the number of non-constant parameters.
    # if continuous parameter element of the list contains the inf and sup bound of the uniform kernel.
    # if switch parameter then it contains p=-1, p=0, p=+1

    link_info.getKernels(population, weights)

    return kernel


# Here params refers to one particle
# The function changes params in place and returns the probability (which may be zero)
def perturbParticle(params, priors, kernel, kernel_type, special_cases, link_info):
    np = len( priors )
    prior_prob = 1

    #special_cases = 0

    if special_cases == 1:
        # this is the case where kernel is uniform and all priors are uniform
        # the perturbation can be made "prior aware" so that we sample efficiently
        ind=0
        for n in kernel[0]:
            if priors[n][0] == 4:
                # skip but increment ind
                ind+=1
            else:
                lflag = (params[n] + kernel[2][ind][0]) < priors[n][1]
                uflag = (params[n] + kernel[2][ind][1]) > priors[n][2]

                lower = kernel[2][ind][0]
                upper = kernel[2][ind][1]
                if lflag == True:
                    lower = -(params[n] - priors[n][1])
                if uflag == True:
                    upper = priors[n][2] - params[n]

                delta = 0
                positive = False
                if lflag == False and uflag == False:
                    # proceed as normal
                    delta = rnd.uniform(low=kernel[2][ind][0],high=kernel[2][ind][1])
                else:
                    # decide if the particle is to be perturbed positively or negatively
                    positive = rnd.uniform(0,1) > abs(lower)/(abs(lower)+upper)

                    if positive == True:
                        # theta = theta + U(0, min(prior,kernel) )
                        delta = rnd.uniform(low=0, high=upper)
                    else:
                        # theta = theta + U( max(prior,kernel), 0 )
                        delta = rnd.uniform(low=lower, high=0)

                params[n] = params[n] + delta
                ind+=1


        # perturb all the links jointly together
        params = link_info.perturbLinks(params)

        # this is not the actual value of the pdf but we only require it to be non zero
        return 1.0

    else:
        if kernel_type==1:
            ind=0
            # n refers to the index of the parameter (integer between 0 and np-1)
            # ind is an integer between 0 and len(kernel[0])-1 which enables to determine the kernel to use

            for n in kernel[0]:
                if priors[n][0] == 4:
                    # skip but increment ind
                    ind+=1
                else:
                    params[n] = params[n] + rnd.uniform(low=kernel[2][ind][0],high=kernel[2][ind][1])
                    ind+=1

            # perturb all the links jointly together
            ##params = link_stats.perturbLinks(kernel, params, priors, link_info)
            params = link_info.perturbLinks(params)

        # compute the probability under the uniform prior
        prior_prob=1
        for n in range(np):
            x = 1.0
            
            if priors[n][0]==2: 
                x=statistics.getPdfUniform(priors[n][1],priors[n][2],params[n])
            prior_prob = prior_prob*x

        ##print "p(params):", prior_prob

        ## get the prior for this set of links
        prior_links = link_info.getLinksPriorPdf(params)
        ## print "p(links):", prior_links
        prior_prob = prior_prob*prior_links

        ##print "done perturbation"
        return prior_prob


# Here params and params0 refer to one particle each.
# Auxilliary is a vector size of nparameters
def getPdfParameterKernel(params, params0, priors, kernel, auxilliary, kernel_type, link_info):
    if kernel_type==1:
        prob=1
	# n refers to the index of the parameter (integer between 0 and np-1)
	# ind is an integer between 0 and len(kernel[0])-1 which enables to determine the kernel to use
	ind=0
        for n in kernel[0]:
	    if priors[n][0] == 4:
                # skip but increment ind
                ind += 1
            else:
                kern = statistics.getPdfUniform(params0[n]+kernel[2][ind][0],params0[n]+kernel[2][ind][1], params[n])
                ##print "n, kern:,", n, kern, params[n], params0[n], kernel[2][ind][0], kernel[2][ind][1]
                prob=prob*kern
                ind += 1

        # get the probability of this set of links given the previous set of links
        probl = link_info.getLinksKernelPdf(params, params0)
        ## print "prob params / links", prob, probl
        prob=prob*probl
        
        return prob

# Here models and parameters refer to the whole population
def getAuxilliaryInfo(kernel_type, models, parameters, model_objs, kernel ):
    nparticles = len(parameters)
    ret = [0 for i in  range(nparticles)]
    return ret






