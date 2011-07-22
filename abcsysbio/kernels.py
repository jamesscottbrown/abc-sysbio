# module defining all things kernel

# kernel_type = 1 : component-wise uniform kernels
# kernel_type = 2 : component-wise normal kernels

import numpy
from numpy import random as rnd
from scipy.stats import norm
from abcsysbio import statistics

# Compute the kernels WITHIN a model by examining the previous population of particles
# populations, weights refers to particles and weights from previous population for one model
def getKernel(kernel_type, kernel, population, weights):
    pop_size = population.shape[0]
    npar = population.shape[1]
    
    if pop_size == 1:
        print "WARNING: getKernel : only one particle so adaptation is not possible"
        return kernel

    if kernel_type == 1:
        tmp=list()
        # component-wise uniform kernels
        for param in range(npar):
            minimum=min(population[:,param])
            maximum=max(population[:,param])
            scale=(maximum-minimum)
            tmp.append([-scale/2.0,scale/2.0])
            #kernel[param][1]=-scale/2.0
            #kernel[param][2]=scale/2.0
        kernel[2]=tmp    
    elif kernel_type == 2:
        # component-wise normal kernels
        tmp=list()
        for param in range(npar):
            # optimal kernel
            #s2 = numpy.var( population[:,param] ) # should be weighted
            s2w = statistics.wtvar(population[:,param], weights, method = "R")
            #print "vars:", s2, s2w

            #kernel[param][2]=2*s2w
            tmp.append(2*s2w)
            kernel[2]=tmp
    return kernel

# Here params refers to one particle
# The function changes params in place and returns the probability (which may be zero)
def perturbParticle(params, priors, kernel, kernel_type):
    np = len( priors )
    prior_prob = 1
    print 'kernel_type:', kernel_type
    for n in range(0,np):
        if not(priors[n][0]==0):
    
            if kernel_type==1:
                params[n] = params[n] + rnd.uniform(low=kernel[2][n][0],high=kernel[2][n][1])

            if kernel_type==2: 
                params[n] = params[n] + rnd.normal(kernel[2][n][0], numpy.sqrt(kernel[2][n][1]) )

            x = 1.0
            if priors[n][0]==1: 
                x=statistics.getPdfGauss(priors[n][1], numpy.sqrt(priors[n][2]), params[n])

            if priors[n][0]==2: 
                x=statistics.getPdfUniform(priors[n][1],priors[n][2],params[n])

            if priors[n][0]==3: 
                x=statistics.getPdfLognormal(priors[n][1],priors[n][2],params[n])

            prior_prob = prior_prob*x

    return prior_prob

# Here params and params0 refer to one particle each.
# Auxilliary is a vector size of nparameters
def getPdfParameterKernel(params, params0, priors, kernel, auxilliary, kernel_type):
    prob = 1
    print 'kernel_type:', kernel_type
    # loop over component-wise
    for n in range(0,len(priors)):
        kern = 0.0

        if not(priors[n][0]==0):     
            if kernel_type==1: 
                scale1 = params0[n] + kernel[2][n][0] # - range/2
                scale2 = params0[n] + kernel[2][n][1] # + range/2
                kern = statistics.getPdfUniform(scale1,scale2, params[n])
               
            if kernel_type==2:
                mean = params0[n]
                scale = numpy.sqrt(kernel[2][n])
                #CDF2 = norm.cdf(priors[n][2],mean,scale)
                #CDF1 = norm.cdf(priors[n][1],mean,scale)
                #CDF = (CDF2-CDF1)**(-1)
                kern = statistics.getPdfGauss(mean,scale,params[n])
                kern = kern/auxilliary[n]

        else: kern=1.0
        prob=prob*kern

    #print prob
    return prob

# Here models and parameters refer to the whole population
def getAuxilliaryInfo(kernel_type, models, parameters, model_objs, kernel ):

    nparticles = len(parameters)
    ret = []
    print 'kernel_type:', kernel_type
    if kernel_type == 2:
        # component-wise Normal kernels
        for k in range(nparticles):
            this_prior = model_objs[models[k]].prior
            this_kernel = kernel[models[k]]

            nparam = len(this_prior)
            ret.append( [1.0 for n in range(nparam)] )

            for n in range(nparam):
                if not(this_prior[n][0]==0):
                    mean = parameters[k][n]
                    scale = numpy.sqrt(this_kernel[2][n])
                    print 'scale:', scale
                    ret[k][n] = norm.cdf(this_prior[n][2],mean,scale) - norm.cdf(this_prior[n][1],mean,scale)
    else:
        ret = [0 for i in  range(nparticles)]
    
    return ret
