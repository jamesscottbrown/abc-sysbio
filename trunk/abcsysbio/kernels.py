# module defining all things kernel

from numpy import random as rnd
from scipy.stats import norm
from abcsysbio import statistics

################## compute the kernels by examining the previous population of particles
def getKernel(kern,population):
    #print "\t\t***getKernel"

    pop_size = population.shape[0]
    #print population
    #print "\t\tpop_size, params", pop_size, len(population[0])
    
    if pop_size == 1:
        return kern

    for param in range(len(population[0])):
        #print "\t\tdist:", population[:,param]
        #print "\t\told kernel", kern[param]
        if kern[param][0]==1: 
            # Uniform kernels
            minimum=min(population[:,param])
            maximum=max(population[:,param])
            scale=(maximum-minimum)
            kern[param][1]=-scale/2.0
            kern[param][2]=scale/2.0
        if kern[param][0]==2:
            # Gaussian kernels
            var = var( population[:,param] ) # should be weighted
            kern[param][2]=2*sqrt(var)

        #print "\t\tnew kernel", kern[param]
    return kern

def perturbParticle(params, priors, kernel):
    np = len( priors )
    prior_prob = 1

    for n in range(0,np):
         
        if kernel[n][0]==1:
            params[n] = params[n] + rnd.uniform(low=kernel[n][1],high=kernel[n][2])

        if kernel[n][0]==2: 
            params[n] = params[n] + rnd.normal(kernel[n][1],kernel[n][2])

        x = 1.0
        if priors[n][0]==2: 
            x=statistics.getPdfUniform(priors[n][1],priors[n][2],params[n])

        if priors[n][0]==3: 
            x=statistics.getPdfLognormal(priors[n][1],priors[n][2],params[n])
            
        prior_prob = prior_prob*x

    return prior_prob

def getPdfParameterKernel(params, params0, priors, kernel):
    prob = 1
    for n in range(0,len(priors)):
        kern = 0

        if not(priors[n][0]==0):
                       
            if kernel[n][0]==1: 
                scale = (kernel[n][2]-kernel[n][1])/2.0
                scale1 = max( (params0[n] - scale), priors[n][1])
                scale2 = min( (params0[n] + scale), priors[n][2])
                kern = statistics.getPdfUniform(scale1,scale2, params[n])

            if kernel[n][0]==2:
                mean = params0[n]
                scale = kernel[n][2]
                CDF2 = norm.cdf(priors[n][2],mean,scale)
                CDF1 = norm.cdf(priors[n][1],mean,scale)
                CDF = (CDF2-CDF1)**(-1)
                kern = statistics.getPdfGauss(mean,scale,params[n])
                kern = kern*CDF

        else: kern=1.0
        prob=prob*kern

    return prob

