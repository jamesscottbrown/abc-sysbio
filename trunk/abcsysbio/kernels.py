# module defining all things kernel

# kernel_type = 1 : component-wise uniform kernels
# kernel_type = 2 : component-wise normal kernels
# kernel_type = 3 : multi-variate normal kernels
# kernel_type = 4 : multi-variate normal kernels using the nearest neighbours of the particles
# kernel_type = 4 : multi-variate normal kernels using an "optimal covaraince matrix" 

import numpy
from numpy import random as rnd
from scipy.stats import norm
from abcsysbio import statistics


# kernel is a list of length 3 such that :
# kernel[0] contains the index of the non-constant paramameters
# kernel[1] contains the informations required to build the kernels in function getKernels - those informations are defined in the input file
# kernel[2] contains the kernel (list, matrix or dictionnary) once it has been built


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
	for param in kernel[0]:
            minimum=min(population[:,param])
            maximum=max(population[:,param])
            scale=(maximum-minimum)
            tmp.append([-scale/2.0,scale/2.0])
        kernel[2]=tmp    
    elif kernel_type == 2:
        # component-wise normal kernels
        tmp=list()
	for param in kernel[0]:
            s2w = statistics.wtvar(population[:,param], weights, method = "R")
            tmp.append(2*s2w)
            kernel[2]=tmp

    elif kernel_type == 3:
        # multi-variate normal kernel whose covariance is based on all the previous population
        pop=list()
        for param in kernel[0]:
            pop.append(population[:,param])
        cov = statistics.compute_cov(pop,weights)
        kernel[2]=2*cov
    if kernel_type == 4:
        # multi-variate normal kernel whose covariance is based on the K nearest neighbours of the particle
        k=int(kernel[1])
        D={}
        pop=list()
        # to compute the neighbours, restrain the population to the non constant parameters
        for param in kernel[0]:
            pop.append(population[:,param])
        for n in range(len(pop[0])):
            # compute the index of the neighbours
            kset = statistics.kNearestNeighEuc(n,pop,k)
            # save the coordinate of the particule (with all componants)
            pop_cur=list()            
            for param in range(npar):
                pop_cur.append(population[n,param])
            # construct the list of the neighbours given kset (restrained to the non constant components) and the corresponding weights
            subpop=list()
            subwei=list()
            for param in range(0, len(pop)):
                subpop.append([])
                for j in range(len(kset)):
                    subpop[param].append(pop[param][kset[j]])
            for j in range(len(kset)):
                subwei.append(weights[kset[j]])
            # compute the covariance and write it into the dictionnary
            D[str(pop_cur)]=2*statistics.compute_cov(subpop,subwei)
        kernel[2]=D
    if kernel_type==5:
        # multi-variate normal kernel whose covariance is the OCM
        pop=list()
        for param in kernel[0]:
            pop.append(population[:,param])
        D={}
        for n in range(pop_size):
            pop_cur=list()
            for param in range(npar):
                pop_cur.append(population[n, param])
            D[str(pop_cur)]=statistics.compute_optcovmat(pop, weights,pop_cur)
        kernel[2]=D


    return kernel




# Here params refers to one particle[
# The function changes params in place and returns the probability (which may be zero)
def perturbParticle(params, priors, kernel, kernel_type):
    np = len( priors )
    prior_prob = 1
    if kernel_type==1:
	ind=0
        for n in kernel[0]:
            params[n] = params[n] + rnd.uniform(low=kernel[2][ind][0],high=kernel[2][ind][1])
	    ind+=1

    if kernel_type==2:
	ind=0
        for n in kernel[0]:
            params[n] = rnd.normal(params[n],numpy.sqrt(kernel[2][ind]))
	    ind+=1

    if kernel_type==3:
	mean=list()
	for n in kernel[0]:
            mean.append(params[n])
        tmp = statistics.mvnd_gen(mean, kernel[2])
        ind=0
        for n in kernel[0]: 
            params[n] = tmp[ind]
            ind=ind+1
            
    if (kernel_type==4 or kernel_type==5):
	mean=list()
	for n in kernel[0]:
            mean.append(params[n])
        D=kernel[2]
        tmp = statistics.mvnd_gen(mean,D[str(params)] )
        ind=0
        for n in kernel[0]: 
            params[n] = tmp[ind]
            ind=ind+1
        
    # compute the likelihood
    prior_prob=1
    for n in range(np):
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
    if kernel_type==1:
        prob=1
        for n in range(len(kernel[0])):
            scale = (kernel[2][n][1]-kernel[2][n][0])/2.0
            scale1 = max( (params0[n] - scale), priors[n][1])
            scale2 = min( (params0[n] + scale), priors[n][2])
            kern =statistics.getPdfUniform(scale1,scale2, params[n])
            prob=prob*kern
        return prob

    elif kernel_type==2:
        prob=1
        for n in range(len(kernel[0])): 
            mean = params0[n]
            scale = numpy.sqrt(kernel[2][n])
            kern = statistics.getPdfGauss(mean,scale,params[n])
            kern = kern/auxilliary[n]
            prob=prob*kern
        return prob
    
    elif kernel_type==3:
        p0 = list()
        p = list()
        for n in kernel[0]:
              p0.append(params0[n])
              p.append(params[n])
	kern = statistics.getPdfMultinormal(p0,kernel[2],p)
	kern = kern/auxilliary
	return kern

    elif (kernel_type==4  or kernel_type==5):
        p0 = list()
        p = list()
        D=kernel[2]
        for n in kernel[0]:
              p0.append(params0[n])
              p.append(params[n])
       	kern = statistics.getPdfMultinormal(p0,D[str(params0)],p)
	kern = kern/auxilliary
	return kern
    

# Here models and parameters refer to the whole population
def getAuxilliaryInfo(kernel_type, models, parameters, model_objs, kernel ):
    nparticles = len(parameters)
    ret = []

    for k in range(nparticles):
        this_prior = model_objs[models[k]].prior
        this_kernel = kernel[models[k]]
        nparam = len(this_prior)
        if kernel_type == 2:
            ret.append( [1.0 for n in range(nparam)] )
            for n in range(len(this_kernel[0])):
                if not(this_prior[n][0]==0):
                    mean = parameters[k][n]
                    if len(this_kernel[2])==1:
                        ret[k][n]=0
                    else:
                        scale = numpy.sqrt(this_kernel[2][n])
                        ret[k][n] = norm.cdf(this_prior[n][2],mean,scale) - norm.cdf(this_prior[n][1],mean,scale)
        elif kernel_type==3:
            up=list()
            low=list()
            mean=list()
            for n in this_kernel[0]:
                low.append(this_prior[n][1])
                up.append(this_prior[n][2])
                mean.append(parameters[k][n])
            scale=this_kernel[2]
            ret.append(statistics.mvnormcdf(low,up,mean, scale))
        elif (kernel_type==4 or kernel_type==5):
            up=list()
            low=list()
            mean=list()
            for n in this_kernel[0]:
                low.append(this_prior[n][1])
                up.append(this_prior[n][2])
                mean.append(parameters[k][n])
            cur_part=list()
            for n in range(nparam):
                cur_part.append(parameters[k][n])
            D = this_kernel[2]
            scale=D[str(cur_part)]
            ret.append(statistics.mvnormcdf(low,up,mean, scale))                
        else:
            ret = [0 for i in  range(nparticles)]
    
    return ret






