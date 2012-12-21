
# module defining all things kernel

import sys
import numpy
from numpy import random as rnd
from scipy.stats import norm
from abcsysbio import statistics

cstats = statistics.c_statistics()


################### kernel_info ################### 
# kernel_info  is a list with length = no of models.   kernel_info[i] is always passed to the kernel functions in this script.
# kernel_info[i] is a list of length 4:
# kernel_info[i][0] is a list of indices of non-constant parameters.
# kernel_info[i][1] contains a list of upper and lower bounds for each parameter, which are defined by the parameter's prior
                    # and are used for perturbing particles with truncated distributions. 
# kernel_info[i][2] kernel info defined in the input file.
# kernel_info[i][3] contains the kernel once it has been built (using getXKinfo).



################### kernels ################### 
# kernel_type = 1 : component-wise uniform kernels
# kernel_type = 2 : component-wise normal kernels
# kernel_type = 3 : multi-variate normal kernels
# kernel_type = 4 : multi-variate normal kernels using the nearest neighbours of the particles
# kernel_type = 5 : multi-variate normal kernels using an "optimal covariance matrix" 




##############################################################################################################################
                                          ### GET KERNEL INFO  FUNCTIONS   ###

# Compute the kernels WITHIN a model by examining the previous population of particles
# populations, weights refer to particles and weights from previous population for one model

##############################################################################################################################

def getUniformKinfo(model_kinfo, population, weights):      # Get kernel information required for perturbation with uniform CW kernels
    pop_size = population.shape[0]
    npar = population.shape[1]
    tmp=list()

    if pop_size == 0:
        print "Error - No particles for a model!  Program closing; please recheck models"
        sys.exit()
    elif pop_size == 1:
        print "WARNING: getKernel : zero or one particle(s), so adaptation is not possible"
        for param in model_kinfo[0]:                       
            tmp.append([-1, 1])	
    else: 
        for param in model_kinfo[0]:
            #print "Population here", population
            minimum=min(population[:,param])                 #Max and min values for a parameter in the previous population
            maximum=max(population[:,param])
            scale=(maximum-minimum)/2.0
            tmp.append([-scale/2.0,scale/2.0])
    model_kinfo[3]=tmp              # N.B. Bounds are stored in model_kinfo[3] - different from normal CW kernel (see below)
    model_kinfo[2]=[0 for i in model_kinfo[0]]                                   
    return model_kinfo
	# model_kinfo[3] is a list of length the number of non-constant parameters. Each element of the list contains the inf and sup bound of the uniform kernel.



    
def getGaussKinfo(model_kinfo, population, weights):         # Get kernel information required for perturbation with Gaussian CW kernels
    pop_size = population.shape[0]
    npar = population.shape[1]
    tmp=list()

    if pop_size == 0:
        print "Error - No particles for a model!  Program closing; please recheck models"
        sys.exit()
    elif pop_size == 1:
        print "WARNING: getKernel : only one particle so adaptation is not possible"
        tmp=[1 for param in model_kinfo[0]]         
    else:
        for param in model_kinfo[0]:
            s2w = statistics.wtvar(population[:,param], weights, method = "R")
            tmp.append(2*s2w)                                # Variance of previous population
    model_kinfo[2]=tmp      # N.B. Variance is stored in model_kinfo[2] - different from uniform CW kernel (see above)
    model_kinfo[3]=[[0,0] for i in model_kinfo[0]]                                      
    return model_kinfo
	# model_kinfo[2] is a list of length the number of non-constant parameters. Each element of the list contains the variance.

########## AS: GET KERNEL FUNCTIONS - THREE ADDITIONAL FUNCTIONS DEPENDING ON TYPE OF MVN ###########    
#################### AS: KERNEL TYPE 3 GETKERNEL FUNCTION ###################

def getMVNkernel(model_kinfo, population, weights):
    print "getMVNkernel called"
    pop_size = population.shape[0] 
    npar = population.shape[1]
    # multi-variate normal kernel whose covariance is based on all the previous population
    if pop_size == 1:
        print "WARNING: getKernel : only one particle so adaptation is not possible"
        cov=numpy.eye(len(model_kinfo[0]))
        model_kinfo[3]=2*cov
    else:
        pop=list()
        for param in model_kinfo[0]:
            pop.append(population[:,param])
        cov = statistics.compute_cov(pop,weights)
        model_kinfo[3]=2*cov

##        # Check that model_kinfo[3] is positive definite
##        e_val,e_vec = numpy.linalg.eigh(model_kinfo[3]) # Compute eigenvalues and vectors
##        e_valt = numpy.sort(e_val)
##
##        # Algorithm taken straight from R implementation
##        if (e_valt[0]<0): # Test of positive definiteness
##            print "Not pos def"
##            # Begin to create a positive definite matrix
##            e_vecT = zip(*e_vec)
##
##            abs_eval = []
##            for i in range(len(e_val)):
##                abs_eval.append(abs(e_val[i]))
##
##            tol = len(e_val)*2*max(abs_eval)*2.220446e-16 # typical machine precision
##            #print "tol = " + repr(tol)
##            tol_intermediate = tol - e_val
##            for i in range(len(e_val)):
##                tol_intermediate[i] = max(0, tol_intermediate[i])
##            
##
##            diag_vect = numpy.diag(tol_intermediate)
##
##            dA = numpy.dot(numpy.dot(e_vec, diag_vect),numpy.transpose(e_vec))
##
##            model_kinfo[3] = (numpy.array(model_kinfo[3]) + dA).tolist()
##
##            e_val2,e_vec2 = numpy.linalg.eigh(model_kinfo[3]) # Compute eigenvalues and vectors
##            e_val2t = numpy.sort(e_val2)
##            if e_val2t[0]<0: print "Still not positive definite"
            
        
    # model_kinfo[3] is the covaraince matrix of the multivariate normal kernel of size len(kernel[0])*len(kernel[0])

    return model_kinfo


def getMVNkernelNN(model_kinfo, population, weights):
    # multi-variate normal kernel whose covariance is based on the K nearest neighbours of the particle
    pop_size = population.shape[0] 
    npar = population.shape[1]
    k=int(model_kinfo[2])
    D={}
    pop=list()
    if pop_size == 1:
        print "WARNING: getKernel : only one particle so adaptation is not possible"
        pop_cur=list()
        for param in range(npar):
            pop_cur.append(population[0,param])
            D[str(pop_cur)]=2*numpy.eye(len(model_kinfo[0])) 
    else:
        # to compute the neighbours, restrain the population to the non constant parameters
        for param in model_kinfo[0]:
            pop.append(population[:,param])
        for n in range(pop_size):
            # compute the index of the neighbours
            kset = cstats.kNearestNeighEuc(n,pop,k)
            # save the coordinate of the particule (with all components)
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

##            # Algorithm taken straight from R implementation
##            # Check that D[str(pop_cur)] is positive definite
##            e_val,e_vec = numpy.linalg.eigh(D[str(pop_cur)]) # Compute eigenvalues and vectors
##
##            if (numpy.sort(e_val)[0]<0): # Test of positive definiteness
##                # Begin to create a positive definite matrix
##                e_vecT = zip(*e_vec)
##
##                abs_eval = []
##                for i in range(len(e_val)):
##                    abs_eval.append(abs(e_val[i]))
##
##                tol = len(e_val)*2*max(abs_eval)*2.220446e-16 # typical machine precision
##                #print "tol = " + repr(tol)
##                tol_intermediate = tol - e_val
##                for i in range(len(e_val)):
##                    tol_intermediate[i] = max(0, tol_intermediate[i])
##                
##
##                diag_vect = numpy.diag(tol_intermediate)
##
##                dA = numpy.dot(numpy.dot(e_vec, diag_vect),numpy.transpose(e_vec))
##
##                D[str(pop_cur)] = (numpy.array(D[str(pop_cur)]) + dA).tolist()

            
    model_kinfo[3]=D
    # model_kinfo[3] is a dictionary with pop_size keys. Each key is string(p) where p is a particle (with nparam dimension) of the previous population. The element of the dictionnary for a given key is a covaraince matrix of size len(kernel[0])*len(kernel[0])
    return model_kinfo

def getMVNkernelOpt(model_kinfo, population, weights):
    # multi-variate normal kernel whose covariance is the OCM
    pop_size = population.shape[0] 
    npar = population.shape[1]

    pop=list()
    D={}
    if pop_size == 1:
        print "WARNING: getKernel : only one particle so adaptation is not possible"
        pop_cur=list()
        for param in range(npar):
            pop_cur.append(population[0,param])
            D[str(pop_cur)]=2*numpy.eye(len(model_kinfo[0])) 
    else:
        for param in model_kinfo[0]:
            pop.append(population[:,param])
        for n in range(pop_size):
            pop_cur=list()
            for param in range(npar):
                pop_cur.append(population[n, param])
            D[str(pop_cur)]=statistics.compute_optcovmat(pop, weights,pop_cur)

##            # Algorithm taken straight from R implementation
##            # Check that D[str(pop_cur)] is positive definite
##            e_val,e_vec = numpy.linalg.eigh(D[str(pop_cur)]) # Compute eigenvalues and vectors
##
##            if (numpy.sort(e_val)[0]<0): # Test of positive definiteness
##                # Begin to create a positive definite matrix
##                e_vecT = zip(*e_vec)
##
##                abs_eval = []
##                for i in range(len(e_val)):
##                    abs_eval.append(abs(e_val[i]))
##
##                tol = len(e_val)*2*max(abs_eval)*2.220446e-16 # typical machine precision
##                #print "tol = " + repr(tol)
##                tol_intermediate = tol - e_val
##                for i in range(len(e_val)):
##                    tol_intermediate[i] = max(0, tol_intermediate[i])
##
##                diag_vect = numpy.diag(tol_intermediate)
##
##                dA = numpy.dot(numpy.dot(e_vec, diag_vect),numpy.transpose(e_vec))
##
##                D[str(pop_cur)] = (numpy.array(D[str(pop_cur)]) + dA).tolist()
            
    model_kinfo[3]=D
    # model_kinfo[3] is a dictionnary with pop_size keys. Each key is string(p) where p is a particle (with nparam dimension) of the previous population. The element of the dictionnary for a given key is a covaraince matrix of size len(kernel[0])*len(kernel[0])
    return model_kinfo


##############################################################################################################################
                                          ### PERTURBATION    FUNCTIONS   ###

# Here params refers to one particle
# The function changes params in place and returns the probability (which may be zero)


##############################################################################################################################

def perturbParticleCW(param_values, model_kfunctions, model_kinfo):
    np = len(param_values)
    ind = 0
    for par in model_kinfo[0]:  
        sampled_number=model_kfunctions[1][ind]( a_p=param_values[par]+model_kinfo[3][ind][0], b_p=param_values[par]+model_kinfo[3][ind][1], mean_p=param_values[par], sd_p=numpy.sqrt(model_kinfo[2][ind]), lower_p=model_kinfo[1][ind][0], upper_p=model_kinfo[1][ind][1] )
        param_values[par] = sampled_number
        ind += 1


########## AS: PERTURBATION FUNCTIONS - FOUR ADDITIONAL FUNCTIONS DEPENDING ON WHETHER PRIORS ARE TRUNCATED AND TYPE OF MVN ###########    
#################### AS 1: NON-TRUNCATED PERTURB FUNCTION ###################
    
def perturbParticleMVN(params_values, model_kfunctions, model_kinfo):
#kernel, prior): Old version

    mean=list()
    for n in model_kinfo[0]: # model_kinfo[0] = list of non-constant parameters
        mean.append(params_values[n])
    kernel2Vector = numpy.reshape(model_kinfo[3], -1).tolist() # model_kinfo[3] = covariance matrix
    
    tmp = cstats.samplingMultivariateGaussReg(mean, kernel2Vector, len(mean))
    ind=0
    for n in model_kinfo[0]: 
        params_values[n] = tmp[ind] # New perturbed positions
        ind=ind+1

#################### AS 2: TRUNCATED PERTURB FUNCTION #######################

def perturbParticleMVNT(params_values, model_kfunctions, model_kinfo):

    mean=list()
    for n in model_kinfo[0]: # model_kinfo[0] = list of non-constant parameters
        mean.append(params_values[n])
    kernel2Vector = numpy.reshape(model_kinfo[3], -1).tolist() # model_kinfo[3] = covariance matrix

    lower = []
    upper = []
    for i in range(len(model_kinfo[0])):
        lower.append(model_kinfo[1][i][0])
        upper.append(model_kinfo[1][i][1])
    tmp = cstats.sampleMultivariateGaussTrunc(mean, kernel2Vector, len(mean), lower, upper)
    ind=0
    for n in model_kinfo[0]: 
        params_values[n] = tmp[ind] # New perturbed positions
        ind=ind+1

#################### AS 3: NON-TRUNCATED LOCAL PERTURB FUNCTION ###################
    
def perturbParticleMVN_local(params_values, model_kfunctions, model_kinfo):

    mean=list()
    for n in model_kinfo[0]: # model_kinfo[0] = list of non-constant parameters
        mean.append(params_values[n])
    D=model_kinfo[3]
    kernel2Vector = numpy.reshape(D[str(params_values)], -1).tolist() # model_kinfo[3] = dictionary of covariance matrices
    tmp = cstats.samplingMultivariateGaussReg(mean, kernel2Vector, len(mean))
    ind=0
    for n in model_kinfo[0]: 
        params_values[n] = tmp[ind]
        ind=ind+1

#################### AS 4: TRUNCATED LOCAL PERTURB FUNCTION ###################
    
def perturbParticleMVNT_local(params_values, model_kfunctions, model_kinfo):

    mean=list()
    for n in model_kinfo[0]: # model_kinfo[0] = list of non-constant parameters
        mean.append(params_values[n])
    D=model_kinfo[3]
    kernel2Vector = numpy.reshape(D[str(params_values)], -1).tolist() # model_kinfo[3] = dictionary of covariance matrices

    lower = zip(*model_kinfo[1])[0]
    upper = zip(*model_kinfo[1])[1]
        
    tmp = cstats.sampleMultivariateGaussTrunc(mean, kernel2Vector, len(mean), lower, upper)
    ind=0
    for n in model_kinfo[0]: 
        params_values[n] = tmp[ind]
        ind=ind+1



##############################################################################################################################
                                          ### GET PDF PARAMETER FUNCTIONS   ###

# Here params and params0 refer to one particle each.
# Auxilliary is a vector size of nparameters.  BP 13/2 - AUXILLIARY HAS BEEN REMOVED!!!!!!!!!!!!!!!!!!!!!!!!!!

##############################################################################################################################



def getPdfParameterCW (param_values, param_values0, model_kfunctions, model_kinfo, *args):
    prob = 1
    ind = 0
    for par in model_kinfo[0]:
        kern = model_kfunctions[2][ind]( x_p=param_values[par], a_p=param_values0[par]+model_kinfo[3][ind][0], b_p=param_values0[par]+model_kinfo[3][ind][1], mean_p=param_values0[par], sd_p=numpy.sqrt(model_kinfo[2][ind]), lower_p=model_kinfo[1][ind][0], upper_p=model_kinfo[1][ind][1] )
        prob=prob*kern
        ind+=1
    return prob

########## AS: PDF FUNCTIONS - FOUR ADDITIONAL FUNCTIONS DEPENDING ON WHETHER PRIORS ARE TRUNCATED AND IF GLOBAL/LOCAL ###########    
#################### AS 1: NON-TRUNCATED MVN PDF FUNCTION ###################
def getPdfParameterKernelMVN(params_values, params_values0, model_kfunctions, model_kinfo, *args):
    p0 = list()
    p = list()
    for n in model_kinfo[0]:
        p0.append(params_values0[n])
        p.append(params_values[n])
    kernel2Vector = numpy.reshape(model_kinfo[3], -1).tolist() # model_kinfo[3] = covariance matrix
    kern = cstats.pdfMultivariateGaussReg(p0, kernel2Vector, len(p0), p)
    return kern

#################### AS 2: TRUNCATED MVN PDF FUNCTION ###################
def getPdfParameterKernelMVNT(params_values, params_values0, model_kfunctions, model_kinfo, *args):
    #print "length of params_values: " + repr(len(params_values))
    p0 = list()
    p = list()
    for n in model_kinfo[0]:
        p0.append(params_values0[n])
        p.append(params_values[n])
    kernel2Vector = numpy.reshape(model_kinfo[3], -1).tolist() # model_kinfo[3] = covariance matrix
    lower = zip(*model_kinfo[1])[0]
    upper = zip(*model_kinfo[1])[1]
    kern = cstats.pdfMultivariateGaussTrunc(p0, kernel2Vector, p, lower, upper, len(p0))
    return kern

#################### AS 3: NON-TRUNCATED LOCAL PDF FUNCTION ###################
def getPdfParameterKernelMVN_local(params_values, params_values0, model_kfunctions, model_kinfo, *args):
    p0 = list()
    p = list()
    D=model_kinfo[3]
    for n in model_kinfo[0]:
        p0.append(params_values0[n])
        p.append(params_values[n])
    kernel2Vector = numpy.reshape(D[str(params.values0)], -1).tolist() # model_kinfo[3] = covariance matrix dict
    kern = cstats.pdfMultivariateGaussReg(p0, kernel2Vector, len(p0), p)
    return kern

#################### AS 4: TRUNCATED LOCAL PDF FUNCTION ###################
def getPdfParameterKernelMVNT_local(params_values, params_values0, model_kfunctions, model_kinfo, *args):
    p0 = list()
    p = list()
    D=model_kinfo[3]
    for n in model_kinfo[0]:
        p0.append(params_values0[n])
        p.append(params_values[n])
    kernel2Vector = numpy.reshape(D[str(params.values0)], -1).tolist() # model_kinfo[3] = covariance matrix dict
    lower = zip(*model_kinfo[1])[0]
    upper = zip(*model_kinfo[1])[1]
    kern = cstats.pdfMultivariateGaussTrunc(p0, kernel2Vector, p, lower, upper, len(p0))
    return kern





