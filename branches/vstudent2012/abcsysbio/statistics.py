# statistical functions

from numpy import *
from numpy import random as rnd
from numpy import linalg as LA
import scipy
import math
import scipy.stats.mvn
from ctypes import *
import os


########### Contents of statistical functions #############
## NB:  These contain a mixture of pure python and c-type python functions
##      The c-type functions are flagged with a '(CT)' flag in the comments
##      c-type functions are defined as class functions in a 'statistics' class and should be called as such

### S: Sampling functions
        ## S.Reg: Regular (Non-truncated) sampling
		## S.Reg.0:     Multinomial (CT)
		## S.Reg.1:     Uniform (CT)
		## S.Reg.2:     Normal (CT)
		## S.Reg.3:     Multivariate Normal (CT)
		
	## S.Trunc: Truncated sampling
		## S.Trunc.1:   Uniform (CT)
		## S.Trunc.2:   Normal (CT)
		## S.Trunc.3a:  Multivariate Normal (non-optimal) (CT)
		## S.Trunc.3b:  Multivariate Normal (MCMC - optimal) (CT)

### P: PDF functions		
	## P.Reg: Regular (Non-truncatd) PDF
		## P.Reg.1:     Uniform (CT)
		## P.Reg.2:     Normal (CT)
		## P.Reg.3:     Multivariate Normal (CT)

	## P.Trunc: Truncated PDF
		## P.Trunc.1:   Uniform (CT)
		## P.Trunc.2:   Normal (CT)
		## P.Trunc.3a:  Multivariate Normal (Genz) (CT)

        ## P.Log: Lognormal PDFs
                ## P.Log.1      Regular lognormal PDF
                ## P.Log.2      Truncated lognormal PDF

### M: Miscellaneous functions
                ## M.1:         Compute weighted variance
                ## M.2:         Compute k-nearest neighbours
                ## M.3:         Compute weighted covariance
                ## M.4:         Compute optimal covariance matrix
                ## M.5:         Compute truncation factor for MVN
                ## M.6:         Integral of multivariate normal for rectangular region

### O: Old Python functions (replaced with c-type functions above)
                ##              Multinomial sampling, MVN sampling, PDF uniform,
                ##              PDF Normal, PDF MVN

#################################################################
class c_statistics():

    def __init__(self):
        dll = os.path.join(os.path.split(os.path.realpath(__file__))[0],'libstats.so.1.0') 
        #print "DLL", dll
        self.libmylib = CDLL(dll) 
        self.libmylib.setRandomGenerator()
        


####### S.Reg.0: Multinomial sampling (CT) #######
    def sample_multinomial(self,weight_p, nweight_p):
        arr_type = nweight_p*c_double
        weight = arr_type()
        #weight = arr_type(weight_p)
        for i in range(nweight_p):#try and find an alternative to the loop
            weight[i-1] = weight_p[i]#c goes 0,1,2,3 whilst python 1,2,3,4
        nweight=c_int(nweight_p)
        i = c_int(0)
        i = self.libmylib.sample_multinomial(byref(weight), nweight)
        return i
        #print "multinomial sampling:", weight[i]


####### S.Reg.1: Uniform sampling - Regular (CT) #######
    def sampleUniformReg(self,a_p, b_p, **kwargs):
        a=c_double(a_p)
        b=c_double(b_p)
        sample=c_double(0)
        self.libmylib.sampleUniform( a, b, byref(sample) )
        #print "sample from uniform", sample.value
        return sample.value


####### S.Reg.2: Normal sampling - Regular (CT) #######
    def sampleGaussReg(self, mean_p, sd_p, **kwargs):
        mean = c_double(mean_p)
        sd = c_double(sd_p)
        sample = c_double(0)

        self.libmylib.sampleGauss( mean, sd, byref(sample) )
        #print "sample from normal", sample.value
        return sample.value


####### S.Reg.3: Multivariate normal sampling - Regular (CT) #######
    def samplingMultivariateGaussReg(self, mean_p, A_p, sizeA_p):
        sizeA=c_int(sizeA_p)
        arr_type = sizeA_p*c_double
        arr_type2 = sizeA_p*sizeA_p*c_double
        mean = arr_type()
        for i in range(sizeA_p):
            mean[i] =mean_p[i]
        sample = arr_type()
        for i in range(sizeA_p):
            sample[i] =0
        A= arr_type2()
        for i in range(sizeA_p*sizeA_p):
            A[i]=A_p[i]
        self.libmylib.samplingMultivariateGauss(byref(mean), byref(A), sizeA, byref(sample))
        result=[]
        for i in range(sizeA_p):
            result.append(sample[i])
        return result


####### S.Trunc.1: Uniform sampling - Truncated (CT) #######
    def sampleUniformTrunc(self,a_p, b_p, lower_p, upper_p, **kwargs):
        a=c_double(a_p)
        b=c_double(b_p)
        lower=c_double(lower_p)
        upper=c_double(upper_p)
        sample=c_double(0)
        self.libmylib.sampleTruncatedUniform( a, b, byref(sample), lower, upper )
        #print "sample from truncated uniform", sample.value
        return sample.value


####### S.Trunc.2: Normal sampling - Truncated (CT) #######
    def sampleGaussTrunc(self, mean_p, sd_p, lower_p, upper_p, **kwargs):
        mean = c_double(mean_p)
        sd = c_double(sd_p)
        sample = c_double(0)
        lower=c_double(lower_p)
        upper=c_double(upper_p)
        self.libmylib.sampleTruncatedGauss( mean, sd, byref(sample), lower, upper )
        #print "sample from normal", sample.value
        return sample.value


####### S.Trunc.3a: Multivariate normal sampling - Truncated, non-optimal (CT) #######
    def sampleMultivariateGaussTrunc(self, mean_p, A_p, sizeA_p,lower_p, upper_p):	
        sizeA=c_int(sizeA_p)
        arr_type = sizeA_p*c_double
        arr_type2 = sizeA_p*sizeA_p*c_double
        mean = arr_type()
        lower=arr_type()
        upper=arr_type()
        for i in range(sizeA_p):
            mean[i] =mean_p[i]
            lower[i]=lower_p[i]
            upper[i]=upper_p[i]
        sample = arr_type()
        for i in range(sizeA_p):
            sample[i] =0

        A= arr_type2()
        for i in range(sizeA_p*sizeA_p):
            A[i]=A_p[i]	
        self.libmylib.sampleMultivariateTGauss(byref(mean), byref(A), sizeA, byref(sample), byref(lower), byref(upper))
        result=[]
        for i in range(sizeA_p):
            result.append(sample[i])
        return result


####### S.Trunc.3b: Multivariate normal sampling - Truncated, MCMC (CT) #######
    def MCMCsampleMultivariateGaussTrunc(self, A_p, mean_p, lower_p, upper_p, sizeA_p):
        sizeA=c_int(sizeA_p)
        arr_type = sizeA_p*c_double
        arr_type2 = sizeA_p*sizeA_p*c_double
        mean = arr_type()
        lower=arr_type()
        upper=arr_type()
        for i in range(sizeA_p):
            mean[i] =mean_p[i]
            lower[i]=lower_p[i]
            upper[i]=upper_p[i]
        sample = arr_type()
        result=arr_type()
        for i in range(sizeA_p):
            sample[i] =0
            result[i]=0
        A= arr_type2()
        for i in range(sizeA_p*sizeA_p):
            A[i]=A_p[i]	
        upperd=c_double(0)
        lowerd=c_double(0)
                    
        self.libmylib.MCMCsampleMultivariateTGauss(byref(A), byref(mean), byref(lower), byref(upper), sizeA, byref(result), byref(sample), byref(upperd), byref(lowerd))
        final_result=[]
        for i in range(sizeA_p):
            final_result.append(result[i])
        return final_result	


####### P.Reg.1: Uniform PDF - Regular (CT) #######
    def pdfUniformReg(self, x_p, a_p, b_p, **kwargs):
        x=c_double(x_p)
        a=c_double(a_p)
        b=c_double(b_p)
        density=c_double(0)
        self.libmylib.pdfUniform( x, a, b, byref(density) )
        #print "density from uniform", density.value
        return density.value


####### P.Reg.2: Normal PDF - Regular (CT) #######
    def pdfGaussReg(self, x_p,mean_p,sd_p, **kwargs):
        mean=c_double(mean_p)
        sd=c_double(sd_p)
        x=c_double(x_p)
        density=c_double(0)

        self.libmylib.pdfGauss( mean, sd, x, byref(density))
        #print "density from normal", density.value
        return(density.value)


####### P.Reg.3: Multivariate normal PDF - Regular (CT) #######
    def pdfMultivariateGaussReg(self, mean_p, A_p,sizeA_p, x_p):
        sizeA=c_int(sizeA_p)
        arr_type = sizeA_p*c_double
        mean = arr_type()
        x=arr_type()
        for i in range(sizeA_p):
            mean[i] =mean_p[i]
            x[i]=x_p[i]
        arr_type2 = sizeA_p*sizeA_p*c_double
        A= arr_type2()
        for i in range(sizeA_p*sizeA_p):
            A[i]=A_p[i]
        density = c_double(0)

        self.libmylib.pdfMultivariateGauss(byref(mean), byref(A), sizeA, byref(x), byref(density))
        #print "density from MVN", density.value
        return density.value


####### P.Trunc.1: Uniform PDF - Truncated (CT) #######
    def pdfUniformTrunc(self,x_p, a_p, b_p, lower_p, upper_p, **kwargs):
        x=c_double(x_p)
        a=c_double(a_p)
        b=c_double(b_p)
        lower=c_double(lower_p)
        upper=c_double(upper_p)
        density=c_double(0)

        self.libmylib.pdfTruncatedUniform(x, a, b, lower, upper, byref(density))
        #print "density from truncated uniform", density.value
        return density.value

####### P.Trunc.2: Normal PDF - Truncated (CT) #######
    def pdfGaussTrunc(self,mean_p, sd_p, x_p, lower_p, upper_p, **kwargs):
        mean=c_double(mean_p)
        sd=c_double(sd_p)
        x=c_double(x_p)
        density=c_double(0)
        lower=c_double(lower_p)
        upper=c_double(upper_p)

        self.libmylib.pdfTruncatedGauss(mean, sd, x, byref(density), lower, upper)
        #print "density from truncated normal", density.value
        return density.value


####### P.Trunc.3: Multivariate normal PDF - Truncated (CT) #######
    def pdfMultivariateGaussTrunc(self, mean_p, A_p, x_p, lower_p, upper_p, sizeA_p):
        sizeA=c_int(sizeA_p)
        arr_type = sizeA_p*c_double
        arr_type2 = sizeA_p*sizeA_p*c_double
        mean = arr_type()
        lower=arr_type()
        upper=arr_type()
        x=arr_type()
        for i in range(sizeA_p):
            mean[i] =mean_p[i]
            lower[i]=lower_p[i]
            upper[i]=upper_p[i]
            x[i]=x_p[i]
        A= arr_type2()
        for i in range(sizeA_p*sizeA_p):
            A[i]=A_p[i]
        density=c_double(0)


        ######### Accuracy tuning of truncated MVN pdf - Start #########

        ## Three parameters to tune - errorlimit, MCconfidence, and NMax 
        ## Eg. errorlimit=0.05, MCconfidence=1.28155, NMax=2 implies the following: 
        ## There is a ~90% probability that a calculated pdf value has an error <= 5%; 
        ## however a maximum of 2 iterations (within the C code) is allowed

        errorlimit=c_double(0.05) # Expected error of the truncated probability distribution function (- used to normalise the pdf)
        MCconfidence=c_double(1.28155) # Std MC confidence factor
        NMax=c_int(2) 

        ######### Accuracy tuning of truncated MVN pdf - End   #########
        
                    
        self.libmylib.pdfTMultivariateGauss(byref(mean), byref(A), sizeA, byref(x), byref(lower), byref(upper), byref(density),
                                            errorlimit, MCconfidence, NMax)
        #print "density from truncated MVN", density
        return density.value



####### P.Log.1: Lognormal PDF - Regular #######
def getPdfLognormal(mean,sigma,parameter):
    x = exp(-0.5*(log(parameter)-mean)*(log(parameter)-mean)/(sigma*sigma) )
    x = x/( parameter*sigma*sqrt(2*pi) )
    #print "density from lognormal", x
    return x

####### P.Log.2: Lognormal PDF - Truncated (CT) #######
## Not coded - not needed as no Lognormal kernel used



####### M: Miscellaneous functions #######
####### M.1: compute weighted variance - http://adorio-research.org/wordpress/?p=259 #######
def wtvar(X, W, method = "R"):
    sumW = sum(W)
    if method == "nist":
        xbarwt = sum([w * x for w,x in zip(W, X)])/sumW    # fixed.2009.03.07, divisor added.
        Np = sum([ 1 if (w != 0) else 0 for w in W])
        D = sumW * (Np-1.0)/Np
        return sum([w * (x - xbarwt)**2 for w,x in zip(W,X)])/D
    else: # default is R 
        sumW2 = sum([w **2 for w in W])
        xbarwt = sum([(w * x)  for (w,x) in zip(W, X)])/sumW
        return sum([(w * (x - xbarwt)**2) for (w,x) in zip(W, X)])* sumW/(sumW**2 - sumW2)


####### M.2: Compute the k nearest neighbors of a point inside a set S of points using the euclidian distance #######
def kNearestNeighEuc(ind,S,k):
    n = len(S[0])
    SS = array(S)
    aa=zeros((len(S),n))
    for param in range(len(S)):
        aa[param,:]=S[param][ind]*ones((1,n))
    kmin = list()
    dist=sum((aa-SS)**2,2)
    M=max(dist)
    for i in range(min(k,n)):
        im = argmin(dist)
        kmin.append(im)
        dist[im]=M
    return kmin


####### M.3: Compute the weighted covariance #######
def compute_cov(population,w):
    dimvar = len(population)
    nbsamp = len(population[0])
    cov = zeros( [dimvar,dimvar], float)
    m=list()
    for k in range(dimvar):
        m.append(0)
        for i in range(nbsamp):
            m[k] += w[i]*population[k][i]
        m[k]=m[k]/sum(w)
    for i in range(nbsamp):
        for j in range(dimvar):
            for k in range(j):
                cov[j,k] += w[i]*(population[j][i] - m[j])*(population[k][i] - m[k])
    cov =cov + transpose(cov)
    for i in range(nbsamp):
        for j in range(dimvar):
            cov[j,j] += w[i]*(population[j][i] - m[j])*(population[j][i] - m[j])    
    for j in range(dimvar):
        for k in range(dimvar):
            cov[j,k]=cov[j,k]/sum(w)
    return cov


####### M.4: compute the optimal covariance matrix #######
def compute_optcovmat(population,w,p):
    dimvar = len(population)
    nbsamp = len(population[0]) 
    C = zeros( [dimvar,dimvar], float)
    for n in range(nbsamp):
        for i in range(dimvar):
            for j in range(i):
                C[i,j]+=w[n]*(population[i][n]-p[i])*(population[j][n]-p[j])
    C= C + transpose(C)
    for i in range(dimvar):
        C[i,i]+=w[n]*(population[i][n]-p[i])*(population[i][n]-p[i])    
    for j in range(dimvar):
        for k in range(dimvar):
            C[j,k]=C[j,k]/sum(w)    
    return C


####### M.5: Compute the truncation factor for a multivariate normal kernel #######
def mvstdnormcdf(lower,upper,corrcoef, **kwds):
    n = len(lower)
    #don't know if converting to array is necessary,
    #but it makes ndim check possible
    lower = array(lower)
    upper = array(upper)
    corrcoef = array(corrcoef)
    
    correl = zeros(n*(n-1)/2.0)  #dtype necessary?
    
    if (lower.ndim != 1) or (upper.ndim != 1):
        raise ValueError, 'can handle only 1D bounds'
    if len(upper) != n:
        raise ValueError, 'bounds have different lengths'
    if n==2 and corrcoef.size==1:
        correl = corrcoef
        #print 'case scalar rho', n
    elif corrcoef.ndim == 1 and len(corrcoef) == n*(n-1)/2.0:
        #print 'case flat corr', corrcoeff.shape
        correl = corrcoef
    elif corrcoef.shape == (n,n):
        #print 'case square corr',  correl.shape
        for ii in range(n):
            for jj in range(ii):
                correl[ jj + ((ii-2)*(ii-1))/2] = corrcoef[ii,jj]
    else:
        raise ValueError, 'corrcoef has incorrect dimension'

    if not 'maxpts' in kwds:
        if n >2:
            kwds['maxpts'] = 10000*n

    lowinf = isneginf(lower)
    uppinf = isposinf(upper)
    infin = 2.0*ones(n)
    
    putmask(infin,lowinf,0)# infin.putmask(0,lowinf)
    putmask(infin,uppinf,1) #infin.putmask(1,uppinf)
    #this has to be last
    putmask(infin,lowinf*uppinf,-1)
    error, cdfvalue, inform = scipy.stats.mvn.mvndst(lower,upper,infin,correl,**kwds)
    if inform:
        print 'something wrong', inform, error, cdfvalue
       # print 'lower=', lower
       # print 'upper=', upper
       # print 'corrcoef=', corrcoef
    return cdfvalue

####### M.6: Integral of multivariate normal for rectangular region #######
def mvnormcdf(lower, upper, mu, cov, **kwds):
    '''find integral of multivariate normal for rectangular region'''
    
    lower = array(lower)
    upper = array(upper)
    cov = array(cov)
    stdev = sqrt(diag(cov))
    lower = (lower - mu)/stdev
    upper = (upper - mu)/stdev
    divrow = atleast_2d(stdev)
    corr = cov/divrow/divrow.T
    #v/sqrt(atleast_2d(diag(covv)))/sqrt(atleast_2d(diag(covv))).T

    return mvstdnormcdf(lower, upper, corr, **kwds)


####### O: Old Python functions #######
################## multinomial sampling
def w_choice(item,weight):
    n = rnd.random_sample()
    for i in range(0,len(weight)): 
        if n < weight[i]:
            break
        n = n - weight[i] 
    return i ## AS: The probabilities are equal to the weights

######### sample according to a multivariate normal
def mvnd_gen(mean, cov):
    a=list(rnd.normal(0,1,len(mean)))
    lambdas, vect = LA.eig(cov)
    print lambdas
    tmp=mat(vect)*mat(diag(sqrt(lambdas)))*transpose(mat(a))
    res=list()
    for i in range(len(mean)):
       res.append(mean[i]+tmp[i,0]) 
    return res

################## compute the pdf of uniform distribution
def getPdfUniform(scale1,scale2,parameter):
    if ((parameter>scale2) or (parameter<scale1)):
        return 0.0
    else:
        return 1/(scale2-scale1)

################## compute the pdf of gauss distribution
def getPdfGauss(mean,scale,parameter):
    x = exp( -0.5*(parameter-mean)*(parameter-mean)/(scale*scale) )
    x = x/( scale*sqrt(2*pi) )
    return x

############ compute the pdf of a multinormal distribution
def getPdfMultinormal(x,cov,y):
    a=0
    k = len(x) 
    inv=LA.inv(cov)
    for i in range(k):
      for j in range(k):
          a+=inv[i,j]*(x[i]-y[i])*(x[j]-y[j])
    det = LA.det(cov)
    return exp(-1.0*a/2.0)/(sqrt((2*pi)**k * det))







