# statistical functions

from numpy import *
from numpy import random as rnd
from numpy import linalg as LA
import scipy
import math
import scipy.stats.mvn

################## multinomial sampling
def w_choice(item,weight):
    n = rnd.random_sample()
    for i in range(0,len(weight)):
        if n < weight[i]:
            break
        n = n - weight[i]
    return i

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

################ compute the pdf of lognormal distribution
def getPdfLognormal(mean,sigm,parameter):
    x = exp(-0.5*(log(parameter)-mean)*(log(parameter)-mean)/(sigm*sigm) )
    x = x/( parameter*sigm*sqrt(2*pi) )
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

################ compute weighted variance
################ http://adorio-research.org/wordpress/?p=259
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


################# compute the truncation factor for a multivariate normal kernel
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
    
    

####### Compute the k nearest neighbors of a point inside a set S of points using the euclidian distance
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


########### Compute the weighted covariance
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



#### compute the optimal covariance matrix
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
