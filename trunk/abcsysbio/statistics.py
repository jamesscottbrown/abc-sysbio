# statistical functions

from numpy import *
from numpy import random as rnd

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
