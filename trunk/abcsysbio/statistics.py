# statistical functions

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
    x=((math.sqrt(2*math.pi)*scale)**(-1))*math.exp(-1.0*((parameter-mean)**2)/(2.0*scale**2))
    return x


################ compute the pdf of lognormal distribution
def getPdfLognormal(mean,sigm,parameter):
    x=((sqrt(2*math.pi)*sigm*parameter)**(-1))*exp(-0.5*(((math.ln(parameter)-sigm)**2)/sigm)**2)
    return x

