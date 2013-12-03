from numpy import *
from numpy.linalg import *
# import numdifftools as nd

twidth = 0.5
total_time = 100
total_points = 100/twidth
npulse = 10

# define moving average using convolution
def moving_average(x, k):
    window = ones(int(k))/float(k)
    return convolve(x, window, 'same')

# detect maxima and minima
def detect_pulses( x ):
    xa = moving_average(x, 9)
    xd = diff(xa)
   
    statpts = (where(xd[:-1] * xd[1:] < 0))[0]

    maxima = []
    for i in range(len(statpts)):
  
        if statpts[i] > 0:
            if (xd[ statpts[i]+1 ] -  xd[ statpts[i] ] < 0):
                maxima.append( statpts[i] )

    minima = []
    for i in range(len(statpts)):
  
        if statpts[i] > 0:
            if (xd[ statpts[i]+1 ] -  xd[ statpts[i] ] > 0):
                minima.append( statpts[i] )

    # convert maxima and minima back to normal coordinates
    xmaxima = add(maxima,1)
    xminima = add(minima,1)

    return [xmaxima, xminima]


def distance(data1, data2, parameters, model):
    # data1 is simulated, and has shape npoints x beta
    # model is the model number

    sig = data1[:,0]

    xmaxima, xminima = detect_pulses(sig)

    if len(xmaxima) == 0 or len(xminima) == 0:
        return [ None ]

    av_sig = mean( sig )
    av_min = mean( sig[xminima] )
    min_max = min( sig[xmaxima] )
    av_dmax = mean(diff(xmaxima))
    av_dmin = mean(diff(xminima))

    d1 = (len(xmaxima) - npulse)*(len(xmaxima) - npulse)
    d2 = (len(xminima) - npulse)*(len(xminima) - npulse)
    d3 = av_min
    d4 = 1/min_max
    d5 = 1/av_sig 
    d6 = (av_dmax - (total_points/npulse))*(av_dmax - (total_points/npulse))
    d7 = (av_dmin - (total_points/npulse))*(av_dmin - (total_points/npulse))
    
    ## print d1, d2, d3, d4, d5, d6, d7
    return [d1, d2, d3, d4, d5, d6, d7]

    
