from numpy import *
import savitzky_golay

# define moving average using convolution
# k must be odd numbered integer
def moving_average_edge_corr(x, k):
    if k % 2 != 1 or k < 1:
        raise TypeError("window_size size must be a positive odd number")

    window = ones(int(k))/float(k)

    padded = zeros([len(x)+(k-1)])
    padded[0:k/2] = x[0:k/2]
    padded[k/2:k/2 + len(x)] = x[:]
    padded[k/2 + len(x):len(padded)] =  x[len(x)-k/2:len(x)]
    
    ret = convolve(padded, window, 'valid')
    return ret

# define moving average using convolution
def moving_average_basic(x, k):
    window = ones(int(k))/float(k)
    return convolve(x, window, 'same')

def detect_pulses_full( x, times ):
    nt = len(times)

    # SG
    xs = savitzky_golay.savitzky_golay(x, window_size=11, order=1, deriv=0)
    xd = savitzky_golay.savitzky_golay(x, window_size=11, order=1, deriv=1)

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

    # fit line sig = p0 + p1*time
    A = vstack([times, ones(nt)]).T
    p1, p0 = linalg.lstsq(A, x )[0]
    fit = p0 + p1*times
    residuals = xa-fit

    # calculate mad of residuals
    #residuals = xa-fit
    #md = median(residuals)
    #mad = median( residuals - md )/0.67

    # apply signal to noise check on the peaks, in the de-trended data
    w = 1
    thresh = 2

    peaks = []
    troughs = []
    for i in range(len(xmaxima)):
        loc = xmaxima[i]
        wl = loc-w
        wu = loc+w
        if wl < 0: wl = 0
        if wu > len(x)-1: wu = len(x)-1
        
        sig = mean(residuals[ wl:wu ])
        if sig > thresh:
            #print sig
            peaks.append(loc)

    for i in range(len(xminima)):
        loc = xminima[i]
        wl = loc-w
        wu = loc+w
        if wl < 0: wl = 0
        if wu > len(x)-1: wu = len(x)-1
        
        sig = mean(residuals[ wl:wu ])
        if -sig > thresh:
            #print sig
            troughs.append(loc)

    return [xmaxima, xminima, xs, xa, peaks, troughs, p1, fit]

def detect_pulses_basic( x, times, debug=False ):
    nt = len(times)
    
    # SG
    if debug == True:
        xs = savitzky_golay.savitzky_golay(x, window_size=11, order=1, deriv=0)
    else:
        xs = 0

    xd = savitzky_golay.savitzky_golay(x, window_size=11, order=1, deriv=1)

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

    return [xmaxima, xminima, xs, xd]
