import scipy
import numpy
import math
import sys

def euclidianDistance(data1, data2):
    """Returns the euclidian distance between two data sets.
    Data sets must have the same dimensions.

    ***** args *****
    
    data1, data2:
            numpy array objects with the same dimensions.

    """
    
    distance = -1
    #data1 and data2 are two numpy arrays
    if(numpy.shape(data1)!=numpy.shape(data2)):
        print "\ndata sets have different dimensions\n"
        sys.exit()
    else:
        z = (data1 - data2)**2
        distance = math.sqrt(scipy.sum(z))
	#distance = scipy.sum(z)/30

    if distance < 0:
        return None
    else:
        return distance
        
