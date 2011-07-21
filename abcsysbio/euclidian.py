import numpy
import sys

def euclidianDistance(data1, data2, parameters, model):
    """Returns the euclidian distance between two data sets.
    Data sets must have the same dimensions.

    ***** args *****
    
    data1, data2:
            numpy array objects with the same dimensions.

    """
    
    distance = -1
    #data1 and data2 are two numpy arrays
    if(numpy.shape(data1)!=numpy.shape(data2)):
        print "\neuclidianDistance: data sets have different dimensions\n"
        sys.exit()
    else:
        z = (data1 - data2)*(data1 - data2)
        distance = numpy.sqrt(numpy.sum(z))
	
    if distance < 0:
        return [None]
    else:
        return [distance]
        
