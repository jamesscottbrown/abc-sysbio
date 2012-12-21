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
        print "data 1",data1, "data2",data2
        print "\neuclidianDistance: data sets have different dimensions\n"
        sys.exit()
    else:
        z = (data1 - data2)*(data1 - data2)
        distance = numpy.sqrt(numpy.sum(z))

    if distance < 0:
        return [None]
    else:
        return [distance]



##########################################################
def euclidianLogDistance(data1, data2, parameters, model):
    """Returns the euclidian distance between the log  of two data sets.
    Data sets must have the same dimensions.
    ***** args *****
        data1, data2:
            numpy array objects with the same dimensions.
    """
    distance_log = -1
    #data1 and data2 are two numpy arrays
    if(numpy.shape(data1)!=numpy.shape(data2)):
        print "data 1",data1, "data2",data2
        print "\neuclidianDistance: data sets have different dimensions\n"
        sys.exit()
    else:
        z = (numpy.log(data1+numpy.power(10.0,-5)) - numpy.log(data2+numpy.power(10.0,-5))) * (numpy.log(data1+numpy.power(10.0,-5)) - numpy.log(data2+numpy.power(10.0,-5)))
        distance_log = numpy.sqrt(numpy.sum(z))
        #print data1

   
    if distance_log < 0:
        return [None]
    else:
        return [distance_log]
