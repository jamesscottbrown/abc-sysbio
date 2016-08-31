import numpy
import sys


def euclidian_distance(data1, data2, parameters, model):
    """
    Returns the Euclidian distance between two data sets (data1 and data2), which must have the same dimensions.

    The parameters and model arguments are ignored: they are included as all distance functions have the same signature,
    and custom functions may need this information.

    Parameters
    ----------
    data1 - first dataset
    data2 - second dataset
    parameters - not used
    model - not used


    """

    # data1 and data2 are two numpy arrays
    if numpy.shape(data1) != numpy.shape(data2):
        print "\neuclidian_distance: data sets have different dimensions\n"
        sys.exit()
    else:
        z = (data1 - data2) * (data1 - data2)
        distance = numpy.sqrt(numpy.sum(z))

    if distance < 0:
        return [None]
    else:
        return [distance]
