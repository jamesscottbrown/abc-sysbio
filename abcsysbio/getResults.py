import re, math
import numpy

import matplotlib
from matplotlib.ticker import FormatStrFormatter
from pylab import *

from matplotlib.backends.backend_pdf import PdfPages

import abcsmc
from abcsmc import transform_data_for_fitting


# weighted histogramming
def bin_data(d, w, nbins):
    d_max = numpy.max(d)
    d_min = numpy.min(d) - 1e-6  # ensures that the lowest entry is included in the first bin
    bin_width = (d_max - d_min) / nbins

    bin_l = numpy.array([d_min + i * bin_width for i in range(nbins)])
    bin_u = numpy.array([d_min + (i + 1) * bin_width for i in range(nbins)])
    bin_c = numpy.array([bin_l[i] + bin_width / 2 for i in range(nbins)])

    count = numpy.zeros([nbins])

    for k in range(len(d)):
        kd = d[k]
        kw = w[k]

        i = 0
        for i in range(nbins):
            if bin_l[i] < kd <= bin_u[i]:
                count[i] = count[i] + kw
                break

    return [bin_c, count]


def MatrixToTextfile(matrix, filename, model, eps):
    """
    Write the part of a three dimensional matrix indexed by model and eps
    to a text file

    ***** args *****
    
    matrix:

            A three-dimensional matrix
            
    filename:

            A string, the name of the file to be written

    model:

            An integer
            
    eps:

            An integer
            
    """

    out_file = open(filename + ".txt", "w")
    for particle in range(len(matrix[model][eps][0])):
        for param in range(len(matrix[model][eps])):
            out_file.write(repr(matrix[model][eps][param][particle]) + " ")
        out_file.write("\n")
    out_file.close()


def printModelDistribution(matrix, eps, filename='model_distribution.txt'):
    """
    Write the contents of a two-dimensional matrix indexed by
    eps to a text file.
    
    ***** args *****

    matrix:

            A two-dimensional matrix

    eps:

            An integer

    ***** kwargs *****
    
    filename:

            String.
            The name of the text file to be written.
    """

    out_file = open(filename, "w")
    for j in range(eps + 1):
        for i in range(len(matrix[0])):
            out_file.write(repr(matrix[j][i]) + " ")
        out_file.write("\n")
    out_file.close()


def plotData(data, filename):
    matplotlib.pyplot.subplot(111)
    clf()
    matplotlib.pyplot.plot(data.timepoints, data.values, 'o')
    xlabel('time')
    ylabel('Unit')
    matplotlib.pyplot.savefig(filename)
    matplotlib.pylab.clf()


def plotTimeSeries2(model, pars, data, beta, filename, traj2, population, plotdata=True):
    """
    Plot simulated trajectories from the model with accepted parameters.
        
    ***** args *****
    
    model:

            model object

    pars:
            2D list of parameters for plotting to be done
    
    data:

            data object
   
    filename:

            Name of the output file to write to.

    plotdata:

            Boolean
            Whether or not to plot the data over the trajectories.

    """

    nsim = len(pars)
    matplotlib.pyplot.subplot(111)
    clf()
    for i in range(nsim):
        for j in range(beta):
            points_sim = traj2[i][j]
            matplotlib.pyplot.plot(data.timepoints, points_sim)

    if plotdata:
        matplotlib.pyplot.plot(data.timepoints, data.values, 'o')
        xlabel('time')
        ylabel('Unit')

    matplotlib.pyplot.savefig(filename)
    matplotlib.pylab.clf()


def plotTimeSeries(model, pars, data, beta, filename, plotdata=True):
    """
    Plot simulated trajectories from the model with accepted parameters.
        
    ***** args *****
    
    model:

            model object

    pars:
            2D list of parameters for plotting to be done
    
    data:

            data object
   
    filename:

            Name of the output file to write to.

    plotdata:

            Boolean
            Whether or not to plot the data over the trajectories.

    """

    nsim = len(pars)
    sims = model.simulate(pars, data.timepoints, nsim, beta=beta)

    matplotlib.pyplot.subplot(111)
    clf()
    for i in range(nsim):
        for j in range(beta):
            points = sims[i, j, :, :]
            points_sim = abcsmc.transform_data_for_fitting(model.fit, points)
            matplotlib.pyplot.plot(data.timepoints, points_sim)

    if plotdata:
        matplotlib.pyplot.plot(data.timepoints, data.values, 'o')
        xlabel('time')
        ylabel('Unit')

    matplotlib.pyplot.savefig(filename)
    matplotlib.pylab.clf()


def getAllHistograms(matrix, weights, population=1, PlotName='AllScatterPlots', model=1):
    """
    Plot weighted histograms.
    
    ***** args *****

    matrix:

            Matrix of data

    weights:

            Weights to assign to data

    *** kwargs *****

    populations:

            Integer, to index matrix with required population

    PlotName:

            String.
            Name for the saved plot

    model:

            Integer, to index matrix with required model.

    """

    matplotlib.pylab.clf()
    npar = len(matrix[int(model) - 1][0])

    # Maximum plots per page is 16
    # If we are over this then require multiple plots
    multi = False
    if npar > 16:
        multi = True

    if not multi:
        # print "******************* DOING SINGLE"
        # In the below max1 refers to the number of rows
        # and max2 refers to the number of columns ie max2 x max1

        # Lets check to see whether we have more than four parameters
        # If so we just make the plots 1 x npar

        max1 = 4.0
        if max1 > npar:
            max1 = npar
            max2 = 1

        dim = math.ceil(npar / max1)

        if dim == 1:
            max1 = 2
            max2 = 2
        elif max1 > dim:
            max2 = dim
        else:
            max2 = max1

        num_plots = math.ceil(dim / max2)

        for p in range(int(num_plots)):
            start = p * max1 ** 2
            end = p * max1 * max2 + max1 * max2
            for i in range(int(start), int(end)):
                if i >= len(matrix[int(model) - 1][0]):
                    break
                matplotlib.pyplot.subplot(max1, max2, i - start + 1)
                subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.5)
                x = matrix[int(model) - 1][int(population) - 1][i]
                w = weights[int(model) - 1][int(population) - 1][i]
                bins = 20.0
                if not (len(x) == 0):
                    cla()

                    histogram_x, histogram_y = bin_data(x, w, int(bins))

                    max_x = max(histogram_x)
                    min_x = min(histogram_x)
                    range_x = max_x - min_x

                    matplotlib.pyplot.bar(histogram_x, histogram_y, color='#1E90FF', width=range_x / bins, align='center')
                    xlabel('parameter ' + repr(i + 1), size='xx-small')

                    xmin, xmax = xlim()
                    ymin, ymax = ylim()

                    ax = gca()
                    ay = gca()

                    if (xmax - xmin) < 0.1 or (xmax - xmin) >= 1000:
                        x_formatter = FormatStrFormatter('%0.1e')
                    else:
                        x_formatter = FormatStrFormatter('%0.2f')
                    ax.xaxis.set_major_formatter(x_formatter)

                    y_formatter = FormatStrFormatter('%i')
                    axis([xmin, xmax, ymin, ymax])
                    yticks((ymin, (ymin + ymax) / 2.0, ymax), size='xx-small')
                    xticks((xmin, (xmin + xmax) / 2.0, xmax), size='xx-small')

            savefig(PlotName)
            matplotlib.pylab.clf()
            matplotlib.pyplot.subplot(111)

    else:

        s_num = 1
        p_num = 1
        for i in range(npar):
            matplotlib.pyplot.subplot(4, 4, s_num)
            subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.5)
            x = matrix[int(model) - 1][int(population) - 1][i]
            w = weights[int(model) - 1][int(population) - 1][i]

            cla()
            bins = 20
            histogram_x, histogram_y = bin_data(x, w, int(bins))

            max_x = max(histogram_x)
            min_x = min(histogram_x)
            range_x = max_x - min_x

            matplotlib.pyplot.bar(histogram_x, histogram_y, color='#1E90FF', width=range_x / bins, align='center')
            xlabel('parameter ' + repr(i + 1), size='xx-small')

            xmin, xmax = xlim()
            ymin, ymax = ylim()

            ax = gca()
            ay = gca()

            if (xmax - xmin) < 0.1 or (xmax - xmin) >= 1000:
                x_formatter = FormatStrFormatter('%0.1e')
            else:
                x_formatter = FormatStrFormatter('%0.2f')
            ax.xaxis.set_major_formatter(x_formatter)

            y_formatter = FormatStrFormatter('%i')
            axis([xmin, xmax, ymin, ymax])
            yticks((ymin, (ymin + ymax) / 2.0, ymax), size='xx-small')
            xticks((xmin, (xmin + xmax) / 2.0, xmax), size='xx-small')

            s_num += 1
            if s_num == 17:
                savefig(PlotName + '_' + repr(p_num))
                s_num = 1
                p_num += 1
                clf()

        savefig(PlotName + '_' + repr(p_num))
        matplotlib.pylab.clf()
        matplotlib.pyplot.subplot(111)


def getAllScatterPlots(matrix, weights, populations=(1,), PlotName='AllScatterPlots', model=1):
    """
    Plot scatter plots and histograms of data given in matrix.
    Used to plot posterior parameter distributions.
    ***** args *****

    matrix:

            Matrix of data.

    *** kwargs *****
    
    populations:

            Ordered tuple of integers.
            Determines which data will be plotted.
            Used to index the matrix.
            
    PlotName:

            String
            Name for the saved plot

    model:

            Integer
            Determines which data will be plotted.
            Used to index the matrix.

    """

    matplotlib.pylab.clf()
    dim = len(matrix[int(model) - 1][0])

    my_colors = ['#000000', '#003399', '#3333FF', '#6666FF', '#990000', '#CC0033', '#FF6600', '#FFCC00', '#FFFF33',
                 '#33CC00', '#339900', '#336600']

    if len(populations) > len(my_colors):
        q = int(math.ceil(len(populations) / len(my_colors)))

        for slopes in range(q):
            my_colors.extend(my_colors)

    max1 = 4.0
    if dim <= max1:
        permutation = numpy.zeros([dim ** 2, 2])
        k = 0
        for i in range(1, dim + 1):
            for j in range(1, dim + 1):
                permutation[k][0] = i
                permutation[k][1] = j
                k += 1

        bin_b = 20.0
        i2 = 0
        for i in range(len(permutation)):
            matplotlib.pyplot.subplot(dim, dim, i + 1)
            subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.5)
            w = weights[int(model) - 1][int(populations[len(populations) - 1]) - 1][int(permutation[i][0]) - 1]

            for j in range(len(populations)):

                x = matrix[int(model) - 1][int(populations[j]) - 1][int(permutation[i][0]) - 1]
                y = matrix[int(model) - 1][int(populations[j]) - 1][int(permutation[i][1]) - 1]
                g = (j + 1) * ((len(populations) * 1.5) ** (-1))
                if permutation[i][0] == permutation[i][1]:
                    if j == len(populations) - 1:
                        cla()
                        if not (len(x) == 0):
                            i2 += 1
                            histogram_x, histogram_y = bin_data(x, w, int(bin_b))

                            max_x = max(histogram_x)
                            min_x = min(histogram_x)
                            range_x = max_x - min_x
                            matplotlib.pyplot.bar(histogram_x, histogram_y, width=range_x / bin_b, color=my_colors[j],
                                                  align='center')
                            xlabel('parameter ' + repr(i2), size='xx-small')

                else:
                    if not (len(x) == 0):
                        scatter(x, y, s=10, marker='o', c=my_colors[j], edgecolor=my_colors[j])
                        ylabel('parameter ' + repr(int(permutation[i][1])), size='xx-small')
                xlabel('parameter ' + repr(int(permutation[i][0])), size='xx-small')

                xmin, xmax = xlim()
                ymin, ymax = ylim()

                ax = gca()
                ay = gca()
                if (xmax - xmin) < 0.1 or (xmax - xmin) >= 1000:
                    x_formatter = FormatStrFormatter('%0.1e')
                else:
                    x_formatter = FormatStrFormatter('%0.2f')
                ax.xaxis.set_major_formatter(x_formatter)

                if (ymax - ymin) < 0.1 or (ymax - ymin) >= 1000:
                    y_formatter = FormatStrFormatter('%0.1e')
                else:
                    y_formatter = FormatStrFormatter('%0.2f')
                ay.yaxis.set_major_formatter(y_formatter)

                if permutation[i][0] == permutation[i][1]:
                    y_formatter = FormatStrFormatter('%i')

                axis([xmin, xmax, ymin, ymax])
                xticks((xmin, (xmin + xmax) / 2.0, xmax), size='xx-small')
                yticks((ymin, (ymin + ymax) / 2.0, ymax), size='xx-small')

        savefig(PlotName)
        matplotlib.pylab.clf()
        matplotlib.pyplot.subplot(111)

    else:
        i2 = 0
        permutation = zeros([dim * (dim + 1) / 2, 2])
        k = 0
        for i in range(1, dim + 1):
            for j in range(i, dim + 1):
                permutation[k][0] = i
                permutation[k][1] = j
                k += 1

        bin_b = 20.0

        dim = math.ceil(len(permutation) / max1)
        num_plots = math.ceil(dim / max1)

        for p in range(int(num_plots)):

            start = p * max1 ** 2
            end = p * max1 ** 2 + max1 ** 2
            for i in range(int(start), int(end)):
                if i >= len(permutation):
                    break
                matplotlib.pyplot.subplot(max1, max1, i - start + 1)
                subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.5)
                w = weights[int(model) - 1][int(populations[len(populations) - 1]) - 1][int(permutation[i][0]) - 1]
                for j in range(len(populations)):
                    x = matrix[int(model) - 1][int(populations[j]) - 1][int(permutation[i][0]) - 1]
                    y = matrix[int(model) - 1][int(populations[j]) - 1][int(permutation[i][1]) - 1]
                    g = (j + 1) * ((len(populations) * 1.5) ** (-1))
                    if permutation[i][0] == permutation[i][1]:
                        if j == len(populations) - 1:
                            cla()
                            if not (len(x) == 0):
                                i2 += 1

                                histogram_x, histogram_y = bin_data(x, w, int(bin_b))

                                max_x = max(histogram_x)
                                min_x = min(histogram_x)
                                range_x = max_x - min_x
                                matplotlib.pyplot.bar(histogram_x, histogram_y, color=my_colors[j], width=range_x / bin_b,
                                                      align='center')
                                xlabel('parameter ' + repr(i2), size='xx-small')

                    else:
                        if not (len(x) == 0):
                            scatter(x, y, s=10, marker='o', c=my_colors[j], edgecolor=my_colors[j])
                            ylabel('parameter ' + repr(int(permutation[i][1])), size='xx-small')
                            xlabel('parameter ' + repr(int(permutation[i][0])), size='xx-small')

                    xmin, xmax = xlim()
                    ymin, ymax = ylim()

                    ax = gca()
                    ay = gca()
                    if (xmax - xmin) < 0.1 or (xmax - xmin) >= 1000:
                        x_formatter = FormatStrFormatter('%0.1e')
                    else:
                        x_formatter = FormatStrFormatter('%0.2f')
                    ax.xaxis.set_major_formatter(x_formatter)

                    if (ymax - ymin) < 0.1 or (ymax - ymin) >= 1000:
                        y_formatter = FormatStrFormatter('%0.1e')
                    else:
                        y_formatter = FormatStrFormatter('%0.2f')
                    ay.yaxis.set_major_formatter(y_formatter)

                    if permutation[i][0] == permutation[i][1]:
                        y_formatter = FormatStrFormatter('%i')

                    axis([xmin, xmax, ymin, ymax])
                    xticks((xmin, (xmin + xmax) / 2.0, xmax), size='xx-small')
                    yticks((ymin, (ymin + ymax) / 2.0, ymax), size='xx-small')

            savefig(PlotName + "_" + repr(p))
            matplotlib.pylab.clf()
            matplotlib.pyplot.subplot(111)


def getScatterPlot(matrix, parameter, populations=(1,), PlotName='ScatterPlot', model=1):
    """
    Plot a single scatter plot of accepted parameters.
    ***** args *****

    matrix:

            Matrix of accepted parameters.

    parameter:

            List of integers of length 2.
            Used to index the matrix,
            to determine which two parameters to plot against each other.

    ***** kwargs *****
    
    populations:

            Tuple of integers.
            Used to index the matrix.
            Accepted parameters from these populations will be plotted.

    PlotName:

            String
            Name for the saved plot

    model:

            Integer
            Number for the model from which accepted parameters should be plotted.

    """

    matplotlib.pylab.clf()
    matplotlib.pyplot.subplot(111)

    for j in range(len(populations)):
        g = (j + 1) * ((len(populations) * 1.5) ** (-1))
        x = matrix[int(model) - 1][int(populations[j]) - 1][int(parameter[0]) - 1]
        y = matrix[int(model) - 1][int(populations[j]) - 1][int(parameter[1]) - 1]
        if not (len(x) == 0):
            scatter(x, y, s=10, c=repr(g), edgecolor=repr(g))
    savefig(PlotName)
    matplotlib.pylab.clf()


def getModelDistribution(matrix, epsilon, rate, PlotName='ModelDistribution'):
    """
    Plot a histogram of the posterior distributions of the models
    ***** args *****

    matrix:

            Matrix containing the model distributions
            after each population.

    epsilon:

            Epsilon used for each population
            (To be displayed above the plot)

    rate:

            Acceptance rate for each population
            (To be displayed above the plot)

    ***** kwargs *****

    PlotName:

            Prefix for file name for the saved plot.
            If the saved plot runs over multiple pages,
            each page will be saved to a separate file, with this prefix.

    """

    matplotlib.pylab.clf()
    max1 = 4.0  # must be float or double, but no integer
    if max1 > len(matrix):
        max1 = len(matrix)
        max2 = 1

    dim = math.ceil(len(matrix) / max1)

    if dim == 1:
        max1 = 2
        max2 = 2
    elif max1 > dim:
        max2 = dim
    else:
        max2 = max1

    num_plots = math.ceil(dim / max2)

    for p in range(int(num_plots)):

        start = p * max1 ** 2
        end = p * max1 * max2 + max1 * max2
        for i in range(int(start), int(end)):
            if i >= len(matrix): 
                break
            matplotlib.pyplot.subplot(max1, max2, i - start + 1)
            subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.8)
            left = arange(1, matrix.shape[1] + 1, 1)
            height = matrix[i]
            bar(left, height, width=1.0, color='#1E90FF', align='center')
            xmin = 0
            xmax = matrix.shape[1] + 1
            ymin, ymax = ylim()
            axis([math.floor(xmin), math.ceil(xmax), 0, ymax + ymax * 0.1])
            yticks(size='xx-small')
            xticks(left, size='xx-small')
            title("(" + repr(i + 1) + ") " + str(epsilon[i]) + "\n" + str(rate[i]), size='xx-small')
            savefig(PlotName + '_' + repr(p + 1))
        matplotlib.pylab.clf()
        matplotlib.pyplot.subplot(111)
