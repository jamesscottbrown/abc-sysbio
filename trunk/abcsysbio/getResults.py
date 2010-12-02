import re, math
import numpy

import matplotlib
from matplotlib.ticker import FormatStrFormatter
from pylab import *

# need to get howToFitData
import abcsysbio.abcSMC_model

# weighted histogramming 
def bin_data(d, w, nbins):

    d_max = numpy.max(d)
    d_min = numpy.min(d) - 1e-6 # ensures that the lowest entry is included in the first bin
    bin_width = (d_max - d_min)/nbins

    bin_l = numpy.array( [ d_min + i*bin_width for i in range(nbins) ] )
    bin_u = numpy.array( [ d_min + (i+1)*bin_width for i in range(nbins) ] )
    bin_c = numpy.array( [ bin_l[i] + bin_width/2 for i in range(nbins) ] )

    count = numpy.zeros( [ nbins ] )

    for k in range(len(d)):
        kd = d[k]
        kw = w[k]

        i = 0
        for i in range(nbins):
            if kd > bin_l[i] and kd <= bin_u[i]:
                count[i] = count[i] + kw
                break

    return [ bin_c, count ]


def MatrixToTextfile(matrix,filename,model,eps):

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
    
    out_file=open(filename+".txt","w")
    for j in range(0,len(matrix[model][eps][0])):###for each parameter of the choosen model
        for i in range(0,len(matrix[model][eps])):###for each value of the parameter of the choosen model
            out_file.write(repr(matrix[model][eps][i][j])+" ")
        out_file.write("\n") 
    out_file.close()



def printModelDistribution(matrix,eps,filename='model_distribution.txt'):

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
    

    out_file=open(filename,"w")
    for j in range(0, eps+1):
        for i in range(0,len(matrix[0])):
            out_file.write(repr(matrix[j][i])+" ")
        out_file.write("\n")
    out_file.close()



def plotTimeSeries(ModelName,matrix,integrationType,InitValues,timepoints,fit,populationNumber=1,dt=0.1,amount=1,model=1,filename='TimeSeries',data=None,plotdata=False):

    """
    Plot simulated trajectories from the model with accepted parameters.
        
    ***** args *****
    
    ModelName:

            Name of the Python file representing the model.
    
    matrix:

            Matrix of accepted parameters

    integrationType:

            Integration type for simulating the model.
            One of ODE, SDE or Gillespie

    InitValues:

            Initial values for the simulation

    timepoints:

            Times at which the simulation should be made.

    fit:

            Fitting instructions for the trajectories from the model.

    ***** kwargs *****
    
    populationNumber:

            Index for the matrix, which population to plot
          
    dt:

            Internal timestep for the simulation.

    amount:

            Number of trajectories to plot.

    model:

            Index for the matrix, which model to plot parameters for.

    filename:

            Name of the output file to write to.

    data:

            Data to plot over the trajectories.

    plotdata:

            Boolean
            Whether or not to plot the data over the trajectories.

    """
        
    module = __import__(ModelName)

    o=re.compile('ODE')
    s=re.compile('SDE')
    g=re.compile('Gillespie')

    matplotlib.pyplot.subplot(111)
    
    time_start=timepoints[0]
    time_end=timepoints[-1]
    step=(time_end-time_start)/500

    timepoints_new=arange(time_start,time_end+step,step)

    clf()
    for i in range(0,amount):
        parameter=[]
        for j in range(0,len(matrix[int(model)][populationNumber-1])):
            parameter.append(matrix[int(model)][int(populationNumber)-1][j][i])

        if o.match(integrationType):
            from abcsysbio import abcodesolve
            data_sim=abcodesolve.abcodeint(module, InitValues, timepoints_new,parameter, dt)
        if s.match(integrationType):
            from abcsysbio import sdeint
            data_sim=sdeint.sdeint(module,InitValues,parameter,timepoints_new,dt=0.01)
        if g.match(integrationType):
            from abcsysbio import GillespieAlgorithm    
            data_sim=GillespieAlgorithm.GillespieInt(module, InitValues,parameter, timepoints_new)

        points_sim = abcsysbio.abcSMC_model.howToFitData(fit,data_sim)

        matplotlib.pyplot.plot(timepoints_new,points_sim)
    if plotdata==True: 
        matplotlib.pyplot.plot(timepoints,data,'o')
        xlabel('time')
        ylabel('Unit')
   
    
    matplotlib.pyplot.savefig(filename)
    matplotlib.pylab.clf()

def getAllHistograms(matrix,weights,population=1,PlotName='AllScatterPlots', model=1):

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
    dim=len(matrix[int(model)-1][0])

    max1=4.0 #must be float or double, but no integer
    if max1>len(matrix[int(model)-1][0]): 
        max1=len(matrix[int(model)-1][0])
        max2=1

    dim=math.ceil(len(matrix[int(model)-1][0])/max1)

    if dim == 1:
        max1=2
        max2=2
    elif (max1 > dim):
        max2=dim
    else: max2=max1

    numOfPlots=math.ceil(dim/max2)

    for p in range(0,int(numOfPlots)):
        start=p*max1**2
        end=p*max1*max2+max1*max2
        for i in range(int(start),int(end)):
            if i>=len(matrix[int(model)-1][0]): break
            matplotlib.pyplot.subplot(max1,max2,i-start+1)
            subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=0.6, hspace=0.5)    
            x=matrix[int(model)-1][int(population)-1][i]
            w=weights[int(model)-1][int(population)-1][i]
            bins=20.0
            if not(len(x)==0):
                cla()
                #hist=Histogram.WeightedHistogram(x,w,int(bins))
                #length=hist.__len__()
                #histogramX=[]
                #histogramY=[]
                #for k in range(0, length):
                #    histogramX.append(hist.__getitem__(k)[0])
                #    histogramY.append(hist.__getitem__(k)[1])

                histogramX=[]
                histogramY=[]
                histogramX, histogramY = bin_data(x, w, int(bins))
                
                maxX=max(histogramX)
                minX=min(histogramX)
                xrange=maxX-minX
                
                matplotlib.pyplot.bar(histogramX,histogramY,color='#1E90FF',width=xrange/bins,align='center')
                xlabel('parameter '+repr(i+1),size='xx-small')
            
                xmin,xmax=xlim()
                ymin,ymax=ylim()
                
                ax = gca()
                ay = gca()

                if ((xmax-xmin)<0.1 or (xmax-xmin)>=1000): xFormatter = FormatStrFormatter('%0.1e')
                else: xFormatter = FormatStrFormatter('%0.2f')
                ax.xaxis.set_major_formatter(xFormatter)    


                yFormatter = FormatStrFormatter('%i')    
                axis([xmin, xmax, ymin, ymax])
                yticks((ymin, (ymin+ymax)/2.0, ymax),size = 'xx-small')
                xticks((xmin, (xmin+xmax)/2.0, xmax),size='xx-small')




        savefig(PlotName)
        matplotlib.pylab.clf()
        matplotlib.pyplot.subplot(111)



def getAllScatterPlots(matrix,weights,populations=(1,),PlotName='AllScatterPlots', model=1):

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
    dim=len(matrix[int(model)-1][0])

  
    myColors=['#000000','#003399','#3333FF','#6666FF','#990000','#CC0033','#FF6600','#FFCC00','#FFFF33','#33CC00','#339900','#336600']

    if len(populations)>len(myColors):
        q=int(math.ceil(len(populations)/len(myColors)))
        
        for slopes in range(0,q):
            myColors.extend(myColors)

    max1=4.0
    if dim<=max1:
        permutation=numpy.zeros([dim**2,2])
        k=0
        for i in range(1,dim+1):
            for j in range(1,dim+1):
                permutation[k][0]=i
                permutation[k][1]=j
                k=k+1

        binB=20.0
        i2=0
        for i in range(0,len(permutation)):
            matplotlib.pyplot.subplot(dim,dim,i+1)
            subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=0.6, hspace=0.5)
            w=weights[int(model)-1][int(populations[len(populations)-1])-1][int(permutation[i][0])-1]
          
            for j in range(0,len(populations)):


                x=matrix[int(model)-1][int(populations[j])-1][int(permutation[i][0])-1]
                y=matrix[int(model)-1][int(populations[j])-1][int(permutation[i][1])-1]
                g=(j+1)*((len(populations)*1.5)**(-1))
                if permutation[i][0]==permutation[i][1]:
                    if j==len(populations)-1:
                        cla()
                        if not(len(x)==0):
                            i2=i2+1
                            #hist=Histogram.WeightedHistogram(x,w,int(binB))
                            #length=hist.__len__()
                            #histogramX=[]
                            #histogramY=[]

                            #for k in range(0, length):
                            #    histogramX.append(hist.__getitem__(k)[0])
                            #    histogramY.append(hist.__getitem__(k)[1])

                            histogramX=[]
                            histogramY=[]
                            histogramX, histogramY = bin_data(x, w, int(binB))
                            
                            maxX=max(histogramX)
                            minX=min(histogramX)
                            xrange=maxX-minX
                            matplotlib.pyplot.bar(histogramX,histogramY,width=xrange/binB,color=myColors[j],align='center')
                            xlabel('parameter '+repr(i2),size='xx-small')
                 
                else:
                    if not(len(x)==0):
                        scatter(x,y,s=10,marker='o',c=myColors[j],edgecolor=myColors[j])
                        ylabel('parameter '+repr(int(permutation[i][1])),size='xx-small')
                xlabel('parameter '+repr(int(permutation[i][0])),size='xx-small')
                
                xmin,xmax=xlim()
                ymin,ymax=ylim()

                ax = gca()
                ay = gca()
                if ((xmax-xmin)<0.1 or (xmax-xmin)>=1000): xFormatter = FormatStrFormatter('%0.1e')
                else: xFormatter = FormatStrFormatter('%0.2f')
                ax.xaxis.set_major_formatter(xFormatter)
                
                if ((ymax-ymin)<0.1 or (ymax-ymin)>=1000): yFormatter = FormatStrFormatter('%0.1e')
                else: yFormatter = FormatStrFormatter('%0.2f')
                ay.yaxis.set_major_formatter(yFormatter)

                if permutation[i][0]==permutation[i][1]: yFormatter = FormatStrFormatter('%i')
                
                axis([xmin, xmax, ymin, ymax])
                xticks((xmin, (xmin+xmax)/2.0, xmax),size='xx-small')
                yticks((ymin, (ymin+ymax)/2.0, ymax),size = 'xx-small')
       
               

               
        savefig(PlotName)
        matplotlib.pylab.clf()
        matplotlib.pyplot.subplot(111)

    else:
        i2=0
        permutation=zeros([dim*(dim+1)/2,2])
        k=0
        for i in range(1,dim+1):
            for j in range(i,dim+1):
                permutation[k][0]=i
                permutation[k][1]=j
                k=k+1
        
        binB=20.0

        dim=math.ceil(len(permutation)/max1)
        numOfPlots=math.ceil(dim/max1)
       
        for p in range(0,int(numOfPlots)):

            start=p*max1**2
            end=p*max1**2+max1**2
            for i in range(int(start),int(end)):
                if i>=len(permutation): break
                matplotlib.pyplot.subplot(max1,max1,i-start+1)
                subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=0.6, hspace=0.5)
                w=weights[int(model)-1][int(populations[len(populations)-1])-1][int(permutation[i][0])-1]
                for j in range(0,len(populations)):
                    x=matrix[int(model)-1][int(populations[j])-1][int(permutation[i][0])-1]
                    y=matrix[int(model)-1][int(populations[j])-1][int(permutation[i][1])-1]
                    g=(j+1)*((len(populations)*1.5)**(-1))
                    if permutation[i][0]==permutation[i][1]:
                        if j==len(populations)-1:
                            cla()
                            if not(len(x)==0):
                                i2=i2+1
                                #hist=Histogram.WeightedHistogram(x,w,int(binB))
                                #length=hist.__len__()
                                #histogramX=[]
                                #histogramY=[]
                                #for k in range(0, length):
                                #    histogramX.append(hist.__getitem__(k)[0])
                                #    histogramY.append(hist.__getitem__(k)[1])

                                histogramX=[]
                                histogramY=[]
                                histogramX, histogramY = bin_data(x, w, int(binB))
                            
                                maxX=max(histogramX)
                                minX=min(histogramX)
                                xrange=maxX-minX
                                matplotlib.pyplot.bar(histogramX,histogramY,color=myColors[j],width=xrange/binB,align='center')
                                xlabel('parameter '+repr(i2),size='xx-small')

                    #        matplotlib.pyplot.hist(matrix[int(model)-1][j][int(permutation[i][0])-1],bins=20,facecolor='white')
                    #        xlabel('parameter '+repr(int(permutation[i][0])),size='xx-small')
                    
                    else:
                        if not(len(x)==0):
                            scatter(x,y,s=10,marker='o',c=myColors[j],edgecolor=myColors[j])
                            ylabel('parameter '+repr(int(permutation[i][1])),size='xx-small')
                            xlabel('parameter '+repr(int(permutation[i][0])),size='xx-small')
                
                    xmin,xmax=xlim()
                    ymin,ymax=ylim()
                   
                    ax = gca()
                    ay = gca()
                    if ((xmax-xmin)<0.1 or (xmax-xmin)>=1000): xFormatter = FormatStrFormatter('%0.1e')
                    else: xFormatter = FormatStrFormatter('%0.2f')
                    ax.xaxis.set_major_formatter(xFormatter)

                    if ((ymax-ymin)<0.1 or (ymax-ymin)>=1000): yFormatter = FormatStrFormatter('%0.1e')
                    else: yFormatter = FormatStrFormatter('%0.2f')
                    ay.yaxis.set_major_formatter(yFormatter)

                    if permutation[i][0]==permutation[i][1]: yFormatter = FormatStrFormatter('%i')
                                 
                    axis([xmin, xmax, ymin, ymax])
                    xticks((xmin, (xmin+xmax)/2.0, xmax),size='xx-small')
                    yticks((ymin, (ymin+ymax)/2.0, ymax),size = 'xx-small')
       
               
            savefig(PlotName+"_"+repr(p))
            matplotlib.pylab.clf()
            matplotlib.pyplot.subplot(111)



def getScatterPlot(matrix,parameter,populations=(1,),PlotName='ScatterPlot',model=1):

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

    for j in range(0,len(populations)):
        g=(j+1)*((len(populations)*1.5)**(-1))
        x=matrix[int(model)-1][int(populations[j])-1][int(parameter[0])-1]
        y=matrix[int(model)-1][int(populations[j])-1][int(parameter[1])-1]   
        if not(len(x)==0):
            scatter(x,y,s=10,c=repr(g),edgecolor=repr(g))
    savefig(PlotName)
    matplotlib.pylab.clf()



def getModelDistribution(matrix,epsilon,rate,PlotName='ModelDistribution'):

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
    max1=4.0 #must be float or double, but no integer
    if max1>len(matrix): 
        max1=len(matrix)
        max2=1

    dim=math.ceil(len(matrix)/max1)

    if dim == 1:
        max1=2
        max2=2
    elif (max1 > dim):
        max2=dim
    else: max2=max1

    numOfPlots=math.ceil(dim/max2)

    for p in range(0,int(numOfPlots)):
        
        start=p*max1**2
        end=p*max1*max2+max1*max2
        for i in range(int(start),int(end)):
            if i>=len(matrix): break
            matplotlib.pyplot.subplot(max1,max2,i-start+1)
            subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=0.6, hspace=0.8)
            left=arange(1,matrix.shape[1]+1,1)
            height=matrix[i]
            bar(left,height,width=1.0,color='#1E90FF',align='center')
            xmin=0
            xmax=matrix.shape[1]+1
            ymin,ymax=ylim()
            axis([math.floor(xmin), math.ceil(xmax),0,ymax + ymax*0.1])
            yticks(size='xx-small')
            xticks(left,size='xx-small')
            #yticks(arange(0,ymax+2,ymax*0.4),size='xx-small')

            title("("+repr(i+1)+") "+str(epsilon[i])+"\n"+str(rate[i]),size='xx-small')
            savefig(PlotName+'_'+repr(p+1))
        matplotlib.pylab.clf()
        matplotlib.pyplot.subplot(111)
