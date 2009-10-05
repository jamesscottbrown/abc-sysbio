import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from pylab import *

import abc
import numpy
import sys


def printOptions():
    
    print "\nList of possible options:"
    print "-d\t--diagnostic\t no printing of diagnostic plots"
    print "-s\t--save\t\t no backup after each population"
    print "-t\t--timeseries\t no plotting of simulation results after each population"
    print "-p\t--plotdata\t no plotting of given data points""-p"
    print "-i\t--infile\t declaration of the input file. This input file has to be provided to run the program!"
    print "-y\t--synthetic\tPrint a synthetic data set at the timepoints given in the input file, using initial conditions and constant parameters."
    print "-S\t--simulate\tSimulate the model over the range of timepoints, using initial conditions and constant parameters."
    print "-f\t--fulloutput\tprint epsilon, sampling steps and acceptence rates after each population."
    print "-h\t--help\tPrint this list of options."

    print "\n"

    
def main():
    
    diagnostic=True
    pickling=True
    file_exist=False
    plot=True
    plotTimeSeries=True
    simulate=False
    synthetic=False
    full=False

    for i in range(1,len(sys.argv)):

        if sys.argv[i].startswith('--'):
            option = sys.argv[i][2:]
            if option == 'help':
                printOptions()
                sys.exit()
            elif option == 'diagnostic': diagnostic=False
            elif option == 'save': pickling=False
            elif option == 'timeseries': plotTimeSeries=False
            elif option == 'plotdata': plot=False
            elif option == 'simulate': simulate=True
            elif option == 'fulloutput': full=True
            elif option == 'synthetic':
                simulate=True
                synthetic=True
            elif option == 'infile': 
                input_file=sys.argv[i+1]
                file_exist=True
            elif not(sys.argv[i-1][2:] == 'infile'): 
                print "\nunknown option "+sys.argv[i]
                printOptions()
                sys.exit()
            
        
        elif sys.argv[i].startswith('-'):
            option = sys.argv[i][1:]
            if option == 'h':
                printOptions()
                sys.exit()
            elif option == 'd': diagnostic=False
            elif option == 's': pickling=False
            elif option == 't': plotTimeSeries=False
            elif option == 'p': plot=False
            elif option == 'S': simulate=True
            elif option == 'f': full=True
            elif option == 'y':
                simulate=True
                synthetic=True
            elif option == 'i': 
                input_file=sys.argv[i+1]
                file_exist=True
            elif not(sys.argv[i-1][2:] == 'i'): 
                print "\nunknown option "+sys.argv[i]
                printOptions()
                sys.exit()
        elif not((sys.argv[i-1][2:] == 'infile') or (sys.argv[i-1][1:] == 'i')):
            print "\nunknown expression \""+sys.argv[i]+"\""
            printOptions()
            sys.exit()
                
    if (file_exist==True and simulate==False):            
        restart,name,data,timepoints,numOutput,epsilon,InitValues,integrationType,modelWeight,prior,kernel,sampleFromPrior,source,fit,beta,dt, rtol, atol, constKernels=abc.parseInfo.getAlgorithmInfo(input_file)
        
        abc.importSBML(source,integrationType,name)
   
        if plot==True:
            plt.plot(timepoints,data,'o')
            plt.savefig('_data')
            cla()

    
#####run ABCSMC

        modelDistribution,inferred_parameter,rate=abc.abcSMC_model.abcSimulator(name,data, timepoints, numOutput, epsilon, InitValues,integrationType,modelWeight,prior,kernel,dt=dt,fit=fit,sampleFromPrior=sampleFromPrior,beta=beta,restart=restart,diagnostic=diagnostic,pickling=pickling,plotDataSeries=plotTimeSeries, full=full,constKernels=constKernels)

#####if -S or --simulate option applies, simulate data at timepoints in infile using initialValues in infile.
#####plots are made - do we want to do anythin else e.g. save the data to a text file?
        
    elif (file_exist==True and simulate==True):
        name, timepoints, InitValues, integrationType, ConstantParameters, source, dt, rtol, atol =abc.parseInfo.getInfoForSimulation(input_file)

        modelCorrect, message=abc.checkInputArguments.checkInputSimulation(name, timepoints, InitValues, integrationType, ConstantParameters, source, dt)

        if modelCorrect == False:
            print message
            sys.exit()

        abc.importSBML(source,integrationType,name)

        if synthetic!=True:
            timepoints=numpy.arange(timepoints[0], timepoints[-1]+dt, dt)

        for i in range(0, len(name)):
            data=abc.abcSMC_model.simulateData(name,i,integrationType,InitValues,timepoints,ConstantParameters[i],dt)

            if synthetic==True:
                print name[i]
                print data
                plt.plot(timepoints, data, 'o')
                plt.savefig("Synthetic_Data_"+name[i])
                cla()
            else: 
                plt.plot(timepoints, data)
                plt.savefig("Simulated_Data_"+name[i])
                cla()

    else:
        print "\nNo input_file is given!\nUse: \n\t-i 'inputfile' \nor: \n\t--infile 'inputfile' \n"
####plot data with final parameter
#abc.getResults.plotTimeSeries(name[0],inferred_parameter,integrationType[0],InitValues[0],timepoints,populationNumber=5,amount=10,model=1)
    

main()
