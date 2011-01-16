import os, sys, pickle
import numpy

from getResults import getAllScatterPlots
from getResults import getAllHistograms
from getResults import plotTimeSeries
from getResults import getModelDistribution

class input_output:
    
    def __init__(self, folder, restart):
        self.folder = folder

        if restart==True:
            self.folder = self.folder + '_restart'

    ################write rates, distances, trajectories    
    def write_data(self, population, results, timing, models, diagnostic=True):
        # results abcsmc_results class
        rate_file = open(self.folder+'/rates.txt',"a")
        rate_file.write(repr(population+1)+"\t"+repr(results.sampled)+"\t"+repr(results.rate))
        rate_file.write("\t" + str(round(timing,2)) + " secs\n")
        rate_file.close()


        distance_file=open(self.folder + '/distance_Population' + repr(population+1) + '.txt',"a")
        for i in range(len(results.distances)): 
            for j in range(len(results.distances[i])):
                print >>distance_file, i+1, j, results.distances[i][j], results.models[i]
        distance_file.close()

        traj_file=open(self.folder + '/traj_Population' + repr(population+1) + '.txt',"a")
        for i in range(len(results.trajectories)):
            arr = results.trajectories[i]
            nrow, ncol = numpy.shape( arr )
            for ic in range(ncol):
                for ir in range(nrow):
                    print >>traj_file, arr[ir,ic],
                print >>traj_file, ""
        traj_file.close()

        model_file = open(self.folder + '/ModelDistribution.txt',"a")
        for m in results.margins: 
            print >>model_file, m,
        model_file.close()

        nmodels = len(models)
        for mod in range(nmodels):
            print self.folder + '/results_' + models[mod].name
            try:
                print self.folder + '/results_' + models[mod].name
                os.chdir(self.folder + '/results_' + models[mod].name)
                print os.getcwd()
                os.mkdir("Population_" + repr(population+1))
                os.chdir("../..")
            except:
                print "\nCan not create the folder Population_"+repr(population+1)+"!\n"
                sys.exit()

        #if nmodels>1 and diagnostic==True:
        #    getModelDistribution(results.margins,epsilon,rate,PlotName=folder+'/ModelDistribution')

        population_mod=[]
        weights_mod=[]
        for mod in range(nmodels):
            weight_file=open(self.folder+'/results_'+ models[mod].name + '/Population_'+repr(population+1)+'/data_Weights'+repr(population+1)+".txt","w")
            param_file = open(self.folder+'/results_'+ models[mod].name + '/Population_'+repr(population+1)+'/data_Population'+repr(population+1)+".txt","w") 
            PlotName = self.folder+'/results_' + models[mod].name + '/Population_'+repr(population+1) + '/ScatterPlots_Population' + repr(population+1)
            PlotName2 = self.folder+'/results_' + models[mod].name + '/Population_'+repr(population+1) + '/weightedHistograms_Population' + repr(population+1)

            nparticles = len(results.weights)
            for g in range(nparticles):
                print g
                if( results.models[g] == mod ):
                    print g, results.models[g]

                    for k in range(len(results.parameters[g])):
                        print >>param_file, results.parameters[g][k],
                    print >>param_file, ""

            weight_file.close()
            param_file.close()
            
            #if diagnostic==True:
            #    population_mod.append([])
            #    weights_mod.append([])
            #    for eps in range(len(parameters[0])):
            #        population_mod[mod].append([])
            #        weights_mod[mod].append([])

            #getAllScatterPlots(population_mod,weights_mod,populations=numpy.arange(1,t+2,1),PlotName=PlotName,model=mod+1)
            #getAllHistograms(population_mod,weights_mod,population=t+1,PlotName=PlotName2, model=mod+1)
            

    ################make some nice plots
    def make_analysis_plots(folder, diagnostic, plotDataSeries, timepoints, fit, dt, data, t, nmodel, ModelName, integrationType,
                            InitValues, priors, population, numbers, modelDistribution, weights, epsilon, rate ):            

        if plotDataSeries==True:
            for mod in range(0,nmodel):
                x=2
                if numbers[mod][t]<x: x=int(numbers[mod][t])
                plotTimeSeries(ModelName[mod],
                               population,
                               integrationType[mod],
                               InitValues[mod],
                               timepoints,
                               fit[mod],
                               populationNumber=t+1,
                               dt=dt,
                               amount=x,
                               model=mod,
                               filename=folder+'/results_'+ModelName[mod]+'/Population_'+repr(t+1)+'/Timeseries_Population'+repr(t+1),
                               data=data,
                               plotdata=True)

    ################create output folders
    def create_output_folders(self, modelnames, numOutput, pickling):
        try:
            os.mkdir(self.folder)
            os.chdir(self.folder)
        except:
            print "\nThe folder "+ self.folder +" already exists!\n"
            sys.exit()

        for mod in modelnames:
            try:
                os.mkdir('results_' + mod)
            except:
                print "\nThe folder "+ self.folder + "/results_" + mod +" already exists!\n"
                sys.exit()

        os.chdir('..')

        if pickling==True:
            try:
                os.chdir(self.folder)
                os.mkdir('copy')
                os.chdir('..')
            except:
                print "\nThe folder \'copy\' already exists!\n"
                sys.exit()
    
            out_file=open(self.folder+'/copy/algorithm_parameter.dat',"w")
            pickle.dump(numOutput, out_file)
            out_file.close()
    ###

    ################read the stored data
    def read_pickled(self, location):
        # pickle numbers selected model of previous population
        # pickle population of selected model of previous population pop_pickled[selected_model][n][vectpos]
        # pickle weights of selected model of previous population weights_pickled[selected_model][n][vectpos]

        try:
            in_file=open(location + '/copy/model_last.dat',"r")
            model_pickled=pickle.load(in_file)
            in_file.close()
        except:
            print "\nCan not find file \'model_last.dat\' in folder \'copy\'!\n"
            sys.exit()

        try:
            in_file=open(location + '/copy/weights_last.dat',"r")
            weights_pickled=pickle.load(in_file)
            in_file.close()
        except:
            print "\nCan not find file \'weights_last.dat\' in folder \'copy\'!\n"
            sys.exit()

        try:
            in_file=open(location + '/copy/params_last.dat',"r")
            parameters_pickled=pickle.load(in_file)
            in_file.close()
        except:
            print "\nCan not find file \'params_last.dat\' in folder \'copy\'!\n"
            sys.exit()

        try:
            in_file=open(location + '/copy/margins_last.dat',"r")
            margins_pickled=pickle.load(in_file)
            in_file.close()
        except:
            print "\nCan not find file \'margins_last.dat\' in folder \'copy\'!\n"
            sys.exit()
    
        try:
            in_file=open(location + '/copy/kernels_last.dat',"r")
            kernel=pickle.load(in_file)
            in_file.close()
        except:
            print "\nCan not find file \'kernels_last.dat\' in folder \'copy\'!\n"
            sys.exit()

        #print "\n\n\n Reading previous population"
        #print "model_pickled", model_pickled, "\n\n\n"
        #print "weights_pickled", weights_pickled, "\n\n\n"
        #print "parameters_pickled", parameters_pickled, "\n\n\n"
        #print "margins_pickled", margins_pickled, "\n\n\n"

        return [model_pickled, weights_pickled, parameters_pickled, margins_pickled, kernel]
    ###

    ################write the stored data
    def write_pickled(self, nmodel, model_prev, weights_prev, parameters_prev, margins_prev, kernel):

        out_file=open(self.folder+'/copy/model_last.dat',"w")
        x = model_prev[:]
        pickle.dump(x,out_file)
        out_file.close()

        out_file=open(self.folder+'/copy/weights_last.dat',"w")
        x = weights_prev[:]
        pickle.dump(x,out_file)
        out_file.close()

        out_file=open(self.folder+'/copy/params_last.dat',"w")
        x= parameters_prev
        pickle.dump(x,out_file)
        out_file.close()
    
        out_file=open(self.folder+'/copy/margins_last.dat',"w")
        x = margins_prev[:]
        pickle.dump(x,out_file)
        out_file.close()
 
        out_file=open(self.folder+'/copy/kernels_last.dat',"w")
        x=[]
        for mod in range(0,nmodel):
            x.append(kernel[mod])
        pickle.dump(x,out_file)
        out_file.close()
    ###

    
