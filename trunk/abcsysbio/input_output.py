import os, sys, pickle
import numpy

from getResults import getAllScatterPlots
from getResults import getAllHistograms
from getResults import plotTimeSeries2
from getResults import getModelDistribution
from getResults import plotData

class input_output:
    
    def __init__(self, folder, restart, diagnostic, plotDataSeries, havedata=True):
        self.folder = folder
        self.diagnostic = diagnostic
        self.plotDataSeries = plotDataSeries
        self.havedata = havedata

        # Hold all data here for plotting purposes.
        # May want to remove this as could get large
        self.all_results = []
        
        if restart==True:
            self.folder = self.folder + '_restart'

    def plot_data(self, data):
        if self.havedata == True: plotData( data, self.folder+'/_data' )


    ################write rates, distances, trajectories    
    def write_data(self, population, results, timing, models, data):

        # results abcsmc_results class
        self.all_results.append( results )

        beta = len(results.trajectories[0])
        
        rate_file = open(self.folder+'/rates.txt',"a")
        print >>rate_file, population+1, results.epsilon, results.sampled, results.rate, round(timing,2)
        rate_file.close()

        # distances are stored as [nparticle][nbeta][d1, d2, d3 .... ]
        distance_file=open(self.folder + '/distance_Population' + repr(population+1) + '.txt',"a")
        for i in range(len(results.distances)):
            for j in range(len(results.distances[i])):
                print >>distance_file, i+1, j, results.distances[i][j], results.models[i]
        distance_file.close()

        # trajectories are stored as [nparticle][nbeta][ species ][ times ]
        traj_file=open(self.folder + '/traj_Population' + repr(population+1) + '.txt',"a")
        for i in range(len(results.trajectories)):
            for j in range(len(results.trajectories[i])): 
                arr = results.trajectories[i][j]
                nrow, ncol = numpy.shape( arr )
                for ic in range(ncol):
                    print >>traj_file, i, j, results.models[i], ic, 
                    for ir in range(nrow):
                        print >>traj_file, arr[ir,ic],
                    print >>traj_file, ""
        traj_file.close()

        if len(results.margins) > 1:
            model_file = open(self.folder + '/ModelDistribution.txt',"a")
            for m in results.margins: 
                print >>model_file, m,
            print >>model_file, ""
            model_file.close()

        nmodels = len(models)
        for mod in range(nmodels):
            #print self.folder + '/results_' + models[mod].name
            try:
                #print self.folder + '/results_' + models[mod].name
                os.chdir(self.folder + '/results_' + models[mod].name)
                #print os.getcwd()
                os.mkdir("Population_" + repr(population+1))
                os.chdir("../..")
            except:
                print "\nCan not create the folder Population_"+repr(population+1)+"!\n"
                sys.exit()

        # count number of particles in each model so that we can skip empty models
        counts = numpy.zeros([nmodels])
        nparticles = len(results.weights)
        for np in range(nparticles):
            counts[ results.models[np] ] = counts[ results.models[np] ] + 1
            

        # print out particles and weights if there are particles
        for mod in range(nmodels):
            if counts[mod] > 0:
                weight_file=open(self.folder+'/results_'+ models[mod].name + '/Population_'+repr(population+1)+'/data_Weights'+repr(population+1)+".txt","w")
                param_file = open(self.folder+'/results_'+ models[mod].name + '/Population_'+repr(population+1)+'/data_Population'+repr(population+1)+".txt","w") 
         
                nparticles = len(results.weights)
                for g in range(nparticles):
                    if( results.models[g] == mod ):
                        for k in range(len(results.parameters[g])):
                            print >>param_file, results.parameters[g][k],
                        print >>param_file, ""

                        print >>weight_file, results.weights[g]

                weight_file.close()
                param_file.close()


        # do diagnostics such as scatter plots, histograms and model distribution
        npop = len(self.all_results)
        if self.diagnostic == True:
            
            if nmodels > 1:
                # create matrix [npop][nmodel]
                m = numpy.zeros( [ len(self.all_results), nmodels ] )
                r = []
                e = []
                for i in range( len(self.all_results) ):
                    m[i,:] = self.all_results[i].margins
                    r.append( self.all_results[i].rate )
                    e.append( self.all_results[i].epsilon )
            
                getModelDistribution(m,e,r,PlotName=self.folder+'/ModelDistribution')

            # for scatter plots and histograms we require container [model][population][parameter][values]
            population_mod=[]
            weights_mod=[]
            
            for mod in range(nmodels):
                population_mod.append([])
                weights_mod.append([])

                if counts[mod] > 0:
                    PlotName = self.folder+'/results_' + models[mod].name + '/Population_'+repr(population+1) + '/ScatterPlots_Population' + repr(population+1)
                    PlotName2 = self.folder+'/results_' + models[mod].name + '/Population_'+repr(population+1) + '/weightedHistograms_Population' + repr(population+1)
                
                    for eps in range( npop ):
                        population_mod[mod].append([])
                        weights_mod[mod].append([])

                        non_const = 0
                        for param in range(models[mod].nparameters):
                            if not( models[mod].prior[param][0]==0 ):
                                population_mod[mod][eps].append([])
                                weights_mod[mod][eps].append([])

                                nparticles = len(self.all_results[eps].weights)
                                for np in range( nparticles ):
                                    if self.all_results[eps].models[np] == mod:
                                        # print mod, eps, param, self.all_results[eps].models[np], self.all_results[eps].parameters[np][param]
                                        population_mod[mod][eps][non_const].append( self.all_results[eps].parameters[np][param] )
                                        weights_mod[mod][eps][non_const].append( self.all_results[eps].weights[np] )

                                non_const = non_const + 1

                    getAllScatterPlots(population_mod, weights_mod, populations=numpy.arange(1,population+2),PlotName=PlotName,model=mod+1)
                    getAllHistograms(population_mod, weights_mod, population=population+1, PlotName=PlotName2, model=mod+1)

            if self.plotDataSeries == True:
                for mod in range(nmodels):
                    # get the first n of the accepted particles for this model
                    pars = []
                    traj2 = []
                    n = 10
                    count = 0
                    nparticles = len(results.weights)
                    for np in range(nparticles):
                        if results.models[np] == mod and count < n:
                            pars.append( results.parameters[np] )
                            traj2.append(results.trajectories[np])
                            count = count + 1

                    if len(pars) > 0:
                        filename = self.folder + '/results_' + models[mod].name + '/Population_' + repr(npop) + '/Timeseries_Population' + repr(npop)
                      #  plotTimeSeries(models[mod],pars,data,beta,filename,plotdata=self.havedata)
                        plotTimeSeries2(models[mod],pars,data,beta,filename,traj2,population,plotdata=self.havedata)


    ################ writes trajectories and parameters from simulations    
    def write_data_simulation(self, population, results, timing, models, data):

        # results abcsmc_results class
        self.all_results.append( results )

        nparticles = len(results.trajectories)
        beta = len(results.trajectories[0])

        # trajectories are stored as [nparticle][nbeta][ species ][ times ]
        traj_file=open(self.folder + '/trajectories.txt',"a")
        for i in range(nparticles):
            for j in range(beta): 
                arr = results.trajectories[i][j]
                nrow, ncol = numpy.shape( arr )
                for ic in range(ncol):
                    print >>traj_file, i, j, results.models[i], ic, 
                    for ir in range(nrow):
                        print >>traj_file, arr[ir,ic],
                    print >>traj_file, ""
        traj_file.close()

        # dump out all the parameters
        param_file=open(self.folder + '/particles.txt',"a")
        for i in range(nparticles):
            print >>param_file, i, results.models[i],
            for j in results.parameters[i]:
                print >>param_file, j,
            print >>param_file, ""
        param_file.close()

        # do timeseries plots
        npop = len(self.all_results)
        nmodels = len(models)

        # separate timeseries for each model
        for mod in range(nmodels):
            # get the first n of the accepted particles for this model
            pars = []
            traj2 = []
            n = nparticles
            count = 0
            for np in range(nparticles):
                if results.models[np] == mod and count < n:
                    pars.append( results.parameters[np] )
                    traj2.append(results.trajectories[np])

                    count = count + 1

            if len(pars) > 0:
                filename = self.folder + '/' + models[mod].name + '_timeseries'
               #  plotTimeSeries(models[mod],pars,data,beta,filename,plotdata=False)
                plotTimeSeries2(models[mod],pars,data,beta,filename,traj2,population,plotdata=False)

                        
    ################create output folders
    def create_output_folders(self, modelnames, numOutput, pickling, simulation):

        if simulation == True: pickling = False

        try:
            os.mkdir(self.folder)
            os.chdir(self.folder)
        except:
            print "\nThe folder "+ self.folder +" already exists!\n"
            sys.exit()

        if simulation == False:
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

    
