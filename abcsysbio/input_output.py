import os, sys, pickle, string
import numpy

from getResults import getAllScatterPlots
from getResults import getAllHistograms
from getResults import plotTimeSeries
from getResults import getModelDistribution
from getResults import plotData

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


class input_output:
    
    def __init__(self, folder, restart, diagnostic, plotDataSeries, havedata=True):
        self.folder = folder
        self.diagnostic = diagnostic
        self.plotDataSeries = plotDataSeries
        self.havedata = havedata
      
        # Hold all data here for plotting purposes.
        # May want to remove this as could get large
        # self.all_results = []
        
        if restart==True:
            self.folder = self.folder + '_restart'

    def plot_data(self, data):
        if self.havedata == True: plotData( data, self.folder+'/_data' )


    ################write rates, distances, trajectories    
    def write_data(self, population, results, timing, models, data):

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
                link_file = open(self.folder+'/results_'+ models[mod].name + '/Population_'+repr(population+1)+'/data_Links'+repr(population+1)+".txt","w") 
                link_sum_file = open(self.folder+'/results_'+ models[mod].name + '/Population_'+repr(population+1)+'/data_Links_sum'+repr(population+1)+".txt","w")

                # count the number of links in this model
                nl = 0
                npar = models[ mod ].nparameters
                for i in range(npar):
                    if models[ mod ].prior[i][0] == 4:
                        nl += 1
                ##print "number of links in this model", nl
                counter = {}
         
                nparticles = len(results.weights)
                for g in range(nparticles):
                    if( results.models[g] == mod ):

                        # print parameter
                        for k in range(len(results.parameters[g])):
                            print >>param_file, results.parameters[g][k],
                        print >>param_file, ""

                        # print weights
                        print >>weight_file, results.weights[g]

                        # print weights and links
                        print >>link_file, results.weights[g], 
                        for k in range(len(results.parameters[g])):
                            if models[mod].prior[k][0] == 4:
                                print >>link_file, results.parameters[g][k],
                        print >>link_file, ""


                        # place links into dictionary
                        rep = ""
                        for k in range(len(results.parameters[g])):
                            if models[ mod ].prior[k][0] == 4:
                                # add one to representation so they are positive
                                rep += repr( results.parameters[g][k]+1 )

                        ## print rep

                        if not rep in counter:
                            counter[rep] = results.weights[g]
                        else:
                            counter[rep] += results.weights[g]

                # print links summary
                ##print counter
                for keys in sorted(counter, key=counter.get, reverse=True):
                    m = [int(keys[i])-1 for i in range(nl)]
                    for j in range(nl):
                        print >>link_sum_file, m[j],
                    print >>link_sum_file, sum(numpy.abs(m)), counter[keys]   

                weight_file.close()
                param_file.close()
                link_file.close()
                link_sum_file.close()
                

        # do diagnostics such as scatter plots, histograms and model distribution
        if self.diagnostic == True:
            
            #if nmodels > 1:
            #    # create matrix [npop][nmodel]
            #    m = numpy.zeros( [ npop, nmodels ] )
            #    r = []
            #    e = []
            #    for i in range( npop ):
            #        m[i,:] = self.all_results[i].margins
            #        r.append( self.all_results[i].rate )
            #        e.append( self.all_results[i].epsilon )
            
            #    getModelDistribution(m,e,r,PlotName=self.folder+'/ModelDistribution')

            # for scatter plots and histograms we require container [model][population][parameter][values]
            population_mod=[]
            weights_mod=[]
            
            for mod in range(nmodels):
                population_mod.append([])
                weights_mod.append([])

                if counts[mod] > 0:
                    #PlotName = self.folder+'/results_' + models[mod].name + '/Population_'+repr(population+1) + '/ScatterPlots_Population' + repr(population+1)
                    PlotName2 = self.folder+'/results_' + models[mod].name + '/Population_'+repr(population+1) + '/weightedHistograms_Population' + repr(population+1)

                    population_mod[mod].append([])
                    weights_mod[mod].append([])

                    non_const = 0
                    for param in range(models[mod].nparameters):
                        if not( models[mod].prior[param][0]==0 ):
                            population_mod[mod][0].append([])
                            weights_mod[mod][0].append([])

                            nparticles = len(results.weights)
                            for np in range( nparticles ):
                                if results.models[np] == mod:
                                    ##print mod, param, results.models[np], results.parameters[np][param]
                                    population_mod[mod][0][non_const].append( results.parameters[np][param] )
                                    weights_mod[mod][0][non_const].append( results.weights[np] )

                            non_const = non_const + 1

                    #getAllScatterPlots(population_mod, weights_mod, populations=numpy.arange(1,population+2),PlotName=PlotName,model=mod+1)
                    getAllHistograms(population_mod, weights_mod, population=1, PlotName=PlotName2, model=mod+1)

            if self.plotDataSeries == True:
                for mod in range(nmodels):
                    # get the first n of the accepted particles for this model
                    pars = []
                    n = 10
                    count = 0
                    nparticles = len(results.weights)
                    for np in range(nparticles):
                        if results.models[np] == mod and count < n:
                            pars.append( results.parameters[np] )
                            count = count + 1

                    if len(pars) > 0:
                        filename = self.folder + '/results_' + models[mod].name + '/Population_' + repr(population+1) + '/Timeseries_Population' + repr(population+1)
                        #plotTimeSeries(models[mod],pars,data,beta,filename,plotdata=self.havedata)
                        
                        pp = PdfPages( filename + ".pdf" )

                        # trajectories are stored as list [nparticle][nbeta][ species ][ times ] not numpy array

                        # here we assume beta=1
                        for i in range(len(results.trajectories)):
                            # print "printing traj", i
                            if( results.models[i] == mod):
                            #if i < 500:
                                arr = results.trajectories[i][0]
                                nrow, ncol = numpy.shape( arr )
                                # print nrow, ncol
                                for ic in range(ncol):
                                    plt.plot(arr[:,ic], label='sp '+repr(ic))

                                #dists = [ round( results.distances[i][0][k],3) for k in range(len( results.distances[i][0])) ]
                                #labs = [ repr(x) for x in dists ]
                                #plt.text(10, 0.8*max(arr), string.join(labs, ", ") )

                                plt.title("particle " +repr(i) )
                                legend = plt.legend(loc='upper left', shadow=False)
                                # Set the fontsize
                                for label in legend.get_texts():
                                    label.set_fontsize('small')

                                for label in legend.get_lines():
                                    label.set_linewidth(0.5)

                                pp.savefig()
                                plt.close()

                        pp.close()




    ################ writes trajectories and parameters from simulations    
    def write_data_simulation(self, population, results, timing, models, data):

        # results abcsmc_results class
        # self.all_results.append( results )

        nparticles = len(results.trajectories)
        beta = len(results.trajectories[0])

        # trajectories are stored as [nparticle][nbeta][ species ][ times ]
        traj_file=open(self.folder + '/trajectories.txt',"a")
        for i in range(nparticles):
            ##print results.trajectories[i][0]
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

        # dump out the links
        link_file = open(self.folder+'/links.txt',"a")
        for i in range(nparticles):
            mod = results.models[i]
            print >>link_file, i, results.models[i],
            for j in range(len(results.parameters[i])):
                if models[ mod ].prior[j][0] == 4:
                    print >>link_file, results.parameters[i][j],  
            print >>link_file, ""
        link_file.close()

        # get links summary but assuming only one model
        # count the number of links in this model
        nl = 0
        npar = models[ 0 ].nparameters
        for i in range(npar):
            if models[ 0 ].prior[i][0] == 4:
                nl += 1

        ## print "number of links in this model", nl
        
        weights = [1/float(nparticles) for i in range(nparticles)]
        counter = {}
        for i in range(nparticles):
            mod = results.models[i]
            rep = ""
            for j in range(len(results.parameters[i])):
                if models[ mod ].prior[j][0] == 4:
                    # add one to representation so they are positive
                    rep += repr( results.parameters[i][j]+1 )

            ## print rep
            if not rep in counter:
                counter[rep] = weights[i]
            else:
                counter[rep] += weights[i]

        # print links summary
        ## print counter
        link_sum_file = open(self.folder+'/link_summary.txt',"a")
        for keys in sorted(counter, key=counter.get, reverse=True):
            m = [int(keys[i])-1 for i in range(nl)]
            for j in range(nl):
                print >>link_sum_file, m[j],
            print >>link_sum_file, sum(numpy.abs(m)), counter[keys]   
        link_sum_file.close()
        
        # do timeseries plots
        npop = 1
        nmodels = len(models)

        # separate timeseries for each model
        for mod in range(nmodels):
            # get the first n of the accepted particles for this model
            pars = []
            n = nparticles
            count = 0
            for np in range(nparticles):
                if results.models[np] == mod and count < n:
                    pars.append( results.parameters[np] )
                    count = count + 1

            if len(pars) > 0:
                filename = self.folder + '/' + models[mod].name + '_timeseries'
                plotTimeSeries(models[mod],pars,data,beta,filename,plotdata=False)
                
                        
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

    
