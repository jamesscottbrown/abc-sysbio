# standard imports
import re, os, sys, pickle, time, math

import numpy
from numpy import random as rnd
from scipy.stats import norm

from abcodesolve        import abcodeint
from sdeint             import sdeint
from GillespieAlgorithm import GillespieInt
from realize            import realize_particles

from getResults import MatrixToTextfile
from getResults import getAllScatterPlots
from getResults import getAllHistograms
from getResults import plotTimeSeries
from getResults import printModelDistribution
from getResults import getModelDistribution
from evaluate import evaluateDistance

import abcsysbio

__LARGE_DISTANCE__ = 1000000
__XMAX__= 2000

def w_choice(item,weight):#http://snippets.dzone.com/posts/show/732 - modified for vectors
    n = rnd.random_sample()
    for i in range(0,len(weight)):
        if n < weight[i]:
            break
        n = n - weight[i]
    return i


def w_choice2(item,weight):#http://snippets.dzone.com/posts/show/732 - modified for vectors
    n = rnd.random_sample()
    for i in range(0,len(weight)):
        if n < weight[i]:
            break
        n = n - weight[i]
    return item[i]

# calculates the index of the model, given a result index
def _getModelNumber(index, numberOfModels):
    sumIndex = 0
    for i in range(len(numberOfModels)):
        if(sumIndex > index):
            return i-1
        sumIndex += numberOfModels[i]
        
    return len(numberOfModels)-1

################## compute the kernels
def getKernel(kern,population):

    for param in range(0, len(population)):
        minimum=min(population[param])
        maximum=max(population[param])
        scale=(maximum-minimum)
        if kern[param][0]==1: 
            kern[param][1]=-scale/2.0
            kern[param][2]=scale/2.0
        if kern[param][0]==2:
            kern[param][2]=scale/2.0

    return kern


##################compute the pdf of uniform distribution
def getPdfUniform(scale1,scale2,parameter):
    if ((parameter>scale2) or (parameter<scale1)):
        return 0.0
    else:
        return 1/(scale2-scale1)


##################compute the pdf of gauss distribution
def getPdfGauss(mean,scale,parameter):
    x=((math.sqrt(2*math.pi)*scale)**(-1))*math.exp(-1.0*((parameter-mean)**2)/(2.0*scale**2))

    return x


################compute the pdf of lognormal distribution
def getPdfLognormal(mean,sigm,parameter):
    x=((sqrt(2*math.pi)*sigm*parameter)**(-1))*exp(-0.5*(((math.ln(parameter)-sigm)**2)/sigm)**2)
    
    return x


###########simulate data with sampled parameter set for the sampled model
def simulateData(ModelName,selected_model,integrationType,InitValues,timepoints,sampleParameters,dt,rtol=None,atol=None):

    o=re.compile('ODE')
    s=re.compile('SDE')
    g=re.compile('Gillespie')
  
    module = __import__(ModelName[selected_model])
    print ModelName[selected_model], module

    if o.match(integrationType[selected_model]):
        samplePoints=abcodeint(module, InitValues[selected_model], timepoints, sampleParameters, dt, rtol, atol)
    if s.match(integrationType[selected_model]): 
        samplePoints=sdeint(module,InitValues[selected_model],sampleParameters,timepoints,dt)
    if g.match(integrationType[selected_model]): 
        samplePoints=GillespieInt(module, InitValues[selected_model], sampleParameters,timepoints)

    return samplePoints

###########how to fit variable to the given data (if fit is None, data for all variables are available in order of the model)
def howToFitData(fitting_instruction,samplePoints):

    if fitting_instruction!=None:
        points=numpy.zeros([len(samplePoints), len(fitting_instruction)])
        for i in range(0, len(fitting_instruction)):
            points[:,i] = eval(fitting_instruction[i])
           
    else:
        points=samplePoints


    return points

def sample_particle(nparticle, selected_model, margins_prev, model_prev, weights_prev ):
    u = rnd.uniform(low=0, high=margins_prev[selected_model])
    F = 0

    for i in range(0,nparticle):
        if int(model_prev[i]) == int(selected_model) :
            F = F + weights_prev[i]

            if(F > u):
                break

    return i

def perturbParticle(params, priors, kernel):
    np = len( priors )
    prior_prob = 1
            
    for n in range(0,np):
         
        if kernel[n][0]==1:
            params[n] = params[n] + rnd.uniform(low=kernel[n][1],high=kernel[n][2])

        if kernel[n][0]==2: 
            params[n] = params[n] + rnd.normal(kernel[n][1],kernel[n][2])

        x = 1.0
        if priors[n][0]==2: 
            x=getPdfUniform(priors[n][1],priors[n][2],params[n])

        if priors[n][0]==3: 
            x=getPdfLognormal(priors[n][1],priors[n][2],params[n])
            
        prior_prob = prior_prob*x

    return prior_prob
    
############sample parameters
def sampleTheParameter(t, nparticle, priors, kernel, selected_model, margins_prev, model_prev, weights_prev, parameters_prev, restart, perturbfn ):
    np = len( priors )
    sampleParameters = [ 0 for i in range(0,np) ]
    
    if(restart==False and t==0):
        
        # sample from prior
        for n in range(0,np):
            if priors[n][0] == 0: sampleParameters[n]=priors[n][1]
            #if priors[n][0] == 1: sampleParameters[n]= sample from Jeffrey's prior
            if priors[n][0] == 2: sampleParameters[n]=rnd.uniform(low=priors[n][1],high=priors[n][2])
            if priors[n][0] == 3: sampleParameters[n]=rnd.lognormal(priors[n][1],priors[n][2])
        
    else:

        prior_prob = -1
        while prior_prob <= 0 :

            # sample putative particle from previous population
            p = sample_particle(nparticle, selected_model, margins_prev, model_prev, weights_prev )
            for nn in range(np):
                sampleParameters[nn] = parameters_prev[ p ][nn]
          
            prior_prob = perturbfn( sampleParameters, priors, kernel )

    return sampleParameters

############sample model
def sampleTheModel(t, nmodel, margins_prev, modelDistrib0, dead_models, restart, modelKernel ):
    
    if nmodel == 1:
        return 0

    if(restart == False and t == 0):
        # sample from prior
        selected_model=w_choice( range(0,nmodel), modelDistrib0 )
        return selected_model

    else:
        sampled_model=w_choice( range(0,nmodel), margins_prev )

        if( len(dead_models) == nmodel-1 ):
            return sampled_model
        
        # perturb model
        u = rnd.uniform(low=0, high=1)
        if u > modelKernel:

            # sample randomly from other (non dead) models
            not_available = dead_models[:]
            not_available.append(sampled_model)
            ss = set( not_available )
            s = set( range(0,nmodel) )

            #print "perturbing model ", set( range(0,nmodel) ), ss, s-ss, list(s-ss), numpy.array( list(s-ss) )

            ar = numpy.array( list(s-ss) )
            rnd.shuffle( ar )
            perturbed_model = ar[0]
            #print "perturbation ", sampled_model, "->", perturbed_model
            return perturbed_model
        else:
            #print "perturbation ", sampled_model, "->", sampled_model
            return sampled_model         

def getPdfModelKernel(m,m0,nmodel,dead_models):
    ndead = len(dead_models)
   
    if(ndead == nmodel-1):
        return 1.0
    else:

        if(m == m0):
            return 0.7
        else:
            return 0.3/(nmodel-ndead)

def getPdfParameterKernel(params, params0, priors, kernel):

    prob = 1
    for n in range(0,len(priors)):
        kern = 0

        if not(priors[n][0]==0):
                       
            if kernel[n][0]==1: 
                scale = (kernel[n][2]-kernel[n][1])/2.0
                scale1 = max( (params0[n] - scale), priors[n][1])
                scale2 = min( (params0[n] + scale), priors[n][2])
                kern = getPdfUniform(scale1,scale2, params[n])

            if kernel[n][0]==2:
                mean = params0[n]
                scale = kernel[n][2]
                CDF2 = norm.cdf(priors[n][2],mean,scale)
                CDF1 = norm.cdf(priors[n][1],mean,scale)
                CDF = (CDF2-CDF1)**(-1)
                kern = getPdfGauss(mean,scale,params[n])
                kern = kern*CDF

        else: kern=1.0
        prob=prob*kern

    return prob


def computeParticleWeights(t, sampleParameters, selected_model, nmodel, numOutput, priors, kernel, model_prev,
                           margins_prev, weights_prev, parameters_prev, BETA, dead_models, restart, kernelpdffn):

    if(restart==False and t == 0):
        return BETA
    else:
        # calculate model prior probility = 1/nmodel
        mprob = 1/float(nmodel)

        # particle prior probability
        pprob = 1
        for n in range(0,len(priors)):
            x = 1.0
            if priors[n][0]==0:
                x=1
            if priors[n][0]==2: 
                x=getPdfUniform(priors[n][1],priors[n][2],sampleParameters[n])
            if priors[n][0]==3: 
                x=getPdfLognormal(priors[n][1],priors[n][2],sampleParameters[n])

            pprob = pprob*x

        numer = BETA * mprob * pprob
        
        denom_m = 0
        for i in range(nmodel):
            denom_m = denom_m + margins_prev[i]*getPdfModelKernel(selected_model, i, nmodel, dead_models)
        denom = 0
        #print "Calculating denom\t", selected_model, sampleParameters
        for j in range(numOutput):

            if(int(selected_model) == int(model_prev[j]) ):
                #print "\t", j, model_prev[j], weights_prev[j], parameters_prev[j]
                denom = denom + weights_prev[j] * kernelpdffn(sampleParameters, parameters_prev[j], priors, kernel )

        #print "numer/denom_m/denom/m(t-1) : ", numer,denom_m, denom, margins_prev[selected_model]

        return numer/(denom_m*denom/margins_prev[selected_model])

def normalizeWeights(numOutput, weights_curr):
    n = sum( weights_curr )
    for i in range(numOutput):
        weights_curr[i] = weights_curr[i]/float(n)

    return weights_curr

def modelMarginals(nmodel, numOutput, model_curr, weights_curr):
    ret  = [0 for i in range(nmodel)]

    for i in range(nmodel):
        for j in range(numOutput):
            if int(model_curr[j]) == int(i):
                ret[i] = ret[i] + weights_curr[j]
            
    return ret

################create output folders
def create_output_folders(fname, Nmodel, ModelName, numOutput, restart, pickling):
    if restart==True:
        folder=fname + '_restart'
    else:
        folder=fname
    try:
        os.mkdir(folder)
        os.chdir(folder)
    except:
        print "\nThe folder "+ folder +" already exists!\n"
        sys.exit()
    for mod in range(0,len(ModelName)):
        try:
            os.mkdir('results_'+ModelName[mod])
        except:
            print "\nThe folder "+ folder + "/results_" + ModelName[mod]+" already exists!\n"
            sys.exit()

    os.chdir('..')

            
    #rate_file=open(folder+'/rates.txt',"w")

        
    if pickling==True:
        try:
            os.chdir(folder)
            os.mkdir('copy')
            os.chdir('..')
        except:
            print "\nThe folder \'copy\' already exists!\n"
            sys.exit()
    
        out_file=open(folder+'/copy/algorithm_parameter.dat',"w")
        pickle.dump(numOutput,out_file)
        out_file.close()

    return folder

################read the stored data
def read_pickled(fname):
    # pickle numbers selected model of previous population
    # pickle population of selected model of previous population pop_pickled[selected_model][n][vectpos]
    # pickle weights of selected model of previous population weights_pickled[selected_model][n][vectpos]

    try:
        in_file=open(fname + '/copy/model_last.dat',"r")
        model_pickled=pickle.load(in_file)
        in_file.close()
    except:
        print "\nCan not find file \'model_last.dat\' in folder \'copy\'!\n"
        sys.exit()

    try:
        in_file=open(fname + '/copy/weights_last.dat',"r")
        weights_pickled=pickle.load(in_file)
        in_file.close()
    except:
        print "\nCan not find file \'weights_last.dat\' in folder \'copy\'!\n"
        sys.exit()

    try:
        in_file=open(fname + '/copy/params_last.dat',"r")
        parameters_pickled=pickle.load(in_file)
        in_file.close()
    except:
        print "\nCan not find file \'params_last.dat\' in folder \'copy\'!\n"
        sys.exit()

    try:
        in_file=open(fname + '/copy/margins_last.dat',"r")
        margins_pickled=pickle.load(in_file)
        in_file.close()
    except:
        print "\nCan not find file \'margins_last.dat\' in folder \'copy\'!\n"
        sys.exit()
    
    try:
        in_file=open(fname + '/copy/kernels_last.dat',"r")
        kernel=pickle.load(in_file)
        in_file.close()
    except:
        print "\nCan not find file \'kernels_last.dat\' in folder \'copy\'!\n"
        sys.exit()

    print "\n\n\n Reading previous population"

    print "model_pickled", model_pickled, "\n\n\n"
    print "weights_pickled", weights_pickled, "\n\n\n"
    print "parameters_pickled", parameters_pickled, "\n\n\n"
    print "margins_pickled", margins_pickled, "\n\n\n"

    return [model_pickled, weights_pickled, parameters_pickled, margins_pickled, kernel]

################write the stored data
def write_pickled(folder, nmodel, numOutput, priors, model_prev, weights_prev, parameters_prev, margins_prev, kernel):

    out_file=open(folder+'/copy/model_last.dat',"w")
    x = model_prev[:]
    pickle.dump(x,out_file)
    out_file.close()

    out_file=open(folder+'/copy/weights_last.dat',"w")
    x = weights_prev[:]
    pickle.dump(x,out_file)
    out_file.close()

    out_file=open(folder+'/copy/params_last.dat',"w")
    x= parameters_prev
    pickle.dump(x,out_file)
    out_file.close()
    
    out_file=open(folder+'/copy/margins_last.dat',"w")
    x = margins_prev[:]
    pickle.dump(x,out_file)
    out_file.close()
 
    out_file=open(folder+'/copy/kernels_last.dat',"w")
    x=[]
    for mod in range(0,nmodel):
        x.append(kernel[mod])
    pickle.dump(x,out_file)
    out_file.close()

################make some nice plots
def make_analysis_plots(folder, diagnostic, plotDataSeries, timepoints, fit, dt, data, t, nmodel, ModelName, integrationType,
                        InitValues, priors, population, numbers, modelDistribution, weights, epsilon, rate ):
    
    for mod in range(0,nmodel):
        try:
            os.chdir(folder+'/results_'+ModelName[mod])
            os.mkdir("Population_"+repr(t+1))
            os.chdir("../..")
        except:
            print "\nCan not create the folder Population_"+repr(t+1)+"!\n"
            sys.exit()


    for mod in range(0,nmodel):
        weight_file=open(folder+'/results_'+ModelName[mod]+'/Population_'+repr(t+1)+'/data_Weights'+repr(t+1)+".txt","w")
        for g in range(0,len(weights[mod][t])):
            weight_file.write(repr(weights[mod][t][g])+"\n")
        weight_file.close()

            
    population_mod=[]
    weights_mod=[]
    for mod in range(0,nmodel):
        PlotName=folder+'/results_'+ModelName[mod]+'/Population_'+repr(t+1)+'/ScatterPlots_Population'+repr(t+1)
        filename1=folder+'/results_'+ModelName[mod]+'/Population_'+repr(t+1)+'/data_Population'+repr(t+1)
        filename2=folder+'/results_'+ModelName[mod]+'/Population_'+repr(t+1)+'/data_Weights'+repr(t+1)
        PlotName2=folder+'/results_'+ModelName[mod]+'/Population_'+repr(t+1)+'/weightedHistograms_Population'+repr(t+1)
           
        MatrixToTextfile(population,filename=filename1,model=mod,eps=t)
        if diagnostic==True:
                
            population_mod.append([])
            weights_mod.append([])
            for eps in range(0,len(population[0])):
                population_mod[mod].append([])
                weights_mod[mod].append([])
                for param in range(0,len(priors[mod])):
                    if not(priors[mod][param][0]==0):
                        population_mod[mod][eps].append(population[mod][eps][param])
                        weights_mod[mod][eps].append(weights[mod][eps])

            getAllScatterPlots(population_mod,weights_mod,populations=numpy.arange(1,t+2,1),PlotName=PlotName,model=mod+1)
            getAllHistograms(population_mod,weights_mod,population=t+1,PlotName=PlotName2, model=mod+1)

    if nmodel>1:
        printModelDistribution(modelDistribution,eps=t,filename=folder+'/ModelDistribution.txt')
        if diagnostic==True:
            getModelDistribution(modelDistribution,epsilon,rate,PlotName=folder+'/ModelDistribution')


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




#######################################################################
# The ABC SMC algorithm for model selection                           #
#######################################################################

def abcSimulator(ModelName,data,timepoints,numOutput, epsilon, InitValues,integrationType,modelWeights,priors,kernel,fit=None,
                 sampleFromPrior=False,beta=1,dt=0.01,restart=False,source=None,diagnostic=True,pickling=True,plotDataSeries=True,
                 full=False, rtol=None, atol=None, constKernels=False,
                 runmode=0, timing=False, nthread=None, nblock=None, fname="_results_", mt_data="MersenneTwister.dat",app_file="app.txt",
                 custom_kernel=False, custom_distance=False, modelKernel=0.7):   
    #ToDo: Allow for different initialization values?!

    """

    The simulator for the abc algorithm:

    ***** args *****

    ModelName:
              a list of strings, where each string represents a name for a model. 
              The length of this list is determind by the number of models that 
              are investigated.
    
              Example: 
              To investigate 3 models with the names "model1", "model2" and "model3": 
              ModelName = ('model1','model2','model3')

              To investigate one model with the name "MyInterestingModel":
              ModelName = ('MyInterestingModel',)

    data:
              a list of floats that coresponds to the experimental data. 
              The length of this list is determind by the length of the argument timepoints.
              
    timepoints:
              a list of floats.
              The length must not differ to the length of data. Each timepoint describes an
              external timepoint, where data are available. 

    numOutput:
              an integer number.
              numOutput corresponds to the population size, i.e. how many particles must be 
              accepted in each population.

    epsilon: 
              a 2D list of floats.
              Epsilon corresponds to the maximum distance between experimental and simmulated 
              data. The length of this list corresponds to the number of populations. If the 
              second dimension > 1, the user needs to define its own distance function!

              Example:
              maximum distance for ABC Rejection (only one population):
              epsilon = [[1.5,]]
              maximum distance for ABC SMC with five populations:
              epsilon = [[3.5,3.0,2.5,2.0,1.5]]
              maximum distance for ABC SMC with five populations, but more complex distance function:
              epsilon = [[3.5,3.0,2.5,2.0,1.5],[35,30,25,20,15],[1.0,0.8,0.6,0.4,0.2]]

    InitValues: 
              a 2D list of floats.
              The first dimension describes each model, the second dimension describes the 
              initial condition for each species in the particular model. 
              
              Example:
              2 models, the first model has 3 species, the second model has 4 species:
              InitValues = [ [1.9 , 20 , 32.1] , [2.0 , 1.3 , 2.4 , 3.4] ]
              1 model with only one species:
              InitValues = [ [2.0] ,]

    IntegrationType:
              a list of predefined strings.
              The length of this list is determind by the number of models to be investigated.

              Possible strings are:
              'ODE'       --- use the odeint (Scipy) to simulate a deterministic system
              'SDE'       --- use the sdeint (abc) to simulate the system with SDEs
              'Gillespie' --- use the GillespieAlgorithm (abc) to simulate a stochastic system

              Example:
              1 model to be simulated with SDEs
              IntegrationType = ('SDE' ,)
              3 models, first two deterministic, third stochastic via Gillespie
              IntegrationType = ('ODE' , 'ODE' , 'Gillespie')

    modelWeights:
              a list of floats.
              This list describes the discrete prior distribution of the models.The length of 
              the list is determined by the number of models to be investigated. The entries 
              of the list display the weights for each model, starting from a uniform 
              distribution. These weights do not need to be normalized to 1.

              Example:
              1 model
              modelWeights = (1,)
              5 models with purely uniform distribution
              modelWeights = (1,1,1,1,1)
              5 different weighted models
              modelWeights = (2,1,4,3,7)

    priors: 
              a 3D list.
              The first dimension represents the number of models, the second dimension
              represents the number of parameter for a particular model and the third dimension
              represents the distribution for this parameter and has a constant length of 3.
              The first entry for each parameter is an integer number that stands for a
              specific distribution.
              
              Implemented distributions are:
              0   ---   constant parameter.
                        Example: constant parameter with value 12.3
                        [0 , 12.3 , x] , where x can be any number

              1   ---   uniform distribution. 
                        Example: uniform distribution in the range 0.1 to 50
                        [1 , 0.1 , 50]

              2   ---   lognormal distribution.
                        Example: lognormal distribution with mean 3 and varianz 1.5
                        [2 , 3 , 1.5]

              Example:
              1 model with 3 parameter, the first two parameter have uniform prior between 0 and
              5, the third parameter has lognormal prior with mean 1 and varianz 3.
              [ [ [1,0,5],[1,0,5],[2,1,3] ], ]
              2 models where the first model has 2 parameter (the first is constant 3 and the 
              second is uniform between 0 and 1) and the second model has 1 lognormal parameter
              with 0 mean and varianz 0.5
              [ [ [0,3,0],[1,0,1] ] , [ [2,0,0.5],] ]

    
    kernel:
              a 3D list.
              This list has the same structure as the 'priors', but it has different 
              distributions implemented.

              Implemented distributions are:
              0   ---   uniform distribution. 
                        Example: uniform distribution in the range -2 to 2
                        [0 , -2 , 2]
              
              1   ---   gaussian distribution.
                        Example : gaussian distribution with mean 0 and varianz 0.5
                        [1 , 0 , 0.5]


    ***** kwargs *****
    
    
    fit:      
              a 2D list of strings.
              This list contains the fitting instructions and therefore defines how to fit the
              experimental data to the systems variables. The first dimension represents the 
              models to be investigated, the second dimension represents the amount of data.
              
              Example:
              1 model with 7 species, 3 data series
              fit = [ ['species1+species2' , 'species5' , '3*species4-species7'], ]
              2 models with each 7 species, 2 data series
              fit = [ ['species1' , 'species2*species3'] , ['species1' , 'species2*species4'] ]
              
              The default value is 'None', i.e. the order of data series corresponds exactly 
              to the order of species.

    beta:
              an integer number.
              It represents the number how often a system is simulated with a sampled parameter
              set. For deterministic systems beta is 1, for stochastic systems beta should be 
              larger than 1.

    dt:
              a float number.
              dt is the internal time step for SDE systems. The ODE solver also uses dt as an 
              additional internal time step.

    restart:
              a boolean variable.
              Start the abc simulator from a given population.
              The default value is 'False'.

    source:   a list of strings.
              This length of the list is equal to the number of models to be investigated.
              The entries describe the SBML sources for each model. 
              The default value is 'None'.


    diagnostic:
              a boolean variable.
              Create diagnostic plots after each population.
              Diagnostic plots are: 
              - model distribution (if more than 1 model to analize)
              - all combinations of scatter plots for parameters to be inferred
              The default value is 'False'.

    pickling:
              a boolean variable.
              Make a backup after each population that allows to start the abc simulator 
              from a given population (see also 'restart').
              The dafault value is 'False'.

    plotDataSeries:
              a boolean variable.
              Compute and plot after each population the timeseries for the first 10 sampled
              parameters inclusive data points.
              The dafault value is 'True'.            

    full:
              a boolean variable.
              Print status while running the abc simulator.
              The default value is 'False'.

    constKernels:
              a boolean variable.
              Do not correct the pertubation kernels in respect to the previous
              population.
              The default value is 'False'.

    """

    #
    # Define all variable and containers here
    #
    npop = len(epsilon[0])
    nmodel = len(ModelName)

    # model prior distribution is uniform over number of models
    modelDistrib0=numpy.zeros([nmodel])
    for mod in range(nmodel):
        modelDistrib0[mod] = 1/float(nmodel)

    # modelDistribution contains the marginal model distribution indexed via [pop][model]
    modelDistribution=numpy.zeros([npop,nmodel])

    # numbers contains the marginal model distribution indexed via [model][pop]
    numbers=numpy.zeros([nmodel, npop])
     
    # population contains the parameters indexed via [model][pop][parameter]
    # weights contains the weights indexed via [model][pop]
    population=[]
    weights=[]

    for j in range(0,nmodel):
        population.append([])
        weights.append([])
        for eps in range(0,npop):
            population[j].append([])
            weights[j].append([])
            for param in range(0,len(priors[j])):
                population[j][eps].append([])                              

    # hits, sampling, rate are arrays indexed by population 
    hits=[]
    sampling=[]
    rate=[]
    
    for j in range(0,npop):
        hits.append(0.0)
        sampling.append(0.0)
        rate.append(0.0)

    #
    # containers for joint space (m, theta) inference
    #
    model_prev      = [0  for i in range(0,numOutput)]
    weights_prev    = [0  for i in range(0,numOutput)]
    parameters_prev = [[] for i in range(0,numOutput)]
    margins_prev    = [0  for i in range(0,nmodel)] 

    model_curr      = [0  for i in range(0,numOutput)]
    weights_curr    = [0  for i in range(0,numOutput)]
    parameters_curr = [[] for i in range(0,numOutput)]
    margins_curr    = [0  for i in range(0,nmodel)] 

    dead_models = []

    
    # Create results folder for each model and create the pickling folder
    folder = create_output_folders(fname, nmodel, ModelName, numOutput, restart, pickling) 

    # Load the pickled data
    if restart==True:
        model_pickled, weights_pickled, parameters_pickled, margins_pickled, kernel =  read_pickled( fname )
        
        # Load the data into the previous population containers
        for p in range(0,numOutput):
            model_prev[p] = model_pickled[p]
            weights_prev[p] = weights_pickled[p]

            parameters_prev[p] = parameters_pickled[p][ : ]

        for m in range(0,nmodel):
            margins_prev[m] = margins_pickled[m]

        # Check for dead models
        dead_models = []
        for j in range(nmodel):
            if margins_prev[j] < 1e-6:
                dead_models.append(j)

    #
    # Read, compile cuda code. Initialize RNGs 
    #
    if runmode == 2:
        
        o=re.compile('ODE')
        s=re.compile('SDE')
        g=re.compile('Gillespie')

        # First check that all models are either ode, or mjp
        int_type = 0
        if s.match(integrationType[0]):
            print "\n cuda SDE integration not implemented \n"
            sys.exit()
        elif g.match(integrationType[0]):
            int_type = 3
        elif o.match(integrationType[0]):
            int_type = 1

        for j in range(1,nmodel):
            if integrationType[0] != integrationType[j]:
                print "\n cuda supports same integration type only \n"
                sys.exit()
                
        model_np = []
        model_nx = []

        # pmax, xmax
        for j in range(0,nmodel):
            model_np.append( len(priors[j]) )
            model_nx.append( len(InitValues[j]) )

        pmax = max(model_np)
        xmax = max(model_nx)

        import abccuda

        if int_type == 1:
            import abccuda.lsoda

            # get the location of the cuLsoda_all.cu file
            s = abccuda.__file__
            ls_cu = s.replace('__init__.pyc','cuLsoda_all.cu')
            print '#### using cuLsoda_all.cu at', ls_cu
            
            ptrs = abccuda.lsoda.compile_lsoda( app_file, ls_cu, nthread, nblock, pmax, xmax, options=[], write=False )
        
        if int_type == 3:
            import abccuda.gillespie
            
            # get the location of the MersenneTwister.cu file
            s = abccuda.__file__
            mt_cu = s.replace('__init__.pyc','MersenneTwister.cu')
            print '#### using MersenneTwister.cu at', mt_cu

            # get the location of the MersenneTwister.dat file
            s = abccuda.__file__
            mt_data = s.replace('__init__.pyc','MersenneTwister.dat')
            print '#### using MersenneTwister.dat at', mt_data

            ptrs = abccuda.gillespie.compile_gillespie ( app_file, mt_cu, mt_data, nthread, nblock, pmax, xmax, options=[], write=False )
    
    #
    # Using cuda-sim
    # Read, compile cuda code. Initialize RNGs 
    #
    if runmode == 3:
        
        o=re.compile('ODE')
        s=re.compile('SDE')
        g=re.compile('Gillespie')
        
        import cudasim
        cudaModel = []
        
        # initialize and compile all CUDA kernels; assume that CUDA kernels
        # given in appfile= are in the same order as in the input file
        # assume that CUDA kernel-files are separated with commas
        cudaKernels = app_file.split(",")
        print cudaKernels
        
        for i in range(nmodel):
        
            # check integration type
            if s.match(integrationType[i]):
                int_type = 2
            elif g.match(integrationType[i]):
                int_type = 3
            elif o.match(integrationType[i]):
                int_type = 1
            
            if int_type == 1:
                import cudasim.Lsoda as Lsoda
                codeFile = cudaKernels[i]
                cudaModel.append(Lsoda.Lsoda(timepoints, codeFile, dt=dt))
            elif int_type == 2:
                import cudasim.EulerMaruyama as EulerMaruyama
                codeFile = cudaKernels[i]
                cudaModel.append(EulerMaruyama.EulerMaruyama(timepoints, codeFile, dt=dt, beta=beta))
            else:
                import cudasim.Gillespie as Gillespie
                codeFile = cudaKernels[i]
                cudaModel.append(Gillespie.Gillespie(timepoints, codeFile, dt=dt, beta=beta))
        
    #
    # Set distance function
    #
    if custom_distance == False:
        from abcsysbio import euclidian
        distancefn = euclidian.euclidianDistance
    else:
        import customABC
        distancefn = customABC.distance

    #
    # Set kernel function
    #
    if custom_kernel == False:
        kernelfn = getKernel
        kernelpdffn = getPdfParameterKernel
        perturbfn = perturbParticle

    else:
        import customABC
        kernelfn = customABC.getKernel
        kernelpdffn = customABC.getPdfParameterKernel
        perturbfn = customABC.perturbParticle
    #
    # Start of ABC-SMC algoritm
    #
    if timing == True:
        start_time = time.time()

    #
    # Loop over populations
    #
    for t in range(0,npop):
	i=0
        distance_file=open(folder+'/distance_Population'+repr(t+1)+'.txt',"a")
        if runmode > 1: traj_file=open(folder+'/traj_Population'+repr(t+1)+'.txt',"a")

        beta_values = numpy.zeros( [numOutput] )

        #
        # Loop over accepted particles
        #
        
        
        startTime = time.time()
        while(i<numOutput):
            #param_file = open(folder+'/param'+str(t)+'.txt',"a")
        
            if runmode == 0 or runmode == 1:
                sampling[t]=sampling[t]+1
                if(full == True and sampling[t]%100 == 0):
                    print "sampled: ", sampling[t] 
                    sys.stdout.flush()

                #
                # Sample and perturb the model
                #
                selected_model = sampleTheModel(t, nmodel, margins_prev,  modelDistrib0, dead_models, restart, modelKernel)
            
                #
                # Sample and perturb the parameters given the model
                #
                sampleParameters = sampleTheParameter(t, numOutput, priors[selected_model], kernel[selected_model],
                                                      selected_model, margins_prev, model_prev, weights_prev, parameters_prev, restart, perturbfn ) 

                #
                # Realize the particles
                #
                s=re.compile('SDE')
                g=re.compile('Gillespie')
                if s.search(integrationType[selected_model]) or g.search(integrationType[selected_model]):
                    repeat=beta
                else: 
                    repeat=1

                dist=False
                for r in range(0,repeat):

                    if runmode == 0 :
                        samplePoints=simulateData(ModelName,selected_model,integrationType,InitValues,timepoints,sampleParameters,dt,rtol,atol)
                        points=howToFitData(fit[selected_model],samplePoints)
                        distance=distancefn(points,data)
                    else :
                        distance = realize_particles(  ModelName,
                                                       selected_model,
                                                       integrationType,
                                                       InitValues[selected_model],
                                                       timepoints,
                                                       sampleParameters,
                                                       epsilon,
                                                       t,
                                                       fit[selected_model],
                                                       data,
                                                       distancefn,
                                                       howToFitData,
                                                       dt,
                                                       rtol,
                                                       atol,
                                                       __LARGE_DISTANCE__,
                                                       __XMAX__)

                        #print selected_model, sampleParameters, r, distance 

                    if(distance<=epsilon[0][t] and distance>=0):
                        dist=True
                        distance_file.write(repr(int(hits[t]+1))+'\t'+repr(distance)+'\t'+repr(selected_model+1)+'\n')
                        distance_file.flush()
                        beta_values[i] = beta_values[i] + 1

                # end of loop over B
                
                if dist==True:
                    hits[t]=hits[t]+1
                    if full == True and repeat == 1:
                        print "hits: "+repr(hits[t])+"\n"
                        sys.stdout.flush()
                    elif full == True and repeat > 1:
                        print "hits: "+repr(hits[t])+" : B = "+repr(beta_values[i])+"\n"
                        sys.stdout.flush()

                    rate[t]=hits[t]/sampling[t]

                    #
                    # Fill accepted values
                    #
                    model_curr[i] = selected_model

                    for p in range(0, len(priors[selected_model]) ):
                        parameters_curr[i].append(sampleParameters[p])

                    #
                    # Count how often a specific model was accepted
                    #
                    numbers[selected_model][t]=int(numbers[selected_model][t])+1
                    i=i+1

            if runmode == 2:
                nt = nthread * nblock 
                sampling[t]=sampling[t] + nt
                if full == True:
                    print "population", t, ", cuda batch:", sampling[t], i

                # sample nt models
                batch_model = numpy.zeros( [nt], dtype=numpy.int32 )
                for k in range(nt):
                    batch_model[k] = sampleTheModel(t, nmodel, margins_prev,  modelDistrib0, dead_models, restart, modelKernel)

                # sample nt parameters given the models
                batch_param = numpy.zeros( [nt,pmax], dtype=numpy.float32 )
                for k in range(nt):
                    m = batch_model[k]
                    params = sampleTheParameter(t, numOutput, priors[m], kernel[m],
                                                m, margins_prev, model_prev, weights_prev, parameters_prev, restart, perturbfn ) 

                    for kk in range( model_np[m] ):
                        batch_param[k,kk] = params[kk]

                # Loop over B
                accept = [0 for k in range(0,nt)]
                batch_beta = [0 for k in range(0,nt)]
                for r in range(beta):
                   
                    #
                    # run the ode solver
                    if int_type == 1:
                        simulated_data = abccuda.lsoda.run_lsoda( ptrs, nthread, nblock, nmodel, pmax, xmax, model_nx, 
                                                                  batch_model, batch_param, timepoints, InitValues, dt )
                    #
                    # run the gillespie
                    elif int_type == 3:
                        simulated_data = abccuda.gillespie.run_gillespie( ptrs, nthread, nblock, nmodel, pmax, xmax,
                                                                          batch_model, batch_param, timepoints, InitValues )

                    # calculate distances and fill beta_values
                    for k in range(nt):
                        m = batch_model[k]
                        d = simulated_data[k, :, 0:model_nx[m] ]
                    
                        points=howToFitData( fit[m], d )
                        distance=distancefn(points,data)

                        if custom_kernel == True:
                            print t,",",
                            for kk in range( model_np[m] ):
                                print batch_param[k,kk],",",
                            print distance
                        

                        dist=evaluateDistance(distance,epsilon,t)
                        
                        if dist == True:
                            distance_file.write(repr(int(hits[t]+1))+'\t'+repr(distance)+'\t'+repr(m+1)+'\n')
                            accept[k] = 1
                            batch_beta[k] = batch_beta[k] + 1

                            # output trajectory
                            nrow, ncol = numpy.shape( points )
                            for ic in range(ncol):
                                for ir in range(nrow):
                                    print >>traj_file, points[ir,ic],
                                print >>traj_file, ""

                # end of loop over beta

                # copy accepted particles over
                for k in range(nt):
                    if accept[k] == 1:
                        hits[t]=hits[t]+1
                        rate[t]=hits[t]/sampling[t]
                        
                        if i < numOutput:
                            m = batch_model[k]
                            model_curr[i] = m
                            
                            for p in range(model_np[m]):
                                parameters_curr[i].append( batch_param[k,p] )

                            beta_values[i] = batch_beta[k]
                            
                            #print parameters_curr[i]
                            numbers[m][t]=numbers[m][t] + 1
                            i=i+1
            
            if runmode == 3:
                # ToDo: Some more sophisticated way to chose batchSize
                batchSize = 25000/beta
                
                if(full == True):
                    print "population " +str(t) + ", sampled:", sampling[t], "accepted:", i  
                
                # sample models; save how often specific models were sampled
                batchModel = numpy.zeros(nmodel, dtype=numpy.int32 )
                for k in range(batchSize):
                    batchModel[sampleTheModel(t, nmodel, margins_prev,  modelDistrib0, dead_models, restart, modelKernel)] += 1
                    
                # sample parameters for each model seperately and simulate separately
                allResults = []
                allParameters = []
                allResultsIndices = range(batchSize)
                for k in range(nmodel):
                    if batchModel[k] == 0: continue

                    modelParameters = numpy.zeros((batchModel[k],len(priors[k])), dtype=numpy.float32 )
                    
                    # sample batchModel[k] many parameters for model k
                    for l in range(batchModel[k]):
                        modelParameters[l] = sampleTheParameter(t, numOutput, priors[k], kernel[k],
                            k, margins_prev, model_prev, weights_prev, parameters_prev, restart, perturbfn)
                        #for bla in range(len(modelParameters[l])):
                            #param_file.write(str(modelParameters[l][bla]) + ' ')
                        #param_file.write("\n")
                    
                    # create batchModel[k] many initValues
                    # ToDo allow more than one init Value!
                    species = numpy.zeros((batchModel[k],len(InitValues[k])), dtype=numpy.float32)
                    for l in range(batchModel[k]):
                        species[l,:] = InitValues[k][:]
                    
                    # simulate
                    result = cudaModel[k].run(modelParameters, species, seed=rnd.randint(0,4294967295))
                    
                    # put results from different models together
                    for l in range(len(result)):
                        allResults.append(result[l])
                        allParameters.append(modelParameters[l])
                
                # randomly choose results and test for distance
                rnd.shuffle(allResultsIndices)
                #param_file.close()
                
                print "start distance calculation for batch.. ",
                for k in range(batchSize):
                    # cancel if numOutput has been reached
                    if(i >= numOutput):
                        break
                    
                    tempIndex = allResultsIndices[k]
                    tempResult = allResults[tempIndex]
                    tempModelIndex = _getModelNumber(tempIndex,batchModel)
                    tempParameters = allParameters[tempIndex]
                    sampling[t] += 1
                    
                    # loop over beta
                    foundHit = 0
                    for l in range(len(tempResult)):
                        samplePoints = howToFitData(fit[tempModelIndex],tempResult[l])
                        distance = distancefn(samplePoints, data)
                        
                        dist=evaluateDistance(distance,epsilon,t)

                        if dist == True:
                            distance_file.write(repr(int(hits[t]+1))+", "+ str(foundHit) +'\t'+repr(distance)+'\t'+repr(tempModelIndex)+'\n')
                            foundHit += 1
                            
                            # output trajectory
                            nrow, ncol = numpy.shape( samplePoints )
                            for ic in range(ncol):
                                for ir in range(nrow):
                                    print >>traj_file, samplePoints[ir,ic],
                                print >>traj_file, ""
                        
                    if (foundHit>0):
                        hits[t] += 1
                        rate[t] = hits[t]/sampling[t]

                        model_curr[i] = tempModelIndex
                        parameters_curr[i] = tempParameters
                        #input(parameters_curr[i])
                        beta_values[i] = foundHit
                        numbers[tempModelIndex][t] += 1
                        
                        i += 1
                
            # end of loop over particles
        distance_file.close()
        if runmode > 1 : traj_file.close()

        #
        # Compute the weights
        #  
        for i in range(numOutput):
            selected_model = model_curr[i]
            weights_curr[i] = computeParticleWeights(t, parameters_curr[i], selected_model, nmodel, numOutput, priors[selected_model], kernel[selected_model],
                                                     model_prev, margins_prev, weights_prev, parameters_prev, beta_values[i], dead_models, restart, kernelpdffn)

            #print selected_model, parameters_curr[i], weights_curr[i]

        #
        # Normalize the weights
        #
        weights_curr = normalizeWeights(numOutput, weights_curr)

        #
        # Get the marginal model distribution
        #  
        margins_curr = modelMarginals(nmodel, numOutput, model_curr, weights_curr)

        #
        # Prepare for next population
        # 
        margins_prev = margins_curr[ : ]
        weights_prev = weights_curr[ : ]
        model_prev = model_curr[ : ]
        parameters_prev = []
        for i in range(numOutput):
            parameters_prev.append(parameters_curr[i][ : ])

        model_curr      = [0  for j in range(0,numOutput)]
        weights_curr    = [0  for j in range(0,numOutput)]
        parameters_curr = [[] for j in range(0,numOutput)]
        margins_curr    = [0  for j in range(0,nmodel)] 

        #
        # Check for dead models
        #
        dead_models = []
        for j in range(nmodel):
            if margins_prev[j] < 1e-6:
                dead_models.append(j)


        #
        # Fill global containers
        #
        for mod in range(nmodel):
            for part in range(numOutput):
                if int(model_prev[part]) == mod:
                    weights[mod][t].append( weights_prev[part] )

                    # Marginals
                    modelDistribution[t][mod] = modelDistribution[t][mod] + weights_prev[part]
                    numbers[mod][t] = numbers[mod][t] + weights_prev[part]

                    for param in range(len(priors[mod])):
                        population[mod][t][param].append( parameters_prev[part][param] )

        #
        # Compute kernels
        #
        if constKernels==False: ### numbers[lastpopulation]
            for mod in range(0, nmodel):
                if margins_prev[mod] > 0.1:  ### only compute corrected kernels if the model was selected at least 10%
                    kernel[mod] = kernelfn(kernel[mod], population[mod][t])

        if pickling==True:
            write_pickled(folder, nmodel, numOutput, priors, model_prev, weights_prev, parameters_prev, margins_prev, kernel)

        if full == True:
            print "number of sampling steps in population",t+1,":", sampling[t]
            print "model marginals for iteration",t+1,":",margins_prev
            print "acceptance rate:",rate[t]
            if(len(dead_models) > 0): print "dead_models:", dead_models

        rate_file = open(folder+'/rates.txt',"a")
        rate_file.write(repr(t+1)+"\t"+repr(sampling[t])+"\t"+repr(rate[t]))
        rate_file.write("\t" + str(round((time.time()-startTime)/60,2)) + " min\n")
        rate_file.close()

        if timing == True:
            print "####abcSimulator:runmode/pop/time:", runmode, t+1, time.time() - start_time

        #
        # Print histograms and scatter plots (analytical plots)
        #
        
        make_analysis_plots(folder, diagnostic, plotDataSeries, timepoints, fit, dt, data, t, nmodel, ModelName,
                            integrationType, InitValues, priors, population, numbers, modelDistribution, weights, epsilon[0], rate )            
        restart=False

        # end of loop over populations

    if timing == True:
        print "####abcSimulator:runmode/time:", runmode, time.time() - start_time
        
    rate_file.close()

    return(modelDistribution,population,rate)
