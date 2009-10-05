import re
from numpy import *
from scipy.integrate.odepack import odeint
from abc.sdeint import *
from abc.getResults import *
from scipy.stats import *
from abc import GillespieAlgorithm
from abc import abcodesolve
from math import *
from numpy import random as rnd
from abc import euclidian
from abc.checkInputArguments import *
import os
import pickle


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


################## compute the kernels
def getKernel(kern,population):

    for param in range(0, len(population)):
        minimum=min(population[param])
        maximum=max(population[param])
        scale=(maximum-minimum)/4.0
        if kern[param][0]==1: 
            kern[param][1]=-scale/2.0
            kern[param][2]=scale/2.0
        if kern[param][0]==2:
            kern[param][2]=scale/2.0

    return kern


##################compute the pdf of uniform distribution
def getPdfUniform(scale1,scale2,parameter):
    if ((parameter>scale2) or (parameter<scale1)): x=0.0
    else:
        scale=scale2-scale1
        x=(scale*2)**(-1)

    return x


##################compute the pdf of gauss distribution
def getPdfGauss(mean,scale,parameter):
    x=((math.sqrt(2*math.pi)*scale)**(-1))*math.exp(-1.0*((parameter-mean)**2)/(2.0*scale**2))

    return x


################compute the pdf of lognormal distribution
def getPdfLognormal(mean,sigm,parameter):
    x=((sqrt(2*math.pi)*sigm*parameter)**(-1))*exp(-0.5*(((math.ln(parameter)-sigm)**2)/sigm)**2)
    
    return x

################ sample parameters from prior
def sampleParametersFromPrior(priors):
    sampleParameters=zeros([len(priors)])
    for n in range(0,len(priors)):
        if priors[n][0]==0: sampleParameters[n]=priors[n][1]
    #   if priors[n][0]==1: sampleParameters[n]= sample from Jeffrey's prior
        if priors[n][0]==2: sampleParameters[n]=rnd.uniform(low=priors[n][1],high=priors[n][2])
        if priors[n][0]==3: sampleParameters[n]=rnd.lognormal(priors[n][1],priors[n][2])
        
    return sampleParameters

################sample parameters from previous population with weights
def sampleParametersFromPrevious(numbers,population,weights,kernel,priors,sampleParameters):

    xEnd=array(0.0)
    while(xEnd<=0):
        xAll=array(1.0)
        for n in range(0,len(priors)):
            if not priors[n][0]==0:
                line=w_choice(population[n],weights[n])
                break
        ###add pertubation to parameter:
        for n in range(0,len(priors)):
            if priors[n][0]==0: 
                sampleParameters[n]=priors[n][1]
                x=array(1.0)
            else:
                if kernel[n][0]==1:
                    sampleParameters[n]=population[n][line]+rnd.uniform(low=kernel[n][1],high=kernel[n][2])
                if kernel[n][0]==2: 
                    sampleParameters[n]=population[n][line]+rnd.normal(kernel[n][1],kernel[n][2])

                if priors[n][0]==2: 
                    x=getPdfUniform(priors[n][1],priors[n][2],sampleParameters[n])
                if priors[n][0]==3: 
                    x=getPdfLognormal(priors[n][1],priors[n][2],sampleParameters[n])

            xAll=xAll*x
        xEnd=xEnd+xAll
    
    return sampleParameters



############sample parameters
def sampleTheParameter(priors,sampleParameters,numbers,t,population,weights,kernel,restart):
    if ((restart==False and t==0) or numbers==0):
        sampleParameters=sampleParametersFromPrior(priors)

    else: 
        sampleParameters=sampleParametersFromPrevious(numbers,population,weights,kernel,priors,sampleParameters)
    return sampleParameters


###########simulate data with sampled parameter set for the sampled model
def simulateData(ModelName,selected_model,integrationType,InitValues,timepoints,sampleParameters,dt,rtol=None,atol=None):

    o=re.compile('ODE')
    s=re.compile('SDE')
    g=re.compile('Gillespie')
  
    module = __import__(ModelName[selected_model])


    if o.match(integrationType[selected_model]):
        samplePoints=abcodesolve.abcodeint(module, InitValues[selected_model], timepoints, sampleParameters, dt, rtol, atol)
    if s.match(integrationType[selected_model]): 
        samplePoints=sdeint(module,InitValues[selected_model],sampleParameters,timepoints,dt)
    if g.match(integrationType[selected_model]): 
        samplePoints=GillespieAlgorithm.GillespieInt(module, InitValues[selected_model], sampleParameters,timepoints)

    return samplePoints


###########how to fit variable to the given data (if fit is None, data for all variables are available in order of the model)
def howToFitData(fitting_instruction,samplePoints):

    if fitting_instruction!=None:
        points=zeros([len(samplePoints), len(fitting_instruction)])
        for i in range(0, len(fitting_instruction)):
            points[:,i] = eval(fitting_instruction[i])
           
    else:
        points=samplePoints


    return points

############compute weights of particles
def computeParticleWeights(priors,sampleParameters,numbers1,numbers0,t,kernel,population1,population0,weights1,weights0,BETA,restart):
    if(((restart==False and t>0) or restart==True) and numbers0>0):
        xAll=1.0
        for n in range(0,len(priors)):
            if priors[n][0]==0:
                x=1
            if priors[n][0]==2: 
                x=getPdfUniform(priors[n][1],priors[n][2],sampleParameters[n])
            if priors[n][0]==3: 
                x=getPdfLognormal(priors[n][1],priors[n][2],sampleParameters[n])
            xAll=xAll*x
        sum=0
           
       
       
        for p in range(0,int(numbers0)):
            kernAll=0.0
            for n in range(0,len(priors)):
                if not(priors[n][0]==0):

                    if kernel[n][0]==1: 
                        scale=(kernel[n][2]-kernel[n][1])/2.0
                        scale1=max((population0[n][p]-scale),priors[n][1])
                        scale2=min((population0[n][p]+scale),priors[n][2])
                        kern=getPdfUniform(scale1,scale2,population1[n][int(numbers1)])

                    if kernel[n][0]==2:
                        mean=population0[n][p]
                        scale=kernel[n][2]
                        CDF2 = norm.cdf(priors[n][2],mean,scale)
                        CDF1 = norm.cdf(priors[n][1],mean,scale)
                        CDF=(CDF2-CDF1)**(-1)
                        kern=getPdfGauss(mean,scale,population1[n][int(numbers1)])
                        kern=kern*CDF

                else: kern=1.0
                kernAll=kernAll+kern
            sum=sum+weights0[0][p]*kernAll
        for n in range(0,len(priors)):
            weights1[n].append(xAll*BETA/sum)



    return weights1
     

################normalize the weights
def normalizeWeights(len_ModelName,priors,t,numbers,weights,restart):
    for mod in range(0,len_ModelName):
        if numbers[mod][t]>0:
            for h in range(0,len(priors[mod])):
                weight=0
                for j in range(0,int(numbers[mod][t])):
                    if (t==0 or numbers[mod][t-1]==0) and restart==False: weights[mod][t][h].append(1.0)
                    weight=weight+weights[mod][t][h][j]

                for j in range(0,int(numbers[mod][t])):
                    weights[mod][t][h][j]=weights[mod][t][h][j]/float(weight)
    return weights






#######################################################################
# The ABC SMC algorithm for model selection                           #
#######################################################################

def abcSimulator(ModelName,data,timepoints,numOutput, epsilon, InitValues,integrationType,modelWeights,priors,kernel,fit=None,sampleFromPrior=False,beta=1,dt=0.01,restart=False,source=None,diagnostic=True,pickling=True,plotDataSeries=True,full=False, rtol=None, atol=None, constKernels=False):   

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
              a list of floats.
              Epsilon corresponds to the maximum distance between experimental and simmulated 
              data. The length of this list corresponds to the number of populations.

              Example:
              maximum distance for ABC Rejection (only one population):
              epsilon = (1.5,)
              maximum distance for ABC SMC with five populations:
              epsilon = (3.5,3.0,2.5,2.0,1.5)

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
    
    modelCorrect,message=checkInputABC(ModelName,
                                       data,
                                       timepoints,
                                       numOutput,
                                       epsilon,
                                       InitValues,
                                       integrationType,
                                       modelWeights,
                                       priors,
                                       kernel,
                                       fit,
                                       restart,
                                       source,
                                       constKernels)

    if modelCorrect==False:
        print message
        sys.exit()

    modelDistribution=zeros([len(epsilon),len(ModelName)])


###########################create a result folder for each model

    if restart==True:
        folder='_restart_results_'
    else:
        folder='_results_'
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
            print "\nThe folder "+ folder + ModelName[mod]+" already exists!\n"
            sys.exit()

    os.chdir('..')

            
    rate_file=open(folder+'/rates.txt',"w")

        
    if pickling==True:
        try:
            os.chdir(folder)
            os.mkdir('copy')
            os.chdir('..')
        except:
            print "\nThe folder \'copy\' already exists!\n"
            sys.exit()
    

##########################save algorithm parameter via pickling
    if pickling==True:
        out_file=open(folder+'/copy/algorithm_parameter.dat',"w")
        pickle.dump(numOutput,out_file)
        out_file.close()


# population=empty([len(ModelName),len(epsilon),parameternumber,numOutput for each parameter])
    population=[]
    weights=[]

    ModelNumber=arange(0,len(ModelName),1)
    numbers=zeros([len(ModelName),len(epsilon)])
    modelDistrib0=zeros([len(ModelName)])
    modelDistrib=zeros([len(ModelName)])

    hits=[]
    sampling=[]
    rate=[]
    
    for j in range(0,len(ModelName)):
        population.append([])
        weights.append([])
        for eps in range(0,len(epsilon)):
            population[j].append([])
            weights[j].append([])
            for param in range(0,len(priors[j])):
                population[j][eps].append([])
                weights[j][eps].append([])                              

    for j in range(0,len(epsilon)):
        hits.append(0.0)
        sampling.append(0.0)
        rate.append(0.0)

################normalize the prior model distribution
    weight=0
    if restart==True:

        # pickle numbers selected model of previous population
        # pickle population of selected model of previous population pop_pickled[selected_model][n][vectpos]
        # pickle weights of selected model of previous population weights_pickled[selected_model][n][vectpos]
        
        try:
            in_file=open('_results_/copy/numlast.dat',"r")
            num_pickled=pickle.load(in_file)
            in_file.close()
        except:
            print "\nCan not find file \'numlast.dat\' in folder \'copy\'!\n"
            sys.exit()
        
        try:
            in_file=open('_results_/copy/poplast.dat',"r")
            pop_pickled=pickle.load(in_file)
            in_file.close()
        except:
            print "\nCan not find file \'poplast.dat\' in folder \'copy\'!\n"
            sys.exit()

        try:
            in_file=open('_results_/copy/weightslast.dat',"r")
            weights_pickled=pickle.load(in_file)
            in_file.close()
        except:
            print "\nCan not find file \'weightslast.dat\' in folder \'copy\'!\n"
            sys.exit()
        

        try:
            in_file=open('_results_/copy/kernelslast.dat',"r")
            kernel=pickle.load(in_file)
            in_file.close()
        except:
            print "\nCan not find file \'kernelslast.dat\' in folder \'copy\'!\n"
            sys.exit()


    for mod in range(0,len(modelWeights)):
        weight=weight+modelWeights[mod]
    for mod in range(0,len(ModelName)):
        modelDistrib0[mod]=modelWeights[mod]*weight**(-1)


################start model selection with external loop
    for t in range(0,len(epsilon)):
	i=0
        distance_file=open(folder+'/distance_Population'+repr(t+1)+'.txt',"w")

#################count how many particles are accepted
	while(i<numOutput):
            sampling[t]=sampling[t]+1

################sample model from prior or intermediate distribution         
            if(sampleFromPrior==True or (t==0 and restart==False)): selected_model=w_choice(ModelNumber,modelDistrib0)
            else: #sample from previous
                if restart==True: #restart and sample from previous
                    weight=0
                    for mod in range(0,len(num_pickled)):
                        weight=weight+num_pickled[mod]
                    for mod in range(0,len(ModelName)):
                        modelDistrib[mod]=num_pickled[mod]*weight**(-1)
                
                if restart==False: #normal run (or t>0 for restarted run) and sample from previous
                    weight=0
                    for mod in range(0,len(ModelName)):
                        weight=weight+modelDistribution[t-1][mod]
                    for mod in range(0,len(ModelName)):
                        modelDistrib[mod]=modelDistribution[t-1][mod]/weight
                
                selected_model=w_choice(ModelNumber,modelDistrib) # sampleFromPrior is False AND t is larger than 0

            dim=len(priors[selected_model])
            sampleParameters=[]
            for z in range(0,dim):
                sampleParameters.append(0)

#############sample parameters from prior or from previous population
            dist=False
            BETA=0
            if restart==False:
                sampleParameters=sampleTheParameter(priors[selected_model],
                                                    sampleParameters,
                                                    numbers[selected_model][t-1],
                                                    t,
                                                    population[selected_model][t-1],
                                                    weights[selected_model][t-1],
                                                    kernel[selected_model],
                                                    restart)
            else:
                sampleParameters=sampleTheParameter(priors[selected_model],
                                                    sampleParameters,
                                                    num_pickled[selected_model],
                                                    t,
                                                    pop_pickled[selected_model],
                                                    weights_pickled[selected_model],
                                                    kernel[selected_model],
                                                    restart)

############simulate data with sampled parameter set for the sampled model
            s=re.compile('SDE')
            g=re.compile('Gillespie')
            if s.search(integrationType[selected_model]) or g.search(integrationType[selected_model]):
                repeat=beta
            else: 
                repeat=1
            for r in range(0,repeat):
                samplePoints=simulateData(ModelName,selected_model,integrationType,InitValues,timepoints,sampleParameters,dt,rtol,atol)

###########how to fit variable to the given data (if fit is None, data for all variables are available in order of the model)
                points=howToFitData(fit[selected_model],samplePoints)

###########compute the distance and compare it with epsilon
                distance=euclidian.euclidianDistance(points,data)
                if(distance<=epsilon[t]):
                    distance_file.write(repr(int(hits[t]+1))+'\t'+repr(distance)+'\t'+repr(selected_model+1)+'\n')
                    dist=True
                    BETA=BETA+1
                    hits[t]=hits[t]+1
                    if full == True:
                        print "hits: "+repr(hits[t])+"\n"
                    rate[t]=hits[t]/sampling[t]
            
            if dist==True:
                for param in range(0,dim):
                    population[selected_model][t][param].append(sampleParameters[param])

############compute weights of particles
                if restart==False:
                    weights[selected_model][t]=computeParticleWeights(priors[selected_model],
                                                                      sampleParameters,
                                                                      numbers[selected_model][t],
                                                                      numbers[selected_model][t-1],
                                                                      t,
                                                                      kernel[selected_model],
                                                                      population[selected_model][t],
                                                                      population[selected_model][t-1],
                                                                      weights[selected_model][t],
                                                                      weights[selected_model][t-1],
                                                                      BETA,
                                                                      restart)
                    
                else:
                    weights[selected_model][t]=computeParticleWeights(priors[selected_model],
                                                                      sampleParameters,
                                                                      numbers[selected_model][t],
                                                                      num_pickled[selected_model],
                                                                      t,
                                                                      kernel[selected_model],
                                                                      population[selected_model][t],
                                                                      pop_pickled[selected_model],
                                                                      weights[selected_model][t],
                                                                      weights_pickled[selected_model],
                                                                      BETA,
                                                                      restart)


########count how often a specific model was accepted
                numbers[selected_model][t]=int(numbers[selected_model][t])+1
                i=i+1


        distance_file.close()
##############normalize the weights
        weights=normalizeWeights(len(ModelName),priors,t,numbers,weights,restart)

########compute new kernels in respect to the previous population
        if constKernels==False: ### numbers[lastpopulation]
            for mod in range(0, len(ModelName)):
                if numbers[mod][t]>numOutput*0.1:  ### only compute corrected kernels if the model was selected at least 10%
                    kernel[mod]=getKernel(kernel[mod],population[mod][t])

      
#############get the model distribution           
        for mod in range(0,len(ModelName)):
            modelDistribution[t][mod]=int(numbers[mod][t])



#############save last population and last model distribution in file via pickling
        if pickling==True:
            out_file=open(folder+'/copy/poplast.dat',"w")
            x=[]
            for mod in range(0,len(ModelName)):
                x.append(population[mod][t])
            pickle.dump(x,out_file)
            out_file.close()

            out_file=open(folder+'/copy/numlast.dat',"w")
            x=[]
            for mod in range(0,len(numbers)):
                x.append(numbers[mod][t])
            pickle.dump(x,out_file)
            out_file.close()

            out_file=open(folder+'/copy/weightslast.dat',"w")
            x=[]
            for mod in range(0,len(ModelName)):
                x.append(weights[mod][t])
            pickle.dump(x,out_file)
            out_file.close()

            out_file=open(folder+'/copy/kernelslast.dat',"w")
            x=[]
            for mod in range(0,len(ModelName)):
                x.append(kernel[mod])
            pickle.dump(x,out_file)
            out_file.close()

############print numerical parameters and measures
        if full == True:
            print "number of sampling steps in population",t+1,":", sampling[t]
            print "accepted model for iteration",t+1,":",hits[t]
            print "acceptence rate:",rate[t],"\n"

        rate_file.write(repr(t+1)+"\t"+repr(sampling[t])+"\t"+repr(rate[t])+"\n")


###########print histograms and scatter plots (analytical plots)
#WeightedHistogram(|data|, |weights|, |bins|, |range|=None)


        
        for mod in range(0,len(ModelName)):
            try:
                os.chdir(folder+'/results_'+ModelName[mod])
                os.mkdir("Population_"+repr(t+1))
                os.chdir("../..")
            except:
                print "\nCan not creat the folder Population_"+repr(t+1)+"!\n"
                sys.exit()


        for mod in range(0,len(ModelName)):
            weight_file=open(folder+'/results_'+ModelName[mod]+'/Population_'+repr(t+1)+'/data_Weights'+repr(t+1)+".txt","w")
            for g in range(0,len(weights[mod][t][0])):
                for j in range(0, len(weights[mod][t])): #for each entry in the first parameter
                    weight_file.write(repr(weights[mod][t][j][g])+"\t")
                weight_file.write("\n")
            weight_file.close()

            
        population_mod=[]
        weights_mod=[]
        for mod in range(0,len(ModelName)):
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
                            weights_mod[mod][eps].append(weights[mod][eps][param])

                getAllScatterPlots(population_mod,weights_mod,populations=arange(1,t+2,1),PlotName=PlotName,model=mod+1)
                getAllHistograms(population_mod,weights_mod,population=t+1,PlotName=PlotName2, model=mod+1)

        if len(ModelName)>1:
            printModelDistribution(modelDistribution,eps=t,filename=folder+'/ModelDistribution.txt')
            if diagnostic==True:
                getModelDistribution(modelDistribution,epsilon,rate,PlotName=folder+'/ModelDistribution')



        if plotDataSeries==True:
            for mod in range(0,len(ModelName)):
                x=10
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

            
        restart=False
       

    rate_file.close()

    return(modelDistribution,population,rate)
