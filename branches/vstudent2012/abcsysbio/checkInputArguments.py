import pickle
import re
import numpy



def checkInputABC(info_new, fname, custom_distance, design ):
                  
    """
    Check that the information in the input file is consistent with each other and with the model,
    and that it is in the format required to run the abc-SMC algorithm.
    Return boolean, string (empty if boolean is True)

    Takes as input an algorithm_info object, the output folder name, whether custom distance is specified, and whether the mode is design
    
    """

    restart =            info_new.restart    #True/False   
    ModelName =         info_new.name
    subModelName =      info_new.submodelname  
    data =              info_new.data   
    timepoints =        info_new.times
    numOutput =         info_new.particles
    epsilon =           info_new.epsilon
    integrationType =   info_new.type
    kernel =            info_new.kernel        
    modelKernel =       info_new.modelkernel   
    modelWeights =      info_new.modelprior    
    priors =            info_new.gprior    
    x0priors =          info_new.x0prior    
    source =            info_new.source
    fit =               info_new.fit
    beta =              info_new.beta
    dt =                info_new.dt
    rtol =              info_new.rtol
    atol =              info_new.atol
    localpriors =       info_new.localprior
    
    ### check general properties of the given arguments that are independent of the model


    # Model names must be given for every model:
    for i in range(0, len(ModelName)):
        if ModelName[i] == "":
            return False, "\nPlease do not give empty strings for model names!\n"

    # To accommodate submodels per model:
    for i in range(0, len(subModelName)):
        for j in range(0,len(subModelName[i])):
            if subModelName[i][j] == "":
                return False, "\nPlease do not give empty strings for submodel names!\n"

    # Number of models = number of types:
    if not len(ModelName)==len(integrationType): 
        return False,"\nPlease provide the same number of models and integration types!\n"

    # Number of models = number of weights:
    if not len(ModelName)==len(modelWeights): 
        return False,"\nPlease provide the same number of models and model weights!\n"

    # Fit instruction or 'None' for each model:
    if not len(ModelName)==len(fit):
        return False,"\nPlease provide a fit instruction (or None) for each model. If the fit instruction is None all data will be fitted to the model data.\n"

    # If all parameters and initial conditions are constant, parameter inference cannot occur!
    for i in range(0,len(priors)):
        count = 0
        for j in range(0,len(priors[i])):
                if priors[i][j][0] != 0:
                    count = count + 1
        for j in range(0,len(x0priors[i])):
            for k in range(0,len(x0priors[i][j])):
                if x0priors[i][j][k][0] != 0:
                    count = count + 1
        for j in range(0,len(localpriors[i])):
            for k in range(0,len(localpriors[i][j])):
                if localpriors[i][j][k][0] != 0:
                    count = count + 1
        if count == 0:
            return False,"\nThere must be at least one non-constant parameter per model for parameter inference and model selection.\n"

    # Timpoints must be the same length as variable values list in each dataset:
    if design == False:
        for i in (0,len(data)/2):
            if not len(data[i])==len(timepoints[i]):
                return False,"\nPlease provide data that correspond to the length of your timepoints!\n"

    # If fit is not defined for a given submodel, the number of initial values must be the same as the number of variables in the data set:
    for i in range(0,len(x0priors)):
        for j in range(0,len(x0priors[i])):
            if fit[i][j] == None:
                if len(x0priors[i][j]) != len(data[j][0]):
                    return False,"\nIf fit is not defined per submodel, please ensure that the number of initial values matches the number of variables in data.\n"

    # checks that prior distributions are provided for every model
    if not len(ModelName)==len(priors):
        return False,"\nPlease provide prior distributions for each model!\n"




    # Initial values for each submodel:
    totalIniValSets = 0
    totalModels = 0

    for i in range(0,len(x0priors)):
        totalIniValSets = totalIniValSets + len(x0priors[i])
    for i in range(0,len(subModelName)):
        totalModels=totalModels + len(subModelName[i])

    if not totalModels==totalIniValSets:
        return False,"\nPlease provide initial values for each submodel!\n"


    # Check that a kernel number is supplied:
    if not modelKernel>0.0:
        return False,"\nPlease provide a model Kernel larger than 0!\n"


    # User must supply a custom distance function when using multiple epsilon schedules:
    if (len(epsilon)>1 and  custom_distance==False):
        return False,"\nPlease provide a custom distance function when you specify more than 1 epsilon schedule!\n"


    # The number of submodels must be the same for each model:
    for i in range (1,len(subModelName)):
        if not len(subModelName[0])==len(subModelName[i]):
            return False, "\nPlease provide the same number of submodels, in the same order, for each model!\n"


    # The number of submodels = number of datasets (indirectly ensures that there are not more datasets than source files!):
    if not len(subModelName[0])==len(data):
        return False,"\nPlease provide the same number of data files as submodels!\n"
        
        
########## checking further properties independent of the model #########################################################

    sde=re.compile('SDE')
    ode=re.compile('ODE')
    gillespie=re.compile('Gillespie')
            
    for mod in range(0,len(ModelName)):
        # Integrations types must be legitimate
        string=integrationType[mod]
        if (not sde.search(string)) and (not ode.search(string)) and (not gillespie.search(string)):
            return False,"\nThe integration type for model "+ModelName[mod]+" does not exist!\n"        

        # Check for correct definition of parameter distributions:
        for param in range(0,len(priors[mod])):
            if len(priors[mod][param])!=3 and (len(priors[mod][param])!=4 and priors[mod][param][0]!=4):
                return False, "\nThe prior distribution of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is incorrectly defined!\n"

            # Checks prior distributions are legitimate:
            if not (priors[mod][param][0]==0 or priors[mod][param][0]==1 or priors[mod][param][0]==2 or priors[mod][param][0]==3 or priors[mod][param][0]==4):
                return False, "\nThe prior distribution of parameter "+repr(param+1)+" in model "+ModelName[mod]+" does not exist!\n"

            # There must be a prior distribution specified for each of the inivals for each submodel:
            for i in range(len(x0priors[mod])):
                for j in range(len(x0priors[mod][i])):
                    if not (x0priors[mod][i][j][0]==0 or x0priors[mod][i][j][0]==1 or x0priors[mod][i][j][0]==2 or x0priors[mod][i][j][0]==3 or x0priors[mod][i][j][0]==4):
                        return False, "\nThe prior distribution of initial condition "+repr(j+1)+" in "+subModelName[mod][i]+" does not exist!\n"      


            # Uniform distributions must have legitimate parameters:
            if priors[mod][param][0]==2:
                if not priors[mod][param][1]<priors[mod][param][2]:
                    return False, "\nThe range of the uniform prior distribution of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is incorrectly defined!\n"

            # Checks the definition of parameters of prior distributions:
            if priors[mod][param][0]==3:
                if not (priors[mod][param][1]>=0 or priors[mod][param][2]>=0):
                    return False, "\nThe mean or scale of the lognormal prior distribution of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is incorrectly defined!\n"


        ##if kernel[mod][param][0]==1 and not priors[mod][param][0]==0:
            ##    if constKernels==False:
            ##        if not kernel[mod][param][1]<kernel[mod][param][2]:
            ##            return False, "\nThe range of the uniform pertubation kernel of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is wrong defined!\n"
        
        ##if kernel[mod][param][0]==2 and not priors[mod][param][0]==0:
            ##    if constKernels==False:
            ##        if not kernel[mod][param][2]>0:
            ##            return False, "\nThe variance of the gaussian pertubation kernel of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is wrong defined!\n"



############## check model specific properties (comparing with SBML model) ######################################

    if not source==None:
        import libsbml

        # Check number of submodels = number of SBML files associated with this model:
        totalSource = len(source)*len(source[0])
        if not totalSource==totalModels: 
            return False,"\nPlease provide the same number of submodel sources and submodel names!\n"
        

        # Iterate through each source file via source array to get pararmeter and species information:
        reader=libsbml.SBMLReader()
        for mod in range (0,len(source)):
            for l in range(0,len(source[mod])):
        
                document=reader.readSBML(source[mod][l])
                model=document.getModel()

                numSpecies=model.getNumSpecies()
                numGlobalParameters=model.getNumParameters()
                listOfParameter=[]

                NumCompartments=model.getNumCompartments()   
            
                for i in range(0,NumCompartments):
                    if model.getCompartment(i).isSetVolume():
                        numGlobalParameters=numGlobalParameters+1
                        listOfParameter.append(model.getListOfCompartments()[i])

                for i in range(0, numGlobalParameters-NumCompartments):
                    listOfParameter.append(model.getParameter(i))


                numLocalParameters=0
                for i in range(0,model.getNumReactions()):
                    local=model.getReaction(i).getKineticLaw().getNumParameters()
                    numLocalParameters=numLocalParameters+local
                    for k in range(0,local):
                        listOfParameter.append(model.getListOfReactions()[i].getKineticLaw().getParameter(k))
                              

                numParameters=numLocalParameters+numGlobalParameters

                species = model.getListOfSpecies()
                #for k in range(0, len(species)):
                    #if (species[k].getConstant() == True):
                    #    numParameters=numParameters+1
                    #    numSpecies=numSpecies-1

                listOfRules = model.getListOfRules()
                for k in range (0, len(listOfParameter)):
                    if listOfParameter[k].getConstant()==False:
                        for j in range(0, len(listOfRules)):
                            if listOfRules[j].isRate():
                                if listOfParameter[k].getId()==listOfRules[j].getVariable():
                                    numSpecies=numSpecies+1
                                    numParameters=numParameters-1


                # Checks parameters for each submodel against the proposed parameters of the model.
                # Checks that all of the submodels have the same number of parameters. 
                if not len(info_new.pindex[mod][0][0])==numParameters:
                    return False,"\nThe number of given prior distributions for "+ModelName[mod]+" is not correct! This may be because you have supplied submodels with different parameter sets. Parameter sets for all submodels within a model MUST be the same!\n"


                if not len(x0priors[mod][l])==numSpecies:
                    return False,"\nPlease provide an initial value for each species in submodel "+subModelName[mod][l]+"!\n"


    ### check arguments connected to pickling
    if restart==True:

        try:
            in_file=open(fname + '/copy/algorithm_parameter.dat',"r")
            numOutput_pickled=pickle.load(in_file)
            in_file.close()
        except:
            return False,"\nCan not find file \'algorithm_parameter.dat\' in folder \'copy\'!\n"
        
        if not numOutput<=numOutput_pickled:
            return False, "\nRunning the abc algorithm from a previous point is not possible with a larger population size!\n"



    return True,""


##########################################################################################################
##def checkInputSimulation(name, timepoints, InitValues, IntegrationType, ConstantParameters, source, dt):
def checkInputSimulation(info_new , fname):
    """
    Check that the information in the input file is consistent with each other and with the model,
    and that it is in the format required to simulate the model.
    Return boolean, string (empty if boolean is True)
    
    """
    ModelName =         info_new.name
    subModelName =      info_new.submodelname
    timepoints =        info_new.times
    IntegrationType =   info_new.type
    modelWeights =      info_new.modelprior
    priors =            info_new.prior
    x0priors =          info_new.x0prior
    source =            info_new.source
    fit =               info_new.fit
    beta =              info_new.beta
    dt =                info_new.dt


    # Check all models have names:
    for i in range(0, len(ModelName)):
        if ModelName[i] == "":
            return False, "\nPlease do not give empty strings for model names!\n"

    # Check all submodels have names:
    for i in range(0, len(subModelName)):
        for j in range(0,len(subModelName[i])):
            if subModelName[i][j] == "":
                return False, "\nPlease do not give empty strings for submodel names!\n"

    # Number of models = number of types:
    if not len(ModelName)==len(IntegrationType): 
        return False,"\nPlease provide the same number of models and integration types!\n"



########### check model specific properties (comparing with SBML model) #########################################

    if not source==None:
        import libsbml

        totalModels = 0

        for i in range(0,len(subModelName)):
            totalModels=totalModels + len(subModelName[i])

        totalSource = len(source)*len(source[0])
        if not totalSource==totalModels: 
            return False,"\nPlease provide the same number of submodel sources and submodel names!\n"


        # Collect parameter and species information from model sets:
        
        reader=libsbml.SBMLReader()
        for mod in range(0,len(source)):
            for l in range(0,len(source[mod])):    

                document=reader.readSBML(source[mod][l])
                model=document.getModel()

                numSpecies=model.getNumSpecies()
                numGlobalParameters=model.getNumParameters()

                NumCompartments=model.getNumCompartments()   
                for i in range(0,NumCompartments):
                    if model.getCompartment(i).isSetVolume():
                        numGlobalParameters=numGlobalParameters+1

                numLocalParameters=0
                for i in range(0,model.getNumReactions()):
                    numLocalParameters=numLocalParameters+model.getReaction(i).getKineticLaw().getNumParameters()

                numParameters=numLocalParameters+numGlobalParameters

                species = model.getListOfSpecies()
                #for k in range(0, len(species)):
                #    if (species[k].getConstant() == True):
                #        numParameters=numParameters+1
                #        numSpecies=numSpecies-1

                if not info_new.nparameters[mod] == numParameters:
                    return False,"\nThe number of given parameters for "+ModelName[mod]+" is not correct! This may be because you have supplied submodels with different parameter sets. Parameter sets for all submodels within a model MUST be the same!\n"
            
                if not len(info_new.x0prior[mod][l]) == numSpecies:
                    return False,"\nPlease provide an initial value for each species in submodel "+subModelName[mod][l]+"!\n"

    # Timepoints must be provided:
    if len(timepoints)==0:
        return False, "\nPlease give timepoints at which to return simulated data points!\n"

    # Check integration types:
    sde=re.compile('SDE')
    ode=re.compile('ODE')
    gillespie=re.compile('Gillespie')
    for mod in range(0,len(ModelName)):

        string=IntegrationType[mod]
        if (not sde.search(string)) and (not ode.search(string)) and (not gillespie.search(string)):
            return False,"\nThe integration type for "+ModelName[mod]+" does not exist!\n"


    return True,""




