import pickle
import re

def checkInputABC(info_new , fname, custom_distance, design ):
                  
    """
    Check that the information in the input file is consistent with each other and with the model,
    and that it is in the format required to run the abc-SMC algorithm.
    Return boolean, string (empty if boolean is True)

    Takes as input an algorithm_info object, the output folder name, and whether custom distance is specified
    
    """

    restart = 	        info_new.restart
    ModelName = 	info_new.name
    data = 		info_new.data
    timepoints = 	info_new.times
    numOutput = 	info_new.particles
    epsilon = 	        info_new.epsilon
    integrationType =   info_new.type
    modelWeights = 	info_new.modelprior
    priors = 	        info_new.prior
    x0priors =          info_new.x0prior
    kernel = 	        info_new.kernel
    source = 	        info_new.source
    fit = 		info_new.fit
    beta = 		info_new.beta
    dt = 		info_new.dt
    rtol = 		info_new.rtol
    atol = 		info_new.atol
    modelKernel = 	info_new.modelkernel

    
    ### check general properties of the given arguments that are independent of the model

    for i in range(0, len(ModelName)):
        if ModelName[i] == "":
            return False, "\nPlease do not give empty strings for model names!\n"

    if not len(ModelName)==len(integrationType): 
	return False,"\nPlease provide the same amount of model sources and integration types!\n"
    
    if not len(ModelName)==len(modelWeights): 
        return False,"\nPlease provide the same amount of model sources and model weights!\n"

    if not len(ModelName)==len(fit):
        return False,"\nPlease provide a fit instruction (or None) for each model. If the fit instruction is None all data will be fitted to the model data.\n"
    
    if design == False:
        if not len(data)==len(timepoints):
            return False,"\nPlease provide data that correspond to the length of your timepoints!\n"

    if not len(ModelName)==len(priors):
	return False,"\nPlease provide prior distributions for each model!\n"
    
    if not len(ModelName)==len(x0priors):
	return False,"\nPlease provide initial values for each model!\n"
    
    if not modelKernel>0.0:
	return False,"\nPlease provide a model Kernel larger than 0!\n"
    
    if (len(epsilon)>1 and  custom_distance==False):
        return False,"\nPlease provide a custom distance function when you specify more than 1 epsilon schedule!\n"

    ### check model specific properties (comparing with SBML model)
    if not source==None:
        import libsbml

        if not len(source)==len(ModelName): 
            return False,"\nPlease provide the same amount of model sources and model names!\n"
        
	reader=libsbml.SBMLReader()
	for mod in range(0,len(source)):
	    
	    document=reader.readSBML(source[mod])
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
                    
            
            if not len(priors[mod])==numParameters:
                return False,"\nThe number of given prior distributions for model "+ModelName[mod]+" is not correct!\n"
            
            if not len(x0priors[mod])==numSpecies:
                return False,"\nPlease provide an initial value for each species in model "+ModelName[mod]+"!\n"

	  # number species:  fit-data-variables



    ### checking further properties independent of the model
    sde=re.compile('SDE')
    ode=re.compile('ODE')
    gillespie=re.compile('Gillespie')

            
    for mod in range(0,len(ModelName)):

	string=integrationType[mod]
	if (not sde.search(string)) and (not ode.search(string)) and (not gillespie.search(string)):
	    return False,"\nThe integration type for model "+ModelName[mod]+" does not exist!\n"

        for ic in range(len(x0priors[mod])):
             if not (x0priors[mod][ic][0]==0 or x0priors[mod][ic][0]==1 or x0priors[mod][ic][0]==2 or x0priors[mod][ic][0]==3):
		return False, "\nThe prior distribution of initial condition "+repr(ic+1)+" in model "+ModelName[mod]+" does not exist!\n"           

	for param in range(0,len(priors[mod])):
	    if not len(priors[mod][param])==3:
		return False, "\nThe prior distribution of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is wrong defined!\n"
	    
	    if not (priors[mod][param][0]==0 or priors[mod][param][0]==1 or priors[mod][param][0]==2 or priors[mod][param][0]==3):
		return False, "\nThe prior distribution of parameter "+repr(param+1)+" in model "+ModelName[mod]+" does not exist!\n"
            
	    if priors[mod][param][0]==2:
		if not priors[mod][param][1]<priors[mod][param][2]:
		    return False, "\nThe range of the uniform prior distribution of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is wrong defined!\n"
	    
	    if priors[mod][param][0]==3:
		if not (priors[mod][param][1]>=0 or priors[mod][param][2]>=0):
		    return False, "\nThe mean or scale of the lognormal prior distribution of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is wrong defined!\n"
	    
	    ##if kernel[mod][param][0]==1 and not priors[mod][param][0]==0:
            ##    if constKernels==False:
            ##        if not kernel[mod][param][1]<kernel[mod][param][2]:
            ##            return False, "\nThe range of the uniform pertubation kernel of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is wrong defined!\n"
	    
	    ##if kernel[mod][param][0]==2 and not priors[mod][param][0]==0:
            ##    if constKernels==False:
            ##        if not kernel[mod][param][2]>0:
            ##            return False, "\nThe variance of the gaussian pertubation kernel of parameter "+repr(param+1)+" in model "+ModelName[mod]+" is wrong defined!\n"


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

##def checkInputSimulation(name, timepoints, InitValues, IntegrationType, ConstantParameters, source, dt):
def checkInputSimulation(info_new , fname):
    """
    Check that the information in the input file is consistent with each other and with the model,
    and that it is in the format required to simulate the model.
    Return boolean, string (empty if boolean is True)
    
    """

    name = 	info_new.name
    timepoints = 	info_new.times
    IntegrationType =   info_new.type
    modelWeights = 	info_new.modelprior
    priors = 	        info_new.prior
    x0priors =          info_new.x0prior
    source = 	        info_new.source
    fit = 		info_new.fit
    beta = 		info_new.beta
    dt = 		info_new.dt

    for i in range(0, len(name)):
        if name[i] == "":
            return False, "\nPlease do not give empty strings for model names!\n"

    if not len(name)==len(IntegrationType):
        return False, "\nPlease provide the same number of model names and and integration types!\n"

    ### check model specific properties (comparing with SBML model)
    if not source==None:
        import libsbml

        if not len(source)==len(name): 
            return False,"\nPlease provide the same amount of model sources and model names!\n"
        
	reader=libsbml.SBMLReader()
	for mod in range(0,len(source)):
	    
	    document=reader.readSBML(source[mod])
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
                return False,"\nThe number of given parameters for model "+name[mod]+" is not correct!\n"
            
            if not len(info_new.x0prior[mod]) == numSpecies:
                return False,"\nPlease provide an initial value for each species in model "+name[mod]+"!\n"

    if len(timepoints)==0:
        return False, "\nPlease give timepoints at which to return simulated data points!\n"

    sde=re.compile('SDE')
    ode=re.compile('ODE')
    gillespie=re.compile('Gillespie')
            
    for mod in range(0,len(name)):

	string=IntegrationType[mod]
	if (not sde.search(string)) and (not ode.search(string)) and (not gillespie.search(string)):
	    return False,"\nThe integration type for model "+ModelName[mod]+" does not exist!\n"


    return True,""
