import sys

def getSpeciesValue(species):
    """
    Return the initial amount of a species.
    If species.isSetInitialAmount() == True, return the initial amount.
    Otherwise, return the initial concentration.

    ***** args *****
    
    species:    a libsbml.Species object
    
    """

    if species.isSetInitialAmount():
        return species.getInitialAmount()
    else:
        return species.getInitialConcentration()



def generateTemplate(source, filename, sumname, length, dataname=None):

    """

    Generate a model summary file (model_summary.txt) and a template file (filename) from one or more SBML source files.
    
    
    ***** args *****
    
    source:    a list of strings.
               Each entry describes a SBML file. 


    ***** kwargs *****
    
    filename:  a string.
               The name of the template to be generated.

    sumnname:  a string.
               The name of the summary to be generated.

    dataname:  a string.
               The name of a datafile.

    """

   # get source list in correct order to keep submodel sets together in the correct models

    nmodels = len(source)/length
    
    sourceList = []
    for i in range(0,nmodels):
        for j in range(0,length):
            sourceList.append(source[j*nmodels+i])
    print sourceList

    source = sourceList

    # name files to which data is being written
    out_file=open(filename,"w")
    sum_file=open(sumname,"w")


    # ################ READ IN DATA FROM TEXT FILES ACCORDING TO COMMAND-LINE INPUT (LIST OF STRINGS) ############################

    have_data = False
    times = [[]]
    vars = [[]]
    nvar = [0 for i in range(0,length)]

    first = True
    if dataname != None:
        have_data = True
        for group in range(0,length):        
            df = open(dataname[group],'r')         
            if first != True:
                times.append([])
                vars.append([])
                first = True
            for line in df:
                strs = str(line).split(' ')
                vals = [float(i) for i in strs]
            
                if first==True:
                    for j in range(1,len(vals)):
                        vars[group].append([])		
                    first=False
                    nvar[group] = len(vals)-1		
                
                times[group].append(vals[0])		
                for j in range(1,len(vals)):
                    vars[group][j-1].append(vals[j])	


    # write output file to be used as input for run-abc-sysbio
    out_file.write("<input>\n\n")
    out_file.write("######################## number of models\n\n")
    out_file.write("# Number of models for which details are described in this input file\n")
    out_file.write("<modelnumber> "+repr(len(source)/length)+ " </modelnumber>\n\n")
    out_file.write("<submodelnumber> "+repr(length)+ " </submodelnumber>\n\n")

    out_file.write("######################## restart\n\n")
    out_file.write("# Restart from previous (pickled) population?\n")
    out_file.write("<restart> False </restart>\n\n")

    out_file.write("######################## epsilon schedule\n\n")
    out_file.write("# Automatic epsilon schedule. Provide a vector of final epsilons and the alpha (defaults to 0.9)\n")
    out_file.write("<autoepsilon>\n") 
    out_file.write("<finalepsilon> 1.0 </finalepsilon>\n")
    out_file.write("<alpha> 0.9 </alpha>\n")
    out_file.write("</autoepsilon>\n\n")

    out_file.write("# OR\n")
    out_file.write("# Series of epsilons. (Whitespace delimited list)\n")
    out_file.write("# Multiple epsilon schedules can be specified by giving additional vectors enclosed in <e2> </e2>, <e3> </e3> etc\n")
    out_file.write("# NOTE: the parser always reads them in order and ignores the tag value\n")
    out_file.write("<!-- <epsilon> -->\n")
    out_file.write("<!-- <e1> 5.0 3.0 2.0 1.0 </e1> -->\n")
    out_file.write("<!-- </epsilon> -->\n")
    out_file.write("\n")
        
    out_file.write("######################## particles\n\n")
    out_file.write("<particles> 100 </particles>\n\n")

    out_file.write("######################## beta\n\n")
    out_file.write("# Beta is the number of times to simulate each sampled parameter set.\n# This is only applicable for models simulated using Gillespie and SDE\n")
    out_file.write("<beta> 1 </beta>\n\n")

    out_file.write("######################## dt\n\n")
    out_file.write("# Internal timestep for solver.\n# Make this small for a stiff model.\n")
    out_file.write("<dt> 0.01 </dt>\n\n")

    out_file.write("######################## perturbation kernels : OPTIONAL (default uniform)\n\n")
    out_file.write("# The pertubation kernels are computed with respect to the previous parameter distribution\n")
    out_file.write("# Currently adaptive and non-adaptive uniform and normal perturbation kernels can be implemented\n")
    out_file.write("# To use adaptive kernels, use just <kernel> uniform </kernel> or <kernel> normal </kernel>\n")
    out_file.write("# To use non-adaptive kernels, use <kernel> uniform lowerbound upperbound </kernel> or <kernel> normal mean variance </kernel>\n")
    out_file.write("<kernel> uniform </kernel>\n\n")

    out_file.write("######################## model kernel : OPTIONAL (default 0.7)\n\n")
    out_file.write("# Probability of perturbing the sampled model (ignored when modelnumber = 1)\n")
    out_file.write("<modelkernel> 0.7 </modelkernel>\n\n")

    out_file.write("######################## ODE solver control : OPTIONAL \n\n")
    out_file.write("# rtol and atol can be specified here. If the model is stiff then setting these to small might help the simulation to run\n")
    out_file.write("#<rtol> </rtol> \n#<atol> </atol>\n\n")

    ########################################## DEFINE DATA SETS CORRESPONDING TO EACH SUBMODEL TYPE (E.G. WT/KO1/KO2 ETC.) ############################

    out_file.write("######################## User-supplied data\n\n")
    out_file.write("<data>\n\n")

    out_file.write("# times: For ABC SMC, times must be a whitespace delimited list\n")
    out_file.write("# In simulation mode these are the timepoints for which the simulations will be output\n\n")
    out_file.write("# variables: For ABC SMC, whitespace delimited lists of concentrations (ODE or SDE) or molecule numbers (Gillespie)\n")
    out_file.write("# Denote your data via tags <v1> </v1> or <var1> </var1> or <v2> </v2> etc. The tags are ignored and the data read in order\n")
    out_file.write("# For simulation these data are ignored\n")
    out_file.write("# See fitting instruction below if the dimensionality of your data sets differ from the dimensionality of your model\n\n")

    sourceCounter = 0

    if have_data == False:
        out_file.write("<times> 0 1 2 3 4 5 6 7 8 9 10 </times>\n\n")
        out_file.write(" <var1> </var1>\n")

    else:
        for i in range (0,len(dataname)):
            out_file.write("<data"+repr(i+1)+">\n")

            out_file.write("<times>");
            for k in times[i]:
                out_file.write(" "+repr(k) )
            out_file.write(" </times>\n\n");
            out_file.write("<variables>\n")

            for k in range(nvar[i]):
                out_file.write("<var"+repr(k+1)+"> ");
                for l in vars[i][k]:
                    out_file.write(" "+repr(l) )
                out_file.write(" </var"+repr(k+1)+">\n");

            out_file.write("</variables>\n")
            out_file.write("</data"+repr(i+1)+">\n\n\n")

    out_file.write("</data>\n\n")

    out_file.write("######################## Models with associated submodels\n\n")
    out_file.write("<models>\n")


##################################### Extract model from SBML files and store in arrays ################################################################

    import libsbml
    reader=libsbml.SBMLReader()

    parameter=[]
    parameterId=[]
    parameterId2=[]
    listOfParameter=[]
    species=[]
    numSpecies=[]
    numGlobalParameters=[]
    listOfRules=[]
    NumCompartments=[]
    NumReactions=[]
    numLocalParameters=[]
    numParameters=[]
    numRules=[]
    numFunctionDefinitions=[]
    numEvents=[]
    paramAsSpecies=[]
    string=[]
    counter=[]

    inival=[]
    speciesId=[]
    conOrVar=[]
    isRate=[]
    getVar=[]

    for i in range (0,len(source)):

        document=reader.readSBML(source[i])
        model=document.getModel()
        
        numSpeciesPerMod=model.getNumSpecies()
        numSpecies.append(numSpeciesPerMod)

        numGlobalParametersPerMod=model.getNumParameters()
        numGlobalParameters.append(numGlobalParametersPerMod)    
        
        parameterPerMod=[]
        parameterIdPerMod = []
        parameterId2PerMod=[]
        listOfParameterPerMod=[]

        inivalPerMod=[]
        speciesIdPerMod=[]
        conOrVarPerMod=[]
        isRatePerMod=[]
        getVarPerMod=[]

        r1=0
        r2=0
        r3=0
        listOfRulesPerMod=model.getListOfRules()
        listOfRules.append(listOfRulesPerMod)
        for k in range(0, len(listOfRules[i])):
                if model.getRule(k).isAlgebraic(): r1=r1+1
                if model.getRule(k).isAssignment(): r2=r2+1
                if model.getRule(k).isRate(): r3=r3+1

        comp=0
        NumCompartmentsPerMod=model.getNumCompartments()
        NumCompartments.append(NumCompartmentsPerMod)   
        for k in range(0,NumCompartments[i]):
                if model.getCompartment(k).isSetVolume():
                        comp=comp+1
                        numGlobalParameters[i]=numGlobalParameters[i]+1
                        parameterPerMod.append(model.getListOfCompartments()[k].getVolume())
                        parameterIdPerMod.append(model.getListOfCompartments()[k].getId())
                        parameterId2PerMod.append('compartment'+repr(k+1))
                        listOfParameterPerMod.append(model.getListOfCompartments()[k])
        
        for k in range(0,numGlobalParameters[i]-comp):
                param=model.getParameter(k)
                parameterPerMod.append(param.getValue())
                parameterIdPerMod.append(param.getId())
                parameterId2PerMod.append('parameter'+repr(k+1))
                listOfParameterPerMod.append(param)
      
        numLocalParametersPerMod=0
        NumReactionsPerMod=model.getNumReactions()
        NumReactions.append(NumReactionsPerMod)
        for k in range(0,NumReactions[i]):
                local=model.getReaction(k).getKineticLaw().getNumParameters()
                numLocalParametersPerMod=numLocalParametersPerMod+local
                numLocalParameters.append(numLocalParametersPerMod)

        for j in range(0,local):
                parameterPerMod.append(model.getListOfReactions()[k].getKineticLaw().getParameter(j).getValue()) 
                parameterIdPerMod.append(model.getListOfReactions()[k].getKineticLaw().getParameter(j).getId())
                x=len(parameterId)-comp
                parameterId2PerMod.append('parameter'+repr(x))
                listOfParameterPerMod.append(model.getListOfReactions()[k].getKineticLaw().getParameter(j))

        parameter.append(parameterPerMod)
        parameterId.append(parameterIdPerMod)
        parameterId2.append(parameterId2PerMod)
        listOfParameter.append(listOfParameterPerMod)

        numParametersPerMod=numLocalParameters[i]+numGlobalParameters[i]
        numParameters.append(numParametersPerMod)
        
        numRulesPerMod=model.getNumRules()
        numRules.append(numRulesPerMod)

        numFunctionDefinitionsPerMod=model.getNumFunctionDefinitions()
        numFunctionDefinitions.append(numFunctionDefinitionsPerMod)

        numEventsPerMod=model.getNumEvents()
        numEvents.append(numEventsPerMod)


# ----------------------------------- Defining initial conditions according to SBML files --------------------------------

        speciesPerMod = model.getListOfSpecies()
        species.append(speciesPerMod)

        paramAsSpeciesPerMod=0

        x=0
        for k in range(0,len(species[i])):
            x=x+1
            inivalPerMod.append(getSpeciesValue(species[i][k]))
            speciesIdPerMod.append(species[i][k].getId())

        inival.append(inivalPerMod)
        speciesId.append(speciesIdPerMod)

        for k in range(0,len(listOfParameter[i])):
             conOrVarPerMod.append(listOfParameter[i][k].getConstant())

        conOrVar.append(conOrVarPerMod)
        
        for j in range(0, len(listOfRules[i])):
             isRatePerMod.append(listOfRules[i][k].isRate())
             getVarPerMod.append(listOfRules[i][j].getVariable())

        isRate.append(isRatePerMod)
        getVar.append(getVarPerMod)

        for k in range(0,len(listOfParameter[i])):
            if conOrVar[i][k]==False:
                for j in range(0, len(listOfRules[i])):
                    if isRate:
                        if parameterId[i][k]==getVar:
                            x=x+1
                            paramAsSpeciesPerMod=paramAsSpeciesPerMod+1

        paramAsSpecies.append(paramAsSpeciesPerMod)



################################# Info for models and submodels ###################################################################################

    # Print model information to the template file and summary file
    counterPerMod=0
    for i in range(0,nmodels):
        out_file.write("######################## Model (submodels listed in section below)\n\n")
        sum_file.write("Model "+repr(i+1)+"\n\n")


        out_file.write("<model"+repr(i+1)+">\n")
        sum_file.write("<model"+repr(i+1)+">\n")
        out_file.write("<name> model"+repr(i+1)+" </name>\n")
        sum_file.write("<name> model"+repr(i+1)+" </name>\n")
        out_file.write("# type: the method used to simulate your model. ODE, SDE or Gillespie.\n")
        out_file.write("<type> SDE </type>\n\n")

        out_file.write("# Priors on initial conditions and parameters:\n")
        out_file.write("# one of \n")
        out_file.write("#       constant, value \n")
        out_file.write("#       normal, mean, variance \n")
        out_file.write("#       uniform, lower, upper \n")
        out_file.write("#       lognormal, mean, variance \n")
        out_file.write("#       trunc_normal, mean, variance, lowerbound, upperbound (either lowerbound or upperbound can be 'i')\n\n")

        if(numGlobalParameters[i]==0): stringPerMod=" (all of them are local parameters)\n"
        elif(numGlobalParameters[i]==1): stringPerMod=" (the first parameter is a global parameter)\n"
        elif(numLocalParameters[i]==0): stringPerMod=" (all of them are global parameters)\n"
        else: stringPerMod=" (the first "+repr(numGlobalParameters[i])+" are global parameter)\n"
        string.append(stringPerMod)

        sum_file.write("Parameter: "+repr(numParameters[i])+string[i])
        sum_file.write("("+repr(paramAsSpecies[i])+" parameter is treated as species)\n")

        out_file.write("<parameters>\n")


        paramSamplingPerMod=0
        paramSamplingPerMod=counterPerMod*length
        n=paramSamplingPerMod

        paramCounter=1
        for k in range(0,numParameters[n]-paramAsSpecies[n]):
            Print = True
            if k<len(listOfParameter[n]):
                if conOrVar[n][k]==False:
                    for j in range(0, len(listOfRules[n])):
                        if isRate:
                            if parameterId[n][k]==getVar: Print = False
            else: Print == True
            if Print ==True:
                sum_file.write("P"+repr(paramCounter)+":\t"+parameterId[n][k]+"\t"+parameterId2[n][k]+"\t("+repr(parameter[n][k])+")\n")
                out_file.write("<parameter"+repr(paramCounter)+">")
                out_file.write(" constant ")
                out_file.write(repr(parameter[n][k])+" </parameter"+repr(paramCounter)+">\n")
                paramCounter=paramCounter+1
        counterPerMod=counterPerMod+1
        counter.append(counterPerMod)    


        sum_file.write("\n############################################################\n\n")

        out_file.write("</parameters>\n\n")


# ---------------------------------------------------------------- Submodels ------------------------------------------------------------------------
        
        # Print submodel information to the input file and summary file within the model section (i.e. submodels from a given model)
        out_file.write("######################## Set of associated submodels\n\n")
        out_file.write("<submodels>\n")

        for j in range(0,length):
                out_file.write("<submodel"+repr(i+1)+"_"+repr(j+1)+">\n")
                sum_file.write("<submodel> "+repr(i+1)+"_"+repr(j+1)+" </submodel>\n")

                out_file.write("<name> model"+repr(i+1)+"_"+repr(j+1)+" </name>\n<source> "+source[j*nmodels+i]+" </source>\n\n")
                sum_file.write("<name> model"+repr(i+1)+"_"+repr(j+1)+" </name>\n<source> "+source[j*nmodels+i]+" </source>\n\n")



                out_file.write("# Fitting information. If fit is None, all species in the model are fitted to the data in the order they are listed in the model.\n")
                out_file.write("# Otherwise, give a whitespace delimited list of fitting instructions the same length as the dimensions of your data.\n")
                out_file.write("# Use speciesN to denote the Nth species in your model. Simple arithmetic operations can be performed on the species from your model.\n")
                out_file.write("# For example, to fit the sum of the first two species in your model to your first variable, write fit: species1+species2\n")
                out_file.write("<fit> None </fit>\n\n")


                out_file.write("# Optional argument. Distance function can be euclidian, log, or custom (to be defined in customABC.py). The default is euclidian, unless 'custd' is given in the command line, in which case custom is the default.\n")
                out_file.write("<!-- <distancefunction>    </distancefunction> -->\n\n")

                out_file.write("# If any of the parameters must be considered as local parameters (i.e. sampled using different distributions for different submodels), specify these below.\n")
                out_file.write("# Ensure that all local parameters are subsequently deleted from the overall model parameter list (above).\n")
                out_file.write("# Local parameters can be specified by giving additional vectors enclosed in <param2> </param2>, <param3> </param3> etc\n")
                out_file.write("# NOTE: the parser always reads them in order and ignores the tag value\n")
                out_file.write("<!-- <localparameters> -->\n")
                out_file.write("<!-- <param1> constant 1.0 </param1> -->\n")
                out_file.write("<!-- <param2> uniform 1.0 10.0 </param2> -->\n")
                out_file.write("<!-- </localparameters> -->\n\n")

                out_file.write("# Optional argument. Distance function can be euclidian, log, or custom (to be defined in customABC.py).\n")
                out_file.write("# The default is euclidian, unless 'custd' is given in the command line, in which case custom is the default.\n")
                out_file.write("<!-- <distancefunction>    </distancefunction> -->\n\n")


# ------------------------------------ Extract initial conditions for each submodel using the arrays ------------------------------------------------

                sum_file.write("number of compartments: "+repr(NumCompartments[i])+"\n")
                sum_file.write("number of reactions: "+repr(NumReactions[i])+"\n")
                sum_file.write("number of rules: "+repr(numRules[i])+"\n")
                if numRules[i]>0:
                        sum_file.write("\t Algebraic rules: "+repr(r1)+"\n")
                        sum_file.write("\t Assignment rules: "+repr(r2)+"\n")
                        sum_file.write("\t Rate rules: "+repr(r3)+"\n\n")
                sum_file.write("number of functions: "+repr(numFunctionDefinitions[i])+"\n")
                sum_file.write("number of events: "+repr(numEvents[i])+"\n\n")

                
                sum_file.write("Species with initial values: "+repr(numSpecies[i])+"\n")


                out_file.write("<initial>\n")

                n=sourceCounter

                x=0
                for k in range(0,len(inival[n])):
                    x=x+1
                    out_file.write(" <ic"+repr(x)+"> constant "+repr(inival[n][k])+" </ic"+repr(x)+">\n")
                    sum_file.write("S"+repr(x)+":\t"+repr(speciesId[n][k])+"\tspecies"+repr(k+1)+"\t("+repr(inival[n][k])+")\n")
                for k in range(0,len(listOfParameter[n])):
                    if conOrVar[n][k]==False:
                        for j in range(0, len(listOfRules[n])):
                            if isRate:
                                if parameterId[n][k]==getVar:
                                    x=x+1
                                    out_file.write(" <ic"+repr(x)+"> constant "+repr(listOfParameter[n][k].getValue())+" </ic"+repr(x)+">\n")
                                    sum_file.write("S"+repr(x)+":\t"+listOfParameter[n][k].getId()+"\tparameter"+repr(k+1-comp)+"\t("+repr(listOfParameter[n][k].getValue())+") (parameter included in a rate rule and therefore treated as species)\n")


                sourceCounter = sourceCounter + 1

                out_file.write("</initial>\n")
                sum_file.write("\n")
                sum_file.write("############################################################\n")
                sum_file.write("\n")

                out_file.write("</submodel"+repr(i+1)+"_"+repr(j+1)+">\n\n")
        out_file.write("</submodels>\n\n")

        out_file.write("</model"+repr(i+1)+">\n\n")

    # Close tags and close files:
    out_file.write("</models>\n\n")
    out_file.write("</input>\n\n")

    out_file.close()
    sum_file.close()




   
