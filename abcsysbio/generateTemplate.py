

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



def generateTemplate(source,filename="input_file_template.txt", sumname="model_summary.txt"):

    """

    Generate a model summary file (model_summary.txt) and a template file (filename) from one or more SBML source files.
    
    
    ***** args *****
    
    source:    a tuple of strings.
               Each tuple entry describes a SBML file. 


    ***** kwargs *****
    
    filename:  a string.
               The name of the template to be generated.

    """

    out_file=open(filename,"w")
    sum_file=open(sumname,"w")

    out_file.write("#Number of models for which details are described in this input file\n\n")
    out_file.write("<modelnumber> "+repr(len(source))+"\n\n")
    out_file.write("#Probability of perturbing the sampled model (ignored when modelnumber = 1)\n\n")
    out_file.write("<model kernel> 0.7 \n\n")
    out_file.write("#Restart from previous (pickled) population?\n\n")
    out_file.write("<restart> False\n\n")
    out_file.write("#Series of epsilons. (Comma-delimited list)\n#If the length of the epsilon series is one and you have only one model\n#you are effectively doing abc Rejection\n\n")                   
    out_file.write("<epsilon> 1.0,0.5,0.2\n\n")
    out_file.write("#Population size\n\n")
    out_file.write("<population size> 1000\n\n")
    out_file.write("#Beta is the number of times to simulate each sampled parameter set.\n#This is only applicable for models simulated using SDE.\n\n")
    out_file.write("<beta> 1\n\n")
    out_file.write("#Internal timestep for solver.\n#Make this small for a stiff model.\n\n")
    out_file.write("<dt> 0.01\n\n")
    out_file.write("#The pertubation kernels are computed in respect to the previous parameter \n#distribution if constant parameter is set 'False'. \n#In this case you only need to provide the type of the pertubation kernel. \n#Otherwise you have to define the whole pertubation kernels.\n\n")
    out_file.write("<constant kernels> False\n\n")
    out_file.write("#rtol and atol can be specified here.\n#If the model is stiff then setting these to small\n#might help the simulation to run.\n#Only applicable for models simulated using ODE.\n\n")
    out_file.write("#<rtol>\n#<atol>\n\n")
    out_file.write("#User-supplied data.\n\n")
    out_file.write("<data>\n\n")
    out_file.write("#times: For abc-SMC, times must be a comma-delimited list starting with 0.\n#For simulation only the first and last timepoints are used.\n#To make a synthetic data set give a comma-delimited list of timepoints at which data points are required.\n\n")
    out_file.write("times:0,1,2,3,4,5,6,7,8,9,10\n\n")
    out_file.write("#variables: For abc-SMC, comma-delimited lists of concentrations (ODE or SDE) or molecule numbers (Gillespie).\n#Denote your data as variable1, variable2, ..., variableN.\n#For simulation or synthetic data sets these data are ignored.\n#See fitting instruction below if the dimensionality of your data sets differ from the dimensionality of your model.\n\n") 
    out_file.write("variable1:\n\n//\n\n")


    import libsbml
    reader=libsbml.SBMLReader()
        
      
    for i in range(0,len(source)):
        sum_file.write("Model "+repr(i+1)+"\n")
        out_file.write("<model"+repr(i+1)+">\n")
        sum_file.write("name: model"+repr(i+1)+"\nsource: "+source[i]+"\n\n")
        out_file.write("name: model"+repr(i+1)+"\nsource: "+source[i]+"\n\n")
        out_file.write("#type: the method used to simulate your model. ODE, SDE or Gillespie.\n\n")
        out_file.write("type: ODE\n\n")

        document=reader.readSBML(source[i])
        model=document.getModel()
        
        numSpecies=model.getNumSpecies()
        numGlobalParameters=model.getNumParameters()    
        
        parameter=[]
        parameterId=[]
        parameterId2=[]
        listOfParameter=[]
        
        r1=0
        r2=0
        r3=0
        listOfRules=model.getListOfRules()
        for k in range(0, len(listOfRules)):
            if model.getRule(k).isAlgebraic(): r1=r1+1
            if model.getRule(k).isAssignment(): r2=r2+1
            if model.getRule(k).isRate(): r3=r3+1

        comp=0
        NumCompartments=model.getNumCompartments()   
        for k in range(0,NumCompartments):
            if model.getCompartment(k).isSetVolume():
                comp=comp+1
                numGlobalParameters=numGlobalParameters+1
                parameter.append(model.getListOfCompartments()[k].getVolume())
                parameterId.append(model.getListOfCompartments()[k].getId())
                parameterId2.append('compartment'+repr(k+1))
                listOfParameter.append(model.getListOfCompartments()[k])
        
        for k in range(0,numGlobalParameters-comp):
            param=model.getParameter(k)
            parameter.append(param.getValue())
            parameterId.append(param.getId())
            parameterId2.append('parameter'+repr(k+1))
            listOfParameter.append(param)
      
        numLocalParameters=0
        NumReactions=model.getNumReactions()
        for k in range(0,NumReactions):
            local=model.getReaction(k).getKineticLaw().getNumParameters()
            numLocalParameters=numLocalParameters+local

            for j in range(0,local):
                parameter.append(model.getListOfReactions()[k].getKineticLaw().getParameter(j).getValue()) 
                parameterId.append(model.getListOfReactions()[k].getKineticLaw().getParameter(j).getId())
                x=len(parameterId)-comp
                parameterId2.append('parameter'+repr(x))
                listOfParameter.append(model.getListOfReactions()[k].getKineticLaw().getParameter(j))

        numParameters=numLocalParameters+numGlobalParameters
        
        species = model.getListOfSpecies()
        for k in range(0, len(species)):
            if (species[k].getConstant() == True):
                numParameters=numParameters+1
                parameter.append(getSpeciesValue(species[k]))
                parameterId.append(species[k].getId())
                parameterId2.append('species'+repr(k+1))
                numSpecies=numSpecies-1

        sum_file.write("number of compartments: "+repr(NumCompartments)+"\n")
        sum_file.write("number of reactions: "+repr(NumReactions)+"\n")
        sum_file.write("number of rules: "+repr(model.getNumRules())+"\n")
        if model.getNumRules()>0:
            sum_file.write("\t Algebraic rules: "+repr(r1)+"\n")
            sum_file.write("\t Assignment rules: "+repr(r2)+"\n")
            sum_file.write("\t Rate rules: "+repr(r3)+"\n\n")
        sum_file.write("number of functions: "+repr(model.getNumFunctionDefinitions())+"\n")
        sum_file.write("number of events: "+repr(model.getNumEvents())+"\n\n")
        

        paramAsSpecies=0
        sum_file.write("Species with initial values: "+repr(numSpecies)+"\n")
        out_file.write("#Initial values: Comma-delimited list.\n#Initial values are taken from the SBML model and are in the same order as the species\n#are given in the model.\n\n")
        out_file.write("initial values: ")
        x=0
        for k in range(0,len(species)):
            if (species[k].getConstant() == False):
                x=x+1
                out_file.write(repr(getSpeciesValue(species[k]))+", ")
                sum_file.write("S"+repr(x)+":\t"+species[k].getId()+"\tspecies"+repr(k+1)+"\t("+repr(getSpeciesValue(species[k]))+")\n")
        for k in range(0,len(listOfParameter)):
            if listOfParameter[k].getConstant()==False:
                for j in range(0, len(listOfRules)):
                    if listOfRules[j].isRate():
                        if parameterId[k]==listOfRules[j].getVariable():
                            x=x+1
                            paramAsSpecies=paramAsSpecies+1
                            out_file.write(repr(listOfParameter[k].getValue())+", ")
                            sum_file.write("S"+repr(x)+":\t"+listOfParameter[k].getId()+"\tparameter"+repr(k+1-comp)+"\t("+repr(listOfParameter[k].getValue())+") (parameter included in a rate rule and therefore treated as species)\n")


        out_file.write("\n")
        sum_file.write("\n")
        
        if(numGlobalParameters==0): string=" (all of them are local parameters)\n"
        elif(numGlobalParameters==1): string=" (the first parameter is a global parameter)\n"
        elif(numLocalParameters==0): string=" (all of them are global parameters)\n"
        else: string=" (the first "+repr(numGlobalParameters)+" are global parameter)\n"

        out_file.write("\n#Fitting information. If fit is None, all species in the model are fitted to the data in the order they are listed in the model.\n#Otherwise, give a comma-delimited list of fitting instrictions the same length as the dimensions of your data.\n#Use speciesN to denote the Nth species in your model. Simple arithmatic operations can be performed on the species from your model.\n#For example, to fit the sum of the first two species in your model to your first variable, write fit: species1+species2\n") 
        out_file.write("\nfit: None\n\n")

        out_file.write("#Parameters:\n#Priors:\n#one of \n#\tconstant, value\n#\tuniform, lower, upper\n#\tlognormal, location, scale\n")
        out_file.write("#Kernels: By default, gaussian with variance half the range of the prior distribution.\n#To force an alternative kernel use the switch in abcScript.py and give one of \n#\tuniform, lower, upper\n#\tgaussian, mean, variance\n\n")
        sum_file.write("Parameter: "+repr(numParameters)+string)
        sum_file.write("("+repr(paramAsSpecies)+" parameter is treated as species)\n")

        counter=0
        for k in range(0,numParameters-paramAsSpecies):
            Print = True
            if k<len(listOfParameter):
                if listOfParameter[k].getConstant()==False:
                    for j in range(0, len(listOfRules)):
                        if listOfRules[j].isRate():
                            if parameterId[k]==listOfRules[j].getVariable(): Print = False
            else: Print == True
            if Print ==True:
                counter=counter+1
                out_file.write("parameter"+repr(counter)+":\n")
                sum_file.write("P"+repr(counter)+":\t"+parameterId[k]+"\t"+parameterId2[k]+"\t("+repr(parameter[k])+")\n")
                out_file.write("\tprior: constant, ")
                out_file.write(repr(parameter[k])+"\n")
                out_file.write("\tkernel: uniform, -"+repr(parameter[k])+","+repr(parameter[k])+"\n")
        out_file.write("//\n\n")
        sum_file.write("\n############################################################\n\n")
            
    out_file.close()
    sum_file.close()

   

   
