

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



def generateTemplate(source, filename, sumname, dataname=None):

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
    


    out_file=open(filename,"w")
    sum_file=open(sumname,"w")

    have_data = False
    times = []
    vars = []
    nvar = 0
    first = True
    if dataname != None:
        have_data = True
        df = open(dataname,'r')
        for line in df:
            strs = str(line).split(' ')
            vals = [float(i) for i in strs]
            
            if first==True:
                for j in range(1,len(vals)):
                    vars.append([])
                first=False
                nvar = len(vals)-1

            times.append(vals[0])

            for j in range(1,len(vals)):
                vars[j-1].append(vals[j])

    #print times
    #print vars
        
    out_file.write("<input>\n\n")
    out_file.write("######################## number of models\n\n")
    out_file.write("# Number of models for which details are described in this input file\n")
    out_file.write("<modelnumber> "+repr(len(source))+ " </modelnumber>\n\n")

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
    out_file.write("<!--</epsilon> -->\n")
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
    out_file.write("# Currently uniform and normal are implemented\n")
    out_file.write("<kernel> uniform </kernel>\n\n")

    out_file.write("######################## model kernel : OPTIONAL (default 0.7)\n\n")
    out_file.write("# Probability of perturbing the sampled model (ignored when modelnumber = 1)\n")
    out_file.write("<modelkernel> 0.7 </modelkernel>\n\n")

    out_file.write("######################## ODE solver control : OPTIONAL \n\n")
    out_file.write("# rtol and atol can be specified here. If the model is stiff then setting these to small might help the simulation to run\n")
    out_file.write("#<rtol> </rtol> \n#<atol> </atol>\n\n")

    out_file.write("######################## User-supplied data\n\n")
    out_file.write("<data>\n")
    out_file.write("# times: For ABC SMC, times must be a whitespace delimited list\n")
    out_file.write("# In simulation mode these are the timepoints for which the simulations will be output\n")
    if have_data == False:
        out_file.write("<times> 0 1 2 3 4 5 6 7 8 9 10 </times>\n\n")
    else:
        out_file.write("<times>");
        for i in times:
            out_file.write(" "+repr(i) )
        out_file.write(" </times>\n\n");
    
    out_file.write("# variables: For ABC SMC, whitespace delimited lists of concentrations (ODE or SDE) or molecule numbers (Gillespie)\n")
    out_file.write("# Denote your data via tags <v1> </v1> or <var1> </var1> or <v2> </v2> etc. The tags are ignored and the data read in order\n")
    out_file.write("# For simulation these data are ignored\n")
    out_file.write("# See fitting instruction below if the dimensionality of your data sets differ from the dimensionality of your model\n")
    out_file.write("<variables>\n")
    if have_data == False:
        out_file.write(" <var1> </var1>\n")
    else:
        for k in range(nvar):
            out_file.write("<var"+repr(k+1)+"> ");
            for i in vars[k]:
                out_file.write(" "+repr(i) )
            out_file.write(" </var"+repr(k+1)+">\n");

    out_file.write("</variables>\n")
    out_file.write("</data>\n\n")

    out_file.write("######################## Models\n\n")
    out_file.write("<models>\n")

    import libsbml
    reader=libsbml.SBMLReader()
        
      
    for i in range(0,len(source)):
        sum_file.write("Model "+repr(i+1)+"\n")
        sum_file.write("name: model"+repr(i+1)+"\nsource: "+source[i]+"\n\n")

        out_file.write("<model"+repr(i+1)+">\n")
        out_file.write("<name> model"+repr(i+1)+" </name>\n<source> "+source[i]+" </source>\n\n")
        out_file.write("# type: the method used to simulate your model. ODE, SDE or Gillespie.\n")
        out_file.write("<type> SDE </type>\n\n")

        out_file.write("# Fitting information. If fit is None, all species in the model are fitted to the data in the order they are listed in the model.\n")
        out_file.write("# Otherwise, give a whitespace delimited list of fitting instrictions the same length as the dimensions of your data.\n")
        out_file.write("# Use speciesN to denote the Nth species in your model. Simple arithmetic operations can be performed on the species from your model.\n")
        out_file.write("# For example, to fit the sum of the first two species in your model to your first variable, write fit: species1+species2\n")
        out_file.write("<fit> None </fit>\n\n")

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
        ##for k in range(0, len(species)):
            ##if (species[k].getConstant() == True):
                ##numParameters=numParameters+1
                ##parameter.append(getSpeciesValue(species[k]))
                ##parameterId.append(species[k].getId())
                ##parameterId2.append('species'+repr(k+1))
                ##numSpecies=numSpecies-1

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

        out_file.write("# Priors on initial conditions and parameters:\n")
        out_file.write("# one of \n")
        out_file.write("#       constant, value \n")
        out_file.write("#       normal, mean, variance \n")
        out_file.write("#       uniform, lower, upper \n")
        out_file.write("#       lognormal, mean, variance \n\n")

        out_file.write("<initial>\n")

        x=0
        for k in range(0,len(species)):
            ##if (species[k].getConstant() == False):
            x=x+1
            #out_file.write(repr(getSpeciesValue(species[k]))+", ")
            out_file.write(" <ic"+repr(x)+"> constant "+repr(getSpeciesValue(species[k]))+" </ic"+repr(x)+">\n")
            sum_file.write("S"+repr(x)+":\t"+species[k].getId()+"\tspecies"+repr(k+1)+"\t("+repr(getSpeciesValue(species[k]))+")\n")
        for k in range(0,len(listOfParameter)):
            if listOfParameter[k].getConstant()==False:
                for j in range(0, len(listOfRules)):
                    if listOfRules[j].isRate():
                        if parameterId[k]==listOfRules[j].getVariable():
                            x=x+1
                            paramAsSpecies=paramAsSpecies+1
                            #out_file.write(repr(listOfParameter[k].getValue())+", ")
                            out_file.write(" <ic"+repr(x)+"> constant "+repr(listOfParameter[k].getValue())+" </ic"+repr(x)+">\n")
                            sum_file.write("S"+repr(x)+":\t"+listOfParameter[k].getId()+"\tparameter"+repr(k+1-comp)+"\t("+repr(listOfParameter[k].getValue())+") (parameter included in a rate rule and therefore treated as species)\n")

        out_file.write("</initial>\n\n")

        sum_file.write("\n")
        
        if(numGlobalParameters==0): string=" (all of them are local parameters)\n"
        elif(numGlobalParameters==1): string=" (the first parameter is a global parameter)\n"
        elif(numLocalParameters==0): string=" (all of them are global parameters)\n"
        else: string=" (the first "+repr(numGlobalParameters)+" are global parameter)\n"

        sum_file.write("Parameter: "+repr(numParameters)+string)
        sum_file.write("("+repr(paramAsSpecies)+" parameter is treated as species)\n")

        out_file.write("<parameters>\n")

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
                sum_file.write("P"+repr(counter)+":\t"+parameterId[k]+"\t"+parameterId2[k]+"\t("+repr(parameter[k])+")\n")
                out_file.write("<parameter"+repr(counter)+">")
                out_file.write(" constant ")
                out_file.write(repr(parameter[k])+" </parameter"+repr(counter)+">\n")
    
        sum_file.write("\n############################################################\n\n")

        out_file.write("</parameters>\n")
        out_file.write("</model"+repr(i+1)+">\n\n")

    out_file.write("</models>\n\n")
    out_file.write("</input>\n\n")

    out_file.close()
    sum_file.close()

   

   
