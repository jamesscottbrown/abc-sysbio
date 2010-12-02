from numpy import *
from libsbml import *
import re
import os
from abcsysbio.relations import *

##To call the parser:
##    SBMLparse.importSBML(source, integrationType, ModelName=None,method=None)
##    All arguments to function must be passed as tuples.
##    If there is only one source to parse it must still be passed as a tuple ('source.xml',)
##    with an integrationType passed as ('Gillespie',)

def write_GillespieFunctions(stoichiometricMatrix, kineticLaw, numSpecies, numReactions, species, parameterId, InitValues, speciesId,name, listOfFunctions, FunctionArgument, FunctionBody, listOfRules, ruleFormula, ruleVariable, listOfEvents, EventCondition, EventVariable, EventFormula):
    """
    Write the python files with the arguments to the Gillespie algorithm taken from the parser
    """

    for i in range(0,len(listOfRules)):
        if listOfRules[i].isRate():
            print "\n Model '"+name+"' contains at least one rate rule. This model can not be parsed and simmulated with the Gillespie algorithm! Please change the simmulation Type! \n"
            sys.exit()

    p=re.compile('\s')
    out_file=open(name+".py","w")

    out_file.write("from abcsysbio.relations import *\n\n#Functions\n")

    for i in range(0, len(listOfFunctions)):
        out_file.write("def ")
        out_file.write(listOfFunctions[i].getId())
        out_file.write("(")
        for j in range(0, listOfFunctions[i].getNumArguments()):
            out_file.write(FunctionArgument[i][j])
            out_file.write(",")
        out_file.write("):\n\n\toutput=")
        out_file.write(FunctionBody[i])
        out_file.write("\n\n\treturn output\n\n") 

    out_file.write("\n#Gillespie Hazards\n\n")

    out_file.write("def Hazards((")

    for i in range(0,numSpecies):
        if (species[i].getConstant() == False):
            out_file.write(speciesId[i])
            out_file.write(",")
  
    out_file.write("),parameter):\n\n")  
    
    for i in range(0,len(parameterId)):
        out_file.write("\t"+parameterId[i]+"=parameter["+repr(i)+"]\n")

    counter = len(parameterId)
    for i in range(0, numSpecies):
        if (species[i].getConstant() == True):
            out_file.write("\t"+speciesId[i]+"=parameter["+repr(counter)+"]\n")
            counter = counter+1
    
    out_file.write("\n")

    for i in range(0,numReactions):
        out_file.write("\tHazard_"+repr(i)+" = "+kineticLaw[i])
        out_file.write("\n")
                        
    out_file.write("\n\treturn(")
    
    for i in range(0,numReactions):
        out_file.write("Hazard_"+repr(i))
        if(not i==(numReactions-1)):
            out_file.write(", ")

    out_file.write(")\n\n")

    out_file.write("#Gillespie Reactions\n\n")

    for i in range(0, numReactions):
        out_file.write("def Reaction"+repr(i)+"((")
        for k in range(0,numSpecies):
            if (species[k].getConstant() == False):
                out_file.write(speciesId[k])
                out_file.write(",")

        out_file.write(")):\n\n")

        for k in range(0,numSpecies):
            if (species[k].getConstant() == False):
                out_file.write("\t"+speciesId[k]+"_new="+speciesId[k]+"+("+str(stoichiometricMatrix[k][i])+")\n")

        out_file.write("\n\treturn(")
        for k in range(0,numSpecies):
            if (species[k].getConstant() == False):
                out_file.write(speciesId[k]+"_new")
                out_file.write(",")
        out_file.write(")\n\n")
        
    out_file.write("#Dictionary of reactions\ndef defaultfunc():\n\tpass\n\ndef Switch():\n\tswitch = {\n")
    for i in range(0, numReactions):
        out_file.write("\t\t"+repr(i)+" : Reaction"+repr(i)+",\n")

    out_file.write("\t\t\"default\": defaultfunc\n\t\t}\n\treturn switch\n\n")

    out_file.write("#Rules and Events\n")

    out_file.write("def rules((")
    for i in range(0,numSpecies):
            if (species[i].getConstant() == False):
                out_file.write(speciesId[i])
                out_file.write(",")
    out_file.write("),(")
    for i in range(0,len(parameterId)):
        out_file.write(parameterId[i])
        out_file.write(',')
    out_file.write("),t):\n\n")

    for i in range(0, len(listOfRules)):
        if listOfRules[i].isAssignment():
            out_file.write("\t")
            out_file.write(ruleVariable[i])
            out_file.write("=")
            out_file.write(ruleFormula[i])
            out_file.write("\n")

    out_file.write("\n\treturn((")
    for i in range(0, numSpecies):
        if (species[i].getConstant() == False):
            out_file.write(speciesId[i])
            out_file.write(",")

    out_file.write("),(")
    for i in range(0,len(parameterId)):
        out_file.write(parameterId[i])
        out_file.write(',')
    out_file.write("))\n\n")

    out_file.write("def events((")
    for i in range(0,numSpecies):
            if (species[i].getConstant() == False):
                out_file.write(speciesId[i])
                out_file.write(",")
    out_file.write("),(")
    for i in range(0,len(parameterId)):
        out_file.write(parameterId[i])
        out_file.write(',')
    out_file.write("),t):\n\n")

    for i in range(0, len(listOfEvents)):
        out_file.write("\tif ")
        out_file.write(mathMLConditionParser(EventCondition[i]))
        out_file.write(":\n")
        listOfAssignmentRules = listOfEvents[i].getListOfEventAssignments()
        for j in range(0, len(listOfAssignmentRules)):
            out_file.write("\t\t")
            out_file.write(EventVariable[i][j])
            out_file.write("=")
            out_file.write(EventFormula[i][j])
            out_file.write("\n")
        out_file.write("\n")

    out_file.write("\n\treturn((")
    for i in range(0, numSpecies):
        if (species[i].getConstant() == False):
            out_file.write(speciesId[i])
            out_file.write(",")

    out_file.write("),(")
    for i in range(0,len(parameterId)):
        out_file.write(parameterId[i])
        out_file.write(',')
    out_file.write("))\n\n")

                   
    out_file.close()


def write_ODEFunctions(stoichiometricMatrix, kineticLaw, species, numReactions,speciesId,listOfParameter, parameterId,parameter,InitValues,name, listOfFunctions, FunctionArgument, FunctionBody, listOfRules, ruleFormula, ruleVariable, listOfEvents, EventCondition, EventVariable, EventFormula):
    """
    Write the python file with ODE functions using the information taken by the parser
    """ 
    p=re.compile('\s')
    #Open the outfile
    out_file=open(name+".py","w")
    #Import the necessaries
    out_file.write("from math import *\nfrom numpy import *\nfrom abcsysbio.relations import *\n\n")
    #The user-defined functions used in the model must be written in the file

    for i in range(0, len(listOfFunctions)):
        out_file.write("def ")
        out_file.write(listOfFunctions[i].getId())
        out_file.write("(")
        for j in range(0, listOfFunctions[i].getNumArguments()):
            out_file.write(FunctionArgument[i][j])
            out_file.write(",")
        out_file.write("):\n\n\toutput=")
        out_file.write(FunctionBody[i])
        out_file.write("\n\n\treturn output\n\n")

    #Write the modelfunction
   
    out_file.write("def modelfunction((")

    numSpecies = len(species)
    
    for i in range(0, len(species)):
        if (species[i].getConstant() == False):
            out_file.write(speciesId[i])
            out_file.write(",")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write(parameterId[i])
                    out_file.write(",")

    out_file.write("),time,parameter=(")
    for i in range(0,len(parameterId)):
        dontPrint = False
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]): 
                    dontPrint=True
        if dontPrint == False:
            out_file.write(repr(parameter[i]))
            out_file.write(",")

    for i in range(0, numSpecies):
        if (species[i].getConstant() == True):
            out_file.write(repr(InitValues[i]))
            out_file.write(",")
    out_file.write(")):\n\n")  
    
    counter=0
    for i in range(0,len(parameterId)):
        dontPrint = False
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]): dontPrint=True
        if dontPrint == False:
            out_file.write("\t"+parameterId[i]+"=parameter["+repr(counter)+"]\n")
            counter=counter+1

    for i in range(0, numSpecies):
        if (species[i].getConstant() == True):
            out_file.write("\t"+speciesId[i]+"=parameter["+repr(counter)+"]\n")
            counter = counter+1
    
    out_file.write("\n")



    #write the derivatives

    for i in range(0,numSpecies):
        if (species[i].getConstant() == False):
            out_file.write("\td_"+speciesId[i]+"=")
            if (species[i].isSetCompartment() == True):
                out_file.write("(")
            for k in range(0,numReactions):
                if(not stoichiometricMatrix[i][k]==0.0):
                    out_file.write("(")
                    out_file.write(repr(stoichiometricMatrix[i][k]))
                    out_file.write(")*(")
                    string=p.sub('',kineticLaw[k])
                    out_file.write(string)
                    out_file.write(")+")
            out_file.write("0")
            if (species[i].isSetCompartment() == True):
                out_file.write(")/")
                mySpeciesCompartment = species[i].getCompartment()
                for j in range(0, len(listOfParameter)):
                    if (listOfParameter[j].getId() == mySpeciesCompartment):
                        out_file.write(parameterId[j])
                        break
            out_file.write("\n")

    for i in range(0,len(listOfRules)):
        if listOfRules[i].isRate() == True:
            out_file.write("\td_")
            out_file.write(ruleVariable[i])
            out_file.write("=")
            out_file.write(ruleFormula[i])
            out_file.write("\n")
            
            
                
    out_file.write("\n\treturn(")
  
    for i in range(0, len(species)):
        if (species[i].getConstant() == False):
            out_file.write("d_"+speciesId[i])
            out_file.write(",")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write("d_"+parameterId[i])
                    out_file.write(",")

    out_file.write(")\n")

    #Write the rules
    out_file.write("\ndef rules((")

    for i in range(0, len(species)):
        if (species[i].getConstant() == False):
            out_file.write(speciesId[i])
            out_file.write(",")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write(parameterId[i])
                    out_file.write(",")
    out_file.write("),(")
    for i in range(0,len(parameterId)):
        dontPrint = False
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]): 
                    dontPrint=True
        if dontPrint == False:
            out_file.write(parameterId[i])
            out_file.write(",")

    for i in range(0, numSpecies):
        if (species[i].getConstant() == True):
            out_file.write(speciesId[i])
            out_file.write(",")

    out_file.write("),time):\n\n")

    #Write the events

    for i in range(0, len(listOfEvents)):
        out_file.write("\tif ")
        out_file.write(mathMLConditionParser(EventCondition[i]))
        out_file.write(":\n")
        listOfAssignmentRules = listOfEvents[i].getListOfEventAssignments()
        for j in range(0, len(listOfAssignmentRules)):
            out_file.write("\t\t")
            out_file.write(EventVariable[i][j])
            out_file.write("=")
            out_file.write(EventFormula[i][j])
            out_file.write("\n")
        out_file.write("\n")

    out_file.write("\n")

    #write the rules

    for i in range(0, len(listOfRules)):
        if listOfRules[i].isAssignment():
            out_file.write("\t")
            out_file.write(ruleVariable[i])
            out_file.write("=")
            out_file.write(mathMLConditionParser(ruleFormula[i]))
            out_file.write("\n")


    out_file.write("\n\treturn((")
    for i in range(0, numSpecies):
        if (species[i].getConstant() == False):
            out_file.write(speciesId[i])
            out_file.write(",")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write(parameterId[i])
                    out_file.write(",")
    out_file.write("),(")

    for i in range(0,len(parameterId)):
        dontPrint = False
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]): 
                    dontPrint=True
        if dontPrint == False:
            out_file.write(parameterId[i])
            out_file.write(",")

    for i in range(0, numSpecies):
        if (species[i].getConstant() == True):
            out_file.write(speciesId[i])
            out_file.write(",")
    out_file.write("))\n\n")               
    out_file.close()

def write_SDEFunctions(stoichiometricMatrix, kineticLaw, species, numReactions,speciesId, listOfParameter,parameterId,parameter,InitValues,method,name, listOfFunctions, FunctionArgument, FunctionBody,listOfRules, ruleFormula, ruleVariable, listOfEvents, EventCondition, EventVariable, EventFormula):
    """
    Write the python file with SDE functions using the information taken by the parser
    """ 
    p=re.compile('\s')
    out_file=open(name+".py","w")
    out_file.write("from math import sqrt\nfrom numpy import random\nfrom abcsysbio.relations import *\n\n")

    for i in range(0, len(listOfFunctions)):
        out_file.write("def ")
        out_file.write(listOfFunctions[i].getId())
        out_file.write("(")
        for j in range(0, listOfFunctions[i].getNumArguments()):
            out_file.write(FunctionArgument[i][j])
            out_file.write(",")
        out_file.write("):\n\n\toutput=")
        out_file.write(FunctionBody[i])
        out_file.write("\n\n\treturn output\n\n")
        
    out_file.write("def modelfunction((")

    numSpecies = len(species)

    for i in range(0, len(species)):
        if (species[i].getConstant() == False):
            out_file.write(speciesId[i])
            out_file.write(",")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write(parameterId[i])
                    out_file.write(",")

    out_file.write(")=(")

    for i in range(0, len(species)):
        if (species[i].getConstant() == False):
            out_file.write(repr(InitValues[i]))
            out_file.write(",")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write(repr(parameter[i]))
                    out_file.write(",")

    out_file.write("),dt=0,parameter=(")

    for i in range(0,len(parameterId)):
        dontPrint = False
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]): 
                    dontPrint=True
        if dontPrint == False:
            out_file.write(repr(parameter[i]))
            out_file.write(",")

    for i in range(0, numSpecies):
        if (species[i].getConstant() == True):
            out_file.write(repr(InitValues[i]))
            out_file.write(",")
 
    out_file.write("),time=0):\n\n")

    counter=0
    for i in range(0,len(parameterId)):
        dontPrint = False
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]): dontPrint=True
        if dontPrint == False:
            out_file.write("\t"+parameterId[i]+"=parameter["+repr(counter)+"]\n")
            counter=counter+1

    for i in range(0, numSpecies):
        if (species[i].getConstant() == True):
            out_file.write("\t"+speciesId[i]+"=parameter["+repr(counter)+"]\n")
            counter = counter+1

    out_file.write("\n")
	

    
    out_file.write("\n")
    
    for i in range(0,numSpecies):
        if (species[i].getConstant() == False):
            out_file.write("\td_"+speciesId[i]+"=")
            if (species[i].isSetCompartment() == True):
                out_file.write("(")
            for k in range(0,numReactions):
                if(not stoichiometricMatrix[i][k]==0.0):
                    out_file.write("(")
                    out_file.write(repr(stoichiometricMatrix[i][k]))
                    out_file.write(")*(")
                    string=p.sub('',kineticLaw[k])
                    out_file.write(string)
                    out_file.write(")+")
            out_file.write("0")
            if (species[i].isSetCompartment() == True):
                out_file.write(")/")
                mySpeciesCompartment = species[i].getCompartment()
                for j in range(0, len(listOfParameter)):
                    if (listOfParameter[j].getId() == mySpeciesCompartment):
                        out_file.write(parameterId[j])
                        break
            out_file.write("\n")

    for i in range(0,len(listOfRules)):
        if listOfRules[i].isRate() == True:
            out_file.write("\td_")
            out_file.write(ruleVariable[i])
            out_file.write("=")
            out_file.write(ruleFormula[i])
            out_file.write("\n")
            

##################################################
#noise terms
##################################################

    out_file.write("\n")

    if(method==1):

        for i in range(0,numSpecies):
            if (species[i].getConstant() == False):
                out_file.write("\tnoise_"+speciesId[i]+"=")
                for k in range(0,numReactions):
                    if(not stoichiometricMatrix[i][k]==0.0):
                        out_file.write("(" + repr(stoichiometricMatrix[i][k]))
                        out_file.write(")*sqrt(")
                        string=p.sub('',kineticLaw[k])
                        out_file.write(string)
                        out_file.write(")*")
                        out_file.write("random.normal(0.0,sqrt(dt))")
                        out_file.write("+")
                out_file.write("0\n")

        for i in range(0,len(listOfRules)):
            if listOfRules[i].isRate() == True:
                out_file.write("\tnoise_")
                out_file.write(ruleVariable[i])
                out_file.write("= sqrt(")
                out_file.write(ruleFormula[i])
                out_file.write(" ) * random.normal(0.0,sqrt(dt))")
                out_file.write("\n")
                

    if(method==2):

        for i in range(0,numSpecies):
            if (species[i].getConstant() == False):
                out_file.write("\tnoise_"+speciesId[i]+"=")
                out_file.write("random.normal(0.0,sqrt(dt))\n")
                
        for i in range(0,len(listOfRules)):
            if listOfRules[i].isRate() == True:
                out_file.write("\tnoise_")
                out_file.write(ruleVariable[i])
                out_file.write("= ")
                out_file.write("random.normal(0.0,sqrt(dt))\n")


    if(method==3):
        
        for i in range(0,numSpecies):
            if (species[i].getConstant() == False):
                out_file.write("\tnoise_"+speciesId[i]+"=")
                for k in range(0,numReactions):
                    if(not stoichiometricMatrix[i][k]==0.0):
                        out_file.write("(")
                        out_file.write(repr(stoichiometricMatrix[i][k]))
                        out_file.write(")*(")
                        string=p.sub('',kineticLaw[k])
                        out_file.write(string)
                        out_file.write(")*")
                        out_file.write("random.normal(0.0,sqrt(dt))")
                        out_file.write("+")
                out_file.write("0\n")

        for i in range(0,len(listOfRules)):
            if listOfRules[i].isRate() == True:
                out_file.write("\tnoise_")
                out_file.write(ruleVariable[i])
                out_file.write("= (")
                out_file.write(ruleFormula[i])
                out_file.write(" ) * random.normal(0.0,sqrt(dt))")
                out_file.write("\n")
                 

    
    out_file.write("\n\treturn((")
                
    for i in range(0, len(species)):
        if (species[i].getConstant() == False):
            out_file.write("d_"+speciesId[i])
            out_file.write(",")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write("d_"+parameterId[i])
                    out_file.write(",")

    out_file.write("),(")
    for i in range(0,numSpecies):
        if (species[i].getConstant() == False):
            out_file.write("noise_"+speciesId[i])
            out_file.write(", ")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write("noise_"+parameterId[i])
                    out_file.write(",")

    out_file.write("))\n\n")

    
                
    #Write the assignment rules
    out_file.write("\ndef rules((")

    
    for i in range(0, len(species)):
        if (species[i].getConstant() == False):
            out_file.write(speciesId[i])
            out_file.write(",")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write(parameterId[i])
                    out_file.write(",")
    out_file.write("),(")
    for i in range(0,len(parameterId)):
        dontPrint = False
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]): 
                    dontPrint=True
        if dontPrint == False:
            out_file.write(parameterId[i])
            out_file.write(",")

    for i in range(0, numSpecies):
        if (species[i].getConstant() == True):
            out_file.write(speciesId[i])
            out_file.write(",")

    out_file.write("),time):\n\n")

	#Write the events

    for i in range(0, len(listOfEvents)):
        out_file.write("\tif ")
        out_file.write(mathMLConditionParser(EventCondition[i]))
        out_file.write(":\n")
        listOfAssignmentRules = listOfEvents[i].getListOfEventAssignments()
        for j in range(0, len(listOfAssignmentRules)):
            out_file.write("\t\t")
            out_file.write(EventVariable[i][j])
            out_file.write("=")
            out_file.write(EventFormula[i][j])
            out_file.write("\n")
        out_file.write("\n")

    out_file.write("\n")

    for i in range(0, len(listOfRules)):
        if listOfRules[i].isAssignment():
            out_file.write("\t")
            out_file.write(ruleVariable[i])
            out_file.write("=")
            out_file.write(ruleFormula[i])
            out_file.write("\n")

    out_file.write("\n\treturn((")
    for i in range(0, numSpecies):
        if (species[i].getConstant() == False):
            out_file.write(speciesId[i])
            out_file.write(",")
    for i in range(0, len(listOfParameter)):
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]):
                    out_file.write(parameterId[i])
                    out_file.write(",")
    out_file.write("),(")

    for i in range(0,len(parameterId)):
        dontPrint = False
        if (listOfParameter[i].getConstant() == False):
            for k in range(0, len(listOfRules)):
                if (listOfRules[k].isRate() and ruleVariable[k] == parameterId[i]): 
                    dontPrint=True
        if dontPrint == False:
            out_file.write(parameterId[i])
            out_file.write(",")

    for i in range(0, numSpecies):
        if (species[i].getConstant() == True):
            out_file.write(speciesId[i])
            out_file.write(",")
    out_file.write("))\n\n")               
    out_file.close()


################################################################################
# Function to get initial amount given a species and an algorithm type         #
# Need to pass to this a libsbml species object and a type an integration type #
################################################################################

def getSpeciesValue(species, intType):
    if species.isSetInitialAmount() and species.isSetInitialConcentration():
        if intType==ODE or intType==SDE:
            return species.getInitialConcentration()
        else: #implies intType = Gillespie
            return species.getInitialAmount()

    if species.isSetInitialAmount():
        return species.getInitialAmount()
    else:
        return species.getInitialConcentration()

##########################################
#Rename all parameters and species       #
##########################################
def rename(node,name,new_name):
    typ = node.getType()
    
    if (typ==AST_NAME or typ==AST_NAME_TIME):
        nme = node.getName()
        if nme == name:
            node.setName(new_name)

    for n in range(0,node.getNumChildren()):
        rename(node.getChild(n),name,new_name)
    return node


##############
# The PARSER #
##############

def importSBML(source,integrationType,ModelName=None,method=None):

    """
    ***** args *****
    source:
                  a list of strings.
                  Each tuple entry describes a SBML file to be parsed.

    integrationType:
                  a list of strings.
                  The length of this tuple is determined by the number of SBML 
                  files to be parsed. Each entry describes the simulation algorithm.
                  Possible algorithms are:
                  ODE         ---   for deterministic systems; solved with odeint (scipy)
                  SDE         ---   for stochastic systems; solved with sdeint (abc)
                  Gillespie   ---   for staochastic systems; solved with GillespieAlgorithm (abc)

    ***** kwargs *****
    ModelName:
                  a list of strings.
                  ModelName describes the names of the parsed model files.

    method:
                  an integer number.
                  Type of noise in a stochastic system.
                  (Only implemented for stochastic systems solved with sdeint.)
                  Possible options are:
                  1 --- default
                  2 --- Ornstein-Uhlenbeck
                  3 --- geometric Brownian motion

    """

    #regular expressions for detecting integration types
    g=re.compile('Gillespie')
    o=re.compile('ODE')
    s=re.compile('SDE')

    #check that you have appropriate lengths of integration types and sources
    #(need equal lengths)
    if(not(len(source)==len(integrationType))):
        print "\nError: Number of sources is not the same as number of integrationTypes!\n"
    #check that you have model names,
    #if not the models will be named model1, model2, etc
    else:
        if(ModelName==None):
            ModelName=[]
            for x in range(0,len(source)):
                ModelName.append("model"+repr(x+1))

        #if no method is specified and the integrationType is "SDE"
        #the method type defaults to 1
        for models in range(0,len(source)):
            intType = integrationType[models]
            if method==None:
                if s.match(integrationType[models]):
                    method=[]
                    for x in range(0, len(source)):
                        method.append(1)

        #All the below should happen for each model
        #Arguments to pass to the writing functions:
            #species IDs
            #species concentrations (initial values from model)
            #reactions in the form of kinetic law list
            #stoichiometric matrix
            #parameters
            #values of parameters
            #name of output file
            #list of functions if they need to be defined at the top of the written .py file

                        
            #I think that we can pass parameters directly to the writing functions, non?
            parameterId=[]
            parameterId2=[]
            parameter=[]
            listOfParameter=[]
            #Likewise species?
            speciesId=[]
            speciesId2=[]
            species=[]

##    r=re.compile('.mod')
##    if(r.search(source)):
##        old_source=source
##        source=r.sub(".xml",old_source)
##        call='python mod2sbml.py '+old_source+' > '+ source
##        os.system(call)

            #Get the model
            reader=SBMLReader()
            document=reader.readSBML(source[models])
            model=document.getModel()
            
            #get basic model properties
            numSpeciesTypes=model.getNumSpeciesTypes()
            numSpecies=model.getNumSpecies()
            numReactions=model.getNumReactions()
            numGlobalParameters=model.getNumParameters()
            numFunctions=model.getNumFunctionDefinitions()
            
            stoichiometricMatrix=empty([numSpecies, numReactions])


 
#################################################################################################
# get compartment volume/size - if it is set, pass as parameter with corresponding Id and value #
#################################################################################################

            listOfCompartments = model.getListOfCompartments()  
            comp=0
            for i in range(0, len(listOfCompartments)):
            #    listOfCompartments[i].setId('compartment'+repr(i+1))
                if listOfCompartments[i].isSetVolume():
                    comp=comp+1
                    parameterId.append(listOfCompartments[i].getId())
                    parameterId2.append('compartment'+repr(i+1))
                    parameter.append(listOfCompartments[i].getVolume())
                    listOfParameter.append(model.getCompartment(i))


#########################
# get global parameters #
#########################

            for i in range(0,numGlobalParameters):
                parameterId.append(model.getParameter(i).getId())
                parameterId2.append('parameter'+repr(i+1))
                parameter.append(model.getParameter(i).getValue())
                listOfParameter.append(model.getParameter(i))


###############
# get species #
###############

            #Empty matrix to hold reactants
            reactant=[]
            #Empty matrix to hold products
            product=[]
            #Empty matrix to hold Species Ids
            #Empty matrix to hold the InitValues used going forward
            InitValues=[]
            S1 = []
            S2 = []

            #Get a list of species
            listOfSpecies = model.getListOfSpecies()
            #Make the matrices long enough
            for k in range(0, len(listOfSpecies)):
                species.append(listOfSpecies[k])
                speciesId.append(listOfSpecies[k].getId())
                speciesId2.append('species'+repr(k+1))
                #get the initial value
                #Need to fix this part
                #So that it will take getInitialConcentration
                #or getInitialValue as appropriate
                InitValues.append(getSpeciesValue(listOfSpecies[k],intType))
                #I'm not really sure what this part is doing
                #Hopefully it will become more clear later
                S1.append(0.0)
                S2.append(0.0)
                #placeholder in reactant matrix for this species
                reactant.append(0)
                #placeholder in product matrix for this species
                product.append(0)

###############################
# analyse the model structure #
###############################

            reaction=[]
            numReactants=[]
            numProducts=[]
            kineticLaw=[]
            numLocalParameters=[]

            #Get the list of reactions
            listOfReactions = model.getListOfReactions()

            #For every reaction    
            for i in range(0, len(listOfReactions)):
                #What does this part do?
                for a in range(0, len(species)):
                    #what do S1 and S2 represent?
                    #S1 is something to do with stoichimetry of reactants
                    #At the moment S1 and S2 are as long as len(species)
                    S1[a]=0.0
                    #S2 is something to do with stoichiometry of products
                    S2[a]=0.0
        
                numReactants.append(listOfReactions[i].getNumReactants())
                numProducts.append(listOfReactions[i].getNumProducts())
                
                kineticLaw.append(listOfReactions[i].getKineticLaw().getFormula())
                numLocalParameters.append(listOfReactions[i].getKineticLaw().getNumParameters())

                for j in range(0, numReactants[i]):
                    reactant[j]=listOfReactions[i].getReactant(j)
        
                    for k in range(0,len(species)):
                        if (reactant[j].getSpecies()==species[k].getId()):
                            S1[k]=reactant[j].getStoichiometry()

                    
                for l in range(0,numProducts[i]):
                    product[l]=listOfReactions[i].getProduct(l)
        
                    for k in range(0,len(species)):
                        if (product[l].getSpecies()==species[k].getId()):
                            S2[k]=product[l].getStoichiometry()

                for m in range(0, len(species)):
                    stoichiometricMatrix[m][i]=-S1[m]+S2[m]

                for n in range(0,numLocalParameters[i]):
                    parameterId.append(listOfReactions[i].getKineticLaw().getParameter(n).getId())
                    parameterId2.append('parameter'+repr(len(parameterId)-comp))
                    parameter.append(listOfReactions[i].getKineticLaw().getParameter(n).getValue())
                    listOfParameter.append(listOfReactions[i].getKineticLaw().getParameter(n))

                    name=listOfReactions[i].getKineticLaw().getParameter(n).getId()
                    new_name='parameter'+repr(len(parameterId)-comp)
                    node=model.getReaction(i).getKineticLaw().getMath()
                    new_node=rename(node,name,new_name)
                    kineticLaw[i]=formulaToString(new_node)

                for n in range(0,comp):
                    
                    name=parameterId[n]
                    new_name='compartment'+repr(n+1)
                    node=model.getReaction(i).getKineticLaw().getMath()
                    new_node=rename(node,name,new_name)
                    kineticLaw[i]=formulaToString(new_node)

#####################
# analyse functions #
#####################

            #Get the list of functions

            listOfFunctions = model.getListOfFunctionDefinitions()


            FunctionArgument=[]
            FunctionBody=[]
                
            for fun in range(0,len(listOfFunctions)):
                FunctionArgument.append([])
                for funArg in range(0, listOfFunctions[fun].getNumArguments()):
                    FunctionArgument[fun].append(formulaToString(listOfFunctions[fun].getArgument(funArg)))

                FunctionBody.append(formulaToString(listOfFunctions[fun].getBody()))

            for fun in range(0, len(listOfFunctions)):
                for funArg in range(0,listOfFunctions[fun].getNumArguments()):
                    name=FunctionArgument[fun][funArg]
                    node=listOfFunctions[fun].getBody()
                    new_node=rename(node,name,"a"+repr(funArg+1))
                    FunctionBody[fun]=formulaToString(new_node)
                    FunctionArgument[fun][funArg]='a'+repr(funArg+1)
                   
        
#################
# analyse rules #
#################

            #Get the list of rules
            ruleFormula=[]
            ruleVariable=[]
            listOfRules = model.getListOfRules()
            for ru in range(0,len(listOfRules)):
                ruleFormula.append(listOfRules[ru].getFormula())
                ruleVariable.append(listOfRules[ru].getVariable())


##################
# analyse events #
##################
   
            listOfEvents = model.getListOfEvents()

            EventCondition=[]
            EventVariable=[]
            EventFormula=[]
           # listOfAssignmentRules=[]

            for eve in range(0,len(listOfEvents)):
                EventCondition.append(formulaToString(listOfEvents[eve].getTrigger().getMath()))
                listOfAssignmentRules=listOfEvents[eve].getListOfEventAssignments()
                EventVariable.append([])
                EventFormula.append([])
                for ru in range(0, len(listOfAssignmentRules)):
                   EventVariable[eve].append(listOfAssignmentRules[ru].getVariable())
                   EventFormula[eve].append(formulaToString(listOfAssignmentRules[ru].getMath()))


          
########################################################################
#rename parameters and species in reactions, events, rules             #
########################################################################

            NAMES=[[],[]]
            NAMES[0].append(parameterId)
            NAMES[0].append(parameterId2)
            NAMES[1].append(speciesId)
            NAMES[1].append(speciesId2)

            for nam in range(0,2):
                
                for i in range(0, len(NAMES[nam][0])):
                    name=NAMES[nam][0][i]
                    new_name=NAMES[nam][1][i]
             
                    for k in range(0, numReactions):
                        node=model.getReaction(k).getKineticLaw().getMath()
                        new_node=rename(node,name,new_name)
                        kineticLaw[k]=formulaToString(new_node)
                        
                    for k in range(0,len(listOfRules)):
                        node=listOfRules[k].getMath()
                        new_node=rename(node,name,new_name)
                        ruleFormula[k]=formulaToString(new_node)
                        if ruleVariable[k]==name: ruleVariable[k]=new_name


                    for k in range(0,len(listOfEvents)):
                        node=listOfEvents[k].getTrigger().getMath()
                        new_node=rename(node,name,new_name)
                        EventCondition[k]=formulaToString(new_node)
                        listOfAssignmentRules=listOfEvents[k].getListOfEventAssignments()
                        for cond in range(0, len(listOfAssignmentRules)):
                            node=listOfAssignmentRules[cond].getMath()
                            new_node=rename(node,name,new_name)
                            EventFormula[k][cond]=formulaToString(new_node)
                            if EventVariable[k][cond]==name: EventVariable[k][cond]=new_name
                        

          
##########################
# call writing functions #
##########################

            s=re.compile('SDE')
            if o.match(integrationType[models]):
                write_ODEFunctions(stoichiometricMatrix, kineticLaw, species, numReactions,speciesId2,listOfParameter, parameterId2,parameter,InitValues,ModelName[models], listOfFunctions,FunctionArgument,FunctionBody,listOfRules, ruleFormula,ruleVariable, listOfEvents,EventCondition,EventVariable,EventFormula)
            if s.match(integrationType[models]):
                write_SDEFunctions(stoichiometricMatrix, kineticLaw, species, numReactions,speciesId2,listOfParameter, parameterId2,parameter,InitValues,method[models],ModelName[models], listOfFunctions,FunctionArgument,FunctionBody, listOfRules, ruleFormula, ruleVariable, listOfEvents, EventCondition, EventVariable, EventFormula)
            if g.match(integrationType[models]):
                  write_GillespieFunctions(stoichiometricMatrix, kineticLaw, numSpecies, numReactions, species, parameterId2, InitValues, speciesId2,ModelName[models], listOfFunctions,FunctionArgument,FunctionBody, listOfRules, ruleFormula, ruleVariable, listOfEvents, EventCondition, EventVariable, EventFormula)
