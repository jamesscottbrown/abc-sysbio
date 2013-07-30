from libsbml import *
from abcsysbio.relations import *
import os
import re
from Writer import Writer

class SDEPythonWriter(Writer):
    def __init__(self, sbmlFileName, modelName="", inputPath="", outputPath=""):
        Writer.__init__(self, sbmlFileName, modelName, inputPath, outputPath)
        self.out_file = open(os.path.join(outputPath,self.parsedModel.name+".py"),"w")
        
    def write(self, method):
        p=re.compile('\s')
        self.out_file.write("from math import sqrt\nfrom numpy import random\nfrom abcsysbio.relations import *\n\n")

        self.out_file.write("def trunc_sqrt(x):\n\tif x < 0: return 0\n\telse: return sqrt(x)\n\n")

        for i in range(0, len(self.parsedModel.listOfFunctions)):
            self.out_file.write("def ")
            self.out_file.write(self.parsedModel.listOfFunctions[i].getId())
            self.out_file.write("(")
            for j in range(0, self.parsedModel.listOfFunctions[i].getNumArguments()):
                self.out_file.write(self.parsedModel.functionArgument[i][j])
                self.out_file.write(",")
            self.out_file.write("):\n\n\toutput=")
            self.out_file.write(self.parsedModel.functionBody[i])
            self.out_file.write("\n\n\treturn output\n\n")
            
        self.out_file.write("def modelfunction((")
    
        for i in range(0, len(self.parsedModel.species)):
            ##if (self.parsedModel.species[i].getConstant() == False):
            self.out_file.write(self.parsedModel.speciesId[i])
            self.out_file.write(",")
        for i in range(0, len(self.parsedModel.listOfParameter)):
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]):
                        self.out_file.write(self.parsedModel.parameterId[i])
                        self.out_file.write(",")
    
        self.out_file.write(")=(")
    
        for i in range(0, len(self.parsedModel.species)):
            ##if (self.parsedModel.species[i].getConstant() == False):
            self.out_file.write(repr(self.parsedModel.initValues[i]))
            self.out_file.write(",")
        for i in range(0, len(self.parsedModel.listOfParameter)):
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]):
                        self.out_file.write(repr(self.parsedModel.parameter[i]))
                        self.out_file.write(",")
    
        self.out_file.write("),dt=0,parameter=(")
    
        for i in range(0,len(self.parsedModel.parameterId)):
            dontPrint = False
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]): 
                        dontPrint=True
            if dontPrint == False:
                self.out_file.write(repr(self.parsedModel.parameter[i]))
                self.out_file.write(",")
    
        ##for i in range(0, self.parsedModel.numSpecies):
            ##if (self.parsedModel.species[i].getConstant() == True):
            ##    self.out_file.write(repr(self.parsedModel.initValues[i]))
            ##    self.out_file.write(",")
     
        self.out_file.write("),time=0):\n\n")
    
        counter=0
        for i in range(0,len(self.parsedModel.parameterId)):
            dontPrint = False
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]): dontPrint=True
            if dontPrint == False:
                self.out_file.write("\t"+self.parsedModel.parameterId[i]+"=parameter["+repr(counter)+"]\n")
                counter=counter+1
    
        ##for i in range(0, self.parsedModel.numSpecies):
            ##if (self.parsedModel.species[i].getConstant() == True):
            ##    self.out_file.write("\t"+self.parsedModel.speciesId[i]+"=parameter["+repr(counter)+"]\n")
            ##    counter = counter+1
    
        self.out_file.write("\n")
        
    
        
        self.out_file.write("\n")
        
        for i in range(0,self.parsedModel.numSpecies):
            ##if (self.parsedModel.species[i].getConstant() == False):
            self.out_file.write("\td_"+self.parsedModel.speciesId[i]+"=")
            if (self.parsedModel.species[i].isSetCompartment() == True):
                self.out_file.write("(")
            for k in range(0,self.parsedModel.numReactions):
                if(not self.parsedModel.stoichiometricMatrix[i][k]==0.0):
                    self.out_file.write("(")
                    self.out_file.write(repr(self.parsedModel.stoichiometricMatrix[i][k]))
                    self.out_file.write(")*(")
                    string=p.sub('',self.parsedModel.kineticLaw[k])
                    self.out_file.write(string)
                    self.out_file.write(")+")
            self.out_file.write("0")
            if (self.parsedModel.species[i].isSetCompartment() == True):
                self.out_file.write(")/")
                mySpeciesCompartment = self.parsedModel.species[i].getCompartment()
                for j in range(0, len(self.parsedModel.listOfParameter)):
                    if (self.parsedModel.listOfParameter[j].getId() == mySpeciesCompartment):
                        self.out_file.write(self.parsedModel.parameterId[j])
                        break
            self.out_file.write("\n")
    
        for i in range(0,len(self.parsedModel.listOfRules)):
            if self.parsedModel.listOfRules[i].isRate() == True:
                self.out_file.write("\td_")
                self.out_file.write(self.parsedModel.ruleVariable[i])
                self.out_file.write("=")
                self.out_file.write(self.parsedModel.ruleFormula[i])
                self.out_file.write("\n")
                
    
    ##################################################
    #noise terms
    ##################################################
    
        self.out_file.write("\n")

        # check for columns of the stochiometry matrix with more than one entry
        randomVariables = ["*random.normal(0.0,sqrt(dt))"] * self.parsedModel.numReactions
        for k in range(0,self.parsedModel.numReactions):
            countEntries = 0
            for i in range(0,self.parsedModel.numSpecies):
                if(self.parsedModel.stoichiometricMatrix[i][k] != 0.0): countEntries += 1
                
            # define specific randomVariable
            if countEntries > 1:
                self.out_file.write("\trand"+repr(k)+" = random.normal(0.0,sqrt(dt))\n")
                randomVariables[k] = "*rand" + repr(k)
                    
    
        if(method==1):
    
            for i in range(0,self.parsedModel.numSpecies):
                ##if (self.parsedModel.species[i].getConstant() == False):
                self.out_file.write("\tnoise_"+self.parsedModel.speciesId[i]+"=")
                for k in range(0,self.parsedModel.numReactions):
                    if(not self.parsedModel.stoichiometricMatrix[i][k]==0.0):
                        self.out_file.write("(" + repr(self.parsedModel.stoichiometricMatrix[i][k]))
                        self.out_file.write(")*trunc_sqrt(")
                        string=p.sub('',self.parsedModel.kineticLaw[k])
                        self.out_file.write(string)
                        self.out_file.write(")")
                        self.out_file.write(randomVariables[k])
                        self.out_file.write("+")
                self.out_file.write("0\n")
    
            for i in range(0,len(self.parsedModel.listOfRules)):
                if self.parsedModel.listOfRules[i].isRate() == True:
                    self.out_file.write("\tnoise_")
                    self.out_file.write(self.parsedModel.ruleVariable[i])
                    self.out_file.write("= trunc_sqrt(")
                    self.out_file.write(self.parsedModel.ruleFormula[i])
                    self.out_file.write(")")
                    self.out_file.write(randomVariables[k])
                    self.out_file.write("\n")
                    
    
        if(method==2):
    
            for i in range(0,self.parsedModel.numSpecies):
                ##if (self.parsedModel.species[i].getConstant() == False):
                self.out_file.write("\tnoise_"+self.parsedModel.speciesId[i]+"=")
                self.out_file.write("random.normal(0.0,sqrt(dt))\n")
                    
            for i in range(0,len(self.parsedModel.listOfRules)):
                if self.parsedModel.listOfRules[i].isRate() == True:
                    self.out_file.write("\tnoise_")
                    self.out_file.write(self.parsedModel.ruleVariable[i])
                    self.out_file.write("= ")
                    self.out_file.write("random.normal(0.0,sqrt(dt))\n")
    
    
        if(method==3):
            
            for i in range(0,self.parsedModel.numSpecies):
                ##if (self.parsedModel.species[i].getConstant() == False):
                self.out_file.write("\tnoise_"+self.parsedModel.speciesId[i]+"=")
                for k in range(0,self.parsedModel.numReactions):
                    if(not self.parsedModel.stoichiometricMatrix[i][k]==0.0):
                        self.out_file.write("(")
                        self.out_file.write(repr(self.parsedModel.stoichiometricMatrix[i][k]))
                        self.out_file.write(")*(")
                        string=p.sub('',self.parsedModel.kineticLaw[k])
                        self.out_file.write(string)
                        self.out_file.write(")*")
                        self.out_file.write("random.normal(0.0,sqrt(dt))")
                        self.out_file.write("+")
                self.out_file.write("0\n")
    
            for i in range(0,len(self.parsedModel.listOfRules)):
                if self.parsedModel.listOfRules[i].isRate() == True:
                    self.out_file.write("\tnoise_")
                    self.out_file.write(self.parsedModel.ruleVariable[i])
                    self.out_file.write("= (")
                    self.out_file.write(self.parsedModel.ruleFormula[i])
                    self.out_file.write(" ) * random.normal(0.0,sqrt(dt))")
                    self.out_file.write("\n")
                     
    
        
        self.out_file.write("\n\treturn((")
                    
        for i in range(0, len(self.parsedModel.species)):
            ##if (self.parsedModel.species[i].getConstant() == False):
            self.out_file.write("d_"+self.parsedModel.speciesId[i])
            self.out_file.write(",")
        for i in range(0, len(self.parsedModel.listOfParameter)):
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]):
                        self.out_file.write("d_"+self.parsedModel.parameterId[i])
                        self.out_file.write(",")
    
        self.out_file.write("),(")
        for i in range(0,self.parsedModel.numSpecies):
            ##if (self.parsedModel.species[i].getConstant() == False):
            self.out_file.write("noise_"+self.parsedModel.speciesId[i])
            self.out_file.write(", ")
        for i in range(0, len(self.parsedModel.listOfParameter)):
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]):
                        self.out_file.write("noise_"+self.parsedModel.parameterId[i])
                        self.out_file.write(",")
    
        self.out_file.write("))\n\n")
    
        
                    
        #Write the assignment rules
        self.out_file.write("\ndef rules((")
    
        
        for i in range(0, len(self.parsedModel.species)):
            ##if (self.parsedModel.species[i].getConstant() == False):
            self.out_file.write(self.parsedModel.speciesId[i])
            self.out_file.write(",")
        for i in range(0, len(self.parsedModel.listOfParameter)):
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]):
                        self.out_file.write(self.parsedModel.parameterId[i])
                        self.out_file.write(",")
        self.out_file.write("),(")
        for i in range(0,len(self.parsedModel.parameterId)):
            dontPrint = False
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]): 
                        dontPrint=True
            if dontPrint == False:
                self.out_file.write(self.parsedModel.parameterId[i])
                self.out_file.write(",")
    
        ##for i in range(0, self.parsedModel.numSpecies):
            ##if (self.parsedModel.species[i].getConstant() == True):
            ##    self.out_file.write(self.parsedModel.speciesId[i])
            ##    self.out_file.write(",")
    
        self.out_file.write("),time):\n\n")
    
        #Write the events
    
        for i in range(0, len(self.parsedModel.listOfEvents)):
            self.out_file.write("\tif ")
            self.out_file.write(mathMLConditionParser(self.parsedModel.eventCondition[i]))
            self.out_file.write(":\n")
            listOfAssignmentRules = self.parsedModel.listOfEvents[i].getListOfEventAssignments()
            for j in range(0, len(listOfAssignmentRules)):
                self.out_file.write("\t\t")
                self.out_file.write(self.parsedModel.eventVariable[i][j])
                self.out_file.write("=")
                self.out_file.write(self.parsedModel.eventFormula[i][j])
                self.out_file.write("\n")
            self.out_file.write("\n")
    
        self.out_file.write("\n")
    
        for i in range(0, len(self.parsedModel.listOfRules)):
            if self.parsedModel.listOfRules[i].isAssignment():
                self.out_file.write("\t")
                self.out_file.write(self.parsedModel.ruleVariable[i])
                self.out_file.write("=")
                self.out_file.write(self.parsedModel.ruleFormula[i])
                self.out_file.write("\n")
    
        self.out_file.write("\n\treturn((")
        for i in range(0, self.parsedModel.numSpecies):
            ##if (self.parsedModel.species[i].getConstant() == False):
            self.out_file.write(self.parsedModel.speciesId[i])
            self.out_file.write(",")
        for i in range(0, len(self.parsedModel.listOfParameter)):
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]):
                        self.out_file.write(self.parsedModel.parameterId[i])
                        self.out_file.write(",")
        self.out_file.write("),(")
    
        for i in range(0,len(self.parsedModel.parameterId)):
            dontPrint = False
            if (self.parsedModel.listOfParameter[i].getConstant() == False):
                for k in range(0, len(self.parsedModel.listOfRules)):
                    if (self.parsedModel.listOfRules[k].isRate() and self.parsedModel.ruleVariable[k] == self.parsedModel.parameterId[i]): 
                        dontPrint=True
            if dontPrint == False:
                self.out_file.write(self.parsedModel.parameterId[i])
                self.out_file.write(",")
    
        ##for i in range(0, self.parsedModel.numSpecies):
            ##if (self.parsedModel.species[i].getConstant() == True):
                ##self.out_file.write(self.parsedModel.speciesId[i])
                ##self.out_file.write(",")
        self.out_file.write("))\n\n")               
        self.out_file.close()
