from libsbml import *
from abcsysbio.relations import *
from Writer import Writer
import os
import re

class CWriter(Writer):
    def __init__(self, sbmlFileName, modelName="", inputPath="", outputPath=""):
        Writer.__init__(self, sbmlFileName, modelName, inputPath, outputPath)
        
        self.hppOutputFile = open(os.path.join(outputPath,self.parsedModel.name+".hpp"), "w")
        self.cppOutputFile = open(os.path.join(outputPath,self.parsedModel.name+".cpp"), "w")
        
    def write(self):
        self.writeCheader()
        self.writeCsourceCode()
            
    def writeCheader(self):
            self.hppOutputFile.write("#ifndef ")
            self.hppOutputFile.write(self.parsedModel.name.upper())
            self.hppOutputFile.write("_HPP_\n")
            
            self.hppOutputFile.write("#define ")
            self.hppOutputFile.write(self.parsedModel.name.upper())
            self.hppOutputFile.write("_HPP_\n")
            
            self.hppOutputFile.write("""
		
		#include <vector>
		#include <iostream>
		#include "newmat.h"
		#include "newmatio.h"
		#include "newmatap.h"
    		class ChildModel {
              public: 
            
              /**
               * Number of reactions of the model
               */
              int NREACTIONS;
              
              /**
               * Number of species of the model
               */
              int NSPECIES;
              
              /**
               * Stoichiometric Matrix of the system (the rows represent the species and the columns the reactions)
               */
              Matrix* pstoichiometricMatrix;
            
              ChildModel(int i);
              void init();
            
             /**
               * Virtual method (ie method defined in the child class) setting the values of the stoichiometric matrix
               *
               * @param void
               * @return void
               */
              void getStoichiometricMatrix();
              
              /**
               * Virtual method computing the hazards of the different reactions for a given concentration of species (yi) and some parameter values
               *
               * @param double concentrations[] Array of size NSPECIES containing the concentrations of the species for which we want to compute the hazards
               * @param double parameters[] Array containing the parameter's values for which we want to compute the hazards (the number of parameters depend on the model and doesn't have to be the number of reactions)
               */
              ColumnVector getHazards(const double concentrations[],
            				  const double parameters[]);
              
              /**
               * Virtual method modifying the concentrations and parameters depending on some criteria defined by the SBML
               *
               * @param double concentrations[] Array of size NSPECIES containing the concentrations of the species
               * @param double parameters[] Array containing the parameter's values
               */
              void applyRulesAndEvents(double concentrations[],
            				   double parameters[], double time);
            """)
            
            for i in range(0, len(self.parsedModel.listOfFunctions)):
               
                self.hppOutputFile.write("double ")
                string = self.parsedModel.listOfFunctions[i].getId()
                string = re.sub('_', '', string)
                self.hppOutputFile.write(string)
                self.hppOutputFile.write("(")
                self.hppOutputFile.write("\tdouble\t" + self.parsedModel.functionArgument[i][0])
                
                for j in range(1, self.parsedModel.listOfFunctions[i].getNumArguments()):
                   self.hppOutputFile.write(", ")
                   self.hppOutputFile.write("double " + self.parsedModel.functionArgument[i][j])
                self.hppOutputFile.write(");\n")
            
            self.hppOutputFile.write('\n};\n')
            self.hppOutputFile.write('#endif /*')
            self.hppOutputFile.write(self.parsedModel.name.upper())
            self.hppOutputFile.write('_HPP_ */\n')
            
    def writeCsourceCode(self):
            p1 = re.compile('species(\d+)')
            p2 = re.compile('parameter(\d+)')
        
            #self.cppOutputFile.write('#include "' + self.parsedModel.name + '.hpp"\n')
            self.cppOutputFile.write('#include "ChildModel.hpp"\n')
            self.cppOutputFile.write('#include <cmath>\n')
            self.writeModelConstructor()
            self.writeUserDefinedFunctions()
            self.writeStoichiometricMatrix()
            self.writeGetHazardFunction(p1,p2)
            self.writeRulesAndEvents(p1,p2)
        
    def  writeModelConstructor(self):
    
            self.cppOutputFile.write("\nChildModel::ChildModel(int i){")
            self.cppOutputFile.write("\n\tNSPECIES = " + str(self.parsedModel.numSpecies) + ";")                                    
            self.cppOutputFile.write("\n\tNREACTIONS = " + str(self.parsedModel.numReactions) + ";")                                  
            self.cppOutputFile.write("\n\tpstoichiometricMatrix = new Matrix(NSPECIES,NREACTIONS);")
            self.cppOutputFile.write("\n\t(*pstoichiometricMatrix) = 0.0;");        
            self.cppOutputFile.write("\n\tgetStoichiometricMatrix();");
            self.cppOutputFile.write("\n}");
    
    def  writeUserDefinedFunctions(self):
    
            #The user-defined functions used in the model must be written in the file
       
            for i in range(0, len(self.parsedModel.listOfFunctions)):
                self.cppOutputFile.write("double ChildModel::")
                string = self.parsedModel.listOfFunctions[i].getId()
                string = re.sub('_', '', string)
                self.cppOutputFile.write(string)
                self.cppOutputFile.write("(")
    
                self.cppOutputFile.write("double  " + self.parsedModel.functionArgument[i][0])
                for j in range(1, self.parsedModel.listOfFunctions[i].getNumArguments()):
                    self.cppOutputFile.write(",")
                    self.cppOutputFile.write(" double  " + self.parsedModel.functionArgument[i][j])
                self.cppOutputFile.write("){\n\n\t\tdouble output=")
                self.cppOutputFile.write(self.parsedModel.functionBody[i] + ";")
                self.cppOutputFile.write("\n\n\t\treturn output;\n\t}\n")
    
    def  writeStoichiometricMatrix(self):
    
            self.cppOutputFile.write("\n\n\tvoid ChildModel::getStoichiometricMatrix() {")
        
            for i in range(0, self.parsedModel.numReactions):
                for k in range(0, self.parsedModel.numSpecies):
                    ##if (self.parsedModel.species[k].getConstant() == False):
                    self.cppOutputFile.write("\n\t\t (*pstoichiometricMatrix)(" + repr(k) + "+1," + repr(i) + "+1)= " + str(self.parsedModel.stoichiometricMatrix[k][i]) + ";")    
            self.cppOutputFile.write("\n\t}")
        
    def writeGetHazardFunction(self,p1,p2):
            self.cppOutputFile.write("\n\n\tColumnVector ChildModel::getHazards(const double concentrations[],const double parameters[]) {")
            self.cppOutputFile.write("\n\t\tColumnVector hazards(NREACTIONS);\n")
            for i in range(0, self.parsedModel.numReactions): 
                string = self.parsedModel.kineticLaw[i];
                string = re.sub('_', '', string)
                string = p1.sub(r"concentrations[\g<1>-1]", string);
                string = p2.sub(r"parameters[\g<1>]", string);
                string = re.sub("compartment1", "parameters[0]", string);
                self.cppOutputFile.write("\n\t\thazards(" + repr(i) + "+1) = " + string)
                self.cppOutputFile.write(";\n")
            self.cppOutputFile.write("\t\treturn hazards;\n")
            self.cppOutputFile.write("\t}\n")
    
    def writeRulesAndEvents(self,p1,p2):
        
              #Write the rules and events
            self.cppOutputFile.write("\n\tvoid ChildModel::applyRulesAndEvents(double concentrations[], double parameters[], double time) {\n")
            
            self.writeEvents(p1,p2)
            self.writeRules(p1,p2)
            
            self.cppOutputFile.write("\n\t}\n")
            
    def writeEvents(self,p1,p2):
            #Write the events
            
            for i in range(0, len(self.parsedModel.listOfEvents)):
                self.cppOutputFile.write("\t\tif ")
                string = mathMLConditionParser(self.parsedModel.eventCondition[i])
                string = re.sub(',', '>=', string)
                string = re.sub("geq", " ", string)
                self.cppOutputFile.write(string)
                self.cppOutputFile.write("{\n")
                listOfAssignmentRules = self.parsedModel.listOfEvents[i].getListOfEventAssignments()
                
                for j in range(0, len(listOfAssignmentRules)):
                    self.cppOutputFile.write("\t\t\t")
        
                    string = self.parsedModel.eventVariable[i][j]
                    string = re.sub('_', '', string)
                    string = p1.sub(r"concentrations[\g<1>-1]", string);
                    string = p2.sub(r"parameters[\g<1>]", string);
                    string = re.sub("compartment1", "parameters[0]", string);
                    self.cppOutputFile.write(string)
                    
                    self.cppOutputFile.write("=")
                    
                    string = self.parsedModel.eventFormula[i][j]
                    string = re.sub('_', '', string)
                    string = p1.sub(r"concentrations[\g<1>-1]", string); 
                    string = p2.sub(r"parameters[\g<1>]", string);
                    string = re.sub("compartment1", "parameters[0]", string);
                    self.cppOutputFile.write(string)
        
                    self.cppOutputFile.write(";\n\t\t}\n")
                self.cppOutputFile.write("\n")
        
            self.cppOutputFile.write("\n")
            
    def writeRules(self,p1,p2):
            #write the rules
        
            for i in range(0, len(self.parsedModel.listOfRules)):
                if self.parsedModel.listOfRules[i].isAssignment():
                    self.cppOutputFile.write("\t\t")
                    string = self.parsedModel.ruleVariable[i];
                    string = re.sub('_', '', string)
                    string = p1.sub(r"concentrations[\g<1>-1]", string); 
                    string = p2.sub(r"parameters[\g<1>]", string);
                    string = re.sub("compartment1", "parameters[0]", string);
                    self.cppOutputFile.write(string)
        
                    self.cppOutputFile.write("=")
        
                    string = mathMLConditionParser(self.parsedModel.ruleFormula[i])
                    string = re.sub('_', '', string)
                    string = p1.sub(r"concentrations[\g<1>-1]", string); 
                    string = p2.sub(r"parameters[\g<1>]", string);
                    string = re.sub("compartment1", "parameters[0]", string);
                    self.cppOutputFile.write(string)
        
                    self.cppOutputFile.write(";\n")
