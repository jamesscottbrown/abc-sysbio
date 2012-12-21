from Parser import Parser
from libsbml import *
from numpy import *
import re
import os
from abcsysbio.relations import *
from ParsedModel import ParsedModel

class CandPythonParser(Parser):
    def __init__(self, sbmlFileName, modelName, integrationType, method, inputPath="", outputPath=""):
        Parser.__init__(self, sbmlFileName, modelName, integrationType, method, inputPath, outputPath)
        
    def parse(self):
        Parser.parse(self)
        self.getGlobalParameters()
        self.getSpecies()
        self.analyseModelStructure()
        self.analyseFunctions()
        self.analyseRules()
        self.analyseEvents()
        self.renameEverything()
            
    def getGlobalParameters(self):
        Parser.getGlobalParameters(self)
        for i in range(0,self.writer.parsedModel.numGlobalParameters):
            self.writer.parsedModel.parameterId.append("parameter"+repr(i+1))
            
    def getSpecies(self):
        Parser.getSpecies(self)
        for k in range(0, len(self.listOfSpecies)):
            self.writer.parsedModel.speciesId.append("species"+repr(k+1))
                
    def analyseModelStructure(self):
        Parser.analyseModelStructure(self)
        for i in range(0, len(self.listOfReactions)):
            for n in range(0, self.numLocalParameters[i]):
                self.parameterId.append(self.listOfReactions[i].getKineticLaw().getParameter(n).getId())
		self.writer.parsedModel.parameterId.append("parameter"+repr(len(self.parameterId)-self.comp))
        	name = self.listOfReactions[i].getKineticLaw().getParameter(n).getId()
                new_name = 'parameter' + repr(len(self.parameterId) - self.comp)
                node = self.sbmlModel.getReaction(i).getKineticLaw().getMath()
                new_node = self.rename(node, name, new_name)
                self.writer.parsedModel.kineticLaw[i] = formulaToString(new_node)
    
