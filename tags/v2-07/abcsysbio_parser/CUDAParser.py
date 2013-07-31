from Parser import Parser
from libsbml import formulaToString
import re

class CUDAParser(Parser):
    def __init__(self, sbmlFileName, modelName, integrationType, method, inputPath="", outputPath=""):
        self.mathPython = []
        self.mathCuda = []
        Parser.__init__(self, sbmlFileName, modelName, integrationType, method, inputPath="", outputPath="")
        
    def parse(self):
        Parser.parse(self)
        self.getGlobalParameters()
        self.getSpecies()
        self.analyseModelStructure()
        self.analyseFunctions()
        self.analyseRules()
        self.analyseEvents()
        self.renameMathFunctions()
        self.renameEverything()
        self.writer.parsedModel.numGlobalParameters = self.writer.parsedModel.numGlobalParameters + 1
        
    def getGlobalParameters(self):
        Parser.getGlobalParameters(self)
        for i in range(0, self.writer.parsedModel.numGlobalParameters):
            if ((len(self.writer.parsedModel.parameterId) - self.comp) < 9):
                self.writer.parsedModel.parameterId.append("parameter0" + repr(i + 1))
            else:
                self.writer.parsedModel.parameterId.append("parameter" + repr(i + 1))
    def getSpecies(self):
        Parser.getSpecies(self)
        for k in range(0, len(self.listOfSpecies)):
            if ((len(self.writer.parsedModel.speciesId) - self.comp) < 9):
                self.writer.parsedModel.speciesId.append("species0" + repr(k + 1))
            else:
                self.writer.parsedModel.speciesId.append("species" + repr(k + 1))
            
    def analyseModelStructure(self):
        Parser.analyseModelStructure(self)
        for i in range(0, len(self.listOfReactions)):
            for n in range(0, self.numLocalParameters[i]):
                self.parameterId.append(self.listOfReactions[i].getKineticLaw().getParameter(n).getId())
		if ((len(self.writer.parsedModel.parameterId) - self.comp) < 10):
                    self.writer.parsedModel.parameterId.append("parameter0" + repr(len(self.parameterId) - self.comp))
                else:
                    self.writer.parsedModel.parameterId.append("parameter" + repr(len(self.parameterId) - self.comp))
        	name = self.listOfReactions[i].getKineticLaw().getParameter(n).getId()
                new_name = 'parameter' + repr(len(self.parameterId) - self.comp)
                node = self.sbmlModel.getReaction(i).getKineticLaw().getMath()
                new_node = self.rename(node, name, new_name)
                self.writer.parsedModel.kineticLaw[i] = formulaToString(new_node)

    def renameMathFunctions(self):
        self.mathPython.append('log10')
        self.mathPython.append('acos')
        self.mathPython.append('asin')
        self.mathPython.append('atan')
        self.mathPython.append('exp')
        self.mathPython.append('sqrt')
        self.mathPython.append('pow')
        self.mathPython.append('log')
        self.mathPython.append('sin')
        self.mathPython.append('cos')
        self.mathPython.append('ceil')
        self.mathPython.append('floor')
        self.mathPython.append('tan')
        self.mathPython.append('time')
        
        self.mathCuda.append('log10')
        self.mathCuda.append('acos')
        self.mathCuda.append('asin')
        self.mathCuda.append('atan')
        self.mathCuda.append('exp')
        self.mathCuda.append('sqrt')
        self.mathCuda.append('mpow')
        self.mathCuda.append('log')
        self.mathCuda.append('sin')
        self.mathCuda.append('cos')
        self.mathCuda.append('ceil')
        self.mathCuda.append('floor')
        self.mathCuda.append('tan')
    
    def renameEverything(self):
        Parser.renameEverything(self)
        for nam in range(0, len(self.mathPython)): 
            for k in range(0, len(self.writer.parsedModel.kineticLaw)):
                if re.search(self.mathPython[nam], self.writer.parsedModel.kineticLaw[k]):
                    s = self.writer.parsedModel.kineticLaw[k]
                    s = re.sub(self.mathPython[nam], self.mathCuda[nam], s)
                    self.writer.parsedModel.kineticLaw[k] = s
    
            for k in range(0, len(self.writer.parsedModel.ruleFormula)):
                if re.search(self.mathPython[nam], self.writer.parsedModel.ruleFormula[k]):
                    s = self.writer.parsedModel.ruleFormula[k]
                    s = re.sub(self.mathPython[nam], self.mathCuda[nam], s)
                    self.writer.parsedModel.ruleFormula[k] = s
    
            for k in range(0, len(self.writer.parsedModel.eventFormula)):
                for cond in range(0, len(self.listOfAssignmentRules)):
                    if re.search(self.mathPython[nam], self.writer.parsedModel.eventFormula[k][cond]):
                        s = self.writer.parsedModel.eventFormula[k][cond]
                        s = re.sub(self.mathPython[nam], self.mathCuda[nam], s)
                        self.writer.parsedModel.eventFormula[k][cond] = s
    
            for k in range(0, len(self.writer.parsedModel.eventCondition)):
                if re.search(self.mathPython[nam], self.writer.parsedModel.eventCondition[k]):
                    s = self.writer.parsedModel.eventCondition[k]
                    s = re.sub(self.mathPython[nam], self.mathCuda[nam], s)
                    self.writer.parsedModel.eventCondition[k] = s
    
            for k in range(0, len(self.writer.parsedModel.functionBody)):
                if re.search(self.mathPython[nam], self.writer.parsedModel.functionBody[k]):
                    s = self.writer.parsedModel.functionBody[k]
                    s = re.sub(self.mathPython[nam], self.mathCuda[nam], s)
                    self.writer.parsedModel.functionBody[k] = s
    
            for fun in range(0, len(self.writer.parsedModel.listOfFunctions)):
                for k in range(0, len(self.writer.parsedModel.functionArgument[fun])):
                    if re.search(self.mathPython[nam], self.writer.parsedModel.functionArgument[fun][k]):
                        s = self.writer.parsedModel.functionArgument[fun][k]
                        s = re.sub(self.mathPython[nam], self.mathCuda[nam], s)
                        self.writer.parsedModel.functionArgument[fun][k] = s
        
