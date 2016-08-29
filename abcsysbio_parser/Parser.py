from numpy import *
from libsbml import *
import re
import os
from abcsysbio.relations import *
from CWriter import CWriter
from SDEPythonWriter import SDEPythonWriter
from ODEPythonWriter import ODEPythonWriter
from GillespiePythonWriter import GillespiePythonWriter
from SDECUDAWriter import SdeCUDAWriter
from ODECUDAWriter import OdeCUDAWriter
from GillespieCUDAWriter import GillespieCUDAWriter


class Parser:
    def __init__(self, sbmlFileName, modelName, integrationType, method, inputPath="", outputPath=""):
        c = re.compile('C', re.IGNORECASE)
        py = re.compile('Python', re.I)
        cuda = re.compile('CUDA', re.I)

        gil = re.compile('Gillespie', re.I)
        ode = re.compile('ODE', re.I)
        sde = re.compile('SDE', re.I)
        euler = re.compile('Euler', re.I)
        heun = re.compile('Heun', re.I)
        milstein = re.compile('Milstein', re.I)

        if cuda.search(integrationType):
            if gil.search(integrationType):
                self.writer = GillespieCUDAWriter(sbmlFileName, modelName, inputPath, outputPath)
            elif ode.search(integrationType):
                self.writer = OdeCUDAWriter(sbmlFileName, modelName, inputPath, outputPath)
            elif sde.search(integrationType):
                self.writer = SdeCUDAWriter(sbmlFileName, modelName, inputPath, outputPath)

        elif c.search(integrationType):
            self.writer = CWriter(sbmlFileName, modelName, inputPath, outputPath)

        elif py.search(integrationType):
            if gil.search(integrationType):
                self.writer = GillespiePythonWriter(sbmlFileName, modelName, inputPath, outputPath)
            elif ode.search(integrationType):
                self.writer = ODEPythonWriter(sbmlFileName, modelName, inputPath, outputPath)
            elif sde.search(integrationType):
                self.writer = SDEPythonWriter(sbmlFileName, modelName, inputPath, outputPath)

        reader = SBMLReader()
        document = reader.readSBML(inputPath + sbmlFileName)
        self.sbmlModel = document.getModel()

        self.parameterId = []

        self.listOfSpecies = []  # Used by the child
        self.speciesId = []

        self.product = []
        self.reactant = []
        self.S1 = []
        self.S2 = []

        self.listOfReactions = []  # Used by the child
        self.listOfAssignmentRules = []
        self.numLocalParameters = []  # Used by the child

        self.comp = 0
        self.parse()
        if (py.search(integrationType) or cuda.search(integrationType)) and sde.search(integrationType):
            self.writer.write(method)
        else:
            self.writer.write()

    def parse(self):
        self.getBasicModelProperties()
        self.writer.parsedModel.stoichiometricMatrix = empty(
            [self.writer.parsedModel.numSpecies, self.writer.parsedModel.numReactions])
        self.getCompartmentVolume()

    def getBasicModelProperties(self):
        self.writer.parsedModel.numSpecies = self.sbmlModel.getNumSpecies()
        self.writer.parsedModel.numReactions = self.sbmlModel.getNumReactions()
        self.writer.parsedModel.numGlobalParameters = self.sbmlModel.getNumParameters()

    def getCompartmentVolume(self):
        listOfCompartments = self.sbmlModel.getListOfCompartments()

        for i in range(len(listOfCompartments)):
            if listOfCompartments[i].isSetVolume():
                self.comp += 1
                self.parameterId.append(listOfCompartments[i].getId())
                self.writer.parsedModel.parameterId.append('compartment' + repr(i + 1))
                self.writer.parsedModel.parameter.append(listOfCompartments[i].getVolume())
                self.writer.parsedModel.listOfParameter.append(self.sbmlModel.getCompartment(i))

    def getGlobalParameters(self):
        # Differs between CUDA and Python/C
        for i in range(self.writer.parsedModel.numGlobalParameters):
            self.parameterId.append(self.sbmlModel.getParameter(i).getId())
            self.writer.parsedModel.parameter.append(self.sbmlModel.getParameter(i).getValue())
            self.writer.parsedModel.listOfParameter.append(self.sbmlModel.getParameter(i))

    def getSpecies(self):
        # Differs between CUDA and Python/C
        self.listOfSpecies = self.sbmlModel.getListOfSpecies()

        for k in range(len(self.listOfSpecies)):
            self.writer.parsedModel.species.append(self.listOfSpecies[k])
            self.speciesId.append(self.listOfSpecies[k].getId())
            self.S1.append(0.0)
            self.S2.append(0.0)
            self.reactant.append(0)
            self.product.append(0)
            self.writer.parsedModel.initValues.append(
                self.getSpeciesValue(self.listOfSpecies[k]))  # Only used by the python writer

    def analyseModelStructure(self):
        # Differs between CUDA and Python/C
        reaction = []
        numReactants = []
        numProducts = []

        self.listOfReactions = self.sbmlModel.getListOfReactions()

        for i in range(len(self.listOfReactions)):
            for a in range(len(self.writer.parsedModel.species)):
                self.S1[a] = 0.0
                self.S2[a] = 0.0

            numReactants.append(self.listOfReactions[i].getNumReactants())
            numProducts.append(self.listOfReactions[i].getNumProducts())
            self.writer.parsedModel.kineticLaw.append(self.listOfReactions[i].getKineticLaw().getFormula())
            self.numLocalParameters.append(self.listOfReactions[i].getKineticLaw().getNumParameters())

            for j in range(numReactants[i]):
                self.reactant[j] = self.listOfReactions[i].getReactant(j)

                for k in range(len(self.writer.parsedModel.species)):
                    if self.reactant[j].getSpecies() == self.writer.parsedModel.species[k].getId():
                        self.S1[k] = self.reactant[j].getStoichiometry()

            for l in range(numProducts[i]):
                self.product[l] = self.listOfReactions[i].getProduct(l)

                for k in range(len(self.writer.parsedModel.species)):
                    if self.product[l].getSpecies() == self.writer.parsedModel.species[k].getId():
                        self.S2[k] = self.product[l].getStoichiometry()

            for m in range(len(self.writer.parsedModel.species)):
                self.writer.parsedModel.stoichiometricMatrix[m][i] = -self.S1[m] + self.S2[m]

            for n in range(self.numLocalParameters[i]):
                self.writer.parsedModel.parameter.append(
                    self.listOfReactions[i].getKineticLaw().getParameter(n).getValue())
                self.writer.parsedModel.listOfParameter.append(self.listOfReactions[i].getKineticLaw().getParameter(n))

            for n in range(self.comp):
                compartment_name = self.parameterId[n]
                new_name = 'compartment' + repr(n + 1)
                node = self.sbmlModel.getReaction(i).getKineticLaw().getMath()
                new_node = self.rename(node, compartment_name, new_name)
                self.writer.parsedModel.kineticLaw[i] = formulaToString(new_node)

    def analyseFunctions(self):
        sbmlListOfFunctions = self.sbmlModel.getListOfFunctionDefinitions()

        for fun in range(len(sbmlListOfFunctions)):
            self.writer.parsedModel.listOfFunctions.append(sbmlListOfFunctions[fun])
            self.writer.parsedModel.functionArgument.append([])
            self.writer.parsedModel.functionBody.append(
                formulaToString(self.writer.parsedModel.listOfFunctions[fun].getBody()))

            for funArg in range(self.writer.parsedModel.listOfFunctions[fun].getNumArguments()):
                self.writer.parsedModel.functionArgument[fun].append(
                    formulaToString(self.writer.parsedModel.listOfFunctions[fun].getArgument(funArg)))
                old_name = self.writer.parsedModel.functionArgument[fun][funArg]
                node = self.writer.parsedModel.listOfFunctions[fun].getBody()
                new_node = self.rename(node, old_name, "a" + repr(funArg + 1))
                self.writer.parsedModel.functionBody[fun] = formulaToString(new_node)
                self.writer.parsedModel.functionArgument[fun][funArg] = "a" + repr(funArg + 1)

    def analyseRules(self):
        self.writer.parsedModel.listOfRules = self.sbmlModel.getListOfRules()
        for rule in range(len(self.writer.parsedModel.listOfRules)):
            self.writer.parsedModel.ruleFormula.append(self.writer.parsedModel.listOfRules[rule].getFormula())
            self.writer.parsedModel.ruleVariable.append(self.writer.parsedModel.listOfRules[rule].getVariable())

    def analyseEvents(self):
        self.writer.parsedModel.listOfEvents = self.sbmlModel.getListOfEvents()
        for event in range(len(self.writer.parsedModel.listOfEvents)):
            self.writer.parsedModel.eventCondition.append(
                formulaToString(self.writer.parsedModel.listOfEvents[event].getTrigger().getMath()))
            self.listOfAssignmentRules = self.writer.parsedModel.listOfEvents[event].getListOfEventAssignments()
            self.writer.parsedModel.eventVariable.append([])
            self.writer.parsedModel.eventFormula.append([])

            for rule in range(len(self.listOfAssignmentRules)):
                self.writer.parsedModel.eventVariable[event].append(self.listOfAssignmentRules[rule].getVariable())
                self.writer.parsedModel.eventFormula[event].append(
                    formulaToString(self.listOfAssignmentRules[rule].getMath()))

    def renameEverything(self):

        NAMES = [[], []]
        NAMES[0].append(self.parameterId)
        NAMES[0].append(self.writer.parsedModel.parameterId)
        NAMES[1].append(self.speciesId)
        NAMES[1].append(self.writer.parsedModel.speciesId)

        for nam in range(2):

            for i in range(len(NAMES[nam][0])):
                old_name = NAMES[nam][0][i]
                new_name = NAMES[nam][1][i]

                for k in range(self.writer.parsedModel.numReactions):
                    node = self.sbmlModel.getReaction(k).getKineticLaw().getMath()
                    new_node = self.rename(node, old_name, new_name)
                    self.writer.parsedModel.kineticLaw[k] = formulaToString(new_node)

                for k in range(len(self.writer.parsedModel.listOfRules)):
                    node = self.writer.parsedModel.listOfRules[k].getMath()
                    new_node = self.rename(node, old_name, new_name)
                    self.writer.parsedModel.ruleFormula[k] = formulaToString(new_node)
                    if self.writer.parsedModel.ruleVariable[k] == old_name:
                        self.writer.parsedModel.ruleVariable[k] = new_name

                for k in range(len(self.writer.parsedModel.listOfEvents)):
                    node = self.writer.parsedModel.listOfEvents[k].getTrigger().getMath()
                    new_node = self.rename(node, old_name, new_name)
                    self.writer.parsedModel.eventCondition[k] = formulaToString(new_node)
                    self.listOfAssignmentRules = self.writer.parsedModel.listOfEvents[k].getListOfEventAssignments()

                    for cond in range(len(self.listOfAssignmentRules)):
                        node = self.listOfAssignmentRules[cond].getMath()
                        new_node = self.rename(node, old_name, new_name)
                        self.writer.parsedModel.eventFormula[k][cond] = formulaToString(new_node)
                        if self.writer.parsedModel.eventVariable[k][cond] == old_name:
                            self.writer.parsedModel.eventVariable[k][cond] = new_name

    def rename(self, node, old_name, new_name):
        typ = node.getType()

        if typ == AST_NAME or typ == AST_NAME_TIME:
            nme = node.getName()
            if nme == old_name:
                node.setName(new_name)

        for n in range(node.getNumChildren()):
            self.rename(node.getChild(n), old_name, new_name)
        return node

    def getSpeciesValue(self, specie):
        if specie.isSetInitialAmount() and specie.isSetInitialConcentration():
            return specie.getInitialConcentration()  # The initial values are only used in ODE and SDE solvers so we take the concentration (if it was used in gillespie we would have taken the value)
        if specie.isSetInitialAmount():
            return specie.getInitialAmount()
        else:
            return specie.getInitialConcentration()
