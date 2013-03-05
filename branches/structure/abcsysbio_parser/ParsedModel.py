class ParsedModel:
    def __init__(self):
        self.species = []
        self.speciesId = []
        self.numSpecies = 0
        
        self.listOfParameter = []
        self.parameter = []
        self.parameterId = []
        self.numGlobalParameters = 0
        
        self.listOfFunctions = []
        self.functionArgument = []
        self.functionBody = []
        
        self.kineticLaw = []
        self.numReactions = 0
        
        self.listOfRules = []
        self.ruleFormula = []
        self.ruleVariable = []
        
        self.listOfEvents = []
        self.eventCondition = []
        self.eventVariable = []
        self.eventFormula = []
        
        self.initValues = []#Only used by the Python writer
        self.name = ""
        self.stoichiometricMatrix = None