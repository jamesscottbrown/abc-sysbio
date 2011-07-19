from CUDAParser import CUDAParser

class SdeAndGillespieCUDAParser(CUDAParser):
    def __init__(self, sbmlFileName, modelName, integrationType, method, inputPath="", outputPath=""):
        CUDAParser.__init__(self, sbmlFileName, modelName, integrationType, method, inputPath="", outputPath="")
        
    def renameMathFunctions(self):
        CUDAParser.renameMathFunctions(self)
        self.mathCuda.append('t')
