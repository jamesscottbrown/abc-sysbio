from ODEPythonWriter import ODEPythonWriter
from GillespiePythonWriter import GillespiePythonWriter
from SDEPythonWriter import SDEPythonWriter
from ODECUDAWriter import OdeCUDAWriter
from SDECUDAWriter import SdeCUDAWriter
from GillespieCUDAWriter import GillespieCUDAWriter
#from CWriter import CWriter
from CandPythonParser import CandPythonParser
from SDEAndGillespieCUDAParser import SdeAndGillespieCUDAParser
from ODECUDAParser import OdeCUDAParser
import re

def ParseAndWrite(source, integrationType, modelName = None, inputPath = "", outputPath = "", method = None):
    
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
    modelName:
                  a list of strings.
                  modelName describes the names of the parsed model files.

    method:
                  an integer number.
                  Type of noise in a stochastic system.
                  (Only implemented for stochastic systems solved with sdeint.)
                  Possible options are:
                  1 --- default
                  2 --- Ornstein-Uhlenbeck
                  3 --- geometric Brownian motion

    """

     #regular expressions for detecting integration types and integration language
    c=re.compile('C', re.IGNORECASE)
    py=re.compile('Python', re.I)
    cuda=re.compile('CUDA', re.I)
    
    ode=re.compile('ODE', re.I)
    sde=re.compile('SDE', re.I)
    heun=re.compile('Heun', re.I)
    milstein=re.compile('Milstein', re.I)
    gil = re.compile('Gillespie', re.I)

    #check that you have appropriate lengths of integration types and sources
    #(need equal lengths)
    if(not(len(source)==len(integrationType))):
        print "\nError: Number of sources is not the same as number of integrationTypes!\n"
    #check that you have model names,
    #if not the models will be named model1, model2, etc
    else:
        if(modelName==None):
            modelName=[]
            for x in range(0,len(source)):
                modelName.append("model"+repr(x+1))
        else:
            for x in range(0,len(modelName)):
                if(modelName[x]==""):
                    modelName[x]="model"+repr(x+1)

        #if no method is specified and the integrationType is "SDE"
        #the method type defaults to 1
        for model in range(0,len(source)):
            if cuda.search(integrationType[model]):
                if(not(sde.search(integrationType[model]) or gil.search(integrationType[model]) or ode.search(integrationType[model]))):
                    print "\nError: an integration type is required for CUDA"
                elif (sde.search(integrationType[model])):
                    if(heun.search(integrationType[model]) or milstein.search(integrationType[model])):
                        print "\nError: Only Euler is available in Cuda"
                    else:
                        if(method==None or method[model]==""):
                            parser = SdeAndGillespieCUDAParser(source[model], modelName[model], "CUDA SDE", 1, inputPath, outputPath)
                        else:
                            parser = SdeAndGillespieCUDAParser(source[model], modelName[model], "CUDA SDE", method[model], inputPath, outputPath)
                elif(gil.search(integrationType[model])):
                    parser = SdeAndGillespieCUDAParser(source[model], modelName[model], integrationType[model], None, inputPath, outputPath)
                else:
                    parser = OdeCUDAParser(source[model], modelName[model], integrationType[model], None, inputPath, outputPath)
                    
            elif c.search(integrationType[model]):
                if (sde.search(integrationType[model])):
                    if (not (method==None or method==1)):
                        print "\nError: Only the method 1 of SDE resolution can be used in C"
                    else:
                        parser = CandPythonParser(source[model],modelName[model], "C", None, inputPath, outputPath)
                else:
                    parser = CandPythonParser(source[model],modelName[model], "C", None, inputPath, outputPath)
            
            elif py.search(integrationType[model]):
                if(integrationType==None):
                    print "\nError: an integration type is required for Python"
                elif (sde.search(integrationType[model])):
                    if(heun.search(integrationType[model]) or milstein.search(integrationType[model])):
                        print "\nError: Only Euler is available in Python"
                    else:
                        if(method==None or method[model]==""):
			    parser = CandPythonParser(source[model], modelName[model], "Python SDE", 1, inputPath, outputPath)
                        else:
                            parser = CandPythonParser(source[model], modelName[model], "Python SDE", method[model], inputPath, outputPath)
                else:
                    parser = CandPythonParser(source[model], modelName[model], integrationType[model], None, inputPath, outputPath)
           
        
                        
                        
                        
                        
                        
                        
                        
                        
                        
