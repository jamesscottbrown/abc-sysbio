from ODEPythonWriter import ODEPythonWriter
from GillespiePythonWriter import GillespiePythonWriter
from SDEPythonWriter import SDEPythonWriter
from ODECUDAWriter import OdeCUDAWriter
from SDECUDAWriter import SdeCUDAWriter
from GillespieCUDAWriter import GillespieCUDAWriter
from CWriter import CWriter
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
                  Gillespie   ---   for stochastic systems; solved with GillespieAlgorithm (abc)

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

   
    # create lists from the matrices that now exist in order to get the input in the same form as it used to be before submodels were incorporated
    integrationTypelist = []
    sourcelist = []

    for i in range(0, len(source)):
        for j in range(0, len(source[i])):
            sourcelist.append(source[i][j])
            integrationTypelist.append(integrationType[i])



    if (modelName != None):
        namelist = []
        for i in range(0, len(source)):
            for j in range(0, len(source[i])):
                namelist.append(modelName[i][j])
        nmodel = len(modelName)
        nsubmodel = len(namelist)/len(modelName)



    #check that you have appropriate lengths of integration types and sources
    #(need equal lengths)
    if(not(len(source)==len(integrationType))):
        print "\nError: Number of models is not the same as number of integrationTypes!\n" 
    
    #check that there are model names,
    #if not the models will be named model1, model2, etc
    else:
        nmodel = len(integrationType)
        nsubmodel = len(sourcelist)/len(integrationType)
        if(modelName==None):
            namelist = []
            for x in range(0,nmodel):
                for j in range(0,nsubmodel):
                    namelist.append("model"+repr(x+1)+"_"+repr(j+1))

        else:                         
            for x in range(0,nmodel):
                for j in range(0,nsubmodel):
                    if(namelist[x*nsubmodel+j]==""):
                        namelist[x*nsubmodel+j]="model"+repr(x+1)+"_"+repr(j+1)




        #if no method is specified and the integrationType is "SDE"
        #the method type defaults to 1

        for model in range(0,len(sourcelist)):
            if cuda.search(integrationTypelist[model]):
                if(not(sde.search(integrationTypelist[model]) or gil.search(integrationTypelist[model]) or ode.search(integrationTypelist[model]))):
                    print "\nError: an integration type is required for CUDA"
                elif (sde.search(integrationTypelist[model])):
                    if(heun.search(integrationTypelist[model]) or milstein.search(integrationTypelist[model])):
                        print "\nError: Only Euler is available in Cuda"
                    else:
                        if(method==None or method[model]==""):
                            parser = SdeAndGillespieCUDAParser(sourcelist[model], namelist[model], "CUDA SDE", 1, inputPath, outputPath)
                        else:
                            parser = SdeAndGillespieCUDAParser(sourcelist[model], namelist[model], "CUDA SDE", method[model], inputPath, outputPath)
                elif(gil.search(integrationTypelist[model])):
                    parser = SdeAndGillespieCUDAParser(sourcelist[model], namelist[model], integrationTypelist[model], None, inputPath, outputPath)
                else:
                    parser = OdeCUDAParser(sourcelist[model], namelist[model], integrationTypelist[model], None, inputPath, outputPath)



                    
            elif c.search(integrationTypelist[model]):
                if (sde.search(integrationTypelist[model])):
                    if (not (method==None or method==1)):
                        print "\nError: Only the method 1 of SDE resolution can be used in C"
                    else:
                        parser = CandPythonParser(sourcelist[model],namelist[model], "C", None, inputPath, outputPath)
                else:
                    parser = CandPythonParser(sourcelist[model],namelist[model], "C", None, inputPath, outputPath)
 


           
            elif py.search(integrationTypelist[model]):
                if(integrationTypelist==None):
                    print "\nError: an integration type is required for Python"
                elif (sde.search(integrationTypelist[model])):
                    if(heun.search(integrationTypelist[model]) or milstein.search(integrationTypelist[model])):
                        print "\nError: Only Euler is available in Python"
                    else:
                        if(method==None or method[model]==""):
                            parser = CandPythonParser(sourcelist[model], namelist[model], "Python SDE", 1, inputPath, outputPath)
                        else:
                            parser = CandPythonParser(sourcelist[model], namelist[model], "Python SDE", method[model], inputPath, outputPath)
                else:
                    parser = CandPythonParser(sourcelist[model], namelist[model], integrationTypelist[model], None, inputPath, outputPath)
           


                        
                        
                        
                        
                        
                        
                        
