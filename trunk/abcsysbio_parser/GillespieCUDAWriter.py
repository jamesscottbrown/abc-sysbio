from libsbml import *
from numpy import *
from abcsysbio.relations import *
import os
import re
from Writer import Writer

## replace the species and parameters recursively
##
## replace
## pq = re.compile(speciesId[q])
## string=pq.sub('y['+repr(q)+']' ,string)
## with
## string = rep(string, speciesId[q],'y['+repr(q)+']')

def rep(str,find,replace):

    ex = find+"[^0-9]"
    ss = str;
    while re.search(ex,ss) != None:
        res = re.search(ex,ss)
        ss = ss[0:res.start()] + replace + " " + ss[res.end()-1:]

    ex = find+"$"
    if re.search(ex,ss) != None:
        res = re.search(ex,ss)
        ss = ss[0:res.start()] + replace + " " + ss[res.end():]
 
    return ss;

class GillespieCUDAWriter(Writer):
    def __init__(self, sbmlFileName, modelName="", inputPath="", outputPath=""):
       Writer.__init__(self, sbmlFileName, modelName, inputPath, outputPath)        
       self.out_file=open(os.path.join(outputPath,self.parsedModel.name+".cu"),"w")
        
    def mathMLConditionParserCuda(self, mathMLstring):
        """
        Replaces and and or with and_ and or_ in a MathML string.
        Returns the string with and and or replaced by and_ and or_
    
        ***** args *****
    
        mathMLstring:
    
                A mathMLstring
    
        """
        
        andString = re.compile("and")
        orString = re.compile("or")
        mathMLstring = andString.sub("and_", mathMLstring)
        mathMLstring = orString.sub("or_", mathMLstring)

        return mathMLstring
    
    def write(self):
        p=re.compile('\s')
        
        self.out_file.write("#define NSPECIES " + str(self.parsedModel.numSpecies) + "\n")
        self.out_file.write("#define NPARAM " + str(self.parsedModel.numGlobalParameters) + "\n")
        self.out_file.write("#define NREACT " + str(self.parsedModel.numReactions) + "\n")
        self.out_file.write("\n")
    
        numEvents = len(self.parsedModel.listOfEvents)
        numRules = len(self.parsedModel.listOfRules)
        num = numEvents+numRules
        if num>0:
            self.out_file.write("#define leq(a,b) a<=b\n")
            self.out_file.write("#define neq(a,b) a!=b\n")
            self.out_file.write("#define geq(a,b) a>=b\n")
            self.out_file.write("#define lt(a,b) a<b\n")
            self.out_file.write("#define gt(a,b) a>b\n")
            self.out_file.write("#define eq(a,b) a==b\n")
            self.out_file.write("#define and_(a,b) a&&b\n")
            self.out_file.write("#define or_(a,b) a||b\n")
        
    
        for i in range(0,len(self.parsedModel.listOfFunctions)):
            self.out_file.write("__device__ float "+self.parsedModel.listOfFunctions[i].getId()+"(")
            for j in range(0, self.parsedModel.listOfFunctions[i].getNumArguments()):
                self.out_file.write("float "+self.parsedModel.functionArgument[i][j])
                if(j<( self.parsedModel.listOfFunctions[i].getNumArguments()-1)):
                    self.out_file.write(",")
            self.out_file.write("){\n    return ")
            self.out_file.write(self.parsedModel.functionBody[i])
            self.out_file.write(";\n}\n")
            self.out_file.write("")
    
        self.out_file.write("\n\n__constant__ int smatrix[]={\n")
        for i in range(0,len(self.parsedModel.stoichiometricMatrix[0])):
            for j in range(0,len(self.parsedModel.stoichiometricMatrix)):
                self.out_file.write("    "+repr(self.parsedModel.stoichiometricMatrix[j][i]))
                if (not(i==(len(self.parsedModel.stoichiometricMatrix)-1) and (j==(len(self.parsedModel.stoichiometricMatrix[0])-1)))):
                    self.out_file.write(",")
            self.out_file.write("\n")
    
    
        self.out_file.write("};\n\n\n")
    
        
        #stoichiometry function moved to Gillespie.py
        #self.out_file.write("__device__ void stoichiometry(int *y, int r, int tid){\n")
        #self.out_file.write("    for(int i=0; i<"+repr(len(species))+"; i++){\n        y[i]+=smatrix[r*"+repr(len(species))+"+ i];\n    }\n}\n\n\n")
        
        self.out_file.write("__device__ void hazards(int *y, float *h, float t, int tid){")
        # wirte rules and events 
        for i in range(0,len(self.parsedModel.listOfRules)):
            if self.parsedModel.listOfRules[i].isRate() == True:
                self.out_file.write("    ")
                if not(self.parsedModel.ruleVariable[i] in speciesId):
                    self.out_file.write(self.parsedModel.ruleVariable[i])
                else:
                    string = "y["+repr(self.parsedModel.speciesId.index(self.parsedModel.ruleVariable[i]))+"]"
                    self.out_file.write(string)
                self.out_file.write("=")
    
                string = self.parsedModel.ruleFormula[i]
                for q in range(0,len(self.parsedModel.speciesId)):
                    #pq = re.compile(self.parsedModel.speciesId[q])
                    #string=pq.sub('y['+repr(q)+']' ,string)
                    string = rep(string, self.parsedModel.speciesId[q],'y['+repr(q)+']')
                for q in range(0,len(self.parsedModel.parameterId)):
                    if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                        flag = False
                        for r in range(0,len(self.parsedModel.eventVariable)):
                            if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                flag = True
                        if flag==False:
                            pq = re.compile(self.parsedModel.parameterId[q])
                            string=pq.sub('tex2D(param_tex,'+repr(q)+',tid)' ,string)
    
                self.out_file.write(string)
                self.out_file.write(";\n")
                
        for i in range(0,len(self.parsedModel.listOfEvents)):
            self.out_file.write("    if( ")
            self.out_file.write(self.mathMLConditionParserCuda(self.parsedModel.eventCondition[i]))
            self.out_file.write("){\n")
            listOfAssignmentRules = self.parsedModel.listOfEvents[i].getListOfEventAssignments()
            for j in range(0, len(listOfAssignmentRules)):
                self.out_file.write("        ")
                #self.out_file.write("float ")
                if not(self.parsedModel.eventVariable[i][j] in self.parsedModel.speciesId):
                    self.out_file.write(self.parsedModel.eventVariable[i][j])
                else:
                    string = "y["+repr(self.parsedModel.speciesId.index(self.parsedModel.eventVariable[i][j]))+"]"
                    self.out_file.write(string)           
                self.out_file.write("=")
                
                string = self.parsedModel.eventFormula[i][j]
                for q in range(0,len(self.parsedModel.speciesId)):
                    #pq = re.compile(self.parsedModel.speciesId[q])
                    #string=pq.sub('y['+repr(q)+']' ,string)
                    string = rep(string, self.parsedModel.speciesId[q],'y['+repr(q)+']')
                for q in range(0,len(self.parsedModel.parameterId)):
                    if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                        flag = False
                        for r in range(0,len(self.parsedModel.eventVariable)):
                            if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                flag = True
                        if flag==False:
                            pq = re.compile(self.parsedModel.parameterId[q])
                            string=pq.sub('tex2D(param_tex,'+repr(q)+',tid)' ,string)
    
                self.out_file.write(string)
                self.out_file.write(";\n")
            self.out_file.write("    }\n")
    
        self.out_file.write("\n")
    
        for i in range(0, len(self.parsedModel.listOfRules)):
            if self.parsedModel.listOfRules[i].isAssignment():
                self.out_file.write("    ")
                if not(self.parsedModel.ruleVariable[i] in self.parsedModel.speciesId):
                    self.out_file.write("float ")
                    self.out_file.write(self.parsedModel.ruleVariable[i])
                else:
                    string = "y["+repr(self.parsedModel.speciesId.index(self.parsedModel.ruleVariable[i]))+"]"
                    self.out_file.write(string)
                self.out_file.write("=")
     
                string = self.mathMLConditionParserCuda(self.parsedModel.ruleFormula[i])
                for q in range(0,len(self.parsedModel.speciesId)):
                    pq = re.compile(self.parsedModel.speciesId[q])
                    string=pq.sub("y["+repr(q)+"]" ,string)
                for q in range(0,len(self.parsedModel.parameterId)):
                    if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                        flag = False
                        for r in range(0,len(self.parsedModel.eventVariable)):
                            if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                flag = True
                        if flag==False:
                            pq = re.compile(self.parsedModel.parameterId[q])
                            x = "tex2D(param_tex,"+repr(q)+",tid)"
                            string=pq.sub(x,string)
                self.out_file.write(string)
                self.out_file.write(";\n")
        self.out_file.write("\n")
    
    
    
        for i in range(0,self.parsedModel.numReactions):
            self.out_file.write("    h["+repr(i)+"] = ")
    
            string = self.parsedModel.kineticLaw[i]
            for q in range(0,len(self.parsedModel.speciesId)):
                #pq = re.compile(self.parsedModel.speciesId[q])
                #string=pq.sub('y['+repr(q)+']' ,string)
                string = rep(string, self.parsedModel.speciesId[q],'y['+repr(q)+']')
            for q in range(0,len(self.parsedModel.parameterId)):
                if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                    flag = False
                    for r in range(0,len(self.parsedModel.eventVariable)):
                        if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                            flag = True
                    if flag==False:
                        pq = re.compile(self.parsedModel.parameterId[q])
                        string=pq.sub('tex2D(param_tex,'+repr(q)+',tid)' ,string)
       
            string=p.sub('',string)
            self.out_file.write(string+";\n")
            
        self.out_file.write("\n")
        self.out_file.write("}\n\n")
