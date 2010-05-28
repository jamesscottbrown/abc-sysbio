import re
import sys
from numpy import *

####################read the input file
def getAlgorithmInfo(filename):

    """
    Parse the user-provided input file and return all information required to run the abc-SMC algorithm.
    
    ***** args *****
    filename:  a string that describes the source of the user 
               provided input file.
    
    """

    modnum=re.compile('<modelnumber>')
    numout=re.compile('<population size>')
    eps=re.compile('<epsilon\d*>')
    dat=re.compile('<data>')
    datfile=re.compile('<datafile>')
    flowdatfile=re.compile('<flowdatafile>')
    res=re.compile('<restart>\s*True')
    beta=re.compile('<beta>')
    times=re.compile('times:')
    var=re.compile('variable\d+:')
    numvar=re.compile('number of variables:')
    s=re.compile('\/\/\n')
    model=re.compile('model\d+:')
    timestep=re.compile('<dt>')
    rtolRE=re.compile('<rtol>')
    atolRE=re.compile('<atol>')
    constKernels = re.compile('<consant kernels>')
    modelKernel = re.compile('<model kernel>')
    
    plus=re.compile('\+')
    minus=re.compile('\-')
    multiplikation=re.compile('\*')

    n=re.compile('name:')
    q=re.compile('source:')
    t=re.compile('type:')
    init=re.compile('initial values:')

    REtrue=re.compile('[true]|[True]')
    
    #f is the regular expression to find the fit information in the input file
    f=re.compile('fit:')
    #noneRE is the regular expression to find if fit is none.
    #This is not an ideal solution
    noneRE=re.compile('None')
    #commaRE finds , surrounded by zero or more blank spaces
    commaRE=re.compile('\s*,\s*')
    
    modweight=re.compile('<model weights>')
    
    newLine=re.compile('\n')
   # digit=re.compile('\,*\s*(.+?)\s*\,')
    digit=re.compile('-?\d+\.*\d*e?\+?-?\d*') # change the format to identify a digit (int and floats)
    dig=re.compile('\d+')
    
    param=re.compile('parameter\d')
    const=re.compile('constant')
    uni=re.compile('uniform')
    gauss=re.compile('gauss')
    logn=re.compile('lognormal')
    pri=re.compile('prior:')
    kern=re.compile('kernel:')

    commentRE=re.compile('^\#+')
    
    name=[]
    source=[]
    initValues=[]
    integrationType=[]
    prior=[]
    kernel=[]
    epsilon=[]
    modelWeight=[]
    timepoints=[]
    data=[]
    fit = []


######beta
    BETA=1
    try:
        f_info=open(filename)
        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif beta.match(line):
                try:
                    BETA=int(beta.sub('',line))
                except:
                    print "\n Please provide an integer value for beta <beta> in file '"+filename+"'\n"
                    sys.exit()
                break;
        f_info.close()
    except:
        print "\n can not open '"+filename+"'\n"
        sys.exit()


######restart?
    restart=False
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif res.match(line) and REtrue.search(line):
            restart=True
            break;
    f_info.close()

######how many models
    found_modelnumber=False
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif modnum.match(line):
            try:
                modelnumber=int(modnum.sub('',line))
            except:
                print "\n Please provide an integer value for modelnumber <modelnumber> in file '"+filename+"'\n"
                sys.exit()
            found_modelnumber=True
            break;
    f_info.close()
    
    if found_modelnumber==False:
        print "\n No value for modelnumber  <modelnumber> is given in file '"+filename+"'\n"
        sys.exit()


#####model kernel
    modelkernel=0.7
    if modelnumber>1:
        f_info=open(filename)
        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif modelKernel.match(line):
                try:
                    modelkernel=float(modelKernel.sub('',line))
                except:
                    print "\n Please provide an float value for model kernel <model kernel> in file '"+filename+"'\n"
                    sys.exit()
                break;
        f_info.close()




#####get epsilon
    foundEps = -1
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif eps.match(line):
            if foundEps ==(-1):
                foundEps = 0
                epsilon.append([])
                sequence=eps.sub('',line)
                seq=digit.findall(sequence)
                for x in range(0,len(seq)):
                    try:
                        epsilon[foundEps].append(float(seq[x]))
                    except:
                        print "\nMaximum distance values <epsilon> have a wrong format in file '"+filename+"'\n"
                        sys.exit()
            else:
                foundEps = foundEps + 1
                epsilon.append([])
                sequence=eps.sub('',line)
                seq=digit.findall(sequence)
                for x in range(0,len(seq)):
                    try:
                        epsilon[foundEps].append(float(seq[x]))
                    except:
                        print "\nMaximum distance values <epsilon> have a wrong format in file '"+filename+"'\n"
                        sys.exit()


    f_info.close()


#####get dt
    found_dt=False
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif timestep.match(line):
            try:
                dt=float(timestep.sub('',line))
            except:
                print "\nValue for internal timestep <dt> has a wrong format in file '"+filename+"'\n"
                sys.exit()
            found_dt=True
            break
    f_info.close()
    
    if found_dt==False:
        print "\n No value for the internal timestep <dt> is given in file ' "+filename+"'\n"
        sys.exit()


######use constant kernels?
    constKern=False
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif constKernels.match(line) and REtrue.search(line):
            constKern=True
            break;
    f_info.close()


#####get rtol and atol
    rtol=None
    atol=None
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif rtolRE.match(line):
            try:
                rtol=float(rtolRE.sub('',line))
            except:
                print "\nValue for <rtol> has a wrong format in file '"+filename+"'\n"
                sys.exit()
        elif atolRE.match(line):
            try:
                atol=float(atolRE.sub('',line))
            except:
                print "\nValue for <atol> has a wrong format in file '"+filename+"'\n"
                sys.exit()
    f_info.close()


#####get model weights
    weightsFound=False
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif modweight.match(line):
            weightsFound=True
            sequence=modweight.sub('',line)
            seq=digit.findall(sequence)
            for x in range(0,len(seq)):
                try:
                    modelWeight.append(float(seq[x]))
                except:
                    print "\nValues for model weights <modelweights> have a wrong format in file '"+filename+"'\n"
                    sys.exit()
            break;
    f_info.close()

    if weightsFound==False:
        for x in range(0,modelnumber):
            modelWeight.append(1.0)


#######get num Output (population size)
    found_numOutput=False
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif numout.match(line):
            try:
                numOutput=int(numout.sub('',line))
                found_numOut=True
                break;
            except:
                print "\n Please provide an integer value for the population size <population size> in file '"+filename+"'\n"
    f_info.close()

    if numOutput==False:
        print "\n No value for the population size <population size> is given in file '"+filename+"'\n"
        sys.exit()

#########get flow cytometry type data
    f_info=open(filename)

    for line in f_info.readlines():
        if commentRE.match(line): continue
        else:

            if flowdatfile.match(line):
                store = []

                try:
                    datafile = flowdatfile.sub('',line).strip()
                    print 'reading datafile:', datafile
                    
                    tf = open(datafile,'r')

                except:
                    print "\n Please provide a valid file for <flowdatafile> in file '"+filename+"'\n"

                currenttime = None
                data = []
                store_m = []
                for dline in tf:
                    fields = dline.strip().split(' ');

                    tinfile = float(fields[0])

                    if tinfile > currenttime:
                        # add current store_m to data
                        if currenttime != None :
                            data.append( store_m[:] )

                        # new time point
                        timepoints.append( tinfile )
                        store_m = []

                        currenttime = tinfile
                    
                    elif tinfile < currenttime:
                        print "\n Times in flowdata are non increasing \n"
                        sys.exit()
                    
                    store_m.append( fields[1:] )

                # at the end require append of current matrix
                data.append( store_m[:] )

                for i in range(len(timepoints)):
                    for j in range(3):
                        print timepoints[i], data[i][j][:]
                
                sys.exit()
                break;

#########get data
    f_info=open(filename)

    stop=False
    dataFound=False
    dataFromFile = []
    maskFromFile = []
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif(stop==True): break
        else:

            if dat.match(line):
                dataFound=True
        
            if (s.match(line) and dataFound==True):
                stop=True
            
            if (stop==False and dataFound==True):
                if times.match(line):
                    sequence=times.sub('',line)
                    seq=digit.findall(sequence)
                    for x in range(0,len(seq)):
                        timepoints.append(float(seq[x]))
                if var.match(line):
                    sequence=var.sub('',line)

                    # get the data values
                    seq = re.split(',',sequence.strip())
                    
                    # detect NA values
                    mask = [ re.match("\s*NA\s*", seq[i]) != None for i in range(len(seq)) ]
                    
                    # set NA values to zero
                    for i in range(len(seq)):
                        if( mask[i] == True): seq[i] = 0
                    
                    dataFromFile.append(seq)
                    maskFromFile.append(mask)

            if datfile.match(line):
                store = []
                try:
                    datafile = datfile.sub('',line).strip()
                    print 'reading datafile:', datafile
                    
                    tf = open(datafile,'r')

                except:
                    print "\n Please provide a valid file for <datafile> in file '"+filename+"'\n"
                    
                for dline in tf:
                    fields = dline.strip().split(' ');
                    timepoints.append(float(fields[0]))

                    tstore = []
                    for i in range(1,len(fields)):
                        tstore.append(fields[i])

                    store.append( tstore[:] )
                    
                nvars = len(store[0])
                ntimes = len(timepoints)
            
                for nv in range(nvars):
                    species = [ store[nt][nv] for nt in range(ntimes) ]

                    # detect NA values
                    mask = [ re.match("\s*NA\s*", species[nt]) != None for nt in range(ntimes) ]

                    # set NA values to zero
                    for i in range(ntimes):
                        if( mask[i] == True): species[i] = 0

                    dataFromFile.append(species)
                    maskFromFile.append(mask)

                dataFound=True
                break;
                
            
                

    f_info.close()
    data_unmasked = zeros([len(timepoints),len(dataFromFile)])
    data_mask = zeros([len(timepoints),len(dataFromFile)])
    for x in range(0, len(dataFromFile)):
        try:
            data_unmasked[:,x]=dataFromFile[x]
            data_mask[:,x]=maskFromFile[x]
        except:
            print "\n<data> is wrong defined: length of 'times' and 'variables' are not the same!\n"

    data = ma.array(data_unmasked, mask = data_mask)

    #print data
    #for i in range(len(timepoints)):
    #    print timepoints[i], data[i,:].transpose()
    #sys.exit()

######sample from prior
    sampleFromPrior=True


####get model properties
    for i in range(0,modelnumber):
    
        initValues.append([])
        stop=False
        rightModel=False
        model=re.compile('<model%s>' %(i+1))
        f_info=open(filename)
        paramFound=False
        k=0
        prior0=[]
        kernel0=[]

        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif stop==True:
                break
            else:
                if s.match(line) and rightModel==True:
                    stop=True
                if stop==False and rightModel==True:
                
                    if n.match(line): 
                        sequence=n.sub('',line)
                        name.append(sequence.strip())
                    if q.match(line): 
                        sequence=q.sub('',line)
                        source.append(sequence.strip())
                    if t.match(line): 
                        sequence=t.sub('',line)
                        integrationType.append(sequence.strip())            
                
                    if init.match(line):
                        sequence=init.sub('',line)
                        seq=digit.findall(sequence)
                        for x in range(0,len(seq)):
                            try:
                                initValues[i].append(float(seq[x]))
                            except:
                                print "\nInitial values for model "+name[i]+" have a wrong format in file '"+filename+"'\n"
                                sys.exit()

                    
                    if paramFound==True and pri.search(line):
                        sequence=pri.sub('',line)
                        prior0.append([0,0,0])
                        kernel0.append([0,0,0])
                        if const.search(sequence):
                            prior0[k][0]=0
                            seq0=const.sub('',sequence)
                            seq=digit.findall(seq0)
                            try:
                                prior0[k][1]=float(seq[0])
                            except:
                                print "\nValue of the prior for model "+name[i]+" have a wrong format in file '"+filename+"'\n"
                                sys.exit()
                        if uni.search(sequence): 
                            prior0[k][0]=2
                            seq0=uni.sub('',sequence)
                            seq=digit.findall(seq0)
                            try:
                                prior0[k][1]=float(seq[0])
                                prior0[k][2]=float(seq[1])
                            except:
                                print "\nValue of the prior for model "+name[i]+" have a wrong format in file '"+filename+"'\n"
                                sys.exit()
                        if logn.search(sequence): 
                            prior0[k][0]=3
                            seq0=logn.sub('',sequence)
                            seq=digit.findall(seq0)
                            prior0[k][1]=float(seq[0])
                            prior0[k][2]=float(seq[1])
                        if prior0[k][0]==0: 
                            paramFound=False
                            k=k+1
                
                    if paramFound==True and kern.search(line):
                        sequence=kern.sub('',line)
                        if uni.search(sequence): 
                            kernel0[k][0]=1
                            seq=digit.findall(sequence)
                            kernel0[k][1]=float(seq[0])
                            kernel0[k][2]=float(seq[1])
                        if gauss.search(sequence): 
                            kernel0[k][0]=2
                            seq=digit.findall(sequence)
                            kernel0[k][1]=float(seq[0])
                            kernel0[k][2]=float(seq[1])
                        k=k+1
                        paramFound=False
                          
                    if param.match(line):
                        paramFound=True
                
                if model.match(line):
                    rightModel=True

        prior.append(prior0)
        kernel.append(kernel0)

        f_info.close()

#######get fitting rules
    speciesIntRE = re.compile('species\d+\n|species\d+\s+|species\d+\+|species\d+\-|species\d+\*|species\d+$')
    speciesRE = re.compile('species')
    operatorRE = re.compile('\+|\-|\*')
    speciesZeroRE = re.compile('species0')
    f_info=open(filename)

    for line in f_info.readlines():
        ##
        ##The fitting rules need to be written differently
        ##if samplePoints is one-dimensional
        ##or two-dimensional
        ##I really need a two-dimensional model
        ##to play with
        ##This information could be obtained from the length of "Initial Values"
        ##
        index = 0
        if commentRE.match(line): continue
        elif f.match(line):
            line = f.sub('', line)
            if noneRE.search(line):
                fit.append(None)
                index = index + 1
            else:
                fit.append([])
                dimensionsOfModel = len(initValues[index])
                fit[-1]=line.strip().split(',')
                for i in range(0, len(fit[-1])):
                    fit[-1][i] = fit[-1][i].strip()
                    if speciesZeroRE.search(fit[-1][i]):
                        print "\nIncorrect species indexing!\n"
                        sys.exit()
                    if dimensionsOfModel == 1:
                        fit[-1][i]=speciesIntRE.sub('samplePoints', fit[-1][i])
                    else:
                        for m in speciesIntRE.finditer(fit[-1][i]):
                            mySpecies = str(m.group())
                            mySpecies = operatorRE.sub('', mySpecies)                          
                            mySpeciesIndex = speciesRE.sub('', mySpecies)
                            mySpecies = re.compile(mySpecies)
                            mySpeciesIndex = int(mySpeciesIndex)
                            mySpeciesIndex = mySpeciesIndex-1
                            mySpeciesIndex = str(mySpeciesIndex)
                            fit[-1][i]=mySpecies.sub('samplePoints[:,'+mySpeciesIndex+']', fit[-1][i])

    f_info.close()
        

    return restart,name,data,timepoints,numOutput,epsilon,initValues,integrationType,modelWeight,prior,kernel,sampleFromPrior,source,fit,BETA,dt, rtol, atol,constKern,modelkernel




def getInfoForSimulation(filename):

    """
    Parse the user-provided input file and return all information required to simulate the model.
    
    ***** args *****
    filename:  a string that describes the source of the user 
               provided input file.
    
    """
    
    modnum=re.compile('<modelnumber>')
    numout=re.compile('<population size>')
    dat=re.compile('<data>')
    times=re.compile('times:')
    s=re.compile('\/\/\n')
    model=re.compile('model\d+:')
    timestep=re.compile('<dt>')
    rtolRE=re.compile('<rtol>')
    atolRE=re.compile('<atol>')

    
    plus=re.compile('\+')
    minus=re.compile('\-')
    multiplikation=re.compile('\*')

    n=re.compile('name:')
    q=re.compile('source:')
    t=re.compile('type:')
    init=re.compile('initial values:')

    REtrue=re.compile('[true]|[True]')
    
    #f is the regular expression to find the fit information in the input file
    f=re.compile('fit:')
    #noneRE is the regular expression to find if fit is none.
    #This is not an ideal solution
    noneRE=re.compile('None')
    #commaRE finds , surrounded by zero or more blank spaces
    commaRE=re.compile('\s*,\s*')
    
   
    newLine=re.compile('\n')
    digit=re.compile('-?\d+\.*\d*e?\+?-?\d*') # change the format to identify a digit (int and floats)
    dig=re.compile('\d+')
    
    param=re.compile('parameter\d')
    const=re.compile('constant')
    uni=re.compile('uniform')
    gauss=re.compile('gauss')
    logn=re.compile('lognormal')
    pri=re.compile('prior:')
    kern=re.compile('kernel:')

    name=[]
    source=[]
    initValues=[]
    integrationType=[]
    parametersConstant=[]
    timepoints=[]
    fit = []

    commentRE=re.compile('^\#+')

    s=re.compile('\/\/\n')

#########get timepoints
    f_info=open(filename)
    dataFound=False

    stop=False
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif(stop==True): break
        else:
            if (s.match(line) and dataFound==True):
                stop=True
            
            if (stop==False and dataFound==True):
                if times.match(line):
                    sequence=times.sub('',line)
                    seq=digit.findall(sequence)
                    for x in range(0,len(seq)):
                        timepoints.append(float(seq[x]))
            if dat.match(line):
                dataFound=True
        
    f_info.close()
 

######how many models
    modnum=re.compile('<modelnumber>')
    found_modelnumber=False
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif modnum.match(line):
            try:
                modelnumber=int(modnum.sub('',line))
            except:
                print "\n Please provide an integer value for modelnumber <modelnumber> in file '"+filename+"'\n"
                sys.exit()
            found_modelnumber=True
            break;
    f_info.close()
    
    if found_modelnumber==False:
        print "\n No value for modelnumber  <modelnumber> is given in file '"+filename+"'\n"
        sys.exit()

    

    ####get model properties
    for i in range(0,modelnumber):
    
        initValues.append([])
        stop=False
        rightModel=False
        model=re.compile('<model%s>' %(i+1))
        f_info=open(filename)
        paramFound=False
        k=0
        parameter0=[]

        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif stop==True:
                break
            else:
                if s.match(line) and rightModel==True:
                    stop=True
                if stop==False and rightModel==True:
                
                    if n.match(line): 
                        sequence=n.sub('',line)
                        name.append(sequence.strip())
                    if q.match(line): 
                        sequence=q.sub('',line)
                        source.append(sequence.strip())    
                    if t.match(line): 
                        sequence=t.sub('',line)
                        integrationType.append(sequence.strip())                          
                    if init.match(line):
                        sequence=init.sub('',line)
                        seq=digit.findall(sequence)
                        for x in range(0,len(seq)):
                            try:
                                initValues[i].append(float(seq[x]))
                            except:
                                print "\nInitial values for model "+name[i]+" have a wrong format in file '"+filename+"'\n"
                                sys.exit()
                                
                    if paramFound==True and pri.search(line):
                        sequence=pri.sub('',line)
                        if const.search(sequence):
                            seq0=const.sub('',sequence)
                            seq=digit.findall(seq0)
                            try:
                                parameter0.append(float(seq[0]))
                            except:
                                print "\nValue of the prior for model "+name[i]+" have a wrong format in file '"+filename+"'\n"
                                sys.exit()
                        if uni.search(sequence):
                            print "\nYou can only supply constant priors for a simulation.\n"
                            sys.exit()
                        if logn.search(sequence): 
                            print "\nYou can only supply constant priors for a simulation.\n"
                            sys.exit()
##                        if prior0[k][0]==0: 
##                            paramFound=False
##                            k=k+1
            
                    if param.match(line):
                        paramFound=True
                
                if model.match(line):
                    rightModel=True

        parametersConstant.append(parameter0)

        f_info.close()

    

    

    #####get dt
    timestep=re.compile('<dt>')
    found_dt=False
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif timestep.match(line):
            try:
                dt=float(timestep.sub('',line))
            except:
                print "\nValue for internal timestep <dt> has a wrong format in file '"+filename+"'\n"
                sys.exit()
            found_dt=True
            break
    f_info.close()

    #####get rtol and atol
    rtol=None
    atol=None
    f_info=open(filename)
    for line in f_info.readlines():
        if commentRE.match(line): continue
        elif rtolRE.match(line):
            try:
                rtol=float(rtolRE.sub('',line))
            except:
                print "\nValue for <rtol> has a wrong format in file '"+filename+"'\n"
                sys.exit()
        elif atolRE.match(line):
            try:
                atol=float(atolRE.sub('',line))
            except:
                print "\nValue for <atol> has a wrong format in file '"+filename+"'\n"
                sys.exit()
    f_info.close()

    return name, timepoints, initValues, integrationType, parametersConstant, source, dt, rtol, atol
