# Algorithm information
import re, sys, numpy

class algorithm_info:
    """
    A class to parse the user-provided input file and return all information required to run the abc-SMC algorithm. 

    """

    def __init__(self, filename):
        #print filename

        self.name=[]
        self.source=[]
        self.initValues=[]
        self.integrationType=[]
        self.prior=[]
        self.kernel=[]
        self.epsilon=[]
        self.modelWeight=[]
        self.timepoints=[]
        self.data=[]
        self.fit = []
        self.modelnumber = 0
        self.numOutput = 0
        self.dt = 0
        self.BETA=1
        self.restart=False
        self.modelkernel=0.7
        self.constKern=False
        self.rtol=None
        self.atol=None
        self.sampleFromPrior=False

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
    
        # f is the regular expression to find the fit information in the input file
        f=re.compile('fit:')
        # noneRE is the regular expression to find if fit is none.
        # This is not an ideal solution
        noneRE=re.compile('None')
        # commaRE finds , surrounded by zero or more blank spaces
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
    
        # beta
        try:
            f_info=open(filename)
            for line in f_info.readlines():
                if commentRE.match(line): continue
                elif beta.match(line):
                    try:
                        self.BETA=int(beta.sub('',line))
                    except:
                        print "\n Please provide an integer value for beta <beta> in file '"+filename+"'\n"
                        sys.exit()
                    break;
            f_info.close()
        except:
            print "\n can not open '"+filename+"'\n"
            sys.exit()


        # restart?
        f_info=open(filename)
        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif res.match(line) and REtrue.search(line):
                self.restart=True
                break;
        f_info.close()

        # how many models
        found_modelnumber=False
        f_info=open(filename)
        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif modnum.match(line):
                try:
                    self.modelnumber=int(modnum.sub('',line))
                except:
                    print "\n Please provide an integer value for modelnumber <modelnumber> in file '"+filename+"'\n"
                    sys.exit()
                found_modelnumber=True
                break;
        f_info.close()
    
        if found_modelnumber==False:
            print "\n No value for modelnumber  <modelnumber> is given in file '"+filename+"'\n"
            sys.exit()


        # model kernel
        if self.modelnumber>1:
            f_info=open(filename)
            for line in f_info.readlines():
                if commentRE.match(line): continue
                elif modelKernel.match(line):
                    try:
                        self.modelkernel=float(modelKernel.sub('',line))
                    except:
                        print "\n Please provide an float value for model kernel <model kernel> in file '"+filename+"'\n"
                        sys.exit()
                    break;
            f_info.close()

        # get epsilon
        foundEps = -1
        f_info=open(filename)
        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif eps.match(line):
                if foundEps ==(-1):
                    foundEps = 0
                    self.epsilon.append([])
                    sequence=eps.sub('',line)
                    seq=digit.findall(sequence)
                    for x in range(0,len(seq)):
                        try:
                            self.epsilon[foundEps].append(float(seq[x]))
                        except:
                            print "\nMaximum distance values <epsilon> have a wrong format in file '"+filename+"'\n"
                            sys.exit()
                else:
                    foundEps = foundEps + 1
                    self.epsilon.append([])
                    sequence=eps.sub('',line)
                    seq=digit.findall(sequence)
                    for x in range(0,len(seq)):
                        try:
                            self.epsilon[foundEps].append(float(seq[x]))
                        except:
                            print "\nMaximum distance values <epsilon> have a wrong format in file '"+filename+"'\n"
                            sys.exit()


        f_info.close()


        # get dt
        found_dt=False
        f_info=open(filename)
        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif timestep.match(line):
                try:
                    self.dt=float(timestep.sub('',line))
                except:
                    print "\nValue for internal timestep <dt> has a wrong format in file '"+filename+"'\n"
                    sys.exit()
                found_dt=True
                break
        f_info.close()
    
        if found_dt==False:
            print "\n No value for the internal timestep <dt> is given in file ' "+filename+"'\n"
            sys.exit()


        # use constant kernels?
        f_info=open(filename)
        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif constKernels.match(line) and REtrue.search(line):
                self.constKern=True
                break;
        f_info.close()


        # get rtol and atol
        f_info=open(filename)
        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif rtolRE.match(line):
                try:
                    self.rtol=float(rtolRE.sub('',line))
                except:
                    print "\nValue for <rtol> has a wrong format in file '"+filename+"'\n"
                    sys.exit()
            elif atolRE.match(line):
                try:
                    self.atol=float(atolRE.sub('',line))
                except:
                    print "\nValue for <atol> has a wrong format in file '"+filename+"'\n"
                    sys.exit()
        f_info.close()


        # get model weights
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
                        self.modelWeight.append(float(seq[x]))
                    except:
                        print "\nValues for model weights <modelweights> have a wrong format in file '"+filename+"'\n"
                        sys.exit()
                break;
        f_info.close()

        if weightsFound==False:
            for x in range(0,self.modelnumber):
                self.modelWeight.append(1.0/float(self.modelnumber))


        # get num Output (population size)
        found_numOutput=False
        f_info=open(filename)
        for line in f_info.readlines():
            if commentRE.match(line): continue
            elif numout.match(line):
                try:
                    self.numOutput=int(numout.sub('',line))
                    found_numOut=True
                    break;
                except:
                    print "\n Please provide an integer value for the population size <population size> in file '"+filename+"'\n"
        f_info.close()

        if self.numOutput==False:
            print "\n No value for the population size <population size> is given in file '"+filename+"'\n"
            sys.exit()

        # get flow cytometry type data
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
                    self.data = []
                    store_m = []
                    for dline in tf:
                        fields = dline.strip().split(' ');

                        tinfile = float(fields[0])
                    
                        if tinfile > currenttime:
                            # add current store_m to data
                            if currenttime != None :
                                self.data.append( store_m[:] )

                            # new time point
                            self.timepoints.append( tinfile )
                            store_m = []

                            currenttime = tinfile
                    
                        elif tinfile < currenttime:
                            print "\n Times in flowdata are non increasing \n"
                            sys.exit()
                    
                        store_m.append( fields[1:] )

                    # at the end require append of current matrix
                    self.data.append( store_m[:] )

                    for i in range(len(self.timepoints)):
                        for j in range(3):
                            print self.timepoints[i], self.data[i][j][:]
                
                    sys.exit()
                    break;

        # get data
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
                            self.timepoints.append(float(seq[x]))
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
                        self.timepoints.append(float(fields[0]))

                        tstore = []
                        for i in range(1,len(fields)):
                            tstore.append(fields[i])

                        store.append( tstore[:] )
                    
                    nvars = len(store[0])
                    ntimes = len(self.timepoints)
            
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
        data_unmasked = numpy.zeros([len(self.timepoints),len(dataFromFile)])
        data_mask = numpy.zeros([len(self.timepoints),len(dataFromFile)])
        for x in range(0, len(dataFromFile)):
            try:
                data_unmasked[:,x]=dataFromFile[x]
                data_mask[:,x]=maskFromFile[x]
            except:
                print "\n<data> is wrong defined: length of 'times' and 'variables' are not the same!\n"

        self.data = numpy.ma.array(data_unmasked, mask = data_mask)

        # print data
        # for i in range(len(self.timepoints)):
        #    print self.timepoints[i], data[i,:].transpose()
        # sys.exit()

        # get model properties
        for i in range(0,self.modelnumber):
    
            self.initValues.append([])
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
                            self.name.append(sequence.strip())
                        if q.match(line): 
                            sequence=q.sub('',line)
                            self.source.append(sequence.strip())
                        if t.match(line): 
                            sequence=t.sub('',line)
                            self.integrationType.append(sequence.strip())            
                
                        if init.match(line):
                            sequence=init.sub('',line)
                            seq=digit.findall(sequence)
                            for x in range(0,len(seq)):
                                try:
                                    self.initValues[i].append(float(seq[x]))
                                except:
                                    print "\nInitial values for model "+self.name[i]+" have a wrong format in file '"+filename+"'\n"
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
                                    print "\nValue of the prior for model "+self.name[i]+" have a wrong format in file '"+filename+"'\n"
                                    sys.exit()
                            if uni.search(sequence): 
                                prior0[k][0]=2
                                seq0=uni.sub('',sequence)
                                seq=digit.findall(seq0)
                                try:
                                    prior0[k][1]=float(seq[0])
                                    prior0[k][2]=float(seq[1])
                                except:
                                    print "\nValue of the prior for model "+self.name[i]+" have a wrong format in file '"+filename+"'\n"
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

            self.prior.append(prior0)
            self.kernel.append(kernel0)

            f_info.close()
        
        # get fitting rules
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
                    self.fit.append(None)
                    index = index + 1
                else:
                    self.fit.append([])
                    dimensionsOfModel = len(self.initValues[index])
                    self.fit[-1]=line.strip().split(',')
                    for i in range(0, len(self.fit[-1])):
                        self.fit[-1][i] = self.fit[-1][i].strip()
                        if speciesZeroRE.search(self.fit[-1][i]):
                            print "\nIncorrect species indexing!\n"
                            sys.exit()
                        if dimensionsOfModel == 1:
                            self.fit[-1][i]=speciesIntRE.sub('samplePoints', self.fit[-1][i])
                        else:
                            for m in speciesIntRE.finditer(self.fit[-1][i]):
                                mySpecies = str(m.group())
                                mySpecies = operatorRE.sub('', mySpecies)                          
                                mySpeciesIndex = speciesRE.sub('', mySpecies)
                                mySpecies = re.compile(mySpecies)
                                mySpeciesIndex = int(mySpeciesIndex)
                                mySpeciesIndex = mySpeciesIndex-1
                                mySpeciesIndex = str(mySpeciesIndex)
                                self.fit[-1][i]=mySpecies.sub('samplePoints[:,'+mySpeciesIndex+']', self.fit[-1][i])

        f_info.close()
        
    
    #return restart,name,data,timepoints,numOutput,epsilon,initValues,integrationType,modelWeight,prior,kernel,sampleFromPrior,source,fit,BETA,dt, rtol, atol,constKern,modelkernel
        
    
