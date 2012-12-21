# Algorithm information

import re, sys, numpy, copy

from xml.dom import minidom

# implemented priors
re_prior_const=re.compile('constant')
re_prior_uni=re.compile('uniform')
re_prior_normal=re.compile('normal')
re_prior_logn=re.compile('lognormal') 
re_prior_trn=re.compile('trunc_normal')


# implemented kernels
re_kernel_uniform=re.compile('uniform')
re_kernel_normal=re.compile('normal')
re_kernel_mvnormal=re.compile('multiVariateNormal')
re_kernel_mvnormalKN=re.compile('multiVariateNormalKNeigh')
re_kernel_mvnormalOCM=re.compile('multiVariateNormalOCM')
re_kernel_nonadapt_uniform=re.compile('nonadapt_uniform\s\d+\.\?\d*\s\d+\.?\d*')
re_kernel_nonadapt_normal=re.compile('nonadapt_normal\s\d+\.?\d*')

# True/False
re_true=re.compile('True')
re_none=re.compile('None')

def parse_required_single_value( node, tagname, message, cast ):
    try:
        data = node.getElementsByTagName(tagname)[0].firstChild.data
    except:
        print message
        sys.exit()

    ret = 0
    try:
        ret = cast( data )
    except:
        print message
        sys.exit()

    return(ret)

def parse_required_vector_value( node, tagname, message, cast, j ): 
    try:
        data = node.getElementsByTagName(tagname)[j].firstChild.data 
    except:
        print message
        sys.exit()

    tmp = str( data ).split()
    ret = []
    try:
        ret = [ cast(i) for i in tmp ]
    except:
        print message
        sys.exit()

    if len(ret) == 0:
        print message
        sys.exit()

    return(ret)

def process_prior( tmp ):
    prior_tmp = [0,0,0,0]

    if re_prior_const.match( tmp[0] ):
        prior_tmp[0] = 0
        try:
            prior_tmp[1] = float( tmp[1] )
        except:
            print "\nValue of the prior for model ", self.name[self.nmodels-1], "has the wrong format:", tmp[1]
            sys.exit()
                                
    elif re_prior_normal.match( tmp[0] ):
        prior_tmp[0] = 1
        try:
            prior_tmp[1] = float( tmp[1] )
            prior_tmp[2] = float( tmp[2] )
        except:
            print "\nValue of the prior for model ", self.name[self.nmodels-1], "has the wrong format:", tmp[1]
            sys.exit()

    elif re_prior_uni.match( tmp[0] ):
        prior_tmp[0] = 2
        try:
            prior_tmp[1] = float( tmp[1] )
            prior_tmp[2] = float( tmp[2] )
        except:
            print "\nValue of the prior for model ", self.name[self.nmodels-1], "has the wrong format:", tmp[1]
            sys.exit()
                                
    elif re_prior_logn.match( tmp[0] ):
        prior_tmp[0] = 3
        try:
            prior_tmp[1] = float( tmp[1] )
            prior_tmp[2] = float( tmp[2] )
        except:
            print "\nValue of the prior for model ", self.name[self.nmodels-1], "has the wrong format:", tmp[1]
            sys.exit()
            
    elif re_prior_trn.match( tmp[0] ):      # TRUNCATED NORMAL PRIOR BP 29/2
        prior_tmp[0] = 4
        try:
            prior_tmp[1] = float( tmp[1] )
            prior_tmp[2] = float( tmp[2] )

            if tmp[3] == "i" and tmp[4] != "i":
                prior_tmp[3] = ([ neg_inifity, float(tmp[4]) ])
            elif tmp[3] != "i" and tmp[4] == "i":
                prior_tmp[3] = ([ float(tmp[3]), pos_infinity ])
            else:
                prior_tmp[3] = ( [float( tmp[3] ), float( tmp[4] )] )     # BP 29/2 - Upper and lower bounds 
        except:
            print "\nValue of the prior for model ", self.name[self.nmodels-1], "has the wrong format:", tmp[1]
            sys.exit()

    else:
        print "\nSupplied parameter prior ", tmp[0], " unsupported"
        sys.exit()

    return prior_tmp

def parse_fitting_information( node ):
    fitref = node.getElementsByTagName('fit')[0]
    tmp = str( fitref.firstChild.data ).split()
    ret = []

    if len(tmp) == 1 and re_none.match( tmp[0] ):
        return None
    else:
        for i in tmp:
            # replace species with samplePoints            
            ttmp = re.sub('species','samplePoints', i )

            # find all instances of samplePoints[ ] and extract the list of species numbers
            sp_strs = re.findall("samplePoints([0-9]+)",ttmp)
            sp_nums = [int(j) for j in sp_strs]
            sp_nums.sort()
            sp_nums.reverse()

            # loop over the species numbers and replace
            for n in sp_nums:
                ttmp = re.sub('ts'+str(n),'ts[:,'+str(n-1)+']',ttmp)

            ret.append( ttmp )

        return( ret )

class algorithm_info:
    """
    A class to parse the user-provided input file and return all information required to run the abc-SMC algorithm.
    
    """ 
    
    def __init__(self, filename, mode):
        xmldoc = minidom.parse(filename)
        self.mode = mode
        ### mode is 0  inference, 1 simulate, 2 design

        self.modelnumber = 0
        self.restart = False
        self.particles = 0
        self.beta = 0
        self.dt = 0
        self.epsilon = []
        self.final_epsilon = []
        self.alpha = 0.9
        self.times = [[]]   
        self.ntimes = [[]]
        self.data = [[]]    
        
        
        self.nmodels = 0
        self.nsubmodels = 0 
        self.nparameters = []
        self.nspecies = [[]] 
        self.name = []  
        self.submodelname = [[]]  
        self.source = [[]]  
        self.type = []
        self.prior = []
        self.localprior = [[]] 
        self.x0prior = [[]]  
        self.fit = [[]]  
        self.logp = []
        self.dist = [[]] 

        self.modelkernel = 0.7
        self.kernel = [1,0]
        self.modelprior = []
        self.rtol = 1e-5
        self.atol = 1e-5

        ##################################################
        ## Required arguments

        ### get number of models
        self.modelnumber = parse_required_single_value( xmldoc, "modelnumber", "Please provide an integer value for <modelnumber>", int )
        
        ### get number of submodels 
        self.submodelnumber = parse_required_single_value( xmldoc, "submodelnumber", "Please provide an integer value for <submodelnumber>", int )

        ### get number of particles
        self.particles = parse_required_single_value( xmldoc, "particles", "Please provide an integer value for <particles>", int )

        ### get beta value
        self.beta = parse_required_single_value( xmldoc, "beta", "Please provide an integer value for <beta>", int )
        
        ### get dt
        self.dt = parse_required_single_value( xmldoc, "dt", "Please provide a float value for <dt>", float )

        ### get epsilon
        if self.mode != 1:
            # automated epsilon takes priority
            if( len( xmldoc.getElementsByTagName('autoepsilon') ) > 0 ):
                # found automated epsilon
                epsref = xmldoc.getElementsByTagName('autoepsilon')[0]
                self.final_epsilon = parse_required_vector_value( epsref, "finalepsilon", "Please provide a whitespace separated list of values for <autoepsilon><finalepsilon>" , float, 0 ) 
                try:
                    self.alpha = parse_required_single_value( epsref, "alpha", "Please provide a float value for <autoepsilon><alpha>" , float)
                except:
                    null=0
            else:
                # first do a scan to get the number of epsilon series and the number of epsilon in the series
                neps1 = 0
                neps2 = 0
                epsref = xmldoc.getElementsByTagName('epsilon')[0]
                for e in epsref.childNodes:
                    if e.nodeType == e.ELEMENT_NODE:
                        neps1 += 1
                        neps2 = len( str( e.firstChild.data ).split() )

                # create matrix
                self.epsilon = numpy.zeros([neps1,neps2])
                i1 = 0
                for e in epsref.childNodes:
                    if e.nodeType == e.ELEMENT_NODE:
                        tmp = str( e.firstChild.data ).split()

                        for i in range(neps2):
                            self.epsilon[i1,i] = float( tmp[i] )

                        i1 += 1  

        ### get data attributes
        dataref = xmldoc.getElementsByTagName('data')[0]
        timesets = xmldoc.getElementsByTagName('times') 
  
        # times
        nvar = [[]] 
        first = True 
        for i in range(0,len(timesets)): 
            if first != True: 
                self.times.append([]) 
                self.ntimes.append([]) 
                self.data.append([]) 
                nvar.append([]) 
            first = False 
            self.times[i]=parse_required_vector_value(dataref,"times","Please provide a whitespace separated list of values for <data><times>",float , i) 

            self.ntimes[i] = len(self.times[i])

            # variables
            if self.mode == 0: 
                # first do a scan to get the number of timeseries
                nvar[i] = 0 
                varref = dataref.getElementsByTagName('variables')[i]  
                for v in varref.childNodes:
                    if v.nodeType == v.ELEMENT_NODE:
                        nvar[i] += 1 

                # create matrix and mask
                data_unmasked = numpy.zeros([self.ntimes[i],nvar[i]]) 
                data_mask = numpy.zeros([self.ntimes[i],nvar[i]], dtype=numpy.int32) 
                nvar[i] = 0 
                for v in varref.childNodes:
                    if v.nodeType == v.ELEMENT_NODE:
                        tmp = str( v.firstChild.data ).split()
                        
                        for j in range(self.ntimes[i]): 
                        # Search for NA
                            if re.match("\s*NA\s*", tmp[j]) != None:
                                data_mask[j,nvar[i]] = 1 
                                tmp[j] = 0

                            data_unmasked[j,nvar[i]] = float( tmp[j] ) 

                        nvar[i] += 1         

                # create masked data
                self.data[i] = numpy.ma.array(data_unmasked, mask = data_mask) 

        ### get model attributes
        modelref = xmldoc.getElementsByTagName('models')[0]  
        for m in modelref.childNodes:
            if m.nodeType == m.ELEMENT_NODE:
                self.nmodels += 1
                if self.nmodels > 1:     
                    self.submodelname.append([]) 
                    self.localprior.append([])
                    self.source.append( []) 
                    self.x0prior.append([]) 
                    self.fit.append([]) 
                    self.nspecies.append( [] ) 
                    self.dist.append([]) 
                self.name.append(str(m.getElementsByTagName('name')[0].firstChild.data).strip() )
                self.prior.append([])
                self.type.append( str(m.getElementsByTagName('type')[0].firstChild.data).strip() )
                try:
                    tmp = str( m.getElementsByTagName('logp')[0].firstChild.data ).strip()
                    if re_true.match( tmp ):
                        self.logp.append( True )
                    else:
                        self.logp.append( False )
                except:
                    self.logp.append( False )
                
                
                nparameter = 0
                paramref = xmldoc.getElementsByTagName('parameters')[self.nmodels-1]
                for p in paramref.childNodes: 
                    if p.nodeType == p.ELEMENT_NODE:
                        nparameter += 1
                        prior_tmp = [0,0,0]
                        tmp = str( p.firstChild.data ).split()
                        self.prior[self.nmodels-1].append( process_prior( tmp ) )

                if nparameter == 0:
                    print "\nNo parameters specified in model ", self.name[self.nmodels-1]
                    sys.exit()
        
                submodelref = xmldoc.getElementsByTagName('submodels')[self.nmodels-1]
                self.nsubmodels = 0
                
                for n in submodelref.childNodes:
                    if n.nodeType == n.ELEMENT_NODE:
                        self.nsubmodels += 1
                        self.localprior[self.nmodels-1].append([])
                        self.x0prior[self.nmodels-1].append([])
                        self.submodelname[self.nmodels-1].append(str(n.getElementsByTagName('name')[0].firstChild.data).strip() )
                        self.source[self.nmodels-1].append( str(n.getElementsByTagName('source')[0].firstChild.data).strip() )
                        self.fit[self.nmodels-1].append( parse_fitting_information( n )  )
                        
                        #initref = m.getElementsByTagName('initialvalues')[0]
                        #tmp = str( initref.firstChild.data ).split()
                        #self.init.append( [ float(i) for i in tmp ] )
                        #self.nspecies.append( len( self.init[self.nmodels-1] ) )

                        if( len(n.getElementsByTagName('localparameters')) > 0 ):
                            paramref = n.getElementsByTagName('localparameters')[0]
                            for p in paramref.childNodes: 
                                if p.nodeType == p.ELEMENT_NODE:
                                    nparameter += 1
                                    prior_tmp = [0,0,0]
                                    tmp = str(p.firstChild.data).split()
                                    self.localprior[self.nmodels-1][self.nsubmodels-1].append( process_prior( tmp ) )
                        
                        if( len(n.getElementsByTagName('distancefunction')) > 0 ):
                            dist = n.getElementsByTagName('distancefunction')[0]
                            
                            self.dist[self.nmodels-1].append(str(n.getElementsByTagName('distancefunction')[0].firstChild.data).strip() )
                        else:
                            self.dist[self.nmodels-1].append('default')
                        
                        
                        ninit = 0
                        initref = n.getElementsByTagName('initial')[0]
                        for inn in initref.childNodes:
                            if inn.nodeType == inn.ELEMENT_NODE:
                                ninit += 1
                                prior_tmp = [0,0,0]
                                tmp = str( inn.firstChild.data ).split()
                                self.x0prior[self.nmodels-1][self.nsubmodels-1].append( process_prior( tmp ) )
                    
                        if ninit == 0:
                            print "\nNo initial conditions specified in model ", self.name[self.nmodels-1]
                            sys.exit()
                        self.nspecies[self.nmodels-1].append( ninit ) 
                
                self.nparameters.append( nparameter )
                

        # create full-length prior and create pindex to identify which parameters belong to each submodel
        self.gprior = [[] for i in range(self.nmodels)]
        self.pindex = [[[[],[]] for i in range(0,len(self.submodelname[j]))] for j in range(0,self.nmodels)]
        self.gparameters = []
        for j in range(0,self.nmodels):
            self.gparameters.append(len(self.prior[j])) 
            for i in range(0,len(self.submodelname[j])):
                self.pindex[j][i][0] = range(self.gparameters[j])
                for x in self.localprior[j][i]:
                    self.prior[j].append( x[:] )
                    self.pindex[j][i][0].append(len(self.prior[j])-1)
            self.gprior[j] = copy.deepcopy(self.prior[j])
            for i in range(0,len(self.submodelname[j])):
                for x in self.x0prior[j][i]:
                    self.prior[j].append( x[:] )
                    self.pindex[j][i][1].append(len(self.prior[j])-1)
        
        if self.nmodels == 0:
            print "\nNo models specified"
            sys.exit()
            

        ## Optional arguments

        ### get atol
        try:
            data = xmldoc.getElementsByTagName('atol')[0].firstChild.data
            self.atol = float(data)
        except:
            null = 0
        
        ### get rtol
        try:
            data = xmldoc.getElementsByTagName('rtol')[0].firstChild.data
            self.rtol = float(data)
        except:
            null = 0    

        ### get restart
        try:
            tmp = str( xmldoc.getElementsByTagName('restart')[0].firstChild.data ).strip()
            if re_true.match( tmp ):
                self.restart = True
        except:
            null = 0

        ### get model kernel
        try:
            data = xmldoc.getElementsByTagName('modelkernel')[0].firstChild.data
            try:
                self.modelkernel = float(data)
            except:
                print "\n#################\n<modelkernel> must be a float so I am going to ignore your argument"
                
            if self.modelkernel > 1.0:
                print "\n#################\n<modelkernel> must be <= 1.0  so I am going to ignore your argument"
                self.modelkernel = 0.7
        except:
            null = 0

        ### get kernel
        #try:
        data = str(xmldoc.getElementsByTagName('kernel')[0].firstChild.data).strip()
        if re_kernel_uniform.match( data ):
            self.kernel[0] = 1
        elif re_kernel_normal.match( data ):
            self.kernel[0] = 2
        elif re_kernel_mvnormal.match( data ):
            self.kernel[0] = 3
        elif re_kernel_mvnormalKN.match( data ):
            self.kernel[0] = 4
        elif re_kernel_mvnormalOCM.match( data ):
            self.kernel[0] = 5
        elif re_kernel_nonadapt_uniform.match( data ):      # BP 7/3
            self.kernel[0] = 7
            tmp = str(data).split()
            self.kernel[1] = [float(tmp[1]), float(tmp[2])]
        elif re_kernel_nonadapt_normal.match( data ):       # BP 7/3
            self.kernel[0]= 8
            tmp = str(data).split()
            self.kernel[1] = float(tmp[1])
        else:
            print "\n#################\n<kernel> must be one of uniform, normal, multivariateNormal, multivariateNormalKNeigh and multivariateNormalOCM  so I am going to ignore your argument"
        #except:
         #   null = 0


        ### get model priors
        self.modelprior = [1/float(self.nmodels) for i in range(self.nmodels) ]
        try:
            data = xmldoc.getElementsByTagName("modelprior")[0].firstChild.data
            tmp = str( data ).split()
            #print tmp
            ret = []
            try:
                 ret = [ float(i) for i in tmp ]
            except:
                print "\n#################\n<modelprior> must be a vector of floats so I am going to ignore your argument"

            if sum(ret) != 1.0 or len(ret) != self.nmodels:
                print "\n#################\n<modelprior> must sum to one and be the same length as the number of models so I am going to ignore your argument"
            else:
                self.modelprior = ret[:]
        except:
            null = 0

        
    def print_info(self):
        print "\nALGORITHM INFO"
        print "modelnumber:", self.modelnumber
        print "submodelnumber:", self.submodelnumber  
        print "restart:", self.restart
        print "particles:", self.particles
        print "beta:", self.beta
        print "dt:", self.dt
        if self.mode != 1:
            if len(self.final_epsilon) == 0:
                print "manual epsilon:" 
                for i in range(self.epsilon.shape[0]):
                    print "\t", 
                    for j in range(self.epsilon.shape[1]):
                        print "", self.epsilon[i,j],
                    print ""
            else:
                print "auto epsilon:" 
                print "\t", self.final_epsilon
                print "\talpha:", self.alpha

            print "kernel:", self.kernel
            print "model kernel:", self.modelkernel
        print "model prior:", self.modelprior

        print "\nDATA:"    
        for i in range(0,len(self.times)):  
            print "\ndataset", i+1, ":"
            print "\ttimes:", self.times[i]    
            if self.mode == 0:
                print "\tvars:"
                for k in range(len(self.data[i][0,:])):
                    print "\t"    
                    for j in range(self.ntimes[i]):    
                        print "", self.data[i][j,k],    
                    print ""    
        ####
        ### models
        print "\nMODELS:", self.nmodels
        for i in range(self.nmodels):
            print "\t", "name:", self.name[i]
            print "\t", "npar:", self.nparameters[i]
            print "\t", "prior:", self.prior[i][0:self.gparameters[i]]
            print "\t", "type:", self.type[i]
            print "\t", "logp:", self.logp[i]
            print "\t", "SUBMODELS:", self.nsubmodels
            ### submodels
            for j in range(self.nsubmodels):
                print "\t\t", "name:", self.submodelname[i][j]
                print "\t\t", "source:", self.source[i][j]
                print "\t\t", "local prior:", self.localprior[i][j]
                print "\t\t", "nspecies:", self.nspecies[i][j]
                print "\t\t", "fit:", self.fit[i][j]
                print "\t\t", "init:", self.x0prior[i][j]
                print "\t\tdistancefn:", self.dist[i][j],"\n"
                
