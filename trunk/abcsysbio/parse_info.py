# Algorithm information

import re, sys, numpy

from xml.dom import minidom

# implemented priors
re_prior_const=re.compile('constant')
re_prior_uni=re.compile('uniform')
re_prior_normal=re.compile('normal')
re_prior_logn=re.compile('lognormal') 

# implemented kernels
re_kernel_uniform=re.compile('uniform')
re_kernel_normal=re.compile('normal')

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

def parse_required_vector_value( node, tagname, message, cast ):
    try:
        data = node.getElementsByTagName(tagname)[0].firstChild.data
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
    
    def __init__(self, filename):
        xmldoc = minidom.parse(filename)

        self.modelnumber = 0
        self.restart = False
        self.particles = 0
        self.beta = 0
        self.dt = 0
        self.epsilon = []
        self.times = []
        self.ntimes = 0
        self.data = []
        
        
        self.nmodels = 0
        self.nparameters = []
        self.nspecies = []
        self.name = []
        self.source = []
        self.init = []
        self.type = []
        self.prior = []
        self.fit = []

        self.modelkernel = 0.7
        self.kernel = 0
        self.modelprior = []
        self.rtol = 1e-5
        self.atol = 1e-5

        ##################################################
        ## Required arguments

        ### get number of models
        self.modelnumber = parse_required_single_value( xmldoc, "modelnumber", "Please provide an integer value for <modelnumber>", int )

        ### get number of particles
        self.particles = parse_required_single_value( xmldoc, "particles", "Please provide an integer value for <particles>", int )

        ### get beta value
        self.beta = parse_required_single_value( xmldoc, "beta", "Please provide an integer value for <beta>", int )
        
        ### get dt
        self.dt = parse_required_single_value( xmldoc, "dt", "Please provide an float value for <dt>", float )

        ### get epsilon
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

        ##self.epsilon = parse_required_vector_value( xmldoc, "epsilon", "Please provide a whitespace separated list of values for <epsilon>" , float )

        ### get data attributes
        dataref = xmldoc.getElementsByTagName('data')[0]
        # times
        self.times = parse_required_vector_value( dataref, "times", "Please provide a whitespace separated list of values for <data><times>" , float )
        self.ntimes = len(self.times)

        # variables
        # first do a scan to get the number of timeseries
        nvar = 0
        varref = dataref.getElementsByTagName('variables')[0]
        for v in varref.childNodes:
            if v.nodeType == v.ELEMENT_NODE:
                nvar += 1

        # create matrix
        self.data = numpy.zeros([self.ntimes,nvar])
        nvar = 0
        for v in varref.childNodes:
            if v.nodeType == v.ELEMENT_NODE:
                tmp = str( v.firstChild.data ).split()

                for i in range(self.ntimes):
                    self.data[i,nvar] = float( tmp[i] )

                nvar += 1        
                
                
        ### get model attributes
        modelref = xmldoc.getElementsByTagName('models')[0]
        for m in modelref.childNodes:
            if m.nodeType == v.ELEMENT_NODE:
                self.nmodels += 1
                self.prior.append([])

                self.name.append( str(m.getElementsByTagName('name')[0].firstChild.data).strip() )
                self.source.append( str(m.getElementsByTagName('source')[0].firstChild.data).strip() )
                self.type.append( str(m.getElementsByTagName('type')[0].firstChild.data).strip() )
                
                self.fit.append( parse_fitting_information( m )  )

                initref = m.getElementsByTagName('initialvalues')[0]
                tmp = str( initref.firstChild.data ).split()
                self.init.append( [ float(i) for i in tmp ] )
                self.nspecies.append( len( self.init[self.nmodels-1] ) )

                nparameter = 0
                paramref = m.getElementsByTagName('parameters')[0]
                for p in paramref.childNodes:
                    if p.nodeType == p.ELEMENT_NODE:
                        nparameter += 1
                        prior_tmp = [0,0,0]
                        tmp = str( p.firstChild.data ).split()
                        
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
                        else:
                            print "\nIn model", self.name[self.nmodels-1], " supplied parameter prior ", tmp[0], " unsupported"
                            sys.exit()
                        
                        self.prior[self.nmodels-1].append( prior_tmp )

                if nparameter == 0:
                    print "\nNo parameters specified in model ", self.name[self.nmodels-1]
                    sys.exit()
                self.nparameters.append( nparameter )
                
        if self.nmodels == 0:
            print "\nNo models specified"
            sys.exit()

        ##################################################
        ## Optional arguments

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
        try:
            data = str(xmldoc.getElementsByTagName('kernel')[0].firstChild.data).strip()
            if re_kernel_uniform.match( data ):
                self.kernel = 0
            elif re_kernel_normal.match( data ):
                self.kernel = 1
            else:
                print "\n#################\n<kernel> must be one of uniform, normal  so I am going to ignore your argument"
        except:
            null = 0


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
        print "restart:", self.restart
        print "particles:", self.particles
        print "beta:", self.beta
        print "dt:", self.dt
        print "epsilon:" 
        for i in range(self.epsilon.shape[0]):
            print "\t", 
            for j in range(self.epsilon.shape[1]):
                print "", self.epsilon[i,j],
            print ""
        print "kernel:", self.kernel
        print "model kernel:", self.modelkernel
        print "model prior:", self.modelprior
        
        print "DATA:"
        print "\ttimes:", self.times
        print "\tvars:" 
        for i in range(len(self.data[0,:])):
            print "\t", 
            for j in range(self.ntimes):
                print "", self.data[j,i],
            print ""
        
        print "MODELS:", self.nmodels
        for i in range(self.nmodels):
            print "\t", "npar:", self.nparameters[i]
            print "\t", "nspecies:", self.nspecies[i]
            print "\t", "name:", self.name[i]
            print "\t", "source:", self.source[i]
            print "\t", "type:", self.type[i]
            print "\t", "fit:", self.fit[i]
            print "\t", "init:", self.init[i]
            print "\t", "prior:", self.prior[i]
            print "\n"

###x = algorithm_info('input.xml')
###x.print_info()
