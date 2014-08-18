# sampling link functions

import numpy
from numpy import random as rnd
#import statsmodels.api as sm

# define class to hold link information
class link_stats:

    def __init__(self, model, linkp, nfixed, adapt, regf, enum_file):

        self.linkp = linkp
        self.nfixed = nfixed
        self.nparameters = model.nparameters
        self.adaptive = adapt
        self.regf = regf
        self.prior_size = -1
        self.enum = False
        self.enum_file = enum_file

        # prior is either free or enumerated
        if enum_file != "":
            self.enum = True
        else:
            self.enum = False
            if nfixed != -1:
                print "Free prior and fixed links no longer implemented"
                exit()
                

        # store the kernel internally to link_stats
        self.kernel = []

        self.nlinks = 0
        for n in range(self.nparameters):
            if model.prior[n][0] == 4:
                self.nlinks += 1

        
        # get the types of discrete priors and their locations
        self.prior_type = [0 for i in range(self.nlinks)]
        self.p0s = [0 for i in range(self.nlinks)]
        self.locations = [0 for i in range(self.nlinks)]
        count = 0
        for n in range(self.nparameters):		
            if model.prior[n][0] == 4:
                self.locations[count] = n
                self.prior_type[count] = model.prior[n][1]
                self.p0s[count] = model.prior[n][2]
                count += 1

        # get the enumerated prior
        if self.enum == True:
            pfile = self.enum_file
            self.prior_size = 0
            self.prior_data = []
            self.prior_dict = dict()

            of = open(pfile,'r')
            for l in of:
                fields = l.split(" ")
                code = fields[0]
                nedges = int(fields[1])
                edges = [int(i) for i in fields[2:(2+self.nlinks)]]

                # print code, nedges, edges
                
                # if this is a free edges model keep all of the prior
                if self.nfixed == -1:
                    self.prior_size = self.prior_size + 1
                    self.prior_data.append( edges )
                    self.prior_dict[ code ] = 0

                # maximum number of edges model
                if self.nfixed < -1 and nedges <= abs(self.nfixed) :
                    self.prior_size = self.prior_size + 1
                    self.prior_data.append( edges )
                    self.prior_dict[ code ] = 0

                # fixed number model
                if self.nfixed > 0 and nedges == self.nfixed :
                    self.prior_size = self.prior_size + 1
                    self.prior_data.append( edges )
                    self.prior_dict[ code ] = 0

            #print self.prior_data
            #print self.prior_dict
            #exit()

        print "EDGE PRIORS\n"
        if self.enum == True:
            print "\tENUMERATED EDGE PRIOR"
            print "\t\tenumeration file :", self.enum_file
            print "\t\tprior size       :", self.prior_size
            
            if self.nfixed != -1:
                print "\nREDUCED MODEL SPACE"
                print "\tnfixed           :", self.nfixed
        else:
            print "\tFREE EDGE PRIOR"

        print "EDGE KERNELS\n"
        if self.adaptive > 0:
            print "\tADAPTIVE EDGE KERNELS"
            print "\t\tregularisation   :", self.regf
        else:
            print "\tDEPENDENT EDGE KERNELS"
            print "\t\tlinkp            :", self.linkp
        

            
            

        
    ########################################################3
    # THIS IS WHERE WE CONTROL ALL THE OPTIONS
    def sample_links_from_prior(self, reti):
        if self.enum == 0:
            return link_stats.sample_links_from_prior_free(self, reti)
        else:
            return link_stats.sample_links_from_prior_fixed_em(self, reti)


    def getLinksPriorPdf(self, params ):
         if self.enum == 0:
             return link_stats.getLinksPriorPdf_free(self, params)
         else:
             return link_stats.getLinksPriorPdf_fixed_em(self, params)

    def perturbLinks(self, params):
        if self.adaptive == 0:
            return link_stats.perturbLinks_1(self, params)
        else:
            return link_stats.perturbLinks_1_adaptive(self, params)


    def getLinksKernelPdf(self, params, params0):
        if self.adaptive == 0:
            return link_stats.getLinksKernelPdf_1(self, params, params0) ## looks at individual links
        else:
            return link_stats.getLinksKernelPdf_1_adaptive(self, params, params0)
        

    # get kernels for the links
    def getKernels(self, population, weights):
        
        if self.adaptive == 1:
            tmp=list()
            for i in range(self.nlinks):
                # extract the link info from the data
                xx = numpy.array( population[:, self.locations[i] ] )

                #pm1 = len( numpy.where(xx == -1)[0] )/float(len(xx))
                #p0  = len( numpy.where(xx ==  0)[0] )/float(len(xx))
                #pp1 = len( numpy.where(xx == +1)[0] )/float(len(xx))
                pm1 = numpy.sum( weights[ numpy.where(xx == -1)[0] ] )
                p0  = numpy.sum( weights[ numpy.where(xx ==  0)[0] ] )
                pp1 = numpy.sum( weights[ numpy.where(xx == +1)[0] ] )
                tmp.append( [pm1, p0, pp1] )

            self.kernel = tmp 

            # regularise
            link_stats.regularise_kernel(self)
    
    # This function takes the previous population link probabilities and
    # tweaks them by adding a uniform component in order to provide some level of robustness
    def regularise_kernel(self):

        for i in range(self.nlinks):
            kdist = self.kernel[i]
            ret = kdist[:]
        
            if self.prior_type[i] == 3:
                # v = [-1,0,+1]
                ret[0] += self.regf/3.0
                ret[1] += self.regf/3.0
                ret[2] += self.regf/3.0

            if self.prior_type[i] == -2:
                # v = [-1,0]
                ret[0] += self.regf/2.0
                ret[1] += self.regf/2.0

            if self.prior_type[i] == 2:
                # v = [0,1]
                ret[1] += self.regf/2.0
                ret[2] += self.regf/2.0

            sr = sum(ret)
            ret[0] = ret[0]/sr
            ret[1] = ret[1]/sr
            ret[2] = ret[2]/sr
        
            print "kernel reg:", self.regf, kdist, ret

            self.kernel[i] = ret

        #print kernel[2]

    ######################################################
    # SAMPLING FROM THE PRIOR
    def sample_links_from_prior_fixed_em(self, reti):
   
        # sample a number from 1 to prior_size
        u = numpy.random.randint(0,self.prior_size,1)

        links = [ self.prior_data[u][i] for i in range(self.nlinks)]
        
        # print 'sample from prior:', links
        # loop over and fill links from the precomputed prior

        for i in range(self.nlinks):
            reti[ self.locations[i] ] = links[i]

        return reti

    def sample_links_from_prior_free(self, reti):
  
        links = [0 for i in range(self.nlinks)]
        for i in range(self.nlinks):
            p0 = self.p0s[i]
            #print p0
            #print self.prior_type[i]

            if self.prior_type[i] == 3:
                v = [-1,0,1]
                str_prior = [(1-p0)/2, p0, (1-p0)/2]
                links[i] = v[ numpy.where(rnd.multinomial(n=1,pvals=str_prior,size=1)[0]==1)[0][0]]

            if self.prior_type[i] == -2:
                v = [-1,0]
                str_prior = [(1-p0), p0 ]
                links[i] = v[ numpy.where(rnd.multinomial(n=1,pvals=str_prior,size=1)[0]==1)[0][0]]

            if self.prior_type[i] == 2:
                v = [0,1]
                str_prior = [p0, (1-p0) ]
                #print str_prior, v[ numpy.where(rnd.multinomial(n=1,pvals=str_prior,size=1)[0]==1)[0][0]]
                links[i] = v[ numpy.where(rnd.multinomial(n=1,pvals=str_prior,size=1)[0]==1)[0][0]]

        # print 'sample from prior:', links
        # loop over and fill links
        for i in range(self.nlinks):
            reti[ self.locations[i] ] = links[i]

        return reti

    ######################################################
    # PRIOR PDF
    def getLinksPriorPdf_fixed_em(self, params):

        # In this special case we have constructed prior
        # such that all hard work has been done
        prior_prob = 1/float(self.prior_size)

        return prior_prob

    def getLinksPriorPdf_free(self, params):

        links = [0 for i in range(self.nlinks)]

        prior_prob = 1.0
        for i in range(self.nlinks):
            links[i] = params[ self.locations[i] ]

            p0 = self.p0s[i]

            x = 1.0
            if self.prior_type[i] == 3:
                if links[i] == 0 : x = p0
                else :  x = (1-p0)/2
            else:
                if links[i] == 0: x = p0
                else:  x = 1-p0

            prior_prob = prior_prob*x

        ## print "prior prob:", prior_prob
        return prior_prob


    ######################################################
    # PERTURBATION KERNELS

    ######################################################
    # Here every link is perturbed independently
    def perturbLinks_1(self, params):
       
        links = [0 for i in range(self.nlinks)]
        newlinks = [0 for i in range(self.nlinks)]

        # sample a network valid under the prior
        prior_prob = 0
        while prior_prob == 0:
            for i in range(self.nlinks):
                links[i] = params[ self.locations[i] ]

                s = 0
                if self.prior_type[i] == 3:
                    s = set([-1, 0, 1])
                elif  self.prior_type[i] == 2:
                    s = set([0, 1])
                else:
                    s = set([-1, 0])

                # change the link with probability linkp
                u = rnd.uniform(0,1)
                if u < self.linkp:
                    ss = set( [ links[i] ] )
                    ar = numpy.array( list(s-ss) )
                    rnd.shuffle( ar )
                    newlinks[i] = ar[0]
                else:
                    newlinks[i] = links[i]

            if self.enum == True:
                # code newlinks and check the prior dictionary
                code = "".join([str(k+1) for k in newlinks])
                # print newlinks, code
                if code in self.prior_dict:
                    prior_prob = 1
                else:
                    prior_prob = 0
            else:
                # everything exists under the prior
                prior_prob = 1

        ##print links, "->", newlinks 

        for i in range(self.nlinks):
            params[ self.locations[i] ] = newlinks[i]

        return params

    def getLinksKernelPdf_1(self, params, params0):
        prob=1
        
        for i in range(self.nlinks):
            # p( t | t -1 ) =  1-linkp (t == t-1), linkp (t != t-1 )
            kern = 0;
            if params[ self.locations[i] ] == params0[ self.locations[i] ]:
                kern = 1-self.linkp
            else:
                kern = self.linkp
           
            prob=prob*kern
        ##print "pdf kernel:", prob
        return prob

    ######################################################
    
    # This is an independent kernel, effectively sampling links from the previous population
    # There is no dependence on the current particle
    def perturbLinks_1_adaptive(self, params):
        #print "perturbing links"

        links = [0 for i in range(self.nlinks)]
        newlinks = [0 for i in range(self.nlinks)]

        ss = [-1, 0, +1]

        # sample a network valid under the prior
        prior_prob = 0
        while prior_prob == 0:

            for i in range(self.nlinks):
                links[i] = params[ self.locations[i] ]
                kdist = self.kernel[i]

                newlinks[i] = ss[ numpy.where( numpy.random.multinomial(1, kdist ) == 1 )[0] ]

            if self.enum == True:
                # code newlinks and check the prior dictionary
                code = "".join([str(k+1) for k in newlinks])
                # print newlinks, code
                if code in self.prior_dict:
                    prior_prob = 1
                else:
                    prior_prob = 0
            else:
                # everything exists under the prior
                prior_prob = 1
                    
        #print links, "->", newlinks 

        for i in range(self.nlinks):
            params[ self.locations[i] ] = newlinks[i]

        return params

    # adaptive (independent of params0)
    def getLinksKernelPdf_1_adaptive(self, params, params0):
        prob=1
        
        for i in range(self.nlinks):
            kdist = self.kernel[i]
            
            # p( t | t -1 ) =  1-linkp (t == t-1), linkp (t != t-1 )
            #kern = 0;
            #if params[ self.locations[i] ] == params0[ self.locations[i] ]:
            #    kern = 1-self.linkp
            #else:
            #    kern = self.linkp

            # Here the probability just depends on the kernel
            ind = params[ self.locations[i] ] + 1 # 0, 1, 2
            kern = kdist[ ind ]
           
            prob=prob*kern
        ##print "pdf kernel:", prob
        return prob

    ######################################################

    # perturb by sampling from the enumerated prior
    def perturbLinks_prior_em(self,params):
        
        # change the links with probability linkp
        u = rnd.uniform(0,1)
        if u < self.linkp:
            params = link_stats.sample_links_from_prior_fixed_em(self, params)
        
        return params

    def getLinksKernelPdf_2(self, params, params0):
    
        nochange = True
        prob = (1-self.linkp)

        for i in range(self.nlinks):
            if not (params[ self.locations[i] ] == params0[ self.locations[i] ]):
                nochange = False

        if nochange == False:
            prob = self.linkp
        
        ## print "pdf kernel:", prob
        return prob



