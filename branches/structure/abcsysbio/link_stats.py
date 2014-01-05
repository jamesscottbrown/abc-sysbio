# sampling link functions

import numpy
from numpy import random as rnd

# define class to hold link information
class link_stats:

    def __init__(self, model, linkp, nfixed, adapt, regf):

        self.linkp = linkp
        self.nfixed = nfixed
        self.nparameters = model.nparameters
        self.adaptive = adapt
        self.regf = regf

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

        # calculate where the links occur in the kernel, non constant, parameter list
        self.klocations = [0 for i in range(self.nlinks)]
        count_non_const = 0
        count = 0
        for n in range(self.nparameters):
            if model.prior[n][0] != 0:
                
                if model.prior[n][0] == 4:
                    self.klocations[count] = count_non_const
                    count += 1

                count_non_const +=1

    # This function takes the previous population link probabilities and
    # tweaks them by adding a uniform component in order to provide some level of robustness
    def regularise_kernel(self, kernel ):

        for i in range(self.nlinks):
            kdist = kernel[2][ self.klocations[i] ]
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

            kernel[2][ self.klocations[i] ] = ret

        #print kernel[2]

    def sample_links_from_prior(self, reti):
        if self.nfixed == -1:
            return link_stats.sample_links_from_prior_free(self, reti)
        else:
            return link_stats.sample_links_from_prior_fixed(self, reti)

    def getLinksPriorPdf(self, params ):
         if self.nfixed == -1:
             return link_stats.getLinksPriorPdf_free(self, params)
         else:
             return link_stats.getLinksPriorPdf_fixed(self, params)

    def perturbLinks(self, params, kernel):
        if self.adaptive == 0:
            if self.nfixed == -1:
                return link_stats.perturbLinks_1(self, params)
            else:
                return link_stats.perturbLinks_prior(self,params)
        else:
            if self.nfixed == -1:
                return link_stats.perturbLinks_1_adaptive(self, params, kernel)
            else:
                print "adaptive and fixed links not implemented"
                exit()

    def getLinksKernelPdf(self, params, params0, kernel ):
        if self.adaptive == 0:
            if self.nfixed == -1:
                return link_stats.getLinksKernelPdf_1(self, params, params0) ## looks at individual links
            else:
                return link_stats.getLinksKernelPdf_2(self, params, params0) ## looks at all the links together
        else:
            if self.nfixed == -1:
                return link_stats.getLinksKernelPdf_1_adaptive(self, params, params0, kernel)
            else:
                print "adaptive and fixed links not implemented"
                exit()
            

    ######################################################
    # FIXED EDGES
    def sample_links_from_prior_fixed(self, reti):
   
        ##print "prior before:", reti

        accept = False
        while(accept==False):
            links = [0 for i in range(self.nlinks)]
            for i in range(self.nlinks):
                if self.prior_type[i] == 3:
                    v = [-1,0,1]
                    str_prior = [1/float(3), 1/float(3), 1/float(3) ]
                    links[i] = v[ numpy.where(rnd.multinomial(n=1,pvals=str_prior,size=1)[0]==1)[0][0]]

                elif self.prior_type[i] == -2:
                   v = [-1,0]
                   str_prior = [0.5, 0.5]
                   links[i] = v[ numpy.where(rnd.multinomial(n=1,pvals=str_prior,size=1)[0]==1)[0][0]]

                else:
                   v = [0,1]
                   str_prior = [0.5, 0.5]
                   links[i] = v[ numpy.where(rnd.multinomial(n=1,pvals=str_prior,size=1)[0]==1)[0][0]]

            total_links = sum( numpy.abs(links) )

            # if nlinks is greater than zero then we can the exact number
            # in nlinks is less than zero we want less than
            # ie nlinks = 3 means nlinks == 3 but nlinks == -3 means nlinks <= 3
            if( self.nfixed > 0 and total_links == self.nfixed ):
                accept = True
            elif( self.nfixed < 0 and total_links <= abs(self.nfixed) ):
                accept = True

        ## print 'sample from prior:', links
        # loop over and fill links
        for i in range(self.nlinks):
            reti[ self.locations[i] ] = links[i]

        return reti

    
    def getLinksPriorPdf_fixed(self, params):

        links = [0 for i in range(self.nlinks)]
        for i in range(self.nlinks):
            links[i] = params[ self.locations[i] ]

        total_links = sum( numpy.abs(links) )

        prior_prob = 1.0
        if self.nfixed > 0:
            if total_links == self.nfixed:
                prior_prob = 1.0 # should really be 1/(total number of possible networks)
            else:
                prior_prob = 0.0

        if self.nfixed < 0:
            if total_links <= abs(self.nfixed):
                prior_prob = 1.0 # should really be 1/(total number of possible networks)
            else:
                prior_prob = 0.0

        return prior_prob
   
    ######################################################
    # FREE EDGES
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

       
    ####################### KERNELS

    # This is an independent kernel, effectively sampling from the prior distribution
    def perturbLinks_prior(self,params):
        
        # change the links with probability linkp
        u = rnd.uniform(0,1)
        if u < self.linkp:
               params = link_stats.sample_links_from_prior_fixed(self, params)
        
        return params

    # This is an independent kernel, effectively sampling links from the previous population
    # There is no dependence on the current particle
    def perturbLinks_1_adaptive(self, params, kernel):

        links = [0 for i in range(self.nlinks)]
        newlinks = [0 for i in range(self.nlinks)]

        ss = [-1, 0, +1]

        for i in range(self.nlinks):
            links[i] = params[ self.locations[i] ]
            kdist = kernel[2][ self.klocations[i] ]

            # change the link with probability linkp
            #u = rnd.uniform(0,1)
            #if u < self.linkp:
            #    #newlinks[i] = rnd.choice( [-1, 0, +1], p=kdist)
            #    newlinks[i] = ss[ numpy.where( numpy.random.multinomial(1, kdist ) == 1 )[0] ]
            #else:
            #    newlinks[i] = links[i]

            newlinks[i] = ss[ numpy.where( numpy.random.multinomial(1, kdist ) == 1 )[0] ]

        #print links, "->", newlinks 

        for i in range(self.nlinks):
            params[ self.locations[i] ] = newlinks[i]

        return params

    # Here every link is perturbed independently
    def perturbLinks_1(self, params):
       
        links = [0 for i in range(self.nlinks)]
        newlinks = [0 for i in range(self.nlinks)]

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

        ##print links, "->", newlinks 

        for i in range(self.nlinks):
            params[ self.locations[i] ] = newlinks[i]

        return params

    # no constraints and adaptive
    def getLinksKernelPdf_1_adaptive(self, params, params0, kernel):
        prob=1
        
        for i in range(self.nlinks):
            kdist = kernel[2][ self.klocations[i] ]
            
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

    # no constraints
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

    # constraints
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



## def perturbLinks_2(kernel, params, priors, link_info):
##     # Here only one link is perturbed independently
    
##     # change one link with probability linkp
##     u = rnd.uniform(0,1)
##     if u < link_info.linkp :
##         # choose a link to perturb uniformly
##         v = [i for i in range(link_info.nlinks) ]
##         p = [1/float(nlinks) for i in range(links_info.nlinks) ]
##         linknum = v[ numpy.where(rnd.multinomial(n=1,pvals=p,size=1)[0]==1)[0][0] ]

##         #print "linknum:", linknum
##         # loop over to modify link num
##         count = 0
##         for n in kernel[0]:
##             if priors[n][0] == 4:
##                 if linknum == count:
##                     # perturb this link

##                     s = 0
##                     if priors[n][1] == 3:
##                         s = set([-1, 0, 1])

##                     elif priors[n][1] == 2:
##                         s = set([0, 1])

##                     else:
##                         s = set([-1, 0])

##                     ss = set( [params[n]] )
##                     ar = numpy.array( list(s-ss) )
##                     rnd.shuffle( ar )
##                     ret = ar[0]

##                     #print linknum, params[n], "->", ret 
##                     params[n] = ret

##                 count += 1
        
##     return params



