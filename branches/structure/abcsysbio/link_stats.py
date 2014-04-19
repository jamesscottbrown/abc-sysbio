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

        if enum_file != "":
            self.enum = True
        
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
            prior_raw_data = numpy.genfromtxt( pfile, dtype=numpy.int32, delimiter=" ")

            #print self.prior_data
            #for i in range(self.prior_size):
            #    print self.prior_data[i,:]

            # extract the models within the prior that we need
            new_prior = []
            prior_count = 0
            for i in range(numpy.shape(prior_raw_data)[0]):

                # if this is a free edges model keep all of the prior
                if self.nfixed == -1:
                    prior_count = prior_count + 1
                    new_prior.append( prior_raw_data[i,:])

                # maximum number of edges model
                if self.nfixed < -1 and prior_raw_data[i,1] <= abs(self.nfixed) :
                    prior_count = prior_count + 1
                    new_prior.append( prior_raw_data[i,:])

                # fixed number model
                if self.nfixed > 0 and prior_raw_data[i,1] == self.nfixed :
                    prior_count = prior_count + 1
                    new_prior.append( prior_raw_data[i,:])

            self.prior_size = prior_count
            self.prior_data = new_prior
            
        if self.adaptive > 0:
            print "FREE EDGE PRIORS AND ADAPTIVE EDGE KERNELS"
            print "\tregularisation   :", self.regf

        if self.nfixed != -1:
            print "CONSTRAINED EDGE PRIORS AND INDEPENDENT PRIOR EDGE KERNELS"
            print "\tnfixed           :", self.nfixed
            print "\tlinkp            :", self.linkp

            # ENUMERATED PRIOR CURRENTLY ONLY WORKS WITH CONSTRAINED EDGES
            if self.enum == True:
                print "\tENUMERATED EDGE PRIOR"
                print "\t\tenumeration file :", self.enum_file
                print "\t\tprior size       :", self.prior_size
                print "\t\tprior size       :", len(self.prior_data)

        
    ########################################################3
    # THIS IS WHERE WE CONTROL ALL THE OPTIONS
    def sample_links_from_prior(self, reti):
        if self.nfixed == -1:
            return link_stats.sample_links_from_prior_free(self, reti)
        else:
            if self.enum == 0: return link_stats.sample_links_from_prior_fixed(self, reti)
            else : return link_stats.sample_links_from_prior_fixed_em(self, reti)

    def getLinksPriorPdf(self, params ):
         if self.nfixed == -1:
             return link_stats.getLinksPriorPdf_free(self, params)
         else:
             if self.enum == 0: return link_stats.getLinksPriorPdf_fixed(self, params)
             else : return link_stats.getLinksPriorPdf_fixed_em(self, params)

    def perturbLinks(self, params):
        if self.adaptive == 0:
            if self.nfixed == -1:
                return link_stats.perturbLinks_1(self, params)
            else:
                if self.enum == 0: return link_stats.perturbLinks_prior(self, params)
                else: return link_stats.perturbLinks_prior_em(self, params)
        else:
            if self.nfixed == -1:
                if self.adaptive == 1: return link_stats.perturbLinks_1_adaptive(self, params)
                #if self.adaptive == 2: link_stats.perturbLinks_1_logistic(self, params)
            else:
                print "adaptive and fixed links not implemented"
                exit()

    def getLinksKernelPdf(self, params, params0):
        if self.adaptive == 0:
            if self.nfixed == -1:
                return link_stats.getLinksKernelPdf_1(self, params, params0) ## looks at individual links
            else:
                return link_stats.getLinksKernelPdf_2(self, params, params0) ## looks at all the links together
        else:
            if self.nfixed == -1:
                if self.adaptive == 1: return link_stats.getLinksKernelPdf_1_adaptive(self, params, params0)
                #if self.adaptive == 2: return link_stats.getLinksKernelPdf_1_logistic(self, params, params0)
            else:
                print "adaptive and fixed links not implemented"
                exit()
            
       
        

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

        if self.adaptive == 2:
            # extract all the link info from the data
            x = numpy.array( population[:, self.locations ] )
            nparticles = numpy.shape(x)[0]

            tmp=list()
            tmp.append([])
            tmp.append([])
            tmp.append([])
            
            # estimate probabilities of S1
            type = 0
            obs = x[:,0]
            pm1 = len( numpy.where(obs == -1)[0] )/float(nparticles)
            p0  = len( numpy.where(obs ==  0)[0] )/float(nparticles)
            pp1 = len( numpy.where(obs == +1)[0] )/float(nparticles)
            tmp[0].append(self.prior_type[0])
            tmp[1].append(type)
            tmp[2].append([pm1, p0, pp1])
          
            # need to do nlink-1 regressions S_i | S_(i-1), S_(i-2), S_(i-3) .... 0
            for i in range(1,self.nlinks):
                obs = x[:,i] # should be 0, 1
                expl = x[:,range(0,i)]
                expl = sm.add_constant(expl, prepend=False)
                # print i+1, numpy.arange(0,i+1), numpy.shape(obs),  numpy.shape(expl)

                #print "obs:"
                #for j in range(len(obs) ):
                #    print obs[j], "\t", expl[j,:]
                #pass 

                try:
                    glm_binom = sm.GLM(obs, expl, family=sm.families.Binomial())
                    res = glm_binom.fit()
                    # print "logistic regression params:", res.params

                    tmp[0].append( self.prior_type[0] )
                    tmp[1].append( 1 )
                    tmp[2].append( res.params )
                    
                except:
                    type = 0
                    pm1 = len( numpy.where(obs == -1)[0] )/float(nparticles)
                    p0  = len( numpy.where(obs ==  0)[0] )/float(nparticles)
                    pp1 = len( numpy.where(obs == +1)[0] )/float(nparticles)
                    tmp[0].append(self.prior_type[0])
                    tmp[1].append(type)
                    tmp[2].append([pm1, p0, pp1])

            self.kernel = tmp 
            print self.kernel

    
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
    # FIXED EDGES WITH ENUMERATED PRIOR
    def sample_links_from_prior_fixed_em(self, reti):
   
        # sample a number from 1 to prior_size
        u = numpy.random.randint(0,self.prior_size,1)

        # the links start at third column
        links = [ self.prior_data[u][i+2] for i in range(self.nlinks)]
        
        # print 'sample from prior:', links
        # loop over and fill links from the precomputed prior

        for i in range(self.nlinks):
            reti[ self.locations[i] ] = links[i]

        return reti

    
    def getLinksPriorPdf_fixed_em(self, params):

        # In this special case we have constructed prior
        # such that all hard work has been done
        prior_prob = 1/float(self.prior_size)

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

    def perturbLinks_prior_em(self,params):
        
        # change the links with probability linkp
        u = rnd.uniform(0,1)
        if u < self.linkp:
            params = link_stats.sample_links_from_prior_fixed_em(self, params)
        
        return params

    # This is an independent kernel, effectively sampling links from the previous population
    # There is no dependence on the current particle
    def perturbLinks_1_adaptive(self, params):

        links = [0 for i in range(self.nlinks)]
        newlinks = [0 for i in range(self.nlinks)]

        ss = [-1, 0, +1]

        for i in range(self.nlinks):
            links[i] = params[ self.locations[i] ]
            kdist = self.kernel[i]

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

    # This is an independent kernel, effectively sampling links from the previous population
    # There is no dependence on the current particle
    def perturbLinks_1_logistic(self, params):
        # print "\n#### sampling from logistic model"

        links = [params[ self.locations[i]] for i in range(self.nlinks)]
        newlinks = [0 for i in range(self.nlinks)]

        ss = [-1, 0, +1]

        # sample the first link according to p
        kdist = self.kernel[2][0]
        newlinks[0] = ss[ numpy.where( numpy.random.multinomial(1, kdist ) == 1 )[0] ]
      
        # here i goes from 1 to nlinks
        for i in range(1,self.nlinks):

            prior = self.kernel[0][i]
            type = self.kernel[1][i]
            # print prior, type
            
            if type == 0:
                # just standard sample from p
                kdist = self.kernel[2][i]
                newlinks[i] = ss[ numpy.where( numpy.random.multinomial(1, kdist ) == 1 )[0] ]

            elif type == 1:
                # here we are fitting link S_(i) | S_(i-1), S_(i-2) .... 
                #  S1 | S0
                #  S2 | S1, S0
                #  S3 | S2, S1, S0
                
                # sample from single logistic regression
                beta  = numpy.matrix( self.kernel[2][i] )
                # there are i regressors + constant
                # create a column vector, add the previous samples and the constant term
                x = numpy.matrix( numpy.zeros([i+1, 1]) )
                x[0,:] =  newlinks[0:i]
                x[1,:] =  1
                #print "x:\n", x
                #print numpy.shape(x)
                
                # obtain theoretical p(x) value p(Y = 1 | X = x) = 1/(1 + exp( -( BX ) )
                # print  "beta:", beta, numpy.shape(beta), beta*x
                
                px = 1/(1 + numpy.exp( - (beta*x)[0,0] ) )
                # print px, numpy.where( rnd.multinomial(1, [ 1-px, px ] ) == 1 )[0]

                if prior ==  2 : v = [0, 1]
                if prior == -2 : v = [0, -1]
                newlinks[i] = v[ numpy.where( rnd.multinomial(1, [ 1-px, px ] ) == 1 )[0] ]
                
            else:
                # sample from two logistic regressions
                pass
 
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

    # no constraints and adaptive (independent of params0)
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

    # no constraints and adaptive (independent of params0)
    def getLinksKernelPdf_1_logistic(self, params, params0):

        links = [params[ self.locations[i]] for i in range(self.nlinks)]
        prob=1

        # these need to be compared link by link
        for i in range(self.nlinks):
            prior = self.kernel[0][i]
            type = self.kernel[1][i]
            ## print prior, type
            
            if type == 0:
                # just standard sample from p
                kdist = self.kernel[2][i]
                ind = links[i] + 1 # 0, 1, 2
                kern = kdist[ ind ]
                prob=prob*kern

            elif type == 1:
                # here we are calculating p( S_(i) | S_(i-1), S_(i-2) ....) 
                #  p( S1 | S0 )
                #  p( S2 | S1, S0)
                #  p( S3 | S2, S1, S0)
                
                # sample from single logistic regression
                beta  = numpy.matrix( self.kernel[2][i] )
                # there are i regressors + constant
                # create a column vector, add the previous samples and the constant term
                x = numpy.matrix( numpy.zeros([i+1, 1]) )
                x[0,:] =  links[0:i]
                x[1,:] =  1
                #print "x:\n", x
                #print numpy.shape(x)
                
                # obtain theoretical p(x) value p(Y = 1 | X = x) = 1/(1 + exp( -( BX ) )
                px = 1/(1 + numpy.exp( - (beta*x)[0,0] ) )
                
                if links[i] == 0:
                    prob=prob*(1-px)
                else:
                    prob=prob*px
         
        #print "pdf kernel:", prob
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



