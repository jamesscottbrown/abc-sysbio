# sampling link functions

import numpy
from numpy import random as rnd

# define class to hold link information
class link_stats:

    def __init__(self, model, linkp, nfixed):

        self.linkp = linkp
        self.nfixed = nfixed
        self.nparameters = model.nparameters

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

            if( total_links == self.nfixed ):
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
        if total_links == self.nfixed:
            prior_prob = 1.0
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
            ##print prior_prob,

        print ""
        return prior_prob

       
####################### KERNELS

def perturbLinks_1(kernel, params, priors, link_info):
    # Here every link is perturbed independently

    ##print "link:",
    for n in kernel[0]:
        if priors[n][0] == 4:
                    
            # ret is the returned link value
            ret = params[n]

            s = 0
            if priors[n][1] == 3:
                s = set([-1, 0, 1])

            elif priors[n][1] == 2:
                s = set([0, 1])

            else:
                s = set([-1, 0])
    
            # change the link with probability linkp
            u = rnd.uniform(0,1)
            if u < link_info.linkp:
                ss = set( [params[n]] )
                ar = numpy.array( list(s-ss) )
                rnd.shuffle( ar )
                ret = ar[0]

            ##print params[n], "->", ret, 
            params[n] = ret

    ##print ""
    return params

def perturbLinks_2(kernel, params, priors, link_info):
    # Here only one link is perturbed independently
    
    # change one link with probability linkp
    u = rnd.uniform(0,1)
    if u < link_info.linkp :
        # choose a link to perturb uniformly
        v = [i for i in range(link_info.nlinks) ]
        p = [1/float(nlinks) for i in range(links_info.nlinks) ]
        linknum = v[ numpy.where(rnd.multinomial(n=1,pvals=p,size=1)[0]==1)[0][0] ]

        #print "linknum:", linknum
        # loop over to modify link num
        count = 0
        for n in kernel[0]:
            if priors[n][0] == 4:
                if linknum == count:
                    # perturb this link

                    s = 0
                    if priors[n][1] == 3:
                        s = set([-1, 0, 1])

                    elif priors[n][1] == 2:
                        s = set([0, 1])

                    else:
                        s = set([-1, 0])

                    ss = set( [params[n]] )
                    ar = numpy.array( list(s-ss) )
                    rnd.shuffle( ar )
                    ret = ar[0]

                    #print linknum, params[n], "->", ret 
                    params[n] = ret

                count += 1
        
    return params

## get probability of all links given another set of links and the the kernel
def getLinksKernelPdf_1(kernel, params, params0, priors, link_info):
    prob=1
    ind=0
    for n in kernel[0]:
        if priors[n][0] == 4:
            # p( t | t -1 ) =  1-linkp (t == t-1), linkp (t != t-1 )
            kern = 0;
            if params[n] == params0[n]:
                kern = 1-link_info.linkp
            else:
                kern = link_info.linkp
           
            prob=prob*kern
    ##print "pdf kernel:", prob
    return prob

def getLinksKernelPdf_2(kernel, params, params0, priors, link_info):
    
    nochange = True
    prob = (1-link_info.linkp)
    for n in kernel[0]:
        if priors[n][0] == 4:
            if not (params[n] == params0[n]):
                nochange = False

    if nochange == False:
        prob = link_info.linkp
    
    ##print "pdf kernel:", prob
    return prob




## calculate the probability of all the links under the prior
##def getLinksPriorPdf_fixed(params, priors):  
##  return prior_prob

# set the defaults
def perturbLinks(kernel, params, priors, linkp):
    return perturbLinks_1(kernel, params, priors, linkp)
    ##return perturbLinks_2(kernel, params, priors, linkp, nlinks=2)

def getLinksKernelPdf(kernel, params, params0, priors, linkp):
    return getLinksKernelPdf_1(kernel, params, params0, priors, linkp)
    ##return getLinksKernelPdf_2(kernel, params, params0, priors, linkp)


