# writer functions

import numpy as np
import string, re

class input_writer:
    def __init__(self, network, xmlfile, times, model_name):

        self.xmlfile = open(xmlfile,'w')

        print >>self.xmlfile, "<input>"
        print >>self.xmlfile, "<modelnumber> 1 </modelnumber>"
        print >>self.xmlfile, "<restart> False </restart>"
        print >>self.xmlfile, "<autoepsilon>"
        print >>self.xmlfile, "<finalepsilon> 100.0 </finalepsilon>"
        print >>self.xmlfile, "<alpha> 0.9 </alpha>"
        print >>self.xmlfile, "</autoepsilon>"
        print >>self.xmlfile, "<particles> 1000 </particles>"
        print >>self.xmlfile, "<beta> 1 </beta>"
        print >>self.xmlfile, "<dt> 0.01 </dt>"
        print >>self.xmlfile, "<kernel> uniform </kernel>"
        print >>self.xmlfile, "<modelkernel> 0.7 </modelkernel>"
 
        print >>self.xmlfile, "<data>\n <times>"
        for t in np.arange(times[0],times[1],times[2]):
            print >>self.xmlfile, t, " ",
        print >>self.xmlfile, "</times>"
        print >>self.xmlfile, "<variables> <var1> </var1> </variables>"
        print >>self.xmlfile, "</data>"


        print >>self.xmlfile, "<models>"
        print >>self.xmlfile, "<model1>\n"
        print >>self.xmlfile, "<name> "+ model_name + " </name>"
        print >>self.xmlfile, "<source> </source>"
        print >>self.xmlfile, "<type> ODE </type>"
        print >>self.xmlfile, "<fit> species6 </fit>"
        print >>self.xmlfile, "<logp> True </logp>\n"

        print >>self.xmlfile, "<initial>"
        for s in network.species:
            print >>self.xmlfile,  "<"+ s +"> constant 1 </"+ s +">"
        print >>self.xmlfile,"</initial>\n"

        
        print >>self.xmlfile, "<parameters>"
        for p in network.params:
            print >>self.xmlfile,  "<"+ p +"> uniform -3 2 </"+ p +">"
        print >>self.xmlfile,"</parameters>\n"

        print >>self.xmlfile, "</model1>"
        print >>self.xmlfile, "</models>"
        print >>self.xmlfile, "</input>"

        
        self.xmlfile.close()


class base_code_writer:
    def __init__(self, network, outfile):
        self.fout = open(outfile,'w')

        self.network = network
       
    def write_laws(self):

        for l in self.network.laws:
            print >>self.fout, "__device__ double "+l.name+"(",

            ni = len(l.inputs)
            for i in range(ni):
                print >>self.fout, "double "+l.inputs[i],
                if i < ni-1:
                    print >>self.fout, ",",
            print >>self.fout, "){"
            print >>self.fout, "\t return "+ l.output + ";"
            print >>self.fout, "}\n"


    def write_header(self):
        print >>self.fout, "#define NSPECIES " + repr(self.network.nspecies)
        print >>self.fout, "#define NPARAM " + repr(self.network.npar)
        print >>self.fout, "#define NREACT " + repr(self.network.sdim[0])

        print >>self.fout, "\n\n"
        print >>self.fout, "#define geq(a,b) a>=b"
        #print >>self.fout, "#define par(a) tex2D(param_tex,a,tid)"
        print >>self.fout, "\n"
        print >>self.fout, "// if v is -1 0 +1 then switchp(v) is 0 0 1"
        print >>self.fout, "#define switchp(v) 0.5*v*(v+1)"
        print >>self.fout, "// if v is -1 0 +1 then switchn(v) is 1 0 0"
        print >>self.fout, "#define switchn(v) 0.5*v*(v-1)"
        print >>self.fout, "\n"

    def write_params(self):
        for p in range(len(self.network.params)):
            print >>self.fout, "#define "+  self.network.params[p] + " tex2D(param_tex,"+repr(p)+",tid)"
        print >>self.fout, "\n"

    def finish(self):
        self.fout.close()

# helper function to replace the species in the rates
def fill_species(network,rate):

    # loop over all species in the network
    # replace the occurances in rate
    ret = rate
    #print "\n"
    #print "\t", rate
    for i in range(network.nspecies):
        cs = "y["+repr(i)+"]";
        #print "\t", i, "replace", network.species[i], "with", cs,

        # old way
        #ret = string.replace(ret, network.species[i], cs)

        # new way
        #p = re.compile("(?<!-)"+network.species[i]+"(?!-)")
        #ret = p.sub(cs,ret);

        # new way 2
        p = re.compile("(?<=\s)"+network.species[i]+"(?=\s)")
        ret = p.sub(cs,ret);

        #print "to get", ret

    return ret

class mjp_code_writer(base_code_writer):

    def __init__(self, network, outfile):         

        # initialise using the base class
        base_code_writer.__init__(self, network, outfile)
        base_code_writer.write_header(self)
        base_code_writer.write_params(self)
        base_code_writer.write_laws(self)
       
        print >>self.fout, "__constant__ int smatrix[]={"
        print >>self.fout, "//",
        for i in range(network.nspecies):
             print >>self.fout, "\t"+network.species[i],
        print >>self.fout, ""

        for i in range(network.sdim[0]):
            for j in range(network.sdim[1]):
                print >>self.fout, "\t"+repr(network.smatrix[i,j])+",",
            print >>self.fout, ""
        print >>self.fout, "};"

        print >>self.fout, "\n\n"

        print >>self.fout, "__device__ void hazards(int *y, float *h, float t, int tid){"
        #print >>self.fout, "\tif( t >= 50 && t <= 51){"
        #print >>self.fout, "\t\ty[0]=60;"
        #print >>self.fout, "\t}"

        print >>self.fout, "\t// reactions //\n"
        for i in range(network.nreactions):
            print >>self.fout, "\t// R"+repr(i)+":",
            network.total_reactions[i].fdisplay(self.fout)
            print >>self.fout, "\th["+repr(i)+"] = ",
            print >>self.fout, fill_species(network, network.total_reactions[i].rate ),
            print >>self.fout, ";\n"

        print >>self.fout, "}"  

        base_code_writer.finish(self)


class ode_code_writer(base_code_writer):

    def __init__(self, network, outfile):         

        # initialise using the base class
        base_code_writer.__init__(self, network, outfile)
        base_code_writer.write_header(self)
        base_code_writer.write_params(self)
        base_code_writer.write_laws(self)
        
       
        print >>self.fout, "struct myFex{"
        print >>self.fout, "__device__ void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/){"
        print >>self.fout, "\tint tid = blockDim.x * blockIdx.x + threadIdx.x;"
        #print >>self.fout, "\tif( t[0] >= 100){"
        #print >>self.fout, "\t\ty[0]=0.6;"
        #print >>self.fout, "\t}"

        # loop over the stoichiometry matrix
        smt = np.transpose( self.network.smatrix )
       
        reaction_count = 0
        for i in range(network.nspecies):
            print >>self.fout, "\tydot["+repr(i)+"] = ",

            # get the non zero contributions from the S matrix
            r = np.nonzero( smt[i,:] )[0]
            nr = len( r )
            for j in range(nr):
                print >>self.fout, "("+repr(smt[i,r[j]] )+")*"+ fill_species(network,network.total_reactions[r[j]].rate),
 
                if j != nr-1:
                    print >>self.fout, " + ",
                
            print >>self.fout, ";"

        print >>self.fout, "\n"
        print >>self.fout, "\t}"
        print >>self.fout, "};"
        
        jac_string = """
struct myJex{
       __device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd/*, void *otherData*/){
                return; 
       }
};
"""

        print >>self.fout, jac_string
        
        base_code_writer.finish(self)

        
