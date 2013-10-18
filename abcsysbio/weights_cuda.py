
import pycuda
import numpy as np
import pycuda.driver as driver
from pycuda.compiler import SourceModule
from abcsysbio import statistics

def fp(v):
    # print the list, v, as a string
    return str(v).strip('[]')


class weights_cuda:

    # this only works when there is one model. prior and kernels are for the single model

    def __init__(self, nparticles, nparam, prior, kernel, link_info, parameters_curr, parameters_prev, weights_prev, weights_out):

        if nparticles % 64 != 0:
            print "number of particles must be divisible by 64"
            exit()


        # get some required information
        uniform_params = []
        links_params = []
        kernel_vals = []

        for n in range(nparam):
            if prior[n][0] == 2: 
                uniform_params.append(n)
            if prior[n][0] == 4: 
                links_params.append(n)

        ind = 0 # here ind refers to the non constant parameters (but includes links)
        for n in kernel[0]:
	    if prior[n][0] == 4:
                # skip but increment ind
                ind += 1
            else:
                kernel_vals.append( kernel[2][ind][1] )
                ind += 1


        predef  = ""
        predef += "\n#define NPARTICLES "+repr(nparticles)
        predef += "\n#define NPARAM "+repr(nparam)
        predef += "\n#define NUNIF "+repr( len(uniform_params) )
        predef += "\n__device__ const int uniform_pars["+repr(len(uniform_params))+"] = {"+fp(uniform_params)+"};"
        predef += "\n__device__ const double unif_kern["+repr(len(kernel_vals))+"] = {"+fp(kernel_vals)+"};"

        if link_info.nlinks > 0:
            predef += "\n#define LINKP "+repr(link_info.linkp)
            predef += "\n#define NLINKS "+repr( len(links_params) )
            predef += "\n__device__ const int link_locations["+repr(len(links_params))+"] = {"+fp(links_params)+"};"
       

        parameter_unif_pdf_code="""

        __device__ double dunif(double x, double a, double b)
        {
            if( x >= a && x <= b) return 1/(b-a);
            else return 0;
        }

        """
        
        parameter_kernel_pdf_code="""
        
        // this function take pointers to the first elements of the parameter arrays
        // j is previous particle, i is current particle
        // this function is for uniform kernels only
        __device__ double parameter_kernel_pdf(double* par_j, double* par_i)
        {
            double p = 1;
            for(int k=0; k<NUNIF; k++){
                p *= dunif( par_i[ uniform_pars[k]], par_j[ uniform_pars[k]] - unif_kern[k], par_j[ uniform_pars[k]] + unif_kern[k] );
            }
            return p;
        }

        """

        link_kernel_pdf_code_1="""
        
        // this functions take pointers to the first elements of the parameter arrays
        __device__ double link_kernel_pdf(double* par_j, double* par_i)
        {
            // this function is for a free edge model. We just at the change of each link
            double prob = 1;
            double s = 0;

            for(int k=0; k<NLINKS; k++){
                s = ( par_j[ link_locations[k] ] - par_i[ link_locations[k] ] );
                if( s < 1e-3 ) prob = prob*(1-LINKP);
                else  prob = prob*LINKP;
            }
            return prob;
        }

        """

        link_kernel_pdf_code_2="""
        
        // this functions take pointers to the first elements of the parameter arrays
        __device__ double link_kernel_pdf(double* par_j, double* par_i)
        {
            // this function is for a constrained edge model. We just look to see if there is any difference in the edges
            double s = 0;
            for(int k=0; k<NLINKS; k++){
                s += par_j[ link_locations[k] ] - par_i[ link_locations[k] ];
            }

            if( s < 1e-3 ) return (double) 1-LINKP;
            else return (double) LINKP;
        }

        """

        link_kernel_pdf_code_null="""
        
        // this function is a place holder for when there are no free links
        __device__ double link_kernel_pdf(double* par_j, double* par_i)
        {
           return 1.0;
        }

        """
        

        weight_calculation_code="""
        
        // This kernel calculate the denominator of the weight assuming there is only one model
        __global__ void calculate_weights(double* weight_prev, double* param_prev, double* param_curr, double* weight)
        {
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;

            // weight_prev, param_prev, param_curr could go into shared memory to speed things up

            double w = 0;
            for(int j=0; j<NPARTICLES; ++j){
                w += weight_prev[j]*parameter_kernel_pdf( &(param_prev[NPARAM*j]) , &(param_curr[NPARAM*tid]) )*link_kernel_pdf( &(param_prev[NPARAM*j]), &(param_curr[NPARAM*tid]) );
            }

            //printf("weight : %d %f \\n", tid, w);
            weight[ tid ] = w;
        }

        """

        # Build the code
        
        if link_info.nfixed == -1:
            print "weights_cuda : running free link model"
            link_kernel_pdf_code = link_kernel_pdf_code_1  ## looks at individual links
        else:
            print "weights_cuda : running constrained link model"
            link_kernel_pdf_code = link_kernel_pdf_code_2  ## looks at all the links together

        
        if link_info.nlinks > 0:
            all_code = predef + parameter_unif_pdf_code + parameter_kernel_pdf_code + link_kernel_pdf_code + weight_calculation_code
        else:
            all_code = predef + parameter_unif_pdf_code + parameter_kernel_pdf_code + link_kernel_pdf_code_null + weight_calculation_code
        

        if False:
            of = open("full_weight_kernel.cu","w")
            print >>of, all_code

        compiled = pycuda.compiler.SourceModule( all_code, nvcc="nvcc" )
        kernel = compiled.get_function("calculate_weights")
     


        # Define the local  input arrays
        w_p  = np.zeros( [nparticles], dtype=np.float64)
        p_p  = np.zeros( [nparticles*nparam], dtype=np.float64)
        p_c  = np.zeros( [nparticles*nparam], dtype=np.float64)
        wout = np.zeros( [nparticles], dtype=np.float64)

        # Fill the values
        count_2D = 0
        for k in range(nparticles):
            w_p[k] = weights_prev[k]

            for j in range(nparam):
                p_p[count_2D] = parameters_prev[k][j]
                p_c[count_2D] = parameters_curr[k][j]
                count_2D += 1
            
    
        # allocate on device
        d_w_p  = driver.mem_alloc(w_p.size * w_p.dtype.itemsize)
        d_p_p  = driver.mem_alloc(p_p.size * p_p.dtype.itemsize)
        d_p_c  = driver.mem_alloc(p_c.size * p_c.dtype.itemsize)
        d_wout = driver.mem_alloc(wout.size * wout.dtype.itemsize)
       
        # copy to device
        driver.memcpy_htod(d_w_p, w_p)
        driver.memcpy_htod(d_p_p, p_p)
        driver.memcpy_htod(d_p_c, p_c)
        driver.memcpy_htod(d_wout, wout)
        
        # run
        threads = 64
        blocks = nparticles/threads
        kernel( d_w_p, d_p_p, d_p_c, d_wout, block=(threads,1,1), grid=(blocks,1) )
        driver.memcpy_dtoh(wout, d_wout)

        # calculate weight given the prior and the denominator
        for k in range(nparticles):
            # particle prior probability given by prior(theta) * prior (l)
            pr = 1
            for n in range(nparam):    
                if prior[n][0]==2: 
                    pr = pr * statistics.getPdfUniform(prior[n][1], prior[n][2], parameters_curr[k][n])
                    
            ## multiply by the prior for this set of links
            prior_links = link_info.getLinksPriorPdf( parameters_curr[k] )
            pr = pr*prior_links

            ##print pr, wout[k]
            
            weights_out[k] = pr/ wout[k]
             


## weights_cuda( 100, 10, [], [], 0.7, [], [], [], [] )
           
