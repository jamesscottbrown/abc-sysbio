struct myFex
{
        __device__ void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/)
        {
          int tid = blockDim.x * blockIdx.x + threadIdx.x; 
          
	  ydot[0]= tex2D(param_tex, 0, tid) - y[0]*tex2D(param_tex, 1, tid); 
	}
};

struct myJex
{
        __device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd/*, void *otherData*/)
        {
          return;
        }
};
