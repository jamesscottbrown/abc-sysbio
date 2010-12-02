// SIRmodels.cu

struct myFex
{
	__device__ void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/)
	{
	  int tid = blockDim.x * blockIdx.x + threadIdx.x;
	  int model = tex1Dfetch(model_tex, tid);
	   
	  if(model == 0){
	    ydot[0] =  tex2D(param_tex, 4, tid)      - tex2D(param_tex, 1, tid)*y[0]      - tex2D(param_tex, 3, tid)*y[0]*y[1];
	    ydot[1] = -tex2D(param_tex, 1, tid)*y[1] + tex2D(param_tex, 3, tid)*y[0]*y[1] - tex2D(param_tex, 2, tid)*y[1];
	    ydot[2] = -tex2D(param_tex, 1, tid)*y[2] + tex2D(param_tex, 2, tid)*y[1];
	  }
	  
	  if(model == 1){
	    ydot[0] =  tex2D(param_tex, 4, tid)      - tex2D(param_tex, 1, tid)*y[0]      - tex2D(param_tex, 3, tid)*y[0]*y[2];
	    ydot[1] = -tex2D(param_tex, 1, tid)*y[1] + tex2D(param_tex, 3, tid)*y[0]*y[2] - tex2D(param_tex, 5, tid)*y[1];
	    ydot[2] = -tex2D(param_tex, 1, tid)*y[2] + tex2D(param_tex, 5, tid)*y[1]      - tex2D(param_tex, 2, tid)*y[2];
	    ydot[3] = -tex2D(param_tex, 1, tid)*y[3] + tex2D(param_tex, 2, tid)*y[2];
	  }

	  if(model == 2){
	    ydot[0] =  tex2D(param_tex, 4, tid)      - tex2D(param_tex, 1, tid)*y[0]      - tex2D(param_tex, 3, tid)*y[0]*y[1] + tex2D(param_tex, 5, tid)*y[2];
	    ydot[1] = -tex2D(param_tex, 1, tid)*y[1] + tex2D(param_tex, 3, tid)*y[0]*y[1] - tex2D(param_tex, 2, tid)*y[1];
	    ydot[2] = -tex2D(param_tex, 1, tid)*y[2] + tex2D(param_tex, 2, tid)*y[1]      - tex2D(param_tex, 5, tid)*y[2];
	  }
	}
};

struct myJex
{
	__device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd/*, void *otherData*/)
	{
	  return;
	}
};
