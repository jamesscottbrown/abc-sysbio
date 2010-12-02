#define NRMAX 2

__constant__ int smatrix[] = {1, -1};

__device__ void stoichiometry(int *x, int r, int tid){
  x[0] = x[0] + smatrix[r];
}

__device__ void hazards(int *x, float *h, int tid){
  
  if( tex1Dfetch(model_tex, tid) == 0 ){
    h[0] = tex2D(param_tex, 0, tid);
    h[1] = x[0] > 0 ? tex2D(param_tex, 1, tid)*x[0] : 0;
  }
}

__device__ void rules(int *x, float t){ }
