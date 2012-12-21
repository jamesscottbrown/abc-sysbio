#define NSPECIES 6
#define NPARAM 6
#define NREACT 12

__device__ float function_1(float a1,float a2,float a3,float a4){
    return a1 / (1 + __powf(a2, a3)) + a4;
}

struct myFex{
    __device__ void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/)
    {
        int tid = blockDim.x * blockIdx.x + threadIdx.x;



        ydot[5]=(1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[4])-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[5]))/tex2D(param_tex,0,tid);
        ydot[4]=(1.0*(tex2D(param_tex,0,tid)*function_1(tex2D(param_tex,4,tid),y[3],tex2D(param_tex,2,tid),tex2D(param_tex,1,tid)))-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,5,tid)*y[4]))/tex2D(param_tex,0,tid);
        ydot[3]=(1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[2])-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[3]))/tex2D(param_tex,0,tid);
        ydot[2]=(1.0*(tex2D(param_tex,0,tid)*function_1(tex2D(param_tex,4,tid),y[1],tex2D(param_tex,2,tid),tex2D(param_tex,1,tid)))-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,5,tid)*y[2]))/tex2D(param_tex,0,tid);
        ydot[1]=(1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[0])-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[1]))/tex2D(param_tex,0,tid);
        ydot[0]=(1.0*(tex2D(param_tex,0,tid)*function_1(tex2D(param_tex,4,tid),y[5],tex2D(param_tex,2,tid),tex2D(param_tex,1,tid)))-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,5,tid)*y[0]))/tex2D(param_tex,0,tid);

    }
};


 struct myJex{
    __device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd/*, void *otherData*/){
        return; 
    }
};