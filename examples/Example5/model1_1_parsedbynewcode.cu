#define NSPECIES 6
#define NPARAM 6
#define NREACT 12

//Code for texture memory
__device__ float function_1(float a1,float a2,float a3,float a4){
    return a1 / (1 + __powf(a2, a3)) + a4;
}


__device__ void step(float *y, float t, unsigned *rngRegs, int tid){

    float d_y5= DT * ((1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[4])-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[5]))/tex2D(param_tex,0,tid));
    float d_y4= DT * ((1.0*(tex2D(param_tex,0,tid)*function_1(tex2D(param_tex,4,tid),y[3],tex2D(param_tex,2,tid),tex2D(param_tex,1,tid)))-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,5,tid)*y[4]))/tex2D(param_tex,0,tid));
    float d_y3= DT * ((1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[2])-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[3]))/tex2D(param_tex,0,tid));
    float d_y2= DT * ((1.0*(tex2D(param_tex,0,tid)*function_1(tex2D(param_tex,4,tid),y[1],tex2D(param_tex,2,tid),tex2D(param_tex,1,tid)))-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,5,tid)*y[2]))/tex2D(param_tex,0,tid));
    float d_y1= DT * ((1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[0])-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[1]))/tex2D(param_tex,0,tid));
    float d_y0= DT * ((1.0*(tex2D(param_tex,0,tid)*function_1(tex2D(param_tex,4,tid),y[5],tex2D(param_tex,2,tid),tex2D(param_tex,1,tid)))-1.0*(tex2D(param_tex,0,tid)*tex2D(param_tex,5,tid)*y[0]))/tex2D(param_tex,0,tid));


    d_y0 += ((1.0*sqrt(tex2D(param_tex,0,tid)*function_1(tex2D(param_tex,4,tid),y[5],tex2D(param_tex,2,tid),tex2D(param_tex,1,tid)))*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(tex2D(param_tex,0,tid)*tex2D(param_tex,5,tid)*y[0])*randNormal(rngRegs,sqrt(DT)))/tex2D(param_tex,0,tid));
    d_y1 += ((1.0*sqrt(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[0])*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[1])*randNormal(rngRegs,sqrt(DT)))/tex2D(param_tex,0,tid));
    d_y2 += ((1.0*sqrt(tex2D(param_tex,0,tid)*function_1(tex2D(param_tex,4,tid),y[1],tex2D(param_tex,2,tid),tex2D(param_tex,1,tid)))*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(tex2D(param_tex,0,tid)*tex2D(param_tex,5,tid)*y[2])*randNormal(rngRegs,sqrt(DT)))/tex2D(param_tex,0,tid));
    d_y3 += ((1.0*sqrt(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[2])*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[3])*randNormal(rngRegs,sqrt(DT)))/tex2D(param_tex,0,tid));
    d_y4 += ((1.0*sqrt(tex2D(param_tex,0,tid)*function_1(tex2D(param_tex,4,tid),y[3],tex2D(param_tex,2,tid),tex2D(param_tex,1,tid)))*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(tex2D(param_tex,0,tid)*tex2D(param_tex,5,tid)*y[4])*randNormal(rngRegs,sqrt(DT)))/tex2D(param_tex,0,tid));
    d_y5 += ((1.0*sqrt(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[4])*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(tex2D(param_tex,0,tid)*tex2D(param_tex,3,tid)*y[5])*randNormal(rngRegs,sqrt(DT)))/tex2D(param_tex,0,tid));

    y[0] += d_y0;
    y[1] += d_y1;
    y[2] += d_y2;
    y[3] += d_y3;
    y[4] += d_y4;
    y[5] += d_y5;
}
//Code for shared memory
__device__ float function_1(float a1,float a2,float a3,float a4){
    return a1 / (1 + __powf(a2, a3)) + a4;
}


__device__ void step(float *parameter, float *y, float t, unsigned *rngRegs){

    float d_y0= DT * ((1.0*(parameter[0]*function_1(parameter[4],y[5],parameter[2],parameter[1]))-1.0*(parameter[0]*parameter[5]*y[0]))/parameter[0]);
    float d_y1= DT * ((1.0*(parameter[0]*parameter[3]*y[0])-1.0*(parameter[0]*parameter[3]*y[1]))/parameter[0]);
    float d_y2= DT * ((1.0*(parameter[0]*function_1(parameter[4],y[1],parameter[2],parameter[1]))-1.0*(parameter[0]*parameter[5]*y[2]))/parameter[0]);
    float d_y3= DT * ((1.0*(parameter[0]*parameter[3]*y[2])-1.0*(parameter[0]*parameter[3]*y[3]))/parameter[0]);
    float d_y4= DT * ((1.0*(parameter[0]*function_1(parameter[4],y[3],parameter[2],parameter[1]))-1.0*(parameter[0]*parameter[5]*y[4]))/parameter[0]);
    float d_y5= DT * ((1.0*(parameter[0]*parameter[3]*y[4])-1.0*(parameter[0]*parameter[3]*y[5]))/parameter[0]);


    d_y0+= ((1.0*sqrt(parameter[0]*function_1(parameter[4],y[5],parameter[2],parameter[1]))*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(parameter[0]*parameter[5]*y[0])*randNormal(rngRegs,sqrt(DT)))/parameter[0]);
    d_y1+= ((1.0*sqrt(parameter[0]*parameter[3]*y[0])*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(parameter[0]*parameter[3]*y[1])*randNormal(rngRegs,sqrt(DT)))/parameter[0]);
    d_y2+= ((1.0*sqrt(parameter[0]*function_1(parameter[4],y[1],parameter[2],parameter[1]))*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(parameter[0]*parameter[5]*y[2])*randNormal(rngRegs,sqrt(DT)))/parameter[0]);
    d_y3+= ((1.0*sqrt(parameter[0]*parameter[3]*y[2])*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(parameter[0]*parameter[3]*y[3])*randNormal(rngRegs,sqrt(DT)))/parameter[0]);
    d_y4+= ((1.0*sqrt(parameter[0]*function_1(parameter[4],y[3],parameter[2],parameter[1]))*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(parameter[0]*parameter[5]*y[4])*randNormal(rngRegs,sqrt(DT)))/parameter[0]);
    d_y5+= ((1.0*sqrt(parameter[0]*parameter[3]*y[4])*randNormal(rngRegs,sqrt(DT))-1.0*sqrt(parameter[0]*parameter[3]*y[5])*randNormal(rngRegs,sqrt(DT)))/parameter[0]);

    y[0] += d_y0;
    y[1] += d_y1;
    y[2] += d_y2;
    y[3] += d_y3;
    y[4] += d_y4;
    y[5] += d_y5;
}
