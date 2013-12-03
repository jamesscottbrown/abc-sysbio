#define NSPECIES 4
#define NPARAM 22
#define NREACT 12

#include <stdio.h>

#define geq(a,b) a>=b

__device__ double switchp( float v ){
return (double) 0.5*v*(v+1);
}
// if v is -1 0 +1 then switchn(v) is 1 0 0
__device__ double switchn( float v ){
        return (double) 0.5*v*(v-1);
}

#define k1_1 tex2D(param_tex,0,tid)
#define k1_2 tex2D(param_tex,1,tid)

#define Rm_1 tex2D(param_tex,2,tid)
#define Rm_2 tex2D(param_tex,3,tid)

#define sigma_1s tex2D(param_tex,4,tid)
#define sigma_2s tex2D(param_tex,5,tid)

#define sigma_1t1 tex2D(param_tex,6,tid)
#define sigma_2t1 tex2D(param_tex,7,tid)

#define delta_1s tex2D(param_tex,8,tid)
#define delta_2s tex2D(param_tex,9,tid)

#define delta_1t1 tex2D(param_tex,10,tid)
#define delta_2t1 tex2D(param_tex,11,tid)

#define tlA tex2D(param_tex,12,tid)
#define tlB tex2D(param_tex,13,tid)

#define degA tex2D(param_tex,14,tid)
#define degB tex2D(param_tex,15,tid)

#define degmA tex2D(param_tex,16,tid)
#define degmB tex2D(param_tex,17,tid)

#define S1 tex2D(param_tex,18,tid)
#define S2 tex2D(param_tex,19,tid)
#define S3 tex2D(param_tex,20,tid)
#define S4 tex2D(param_tex,21,tid)

//Code for texture memory

// promoter 1 
__device__ double fa( double A , double B , double tid ){
	return Rm_1*(k1_1 + switchp(S1)*sigma_1s*A*A + switchp(S3)*sigma_1t1*B*B)/(1 + k1_1 + switchp(S1)*sigma_1s*A*A + switchn(S1)*delta_1s*A*A + switchp(S3)*sigma_1t1*B*B + switchn(S3)*delta_1t1*B*B );
}

// promoter 2
__device__ double fb( double A , double B , double tid ){
	return Rm_2*(k1_2 + switchp(S2)*sigma_2s*B*B + switchp(S4)*sigma_2t1*A*A)/(1 + k1_2 + switchp(S2)*sigma_2s*B*B + switchn(S2)*delta_2s*B*B + switchp(S4)*sigma_2t1*A*A + switchn(S4)*delta_2t1*A*A);
}

__device__ void step(float *y, float t, unsigned *rngRegs, int tid){

	float W0 = randNormal(rngRegs,sqrt(DT));
	float W1 = randNormal(rngRegs,sqrt(DT));
	float W2 = randNormal(rngRegs,sqrt(DT));
	float W3 = randNormal(rngRegs,sqrt(DT));
	float W4 = randNormal(rngRegs,sqrt(DT));
	float W5 = randNormal(rngRegs,sqrt(DT));
	float W6 = randNormal(rngRegs,sqrt(DT));
	float W7 = randNormal(rngRegs,sqrt(DT));

	float dy0 = DT *(  (1.0)* fa( y[1] , y[3] , tid )   +  (-1.0)* degmA * y[0]  );
	float dy1 = DT *(  (1.0)* tlA * y[0]   +  (-1.0)* degA * y[1]  );
	float dy2 = DT *(  (1.0)* fb( y[1] , y[3] , tid )   +  (-1.0)* degmB * y[2]  );
	float dy3 = DT *(  (1.0)* tlB * y[2]   +  (-1.0)* degB * y[3]  );


	dy0 += (1.0)*sqrt( fa( y[1] , y[3] , tid ) )*W0  +  (-1.0)*sqrt( degmA * y[0] )*W1 ;
	dy1 += (1.0)*sqrt( tlA * y[0] )*W2  +  (-1.0)*sqrt( degA * y[1] )*W3 ;
	dy2 += (1.0)*sqrt( fb( y[1] , y[3] , tid ) )*W4  +  (-1.0)*sqrt( degmB * y[2] )*W5 ;
	dy3 += (1.0)*sqrt( tlB * y[2] )*W6  +  (-1.0)*sqrt( degB * y[3] )*W7 ;

	y[0] += dy0;
	y[1] += dy1;
	y[2] += dy2;
	y[3] += dy3;
}
