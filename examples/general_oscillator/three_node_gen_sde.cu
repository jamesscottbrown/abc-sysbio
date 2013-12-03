#define NSPECIES 6
#define NPARAM 42
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
#define k1_3 tex2D(param_tex,2,tid)

#define Rm_1 tex2D(param_tex,3,tid)
#define Rm_2 tex2D(param_tex,4,tid)
#define Rm_3 tex2D(param_tex,5,tid)

#define sigma_1s tex2D(param_tex,6,tid)
#define sigma_2s tex2D(param_tex,7,tid)
#define sigma_3s tex2D(param_tex,8,tid)

#define sigma_1t1 tex2D(param_tex,9,tid)
#define sigma_2t1 tex2D(param_tex,10,tid)
#define sigma_3t1 tex2D(param_tex,11,tid)
#define sigma_1t2 tex2D(param_tex,12,tid)
#define sigma_2t2 tex2D(param_tex,13,tid)
#define sigma_3t2 tex2D(param_tex,14,tid)

#define delta_1s tex2D(param_tex,15,tid)
#define delta_2s tex2D(param_tex,16,tid)
#define delta_3s tex2D(param_tex,17,tid)

#define delta_1t1 tex2D(param_tex,18,tid)
#define delta_2t1 tex2D(param_tex,19,tid)
#define delta_3t1 tex2D(param_tex,20,tid)
#define delta_1t2 tex2D(param_tex,21,tid)
#define delta_2t2 tex2D(param_tex,22,tid)
#define delta_3t2 tex2D(param_tex,23,tid)

#define tlA tex2D(param_tex,24,tid)
#define tlB tex2D(param_tex,25,tid)
#define tlC tex2D(param_tex,26,tid)

#define degA tex2D(param_tex,27,tid)
#define degB tex2D(param_tex,28,tid)
#define degC tex2D(param_tex,29,tid)

#define degmA tex2D(param_tex,30,tid)
#define degmB tex2D(param_tex,31,tid)
#define degmC tex2D(param_tex,32,tid)

#define S1 tex2D(param_tex,33,tid)
#define S2 tex2D(param_tex,34,tid)
#define S3 tex2D(param_tex,35,tid)
#define S4 tex2D(param_tex,36,tid)
#define S5 tex2D(param_tex,37,tid)
#define S6 tex2D(param_tex,38,tid)
#define S7 tex2D(param_tex,39,tid)
#define S8 tex2D(param_tex,40,tid)
#define S9 tex2D(param_tex,41,tid)

//Code for texture memory

// promoter 1 
__device__ double fa( double A , double B , double C , double tid ){
	return Rm_1*(k1_1 + switchp(S1)*sigma_1s*A*A + switchp(S5)*sigma_1t1*B*B + switchp(S7)*sigma_1t2*C*C)/(1 + k1_1 + switchp(S1)*sigma_1s*A*A + switchn(S1)*delta_1s*A*A + switchp(S5)*sigma_1t1*B*B + switchn(S5)*delta_1t1*B*B + switchp(S7)*sigma_1t2*C*C + switchn(S7)*delta_1t2*C*C );
}

// promoter 2
__device__ double fb( double A , double B , double C , double tid ){
	return Rm_2*(k1_2 + switchp(S2)*sigma_2s*B*B + switchp(S6)*sigma_2t1*C*C + switchp(S8)*sigma_2t2*A*A)/(1 + k1_2 + switchp(S2)*sigma_2s*B*B + switchn(S2)*delta_2s*B*B + switchp(S6)*sigma_2t1*C*C + switchn(S6)*delta_2t1*C*C + switchp(S8)*sigma_2t2*A*A + switchn(S8)*delta_2t2*A*A);
}

// promoter 3 
__device__ double fc( double A , double B , double C , double tid ){
	return Rm_3*(k1_3 + switchp(S3)*sigma_3s*C*C + switchp(S4)*sigma_3t1*A*A + switchp(S9)*sigma_3t2*B*B)/(1 + k1_3 + switchp(S3)*sigma_3s*C*C + switchn(S3)*delta_3s*C*C + switchp(S4)*sigma_3t1*A*A + switchn(S4)*delta_3t1*A*A + switchp(S9)*sigma_3t2*B*B + switchn(S9)*delta_3t2*B*B);
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
	float W8 = randNormal(rngRegs,sqrt(DT));
	float W9 = randNormal(rngRegs,sqrt(DT));
	float W10 = randNormal(rngRegs,sqrt(DT));
	float W11 = randNormal(rngRegs,sqrt(DT));


	float dy0 = DT *(  (1.0)* fa( y[1] , y[3] , y[5] , tid )   +  (-1.0)* degmA * y[0]  );
	float dy1 = DT *(  (1.0)* tlA * y[0]   +  (-1.0)* degA * y[1]  );
	float dy2 = DT *(  (1.0)* fb( y[1] , y[3] , y[5] , tid )   +  (-1.0)* degmB * y[2]  );
	float dy3 = DT *(  (1.0)* tlB * y[2]   +  (-1.0)* degB * y[3]  );
	float dy4 = DT *(  (1.0)* fc( y[1] , y[3] , y[5] , tid )   +  (-1.0)* degmC * y[4]  );
	float dy5 = DT *(  (1.0)* tlC * y[4]   +  (-1.0)* degC * y[5]  );


	dy0 += (1.0)*sqrt( fa( y[1] , y[3] , y[5] , tid ) )*W0  +  (-1.0)*sqrt( degmA * y[0] )*W1 ;
	dy1 += (1.0)*sqrt( tlA * y[0] )*W2  +  (-1.0)*sqrt( degA * y[1] )*W3 ;
	dy2 += (1.0)*sqrt( fb( y[1] , y[3] , y[5] , tid ) )*W4  +  (-1.0)*sqrt( degmB * y[2] )*W5 ;
	dy3 += (1.0)*sqrt( tlB * y[2] )*W6  +  (-1.0)*sqrt( degB * y[3] )*W7 ;
	dy4 += (1.0)*sqrt( fc( y[1] , y[3] , y[5] , tid ) )*W8  +  (-1.0)*sqrt( degmC * y[4] )*W9 ;
	dy5 += (1.0)*sqrt( tlC * y[4] )*W10  +  (-1.0)*sqrt( degC * y[5] )*W11 ;


	y[0] += dy0;
	y[1] += dy1;
	y[2] += dy2;
	y[3] += dy3;
	y[4] += dy4;
	y[5] += dy5;


}
