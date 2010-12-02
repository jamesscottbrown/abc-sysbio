from math import *
from numpy import *
from abcsysbio.relations import *

def modelfunction((m1,p1,m2,p2,m3,p3,),t,parameter=(1,2,5,1000,)):

	alpha0=parameter[0]
	n=parameter[1]
	beta=parameter[2]
	alpha=parameter[3]
	
	d_m1=-m1+alpha/(1+p3**n)+alpha0
	d_p1=-beta*(p1-m1)
	d_m2=-m2+alpha/(1+p1**n)+alpha0
	d_p2=-beta*(p2-m2)
	d_m3=-m3+alpha/(1+p2**n)+alpha0
	d_p3=-beta*(p3-m3)


	return(d_m1,d_p1,d_m2,d_p2,d_m3,d_p3, )


def rules((m1,p1,m2,p2,m3,p3,),parameter,t):



	return((m1,p1,m2,p2,m3,p3, ),parameter)

