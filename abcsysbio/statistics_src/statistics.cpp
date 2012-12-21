//============================================================================
// Name        : statistics.cpp
// Author      : Aaron Sim
// Version     :
// Copyright   : Aaron Sim
// Description : Statistics library in C for abcsysbio
//============================================================================

//This file contains a multinomial sampling function, sampling functions and probability density functions for uniform, gaussian, and truncated normal functions in both univariate and multivariate cases. It uses GSL library as it a free and computationally efficient.

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

#include <sys/time.h>
#include <time.h>
#include <sys/types.h>

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//tell C++ compiler that the function is compiled in C style
extern "C"{

extern const gsl_rng_type *gsl_rng_default;
extern unsigned long int gsl_rng_default_seed;
gsl_rng * r;
const gsl_rng_type * T;

//set the random Generator function // used gettimeofday as clock_gettime has problems in compiling in Mac OS X
void setRandomGenerator() {
  const gsl_rng_type * T;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  //timespec ts;
  //clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  //gsl_rng_set(r, (int) ts.tv_nsec);
  
  struct timeval tv;
  gettimeofday(&tv, NULL);
  gsl_rng_set(r, (int) tv.tv_usec);
}



///MULTINOMIAL SAMPLING: "weight" is a list containing "nweight" number of weights from which one is randomly chosen.
int sample_multinomial (double* weight, int nweight){

	double weight_chooser =  gsl_rng_uniform(r);
	double totalweight = 0;
	unsigned int i;
	double a=0;
	for (i=0; i < nweight; ++i){
		 totalweight = totalweight + weight[i];
	}
	for (i=0; i < nweight; ++i){
		if(weight[i] != 0)  {
			a += weight[i]/totalweight;
		}
		if (weight_chooser < a){
			break;
		 }
	}
	return i;
}


///////////////SAMPLING UNIFORM DISTRIBUTION: "sample" is a random number generated from a distribution U(a,b)
int sampleUniform (double a, double b, double* sample){
	(*sample)= gsl_ran_flat(r, a, b);
//	cout<< "s:"<<sample << endl;
	return 0;
}


///////////// SAMPLING GAUSSIAN DISTRIBUTION: "sample" is a random number generated from a distribution N(mean,sd^2)

  int sampleGauss (double mean, double sd, double* sample){
    (*sample) = mean + gsl_ran_gaussian(r,sd);//we don't want it with mean 0
    return 0;
  }

///////////////SAMPLING TRUNCATED UNIFORM DISTRIBUTION:"sample" is a random number generated from a distribution U(a,b) truncated at "lower" and "upper".
int sampleTruncatedUniform (double a, double b, double* sample, double lower, double upper){
	//if (a=-0){
	//	a=0;
	//}
	if (a>lower && b<upper){
  		(*sample)= gsl_ran_flat(r, a, b);
	}
	else if(a<=lower && b<=upper){
		(*sample)= gsl_ran_flat(r, lower, b);
	}
 	else if (a<=lower && b>=upper){
		(*sample)= gsl_ran_flat(r, lower, upper);
	}
	else if (a>=lower && b>=upper){
  		(*sample)= gsl_ran_flat(r, a, upper);
	}else{
		cout << "oops!" << endl;
	}
	return 0;
}

///////////// SAMPLING TRUNCATED GAUSSIAN DISTRIBUTION: "sample" is a random number generated from a distribution N(mean,sd^2)truncated at "lower" and "upper".
int sampleTruncatedGauss(double mean, double sd, double* sample, double lower, double upper){
    //double x= mean + gsl_ran_gaussian(r,sd);//we don't want it with mean 0
    //if (j<x<k){
	//(*sample)=x
    	//return 0;
    //}
	//we don't want it with mean 0
	double g=-1.0;
	double u=0.0;
	double z=0.0;
	double lowerd=(lower-mean)/sd;
	double upperd=(upper-mean)/sd;
	//cout << "inside C" << endl;
	while(g<u){
		z=gsl_ran_flat(r,lowerd,upperd);
		//cout<<"z:"<<z<<endl;
		if (lowerd<=0 && 0<=upperd){
			g=exp(-(z*z)/2);
			//cout << "case1" << endl;
		}
		else if (upperd<0){
			g=exp(((upperd*upperd)-(z*z))/2);
			//cout << "case2" << endl;
		}
		else if (0<lowerd){
			g=exp(((lowerd*lowerd)-(z*z))/2);
			//cout << "case3" << endl;
		}
		u=gsl_ran_flat(r,0,1);
	}
	(*sample)=mean + sd*z;
	return 0;
}

///////////PDF OF UNIFORM DISTRIBUTION:"density" is the probability of point "x" being generated from a distribution U(a,b)
int pdfUniform( double x, double a, double b, double* density){

  (*density) =  gsl_ran_flat_pdf(x, a, b);
  return 0;
}


//////////PDF OF GAUSSIAN DISTRIBUTION:"density" is the probability of point "x" being generated from a distribution N(mean,sd^2)

int pdfGauss (double mean, double sd, double x, double* density){
	(*density)= gsl_ran_gaussian_pdf (x-mean, sd);
	return 0;
}

//////////PDF OF UNIFORM DISTRIBUTION + PRIORS:"density" is the probability of point "x" being generated from a distribution U(a,b) truncated at "lower" and "upper".

int pdfTruncatedUniform (double x, double a, double b, double lower, double upper, double* density){
	if (x>=lower && x<=upper){
		if (a>lower && b<upper){
  			(*density) =  gsl_ran_flat_pdf(x, a, b);
		}
		else if (a<=lower && b<=upper){
			(*density) =  gsl_ran_flat_pdf(x, lower, b);
		}
 		else if (a<=lower && b>=upper){
			(*density) =  gsl_ran_flat_pdf(x, lower, upper);
		}
		else if (a>=lower && b>=upper){
  			(*density) =  gsl_ran_flat_pdf(x, a, upper);
			//cout << "case4" << endl;
		}
 	}
	else {
		cout << "oops!" << endl;
		(*density)=0;
	}
	return 0;
}

/////////PDF OF GAUSSIAN DISTRIBUTION + PRIORS:"density" is the probability of point "x" being generated from a distribution N(mean,sd^2)truncated at "lower" and "upper".

int pdfTruncatedGauss (double mean, double sd, double x, double* density, double lower, double upper){
  if (lower<x<upper){
      if (upper == 10001.0) {
          (*density)= (1/sd)*(gsl_ran_gaussian_pdf ((x-mean)/sd, sd))/(1 -(gsl_cdf_gaussian_P((lower-mean)/sd, sd)));
          return 0;
      }
      else if (lower == -10001.0) {
          (*density)= (1/sd)*(gsl_ran_gaussian_pdf ((x-mean)/sd, sd))/(gsl_cdf_gaussian_P((upper-mean)/sd, sd));
          return 0;
      }
	(*density)= (1/sd)*(gsl_ran_gaussian_pdf ((x-mean)/sd, sd))/((gsl_cdf_gaussian_P((upper-mean)/sd, sd))-(gsl_cdf_gaussian_P((lower-mean)/sd, sd)));
	return 0;
 }
  else{
	(*density)=0;
	return 0;
 }
}

/////////////SAMPLING FOR MULTIVARIATE NORMAL:"sample" is a vector obtained by sampling from a multivariate distribution N(mean,A) where the size of the vectors/matrix is sizeA.

int samplingMultivariateGauss(double* mean, double* A, int sizeA, double* sample){

	int i, j;
	j=0;

	//create gsl vector mean and gsl matrix
	gsl_vector* gsl_mean=gsl_vector_alloc(sizeA);
	gsl_matrix* gsl_A= gsl_matrix_alloc(sizeA,sizeA);

	for(i=0; i<sizeA*sizeA; i++){
		if(i%sizeA==0 && i!=0){
			j++;
			}
		gsl_matrix_set(gsl_A,i%sizeA,j,A[i]);
	}
	for (i=0; i<sizeA; i++){
		gsl_vector_set(gsl_mean,i,mean[i]);
	}

	//if matrix is positive definite cholesky decomposition
	int error = gsl_linalg_cholesky_decomp(gsl_A);
	if( error == GSL_EDOM ){
		cout<<"matrix is not positive definite\n";
		exit(1);
		}

	int a;
	gsl_vector* outputVector=gsl_vector_alloc(sizeA);
	for (a=0; a<sizeA; a++){ //we want n independant standarn normal variates and saving them in outputVector
		gsl_vector_set( outputVector, a, gsl_ran_gaussian(r,1) );
	}

	//multiplication between lower triangle of cholesky matrix and outputVector
	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, gsl_A, outputVector);
	gsl_matrix_free(gsl_A);
	gsl_vector_add(outputVector, gsl_mean);
	for(i=0;i<sizeA;i++){
		sample[i]=gsl_vector_get (outputVector,i);
	}
	gsl_vector_free(gsl_mean);
	gsl_vector_free(outputVector);

	return 0;
}

///////////PDF FOR MULTIVARIATE NORMAL:"density" is the probability of vector"x" being generated from a multivariate distribution N(mean,A)

int pdfMultivariateGauss(double* mean, double* A, int sizeA, double* x,double* density){

	gsl_vector* gsl_mean=gsl_vector_alloc(sizeA);
	gsl_vector* gsl_x=gsl_vector_alloc(sizeA);
	gsl_matrix* gsl_A= gsl_matrix_alloc(sizeA,sizeA);
	gsl_matrix* gsl_Ainv= gsl_matrix_alloc(sizeA,sizeA);
	gsl_permutation *p = gsl_permutation_alloc(sizeA);

	int i,j;
	j=0;
	for(i=0; i<sizeA*sizeA; i++){

		if(i%sizeA==0 && i!=0){
			j++;
			}
		gsl_matrix_set(gsl_A,i%sizeA,j,A[i]);
	}

	for (i=0; i<sizeA; i++){
		gsl_vector_set(gsl_mean,i,mean[i]);
		gsl_vector_set(gsl_x,i,x[i]);
	}

	int a;
	gsl_linalg_LU_decomp( gsl_A, p, &a );// LU decompostion
	gsl_linalg_LU_invert( gsl_A, p, gsl_Ainv );//get the inverse
	double det = gsl_linalg_LU_det( gsl_A, a);//get the determinant
	gsl_matrix_free(gsl_A);
	gsl_permutation_free(p);

	gsl_vector* xm = gsl_vector_alloc(sizeA);//define a new vector
	gsl_vector_memcpy( xm, gsl_x);//copy elements of vector x into xm
	gsl_vector_sub( xm, gsl_mean );//subtract the elements in mean from the elements of vector x
	gsl_vector* ym = gsl_vector_alloc(sizeA);
	gsl_blas_dsymv(CblasUpper,1.0,gsl_Ainv,xm,0.0,ym);//matrix vector multiplication : Y = A^-1 * X

	gsl_matrix_free(gsl_Ainv);
	double kern; //you could also just directly use density
	gsl_blas_ddot( xm, ym, &kern);//scalar product x^T y for the vectors xm and ym, returning the result in &ay : X^t * Y
	(*density) = exp(-0.5*kern)/sqrt( pow((2*M_PI),sizeA)*det );


	gsl_vector_free(gsl_mean);
	gsl_vector_free(gsl_x);
	gsl_vector_free(xm);
	gsl_vector_free(ym);

	return 0;
}

///////////SAMPLING FOR MULTIVARIATE TRUNCATED NORMAL:"chosen_sample" is a vector obtained by sampling from a multivariate distribution N(mean,A) where the size of the vectors/matrix is sizeA and which has been truncated at vectors "lower" and "upper".

int sampleMultivariateTGauss(double* mean, double* A, int sizeA, double* sample, double* lower, double* upper){

	int i,j;
	j=0;
	while (j!=sizeA){
		j=0;
		//"sample" must be between "lower" and "upper" in order to be accepted
		int x=samplingMultivariateGauss(mean, A, sizeA, sample);
		for (i=0; i<sizeA; i++){
			if (lower[i]<sample[i] && sample[i]<upper[i]){
				//chosen_sample[i]=sample[i];
				j++;
			}
			else {
				j=0;
			}
		}
	}
	return 0;
}

///////////MCMC SAMPLING FOR MULTIVARIATE TRUNCATED NORMAL

int MCMCsampleMultivariateTGauss(double* A, double* mean, double* lower, double* upper, int sizeA, double* result, double* sample, double* upperd, double* lowerd){

	gsl_vector* gsl_mean=gsl_vector_alloc(sizeA);
	gsl_vector* gsl_sample=gsl_vector_alloc(sizeA);
	gsl_vector* gsl_mean_i=gsl_vector_alloc(sizeA-1);
	gsl_vector* gsl_sample_i=gsl_vector_alloc(sizeA-1);
	gsl_vector* outputVector1=gsl_vector_alloc(sizeA-1);
	gsl_vector* vector_row=gsl_vector_alloc(sizeA);
	gsl_vector* outputVector2=gsl_vector_alloc(sizeA-1);
	gsl_vector* vector_i=gsl_vector_alloc(sizeA-1);

	gsl_permutation* p = gsl_permutation_alloc(sizeA-1);
	gsl_matrix* gsl_A= gsl_matrix_alloc(sizeA,sizeA);
	gsl_matrix* gsl_A_i= gsl_matrix_alloc(sizeA-1,sizeA-1);
	gsl_matrix* gsl_Ainv_i= gsl_matrix_alloc(sizeA-1,sizeA-1);

	int i,j,k;
	j=0;
	//set the matrix
	for(i=0; i<sizeA*sizeA; i++){
		if(i%sizeA==0 && i!=0){
			j++;
			}
		gsl_matrix_set(gsl_A,i%sizeA,j,A[i]);
	}

	//sample from above multiv
	sampleMultivariateTGauss(mean, A, sizeA, sample, lower, upper);

	for (i=0; i<sizeA; i++){
		gsl_vector_set(gsl_mean,i,mean[i]);//set the mean vector
		gsl_vector_set(gsl_sample,i,sample[i]);
	}
	for (i=0; i<sizeA; i++){
		gsl_matrix_swap_columns (gsl_A, i, sizeA-1);
		gsl_matrix_swap_rows (gsl_A, i, sizeA-1);

		gsl_vector_swap_elements (gsl_mean,i, sizeA-1);
		gsl_vector_swap_elements (gsl_sample,i, sizeA-1);

		for (j=0; j<(sizeA-1); j++){
			gsl_vector_set (gsl_mean_i,j, gsl_vector_get(gsl_mean, j));//set mean-i
			gsl_vector_set (gsl_sample_i,j, gsl_vector_get(gsl_sample, j));//set sample-i
			for (k=0; k<(sizeA-1); k++){
				gsl_matrix_set(gsl_A_i,j,k,gsl_matrix_get (gsl_A, j,k));//set matrix -i row&column
			}
		}

		//set matrix and vectors back to normal order
		gsl_matrix_swap_columns (gsl_A, sizeA-1, i);
		gsl_matrix_swap_rows (gsl_A, sizeA-1, i);
		gsl_vector_swap_elements (gsl_mean, sizeA-1, i);
		gsl_vector_swap_elements (gsl_sample, sizeA-1, i);

		//select i-th column of matrix and copy it to vector
		//cout<<"i="<<i<<endl;
		gsl_matrix_get_col(vector_row,gsl_A,i);
		//cout<<"here2"<<endl;
		gsl_vector_swap_elements (vector_row,i, sizeA-1);
		for (j=0;j<(sizeA-1);j++){
			gsl_vector_set(vector_i,j,gsl_vector_get(vector_row,j));//vector -i
		}


		//gsl_matrix_transpose_memcpy (vector_i_trans, vector_i);//get transpose of vector-i: iy is equivalent to vector-i
		//gsl_linalg_LU_decomp( gsl_A, p, &a );// LU decompostion
		//gsl_linalg_LU_invert( gsl_A, p, gsl_Ainv );//get the inverse
		int a;
		gsl_linalg_LU_decomp( gsl_A_i, p, &a );// LU decompostion
		gsl_linalg_LU_invert( gsl_A_i, p, gsl_Ainv_i );//get the inverse


		//gsl_matrix_mul_elements(gsl_Ainv_i,vector_i_trans);//multiply 2 matrices
		gsl_vector_sub(gsl_sample_i,gsl_mean_i);//subtract vectors
		gsl_blas_dsymv(CblasUpper, 1.0,gsl_Ainv_i,gsl_sample_i,0.0,outputVector1);//matrix vector multiplication
		gsl_blas_dsymv(CblasUpper, 1.0,gsl_Ainv_i,vector_i,0.0,outputVector2);//matrix vector multiplication
		int b=gsl_vector_mul(outputVector1, vector_i);//multiply vectors
		int c=gsl_vector_mul(outputVector2, vector_i);
		double diag_cov=gsl_matrix_get(gsl_A,i,i);
		//cout<<"diag_cov"<<diag_cov<<endl;
		//cout<<"b"<<b<<endl;
		//cout<<"c"<<c<<endl;
		double mu=mean[i]+b;
		//cout<<"mu"<<mu<<endl;
		double sd=diag_cov-c;
		//cout<<"sd"<<sd<<endl;
		double g=-1.0;
		double u=0.0;
		double z=0.0;

		lowerd[i]=lower[i]-mu;
	 	upperd[i]=upper[i]-mu;
		while(g<u){
			z=gsl_ran_flat(r,lowerd[i],upperd[i]);
				if (lowerd[i]<0 && 0<upperd[i]){
					g=exp(-z*z/2);
				}
				if (upperd[i]<0){
					g=exp(((upperd[i]*upperd[i])-z*z)/2);
				}
				if (0<lowerd[i]){
					g=exp(((lowerd[i]*lowerd[i])-upperd[i]*upperd[i])/2);
				}
				u=gsl_ran_flat(r,0,1);
		}
		//cout<<"z"<<z<<endl;
		result[i]=mu +z*sd;
		//cout<<"result"<< " "<<i<< " " << result[i]<<endl;

		//here is should fill sample with z.

	}
	gsl_vector_free(gsl_mean);
	gsl_vector_free(gsl_sample);
	gsl_vector_free(gsl_mean_i);
	gsl_vector_free(gsl_sample_i);
	gsl_vector_free(vector_row);
	gsl_vector_free(vector_i);
	gsl_vector_free(outputVector1);
	gsl_vector_free(outputVector2);

	gsl_matrix_free(gsl_A_i);
	gsl_matrix_free(gsl_Ainv_i);
	gsl_matrix_free(gsl_A);
	gsl_permutation_free(p);
	return 0;
}


///////////PDF FOR MULTIVARIATE TRUNCATED NORMAL - using Genz (1993) algorithm for distribution function truncated MVN

int pdfTMultivariateGauss(double* mean, double* A, int sizeA, double* x, double* lower, double* upper,
							double* density, double errorlimit=0.005, double MCconfidence=2.5, double NMax=100) {
	// Assume for now that the truncation limits are finite
	// Arguments : (mean, covariance, no. of variables, quantile, vector of lower bounds, vector of upper bounds,
	// 				output pdf, prob error, Monte Carlo confidence factor, Limit of number of iterations)

	// Allocate the memory
	gsl_vector* gsl_mean=gsl_vector_alloc(sizeA);
	gsl_vector* gsl_x=gsl_vector_alloc(sizeA);
	gsl_vector* gsl_lower=gsl_vector_alloc(sizeA);
	gsl_vector* gsl_upper=gsl_vector_alloc(sizeA);
	gsl_matrix* gsl_A= gsl_matrix_alloc(sizeA,sizeA);
	gsl_matrix* gsl_Ainv= gsl_matrix_alloc(sizeA,sizeA);

	// Initialise the covariance matrix and simultaneously make a duplicate
	int i,j;
	j=0;
	double A_temp[sizeA*sizeA];
	for(i=0; i<sizeA*sizeA; i++){
		A_temp[i] = A[i];
		if(i%sizeA==0 && i!=0){
			j++;
			}
		gsl_matrix_set(gsl_A,i%sizeA,j,A[i]); // Initialise covariance matrix
	}

	// Initialise the mean, lower and upper vectors
	for (i=0; i<sizeA; i++){
		gsl_vector_set(gsl_mean,i,mean[i]);
		gsl_vector_set(gsl_x,i,x[i]);
		gsl_vector_set(gsl_lower,i,lower[i]);
		gsl_vector_set(gsl_upper,i,upper[i]);
	}

	pdfMultivariateGauss(mean, A_temp, sizeA, x, density); // Calculate the pdf of untruncated version first - then divide by intsum below

	// Compute the lower triangular Cholesky factor C for A
	// ***** NEED TO CHECK IF COVARIANCE IS POSITIVE DEFINITE - checkinputArguments.py *****
	int error = gsl_linalg_cholesky_decomp(gsl_A);
	if( error == GSL_EDOM ){
		cout<<"matrix is not positive definite\n";
		exit(1);
		}

	// Shift lower and upper limits and the queried x to the mean
	gsl_vector_sub(gsl_lower, gsl_mean);
	gsl_vector_sub(gsl_upper, gsl_mean);
	gsl_vector_sub(gsl_x, gsl_mean);

	// Initialise variables
	double intsum=0;
	double varsum=0;
	int N = 0;
	double d[sizeA];
	double e[sizeA];
	double f[sizeA];
	double delta_temp;
	double errorterm=errorlimit+1; // Initalise errorterm to be greater than errorlimit
	double w[sizeA-1];
	if (gsl_vector_get(gsl_lower, 0)==-10001) {
		d[0] = 0; // limit if lower limit is not truncated; -10000001 is a placeholder for -inf
	} else {
		d[0] = gsl_cdf_gaussian_P((gsl_vector_get(gsl_lower, 0))/gsl_matrix_get(gsl_A,0,0), 1);
	}
	if (gsl_vector_get(gsl_upper, 0)==10001) {
		e[0] = 1; // limit if upper limit is not truncated; -10000001 is a placeholder for +inf
	} else {
		e[0] = gsl_cdf_gaussian_P((gsl_vector_get(gsl_upper, 0))/gsl_matrix_get(gsl_A,0,0), 1);
	}
	f[0] = e[0]-d[0];


	//setRandomGenerator();
	while(errorterm >= errorlimit && N != NMax) {
		// Generate uniform random variables w1, w2,..., w_{m-1}
		for (i=0; i<(sizeA-1); i++){
			w[i] = gsl_ran_flat(r, 0, 1);
		}

		// Iterate the MC loop
		double y[sizeA-1];
		for (i=1; i<sizeA; i++) {
			y[i-1] = gsl_cdf_gaussian_Pinv(d[i-1]+w[i-1]*(e[i-1]-d[i-1]), 1);
			double temp_sum = 0;
			for (j=0; j<i; j++){
				temp_sum += gsl_matrix_get(gsl_A, i, j)*y[j];
			}
			if (gsl_vector_get(gsl_lower, i)==-10001) {
				d[i] = 0; // limit if lower limit is not truncated; -1000 is a placeholder for -inf
			} else {
				d[i] = gsl_cdf_gaussian_P((gsl_vector_get(gsl_lower, i)-temp_sum)/gsl_matrix_get(gsl_A, i, i), 1);
			}
			if (gsl_vector_get(gsl_upper, i)==10001) {
				e[i] = 1; // limit if upper limit is not truncated; -1000 is a placeholder for +inf
			} else {
				e[i] = gsl_cdf_gaussian_P((gsl_vector_get(gsl_upper, i)-temp_sum)/gsl_matrix_get(gsl_A, i, i), 1);
			}
			f[i] = (e[i] - d[i])*f[i-1];
		}
		N = N + 1;
		delta_temp = (f[sizeA-1]-intsum)/N;
		intsum = intsum + delta_temp;
		varsum = ((N-2)*varsum)/N + delta_temp*delta_temp;
		errorterm = MCconfidence*sqrt(varsum);
	}

	// Calculate pdf depending on whether it lies within the non-truncated region

	gsl_vector_sub(gsl_lower, gsl_x);
	gsl_vector_sub(gsl_upper, gsl_x);
	if ((gsl_vector_max(gsl_lower)>0) || (gsl_vector_min(gsl_upper)<0)) {
		*density = 0;
	} else {
		*density = *density/intsum;
	}

	// Release the memories
	gsl_vector_free(gsl_mean);
	gsl_vector_free(gsl_x);
	gsl_vector_free(gsl_lower);
	gsl_vector_free(gsl_upper);
	gsl_matrix_free(gsl_A);
	gsl_matrix_free(gsl_Ainv);

	return 0;
}
/*
int main() {
	return 0;
}
*/
//end extern C
}









