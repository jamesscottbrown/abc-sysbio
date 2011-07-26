/*
 * ODESolver.cpp
 *
 *  Created on: Mar 10, 2011
 *      Author: student1
 */
#include "ODESolver.hpp"
#include "ChildModel.hpp"
#include <iostream>
#include <vector>
using namespace std;

/**
 * EXPLAINATION FROM GSL DOCUMENTATION :
 *
 * The routines solve the general n-dimensional first-order system,

     dy_i(t)/dt = f_i(t, y_1(t), ..., y_n(t))

for i = 1, \dots, n. The stepping functions rely on the vector of derivatives f_i and the Jacobian matrix, J_{ij} = df_i(t,y(t)) / dy_j.
A system of equations is defined using the gsl_odeiv_system datatype.

— Data Type: gsl_odeiv_system

    This data type defines a general ODE system with arbitrary parameters.

    int (* function) (double t, const double y[], double dydt[], void * params)
        This function should store the vector elements f_i(t,y,params) in the array dydt, for arguments (t,y) and parameters params. The function should return GSL_SUCCESS if the calculation was completed successfully. Any other return value indicates an error.

    int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params);
        This function should store the vector of derivative elements df_i(t,y,params)/dt in the array dfdt and the Jacobian matrix J_{ij} in the array dfdy, regarded as a row-ordered matrix J(i,j) = dfdy[i * dimension + j] where dimension is the dimension of the system. The function should return GSL_SUCCESS if the calculation was completed successfully. Any other return value indicates an error.

        Some of the simpler solver algorithms do not make use of the Jacobian matrix, so it is not always strictly necessary to provide it (the jacobian element of the struct can be replaced by a null pointer for those algorithms). However, it is useful to provide the Jacobian to allow the solver algorithms to be interchanged—the best algorithms make use of the Jacobian.

    size_t dimension;
        This is the dimension of the system of equations.

    void * params
        This is a pointer to the arbitrary parameters of the system.
 */


extern "C" {

  void maketime (vector<double>& vtimepoints, double* timepoints, int cntimepoints){
    for (int k = 0; k < cntimepoints; k++ ){	
      vtimepoints.push_back(timepoints[k]);
    }
  }

 int MainC(double* initialValues, double* parameters, int cbeta, double* timepoints, int cntimepoints, int NPARAMETERS, int NSPECIES, double initstep, double absoluteError, double relativeError, double* output){
   
    vector<double> vtimepoints;
    maketime(vtimepoints, timepoints, cntimepoints);

    ChildModel model(0);
    ODESolver solver(initialValues, parameters, &model, cbeta, vtimepoints,"ODEtest.txt", initstep, absoluteError, relativeError);				
    // Loop over output from the solver
    int count=0;
    for (int i=0; i<cbeta; i++){
      	for(int k=1; k<NSPECIES+1; k++){
		for(int j=0; j<cntimepoints; j++){
			output[count] = solver.output[i][j][k];
	  		
			count++;
	    	}
      	}
    }
    return 0;  
  }
  
}

int func(double t, const double y[], double f[], void *p) {
	ODESolver* solver = (ODESolver*) p;
	ChildModel* model = solver->model;
	double* v = (*solver).reactions(y);
	for (int j = 0; j < (*model).NSPECIES; j++) {
		f[j] = v[j];
	}
	delete [] v;
	return GSL_SUCCESS;
}

int jac(double t, const double concentrations[], double *dfdy, double dfdt[], void *p) {
	return GSL_SUCCESS;
}

ODESolver::ODESolver(double* ainitialValues, double* aparameters, ChildModel* mmodel,
		int inbOfIterations, vector<double> vtimepoints, string sfilename,
		double dinitialStepSize, double dabsoluteError,
		double drelativeError) :
	Solver(ainitialValues, aparameters, mmodel, inbOfIterations, vtimepoints, sfilename) {
	initialStepSize = dinitialStepSize;
	absoluteError = dabsoluteError;
	relativeError = drelativeError;
	run();
}

ODESolver::~ODESolver() {

}

void ODESolver::solve() {
	setInitialValues();
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;

	//We first define the co;ponents of the system :

	//stepping function which advance a solution from time t to t+h for a fixed step-size h and estimate the resulting local error
	gsl_odeiv_step * s = gsl_odeiv_step_alloc(T, (*model).NSPECIES);

	//control function examines the proposed change to the solution produced by a stepping
	//function and attempts to determine the optimal step-size for a user-specified level of error
	gsl_odeiv_control * c = gsl_odeiv_control_y_new(absoluteError, relativeError);

	//the evolution function which combines the results of a stepping function and control
	//function to reliably advance the solution forward over an interval (t_0, t_1).
	//If the control function signals that the step-size should be decreased the evolution function
	//backs out of the current step and tries the proposed smaller step-size. This process is continued until an acceptable step-size is found.
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc((*model).NSPECIES);

	//the system itself
	gsl_odeiv_system sys = { func, jac, (*model).NSPECIES, this };

	double t = 0.0;
	

	for (counter = 0; counter < timepoints.size(); counter++) {
		double ti = timepoints[counter];
		(*model).applyRulesAndEvents(concentrations, parameters, ti);
		while (t < ti) {
			//This makes the system evolve
			int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, ti, &initialStepSize,
					concentrations);
			if (status != GSL_SUCCESS)
				break;
		}
		storePoints(ti);
	}

	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);
}

double* ODESolver::reactions(const double local_concentrations[]) {

	double* v = new double[(*model).NSPECIES];

	ColumnVector hazards = (*model).getHazards(local_concentrations, parameters);

	int j;
	int k;
	for (k = 0; k < (*model).NSPECIES; k++) {
		v[k]=0;
		for (j = 0; j < (*model).NREACTIONS; j++) {
			v[k] += (*model->pstoichiometricMatrix)(k + 1, j + 1)
					* hazards(j + 1);
		}
	}
	return v;
}
