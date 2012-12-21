/*
 * EulerSDESolver.cpp
 *
 *  Created on: Mar 15, 2011
 *      Author: student1 and student3
 */
#include "EulerSDESolver.hpp"
#include "ChildModel.hpp"
#include <iostream>
#include <vector>

using namespace std;

extern "C" {

  void maketime (vector<double>& vtimepoints, double* timepoints, int cntimepoints){
    for (int k = 0; k < cntimepoints; k++ ){	
      vtimepoints.push_back(timepoints[k]);
    }
  }

 int MainC(double* initialValues, double* parameters, int cbeta, double* timepoints, double dt, int cntimepoints, int NPARAMETERS, int NSPECIES, double* output){

    vector<double> vtimepoints;
    maketime(vtimepoints, timepoints, cntimepoints);

    ChildModel model(0);
    EulerSDESolver solver(initialValues, parameters, &model, cbeta, vtimepoints,"Eulertest.txt",dt);				
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


EulerSDESolver::EulerSDESolver(double* ainitialValues, double* aparameters,
		ChildModel* mmodel, int inbOfIterations, vector<double> vtimepoints,
		string sfilename, double ddt) :
	SDESolver(ainitialValues, aparameters, mmodel, inbOfIterations,
			vtimepoints, sfilename, ddt) {
	run();
}

EulerSDESolver::~EulerSDESolver() {

}

vector<double> EulerSDESolver::getNoise() {
	//We only need one noise per reaction for this algorithm
	vector<double> noise(model->NREACTIONS, 0);
	int j;
	for (j = 0; j < model->NREACTIONS; j++) {
		noise[j] = gsl_ran_gaussian(r, sqrt(dt));
	}
	return noise;
}

void EulerSDESolver::getEulerNoiseTerm(const ColumnVector& hazards,
		const vector<double>& noise) {
	//For each specie the term given by the euler approximation is the sum over the reactions of the diffusion matrix for this specie and this reaction multiplied by the noise of this reaction
	Matrix diffusionMatrix = getDiffusionMatrix(hazards);
	int k;
	int j;
	for (k = 0; k < model->NSPECIES; k++) {
		dyEuler[k] = 0;
		for (j = 0; j < model->NREACTIONS; j++) {
			dyEuler[k] += diffusionMatrix(k + 1, j + 1) * noise[j];
		}
	}
}

int EulerSDESolver::checkY(double* candidateDy) {
	ColumnVector candidateY(model->NSPECIES);
	int stop = 1;
	int k;
	//We create a vector containing the values of the concentrations if we kept the candidateDy
	for (k = 0; k < model->NSPECIES; k++) {
		candidateY(k + 1) = candidateDy[k] + concentrations[k];
	}
	//And we check if all the values are positive, if they're not we return 0 meaning that the candidateDy is not valid.
	if (candidateY.minimum() < 0) {
		stop = 0;
	}
	return stop;
}

void EulerSDESolver::changeConc() {
	//We first check if the candidate dy is valid and as long as it is not we keep on redoing the step
	while (checkY(dy) == 0) {
		step();
	}
	//Then we increment the concentrations by the dy
	for (int k = 0; k < model->NSPECIES; k++) {
		concentrations[k] = concentrations[k] + dy[k];
	}
}

void EulerSDESolver::step() {
	// 1- Get a random noise for each reaction
	vector<double> noise = getNoise();

	// 2- Compute the hazards for the current concentrations and parameters
	ColumnVector hazards = model->getHazards(concentrations, parameters);

	// 3- Get the ODE part of the derivatives as well as the Euler part
	dyODE = getODEterm(hazards);
	getEulerNoiseTerm(hazards, noise);

	// 4- the sum of the two parts give us the derivative value
	int k;
	for (k = 0; k < model->NSPECIES; k++) {
		dy[k] = dyODE[k] + dyEuler[k];
	}
}

