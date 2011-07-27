/*
 * HeunSDESolver.cpp
 *
 *  Created on: Mar 14, 2011
 *      Author: student1
 */
#include "HeunSDESolver.hpp"
#include "ChildModel.hpp"
#include <iostream>
#include <vector>

using namespace std;

extern "C" {

  void maketime (vector<double>& vtimepoints, double* timepoints, int cntimepoints){
    for (unsigned int k = 0; k < cntimepoints; k++ ){	
      vtimepoints.push_back(timepoints[k]);
    }
  }

 int MainC(double* initialValues, double* parameters, int cbeta, double* timepoints, double dt, int cntimepoints, int NPARAMETERS, int NSPECIES, double* output){
    vector<double> vtimepoints;
    maketime(vtimepoints, timepoints, cntimepoints);

    ChildModel model(0);
    HeunSDESolver solver(initialValues, parameters, &model, cbeta, vtimepoints,"Heuntest.txt",dt);				
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
HeunSDESolver::HeunSDESolver(double* ainitialValues, double* aparameters,
		ChildModel* mmodel, int inbOfIterations, vector<double> vtimepoints,
		string sfilename, double ddt) :
	SDESolver(ainitialValues, aparameters, mmodel, inbOfIterations,
			vtimepoints, sfilename, ddt) {

	intermediate_y = new double[model->NSPECIES];
	intermediate_dy = new double[model->NSPECIES];
	intermediate_dyODE = new double[model->NSPECIES];
	intermediate_dyEulerNoise = new double[model->NSPECIES];

	run();
}

HeunSDESolver::~HeunSDESolver() {
	delete[] intermediate_y;
	delete[] intermediate_dy;
	delete[] intermediate_dyODE;
	delete[] intermediate_dyEulerNoise;
}

double* HeunSDESolver::getEulerNoiseTerm(const ColumnVector& hazards,
		const vector<double>& noise) {
	//For each specie the term given by the euler approximation is the sum over the reactions of the diffusion matrix for this specie and this reaction multiplied by the noise of this reaction
	double* local_dyEuler = new double[model->NSPECIES];

	Matrix diffusionMatrix = getDiffusionMatrix(hazards);
	int k;
	int j;
	for (k = 0; k < model->NSPECIES; k++) {
		local_dyEuler[k] = 0;
		for (j = 0; j < model->NREACTIONS; j++) {
			local_dyEuler[k] += diffusionMatrix(k + 1, j + 1) * noise[j];
		}
	}

	return local_dyEuler;
}

vector<double> HeunSDESolver::getNoise() {
	//We only need one noise per reaction for this algorithm
	vector<double> noise(model->NREACTIONS, 0);
	int j;
	for (j = 0; j < model->NREACTIONS; j++) {
		noise[j] = gsl_ran_gaussian(r, sqrt(dt));
	}
	return noise;
}

int HeunSDESolver::checkY(double* candidateDy) {

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

void HeunSDESolver::changeConc() {
	//We first check if the candidate dy is valid and as long as it is not we keep on redoing the step
	while (checkY(dy) == 0) {
		step();
	}
	//Then we increment the concentrations by the dy
	for (int k = 0; k < model->NSPECIES; k++) {
		concentrations[k] = concentrations[k] + dy[k];
	}
}

void HeunSDESolver::step() {
	/**
	 * 1- Compute the Euler approximation for the current concentrations
	 *
	 * 2- Use the concentrations obtained (called intermediate_y) to compute the Euler approximation from this new point CAREFUL : using the same noise!!
	 *
	 * 3- Take the mean of the two sets of derivatives obtained to have our final set of derivatives
	 */
	vector<double> noise(model->NREACTIONS, 0);

	ColumnVector hazards = model->getHazards(concentrations, parameters);
	//We get the ODE term for the current concentrations
	intermediate_dyODE = getODEterm(hazards);
	//We initialise intermediate_dy to be sure its sum with concentrations will give negative values so the we're sure to enter the while loop and do the operations in it at least once
	int k;
	for (k = 0; k < model->NSPECIES; k++) {
		intermediate_dy[k] = -concentrations[k] - 1;
	}

	//As long as the intermediate_dy is not such that its sum with concentrations only gives positive values we keep doing the following operations
	while (checkY(intermediate_dy) == 0) {
		//We get randoms numbers for the noise
		noise = getNoise();
		//We get the Euler noise term for the current concentrations
		intermediate_dyEulerNoise = getEulerNoiseTerm(hazards, noise);
		//The derivative value at this point (which is just an intermediate value) is the sum of ODe and Euler
		for (k = 0; k < model->NSPECIES; k++) {
			intermediate_dy[k] = intermediate_dyODE[k]
					+ intermediate_dyEulerNoise[k];
		}

	}
	//And the intermediate concentrations are the sum of the intermediate derivatives and the concentrations
	for (k = 0; k < model->NSPECIES; k++) {
		intermediate_y[k] = intermediate_dy[k] + concentrations[k];
	}
	//We compute the hazards for these intermediate concentrations
	ColumnVector intermediateHazards = model->getHazards(intermediate_y,
			parameters);

	//The final derivatives are the mean of the derivatives for the current concentrations and the derivatives for the intermediate concentrations using in both the same noise
	for (k = 0; k < model->NSPECIES; k++) {
		dyODE[k] = 0.5 * (getODEterm(hazards)[k] + getODEterm(
				intermediateHazards)[k]);
		dyEuler[k] = 0.5 * (getEulerNoiseTerm(hazards, noise)[k]
				+ getEulerNoiseTerm(intermediateHazards, noise)[k]);
		dy[k] = dyODE[k] + dyEuler[k];
	}
}
