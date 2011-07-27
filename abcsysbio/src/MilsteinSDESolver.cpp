/*
 * MilsteinSDESolver.cpp
 *
 *  Created on: Mar 14, 2011
 *      Author: student1
 */
#include "MilsteinSDESolver.hpp"
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
    MilsteinSDESolver solver(initialValues, parameters, &model, cbeta, vtimepoints,"Milsteintest.txt",dt);				
    
// Loop over output from the solver
    int count=0;
    for (int i=0; i<cbeta; i++){
      	for(int k=1; k<NSPECIES+1; k++){
		for(int j=0; j<cntimepoints; j++){
			output[count] = solver.output[i][j][k];
	  		//cout << count << "   " << solver.output[i][j][k] << "   "<< output[count] << endl;
	    		//cout << i << "\t" << j << "\t" << k << endl;
			count++;
	    }
      }
    }
    return 0;  
  }
  
}

MilsteinSDESolver::MilsteinSDESolver(double* ainitialValues,
		double* aparameters, ChildModel* mmodel, int inbOfIterations,
		vector<double> vtimepoints, string sfilename, double ddt) :
	SDESolver(ainitialValues, aparameters, mmodel, inbOfIterations,
			vtimepoints, sfilename, ddt) {

	dyMilstein = new double[model->NSPECIES];
	adaptedYj = new double[model->NSPECIES];
	lastStepConcentrations = new double[model->NSPECIES];

	run();
}

MilsteinSDESolver::~MilsteinSDESolver() {
	delete[] dyMilstein;
	delete[] adaptedYj;
	delete[] lastStepConcentrations;
}

void MilsteinSDESolver::getEulerNoiseTerm(const ColumnVector& hazards,
		const vector<vector<double> >& noise) {
	//For each specie the term given by the euler approximation is the sum over the reactions of the diffusion matrix
	//for this specie and this reaction multiplied by the noise of this reaction
	Matrix diffusionMat = getDiffusionMatrix(hazards);

	int k;
	int j;
	for (k = 0; k < model->NSPECIES; k++) {
		dyEuler[k] = 0;
		for (j = 0; j < model->NREACTIONS; j++) {
			dyEuler[k] += diffusionMat(k + 1, j + 1) * noise[j][0];
		}
	}
}

double MilsteinSDESolver::getDoubleStochasticIntegral(int j1, int j2,
		const vector<vector<double> >& noise) {

	double I = 0.0;
	if (j1 != j2) {
		I = (0.5 * noise[j1][0] * noise[j2][0] + sqrt(1.0 / 12.0 - 1.0 / (2
				* pow(M_PI, 2.0))) * (noise[j1][1] * sqrt(dt) * noise[j2][0]
				- noise[j2][1] * sqrt(dt) * noise[j1][0])) + dt / 2
				* (noise[j1][2] * (sqrt(2) * noise[j1][0] / sqrt(dt)
						+ noise[j2][3]) - noise[j2][2] * (sqrt(2)
						* noise[j2][0] / sqrt(dt) + noise[j2][3]));
	} else {
		I = 1.0 / 2.0 * (pow(noise[j1][0], 2) - dt);
	}

	return I;
}

void MilsteinSDESolver::getMilsteinNoiseTerm(const Matrix& diffusionMatrix,
		const ColumnVector& driftVector, const double yColumn[], const vector<
				vector<double> >& noise) {

	int k;
	int j1;
	int j2;
	int k1;
	double I;
	ColumnVector adaptedHazards(model->NREACTIONS);
	Matrix adaptedDiffusionMatrix(model->NSPECIES, model->NREACTIONS);
	// For each specie the Milstein part of the derivative is the double sum over all the reactions of the double
	// stochastic integral times the difference of the diffusionMatrix for adapted concentrations minus the diffusion matrix for the current concentrations
	for (k = 0; k < model->NSPECIES; k++) {
		dyMilstein[k] = 0;

		for (j1 = 0; j1 < model->NREACTIONS; j1++) {
			//We first compute the adapted concentrations for the reaction j1
			for (k1 = 0; k1 < model->NSPECIES; k1++) {
				adaptedYj[k1] = yColumn[k1] + driftVector(k1 + 1) * dt
						+ diffusionMatrix(k1 + 1, j1 + 1) * sqrt(dt);
			}
			//To enable us to compute the adaptedHazards and thus adapted diffusion matrix
							//NOTE: If the previous step choose the "wrong" noise in the stochastic integral the adapted hazards
							//can get negative, and thus dyMilstein won't be defined
							//In such a case computing the step again won't help because the previous calculus is not stochastic
							//so we need to go back to 2 steps before. Hence the change in the checkY and changeConc compared
							//to the other algorithms
			adaptedHazards = model->getHazards(adaptedYj, parameters);
			adaptedDiffusionMatrix = getDiffusionMatrix(adaptedHazards);

			for (j2 = 0; j2 < model->NREACTIONS; j2++) {
				//Then we get the value of the stochastic integral for each j2 and compute the value of our derivative
							//NOTE : the stochasticity of the algorithm is completely contained in this stochastic integral,
				I = getDoubleStochasticIntegral(j1, j2, noise);

				dyMilstein[k] += 1.0 / sqrt(dt) * (adaptedDiffusionMatrix(
						k + 1, j2 + 1) - diffusionMatrix(k + 1, j2 + 1)) * I;
			}
		}
	}
}

vector<vector<double> > MilsteinSDESolver::getNoise() {
	//In this algorithm we need 4 noise terms
	vector<vector<double> > noise(model->NREACTIONS, vector<double> (4, 0));
	int j;
	for (j = 0; j < model->NREACTIONS; j++) {
		noise[j][0] = gsl_ran_gaussian(r, sqrt(dt));
		noise[j][1] = gsl_ran_gaussian(r, 1);
		noise[j][2] = gsl_ran_gaussian(r, 1);
		noise[j][3] = gsl_ran_gaussian(r, 1);
	}
	return noise;
}

int MilsteinSDESolver::checkY(double* candidateDy) {
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
	} else {
		for (k = 0; k < model->NSPECIES; k++) {
			//We also check if the values are defined, if not we change the concentrations back to their
			//previous value (2 steps ago) and return 0 to say that the candidate is not valid
			if (dy[k] != dy[k]) {
				concentrations[k] = lastStepConcentrations[k];
				stop = 0;
			}
		}
	}
	return stop;
}

void MilsteinSDESolver::changeConc() {
	//We first check if the candidate dy is valid and as long as it is not we keep on redoing the step
	while (checkY(dy) == 0) {
		step();
	}
	//Then we store the value of the concentrations to lastStepConcentrations and increment the concentrations by the dy
	for (int k = 0; k < model->NSPECIES; k++) {
		lastStepConcentrations[k] = concentrations[k];
		concentrations[k] = concentrations[k] + dy[k];
	}
}

void MilsteinSDESolver::step() {
	vector<vector<double> > noise = getNoise();

	ColumnVector hazards = model->getHazards(concentrations, parameters);

	Matrix diffusionMatrix = getDiffusionMatrix(hazards);

	ColumnVector driftVector = getDriftVector(hazards);

	dyODE = getODEterm(hazards);
	getEulerNoiseTerm(hazards, noise);
	getMilsteinNoiseTerm(diffusionMatrix, driftVector, concentrations, noise);

	int k;
	for (k = 0; k < model->NSPECIES; k++) {
		dy[k] = dyODE[k];
		dy[k] += dyEuler[k];
		dy[k] += dyMilstein[k];
	}
}
