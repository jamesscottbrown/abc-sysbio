/*
 * GillespieSolver.cpp
 *
 *  Created on: Mar 15, 2011
 *      Author: student1
 */
#include "GillespieSolver.hpp"
#include "ChildModel.hpp"
#include <iostream>
#include <vector>
using namespace std;

extern "C" {

  void maketime (vector<double>& vtimepoints, double* timepoints, int cntimepoints){
    for ( int k = 0; k < cntimepoints; k++ ){	
      vtimepoints.push_back(timepoints[k]);
    }
  }

 void MainC(double* initialValues, double* parameters, int cbeta, double* timepoints, int cntimepoints, int NPARAMETERS_ARG, int NSPECIES_ARG, double* output){
    int NPARAMETERS_INTERNAL=NPARAMETERS_ARG;
    int NSPECIES_INTERNAL=NSPECIES_ARG;
    vector<double> vtimepoints;
    maketime(vtimepoints, timepoints, cntimepoints);   

    ChildModel model(0);
  
    GillespieSolver solver(initialValues, parameters, &model, cbeta, vtimepoints,"Gillespietest.txt");				
   
    int count=0;
    for (int i=0; i<cbeta; i++){
      for(int k=1; k<NSPECIES_INTERNAL+1; k++){ // skip the time values
	for(int j=0; j<cntimepoints; j++){
	  output[count] = (double)solver.output[i][j][k];
	  //cout << i << "," << j << "," << k << "," << count << "\t" << output[count] << endl;
	  count++;
	}
      }
    }
 }
  

}

GillespieSolver::GillespieSolver(double* ainitialValues, double* aparameters,
		ChildModel* mmodel, int inbOfIterations, vector<double> vtimepoints, string sfilename) :
	Solver(ainitialValues, aparameters, mmodel, inbOfIterations, vtimepoints, sfilename) {
 
	setRandomGenerator();
	run();
}

GillespieSolver::~GillespieSolver() {

}

void GillespieSolver::setRandomGenerator() {
	extern const gsl_rng_type *gsl_rng_default;

	const gsl_rng_type * T;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	timespec ts;

	#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
	//cout << "######## MAC #######" << endl;
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	ts.tv_sec = mts.tv_sec;
	ts.tv_nsec = mts.tv_nsec;

        #else
	clock_gettime(CLOCK_REALTIME, &ts);
        #endif

	gsl_rng_set(r, (int) ts.tv_nsec);
}

int GillespieSolver::chooseHazard(ColumnVector hazards) {
	//Chooses a reaction by taking a uniformaly distributed random number and comparing its value to the normalised hazard of the reactions
	double hazard_chooser = gsl_rng_uniform(r);

	int i;
	double totalhazard = 0;
	double a = 0;
	for (i = 0; i < model->NREACTIONS; ++i) {
		totalhazard = totalhazard + hazards(i + 1);
	}
	for (i = 0; i < model->NREACTIONS; ++i) {
		if (hazards(i + 1) != 0) {
			a += hazards(i + 1) / totalhazard;
		}
		if (hazard_chooser < a) {
			break;
		}
	}
	return i;
}

void GillespieSolver::changeConcentration(int chosenHazard) {
	int i;
	
	for (i = 0; i < model->NSPECIES; i++) {

	  concentrations[i] = concentrations[i] + (*model->pstoichiometricMatrix)(i+1, chosenHazard+1);
	}
}

double GillespieSolver::changeTime(ColumnVector hazards, double time) {
	int i;
	double totalhazard = 0;
	for (i = 0; i < model->NREACTIONS; ++i) {
		totalhazard = totalhazard + hazards(i + 1);
	}
	double b = gsl_ran_exponential(r, 1/totalhazard);

	time = time + b;
	return time;
}

void GillespieSolver::solve() {
	//initialising the values
	ColumnVector hazards = model->getHazards(concentrations, parameters);
	int flag = 1;
	double time = 0;
	double t = timepoints[0];
	counter = 0;
	setInitialValues();
	
	//As long as t is in the timepoints array
	while (flag==1) {
		//As long as time is less than the current timepoint t
		while (time <= t) {
			//1- We change the concentrations and parameters depending on the time and rules and events
			model->applyRulesAndEvents(concentrations, parameters, time);

			//cout << "conc:\t" << concentrations[0] << endl; 

			//2- We compute the hazards for the current concentrations and parameters
			hazards = model->getHazards(concentrations, parameters);
			
			// Check they are positive
			for(int k=0;k<model->NREACTIONS; ++k) hazards(k+1) = hazards(k+1) < 0 ? 0 : hazards(k+1);

			//3- We choose one reaction for this step
			int chosenHazard = chooseHazard(hazards);
			
			//4- We change the concentration given the chosen reaction
			changeConcentration(chosenHazard);

			//5- We change the time for the next step
			time = changeTime(hazards, time);
		}
		//once we reached the timestep t we store the computed value and we increment the counter
		storePoints(t);
		counter++;
		
		if (counter >= timepoints.size()) {
			flag = 0;
		} else {
			t = timepoints[counter];
		}
	}
}
