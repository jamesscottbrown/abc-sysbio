/*
 * SDESolver.cpp
 *
 *  Created on: Mar 14, 2011
 *      Author: student1
 */
#include "SDESolver.hpp"

SDESolver::SDESolver(double* ainitialValues, double* aparameters,
		ChildModel* mmodel, int inbOfIterations, vector<double> vtimepoints,
		string sfilename, double ddt) :
	Solver(ainitialValues, aparameters, mmodel, inbOfIterations, vtimepoints,
			sfilename) {
	dt = ddt;
	maxT = (int) timepoints[timepoints.size() - 1];

	dy = new double[model->NSPECIES];
	dyODE = new double[model->NSPECIES];
	dyEuler = new double[model->NSPECIES];
	setRandomGenerator();
}

SDESolver::~SDESolver() {
	delete [] dy;
	delete [] dyODE;
	delete [] dyEuler;
}

ColumnVector SDESolver::getDriftVector(const ColumnVector& hazards) {
	ColumnVector driftVector(model->NREACTIONS);

	driftVector = (*model->pstoichiometricMatrix) * hazards;
	return driftVector;
}

Matrix SDESolver::getDiffusionMatrix(const ColumnVector& hazards) {

	DiagonalMatrix diagSqrtH(model->NREACTIONS);
	diagSqrtH = 0;

	Matrix diffusionMatrix(model->NSPECIES, model->NREACTIONS);

	for (int i = 1; i <= model->NREACTIONS; i++) {
		diagSqrtH(i) = sqrt(hazards(i));
	}
	diffusionMatrix = (*model->pstoichiometricMatrix) * diagSqrtH;

	return diffusionMatrix;
}

double* SDESolver::getODEterm(const ColumnVector& hazards) {
	double* local_dyODE = new double[model->NSPECIES];
	ColumnVector driftVector = getDriftVector(hazards);
	int k;
	for (k = 0; k < model->NSPECIES; k++) {
		local_dyODE[k] = driftVector(k + 1) * dt;
	}
	return local_dyODE;
}

void SDESolver::setRandomGenerator() {
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

int SDESolver::storeWantedPoints() {
	if (counter >= timepoints.size()) {
		return 1;
	}
	if ((timepoints[counter] - t) < dt) {
		storePoints(timepoints[counter]);
		counter++;
	}
	return 0;
}

void SDESolver::step() {

}

void SDESolver::solve() {
	setInitialValues();
	stop = 0;
	counter = 0;
	t = -dt;
	int stepnumber = 1;
	//while stepnumber is within the number of time steps so while t<=maxT
	while (t+dt<maxT) {
		//first we change the time to the next value
		t = t+dt;
		//Then we apply the rules and event to change the parametes and concentrations in case it's needed
		model->applyRulesAndEvents(concentrations, parameters, t);
		//Then we compute the step
		step();
		//and we check the validity of the dy given by the step, re-do the step if need be and increment the concentrations
		changeConc();
		//finally we check if the counter is not too high and we store the data in the output array if the time matches one of the timepoints
		stop = storeWantedPoints();
		if (stop == 1) {
			break;
		}
		stepnumber++;
	}
}
