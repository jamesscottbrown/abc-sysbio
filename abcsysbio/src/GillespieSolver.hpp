/*
 * GillespieSolver.hpp
 *
 *  Created on: Mar 15, 2011
 *      Author: student1
 */

#ifndef GILLESPIESOLVER_HPP_
#define GILLESPIESOLVER_HPP_

#include "Solver.hpp"
#include "ChildModel.hpp"
#include <cmath>

#ifdef __MACH__
#include <mach/mach.h>
#include <mach/clock.h>
#endif

class GillespieSolver: public Solver {
public:
	/**
	 * Random generator initialiser
	 */
	gsl_rng * r;

	/**
	 * Constructor - Runs the simulations
	 *
	 * @param double* ainitialValues array of size NSPECIES containing the initial values wanted to solve the system
	 * @param aparameters double* array of the parameters values wanted to solve the system
	 * @param Model* mmodel model we want to solve
	 * @param int inbOfIterations number of time we want to simulate the system
	 * @param vector<double> vtimepoints vector of the timepoints when we want the data to be computed and stored
	 * @param string sfilename name of the file we want to print the output
	 *
	 * @return void
	 */
	GillespieSolver(double* ainitialValues, double* aparameters, ChildModel* mmodel,
			int inbOfIterations, vector<double> timepoints, string sfilename);

	/**
	 * Generator
	 */
	~GillespieSolver();

	/**
	 * Method locking at all the hazards and chooses one (the more the hazard of one reaction is high the higher the probability for it to be choosen)
	 *
	 * @param ColumnVector hazards the hazards of each reaction
	 * @return int the number of the chosen hazards (or reaction...)
	 */
	int chooseHazard(ColumnVector hazards);

	/**
	 * Initialises the random generator
	 */
	void setRandomGenerator();

	/**
	 * method changing the concentration given the hazard that's been chosen
	 *
	 * @param int chosenHazard the number of the hazard or reaction that's been chosen
	 * @return void
	 */
	void changeConcentration(int chosenHazard);

	/**
	 * Method incrementing the time with an random number from an exponential distribution of lambe the sum of the hazards
	 *
	 * @param ColumnVector Column Vector of size NREACTION with the hazards the hazards of the different functions
	 * @param double time the current time to be incremented
	 */
	double changeTime(ColumnVector hazards, double time);

	/**
	 * Method doing one simulation of the system
	 *
	 * @param void
	 * @return void
	 */
	void solve();
};
#endif /* GILLESPIESOLVER_HPP_ */
