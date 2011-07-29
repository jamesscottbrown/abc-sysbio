/*
 * SDESolver.hpp
 * Defines the base class SDESolver, that inheritates from the base class Solver and that deals with the methods and attributes defining any SDESolver whichever algorithm is used to solve it
 *  Created on: Mar 14, 2011
 *      Author: student1
 */

#ifndef SDESOLVER_HPP_
#define SDESOLVER_HPP_

#include "Solver.hpp"
#include <cmath>

#ifdef __MACH__
#include <mach/mach.h>
#include <mach/clock.h>
#endif

class SDESolver: public Solver {
public:
	/**
	 * time step wanted to solve the system
	 */
	double dt;

	/**
	 * Maximum time asked by the user (it is the higher number in the timepoints array)
	 */
	int maxT;

	/**
	 * Array of size NSPECIES - Derivative of the concentrations at the current moment
	 */
	double* dy;

	/**
	 * Array of size NSPECIES containing the ODE part of the derivatives of the concentrations
	 */
	double* dyODE;

	/**
	 * Array of size NSPECIES containing the Euler approximation part of the derivatives of the concentrations
	 */
	double* dyEuler;

	/**
	 * current time of the simulation
	 */
	double t;

	/**
	 * Flag to stop the simulation when we reach the end of the timepoints wanted.
	 */
	int stop;

	/**
	 * random generator initialisator
	 */
	gsl_rng * r;

	/**
	 * Constructor - initialises the attributes specific to an SDE solver and the random generator
	 *
	 * @param double* ainitialValues array of size NSPECIES containing the initial values wanted to solve the system
	 * @param aparameters double* array of the parameters values wanted to solve the system
	 * @param Model* mmodel model we want to solve
	 * @param int inbOfIterations number of time we want to simulate the system
	 * @param vector<double> vtimepoints vector of the timepoints when we want the data to be computed and stored
	 * @param string sfilename name of the file we want to print the output
	 * @param double ddt time step wanted
	 * @return void
	 */
	SDESolver(double* ainitialValues, double* aparameters, ChildModel* mmodel,
			int inbOfIterations, vector<double> vtimepoints, string sfilename,
			double ddt);

	/**
	 * Destructor
	 */
	~SDESolver();

	/**
	 * Method returning the drift vector
	 *
	 * @param ColumnVector hazards Array of size NREACTIONS containing the hazards of all the reactions of the system
	 * @return ColumnVector Drift vector
	 */
	ColumnVector getDriftVector(const ColumnVector& hazards);

	/**
	 * Method returning the diffusion matrix
	 *
	 * @param ColumnVector hazards Array of size NREACTIONS containing the hazards of all the reactions of the system
	 * @return Matrix the diffusion Matrix
	 */
	Matrix getDiffusionMatrix(const ColumnVector& hazards);

	/**
	 * Method returning the ODE part of the derivative of the concentration (overwrites dyODE)
	 *
	 * @param ColumnVector hazards Array of size NREACTIONS containing the hazards of all the reactions of the system
	 * @return  Array of size NSPECIES containing the ODE part of the derivatives of the concentrations used to compute the hazards
	 */
	double* getODEterm(const ColumnVector& hazards);

	/**
	 * Virtual method computing one step of the simulation
	 * (differs given the algorithm implemented)
	 *
	 * @param void
	 * @return void
	 */
	virtual void step();

	/**
	 * Virtual method checking at the end of each step the validity of the output (sometimes the steps has to be done again to have a good value because of the noise part)
	 * (differs given the algorithm implemented)
	 *
	 * @param double* candidateDy the derivative we want to check the validity of
	 * @return 0 if the step has to be done again and 1 if the candidateDy is fine and can be used to increment the concentrations
	 */
	virtual int checkY(double* candidateDy)=0;

	/**
	 * Method filling the times array
	 *
	 * @param int maxT maximum times to reach during the simulation (given by the maximum of the user input timepoints vector)
	 * @return void
	 */
	void setTimes(int maxT);

	/**
	 * Increment the concentrations at the end of one step
	 *
	 * @param void
	 * @return void
	 */
	virtual void changeConc()=0;

	/**
	 * Initialises the random generator
	 *
	 * @param void
	 * @return void
	 */
	void setRandomGenerator();

	/**
	 * Stores the points in the output array when their time fits with the timepoints
	 *
	 * @param void
	 * @return int Stop flag saying if the simulation should go on or stop (returns 1 if the counter is superior or equals to the number of timepoints meaning we reached the end of the timepoints)
	 */
	int storeWantedPoints();

	/**
	 * Virtual method doing one simulation of the system
	 *
	 * @param void
	 * @return void
	 */
	void solve();
};

#endif /* SDESOLVER_HPP_ */
