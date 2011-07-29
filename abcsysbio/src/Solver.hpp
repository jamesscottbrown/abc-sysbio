/*
 * Solveur.hpp
 * Base class defining a Solver
 *
 *  Created on: Mar 8, 2011
 *      Author: student1
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_
#include "ChildModel.hpp"
#include <fstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "newmat.h"
#include "newmatio.h"
#include "newmatap.h"

#ifdef __MACH__
#include <mach/mach.h>
#include <mach/clock.h>
#endif

class Solver {
public:
	/**
	 * Array of size NSPECIES containing the initial concentrations of the system to solve (doubles)
	 */
	double* initialValues;

	/**
	 * Array of doubles containing the initial values of the parameters of the system to solve. The size depend on the model to solve.
	 * It remains constant during the whole life of the object of type Solver
	 */
	double* initialParameters;

	/**
	 * Array containing the values of the parameter while the system is being solved (in case some rule or event would change them at some point)
	 */
	double* parameters;

	/**
	 * Array of size NSPECIES containing the values of the concentrations of the different species during the simulation
	 */
	double* concentrations;

	/**
	 * Model that we want to solve ie an instance of a child of the base class Model that corresponds to the biological model given by the SBML file
	 */
	ChildModel* model;

	/**
	 * Number of iteration of the same simulation that the user wants to do
	 */
	int nbOfIterations;

	/**
	 * vector containing the times at which the users wants to have the data
	 */
	vector<double> timepoints;

	/**
	 * Number of the current iteration running
	 */
	int iterationNumber;

	/**
	 * Number of the current step of the simulation running
	 */
	unsigned int counter;

	/**
	 * Tri-dimensional array where the data simulated are stored and that then can be used as output.
	 * First dimension : size nbOfIterations, contains the iteration number
	 * Second dimension : size of the timepoints array , contains the timepoint
	 * third dimension : size NSPECIES, contains the concentrations of the different species at the related timepoint and iteration
	 */
	double*** output;

	/**
	 * output stream to rpint the data in a file after simulation. (TODO : get rid of it once all the plots made and all)
	 */
	ofstream out;

	/**
	 * Constructor - initialises all the attributes
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
	Solver(double* ainitialValues, double* aparameters, ChildModel* mmodel, int inbOfIterations, vector<double> vtimepoints, string sfilename);

	/**
	 * Destructor
	 *
	 * @param void
	 * @return void
	 */
	virtual ~Solver();

	/**
	 * Method allowing to fill the arrays concentrations and parameters with their initial values thus (re)initialising the system
	 *
	 * @param void
	 * @return void
	 */
	void setInitialValues();

	/**
	 * Method printing the values stored in the output array in the file defined by filename by the user at the end of all the simulations
	 *
	 * @param void
	 * @return void
	 */
	void storeDataInFile();

	/**
	 * Method storing the current concentrations, timepoint and iteration number in the output array
	 *
	 * @param double t timepoints that we wish to store
	 * @return void
	 */
	void storePoints(double t);

	/**
	 * virtual method doing one simulation of the system
	 *
	 * @param void
	 * @return void
	 */
	virtual void solve() = 0;

	/**
	 * Method running the simulation as many times as asked by the user and dealing with the printing in the file
	 *
	 * @param void
	 * @return void
	 */
	void run();
};

#endif /* SOLVEUR_HPP_ */

