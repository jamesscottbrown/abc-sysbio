/*
 * ODESolver.hpp
 * Class defining a Solveur using ODEs to solve the system
 *  Created on: Mar 10, 2011
 *      Author: student1
 */

#ifndef ODESOLVER_HPP_
#define ODESOLVER_HPP_

#include "Solver.hpp"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>


class ODESolver: public Solver {
public:

	/**
	 * The initial step size used by the solver, but it is then changed automatically by the solver to fit best the system
	 */
	double initialStepSize;

	/**
	 * The controls implemented will keep the local error on each step within an absolute error of absoluteError and relative error of relativeError with respect to the solution
	 */
	double absoluteError;

	/**
	 * The controls implemented will keep the local error on each step within an absolute error of absoluteError and relative error of relativeError with respect to the solution
	 */
	double relativeError;

	/**
	 * Constructor - Runs the simulations
	 *
	 * @param double* ainitialValues array of size NSPECIES containing the initial values wanted to solve the system
	 * @param aparameters double* array of the parameters values wanted to solve the system
	 * @param Model* mmodel model we want to solve
	 * @param int inbOfIterations number of time we want to simulate the system
	 * @param vector<double> vtimepoints vector of the timepoints when we want the data to be computed and stored
	 * @param string sfilename name of the file we want to print the output
	 * @param double dinitialStepSize The initial step size used by the solver, but it is then changed automatically by the solver to fit best the system
	 * @param double dabsoluteError The controls implemented will keep the local error on each step within an absolute error of absoluteError and relative error of relativeError with respect to the solution
	 * @param double drelativeError The controls implemented will keep the local error on each step within an absolute error of absoluteError and relative error of relativeError with respect to the solution
	 * @return void
	 */
	ODESolver(double* ainitialValues, double* aparameters, ChildModel* mmodel,
			int inbOfIterations, vector<double> vtimepoints, string sfilename,
			double dinitialStepSize = 1e-6, double dabsoluteError = 1e-6,
			double drelativeError = 0.0);

	/**
	 * Destructor
	 */
	~ODESolver();

	/**
	 * Method doing one simulation of the system
	 *
	 * @param void
	 * @return void
	 */
	void solve();

	/**
	 * Method giving the values of each reaction to gsl when the solver needs it
	 */
	double* reactions(const double concentrations[]);
};

#endif /* ODESOLVER_HPP_ */
