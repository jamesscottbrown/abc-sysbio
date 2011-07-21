/*
 * EulerSDESolver.hpp
 * Class defining an SDE Solver using the Euler approximation
 *
 *  Created on: Mar 15, 2011
 *      Author: student1
 */

#ifndef EULERSDESOLVER_HPP_
#define EULERSDESOLVER_HPP_
#include "SDESolver.hpp"

class EulerSDESolver: public SDESolver {
public:

	/**
	 * Constructor - Runs the simulations
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
	EulerSDESolver(double* ainitialValues, double* aparameters, ChildModel* mmodel,
			int inbOfIterations, vector<double> vtimepoints, string filename,
			double ddt);

	/**
	 * Destructor
	 */
	~EulerSDESolver();

	/**
	 * Method returning a vector of size NREACTIONS giving the noise for each reaction that will be used for this step
	 *
	 * @param void
	 * @return vector<double> noise vector of size NREACTIONS giving the noise for each reaction that will be used for this step
	 */
	vector<double> getNoise();

	/**
	 * Method returning the Euler approximation part of the derivative of the concentration
	 *
	 * @param ColumnVector hazards Array of size NREACTIONS containing the hazards of all the reactions of the system
	 * @param vector<double>& noise vector of size NREACTIONS containing the noises associated at each reaction
	 * @return void
	 */
	void getEulerNoiseTerm(const ColumnVector& hazards,
			const vector<double>& noise);

	/**
	 * Method checking at the end of each step the validity of the output (sometimes the steps has to be done again to have a good value because of the noise part)
	 *
	 * @param double* candidateDy the derivative we want to check the validity of
	 * @return 0 if the step has to be done again and 1 if the candidateDy is fine and can be used to increment the concentrations
	 */
	int checkY(double* candidateDy);

	/**
	 * Increment the concentrations at the end of one step
	 *
	 * @param void
	 * @return void
	 */
	void changeConc();

	/**
	 * Method doing one simulation of the system
	 *
	 * @param void
	 * @return void
	 */
	void step();
};

#endif /* EULERSDESOLVER_HPP_ */
