/*
 * MilsteinSDESolver.hpp
 *
 *  Created on: Mar 14, 2011
 *      Author: student1
 */

#ifndef MILSTEINSDESOLVER_HPP_
#define MILSTEINSDESOLVER_HPP_

#include "SDESolver.hpp"

class MilsteinSDESolver: public SDESolver {
public:

	/**
	 * Array of size NSPECIES containing the Milstein approximation part of the derivatives of the concentrations
	 */
	double* dyMilstein;

	/**
	 * Array of size NSPECIES containing intermediate concentrations used in the algorithm (adaptedYj[k] = concentrations[k]+driftVector[k]*dt+DiffusionMatrix[k,j]*sqrt(dt))
	 */
	double* adaptedYj;

	/**
	 * Array of size NSPECIES containing the concentrations at the step before so as to be able to go back ni case the adaptedHazard gets negativ at the next step
	 */
	double* lastStepConcentrations;

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
	MilsteinSDESolver(double* ainitialValues, double* aparameters,
			ChildModel* mmodel, int inbOfIterations, vector<double> vtimepoints,
			string sfilename, double ddt);

	/**
	 * Destructor
	 */
	~MilsteinSDESolver();

	/**
	 * Method returning a bi dimensionnal vector of size NREACTIONS and 4 giving the 4 noises for each reaction that will be used for this step
	 *
	 * @param void
	 * @return vector<vector<double> > noise vector of size NREACTIONS and 4 giving the noises for each reaction that will be used for this step
	 */
	vector<vector<double> > getNoise();

	/**
	 * Method returning the Euler approximation part of the derivative of the concentration
	 *
	 * @param ColumnVector hazards Array of size NREACTIONS containing the hazards of all the reactions of the system
	 * @param vector<vector<double> >& noise vector of size NREACTIONS containing the noises associated at each reaction
	 * @return void
	 */
	void getEulerNoiseTerm(const ColumnVector& hazards, const vector<vector<
			double> >& noise);

	/**
	 * Method returning the stochastic integral over the two reaction j1 and j2
	 *
	 * @param int j1 first reaction for which we want to compute the integral
	 * @param int j2 second reaction
	 * @param vector<vector<double> >& noise table with the 4 noises for each reaction
	 *
	 * @return double approximated value of the stochastic integral
	 */
	double getDoubleStochasticIntegral(int j1, int j2, const vector<vector<
			double> >& noise);

	/**
	 * Method returning the Milstein approximation part of the derivative of the concentration
	 *
	 * @param Matrix diffusionMatrix The Diffusion Matrix
	 * @param ColumnVector driftVector The Drift Vector
	 * @param vector<vector<double> >& noise vector of size NREACTIONS containing the noises associated at each reaction
	 * @return void
	 */
	void getMilsteinNoiseTerm(const Matrix& diffusionMatrix,
			const ColumnVector& driftVector, const double yColumn[],
			const vector<vector<double> >& noise);

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

	void step();
};

#endif /* MILSTEINSDESOLVER_HPP_ */
