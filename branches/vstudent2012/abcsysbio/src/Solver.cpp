/*
 * Solveur.cpp
 *
 *  Created on: Mar 8, 2011
 *      Author: student1
 */
#include "Solver.hpp"

Solver::Solver(double* ainitialValues, double* aparameters, ChildModel* mmodel,
		int inbOfIterations, vector<double> vtimepoints, string sfilename) {
	initialValues = ainitialValues;
	initialParameters = aparameters;
	concentrations = new double[mmodel->NSPECIES];
	model = mmodel;

	nbOfIterations = inbOfIterations;
	timepoints = vtimepoints;
	//out.open(sfilename.c_str(), ios::out);
	output = new double**[nbOfIterations];

	for (int i = 0; i < nbOfIterations; i++) {
		output[i] = new double*[timepoints.size()];
		for (unsigned int j = 0; j < timepoints.size(); j++) {
			output[i][j] = new double[model->NSPECIES + 1];
		}
	}
	setInitialValues();
}

Solver::~Solver() {
	delete[] concentrations;
	for (int i = 0; i < nbOfIterations; i++) {
		for (unsigned int j = 0; j < timepoints.size(); j++) {
			delete[] output[i][j];
		}
		delete[] output[i];
	}
	delete[] output;
	
}

void Solver::storePoints(double t) {
	output[iterationNumber][counter][0] = t;

	int x;
	for (x = 0; x < model->NSPECIES; x++) {
		output[iterationNumber][counter][x + 1] = (double)concentrations[x];
	}
}

void Solver::storeDataInFile() {
	int k;
	for (unsigned int ii = 0; ii < timepoints.size(); ++ii) {
		out << output[0][ii][0] << ";";
		for (k = 0; k < model->NSPECIES; k++) {
			out << ";" << output[0][ii][k + 1];
		}
		out << endl;
	}
	out.close();
}

void Solver::run() {
	for (iterationNumber = 0; iterationNumber < nbOfIterations; ++iterationNumber) {
		solve();
	}
	//storeDataInFile();
}

void Solver::setInitialValues() {
	for (int k = 0; k < model->NSPECIES; k++) {
		concentrations[k] = initialValues[k];
	}
	parameters = initialParameters;
}

