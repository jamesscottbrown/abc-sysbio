/*
 * Model.hpp
 * Base Class defining a biological model
 * The children will be defined by the parser and take their values from the SBML model
 *
 *  Created on: Mar 7, 2011
 *      Author: student1
 */

#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <vector>
#include <iostream>
#include "newmat.h"
#include "newmatio.h"
#include "newmatap.h"

class Model {
public:

  //int blank;
  //Matrix* otherMatrix; 
	/**
	 * Number of reactions of the model
	 */
	int NREACTIONS;

	/**
	 * Number of species of the model
	 */
	int NSPECIES;

	/**
	 * Stoichiometric Matrix of the system (the rows represent the species and the columns the reactions)
	 */
	Matrix* pstoichiometricMatrix;

	/**
	 * Constructor
	 *
	 * @param void
	 * @return void
	 */
	Model(int i);

	/**
	 * Destructor
	 */
	virtual ~Model();

	/**
	 * Virtual method (ie method defined in the child class) setting the values of the stoichiometric matrix
	 *
	 * @param void
	 * @return void
	 */
	virtual void getStoichiometricMatrix() =0;

	/**
	 * Virtual method computing the hazards of the different reactions for a given concentration of species (yi) and some parameter values
	 *
	 * @param double concentrations[] Array of size NSPECIES containing the concentrations of the species for which we want to compute the hazards
	 * @param double parameters[] Array containing the parameter's values for which we want to compute the hazards (the number of parameters depend on the model and doesn't have to be the number of reactions)
	 */
	virtual ColumnVector getHazards(const double concentrations[],
			const double parameters[]) =0;

	/**
	 * Virtual method modifying the concentrations and parameters depending on some criteria defined by the SBML
	 *
	 * @param double concentrations[] Array of size NSPECIES containing the concentrations of the species
	 * @param double parameters[] Array containing the parameter's values
	 */
	virtual void applyRulesAndEvents(double concentrations[],
			double parameters[], double time) =0;


};

#endif /* MODEL_HPP_ */
