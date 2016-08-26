import pickle
import re
from PriorType import PriorType


def checkInputABC(info_new, fname, custom_distance, design):
    """
    Check that the information in the input file is consistent with each other and with the model,
    and that it is in the format required to run the abc-SMC algorithm.
    Return boolean, string (empty if boolean is True)

    Takes as input an algorithm_info object, the output folder name, and whether custom distance is specified
    
    """

    restart = info_new.restart
    model_name = info_new.name
    data = info_new.data
    timepoints = info_new.times
    num_particles = info_new.particles
    epsilon = info_new.epsilon
    integration_type = info_new.type
    model_weights = info_new.modelprior
    priors = info_new.prior
    x0priors = info_new.x0prior
    source = info_new.source
    fit = info_new.fit
    model_kernel = info_new.modelkernel

    # check general properties of the given arguments that are independent of the model
    for i in range(0, len(model_name)):
        if model_name[i] == "":
            return False, "\nPlease do not give empty strings for model names!\n"

    if not len(model_name) == len(integration_type):
        return False, "\nPlease provide the same amount of model sources and integration types!\n"

    if not len(model_name) == len(model_weights):
        return False, "\nPlease provide the same amount of model sources and model weights!\n"

    if not len(model_name) == len(fit):
        return False, "\nPlease provide a fit instruction (or None) for each model. If the fit instruction is None all data will be fitted to the model data.\n"

    if not design:
        if not len(data) == len(timepoints):
            return False, "\nPlease provide data that correspond to the length of your timepoints!\n"

    if not len(model_name) == len(priors):
        return False, "\nPlease provide prior distributions for each model!\n"

    if not len(model_name) == len(x0priors):
        return False, "\nPlease provide initial values for each model!\n"

    if not model_kernel > 0.0:
        return False, "\nPlease provide a model Kernel larger than 0!\n"

    if len(epsilon) > 1 and (custom_distance is False):
        return False, "\nPlease provide a custom distance function when you specify more than 1 epsilon schedule!\n"

    # check model specific properties (comparing with SBML model)
    if source is not None:
        import libsbml

        if not len(source) == len(model_name):
            return False, "\nPlease provide the same amount of model sources and model names!\n"

        reader = libsbml.SBMLReader()
        for mod in range(0, len(source)):

            document = reader.readSBML(source[mod])
            model = document.getModel()

            num_species = model.getNumSpecies()
            num_global_parameters = model.getNumParameters()
            list_of_parameters = []

            num_compartments = model.getNumCompartments()

            for i in range(0, num_compartments):
                if model.getCompartment(i).isSetVolume():
                    num_global_parameters += 1
                    list_of_parameters.append(model.getListOfCompartments()[i])

            for i in range(0, num_global_parameters - num_compartments):
                list_of_parameters.append(model.getParameter(i))

            num_local_parameters = 0
            for i in range(0, model.getNumReactions()):
                local = model.getReaction(i).getKineticLaw().getNumParameters()
                num_local_parameters = num_local_parameters + local
                for k in range(0, local):
                    list_of_parameters.append(model.getListOfReactions()[i].getKineticLaw().getParameter(k))

            num_parameters = num_local_parameters + num_global_parameters

            list_of_rules = model.getListOfRules()
            for k in range(0, len(list_of_parameters)):
                if not list_of_parameters[k].getConstant():
                    for j in range(0, len(list_of_rules)):
                        if list_of_rules[j].isRate():
                            if list_of_parameters[k].getId() == list_of_rules[j].getVariable():
                                num_species += 1
                                num_parameters -= 1

            if not len(priors[mod]) == num_parameters:
                return False, "\nThe number of given prior distributions for model " + model_name[
                    mod] + " is not correct!\n"

            if not len(x0priors[mod]) == num_species:
                return False, "\nPlease provide an initial value for each species in model " + model_name[mod] + "!\n"

    # checking further properties independent of the model
    sde = re.compile('SDE')
    ode = re.compile('ODE')
    gillespie = re.compile('Gillespie')

    for mod in range(0, len(model_name)):

        string = integration_type[mod]
        if (not sde.search(string)) and (not ode.search(string)) and (not gillespie.search(string)):
            return False, "\nThe integration type for model " + model_name[mod] + " does not exist!\n"

        for ic in range(len(x0priors[mod])):
            if x0priors[mod][ic][0] not in [PriorType.constant, PriorType.normal, PriorType.uniform, PriorType.lognormal]:
                return False, "\nThe prior distribution of initial condition " + repr(ic + 1) + " in model " + \
                       model_name[mod] + " does not exist!\n"

        for param in range(0, len(priors[mod])):
            if not len(priors[mod][param]) == 3:
                return False, "\nThe prior distribution of parameter " + repr(param + 1) + " in model " + model_name[
                    mod] + " is wrong defined!\n"

            if priors[mod][param][0] not in [PriorType.constant, PriorType.normal, PriorType.uniform, PriorType.lognormal]:
                return False, "\nThe prior distribution of parameter " + repr(param + 1) + " in model " + model_name[
                    mod] + " does not exist!\n"

            if priors[mod][param][0] == PriorType.uniform:
                if not priors[mod][param][1] < priors[mod][param][2]:
                    return False, "\nThe range of the uniform prior distribution of parameter " + repr(
                        param + 1) + " in model " + model_name[mod] + " is wrong defined!\n"

            if priors[mod][param][0] == PriorType.lognormal:
                if not (priors[mod][param][1] >= 0 or priors[mod][param][2] >= 0):
                    return False, "\nThe mean or scale of the lognormal prior distribution of parameter " + repr(
                        param + 1) + " in model " + model_name[mod] + " is wrong defined!\n"

    # check arguments connected to pickling
    if restart:

        try:
            in_file = open(fname + '/copy/algorithm_parameter.dat', "r")
            num_particles_pickled = pickle.load(in_file)
            in_file.close()
        except IOError:
            return False, "\nCan not find file \'algorithm_parameter.dat\' in folder \'copy\'!\n"

        if not num_particles <= num_particles_pickled:
            return False, "\nRunning the abc algorithm from a previous point is not possible with a larger population size!\n"

    return True, ""


def checkInputSimulation(info_new):
    """
    Check that the information in the input file is consistent with each other and with the model,
    and that it is in the format required to simulate the model.
    Return boolean, string (empty if boolean is True)
    
    """

    name = info_new.name
    timepoints = info_new.times
    integration_type = info_new.type
    source = info_new.source

    for i in range(0, len(name)):
        if name[i] == "":
            return False, "\nPlease do not give empty strings for model names!\n"

    if not len(name) == len(integration_type):
        return False, "\nPlease provide the same number of model names and and integration types!\n"

    # check model specific properties (comparing with SBML model)
    if source is not None:
        import libsbml

        if not len(source) == len(name):
            return False, "\nPlease provide the same amount of model sources and model names!\n"

        reader = libsbml.SBMLReader()
        for mod in range(0, len(source)):

            document = reader.readSBML(source[mod])
            model = document.getModel()

            num_species = model.getNumSpecies()
            num_global_parameters = model.getNumParameters()

            num_compartments = model.getNumCompartments()
            for i in range(0, num_compartments):
                if model.getCompartment(i).isSetVolume():
                    num_global_parameters += 1

            num_local_parameters = 0
            for i in range(0, model.getNumReactions()):
                num_local_parameters = num_local_parameters + model.getReaction(i).getKineticLaw().getNumParameters()

            num_parameters = num_local_parameters + num_global_parameters

            if not info_new.nparameters[mod] == num_parameters:
                return False, "\nThe number of given parameters for model " + name[mod] + " is not correct!\n"

            if not len(info_new.x0prior[mod]) == num_species:
                return False, "\nPlease provide an initial value for each species in model " + name[mod] + "!\n"

    if len(timepoints) == 0:
        return False, "\nPlease give timepoints at which to return simulated data points!\n"

    sde = re.compile('SDE')
    ode = re.compile('ODE')
    gillespie = re.compile('Gillespie')

    for mod in range(0, len(name)):

        string = integration_type[mod]
        if (not sde.search(string)) and (not ode.search(string)) and (not gillespie.search(string)):
            return False, "\nThe integration type for model " + str(mod + 1) + " does not exist!\n"

    return True, ""
