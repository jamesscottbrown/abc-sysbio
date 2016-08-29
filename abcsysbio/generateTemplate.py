def getSpeciesValue(species):
    """
    Return the initial amount of a species.
    If species.isSetInitialAmount() == True, return the initial amount.
    Otherwise, return the initial concentration.

    ***** args *****
    
    species:    a libsbml.Species object
    
    """

    if species.isSetInitialAmount():
        return species.getInitialAmount()
    else:
        return species.getInitialConcentration()


def generateTemplate(source, filename, sumname, dataname=None):
    """

    Generate a model summary file (model_summary.txt) and a template file (filename) from one or more SBML source files.
    
    
    ***** args *****
    
    source:    a list of strings.
               Each entry describes a SBML file. 


    ***** kwargs *****
    
    filename:  a string.
               The name of the template to be generated.

    sumnname:  a string.
               The name of the summary to be generated.

    dataname:  a string.
               The name of a datafile.

    """

    out_file = open(filename, "w")
    sum_file = open(sumname, "w")

    have_data = False
    times = []
    variables = []
    nvar = 0
    first = True
    if dataname is not None:
        have_data = True
        df = open(dataname, 'r')
        for line in df:
            strs = str(line).split(' ')
            vals = [float(i) for i in strs]

            if first:
                for j in range(1, len(vals)):
                    variables.append([])
                first = False
                nvar = len(vals) - 1

            times.append(vals[0])

            for j in range(1, len(vals)):
                variables[j - 1].append(vals[j])

    # print times
    # print vars

    out_file.write("<input>\n\n")
    out_file.write("######################## number of models\n\n")
    out_file.write("# Number of models for which details are described in this input file\n")
    out_file.write("<modelnumber> " + repr(len(source)) + " </modelnumber>\n\n")

    out_file.write("######################## restart\n\n")
    out_file.write("# Restart from previous (pickled) population?\n")
    out_file.write("<restart> False </restart>\n\n")

    out_file.write("######################## epsilon schedule\n\n")
    out_file.write("# Automatic epsilon schedule. Provide a vector of final epsilons and the alpha (defaults to 0.9)\n")
    out_file.write("<autoepsilon>\n")
    out_file.write("<finalepsilon> 1.0 </finalepsilon>\n")
    out_file.write("<alpha> 0.9 </alpha>\n")
    out_file.write("</autoepsilon>\n\n")

    out_file.write("# OR\n")
    out_file.write("# Series of epsilons. (Whitespace delimited list)\n")
    out_file.write(
        "# Multiple epsilon schedules can be specified by giving additional vectors enclosed in <e2> </e2>, <e3> </e3> etc\n")
    out_file.write("# NOTE: the parser always reads them in order and ignores the tag value\n")
    out_file.write("<!-- <epsilon> -->\n")
    out_file.write("<!-- <e1> 5.0 3.0 2.0 1.0 </e1> -->\n")
    out_file.write("<!--</epsilon> -->\n")
    out_file.write("\n")

    out_file.write("######################## particles\n\n")
    out_file.write("<particles> 100 </particles>\n\n")

    out_file.write("######################## beta\n\n")
    out_file.write(
        "# Beta is the number of times to simulate each sampled parameter set.\n# This is only applicable for models simulated using Gillespie and SDE\n")
    out_file.write("<beta> 1 </beta>\n\n")

    out_file.write("######################## dt\n\n")
    out_file.write("# Internal timestep for solver.\n# Make this small for a stiff model.\n")
    out_file.write("<dt> 0.01 </dt>\n\n")

    out_file.write("######################## perturbation kernels : OPTIONAL (default uniform)\n\n")
    out_file.write("# The pertubation kernels are computed with respect to the previous parameter distribution\n")
    out_file.write("# Currently uniform and normal are implemented\n")
    out_file.write("<kernel> uniform </kernel>\n\n")

    out_file.write("######################## model kernel : OPTIONAL (default 0.7)\n\n")
    out_file.write("# Probability of perturbing the sampled model (ignored when modelnumber = 1)\n")
    out_file.write("<modelkernel> 0.7 </modelkernel>\n\n")

    out_file.write("######################## ODE solver control : OPTIONAL \n\n")
    out_file.write(
        "# rtol and atol can be specified here. If the model is stiff then setting these to small might help the simulation to run\n")
    out_file.write("#<rtol> </rtol> \n#<atol> </atol>\n\n")

    out_file.write("######################## User-supplied data\n\n")
    out_file.write("<data>\n")
    out_file.write("# times: For ABC SMC, times must be a whitespace delimited list\n")
    out_file.write("# In simulation mode these are the timepoints for which the simulations will be output\n")
    if not have_data:
        out_file.write("<times> 0 1 2 3 4 5 6 7 8 9 10 </times>\n\n")
    else:
        out_file.write("<times>")
        for i in times:
            out_file.write(" " + repr(i))
        out_file.write(" </times>\n\n")

    out_file.write(
        "# variables: For ABC SMC, whitespace delimited lists of concentrations (ODE or SDE) or molecule numbers (Gillespie)\n")
    out_file.write(
        "# Denote your data via tags <v1> </v1> or <var1> </var1> or <v2> </v2> etc. The tags are ignored and the data read in order\n")
    out_file.write("# For simulation these data are ignored\n")
    out_file.write(
        "# See fitting instruction below if the dimensionality of your data sets differ from the dimensionality of your model\n")
    out_file.write("<variables>\n")
    if not have_data:
        out_file.write(" <var1> </var1>\n")
    else:
        for k in range(nvar):
            out_file.write("<var" + repr(k + 1) + "> ")
            for i in variables[k]:
                out_file.write(" " + repr(i))
            out_file.write(" </var" + repr(k + 1) + ">\n")

    out_file.write("</variables>\n")
    out_file.write("</data>\n\n")

    out_file.write("######################## Models\n\n")
    out_file.write("<models>\n")

    import libsbml
    reader = libsbml.SBMLReader()

    for i in range(len(source)):
        sum_file.write("Model " + repr(i + 1) + "\n")
        sum_file.write("name: model" + repr(i + 1) + "\nsource: " + source[i] + "\n\n")

        out_file.write("<model" + repr(i + 1) + ">\n")
        out_file.write("<name> model" + repr(i + 1) + " </name>\n<source> " + source[i] + " </source>\n\n")
        out_file.write("# type: the method used to simulate your model. ODE, SDE or Gillespie.\n")
        out_file.write("<type> SDE </type>\n\n")

        out_file.write(
            "# Fitting information. If fit is None, all species in the model are fitted to the data in the order they are listed in the model.\n")
        out_file.write(
            "# Otherwise, give a whitespace delimited list of fitting instrictions the same length as the dimensions of your data.\n")
        out_file.write(
            "# Use speciesN to denote the Nth species in your model. Simple arithmetic operations can be performed on the species from your model.\n")
        out_file.write(
            "# For example, to fit the sum of the first two species in your model to your first variable, write fit: species1+species2\n")
        out_file.write("<fit> None </fit>\n\n")

        document = reader.readSBML(source[i])
        model = document.getModel()

        num_species = model.getNumSpecies()
        num_global_parameters = model.getNumParameters()

        parameter = []
        parameter_id = []
        parameter_id2 = []
        list_of_parameters = []

        r1 = 0
        r2 = 0
        r3 = 0
        list_of_rules = model.getListOfRules()
        for k in range(len(list_of_rules)):
            if model.getRule(k).isAlgebraic():
                r1 += 1
            if model.getRule(k).isAssignment():
                r2 += 1
            if model.getRule(k).isRate():
                r3 += 1

        comp = 0
        num_compartments = model.getNumCompartments()
        for k in range(num_compartments):
            if model.getCompartment(k).isSetVolume():
                comp += 1
                num_global_parameters += 1
                parameter.append(model.getListOfCompartments()[k].getVolume())
                parameter_id.append(model.getListOfCompartments()[k].getId())
                parameter_id2.append('compartment' + repr(k + 1))
                list_of_parameters.append(model.getListOfCompartments()[k])

        for k in range(num_global_parameters - comp):
            param = model.getParameter(k)
            parameter.append(param.getValue())
            parameter_id.append(param.getId())
            parameter_id2.append('parameter' + repr(k + 1))
            list_of_parameters.append(param)

        num_local_parameters = 0
        num_reactions = model.getNumReactions()
        for k in range(num_reactions):
            local = model.getReaction(k).getKineticLaw().getNumParameters()
            num_local_parameters = num_local_parameters + local

            for j in range(local):
                parameter.append(model.getListOfReactions()[k].getKineticLaw().getParameter(j).getValue())
                parameter_id.append(model.getListOfReactions()[k].getKineticLaw().getParameter(j).getId())
                x = len(parameter_id) - comp
                parameter_id2.append('parameter' + repr(x))
                list_of_parameters.append(model.getListOfReactions()[k].getKineticLaw().getParameter(j))

        num_parameters = num_local_parameters + num_global_parameters

        species = model.getListOfSpecies()

        sum_file.write("number of compartments: " + repr(num_compartments) + "\n")
        sum_file.write("number of reactions: " + repr(num_reactions) + "\n")
        sum_file.write("number of rules: " + repr(model.getNumRules()) + "\n")
        if model.getNumRules() > 0:
            sum_file.write("\t Algebraic rules: " + repr(r1) + "\n")
            sum_file.write("\t Assignment rules: " + repr(r2) + "\n")
            sum_file.write("\t Rate rules: " + repr(r3) + "\n\n")
        sum_file.write("number of functions: " + repr(model.getNumFunctionDefinitions()) + "\n")
        sum_file.write("number of events: " + repr(model.getNumEvents()) + "\n\n")

        param_as_species = 0
        sum_file.write("Species with initial values: " + repr(num_species) + "\n")

        out_file.write("# Priors on initial conditions and parameters:\n")
        out_file.write("# one of \n")
        out_file.write("#       constant, value \n")
        out_file.write("#       normal, mean, variance \n")
        out_file.write("#       uniform, lower, upper \n")
        out_file.write("#       lognormal, mean, variance \n\n")

        out_file.write("<initial>\n")

        x = 0
        for k in range(len(species)):
            x += 1
            out_file.write(
                " <ic" + repr(x) + "> constant " + repr(getSpeciesValue(species[k])) + " </ic" + repr(x) + ">\n")
            sum_file.write("S" + repr(x) + ":\t" + species[k].getId() + "\tspecies" + repr(k + 1) + "\t(" + repr(
                getSpeciesValue(species[k])) + ")\n")
        for k in range(len(list_of_parameters)):
            if not list_of_parameters[k].getConstant():
                for j in range(len(list_of_rules)):
                    if list_of_rules[j].isRate():
                        if parameter_id[k] == list_of_rules[j].getVariable():
                            x += 1
                            param_as_species += 1
                            out_file.write(
                                " <ic" + repr(x) + "> constant " + repr(list_of_parameters[k].getValue()) + " </ic" + repr(
                                    x) + ">\n")
                            sum_file.write("S" + repr(x) + ":\t" + list_of_parameters[k].getId() + "\tparameter" + repr(
                                k + 1 - comp) + "\t(" + repr(list_of_parameters[
                                                                 k].getValue()) + ") (parameter included in a rate rule and therefore treated as species)\n")

        out_file.write("</initial>\n\n")

        sum_file.write("\n")

        if num_global_parameters == 0:
            string = " (all of them are local parameters)\n"
        elif num_global_parameters == 1:
            string = " (the first parameter is a global parameter)\n"
        elif num_local_parameters == 0:
            string = " (all of them are global parameters)\n"
        else:
            string = " (the first " + repr(num_global_parameters) + " are global parameter)\n"

        sum_file.write("Parameter: " + repr(num_parameters) + string)
        sum_file.write("(" + repr(param_as_species) + " parameter is treated as species)\n")

        out_file.write("<parameters>\n")

        counter = 0
        for k in range(num_parameters - param_as_species):
            do_print = True
            if k < len(list_of_parameters):
                if not list_of_parameters[k].getConstant():
                    for j in range(len(list_of_rules)):
                        if list_of_rules[j].isRate():
                            if parameter_id[k] == list_of_rules[j].getVariable():
                                do_print = False

            if do_print:
                counter += 1
                sum_file.write("P" + repr(counter) + ":\t" + parameter_id[k] + "\t" + parameter_id2[k] + "\t(" + repr(
                    parameter[k]) + ")\n")
                out_file.write("<parameter" + repr(counter) + ">")
                out_file.write(" constant ")
                out_file.write(repr(parameter[k]) + " </parameter" + repr(counter) + ">\n")

        sum_file.write("\n############################################################\n\n")

        out_file.write("</parameters>\n")
        out_file.write("</model" + repr(i + 1) + ">\n\n")

    out_file.write("</models>\n\n")
    out_file.write("</input>\n\n")

    out_file.close()
    sum_file.close()
