import libsbml


def get_species_value(species):
    """
    Return the initial amount of a species.
    If species.isSetInitialAmount() == True, return the initial amount.
    Otherwise, return the initial concentration.

    Parameters
    ----------
    species : a libsbml.Species object

    """

    if species.isSetInitialAmount():
        return species.getInitialAmount()
    else:
        return species.getInitialConcentration()


def generate_template(source, filename, sumname, dataname=None):
    """

    Generate a model summary file (model_summary.txt) and a template file (filename) from one or more SBML source files.
    
    
    Parameters
    ----------
    source :     a list of strings, each the name of a SBML file.
    filename :   a string containing name of the template to be generated.
    sumname :   a string containing the name of the summary to be generated.
    dataname :   a string containing name of a datafile.

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

    out_file.write("""
<input>


######################## number of models
# Number of models for which details are described in this input file
""")

    out_file.write("<modelnumber> %s </modelnumber>" % repr(len(source)))

    out_file.write("""
######################## restart
# Restart from previous (pickled) population?
<restart> False </restart>

######################## epsilon schedule
# Automatic epsilon schedule. Provide a vector of final epsilons and the alpha (defaults to 0.9)
<autoepsilon>
<finalepsilon> 1.0 </finalepsilon>
<alpha> 0.9 </alpha>
</autoepsilon>

# OR
# Series of epsilons. (Whitespace delimited list)
# Multiple epsilon schedules can be specified by giving additional vectors enclosed in <e2> </e2>, <e3> </e3> etc
# NOTE: the parser always reads them in order and ignores the tag value
<!-- <epsilon> -->
<!-- <e1> 5.0 3.0 2.0 1.0 </e1> -->
<!--</epsilon> -->


######################## particles
<particles> 100 </particles>

######################## beta
# Beta is the number of times to simulate each sampled parameter set.
# This is only applicable for models simulated using Gillespie and SDE
<beta> 1 </beta>

######################## dt
# Internal timestep for solver.
# Make this small for a stiff model.
<dt> 0.01 </dt>

######################## perturbation kernels : OPTIONAL (default uniform)
# The pertubation kernels are computed with respect to the previous parameter distribution
# Currently uniform and normal are implemented
<kernel> uniform </kernel>

######################## model kernel : OPTIONAL (default 0.7)
# Probability of perturbing the sampled model (ignored when modelnumber = 1)
<modelkernel> 0.7 </modelkernel>

######################## ODE solver control : OPTIONAL
# rtol and atol can be specified here. If the model is stiff then making these small might help the simulation to run
#<rtol> </rtol>
#<atol> </atol>

######################## User-supplied data
<data>
# times: For ABC SMC, times must be a whitespace delimited list
# In simulation mode these are the timepoints for which the simulations will be output
""")    

    if not have_data:
        out_file.write("<times> 0 1 2 3 4 5 6 7 8 9 10 </times>\n\n")
    else:
        out_file.write("<times>")
        for i in times:
            out_file.write(" " + repr(i))
        out_file.write(" </times>\n\n")

    out_file.write("""
# variables: For ABC SMC, whitespace delimited lists of concentrations (ODE or SDE) or molecule numbers (Gillespie)
# Denote your data via tags <v1> </v1> or <var1> </var1> or <v2> </v2> etc.
# The tags are ignored and the data read in order
# For simulation these data are ignored
# See fitting instruction below if the dimensionality of your data sets differ from the dimensionality of your model
<variables>
""")

    if not have_data:
        out_file.write(" <var1> </var1>\n")
    else:
        for k in range(nvar):
            out_file.write("<var" + repr(k + 1) + "> ")
            for i in variables[k]:
                out_file.write(" " + repr(i))
            out_file.write(" </var" + repr(k + 1) + ">\n")

    out_file.write("""
</variables>
</data>

######################## Models
<models>
""")

    reader = libsbml.SBMLReader()

    for i in range(len(source)):
        sum_file.write("Model " + repr(i + 1) + "\n")
        sum_file.write("name: model" + repr(i + 1) + "\nsource: " + source[i] + "\n\n")

        out_file.write("<model" + repr(i + 1) + ">\n")
        out_file.write("<name> model" + repr(i + 1) + " </name>\n<source> " + source[i] + " </source>\n\n")
        out_file.write("# type: the method used to simulate your model. ODE, SDE or Gillespie.\n")
        out_file.write("<type> SDE </type>\n\n")

        out_file.write("""
# Fitting information.
# If fit is None, all species in the model are fitted to the data in the order they are listed in the model.
# Otherwise, give a whitespace delimited list of fitting instrictions the same length as the dimensions of your data.
# Use speciesN to denote the Nth species in your model. Simple arithmetic operations can be performed on these.
# For example, to fit the sum of the first two species in the model to your first variable, write fit: species1+species2
<fit> None </fit>

""")

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
                " <ic" + repr(x) + "> constant " + repr(get_species_value(species[k])) + " </ic" + repr(x) + ">\n")
            sum_file.write("S" + repr(x) + ":\t" + species[k].getId() + "\tspecies" + repr(k + 1) + "\t(" + repr(
                get_species_value(species[k])) + ")\n")
        for k in range(len(list_of_parameters)):
            if not list_of_parameters[k].getConstant():
                for j in range(len(list_of_rules)):
                    if list_of_rules[j].isRate():
                        if parameter_id[k] == list_of_rules[j].getVariable():
                            x += 1
                            param_as_species += 1
                            out_file.write(
                                " <ic" + repr(x) + "> constant " + repr(list_of_parameters[k].getValue()) + " </ic" +
                                repr(x) + ">\n")
                            sum_file.write("S" + repr(x) + ":\t" + list_of_parameters[k].getId() + "\tparameter" +
                                           repr(k + 1 - comp) + "\t(" + repr(list_of_parameters[k].getValue()) +
                                           ") (parameter included in a rate rule and therefore treated as species)\n")

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
                out_file.write("<parameter%s realName='%s'>" % (counter, parameter_id[k]))
                out_file.write(" constant ")
                out_file.write(repr(parameter[k]) + " </parameter" + repr(counter) + ">\n")

        sum_file.write("\n############################################################\n\n")

        out_file.write("</parameters>\n")
        out_file.write("</model" + repr(i + 1) + ">\n\n")

    out_file.write("</models>\n\n")
    out_file.write("</input>\n\n")

    out_file.close()
    sum_file.close()
