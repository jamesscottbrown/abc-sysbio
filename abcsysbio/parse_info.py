# Algorithm information

import re, sys, numpy

from xml.dom import minidom
from KernelType import KernelType
from PriorType import PriorType
from Prior import Prior

# implemented priors
re_prior_const = re.compile('constant')
re_prior_uni = re.compile('uniform')
re_prior_normal = re.compile('normal')
re_prior_logn = re.compile('lognormal')
re_prior_catg = re.compile('categorical')

# implemented kernels
re_kernel_uniform = re.compile('uniform')
re_kernel_normal = re.compile('normal')
re_kernel_mvnormal = re.compile('multiVariateNormal')
re_kernel_mvnormalKN = re.compile('multiVariateNormalKNeigh')
re_kernel_mvnormalOCM = re.compile('multiVariateNormalOCM')

# True/False
re_true = re.compile('True')
re_none = re.compile('None')


def parse_required_single_value(node, tagname, message, cast):
    """
    Given a tag, try to find the child element with given tagname, and cast it contents to a given type.
    If this fails, exit.

    Parameters
    ----------
    node : xml.dom.minimdom.Document instance
    tagname : string containing name of tag (e.g. 'modelnumber')
    message : message to print on failure
    cast : type to cast value to

    Returns
    -------
    value of type cast

    """
    try:
        data = node.getElementsByTagName(tagname)[0].firstChild.data
        ret = cast(data)
    except (IndexError, ValueError):
        sys.exit(message)

    return ret


def parse_required_vector_value(node, tag_name, message, cast):
    """
     Given a tag, try to find the child element with given tagname, and split its contents into a list, and cast each
     element to a given type. If this fails, exit.

    Parameters
    ----------
     node : xml.dom.* object
     tag_name : string containing name of tag (e.g. 'modelnumber')
     message : message to print on failure
     cast : type to cast value to

    Returns
    -------
    list of values of type cast

    """
    try:
        data = node.getElementsByTagName(tag_name)[0].firstChild.data
        tmp = str(data).split()
        ret = [cast(i) for i in tmp]
    except (IndexError, ValueError):
        sys.exit(message)

    if len(ret) == 0:
        sys.exit(message)

    return ret


def process_prior(tmp, model_num):
    """
    Check that prior has been specified correctly.

    Parameters
    ----------
    tmp : list of strings specifying a prior (e.g. ['uniform', '0', '10'] or ['uniform', '1.0'])
    model_num : number of model (used in error messages)

    Returns
    -------
    A NamedTuple representing the prior.
    The "type" field represents the prior type, and determines which other fields exist.

    """

    if re_prior_const.match(tmp[0]):
        try:
            prior_params = Prior(type=PriorType.constant, value=float(tmp[1]))
        except ValueError:
            sys.exit("\nValue of the prior for model %s (counting from 1) has wrong format: %s" % (model_num, tmp[1]))

    elif re_prior_normal.match(tmp[0]):
        try:
            prior_params = Prior(type=PriorType.normal, mean=float(tmp[1]), variance=float(tmp[2]))
        except ValueError:
            sys.exit("\nValue of the prior for model %s (counting from 1) has wrong format: %s" % (model_num, tmp[1]))

    elif re_prior_uni.match(tmp[0]):
        try:
            prior_params = Prior(type=PriorType.uniform, lower_bound=float(tmp[1]), upper_bound=float(tmp[2]))
        except ValueError:
            sys.exit("\nValue of the prior for model %s (counting from 1) has wrong format: %s" % (model_num, tmp[1]))

    elif re_prior_logn.match(tmp[0]):
        try:
            prior_params = Prior(type=PriorType.lognormal, mu=float(tmp[1]), sigma=float(tmp[2]))
        except ValueError:
            sys.exit("\nValue of the prior for model %s (counting from 1) has wrong format: %s" % (model_num, tmp[1]))
    elif re_prior_catg.match(tmp[0]):
        try:
            prior_params = Prior(type=PriorType.lognormal, p=[float(x) for x in tmp[1:]])
        except ValueError:
            sys.exit("\nValue of the prior for model %s (counting from 1) has wrong format: %s" % (model_num, tmp[1]))

    else:
        sys.exit("\nSupplied parameter prior %s unsupported" % tmp[0])

    return prior_params


def parse_fitting_information(node):
    """
    Given an element whose child is a <fit> tag, construct the 'fit' strings used in HowToFitData.

    Parameters
    ----------
    node : some kind of xml.dom.* object

    Returns
    -------
    a list, each element of which is a string representing what function of the simulation results is intended to be
    fitted to the corresponding experimental measurement
    """

    fit_tag = node.getElementsByTagName('fit')[0]
    raw_fiting_strings = str(fit_tag.firstChild.data).split()

    if len(raw_fiting_strings) == 1 and re_none.match(raw_fiting_strings[0]):
        return None
    else:
        fitting_strings = []
        for i in raw_fiting_strings:
            # replace species with sample_points
            fitting_string = re.sub('species', 'sample_points', i)

            # find all instances of sample_points[ ] and extract the list of species numbers
            sp_strings = re.findall("sample_points([0-9]+)", fitting_string)
            sp_nums = [int(j) for j in sp_strings]
            sp_nums.sort()
            sp_nums.reverse()

            # loop over the species numbers and replace
            for n in sp_nums:
                fitting_string = re.sub('ts' + str(n), 'ts[:,' + str(n - 1) + ']', fitting_string)

            fitting_strings.append(fitting_string)

        return fitting_strings


class AlgorithmInfo:
    """
    A class to parse the user-provided input file and return all information required to run the abc-SMC algorithm.
    
    """

    def __init__(self, filename, mode):
        xmldoc = minidom.parse(filename)
        self.mode = mode
        # mode is 0  inference, 1 simulate, 2 design

        self.modelnumber = 0
        self.restart = False
        self.particles = 0
        self.beta = 0
        self.dt = 0
        self.epsilon = []
        self.final_epsilon = []
        self.alpha = 0.9
        self.times = []
        self.ntimes = 0
        self.data = []

        self.nmodels = 0
        self.nparameters = []
        self.nspecies = []
        self.name = []
        self.source = []
        self.type = []
        self.prior = []
        self.x0prior = []
        self.fit = []
        self.logp = []

        self.modelkernel = 0.7
        self.kernel = KernelType.component_wise_uniform
        self.modelprior = []
        self.rtol = 1e-5
        self.atol = 1e-5

        # Required arguments
        self.modelnumber = parse_required_single_value(xmldoc, "modelnumber",
                                                       "Please provide an integer value for <modelnumber>", int)

        self.particles = parse_required_single_value(xmldoc, "particles",
                                                     "Please provide an integer value for <particles>", int)

        self.beta = parse_required_single_value(xmldoc, "beta", "Please provide an integer value for <beta>", int)

        self.dt = parse_required_single_value(xmldoc, "dt", "Please provide an float value for <dt>", float)

        # Get epsilon values, unless in simulate mode
        if self.mode != 1:

            # automated epsilon takes priority
            autoepsilon_tags = xmldoc.getElementsByTagName('autoepsilon')
            if len(autoepsilon_tags) > 0:
                self.final_epsilon = parse_required_vector_value(autoepsilon_tags[0], "finalepsilon",
                                                                 "Please provide a whitespace separated list of " +
                                                                 "values for <autoepsilon><finalepsilon>", float)
                self.alpha = parse_required_single_value(autoepsilon_tags[0], "alpha",
                                                         "Please provide a float value for <autoepsilon><alpha>",
                                                         float)
            else:
                # do a first pass to find the number of epsilon values we need to store
                num_schedules = 0
                num_epsilons_per_schedule = 0
                epsilon_tag = xmldoc.getElementsByTagName('epsilon')[0]

                for e in epsilon_tag.childNodes:
                    if e.nodeType == e.ELEMENT_NODE:
                        num_schedules += 1
                        num_epsilons_per_schedule = len(str(e.firstChild.data).split())

                # do a second pass to store these values
                self.epsilon = numpy.zeros([num_schedules, num_epsilons_per_schedule])
                schedule = 0
                for e in epsilon_tag.childNodes:
                    if e.nodeType == e.ELEMENT_NODE:
                        tmp = str(e.firstChild.data).split()

                        for population in range(num_epsilons_per_schedule):
                            self.epsilon[schedule, population] = float(tmp[population])

                        schedule += 1

        # Get data attributes
        dataref = xmldoc.getElementsByTagName('data')[0]

        self.times = parse_required_vector_value(dataref, "times",
                                                 "<data><times> requires a whitespace separated list of values", float)
        self.ntimes = len(self.times)

        # variables
        if self.mode == 0:
            # first do a scan to get the number of timeseries
            num_vars = 0
            var_tag = dataref.getElementsByTagName('variables')[0]
            for v in var_tag.childNodes:
                if v.nodeType == v.ELEMENT_NODE:
                    num_vars += 1

            # create matrix and mask
            data_unmasked = numpy.zeros([self.ntimes, num_vars])
            data_mask = numpy.zeros([self.ntimes, num_vars], dtype=numpy.int32)
            num_vars = 0
            for v in var_tag.childNodes:
                if v.nodeType == v.ELEMENT_NODE:
                    tmp = str(v.firstChild.data).split()

                    for population in range(self.ntimes):
                        # Search for NA
                        if re.match("\s*NA\s*", tmp[population]) is not None:
                            data_mask[population, num_vars] = 1
                            tmp[population] = 0

                        data_unmasked[population, num_vars] = float(tmp[population])

                    num_vars += 1

            # create masked data
            self.data = numpy.ma.array(data_unmasked, mask=data_mask)

        # get model attributes
        model_tag = xmldoc.getElementsByTagName('models')[0]
        for m in model_tag.childNodes:
            if m.nodeType == m.ELEMENT_NODE:
                self.nmodels += 1
                self.prior.append([])
                self.x0prior.append([])

                self.name.append(str(m.getElementsByTagName('name')[0].firstChild.data).strip())
                self.source.append(str(m.getElementsByTagName('source')[0].firstChild.data).strip())
                self.type.append(str(m.getElementsByTagName('type')[0].firstChild.data).strip())

                self.fit.append(parse_fitting_information(m))

                try:
                    tmp = str(m.getElementsByTagName('logp')[0].firstChild.data).strip()
                    if re_true.match(tmp):
                        self.logp.append(True)
                    else:
                        self.logp.append(False)
                except IndexError:
                    self.logp.append(False)

                num_params = 0
                param_tag = m.getElementsByTagName('parameters')[0]
                for p in param_tag.childNodes:
                    if p.nodeType == p.ELEMENT_NODE:
                        num_params += 1
                        tmp = str(p.firstChild.data).split()
                        self.prior[self.nmodels - 1].append(process_prior(tmp, self.nmodels))

                num_initial_conditions = 0
                initial_tag = m.getElementsByTagName('initial')[0]
                for inn in initial_tag.childNodes:
                    if inn.nodeType == inn.ELEMENT_NODE:
                        num_initial_conditions += 1
                        tmp = str(inn.firstChild.data).split()
                        self.x0prior[self.nmodels - 1].append(process_prior(tmp, self.nmodels))

                if num_params == 0:
                    sys.exit("\nNo parameters specified in model %s" % self.name[self.nmodels - 1])
                if num_initial_conditions == 0:
                    sys.exit("\nNo initial conditions specified in model %s" % self.name[self.nmodels - 1])
                self.nparameters.append(num_params)
                self.nspecies.append(num_initial_conditions)

        if self.nmodels == 0:
            sys.exit("\nNo models specified")

        ##################################################
        # Optional arguments

        # get atol
        try:
            data = xmldoc.getElementsByTagName('atol')[0].firstChild.data
            self.atol = float(data)
        except (IndexError, ValueError):
            pass

        # get rtol
        try:
            data = xmldoc.getElementsByTagName('rtol')[0].firstChild.data
            self.rtol = float(data)
        except (IndexError, ValueError):
            pass

        # get restart
        try:
            tmp = str(xmldoc.getElementsByTagName('restart')[0].firstChild.data).strip()
            if re_true.match(tmp):
                self.restart = True
        except IndexError:
            pass

        # get model kernel
        try:
            data = xmldoc.getElementsByTagName('modelkernel')[0].firstChild.data
            try:
                self.modelkernel = float(data)
            except ValueError:
                print "\n#################\n<modelkernel> must be a float so I am going to ignore your argument"

            if self.modelkernel > 1.0:
                print "\n#################\n<modelkernel> must be <= 1.0  so I am going to ignore your argument"
                self.modelkernel = 0.7
        except IndexError:
            pass

        # get kernel
        try:
            data = str(xmldoc.getElementsByTagName('kernel')[0].firstChild.data).strip()
            if re_kernel_uniform.match(data):
                self.kernel = KernelType.component_wise_uniform
            elif re_kernel_normal.match(data):
                self.kernel = KernelType.component_wise_normal
            elif re_kernel_mvnormal.match(data):
                self.kernel = KernelType.multivariate_normal
            elif re_kernel_mvnormalKN.match(data):
                self.kernel = KernelType.multivariate_normal_nn
            elif re_kernel_mvnormalOCM.match(data):
                self.kernel = KernelType.multivariate_normal_ocm
            else:
                print "\n#################"
                print "<kernel> must be one of uniform, normal, multivariateNormal, multivariateNormalKNeigh or " + \
                      "multivariateNormalOCM  so I will ignore your argument"
        except IndexError:
            pass

        # get model priors
        self.modelprior = [1 / float(self.nmodels)] * self.nmodels
        try:
            data = xmldoc.getElementsByTagName("modelprior")[0].firstChild.data
            tmp = str(data).split()

            ret = []
            try:
                ret = [float(population) for population in tmp]
            except ValueError:
                print "\n#################\n<modelprior> must be a vector of floats so I will ignore your argument"

            if sum(ret) != 1.0 or len(ret) != self.nmodels:
                print "\n#################"
                print "<modelprior> must sum to one and be the same length as the number of models so I will ignore " \
                      "your argument"
            else:
                self.modelprior = ret[:]
        except IndexError:
            pass

    def print_info(self):
        """
        Pretty-print the attributes of this algorithm_info object.
        """

        print "\nALGORITHM INFO"
        print "modelnumber:", self.modelnumber
        print "restart:", self.restart
        print "particles:", self.particles
        print "beta:", self.beta
        print "dt:", self.dt
        if self.mode != 1:
            if len(self.final_epsilon) == 0:
                print "manual epsilon:"
                for i in range(self.epsilon.shape[0]):
                    print "\t",
                    for j in range(self.epsilon.shape[1]):
                        print "", self.epsilon[i, j],
                    print ""
            else:
                print "auto epsilon:"
                print "\t", self.final_epsilon
                print "\talpha:", self.alpha

            print "kernel:", self.kernel
            print "model kernel:", self.modelkernel
        print "model prior:", self.modelprior

        print "DATA:"
        print "\ttimes:", self.times
        if self.mode == 0:
            print "\tvars:"
            for i in range(len(self.data[0, :])):
                print "\t",
                for j in range(self.ntimes):
                    print "", self.data[j, i],
                print ""

        print "MODELS:", self.nmodels
        for i in range(self.nmodels):
            print "\t", "npar:", self.nparameters[i]
            print "\t", "nspecies:", self.nspecies[i]
            print "\t", "name:", self.name[i]
            print "\t", "source:", self.source[i]
            print "\t", "type:", self.type[i]
            print "\t", "fit:", self.fit[i]
            print "\t", "init:", self.x0prior[i]
            print "\t", "prior:", self.prior[i]
            print "\t", "logp:", self.logp[i]
            print "\n"
