import numpy
import os
import re
from ctypes import *
from Prior import *


# compilation step to create share object for a correct solver and model
def compile(name, integration):
    integ = integration + "Solver"
    libname = "./lib" + name + ".so.1.0"

    abc_gsl_lib = os.getenv("GSL_LIB")
    abc_gsl_inc = os.getenv("GSL_INC")

    if abc_gsl_lib is None:
        abc_gsl_lib = "/usr/local/lib"
    if abc_gsl_inc is None:
        abc_gsl_inc = "/usr/local/include"

    abc_nm_lib = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'src/newmat11/')
    abc_nm_inc = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'src/newmat11/')
    abc_src_dir = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'src/')

    command = "make -f " + abc_src_dir + "makefile --quiet "
    command = command + " MODEL=" + name + " SOLVER=" + integ + " LIBNAME=" + libname + " "
    command = command + "ABC_GSL_LIB=" + abc_gsl_lib + " "
    command = command + "ABC_GSL_INC=" + abc_gsl_inc + " "
    command = command + "ABC_NM_LIB=" + abc_nm_lib + " "
    command = command + "ABC_NM_INC=" + abc_nm_inc + " "
    command = command + "ABC_SRC_DIR=" + abc_src_dir + " "

    os.system(command)
    return CDLL(libname)


class Model:
    # instantiation
    def __init__(self, name, nspecies, nparameters, prior, x0prior, source, integration, fit, dt, beta, initstep,
                 relative_error, absolute_error, logp):
        gil = re.compile('Gillespie', re.I)
        ode = re.compile('ODE', re.I)
        euler = re.compile('Euler', re.I)
        heun = re.compile('Heun', re.I)
        milstein = re.compile('Milstein', re.I)

        solver_name = ""
        if gil.search(integration):
            solver_name = 'Gillespie'
        elif ode.search(integration):
            solver_name = 'ODE'
        elif euler.search(integration):
            solver_name = 'EulerSDE'
        elif heun.search(integration):
            solver_name = 'HeunSDE'
        elif milstein.search(integration):
            solver_name = 'MilsteinSDE'
        else:
            print "C model : unrecognised integrator : ", solver_name

        self.name = name
        self.nspecies = nspecies

        # combine the parameters with the species
        self.kparameters = nparameters
        self.nparameters = nparameters + nspecies

        self.prior = prior
        self.prior.extend(x0prior)

        self.source = source
        self.integration = solver_name
        self.fit = fit
        self.dt = dt
        self.beta = beta
        self.initstep = initstep
        self.relativeError = relative_error
        self.absoluteError = absolute_error
        self.logp = logp

        self.lib = compile(self.name, self.integration)

        if self.integration == 'ODE':
            self.simulate = self.simulate_ode
        if self.integration == 'EulerSDE':
            self.simulate = self.simulate_sde
        if self.integration == 'MilsteinSDE':
            self.simulate = self.simulate_sde
        if self.integration == 'HeunSDE':
            self.simulate = self.simulate_sde
        elif self.integration == 'Gillespie':
            self.simulate = self.simulate_gillespie

    def simulate_ode(self, p, t, n, beta):
        # must return a structure of the form [n][ntimepoints][nspecies]
        # n is the number of models
        ntimepoints = len(t)
        ret = numpy.zeros([n, beta, ntimepoints, self.nspecies])
        # set output from cpp file to python
        output_arr_type = beta * (self.nspecies + 1) * ntimepoints * c_double
        output = output_arr_type()

        # set timepoints; ctypes
        tim_arr_type = ntimepoints * c_double
        ctime = tim_arr_type()
        for i in range(ntimepoints):
            ctime[i] = t[i]

        # set other ctypes
        c_beta = c_int(beta)
        c_ntimepoints = c_int(ntimepoints)
        c_nparameters = c_int(self.nparameters)
        c_nspecies = c_int(self.nspecies)
        c_initstep = c_double(self.initstep)
        c_absolute_error = c_double(self.absoluteError)
        c_relative_error = c_double(self.relativeError)
        for ni in range(n):
            # set parameters; ctypes
            par_arr_type = self.nparameters * c_double
            cparam = par_arr_type()
            for i in range(self.kparameters):
                if not self.logp:
                    cparam[i] = p[ni][i]
                else:
                    cparam[i] = numpy.power(10, p[ni][i])

            # set initial values; ctypes
            init_arr_type = self.nspecies * c_double
            cinit = init_arr_type()
            j = 0
            for i in range(self.kparameters, self.nparameters):
                cinit[j] = p[ni][i]
                j += 1

            dat = self.lib.MainC(byref(cinit), byref(cparam), c_beta, byref(ctime), c_ntimepoints, c_nparameters,
                                 c_nspecies, c_initstep, c_absolute_error, c_relative_error, byref(output))

            count = 0
            for j in range(beta):
                for h in range(self.nspecies):
                    for k in range(ntimepoints):
                        ret[ni][j][k][h] = output[count]
                        count += 1
        return ret

    def simulate_sde(self, p, t, n, beta):
        # must return a structure of the form [n][ntimepoints][nspecies]
        # where p is an 2D array of parameters, t is a list of timepoints, n is the number of models and beta is number of model 
        ntimepoints = len(t)
        ret = numpy.zeros([n, beta, ntimepoints, self.nspecies])
        # set output from cpp file to python
        output_arr_type = beta * (self.nspecies + 1) * ntimepoints * c_double
        output = output_arr_type()

        # set timepoints; ctypes
        tim_arr_type = ntimepoints * c_double
        c_time = tim_arr_type()
        for i in range(ntimepoints):
            c_time[i] = t[i]
        # set other ctypes
        c_beta = c_int(beta)
        c_ntimepoints = c_int(ntimepoints)
        c_dt = c_double(self.dt)
        c_nparameters = c_int(self.nparameters)
        c_npsecies = c_int(self.nspecies)

        for ni in range(n):
            # set parameters; ctypes
            par_arr_type = self.nparameters * c_double
            cparam = par_arr_type()
            for i in range(self.kparameters):
                if not self.logp:
                    cparam[i] = p[ni][i]
                else:
                    cparam[i] = numpy.power(10, p[ni][i])

            # set initial values; ctypes
            init_arr_type = self.nspecies * c_double
            cinit = init_arr_type()
            j = 0
            for i in range(self.kparameters, self.nparameters):
                cinit[j] = p[ni][i]
                j += 1

            # double* initialValues, double* parameters, int beta, double* timepoints, int ntimepoints, double dt, NPARAMETERS, NSPECIES
            dat = self.lib.MainC(byref(cinit), byref(cparam), c_beta, byref(c_time), c_dt, c_ntimepoints, c_nparameters,
                                 c_npsecies, byref(output))

            count = 0
            for j in range(beta):
                for h in range(self.nspecies):
                    for k in range(ntimepoints):
                        ret[ni][j][k][h] = output[count]
                        count += 1
        return ret

    def simulate_gillespie(self, p, t, n, beta):
        # must return a structure of the form [n][ntimepoints][nspecies]
        ntimepoints = len(t)
        ret = numpy.zeros([n, beta, ntimepoints, self.nspecies])

        # set output from cpp file to python
        output_arr_type = beta * (self.nspecies + 1) * ntimepoints * c_double
        output = output_arr_type()

        # set timepoints; ctypes
        tim_arr_type = ntimepoints * c_double
        ctime = tim_arr_type()
        for i in range(ntimepoints):
            ctime[i] = t[i]

        # set other ctypes
        c_beta = c_int(beta)
        c_ntimepoints = c_int(ntimepoints)
        c_nparameters = c_int(self.nparameters)
        c_nspecies = c_int(self.nspecies)

        for ni in range(n):
            # set parameters; ctypes
            par_arr_type = self.nparameters * c_double
            cparam = par_arr_type()
            for i in range(self.kparameters):
                if not self.logp:
                    cparam[i] = p[ni][i]
                else:
                    cparam[i] = numpy.power(10, p[ni][i])

            # set initial values; ctypes
            init_arr_type = self.nspecies * c_double
            cinit = init_arr_type()
            j = 0
            for i in range(self.kparameters, self.nparameters):
                cinit[j] = p[ni][i]
                j += 1

            self.lib.MainC(byref(cinit), byref(cparam), c_beta, byref(ctime), c_ntimepoints, c_nparameters, c_nspecies,
                           byref(output))

            count = 0
            for j in range(beta):
                for h in range(self.nspecies):
                    for k in range(ntimepoints):
                        ret[ni][j][k][h] = output[count]
                        count += 1
        return ret
