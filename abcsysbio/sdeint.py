from numpy import *
import math


def sdeint(func, init_values, parameter, timepoints, dt=0.01):
    """
    ***** args *****
    
    func:        a python function that defines the stochastic system.
                    This function takes the species concentration, the systems 
                parameter, the actual integration time and the internal time step
                as arguments. It returns the derivatives of the species 
                concentrations and the noise terms for each species.
    
    InitValues:  a tuple of floats.
                    InitValues contains the initial concentrations for each species.
    
    parameter:   a tuple of floats.
                    It contains the values for all systems parameter.
    
    timepoints:  a tuple of floats.
                    This tuple contains all external time points. It has to be sorted 
                from low to high.
    
    ***** kwargs *****
    
    dt:          a float number.
                    dt describes the internal time step.
    
    
    """

    max_t = timepoints[len(timepoints) - 1]
    times = arange(0.0, max_t + dt, dt)
    length = len(times)

    dim = len(init_values)
    solutions = zeros([length, dim])
    solutions_out = zeros([len(timepoints), dim])
    init_values, parameter = func.rules(init_values, parameter, 0)  # change this to t=0
    solutions[0] = init_values

    n = 0
    for i in range(1, length):
        new, w = func.modelfunction(solutions[i - 1], dt, parameter, time=times[i])

        for k in range(0, dim):
            solutions[i][k] = solutions[i - 1][k] + new[k] * dt + w[k]

        solutions[i], parameter = func.rules(solutions[i], parameter, times[i])

        # if any concentration has gone negative, terminate the simulation
        for k in range(0, dim):
            if solutions[i][k] < 0.0:
                solutions[i][k] = 0.0
                return solutions_out

        if n >= len(timepoints):
            return zeros([len(timepoints), dim])
        if (timepoints[n] - times[i]) < 0.000000001:
            solutions_out[n] = solutions[i]
            n += 1

    return solutions_out


def sdeint_onestep(func, current_concentrations, t1, t2, parameters, dt=0.01):
    dim = len(current_concentrations)
    x = zeros([1, dim])

    time = t1
    next_time = min(time + dt, t2)

    x[0], parameters = func.rules(current_concentrations, parameters, t=time)

    while 1:
        dx, w = func.modelfunction(x[0], next_time - time, parameters, time=next_time)

        for k in range(0, dim):
            x[0][k] += dx[k] * (next_time - time) + w[k]

        x[0], parameters = func.rules(x[0], parameters, next_time)

        for k in range(0, dim):
            if x[0][k] < 0.0:
                x[0][k] = 0.0
                return [x[0], next_time, False]

        if (t2 - next_time) < 0.000000001:
            return [x[0], next_time, True]

        time = next_time
        next_time = min(time + dt, t2)

    return [x[0], time, False]
