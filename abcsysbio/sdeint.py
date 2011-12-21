from numpy import *
import math


def sdeint(func,InitValues,parameter, timepoints,dt=0.01):
    
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
    
    maxT=timepoints[len(timepoints)-1]
    times = arange(0.0, maxT+dt,dt)    
    length = len(times)   
    
    dim=len(InitValues)
    solutions=zeros([length,dim])
    solutions_out=zeros([len(timepoints),dim])
    #InitValues, parameter = func.rules(InitValues, parameter, times[0])
    InitValues, parameter = func.rules(InitValues, parameter,0) # change this to t=0
    solutions[0]=InitValues
    #solutions_out[0]=InitValues
    
    n=0
    for i in range(1,length):
        new,W=func.modelfunction(solutions[i-1],dt,parameter,time=times[i])
        
        for k in range(0,dim):
            solutions[i][k]=solutions[i-1][k]+new[k]*dt+W[k]
            
        solutions[i],parameter=func.rules(solutions[i],parameter,times[i])
            
        if(solutions[i][k]<0.0):
            solutions[i][k]=0.0
            #print "\nSDE simulation failed. Try smaller timestep.\n"
            return solutions_out
        
        if(n>=len(timepoints)):
            return zeros([len(timepoints),dim])
        if((timepoints[n]-times[i])<0.000000001):
            solutions_out[n]=solutions[i]
            n=n+1

    return solutions_out

def sdeint_onestep(func, current_concentrations, t1, t2, parameters, dt=0.01):

    dim = len(current_concentrations)
    X = zeros([1,dim])

    time = t1
    next_time = min(time+dt,t2)

    X[0], parameters = func.rules(current_concentrations, parameters, t=time)

    while 1:
        dX, W = func.modelfunction(X[0], next_time - time, parameters, time=next_time)

        #print 's', current_concentrations
        #print 'new', dX[0]
        #print 'W', W[0]
        
        for k in range(0,dim):
            X[0][k] += dX[k]*(next_time-time) + W[k]

        X[0], parameters = func.rules(X[0], parameters, next_time)
        
        for k in range(0,dim):
            if(X[0][k]<0.0):
                X[0][k]=0.0
                #print "\nSDE simulation failed. Try smaller timestep.\n"
                return [ X[0], next_time, False ]

        if ((t2 - next_time)<0.000000001):
            return [ X[0], next_time, True ]

        time = next_time
        next_time = min(time+dt,t2)

    return [ X[0], time, False ]
