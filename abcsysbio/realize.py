# combines simulations and distance calculation
# simulation stops when distance is increased
# much more efficient

import re, numpy

from GillespieAlgorithm import GillespieInt_onestep
from abcodesolve import abcodeint_onestep
from sdeint import sdeint_onestep
from evaluate import evaluateDistance

#############
def realize_particles( ModelName, selected_model, integrationType, initValues, timepoints, parameters,
                       epsilon, t, fit, data, distancefn, fitfn, dt, rtol, atol, large_dist, x_max ):

    o=re.compile('ODE')
    s=re.compile('SDE')
    g=re.compile('Gillespie')

    mtyp = 0
    if( o.match(integrationType[selected_model]) ):
        mtyp = 1
    elif( s.match(integrationType[selected_model]) ): 
        mtyp = 2
    elif( g.match(integrationType[selected_model]) ):
        mtyp = 3
    else:
        print "\n undefined integration type in realize_particles\n"
    
    mdl = __import__(ModelName[selected_model])

    #Set initial time
    time = timepoints[0]

    #set initial concentrations
    concentrations = numpy.zeros([len(timepoints),len(initValues)])
    current_concentrations = tuple(initValues)
    
    if(mtyp == 3) : current_concentrations, parameters=mdl.events(current_concentrations, parameters,time)
    current_concentrations, parameters=mdl.rules(current_concentrations, parameters,time)
        
    concentrations[0] = current_concentrations
    
    counter = 1
    distance = large_dist

    while 1:
        if(mtyp==1):
            current_concentrations, time, flag = abcodeint_onestep(mdl, current_concentrations, time, timepoints[counter], parameters, dt, rtol, atol )
        if(mtyp==2):
            current_concentrations, time, flag = sdeint_onestep(mdl, current_concentrations, time, timepoints[counter], parameters, dt )
        if(mtyp==3):
            current_concentrations, time, flag = GillespieInt_onestep(mdl, current_concentrations, time, timepoints[counter], parameters, x_max)

        #print "Main loop:", current_concentrations, time,  timepoints[counter], parameters

        if flag == False:
            distance = large_dist
            break

        concentrations[counter] = current_concentrations
        points = fitfn(fit,concentrations[0:counter+1])
        distance = distancefn(points,data[0:counter+1])

        dist=evaluateDistance(distance,epsilon,t)
        if dist == False:
            break

        counter = counter + 1
        if counter >= len(timepoints):
            break

    return distance
