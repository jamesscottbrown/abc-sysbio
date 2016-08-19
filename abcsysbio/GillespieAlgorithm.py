import numpy

#Gillespie solver
#Needs to take as parameters initial time (assume zero), maximum time,
#values for parameters (passed to the function Hazards) and initial values of
#species (passed to the first row of the concentrations[][] matrix during the first iteration)

def select_hazard(total_hazard,list_of_hazards):
    """chose a hazard from the list by by hazard weights.

    ***** args *****

    total hazard:

            Sum of hazards

    list of hazards:

            List of hazards, one for each reaction

    """

    hazard_chooser = numpy.random.random_sample()
    for i in range(0,len(list_of_hazards)):
        x=list_of_hazards[i]/total_hazard #normal hazards
        if hazard_chooser<x:
            break
        hazard_chooser=hazard_chooser-x
    return i




def GillespieInt(func, initValues, parameters, outputtimes):
        """Simulate the function func from the initial values initValues using the parameters given by parameters.
        Return values at the times in outputtimes.

        ***** args *****
        
        func:

                a function to integrate, generated by parseInfo.py from an SBML model.
                
        initValues:

                a list of floats representing the initial values of the species whose trajectories are described by func.
                
        parameters:
        
                a tuple of parameters to pass to func
                
        outputtimes:
        
                a list of times at which to return data points.

        """
        
        #Set initial time
        time = 0 #outputtimes[0] # change this to t=0
        #set initial concentrations
        concentrations = numpy.zeros([len(outputtimes),len(initValues)])
        current_concentrations = tuple(initValues)
        current_concentrations, parameters=func.events(current_concentrations, parameters,time)
        current_concentrations, parameters=func.rules(current_concentrations, parameters,time)
        
        concentrations[0] = current_concentrations
        #Initialise outputtimes
        counter = 0 # 1 # change this to 0
        #Get the dictionary for choosing a function        
        switch = func.Switch()

        flag = True 

        while flag:
                current_concentrations, parameters=func.events(current_concentrations, parameters,time)
                current_concentrations, parameters=func.rules(current_concentrations, parameters,time)
                list_of_hazards = func.Hazards(current_concentrations, parameters)
                total_hazard = sum(list_of_hazards)
                #Check that total hazard is not less than zero.

                if total_hazard <= 0.0:
                        time = outputtimes[-1]
                else:
                        URN = numpy.random.uniform()
                        new_time = (-1.0/total_hazard) * numpy.log(URN)
                        time += new_time

                while time >= outputtimes[counter]:
                        concentrations[counter] = current_concentrations
                        counter += 1
                        if counter >= len(outputtimes):
                                flag = False
                                break

                if total_hazard <= 0.0:
                    pass
                else:
                    hazard_chosen=select_hazard(total_hazard,list_of_hazards)

                    try:
                        current_concentrations = tuple(switch[hazard_chosen](tuple(current_concentrations))) # Switch based on value of hazard_chosen
                    except KeyError:
                        switch["default"]()
                                                                  
        return concentrations


def GillespieInt_onestep(func, current_concentrations, t1, t2, parameters, xmax):
    #print current_concentrations, t1,  t2, parameters

    time = t1
    switch = func.Switch()

    while 1:
        current_concentrations, parameters=func.events(current_concentrations, parameters, time)
        current_concentrations, parameters=func.rules(current_concentrations, parameters, time)

        # print current_concentrations, time

        list_of_hazards = func.Hazards(current_concentrations, parameters)
        total_hazard = sum(list_of_hazards)

        if sum( current_concentrations ) > xmax :
            return[ current_concentrations, time, False ]

        if total_hazard <= 0.0:
            return[ current_concentrations, time, False ]

        URN = numpy.random.uniform()
        new_time = (-1.0/total_hazard) * numpy.log(URN)
        time += new_time

        hazard_chosen=select_hazard(total_hazard,list_of_hazards)

        try:
            current_concentrations = tuple(switch[hazard_chosen](tuple(current_concentrations))) # Switch based on value of hazard_chosen
        except KeyError:
            switch["default"]()

        if time >= t2:
            return [ current_concentrations, time, True ]
