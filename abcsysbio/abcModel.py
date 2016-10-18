class AbcModel:

    """ABCModel class.
    The abcsmc class requires model objects; if cudasim is used, these will be the appopriate simulator class.
    This provides as alternative for when cudasim is not being used.
    """

    def __init__(self,
                 name,  # just a string
                 simulationFn,
                 distanceFn,
                 prior,  # array of length nparameters
                 nparameters,  # including initial conditions, etc
                 parameterNames=None,
                 ):
        self.name = name
        self.simulationFn = simulationFn
        self.distanceFn   = distanceFn
        self.nparameters  = nparameters
        self.prior        = prior  # this is stupid, should be an array, so really should be called priors!
        if parameterNames is None:
            parameterNames = ['P%d' % x for x in range(self.nparameters)]
        self.parameterNames = parameterNames



    # TODO: remove some of these args?
    def simulate(self, params, timepoints, num_simulations, beta):
        simulatedData = apply(self.simulationFn, (params,))
        return simulatedData

    def distance(self, simulatedData, targetData, params, _unusedModel):
        d = apply(self.distanceFn, (simulatedData, targetData, params, self))
        return d