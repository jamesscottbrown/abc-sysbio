import numpy as np


class EpsilonSchedule(object):

    """Class to construct a fixed/non-adaptive epsilon schedule"""

    def __init__(self, tol_type, tmin, tmax,nt):
        """Init innit.
        Input:
            tol_type: specify const, linear,exp or log; default is exp
            tmin: minimum threshold for metric
            tmax: maximum threshold for metric
            nt: number of iterations
        """
        self.tol_type = tol_type
        self.nt = nt
        self.tmin = tmin
        self.tmax = tmax
        self.tol = self.set_tolerance()

    def set_tolerance(self):
        """Method to set tolerance type either const, linear,exp or log."""
        if self.tol_type == "const":
            return self.const_tol()
        elif self.tol_type == "linear":
            return self.linear_tol()
        elif self.tol_type == "exp":
            return self.exp_tol()
        elif self.tol_type == "log":
            return self.log_tol()
        else:
            print "Specify either const, linear, exp or log for tolerance class"

    def linear_tol(self):
        """Linearly decreasing tolerance level."""
        return np.linspace(self.tmax, self.tmin, num=self.nt)

    def log_tol(self):
        """Log decreasing tolerance level."""
        return np.logspace(self.tmax, self.tmin, num=self.nt)

    def const_tol(self):
        """Constant tolerance level for every iteration."""
        return np.ones(self.nt)*self.tmin

    def exp_tol(self):
        """Exponentially decreasing tolerance level."""
        return np.logspace(np.log10(self.tmax), np.log10(self.tmin), num=self.nt)