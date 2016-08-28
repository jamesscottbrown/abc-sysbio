from collections import namedtuple

Prior = namedtuple('prior', ['type', 'value', 'mean', 'variance', 'lower_bound', 'upper_bound', 'mu', 'sigma'])
Prior.__new__.__defaults__ = (None,) * len(Prior._fields)
