from enum import Enum


class KernelType(Enum):
    component_wise_uniform = 1
    component_wise_normal = 2
    multivariate_normal = 3
    multivariate_normal_nn = 4
    multivariate_normal_ocm = 5
