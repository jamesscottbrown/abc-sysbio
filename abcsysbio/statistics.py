# statistical functions

from numpy import *
from numpy import random as rnd
from numpy import linalg as la
import scipy
import math
import scipy.stats.mvn


def w_choice(weight):
    """
    Sample from the categorical distribution with probabilities given by weight.

    Parameters
    ----------
    weight : list of probability for each category

    Returns
    -------

    """
    n = rnd.random_sample()
    for i in range(len(weight)):
        if n < weight[i]:
            return i
        n = n - weight[i]
    return len(weight) - 1


def get_pdf_uniform(min_val, max_val, x):
    """
    Evaluate the P(x) for x ~ U(min_val, max_val)

    Parameters
    ----------
    x : value at which to evaluate the p.d.f
    max_val : max value of uniform distribution
    min_val : min value of uniform distribution
    """
    if (x > max_val) or (x < min_val):
        return 0.0
    else:
        return 1 / (max_val - min_val)


def get_pdf_gauss(m, scale, x):
    """
    Evaluate P(x) for Gaussian distribution, x ~ N(m, scale^2)

    Parameters
    ----------
    x : value at which to evaluate the p.d.f
    scale : standard deviation of the Gaussian
    m : mean of the Gaussian
    """
    x = exp(-0.5 * (x - m) * (x - m) / (scale * scale))
    x = x / (scale * sqrt(2 * pi))
    return x


def get_pdf_lognormal(m, sigma, x):
    """
    Evaluate P(x) for lognormal distribution, x ~ ln N(m, sigma^2)

    N.B. If X is log-normal, then Y = ln(X) is normally distrbuted.

    Parameters
    ----------
    x : value at which to evaluate the p.d.f
    sigma : standard deviation of the associated normal
    m : mean of the associated normal
    """
    x = exp(-0.5 * (log(x) - m) * (log(x) - m) / (sigma * sigma))
    x = x / (x * sigma * sqrt(2 * pi))
    return x


# compute the pdf of a multinormal distribution
def get_pdf_multinormal(x, covariances, m):
    """
    Evaluate P(x) for multivariate normal distribution, x ~ N(means, covariances)

    Parameters
    ----------
    x : value at which to evaluate the p.d.f
    covariances : covariance matrix
    m : mean vector

    Returns
    -------

    """
    a = 0
    k = len(x)
    inv = la.inv(covariances)
    for i in range(k):
        for j in range(k):
            a += inv[i, j] * (x[i] - m[i]) * (x[j] - m[j])
    det = la.det(covariances)
    return exp(-1.0 * a / 2.0) / (sqrt((2 * pi) ** k * det))


def wtvar(x, weights, method="R"):
    """
    Compute the weighted variances for measurements x.

    This function is based on http://adorio-research.org/wordpress/?p=259

    Parameters
    ----------
    x : measurements
    weights : weigths for each measurement
    method : 'R' (Default) or 'nist'

    Returns
    -------
    the weighted variance of the measurements
    """
    sum_w = sum(weights)
    if method == "nist":
        x_bar_wt = sum([weights * x for weights, x in zip(weights, x)]) / sum_w  # fixed.2009.03.07, divisor added.
        np = sum([1 if (weights != 0) else 0 for weights in weights])
        d = sum_w * (np - 1.0) / np
        return sum([weights * (x - x_bar_wt) ** 2 for weights, x in zip(weights, x)]) / d
    else:
        sum_w2 = sum([weights ** 2 for weights in weights])
        x_bar_wt = sum([(weights * x) for (weights, x) in zip(weights, x)]) / sum_w
        return sum([(weights * (x - x_bar_wt) ** 2) for (weights, x) in zip(weights, x)]) * sum_w / (
            sum_w ** 2 - sum_w2)


def mvnd_gen(m, c):
    """
    Draw a sample from a multivariate normal distribution.

    Parameters
    ----------
    m :  mean vector
    c : covariance

    Returns
    -------
    a sample from the distribution
    """
    a = list(rnd.normal(0, 1, len(m)))
    lambdas, vect = la.eig(c)
    print lambdas
    tmp = mat(vect) * mat(diag(sqrt(lambdas))) * transpose(mat(a))
    res = list()
    for i in range(len(m)):
        res.append(m[i] + tmp[i, 0])
    return res


def mvstdnormcdf(lower, upper, corr_coef, **kwds):
    """
    standardized multivariate normal cumulative distribution function

    This is a wrapper for scipy.stats.kde.mvn.mvndst which calculates
    a rectangular integral over a standardized multivariate normal
    distribution.

    This function assumes standardized scale, that is the variance in each dimension
    is one, but correlation can be arbitrary, covariance = correlation matrix

    From statsmodels.sandbox.distributions.extras

    Parameters
    ----------
    lower, upper : array_like, 1d
       lower and upper integration limits with length equal to the number
       of dimensions of the multivariate normal distribution. It can contain
       -np.inf or np.inf for open integration intervals
    corrcoef : float or array_like
       specifies correlation matrix in one of three ways, see notes
    optional keyword parameters to influence integration
        * maxpts : int, maximum number of function values allowed. This
             parameter can be used to limit the time. A sensible
             strategy is to start with `maxpts` = 1000*N, and then
             increase `maxpts` if ERROR is too large.
        * abseps : float absolute error tolerance.
        * releps : float relative error tolerance.

    Returns
    -------
    cdfvalue : float
        value of the integral


    Notes
    -----
    The correlation matrix corrcoef can be given in 3 different ways
    If the multivariate normal is two-dimensional than only the
    correlation coefficient needs to be provided.
    For general dimension the correlation matrix can be provided either
    as a one-dimensional array of the upper triangular correlation
    coefficients stacked by rows, or as full square correlation matrix

    See Also
    --------
    mvnormcdf : cdf of multivariate normal distribution without
        standardization

    Examples
    --------

    >> print mvstdnormcdf([-np.inf, -np.inf], [0.0, np.inf], 0.5)
    0.5
    >> corr = [[1.0, 0, 0.5], [0, 1, 0], [0.5, 0,1]]
    >> print mvstdnormcdf([-np.inf, -np.inf, -100.0], [0.0, 0.0, 0.0], corr, abseps=1e-6)
    0.166666399198
    >> print mvstdnormcdf([-np.inf, -np.inf, -100.0],[0.0, 0.0, 0.0], corr, abseps=1e-8)
    something wrong completion with ERROR > EPS and MAXPTS function values used;
                        increase MAXPTS to decrease ERROR; 1.048330348e-006
    0.166666546218
    >> print mvstdnormcdf([-np.inf, -np.inf, -100.0], [0.0, 0.0, 0.0], corr, maxpts=100000, abseps=1e-8)
    0.166666588293
    """

    n = len(lower)

    lower = array(lower)
    upper = array(upper)
    corr_coef = array(corr_coef)

    correl = zeros(n * (n - 1) / 2.0)  # dtype necessary?

    if (lower.ndim != 1) or (upper.ndim != 1):
        raise ValueError('can handle only 1D bounds')

    if len(upper) != n:
        raise ValueError('bounds have different lengths')

    if n == 2 and corr_coef.size == 1:
        correl = corr_coef
    elif corr_coef.ndim == 1 and len(corr_coef) == n * (n - 1) / 2.0:
        correl = corr_coef
    elif corr_coef.shape == (n, n):
        for ii in range(n):
            for jj in range(ii):
                correl[jj + ((ii - 2) * (ii - 1)) / 2] = corr_coef[ii, jj]
    else:
        raise ValueError('corrcoef has incorrect dimension')

    if 'maxpts' not in kwds:
        if n > 2:
            kwds['maxpts'] = 10000 * n

    lowinf = isneginf(lower)
    uppinf = isposinf(upper)
    infin = 2.0 * ones(n)

    putmask(infin, lowinf, 0)  # infin.putmask(0,lowinf)
    putmask(infin, uppinf, 1)  # infin.putmask(1,uppinf)
    # this has to be last
    putmask(infin, lowinf * uppinf, -1)
    error, cdfvalue, inform = scipy.stats.mvn.mvndst(lower, upper, infin, correl, **kwds)
    if inform:
        print 'something wrong', inform, error, cdfvalue
    return cdfvalue


def mvnormcdf(lower, upper, mu, c, **kwds):
    """
    multivariate normal cumulative distribution function

    This is a wrapper for scipy.stats.kde.mvn.mvndst which calculates
    a rectangular integral over a multivariate normal distribution.

    From statsmodels.sandbox.distributions.extras

    Parameters
    ----------
    lower, upper : array_like, 1d
       lower and upper integration limits with length equal to the number
       of dimensions of the multivariate normal distribution. It can contain
       -np.inf or np.inf for open integration intervals
    mu : array_lik, 1d
       list or array of means
    c : array_like, 2d
       specifies covariance matrix
    optional keyword parameters to influence integration
        * maxpts : int, maximum number of function values allowed. This
             parameter can be used to limit the time. A sensible
             strategy is to start with `maxpts` = 1000*N, and then
             increase `maxpts` if ERROR is too large.
        * abseps : float absolute error tolerance.
        * releps : float relative error tolerance.

    Returns
    -------
    cdfvalue : float
        value of the integral


    Notes
    -----
    This function normalizes the location and scale of the multivariate
    normal distribution and then uses `mvstdnormcdf` to call the integration.

    See Also
    --------
    mvstdnormcdf : location and scale standardized multivariate normal cdf
    """

    lower = array(lower)
    upper = array(upper)
    c = array(c)
    stdev = sqrt(diag(c))
    lower = (lower - mu) / stdev
    upper = (upper - mu) / stdev
    divrow = atleast_2d(stdev)
    corr = c / divrow / divrow.T

    return mvstdnormcdf(lower, upper, corr, **kwds)


def k_nearest_neighbours(ind, s, k):
    """
    Compute the k nearest neighbors of a point inside a set S of points using the Euclidian distance.

    Parameters
    ----------
    ind : the index of the point within s whose nearest-neighbours are wanted
    s : positions of points, s[param, point]
    k : the number of nearest-neighbours to identify

    Returns
    -------

    """
    n = len(s[0])
    ss = array(s)
    aa = zeros((len(s), n))

    for param in range(len(s)):
        aa[param, :] = s[param][ind] * ones((1, n))

    dist = sum((aa - ss) ** 2, 2)
    m = max(dist)

    k_min = list()
    for i in range(min(k, n)):
        im = argmin(dist)
        k_min.append(im)
        dist[im] = Inf
    return k_min


def compute_cov(x, weights):
    """
    Compute the weighted covariance matrix for a set of measurements, by first calculating the weighted mean.

    Parameters
    ----------
    x : measurements
    weights : weights

    Returns
    -------

    """
    num_dimensions = len(x)
    num_samples = len(x[0])

    # Calculate weighted mean of the samples
    m = list()
    for d in range(num_dimensions):
        m.append(0)
        for sample in range(num_samples):
            m[d] += weights[sample] * x[d][sample]
        m[d] = m[d] / sum(weights)

    return compute_optcovmat(x, weights, m)


def compute_optcovmat(x, weights, m):
    """
    Compute the weighted covariance matrix for a set of measurements, given means.

    Parameters
    ----------
    x : measurements
    weights : weights
    m : mean

    Returns
    -------

    """
    num_dimensions = len(x)
    num_samples = len(x[0])
    c = zeros([num_dimensions, num_dimensions], float)

    # Fill in upper-half of c
    for sample in range(num_samples):
        for d1 in range(num_dimensions):
            for d2 in range(d1):
                c[d1, d2] += weights[sample] * (x[d1][sample] - m[d1]) * (x[d2][sample] - m[d2])

    # Fill in lower-half by symmetry
    c = c + transpose(c)

    # Fill in diagonal
    for sample in range(num_samples):
        for d1 in range(num_dimensions):
            c[d1, d1] += weights[sample] * (x[d1][sample] - m[d1]) * (x[d1][sample] - m[d1])

    # Divide every element by the total weight
    for d1 in range(num_dimensions):
        for d2 in range(num_dimensions):
            c[d1, d2] = c[d2, d2] / sum(weights)
    return c
