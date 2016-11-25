"""
A set of functions for data fitting.

GML November 2016
"""

import inspect
import numpy as np
import scipy as sc
import scipy.optimize as optim


def get_errors(cov):
    """
    Find errors from covariance matrix

    Parameters
    ----------
    cov : np.ndarray
        Covariance matrix of the fit parameters.

    Returns
    -------
    err : 1-dim np.ndarray
        Errors asociated to the fit parameters.
    """
    return np.sqrt(np.diag(cov))


def get_from_name(name, glob=globals(), loc={}):
    """
    Checks whether the variable exists in current scope.

    Parameters
    ----------
    name : string
        Function name
    glob : dictionary
        Holds global variables
    loc : dictionary
        Holds local variables

    Returns
    -------
    obj : object
        Object pointed by variable *name*. Raises an error if not found.
    """
    if name not in glob and name not in loc:
        raise KeyError("Function {} is not defined".format(name))
    return glob[name] if name in glob else loc[name]


def gauss(x, amp, mu, sigma):
    if sigma<0:
        return np.inf
    return amp/(2*np.pi)**.5/sigma * np.exp(-0.5*(x-mu)**2./sigma**2.)


def polynom(x, *coeffs):
    return coeffs[0] + x * polynom(x, *coeffs[1:]) if len(coeffs) else 0.


def expo(x, const, mean):
    return const * np.exp(x/mean)


def power(x, const, pow_):
    return const * np.power(x, pow_)


def build_gauss(amp, mu, sigma):
    return lambda x: gauss(x, amp, mu, sigma)


def build_polynom(*coeffs):
    return lambda x: polynom(x, *coeffs)


def build_expo(const, mean):
    return lambda x: expo(x, const, mean)


def build_power(const, pow_):
    return lambda x: power(x, const, pow_)


def build(func, *args):
    f = get_from_name(func)
    return lambda x: f(x, *args)


def fit_gauss(x, y, seed=()):
    return sc.optimize.curve_fit(gauss, x, y, seed)


def fit_polynom(x, y, seed=()):
    return optim.curve_fit(polynom, x, y, seed)


def fit_expo(x, y, seed=()):
    return optim.curve_fit(expo, x, y, seed)


def fit_power(x, y, seed=()):
    return optim.curve_fit(power, x, y, seed)


def fit_simple(func, x, y, seed=()):
    """
    Fit x, y data to an specific already defined python function.

    Parameters
    ----------
    func : string
        Function name
    x, y : iterables
        Data sets to be fitted.
    seed : sequence
        Initial estimation of the fit parameters. Either all or none of them
        must be given.

    Returns
    -------
    fitted_fun : callable
        Callable function of 1 variable with the resulting parameters.
    vals : np.ndarray
        Optimal values for the parameters so that it minimizes the sum of the
        squares of func(x)-y
    err : 1-dim np.ndarray
        The estimated errors of vals.
    """
    fit_fun = get_from_name(func)
    vals, cov = optim.curve_fit(fit_fun, x, y, seed)
    vals, cov = fit_fun(x, y, seed)
    f = build(func, *vals)
    return f, vals, get_errors(cov)


def fit(func, x, y, seed=(), **kwargs):
    """
    Fit x, y data to a generic relation of already defined python functions.

    Parameters
    ----------
    func : string
        Generic relation to already defined functions. Because polynom is
        degree-free, it must be the last one in func.
        e.g. gauss + power / polynom.
    x, y : iterables
        Data sets to be fitted.
    seed : sequence
        Initial estimation of the fit parameters. Either all or none of them
        must be given.

    Returns
    -------
    fitted_fun : callable
        Callable function of 1 variable with the resulting parameters.
    vals : np.ndarray
        Optimal values for the parameters so that it minimizes the sum of the
        squares of func(x)-y
    err : 1-dim np.ndarray
        The estimated errors of vals.
    """
    str_fun = "lambda x, *args: "
    start, end = 0, -1
    local = locals()
    for token in func.split(" "):
        if token in ["+", "-", "*", "/", "**"]:
            str_fun += " {} ".format(token)
        else:
            f = get_from_name(token, globals(), locals())
            end = start + len(inspect.getargspec(f).args) - 1
            if start == end:
                end = ""
            str_fun += token + "(x, *args[{}:{}])".format(start, end)
        start = end
    exec("local['fit_fun'] = {}".format(str_fun)) in globals(), locals()
    fit_fun = local["fit_fun"]
    vals, cov = optim.curve_fit(fit_fun, x, y, seed, **kwargs)
    return lambda x: fit_fun(x, *vals), vals, get_errors(cov)
