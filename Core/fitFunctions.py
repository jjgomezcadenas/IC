"""
A set of functions for data fitting.

GML November 2016
"""

import inspect
import numpy as np
import scipy as sc
import scipy.optimize as optim


def in_range(data, minval=-np.inf, maxval=np.inf):
    """
    Find values in range.

    Parameters
    ---------
    data : np.ndarray
        Data set of arbitrary dimension.
    minval : int or float, optional
        Range minimum. Defaults to -inf.
    maxval : int or float, optional
        Range maximum. Defaults to +inf.

    Returns
    -------
    selection : np.ndarray
        Boolean array with the same dimension as the input. Contains True
        for those values of data in the input range and False for the others.
    """
    return (minval <= data) & (data <= maxval)


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


def profileX(xdata, ydata, nbins, xrange=None, yrange=None, drop_nan=True):
    """
    Compute the x-axis binned average of a dataset.

    Parameters
    ----------
    xdata, ydata : 1-dim np.ndarray
        x and y coordinates from a dataset.
    nbins : int
        Number of divisions in the x axis.
    xrange : tuple of ints/floats or None, optional
        Range over the x axis. Defaults to dataset extremes.
    yrange : tuple of ints/floats or None, optional
        Range over the y axis. Defaults to dataset extremes.
    drop_nan : bool, optional
        Exclude empty bins. Defaults to True.

    Returns
    -------
    x_out : 1-dim np.ndarray.
        Bin centers.
    y_out : 1-dim np.ndarray
        Data average for each bin.
    y_err : 1-dim np.ndarray
        Average error for each bin.
    """
    xmin, xmax = (np.min(xdata), np.max(xdata)) if xrange is None else xrange
    ymin, ymax = (np.min(ydata), np.max(ydata)) if yrange is None else yrange

    x_out = np.linspace(xmin, xmax, nbins+1)
    y_out = np.empty(nbins)
    y_err = np.empty(nbins)
    dx = np.diff(x_out)[0]

    selection = in_range(xdata, xmin, xmax) & in_range(ydata, ymin, ymax)
    xdata, ydata = xdata[selection], ydata[selection]
    for i in range(nbins):
        bin_data = np.extract(in_range(xdata, x_out[i], x_out[i+1]), ydata)
        y_out[i] = np.mean(bin_data)
        y_err[i] = np.std(bin_data) / bin_data.size**0.5
    x_out += dx / 2.
    x_out = x_out[:-1]
    if drop_nan:
        selection = ~(np.isnan(y_out) | np.isnan(y_err))
        x_out = x_out[selection]
        y_out = y_out[selection]
        y_err = y_err[selection]
    return x_out, y_out, y_err


def profileY(xdata, ydata, nbins, yrange=None, xrange=None, drop_nan=True):
    """
    Compute the y-axis binned average of a dataset.

    Parameters
    ----------
    xdata, ydata : 1-dim np.ndarray
        x and y coordinates from a dataset.
    nbins : int
        Number of divisions in the y axis.
    yrange : tuple of ints/floats or None, optional
        Range over the y axis. Defaults to dataset extremes.
    xrange : tuple of ints/floats or None, optional
        Range over the x axis. Defaults to dataset extremes.
    drop_nan : bool, optional
        Exclude empty bins. Defaults to True.

    Returns
    -------
    x_out : 1-dim np.ndarray.
        Bin centers.
    y_out : 1-dim np.ndarray
        Data average for each bin.
    y_err : 1-dim np.ndarray
        Average error for each bin.
    """
    xmin, xmax = (np.min(xdata), np.max(xdata)) if xrange is None else xrange
    ymin, ymax = (np.min(ydata), np.max(ydata)) if yrange is None else yrange

    x_out = np.linspace(ymin, ymax, nbins+1)
    y_out = np.empty(nbins)
    y_err = np.empty(nbins)
    dx = np.diff(x_out)[0]

    selection = in_range(xdata, xmin, xmax) & in_range(ydata, ymin, ymax)
    xdata, ydata = xdata[selection], ydata[selection]
    for i in range(nbins):
        bin_data = np.extract(in_range(ydata, x_out[i], x_out[i+1]), xdata)
        y_out[i] = np.mean(bin_data)
        y_err[i] = np.std(bin_data) / bin_data.size**0.5
    x_out += dx / 2.
    x_out = x_out[:-1]
    if drop_nan:
        selection = ~(np.isnan(y_out) | np.isnan(y_err))
        x_out = x_out[selection]
        y_out = y_out[selection]
        y_err = y_err[selection]
    return x_out, y_out, y_err


def profileXY(xdata, ydata, zdata, nbinsx, nbinsy,
              xrange=None, yrange=None, zrange=None, drop_nan=True):
    """
    Compute the xy-axis binned average of a dataset.

    Parameters
    ----------
    xdata, ydata, zdata : 1-dim np.ndarray
        x, y, z coordinates from a dataset.
    nbinsx, nbinsy : int
        Number of divisions in each axis.
    xrange : tuple of ints/floats or None, optional
        Range over the x axis. Defaults to dataset extremes.
    yrange : tuple of ints/floats or None, optional
        Range over the y axis. Defaults to dataset extremes.
    zrange : tuple of ints/floats or None, optional
        Range over the z axis. Defaults to dataset extremes.
    drop_nan : bool, optional
        Exclude empty bins. Defaults to True.

    Returns
    -------
    x_out : 1-dim np.ndarray.
        Bin centers in the x axis.
    y_out : 1-dim np.ndarray.
        Bin centers in the y axis.
    z_out : 1-dim np.ndarray
        Data average for each bin.
    z_err : 1-dim np.ndarray
        Average error for each bin.
    """
    xmin, xmax = (np.min(xdata), np.max(xdata)) if xrange is None else xrange
    ymin, ymax = (np.min(ydata), np.max(ydata)) if yrange is None else yrange
    zmin, zmax = (np.min(zdata), np.max(zdata)) if zrange is None else zrange

    x_out = np.linspace(xmin, xmax, nbinsx+1)
    y_out = np.linspace(ymin, ymax, nbinsy+1)
    z_out = np.empty((nbinsx, nbinsy))
    z_err = np.empty((nbinsx, nbinsy))
    dx = np.diff(x_out)[0]
    dy = np.diff(y_out)[0]

    selection = (in_range(xdata, xmin, xmax) &
                 in_range(ydata, ymin, ymax) &
                 in_range(zdata, zmin, zmax))
    xdata, ydata, zdata = xdata[selection], ydata[selection], zdata[selection]
    for i in range(nbinsx):
        for j in range(nbinsy):
            selection = (in_range(xdata, x_out[i], x_out[i+1]) &
                         in_range(ydata, y_out[j], y_out[j+1]))
            bin_data = np.extract(selection, zdata)
            z_out[i,j] = np.nanmean(bin_data) if bin_data.size else 0.
            z_err[i,j] = np.nanstd(bin_data) / bin_data.size**0.5 if bin_data.size else 0.
    x_out += dx / 2.
    y_out += dy / 2.
    x_out = x_out[:-1]
    y_out = y_out[:-1]
    if drop_nan:
        selection = (np.isnan(z_out) | np.isnan(z_err))
        z_out[selection] = 0
        z_err[selection] = 0
    return x_out, y_out, z_out, z_err

def projectionX(xdata, ydata, nbins, xrange=None, yrange=None):
    """
    Compute the projection of a dataset over the X axis.

    Parameters
    ----------
    xdata, ydata : 1-dim np.ndarray
        x and y coordinates from a dataset.
    nbins : int
        Number of divisions in the x axis.
    xrange : tuple of ints/floats or None, optional
        Range over the x axis. Defaults to dataset extremes.
    yrange : tuple of ints/floats or None, optional
        Range over the y axis. Defaults to dataset extremes.

    Returns
    -------
    x_out : 1-dim np.ndarray.
        Bin centers.
    y_out : 1-dim np.ndarray
        Number of points within each bin.
    """
    xmin, xmax = (np.min(xdata), np.max(xdata)) if xrange is None else xrange
    ymin, ymax = (np.min(ydata), np.max(ydata)) if yrange is None else yrange

    x_out = np.linspace(xmin, xmax, nbins+1)
    y_out = np.empty(nbins)
    dx = np.diff(x_out)[0]

    selection = in_range(xdata, xmin, xmax) & in_range(ydata, ymin, ymax)
    xdata, ydata = xdata[selection], ydata[selection]
    for i in range(nbins):
        bin_data = np.extract(in_range(xdata, x_out[i], x_out[i+1]), ydata)
        y_out[i] = bin_data.size
    x_out += dx / 2.
    x_out = x_out[:-1]
    return x_out, y_out


def projectionY(xdata, ydata, nbins, yrange=None, xrange=None):
    """
    Compute the projection of a dataset over the X axis.

    Parameters
    ----------
    xdata, ydata : 1-dim np.ndarray
        x and y coordinates from a dataset.
    nbins : int
        Number of divisions in the x axis.
    yrange : tuple of ints/floats or None, optional
        Range over the y axis. Defaults to dataset extremes.
    xrange : tuple of ints/floats or None, optional
        Range over the x axis. Defaults to dataset extremes.

    Returns
    -------
    x_out : 1-dim np.ndarray.
        Bin centers.
    y_out : 1-dim np.ndarray
        Number of points within each bin.
    """
    xmin, xmax = (np.min(xdata), np.max(xdata)) if xrange is None else xrange
    ymin, ymax = (np.min(ydata), np.max(ydata)) if yrange is None else yrange

    x_out = np.linspace(ymin, ymax, nbins+1)
    y_out = np.empty(nbins)
    dx = np.diff(x_out)[0]

    selection = in_range(xdata, xmin, xmax) & in_range(ydata, ymin, ymax)
    xdata, ydata = xdata[selection], ydata[selection]
    for i in range(nbins):
        bin_data = np.extract(in_range(ydata, x_out[i], x_out[i+1]), xdata)
        y_out[i] = bin_data.size
    x_out += dx / 2.
    x_out = x_out[:-1]
    return x_out, y_out
