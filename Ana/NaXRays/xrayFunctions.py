#
# functions to study the X-ryas of Na22 in IC
#
# author: JA Hernando
# date: 21/11/2016

# from __future__ import print_function

import numpy as np
from scipy.optimize import least_squares

import Core.fitFunctions as fitf

import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')


KALPHA = 29.7  # keV
KBETA = 33.8  # keV
QBB = 2458.  # keV
FWHM = 2.*np.sqrt(2.*np.log(2.))

TSCALE = 0.025  # ns time slice time


# ----------------------------------
# Generic Utilities - TO Export
# ----------------------------------


def np_profile(xs, ys, xedges):
    """ compute the profile of ys vs xs
    Parameters:
    -----------
    xs: iterable
        the x values, the axis of the bins
    ys: iterable
        the ys values, the axis to compute the profile
    xedges: iterable - ordered list
        of the x-values of the bins
    Returns:
    --------
    cxs: iterable
        the center of the profile bins
    cys: iteramble
        the average of ys in the profile bins
    ceys: the errors on y
    """
    cxs, cys, ceys = [], [], []
    for x0, x1 in zip(xedges[:-1], xedges[1:]):
        sel = (xs >= x0) & (xs < x1)
        selxs = xs[sel]
        selys = ys[sel]
        cxs.append(np.mean(selxs))
        cys.append(np.mean(selys))
        ceys.append(np.std(selys)/np.sqrt(len(selys)))
    npa = np.array
    return (npa(cxs), npa(cys), npa(ceys))


def fit(ps0, xs, ys, eys, fun, bounds=None):
    """ Wrapper to least_squares.
    Inputs:
        ps0: the initial guess parameters,
        xs, ys: np arrays with the x, y data,
        fun: is a function of the ps.
    Outputs:
        result with the result of the fit. It has as members:
        chi2: the chi2 of the fit
        fun: the function with the optimized parameters,
        x: the optimized parameters
        cov: the cov. matrix of the optimized paramters
    """
    def errfunc(ps):
        return (ys-fun(ps, xs))/eys

    def func(ps):
        def ffun(xs):
            return fun(ps, xs)
        return ffun

    if (bounds):
        result = least_squares(errfunc, ps0, bounds=bounds)
    else:
        result = least_squares(errfunc, ps0)

    chi2 = -1.
    cov = []
    if (result.success):
        chi2 = np.sum(result.fun**2)/(1.*(len(ys)-len(ps0)))
        cov = (np.linalg.inv(np.dot(result.jac.T, result.jac)))*chi2
    result.chi2 = chi2
    result.cov = cov
    result.fun = func(result.x)

    return result


def fsline(ps, xs):
    return ps[0]+ps[1]*xs


def fparabola(ps, xs):
    return ps[0]+ps[1]*xs+ps[2]*xs*xs


def fxrays(ps, xs, iray=2):
    n0, f, mu, sigma = ps[0], ps[1], ps[2], ps[3]
    mu2 = (KBETA/KALPHA)*mu
    sigma2 = sigma*np.sqrt(KALPHA/KBETA)
    factor1 = (1./(np.sqrt(2*np.pi)*sigma))
    factor2 = (1./(np.sqrt(2*np.pi)*sigma2))
    g1 = factor1*n0*(1-f)*np.exp(-(xs-mu)*(xs-mu)/(2*sigma*sigma))
    g2 = factor2*n0*f*np.exp(-(xs-mu2)*(xs-mu2)/(2*sigma2*sigma2))
    if (iray == 0):
        return g1
    if (iray == 1):
        return g2
    return g1+g2


# ----------------------------
# Energy resolution
# ----------------------------


def energy_resolution(enes, range, bins=50):
    """ fit the list of energy values to 2 Xrays.
    Input:
        enes: list of energy values
        range: tupe with the lower, upper bound of energy
        bins: number of bins of the histogram
    Output:
        fitresult (with chi2, x, fun, cov attributes)
        xres (resolution FWHM at Kalpha)
        bbres (resolution FWHH at Qbb)
    """

    ys, xedges = np.histogram(enes, bins=bins, range=range)
    xs = 0.5*(xedges[:-1]+xedges[1:])
    eys = 1.+np.sqrt(ys)
    xbounds = ((0, 0., 4000., 100.), (400000., 1., 6000., 800.))
    ps0 = (200000., 0.2, 5500., 600.)
    fitresult = fit(ps0, xs, ys, eys, fxrays, bounds=xbounds)

    chi2, pshat, cov = fitresult.chi2, fitresult.x, fitresult.cov
    mu, sig = pshat[2], pshat[3]
    smu, ssig = np.sqrt(cov[2][2]), np.sqrt(cov[3][3])
    rat = sig/mu
    erat = np.sqrt((ssig*ssig)/(mu*mu) + (rat*rat*smu*smu)/(mu*mu))
    rat = float(100.*rat*FWHM)
    erat = float(100.*erat*FWHM)
    factor = np.sqrt(KALPHA/QBB)
    xres = (rat, erat)
    bbres = (factor*rat, factor*erat)

    return fitresult, xres, bbres

# -----------------------
# attachment
# -----------------------


def attachment(times, enes, bins, range):
    """ Computes the attachment
    Input:
        times: list of times
        enes: list of energies
        bins: int of the profile
        range: tupe with range on times
    Output:
        fitresult: to a straight line
        tuple: with the x, y , error-y of the bins of the profile
        mm: value and error of the slope (pes/us)
        tau: value and error of the life time (ms)
    """
    edges = np.linspace(range[0], range[1], bins)
    cxs, cys, ceys = np_profile(times, enes, edges)
    ps0 = (5000.,  -1.)
    fitresult = fit(ps0, cxs, cys, ceys, fsline)

    pshat, cov = fitresult.x, fitresult.cov
    n0, m = pshat[0], pshat[1]
    en0, em = np.sqrt(cov[0][0]), np.sqrt(cov[1][1])
    tau = -1.*n0/m
    stau = np.sqrt(en0*en0/(m*m)+tau*tau*em*em/(m*m))

    mm = (m, em)
    mtau = (tau*1.e-3, stau*1.e-3)

    print('mu = ', mm)
    print('tau = ', mtau)
    return fitresult, (cxs, cys, ceys), mm, mtau


def correct_attachment(times, enes, mu):
    enes1 = enes - times*mu
    return enes1


# ----------------------------
# Geometrical corrections
# ----------------------------

def radius_correction(rads, enes, bins, range):
    """ compute the readius correction
    Input:
        rads: list of the values of the radius
        enes: list of the values of the energy
        bins: int bins of the profile
        range: tuple of the range in radius (lower, upper) bounds
    Output:
        fitresult: result of the fit to a 2d polynomial
        tuple: with the x, y, error-y value sof the profile
    """

    def remove_corner(zi):
        return True
        # r, q = zi[0], zi[1]
        # if (q > 7000.):
        #     return False
        # return (q < 7000.-(2000./(220-150.))*(r-150.))

    zz = zip(rads, enes)
    zf = [zi for zi in zz if remove_corner(zi)]

    frads = np.array([zi[0] for zi in zf])
    fenes = np.array([zi[1] for zi in zf])

    edges = np.linspace(range[0], range[1], bins)
    cxs, cys, ceys = np_profile(frads, fenes, edges)

    ps0 = (5000., -1., -1.)
    fitresult = fit(ps0, cxs, cys, ceys, fparabola)
    pshat = fitresult.x
    print(' radius - fit result ', pshat)

    return fitresult, (cxs, cys, ceys)


def correct_energy_by_radius(enes, rads, mu, c):
    enes1 = enes-m*rads-c*rads*rads
    return enes1
