# functions to study the X-ryas of Na22 in IC
#
# author: JA Hernando
# date: 21/11/2016

import numpy as np
import xrayFunctions as xrayf

import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')


LDIM = 240.


KALPHA = 29.7  # keV
KBETA = 33.8  # keV
QBB = 2458.  # keV
FWHM = 2.*np.sqrt(2.*np.log(2))

TSCALE = 0.025  # time in ms of the slices

# -----------------------
# Plots
# -----------------------


def plot_xybarycenter(pan, bins=60, ldim=LDIM):
    fig, axs = plt.subplots(2, 2, figsize=(4*3, 3*3))
    axs[0][0].scatter(pan['x'].values, pan['y'].values)
    axs[0][0].set_xlabel('x (mm)')
    axs[0][0].set_ylabel('y (mm)')
    axs[0][0].set_xlim((-ldim, ldim))
    axs[0][0].set_ylim((-ldim, ldim))
    axs[0][1].hist(pan['y'].values, bins=bins, range=(-ldim, ldim))
    axs[0][1].set_xlabel(' y (mm)')
    axs[1][0].hist(pan['x'].values, bins=bins, range=(-ldim, ldim))
    axs[1][0].set_xlabel('x (mm)')
    xs, ys = pan['x'].values, pan['y'].values
    cc, cx, cy, cimg = axs[1][1].hist2d(xs, ys, bins=bins/2, cmap=plt.cm.jet)
    axs[1][1].grid(False)
    axs[1][1].set_xlabel('x (mm)')
    axs[1][1].set_ylabel('y (mm)')
    fig.colorbar(cimg, ax=axs[1][1])
    return fig


def plot_energy_resolution(enes, bins, range, textpos):
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    ys, edges, cc = ax.hist(enes, bins=bins, range=range,
                            color='blue', alpha=0.4)
    fresult, xres, bbres = xrayf.energy_resolution(enes, bins=bins,
                                                   range=range)
    xs = 0.5*(edges[1:]+edges[:-1])
    pshat = fresult.x
    fun = xrayf.fxrays
    ax.plot(xs, fun(pshat, xs), color='blue', linewidth=2)
    ax.plot(xs, fun(pshat, xs, iray=0), 'g:', linewidth=2)
    ax.plot(xs, fun(pshat, xs, iray=1), 'r:', linewidth=2)
    text = r'FWHM Qbb = {0:2.1f} $\pm$ {1:2.1f} %'.format(bbres[0],
                                                          bbres[1])
    ax.text(textpos[0], textpos[1], text, fontsize=14)
    return fig, (fresult, xres, bbres)


def plot_attachment(xs, ys, bins, range, textpos):
    fig, axs = plt.subplots(2, 1, figsize=(8, 6))

    axs[0].scatter(xs, ys, alpha=0.5)
    axs[0].set_xlim(range)
    axs[0].set_xlabel(u't ($\mu$s)')
    axs[0].set_ylabel('pes')

    fresult, cc, mm, mtau = xrayf.attachment(xs, ys, bins=bins, range=range)
    fun = fresult.fun
    
    cxs, cys, ceys = cc[0], cc[1], cc[2]
    axs[0].errorbar(cxs, cys, yerr=ceys)
    axs[1].errorbar(cxs, cys, yerr=ceys)
    xxs = np.linspace(range[0], range[1], 100.)
    axs[1].plot(xxs, fun(xxs), '-')
    axs[1].set_xlim(range)
    # axs[1].set_ylim((3000., 6000.))
    axs[1].set_xlabel(u'time ($\mu$s)')
    axs[1].set_ylabel('pes')

    text = u'$\\tau$ = {0:.3f} $\pm$ {1:.3f} ms'.format(mtau[0], mtau[1])
    axs[1].text(textpos[0], textpos[1], text, fontsize=14)
    fig.tight_layout()
    return fig, (fresult, cc, mm, mtau)


def plot_radius_correction(rads, enes, bins=50, range=(0., 220.)):

    fig, axs = plt.subplots(2, 1, figsize=(8, 6))

    axs[0].scatter(rads, enes, c='blue', alpha=0.5)
    axs[0].set_xlim(-10., 240.)
    axs[0].set_xlabel('radius (mm)')
    axs[0].set_ylabel(' pes')

    fitresult, cc = xrayf.radius_correction(rads, enes, bins, range)

    cxs, cys, ceys = cc
    axs[0].errorbar(cxs, cys, yerr=ceys, c='black')
    axs[1].errorbar(cxs, cys, yerr=ceys)
    axs[1].set_xlabel('radius (mm)')
    axs[1].set_ylabel(' pes')

    psfun = fitresult.fun
    pshat = fitresult.x
    print(' fit result ', pshat)
    xxs = np.linspace(10., 200., 10)
    axs[1].plot(xxs, psfun(xxs), c='black')
    fig.tight_layout()
    return fig, (fitresult, cc)
