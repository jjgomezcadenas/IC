#
# control plots of TAMARA
#
# author: JA Hernando
# date: 21/11/2016

import numpy as np

import qxyFunctions as qxyf
import Calib.calib as cb   # TODO move functions from here to its place
from Core.Bridges import Signal

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
matplotlib.style.use('ggplot')


LDIM = 250.

# TAM = True
DOR = False  # Ana PMaps from Dorothea
if (DOR):
    print('nana: PMaps from DOROTHEA')

UNKNOWN, S1, S2 = Signal.UNKNOWN, Signal.S1, Signal.S2

TSCALE = 0.025
SIAXIS, TIMEAXIS = 0, 1
if (DOR):
    TSCALE = 1.  # Dorothea PMAPS
    SIAXIS, TIMEAXIS = 1, 0  # Dorothea PMaps

# -------------------------
# General
# --------------------------


def plot_pan(pan, labels=[], bins=100):
    if (len(labels) <= 0):
        labels = pan.columns
    nx, ny, figsize = cb.plt_subplots(len(labels))
    fig, axs = plt.subplots(nx, ny, figsize=figsize)
    axs = axs.ravel()
    for i, label in enumerate(labels):
        axs[i].hist(pan[label].values, bins=bins)
        axs[i].set_xlabel(label)
    fig.tight_layout()
    return fig


# --------------------------
# Plot WF functions
# ---------------------------


def polo_pmt_wfs(rfws, cwfs=[], blrs=[], title='pmt'):
    """ plots a list of wfs
    """
    n = len(rfws)
    nx, ny, figsize = cb.plt_subplots(n)
    fig, axs = plt.subplots(nx, ny, figsize=figsize)
    axs = axs.ravel()
    for i, ax in enumerate(axs):
        if (len(rfws) > i):
            ax.plot(rfws[i], color='red', alpha=0.5)
        if (len(cwfs) > i):
            ax.plot(cwfs[i], color='blue', alpha=0.5)
        if (len(blrs) > i):
            ax.plot(blrs[i], color='black', alpha=0.5)
        ax.set_xlabel('time')
        ax.set_ylabel('adc')
        ax.set_title(title+' ['+str(i)+']')
    fig.tight_layout()
    return fig


# ----------------------
# Plot of pmaps
# ---------------------


def plot_pmap_wfs(pmap, title=''):
    """ plot the wfs of the peaks of a pmap
    """
    fig, axs = plt.subplots(2, 1, figsize=(8, 6))
    axs[0].set_title('cathode '+title)
    axs[0].set_xlim(0., 800)
    axs[1].set_xlim(0., 800.)
    axs[1].set_title('anode '+title)
    xcolors = {UNKNOWN: 'g', S1: 'b*', S2: 'B'}
    for peak in pmap.peaks:
        xcol = xcolors[peak.signal]
        axs[0].plot(TSCALE*peak.times, peak.cathode, xcol)
        ctimes = TSCALE*peak.times[0] + np.arange(peak.anode.shape[TIMEAXIS])
        axs[1].plot(ctimes, np.nansum(peak.anode, SIAXIS), xcol)
    fig.tight_layout()
    return fig


def plot_peak(peak, xs, ys, title='', q0=0, bary=True,
              xlim=(-LDIM, LDIM), ylim=(-LDIM, LDIM)):
    """ plot the peak wfs (anode and cathode) and SiPM plane
    """
    fig = plt.figure(figsize=(4*3, 3*3))
    gs = gridspec.GridSpec(4, 5)
    ax1 = plt.subplot(gs[0:2, 0:2])
    ctimes = peak.times[0]*TSCALE + np.arange(len(peak.cathode))*TSCALE
    cmax = np.max(peak.cathode)
    amax = np.max(np.nansum(peak.anode, axis=SIAXIS))
    asamples = peak.anode.shape[TIMEAXIS]
    atimes = peak.times[0]*TSCALE + np.arange(asamples)
    ax1.plot(ctimes, peak.cathode/cmax)
    if (not DOR):
        atimes += 1
    ax1.plot(atimes, np.nansum(peak.anode, axis=SIAXIS)/amax)
    ax1.set_xlabel(u'time $\mu$s')
    ax1.set_ylabel('charge (normalized)')
    ax2 = plt.subplot(gs[2:, 0:2])
    nsipms = peak.anode.shape[SIAXIS]
    ax2.scatter(range(nsipms), np.nansum(peak.anode, axis=TIMEAXIS))
    ax2.set_xlabel('sipm index')
    ax2.set_ylabel('charge ps')
    ax3 = plt.subplot(gs[:, 2:])
    qs = np.nansum(peak.anode, axis=TIMEAXIS)
    cqs, cxs, cys = qxyf.qxy_slice(qs, xs, ys, q0=q0, xlim=xlim, ylim=ylim)
    ct = ax3.scatter(cxs, cys, c=cqs, cmap=plt.cm.jet, s=50)
    if (bary):
        point, cov = qxyf.qxy_point(qs, xs, ys, q0=q0, xlim=xlim, ylim=ylim)
        cq, cx, cy = point
        ex, ey = np.sqrt(cov[0]), np.sqrt(cov[1])
        ax3.errorbar(cx, cy, xerr=ex, yerr=ey, c='black', marker='*', ms=20)
    if (title != ''):
        ax3.set_title(title)
    ax3.set_xlabel('x (mm)')
    ax3.set_ylabel('y (mm)')
    fig.colorbar(ct, ax=ax3)
    fig.tight_layout()
    return fig
