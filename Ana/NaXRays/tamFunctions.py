#
# functions for TAMARA: to produce pMaps from wfs and seelct peaks in PMaps
#
# author: JA Hernando
# date: 21/11/2016

import operator
import numpy as np

from Core.Bridges import Signal

import qxyFunctions as qxyf


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

# peaks selection
Q0 = 4.0  # min charge to start a peak (previous 4.)
S1WIDTH = 20.   # 0.5 us
S1QTOT = 20.0   # 5 pes in S1 (previus 10. )
S2QIMIN = 4.0  # 5 pes in 25 ns slice for S1 (previous 20.)
S2SEPARATION = 80.  # 2 us!
TWIDTH = 40.  # 25 ns to 1 us
Q0BARY = 1.2
# Minimum charge in SiPMs to enter into the barycenter calculation


# --------------------------
# WF functions
# --------------------------


def wfs_cal(wfs, cals=1.):
    """ add the array of wfs, normalize to calib
    """
    if (isinstance(cals, float)):
        cals = np.ones(len(wfs))*cals
    iwfs = [wf/abs(cal) for wf, cal in zip(wfs, cals)]
    return iwfs


def wfs_onewf(wfs, cals=1.):
    """ add the array of wfs, normalize to calib
    """
    iwfs = wfs_cal(wfs, cals)
    return reduce(operator.add, iwfs)


def sipmwf_basewfs(wfs, irange=(0, 90)):
    """ compute the baseline
    """
    nsample = len(wfs[0])
    base = np.mean(wfs[:, irange[0]:irange[1]], axis=1)
    # base = na.wfs_aveonebin(wfs, irange=(irange[0], irange[1]))
    basewf = np.array([bi*np.ones(nsample) for bi in base])
    return basewf


def sipmwf_cwfs(wfs, cals=16.):
    """ remove the SiPM baseline and apply calibration
    """
    basewf = sipmwf_basewfs(wfs)
    cwfs = (wfs-basewf)/cals
    return cwfs


# --------------------------
# PMaps functions from wfs
# --------------------------

#  Classify peack from the WF

def peaks_find(wf, q0=Q0, separation=S2SEPARATION):
    """ returns a list of intervals where the charge is over a threshold (q0min)
    it adds interval if the distance between them is smaller than (separation)
    """
    ss = []
    index = wf > q0
    it0 = False
    i0 = -10*separation
    for i, it in enumerate(index):
        if (it and (not it0)):
            if ((i-i0) > separation):
                si = [i]
            else:
                si = ss.pop()
                si += range(si[-1]+1, i+1)
        elif (it and it0):
            si.append(i)
        elif ((not it) and it0):
            ss.append(si)
            i0 = i-1
        it0 = it
    return ss


def peaks_classify(wf, ss, width=S1WIDTH, s1qtotmin=S1QTOT, s2qmin=S2QIMIN):
    """ classifies the peaks in 0,1,2 for S1, S2 or unknown (0).
    wf is the wavefunction
    ss is the list of the interval of the Peacks
    width is the maximum distance to consider a S1
    s1qtotmin is the minimun total charge in the S1 peak
    s2qmin is the minimum charge of a bin int he S2 peak
    returns a list with (0,1,2) types in the same order of the ss intervals
    """
    ssi = [UNKNOWN, ]*len(ss)
    for i, si in enumerate(ss):
        qtot = np.sum(wf[si])
        qmax = np.max(wf[si])
        if (len(si) < width):
            if (qtot > s1qtotmin):
                ssi[i] = S1
        elif (qmax > s2qmin):
            ssi[i] = S2
    return ssi


def peaks_anode(ss, siwfs, twidth=TWIDTH):
    """ it returns the anode information of the peaks given by the ss list
    """
    qs = []
    for si in ss:
        i0, i1 = si[0], si[-1]
        j0, j1 = max(0, int(i0/twidth)), min(int(i1/twidth), 800)
        iq = siwfs[:, j0:j1]
        qs.append(iq)
    return qs


def pmap_summary(pmap):
    s = 'pmap: n '+str(len(pmap.peaks))
    s += ', types '+str([peak.signal for peak in pmap.peaks])
    s += ', times '+str([peak.tmin+0.5*peak.width for peak in pmap.peaks])
    s += ', cathode '+str([peak.cathode_integral for peak in pmap.peaks])
    s += ', anode '+str([peak.anode_integral for peak in pmap.peaks])
    return s


# Selection of Peaks and PMaps
# -----------------------------


def peak_has_energy(peak, range, scale=1.):
    ene = scale*peak.cathode_integral
    return (ene >= range[0]) and (ene <= range[1])


def peak_has_time(peak, range, scale=1.):
    t0, t1 = scale*peak.tmin, scale*peak.tmax
    return ((t0 >= range[0]) and (t1 <= range[1]))


def peak_has_type(peak, ptype):
    return (peak.signal == ptype)


def s1_has_time(peak, range, scale=1):
    if (peak.signal != S1):
        return False
    return peak_has_time(peak, range=range, scale=scale)


def s2_has_energy(peak, range, scale=1.):
    if (peak.signal != S2):
        return False
    return peak_has_energy(peak, range=range, scale=scale)


def pmap_is_gold(pmap):
    ns1 = len([peak for peak in pmap.peaks if (peak.signal == S1)])
    ns2 = len([peak for peak in pmap.peaks if (peak.signal == S2)])
    return ((ns1 == 1) and (ns2 >= 1))


def pmap_is_platinum(pmap):
    ns1 = len([peak for peak in pmap.peaks if (peak.signal == S1)])
    ns2 = len([peak for peak in pmap.peaks if (peak.signal == S2)])
    return ((ns1 == 1) and (ns2 == 1))


def pmap_peaks_with_condition(pmap, condition):
    return [peak for peak in pmap.peaks if condition(peak)]


def pmaps_peaks_with_condition(pmaps, condition):
    ns, peaks = [], []
    for pmap in pmaps:
        xpeaks = pmap_peaks_with_condition(pmap, condition)
        ns.append(len(xpeaks))
        if (len(xpeaks) >= 0):
            peaks += xpeaks
    return ns, peaks


#  Functions for Peaks
# ----------------------------

def peak_barycenter(peak, xs, ys, q0bary=Q0BARY):
    qs = np.nansum(peak.anode, axis=TIMEAXIS)
    wd = peak.anode.shape[TIMEAXIS]
    return qxyf.qxy_point(qs, xs, ys, q0=q0bary*wd)
