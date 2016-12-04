"""
Functions to find peaks, S12 selection etc.
JJGC and GML December 2016

"""
from __future__ import print_function

import math
import numpy as np
import pandas as pd


def pmt_sum(CWF, adc_to_pes):
    """
    input: A CWF list or array
           a vector with the adc_to_pes values (must be positive)
    returns: the sum of CWF, in pes
    """

    NPMT = len(CWF)
    NWF = len(CWF[0])

    csum = np.zeros(NWF, dtype=np.double)
    for j in range(NPMT):
        csum += CWF[j]*1./adc_to_pes[j]
    return csum


def wfdf(time,energy_pes):
    """
    takes two vectors (time, energy) and returns a data frame representing a waveform
    """
    swf = {}
    swf['time_mus'] = time/units.mus
    swf['ene_pes'] = energy_pes
    return pd.DataFrame(swf)


def wf_thr(wf,threshold=0):
    """
    return a zero supressed waveform (more generally, the vaules of wf above threshold)
    """
    return wf.loc[lambda df: df.ene_pes.values >threshold, :]


def find_peaks(wfzs, stride=4, lmin=8):
    """
    Find peaks.
    do not interrupt the peak if next sample comes within stride
    accept the peak only if larger than lmin samples
    """
    T = wfzs['time_mus'].values
    P = wfzs['ene_pes'].values
    I = wfzs.index.values

    S12 = {}
    pulse_on = 1
    j=0

    S12[0] = []
    S12[0].append([T[0],P[0],I[0]])

    for i in range(1,len(wfzs)) :
        if wfzs.index[i]-stride > wfzs.index[i-1]:  #new s12
            j+=1
            S12[j] = []
            S12[j].append([T[i],P[i],I[i]])
        else:
            S12[j].append([T[i],P[i],I[i]])

    S12L=[]
    for i in S12.keys():
        if len(S12[i]) > lmin:
            S12L.append(pd.DataFrame(S12[i], columns=['time_mus','ene_pes','index']))
    return S12L
