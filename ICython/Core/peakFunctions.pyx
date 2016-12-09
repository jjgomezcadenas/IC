"""
Cython version of some peak functions
JJGC December, 2016
"""
cimport numpy as np
import numpy as np
from scipy import signal

cpdef calibrated_pmt_sum(double [:, :] CWF,
                         double [:] adc_to_pes,
                         int n_MAU=200,
                         double thr_MAU=5):
    """
    Computes the ZS calibrated sum of the PMTs
    after correcting the baseline with a MAU to suppress low frequency noise.
    input:
    CWF:    Corrected waveform (passed by BLR)
    adc_to_pes: a vector with calibration constants
    n_MAU:  length of the MAU window
    thr_MAU: treshold above MAU to select sample

    """

    cdef int j, k
    cdef int NPMT = CWF.shape[0]
    cdef int NWF = CWF.shape[1]
    cdef double [:] MAU = np.array(np.ones(n_MAU),
                                   dtype=np.double)*(1./float(n_MAU))


    # CWF if above MAU threshold
    cdef double [:, :] pmt_thr = np.zeros((NPMT,NWF), dtype=np.double)
    cdef double [:] csum = np.zeros(NWF, dtype=np.double)
    cdef double [:] MAU_pmt = np.zeros(NWF, dtype=np.double)

    for j in range(NPMT):
        # MAU for each of the PMTs, following the waveform
        MAU_pmt = signal.lfilter(MAU,1,CWF[j,:])

        for k in range(NWF):
            if CWF[j,k] > MAU_pmt[k] + thr_MAU:
                pmt_thr[j,k] = CWF[j,k]

    for j in range(NPMT):
        for k in range(NWF):
            csum[k] += pmt_thr[j, k]*1./adc_to_pes[j]
    return np.asarray(csum)


cpdef wfzs(double [:] wf, double threshold=0):
    """
    takes a waveform wf and returns the values of the wf above threshold:
    if the input waveform is of the form [e1,e2,...en],
    where ei is the energy of sample i,
    then then the algorithm returns a vector [e1,e2...ek],
    where k <=n and ei > threshold and
    a vector of indexes [i1,i2...ik] which label the position
    of the zswf of [e1,e2...ek]
    For example if the input waveform is:
    [1,2,3,5,7,8,9,9,10,9,8,5,7,5,6,4,1] and the trhesold is 5
    then the algoritm returns
    a vector of amplitudes [7,8,9,9,10,9,8,7,6] and a vector of indexes
    [4,5,6,7,8,9,10,12,14]

    """
    cdef int len_wf = wf.shape[0]
    cdef double [:] wfzs_e = np.zeros(len_wf, dtype=np.double)
    cdef int [:] wfzs_i = np.zeros(len_wf, dtype=np.int32)

    cdef int i,j
    j=0
    for i in range(len_wf):
        if wf[i] > threshold:
            wfzs_e[j] = wf[i]
            wfzs_i[j] = i
            j+=1

    cdef double [:] wfzs_ene = np.zeros(j, dtype=np.double)
    cdef int [:] wfzs_indx = np.zeros(j, dtype=np.int32)

    for i in range(j):
        wfzs_ene[i] =  wfzs_e[i]
        wfzs_indx[i] = wfzs_i[i]

    return np.asarray(wfzs_ene), np.asarray(wfzs_indx)


cpdef time_from_index(int [:] indx):
    """
    returns the times (in ns) corresponding to the indexes in indx
    """
    cdef int len_indx = indx.shape[0]
    cdef double [:] tzs = np.zeros(len_indx, dtype=np.double)

    cdef int i
    cdef double step = 25 #ns
    for i in range(len_indx):
        tzs[i] = step * float(indx[i])

    return np.asarray(tzs)


cpdef find_S12(double [:] wfzs, int [:] index,
               double tmin=0, double tmax=1e+6,
               int stride=4):
    """
    Find S1/S2 peaks.
    input:
    wfzs:   a vector containining the zero supressed wf
    indx:   a vector of indexes
    returns a dictionary

    do not interrupt the peak if next sample comes within stride
    accept the peak only if within [lmin, lmax)
    accept the peak only if within [tmin, tmax)
    """

    cdef double [:] P = wfzs
    cdef double [:] T = time_from_index(index)

    assert(len(wfzs) == len(index))

    cdef dict S12 = {}
    cdef int pulse_on = 1
    cdef int i, j

    cdef list s12 = []
    S12[0] = s12
    S12[0].append([T[0],P[0]])

    j=0
    for i in range(1,len(wfzs)) :

        if T[i] > tmax:
            break

        if T[i] < tmin:
            continue

        if index[i] - stride > index[i-1]:  #new s12
            j+=1
            s12 = []
            S12[j] = s12
            S12[j].append([T[i],P[i]])
        else:
            S12[j].append([T[i],P[i]])

    return S12
