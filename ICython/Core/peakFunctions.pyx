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
    return csum
