"""
Cython version of some peak functions
JJGC December, 2016

calibrated_pmt_sum computes the ZS calibrated sum of the PMTs
after correcting the baseline with a MAU to suppress low frequency noise.
input:
CWF:    Corrected waveform (passed by BLR)
adc_to_pes: a vector with calibration constants
n_MAU:  length of the MAU window
thr_MAU: treshold above MAU to select sample

"""
cimport numpy as np
import numpy as np
from scipy import signal

cpdef calibrated_pmt_sum(double [:, :] CWF,
                         double [:] adc_to_pes,
                         int n_MAU=*,
                         double thr_MAU=*)
