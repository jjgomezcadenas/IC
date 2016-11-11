
import numpy as np
cimport numpy as np
from scipy import signal as SGN

cpdef test():
    print('hello')

cpdef BLR(np.ndarray[np.int16_t, ndim=1] signal_daq, float coef,
           int nm, float thr1, float thr2, float thr3):
    """
    Deconvolution offline of the DAQ signal
    """
    cdef int len_signal_daq
    len_signal_daq = signal_daq.shape[0]

    cdef np.ndarray[np.float64_t, ndim=1] signal_i = signal_daq.astype(float)
    cdef np.ndarray[np.float64_t, ndim=1] B_MAU = (1./nm)*np.ones(nm, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] MAU = np.zeros(len_signal_daq, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] acum = np.zeros(len_signal_daq, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] signal_r = np.zeros(len_signal_daq, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] pulse_ = np.zeros(len_signal_daq, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] wait_ = np.zeros(len_signal_daq, dtype=np.float64)

    cdef float BASELINE,upper,lower

    MAU[0:nm] = SGN.lfilter(B_MAU,1, signal_daq[0:nm])

    acum[nm] =  MAU[nm]
    BASELINE = MAU[nm-1]

#----------

# While MAU inits BLR is switched off, thus signal_r = signal_daq

    signal_r[0:nm] = signal_daq[0:nm]

    # MAU has computed the offset using nm samples
    # now loop until the end of DAQ window

    cdef int k,j
    cdef float trigger_line, pulse_on, wait_over, offset
    cdef float part_sum

    pulse_on=0
    wait_over=0
    offset = 0

    for k in range(nm,len_signal_daq):
        trigger_line = MAU[k-1] + thr1
        pulse_[k] = pulse_on
        wait_[k] = wait_over

        # condition: raw signal raises above trigger line and
        # we are not in the tail
        # (wait_over == 0)

        if signal_daq[k] > trigger_line and wait_over == 0:

            # if the pulse just started pulse_on = 0.
            # In this case compute the offset as value
            #of the MAU before pulse starts (at k-1)

            if pulse_on == 0: # pulse just started
                #offset computed as the value of MAU before pulse starts
                offset = MAU[k-1]
                pulse_on = 1

            #Pulse is on: Freeze the MAU
            MAU[k] = MAU[k-1]
            signal_i[k] = MAU[k-1]  #signal_i follows the MAU

            #update recovered signal, correcting by offset

            acum[k] = acum[k-1] + signal_daq[k] - offset;
            signal_r[k] = signal_daq[k] + coef*acum[k]

            #signal_r[k] = signal_daq[k] + signal_daq[k]*(coef/2) + coef*acum[k-1]
            #acum[k] = acum[k-1] + signal_daq[k] - offset

        else:  #no signal or raw signal has dropped below threshold
        # but raw signal can be negative for a while and still contribute to the
        # reconstructed signal.

            if pulse_on == 1: #reconstructed signal still on
                # switch the pulse off only when recovered signal
                # drops below threshold
                # lide the MAU, still frozen.
                # keep recovering signal

                MAU[k] = MAU[k-1]
                signal_i[k] = MAU[k-1]

                acum[k] = acum[k-1] + signal_daq[k] - offset;
                signal_r[k] = signal_daq[k] + coef*acum[k]

                #signal_r[k] = signal_daq[k] + signal_daq[k]*(coef/2) + coef*acum[k-1]
                #acum[k] = acum[k-1] + signal_daq[k] - offset

                #if the recovered signal drops before trigger line
                #rec pulse is over!
                if signal_r[k] < trigger_line + thr2:
                    wait_over = 1  #start tail compensation
                    pulse_on = 0   #recovered pulse is over

            else:  #recovered signal has droped below trigger line
            #need to compensate the tail to avoid drifting due to erros in
            #baseline calculatoin

                if wait_over == 1: #compensating pulse
                    # recovered signal and raw signal
                    #must be equal within a threshold
                    # otherwise keep compensating pluse

                    if signal_daq[k-1] < signal_r[k-1] - thr3:
                        # raw signal still below recovered signal
                        # keep compensating pulse
                        # is the recovered signal near offset?
                        upper = offset + (thr3 + thr2)
                        lower = offset - (thr3 + thr2)

                        if lower < signal_r[k-1] < upper:
                            # we are near offset, activate MAU.

                            signal_i[k] = signal_r[k-1]
                            part_sum = 0.
                            for j in range(k-nm, k):
                                part_sum += signal_i[j]/nm
                            MAU[k] = part_sum

                            #MAU[k] = np.sum(signal_i[k-nm:k])*1./nm

                        else:
                            # rec signal not near offset MAU frozen
                            MAU[k] = MAU[k-1]
                            signal_i[k] = MAU[k-1]

                        # keep adding recovered signal
                        acum[k] = acum[k-1] + signal_daq[k] - offset;
                        signal_r[k] = signal_daq[k] + coef*acum[k]
                        #signal_r[k] = signal_daq[k] + signal_daq[k]*(coef/2) + coef*acum[k-1]
                        #acum[k] = acum[k-1] + signal_daq[k] - offset

                    else:  # raw signal above recovered signal: we are done

                        wait_over = 0
                        acum[k] = MAU[k-1]
                        signal_r[k] = signal_daq[k]
                        signal_i[k] = signal_r[k]
                        part_sum = 0.
                        for j in range(k-nm, k):
                            part_sum += signal_i[j]/nm
                        MAU[k] = part_sum
                        #MAU[k] = np.sum(signal_i[k-nm:k])*1./nm

                else: #signal still not found

                    #update MAU and signals
                    part_sum = 0.
                    for j in range(k-nm, k):
                        part_sum += signal_i[j]/nm
                    MAU[k] = part_sum
                    #MAU[k] = np.sum(signal_i[k-nm:k]*1.)/nm
                    acum[k] = MAU[k-1]
                    signal_r[k] = signal_daq[k]
                    signal_i[k] = signal_r[k]
    #energy = np.dot(pulse_f,(signal_r-BASELINE))

    signal_r = signal_r - BASELINE

    #return  signal_r.astype(int)
    return  signal_r.astype(int), MAU, pulse_, wait_

cpdef deconvolve_signal_acum(np.ndarray[np.int16_t, ndim=1] signal_i,
                             int n_baseline=500,
                             float coef_clean=2.905447E-06,
                             float coef_blr=1.632411E-03,
                             float noise_rms = 0.9,
                             float thr_trigger=5, float thr_acum=800,
                             float coeff_acum = 0.9995):

    """
    The accumulator approach by Master VHB

    """

    cdef float coef = coef_blr
    cdef int nm = n_baseline

    cdef int len_signal_daq = len(signal_i)
    cdef np.ndarray[np.float64_t, ndim=1] signal_r = np.zeros(len_signal_daq,
                                                              dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] acum = np.zeros(len_signal_daq,
                                                          dtype=np.double)

    cdef np.ndarray[np.float64_t, ndim=1] signal_daq = signal_i.astype(float)
    cdef int j
    cdef float baseline = 0.
    for j in range(0,nm):
        baseline += signal_daq[j]
    baseline /= nm
    cdef float trigger_line
    trigger_line = thr_trigger*noise_rms

    signal_daq =  baseline - signal_daq

    b_cf, a_cf = SGN.butter(1, coef_clean, 'high', analog=False);
    signal_daq = SGN.lfilter(b_cf,a_cf,signal_daq)

    # BLR
    signal_r[0:nm] = signal_daq[0:nm]
    cdef int k
    for k in range(nm,len_signal_daq):
        # condition: raw signal raises above trigger line
        if (signal_daq[k] > trigger_line) or (acum[k-1] > thr_acum):

            signal_r[k] = signal_daq[k] + signal_daq[k]*(coef/2.0) + coef*acum[k-1]
            acum[k] = acum[k-1] + signal_daq[k]

        else:
            signal_r[k] = signal_daq[k]
            # deplete the accumulator before or after the signal to avoid runoffs
            if (acum[k-1]>0):
                acum[k]=acum[k-1]*coeff_acum

    return signal_r.astype(int), acum
