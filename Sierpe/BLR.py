"""
BLR algorithm
VH, RE, JJGC
May-October, 2016

Changelog:
13/10: added reco_signal to describe the core of BLR. Include latest recipe from VH

"""
from __future__ import print_function
from Util import *
from LogConfig import *
from scipy import signal as SGN

import FEParam as FP
import FEE2 as FE

def reco_signal(k, signal_daq,coef,acum,offset):
    """
    The core of the BLR algorithm. Obtains the reconstructed signal by
    adding the raw signal multiplied by the accumulator times the deconvolution coefficients
    """
    signal_reco = signal_daq[k] + signal_daq[k]*(coef/2) + coef*acum[k-1]
    accum = acum[k-1] + signal_daq[k] - offset
    return signal_reco, accum

class SBLR:
    """
    Signals BLR: a simple class to hold of the signals relevant for BLR
    """
    def __init__(self, len_signal_daq):
        """
        len_signal_daq: length of the signal to deconvolute


        """

        self.MAU = np.zeros(len_signal_daq, dtype=np.double)
        self.acum = np.zeros(len_signal_daq, dtype=np.double)
        self.signal_r = np.zeros(len_signal_daq, dtype=np.double)
        self.pulse_on = np.zeros(len_signal_daq, dtype=np.double)
        self.wait_over = np.zeros(len_signal_daq, dtype=np.double)
        self.offset = np.zeros(len_signal_daq, dtype=np.double)
        self.BASELINE = 0


def BLR(signal_daq, coef, mau_len=250, thr1 = 3*FP.NOISE_ADC, thr2 = 0,
        thr3 = FP.NOISE_ADC, log='INFO'):

    """
    Deconvolution offline of the DAQ signal using a MAU
    moving window-average filter of a vector data. See notebook
    y(n) = (1/WindowSize)(x(n) + x(n-1) + ... + x(n-windowSize))
    in a filter operation filter(b,a,x):
    b = (1/WindowSize)*ones(WindowSize) = (1/WS)*[1,1,1,...]: numerator
    a = 1 : denominator
    y = filter(b,a,x)
    y[0] = b[0]*x[0] = (1/WS) * x[0]
    y[1] = (1/WS) * (x[0] + x[1])
    y[WS-1] = mean(x[0:WS])
    y[WS] = mean(x[1:WS+1])
    and so on
    """
    lg = 'logging.'+log
    logger.setLevel(eval(lg))
    len_signal_daq = len(signal_daq)
    sblr = SBLR(len_signal_daq)

    signal_i = np.copy(signal_daq) #uses to update MAU while procesing signal
    nm = mau_len
    B_MAU = (1./nm)*np.ones(nm)

#   MAU averages the signal in the initial tranch
#    allows to compute the baseline of the signal

    sblr.MAU[0:nm] = SGN.lfilter(B_MAU,1, signal_daq[0:nm])
    sblr.acum[nm] =  sblr.MAU[nm]
    sblr.BASELINE = sblr.MAU[nm-1]

    logging.debug("""-->BLR:
                     MAU_LEN={}
                     thr1 = {}, thr2 = {}, thr3 = {} =""".format(
                     mau_len, thr1, thr2, thr3))
    logging.debug("n = {}, acum[n] = {} BASELINE ={}".format(nm, sblr.acum[nm],sblr.BASELINE))

#----------

# While MAU inits BLR is switched off, thus signal_r = signal_daq

    sblr.signal_r[0:nm] = signal_daq[0:nm]
    pulse_on=0
    wait_over=0
    offset = 0

    # MAU has computed the offset using nm samples
    # now loop until the end of DAQ window

    logging.debug("nm = {}".format(nm))

    for k in range(nm,len_signal_daq):
        trigger_line = sblr.MAU[k-1] + thr1
        sblr.pulse_on[k] = pulse_on
        sblr.wait_over[k] = wait_over
        sblr.offset[k] = offset

        # condition: raw signal raises above trigger line and
        # we are not in the tail
        # (wait_over == 0)
        if signal_daq[k] > trigger_line and wait_over == 0:

            # if the pulse just started pulse_on = 0.
            # In this case compute the offset as value
            #of the MAU before pulse starts (at k-1)

            if pulse_on == 0: # pulse just started
                #offset computed as the value of MAU before pulse starts
                offset = sblr.MAU[k-1]
                pulse_on = 1

            #Pulse is on: Freeze the MAU
            sblr.MAU[k] = sblr.MAU[k-1]
            signal_i[k] = sblr.MAU[k-1]  #signal_i follows the MAU

            #update recovered signal, correcting by offset
            sblr.signal_r[k], sblr.acum[k] = reco_signal(k, signal_daq,coef,sblr.acum,offset)

        else:  #no signal or raw signal has dropped below threshold

        # but raw signal can be negative for a while and still contribute to the
        # reconstructed signal.

            if pulse_on == 1: #reconstructed signal still on
                # switch the pulse off only when recovered signal
                #drops below threshold
                #slide the MAU, still frozen.
                # keep recovering signal

                sblr.MAU[k] = sblr.MAU[k-1]
                signal_i[k] = sblr.MAU[k-1]
                sblr.signal_r[k], sblr.acum[k] = reco_signal(k, signal_daq,coef,sblr.acum,offset)
                # sblr.acum[k] = sblr.acum[k-1] + signal_daq[k] - offset;
                # sblr.signal_r[k] = signal_daq[k] + coef*sblr.acum[k]

                #if the recovered signal drops before trigger line
                #rec pulse is over!
                if sblr.signal_r[k] < trigger_line + thr2:
                    wait_over = 1  #start tail compensation
                    pulse_on = 0   #recovered pulse is over


            else:  #recovered signal has droped below trigger line
            #need to compensate the tail to avoid drifting due to erros in
            #baseline calculatoin

                if wait_over == 1: #compensating pulse
                    # recovered signal and raw signal
                    #must be equal within a threshold
                    # otherwise keep compensating pluse

                    if signal_daq[k-1] < sblr.signal_r[k-1] - thr3:
                        # raw signal still below recovered signal
                        # keep compensating pulse
                        # is the recovered signal near offset?
                        upper = offset + (thr3 + thr2)
                        lower = offset - (thr3 + thr2)

                        if lower < sblr.signal_r[k-1] < upper:
                            # we are near offset, activate MAU.

                            signal_i[k] = sblr.signal_r[k-1]
                            sblr.MAU[k] = np.sum(signal_i[k-nm:k])*1./nm

                        else:
                            # rec signal not near offset MAU frozen
                            sblr.MAU[k] = sblr.MAU[k-1]
                            signal_i[k] = sblr.MAU[k-1]

                        # keep adding recovered signal
                        sblr.signal_r[k], sblr.acum[k] = reco_signal(k, signal_daq,coef,sblr.acum,offset)
                        # sblr.acum[k] = sblr.acum[k-1] + signal_daq[k] - sblr.MAU[k]
                        # sblr.signal_r[k] = signal_daq[k] + coef*sblr.acum[k]

                    else:  # raw signal above recovered signal: we are done

                        wait_over = 0
                        sblr.acum[k] = sblr.MAU[k-1]
                        sblr.signal_r[k] = signal_daq[k]
                        signal_i[k] = sblr.signal_r[k]
                        sblr.MAU[k] = np.sum(signal_i[k-nm:k])*1./nm


                else: #signal still not found

                    #update MAU and signals
                    sblr.MAU[k] = np.sum(signal_i[k-nm:k]*1.)/nm
                    sblr.acum[k] = sblr.MAU[k-1]
                    sblr.signal_r[k] = signal_daq[k]
                    signal_i[k] = sblr.signal_r[k]
    #energy = np.dot(pulse_f,(signal_r-BASELINE))

    sblr.signal_r = sblr.signal_r - sblr.BASELINE
    return  sblr

def accumulator_coefficients(CA,NPMT,len_WF):
    """
    Compute the accumulator coefficients for DBLR
    It computes the inverse function of the HPF and takes
    the accumulator as the value of the function anywhere
    but the first bin (the inverse is a step function with
    constant value equal to the accumulator)
    CA are the values of the capacitances defining the filter
    (1/(2CR)) for each PMT
    """
    coef_acc =np.zeros(NPMT, dtype=np.double)

    signal_t = np.arange(0.0, len_WF*1., 1., dtype=np.double)

    for j in range(NPMT):
        fee = FE.FEE(C=CA[j],R= FP.R, f=FP.freq_LPF, RG=FP.V_GAIN)
        signal_inv_daq = fee.InverseSignalDAQ(signal_t)  #inverse function
        coef_acc[j] = signal_inv_daq[10] #any index is valid, function is flat

    return coef_acc


def discharge_acum(length_d=5000, tau=2500, compress=0.005):
    t_discharge = np.arange(0,length_d,1,dtype=np.double)
    exp =  np.exp(-(t_discharge-length_d/2)/tau)
    cf = 1./(1. + exp)
    discharge_curve = compress*(1. - cf) + (1. - compress)
    return discharge_curve


def deconvolve_signal_acum(signal_i, n_baseline=500, noise_rms= 0.8,
                      coef_clean=2.905447E-06, coef_blr=1.632411E-03,
                      thr_trigger=5, thr_acum=1000, coeff_acum = 0.9995,
                      acum_discharge_length = 5000, acum_tau=2500, acum_compress=0.0025,
                      filter_c=True):

    """
    The accumulator approach by Master VHB
    decorations by JJGC
    """

    coef = coef_blr
    nm = n_baseline
    len_signal_daq = len(signal_i)

    signal_r = np.zeros(len_signal_daq, dtype=np.double)
    acum = np.zeros(len_signal_daq, dtype=np.double)
    pulse_on = np.zeros(len_signal_daq, dtype=np.int8)
    j_reg = np.zeros(len_signal_daq, dtype=np.double)
    trigger = np.zeros(len_signal_daq, dtype=np.double)

    # signal_daq in floats
    signal_daq = signal_i.astype(float)

    # compute baseline and noise
    baseline = np.mean(signal_daq[0:nm])
    noise_rms = np.std(signal_daq[0:nm],ddof=1)

    # change sign and subtract baseline
    signal_daq =  baseline - signal_daq

    # clean the signal_daq
    if filter_c==True:
        b_cf, a_cf = SGN.butter(1, coef_clean, 'high', analog=False);
        signal_daq = SGN.lfilter(b_cf,a_cf,signal_daq)

    # compute discharge curve
    discharge_curve = discharge_acum(length_d=acum_discharge_length,
                                     tau=acum_tau,
                                     compress=acum_compress)


    # signal_r equals to signal_d (baseline suppressed and change signed)
    # for the first nm samples
    signal_r[0:nm] = signal_daq[0:nm]
    p_on = 0
    j=0
    # print ("baseline = {}, noise (LSB_rms) = {} MAU[nm] ={} ".format(
    #        baseline, noise_rms, MAU[nm]))

    trigger_line = thr_trigger*noise_rms  # fixed trigger line

    for k in range(nm,len_signal_daq):
        pulse_on[k] = p_on
        trigger[k] = trigger_line

        # update recovered signal
        signal_r[k] = signal_daq[k] + signal_daq[k]*(coef/2.0) + coef*acum[k-1]

        # condition: raw signal raises above trigger line
        # once signal raises above trigger line condition is on until
        # accumulator drops below thr_acum
        if (signal_daq[k] > trigger_line) or (acum[k-1] > thr_acum):
            if p_on == 0:
                p_on = 1

            # update accumulator
            acum[k] = acum[k-1] + signal_daq[k]
        else:
            if p_on == 1:
                p_on = 0
                j = 0

            # discharge acumulator
            if acum[k-1]>1:
                acum[k] = acum[k-1]*discharge_curve[j]
                if j<acum_discharge_length-1:
                    j=j+1
                else:
                    j=acum_discharge_length-1
            else:
                acum[k]=0
                j=0
        j_reg[k]=j

    BLR={}
    BLR['acum'] = acum
    BLR['pulse_on'] = pulse_on
    BLR['signal_daq'] = signal_daq
    BLR['signal_r'] = signal_r
    BLR['j_reg'] = j_reg
    BLR['trigger'] = trigger

    return pd.DataFrame(BLR)


def DBLR(pmtrd, event_number, coeff_acc, mau_len=250,
         thr1 = FP.NOISE_ADC, thr2=0, thr3 = FP.NOISE_ADC, log='INFO'):
    """
    Peform Base line Restoration
    coeff_acc is an array with the coefficients of the accumulator
    Threshold 1 is used to decide when raw signal raises up above trigger line
    Threshold 2 is used to decide when reconstructed signal is above trigger line
    Threshold 3 is used to compare Raw and Rec signal
    """
    perform_blr = lambda wf,coef: BLR(FP.ceiling - wf, coef, mau_len, thr1, thr2, thr3, log)
    return map( perform_blr, pmtrd[event_number], coeff_acc )
