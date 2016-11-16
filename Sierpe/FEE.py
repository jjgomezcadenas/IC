"""
Classes and functions describing the electronics of the
PMT plane FEE.
(full model)
VH, JJGC, November, 2016
"""

from __future__ import print_function

import numpy as np
from scipy import signal
import system_of_units as units

# globals describing FEE
PMT_GAIN = 1.7e6
FEE_GAIN = 582.237*units.ohm
DAQ_GAIN = 1.25
NBITS = 12
LSB = 2.0*units.V/2**NBITS/DAQ_GAIN
NOISE_I = LSB/(FEE_GAIN*DAQ_GAIN)
NOISE_DAQ = 0.313*units.mV

C2 = 8*units.nF
C1 = 2714*units.nF
R1 = 1567*units.ohm
Zin = 62*units.ohm
t_sample = 25*units.ns
f_sample = 1./t_sample
f_mc = 1./(1*units.ns)
f_LPF1 = 3*units.MHZ
f_LPF2 = 10*units.MHZ
ADC_TO_PES = 20  # nominal factor, comes out from spe area
OFFSET = 2500  # offset adc


def i_to_adc():
    """
    current to adc counts
    """
    return FEE_GAIN/LSB


def i_to_v():
    """
    current to voltage
    """
    return FEE_GAIN


def v_to_adc():
    """
    voltage to adc
    """
    return 1./LSB


class SPE:
    """
    Represents a single photo-electron in the PMT
    """

    def __init__(self, pmt_gain=PMT_GAIN, x_slope=5.*units.ns,
                 x_flat=1.*units.ns):

        self.pmt_gain = pmt_gain
        self.x_slope = x_slope
        self.x_flat = x_flat
        self.spe_base = self.x_slope + self.x_flat
        self.spe_length = 2*self.x_slope + self.x_flat

        # current
        self.A = self.pmt_gain*units.eplus/self.spe_base

        time_step = 1.*units.ns
        self.t = np.arange(0, self.spe_length, time_step)
        nns = int(self.x_slope/time_step)
        nnf = int(self.x_flat/time_step)
        rise = np.linspace(0, self.A, num=nns)
        fall = np.linspace(self.A, 0, num=nns)
        flat = self.A*np.ones(nnf)
        self.spe = np.concatenate((rise, flat, fall))

    def __str__(self):
        """
        output the class to string
        """
        s = """
        (PMT gain = {0:5.2g}, amplitude = {1:5.2g} muA
         slope = {2:5.2f} ns, flat = {3:5.2f} ns)
        """.format(self.pmt_gain, self.A/units.muA,
                   self.x_slope/units.ns,
                   self.x_flat/units.ns)
        return s

    def __repr__(self):
        """
        class representation
        """
        return self.__str__()


def spe_pulse(spe, t0=100*units.ns, tmax=1e+6*units.ns,
              time_step=1*units.ns):
    """
    input: an instance of class SPE
    Returns a SPE pulse at time t0
    with baseline extending in steps of time_step from 0 to tmax
    determined by DELTA_L
    """
    n = int(t0/time_step)
    nmax = int(tmax/time_step)

    DELTA = np.zeros(nmax)   # Dirac delta of size DELTA_L
    DELTA[n] = 1
    # step = time_step/units.ns
    # spe_pulse_t =np.arange(0, len(DELTA) + len(self.spe) -1, step)
    spe_pulse = signal.convolve(DELTA, spe.spe)

    return spe_pulse


def spe_pulse_train(spe,
                    signal_start=2000*units.ns,
                    signal_length=5000*units.ns, daq_window=20*units.mus,
                    time_step=1*units.ns):
    """
    Input: an instance of class SPE
    Returns a train of SPE pulses between signal_start
    and start+length in daq_window separated by tstep
    """
    nmin = int(signal_start/time_step)
    nmax = int((signal_start + signal_length)/time_step)
    NMAX = int(daq_window/time_step)
    # step = time_step/units.ns

    DELTA = np.zeros(NMAX)
    DELTA[nmin:nmax+1] = 1
    # spe_pulse_t =np.arange(0,len(DELTA) + len(self.spe) -1,step)
    spe_pulse = signal.convolve(DELTA, spe.spe)

    return spe_pulse


def spe_pulse_from_vector(spe, cnt):
    """
    input: an instance of spe
    Returns a train of SPE pulses corresponding to vector cnt
    """

    spe_pulse = signal.convolve(cnt[0:-len(spe.spe)+1], spe.spe)
    return spe_pulse


class FEE:
    """
    Complete model of Front-end electronics.
    """

    def __init__(self, gain=FEE_GAIN,
                 c2=C2, c1=C1, r1=R1, zin=Zin, fsample=f_sample,
                 flpf1=f_LPF1, flpf2=f_LPF2,
                 noise_FEEPMB_rms=NOISE_I, noise_DAQ_rms=NOISE_DAQ, lsb=LSB):

        self.R1 = r1
        self.Zin = zin
        self.C2 = c2
        self.C1 = c1
        self.GAIN = gain
        self.A1 = self.R1 * self.Zin/(self.R1 + self.Zin)  # ohms
        self.A2 = gain/self.A1  # ohms/ohms = []
        self.R = self.R1 + self.Zin
        self.Cr = 1. + self.C1/self.C2
        self.C = self.C1/self.Cr
        self.ZC = self.Zin/self.Cr

        self.f_sample = fsample
        self.freq_LHPF = 1./(self.R * self.C)
        self.freq_LPF1 = flpf1*2*np.pi
        self.freq_LPF2 = flpf2*2*np.pi

        self.freq_LHPFd = self.freq_LHPF/(self.f_sample*np.pi)
        self.freq_LPF1d = self.freq_LPF1/(self.f_sample*np.pi)
        self.freq_LPF2d = self.freq_LPF2/(self.f_sample*np.pi)

        self.noise_FEEPMB_rms = noise_FEEPMB_rms
        self.LSB = lsb
        self.voltsToAdc = self.LSB/units.volt
        self.DAQnoise_rms = noise_DAQ_rms

    def __str__(self):
        """
        output the class to string
        """
        s = """
        (C1 = {0:7.1f} nf,
         C2 = {1:7.1f} nf,
         R1 = {2:7.1f} ohm,
         Zin = {3:7.1f} ohm,
         gain = {4:7.1f} ohm,
         f_sample = {5:7.1f} MHZ,
         freq_LHPF = {6:7.2f} kHz,
         freq_LPF1 = {7:7.2f} MHZ,
         freq_LPF2 = {8:7.2f} MHZ,
         freq_LHPFd = {9:8.5f},
         freq_LPF1d = {10:7.2f},
         freq_LPF2d = {11:7.2f},
         noise_FEEPMB_rms = {12:7.2f} muA,
         LSB = {13:7.2g} mV,
         volts to adc = {14:7.2g},
         DAQnoise_rms = {15:7.2g}
        )
        """.format(self.C1/units.nF,
                   self.C2/units.nF,
                   self.R1/units.ohm,
                   self.Zin/units.ohm,
                   self.GAIN/units.ohm,
                   self.f_sample/units.MHZ,
                   self.freq_LHPF/(units.kHz*2*np.pi),
                   self.freq_LPF1/(units.MHZ*2*np.pi),
                   self.freq_LPF2/(units.MHZ*2*np.pi),
                   self.freq_LHPFd,
                   self.freq_LPF1d,
                   self.freq_LPF2d,
                   self.noise_FEEPMB_rms/units.muA,
                   self.LSB/units.mV,
                   self.voltsToAdc,
                   self.DAQnoise_rms/units.mV)
        return s

    def __repr__(self):
        """
        class representation
        """
        return self.__str__()


def noise_adc(fee, signal_in_adc):
    """
    Equivalent Noise of the DAQ added at the output
    of the system
    input: a signal (in units of adc counts)
           an instance of FEE class
    output: a signal with DAQ noise added
    """
    noise_daq = fee.DAQnoise_rms*v_to_adc()
    return signal_in_adc + np.random.normal(0,
                                            noise_daq,
                                            len(signal_in_adc))


def filter_fee(feep):
    """
    input: an instance of class FEE
    output: buttersworth parameters of the equivalent FEE filter
    """

    # high pass butterswoth filter ~1/RC
    b1, a1 = signal.butter(1,
                           feep.freq_LHPFd,
                           'high', analog=False)
    b2, a2 = signal.butter(1,
                           feep.freq_LHPFd,
                           'low', analog=False)

    b0 = b2*feep.ZC + b1*feep.A1  # in ohms
    a0 = a1

    # LPF order 1
    b1l, a1l = signal.butter(1,
                             feep.freq_LPF1d, 'low',
                             analog=False)
    # LPF order 4
    b2l, a2l = signal.butter(4,
                             feep.freq_LPF2d, 'low',
                             analog=False)
    # convolve HPF, LPF1
    a_aux = np.convolve(a0, a1l, mode='full')
    b_aux = np.convolve(b0, b1l, mode='full')
    # convolve HPF+LPF1, LPF2
    a = np.convolve(a_aux, a2l, mode='full')
    b_aux2 = np.convolve(b_aux, b2l, mode='full')
    b = feep.A2*b_aux2  # in ohms

    return b, a


def filter_cleaner(feep):
    """
    cleans the input signal
    """
    freq_zero = 1./(feep.R1*feep.C1)
    freq_zerod = freq_zero/(feep.f_sample*np.pi)
    b, a = signal.butter(1, freq_zerod, 'high', analog=False)

    return b, a


def signal_v_fee(feep, signal_i):
    """
    input: signal_i = signal current (i = A)
           instance of class FEE
    output: signal_v (in volts) with effect FEE

    ++++++++++++++++++++++++++++++++++++++++++++++++
    +++++++++++ PMT+FEE NOISE ADDED HERE +++++++++++
    ++++++++++++++++++++++++++++++++++++++++++++++++

    """
    if (feep.noise_FEEPMB_rms == 0.0):
        noise_FEEin = np.zeros(len(signal_i))
    else:
        noise_FEEin = np.random.normal(0,
                                       feep.noise_FEEPMB_rms,
                                       len(signal_i))

    # Equivalent Noise of the FEE + PMT BASE added at the input
    # of the system to get the noise filtering effect

    b, a = filter_fee(feep)  # b in ohms
    # filtered signal in I*R = V
    return signal.lfilter(b, a, signal_i + noise_FEEin)


def signal_clean(feep, signal_fee):
    """
    input: signal_fee = adc, convoluted
           instance of class FEE
    output: signal_c cleaning filter passed

    ++++++++++++++++++++++++++++++++++++++++++++++++
    +++++++++++ PMT+FEE NOISE ADDED HERE +++++++++++
    ++++++++++++++++++++++++++++++++++++++++++++++++

    """
    b, a = filter_cleaner(feep)
    return signal.lfilter(b, a, signal_fee)


def daq_decimator(f_sample1, f_sample2, signal_in):
    """
    downscales the signal vector according to the
    scale defined by f_sample1 (1 GHZ) and
    f_sample2 (40 Mhz).
    Includes anti-aliasing filter
    """

    scale = int(f_sample1/f_sample2)
    return signal.decimate(signal_in, scale, ftype='fir')
