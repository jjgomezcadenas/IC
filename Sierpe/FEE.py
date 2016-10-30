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


class SPE:
    """
    Represents a single photo-electron in the PMT
    """

    def __init__(self, pmt_gain=4.5e6, x_slope=5.*units.ns,
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
        (PMT gain = {0:5.2g}, slope = {1:5.2f} ns, flat = {2:5.2f} ns)
        """.format(self.pmt_gain, self.x_slope/units.ns,
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

    def __init__(self, gain=582.2*units.ohm,
                 C2=8*units.nF, C1=2714*units.nF,
                 R1=1567*units.ohm, Zin=62*units.ohm,
                 f_sample=1./(25*units.ns), f_LPF1=3E6*units.hertz,
                 f_LPF2=10E6*units.hertz,
                 noise_FEEPMB_rms=0.20*units.muA,
                 nbits=12, DAQnoise_rms=0.313*units.mV):

        self.R1 = R1
        self.Zin = Zin
        self.C2 = C2
        self.C1 = C1
        self.GAIN = gain
        self.A1 = R1*Zin/(R1+Zin)
        self.A2 = gain/self.A1
        self.R = self.R1+self.Zin
        self.Cr = 1+self.C1/self.C2
        self.C = self.C1/self.Cr
        self.ZC = self.Zin/self.Cr

        self.f_sample = f_sample
        self.freq_LHPF = 1./(self.R*self.C)
        self.freq_LPF1 = f_LPF1*2*np.pi
        self.freq_LPF2 = f_LPF2*2*np.pi

        self.freq_LHPFd = self.freq_LHPF/(self.f_sample*np.pi)
        self.freq_LPF1d = self.freq_LPF1/(self.f_sample*np.pi)
        self.freq_LPF2d = self.freq_LPF2/(self.f_sample*np.pi)

        self.noise_FEEPMB_rms = noise_FEEPMB_rms

        self.NBITS = nbits
        self.LSB = 2*units.volt/2**self.NBITS/1.25
        self.voltsToAdc = self.LSB/units.volt
        self.DAQnoise_rms = DAQnoise_rms

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
         NIBTS = {13:d},
         LSB = {14:7.2g} mV,
         volts to adc = {15:7.2g},
         DAQnoise_rms = {16:7.2g}
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
                   self.NBITS, self.LSB/units.mV,
                   self.voltsToAdc,
                   self.DAQnoise_rms/units.mV)
        return s

    def __repr__(self):
        """
        class representation
        """
        return self.__str__()


def i_to_adc(fee):
    """
    input: instances of classes SimpleFEE, SimpleDAQ
    outputs: current to adc counts
    """
    return fee.GAIN/fee.voltsToAdc


def i_to_v(fee):
    """
    input: instance of class FEE
    output: current to voltage
    """
    return fee.GAIN


def noise_adc(fee, signal_in):
    """
    Equivalent Noise of the DAQ added at the output
    of the system
    input: a signal (in units of voltage)
           an instance of FEE class
    output: a signal with DAQ noise added
    """

    return signal_in + np.random.normal(0,
                                        fee.DAQnoise_rms,
                                        len(signal_in))


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

    b0 = b2*feep.ZC + b1*feep.A1  # in ohm
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
    b = feep.A2*b_aux2  # b has no units

    return b, a


def filter_cleaner(feep):
    """
    cleans the input signal
    """
    freq_zero = 1./(feep.R1*feep.C1)
    freq_zerod = freq_zero/(feep.f_sample*np.pi)
    b, a = signal.butter(1, freq_zerod, 'high', analog=False)

    return b, a


def signal_fee(feep, signal_in):
    """
    input: instance of class feep and a signal
    outputs: signal convolved with effect FEE

    ++++++++++++++++++++++++++++++++++++++++++++++++
    +++++++++++ PMT+FEE NOISE ADDED HERE +++++++++++
    ++++++++++++++++++++++++++++++++++++++++++++++++

    """
    if (feep.noise_FEEPMB_rms == 0.0):
        noise_FEEin = np.zeros(len(signal_in))
    else:
        noise_FEEin = np.random.normal(0,
                                       feep.noise_FEEPMB_rms,
                                       len(signal_in))

    # Equivalent Noise of the FEE + PMT BASE added at the input
    # of the system to get the noise filtering effect

    b, a = filter_fee(feep)
    return signal.lfilter(b, a, signal_in+noise_FEEin)
