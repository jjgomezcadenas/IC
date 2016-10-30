"""
Classes and functions describing the electronics of the
PMT plane FEE (simple model)
VH, JJGC, November, 2016
"""

from __future__ import print_function

import numpy as np
from scipy import signal
import system_of_units as units


class SimpleFee:
    """
    Simplified model of Front-end electronics.
    """

    def __init__(self, gain=582.2*units.ohm, C2=8*units.nF,
                 R1=1567*units.ohm, Zin=62*units.ohm,
                 f_sample=1./(25*units.ns), f_LPF1=3E6*units.hertz,
                 f_LPF2=10E6*units.hertz, noise_rms=0.3*units.mV):

        self.R1 = R1
        self.Zin = Zin
        self.C2 = C2
        self.GAIN = gain
        self.R = self.R1+self.Zin
        self.f_sample = f_sample
        self.freq_HPF = 1./(self.R*self.C2)
        self.freq_LPF1 = f_LPF1*2*np.pi
        self.freq_LPF2 = f_LPF2*2*np.pi
        self.freq_HPFd = self.freq_HPF/(self.f_sample*np.pi)
        self.freq_LPF1d = self.freq_LPF1/(self.f_sample*np.pi)
        self.freq_LPF2d = self.freq_LPF2/(self.f_sample*np.pi)

        self.noise_rms = noise_rms

    def __str__(self):
        """
        output the class to string
        """
        s = """
        (C2 = {0:7.1f} nf, R1 = {1:7.1f} ohm, Zin = {2:7.1f} ohm,
         gain = {3:7.1f} ohm, f_sample = {4:7.1f} MH,
         self.freq_HPF = {5:7.2f} kHz, self.freq_LPF1 = {6:7.2f} MHZ,
         self.freq_LPF2 = {7:7.2f} MHZ,
         self.freq_HPFd = {8:8.5f}  self.freq_LPF1d = {9:7.2f}
         self.freq_LPF2d = {10:7.2f},
         noise_rms = {11:7.2f} mV)
        """.format(self.C2/units.nF, self.R1/units.ohm,
                   self.Zin/units.ohm, self.GAIN/units.ohm,
                   self.f_sample/units.MHZ,
                   self.freq_HPF/(units.kHz*2*np.pi),
                   self.freq_LPF1/(units.MHZ*2*np.pi),
                   self.freq_LPF2/(units.MHZ*2*np.pi),
                   self.freq_HPFd, self.freq_LPF1d, self.freq_LPF2d,
                   self.noise_rms/units.mV)
        return s

    def __repr__(self):
        """
        class representation
        """
        return self.__str__()


class SimpleDAQ:
    """
    Simplified model of DAQ
    """

    def __init__(self, nbits=12):

        self.NBITS = nbits
        self.LSB = 2*units.volt/2**self.NBITS
        self.voltsToAdc = self.LSB/(units.volt*1.25)

    def __str__(self):
        """
        output the class to string
        """
        s = """
        (NIBTS = {0:d}, LSB = {1:7.2g} volts/adc, volts to adc = {2:7.2g})
        """.format(self.NBITS, self.LSB/units.volt, self.voltsToAdc)
        return s

    def __repr__(self):
        """
        class representation
        """
        return self.__str__()

def filter_sfee(sfe):
    """
    input: an instance of class SimpleFee
    output: buttersworth parameters of the equivalent FEE filter
    """
    # high pass butterswoth filter ~1/RC
    b0, a0 = signal.butter(1, sfe.freq_HPFd, 'high', analog=False)
    # LPF order 1
    b1, a1 = signal.butter(1, sfe.freq_LPF1d, 'low', analog=False)
    # LPF order 4
    b2, a2 = signal.butter(4, sfe.freq_LPF2d, 'low', analog=False)
    # convolve HPF, LPF1
    a_aux = np.convolve(a0, a1, mode='full')
    b_aux = np.convolve(b0, b1, mode='full')
    # convolve HPF+LPF1, LPF2
    a = np.convolve(a_aux, a2, mode='full')
    b_aux2 = np.convolve(b_aux, b2, mode='full')
    b = sfe.GAIN*b_aux2
    return b, a


def filter_sfee_lpf(sfe):
    """
    input: an instance of class SimpleFee
    output: buttersworth parameters of the equivalent LPT FEE filter
    """
    # LPF order 1
    b1, a1 = signal.butter(1, sfe.freq_LPF1d, 'low', analog=False)
    # LPF order 4
    b2, a2 = signal.butter(4, sfe.freq_LPF2d, 'low', analog=False)
    # convolve LPF1, LPF2
    a = np.convolve(a1, a2, mode='full')
    b_aux = np.convolve(b1, b2, mode='full')
    b = sfe.GAIN*b_aux
    return b, a


def filter_sfee_hpf(sfe):
    """
    input: an instance of class SimpleFee
    output: buttersworth parameters of the HPF FEE filter
    """
    # high pass butterswoth filter ~1/RC
    b0, a = signal.butter(1, sfe.freq_HPFd, 'high', analog=False)
    b = sfe.GAIN*b0
    return b, a

def signal_fee(sfe, signal_in):
    """
    input: instance of class sfe and a signal
    outputs: signal convolved with effect FEE
    """
    b, a = filter_fee(sfe)
    return signal.lfilter(b/sfe.GAIN, a, signal_in)


def signal_fee_lpf(sfe, signal_in):
    """
    input: instance of class sfe and a signal
    outputs: signal convolved with LPF
    """
    b, a = filter_fee_lpf(sfe)
    return signal.lfilter(b/sfe.GAIN, a, signal_in)


def signal_fee_hpf(sfe, signal_in):
    """
    input: instance of class sfe and a signal
    outputs: signal convolved with HPF
    """
    b, a = filter_fee_hpf(sfe)
    return signal.lfilter(b/sfe.GAIN, a, signal_in)


class FeePmt:
    """
    Complete model of Front-end electronics.
    """

    def __init__(self, gain=582.2*units.ohm,
                 C2=8*units.nF, C1=2714*units.nF,
                 R1=1567*units.ohm, Zin=62*units.ohm,
                 f_sample=1./(25*units.ns), f_LPF1=3E6*units.hertz,
                 f_LPF2=10E6*units.hertz,
                 noise_FEEPMB_rms=0.20*units.muA):

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

        self.LSB = (2.0/2**12/1.25)*units.V
        # LSB measured at FEE output

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
         self.freq_LHPF = {6:7.2f} kHz,
         self.freq_LPF1 = {7:7.2f} MHZ,
         self.freq_LPF2 = {8:7.2f} MHZ,
         self.freq_LHPFd = {9:8.5f},
         self.freq_LPF1d = {10:7.2f},
         self.freq_LPF2d = {11:7.2f},
         self.noise_FEEPMB_rms = {12:7.2f} muA,
         self.LSB = {13:7.2f} mV)
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
                   self.LSB/units.mV)
        return s

    def __repr__(self):
        """
        class representation
        """
        return self.__str__()


class FeeDAQ:
    """
    Simplified model of DAQ + NOISE
    """

    def __init__(self, nbits=12, DAQnoise_rms=0.313*units.mV):

        self.NBITS = nbits
        self.LSB = 2*units.volt/2**self.NBITS
        self.voltsToAdc = self.LSB/(units.volt*1.25)
        self.DAQnoise_rms = DAQnoise_rms

    def __str__(self):
        """
        output the class to string
        """
        s = """
        (NIBTS = {0:d}, LSB = {1:7.2g} volts/adc,
        volts to adc = {2:7.2g},
         noise = {3:7.2g})
        """.format(self.NBITS, self.LSB/units.volt,
                   self.voltsToAdc, self.DAQnoise_rms/units.mV)
        return s

    def __repr__(self):
        """
        class representation
        """
        return self.__str__()
