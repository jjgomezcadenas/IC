"""
Parameters defining NEW PMT + FEE
"""
import math
import numpy as np

import system_of_units as units


np.set_printoptions(threshold='nan')

# offset = 700
offset = 2500
ceiling = 4096
n_PMT = 12
PMT_GAIN = 1.7e6
MAU_WindowSize = 512
# V_GAIN = 570*units.ohm*0.5 #Gain current --> V = IR Gain = R = 250 Ohm
V_GAIN = 250*units.ohm
C = 6.2*units.nF  # decoupling capacitor in pF
R = 2350*units.ohm  # resistor in Ohms
C12 = np.array([6.02975448, 6.22547194, 6.0671337, 6.22159457,
                6.29999787, 6.09892384, 6.18289435, 6.21775591,
                6.19306671, 6.30518792, 6.20359891,  6.31231192])*units.nF

AC = np.array([0.001652348, 0.00163149, 0.001631774, 0.00162948,
               0.001582021, 0.001594688, 0.001624928, 0.001642566,
               0.001626805, 0.001620873, 0.001600633, 0.001611269])

time_step = 1*units.ns
time_DAQ = 25*units.ns
time_bin = time_DAQ/time_step
f_sample = 1.0/time_step  # this is the inverse of the time step 1 ns
f_sample_DAQ = 1.0/time_DAQ

freq_HPF = 1.0/(2*math.pi*R*C)  # HPF filter
freq_LPF = 3E6*units.hertz

NBITS = 12
NBITS_FRAC = 11
NBITS_acum = 32
NBITS_FRAC_acum = 20

LSB = 2*units.volt/2**NBITS  # value of LSB
LSB_OUT = LSB/1.25  # effective value to convert mV to ADC counts
voltsToAdc = LSB_OUT

NOISE_FEE_rms = 0.3*units.mV
NOISE_FEE_adc = NOISE_FEE_rms/voltsToAdc
MAU_thr = LSB/units.volt  # MAU operates in ADC counts

NOISE_FEE = NOISE_FEE_rms
NOISE_ADC = NOISE_FEE_adc


def pmt_gain():
    return PMT_GAIN


def sampling_time():
    return time_step


def sampling_DAQ():
    return time_DAQ


def spe_i_to_v():
    return V_GAIN


def spe_i_to_adc():
    return V_GAIN/LSB_OUT


def decoupling_capacitor():
    return C


def decoupling_resitor():
    return R


def f_HPF():
    return freq_HPF


def f_LPF():
    return freq_LPF


def W_HPF_fine():
    return freq_HPF/f_sample


def W_HPF_daq():
    return freq_HPF/f_sample_DAQ


# def W_HPF_mc():
#     return freq_HPF/f_sample_MC


def W_LPF_fine():
    return freq_LPF/f_sample


def W_LPF_daq():
    return freq_LPF/f_sample_DAQ


# def W_LPF_mc():
#     return freq_LPF/f_sample_MC


def print_FEE():
    print """
    NEW FEE: DEFAULT PARAMETERS
    PMT gain = %7.2g
    sampling time: (fine) = %7.2f ns (DAQ) = %7.2f ns
    decoupling capacitor = %7.2f nF
    decoupling resistor = %7.2f ohm
    HPF frequency = %7.2f Hz  W_HPF_fine = %7.2g W_HPF_daq = %7.2g
    LPF frequency = %7.2f Hz  W_LPF_fine = %7.2g W_LPF_daq = %7.2g
    noise = %7.2f mV
    noise (adc) = %7.2f
    vots to adc factor = %7.2f
    """.format(pmt_gain(), sampling_time()/units.ns, sampling_DAQ()/units.ns,
               decoupling_capacitor()/units.nF, decoupling_resitor()/units.ohm,
               f_HPF()/units.hertz, W_HPF_fine(), W_HPF_daq(),
               f_LPF()/units.hertz, W_LPF_fine(), W_LPF_daq(),
               NOISE_FEE_rms/units.mV, NOISE_FEE_rms/voltsToAdc,
               units.mV/voltsToAdc)
    print "decoupling capacitors for energy plane = %s".format(C12/units.nF)


def pulse_area(pulse):
    return np.sum(pulse)


def pulse_area_positive(pulse):
    return np.sum(pulse[np.where(pulse > 0)])


def pulse_area_threshold(pulse, thr):
    return np.sum(pulse[np.where(pulse > thr)])


def pulse_len_positive(pulse):
    return len(pulse[np.where(pulse > 0)])
