import Sierpe.FEE as FE
import Core.system_of_units as units
from nose.tools import *
import numpy as np

def test_fee_params():
    """
    Check the values of the FEE params
    """
    assert_equal(FE.PMT_GAIN, 1.7e6)
    assert_equal(FE.FEE_GAIN, 582.237*units.ohm)
    assert_equal(FE.DAQ_GAIN, 1.25)
    assert_equal(FE.NBITS, 12)
    assert_equal(FE.LSB, 2.0*units.V/2**FE.NBITS/FE.DAQ_GAIN)
    assert_equal(FE.NOISE_I, FE.LSB/(FE.FEE_GAIN*FE.DAQ_GAIN))
    assert_equal(FE.NOISE_DAQ, 0.313*units.mV)
    assert_equal(FE.C2, 8*units.nF)
    assert_equal(FE.C1, 2714*units.nF)
    assert_equal(FE.R1, 1567*units.ohm)
    assert_equal(FE.Zin, 62*units.ohm)
    assert_equal(FE.t_sample, 25*units.ns)
    assert_equal(FE.f_sample, 1./FE.t_sample)
    assert_equal(FE.f_mc, 1./(1*units.ns))
    assert_equal(FE.f_LPF1, 3*units.MHZ)
    assert_equal(FE.f_LPF2, 10*units.MHZ)
    assert_equal(FE.ADC_TO_PES, 20)  # nominal factor, comes out from spe area
    assert_equal(FE.OFFSET, 2500)  # offset adc

def test_spe():
    """
    Check the values of the FEE params
    """
    spe = FE.SPE()
    fee = FE.FEE(noise_FEEPMB_rms=0*units.mA, noise_DAQ_rms=0)
    assert_equal(FE.PMT_GAIN, 1.7e6)
    spe_i = FE.spe_pulse(spe,t0=100*units.ns, tmax=200*units.ns)
    spe_v = FE.signal_v_fee(fee, spe_i)
    spe_adc = FE.daq_decimator(1000.*units.MHZ, 40*units.MHZ, spe_v*FE.v_to_adc())
    adc_to_pes = np.sum(spe_adc)
    assert_greater(adc_to_pes, 18)
    assert_less(adc_to_pes, 22)
