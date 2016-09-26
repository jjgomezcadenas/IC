"""
Sensors Response
JJGC August 2016


Simulates the response of the energy plane (HPF, LPF, noise) and the response of the
SiPM plane (not yet implemented)
 
"""

from __future__ import print_function
from Util import *
from LogConfig import *
import FEParam as FP
import SPE as SP
import FEE2 as FE

#import tables
from FEE2 import down_scale_signal_


"""
Code
"""

def FEE_param_table(fee_table):
    """
    Stores the parameters of the EP FEE simulation 
    """
    row = fee_table.row
    row['offset'] = FP.offset
    row['pmt_gain'] = FP.PMT_GAIN
    row['V_gain'] = FP.V_GAIN
    row['R'] = FP.R
    row['C12'] = FP.C12
    row['CO12'] = FP.C12 # to be rewritten by ISIDORA
    row['time_step'] = FP.time_step
    row['time_daq'] = FP.time_DAQ
    row['freq_LPF'] = FP.freq_LPF
    row['freq_HPF'] = 1./(2*pi*FP.R*FP.C)
    row['LSB'] = FP.LSB
    row['volts_to_adc'] = FP.voltsToAdc/volt
    row['noise_fee_rms'] = FP.NOISE_FEE
    row['noise_adc'] = FP.NOISE_ADC
    
    row.append()
    

def energy_pes(event_number, sensord):
    """
    Sum the WFs of PMTs and SiPMs (MC) and store the total energy in PES
    """     
    rdata = []

    for j in range(sensord.shape[1]):
        swf = sensord[event_number, j]
        ene = np.sum(swf)
        rdata.append(ene)
        
    return np.array(rdata) 

def simulate_sipm_response(event_number,sipmrd_):
    """
    For the moment use a dummy rutne that simply copies the sipm EARRAY
    """
    rdata = []

    for j in range(sipmrd_.shape[1]):
        rdata.append(sipmrd_[event_number, j])
    return np.array(rdata)


def simulate_pmt_response(event_number,pmtrd_):
    """
    Sensor Response
    Given a signal in PE (photoelectrons in bins of 1 ns) and the response function of 
    for a single photoelectron (spe) and front-end electronics (fee)
    this function produces the PMT raw data (adc counts bins 25 ns)

    pmtrd_ dataset that holds the PMT PE data for each PMT
    pmtrd25 dataset to be created with adc counts, bins 25 ns 
    after convoluting with electronics
    """
  
    rdata = []

    for j in range(pmtrd_.shape[1]):
        logger.debug("-->PMT number ={}".format(j))
                
        pmt = pmtrd_[event_number, j] #waveform for event event_number, PMT j
        
        fee = FE.FEE(C=FP.C12[j],R= FP.R, f=FP.freq_LPF, RG=FP.V_GAIN) 
        spe = SP.SPE(pmt_gain=FP.PMT_GAIN,x_slope = 5*ns,x_flat = 1*ns)
    
        signal_PMT = spe.SpePulseFromVectorPE(pmt) #PMT response

        #Front end response to PMT pulse (in volts)
        signal_fee = fee.FEESignal(signal_PMT, noise_rms=FP.NOISE_FEE) 

        #Signal out of DAQ
        #positive signal convention
        #signal_daq = fee.daqSignal(signal_fee, noise_rms=0) - FP.offset
        #negative signals convention!

        signal_daq = FP.offset -fee.daqSignal(signal_fee, noise_rms=0) 
    
        rdata.append(signal_daq)
    return np.array(rdata)

def decimate_signal(event_number,pmtrd_):
    """
    Decimates the MCRD signal to produce TWF (pes, bins 25 ns)
    """
  
    rdata = []

    for j in range(pmtrd_.shape[1]):
        logger.debug("-->PMT number ={}".format(j))
                
        pmt = pmtrd_[event_number, j] #waveform for event event_number, PMT j
        twf = down_scale_signal_(pmt, int(FP.time_DAQ))
        
        
        rdata.append(twf)
    return np.array(rdata)

def rebin_pmt_array(a,stride):
    """
    rebins pmt array a according to stride
    there is a rebin_array in util which uses
    lenb = (len(a))/int(stride)
    this version uses (lean(a)+1) to correct from the fact that the
    MCRD is 599999 channels (should be 600000)
    """
    
    lenb = (len(a)+1)/int(stride)
    b = np.zeros(lenb)
    j=0
    for i in range(lenb):
        b[i] = np.sum(a[j:j+stride])
        j+= stride
    return b

def rebin_signal(event_number,pmtrd_, stride):
    """
    rebins the MCRD signal to produce TWF (pes, bins 25 ns)
    """
    
    rdata = []

    for j in range(pmtrd_.shape[1]):
        logger.debug("-->PMT number ={}".format(j))
                
        pmt = pmtrd_[event_number, j] #waveform for event event_number, PMT j
        twf = rebin_pmt_array(pmt, stride)
        
        rdata.append(twf)
    return np.array(rdata)

  




