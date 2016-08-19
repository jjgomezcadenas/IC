"""
Tables defining the DM
"""
import tables
class DetectorGeometry(tables.IsDescription):
    """
    Stores geometry information for the detector
    """
    x_det = tables.Float64Col(shape=2, pos=1) #xmin, xmax
    y_det = tables.Float64Col(shape=2, pos=2) #ymin, ymax
    z_det = tables.Float64Col(shape=2, pos=3) #zmin, zmax
    r_det = tables.Float64Col(pos=4) # radius


class DataPMT(tables.IsDescription):
    """
    Stores metadata information for the PMTs
    (position, gain, calibration-constant, mask)
    """
    channel = tables.Int16Col(pos=1) #electronic channel
    active = tables.Int16Col(pos=2) # 1 if active. 0 if dead
    position = tables.Float64Col(shape=3, pos=3)
    gain =tables.Float64Col(pos=4)
    adc_to_pes =tables.Float64Col(pos=5)

class DataSiPM(tables.IsDescription):
    """
    Stores metadata information for the SiPMs
    (position, gain, calibration-constant, mask)
    """
    channel = tables.Int16Col(pos=1) #electronic channel
    active = tables.Int16Col(pos=2) # 1 if active. 0 if dead
    position = tables.Float64Col(shape=3, pos=3)
    gain =tables.Float64Col(pos=4)
    adc_to_pes =tables.Float64Col(pos=5)

class MCTrack(tables.IsDescription):
    """
    Stores the parameters used by the simulation as metadata
    using Pytables
    """
    event_indx = tables.Int16Col(pos=1) 
    mctrk_indx = tables.Int16Col(pos=2) 
    particle_name = tables.StringCol(10,pos=3)  #displaces the baseline (e.g, 700)
    pdg_code = tables.Int16Col(pos=4)   # number of PMTs (12) 
    initial_vertex =tables.Float64Col(shape=3, pos=5)
    final_vertex =tables.Float64Col(shape=3, pos=6)
    momentum =tables.Float64Col(shape=3, pos=7)
    energy =tables.Float64Col(pos=8)
    nof_hits = tables.Int16Col(pos=9) 
    hit_indx = tables.Int16Col(pos=10)
    hit_position = tables.Float64Col(shape=3, pos=11)
    hit_time =tables.Float64Col(pos=12)
    hit_energy =tables.Float64Col(pos=13)


class FEE(tables.IsDescription):
    """
    Stores the parameters used by the EP simulation as metadata
    """
    
    offset = tables.Int16Col(pos=1)  #displaces the baseline (e.g, 700)
    pmt_gain =tables.Float32Col(pos=2)  #Gain of PMT (4.5e6)
    V_gain =tables.Float32Col(pos=3)  #FE gain (250*ohm)
    R = tables.Float32Col(pos=4) # resistor in Ohms (2350*ohm)
    C12 = tables.Float32Col(shape=12,pos=5) #6.2*nF  decoupling capacitor in pF
    CO12 = tables.Float32Col(shape=12,pos=6) #Accumulator coefficients
    time_step=tables.Float32Col(pos=7) #1*ns input MC time bins
    time_daq=tables.Float32Col(pos=8) #25*ns DAQ time 
    freq_LPF=tables.Float32Col(pos=9) #3E6*hertz
    freq_HPF=tables.Float32Col(pos=10) #1/2piRC
    LSB = tables.Float32Col(pos=11)    # Least Significant Bit 2*volt/2**NBITS, 
    volts_to_adc = tables.Float32Col(pos=12) # conversion from volts to adc counts
    noise_fee_rms = tables.Float32Col(pos=13) # noise FEE in volts
    noise_adc = tables.Float32Col(pos=14) # noise FEE in ADC counts

