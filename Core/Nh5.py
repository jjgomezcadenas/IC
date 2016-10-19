"""
Tables defining the DM
"""
import tables
class DetectorGeometry(tables.IsDescription):
    """
    Stores geometry information for the detector
    """
    x_det = tables.Float32Col(shape=2, pos=1) #xmin, xmax
    y_det = tables.Float32Col(shape=2, pos=2) #ymin, ymax
    z_det = tables.Float32Col(shape=2, pos=3) #zmin, zmax
    r_det = tables.Float32Col(pos=4) # radius


class DataPMT(tables.IsDescription):
    """
    Stores metadata information for the PMTs
    (position, gain, calibration-constant, mask)
    """
    channel = tables.Int16Col(pos=1) #electronic channel
    active = tables.Int16Col(pos=2) # 1 if active. 0 if dead
    position = tables.Float32Col(shape=3, pos=3)
    gain =tables.Float32Col(pos=4)   #for PMT gain is the accumulator coeff
    adc_to_pes =tables.Float32Col(pos=5)

class DataSiPM(tables.IsDescription):
    """
    Stores metadata information for the SiPMs
    (position, gain, calibration-constant, mask)
    """
    channel = tables.Int16Col(pos=1) #electronic channel
    active = tables.Int16Col(pos=2) # 1 if active. 0 if dead
    position = tables.Float32Col(shape=3, pos=3)
    gain =tables.Float32Col(pos=4)
    adc_to_pes =tables.Float32Col(pos=5)

class MCTrack(tables.IsDescription):
    """
    Stores the parameters used by the simulation as metadata
    using Pytables
    """
    event_indx = tables.Int16Col(pos=1)
    mctrk_indx = tables.Int16Col(pos=2)
    particle_name = tables.StringCol(10,pos=3)  #displaces the baseline (e.g, 700)
    pdg_code = tables.Int16Col(pos=4)   # number of PMTs (12)
    initial_vertex =tables.Float32Col(shape=3, pos=5)
    final_vertex =tables.Float32Col(shape=3, pos=6)
    momentum =tables.Float32Col(shape=3, pos=7)
    energy =tables.Float32Col(pos=8)
    nof_hits = tables.Int16Col(pos=9)
    hit_indx = tables.Int16Col(pos=10)
    hit_position = tables.Float32Col(shape=3, pos=11)
    hit_time =tables.Float32Col(pos=12)
    hit_energy =tables.Float32Col(pos=13)

class SENSOR_WF(tables.IsDescription):
    """
    Describes a true waveform (zero supressed)
    """
    event = tables.UInt32Col(pos=0)
    ID = tables.UInt32Col(pos=1)
    time_mus = tables.Float32Col(pos=2)
    ene_pes = tables.Float32Col(pos=3)

class FEE(tables.IsDescription):
    """
    Stores the parameters used by the EP simulation as metadata
    """

    offset = tables.Int16Col(pos=1)  #displaces the baseline (e.g, 700)
    pmt_gain =tables.Float32Col(pos=2)  #Gain of PMT (4.5e6)
    V_gain =tables.Float32Col(pos=3)  #FE gain (250*ohm)
    R = tables.Float32Col(pos=4) # resistor in Ohms (2350*ohm)
    C12 = tables.Float32Col(shape=12,pos=5) #6.2*nF  decoupling capacitor in nF
    CR = tables.Float32Col(shape=12,pos=6) # calibration constants RAW
    CB = tables.Float32Col(shape=12,pos=7) # calibration constants BLR
    AC = tables.Float32Col(shape=12,pos=8) #Accumulator coefficients
    time_step=tables.Float32Col(pos=9) #1*ns input MC time bins
    time_daq=tables.Float32Col(pos=10) #25*ns DAQ time
    freq_LPF=tables.Float32Col(pos=11) #3E6*hertz
    freq_HPF=tables.Float32Col(pos=12) #1/2piRC
    LSB = tables.Float32Col(pos=13)    # Least Significant Bit 2*volt/2**NBITS,
    volts_to_adc = tables.Float32Col(pos=14) # conversion from volts to adc counts
    noise_fee_rms = tables.Float32Col(pos=15) # noise FEE in volts
    noise_adc = tables.Float32Col(pos=16) # noise FEE in ADC counts
