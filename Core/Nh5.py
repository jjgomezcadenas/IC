"""
Tables defining the DM
"""
import tables

class RunInfo(tables.IsDescription):
    run_number = tables.Int32Col(shape=(), pos=0)


class EventInfo(tables.IsDescription):
    evt_number = tables.Int32Col(shape=(), pos=0)
    timestamp = tables.UInt64Col(shape=(), pos=1)


class DetectorGeometry(tables.IsDescription):
    """
    Stores geometry information for the detector
    """
    x_det = tables.Float32Col(shape=2, pos=1)  # xmin, xmax
    y_det = tables.Float32Col(shape=2, pos=2)  # ymin, ymax
    z_det = tables.Float32Col(shape=2, pos=3)  # zmin, zmax
    r_det = tables.Float32Col(pos=4)  # radius


class DataSensor(tables.IsDescription):
    """
    Stores metadata information for the SiPMs
    (position, gain, calibration-constant, mask)
    """
    channel = tables.Int32Col(pos=0)  # electronic channel
    position = tables.Float32Col(shape=3, pos=1)
    coeff = tables.Float64Col(pos=2)
    adc_to_pes = tables.Float32Col(pos=3)
    noise_rms = tables.Float32Col(pos=4)


class MCTrack(tables.IsDescription):
    """
    Stores the parameters used by the simulation as metadata
    using Pytables
    """
    event_indx = tables.Int16Col(pos=1)
    mctrk_indx = tables.Int16Col(pos=2)
    particle_name = tables.StringCol(10, pos=3)
    pdg_code = tables.Int16Col(pos=4)
    initial_vertex = tables.Float32Col(shape=3, pos=5)
    final_vertex = tables.Float32Col(shape=3, pos=6)
    momentum = tables.Float32Col(shape=3, pos=7)
    energy = tables.Float32Col(pos=8)
    nof_hits = tables.Int16Col(pos=9)
    hit_indx = tables.Int16Col(pos=10)
    hit_position = tables.Float32Col(shape=3, pos=11)
    hit_time = tables.Float32Col(pos=12)
    hit_energy = tables.Float32Col(pos=13)


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
    offset = tables.Int16Col(pos=1)  # displaces the baseline (e.g, 700)
    ceiling = tables.Int16Col(pos=2)  # adc top count (4096)
    pmt_gain = tables.Float32Col(pos=3)  # Gain of PMT (4.5e6)
    V_gain = tables.Float32Col(pos=4)  # FE gain (250*ohm)
    R = tables.Float32Col(pos=5)  # resistor in Ohms (2350*ohm)
    time_step = tables.Float32Col(pos=6)  # 1*ns input MC time bins
    time_daq = tables.Float32Col(pos=7)  # 25*ns DAQ time
    freq_LPF = tables.Float32Col(pos=8)  # 3E6*hertz
    freq_HPF = tables.Float32Col(pos=9)  # 1/2piRC
    LSB = tables.Float32Col(pos=10)  # Least Significant Bit 2*volt/2**NBITS
    volts_to_adc = tables.Float32Col(pos=11)  # volts to adc counts
    noise_fee_rms = tables.Float32Col(pos=12)  # noise FEE in volts
    noise_adc = tables.Float32Col(pos=13)   # noise FEE in ADC counts
    C12 = tables.Float32Col(shape=12, pos=14)  # 6.2*nF decoupling capacitor
    AC = tables.Float32Col(shape=12, pos=15)  # Accumulator coefficients
    CR = tables.Float32Col(shape=12, pos=16)  # calibration constants RAW
    CB = tables.Float32Col(shape=12, pos=17)  # calibration constants BLR


class PMAP(tables.IsDescription):
    """
    Store for a PMap
    """
    event = tables.Int32Col(pos=0)
    peak = tables.UInt8Col(pos=1)
    signal = tables.StringCol(2, pos=2)
    time = tables.Float32Col(pos=3)
    ToT = tables.UInt16Col(pos=4)
    cathode = tables.Float32Col(pos=5)
    anode = tables.Float32Col(pos=6, shape=(1792,))
