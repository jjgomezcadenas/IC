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
    OFFSET = tables.Int16Col(pos=1)  # displaces the baseline (e.g, 700)
    CEILING = tables.Int16Col(pos=2)  # adc top count (4096)
    PMT_GAIN = tables.Float32Col(pos=3)  # Gain of PMT (4.5e6)
    FEE_GAIN = tables.Float32Col(pos=4)  # FE gain (250*ohm)
    R1 = tables.Float32Col(pos=5)  # resistor in Ohms (2350*ohm)
    C1 = tables.Float32Col(pos=6)  # Capacitor C1 in nF
    C2 = tables.Float32Col(pos=7)  # Capacitor C2 in nF
    ZIN = tables.Float32Col(pos=8)  # equivalent impedence
    DAQ_GAIN = tables.Float32Col(pos=9)
    NBITS = tables.Float32Col(pos=10)  # number of bits ADC
    LSB = tables.Float32Col(pos=11)  # LSB (adc count)
    NOISE_I = tables.Float32Col(pos=12)  # Noise at the input
    NOISE_DAQ = tables.Float32Col(pos=13) # Noise at DAQ
    t_sample = tables.Float32Col(pos=14) # sampling time
    f_sample = tables.Float32Col(pos=15) # sampling frequency
    f_mc = tables.Float32Col(pos=16) # sampling frequency in MC (1ns)
    f_LPF1 = tables.Float32Col(pos=17)  # LPF
    f_LPF2 = tables.Float32Col(pos=18)  # LPF
    coeff_c = tables.Float64Col(shape=12, pos=19)  # cleaning coeff
    coeff_blr = tables.Float64Col(shape=12, pos=20)  # COEFF BLR
    adc_to_pes = tables.Float32Col(shape=12, pos=21)  # CALIB CONST
    pmt_noise_rms = tables.Float32Col(shape=12, pos=22)  # rms noise


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
