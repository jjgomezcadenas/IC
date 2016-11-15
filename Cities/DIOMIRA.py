"""
DIOMIRA
JJGC August-October 2016
GML October 2016

What DIOMIRA does:
1) Reads a MCRD file containing MC waveforms for the 12 PMTs of the EP.
   Each waveform contains number of PEs in bins of 1 ns.
2) Convolves the PE waveform with the response of the FEE electronics.
3) Decimates the waveform, simulating the effect of the DAQ sampling
(25 ns bins)
4) Writes a RWF file with the new data and adds the FEE simulation parameters
as metadata
"""

from __future__ import print_function
import sys
import numpy as np
import tables
from scipy import signal as SGN
from time import time

import system_of_units as units
from LogConfig import logger
from Configure import configure, define_event_loop
from Nh5 import FEE, SENSOR_WF

import FEParam as FP
import SPE as SP
import FEE2 as FE

import wfmFunctions as wfm
import coreFunctions as cf
import tblFunctions as tbl
import sensorFunctions as snf

from RandomSampling import NoiseSampler as SiPMsNoiseSampler
"""

DIOMIRA
ChangeLog:

26.9

Changed types of PMTRWF, SIPMRWF and PMTTWF to Float32 for
    (compatibility with ART/GATE)

Do not store EPMT and ESIPM (can be computed on the fly)

Change sign of pmtrwf to negative (as produced by the DAQ)

28.9 add cython

29.9 changed the way true waveforms are treated.
before: --> full waveform (in bins of 25 ns) in an EArray
now--> ZS waveform rebinned at 1 mus in a Table
(adavantages:) faster processing less space

01.10 moved loop definition to Configure.py and added index to TWF table

11.10 introduced SiPM noise. SiPMs' true waveforms stored under TWF group.
Some variables, classes and functions renamed for clarity.

12.10 ZS functions to store the SiPMs

13.10 Reutilization of functions and some duplicities removed.
      SiPMs' waveforms to be stored without ZS.

14.10 JJ: Store all waveforms as Int16 (and as adc counts)

17.10: JJ PEP8 compliant and pyflakes compliant
(adjust to PEP8 syntax, eliminate import* etc.)

18.10 GML, add soft cut to eliminate noise below 0.5 pes in sipm plane

19.10 JJ: write calibration constants to FEE table!

20.10: JJ, store BLR with positive signal and baseline subtracted

20.10: GML, overwrite calibration constants in DataPMT with values
from FEE table. PRE-RELEASE

"""


def FEE_param_table(fee_table):
    """
    Stores the parameters of the EP FEE simulation
    """
    row = fee_table.row
    row["offset"] = FP.offset
    row["ceiling"] = FP.ceiling
    row["pmt_gain"] = FP.PMT_GAIN
    row["V_gain"] = FP.V_GAIN
    row["R"] = FP.R
    row["C12"] = FP.C12
    row["CR"], row["CB"] = calibration_constants_from_spe()
    row["AC"] = FP.AC
    row["time_step"] = FP.time_step
    row["time_daq"] = FP.time_DAQ
    row["freq_LPF"] = FP.freq_LPF
    row["freq_HPF"] = 1./(2*np.pi*FP.R*FP.C)
    row["LSB"] = FP.LSB
    row["volts_to_adc"] = FP.voltsToAdc/units.volt
    row["noise_fee_rms"] = FP.NOISE_FEE
    row["noise_adc"] = FP.NOISE_ADC

    row.append()
    fee_table.flush()


def save_pmt_cal_consts(pmt_table, consts):
    """
    Overwrite PMT cal constats in table.
    """
    adc_to_pes = pmt_table.cols.adc_to_pes
    for i, const in enumerate(consts):
        adc_to_pes[i] = const


def simulate_sipm_response(event_number, sipmrd_, sipms_noise_sampler):
    """
    Add noise with the NoiseSampler class and return the noisy waveform.
    """
    return sipmrd_[event_number] + sipms_noise_sampler.Sample()


def simulate_pmt_response(event_number, pmtrd_, blr_mau=500):
    """
    Input:
     1) extensible array pmtrd_ (events, sensors, waveform)
     2) event_number

    returns:
    array of raw waveforms (RWF), obtained by convoluting pmtrd_ with the PMT
    front end electronics (LPF, HPF)
    array of BLR waveforms (only convolution with LPF)
    """

    RWF = []
    BLRX = []

    for j in range(pmtrd_.shape[1]):
        logger.debug("-->PMT number ={}".format(j))

        pmt = pmtrd_[event_number, j]  # waveform for event event_number, PMT j
        fee = FE.FEE(PMTG=FP.PMT_GAIN, C=FP.C12[j], R=FP.R, f=FP.freq_LPF,
                     RG=FP.V_GAIN)  # instantiate FEE class
        # instantiate single photoelectron class
        spe = SP.SPE(pmt_gain=FP.PMT_GAIN, x_slope=5*units.ns,
                     x_flat=1*units.ns)

        # waveform "pmt" is passed to spe, output is a signal current
        signal_PMT = spe.SpePulseFromVectorPE(pmt)  # PMT response

        # Front end response to PMT pulse (in volts)
        signal_fee, signal_blr = fee.FEESignal(signal_PMT, FP.NOISE_FEE)

        # daq response (decimation)
        signal_daq = FP.offset - fee.daqSignal(signal_fee, noise_rms=0)

        signal_daq_blr = (FP.ceiling - FP.offset +
                          fee.daqSignal(signal_blr, noise_rms=0))
        nm = blr_mau
        MAU = np.zeros(nm, dtype=np.double)
        B_MAU = (1./nm)*np.ones(nm, dtype=np.double)

        MAU[0:nm] = SGN.lfilter(B_MAU, 1, signal_daq_blr[0:nm])
        BASELINE = MAU[nm-1]

        RWF.append(signal_daq)
        BLRX.append(signal_daq_blr - BASELINE)

    return np.array(RWF), np.array(BLRX)


def calibration_constants_from_spe(start_pulse=100*units.ns,
                                   end_pulse=500*units.ns):
    """
    Computes calibration constants from the are of a SPE
    """
    spe = SP.SPE()
    cr = []
    cb = []

    for pmt, C in enumerate(FP.C12):
        fee = FE.FEE(PMTG=FP.PMT_GAIN, C=C, R=FP.R, f=FP.freq_LPF,
                     RG=FP.V_GAIN)
        # PMT response to a single photon (single pe current pulse)
        signal_t, signal_PE = spe.SpePulse(start_pulse, tmax=end_pulse)
        # effect of FEE
        signal_fee, signal_blr = fee.FEESignal(signal_PE, FP.NOISE_FEE)
        # effect of DAQ
        signal_daq = fee.daqSignal(signal_fee, noise_rms=0)
        signal_daq_blr = fee.daqSignal(signal_blr, noise_rms=0)
        area = np.sum(signal_daq)
        area_blr = np.sum(signal_daq_blr)
        logger.debug("PMT {}: cc {}, cc blr {}".format(pmt, area, area_blr))
        cr.append(area)
        cb.append(area_blr)

    return cr, cb


def DIOMIRA(argv):
    """
    Diomira driver
    """
    DEBUG_LEVEL, INFO, CFP = configure(argv[0], argv[1:])

    if INFO:

        print("""
        DIOMIRA:
         1. Reads a MCRD file produced by art/centella, which stores MCRD
        waveforms for PMTs (bins of 1 ns) and SiPMs (bins of 1 mus)
        2. Simulates the response of the energy plane and outputs both RWF
        and TWF
        3. Simulates the response of the tracking plane in the SiPMs and
        outputs SiPM RWF
        4. Add a table describing the FEE parameters used for simulation
        5. Copies the tables on geometry, detector data and MC
        """)
        FP.print_FEE()

    PATH_IN = CFP["PATH_IN"]
    PATH_OUT = CFP["PATH_OUT"]
    FILE_IN = CFP["FILE_IN"]
    FILE_OUT = CFP["FILE_OUT"]
    PATH_DB = CFP["PATH_DB"]
    FIRST_EVT = CFP["FIRST_EVT"]
    LAST_EVT = CFP["LAST_EVT"]
    RUN_ALL = CFP["RUN_ALL"]
    COMPRESSION = CFP["COMPRESSION"]
    NOISE_CUT = CFP["NOISE_CUT"]
    NEVENTS = LAST_EVT - FIRST_EVT

    logger.info("Debug level = {}".format(DEBUG_LEVEL))
    logger.info("Input path = {}; output path = {}".format(PATH_IN, PATH_OUT))
    logger.info("File_in = {} file_out = {}".format(FILE_IN, FILE_OUT))
    logger.info("Path to database = {}".format(PATH_DB))
    logger.info("First event = {} last event = {} "
                "# events requested = {}".format(FIRST_EVT, LAST_EVT, NEVENTS))
    logger.info("Compression library/level = {}".format(COMPRESSION))
    logger.info("Noise cut = {} pes ".format(NOISE_CUT))

    # open the input file
    with tables.open_file("{}/{}".format(PATH_IN, FILE_IN), "r") as h5in:
        # access the PMT raw data in file
        pmtrd_ = h5in.root.pmtrd
        sipmrd_ = h5in.root.sipmrd
        # pmtrd_.shape = (nof_events, nof_sensors, wf_length)

        NPMT = pmtrd_.shape[1]
        NSIPM = sipmrd_.shape[1]
        PMTWL = pmtrd_.shape[2]
        # PMTWL_FEE = int((PMTWL+1)/FP.time_DAQ) #old format
        PMTWL_FEE = int(PMTWL/FP.time_DAQ)
        SIPMWL = sipmrd_.shape[2]
        NEVENTS_DST = pmtrd_.shape[0]

        logger.info("nof PMTs = {} nof  SiPMs = {} "
                    "nof events in input DST = {} ".format(NPMT, NSIPM,
                                                           NEVENTS_DST))
        logger.info("lof PMT WF = {} lof SiPM WF (MC) = {} "
                    "lof PMT WF (FEE) = {}".format(PMTWL, SIPMWL, PMTWL_FEE))

        # access the geometry and the sensors metadata info
        geom_t = h5in.root.Detector.DetectorGeometry
        pmt_t = h5in.root.Sensors.DataPMT
        blr_t = h5in.root.Sensors.DataBLR
        sipm_t = h5in.root.Sensors.DataSiPM
        mctrk_t = h5in.root.MC.MCTracks
        sipmdf = snf.read_data_sensors(sipm_t)

        # Create instance of the noise sampler
        noise_sampler_ = SiPMsNoiseSampler(PATH_DB+"/NoiseSiPM_NEW.h5",
                                           SIPMWL, True)
        sipms_thresholds_ = NOISE_CUT * np.array(sipmdf["adc_to_pes"])

        # open the output file
        with tables.open_file("{}/{}".format(PATH_OUT, FILE_OUT), "w",
                              filters=tbl.filters(COMPRESSION)) as h5out:

            # create a group to store MC data
            mcgroup = h5out.create_group(h5out.root, "MC")
            # copy the mctrk table
            mctrk_t.copy(newparent=mcgroup)

            # create a group  to store geom data
            detgroup = h5out.create_group(h5out.root, "Detector")
            # copy the geom table
            geom_t.copy(newparent=detgroup)

            # create a group  store sensor data
            sgroup = h5out.create_group(h5out.root, "Sensors")
            # copy the pmt table
            pmt_t.copy(newparent=sgroup)
            blr_t.copy(newparent=sgroup)
            # copy the sipm table
            sipm_t.copy(newparent=sgroup)

            # create a table to store Energy plane FEE, hang it from MC group
            fee_table = h5out.create_table(mcgroup, "FEE", FEE,
                                           "EP-FEE parameters",
                                           tbl.filters("NOCOMPR"))

            # create a group to store True waveform data
            twfgroup = h5out.create_group(h5out.root, "TWF")
            # create a table to store true waveform (zs, rebinned)
            pmt_twf_table = h5out.create_table(twfgroup, "PMT", SENSOR_WF,
                                               "Store for PMTs TWF",
                                               tbl.filters(COMPRESSION))

            sipm_twf_table = h5out.create_table(twfgroup, "SiPM", SENSOR_WF,
                                                "Store for SiPM TWF",
                                                tbl.filters(COMPRESSION))

            # and index in event column
            pmt_twf_table.cols.event.create_index()
            sipm_twf_table.cols.event.create_index()

            # fill FEE table
            FEE_param_table(fee_table)
            pmt_t_copy = h5out.root.Sensors.DataPMT
            blr_t_copy = h5out.root.Sensors.DataBLR
            save_pmt_cal_consts(pmt_t_copy, fee_table.cols.CR[0])
            save_pmt_cal_consts(blr_t_copy, fee_table.cols.CB[0])

            # create a group to store RawData
            h5out.create_group(h5out.root, "RD")

            # create an extensible array to store the RWF waveforms
            pmtrwf = h5out.create_earray(h5out.root.RD, "pmtrwf",
                                         atom=tables.Int16Atom(),
                                         shape=(0, NPMT, PMTWL_FEE),
                                         expectedrows=NEVENTS_DST)

            pmtblr = h5out.create_earray(h5out.root.RD, "pmtblr",
                                         atom=tables.Int16Atom(),
                                         shape=(0, NPMT, PMTWL_FEE),
                                         expectedrows=NEVENTS_DST)

            sipmrwf = h5out.create_earray(h5out.root.RD, "sipmrwf",
                                          atom=tables.Int16Atom(),
                                          shape=(0, NSIPM, SIPMWL),
                                          expectedrows=NEVENTS_DST)
            # LOOP
            first_evt, last_evt, prind_mod = define_event_loop(FIRST_EVT,
                                                               LAST_EVT,
                                                               NEVENTS,
                                                               NEVENTS_DST,
                                                               RUN_ALL)
            t0 = time()
            for i in range(first_evt, last_evt):
                if not i%print_mod:
                    logger.info("-->event number = {}".format(i))

                # supress zeros in MCRD and rebin the ZS function in 1 mus bins
                rebin = int(units.mus/units.ns)

                trueSiPM = wfm.zero_suppression(sipmrd_[i], 0.)

                # dict_map applies a function to the dictionary values
                truePMT = cf.dict_map(lambda df: wfm.rebin_df(df, rebin),
                                      wfm.zero_suppression(pmtrd_[i],
                                      0., to_mus=int(units.ns/units.ms)))

                # store in table
                tbl.store_wf(i, pmt_twf_table, truePMT)
                tbl.store_wf(i, sipm_twf_table, trueSiPM)

                # simulate PMT response and return an array with RWF;BLR
                # convert to float, append to EVector

                dataPMT, blrPMT = simulate_pmt_response(i, pmtrd_)
                pmtrwf.append(dataPMT.astype(int).reshape(1, NPMT, PMTWL_FEE))

                pmtblr.append(blrPMT.astype(int).reshape(1,
                                                         NPMT,
                                                         PMTWL_FEE))

                # simulate SiPM response and return an array with RWF
                # convert to float, zero suppress and dump to table
                dataSiPM = simulate_sipm_response(i, sipmrd_, noise_sampler_)
                dataSiPM = wfm.to_adc(dataSiPM, sipmdf)
                dataSiPM = wfm.noise_suppression(dataSiPM, sipms_thresholds_)

                sipmrwf.append(dataSiPM.astype(int).reshape(1, NSIPM, SIPMWL))

            pmtrwf.flush()
            sipmrwf.flush()
            pmtblr.flush()

            t1 = time()
            dt = t1 - t0
            print("DIOMIRA has run over {} events in {} seconds".format(i+1,
                                                                        dt))
    print("Leaving Diomira. Safe travels!")

if __name__ == "__main__":
    from cities import diomira
    print(diomira)
    DIOMIRA(sys.argv)
