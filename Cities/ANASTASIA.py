"""
ANASTASIA
GML October 2016

What ANASTASIA does:
1) Reads a hdf5 file containing the PMT's CWF and the SiPMs' RWF in ADC counts.
2) Creates a single "big" PMT summing up PMTs' waveforms.
3) Subtracts the SiPMs' baseline.
4) Applies zero-suppression to both the big PMT and the individual SiPMs.
5) Converts the waveforms from adc to pes.
6) Writes the ZS waveforms in the same file as earrays.
"""

from __future__ import print_function

import sys
from time import time

import numpy as np
import tables as tb
import scipy as sc
import scipy.signal

from LogConfig import logger
from Configure import configure, define_event_loop

import sensorFunctions as snf
import wfmFunctions as wfm

from RandomSampling import NoiseSampler as SiPMsNoiseSampler
"""

ANASTASIA
ChangeLog:

14.10 First version.

18.10 Big PMT and ZS methods implemented.

21.10 Several fixes. PRE-RELEASE.

31.10 Baseline subtraction for SiPMs introduced.

"""


def ANASTASIA(argv):
    """
    ANASTASIA driver
    """
    DEBUG_LEVEL, INFO, CFP = configure(argv[0], argv[1:])

    if INFO:
        print(__doc__)

    PATH_IN = CFP["PATH_IN"]
    FILE_IN = CFP["FILE_IN"]
    PATH_DB = CFP["PATH_DB"]
    FIRST_EVT = CFP["FIRST_EVT"]
    LAST_EVT = CFP["LAST_EVT"]
    RUN_ALL = CFP["RUN_ALL"]
    NEVENTS = LAST_EVT - FIRST_EVT

    PMT_NOISE_CUT_RAW = CFP["PMT_NOISE_CUT_RAW"]
    PMT_NOISE_CUT_BLR = CFP["PMT_NOISE_CUT_BLR"]
    SIPM_ZS_METHOD = CFP["SIPM_ZS_METHOD"]
    SIPM_NOISE_CUT = CFP["SIPM_NOISE_CUT"]

    logger.info("Debug level = {}".format(DEBUG_LEVEL))
    logger.info("input file = {}/{}".format(PATH_IN, FILE_IN))
    logger.info("path to database = {}".format(PATH_DB))
    logger.info("First event = {} last event = {} "
                "# events requested = {}".format(FIRST_EVT, LAST_EVT, NEVENTS))
    logger.info("ZS method PMTS RAW = {}. "
                "Cut value = {}".format("RMS_CUT", PMT_NOISE_CUT_RAW))
    logger.info("ZS method PMTS BLR = {}. "
                "Cut value = {}".format("ABSOLUTE", PMT_NOISE_CUT_BLR))
    logger.info("ZS method SIPMS = {}. "
                "Cut value = {}".format(SIPM_ZS_METHOD, SIPM_NOISE_CUT))

    with tb.open_file("{}/{}".format(PATH_IN, FILE_IN), "r+") as h5in:
        pmtblr = h5in.root.RD.pmtblr
        pmtcwf = h5in.root.RD.pmtcwf
        sipmrwf = h5in.root.RD.sipmrwf
        pmtdata = h5in.root.Sensors.DataPMT
        blrdata = h5in.root.Sensors.DataBLR
        sipmdata = h5in.root.Sensors.DataSiPM
        pmtdfraw = snf.read_data_sensors(pmtdata)
        pmtdfblr = snf.read_data_sensors(blrdata)
        sipmdf = snf.read_data_sensors(sipmdata)

        NEVT, NPMT, PMTWL = pmtcwf.shape
        NEVT, NSIPM, SIPMWL = sipmrwf.shape

        logger.info("# events in DST = {}".format(NEVT))
        logger.info("#PMTs = {} #SiPMs = {}".format(NPMT, NSIPM))
        logger.info("PMT WFL = {} SiPM WFL = {}".format(PMTWL, SIPMWL))

        # Calibration constants and their average
        pmt_cal_consts_raw = abs(pmtdfraw["adc_to_pes"].reshape(NPMT, 1))
        pmt_cal_consts_blr = abs(pmtdfblr["adc_to_pes"].reshape(NPMT, 1))
        pmt_ave_consts_raw = np.mean(pmt_cal_consts_raw)
        pmt_ave_consts_blr = np.mean(pmt_cal_consts_blr)

        # FEE noise in ADC
        noise_adc = 0.789#h5in.root.MC.FEE.col("noise_adc")[0]

        # Create instance of the noise sampler and compute noise thresholds
        sipms_noise_sampler_ = SiPMsNoiseSampler(PATH_DB+"/NoiseSiPM_NEW.h5",
                                                 SIPMWL)

        # Increate thresholds by 1% for safety
        pmts_noise_threshold_raw_ = (PMT_NOISE_CUT_RAW * NPMT /
                                     pmt_ave_consts_raw * 1.01)
        pmts_noise_threshold_blr_ = (PMT_NOISE_CUT_BLR * noise_adc /
                                     pmt_ave_consts_blr * NPMT**0.5 * 1.01)

        if SIPM_ZS_METHOD == "FRACTION":
            sipms_thresholds_ = sipms_noise_sampler_.ComputeThresholds(
                                SIPM_NOISE_CUT, sipmdf['adc_to_pes'])
        else:
            sipms_thresholds_ = np.ones(NSIPM) * SIPM_NOISE_CUT

        if "/ZS" not in h5in:
            h5in.create_group(h5in.root, "ZS")
        if "/ZS/PMT" in h5in:
            h5in.remove_node("/ZS", "PMT")
        if "/ZS/BLR" in h5in:
            h5in.remove_node("/ZS", "BLR")
        if "/ZS/SiPM" in h5in:
            h5in.remove_node("/ZS", "SiPM")

        # Notice the Int16, not Float32! bad for compression
        pmt_zs_ = h5in.create_earray(h5in.root.ZS, "PMT",
                                     atom=tb.Int16Atom(),
                                     shape=(0, 1, PMTWL),
                                     expectedrows=NEVT)

        pmt_zs_blr_ = h5in.create_earray(h5in.root.ZS, "BLR",
                                         atom=tb.Int16Atom(),
                                         shape=(0, 1, PMTWL),
                                         expectedrows=NEVT)

        sipm_zs_ = h5in.create_earray(h5in.root.ZS, "SiPM",
                                      atom=tb.Int16Atom(),
                                      shape=(0, NSIPM, SIPMWL),
                                      expectedrows=NEVT)

        first_evt, last_evt = define_event_loop(FIRST_EVT, LAST_EVT,
                                                NEVENTS, NEVT, RUN_ALL)

        t0 = time()
        for i in range(first_evt, last_evt):
            logger.info("-->event number ={}".format(i))

            # Integrate PMT plane in pes (not in time!)
            pmtcwf_int_pes = (pmtcwf[i] / pmt_cal_consts_raw).sum(axis=0)
            pmtblr_int_pes = (pmtblr[i] / pmt_cal_consts_blr).sum(axis=0)

            # suppress_wf puts zeros where the wf is below the threshold
            pmtcwf_int_pes = wfm.suppress_wf(pmtcwf_int_pes,
                                             pmts_noise_threshold_raw_)
            pmtblr_int_pes = wfm.suppress_wf(pmtblr_int_pes,
                                             pmts_noise_threshold_blr_)

            pmt_zs_.append(pmtcwf_int_pes.reshape(1, 1, PMTWL))
            pmt_zs_blr_.append(pmtblr_int_pes.reshape(1, 1, PMTWL))

            SiPMdata = wfm.subtract_baseline(sipmrwf[i], None)
            SiPMdata = wfm.noise_suppression(SiPMdata, sipms_thresholds_)
            #SiPMdata = wfm.to_pes(SiPMdata, sipmdf)

            sipm_zs_.append(SiPMdata.reshape(1, NSIPM, SIPMWL))
        t1 = time()
        dt = t1-t0

        print("ANASTASIA has run over {} events in {} seconds".format(i, dt))
    print("Leaving ANASTASIA. Safe travels!")


if __name__ == "__main__":
    from cities import anastasia
    print(anastasia)
    ANASTASIA(sys.argv)
