"""
ANASTASIA
GML October 2016

What ANASTASIA does:
1) Reads a hdf5 file containing the PMT's CWF and the SiPMs' RWF in ADC counts.
2) Subtracts the SiPMs' baseline.
3) Applies zero-suppression to both the big PMT and the individual SiPMs.
4) Writes the ZS waveforms in the same file as earrays.
"""

from __future__ import print_function

import sys
from time import time

import numpy as np
import tables as tb

from Core.LogConfig import logger
from Core.Configure import configure, define_event_loop

import Core.sensorFunctions as snf
import Core.wfmFunctions as wfm

from Core.RandomSampling import NoiseSampler as SiPMsNoiseSampler
"""

ANASTASIA
ChangeLog:

14.10 First version.

18.10 Big PMT and ZS methods implemented.

21.10 Several fixes. PRE-RELEASE.

31.10 Baseline subtraction for SiPMs introduced.

10.11 Waveforms stay in adc counts. All PMTs are now stored.
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
    FIRST_EVT = CFP["FIRST_EVT"]
    LAST_EVT = CFP["LAST_EVT"]
    RUN_ALL = CFP["RUN_ALL"]
    NEVENTS = LAST_EVT - FIRST_EVT

    # Increate thresholds by 1% for safety
    PMT_NOISE_CUT_RAW = CFP["PMT_NOISE_CUT_RAW"] * 1.01
    PMT_NOISE_CUT_BLR = CFP["PMT_NOISE_CUT_BLR"] * 1.01
    SIPM_ZS_METHOD = CFP["SIPM_ZS_METHOD"]
    SIPM_NOISE_CUT = CFP["SIPM_NOISE_CUT"]

    logger.info("Debug level = {}".format(DEBUG_LEVEL))
    logger.info("input file = {}/{}".format(PATH_IN, FILE_IN))
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
        sipmdata = h5in.root.Sensors.DataSiPM
        sipmdf = snf.read_data_sensors(sipmdata)

        NEVT, NPMT, PMTWL = pmtcwf.shape
        NEVT, NSIPM, SIPMWL = sipmrwf.shape

        logger.info("# events in DST = {}".format(NEVT))
        logger.info("#PMTs = {} #SiPMs = {}".format(NPMT, NSIPM))
        logger.info("PMT WFL = {} SiPM WFL = {}".format(PMTWL, SIPMWL))

        # Create instance of the noise sampler and compute noise thresholds
        sipms_noise_sampler_ = SiPMsNoiseSampler(SIPMWL)

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
                                     shape=(0, NPMT, PMTWL),
                                     expectedrows=NEVT)

        blr_zs_ = h5in.create_earray(h5in.root.ZS, "BLR",
                                     atom=tb.Int16Atom(),
                                     shape=(0, NPMT, PMTWL),
                                     expectedrows=NEVT)

        sipm_zs_ = h5in.create_earray(h5in.root.ZS, "SiPM",
                                      atom=tb.Int16Atom(),
                                      shape=(0, NSIPM, SIPMWL),
                                      expectedrows=NEVT)

        first_evt, last_evt, print_mod = define_event_loop(FIRST_EVT, LAST_EVT,
                                                           NEVENTS,
                                                           NEVT, RUN_ALL)

        t0 = time()
        for i in range(first_evt, last_evt):
            if not i%print_mod:
                logger.info("-->event number = {}".format(i))

            pmtzs = wfm.noise_suppression(pmtcwf[i], PMT_NOISE_CUT_RAW)
            blrzs = wfm.noise_suppression(pmtblr[i], PMT_NOISE_CUT_BLR)

            pmt_zs_.append(pmtzs[np.newaxis])
            blr_zs_.append(blrzs[np.newaxis])

            sipmzs = sipmrwf[i]
            if "/MC" not in h5in:
                sipmzs = wfm.subtract_baseline(sipmzs, None)
            sipmzs = wfm.noise_suppression(sipmzs, sipms_thresholds_)
            sipm_zs_.append(sipmzs[np.newaxis])

        t1 = time()
        dt = t1-t0

        print("ANASTASIA has run over {} events in {} seconds".format(i, dt))
    print("Leaving ANASTASIA. Safe travels!")


if __name__ == "__main__":
    from cities import anastasia
    print(anastasia)
    ANASTASIA(sys.argv)
