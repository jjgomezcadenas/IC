"""
ISIDORA
JJGC Agusut 2016

What ISIDORA does:
1) Reads a RWF file written by DIOMIRA
2) Performs DBLR
3) Write the corrected waveforms (CWF) to the file as new Evectors.

Change Log:
18/10 JJG added cython module. PRE-RELEASE

11.11 JJGC: A major refactoring of the code, now based in a much improved
deconv algorithm
"""

from __future__ import print_function

import sys
from time import time
import numpy as np
import tables

import Core.system_of_units as units
from Core.LogConfig import logger
from Core.Configure import configure, define_event_loop

import Sierpe.FEE as FE
import ICython.cBLR as cblr
import Database.loadDB as DB


def DBLR(pmtrwf, event_number, n_baseline=500, thr_trigger=5,
         thr_acum=2000,
         acum_discharge_length = 5000,
         acum_tau=2500,
         acum_compress=0.01):
    """
    Peform Base line Restoration
    """
    DataPMT = DB.DataPMT()
    NPMT = pmtrwf.shape[1]
    CWF = {}
    ACUM = {}

    for pmt in range(NPMT):
        signal_r, acum = cblr.\
          deconvolve_signal_acum(pmtrwf[event,pmt],
                                 n_baseline=n_baseline,
                                 coef_clean=DataPMT.coeff_c[pmt],
                                 coef_blr=DataPMT.coeff_blr[pmt],
                                 thr_trigger=thr_trigger,
                                 thr_acum=thr_acum,
                                 acum_discharge_length=acum_discharge_length,
                                 acum_tau=acum_tau,
                                 acum_compress=acum_compress)
        CWF[pmt] = signal_r
        ACUM[pmt] = acum
        return CWF, ACUM


def ISIDORA(argv):
    DEBUG_LEVEL, INFO, CFP = configure(argv[0], argv[1:])

    if INFO:

        print("""
        ISIDORA:
        1. Reads an Nh5 file produced by DIOMIRA, which stores the
            raw waveforms (RWF) for the PMTs and SiPMs waveforms, as well as
            data on geometry, sensors and MC. The RDWF of the PMTs
            show negative swing due to the HPF of the EP FEE electronics

        2. Performs DBLR on the PMT RWF and produces corrected waveforms (CWF).

        3. Adds the CWF and ancilliary info to the DST

        4. Computes the energy of the PMTs per each event and writes to DST

        """)

    PATH_IN = CFP["PATH_IN"]
    FILE_IN = CFP["FILE_IN"]
    FIRST_EVT = CFP["FIRST_EVT"]
    LAST_EVT = CFP["LAST_EVT"]
    RUN_ALL = CFP["RUN_ALL"]
    N_BASELINE = CFP["N_BASELINE"]
    THR_TRIGGER = CFP["THR_TRIGGER"]
    THR_ACUM = CFP["THR_ACUM"]
    ACUM_DISCHARGE_LENGTH = CFP["ACUM_DISCHARGE_LENGTH"]
    ACUM_TAU = CFP["ACUM_TAU"]
    ACUM_COMPRESS = CFP["ACUM_COMPRESS"]
    NEVENTS = LAST_EVT - FIRST_EVT

    logger.info("Debug level = {}".format(DEBUG_LEVEL))
    logger.info("Input path ={}; file_in ={} ".format(PATH_IN, FILE_IN))
    logger.info("First event = {} last event = {} "
                "# events requested = {}".format(FIRST_EVT, LAST_EVT, NEVENTS))
    logger.info("Baseline calculation length = {}"
                "n_sigma for trigger = {}".format(N_BASELINE, THR_TRIGGER)

    logger.info("""Accumulator Parameters:
                accumulator threshold = {}
                length for discharge = {}
                tau for discharge = {}
                compression factor = {}
                """.format(THR_ACUM, ACUM_DISCHARGE_LENGTH,
                           ACUM_TAU, ACUM_COMPRESS))

    # open the input file in mode append
    with tables.open_file("{}/{}".format(PATH_IN, FILE_IN), "a") as h5in:
        # access the PMT raw data in file
        pmtrd_ = h5in.root.RD.pmtrwf     # PMT raw data must exist

        # pmtrd_.shape = (nof_events, nof_sensors, wf_length)
        NEVENTS_DST, NPMT, PMTWL = pmtrd_.shape

        logger.info("nof PMTs = {} WF side = {} ".format(NPMT, PMTWL))
        logger.info("nof events in input DST = {} ".format(NEVENTS_DST))

        # create an extensible array to store the CWF waveforms
        # if it exists remove and create again
        if "/RD/pmtcwf" in h5in:
            h5in.remove_node("/RD", "pmtcwf")

        pmtcwf = h5in.create_earray(h5in.root.RD, "pmtcwf",
                                    atom=tables.Int16Atom(),
                                    shape=(0, NPMT, PMTWL),
                                    expectedrows=NEVENTS_DST)
        if "/RD/pmtacum" in h5in:
            h5in.remove_node("/RD", "pmtacum")

        pmtacum = h5in.create_earray(h5in.root.RD, "pmtacum",
                                    atom=tables.Int16Atom(),
                                    shape=(0, NPMT, PMTWL),
                                    expectedrows=NEVENTS_DST)

        # LOOP
        first_evt, last_evt, print_mod = define_event_loop(FIRST_EVT, LAST_EVT,
                                                           NEVENTS,
                                                           NEVENTS_DST,
                                                           RUN_ALL)

        t0 = time()
        for i in range(first_evt, last_evt):
            if not i % print_mod:
                logger.info("-->event number = {}".format(i))

            signal_r, acum = DBLR(pmtrd_, i,
                             n_baseline=N_BASELINE,
                             thr_trigger=THR_TRIGGER,
                             thr_acum=THR_ACUM,
                             acum_discharge_length=ACUM_DISCHARGE_LENGTH,
                             acum_tau=ACUM_TAU,
                             acum_compress=ACUM_COMPRESS)

            # append to pmtcwf
            pmtcwf.append(np.array(signal_r).reshape(1, NPMT, PMTWL))
            # append to pmtacum
            pmtacum.append(np.array(acum).reshape(1, NPMT, PMTWL))

        t1 = time()
        dt = t1 - t0
        pmtcwf.flush()
        pmtacum.flush()

        print("ISIDORA has run over {} events in {} seconds".format(i+1, dt))
    print("Leaving ISIDORA. Safe travels!")

if __name__ == "__main__":
    from cities import isidora
    print(isidora)
    ISIDORA(sys.argv)
