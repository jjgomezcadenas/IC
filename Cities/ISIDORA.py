"""
ISIDORA
JJGC Agusut 2016

What ISIDORA does:
1) Reads a RWF file written by DIOMIRA
2) Performs DBLR
3) Write the corrected waveforms (CWF) to the file as new Evectors.

Change Log:
18/10 JJG added cython module
"""

from __future__ import print_function

import sys
from time import time
import numpy as np
import tables

import system_of_units as units
from LogConfig import logger
from Configure import configure, define_event_loop

import FEParam as FP
import coreFunctions as cf

CYTHON = True
from BLR import accumulator_coefficients, BLR


def DBLR(pmtrd, event_number, coeff_acc, mau_len=250,
         thr1 = FP.NOISE_ADC, thr2=0, thr3 = FP.NOISE_ADC):
    """
    Peform Base line Restoration
    coeff_acc is an array with the coefficients of the accumulator
    Threshold 1 is used to decide when raw signal raises up above trigger line
    Threshold 2 is used to decide when reconstructed signal is above trigger line
    Threshold 3 is used to compare Raw and Rec signal
    """

    len_WF = pmtrd.shape[2]
    NPMT = pmtrd.shape[1]
    BLRS =[]

    for j in range(NPMT):
        sgn_raw = FP.ceiling - pmtrd[event_number, j]
        sgn_rec, MAU, pulse_on, wait_over = cBLR(sgn_raw,
                                                 coeff_acc[j],
                                                 mau_len, thr1, thr2, thr3)

        BLRS.append((sgn_rec, MAU, pulse_on, wait_over))

    return zip(*BLRS)


def ISIDORA(argv):
    DEBUG_LEVEL, INFO, CFP = configure(argv[0],argv[1:])

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

    PATH_IN = CFP['PATH_IN']
    FILE_IN = CFP['FILE_IN']
    FIRST_EVT = CFP['FIRST_EVT']
    LAST_EVT = CFP['LAST_EVT']
    RUN_ALL = CFP['RUN_ALL']
    COEF = CFP['COEF']
    CA = cf.farray_from_string(CFP['CA'])*units.nF
    AC = cf.farray_from_string(CFP['AC'])
    MAU_LEN = CFP['MAU_LEN']
    NSIGMA1 = CFP['NSIGMA1']
    NSIGMA2 = CFP['NSIGMA2']
    NSIGMA3 = CFP['NSIGMA3']

    NEVENTS = LAST_EVT - FIRST_EVT

    logger.info('Debug level = {}'.format(DEBUG_LEVEL))

    logger.info("input path ={}; file_in ={} ".format(
        PATH_IN,FILE_IN))

    logger.info("""first event = {} last event = {}
                    nof events requested = {} """.format(FIRST_EVT,
                                                         LAST_EVT,
                                                         NEVENTS))

    logger.info("""MAU length = {}
                   n_sigma1 = {} n_sigma2 = {} """.format(MAU_LEN,
                                                          NSIGMA1,
                                                          NSIGMA2))
    logger.info("CA  = {} nF ".format(CA/nF))
    logger.info("Accumulator Coefficients = {}  ".format(AC))

    # open the input file in mode append
    with tables.open_file("{}/{}".format(PATH_IN, FILE_IN), "a") as h5in:
        # access the PMT raw data in file
        pmtrd_ = h5in.root.RD.pmtrwf     # PMT raw data must exist

        #pmtrd_.shape = (nof_events, nof_sensors, wf_length)
        NEVENTS_DST, NPMT, PMTWL = pmtrd_.shape

        logger.info("nof PMTs = {} WF side = {} ".format(NPMT, PMTWL))
        logger.info("nof events in input DST = {} ".format(NEVENTS_DST))

        # create an extensible array to store the CWF waveforms
        # if it exists remove and create again
        if '/RD/pmtcwf' in h5in:
            pmtcwf = h5in.root.RD.pmtcwf
            h5in.remove_node("/RD","pmtcwf")

        pmtcwf = h5in.create_earray(h5in.root.RD, "pmtcwf",
                                atom=tables.Int16Atom(),
                                shape=(0, NPMT, PMTWL),
                                expectedrows=NEVENTS_DST)

        # create a group to store BLR configuration
        # mau, acummulator, pulse_on and wait_over stored for pmt 0
        # baseline stored for all PMTs.
        if not '/BLR' in h5in:
            rgroup = h5in.create_group(h5in.root, "BLR")

        if '/BLR/mau' in h5in:
            mau = h5in.root.BLR.mau
            h5in.remove_node("/BLR","mau")
        mau = h5in.create_earray(h5in.root.BLR, "mau",
                                atom=tables.Int16Atom(),
                                shape=(0, PMTWL),
                                expectedrows=NEVENTS_DST)

        if '/BLR/pulse_on' in h5in:
            pulse_on = h5in.root.BLR.pulse_on
            h5in.remove_node("/BLR","pulse_on")
        pulse_on = h5in.create_earray(h5in.root.BLR, "pulse_on",
                                atom=tables.Int16Atom(),
                                shape=(0, PMTWL),
                                expectedrows=NEVENTS_DST)

        if '/BLR/wait_over' in h5in:
            wait_over = h5in.root.BLR.wait_over
            h5in.remove_node("/BLR","wait_over")
        wait_over = h5in.create_earray(h5in.root.BLR, "wait_over",
                                atom=tables.Int16Atom(),
                                shape=(0, PMTWL),
                                expectedrows=NEVENTS_DST)

        # compute the accumulator coefficients from the nominal values
        # of capacitances (if COEF =0) or use previously computed coeff
        coeff_acc = accumulator_coefficients(CA,
                                             NPMT, PMTWL) if COEF == 0 else AC

        # LOOP
        first_evt, last_evt = define_event_loop(FIRST_EVT, LAST_EVT,
                                                NEVENTS, NEVENTS_DST, RUN_ALL)

        t0 = time()
        for i in range(first_evt,last_evt):

            logger.info("-->event number ={}".format(i))
            signal_r, xmau, pulse_, wait_ = DBLR(pmtrd_, i, coeff_acc,
                                                 mau_len=MAU_LEN,
                                                 thr1 = NSIGMA1*FP.NOISE_ADC,
                                                 thr2=NSIGMA2*FP.NOISE_ADC,
                                                 thr3 = NSIGMA3*FP.NOISE_ADC)
            # append to pmtcwd
            pmtcwf.append(np.array(signal_r).reshape(1, NPMT, PMTWL))

            # append BLR variables
            mau.append(xmau[0].reshape(1, PMTWL))
            pulse_on.append (pulse_[0].reshape(1, PMTWL))
            wait_over.append(wait_[0].reshape(1, PMTWL))

        t1 = time()
        pmtcwf.flush()
        mau.flush()
        pulse_on.flush()
        wait_over.flush()

        print("ISIDORA has run over {} events in {} seconds".format(i+1,
                                                                    t1-t0))
    print("Leaving ISIDORA. Safe travels!")

if __name__ == '__main__':
    #ISIDORA(sys.argv)
