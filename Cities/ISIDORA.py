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
import tables as tb

from Core.LogConfig import logger
from Core.Configure import configure, define_event_loop, print_configuration
from Core.Nh5 import DECONV_PARAM
import Core.tblFunctions as tbl

import ICython.Sierpe.cBLR as cblr
import Database.loadDB as DB


def DBLR(pmtrwf, n_baseline=500, thr_trigger=5,
         discharge_length=5000,
         acum_tau=2500,
         acum_compress=0.01):
    """
    Peform Base line Restoration
    """
    DataPMT = DB.DataPMT()
    NPMT, PMTWL = pmtrwf.shape
    CWF = np.empty(pmtrwf.shape)
    ACUM = np.empty(pmtrwf.shape)
    BSL = np.empty(pmtrwf.shape[0])
    BSLE = np.empty(pmtrwf.shape[0])
    BSLN = np.empty(pmtrwf.shape[0])

    for pmt in range(NPMT):
        # VH
        signal_r, acum, baseline, baseline_end, noise_rms = cblr.\
          deconvolve_signal_acum(pmtrwf[pmt],
                                 n_baseline=500,
                                 coef_clean=DataPMT.coeff_c[pmt],
                                 coef_blr=DataPMT.coeff_blr[pmt],
                                 thr_trigger=thr_trigger,
                                 acum_discharge_length=discharge_length)

        # JJ
        # signal_r, acum, baseline, baseline_end, noise_rms = cblr.\
        #   deconvolve_signal_acum_v2(pmtrwf[pmt],
        #                             n_baseline=500,
        #                             coef_clean=DataPMT.coeff_c[pmt],
        #                             coef_blr=DataPMT.coeff_blr[pmt],
        #                             thr_trigger=thr_trigger,
        #                             acum_tau=acum_tau,
        #                             acum_compress=acum_compress,
        #                             acum_discharge_length=discharge_length)

        # signal_r, acum = cblr.deconvolve_signal_acum(
        #                  pmtrwf[pmt],
        #                  n_baseline=n_baseline,
        #                  coef_clean=DataPMT.coeff_c[pmt],
        #                  coef_blr=DataPMT.coeff_blr[pmt],
        #                  thr_trigger=thr_trigger,
        #                  thr_acum=thr_acum,
        #                  acum_discharge_length=acum_discharge_length,
        #                  acum_tau=acum_tau,
        #                  acum_compress=acum_compress)
        CWF[pmt] = signal_r
        ACUM[pmt] = acum
        BSL[pmt] = baseline
        BSLE[pmt] = baseline_end
        BSLN[pmt] = noise_rms

    return CWF, ACUM, BSL, BSLE, BSLN


def ISIDORA(argv=sys.argv):
    CFP = configure(argv)

    if CFP["INFO"]:

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

    N_BASELINE = CFP["N_BASELINE"]
    THR_TRIGGER = CFP["THR_TRIGGER"]
    ACUM_DISCHARGE_LENGTH = CFP["ACUM_DISCHARGE_LENGTH"]
    ACUM_TAU = CFP["ACUM_TAU"]
    ACUM_COMPRESS = CFP["ACUM_COMPRESS"]
    COMPRESSION = CFP["COMPRESSION"]

    # open the input file in mode append
    with tb.open_file(CFP["FILE_IN"], "a",
                      filters=tbl.filters(COMPRESSION)) as h5in:
        # access the PMT raw data in file
        pmtrd_ = h5in.root.RD.pmtrwf  # PMT raw data must exist
        NEVENTS_DST, NPMT, PMTWL = pmtrd_.shape

        print_configuration({"# PMT": NPMT, "PMT WL": PMTWL,
                             "# events in DST": NEVENTS_DST})

        # create an extensible array to store the CWF waveforms
        # if it exists remove and create again
        if "/RD/pmtcwf" in h5in:
            h5in.remove_node("/RD", "pmtcwf")

        pmtcwf = h5in.create_earray(h5in.root.RD, "pmtcwf",
                                    atom=tb.Int16Atom(),
                                    shape=(0, NPMT, PMTWL),
                                    expectedrows=NEVENTS_DST,
                                    filters=tbl.filters(COMPRESSION))
        # if "/RD/pmtacum" in h5in:
        #     h5in.remove_node("/RD", "pmtacum")
        #
        # pmtacum = h5in.create_earray(h5in.root.RD, "pmtacum",
        #                              atom=tb.Int16Atom(),
        #                              shape=(0, NPMT, PMTWL),
        #                              expectedrows=NEVENTS_DST)

        if "/Deconvolution" not in h5in:
            h5in.create_group(h5in.root, "Deconvolution")
        if "/Deconvolution/Parameters" in h5in:
            h5in.remove_node("/Deconvolution", "Parameters")
        if "/Deconvolution/BL" in h5in:
            h5in.remove_node("/Deconvolution", "BL")

        deconv_table = h5in.create_table(h5in.root.Deconvolution,
                                         "Parameters",
                                         DECONV_PARAM,
                                         "Deconvolution parameters",
                                         tbl.filters("NOCOMPR"))
        tbl.store_deconv_table(deconv_table, CFP)

        bl_array = h5in.create_earray(h5in.root.Deconvolution, "BL",
                                      atom=tb.Int16Atom(),
                                      shape=(0, NPMT, 3),
                                      expectedrows=NEVENTS_DST,
                                      filters=tbl.filters(COMPRESSION))
        # LOOP
        t0 = time()
        for i in define_event_loop(CFP, NEVENTS_DST):
            data = DBLR(pmtrd_[i],
                        n_baseline=N_BASELINE,
                        thr_trigger=THR_TRIGGER,
                        discharge_length=ACUM_DISCHARGE_LENGTH,
                        acum_tau=ACUM_TAU,
                        acum_compress=ACUM_COMPRESS)
            signal_r, bl_data = data[0], data[2:]

            # append to pmtcwf
            pmtcwf.append(signal_r.reshape(1, NPMT, PMTWL))
            bl_array.append(np.array(bl_data).T[np.newaxis])
            # append to pmtacum
            # pmtacum.append(acum.reshape(1, NPMT, PMTWL))

        t1 = time()
        dt = t1 - t0
        pmtcwf.flush()
        bl_array.flush()
        # pmtacum.flush()

        print("ISIDORA has run over {} events in {} seconds".format(i+1, dt))
    print("Leaving ISIDORA. Safe travels!")


if __name__ == "__main__":
    from cities import isidora
    print(isidora)
    ISIDORA(sys.argv)
