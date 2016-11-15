"""
DOROTHEA
JJGC August-October 2016
GML October 2016

What DOROTHEA does:
1) Reads a hdf5 file containing ZS waveform for all PMTs and SiPMs in adc.
2) Converts waveforms to pes and creates a summed PMT.
2) Finds the peaks in the summed waveform and links them with those found in
the SiPMs in the same time window.
3) Rebins the PMT waveforms to 1 mus.
4) Writes pmaps into a new file.
"""

from __future__ import print_function
import sys
import math
import numpy as np
import tables
from time import time

import system_of_units as units
from LogConfig import logger
from Configure import configure, define_event_loop
from HLObjects import Signal, Peak, PMap
from Nh5 import PMAP

import tblFunctions as tbl


"""

DOROTHEA
ChangeLog:

24.10 First version

10.11 Conversion to pes and creation of the summed PMT.
"""


def classify_signal(slices, foundS2):
    if len(slices)>1:
        sig = Signal.S2
    elif not foundS2:
        sig = Signal.S1
    else:
        sig = Signal.UNKNOWN
    return sig


def build_pmap(pmtwf, sipmwfs, stride=40):
    """
    Finds any peak in the waveform and rebins it.
    """
    to_mus = 25*stride*units.ns/units.mus

    nbins = int(math.ceil(len(pmtwf)*1.0/stride))
    pmap = PMap()
    ene_pmt = []
    ene_sipms = []
    time_over_thrs = []
    tmin = float("inf")
    foundS2 = False

    for i in range(nbins):
        low = i * stride
        upp = low + stride
        slice_ = pmtwf[low:upp]
        e = slice_.sum()
        q = sipmwfs[:, i].flatten()
        t = i

        # Non-empty slice, append it and carry on
        if e > 0.:
            if t<tmin:
                tmin = t
            ene_pmt.append(e)
            # q = np.concatenate((q,np.zeros(3)))
            ene_sipms.append(q)
            time_over_thrs.append(np.nonzero(slice_)[0].size)

        # Empty slice. Everything accumulated so far is a peak.
        # It will be S1-like if it is a short peak
        elif len(ene_pmt) > 0:
            sigtype = classify_signal(ene_pmt, foundS2)
            if sigtype==Signal.S2:
                foundS2 = True
            tmax = t
            peak = Peak(np.arange(tmin, tmax)*to_mus,
                        ene_pmt, ene_sipms,
                        time_over_thrs, sigtype)
            pmap.peaks.append(peak)

            rebin_wf = []
            tmin = float("inf")
            ene_pmt = []
            ene_sipms = []
            time_over_thrs = []

    return pmap


def DOROTHEA(argv):
    """
    DOROTHEA driver
    """
    DEBUG_LEVEL, INFO, CFP = configure(argv[0], argv[1:])

    if INFO:
        print(__doc__)

    PATH_IN = CFP["PATH_IN"]
    PATH_OUT = CFP["PATH_OUT"]
    FILE_IN = CFP["FILE_IN"]
    FILE_OUT = CFP["FILE_OUT"]
    PATH_DB = CFP["PATH_DB"]
    FIRST_EVT = CFP["FIRST_EVT"]
    LAST_EVT = CFP["LAST_EVT"]
    RUN_ALL = CFP["RUN_ALL"]
    COMPRESSION = CFP["COMPRESSION"]
    NEVENTS = LAST_EVT - FIRST_EVT

    logger.info("Debug level = {}".format(DEBUG_LEVEL))
    logger.info("Input path = {}; output path = {}".format(PATH_IN, PATH_OUT))
    logger.info("File_in = {} file_out = {}".format(FILE_IN, FILE_OUT))
    logger.info("Path to database = {}".format(PATH_DB))
    logger.info("First event = {} last event = {} "
                "# events requested = {}".format(FIRST_EVT, LAST_EVT, NEVENTS))
    logger.info("Compression library/level = {}".format(COMPRESSION))

    # open the input file
    with tables.open_file("{}/{}".format(PATH_IN, FILE_IN), "r") as h5in:

        # access the PMT ZS data in file
        pmtzs_ = h5in.root.ZS.PMT
        blrzs_ = h5in.root.ZS.BLR
        sipmzs_ = h5in.root.ZS.SiPM

        NEVT, NPMT, PMTWL = pmtzs_.shape
        NEVT, NSIPM, SIPMWL = sipmzs_.shape

        logger.info("# events in DST: {}".format(NEVT))
        logger.info("# PMTs = {}, # SiPMs = {} ".format(NPMT, NSIPM))
        logger.info("PMT WFL = {}, SiPM WFL = {}".format(PMTWL, SIPMWL))

        # access the geometry and the sensors metadata info
        geom_t = h5in.root.Detector.DetectorGeometry
        pmt_t = h5in.root.Sensors.DataPMT
        blr_t = h5in.root.Sensors.DataBLR
        sipm_t = h5in.root.Sensors.DataSiPM

        pmtdf = tbl.read_sensors_table(pmt_t)
        blrdf = tbl.read_sensors_table(blr_t)
        sipmdf = tbl.read_sensors_table(sipm_t)

        pmt_to_pes = abs(1.0 / pmtdf.adc_to_pes.reshape(NPMT, 1))
        blr_to_pes = abs(1.0 / blrdf.adc_to_pes.reshape(NPMT, 1))
        sipm_to_pes = abs(1.0 / sipmdf.adc_to_pes.reshape(NSIPM, 1))

        # open the output file
        with tables.open_file("{}/{}".format(PATH_OUT, FILE_OUT), "w",
                              filters=tbl.filters(COMPRESSION)) as h5out:

            # create groups and copy MC data to the new file
            if "/MC" in h5in:
                mcgroup = h5out.create_group(h5out.root, "MC")
                h5in.root.MC.MCTracks.copy(newparent=mcgroup)
                h5in.root.MC.FEE.copy(newparent=mcgroup)

            detgroup = h5out.create_group(h5out.root, "Detector")
            geom_t.copy(newparent=detgroup)

            sgroup = h5out.create_group(h5out.root, "Sensors")
            pmt_t.copy(newparent=sgroup)
            blr_t.copy(newparent=sgroup)
            sipm_t.copy(newparent=sgroup)

            pmapsgroup = h5out.create_group(h5out.root, "PMAPS")

            # create a table to store pmaps (rebined, linked, zs wfs)
            pmaps_ = h5out.create_table(pmapsgroup, "PMaps", PMAP,
                                        "Store for PMaps",
                                        tbl.filters(COMPRESSION))

            pmaps_blr_ = h5out.create_table(pmapsgroup, "PMapsBLR", PMAP,
                                            "Store for PMaps made with BLR",
                                            tbl.filters(COMPRESSION))

            # add index in event column
            pmaps_.cols.event.create_index()
            pmaps_blr_.cols.event.create_index()

            # LOOP
            first_evt, last_evt, print_mod = define_event_loop(FIRST_EVT,
                                                               LAST_EVT,
                                                               NEVENTS, NEVT,
                                                               RUN_ALL)
            t0 = time()
            for i in range(first_evt, last_evt):
                if not i%print_mod:
                    logger.info("-->event number = {}".format(i))

                pmtwf = np.sum(pmtzs_[i] * pmt_to_pes, axis=0)
                blrwf = np.sum(blrzs_[i] * blr_to_pes, axis=0)
                sipmwfs = sipmzs_[i] * sipm_to_pes

                pmap = build_pmap(pmtwf, sipmwfs)
                tbl.store_pmap(pmap, pmaps_, i)

                pmap_blr = build_pmap(blrwf, sipmwfs)
                tbl.store_pmap(pmap_blr, pmaps_blr_, i)

            t1 = time()
            dt = t1 - t0
            print("DOROTHEA has run over {} events in {} seconds".format(i+1,
                                                                         dt))
    print("Leaving DOROTHEA. Safe travels!")

if __name__ == "__main__":
    from cities import dorothea
    print(dorothea)
    DOROTHEA(sys.argv)
