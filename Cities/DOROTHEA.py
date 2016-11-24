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
import tables as tb
from time import time

import Core.system_of_units as units
from Core.LogConfig import logger
from Core.Configure import configure, define_event_loop, print_configuration
from Core.Bridges import Signal, Peak, PMap
from Core.Nh5 import PMAP

import Core.tblFunctions as tbl
import Database.loadDB as DB


"""

DOROTHEA
ChangeLog:

24.10 First version

10.11 Conversion to pes and creation of the summed PMT.

16.11 Using new database utility.
"""


def classify_peaks(pmap, **options):
    """
    Classify peaks according to given criteria.
    """
    foundS1 = False
    foundS2 = False
    s1_min_int = options.get("MIN_S1_INTEGRAL", 0.)
    s1_max_tot = options.get("MAX_S1_ToT", 40)
    s2_min_wid = options.get("MIN_S2_WIDTH", 0)
    s2_min_hei = options.get("MIN_S2_HEIGHT", 0.)
    s2_min_int = options.get("MIN_S2_INTEGRAL", 0.)
    for peak in pmap:
        peak.signal = Signal.UNKNOWN
        if (len(peak) >= s2_min_wid and peak.cathode_integral > s2_min_int and
           peak.peakmax[1] > s2_min_hei):
            peak.signal = Signal.S2
            foundS2 = True
        elif len(peak) == 1 and peak.tothrs.sum() < s1_max_tot and not foundS2:
            if not foundS1 or peak.cathode_integral > s1_min_int:
                peak.signal = Signal.S1
                foundS1 = True

    if options.get("FILTER_OUTPUT", False):
        filter_peaks(pmap)


def filter_peaks(pmap):
    """
    Remove unknowns.
    """
    for i in reversed(range(len(pmap.peaks))):
        if pmap.peaks[i].signal == Signal.UNKNOWN:
            pmap.peaks.pop(i)


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

    for i in range(nbins):
        low = i * stride
        upp = low + stride
        slice_ = pmtwf[low:upp]
        e = slice_.sum()
        q = sipmwfs[:, i].flatten()
        t = i

        # Non-empty slice, append it and carry on
        if e > 0.:
            if t < tmin:
                tmin = t
            ene_pmt.append(e)
            # q = np.concatenate((q,np.zeros(3)))
            ene_sipms.append(q)
            time_over_thrs.append(np.nonzero(slice_)[0].size)

        # Empty slice. Everything accumulated so far is a peak.
        # It will be S1-like if it is a short peak
        elif len(ene_pmt) > 0:
            tmax = t
            peak = Peak(np.arange(tmin, tmax)*to_mus,
                        ene_pmt, ene_sipms, time_over_thrs)
            pmap.peaks.append(peak)

            tmin = float("inf")
            ene_pmt = []
            ene_sipms = []
            time_over_thrs = []

    return pmap


def DOROTHEA(argv=sys.argv):
    """
    DOROTHEA driver
    """
    CFP = configure(argv)

    if CFP["INFO"]:
        print(__doc__)

    COMPRESSION = CFP["COMPRESSION"]

    # open the input file
    with tb.open_file(CFP["FILE_IN"], "r") as h5in:
        # access the PMT ZS data in file
        pmtzs_ = h5in.root.ZS.PMT
        blrzs_ = h5in.root.ZS.BLR
        sipmzs_ = h5in.root.ZS.SiPM

        NEVT, NPMT, PMTWL = pmtzs_.shape
        NEVT, NSIPM, SIPMWL = sipmzs_.shape

        print_configuration({"# PMT": NPMT, "PMT WL": PMTWL,
                             "# SiPM": NSIPM, "SIPM WL": SIPMWL,
                             "# events in DST": NEVT})

        pmtdf = DB.DataPMT()
        sipmdf = DB.DataSiPM()

        pmt_to_pes = abs(1.0 / pmtdf.adc_to_pes.reshape(NPMT, 1))
        sipm_to_pes = abs(1.0 / sipmdf.adc_to_pes.reshape(NSIPM, 1))

        # open the output file
        with tb.open_file(CFP["FILE_OUT"], "w",
                          filters=tbl.filters(COMPRESSION)) as h5out:

            # create groups and copy MC data to the new file
            if "/MC" in h5in:
                mcgroup = h5out.create_group(h5out.root, "MC")
                twfgroup = h5out.create_group(h5out.root, "TWF")

                h5in.root.MC.MCTracks.copy(newparent=mcgroup)
                h5in.root.MC.FEE.copy(newparent=mcgroup)
                h5in.root.TWF.PMT.copy(newparent=twfgroup)
                h5in.root.TWF.SiPM.copy(newparent=twfgroup)

            rungroup = h5out.create_group(h5out.root, "Run")
            h5in.root.Run.runInfo.copy(newparent=rungroup)
            h5in.root.Run.events.copy(newparent=rungroup)

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
            t0 = time()
            for i in define_event_loop(CFP, NEVT):
                pmtwf = np.sum(pmtzs_[i] * pmt_to_pes, axis=0)
                blrwf = np.sum(blrzs_[i] * pmt_to_pes, axis=0)
                sipmwfs = sipmzs_[i] * sipm_to_pes

                pmap = build_pmap(pmtwf, sipmwfs)
                classify_peaks(pmap, **CFP)
                tbl.store_pmap(pmap, pmaps_, i)

                pmap_blr = build_pmap(blrwf, sipmwfs)
                classify_peaks(pmap_blr, **CFP)
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
