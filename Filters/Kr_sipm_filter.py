"""
Kr_sipm_filter takes a bunch of files and applies a filter on them,
rejecting those that are not suitable and merging the remaining into
a single new file.
"""
from __future__ import print_function

import sys
import tables as tb
import numpy as np

import tblFunctions as tbl


def print_usage():
    print("\nUsage: python Kr_sipm_filter.py [outputfile] [inputfiles]\n")
    sys.exit()


def create_new_file(h5out, h5in, nfiles):
    NEVT, NPMT, PMTWL = h5in.root.RD.pmtrwf.shape
    NEVT, NSIPM, SIPMWL = h5in.root.RD.sipmrwf.shape
    NEVT *= nfiles

    rungroup = h5out.create_group(h5out.root, "Run")
    h5in.root.Run.runInfo.copy(newparent=rungroup)
    evt_out = h5out.create_earray(h5out.root.Run, "event_number",
                                  atom=tb.Int16Atom(),
                                  shape=(0,),
                                  expectedrows=NEVT)

    detgroup = h5out.create_group(h5out.root, "Detector")
    h5in.root.Detector.DetectorGeometry.copy(newparent=detgroup)

    snsgroup = h5out.create_group(h5out.root, "Sensors")
    h5in.root.Sensors.DataPMT.copy(newparent=snsgroup)
    h5in.root.Sensors.DataBLR.copy(newparent=snsgroup)
    h5in.root.Sensors.DataSiPM.copy(newparent=snsgroup)

    h5out.create_group(h5out.root, "RD")
    pmt_out = h5out.create_earray(h5out.root.RD, "pmtrwf",
                                  atom=tb.Int16Atom(),
                                  shape=(0, NPMT, PMTWL),
                                  expectedrows=NEVT)

    blr_out = h5out.create_earray(h5out.root.RD, "pmtblr",
                                  atom=tb.Int16Atom(),
                                  shape=(0, NPMT, PMTWL),
                                  expectedrows=NEVT)

    sipm_out = h5out.create_earray(h5out.root.RD, "sipmrwf",
                                   atom=tb.Int16Atom(),
                                   shape=(0, NSIPM, SIPMWL),
                                   expectedrows=NEVT)
    return evt_out, pmt_out, blr_out, sipm_out


def filter_events(h5in):
    if "/RD/sipmrwf" not in h5in:
        return np.array([])
    else:
        return np.arange(h5in.root.RD.sipmrwf.shape[0])


def Kr_sipm_filter(outputfilename, *inputfilenames, **options):
    COMPRESSION = options.get("compression", "ZLIB4")
    print("Using compression mode", COMPRESSION)

    with tb.open_file(outputfilename, "w",
                      filters=tbl.filters(COMPRESSION)) as h5out:
        create_file = True
        n_events_in = 0
        n_events_out = 0
        for i, filename in enumerate(inputfilenames):
            print("Opening", filename, end="... ")
            with tb.open_file(filename, "r") as h5in:
                filtered_events = filter_events(h5in)
                n_events_in += h5in.root.RD.pmtrwf.shape[0]
                n_events_out += len(filtered_events)
                if create_file and filtered_events.size:
                    evt_out, pmt_out, blr_out, sipm_out = create_new_file(
                                                          h5out, h5in,
                                                          len(inputfilenames))
                    create_file = False
                for evt in filtered_events:
                    evt_out.append(h5in.root.Run.event_number[evt][np.newaxis])
                    pmt_out.append(h5in.root.RD.pmtrwf[evt][np.newaxis])
                    blr_out.append(h5in.root.RD.pmtblr[evt][np.newaxis])
                    sipm_out.append(h5in.root.RD.sipmrwf[evt][np.newaxis])
            print("Done")
        pmt_out.flush()
        blr_out.flush()
        sipm_out.flush()
    print("# events in = {}, # events out = {}, ratio = {}"
          "".format(n_events_in, n_events_out,
                    n_events_out * 1.0 / n_events_in))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print_usage()
    Kr_sipm_filter(sys.argv[1], *sys.argv[2:])
