"""
MC_merger takes a bunch of files and produces another one with the same
events.
"""
from __future__ import print_function

import sys
import tables as tb
import numpy as np

import wfmFunctions as wfm
import tblFunctions as tbl
from Nh5 import EventInfo


def print_usage():
    print("\nUsage: python MC_merger.py [outputfile] [inputfiles]\n")
    sys.exit()


def create_new_file(h5out, h5in, nfiles):
    NEVT, NPMT, PMTWL = h5in.root.pmtrd.shape
    NEVT, NSIPM, SIPMWL = h5in.root.sipmrd.shape
    NEVT *= nfiles

    rungroup = h5out.create_group(h5out.root, "Run")
    h5in.root.Run.runInfo.copy(newparent=rungroup)
    evt_out = h5out.create_table(h5out.root.Run, "events", EventInfo,
                                 "Events information",
                                 tbl.filters("NOCOMPR"))

    mcgroup = h5out.create_group(h5out.root, "MC")
    h5in.root.MC.MCTracks.copy(newparent=mcgroup)

    detgroup = h5out.create_group(h5out.root, "Detector")
    h5in.root.Detector.DetectorGeometry.copy(newparent=detgroup)

    snsgroup = h5out.create_group(h5out.root, "Sensors")
    h5in.root.Sensors.DataPMT.copy(newparent=snsgroup)
    h5in.root.Sensors.DataBLR.copy(newparent=snsgroup)
    h5in.root.Sensors.DataSiPM.copy(newparent=snsgroup)

    pmt_out = h5out.create_earray(h5out.root, "pmtrd",
                                  atom=tb.Int16Atom(),
                                  shape=(0, NPMT, PMTWL),
                                  expectedrows=NEVT)

    sipm_out = h5out.create_earray(h5out.root, "sipmrd",
                                   atom=tb.Int16Atom(),
                                   shape=(0, NSIPM, SIPMWL),
                                   expectedrows=NEVT)
    return evt_out, pmt_out, sipm_out


def MC_merger(outputfilename, *inputfilenames, **options):
    COMPRESSION = options.get("compression", "ZLIB4")
    print("Using compression mode", COMPRESSION)

    with tb.open_file(outputfilename, "w",
                      filters=tbl.filters(COMPRESSION)) as h5out:
        create_file = True
        n_events = 0
        for i, filename in enumerate(inputfilenames):
            print("Opening", filename, end="... ")
            sys.stdout.flush()
            try:
                with tb.open_file(filename, "r") as h5in:
                    NEVT = h5in.root.pmtrd.shape[0]
                    n_events += NEVT
                    if create_file:
                        output_data = create_new_file(h5out, h5in,
                                                      len(inputfilenames))
                        evt_out, pmt_out, sipm_out = output_data
                        create_file = False

                    run = h5in.root.Run.events.cols
                    evtrow = evt_out.row
                    pmtrd = h5in.root.pmtrd
                    sipmrd = h5in.root.sipmrd
                    for evt in range(NEVT):
                        evtrow["evt_number"] = run.evt_number[evt]
                        evtrow["timestamp"] = run.timestamp[evt]
                        evtrow.append()

                        pmt_out.append(pmtrd[evt][np.newaxis])
                        sipm_out.append(sipmrd[evt][np.newaxis])

                    evt_out.flush()
                print("OK")
            except:
                print("Error")
        pmt_out.flush()
        sipm_out.flush()
    print("# events = {}".format(n_events))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print_usage()
    MC_merger(sys.argv[1], *sys.argv[2:])
