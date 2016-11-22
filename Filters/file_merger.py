"""
Merge a bunch of files into a single one reproducing the structure.

@Author: G. Martinez Lema
"""
from __future__ import print_function

import sys
import importlib
import argparse
import tables as tb
import numpy as np

import Core.tblFunctions as tbl
from Core.Nh5 import EventInfo, SENSOR_WF, PMAP
from Core.Configure import read_config_file, filter_options


def init_filter(filtername, **options):
    """
    Creates an instance of a filter.

    Parameters
    ----------
    filtername : string
        Name of the filter to be used. Must be a class with the __call__
        method or a function taking keyword arguments and returnin another
        function. The filter must be contained in a module with the same name.

    Optional keyword arguments for the filter configuration can be passed
    through the options keyword argument.

    Returns
    -------
    filterinstance : class or function
        Callable object taking an hdf5 file and event number and
        returning a boolean.
    """
    module = importlib.import_module("Filters." + filtername)
    filter_ = getattr(module, filtername)
    globals()[filtername] = module
    print("Initializing filter {}".format(filtername))
    return filter_(**filter_options(options, filtername))


def create_new_file(outputfilename, inputfilename, **options):
    """
    Create a new file copying the structure of the inputfile.

    Parameters
    ----------
    outputfilename : string
        Name of the output file (path included)
    inputfilename : string
        Name of the input file (path included)

    Optional keyword arguments for the file configuration can be passed
    through the options keyword argument. Options:
    compression : string
        Compression library and level to be used. Available options in
        tblFunctions.filters.

    Returns
    -------
    newfile : tables.File instance
        The (already open) new file with the copied structure.
    """
    COMPRESSION = options.get("compression", "ZLIB4")
    print("Using compression mode", COMPRESSION)
    h5out = tb.open_file(outputfilename, "w", tbl.filters(COMPRESSION))

    with tb.open_file(inputfilename, "r") as h5in:
        NEVT = h5in.root.Run.events.cols.evt_number[:].size
        NEVT *= options.get("nfiles", 1)

        rungroup = h5out.create_group(h5out.root, "Run")
        h5in.root.Run.runInfo.copy(newparent=rungroup)
        h5out.create_table(h5out.root.Run, "events", EventInfo,
                           "Events information",
                           tbl.filters("NOCOMPR"))

        if "/MC" in h5in:
            mcgroup = h5out.create_group(h5out.root, "MC")
            h5in.root.MC.MCTracks.copy(newparent=mcgroup)
            if "/MC/FEE" in h5in:
                h5in.root.MC.FEE.copy(newparent=mcgroup)

        if "/Detector" in h5in:
            detgroup = h5out.create_group(h5out.root, "Detector")
            h5in.root.Detector.DetectorGeometry.copy(newparent=detgroup)

        if "/Sensors" in h5in:
            snsgroup = h5out.create_group(h5out.root, "Sensors")
            h5in.root.Sensors.DataPMT.copy(newparent=snsgroup)
            h5in.root.Sensors.DataBLR.copy(newparent=snsgroup)
            h5in.root.Sensors.DataSiPM.copy(newparent=snsgroup)

        if "/pmtrd" in h5in:
            _, NPMT, PMTWL = h5in.root.pmtrd.shape
            h5out.create_earray(h5out.root, "pmtrd",
                                atom=tb.Int16Atom(),
                                shape=(0, NPMT, PMTWL),
                                expectedrows=NEVT)

            _, NSIPM, SIPMWL = h5in.root.sipmrd.shape
            h5out.create_earray(h5out.root, "sipmrd",
                                atom=tb.Int16Atom(),
                                shape=(0, NSIPM, SIPMWL),
                                expectedrows=NEVT)

        if "/RD" in h5in:
            rdgroup = h5out.create_group(h5out.root, "RD")
            _, NPMT, PMTWL = h5in.root.RD.pmtrwf.shape
            h5out.create_earray(rdgroup, "pmtrwf",
                                atom=tb.Int16Atom(),
                                shape=(0, NPMT, PMTWL),
                                expectedrows=NEVT)

            h5out.create_earray(rdgroup, "pmtblr",
                                atom=tb.Int16Atom(),
                                shape=(0, NPMT, PMTWL),
                                expectedrows=NEVT)

            _, NSIPM, SIPMWL = h5in.root.RD.sipmrwf.shape
            h5out.create_earray(rdgroup, "sipmrwf",
                                atom=tb.Int16Atom(),
                                shape=(0, NSIPM, SIPMWL),
                                expectedrows=NEVT)
            if "/RD/pmtcwf" in h5in:
                h5out.create_earray(rdgroup, "pmtcwf",
                                    atom=tb.Int16Atom(),
                                    shape=(0, NPMT, PMTWL),
                                    expectedrows=NEVT)

        if "/TWF" in h5in:
            twfgroup = h5out.create_group(h5out.root, "TWF")
            pmt_twf_table = h5out.create_table(twfgroup, "PMT", SENSOR_WF,
                                               "Store for PMTs TWF",
                                               tbl.filters(COMPRESSION))

            sipm_twf_table = h5out.create_table(twfgroup, "SiPM", SENSOR_WF,
                                                "Store for SiPM TWF",
                                                tbl.filters(COMPRESSION))
            pmt_twf_table.cols.event.create_index()
            sipm_twf_table.cols.event.create_index()

        if "/BLR" in h5in:
            h5out.create_group(h5out.root, "BLR")

            h5out.create_earray(h5out.root.BLR, "mau",
                                atom=tb.Int16Atom(),
                                shape=(0, PMTWL),
                                expectedrows=NEVT)

            h5out.create_earray(h5out.root.BLR, "pulse_on",
                                atom=tb.Int16Atom(),
                                shape=(0, PMTWL),
                                expectedrows=NEVT)

            h5out.create_earray(h5out.root.BLR, "wait_over",
                                atom=tb.Int16Atom(),
                                shape=(0, PMTWL),
                                expectedrows=NEVT)

        if "/ZS" in h5in:
            zsgroup = h5out.create_group(h5out.root, "ZS")
            h5out.create_earray(zsgroup, "PMT",
                                atom=tb.Int16Atom(),
                                shape=(0, NPMT, PMTWL),
                                expectedrows=NEVT)

            h5out.create_earray(zsgroup, "BLR",
                                atom=tb.Int16Atom(),
                                shape=(0, NPMT, PMTWL),
                                expectedrows=NEVT)

            h5out.create_earray(zsgroup, "SiPM",
                                atom=tb.Int16Atom(),
                                shape=(0, NSIPM, SIPMWL),
                                expectedrows=NEVT)

        if "/PMAPS" in h5in:
            pmapgroup = h5out.create_group(h5out.root, "PMAPS")
            h5out.create_table(pmapgroup, "PMaps", PMAP,
                               "Store for PMaps", tbl.filters(COMPRESSION))

            h5out.create_table(pmapgroup, "PMapsBLR", PMAP,
                               "Store for PMaps made with BLR",
                               tbl.filters(COMPRESSION))

    return h5out


def file_merger(outputfilename, discardedfilename, *inputfilenames, **options):
    options["nfiles"] = len(inputfilenames)

    filters = options.get("FILTERS", [])
    if not isinstance(filters, list):
        filters = [filters]
    filters = map(lambda f: init_filter(f, **options), filters)

    h5out = create_new_file(outputfilename, inputfilenames[0], **options)
    evtrow_out = h5out.root.Run.events.row

    if "/pmtrd" in h5out:
        pmtrd_out = h5out.root.pmtrd
        sipmrd_out = h5out.root.sipmrd

    if "/RD" in h5out:
        pmtrwf_out = h5out.root.RD.pmtrwf
        pmtblr_out = h5out.root.RD.pmtblr
        sipmrwf_out = h5out.root.RD.sipmrwf
        if "/RD/pmtcwf" in h5out:
            pmtcwf_out = h5out.root.RD.pmtcwf

    if "/TWF" in h5out:
        pmttwf_out = h5out.root.TWF.PMT
        sipmtwf_out = h5out.root.TWF.SiPM

    if "/BLR" in h5out:
        mau_out = h5out.root.BLR.mau
        pulse_on_out = h5out.root.BLR.pulse_on
        wait_over_out = h5out.root.BLR.wait_over

    if "/ZS" in h5out:
        pmtzs_out = h5out.root.ZS.PMT
        blrzs_out = h5out.root.ZS.BLR
        sipmzs_out = h5out.root.ZS.SiPM

    if "/PMAPS" in h5out:
        pmaps_out = h5out.root.PMAPS.PMaps
        pmaps_blr_out = h5out.root.PMAPS.PMapsBLR

    dump_unselected = discardedfilename is not None
    if dump_unselected:
        h5dis = create_new_file(discardedfilename, inputfilenames[0],
                                **options)
        evtrow_dis = h5dis.root.Run.events.row

        if "/pmtrd" in h5dis:
            pmtrd_dis = h5dis.root.pmtrd
            sipmrd_dis = h5dis.root.sipmrd

        if "/RD" in h5dis:
            pmtrwf_dis = h5dis.root.RD.pmtrwf
            pmtblr_dis = h5dis.root.RD.pmtblr
            sipmrwf_dis = h5dis.root.RD.sipmrwf
            if "/RD/pmtcwf" in h5dis:
                pmtcwf_dis = h5dis.root.RD.pmtcwf

        if "/TWF" in h5dis:
            pmttwf_dis = h5dis.root.TWF.PMT
            sipmtwf_dis = h5dis.root.TWF.SiPM

        if "/BLR" in h5dis:
            mau_dis = h5dis.root.BLR.mau
            pulse_on_dis = h5dis.root.BLR.pulse_on
            wait_over_dis = h5dis.root.BLR.wait_over

        if "/ZS" in h5dis:
            pmtzs_dis = h5dis.root.ZS.PMT
            blrzs_dis = h5dis.root.ZS.BLR
            sipmzs_dis = h5dis.root.ZS.SiPM

        if "/PMAPS" in h5dis:
            pmaps_dis = h5dis.root.PMAPS.PMaps
            pmaps_blr_dis = h5dis.root.PMAPS.PMapsBLR

    n_events_in = 0
    n_events_out = 0
    n_events_dis = 0
    for i, filename in enumerate(inputfilenames):
        print("Opening", filename, end="... ")
        sys.stdout.flush()
        try:
            with tb.open_file(filename, "r") as h5in:
                run = h5in.root.Run.events.cols
                NEVT = run.evt_number[:].size
                n_events_in += NEVT

                if "/pmtrd" in h5in:
                    pmtrd_in = h5in.root.pmtrd
                    sipmrd_in = h5in.root.sipmrd

                if "/RD" in h5in:
                    pmtrwf_in = h5in.root.RD.pmtrwf
                    pmtblr_in = h5in.root.RD.pmtblr
                    sipmrwf_in = h5in.root.RD.sipmrwf
                    if "/RD/pmtcwf" in h5in:
                        pmtcwf_in = h5in.root.RD.pmtcwf

                if "/TWF" in h5in:
                    pmttwf_in = h5in.root.TWF.PMT
                    sipmtwf_in = h5in.root.TWF.SiPM

                if "/BLR" in h5in:
                    mau_in = h5in.root.BLR.mau
                    pulse_on_in = h5in.root.BLR.pulse_on
                    wait_over_in = h5in.root.BLR.wait_over

                if "/ZS" in h5in:
                    pmtzs_in = h5in.root.ZS.PMT
                    blrzs_in = h5in.root.ZS.BLR
                    sipmzs_in = h5in.root.ZS.SiPM

                if "/PMAPS" in h5in:
                    pmaps_in = h5in.root.PMAPS.PMaps
                    pmaps_blr_in = h5in.root.PMAPS.PMapsBLR

                for evt in range(NEVT):
                    if all([filter_(h5in, evt) for filter_ in filters]):
                        evtrow_out["evt_number"] = run.evt_number[evt]
                        evtrow_out["timestamp"] = run.timestamp[evt]
                        evtrow_out.append()

                        if "/pmtrd" in h5out:
                            pmtrd_out.append(pmtrd_in[evt][np.newaxis])
                            sipmrd_out.append(sipmrd_in[evt][np.newaxis])

                        if "/RD" in h5out:
                            pmtrwf_out.append(pmtrwf_in[evt][np.newaxis])
                            pmtblr_out.append(pmtblr_in[evt][np.newaxis])
                            sipmrwf_out.append(sipmrwf_in[evt][np.newaxis])
                            if "/RD/pmtcwf" in h5out:
                                pmtcwf_out.append(pmtcwf_in[evt][np.newaxis])

                        if "/TWF" in h5out:
                            wf = tbl.read_wf_table(pmttwf_in, evt)
                            tbl.store_wf(pmttwf_out, evt, wf)
                            wf = tbl.read_wf_table(sipmtwf_in, evt)
                            tbl.store_wf(sipmtwf_out, evt, wf)

                        if "/BLR" in h5out:
                            mau_out.append(mau_in[evt][np.newaxis])
                            pulse_on_out.append(pulse_on_in[evt][np.newaxis])
                            wait_over_out.append(wait_over_in[evt][np.newaxis])

                        if "/ZS" in h5out:
                            pmtzs_out.append(pmtzs_in[evt][np.newaxis])
                            blrzs_out.append(blrzs_in[evt][np.newaxis])
                            sipmzs_out.append(sipmzs_in[evt][np.newaxis])

                        if "/PMAPS" in h5in:
                            pmap = tbl.read_pmap(pmaps_in, evt)
                            tbl.store_pmap(pmap, pmaps_out, evt)
                            pmap = tbl.read_pmap(pmaps_blr_in, evt)
                            tbl.store_pmap(pmap, pmaps_blr_out, evt)

                        n_events_out += 1
                    else:
                        evtrow_dis["evt_number"] = run.evt_number[evt]
                        evtrow_dis["timestamp"] = run.timestamp[evt]
                        evtrow_dis.append()

                        if "/pmtrd" in h5dis:
                            pmtrd_dis.append(pmtrd_in[evt][np.newaxis])
                            sipmrd_dis.append(sipmrd_in[evt][np.newaxis])

                        if "/RD" in h5dis:
                            pmtrwf_dis.append(pmtrwf_in[evt][np.newaxis])
                            pmtblr_dis.append(pmtblr_in[evt][np.newaxis])
                            sipmrwf_dis.append(sipmrwf_in[evt][np.newaxis])
                            if "/RD/pmtcwf" in h5dis:
                                pmtcwf_dis.append(pmtcwf_in[evt][np.newaxis])

                        if "/TWF" in h5dis:
                            wf = tbl.read_wf_table(pmttwf_in, evt)
                            tbl.store_wf(pmttwf_dis, evt, wf)
                            wf = tbl.read_wf_table(sipmtwf_in, evt)
                            tbl.store_wf(sipmtwf_dis, evt, wf)

                        if "/BLR" in h5dis:
                            mau_dis.append(mau_in[evt][np.newaxis])
                            pulse_on_dis.append(pulse_on_in[evt][np.newaxis])
                            wait_over_dis.append(wait_over_in[evt][np.newaxis])

                        if "/ZS" in h5dis:
                            pmtzs_dis.append(pmtzs_in[evt][np.newaxis])
                            blrzs_dis.append(blrzs_in[evt][np.newaxis])
                            sipmzs_dis.append(sipmzs_in[evt][np.newaxis])

                        if "/PMAPS" in h5in:
                            pmap = tbl.read_pmap(pmaps_in, evt)
                            tbl.store_pmap(pmap, pmaps_dis, evt)
                            pmap = tbl.read_pmap(pmaps_blr_in, evt)
                            tbl.store_pmap(pmap, pmaps_blr_dis, evt)

                        n_events_dis += 1
                h5out.root.Run.events.flush()
                if dump_unselected:
                    h5dis.root.Run.events.flush()
            print("OK")
        except:
            print("Error")
    ratio_out = n_events_out * 100. / n_events_in
    ratio_dis = n_events_dis * 100. / n_events_in
    print("# events in = {}".format(n_events_in))
    print("# events accepted = {} ({:.2f}%)".format(n_events_out, ratio_out))
    print("# events discarded = {} ({:.2f}%)".format(n_events_dis, ratio_dis))
    h5out.close()
    if dump_unselected:
        h5dis.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", metavar="ifile", type=str, nargs="+",
                        help="input files to be merged", required=True)
    parser.add_argument("-o", metavar="ofile", type=str,
                        help="output file", required=True)
    parser.add_argument("-d", metavar="dfile", type=str,
                        help="output file with discarded events")
    parser.add_argument("-c", metavar="cfile", type=str,
                        help="configuration file")

    args = parser.parse_args()
    options = read_config_file(args.c) if args.c else {}
    file_merger(args.o, args.d, *args.i, **options)
