
# TAMARA
#
# TAMARA is a main functions that reads RW (Corrected RW for PMTs) and SiPM WF
# compute the 'PMaps' (a la Hernnado), selects PMAPs and output a DataFrame
# with the S1 & S2 information
#
# Used for the analysis presented at Canfranc 1/12/16
#
# author: JA Hernando
# date: 21/11/2016
#

from __future__ import print_function

import functools as ft
import pandas as pd
import tables as tb
import numpy as np

import Database.loadDB as db
from Core.Bridges import Peak, PMap

import tamFunctions as tamf

# --------------------------------------
# Parameters
# --------------------------------------

# Input Oputput
idatapath = '/Users/hernando/Investigacion/NEXT/Data/NEWStar/Na/r2948/'
ifilename = 'run_2948_idfile_RWF.h5'
idfiles = [str(i) for i in range(0, 10)]
odatapath = 'data/'
ofilename = 'r2948_tam_s1trg'
nevents = -1
ifiles = [idatapath+ifilename.replace('idfile', id) for id in idfiles]
ofile = odatapath+ofilename

# Pmaps configuration
q0cat = 4.  # pes
q0ano = 1.2  # pes


# PMaps selection
tscale = 0.025
s1time = (95., 105.)
s2ene = (0., 1.e6)
select_pmap = tamf.pmap_is_gold
select_s1 = ft.partial(tamf.s1_has_time, range=s1time, scale=tscale)
select_s2 = ft.partial(tamf.s2_has_energy, range=s2ene, scale=1.)

#  select_s1 = tambda s1: tamf.s1_has_time(s1, range=s1time, scale=tscake)
#  select_s2 = lambda s2: tamf.s2_has_energy(s2, range=s2ene, scale = 1.)

# Set options into the opts dictionary
opts = {}
opts['q0cat'] = q0cat
opts['q0ano'] = q0ano
opts['select_pmap'] = select_pmap
opts['select_s1'] = select_s1
opts['select_s2'] = select_s2


# ----------------------------
#  Driver
# ----------------------------


def run(ifiles, ofile, nevents, **opts):

    pmaps, pan = [], []
    for i, ifile in enumerate(ifiles):
        print('reading ', ifile)
        ipmaps, ipan = run_file(ifile, nevents, **opts)
        pmaps += ipmaps
        if (len(pan) == 0):
            pan = ipan
        else:
            pan = pan.append(ipan)
    print('saving pan into ', ofile)
    print('totan entries ', len(pan))
    pan.to_csv(ofile+'.csv')
    pan.to_hdf(ofile+'.h5', 'df', mode='w', format='table', data_columns=True)
    return pmaps, pan


def run_file(ifile, nevents=-1, **opts):

    pmaps = get_pmaps(ifile, nevents, **opts)

    fselect_pmap = opts.pop('select_pmap', select_pmap)
    cpmaps = select_pmaps(pmaps, fselect_pmap)

    fselect_s1 = opts.pop('select_s1', select_s1)
    fselect_s2 = opts.pop('select_s2', select_s2)
    pan = ana_pmaps(cpmaps, fselect_s1, fselect_s2, **opts)

    return cpmaps, pan


def get_pmaps(ifile, nevents, **opts):
    """ From an input file, read nevents and produce 'PMaps'
    """

    h5f = tb.open_file(ifile)
    pmtcwfs = h5f.root.RD.pmtcwf
    sipmrwfs = h5f.root.RD.sipmrwf

    pmtdb = db.DataPMT()

    if (nevents == -1):
        nevents = len(pmtcwfs)

    pmtcals = pmtdb['adc_to_pes'].values

    pmaps = []
    for ievt in range(nevents):
        wf = tamf.wfs_onewf(pmtcwfs[ievt], cals=pmtcals)
        siwfs = tamf.sipmwf_cwfs(sipmrwfs[ievt])
        ok, pmap = create_pmap(wf, siwfs, **opts)
        if (not ok):
            continue
        pmap.ievt = ievt
        pmaps.append(pmap)

    print('n events read ', nevents)
    print('n created pmaps ', len(pmaps))
    return pmaps


def create_pmap(wf, siwfs, **opts):
    """ create a 'PMap' from the corrected sum wf of all PMTs and the siwfs
    """

    ss = tamf.peaks_find(wf)
    if (len(ss) <= 0):
        print('No peaks!')
        return 0, None
    ssi = tamf.peaks_classify(wf, ss)
    qs = tamf.peaks_anode(ss, siwfs)
    # print(' pmap type ', ssi)
    # print(' pmap size ', [len(iss) for iss in ss])
    peaks = [Peak(times, wf[times[0]: times[-1]+1], qi, [], peaktype=itype)
             for itype, times, qi in zip(ssi, ss, qs)]
    pmap = PMap(peaks=peaks)
    # print(tamf.pmap_summary(pmap))
    return True, pmap


def select_pmaps(pmaps, select_pmap):
    """ select pmaps using the select_pmap conditional function
    """

    cpmaps = [pmap for pmap in pmaps if select_pmap(pmap)]
    print('n pmaps', len(pmaps))
    print('selected pmaps ', len(cpmaps))
    return cpmaps


def ana_pmaps(pmaps, select_s1, select_s2, **opts):
    """ produce a DataFrame with the S1 & S2 peaks information
    forma list of pmaps. Select the S1 and S2 with the select_s1 and
    select_s2 conditional functions
    """

    sipmdb = db.DataSiPM()
    xs = sipmdb['X']
    ys = sipmdb['Y']
    q0 = opts.pop('q0ano', q0ano)

    pan = {}
    labels = ['evt', 'ns1', 'ns2', 'is1', 'is2',
              't0', 't0wd', 't0cat', 't0ano', 't0icat',
              't', 'twd', 'tcat', 'tano', 'ticat',
              'x', 'y', 'ex', 'ey', 'exy']
    for label in labels:
        pan[label] = []

    for pmap in pmaps:
        s1s = tamf.pmap_peaks_with_condition(pmap, select_s1)
        ns1 = len(s1s)
        s2s = tamf.pmap_peaks_with_condition(pmap, select_s2)
        ns2 = len(s2s)
        for is1, s1 in enumerate(s1s):
            for is2, s2 in enumerate(s2s):
                pan['evt'].append(pmap.ievt)
                pan['ns1'].append(ns1)
                pan['ns2'].append(ns2)
                pan['is1'].append(is1)
                pan['is2'].append(is2)
                pan['t0'].append(s1.tmin+0.5*s1.width)
                pan['t0wd'].append(s1.width)
                pan['t0cat'].append(s1.cathode_integral)
                pan['t0ano'].append(s1.anode_integral)
                pan['t0icat'].append(s1.peakmax[1])
                pan['t'].append(s2.tmin+0.5*s2.width)
                pan['twd'].append(s2.width)
                pan['tcat'].append(s2.cathode_integral)
                pan['tano'].append(s2.anode_integral)
                pan['ticat'].append(s2.peakmax[1])
                point, epoint = tamf.peak_barycenter(s2, xs, ys, q0)
                pan['x'].append(point[1])
                pan['y'].append(point[2])
                pan['ex'].append(epoint[0])
                pan['ey'].append(epoint[1])
                pan['exy'].append(epoint[2])

    dpan = pd.DataFrame(pan)
    print('average s1', np.mean(dpan['ns1']))
    print('average s2', np.mean(dpan['ns2']))
    print('entries', len(dpan))
    return dpan


if __name__ == 'main':
    run(ifiles, ofile, nevents, **opts)
