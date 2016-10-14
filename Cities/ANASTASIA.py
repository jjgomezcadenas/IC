"""
ANASTASIA
GML October 2016

What ANASTASIA does:
1) Reads a hdf5 file containing the PMT's CWF and the SiPMs' RWF expressed in ADC counts.
2) Expresses the waveforms in pes.
3) Applies zero-suppression.
4) Writes back in the input file the ZS waveforms.
"""

from __future__ import print_function
from Util import *
from LogConfig import *
from Configure import *
from Nh5 import SENSOR_WF
from FEParam import NOISE_ADC

import tables
from time import time
import wfmFunctions as wfm
import pandas as pd
#------

"""

ANASTASIA
ChangeLog:

14.10 First version.

"""

def scale_to_pes( sens_wf, sensdf ):
    '''
        Transform the ene_pes field to pes for each sensor.
    '''
    return { key : pd.DataFrame( {'time_mus' : df.time_mus, 'ene_pes' : df.ene_pes * sensdf['adc_to_pes'][key]} ) for key,df in sens_wf.items() }

def ANASTASIA(argv):
    """
    ANASTASIA driver
    """
    DEBUG_LEVEL, INFO, CYTHON, CFP = configure(argv[0],argv[1:])

    if INFO:
        print(__doc__)

    PATH_IN   = CFP['PATH_IN']
    PATH_OUT  = CFP['PATH_OUT']
    FILE_IN   = CFP['FILE_IN']
    PATH_DB   = CFP['PATH_DB']
    FIRST_EVT = CFP['FIRST_EVT']
    LAST_EVT  = CFP['LAST_EVT']
    RUN_ALL   = CFP['RUN_ALL']
    CLIB      = CFP['CLIB']
    CLEVEL    = CFP['CLEVEL']
    NEVENTS   = LAST_EVT - FIRST_EVT

    ZS_METHOD_SIPMS = CFP['ZS_METHOD_SIPMS']
    NOISE_CUT_PMTS  = CFP['NOISE_CUT_PMTS']
    NOISE_CUT_SIPMS = CFP['NOISE_CUT_SIPMS']

    logger.info('Debug level = {}'.format(DEBUG_LEVEL))
    logger.info("input file = {}/{}".format(PATH_IN,FILE_IN))
    logger.info("path to database = {}".format(PATH_DB))
    logger.info("first event = {}; last event = {} nof events requested = {} ".format(FIRST_EVT,LAST_EVT,NEVENTS))
    logger.info("Compression library = {} Compression level = {} ".format(CLIB,CLEVEL))
    logger.info("ZS method PMTS  = {}. Cut value = {}".format("RMS CUT",NOISE_CUT_PMTS))
    logger.info("ZS method SIPMS = {}. Cut value = {}".format(ZS_METHOD,NOISE_CUT_SIPMS))

    # open the input file
    with tables.open_file("{}/{}".format(PATH_IN,FILE_IN), "r+") as h5in:
        # access the PMT raw data in file

        pmttwf, sipmtwf, pmtrwf, sipmrwf, pmtdf, sipmdf, gdf = tbl.get_vectors(h5f)

        NEVT, NPMT , PMTWL  = pmtrwf.shape
        NEVT, NSIPM, SIPMWL = sipmrwf.shape

        logger.info("#PMTs = {}; #SiPMs = {}; #events in DST = {}".format(NPMT,NSIPM,NEVT))
        logger.info("PMT WFL = {}; SiPM WFL = {}".format(PMTWL,SIPMWL))

        # Create instance of the noise sampler and compute noise thresholds
        sipms_noise_sampler_    = SiPMsNoiseSampler(PATH_DB+"/NoiseSiPM_NEW.dat",sipmdf,SIPMWL,False)
        pmts_noise_thresholds_  = np.ones(NPMT) * NOISE_ADC * NOISE_CUT_PMTS
        sipms_noise_thresholds_ = sipms_noise_sampler_.ComputeThresholds(NOISE_CUT_SIPMS) if ZS_METHOD == 'FRACTION' else np.ones(NSIPM) * NOISE_CUT_SIPMS

        zsgroup = h5in.create_group(h5in.root, "ZS")
        pmt_zs_table  = h5in.create_table( zsgroup, "PMT", SENSOR_WF, "Store for PMTs ZSWF",
                                           tables.Filters(complib=CLIB, complevel=CLEVEL) )

        sipm_zs_table = h5in.create_table( zsgroup, "SiPM", SENSOR_WF, "Store for SiPMs ZSWF",
                                           tables.Filters(complib=CLIB, complevel=CLEVEL) )

        #and index in event column
        pmt_zs_table.cols.event.create_index()
        sipm_zs_table.cols.event.create_index()

        first_evt, last_evt = define_event_loop(FIRST_EVT,LAST_EVT,NEVENTS,NEVENTS_DST,RUN_ALL)

        t0 = time()
        for i in range(first_evt,last_evt):
            logger.info("-->event number ={}".format(i))

            dataPMT = wfm.sensor_wise_zero_suppresion(pmtrwf[i],pmt_noise_thresholds_)
            tbl.store_wf( i, pmt_zs_table, scale_to_pes(dataPMT,pmtdf) )

            dataSiPM = wfm.sensor_wise_zero_suppresion(sipmrwf[i],sipms_noise_thresholds_)
            tbl.store_wf( i, sipm_zs_table, scale_to_pes(dataSiPM,sipmdf) )

        t1 = time()

        print("ANASTASIA has run over {} events in {} seconds".format(i, t1-t0))
    print("Leaving ANASTASIA. Safe travels!")

if __name__ == '__main__':
    #import cProfile

    #cProfile.run('ANASTASIA(sys.argv)', sort='time')
    ANASTASIA(sys.argv)
