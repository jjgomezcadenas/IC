"""
ANASTASIA
GML October 2016

What ANASTASIA does:
1) Reads a hdf5 file containing the PMT's CWF and the SiPMs' RWF in ADC counts.
2) Applies zero-suppression.
3) Expresses the waveforms in pes.
4) Writes back in the input file the ZS waveforms as tables.
"""

from __future__ import print_function
from Util import *
from LogConfig import *
from Configure import *
from Nh5 import SENSOR_WF
from FEParam import NOISE_ADC

import tables
from time import time
import sensorFunctions as snf
import wfmFunctions as wfm
import tblFunctions as tbl
import pandas as pd

from RandomSampling import NoiseSampler as SiPMsNoiseSampler
#------

'''

ANASTASIA
ChangeLog:

14.10 First version.

'''

def scale_to_pes( sens_wf, sensdf ):
    '''
        Transform the ene_pes field to pes for each sensor.
    '''
    return { key : pd.DataFrame( {'time_mus' : df.time_mus, 'ene_pes' : -df.ene_pes / sensdf['adc_to_pes'][key]} ) for key,df in sens_wf.items() }

def ANASTASIA(argv):
    '''
        ANASTASIA driver
    '''
    DEBUG_LEVEL, INFO, CYTHON, CFP = configure(argv[0],argv[1:])

    if INFO:
        print(__doc__)

    PATH_IN   = CFP['PATH_IN']
    FILE_IN   = CFP['FILE_IN']
    PATH_OUT  = CFP['PATH_OUT']
    FILE_OUT  = CFP['FILE_OUT']
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
    logger.info("ZS method SIPMS = {}. Cut value = {}".format(ZS_METHOD_SIPMS,NOISE_CUT_SIPMS))

    # open the input file
    with tables.open_file("{}/{}".format(PATH_IN,FILE_IN), "r+") as h5in:
        # access the PMT raw data in file

        pmtcwf   = h5in.root.RD.pmtcwf
        sipmrwf  = h5in.root.RD.sipmrwf
        mcdata   = h5in.root.MC.MCTracks
        geodata  = h5in.root.Detector.DetectorGeometry
        pmtdata  = h5in.root.Sensors.DataPMT
        sipmdata = h5in.root.Sensors.DataSiPM
        pmtdf    = snf.read_data_sensors(pmtdata)
        sipmdf   = snf.read_data_sensors(sipmdata)

        NEVT, NPMT , PMTWL  = pmtcwf.shape
        NEVT, NSIPM, SIPMWL = sipmrwf.shape

        logger.info("#PMTs = {}; #SiPMs = {}; #events in DST = {}".format(NPMT,NSIPM,NEVT))
        logger.info("PMT WFL = {}; SiPM WFL = {}".format(PMTWL,SIPMWL))

        pmt_adc_consts = -1.0/np.array(pmtdf['adc_to_pes']).reshape(NPMT,1)
        # Create instance of the noise sampler and compute noise thresholds
        sipms_noise_sampler_    = SiPMsNoiseSampler(PATH_DB+"/NoiseSiPM_NEW.dat",sipmdf,SIPMWL)
        pmts_noise_threshold_   = NOISE_ADC * np.mean(pmt_adc_consts) * NOISE_CUT_PMTS
        sipms_noise_thresholds_ = sipms_noise_sampler_.ComputeThresholds(NOISE_CUT_SIPMS,sipmdf = sipmdf) if ZS_METHOD_SIPMS == 'FRACTION' else np.ones(NSIPM) * NOISE_CUT_SIPMS

        with tables.open_file("{}/{}".format(PATH_OUT,FILE_OUT), "w") as h5out:
            # Create groups
            h5out.create_group(h5out.root, "MC")
            h5out.create_group(h5out.root, "Detector")
            h5out.create_group(h5out.root, "Sensors")
            h5out.create_group(h5out.root, "ZS")

            # Copy data from input
            mcdata  .copy(newparent=h5out.root.MC)
            geodata .copy(newparent=h5out.root.Detector)
            pmtdata .copy(newparent=h5out.root.Sensors)
            sipmdata.copy(newparent=h5out.root.Sensors)


            # Create tables for PMTs and SiPMs ZSWF
            pmt_zs_table  = h5in.create_table( h5out.root.ZS, "PMT", SENSOR_WF, "Store for PMTs ZSWF",
                                               tables.Filters(complib=CLIB, complevel=CLEVEL) )

            sipm_zs_table = h5in.create_table( h5out.root.ZS, "SiPM", SENSOR_WF, "Store for SiPMs ZSWF",
                                               tables.Filters(complib=CLIB, complevel=CLEVEL) )

            # Add index in event column
            pmt_zs_table .cols.event.create_index()
            sipm_zs_table.cols.event.create_index()

            first_evt, last_evt = define_event_loop(FIRST_EVT,LAST_EVT,NEVENTS,NEVT,RUN_ALL)

            t0 = time()
            for i in range(first_evt,last_evt):
                logger.info("-->event number ={}".format(i))

                # Integrate PMT plane in pes (not in time!)
                pmt_int_pes = (pmtcwf[i] * pmt_adc_consts).sum(axis=0)

                # sensor_wise_zero_suppresion returns a dictionary holding
                # the surviving sensors dataframes. Time in mus, ene in pes (PMTs) or adc (SiPMs).
                # They are converted just before storing them in the table.
                dataPMT  = wfm.zs_wf(pmt_int_pes,pmts_noise_threshold_,25*ns/mus)
                tbl.store_wf( i, pmt_zs_table, { 12 : dataPMT } )

                tmin, tmax = dataPMT.time_mus.min(),dataPMT.time_mus.max()

                dataSiPM = wfm.sensor_wise_zero_suppression(sipmrwf[i],sipms_noise_thresholds_,1.0)
                dataSiPM = wfm.in_window(dataSiPM,tmin,tmax)
                tbl.store_wf( i, sipm_zs_table, scale_to_pes(dataSiPM,sipmdf) )

            t1 = time()

        print("ANASTASIA has run over {} events in {} seconds".format(i, t1-t0))
    print("Leaving ANASTASIA. Safe travels!")

if __name__ == '__main__':
    #import cProfile

    #cProfile.run('ANASTASIA(sys.argv)', sort='time')
    ANASTASIA(sys.argv)
