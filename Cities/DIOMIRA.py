"""
DIOMIRA
JJGC August-October 2016

What DIOMIRA does:
1) Reads a MCRD file containing MC waveforms for the 12 PMTs of the EP.
   Each waveform contains number of PEs in bins of 1 ns.
2) Convolves the PE waveform with the response of the FEE electronics.
3) Decimates the waveform, simulating the effect of the DAQ sampling (25 ns bins)
4) Writes a RWF file with the new data and adds the FEE simulation parameters as metadata
"""

from __future__ import print_function
from Util import *
from LogConfig import *
from Configure import *
from Nh5 import *
import FEParam as FP
import SPE as SP
import FEE2 as FE

import tables
from time import time
import wfmFunctions as wfm
import sensorFunctions as sns
import pandas as pd

from RandomSampling import NoiseSampler as SiPMsNoiseSampler
#------


"""

DIOMIRA
ChangeLog:

26.9

Changed types of PMTRWF, SIPMRWF and PMTTWF to Float32 for
    (compatibility with ART/GATE)

Do not store EPMT and ESIPM (can be computed on the fly)

Change sign of pmtrwf to negative (as produced by the DAQ)

28.9 add cython

29.9 changed the way true waveforms are treated.
before: --> full waveform (in bins of 25 ns) in an EArray
now--> ZS waveform rebinned at 1 mus in a Table
(adavantages:) faster processing less space

01.10 moved loop definition to Configure.py and added index to TWF table

11.10 introduced SiPM noise. SiPMs' true waveforms stored under TWF group.
Some variables, classes and functions renamed for clarity.

"""
def FEE_param_table(fee_table):
    """
    Stores the parameters of the EP FEE simulation
    """
    row = fee_table.row
    row['offset'] = FP.offset
    row['pmt_gain'] = FP.PMT_GAIN
    row['V_gain'] = FP.V_GAIN
    row['R'] = FP.R
    row['C12'] = FP.C12
    row['AC'] = FP.AC
    row['time_step'] = FP.time_step
    row['time_daq'] = FP.time_DAQ
    row['freq_LPF'] = FP.freq_LPF
    row['freq_HPF'] = 1./(2*pi*FP.R*FP.C)
    row['LSB'] = FP.LSB
    row['volts_to_adc'] = FP.voltsToAdc/volt
    row['noise_fee_rms'] = FP.NOISE_FEE
    row['noise_adc'] = FP.NOISE_ADC

    row.append()


def store_wf(event, table, WF):
    """
    Store WF in table
    """
    row = table.row
    for isens,wf in WF.items():
        for t,e in zip(wf.time_mus, wf.ene_pes):
            row['event'] = event
            row['ID'] = isens
            row['time_mus'] = t
            row['ene_pes'] = e
            row.append()
    table.flush()

def rebin_twf(t, e, stride = 40):
    """
    rebins the a waveform according to stride
    The input waveform is a vector such that the index expresses time bin and the
    contents expresses energy (e.g, in pes)
    The function returns a DataFrame. The time bins and energy are rebinned according to stride
    """

    n = len(t)/int(stride)
    r = len(t)%int(stride)

    lenb = n
    if r > 0:
        lenb = n+1

    T = np.zeros(lenb,dtype=np.float32)
    E = np.zeros(lenb,dtype=np.float32)

    j=0
    for i in range(n):
        E[i] = np.sum(e[j:j+stride])
        T[i] = np.mean(t[j:j+stride])
        j+= stride

    if r > 0:
        E[n] = np.sum(e[j:])
        T[n] = np.mean(t[j:])

    return T,E


def pmt_twf_signal(event_number,pmtrd, stride):
    """
    1) takes pmtrd
    2) Performs ZS
    3) Rebins resulting wf according to stride
    """

    rdata = {}

    for j in range(pmtrd.shape[1]):
        logger.debug("-->PMT number ={}".format(j))

        energy_pes = pmtrd[event_number, j] #waveform for event event_number, PMT j
        time_mus = np.arange(pmtrd.shape[2])*ns/mus

        twf_zs = wfm.wf_thr(wfm.wf2df(time_mus,energy_pes),0.5)
        time_mus, ene_pes = rebin_twf(twf_zs.time_mus.values,twf_zs.ene_pes.values,stride)
        if not time_mus.any(): continue
        twf = wfm.wf2df(time_mus, ene_pes)

        logger.debug("-->len(twf) ={}".format(len(twf)))

        rdata[j] = twf
    return rdata


def sipm_twf_signal(event_number,sipmrd):
    '''
        Removes zeros from SiPM RD.
    '''
    out = {}
    for index,wfm in enumerate(sipmrd[event_number]):
        time_mus = np.where( wfm > 0. )[0]
        if not time_mus.any(): continue
        ene_pes = wfm[time_mus]
        out[index] = pd.DataFrame( {'time_mus':time_mus, 'ene_pes':ene_pes} )
    return out


def simulate_sipm_response(event_number,sipmrd_,sipms_noise_sampler):
    """
    Add noise with the SiPMNoiseSampler class and return the noisy waveform.
    """
    return sipmrd_[event_number] + sipms_noise_sampler.Sample()


def simulate_pmt_response(event_number,pmtrd_):
    """
    Sensor Response
    Given a signal in PE (photoelectrons in bins of 1 ns) and the response function of
    for a single photoelectron (spe) and front-end electronics (fee)
    this function produces the PMT raw data (adc counts bins 25 ns)

    pmtrd_ dataset that holds the PMT PE data for each PMT
    pmtrd25 dataset to be created with adc counts, bins 25 ns
    after convoluting with electronics
    """

    rdata = []

    for j in range(pmtrd_.shape[1]):
        logger.debug("-->PMT number ={}".format(j))

        pmt = pmtrd_[event_number, j] #waveform for event event_number, PMT j

        fee = FE.FEE(C=FP.C12[j],R= FP.R, f=FP.freq_LPF, RG=FP.V_GAIN)
        spe = SP.SPE(pmt_gain=FP.PMT_GAIN,x_slope = 5*ns,x_flat = 1*ns)

        signal_PMT = spe.SpePulseFromVectorPE(pmt) #PMT response

        #Front end response to PMT pulse (in volts)
        signal_fee = fee.FEESignal(signal_PMT, noise_rms=FP.NOISE_FEE)

        #Signal out of DAQ
        #positive signal convention
        #signal_daq = fee.daqSignal(signal_fee, noise_rms=0) - FP.offset
        #negative signals convention!

        signal_daq = FP.offset -fee.daqSignal(signal_fee, noise_rms=0)

        rdata.append(signal_daq)
    return np.array(rdata)


def DIOMIRA(argv):
    """
    Diomira driver
    """
    DEBUG_LEVEL, INFO, CYTHON, CFP = configure(argv[0],argv[1:])

    if INFO:

        print("""
        DIOMIRA:
         1. Reads an MCRD file produced by art/centella, which stores MCRD
         waveforms for PMTs (bins of 1 ns)
        and the SiPMs (bins of 1 mus)


        2. Simulates the response of the energy plane in the PMTs MCRD,
        and produces both RWF and TWF:
        see: http://localhost:8931/notebooks/Nh5-Event-Model.ipynb#Reconstructed-Objects


        3. Simulates the response of the tracking plane in the SiPMs MCRD and outputs
            SiPM RWF (not yet implemented, for the time being simply copy the MCRD)

        4. Add a table describing the FEE parameters used for simulation

        5. Copies the tables on geometry, detector data and MC


        """)
        FP.print_FEE()

    PATH_IN =CFP['PATH_IN']
    PATH_OUT =CFP['PATH_OUT']
    FILE_IN =CFP['FILE_IN']
    FILE_OUT =CFP['FILE_OUT']
    PATH_DB =CFP['PATH_DB']
    FIRST_EVT =CFP['FIRST_EVT']
    LAST_EVT =CFP['LAST_EVT']
    RUN_ALL =CFP['RUN_ALL']
    CLIB =CFP['CLIB']
    CLEVEL =CFP['CLEVEL']
    NOISE_CUT_FRACTION = CFP['NOISE_CUT_FRACTION']
    NEVENTS = LAST_EVT - FIRST_EVT

    print('Debug level = {}'.format(DEBUG_LEVEL))

    print("input path ={}; output path = {}; file_in ={} file_out ={}".format(
        PATH_IN,PATH_OUT,FILE_IN, FILE_OUT))

    print("first event = {} last event = {} nof events requested = {} ".format(
        FIRST_EVT,LAST_EVT,NEVENTS))

    print("Compression library = {} Compression level = {} ".format(
        CLIB,CLEVEL))

    # open the input file
    with tables.open_file("{}/{}".format(PATH_IN,FILE_IN), "r") as h5in:
        # access the PMT raw data in file

        pmtrd_ = h5in.root.pmtrd
        sipmrd_ = h5in.root.sipmrd

        #pmtrd_.shape = (nof_events, nof_sensors, wf_length)

        NPMT = pmtrd_.shape[1]
        NSIPM = sipmrd_.shape[1]
        PMTWL = pmtrd_.shape[2]
        #PMTWL_FEE = int((PMTWL+1)/FP.time_DAQ) #old format
        PMTWL_FEE = int(PMTWL/FP.time_DAQ)
        SIPMWL = sipmrd_.shape[2]
        NEVENTS_DST = pmtrd_.shape[0]

        print("nof PMTs = {} nof  SiPMs = {} nof events in input DST = {} ".format(
        NPMT,NSIPM,NEVENTS_DST))

        print("lof SiPM WF = {} lof PMT WF (MC) = {} lof PMT WF (FEE) = {}".format(
        PMTWL,SIPMWL,PMTWL_FEE))


        #access the geometry and the sensors metadata info

        geom_t = h5in.root.Detector.DetectorGeometry
        pmt_t = h5in.root.Sensors.DataPMT
        sipm_t = h5in.root.Sensors.DataSiPM
        mctrk_t = h5in.root.MC.MCTracks

        # Map of the SiPMs' sensorID to the index used by tables
        index_map = { sipm_t[i][0] : i for i in range(sipm_t.shape[0]) }
        # Create instance of the noise sampler
        sipms_noise_sampler_ = SiPMsNoiseSampler(PATH_DB+"/NoiseSiPM_NEW.dat",index_map,SIPMWL,True)
        sipms_noise_thresholds_ = SiPMsNoiseSampler.ComputeThresholds(NOISE_CUT_FRACTION)
        # open the output file
        with tables.open_file("{}/{}".format(PATH_OUT,FILE_OUT), "w",
            filters=tables.Filters(complib=CLIB, complevel=CLEVEL)) as h5out:

            # create a group to store MC data
            mcgroup = h5out.create_group(h5out.root, "MC")
            # copy the mctrk table
            mctrk_t.copy(newparent=mcgroup)

            # create a group  to store geom data
            detgroup = h5out.create_group(h5out.root, "Detector")
            # copy the geom table
            geom_t.copy(newparent=detgroup)

            # create a group  store sensor data
            sgroup = h5out.create_group(h5out.root, "Sensors")
            # copy the pmt table
            pmt_t.copy(newparent=sgroup)
            # copy the sipm table
            sipm_t.copy(newparent=sgroup)

            # create a table to store Energy plane FEE data and hang it from MC group
            fee_table = h5out.create_table(mcgroup, "FEE", FEE,
                          "EP-FEE parameters",tables.Filters(0))


            # create a group to store True waveform data
            twfgroup = h5out.create_group(h5out.root, "TWF")

            # create a table to store true waveform (zs, rebinned)
            pmt_twf_table = h5out.create_table( twfgroup, "PMT", SENSOR_WF, "Store for PMTs TWF",
                                                tables.Filters(complib=CLIB, complevel=CLEVEL) )

            sipm_twf_table = h5out.create_table( twfgroup, "SiPM", SENSOR_WF, "Store for SiPMs TWF",
                                                 tables.Filters(complib=CLIB, complevel=CLEVEL) )

            #and index in event column
            pmt_twf_table.cols.event.create_index()
            sipm_twf_table.cols.event.create_index()

            # fill FEE table
            FEE_param_table(fee_table)

            # create a group to store RawData
            rgroup = h5out.create_group(h5out.root, "RD")

            # create an extensible array to store the RWF waveforms
            pmtrwf = h5out.create_earray(h5out.root.RD, "pmtrwf",
                                    atom=tables.Float32Atom(),
                                    shape=(0, NPMT, PMTWL_FEE),
                                    expectedrows=NEVENTS_DST)

            sipm_rwf_table = h5out.create_table( rgroup, "sipmrwf", SENSOR_WF, "Store for SiPMs RWF",
                                                 tables.Filters(complib=CLIB, complevel=CLEVEL) )

            #LOOP
            first_evt, last_evt = define_event_loop(FIRST_EVT,LAST_EVT,NEVENTS,NEVENTS_DST,RUN_ALL)

            t0 = time()
            for i in range(first_evt,last_evt):
                logger.info("-->event number ={}".format(i))

                # supress zeros in MCRD and rebins the ZS function in 1 mus bins

                rebin = int(1*mus/1*ns)  #rebins zs function in 1 mus bin

                #list with zs twf
                truePMT  =  pmt_twf_signal(i,pmtrd_, rebin)
                trueSiPM = sipm_twf_signal(i,sipmrd_)

                #store in table
                store_wf(i, pmt_twf_table, truePMT)
                store_wf(i, sipm_twf_table, trueSiPM)

                #simulate PMT response and return an array with RWF
                dataPMT = simulate_pmt_response(i,pmtrd_)

                #convert to float
                dataPMT.astype(float)

                #append to EVECTOR
                pmtrwf.append(dataPMT.reshape(1, NPMT, PMTWL_FEE))


                #simulate SiPM response and return an array with RWF
                #convert to float, append to EVector

                dataSiPM = simulate_sipm_response(i,sipmrd_,sipms_noise_sampler_)
                dataSiPM.astype(float)

                zs_wfms = sns.sensor_wise_zero_suppresion(dataSiPM,thresholds)

                store_wf( i, sipm_rwf_table, zs_wfms )

            t1 = time()
            pmtrwf.flush()

            print("DIOMIRA has run over {} events in {} seconds".format(i, t1-t0))
    print("Leaving Diomira. Safe travels!")

if __name__ == '__main__':
    #import cProfile

    #cProfile.run('DIOMIRA(sys.argv)', sort='time')
    DIOMIRA(sys.argv)
