"""
ISIDORA
JJGC Agusut 2016

What ISIDORA does:
1) Reads a RWF file written by DIOMIRA 
2) Performs DBLR
3) Write the corrected waveforms (CWF) to the file as new Evectors. 
4) Computes the energy of the CWF and adds it to the file
"""

from __future__ import print_function
from Util import *
from PlotUtil import *
from Nh5 import *
from cities import isidora

import numpy as np

import FEParam as FP
import SPE as SP
import FEE2 as FE

import tables
import pandas as pd
import logging
import getopt



"""
Code
"""


def energy_pes(event_number, sensord):
    """
    Sum the WFs of PMTs and SiPMs (MC) and store the total energy in PES
    """     
    rdata = []

    for j in range(sensord.shape[1]):
        swf = sensord[event_number, j]
        ene = np.sum(swf)
        rdata.append(ene)
        
    return np.array(rdata) 

def BLR(event_number,pmtrd25):
    """
    Peform Base line Restoration
    """
    ene_pmt =np.zeros(len(CP.CFGP['PMTS']), dtype=np.int32)
    ipmt = 0
    PMTWF ={}
    for j in CP.CFGP['PMTS']:
        logging.debug("-->PMT number ={}".format(j))

        pmtrd = pmtrd_[event_number, j] #waveform for event event_number, PMT j
        
        if CP.CPLT['plot_RAW'] == True:
            print("RAW signal")
            plot_signal(CP.signal_t,pmtrd, title = 'signal RAW', 
                    signal_start=0, signal_end=CP.CFGP['LEN_WVF25'], 
                    units='adc')

        #Deconvolution
        fee = FE.FEE(C=FP.C12[j],R= FP.R, f=FP.freq_LPF, RG=FP.V_GAIN)
        signal_inv_daq = fee.InverseSignalDAQ(CP.signal_t)  #inverse function
        coef = signal_inv_daq[10]  #accumulator coefficient
        

        if CP.CPLT['plot_signal_inv'] :   
            plot_signal(CP.signal_t/ns,signal_inv_daq,
                title = 'Inverse DAQ', 
                signal_start=0*ns, signal_end=10*ns, 
                units='')
            print("inverse coef fee: = {}".format(coef))

        # print "calling MauDeconv"
   
        signal_blr, eadc = DB.BLR(pmtrd, coef, n_sigma = CP.CFGP['NSIGMA'], 
                            NOISE_ADC=FP.NOISE_ADC, thr2=FP.NOISE_ADC/4., thr3 = FP.NOISE_ADC/2.,
                            plot=CP.CPLT['plot_BLR'])
        if CP.CPLT['plot_DEC'] == True:
            print("DBLR signal")
            plot_signal(CP.signal_t,signal_blr, title = 'signal BLR', 
                    signal_start=0, signal_end=CP.CFGP['LEN_WVF25'], 
                    units='adc')

        ene_pmt[ipmt] = eadc
        PMTWF[ipmt]=signal_blr
        ipmt+=1
    return pd.Series(ene_pmt), pd.DataFrame(PMTWF)
        
def usage():
    """
    Usage of program
    """
    print("""
        Usage: python (run) ISIDORA [args]
        where args are:
         -h (--help) : this text
         -i (--info) : print a text describing the invisible city of DIOMIRA
         -d (--debug) : can be set to 'DEBUG','INFO','WARNING','ERROR'
         -c (--cfile) : full path to a configuration file
         
         example of configuration file 

        #Header is not read
        Names of parameters (comma separated)
        Values of parameters (comma separated)
        
        The parameters for DIOMIRA are:
        PATH_IN,PATH_OUT,FILE_IN,FILE_OUT,FIRST_EVT,LAST_EVT,RUN_ALL 

        The parameters are self-explaining. 
        RUN_ALL is used to decide whether to run all the events in the file 
        in case that the total number of events requested (LAST_EVT-FIRST_EVT) 
        exceeds the number of events in the DST file. If RUN_ALL is set to 1 (True), 
        the script will run over all elements in the DST, 
        otherwise it will exit with a warning.

        """)
def configure(argv):
    """
    reads arguments from the command line and configures job
    """
    
    global DEBUG, PATH_IN, PATH_OUT, FILE_IN, FILE_OUT
    global  FIRST_EVT, LAST_EVT,NEVENTS, RUN_ALL, INFO
    
    DEBUG='INFO'
    INFO = True
    cfile =''
    try:
        opts, args = getopt.getopt(argv, "hidc:", ["help","info","debug","cfile"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-d", "--debug"):
            DEBUG = arg
        elif opt in ("-i", "--info"):
            INFO = arg
        elif opt in ("-c", "--cfile"):
            cfile = arg
 
    if DEBUG == 'ERROR' or DEBUG == 'error' or DEBUG == 'e':
        logging.basicConfig(level=logging.ERROR)
    elif DEBUG == 'DEBUG' or DEBUG == 'debug' or DEBUG == 'd':
        logging.basicConfig(level=logging.DEBUG)
    elif DEBUG == 'WARNING' or DEBUG == 'warning' or DEBUG == 'w':
        logging.basicConfig(level=logging.WARNING)
    elif DEBUG == 'INFO' or DEBUG == 'info' or DEBUG == 'i':
        logging.basicConfig(level=logging.INFO)
    else:
        print("value of debug option not defined")
        usage()

    if cfile == '':
        print("Path to configuration file not given")
        usage()
        sys.exit()

    CFP =pd.read_csv(cfile,skiprows=1)
    print("""
        Configuration parameters \n 
        {}
        """.format(CFP))

    PATH_IN=CFP['PATH_IN'][0] 
    PATH_OUT=CFP['PATH_OUT'][0]
    FILE_IN=CFP['FILE_IN'][0]
    FILE_OUT=CFP['FILE_OUT'][0]
    FIRST_EVT=CFP['FIRST_EVT'][0]
    LAST_EVT=CFP['LAST_EVT'][0]
    RUN_ALL=CFP['RUN_ALL'][0]
    NEVENTS = LAST_EVT -  FIRST_EVT

        
if __name__ == '__main__':
    configure(sys.argv[1:])

    if INFO:
        print(diomira)

    wait()
    
    print("""
        DIOMIRA:
        1. Reads an Nh5 file produced by art/centella, which stores the
            pre-raw data (PRD) for the PMTs and SiPMs waveforms, as well as
            data on geometry, sensors and MC.

        2. Simulates the response of the energy plane in the PMTs PRD, and
            outputs PMT Raw-Data (RD), e.g., waveforms in bins of 25 ns
            which correspond to the output of the EP FEE.

        3. Re-formats the data on detector geometry and sensors as pandas
            dataframes (PDF) and Series (PS)

        4. Add a PS describing the FEE parameters used for simulation

        5. Copies the MC table and the SiPM PRD (the simulation of SiPM response
        is pending)

        """)
    FP.print_FEE()
    wait()


    print("input path ={}; output path = {}; file_in ={} file_out ={}".format(
        PATH_IN,PATH_OUT,FILE_IN, FILE_OUT))

    print("first event = {} last event = {} nof events requested = {} ".format(
        FIRST_EVT,LAST_EVT,NEVENTS))

    # open the input file 
    with tables.open_file("{}/{}".format(PATH_IN,FILE_IN), "r+") as h5in: 
        # access the PMT raw data in file 
        pmtrd_ = h5in.root.pmtrd
        sipmrd_ = h5in.root.sipmrd

        #pmtrd_.shape = (nof_events, nof_sensors, wf_length)
        NPMT = pmtrd_.shape[1]
        NSIPM = sipmrd_.shape[1]
        PMTWL = pmtrd_.shape[2] 
        PMTWL_FEE = int((PMTWL+1)/FP.time_DAQ)
        SIPMWL = sipmrd_.shape[2]
        NEVENTS_DST = pmtrd_.shape[0]

        print("nof PMTs = {} nof  SiPMs = {} nof events in input DST = {} ".format(
        NPMT,NSIPM,NEVENTS_DST))

        print("lof SiPM WF = {} lof PMT WF (MC) = {} lof PMT WF (FEE) = {}".format(
        PMTWL,SIPMWL,PMTWL_FEE))

        wait()

        #access the geometry and the sensors metadata info

        geom_t = h5in.root.Detector.DetectorGeometry
        pmt_t = h5in.root.Sensors.DataPMT
        sipm_t = h5in.root.Sensors.DataSiPM
        mctrk_t = h5in.root.MC.MCTracks

        # #return a pandas DF for sensors and geometry
        # pmtdf = read_data_sensors(pmt_t)
        # sipmdf = read_data_sensors(sipm_t)
        # geodf = read_data_geom(geom_t)

        
        # open the output file 
        with tables.open_file("{}/{}".format(PATH_OUT,FILE_OUT), "w",
            filters=tables.Filters(complib="blosc", complevel=9)) as h5out:
 
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
                          "EP-FEE parameters",
                           tables.Filters(0))

            # fill table
            FEE_param_table(fee_table)

            # create a group to store RawData
            rgroup = h5out.create_group(h5out.root, "RD")
            
            # create an extensible array to store the RWF waveforms
            pmtrd = h5out.create_earray(h5out.root.RD, "pmtrd", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NPMT, PMTWL_FEE), 
                                    expectedrows=NEVENTS_DST)
            # pmtrd = h5out.create_earray(h5out.root.RD, "pmtrd", 
            #                         atom=tables.IntAtom(), 
            #                         shape=(0, NPMT, PMTWL), 
            #                         expectedrows=NEVENTS_DST)

            sipmrd = h5out.create_earray(h5out.root.RD, "sipmrd", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NSIPM, SIPMWL), 
                                    expectedrows=NEVENTS_DST)

            #create an extensible array to store the energy in PES of PMTs 
            epmt = h5out.create_earray(h5out.root.RD, "epmt", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NPMT), 
                                    expectedrows=NEVENTS_DST)

            # create an extensible array to store the energy in PES of SiPMs 
            esipm = h5out.create_earray(h5out.root.RD, "esipm", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NSIPM), 
                                    expectedrows=NEVENTS_DST)

            
            if NEVENTS > NEVENTS_DST and RUN_ALL == False:
                print("""
                Refusing to run: you have requested
                FIRST_EVT = {}
                LAST_EVT  = {}
                Thus you want to run over {} events
                but the size of the DST is {} events.
                Please change your choice or select RUN_ALL = TRUE
                to run over the whole DST when this happens
                """.format(FIRST_EVT,LAST_EVT,NEVENTS,NEVENTS_DST))
                sys.exit(0)
            elif  NEVENTS > NEVENTS_DST and RUN_ALL == True:
                FIRST_EVT = 0
                LAST_EVT = NEVENTS_DST
                NEVENTS = NEVENTS_DST

            #create a 2d array to store the true energy of the PMTs (pes)
            ene_pmt = np.zeros((NEVENTS,NPMT))

            #create a 2d array to store the true energy of the SiPMs (pes)
            ene_sipm = np.zeros((NEVENTS,NSIPM))

            for i in range(FIRST_EVT,LAST_EVT):
                print("-->event number ={}".format(i))
                logging.info("-->event number ={}".format(i))

                #simulate PMT response and return an array with new WF
                dataPMT = simulate_pmt_response(i,pmtrd_)
                #dataPMT = copy_pmt(i,pmtrd_)
                
                #append to PMT EARRAY
                pmtrd.append(dataPMT.reshape(1, NPMT, PMTWL_FEE))
                #pmtrd.append(dataPMT.reshape(1, NPMT, PMTWL))
                   
                #simulate SiPM response and return an array with new WF
                dataSiPM = simulate_sipm_response(i,sipmrd_)
                
                #append to SiPM EARRAY
                sipmrd.append(dataSiPM.reshape(1, NSIPM, SIPMWL))

                #fill ene_pmt vector
                enePMT = energy_pes(i, pmtrd_)
                #append to epmt EARRAY
                epmt.append(enePMT.reshape(1, NPMT))

                #fill ene_sipm vector
                eneSIPM = energy_pes(i, sipmrd_)
                esipm.append(eneSIPM.reshape(1, NSIPM))

            pmtrd.flush()
            sipmrd.flush()
            epmt.flush()
            esipm.flush()

            #store the DF with sensors and geom info
            # store_export = pd.HDFStore(PATH_OUT+FILE)
            # store_export.append('Sensors/DataPMT', pmtdf, data_columns=pmtdf.columns)
            # store_export.append('Sensors/DataSiPM', sipmdf, data_columns=pmtdf.columns)
            # store_export.append('Detector/DetectorGeometry', geodf)
            # store_export.append('EnergyPlane/FEE', feedf)
            # store_export.close()

    print("done!")

    

        
     