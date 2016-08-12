"""
DIOMIRA
JJGC Agusut 2016

What SIERPE does:
1) Reads a Nh5 PRD (pre-raw data) containing MC waveforms for the 12 PMTs of the EP.
   Each waveform contains number of PEs in bins of 1 ns.
2) Convolves the PE waveform with the response of the FEE electronics.
3) Decimates the waveform, simulating the effect of the DAQ sampling (25 ns bins)
4) Writes a Nh5 dst with the new data 
5) Adds to Nh5 a new table with metadata
"""

from __future__ import print_function
from Util import *
from PlotUtil import *
from Nh5 import *
from cities import diomira

#import SierpeConfig as CP
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
    row['time_step'] = FP.time_step
    row['time_daq'] = FP.time_DAQ
    row['freq_LPF'] = FP.freq_LPF
    row['freq_HPF'] = 1./(2*pi*FP.R*FP.C)
    row['LSB'] = FP.LSB
    row['volts_to_adc'] = FP.voltsToAdc/volt
    row['noise_fee_rms'] = FP.NOISE_FEE
    row['noise_adc'] = FP.NOISE_ADC
    
    row.append()
    

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

def simulate_sipm_response(event_number,sipmrd_):
    """
    For the moment use a dummy rutne that simply copies the sipm EARRAY
    """
    rdata = []

    for j in range(sipmrd_.shape[1]):
        logging.debug("-->SiPM number ={}".format(j))
        rdata.append(sipmrd_[event_number, j])
    return np.array(rdata)

def copy_pmt(event_number,pmtrd_):
    """
    For the moment use a dummy rutne that simply copies the sipm EARRAY
    """
    rdata = []

    for j in range(pmtrd_.shape[1]):
        rdata.append(pmtrd_[event_number, j])
    return np.array(rdata)

def simulate_pmt_response(event_number,pmtrd_):
    """
    Sensor Response
    Given a signal in PE (photoelectrons in bins of 1 ns) and the response function of 
    for a single photoelectron (spe) and front-end electronics (fee)
    this function produces the PMT raw data (adc counts bins 25 ns)

    pmtrd_ dataset that holds the PMT PE data for each PMT
    pmtrd25 dataset to be created with adc counts, bins 25 ns after convoluting with electronics
    """
  
    rdata = []

    for j in range(pmtrd_.shape[1]):
        logging.debug("-->PMT number ={}".format(j))
                
        pmt = pmtrd_[event_number, j] #waveform for event event_number, PMT j
        
        fee = FE.FEE(C=FP.C12[j],R= FP.R, f=FP.freq_LPF, RG=FP.V_GAIN) 
        spe = SP.SPE(pmt_gain=FP.PMT_GAIN,x_slope = 5*ns,x_flat = 1*ns)
    
        signal_PMT = spe.SpePulseFromVectorPE(pmt) #PMT response

        #Front end response to PMT pulse (in volts)
        signal_fee = fee.FEESignal(signal_PMT, noise_rms=FP.NOISE_FEE) 

        #Signal out of DAQ
        signal_daq = fee.daqSignal(signal_fee, noise_rms=0)

        rdata.append(signal_daq)
    return np.array(rdata)

def usage():
    """
    Usage of program
    """
    print("""
        Usage: python (run) DIOMIRA [args]
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

    

        
     