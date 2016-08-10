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

"""
Configuration parameters
"""


logging.basicConfig(level=logging.INFO)


PATH_IN ="/Users/jjgomezcadenas/Documents/Development/NEXT/data/Waveforms/"
PATH_OUT ="/Users/jjgomezcadenas/Documents/Development/NEXT/data/Waveforms/25ns/"

FILE = "WF_Tl_0.h5"

FIRST_EVT = 0
LAST_EVT = 10
NEVENTS = LAST_EVT - FIRST_EVT
RUN_ALL = True


"""
Code
"""

def get_column(pmta,ic):
    """
    access column ic of table pmta and returns column as an array
    """
    col =[]
    for i in range(pmta.shape[0]):
        col.append(pmta[i][ic])
    return np.array(col)
 
def read_data_sensors(sensor_table):
    """
    reads the sensors table and returns a data frame
    """
    pmta = sensor_table.read()
    PMT={}
    PMT['channel'] = get_column(pmta,0)
    PMT['active'] = get_column(pmta,1)
    PMT['x'] = get_column(pmta,2).T[0]
    PMT['y'] = get_column(pmta,2).T[1]
    PMT['gain'] = get_column(pmta,3)
    PMT['adc_to_pes'] = get_column(pmta,4)
        
    return pd.DataFrame(PMT)


def FEE_param():
    """
    Stores the parameters of the EP FEE simulation as a pd Series
    """
    fp = pd.Series([FP.offset,FP.PMT_GAIN,FP.V_GAIN,FP.R,FP.time_step,FP.time_DAQ,
                    FP.freq_LPF,1./(2*pi*FP.R*FP.C),FP.LSB,FP.voltsToAdc/volt,
                    FP.NOISE_FEE,FP.NOISE_ADC], 
                    index=['offset','pmt_gain','V_gain','R',
                                'time_step','time_daq','freq_LPF',
                                'freq_HPF','LSB','volts_to_adc',
                                'noise_fee_rms','noise_adc'])
    return fp

def read_data_geom(geom_t):
    """
    Reads the geom data en returns a PD Series
    """
        
    ga = geom_t.read()
    G ={}
    G = pd.Series([ga[0][0][0],ga[0][0][1],ga[0][1][0],ga[0][1][1],
                    ga[0][2][0],ga[0][2][1],ga[0][3]],
                    index=['xdet_min','xdet_max','ydet_min','ydet_max',
                            'zdet_min','zdet_max','R'])
    return G

def copy_sipm(event_number,sipmrd_):
    """
    Copies the sipm EARRAY
    """
    rdata = []

    for j in range(sipmrd_.shape[1]):
        logging.debug("-->SiPM number ={}".format(j))
        rdata.append(sipmrd_[event_number, j])
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

        
if __name__ == '__main__':
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


    print("input path ={}; output path = {}; file name ={}".format(
        PATH_IN,PATH_OUT,FILE))

    print("first event = {} last event = {} number of events = {} ".format(
        FIRST_EVT,LAST_EVT,NEVENTS))

    # open the input file 
    with tables.open_file(PATH_IN+FILE, "r") as h5in: 
        # access the PMT raw data in file 
        pmtrd_ = h5in.root.pmtrd
        sipmrd_ = h5in.root.sipmrd

        #pmtrd_.shape = (nof_events, nof_sensors, wf_length)

        #access the geometry and the sensors metadata info

        geom_t = h5in.root.Detector.DetectorGeometry
        pmt_t = h5in.root.Sensors.DataPMT
        sipm_t = h5in.root.Sensors.DataSiPM
        mctrk_t = h5in.root.MC.MCTracks

        #return a pandas DF for sensors and geometry
        pmtdf = read_data_sensors(pmt_t)
        sipmdf = read_data_sensors(sipm_t)
        geodf = read_data_geom(geom_t)

        #FEE simulation data
        feedf =FEE_param()
        

        # open the output file 
        with tables.open_file(PATH_OUT+FILE, "w",
            filters=tables.Filters(complib="blosc", complevel=9)) as h5out:
 
            # create a group and a table to store MC data
            mcgroup = h5out.create_group(h5out.root, "MC")

            # copy the mctrk table
            mctrk_t.copy(newparent=mcgroup)
            
            # create an extensible array to store the waveforms
            pmtrd = h5out.create_earray(h5out.root, "pmtrd", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, pmtrd_.shape[1], 
                                        int((pmtrd_.shape[2]+1)/FP.time_DAQ)), 
                                    expectedrows=pmtrd_.shape[0])

            sipmrd = h5out.create_earray(h5out.root, "sipmrd", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, sipmrd_.shape[1], sipmrd_.shape[2]), 
                                    expectedrows=sipmrd_.shape[0])

            
            if NEVENTS > pmtrd_.shape[0] and RUN_ALL == False:
                print("""
                Refusing to run: you have requested
                FIRST_EVT = {}
                LAST_EVT  = {}
                Thus you want to run over {} events
                but the size of the DST is {}
                Please change your choice or select RUN_ALL = TRUE
                to run over the whole DST when this happens
                """.format(FIRST_EVT,LAST_EVT,NEVENTS,pmtrd_.shape[0]))
                sys.exit(0)
            elif  NEVENTS > pmtrd_.shape[0] and RUN_ALL == True:
                FIRST_EVT = 0
                LAST_EVT = pmtrd_.shape[0]
            
            for i in range(FIRST_EVT,LAST_EVT):
                print("-->event number ={}".format(i))
                logging.info("-->event number ={}".format(i))
                
                dataPMT = simulate_pmt_response(i,pmtrd_)
                
                pmtrd.append(dataPMT.reshape(1, pmtrd_.shape[1], 
                    int((pmtrd_.shape[2]+1)/FP.time_DAQ)))

                dataSiPM = copy_sipm(i,sipmrd_)
                
                sipmrd.append(dataSiPM.reshape(1, sipmrd_.shape[1], sipmrd_.shape[2]))

            pmtrd.flush()
            sipmrd.flush()

            #store the DF with sensors and geom info
            store_export = pd.HDFStore(PATH_OUT+FILE)
            store_export.append('Sensors/DataPMT', pmtdf, data_columns=pmtdf.columns)
            store_export.append('Sensors/DataSiPM', sipmdf, data_columns=pmtdf.columns)
            store_export.append('Detector/DetectorGeometry', geodf)
            store_export.append('EnergyPlane/FEE', feedf)
            store_export.close()

    print("done!")

    

        
     