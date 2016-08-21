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
import Configure as CF

import numpy as np

import FEParam as FP
import SPE as SP
import FEE2 as FE
import BLR

import tables
import pandas as pd
import logging


"""
Code
"""


def accumulator_coefficients(pmtrd_,CA):
    """
    Compute the accumulator coefficients for DBLR
    It computes the inverse function of the HPF and takes
    the accumulator as the value of the function anywhere
    but the first bin (the inverse is a step function with
    constant value equal to the accumulator)
    CA are the values of the capacitances defining the filter
    (1/(2CR)) for each PMT
    """
    len_WF = pmtrd_.shape[2]
    NPMT = pmtrd_.shape[1]
    
    coef_acc =np.zeros(NPMT, dtype=np.double)

    signal_t = np.arange(0.0, len_WF*1., 1., dtype=np.double)

    for j in range(NPMT):
        
        fee = FE.FEE(C=CA[j],R= FP.R, f=FP.freq_LPF, RG=FP.V_GAIN)
        signal_inv_daq = fee.InverseSignalDAQ(signal_t)  #inverse function
        coef_acc[j] = signal_inv_daq[10] #any index is valid, function is flat
        
    return coef_acc

def DBLR(pmtrd_, event_number, coeff_acc, mau_len=250,
         thr1 = 20., thr2=0., thr3 = 5., log='INFO', pmt_range=0):
    """
    Peform Base line Restoration
    coeff_acc is an array with the coefficients of the accumulator
    Threshold 1 is used to decide when raw signal raises up above trigger line
    Threshold 2 is used to decide when reconstructed signal is above trigger line
    Threshold 3 is used to compare Raw and Rec signal
    """
    
    len_WF = pmtrd_.shape[2]
    NPMT = pmtrd_.shape[1]
    ene_pmt =np.zeros(NPMT, dtype=np.int64)
   
    PMTWF =[]
    
    pmts = pmt_range
    if pmt_range == 0:
        pmts = range(NPMT)

    for j in pmts:
        pmtrd = pmtrd_[event_number, j] #waveform for event event_number, PMT j
        
        pmtwf, ene_pmt[j] = BLR.BLR(pmtrd, coeff_acc[j], mau_len, thr1, thr2, thr3, log)
        PMTWF.append(pmtwf)
       
    return ene_pmt, np.array(PMTWF)
                                        

if __name__ == '__main__':
    INFO, CFP = CF.configure(sys.argv[0],sys.argv[1:])

    if INFO:
        print(isidora)

    wait()
    
    print("""
        ISIDORA:
        1. Reads an Nh5 file produced by DIOMIRA, which stores the
            raw waveforms (RWF) for the PMTs and SiPMs waveforms, as well as
            data on geometry, sensors and MC. The RDWF of the PMTs
            show negative swing due to the HPF of the EP FEE electronics

        2. Performs DBLR on the PMT RWF and produces corrected waveforms (CWF).

        3. Adds the CWF to the DST 

        4. Computes the energy of the PMTs per each event and writes to DST

        """)

    
    PATH_IN=CFP['PATH_IN'] 
    FILE_IN=CFP['FILE_IN']
    FIRST_EVT=CFP['FIRST_EVT']
    LAST_EVT=CFP['LAST_EVT']
    RUN_ALL=CFP['RUN_ALL']
    CA=farray_from_string(CFP['CA'])*nF 
    MAU_LEN=CFP['MAU_LEN']
    THR1=CFP['THR1'] 
    THR2=CFP['THR2'] 
    THR3=CFP['THR3']

    NEVENTS = LAST_EVT -  FIRST_EVT


    print("input path ={}; file_in ={} ".format(
        PATH_IN,FILE_IN))

    print("first event = {} last event = {} nof events requested = {} ".format(
        FIRST_EVT,LAST_EVT,NEVENTS))

    print("MAU length = {} THR1 = {} THR2 = {} THR3 = {} ".format(
        MAU_LEN,THR1,THR2,THR3))
    print("CA (nf) = {}  ".format(CA/nF))
    

    # open the input file in mode append 
    with tables.open_file("{}/{}".format(PATH_IN,FILE_IN), "a") as h5in: 
        # access the PMT raw data in file 
        pmtrd_ = h5in.root.RD.pmtrd

        #pmtrd_.shape = (nof_events, nof_sensors, wf_length)    
        
        NPMT = pmtrd_.shape[1]
        PMTWL = pmtrd_.shape[2] 
        NEVENTS_DST = pmtrd_.shape[0]

        print("nof PMTs = {} nof events in input DST = {} ".format(
        NPMT,NEVENTS_DST))

        print("lof PMT WF (MC) = {} ".format(
        PMTWL))

        wait()
            
        # create an extensible array to store the CWF waveforms
        pmtcwf = h5in.create_earray(h5in.root.RD, "pmtcwf", 
                                    atom=tables.FloatAtom(), 
                                    shape=(0, NPMT, PMTWL), 
                                    expectedrows=NEVENTS_DST)
            

        #create an extensible array to store the energy in adc counts of CWF 
        ecwf = h5in.create_earray(h5in.root.RD, "ecwf", 
                                    atom=tables.FloatAtom(), 
                                    shape=(0, NPMT), 
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

        #create a 2d array to store the  energy of the CWF
        ene_cwf = np.zeros((NEVENTS,NPMT))

            
        for i in range(FIRST_EVT,LAST_EVT):
            print("-->event number ={}".format(i))
            logging.info("-->event number ={}".format(i))

            #Perform DBLR
            coeff_acc = accumulator_coefficients(pmtrd_,CA)

            eneCWF, pmtCWF = DBLR(pmtrd_, i, coeff_acc, mau_len=MAU_LEN,
                                  thr1 = THR1, thr2=THR2, thr3 = THR3, 
                                  log='INFO', pmt_range=0)
            
                
            #append to PMT EARRAY
            pmtcwf.append(pmtCWF.reshape(1, NPMT, PMTWL))
            #append to ecwf EARRAY
            ecwf.append(eneCWF.reshape(1, NPMT))
                
                
        pmtcwf.flush()
        ecwf.flush()
        

    print("Leaving Isidora. Safe travels!")

    

        
     