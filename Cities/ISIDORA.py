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
def accumulator_coefficients(NPMT,CA):
    """
    Compute the accumulator coefficients
    """
    import FEParam as FP
    import FEE2 as FE
    
    coef_acc =np.zeros(NPMT, dtype=np.double)
    for j in range(NPMT):
        
        fee = FE.FEE(C=CA[j],R= FP.R, f=FP.freq_LPF, RG=FP.V_GAIN)
        signal_inv_daq = fee.InverseSignalDAQ(signal_t)  #inverse function
        coef_acc[j] = signal_inv_daq[10] #any index is valid, function is flat
        
    return coef_acc

def DBLR(pmtrd_, event_number, CA, pmt_range=0, MAU_LEN=250,
         thr1 = 20., thr2=0., thr3 = 5., plot=False , log='INFO'):
    """
    Peform Base line Restoration
    CA is an array with the coefficients of the accumulator
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
        
        pmtwf, ene_pmt[j] = BLR(pmtrd, CA[j], MAU_LEN, thr1, thr2, thr3, plot,log)
        PMTWF.append(pmtwf)
       
    return ene_pmt, np.array(PMTWF)
                                        
def DBLR(pmtrd_,signal_t, event_number=0, 
         CA=FP.C12, pmt_range=0, 
         mau_len=FP.MAU_WindowSize,
         thr1 = 20., thr2=0., thr3 = 5., 
         plot='True', log='DEBUG'):
    """
    Peform Base line Restoration
    CA is an array with the values of the capacitances for the PMTs
    Threshold 1 is used to decide when raw signal raises up above trigger line
    Threshold 2 is used to decide when reconstructed signal is above trigger line
    Threshold 3 is used to compare Raw and Rec signal
    """
    import FEParam as FP
    import FEE2 as FE
    
    len_WF = pmtrd_.shape[2]
    NPMT = pmtrd_.shape[1]
    ene_pmt =np.zeros(NPMT, dtype=np.double)
    coef_pmt =np.zeros(NPMT, dtype=np.double)
    
    PMTWF =[]
    
    pmts = pmt_range
    if pmt_range == 0:
        pmts = range(NPMT)
    for j in pmts:

        pmtrd = pmtrd_[event_number, j] #waveform for event event_number, PMT j
        
        #Deconvolution
        fee = FE.FEE(C=CA[j],R= FP.R, f=FP.freq_LPF, RG=FP.V_GAIN)
        signal_inv_daq = fee.InverseSignalDAQ(signal_t)  #inverse function
        coef = signal_inv_daq[10]  #accumulator coefficient
       
        signal_blr, eadc = BLR(pmtrd, coef, j, MAU_LEN, thr1, thr2, thr3, plot,log)
    
        ene_pmt[j] = eadc
        coef_pmt[j] = coef
        PMTWF.append(signal_blr)
       
    return np.array(PMTWF), ene_pmt, coef_pmt
        
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

    print("input path ={}; file_in ={} ".format(
        PATH_IN,FILE_IN))

    print("first event = {} last event = {} nof events requested = {} ".format(
        FIRST_EVT,LAST_EVT,NEVENTS))

    # open the input file in mode append 
    with tables.open_file("{}/{}".format(PATH_IN,FILE_IN), "a") as h5in: 
        # access the PMT raw data in file 
        pmtrd_ = h5in.root.pmtrd
        sipmrd_ = h5in.root.sipmrd

        #pmtrd_.shape = (nof_events, nof_sensors, wf_length)

        #define a time vector needed for deconvolution
        signal_t = np.arange(0.0, pmtrd_.shape[2]*1., 1., dtype=np.double)

        
        NPMT = pmtrd_.shape[1]
        NSIPM = sipmrd_.shape[1]
        PMTWL = pmtrd_.shape[2] 
        SIPMWL = sipmrd_.shape[2]
        NEVENTS_DST = pmtrd_.shape[0]

        print("nof PMTs = {} nof  SiPMs = {} nof events in input DST = {} ".format(
        NPMT,NSIPM,NEVENTS_DST))

        print("lof SiPM WF = {} lof PMT WF (MC) = {} ".format(
        PMTWL,SIPMWL))

        wait()
            
        # create an extensible array to store the CWF waveforms
        pmtcwf = h5in.create_earray(h5in.root.RD, "pmtcwf", 
                                    atom=tables.FloatAtom(), 
                                    shape=(0, NPMT, PMTWL), 
                                    expectedrows=NEVENTS_DST)
            

        #create an extensible array to store the energy in adc counts of CWF 
        ecwf = h5out.create_earray(h5out.root.RD, "ecwf", 
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
            pmtCWF, eneCWF, coefBLR = DBLR(pmtrd_,signal_t, 
                                        event_number =i, CA =C12, 
                                        pmt_range=0, mau_len=MAU_LEN,
                                        thr1 = THR1, thr2=THR2, thr3 =THR3, 
                                        plot=PLOT, log=LOG):

            
                
            #append to PMT EARRAY
            pmtcwf.append(pmtCWF.reshape(1, NPMT, PMTWL))
            #append to ecwf EARRAY
            ecwf.append(eneCWF.reshape(1, NPMT))
                
                
        pmtcwf.flush()
        ecwf.flush()
        

    print("done!")

    

        
     