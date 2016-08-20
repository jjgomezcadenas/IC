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

def BLR(signal_daq, coef, mau_len, thr1, thr2, thr3, log):

    """
    Deconvolution offline of the DAQ signal using a MAU
    moving window-average filter of a vector data. See notebook 
    y(n) = (1/WindowSize)(x(n) + x(n-1) + ... + x(n-windowSize))
    in a filter operation filter(b,a,x):
    b = (1/WindowSize)*ones(WindowSize) = (1/WS)*[1,1,1,...]: numerator
    a = 1 : denominator
    y = filter(b,a,x)
    y[0] = b[0]*x[0] = (1/WS) * x[0]
    y[1] = (1/WS) * (x[0] + x[1])
    y[WS-1] = mean(x[0:WS])
    y[WS] = mean(x[1:WS+1])
    and so on
    """
    
    logl ='logging.'+log
    logging.basicConfig(level=eval(logl))
    

    len_signal_daq = len(signal_daq)
    MAU = np.zeros(len_signal_daq, dtype=np.double)
    acum = np.zeros(len_signal_daq, dtype=np.double)
    signal_r = np.zeros(len_signal_daq, dtype=np.double)
    pulse_f = np.zeros(len(signal_daq), dtype=np.double)
    pulse_ff = np.zeros(len(signal_daq), dtype=np.double)
    pulse_t = np.zeros(len(signal_daq), dtype=np.double)
    pulse_w = np.zeros(len(signal_daq), dtype=np.double)
    
    signal_i = np.copy(signal_daq) #uses to update MAU while procesing signal

    nm = mau_len
    B_MAU = (1./nm)*np.ones(nm)

#   MAU averages the signal in the initial tranch 
#    allows to compute the baseline of the signal  
    
    MAU[0:nm] = SGN.lfilter(B_MAU,1, signal_daq[0:nm])
    acum[nm] =  MAU[nm]
    BASELINE = MAU[nm-1]

    logging.debug("""-->BLR: 
                     PMT number = {}
                     MAU_LEN={}
                     thr1 = {}, thr2 = {}, thr3 = {} =""".format(
                     pmt, mau_len, thr1, thr2, thr3))
    logging.debug("n = {}, acum[n] = {} BASELINE ={}".format(nm, acum[nm],BASELINE))

#----------

# While MAU inits BLR is switched off, thus signal_r = signal_daq 

    signal_r[0:nm] = signal_daq[0:nm] 
    pulse_on=0
    wait_over=0
    offset = 0
    
    # MAU has computed the offset using nm samples
    # now loop until the end of DAQ window

    logging.debug("nm = {}".format(nm))
    
    for k in range(nm,len_signal_daq): 

        trigger_line = MAU[k-1] + thr1
        pulse_t[k] = trigger_line
        pulse_f[k] = pulse_on
        pulse_w[k] = wait_over 
        pulse_ff[k] = signal_daq[k] - signal_r[k]
        
        # condition: raw signal raises above trigger line and 
        # we are not in the tail
        # (wait_over == 0)
        if signal_daq[k] > trigger_line and wait_over == 0:

            # if the pulse just started pulse_on = 0.
            # In this case compute the offset as value 
            #of the MAU before pulse starts (at k-1)

            if pulse_on == 0: # pulse just started
                #offset computed as the value of MAU before pulse starts
                offset = MAU[k-1]  
                pulse_on = 1 
                
            #Pulse is on: Freeze the MAU
            MAU[k] = MAU[k-1]  
            signal_i[k] =MAU[k-1]  #signal_i follows the MAU
            
            #update recovered signal, correcting by offset
            acum[k] = acum[k-1] + signal_daq[k] - offset;
            signal_r[k] = signal_daq[k] + coef*acum[k] 
                  
            
        else:  #no signal or raw signal has dropped below threshold
                      
        # but raw signal can be negative for a while and still contribute to the
        # reconstructed signal.

            if pulse_on == 1: #reconstructed signal still on
                # switch the pulse off only when recovered signal 
                #drops below threshold
                #slide the MAU, still frozen. 
                # keep recovering signal
                
                MAU[k] = MAU[k-1] 
                signal_i[k] =MAU[k-1]
                acum[k] = acum[k-1] + signal_daq[k] - offset;
                signal_r[k] = signal_daq[k] + coef*acum[k] 
                
                
                #if the recovered signal drops before trigger line 
                #rec pulse is over!
                if signal_r[k] < trigger_line + thr2:
                    wait_over = 1  #start tail compensation
                    pulse_on = 0   #recovered pulse is over
                      

            else:  #recovered signal has droped below trigger line
            #need to compensate the tail to avoid drifting due to erros in 
            #baseline calculatoin

                if wait_over == 1: #compensating pulse
                    # recovered signal and raw signal 
                    #must be equal within a threshold
                    # otherwise keep compensating pluse

                        
                    if signal_daq[k-1] < signal_r[k-1] - thr3:
                        # raw signal still below recovered signal 
                        # keep compensating pulse
                        # is the recovered signal near offset?
                        upper = offset + (thr3 + thr2)
                        lower = offset - (thr3 + thr2)
                        
                        if signal_r[k-1] > lower and signal_r[k-1] < upper:
                            # we are near offset, activate MAU. 
                            
                            signal_i[k] = signal_r[k-1]
                            MAU[k] = np.sum(signal_i[k-nm:k])*1./nm
                            
                                      
                        else: 
                            # rec signal not near offset MAU frozen  
                            MAU[k] = MAU[k-1]
                            signal_i[k] = MAU[k-1]
                            

                        # keep adding recovered signal  
                        acum[k] = acum[k-1] + signal_daq[k] - MAU[k]
                        signal_r[k] = signal_daq[k] + coef*acum[k]
                        
                    else:  # raw signal above recovered signal: we are done 
                        
                        wait_over = 0
                        acum[k] = MAU[k-1]
                        signal_r[k] = signal_daq[k]
                        signal_i[k] = signal_r[k]
                        MAU[k] = np.sum(signal_i[k-nm:k])*1./nm
                        
                            
                else: #signal still not found
                    
                    #update MAU and signals
                    MAU[k] = np.sum(signal_i[k-nm:k]*1.)/nm   
                    acum[k] = MAU[k-1]
                    signal_r[k] = signal_daq[k]
                    signal_i[k] = signal_r[k]  
                                                                                                       
    energy = np.dot(pulse_f,(signal_r-BASELINE)) 
                       
    return  signal_r-BASELINE, energy

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
        
        pmtwf, ene_pmt[j] = BLR(pmtrd, coeff_acc[j], mau_len, thr1, thr2, thr3, log)
        PMTWF.append(pmtwf)
       
    return ene_pmt, np.array(PMTWF)
                                        

def usage():
    """
    Usage of program
    """
    print("""
        Usage: python (run) ISIDORA [args]
        where args are:
         -h (--help) : this text
         -i (--info) : print a text describing the invisible city of ISIDORA
         -d (--debug) : can be set to 'DEBUG','INFO','WARNING','ERROR'
         -c (--cfile) : full path to a configuration file
         
         example of configuration file 

        # comment line  
        Names of parameters (comma separated)
        Values of parameters (comma separated)
        
        The parameters for ISIDORA are:

        PATH_IN = path to DST file
        FILE_IN = name of DST file
        FIRST_EVT,LAST_EVT,RUN_ALL,

        RUN_ALL is used to decide whether to run all the events in the file
        in case that the total number of events requested (LAST_EVT-FIRST_EVT) 
        exceeds the number of events in the DST file. If RUN_ALL is set to 1 (True), 
        the script will run over all elements in the DST, 
        otherwise it will exit with a warning.

        CA = Nominal (measured) values of the capacitors defining the filter
        MAU_LEN = length of MAU
        THR1, THR2, THR3, are trhesholds for the DBLR algorithm 
        Threshold 1 is used to decide when raw signal raises up above trigger line
        Threshold 2 is used to decide when reconstructed signal is above trigger line
        Threshold 3 is used to compare Raw and Rec signal
        

        """)
def configure(argv):
    """
    reads arguments from the command line and configures job
    """
    
    global DEBUG, PATH_IN, FILE_IN,INFO 
    global  FIRST_EVT, LAST_EVT, NEVENTS, RUN_ALL
    global  CA, MAU_LEN, THR1, THR2, THR3
    
    DEBUG='INFO'
    cfile =''
    INFO = True
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

    lg = 'logging.'+DEBUG
    logging.basicConfig(level=eval(lg))
 
    if cfile == '':
        print("Path to configuration file not given")
        usage()
        sys.exit()

    CFP =pd.read_csv(cfile,comment="#")
    print("""
        Configuration parameters \n 
        {}
        """.format(CFP))
    
    PATH_IN=CFP['PATH_IN'][0] 
    FILE_IN=CFP['FILE_IN'][0]
    FIRST_EVT=CFP['FIRST_EVT'][0]
    LAST_EVT=CFP['LAST_EVT'][0]
    RUN_ALL=CFP['RUN_ALL'][0]
    CA=CFP['CA'][0]
    MAU_LEN=CFP['MAU_LEN'][0]
    THR1=CFP['THR1'][0] 
    THR2=CFP['THR2'][0] 
    THR3=CFP['THR3'][0]

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

    

        
     