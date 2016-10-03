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
from LogConfig import *
from Configure import *

from BLR import accumulator_coefficients,DBLR
import FEParam as FP
import tables
from time import time 
#import pandas as pd


"""
Code
"""
def ISIDORA(argv):
    DEBUG_LEVEL, INFO, CYTHON, CFP = configure(argv[0],argv[1:])
    
    if INFO:
        
    
        print("""
        ISIDORA:
        1. Reads an Nh5 file produced by DIOMIRA, which stores the
            raw waveforms (RWF) for the PMTs and SiPMs waveforms, as well as
            data on geometry, sensors and MC. The RDWF of the PMTs
            show negative swing due to the HPF of the EP FEE electronics

        2. Performs DBLR on the PMT RWF and produces corrected waveforms (CWF).

        3. Adds the CWF and ancilliary info to the DST 

        4. Computes the energy of the PMTs per each event and writes to DST

        """)

    
    PATH_IN=CFP['PATH_IN'] 
    FILE_IN=CFP['FILE_IN']
    FIRST_EVT=CFP['FIRST_EVT']
    LAST_EVT=CFP['LAST_EVT']
    RUN_ALL=CFP['RUN_ALL']
    COEF=CFP['COEF']
    CA=farray_from_string(CFP['CA'])*nF 
    AC=farray_from_string(CFP['AC']) 
    MAU_LEN=CFP['MAU_LEN']
    NSIGMA1=CFP['NSIGMA1'] 
    NSIGMA2=CFP['NSIGMA2'] 
    NSIGMA3=CFP['NSIGMA3'] 

    NEVENTS = LAST_EVT -  FIRST_EVT

    print('Debug level = {}'.format(DEBUG_LEVEL))

    print("input path ={}; file_in ={} ".format(
        PATH_IN,FILE_IN))

    print("first event = {} last event = {} nof events requested = {} ".format(
        FIRST_EVT,LAST_EVT,NEVENTS))

    print("MAU length = {} n_sigma1 = {} n_sigma2 = {} ".format(
        MAU_LEN,NSIGMA1,NSIGMA2))
    print("CA  = {} nF ".format(CA/nF))
    print("Accumulator Coefficients = {}  ".format(AC))
    

    # open the input file in mode append 
    with tables.open_file("{}/{}".format(PATH_IN,FILE_IN), "a") as h5in: 
        # access the PMT raw data in file 
        pmtrd_ = h5in.root.RD.pmtrwf     #PMT raw data must exist

        #pmtrd_.shape = (nof_events, nof_sensors, wf_length)    
        
        NPMT = pmtrd_.shape[1]
        PMTWL = pmtrd_.shape[2] 
        NEVENTS_DST = pmtrd_.shape[0]

        print("nof PMTs = {} nof events in input DST = {} ".format(
        NPMT,NEVENTS_DST))

        print("lof PMT WF (MC) = {} ".format(
        PMTWL))

        #wait()
            
        # create an extensible array to store the CWF waveforms
        # if it exists remove and create again
        pmtcwf =0
        try:
            pmtcwf = h5in.root.RD.pmtcwf
            h5in.remove_node("/RD","pmtcwf")
            pmtcwf = h5in.create_earray(h5in.root.RD, "pmtcwf", 
                                    atom=tables.Float32Atom(), 
                                    shape=(0, NPMT, PMTWL), 
                                    expectedrows=NEVENTS_DST)
        except tables.exceptions.NodeError:
            pmtcwf = h5in.create_earray(h5in.root.RD, "pmtcwf", 
                                    atom=tables.Float32Atom(), 
                                    shape=(0, NPMT, PMTWL), 
                                    expectedrows=NEVENTS_DST)

        # create a group to store BLR configuration 
        # mau, acummulator, pulse_on and wait_over stored for pmt 0
        # baseline stored for all PMTs.

        rgroup = 0
        try:
            rgroup = h5in.root.BLR
            
        except tables.exceptions.NodeError:
            rgroup = h5in.create_group(h5in.root, "BLR")
        
        mau = 0
        acum = 0
        pulse_on = 0
        wait_over = 0
        baseline = 0

        try:
            mau = h5in.root.BLR.mau
            h5in.remove_node("/BLR","mau")
            mau = h5in.create_earray(h5in.root.BLR, "mau", 
                                    atom=tables.Float32Atom(), 
                                    shape=(0, PMTWL), 
                                    expectedrows=NEVENTS_DST)
        except tables.exceptions.NodeError:
            mau = h5in.create_earray(h5in.root.BLR, "mau", 
                                    atom=tables.Float32Atom(), 
                                    shape=(0, PMTWL), 
                                    expectedrows=NEVENTS_DST)
        try:
            pulse_on = h5in.root.BLR.pulse_on
            h5in.remove_node("/BLR","pulse_on")
            pulse_on = h5in.create_earray(h5in.root.BLR, "pulse_on", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, PMTWL), 
                                    expectedrows=NEVENTS_DST)
        except tables.exceptions.NodeError:
            pulse_on = h5in.create_earray(h5in.root.BLR, "pulse_on", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, PMTWL), 
                                    expectedrows=NEVENTS_DST)
        try:
            wait_over = h5in.root.BLR.wait_over
            h5in.remove_node("/BLR","wait_over")
            wait_over = h5in.create_earray(h5in.root.BLR, "wait_over", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, PMTWL), 
                                    expectedrows=NEVENTS_DST)
        except tables.exceptions.NodeError:
            wait_over = h5in.create_earray(h5in.root.BLR, "wait_over", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, PMTWL), 
                                    expectedrows=NEVENTS_DST)
        try:
            acum  = h5in.root.BLR.acum 
            h5in.remove_node("/BLR","acum")
            acum  = h5in.create_earray(h5in.root.BLR, "acum", 
                                    atom=tables.Float32Atom(), 
                                    shape=(0, PMTWL), 
                                    expectedrows=NEVENTS_DST)
        except tables.exceptions.NodeError:
            acum  = h5in.create_earray(h5in.root.BLR, "acum", 
                                    atom=tables.Float32Atom(), 
                                    shape=(0, PMTWL), 
                                    expectedrows=NEVENTS_DST)
        try:
            baseline  = h5in.root.BLR.baseline 
            h5in.remove_node("/BLR","baseline")
            baseline  = h5in.create_earray(h5in.root.BLR, "baseline", 
                                    atom=tables.Float32Atom(), 
                                    shape=(0, NPMT), 
                                    expectedrows=NEVENTS_DST)
        except tables.exceptions.NodeError:
            baseline  = h5in.create_earray(h5in.root.BLR, "baseline", 
                                    atom=tables.Float32Atom(), 
                                    shape=(0, NPMT), 
                                    expectedrows=NEVENTS_DST)

            
        #LOOP
        first_evt, last_evt = define_event_loop(FIRST_EVT,LAST_EVT,NEVENTS,NEVENTS_DST,RUN_ALL)
            
        t0 = time()
        for i in range(first_evt,last_evt):
            
            logger.info("-->event number ={}".format(i))

            #DBLR
            coeff_acc = AC
            if COEF == 0:
                coeff_acc = accumulator_coefficients(pmtrd_,CA)


            BLRS = DBLR(pmtrd_, i, coeff_acc, mau_len=MAU_LEN,
                        thr1 = NSIGMA1*FP.NOISE_ADC, thr2=NSIGMA2*FP.NOISE_ADC, 
                        thr3 = NSIGMA3*FP.NOISE_ADC, log=DEBUG_LEVEL)

            
            pmtCWF = []
            for blr in BLRS:
                pmtCWF.append(blr.signal_r)
                
            #append to pmtcwd
            pmtcwf.append(np.array(pmtCWF).reshape(1, NPMT, PMTWL))

            # append BLR variables
            
            
            BASELINE = []
            for blr in BLRS:
                BASELINE.append(blr.BASELINE)
            baseline.append(np.array(BASELINE).reshape(1, NPMT))

            mau_pmt0 = BLRS[0].MAU
            mau.append(mau_pmt0.reshape(1, PMTWL))

            pulse_on0 = BLRS[0].pulse_on
            pulse_on.append(pulse_on0.reshape(1, PMTWL))

            wait_over0 = BLRS[0].wait_over
            wait_over.append(wait_over0.reshape(1, PMTWL))

            acum0  = BLRS[0].acum
            acum.append(acum0.reshape(1, PMTWL))
                
        t1 = time()
        pmtcwf.flush()
        mau.flush()
        pulse_on.flush()
        wait_over.flush()
        acum.flush()
        baseline.flush()
        
        print("ISIDORA has run over {} events in {} seconds".format(i, t1-t0))
    print("Leaving ISIDORA. Safe travels!")
                                  

if __name__ == '__main__':
    ISIDORA(sys.argv)