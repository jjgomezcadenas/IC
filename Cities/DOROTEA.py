"""
DOROTEA
JJGC Agusut 2016

DOROTEA is an alternative version of DIOMIRA. It does the same than DIOMIRA
but adds the RWFs to the existing file with the MCRD, rather than creating a new file 
1) Reads a MCRD file containing MC waveforms for the 12 PMTs of the EP.
   Each waveform contains number of PEs in bins of 1 ns.
2) Convolves the PE waveform with the response of the FEE electronics.
3) Decimates the waveform, simulating the effect of the DAQ sampling (25 ns bins)
4) Writes a RWF file with the new data and adds the FEE simulation parameters as metadata
"""

from __future__ import print_function
from SensorsResponse import *
from PlotUtil import *
import Configure as CF
from Nh5 import *
from cities import dorotea
from LogConfig import *

"""
Code
"""

        
if __name__ == '__main__':
    INFO, CFP = CF.configure(sys.argv[0],sys.argv[1:])

    if INFO:
        print(dorotea)

    wait()
    
    print("""
        DOROTEA:
        1. Reads an MCRD file produced by NEXUS, which stores TWF1ns waveforms 
            for the PMTs and SiPMs as well as data on geometry, sensors and MC.

        2. Simulates the response of the energy plane in the PMTs  and
            outputs PMT TWF25ns e.g., waveforms in bins of 25 ns (and in adc counts)
            which correspond to the output of the EP FEE.

        3. Adds the RWF to the file. It also adds vectors
            holding the waveform energy in pes (per PMT/SiPM)

        4. Add a table describing the FEE parameters used for simulation


        """)
    FP.print_FEE()
    wait()

    PATH_IN =CFP['PATH_IN']
    FILE_IN =CFP['FILE_IN']
    FIRST_EVT =CFP['FIRST_EVT']
    LAST_EVT =CFP['LAST_EVT']
    RUN_ALL =CFP['RUN_ALL']
    NEVENTS = LAST_EVT - FIRST_EVT

    print("input path ={};  file_in ={}".format(
        PATH_IN,FILE_IN))

    print("first event = {} last event = {} nof events requested = {} ".format(
        FIRST_EVT,LAST_EVT,NEVENTS))

    # open the input file 
    with tables.open_file("{}/{}".format(PATH_IN,FILE_IN), "a") as h5in: 
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
        logger.info("nof PMTs = {} nof  SiPMs = {} nof events in input DST = {} ".format(
        NPMT,NSIPM,NEVENTS_DST))

        print("lof SiPM WF = {} lof PMT WF (MC) = {} lof PMT WF (FEE) = {}".format(
        PMTWL,SIPMWL,PMTWL_FEE))

        logger.info("lof SiPM WF = {} lof PMT WF (MC) = {} lof PMT WF (FEE) = {}".format(
        PMTWL,SIPMWL,PMTWL_FEE))

        wait()

        
        logger.info("Check it table FEE exits. Otherwise create it")
        
        fee_table = 0
        mcgroup = h5in.root.MC
        try:
            fee_table = h5in.root.MC.FEE
            
        except tables.exceptions.NodeError:
            print("creating a table for FEE data")
            logger.info("Check it table FEE exits. Otherwise create it")
            fee_table = h5in.create_table(mcgroup, "FEE", FEE,
                          "EP-FEE parameters",
                           tables.Filters(0))

        logger.info("Fill FEE table")
        FEE_param_table(fee_table)


        logger.info("Check it group RD exits. Otherwise create it")
        rgroup = 0
        try:
            rgroup = h5in.root.RD
        except tables.exceptions.NodeError:
            rgroup = h5in.create_group(h5in.root, "RD")
            

        
        logger.info("""create extensible arrays to store the RWF waveforms
                        delete any existing (previous) array""")
        pmtrwf = 0
        try:
            pmtrwf = h5in.root.RD.pmtrwf
            h5in.remove_node("/RD","pmtrwf")
            pmtrwf = h5in.create_earray(h5in.root.RD, "pmtrwf", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NPMT, PMTWL_FEE), 
                                    expectedrows=NEVENTS_DST)

        except tables.exceptions.NodeError:
            pmtrwf = h5in.create_earray(h5in.root.RD, "pmtrwf", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NPMT, PMTWL_FEE), 
                                    expectedrows=NEVENTS_DST)
            
        sipmrwf = 0
        try:
            sipmrwf = h5in.root.RD.sipmrwf
            h5in.remove_node("/RD","sipmrwf")
            sipmrwf = h5in.create_earray(h5in.root.RD, "sipmrwf", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NSIPM, SIPMWL), 
                                    expectedrows=NEVENTS_DST)
        except tables.exceptions.NodeError:
            sipmrwf = h5in.create_earray(h5in.root.RD, "sipmrwf", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NSIPM, SIPMWL), 
                                    expectedrows=NEVENTS_DST)

        epmt = 0
        try:
            epmt = h5in.root.RD.epmt
            h5in.remove_node("/RD","epmt")
            epmt = h5in.create_earray(h5in.root.RD, "epmt", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NPMT), 
                                    expectedrows=NEVENTS_DST)
        except tables.exceptions.NodeError:
            pmt = h5in.create_earray(h5in.root.RD, "epmt", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NPMT), 
                                    expectedrows=NEVENTS_DST)

        logger.info("""create extensible arrays to store the true
            energy (pes) of PMTs and SiPMs. Delete existing previous arrays""") 

        esipm = 0
        try:
            esipm = h5in.root.RD.esipm
            h5in.remove_node("/RD","esipm")
            esipm = h5in.create_earray(h5in.root.RD, "esipm", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, NSIPM), 
                                    expectedrows=NEVENTS_DST)
        except tables.exceptions.NodeError:
            esipm = h5in.create_earray(h5in.root.RD, "esipm", 
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


        for i in range(FIRST_EVT,LAST_EVT):
            print("-->event number ={}".format(i))
            logger.info("-->event number ={}".format(i))

            #simulate PMT response and return an array with new WF
            dataPMT = simulate_pmt_response(i,pmtrd_)
                
            #append to PMT EARRAY
            pmtrwf.append(dataPMT.reshape(1, NPMT, PMTWL_FEE))
                   
            #simulate SiPM response and return an array with new WF
            dataSiPM = simulate_sipm_response(i,sipmrd_)
                
            #append to SiPM EARRAY
            sipmrwf.append(dataSiPM.reshape(1, NSIPM, SIPMWL))

            #fill ene_pmt vector
            enePMT = energy_pes(i, pmtrd_)
            #append to epmt EARRAY
            epmt.append(enePMT.reshape(1, NPMT))

            #fill ene_sipm vector
            eneSIPM = energy_pes(i, sipmrd_)
            esipm.append(eneSIPM.reshape(1, NSIPM))

            pmtrwf.flush()
            sipmrwf.flush()
            epmt.flush()
            esipm.flush()


    print("Leaving Dorothea. Safe travel!")

    

        
     