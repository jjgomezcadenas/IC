"""
DOROTEA
The cities and desire, chaper I
JJGC August 2016

What DOROTEA does:
1) Reads the pyTNT (pyTablesNexT) DST 
2) Analizes MCTrack data

"""
from __future__ import print_function
from Util import *
from PlotUtil import *

#import SierpeConfig as CP
import numpy as np

import tables
import logging

"""
Configuration parameters
"""

logging.basicConfig(level=logging.INFO)

DATASET = ["WF_Kr_0.h5",
           "WF_Na_0.h5",
           "WF_Tl_0.h5"]

PATH_IN="/Users/jjgomezcadenas/Documents/Development/NEXT/data/Waveforms/"

FRST_EVT = 0
LST_EVT = 3

"""
Code
"""


if __name__ == '__main__':
    
    print("""DOROTEA
        This program reads the pyTNT DST and analyzes MCTrack data 
        """)
    
    file = DATASET[0]

    print("input path ={} file name ={}".format(
        PATH_IN,file))

    print("number of events = {} ".format(
        LST_EVT, FRST_EVT))

    # open the input file 
    with tables.open_file(PATH_IN+file, "r") as h5in:
        table = h5in.root.MC.MCTracks

        for row in table.iterrows():
            if row['hit_indx'] == 0

       
pressure = [x['pressure'] for x in table.iterrows() if x['TDCcount'] > 3 and 20 <= x['pressure'] < 50]
>>> pressure
            for i in range(RUNP['frst_evt'],RUNP['lst_evt']):
                print("-->event number ={}".format(i))
                logging.info("-->event number ={}".format(i))
                
                rdata = simulate_pmt_response(i,pmtrd_)
                # fill raw data in
                pmtrd25.append(rdata.reshape(1, len(CFGP['PMTS']), CFGP['LEN_WVF25']))

                if CPLT['INTER'] == True:
                    wait()

            pmtrd25.flush()

    print("done!")

    

        
     