"""
DIOMIRA
Leaving there and proceeding for three days toward the east, 
you reach Diomira, a city with sixty silver domes, bronze statues of all the gods, 
streets paved with lead, a crystal theater, 
a golden cock that crows every morning on a tower. 
All these beauties will already be familiar to the visitor, 
who has seen them also in other cities. 
But the special quality of this city for the man who arrives there on a September evening, 
when the days are growing shorter 
and the multicolored lamps are lighted all at once at the doors of the food stalls 
and from a terrace a woman's voice cries ooh!, 
is that he feels envy toward those who now believe they 
have once before lived an evening identical to this and who think they were happy, that time.

DIOMIRA reads the pyTNT Monte Carlo DST and produces maps of PMT and SiPM energy
JJGC August 2016


"""
from __future__ import print_function
from Util import *
from PlotUtil import *

#import SierpeConfig as CP
import numpy as np

import tables as tb
import pandas as pd
import logging

import sys
import getopt

"""
Configuration parameters
"""

logging.basicConfig(level=logging.INFO)


"""
Code
"""

def fill_data_(row):
    data ={}
    pos = row['position']
    data['channel'] =row['channel']
    data['x']= pos[0]
    data['y'] =pos[1]
    return data
    
def read_data_sensors(pfile):
    """
    reads the PMT and SIPM metadata and return data frames
    """
    with tb.open_file(pfile, "r") as h5in:
        #geom_t = h5in.root.Detector.DetectorGeometry
        pmt_t = h5in.root.Sensors.DataPMT
        sipm_t = h5in.root.Sensors.DataSiPM
        PMT={}
        SIPM={}
        i=0
        for row in pmt_t.iterrows():
            PMT[i]=fill_data_(row)
            i+=1
        i=0
        for row in sipm_t.iterrows():
            SIPM[i]=fill_data_(row)
            i+=1
    return pd.DataFrame(PMT),pd.DataFrame(SIPM)

def sensor_energy(pfile,event_number,pmtdf,sipmdf):
    """
    Open the pyTNT, loads the pmtrd_ extensible array with 
    PMT waveforms and computes the energy of the waveforms of each PMT and SiPM.
    """
    with tb.open_file(pfile, "r") as h5in:
        pmtrd_ = h5in.root.pmtrd
        sipmrd_ = h5in.root.sipmrd
        enePMT = np.zeros(len(pmtdf.columns), dtype=np.float64)
        eneSIPM = np.zeros(len(sipmdf.columns), dtype=np.float64)
        for j in range(len(pmtdf.columns)):
            pmtrd = pmtrd_[event_number, j]
            enePMT[j] = sum(pmtrd)
        for j in range(len(sipmdf.columns)):
            sipmrd = sipmrd_[event_number, j]
            eneSIPM[j] = sum(sipmrd)
        
    pmtdf.ix['pes'] =enePMT 
    sipmdf.ix['pes'] =eneSIPM 

def plot_sensor(sdf, rd):
    """
    plots a map of PMTs/SiPMs
    """
    plt.figure(figsize=(10,10))
    ax = plt.subplot(aspect='equal')
    x =[]
    y =[]
    r =[]
    col = []
    for i in sdf.columns:
        sensor = sdf[i]
        x.append(sensor['x'])
        y.append(sensor['y'])
        r.append(rd)
        col.append(sensor['pes'])
    circles(x, y, r, c=col, alpha=0.5, ec='none')
    plt.colorbar()
    plt.xlim(-198,198)  #dimensions of NEW
    plt.ylim(-198,198)
    plt.show()

def main(argv):
    path='/Users/jjgomezcadenas/Documents/Development/NEXT/data/Waveforms/'
    ffile="WF_Tl_0.h5"
    pfile = path+ffile
    event = int(0)
    try:
        opts, args = getopt.getopt(argv, "hp", ["help", "path"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-p", "--path"):
            pfile = arg
        

    print("reading file ={}".format(pfile))
    print("analyzing event ={}".format(event))
    pmtdf,sipmdf =read_data_sensors(pfile)
    sensor_energy(pfile,event,pmtdf,sipmdf)
    plot_sensor(pmtdf,10)
    plot_sensor(sipmdf,2)

def usage():
    print("""usage: python DIOMIRA.py --help (-h) --path (-p) 
        --help (-h): this text
        --path (-p): full path to your pyTNT MC file:
        
        """)
    
    
if __name__ == '__main__':
    main(sys.argv[1:])
    
    
    

    

        
     