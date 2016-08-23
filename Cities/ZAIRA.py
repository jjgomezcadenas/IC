"""
ZAIRA 
JMB August 2016

What ZAIRA does:
1) Read a file containing MC waveforms (25 ns)
2) Energy histogram
3) Plot real&measured energy
4) Correlation plots between coordinates and energy
5) Add new tables (?)
"""

from __future__ import print_function

from Util import *

import matplotlib.pyplot as plt
import pandas as pd
import sys
import tables as tb
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import getopt

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

def energyTracks(mctracks):
    """
    Sum of energy of all tracks of each event
    """
    track_energy = {}
    for row in mctracks.iterrows():
        evt = row['event_indx']
        mctrk = row['mctrk_indx']
        if evt not in track_energy:
            track_energy[evt] = {}
            track_energy[evt][mctrk] = row['energy']
        else:
            if mctrk not in track_energy[evt]:
                track_energy[evt][mctrk] = row['energy']
    return map(lambda trks: sum(trks.values()),track_energy.values())


def getMaxEnergyTrack(tracks):
    return max(tracks.values(), key=lambda a: a[0])

def energyPositionTracks(mctracks):
    """
    Find track with max energy for each event and returns energy,x,y,z
    """
    evt = 0
    energies = []
    xs = []
    ys = []
    zs = []
    
    track_energy = {}
    
    for row in mctracks.iterrows():
        if(row['event_indx'] != evt):
            maxEnergyTrk = getMaxEnergyTrack(track_energy[evt])
            track_energy = {}
            energies.append(maxEnergyTrk[0])
            xs.append(maxEnergyTrk[1][0])
            ys.append(maxEnergyTrk[1][1])
            zs.append(maxEnergyTrk[1][2])
            
        evt = row['event_indx']
        mctrk = row['mctrk_indx']
        if evt not in track_energy:
            track_energy[evt] = {}
            track_energy[evt][mctrk] = (row['energy'],row['initial_vertex'])
        else:
            if mctrk not in track_energy[evt]:
                track_energy[evt][mctrk] = (row['energy'],row['initial_vertex'])
                
    #Don't loose last event
    maxEnergyTrk = getMaxEnergyTrack(track_energy[evt])
    track_energy = {}
    energies.append(maxEnergyTrk[0])
    xs.append(maxEnergyTrk[1][0])
    ys.append(maxEnergyTrk[1][1])
    zs.append(maxEnergyTrk[1][2])
    
    return np.array(energies), np.array(xs), np.array(ys), np.array(zs)

def energyHist(energies):
    n, bins, patches = plt.hist(energies, 50, normed=0, facecolor='green', alpha=0.75)
    
    plt.xlabel('E (MeV)')
    plt.ylabel('Events')
    plt.title(r'Events energy')
    #plt.axis([40, 160, 0, 0.03])
    plt.grid(True)
    plt.show()

def energyCorrelation(energies, pmtEnergies):
    plt.scatter(energies, pmtEnergies, s=80, marker="+")
    plt.xlabel('True Energy (MeV)')
    plt.ylabel('Q in PMTs (pe)')
    plt.title('Correlation between true energy and PMT energy')
    plt.show()

def positionCorrelation(pmtEnergies,xs,ys,zs):
    plt.figure(figsize=(8, 20))
    plt.subplot(311)
    plt.xlabel('x (mm)')
    plt.ylabel('Q (pe)')
    plt.title('Correlation between x coordinate and PMT energy')
    plt.grid(True)
    plt.scatter(xs,pmtEnergies, s=80, marker="+")
    plt.subplot(312)
    plt.xlabel('y (mm)')
    plt.ylabel('Q (pe)')
    plt.title('Correlation between y coordinate and PMT energy')
    plt.grid(True)
    plt.scatter(ys,pmtEnergies, s=80, marker="+")
    plt.subplot(313)
    plt.xlabel('z (mm)')
    plt.ylabel('Q (pe)')
    plt.title('Correlation between z coordinate and PMT energy')
    plt.grid(True)
    plt.scatter(zs,pmtEnergies, s=80, marker="+")
    plt.show()
    
def positionCorrelation3d(pmtEnergies,xs,ys,zs):
    fig = plt.figure(1);
    fig.set_figheight(10.);
    fig.set_figwidth(11.);
        
    ax1 = fig.add_subplot(111,projection='3d');
    
    s1 = ax1.scatter(xs,ys,zs,marker='s',linewidth=0.5,alpha=0.6,s=5,c=np.log10(pmtEnergies),cmap=plt.get_cmap('rainbow'));

    ax1.set_xlabel('x (mm)')
    ax1.set_ylabel('y (mm)')
    ax1.set_zlabel('z (mm)')
    plt.title('Distribution of events (energy in color code)')
    plt.grid(True)
    plt.show()

def usage():
    """
    Usage of program
    """
    print("""
        Usage: python (run) ZAIRA [args]
        where args are:
         -h (--help) : this text
         -i (--info) : print a text describing the invisible city of ZAIRA
         -c (--cfile) : full path to a configuration file
         
         example of configuration file 

        #Header is not read
        Names of parameters (comma separated)
        Values of parameters (comma separated)
        
        The parameters for DIOMIRA are:
        PATH_IN,PATH_OUT,FILE_IN,FILE_OUT,FIRST_EVT,LAST_EVT

        The parameters are self-explaining. 

        """)

def configure(argv):
    """
    read arguments from the command line and configures job
    """

    global PATH_IN, PATH_OUT, FILE_IN, FILE_OUT, FIRST_EVT, LAST_EVT
    cfile = ''

    try:
        opts, args = getopt.getopt(argv, "hidc:", ["help","info","cfile"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            print("help")
            sys.exit()
        elif opt in ("-i", "--info"):
            INFO = arg
            print("info")
        elif opt in ("-c", "--cfile"):
            cfile = arg
            print("cfile")
    

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


def zaira():
    h5in = tb.open_file(PATH_IN+'/'+FILE_IN, "r+")
    pmtrd_ = h5in.root.pmtrd
    sipmrd_ = h5in.root.sipmrd
    mctracks_ = h5in.root.MC.MCTracks

    NEVENTS = pmtrd_.shape[0]
    print("Number of events in file: {}".format(NEVENTS))
    
    #Print some values
    print("PMT energies for first 10 events:\n")
    for i in xrange(0, 10):
        print(energy_pes(i,pmtrd_))
        
    #Compute pmt & tracks total energies
    pmtEnergies = map(lambda i: energy_pes(i, pmtrd_).sum(), xrange(NEVENTS))
    totalEnergies = energyTracks(mctracks_)

    tkrEnergies, xs, ys, zs = energyPositionTracks(mctracks_)
    
    #Plot histogram
    energyHist(totalEnergies)
    
    #Energy correlation
    energyCorrelation(totalEnergies,pmtEnergies)
    
    #Position correlation
    positionCorrelation(pmtEnergies,xs,ys,zs)
    positionCorrelation3d(pmtEnergies,xs,ys,zs)
    
    #h5in.close()

    print("Leaving Zaira. Safe travel!")

if __name__ == '__main__':
    configure(sys.argv[1:])
    zaira()
