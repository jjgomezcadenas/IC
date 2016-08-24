"""
MCWFControlPlots
JMB August 2016

What MCWFControlPlots does:
1) Read a file containing MC waveforms (25 ns)
2) Energy histogram
3) Plot real&measured energy
4) Correlation plots between coordinates and energy
"""

from __future__ import print_function

from Util import *

import matplotlib.pyplot as plt
import pandas as pd
import sys
import tables as tb
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import Configure as CF

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



def mcwfcontrolplots(CFP):
    PATH_IN=CFP['PATH_IN']
    FILE_IN=CFP['FILE_IN']
    FIRST_EVT=CFP['FIRST_EVT']
    LAST_EVT=CFP['LAST_EVT']

    h5in = tb.open_file(PATH_IN+'/'+FILE_IN, "r+")
    pmtcwf_ = h5in.root.RD.pmtcwf
    sipmrwf_ = h5in.root.RD.sipmrwf
    mctracks_ = h5in.root.MC.MCTracks

    NEVENTS = pmtcwf_.shape[0]
    print("Number of events in file: {}".format(NEVENTS))
    
    #Print some values
    print("PMT energies for first 10 events:\n")
    for i in xrange(0, 10):
        print(energy_pes(i,pmtcwf_))
        
    #Compute pmt & tracks total energies
    pmtEnergies = map(lambda i: energy_pes(i, pmtcwf_).sum(), xrange(NEVENTS))
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
    INFO, CFP = CF.configure(sys.argv[0], sys.argv[1:])
    mcwfcontrolplots(CFP)
