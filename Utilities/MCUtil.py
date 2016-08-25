import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

def energy_tracks(mctracks):
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

def max_energy_track(tracks):
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
            maxEnergyTrk = max_energy_track(track_energy[evt])
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
    maxEnergyTrk = max_energy_track(track_energy[evt])
    track_energy = {}
    energies.append(maxEnergyTrk[0])
    xs.append(maxEnergyTrk[1][0])
    ys.append(maxEnergyTrk[1][1])
    zs.append(maxEnergyTrk[1][2])
    
    return np.array(energies), np.array(xs), np.array(ys), np.array(zs)


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
