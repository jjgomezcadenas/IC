
# coding: utf-8

# In[7]:

from __future__ import print_function

import matplotlib.pyplot as plt
import numpy             as np
import h5py
import sys

# Add the Utilties directory to the system path so the file sipm_param can be imported.
sys.path.append("utilities")
from sipm_param import *
TOTALEVTS = 100

# Get the arguments.
args = sys.argv
EVT_START = int(args[1])
EVT_END = int(args[2])

# Define the 3D grid extent and range limits.

# In[8]:

# Range limits from extraction 
NX = 200
NY = 200
NZ = 200

# Projection voxel sizes in mm
vSizeX = 2
vSizeY = 2

# The slice width, in Geant4 voxels
slice_width = 2.

# Range limit in x and y for slices (assuming a square range), in voxels
RNG_LIM = 200 # 196 # less than 200 because in NEW EL plane is smaller than SiPM plane

# SiPM plane geometry definition
nsipm = 20             # number of SiPMs in response map
sipm_pitch = 10.       # distance between SiPMs
sipm_edge_width = 5.   # distance between SiPM and edge of board

# -------------------------------------------------------------------------------------
xlen = 2*sipm_edge_width + (nsipm-1)*sipm_pitch       # (mm) side length of rectangle
ylen = 2*sipm_edge_width + (nsipm-1)*sipm_pitch       # (mm) side length of rectangle
wbin = 2.0                                            # (mm) bin width

# Compute the positions of the SiPMs
pos_x = np.ones(nsipm**2)*sipm_edge_width +          (np.ones(nsipm*nsipm)*range(nsipm**2) % nsipm)*sipm_pitch
pos_y = np.ones(nsipm**2)*sipm_edge_width +          np.floor(np.ones(nsipm*nsipm)*range(nsipm**2) / nsipm)*sipm_pitch


# Load Data

# In[9]:

sg = h5py.File('in_data/dnn_MAGBOX_oldPhys_600k_si_v2x2x2_r200x200x200.h5')
bg = h5py.File('in_data/dnn_MAGBOX_oldPhys_600k_bg_v2x2x2_r200x200x200.h5')


# Create event slicer

# In[10]:

def slice_evt(hfile,nevt,zwidth):
    """
    Create slices for the specified event.
    hfile: the HDF5 files containing the events
    nevt: the event number to slice
    zwidth: the slice width in mm
    
    returns: [energies, slices]
    where energies is a list of the energies in each slice
    and slices is a matrix of size [Nslices,NY,NX] containing normalized slices
    """
    
    # Get the event from the file.
    htrk = hfile['trk{0}'.format(nevt)]
    
    # Get the z-range.
    zmin = np.min(htrk[2]); zmax = np.max(htrk[2])
    
    # Create slices of width zwidth beginning from zmin.
    nslices = int(np.ceil((zmax - zmin)/zwidth)) # make sure works same as in math
    #print("{0} slices for event {1}".format(nslices,nevt))
    
    slices = np.zeros([nslices,NY,NX],dtype=np.float32)
    energies = np.zeros(nslices, dtype=np.float32)
    for x,y,z,e in zip(htrk[0],htrk[1],htrk[2],htrk[3]):
        
        # Add the energy at (x,y,z) to the (x,y) value of the correct slice.
        islice = int((z - zmin)/zwidth)
        if(islice == nslices): islice -= 1
        slices[islice][y][x] += e
        energies[islice] += e
    
    # Normalize the slices to the total energy
    etot = np.sum(energies)
    for s in range(nslices):
        slices[s] /= etot
#        slices[s] /= energies[s]
        
    # Return the list of slices and energies.
    return [energies, slices]


# Define function to write and save data for neural net

# In[11]:

def save_data(fread,fwrite):
    xrng = []; yrng = []   # x- and y-ranges
    nspevt = []            # number of slices per event
    slices_x = []; slices_y = []; slices_e = []   # slice arrays
    Ntrks = EVT_END - EVT_START #min(len(fread),TOTALEVTS) # this = 5k here
    #if Ntrks !=  TOTALEVTS: 
        #print('WARNING: Not enough data. TOTALEVTS > number events in data file...  Preparing: ' + str(Ntrks) + 'events...')
        #print('Preparing ' + str(Ntrks) + ' events...')
    for ee in range(Ntrks):

        if(Ntrks < 100 or (ee + 1) % int(Ntrks/100) == 0):
            print("Slicing event {0}".format(ee+EVT_START))

        # Slice the event.
        en,sl = slice_evt(fread,ee+EVT_START,slice_width)
        nslices = len(en)
        nspevt.append(nslices)
        
        # Get information about each slice.
        for ss in range(nslices):

            # Don't include 0-energy slices.
            if(en[ss] < 0):
                print('skipped_slice')
                continue

            # Get lists of the nonzero x,y,z indices and E values.
            cslice = sl[ss]
            nzy,nzx = np.nonzero(cslice)
            nze = cslice[np.nonzero(cslice)]

            # Extract several quantities of interest.
            xmin = np.min(nzx); xmax = np.max(nzx)
            ymin = np.min(nzy); ymax = np.max(nzy)
            xrng.append(xmax - xmin + 1)
            yrng.append(ymax - ymin + 1)

            # Save the slice if within range.
            if((xmax - xmin) >= RNG_LIM-1 or (ymax - ymin) >= RNG_LIM-1):
                print("Range of {0} for event {1} slice {2}, energy {3}; slice not included".format(xmax-xmin,ee+EVT_START,ss,en[ss]))
            else:

                # Center the slices about RNG_LIM/2.
                x0 = int((xmin + xmax)/2. - RNG_LIM/2.)
                y0 = int((ymin + ymax)/2. - RNG_LIM/2.)
                nzx -= x0; nzy -= y0

                # Create the slice array.
                snum = len(slices_x)
                slices_x.append(nzx); slices_y.append(nzy); slices_e.append(nze)
                carr = np.array([nzx, nzy, nze])

                # Create the corresponding SiPM map.
                sipm_map = np.zeros(nsipm*nsipm,dtype=np.float32)
                for xpt,ypt,ept in zip(nzx,nzy,nze):

                    # Compute the distances and probabilities.  Add the probabilities to the sipm map.
                    rr = np.array([np.sqrt((xi - xpt)**2 + (yi - ypt)**2) for xi,yi in zip(pos_x,pos_y)])
                    probs = 0.5*(sipm_par(0, rr) + sipm_par(1, rr))
                    sipm_map += probs*ept

                # Normalize the probability map, and set sigma = 1.
                #sipm_map -= np.mean(sipm_map)
                #sipm_map /= np.std(sipm_map)

                # Save the slice and the SiPM map to an HDF5 file.
                fwrite.create_dataset("slice{0}".format(snum),data=carr,     dtype='float32')
                fwrite.create_dataset("sipm{0}".format(snum), data=sipm_map, dtype='float32')
    fwrite.create_dataset("nspevt", (Ntrks,), data=np.array(nspevt))
    fwrite.close()


# Create one file for signal, one for background

# In[12]:

sg_out = h5py.File("out_data/signal_2x2x2_{0}_to_{1}.h5".format(EVT_START,EVT_END-1)    , 'w')
bg_out = h5py.File("out_data/background_2x2x2_{0}_to_{1}.h5".format(EVT_START,EVT_END-1), 'w')
print('---generating signal---')
save_data(sg,sg_out)
print('---generating background---')
save_data(bg,bg_out)
print('---finished.---')

