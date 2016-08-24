from __future__ import print_function
import random
import tables
import numpy as np
from sipm_param import *

NSIPM = 8
NUM_ELPT = 1     # max num points observed in EL per event


sipm_pitch = 10.0     # SiPM pitch in mm
sipm_edge_width = 5.0 # width of edge of dice board in mm
ze = 10.0             # distance between SiPM plane and EL gap
grid_space = 2.0      # grid spacing in mm
d_gap = 5.0           # length of EL gap
n_tbins = 2           # number of time bins collected as electron crosses EL

N = 1   # Number of photons (probably not necessary due to normalization)

nevts = 50000 # num events for each num EL pt

max_xy = (NSIPM-1)*sipm_pitch + 2*sipm_edge_width # maximum x and y value (80 mm)
max_p = max_xy /grid_space                        # number of points per line (40)

# initialize
sipm_res = np.zeros((nevts,NSIPM**2),dtype=np.float64) #SiPM responses
x = np.empty((NUM_ELPT,nevts),dtype = np.float64)      #xcoord for evt
y = np.empty((NUM_ELPT,nevts),dtype = np.float64)      #ycoord for evt
E = np.empty((NUM_ELPT,nevts),dtype = np.float64)      #energy fraction for evt


#Set up the SiPM positions
sipm_pos_x = []
for i in range(NSIPM): sipm_pos_x.extend([sipm_edge_width + i*sipm_pitch]*NSIPM)
sipm_pos_y = [sipm_edge_width + i*sipm_pitch for i in range(NSIPM)]*NSIPM

#Generate energy fractions for each point from uniform distribution so that they sum to 1
E = np.random.uniform(0,1,(NUM_ELPT,nevts))
E = E / E.sum(axis=0)

#Generate x,y coords for each elpt
for i in range(NUM_ELPT):
    elpt = np.random.randint(0,1600,nevts)
    x[i] = (elpt % max_p)*grid_space + 1
    y[i] = (np.floor(elpt/max_p))*grid_space + 1

#Generate SiPM responses
posxy = zip(sipm_pos_x,sipm_pos_y)
p = 0
#for each EL point (NUM_ELPT)...
for e in E:
    xy = zip(x[p],y[p])

    evt = 0
    # for each event (xi,yi) are its coordinates...
    for xi,yi in xy:
        idx = 0

        # for each sipm get sipm_par.
        for posx,posy in posxy:

            #calculate distance between this sipm and xi,yi
            r = np.sqrt((posx - xi)**2 + (posy - yi)**2)

            #add this response to previous responses weighted by this event's energy fraction
            sipm_res[evt,idx] = sipm_res[evt,idx] +  sipm_par(0,r)*e[evt]
            if n_tbins == 2:
                sipm_res[evt,idx] = sipm_res[evt,idx] +  sipm_par(1,r)*e[evt]

            idx += 1 #iterate thru SIPMs
        evt += 1 #iterate thru events
    p += 1 #iterate thru points

#multiply by number of photons, not neccessary when assuming infinite photons (set N = 1)
if n_tbins == 2: sipm_res = sipm_res * N / 2 #average time of 2 bins
else: sipm_res *= N

#Normalize for the DNN
mean = np.mean(sipm_res)
std = np.std(sipm_res)
sipm_res = (sipm_res - mean)/std

#Save responses in a table
# Store "x" in a chunked array with level 5 BLOSC compression...
f = tables.open_file('resp.h', 'w')

filters = tables.Filters(complib='blosc', complevel=9)

groupname = 'sim_' + str(NUM_ELPT) + 'pt'
group = f.create_group(f.root, groupname, 'Group for ' + str(NUM_ELPT) + ' ELPTs')

#x
atom = tables.Atom.from_dtype(x.dtype)
x1 = f.create_earray(group, 'x', atom, (0,nevts), filters=filters) # earrays are extensible along first index
                                            # (data for an additional eltp point
                                            # can be added if deseired)                                                                                                                                                                               
#y
atom = tables.Atom.from_dtype(y.dtype)
y1 = f.create_earray(group, 'y', atom, (0,nevts), filters=filters)


#sipm
atom = tables.Atom.from_dtype(sipm_res.dtype)  #sipm_res1 probably doesn't need to be extendable
sipm_res1 = f.create_earray(group, 'sipm_resp', atom, (0,nevts,NSIPM**2), filters=filters)


for p in range(NUM_ELPT):
    x1.append([x[p]])
    y1.append([y[p]])
sipm_res1.append([sipm_res])

print(f)

f.close()
