from __future__ import print_function
import random
import tables
import numpy as np
from sipm_param import *

MOD = False #set mod to true to add a new group to an existing .h file
filename = 'forrelu.h'
NSIPM = 8
NUM_ELPT = 2  # max num points observed in EL in same timestep
#Vis = True

sipm_pitch = 1.0     # SiPM pitch in mm
sipm_edge_width = 5.0 # width of edge of dice board in mm
ze = 10.0             # distance between SiPM plane and EL gap
grid_space = 2.0      # grid spacing in mm
d_gap = 5.0           # length of EL gap
n_tbins = 2           # number of time bins collected as electron crosses EL

N = 1   # Number of photons (probably not necessary due to normalization)

nevts = 100000 # num events for each num EL pt

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
    #elpt = np.random.randint(0,1600,nevts)
    #x[i] = (elpt % max_p)*grid_space + 1
    #y[i] = (np.floor(elpt/max_p))*grid_space + 1
    x[i] = np.random.uniform(0,80,nevts)
    y[i] = np.random.uniform(0,80,nevts)

if NUM_ELPT == 2:
    x[0] = np.random.uniform(0,max_xy/2,nevts)
    x[1] = np.random.uniform(max_xy/2, max_xy, nevts)

sipm_pos = np.array([sipm_pos_x,sipm_pos_y])
# generate responses of sipms
p = 0
#for each EL point (NUM_ELPT)...
for e in E:
    xy = zip(x[p],y[p]) #get coord pair for this elpt

    evt = 0
    # for each event (xi,yi) are its coordinates...
    for xi,yi in xy:
        r = np.sqrt((sipm_pos[0] - xi)**2 + (sipm_pos[1] - yi)**2)
        #add this response to previous responses weighted by this event's energy fraction
        sipm_res[evt] = sipm_res[evt] + (sipm_par(0,r) + sipm_par(1,r))*e[evt]
        evt += 1 #iterate thru events

    p += 1 #iterate thru points

#multiply by number of photons, not neccessary when assuming infinite photons (set N = 1)
sipm_res = sipm_res * N/2

#Normalize for the DNN
mean = np.mean(sipm_res)
mean = 0
std = np.std(sipm_res)
sipm_res = (sipm_res - mean)/std

#Save responses in a table
if MOD: f = tables.open_file(filename, 'r+') #modify
else:  f = tables.open_file(filename, 'w')   #else write

filters = tables.Filters(complib='blosc', complevel=9, shuffle=False)

groupname = 'sim_' + str(NUM_ELPT) + 'pt'
group = f.create_group(f.root, groupname, 'Group for ' + str(NUM_ELPT) + ' ELPTs')
###Initialize arrays:
#x
atom = tables.Atom.from_dtype(x.dtype)
x1 = f.create_earray(group, 'x', atom, (0,nevts), filters=filters) # earrays are extensible along first index
                                            # (data for an additional eltp point
                                            # can be added if deseired)
#y
atom = tables.Atom.from_dtype(y.dtype)
y1 = f.create_earray(group, 'y', atom, (0,nevts), filters=filters)


#E
atom = tables.Atom.from_dtype(E.dtype)
E1 = f.create_earray(group, 'E', atom, (0,nevts), filters=filters)

#sipm
atom = tables.Atom.from_dtype(sipm_res.dtype)  #sipm_res1 probably doesn't need to be extendable
sipm_res1 = f.create_earray(group, 'sipm_resp', atom, (0,nevts,NSIPM**2), filters=filters)

#Apppend to arrays:
for p in range(NUM_ELPT):
    x1.append([x[p]])
    y1.append([y[p]])
    E1.append([E[p]])
sipm_res1.append([sipm_res])

print(f)

f.close()
