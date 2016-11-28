# Modified version of New_kr_pmaps_reader
# NEW_pmaps_reader reads pmaps and constructs true labels
# and corresponding sipm maps for the DNN
from __future__ import print_function
import tables as tb
import numpy as np

extend_training_set = False
gen_all = True

# -------------------------------------------------------------------------------------------------
# Input parameters
#
etype = 'bg'                  # 'si' or 'bg': set depending on type of events in file below (dfile) 
dfile = '/home/jrenner/data/SE/hdf5_NEXT_v0_08_06_NEW_se_1M_combined.h5'
#dfile = '/home/jrenner/data/0vbb/hdf5_NEXT_v0_08_06_NEW_bb_1M_combined.h5'
sipm_param_file = '/home/jrenner/data/sipm_param/ReproducedFull.h5'

nevts = 10000                            # number of events to read
tbin = 5                                 # bin size in microseconds
max_slices = 30                          # maximum number of slices per map
# -------------------------------------------------------------------------------------------------

# Open the hdf5 file containing the PMaps.
fpmap = tb.open_file(dfile,'r')
pmaps = fpmap.root.PMaps.PMaps
reco = fpmap.root.Reco.Reco

# Get the map of XY values vs. ID.
f_ids = tb.open_file(sipm_param_file)
s_ids = np.array(f_ids.root.Sensors.XY)
f_ids.close()

# Construct a dictionary that will map sensor id to a x,y coord pair
id_to_coords = {}
for ID, x, y in zip(s_ids[12:, 0], s_ids[12:, 1], s_ids[12:, 2]):
    id_to_coords[np.int32(ID)] = np.array([-x, y])  # -x b/c old IDs!

# Loop over all events.
rnum = 0                    # row number in table iteration
ev = 0                      # processed event number (starts at 0 and runs to nevts-1)
evtnum = pmaps[0]['event']  # event number from file
maps = []
energies = []
while(rnum < pmaps.nrows and ev < nevts):

    print("Attempting to process event ", evtnum, "...")

    # Attempt to get all times, cathode energies, and anode values
    #  for one event.
    times = []
    cathodes = []
    anodes = []
    while(rnum < pmaps.nrows and pmaps[rnum]['event'] == evtnum):
        if(pmaps[rnum]['signal'] == b'S2'):
            times.append(pmaps[rnum]['time'])
            cathodes.append(pmaps[rnum]['cathode'])
            anodes.append(pmaps[rnum]['anode'])
        rnum += 1
 
    # If we had an S2 for this event, add these values to an SiPM map.
    if(len(times) > 0):

        t0 = times[0]
        tf = times[-1]
        nslices = int((tf - t0) / tbin)
        if((1.0*(tf - t0) / tbin) - nslices > 0):     # if any remainder, add 1 slice
            nslices += 1
      
        # Warn and note if the event contains more than max_slices bins. 
        if(nslices >= max_slices):
            print("WARNING: event with ", nslices, " slices too large for max bins =", max_slices, ". Increase bin size.")
    
        else:
    
            # Calculate a slice offset so that not all slices begin at 0.
            offset = np.random.randint(0,max_slices-nslices) 
        
            # Create SiPM maps for this event.
            tmap = np.zeros((48, 48, max_slices), dtype=np.float32)
            eslices = np.zeros(max_slices)
            for time,cathode,amps in zip(times, cathodes, anodes):
        
                print("-- Processing slice with time ", time, " and energy ", cathode)
        
                # Calculate the bin number.
                bb = int((time - t0) / tbin)
        
                # Add the energy to the energy array.
                eslices[bb + offset] += cathode
        
                # Add the amplitudes to the map.
                for sipm, amp in enumerate(amps):
                    if amp != 0.0:
                        ID = sipm % 64 + 1000 * (int(sipm / 64) + 1)
                        [i, j] = (id_to_coords[ID] + 235) / 10
                        tmap[np.int8(i), np.int8(j), bb + offset] += amp
        
            # Weight the slices assuming a total slice energy of 1.
            etot = np.sum(eslices)
            for nsl,esl in enumerate(eslices):
                tmap[:,:,nsl] *= 1.0*esl/etot
        
            # Add the event to the lists of maps and energies if it contains some information.
            maxval = np.max(tmap)
            if(maxval > 0):
                tmap /= maxval     # normalize so maximum value is = 1
                maps.append(tmap)
                energies.append(eslices)
                ev += 1
                print("** Added map for event ", ev, "; file event number ", evtnum, "\n")

    # Set to the next event.
    if(rnum < pmaps.nrows):
        evtnum = pmaps[rnum]['event']

# Save the data
print('...Saving data...')
maps = np.array(maps)
energies = np.array(energies)
print("Maps of length = ", len(maps))

f = tb.open_file('NEW_training_MC_{0}.h5'.format(etype), 'w')
filters = tb.Filters(complib='blosc', complevel=9, shuffle=False)
atom_m = tb.Atom.from_dtype(maps.dtype)
maparray = f.create_earray(f.root, 'maps', atom_m, (0, 48, 48, max_slices), filters=filters)
atom_e = tb.Atom.from_dtype(energies.dtype)
earray = f.create_earray(f.root, 'energies', atom_e, (0, max_slices), filters=filters)

for i in range(len(maps)):
    maparray.append([maps[i]])
    earray.append([energies[i]])

print('...Save complete.')
print(f)
f.close()
