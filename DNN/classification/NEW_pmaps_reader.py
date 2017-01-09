# Modified version of New_kr_pmaps_reader
# NEW_pmaps_reader reads pmaps and constructs true labels
# and corresponding sipm maps for the DNN
from __future__ import print_function
import tables as tb
import numpy as np
import h5py
import sys


# -------------------------------------------------------------------------------------------------
# Input parameters
#

# Process files with base dfile and ending _(num).h5 where (num) runs from 0 to fnum-1.
#dfile = '/home/jrenner/data/SE/hdf5_NEXT_v0_08_06_NEW_se_1M_combined.h5'
#dfile = '/home/jrenner/data/0vbb/bb_1M_v0_08_07/hdf5_NEXT_v0_08_06_NEW_bb_1M_combined'
#dfile = '/home/jrenner/data/0vbb/bb_1M_v0_08_07/hdf5_NEXT_NEW_bb_1M_v0_08_07'
dfile = '/home/jrenner/data/SE/se_1M_v0_08_07/hdf5_NEXT_NEW_se_1M_v0_08_07'
#fnum = 50

sipm_param_file = '/home/jrenner/data/sipm_param/ReproducedFull.h5'

nevts = 100000                              # number of events to read
tbin = 2                                 # bin size in microseconds
max_slices = 60                          # maximum number of slices per map
# -------------------------------------------------------------------------------------------------

# Get input arguments
iargs = sys.argv
if(len(iargs) < 4):
    print("python NEW_pmaps_reader.py <si_or_bg> <start_file> <end_file+1>")
    exit()
etype = iargs[1]          # 'si' or 'bg': set depending on type of events in file below (dfile)
f_start = int(iargs[2])
f_end = int(iargs[3])

# Get the map of XY values vs. ID.
f_ids = tb.open_file(sipm_param_file)
s_ids = np.array(f_ids.root.Sensors.XY)
f_ids.close()

# Open the 20x20 window table.
fwtbl = h5py.File("wtbl.h5",'r')
wtbl = fwtbl['wtbl']

# Construct a dictionary that will map sensor id to a x,y coord pair
id_to_coords = {}
for ID, x, y in zip(s_ids[12:, 0], s_ids[12:, 1], s_ids[12:, 2]):
    id_to_coords[np.int32(ID)] = np.array([-x, y])  # -x b/c old IDs!

# Loop over all files.
maps = []
energies = []
for fn in range(f_start,f_end):

    # Open the hdf5 file containing the PMaps.
    fname = "{0}_{1}.h5".format(dfile,fn)
    print("Opening file: {0}".format(fname))
    fpmap = tb.open_file(fname,'r')
    pmaps = fpmap.root.PMaps.PMaps
    print(pmaps)

    # Loop over all events.
    rnum = 0                    # row number in table iteration
    ev = 0                      # processed event number (starts at 0 and runs to nevts-1)
    evtnum = pmaps[0]['event']  # event number from file
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
            print("* Event contains", nslices, "slices.")
            if(nslices >= max_slices):
                print("WARNING: event with ", nslices, " slices too large for max bins =", max_slices, ". Increase bin size.")
        
            else:
        
                # Calculate a slice offset so that not all slices begin at 0.
                offset = 0 #np.random.randint(0,max_slices-nslices) 
            
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

                # -----------------------------------------------------------
                # Attempt to confine the event to a 20x20 window.
                tmap20 = np.zeros((20, 20, max_slices), dtype=np.float32)
                tsproj = np.zeros((48, 48), dtype=np.float32)
                tsipm_sum = 0.
                for isipm in range(48):
                    for jsipm in range(48):
                        sval = np.sum(tmap[isipm, jsipm, :])
                        tsproj[isipm, jsipm] = sval
                        tsipm_sum += sval
                print(tsproj)

                # Get the index of the max SiPM and fill the corresponding 20x20 maps.
                imax_sum = np.argmax(tsproj)
                print("** Max is",imax_sum)
                for nslice in range(len(eslices)):
                    for s20id in range(len(wtbl[imax_sum])):
                        s48id = wtbl[imax_sum][s20id]
                        i20 = int(s20id / 20)
                        j20 = s20id % 20
                        i48 = int(s48id / 48)
                        j48 = s48id % 48
                        tmap20[i20][j20][nslice] = tmap[i48][j48][nslice]

                # Ensure that the charge in the 20x20 window is sufficient.
                chg_frac = np.sum(tmap20) / tsipm_sum
                if(chg_frac > 0.9):                         
 
                    # Weight the slices assuming a total slice energy of 1.
                    etot = np.sum(eslices)
                    for nsl,esl in enumerate(eslices):
                        nval = np.sum(tmap20[:,:,nsl])
                        if(nval > 0):
                            tmap20[:,:,nsl] *= 1.0/nval
                            tmap20[:,:,nsl] *= 1.0*esl/etot
            
                    # Add the event to the lists of maps and energies if it contains some information.
                    maxval = np.max(tmap)
                    if(maxval > 0):
                        #tmap /= maxval     # normalize so maximum value is = 1
                        maps.append(tmap20)
                        energies.append(eslices)
                        ev += 1
                        print("** Added map for event ", ev, ", file", fn, "; file event number ", evtnum, "\n")
                else:
                    print("--> Event did not make charge fraction cut with charge",chg_frac)
    
        # Set to the next event.
        if(rnum < pmaps.nrows):
            evtnum = pmaps[rnum]['event']

    fpmap.close()

# Save the data
print('...Saving data...')
maps = np.array(maps)
energies = np.array(energies)
print("Maps of length = ", len(maps))

f = tb.open_file('NEW_training_MC_{0}_{1}_to_{2}.h5'.format(etype,f_start,f_end), 'w')
filters = tb.Filters(complib='blosc', complevel=9, shuffle=False)
atom_m = tb.Atom.from_dtype(maps.dtype)
maparray = f.create_earray(f.root, 'maps', atom_m, (0, 20, 20, max_slices), filters=filters)
atom_e = tb.Atom.from_dtype(energies.dtype)
earray = f.create_earray(f.root, 'energies', atom_e, (0, max_slices), filters=filters)

for i in range(len(maps)):
    maparray.append([maps[i]])
    earray.append([energies[i]])

print('...Save complete.')
print(f)
f.close()
