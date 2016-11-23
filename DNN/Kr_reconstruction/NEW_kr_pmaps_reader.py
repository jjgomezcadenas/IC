# NEW_pmaps_reader reads pmaps and constructs true labels
# and corresponding sipm maps for the DNN
from __future__ import print_function
import tables as tb
import numpy as np

extend_training_set = False
gen_all = True

montecarlo = tb.open_file(
             '/home/jrenner/data/kr/hdf5_NEXT_v0_08_06_NEW_kr_1pt5M.h5',
             'r')
tracks = montecarlo.root.MC.MCTracks
pmaps = montecarlo.root.PMaps.PMaps
reco = montecarlo.root.Reco.Reco

# Can't load too many events into memory: load max 1e6
# Note that is more than 1e6 tracks

# Start reading file at start_event
start_event = 0

# End_event non-inclusive
if gen_all:
    end_event = np.max(tracks[:]['event_indx']) + 1
else:
    end_event = min(start_event + int(1e4),
                    np.max(tracks[:]['event_indx']) + 1)

start_found = False
i_start = -1
i_end = -1

# Find the indices for start and end events
for i, ev in enumerate(tracks[:]['event_indx']):
    if ev == start_event and not start_found:
        i_start = i
        start_found = True
    elif ev == end_event:
        i_end = i
        break

# Construct true labels from montecarlo track information
event_indx = np.array(tracks[i_start:i_end]['event_indx'], dtype=np.int32)
hit_energy = np.array(tracks[i_start:i_end]['hit_energy'], dtype=np.float32)
hit_pos = np.array(tracks[i_start:i_end]['hit_position'], dtype=np.float32)

lbl = np.zeros((np.max(event_indx) + 1, 3), dtype=np.float32)

disc = []  # list of event indices to discard

tot_energy = 0

# total energy for any event should be .041
print('---')
print('Total event energy: ' + str(
     np.sum(hit_energy[np.where(event_indx == 10)[0]])))
print('---')
print('')
print('...Finding objective truth...')

no_label = []

# Calculate "true" labels for dnn
for i, ev in enumerate(event_indx):

    # Compute average hit_pos, weighted by hit_energy
    lbl[ev] = lbl[ev] + hit_energy[i] * hit_pos[i]
    tot_energy += hit_energy[i]

    # Will throw index ev is last event number
    try:
        if event_indx[i + 1] > ev:
            lbl[ev] /= tot_energy
            tot_energy = 0
            if event_indx[i + 1] != ev + 1:
                for i in range(ev + 1, event_indx[i + 1]):
                    no_label.append(i)
    except IndexError:
        if ev == np.max(event_indx):
            lbl[ev] /= tot_energy
            tot_energy = -99999

        # This should not happen, so throw error
        else:
            print('ERROR: ' + str(event_indx[i + 1]))
print('...Found.')

# Now create the corresponding maps
if not gen_all:
    ip_start = pmaps.get_where_list('event == start_event', sort=True)[0]
    print('---PMAP start found...')
    ip_end = pmaps.get_where_list('event == end_event', sort=True)[0]
    print('---PMAP end found...')

    print('---Loading events...')
    event_indx = np.array(pmaps[ip_start:ip_end]['event'], dtype=np.int32)
    amp_list = np.array(pmaps[ip_start:ip_end]['anode'], dtype=np.float32)
    montecarlo.close()
    print('Load Complete.')

f_ids = tb.open_file('ReproducedFull.h5')
s_ids = np.array(f_ids.root.Sensors.XY)
f_ids.close()

# Construct a dictionary that will map sensor id to a x,y coord pair
id_to_coords = {}
for ID, x, y in zip(s_ids[12:, 0], s_ids[12:, 1], s_ids[12:, 2]):
    id_to_coords[np.int32(ID)] = np.array([-x, y])  # -x b/c old IDs!

if not gen_all:
    print('...Populating SiPM maps...')
    # Populate SiPM maps for each event
    maps = np.zeros((np.max(event_indx) + 1, 48, 48), dtype=np.float32)
    for ev, amps in zip(event_indx, amp_list):
        for sipm, amp in enumerate(amps):
            if amp != 0.0:
                ID = sipm % 64 + 1000 * (int(sipm / 64) + 1)
                [i, j] = (id_to_coords[ID] + 235) / 10
                maps[ev, np.int8(i), np.int8(j)] += amp


if gen_all:
    prevprint = -1

    maps = np.zeros((np.max(event_indx) + 1, 48, 48), dtype=np.float32)
    print('...Populating SiPM maps...')

    # Populate SiPM maps for each event
    for row in pmaps.iterrows():
        for sipm, amp in enumerate(row['anode']):
            if amp != 0.0:
                ID = sipm % 64 + 1000 * (int(sipm / 64) + 1)
                [i, j] = (id_to_coords[ID] + 235) / 10
                maps[row['event'], np.int8(i), np.int8(j)] += amp
        if (row['event'] + 1) % 10000 == 0:
            if prevprint != row['event']:
                print(str(float(row['event'] + 1) / 1e6 * 100) + '% complete.')
                prevprint = row['event']


# Delete events that have blank maps or no truth label
print('...Deleting blank SiPM maps...')
empty_maps = []

for ev, m in enumerate(maps):
    if np.max(m) == 0.0:
        empty_maps.append(ev)

delete_events = list(set(empty_maps).union(no_label))

maps = np.delete(maps, delete_events, axis=0)
lbl = np.delete(lbl, delete_events, axis=0)
print(str(len(delete_events)) + ' events deleted.')
print('Training set dimensions: ')
print(maps.shape, lbl.shape)

# Check no blank maps
for m in maps:
    if np.max(m) == 0:
        print('BROKEN')
        break

# Save the data
print('...Saving data...')

if extend_training_set:
    f = tb.open_file('NEW_training_diff_MC_1pt5M2.h', 'r+')
    tmaps = f.root.maps
    tcoords = f.root.coords
else:
    f = tb.open_file('NEW_training_diff_MC_1pt5M2.h', 'w')
    filters = tb.Filters(complib='blosc', complevel=9, shuffle=False)
    atom = tb.Atom.from_dtype(maps.dtype)
    tmaps = f.create_earray(f.root, 'maps', atom, (0, 48, 48), filters=filters)
    atom = tb.Atom.from_dtype(lbl.dtype)
    tcoords = f.create_earray(f.root, 'coords', atom, (0, 3), filters=filters)

for i in range(len(lbl)):
    tmaps.append([maps[i]])
    tcoords.append([lbl[i]])

print('...Save complete.')
print(f)
f.close()
