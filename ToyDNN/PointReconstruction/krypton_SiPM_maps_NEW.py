
from __future__ import print_function
import tables as tb
import numpy as np

nevts = 500000
EL_rad = 198
N_ELpts = 1


def gen_polar_hits(num, N_ELpts, EL_radius):
    """
    Generate hits in polar coordinates. Using polar coords for convenience,
    since hits must be within EL radius
    """
    yp_ = np.empty((num, N_ELpts * 2), np.float32)
    yp_[:, 0: N_ELpts] = np.random.uniform(0, EL_radius, (num, N_ELpts))
    yp_[:, N_ELpts: 2 * N_ELpts] = np.random.uniform(0, 2 * np.pi,
                                                     (num, N_ELpts))
    return yp_


def cartesian_convert(r, theta, N_ELpts):
    """
    Convert polar coords to cartesian
    """
    yc_ = np.empty((r.shape[0], r.shape[1] * 2), np.float32)
    for i in range(N_ELpts):
        yc_[:, i] = np.cos(theta[:, i]) * r[:, i]  # xcoords
        yc_[:, N_ELpts + i] = np.sin(theta[:, i]) * r[:, i]  # ycoords
    return yc_


# Generate random events in EL plane
polar_coords = gen_polar_hits(nevts, N_ELpts, EL_rad)
coords = cartesian_convert(polar_coords[:, :N_ELpts],
                           polar_coords[:, N_ELpts:], N_ELpts)

# Round for easy extraction from the light table
rcoords = np.array(np.round(coords), dtype=np.int)

# Open the light table
f = tb.open_file('ReproducedFull.h5', 'r')
table = f.root.Probabilities.data
grid = np.array(table[:]['grid_xy'], dtype=np.int)
probs = np.array(table[:]['sens_prob'], dtype=np.float32)
ids = np.array(table[:]['sens_id'], dtype=np.int)
id_pos = np.array(f.root.Sensors.XY[12:])

# *Note table where will not work, so replace ids with positions
# **takes 30 seconds**
pos = np.ones((len(ids), 2), dtype=np.int) * -9999  # -9999 --> PMTs
for p in id_pos:
    pos[np.where(ids == p[0])[0]] = p[1:]

maps = np.zeros((nevts, 48, 48), dtype=np.float32)

# For each time slice: collect non zero sipm responses, populate maps
for z, xy in enumerate(rcoords):

    # Collect non zero sipm responses
    i = np.where((grid == xy).all(axis=1))[0]
    sli_pos = pos[i]
    sli_idx = np.array((sli_pos + 235) / 10, dtype=np.int)
    sli_probs = probs[i]

    # Populate sipm maps for this time slice
    for xidx, yidx, probj in zip(sli_idx[:, 0], sli_idx[:, 1], sli_probs):
        if xidx != -977:
            maps[z, xidx, yidx] = probj

    if (z + 1) % 1000 == 0:
        print(str(z/float(nevts)) + ' complete...')

# Save to file
f = tb.open_file('NEW_krypton_maps_500000b.h', 'w')
filters = tb.Filters(complib='blosc', complevel=9, shuffle=False)

atom = tb.Atom.from_dtype(maps.dtype)
tmaps = f.create_earray(f.root, 'maps',   atom, (0, 48, 48), filters=filters)

atom = tb.Atom.from_dtype(coords.dtype)
tcoords = f.create_earray(f.root, 'coords', atom, (0, 2), filters=filters)

for i in range(nevts):
    tmaps.append([maps[i]])
    tcoords.append([coords[i]])  # Save unrounded coords!

print(f)
f.close()
