import numpy  as np


#--- Parameters for toy_reconstruction.py
# Number of EL points to generate per event.
N_ELpts = 1

# Geometry definition.
nsipm = 8
sipm_pitch = 10.       # distance between SiPMs
sipm_edge_width = 5.   # distance between SiPM and edge of board

# Variables for computing an EL point location.
xlen = 2*sipm_edge_width + 7*sipm_pitch       # (mm) side length of rectangle
ylen = 2*sipm_edge_width + 7*sipm_pitch       # (mm) side length of rectangle
if xlen == ylen: max_xy = xlen
wbin = 2.0                                    # (mm) bin width

# Compute the positions of the SiPMs.
pos_x = np.ones(nsipm**2)*sipm_edge_width + (np.ones(nsipm*nsipm)*range(nsipm**2) % nsipm)*sipm_pitch
pos_y = np.ones(nsipm**2)*sipm_edge_width + np.floor(np.ones(nsipm*nsipm)*range(nsipm**2) / nsipm)*sipm_pitch


# Below are a few switches that determine how the data used/created here are loaded/stored.
Save_model = False     # Save the neural net (including the learned weights) for future training/teseting
Save_hist  = False      # Set to True to save the history of the model
MOD_hist   = False     # if MOD_hist, script will modify existing HISTORY .h file, adding a node in pytable
                       # with loss history
                       # if not MOD_hist, will create a new .h file and do the same thing
                       # **Note if Save_hist is False, MOD_hist does not matter

Load_MC = True    # set to true to load existing SiPM response data

# Load events
if Load_MC: Save_MC = False
# Generate events
else:
    Nevts_train = 10000
    Nevts_valid = 1000
    Nevts = Nevts_train + Nevts_valid
    Save_MC = True     # (do not save data that is already saved)

if Save_MC:
    # Set to true to create a new group (for a new N_ELpts) in an existing .h file
    # Set to false to create a new .h file entirely (will overwrite existing file)
    MOD_data = False
MC_filename =  'TMC_data.h'

# set to true to implement early stopping
ES = False
if ES == True:
    # must specify number of total layers now (not -1)
    N_layers = -1
    ES_filepath =  'ES_' + str(N_ELpts) + 'hits' + str(N_layers) + 'lay.h'
