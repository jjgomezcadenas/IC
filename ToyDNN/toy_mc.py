import tables as tb
import numpy  as np
import sys

# Add the Utilties directory to the system path so the file sipm_param can be imported.
sys.path.append("../Utilities")
from sipm_param import *
from toy_inputs import *


# toy_mc.py provides functions for generating and saving, or loading
# montecarlo data for the toy reconstruction.

def load_data(MC_filename):
    f = tb.open_file(MC_filename, 'r')

    # access correct group
    if N_ELpts  == 1: dat_tree = f.root.sim_1pt
    elif N_ELpts ==2: dat_tree = f.root.sim_2pt
    elif N_ELpts ==3: dat_tree = f.root.sim_3pt
    elif N_ELpts ==4: dat_tree = f.root.sim_4pt

    print('-- Loading data from file: ' + MC_filename + ' ...')
    x_train = dat_tree.xtrain[:]
    y_train = dat_tree.ytrain[:]
    x_valid = dat_tree.xvalid[:]
    y_valid = dat_tree.yvalid[:]

    Nevts_train = y_train.shape[0]
    Nevts_valid = y_valid.shape[0]
    Nevts = Nevts_train + Nevts_valid

    print(dat_tree)
    f.close()
    print('-- Load complete')
    return x_train,y_train,x_valid,y_valid,Nevts_train,Nevts_valid

def generate_data():
    print('-- Generating hits and response data ... ')
    # Set up the training sets.
    x_train = np.zeros([Nevts_train,nsipm*nsipm]); x_valid = np.zeros([Nevts_valid,nsipm*nsipm])
    y_train = np.zeros([Nevts_train,2*N_ELpts]);   y_valid = np.zeros([Nevts_valid,2*N_ELpts])

    # generate and sort ypts
    Xpts = np.sort(np.random.uniform(0,xlen,(Nevts,N_ELpts)),axis=1)
    Ypts = np.random.uniform(0,ylen,(Nevts,N_ELpts))

    for nn in range(Nevts):

        if(nn % int(Nevts/10) == 0):
            print "-- Event {0} of {1} ...".format(nn,Nevts)

        # For each event, generate a number of random EL points and compute the probabilities of detection for each SiPM.
        sipm_map = np.zeros(nsipm*nsipm)
        for npt in range(N_ELpts):

            # Generate the point.
            #xpt = rd.random()*xlen
            xpt = Xpts[nn,npt]
            ypt = Ypts[nn,npt]

            # Compute the distances and probabilities.  Add the probabilities to the sipm map.
            rr = np.array([np.sqrt((xi - xpt)**2 + (yi - ypt)**2) for xi,yi in zip(pos_x,pos_y)])
            probs = 0.5*(sipm_par(0, rr) + sipm_par(1, rr))
            sipm_map += probs

            # Fill the y matrix with the generated points.
            if(nn >= Nevts_train):
                y_valid[nn-Nevts_train,2*npt] = xpt/xlen
                y_valid[nn-Nevts_train,2*npt+1] = ypt/ylen
            else:
                y_train[nn,2*npt] = xpt/xlen
                y_train[nn,2*npt+1] = ypt/ylen

        # Normalize the probability map, and set sigma = 1.
        sipm_map -= np.mean(sipm_map)
        sipm_map /= np.std(sipm_map)

        # Fill the x matrix for the generated points.
        if(nn >= Nevts_train):
            x_valid[nn-Nevts_train] = sipm_map
        else:
            x_train[nn] = sipm_map

    return x_train,y_train,x_valid,y_valid

def save_data(MC_filename,x_train,y_train,x_valid,y_valid):
    if MOD_data:
        f = tb.open_file(MC_filename, 'r+')    # modify
        print('-- Modifying '  +  MC_filename  + ' to contain new montecarlo data ...')
    else:
        f = tb.open_file(MC_filename, 'w')      # else write
        print('-- Writing new ' +  MC_filename   + ' to contain new montecarlo data ...')


    # compress tables arrays with blosc, but don't shuffle
    filters = tb.Filters(complib='blosc', complevel=9, shuffle=False)

    groupname = 'sim_' + str(N_ELpts) + 'pt'
    group = f.create_group(f.root, groupname, 'Group for ' + str(N_ELpts) + ' ELPTs')

    # ** The arrays used here are extensible (EArrays), this means that we can add more data later if we want
    #    so that in the future we can train on a larger training set (or validate with larger valid set).
    #    To add one SiPM response map to the training data you would open the table in the .h file,
    #    access the correct group, and then maps to xdat.
    #    For example, to append one map to the training data with 3 E_LPts you would write:
    #    f.root.sim_3pt.xtrain.append([map]) when map is a numpy ndarray with one dimension of size nsipm^2

    atom = tb.Atom.from_dtype(x_train.dtype)
    xtrain = f.create_earray(group, 'xtrain', atom, (0,nsipm**2),  filters=filters, expectedrows=Nevts_train)

    atom = tb.Atom.from_dtype(y_train.dtype)
    ytrain = f.create_earray(group, 'ytrain', atom, (0,2*N_ELpts), filters=filters, expectedrows=Nevts_train)

    atom = tb.Atom.from_dtype(x_valid.dtype)
    xvalid = f.create_earray(group, 'xvalid', atom, (0,nsipm**2),  filters=filters, expectedrows=Nevts_valid)

    atom = tb.Atom.from_dtype(y_valid.dtype)
    yvalid = f.create_earray(group, 'yvalid', atom, (0,2*N_ELpts), filters=filters, expectedrows=Nevts_valid)

    print('-- Saving training data ...')

    #Now put the data in the arrays
    for ev in range(Nevts_train):
        xtrain.append([x_train[ev]])
        ytrain.append([y_train[ev]])

    print('-- Saving validation data ...')

    for ev in range(Nevts_valid):
        xvalid.append([x_valid[ev]])
        yvalid.append([y_valid[ev]])

    print(f)
    f.close()
