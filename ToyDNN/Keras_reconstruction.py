
# coding: utf-8

# # Reconstruction with Keras
#

# In[ ]:

from keras.models       import Sequential
from keras.layers       import Dense, Activation, Dropout
from keras.optimizers   import SGD
from keras              import callbacks

from matplotlib.patches import Ellipse

import matplotlib.pyplot as plt
import numpy  as np
import random as rd
import tables as tb
import math
import sys


# Add the Utilties directory to the system path so the file sipm_param can be imported.
sys.path.append("../Utilities")
from sipm_param import *

MOD = True       # if MOD, script will modify existing HISTORY .h file, adding a node in pytable
                 # with loss history
                 # if not MOD, will create a new .h file and do the same thing

Save = True      # Set to True to save the model


# The parameterization is done in two time bins and can be accessed through the function sipm_par(tbin, r) defined in Utilities/sipm_param.py.  Here tbin is the time bin and r is the radial distance from the EL point at which the SiPM of interest is located.

# In[ ]:

# Geometry definition.
nsipm = 8
sipm_pitch = 10.       # distance between SiPMs
sipm_edge_width = 5.   # distance between SiPM and edge of board

# Number of EL points to generate per event.
N_ELpts = 1

ES = False            # set to true to implement early stopping
if ES == True:
    N_layers = 4       # must specify number of total layers now
    ES_filepath =  'ES_' + str(N_ELpts) + 'hits' + str(N_layers) + 'lay.h'
    print(ES_filepath)

# Variables for computing an EL point location.
xlen = 2*sipm_edge_width + 7*sipm_pitch       # (mm) side length of rectangle
ylen = 2*sipm_edge_width + 7*sipm_pitch       # (mm) side length of rectangle
if xlen == ylen: max_xy = xlen
wbin = 2.0                                    # (mm) bin width

# Compute the positions of the SiPMs.
pos_x = np.ones(nsipm**2)*sipm_edge_width + (np.ones(nsipm*nsipm)*range(nsipm**2) % nsipm)*sipm_pitch
pos_y = np.ones(nsipm**2)*sipm_edge_width + np.floor(np.ones(nsipm*nsipm)*range(nsipm**2) / nsipm)*sipm_pitch


# In[ ]:

# Number of events to generate.
Nevts_train = 100000
Nevts_test = 10000

# Set up the training sets.
x_train = np.zeros([Nevts_train,nsipm*nsipm]); x_test = np.zeros([Nevts_test,nsipm*nsipm])
y_train = np.zeros([Nevts_train,2*N_ELpts]); y_test = np.zeros([Nevts_test,2*N_ELpts])

# Generate the events.
Nevts = Nevts_train + Nevts_test

# generate and sort ypts
Xpts = np.sort(np.random.uniform(0,xlen,(Nevts,N_ELpts)),axis=1)

for nn in range(Nevts):

    if(nn % int(Nevts/10) == 0):
        print "-- Event {0} of {1} ...".format(nn,Nevts)

    # For each event, generate a number of random EL points and compute the probabilities of detection for each SiPM.
    sipm_map = np.zeros(nsipm*nsipm)
    for npt in range(N_ELpts):

        # Generate the point.
        #xpt = rd.random()*xlen
        xpt = Xpts[nn,npt]
        ypt = rd.random()*ylen

        # Compute the distances and probabilities.  Add the probabilities to the sipm map.
        rr = np.array([math.sqrt((xi - xpt)**2 + (yi - ypt)**2) for xi,yi in zip(pos_x,pos_y)])
        probs = 0.5*(sipm_par(0, rr) + sipm_par(1, rr))
        sipm_map += probs

        # Fill the y matrix with the generated points.
        if(nn >= Nevts_train):
            y_test[nn-Nevts_train][2*npt] = xpt/xlen
            y_test[nn-Nevts_train][2*npt+1] = ypt/ylen
        else:
            y_train[nn][2*npt] = xpt/xlen
            y_train[nn][2*npt+1] = ypt/ylen

    # Normalize the probability map, and set sigma = 1.
    sipm_map -= np.mean(sipm_map)
    sipm_map /= np.std(sipm_map)

    # Fill the x matrix for the generated points.
    if(nn >= Nevts_train):
        x_test[nn-Nevts_train] = sipm_map
    else:
        x_train[nn] = sipm_map



# In[ ]:

# construct a 5 layer network
def lay5():
    model = Sequential()
    model.add(Dense(output_dim=2048, input_dim=nsipm*nsipm))
    model.add(Activation("relu"))
    model.add(Dense(output_dim=1024))
    model.add(Activation("relu"))
    model.add(Dense(output_dim=512))
    model.add(Activation("relu"))
    model.add(Dense(output_dim=256))
    model.add(Activation("relu"))
    #model.add(Dropout(.25))
    model.add(Dense(output_dim=2*N_ELpts))
    model.add(Activation("sigmoid"))
    model.compile(loss='mse', optimizer=SGD(lr=1.0, momentum=0.9, nesterov=True))
    model.summary()
    N_layers = 5
    return model,N_layers

# construct a 4 layer network
def lay4():
    model = Sequential()
    model.add(Dense(output_dim=2048, input_dim=nsipm*nsipm))
    model.add(Activation("relu"))
    model.add(Dense(output_dim=1024))
    model.add(Activation("relu"))
    model.add(Dense(output_dim=512))
    model.add(Activation("relu"))
    #model.add(Dropout(.25))
    model.add(Dense(output_dim=2*N_ELpts))
    model.add(Activation("sigmoid"))
    model.compile(loss='mse', optimizer=SGD(lr=1.0, momentum=0.9, nesterov=True))
    model.summary()
    N_layers = 4
    return model,N_layers

# etc ..
def lay3():
    model = Sequential()
    model.add(Dense(output_dim=2048, input_dim=nsipm*nsipm))
    model.add(Activation("relu"))
    model.add(Dense(output_dim=1024))
    model.add(Activation("relu"))
    #model.add(Dropout(.25))
    model.add(Dense(output_dim=2*N_ELpts))
    model.add(Activation("sigmoid"))
    model.compile(loss='mse', optimizer=SGD(lr=1.0, momentum=0.9, nesterov=True))
    model.summary()
    N_layers = 3
    return model,N_layers

def lay2():
    model = Sequential()
    model.add(Dense(output_dim=2048, input_dim=nsipm*nsipm))
    model.add(Activation("relu"))
    #model.add(Dropout(.25))
    model.add(Dense(output_dim=2*N_ELpts))
    model.add(Activation("sigmoid"))
    model.compile(loss='mse', optimizer=SGD(lr=1.0, momentum=0.9, nesterov=True))
    model.summary()
    N_layers = 2
    return model,N_layers

def lay1():
    model = Sequential()
    model.add(Dense(output_dim=2*N_ELpts, input_dim=nsipm**2))
    model.add(Activation("sigmoid"))
    model.compile(loss='mse', optimizer=SGD(lr=1.0, momentum=0.9, nesterov=True))
    model.summary()
    N_layers = 1
    return model,N_layers


# Call layered a network and begin training

# In[ ]:

# call desired network
model,N_layers = lay1()

# stop early and save best model
if ES and Save:
    callbacks = [
        callbacks.EarlyStopping(monitor='val_loss', patience=3, mode='min'), # stop training if val_loss
                                                                             # stops decreasing for 3 epochs

        callbacks.ModelCheckpoint(ES_filepath, monitor='val_loss', save_best_only=True, mode='min')] # save best model
# stop early do not save
elif ES:
    callbacks =[callbacks.EarlyStopping(monitor='val_loss', patience=3, mode='min')]

# do not stop early, decide whether to save later
else: callbacks =[]

hist = model.fit(x_train, y_train, nb_epoch=20, batch_size=100,
                 validation_data=(x_test,y_test),
                 verbose=2, callbacks=callbacks);


# Plot training and validation error over epochs

# In[ ]:

val   = plt.plot(range(0,len(hist.history['loss'])),
                 np.sqrt(np.array(hist.history['val_loss'])/N_ELpts)*max_xy, 'b', label='Val error')
train = plt.plot(range(0,len(hist.history['loss'])),
                 np.sqrt(np.array(hist.history['loss'])    /N_ELpts)*max_xy, 'r', label='Train error')
plt.legend(loc='upper right')
plt.xlabel('Epochs')
plt.ylabel('Mean error per hit (mm)')
plt.title(str(N_ELpts) + ' hit(s), ' + str(N_layers) +' layers')
plt.show()


# Plot histogram of error

# In[ ]:

# Get model's predictions, plot the error in a histogram
predictions = model.predict(x_test)
st = (predictions - y_test)**2
E = np.empty((Nevts_test,N_ELpts))
for pt in range(N_ELpts):
    E[:,pt] = np.sqrt(st[:,2*pt] + st[:,2*pt+1])*80

# plot the histogram, of average error per point in distance
n, bins, patches = plt.hist(E.flatten(), 50, color='red', normed=1, alpha=0.75)
plt.xlim(0,30)
plt.xlabel('Error/Hit (mm)')
plt.grid(True)
plt.title(str(N_ELpts) + ' hit(s), ' + str(N_layers) +' layers')
plt.show()


# Plot error as a function of the distance from 40,40 of the EL hit

# In[ ]:

from collections import defaultdict
dd=defaultdict(int) # defaultdict whose keys are distance from 40,40, whose values are error
dn=defaultdict(int) # defaultdict whose keys are distance from 40,40, whose values are the number of evts
                    # with that distance

for ev in range(Nevts_test):
    for pt in range(N_ELpts):
        dist = np.sqrt((y_test[ev,2*pt]*xlen-40)**2 + (y_test[ev,2*pt+1]*ylen-40)**2)
        dd[dist] += E[ev,pt]
        dn[dist] += 1

for key in dd: dd[key] /= dn[key]

# this line is pretty meaningless since I am picking the degree polynomial --------->
bestfit = plt.plot(np.unique(dd.keys()), np.poly1d(np.polyfit(dd.keys(),dd.values(), 3))(np.unique(dd.keys())))

# bin error for plotting
bins = defaultdict(int)
no   = defaultdict(int)
for key in dd:
    bins[round(key)] += dd[key]
    no[round(key)]   += 1
for key in bins: bins[key] /= no[key]
plt.scatter(bins.keys(),bins.values())

plt.xlabel('Radial distance of EL hit from center of EL plane')
plt.ylabel('Average error per hit (mm)')
plt.title(str(N_ELpts) + ' hit(s), ' + str(N_layers) +' layers')
plt.grid(True)
plt.show()



# Plot a training event

# In[ ]:

# Plot one event.
plt_train = False
pevt = np.random.randint(0,Nevts_test)
fig = plt.figure();
ax1 = fig.add_subplot(111);
ax1.axis([0, xlen, 0, ylen]);

if(plt_train):
    xarr = x_train[pevt]
    yarr = y_train[pevt]
else:
    xarr = x_test[pevt]
    yarr = y_test[pevt]

# Create circles and plot them according to the probabilities.
probs = (xarr - min(xarr))
probs /= max(probs)
for x,y,p in zip(pos_x, pos_y, probs):

    #print "Placing sipm at ({0},{1}) with prob {2}".format(x,y,p);

    # Set up the location; note we must invert y due to a pi rotation
    #  about the x-axis.
    r = Ellipse(xy=(x,y), width=2., height=2.);
    r.set_facecolor('0');
    r.set_alpha(0.02 + 0.98*p);
    ax1.add_artist(r);

# Place large blue circles for actual EL points.
for npt in range(len(yarr)/2):
    xpt = yarr[2*npt]*xlen
    ypt = yarr[2*npt+1]*ylen
    mrk = Ellipse(xy=(xpt,ypt), width=4., height=4.);
    mrk.set_facecolor('b');
    ax1.add_artist(mrk);

# Place small red circles for predicted EL points.
for npt in range(len(predictions[pevt])/2):
    xpt = predictions[pevt][2*npt]*xlen
    ypt = predictions[pevt][2*npt+1]*ylen
    mrk = Ellipse(xy=(xpt,ypt), width=2., height=2.);
    mrk.set_facecolor('r');
    ax1.add_artist(mrk);

# Show the plot.
plt.xlabel("x (mm)");
plt.ylabel("y (mm)");
plt.show()


# Compare with other machine learning tools:

# In[ ]:

from sklearn.linear_model import LinearRegression
linreg = LinearRegression()
linreg_model = linreg.fit(x_train, y_train)
print('----training now----')
print('LinReg Error in mm/hit: '  + str(
        np.sqrt(np.sum((linreg_model.predict(x_test) - y_test)**2/N_ELpts/y_test.shape[0])) * 80) )


# Save the model if not already saved

# In[ ]:

if not ES and Save:
    model_filename = str(N_ELpts) + 'hit' + str(N_layers) +'lay.h'
    model.save(model_filename)


# Also save the history for DNNcomparison

# In[ ]:

hist_filename = 'hist_' + str(N_ELpts) + 'hit.h'

if Save and not ES: # Can get rid of 'and not ES', but at the moment i dont want to keep the history
                    # for comparison if each of the histories will have different numbers of epochs

    if MOD: f = tb.open_file(hist_filename, 'r+')
    else:   f = tb.open_file(hist_filename, 'w' )

    filters = tb.Filters(complib='blosc', complevel=9, shuffle=False) # define tables filters
    arrayname = str(N_layers) + 'lay'                                 # define new array name, (could be carray)

    # put history in ndarray
    val_err_hist = np.sqrt(np.array(hist.history['val_loss'])/N_ELpts)*max_xy

    # put ndarray in tables earray
    atom = tb.Atom.from_dtype(val_err_hist.dtype)
    err_hist = f.create_earray(f.root,arrayname, atom, (0,val_err_hist.shape[0]), filters=filters)
    err_hist.append([val_err_hist])



# In[ ]:

print(f)
f.close()
