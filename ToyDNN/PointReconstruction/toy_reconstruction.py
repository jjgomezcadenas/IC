
# # Main script for toy reconstruction with Keras
from sklearn.linear_model       import LinearRegression
from toy_inputs                 import *
import matplotlib.pyplot as plt
import numpy             as np
import tables            as tb
import toy_mc
import toy_nets
import toy_plots

# Now load or generate the hits and SiPM response data
# ** All of the parameters listed above must be the same as those used to create the data in f.root.group ** #
# ** This can be changed if needed so that the above parameters are also saved/loaded in the pytable, but ** #
# ** for the moment this seems unnecessary.                                                               ** #
if Load_MC:
    x_train,y_train,x_valid,y_valid,Nevts_train,Nevts_valid = toy_mc.load_data(MC_filename)
    Nevts = Nevts_train + Nevts_valid
else:
    Nevts = Nevts_train + Nevts_valid
    x_train,y_train,x_valid,y_valid = toy_mc.generate_data()
if Save_MC: toy_mc.save_data(MC_filename,x_train,y_train,x_valid,y_valid)


# Call a network and begin training
model,N_layers = toy_nets.cnn()

# stop early and save best model
if ES and Save_model:
    callbacks = [
        callbacks.EarlyStopping(monitor='val_loss',      # stop training if val_loss
                                patience=5, mode='min'), # stops decreasing for 5 epochs
        callbacks.ModelCheckpoint(ES_filepath,           # save best model
                                monitor='val_loss',
                                save_best_only=True,
                                mode='min')]
# stop early do not save
elif ES: callbacks =[callbacks.EarlyStopping(monitor='val_loss', patience=5, mode='min')]

# do not stop early, decide whether to save later *use this when comparing error of different numbers of layers
else: callbacks =[]

# finally, train the network
if N_layers == 'cnn':
    hist = model.fit(x_train.reshape((Nevts_train,1,nsipm,nsipm)), y_train, nb_epoch=20, batch_size=100,
                 validation_data=(x_valid.reshape((Nevts_valid,1,nsipm,nsipm)),y_valid),
                 verbose=2, callbacks=callbacks);
else:
    hist = model.fit(x_train, y_train, nb_epoch=20, batch_size=100,
                 validation_data=(x_valid,y_valid),
                 verbose=2, callbacks=callbacks);


# Study the network a little with some plots:
#
# Plot of training and validation error over epochs
toy_plots.plot_net_history(hist.history)

# Get the trained model's predictions for the cv dataset
if N_layers == 'cnn':
    predictions = model.predict(x_valid.reshape((Nevts_valid,1,nsipm,nsipm)))
else:
    predictions = model.predict(x_valid)

# Error in x and in y as a function of the x and y coordinate of the hit
toy_plots.plot_xy_error(predictions,y_valid)

# Compute sqauared error of predictions
se = (predictions - y_valid)**2
E = np.empty((Nevts_valid,N_ELpts))
for pt in range(N_ELpts):
    # compute error in distance of predictions
    E[:,pt] = np.sqrt(se[:,2*pt] + se[:,2*pt+1])*max_xy

# plot the histogram, of average error per point in distance
toy_plots.plot_histogram(2000,3)

# Plot error as a function of the distance of the EL hit from the center of the EL plane (40,40)
toy_plots.plot_radial_error()


# Plot a sipm map with a label and prediction
toy_plots.plot_event()

# Compare accuracy with linear regression (sklearn)
def linreg(xtrain,ytrain,xvalid,yvalid,N_ELpts,max_xy):
    linreg = LinearRegression()
    linreg_model = linreg.fit(xtrain, ytrain)
    print('----training now----')
    print('LinReg Error in mm/hit: '  + str(
            np.sqrt(np.sum((linreg_model.predict(xvalid) - yvalid)**2/N_ELpts/yvalid.shape[0])) * max_xy) )
linreg(x_train,y_train,x_valid,y_valid,N_ELpts,max_xy)


# Save the model if not already saved:
# Note if ES then the model should already be saved!
if not ES and Save_model:
    model_filename = str(N_ELpts) + 'hit' + str(N_layers) +'lay.h'
    model.save(model_filename)

# Save the model's training/valid error history for network comparison
def save_history(hist_filename,MOD):
    if MOD: f = tb.open_file(hist_filename, 'r+')
    else:   f = tb.open_file(hist_filename, 'w' )

    filters = tb.Filters(complib='blosc', complevel=9, shuffle=False) # define tables filters
    if N_layers == 'cnn': arrayname = 'cnn'
    else: arrayname = str(N_layers) + 'lay'                           # define new array name, (could be carray)

    # put history in ndarray
    val_err_hist = np.sqrt(np.array(hist.history['val_loss'])/N_ELpts)*max_xy

    # put ndarray in tables earray
    atom = tb.Atom.from_dtype(val_err_hist.dtype)
    err_hist = f.create_earray(f.root,arrayname, atom, (0,val_err_hist.shape[0]), filters=filters)
    err_hist.append([val_err_hist])
    print(f)
    f.close()

if Save_hist: save_history('hist_' + str(N_ELpts) + 'hit.h', MOD_hist)
