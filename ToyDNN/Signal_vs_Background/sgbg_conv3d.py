
# sgbg_3dconv will use magbox SiPM response data to classify tracks into
# signal or background classes
from __future__ import print_function
from keras.optimizers import Nadam
from keras.models import Model
from keras.layers import Input, Dense, Convolution3D, Flatten
from keras import callbacks
import tables as tb
import numpy as np
import h5py

load = True  # if load, will load saved prepared and increase speed


# Pad and prepare data for dnn
def create_clsf_dataset(num_slices, fs, fb):
    """
    create_clsf_dataset creates x and y_ for the big classfier dnn. x and y_
    can be for training, validating, or testing

    num_slices is the number of slices that will be given to each event.
    events with more slices are discarded, events with fewer are padded with
    zeros.

    fs is the file with signal data, sb is the file with background data.
    """
    temp_num_evts = len(fs['nspevt'])
    num_sipms = fs['sipm0'].shape[0]
    x = []
    y_ = []

    # keep track of slice numbers
    sg_snum = 0
    bg_snum = 0
    print('---Generating ' + str(temp_num_evts) + ' events')
    for i in range(temp_num_evts):

        if (i + 1) % (temp_num_evts / 10) == 0:
            print('---Event: ' + str(i + 1) + '/' + str(temp_num_evts))

        # ---Signal
        slice_len = int(fs['nspevt'][i])

        # discard events with more than max_slices
        if slice_len <= num_slices:
            y_.append(np.array([1]))  # signal

            # distribute 0-padded slices randomly before and after track
            if num_slices - slice_len == 0:
                rand = 0
            else:
                rand = np.random.randint(num_slices - slice_len)
            xevt = np.zeros((num_slices, num_sipms), dtype=np.float32)
            for s in range(slice_len):
                xevt[rand + s] = fs['sipm{0}'.format(s + sg_snum)]
            x.append(xevt)

        sg_snum += slice_len

        # ---Background
        slice_len = int(fb['nspevt'][i])

        # discard events with more than max_slices
        if slice_len <= num_slices:
            y_.append(np.array([0]))  # background

            # distribute 0-padded slices randomly before and after track
            if num_slices - slice_len == 0:
                rand = 0
            else:
                rand = np.random.randint(num_slices - slice_len)

            xevt = np.zeros((num_slices, num_sipms), dtype=np.float32)
            for s in range(slice_len):
                xevt[rand + s] = fb['sipm{0}'.format(s + sg_snum)]
            x.append(xevt)

        bg_snum += slice_len
    print('---finished.')
    return(np.array(x), np.array(y_))

if load:
    print('---Loading Data')
    dnn_dat = tb.open_file('dnn_data/prepared_events_200k.h', 'r')
    print(dnn_dat)
    x_tv = np.array(dnn_dat.root.x)
    print(sum(sum(x_tv[0])))
    x_tv /= np.max(x_tv)
    print(np.min(x_tv), np.max(x_tv))
    y_tv = np.array(dnn_dat.root.y)
    print('Percent Signal Events: ' + str(float(sum(y_tv)) / float(len(y_tv))))
    nslices = x_tv.shape[1]
    dnn_dat.close()
    print('---Load Complete.')

else:
    # Load and pad SiPM response data,
    # so that each event has same number of time slices!
    sg = h5py.File('dnn_data/signal_2x2x2_200k.h5', 'r')
    bg = h5py.File('dnn_data/background_2x2x2_200k.h5', 'r')

    nslices = 30
    x_tv, y_tv = create_clsf_dataset(nslices, sg, bg)
    x_tv /= np.max(x_tv)
    print(np.min(x_tv), np.max(x_tv))
    dnn_dat = tb.open_file('dnn_data/prepared_events_200k.h', 'w')
    filters = tb.Filters(complib='blosc', complevel=9, shuffle=False)
    atomx = tb.Atom.from_dtype(x_tv.dtype)
    atomy = tb.Atom.from_dtype(y_tv.dtype)
    x = dnn_dat.create_earray(dnn_dat.root, 'x', atomx,
                              (0, nslices, x_tv.shape[2]),
                              filters=filters, expectedrows=len(x_tv))
    y = dnn_dat.create_earray(dnn_dat.root, 'y', atomy, (0, 1),
                              filters=filters, expectedrows=len(y_tv))

    print('---Saving Data')
    for ev in range(len(y_tv)):
        x.append([x_tv[ev]])
        y.append([y_tv[ev]])
    dnn_dat.close()
    print('---Save Complete.')

    sg.close()
    bg.close()

# Set up a DNN for 3D convolution
nsipm = int(np.sqrt(x_tv.shape[2]))
div = 1

inputs = Input(shape=(1, nslices, nsipm, nsipm))
cinputs = Convolution3D(8 / div, 3, 3, 3, border_mode='same',
                        subsample=(1, 1, 1), activation='relu')(inputs)
cinputs = Convolution3D(32 / div, 3, 3, 3, border_mode='valid',
                        subsample=(3, 1, 1), activation='relu')(cinputs)
cinputs = Convolution3D(64 / div, 3, 3, 3, border_mode='valid',
                        subsample=(1, 2, 2), activation='relu')(cinputs)
cinputs = Convolution3D(64 / div, 2, 2, 2, border_mode='valid',
                        subsample=(2, 2, 2), activation='relu')(cinputs)
cinputs = Convolution3D(128 / div, 2, 2, 2, border_mode='valid',
                        subsample=(2, 2, 2), activation='relu')(cinputs)
f1 = Flatten()(cinputs)
f1 = Dense(output_dim=512, activation='relu')(f1)
# f1 = Dropout(.3)(f1)

inc_output = Dense(output_dim=1, activation='sigmoid')(f1)
incep = Model(inputs, inc_output)

incep.compile(loss='binary_crossentropy',
              optimizer=Nadam(lr=0.002, beta_1=0.9, beta_2=0.999,
                              epsilon=1e-08, schedule_decay=0.004),
              metrics=['accuracy'])

incep.summary()

# Train the network
Ntrain = 360000
callbacks = [callbacks.EarlyStopping(monitor='val_loss', patience=20,
                                     mode='min'),
             callbacks.ModelCheckpoint('models/conv3d_200k.h',
                                       monitor='val_loss', save_best_only=True,
                                       mode='min'),
             # callbacks.ReduceLROnPlateau(monitor='val_acc', mode='max',
             # factor=0.5, patience=5)
             ]

# Reshape the network for training
x_tv = np.reshape(x_tv, (len(x_tv), 1, nslices, nsipm, nsipm))
hist = incep.fit(x_tv[:Ntrain], y_tv[:Ntrain],
                 shuffle=True, nb_epoch=40, batch_size=500, verbose=1,
                 callbacks=callbacks,
                 validation_data=(x_tv[Ntrain:], y_tv[Ntrain:]))

# Save History
f = tb.open_file('histories/conv3d_200k.h', 'w')
filters = tb.Filters(complib='blosc', complevel=9, shuffle=False)
val_err_hist = np.array(hist.history['val_loss'])
atom = tb.Atom.from_dtype(val_err_hist.dtype)
err_hist = f.create_earray(f.root, 'conv3d', atom, (0, val_err_hist.shape[0]),
                           filters=filters)
err_hist.append([val_err_hist])
f.close()
