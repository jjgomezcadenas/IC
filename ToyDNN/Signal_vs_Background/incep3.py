
# Sg_vs_Bg_dnn will use magbox SiPM response data to classify tracks into signal or background classes
from   __future__         import print_function
from   keras.optimizers   import SGD, Nadam
from   keras.models       import Sequential, Model
from   keras.layers       import Input, Dense, Activation, Convolution2D, AveragePooling2D, MaxPooling2D, merge, Reshape, Flatten, Dropout
from   keras              import callbacks
from   matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import tables            as tb
import numpy             as np
import copy
import h5py
import sys

load = True # if load, will load saved prepared and increase speed
inception = True
one_hot_encode = False # Do not set to true if inception

# Set max number of slices allowed
# Pad and prepare data for dnn
def create_clsf_dataset(num_slices,fs,fb):
    """
    create_clsf_dataset creates x and y_ for the big classfier dnn. x and y_
    can be for training, validating, or testing

    num_slices is the number of slices that will be given to each event.
    events with more slices are discarded, events with fewer are padded with zeros

    fs is the file with signal data, sb is the file with background data.
    """
    temp_num_evts = len(fs['nspevt'])
    num_sipms = fs['sipm0'].shape[0]
    x = []; y_ = []

    # keep track of slice numbers
    sg_snum=0; bg_snum=0
    print('---Generating ' + str(temp_num_evts) + ' events')
    for i in range(temp_num_evts):

        if (i+1) % (temp_num_evts/10) == 0 :
            print('---Event: ' + str(i+1) + '/' + str(temp_num_evts))

        # ---Signal
        slice_len = fs['nspevt'][i]

        # discard events with more than max_slices
        if slice_len <= num_slices:
            y_.append(np.array([1])) # signal

            # distribute 0-padded slices randomly before and after track
            if num_slices - slice_len == 0: rand = 0
            else: rand = np.random.randint(num_slices - slice_len)
            xevt = np.zeros((num_slices, num_sipms),dtype=np.float32)
            for s in range(slice_len): xevt[rand + s] = fs['sipm{0}'.format(s + sg_snum)]
            x.append(xevt)

        sg_snum += slice_len

        # ---Background
        slice_len = fb['nspevt'][i]

        # discard events with more than max_slices
        if slice_len <= num_slices:
            y_.append(np.array([0])) # background

            # distribute 0-padded slices randomly before and after track
            if num_slices - slice_len == 0: rand = 0
            else: rand = np.random.randint(num_slices - slice_len)
            xevt = np.zeros((num_slices, num_sipms),dtype=np.float32)
            for s in range(slice_len): xevt[rand + s] = fb['sipm{0}'.format(s + sg_snum)]
            x.append(xevt)

        bg_snum += slice_len
    print('---finished.')
    return(np.array(x),np.array(y_))

if load:
    print('---Loading Data')
    dnn_dat = tb.open_file('dnn_data/prepared_events.h', 'r'); print(dnn_dat)
    x_tv    = np.array(dnn_dat.root.x)
    y_tv    = np.array(dnn_dat.root.y)
    nslices = x_tv.shape[1]
    dnn_dat.close()
    print('---Load Complete.')


else:
    # Load and pad SiPM response data, so that each event has same number of time slices!
    sg = h5py.File('dnn_data/signal_2x2x2.h5',    'r')
    bg = h5py.File('dnn_data/background_2x2x2.h5','r')

    nslices = 30
    x_tv,y_tv = create_clsf_dataset(nslices,sg,bg)
    dnn_dat   = tb.open_file('dnn_data/prepared_events.h', 'w')
    filters   = tb.Filters(complib='blosc', complevel=9, shuffle=False)
    atomx     = tb.Atom.from_dtype(x_tv.dtype)
    atomy     = tb.Atom.from_dtype(y_tv.dtype)
    x         = dnn_dat.create_earray(dnn_dat.root, 'x', atomx,
                                     (0,nslices,x_tv.shape[2]),
                                     filters=filters, expectedrows=len(x_tv))
    y         = dnn_dat.create_earray(dnn_dat.root, 'y', atomy,
                                     (0,1),
                                     filters=filters, expectedrows=len(y_tv))

    print('---Saving Data')
    for ev in range(len(y_tv)):
        x.append([x_tv[ev]])
        y.append([y_tv[ev]])
    dnn_dat.close()
    print('---Save Complete.')

    sg.close()
    bg.close()

if one_hot_encode:
    # make labels one hot encoded
    yh_tv = np.zeros((y_tv.shape[0],2))
    ev=0
    for i in y_tv:
        yh_tv[ev,i] = 1
        ev+=1

# Set up a DNN (no inception)
if not inception:
    nsipm = int(np.sqrt(x_tv.shape[2]))

    # Construct a DNN.
    model = Sequential()

    # Reshape the input
    model.add(Reshape((nslices, nsipm, nsipm), input_shape=(nslices, nsipm**2)))
    model.add(Flatten())
    model.add(Dense(output_dim=4096))
    model.add(Activation("relu"))
    model.add(Dense(output_dim=1024))
    model.add(Activation("relu"))
    model.add(Dense(output_dim=256))
    model.add(Activation("relu"))
    model.add(Dropout(.2))

    if one_hot_encode:
        print('one hot encoded output')
        model.add(Dense(output_dim=2))
        model.add(Activation("softmax"))
        # Setup
        #model.compile(loss='categorical_crossentropy',
        #              optimizer=Nadam(lr=0.005, beta_1=0.9, beta_2=0.999, epsilon=1e-02, schedule_decay=0.004),
        #              metrics=['accuracy'])
        model.compile(loss='categorical_crossentropy',
                      optimizer=SGD(lr=0.01, momentum=0.9, nesterov=True),
                      metrics=['accuracy'])

    else:
        model.add(Dense(output_dim=1))
        model.add(Activation("sigmoid"))

        # Setup
        #model.compile(loss='binary_crossentropy',
        #              optimizer=Nadam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-02, schedule_decay=0.004),
        #              metrics=['accuracy'])
        model.compile(loss='binary_crossentropy',
                      optimizer=SGD(lr=0.005, momentum=0.9, nesterov=True),
                      metrics=['accuracy'])
    model.summary()


# Set up a dnn with inception modules
if inception:
    nsipm  = int(np.sqrt(x_tv.shape[2]))
    div    = 1
    inputs = Input(shape=(nslices,nsipm,nsipm))

    #cinputs = Convolution2D(30  /div, 1, 1, border_mode='same', subsample=(1,1), activation='relu')(inputs)
    #cinputs = Convolution2D(512 /div, 3, 3, border_mode='valid', subsample=(1,1), activation='relu')(inputs)

    # inception 1
    c1  = Convolution2D(64 /div, 1, 1, border_mode='same', subsample=(1,1), activation='relu')(inputs)
    c3r = Convolution2D(96 /div, 1, 1, border_mode='same', subsample=(1,1), activation='relu')(inputs)
    c3  = Convolution2D(128/div, 3, 3, border_mode='same', subsample=(1,1), activation='relu')(c3r)
    c5r = Convolution2D(16 /div, 1, 1, border_mode='same', subsample=(1,1), activation='relu')(inputs)
    c5  = Convolution2D(32 /div, 5, 5, border_mode='same', subsample=(1,1), activation='relu')(c5r)
    mp  = MaxPooling2D (pool_size=(2, 2), strides=(1,1), border_mode='same', dim_ordering='default')(inputs)
    mpr = Convolution2D(32/ div, 1, 1, border_mode='same', subsample=(1,1), activation='relu')(mp)
    m1  = merge([c1,c3,c5,mpr], mode='concat', concat_axis=1)

    #m1 = MaxPooling2D (pool_size=(2, 2), strides=(2,2), border_mode='same', dim_ordering='default')(m1)

    # inception 2
    c12  = Convolution2D(128/div, 1, 1, border_mode='same', subsample=(1,1), activation='relu')(m1)
    c32r = Convolution2D(128/div, 1, 1, border_mode='same', subsample=(1,1), activation='relu')(m1)
    c32  = Convolution2D(192/div, 3, 3, border_mode='same', subsample=(1,1), activation='relu')(c32r)
    c52r = Convolution2D(32 /div, 1, 1, border_mode='same', subsample=(1,1), activation='relu')(m1)
    c52  = Convolution2D(96 /div, 5, 5, border_mode='same', subsample=(1,1), activation='relu')(c52r)
    mp2  = MaxPooling2D (pool_size=(2, 2), strides=(1,1), border_mode='same', dim_ordering='default')(m1)
    mp2r = Convolution2D(64 /div, 1, 1, border_mode='same', subsample=(1,1), activation='relu')(mp2)
    m2   = merge([c12,c32,c52,mp2r], mode='concat', concat_axis=1)


    m2 = MaxPooling2D (pool_size=(2, 2), strides=(2,2), border_mode='same', dim_ordering='default')(m2)

    # inception 3
    c13  = Convolution2D(192/div, 1, 1, border_mode='same',subsample=(1,1), activation='relu')(m2)
    c33r = Convolution2D(96 /div, 1, 1, border_mode='same',subsample=(1,1), activation='relu')(m2)
    c33  = Convolution2D(208/div, 3, 3, border_mode='same',subsample=(1,1), activation='relu')(c33r)
    c53r = Convolution2D(16 /div, 1, 1, border_mode='same',subsample=(1,1), activation='relu')(m2)
    c53  = Convolution2D(48 /div, 5, 5, border_mode='same',subsample=(1,1), activation='relu')(c53r)
    mp3  = MaxPooling2D (pool_size=(2, 2), strides=(1,1), border_mode='same', dim_ordering='default')(m2)
    mp3r = Convolution2D(64 /div, 1, 1, border_mode='same',subsample=(1,1), activation='relu')(mp3)
    m3   = merge([c13,c33,c53,mp3r], mode='concat', concat_axis=1)

    # dense
    #d1 = Dense(output_dim=2048, activation='relu')(fl)
    a1 = AveragePooling2D(pool_size=(2,2), strides=(2,2), border_mode='same')(m3)
    fc1 = Convolution2D(128 /div, 1, 1, border_mode='same',subsample=(1,1), activation='relu')(a1)
    fl = Flatten()(fc1)
    d2 = Dense(output_dim=1024 , activation='relu')(fl)
    #d2 = Dropout(.7)(d2)

    inc_output = Dense(output_dim=1,activation='sigmoid')(d2)
    incep = Model(inputs,inc_output)

    #incep.compile(loss='binary_crossentropy',
    #              optimizer=SGD(lr=0.001, momentum=0.9, nesterov=True),
    #              metrics=['accuracy'])

    #incep.compile(loss='binary_crossentropy',
    #                  optimizer=Nadam(lr=0.0001, beta_1=0.9, beta_2=0.999, epsilon=1e-02, schedule_decay=0.004),
    #                  metrics=['accuracy'])

    incep.compile(loss='binary_crossentropy',
                  optimizer=Nadam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, schedule_decay=0.004),
                  metrics=['accuracy'])

    incep.summary()

# Train the network
Ntrain    = 36000
callbacks = [callbacks.EarlyStopping(monitor='val_acc', patience=20, mode='max'),
             callbacks.ModelCheckpoint('models/real_inception.h', monitor='val_acc', save_best_only=True, mode='max'),
             #callbacks.ReduceLROnPlateau(monitor='val_acc', mode='max', factor=0.5, patience=5)
            ]

if inception:
    # Reshape the network for training
    x_tv = np.reshape(x_tv,(len(x_tv),nslices,nsipm,nsipm))
    hist=incep.fit(x_tv[:Ntrain], y_tv[:Ntrain],
                   shuffle=True, nb_epoch=40, batch_size=50, verbose=1, callbacks=callbacks,
                   validation_data=(x_tv[Ntrain:],y_tv[Ntrain:]))
elif one_hot_encode:
    hist=model.fit(x_tv[:Ntrain], yh_tv[:Ntrain], shuffle=True, nb_epoch=10, batch_size=50, validation_data=(x_tv[Ntrain:], yh_tv[Ntrain:]))
else:
    hist=model.fit(x_tv[:Ntrain], y_tv [:Ntrain], shuffle=True, nb_epoch=10, batch_size=50, validation_data=(x_tv[Ntrain:], y_tv [Ntrain:]),verbose=2)

# Save History
f = tb.open_file('histories/fixed2.h', 'w')
filters = tb.Filters(complib='blosc', complevel=9, shuffle=False) # define tables filters
val_err_hist = np.array(hist.history['val_loss'])
atom = tb.Atom.from_dtype(val_err_hist.dtype)
err_hist = f.create_earray(f.root,'real_inception', atom, (0,val_err_hist.shape[0]), filters=filters)
err_hist.append([val_err_hist])
f.close()
