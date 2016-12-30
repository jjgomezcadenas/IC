from __future__  import print_function

from matplotlib.patches         import Ellipse
import matplotlib.pyplot as plt
import numpy  as np
import random as rd
import tables as tb
import keras.backend.tensorflow_backend as K

from scipy.stats import threshold

from keras.models               import Model, load_model
from keras.layers               import Input, Dense, MaxPooling3D, Convolution3D, Activation, Dropout
from keras.layers.normalization import BatchNormalization
from keras.optimizers           import SGD, Adam, Nadam         
#from keras.callbacks            import ReduceLROnPlateau !!!!!
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.layers.core          import Flatten
from keras                      import callbacks
from keras.regularizers         import l2, activity_l2

# Load the data.
Ntrain = 48000     # number of training events per sample
Ntot = 50000
lmodel = False
Nepochs = 2000

# Signal events.
s_dat = tb.open_file('/home/jrenner/data/classification/NEW_training_MC_si_nst_nonorm.h5', 'r')
print(s_dat)
s_array = np.array(s_dat.root.maps)
x_t = s_array[:Ntrain]
x_v = s_array[Ntrain:Ntot]
y_t = np.ones([Ntrain, 1])
y_v = np.ones([Ntot-Ntrain, 1])

s_earray = np.array(s_dat.root.energies)

# Background events.
b_dat = tb.open_file('/home/jrenner/data/classification/NEW_training_MC_bg_nst_nonorm.h5', 'r')
print(b_dat)
b_array = np.array(b_dat.root.maps)
print("Concatenating datasets...")
x_t = np.concatenate([x_t, b_array[:Ntrain]])
x_v = np.concatenate([x_v, b_array[Ntrain:Ntot]])
y_bt = np.zeros([Ntrain, 1])
y_t = np.concatenate([y_t, y_bt])
y_bv = np.zeros([Ntot-Ntrain, 1])
y_v = np.concatenate([y_v, y_bv])

b_earray = np.array(b_dat.root.energies)

# Normalize
#mval = max(np.max(s_array),np.max(b_array))
#print("Normalizing with max value of", mval)
#x_t /= mval
#x_v /= mval

# Include the final dimension (single-channel).
print("Reshaping...")
#x_t = np.expand_dims(x_t, axis=1)
#x_v = np.expand_dims(x_v, axis=1)
x_t = np.reshape(x_t, (len(x_t), 48, 48, 30, 1))
x_v = np.reshape(x_v, (len(x_v), 48, 48, 30, 1))
print("Prepared", len(x_t), "training events and", len(x_v), "validation events.")

# Remove noise.
nthr = 5.0e-5
inoise = x_t < nthr
x_t[inoise] = 0

inoise = x_v < nthr
x_v[inoise] = 0

# Set up a DNN for 3D convolution
with K.tf.device('/gpu:1'):
    
    if(lmodel):
        incep = load_model('models/smallnet.h5')
    else:
        K.set_session(K.tf.Session(config=K.tf.ConfigProto(allow_soft_placement=True, log_device_placement=True)))

        inputs = Input(shape=(48, 48, 30, 1))
        cinputs = Convolution3D(64, 6, 6, 6, border_mode='same', subsample=(3, 3, 5), activation='relu',init='lecun_uniform', W_regularizer=l2(0.0001))(inputs)
        cinputs = MaxPooling3D(pool_size=(3, 3, 3), strides=(2, 2, 3), border_mode='same', dim_ordering='default')(cinputs)
        cinputs = BatchNormalization(epsilon=1e-05, mode=0, axis=4, momentum=0.99, weights=None, beta_init='zero', gamma_init='one', gamma_regularizer=None, beta_regularizer=None)(cinputs)
        cinputs = Convolution3D(16, 1, 1, 1, border_mode='same', subsample=(1, 1, 1), activation='relu',init='lecun_uniform', W_regularizer=l2(0.0001))(cinputs)
        cinputs = Convolution3D(32, 2, 2, 2, border_mode='same', subsample=(1, 1, 1), activation='relu',init='lecun_uniform', W_regularizer=l2(0.0001))(cinputs)
        cinputs = BatchNormalization(epsilon=1e-05, mode=0, axis=4, momentum=0.99, weights=None, beta_init='zero', gamma_init='one', gamma_regularizer=None, beta_regularizer=None)(cinputs)
        cinputs = MaxPooling3D(pool_size=(2, 2, 2), strides=(2, 2, 2), border_mode='same', dim_ordering='default')(cinputs)
        f1 = Flatten()(cinputs)
        f1 = Dense(output_dim=32, activation='relu', init='lecun_uniform', W_regularizer=l2(0.0001))(f1)
        f1 = Dropout(.6)(f1)

        inc_output = Dense(output_dim=1, activation='sigmoid',init='normal', W_regularizer=l2(0.0002))(f1)
        incep = Model(inputs, inc_output)

        incep.compile(loss='binary_crossentropy', optimizer=Nadam(lr=0.0005, beta_1=0.9, beta_2=0.999, epsilon=1e-08, schedule_decay=0.005), metrics=['accuracy'])
        lcallbacks = [callbacks.ModelCheckpoint('models/conv3d_classifier.h', monitor='val_loss', save_best_only=True, mode='min')]

    incep.summary()
 
# Train.
hist = incep.fit(x_t, y_t, shuffle=True, nb_epoch=Nepochs, batch_size=100, verbose=1, validation_data=(x_v, y_v), callbacks=lcallbacks)

# Save the model
incep.save('models/nonoise_long.h5')
