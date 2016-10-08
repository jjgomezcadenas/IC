from keras.models               import Sequential
from keras.layers               import Dense, Activation, Dropout
from keras.optimizers           import SGD, Adam, Nadam
from keras                      import callbacks
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.layers.core          import Flatten

from toy_inputs import *


# Define some sample networks:
# Convolutional
def cnn():
    model = Sequential()

    # **Worth taking into consideration that our image size is tiny (8x8), convolution may work much better for
    # **with 1792 sipms

    # kernal size is 3x3, 32 filters, padding is same.
    # Same padding works better, this is probably because same padding makes it easier for network no to retain as
    # much information as possible around the edges.
    model.add(Convolution2D(32,3,3,border_mode='same',input_shape=(1, nsipm, nsipm)))
    model.add(Activation('relu'))
    model.add(Convolution2D(32,3,3,border_mode='same', input_shape=(32, nsipm, nsipm)))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Flatten())
    model.add(Dense(output_dim=128))
    model.add(Activation('relu'))
    model.add(Dense(output_dim=64))
    model.add(Activation('relu'))
    model.add(Dense(output_dim=2*N_ELpts))
    model.add(Activation('sigmoid'))

    # Nadam optimizer is a safe choice at least for deep networks. It is adam optimizer with Nesterov
    # Momentum. Nesterov Momentum takes into account future expected future gradient gradient, unlike traditional Mom.
    model.compile(loss='mse', optimizer=Nadam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, schedule_decay=0.004))
    N_layers = 'cnn'
    return model,N_layers

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
    model.add(Dropout(.25))    # **Note this is drop rate, not keep rate (like in tf)
    model.add(Dense(output_dim=2*N_ELpts))
    model.add(Activation("sigmoid"))

    # Nadam optimizer is probabaly a safe choice at least for our deep networks. It is adam optimizer with Nesterov
    # Momentum. Nesterov Momentum takes into account future expected future gradient gradient, unlike traditional Mom.
    model.compile(loss='mse', optimizer=Nadam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, schedule_decay=0.004))
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
    model.add(Dropout(.2))
    model.add(Dense(output_dim=2*N_ELpts))
    model.add(Activation("sigmoid"))
    model.compile(loss='mse', optimizer=Nadam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, schedule_decay=0.004))
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
    #model.add(Dropout(.15))
    model.add(Dense(output_dim=2*N_ELpts))
    model.add(Activation("sigmoid"))
    #model.compile(loss='mse', optimizer=SGD(lr=1.0, momentum=0.9, nesterov=True))
    model.compile(loss='mse', optimizer=Nadam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, schedule_decay=0.004))

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
    model.compile(loss ='mse',  optimizer=SGD(lr=1.0, momentum=0.9, nesterov=True))
    #model.compile(loss='mse', optimizer=Nadam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, schedule_decay=0.004))
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
