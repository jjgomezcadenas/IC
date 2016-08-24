"""

Configuration for different nets.

"""

import h5py
import numpy as np
import tensorflow as tf

from inputs import *
# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Helper methods.
# -----------------------------------------------------------------------------
def weight_variable(shape):
  initial = tf.truncated_normal(shape, stddev=0.01)
  return tf.Variable(initial)

def bias_variable(shape):
  initial = tf.constant(0.01, shape=shape)
  return tf.Variable(initial)

def conv2d(x, W):
  return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

def conv3d(x, W):
  return tf.nn.conv3d(x, W, strides=[1, 1, 1, 1, 1], padding='SAME')

def max_pool_2x2x2(x):
  return tf.nn.max_pool3d(x, ksize=[1, 2, 2, 2, 1],
                           strides=[1, 2, 2, 2, 1], padding='SAME')

def max_pool_2x2(x):
  return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],
                        strides=[1, 2, 2, 1], padding='SAME')

def max_pool_nrs(x):
  return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],
                        strides=[1, 1, 1, 1], padding='SAME')


def avg_pool_2x2(x):
  return tf.nn.avg_pool(x, ksize=[1, 2, 2, 1],
                        strides=[1, 2, 2, 1], padding='SAME')

# -----------------------------------------------------------------------------
# Set up the neural network.
# -----------------------------------------------------------------------------
def getNet(name,x_input):

    if(name == "MNISTbasic"):
        return MNISTbasic_net(x_input)
    if(name == "MNISTadv"):
        return MNISTadv_net(x_input)
    if(name == "MNISTadv_nodrop"):
        return MNISTadv_nodrop_net(x_input)


# -----------------------------------------------------------------------------
# DNN definitions
# -----------------------------------------------------------------------------
def MNISTbasic_net(x_input):

    # Softmax readout: 2 classifications
    W = weight_variable([nsensors, ngrid])
    b = bias_variable([ngrid])
    y = tf.nn.softmax(tf.matmul(x_input, W) + b)

    return y


# Function to write for part 2.
def MNISTadv_net(x_input):

    # Resize the array to pdim x pdim x 3
    x_image = tf.reshape(x_input, [-1,nsipm,nsipm,1])

    # First convolutional layer
    W_conv1 = weight_variable([4, 4, 1, 32])
    b_conv1 = bias_variable([32])
    h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)
    h_pool1 = max_pool_2x2(h_conv1)

    # Second convolutional layer
    W_conv2 = weight_variable([2, 2, 32, 64])
    b_conv2 = bias_variable([64])
    h_conv2 = tf.nn.relu(conv2d(h_pool1, W_conv2) + b_conv2)
    h_pool2 = max_pool_2x2(h_conv2)

    # Densely connected layer
    new_dim = int(nsipm / 4)
    W_fc1 = weight_variable([new_dim * new_dim * 64, 1024])
    b_fc1 = bias_variable([1024])
    h_pool2_flat = tf.reshape(h_pool2, [-1, new_dim*new_dim*64])
    h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1)

    # Dropout
    keep_prob = 0.4  #tf.placeholder("float")
    h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

    # Softmax readout
    W_fc2 = weight_variable([1024, ngrid])
    b_fc2 = bias_variable([ngrid])
    # y = tf.nn.softmax(tf.matmul(h_fc1_drop, W_fc2) + b_fc2)
    y = tf.matmul(h_fc1_drop, W_fc2) + b_fc2
    return y

# Function to write for part 2.
def MNISTadv_nodrop_net(x_input):

    # Resize the array to pdim x pdim x 3
    x_image = tf.reshape(x_input, [-1,nsipm,nsipm,1])

    # First convolutional layer
    W_conv1 = weight_variable([4, 4, 1, 32])
    b_conv1 = bias_variable([32])
    h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)
    h_pool1 = max_pool_2x2(h_conv1)

    # Second convolutional layer
    W_conv2 = weight_variable([2, 2, 32, 64])
    b_conv2 = bias_variable([64])
    h_conv2 = tf.nn.relu(conv2d(h_pool1, W_conv2) + b_conv2)
    h_pool2 = max_pool_2x2(h_conv2)

    # Densely connected layer
    new_dim = int(nsipm / 4)
    W_fc1 = weight_variable([new_dim * new_dim * 64, 1024])
    b_fc1 = bias_variable([1024])
    h_pool2_flat = tf.reshape(h_pool2, [-1, new_dim*new_dim*64])
    h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1)

    # Dropout
    keep_prob = 1.0  #tf.placeholder("float")
    h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

    # Softmax readout
    W_fc2 = weight_variable([1024, ngrid])
    b_fc2 = bias_variable([ngrid])
    # y = tf.nn.softmax(tf.matmul(h_fc1_drop, W_fc2) + b_fc2)
    y = tf.matmul(h_fc1_drop, W_fc2) + b_fc2
    return y
