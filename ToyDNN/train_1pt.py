"""
dnntrainreg.py

Trains the DNN analysis for the specified configuration

"""

import h5py
import numpy as np
import tensorflow as tf
import os
import logging
import tables
import sys
import neuralnets
from inputs import *

# Ensure the appropriate directory structure exists.
if(not os.path.isdir(rdir)): os.mkdir(rdir)
if(not os.path.isdir("{0}/{1}".format(rdir,rname))): os.mkdir("{0}/{1}".format(rdir,rname))
if(not os.path.isdir("{0}/{1}/acc".format(rdir,rname))): os.mkdir("{0}/{1}/acc".format(rdir,rname))

# Create the logger object.
if(log_to_file):
    logging.basicConfig(filename="{0}/{1}/{2}.log".format(rdir,rname,rname),format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',level=logging.DEBUG)
else:
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',level=logging.DEBUG)
logging.info("Params:\n ntrain_evts = {0}\n num_epochs = {1}\n epoch_blk_size = {2}\n dtblk_size = {3}\n batch_size = {4}\n nval_evts = {5}".format(ntrain_evts,num_epochs,epoch_blk_size,dtblk_size,batch_size,nval_evts))

# Checks on parameters.
if(ntrain_evts % dtblk_size != 0):
    logging.error("ERROR: ntrain_evts must be evenly divisible by dtblk_size..."); exit()
if(num_epochs % epoch_blk_size != 0):
    logging.error("ERROR: num_epochs must be evenly divisible by epoch_blk_size..."); exit()
if(ntrain_evts % batch_size != 0):
    logging.error("ERROR: ntrain_evts must be evenly divisible by batch_size..."); exit()
if(nval_evts % batch_size != 0):
    logging.error("ERROR: nval_evts must be evenly divisible by batch_size..."); exit()

# Constructed file names.
fn_data = "{0}/{1}.h".format(datdir,dname)
fn_saver = "{0}/{1}/tfmdl_{2}.ckpt".format(rdir,rname,rname)   # for saving trained network
fn_acc = "{0}/{1}/acc/accuracy_{2}.dat".format(rdir,rname,rname)
fn_prob = "{0}/{1}/acc/prob_{2}".format(rdir,rname,rname)

# ---------------------------------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------------------------------

# Evaluate the performance.
def eval_performance(fsummary,epoch,sess,loss,y_out,dat_train,dat_test,lbl_train,lbl_test):

    logging.info(" \n --- Calling eval_performance\n")

    # ----------------------------
    # Evaluate the training data.
    # ----------------------------
    f_distance_train = open("{0}_train_ep{1}.dat".format(fn_prob,epoch),"w")
    avg_dist_tr = 0.; lval_tr = 0.
    nevt = 0; nbatches = 0
    while(nevt < nval_evts):

        # Run the classification.
        ltemp,ytemp = sess.run([loss,y_out],feed_dict={x: dat_train[nevt:nevt+batch_size], y_: lbl_train[nevt:nevt+batch_size]})
        evtno = 0
        for appx in ytemp:

            # Get actual X, Y labels
            exact  = lbl_train[nevt+evtno]

            # Calculate distance between (X,Y) and its approximation
            dist = np.sqrt(((exact[0]-appx[0])*max_xy)**2 + ((exact[1] - appx[1])*max_xy)**2)
            avg_dist_tr += dist

            # Write the probabilities to file.
            f_distance_train.write("{0} {1}".format(exact," ".join(map(str,appx))))

            evtno += 1

        lval_tr += ltemp
        nevt += batch_size; nbatches += 1

    avg_dist_tr /= nval_evts
    f_distance_train.close()


    lval_tr /= nbatches

    # ---------------------------
    # Evaluate the test data.
    # ---------------------------
    f_distance_test = open("{0}_test_ep{1}.dat".format(fn_prob,epoch),"w")
    avg_dist_te = 0; lval_te = 0.
    nevt = 0; nbatches = 0
    while(nevt < nval_evts):

        # Run the regression.
        ltemp,ytemp = sess.run([loss,y_out],feed_dict={x: dat_test[nevt:nevt+batch_size], y_: lbl_test[nevt:nevt+batch_size]})
        evtno = 0
        for appx in ytemp:
            # Get actual X, Y labels
            exact = lbl_train[nevt+evtno]

            # Calculate distance between (X,Y) and its approximation
            dist = np.sqrt(((exact[0]-appx[0])*max_xy)**2 + ((exact[1] - appx[1])*max_xy)**2)
            avg_dist_te += dist

            # Write the probabilities to the file.
            f_distance_test.write("{0} {1}\n".format(exact," ".join(map(str,appx))))

            evtno += 1

        lval_te += ltemp
        nevt += batch_size; nbatches += 1

    avg_dist_te /= nval_evts
    f_distance_test.close()

    lval_te /= nbatches

    # Write to the final summary file.
    fsummary.write("{0} {1} {2} {3} {4}\n".format(epoch,avg_dist_tr,lval_tr,avg_dist_te,lval_te))

# Read in all the data from evt_start to evt_end-1.
def read_data(fdata,evt_start,evt_end):

    logging.info("-- read_data: Reading events from {0} to {1} with {2} sensors and {3} grid points".format(evt_start,evt_end,nsensors,ngrid))

    # Set up the data arrays.
    #dat_sensors = np.zeros([nevts,nsensors]); lbl_sensors = np.zeros([nevts,ngrid])

    # Read in all events from the data file.
    dat_tree = fdata.root.sim_1pt

    #get silicon PM maps and x,y coord labels for these events ##VECTORIZED
    dat_sensors = dat_tree.sipm_resp[0][evt_start:evt_end]
    lbl_sensors = ((np.array([dat_tree.x[0][evt_start:evt_end],dat_tree.y[0][evt_start:evt_end]])).transpose() - 40)/80

    ##VECTORIZED^^
    """

    while(nmap < nevts):
        # Get the elements in this entry.
        #dat_tree.GetEntry(nmap)

        # Fill the data array with the sipm probabilities.
        ss = 0
        for p in dat_tree.sipm_resp[0][nmap]:
            dat_sensors[nmap][ss] = p
            ss += 1


        # Set the label to the active EL point, and normalize it [-0.5, 0.5]
        lbl_sensors[nmap][0] = (dat_tree.x[0][nmap] - (max_xy/2))/max_xy
        lbl_sensors[nmap][1] = (dat_tree.y[0][nmap] - (max_xy/2))/max_xy


        #lbl_sensors[nmap][dat_tree.elpt] = 1

        nmap += 1
    """

    # Return the data and labels.
    return (dat_sensors,lbl_sensors)

# Set up the neural network.
def net_setup():

    logging.info("\n\n-- net_setup():  SETTING UP NETWORK --")

    logging.info("Creating placeholders for input and output variables...")
    x_input = tf.placeholder(tf.float32, [batch_size, nsensors]) # npix])
    y_ = tf.placeholder(tf.float32, [batch_size, ngrid])

    y_out = nets.neuralnets.getNet(net_name,x_input)

    # Set up for training
    logging.info("Setting up tf training variables...")
    #cross_entropy = -tf.reduce_sum(y_*tf.log(y_out + 1.0e-9))

    #loss = tf.reduce_mean(cross_entropy, name='ropy_mean')
    loss = tf.nn.l2_loss(y_ - y_out, name = 'MSE')
    gstep = tf.Variable(0, trainable=False)
    lrate = tf.train.exponential_decay(opt_lr, gstep,
                                           opt_ndecayepochs*batches_per_epoch,
                                           opt_decaybase, staircase=True)
    #train_step = tf.train.MomentumOptimizer(learning_rate=opt_lr,momentum=opt_mom).minimize(cross_entropy)
    #train_step = tf.train.AdamOptimizer(learning_rate=lrate,epsilon=opt_eps).minimize(cross_entropy,global_step=gstep)
    #train_step = tf.train.GradientDescentOptimizer(0.3).minimize(cross_entropy)
    train_step =  tf.train.AdamOptimizer(learning_rate=lrate,epsilon=opt_eps).minimize(loss,global_step=gstep)
    logging.info("Setting up session...")
    sess = tf.Session()
    init_op = tf.initialize_all_variables()
    sess.run(init_op)

    # Create a saver to save the DNN.
    saver = tf.train.Saver()

    #  previously trained data.
    if(not train_init):
        logging.info("Restoring previously trained net from file {0}".format(fn_saver))
        saver.restore(sess,fn_saver)

    return (sess,train_step,loss,x_input,y_,y_out,saver,gstep)

# -----------------------------------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------------------------------

# Set up the DNN.
(sess,train_step,loss,x,y_,y_out,saver,gstep) = net_setup()



# Open the necessary files.
f_acc = open(fn_acc,'w')
f_dat = tables.open_file(fn_data, 'r')


# Read in a validation set for short checks on accuracy.
dat_val = np.zeros([batch_size,nsensors]); lbl_val = np.zeros([batch_size,ngrid])
(dat_val[:],lbl_val[:]) = read_data(f_dat,ntrain_evts,ntrain_evts+batch_size)



# Set up the arrays for the training and validation evaluation datasets.
dat_train_eval = []; lbl_train_eval = []
dat_test_eval = []; lbl_test_eval = []

# Iterate over all epoch blocks.
for eblk in range(num_epoch_blks):

    logging.info("\n\n**EPOCH BLOCK {0}".format(eblk))

    # Iterate over data blocks.
    for dtblk in range(num_dt_blks):

        logging.info("- DATA BLOCK {0}".format(dtblk))

        # Read in the data.
        if(num_dt_blks > 1 or eblk == 0):
            evt_start = dtblk*dtblk_size
            evt_end = (dtblk+1)*dtblk_size
            dat_train = np.zeros([dtblk_size,nsensors])
            lbl_train = np.zeros([dtblk_size,ngrid])
            #gc.collect()  # force garbage collection to free memory

            #print('done2! ---')
            #sys.stdout.flush()
            (dat_train[0:dtblk_size],lbl_train[0:dtblk_size]) = read_data(f_dat,evt_start,evt_end)
            #print('done2!')
            #sys.stdout.flush()
        # Iterate over epochs within the block.
        for ep in range(epoch_blk_size):

            logging.info("-- EPOCH {0} of block size {1}".format(ep,epoch_blk_size))

            # Shuffle the data.
            logging.info("--- Shuffling data...")
            perm = np.arange(len(dat_train))
            np.random.shuffle(perm)
            dat_train = dat_train[perm]
            lbl_train = lbl_train[perm]

            # Train the NN in batches.
            for bnum in range(batches_per_epoch):

                logging.info("--- Training batch {0} of {1}".format(bnum,batches_per_epoch))

                batch_xs = dat_train[bnum*batch_size:(bnum + 1)*batch_size,:]
                batch_ys = lbl_train[bnum*batch_size:(bnum + 1)*batch_size,:]
                _, loss_val = sess.run([train_step, loss], feed_dict={x: batch_xs, y_: batch_ys})
                logging.info("--- [Step {0}] Got loss value of {1}".format(sess.run(gstep),loss_val))

            # Run a short accuracy check.
            acc_train = 0.; acc_test = 0.
            ltemp,ytemp = sess.run([loss,y_out],feed_dict={x: dat_train[0:batch_size], y_: lbl_train[0:batch_size]})
            for yin,yout in zip(lbl_train[0:batch_size],ytemp):
                #if(np.argmax(yin) == np.argmax(yout)): acc_train += 1
                acc_train += ((yin - yout)*max_xy)**2
            acc_train /= batch_size
            ltemp,ytemp = sess.run([loss,y_out],feed_dict={x: dat_val, y_: lbl_val})
            for yin,yout in zip(lbl_val,ytemp):
                #print(yin)
                #(yout)
                #sys.stdout.flush()

                #if(np.argmax(yin) == np.argmax(yout)): acc_test += 1
                acc_test += ((yin - yout)*max_xy)**2
            acc_test /= batch_size
            logging.info("--- Training MSE = {0}; Test MSE = {1}".format(acc_train,acc_test))

    # Calculate the number of epochs run.
    epoch = eblk*epoch_blk_size
    logging.info("Checking accuracy after {0} epochs".format(epoch+1))

    # Read in the data to be used in the accuracy check.
    if(len(dat_train_eval) == 0):
        dat_train_eval = np.zeros([nval_evts,nsensors]); lbl_train_eval = np.zeros([nval_evts,ngrid])
        (dat_train_eval[:],lbl_train_eval[:]) = read_data(f_dat,0,nval_evts)

    if(len(dat_test_eval) == 0):
        dat_test_eval = np.zeros([nval_evts,nsensors]); lbl_test_eval = np.zeros([nval_evts,ngrid])
        (dat_test_eval[:],lbl_test_eval[:]) = read_data(f_dat,ntrain_evts,ntrain_evts+nval_evts)

    # Run the accuracy check.
    eval_performance(f_acc,epoch,sess,loss,y_out,dat_train_eval,dat_test_eval,lbl_train_eval,lbl_test_eval)

    # Save the trained model.
    logging.info("Saving trained model to: {0}".format(fn_saver))
    save_path = saver.save(sess, fn_saver)

# Close the relevant files.
f_acc.close()
f_dat.close()
