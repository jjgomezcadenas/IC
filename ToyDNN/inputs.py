import logging

# ---------------------------------------------------------------------------------------------------------
# USER INPUTS
# ---------------------------------------------------------------------------------------------------------

# Directory and run names
datdir = "/Users/alej/Desktop/Valencia/DNNtoy/indata"                   # the base data directory
dname = "resp"         # the data name
rdir = "/Users/alej/Desktop/Valencia/DNNtoy/outdata"             # the run directory
rname = "sim_1_elpt"

# Net configuration parameters
net_name = "MNISTadv_nodrop"                    # name of the neural net described in nets/neuralnets.py
train_init = True                          # if true, train from net with standard pre-training; if false, read in a previously trained net

# Response map parameters
nsipm = 8
ngrid = 2
max_xy = 80

# Parameters describing training intervals and number of events for training and validation
ntrain_evts = 40000    # number of training evts per dataset
nval_evts = 5000       # number of validation events
num_epochs = 80        # total number of epochs to train
epoch_blk_size = 1     # number of epochs to run per block (before reading new dataset); set equal to num_epochs unless data to be read in multiple blocks
dtblk_size = ntrain_evts      # number of signal and background events per training block
batch_size = 250       # training batch size

# Training optimizer parameters
opt_lr = 1.0e-3       # optimizer learning rate
opt_eps = 1.0e-7      # optimizer epsilon (for AdamOptimizer)
opt_mom = 0.9          # optimizer momentum
opt_decaybase = 0.1    # multiplicative factor for learning rate decay
opt_ndecayepochs = 50  # decay interval: apply decay by a factor of opt_decaybase every opt_ndecayepochs epochs

# Plotting and logging parameters
log_to_file = True         # set to True to output log information to a file rather than the console
logging_lvl = logging.INFO  # logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL
#plt_show = False       # show plots on-screen for dnnplot
#plt_imgtype = "png"    # image type to which to save plots

# END USER INPUTS
# ------------------------------------------------------------------------------------------

# Calculated parameters (based on user inputs)
batches_per_epoch = int(dtblk_size/batch_size)
num_epoch_blks = num_epochs / epoch_blk_size
num_dt_blks = ntrain_evts / dtblk_size
nsensors = nsipm * nsipm
