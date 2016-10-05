# toy_plots.py holds the plotting functions for analysis of the
# toy reconstruction!
from collections         import defaultdict
from matplotlib.patches  import Ellipse

from toy_inputs          import *

import matplotlib.pyplot as plt
import numpy as np


# training/cv loss error over epochs of training
def plot_net_history(history):
    val   = plt.plot(range(0,len(history['loss'])),
                     np.sqrt(np.array(history['val_loss'])/N_ELpts)*max_xy, 'b', label='Val error')
    train = plt.plot(range(0,len(history['loss'])),
                     np.sqrt(np.array(history['loss'])    /N_ELpts)*max_xy, 'r', label='Train error')
    plt.legend(loc='upper right')
    plt.xlabel('Epochs')
    plt.ylabel('Mean error per hit (mm)')
    plt.title(str(N_ELpts) + ' hit(s), ' + str(N_layers) +' layers')
    plt.grid(True)
    plt.show()

# Error in x and in y as a function of the x and y coordinate of the hit
def plot_xy_error(predictions,y_):
    xerr = abs(predictions[:,0] - y_[:,0])*max_xy
    yerr = abs(predictions[:,1] - y_[:,1])*max_xy

    # could replace these best fit lines with something average data points
    # x, y error as a function of x, y  EL hit coordinate
    #plt.scatter(y_[0:10,0]*max_xy,xerr[0:10],color='b') # sample points
    #plt.scatter(y_[0:10,1]*max_xy,yerr[0:10],color='g') # sample points
    bestfitx = plt.plot(np.unique(y_[:,0]*max_xy), np.poly1d(np.polyfit(y_[:,0]*max_xy,xerr, 8))(np.unique(y_[:,0]*max_xy)),label='x')
    bestfitx = plt.plot(np.unique(y_[:,1]*max_xy), np.poly1d(np.polyfit(y_[:,1]*max_xy,yerr, 8))(np.unique(y_[:,1]*max_xy)),label='y')
    plt.legend(loc='upper right')
    plt.xlabel('EL hit coordinate')
    plt.ylabel('EL hit error (mm)')
    plt.show()

# plot the histogram, of average error per point in distance
def plot_histogram(n_bins,xlim):
    n, bins, patches = plt.hist(E.flatten(), n_bins, color='red', normed=1, alpha=0.75)
    plt.xlim(0,xlim)
    plt.xlabel('Error/Hit (mm)')
    plt.grid(True)
    plt.title(str(N_ELpts) + ' hit(s), ' + str(N_layers) +' layers')
    plt.show()

# Plot error as a function of the distance of the EL hit from the center of the EL plane (40,40)
def plot_radial_error():
    dd=defaultdict(int) # defaultdict whose keys are distance from 40,40, whose values are error
    dn=defaultdict(int) # defaultdict whose keys are distance from 40,40, whose values are the number of evts
                        # with that distance

    for ev in range(Nevts_valid):
        for pt in range(N_ELpts):
            dist = np.sqrt((y_valid[ev,2*pt]*xlen-40)**2 + (y_valid[ev,2*pt+1]*ylen-40)**2)
            dd[dist] += E[ev,pt]
            dn[dist] += 1

    for key in dd: dd[key] /= dn[key]

    # this line is pretty meaningless since I am picking the degree polynomial --------->
    bestfit = plt.plot(np.unique(dd.keys()), np.poly1d(np.polyfit(dd.keys(),dd.values(), 5))(np.unique(dd.keys())))

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

# Plot a sipm map with a label and prediction
def plot_event():
    plt_train = False
    pevt = np.random.randint(0,Nevts_valid)
    fig = plt.figure();
    ax1 = fig.add_subplot(111);
    ax1.axis([0, xlen, 0, ylen]);

    if(plt_train):
        xarr = x_train[pevt]
        yarr = y_train[pevt]
    else:
        xarr = x_valid[pevt]
        yarr = y_valid[pevt]

    # Create circles and plot them according to the probabilities.
    probs = (xarr - min(xarr))
    probs /= max(probs)
    for x,y,p in zip(pos_x, pos_y, probs):
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
