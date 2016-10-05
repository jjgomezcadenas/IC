

# Here execute a Corona algorithm so that we may compare its accuracy to the that of a neural network
import copy
import numpy as np
import tables as tb
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse




def corona(m,coords):
    """
    corona takes as input a map
    outputs a list of tuples, each tuple contains the coordinates of a EL hit
    called by toy_corona
    """
    x=0
    y=0
    hit_energy = 0
    mt = copy.deepcopy(m)
    (xi,yi) = np.unravel_index(mt.argmax(), mt.shape)

    if mt[xi,yi] > 4.1:

        for h in range(-3,4):
            for v in range(-3,4):
                try:
                    x += (xi+h) * mt[xi+h,yi+v]
                    y += (yi+v) * mt[xi+h,yi+v]
                    hit_energy += mt[xi+h,yi+v]
                    mt[xi+h,yi+v] = 0

                except: pass

        x /= hit_energy
        y /= hit_energy

        coords.append((x,y))
        corona(mt,coords)

    else: return
    return

def toy_corona_error(filepath,plot=False):
    """
    input: a data file path; if plot==True, plots error as a function of coord
    returns: error of corona algorithm
    uses: corona
    """
    # Load SiPM maps
    f    = tb.open_file(filepath,'r')
    data = f.root.sim_1pt
    maps        = data.xvalid[:]
    real_coords = data.yvalid[:]
    for m in maps: m -= np.min(m) # make minimum=0 (min is negative)

    E  = np.empty((maps.shape[0]))
    xe = np.empty((maps.shape[0]))
    xc = np.empty((maps.shape[0]))
    ye = np.empty((maps.shape[0]))
    yc = np.empty((maps.shape[0]))

    e = 0
    for m,c in zip(maps,real_coords):
        coords =[]
        m    = m.reshape((8,8)).T
        corona(m,coords)
        E[e] = np.sqrt(np.sum(((np.array(coords[0]))*10+5 - c*80)**2))
        xe[e] = abs(coords[0][0]*10+5 - c[0]*80)
        xc[e] = c[0]*80
        ye[e] = abs(coords[0][1]*10+5 - c[1]*80)
        yc[e] = c[1]*80
        e += 1

    print('Mean error with 1 hit (mm): ' + str(np.mean(E)))

    if plot:
        # x, y error as a function of x, y  EL hit coordinate
        plt.scatter(xc[0:50],xe[0:50],color='b') # sample points
        plt.scatter(yc[0:50],ye[0:50],color='g') # sample points
        bestfitx = plt.plot(np.unique(xc), np.poly1d(np.polyfit(xc,xe, 8))(np.unique(xc)),label='x')
        bestfity = plt.plot(np.unique(yc), np.poly1d(np.polyfit(yc,ye, 8))(np.unique(yc)),label='y')
        plt.legend(loc='upper right')
        plt.xlabel('EL hit coordinate')
        plt.ylabel('EL hit error (mm)')
        plt.show()

    return(xc,xe,yc,ye)
