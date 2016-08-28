"""
A utility module for plots with matplotlib
"""
from Util import *
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import numpy as np
import os
import sys

from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

def hbins(x, nsigma=5, nbins=10):
  xmin =np.average(x) - nsigma*np.std(x)
  xmax =np.average(x) + nsigma*np.std(x)
  bins = np.linspace(xmin, xmax, nbins)
  return bins
  
def HSimple1(x,nbins,title='hsimple',xlabel = '', ylabel = 'Frequency', 
             save=False,filename='hsimple.png', filepath='./'):
  plt.hist(x, nbins, histtype='bar', alpha=0.75)
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  
  if save:
    pathfile = filepath+filename
    print "saving histogram %s in %s"%(filename, pathfile)
    plt.savefig(pathfile, bbox_inches='tight')
    plt.close()
  else:
    plt.figure()


def plot_signal(signal_t,signal, 
                title = 'signal', signal_start=0, signal_end=1e+4, units=''):

  ax1 = plt.subplot(1,1,1)
  ax1.set_xlim([signal_start, signal_end])
  SetPlotLabels(xlabel='t (ns)', ylabel='signal (%s)'%units)
  plt.title(title)
  plt.plot(signal_t, signal)
  plt.show()

def plts(signal, signal_start=0, signal_end=1e+4, offset=5):
  """
  A more concise version
  """

    ax1 = plt.subplot(1,1,1)
    ymin =np.amin(signal[signal_start:signal_end]) - offset
    ymax =np.amax(signal[signal_start:signal_end]) + offset
    ax1.set_xlim([signal_start, signal_end])
    ax1.set_ylim([ymin, ymax])
    plt.plot(signal)




def SetPlotLabels(xlabel="", ylabel="",grid=True):
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if grid == True:
    plt.grid(which='both', axis='both')

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circle in data unit.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    plt.colorbar()

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None
    if 'fc' in kwargs: kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs: kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs: kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs: kwargs.setdefault('linewidth', kwargs.pop('lw'))

    patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(x, y, s)]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        collection.set_array(np.asarray(c))
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    if c is not None:
        plt.sci(collection)
    return collection






