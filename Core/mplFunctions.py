"""
A utility module for plots with matplotlib
"""
from __future__ import print_function
from Util import *
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import numpy as np
import wfmFunctions as wfm
import tblFunctions as tbl

from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

# histograms, signals and shortcuts

def hbins(x, nsigma=5, nbins=10):
  """
  Given an array x, hbins returns the number of bins
  in an interval of  [<x> - nsigma*std(x), <x> + nsigma*std(x)]
  """
  xmin =np.average(x) - nsigma*np.std(x)
  xmax =np.average(x) + nsigma*np.std(x)
  bins = np.linspace(xmin, xmax, nbins)
  return bins

def histo(x,nbins,title='hsimple',xlabel = '', ylabel = 'Frequency'):
  """
  histograms
  """

  plt.hist(x, nbins, histtype='bar', alpha=0.75)
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)


def HSimple1(x,nbins,title='hsimple',xlabel = '', ylabel = 'Frequency',
             save=False,filename='hsimple.png', filepath='./'):
  """
  an interface for plt.hist with some decorations and default options
  """

  plt.hist(x, nbins, histtype='bar', alpha=0.75)
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)

  if save:
    pathfile = filepath+filename
    print("saving histogram %s in %s"%(filename, pathfile))
    plt.savefig(pathfile, bbox_inches='tight')
    plt.close()
  else:
    plt.figure()


def plts(signal, signal_start=0, signal_end=1e+4, offset=5):
  """
  plots a signal in a give interval, control offset by hand
  """
  ax1 = plt.subplot(1,1,1)
  ymin =np.amin(signal[signal_start:signal_end]) - offset
  ymax =np.amax(signal[signal_start:signal_end]) + offset
  ax1.set_xlim([signal_start, signal_end])
  ax1.set_ylim([ymin, ymax])
  plt.plot(signal)

def plot_signal(signal_t,signal,
                title = 'signal', signal_start=0, signal_end=1e+4, units=''):
  """
  Given a series signal (t, signal), plot the signal
  """

  ax1 = plt.subplot(1,1,1)
  ax1.set_xlim([signal_start, signal_end])
  SetPlotLabels(xlabel='t (ns)', ylabel='signal (%s)'%units)
  plt.title(title)
  plt.plot(signal_t, signal)
  plt.show()

def SetPlotLabels(xlabel="", ylabel="",grid=True):
  """
  Short cut to set labels in plots
  """
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if grid == True:
    plt.grid(which='both', axis='both')


#### Circles!
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

# waveforms, sensors, tracks

# waveforms

def plot_waveforms(pmtwfdf, maxlen=0, zoom = False, window_size = 800):
    """
    Takes as input a df storing the PMT wf and plots the 12 PMT WF
    """

    plt.figure(figsize=(12,12))

    len_pmt = len(pmtwfdf[0])

    if maxlen > 0:
        len_pmt = maxlen
    for i in range(12):
        first, last = define_window(pmtwfdf[i],window_size) if zoom else (0, len(pmtwfdf[i]))
        ax1 = plt.subplot(3,4,i+1)
        #ax1.set_xlim([0, len_pmt])
        SetPlotLabels(xlabel='samples', ylabel='adc')
        plt.plot(pmtwfdf[i][first:last])


    plt.show()


def scan_waveforms(pmtea,list_of_events=[0]):
    """
    Takes the earray pmtea and a list of events and scan the waveforms
    """

    for event in list_of_events:
        plot_waveforms(wfm.get_waveforms(pmtea,event_number=event))
        wait()

def define_window( wf, window_size ):
    peak = np.argmax(abs(wf-np.mean(wf)))
    return max(0,peak-window_size),min(len(wf),peak+window_size)

def overlap_waveforms(wfset,event,zoom = True, window_size = 800):
    '''
        Draw all waveforms together.
    '''
    wfs = wfset[event]
    first, last = define_window(wfs[0],window_size) if zoom else (0, wfs.shape[1])
    for wf in wfs:
        plt.plot(wf[first:last])

def compare_raw_blr( pmtrwf, pmtblr, evt = 0, zoom = True, window_size = 800 ):
    '''
        Compare PMT RWF and BLR WF. Option zoom takes a window around the peak
        of size window_size.
    '''
    plt.figure( figsize = (12,12) )
    for i,(raw,blr) in enumerate(zip(pmtrwf[evt],pmtblr[evt])):
        first, last = define_window(raw,window_size) if zoom else (0, pmtrwf.shape[2])
        splot = plt.subplot(3,4,i+1)
        plt.plot(raw[first:last])
        plt.plot(blr[first:last])

def compare_corr_raw( pmtcwf, pmtblr, evt = 0, zoom = True, window_size = 800 ):
    '''
        Compare PMT CWF and RWF (or BLR). Option zoom takes a window around the peak
        of size window_size.
    '''
    transform = lambda wf: 4096 - wf - (4096-2500)
    pmtblr = map( transform, pmtblr )
    plt.figure( figsize = (12,12) )
    for i,(raw,blr) in enumerate(zip(pmtcwf[evt],pmtblr[evt])):
        first, last = define_window(raw,window_size) if zoom else (0, pmtcwf.shape[2])
        splot = plt.subplot(3,4,i+1)
        plt.plot(raw[first:last])
        plt.plot(blr[first:last])

def plot_pmtwf(PMTWF):
    """
    Plots pmtwf
    """

    pmtwf = PMTWF[0]
    plt.plot(pmtwf['time_mus']/mus,pmtwf['ene_pes'])
    ene = pmtwf['ene_pes'].values/12.
    time = pmtwf['time_mus'].values/mus
    plt.xlabel('t (mus)')
    plt.ylabel('E (pes)')
    plt.show()
    plt.figure(figsize=(12,12))
    for i in range(1,len(PMTWF)):
        ax1 = plt.subplot(3,4,i)
        pmtwf = PMTWF[i]
        plt.plot(pmtwf['time_mus']/mus,pmtwf['ene_pes'])
        plt.plot(time,ene)

    plt.show()


def plot_sensor(geom_df,sensor_df, energy_df, event=0, radius=10):
    """
    plots the energy of the sensors
    """
    x =sensor_df['x'].values
    y =sensor_df['y'].values
    r =np.ones(len(sensor_df['x'].values))*radius
#    col = energy_df[event].values ### BUG! we were taking columns

 #   col = energy_df.iloc[[event]].values.flatten() JMB fix
    col = energy_df.ix[event].values # another fix more concise

    plt.figure(figsize=(10,10))
    ax = plt.subplot(aspect='equal')
    circles(x, y, r, c=col, alpha=0.5, ec='none')
    plt.colorbar()
    #xlim(-198,198)  #one should use geom info
    #ylim(-198,198)
    xlim(geom_df['xdet_min'],geom_df['xdet_max'])
    ylim(geom_df['ydet_min'],geom_df['ydet_max'])
    return col

def plot_ene_pmt(geom_df,sensor_df, epmt, event_number=0, radius=10):
    """
    plots the reconstructed energy of the PMTs
    energy_se is a series describing the reconstructed energy
    in each PMT
    """
    x =sensor_df['x'].values
    y =sensor_df['y'].values
    r =np.ones(len(sensor_df['x'].values))*radius
    col = epmt[event_number]

    plt.figure(figsize=(10,10))
    ax = plt.subplot(aspect='equal')
    circles(x, y, r, c=col, alpha=0.5, ec='none')
    plt.colorbar()
    #xlim(-198,198)  #one should use geom info
    #ylim(-198,198)
    xlim(geom_df['xdet_min'],geom_df['xdet_max'])
    ylim(geom_df['ydet_min'],geom_df['ydet_max'])
    return col

def plot_best(sipmrwf,sipmtwf, sipmdf, evt = 0):
    '''
        Plot the noisy waveform of the SiPM with greatest charge and superimpose the true waveform.
    '''
    plt.figure( figsize = (10,8) )
    #Find SiPM with greatest peak
    maxsipm = np.unravel_index(sipmrwf[evt].argmax(),sipmrwf[evt].shape)[0]
    print("SiPM with greatest peak is at index {} with ID {}".format(maxsipm,sipmdf.ix[maxsipm].channel))

    # Plot noisy waveform in red and noiseless waveform in blue
    true_times, true_amps = tbl.read_wf(sipmtwf,evt,maxsipm)
    plt.plot(sipmrwf[evt,maxsipm,:])
    plt.plot(true_times,np.array(true_amps)*sipmdf['adc_to_pes'][maxsipm])
    plt.xlabel('time ($\mu$s)');plt.ylabel('Energy (adc)');

def plot_best_group(sipmrwf,sipmtwf,sipmdf,evt = 0, nsipms = 8, ncols = 3):
    '''
        Plot the noisy (red) and true (blue) waveforms of the nsipms SiPMs with greatest charge.
    '''
    plt.figure( figsize = (10,8) )
    #Find SiPM with greatest peak
    sipms = sorted( enumerate(sipmrwf[evt]), key = lambda x: max(x[1]), reverse = True )[:nsipms]
    nrows = int(ceil(nsipms*1.0/ncols))

    for i,(sipm_index, sipm_wf) in enumerate(sipms):
        try:
            true_times, true_amps = tbl.read_wf(sipmtwf,evt,sipm_index)
        except:
            true_times, true_amps = np.arange(len(sipm_wf)),np.zeros(len(sipm_wf))
        plt.subplot(nrows,ncols,i+1)
        plt.plot(sipm_wf)
        plt.plot(true_times,np.array(true_amps)*sipmdf['adc_to_pes'][sipm_index])
        plt.xlabel('time ($\mu$s)');plt.ylabel('Energy (adc)');
    plt.tight_layout()

def plot_track(geom_df,mchits_df,vox_size=10, zoom = False):
    """
    plot the hits of a mctrk. Adapted from JR plotting functions
    notice that geom_df and mchits_df are pandas objects defined above
    if zoom = True, the track is zoomed in (min/max dimensions are taken from track)
    if zoom = False, detector dimensions are used
    """
    from mpl_toolkits.mplot3d import Axes3D


    grdcol = 0.99

    varr_x = mchits_df['x'].values*vox_size
    varr_y = mchits_df['y'].values*vox_size
    varr_z = mchits_df['z'].values*vox_size
    varr_c = mchits['energy'].values/keV

    min_x = geom_df['xdet_min']*vox_size
    max_x = geom_df['xdet_max']*vox_size
    min_y = geom_df['ydet_min']*vox_size
    max_y =geom_df['ydet_max']*vox_size
    min_z = geom_df['zdet_min']*vox_size
    max_z = geom_df['zdet_max']*vox_size
    emin=0
    emax=np.max(varr_c)

    if zoom == True:
        min_x = np.min(varr_x)
        max_x = np.max(varr_x)
        min_y = np.min(varr_y)
        max_y = np.max(varr_y)
        min_z = np.min(varr_z)
        max_z = np.max(varr_z)
        emin=np.min(varr_c)

    # Plot the 3D voxelized track.
    fig = plt.figure(1)
    fig.set_figheight(6.)
    fig.set_figwidth(8.)


    ax1 = fig.add_subplot(111,projection='3d');
    s1 = ax1.scatter(varr_x,varr_y,varr_z,marker='s',linewidth=0.5,
                         s=2*vox_size,c=varr_c,cmap=plt.get_cmap('rainbow'),
                         vmin=emin,vmax=emax)

    # this disables automatic setting of alpha relative of distance to camera
    s1.set_edgecolors = s1.set_facecolors = lambda *args:None;


    print(" min_x ={} max_x ={}".format(min_x,max_x))
    print(" min_y ={} max_y ={}".format(min_y,max_y))
    print(" min_z ={} max_z ={}".format(min_z,max_z))
    print("min_e ={} max_e ={}".format(emin, emax))

    ax1.set_xlim([min_x, max_x])
    ax1.set_ylim([min_y, max_y])
    ax1.set_zlim([min_z, max_z])

    #    ax1.set_xlim([0, 2 * vox_ext]);
    #    ax1.set_ylim([0, 2 * vox_ext]);
    #    ax1.set_zlim([0, 2 * vox_ext]);
    ax1.set_xlabel("x (mm)");
    ax1.set_ylabel("y (mm)");
    ax1.set_zlabel("z (mm)");
    ax1.set_title("");

    lb_x = ax1.get_xticklabels();
    lb_y = ax1.get_yticklabels();
    lb_z = ax1.get_zticklabels();
    for lb in (lb_x + lb_y + lb_z):
        lb.set_fontsize(8);

    ax1.w_xaxis.set_pane_color((1.0,1.0,1.0,1.0));
    ax1.w_yaxis.set_pane_color((1.0,1.0,1.0,1.0));
    ax1.w_zaxis.set_pane_color((1.0,1.0,1.0,1.0));
    ax1.w_xaxis._axinfo.update({'grid' : {'color': (grdcol, grdcol, grdcol, 1)}});
    ax1.w_yaxis._axinfo.update({'grid' : {'color': (grdcol, grdcol, grdcol, 1)}});
    ax1.w_zaxis._axinfo.update({'grid' : {'color': (grdcol, grdcol, grdcol, 1)}});

    cb1 = plt.colorbar(s1);
    cb1.set_label('Energy (keV)');

    plt.show()

def plot_track_projections(geom_df,mchits_df,vox_size=10, zoom = False):
    """
    plot the projections of an MC track. Adapted from function plot_track above
    notice that geom_df and mchits_df are pandas objects defined above
    if zoom = True, the track is zoomed in (min/max dimensions are taken from track)
    if zoom = False, detector dimensions are used

    For now, it is assumed that vox_sizeX = vox_sizeY = vox_sizeZ

    """
    from mpl_toolkits.mplot3d import Axes3D

    grdcol = 0.99

    vox_sizeX = vox_size
    vox_sizeY = vox_size
    vox_sizeZ = vox_size

    varr_x = mchits_df['x'].values*vox_size
    varr_y = mchits_df['y'].values*vox_size
    varr_z = mchits_df['z'].values*vox_size
    varr_c = mchits['energy'].values/keV

    min_x = geom_df['xdet_min']*vox_size
    max_x = geom_df['xdet_max']*vox_size
    min_y = geom_df['ydet_min']*vox_size
    max_y =geom_df['ydet_max']*vox_size
    min_z = geom_df['zdet_min']*vox_size
    max_z = geom_df['zdet_max']*vox_size
    emin=0
    emax=np.max(varr_c)

    if zoom == True:
        min_x = np.min(varr_x)
        max_x = np.max(varr_x)
        min_y = np.min(varr_y)
        max_y = np.max(varr_y)
        min_z = np.min(varr_z)
        max_z = np.max(varr_z)
        emin=np.min(varr_c)

    # Plot the 2D projections.
    fig = plt.figure(1);
    fig.set_figheight(5.);
    fig.set_figwidth(20.);

    # Create the x-y projection.
    ax1 = fig.add_subplot(131);
    hxy, xxy, yxy = np.histogram2d(varr_y, varr_x, weights=varr_c, normed=False,
                                   bins=((1.05*max_y - 0.95*min_y)/vox_sizeY, (1.05*max_x - 0.95*min_x)/vox_sizeX),
                                   range=[[0.95*min_y,1.05*max_y],[0.95*min_x,1.05*max_x]])

    extent1 = [yxy[0], yxy[-1], xxy[0], xxy[-1]]
    sp1 = ax1.imshow(hxy, extent=extent1, interpolation='none', aspect='auto', origin='lower')
    ax1.set_xlabel("x (mm)")
    ax1.set_ylabel("y (mm)")
    cbp1 = plt.colorbar(sp1);
    cbp1.set_label('Energy (keV)');

    # Create the y-z projection.
    ax2 = fig.add_subplot(132);
    hyz, xyz, yyz = np.histogram2d(varr_z, varr_y, weights=varr_c, normed=False,
                                   bins=((1.05*max_z - 0.95*min_z)/vox_sizeZ, (1.05*max_y - 0.95*min_y)/vox_sizeY),
                                   range=[[0.95*min_z,1.05*max_z],[0.95*min_y,1.05*max_y]])
    extent2 = [yyz[0], yyz[-1], xyz[0], xyz[-1]]
    sp2 = ax2.imshow(hyz, extent=extent2, interpolation='none', aspect='auto', origin='lower')
    ax2.set_xlabel("y (mm)")
    ax2.set_ylabel("z (mm)")
    cbp2 = plt.colorbar(sp2);
    cbp2.set_label('Energy (keV)');

    # Create the x-z projection.
    ax3 = fig.add_subplot(133);
    hxz, xxz, yxz = np.histogram2d(varr_z, varr_x, weights=varr_c, normed=False,
                                   bins=((1.05*max_z - 0.95*min_z)/vox_sizeZ, (1.05*max_x - 0.95*min_x)/vox_sizeX),
                                   range=[[0.95*min_z,1.05*max_z],[0.95*min_x,1.05*max_x]])
    extent3 = [yxz[0], yxz[-1], xxz[0], xxz[-1]]
    sp3 = ax3.imshow(hxz, extent=extent3, interpolation='none', aspect='auto', origin='lower')
    ax3.set_xlabel("x (mm)")
    ax3.set_ylabel("z (mm)")
    cbp3 = plt.colorbar(sp3);
    cbp3.set_label('Energy (keV)');


    plt.show()
