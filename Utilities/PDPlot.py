"""
Panda Plotting functions
JJG September 2016
"""
#from __future__ import print_function
import pandas as pd
from PlotUtil import *



def plot_sensor(geom_df,sensor_df, energy_df, event=0, radius=10):
    """
    plots the energy of the sensors (in pes)
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

def plot_waveforms(pmtwfdf):
    """
    Takes as input a df storing the PMT wf and plots the 12 PMT WF
    """
    
    plt.figure(figsize=(12,12))
    
    len_pmt = len(pmtwfdf[0])
    for i in range(12):
        ax1 = plt.subplot(3,4,i+1)
        ax1.set_xlim([0, len_pmt])
        SetPlotLabels(xlabel='t (ns)', ylabel='adc')
        plt.plot(pmtwfdf[i])

    
    plt.show()

def scan_waveforms(pmtea,list_of_events=[0]):
    """
    Takes the earray pmtea and a list of events and scan the waveforms
    """
    
    for event in list_of_events:
        plot_waveforms(get_waveforms(pmtea,event_number=event))
        wait()

