"""
Utilities
"""
#from __future__ import print_function
import pandas as pd
from PlotUtil import *

def read_data_geom(geom_t):
    """
    Reads the geom table en returns a PD Series
    """
        
    ga = geom_t.read()
    G = pd.Series([ga[0][0][0],ga[0][0][1],ga[0][1][0],ga[0][1][1],
                    ga[0][2][0],ga[0][2][1],ga[0][3]],
                    index=['xdet_min','xdet_max','ydet_min','ydet_max',
                            'zdet_min','zdet_max','R'])
    return G

def read_data_FEE(fee_t):
    """
    Reads the FEE table en returns a PD Series for the simulation parameters 
    and a PD series for the values of the capacitors used in the simulation
    """
        
    fa = fee_t.read()
    F = pd.Series([fa[0][0],fa[0][1],fa[0][2],fa[0][3],fa[0][5],fa[0][6],fa[0][7],fa[0][8],fa[0][9],fa[0][10],
                   fa[0][11],fa[0][12]],
                    index=['offset','pmt_gain','V_gain','R',"time_step",
                           "time_daq",
                            "freq_LPF",
                            "freq_HPF",
                            "LSB",
                            "volts_to_adc",
                            "noise_fee_rms",
                            "noise_adc"])
    C =pd.Series([fa[0][4]],index=['C12'])
    return F,C

def get_column_(pmta,ic):
    """
    access column ic of table pmta and returns column as an array
    """
    col =[]
    for i in range(pmta.shape[0]):
        col.append(pmta[i][ic])
    return np.array(col)
 
def read_data_sensors(sensor_table):
    """
    reads the sensors table and returns a data frame
    """
    pmta = sensor_table.read()
    PMT={}
    PMT['channel'] = get_column_(pmta,0)
    PMT['active'] = get_column_(pmta,1)
    PMT['x'] = get_column_(pmta,2).T[0]
    PMT['y'] = get_column_(pmta,2).T[1]
    PMT['gain'] = get_column_(pmta,3)
    PMT['adc_to_pes'] = get_column_(pmta,4)
        
    return pd.DataFrame(PMT)

def read_energy_sensors(energy_v):
    """
    reads the sensors energy and returns a data frame
    """        
    return pd.DataFrame(energy_v.read())

def get_mctrks(mctrk,event_number = 0):
    """
    return all the mc trks in an event
    """
    mcparticle ={}
    mc_name = []
    mc_pdg = []
    mc_vxi = []
    mc_vxf= []
    mc_nhits= []
    mc_energy= []
    for row in mctrk.iterrows():
        if row['hit_indx'] == 0 and row['event_indx'] == event_number:
            mc_name.append(row['particle_name'])
            mc_pdg.append(row['pdg_code'])
            mc_vxi.append(row['initial_vertex'])
            mc_vxf.append(row['final_vertex'])
            mc_nhits.append(row['nof_hits'])
            mc_energy.append(row['energy']) 
                
    mcparticle['name'] = mc_name
    mcparticle['pdg'] = mc_pdg
    mcparticle['vxi'] = mc_vxi
    mcparticle['vxf'] = mc_vxf
    mcparticle['nhits'] = mc_nhits
    mcparticle['energy'] = mc_energy
                
    return pd.DataFrame(mcparticle)

def get_mchits(mctrk,event_number = 0, particle_number =0):
    """
    return the mc hits of a mc particle in an event
    Takes the pointer to the table (not to the mctrk data frame)
    """
    mchits ={}
    hit_x = []
    hit_y = []
    hit_z = []
    hit_time = []
    hit_energy = []
    
    for row in mctrk.iterrows():
        if row['mctrk_indx'] == particle_number and row['event_indx'] == event_number:
            hit_x.append(row['hit_position'][0])
            hit_y.append(row['hit_position'][1])
            hit_z.append(row['hit_position'][2])
            hit_time.append(row['hit_time'])
            hit_energy.append(row['hit_energy'])
            
                
    mchits['x'] =hit_x
    mchits['y'] =hit_y
    mchits['z'] =hit_z
    mchits['time'] = hit_time
    mchits['energy'] = hit_energy
                
    return pd.DataFrame(mchits)

def plot_sensor(geom_df,sensor_df, energy_df, event=0, radius=10):
    """
    plots the energy of the sensors (in pes)
    """
    x =sensor_df['x'].values
    y =sensor_df['y'].values
    r =np.ones(len(sensor_df['x'].values))*radius
#    col = energy_df[event].values
    col = energy_df.iloc[[event]].values.flatten()
    
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

