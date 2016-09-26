"""
Panda Functions
JJGC, September 2016
"""
#from __future__ import print_function
import pandas as pd
import numpy as np

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

def get_energy_sensors(energy_v):
    """
    reads the sensors energy and returns a data frame
    """        
    return pd.DataFrame(energy_v.read())

def get_waveforms(pmtea,event_number=0):
    """
    Takes the earray pmtea and returns a DF for event_number
    """
    
    PMTWF ={}
    NPMT = pmtea.shape[1]
    
    for j in range(NPMT):
        PMTWF[j] = pmtea[event_number, j] #waveform for event event_number, PMT j
       
    return pd.DataFrame(PMTWF)

def get_waveforms_and_energy(pmtea,event_number=0):
    """
    Takes the earray pmtea and returns a DF for the wf
    and a Series with the sum of the energies for event_number
    """
    
    PMTWF ={}
    EPMT = []
    NPMT = pmtea.shape[1]
    
    for j in range(NPMT):
        PMTWF[j] = pmtea[event_number, j] #waveform for event event_number, PMT j
        epmt = np.sum(PMTWF[j])
        EPMT.append(epmt)
    return pd.DataFrame(PMTWF), pd.Series(EPMT)

def get_energy(pmtea,event_list=[0]):
    """
    Takes the earray pmtea and a list of events and returns a DF
    with the sum of the energies for event_number
    """
    
    NPMT = pmtea.shape[1]
    epmt = np.zeros(NPMT)
    EPMT=[]
    
    for i in event_list:
        for j in range(NPMT):
            epmt[j] = np.sum(pmtea[i, j])
        EPMT.append(epmt)
        
    return pd.DataFrame(EPMT)

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
