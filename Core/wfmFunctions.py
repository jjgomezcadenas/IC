"""
Waveform Functions
JJGC, September-October 2016

ChangeLog
12/10: change from time_ns to time_mus
"""
#from __future__ import print_function
import pandas as pd
import numpy as np
import FEParam as FP
from system_of_units import *
from Util import dict_map
from math import *
from LogConfig import *


def to_adc( wfs, sensdf ):
    '''
        Scale waveform in pes to adc.
    '''
    return wfs * sensdf['adc_to_pes'].reshape(wfs.shape[0],1)

def to_pes( wfs, sensdf ):
    '''
        Scale waveform in adc to pes.
    '''
    return wfs / sensdf['adc_to_pes'].reshape(wfs.shape[0],1)


def rebin_twf(t, e, stride = 40):
    """
    rebins the a waveform according to stride
    The input waveform is a vector such that the index expresses time bin and the
    contents expresses energy (e.g, in pes)
    The function returns the ned times and energies
    """
    n = int(ceil(len(t)/float(stride)))

    T = np.zeros(n,dtype=np.float32)
    E = np.zeros(n,dtype=np.float32)

    j=0
    for i in range(n):
        E[i] = np.sum(e[j:j+stride])
        T[i] = np.mean(t[j:j+stride])
        j+= stride

    return T,E

def rebin_df(df,stride=40):
    '''
        applies the rebin_wf function to a dataframe.
    '''
    return wf2df(*rebin_twf(*df2wf(df), stride = stride))

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

def wfdf(time_mus,energy_pes,indx):
    """
    takes three vectors (time, energy and indx) and returns a data frame
    representing a waveform
    """
    return pd.DataFrame({'time_mus':time_mus,'ene_pes':energy_pes,'indx':indx})

def wf2df(time_mus,energy_pes, dropnan=False):
    """
    takes two vectors (time, energy) and returns a data frame representing a waveform
    """
    if dropnan == False:
        return pd.DataFrame({'time_mus':time_mus,'ene_pes':energy_pes})
    else:
        return pd.DataFrame({'time_mus':time_mus,'ene_pes':energy_pes}).dropna()

def df2wf(df):
    '''
        takes a data frame and returns the array of times and the array of energies.
    '''
    return df['time_mus'], df['ene_pes']

def add_cwf(cwfdf,pmtDF):
    """
    input: cwfdf: each colum is the wf for one PMT.
    output: swf is a data frame with two columns:
    time_mus = counts the time in ns
    ene_pes = conths the energy in pes
    """
    wf =0
    NPMT = len(pmtDF)
    for i in range(NPMT):
        adc_to_pes = pmtDF['adc_to_pes'][i]
        wf += cwfdf[i].values/adc_to_pes

    return wfdf(np.array(range(len(wf)))*FP.time_DAQ, wf,np.array(range(len(wf))))

def wf_thr(wf,threshold=1):
    """
    return a zero supressed waveform (more generally, the vaules of wf above threshold)
    """
    return wf.loc[lambda df: df.ene_pes.values >threshold, :]

def zs_wf(waveform,threshold,to_mus=None):
    '''
        get a zero-supressed wf.
    '''
    t = np.argwhere(waveform>threshold).flatten()
    if not t.size: return None
    return wf2df( t if to_mus is None else t*to_mus,waveform[t] )

def sensor_wise_zero_suppression(data,thresholds, to_mus=None):
    '''
        takes an array of waveforms, applies the corresponding threshold to
        each row and returns a dictionary with the data frames of the survivors.
    '''
    # If threshold is a single value, transform it into an array
    if not hasattr(thresholds, '__iter__'): thresholds = np.ones( data.shape[0] ) * thresholds
    return { i : df for i,df in enumerate(map(zs_wf,data,thresholds)) if df is not None }

def suppress_wf(wf,th):
    '''
        Put zeros where wf is below th.
    '''
    wf[wf<=th] = 0

def noise_suppression(data,thresholds):
    '''
        takes an array of waveforms, applies the corresponding threshold to
        each row and returns a dictionary with the data frames of the survivors.
    '''
    suppressed_data = np.copy(data)
    if not hasattr(thresholds, '__iter__'):
        thresholds = np.ones( data.shape[0] ) * thresholds
    map( suppress_wf, suppressed_data, thresholds)
    return suppressed_data

def in_window( data, tmin, tmax ):
    '''
        Filters out data outside specified window.
    '''
    filter_df = lambda df: df[ (df.time_mus >= tmin) & (df.time_mus <= tmax) ]
    return dict_map( filter_df, data )

def find_S12(swf, stride=40):
    """
    Find S1 or S2 signals. The input is a zero-supressed WF. The stride defines the contiguity criterium.
    The stride is applied to the indexes which keep the ordering of the original (non-zs) WF.
    For example, with a stride of 40 (corresponding to steps of 1 mus for a DAQ timing of 25 ns) index 1
    and index 39 are in the same S12.
    """
    T = swf['time_mus'].values
    P = swf['ene_pes'].values
    I = swf['indx'].values

    S12 = {}
    pulse_on = 1
    j=0

    S12[0] = []
    S12[0].append([T[0],P[0],I[0]])

    for i in range(1,len(swf)) :
        if swf.index[i]-stride > swf.index[i-1]:  #new s12
            j+=1
            S12[j] = []
            S12[j].append([T[i],P[i],I[i]])
        else:
            S12[j].append([T[i],P[i],I[i]])

    S12L=[]
    for i in S12.keys():
        S12L.append(pd.DataFrame(S12[i], columns=['time_mus','ene_pes','indx']))
    return S12L

def rebin_waveform(swf, stride = 40):
    """
    rebins the a waveform according to stride
    The input waveform is a vector such that the index expresses time bin and the
    contents expresses energy (e.g, in pes)
    The function returns a DataFrame. The time bins and energy are rebinned according to stride
    """

    t = swf['time_mus'].values
    e = swf['ene_pes'].values
    I = swf['indx'].values
    n = len(swf)/int(stride)
    r = len(swf)%int(stride)

    lenb = n
    if r > 0:
        lenb = n+1

    T = np.zeros(lenb)
    E = np.zeros(lenb)
    II = np.zeros(lenb, dtype=int)

    j=0
    for i in range(n):
        E[i] = np.sum(e[j:j+stride])
        T[i] = np.mean(t[j:j+stride])
        II[i] = I[(j+stride)/2]
        j+= stride

    if r > 0:
        E[n] = np.sum(e[j:])
        T[n] = np.mean(t[j:])
        II[n] = I[(len(swf) - j/2)]


    rbw={}
    rbw['ene_pes'] = E
    rbw['time_mus'] = T
    rbw['indx'] = II
    return pd.DataFrame(rbw)

def find_t0(s1):
    """
    returns t0 as the peak of the S1 signal
    """
    emax = np.amax(s1.ene_pes.values)
    return s1.loc[lambda df: df.ene_pes.values ==emax, :]

def s2_energy(s2):
    """
    total energy in pes
    """
    return np.sum(s2.ene_pes.values)

def s2_length(s2):
    """
    s2 length in ns
    """
    t = s2.time_mus.values
    return t[-1] - t[0]

def pmt_wf(cwfdf,pmtDF):
    """
    input: cwfdf: each colum is the wf for one PMT.
    output: wf is a list of waveform data frames
    time_mus = counts the time in ns
    ene_pes = conths the energy in pes
    returns a data frame with one data frame per PMT
    each pmt DF expresses the waveform in the PMT
    """

    NPMT = len(pmtDF)
    PMTWF = []
    for i in range(NPMT):
        adc_to_pes = pmtDF['adc_to_pes'][i]
        wf = cwfdf[i].values/adc_to_pes
        PMTWF.append(wfdf(np.array(range(len(wf)))*FP.time_DAQ, wf,np.array(range(len(wf)))))


    return PMTWF

def pmaps_EP(pmtcwf,pmtDF,list_of_events=[0], thr = 1, stride=40):
    """
    computes pmaps in the EP plane
    Returns a list of S2PMAPS (one entry per event) and a t0 DataFrame
    S2PMAPS is a list of PMAPS of one entry per S2 found in the sum of PMTs
    Each EP PMAP is the collection of 12+1 S2s
    """
    NPMT = len(pmtDF)
    S2PMAP = []
    for event in list_of_events:
        # compute the sum function (swf), supress zeros (swf_zs) and finds s12 for swf_zs
        cwfdf = get_waveforms(pmtcwf,event_number=event)
        swf = add_cwf(cwfdf,pmtDF)
        swf_zs = wf_thr(swf,threshold =thr)
        s12 = find_S12(swf_zs, stride=stride)

        is2=0
        t0 = -999
        if len(s12) > 1:  #if s1 exists is s12[0]
            is2=1  #s2 starts in index 1 if s1 exists in 0 otherwise
            s1 = s12[0]
            t0 = find_t0(s1)

        S2L = []
        for s2 in s12[is2:]: #loop over s2 found is swf_zs
            PMTWF = pmt_wf(cwfdf,PMT) #wf for each of the PMTs
            #scan_pmtwf_s2(PMTWF,s2)

            PMAP = []
            s2rb = rebin_waveform(s2, stride = stride)
            PMAP.append(s2rb)
            for i in range(NPMT):
                pmtwf = PMTWF[i]
                pmtdf = wfdf(pmtwf.time_mus.values[s2.indx.values[0]:s2.indx.values[-1]+1],
                                 pmtwf.ene_pes.values[s2.indx.values[0]:s2.indx.values[-1]+1],
                                 pmtwf.indx.values[s2.indx.values[0]:s2.indx.values[-1]+1])
                pmtrb = rebin_waveform(pmtdf, stride = stride)

                PMAP.append(pmtrb)
            S2L.append(PMAP)
        S2PMAP.append(S2L)
    return t0, S2PMAP

    ### PMAPS

class PMAP:
    """
    A simple class to hold the EP and TP pmap info
    """
    def __init__(self, t0):
        """
        inits the class with t0

        """

        self.t0 = t0
        self.s2PMAP  = []
        self.epPMAP  = []
        self.sipmS2P  = []

    def add_pmap(self,s2pmap, sipms2p, epmap):
        self.s2PMAP.append(s2pmap)
        self.sipmS2P.append(sipms2p)
        self.epPMAP.append(epmap)


    def nof_s2(self):
        return len(self.s2PMAP)


def sipm_panel(sipmrwf, SIPMDF, event_number=0):
    """
    Organize the SiPM as a PD panel, that is a collection of PD DataFrames

    1. items = number of sipm
    2. One DataFrame per SiPM
    """
    sipmwf = sipmrwf[event_number]
    SIPM = {}
    NSIPM = sipmwf.shape[0]
    sipmwl = sipmwf.shape[1]
    for i in range(NSIPM):
        adc_to_pes = SIPMDF.adc_to_pes[i]
        energy_pes = sipmwf[i]/adc_to_pes
        time_mus = np.array(range(sipmwl))*1e+3 #steps are mus
        indx = np.ones(sipmwl)*i
        SIPM[i] = wf_mus(wfdf(time_mus,energy_pes,indx))
    return pd.Panel(SIPM)

def sipm_s2(sipmdf, s2df):
    """
    Takes a sipm DF and an s2df
    Returns a DF with the sipm values in the range specified by s2
    """
    s2ti = s2df.time_mus.values[0]
    s2tf = s2df.time_mus.values[-1]
    dfl = sipmdf.loc[lambda df: df.time_mus.values >= s2ti, :]
    dfu = dfl.loc[lambda df: df.time_mus.values < s2tf, :]
    return dfu

def sipmp_s2(sipmp, s2df, thr=0.5):
    """
    Takes a sipm panel and a s2df
    Returns a sipm panel with a collection of sipm DF such that:
    1. the range of the sipm is specified by s2
    2. the sipm energy are above threshold.
    """
    SIPM={}
    j=0
    for i in sipmp.items:
        sipm = sipmp[i]
        sipms2 = sipm_s2(sipm, s2df)
        if np.sum(sipms2).ene_pes > thr:
            SIPM[j] = sipms2
            j+=1
    return pd.Panel(SIPM)

def sipm_hit_index(sipmp):
    """
    Store the indexes of the (sipm number)
    """
    hi =[]
    for i in sipmp.items:
        sipm = sipmp[i]
        hi.append(sipm.indx.values[0])
    return pd.Series(hi)

def sipmps2p_energy(sipms2p):
    """
    Takes a sipms2p as input
    Returns a DataFrame with the index and the energy in the SiPM as columns:

    """

    SIPM=[]
    for i in sipms2p.items:
        swf = {}
        swf['ene_pes'] = np.sum(sipms2p[i].ene_pes.values)
        swf['indx'] = sipms2p[i].indx.values[0]
        SIPM.append(swf)

    return pd.DataFrame(SIPM)

def sipm_s2_panel(sipmrwf, SIPMDF, s2df, thr_min=0.5, thr_s2 =1, event_number=0):
    """
    Takes the sipmrwf and a s2df
    Returns a sipm panel with a collection of sipm DF such that:
    1. the range of the sipm is specified by s2
    2. the sipm energy are above threshold.
    """

    sipmwf = sipmrwf[event_number]
    SIPM = {}
    NSIPM = sipmwf.shape[0]
    sipmwl = sipmwf.shape[1]

    j=0

    for i in range(NSIPM):

        adc_to_pes = SIPMDF.adc_to_pes[i]
        energy_pes = sipmwf[i]/adc_to_pes
        if np.sum(energy_pes) < thr_min:  #only worry about SiPM with energy above threshold
            continue

        time_mus = np.array(range(sipmwl))*1e+3 #steps are mus
        indx = np.ones(sipmwl)*i
        sipm = wf_mus(wfdf(time_mus,energy_pes,indx))
        sipms2 = sipm_s2(sipm, s2df)
        if np.sum(sipms2).ene_pes > thr_s2:
            SIPM[j] = sipms2
            j+=1
    return pd.Panel(SIPM)
