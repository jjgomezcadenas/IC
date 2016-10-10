def pmt_alpha(pmtrwf,pmtdf, geomdf, thr=7*pes, t_trigger = 600, log='INFO', plot=False, event_list=[0]):
    """
    alpha analysis based on PMTs
    """
    lg = 'logging.'+DEBUG
    logger.setLevel(eval(lg))

    evl = len(event_list)
    t0 = np.zeros(evl, dtype=np.float32)
    t = np.zeros(evl, dtype=np.float32)
    xb = np.zeros(evl, dtype=np.float32)
    yb = np.zeros(evl, dtype=np.float32)
    s2e = np.zeros(evl, dtype=np.float32)
    s2l = np.zeros(evl, dtype=np.float32)
    ns1 = np.zeros(evl, dtype=np.int32)
    ns2 = np.zeros(evl, dtype=np.int32)

    for event in event_list:
        logger.info('event = {}'.format(event))

        PMT, BSL  = waveform_panel(pmtrwf,pmtdf,event=event)
        if plot:
            plot_PPMT(PMT, tmin=0, tmax=1000, emin = -10, emax = 20, option='all')
            plt.show()
            wait()
            plot_PPMT(PMT, tmin=0, tmax=1000, emin = -10, emax = 150, option='sum')
            plt.show()
            wait()

        s12 = find_S12(wf_thr(sPMT(PMT),threshold=thr))

        logger.debug('length of s12 = {}'.format(len(s12)))
        S1 = []
        S2 = []
        for s in s12:
            logger.debug('evaluating s in s12: s ={}'.format(s.describe()))
            logger.debug('tmax ={}'.format(s.describe().time_mus.max()))
            if s.describe().time_mus.max() < t_trigger: #s1
                S1.append(s)
            else:
                S2.append(s)

        logger.debug('length of S1 list = {}'.format(len(S1)))
        logger.debug('length of S2 list = {}'.format(len(S2)))

        if (len(S1) == 0):
            logger.warning("S1 not found, ignore event")
            t0[event] = -999
            t[event] = -999
            xb[event] = -999
            yb[event] = -999
            s2e[event] = -999
            s2l[event] = -999
            ns1[event] = 0
            ns2[event] = len(S2)
            continue
        if (len(S2) == 0):
            logger.warning("S2 not found, ignore event")
            t0[event] = -999
            t[event] = -999
            xb[event] = -999
            yb[event] = -999
            s2e[event] = -999
            s2l[event] = -999
            ns1[event] = len(S1)
            ns2[event] = 0
            continue

        ns1[event] = len(S1)
        ns2[event] = len(S2)

        s1 = S1[0]

        if len(S1) > 1:
            cmax = 0
            i=0
            imax = 0
            for s in S1:
                if s.describe().time_mus.count() > cmax:
                    cmax = s.describe().time_mus.count()
                    imax = i
            i+=1
            s1 = S1[imax]

        logger.debug('found s1 = {}'.format(s1.describe()))



        s2 = S2[0]
        es2 = s12_energy(s2)

        if len(S2) > 1:
            emax = 0
            i=0
            imax = 0
            for s in S2:
                es2 = s12_energy(s)
                if es2 > emax:
                    emax =  es2
                    imax = i
            i+=1
            s2 = S2[imax]

        epmt = energy_sum(PMT, thr=0)

        logger.debug('found s2 = {}'.format(s2.describe()))

        t0[event] = find_t0(s1)
        t[event] = find_t(s1,s2)
        xb[event], yb[event] = pmt_barycenter(pmtdf, epmt)
        s2l[event] = s12_length(s2)


def get_vectors(h5f):
    """
    input: file pointer
    returns: data vectors
    """
    pmtrwf = h5f.root.RD.pmtrwf
    sipmrwf = h5f.root.RD.sipmrwf
    geom_t = h5f.root.Detector.DetectorGeometry
    pmt_t = h5f.root.Sensors.DataPMT
    sipm_t = h5f.root.Sensors.DataSiPM
    gdf = snf.read_data_geom(geom_t)
    pmtdf = snf.read_data_sensors(pmt_t)
    sipmdf = snf.read_data_sensors(sipm_t)
    return pmtrwf,sipmrwf,pmtdf,sipmdf,gdf

def get_pmt_vectors(h5f):
    """
    input: file pointer
    returns: data vectors
    """
    pmtrwf = h5f.root.RD.pmtrwf
    geom_t = h5f.root.Detector.DetectorGeometry
    pmt_t = h5f.root.Sensors.DataPMT
    gdf = snf.read_data_geom(geom_t)
    pmtdf = snf.read_data_sensors(pmt_t)
    return pmtrwf,pmtdf,gdf

def wfdf(time,energy_pes,indx):
    """
    takes three vectors (time, energy and indx) and returns a data frame representing a waveform
    """
    swf = {}
    swf['time_mus'] = time/mus
    swf['ene_pes'] = energy_pes
    swf['indx'] = indx
    return pd.DataFrame(swf)

def waveform_panel(pmtrwf,pmtdf,mau_len = 500, calib_constat =True, adc_to_pes=20,
                   type = 'PMT', daq_ceiling=4096, event=0):
    """
    input: sensor (pmt or sipm) data vector, sensor data frame (position, calibration)
    returns: a panel holding waveforms for all sensors, and a series for the baselines
    """
    PMT = {}
    nm = mau_len
    B_MAU = (1./nm)*np.ones(nm)
    pmt_len = pmtrwf.shape[2]
    NPMT = pmtrwf.shape[1]
    MAU = np.zeros(nm)
    BSL = {}

    time_ns = np.arange(pmt_len)*mus
    indx = np.arange(pmt_len)

    if type == 'PMT':
        time_ns = np.arange(pmt_len)*FP.time_DAQ

    ene_sum = 0
    for j in range(NPMT):

        if calib_constat == True:
            adc_to_pes = abs(pmtdf['adc_to_pes'][j])

        signal_daq = pmtrwf[event,j]
        if type == 'PMT':
            signal_daq = daq_ceiling - pmtrwf[event,j]

        MAU[0:nm] = SGN.lfilter(B_MAU,1, signal_daq[0:nm])
        BASELINE = MAU[nm-1]

        ene_pes = (signal_daq - BASELINE)/adc_to_pes
        if type == 'PMT':
            ene_sum += ene_pes

        PMT[j] = wfdf(time_ns,ene_pes,indx)
        BSL[j] = BASELINE
    PMT[j+1] = wfdf(time_ns,ene_sum,indx)
    return pd.Panel(PMT),pd.Series(BSL)
def plot_PPMT(pmt_panel, tmin=0, tmax=1200, emin = 0, emax = 10000, option='sum'):
    """
    Plots pmtwf
    """
    plt.figure(figsize=(10,10))

    if option == 'sum':
        ax1 = plt.subplot(1,1,1)
        ax1.set_xlim([tmin, tmax])
        ax1.set_ylim([emin, emax])
        indx = pmt_panel.items[-1]
        pmtwf = pmt_panel[indx]
        plt.plot(pmtwf['time_mus'],pmtwf['ene_pes'])
    else:

        for i in pmt_panel.items[0:-1]:
            ax1 = plt.subplot(3,4,int(i)+1)
            ax1.set_xlim([tmin, tmax])
            ax1.set_ylim([emin, emax])

            pmtwf = pmt_panel[i]
            plt.plot(pmtwf['time_mus'],pmtwf['ene_pes'])

    plt.show()
def wf_thr(wf,threshold=0):
    """
    return a zero supressed waveform (more generally, the vaules of wf above threshold)
    """
    return wf.loc[lambda df: df.ene_pes.values >threshold, :]
def energy_sum(sensor_panel, thr=0):
    """
    Sum the WFs of PMTs and SiPMs (MC) and store the total energy in PES
    """
    EPES = []

    for i in sensor_panel.items[0:-1]:
        pmtwf = sensor_panel[i]
        EPES.append(np.sum(pmtwf.ene_pes.values[np.where(pmtwf.ene_pes.values>thr)]))
    return pd.Series(EPES)
def plot_sensors(geom_df,sensor_df, energy, radius=10):
    """
    plots the energy of the sensors
    """
    x =sensor_df['x'].values
    y =sensor_df['y'].values
    r =np.ones(len(sensor_df['x'].values))*radius

    plt.figure(figsize=(10,10))
    ax = plt.subplot(aspect='equal')
    mpl.circles(x, y, r, c=energy, alpha=0.5, ec='none')
    plt.colorbar()

    plt.xlim(geom_df['xdet_min'],geom_df['xdet_max'])
    plt.ylim(geom_df['ydet_min'],geom_df['ydet_max'])
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
    returns t0
    """
    emax = np.amax(s1.ene_pes.values)
    return s1.loc[lambda df: df.ene_pes.values ==emax, :]
def s12_energy(s12):
    """
    total energy in pes
    """
    return np.sum(s12.ene_pes.values)
def s12_length(s12):
    """
    s2 length in mus
    """

    return s12.describe().time_mus['max'] - s12.describe().time_mus['min']
def s12_peak(s2):
    """
    s2 peak in mus
    """

    return s12.describe().time_mus['max'], s2.describe().ene_pes['max']
def find_t(s1,s2):
    """
    returns the time of the interaction
    """
    t0 = find_t0(s1).time_mus.values[0]
    ts2,es2 = s12_peak(s2)
    return ts2 - t0
