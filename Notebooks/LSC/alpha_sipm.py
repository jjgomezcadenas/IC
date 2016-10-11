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

def sipm_s2_panel(sipmp, s2, thr_min=0.5, thr_s2 =1, event_number=0):
    """
    Takes a sipmp
    Returns a sipm panel with a collection of sipm DF such that:
    1. the range of the sipm is specified by s2
    2. the sipm energy are above threshold.
    """

    j=0
    SIPM={}
    ESIPM=[]
    for i in (sipmp.items):
        sipm = sipmp[i]
        ESIPM.append(np.sum(sipm.ene_pes))

        if np.sum(sipm.ene_pes) < thr_min:  #only worry about SiPM with energy above threshold
            continue

        sipms2 = sipm_s2(sipm, s2)
        if np.sum(sipms2).ene_pes > thr_s2:
            SIPM[j] = sipms2
            j+=1
    return pd.Panel(SIPM), np.array(ESIPM)
    
def sipm_s2_energy(sipmp, s2, thr_min=0.5, thr_s2 =0.5):
    """
    Takes a sipmp
    Returns a sipm panel where every member is a SiPM which has energy above trheshold:
    The energy of the SiPM corresponds to the S2 width

    """

    j=0
    SIPM={}
    for i in (sipmp.items):
        sipm = sipmp[i]
        etot = np.sum(sipm.ene_pes)
        if etot < thr_min:  #only worry about SiPM with energy above threshold
            continue

        sipms2 = sipm_s2(sipm, s2)
        es2 = np.sum(sipms2).ene_pes
        if es2 > thr_s2:
            DATA = []
            DATA.append(i)
            DATA.append(etot)
            DATA.append(es2)
            SIPM[j] = DATA
            j+=1
    return pd.DataFrame(data=SIPM.values(), index=SIPM.keys(), columns=['sipm_indx','etot_pes','es2_pes'])
