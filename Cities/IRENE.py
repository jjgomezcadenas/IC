
from __future__ import print_function
import numpy as np
import tables
from time import time

import tables as tb
import numpy as np
import Core.mplFunctions as mpl
import Core.wfmFunctions as wfm
from Database import loadDB
import ICython.Sierpe.BLR as blr
import ICython.Core.peakFunctions as cpf
import Core.peakFunctions as pf
from ICython.Core.system_of_units import SystemOfUnits

#import Core.system_of_units as units
from Core.LogConfig import logger
from Core.Configure import configure, define_event_loop, print_configuration
import Core.tblFunctions as tbl

units = SystemOfUnits()

class S12(tb.IsDescription):
    """
    Store for a S1/S2
    The table maps a S12:
    peak is the index of the S12 dictionary, running over the number of peaks found
    time and energy of the peak
    """
    event = tb.Int32Col(pos=0)
    peak = tb.UInt8Col(pos=1)  # peak number
    time = tb.Float32Col(pos=2) # time in ns
    ene = tb.Float32Col(pos=3) # energy in pes

class S2Si(tb.IsDescription):
    """
    Store for a S2Si
    The table maps a S2Si
    peak is the same than the S2 peak
    nsipm gives the SiPM number
    only energies are stored (times are defined in S2)
    """
    event = tb.Int32Col(pos=0)
    peak = tb.UInt8Col(pos=1)  # peak number
    nsipm = tb.Int16Col(pos=2)  # sipm number
    nsample = tb.Int16Col(pos=3) # sample number
    ene = tb.Float32Col(pos=4) # energy in pes


class Irene:
    """
    The city of IRENE performs a fast processing directly
    from raw data (pmtrwf and sipmrwf) to PMAPS.
    It is optimized for speed (use of CYTHON functions) and intended
    for fast processing of data

    """
    def __init__(self, run_number):
        """
        Inits the machine with the run number
        loads the data base to access calibration and geometry
        sets all switches to default value (False most of the time)
        """
        DataPMT = loadDB.DataPMT(run_number)
        DataSiPM = loadDB.DataSiPM(run_number)

        self.adc_to_pes = abs(DataPMT.adc_to_pes.values).astype(np.double)
        self.sipm_adc_to_pes = DataSiPM.adc_to_pes.values.astype(np.double)
        self.coeff_c = DataPMT.coeff_c.values.astype(np.double)
        self.coeff_blr = DataPMT.coeff_blr.values.astype(np.double)
        self.xs = DataSiPM.X.values
        self.ys = DataSiPM.Y.values

        self.setFiles = False
        self.setPmapStore = False

        # BLR default values (override with set_BLR)
        self.n_baseline = 28000
        self.thr_trigger = 5*units.adc

        # MAU default values (override with set_MAU)
        self.n_MAU =  100
        self.thr_MAU = thr_MAU=3*units.adc

        # CSUM default values (override with set_CSUM)
        self.thr_csum= 1*units.pes

        self.setS1 = False
        self.setS2 = False
        self.setSiPM = False

        self.plot_csum = False
        self.plot_s1 = False
        self.plot_s2 = False
        self.plot_sipm = False
        self.plot_sipmzs = False
        self.plot_simap = False

        self.nprint = 1000000


    def set_plot(self, plot_csum=False, plot_s1=False, plot_s2=False,
                 plot_sipm=False, plot_sipmzs=False, plot_simap=False):
        """
        decides what to plot
        """

        self.plot_csum = plot_csum
        self.plot_s1 = plot_s1
        self.plot_s2 = plot_s2
        self.plot_sipm = plot_sipm
        self.plot_sipmzs = plot_sipmzs
        self.plot_simap = plot_simap

    def set_print(self, nprint=10):
        """
        print frequency
        """
        self.nprint = nprint

    def set_input_files(self,path, input_files):
        """
        Sets the input files
        """
        self.path = path
        self.input_files = input_files
        self.setFiles = True

    def set_pmap_store(self, path, pmap_file, compression='ZLIB4'):
        """
        Sets the input files
        """
        filename = path + pmap_file
        self.pmapFile = tb.open_file(filename, "w",
                          filters=tbl.filters(compression))

        pmapsgroup = self.pmapFile.create_group(self.pmapFile.root,
                                                "PMAPS")

        # create tables to store pmaps
        self.s1t = self.pmapFile.create_table(pmapsgroup,
                                                 "S1", S12,
                                                 "S1 Table",
                                                 tbl.filters(compression))

        self.s2t = self.pmapFile.create_table(pmapsgroup,
                                                 "S2", S12,
                                                 "S2 Table",
                                                 tbl.filters(compression))

        self.s2sit = self.pmapFile.create_table(pmapsgroup,
                                                 "S2Si", S2Si,
                                                 "S2Si Table",
                                                 tbl.filters(compression))
        self.s1t.cols.event.create_index()
        self.s2t.cols.event.create_index()
        self.s2sit.cols.event.create_index()

        self.setPmapStore = True

    def set_BLR(self, n_baseline=38000, thr_trigger=5*units.adc):
        """
        Parameters of the BLR
        """
        self.n_baseline = n_baseline
        self.thr_trigger = thr_trigger

    def set_MAU(self, n_MAU=100, thr_MAU=3*units.adc):
        """
        Parameters of the MAU used to remove low frequency noise
        """
        self.n_MAU =  n_MAU
        self.thr_MAU = thr_MAU

    def set_CSUM(self, thr_csum=1*units.pes):
        """
        Parameter for ZS in the calibrated sum
        """
        self.thr_csum= thr_csum

    def set_S1(self, tmin=0*units.mus, tmax=590*units.mus,
               stride=4, lmin=4, lmax=20):
        """
        Parameters for S1 search
        """
        self.tmin_s1 = tmin
        self.tmax_s1 = tmax
        self.stride_s1 = stride
        self.lmin_s1 = lmin
        self.lmax_s1 = lmax
        self.setS1 = True

    def set_S2(self, tmin=590*units.mus, tmax=620*units.mus,
               stride=40, lmin=100, lmax=1000000):
        """
        Parameters for S2 search
        """
        self.tmin_s2 = tmin
        self.tmax_s2 = tmax
        self.stride_s2 = stride
        self.lmin_s2 = lmin
        self.lmax_s2 = lmax
        self.setS2 = True

    def set_SiPM(self, thr_zs=5*units.pes, thr_sipm_s2=50*units.pes):
        """
        Parameters for SiPM analysis
        """
        self.thr_zs = thr_zs
        self.thr_sipm_s2 = thr_sipm_s2
        self.setSiPM = True


    def plot_ene_sipm(self,S2Si, radius=3):
        """
        plots the reconstructed energy of the SiPMs
        input: sipm dictionary
        """

        r = np.ones(len(self.xs)) * radius
        col = np.zeros(len(self.xs))
        for i in S2Si.keys():
            sipml = S2Si[i]
            for sipm in sipml:
                sipm_n = sipm[0]
                sipm_wf = sipm[1]
                sipm_e = np.sum(sipm_wf)
                col[sipm_n] = sipm_e

        plt.figure(figsize=(10, 10))
        plt.subplot(aspect="equal")
        mpl.circles(self.xs, self.ys, r, c=col, alpha=0.5, ec="none")
        plt.colorbar()

        plt.xlim(-198, 198)
        plt.ylim(-198, 198)

    def store_pmaps(self, evt, S1, S2, S2Si):
        """
        Store PMAPS
        """
        #S1
        row = self.s1t.row

        for i in S1.keys():
            time = S1[i][0]
            ene = S1[i][1]
            assert(len(time) == len(ene))
            for j in range(len(time)):
                row["event"] = evt
                row["peak"] = i
                row["time"] = time[j]
                row["ene"] = ene[j]
                row.append()

        self.s1t.flush()

        #S2
        row = self.s2t.row

        for i in S2.keys():
            time = S2[i][0]
            ene = S2[i][1]
            assert(len(time) == len(ene))
            for j in range(len(time)):
                row["event"] = evt
                row["peak"] = i
                row["time"] = time[j]
                row["ene"] = ene[j]
                row.append()

        self.s2t.flush()

        # S2
        row = self.s2sit.row

        for i in S2Si.keys():
            sipml = S2Si[i]
            for sipm in sipml:
                nsipm = sipm[0]
                ene = sipm[1]
                for j in range(len(ene)):
                    if ene[j] > 0:
                        row["event"] = evt
                        row["peak"] = i
                        row["nsipm"] = nsipm
                        row["nsample"] = j
                        row["ene"] = ene[j]
                        row.append()

        self.s2t.flush()


    def run(self, nmax, store_pmaps=False):
        """
        Run the machine
        nmax is the max number of events to run
        store_pmaps decides whether to store pmaps or not
        """
        n_events_tot = 0

        if self.setFiles == False:
            raise IOError('must set files before running')
        if self.path =='':
            raise IOError('path is empty')
        if len(self.input_files) == 0:
            raise IOError('input file list is empty')


        if self.setS1 == False:
            raise IOError('must set S1 parameters before running')
        if self.setS2 == False:
            raise IOError('must set S2 parameters before running')
        if self.setSiPM == False:
            raise IOError('must set Sipm parameters before running')


        print("""
                 IRENE will run a max of {} events
                 Storing PMAPS (1=yes/0=no)
                 Input path ={}
                 Input Files ={}
                          """.format(nmax, store_pmaps,
                                     self.path,self.input_files))

        print("""
                 S1 parameters
                 tmin = {} mus tmax = {} mus stride = {}
                 lmin = {} lmax = {}
                          """.format(self.tmin_s1/units.mus,
                                     self.tmax_s1/units.mus,
                                     self.stride_s1,
                                     self.lmin_s1,
                                     self.lmax_s1))
        print("""
                 S2 parameters
                 tmin = {} mus tmax = {} mus stride = {}
                 lmin = {} lmax = {}
                          """.format(self.tmin_s2/units.mus,
                                     self.tmax_s2/units.mus,
                                     self.stride_s2,
                                     self.lmin_s2,
                                     self.lmax_s2))
        print("""
                 S2Si parameters
                 threshold min charge per SiPM = {} pes
                 threshold min charge in S2 = {} pes
                          """.format(self.thr_zs,
                                     self.thr_sipm_s2))

        # loop over input files
        first=False
        for ffile in self.input_files:

            print("Opening", ffile, end="... ")
            filename = self.path + ffile

            try:
                with tb.open_file(filename, "r+") as h5in:
                    pmtrwf = h5in.root.RD.pmtrwf
                    sipmrwf = h5in.root.RD.sipmrwf

                    if first == False:

                        self.NEVT, self.NPMT, self.PMTWL = pmtrwf.shape
                        NEVT, self.NSIPM, self.SIPMWL = sipmrwf.shape

                        print("""
                        Number of events in file = {}
                        Number of PMTs = {}
                        PMTWL = {}
                        Number of SiPMs = {}
                        SiPMWL = {}
                          """.format(self.NEVT, self.NPMT, self.PMTWL,
                                     self.NSIPM, self.SIPMWL))

                        self.signal_t = np.arange(0., self.PMTWL * 25, 25)
                        first = True

                    for evt in range(self.NEVT):
                        # deconvolve
                        CWF = blr.deconv_pmt(pmtrwf[evt],
                                             self.coeff_c, self.coeff_blr,
                                             n_baseline=self.n_baseline,
                                             thr_trigger=self.thr_trigger)

                        # calibrated PMT sum
                        csum = cpf.calibrated_pmt_sum(CWF,
                                                      self.adc_to_pes,
                                                      n_MAU=self.n_MAU,
                                                      thr_MAU=self.thr_MAU)
                        if self.plot_csum:
                            mpl.plot_signal(self.signal_t/units.mus,
                                            csum, title="calibrated sum, ZS",
                                            signal_start=0, signal_end=1200,
                                            ymax = 200,
                                            t_units='mus', units="pes")

                            plt.show()
                            raw_input('->')


                        # Supress samples below threshold (in pes)
                        wfzs_ene, wfzs_indx = cpf.wfzs(csum,
                                                       threshold=self.thr_csum)

                        #wfzs_t = cpf.time_from_index(wfzs_indx)

                        # find S1 and S2
                        S1 = cpf.find_S12(wfzs_ene, wfzs_indx,
                                         tmin=self.tmin_s1,
                                         tmax=self.tmax_s1,
                                         lmin=self.lmin_s1, lmax=self.lmax_s1,
                                         stride=self.stride_s1,
                                         rebin=False)

                        S2 = cpf.find_S12(wfzs_ene, wfzs_indx,
                                         tmin=self.tmin_s2,
                                         tmax=self.tmax_s2,
                                         lmin=self.lmin_s2, lmax=self.lmax_s2,
                                         stride=self.stride_s2,
                                         rebin=True,
                                         rebin_stride=self.stride_s2)

                        if self.plot_s1:
                            pf.scan_S12(S1)
                        if self.plot_s2:
                            pf.scan_S12(S2)

                        # sipms

                        if self.plot_sipm:
                            mpl.plot_sipm(sipmrwf[evt], nmin=0, nmax=16, x=4, y=4)

                        sipmzs = cpf.signal_sipm(sipmrwf[evt],
                                                 self.sipm_adc_to_pes,
                                                 thr=self.thr_zs,
                                                 n_MAU=self.n_MAU)
                        if self.plot_sipmzs:
                            mpl.plot_sipm(sipmzs, nmin=0, nmax=16, x=4, y=4)
                            plt.show()
                            raw_input('->')

                        SIPM = cpf.select_sipm(sipmzs)
                        S2Si = pf.sipm_S2_dict(SIPM, S2, thr=self.thr_sipm_s2)

                        if store_pmaps == True:
                            if self.setPmapStore == False:
                                raise IOError('must set PMAPS before storing')
                            self.store_pmaps(n_events_tot, S1,S2,S2Si)

                        if self.plot_simap:
                            self.plot_ene_sipm(sipmd)
                            plt.show()
                            raw_input('->')

                        n_events_tot +=1
                        if n_events_tot%self.nprint == 0:
                            print('event in file = {}, total = {}'.\
                              format(evt, n_events_tot))

                        if n_events_tot >= nmax:
                            print('reached maximum number of events (={})'.format(nmax))
                            return nmax


            except:
                print('error')
                raise

        if store_pmaps == True:
            self.pmapFile.close()
        return n_events_tot

if __name__ == "__main__":
    input_files =['run_3112.gdcsnext.000.next1el_3112.root.h5',
'run_3112.gdcsnext.001.next1el_3112.root.h5',
'run_3112.gdcsnext.002.next1el_3112.root.h5',
'run_3112.gdcsnext.003.next1el_3112.root.h5',
'run_3112.gdcsnext.004.next1el_3112.root.h5',
'run_3112.gdcsnext.005.next1el_3112.root.h5',
'run_3112.gdcsnext.006.next1el_3112.root.h5',
'run_3112.gdcsnext.007.next1el_3112.root.h5',
'run_3112.gdcsnext.008.next1el_3112.root.h5',
'run_3112.gdcsnext.009.next1el_3112.root.h5',
'run_3112.gdcsnext.010.next1el_3112.root.h5',
'run_3112.gdcsnext.011.next1el_3112.root.h5',
'run_3112.gdcsnext.012.next1el_3112.root.h5',
'run_3112.gdcsnext.013.next1el_3112.root.h5',
'run_3112.gdcsnext.014.next1el_3112.root.h5',
'run_3112.gdcsnext.015.next1el_3112.root.h5',
'run_3112.gdcsnext.016.next1el_3112.root.h5',
'run_3112.gdcsnext.017.next1el_3112.root.h5',
'run_3112.gdcsnext.018.next1el_3112.root.h5',
'run_3112.gdcsnext.019.next1el_3112.root.h5',
'run_3112.gdcsnext.020.next1el_3112.root.h5',
'run_3112.gdcsnext.021.next1el_3112.root.h5',
'run_3112.gdcsnext.022.next1el_3112.root.h5',
'run_3112.gdcsnext.023.next1el_3112.root.h5',
'run_3112.gdcsnext.024.next1el_3112.root.h5',
'run_3112.gdcsnext.025.next1el_3112.root.h5',
'run_3112.gdcsnext.026.next1el_3112.root.h5',
'run_3112.gdcsnext.027.next1el_3112.root.h5',
'run_3112.gdcsnext.028.next1el_3112.root.h5',
'run_3112.gdcsnext.029.next1el_3112.root.h5',
'run_3112.gdcsnext.030.next1el_3112.root.h5',
'run_3112.gdcsnext.031.next1el_3112.root.h5',
'run_3112.gdcsnext.032.next1el_3112.root.h5',
'run_3112.gdcsnext.033.next1el_3112.root.h5',
'run_3112.gdcsnext.034.next1el_3112.root.h5',
'run_3112.gdcsnext.035.next1el_3112.root.h5',
'run_3112.gdcsnext.036.next1el_3112.root.h5',
'run_3112.gdcsnext.037.next1el_3112.root.h5',
'run_3112.gdcsnext.038.next1el_3112.root.h5',
'run_3112.gdcsnext.039.next1el_3112.root.h5',
'run_3112.gdcsnext.040.next1el_3112.root.h5',
'run_3112.gdcsnext.041.next1el_3112.root.h5',
'run_3112.gdcsnext.042.next1el_3112.root.h5',
'run_3112.gdcsnext.043.next1el_3112.root.h5',
'run_3112.gdcsnext.044.next1el_3112.root.h5',
'run_3112.gdcsnext.045.next1el_3112.root.h5',
'run_3112.gdcsnext.046.next1el_3112.root.h5',
'run_3112.gdcsnext.047.next1el_3112.root.h5',
'run_3112.gdcsnext.048.next1el_3112.root.h5',
'run_3112.gdcsnext.049.next1el_3112.root.h5',
'run_3112.gdcsnext.050.next1el_3112.root.h5',
'run_3112.gdcsnext.051.next1el_3112.root.h5',
'run_3112.gdcsnext.052.next1el_3112.root.h5',
'run_3112.gdcsnext.053.next1el_3112.root.h5',
'run_3112.gdcsnext.054.next1el_3112.root.h5',
'run_3112.gdcsnext.055.next1el_3112.root.h5',
'run_3112.gdcsnext.056.next1el_3112.root.h5',
'run_3112.gdcsnext.057.next1el_3112.root.h5',
'run_3112.gdcsnext.058.next1el_3112.root.h5',
'run_3112.gdcsnext.059.next1el_3112.root.h5',
'run_3112.gdcsnext.060.next1el_3112.root.h5',
'run_3112.gdcsnext.061.next1el_3112.root.h5',
'run_3112.gdcsnext.062.next1el_3112.root.h5',
'run_3112.gdcsnext.063.next1el_3112.root.h5',
'run_3112.gdcsnext.064.next1el_3112.root.h5',
'run_3112.gdcsnext.065.next1el_3112.root.h5',
'run_3112.gdcsnext.066.next1el_3112.root.h5',
'run_3112.gdcsnext.067.next1el_3112.root.h5',
'run_3112.gdcsnext.068.next1el_3112.root.h5',
'run_3112.gdcsnext.069.next1el_3112.root.h5',
'run_3112.gdcsnext.070.next1el_3112.root.h5',
'run_3112.gdcsnext.071.next1el_3112.root.h5',
'run_3112.gdcsnext.072.next1el_3112.root.h5',
'run_3112.gdcsnext.073.next1el_3112.root.h5',
'run_3112.gdcsnext.074.next1el_3112.root.h5',
'run_3112.gdcsnext.075.next1el_3112.root.h5',
'run_3112.gdcsnext.076.next1el_3112.root.h5',
'run_3112.gdcsnext.077.next1el_3112.root.h5',
'run_3112.gdcsnext.078.next1el_3112.root.h5',
'run_3112.gdcsnext.079.next1el_3112.root.h5',
'run_3112.gdcsnext.080.next1el_3112.root.h5',
'run_3112.gdcsnext.081.next1el_3112.root.h5',
'run_3112.gdcsnext.082.next1el_3112.root.h5',
'run_3112.gdcsnext.083.next1el_3112.root.h5',
'run_3112.gdcsnext.084.next1el_3112.root.h5',
'run_3112.gdcsnext.085.next1el_3112.root.h5',
'run_3112.gdcsnext.086.next1el_3112.root.h5',
'run_3112.gdcsnext.087.next1el_3112.root.h5',
'run_3112.gdcsnext.088.next1el_3112.root.h5',
'run_3112.gdcsnext.089.next1el_3112.root.h5',
'run_3112.gdcsnext.090.next1el_3112.root.h5',
'run_3112.gdcsnext.091.next1el_3112.root.h5',
'run_3112.gdcsnext.092.next1el_3112.root.h5',
'run_3112.gdcsnext.093.next1el_3112.root.h5',
'run_3112.gdcsnext.094.next1el_3112.root.h5',
'run_3112.gdcsnext.095.next1el_3112.root.h5',
'run_3112.gdcsnext.096.next1el_3112.root.h5',
'run_3112.gdcsnext.097.next1el_3112.root.h5',
'run_3112.gdcsnext.098.next1el_3112.root.h5',
'run_3112.gdcsnext.099.next1el_3112.root.h5',
'run_3112.gdcsnext.100.next1el_3112.root.h5',
'run_3112.gdcsnext.101.next1el_3112.root.h5',
'run_3112.gdcsnext.102.next1el_3112.root.h5',
'run_3112.gdcsnext.103.next1el_3112.root.h5',
'run_3112.gdcsnext.104.next1el_3112.root.h5',
'run_3112.gdcsnext.105.next1el_3112.root.h5',
'run_3112.gdcsnext.106.next1el_3112.root.h5',
'run_3112.gdcsnext.107.next1el_3112.root.h5',
'run_3112.gdcsnext.108.next1el_3112.root.h5',
'run_3112.gdcsnext.109.next1el_3112.root.h5',
'run_3112.gdcsnext.110.next1el_3112.root.h5',
'run_3112.gdcsnext.111.next1el_3112.root.h5',
'run_3112.gdcsnext.112.next1el_3112.root.h5',
'run_3112.gdcsnext.113.next1el_3112.root.h5',
'run_3112.gdcsnext.114.next1el_3112.root.h5',
'run_3112.gdcsnext.115.next1el_3112.root.h5',
'run_3112.gdcsnext.116.next1el_3112.root.h5',
'run_3112.gdcsnext.117.next1el_3112.root.h5',
'run_3112.gdcsnext.118.next1el_3112.root.h5',
'run_3112.gdcsnext.119.next1el_3112.root.h5',
'run_3112.gdcsnext.120.next1el_3112.root.h5',
'run_3112.gdcsnext.121.next1el_3112.root.h5',
'run_3112.gdcsnext.122.next1el_3112.root.h5',
'run_3112.gdcsnext.123.next1el_3112.root.h5',
'run_3112.gdcsnext.124.next1el_3112.root.h5',
'run_3112.gdcsnext.125.next1el_3112.root.h5',
'run_3112.gdcsnext.126.next1el_3112.root.h5',
'run_3112.gdcsnext.127.next1el_3112.root.h5',
'run_3112.gdcsnext.128.next1el_3112.root.h5',
'run_3112.gdcsnext.129.next1el_3112.root.h5',
'run_3112.gdcsnext.130.next1el_3112.root.h5',
'run_3112.gdcsnext.131.next1el_3112.root.h5',
'run_3112.gdcsnext.132.next1el_3112.root.h5',
'run_3112.gdcsnext.133.next1el_3112.root.h5',
'run_3112.gdcsnext.134.next1el_3112.root.h5',
'run_3112.gdcsnext.135.next1el_3112.root.h5',
'run_3112.gdcsnext.136.next1el_3112.root.h5',
'run_3112.gdcsnext.137.next1el_3112.root.h5',
'run_3112.gdcsnext.138.next1el_3112.root.h5',
'run_3112.gdcsnext.139.next1el_3112.root.h5',
'run_3112.gdcsnext.140.next1el_3112.root.h5',
'run_3112.gdcsnext.141.next1el_3112.root.h5',
'run_3112.gdcsnext.142.next1el_3112.root.h5',
'run_3112.gdcsnext.143.next1el_3112.root.h5',
'run_3112.gdcsnext.144.next1el_3112.root.h5',
'run_3112.gdcsnext.145.next1el_3112.root.h5',
'run_3112.gdcsnext.146.next1el_3112.root.h5',
'run_3112.gdcsnext.147.next1el_3112.root.h5',
'run_3112.gdcsnext.148.next1el_3112.root.h5',
'run_3112.gdcsnext.149.next1el_3112.root.h5',
'run_3112.gdcsnext.150.next1el_3112.root.h5',
'run_3112.gdcsnext.151.next1el_3112.root.h5',
'run_3112.gdcsnext.152.next1el_3112.root.h5',
'run_3112.gdcsnext.153.next1el_3112.root.h5',
'run_3112.gdcsnext.154.next1el_3112.root.h5',
'run_3112.gdcsnext.155.next1el_3112.root.h5',
'run_3112.gdcsnext.156.next1el_3112.root.h5',
'run_3112.gdcsnext.157.next1el_3112.root.h5',
'run_3112.gdcsnext.158.next1el_3112.root.h5',
'run_3112.gdcsnext.159.next1el_3112.root.h5',
'run_3112.gdcsnext.160.next1el_3112.root.h5',
'run_3112.gdcsnext.161.next1el_3112.root.h5',
'run_3112.gdcsnext.162.next1el_3112.root.h5',
'run_3112.gdcsnext.163.next1el_3112.root.h5',
'run_3112.gdcsnext.164.next1el_3112.root.h5',
'run_3112.gdcsnext.165.next1el_3112.root.h5',
'run_3112.gdcsnext.166.next1el_3112.root.h5',
'run_3112.gdcsnext.167.next1el_3112.root.h5',
'run_3112.gdcsnext.168.next1el_3112.root.h5',
'run_3112.gdcsnext.169.next1el_3112.root.h5',
'run_3112.gdcsnext.170.next1el_3112.root.h5',
'run_3112.gdcsnext.171.next1el_3112.root.h5',
'run_3112.gdcsnext.172.next1el_3112.root.h5',
'run_3112.gdcsnext.173.next1el_3112.root.h5',
'run_3112.gdcsnext.174.next1el_3112.root.h5',
'run_3112.gdcsnext.175.next1el_3112.root.h5',
'run_3112.gdcsnext.176.next1el_3112.root.h5',
'run_3112.gdcsnext.177.next1el_3112.root.h5',
'run_3112.gdcsnext.178.next1el_3112.root.h5',
'run_3112.gdcsnext.179.next1el_3112.root.h5',
'run_3112.gdcsnext.180.next1el_3112.root.h5',
'run_3112.gdcsnext.181.next1el_3112.root.h5',
'run_3112.gdcsnext.182.next1el_3112.root.h5',
'run_3112.gdcsnext.183.next1el_3112.root.h5',
'run_3112.gdcsnext.184.next1el_3112.root.h5',
'run_3112.gdcsnext.185.next1el_3112.root.h5',
'run_3112.gdcsnext.186.next1el_3112.root.h5',
'run_3112.gdcsnext.187.next1el_3112.root.h5',
'run_3112.gdcsnext.188.next1el_3112.root.h5',
'run_3112.gdcsnext.189.next1el_3112.root.h5',
'run_3112.gdcsnext.190.next1el_3112.root.h5',
'run_3112.gdcsnext.191.next1el_3112.root.h5',
'run_3112.gdcsnext.192.next1el_3112.root.h5',
'run_3112.gdcsnext.193.next1el_3112.root.h5',
'run_3112.gdcsnext.194.next1el_3112.root.h5',
'run_3112.gdcsnext.195.next1el_3112.root.h5',
'run_3112.gdcsnext.196.next1el_3112.root.h5',
'run_3112.gdcsnext.197.next1el_3112.root.h5',
'run_3112.gdcsnext.198.next1el_3112.root.h5',
'run_3112.gdcsnext.199.next1el_3112.root.h5',
'run_3112.gdcsnext.200.next1el_3112.root.h5',
'run_3112.gdcsnext.201.next1el_3112.root.h5',
'run_3112.gdcsnext.202.next1el_3112.root.h5',
'run_3112.gdcsnext.203.next1el_3112.root.h5',
'run_3112.gdcsnext.204.next1el_3112.root.h5',
'run_3112.gdcsnext.205.next1el_3112.root.h5',
'run_3112.gdcsnext.206.next1el_3112.root.h5',
'run_3112.gdcsnext.207.next1el_3112.root.h5',
'run_3112.gdcsnext.208.next1el_3112.root.h5',
'run_3112.gdcsnext.209.next1el_3112.root.h5',
'run_3112.gdcsnext.210.next1el_3112.root.h5',
'run_3112.gdcsnext.211.next1el_3112.root.h5',
'run_3112.gdcsnext.212.next1el_3112.root.h5',
'run_3112.gdcsnext.213.next1el_3112.root.h5',
'run_3112.gdcsnext.214.next1el_3112.root.h5',
'run_3112.gdcsnext.215.next1el_3112.root.h5',
'run_3112.gdcsnext.216.next1el_3112.root.h5',
'run_3112.gdcsnext.217.next1el_3112.root.h5',
'run_3112.gdcsnext.218.next1el_3112.root.h5',
'run_3112.gdcsnext.219.next1el_3112.root.h5',
'run_3112.gdcsnext.220.next1el_3112.root.h5',
'run_3112.gdcsnext.221.next1el_3112.root.h5',
'run_3112.gdcsnext.222.next1el_3112.root.h5',
'run_3112.gdcsnext.223.next1el_3112.root.h5',
'run_3112.gdcsnext.224.next1el_3112.root.h5',
'run_3112.gdcsnext.225.next1el_3112.root.h5',
'run_3112.gdcsnext.226.next1el_3112.root.h5',
'run_3112.gdcsnext.227.next1el_3112.root.h5',
'run_3112.gdcsnext.228.next1el_3112.root.h5',
'run_3112.gdcsnext.229.next1el_3112.root.h5',
'run_3112.gdcsnext.230.next1el_3112.root.h5',
'run_3112.gdcsnext.231.next1el_3112.root.h5',
'run_3112.gdcsnext.232.next1el_3112.root.h5',
'run_3112.gdcsnext.233.next1el_3112.root.h5',
'run_3112.gdcsnext.234.next1el_3112.root.h5',
'run_3112.gdcsnext.235.next1el_3112.root.h5',
'run_3112.gdcsnext.236.next1el_3112.root.h5',
'run_3112.gdcsnext.237.next1el_3112.root.h5',
'run_3112.gdcsnext.238.next1el_3112.root.h5',
'run_3112.gdcsnext.239.next1el_3112.root.h5',
'run_3112.gdcsnext.240.next1el_3112.root.h5',
'run_3112.gdcsnext.241.next1el_3112.root.h5',
'run_3112.gdcsnext.242.next1el_3112.root.h5',
'run_3112.gdcsnext.243.next1el_3112.root.h5',
'run_3112.gdcsnext.244.next1el_3112.root.h5',
'run_3112.gdcsnext.245.next1el_3112.root.h5',
'run_3112.gdcsnext.246.next1el_3112.root.h5',
'run_3112.gdcsnext.247.next1el_3112.root.h5',
'run_3112.gdcsnext.248.next1el_3112.root.h5',
'run_3112.gdcsnext.249.next1el_3112.root.h5',
'run_3112.gdcsnext.250.next1el_3112.root.h5',
'run_3112.gdcsnext.251.next1el_3112.root.h5',
'run_3112.gdcsnext.252.next1el_3112.root.h5',
'run_3112.gdcsnext.253.next1el_3112.root.h5',
'run_3112.gdcsnext.254.next1el_3112.root.h5',
'run_3112.gdcsnext.255.next1el_3112.root.h5',
'run_3112.gdcsnext.256.next1el_3112.root.h5',
'run_3112.gdcsnext.257.next1el_3112.root.h5',
'run_3112.gdcsnext.258.next1el_3112.root.h5',
'run_3112.gdcsnext.259.next1el_3112.root.h5',
'run_3112.gdcsnext.260.next1el_3112.root.h5',
'run_3112.gdcsnext.261.next1el_3112.root.h5',
'run_3112.gdcsnext.262.next1el_3112.root.h5',
'run_3112.gdcsnext.263.next1el_3112.root.h5',
'run_3112.gdcsnext.264.next1el_3112.root.h5',
'run_3112.gdcsnext.265.next1el_3112.root.h5',
'run_3112.gdcsnext.266.next1el_3112.root.h5',
'run_3112.gdcsnext.267.next1el_3112.root.h5',
'run_3112.gdcsnext.268.next1el_3112.root.h5',
'run_3112.gdcsnext.269.next1el_3112.root.h5',
'run_3112.gdcsnext.270.next1el_3112.root.h5',
'run_3112.gdcsnext.271.next1el_3112.root.h5',
'run_3112.gdcsnext.272.next1el_3112.root.h5',
'run_3112.gdcsnext.273.next1el_3112.root.h5',
'run_3112.gdcsnext.274.next1el_3112.root.h5',
'run_3112.gdcsnext.275.next1el_3112.root.h5',
'run_3112.gdcsnext.276.next1el_3112.root.h5',
'run_3112.gdcsnext.277.next1el_3112.root.h5',
'run_3112.gdcsnext.278.next1el_3112.root.h5',
'run_3112.gdcsnext.279.next1el_3112.root.h5',
'run_3112.gdcsnext.280.next1el_3112.root.h5',
'run_3112.gdcsnext.281.next1el_3112.root.h5',
'run_3112.gdcsnext.282.next1el_3112.root.h5',
'run_3112.gdcsnext.283.next1el_3112.root.h5',
'run_3112.gdcsnext.284.next1el_3112.root.h5',
'run_3112.gdcsnext.285.next1el_3112.root.h5',
'run_3112.gdcsnext.286.next1el_3112.root.h5',
'run_3112.gdcsnext.287.next1el_3112.root.h5',
'run_3112.gdcsnext.288.next1el_3112.root.h5',
'run_3112.gdcsnext.289.next1el_3112.root.h5',
'run_3112.gdcsnext.290.next1el_3112.root.h5',
'run_3112.gdcsnext.291.next1el_3112.root.h5',
'run_3112.gdcsnext.292.next1el_3112.root.h5',
'run_3112.gdcsnext.293.next1el_3112.root.h5',
'run_3112.gdcsnext.294.next1el_3112.root.h5',
'run_3112.gdcsnext.295.next1el_3112.root.h5',
'run_3112.gdcsnext.296.next1el_3112.root.h5',
'run_3112.gdcsnext.297.next1el_3112.root.h5',
'run_3112.gdcsnext.298.next1el_3112.root.h5',
'run_3112.gdcsnext.299.next1el_3112.root.h5',
'run_3112.gdcsnext.300.next1el_3112.root.h5',
'run_3112.gdcsnext.301.next1el_3112.root.h5',
'run_3112.gdcsnext.302.next1el_3112.root.h5',
'run_3112.gdcsnext.303.next1el_3112.root.h5',
'run_3112.gdcsnext.304.next1el_3112.root.h5',
'run_3112.gdcsnext.305.next1el_3112.root.h5',
'run_3112.gdcsnext.306.next1el_3112.root.h5',
'run_3112.gdcsnext.307.next1el_3112.root.h5',
'run_3112.gdcsnext.308.next1el_3112.root.h5',
'run_3112.gdcsnext.309.next1el_3112.root.h5',
'run_3112.gdcsnext.310.next1el_3112.root.h5',
'run_3112.gdcsnext.311.next1el_3112.root.h5',
'run_3112.gdcsnext.312.next1el_3112.root.h5',
'run_3112.gdcsnext.313.next1el_3112.root.h5',
'run_3112.gdcsnext.314.next1el_3112.root.h5',
'run_3112.gdcsnext.315.next1el_3112.root.h5',
'run_3112.gdcsnext.316.next1el_3112.root.h5',
'run_3112.gdcsnext.317.next1el_3112.root.h5',
'run_3112.gdcsnext.318.next1el_3112.root.h5',
'run_3112.gdcsnext.319.next1el_3112.root.h5',
'run_3112.gdcsnext.320.next1el_3112.root.h5',
'run_3112.gdcsnext.321.next1el_3112.root.h5',
'run_3112.gdcsnext.322.next1el_3112.root.h5',
'run_3112.gdcsnext.323.next1el_3112.root.h5',
'run_3112.gdcsnext.324.next1el_3112.root.h5',
'run_3112.gdcsnext.325.next1el_3112.root.h5',
'run_3112.gdcsnext.326.next1el_3112.root.h5',
'run_3112.gdcsnext.327.next1el_3112.root.h5',
'run_3112.gdcsnext.328.next1el_3112.root.h5',
'run_3112.gdcsnext.329.next1el_3112.root.h5',
'run_3112.gdcsnext.330.next1el_3112.root.h5',
'run_3112.gdcsnext.331.next1el_3112.root.h5',
'run_3112.gdcsnext.332.next1el_3112.root.h5',
'run_3112.gdcsnext.333.next1el_3112.root.h5',
'run_3112.gdcsnext.334.next1el_3112.root.h5',
'run_3112.gdcsnext.335.next1el_3112.root.h5',
'run_3112.gdcsnext.336.next1el_3112.root.h5',
'run_3112.gdcsnext.337.next1el_3112.root.h5',
'run_3112.gdcsnext.338.next1el_3112.root.h5',
'run_3112.gdcsnext.339.next1el_3112.root.h5',
'run_3112.gdcsnext.340.next1el_3112.root.h5',
'run_3112.gdcsnext.341.next1el_3112.root.h5',
'run_3112.gdcsnext.342.next1el_3112.root.h5',
'run_3112.gdcsnext.343.next1el_3112.root.h5',
'run_3112.gdcsnext.344.next1el_3112.root.h5',
'run_3112.gdcsnext.345.next1el_3112.root.h5',
'run_3112.gdcsnext.346.next1el_3112.root.h5',
'run_3112.gdcsnext.347.next1el_3112.root.h5',
'run_3112.gdcsnext.348.next1el_3112.root.h5',
'run_3112.gdcsnext.349.next1el_3112.root.h5',
'run_3112.gdcsnext.350.next1el_3112.root.h5',
'run_3112.gdcsnext.351.next1el_3112.root.h5',
'run_3112.gdcsnext.352.next1el_3112.root.h5',
'run_3112.gdcsnext.353.next1el_3112.root.h5',
'run_3112.gdcsnext.354.next1el_3112.root.h5',
'run_3112.gdcsnext.355.next1el_3112.root.h5',
'run_3112.gdcsnext.356.next1el_3112.root.h5',
'run_3112.gdcsnext.357.next1el_3112.root.h5',
'run_3112.gdcsnext.358.next1el_3112.root.h5',
'run_3112.gdcsnext.359.next1el_3112.root.h5',
'run_3112.gdcsnext.360.next1el_3112.root.h5',
'run_3112.gdcsnext.361.next1el_3112.root.h5',
'run_3112.gdcsnext.362.next1el_3112.root.h5',
'run_3112.gdcsnext.363.next1el_3112.root.h5',
'run_3112.gdcsnext.364.next1el_3112.root.h5',
'run_3112.gdcsnext.365.next1el_3112.root.h5',
'run_3112.gdcsnext.366.next1el_3112.root.h5',
'run_3112.gdcsnext.367.next1el_3112.root.h5',
'run_3112.gdcsnext.368.next1el_3112.root.h5',
'run_3112.gdcsnext.369.next1el_3112.root.h5',
'run_3112.gdcsnext.370.next1el_3112.root.h5',
'run_3112.gdcsnext.371.next1el_3112.root.h5',
'run_3112.gdcsnext.372.next1el_3112.root.h5',
'run_3112.gdcsnext.373.next1el_3112.root.h5',
'run_3112.gdcsnext.374.next1el_3112.root.h5',
'run_3112.gdcsnext.375.next1el_3112.root.h5',
'run_3112.gdcsnext.376.next1el_3112.root.h5',
'run_3112.gdcsnext.377.next1el_3112.root.h5',
'run_3112.gdcsnext.378.next1el_3112.root.h5',
'run_3112.gdcsnext.379.next1el_3112.root.h5',
'run_3112.gdcsnext.380.next1el_3112.root.h5',
'run_3112.gdcsnext.381.next1el_3112.root.h5',
'run_3112.gdcsnext.382.next1el_3112.root.h5',
'run_3112.gdcsnext.383.next1el_3112.root.h5',
'run_3112.gdcsnext.384.next1el_3112.root.h5',
'run_3112.gdcsnext.385.next1el_3112.root.h5',
'run_3112.gdcsnext.386.next1el_3112.root.h5',
'run_3112.gdcsnext.387.next1el_3112.root.h5',
'run_3112.gdcsnext.388.next1el_3112.root.h5',
'run_3112.gdcsnext.389.next1el_3112.root.h5',
'run_3112.gdcsnext.390.next1el_3112.root.h5',
'run_3112.gdcsnext.391.next1el_3112.root.h5',
'run_3112.gdcsnext.392.next1el_3112.root.h5',
'run_3112.gdcsnext.393.next1el_3112.root.h5',
'run_3112.gdcsnext.394.next1el_3112.root.h5',
'run_3112.gdcsnext.395.next1el_3112.root.h5',
'run_3112.gdcsnext.396.next1el_3112.root.h5',
'run_3112.gdcsnext.397.next1el_3112.root.h5',
'run_3112.gdcsnext.398.next1el_3112.root.h5',
'run_3112.gdcsnext.399.next1el_3112.root.h5',
'run_3112.gdcsnext.400.next1el_3112.root.h5']
    path='/Users/jjgomezcadenas/Documents/Development/NEXT/icdata/LSC/run3112/'
    fpp = Irene(run_number=3112)

    # set the state of the machine
    fpp.set_input_files(path, input_files)
    fpp.set_pmap_store(path, 'pmaps_0_400.h5', compression='ZLIB4')
    fpp.set_print(nprint=10)

    fpp.set_BLR(n_baseline=38000, thr_trigger=5*units.adc)
    fpp.set_MAU(n_MAU=100, thr_MAU=3*units.adc)
    fpp.set_CSUM(thr_csum=1*units.pes)
    fpp.set_S1(tmin=10*units.mus, tmax=590*units.mus, stride=4, lmin=8, lmax=20)
    fpp.set_S2(tmin=590*units.mus, tmax=1190*units.mus, stride=40, lmin=100,
           lmax=100000)
    fpp.set_SiPM(thr_zs=5*units.pes, thr_sipm_s2=25*units.pes)

    t0 = time()
    nevt = fpp.run(nmax=100000, store_pmaps=True)
    t1 = time()
    dt = t1 - t0

    print("run {} evts in {} s, time/event = {}".format(nevt, dt, dt/nevt))
