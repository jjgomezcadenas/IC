
from Centella.AAlgo import AAlgo
from Centella.physical_constants import *
from Centella.system_of_units import *
from Util import *
import numpy as np
import tables 
from time import time

from ROOT import gSystem
gSystem.Load("$GATE_DIR/lib/libGATE")

from ROOT import gate

from math import sqrt

class DetectorGeometry(tables.IsDescription):
    """
    Stores geometry information for the detector
    """
    x_det = tables.Float64Col(shape=2, pos=1) #xmin, xmax
    y_det = tables.Float64Col(shape=2, pos=2) #ymin, ymax
    z_det = tables.Float64Col(shape=2, pos=3) #zmin, zmax
    r_det = tables.Float64Col(pos=4) # radius


class DataPMT(tables.IsDescription):
    """
    Stores metadata information for the PMTs
    (position, gain, calibration-constant, mask)
    """
    channel = tables.Int16Col(pos=1) #electronic channel
    active = tables.Int16Col(pos=2) # 1 if active. 0 if dead
    position = tables.Float64Col(shape=3, pos=3)
    gain =tables.Float64Col(pos=4)
    adc_to_pes =tables.Float64Col(pos=5)

class DataSiPM(tables.IsDescription):
    """
    Stores metadata information for the SiPMs
    (position, gain, calibration-constant, mask)
    """
    channel = tables.Int16Col(pos=1) #electronic channel
    active = tables.Int16Col(pos=2) # 1 if active. 0 if dead
    position = tables.Float64Col(shape=3, pos=3)
    gain =tables.Float64Col(pos=4)
    adc_to_pes =tables.Float64Col(pos=5)

class MCTrack(tables.IsDescription):
    """
    Stores the parameters used by the simulation as metadata
    using Pytables
    """
    event_indx = tables.Int16Col(pos=1) 
    mctrk_indx = tables.Int16Col(pos=2) 
    particle_name = tables.StringCol(10,pos=3)  #displaces the baseline (e.g, 700)
    pdg_code = tables.Int16Col(pos=4)   # number of PMTs (12) 
    initial_vertex =tables.Float64Col(shape=3, pos=5)
    final_vertex =tables.Float64Col(shape=3, pos=6)
    momentum =tables.Float64Col(shape=3, pos=7)
    energy =tables.Float64Col(pos=8)
    nof_hits = tables.Int16Col(pos=9) 
    hit_indx = tables.Int16Col(pos=10)
    hit_position = tables.Float64Col(shape=3, pos=11)
    hit_time =tables.Float64Col(pos=12)
    hit_energy =tables.Float64Col(pos=13)

class MCWF(AAlgo):

    def __init__(self,param=False,level = 1,label="",**kargs):

        """
        """
        
        self.name='MCWaveform'
        
        AAlgo.__init__(self,param,level,self.name,0,label,kargs)
        

    def initialize(self):

        """    
        tables.EArray(parentnode, name, atom=None, shape=None, 
        title='', filters=None, expectedrows=None, 
        chunkshape=None, byteorder=None, _log=True)[source]    
        """
        
        self.m.log(1,'+++Init method of MCWaveform algorithm+++')
        self.NPMTS = int(12) 
        self.LEN_PMT = int(599999)
        self.NSIPM = int(1792) 
        self.LEN_SIPM = int(600)
        path="/Users/jjgomezcadenas/Documents/Development/NEXT/data/Waveforms/"
        file="WF_Tl_0.h5"
        
        self.NEVENTS = self.logman["CNTJob"].ints["NEVENTS"] 
        self.h5f = tables.open_file(path+file, "w",
                                    filters=tables.Filters(complib="blosc", complevel=9))

        self.pmtrd = self.h5f.create_earray(self.h5f.root, "pmtrd", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, self.NPMTS, self.LEN_PMT), 
                                    expectedrows=self.NEVENTS)

        self.sipmrd = self.h5f.create_earray(self.h5f.root, "sipmrd", 
                                    atom=tables.IntAtom(), 
                                    shape=(0, self.NSIPM, self.LEN_SIPM), 
                                    expectedrows=self.NEVENTS)

        group = self.h5f.create_group(self.h5f.root, "Detector")
        self.geom_table = self.h5f.create_table(group, "DetectorGeometry", DetectorGeometry,
                          "DetectorGeometry",
                           tables.Filters(0))

        group = self.h5f.create_group(self.h5f.root, "Sensors")
        self.pmt_table = self.h5f.create_table(group, "DataPMT", DataPMT,
                          "DataPMT",
                           tables.Filters(0))
        self.sipm_table = self.h5f.create_table(group, "DataSiPM", DataSiPM,
                          "DataSiPM",
                           tables.Filters(0))

        group = self.h5f.create_group(self.h5f.root, "MC")
        self.MCTrack_table = self.h5f.create_table(group, "MCTracks", MCTrack, "MCTracks",
                           tables.Filters(0))
        # self.h5f = tables.open_file("pmtrd.h5", "a", 
        #                     filters=tables.Filters(complib="blosc", complevel=9))
        #self.pmt = self.h5f.create_table(self.h5f.root, "pmt", PMTRD)
            
        self.t = time()
        self.n_evt = 0
        
        return

    def execute(self,event=""):

        """            
        """

        t1 = time()
        print "event id = %d, time =%7f.3 sec"%(event.GetEventID(),t1 -self.t)

        if self.n_evt == 0:
            print "instantiating geometry"
            gm = self.run.GetGeometry()
    
            print """
            min detector x coordinate = {}
            max detector x coordinate = {}
            min detector y coordinate = {}
            max detector y coordinate = {}
            min detector z coordinate = {}
            max detector z coordinate = {}
            max detector radius = {}
            """.format(gm.GetXmin(),gm.GetXmax(),gm.GetYmin(),gm.GetYmax(),
            gm.GetZmin(),gm.GetZmax(),gm.GetRmax())

            geom_ = self.geom_table.row
            x = np.array([gm.GetXmin(),gm.GetXmax()])
            y = np.array([gm.GetYmin(),gm.GetYmax()])
            z = np.array([gm.GetZmin(),gm.GetZmax()])
            geom_["x_det"] = x
            geom_["y_det"] = y
            geom_["z_det"] = z
            geom_["r_det"] = gm.GetRmax()
            geom_.append()

        
            pmt_id =[]
            sipm_id =[]
            for sensor in self.run.GetGeometry().GetSensors():
                if sensor and sensor.GetID()<1000:  
                    pmt_id.append(sensor.GetID())
                elif sensor and sensor.GetID() >1000:
                    sipm_id.append(sensor.GetID())

            print " pmt id = {}".format(pmt_id)
            print " sipm id = {}".format(sipm_id)

            pmt_ = self.pmt_table.row
            sipm_ = self.sipm_table.row


            for i in pmt_id:
                pmt =  gm.GetSensor(i)    
                print """
                PMTs
                id = {}
                pos = {}, {}, {}
                """.format(pmt.GetID(),
                    pmt.GetPosition().x(), pmt.GetPosition().y(), 
                    pmt.GetPosition().z())

                pos =  np.array([pmt.GetPosition().x(),pmt.GetPosition().y(),
                                 pmt.GetPosition().z()])
                pmt_["channel"] = pmt.GetID()
                pmt_["active"] = 1 # change if PMT becomes dead
                pmt_["position"] = pos
                pmt_["gain"] = 4.5e+6 #nominal value, change as needed.
                pmt_["adc_to_pes"] = 20 #nominal value, change as needed
                pmt_.append()
        

            for i in sipm_id:
                pmt =  gm.GetSensor(i)    
                print """
                SiPMs
                id = {}
                pos = {}, {}, {}
                """.format(pmt.GetID(),
                    pmt.GetPosition().x(), pmt.GetPosition().y(), 
                    pmt.GetPosition().z())

                pos =  np.array([pmt.GetPosition().x(),pmt.GetPosition().y(),
                                 pmt.GetPosition().z()])

                sipm_["channel"] = pmt.GetID()
                sipm_["active"] = 1 # change if SiPM becomes dead
                sipm_["position"] = pos
                sipm_["gain"] = 1 #nominal value, change as needed.
                sipm_["adc_to_pes"] = 1 #nominal value, change as needed
                sipm_.append()
            wait()


        
        mctrk_ = self.MCTrack_table.row
        vtx =np.zeros(3,dtype=np.float64)
        mctrks = event.GetMCTracks()
        it = 0
        for mctrk in mctrks:
            print "mctrk indx ={}".format(it)

            #get MCParticle generating the track
            pt =  mctrk.GetParticle() 

            print """
                track number {}
                PDG code = {}
                particle name = {}
                initial vertex xyz = {}, {}, {}
                final vertex xyz = {}, {}, {}
                initial momentum = {}, {}, {}
                final momentum = {}, {}, {}

            """.format(it+1,pt.GetPDG(), 
                       pt.GetLabel(),
                       pt.GetInitialVtx().x(),pt.GetInitialVtx().y(),pt.GetInitialVtx().z(),
                       pt.GetFinalVtx().x(),pt.GetFinalVtx().y(),pt.GetFinalVtx().z(),
                       pt.GetInitialMom().x(),pt.GetInitialMom().y(),pt.GetInitialMom().z(),
                       pt.GetFinalMom().x(),pt.GetFinalMom().y(),pt.GetFinalMom().z())

            #get track energy
            ene = mctrk.GetEnergy()
            print "energy of track = {}".format(ene)

            #const std::vector<gate::BHit*>& MCTrack::GetHits() 
            hits = mctrk.GetHits()
            ih = 0
            for hit in hits:
                print "event indx = {} mctrk indx ={} hit indx = {}".format(self.n_evt,it,ih)
                mctrk_["event_indx"] = self.n_evt
                mctrk_["mctrk_indx"] = it
                vtx[0] = pt.GetInitialVtx().x()
                vtx[1] = pt.GetInitialVtx().y()
                vtx[2] = pt.GetInitialVtx().z()
            
                mctrk_["particle_name"] = pt.GetLabel()
                mctrk_["pdg_code"] =  pt.GetPDG()
                mctrk_["initial_vertex"] = vtx

                vtx[0] = pt.GetFinalVtx().x()
                vtx[1] = pt.GetFinalVtx().y()
                vtx[2] = pt.GetFinalVtx().z()

                mctrk_["final_vertex"] = vtx

                vtx[0] = pt.GetInitialMom().x()
                vtx[1] = pt.GetInitialMom().y()
                vtx[2] = pt.GetInitialMom().z()

                mctrk_["momentum"] = vtx
                mctrk_["energy"] = mctrk.GetEnergy()
                vtx[0] = hit.GetPosition().x()
                vtx[1] = hit.GetPosition().y()
                vtx[2] = hit.GetPosition().z()
                mctrk_["nof_hits"] = len(hits)
                mctrk_["hit_indx"] = ih
                mctrk_["hit_position"] = vtx
                mctrk_["hit_time"] = hit.GetTime()
                mctrk_["hit_energy"] = hit.GetAmplitude()

                print """
                hit number = {}
                hit id = {}
                hit position: x = {} y = {} z = {}
                hit time = {}
                hit energy = {}
            """.format(
                ih, hit.GetID(),
                hit.GetPosition().x(),hit.GetPosition().y(),hit.GetPosition().z(),
                hit.GetTime(),
                hit.GetAmplitude()
                )
                mctrk_.append()
                ih+=1
            it+=1 

        sipms = event.GetHits(gate.SIPM)
        print "number of sipms = {}".format(len(sipms))
        #wait()

        rdata = []
        for i in range(len(sipms)): 
            #print "SIPM number ",i
            

            wf = sipms[i].GetWaveform().GetData()
            wa = np.zeros(len(wf), dtype=np.int16)
            
            if i==0:
                print "SIPM: length of buffer ={}".format(len(wf))
            for j in range(len(wf)):
                wa[j] = wf[j].second
                
            rdata.append(wa)

        self.sipmrd.append(np.array(rdata).reshape(1, self.NSIPM, self.LEN_SIPM))

        #wait()


        
        pmts = event.GetHits(gate.PMT)
        print "len(pmts) = ", len(pmts)

        # len_signal = len(pmts[0].GetWaveform().GetData())
        # print "len signal =", len_signal
        

        # # ht = np.arange(0.0, len_signal*1 +0, 1, dtype=np.int)
        # print(self.h5f)
        #wait()
        
        rdata = []
        for i in range(self.NPMTS):
        #for i in range(1): 
            print "PMT number ",i
            
            wf = pmts[i].GetWaveform().GetData()
            wa = np.zeros(self.LEN_PMT, dtype=np.int8)
            
            if i==0:
                print "PMT: length of buffer ={}".format(len(wf))
            for j in range(len(wf)):
                wa[j] = wf[j].second
                
            rdata.append(wa)

        self.pmtrd.append(np.array(rdata).reshape(1, self.NPMTS, self.LEN_PMT))
        
        self.n_evt +=1
        return True
        
    def finalize(self):
        

        self.m.log(1,'+++End method of MCWaveform algorithm+++')
        self.pmtrd.flush()
        self.sipmrd.flush()
        self.geom_table.flush()
        self.pmt_table.flush()
        self.sipm_table.flush()
        self.MCTrack_table.flush()
        self.h5f.close()

        return

    