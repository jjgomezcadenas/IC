"""
MCWFControlPlots
JMB August 2016

What MCWFControlPlots does:
1) Read a file containing MC waveforms (25 ns)
2) Energy histogram
3) Plot real&measured energy
4) Correlation plots between coordinates and energy
"""

from __future__ import print_function

from Util import *

import sys
import tables as tb

import Configure as CF

from MCUtil import *
from PlotUtil import *

def mcwfcontrolplots(CFP):
    PATH_IN=CFP['PATH_IN']
    FILE_IN=CFP['FILE_IN']
    FIRST_EVT=CFP['FIRST_EVT']
    LAST_EVT=CFP['LAST_EVT']

    h5in = tb.open_file(PATH_IN+'/'+FILE_IN, "r+")
    pmtcwf_ = h5in.root.RD.pmtcwf
    sipmrwf_ = h5in.root.RD.sipmrwf
    mctracks_ = h5in.root.MC.MCTracks

    NEVENTS = pmtcwf_.shape[0]
    print("Number of events in file: {}".format(NEVENTS))
    
    #Print some values
    print("PMT energies for first 10 events:\n")
    for i in xrange(0, 10):
        print(energy_pes(i,pmtcwf_))
        
    #Compute pmt & tracks total energies
    pmtEnergies = map(lambda i: energy_pes(i, pmtcwf_).sum(), xrange(NEVENTS))
    totalEnergies = energy_tracks(mctracks_)

    tkrEnergies, xs, ys, zs = energyPositionTracks(mctracks_)
    
    #Plot histogram
    HSimple1(totalEnergies,50,'Events energy','E (MeV)', 'Events')
    
    #Energy correlation
    energyCorrelation(totalEnergies,pmtEnergies)
    
    #Position correlation
    positionCorrelation(pmtEnergies,xs,ys,zs)
    positionCorrelation3d(pmtEnergies,xs,ys,zs)
    
    h5in.close()

if __name__ == '__main__':
    INFO, CFP = CF.configure(sys.argv[0], sys.argv[1:])
    mcwfcontrolplots(CFP)
