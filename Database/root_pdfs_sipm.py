from ROOT import *
from os import mkdir
from sys import argv
import numpy as np

def Usage():
    print '''Usage: python {} [rootfile]'''.format(__file__)
    exit()

if len(argv)<2: Usage()

# Open files
sqlfile = open(argv[1] + '.sql','w')
F = TFile(argv[1],'readonly')

# Get histos
HS = {}
for key in F.GetDirectory("PDF1").GetListOfKeys():
    ID = int(''.join( filter( str.isdigit, key.GetName() ) ))
    HS[ID] = key.ReadObj()

Hdummy = HS.values()[0]

# Set energies
sqlfile.write('TRUNCATE TABLE NoiseBins;\n')
energies = map(Hdummy.GetBinCenter, range(1,Hdummy.GetNbinsX()+1))
for idx,e in enumerate(energies):
    sqlfile.write('INSERT INTO NoiseBins Values ({0},{1});\n'.format(idx,e))

# Normalize and set probabilities
sqlfile.write('TRUNCATE TABLE NoiseSiPM;\n')
nbins = Hdummy.GetNbinsX()
for ID,H in sorted(HS.items()):
    if H.GetNbinsX() != nbins:
        print "ERROR, nbins in {} is {}, different from {}".format(ID,H.GetNbinsX(),nbins)
        break
    data = np.array(map( H.GetBinContent, range(1,H.GetNbinsX()+1) ))
    datanorm = data/np.sum(data) if data.any() else data
    values = ', '.join(map(str,datanorm))
    sqlfile.write('INSERT INTO NoiseSiPM Values ({});\n'.format(values))

sqlfile.close()
