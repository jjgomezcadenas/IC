from ROOT import *
import MySQLdb
from os import mkdir
from sys import argv
import numpy as np

# TODO: This should be read somehow from the file
NEWMAXRUN = 2904

dbIC = MySQLdb.connect(host="localhost", user='USER',passwd='PASS',db="ICNEWDB")
cursorIC = dbIC.cursor()

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
cursorIC.execute('select max(MinRun) from SipmNoiseBins;')
maxRunIC = cursorIC.fetchall()[0][0]
sqlfile.write('UPDATE SipmNoiseBins SET MaxRun={0} where MinRun={1} and MaxRun is NULL;\n'.format(NEWMAXRUN-1,maxRunIC))

energies = map(Hdummy.GetBinCenter, range(1,Hdummy.GetNbinsX()+1))
for idx,e in enumerate(energies):
    sqlfile.write('INSERT INTO SipmNoiseBins Values ({0},NULL,{1},{2});\n'.format(NEWMAXRUN,idx,e))

# Normalize and set probabilities
cursorIC.execute('select max(MinRun) from SipmNoise;')
maxRunIC = cursorIC.fetchall()[0][0]
sqlfile.write('UPDATE SipmNoise SET MaxRun={0} where MinRun={1} and MaxRun is NULL;\n'.format(NEWMAXRUN-1,maxRunIC))

nbins = Hdummy.GetNbinsX()
for ID,H in sorted(HS.items()):
    if H.GetNbinsX() != nbins:
        print "ERROR, nbins in {} is {}, different from {}".format(ID,H.GetNbinsX(),nbins)
        break
    data = np.array(map( H.GetBinContent, range(1,H.GetNbinsX()+1) ))
    datanorm = data/np.sum(data) if data.any() else data
    values = ', '.join(map(str,datanorm))
    sqlfile.write('INSERT INTO SipmNoise Values ({0},NULL,{1},{2});\n'.format(NEWMAXRUN,ID,values))

sqlfile.close()
