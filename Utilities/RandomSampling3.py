'''
    Defines a class for random sampling
'''
from __future__ import print_function
import numpy as np
from time import time

class SiPMNoiseSampler:
    def __init__(self,filename,index_map):
        print("Initializing SiPMNoiseSampler",end=" ")
        data = np.loadtxt(filename)
        self.xbins = data[0,1:]
        self.IDs   = data[1:,0]
        self.d     = np.diff(self.xbins)[0]*0.5
        self.probs = [None]*len(self.IDs)
        for i,ID in enumerate(self.IDs):
            index = index_map.get(int(ID),None)
            if index is None: continue
            self.probs[index] = data[1+i,1:]/sum(data[1+i,1:])
        while None in self.probs:
            self.probs.pop(self.probs.index(None))
        print("DONE")

    def Sample( self, size = 1, smear = True ):
        x = np.array( [ np.random.choice(self.xbins, size=size, p=prob) for prob in self.probs ] )
        return x if not smear else x + np.random.uniform(-self.d,self.d,size=x.shape[0])
