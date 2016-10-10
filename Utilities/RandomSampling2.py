'''
    Defines a class for random sampling
'''
from __future__ import print_function
import numpy as np
import scipy.stats as st
import scipy.interpolate as interpolator
from time import time


class NoiseSiPM(st.rv_continuous):
    def SetCDF(self,x,y):
        x =  np.array(x) * 1.0
        y = (np.array(y) * 1.0).cumsum()
        self.d = np.diff(x)[0]
        self.y = y/y[-1]
        self.a = min(x) - self.d/2
        self.b = max(x)/self.d + 0.5
        self.cdf_ = lambda x: self.y[int((x-self.a)/self.d)]
        print(self.a,self.b,self.d)

    def _cdf(self,x):
        return self.cdf_(x)

    def RVS(self,*args,**kwargs):
        return self.rvs(*args,**kwargs) * self.d

class SiPMNoiseSampler:
    def __init__(self,filename,index_map):
        print("Initializing SiPMNoiseSampler",end=" ")
        data = np.loadtxt(filename)
        self.xbins = data[0,1:]
        self.IDs   = data[1:,0]
        self.probs = { ID : data[1+i,1:] for i,ID in enumerate(self.IDs) }

        xmin = min(self.xbins)
        xmax = max(self.xbins)
        self.RVSs = [ NoiseSiPM(a=xmin,b=xmax,name=str(ID)) for ID in index_map.keys() ]
        for ID,index in index_map.items():
            self.RVSs[index].SetCDF(self.xbins,self.probs[ID])
        print("DONE")

    def Sample( self, size = 1 ):
        return np.array( [ pdf.RVS(size=size) for pdf in self.RVSs ] )

if __name__ == '__main__':
    x = np.array(range(10))
    y = np.array((range(10)))

    a = NoiseSiPM(name = 'oli')
    a.SetCDF(x,y)

    print(a.cdf(8))
    print(a.rvs())
