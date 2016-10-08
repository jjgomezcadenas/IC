'''
    Defines a class for random sampling
'''

import numpy as np
import scipy.stats as st
from time import time


class NoiseSiPM(st.rv_continuous):
    def SetCDF(self,x,y):
        x =  np.array(x) * 1.0
        y = (np.array(y) * 1.0).cumsum()
        self.a = min(x)
        self.b = max(x)
        self.cdf_ = interpolator.interp1d(x,y/y[-1],kind='cubic')#,bounds_error=False)

    def _cdf(self,x):
        return self.cdf_(x)

class SiPMNoiseSampler:
    def __init__(self,filename,index_map):
        data = np.loadtxt(filename)
        self.xbins = data[0,1:]
        self.IDs   = data[1:,0]
        self.probs = { ID : data[1+i,1:] for i,ID in enumerate(self.IDs) }

        xmin = min(self.xbins)
        xmax = max(self.xbins)
        self.RVSs = [ NoiseSiPM(0,a=xmin,b=xmax,name=str(ID)) for ID in index_map.keys() ]
        for ID,index in index_map.items():
            self.RVSs[index].SetCDF(self.xbins,self.probs[ID])

    def Sample( self, size = 1 ):
        print("GO")
        self.t0 = time()
        for pdf in self.RVSs:
            self.t0 = time()
            x = pdf.rvs(size=size)
            print(time() - self.t0)
            # a = np.array( [ pdf.rvs   (size=size) for pdf in self.RVSs ] )
        print(time() - self.t0)
        return a
        print('ok')
        return a

if __name__ == '__main__':
    x = np.array(range(10))
    y = np.array((range(10)))

    a = NoiseSiPMPDF(name = 'oli')
    a.SetCDF(x,y)

    print a.cdf(8)
    print a.rvs()
