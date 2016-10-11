'''
    Defines a class for random sampling
'''
from __future__ import print_function
import numpy as np

class SiPMsNoiseSampler:
    def __init__(self,filename,index_map,sample_size = 1,smear = True):
        print("Initializing SiPMNoiseSampler",end=" ")
        data = np.loadtxt(filename )

        self.xbins = data[0,1:] # xbins are stored in the first line. First value is dummy.
        self.dx    = np.diff(self.xbins)[0]*0.5 # half of the bin size
        data = data[ np.where(map(index_map.__contains__,data[:,0]))[0] ] #Remove masked channels

        self.probs = np.apply_along_axis( lambda probs: probs/np.sum(probs), 1, data[:,1:] ) # normalize probabilities
        self.nsamples = sample_size
        self.nsensors = len(self.probs)
        self.output_shape = (self.nsensors,self.nsamples)

        self._sample_sensor = lambda probs: np.random.choice( self.xbins, size = self.nsamples, p=probs )
        self._discrete_sampler = lambda: np.apply_along_axis( self._sample_sensor, 1, self.probs )
        self._continuous_sampler = lambda: self._discrete_sampler() + np.random.uniform(-self.dx,self.dx,size=self.output_shape)

        self._sampler = self._continuous_sampler if smear else self._discrete_sampler
        print("OK")

    def Sample( self ):
        return self._sampler()
