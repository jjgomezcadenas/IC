'''
    Defines a class for random sampling
'''
from __future__ import print_function
import numpy as np

class NoiseSampler:
    def __init__(self,filename,sipmdf,sample_size = 1,smear = True):
        print("Initializing NoiseSampler...",end=" ")
        data = np.loadtxt(filename )

        self.xbins = data[0,1:] # xbins are stored in the first line. First value is dummy.
        self.dx    = np.diff(self.xbins)[0]*0.5 # half of the bin size
        data = data[ np.where(map(sipmdf['channel'].__contains__,data[:,0]))[0] ] #Remove masked channels

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

    def ComputeThresholds( self, noise_cut = 0.99 ):
        '''
            Find the number of pes at which each noise distribution leaves behind
            the noise_cut fraction of its population.
        '''
        return np.array( [ self.xbins[np.argwhere( probs > noise_cut )[0][0]] for i,probs in enumerate(np.apply_along_axis( np.cumsum, 1, self.probs )) ] )
