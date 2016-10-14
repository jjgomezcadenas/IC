'''
    Defines a class for random sampling
'''
from __future__ import print_function
import numpy as np

class NoiseSampler:
    def __init__(self,filename,sipmdf,sample_size = 1,smear = True, pes2adc = None):
        '''
            filename    -> (string) txt file containing a matrix of noise distributions
            sipmdf      -> (data frame) with the SiPMs info
            sample_size -> (int) number of samples per call
            smear       -> (bool) perform continuous sampling
            pes2adc     -> (float or array) overrule adc conversion from sipmdf
        '''
        print("Initializing NoiseSampler...",end=" ")

        # Read data, backup xbins, remove masked and xbins from matrix
        # and normalize probabilities
        data  = np.loadtxt(filename)
        xbins = data[0,1:]
        data  = data[ np.where(map(sipmdf['channel'].__contains__,data[:,0]))[0] ]
        probs = np.apply_along_axis( lambda ps: ps/np.sum(ps), 1, data[:,1:] )

        self.nsamples = sample_size
        self.nsensors = len(probs)
        self.output_shape = (self.nsensors,self.nsamples)

        # Convert xbins to adc and compute a matrix with half of the bin size
        if pes2adc is None:
            pes2adc = 1.0/sipmdf['adc_to_pes']
        elif not hasattr(pes2adc,__iter__):
            pes2adc = np.ones(self.nsensors) * pes2adc

        xbins   = np.tile(xbins,(self.nsensors,1)) * np.array(pes2adc).reshape(self.nsensors,1)
        self.dx = np.tile(np.diff(self.xbins,axis=1)[:,0]*0.5, self.nsamples)

        # data array has shape (nsensors,2*self.nsamples) for efficiency
        self.data  = np.concatenate((self.xbins,self.probs),axis=1)

        # Sampling functions
        self._sample_sensor      = lambda data: np.random.choice( data[:data.shape/2], size = self.nsamples, p=data[:data.shape/2] )
        self._discrete_sampler   = lambda: np.apply_along_axis( self._sample_sensor, 1, self.data )
        self._continuous_sampler = lambda: self._discrete_sampler() + np.random.uniform(-self.dx,self.dx)

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
