'''
    Defines a class for random sampling
'''
from __future__ import print_function
import numpy as np

class NoiseSampler:
    def __init__(self,filename,sipmdf,sample_size = 1,smear = True):
        '''
            filename    -> (string) txt file containing a matrix of noise distributions
            sipmdf      -> (data frame) with the SiPMs info
            sample_size -> (int) number of samples per call
            smear       -> (bool) perform continuous sampling
        '''
        print("Initializing NoiseSampler...",end=" ")

        # Read data, take xbins and compute (half of) bin size.
        data  = np.loadtxt(filename)
        self.xbins = data[0,1:]
        self.dx    = np.diff(self.xbins)[0] * 0.5

        # Remove masked channels and normalize probabilities
        data  = data[ np.where(map(sipmdf['channel'].values.__contains__,data[:,0])) ]
        self.probs = np.apply_along_axis( lambda ps: ps/np.sum(ps), 1, data[:,1:] )

        self.nsamples = sample_size
        self.nsensors = len(self.probs)
        self.output_shape = (self.nsensors,self.nsamples)

        # Sampling functions
        self._sample_sensor      = lambda probs: np.random.choice( self.xbins, size=self.nsamples, p=probs )
        self._discrete_sampler   = lambda: np.apply_along_axis( self._sample_sensor, 1, self.probs )
        self._continuous_sampler = lambda: self._discrete_sampler() + np.random.uniform(-self.dx,self.dx)

        self._sampler = self._continuous_sampler if smear else self._discrete_sampler
        print("OK")

    def Sample( self ):
        return self._sampler()

    def ComputeThresholds( self, noise_cut = 0.99, **kwargs ):
        '''
            Find the number of pes at which each noise distribution leaves behind
            the noise_cut fraction of its population. Options passed through kwargs:
            - pes_to_adc : float or nsensor-sized array of floats
            - sipmdf     : data frame with the adc_to_pes conversion
            - if neither of the above is present, the thresholds are given in pes.
        '''
        # If any of the options is present, perform pes-to-adc conversion
        if 'sipmdf' in kwargs:
            pes_to_adc = kwargs['sipmdf']['adc_to_pes']
        elif 'pes_to_adc' in kwargs:
            pes_to_adc = kwargs['pes_to_adc']
            if not hasattr(pes_to_adc,__iter__):
                pes_to_adc *= np.ones(self.nsensors)

        return np.array( [ self.xbins[np.argwhere( probs > noise_cut )[0][0]] for i,probs in enumerate(np.apply_along_axis( np.cumsum, 1, self.probs )) ] ) * np.array(pes_to_adc)
