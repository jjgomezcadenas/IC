"""
Defines a class for random sampling.
"""
from __future__ import print_function

import numpy as np


class NoiseSampler:
    def __init__(self, filename, sipmdf, sample_size=1, smear=True):
        """
        Samples a histogram as if it was a PDF.

        Parameters
        ----------
        filename : string
            Path to the txt file containing the matrix of noise distributions.
        sipmdf: pandas DataFrame
            Contains the SiPMs info in DF format.
        sample_size: int
            Number of samples per sensor and call.
        smear: bool
            Flag to choose between performing discrete or continuous sampling.

        Attributes
        ---------
        xbins : numpy.ndarray
            Contains the the bins centers in pes.
        dx: float
            Half of the bin size.
        probs: numpy.ndarray
            Matrix holding the probability for each sensor at each bin.
        nsamples: int
            Number of samples per sensor taken at each call.
        nsensors: int
            Number of sensors to simulate.
        output_shape: tuple of ints
            Shape of the output array
        """

        # Read data, take xbins and compute (half of) bin size.
        data = np.loadtxt(filename)
        self.xbins = data[0, 1:]
        self.dx = np.diff(self.xbins)[0] * 0.5

        # Remove masked channels and normalize probabilities
        self._df = sipmdf
        channels = self._df["channel"]
        data = data[np.where(map(channels.values.__contains__, data[:, 0]))]
        self.probs = np.apply_along_axis(lambda ps: ps/np.sum(ps),
                                         1, data[:, 1:])

        self.nsamples = sample_size
        self.nsensors = len(self.probs)
        self.output_shape = (self.nsensors, self.nsamples)

        # Sampling functions
        self._sample_sensor = lambda probs: np.random.choice(
                                            self.xbins,
                                            size=self.nsamples,
                                            p=probs)
        self._discrete_sampler = lambda: np.apply_along_axis(
                                         self._sample_sensor,
                                         1, self.probs)
        self._continuous_sampler = lambda: (self._discrete_sampler() +
                                            np.random.uniform(-self.dx,
                                                              self.dx))

        self._sampler = (self._continuous_sampler if smear
                         else self._discrete_sampler)

    def Sample(self):
        """
        Return a sample of each distribution.
        """
        return self._sampler()

    def ComputeThresholds(self, noise_cut=0.99, pes_to_adc=None):
        """
        Find the number of pes at which each noise distribution leaves behind
        the a given fraction of its population.

        Parameters
        ----------
        noise_cut : float
            Fraction of the distribution to be left behind. Default is 0.99.
        pes_to_adc : float or array of floats, optional
            Constant(s) for adc to pes conversion (default None).
            If not present, constants are taken from DataFrame given to the
            constructor.

        Returns
        -------
        cuts: array of floats
            Cuts in adc.
        """
        def findcut(probs):
            return self.xbins[probs > noise_cut][0]

        if pes_to_adc is None:
            pes_to_adc = self._df["adc_to_pes"]

        cumprobs = np.apply_along_axis(np.cumsum, 1, self.probs)
        return np.apply_along_axis(findcut, 1, cumprobs) * pes_to_adc
