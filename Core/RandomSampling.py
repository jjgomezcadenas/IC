"""
Defines a class for random sampling.
"""
from __future__ import print_function

import tables as tb
import numpy as np


class NoiseSampler:
    def __init__(self, filename, sample_size=1, smear=True):
        """
        Samples a histogram as if it was a PDF.

        Parameters
        ----------
        filename : string
            Path and name to the hdf5 file containing the noise distributions.
        sample_size: int
            Number of samples per sensor and call.
        smear: bool
            Flag to choose between performing discrete or continuous sampling.

        Attributes
        ---------
        baselines : array of floats
            Pedestal for each SiPM.
        xbins : numpy.ndarray
            Contains the the bins centers in pes.
        dx: float
            Half of the bin size.
        probs: numpy.ndarray
            Matrix holding the probability for each sensor at each bin.
        nsamples: int
            Number of samples per sensor taken at each call.
        """
        self.nsamples = sample_size

        # Read data, take xbins, compute (half of) bin size and normalize
        # probabilities.
        h5in = tb.open_file(filename)
        data = h5in.root.data

        self.baselines = np.copy(data.attrs.baselines)
        self.baselines = self.baselines.reshape(self.baselines.shape[0], 1)

        self.xbins = np.copy(data.attrs.bins)
        self.dx = np.diff(self.xbins)[0] * 0.5

        def norm(ps):
            return ps/np.sum(ps) if ps.any() else ps

        self.probs = np.apply_along_axis(norm, 1, data[:])
        h5in.close()

        # Sampling functions
        def _sample_sensor(probs):
            if not probs.any():
                return np.zeros(self.nsamples)
            return np.random.choice(self.xbins, size=self.nsamples, p=probs)


        _discrete_sampler = lambda: np.apply_along_axis(_sample_sensor,
                                                        1, self.probs)
        _continuous_sampler = lambda: (_discrete_sampler() +
                                       np.random.uniform(-self.dx, self.dx))

        self._sampler = _continuous_sampler if smear else _discrete_sampler

    def Sample(self):
        """
        Return a sample of each distribution.
        """
        return self._sampler() + self.baselines

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
            If not present, the thresholds are given in pes.

        Returns
        -------
        cuts: array of floats
            Cuts in adc or pes.
        """
        def findcut(probs):
            if not probs.any():
                return np.inf
            return self.xbins[probs > noise_cut][0]

        if pes_to_adc is None:
            pes_to_adc = np.ones(self.probs.shape[0])
        pes_to_adc.reshape(self.probs.shape[0], 1)

        cumprobs = np.apply_along_axis(np.cumsum, 1, self.probs)
        return np.apply_along_axis(findcut, 1, cumprobs) * pes_to_adc
