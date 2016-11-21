"""
Select only those events with good baseline recovery.
"""
from __future__ import print_function

import sys


class Baseline:
    """
    Select only those events with a good baseline recovery.

    Parameters
    ----------
    max_adc : int
        Minimum adc at the end of the waveform to consider the waveform valid.
    n_samples : int, optional
        Number of samples at the end of the waveform to measure waveform
        lifting. Default is 10.
    """

    def __init__(self, **opts):
        if "max_adc" not in opts:
            print("max_adc option not found")
            sys.exit(2)

        self.max_adc = opts["max_adc"]
        self.n_samples = opts.get("n_samples", 10)

    def __call__(self, f, i):
        means = np.mean(f.root.RD.pmtcwf[i,:,-n_samples:], axis=1)
        return np.all(means < self.max_adc)
