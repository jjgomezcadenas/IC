"""
Select only those events with good baseline recovery.
"""
from __future__ import print_function

import sys

import numpy as np

import Core.wfmFunctions as wfm


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
    apply_on : string, optional
        Apply filter on different data types: CWF (defult), BLR
    """

    def __init__(self, **opts):
        if "max_adc" not in opts:
            print("max_adc option not found")
            sys.exit(2)

        self.max_adc = opts["max_adc"]
        self.n_samples = opts.get("n_samples", 10)
        self.wftype = opts.get("apply_on", "CWF")
        if self.wftype not in ["CWF", "BLR"]:
            raise ValueError("Wrong value for argument apply_on: {}\n"
                             "Available options are 'CWF' and 'BLR'")

    def __call__(self, f, i):
        if self.wftype == "CWF":
            data = f.root.RD.pmtcwf[i]
        else:
            data = wfm.subtract_baseline(f.root.RD.pmtblr[i])
        means = np.mean(data[:,-self.n_samples:], axis=1)
        return np.all(means < self.max_adc)
