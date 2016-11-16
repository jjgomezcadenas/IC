"""
Select only those events with at least one SiPM with some minimum charge.
"""
from __future__ import print_function

import sys
import numpy as np

import Core.wfmFunctions as wfm


class min_charge:
    """
    Select only those events with at least one SiPM with some minimum charge.

    Parameters
    ----------
    min_signal : int or float
        Minimum amount of signal in the SiPMs.

    sipm_selection : list of integers, optional
        List of sipm indices in which to apply selection. Default is empty.

    mau_len : int, optional
        Window size to compute the baseline. Default is 200.

    integrate_wf : bool, optional
        Flag to apply selection on the integrated wf (True)
        or per sample (False)
    """

    def __init__(self, **opts):
        if "min_signal" not in opts:
            print("min_signal option not found")
            sys.exit(2)

        self.min_signal = opts["min_signal"]
        self.mau_len = opts.get("mau_len", 200)
        self.integrate_wf = opts.get("integrate_wf", False)

        self.do_selection = "sipm_selection" in opts
        if self.do_selection:
            self.selection = np.array(opts["sipm_selection"], ndmin=1)

    def __call__(self, f, i):
        sipms = f.root.RD.sipmrwf

        wfs = sipms[i][self.selection] if self.do_selection else sipms[i]
        wfs = wfm.subtract_baseline(wfs, self.mau_len)
        if self.integrate_wf:
            wfs = np.sum(wfs, axis=1)
        return np.max(wfs) > self.min_signal
