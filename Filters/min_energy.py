"""
Select only those events with energy above some threshold.
"""
from __future__ import print_function

import sys
import numpy as np

import wfmFunctions as wfm


class min_energy:
    """
    Select only those events with integrated energy above some threshold.

    Parameters
    ----------
    min_signal : int or float
        Minimum amount of signal in the PMTs.

    apply_on : string, optional
        WF to be used for computing signal. Options are: "ZS" (default)
        "RWF" or "CWF".
    """

    def __init__(self, **opts):
        if "min_signal" not in opts:
            print("min_signal option not found")
            sys.exit(2)

        self.min_signal = opts["min_signal"]
        self.wftype = opts.get("apply_on", "ZS")

    def __call__(self, f, i):
        if self.wftype == "ZS":
            pmts = f.root.ZS.PMT
        elif self.wftype == "CWF":
            pmts = f.root.RD.pmtcwf
        elif self.wftype == "RWF":
            pmts = f.root.RD.pmtrwf

        ene = pmts[i].sum()
        return ene > self.min_signal
