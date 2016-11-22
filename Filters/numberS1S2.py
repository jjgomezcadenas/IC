"""
Select only those events with N S1 and M S2. UNKNOWN signals are ignored.
"""
from __future__ import print_function

import sys

import Core.Bridges as bdg
import Core.tblFunctions as tbl


class NumberS1S2:
    """
    Select only those events with N S1 and M S2. UNKNOWN signals are ignored.

    Parameters
    ----------
    nS1 : int, optional
        Number of S1 signals
    nS2 : int, optional
        Number of S2 signals
    """

    def __init__(self, **opts):
        if "nS1" not in opts and "nS2" not in opts:
            print("options nS1, nS2 not found."
                  " At least one of them must be given")
            sys.exit(2)

        self.nS1 = opts["nS1"]
        self.nS2 = opts["nS2"]

    def __call__(self, f, i):
        pmap = tbl.read_pmap(f.root.PMAPS.PMaps, i)
        if len(pmap.get(Signal.S1)) != self.nS1:
            return False
        if len(pmap.get(Signal.S2)) != self.nS2:
            return False
        return True
