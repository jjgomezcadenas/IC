"""
Example showing how to store time series of lists of variable length.

Author: Francesc Alted
Date: 2016-08-31
"""

from __future__ import print_function

from time import time
import random
import numpy as np
import tables

NEVENTS = 1000*10    # The number of events
FILENAME = "s12.h5"  # The name of the output file

# Use one of the compression param set below
#FILTERS = tables.Filters(complevel=0)   # no compression
#FILTERS = tables.Filters(complevel=1, complib="zlib")  # zlib
FILTERS = tables.Filters(complevel=5, complib="blosc")  # blosc


class S12(tables.IsDescription):
    event = tables.UInt32Col(pos=0)
    kind = tables.UInt32Col(pos=1)
    t = tables.Int32Col(pos=2)
    value = tables.Float32Col(pos=3)


def store(event, t0, s, kind, table):
    """Store S1 or S2 time series in `table`"""
    row = table.row
    for i, val in enumerate(s):
        row['event'] = event
        row['kind'] = kind
        row['t'] = t0 + i
        row['value'] = val
        row.append()
    table.flush()


def store_s1(event, t, table):
    """Build and store S1 time series in `table`"""
    # Build a S1
    np_s1 = random.randint(5, 50)  # S1 events have between 5 and 50 elements
    s1 = np.random.randn(np_s1) * 10 + 10
    # Store S1
    store(event, t, s1, 1, table)


def store_s2(event, t, table):
    """Build and store S2 time series in `table`"""
    # Build a S2
    np_s2 = random.randint(50, 500)  # S2 events have between 5 and 50 elements
    s2 = np.random.randn(np_s2) * 10 + 10
    # Store S2
    store(event, t, s2, 2, table)


def create_s12_table(nevents):
    """Create a table w/ events with different amounts of S1 and S2 objects on it."""
    f = tables.open_file(FILENAME, "w", filters=FILTERS)
    table = f.create_table(f.root, "t12", S12, "Store for S12")

    # Fill some S12 in table
    t = 0
    for i in range(NEVENTS):
        # Store between 1 and 3 events S1
        for j in range(random.randint(1, 3)):
            t += 10
            store_s1(i, t, table)
        # Store between 1 and 5 events S2
        for j in range(random.randint(1, 5)):
            t += 100
            store_s2(i, t, table)
    f.close()


def read_s12(nevents, s12):
    f = tables.open_file(FILENAME)
    table = f.root.t12
    t0 = time()
    sum_ = 0
    for row in table.where("(kind == s12) & (event < nevents)"):
        sum_ += row['value']
    print("Time for querying S%d objects: %.3f" % (s12, (time() - t0)))
    f.close()


def create_index():
    f = tables.open_file(FILENAME, "a")
    table = f.root.t12
    table.cols.event.create_index()
    f.close()


# Create the S12 table
print("Creating S12 table...")
create_s12_table(NEVENTS)

# Read the S1 objects in the first 3 events
read_s12(nevents=3, s12=1)

# Read the S2 objects in the first 2 events
read_s12(nevents=2, s12=2)

# Create an index for improved query times
create_index()
print("After indexing...")

# Read again and see time improvements
read_s12(nevents=3, s12=1)
read_s12(nevents=2, s12=2)
