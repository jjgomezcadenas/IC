"""
Example of how concatenate a number of files into another one.
Here EArrays are used, but the concept is similar with Tables.
"""

import numpy as np
import tables

NFILES = 10

raw_data = np.arange(1e6, dtype=np.int16)

# Create files (RAW data)
for i in range(NFILES):
    with tables.open_file("raw%d.h5"%i, "w") as f:
        f.create_earray(f.root, 'arr', obj=raw_data)

# Create destination file
f_dest = tables.open_file("dst.h5", "w")
# And dataset for concatenation
earray_dest = f_dest.create_earray(f_dest.root, "dest", tables.Int16Atom(), shape=(0,))
# Concatenate files
for i in range(NFILES):
    with tables.open_file("raw%d.h5" % i) as f:
        data_node = f.get_node('/arr')
        earray_dest.append(data_node[:])
f_dest.close()
