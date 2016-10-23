import sys, os
from tables import *

def checkFiles(dirname):
    for filename in os.listdir(dirname):
        if '.h5' in filename:
            h5file = open_file(filename, "a")
            pmt = h5file.get_node("/", "pmtrd")
            print 'File: {0} -> Chunk size: {1}'.format(filename,pmt.chunkshape)
            h5file.close()


if __name__ == "__main__":
    times = checkFiles(sys.argv[1])
