import numpy             as np
import h5py

tot_evts = 20000   # total number of events to slice
blk_evts = 1000    # number of events per block

if(tot_evts % blk_evts != 0):
    print "ERROR: tot_evts must be evenly divisible by blk_evts"
    exit()

nfiles = tot_evts / blk_evts
cmb_type = "signal"

h5_out = h5py.File("out_data/{0}_2x2x2_combined.h5".format(cmb_type), 'w')

nslice_arr = []
nslice = 0
for nf in range(nfiles):

    estart = nf*blk_evts
    eend = estart + blk_evts

    # Open the file.
    fname = "out_data/{0}_2x2x2_{1}_to_{2}.h5".format(cmb_type,estart,eend-1)
    h5_in = h5py.File(fname,'r')

    # Concatenate the slice number array.
    tnslice_arr = h5_in["nspevt"]
    nslice_arr = np.concatenate((nslice_arr,tnslice_arr))
    
    # Write all slices and sipm maps.
    nslices_in_file = (len(h5_in)-1)/2
    print "Combining file {0} with total of {1} slices".format(fname,nslices_in_file)
    for ns in range(nslices_in_file):
        #h5_out.create_dataset("slice{0}".format(nslice),data=h5_in["slice{0}".format(ns)],dtype='float32')
        h5_out.create_dataset("sipm{0}".format(nslice), data=h5_in["sipm{0}".format(ns)], dtype='float32')
        nslice += 1

    h5_in.close()
h5_out.create_dataset("nspevt",data=nslice_arr) 
h5_out.close()
