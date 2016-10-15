import os

tot_evts = 20000   # total number of events to slice
blk_evts = 1000    # number of events per block

if(tot_evts % blk_evts != 0):
    print "ERROR: tot_evts must be evenly divisible by blk_evts"
    exit()

nfiles = tot_evts / blk_evts
for nf in range(nfiles):

    estart = nf*blk_evts
    eend = estart + blk_evts
    call_str = "python SiPM_gen.py {0} {1} &> out_SiPM_gen_{2}_{3}.dat &".format(estart,eend,estart,eend-1)
    print "Calling {0}".format(call_str)
    os.system(call_str)
