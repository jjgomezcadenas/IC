#
# Functions to deal with points in the anode
#

import numpy as np

# size of the chamber in x or y
LDIM = 250.


def qxy_slice(qs, xs, ys, q0=0., xlim=(-LDIM, LDIM), ylim=(-LDIM, LDIM)):
    def xfilter(zi):
        if (zi[0] < q0):
            return False
        if ((zi[1] < xlim[0]) or (zi[1] > xlim[1])):
            return False
        if ((zi[2] < ylim[0]) or (zi[2] > ylim[1])):
            return False
        return True
    zz = zip(qs, xs, ys)
    zf = filter(xfilter, zz)
    cqs = np.array([zi[0] for zi in zf])
    cxs = np.array([zi[1] for zi in zf])
    cys = np.array([zi[2] for zi in zf])
    return cqs, cxs, cys


def qxy_point(qs, xs, ys, q0=0., xlim=(-LDIM, LDIM), ylim=(-LDIM, LDIM)):
    cqs, cxs, cys = qxy_slice(qs, xs, ys, q0=q0, xlim=xlim, ylim=ylim)
    if (len(cqs) <= 0):
        return (-LDIM, -LDIM, -LDIM), (LDIM*LDIM, LDIM*LDIM, LDIM*LDIM)
    tqi = np.sum(cqs)
    txi = np.dot(cxs, cqs)/tqi
    tyi = np.dot(cys, cqs)/tqi
    dx = cxs-txi
    sxi = np.dot(dx*dx, cqs)/tqi
    dy = cys-tyi
    syi = np.dot(dy*dy, cqs)/tqi
    sxyi = np.dot(dx*dy, cqs)/tqi
    return (tqi, txi, tyi), (sxi, syi, sxyi)


def qxy_track(qs, xs, ys, q0=0., xlim=(-LDIM, LDIM), ylim=(-LDIM, LDIM)):
    n = qs.shape[1]
    track = []
    for i in range(n):
        qsi = qs[:, i]
        point, cov = qxy_point(qsi, xs, ys, q0=q0, xlim=xlim, ylim=ylim)
        track.append(point)
    ctrack = [p for p in track if (p[0] > 0.)]
    return ctrack
