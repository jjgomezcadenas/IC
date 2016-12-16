#
# functions to study the Na22 in IC
#
# author: JA Hernando
# date: 21/11/2016

import numpy as np
import Calib.calib as cb   # TODO move functions from here to its place

import qxyFunctions as qxyf

import matplotlib
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
matplotlib.style.use('ggplot')


LDIM = 250.


def polo_qxy_slice(qs, xs, ys, q0=0., xlim=(-LDIM, LDIM), ylim=(-LDIM, LDIM)):
    fig, ax = plt.subplots(figsize=(6, 5))
    cqs, cxs, cys = qxyf.qxy_slice(qs, xs, ys, q0=q0, xlim=xlim, ylim=ylim)
    img = ax.scatter(cxs, cys, c=cqs, s=50, alpha=0.5, cmap=plt.cm.jet)
    fig.colorbar(img, ax=ax)
    return fig


def polo_qxy_slices(qs, xs, ys, q0=0., xlim=(-LDIM, LDIM), ylim=(-LDIM, LDIM)):
    nsamples = qs.shape[1]
    nx, ny, figsize = cb.plt_subplots(nsamples)
    fig, axs = plt.subplots(nx, ny, figsize=figsize)
    axs = axs.ravel()
    for i in range(nsamples):
        qsi = qs[:, i]
        cqs, cxs, cys = qxyf.qxy_slice(qsi, xs, ys, q0=q0, xlim=xlim, ylim=ylim)
        # print('slice ', i, ' size ', len(cqs))
        if (len(cqs) <= 0):
            continue
        img = axs[i].scatter(cxs, cys, c=cqs, s=50, alpha=0.5, cmap=plt.cm.jet)
        axs[i].set_xlim(xlim)
        axs[i].set_ylim(ylim)
        fig.colorbar(img, ax=axs[i])
    return fig


def polo_qxy_slices3d(qs, xs, ys, q0=0., xlim=(-LDIM, LDIM),
                      ylim=(-LDIM, LDIM)):
    nsamples = qs.shape[1]
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    for i in range(nsamples):
        qsi = qs[:, i]
        cqs, cxs, cys = qxyf.qxy_slice(qsi, xs, ys, q0=q0, xlim=xlim, ylim=ylim)
        czs = (i+1)*np.ones(len(cqs))
        img = ax.scatter(cxs, cys, czs, c=cqs, s=50, alpha=0.5)
    fig.colorbar(img, ax=ax)
    return fig


def polo_qxy_track(track):
    xs = [hit[1] for hit in track]
    ys = [hit[2] for hit in track]
    qs = [hit[0] for hit in track]
    fig, ax = plt.subplots(figsize=(8, 6))
    ct = ax.scatter(xs, ys, c=qs, alpha=0.5, s=200, cmap=plt.cm.jet)
    # axs[0].set_xlim(xlim)
    # axs[0].set_ylim(ylim)
    ax.plot(xs, ys)
    fig.colorbar(ct, ax=ax)
    return fig


def polo_qxy_track3d(track):
    xs = [hit[1] for hit in track]
    ys = [hit[2] for hit in track]
    qs = [hit[0] for hit in track]
    fig, axs = plt.subplots()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    zs = range(len(xs))
    img = ax.scatter(xs, ys, zs, c=qs, alpha=0.5, s=200)
    ax.plot(xs, ys, zs)
    fig.colorbar(img, ax=ax)
    return fig
