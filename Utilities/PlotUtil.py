"""
A utility module for plots with matplotlib
"""
from Util import *
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import numpy as np
import os
import sys

def hbins(x, nsigma=5, nbins=10):
  xmin =np.average(x) - nsigma*np.std(x)
  xmax =np.average(x) + nsigma*np.std(x)
  bins = np.linspace(xmin, xmax, nbins)
  return bins
  
def HSimple1(x,nbins,title='hsimple',xlabel = '', ylabel = 'Frequency', 
             save=False,filename='hsimple.png', filepath='./'):
  plt.hist(x, nbins, histtype='bar', alpha=0.75)
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  
  if save:
    pathfile = filepath+filename
    print "saving histogram %s in %s"%(filename, pathfile)
    plt.savefig(pathfile, bbox_inches='tight')
    plt.close()
  else:
    plt.figure()


def plot_signal(signal_t,signal, 
                title = 'signal', signal_start=0, signal_end=1e+4, units=''):

  ax1 = plt.subplot(1,1,1)
  ax1.set_xlim([signal_start, signal_end])
  SetPlotLabels(xlabel='t (ns)', ylabel='signal (%s)'%units)
  plt.title(title)
  plt.plot(signal_t, signal)
  plt.show()

def plot_signal2(signal_t,signal, 
                title = 'signal', signal_start=0, signal_end=1e+7,
                signal_min=-1e+7, signal_max=1e+7, 
                units=''):

  ax1 = plt.subplot(1,1,1)
  ax1.set_xlim([signal_start, signal_end])
  ax1.set_ylim([signal_min, signal_max])
  SetPlotLabels(xlabel='t (ns)', ylabel='signal (%s)'%units)
  plt.title(title)
  plt.plot(signal_t, signal)
  plt.figure()

def pulse_plot(pulse_time, pulse_volt):
  """
  Plots pulse
  """
  plt.plot(pulse_time, pulse_volt)
  plt.show()

def SetPlotLabels(xlabel="", ylabel="",grid=True):
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if grid == True:
    plt.grid(which='both', axis='both')


def PrintSignal(title,signal_t,signal):

  print "%s len: signal_t =%d, signal = %d "%(title,len(signal_t), len(signal))
  #print "signal_t =", signal_t
  #print "signal =", signal
  
def PlotPE(cbin,cnt,len_signal_max, s1_l,s1_r,s2_l,s2_r):
    ax1 = plt.subplot(3,1,1)
    ax1.set_xlim([0, len_signal_max/mus])
    SetPlotLabels(xlabel='t (mus)', ylabel='signal (PES)')
    
    plt.plot(cbin/mus, cnt)

    ax2 = plt.subplot(3,1,2)
    ax2.set_xlim([s1_l/mus, s2_r/mus])
    SetPlotLabels(xlabel='t (mus)', ylabel='signal (PES)')
    
    plt.plot(cbin/mus, cnt)
   
    ax3 = plt.subplot(3,1,3)
    ax3.set_xlim([s2_l, s2_r])
    SetPlotLabels(xlabel='t (mus)', ylabel='signal (PES)')
    
    plt.plot(cbin/ns, cnt)
    plt.show()

