import FEParam as FP
from Util import *
from scipy import signal as SGN
import numpy as np

# spe is parameterized as an isosceles trapezoid of diagonal 
#   xxxxx-xxxxx      
#           --- 
#         -     -
#       -         -
#     -             -
#   -                 -
#   1 2 3 4 515 4 3 2 1
# x_slope = 5 ns, x_flat = 1 ns

###########################################################
class SPE:

  def __init__(self,pmt_gain=4.5e6,x_slope = 5*ns,x_flat = 1*ns):
    """
    Defines a spe 
    """
    self.pmt_gain=pmt_gain
    self.x_slope=x_slope
    self.x_flat=x_flat
    
    self.spe_base = self.x_slope + self.x_flat
    self.spe_length = 2*self.x_slope + self.x_flat

    self.A=self.pmt_gain*eplus/self.spe_base  #current  
    self.V = self.A*FP.spe_i_to_v()
    self.nADC = int(self.A*FP.spe_i_to_adc())

    self.t = np.arange(0,self.spe_length,FP.time_step)
    nns = int(self.x_slope/FP.time_step)
    nnf = int(self.x_flat/FP.time_step)
    rise = np.linspace(0,self.A, num=nns)
    fall = np.linspace(self.A,0, num=nns)
    flat = self.A*np.ones(nnf)
    self.spe=np.concatenate((rise,flat,fall))
    
        
  def Spe(self):
    """
    Returns a SPE
    """
    return self.t,self.spe

  def SpePulse(self,t0, tmax=1e+6*ns):
    """
    Returns a SPE pulse at time t0
    with baseline extending in steps of time_step from 0 to tmax determined by DELTA_L
    """
    n = int(t0/FP.time_step)
    nmax = int(tmax/FP.time_step)
    if n >= nmax:
      print "error in SpePulse n = %d nmax = %d "%(
        n, nmax)
      sys.exit()

    
    DELTA=np.zeros(nmax)   #Dirac delta of size DELTA_L
    DELTA[n]=1
    step = FP.time_step/ns
    spe_pulse_t =np.arange(0,len(DELTA) + len(self.spe) -1,step)
    spe_pulse = SGN.convolve(DELTA, self.spe)
    
    return spe_pulse_t,spe_pulse

  def SpePulseTrain(self,tlow,tup,tmax=1e+6*ns):
    """
    Returns a train of SPE pulses between tlow and tup separated by tstep
    """
    nmin = int(tlow/FP.time_step)
    nmax = int(tup/FP.time_step)
    NMAX = int(tmax/FP.time_step)
    step = FP.time_step/ns

    if nmax >= NMAX:
      print "error in SpePulse train nmax = %d NMAX = %d "%(
        nmax, NMAX)
      sys.exit()

    DELTA=np.zeros(NMAX)
    DELTA[nmin:nmax+1] = 1
    spe_pulse_t =np.arange(0,len(DELTA) + len(self.spe) -1,step)
    spe_pulse = SGN.convolve(DELTA, self.spe)
      
    return spe_pulse_t,spe_pulse

  def SpePulseFromVectorPE(self,cnt):
    """
    Returns a train of SPE pulses corresponding to vector cnt
    """
    
    spe_pulse = SGN.convolve(cnt[0:-len(self.spe)+1], self.spe)
      
    return spe_pulse

  def __str__(self):
        
    s= """
        SPE:
        Gain = %7.2f I (muA) = %7.2f, V (mV) = %7.2f, ADC = %d
      """%(self.pmt_gain, self.A/microampere, self.V/mV, self.nADC )
    return s


if __name__ == '__main__':
  from PlotUtil import *

  signal_start=2000*ns
  signal_length=200*microsecond
  daq_window = 1*millisecond
  
  spe = SPE()
  print spe

  signal_end = signal_start + signal_length
  
  signal_t, signal_PE = spe.SpePulseTrain(signal_start,signal_end,daq_window)

  plot_signal(signal_t/ns,signal_PE/muA, 
                  title = 'Input Signal: PE Train', 
                  signal_start=0*ns, signal_end=len(signal_t)/ns, 
                  units='muA')

  plt.show()

