import FEParam as FP
from Util import *
from scipy import signal
import SPE as SP

 

###########################################################

def FilterCutoffFrequency(R=4700, C=8e-9):
	"""
	Takes R (in Ohms) and C (in Farad) and returns frequency in Htz
	"""
	
	return 1./(2*pi*R*C)

class Filter:

 	def __init__(self,type='high',fc=5E3,fs= 1e+9):
		"""
		Defines a Butterworth HPF (high pass frequency) or LPF (low pass frequencey) filter 
		the default sampling frequencey is 1 GHz (inverse of 1 ns time)
		type may be equal to hig or low
		"""
		self.fc = fc
		self.fs = fs
		self.W = 2*self.fc/self.fs
		self.type = type
		self.b, self.a = signal.butter(1, self.W, btype=self.type, analog=False, output='ba')
		self.ba, self.aa = signal.butter(1, 2*self.fc, btype=self.type, analog=True, output='ba')

	def FilterCoef(self):
		"""
		Returns the filter coefficients
		"""
  		return self.b,self.a

  	def FilterAnalogCoef(self):
		"""
		Returns the filter coefficients
		"""
  		return self.ba,self.aa

	def FilterPulse(self,pulse):
		"""
		Filters a pulse
		"""
  		return signal.lfilter(self.b,self.a, pulse)

  	def FilterResponse(self):
  		"""
  		Gives the response of the filter y frequency-amplitude
  		"""

  		self.w, self.h = signal.freqs(self.ba, self.aa)
  		return self.w, self.h


	def __str__(self):
        
		s= """
		Filter:
		fc = %7.2f Hz, fs = %7.2f, W = %7.2f Hz type = %s  
		"""%(self.fc, self.fs, self.W, self.type)
		return s

	

def FilterScan(spe,C):
  	spe_pulse_t, spe_pulse = spe.SpePulse(100*ns)
	print "Capacitance of filter (nF) = %7.2f"%(C/nF)

	Ra = np.logspace(0., 3.0, num=3)
    
  	for R in Ra:
  		print "Resistance of filter (O) = %7.2f"%(R)

  		freq_HPF = FilterCutoffFrequency(R=R*ohm, C=C)
  		print "cutoff frequence (kHz) = %7.2f"%(freq_HPF/1e+3)
  		hpf = Filter(type='high',fc=freq_HPF,fs=f_sample)
  		w,h = hpf.FilterResponse()
  		plt.subplot(3,1,1)
		plt.xscale('log')
		plt.title('HPF  frequency response')
		SetPlotLabels(xlabel='Frequency [radians / second]', ylabel='Amplitude [dB]')

		plt.margins(0, 0.1)
		plt.axvline(freq_HPF, color='green') # cutoff frequency
		plt.plot(w, 20 * np.log10(abs(h)))

		spe_hp_pulse = hpf.FilterPulse(spe_pulse)

		ax1 = plt.subplot(3,1,2)
		ax1.set_xlim([0, 200])
		SetPlotLabels(xlabel='t (ns)', ylabel='I (muA)')
  		plt.title('SPE pulse')
		plt.plot(spe_pulse_t/ns, spe_pulse/microampere)

		ax2 = plt.subplot(3,1,3)
		ax2.set_xlim([0, 200])
		SetPlotLabels(xlabel='t (ns)', ylabel='I (muA)')
  		plt.title('SPE pulse after HPF')
		plt.plot(spe_pulse_t, spe_hp_pulse/microampere)

		plt.show()

def FilterCharacterization(spe):
	spe_pulse_t, spe_pulse = spe.SpePulse(100*ns)

  	hpf = Filter(type='high',fc=FP.freq_HPF,fs=FP.f_sample)
  	lpf = Filter(type='low',fc=FP.freq_LPF,fs=FP.f_sample)

  	w,h = hpf.FilterResponse()

  	plt.subplot(2,1,1)
	plt.xscale('log')
	plt.title('HPF  frequency response')
	SetPlotLabels(xlabel='Frequency [radians / second]', ylabel='Amplitude [dB]')

	plt.margins(0, 0.1)
	plt.axvline(FP.freq_HPF, color='green') # cutoff frequency
	plt.plot(w, 20 * np.log10(abs(h)))

	
	w,h = lpf.FilterResponse()
	plt.subplot(2,1,2)

  	
	plt.xscale('log')
	plt.title('LPF  frequency response')
	SetPlotLabels(xlabel='Frequency [radians / second]', ylabel='Amplitude [dB]')

	plt.margins(0, 0.1)
	plt.axvline(FP.freq_LPF, color='green') # cutoff frequency
	plt.plot(w, 20 * np.log10(abs(h)))

	plt.show()

  	spe_hp_pulse = hpf.FilterPulse(spe_pulse)
  	spe_hp_lp_pulse = lpf.FilterPulse(spe_hp_pulse)

  	ax1 = plt.subplot(3,1,1)
	ax1.set_xlim([0, 200])
	SetPlotLabels(xlabel='t (ns)', ylabel='I (muA)')
  	plt.title('SPE pulse')
	plt.plot(spe_pulse_t/ns, spe_pulse/microampere)

	ax2 = plt.subplot(3,1,2)
	ax2.set_xlim([0, 200])
	SetPlotLabels(xlabel='t (ns)', ylabel='I (muA)')
  	plt.title('SPE pulse after HPF')
	plt.plot(spe_pulse_t, spe_hp_pulse/microampere)

	ax3 = plt.subplot(3,1,3)
	ax3.set_xlim([0, 200])
	SetPlotLabels(xlabel='t (ns)', ylabel='I (muA)')
  	plt.title('SPE pulse after HPF + LPF')
	plt.plot(spe_pulse_t, spe_hp_lp_pulse/microampere)

	plt.show()

def FilterTrain(spe,ti=100*ns,tf=200*ns, tm=50*ns):
	spe_pulse_t, spe_pulse = spe.SpePulseTrain(ti,tf)

	print "pulse train: ti = %7.2f ns, tf = %7.2f ns, area = %7.2g muA"%(ti,tf,spe_pulse.sum())
	
  	hpf = Filter(type='high',fc=FP.freq_HPF,fs=FP.f_sample)
  	lpf = Filter(type='low',fc=FP.freq_LPF,fs=FP.f_sample)

  	w,h = hpf.FilterResponse()

  	plt.subplot(2,1,1)
	plt.xscale('log')
	plt.title('HPF  frequency response')
	SetPlotLabels(xlabel='Frequency [radians / second]', ylabel='Amplitude [dB]')

	plt.margins(0, 0.1)
	plt.axvline(FP.freq_HPF, color='green') # cutoff frequency
	plt.plot(w, 20 * np.log10(abs(h)))

	
	w,h = lpf.FilterResponse()
	plt.subplot(2,1,2)

  	
	plt.xscale('log')
	plt.title('LPF  frequency response')
	SetPlotLabels(xlabel='Frequency [radians / second]', ylabel='Amplitude [dB]')

	plt.margins(0, 0.1)
	plt.axvline(FP.freq_LPF, color='green') # cutoff frequency
	plt.plot(w, 20 * np.log10(abs(h)))

	plt.show()

  	spe_hp_pulse = hpf.FilterPulse(spe_pulse)
  	spe_hp_lp_pulse = lpf.FilterPulse(spe_hp_pulse)
  	area_hp = np.sum(spe_hp_pulse[np.where(spe_hp_pulse>0)])
  	area_hp_lp = np.sum(spe_hp_lp_pulse[np.where(spe_hp_lp_pulse>0)])


  	print "pulse area after HPF = %7.2g muA, pulse area after HPF + LPF  = %7.2g muA"%(
  		area_hp,area_hp_lp)

  	ax1 = plt.subplot(3,1,1)
	ax1.set_xlim([0, tf+tm])
	SetPlotLabels(xlabel='t (ns)', ylabel='I (muA)')
  	plt.title('SPE pulse')
	plt.plot(spe_pulse_t/ns, spe_pulse/microampere)

	ax2 = plt.subplot(3,1,2)
	ax2.set_xlim([0, tf+tm])
	SetPlotLabels(xlabel='t (ns)', ylabel='I (muA)')
  	plt.title('SPE pulse after HPF')
	plt.plot(spe_pulse_t, spe_hp_pulse/microampere)

	ax3 = plt.subplot(3,1,3)
	ax3.set_xlim([0, tf+tm])
	SetPlotLabels(xlabel='t (ns)', ylabel='I (muA)')
  	plt.title('SPE pulse after HPF + LPF')
	plt.plot(spe_pulse_t, spe_hp_lp_pulse/microampere)

	plt.show()



if __name__ == '__main__':
	from PlotUtil import *

	spe = SP.SPE(pmt_gain=FP.PMT_GAIN,x_slope = 5*ns,x_flat = 1*ns)
	#FilterScan(spe,10E-9)
	FilterCharacterization(spe)
	#FilterTrain(spe,ti=100*ns,tf=200*ns, tm=50*ns)
  
  	

