from Util import *
from PlotUtil import *
import FEParam as FP
import SPE as SP
from scipy import signal as SGN
import numpy as np
from deModel import arrayFixedInt
from deModel import DeFixedInt


def DownScaleSignal(signal_t, signal, scale):
 	"""
 	downscales the time and signal vectors
 	"""
	signal_d=SGN.decimate(signal,scale,ftype='fir')
	signal_t_d = SGN.decimate(signal_t,scale,ftype='fir')

	#pulse_plot(signal_t_d, signal_d)
	return signal_t_d, signal_d


##########################################################
class Filter:

 	def __init__(self,ftype='high',fc=5E3,fs= 1e+9):
		"""
		Defines a Butterworth HPF (high pass frequency) or LPF (low pass frequencey) filter 
		the default sampling frequencey is 1 GHz (inverse of 1 ns time)
		type may be equal to hig or low
		"""
		self.fc = fc
		self.fs = fs
		self.W = 2*self.fc/self.fs
		self.ftype = ftype
		self.b, self.a = SGN.butter(1, self.W, btype=self.ftype, analog=False, output='ba')
		self.ba, self.aa = SGN.butter(1, 2*self.fc, btype=self.ftype, analog=True, output='ba')

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
  		return SGN.lfilter(self.b,self.a, pulse)

  	def FilterResponse(self):
  		"""
  		Gives the response of the filter y frequency-amplitude
  		"""

  		self.w, self.h = SGN.freqs(self.ba, self.aa)
  		return self.w, self.h


	def __str__(self):
        
		s= """
		Filter:
		fc = %7.2f Hz, fs = %7.2f, W = %7.2f Hz type = %s  
		"""%(self.fc, self.fs, self.W, self.ftype)
		return s

	

###########################################################

class FEE:
	"""
	Emulates the effect of the PMT + FEE:
	1) PMT gain G
	2) A HPF, due to the FEE decoupling capacitor C and the associated series resitor R
	3) A LPF that shapes the signal, with a frequencey f and a noise
	4) A resitor gain RG to transform current in voltage
	"""

 	def __init__(self,PMTG=4.5e6, C=6.75*nF,R= 2350*ohm, 
 				 f=2E6*hertz, fn=2E5*hertz, RG=250*ohm):

 		self.PMTG = PMTG
 		self.C = C
 		self.R = R
 		self.f_LPF = f
 		self.fn_LPF = fn
 		self.RG = RG
 		self.f_HPF=(1./(2*pi*R*C))
 		self.hpf = Filter(ftype='high',fc=self.f_HPF,fs=FP.f_sample)
		self.lpf = Filter(ftype='low',fc=self.f_LPF,fs=FP.f_sample)
		self.hpfr = Filter(ftype='high',fc=self.f_HPF,fs=FP.f_sample_DAQ)


	def Filter(self,signal):
		"""
		for an input signal in pes, returns the effect of the electronics
		"""
		signal_hp = self.hpf.FilterPulse(signal)
  		signal_hp_lp = self.lpf.FilterPulse(signal_hp)
  	
  		return signal_hp, signal_hp_lp

  	def FilterInverse(self, signal):
  		"""
  		Returns the inverse of the HPF filter in bins of 1ns
  		"""
  		b_HPr,a_HPr = self.hpf.FilterCoef()
		h_t_inv = SGN.lfilter(a_HPr,b_HPr, signal)
		return h_t_inv

	def FilterInverseDAQ(self, signal):
		"""
		Returns the inverse of the HPF filter in bins of time_DAQ (25 ns)
		"""
  		b_HPr,a_HPr = self.hpfr.FilterCoef()
		h_t_inv = SGN.lfilter(a_HPr,b_HPr, signal)
		return h_t_inv

	def FEESignal(self,signal_current, noise_rms=0.3*mV):
 		"""
 		filters the input signal according to the filters and transforms it in volts
 		"""
 		signal_hp, signal_hp_lp = self.Filter(signal_current)
 		noise = self.FEENoise(len(signal_current), noise_rms)
 		return signal_hp_lp*self.RG + noise


 	def DAQSignal(self,signal_t, signal_fee, noise_rms=0.3*mV):
 		"""
 		downscale the signal after the FEE
 		
 		"""

		signal_t_d, signal_d = DownScaleSignal(signal_t, signal_fee, int(FP.time_DAQ))
		signal_daq = signal_d/FP.voltsToAdc
		noise = self.FEENoise(len(signal_daq), noise_rms/FP.voltsToAdc)
 		return signal_t_d, signal_daq + noise

 	def FEENoise(self, signal_length, noise_rms=0.3*mV):
 		"""
 		filters the input signal according to the filters and transforms it in volts
 		"""
 		if noise_rms > 0:
 			noise = np.random.normal(0, noise_rms, signal_length)
 			return noise
 		else:
 			return np.zeros(signal_length)


 	def VSignal(self,signal_current):
 		"""
 		Takes a current input signal and transforms it to volts
 		"""
 		return signal_current*self.RG

 	def InverseSignal(self, signal_t):
 		"""
 		Computes the inverse signal for deconvolution 
 		"""
		pulse = np.zeros(len(signal_t))
		pulse[0]=1
	
		signal_fee_inv = self.FilterInverse(pulse)
		return signal_fee_inv

	def InverseSignalDAQ(self, signal_t):
		"""
		Computes the inverse signal in bins of DAQ
		"""
		pulse = np.zeros(len(signal_t))
		pulse[0]=1
	
		signal_daq_inv = self.FilterInverseDAQ(pulse)
		return signal_daq_inv



 	def __str__(self):
        
		s= """
		NEW FEE
  		PMT gain = %7.2g
  		decoupling capacitor = %7.2f nF
  		decoupling resistor = %7.2f ohm
  		resitor gain = %7.2f ohm
  		HPF frequency = %7.2g Hz  W_HPF = %7.2g 
  		LPF frequency = %7.2g Hz  W_LPF = %7.2g
  		LPF frequency noise = %7.2g Hz  
  """%(self.PMTG,self.C/nF,self.R/ohm, self.RG/ohm,
  		self.f_HPF/hertz,self.hpf.W,
  		self.f_LPF/hertz, self.lpf.W,
  		self.fn_LPF/hertz)
		
		return s

###########################################################

class DAQ:
	"""
	Emulates the effect of the DAQ:
	1) downsamples the signals to the sampling frequence of DAQ (400 MHz)
	2) Quantizes the signal according to the number of bits in the DAQ
	"""

 	def __init__(self,NBITS=12, NBITS_FRAC=11,time_sample=25*ns, 
 				 LSB = 2*volt):
 		
 		self.NBITS=NBITS
 		self.NBITS_FRAC=NBITS_FRAC
 		self.time_sample = time_sample
 		self.f_sample=1./self.time_sample
 		self.LSB = LSB
 		lsb = self.LSB/2**self.NBITS
 		self.voltsToAdc = lsb/1.25     # 1.25 magic from experts
 		

 	def DAQSignal(self,signal_t, signal_fee, noise_rms=0.3*mV):
 		"""
 		downscale the signal after the FEE, quantizes and 
 		transforms to ADC counts. Returns back in floating point
 		
 		"""
 		
 		#pulse_plot(signal_t, signal_fee)
 		sf = int(FP.f_sample/self.f_sample)
		signal_t_d, signal_d = DownScaleSignal(signal_t, signal_fee, sf)
		signal_daq1 = signal_d/self.voltsToAdc

		#signal_daq_fi  = arrayFixedInt(self.NBITS,self.NBITS_FRAC, signal_daq1)
		#signal_daq = fpi_array_fValue(signal_daq_fi)
		
		signal_daq = signal_daq1
		
		if noise_rms > 0:
 			noise = np.random.normal(0, noise_rms/self.voltsToAdc, len(signal_daq))
 			return signal_t_d, signal_daq + noise
 		else:
 			return signal_t_d, signal_daq
		

	def __str__(self):

		s= """
		NEW DAQ
		number of bits = %d, number of bits fraction = %d
		sampling time = %7.2f ns frequencey = %7.2g hertz		  
		"""%(self.NBITS,self.NBITS_FRAC,
			self.time_sample/ns, 
			self.f_sample/hertz)
		return s



def DeconvSimple(signal,signal_inv):
	"""
	Deconvolution of the fine-grained fee signal (no DAQ)
	no noise
	using true start and end of signals
	"""

	coef = signal_inv[100]
	
	print "coef = %7.4g"%coef

	acum = np.zeros(len(signal))

	acum[0]=coef*signal[0]
	for n in xrange(1,len(signal)):
		acum[n] = acum[n-1] + signal[n]

	signal_r = signal + coef*acum

	return signal_r



def DeconvDAQ(signal_PE_daq, signal_t_daq,signal_daq,signal_daq_inv):
	"""
	Deconvolution of the DAQ signal
	no noise
	using true start and end of signals
	"""

	# coef=fi(h_t_inv(100),1,NBITS_acum,NBITS_FRAC_acum)
	# coef = signal_daq_inv[100]
	
	coef = DeFixedInt(FP.NBITS_acum, FP.NBITS_FRAC_acum,signal_daq_inv[100])
	#coef = cf.value

	print "DeconvDAQ: coef = %s "%fpi_rep(coef)
	
	acum  = arrayFixedInt(FP.NBITS_acum, FP.NBITS_FRAC_acum, np.zeros(len(signal_daq)))
	acum[0]=coef*signal_daq[0]

	for n in xrange(1,len(signal_daq)):
		acum[n] = acum[n-1] + signal_daq[n]

	signal_r  = arrayFixedInt(FP.NBITS_acum, FP.NBITS_FRAC_acum, np.zeros(len(signal_daq)))

	for n in xrange(0,len(signal_r)):
		signal_r[n] = signal_daq[n] + coef*acum[n]

	plot_signal(signal_t_daq/ns,signal_PE_daq,
                title = 'PE after DAQ', 
                signal_start=0*ns, signal_end=len(signal_t_daq)*FP.time_DAQ, 
                units='adc counts')

	plot_signal(signal_t_daq/ns,fpi_array_fValue(signal_r),
                title = 'Deconv DAQ', 
                signal_start=0*ns, signal_end=len(signal_t_daq)*FP.time_DAQ, 
                units='adc counts')

	# print "DAQ signal"
	# print len(signal_r)
	# print fpi_array_fValue(signal_r)

	plt.show()
	
	wait()
	thr = 20.
	Q_i = FP.pulse_area_positive(signal_PE_daq)
	Q_r = fpi_sum(signal_r,thr)

	print """
		Deconv DAQ
		Q (input pulse) = %7.4g
		Q (recovered pulse) = %7.4g
		rms = %7.4g
		"""%(Q_i,Q_r, abs(Q_i -Q_r)/Q_i)

def DeconvDAQMAU(signal_PE_daq, signal_t_daq,signal_daq,signal_daq_inv):
	"""
	Deconvolution of the DAQ signal
	using MAUS
	"""

	return 0


def fpi_rep(fpi):
	s= """
	FPI for number %10.4f, 
	rep = %s, int value = %d
	nbits = %d nfrac = %d width = %d
	"""%(fpi.fValue,fpi.rep,fpi.value,
		fpi.intWidth,fpi.fractWidth,fpi.width)
	return s


def fpi_array_value(fpi_array):
	array_value = np.zeros(len(fpi_array))
	for i in xrange(0,len(fpi_array)):
			array_value[i]=fpi_array[i].value
	return array_value

def fpi_array_fValue(fpi_array):
	array_fValue = np.zeros(len(fpi_array))
	for i in xrange(0,len(fpi_array)):
			array_fValue[i]=fpi_array[i].fValue
	return array_fValue


def fpi_sum(pulse, thr):
	print "adding terms above thr =%7.2f"%thr
	fsum = pulse[0]
	
	for i in xrange(1, len(pulse)):	
		if pulse[i].fValue > thr:
			fsum+=pulse[i]


	#fsum =np.sum(pulse[np.where(pulse>thr)])
	#fsum = np.sum(fpi_array)  #sum in fixed point precision
	return fsum.fValue # returns the floating point value


def testFPI():

	a = DeFixedInt(FP.NBITS, FP.NBITS_FRAC, 175.5000)
	b = DeFixedInt(FP.NBITS, FP.NBITS_FRAC, 551.2831)
	c = DeFixedInt(FP.NBITS_acum, FP.NBITS_FRAC_acum, 0.0005)
	
	d = a + b*c
	x = 175.5000 + 551.2831*0.0005

	print """
	 a = %s
	 b = %s
	 c = %s
	 d = %s

	"""%(fpi_rep(a),fpi_rep(b),fpi_rep(c),fpi_rep(d))

	print " diff = %7.2g"%((x-d.fValue)/x)

	wait()

def DoDeconvSimple(signal_PE_v, signal_fee,signal_fee_inv, plot_signal_r):

	signal_r = DeconvSimple(signal_fee,signal_fee_inv)

	if plot_signal_r == 1 or plot_signal_r == 2 or plot_signal_r == 3: 
		plot_signal(signal_t/ns,signal_r/mV,
                	title = 'Deconv simple', 
                	signal_start=0, signal_end=len(signal_t), 
                	units='mV')

	if plot_signal_r == 2 or plot_signal_r == 3:
		plot_signal2(signal_t/ns,signal_r/mV,
                	title = 'Deconv simple', 
                	signal_start=5400, signal_end=5700,
                	signal_min=0, signal_max=10, 
                	units='mV')
	if plot_signal_r == 3:
		plot_signal2(signal_t/ns,signal_r/mV,
                	title = 'Deconv simple', 
                	signal_start=5700, signal_end=7000,
                	signal_min=0, signal_max=2, 
                	units='mV')

	Q_i = FP.pulse_area_positive(signal_PE_v/mV)
	Q_r = FP.pulse_area_threshold(signal_r/mV, 3.0*noise_rms/mV)

	print """
		Deconv Simple FEE
		Q (input pulse) = %7.5g
		Q (recovered pulse) = %7.5g (3*noise_rms cut)
		rms = %7.5g
		"""%(Q_i,Q_r, abs(Q_i -Q_r)/Q_i)
	plt.show()

def DoDeconvSimpleDAQ(signal_PE_adc, signal_daq,signal_inv_daq, plot_signal_r_DAQ):

	signal_r_DAQ = DeconvSimple(signal_daq,signal_inv_daq)

	if plot_signal_r_DAQ == 1 or plot_signal_r_DAQ == 2 or plot_signal_r_DAQ == 3: 
		plot_signal(signal_t_daq/ns,signal_r_DAQ,
                	title = 'Deconv simple DAQ', 
                	signal_start=0, signal_end=len(signal_t), 
                	units='adc counts')

	if plot_signal_r_DAQ == 2 or plot_signal_r_DAQ == 3:
		plot_signal2(signal_t_daq/ns,signal_r_DAQ,
                	title = 'Deconv simple', 
                	signal_start=5400, signal_end=5700,
                	signal_min=0, signal_max=20, 
                	units='adc counts')

	if plot_signal_r_DAQ == 3:
		plot_signal2(signal_t_daq/ns,signal_r_DAQ,
                	title = 'Deconv simple', 
                	signal_start=5700, signal_end=7000,
                	signal_min=0, signal_max=1, 
                	units='adc counts')

	time_bin = FP.time_DAQ/FP.time_step

	Q_i = FP.pulse_area_positive(signal_PE_adc)
	Q_r = time_bin*FP.pulse_area_threshold(signal_r_DAQ, 5*noise_DAQ)

	print """
		Deconv Simple DAQ
		Q (input pulse) = %7.5g
		Q (recovered pulse) = %7.5g (5*noise_rms cut)
		rms = %7.5g
		"""%(Q_i,Q_r, abs(Q_i -Q_r)/Q_i)


	plt.show()
if __name__ == '__main__':

	#testFPI()

	single_pe =False
	signal_start=2000*ns
	signal_length=5000*ns
	daq_window = 20*microsecond
	noise_rms=0.3*mV
	noise_DAQ = noise_rms/FP.voltsToAdc

	plot_signal_pe_muA = True
	plot_signal_pe_mV = True
	plot_signal_pe_adc = True
	plot_signal_noise = True
	plot_signal_fee = True
	plot_signal_daq = True
	plot_signal_r = 1
	plot_signal_fee_inv = True
	plot_signal_inv_daq = True
	plot_signal_r_DAQ = 1

	deconv_simple = True
	deconv_simple_DAQ = True
	

	spe = SP.SPE()
	print spe

	fee = FEE()
	print fee

	daq = DAQ()
	print daq

	FP.print_FEE()

	print "noise rms = %7.2f mV (%7.2f adc counts)"%(noise_rms/mV, noise_DAQ)

	#wait()

	
	# input signal in PEs (current)

	signal_end = signal_start + signal_length
  
  	signal_t = 0
  	signal_PE = 0
  	
  	if single_pe == True:
  		signal_t, signal_PE =spe.SpePulse(100*ns, tmax=1500*ns)
		noise_rms=0.

	signal_t, signal_PE = spe.SpePulseTrain(signal_start,signal_end,daq_window)
	

	if plot_signal_pe_muA:
		plot_signal(signal_t/ns,signal_PE/muA, 
                	title = 'Input Signal: PE Train', 
                	signal_start=0*ns, signal_end=len(signal_t)/ns, 
                	units='muA')

	#input signal in mV
	signal_PE_v = fee.VSignal(signal_PE)
	signal_PE_adc = fee.VSignal(signal_PE)/FP.voltsToAdc

	if plot_signal_pe_mV:
		plot_signal(signal_t/ns,signal_PE_v/mV, 
        	        title = 'Input Signal: PE Train', 
            	    signal_start=0*ns, signal_end=len(signal_t)/ns, 
                	units='mV')

	if plot_signal_pe_adc:
		plot_signal(signal_t/ns,signal_PE_adc, 
        	        title = 'Input Signal: PE Train', 
            	    signal_start=0, signal_end=len(signal_t), 
                	units='adc counts')

	#noise
	noise =fee.FEENoise(len(signal_t), noise_rms=noise_rms)

	if plot_signal_noise:
		plot_signal(signal_t/ns,noise/mV, 
                	title = 'FEE Noise', 
               	 	signal_start=0*ns, signal_end=len(signal_t)/ns, 
                	units='mV')

	# effect of electronics: pass filters and output voltage
	signal_fee = fee.FEESignal(signal_PE, noise_rms=noise_rms)

	
	if plot_signal_fee:
		plot_signal(signal_t/ns,signal_fee/mV, 
                	title = 'Signal after FEE (filtered)', 
                	signal_start=0*ns, signal_end=len(signal_t)/ns, 
                	units='mV')
	#inverse
	signal_fee_inv = fee.InverseSignal(signal_t)


	if plot_signal_fee_inv:
		plot_signal(signal_t/ns,signal_fee_inv,
                title = 'Inverse FEE', 
                signal_start=0*ns, signal_end=10*ns, 
                units='')

	#Signal out of DAQ

	signal_t_daq, signal_daq = fee.DAQSignal(signal_t, signal_fee, noise_rms=0)

	if plot_signal_daq:
		plot_signal(signal_t_daq/ns,signal_daq, 
                	title = 'Signal after DAQ ', 
                	signal_start=0, signal_end=len(signal_t), 
                	units='adc counts')

	#inverse
	signal_inv_daq = fee.InverseSignalDAQ(signal_t)  #in bins of DAQ (25 ns)

	if plot_signal_inv_daq:
		plot_signal(signal_t/ns,signal_inv_daq,
                title = 'Inverse DAQ', 
                signal_start=0*ns, signal_end=10*ns, 
                units='')

	#------deconv simple
	if deconv_simple:
	
		print "++++++Deconv simple"
		DoDeconvSimple(signal_PE_v,signal_fee,signal_fee_inv,plot_signal_r)

	
	#---DECONV Simple after DAQ 
	if deconv_simple_DAQ:

		print "++++++Deconv simple DAQ"
		DoDeconvSimpleDAQ(signal_PE_adc, signal_daq,signal_inv_daq, plot_signal_r_DAQ)

	

