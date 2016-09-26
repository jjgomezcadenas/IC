# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 00:24:52 2016

@author: viherbos
"""

import numpy as np
import matplotlib.pyplot as plt
from BLR import BLR
from panel_to_hdf5 import read_panel_hdf5
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import pandas as pd
import FEParam as FP



thr_ener_delta = 50



def read_panel_hdf5(in_path, PMT, event):

	#Reads PMT signal realted to a given event from HDF5 file

	a = pd.read_hdf(in_path,'data')

	b = np.array(a[PMT,:,event])

	return b


def find_coeff_II(  LIMIT_L, LIMIT_H, PULSE_R,
			hdf5_file, Channel, event_range):
	
	# LIMIT_L / LIMIT_H : Interesting part of the signal
	# hdf5_file, Channel, event  : Storage file, Channel and Event

	thr3_a = 0.001*FP.NOISE_ADC

	# No fast return to baseline
	
	############## AUXILIARY TOOLS #########################################
	
	def BLR_pulse(x, coef, f, X_x):
		# Applies BLR reconstruction function to f signal using coef
		signal = BLR(signal_daq=(4096-f[X_x].flatten().astype(float)), 
									coef=coef, mau_len=250,	thr1 = 1.5*FP.NOISE_ADC, thr2 = 0, thr3= thr3_a)
		return signal.signal_r[x]
		

	def find_baseline(x):
	    # Finds baseline in a given sequence 
		length_signal = np.size(x)
		baseline = x.sum() / length_signal
		
		return baseline	
	
	########################################################################	
	

	signal_p = read_panel_hdf5(hdf5_file, Channel, event_range)
	signal = signal_p[:,0]


	# Mean of all the signals to minimize random noise due to LED fluctuations
	for m in event_range[2:]:
		
		signal = signal_p[:,m] + signal 
	
	signal = signal / np.size(event_range) 		


	# Signal baseline
	baseline 	   = 4096-find_baseline(signal[0:PULSE_R])					
	print ('BASELINE= ', baseline)	

    # X axis data						
	X = np.array(range(0,np.size(signal)))
	
	# Function for LS minimization
	func = lambda coef,x: BLR_pulse(x, coef, signal, X)
	# Error function for LS minimization
	ErrorFunc = lambda coef,x: np.std(func(coef,x),ddof=0)

	p0 = 1.4E-3
	coeff,cov,infodict,mesg,ier = leastsq(ErrorFunc,p0,
										args=(range(LIMIT_L,LIMIT_H,1)),full_output=True)
	ss_err=(infodict['fvec']**2).sum()
	

# We look for no negative lobe in the failing edge of the pulse and 
# a minimum total error in the signal after the falling edge and a stab. period for LED
# We take the standard deviation (or variance) of the tail of the reconstructed signal 
# as the parameter to minimiza because the final value is unknown (LED fluctuation error). 
# When the std is minimum the tail is almost flat which reproduces
# the ideal behavior of the BLR.
# WARNING: For long decays due to errors in the DBLR we seek a minimum std at the 
# beginning of the tail just after the failing edge and the LED transient


	out1 = BLR(signal_daq=(4096-signal[X].flatten().astype(float)),
							coef=coeff,
							mau_len=250,
							thr1 = 1.5*FP.NOISE_ADC, thr2 = 0, 
							thr3 = thr3_a)

	out2 = BLR(signal_daq=(4096-signal[X].flatten().astype(float)),
							coef=coeff*1.1,
							mau_len=250,
							thr1 = 1.5*FP.NOISE_ADC, thr2 = 0, 
							thr3 = thr3_a)

	out3 = BLR(signal_daq=(4096-signal[X].flatten().astype(float)),
							coef=coeff*0.9,
							mau_len=250,
							thr1 = 1.5*FP.NOISE_ADC, thr2 = 0, 
							thr3 = thr3_a)

	# This is to keep an eye on the evolution of the error (+-2.5% change in coeff)

	plt.figure()
	plt.plot(X, 4096-signal)			
	plt.plot(X, out1.signal_r, 'r--', linewidth=1)
	plt.plot(X, out2.signal_r, 'b--', linewidth=2)
	plt.plot(X, out3.signal_r, 'b--', linewidth=2)
	plt.show()

	thr_ener = thr_ener_delta

	ener1=(out1.signal_r>thr_ener)*out1.signal_r
	ener2=(out2.signal_r>thr_ener)*out2.signal_r
	ener3=(out3.signal_r>thr_ener)*out3.signal_r

	print (100*(ener2.sum()-ener1.sum())/ener1.sum(), ener1.sum(), 100*(ener3.sum()-ener1.sum())/ener1.sum())

	# This is to get an idea of the coeff error effect on energy calculations (10% error on coeef --> 1% energy error)
	# This coeff error is absolutely linear so it is calibrable away

	return coeff


def main():

	coeff=np.array(np.zeros((24,5)))
   	
   	# WATCH OUT!!!    LIMIT_L --> Pulse finished and LED fluctuations extinguished
	# 				  LIMIT_H --> Different strategies depending on pulse length.
	#							for short (<=10u) almost no DC error so minimize std
	#							of the whole tail. For long pulses the DC error introduces
	#							a long decay -> make decay as flat (linear) as posible
	#							minimizing std of the "bump" zone

	# for i in range(2,24,1): 
	#    	LIMIT_L       = 1870
	#    	LIMIT_H       = 8000
	#    	PULSE_R  	 = 1738
	# 	PULSE_L 	 = 40
	#    	hdf5_file 	 = 'F:\DATOS_DAC\CALIBRATION\cal_1u.h5.z'
	#    	Channel 		     = i
	#    	event_range      = range(0,500,1)

	# 	coeff[i,0] = find_coeff_II( LIMIT_L, LIMIT_H, PULSE_R, PULSE_L,
	# 								hdf5_file, Channel, event_range )	
	# 	print (coeff)


	# for i in range(2,24,1): 
	#    	LIMIT_L       = 1980
	#    	LIMIT_H       = 8000
	#    	PULSE_R  	 = 1735
	#		PULSE_L 	 = 1600
	#    	hdf5_file 	 = 'F:\DATOS_DAC\CALIBRATION\cal_2u5.h5.z'
	#    	Channel 		     = i
	#    	event_range      = range(0,500,1)

	# 	coeff[i,1] = find_coeff_II( LIMIT_L, LIMIT_H, PULSE_R, PULSE_L,
	# 								hdf5_file, Channel, event_range )	
	# 	print (coeff)


	# for i in range(2,24,1): 
	#    	LIMIT_L       = 2230
	#    	LIMIT_H       = 8000
	#    	PULSE_R  	 = 1736
	#    	PULSE_L 	 = 400
	#    	hdf5_file 	 = 'F:\DATOS_DAC\CALIBRATION\cal_10u.h5.z'
	#    	Channel 		     = i
	#    	event_range      = range(0,500,1)

	# 	coeff[i,2] = find_coeff_II( LIMIT_L, LIMIT_H, PULSE_R, PULSE_L,
	# 								hdf5_file, Channel, event_range )	
	# 	print (coeff)
	   	

	for i in range(23,24,1): 
	   	LIMIT_L       = 3850
	   	LIMIT_H       = 4250
	   	PULSE_R  	 = 1736
		PULSE_L 	 = 2000
	   	hdf5_file 	 = '/mnt/WINDOWS_ntfs/DATOS_DAC/CALIBRATION/cal_50u.h5.z'
	   	Channel 		     = i
	   	event_range      = range(0,500,1)

		coeff[i,3] = find_coeff_II( LIMIT_L, LIMIT_H, PULSE_R,
									hdf5_file, Channel, event_range )	
		print (coeff)


	# for i in range(2,24,1): 
	#    	LIMIT_L       = 5820
	#    	LIMIT_H       = 16000
	#    	PULSE_R  	 = 1736
	# 	PULSE_L 	 = 4000
	#    	hdf5_file 	 = 'F:\DATOS_DAC\CALIBRATION\cal_100u.h5.z'
	#    	Channel 		     = i
	#    	event_range      = range(0,500,1)

	# 	coeff[i,4] = find_coeff_II( LIMIT_L, LIMIT_H, PULSE_R, PULSE_L,
	# 								hdf5_file, Channel, event_range )	
	#	print (coeff)


	#aux_frame = pd.DataFrame(coeff)
	

if __name__ == "__main__":
	main()
		


