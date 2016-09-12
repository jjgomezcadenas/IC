# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 00:24:52 2016

@author: viherbos
"""

import numpy as np
import matplotlib.pyplot as plt
import DBLR
from panel_to_hdf5 import read_panel_hdf5
from scipy.optimize import curve_fit







def find_coeff(  LIMIT_L, LIMIT_H, PULSE_R, PULSE_L,
			pulse_height,
			hdf5_file, PMT, event,
			amplitude_range, delta,
             noise_sigma,
			draw):
	
	# LIMIT_L / LIMIT_H : Interesting part of the signal
	# pulse_height      : Pulse Amplitude
	# hdf5_file, PMT, event  : Storage file, PMT and Event
	# amplitude_range, delta : Range and Delta for the sweep loop of amplitude
	# draw              : Draw all the signals in the loop for debugging 
	
	
	############## AUXILIARY TOOLS #########################################
	
	def BLR_pulse(x, coef, f):
		# Applies BLR reconstruction function to f signal using coef
		
		#f = read_panel_hdf5(hdf5_file, PMT, event)
		signal_r,energy = DBLR.BLR(signal_daq=(4096-f[x].flatten().astype(float)),
							coef=coef,
							n_sigma = noise_sigma, NOISE_ADC=0.75,
							thr1 = 0, thr2 = 0, thr3 = 0, plot = False)
		return signal_r
		
	def pulse(x, pulse_init, pulse_length,baseline,amplitude,error):
		# Generates a pulse signal with all the parameters
		signal = np.array(np.zeros(x.shape[0]))
		signal[pulse_init:(pulse_init+pulse_length)]=1
	
		signal = signal*amplitude+baseline
		
		signal[(pulse_init+pulse_length):]=signal[(pulse_init+pulse_length):]+error
	
		return signal
	
	def find_baseline(x):
	    # Finds baseline in a given sequence 
		length_signal = np.size(x)
		baseline = x.sum() / length_signal
		
		return baseline	
	
	########################################################################	
	
	j=0
	
	signal	   = read_panel_hdf5(hdf5_file, PMT, event)		
	baseline 	   = 4096-find_baseline(signal[0:PULSE_R])
				
	coef_array    = np.array(np.zeros(amplitude_range/delta))    
	std_array 	  = np.array(np.zeros(amplitude_range/delta))    
	std_array_X   = np.array(np.zeros(amplitude_range/delta))
					
					
	print ('BASELINE= ', baseline)			
    	
	X = np.array(range(0,np.size(signal)))
	# X axis data					

	for i in np.arange(pulse_height,pulse_height+amplitude_range,delta):
		
		
		pulse_signal = pulse(  x=X,
						pulse_init   = PULSE_R,
						pulse_length = PULSE_L,
						baseline     = baseline,
						amplitude    = i,
						error        = 0
						)
		# Ideal response signal defined
				
		X_x = X[LIMIT_L:LIMIT_H]
		f=pulse_signal[LIMIT_L:LIMIT_H]
	
		p0 = [1.4E-3]
		coeff, var_matrix = curve_fit(lambda x, coef: BLR_pulse(x, coef, signal), 
							X_x, f, p0=p0, bounds=([1E-3], [1.8E-3]))
		perr = np.sqrt(np.diag(var_matrix))
		
				
		Y_fit = BLR_pulse(X_x,coeff,signal)
		
				
		error_total = Y_fit[PULSE_R+PULSE_L-LIMIT_L:LIMIT_H-LIMIT_L]	
	
		standard_dev = error_total.std(ddof=1)		
		
		print (i,'  ', coeff[0],'  ', perr[0],' ',standard_dev)
	
# We look for no negative lobe in the failing edge of the pulse and 
# a minimum total error in the whole signal after the falling edge
# We take the standard deviation of the tail of the reconstructed signal 
# as the parameter to optimize vs the signal amplitude. 
# When the std is minimum the tail is almost flat which reproduces
# the ideal behavior of the BLR 

		std_array[j]   = standard_dev
		std_array_X[j] = i
		coef_array[j]  = coeff[0]
		
		j=j+1
	
		if (draw==True):
			plt.figure(j)
			plt.plot(X_x,f)
			plt.plot(X_x, Y_fit, 'r--', linewidth=1)
			plt.show()
		# Plots all the signals for debugging purposes


	plt.figure(j+1)
	plt.plot(std_array_X,std_array)
	coeff = coef_array[std_array.argmin()]

	print (coeff, std_array_X[std_array.argmin()], std_array[std_array.argmin()])

	plt.figure(j+2)
	plt.plot(X,pulse(  x=X,
				pulse_init   = PULSE_R,
				pulse_length = PULSE_L,
				baseline     = baseline,
				amplitude    = std_array_X[std_array.argmin()],
				error        = 0
				))
				
	plt.plot(X, BLR_pulse(X,coef_array[std_array.argmin()],signal), 'r--', linewidth=1)
	#plt.show()
	
	return coeff


def main():

    #Params
	
    draw = False
 
    LIMIT_L       = 19000
    LIMIT_H       = 22500
    PULSE_R       = 20142
    PULSE_L       = 1200
    pulse_height  = 545
    hdf5_file 	 = '2052.h5.z'
    PMT 		     = 0
    event          = 0
    amplitude_range = 2
    delta          = 0.1
    noise_sigma    = 4


    # LIMIT_L       = 0
    # LIMIT_H       = 6400
    # PULSE_R       = 3340
    # PULSE_L       = 400
    # pulse_height  = 300
    # hdf5_file 	 = '2271.h5.z'
    # PMT 		     = 0
    # event          = 5
    # amplitude_range = 100
    # delta          = 1
    # noise_sigma    = 4

    coeff=find_coeff(  LIMIT_L, LIMIT_H, PULSE_R, PULSE_L,
			pulse_height,
			hdf5_file, PMT, event,
			amplitude_range, delta,
             noise_sigma,
			draw)			
    plt.show() 

if __name__ == "__main__":
	main()
		


