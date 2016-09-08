# -*- coding: utf-8 -*-
"""
Created on Thu Jun 02 23:44:40 2016

@author: viherbos
"""
import numpy as np
import pandas as pd
from panel_to_hdf5 import get_panel_hdf5
from fit_library import gauss2_fit


def find_SPE_h5(base_path, n_files, start, end, bins, guess, n_figure):
	# This is a basic function build to find the value of the SPE.
	# PARAMETERS:
	#base_path: Base path including the name of the file, number not included (see main)
	#n_files:   Number of files in the dataset
	#start:     Starting Point of the integral
	#end:       End point of the integral
	#bins:      Number of bins of the histogram
	#guess:     Initial guess for the SPE
	#n_figure:  figure number (in case many windows are required)

	# The function returns the fit parameters:
	# a[0] = A1 // a[1] = MU1 // a[2] = SIGMA1 //
	# a[3] = A2 // a[4] = MU2 // a[5] = SIGMA2
	# b --> Error (Sigma) of the fit parameters


	#GRAPHICS WINDOW
	x_text = 'ADC_COUNTS (LSBs)'
	y_text = 'Hits'
	title_text = 'INTEGRATED SINGLE PHOTOELECTRON VALUE'

	integral_r = np.zeros(n_files,float)

	# Party time !!!
	g = get_panel_hdf5(base_path)

	for x in range(1, n_files+1):

		#f = read_panel_hdf5(base_path,0,x-1)
		f = g[0,:,x-1]
		media=np.mean(f)
		f = f[start:end] - media
		#Integrate the SPE (beware of the sampling period)
		integral_r[x-1]=np.sum(f)
		print media

	# Fit two gaussians (one for 0 and the other for the SPE (Initial guess required)
	a,b = gauss2_fit(integral_r, x_text, y_text, title_text, bins, [0,guess], n_figure, 1)

	return a,b


def main():
	base_path = 'spe_1230_2046.h5.z'
	n_files   = 500
	#INTEGRATING RANGE
	start=(43.3*40)
	end  =(43.8*40)
	bins = 200
	guess = -20
	n_figure = 1

	a,b=find_SPE_h5(base_path, n_files, start, end, bins, guess, n_figure=1)

if __name__ == "__main__":
	main()