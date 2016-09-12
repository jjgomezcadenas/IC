# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 15:18:42 2016

@author: viherbos
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as pl

import os


def panel_to_hdf5(in_path, n_PMT, n_events, n_samples, out_path):

	# Reads signals from each event and PMT stored in text files
	# Creates HDF5 based on a single PANEL node where Item=PMT /
	# minor_axis(column)=Event / major_axis=sample

	if os.path.exists(out_path):
		os.remove(out_path)

	store = pd.HDFStore(out_path, complevel=9, complib='zlib')


	data = pd.Panel(dtype='int32',
			     items=range(0,n_PMT),
			     minor_axis=range(0,n_events),
                      major_axis=range(0,n_samples))
	# items: axis 0, each item corresponds to a DataFrame contained inside
	# major_axis: axis 1, it is the index (rows) of each of the DataFrames
	# minor_axis: axis 2, it is the columns of each of the DataFrames


	#The following reads all the data and stores it in a Panel
	for x in range(1, n_PMT+1):
		# Read all the files in the dataset
		for y in range(1, n_events+1):
			num = (x-1)*n_events+y
			path=''.join([in_path,str(num),'.txt'])
			print path
			# File to read
			g=pd.read_csv(path, names=[str(y-1)], header=None, dtype='int32')
			# DataFrame that stores read data and gives columns names.
			# Each column is an event
			data.loc[x-1,:,y-1]=g[str(y-1)]
			#Store data colum in its place

	store.put('data',data)
	# Dumps data to file
	store.close()


def read_panel_hdf5(in_path, PMT, event):

	#Reads PMT signal realted to a given event from HDF5 file

	a = pd.read_hdf(in_path,'data')

	b = np.array(a[PMT,:,event])

	return b


def get_panel_hdf5(in_path):

	#Reads PMT signal related to a given event from HDF5 file

	a = pd.read_hdf(in_path,'data')

	return a


def main():

	#SPE
	# in_path  = 'F:/DATOS_DAC/spe_1230/2046/pmt_0_trace_evt_'
	# out_path = 'spe_1230_2046.h5.z'
	# n_PMT   = 1
	# n_events = 5000
	# n_samples = 3200

	in_path  = 'D:/DATOS_DAC/2280/pmt_1_trace_evt_'
	out_path = 'spe_2280.h5.z'
	n_PMT   = 1
	n_events = 2000
	n_samples = 3200

	#Linearity
	# in_path  = 'D:/DATOS_DAC/2271/pmt_1_trace_evt_'
	# out_path = '2271.h5.z'
	# n_PMT   = 3
	# n_events = 50
	# n_samples = 6400

#	in_path  = 'D:/DATOS_DAC/Canfranc/Calibracion/run_1755'
#	out_path = 'Calibracion_1.h5.z'
#	n_PMT   = 12
#	n_events = 3
#	n_samples = 128000


	panel_to_hdf5(in_path, n_PMT, n_events, n_samples, out_path)

#	res = read_panel_hdf5(out_path, 0, 2)
#
#	pl.plot(res)
#	pl.show()

if __name__ == "__main__":
	main()