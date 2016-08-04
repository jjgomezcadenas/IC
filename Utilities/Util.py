"""
Utilities
"""
#from __future__ import print_function
from math import *
from system_of_units import *
import numpy as np

def wait():
	raw_input("Press a key...")


def drange(start, stop, step):
	"""
	a range of doubles
	"""
	r = start
	while r < stop:
		yield r
		r += step

def lin(x,x0,y0,x1,y1):
	"""
	lineal extrapolation
	"""
	
	return y0 + (x-x0)*(y1-y0)/(x1-x0)

def inRange(x,xmin,xmax):
	if x >= xmin and x <= xmax:
		return True
	else:
		return False

def differenceRatio(x1,x2):
	return 2*abs(x1-x2)/(x1+x2)

def ResFWHM(sigma,E):
	return 100*2.35*sigma/E

def MeanWOL(x,nmean=5, nsigma=5):
	"""
	Mean Without OutLayers
	"""
	x_mean = np.average(x[0:nmean])
	x_std = np.std(x[0:nmean])
	xmin = x_mean - 5*x_std
	xmax = x_mean + 5*x_std

	x2 = []
	for y in np.nditer(x):
		if inRange(y,xmin,xmax) == True:
			x2.append(y)
            
	return(np.array(x2))


