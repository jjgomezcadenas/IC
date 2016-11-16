# -*- coding: utf-8 -*-
"""
Created on Thu Jun 02 23:44:40 2016

@author: viherbos
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



# Fitting functions definition
def line(x, A, B):
    return A*x + B

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x-mu)**2/(2.*sigma**2))

def gauss2(x, *param):
    return param[0] * np.exp(-(x-param[1])**2/(2.*param[2]**2)) + \
           param[3] * np.exp(-(x-param[4])**2/(2.*param[5]**2))

def gauss3(x, *param):
    return param[0] * np.exp(-(x-param[1])**2/(2.*param[2]**2)) + \
           param[3] * np.exp(-(x-param[4])**2/(2.*param[5]**2)) + \
           param[6] * np.exp(-(x-param[7])**2/(2.*param[8]**2))


# The following functions implement several fitting flavors

def line_fit(f,X,f_sigma,x_text,y_text,title_text,n_figure,graph_sw):

# Makes a linear fit for n points (X input vector). 
# f is the mean of the measured data points
# f_sigma is the standard deviation (ddof=1) of the measured data points
# The rest are attributes for the plotting windows (graph_sw = 1 to plot)
# returns coeff (A,B), perr - error for the fit param, 
#         XI2_r --> Squared CHI reduced (Goodness of fit) 

    p0 = [1,(f[1]-f[0])/(X[1]-X[0])]
    coeff, var_matrix = curve_fit(line, X, f,p0=p0)

    #Parameters error estimation (sigma). See numpty user guide
    perr = np.sqrt(np.diag(var_matrix))

    Y_fit = line(X,coeff[0],coeff[1])

    XI2 = np.sum(((Y_fit-f)**2.)/(f_sigma**2.))
    XI2_r = XI2/(len(X)-2)

    if (graph_sw==1):
    # Draws figure with all the properties

        plt.figure(n_figure)
        plt.plot(X, Y_fit, 'r--', linewidth=1)
        plt.errorbar(X, f, fmt='b*', yerr=f_sigma)
        plt.xlabel(x_text)
        plt.ylabel(y_text)
        plt.title(title_text)
        plt.figtext(0.2,0.75, ('CHI2_r = %0.3f ' % (XI2_r)))
        plt.show(block=False)
        #Fit parameters
    print ('Fitted A = ', coeff[0], '( Error_std=', perr[0],')')
    print ('Fitted B = ', coeff[1], '( Error_std=', perr[1],')')

    return coeff, perr, XI2_r


def gauss1_fit(f,x_text,y_text,title_text,bins,n_figure,graph_sw):
	
# Makes a gaussian fit for histogram of the measurements. 
# f is the measurements vector
# The rest are attributes for the plotting windows (graph_sw = 1 to plot)
# returns coeff (see gaussian definition), perr - error for the fit param
	
    hist, bin_edges = np.histogram(f, bins=bins)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

    p0 = [1, np.mean(f), np.std(f)]
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

    #Gets fitted function and residues
    hist_fit = gauss(bin_centres, coeff[0],coeff[1],coeff[2])

    #Parameters error eestimation (sigma). See numpty user guide
    perr = np.sqrt(np.diag(var_matrix))

    if (graph_sw==1):
    # Draws figure with all the properties
	
        plt.figure(n_figure)
        plt.hist(f, bins, facecolor='green')
        plt.plot(bin_centres, hist_fit, 'r--', linewidth=1)
        plt.grid(True)
        plt.xlabel(x_text)
        plt.ylabel(y_text)
        plt.title(title_text)
        plt.figtext(0.2,0.75, ('Fitted MU = %0.3f ( Error_std = %0.3f)' % (coeff[1] , perr[1])))
        plt.figtext(0.2,0.7, ('Fitted SIGMA = %0.3f ( Error_std = %0.3f)' % (coeff[2] , perr[2])))
        plt.show(block=False)
        #Fit parameters
        print ('Fitted A = ', coeff[0], '( Error_std=', perr[0],')')
        print ('Fitted MU = ', coeff[1], '( Error_std=', perr[1],')')
        print ('FItted SIGMA = ', coeff[2], '( Error_std=', perr[2],')')

    return coeff, perr


def gauss2_fit(f,x_text,y_text,title_text,\
               bins,mu_guess,\
               n_figure,graph_sw):

# Makes a two gaussian fit for histogram of the measurements. 
# f is the measurements vector
# mu_guess: (required for precise fitting)
# The rest are attributes for the plotting windows (graph_sw = 1 to plot)
# returns coeff (see gaussian definition), perr - error for the fit param

    hist, bin_edges = np.histogram(f, bins=bins)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

    p0 = (1, mu_guess[0], 1, 1, mu_guess[1], 1)
    coeff, var_matrix = curve_fit(gauss2, bin_centres, hist, p0=p0)

    #Gets fitted function and residues
    hist_fit = gauss2(bin_centres,*coeff)

    #Parameters error eestimation (sigma). See numpty user guide
    perr = np.sqrt(np.diag(var_matrix))

    if (graph_sw==1):
    # Draws figure with all the properties
        plt.figure(n_figure)
        plt.hist(f, bins, facecolor='green')
        plt.plot(bin_centres, hist_fit, 'r--', linewidth=1)
        plt.grid(True)
        plt.xlabel(x_text)
        plt.ylabel(y_text)
        plt.title(title_text)
        plt.figtext(0.2,0.8, ('Fitted MU1 = %0.3f ( Error_std = %0.3f)' % (coeff[1] , perr[1])))
        plt.figtext(0.2,0.75, ('Fitted MU2 = %0.3f ( Error_std = %0.3f)' % (coeff[4] , perr[4])))
        plt.show(block=False)
        #Fit parameters
        print ('Fitted A1 = ', coeff[0], '( Error_std=', perr[0],')')
        print ('Fitted MU1 = ', coeff[1], '( Error_std=', perr[1],')')
        print ('FItted SIGMA1 = ', coeff[2], '( Error_std=', perr[2],')')
        print ('Fitted A2 = ', coeff[3], '( Error_std=', perr[3],')')
        print ('Fitted MU2 = ', coeff[4], '( Error_std=', perr[4],')')
        print ('Fitted SIGMA2 = ', coeff[5], '( Error_std=', perr[5],')')

    return coeff, perr


def gauss3_fit(f,x_text,y_text,title_text,\
               bins,mu_guess,\
               n_figure,graph_sw):
# Makes a three gaussian fit for histogram of the measurements. 
# f is the measurements vector
# mu_guess: (required for precise fitting)
# The rest are attributes for the plotting windows (graph_sw = 1 to plot)
# returns coeff (see gaussian definition), perr - error for the fit param
															
    hist, bin_edges = np.histogram(f, bins=bins)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

    p_guess = (1, mu_guess[0], 1, 1, mu_guess[1], 1, 1, mu_guess[2], 1)
    coeff, var_matrix = curve_fit(gauss3, bin_centres, hist, p0=p_guess)

    #Gets fitted function and residues
    hist_fit = gauss3(bin_centres,*coeff)

    #Parameters error estimation (sigma). See numpty user guide
    perr = np.sqrt(np.diag(var_matrix))

    if (graph_sw==1):
    # Draws figure with all the properties

        plt.figure(n_figure)
        plt.hist(f, bins, facecolor='green')
        plt.plot(bin_centres, hist_fit, 'r--', linewidth=1)
        plt.grid(True)
        plt.xlabel(x_text)
        plt.ylabel(y_text)
        plt.title(title_text)
        plt.figtext(0.2,0.8, ('Fitted MU1 = %0.3f ( Error_std = %0.3f)' % (coeff[1] , perr[1])))
        plt.figtext(0.2,0.75, ('Fitted MU2 = %0.3f ( Error_std = %0.3f)' % (coeff[4] , perr[4])))
        plt.figtext(0.2,0.7, ('Fitted MU3 = %0.3f ( Error_std = %0.3f)' % (coeff[7] , perr[7])))
        plt.show(block=False)
        #Fit parameters
        print ('Fitted A1 = ', coeff[0], '( Error_std=', perr[0],')')
        print ('Fitted MU1 = ', coeff[1], '( Error_std=', perr[1],')')
        print ('FItted SIGMA1 = ', coeff[2], '( Error_std=', perr[2],')')
        print ('Fitted A2 = ', coeff[3], '( Error_std=', perr[3],')')
        print ('Fitted MU2 = ', coeff[4], '( Error_std=', perr[4],')')
        print ('FItted SIGMA2 = ', coeff[5], '( Error_std=', perr[5],')')
        print ('Fitted A3 = ', coeff[6], '( Error_std=', perr[6],')')
        print ('Fitted MU3 = ', coeff[7], '( Error_std=', perr[7],')')
        print ('Fitted SIGMA3 = ', coeff[7], '( Error_std=', perr[8],')')
    return coeff, perr