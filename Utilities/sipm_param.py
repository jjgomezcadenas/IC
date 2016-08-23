"""
sipm_param.py
author: jrenner

Defines the SiPM parameterization functions as:

N(x) = M*sum(c_n*x^n) for n = 0 to n = 9

where x is the distance of the SiPM from some central point of
light emission.

Because the response is characterized over several time bins, we
have several values for M and the coefficients.

"""
import numpy as np

# Number of time bins
n_tbins = 2

# Coefficients from S2 parameterization
M = [1.599, 1.599]
c0 = [7.72708346764e-05, 0.000116782596518]
c1 = [-1.69330613273e-07, 3.05115354927e-06]
c2 = [-1.52173658255e-06, -7.00800605142e-06]
c3 = [-2.4985972302e-07, 6.53907883449e-07]
c4 = [1.12327204397e-07, 8.95230202525e-08]
c5 = [-1.49353264606e-08, -2.27173290582e-08]
c6 = [1.04614146487e-09, 2.00740799864e-09]
c7 = [-4.19111362353e-11, -9.21915945523e-11]
c8 = [9.12129133361e-13, 2.20534216312e-12]
c9 = [-8.40089561697e-15, -2.1795164563e-14]

# Maximum radial extent of parameterization
rmax = 20.

# Return the SiPM response for the specified time bin and radial distance.
def sipm_par(tbin,r):

    # Ensure the time bin value is valid.
    if(tbin < 0 or tbin >= n_tbins):
        print "Invalid time bin in sipm_param: returning 0.0 ..."
        return 0.0

    # Calculate the response based on the parametrization.
    vpar = M[tbin]*(c0[tbin] + c1[tbin]*r + c2[tbin]*r**2 + c3[tbin]*r**3 + 
    c4[tbin]*r**4 + c5[tbin]*r**5 + c6[tbin]*r**6 + c7[tbin]*r**7 + 
    c8[tbin]*r**8 + c9[tbin]*r**9)

    # Zero the response for radii too large.
    if(hasattr(vpar, "__len__")):
        ret = np.zeros(len(vpar)); iret = 0
        for rv,pv in zip(r,vpar):
            if(rv < rmax):
                ret[iret] = pv
            iret += 1
        return ret
    else:
        if(r < rmax):
            return vpar
        return 0.0
        
        
