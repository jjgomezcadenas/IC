
# coding: utf-8

# # CoreUtil nb

# ## This nb bundles together core utilitiy functionality such as job configuration

# #### Logger

# In[22]:

# import logging
# import sys
# logger = logging.getLogger()
# logger.handlers[0].stream = sys.stdout
# logger.setLevel(logging.DEBUG)

import logging
import sys

logger = logging.getLogger()
# create console handler
ch = logging.StreamHandler()
logger.addHandler(ch)
logger.handlers[0].stream = sys.stdout
logger.setLevel(logging.INFO)


# #### Pandas is used to read the csv file

# In[25]:

import pandas as pd


def cdf_to_dict(cdf):
    """
    transforms the configuration data frame into a dictionary
    """

    dc ={}
    for k in cdf.keys():
        dc[k] = cdf[k][0]
    return dc
    

def configure(cfile, INFO=False, level='INFO'):

    """
    Configures job
    """
    
    lg = 'logging.'+level
    logger.setLevel(eval(lg))

    if cfile == '':
        print("Path to configuration file not given. Please specify path")
        
    cfp =pd.read_csv(cfile,comment="#")
    CFP = cdf_to_dict(cfp)
    
    logger.info("Configuration Parameters (CFP) dictionary  = {}".format(CFP))
    return INFO, CFP
    
def farray_from_string(sfl):
    """
    returns a np array of floats from a string (sfl)
    representing the array in which the numbers are separated by blanks
    e.g, '1.4 3.6 6.7' -->np.array(1.4,3.6,6.7)
    """
    sarr = sfl.split(' ')
    arr = []
    for x in sarr:
        arr.append(float(x))
    return np.array(arr)

def rebin_array(a, stride):
    """
    rebins the array according to stride 
    """
    lenb = len(a)/int(stride)
    b = np.zeros(lenb)
    j=0
    for i in range(lenb):
        b[i] = np.sum(a[j:j+stride])
        j+= stride
    return b



