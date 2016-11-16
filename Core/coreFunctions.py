"""
Core functions
"""
import numpy as np
import pandas as pd


def wait():
    """
    A simple convenience name for raw_input
    """
    raw_input("Press a key...")


def dict_map(F, D):
    '''
    Apply map to dictionary values maintaining correspondence.
    '''
    return {key: F(val) for key, val in D.iteritems()}


def df_map(F, DF, field):
    '''
    Apply map to some DataFrame field.
    '''
    out = pd.DataFrame(DF)
    out[field] = map(F, out[field])
    return out


def dict_filter(C, D):
    '''
    Apply filter to dictionary values without losing correspondence.
    '''
    return {key: val for key, val in D.iteritems() if C(val)}


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
    j = 0
    for i in range(lenb):
        b[i] = np.sum(a[j:j+stride])
        j += stride
    return b


def drange(start, stop, step):
    """
    a range of doubles
    """
    r = start
    while r < stop:
        yield r
        r += step
