# -*- coding: utf-8 -*-
"""
Created on Fri May  6 13:40:46 2016

@author: mablou
"""

#import ctypes
import numpy as np

_Bayeslib = np.ctypeslib.load_library('libBSS_library.so','.')

_Bayeslib.normpdf.argtypes = (np.ctypeslib.ndpointer(float),
                                      np.ctypeslib.ctypes.c_int,
                                      np.ctypeslib.ctypes.c_double,
                                      np.ctypeslib.ctypes.c_double,
                                      np.ctypeslib.ndpointer(float))

def normpdf(in_array,mean,stddev):
    in_array = np.require(in_array,float)
    out_array = np.empty_like(in_array)
    _Bayeslib.normpdf(in_array,in_array.size,mean,stddev,out_array)
    return out_array

    