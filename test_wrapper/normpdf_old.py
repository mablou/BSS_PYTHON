# -*- coding: utf-8 -*-
"""
Created on Fri May  6 13:40:46 2016

@author: mablou
"""

import ctypes
import numpy as np

_Bayeslib = ctypes.CDLL('libBSS_library.so')

_Bayeslib.normpdf.argtypes = (ctypes.POINTER(ctypes.c_double),
                                  ctypes.c_int,ctypes.c_double, 
                                      ctypes.c_double, 
                                          ctypes.POINTER(ctypes.c_double))
_Bayeslib.restype = ctypes.POINTER(ctypes.c_double)

def normpdf(in_array,mean,stddev):
    global _Bayeslib    
    n_array = len(in_array)
    out_array = np.zeros(n_array)
    array_type = ctypes.c_double * n_array
    toto=_Bayeslib.normpdf(array_type(*in_array),ctypes.c_int(n_array),
                          ctypes.c_double(mean),ctypes.c_double(stddev),
                                                    array_type(*out_array))
    return toto
    