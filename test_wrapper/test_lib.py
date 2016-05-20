# -*- coding: utf-8 -*-
"""
Created on Mon May  9 10:31:25 2016

@author: mablou
"""

import numpy as np

from normpdf import normpdf
#from BSSlib import getMinMax

#x = np.arange(1,100,1,dtype=float)
x = range(200)
y = normpdf(x,25,10)
import matplotlib.pyplot as plt

plt.plot(x,y)

#%%

from BSSlib import getMinMax
x= np.random.random_integers(1,1000,25)
minmax = getMinMax(x)
print minmax.min
print minmax.max