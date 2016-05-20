# -*- coding: utf-8 -*-
"""
Created on Tue May 10 13:44:54 2016

@author: mablou
"""
#%%
import os
import numpy as np
from scipy.io import loadmat
import BSSlib as bs

os.chdir('/Users/mablou/Dropbox/postdoc_IFP/BSS_PYTHON')
#%%
well_data = np.loadtxt('data/Davis_IP_low_def')
seismic = np.loadtxt('data/smallData.ixta')
pdf2d = loadmat('preprocess/pdf2d.mat')
stat_poro = loadmat('preprocess/stat_poro.mat')
prob = loadmat('preprocess/prob.mat')
#%%
#Number of structures
NVario = 2
#Normalized variance for each structure
Variance = [0.95,0.05]
#Type of structures
Vario = [3,8]
alpha = [1,1]
#vario range
scf = [10,10,10]

#We define axes for anisotropy
#No anisotropy in this case
ap = np.array([np.diag([1,1,1])for vario in range(NVario)]).flatten()

variogram = dict((('Nvario',NVario),('vario',Vario),
          ('alpha',alpha),('ap',ap),('scf',scf),('var',Variance)))
#%%
         