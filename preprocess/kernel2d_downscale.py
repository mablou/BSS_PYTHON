# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 09:54:47 2015

@author: blouinm
"""

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import math
from __future__ import division
from pylab import ginput
from scipy import interpolate
import pickle
#definition de fonctions utiles

def iqr_mb(x):
    q75, q25 = np.percentile(x, [75 ,25])
    return q75 - q25

#%%
# CONSTRUCTION DE LA MATRICE TOTALE

#Coordonnees X Y des puits en reference avec la sismique 3d
Inline_log = 641;
Xline_log = 400;

# On charge les donnees de puits 

inversion = np.loadtxt('/home/irsrvshare3/R16/XSE_BLOUIN/Barnett/SEISMIC/smallData.ixta')
coloc = (inversion[:,0]==640) * (inversion[:,1]==401)
IP_LOW = inversion[coloc,3]
t_low = inversion[coloc,2]
#IP_LOW=np.hstack([IP_LOW[0],IP_LOW])

IP = np.loadtxt('/home/irsrvshare3/R16/XSE_BLOUIN/Barnett/WELL/Davis_IP_low_def')[:,1]
IP_500 = np.loadtxt('/home/irsrvshare3/R16/XSE_BLOUIN/Barnett/WELL/Davis_IP_time_500us')[:,1]





t_500 = np.loadtxt('/home/irsrvshare3/R16/XSE_BLOUIN/Barnett/WELL/Davis_IP_time_500us')[:,0]
#t_low = np.hstack([602,t_low])

f=interpolate.interp1d(t_low,IP_LOW,kind='nearest')
IP_LOW = f(t_500)
#%%
#####################################################################
#Définition des vecteurs

min_IP=7000
max_IP=20000
min_PHI=0.10
max_PHI=0.5
# # Definition des vecteurs pour trouver les indices dans simul_bayes
vec_IP=np.linspace(min_IP,max_IP,500)
vec_PHI=np.linspace(min_PHI,max_PHI,500)


#####################################################################
#Definition des largeurs de bande
## Largeur de bande proposïee par Silverman (1986)


   
#==============================================================================
# 
# # Preparation de la grille pour le kernel
#==============================================================================



#%%
#test scipy
from scipy import stats

xmin =np.min(vec_IP)
xmax = np.max(vec_IP)
ymin =np.min(vec_IP)
ymax =np.max(vec_IP)


X, Y = np.mgrid[np.min(vec_IP):np.max(vec_IP):500j,np.min(vec_IP):np.max(vec_IP):500j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([IP_LOW, IP_500])
kernel = stats.gaussian_kde(values,bw_method=1)
Z = np.reshape(kernel(positions).T, X.shape)


fig, ax = plt.subplots(figsize=(5,5))
ax.imshow(np.rot90(Z), #cmap=plt.cm.gist_earth_r,
           extent=[xmin, xmax, ymin, ymax])
ax.plot(IP_LOW,IP_500, 'k.', markersize=4)
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
#ax.set_aspect(12000/0.4)
ax.set_xlabel('IP$_{inversion}$',fontsize=14)
ax.set_ylabel('IP$_{500}$',fontsize=14)
fig.show()
fig.savefig('kernel2d_downscale.png', bbox_inches='tight',dpi=600)
plt.show()


with open('stats_downscale.pickle', 'w') as f:
    pickle.dump([vec_IP,vec_PHI,Z], f)
    


#%%
pick_IP=(13250,13250)
Phi_line=(ymin,ymax)

fig,ax = plt.subplots(figsize=(5,5))

ax.plot(pick_IP,Phi_line,'k',linewidth=4)

ax.imshow(np.rot90(Z), #cmap=plt.cm.gist_earth_r,
           extent=[xmin, xmax, ymin, ymax])
ax.plot(IP_LOW,IP_500, 'k.', markersize=4)
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
#ax.set_aspect(12000/0.4)
ax.set_xlabel('IP$_{inversion}$',fontsize=14)
ax.set_ylabel('IP$_{500}$',fontsize=14)

fig.savefig('kernel2d_downscale_w_pick.png', bbox_inches='tight',dpi=600)
#%%

fig,ax = plt.subplots(figsize=(5,5))

ax.plot(vec_IP,Z[240,:])
ax.set_xlabel('IP$_{500}$',fontsize=14)
ax.yaxis.set_ticks([])
ax.set_ylabel('Probabilit'u'é',fontsize=14)

fig.savefig('downscale_function_ex.png', bbox_inches='tight',dpi=600)

#%%
xmin =np.min(vec_IP)
xmax = np.max(vec_IP)
ymin =np.min(vec_PHI)
ymax =np.max(vec_PHI)

PHI_500 = np.loadtxt('/home/irsrvshare3/R16/XSE_BLOUIN/Barnett/WELL/Davis_NPOR_time_500us')[:,1]
PHI_500=PHI_500[1:]

X, Y = np.mgrid[xmin:xmax:500j,ymin:ymax:500j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([IP_500, PHI_500])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)


fig, ax = plt.subplots(figsize=(5,5))
ax.imshow(np.rot90(Z), #cmap=plt.cm.gist_earth_r,
           extent=[xmin, xmax, ymin, ymax])
ax.plot(IP_500,PHI_500, 'k.', markersize=4)
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.set_aspect(12000/0.4)
ax.set_xlabel('IP$_{500}$',fontsize=14)
ax.set_ylabel('PHI$_{500}$',fontsize=14)
fig.show()
fig.savefig('stats_pdf2d_500us.png', bbox_inches='tight',dpi=600)
plt.show()


with open('stats_pdf2d_500us.pickle', 'w') as f:
    pickle.dump([vec_IP,vec_PHI,Z], f)
    
#%%
PHI_HIGH = np.loadtxt('/home/irsrvshare3/R16/XSE_BLOUIN/Barnett/WELL/Davis_NPOR_high_def')[:,1]
stats.entropy(PHI_HIGH)
stats.entropy(PHI_500)
stats.entropy(PHI)
