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
import pickle

#function to calculate interquartile range

def iqr_mb(x):
    q75, q25 = np.percentile(x, [75 ,25])
    return q75 - q25

#%%
# Load well data Porosity and IP

IP = np.loadtxt('../data/Davis_IP_low_def')[:,1]
PHI = np.loadtxt('../data/Davis_NPOR_low_def')[:,1]

#%%
#First method to construct the kernel with scipy
from scipy import stats

xmin =np.min(vec_IP)
xmax = np.max(vec_IP)
ymin =np.min(vec_PHI)
ymax =np.max(vec_PHI)


X, Y = np.mgrid[np.min(vec_IP):np.max(vec_IP):500j,np.min(vec_PHI):np.max(vec_PHI):500j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([IP, PHI])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)


import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.imshow(np.rot90(Z), #cmap=plt.cm.gist_earth_r,
           extent=[xmin, xmax, ymin, ymax])
ax.plot(IP, PHI, 'k.', markersize=2)
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.set_aspect(12000/0.4)
plt.show()

#%%

#####################################################################




#####################################################################
#Bandwidth defintion
## Using Silverman (1986)

l_IP1 = 0.9 * (min(np.std(IP),iqr_mb(IP)/1.349)) * len(IP)**(-1/5)
l_PHI1 = 0.9 * (min(np.std(PHI),iqr_mb(PHI)/1.349)) * len(PHI)**(-1/5)

## Bowman et Foster
l_IP2=np.std(IP)*len(IP)**(-1/6)
l_PHI2=np.std(PHI)*len(PHI)**(-1/6)


## Thumb rule (Deheurels (1977))
l_IP3=1.06*np.std(IP)* len(IP)**(-1/5)
l_PHI3=1.06*np.std(PHI)* len(PHI)**(-1/5)

print('l_IP')
print(l_IP1)
print(l_IP2)
print(l_IP3)

#Pick a Bandwidth for IP
l_IP=input('Entrez la valeur de l_IP : ')
#%%
print('l_PHI')
print(l_PHI1)
print(l_PHI2)
print(l_PHI3)

#Pick a Bandwidth for PHI
l_PHI=input('Entrez la valeur de l_PHI : ')

#==============================================================================
#
# # Prepare the grid on which to calculate the kernel
#==============================================================================
#%%
min_IP = np.min(IP) - l_IP/2
max_IP = np.max(IP) + l_IP/2


min_PHI = np.min(PHI) - l_PHI
max_PHI = np.max(PHI) + l_PHI

min_IP=7000
max_IP=20000
min_PHI=0.10
max_PHI=0.5

# # Definition of two vectors for kernel resolution and limits
vec_IP=np.linspace(min_IP,max_IP,500)
vec_PHI=np.linspace(min_PHI,max_PHI,500)


grid_x,grid_y = np.meshgrid(vec_IP,vec_PHI)
gril2d = (np.vstack([np.reshape(grid_x,250000),np.reshape(grid_y,250000)])).T

l=len(gril2d)

# KERNEL CALCULATION
r=np.zeros(l)

for q in range(len(IP)):

    g=(1/ l_IP)*(1/ l_PHI)*(1/np.sqrt(2*math.pi))*np.exp(-0.5* \
        (((IP[q]-gril2d[:,0])/ l_IP)**2))*(1/np.sqrt(2*math.pi))*np.exp(-0.5*\
        (((PHI[q]-gril2d[:,1])/l_PHI)**2))

    r=r+g


pdf22=r/len(IP) #On divise par le nombre d'echantillons n (voir KDE def)

## RESHAPE PDF
pdf2d=np.reshape(pdf22,[len(vec_IP),len(vec_PHI)])

# SAVE PDF and variables vectors
with open('stats.pickle', 'w') as f:
    pickle.dump([vec_IP,vec_PHI,pdf2d], f)
#%%
# Start from here to load previously calculated kernel and display it
with open('stats.pickle') as f:
    vec_IP,vec_PHI,pdf2d = pickle.load(f)

fig0,ax0 = plt.subplots(figsize=(5,5))
#ax0.plot(IP,PHI,'o')
ax0.set_xlabel('IP',fontsize=14)
ax0.set_ylabel('Phi',fontsize=14)
#fig0.savefig('relation_scatter.png', bbox_inches='tight',dpi=600)
#
#fig1,ax1 = plt.subplots(figsize=(5,5))
ax0.imshow(np.rot90((pdf2d).T),extent=[7000,20000,0.1,0.5])
#ax0.plot(IP,PHI,'o')
ax0.set_aspect(13000/0.4)
ax0.set_xlim([7000,20000])
ax0.set_ylim([0.1,0.5])
#ax1.set_xlabel('IP',fontsize=14)
#ax1.set_ylabel('Phi',fontsize=14)
fig0.show()
#fig1.savefig('kernel2d.png', bbox_inches='tight',dpi=600)
#%%
pick_IP=(13250,13250)
Phi_line=(0.1,0.5)
fig2,ax2 = plt.subplots(figsize=(5,5))
ax2.imshow(np.rot90((pdf2d).T),extent=[8000,20000,0.1,0.5])
ax2.plot(pick_IP,Phi_line,'k',linewidth=4)
ax2.set_aspect(12000/0.4)
ax2.set_xlabel('IP',fontsize=30)
ax2.set_ylabel('Phi',fontsize=30)
ax2.xaxis.set_ticks([])
ax2.yaxis.set_ticks([])
fig2.show()
fig2.savefig('kernel2d_w_pick_noaxes.png', bbox_inches='tight',dpi=600)

#==============================================================================
#
# # In the case of multiple statistical families in data distribution
## Generate multiple kernels and define a probability of being in each family_probs
# depending on your secondary data (here IP)
#==============================================================================

#%%
def point_inside_polygon(x,y,poly):

    n = len(poly)
    inside =False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if min(p1y,p2y)< y < max(p1y,p2y):
            if x <= max(p1x,p2x):
                if p1y != p2y:
                    xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                if p1x == p2x or x <= xinters:
                    inside = not inside
        p1x,p1y = p2x,p2y

    return inside

#%%


# option ginput et inpolygon pour classer les donnees
############################################################
#Famille 1

print('Family1')
#   pick values in each class
#input('CLASSE '+ str(i))
poly = ginput(20,timeout=0)

poly.append(poly[0])

inside=[]

for xi,yi in zip(IP,PHI):
    inside.append(point_inside_polygon(xi,yi,poly))

inside=np.array(inside)
ax0.plot(IP[inside],PHI[inside],'y*',markersize=15)

IP_family1 = IP[inside]       # IA data in facies
PHI_family1 = PHI[inside]     # Phi data in facies

#%%
############################################################
#Famille 2

print('Family2')
#   pick values in each class
poly = ginput(20,timeout=0)#input('CLASSE '+ str(i))


poly.append(poly[0])

inside=[]
for xi,yi in zip(IP,PHI):
    inside.append(point_inside_polygon(xi,yi,poly))

inside=np.array(inside)
ax0.plot(IP[inside],PHI[inside],'g.',markersize=15)

IP_family2 = IP[inside]       # IA data in facies
PHI_family2 = PHI[inside]     # Phi data in facies
#%%
#####################################################################
#Family1

#Definition des largeurs de bande
## Largeur de bande proposïee par Silverman (1986)

l_IP1 = 0.9 * (min(np.std(IP_family1),iqr_mb(IP_family1)/1.349)) * len(IP_family1)**(-1/5)
l_PHI1 = 0.9 * (min(np.std(PHI_family1),iqr_mb(PHI_family1)/1.349)) * len(PHI_family1)**(-1/5)

## Bowman et Foster
l_IP2=np.std(IP_family1)*len(IP_family1)**(-1/6)
l_PHI2=np.std(PHI_family1)*len(PHI_family1)**(-1/6)


## Regle de pouce (Deheurels (1977))
l_IP3=1.06*np.std(IP_family1)* len(IP_family1)**(-1/5)
l_PHI3=1.06*np.std(PHI_family1)* len(PHI_family1)**(-1/5)

print('l_IP')
print(l_IP1)
print(l_IP2)
print(l_IP3)

#Demande a l'utilisateur la largeur de bande choisie
l_IP_family1=input('Entrez la valeur de l_IP : ')
#%%

print('l_PHI')
print(l_PHI1)
print(l_PHI2)
print(l_PHI3)


#Demande a l'utilisateur la largeur de bande choisie
l_PHI_family1=input('Entrez la valeur de l_PHI : ')

#%%
r=np.zeros(250000)
for q in range(len(IP_family1)):

    g=(1/ l_IP_family1)*(1/ l_PHI_family1)*(1/np.sqrt(2*math.pi))*np.exp(-0.5* \
        (((IP_family1[q]-gril2d[:,0])/ l_IP_family1)**2))*(1/np.sqrt(2*math.pi))*np.exp(-0.5*\
        (((PHI_family1[q]-gril2d[:,1])/l_PHI_family1)**2))

    r=r+g


pdf22_family1=r/len(IP_family1) #On divise par le nombre d'echantillons n (voir KDE def)

## RESHAPE DU PDF 3D EN CUBE
pdf2d_family1=np.reshape(pdf22_family1,[len(vec_IP),len(vec_PHI)])

#%%

#####################################################################
#Family2

#Definition des largeurs de bande
## Largeur de bande proposïee par Silverman (1986)

l_IP1 = 0.9 * (min(np.std(IP_family2),iqr_mb(IP_family2)/1.349)) * len(IP_family2)**(-1/5)
l_PHI1 = 0.9 * (min(np.std(PHI_family2),iqr_mb(PHI_family2)/1.349)) * len(PHI_family2)**(-1/5)

## Bowman et Foster
l_IP2=np.std(IP_family2)*len(IP_family2)**(-1/6)
l_PHI2=np.std(PHI_family2)*len(PHI_family2)**(-1/6)


## Regle de pouce (Deheurels (1977))
l_IP3=1.06*np.std(IP_family2)* len(IP_family2)**(-1/5)
l_PHI3=1.06*np.std(PHI_family2)* len(PHI_family2)**(-1/5)

print('l_IP')
print(l_IP1)
print(l_IP2)
print(l_IP3)

#Demande a l'utilisateur la largeur de bande choisie
l_IP_family2=input('Entrez la valeur de l_IP : ')
#%%

print('l_PHI')
print(l_PHI1)
print(l_PHI2)
print(l_PHI3)

#Demande a l'utilisateur la largeur de bande choisie
l_PHI_family2=input('Entrez la valeur de l_PHI : ')

#%%
r=np.zeros(250000)
for q in range(len(IP_family2)):

    g=(1/ l_IP_family2)*(1/ l_PHI_family2)*(1/np.sqrt(2*math.pi))*np.exp(-0.5* \
        (((IP_family2[q]-gril2d[:,0])/ l_IP_family2)**2))*(1/np.sqrt(2*math.pi))*np.exp(-0.5*\
        (((PHI_family2[q]-gril2d[:,1])/l_PHI_family2)**2))

    r=r+g


pdf22_family2=r/len(IP_family2) #On divise par le nombre d'echantillons n (voir KDE def)

## RESHAPE DU PDF 3D EN CUBE
pdf2d_family2=np.reshape(pdf22_family2,[len(vec_IP),len(vec_PHI)])

#####################################################################################
#%%

sum1 = np.sum(pdf2d_family1,axis=0)/np.sum(pdf2d_family1)
sum2 = np.sum(pdf2d_family2,axis=0)/np.sum(pdf2d_family2)

prob1=sum1/(sum1+sum2)
prob2=1-prob1

fig4,ax4 = plt.subplots(figsize=(5,5))

ax4.plot(vec_IP,prob1,'r-')
ax4.plot(vec_IP,prob2,'b-')
ax4.set_ylabel('Probabilit'u'é',fontsize=14)
ax4.set_xlabel('Imp'u'é''dance acoustique',fontsize=14)
ax4.set_xlim(8000,20000)
ax4.legend({'Famille 1','Famille 2'},fontsize=14)
fig4.savefig('family_probs.png', bbox_inches='tight',dpi=600)

#%%

fig5,ax5 = plt.subplots(figsize=(5,5))
ax5.set_xlabel('AI',fontsize=30)
ax5.set_ylabel(r'$\phi$',fontsize=35)
ax5.imshow(np.rot90((pdf2d).T),extent=[7000,20000,0.1,0.5])
#ax5.plot(IP_family1,PHI_family1,'y*',markersize=12)
#ax5.plot(IP_family2,PHI_family2,'k.',markersize=12)
ax5.set_aspect(13000/0.4)
ax5.set_xlim([7000,20000])
ax5.set_ylim([0.1,0.5])
ax5.set_xticks([])
ax5.set_yticks([])
#ax5.legend({'Family 1','Family 2'},fontsize=20)
fig5.show()

fig5.savefig('kernel.png', bbox_inches='tight',dpi=600)
