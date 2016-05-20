# -*- coding: utf-8 -*-
"""
Created on Mon May  9 13:55:26 2016

@author: mablou
"""
#%%
import numpy as np

_Bayeslib = np.ctypeslib.load_library('libBSS_library.so','.')

#%%
#==============================================================================
#
#           STRUCTURES DEFINITION
#
#==============================================================================
#==============================================================================
# variogram                                            */
# Nvario: number of combined variogram models          */
# vario: model of variogram per variogram model        */
#        1 --> exponential                             */
#        2 --> gaussian                                */
#        3 --> spherical                               */
#        4 --> cardsin                                 */
#        5 --> stable                                  */
#        6 --> gamma                                   */
#        7 --> cubic                                   */
#        8 --> nugget                                  */
#        9 --> power                                   */
# alpha: exponent for the stable and gamma variogram   */
#        per variogram model                           */
# ap: anisotropy axes per variogram model              */
# scf: correlation lengths per variogram model         */
# var: normalized variance per variogram model(sum = 1)*/
# 
#==============================================================================
class vario_mod(np.ctypeslib.ctypes.Structure):
    _fields_ = [
                ('Nvario',np.ctypeslib.ctypes.c_int),
                ('vario',np.ctypeslib.ndpointer(int)),
                ('alpha',np.ctypeslib.ndpointer(float)),
                ('ap',np.ctypeslib.ctypes.ndpointer(float)),
                ('scf',np.ctypeslib.ctypes.ndpointer(float)),
                ('var',np.ctypeslib.ctypes.ndpointer(float))
                ]
                
#==============================================================================
# grid_seismic								
# x: x value at the center of cell			
# y: y value at the center of cell			
# z: z value at the center of cell			
# AI: AI value at the center of cell	   	
# N: N number of points in the grid			
# 
#==============================================================================

class grid_seismic(np.ctypeslib.ctypes.Structure):
    _fields_ = [
                ('x',np.ctypeslib.ndpointer(float)),
                ('y',np.ctypeslib.ndpointer(float)),
                ('z',np.ctypeslib.ctypes.ndpointer(float)),
                ('AI',np.ctypeslib.ctypes.ndpointer(float)),
                ('IS',np.ctypeslib.ctypes.ndpointer(float)),
                ('N', np.ctypeslib.ctypes.c_long)
                ]                
#==============================================================================
# prob_mod															
# fam_prob	: array containing the probability of being in family 1	
# pdf2d		: Joint probability distribution of var1/var2		
# vec_var1	: useful Vector of shape min_var1...dvar1...maxvar1	
# n		: length of vector var1 and var2 must be the same	
# vec_var2	: useful Vector of shape min_var2...dvar2...maxvar2	
# mean		: double array containing mean for each family		
# std		: standard deviation for each family				
# nb_family : number of family in kernel or # of facies 		
#==============================================================================
class prob_mod(np.ctypeslib.ctypes.Structure):
    _fields_ = [
                ('fam_prob',np.ctypeslib.ndpointer(float)),
                ('pdf2d',np.ctypeslib.ndpointer(float)),
                ('pdf3d',np.ctypeslib.ctypes.ndpointer(float)),
                ('pdfupscale',np.ctypeslib.ctypes.ndpointer(float)),
                ('pdfupscale2',np.ctypeslib.ctypes.ndpointer(float)),
                ('vec_var1',np.ctypeslib.ctypes.ndpointer(float)),
                ('vec_var2',np.ctypeslib.ctypes.ndpointer(float)),
                ('vec_var3',np.ctypeslib.ctypes.ndpointer(float)),
                ('mean',np.ctypeslib.ctypes.ndpointer(float)),
                ('std',np.ctypeslib.ctypes.ndpointer(float)),
                ('n', np.ctypeslib.ctypes.c_long),
                ('nb_family', np.ctypeslib.ctypes.c_int)
                ]

#==============================================================================
# /*well data                                               */
# /*nwell: number of wells                                  */
# /*n: number of measurement points per well                */
# /*   i = [0...nwell-1]                                    */
# /*ntype: number of measurement types                      */
# /*code: status of the measurements i=[0...ntype-1]        */
# /*      --> 0 : Gaussian white noise                      */
# /*      --> 1: standard Normal                            */
# /*      --> 2: non standard Normal                        */
# /*      --> 3: lognormal (neperien)                       */
# /*      --> 4: lognormal (log10)                          */
# /*      --> 5: facies                                     */
# /*      --> 6: uniform                                    */
# /*      --> 7: any                                        */
# /*x: X-coordinates of the measurements                    */
# /*   i = [0 ... n[0]-1 n[0] ... n[0]+n[1]-1...sum(n[k])-1]*/
# /*y: Y-coordinates of the measurements                    */
# /*   i = [0 ... n[0]-1 n[0] ... n[0]+n[1]-1...sum(n[k])-1]*/
# /*z: Z-coordinates of the measurements                    */
# /*   i = [0 ... n[0]-1 n[0] ... n[0]+n[1]-1...sum(n[k])-1]*/
# /*var1: values of the measurements                     */
# /*   same kind of indexation, but repeated per type of    */
# /*   measurement                                          */
# /*   type 1 :                                             */
# /*   i = [0 ... n[0]-1 n[0] ... n[0]+n[1]-1...sum(n[k])-1]*/
# /*   type 2 :                                             */
# /*   i=[sum(n[k])... sum(n[k])+n[0]-1 ... 2*(sum(n[k])-1)]*/
#==============================================================================
class welldata_mod(np.ctypeslib.ctypes.Structure):
    _fields_ = [
                ('nwell',np.ctypeslib.ctypes.c_int),
                ('n',np.ctypeslib.ctypes.c_int),
                ('ntype',np.ctypeslib.ctypes.c_int),
                ('code',np.ctypeslib.ctypes.ndpointer(int)),
                ('x',np.ctypeslib.ctypes.ndpointer(float)),
                ('y',np.ctypeslib.ctypes.ndpointer(float)),
                ('z',np.ctypeslib.ctypes.ndpointer(float))
                ('var1',np.ctypeslib.ctypes.ndpointer(float))
                ] 

#==============================================================================
# realization_mod 
# 																
# result: 2D array with 3 first column XYZ + nb_simulations columns				
# n	: n_rows correspond to number of points in the grid + number of well data 	
# rn	: array containing the random sequence of nodes visited					
# mode	: 1: normal 2-->downscale | 3---> 3D pdf |4 ---> both									
#==============================================================================
class realization_mod(np.ctypeslib.ctypes.Structure):
    _fields_ = [
                ('mode',np.ctypeslib.ctypes.c_int),
                ('n',np.ctypeslib.ctypes.c_long),
                ('result',np.ctypeslib.ctypes.ndpointer(float)),
                ('rn',np.ctypeslib.ctypes.ndpointer(float)),
                ] 


#==============================================================================
# pair structure is used to return two values from minMax()
#==============================================================================
class pair(np.ctypeslib.ctypes.Structure):
    _fields_ = [
                ('min',np.ctypeslib.ctypes.c_double),
                ('max',np.ctypeslib.ctypes.c_double)    
                ]
                

 #%%
#==============================================================================
# 
#               FUNCTION ARGUMENTS DEFINITION
#              
#==============================================================================
               

_Bayeslib.BayesianSimulation.argtypes = (welldata_mod,vario_mod,
                                         prob_mod, grid_seismic, 
                                         np.ctypeslib.ctypes.c_int,
                                         np.ctypeslib.ctypes.c_int,
                                         np.ctypeslib.ctypes.POINTER(realization_mod))


                
_Bayeslib.getMinMax.argtypes = (np.ctypeslib.ctypes.POINTER(pair),
                                np.ctypeslib.ndpointer(float),
                                np.ctypeslib.ctypes.c_int)
                                
#%%
#==============================================================================
# 
#               PYTHON FUNCTIONS
#                                
#==============================================================================

def BayesianSimulation(points,grid,variogram,stat_struct,
                           nb_simul=1,mode='normal'):
    
    points = np.require(points,float)
    grid = np.require(grid,float)
    n_tot = len(points) + len(grid)    
    
    realization_out = realization_mod()
    realization_out.mode = mode
    realization_out.n = n_tot
    realization_out.result = np.zeros((n_tot,3+nb_simul))    
    
    hard_data = welldata_mod()
    hard_data.x = points[:,0]
    hard_data.y = points[:,1]
    hard_data.z = points[:,2]
    hard_data.var1 = points[:,3]
    hard_data.n = len(points)
    
    _vario = vario_mod()
    
    
    
    _Bayeslib.BayesianSimulation()
    


def getMinMax(in_array):
    in_array = np.require(in_array,float)
    out_pair = pair()
    _Bayeslib.getMinMax(out_pair,in_array,in_array.size)
    return out_pair        

















                        