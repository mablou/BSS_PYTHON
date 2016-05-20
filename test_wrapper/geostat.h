#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>


#ifndef _GEOSTAT_H
#define _GEOSTAT_H



#define NTAB 32

/*STRUCTURES*/
/*----------*/
/*variogram                                            */
/*Nvario: number of combined variogram models          */
/*vario: model of variogram per variogram model        */
/*       1 --> exponential                             */
/*       2 --> gaussian                                */
/*       3 --> spherical                               */
/*       4 --> cardsin                                 */
/*       5 --> stable                                  */
/*       6 --> gamma                                   */
/*       7 --> cubic                                   */
/*       8 --> nugget                                  */
/*       9 --> power                                   */
/*alpha: exponent for the stable and gamma variogram   */
/*       per variogram model                           */
/*ap: anisotropy axes per variogram model              */
/*scf: correlation lengths per variogram model         */
/*var: normalized variance per variogram model(sum = 1)*/
struct vario_mod {
  int Nvario;
  int *vario;
  double *alpha;
  double *ap;
  double *scf;
  double *var;
};

/*grid_seismic								*/
/*x: x value at the center of cell			*/
/*y: y value at the center of cell			*/
/*z: z value at the center of cell			*/
/*AI: AI value at the center of cell		*/
/*N: N number of points in the grid			*/
struct grid_seismic  {
	double *x;
	double *y;
	double *z;
	double *AI;
	double *IS;
	long int N;
};

/*prob_mod															*/
/*fam_prob	: array containing the probability of being in family 1	*/
/*pdf2d		: Joint probability distribution of var1/var2		*/
/*vec_var1	: useful Vector of shape min_var1...dvar1...maxvar1	*/
/*n		: length of vector var1 and var2 must be the same	*/
/*vec_var2	: useful Vector of shape min_var2...dvar2...maxvar2	*/
/*mean		: double array containing mean for each family		*/
/*std		: standard deviation for each family				*/
/*nb_family : number of family in kernel or # of facies 		*/

struct prob_mod  {
	double *fam_prob;
	double *pdf2d;
	double *pdf3d;
	double *pdfupscale;
	double *pdfupscale2;
	double  *vec_var1;
	double  *vec_var2;
	double *vec_var3;
	double *mean;
	double *std;
	long int n;
	int nb_family;
};

/*grid                                     */
/*NX: number of gridblocks along the X axis*/
/*NY: number of gridblocks along the Y axis*/
/*NZ: number of gridblocks along the Z axis*/
/*DX: gridblock length along the X axis    */
/*DY: gridblock length along the Y axis    */
/*DZ: gridblock length along the Z axis    */
/*Xo: X-cell number of the origin cell     */
/*Yo: Y-cell number of the origin cell     */
/*Zo: Z-cell number of the origin cell     */
struct grid_mod {
  int NX, NY, NZ;
  double DX,DY,DZ;
  double Xo,Yo,Zo;
};


/*well data                                               */
/*nwell: number of wells                                  */
/*n: number of measurement points per well                */
/*   i = [0...nwell-1]                                    */
/*ntype: number of measurement types                      */
/*code: status of the measurements i=[0...ntype-1]        */
/*      --> 0 : Gaussian white noise                      */
/*      --> 1: standard Normal                            */
/*      --> 2: non standard Normal                        */
/*      --> 3: lognormal (neperien)                       */
/*      --> 4: lognormal (log10)                          */
/*      --> 5: facies                                     */
/*      --> 6: uniform                                    */
/*      --> 7: any                                        */
/*x: X-coordinates of the measurements                    */
/*   i = [0 ... n[0]-1 n[0] ... n[0]+n[1]-1...sum(n[k])-1]*/
/*y: Y-coordinates of the measurements                    */
/*   i = [0 ... n[0]-1 n[0] ... n[0]+n[1]-1...sum(n[k])-1]*/
/*z: Z-coordinates of the measurements                    */
/*   i = [0 ... n[0]-1 n[0] ... n[0]+n[1]-1...sum(n[k])-1]*/
/*var1: values of the measurements                     */
/*   same kind of indexation, but repeated per type of    */
/*   measurement                                          */
/*   type 1 :                                             */
/*   i = [0 ... n[0]-1 n[0] ... n[0]+n[1]-1...sum(n[k])-1]*/
/*   type 2 :                                             */
/*   i=[sum(n[k])... sum(n[k])+n[0]-1 ... 2*(sum(n[k])-1)]*/
struct welldata_mod {
  int nwell;
  int n;
  int ntype;
  int *code;
  double *x;
  double *y;
  double *z;
  double *var1;
};



/*kriging mod                                        */
/*n: number of components                             */

struct kriging_mod {
  double n;
  double *vector;
  double variance;
};



/*variogram table                                      */
/*Nvario: number of combined variogram models          */
/*vario: model of variogram per variogram model        */
/*       1 --> exponential                             */
/*       2 --> gaussian                                */
/*       3 --> spherical                               */
/*       4 --> cardsin                                 */
/*       5 --> stable                                  */
/*       6 --> gamma                                   */
/*       7 --> cubic                                   */
/*       8 --> nugget                                  */
/*       9 --> power                                   */
/*alpha: exponent for the stable and gamma variogram   */
/*       per variogram model                           */
/*ap: anisotropy axes per variogram model              */
/*scf: correlation lengths per variogram model         */
/*var: normalized variance per variogram model(sum = 1)*/
struct variotable_mod {
  int number_of_variograms;
  int *Nvario;
  int *vario;
  double *alpha;
  double *ap;
  double *scf;
  double *var;
};


/*realization_mod 																*/
/*result: 2D array with 3 first column XYZ + nb_simulations columns				*/
/*n	: n_rows correspond to number of points in the grid + number of well data 	*/
/*rn	: array containing the random sequence of nodes visited					*/
/*mode	: 1: normal 2-->upscale | 3---> 3D pdf |4 ---> both									*/
struct realization_mod  {
	double *result;
	double *rn;
	int mode;
	long int n;
};

/* structure is used to return two values from minMax() */
struct pair
{
  double min;
  double max;
}; 

/*points   									*/
/*n: number of points						*/
/*x: X-coordinates of the measurements		*/
/*y: Y-coordinates of the measurements		*/
/*z: Z-coordinates of the measurements    	*/
/*var1: values of the measurements			*/

struct points {
  long int n;
  double *x;
  double *y;
  double *z;
  double *var1;
};
 
/*=====================================================*/

/*FUNCTIONS*/
/*---------*/

/*normalization of the anostropy axes                  */
/*ap: anisotropy axes                                  */
/*scf: correlation lengths                             */
/* The returned normalized axes are in ap              */
void axes(double *ap, double *scf, int N);

/*cardsin covariance value for lag h*/
double cardsin(double h);

/*calculation of the covariance value for a distance h */
/*defined by i,j,k                                     */
/*available variogram model:                           */
/* 1 -> exponential                                    */
/* 2 -> gaussian                                       */
/* 3 -> spherical                                      */
/* 4 -> cardsin                                        */
/* 5 -> stable                                         */
/* 6 -> gamma                                          */
/* 7 -> cubic                                          */
/* 8 -> nugget                                         */
/* 9 -> power                                          */
/*variogram: variogram with the structure defined above*/
/*di: distance along the X axis                        */
/*dj: distance along the Y axis                        */
/*dk: distance along the Z axis                        */
/* The returned value is the computed covariance value */
double cov_value(struct vario_mod variogram,double di,double dj,double dk);

/*cubic covariance value for lag h*/
double cubic(double h);

/*exponential covariance value for lag h*/
double exponential(double h);



/*gamma covariance value for lag h and exponent alpha*/
double gammafu(double h,double alpha);

/*gaussian covariance value for lag h*/
double gaussian(double h);

/*computation of the covariance matrix for the well data*/
/*well coordinates have no specific unit                */
/*The dimension of C is given by n                      */
/*C defines the correlation between the first n wells   */
/*C is recorded as a vector so that                     */
/*C[k] = Cij with i = [0...nwell-1], j = [0...nwell-1]  */
/*and k = j+i(i+1)/2                                    */
/*variogram: structure defined above                    */
/*well: structure defined above                         */
void gen_cov_matrix(double *C, struct vario_mod variogram, struct points well, int n);

/*calculates C.x*/
/* C     : symmetric positive-definite matrix recorded */
/*         (per raws) as a vector with only components */
/*         Cij so that j <= i, i = [0...n-1]           */
/* x     : vector, xi avec i = [0...n-1]               */
/* b     : vector, bi avec i = [0...n-1]               */
/* n     : dimension of matrix Cij                     */
/*                                                     */
/* The result vector is returned in b                  */
void mat_vec(double *C, double *x, double *b, int n);

/*2-norm of vector b                 */
/* b : vector                        */
/* n : length of b, bi, i = [0...n-1]*/
/*returns the norm of b              */
double norm(double *b,int n);

/*nugget covariance value for lag h*/
double nugget(double h);

/*power covariance value for lag h and exponent alpha*/
double power(double h,double alpha);

/*calculates bt.b                */
/* b : vector, bi, i = [0...n-1] */
/* n : length of b               */
/*returns the scalar product of b*/
double scal_vec(double *b,int n);


/*kriging                                            */
/*interpolation from ponctual data with simple       */
/*kriging. kriging is used to interpolate standard   */
/*Gaussian data only                                 */
/*inversion method = CONJUGATE GRADIENTS             */
/*INPUT                                              */
/*well: structure with the well data                 */
/*variogram: structure describing the variogram model*/
/*OUTPUT                                             */
/*realout: structure defining the output realization */

void SimpleKriging(struct points data,struct vario_mod variogram,double ptToKrig[3], struct kriging_mod *realout);



/*solves the set of n linear equations Cx = D           */
/* C     : symmetric positive-definite matrix recorded  */
/*         (per raws) as a vector with only components  */
/*         Cij so that j <= i, i = [0...n-1]            */
/* D     : right-hand side vector, Di avec i = [0...n-1]*/
/* n     : dimension of matrix Cij                      */
/*                                                      */
/* The solution vector is returned in D                 */
/* CONJUGATE GRADIENT method                            */
void solve3(double *C, double *D, int n);

/*spherical covariance value for lag h*/
double spherical(double h);

/*stable covariance value for lag h and exponent alpha*/
double stable(double h, double alpha);

/* tb1.b2                                    */
/* b1 : vector                               */
/* b2 : vector                               */
/* n : length of b1 and b2, bi, i = [0...n-1]*/
/*returns the norm the product tb1.b2        */
double vec_vec(double *b1, double *b2,int n);

/*======================================================================================*/
/*	functions added for Baysian simulations						*/
/*======================================================================================*/

/*Function for Bayesian Simulation		*/
/*The algorithm uses a fixed random path	*/
/*Use of this function requires preliminary	*/
/*variogram construction and KDE estimation	*/
/*
 *INPUT:					
 * data: structure with initial conditionning well data 
 * variogram1: structure containing the variogram model parameters for family 1
 * variogram2: structure containing the variogram model parameters for family 2
 * stats: structure containing stats for var1 for each family and the link with var2
 * grid: structure containing the simulation grid populated with var2
 * nb_neighbours: integer designing the number of points used in the kriging routine
 * nb_simul: integer designing the number of stochastic simulations to realize

 *OUTPUT:
 * (*bayesout): structure containing all the realizations for all the points in the grid
 */

void BayesianSimulation (struct welldata_mod data, struct vario_mod variogram, struct prob_mod stats, struct grid_seismic grid,int nb_neighbours,int nb_simul, struct realization_mod *bayesout);

/*function to construct a cumulative density function (cdf)*/
/*uses array of normalized posteriori distribution*/
void buildcdf(double *post,double *cdf,int n);

/*function to calculate the euclidian distance between to points x,y,z*/
double calc_radius(double pt1[3],double pt2[3]);

/*calculates D.x*/
/* D     : a cubic matrix of n*n*n				*/
/* x     : vector, xi avec i = [0...n-1]*/
/* b     : vector, bi avec i = [0...n-1]*/
/* n     : dimension of matrix Cij*/
/* */
/* The solution vector is returned in b*/

void cubic_mat_vec(double *D, double *x, double *b, int n);

/*function thats does the division element by element of 2 vectors*/
/* array 1 = array1 ./ array2					  */
void divide_2vec(double *array1, double *array2, int n);

/*function used to find nearest neighbour*/
/*function called by interp1*/
int findNearestNeighbourIndex( double value, double *x, int len );

/*calculates C.x*/
/* C     : square matrix				*/
/* x     : vector, xi avec i = [0...n-1]*/
/* b     : vector, bi avec i = [0...n-1]*/
/* n     : dimension of matrix Cij*/
/* */
/* The solution vector is returned in b*/

void fullmat_vec(double *C, double *x, double *b, int n);

/*function that return min and max of an array*/
void getMinMax(struct pair *minmax,double *arr, int n);

/*function to determine points used for kriging*/
/*according to the nb_neighbours input by the user*/
/*INPUT:																			*/
/* nb_neighbours: integer representing the number of points max used for kriging 	*/
/* ptToKrig	: coordinates and value of var1 to calculate distance					*/
/* temp_dataset	: conditionning points for the simulation at this step				*/

/*OUTPUT:																			*/
/* idx_list: indexes of points within neighbourhood									*/
void getNeighbourhood (int nb_neighbours,double ptToKrig[3],double *data,long int n,long int *idx_list,long int N);

/*Function to build an array in 1D with regular steps
 *Takes min, spacing and number of values as inputs and the pointer to the returned
 *array					*/
void grid1D (double min, double dx, int nx, double *xx);


/*function used for linear interpolation in 1D*/
/*similar usage as interp1 in Matlab*/
void interp1(double *x, int nx, double *y, double *xx, int nxx, double *yy);

/*function to calculate the marginal probability*/
/*Taking a 2d or a 3d array as input and doing the cumulative sum column by column*/
void marg_calc(double *pdf,int n,double *marginal,int dim);

/*calculates D(TRANSPOSED)*D*/
/* D     : inverse of C output Cholesky decomposed matrix (inv(L))*/
/* n     : dimension of matrix Dij*/
/* */
/* The solution vector is returned in *D */

void matT_mat(double *D, int n);

/*function to calculate the mean of values in an array
 * inputs required are an array of doubles and the number of elements
 * in the array								*/
double mean_function(double array[],int n);

/*function thats does the multiplication element by element of 2 vectors*/
/* array3= array1 .* array2											  */
void multi_2vec(double *array1, double *array2, int n,double *array3);

/*function that returns the sum of an array*/
/*Takes the array and its size n as inputs*/
double sum_function(double *array,int n);

/*function to divide the elements of an array by the sum of its element*/
void normalize(double *array, int n);

/*function to calculate a normal distribution*/
/*inputs required are mean, std dev and a point x*/
void normpdf (double *x, int nx, double mean, double stddev, double *gauss );

/*Perform cholesky decomposition of C			        */
/*And then solve Lx = D									*/

void solve_choldc(double *C, double *D, int n);

/*program to sort x,y,z,var1 values according to radius*/
void sort_ascend(double *radius,long int *idx_list,int n,int nb_neighbours);

/*function thats substract a constant from elements of an array*/
/* array 1 = array1 - cst				  */
void sub_2vec(double *array1, double mean, int n);



#endif
