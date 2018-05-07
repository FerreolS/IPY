/*
 * svd2D.c --
 *
 * Implements basic "svd2D" as in global bioimlib

  Let X be a NxMx3 matrix such that:
  
  P_mn = [X(n,m,1) X(n,m,2)
          X(n,m,2) X(n,m,3)]
          
  is a symmetric matrix. Then the present function computes the eigenvalues
  E(n,m,1) E(n,m,2) and the first eigenvector [V(n,m,1) V(n,m,2)] (the second 
  one being [V(n,m,2) -V(n,m,1)]). Hence the function outputs two matrices E
  and V of size NxMx2.
 *
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2018: Ferréol Soulez
 *
 * This file is part of free software IPY: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * IPY is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *-----------------------------------------------------------------------------
 */
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>


#include <pstdlib.h>
#include <play.h>
#include <yapi.h>

/* Define some macros to get rid of some GNU extensions when not compiling
   with GCC. */
#if ! (defined(__GNUC__) && __GNUC__ > 1)
#   define __attribute__(x)
#   define __inline__
#   define __FUNCTION__        ""
#   define __PRETTY_FUNCTION__ ""
#endif

#define TRUE  1
#define FALSE 0

PLUG_API void y_error(const char *) __attribute__ ((noreturn));

static int same_dims(const long xdims[], const long ydims[]);



/*---------------------------------------------------------------------------*/
/* BUILT-IN FUNCTIONS */
/* E=svd2D_decomp(V,X) */
extern void Y_svd2D_decomp(int argc);

void
Y_svd2D_decomp(int argc)
{
  double *x;
  long dims[Y_DIMSIZE];
  long tdims[Y_DIMSIZE];
  long dimsout[Y_DIMSIZE];
  double *Yv,*Ye ;
  long ntot, i, k, num_of_mat;
  int rtype;
  int fresh;

	

  double n;
  double tmp[3];
  double E[2];
  double U[2];    
  double trace;
  double delta;
        
  
  if (argc != 2) y_error("expecting exactly two argument");
 
  x = ygeta_d(0, &ntot, dims);
  

  if (dims[0]!=3)  y_error("dims(X) should be 3");

  if (dims[3]!=3)  y_error("The last dimension of the input should be equal to 3.");

  num_of_mat = dims[1]*dims[2];   // number lateral entries (i.e. number of svd to compute)


  dimsout[0] = 3;
  dimsout[1] =  dims[1];
  dimsout[2] =  dims[2];
  dimsout[3] =  2;

  
  Yv = NULL;
  Ye = NULL;
  if (argc == 2) {
    /* Get current value of destination. */
    if (yarg_typeid(1) == Y_DOUBLE) {
      Yv = ygeta_any(1, NULL, tdims, &rtype);
      if (!(rtype == Y_DOUBLE && same_dims(tdims, dims))) {
        Yv = NULL;
      }
    }
    if (Yv == NULL) {
      /* Discard destination value because it cannot be reused. */
      ypush_nil();
      yarg_swap(argc, 0);
      yarg_drop(1);
    }
  }
  
  /* Create the destination if needed. */
  fresh = (Yv == NULL);
  if (fresh) {
      Yv = ypush_d(dimsout);
      Ye = ypush_d(dimsout);
    
  }

   #pragma omp parallel for shared(x,Ye,Yv) private(i, k, trace, delta, n, tmp, E, U)
    for(i=0; i < num_of_mat; i++){
    	for (k=0;k<3;k++)   // get the matrix value [X(1,1) X(2,1)=X(1,2), X(2,2)]
        	tmp[k]=x[i+num_of_mat*k];
        	
        if (fabs(tmp[1]) < 1e-15){
    		E[0]=tmp[0];
    		E[1]=tmp[2];
    		U[0]=1.0;
    		U[1]=0.0;
    	}
    	else{
    		trace=tmp[0]+tmp[2];
    		delta=(tmp[0]-tmp[2])*(tmp[0]-tmp[2])+4*tmp[1]*tmp[1];
    		E[0]=0.5*(trace+sqrt(delta));
    		E[1]=0.5*(trace-sqrt(delta));
    		n=sqrt((E[0]-tmp[0])*(E[0]-tmp[0])+tmp[1]*tmp[1]);
    		U[0]=tmp[1]/n;
    		U[1]=(E[0]-tmp[0])/n;
  		}
  		
  		for (k=0;k<2;k++){  // set result
        	Ye[i+num_of_mat*k]=E[k];
  			Yv[i+num_of_mat*k]=U[k];
  		}
    }
}



/*---------------------------------------------------------------------------*/
/* BUILT-IN FUNCTIONS */
/*[X]=svd2D_recomp(E,V) */
extern void Y_svd2D_recomp(int argc);

void
Y_svd2D_recomp(int argc)
{
  double *x;
  long dims[Y_DIMSIZE];
  long tdims[Y_DIMSIZE];
  long dimsout[Y_DIMSIZE];
  double *Yv,*Ye ;
  long ntot, i, k, num_of_mat;
  
  
  double ee[2];
  double vv[2];
  double tmp[3];  
	
      
  
  if (argc != 2) y_error("expecting exactly two argument");
 
  Ye = ygeta_d(1, &ntot, dims);
  

  if (dims[0]!=3)  y_error("dims(X) should be 3");
  
  if (dims[3]!=2)  y_error("The last dimension of the input should be equal to 2.");
  
  num_of_mat = dims[1]*dims[2];   // number lateral entries (i.e. number of svd to compute)
  
  
  if (yarg_typeid(0) == Y_DOUBLE) {
    Yv = ygeta_d(0, NULL, tdims );
    if (!same_dims(tdims, dims)) {
      y_error("E and V should have the same dimension.");
    }
  }
  
  
  dimsout[0] = 3;
  dimsout[1] =  dims[1];
  dimsout[2] =  dims[2];
  dimsout[3] =  3;
  
  
  x  = ypush_d(dimsout);
  
  
  
#pragma omp parallel for shared(Ye,Yv,x) private(i, k, ee,vv,tmp)
  for(i=0; i < num_of_mat; i++){
    for (k=0;k<2;k++){   // get the eigenvalues and eigenvectors
      ee[k]=Ye[i+num_of_mat*k];
      vv[k]=Yv[i+num_of_mat*k];
    }
    tmp[0]=ee[0]*pow(vv[0],2) + ee[1]*pow(vv[1],2);
    tmp[1]=vv[0]*vv[1]*(ee[0]-ee[1]);
    tmp[2]=ee[0]*pow(vv[1],2)+ee[1]*pow(vv[0],2);
    for (k=0;k<3;k++){  // set result
      x[i+num_of_mat*k]=tmp[k];
    }
  }
  
}   	




/*---------------------------------------------------------------------------*/
/* UTILITIES */

static int
same_dims(const long xdims[], const long ydims[])
{
  int i, xrank, yrank;
  if (xdims == ydims) return TRUE;
  xrank = (xdims == NULL ? 0 : xdims[0]);
  yrank = (ydims == NULL ? 0 : ydims[0]);
  if (xrank != yrank) return FALSE;
  for (i = 1; i <= xrank; ++i) {
    if (xdims[i] != ydims[i]) return FALSE;
  }
  return TRUE;
}
