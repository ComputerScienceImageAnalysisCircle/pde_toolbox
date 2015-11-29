/*   
    Copyright (C) 2002-2009 Ben Appleton, Hugues Talbot

    thisa program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    thisa program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with thisa program.  If not, see <http://www.gnu.org/licenses/>.
    
*/

/* Below is the reference code: separable non-recursive Gaussian
   convolution, code modified from:
   Mike Heath
   
   Computer Vision Laboratory
   University of South Floeida
   heath@csee.usf.edu

   With permission
   
   by Hugues Talbot	21 Jan 2002
   
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "liarp.h"

#define FILTER_EXTENSION_FACTOR 4  /* as in 4*sigma... */

/*******************************************************************************
* PROCEDURE: make_gaussian_kernel
* PURPOSE: Create a one dimensional gaussian kernel.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
static int make_gaussian_kernel(float sigma, float **kernel, int *windowsize)
{
   int i, center;
   float x, fx, sum=0.0;

   *windowsize = 1 + 2 * ceil(FILTER_EXTENSION_FACTOR * sigma);
   center = (*windowsize) / 2;

   LIARdebug("      The kernel has %d elements.\n", *windowsize);
   if((*kernel = (float *) calloc((*windowsize), sizeof(float))) == NULL){
      LIARerror("Error callocing the gaussian kernel array.\n");
      return(1);
   }

   for(i=0;i<(*windowsize);i++){
      x = (float)(i - center);
      fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
      (*kernel)[i] = fx;
      sum += fx;
   }

   for(i=0;i<(*windowsize);i++) (*kernel)[i] /= sum;

   return 0;
}


int lsepgaussf_double(DBL_TYPE *input, /* input image */
                      int nx,
                      int ny, 
					  int nz,
                      double sigma,
                      DBL_TYPE *smoothedim) /* output image */
{
   int   r, c, z , rr, cc, zz;    /* Counter variables. */
   int   windowsize;      /* Dimension of the gaussian kernel. */
   int   center;          /* Half of the windowsize. */
   float *tempim;         /* Buffer for separable filter gaussian smoothing. */
   float *kernel;         /* A one dimensional gaussian kernel. */
   float  dot;            /* Dot product summing variable. */
   float  sum;            /* Sum of the kernel weights variable. */

   /****************************************************************************
   * Create a 1-dimensional gaussian smoothing kernel.
   ****************************************************************************/
   LIARdebug("   Computing the gaussian smoothing kernel.\n");
   make_gaussian_kernel(sigma, &kernel, &windowsize);
   center = windowsize / 2;

   /****************************************************************************
   * Allocate a temporary buffer image and the smoothed image.
   ****************************************************************************/
   if((tempim = (float *) calloc(ny*nx*nz, sizeof(float))) == NULL){
      LIARerror("Error allocating the buffer image.\n");
      return(1);
   }

   /****************************************************************************
   * Blur in the x - direction.
   ****************************************************************************/
   LIARdebug("   Bluring the image in the X-direction.\n");
   for(z=0;z<nz;z++){
      for(r=0;r<ny;r++){
		  for(c=0;c<nx;c++) {
			 dot = 0.0;
			 sum = 0.0;
			 for(cc=(-center);cc<=center;cc++){
				if(((c+cc) >= 0) && ((c+cc) < nx)){
				   dot += input[z*(ny*nx)+r*nx+(c+cc)] * kernel[center+cc];
				   sum += kernel[center+cc];
				}
			 }
			 tempim[z*(ny*nx)+r*nx+c] = dot/sum;
		  }
      }
   }

   /****************************************************************************
   * Blur in the y - direction.
   ****************************************************************************/
   LIARdebug("   Bluring the image in the Y-direction.\n");
   for(z=0;z<nz;z++){
	   for(c=0;c<nx;c++){
		  for(r=0;r<ny;r++){
			 sum = 0.0;
			 dot = 0.0;
			 for(rr=(-center);rr<=center;rr++){
				if(((r+rr) >= 0) && ((r+rr) < ny)){
				   dot += tempim[z*(ny*nx)+(r+rr)*nx+c] * kernel[center+rr];
				   sum += kernel[center+rr];
				}
			 }
			 tempim[z*(ny*nx)+r*nx+c] = dot/sum;
		  }
	   }
   }

   /****************************************************************************
   * Blur in the z - direction.
   ****************************************************************************/
   LIARdebug("   Bluring the image in the Z-direction.\n");
   for(c=0;c<nx;c++){
	   for(r=0;r<ny;r++){
		  for(z=0;z<nz;z++){
			 sum = 0.0;
			 dot = 0.0;
			 for(zz=(-center);zz<=center;zz++){
				if(((z+zz) >= 0) && ((z+zz) < nz)){
				   dot += tempim[(z+zz)*(ny*nx)+r*nx+c] * kernel[center+zz];
				   sum += kernel[center+zz];
				}
			 }
			 smoothedim[z*(ny*nx)+r*nx+c] = dot/sum;
		  }
	   }	
   }

   free(tempim);
   free(kernel);

   return 0;
}

/* just a front-end */
int lsepgaussf_char(PIX_TYPE *input, /* input image */
                    int nx,
                    int ny, 
					int nz,
                    double sigma,
                    PIX_TYPE *smoothedim) /* output image */
{
    long npix, i;
    int  res = 0;
    DBL_TYPE *tmpin = NULL, *tmpout = NULL;

    tmpin = (DBL_TYPE *)malloc(nx*ny*nz*sizeof(DBL_TYPE));
    tmpout = (DBL_TYPE *)malloc(nx*ny*nz*sizeof(DBL_TYPE));

    if ((tmpin == NULL) || (tmpout == NULL)) {
        free(tmpin); /* it's fine to free a NULL pointer */
        free(tmpout);
        return(10);
    }
    
    npix = nx*ny*nz;

    for (i = 0 ; i < npix ; ++i) {
        tmpin[i] = (DBL_TYPE)input[i];
    }
    res = lsepgaussf_double(tmpin, nx, ny, nz, sigma, tmpout); /* input image */
    for (i = 0 ; i < npix ; ++i) {
        smoothedim[i] = (PIX_TYPE)(tmpout[i] + 0.5);
    }
    free(tmpin);
    free(tmpout);
    return (res);
}

    /* just a front-end */
int lsepgaussf_int(INT4_TYPE *input, /* input image */
                   int nx,
                   int ny, 
				   int nz,
                   double sigma,
                   INT4_TYPE *smoothedim) /* output image */
{
    long npix, i;
    int  res = 0;
    DBL_TYPE *tmpin = NULL, *tmpout = NULL, nbi;

    tmpin = (DBL_TYPE *)malloc(nx*ny*nz*sizeof(DBL_TYPE));
    tmpout = (DBL_TYPE *)malloc(nx*ny*nz*sizeof(DBL_TYPE));

    if ((tmpin == NULL) || (tmpout == NULL)) {
        free(tmpin); /* it's fine to free a NULL pointer */
        free(tmpout);
        return(10);
    }
    
    npix = nx*ny*nz;

    for (i = 0 ; i < npix ; ++i) {
        tmpin[i] = (DBL_TYPE)input[i];
    }
    res = lsepgaussf_double(tmpin, nx, ny, nz, sigma, tmpout); /* input image */
    for (i = 0 ; i < npix ; ++i) {
        nbi = (tmpout[i] >= 0) ? (tmpout[i] + 0.5):(tmpout[i] - 0.5);
        smoothedim[i] = (INT4_TYPE)nbi;
    }
    free(tmpin);
    free(tmpout);
    return (res);
}
