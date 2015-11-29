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

/*
 * File:		lbimage_utils.c
 *
 * Written by:		Image Analysis Group staff,
 * 			CSIRO Mathematical and Information Sciences.
 *
 * Date:		November 2002
 *
 *
*/

/*********************************************************************************************
 lbimage_utils.c
 ------

  DESCRIPTION:
  A VoiR-compatible wrapper for the utilities in bimage_utils.c

  HISTORY:
  Created by Ben Appleton (25/11/02)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pde_toolbox_bimage.h"
#include "pde_toolbox_LSTB.h"
//#include "pde.h"
#include "liarp.h"
//#include "liarlmts.h"


/* lblur
Blurs an arbitrary dimensional image with a symmetric row of Pascal's triangle
 */
int lblur(
	DBL_TYPE * in,					/* The input image */
	int * dim_buf,					/* The image dimensions */
	int dim_length,					/* The number of image dimesions */
	int blurring,					/* The amount of blurring to perform */
	DBL_TYPE *out					/* The blurred output image */
)
{
	BVECT * dim;
	BIMAGE * in_image;
	int i, num_pixels;

	/* Set up a BVECT to describe the image dimensions */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));

	/* Construct image, copy data (casting to float) */
	in_image = BIMAGE_constructor_double(in, dim);

	/* Blur the image */
	BIMAGE_blur(in_image, in_image, blurring);

	/* Display any LSTB_error or LSTB_debugages */
	lreadLSTBmsgs();

	/* Cast-copy blurred image to output */
	num_pixels = BVECT_prod(dim);
	for (i = 0; i < num_pixels; i++) {
		out[i] = (DBL_TYPE)in_image->buf[i];
	}

	/* Free everything */
	BIMAGE_destructor(in_image);
	BVECT_destructor(dim);

	return 0;
}

/* laosblur
Blurs an arbitrary dimensional image with an AOS (IIR-like) Gaussian approximation
 */
int laosblur(
	DBL_TYPE * in,					/* The input image */
	int * dim_buf,					/* The image dimensions */
	int dim_length,					/* The number of image dimesions */
	DBL_TYPE blurring,					/* The amount of blurring to perform */
	DBL_TYPE *out					/* The blurred output image */
)
{
	BVECT * dim;
	BIMAGE * in_image;
	int i, num_pixels;

	/* Set up a BVECT to describe the image dimensions */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));

	/* Construct image, copy data (casting to float) */
	in_image = BIMAGE_constructor_double(in, dim);

	/* Blur the image */
	BIMAGE_aosblur(in_image, in_image, (float)blurring);

	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Cast-copy blurred image to output */
	num_pixels = BVECT_prod(dim);
	for (i = 0; i < num_pixels; i++) {
		out[i] = (DBL_TYPE)in_image->buf[i];
	}

	/* Free everything */
	BIMAGE_destructor(in_image);
	BVECT_destructor(dim);

	return 0;
}

/* labs_grad
	Computes the absolute gradient of a BIMAGE
*/
int labs_grad(
	DBL_TYPE * in,					/* The input image */
	int * dim_buf,					/* The image dimensions */
	int dim_length,					/* The number of image dimesions */
	DBL_TYPE *out					/* The gradient image */
)
{
	BVECT * dim;
	BIMAGE * in_image;
	int i, num_pixels;

	/* Set up a BVECT to describe the image dimensions */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));

	/* Construct image, copy data (casting to float) */
	in_image = BIMAGE_constructor_double(in, dim);

	/* Call the gradient function */
	BIMAGE_abs_grad(in_image, in_image);

	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Cast-copy blurred image to output */
	num_pixels = BVECT_prod(dim);
	for (i = 0; i < num_pixels; i++) {
		out[i] = (DBL_TYPE)in_image->buf[i];
	}

	/* Free everything */
	BIMAGE_destructor(in_image);
	BVECT_destructor(dim);

	return 0;
}


/* lrad_grad
	Computes the radial gradient of a BIMAGE
*/
int lrad_grad(
	DBL_TYPE * in,					/* The input image */
	int * dim_buf,					/* The image dimensions */
	int dim_length,					/* The number of image dimesions */
	int * centre_buf,				/* Centre for radial gradient */
	DBL_TYPE *out					/* The gradient image */
)
{
	BVECT * dim, * centre;
	BIMAGE * in_image;
	int i, num_pixels;

	/* Set up a BVECT to describe the image dimensions and centre */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));
	centre = BVECT_constructor(dim_length);
	memcpy(centre->buf, centre_buf, dim_length*sizeof(int));

	/* Construct image, copy data (casting to float) */
	in_image = BIMAGE_constructor_double(in, dim);

	/* Call the gradient function */
	BIMAGE_rad_grad(in_image, in_image, centre);

	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Cast-copy blurred image to output */
	num_pixels = BVECT_prod(dim);
	for (i = 0; i < num_pixels; i++) {
		out[i] = (DBL_TYPE)in_image->buf[i];
	}

	/* Free everything */
	BIMAGE_destructor(in_image);
	BVECT_destructor(dim);
	BVECT_destructor(centre);

	return 0;
}


/* - lradialweighting:
	Weight an image by 1/distance^(N-1) (distance to specified centre point, N is image dimension).

*/
int lradial_weighting(
	DBL_TYPE * in,					/* The input image */
	int * dim_buf,					/* The image dimensions */
	int dim_length,					/* The number of image dimesions */
	int * centre_buf,				/* The centre point */
	DBL_TYPE * out					/* The output image */
)
{
	BVECT * dim, * centre;
	BIMAGE * in_image;
	int i, num_pixels;

	/* Set up a BVECT to describe the image dimensions and centre */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));
	centre = BVECT_constructor(dim_length);
	memcpy(centre->buf, centre_buf, dim_length*sizeof(int));

	/* Construct image, copy data (casting to float) */
	in_image = BIMAGE_constructor_double(in, dim);

	/* Call the radial weighting function */
	BIMAGE_radial_weighting(in_image, centre);

	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Cast-copy image to output */
	num_pixels = BVECT_prod(dim);
	for (i = 0; i < num_pixels; i++) {
		out[i] = (DBL_TYPE)in_image->buf[i];
	}

	/* Free everything */
	BIMAGE_destructor(in_image);
	BVECT_destructor(dim);
	BVECT_destructor(centre);
	
	return 0;
}


/* lmetrify
	Converts a gradient image to a metric image.
*/
int lmetrify(
	DBL_TYPE * in,					/* The input image */
	int * dim_buf,					/* The image dimensions */
	int dim_length,					/* The number of image dimesions */
	DBL_TYPE normalisation,			/* The intensity normalisation factor */
	int p,							/* The power applied to the gradient values */
	DBL_TYPE epsilon,				/* Constant length-penalty coefficient */
	DBL_TYPE *out					/* The metric image */
)
{
	BVECT * dim;
	BIMAGE * in_image;
	int i, num_pixels;

	/* Set up a BVECT to describe the image dimensions */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));

	/* Construct image, copy data (casting to float) */
	in_image = BIMAGE_constructor_double(in, dim);

	/* Call the gradient function */
	BIMAGE_metrify(in_image, in_image, normalisation, p, epsilon);

	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Cast-copy metric image to output */
	num_pixels = BVECT_prod(dim);
	for (i = 0; i < num_pixels; i++) {
		out[i] = (DBL_TYPE)in_image->buf[i];
	}

	/* Free everything */
	BIMAGE_destructor(in_image);
	BVECT_destructor(dim);

	return 0;
}


/* ltensor_metric
	Converts a blurred image to a tensor metric image.
*/
int ltensor_metric(
	DBL_TYPE * in,					/* The input image */
	int * dim_buf,					/* The image dimensions */
	int dim_length,					/* The number of image dimesions */
	DBL_TYPE normalisation,			/* The intensity normalisation factor */
	int p,							/* The power applied to the gradient values */
	DBL_TYPE epsilon,				/* Constant length-penalty coefficient */
	DBL_TYPE *abuf,					/* The SPD metric field g = */
	DBL_TYPE *bbuf,					/* [a b] */
	DBL_TYPE *cbuf					/* [b c] */
)
{
	BVECT * dim;
	BIMAGE * in_image;
	METRIC2D * metric;
	int i, num_pixels;

	/* Set up a BVECT to describe the image dimensions */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));

	/* Construct the metric structure */
	metric = (METRIC2D *)malloc(sizeof(METRIC2D));
	metric->a = BIMAGE_constructor(dim);
	metric->b = BIMAGE_constructor(dim);
	metric->c = BIMAGE_constructor(dim);
	/* aniso_ratios not necessary! */

	/* Construct image, copy data (casting to float) */
	in_image = BIMAGE_constructor_double(in, dim);

	/* Call the metric creating function */
	BIMAGE_tensor_metric(in_image, metric, (float)normalisation, p, (float)epsilon);

	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Copy metric image to output */
	num_pixels = BVECT_prod(dim);
	for (i = 0; i < num_pixels; i++) {
		abuf[i] = (DBL_TYPE)metric->a->buf[i];
		bbuf[i] = (DBL_TYPE)metric->b->buf[i];
		cbuf[i] = (DBL_TYPE)metric->c->buf[i];
	}

	/* Free everything */
	BIMAGE_destructor(in_image);
	BVECT_destructor(dim);
	BIMAGE_destructor(metric->a);
	BIMAGE_destructor(metric->b);
	BIMAGE_destructor(metric->c);
	free((void *)metric);

	return 0;
}

/* lline_metric
	Converts a blurred line image to a tensor metric image.
*/
int lline_metric(
	DBL_TYPE * in,					/* The input image */
	int * dim_buf,					/* The image dimensions */
	int dim_length,					/* The number of image dimesions */
	DBL_TYPE normalisation,			/* Normalisation factor */
	DBL_TYPE epsilon,				/* Euclidean length penalty */
	DBL_TYPE *abuf,					/* The SPD metric field g = */
	DBL_TYPE *bbuf,					/* [a b] */
	DBL_TYPE *cbuf					/* [b c] */
)
{
	BVECT * dim;
	BIMAGE * in_image;
	METRIC2D * metric;
	int i, num_pixels;

	/* Set up a BVECT to describe the image dimensions */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));

	/* Construct the metric structure */
	metric = (METRIC2D *)malloc(sizeof(METRIC2D));
	metric->a = BIMAGE_constructor(dim);
	metric->b = BIMAGE_constructor(dim);
	metric->c = BIMAGE_constructor(dim);
	/* aniso_ratios not necessary! */

	/* Construct image, copy data (casting to float) */
	in_image = BIMAGE_constructor_double(in, dim);

	/* Call the metric creating function */
	BIMAGE_line_metric(in_image, metric, (float)normalisation, (float)epsilon);

	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Copy metric image to output */
	num_pixels = BVECT_prod(dim);
	for (i = 0; i < num_pixels; i++) {
		abuf[i] = (DBL_TYPE)metric->a->buf[i];
		bbuf[i] = (DBL_TYPE)metric->b->buf[i];
		cbuf[i] = (DBL_TYPE)metric->c->buf[i];
	}

	/* Free everything */
	BIMAGE_destructor(in_image);
	BVECT_destructor(dim);
	BIMAGE_destructor(metric->a);
	BIMAGE_destructor(metric->b);
	BIMAGE_destructor(metric->c);
	free((void *)metric);

	return 0;
}


/* ldiffusivity
	Converts a gradient image to a diffusivity image.
*/
int ldiffusivity(
	DBL_TYPE * in,					/* The input image */
	int * dim_buf,					/* The image dimensions */
	int dim_length,					/* The number of image dimesions */
	DBL_TYPE lambda,				/* Soft threshold for gradient significance */
	DBL_TYPE *out					/* The diffusivity image */
)
{
	BVECT * dim;
	BIMAGE * in_image;
	int i, num_pixels;

	/* Set up a BVECT to describe the image dimensions */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));

	/* Construct image, copy data (casting to float) */
	in_image = BIMAGE_constructor_double(in, dim);

	/* Call the gradient function */
	BIMAGE_diffusivity(in_image, in_image, lambda);
	
	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Cast-copy diffusivity image to output */
	num_pixels = BVECT_prod(dim);
	for (i = 0; i < num_pixels; i++) {
		out[i] = (DBL_TYPE)in_image->buf[i];
	}

	/* Free everything */
	BIMAGE_destructor(in_image);
	BVECT_destructor(dim);

	return 0;
}
