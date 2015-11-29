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
 * File:		lcontmaxflow.c
 *
 * Written by:		Image Analysis Group staff,
 * 			CSIRO Mathematical and Information Sciences.
 *
 * Date:		July 2003
 *
 *
*/

/*********************************************************************************************
 lcontmaxflow.c
 ------

  DESCRIPTION:
  A continuous maximum flow algorithm

  HISTORY:
  Created by Ben Appleton (14/07/03)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

/************************************* Includes **********************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pde_toolbox_bimage.h"
#include "pde_toolbox_LSTB.h"
#include "pde_toolbox_defs.h"
#include "liarp.h"
//#include "liarlmts.h"

/* Perform greedy multiscale for speed ? */
/* #define GREEDY */
#ifdef GREEDY
/* If P is within P_TOLERANCE of binary it is fixed as source/sink */
	#define P_TOLERANCE 0.05
#endif

/* Prevent over subsampling */
static const int MIN_IMAGE_SIZE = 20;

/************************************* Implementations ***************************************/
/* lcontmaxflow
 	A continuous maximum flow algorithm
*/
int lcontmaxflow(
	DBL_TYPE * g,				/* The isotropic but non-homogeneous metric */
	int * dim_buf,				/* The image dimensions */
	int dim_length,
	DBL_TYPE * dbl_type,		/* Label image of node types */
	DBL_TYPE initial_time,		/* Simulation time at first scale */
	DBL_TYPE settling_time,		/* Simulation time at subsequent scales */
	DBL_TYPE * out,				/* The output */
	int num_scales,				/* The number of scales (>1 for multiscale) */
	pde_hook_func * cbf,		/* callback function e.g for display, can be NULL */
	DBL_TYPE * dbgP,			/* debug pressure (can be NULL) */
	DBL_TYPE * dbgFx,			/* debug velocity components (can be NULL) */
	DBL_TYPE * dbgFy,
	DBL_TYPE * dbgFz,
	int period,					/* frequency of call back (every N step...) */
	int num_threads
)
{
	int scale, i, j, num_pixels;
	BVECT * dim, * temp_dim;
	int return_val = LSTB_FALSE;
	BIMAGE * * image_pyramid;
	BIMAGE * * P_pyramid;
	BIMAGE * * * F_pyramid;
	char * * type_pyramid;

	/* Allocate memory */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));
	temp_dim = BVECT_constructor(dim_length);

	/* Compute the number of pixels for reuse later */
	num_pixels = BVECT_prod(dim);

	/* Issue warnings */
	if (num_scales < 1) {
		printf("Warning: num_scales too low, setting to 1\n");
		num_scales = 1;
	}
	if (BVECT_min(dim) / (1<<(num_scales-1)) < MIN_IMAGE_SIZE) {
		printf("Too much downsampling, ignoring multiscale parameters\n");
		num_scales = 1;
	}

	/* Construct the necessary BIMAGEs */
	image_pyramid = (BIMAGE **)malloc(num_scales * sizeof(BIMAGE *));
	image_pyramid[0] = BIMAGE_constructor_double(g, dim);


	/* Create a pyramid from the metric */
	BVECT_copy(temp_dim, dim);
	for (scale = 1; scale < num_scales; scale++) {
		for (j = 0; j < temp_dim->length; j++) {
			temp_dim->buf[j] /= 2;
		}

		image_pyramid[scale] = BIMAGE_constructor(temp_dim);
		BIMAGE_downsample(image_pyramid[scale-1], image_pyramid[scale], 2, LSTB_TRUE);
	}

	/* Construct image pyramids for the P and F fields */
	P_pyramid = (BIMAGE **)malloc(num_scales * sizeof(BIMAGE *));
	/* Construct the P pyramid */
	BVECT_copy(temp_dim, dim);
	for (scale = 0; scale < num_scales; scale++) {
		P_pyramid[scale] = BIMAGE_constructor(temp_dim);
		for (j = 0; j < temp_dim->length; j++) {
			temp_dim->buf[j] /= 2;
		}
	}

	F_pyramid = (BIMAGE ***)malloc(num_scales * sizeof(BIMAGE **));
	/* For each scale */
	BVECT_copy(temp_dim, dim);
	for (scale = 0; scale < num_scales; scale++) {
		/* Allocate the vector of images for thisa scale  */
		F_pyramid[scale] = (BIMAGE **)malloc(dim->length * sizeof(BIMAGE *));

		/* Allocate each component of the F image at thisa scale */
		for (j = 0; j < temp_dim->length; j++) {
			F_pyramid[scale][j] = BIMAGE_constructor(temp_dim);
		}

		/* Compute the dimensions at the next scale */
		for (j = 0; j < temp_dim->length; j++) {
			temp_dim->buf[j] /= 2;
		}
	}

	/* Construct the type pyramid */
	type_pyramid = (char **)malloc(num_scales * sizeof(char *));
	type_pyramid[0] = (char *)malloc(BVECT_prod(dim) * sizeof(char));
	/* Convert the given type image */
	for (i = 0; i < num_pixels; i++) {
		type_pyramid[0][i] = (char)LSTB_SIGN(dbl_type[i]);
	}
	/* Downsample to produce pyramid */
	BVECT_copy(temp_dim, dim);
	for (scale = 1; scale < num_scales; scale++) {
		for (j = 0; j < temp_dim->length; j++) {
			temp_dim->buf[j] /= 2;
		}

		type_pyramid[scale] = (char *)malloc(BVECT_prod(temp_dim) * sizeof(char));

		/* Downsample */
		label_image_downsample(
			type_pyramid[scale - 1],
			image_pyramid[scale - 1]->dim,
			type_pyramid[scale],
			image_pyramid[scale]->dim,
			2
		);
	}

	/* For each scale */
	for (scale = num_scales - 1; scale >= 0; scale--) {
		/* Initialise at lowest scale */
		if (scale == num_scales - 1) {
			/* Discrete maximum flow */
			maxflow(image_pyramid[scale], type_pyramid[scale], P_pyramid[scale], F_pyramid[scale]);

			/* DEBUGGING - initialise P to be source indicator function (for movies) */
			/*
			int num_pixels = BVECT_prod(P_pyramid[scale]->dim);
			for (i = 0; i < num_pixels; ++i) {
				P_pyramid[scale]->buf[i] = (type_pyramid[scale][i] == MAXFLOW_SOURCE) ? 1 : 0;
			}
			*/
		}

		/* Iterate to convergence */
		return_val = contmaxflow(
			image_pyramid[scale],
			type_pyramid[scale],
			((scale == num_scales - 1) ? ((float)initial_time) : ((float)settling_time)),
			P_pyramid[scale],
			F_pyramid[scale],
			cbf,
			dbgP,
			dbgFx,
			dbgFy,
			dbgFz,
			dim,
			period,
			num_threads
		);

		if (return_val) {
			LSTB_error("Error: contmaxflow dropped out early!\n");
			return return_val;
		}

		/* Upsample P, F for next scale */
		if (scale > 0) {
			BIMAGE_upsample(P_pyramid[scale], P_pyramid[scale - 1], 2, LSTB_TRUE);
			for (i = 0; i < dim->length; i++) {
				BIMAGE_upsample(F_pyramid[scale][i], F_pyramid[scale - 1][i], 2, LSTB_TRUE);
			}

			/* GREEDY mode:
				For a fast but greedy version of the multiscale continuous maximum flow,
				we only consider a small band around the minimal surface from the previous
				scale!  thisa produces an O(N ^ D-1) algorithm instead of O(N^ D+1) I believe.
				Only applied to final scale
			*/
			#ifdef GREEDY
			if (scale == 1) {
				/* Threshold P and compute new source/sink set. */
				int num_pixels = BVECT_prod(P_pyramid[scale - 1]->dim);
				for (i = 0; i < num_pixels; ++i) {
					/* Only modify the types of normal nodes */
					if (type_pyramid[scale - 1][i] == MAXFLOW_NORMAL) {
						if (P_pyramid[scale - 1]->buf[i] > (float)(1.0 - P_TOLERANCE)) {
							type_pyramid[scale - 1][i] = MAXFLOW_SOURCE;
						} else if (P_pyramid[scale - 1]->buf[i] < (float)(P_TOLERANCE)) {
							type_pyramid[scale - 1][i] = MAXFLOW_SINK;
						}
					}
					/* Reset P on source/sink nodes after upsampling */
					if (type_pyramid[scale - 1][i] == MAXFLOW_SOURCE) {
						P_pyramid[scale - 1]->buf[i] = 1.0f;
					} else if (type_pyramid[scale - 1][i] == MAXFLOW_SINK) {
						P_pyramid[scale - 1]->buf[i] = 0.0f;
					}
				}
			}
			#endif
		}
	}

	/* Copy output buffer across */
	for (i = 0; i < num_pixels; i++) {
		out[i] = (DBL_TYPE)P_pyramid[0]->buf[i];
	}

	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Free the metric image pyramid */
	for (i = 0; i < num_scales; i++) {
		BIMAGE_destructor(image_pyramid[i]);
	}
	free((void *)image_pyramid);

	/* Free the P image pyramid */
	for (i = 0; i < num_scales; i++) {
		BIMAGE_destructor(P_pyramid[i]);
	}
	free((void *)P_pyramid);
	
	/* Free the F image pyramid */
	for (i = 0; i < num_scales; i++) {
		/* Free each component at thisa scale */
		for (j = 0; j < dim->length; j++) {
			BIMAGE_destructor(F_pyramid[i][j]);
		}
		free((void *)F_pyramid[i]);
	}
	free((void *)F_pyramid);

	/* Free the type image pyramid */
	for (i = 0; i < num_scales; i++) {
		free((void *)type_pyramid[i]);
	}
	free((void *)type_pyramid);

	/* Free working coordinates */
	BVECT_destructor(dim);
	BVECT_destructor(temp_dim);
	fprintf (stderr,"end lcontmaxflow\n");

	return return_val;
}
