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
 * File:		lmaxflow.c
 *
 * Written by:		Image Analysis Group staff,
 * 			CSIRO Mathematical and Information Sciences.
 *
 * Date:		July 2003
 *
 *
*/

/*********************************************************************************************
 lmaxflow.c
 ------

  DESCRIPTION:
  A discrete maximum flow algorithm (Based on pre-flow push, Goldberg + Tarjan)

  HISTORY:
  Created by Ben Appleton (18/07/03)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

/************************************* Includes **********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pde_toolbox_bimage.h"
#include "pde_toolbox_LSTB.h"
//#include "pde.h"
#include "liarp.h"
//#include "liarlmts.h"

#include "maxflowBOOST.h"

/************************************* Implementations ***************************************/
/* lmaxflow
 	A discrete maximum flow algorithm
*/
int lmaxflow(
	DBL_TYPE * g,				/* The isotropic but non-homogeneous metric */
	int * dim_buf,				/* The image dimensions */
	int dim_length,				/* The number of image dimesions */
	DBL_TYPE * dbl_type,		/* Label image of node types */
	DBL_TYPE * out				/* The output */
)
{
	int i, num_pixels;
	BVECT * dim;
	int return_val;
	BIMAGE * g_bimage, * out_bimage;
	char * type;

	/* Set up a BVECT to describe the image dimensions, and the source vector */
	dim = BVECT_constructor(dim_length);
	memcpy(dim->buf, dim_buf, dim_length*sizeof(int));

	/* Compute the number of pixels for reuse later */
	num_pixels = BVECT_prod(dim);

	/* Allocate a node type array and convert from input type */
	type = (char *)malloc(num_pixels * sizeof(char));
	for (i = 0; i < num_pixels; i++) {
		type[i] = (char)dbl_type[i];
	}

	/* Construct the necessary BIMAGEs */
	g_bimage = BIMAGE_constructor_double(g, dim);
	out_bimage = BIMAGE_constructor(dim);

	/* Compute a maximum flow */
	return_val = maxflow(
		g_bimage,
		type,
		out_bimage,
		((BIMAGE * *)NULL)
	);
/*  thisa function does not return the mincut, and is only for timings */
/*
	return_val = maxflowBOOST(
		g_bimage,
		type,
		out_bimage,
		((BIMAGE * *)NULL)
	);
*/

	/* Copy output buffer across */
	for (i = 0; i < num_pixels; i++) {
		out[i] = (DBL_TYPE)out_bimage->buf[i];
	}

	/* Display any LSTB_error or LSTB_debug messages */
	lreadLSTBmsgs();

	/* Free everything */
	BVECT_destructor(dim);
	BIMAGE_destructor(g_bimage);
	BIMAGE_destructor(out_bimage);
	free((void *)type);

	return return_val;
}
