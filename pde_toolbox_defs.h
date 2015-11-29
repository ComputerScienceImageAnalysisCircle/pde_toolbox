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
 * File:		pde_toolbox_defs.h
 *
 * Written by:		Ben Appleton
 * 			ITEE, The University of Queensland
 *
 * Date:		October 2002
 *
*/

/*********************************************************************************************
 pde_defs.h
 ------

  DESCRIPTION:
  Definitions for pde functions.

  HISTORY:
  Created by Ben Appleton (10/02)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

#ifndef PDE_DEFS_H
#define PDE_DEFS_H

/* My defines, used by my pde functions */
#ifndef LSTB_SIGN
#define LSTB_SIGN(X) (((X) > 0) ? 1 : ( ((X) < 0) ? -1 : 0 ))
#endif

#ifndef LSTB_ABS
#define LSTB_ABS(X) (((X) >= 0) ? (X) : -(X))
#endif

#ifndef LSTB_SQR
#define LSTB_SQR(X) ((X)*(X))
#endif

#ifndef LSTB_MIN
#define LSTB_MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#endif

#ifndef LSTB_MAX
#define LSTB_MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#endif

#define LSTB_TRUE -1
#define LSTB_FALSE 0

#define LSTB_BIGNUM 1000000
#define LSTB_SMALLNUM 0.000001

/* Define the different level set output types  */
#define LSDIRECT -1						/* Copy level set data directly, no reinitialisation */
#define LSFULL 0						/* Perform a full reinitialisation (slow!) */
#define LSREGION 1						/* Output the internal region */
#define LSCONTOUR 2						/* Output the contour only (fast!) */

/* Define the different speed functions */
#define SETHIAN_SPEED_FUNC 0			/* The speed function from Sethian's book */
#define GAC_SPEED_FUNC 1				/* The speed function for Geodesic Active Contours */
#define GVF_SPEED_FUNC 2				/* The speed function used for Gradient Vector Flows */

/* Stopping criteria */
#define NOSTOPPING 0
#define STOPONMETRIC 1
#define STOPONDISTANCE 2

/* Define how far from the border to place the box level set in init_LS_border */
#define INIT_LS_BORDER_DIST 5.0

#define IM_AOS_BLUR_MAX_TIME_STEP 20.0

/* Bits used to store vertex types in maximum flow algorithm */
#define MAXFLOW_SINK -1
#define MAXFLOW_NORMAL 0
#define MAXFLOW_SOURCE 1


/* TEST CODE
	Hard coded weightings for directions
	For use when data dimensions are not directly comparable
*/
/* #define USE_DIRECTION_WEIGHTINGS */
#ifdef USE_DIRECTION_WEIGHTINGS
static const float direction_weighting[] = {1.0f, 1.0f, 1.0f, 1.0f};
#endif

#endif
