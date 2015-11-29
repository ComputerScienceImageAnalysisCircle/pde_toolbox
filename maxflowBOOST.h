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
 * File:		maxflowBOOST.h
 *
 * Written by:	Ben Appleton
 *				ITEE, The University of Queensland
 *
 * Date:		Feb 2005
 *
*/

/*********************************************************************************************
 maxflowBOOST.h
 ------

  DESCRIPTION:
  Interface to the Boost Graph Library's maxflow code

  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

#ifndef MAXFLOWBOOST_H
#define MAXFLOWBOOST_H

#ifdef __cplusplus
/* Force the function to compile as C-callable */
extern "C" {
#endif

int maxflowBOOST(
	BIMAGE * g,							/* The input metric image (no need to radially weight!) */
	char * type,						/* The type of each vertex */
	BIMAGE * P,							/* The pressure function (output) */
	BIMAGE * * F						/* The flows (scaled to floating point) */
);

#ifdef __cplusplus
}
#endif /* __CPLUSPLUS */
#endif

