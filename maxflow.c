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
 * File:		maxflow.c
 *
 * Written by:		Ben Appleton
 *
 * Date:		July 2002
 *
 *
 *
*/

/*********************************************************************************************
 maxflow.c
 ------

  DESCRIPTION:
  An N-D image segmentation tool based on maximum flows.

  REFERENCE:
  Goldberg and Tarjan's pre-flow push algorithm

  HISTORY:
  Created by Ben Appleton (July 2003)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pde_toolbox_bimage.h"
#include "pde_toolbox_LSTB.h"
#include "pde_toolbox_defs.h"

/* Define the index/flag of a non-existent vertex */
#define MAXFLOW_NULL_VERTEX -1

/* Define DEBUG to enable printfs for debugging */
/* #define MAXFLOW_DEBUG */

/* The number of bits to use in quantising floating point */
#define MAXFLOW_QUANTBITS 24

/********************************* Data Structures *********************************/
/* The integer priority queue structure */
typedef struct {
	/* List bases and ends, (indexed by height) */
	int * first;
	int * last;

	/* Arrays of previous/next node indices for each pixel (indexed by pixel) */
	int * prev;
	int * next;

	/* Current height of highest active, source height, and maximum allowable height */
	int curHeight;
	int sourceHeight;
	int maxHeight;
} PQUEUE_INT;

/* Encapsulate all data structures of interest - the class object */
typedef struct {
	/* Indexing structures */
	PQUEUE_INT * pqueue;		/* The priority queue */
	BVECT * dim;				/* The dimensions of our space */
	int * stride;				/* The strides along each axis */
	char * active;				/* The active set indexed by pixel index */
	int * height;				/* The height array, accessed by image index */

	/* Flow-related variables */
	int * cap;					/* Capacities and flows stored by edge index */
	int * flow;
	int * excess;				/* The excess at each vertex */

	/* Lists of source and sink indices */
	int * source_indices;
	int num_source_indices;
	int * sink_indices;
	int num_sink_indices;

	/* Label pixels (source, sink, normal) */
	char * type;
} MAXFLOWDS;


/************************************* FUNCTION PROTOTYPES **************************************/
/* Quantise metric values to integer capacities */
static float quantise_metric(
	BIMAGE * g,							/* Floating point metric (already radially weighted) */
	int * cap							/* Integer capacities (output) */
);

/* Removes the highest active vertex from the queue */
static int queuePullHighestActive(
	MAXFLOWDS * thisa
);

/* Removes the specified vertex from the queue */
static void queuePull(
	int v,
	MAXFLOWDS *thisa
);

/* Adds the specified vertex to the queue (tests for duplicates) */
static void queuePush(
	int v,
	MAXFLOWDS * thisa
);

/* Initialise the queue */
static void queueInit(
	MAXFLOWDS *thisa
);

#if 0
static void queuePrint(
	MAXFLOWDS * thisa
    );
#endif

static void globalRelabel(
	MAXFLOWDS * thisa
);

static void gap_relabel(
	int height,
	MAXFLOWDS * thisa
);

static void compute_lists(
	MAXFLOWDS * thisa
);


/* Main function:
*/
int maxflow(
	BIMAGE * g,							/* The input metric image (no need to radially weight!) */
	char * type,						/* The type of each vertex */
	BIMAGE * P,							/* The pressure function (output) */
	BIMAGE * * F						/* The flows (scaled to floating point) */
)
{
	/* The class object */
	MAXFLOWDS * thisa;

	/* Miscellaneous variables */
	int i, j, k, iterations, num_pixels;
	float scale;							/* Quantisation scale */
	BVECT * coord;
	int height;

	#ifdef USE_DIRECTION_WEIGHTINGS
	printf("Using direction weightings: [%f, %f, %f, %f]\n",
		direction_weighting[0], direction_weighting[1], direction_weighting[2], direction_weighting[3]);
	#endif

	/* Precompute the number of pixels == vertices */
	num_pixels = BVECT_prod(g->dim);

	/* Allocate working coordinates */
	coord = BVECT_constructor(g->dim->length);

	/* Build the class object */
	thisa = (MAXFLOWDS *)malloc(sizeof(MAXFLOWDS));

	/* Allocate memory */
	/* Initially 0 flow/excess */
	thisa->dim = BVECT_constructor(g->dim->length);
	memcpy(thisa->dim->buf, g->dim->buf, g->dim->length * sizeof(int));
	thisa->cap = (int *)malloc(num_pixels * g->dim->length * sizeof(int));
	thisa->flow = (int *)calloc(num_pixels * g->dim->length, sizeof(int));
	thisa->excess = (int *)calloc(num_pixels, sizeof(int));
	thisa->height = (int *)malloc(num_pixels * sizeof(int));
	thisa->active = (char *)calloc(num_pixels, sizeof(char));				/* Initially false */
	thisa->type = (char *)calloc(num_pixels, sizeof(char));					/* Default to MAXFLOW_NORMAL */
	thisa->stride = (int *)malloc(g->dim->length * sizeof(int));
	thisa->type = type;

	/* Precompute the possible strides */
	thisa->stride[0] = 1;
	for (j = 1; j < g->dim->length; j++) {
		thisa->stride[j] = thisa->stride[j-1] * g->dim->buf[j-1];
	}

	/* Initialise the source and sink index list */
	compute_lists(thisa);

	/* Sanity check on input */
	if (thisa->num_source_indices == 0 || thisa->num_sink_indices == 0) {
		LSTB_error("The seed image must specify some source (+1) and sink (-1) vertices!\n");
		return 1;
	}

	/* Construct the priority queue */
	thisa->pqueue = (PQUEUE_INT *)malloc(sizeof(PQUEUE_INT));
	thisa->pqueue->maxHeight = 10 * BVECT_max(g->dim);						/* Heuristic! */
	thisa->pqueue->sourceHeight = thisa->pqueue->maxHeight / 2;
	thisa->pqueue->last = (int *)malloc((thisa->pqueue->maxHeight+1) * sizeof(int));
	thisa->pqueue->first = (int *)malloc((thisa->pqueue->maxHeight+1) * sizeof(int));
	thisa->pqueue->prev = (int *)malloc(num_pixels * sizeof(int));
	thisa->pqueue->next = (int *)malloc(num_pixels * sizeof(int));
	queueInit(thisa);

	/* Quantise metric to integer capacities and convert to edge capacities */
	scale = quantise_metric(g, thisa->cap);

	/* The source should be initialised with maximum excess and the sink with maximally negative excess */
	for (i = 0; i < thisa->num_source_indices; i++) {
		thisa->excess[thisa->source_indices[i]] = 1 << 30;
		thisa->active[thisa->source_indices[i]] = LSTB_TRUE;
	}
	for (i = 0; i < thisa->num_sink_indices; i++) {
		thisa->excess[thisa->sink_indices[i]] = -(1 << 30);
	}

	#ifdef MAXFLOW_DEBUG
		printf("Entering main loop\n");
	#endif

	/* Main loop */
	iterations = 0;
	while(LSTB_TRUE) {
		int u, v;
		int minHeight;

		/* Global relabelling to assign optimum heights (must be run to initialise) */
		if (iterations % (10 * num_pixels) == 0) {
			#ifdef MAXFLOW_DEBUG
				printf("Iteration %i: Global relabelling\n", iterations);
			#endif
			globalRelabel(thisa);
		}

		/* Retrieve the highest active vertex from the queue */
		v = queuePullHighestActive(thisa);

		/* Debugging */
		#ifdef MAXFLOW_DEBUG
			intToBvect(v, coord, thisa->dim);
			printf("%i\t", iterations);
			BVECT_print(coord);
			printf("\n");
			printf("Vertex has height %i\n", thisa->height[v]);
		#endif

		if (v == MAXFLOW_NULL_VERTEX) break;

		/* Extract the coordinate of thisa point */
		intToBvect(v, coord, thisa->dim);

		/* Push flow out of all eligible edges leaving the vertex */
		minHeight = 1 << 30;
		for (j = 0; j < thisa->dim->length; j++) {
			for (k = -1; k <= 1; k += 2) {			/* Backwards or forwards */
				int P, Q, edge;

				#ifdef MAXFLOW_DEBUG
					intToBvect(v, coord, thisa->dim);
				#endif

				/* Does the proposed neighbouring vertex exist? */
				if (coord->buf[j] + k < 0 || coord->buf[j] + k > thisa->dim->buf[j] - 1) continue;

				/* Compute the index of the corresponding vertex */
				u = v + k * thisa->stride[j];

				#ifdef MAXFLOW_DEBUG
					printf("\n");
					printf("Attempting to push flow to ");
					intToBvect(u, coord, thisa->dim);
					BVECT_print(coord);
					printf("\n");
					printf("Vertex has height %i\n", thisa->height[u]);
					if (thisa->active[u]) {
						printf("Vertex is active\n");
					} else {
						printf("Vertex is inactive\n");
					}
					printf("Vertex has type %i\n", thisa->type[u]);
				#endif

				/* Compute the corresponding edge index */
				edge = j + g->dim->length * LSTB_MIN(u, v);

				#ifdef MAXFLOW_DEBUG
					printf("Edge has capacity %i, flow %i\n", thisa->cap[edge], k * thisa->flow[edge]);
				#endif

				/* Q is the maximum amount of flow we can push *out* along thisa edge (in the residual graph) */
				Q = thisa->cap[edge] - k * thisa->flow[edge];
				#ifdef MAXFLOW_DEBUG
					printf("Residual capacity Q = %i\n", Q);
				#endif

				/* Determine P, the potential flow in thisa direction */
				P = LSTB_MIN(Q, thisa->excess[v]);
				#ifdef MAXFLOW_DEBUG
					printf("Possible flow P = %i\n", P);
				#endif

				/* If the edge is feasible and we can push flow out along it, do so */
				if (
					P > 0 && (thisa->height[v] == thisa->height[u] + 1 || (thisa->type[v] == MAXFLOW_SOURCE && thisa->type[u] != MAXFLOW_SOURCE))
				) {
					#ifdef MAXFLOW_DEBUG
						printf("Edge was feasible\n");
					#endif

					/* Push the flow out */
					thisa->flow[edge] += k * P;
					thisa->excess[v] -= P;
					thisa->excess[u] += P;

					/* Enqueue the target vertex if it's not the source or sink */
					if (thisa->type[u] == MAXFLOW_NORMAL) {
						if (!thisa->active[u]) {
							#ifdef MAXFLOW_DEBUG
								printf("Activating secondary vertex\n");
							#endif
							queuePull(u, thisa);
							thisa->active[u] = LSTB_TRUE;
							queuePush(u, thisa);
						}
					}
				}

				/* If thisa edge is not saturated it remains in the residual graph, so update the min-height */
				if (k * thisa->flow[edge] < thisa->cap[edge]) {
					if (thisa->height[u] < minHeight) {
						minHeight = thisa->height[u];
					}
				}

			}
		}

		/* Determine if thisa vertex is still active */
		thisa->active[v] = (thisa->excess[v] > 0 && thisa->type[v] == MAXFLOW_NORMAL);
		/* Assign the point a new height.  Store the old to check for gaps */
		height = thisa->height[v];
		if (thisa->active[v]) {
			thisa->height[v] = minHeight + 1;
		}
		queuePush(v, thisa);

		/* If a gap has formed, perform a gap-relabel */
		if (thisa->height[v] != height) {
			gap_relabel(height, thisa);
		}

		#ifdef MAXFLOW_DEBUG
			if (thisa->active[v]) {
				printf("Primary still active\n");
			} else {
				printf("Primary deactived\n");
			}
			printf("Primary height now %i\n", thisa->height[v]);
		#endif

		iterations++;
	}

	/* Extract cut */
	 globalRelabel(thisa);

	/* Output pressures and velocities */
	if (P != NULL) {
		for (i = 0; i < num_pixels; i++) {
			P->buf[i] = (thisa->height[i] >= thisa->pqueue->sourceHeight ? 1 : 0);
		}
	}
	if (F != NULL) {
		for (j = 0; j < g->dim->length; j++) {
			for (i = 0; i < num_pixels; i++) {
				F[j]->buf[i] = ((float)thisa->flow[j + g->dim->length * i]) / scale;
			}
		}
	}

	/* Free allocated memory */
	BVECT_destructor(thisa->dim);
	free((void *)thisa->stride);
	free((void *)thisa->source_indices);
	free((void *)thisa->sink_indices);
	free((void *)thisa->cap);
	free((void *)thisa->flow);
	free((void *)thisa->excess);
	free((void *)thisa->height);
	free((void *)thisa->active);
	free((void *)thisa->pqueue->first);
	free((void *)thisa->pqueue->last);
	free((void *)thisa->pqueue->prev);
	free((void *)thisa->pqueue->next);
	free((void *)thisa->pqueue);
	free((void *)thisa);

	BVECT_destructor(coord);

	#ifdef MAXFLOW_DEBUG
		printf("Done!\n");
	#endif

	return 0;
}


/* - quantise_metric:
	Quantise metric to integer capacities
*/
static float quantise_metric(
	BIMAGE * g,							/* Floating point metric (already radially weighted) */
	int * cap							/* Integer edges capacities (output) */
)
{
	int i, j, stride, num_pixels;
	float max_g, scale;
	int * c;
	BVECT * coord;

	/* Precompute the number of pixels */
	num_pixels = BVECT_prod(g->dim);

	/* Allocate memory */
	c = (int *)malloc(num_pixels * sizeof(int));
	coord = BVECT_constructor(g->dim->length);

	/* Locate maximum g */
	max_g = 0.0f;
	for (i = 0; i < num_pixels; i++) {
		if (g->buf[i] > max_g) max_g = g->buf[i];
	}

	/* Compute scale parameter */
	scale = ((float)(1 << MAXFLOW_QUANTBITS)) / max_g;

	/* Apply scale parameter to quantise to c */
	for (i = 0; i < num_pixels; i++) {
		c[i] = (int)(g->buf[i] * scale);
	}

	/* Compute edge capacities */
	stride = 1;
	for (j = 0; j < g->dim->length; j++) {
		if (j > 0) stride *= g->dim->buf[j-1];

		for (i = 0, BVECT_zero(coord); i < num_pixels; i++, BVECT_inc(coord, g->dim)) {
			/* Skip non-existent edges */
			if (coord->buf[j] == g->dim->buf[j] - 1) continue;

			/* Assign average capacity to thisa particular edge */
			cap[j + g->dim->length * i] = (c[i] + c[i + stride]) / 2;
			#ifdef USE_DIRECTION_WEIGHTINGS
				cap[j + g->dim->length * i] *= direction_weighting[j];
			#endif
		}
	}

	/* Free allocated memory */
	free((void *)c);
	BVECT_destructor(coord);

	return scale;
}


/* queuePullHighestActive:
	Removes the highest active vertex from the queue
	Returns MAXFLOW_NULL_VERTEX if no active vertices left (complete)
*/
static int queuePullHighestActive(
	MAXFLOWDS * thisa
)
{
	int v;

	/* Determine the highest active vertex v and call queuePull on it */
	if (thisa->pqueue->curHeight >= 0) {
		v = thisa->pqueue->first[thisa->pqueue->curHeight];
		queuePull(v, thisa);
	} else {
		v = MAXFLOW_NULL_VERTEX;
	}

	return v;
}


/* queuePull:
	Removes the specified vertex from the queue
	Must maintain the height of the highest active point
*/
static void queuePull(
	int v,
	MAXFLOWDS *thisa
)
{
	int prev, next, height;

	/* Determine the neighbours (prev, next) of thisa node on its list */
	prev = thisa->pqueue->prev[v];
	next = thisa->pqueue->next[v];

	/* Redirect the neighbours' to skip v */
	if (prev != MAXFLOW_NULL_VERTEX)
		thisa->pqueue->next[prev] = next;
	if (next != MAXFLOW_NULL_VERTEX)
		thisa->pqueue->prev[next] = prev;

	/* Just for safety, set thisa node's neighbour indices to NULLVERTEX */
	thisa->pqueue->next[v] = MAXFLOW_NULL_VERTEX;
	thisa->pqueue->prev[v] = MAXFLOW_NULL_VERTEX;

	/* Update the pointers into both ends of thisa list */
	height = thisa->height[v];
	if (thisa->pqueue->first[height] == v)
		thisa->pqueue->first[height] = next;
	if (thisa->pqueue->last[height] == v)
		thisa->pqueue->last[height] = prev;

	/* Find the next highest active vertex */
	for ( ; thisa->pqueue->curHeight >= 0; thisa->pqueue->curHeight--) {
		int first;

		first = thisa->pqueue->first[thisa->pqueue->curHeight];
		if (first == MAXFLOW_NULL_VERTEX) continue;
		if (!thisa->active[first]) continue;
		break;
	}
}


/* queuePush:
	Push a vertex onto the queue
	Must maintain the height of the highest active point
*/
static void queuePush(
	int v,
	MAXFLOWDS * thisa
)
{
	int height;

	/* Short name for commonly used variable */
	height = thisa->height[v];

	/* If it's higher than the current highest active, update */
	if (thisa->active[v] && height > thisa->pqueue->curHeight) {
		thisa->pqueue->curHeight = height;
	}

	/* Add to list of appropriate height */
	if (thisa->active[v]) {			/* If active, prepend */
		if (thisa->pqueue->first[height] != MAXFLOW_NULL_VERTEX) {
			thisa->pqueue->prev[thisa->pqueue->first[height]] = v;		/* Invertibility */
		}
		thisa->pqueue->next[v] = thisa->pqueue->first[height];
		thisa->pqueue->prev[v] = MAXFLOW_NULL_VERTEX;
		thisa->pqueue->first[height] = v;
		if (thisa->pqueue->last[height] == MAXFLOW_NULL_VERTEX) thisa->pqueue->last[height] = v;
	} else {				/* If !active, append */
		if (thisa->pqueue->last[height] != MAXFLOW_NULL_VERTEX) {
			thisa->pqueue->next[thisa->pqueue->last[height]] = v;		/* Invertibility */
		}
		thisa->pqueue->prev[v] = thisa->pqueue->last[height];
		thisa->pqueue->next[v] = MAXFLOW_NULL_VERTEX;
		thisa->pqueue->last[height] = v;
		if (thisa->pqueue->first[height] == MAXFLOW_NULL_VERTEX) thisa->pqueue->first[height] = v;
	}
}


/* queueInit:
	Sets all queues to empty (first = last = MAXFLOW_NULL_VERTEX)
*/
static void queueInit(
	MAXFLOWDS *thisa
)
{
	int i;

	/* Set all queues to empty */
	for (i = 0; i <= thisa->pqueue->maxHeight; i++) {
		thisa->pqueue->first[i] = MAXFLOW_NULL_VERTEX;
		thisa->pqueue->last[i] = MAXFLOW_NULL_VERTEX;
	}
	thisa->pqueue->curHeight = -1;
}


#if 0 // no longer used, generates a warning...
/* queuePrint:
	Print out the queue - for debugging
*/
static void queuePrint(
	MAXFLOWDS * thisa
)
{
	int height, v;
	BVECT * coord;

	printf("Queue:\n");

	coord = BVECT_constructor(thisa->dim->length);

	for (height = 0; height <= thisa->pqueue->maxHeight; height++) {
		if (thisa->pqueue->first[height] == MAXFLOW_NULL_VERTEX) continue;
		printf("\n%i:\n", height);
		for (v = thisa->pqueue->first[height]; v != MAXFLOW_NULL_VERTEX; v = thisa->pqueue->next[v]) {
			intToBvect(v, coord, thisa->dim);
			BVECT_print(coord);
			printf("\t");
		}
	}
	printf("\n");

	BVECT_destructor(coord);
}

#endif // 0

/* globalRelabel:
	Globally (re\)label by distance function in residual graph
*/
static void globalRelabel(
	MAXFLOWDS * thisa
)
{
	int i, j, k, num_pixels;
	int u, v, edge, height;
	char * onqueue;
	BVECT * coord;

	/* Precompute the number of pixels */
	num_pixels = BVECT_prod(thisa->dim);

	/* Allocate memory */
	coord = BVECT_constructor(thisa->dim->length);
	onqueue = (char *)malloc(num_pixels * sizeof(char));

	/* Wipe the queue first */
	queueInit(thisa);
	memset(onqueue, 0, num_pixels * sizeof(char));

	/* Compute the distances */
	for (i = 0; i < 2; i++) {
		/* Two phases - compute shortest path to sink, then failing that, source */
		if (i == 0) {
			height = 0;
			/* Start with the sink vertices on the queue */
			for (j = 0; j < thisa->num_sink_indices; j++) {
				thisa->height[thisa->sink_indices[j]] = height;
				queuePush(thisa->sink_indices[j], thisa);
				onqueue[thisa->sink_indices[j]] = LSTB_TRUE;
			}
			/* The starting point for the distance function */
			v = thisa->pqueue->first[height];
		} else {
			height = thisa->pqueue->sourceHeight;
			/* Start with the sink vertices on the queue */
			for (j = 0; j < thisa->num_source_indices; j++) {
				thisa->height[thisa->source_indices[j]] = height;
				queuePush(thisa->source_indices[j], thisa);
				onqueue[thisa->source_indices[j]] = LSTB_TRUE;
			}
			/* The starting point for the distance function */
			v = thisa->pqueue->first[height];
		}

		/* Now BFS using the queue */
		while(LSTB_TRUE) {
			/* Extract the coordinate of thisa point */
			intToBvect(v, coord, thisa->dim);

			/* For each neighbour of the current vertex */
			for (j = 0; j < thisa->dim->length; j++) {
				for (k = -1; k <= 1; k += 2) {			/* Backwards or forwards */
					/* Skip non-existent edges */
					if (coord->buf[j] + k < 0 || coord->buf[j] + k > thisa->dim->buf[j] - 1) continue;

					u = v + k * thisa->stride[j];

					/* Ensure that the source is not added to the queue */
					if (thisa->type[u] != MAXFLOW_NORMAL) continue;

					/* Compute the corresponding edge index */
					edge = j + (thisa->dim->length) * LSTB_MIN(u, v);

					/* If thisa edge is not saturated toward u */
					if (-k * thisa->flow[edge] < thisa->cap[edge]) {
						if (!onqueue[u]) {
							onqueue[u] = LSTB_TRUE;
							thisa->height[u] = height + 1;
							queuePush(u, thisa);
						}
					}
				}
			}

			/* Shift to the next vertex, if there is one */
			if (thisa->pqueue->next[v] != MAXFLOW_NULL_VERTEX) {
				v = thisa->pqueue->next[v];
			} else {
				height++;
				if (height == thisa->pqueue->sourceHeight) break;
				if (height == thisa->pqueue->maxHeight) break;
				if (thisa->pqueue->first[height] == MAXFLOW_NULL_VERTEX) break;
				v = thisa->pqueue->first[height];
			}
		}
	}

	/* Free allocated memory */
	BVECT_destructor(coord);
	free((void *)onqueue);
}


/* - gap_relabel:
	Check for a gap at the specified height and relabel if appropriate
*/
static void gap_relabel(
	int height,
	MAXFLOWDS * thisa
)
{
	int i, v;

	/* Only relabel if height is below source */
	if (height >= thisa->pqueue->sourceHeight - 1) return;

	/* Check for a gap */
	if (thisa->pqueue->first[height] != MAXFLOW_NULL_VERTEX) return;

	/* printf("Gap relabelling from height %i\n", height); */

	/* Relabel heights of all higher points */
	for (i = height + 1; i < thisa->pqueue->sourceHeight; i++) {
		/* Consume the specified list */
		while((v = thisa->pqueue->first[i]) != MAXFLOW_NULL_VERTEX) {
			queuePull(v, thisa);
			thisa->height[v] = thisa->pqueue->sourceHeight;
			queuePush(v, thisa);
		}
	}

}


/* - compute_lists:
	Compute the lists of source and sink indices from the vertex types.
	WARNING - Allocates source/sink indices!
*/
static void compute_lists(
	MAXFLOWDS * thisa
)
{
	int i, num_pixels;

	/* Precompute commonly used variables */
	num_pixels = BVECT_prod(thisa->dim);

	/* Scan to obtain the lists */
	/* Two scans (slightly inefficient) */
	thisa->num_source_indices = 0;
	thisa->num_sink_indices = 0;
	for (i = 0; i < num_pixels; i++) {
		if (thisa->type[i] != MAXFLOW_NORMAL) {
			if (thisa->type[i] == MAXFLOW_SOURCE) {
				thisa->num_source_indices++;
			} else {
				thisa->num_sink_indices++;
			}
		}
	}
	thisa->source_indices = (int *)malloc(thisa->num_source_indices * sizeof(int));
	thisa->sink_indices = (int *)malloc(thisa->num_sink_indices * sizeof(int));
	thisa->num_source_indices = 0;
	thisa->num_sink_indices = 0;
	for (i = 0; i < num_pixels; i++) {
		if (thisa->type[i] != MAXFLOW_NORMAL) {
			if (thisa->type[i] == MAXFLOW_SOURCE) {
				thisa->source_indices[thisa->num_source_indices++] = i;
			} else {
				thisa->sink_indices[thisa->num_sink_indices++] = i;
			}
		}
	}
}

