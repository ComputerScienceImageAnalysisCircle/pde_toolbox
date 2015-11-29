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
 * File:		contmaxflow.c
 *
 * Written by:	Ben Appleton
 *				ITEE, The University of Queensland
 *
 * Date:		July 2003
 *
 * Copyright:	Intended for open-source distribution
 *
*/

/*********************************************************************************************
 contmaxflow.c
 ------

  DESCRIPTION:
  Simulates a PDE which converges to a continuous maximum flow.

  HISTORY:
  Created by Ben Appleton (14/7/03)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

//#include <pthread.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pde_toolbox_bimage.h"
#include "pde_toolbox_LSTB.h"
#include "pde_toolbox_defs.h"

#include "liarp.h"
//#include "liarlmts.h"
//#include "ioliar.h"

/*
#define CONTMAXFLOW_DEBUG
*/

/**************************** DEFINITIONS *********************************/
/* The safety factor in simulation timestep */
#define CONTMAXFLOW_TIMESTEP_SAFETY_FACTOR 0.9


/**************************************************************************/
/* Multithreaded? Comment out to compile without threads */

//#define CONTMAXFLOW_MULTITHREADED

/**************************************************************************/
#ifdef CONTMAXFLOW_MULTITHREADED
#include <pthread.h>
/* How many cpus, and hence how many threads to use? */
	#define CONTMAXFLOW_NUM_THREADS global_num_threads
	static int global_num_threads;

	/* If the processing was evenly spread over the image domain
	   we could cut up work into CONTMAXFLOW_NUM_THREADS equal packets
	   and process them independently.  But different work packets will
	   take differing amounts of time, due to differences in the amount
	   of work in each packets and possible differences in the processors'
	   loads.

	   So, whenever a thread requests a new work packet we pass it a
	   constant fraction of the remaining work (subject to min/max
	   bounds on the size of a work packet).  In thisa way we minimise
	   the number of work packets (and corresponding overheads) while
	   attempting to keep all threads/cpus busy right up until the
	   end of computation.

	   The following number tells us what fraction of the `optimal'
	   work chunk we should take.
	   0.1 is very conservative, 0.5 is aggressive, 1.0 is the maximum.
	*/
	#define CONTMAXFLOW_LOAD_BALANCE_FACTOR 0.1

	/* Work packet bounds (in numbers of rows to process) */
	#define CONTMAXFLOW_MIN_WORK_PACKET 100
	#define CONTMAXFLOW_MAX_WORK_PACKET 10000
#endif

/* An array of the RLE representations for each row */
typedef struct {
	int * start;
	int * end;
	int length;
} RLE;

/* Processor thread states */
#define CONTMAXFLOW_UPDATE_P 1
#define CONTMAXFLOW_UPDATE_F 2
#define CONTMAXFLOW_CONSTRAIN_F 3
#define CONTMAXFLOW_WAIT 4
#define CONTMAXFLOW_END 5

/****************************** GLOBALS *******************************************/

/* To simplify thread parameter passing, store data in a global struct */
struct {
	BIMAGE * g;					/* Metric image */
	BIMAGE * P;					/* Potential image */
	BIMAGE * * F;				/* Flow velocity components */
	BVECT * dim;
	BVECT * face_dim;
	RLE * P_rle_array;			/* RLE representation of used nodes */
	RLE * * F_rle_array;		/* RLE representation of used edges */
	float timestep;				/* Simulation timestep */
	int num_rows;
	int num_pixels;
} contmaxflow_data;

#ifdef CONTMAXFLOW_MULTITHREADED
	/* Synchronise changes to state of algorithm */
	static pthread_mutex_t algorithm_state_mutex;

	/* The processor threads */
	static pthread_t * proc_threads;

	/* Control variables: Main thread in waiting mode */
	static pthread_mutex_t main_wakeup_mutex;
	static pthread_cond_t main_wakeup_cv;
	static int num_proc_waiting;

	/* Control variables: Processor thread in waiting mode */
	static pthread_mutex_t proc_wakeup_mutex;
	static pthread_cond_t proc_wakeup_cv;
#endif

/* What is the current job? */
static char proc_thread_todo;

/* Which row we are up to in processing the current iteration */
static int current_row;

/*************************** FUNCTION PROTOTYPES *********************************/

/* High-level function to control multithreading */
static void update(
	void
);

#ifdef CONTMAXFLOW_MULTITHREADED
/* Obtain and process work packets.  Main function
	for threads */
void * thread_function(
	void * arg
);
#endif

/* Low-level (thread-specific) function to update F */
static void thread_update_F(
	int start_row,
	int end_row
);

/* Low-level (thread-specific) function to update P */
static void thread_update_P(
	int start_row,
	int end_row
);

/* Low-level (thread-specific) function to constrain F */
static void thread_constrain_F(
	int start_row,
	int end_row
);

/* RLE representation constructor */
static void construct_RLE_type(
	char * type,					/* The vertex types */
	BVECT * dim,					/* The dimensions of the data */
	RLE * * P_rle_array,			/* The vertex rle arrays */
	RLE * * * F_rle_array			/* The edge rle arrays */
);

/* RLE representation destructor */
static void destroy_RLE_type(
	BVECT * dim,					/* The dimensions of the space */
	RLE * P_rle_array,			/* A pointer to the array of vertex RLEs */
	RLE * * F_rle_array			/* A pointer to the arrays of edge RLEs */
);

/************************ FUNCTION IMPLEMENTATIONS *******************************/

/* - contmaxflow:
*/
int contmaxflow(
	BIMAGE * g,					/* The isotropic but non-homogeneous metric */
	char * type,				/* The type of each vertex - normal, source or sink */
	float simulationtime,		/* The total time to run the simulation */
	BIMAGE * P,					/* The input/output P field */
	BIMAGE * * F,				/* The input/output F field */
	pde_hook_func * cbf,		/* callback e.g for display, can be NULL */
	double * Pbuf,				/* debug pressure buffer, can be NULL */
	double * Fx,				/* debug velocity buffer, can be NULL */
	double * Fy,				/* debug velocity buffer, can be NULL */
	double * Fz,				/* debug velocity buffer, can be NULL */
	BVECT * out_dim,			/* Dimensions of imview window */
	int period,					/* callback called every N */
	int num_threads
)
{
	int i;
	float time;
	int step;
	BVECT * coord;

#ifdef CONTMAXFLOW_MULTITHREADED
	/* Setting global parameter */
	global_num_threads = num_threads;
#endif

	printf("contmaxflow.c - Entering\n");

	#ifdef USE_DIRECTION_WEIGHTINGS
		printf(
			"Using direction weightings: [%f, %f, %f, %f]\n",
			direction_weighting[0],
			direction_weighting[1],
			direction_weighting[2],
			direction_weighting[3]
		);
	#endif

	/* Store pointers to images in global data structure */
	contmaxflow_data.g = g;
	contmaxflow_data.P = P;
	contmaxflow_data.F = F;

	/* Compute commonly used variables */
	contmaxflow_data.num_pixels = BVECT_prod(g->dim);
	contmaxflow_data.timestep = (float)CONTMAXFLOW_TIMESTEP_SAFETY_FACTOR / (float)sqrt((double)g->dim->length*1.0);

	/* Compute the global number of rows in image (axis 0) */
	contmaxflow_data.num_rows = 1;
	for (i = 1; i < g->dim->length; ++i) {
		contmaxflow_data.num_rows *= g->dim->buf[i];
	}

	/* Establish dimensions - reduce work packet overhead */
	contmaxflow_data.dim = g->dim;
	contmaxflow_data.face_dim = BVECT_constructor(g->dim->length);
	memcpy(contmaxflow_data.face_dim->buf, g->dim->buf, g->dim->length * sizeof(int));
	contmaxflow_data.face_dim->buf[0] = 1;

	/* Allocate working coordinate */
	coord = BVECT_constructor(g->dim->length);

	/* Create RLE representations (axis 0) of the normal vertices and edges */
	construct_RLE_type(
		type,
		g->dim,
		&(contmaxflow_data.P_rle_array),
		&(contmaxflow_data.F_rle_array)
	);

	/* Thread initialisation */
#ifdef CONTMAXFLOW_MULTITHREADED
{
	pthread_attr_t attr;

	printf("Setting up processor threads\n");

	/* Initialise mutices and condition variables */
	pthread_mutex_init(&algorithm_state_mutex, NULL);
	pthread_mutex_init(&main_wakeup_mutex, NULL);
	pthread_cond_init(&main_wakeup_cv, NULL);
	pthread_mutex_init(&proc_wakeup_mutex, NULL);
	pthread_cond_init(&proc_wakeup_cv, NULL);

	/* Set up the algorithm state */
	proc_thread_todo = CONTMAXFLOW_WAIT;
	num_proc_waiting = 0;

	/* Create each thread */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	proc_threads = (pthread_t *)malloc(CONTMAXFLOW_NUM_THREADS * sizeof(pthread_t));
	for (i = 0; i < CONTMAXFLOW_NUM_THREADS; ++i) {
		pthread_create(&proc_threads[i], &attr, thread_function, (void *)&i);
	}
	pthread_attr_destroy(&attr);

	/* Block and wait for acknowledgement from threads.  When all have
	  acknowledged and are waiting, continue */
	while(LSTB_TRUE) {
		/* Are all processors waiting yet? */
		pthread_mutex_lock(&main_wakeup_mutex);
		if (num_proc_waiting == CONTMAXFLOW_NUM_THREADS) {
			pthread_mutex_unlock(&main_wakeup_mutex);
			break;
		} else {
			/* Wait for signal from thread */
			pthread_cond_wait(&main_wakeup_cv, &main_wakeup_mutex);
			pthread_mutex_unlock(&main_wakeup_mutex);
		}
	}
}
#endif

	/* The simulation */
	if (cbf) printf("\n");
	for (
		time = 0.0f, step = 0;
		time < simulationtime;
		time += contmaxflow_data.timestep , ++step
	) {
		/* Update the auxiliary potential P */
		proc_thread_todo = CONTMAXFLOW_UPDATE_P;
		update();

		/* Update the flow F */
		proc_thread_todo = CONTMAXFLOW_UPDATE_F;
		update();

		/* Constrain saturated velocities */
		proc_thread_todo = CONTMAXFLOW_CONSTRAIN_F;
		update();

		/* call callback every period */
		if (cbf && (Pbuf || Fx || Fy || Fz)  && (step % period == 0)) {
			int p;
			int stat_high = 0, stat_low = 0;
			double prop;

			for (p = 0 ; p < contmaxflow_data.num_pixels ; ++p) {
				int output_index;

				/* Compute the corresponding output index */
				intToBvect(p, coord, g->dim);
				output_index = bvectToInt(coord, out_dim);

				/* NOTE: If statements in inner loop are inefficient */
				/* Should be recoded as a memcpy()'s for stride 1 segments */
				/* Only for debugging so don't worry about it */
				if (Pbuf)
					Pbuf[output_index] = (double) (P->buf[p]);
				if (Fx)
					Fx[output_index] = (double) (F[0]->buf[p]);
				if (Fy)
					Fy[output_index] = (double) (F[1]->buf[p]);
				if (Fz)
					Fz[output_index] = (double) (F[2]->buf[p]);
				/* Weird statistic: */
				if (P->buf[p] > 0.95)
					stat_high++; /* count the number of pixels close to 1.0 */
				if (P->buf[p] < 0.05)
					stat_low++; /* count the number of pixels close to 0.0 */
			}

			cbf(); /* call the callback */

			/* Progress bar */
			prop = (stat_high + stat_low)/((double)contmaxflow_data.num_pixels);
			printf(
				"\rTime = %2.2f of %2.2f ; [>.95= %d <.05 = %d sum= %d num_pixel= %d ; prop=%g]",
				time,
				simulationtime,
				stat_high,
				stat_low,
				stat_high +
				stat_low,
				contmaxflow_data.num_pixels,
				prop
			);
			fflush(stdout);
			/*  Convergence criteria too easy to reach
			if ((time > 10) && (prop > 0.99))
			    break;
			*/

			/* To make a video of the simulation, save pressures to file */
			/*
			{
				char filename[255];
				int ptype = DBL_CODE;
				long nx, ny, nz;
				int return_val;
				nx = g->dim->buf[0];
				ny = g->dim->buf[1];
				nz = 1;

				sprintf(filename, "contmaxflow%3.3i_z", step);
				return_val = lsave_zimage(
					filename,
					0,
					&nx, &ny, &nz,
					0, 0, 0,
					&ptype,
                	(void **)(&(P->buf)),
					IS_BIGENDIAN);
				if (return_val != 0) printf("lsave_zimage failed!\n");
#ifdef UNDEFINED
				int lsave_zimage(const char *fname, int numcomp, long *nx, long *ny,
                 	long *nz, long *x0, long *y0, long *z0, int *ptype,
                 	void **pixbuff, const endianness byteorder);
#endif
			}
			*/
		}
	}

	/* Thread finalisation */
#ifdef CONTMAXFLOW_MULTITHREADED
{
	/* Tell threads to terminate */
	printf("Shutting down processor threads\n");
	proc_thread_todo = CONTMAXFLOW_END;
	update();

	/* Clean up */
	pthread_mutex_destroy(&algorithm_state_mutex);
	pthread_mutex_destroy(&main_wakeup_mutex);
	pthread_cond_destroy(&main_wakeup_cv);
	pthread_mutex_destroy(&proc_wakeup_mutex);
	pthread_cond_destroy(&proc_wakeup_cv);
	free((void *)proc_threads);
}
#endif

	/* Free all allocated memory */
	BVECT_destructor(coord);
	BVECT_destructor(contmaxflow_data.face_dim);

	destroy_RLE_type(
		g->dim,
		contmaxflow_data.P_rle_array,
		contmaxflow_data.F_rle_array
	);

	/* Fix up the command line */
	if (cbf) printf("                                                      \n");

	printf("contmaxflow.c - Exiting\n");

	return 0;
}


/* - update
	thisa function branches into threads and synchronises them again.
*/
static void update(
)
{
	/* Set up the work for the threads to do before we create them */
	current_row = 0;

	/* Perform requested processing */
#ifdef CONTMAXFLOW_MULTITHREADED
	int i, status;

#ifdef CONTMAXFLOW_DEBUG
	printf("Main: Broadcasting threads to start work\n");
#endif
	/* Broadcast-signal threads to start working */
	num_proc_waiting = 0;
	pthread_mutex_lock(&proc_wakeup_mutex);
	pthread_cond_broadcast(&proc_wakeup_cv);
	pthread_mutex_unlock(&proc_wakeup_mutex);
#ifdef CONTMAXFLOW_DEBUG
	printf("Main: Broadcast complete\n");
#endif

	if (proc_thread_todo == CONTMAXFLOW_END) {
		printf("Main: Waiting for join\n");
		/* Threads are terminating, join */
		for (i = 0; i < CONTMAXFLOW_NUM_THREADS; ++i) {
			pthread_join(proc_threads[i], (void **)&status);
		}
		printf("Main: Join successful\n");
	} else {
		/* Block and wait for signal from threads.  When all are
	  	waiting, iteration is completed */
#ifdef CONTMAXFLOW_DEBUG
		printf("Main: Waiting for processor threads\n");
#endif
		while(LSTB_TRUE) {
			/* Check when threads have finished */
			pthread_mutex_lock(&main_wakeup_mutex);
#ifdef CONTMAXFLOW_DEBUG
			printf("num_proc_waiting = %i\n", num_proc_waiting);
#endif
			if (num_proc_waiting == CONTMAXFLOW_NUM_THREADS) {
#ifdef CONTMAXFLOW_DEBUG
				printf("Main: Processor threads done\n");
#endif
				pthread_mutex_unlock(&main_wakeup_mutex);
				break;
			} else {
				/* Wait for signal from processor thread */
				pthread_cond_wait(&main_wakeup_cv, &main_wakeup_mutex);
#ifdef CONTMAXFLOW_DEBUG
				printf("Main: Processor thread finished\n");
#endif
				pthread_mutex_unlock(&main_wakeup_mutex);
			}
		}
	}
#else
	/* Bypass the threading */
	switch(proc_thread_todo) {
		case CONTMAXFLOW_UPDATE_P:
			thread_update_P(0, contmaxflow_data.num_rows - 1);
			break;
		case CONTMAXFLOW_UPDATE_F:
			thread_update_F(0, contmaxflow_data.num_rows - 1);
			break;
		case CONTMAXFLOW_CONSTRAIN_F:
			thread_constrain_F(0, contmaxflow_data.num_rows - 1);
			break;
		default:
			break;
	}
#endif
}

#ifdef CONTMAXFLOW_MULTITHREADED
/* thread_function:
	thisa is the main function for each thread.  It repeatedly
	requests a work packet and carries it out.  When no more
	work packets are available the thread waits or exits.
*/
void * thread_function(
	void * arg					/* Unused */
    ) {

#ifdef CONTMAXFLOW_DEBUG
	int thread_number = *((int *)arg);
#endif
        
	/* Loop through iterations */
	while(proc_thread_todo != CONTMAXFLOW_END) {
		/* Loop through work packets */
		while(proc_thread_todo != CONTMAXFLOW_WAIT) {
			int start_row, end_row;
			char iteration_finished = LSTB_FALSE;
			int work_left;
			int work_amount;

			/*****************************************************/
			/* Grab a work packet! */
			pthread_mutex_lock(&algorithm_state_mutex);

			/* Any work left to do? */
			if (current_row >= contmaxflow_data.num_rows) {
				iteration_finished = LSTB_TRUE;
			}

			/* Take a fraction of the remaining work
				Apply bounds.
				If exceeds remaining work, truncate.
		 	*/
			work_left = contmaxflow_data.num_rows - current_row;
			work_amount = CONTMAXFLOW_LOAD_BALANCE_FACTOR * work_left / CONTMAXFLOW_NUM_THREADS;
			if (work_amount < CONTMAXFLOW_MIN_WORK_PACKET) {
				work_amount = CONTMAXFLOW_MIN_WORK_PACKET;
			}
			if (work_amount > CONTMAXFLOW_MAX_WORK_PACKET) {
				work_amount = CONTMAXFLOW_MAX_WORK_PACKET;
			}

			/* Determine the corresponding range of rows to process */
			start_row = current_row;
			end_row = current_row + work_amount - 1;
			if (end_row > contmaxflow_data.num_rows - 1) {
				end_row = contmaxflow_data.num_rows - 1;
			}

			/* Now, change current_row to reflect that thisa work has been apportioned. */
			current_row = end_row + 1;

			pthread_mutex_unlock(&algorithm_state_mutex);
			/*****************************************************/

			if (iteration_finished) break;

#ifdef CONTMAXFLOW_DEBUG
			printf("Thread %i processing rows %i-%i\n", thread_number, start_row, end_row);
#endif
			/* Which update are we performing? */
			switch(proc_thread_todo) {
				case CONTMAXFLOW_UPDATE_P:
					thread_update_P(start_row, end_row);
					break;
				case CONTMAXFLOW_UPDATE_F:
					thread_update_F(start_row, end_row);
					break;
				case CONTMAXFLOW_CONSTRAIN_F:
					thread_constrain_F(start_row, end_row);
					break;
				default:
					break;
			}
		}

#ifdef CONTMAXFLOW_DEBUG
		printf("Thread %i halting until next call\n", thread_number);
#endif

		pthread_mutex_lock(&proc_wakeup_mutex);

#ifdef CONTMAXFLOW_DEBUG
		printf("Thread %i incrementing num waiting\n", thread_number);
#endif

		/* Let main thread know we're done, and halt until next signal */
		pthread_mutex_lock(&main_wakeup_mutex);
		++num_proc_waiting;
		pthread_cond_signal(&main_wakeup_cv);
		pthread_mutex_unlock(&main_wakeup_mutex);

		pthread_cond_wait(&proc_wakeup_cv, &proc_wakeup_mutex);
		pthread_mutex_unlock(&proc_wakeup_mutex);

#ifdef CONTMAXFLOW_DEBUG
		printf("Thread %i has woken up\n", thread_number);
#endif
	}

	/* Finished */
	pthread_exit(NULL);
}
#endif

/* - thread_update_F
	Update the flow F according to F_t = -grad P (explicit with given timestep).
*/
static void thread_update_F(
	int start_row,
	int end_row
)
{
	int j, k, stride;
	BVECT * face_coord;

	/* Allocate working coordinates and precompute commonly used variables */
	face_coord = BVECT_constructor(contmaxflow_data.dim->length);

	/* For each edge, update according to pressure gradient */
	stride = 1;
	for (j = 0; j < contmaxflow_data.dim->length; j++) {		/* Along each dimension */
		/* Compute the index step between pixels along thisa axis */
		if (j > 0) stride *= contmaxflow_data.dim->buf[j - 1];

		if (j == 0) {
			/* For all the given set of 0-rows */
			for (
				k = start_row, intToBvect(k, face_coord, contmaxflow_data.face_dim);
				k <= end_row;
				++k, BVECT_inc(face_coord, contmaxflow_data.face_dim)
			) {
				int base_index, run;

				/* Compute base index of row */
				base_index = k * contmaxflow_data.dim->buf[0];

				/* For each run in thisa row */
				for (
					run = 0;
					run < contmaxflow_data.F_rle_array[j][k].length;
					++run
				) {
					int index, start, end;

					start = contmaxflow_data.F_rle_array[j][k].start[run] + base_index;
					end = contmaxflow_data.F_rle_array[j][k].end[run] + base_index;

					/* For all edges in thisa run */
					for (index = start; index <= end; ++index) {
						/* The update */
						contmaxflow_data.F[j]->buf[index] -=
							contmaxflow_data.timestep *
							(contmaxflow_data.P->buf[index + 1] - contmaxflow_data.P->buf[index]);
					}
				}
			}
		} else {
			/* For all 0-rows */
			for (
				k = start_row, intToBvect(k, face_coord, contmaxflow_data.face_dim);
				k <= end_row;
				++k, BVECT_inc(face_coord, contmaxflow_data.face_dim)
			) {
				int base_index, run;

				/* Compute base index of row */
				base_index = k * contmaxflow_data.dim->buf[0];

				/* For each run in thisa row */
				for (
					run = 0;
					run < contmaxflow_data.F_rle_array[j][k].length;
					++run
				) {
					int index, start, end;

					start = contmaxflow_data.F_rle_array[j][k].start[run] + base_index;
					end = contmaxflow_data.F_rle_array[j][k].end[run] + base_index;

					/* For all edges in thisa run */
					for (index = start; index <= end; ++index) {
						/* The update */
						contmaxflow_data.F[j]->buf[index] -=
							contmaxflow_data.timestep *
							(contmaxflow_data.P->buf[index + stride] - contmaxflow_data.P->buf[index]);
					}
				}
			}
		}
	}

	/* Free allocated memory */
	BVECT_destructor(face_coord);
}


/* - thread_update_P
	Update the auxiliary potential P according to P_t = -div F
	(explicit with given timestep).
*/
static void thread_update_P(
	int start_row,					/* The range of rows to update */
	int end_row
)
{
	int i, j, k, stride;
	BVECT * face_coord;

	/* Set up working coordinates */
	face_coord = BVECT_constructor(contmaxflow_data.dim->length);

	/* Accumulate changes to pressure along each dimension
		- uses the additive separation of div
		- ignore source/sink
	*/
	stride = 1;
	for (j = 0; j < contmaxflow_data.dim->length; j++) {		/* Along each dimension */
		/* Compute the index step between pixels along thisa axis */
		if (j > 0) stride *= contmaxflow_data.dim->buf[j - 1];

		if (j == 0) {
			/* For all the given set of 0-rows */
			for (
				k = start_row, intToBvect(k, face_coord, contmaxflow_data.face_dim);
				k <= end_row;
				++k, BVECT_inc(face_coord, contmaxflow_data.face_dim)
			) {
				int base_index, run;

				base_index = k * contmaxflow_data.dim->buf[0];

				/* Loop through each run in thisa row */
				for (
					run = 0;
					run < contmaxflow_data.P_rle_array[k].length;
					++run
				) {
					int index, start, end;

					start = contmaxflow_data.P_rle_array[k].start[run];
					end = contmaxflow_data.P_rle_array[k].end[run];

					/* Explicitly handle boundary cases */
					if (start == 0) {
						start++;
						index = base_index;
						contmaxflow_data.P->buf[index] +=
							-contmaxflow_data.timestep
							* contmaxflow_data.F[j]->buf[index];
					}
					if (end == contmaxflow_data.dim->buf[0] - 1) {
						end--;
						index = base_index + contmaxflow_data.dim->buf[0] - 1;
						contmaxflow_data.P->buf[index] +=
							contmaxflow_data.timestep
							* contmaxflow_data.F[j]->buf[index - 1];
					}

					/* Fast general case */
					start += base_index;
					end += base_index;
					for (index = start; index <= end; ++index) {
						contmaxflow_data.P->buf[index] +=
							contmaxflow_data.timestep
							* (- (contmaxflow_data.F[j]->buf[index] - contmaxflow_data.F[j]->buf[index - 1]));
					}
				}
			}
		} else {
			for (
				k = start_row, intToBvect(k, face_coord, contmaxflow_data.face_dim);
				k <= end_row;
				++k, BVECT_inc(face_coord, contmaxflow_data.face_dim)
			) {
				int base_index, run;

				base_index = k * contmaxflow_data.dim->buf[0];

				/* Check if we're on the j-boundary */
				if (face_coord->buf[j] <= 0) {
					/* Loop through each run in thisa row */
					for (
						run = 0;
						run < contmaxflow_data.P_rle_array[k].length;
						++run
					) {
						int index, start, end;

						start = contmaxflow_data.P_rle_array[k].start[run];
						end = contmaxflow_data.P_rle_array[k].end[run];

						/* Fast general case */
						index = base_index + start;
						for (i = start; i <= end; ++i, ++index) {
							contmaxflow_data.P->buf[index] +=
								- contmaxflow_data.timestep
								* contmaxflow_data.F[j]->buf[index];
						}
					}
				} else if (face_coord->buf[j] >= contmaxflow_data.dim->buf[j] - 1) {
					/* Loop through each run in thisa row */
					for (
						run = 0;
						run < contmaxflow_data.P_rle_array[k].length;
						++run
					) {
						int index, start, end;

						start = contmaxflow_data.P_rle_array[k].start[run];
						end = contmaxflow_data.P_rle_array[k].end[run];

						/* Fast general case */
						index = base_index + start;
						for (i = start; i <= end; ++i, ++index) {
							contmaxflow_data.P->buf[index] += 
								contmaxflow_data.timestep 
								* contmaxflow_data.F[j]->buf[index - stride];
						}
					}
				} else {
					/* Loop through each run in thisa row */
					for (
						run = 0;
						run < contmaxflow_data.P_rle_array[k].length;
						++run
					) {
						int index, start, end;

						start = contmaxflow_data.P_rle_array[k].start[run];
						end = contmaxflow_data.P_rle_array[k].end[run];

						/* Fast general case */
						index = base_index + start;
						for (i = start; i <= end; ++i, ++index) {
							contmaxflow_data.P->buf[index] += 
								contmaxflow_data.timestep
								* (- (contmaxflow_data.F[j]->buf[index] - contmaxflow_data.F[j]->buf[index - stride]));
						}
					}
				}
			}
		}
	}

	/* Free working memory */
	BVECT_destructor(face_coord);
}


/* - thread_constrain_F
	Update the flow F to satisfy |F| <= g, renormalising where necessary.
	Only treats outward flows, ensuring in-place processing is consistent.
*/
static void thread_constrain_F(
	int start_row,					/* The range of rows to update */
	int end_row
)
{
	int i, j, k;
	int * stride;
	float * flow;
	float scale;
	BVECT * coord, * face_coord;

	/* Setup */
	coord = BVECT_constructor(contmaxflow_data.dim->length);
	face_coord = BVECT_constructor(contmaxflow_data.dim->length);

	/* Precompute the strides for each dimension */
	stride = (int *)malloc(contmaxflow_data.dim->length * sizeof(int));
	stride[0] = 1;
	for (i = 1; i < contmaxflow_data.dim->length; i++) {
		stride[i] = stride[i-1] * contmaxflow_data.dim->buf[i - 1];
	}

	/* Allocate memory for the flow vector */
	flow = (float *)malloc(contmaxflow_data.dim->length * sizeof(float));

	/* For each row along the 0-axis */
	for (
		k = start_row, intToBvect(k, face_coord, contmaxflow_data.face_dim);
		k <= end_row;
		++k, BVECT_inc(face_coord, contmaxflow_data.face_dim)
	) {
		int base_index, run;
		char border;

		base_index = k * contmaxflow_data.dim->buf[0];

		/* Border checking? */
		border = LSTB_FALSE;
		for (j = 1; j < contmaxflow_data.dim->length; j++) {
			if (face_coord->buf[j] == 0 || face_coord->buf[j] == contmaxflow_data.dim->buf[j] - 1) {
				border = LSTB_TRUE;
				break;
			}
		}

		if (border) {	/* In thisa case we check boundaries everywhere */
			/* Loop through each run in thisa row */
			for (
				run = 0;
				run < contmaxflow_data.P_rle_array[k].length;
				++run
			) {
				int index, start, end;

				start = contmaxflow_data.P_rle_array[k].start[run];
				end = contmaxflow_data.P_rle_array[k].end[run];
				for (
					i = start, index = base_index + start;
					i <= end;
					++i, ++index
				) {
					float abs_flow;

					/* Find maximum outward flow at thisa node */
					for (j = 0; j < contmaxflow_data.dim->length; j++) {
						float back_flow, forward_flow;

						/* Extract the flow backwards and forwards */
						if (j == 0) {
							if (i == 0) {
								back_flow = 0;
							} else {
								back_flow = - contmaxflow_data.F[j]->buf[index - stride[j]];
							}
							if (i >= contmaxflow_data.dim->buf[j] - 1) {
								forward_flow = 0;
							} else {
								forward_flow = contmaxflow_data.F[j]->buf[index];
							}
						} else {
							if (face_coord->buf[j] == 0) {
								back_flow = 0;
							} else {
								back_flow = - contmaxflow_data.F[j]->buf[index - stride[j]];
							}
							if (face_coord->buf[j] == contmaxflow_data.dim->buf[j] - 1) {
								forward_flow = 0;
							} else {
								forward_flow = contmaxflow_data.F[j]->buf[index];
							}
						}

						/* Take the maximum positive flow */
						flow[j] = LSTB_MAX(back_flow, forward_flow);
						flow[j] = LSTB_MAX(0.0f, flow[j]);
					}

					/* Test for saturation |F| > g */
					abs_flow = 0.0f;
					for (j = 0; j < contmaxflow_data.dim->length; j++) {
						#ifdef USE_DIRECTION_WEIGHTINGS
							flow[j] /= direction_weighting[j];
						#endif
						abs_flow += (float)LSTB_SQR(flow[j]);
					}

					/* If saturated, reduce flow magnitude */
					if (abs_flow > LSTB_SQR(contmaxflow_data.g->buf[index])) {
						abs_flow = (float)sqrt(abs_flow);

						/* Rescale the abs outward flow vector */
						scale = contmaxflow_data.g->buf[index] / abs_flow;
						for (j = 0; j < contmaxflow_data.dim->length; j++) {
							flow[j] *= scale;
							#ifdef USE_DIRECTION_WEIGHTINGS
								flow[j] *= direction_weighting[j];
							#endif
						}

						/* Threshold component-wise against rescaled abs vector */
						for (j = 0; j < contmaxflow_data.dim->length; j++) {
							/* Only touch outward flows */
							if (j == 0) {
								if (i > 0) {
									if (-contmaxflow_data.F[j]->buf[index - stride[j]] > flow[j]) {
										contmaxflow_data.F[j]->buf[index - stride[j]] = -flow[j];
									}
								}
								if (i < contmaxflow_data.dim->buf[j] - 1) {
									if (contmaxflow_data.F[j]->buf[index] > flow[j]) {
										contmaxflow_data.F[j]->buf[index] = flow[j];
									}
								}
							} else {
								if (face_coord->buf[j] > 0) {
									if (-contmaxflow_data.F[j]->buf[index - stride[j]] > flow[j]) {
										contmaxflow_data.F[j]->buf[index - stride[j]] = -flow[j];
									}
								}
								if (face_coord->buf[j] < contmaxflow_data.dim->buf[j] - 1) {
									if (contmaxflow_data.F[j]->buf[index] > flow[j]) {
										contmaxflow_data.F[j]->buf[index] = flow[j];
									}
								}
							}
						}
					}
				}
			}
		} else {
			/* Loop through each run in thisa row */
			for (
				run = 0;
				run < contmaxflow_data.P_rle_array[k].length;
				++run
			) {
				int index, start, end;

				start = contmaxflow_data.P_rle_array[k].start[run];
				end = contmaxflow_data.P_rle_array[k].end[run];

				/* Explicitly handle boundary cases */
				if (start == 0) {
					start++;
					index = base_index;

					/* Find maximum outward flow at thisa node */
					for (j = 0; j < contmaxflow_data.dim->length; j++) {
						float back_flow, forward_flow;

						/* Extract the flow backwards and forwards */
						if (j == 0) {
							back_flow = 0;
							forward_flow = contmaxflow_data.F[j]->buf[index];
						} else {
							if (face_coord->buf[j] == 0) {
								back_flow = 0;
							} else {
								back_flow = - contmaxflow_data.F[j]->buf[index - stride[j]];
							}
							if (face_coord->buf[j] == contmaxflow_data.dim->buf[j] - 1) {
								forward_flow = 0;
							} else {
								forward_flow = contmaxflow_data.F[j]->buf[index];
							}
						}

						/* Take the maximum positive flow */
						flow[j] = LSTB_MAX(back_flow, forward_flow);
						flow[j] = LSTB_MAX(0.0f, flow[j]);
					}
				}
				if (end == contmaxflow_data.dim->buf[0] - 1) {
					end--;
					index = base_index + contmaxflow_data.dim->buf[0] - 1;

					/* Find maximum outward flow at thisa node */
					for (j = 0; j < contmaxflow_data.dim->length; j++) {
						float back_flow, forward_flow;

						/* Extract the flow backwards and forwards */
						if (j == 0) {
							back_flow = - contmaxflow_data.F[j]->buf[index - stride[j]];
							forward_flow = 0;
						} else {
							if (face_coord->buf[j] == 0) {
								back_flow = 0;
							} else {
								back_flow = - contmaxflow_data.F[j]->buf[index - stride[j]];
							}
							if (face_coord->buf[j] == contmaxflow_data.dim->buf[j] - 1) {
								forward_flow = 0;
							} else {
								forward_flow = contmaxflow_data.F[j]->buf[index];
							}
						}

						/* Take the maximum positive flow */
						flow[j] = LSTB_MAX(back_flow, forward_flow);
						flow[j] = LSTB_MAX(0.0f, flow[j]);
					}
				}

				/* Fast general case */
				start += base_index;
				end += base_index;
				for (index = start; index <= end; ++index) {
					float abs_flow;

					/* Find maximum outward flow at thisa node */
					for (j = 0; j < contmaxflow_data.dim->length; ++j) {
						float back_flow, forward_flow;

						/* Extract the flow backwards and forwards */
						back_flow = - contmaxflow_data.F[j]->buf[index - stride[j]];
						forward_flow = contmaxflow_data.F[j]->buf[index];

						/* Take the maximum positive flow */
						flow[j] = LSTB_MAX(back_flow, forward_flow);
						flow[j] = LSTB_MAX(0.0f, flow[j]);
					}

					/* Test for saturation |F| > g */
					abs_flow = 0.0f;
					for (j = 0; j < contmaxflow_data.dim->length; j++) {
#ifdef USE_DIRECTION_WEIGHTINGS
						flow[j] /= direction_weighting[j];
#endif
						abs_flow += (float)LSTB_SQR(flow[j]);
					}

					/* If saturated, reduce flow magnitude */
					if (abs_flow > LSTB_SQR(contmaxflow_data.g->buf[index])) {
						abs_flow = (float)sqrt(abs_flow);

						/* Rescale the abs outward flow vector */
						scale = contmaxflow_data.g->buf[index] / abs_flow;
						for (j = 0; j < contmaxflow_data.dim->length; j++) {
							flow[j] *= scale;
#ifdef USE_DIRECTION_WEIGHTINGS
							flow[j] *= direction_weighting[j];
#endif
						}

						/* Threshold component-wise against rescaled abs vector */
						for (j = 0; j < contmaxflow_data.dim->length; j++) {
							/* Only touch outward flows */
							if (-contmaxflow_data.F[j]->buf[index - stride[j]] > flow[j]) {
								contmaxflow_data.F[j]->buf[index - stride[j]] = -flow[j];
							}
							if (contmaxflow_data.F[j]->buf[index] > flow[j]) {
								contmaxflow_data.F[j]->buf[index] = flow[j];
							}
						}
					}
				}
			}
		}
	}

	/* Free allocated memory */
	BVECT_destructor(coord);
	BVECT_destructor(face_coord);
	free((void *)stride);
	free((void *)flow);
}


/* - construct_RLE_type:
	Construct the RLE representation (axis 0) of the normal vertices and edges
*/
static void construct_RLE_type(
	char * type,					/* The vertex types */
	BVECT * dim,					/* The dimensions of the data */
	RLE * * P_rle_array,			/* The vertex rle arrays */
	RLE * * * F_rle_array			/* The edge rle arrays */
)
{
	BVECT * coord, * face_coord, * face_dim;
	int num_pixels, num_face_pixels;
	int i, j, k;
	int stride;
	RLE * rle;

	/* Allocate working space */
	coord = BVECT_constructor(dim->length);
	face_coord = BVECT_constructor(dim->length);
	face_dim = BVECT_constructor(dim->length);

	/* Allocate a working RLE of maximum length */
	rle = (RLE*)malloc(sizeof(RLE));
	rle->length = BVECT_max(dim);
	rle->start = (int*)malloc(rle->length * sizeof(int));
	rle->end = (int*)malloc(rle->length * sizeof(int));

	/* Determine the face-dimensions of the volume */
	memcpy(face_dim->buf, dim->buf, dim->length * sizeof(int));
	face_dim->buf[0] = 1;

	/* Precompute commonly used variables */
	num_pixels = BVECT_prod(dim);
	num_face_pixels = BVECT_prod(face_dim);

	/* Allocate the base arrays for the RLE linked lists */
	*P_rle_array = (RLE*)malloc(num_face_pixels * sizeof(RLE));
	*F_rle_array = (RLE**)malloc(dim->length * sizeof(RLE *));
	for (j = 0; j < dim->length; ++j) {
		(*F_rle_array)[j] = (RLE*)malloc(num_face_pixels * sizeof(RLE));
	}

	/* Compute RLEs of normal vertices */
	for (k = 0, BVECT_zero(face_coord); k < num_face_pixels; ++k, BVECT_inc(face_coord, face_dim)) {
		int index;
		char run_open;

		index = k * dim->buf[0];

		run_open = LSTB_FALSE;
		rle->length = 0;			/* Initially empty */
		for (i = 0; i < dim->buf[0]; ++i, ++index) {
			if (!run_open && type[index] == MAXFLOW_NORMAL) {
				/* Open up a run */
				run_open = LSTB_TRUE;
				rle->start[rle->length] = i;
			} else if (run_open && type[index] != MAXFLOW_NORMAL) {
				/* Close off a run */
				run_open = LSTB_FALSE;
				rle->end[rle->length] = i - 1;
				++(rle->length);
			}
		}
		/* Close off an unfinished run */
		if (run_open) {
			run_open = LSTB_FALSE;
			rle->end[rle->length] = dim->buf[0] - 1;
			++(rle->length);
		}

		/* Now copy across to the P_rle_array */
		(*P_rle_array)[k].length = rle->length;
		if (rle->length > 0) {
			(*P_rle_array)[k].start = (int*)malloc(rle->length * sizeof(int));
			(*P_rle_array)[k].end = (int*)malloc(rle->length * sizeof(int));
			memcpy((*P_rle_array)[k].start, rle->start, rle->length * sizeof(int));
			memcpy((*P_rle_array)[k].end, rle->end, rle->length * sizeof(int));
		}

		/* DEBUGGING */
		/*
		if (rle->length > 0) {
			BVECT_print(face_coord); printf(":\t");
			for (i = 0; i < rle->length; i++) {
				printf("(%i, %i)", (*P_rle_array)[k].start[i], (*P_rle_array)[k].end[i]);
			}
			printf("\n");
		}
		*/
	}

	/* Compute RLEs of normal edges */
	stride = 1;
	for (j = 0; j < dim->length; ++j) {
		if (j != 0) {
			stride *= dim->buf[j - 1];
		}

		if (j == 0) {
			for (k = 0, BVECT_zero(face_coord); k < num_face_pixels; ++k, BVECT_inc(face_coord, face_dim)) {
				int index;
				char run_open;

				index = k * dim->buf[0];

				run_open = LSTB_FALSE;
				rle->length = 0;			/* Initially empty */
				for (i = 0; i < dim->buf[0] - 1; ++i, ++index) {
					if (!run_open && (type[index] == MAXFLOW_NORMAL || type[index + 1] == MAXFLOW_NORMAL)) {
						/* Open up a run */
						run_open = LSTB_TRUE;
						rle->start[rle->length] = i;
					} else if (run_open && !(type[index] == MAXFLOW_NORMAL || type[index + 1] == MAXFLOW_NORMAL)) {
						/* Close off a run */
						run_open = LSTB_FALSE;
						rle->end[rle->length] = i - 1;
						++(rle->length);
					}
				}
				/* Close off an unfinished run */
				if (run_open) {
					run_open = LSTB_FALSE;
					rle->end[rle->length] = dim->buf[0] - 2;
					++(rle->length);
				}

				/* Now copy across to the F_rle_array */
				(*F_rle_array)[j][k].length = rle->length;
				if (rle->length > 0) {
					(*F_rle_array)[j][k].start = (int*)malloc(rle->length * sizeof(int));
					(*F_rle_array)[j][k].end = (int*)malloc(rle->length * sizeof(int));
					memcpy((*F_rle_array)[j][k].start, rle->start, rle->length * sizeof(int));
					memcpy((*F_rle_array)[j][k].end, rle->end, rle->length * sizeof(int));
				}
			}
		} else {
			for (k = 0, BVECT_zero(face_coord); k < num_face_pixels; ++k, BVECT_inc(face_coord, face_dim)) {
				int index;
				char run_open;

				/* Check bounds */
				if (face_coord->buf[j] == dim->buf[j] - 1) {
					/* Open up a null RLE */
					(*F_rle_array)[j][k].length = 0;

					continue;
				}

				index = k * dim->buf[0];

				run_open = LSTB_FALSE;
				rle->length = 0;			/* Initially empty */
				for (i = 0; i < dim->buf[0]; ++i, ++index) {
					if (!run_open && (type[index] == MAXFLOW_NORMAL || type[index + stride] == MAXFLOW_NORMAL)) {
						/* Open up a run */
						run_open = LSTB_TRUE;
						rle->start[rle->length] = i;
					} else if (run_open && !(type[index] == MAXFLOW_NORMAL || type[index + stride] == MAXFLOW_NORMAL)) {
						/* Close off a run */
						run_open = LSTB_FALSE;
						rle->end[rle->length] = i - 1;
						++(rle->length);
					}
				}
				/* Close off an unfinished run */
				if (run_open) {
					run_open = LSTB_FALSE;
					rle->end[rle->length] = dim->buf[0] - 1;
					++(rle->length);
				}

				/* Now copy across to the F_rle_array */
				(*F_rle_array)[j][k].length = rle->length;
				if (rle->length > 0) {
					(*F_rle_array)[j][k].start = (int*)malloc(rle->length * sizeof(int));
					(*F_rle_array)[j][k].end = (int*)malloc(rle->length * sizeof(int));
					memcpy((*F_rle_array)[j][k].start, rle->start, rle->length * sizeof(int));
					memcpy((*F_rle_array)[j][k].end, rle->end, rle->length * sizeof(int));
				}
			}
		}
	}

	/* Deallocate everything */
	BVECT_destructor(coord);
	BVECT_destructor(face_coord);
	BVECT_destructor(face_dim);
	free((void *)rle->start);
	free((void *)rle->end);
	free((void *)rle);
}


/* - destroy_RLE_type:
	Destroy the RLE arrays of vertex and edge types
*/
static void destroy_RLE_type(
	BVECT * dim,					/* The dimensions of the space */
	RLE * P_rle_array,				/* A pointer to the array of vertex RLEs */
	RLE * * F_rle_array				/* A pointer to the arrays of edge RLEs */
)
{
	BVECT * face_dim;
	int num_face_pixels;
	int i, j;

	/* Allocate working space */
	face_dim = BVECT_constructor(dim->length);

	/* Determine the face-dimensions of the volume */
	memcpy(face_dim->buf, dim->buf, dim->length * sizeof(int));
	face_dim->buf[0] = 1;

	/* Precompute commonly used variables */
	num_face_pixels = BVECT_prod(face_dim);

	/* Free each RLE */
	for (i = 0; i < num_face_pixels; ++i) {
		if (P_rle_array[i].length > 0) {
			free((void *)P_rle_array[i].start);
			free((void *)P_rle_array[i].end);
		}
	}

	free((void *)(P_rle_array));
	for (j = 0; j < dim->length; ++j) {
		for (i = 0; i < num_face_pixels; ++i) {
			if (F_rle_array[j][i].length > 0) {
				free((void *)F_rle_array[j][i].start);
				free((void *)F_rle_array[j][i].end);
			}
		}
		free((void *)F_rle_array[j]);
	}
	free((void *)F_rle_array);

	/* Deallocate everything */
	BVECT_destructor(face_dim);
}
