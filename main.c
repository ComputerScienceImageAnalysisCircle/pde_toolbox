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
* Simple stand-alone wrapper for the continuous-max-flow algorithm
* 3D implementation and modifications added 2013 by Patryk Ząbkiewicz
*
* Camille Couprie, 2008
* Patryk Ząbkiewicz, 2013
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include "pde_toolbox_bimage.h"
#include "pde_toolbox_LSTB.h"
#include "pde_toolbox_defs.h"
#include "liarp.h"

void read_image(char * filename, int * img, int *row, int *col);
void read_int_image(int * img, int row, int col);
void writeimage( DBL_TYPE * image, char * filename, int row, int col, int z);
void write_classical_image( DBL_TYPE * image, char * filename, int row, int col, int z);
long voisin(int indice, int num, int rs, int cs, int zs);
void compute_metric(int * img, double * output, int row, int col, int z);
double max(double x, double y);
void read_3dimage(char * filename, int * img, int row, int col, int z, int header_size);
void writeimage3d( DBL_TYPE * image, char * filename, int row, int col, int z);


// external prototypes
extern int lcontmaxflow(
	DBL_TYPE * g,                   /* The isotropic but non-homogeneous metric */
	int * dim_buf,                  /* The image dimensions (up to 4D) */
	int dim_length,                 /* number of dimensions  of the image */
	DBL_TYPE * dbl_type,            /* Label image of node types */
	DBL_TYPE initial_time,          /* Simulation time at first scale */
	DBL_TYPE settling_time,         /* Simulation time at subsequent scales */
	DBL_TYPE * out,                 /* The output */
	int num_scales,                 /* The number of scales (>1 for multiscale) */
	pde_hook_func * cbf,            /* callback function e.g for display, can be NULL */
	DBL_TYPE * dbgP,                /* debug pressure (can be NULL) */
	DBL_TYPE * dbgFx,               /* debug velocity components (can be NULL) */
	DBL_TYPE * dbgFy,
	DBL_TYPE * dbgFz,
	int period,                     /* frequency of call back (every N step...) */
	int num_threads                 /* multithreaded implementation */
	);

extern int lsepgaussf_double(DBL_TYPE *input, /* input image */
							 int nx,
							 int ny, 
							 int nz,
							 double sigma,
							 DBL_TYPE *smoothedim); /* output image */

// must be somewhat close to zero, but not too much
// used in writeimage, for a kind of 3-level thresholding
#define MYEPSILON 0.05

static int debug=0;

#define dbgprintf if(debug) printf

/* =============================*/
int main(int argc, char **argv)
	/* =============================*/ 
{


	if (argc < 9)
	{
		printf("usage: %s image_to_segment image_source image_sink image_out row col depth header_size [-d]\n", argv[0]);
		return 0;
	}
	if (argc == 10) {
		debug = 1;
	}

	int row = atoi(argv[5]);
	int col= atoi(argv[6]);
	int z = atoi(argv[7]);
	int header_size = atoi(argv[8]);

	//printf("%d\n", argc);

	dbgprintf("Row=%d, Col=%d Depth=%d\n",row,col,z);

	int * img = (int *)malloc(sizeof(int)* row * col * z);
	read_3dimage(argv[1], img, row, col, z, header_size);


	//segmentation of bronchial lumen
	//need to invert the image
	for(int i = 0; i < row*col*z; ++i) {
		img[i] = 255 - img[i];
	}


	DBL_TYPE * g;                                             /* The isotropic but non-homogeneous metric */
	g = (DBL_TYPE *)malloc(sizeof(DBL_TYPE)*row*col*z);

	compute_metric(img, g, row, col, z);

	//double maxi=0;
	//for(int i = 0; i < row*col*z; ++i) {
	//	if(g[i] > maxi) maxi = g[i];
	//}
	//for(int i = 0; i < row*col*z; ++i) {
	//	g[i] = maxi - g[i] + 1.5;
	//}
	for(int i = 0; i < row*col*z; ++i) {
		g[i] = g[i] - 0.05;
	}

	dbgprintf("read 3d images\n");

	const int dim_length = 3; //make it 3d                                       /* 2 dimentionnal problem */
	int dim_buf[dim_length];	                            /* The image dimensions */
	dim_buf[0] = row;
	dim_buf[1] = col;
	dim_buf[2] = z;

	DBL_TYPE * dbl_type;		                            /* Label image of node types :
															sources pixels > 0 ,sink pixels < 0, others = 0 */
	dbl_type = (DBL_TYPE *)malloc(sizeof(DBL_TYPE)*row*col*z);

	/* ----------------Enter here 2 image maps for the source and the sink----------------------------------- */
	int * source = (int *)malloc(sizeof(int)* row * col * z);
	int * sink = (int *)malloc(sizeof(int)* row * col * z);

	read_3dimage(argv[2], source, row, col, z, header_size);
	read_3dimage(argv[3], sink, row, col, z, header_size);

	dbgprintf("DONE: read 3d images\n");

	int img_size = row*col*z;
	int i;
	for (i = 0; i< img_size; ++i)
	{
		if (source[i] == 1) 
			dbl_type[i] = 1;
		else if (sink[i] == 1) 
			dbl_type[i] = -1;
		else dbl_type[i] = 0;
	}

	/*-------------------------------------------------------------------------------------------------------*/

	int num_scales = 3;		                            /* The number of scales (>1 for multiscale) */
	int nb_iter = 2*max(row,col)/pow(2.f, num_scales);
	dbgprintf("Number of iterations = %d\n", nb_iter);
	DBL_TYPE initial_time = nb_iter  * 0.9 * sqrt(2.f);	    /* Simulation time at first scale */
	DBL_TYPE settling_time = initial_time;            	    /* Simulation time at subsequent scales */
	dbgprintf("initial time = settling time = %g\n", initial_time);
	DBL_TYPE * out ;                                          /* The output */
	out = (DBL_TYPE *)malloc(sizeof(DBL_TYPE)*row*col*z);	  

	int period = 5;	                                    /* frequency of call back (every N step...) */
	int num_threads = 2;                                      /* 1 thread for 1 processor*/

	lcontmaxflow(g, dim_buf, dim_length, dbl_type, initial_time, settling_time, out, num_scales, NULL, NULL, NULL, NULL, NULL, period, num_threads);

	writeimage( out, argv[4],  row, col, z);
	//writeimage( dbl_type, "seeds.pgm",  row, col);
	write_classical_image( g, "g.pgm", row, col, z);

	free(g);
	free(img);
	free(out);

	return 0;
}


/* =============================================================== */
void compute_metric(int * img, double * output, int row, int col, int z)
	/* =============================================================== */
{
	int a,b,c,d,e,f, *p;
	long i;
	double hor, vert, zind;
	DBL_TYPE *indbl, *smooth, *q;

	// quick convert to double
	indbl = (DBL_TYPE *)malloc(row*col*z*sizeof(DBL_TYPE));
	q = indbl;
	p = img;
	for (i = 0 ; i < row*col*z ;++i)
		*q++ = *p++;

	// smooth via Gaussian convolution
	smooth = (DBL_TYPE*)malloc(row*col*z*sizeof(DBL_TYPE));
	lsepgaussf_double(indbl, row, col, z ,1.0,smooth);

	dbgprintf("starting computing convolution\n");

	for (i = 0; i< row*col*z;i++)
	{
		a = voisin(i, 1, row, col, z);
		b = voisin(i, 5, row, col, z);
		c = voisin(i, 3, row, col, z);
		d = voisin(i, 7, row, col, z);
		e = voisin(i, 9, row, col, z);
		f = voisin(i, 11, row, col, z);

		//a = -1;
		//      b = -1;
		//      c = -1;
		//      d = -1;
		//e = -1;
		//f = -1;

		//dbgprintf("compute_metric\n");

		if  ((a == -1) || (b == -1))
			hor = 0;
		else
			hor = (smooth[a] -  smooth[b])*(smooth[a] -  smooth[b]); // pow() is very slow

		if  ((c == -1) || (d == -1))
			vert = 0;
		else {
			//vert = pow(img[c] - img[d],2);
			vert = (smooth[c] - smooth[d])*(smooth[c] - smooth[d]);
		}

		if  ((e == -1) || (f == -1))
			zind = 0;
		else {
			zind = (smooth[e] - smooth[f])*(smooth[e] - smooth[f]);
		}

		//dbgprintf("compute_metric output\n");

		// $$ 1/(1+(\nabla I)^2) $$ is good also
		output[i] = 1/(1+sqrt(hor+vert+zind));
		//printf(" %f ", output[i]);
	}

	//dbgprintf("end convolution\n");
}


/* =============================================================== */
long voisin(int indice, int num, int rs, int cs, int zs)
	//retourne l'indice du voisin numero num
	//   6 7 8  9
	//    5 0 1 
	// 11  4 3 2 
	//retourne -1 si le voisin est au l'exterieur de l'image 
	// rs : rowsize : taille des lignes de l'image 
	// cs : colsize : taille des colonnes de l'image 
	// zs : z-index size
	/* =============================================================== */
{
	//dbgprintf("voisin\n");
	int x, y, z;

	x = indice % rs; 
	y = (int)(indice/rs)%cs;
	z = (int)(indice/(rs*cs));

	//dbgprintf("%d %d %d \n", x,y,z);

	switch(num)
	{
	case 1:
		if (x == rs-1) return -1;
		else return x+1+y*rs+z*rs*cs;
		break;
	case 11: // added for 3d
		if (z == 0) return -1;
		else return x+y*rs+(z-1)*rs*cs;
		break;
	case 9: // added for 3d
		if (z == zs-1) return -1;
		else return  x+y*rs+(z+1)*rs*cs;
		break;
	case 8:
		if ((x == rs-1)||(y == 0)) return -1;
		else return x+1+(y-1)*rs+z*rs*cs;
		break;
	case 7:
		if (y == 0) return -1;
		else return x+(y-1)*rs+z*rs*cs;
		break;
	case 6:
		if ((x == 0)||(y == 0)) return -1;
		else return x-1+(y-1)*rs+z*rs*cs;
		break;
	case 5:
		if (x == 0) return -1;
		else return x-1+y*rs+z*rs*cs;
		break;
	case 4:
		if ((x == 0)||(y == cs-1)) return -1;
		else return x-1+(y+1)*rs+z*rs*cs;
		break;
	case 3:
		if (y == cs-1) return -1;
		else return x+(y+1)*rs+z*rs*cs;
		break;
	case 2:
		if ((x == rs-1)||(y == cs-1)) return -1;
		else return x+1+(y+1)*rs+z*rs*cs;
		break;
	}

	// thisa is not possible
	return -2;
}

/* ===========================*/
double max(double x, double y)
	/* ===========================*/
{
	if (x > y) return x;
	return y;
}

void read_3dimage(char * filename, int * img, int row, int col, int z, int header_size)
{
	FILE *fd = NULL;
	fd = fopen(filename,"r");

	char buffer[1000];
	int c;

	c = fread(buffer, sizeof(unsigned char), header_size, fd);

	unsigned char * img2 = (unsigned char *)malloc(sizeof(unsigned char)* row * col * z);

	int ret = fread(img2, sizeof(unsigned char), col * row * z, fd);
	if (ret != row * col * z)
	{
		fprintf(stderr,"fread failed : %d asked ; %d read\n", col * row * z, ret);
	}

	int i;
	for (i = 0; i< col* row * z ; ++i)
	{
		img[i] = (int)(img2[i]);
	}
	fclose(fd);
}

/* ====================================================================== */
void read_image(char * filename, int * img, int * row, int * col)
	/* ====================================================================== */
{
	FILE *fd = NULL;
	fd = fopen(filename,"r");

	char * read;
	int c, d;
	char buffer[1000];
	read = fgets(buffer,1000, fd); // P5: raw byte bw  ; P2: ascii bw  

	if (!read)
	{
		fprintf(stderr, "read_image: fgets returned without reading\n");
	}

	if (buffer[0] != 'P')
	{  
		fprintf(stderr,"read_image : invalid image format\n");
	}

	do
	{ 
		fgets(buffer, 1000, fd);
	}
	while (!isdigit(buffer[0]));

	c = sscanf(buffer, "%d %d \n", &*row, &*col);
	c = fread(&d, sizeof(unsigned char), 4, fd);
	unsigned char * img2 = (unsigned char *)malloc(sizeof(unsigned char)* *row * *col);

	int ret = fread(img2, sizeof(unsigned char), *col * *row , fd);
	if (ret != *row * *col)
	{
		fprintf(stderr,"fread failed : %d asked ; %d read\n", *col * *row, ret);
	}
	int i;
	for (i = 0; i< *col* *row  ;i++)
	{
		img[i] = (int)(img2[i]);
	}
	fclose(fd);
}


/* ====================================================================== */
void read_int_image(int * img, int row, int col)
	/* ====================================================================== */
{ 
	int i; 
	for (i = 0; i< row*col;i++)
	{
		img[i] = 255;
	}
	img[0]=0;   img[9]=0;  img[18]=0;
	img[63]=0;   img[63-9]=0;  img[63-18]=0;

}


/* ==================================== */
void writeimage( DBL_TYPE * image, char * filename, int row, int col, int z)
	/* ==================================== */
{
	FILE * fd;
	int y;
	unsigned char buffer[1];
	fd = fopen(filename,"w+");
	if (fd == NULL)
	{
		printf("erreur fichier\n");
		exit(-1);
	}

	fprintf(fd,"P2 \n");
	fprintf(fd,"%d ", row);
	fprintf(fd,"%d ", col);
	fprintf(fd,"%d \n", z);
	fprintf(fd,"255 \n");

	for( y=0; y<row*col*z;y++)
	{
		//if(y%row == 0)  fprintf(fd, "\n");
		if (image[y] < MYEPSILON) {
			buffer[0] = 0;
			fwrite(buffer , sizeof(unsigned char),  1, fd);
		}
		else {
			if (image[y] > (1-MYEPSILON)) {
				buffer[0] = 255;
				fwrite( buffer , sizeof(unsigned char),  1, fd);
			} else {
				buffer[0] = 128;
				fwrite( buffer , sizeof(unsigned char),  1, fd);
				//	printf( " %d ",(int)(image[y]));
			}
		}
	}
	fflush(fd);
	fclose(fd);
}


/* ==================================== */
void writeimage3d( DBL_TYPE * image, char * filename, int row, int col, int z)
	/* ==================================== */
{
	FILE * fd;
	DBL_TYPE min, max;
	int y;
	unsigned char buffer[1];
	fd = fopen(filename,"w+");
	if (fd == NULL)
	{
		printf("erreur fichier\n");
		exit(-1);
	}

	fprintf(fd,"P2 \n");
	fprintf(fd,"%d ", row);
	fprintf(fd,"%d ", col);
	fprintf(fd,"%d \n", z);
	fprintf(fd,"255 \n");

	min = max = image[0];
	for( y=1; y<row*col*z;y++)
	{
		if (image[y]>max) max=image[y];
		if (image[y]<min) min=image[y];
	}

	dbgprintf("write_classical_image: min=%g, max=%g\n", min, max);
	if (max!=min) {
		for( y=0; y<row*col*z;y++)
		{
			//if(y%row == 0)  fprintf(fd, "\n");
			buffer[0] = (unsigned char)(255*(image[y]-min)/(max-min));
			fwrite( buffer , sizeof(unsigned char),  1, fd);
		}
	} else {
		for( y=0; y<row*col*z;y++)
		{
			//if(y%row == 0)  fprintf(fd, "\n");
			buffer[0] = 0;
			fwrite( buffer , sizeof(unsigned char),  1, fd);
		}
	}


	fflush(fd);
	fclose(fd);
}







/* ==================================== */
/* writes a DOUBLE image, by converting to char first */
void write_classical_image( DBL_TYPE * image, char * filename, int row, int col, int z)
	/* ==================================== */
{
	FILE * fd;
	DBL_TYPE min, max;
	int y;
	unsigned char buffer[1];
	fd = fopen(filename,"w+");
	if (fd == NULL)
	{
		printf("erreur fichier\n");
		exit(-1);
	}

	fprintf(fd,"P2 \n");
	fprintf(fd,"%d ", row);
	fprintf(fd,"%d ", col);
	fprintf(fd,"%d \n", z);
	fprintf(fd,"255 \n");

	min = max = image[0];
	for( y=1; y<row*col*z;y++)
	{
		if (image[y]>max) max=image[y];
		if (image[y]<min) min=image[y];
	}

	dbgprintf("write_classical_image: min=%g, max=%g\n", min, max);

	if (max!=min) {
		for( y=0; y<row*col*z;y++)
		{
			buffer[0] = (unsigned char)(255*(image[y]-min)/(max-min));
			fwrite( buffer , sizeof(unsigned char),  1, fd);
		}
	} else {
		for( y=0; y<row*col*z;y++)
		{
			buffer[0] = 0;
			fwrite( buffer , sizeof(unsigned char),  1, fd);
		}
	}

	fflush(fd);
	fclose(fd);
}
