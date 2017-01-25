#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>

#include <mpi.h>

#include "cuboid.h"
#include "ghosts.h"
#include "interpolate.h"
#include "jacobi.h"
#include "lueqf.h"
#include "mg.h"
#include "restrict.h"
#include "stencil.h"

static inline int min(int a, int b) {
  return (a<b)?a:b;
}

void mg_setup( mg_data **outdata, int maxlevel, int m, int n, int o,
	       int xstart, int xend, int ystart, int yend, int zstart, int zend,
	       int p, int nu1, int nu2, double omega, int size, double* values, 
	       int* xoff, int* yoff, int* zoff, MPI_Comm cart_comm)
{
/* sets multigrid data for each level: 
    - m, n, o
    - x_off, y_off, z_off
    - m_l, n_l, o_l
    - x_ghosts, y_ghosts, z_ghosts
    - left, right, lower, upper, back, front
    - periodic
    - cartesian communicator cart_comm
   allocates space for each level:
    - sbufxy, sbufxz, sbufyz
    - rbufxy, rbufxz, rbufyz
   allocates space for stencil data on level 0
   inits for each level:
    - v, f, r, e
*/
  int level;
  int i, j, k;
  double sum;
  mg_data *data;

  /* Allocating memory for multigrid parameters */
  data = (mg_data*) malloc(maxlevel*sizeof(mg_data));

  for (level=0;level<maxlevel;level++) {
    /* init stencil data */
    if ( level == 0 ) {
      data[0].size = size;
      data[0].values = (double*) malloc( size*sizeof(double) );
      data[0].x_offsets = (int*) malloc( size*sizeof(int) );
      data[0].y_offsets = (int*) malloc( size*sizeof(int) );
      data[0].z_offsets = (int*) malloc( size*sizeof(int) );
      for( i = 0; i < size; i++ ){
	data[0].values[i] = values[i];
	data[0].x_offsets[i] = xoff[i];
	data[0].y_offsets[i] = yoff[i];
	data[0].z_offsets[i] = zoff[i];
      }
      /* setting zero */
      /* WARNING: currently only zeros at (0.0,0.0,0.0) or at (\pi,\pi,\pi) are supported! */
      sum = 0.0;
      for (i=0; i<size; i++) {
	sum += data[0].values[i];
      }
      if (fabs(sum)<1.0e-10) {
	data[0].zero_at_pi3 = false;
      } else {
	data[0].zero_at_pi3 = true;
      }
    }
    if (level<(maxlevel-1) ){
      /* setting zero */
      /* WARNING: currently only zeros at (0.0,0.0,0.0) or at (\pi,\pi,\pi) are supported! */
      data[level].zero_at_pi3 = false;
      data[level+1].values = NULL;
      data[level+1].x_offsets = NULL;
      data[level+1].y_offsets = NULL;
      data[level+1].z_offsets = NULL;
      set_stencil( data, level, 3 );
    }
    
    /* calculate number of ghosts */
    data[level].x_ghosts = data[level].y_ghosts = data[level].z_ghosts = 0;
    for ( i=0; i<data[level].size; i++ ) {
      if ( data[level].x_offsets[i] > data[level].x_ghosts ) data[level].x_ghosts = data[level].x_offsets[i];
      if ( data[level].y_offsets[i] > data[level].y_ghosts ) data[level].y_ghosts = data[level].y_offsets[i];
      if ( data[level].z_offsets[i] > data[level].z_ghosts ) data[level].z_ghosts = data[level].z_offsets[i];
    }
    /* printf("x_ghosts = %d, y_ghosts = %d, z_ghosts = %d\n",data[level].x_ghosts,data[level].y_ghosts,data[level].z_ghosts); */

    /* set periodic flag */
    data[level].periodic = p;

    /* set number of pre and post smoothing steps */
    data[level].nu1 = nu1;
    data[level].nu2 = nu2;

    /* set relaxation coefficient */
    data[level].omega = omega;

    /* set cartesian communicator */
    data[level].cart_comm = cart_comm;

    if (level==0) {
      data[level].m = m + 2*data[level].x_ghosts;
      data[level].n = n + 2*data[level].y_ghosts;
      data[level].o = o + 2*data[level].z_ghosts;
    } else { 
      data[level].m = (data[level-1].m - 2*data[level].x_ghosts)/2 + 2*data[level].x_ghosts;
      data[level].n = (data[level-1].n - 2*data[level].y_ghosts)/2 + 2*data[level].y_ghosts;
      data[level].o = (data[level-1].o - 2*data[level].z_ghosts)/2 + 2*data[level].z_ghosts;
    }

    if (xstart%2==0)
      data[level].x_off = 0;
    else
      data[level].x_off = 1;

    if (ystart%2==0)
      data[level].y_off = 0;
    else
      data[level].y_off = 1;

    if (zstart%2==0)
      data[level].z_off = 0;
    else
      data[level].z_off = 1;

    if (xend>=xstart)
      data[level].m_l = xend - xstart + 1 + 2*data[level].x_ghosts;
    else
      data[level].m_l = 0;
    if (yend>=ystart)
      data[level].n_l = yend - ystart + 1 + 2*data[level].y_ghosts;
    else
      data[level].n_l = 0;
    if (zend>=zstart)
      data[level].o_l = zend - zstart + 1 + 2*data[level].z_ghosts;
    else
      data[level].o_l = 0;

    if (level == 0) {
      /* MPI variables */
      int cart_dims[3], cart_periods[3], cart_coords[3];

      /* Getting processes coordinates */
      MPI_Cart_get(cart_comm,3,cart_dims,cart_periods,cart_coords);
      
      /* Calculating neighbours */
      MPI_Cart_shift(cart_comm, 0, 1, &data[0].left , &data[0].right);
      MPI_Cart_shift(cart_comm, 1, 1, &data[0].lower, &data[0].upper);
      MPI_Cart_shift(cart_comm, 2, 1, &data[0].back,  &data[0].front);
    } else {
      /* Temporary buffer */
      int tmpbuf;

      /* MPI variables */
      int myid;
      MPI_Request reqs[2];
      MPI_Status stats[2];
      
      /* Getting local rank */
      MPI_Comm_rank(cart_comm, &myid);

      if( data[level].periodic ){
        /* Initializing neighbours */
        data[level].left  = myid;
        data[level].right = myid;
        data[level].lower = myid;
        data[level].upper = myid;
        data[level].back  = myid;
        data[level].front = myid;
      }
      else{
        /* Initializing neighbours */
        data[level].left  = MPI_PROC_NULL;
        data[level].right = MPI_PROC_NULL;
        data[level].lower = MPI_PROC_NULL;
        data[level].upper = MPI_PROC_NULL;
        data[level].back  = MPI_PROC_NULL;
        data[level].front = MPI_PROC_NULL;
      }
      
      if (data[level].m_l>0) {
	/* Sending data to left */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].left) && (MPI_PROC_NULL != data[level-1].left) )
	  MPI_Isend((void *) &myid,1,MPI_INT,data[level-1].left,0,cart_comm,
		    &reqs[0]);
	if ( (myid!=data[level-1].right) && (MPI_PROC_NULL != data[level-1].right) )
	  MPI_Irecv((void *) &data[level].right,1,MPI_INT,data[level-1].right,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
	
	/* Sending data to right */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].right) && (MPI_PROC_NULL != data[level-1].right) )
	  MPI_Isend((void *) &myid,1,MPI_INT,data[level-1].right,0,cart_comm,
		    &reqs[0]);
	if ( (myid!=data[level-1].left) && (MPI_PROC_NULL != data[level-1].left) )
	  MPI_Irecv((void *) &data[level].left,1,MPI_INT,data[level-1].left,0,
		    cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
      } else {
	/* Sending data to left */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].left) && (MPI_PROC_NULL != data[level-1].left) )
	  MPI_Isend((void *) &data[level-1].right,1,MPI_INT,data[level-1].left,
		    0,cart_comm,&reqs[0]);
	if ( (myid!=data[level-1].right) && (MPI_PROC_NULL != data[level-1].right) )
	  MPI_Irecv((void *) &tmpbuf,1,MPI_INT,data[level-1].right,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
	
	/* Sending data to right */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].right) && (MPI_PROC_NULL != data[level-1].right) )
	  MPI_Isend((void *) &data[level-1].left,1,MPI_INT,data[level-1].right,
		    0,cart_comm,&reqs[0]);
	if ( (myid!=data[level-1].left) && (MPI_PROC_NULL != data[level-1].left) )
	  MPI_Irecv((void *) &tmpbuf,1,MPI_INT,data[level-1].left,0,
		    cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
      }

      if (data[level].n_l>0) {
	/* Sending data downwards */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].lower) && (MPI_PROC_NULL != data[level-1].lower) )
	  MPI_Isend((void *) &myid,1,MPI_INT,data[level-1].lower,0,cart_comm,
		    &reqs[0]);
	if ( (myid!=data[level-1].upper) && (MPI_PROC_NULL != data[level-1].upper) )
	  MPI_Irecv((void *) &data[level].upper,1,MPI_INT,data[level-1].upper,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
	
	/* Sending data upwards */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].upper) && (MPI_PROC_NULL != data[level-1].upper) )
	  MPI_Isend((void *) &myid,1,MPI_INT,data[level-1].upper,0,cart_comm,
		    &reqs[0]);
	if ( (myid!=data[level-1].lower) && (MPI_PROC_NULL != data[level-1].lower) )
	  MPI_Irecv((void *) &data[level].lower,1,MPI_INT,data[level-1].lower,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
      } else {
	/* Sending data downwards */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].lower) && (MPI_PROC_NULL != data[level-1].lower) )
	  MPI_Isend((void *) &data[level-1].upper,1,MPI_INT,
		    data[level-1].lower,0,cart_comm,&reqs[0]);
	if ( (myid!=data[level-1].upper) && (MPI_PROC_NULL != data[level-1].upper) )
	  MPI_Irecv((void *) &tmpbuf,1,MPI_INT,data[level-1].upper,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
	
	/* Sending data upwards */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].upper) && (MPI_PROC_NULL != data[level-1].upper) )
	  MPI_Isend((void *) &data[level-1].lower,1,MPI_INT,
		    data[level-1].upper,0,cart_comm,&reqs[0]);
	if ( (myid!=data[level-1].lower) && (MPI_PROC_NULL != data[level-1].lower) )
	  MPI_Irecv((void *) &tmpbuf,1,MPI_INT,data[level-1].lower,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
      }

      if (data[level].o_l>0) {
	/* Sending data to back */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].back) && (MPI_PROC_NULL != data[level-1].back) )
	  MPI_Isend((void *) &myid,1,MPI_INT,data[level-1].back,0,cart_comm,
		    &reqs[0]);
	if ( (myid!=data[level-1].front) && (MPI_PROC_NULL != data[level-1].front))
	  MPI_Irecv((void *) &data[level].front,1,MPI_INT,data[level-1].front,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
	
	/* Sending data to front */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].front) && (MPI_PROC_NULL != data[level-1].front) )
	  MPI_Isend((void *) &myid,1,MPI_INT,data[level-1].front,0,cart_comm,
		    &reqs[0]);
	if ( (myid!=data[level-1].back) && (MPI_PROC_NULL != data[level-1].back) )
	  MPI_Irecv((void *) &data[level].back,1,MPI_INT,data[level-1].back,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
      } else {
	/* Sending data to back */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].back) && (MPI_PROC_NULL != data[level-1].back))
	  MPI_Isend((void *) &data[level-1].front,1,MPI_INT,
		    data[level-1].back,0,cart_comm,&reqs[0]);
	if ( (myid!=data[level-1].front) && (MPI_PROC_NULL != data[level-1].front) )
	  MPI_Irecv((void *) &tmpbuf,1,MPI_INT,data[level-1].front,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
	
	/* Sending data to front */
	reqs[0] = MPI_REQUEST_NULL; reqs[1] = MPI_REQUEST_NULL;
	if ( (myid!=data[level-1].front) && (MPI_PROC_NULL != data[level-1].front))
	  MPI_Isend((void *) &data[level-1].back,1,MPI_INT,
		    data[level-1].front,0,cart_comm,&reqs[0]);
	if ( (myid!=data[level-1].back) && (MPI_PROC_NULL != data[level-1].back) )
	  MPI_Irecv((void *) &tmpbuf,1,MPI_INT,data[level-1].back,
		    0,cart_comm,&reqs[1]);
	MPI_Waitall(2,reqs,stats);
      }

      if (data[level].m_l<=0 || data[level].n_l<=0 || data[level].o_l<=0) {
	data[level].left  = myid;
	data[level].right = myid;
	data[level].lower = myid;
	data[level].upper = myid;
	data[level].back  = myid;
	data[level].front = myid;
      }
    }

    if (data[level].m_l>0 && data[level].n_l>0 && data[level].o_l>0) {
      data[level].v = cuboid_alloc(data[level].m_l,data[level].n_l,data[level].o_l);
      data[level].f = cuboid_alloc(data[level].m_l,data[level].n_l,data[level].o_l);
      data[level].r = cuboid_alloc(data[level].m_l,data[level].n_l,data[level].o_l);
      data[level].e = cuboid_alloc(data[level].m_l,data[level].n_l,data[level].o_l);
      data[level].tmp = cuboid_alloc(data[level].m_l,data[level].n_l,data[level].o_l);
      data[level].sbufxy = cuboid_alloc(data[level].m_l,data[level].n_l,data[level].z_ghosts);
      data[level].sbufxz = cuboid_alloc(data[level].m_l,data[level].y_ghosts,data[level].o_l);
      data[level].sbufyz = cuboid_alloc(data[level].x_ghosts,data[level].n_l,data[level].o_l);
      data[level].rbufxy = cuboid_alloc(data[level].m_l,data[level].n_l,data[level].z_ghosts);
      data[level].rbufxz = cuboid_alloc(data[level].m_l,data[level].y_ghosts,data[level].o_l);
      data[level].rbufyz = cuboid_alloc(data[level].x_ghosts,data[level].n_l,data[level].o_l);

      for (i=0;i<data[level].m_l;i++) {
	for (j=0;j<data[level].n_l;j++) {
	  for (k=0;k<data[level].o_l;k++) {
	    data[level].v[i][j][k] = 0.0;
	    data[level].f[i][j][k] = 0.0;
	    data[level].r[i][j][k] = 0.0;
	    data[level].e[i][j][k] = 0.0;
	    data[level].tmp[i][j][k] = 0.0;
	  }
	}
      }
    } else {
      data[level].v = NULL;
      data[level].f = NULL;
      data[level].r = NULL;
      data[level].e = NULL;
      data[level].tmp = NULL;
      data[level].sbufxy = NULL;
      data[level].sbufxz = NULL;
      data[level].sbufyz = NULL;
      data[level].rbufxy = NULL;
      data[level].rbufxz = NULL;
      data[level].rbufyz = NULL;
    }

#ifdef DEBUG
    /* Debugging output */
    {
      int cart_dims[3], cart_periods[3], cart_coords[3];

      /* Getting processes coordinates */
      MPI_Cart_get(cart_comm,3,cart_dims,cart_periods,cart_coords);
      
      printf("[%3d,%3d,%3d]: level = %d, xstart = %3d, xend   = %3d, m_l = %3d\n",
	     cart_coords[0],cart_coords[1],cart_coords[2],level,xstart,xend,
	     data[level].m_l);
      printf("[%3d,%3d,%3d]: level = %d, ystart = %3d, yend   = %3d, n_l = %3d\n",
	     cart_coords[0],cart_coords[1],cart_coords[2],level,ystart,yend,
	     data[level].n_l);
      printf("[%3d,%3d,%3d]: level = %d, zstart = %3d, zend   = %3d, o_l = %3d\n",
	     cart_coords[0],cart_coords[1],cart_coords[2],level,zstart,zend,
	     data[level].o_l);
      printf("[%3d,%3d,%3d]: level = %d, x_off  = %1d, y_off = %1d, z_off = %1d\n",
	     cart_coords[0],cart_coords[1],cart_coords[2],level,
	     data[level].x_off,data[level].y_off,data[level].z_off);
      printf("[%3d,%3d,%3d]: level = %d, left  = %3d, right = %3d\n",
	     cart_coords[0],cart_coords[1],cart_coords[2],level,
	     data[level].left,data[level].right);
      printf("[%3d,%3d,%3d]: level = %d, lower = %3d, upper = %3d\n",
	     cart_coords[0],cart_coords[1],cart_coords[2],level,
	     data[level].lower,data[level].upper);
      printf("[%3d,%3d,%3d]: level = %d, back  = %3d, front = %3d\n",
	     cart_coords[0],cart_coords[1],cart_coords[2],level,
	     data[level].back,data[level].front);
    }
#endif

    /* Calculating data distribution on next level */
    xstart = xstart/2;
    xend = (xend+1)/2 - 1;
    ystart = ystart/2;
    yend = (yend+1)/2 - 1;
    zstart = zstart/2;
    zend = (zend+1)/2 - 1;
  }

  *outdata = data;

  return;
}

void mg_init(double ***v, double ***f, mg_data *data )
{
/* inits first level (0):
    - v
    - f
*/
  int i, j, k;

  /* init initial guess and rhs */
  for (i=0;i<data[0].m_l-2*data[0].x_ghosts;i++) {
    for (j=0;j<data[0].n_l-2*data[0].y_ghosts;j++) {
      for (k=0;k<data[0].o_l-2*data[0].z_ghosts;k++) {
	data[0].v[i+data[0].x_ghosts][j+data[0].y_ghosts][k+data[0].z_ghosts] = v[i][j][k];
	data[0].f[i+data[0].x_ghosts][j+data[0].y_ghosts][k+data[0].z_ghosts] = f[i][j][k];
      }
    }
  }
  
  update_ghosts( data[0].v, data, 0 );
  update_ghosts( data[0].f, data, 0 );

/*       /\* Debugging output *\/ */
/* #ifdef DEBUG */
/*   { */
/*     int cart_dims[3], cart_periods[3], cart_coords[3]; */
    
/*     /\* Getting processes coordinates *\/ */
/*     MPI_Cart_get(data[0].cart_comm,3,cart_dims,cart_periods,cart_coords); */
    
/*     int level; */
/*     for( level = 0; level < maxlevel-1; level++ ){ */
/*       printf("[%3d,%3d,%3d]: level = %d, size = %3d\n", */
/* 	     cart_coords[0],cart_coords[1],cart_coords[2],level,data[level].size); */
/*       printf("[%3d,%3d,%3d]: level = %d, values:\n", */
/* 	     cart_coords[0],cart_coords[1],cart_coords[2],level); */
/*       for( i = 0; i < data[level].size; i++ ) */
/*         printf( "%.10f ", data[level].values[i] ); */
/*       printf("\n[%3d,%3d,%3d]: level = %d, x_offsets:\n", */
/* 	     cart_coords[0],cart_coords[1],cart_coords[2],level); */
/*       for( i = 0; i < data[level].size; i++ ) */
/*         printf( "%3d ", data[level].x_offsets[i] ); */
/*       printf("\n[%3d,%3d,%3d]: level = %d, y_offsets:\n", */
/* 	     cart_coords[0],cart_coords[1],cart_coords[2],level); */
/*       for( i = 0; i < data[level].size; i++ ) */
/*         printf( "%3d ", data[level].y_offsets[i] ); */
/*       printf("\n[%3d,%3d,%3d]: level = %d, z_offsets:\n", */
/* 	     cart_coords[0],cart_coords[1],cart_coords[2],level); */
/*       for( i = 0; i < data[level].size; i++ ) */
/*         printf( "%3d ", data[level].z_offsets[i] ); */
/*       printf( "\n" ); */
/*     } */
/*   } */
/* #endif */
  
  return;
}

void mg_result(double ***v, mg_data *data)
{
  int i, j, k;
 
  for (i=0;i<data[0].m_l-2*data[0].x_ghosts;i++) {
    for (j=0;j<data[0].n_l-2*data[0].y_ghosts;j++) {
      for (k=0;k<data[0].o_l-2*data[0].z_ghosts;k++) {
	v[i][j][k] = data[0].v[i+data[0].x_ghosts][j+data[0].y_ghosts][k+data[0].z_ghosts];
      }
    }
  }

  return;
}

void mg_free(mg_data *data, int maxlevel)
{
  int level;

  for (level=0;level<maxlevel;level++) {
    if (data[level].v != NULL)
      cuboid_free(data[level].v,data[level].m_l,data[level].n_l,data[level].o_l);
    if (data[level].f != NULL)
      cuboid_free(data[level].f,data[level].m_l,data[level].n_l,data[level].o_l);
    if (data[level].r != NULL)
      cuboid_free(data[level].r,data[level].m_l,data[level].n_l,data[level].o_l);
    if (data[level].e != NULL)
      cuboid_free(data[level].e,data[level].m_l,data[level].n_l,data[level].o_l);
    if (data[level].tmp != NULL)
      cuboid_free(data[level].tmp,data[level].m_l,data[level].n_l,data[level].o_l);
    if (data[level].sbufxy != NULL)
      cuboid_free(data[level].sbufxy,data[level].m_l,data[level].n_l,data[level].z_ghosts);
    if (data[level].sbufxz != NULL)
      cuboid_free(data[level].sbufxz,data[level].m_l,data[level].y_ghosts,data[level].o_l);
    if (data[level].sbufyz != NULL)
      cuboid_free(data[level].sbufyz,data[level].x_ghosts,data[level].n_l,data[level].o_l);
    if (data[level].rbufxy != NULL)
      cuboid_free(data[level].rbufxy,data[level].m_l,data[level].n_l,data[level].z_ghosts);
    if (data[level].rbufxz != NULL)
      cuboid_free(data[level].rbufxz,data[level].m_l,data[level].y_ghosts,data[level].o_l);
    if (data[level].rbufyz != NULL)
      cuboid_free(data[level].rbufyz,data[level].x_ghosts,data[level].n_l,data[level].o_l);
    if( data[level].values != NULL )
      free( data[level].values );
    if( data[level].x_offsets != NULL )
      free( data[level].x_offsets );
    if( data[level].y_offsets != NULL )
      free( data[level].y_offsets );
    if( data[level].z_offsets != NULL )
      free( data[level].z_offsets );
  }

  free(data);

  return;
}

double mg_vcycle( mg_data *data, int level, int maxlevel )
{
  int i, j, k;
  double res = 0.0;

  /* Clear old v_2h */
  for (i=0;i<data[level+1].m_l;i++) {
    for (j=0;j<data[level+1].n_l;j++) {
      for (k=0;k<data[level+1].o_l;k++) {
	data[level+1].v[i][j][k] = 0.0;
      }
    }
  }

  /* v_h <- smooth(v_h,f_h,nu1) */
  jacobi( data[level].v, data[level].f, data[level].tmp, data, level, data[level].nu1 );

  /* r_h <- f_h - L * v_h */
  lueqf_res( data[level].v, data[level].f, data[level].r, data, level );
  update_ghosts( data[level].r, data, level );

  /* f_2h <- I_h^2h(r_h) */
  if (data[level+1].m_l>0 && data[level+1].n_l>0 && data[level+1].o_l>0) {
    /*restrict_inj(data[level].r,data[level+1].f, data, level );*/
    restrict_fw(data[level].r,data[level+1].f,data, level); 
    update_ghosts( data[level+1].f, data, level+1 );
  }

  if (level<(maxlevel-1)) {
    if (data[level+1].m_l>0 && data[level+1].n_l>0 && data[level+1].o_l>0) {
      if (level<(maxlevel-2)) {
	/* v_2h <- L^(-1) * f_2h */
	mg_vcycle( data, level + 1, maxlevel );
      } else {
	jacobi( data[level+1].v, data[level+1].f, data[level+1].tmp, data, level+1, data[level].nu1+data[level].nu2 );
      }
    }
  }

  /* e_h <- I_2h^h(v_2h) */
  if (data[level+1].m_l>0 && data[level+1].n_l>0 && data[level+1].o_l>0) {
    interpolate_prepare( data[level].e, data[level+1].v, data, level );
  }

  update_ghosts( data[level].e, data, level );
  interpolate_finish( data[level].e, data, level );
  update_ghosts( data[level].e, data, level );

  /* v_h <- v_h + e_h */
  for (i=0;i<data[level].m_l;i++) {
    for (j=0;j<data[level].n_l;j++) {
      for (k=0;k<data[level].o_l;k++) {
	data[level].v[i][j][k] = data[level].v[i][j][k] + data[level].e[i][j][k];
      }
    }
  }

  /* v_h <- smooth(v_h,f_h,nu2) */
  res = jacobi( data[level].v, data[level].f, data[level].tmp, data, level, data[level].nu2 );

  return(res);
}

double mg(double ***u, double ***f, int maxiter, double tol, int m, int n, int o,
	  int xstart, int xend, int ystart, int yend, int zstart, int zend,
	  int p, int nu1, int nu2, double omega, int size, double* values,
	  int* xoff, int* yoff, int* zoff, MPI_Comm cart_comm, int verbose)
{

  /* MPI variables */
  int myid;

  /* Number of levels */
  int maxlevel;

  /* Pointer to solver data */
  mg_data *data;

  /* Variables for data distribution of the problem */
  int m_l, n_l, o_l;

  /* Loop counters */
  int iter;

  /* Residuals */
  double initres_l, initres, res_l, res;

  /* For time measurement */
  struct timeval start, stop;
  double elapsed;

  /* Getting local rank */
  MPI_Comm_rank(cart_comm, &myid);

  m_l = xend - xstart + 1;
  n_l = yend - ystart + 1;
  o_l = zend - zstart + 1;

  /* Calculating number of levels and checking input size */
  if (p) {
    maxlevel = (int) floor(log((double) (min(min(m,n),o)+1))/log(2.0)) + 1;
    if (min(min(m,n),o)!=(int) pow(2.0,(double) (maxlevel-1))) {
      printf("min(min(m,n),o) = %d != %d = 2^%d\n",min(min(m,n),o),(int) pow(2.0,(double) (maxlevel-1)),maxlevel-1);
      exit(123);
    }
  } else {
    maxlevel = (int) floor(log((double) min(min(m,n),o))/log(2.0)) + 1;
    if (min(min(m,n),o)!=(int) pow(2.0,(double) (maxlevel)) - 1) {
      printf("min(min(m,n),o) = %d != %d = 2^%d - 1\n",min(min(m,n),o),(int) pow(2.0,(double) (maxlevel)) - 1,maxlevel);
      exit(123);
    }
  }

  /* Initializing multigrid solver */
  mg_setup(&data, maxlevel, m, n, o, xstart, xend, ystart, yend, zstart, zend,
  	   p, nu1, nu2, omega, size, values, xoff, yoff, zoff, cart_comm);
  mg_init( u, f, data );
  initres_l = lueqf_res( data[0].v, data[0].f, data[0].tmp, data, 0 );
  MPI_Allreduce(&initres_l,&initres,1,MPI_DOUBLE,MPI_SUM,cart_comm);
  initres = sqrt(initres);
  
  /* Output */
#ifdef FCS_ENABLE_INFO
  if (myid == 0 && verbose > 0) {
    printf("*************** MG with initres = %e ***************\n",initres);
  }
  if (myid == 0 && verbose > 1) {
    printf("\n");
    printf("+-------+-----------+----------+\n");
    printf("| iter. | abs. res. | time     |\n");
    printf("+-------+-----------+----------+\n");
  }
#endif

  for (iter=0;iter<maxiter;iter++) {
    /* Start date */
    gettimeofday(&start,0);

    res_l = mg_vcycle( data, 0, maxlevel );
    MPI_Allreduce(&res_l,&res,1,MPI_DOUBLE,MPI_SUM,cart_comm);
    res = sqrt(res);

    /* Stop date */
    gettimeofday(&stop,0);
    elapsed = (double) (stop.tv_sec - start.tv_sec) + ((double) (stop.tv_usec - start.tv_usec))/1000000.0;

    /* Output */
#ifdef FCS_ENABLE_INFO
    if ((myid == 0 && verbose > 1) || (myid == 0 && verbose > 0 && iter+1 == maxiter)){
      printf("| %5d | %9.3e | %8.6f |\n", iter, initres == 0.0 ? 0.0 : res, elapsed);
    }
#endif

    /* Relative residual */
    /* if (res/initres<tol) */
    /* Absolute residual */
    if (res<tol)
      break;
  }

  /* Output */
#ifdef FCS_ENABLE_INFO
  if (myid == 0 && verbose > 1) {
    printf("+-------+-----------+----------+\n");
    printf("\n");
  }
#endif

  mg_result(u,data);
  mg_free(data,maxlevel);

  return(res);
}
