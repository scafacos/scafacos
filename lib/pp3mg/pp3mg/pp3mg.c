/*
 * pp3mg.c
 *
 * This file contains the functions
 * "pp3mg_init"
 * "pp3mg"
 * "pp3mg_free"
 *
 * Based on Fortran code by Matthias Bolten.
 *
 * Authors: Matthias Bolten, Stephanie Friedhoff
 *
 * Copyright 2009, 2010, 2011, 2012 Matthias Bolten, Stephanie Friedhoff
 * All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include "pp3mg.h"

#include "cuboid.h"
#include "mg.h"
#include "ghosted_grid.h"
#include "interpolation.h"
#include "particle.h"

#define POLYNOMIAL10
#define SIXTHORDER

/*
 *
 * Constants
 *
 */

/* Pi */
#define PP3MG_PI 3.1415926535897932

/*
 *
 * Function definitions
 *
 */

/* Maximum */
#define MAX( a, b ) ((a>b)?a:b)
/* Minimum */
#define MIN( a, b ) ((a<b)?a:b)
/* Square */
#define SQUARE( x ) ((x)*(x))

/* ----------------------------------------------------------------------- */
	
void pp3mg_init( double x_in, double y_in, double z_in, int m_in, int n_in,
		 int o_in, int ghosts_in, int degree_in, int max_particles_in, 
		 int maxiter_in, double tol_in, 
		 enum CHARGE_DISTRIBUTION distribution_in, 
		 enum DISCRETIZATION discretization_in, MPI_Comm mpi_comm,
		 pp3mg_data* data, pp3mg_parameters* params)
{
		
  /* Variables for MPI */
  int mpi_size;
  int mpi_dims[3], mpi_periods[3];
  int mpi_rank;
  int mpi_coords[3];

  /* Copy input data to global variables */
  params->x = x_in;
  params->y = y_in;
  params->z = z_in;
  params->m = m_in;
  params->n = n_in;
  params->o = o_in;
  params->ghosts = ghosts_in;
  params->degree = degree_in;
  params->maxiter = maxiter_in;
  params->tol = tol_in;
  params->distribution = distribution_in;
  params->discretization = discretization_in;
	
  /* Initialize cartesian process grid */
  mpi_dims[0]    = 0;
  mpi_dims[1]    = 0;
  mpi_dims[2]    = 0;
  mpi_periods[0] = 1;
  mpi_periods[1] = 1;	
  mpi_periods[2] = 1;
  MPI_Comm_size( mpi_comm, &mpi_size );
  MPI_Dims_create( mpi_size, 3, mpi_dims );
  MPI_Cart_create( mpi_comm, 3, mpi_dims, mpi_periods, 1, 
		   &(params->mpi_comm_cart) );
  MPI_Comm_rank( params->mpi_comm_cart, &mpi_rank );
  MPI_Cart_coords( params->mpi_comm_cart, mpi_rank, 3, mpi_coords );

  /* Initialize parameters for particle grid */
  params->m_start = ( params->m/mpi_dims[0] ) * mpi_coords[0] 
    + MIN( mpi_coords[0], ( params->m % mpi_dims[0] ) );
  params->m_end   = ( params->m/mpi_dims[0]) * ( mpi_coords[0] + 1 )
    + MIN( mpi_coords[0]+1, ( params->m % mpi_dims[0] ) ) 
    - 1;
  params->n_start = ( params->n/mpi_dims[1] ) * mpi_coords[1]						
    + MIN( mpi_coords[1], ( params->n % mpi_dims[1] ) );
  params->n_end   = ( params->n/mpi_dims[1] ) * ( mpi_coords[1] + 1 )
    + MIN( mpi_coords[1]+1, ( params->n % mpi_dims[1] ) ) 
    - 1;
  params->o_start = ( params->o/mpi_dims[2]) * mpi_coords[2]
    + MIN( mpi_coords[2], ( params->o % mpi_dims[2] ) );
  params->o_end   = ( params->o/mpi_dims[2] ) * ( mpi_coords[2] + 1 )
    + MIN( mpi_coords[2]+1, ( params->o % mpi_dims[2] ) ) 
    - 1;
  params->x_start = 1.0 / mpi_dims[0] * mpi_coords[0] * params->x;
  params->x_end   = 1.0 / mpi_dims[0] * (mpi_coords[0]+1) * params->x;
  params->y_start = 1.0 / mpi_dims[1] * mpi_coords[1] * params->y;
  params->y_end   = 1.0 / mpi_dims[1] * (mpi_coords[1]+1) * params->y;
  params->z_start = 1.0 / mpi_dims[2] * mpi_coords[2] * params->z;
  params->z_end   = 1.0 / mpi_dims[2] * (mpi_coords[2]+1) * params->z;
  params->hx      = params->x / (double) params->m;
  params->hy      = params->y / (double) params->n;
  params->hz      = params->z / (double) params->o;
  params->radius = params->hx;
  if (params->hy < params->radius) params->radius = params->hy;
  if (params->hz < params->radius) params->radius = params->hz;
  params->radius  *= (double) params->ghosts;
  params->radius_2 = params->radius * params->radius;
  params->radius_3 = params->radius_2 * params->radius;
  params->radius_4 = params->radius_2 * params->radius_2;
  params->radius_6 = params->radius_4 * params->radius_2;
  params->radius_7 = params->radius_6 * params->radius;

  /* Allocate storage for particles */
  data->max_particles = max_particles_in;
  data->particles = (pp3mg_particle*) malloc( data->max_particles*sizeof(pp3mg_particle) );
  assert( data->particles != NULL);

  /* -------------------------------------------------------------------
   *
   * Allocate storage for linked list and grids
   *
   * ------------------------------------------------------------------- */
  
  /* Allocate storage for linked list */
  data->hoc = (int***) malloc( (params->m_end-params->m_start+2*(params->ghosts)+1)
			       *sizeof(int**) );
  assert( data->hoc != NULL );
  for(int i = 0; i <= (params->m_end-params->m_start+2*(params->ghosts)); i++ ){
    data->hoc[i] = (int**) malloc( (params->n_end-params->n_start+2*(params->ghosts)+1)
				   *sizeof(int*) );
    assert( data->hoc[i] != NULL );
    for( int j = 0; j <= (params->n_end-params->n_start+2*(params->ghosts)); j++ ){
      data->hoc[i][j] = (int*) malloc( (params->o_end-params->o_start+2*(params->ghosts)+1)
				       *sizeof(int) );
      assert( data->hoc[i][j] != NULL );
    }
  }
  data->ll = (int*) malloc( data->max_particles*sizeof(int) );
  assert( data->ll != NULL );
  
  /* Allocate grids */
  data->f = cuboid_alloc( params->m_end-params->m_start+1, params->n_end-params->n_start+1, 
			  params->o_end-params->o_start+1 );
  
  data->f_ghosted = cuboid_alloc( params->m_end-params->m_start+2*(params->ghosts)+1, 
				  params->n_end-params->n_start+2*(params->ghosts)+1, 
				  params->o_end-params->o_start+2*(params->ghosts)+1 );
  
  data->u = cuboid_alloc( params->m_end-params->m_start+1, params->n_end-params->n_start+1, 
			  params->o_end-params->o_start+1 );
  
  data->u_ghosted = cuboid_alloc( params->m_end-params->m_start+2*(params->ghosts)+1, 
				  params->n_end-params->n_start+2*(params->ghosts)+1, 
				  params->o_end-params->o_start+2*(params->ghosts)+1 );

#ifdef DEBUG
  printf("DEBUG (rank=%d): m = %d, n = %d, n = %d\n",
	 mpi_rank,params->m,params->n,params->o);
  printf("DEBUG (rank=%d): mpi_coords = (%d,%d,%d)\n",
	 mpi_rank,mpi_coords[0],mpi_coords[1],mpi_coords[2]);
  printf("DEBUG (rank=%d): mpi_dims = (%d,%d,%d)\n",
	 mpi_rank,mpi_dims[0],mpi_dims[1],mpi_dims[2]);
  printf("DEBUG (rank=%d): x = %f, hx = %f, x_start = %f x_end = %f\n",
	 mpi_rank,params->x,params->hx,params->x_start,params->x_end);
  printf("DEBUG (rank=%d):y = %f, hy = %f, y_start = %f, y_end = %f\n",
	 mpi_rank,params->y,params->hy,params->y_start,params->y_end);
  printf("DEBUG (rank=%d): z = %f, hz = %f, z_start = %f, z_end = %f\n",
	 mpi_rank,params->z,params->hz,params->z_start,params->z_end);
  printf("DEBUG (rank=%d): radius = %f\n",
	 mpi_rank,params->radius);
  printf("DEBUG (rank=%d): m = %d, m_start = %d, m_end = %d\n",
	 mpi_rank,params->m,params->m_start,params->m_end);
  printf("DEBUG (rank=%d): n = %d, n_start = %d, n_end = %d\n",
	 mpi_rank,params->n,params->n_start,params->n_end);
  printf("DEBUG (rank=%d): o = %d, o_start = %d, o_end = %d\n",
	 mpi_rank,params->o,params->o_start,params->o_end);
#endif  /* DEBUG */
}

/* ----------------------------------------------------------------------- */
	
void pp3mg( double* x, double* y, double* z, double* q, double* e, 
	    double* fx, double* fy, double* fz, int n_local_particles_in,
	    pp3mg_data* data, pp3mg_parameters* params )
{
	
  /* Variables for loops */
  int p, p1, p2;
  int i, j, k;
  int ii, jj, kk;
  int i_start, j_start, k_start;
  int i_end, j_end, k_end;
	
  /* Current cell */
  int cell[3];
	
  /* Distance */
  double r;
  /* Components of distance */
  double rx, ry, rz;
	
  /* Value for correction of energy and forces */
  double val;
	
  /* Parameters for multigrid solver */
  int nu1, nu2;
  double omega;
  int size;
  double* values;
  int* xoff;
  int* yoff;
  int* zoff;
  double ret;
	
  /* Variables for MPI */
  int mpi_periods[3];
  int mpi_coords[3];
  int mpi_dims[3];

  /* -------------------------------------------------------------------
   *
   * Copy input data and initialize variables
   *
   * ------------------------------------------------------------------- */

  /* Copy particles to methods storage */
  if (n_local_particles_in > data->max_particles) {
    data->particles = (pp3mg_particle*) realloc( data->particles, n_local_particles_in*sizeof(pp3mg_particle) );
    data->ll = (int*) realloc( data->ll, n_local_particles_in*sizeof(int));
    data->max_particles = n_local_particles_in;

    if (data->particles == NULL || data->ll == NULL)
	{
	  printf("Realloc failed!");
	  exit(1);
	}
  }

  data->n_local_particles = n_local_particles_in;
  for( p = 0; p < data->n_local_particles; p++ ){
    data->particles[p].x  = x[p];
    data->particles[p].y  = y[p];
    data->particles[p].z  = z[p];
    data->particles[p].q  = q[p];
    data->particles[p].fx = fx[p];
    data->particles[p].fy = fy[p];
    data->particles[p].fz = fz[p];
  }

#ifdef DEBUG
  {
    int mpi_rank;

    MPI_Comm_rank( params->mpi_comm_cart, &mpi_rank );
    printf("DEBUG (rank=%d): n_local_particles = %d\n",
	   mpi_rank,data->n_local_particles);
  }
#endif

  /* Counting particles */
  MPI_Allreduce( (void*) &(data->n_local_particles), (void*) &(data->n_particles), 
      1, MPI_INT, MPI_SUM, params->mpi_comm_cart );			

  mpi_dims[0]    = 0;
  mpi_dims[1]    = 0;
  mpi_dims[2]    = 0;
  MPI_Cart_get( params->mpi_comm_cart, 3, mpi_dims, mpi_periods, mpi_coords );

  /* -------------------------------------------------------------------
   *
   * Initializing linked list and copying ghosts
   *
   * ------------------------------------------------------------------- */
  /* Buildung linked list */
  for( i = 0; i <= (params->m_end-params->m_start+2*params->ghosts); i++ )
    for( j = 0; j <= (params->n_end-params->n_start+2*params->ghosts); j++ )
      for( k = 0; k <= (params->o_end-params->o_start+2*params->ghosts); k++ )
	data->hoc[i][j][k] = -1;

  for( i = 0; i < data->max_particles; i++ )
    data->ll[i] = -1;

  update_linked_list( data->ll, data->hoc, data->particles, 0, data->n_local_particles, 
		      params->x, params->y, params->z, params->m, params->n, params->o, 
		      params->m_start-params->ghosts, params->m_end+params->ghosts,
		      params->n_start-params->ghosts, params->n_end+params->ghosts,
		      params->o_start-params->ghosts,params->o_end+params->ghosts );

  data->n_stored_particles = data->n_local_particles;
  
  update_particle_ghosts( &(data->ll), data->hoc, &(data->particles), params->x, params->y, params->z,
			  params->m, params->n, params->o, params->ghosts,
			  &(data->n_stored_particles), &(data->max_particles),
			  params->m_start, params->m_end, params->n_start,
			  params->n_end, params->o_start, params->o_end,
			  params->mpi_comm_cart );

  /* -------------------------------------------------------------------
   *
   * Assigning charges to the grid
   *
   * ------------------------------------------------------------------- */

  /* Initializing right hand side */
  for( i = 0; i <= (params->m_end-params->m_start); i++ )
    for( j = 0; j <= (params->n_end-params->n_start); j++ )
      for( k = 0; k <= (params->o_end-params->o_start); k++ )
	data->f[i][j][k] = 0.0;

  double r_2;

  double d_2 = 4.0*params->radius*params->radius;
  double d_3 = 2.0*params->radius*d_2;
  double d_4 = d_2*d_2;
  double d_5 = d_3*d_2;
  double d_6 = d_4*d_2;
  double d_7 = d_5*d_2;
  double d_8 = d_4*d_4;
  double d_9 = d_7*d_2;
  double d_10 = d_8*d_2;
  double d_11 = d_9*d_2;
  double d_12 = d_10*d_2;
  double d_13 = d_11*d_2;
  double d_14 = d_12*d_2;
  double d_15 = d_13*d_2;
  double d_16 = d_14*d_2;
  double d_17 = d_15*d_2;

  for( p = 0; p < data->n_stored_particles; p++ ){
    cell[0] = floor( data->particles[p].x / params->x * params->m );
    cell[1] = floor( data->particles[p].y / params->y * params->n );
    cell[2] = floor( data->particles[p].z / params->z * params->o );
    i_start = MAX( cell[0] - params->ghosts, params->m_start );
    i_end = MIN( cell[0] + params->ghosts, params->m_end );
    j_start = MAX( cell[1] - params->ghosts, params->n_start );
    j_end = MIN( cell[1] + params->ghosts, params->n_end );
    k_start = MAX( cell[2] - params->ghosts, params->o_start );
    k_end = MIN( cell[2] + params->ghosts, params->o_end );

    for( i = i_start; i <= i_end; i++ )
      for( j = j_start; j <= j_end; j++ )
	for( k = k_start; k <= k_end; k++ ) {
	  switch (params->distribution) {
	    case polynomial_deg_6:
	      r_2 = SQUARE( i*params->hx - data->particles[p].x ) +
		SQUARE( j*params->hy - data->particles[p].y ) +
		SQUARE( k*params->hz - data->particles[p].z );
	      r = sqrt (r_2);
	      if (r < (params->radius)) {
		double dr = d_2 - 4.0*r_2;
		double dr_2 = dr*dr;
		data->f[i-params->m_start][j-params->n_start][k-params->o_start] +=
		  data->particles[p].q * 315.0 * dr * dr_2 / (8.0 * d_9 * PP3MG_PI);
	      }
	      break;
	    case polynomial_deg_10:
	      r_2 = SQUARE( i*params->hx - data->particles[p].x ) +
		SQUARE( j*params->hy - data->particles[p].y ) +
		SQUARE( k*params->hz - data->particles[p].z );
	      r = sqrt (r_2);
	      if (r < (params->radius)) {
		double dr = d_2 - 4.0*r_2;
		double dr_2 = dr*dr;
		double dr_4 = dr_2*dr_2;
		data->f[i-params->m_start][j-params->n_start][k-params->o_start] +=
		  data->particles[p].q * 9009.0 * dr * dr_4 / (128.0 * d_13 * PP3MG_PI);
	      }
	      break;
	    case polynomial_deg_14:
	      r_2 = SQUARE( i*params->hx - data->particles[p].x ) +
		SQUARE( j*params->hy - data->particles[p].y ) +
		SQUARE( k*params->hz - data->particles[p].z );
	      r = sqrt (r_2);
	      if (r < (params->radius)) {
		double dr = d_2 - 4.0*r_2;
		double dr_2 = dr*dr;
		double dr_4 = dr_2*dr_2;
		double dr_8 = dr_4*dr_4;
		double dr_16 = dr_8*dr_8;
		data->f[i-params->m_start][j-params->n_start][k-params->o_start] +=
		  data->particles[p].q * 109395.0 * dr * dr_2 * dr_4 / (1024.0 * d_17 * PP3MG_PI);
	      }
	      break;
	    case spline_deg_4:
	      r_2 = SQUARE( i*params->hx - data->particles[p].x ) +
		SQUARE( j*params->hy - data->particles[p].y ) +
		SQUARE( k*params->hz - data->particles[p].z );
	      double r_4 = r_2 * r_2;
	      r = sqrt (r_2);
	      if( r < (params->radius / 3.0) )
		data->f[i-params->m_start][j-params->n_start][k-params->o_start] +=
		  data->particles[p].q *
		  ( 27.0 * ( 81.0 * r_4 -
			     54.0 * r_2 * params->radius_2 +
			     11.0 * params->radius_4 ) ) /
		  ( 32.0 * PP3MG_PI * params->radius_7 );
	      else if( r < ( 2.0 * params->radius / 3.0 ) )
		data->f[i-params->m_start][j-params->n_start][k-params->o_start] +=
		  data->particles[p].q *
		  ( 27 * ( (-9) * r_2 +
			   6 * r * params->radius + params->radius_2 ) *
		    ( 27 * r_2 - 42 * r * params->radius +
		      17 * params->radius_2 ) ) /
		  ( 64 * PP3MG_PI * params->radius_7 );
	      else if( r < params->radius )
		data->f[i-params->m_start][j-params->n_start][k-params->o_start] +=
		  data->particles[p].q *
		  ( 2187 * SQUARE(SQUARE( r - params->radius )) ) /
		  ( 64 * PP3MG_PI * params->radius_7 );
	      break;
	    }
	}
  }

  /* Ensuring compatibility condition */
  {
    double sum_local = 0.0;
    double sum = 0.0;

    for (i=0;i<=params->m_end-params->m_start;i++) {
      for (j=0;j<=params->n_end-params->n_start;j++) {
	for (k=0;k<=params->o_end-params->o_start;k++) {
	  sum_local += data->f[i][j][k];
	}
      }
    }
    MPI_Allreduce( &sum_local, &sum, 1, MPI_DOUBLE, MPI_SUM,
		   params->mpi_comm_cart );
#ifdef DEBUG
    {
      int mpi_rank;
      
      MPI_Comm_rank( params->mpi_comm_cart, &mpi_rank );
      printf("DEBUG (rank=%d): sum_local = %f, sum = %f\n",
	     mpi_rank,sum_local,sum);
    }
#endif

    for (i=0;i<=params->m_end-params->m_start;i++) {
      for (j=0;j<=params->n_end-params->n_start;j++) {
	for (k=0;k<=params->o_end-params->o_start;k++) {
	  data->f[i][j][k] = data->f[i][j][k] -
	    sum/(params->m*params->n*params->o);
	}
      }
    }
  }

  if (params->discretization == order_4_compact) {
    /* Preparing right hand side for 4th order compact solver */
    for( i = 0; i <= (params->m_end-params->m_start+2*params->ghosts); i++ )
      for( j = 0; j <= (params->n_end-params->n_start+2*params->ghosts); j++ )
	for( k = 0; k <= (params->o_end-params->o_start+2*params->ghosts); k++ )
	  data->f_ghosted[i][j][k] = 0.0;
    
    for( i = 0; i <= (params->m_end-params->m_start); i++ )
      for( j = 0; j <= (params->n_end-params->n_start); j++ )
	for( k = 0; k <= (params->o_end-params->o_start); k++ )
	  data->f_ghosted[i+params->ghosts][j+params->ghosts][k+params->ghosts] = data->f[i][j][k];
    
    pp3mg_update_ghosts( data->f_ghosted, params->m_end-params->m_start+1,
			 params->n_end-params->n_start+1,
  		       params->o_end-params->o_start+1,
			 params->ghosts, params->mpi_comm_cart );
    
    for( i = 0; i <= (params->m_end-params->m_start); i++ )
      for( j = 0; j <= (params->n_end-params->n_start); j++ )
	for( k = 0; k <= (params->o_end-params->o_start); k++ )
	  data->f[i][j][k] = (1.0/12) *
	    (6 * data->f_ghosted[i+params->ghosts][j+params->ghosts][k+params->ghosts] +
	     data->f_ghosted[i-1+params->ghosts][j+params->ghosts][k+params->ghosts] +
	     data->f_ghosted[i+1+params->ghosts][j+params->ghosts][k+params->ghosts] +
	     data->f_ghosted[i+params->ghosts][j-1+params->ghosts][k+params->ghosts] +
	     data->f_ghosted[i+params->ghosts][j+1+params->ghosts][k+params->ghosts] +
	     data->f_ghosted[i+params->ghosts][j+params->ghosts][k-1+params->ghosts] +
	     data->f_ghosted[i+params->ghosts][j+params->ghosts][k+1+params->ghosts] );
  }
    
  /* -------------------------------------------------------------------
   *
   * Initializing variables for solver
   *
   * ------------------------------------------------------------------- */

  /* Solving linear system */
  for( i = 0; i <= (params->m_end-params->m_start); i++ )
    for( j = 0; j <= (params->n_end-params->n_start); j++ )
      for( k = 0; k <= (params->o_end-params->o_start); k++ )
	data->u[i][j][k] = 0.0;
				
/*   size = 7; */
/*   values = (double*) malloc( size*sizeof( double ) ); */
/*   xoff   = (int*) malloc( size*sizeof( int ) ); */
/*   yoff   = (int*) malloc( size*sizeof( int ) ); */
/*   zoff   = (int*) malloc( size*sizeof( int ) ); */
/*   values[0] = 2.0/(params->hx*params->hx) + 2.0/(params->hy*params->hy) + 2.0/(params->hz*params->hz); xoff[0] = 0; yoff[0] = 0; zoff[0] = 0; */
/*   values[1] = -1.0/(params->hx*params->hx); xoff[1] = -1; yoff[1] =  0; zoff[1] =  0; */
/*   values[2] = -1.0/(params->hx*params->hx); xoff[2] =  1; yoff[2] =  0; zoff[2] =  0; */
/*   values[3] = -1.0/(params->hy*params->hy); xoff[3] =  0; yoff[3] = -1; zoff[3] =  0; */
/*   values[4] = -1.0/(params->hy*params->hy); xoff[4] =  0; yoff[4] =  1; zoff[4] =  0; */
/*   values[5] = -1.0/(params->hz*params->hz); xoff[5] =  0; yoff[5] =  0; zoff[5] = -1; */
/*   values[6] = -1.0/(params->hz*params->hz); xoff[6] =  0; yoff[6] =  0; zoff[6] =  1; */

/*   size = 19; */
/*   values = (double*) malloc( size*sizeof( double ) ); */
/*   xoff   = (int*) malloc( size*sizeof( int ) ); */
/*   yoff   = (int*) malloc( size*sizeof( int ) ); */
/*   zoff   = (int*) malloc( size*sizeof( int ) ); */
	
/*   values[0] = 24.0 / ( 6 * pow( params->hx, 2 ) ); */
/*   for( i = 1; i <= 6; i++ ) */
/*     values[i] = (-2.0) / ( 6 * pow( params->hx, 2 ) ); */
/*   for( i = 7; i <= 18; i++ ) */
/*     values[i] = (-1.0) / ( 6 * pow( params->hx, 2 ) ); */
		
/*   xoff[0]  = 0; */
/*   xoff[1]  = 1; */
/*   xoff[2]  = -1; */
/*   for( i = 3; i <= 6; i++ ) */
/*     xoff[i] = 0; */
/*   xoff[7]  = 1; */
/*   xoff[8]  = -1; */
/*   xoff[9]  = 1; */
/*   xoff[10] = -1; */
/*   xoff[11] = 1; */
/*   xoff[12] = 1; */
/*   xoff[13] = -1; */
/*   xoff[14] = -1; */
/*   for( i = 15; i <= 18; i++ ) */
/*     xoff[i] = 0; */
		
/*   for( i = 0; i <= 2; i++ ) */
/*     yoff[i] = 0; */
/*   yoff[3]  = 1; */
/*   yoff[4]  = -1; */
/*   yoff[5]  = 0; */
/*   yoff[6]  = 0; */
/*   yoff[7]  = -1; */
/*   yoff[8]  = -1; */
/*   yoff[9]  = 1; */
/*   yoff[10] = 1; */
/*   for( i = 11; i <= 14; i++ ) */
/*     yoff[i] = 0; */
/*   yoff[15] = 1; */
/*   yoff[16] = 1; */
/*   yoff[17] = -1; */
/*   yoff[18] = -1; */
	
/*   for( i = 0; i <= 4; i++ ) */
/*     zoff[i] = 0; */
/*   zoff[5]  = 1; */
/*   zoff[6]  = -1; */
/*   for( i = 7; i <= 10; i++ ) */
/*     zoff[i] = 0; */
/*   zoff[11] = 1; */
/*   zoff[12] = -1; */
/*   zoff[13] = 1; */
/*   zoff[14] = -1; */
/*   zoff[15] = 1; */
/*   zoff[16] = -1; */
/*   zoff[17] = 1; */
/*   zoff[18] = -1; */
  
  switch (params->discretization) {
  case order_2:
    /* standard 2-nd order solver */
    size = 7;
    xoff   = (int*) malloc( size*sizeof( int ) );
    yoff   = (int*) malloc( size*sizeof( int ) );
    zoff   = (int*) malloc( size*sizeof( int ) );
    values = (double*) malloc( size*sizeof( double ) );
    xoff[ 0] =  0; yoff[ 0] =  0; zoff[ 0] =  0; values[ 0] = 2.0/(params->hx*params->hx)+2.0/(params->hy*params->hy)+2.0/(params->hz*params->hz);
    xoff[ 1] =  1; yoff[ 1] =  0; zoff[ 1] =  0; values[ 1] = -1.0/(params->hx*params->hx);
    xoff[ 2] = -1; yoff[ 2] =  0; zoff[ 2] =  0; values[ 2] = -1.0/(params->hx*params->hx);
    xoff[ 3] =  0; yoff[ 3] =  1; zoff[ 3] =  0; values[ 3] = -1.0/(params->hy*params->hy);
    xoff[ 4] =  0; yoff[ 4] = -1; zoff[ 4] =  0; values[ 4] = -1.0/(params->hy*params->hy);
    xoff[ 5] =  0; yoff[ 5] =  0; zoff[ 5] =  1; values[ 5] = -1.0/(params->hz*params->hz);
    xoff[ 6] =  0; yoff[ 6] =  0; zoff[ 6] = -1; values[ 6] = -1.0/(params->hz*params->hz);
    
  case order_4_compact:
    /* 4-th order compact solver */
    size = 19;
    xoff   = (int*) malloc( size*sizeof( int ) );
    yoff   = (int*) malloc( size*sizeof( int ) );
    zoff   = (int*) malloc( size*sizeof( int ) );
    values = (double*) malloc( size*sizeof( double ) );
    
    xoff[ 0] =  0; yoff[ 0] =  0; zoff[ 0] =  0; values[ 0] = (2.0-2.0/3.0)*(1.0/(params->hx*params->hx)+1.0/(params->hy*params->hy)+1.0/(params->hz*params->hz));
    xoff[ 1] = -1; yoff[ 1] =  0; zoff[ 1] =  0; values[ 1] = -1.0/(params->hx*params->hx) + (2.0/(params->hx*params->hx)+1.0/(params->hy*params->hy)+1.0/(params->hz*params->hz))/6.0;
    xoff[ 2] =  1; yoff[ 2] =  0; zoff[ 2] =  0; values[ 2] = -1.0/(params->hx*params->hx) + (2.0/(params->hx*params->hx)+1.0/(params->hy*params->hy)+1.0/(params->hz*params->hz))/6.0;
    xoff[ 3] =  0; yoff[ 3] = -1; zoff[ 3] =  0; values[ 3] = -1.0/(params->hy*params->hy) + (1.0/(params->hx*params->hx)+2.0/(params->hy*params->hy)+1.0/(params->hz*params->hz))/6.0;
    xoff[ 4] =  0; yoff[ 4] =  1; zoff[ 4] =  0; values[ 4] = -1.0/(params->hy*params->hy) + (1.0/(params->hx*params->hx)+2.0/(params->hy*params->hy)+1.0/(params->hz*params->hz))/6.0;
    xoff[ 5] =  0; yoff[ 5] =  0; zoff[ 5] = -1; values[ 5] = -1.0/(params->hz*params->hz) + (1.0/(params->hx*params->hx)+1.0/(params->hy*params->hy)+2.0/(params->hz*params->hz))/6.0;
    xoff[ 6] =  0; yoff[ 6] =  0; zoff[ 6] =  1; values[ 6] = -1.0/(params->hz*params->hz) + (1.0/(params->hx*params->hx)+1.0/(params->hy*params->hy)+2.0/(params->hz*params->hz))/6.0;
    xoff[ 7] = -1; yoff[ 7] = -1; zoff[ 7] =  0; values[ 7] = -(1.0/(params->hx*params->hx)+1.0/(params->hy*params->hy))/12.0;
    xoff[ 8] = -1; yoff[ 8] =  1; zoff[ 8] =  0; values[ 8] = -(1.0/(params->hx*params->hx)+1.0/(params->hy*params->hy))/12.0;
    xoff[ 9] =  1; yoff[ 9] = -1; zoff[ 9] =  0; values[ 9] = -(1.0/(params->hx*params->hx)+1.0/(params->hy*params->hy))/12.0;
    xoff[10] =  1; yoff[10] =  1; zoff[10] =  0; values[10] = -(1.0/(params->hx*params->hx)+1.0/(params->hy*params->hy))/12.0;
    xoff[11] = -1; yoff[11] =  0; zoff[11] = -1; values[11] = -(1.0/(params->hx*params->hx)+1.0/(params->hz*params->hz))/12.0;
    xoff[12] = -1; yoff[12] =  0; zoff[12] =  1; values[12] = -(1.0/(params->hx*params->hx)+1.0/(params->hz*params->hz))/12.0;
    xoff[13] =  1; yoff[13] =  0; zoff[13] = -1; values[13] = -(1.0/(params->hx*params->hx)+1.0/(params->hz*params->hz))/12.0;
    xoff[14] =  1; yoff[14] =  0; zoff[14] =  1; values[14] = -(1.0/(params->hx*params->hx)+1.0/(params->hz*params->hz))/12.0;
    xoff[15] =  0; yoff[15] = -1; zoff[15] = -1; values[15] = -(1.0/(params->hy*params->hy)+1.0/(params->hz*params->hz))/12.0;
    xoff[16] =  0; yoff[16] = -1; zoff[16] =  1; values[16] = -(1.0/(params->hy*params->hy)+1.0/(params->hz*params->hz))/12.0;
    xoff[17] =  0; yoff[17] =  1; zoff[17] = -1; values[17] = -(1.0/(params->hy*params->hy)+1.0/(params->hz*params->hz))/12.0;
    xoff[18] =  0; yoff[18] =  1; zoff[18] =  1; values[18] = -(1.0/(params->hy*params->hy)+1.0/(params->hz*params->hz))/12.0;
    
    nu1 = 3;
    nu2 = 3;
    omega = 0.8;
    break;
    
  case order_6:
    /* 6-th order solver */
    size = 19;
    xoff   = (int*) malloc( size*sizeof( int ) );
    yoff   = (int*) malloc( size*sizeof( int ) );
    zoff   = (int*) malloc( size*sizeof( int ) );
    values = (double*) malloc( size*sizeof( double ) );
    
    xoff[ 0] =  0; yoff[ 0] =  0; zoff[ 0] =  0; values[ 0] = (490.0/180.0)*(1.0/(params->hx*params->hx)+1.0/(params->hy*params->hy)+1.0/(params->hz*params->hz));
    xoff[ 1] =  1; yoff[ 1] =  0; zoff[ 1] =  0; values[ 1] = (-270.0/180.0)*(1.0/(params->hx*params->hx));
    xoff[ 2] = -1; yoff[ 2] =  0; zoff[ 2] =  0; values[ 2] = (-270.0/180.0)*(1.0/(params->hx*params->hx));
    xoff[ 3] =  0; yoff[ 3] =  1; zoff[ 3] =  0; values[ 3] = (-270.0/180.0)*(1.0/(params->hy*params->hy));
    xoff[ 4] =  0; yoff[ 4] = -1; zoff[ 4] =  0; values[ 4] = (-270.0/180.0)*(1.0/(params->hy*params->hy));
    xoff[ 5] =  0; yoff[ 5] =  0; zoff[ 5] =  1; values[ 5] = (-270.0/180.0)*(1.0/(params->hz*params->hz));
    xoff[ 6] =  0; yoff[ 6] =  0; zoff[ 6] = -1; values[ 6] = (-270.0/180.0)*(1.0/(params->hz*params->hz));
    xoff[ 7] =  2; yoff[ 7] =  0; zoff[ 7] =  0; values[ 7] = (27.0/180.0)*(1.0/(params->hx*params->hx));
    xoff[ 8] = -2; yoff[ 8] =  0; zoff[ 8] =  0; values[ 8] = (27.0/180.0)*(1.0/(params->hx*params->hx));
    xoff[ 9] =  0; yoff[ 9] =  2; zoff[ 9] =  0; values[ 9] = (27.0/180.0)*(1.0/(params->hy*params->hy));
    xoff[10] =  0; yoff[10] = -2; zoff[10] =  0; values[10] = (27.0/180.0)*(1.0/(params->hy*params->hy));
    xoff[11] =  0; yoff[11] =  0; zoff[11] =  2; values[11] = (27.0/180.0)*(1.0/(params->hz*params->hz));
    xoff[12] =  0; yoff[12] =  0; zoff[12] = -2; values[12] = (27.0/180.0)*(1.0/(params->hz*params->hz));
    xoff[13] =  3; yoff[13] =  0; zoff[13] =  0; values[13] = (-2.0/180.0)*(1.0/(params->hx*params->hx));
    xoff[14] = -3; yoff[14] =  0; zoff[14] =  0; values[14] = (-2.0/180.0)*(1.0/(params->hx*params->hx));
    xoff[15] =  0; yoff[15] =  3; zoff[15] =  0; values[15] = (-2.0/180.0)*(1.0/(params->hy*params->hy));
    xoff[16] =  0; yoff[16] = -3; zoff[16] =  0; values[16] = (-2.0/180.0)*(1.0/(params->hy*params->hy));
    xoff[17] =  0; yoff[17] =  0; zoff[17] =  3; values[17] = (-2.0/180.0)*(1.0/(params->hz*params->hz));
    xoff[18] =  0; yoff[18] =  0; zoff[18] = -3; values[18] = (-2.0/180.0)*(1.0/(params->hz*params->hz));
    
    nu1 = 10;
    nu2 = 10;
    omega = 0.8;
    break;
  }

  /* -------------------------------------------------------------------
   *
   * Solving
   *
   * ------------------------------------------------------------------- */
  ret = mg( data->u, data->f, params->maxiter, params->tol, params->m, params->n, params->o, 
      params->m_start, params->m_end, params->n_start, params->n_end,
      params->o_start, params->o_end,
      1, nu1, nu2, omega, size, values,
      xoff, yoff, zoff, params->mpi_comm_cart, 2);

  free(xoff);
  free(yoff);
  free(zoff);
  free(values);

  /* -------------------------------------------------------------------
   *
   * Interpolating energies to particles
   *
   * ------------------------------------------------------------------- */
  for( p = 0; p < data->n_local_particles; p++ ){
    data->particles[p].e  = 0.0;
    data->particles[p].fx = 0.0;
    data->particles[p].fy = 0.0;
    data->particles[p].fz = 0.0;
  }

  for( i = 0; i <= (params->m_end-params->m_start+2*params->ghosts); i++ )
    for( j = 0; j <= (params->n_end-params->n_start+2*params->ghosts); j++ )
      for( k = 0; k <= (params->o_end-params->o_start+2*params->ghosts); k++ )
	data->u_ghosted[i][j][k] = 0.0;

  for( i = 0; i <= params->m_end-params->m_start; i++ )
    for( j = 0; j <= params->n_end-params->n_start; j++ )
      for( k = 0; k <= params->o_end-params->o_start; k++ )
	data->u_ghosted[i+params->ghosts][j+params->ghosts][k+params->ghosts] = data->u[i][j][k];
				
  pp3mg_update_ghosts( data->u_ghosted, params->m_end-params->m_start+1,
		       params->n_end-params->n_start+1, params->o_end-params->o_start+1,
		       params->ghosts, params->mpi_comm_cart );

  interp_poly( params->degree, data->particles, data->hoc, data->ll, 
	       data->u_ghosted, data->n_local_particles,
	       params->m_start, params->m_end, params->n_start, params->n_end, 
	       params->o_start, params->o_end, params->ghosts, 
	       params->hx, params->hy, params->hz );

  /* ----------------------------------------------------------------------------
   *
   * Near-field correction and energy calculation
   *
   * ---------------------------------------------------------------------------- */

  /* Evaluating pairwise potential and calculating energy */
  for( i = params->ghosts; i <= (params->m_end-params->m_start+params->ghosts); i++ )
    for( j = params->ghosts; j <= (params->n_end-params->n_start+params->ghosts); j++ )
      for( k = params->ghosts; k <= (params->o_end-params->o_start+params->ghosts); k++ ){
	p1 = data->hoc[i][j][k];
	while( p1 >= 0 ){
	  for( ii = i-params->ghosts; ii <= i+params->ghosts; ii++ )
	    for( jj = j-params->ghosts; jj <= j+params->ghosts; jj++ )
	      for( kk = k-params->ghosts; kk <= k+params->ghosts; kk++ ){
		p2 = data->hoc[ii][jj][kk];
		while( p2 >= 0 ){
		  if( p1 != p2 ){
		    double r_2, r_4, r_6;
		    rx = data->particles[p1].x - data->particles[p2].x;
		    ry = data->particles[p1].y - data->particles[p2].y;
		    rz = data->particles[p1].z - data->particles[p2].z;
		    r_2 = SQUARE( rx ) + SQUARE( ry ) + SQUARE( rz );
		    r = sqrt( r_2 );
		    r_4 = r_2 * r_2;
		    r_6 = r_4 * r_2;

		    /* Calculate energy */
		    if( r < params->radius ){
		      switch (params->distribution) {
		      case polynomial_deg_6:
			val = ((((8960.0*r_2 - 11520.0*d_2)*r_2 + 6048.0*d_4)*r_2 - 1680.0*d_6)*r_2 + 315.0*d_8)/(64.0*d_9);
			break;

		      case polynomial_deg_10:
			val = ((((((946176.0*r_2 - 1677312.0*d_2)*r_2 + 1281280.0*d_4)*r_2 - 549120.0*d_6)*r_2 + 144144.0*d_8)*r_2 - 24024.0*d_10)*r_2 + 3003.0*d_12)/(512.0*d_13);
			break;

		      case polynomial_deg_14:
			val = (((((((25740.0/d_17*r_2 - 58344.0/d_15)*r_2 + 58905.0/d_13)*r_2 - 69615.0/(2.0*d_11))*r_2 + 425425.0/(32.0*d_9))*r_2 - 109395.0/(32.0*d_7))*r_2 + 153153.0/(256.0*d_5))*r_2 - 36465.0/(512.0*d_3))*r_2 + 109395.0/(16384.0*2*params->radius);
			break;

		      case spline_deg_4:
			if( r < (params->radius / 3.0 ) )
			  val = ( (-3645) * r_6 +
				  5103 * r_4 * params->radius_2 -
				  3465 * r_2 * params->radius_4 +
				  1673 * params->radius_6 ) /
			    ( 560 * params->radius_7 );
			else if( r < ( 2 * params->radius / 3.0 ) )
			  val = ( 32805 * r * r_6 -
				  102060 * r_6 * params->radius +
				  107163 * r * r_4 * params->radius_2 -
				  28350 * r_4 * params->radius_3 -
				  16065 * r * r_2 * params->radius_4 +
				  9933 * r * params->radius_6 +
				  10 * params->radius_7 ) /
			    ( 3360 * r * params->radius_7 );
			else
			  val = -(  3645 * r * r_6 -
				    20412 * r_6 * params->radius +
				    45927 * r_4 * r  * params->radius_2 -
				    51030 * r_4 * params->radius_3 +
				    25515 * r_2 * r * params->radius_4 -
				    5103 * r * params->radius_6 +
				    338 * params->radius_7 ) /
			    ( 1120 * r * params->radius_7 );
			break;
		      }
			
		      val = val - 1.0 / r;  
		      data->particles[p1].e = data->particles[p1].e -
			1.0 / ( 4 * PP3MG_PI ) * data->particles[p1].q *
			data->particles[p2].q * val;
		    }

		    /* Calculate forces */
		    if( r < params->radius ){
		      switch (params->distribution) {
		      case spline_deg_4:
			if( r < params->radius / 3.0 )
			  val = (-9.0) * r * ( 385 * params->radius_4 -
					   1134 * params->radius_2 * r_2 + 
					   1215 * r_4 ) /
			    ( 280 * params->radius_7 );
			else if( r < ( 2 * params->radius / 3.0 ) )
			  val = ((-5.0) * params->radius_7 -
				 16065 * params->radius_4 * r * r_2 -
				 42525 * params->radius_3 * r_4 +
				 214326 * params->radius_2 * r * r_4 -
				 255150 * params->radius * r_6 +
				 98415 * r * r_6  ) /
			    ( 1680 * params->radius_7 * r_2 );
			else 
			  val = -( (-169.0) * params->radius_7 +
				   25515 * params->radius_4 * r * r_2 -
				   76545 * params->radius_3 * r_4 +
				   91854 * params->radius_2 * r * r_4 -
				   51030 * params->radius * r_6 + 
				   10935 * r * r_6 ) /
			    ( 560 * params->radius_7 * r_2 );
			break;
		      case polynomial_deg_6:
			val = ((((71680.0*r_2-69120.0*d_2)*r_2+24192.0*d_4)*r_2-3360.0*d_6)*r)/(64.0*d_9);
			break;
		      case polynomial_deg_10:
			val = ((((((11354112.0*r_2-16773120.0*d_2)*r_2+10250240.0*d_4)*r_2-3294720.0*d_6)*r_2+576576.0*d_8)*r_2-48048.0*d_10)*r)/(512.0*d_13);
			break;
		      case polynomial_deg_14:
			val = ((((((((6747586560.0*r_2-13382713344.0*d_2)*r_2+11581194240.0*d_4)*r_2-5702860800.0*d_6)*r_2+1742540800.0*d_8)*r_2-336061440.0*d_10)*r_2+39207168.0*d_12)*r_2-2333760.0*d_14)*r)/(16384.0*d_17);
			break;
		      }
			
		      val = val/r + 1.0/r / r_2;
		      data->particles[p1].fx += 
			1.0 / ( 4.0 * PP3MG_PI ) * data->particles[p1].q *
			data->particles[p2].q * 
			rx * val;
		      data->particles[p1].fy +=
			1.0 / ( 4.0 * PP3MG_PI ) * data->particles[p1].q *
			data->particles[p2].q * 
			ry * val;
		      data->particles[p1].fz += 
			1.0 / ( 4.0 * PP3MG_PI ) * data->particles[p1].q *
			data->particles[p2].q * 
			rz * val;
		    }
		  } /* if( p1 != p2) */
		  else{
		    /* Self energy correction */
		    switch (params->distribution) {
		    case polynomial_deg_6:
		      data->particles[p1].e = data->particles[p1].e -
			1.0 / ( 4 * PP3MG_PI ) *
			SQUARE( data->particles[p1].q ) *
			315.0 / ( 64.0 * 2.0 * params->radius );
		      break;

		    case polynomial_deg_10:
		      data->particles[p1].e = data->particles[p1].e -
			1.0 / ( 4 * PP3MG_PI ) *
			SQUARE( data->particles[p1].q ) *
			3003.0 / ( 512.0 * 2.0 * params->radius );
		      break;

		    case polynomial_deg_14:
		      data->particles[p1].e = data->particles[p1].e -
			1.0 / ( 4 * PP3MG_PI ) *
			SQUARE( data->particles[p1].q ) *
			109395.0 / ( 16384.0 * 2.0 * params->radius );
		      break;

		    case spline_deg_4:
		      data->particles[p1].e = data->particles[p1].e -
			1.0 / ( 4 * PP3MG_PI ) *
			SQUARE( data->particles[p1].q ) *
			239 / ( 80.0 * params->radius );
		      break;
		    }
		  }			
		  p2 = data->ll[p2];
		} /* while( p2 >= 0) */
	      } /* for(kk) */
	  p1 = data->ll[p1];
	} /* while(p1>=0) */
      } /* for(k) */

  /* ----------------------------------------------------------------------------
   *
   * Copy output data from internal data structures
   *
   * ---------------------------------------------------------------------------- */
  for( p = 0; p < data->n_local_particles; p++ ){
    e[p]  = data->particles[p].e;
    fx[p] = data->particles[p].fx;
    fy[p] = data->particles[p].fy;
    fz[p] = data->particles[p].fz;
  }
			
				
}

void pp3mg_free( pp3mg_data* data, pp3mg_parameters* params )
{
  int i, j;
	
  MPI_Comm_free (&(params->mpi_comm_cart));
  free( data->particles );
  free( data->ll );
  free( data->f[0][0] );
  free( data->f_ghosted[0][0] );
  free( data->u[0][0] );
  free( data->u_ghosted[0][0] );
	
  for( i = 0; i <= (params->m_end-params->m_start); i++ ){
    free( data->f[i] );
    free( data->u[i] );
  }
	
  for( i = 0; i <= (params->m_end-params->m_start+2*params->ghosts); i++ ){
    for ( j = 0; j <= (params->n_end-params->n_start+2*params->ghosts); j++ ){
      free( data->hoc[i][j] );
    }
    free( data->hoc[i] );
    free( data->f_ghosted[i] );
    free( data->u_ghosted[i] );
  }
  free( data->hoc );
  free( data->f );
  free( data->f_ghosted );
  free( data->u );
  free( data->u_ghosted );
}
