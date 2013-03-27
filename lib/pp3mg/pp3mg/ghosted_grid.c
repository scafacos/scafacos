/* 
 * ghosted_grid.c
 *
 * This file contains the function "update_ghosts" to update
 * "boundary data" with neighboring processors.
 *
 * Based on Fortran code by Matthias Bolten.
 *
 * Authors: Matthias Bolten, Stephanie Friedhoff
 * Created: 2009/06/01
 *
 * Copyright 2009, 2010, 2011, 2012 Matthias Bolten, Stephanie Friedhoff
 * All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

#include "ghosted_grid.h"

void pp3mg_update_ghosts( double*** u, int m, int n, int o, int ghosts, MPI_Comm mpi_comm_cart )
{
  /* MPI variables */
  int mpi_dims[3];
  int mpi_periods[3];
  int mpi_coords[3];
  int mpi_self;
  int mpi_left, mpi_right, mpi_lower, mpi_upper, mpi_back, mpi_front;
  MPI_Request mpi_req[2];
  MPI_Status mpi_stat[2];

  /* Other variables */
  double* buf_left, *buf_right;
  double* buf_lower, *buf_upper;
  double* buf_back, *buf_front;
  int i, j, k;
  int count;

  /* Initializing MPI variables */
  MPI_Comm_rank( mpi_comm_cart, &mpi_self );
  MPI_Cart_get( mpi_comm_cart, 3, mpi_dims, mpi_periods, mpi_coords );
  MPI_Cart_shift( mpi_comm_cart, 0, 1, &mpi_left, &mpi_right );
  MPI_Cart_shift( mpi_comm_cart, 1, 1, &mpi_lower, &mpi_upper );
  MPI_Cart_shift( mpi_comm_cart, 2, 1, &mpi_back, &mpi_front );
  mpi_req[0] = MPI_REQUEST_NULL; 
  mpi_req[1] = MPI_REQUEST_NULL;
	
  /* ----------------------------------------------------------------------------
   *
   * x-direction
   *
   * ---------------------------------------------------------------------------- */

  /* Memory for transfer to left and right neighbors */
  buf_left  = (double*) malloc( ghosts * (n + 2*ghosts) * (o + 2*ghosts) * sizeof( double ) );
  buf_right = (double*) malloc( ghosts * (n + 2*ghosts) * (o + 2*ghosts) * sizeof( double ) );

  /* Sending ghosts to left neighbor */
  if( mpi_left == mpi_self ){
    if( mpi_periods[0] ){	
      for( i = ghosts; i < 2*ghosts; i++ )
	for( j = 0; j < n+2*ghosts; j++ )
	  for( k = 0; k < o+2*ghosts; k++ )
	    u[i+m][j][k] = u[i][j][k];
    }
  }
  else{
    count = 0;
    for( i = ghosts; i < 2*ghosts; i++ )
      for( j = 0; j < n+2*ghosts; j++ )
	for( k = 0; k < o+2*ghosts; k++ ){
	  buf_left[count] = u[i][j][k];
	  count++;					
	}
    MPI_Isend( (void*) buf_left, count, MPI_DOUBLE, mpi_left, 1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( (void*) buf_right, count, MPI_DOUBLE, mpi_right, 1, mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall(2, mpi_req, mpi_stat );
	       
    count = 0;
    for( i = ghosts; i < 2*ghosts; i++ )
      for( j = 0; j < n+2*ghosts; j++ )
	for( k = 0; k < o+2*ghosts; k++ ){
	  u[i+m][j][k] = buf_right[count];
	  count++;					
	}		
  }				
				
  /* Sending ghosts to right neighbor */
  if( mpi_right == mpi_self ){
    if( mpi_periods[0] ){
      for( i = m; i < m+ghosts; i++ )
	for( j = 0; j < n+2*ghosts; j++ )
	  for( k = 0; k < o+2*ghosts; k++ )
	    u[i-m][j][k] = u[i][j][k];
    }
  }
  else{
    count = 0;
    for( i = m; i < m+ghosts; i++ )
      for( j = 0; j < n+2*ghosts; j++ )
	for( k = 0; k < o+2*ghosts; k++ ){
	  buf_right[count] = u[i][j][k];
	  count++;					
	}
    MPI_Isend( (void*) buf_right, count, MPI_DOUBLE, mpi_right, 1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( (void*) buf_left, count, MPI_DOUBLE, mpi_left, 1, mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall(2, mpi_req, mpi_stat );
	       
    count = 0;
    for( i = m; i < m+ghosts; i++ )
      for( j = 0; j < n+2*ghosts; j++ )
	for( k = 0; k < o+2*ghosts; k++ ){
	  u[i-m][j][k] = buf_left[count];
	  count++;					
	}
  }

  /* Freeing memory */
  free( buf_left );
  free( buf_right );
	
  /* ----------------------------------------------------------------------------
   *
   * y-direction
   *
   * ---------------------------------------------------------------------------- */

  /* Memory for transfer to lower and upper neighbors */
  buf_lower = (double*) malloc( (m + 2*ghosts) * ghosts * (o + 2*ghosts) * sizeof( double ) );
  buf_upper = (double*) malloc( (m + 2*ghosts) * ghosts * (o + 2*ghosts) * sizeof( double ) );

  /* Sending ghosts to lower neighbor */
  if( mpi_lower == mpi_self ){
    if( mpi_periods[1] ){
      for( i = 0; i < m+2*ghosts; i++ )
	for( j = ghosts; j < 2*ghosts; j++ )
	  for( k = 0; k < o+2*ghosts; k++ )
	    u[i][j+n][k] = u[i][j][k];
    }
  }
  else{
    count = 0;
    for( i = 0; i < m+2*ghosts; i++ )
      for( j = ghosts; j < 2*ghosts; j++ )
	for( k = 0; k < o+2*ghosts; k++ ){
	  buf_lower[count] = u[i][j][k];
	  count++;					
	}
    MPI_Isend( (void*) buf_lower, count, MPI_DOUBLE, mpi_lower, 1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( (void*) buf_upper, count, MPI_DOUBLE, mpi_upper, 1, mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );
    count = 0;
    for( i = 0; i < m+2*ghosts; i++ )
      for( j = ghosts; j < 2*ghosts; j++ )
	for( k = 0; k < o+2*ghosts; k++ ){
	  u[i][j+n][k] = buf_upper[count];
	  count++;					
	}			
  }


  /* Sending ghosts to upper neighbor */
  if( mpi_upper == mpi_self ){
    if( mpi_periods[1] ){
      for( i = 0; i < m+2*ghosts; i++ )
	for( j = n; j < n+ghosts; j++ )
	  for( k = 0; k < o+2*ghosts; k++ )
	    u[i][j-n][k] = u[i][j][k];			
    }
  }
  else{
    count = 0;
    for( i = 0; i < m+2*ghosts; i++ )
      for( j = n; j < n+ghosts; j++ )
	for( k = 0; k < o+2*ghosts; k++ ){
	  buf_upper[count] = u[i][j][k];
	  count++;					
	}		
    MPI_Isend( (void*) buf_upper, count, MPI_DOUBLE, mpi_upper, 1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( (void*) buf_lower, count, MPI_DOUBLE, mpi_lower, 1, mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );
    count = 0;
    for( i = 0; i < m+2*ghosts; i++ )
      for( j = n; j < n+ghosts; j++ )
	for( k = 0; k < o+2*ghosts; k++ ){
	  u[i][j-n][k] = buf_lower[count];
	  count++;					
	}	
  }

  /* Freeing memory */
  free( buf_lower );
  free( buf_upper );
	
  /* ----------------------------------------------------------------------------
   *
   * z-direction
   *
   * ---------------------------------------------------------------------------- */

  /* Memory for transfer to back and front neighbors */
  buf_front = (double*) malloc( (m + 2*ghosts) * (n + 2*ghosts) * ghosts * sizeof( double ) );
  buf_back  = (double*) malloc( (m + 2*ghosts) * (n + 2*ghosts) * ghosts * sizeof( double ) );
	
  /* Sending ghosts to back neighbor */
  if( mpi_back == mpi_self ){
    if( mpi_periods[2] ){
      for( i = 0; i < m+2*ghosts; i++ )
	for( j = 0; j < n+2*ghosts; j++ )
	  for( k = ghosts; k < 2*ghosts; k++ )
	    u[i][j][k+o] = u[i][j][k];
    }
  }
  else{
    count = 0;
    for( i = 0; i < m+2*ghosts; i++ )
      for( j = 0; j < n+2*ghosts; j++ )
	for( k = ghosts; k < 2*ghosts; k++ ){
	  buf_back[count] = u[i][j][k];
	  count++;					
	}	
    MPI_Isend( (void*) buf_back, count, MPI_DOUBLE, mpi_back, 1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( (void*) buf_front, count, MPI_DOUBLE, mpi_front, 1, mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );
		
    count = 0;
    for( i = 0; i < m+2*ghosts; i++ )
      for( j = 0; j < n+2*ghosts; j++ )
	for( k = ghosts; k < 2*ghosts; k++ ){
	  u[i][j][k+o] = buf_front[count];
	  count++;
	}	
  }

  /* Sending ghosts to front neighbor */
  if( mpi_front == mpi_self ){
    if( mpi_periods[2] ){
      for( i = 0; i < m+2*ghosts; i++ )
	for( j = 0; j < n+2*ghosts; j++ )
	  for( k = o; k < o+ghosts; k++ )
	    u[i][j][k-o] = u[i][j][k];
    }
  }
  else{
    count = 0;
    for( i = 0; i < m+2*ghosts; i++ )
      for( j = 0; j < n+2*ghosts; j++ )
	for( k = o; k < o+ghosts; k++ ){
	  buf_front[count] = u[i][j][k];
	  count++;					
	}
    MPI_Isend( (void*) buf_front, count, MPI_DOUBLE, mpi_front, 1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( (void*) buf_back, count, MPI_DOUBLE, mpi_back, 1, mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );		
    count = 0;
    for( i = 0; i < m+2*ghosts; i++ )
      for( j = 0; j < n+2*ghosts; j++ )
	for( k = o; k < o+ghosts; k++ ){
	  u[i][j][k-o] = buf_back[count];
	  count++;
	}
  }

  /* Freeing memory */
  free( buf_back );
  free( buf_front );
}
