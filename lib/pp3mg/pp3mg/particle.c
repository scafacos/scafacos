/* 
 * particle.c
 *
 * This file contains the functions 
 *   "update_linked_list"     to update the data structure that stores the 
 *                            particles that are contained in a certain grid
 *                            cell.
 *   "update_particle_ghosts" to update boundary particles with neighboring
 *                            processors.
 *
 * Based on Fortran code by Matthias Bolten.
 *
 * Authors: Matthias Bolten, Stephanie Friedhoff
 * Created: 2009/05/28
 *
 * Copyright 2009, 2010, 2011, 2012 Matthias Bolten, Stephanie Friedhoff
 * All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<mpi.h>

#include "pp3mg.h"

#include "particle.h"

void update_linked_list( int* ll, int*** hoc, pp3mg_particle* particles, 
			 int start, int end, double x, double y, double z, int m, int n, int o,
			 int m_start, int m_end, int n_start, int n_end, int o_start, int o_end )
{

	/* Local variables */
	int p;
	int cell[3];

/* 	printf("m_start = %d, m_end = %d\nn_start = %d, n_end = %d\no_start = %d, o_end = %d\n",m_start,m_end,n_start,n_end,o_start,o_end); */
	for( p = start; p < end; p++ ){
		cell[0] = floor( particles[p].x / x * m );    
		cell[1] = floor( particles[p].y / y * n );	
		cell[2] = floor( particles[p].z / z * o );
/* 		printf("p = %d: cell = [%d][%d][%d] -> [%d][%d][%d]\n",p,cell[0],cell[1],cell[2], */
/* 		       cell[0]-m_start,cell[1]-n_start,cell[2]-o_start); */

		ll[p] = hoc[cell[0]-m_start][cell[1]-n_start][cell[2]-o_start];
		hoc[cell[0]-m_start][cell[1]-n_start][cell[2]-o_start] = p;
		
		if( (cell[0] < m_start) || (m_end < cell[0]) ||
		     (cell[1] < n_start) || (n_end < cell[1]) ||
		    (cell[2] < o_start) || (o_end < cell[2]) ) {
		  printf("x = %f, y = %f, z = %f\n",particles[p].x, particles[p].y, particles[p].z);
		  printf( "%d, %d, %d, %d, %d, %d, %d, %d, %d\n",
			  cell[0], cell[1], cell[2], m_start, m_end,
			  n_start, n_end, o_start, o_end );
		}
	}
}

/* ---------------------------------------------------------------------------------------------- */

void update_particle_ghosts( int** ll_addr, int*** hoc, pp3mg_particle** particles_addr,
			     double x, double y, double z, int m, int n, int o, int ghosts, 
			     int* n_stored_particles, int* max_particles,
			     int m_start, int m_end, int n_start, int n_end, 
			     int o_start, int o_end, MPI_Comm mpi_comm_cart )
{

  /* Local variables */
  pp3mg_particle* particles;
  int* ll;
  
  /* Variables for MPI */
  int mpi_self;
  int mpi_left, mpi_right, mpi_lower, mpi_upper, mpi_back, mpi_front;
  int mpi_count;
  int mpi_blockcounts[2];
  MPI_Aint mpi_offsets[2];
  MPI_Request mpi_req[2];
  MPI_Status mpi_stat[2];
  MPI_Datatype mpi_type_particle, mpi_oldtypes[2];
  
  /* Other variables */
  int count, p, start;
  int i, j, k;
  
  /* Initializing MPI variables */
  MPI_Comm_rank( mpi_comm_cart, &mpi_self );
  MPI_Cart_shift( mpi_comm_cart, 0, 1, &mpi_left, &mpi_right );
  MPI_Cart_shift( mpi_comm_cart, 1, 1, &mpi_lower, &mpi_upper );
  MPI_Cart_shift( mpi_comm_cart, 2, 1, &mpi_back, &mpi_front );
  mpi_req[0] = MPI_REQUEST_NULL; 
  mpi_req[1] = MPI_REQUEST_NULL;
  
  /* Creating particle type for MPI */
  mpi_blockcounts[0] = 8;
  mpi_offsets[0] = 0;
  mpi_oldtypes[0] = MPI_DOUBLE;
  MPI_Type_create_struct( 1, mpi_blockcounts, mpi_offsets, mpi_oldtypes, &mpi_type_particle );
  MPI_Type_commit( &mpi_type_particle );
  
  particles = *particles_addr;
  ll = *ll_addr;

  /* --------------------------------------------------------------------------
   *
   * Sending particles to right neighbor
   *
   * -------------------------------------------------------------------------- */
  count = 0;
  for( i = m_end-m_start+1; i <= m_end-m_start+ghosts; i++ )
    for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
      for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	p = hoc[i][j][k];
	while( p >= 0 ){
	  count++;
	  p = ll[p];
	}
      }
  
#ifdef DEBUG
  printf("Rank %d: n_stored_particles = %d, max_particles = %d, count = %d\n",mpi_self,*n_stored_particles,*max_particles,count);
#endif
  start = *max_particles - count;
  while( start <= ( *n_stored_particles + 27*count ) )
    {
      *max_particles = *max_particles*2;
      *particles_addr = (pp3mg_particle*) realloc(particles,*max_particles*sizeof(pp3mg_particle));
      *ll_addr = (int*) realloc(ll,*max_particles*sizeof(int));
      
      if (*particles_addr == NULL || *ll_addr == NULL)
	{
	  printf("Realloc failed!");
	  exit(1);
	}
      else 
	{
	  particles = *particles_addr;

	  ll = *ll_addr;
	  for (int i=*max_particles>>1; i<*max_particles; i++)
	    ll[i] = -1;

	  start = *max_particles - count;
#ifdef DEBUG
	  printf("Rank %d: Reallocated. Now max_particles = %d\n",mpi_self,*max_particles);
#endif
	}
    }

  count = 0;
  /* If at one end, shift particles */
  if( m_end == (m-1) ){
    for( i = m_end-m_start+1; i <= m_end-m_start+ghosts; i++ )
      for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
	for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count] = particles[p];
	    particles[start+count].x  = particles[p].x - x;
	    count++;
	    p = ll[p];
	  }
	}
  }
  else{
    for( i = m_end-m_start+1; i <= m_end-m_start+ghosts; i++ )
      for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
	for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count]  = particles[p];
	    count++;
	    p = ll[p];
	  }
	}
  }
  
  if( mpi_right == mpi_self ){
    /*
    if( periodic ){
    */
    for( p = 0; p < count; p++ ){
      particles[*n_stored_particles+p]  = particles[start+p];
    }
    *n_stored_particles += count;
    /*
    }
    */
  }
  else{
    MPI_Isend( &particles[start], count, mpi_type_particle, mpi_right, 
	       1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( &particles[*n_stored_particles],
	       *max_particles-*n_stored_particles, mpi_type_particle, mpi_left, 1,
	       mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );
    MPI_Get_count( &mpi_stat[1], mpi_type_particle, &mpi_count );
    count = mpi_count;
    *n_stored_particles += count;
  }

  if( *n_stored_particles >= start )
    {
      printf("Buffer too small!\n");
      exit(1);
    }
  
  /* Updating linked list */
  if( count > 0 )
    update_linked_list( ll, hoc, particles,*n_stored_particles-count, *n_stored_particles, 
			x, y, z, m, n, o, 
			m_start-ghosts, m_end+ghosts,
			n_start-ghosts, n_end+ghosts,
			o_start-ghosts, o_end+ghosts );
  
  /* --------------------------------------------------------------------------
   *
   * Sending particles to left neighbor
   *
   * -------------------------------------------------------------------------- */
  
  count = 0;
  for( i = ghosts; i <= 2*ghosts-1; i++ )
    for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
      for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	p = hoc[i][j][k];
	while( p >= 0 ){
	  count++;
	  p = ll[p];
	}
      }
  
#ifdef DEBUG
  printf("Rank %d: n_stored_particles = %d, max_particles = %d, count = %d\n",mpi_self,*n_stored_particles,*max_particles,count);
#endif
  start = *max_particles - count;
  while( start <= ( *n_stored_particles + 27*count ) )
    {
      *max_particles = *max_particles*2;
      *particles_addr = (pp3mg_particle*) realloc(particles,*max_particles*sizeof(pp3mg_particle));
      *ll_addr = (int*) realloc(ll,*max_particles*sizeof(int));
      
      if (*particles_addr == NULL || *ll_addr == NULL)
	{
	  printf("Realloc failed!");
	  exit(1);
	}
      else 
	{
	  particles = *particles_addr;

	  ll = *ll_addr;
	  for (int i=*max_particles>>1; i<*max_particles; i++)
	    ll[i] = -1;

	  start = *max_particles - count;
#ifdef DEBUG
	  printf("Rank %d: Reallocated. Now max_particles = %d\n",mpi_self,*max_particles);
#endif
	}
    }

  count = 0;
  /* If at one end, shift particles */
  if( m_start == 0 ){
    for( i = ghosts; i <= 2*ghosts-1; i++ )
      for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
	for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x + x;
	    particles[start+count].y  = particles[p].y;
	    particles[start+count].z  = particles[p].z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }else{
    for( i = ghosts; i <= 2*ghosts-1; i++ )
      for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
	for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x;
	    particles[start+count].y  = particles[p].y;
	    particles[start+count].z  = particles[p].z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }
  
  if( mpi_left == mpi_self ){
    /*
    if( periodic ){
    */
    for( p = 0; p < count; p++ ){
      particles[*n_stored_particles+p].x  = particles[start+p].x;
      particles[*n_stored_particles+p].y  = particles[start+p].y;
      particles[*n_stored_particles+p].z  = particles[start+p].z;
      particles[*n_stored_particles+p].q  = particles[start+p].q;
      particles[*n_stored_particles+p].e  = particles[start+p].e;
      particles[*n_stored_particles+p].fx = particles[start+p].fx;
      particles[*n_stored_particles+p].fy = particles[start+p].fy;
      particles[*n_stored_particles+p].fz = particles[start+p].fz;
    }
    
    *n_stored_particles += count;
    /*
    }
    */
  }
  else{
    MPI_Isend( &particles[start], count, mpi_type_particle, mpi_left, 
	       1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( &particles[*n_stored_particles],
	       *max_particles-*n_stored_particles, mpi_type_particle, mpi_right, 1,
	       mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );
    MPI_Get_count( &mpi_stat[1], mpi_type_particle, &mpi_count );
    count = mpi_count;
    *n_stored_particles += count;
  }

  if( *n_stored_particles >= start )
    {
      printf("Buffer too small!\n");
      exit(1);
    }
  
  /* Updating linked list */
  if( count > 0 )
    update_linked_list( ll, hoc, particles,*n_stored_particles-count, *n_stored_particles, 
			x, y, z, m, n, o, 
			m_start-ghosts, m_end+ghosts,
			n_start-ghosts, n_end+ghosts,
			o_start-ghosts, o_end+ghosts );
  
  
  /* --------------------------------------------------------------------------
   *
   * Sending particles to upper neighbor
   *
   * -------------------------------------------------------------------------- */
  
  count = 0;
  for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
    for( j = n_end-n_start+1; j <= n_end-n_start+ghosts; j++ )
      for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	p = hoc[i][j][k];
	while( p >= 0 ){
	  count++;
	  p = ll[p];
	}
      }
  
#ifdef DEBUG
  printf("Rank %d: n_stored_particles = %d, max_particles = %d, count = %d\n",mpi_self,*n_stored_particles,*max_particles,count);
#endif
  start = *max_particles - count;
  while( start <= ( *n_stored_particles + 27*count ) )
    {
      *max_particles = *max_particles*2;
      *particles_addr = (pp3mg_particle*) realloc(particles,*max_particles*sizeof(pp3mg_particle));
      *ll_addr = (int*) realloc(ll,*max_particles*sizeof(int));
      
      if (*particles_addr == NULL || *ll_addr == NULL)
	{
	  printf("Realloc failed!");
	  exit(1);
	}
      else 
	{
	  particles = *particles_addr;

	  ll = *ll_addr;
	  for (int i=*max_particles>>1; i<*max_particles; i++)
	    ll[i] = -1;

	  start = *max_particles - count;
#ifdef DEBUG
	  printf("Rank %d: Reallocated. Now max_particles = %d\n",mpi_self,*max_particles);
#endif
	}
    }

  count = 0;
  /* If at one end, shift particles */
  if( n_end == (n-1) ){
    for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
      for( j = n_end-n_start+1; j <= n_end-n_start+ghosts; j++ )
	for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x;
	    particles[start+count].y  = particles[p].y - y;
	    particles[start+count].z  = particles[p].z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }else{
    for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
      for( j = n_end-n_start+1; j <= n_end-n_start+ghosts; j++ )
	for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x;
	    particles[start+count].y  = particles[p].y;
	    particles[start+count].z  = particles[p].z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }
  
  if( mpi_upper == mpi_self ){
    /*
    if( periodic ){
    */
    for( p = 0; p < count; p++ ){
      particles[*n_stored_particles+p].x  = particles[start+p].x;
      particles[*n_stored_particles+p].y  = particles[start+p].y;
      particles[*n_stored_particles+p].z  = particles[start+p].z;
      particles[*n_stored_particles+p].q  = particles[start+p].q;
      particles[*n_stored_particles+p].e  = particles[start+p].e;
      particles[*n_stored_particles+p].fx = particles[start+p].fx;
      particles[*n_stored_particles+p].fy = particles[start+p].fy;
      particles[*n_stored_particles+p].fz = particles[start+p].fz;
    }
    *n_stored_particles += count;
    /*
    }
    */
  }
  else{
    MPI_Isend( &particles[start], count, mpi_type_particle, mpi_upper, 
	       1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( &particles[*n_stored_particles],
	       *max_particles-*n_stored_particles, mpi_type_particle, mpi_lower, 1,
	       mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );
    MPI_Get_count( &mpi_stat[1], mpi_type_particle, &mpi_count );
    count = mpi_count;
    *n_stored_particles += count;
  }
  
  if( *n_stored_particles >= start )
    {
      printf("Buffer too small!\n");
      exit(1);
    }
  
  /* Updating linked list */
  if( count > 0 )
    update_linked_list( ll, hoc, particles,*n_stored_particles-count, *n_stored_particles, 
			x, y, z, m, n, o, 
			m_start-ghosts, m_end+ghosts,
			n_start-ghosts, n_end+ghosts,
			o_start-ghosts, o_end+ghosts );
  
  /* --------------------------------------------------------------------------
   *
   * Sending particles to lower neighbor
   *
   * -------------------------------------------------------------------------- */
  
  count = 0;
  for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
    for( j = ghosts; j <= 2*ghosts-1; j++ )
      for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	p = hoc[i][j][k];
	while( p >= 0 ){
	  count++;
	  p = ll[p];
	}
      }
  
#ifdef DEBUG
  printf("Rank %d: n_stored_particles = %d, max_particles = %d, count = %d\n",mpi_self,*n_stored_particles,*max_particles,count);
#endif
  start = *max_particles - count;
  while( start <= ( *n_stored_particles + 27*count ) )
    {
      *max_particles = *max_particles*2;
      *particles_addr = (pp3mg_particle*) realloc(particles,*max_particles*sizeof(pp3mg_particle));
      *ll_addr = (int*) realloc(ll,*max_particles*sizeof(int));
      
      if (*particles_addr == NULL || *ll_addr == NULL)
	{
	  printf("Realloc failed!");
	  exit(1);
	}
      else 
	{
	  particles = *particles_addr;

	  ll = *ll_addr;
	  for (int i=*max_particles>>1; i<*max_particles; i++)
	    ll[i] = -1;

	  start = *max_particles - count;
#ifdef DEBUG
	  printf("Rank %d: Reallocated. Now max_particles = %d\n",mpi_self,*max_particles);
#endif
	}
    }

  count = 0;
  /* If at one end, shift particles */
  if( n_start == 0 ){
    for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
      for( j = ghosts; j <= 2*ghosts-1; j++ )
	for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x;
	    particles[start+count].y  = particles[p].y + y;
	    particles[start+count].z  = particles[p].z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }else{
    for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
      for( j = ghosts; j <= 2*ghosts-1; j++ )
	for( k = 0; k <= o_end-o_start+2*ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x;
	    particles[start+count].y  = particles[p].y;
	    particles[start+count].z  = particles[p].z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }
  
  if( mpi_lower == mpi_self ){
    /*
    if( periodic ){
    */
    for( p = 0; p < count; p++ ){
      particles[*n_stored_particles+p].x  = particles[start+p].x;
      particles[*n_stored_particles+p].y  = particles[start+p].y;
      particles[*n_stored_particles+p].z  = particles[start+p].z;
      particles[*n_stored_particles+p].q  = particles[start+p].q;
      particles[*n_stored_particles+p].e  = particles[start+p].e;
      particles[*n_stored_particles+p].fx = particles[start+p].fx;
      particles[*n_stored_particles+p].fy = particles[start+p].fy;
      particles[*n_stored_particles+p].fz = particles[start+p].fz;
    }
    *n_stored_particles += count;
    /*
    }
    */
  }
  else{
    MPI_Isend( &particles[start], count, mpi_type_particle, mpi_lower, 
	       1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( &particles[*n_stored_particles],
	       *max_particles-*n_stored_particles, mpi_type_particle, mpi_upper, 1,
	       mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );
    MPI_Get_count( &mpi_stat[1], mpi_type_particle, &mpi_count );
    count = mpi_count;
    *n_stored_particles += count;
  }
  
  if( *n_stored_particles >= start )
    {
      printf("Buffer too small!\n");
      exit(1);
    }
  
  /* Updating linked list */
  if( count > 0 )
    update_linked_list( ll, hoc, particles,*n_stored_particles-count, *n_stored_particles, 
			x, y, z, m, n, o, 
			m_start-ghosts, m_end+ghosts,
			n_start-ghosts, n_end+ghosts,
			o_start-ghosts, o_end+ghosts );
  
  /* --------------------------------------------------------------------------
   *
   * Sending particles to front neighbor
   *
   * -------------------------------------------------------------------------- */
  
  count = 0;
  for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
    for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
      for( k = o_end-o_start+1; k <= o_end-o_start+ghosts; k++ ){
	p = hoc[i][j][k];
	while( p >= 0 ){
	  count++;
	  p = ll[p];
	}
      }
  
#ifdef DEBUG
  printf("Rank %d: n_stored_particles = %d, max_particles = %d, count = %d\n",mpi_self,*n_stored_particles,*max_particles,count);
#endif
  start = *max_particles - count;
  while( start <= ( *n_stored_particles + 27*count ) )
    {
      *max_particles = *max_particles*2;
      *particles_addr = (pp3mg_particle*) realloc(particles,*max_particles*sizeof(pp3mg_particle));
      *ll_addr = (int*) realloc(ll,*max_particles*sizeof(int));
      
      if (*particles_addr == NULL || *ll_addr == NULL)
	{
	  printf("Realloc failed!");
	  exit(1);
	}
      else 
	{
	  particles = *particles_addr;

	  ll = *ll_addr;
	  for (int i=*max_particles>>1; i<*max_particles; i++)
	    ll[i] = -1;

	  start = *max_particles - count;
#ifdef DEBUG
	  printf("Rank %d: Reallocated. Now max_particles = %d\n",mpi_self,*max_particles);
#endif
	}
    }

  count = 0;
  /* If at one end, shift particles */
  if( o_end == (o-1) ){
    for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
      for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
	for( k = o_end-o_start+1; k <= o_end-o_start+ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x;
	    particles[start+count].y  = particles[p].y;
	    particles[start+count].z  = particles[p].z - z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }else{
    for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
      for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
	for( k = o_end-o_start+1; k <= o_end-o_start+ghosts; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x;
	    particles[start+count].y  = particles[p].y;
	    particles[start+count].z  = particles[p].z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }
  
  if( mpi_front == mpi_self ){
    /*
    if( periodic ){
    */
    for( p = 0; p < count; p++ ){
      particles[*n_stored_particles+p].x  = particles[start+p].x;
      particles[*n_stored_particles+p].y  = particles[start+p].y;
      particles[*n_stored_particles+p].z  = particles[start+p].z;
      particles[*n_stored_particles+p].q  = particles[start+p].q;
      particles[*n_stored_particles+p].e  = particles[start+p].e;
      particles[*n_stored_particles+p].fx = particles[start+p].fx;
      particles[*n_stored_particles+p].fy = particles[start+p].fy;
      particles[*n_stored_particles+p].fz = particles[start+p].fz;
    }
    *n_stored_particles += count;
    /*
    }
    */
  }
  else{
    MPI_Isend( &particles[start], count, mpi_type_particle, mpi_front, 
	       1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( &particles[*n_stored_particles],
	       *max_particles-*n_stored_particles, mpi_type_particle, mpi_back, 1,
	       mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );
    MPI_Get_count( &mpi_stat[1], mpi_type_particle, &mpi_count );
    count = mpi_count;
    *n_stored_particles += count;
  }
  
  if( *n_stored_particles >= start )
    {
      printf("Buffer too small!\n");
      exit(1);
    }
  
  /* Updating linked list */
  if( count > 0 )
    update_linked_list( ll, hoc, particles,*n_stored_particles-count, *n_stored_particles, 
			x, y, z, m, n, o, 
			m_start-ghosts, m_end+ghosts,
			n_start-ghosts, n_end+ghosts,
			o_start-ghosts, o_end+ghosts );
  
  /* --------------------------------------------------------------------------
   *
   * Sending particles to back neighbor
   *
   * -------------------------------------------------------------------------- */
  
  count = 0;
  for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
    for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
      for( k = ghosts; k <= 2*ghosts-1; k++ ){
	      p = hoc[i][j][k];
	      while( p >= 0 ){
		count++;
		p = ll[p];
	      }
      }
  
#ifdef DEBUG
  printf("Rank %d: n_stored_particles = %d, max_particles = %d, count = %d\n",mpi_self,*n_stored_particles,*max_particles,count);
#endif
  start = *max_particles - count;
  while( start <= ( *n_stored_particles + 27*count ) )
    {
      *max_particles = *max_particles*2;
      *particles_addr = (pp3mg_particle*) realloc(particles,*max_particles*sizeof(pp3mg_particle));
      *ll_addr = (int*) realloc(ll,*max_particles*sizeof(int));
      
      if (*particles_addr == NULL || *ll_addr == NULL)
	{
	  printf("Realloc failed!");
	  exit(1);
	}
      else 
	{
	  particles = *particles_addr;

	  ll = *ll_addr;
	  for (int i=*max_particles>>1; i<*max_particles; i++)
	    ll[i] = -1;

	  start = *max_particles - count;
#ifdef DEBUG
	  printf("Rank %d: Reallocated. Now max_particles = %d\n",mpi_self,*max_particles);
#endif
	}
    }

  count = 0;
  /* If at one end, shift particles */
  if( o_start == 0 ){
    for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
      for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
	for( k = ghosts; k <= 2*ghosts-1; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x;
	    particles[start+count].y  = particles[p].y;
	    particles[start+count].z  = particles[p].z + z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }
  else{
    for( i = 0; i <= m_end-m_start+2*ghosts; i++ )
      for( j = 0; j <= n_end-n_start+2*ghosts; j++ )
	for( k = ghosts; k <= 2*ghosts-1; k++ ){
	  p = hoc[i][j][k];
	  while( p >= 0 ){
	    particles[start+count].x  = particles[p].x;
	    particles[start+count].y  = particles[p].y;
	    particles[start+count].z  = particles[p].z;
	    particles[start+count].q  = particles[p].q;
	    particles[start+count].e  = particles[p].e;
	    particles[start+count].fx = particles[p].fx;
	    particles[start+count].fy = particles[p].fy;
	    particles[start+count].fz = particles[p].fz;
	    count++;
	    p = ll[p];
	  }
	}
  }
  
  if( mpi_back == mpi_self ){
    /*
    if( periodic ){
    */
    for( p = 0; p < count; p++ ){
      particles[*n_stored_particles+p].x  = particles[start+p].x;
      particles[*n_stored_particles+p].y  = particles[start+p].y;
      particles[*n_stored_particles+p].z  = particles[start+p].z;
      particles[*n_stored_particles+p].q  = particles[start+p].q;
      particles[*n_stored_particles+p].e  = particles[start+p].e;
      particles[*n_stored_particles+p].fx = particles[start+p].fx;
      particles[*n_stored_particles+p].fy = particles[start+p].fy;
      particles[*n_stored_particles+p].fz = particles[start+p].fz;
    }
    *n_stored_particles += count;
    /*
    }
    */
  }else{
    MPI_Isend( &particles[start], count, mpi_type_particle, mpi_back, 
	       1, mpi_comm_cart, &mpi_req[0] );
    MPI_Irecv( &particles[*n_stored_particles],
	       *max_particles-*n_stored_particles, mpi_type_particle, mpi_front, 1,
	       mpi_comm_cart, &mpi_req[1] );
    MPI_Waitall( 2, mpi_req, mpi_stat );
    MPI_Get_count( &mpi_stat[1], mpi_type_particle, &mpi_count );
    count = mpi_count;
    *n_stored_particles += count;
  }
  
  if( *n_stored_particles >= start )
    {
      printf("Buffer too small!\n");
      exit(1);
    }
  
  /* Updating linked list */
  if( count > 0 )
    update_linked_list( ll, hoc, particles,*n_stored_particles-count, *n_stored_particles, 
			x, y, z, m, n, o, 
			m_start-ghosts, m_end+ghosts,
			n_start-ghosts, n_end+ghosts,
			o_start-ghosts, o_end+ghosts );
  MPI_Type_free(&mpi_type_particle);
}
