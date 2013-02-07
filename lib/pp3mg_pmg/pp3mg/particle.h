/*
 * particle.h
 *
 * This file contains definitions of functions found in the file
 * "particle.c".
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


#ifndef _PARTICLE__H_
#define _PARTICLE__H_

void update_linked_list( int* ll, int*** hoc, pp3mg_particle* particles, 
			 int start, int end, double x, double y, double z, int m, int n, int o,
			 int m_start, int m_end, int n_start, int n_end, int o_start, int o_end );

void update_particle_ghosts( int** ll_addr, int*** hoc, pp3mg_particle** particles_addr,
			     double x, double y, double z, int m, int n, int o,
			     int ghosts, int* n_stored_particles, int *max_particles,
			     int m_start, int m_end, int n_start, int n_end, 
			     int o_start, int o_end, MPI_Comm mpi_comm_cart );

#endif  /* ifndef _PARTICLE__H_ */
