/*
 * pp3mg.h
 *
 * This file contains definitions of functions found in the file
 * "pp3mg.c".
 *
 * Based on Fortran code by Matthias Bolten.
 *
 * Authors: Matthias Bolten, Stephanie Friedhoff
 *
 * Copyright 2009, 2010, 2011, 2012 Matthias Bolten, Stephanie Friedhoff
 * All rights reserved.
 *
 */

#ifndef _PP3MG__H_
#define _PP3MG__H_

/*
 *
 * Typedefs 
 *
 */

enum CHARGE_DISTRIBUTION {
  polynomial_deg_6,
  polynomial_deg_10,
  polynomial_deg_14,
  spline_deg_4
};

enum DISCRETIZATION {
  order_2,
  order_4_compact,
  order_6
};

/* Particle */
typedef struct{
  double x;
  double y;
  double z;
  double q;
  double e;
  double fx;
  double fy;
  double fz;
} pp3mg_particle;

/* Parameters */
typedef struct{
  /* Communicator for cartesian process grid */
  MPI_Comm mpi_comm_cart;
	
  /* Physical dimension */
  double x, y, z;

  /* Grid dimension */
  int m, n, o;
	
  /* Boundaries of local part of the grid */
  int m_start, m_end;
  int n_start, n_end;
  int o_start, o_end;

  /* Boundaries of local part of the domain */
  double x_start, x_end;
  double y_start, y_end;
  double z_start, z_end;

  /* Grid spacings */
  double hx, hy, hz;
	
  /* Radius of replacement charge distribution */
  double radius;
  /* various integer powers of the radius */
  double radius_2, radius_3, radius_4, radius_6, radius_7;

  /* Number of ghosts in each direction */
  int ghosts;
	
  /* Degree of interpolation polynomial */
  int degree;

  /* Maximum number of iterations */
  int maxiter;

  /* Relative residiual for solver */
  double tol;

  /* Charge distribution */
  enum CHARGE_DISTRIBUTION distribution;

  /* Discretization */
  enum DISCRETIZATION discretization;
} pp3mg_parameters;

/* Data */
typedef struct{
  /* Maximum number of particles to store on one processor */
  int max_particles;

  /* Total number of particles */
  int n_particles;

  /* Number of local particles (not including particles in ghost cells) */
  int n_local_particles;

  /* Number of locally stored particles (including particles in ghost cells) */
  int n_stored_particles;

  /* Variables for particle store */
  pp3mg_particle* particles;

  /* Linked list */
  int*** hoc;
  int* ll;

  /* Grids */
  double*** f;
  double*** f_ghosted;
  double*** u;
  double*** u_ghosted;
  /*
  double*** ux_ghosted;
  double*** uy_ghosted;
  double*** uz_ghosted;
  */
} pp3mg_data;

/*
 *
 * Functions
 *
 */
void pp3mg_init( double x_in, double y_in, double z_in, int m_in, int n_in,
		 int o_in, int ghosts_in, int degree_in, int max_particles_in, 
		 int maxiter_in, double tol_in,  
		 enum CHARGE_DISTRIBUTION distribution_in, 
		 enum DISCRETIZATION discretization_in, MPI_Comm mpi_comm,
		 pp3mg_data* data, pp3mg_parameters* params);

void pp3mg( double* x, double* y, double* z, double* q, double* e, 
	    double* fx, double* fy, double* fz, int n_local_particles_in,
	    pp3mg_data* data, pp3mg_parameters* params );

void pp3mg_free( pp3mg_data* data, pp3mg_parameters* params );


#endif  /* ifndef _PP3MG__H_ */
