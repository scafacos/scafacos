/*
 * interpolation.h
 *
 * This file contains definitions of functions found in the file
 * "interpolation.c".
 *
 * Based on Fortran code by Matthias Bolten.
 *
 * Authors: Matthias Bolten, Stephanie Friedhoff
 * Created: 2009/06/02
 *
 * Copyright 2009, 2010, 2011, 2012 Matthias Bolten, Stephanie Friedhoff
 * All rights reserved.
 *
 */

#ifndef _INTERPOLATION__H_
#define _INTERPOLATION__H_

void interp_poly( int degree, pp3mg_particle* particles, int*** hoc, int* ll, double*** u, 
		  int n_particles, int m_start, int m_end, int n_start, int n_end, 
		  int o_start, int o_end, int boundary, double hx, double hy, double hz );

#endif  /* ifndef _INTERPOLATION__H_ */
