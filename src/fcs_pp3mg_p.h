/*
  Copyright (C) 2011-2012 Rene Halver

  This file is part of ScaFaCoS.

  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.

  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/



#ifndef FCS_PP3MG_P_INCLUDED
#define FCS_PP3MG_P_INCLUDED


#include "fcs_definitions.h"
#include "fcs_result_p.h"
#include "fcs_interface_p.h"


/**
 * @file fcs_pp3mg_p.h
 * @brief file containing the method specific interface functions
 * for the pp3mg solver (public version)
 * @author Matthias Bolten
 */
typedef struct fcs_pp3mg_parameters_t *fcs_pp3mg_parameters;

/* setter for parameter cells_x */
FCSResult fcs_pp3mg_set_cells_x(FCS handle, fcs_int cells_x);

/* getter for parameter cells_x */
FCSResult fcs_pp3mg_get_cells_x(FCS handle, fcs_int *cells_x);

/* setter for parameter cells_y */
FCSResult fcs_pp3mg_set_cells_y(FCS handle, fcs_int cells_y);

/* getter for parameter cells_y */
FCSResult fcs_pp3mg_get_cells_y(FCS handle, fcs_int *cells_y);

/* setter for parameter cells_z */
FCSResult fcs_pp3mg_set_cells_z(FCS handle, fcs_int cells_z);

/* getter for parameter cells_z */
FCSResult fcs_pp3mg_get_cells_z(FCS handle, fcs_int *cells_z);

/* setter for parameter ghosts */
FCSResult fcs_pp3mg_set_ghosts(FCS handle, fcs_int ghosts);

/* getter for parameter ghosts */
FCSResult fcs_pp3mg_get_ghosts(FCS handle, fcs_int *ghosts);

/* setter for parameter degree */
FCSResult fcs_pp3mg_set_degree(FCS handle, fcs_int degree);

/* getter for parameter degree */
FCSResult fcs_pp3mg_get_degree(FCS handle, fcs_int *degree);

/* setter for parameter max_particles */
FCSResult fcs_pp3mg_set_max_particles(FCS handle, fcs_int max_particles);

/* getter for parameter max_particles */
FCSResult fcs_pp3mg_get_max_particles(FCS handle, fcs_int *degree);

/* setter for parameter maxiter */
FCSResult fcs_pp3mg_set_max_iterations(FCS handle, fcs_int max_iterations);

/* getter for parameter maxiter */
FCSResult fcs_pp3mg_get_max_iterations(FCS handle, fcs_int *max_iterations);

/* setter for parameter tol */
FCSResult fcs_pp3mg_set_tol(FCS handle, fcs_float tol);

/* getter for parameter tol */
FCSResult fcs_pp3mg_get_tol(FCS handle, fcs_float *tol);

/* setter for parameter distribution */
FCSResult fcs_pp3mg_set_distribution(FCS handle, fcs_int disctribution);

/* getter for parameter distribution */
FCSResult fcs_pp3mg_get_distribution(FCS handle, fcs_int *distribution);

/* setter for parameter discretization */
FCSResult fcs_pp3mg_set_discretization(FCS handle, fcs_int discretization);

/* getter for parameter discretization */
FCSResult fcs_pp3mg_get_discretization(FCS handle, fcs_int *discertization);

/* combined setter for all solver parameters */
FCSResult fcs_pp3mg_setup(FCS handle, fcs_int cells_x, fcs_int cells_y, fcs_int cells_z, fcs_int ghosts, fcs_int degree, fcs_int max_particles, fcs_int maxiter, fcs_float tol, fcs_int distribution, fcs_int discretization);

#endif
