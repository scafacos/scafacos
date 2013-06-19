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



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "common/gridsort/gridsort.h"

#include "fcs_pp3mg.h"


#ifdef FCS_ENABLE_DEBUG
# define DEBUG_MOP(_mop_)  do { _mop_; } while (0)
#else
# define DEBUG_MOP(_mop_)  do {  } while (0)
#endif


FCSResult fcs_pp3mg_check(FCS handle)
{
  DEBUG_MOP(printf("fcs_pp3mg_check\n"));

  DEBUG_MOP(printf("fcs_pp3mg_check: done\n"));

	return NULL;
}


FCSResult fcs_pp3mg_init(FCS handle)
{
  FCSResult result;
  fcs_pp3mg_context_t *ctx;

  DEBUG_MOP(printf("fcs_pp3mg_init\n"));

  ctx = malloc(sizeof(fcs_pp3mg_context_t));
  ctx->data = malloc(sizeof(pp3mg_data));
  ctx->parameters = malloc(sizeof(pp3mg_parameters));

  fcs_set_method_context(handle, ctx);
  
  /* 
   * Default parameters:
   * pp3mg_cells_x = pp3mg_cells_y = pp3mg_cells_z = 128
   * h = 1/128 => h^4/(12*5*6)*max(f^(4)) = h^4/360*3840
   */
     
  result = fcs_pp3mg_setup(handle, 128, 128, 128, 6, 3, 10000, 50, 3.9736e-8, 1, 1);
  if (result != NULL) return result;

  DEBUG_MOP(printf("fcs_pp3mg_init: done\n"));

  return NULL;
}


FCSResult fcs_pp3mg_destroy(FCS handle)
{
  fcs_pp3mg_context_t *ctx;

  DEBUG_MOP(printf("fcs_pp3mg_destroy\n"));

  ctx = fcs_get_method_context(handle);

  pp3mg_free(ctx->data,ctx->parameters);

  free(ctx->data);
  free(ctx->parameters);
  free(ctx);

  fcs_set_method_context(handle, NULL);

  DEBUG_MOP(printf("fcs_pp3mg_destroy: done\n"));

  return NULL;
}


FCSResult fcs_pp3mg_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions, fcs_float *charges)
{
  MPI_Comm comm;
  fcs_pp3mg_context_t *ctx;
  fcs_float *box_a, *box_b, *box_c;
  fcs_float x, y, z;
  FCSResult result;

  DEBUG_MOP(printf("fcs_pp3mg_tune\n"));

  comm = fcs_get_communicator(handle);

  box_a = fcs_get_box_a(handle);
  box_b = fcs_get_box_b(handle);
  box_c = fcs_get_box_c(handle);
#define MAX( a, b ) ((a>b)?a:b)
  x = MAX(box_a[0],MAX(box_b[0],box_c[0]));
  y = MAX(box_a[1],MAX(box_b[1],box_c[1]));
  z = MAX(box_a[2],MAX(box_b[2],box_c[2]));

  result = fcs_pp3mg_set_max_particles(handle, local_max_particles);
  if (result != NULL) return result;

  ctx = fcs_get_method_context(handle);

  pp3mg_init(x,y,z,
	     handle->pp3mg_param->m,
	     handle->pp3mg_param->n,
	     handle->pp3mg_param->o,
	     handle->pp3mg_param->ghosts,
	     handle->pp3mg_param->degree,
	     handle->pp3mg_param->max_particles,
	     handle->pp3mg_param->maxiter,
	     handle->pp3mg_param->tol,
	     handle->pp3mg_param->distribution,
	     handle->pp3mg_param->discretization,
	     comm,
	     ctx->data,
	     ctx->parameters);

  DEBUG_MOP(printf("fcs_pp3mg_tune: done\n"));

  return NULL;
}


FCSResult fcs_pp3mg_run(FCS handle, fcs_int local_num_particles, fcs_int local_max_particles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials)
{
  fcs_pp3mg_context_t *ctx;

//  fcs_float e = 0.0;
//  fcs_float *p = NULL;
//  fcs_float *f = NULL;

  DEBUG_MOP(printf("fcs_pp3mg_run\n"));
  
  ctx = fcs_get_method_context(handle);
  
  ctx->last_runtime = MPI_Wtime();
  
  {
    fcs_int i;
    fcs_float *x, *y, *z, *fx, *fy, *fz;
    
    fcs_float box_a[3] = { ctx->parameters->x, 0, 0 };
    fcs_float box_b[3] = { 0, ctx->parameters->y, 0 };
    fcs_float box_c[3] = { 0, 0, ctx->parameters->z };
    fcs_float box_base[3] = { 0, 0, 0 };
    fcs_float lower_bound[3] = { ctx->parameters->x_start, ctx->parameters->y_start, ctx->parameters->z_start };
    fcs_float upper_bound[3] = { ctx->parameters->x_end, ctx->parameters->y_end, ctx->parameters->z_end };
    
    fcs_gridsort_t gridsort;
    
    fcs_int sorted_num_particles;
    fcs_float *sorted_charges, *sorted_positions;
    fcs_float *sorted_field, *sorted_potentials;
    fcs_gridsort_index_t *sorted_indices;
    
    fcs_gridsort_create(&gridsort);
    fcs_gridsort_set_system(&gridsort, box_base, box_a, box_b, box_c, NULL);
    fcs_gridsort_set_bounds(&gridsort, lower_bound, upper_bound);
    fcs_gridsort_set_particles(&gridsort, local_num_particles, local_max_particles, positions, charges);
    fcs_gridsort_sort_forward(&gridsort, 0.0, ctx->parameters->mpi_comm_cart);
    fcs_gridsort_get_real_particles(&gridsort, &sorted_num_particles, &sorted_positions, &sorted_charges, &sorted_indices);
    
    x = (fcs_float*) malloc(sorted_num_particles*sizeof(fcs_float));
    y = (fcs_float*) malloc(sorted_num_particles*sizeof(fcs_float));
    z = (fcs_float*) malloc(sorted_num_particles*sizeof(fcs_float));
    fx = (fcs_float*) malloc(sorted_num_particles*sizeof(fcs_float));
    fy = (fcs_float*) malloc(sorted_num_particles*sizeof(fcs_float));
    fz = (fcs_float*) malloc(sorted_num_particles*sizeof(fcs_float));
    sorted_field = (fcs_float*) malloc(3*sorted_num_particles*sizeof(fcs_float));
    sorted_potentials = (fcs_float*) malloc(sorted_num_particles*sizeof(fcs_float));
    
    for (i=0; i<sorted_num_particles; i++) {
      x[i] = sorted_positions[3*i+0];
      y[i] = sorted_positions[3*i+1];
      z[i] = sorted_positions[3*i+2];
    }
    
    pp3mg(x, y, z, sorted_charges, sorted_potentials, fx, fy, fz, sorted_num_particles, ctx->data, ctx->parameters);
    
    for (i=0; i<sorted_num_particles; i++) {
      sorted_field[3*i+0] = fx[i]/(1.0/(4.0*FCS_PI)*sorted_charges[i]);
      sorted_field[3*i+1] = fy[i]/(1.0/(4.0*FCS_PI)*sorted_charges[i]);
      sorted_field[3*i+2] = fz[i]/(1.0/(4.0*FCS_PI)*sorted_charges[i]);
      sorted_potentials[i] = sorted_potentials[i]/(1.0/(4.0*FCS_PI)*sorted_charges[i]);
    }
    
    fcs_gridsort_sort_backward(&gridsort, sorted_field, sorted_potentials, field, potentials, 1, ctx->parameters->mpi_comm_cart);
    
    if (x) free(x);
    if (y) free(y);
    if (z) free(z);
    if (fx) free(fx);
    if (fy) free(fy);
    if (fz) free(fz);
    if (sorted_field) free(sorted_field);
    if (sorted_potentials) free(sorted_potentials);
    
    fcs_gridsort_free(&gridsort);
    fcs_gridsort_destroy(&gridsort);
    
  }
  ctx->last_runtime = MPI_Wtime() - ctx->last_runtime;

  DEBUG_MOP(printf("fcs_pp3mg_run: done\n"));

  return NULL;
}


/* combined setter function for all pp3mg parameters */
extern FCSResult fcs_pp3mg_setup(FCS handle, fcs_int cells_x, fcs_int cells_y, fcs_int cells_z, fcs_int ghosts, fcs_int degree, fcs_int max_particles, fcs_int max_iterations, fcs_float tol, fcs_int distribution, fcs_int discretization)
{
  FCSResult result;

  result = fcs_pp3mg_set_cells_x(handle, cells_x);
  if (result != NULL) return result;

  result = fcs_pp3mg_set_cells_y(handle, cells_y);
  if (result != NULL) return result;

  result = fcs_pp3mg_set_cells_z(handle, cells_z);
  if (result != NULL) return result;

  result = fcs_pp3mg_set_ghosts(handle, ghosts);
  if (result != NULL) return result;

  result = fcs_pp3mg_set_degree(handle, degree);
  if (result != NULL) return result;

  result = fcs_pp3mg_set_max_particles(handle, max_particles);
  if (result != NULL) return result;

  result = fcs_pp3mg_set_max_iterations(handle, max_iterations);
  if (result != NULL) return result;

  result = fcs_pp3mg_set_tol(handle, tol);
  if (result != NULL) return result;

  result = fcs_pp3mg_set_distribution(handle, distribution);
  if (result != NULL) return result;

  result = fcs_pp3mg_set_discretization(handle, distribution);
  if (result != NULL) return result;

  return NULL;
}

/* setter for parameter cells_x */
FCSResult fcs_pp3mg_set_cells_x(FCS handle, fcs_int cells_x)
{
  handle->pp3mg_param->m = cells_x;
  
  return NULL;
}

/* getter for parameter cells_x */
FCSResult fcs_pp3mg_get_cells_x(FCS handle, fcs_int *cells_x)
{
  *cells_x = handle->pp3mg_param->m;
  
  return NULL;
}

/* setter for parameter cells_y */
FCSResult fcs_pp3mg_set_cells_y(FCS handle, fcs_int cells_y)
{
  handle->pp3mg_param->n = cells_y;
  
  return NULL;
}

/* getter for parameter cells_y */
FCSResult fcs_pp3mg_get_cells_y(FCS handle, fcs_int *cells_y)
{
  *cells_y = handle->pp3mg_param->n;
  
  return NULL;
}

/* setter for parameter cells_z */
FCSResult fcs_pp3mg_set_cells_z(FCS handle, fcs_int cells_z)
{
  handle->pp3mg_param->o = cells_z;
  
  return NULL;
}

/* getter for parameter cells_z */
FCSResult fcs_pp3mg_get_cells_z(FCS handle, fcs_int *cells_z)
{
  *cells_z = handle->pp3mg_param->o;
  
  return NULL;
}

/* setter for parameter ghosts */
FCSResult fcs_pp3mg_set_ghosts(FCS handle, fcs_int ghosts)
{
  handle->pp3mg_param->ghosts = ghosts;
  
  return NULL;
}

/* getter for parameter ghosts */
FCSResult fcs_pp3mg_get_ghosts(FCS handle, fcs_int *ghosts)
{
  *ghosts = handle->pp3mg_param->ghosts;
  
  return NULL;
}

/* setter for parameter degree */
FCSResult fcs_pp3mg_set_degree(FCS handle, fcs_int degree)
{
  handle->pp3mg_param->degree = degree;
  
  return NULL;
}

/* getter for parameter degree */
FCSResult fcs_pp3mg_get_degree(FCS handle, fcs_int *degree)
{
  *degree = handle->pp3mg_param->degree;
  
  return NULL;
}

/* setter for parameter max_particles */
FCSResult fcs_pp3mg_set_max_particles(FCS handle, fcs_int max_particles)
{
  handle->pp3mg_param->max_particles = max_particles;
  
  return NULL;
}

/* getter for parameter max_particles */
FCSResult fcs_pp3mg_get_max_particles(FCS handle, fcs_int *max_particles)
{
  *max_particles = handle->pp3mg_param->max_particles;
  
  return NULL;
}

/* setter for parameter maxiter */
FCSResult fcs_pp3mg_set_max_iterations(FCS handle, fcs_int max_iterations)
{
  handle->pp3mg_param->maxiter = max_iterations;
  
  return NULL;
}

/* getter for parameter maxiter */
FCSResult fcs_pp3mg_get_max_iterations(FCS handle, fcs_int *max_iterations)
{
  *max_iterations = handle->pp3mg_param->maxiter;
  
  return NULL;
}

/* setter for parameter tol */
FCSResult fcs_pp3mg_set_tol(FCS handle, fcs_float tol)
{
  handle->pp3mg_param->tol = tol;
  
  return NULL;
}

/* getter for parameter tol */
FCSResult fcs_pp3mg_get_tol(FCS handle, fcs_float *tol)
{
  *tol = handle->pp3mg_param->tol;
  
  return NULL;
}

/* setter for parameter distribution */
FCSResult fcs_pp3mg_set_distribution(FCS handle, fcs_int distribution)
{
  handle->pp3mg_param->distribution = distribution;
  
  return NULL;
}

/* getter for parameter distribution */
FCSResult fcs_pp3mg_get_distribution(FCS handle, fcs_int *distribution)
{
  *distribution = handle->pp3mg_param->distribution;
  
  return NULL;
}

/* setter for parameter discretization */
FCSResult fcs_pp3mg_set_discretization(FCS handle, fcs_int discretization)
{
  handle->pp3mg_param->discretization = discretization;
  
  return NULL;
}

/* getter for parameter discretization */
FCSResult fcs_pp3mg_get_discretization(FCS handle, fcs_int *discretization)
{
  *discretization = handle->pp3mg_param->discretization;
  
  return NULL;
}

/*
 *
 * Virial stuff
 *
 */
FCSResult fcs_pp3mg_require_virial(FCS handle, fcs_int compute_virial)
{
  return NULL;
}


FCSResult fcs_pp3mg_get_virial(FCS handle, fcs_float *virial)
{
  fcs_int i;

  for (i = 0; i < 9; ++i) virial[i] = 0.0;

  return NULL;
}
