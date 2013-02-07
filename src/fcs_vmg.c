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

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <mpi.h>

#include "fcs_vmg.h"
#include "FCSCommon.h"

FCSResult fcs_vmg_init(FCS handle)
{

  return NULL;
}

FCSResult fcs_vmg_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles,
		       fcs_float* positions, fcs_float* charges)
{
  FCSResult result;
  /*
   * Set default parameters if no parameters were specified
   */
  result = fcs_vmg_set_default(handle);
  if (result)
    return result;

  /*
   * Check parameters
   */
  result = fcs_vmg_check(handle);
  if (result)
    return result;

  /*
   * Setup vmg solver
   */
  fcs_int level;
  fcs_int max_iter;
  fcs_int smoothing_steps;
  fcs_int cycle_type;
  fcs_float precision;
  fcs_int near_field_cells;
  fcs_int interpolation_order;
  fcs_int discretization_order;

  result = fcs_vmg_get_max_level(handle, &level);
  if (result)
    return result;

  fcs_int* periodic = fcs_get_periodicity(handle);

  result  = fcs_vmg_get_max_iterations(handle, &max_iter);
  if (result)
    return result;

  result  = fcs_vmg_get_smoothing_steps(handle, &smoothing_steps);
  if (result)
    return result;

  result  = fcs_vmg_get_cycle_type(handle, &cycle_type);
  if (result)
    return result;

  result  = fcs_vmg_get_precision(handle, &precision);
  if (result)
    return result;

  fcs_float *offset = fcs_get_offset(handle);
  fcs_float* box_a = fcs_get_box_a(handle);

  result  = fcs_vmg_get_near_field_cells(handle, &near_field_cells);
  if (result)
    return result;

  result  = fcs_vmg_get_interpolation_order(handle, &interpolation_order);
  if (result)
    return result;

  result  = fcs_vmg_get_discretization_order(handle, &discretization_order);
  if (result)
    return result;

  MPI_Comm comm = fcs_get_communicator(handle);

  VMG_fcs_setup(level, periodic, max_iter, smoothing_steps,
		cycle_type, precision, offset, box_a[0],
		near_field_cells, interpolation_order,
		discretization_order, comm);

  result = fcs_vmg_library_check(handle);
  if (result)
    return result;

  return NULL;
}

FCSResult fcs_vmg_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles,
		      fcs_float* positions, fcs_float* charges, fcs_float *field,
		      fcs_float *potentials)
{
  VMG_fcs_run(positions, charges, potentials, field, local_particles);

  return NULL;
}

FCSResult fcs_vmg_destroy(FCS handle)
{
  VMG_fcs_destroy();
  return NULL;
}

/**
 * @brief Function to set all vmg parameters with a single call
 *
 * @param handle FCS-object that contains the parameters
 * @param max_level The maximum level of the algorithm, i.e. n_gridpoints = 2^max_level.
 * @param spline_degree The degree of the interpolating B-Splines.
 * @param max_iterations The maximum number of multigrid iterations.
 * @param smoothing_steps Number of pre/postsmoothing steps on each level.
 * @param cycle_type Cycle-number.
 * @param precision Desired precision.
 * @param near_field_cells Number of near field cells.
 * @param interpolation_order Interpolation order.
 * @param discretization_order Discretization order.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_setup(FCS handle, fcs_int max_level,
			       fcs_int max_iterations, fcs_int smoothing_steps,
			       fcs_int cycle_type, fcs_float precision, fcs_int near_field_cells,
                               fcs_int interpolation_order, fcs_int discretization_order)
{
  FCSResult result;
  char fnc_name[] = "fcs_vmg_setup";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  result = fcs_vmg_set_max_level(handle, max_level);
  if (result != NULL)
    return result;

  result = fcs_vmg_set_max_iterations(handle, max_iterations);
  if (result != NULL)
    return result;

  result = fcs_vmg_set_smoothing_steps(handle, smoothing_steps);
  if (result != NULL)
    return result;

  result = fcs_vmg_set_cycle_type(handle, cycle_type);
  if (result != NULL)
    return result;

  result = fcs_vmg_set_precision(handle, precision);
  if (result != NULL)
    return result;

  result = fcs_vmg_set_near_field_cells(handle, near_field_cells);
  if (result != NULL)
    return result;

  result = fcs_vmg_set_interpolation_order(handle, interpolation_order);
  if (result != NULL)
    return result;

  result = fcs_vmg_set_discretization_order(handle, discretization_order);
  if (result != NULL)
    return result;

  return NULL;
}

/**
 * @brief function to set default values if no values are set
 * @param handle the FCS-obect into which the default parameters will be entered
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_vmg_set_default(FCS handle)
{
  MPI_Comm comm;
  int rank;

  const char fnc_name[] = "fcs_vmg_set_default";

  comm = fcs_get_communicator(handle);
  MPI_Comm_rank(comm, &rank);

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  fcs_int max_level;
  fcs_vmg_get_max_level(handle, &max_level);
  if (max_level < 0) {
    max_level = 6;
    if (rank == 0)
      printf("Warning: %s: Parameter %s not set. Set default to %d.\n", fnc_name, "max_level", max_level);
    fcs_vmg_set_max_level(handle, max_level);
  }

  fcs_int max_iter;
  fcs_vmg_get_max_iterations(handle, &max_iter);
  if (max_iter < 0) {
    max_iter = 15;
    if (rank == 0)
      printf("Warning: %s: Parameter %s not set. Set default to %d.\n", fnc_name, "max_iterations", max_iter);
    fcs_vmg_set_max_iterations(handle, max_iter);
  }

  fcs_int smoothing_steps;
  fcs_vmg_get_smoothing_steps(handle, &smoothing_steps);
  if (smoothing_steps < 0) {
    smoothing_steps = 3;
    if (rank == 0)
      printf("Warning: %s: Parameter %s not set. Set default to %d.\n", fnc_name, "smoothing_steps", smoothing_steps);
    fcs_vmg_set_smoothing_steps(handle, smoothing_steps);
  }

  fcs_int cycle_type;
  fcs_vmg_get_cycle_type(handle, &cycle_type);
  if (cycle_type < 0) {
    cycle_type = 1;
    if (rank == 0)
      printf("Warning: %s: Parameter %s not set. Set default to %d.\n", fnc_name, "cycle_type", cycle_type);
    fcs_vmg_set_cycle_type(handle, cycle_type);
  }

  fcs_float precision;
  fcs_vmg_get_precision(handle, &precision);
  if (precision < 0.0) {
    precision = 1.0e-8;
    if (rank == 0)
      printf("Warning: %s: Parameter %s not set. Set default to %e.\n", fnc_name, "precision", precision);
    fcs_vmg_set_precision(handle, precision);
  }

  fcs_int near_field_cells;
  fcs_vmg_get_near_field_cells(handle, &near_field_cells);
  if (near_field_cells < 0) {
    near_field_cells = 4;
    if (rank == 0)
      printf("Warning: %s: Parameter %s not set. Set default to %d.\n", fnc_name, "near_field_cells", near_field_cells);
    fcs_vmg_set_near_field_cells(handle, near_field_cells);
  }

  fcs_int interpolation_order;
  fcs_vmg_get_interpolation_order(handle, &interpolation_order);
  if (interpolation_order < 0) {
    interpolation_order = 5;
    if (rank == 0)
      printf("Warning: %s: Parameter %s not set. Set default to %d.\n", fnc_name, "interpolation_order", interpolation_order);
    fcs_vmg_set_interpolation_order(handle, interpolation_order);
  }

  fcs_int discretization_order;
  fcs_vmg_get_discretization_order(handle, &discretization_order);
  if (discretization_order < 0) {
    discretization_order = 4;
    if (rank == 0)
      printf("Warning: %s: Parameter %s not set. Set default to %d.\n", fnc_name, "discretization_order", discretization_order);
    fcs_vmg_set_discretization_order(handle, discretization_order);
  }

  return NULL;
}

/**
 * @brief Set the maximum level of the multigrid algorithm.
 * This sets the number of grid points to 2^max_level, so basically this
 * parameter ensures that the number of grid points will be a power of two.
 *
 * @param handle FCS-object that contains the parameter
 * @param max_level number of iterations
 *
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_vmg_set_max_level(FCS handle, fcs_int max_level)
{
  char fnc_name[] = "fcs_vmg_set_max_level";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  if (max_level < 3)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "The finest level must be at least 3.");

  handle->vmg_param->max_level = max_level;

  return NULL;
}

/**
 * @brief Get the maximum level of the multigrid algorithm.
 * The number of grid points of the finest multigrid level
 * will be 2^max_level.
 *
 * @param handle FCS-object that contains the parameter
 * @param max_level Number of iterations.
 *
 * @return FCSResult-object containing the return state.
 */
extern FCSResult fcs_vmg_get_max_level(FCS handle, fcs_int *max_level)
{
  char fnc_name[] = "fcs_vmg_get_max_level";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  *max_level = handle->vmg_param->max_level;

  return NULL;
}

/**
 * @brief Set the maximum number of multigrid iterations.
 *
 * @param handle FCS-object that contains the parameter
 * @param max_iterations Number of iterations.
 *
 * @return FCSResult-object containing the return state.
 */
extern FCSResult fcs_vmg_set_max_iterations(FCS handle, fcs_int max_iterations)
{
  char fnc_name[] = "fcs_vmg_set_max_iterations";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  if (max_iterations < 1)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "The maximum number of iterations must be positive.");

  handle->vmg_param->max_iterations = max_iterations;

  return NULL;
}

/**
 * @brief Get the maximum number of multigrid iterations.
 *
 * @param handle FCS-object that contains the parameter.
 * @param max_iterations Number of iterations.
 *
 * @return FCSResult-object containing the return state.
 */
extern FCSResult fcs_vmg_get_max_iterations(FCS handle, fcs_int *max_iterations)
{
  char fnc_name[] = "fcs_vmg_get_max_iterations";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  *max_iterations = handle->vmg_param->max_iterations;

  return NULL;
}

/**
 * @brief Set the number of pre/postsmoothing steps on each level
 *
 * @param handle FCS-object that contains the parameter
 * @param smoothing_steps Number of smoothing steps
 *
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_vmg_set_smoothing_steps(FCS handle, fcs_int smoothing_steps)
{
  char fnc_name[] = "fcs_vmg_set_smoothing_steps";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  if (smoothing_steps < 1)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "The number of smoothing steps must be at least one.");

  handle->vmg_param->smoothing_steps = smoothing_steps;

  return NULL;
}

/**
 * @brief Get the number of pre/postsmoothing steps on each level
 *
 * @param handle FCS-object that contains the parameter
 * @param smoothing_steps Number of pre/postsmoothing steps
 *
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_vmg_get_smoothing_steps(FCS handle, fcs_int *smoothing_steps)
{
  char fnc_name[] = "fcs_vmg_get_smoothing_steps";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  *smoothing_steps = handle->vmg_param->smoothing_steps;

  return NULL;
}

/**
 * @brief Set the cycle-number of the multigrid cycle used. E.g. 1 corresponds
 *        to a V-Cycle and 2 corresponds to a W-Cycle.
 *
 * @param handle FCS-object that contains the parameter.
 * @param cycle_type Gamma.
 *
 * @return FCSReturn-object containing the return state.
 */
extern FCSResult fcs_vmg_set_cycle_type(FCS handle, fcs_int cycle_type)
{
  char fnc_name[] = "fcs_vmg_set_cycle_type";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  if (cycle_type < 1)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Gamma must be at least one.");

  handle->vmg_param->cycle_type = cycle_type;

  return NULL;
}

/**
 * @brief Get the cycle_type-number of the multigrid cycle used. E.g. 1 corresponds
 *        to a V-Cycle and 2 corresponds to a W-Cycle.
 *
 * @param handle FCS-object that contains the parameter.
 * @param cycle_type type of multigrid cycle.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_get_cycle_type(FCS handle, fcs_int* cycle_type)
{
  char fnc_name[] = "fcs_vmg_get_cycle_type";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  *cycle_type = handle->vmg_param->cycle_type;

  return NULL;
}

/**
 * @brief Set the precision of the vmg algorithm. Currently, the relative residual
 *        computed in the discrete L2 norm will be tested against this value.
 *
 * @param handle FCS-object that contains the parameter.
 * @param precision Precision.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_set_precision(FCS handle, fcs_float precision)
{
  char fnc_name[] = "fcs_vmg_set_precision";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  if (precision < 0.0)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "The requested precision must be positive.");

  handle->vmg_param->precision = precision;

  return NULL;
}

/**
 * @brief Get the precision of the vmg algorithm. Currently, the relative residual
 *        computed in the discrete L2 norm will be tested against this value.
 *
 * @param handle FCS-object that contains the parameter.
 * @param precision Precision of the vmg algorithm.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_get_precision(FCS handle, fcs_float* precision)
{
  char fnc_name[] = "fcs_vmg_get_precision";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  *precision = handle->vmg_param->precision;

  return NULL;
}

/**
 * @brief Set the number of near field cells for separating the
 *        near/far field part of the potential.
 *
 * @param handle FCS-object that contains the parameter.
 * @param near_field_cells Near field cells.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_set_near_field_cells(FCS handle, fcs_int near_field_cells)
{
  char fnc_name[] = "fcs_vmg_set_near_field_cells";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  if (near_field_cells < 1)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Number of near field cells  must be positive.");

  handle->vmg_param->near_field_cells = near_field_cells;

  return NULL;
}

/**
 * @brief Get the number of near field cells for separating the
 *        near/far field part of the potential.
 *
 * @param handle FCS-object that contains the parameter.
 * @param near_field_cells Number of near field cells.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_get_near_field_cells(FCS handle, fcs_int* near_field_cells)
{
  char fnc_name[] = "fcs_vmg_get_near_field_cells";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  *near_field_cells = handle->vmg_param->near_field_cells;

  return NULL;
}

/**
 * @brief Set the interpolation order for interpolating
 *        the gridded potential to the particle positions.
 *
 * @param handle FCS-object that contains the parameter.
 * @param interpolation order Interpolation order.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_set_interpolation_order(FCS handle, fcs_int interpolation_order)
{
  char fnc_name[] = "fcs_vmg_set_interpolation_order";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  if (interpolation_order < 3)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Interpolation order must be greater or equal three.");

  handle->vmg_param->interpolation_order = interpolation_order;

  return NULL;
}

/**
 * @brief Get the interpolation order for interpolating
 *        the gridded potential to the particle positions.
 *
 * @param handle FCS-object that contains the parameter.
 * @param interpolation_order Interpolation order.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_get_interpolation_order(FCS handle, fcs_int* interpolation_order)
{
  char fnc_name[] = "fcs_vmg_get_interpolation_order";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  *interpolation_order = handle->vmg_param->interpolation_order;

  return NULL;
}

/**
 * @brief Set the discretization order.
 *        This will affect the discretization error heavily, but a
 *        higher value also implies more computation. Possible values
 *        are 2 or 4.
 *
 * @param handle FCS-object that contains the parameter.
 * @param discretization_order Discretization order.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_set_discretization_order(FCS handle, fcs_int discretization_order)
{
  char fnc_name[] = "fcs_vmg_set_discretization_order";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  if (discretization_order != 2 && discretization_order != 4)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "discretization_order must be 2 or 4.");

  handle->vmg_param->discretization_order = discretization_order;

  return NULL;
}

/**
 * @brief Get the discretization order.
 *        This will affect the discretization error heavily, but a
 *        higher value also implies more computation. Possible values
 *        are 2 or 4.
 *
 * @param handle FCS-object that contains the parameter.
 * @param discretization_order Discretization order.
 *
 * @return FCSResult-object containing the return value.
 */
extern FCSResult fcs_vmg_get_discretization_order(FCS handle, fcs_int* discretization_order)
{
  char fnc_name[] = "fcs_vmg_get_discretization_order";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_VMG)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"vmg\".");

  *discretization_order = handle->vmg_param->discretization_order;

  return NULL;
}

extern FCSResult fcs_vmg_check(FCS handle)
{
  char fnc_name[] = "fcs_vmg_check";
  FCSResult result;

  fcs_int max_level;
  result = fcs_vmg_get_max_level(handle, &max_level);
  if (result != NULL)
    return result;
  if (max_level == -1)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "vmg finest level not set.");
  if (max_level < 3)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg finest level must be greater than 2.");

  fcs_int max_iterations;
  result = fcs_vmg_get_max_iterations(handle, &max_iterations);
  if (result != NULL)
    return result;
  if (max_iterations == -1)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "vmg maximum number of iterations not set.");
  if (max_iterations < 1)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg maximum number of iterations must be positive.");

  fcs_int smoothing_steps;
  result = fcs_vmg_get_smoothing_steps(handle, &smoothing_steps);
  if (result != NULL)
    return result;
  if (smoothing_steps == -1)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "vmg number of smoothing steps not set.");
  if (smoothing_steps < 1)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg number of smoothing steps must be positive.");

  fcs_int cycle_type;
  result = fcs_vmg_get_cycle_type(handle, &cycle_type);
  if (result != NULL)
    return result;
  if (cycle_type == -1)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "vmg cycle number not set.");
  if (cycle_type < 1)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg cycle number must be positive.");

  fcs_float precision;
  result = fcs_vmg_get_precision(handle, &precision);
  if (result != NULL)
    return result;
  if (precision == -1.)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "vmg desired precision not set.");

  fcs_int near_field_cells;
  result = fcs_vmg_get_near_field_cells(handle, &near_field_cells);
  if (result != NULL)
    return result;
  if (near_field_cells == -1)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "vmg number of near field cells not set.");
  if (near_field_cells < 1)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg number of near field cells must be positive.");

  fcs_int interpolation_order;
  result = fcs_vmg_get_interpolation_order(handle, &interpolation_order);
  if (result != NULL)
    return result;
  if (interpolation_order == -1)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "vmg interpolation order not set.");
  if (interpolation_order < 3)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg interpolation order must be greater or equal three.");

  fcs_int discretization_order;
  result = fcs_vmg_get_discretization_order(handle, &discretization_order);
  if (result != NULL)
    return result;
  if (discretization_order == -1)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "vmg discretization order not set.");
  if (discretization_order != 2 && discretization_order != 4)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg discretization order must be 2 or 4.");

  fcs_float* box_a = fcs_get_box_a(handle);
  fcs_float* box_b = fcs_get_box_b(handle);
  fcs_float* box_c = fcs_get_box_c(handle);

  if (!fcs_is_cubic(box_a, box_b, box_c))
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg requires the box to be cubic.");

  if (!fcs_uses_principal_axes(box_a, box_b, box_c))
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg requires the box vectors to be parallel to the principal axes.");

  return NULL;
}

FCSResult fcs_vmg_library_check(FCS handle)
{
  int rval = VMG_fcs_check();
  char fnc_name[] = "fcs_vmg_library_check";

  switch (rval)
    {
    case 0:
      return NULL;
      break;
    case 1:
      return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg requires the number of grid points on each process in every direction to be greater than the number of near field cells.");
      break;
    default:
      return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "vmg unknown error code");
      break;
    }
}

FCSResult fcs_vmg_require_virial(FCS handle, fcs_int flag)
{
  return NULL;
}

FCSResult fcs_vmg_get_virial(FCS handle, fcs_float *virial)
{
  char fnc_name[] =  "fcs_vmg_get_virial";

  if (!handle || !handle->vmg_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name, "null pointer supplied as handle");
  if (!virial)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name, "null pointer supplied for virial");

  fcs_int i;
  for (i=0; i < 9; i++)
    virial[i] = 0.0;

  return NULL;
}
