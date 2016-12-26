/*
  Copyright (C) 2011-2012 Rene Halver
  Copyright (C) 2016 Michael Hofmann

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


#define VMG_CHECK_RETURN_RESULT(_h_, _f_)  do { \
  CHECK_HANDLE_RETURN_RESULT(_h_, _f_); \
  CHECK_METHOD_RETURN_RESULT(_h_, _f_, FCS_METHOD_VMG, "vmg"); \
  } while (0)

#define VMG_CHECK_RETURN_VAL(_h_, _f_, _v_)  do { \
  CHECK_HANDLE_RETURN_VAL(_h_, _f_, _v_); \
  CHECK_METHOD_RETURN_VAL(_h_, _f_, FCS_METHOD_VMG, "vmg", _v_); \
  } while (0)


FCSResult fcs_vmg_init(FCS handle)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  handle->shift_positions = 0;

  handle->destroy = fcs_vmg_destroy;
  handle->set_parameter = fcs_vmg_set_parameter;
  handle->print_parameters = fcs_vmg_print_parameters;
  handle->tune = fcs_vmg_tune;
  handle->run = fcs_vmg_run;

  handle->vmg_param = malloc(sizeof(*handle->vmg_param));
  handle->vmg_param->max_level = -1;
  handle->vmg_param->max_iterations = -1;
  handle->vmg_param->smoothing_steps = -1;
  handle->vmg_param->cycle_type = -1;
  handle->vmg_param->precision = -1.0;
  handle->vmg_param->near_field_cells = -1;
  handle->vmg_param->interpolation_order = -1;
  handle->vmg_param->discretization_order = -1;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_vmg_tune(FCS handle, fcs_int local_particles,
		       fcs_float* positions, fcs_float* charges)
{
  FCSResult result;

  VMG_CHECK_RETURN_RESULT(handle, __func__);

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

  const fcs_int* periodic = fcs_get_periodicity(handle);

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

  const fcs_float *offset = fcs_get_box_origin(handle);
  const fcs_float* box_a = fcs_get_box_a(handle);

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

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_vmg_run(FCS handle, fcs_int local_particles,
		      fcs_float* positions, fcs_float* charges, fcs_float *field,
		      fcs_float *potentials)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  VMG_fcs_run(positions, charges, potentials, field, local_particles);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_vmg_destroy(FCS handle)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  VMG_fcs_destroy();

  free(handle->vmg_param);

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_setup(FCS handle, fcs_int max_level,
			       fcs_int max_iterations, fcs_int smoothing_steps,
			       fcs_int cycle_type, fcs_float precision, fcs_int near_field_cells,
                               fcs_int interpolation_order, fcs_int discretization_order)
{
  FCSResult result;

  VMG_CHECK_RETURN_RESULT(handle, __func__);

  result = fcs_vmg_set_max_level(handle, max_level);
  CHECK_RESULT_RETURN(result);

  result = fcs_vmg_set_max_iterations(handle, max_iterations);
  CHECK_RESULT_RETURN(result);

  result = fcs_vmg_set_smoothing_steps(handle, smoothing_steps);
  CHECK_RESULT_RETURN(result);

  result = fcs_vmg_set_cycle_type(handle, cycle_type);
  CHECK_RESULT_RETURN(result);

  result = fcs_vmg_set_precision(handle, precision);
  CHECK_RESULT_RETURN(result);

  result = fcs_vmg_set_near_field_cells(handle, near_field_cells);
  CHECK_RESULT_RETURN(result);

  result = fcs_vmg_set_interpolation_order(handle, interpolation_order);
  CHECK_RESULT_RETURN(result);

  result = fcs_vmg_set_discretization_order(handle, discretization_order);
  CHECK_RESULT_RETURN(result);

  return FCS_RESULT_SUCCESS;
}

/**
 * @brief function to set default values if no values are set
 * @param handle the FCS-obect into which the default parameters will be entered
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_set_default(FCS handle)
{
  MPI_Comm comm;
  int rank;

  VMG_CHECK_RETURN_RESULT(handle, __func__);

  comm = fcs_get_communicator(handle);
  MPI_Comm_rank(comm, &rank);

  fcs_int max_level;
  fcs_vmg_get_max_level(handle, &max_level);
  if (max_level < 0) {
    max_level = 6;
#ifdef FCS_ENABLE_DEBUG
    if (rank == 0)
      printf("%s: Parameter %s not set. Set default to %d.\n", __func__, "max_level", max_level);
#endif
    fcs_vmg_set_max_level(handle, max_level);
  }

  fcs_int max_iter;
  fcs_vmg_get_max_iterations(handle, &max_iter);
  if (max_iter < 0) {
    max_iter = 15;
#ifdef FCS_ENABLE_DEBUG
    if (rank == 0)
      printf("%s: Parameter %s not set. Set default to %d.\n", __func__, "max_iterations", max_iter);
#endif
    fcs_vmg_set_max_iterations(handle, max_iter);
  }

  fcs_int smoothing_steps;
  fcs_vmg_get_smoothing_steps(handle, &smoothing_steps);
  if (smoothing_steps < 0) {
    smoothing_steps = 3;
#ifdef FCS_ENABLE_DEBUG
    if (rank == 0)
      printf("%s: Parameter %s not set. Set default to %d.\n", __func__, "smoothing_steps", smoothing_steps);
#endif
    fcs_vmg_set_smoothing_steps(handle, smoothing_steps);
  }

  fcs_int cycle_type;
  fcs_vmg_get_cycle_type(handle, &cycle_type);
  if (cycle_type < 0) {
    cycle_type = 1;
#ifdef FCS_ENABLE_DEBUG
    if (rank == 0)
      printf("%s: Parameter %s not set. Set default to %d.\n", __func__, "cycle_type", cycle_type);
#endif
    fcs_vmg_set_cycle_type(handle, cycle_type);
  }

  fcs_float precision;
  fcs_vmg_get_precision(handle, &precision);
  if (precision < 0.0) {
    precision = 1.0e-8;
#ifdef FCS_ENABLE_DEBUG
    if (rank == 0)
      printf("%s: Parameter %s not set. Set default to %e.\n", __func__, "precision", precision);
#endif
    fcs_vmg_set_precision(handle, precision);
  }

  fcs_int near_field_cells;
  fcs_vmg_get_near_field_cells(handle, &near_field_cells);
  if (near_field_cells < 0) {
    near_field_cells = 4;
#ifdef FCS_ENABLE_DEBUG
    if (rank == 0)
      printf("%s: Parameter %s not set. Set default to %d.\n", __func__, "near_field_cells", near_field_cells);
#endif
    fcs_vmg_set_near_field_cells(handle, near_field_cells);
  }

  fcs_int interpolation_order;
  fcs_vmg_get_interpolation_order(handle, &interpolation_order);
  if (interpolation_order < 0) {
    interpolation_order = 5;
#ifdef FCS_ENABLE_DEBUG
    if (rank == 0)
      printf("%s: Parameter %s not set. Set default to %d.\n", __func__, "interpolation_order", interpolation_order);
#endif
    fcs_vmg_set_interpolation_order(handle, interpolation_order);
  }

  fcs_int discretization_order;
  fcs_vmg_get_discretization_order(handle, &discretization_order);
  if (discretization_order < 0) {
    discretization_order = 4;
#ifdef FCS_ENABLE_DEBUG
    if (rank == 0)
      printf("%s: Parameter %s not set. Set default to %d.\n", __func__, "discretization_order", discretization_order);
#endif
    fcs_vmg_set_discretization_order(handle, discretization_order);
  }

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_set_max_level(FCS handle, fcs_int max_level)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  if (max_level < 3)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "The finest level must be at least 3.");

  handle->vmg_param->max_level = max_level;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_get_max_level(FCS handle, fcs_int *max_level)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  *max_level = handle->vmg_param->max_level;

  return FCS_RESULT_SUCCESS;
}

/**
 * @brief Set the maximum number of multigrid iterations.
 *
 * @param handle FCS-object that contains the parameter
 * @param max_iterations Number of iterations.
 *
 * @return FCSResult-object containing the return state.
 */
FCSResult fcs_vmg_set_max_iterations(FCS handle, fcs_int max_iterations)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  if (max_iterations < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "The maximum number of iterations must be positive.");

  handle->vmg_param->max_iterations = max_iterations;

  return FCS_RESULT_SUCCESS;
}

/**
 * @brief Get the maximum number of multigrid iterations.
 *
 * @param handle FCS-object that contains the parameter.
 * @param max_iterations Number of iterations.
 *
 * @return FCSResult-object containing the return state.
 */
FCSResult fcs_vmg_get_max_iterations(FCS handle, fcs_int *max_iterations)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  *max_iterations = handle->vmg_param->max_iterations;

  return FCS_RESULT_SUCCESS;
}

/**
 * @brief Set the number of pre/postsmoothing steps on each level
 *
 * @param handle FCS-object that contains the parameter
 * @param smoothing_steps Number of smoothing steps
 *
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_set_smoothing_steps(FCS handle, fcs_int smoothing_steps)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  if (smoothing_steps < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "The number of smoothing steps must be at least one.");

  handle->vmg_param->smoothing_steps = smoothing_steps;

  return FCS_RESULT_SUCCESS;
}

/**
 * @brief Get the number of pre/postsmoothing steps on each level
 *
 * @param handle FCS-object that contains the parameter
 * @param smoothing_steps Number of pre/postsmoothing steps
 *
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_get_smoothing_steps(FCS handle, fcs_int *smoothing_steps)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  *smoothing_steps = handle->vmg_param->smoothing_steps;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_set_cycle_type(FCS handle, fcs_int cycle_type)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  if (cycle_type < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "Gamma must be at least one.");

  handle->vmg_param->cycle_type = cycle_type;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_get_cycle_type(FCS handle, fcs_int* cycle_type)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  *cycle_type = handle->vmg_param->cycle_type;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_set_precision(FCS handle, fcs_float precision)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  if (precision < 0.0)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "The requested precision must be positive.");

  handle->vmg_param->precision = precision;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_get_precision(FCS handle, fcs_float* precision)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  *precision = handle->vmg_param->precision;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_set_near_field_cells(FCS handle, fcs_int near_field_cells)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  if (near_field_cells < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "Number of near field cells  must be positive.");

  handle->vmg_param->near_field_cells = near_field_cells;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_get_near_field_cells(FCS handle, fcs_int* near_field_cells)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  *near_field_cells = handle->vmg_param->near_field_cells;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_set_interpolation_order(FCS handle, fcs_int interpolation_order)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  if (interpolation_order < 3)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "Interpolation order must be greater or equal three.");

  handle->vmg_param->interpolation_order = interpolation_order;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_get_interpolation_order(FCS handle, fcs_int* interpolation_order)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  *interpolation_order = handle->vmg_param->interpolation_order;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_set_discretization_order(FCS handle, fcs_int discretization_order)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  if (discretization_order != 2 && discretization_order != 4)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "discretization_order must be 2 or 4.");

  handle->vmg_param->discretization_order = discretization_order;

  return FCS_RESULT_SUCCESS;
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
FCSResult fcs_vmg_get_discretization_order(FCS handle, fcs_int* discretization_order)
{
  VMG_CHECK_RETURN_RESULT(handle, __func__);

  *discretization_order = handle->vmg_param->discretization_order;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_vmg_check(FCS handle)
{
  FCSResult result;

  VMG_CHECK_RETURN_RESULT(handle, __func__);

  fcs_int max_level;
  result = fcs_vmg_get_max_level(handle, &max_level);
  CHECK_RESULT_RETURN(result);

  if (max_level == -1)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "vmg finest level not set.");
  if (max_level < 3)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg finest level must be greater than 2.");

  fcs_int max_iterations;
  result = fcs_vmg_get_max_iterations(handle, &max_iterations);
  CHECK_RESULT_RETURN(result);

  if (max_iterations == -1)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "vmg maximum number of iterations not set.");
  if (max_iterations < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg maximum number of iterations must be positive.");

  fcs_int smoothing_steps;
  result = fcs_vmg_get_smoothing_steps(handle, &smoothing_steps);
  CHECK_RESULT_RETURN(result);

  if (smoothing_steps == -1)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "vmg number of smoothing steps not set.");
  if (smoothing_steps < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg number of smoothing steps must be positive.");

  fcs_int cycle_type;
  result = fcs_vmg_get_cycle_type(handle, &cycle_type);
  CHECK_RESULT_RETURN(result);

  if (cycle_type == -1)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "vmg cycle number not set.");
  if (cycle_type < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg cycle number must be positive.");

  fcs_float precision;
  result = fcs_vmg_get_precision(handle, &precision);
  CHECK_RESULT_RETURN(result);

  if (precision == -1.)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "vmg desired precision not set.");

  fcs_int near_field_cells;
  result = fcs_vmg_get_near_field_cells(handle, &near_field_cells);
  CHECK_RESULT_RETURN(result);

  if (near_field_cells == -1)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "vmg number of near field cells not set.");
  if (near_field_cells < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg number of near field cells must be positive.");

  fcs_int interpolation_order;
  result = fcs_vmg_get_interpolation_order(handle, &interpolation_order);
  CHECK_RESULT_RETURN(result);

  if (interpolation_order == -1)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "vmg interpolation order not set.");
  if (interpolation_order < 3)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg interpolation order must be greater or equal three.");

  fcs_int discretization_order;
  result = fcs_vmg_get_discretization_order(handle, &discretization_order);
  CHECK_RESULT_RETURN(result);

  if (discretization_order == -1)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "vmg discretization order not set.");
  if (discretization_order != 2 && discretization_order != 4)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg discretization order must be 2 or 4.");

  const fcs_float* box_a = fcs_get_box_a(handle);
  const fcs_float* box_b = fcs_get_box_b(handle);
  const fcs_float* box_c = fcs_get_box_c(handle);

  if (!fcs_is_cubic(box_a, box_b, box_c))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg requires the box to be cubic.");

  if (!fcs_uses_principal_axes(box_a, box_b, box_c))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg requires the box vectors to be parallel to the principal axes.");

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_vmg_library_check(FCS handle)
{
  int rval = VMG_fcs_check();

  switch (rval)
    {
    case 0:
      return FCS_RESULT_SUCCESS;
      break;
    case 1:
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg requires the number of grid points on each process in every direction to be greater than the number of near field cells.");
      break;
    default:
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "vmg unknown error code");
      break;
    }
}

FCSResult fcs_vmg_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  char *param = *current;
  char *cur = *next;

  *matched = 0;

  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_max_level",            vmg_set_max_level,            FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_max_iterations",       vmg_set_max_iterations,       FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_smoothing_steps",      vmg_set_smoothing_steps,      FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_cycle_type",           vmg_set_cycle_type,           FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_precision",            vmg_set_precision,            FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_near_field_cells",     vmg_set_near_field_cells,     FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_interpolation_order",  vmg_set_interpolation_order,  FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_discretization_order", vmg_set_discretization_order, FCS_PARSE_VAL(fcs_int));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_vmg_print_parameters(FCS handle)
{
  fcs_int level;
  fcs_int max_iter;
  fcs_int smoothing_steps;
  fcs_int cycle_type;
  fcs_float precision;
  fcs_int near_field_cells;
  fcs_int interpolation_order;
  fcs_int discretization_order;

  VMG_CHECK_RETURN_RESULT(handle, __func__);

  fcs_vmg_get_max_level(handle, &level);
  fcs_vmg_get_max_iterations(handle, &max_iter);
  fcs_vmg_get_smoothing_steps(handle, &smoothing_steps);
  fcs_vmg_get_cycle_type(handle, &cycle_type);
  fcs_vmg_get_precision(handle, &precision);
  fcs_vmg_get_near_field_cells(handle, &near_field_cells);
  fcs_vmg_get_interpolation_order(handle, &interpolation_order);
  fcs_vmg_get_discretization_order(handle, &discretization_order);

  printf("vmg max level:            %" FCS_LMOD_INT "d\n", level);
  printf("vmg max iterations:       %" FCS_LMOD_INT "d\n", max_iter);
  printf("vmg smoothing steps:      %" FCS_LMOD_INT "d\n", smoothing_steps);
  printf("vmg cycle_type:           %" FCS_LMOD_INT "d\n", cycle_type);
  printf("vmg precision:            %e\n", precision);
  printf("vmg near field cells:     %" FCS_LMOD_INT "d\n", near_field_cells);
  printf("vmg interpolation degree: %" FCS_LMOD_INT "d\n", interpolation_order);
  printf("vmg discretization order: %" FCS_LMOD_INT "d\n", discretization_order);
  
  return FCS_RESULT_SUCCESS;
}
