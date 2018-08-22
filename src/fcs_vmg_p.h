/*
  Copyright (C) 2011-2012 Rene Halver, Frederik Heber, Julian Iseringhausen

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



#ifndef FCS_VMG_P_INCLUDED
#define FCS_VMG_P_INCLUDED

#include "fcs_definitions.h"
#include "fcs_result_p.h"
#include "fcs_interface_p.h"

/**
 * @file fcs_vmg_p.h
 * @brief file containing all vmg specific functions
 * @author Rene Halver
 */

typedef struct fcs_vmg_parameters_t *fcs_vmg_parameters;

/**
 * @brief Function to set all vmg parameters with a single call
 *
 * @param handle FCS-object that contains the parameters
 * @param max_level The maximum level of the algorithm, i.e. n_gridpoints = 2^max_level.
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
FCSResult fcs_vmg_setup(FCS handle, fcs_int max_level, fcs_int max_iterations,
			       fcs_int smoothing_steps, fcs_int cycle_type,
			       fcs_float precision, fcs_int near_field_cells,
			       fcs_int interpolation_order,
			       fcs_int discretization_order);

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
FCSResult fcs_vmg_set_max_level(FCS handle, fcs_int max_level);

/**
 * @brief Get the maximum level of the multigrid algorithm.
 * The number of grid points of the finest multigrid level
 * will be 2^max_level.
 *
 * @param handle FCS-object that contains the parameter
 * @param max_level Number of iterations on return
 *
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_get_max_level(FCS handle, fcs_int* max_level);

/**
 * @brief Set the maximum number of multigrid iterations.
 *
 * @param handle FCS-object that contains the parameter
 * @param max_iterations Number of iterations.
 *
 * @return FCSResult-object containing the return state.
 */
FCSResult fcs_vmg_set_max_iterations(FCS handle, fcs_int max_iterations);

/**
 * @brief Get the maximum number of multigrid iterations.
 *
 * @param handle FCS-object that contains the parameter.
 * @param max_iterations Number of iterations.
 *
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_get_max_iterations(FCS handle, fcs_int* max_iterations);

/**
 * @brief Set the number of pre/postsmoothing steps on each level
 *
 * @param handle FCS-object that contains the parameter
 * @param smoothing_steps Number of smoothing steps
 *
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_set_smoothing_steps(FCS handle, fcs_int smoothing_steps);

/**
 * @brief Get the number of pre/postsmoothing steps on each level
 *
 * @param handle FCS-object that contains the parameter
 * @param smoothing_steps Number of pre/postsmoothing steps
 *
 * @return FCSReturn-object containing the return state.
 */
FCSResult fcs_vmg_get_smoothing_steps(FCS handle, fcs_int *smoothing_steps);

/**
 * @brief Set the cycle_type-number of the multigrid cycle used. E.g. 1 corresponds
 *        to a V-Cycle and 2 corresponds to a W-Cycle.
 *
 * @param handle FCS-object that contains the parameter.
 * @param cycle_type Gamma.
 *
 * @return FCSReturn-object containing the return state.
 */
FCSResult fcs_vmg_set_cycle_type(FCS handle, fcs_int cycle_type);

/**
 * @brief Get the cycle_type-number of the multigrid cycle used. E.g. 1 corresponds
 *        to a V-Cycle and 2 corresponds to a W-Cycle.
 *
 * @param handle FCS-object that contains the parameter.
 * @param cycle_type type of multigrid cycle.
 *
 * @return FCSResult-object containing the return value.
 */
FCSResult fcs_vmg_get_cycle_type(FCS handle, fcs_int *cycle_type);

/**
 * @brief Set the precision of the vmg algorithm. Currently, the relative residual
 *        computed in the discrete L2 norm will be tested against this value.
 *
 * @param handle FCS-object that contains the parameter.
 * @param precision Precision.
 *
 * @return FCSResult-object containing the return value.
 */
FCSResult fcs_vmg_set_precision(FCS handle, fcs_float precision);

/**
 * @brief Get the precision of the vmg algorithm. Currently, the relative residual
 *        computed in the discrete L2 norm will be tested against this value.
 *
 * @param handle FCS-object that contains the parameter.
 * @param precision Precision of the vmg algorithm.
 *
 * @return FCSResult-object containing the return value.
 */
FCSResult fcs_vmg_get_precision(FCS handle, fcs_float *precision);

/**
 * @brief Set the number of near field cells for separating the
 *        near/far field part of the potential.
 *
 * @param handle FCS-object that contains the parameter.
 * @param near_field_cells Near field cells.
 *
 * @return FCSResult-object containing the return value.
 */
FCSResult fcs_vmg_set_near_field_cells(FCS handle, fcs_int near_field_cells);

/**
 * @brief Get the number of near field cells for separating the
 *        near/far field part of the potential.
 *
 * @param handle FCS-object that contains the parameter.
 * @param near_field_cells Number of near field cells.
 *
 * @return FCSResult-object containing the return value.
 */
FCSResult fcs_vmg_get_near_field_cells(FCS handle, fcs_int *near_field_cells);

/**
 * @brief Set the interpolation order for interpolating
 *        the gridded potential to the particle positions.
 *
 * @param handle FCS-object that contains the parameter.
 * @param interpolation_order Interpolation order.
 *
 * @return FCSResult-object containing the return value.
 */
FCSResult fcs_vmg_set_interpolation_order(FCS handle, fcs_int interpolation_order);

/**
 * @brief Get the interpolation order for interpolating
 *        the gridded potential to the particle positions.
 *
 * @param handle FCS-object that contains the parameter.
 * @param interpolation_order Interpolation order.
 *
 * @return FCSResult-object containing the return value.
 */
FCSResult fcs_vmg_get_interpolation_order(FCS handle, fcs_int *interpolation_order);

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
FCSResult fcs_vmg_set_discretization_order(FCS handle, fcs_int discretization_order);

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
FCSResult fcs_vmg_get_discretization_order(FCS handle, fcs_int *discretization_order);

/**
 * @brief Print runtimes of various vmg subsystems. vmg has to be configured
 *        with --enable-debug-measure-time in order to do so.
 *
 */
void vmg_fcs_print_timer();

#endif
