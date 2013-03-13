/*
  Copyright (C) 2011-2012 Rene Halver, Lukas Arnold, Mathias Winkel

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



#ifndef FCS_PEPC_P_INCLUDED
#define FCS_PEPC_P_INCLUDED

#include "FCSDefinitions.h"
#include "FCSResult_p.h"
#include "FCSInterface_p.h"

/**
 * @file fcs_pepc_p.h
 * @brief file containing all pepc specific functions (public version)
 * @author Rene Halver
 */

typedef struct fcs_pepc_parameters_t *fcs_pepc_parameters;

/**
 * @brief function to set the optional pepc epsilon parameter
 * @param handle FCS-object that is modified
 * @param epsilon epsilon (Plummer potential parameter)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_set_epsilon(FCS handle, fcs_float epsilon);

/**
 * @brief function to get the optional pepc epsilon parameter
 * @param handle FCS-object that contains the parameter
 * @param epsilon epsilon (Plummer potential parameter)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_get_epsilon(FCS handle, fcs_float* epsilon);

/**
 * @brief function to set the optional pepc theta parameter
 * @param handle FCS-object that is modified
 * @param theta multipole acceptance parameter for Barnes-Hut MAC
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_set_theta(FCS handle, fcs_float theta);

/**
 * @brief function to get the optional pepc theta parameter
 * @param handle FCS-object that contains the parameter
 * @param theta multipole acceptance parameter for Barnes-Hut MAC
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_get_theta(FCS handle, fcs_float* theta);

/**
 * @brief function to set the pepc debug level
 * @param handle FCS-object that is modified
 * @param level user-set debug level of pepc solver
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_set_debug_level(FCS handle, fcs_int level);

/**
 * @brief function to get the pepc debug level
 * @param handle FCS-object that contains the parameter
 * @param level user-set debug level of pepc solver
 * @return the debug level
 */
extern FCSResult fcs_pepc_get_debug_level(FCS handle, fcs_int* level);

/**
 * @brief function to set pepcs number of walk threads per MPI rank
 * @param handle FCS-object that contains the parameter
 * @param num_walk_threads number of traversal threads per rank
 */
extern FCSResult fcs_pepc_set_num_walk_threads(FCS handle, fcs_int num_walk_threads);

/**
 * @brief function to get pepcs number of walk threads per MPI rank
 * @param handle FCS-object that contains the parameter
 * @return num_walk_threads number of traversal threads per rank
 */
extern FCSResult fcs_pepc_get_num_walk_threads(FCS handle, fcs_int* num_walk_threads);

/**
 * @brief function for setting pepcs switch for load balancing
 * @param handle FCS-object that contains the parameter
 * @param load_balancing if >0 load balancing is enabled. may only be activated if the frontend does not reorder any particles
 */
extern FCSResult fcs_pepc_set_load_balancing(FCS handle, fcs_int load_balancing);

/**
 * @brief function for getting pepcs switch for load balancing
 * @param handle FCS-object that contains the parameter
 * @param load_balancing if >0 load balancing is enabled. may only be activated if the frontend does not reorder any particles
 */
extern FCSResult fcs_pepc_get_load_balancing(FCS handle, fcs_int* load_balancing);

/**
 * @brief function for setting pepcs switch for adding a dipole correction for periodic systems
 * @param handle FCS-object that contains the parameter
 * @param dipole_correction if >0 the dipole correction is enabled
 */
extern FCSResult fcs_pepc_set_dipole_correction(FCS handle, fcs_int dipole_correction);

/**
 * @brief function for getting pepcs switch for adding a dipole correction for periodic systems
 * @param handle FCS-object that contains the parameter
 * @param dipole_correction if >0 the dipole correction is enabled
 */
extern FCSResult fcs_pepc_get_dipole_correction(FCS handle, fcs_int* dipole_correction);

/**
 * @brief function to set the optional pepc npm (aka np_mult)
 * @param handle FCS-object that is modified
 * @param npm npm (aka np_mult)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_set_npm(FCS handle, fcs_float npm);

/**
 * @brief function to get the optional pepc npm (aka np_mult)
 * @param handle FCS-object that is modified
 * @param npm npm (aka np_mult)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_get_npm(FCS handle, fcs_float* npm);


extern FCSResult fcs_pepc_setup(FCS handle, fcs_float epsilon, fcs_float theta);

void fcs_pepc_setup_f(void *handle, fcs_float epsilon, fcs_float theta, fcs_int level, fcs_int *return_value);

#endif
