/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann

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



#ifndef FCS_WOLF_P_INCLUDED
#define FCS_WOLF_P_INCLUDED


#include "fcs_definitions.h"
#include "fcs_result_p.h"
#include "fcs_interface_p.h"


/**
 * @file fcs_wolf_p.h
 * @brief file containing the method specific interface functions
 * of the wolf solver (public version)
 * @author Michael Hofmann
 */


/**
 * @brief data structure with parameters of the wolf solver
 */
typedef struct fcs_wolf_parameters_t *fcs_wolf_parameters;


/**
 * @brief function to set the cutoff radius
 * @param handle FCS-object
 * @param cutoff cutoff radius
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_wolf_set_cutoff(FCS handle, fcs_float cutoff);


/**
 * @brief function to get the current cutoff radius
 * @param handle FCS-object
 * @param cutoff current cutoff radius
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_wolf_get_cutoff(FCS handle, fcs_float *cutoff);


/**
 * @brief function to set the ewald splitting parameter alpha
 * @param handle FCS-object
 * @param alpha ewald splitting parameter alpha
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_wolf_set_alpha(FCS handle, fcs_float alpha);


/**
 * @brief function to get the current ewald splitting parameter alpha
 * @param handle FCS-object
 * @param alpha current ewald splitting parameter alpha
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_wolf_get_alpha(FCS handle, fcs_float *alpha);


/**
 * @brief function to set all solver parameters
 * @param handle FCS-object
 * @param cutoff cutoff radius (see ::fcs_wolf_set_cutoff)
 * @param alpha ewald splitting parameter alpha (see ::fcs_wolf_set_alpha)
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_wolf_setup(FCS handle, fcs_float cutoff, fcs_float alpha);


/**
 * @brief function to set all solver parameters (Fortran wrapper)
 * @param handle FCS-object
 * @param cutoff cutoff radius (see ::fcs_wolf_set_cutoff)
 * @param alpha ewald splitting parameter alpha (see ::fcs_wolf_set_alpha)
 * @param return_value return state
 */
void fcs_wolf_setup_f(void *handle, fcs_float cutoff, fcs_float alpha, fcs_int *return_value);


#endif
