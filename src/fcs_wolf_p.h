/*
  Copyright (C) 2011,2012,2013 Michael Hofmann

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


#include "FCSDefinitions.h"
#include "FCSResult_p.h"
#include "FCSInterface_p.h"


/**
 * @file fcs_wolf_p.h
 * @brief file containing the method specific interface functions
 * for the wolf solver (public version)
 * @author Michael Hofmann, Rene Halver
 */
typedef struct fcs_wolf_parameters_t *fcs_wolf_parameters;


/**
 * @brief function to set the cutoff radius
 * @param handle FCS-object
 * @param fcs_float cutoff radius
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_wolf_set_cutoff(FCS handle, fcs_float cutoff);


/**
 * @brief function to get the current cutoff radius
 * @param handle FCS-object
 * @param *fcs_float current cutoff radius
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_wolf_get_cutoff(FCS handle, fcs_float *cutoff);


/**
 * @brief function to set the ewald splitting parameter alpha
 * @param handle FCS-object
 * @param fcs_float parameter alpha
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_wolf_set_alpha(FCS handle, fcs_float alpha);


/**
 * @brief function to get the current ewald splitting parameter alpha
 * @param handle FCS-object
 * @param *fcs_float alpha current parameter alpha
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_wolf_get_alpha(FCS handle, fcs_float *alpha);


/**
 * @brief combined setter for all wolf solver parameters
 * @param handle FCS-object
 * @param fcs_float cutoff radius (see fcs_wolf_set_cutoff)
 * @param fcs_float alpha radius (see fcs_wolf_set_alpha)
 * @return FCSResult-object containing th return state
 **/
FCSResult fcs_wolf_setup(FCS handle, fcs_float cutoff, fcs_float alpha);


/**
 * @brief combined setter for all wolf solver parameters (FORTRAN WRAPPER)
 * @param handle void* to be casted into FCS so that the data is stored into it
 * @param cutoff fcs_float the chosen cutoff (see fcs_wolf_set_cutoff)
 * @param return_value fcs_int* return value for FORTRAN interface
 **/
void fcs_wolf_setup_f(void *handle, fcs_float cutoff, fcs_float alpha, fcs_int *return_value);


#endif
