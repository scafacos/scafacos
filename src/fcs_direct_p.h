/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann, Rene Halver

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



#ifndef FCS_DIRECT_P_INCLUDED
#define FCS_DIRECT_P_INCLUDED


#include "FCSDefinitions.h"
#include "FCSResult_p.h"
#include "FCSInterface_p.h"


/**
 * @file fcs_direct_p.h
 * @brief file containing the method specific interface functions
 * for the direct solver (public version)
 * @author Michael Hofmann, Rene Halver
 */
typedef struct fcs_direct_parameters_t *fcs_direct_parameters;


/**
 * @brief function to set the optional cutoff parameter
 * @param handle FCS-object
 * @param cutoff cutoff radius
 *        =0 - no cutoff radius
 *        >0 - compute interactions only inside of the cutoff
 *        <0 - compute interactions only outside of the cutoff
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_set_cutoff(FCS handle, fcs_float cutoff);


/**
 * @brief function to get the optional cutoff parameter
 * @param handle FCS-object
 * @param cutoff current cutoff radius
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_get_cutoff(FCS handle, fcs_float *cutoff);


/**
 * @brief function to set wether the ScaFaCoS near-field solver module should be used to computations with cutoff range
 * @param handle FCS-object
 * @param cutoff_with_near fcs_bool if true, then near-field solver is used for computations with cutoff range
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_set_cutoff_with_near(FCS handle, fcs_bool cutoff_with_near);


/**
 * @brief function to get wether the ScaFaCoS near-field solver module should be used to computations with cutoff range
 * @param handle FCS-object
 * @param cutoff_with_near fcs_bool whether the near-field solver is used for computations with cutoff range
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_get_cutoff_with_near(FCS handle, fcs_bool *cutoff_with_near);


/**
 * @brief function to set the number of image systems used in each (periodic) dimension
 * @param handle FCS-object
 * @param periodic_images array of integers specifying the number of images
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_set_periodic_images(FCS handle, fcs_int *periodic_images);


/**
 * @brief function to set the number of image systems used in each (periodic) dimension
 * @param handle FCS-object
 * @param periodic_images array of integers to store the number of images
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_get_periodic_images(FCS handle, fcs_int *periodic_images);


/**
 * @brief function to set additional input particles (ie, particles for which no results are computed)
 * @param handle FCS-object
 * @param nin_particles number of input particles
 * @param in_positions positions of input particles
 * @param in_charges charges of input particles
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_set_in_particles(FCS handle, fcs_int nin_particles, fcs_float *in_positions, fcs_float *in_charges);


/**
 * @brief combined setter for all direct solver parameters
 * @param handle FCS-object data is stored into
 * @param cutoff the chosen cutoff (see fcs_direct_set_cutoff)
 * @return FCSResult-object containing the return state
 **/
FCSResult fcs_direct_setup(FCS handle, fcs_float cutoff);


/**
 * @brief combined setter for all direct solver parameters (FORTRAN WRAPPER)
 * @param handle void* to be casted into FCS so that the data is stored into it
 * @param cutoff the chosen cutoff (see fcs_direct_set_cutoff)
 * @param return_value return value for FORTRAN interface
 **/
void fcs_direct_setup_f(void *handle, fcs_float cutoff, fcs_int *return_value);


#endif
