/*
  Copyright (C) 2011-2012 Michael Hofmann, Rene Halver

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


#ifndef FCS_DIRECT_INCLUDED
#define FCS_DIRECT_INCLUDED

#include <mpi.h>

#include "fcs_direct_p.h"
#include "FCSResult.h"
#include "FCSInterface.h"
#include "direct/directc.h"


/**
 * @file fcs_direct.h
 * @brief file containing the method specific interface functions
 * for the direct solver (private version)
 * @author Michael Hofmann, Rene Halver
 */


typedef struct fcs_direct_parameters_t
{
  fcs_directc_t directc;

  fcs_int metallic_boundary_conditions;

} fcs_direct_parameters_t;


typedef struct _fcs_direct_context_t
{
  MPI_Comm comm;
  int comm_size, comm_rank;

  fcs_float last_runtime;

} fcs_direct_context_t;


#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief initialization routine for the basic parameters needed by the direct solver
 * @param handle the FCS-object into which the method specific parameters can be entered
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_init(FCS handle);


/**
 * @brief deallocation routine for the custom parameters needed by the direct solver
 * @param handle the FCS-object for which the custom parameter array should be deallocated
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_destroy(FCS handle);


FCSResult fcs_direct_check(FCS handle);


/**
 * @brief tuning method for setting/calculating last parameters, for which positions,
 *  charges, etc. are needed for the direct solver
 * @param handle the FCS-object into which the method specific parameters can be entered
 * @param positons fcs_float* list of positions of particles in form (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions, fcs_float *charges);


/**
 * @brief run method for direct
 * @param handle the FCS-object into which the method specific parameters can be entered
 * @param positons fcs_float* list of positions of particles in form (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @param output FCSOutput* pointer that contains a FCSOutput-object with the results after the run
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles,
                         fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials);


/**
 * @brief a method to calculate direct near-field potentials with the chosen method (without charges of particles)
 *        e.g.: instead of q1q2/r the function return 1/r
 * @param handle pointer to a FCS on which the FCS with the information are saved
 * @param dist fcs_float distance for which the potential
 * @return fcs_float containing the potential value for the given distance (without charges, s. above)
 */
FCSResult fcs_direct_near_field_potential(FCS handle, fcs_float dist);


/**
 * @brief function to enable the virial computation
 * @param handle FCS-object
 * @param compute_virial whether virial should be computed or not
 *        = 0 - no virial
 *        != 0 - compute virial
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_require_virial(FCS handle, fcs_int compute_virial);


/**
 * @brief function to return the virial
 * @param handle FCS-object
 * @param virial array to store the virial
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_direct_get_virial(FCS handle, fcs_float *virial);


FCSResult fcs_direct_set_max_particle_move(FCS handle, fcs_float max_particle_move);


#ifdef __cplusplus
}
#endif


#endif
