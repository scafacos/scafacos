/*
  Copyright (C) 2011-2013 Michael Pippig
  Copyright (C) 2011-2012 Rene Halver
  Copyright (C) 2011 Sebastian Banert

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



#ifndef FCS_P2NFFT_INCLUDED
#define FCS_P2NFFT_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "fcs_p2nfft_p.h"
#include "FCSResult.h"
#include "FCSInterface.h"




typedef struct fcs_p2nfft_parameters_t{
  fcs_int dummy;
} fcs_p2nfft_parameters_t;


fcs_float fcs_p2nfft_compute_near_potential(
    FCS handle, fcs_float dist);
fcs_float fcs_p2nfft_compute_near_field(
   FCS handle, fcs_float dist);
void fcs_p2nfft_compute_near(
    FCS handle, fcs_float dist,
    fcs_float *potential, fcs_float *field);


/**
 * @brief function to check if obligatory p2nfft information is set
 * @param handle FCS-object to be checked if obligatory
 * information needed by p2nfft is set (no plausibility check!)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_p2nfft_check(FCS handle);

/**
 * @brief initialization routine for the basic parameters needed by p2nfft
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_p2nfft_init(FCS handle);

/**
 * @brief tuning method for setting/calculating last parameters, for which positions,
 *  charges, etc. are needed for p2nfft


 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param local_max_particles size of allocated arrays
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_p2nfft_tune(
    FCS handle, fcs_int local_particles, fcs_int local_max_particles,
    fcs_float *positions,  fcs_float *charges);

/**
 * @brief run method for p2nfft
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param local_max_particles size of allocated arrays
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)


 * @param charges fcs_float* list of charges
 * results after the run
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_p2nfft_run(
    FCS handle, fcs_int local_particles, fcs_int local_max_particles,
    fcs_float *positions,  fcs_float *charges,
    fcs_float *field, fcs_float *potentials);

/**
 * @brief a method to calculate p2nfft near-field potentials with the chosen method (without charges of particles)
 *        e.g.: instead of q1q2/r the function return 1/r
 * @param handle pointer to a FCS on which the FCS with the information
 * are saved
 * @param local_particles the amount of particles on local process
 * @param r distance for which the potential
 * @return fcs_float containing the potential value for the given distance (without charges, s. above)
 */
extern FCSResult fcs_p2nfft_near_field_potential (
    FCS handle, fcs_float);

/**
 * @brief clean-up method for p2nfft
 * @param handle the FCS-object, which contains the parameters
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_p2nfft_destroy(FCS handle);

/**
 * @brief function to activate computation of the virial 
 * @param handle FCS-object that contains the parameter
 * @param flag whether or not to compute the virial in the next call
 * to fcs_run().
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_p2nfft_require_virial(FCS handle, fcs_int flag);

/**
 * @brief function to fetch the virial
 * @param handle FCS-object that contains the parameter
 * @param virial pointer to the array where the virial is returned.
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_p2nfft_get_virial(FCS handle, fcs_float *virial);


FCSResult fcs_p2nfft_set_max_particle_move(FCS handle, fcs_float max_particle_move);
FCSResult fcs_p2nfft_set_resort(FCS handle, fcs_int resort);
FCSResult fcs_p2nfft_get_resort(FCS handle, fcs_int *resort);
FCSResult fcs_p2nfft_get_resort_availability(FCS handle, fcs_int *availability);
FCSResult fcs_p2nfft_get_resort_particles(FCS handle, fcs_int *resort_particles);
FCSResult fcs_p2nfft_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm);
FCSResult fcs_p2nfft_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm);
FCSResult fcs_p2nfft_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm);


#endif
