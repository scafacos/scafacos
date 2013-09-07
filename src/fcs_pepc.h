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



#ifndef FCS_PEPC_INCLUDED
#define FCS_PEPC_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "fcs_pepc_p.h"
#include "FCSResult.h"
#include "FCSInterface.h"


/**
 * @file fcs_pepc.h
 * @brief file containing all pepc specific functions
 * @author Lukas Arnold
 */



typedef struct fcs_pepc_parameters_t
{
  /* list of parameters used by pepc */

  /* renormalization epsilon */
  fcs_float epsilon;
  /* MAC theta */
  fcs_float theta;
  /* flag if virial calculation is required */
  fcs_int requirevirial;
  /* number of walk threads per MPI rank */
  fcs_int num_walk_threads;
  /* switch for activating load balancing. may only be set >0 if the frontend does not reorder the particles */
  fcs_int load_balancing;
  /* switch for activating dipole correction for periodic systems */
  fcs_int dipole_correction;
  /* internal memory usage parameter */
  fcs_float npm;
  /* pepc_debug level */
  fcs_int debug_level;

} fcs_pepc_parameters_t;

typedef struct fcs_pepc_internal_t
{
  /* internal data structure for pepc */

  /* virial matrix elements */
  fcs_float virial[9];

  /* pointer to work, used for load balancing */
  fcs_float* work;
  fcs_int work_length;

} fcs_pepc_internal_t;


/**
 * @brief function to check if obligatory pepc information is set
 * @param handle FCS-object to be checked if obligatory
 * information needed by pepc is set (no plausibility check!)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_check(FCS handle);

/**
 * @brief initialization routine for the basic parameters needed by pepc
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_init(FCS handle);

/**
 * @brief tuning routine for the basic parameters needed by pepc
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions,  fcs_float *charges);

/**
 * @brief run method for pepc
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param local_max_particles size of allocated arrays
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)


 * @param charges fcs_float* list of charges
 * @param output FCSOutput* pointer that contains a FCSOutput-object with the
 * results after the run
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles, 
			      fcs_float *positions, 
			      fcs_float *charges, fcs_float *field, fcs_float *potentials);

/**
 * @brief a method to calculate pepc near-field potentials with the chosen method (without charges of particles)
 *        e.g.: instead of q1q2/r the function return 1/r
 * @param handle pointer to a FCS on which the FCS with the information
 * are saved
 * @param local_particles the amount of particles on local process
 * @param r distance for which the potential
 * @return fcs_float containing the potential value for the given distance (without charges, s. above)
 */
extern FCSResult fcs_pepc_near_field_potential (FCS handle, fcs_float);

/**
 * @brief clean-up method for pepc
 * @param handle the FCS-object, which contains the parameters
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_pepc_destroy(FCS handle);

/**
 * @brief function activate virial calculation in pepc
 * @param handle FCS-object that contains the parameter
 * @param choice 1 - user requires virial calculation and can retrieve the result using fcs_pepc_get_virial() later, other values: deactivate this feture
 * @return the debug level
 */
extern FCSResult fcs_pepc_require_virial(FCS handle, fcs_int choice);

/**
 * @brief function to get the virial matrix after force caluclation (i.e. after having set fcs_pepc_require_virial(handle, 1) and calling fcs_pepc_run()
 * @param handle FCS-object that contains the parameter
 * @param virial field for storing virial matrix elements
 * @return the debug level
 */
extern FCSResult fcs_pepc_get_virial(FCS handle, fcs_float virial[9]);

/**
 * @brief pepc init fortran routine that is wrapped by scafacos
 */
extern void pepc_scafacos_initialize(MPI_Fint *comm);

/**
 * @brief pepc finalize fortran routine that is wrapped by scafacos
 */
extern void pepc_scafacos_finalize(MPI_Fint *comm);

/**
 * @brief central pepc fortran routine that is wrapped by scafacos
 */
extern void pepc_scafacos_run(fcs_int *local_particles, fcs_int *total_particles,
			      fcs_float *positions, fcs_float *charges, 
			      fcs_float *efield, fcs_float *potentials, fcs_float *work,
			      fcs_float *virial,
			      fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_int *periodicity, 
			      fcs_int *lattice_corr, fcs_float *eps, fcs_float *theta, 
                              fcs_int *db_level, fcs_int *num_walk_threads, fcs_float *npm );

#endif
