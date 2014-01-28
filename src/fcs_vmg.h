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



/**
 * @file   fcs_vmg.h
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Sat Apr 16 19:14:29 2011
 *
 * @brief  Private declarations for the interface between FCS and vmg libraries.
 *
 */

#ifndef FCS_VMG_INCLUDED
#define FCS_VMG_INCLUDED

#include "fcs_vmg_p.h"
#include "FCSResult.h"
#include "FCSInterface.h"

typedef struct fcs_vmg_parameters_t
{
  fcs_int max_level;
  fcs_int max_iterations;
  fcs_int smoothing_steps;
  fcs_int cycle_type;
  fcs_float precision;
  fcs_int near_field_cells;
  fcs_int interpolation_order;
  fcs_int discretization_order;
}fcs_vmg_parameters_t;

/**
 * @brief function to set default values if no values are set
 * @param handle the FCS-obect into which the default parameters will be entered
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_set_default(FCS handle);

/**
 * @brief function to check if obligatory vmg information is set
 * @param handle FCS-object to be checked if obligatory
 * information needed by vmg is set.
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_check(FCS handle);

/**
 * @brief function to check if the vmg library reports any
 *        configuration errors
 * @param handle FCS-object to be checked if obligatory
 * information needed by vmg is set.
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_library_check(FCS handle);

/**
 * @brief initialization routine for the basic parameters needed by vmg
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_init(FCS handle);

/**
 * @brief tuning method for setting/calculating last parameters, for which positions,
 *  charges, etc. are needed for vmg
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param local_max_particles size of allocated arrays
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles,
			      fcs_float *positions, fcs_float *charges);

/**
 * @brief run method for vmg
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param local_max_particles size of allocated arrays
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles,
			     fcs_float *positions, fcs_float *charges,
			     fcs_float *field, fcs_float *potentials);

/**
 * @brief a method to calculate vmg near-field potentials with the chosen method (without charges of particles)
 *        e.g.: instead of q1q2/r the function return 1/r
 * @param handle pointer to a FCS on which the FCS with the information
 * are saved
 * @param local_particles the amount of particles on local process
 * @param r distance for which the potential
 * @return fcs_float containing the potential value for the given distance (without charges, s. above)
 */
FCSResult fcs_vmg_near_field_potential (FCS handle, fcs_float);

/**
 * @brief clean-up method for vmg
 * @param handle the FCS-object, which contains the parameters
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_destroy(FCS handle);

/**
 * @brief function to activate computation of the virial
 * @param handle FCS-object that contains the parameter
 * @param flag whether or not to compute the virial in the next call
 * to fcs_run().
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_require_virial(FCS handle, fcs_int flag);

/**
 * @brief function to fetch the virial
 * @param handle FCS-object that contains the parameter
 * @param virial pointer to the array where the virial is returned.
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_vmg_get_virial(FCS handle, fcs_float *virial);

FCSResult fcs_vmg_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched);
FCSResult fcs_vmg_print_parameters(FCS handle);

/**
 * @brief External interface definition for setup of the vmg library.
 *
 * @param max_level The finest level (n_gridpoints = 2^max_level).
 * @param periodic Periodicity
 * @param max_iteration The maximum number of multigrid iterations.
 * @param smoothing_steps The number of smoothing steps. Must be positive.
 * @param cycle_type The cycle number (e.g. 1 for V-cycle, 2 for W-cycle)
 * @param precision The desired precision. This number will be tested against the relative residual in the discrete L2-norm.
 * @param box_offset Offset of the box.
 * @param box_size Size of the box.
 * @param near_field_cells Splitting of short/long range part of the potential.
 * @param interpolation_order Interpolation order.
 * @param discretization_order Discretization order.
 * @param comm MPI communicator.
 */
void VMG_fcs_setup(fcs_int max_level, const fcs_int* periodic, fcs_int max_iteration,
			  fcs_int smoothing_steps, fcs_int cycle_type, fcs_float precision,
			  const fcs_float* box_offset, fcs_float box_size, fcs_int near_field_cells,
			  fcs_int interpolation_order, fcs_int discretization_order,
			  MPI_Comm comm);

/**
 * @brief External interface definition for running internal vmg library checks.
 *
 * @return Error code.
 */
int VMG_fcs_check();

/**
 * @brief Run the vmg solver
 *
 * @param x Charge positions
 * @param q Charges
 * @param p Potentials
 * @param f Forces
 * @param num_particles_local Number of particles on this process.
 */
void VMG_fcs_run(fcs_float* x, fcs_float* q, fcs_float* p, fcs_float* f, fcs_int num_particles_local);

 /**
 * @brief Bring the vmg library back to the starting state.
 *
 */
void VMG_fcs_destroy();

#endif /* FCS_VMG_INCLUDED */
