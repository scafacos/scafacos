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



#ifndef _FCS_P3M_H
#define _FCS_P3M_H

#include <mpi.h>
#include "fcs_p3m_p.h"
#include "FCSResult.h"
#include "FCSInterface.h"

/**
 * @brief initialization routine for the basic parameters needed by p3m
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_p3m_init(FCS handle);

/**
 * @brief tuning method for setting/calculating last parameters, for which positions,
 *  charges, etc. are needed for p3m
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param local_max_particles size of allocated arrays
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @return FCSResult-object containing the return state


 */
FCSResult fcs_p3m_tune(FCS handle,
		       fcs_int local_particles, fcs_int local_max_particles,
		       fcs_float *positions,  fcs_float *charges);

/**
 * @brief run method for p3m
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
FCSResult fcs_p3m_run(FCS handle,
		      fcs_int local_particles, fcs_int local_max_particles,
		      fcs_float *positions,  fcs_float *charges,
		      fcs_float *field, fcs_float *potentials);

/**
 * @brief clean-up method for p3m
 * @param handle the FCS-object, which contains the parameters
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_p3m_destroy(FCS handle);

/**
 * @brief function to activate computation of the virial 
 * @param handle FCS-object that contains the parameter
 * @param flag whether or not to compute the virial in the next call
 * to fcs_run().
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_p3m_require_virial(FCS handle, fcs_int flag);

/**
 * @brief function to fetch the virial
 * @param handle FCS-object that contains the parameter
 * @param virial pointer to the array where the virial is returned.
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_p3m_get_virial(FCS handle, fcs_float *virial);

FCSResult fcs_p3m_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched);
FCSResult fcs_p3m_print_parameters(FCS handle);

#endif
