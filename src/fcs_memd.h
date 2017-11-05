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
 * @file   fcs_memd.h
 * @author Florian Fahrenberger <florian.fahrenberger@icp.uni-stuttgart.de>
 * @date   Tue May 10 14:47:00 2011
 *
 * @brief  Private declarations for the interface between FCS and memd libraries.
 *
 */

#ifndef FCS_MEMD_INCLUDED
#define FCS_MEMD_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "fcs_memd_p.h"
#include "FCSResult.h"
#include "FCSInterface.h"

typedef struct fcs_memd_parameters_t
{
    fcs_int dummy;
}fcs_memd_parameters_t;


/**
 * @brief initialization routine for the basic parameters needed by memd
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_memd_init(FCS handle);

/**
 * @brief tuning method for setting/calculating last parameters, for which positions,
 *  charges, etc. are needed for memd
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_memd_tune(FCS handle, fcs_int local_particles, fcs_float *positions, fcs_float *charges);

/**
 * @brief run method for memd
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @param output FCSOutput* pointer that contains a FCSOutput-object with the
 * results after the run
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_memd_run(FCS handle, fcs_int local_particles, fcs_float *positions,  fcs_float *charges, fcs_float *fields, fcs_float *potentials);

/**
 * @brief clean-up method for memd
 * @param handle the FCS-object, which contains the parameters
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_memd_destroy(FCS handle);

FCSResult fcs_memd_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched);
FCSResult fcs_memd_print_parameters(FCS handle);

#endif /* FCS_MEMD_INCLUDED */
