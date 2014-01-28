/*
  Copyright (C) 2011-2012 Pedro Sanchez

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



#ifndef _FCS_MMM2D_H
#define _FCS_MMM2D_H

#include <mpi.h>
#include "fcs_mmm2d_p.h"
#include "FCSResult.h"
#include "FCSInterface.h"

/**
 * @brief initialization routine
 */
FCSResult fcs_mmm2d_init(FCS handle);

/**
 * @brief tuning method
 */
FCSResult fcs_mmm2d_tune(FCS handle,
		       fcs_int local_particles, fcs_int local_max_particles,
		       fcs_float *positions,  fcs_float *charges);

/**
 * @brief run method for mmm2d
 */
FCSResult fcs_mmm2d_run(FCS handle,
			fcs_int local_particles, fcs_int local_max_particles,
			fcs_float *positions,  fcs_float *charges,
			fcs_float *field, fcs_float *potentials);

/**
 * @brief clean-up method for mmm2d
 */
FCSResult fcs_mmm2d_destroy(FCS handle);

FCSResult fcs_mmm2d_require_virial(FCS handle, fcs_int flag);
FCSResult fcs_mmm2d_get_virial(FCS handle, fcs_float *virial);

FCSResult fcs_mmm2d_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched);
FCSResult fcs_mmm2d_print_parameters(FCS handle);

#endif
