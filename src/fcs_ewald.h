/*
  Copyright (C) 2011-2012 Olaf Lenz
  
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



#ifndef _FCS_EWALD_H
#define _FCS_EWALD_H

#include "fcs_ewald_p.h"
#include "FCSResult.h"
#include "FCSInterface.h"

FCSResult fcs_ewald_init(FCS handle);

FCSResult fcs_ewald_destroy(FCS handle);

FCSResult fcs_ewald_tune(FCS handle,
			 fcs_int num_particles,
			 fcs_float *positions, 
			 fcs_float *charges);

FCSResult fcs_ewald_run(FCS handle,
			fcs_int num_particles,
			fcs_float *positions, 
			fcs_float *charges,
			fcs_float *fields,
			fcs_float *potentials);

FCSResult fcs_ewald_set_tolerance(FCS handle, fcs_int tolerance_type, fcs_float tolerance);
FCSResult fcs_ewald_get_tolerance(FCS handle, fcs_int *tolerance_type, fcs_float *tolerance);

FCSResult fcs_ewald_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched);
FCSResult fcs_ewald_print_parameters(FCS handle);

#endif
