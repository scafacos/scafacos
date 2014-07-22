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



#ifndef _FCS_MMM1D_P_H
#define _FCS_MMM1D_P_H

#include "FCSDefinitions.h"
#include "FCSResult_p.h"
#include "FCSInterface_p.h"

/// * @brief file containing all mmm1d specific functions

typedef struct fcs_mmm1d_parameters_t *fcs_mmm1d_parameters;

FCSResult fcs_mmm1d_set_far_switch_radius(FCS handle, fcs_float rad);
FCSResult fcs_mmm1d_get_far_switch_radius(FCS handle, fcs_float *rad);

FCSResult fcs_mmm1d_set_bessel_cutoff(FCS handle, fcs_int cutoff);
FCSResult fcs_mmm1d_get_bessel_cutoff(FCS handle, fcs_int *cutoff);

FCSResult fcs_mmm1d_set_maxPWerror(FCS handle, fcs_float maxPWerror);
FCSResult fcs_mmm1d_get_maxPWerror(FCS handle, fcs_float *maxPWerror);

#endif
