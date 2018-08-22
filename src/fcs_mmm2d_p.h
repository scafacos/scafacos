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



#ifndef _FCS_MMM2D_P_H
#define _FCS_MMM2D_P_H

#include "fcs_definitions.h"
#include "fcs_result_p.h"
#include "fcs_interface_p.h"

/// * @brief file containing all mmm2d specific functions

FCSResult fcs_mmm2d_set_maxPWerror(FCS handle, fcs_float maxPWerror);
FCSResult fcs_mmm2d_get_maxPWerror(FCS handle, fcs_float *maxPWerror);

FCSResult fcs_mmm2d_set_far_cutoff(FCS handle, fcs_float cutoff);
FCSResult fcs_mmm2d_get_far_cutoff(FCS handle, fcs_float *cutoff);

FCSResult fcs_mmm2d_set_dielectric_contrasts(FCS handle, fcs_float delta_top, fcs_float delta_bot);
FCSResult fcs_mmm2d_get_dielectric_contrasts(FCS handle, fcs_float *delta_top, fcs_float *delta_bot);

FCSResult fcs_mmm2d_set_layers_per_node(FCS handle, fcs_int n_layers);
FCSResult fcs_mmm2d_get_layers_per_node(FCS handle, fcs_int *n_layers);

FCSResult fcs_mmm2d_set_skin(FCS handle, fcs_float skin);
FCSResult fcs_mmm2d_get_skin(FCS handle, fcs_float *skin);

FCSResult fcs_mmm2d_require_total_energy(FCS handle, fcs_int flag);
FCSResult fcs_mmm2d_get_total_energy(FCS handle, fcs_float *total_energy);

#endif
