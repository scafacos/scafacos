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



#ifndef _FCS_EWALD_P_H
#define _FCS_EWALD_P_H

#include "FCSDefinitions.h"
#include "FCSResult_p.h"
#include "FCSInterface_p.h"

#ifdef __cplusplus
extern "C" {
#endif
  
FCSResult fcs_ewald_set_kmax(FCS handle, fcs_int kmax);
FCSResult fcs_ewald_set_kmax_tune(FCS handle);
FCSResult fcs_ewald_get_kmax(FCS handle, fcs_int *kmax);

FCSResult fcs_ewald_set_maxkmax(FCS handle, fcs_int maxkmax);
FCSResult fcs_ewald_set_maxkmax_tune(FCS handle);
FCSResult fcs_ewald_get_maxkmax(FCS handle, fcs_int *maxkmax);

FCSResult fcs_ewald_set_r_cut(FCS handle, fcs_float r_cut);
FCSResult fcs_ewald_set_r_cut_tune(FCS handle);
FCSResult fcs_ewald_get_r_cut(FCS handle, fcs_float *r_cut);

FCSResult fcs_ewald_set_alpha(FCS handle, fcs_float alpha);
FCSResult fcs_ewald_set_alpha_tune(FCS handle);
FCSResult fcs_ewald_get_alpha(FCS handle, fcs_float *alpha);

FCSResult fcs_ewald_set_tolerance_field(FCS handle, fcs_float tolerance_field_abs);
FCSResult fcs_ewald_set_tolerance_field_tune(FCS handle);
FCSResult fcs_ewald_get_tolerance_field(FCS handle, fcs_float* tolerance_field_abs);

#ifdef __cplusplus
}
#endif

#endif
