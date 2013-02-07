/*
  Copyright (C) 2011
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "parameters.h"
#include "types.h"
#include "FCSCommon.h"
#include <stdio.h>

///@TODO: condition to force retuning

void mmm1d_set_far_switch_radius_2(void *rd, fcs_float rad2) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  if (!fcs_float_is_equal(rad2, d->far_switch_radius_2))
    d->bessel_calculated = 1;
  d->far_switch_radius_2 = rad2;
  printf("set radius %e\n",d->far_switch_radius_2);
}

void mmm1d_get_far_switch_radius_2(void *rd, fcs_float *rad2) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  *rad2 = d->far_switch_radius_2;
}

void mmm1d_set_bessel_cutoff(void *rd, fcs_int cutoff) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  if (!fcs_float_is_equal(cutoff, d->bessel_cutoff))
    d->bessel_calculated = 1;
  d->bessel_cutoff = cutoff;
}

void mmm1d_get_bessel_cutoff(void *rd, fcs_int *cutoff) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  *cutoff = d->bessel_cutoff;
}

void mmm1d_set_maxPWerror(void *rd, fcs_float maxPWerror) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  if (!fcs_float_is_equal(maxPWerror, d->maxPWerror))
    d->bessel_calculated = 1;
  d->maxPWerror = maxPWerror;
}

void mmm1d_get_maxPWerror(void *rd, fcs_float *maxPWerror) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  *maxPWerror = d->maxPWerror;
}

void mmm1d_set_box_a(void* rd, fcs_float a) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  if (!fcs_float_is_equal(a, d->box_l[0]))
    d->needs_retune = 1;
  d->box_l[0] = a;
}

void mmm1d_set_box_b(void* rd, fcs_float b) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  if (!fcs_float_is_equal(b, d->box_l[1]))
    d->needs_retune = 1;
  d->box_l[1] = b;
}

void mmm1d_set_box_c(void* rd, fcs_float c) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  if (!fcs_float_is_equal(c, d->box_l[2]))
    d->needs_retune = 1;
  d->box_l[2] = c;
}
