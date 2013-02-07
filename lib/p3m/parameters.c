/*
  Copyright (C) 2011,2012 Olaf Lenz
  
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

void ifcs_p3m_set_near_field_flag(void *rd, fcs_int flag) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  d->near_field_flag = flag;
}

void ifcs_p3m_set_box_a(void* rd, fcs_float a) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  if (!fcs_float_is_equal(a, d->box_l[0]))
    d->needs_retune = 1;
  d->box_l[0] = a;
}

void ifcs_p3m_set_box_b(void *rd, fcs_float b) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  if (!fcs_float_is_equal(b, d->box_l[1]))
    d->needs_retune = 1;
  d->box_l[1] = b;
}

void ifcs_p3m_set_box_c(void *rd, fcs_float c) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  if (!fcs_float_is_equal(c, d->box_l[2]))
    d->needs_retune = 1;
  d->box_l[2] = c;
}

void ifcs_p3m_set_r_cut(void *rd, fcs_float r_cut) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  if (!fcs_float_is_equal(r_cut, d->r_cut))
    d->needs_retune = 1;
  d->r_cut = r_cut;
  d->tune_r_cut = 0;
}

void ifcs_p3m_set_r_cut_tune(void *rd) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  d->needs_retune = 1;
  d->tune_r_cut = 1;
}

void ifcs_p3m_get_r_cut(void *rd, fcs_float *r_cut) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  *r_cut = d->r_cut;
}

void ifcs_p3m_set_alpha(void *rd, fcs_float alpha) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  if (!fcs_float_is_equal(alpha, d->alpha))
    d->needs_retune = 1;
  d->alpha = alpha;
  d->tune_alpha = 0;
}

void ifcs_p3m_set_alpha_tune(void *rd) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  d->needs_retune = 1;
  d->tune_alpha = 1;
}

void ifcs_p3m_get_alpha(void *rd, fcs_float *alpha) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  *alpha = d->alpha;
}

void ifcs_p3m_set_grid(void *rd, fcs_int grid) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  if (!(grid == d->grid[0] && grid == d->grid[1] && grid == d->grid[2]))
    d->needs_retune = 1;
  d->grid[0] = grid;
  d->grid[1] = grid;
  d->grid[2] = grid;
  d->tune_grid = 0;
}

void ifcs_p3m_set_grid_tune(void *rd) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  d->needs_retune = 1;
  d->tune_grid = 1;
}

void ifcs_p3m_get_grid(void *rd, fcs_int *grid) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  grid[0] = d->grid[0];
  grid[1] = d->grid[1];
  grid[2] = d->grid[2];
}

void ifcs_p3m_set_cao(void *rd, fcs_int cao) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  if (cao != d->cao)
    d->needs_retune = 1;
  d->cao = cao;
  d->tune_cao = 0;
}

void ifcs_p3m_set_cao_tune(void *rd) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  d->needs_retune = 1;
  d->tune_cao = 1;
}

void ifcs_p3m_get_cao(void *rd, fcs_int *cao) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  *cao = d->cao;
}

void ifcs_p3m_set_tolerance_field(void *rd, fcs_float tolerance_field) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  if (!fcs_float_is_equal(tolerance_field, d->tolerance_field))
    d->needs_retune = 1;
  d->tolerance_field = tolerance_field;
}

void ifcs_p3m_set_tolerance_field_tune(void *rd) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  d->needs_retune = 1;
  d->tolerance_field = P3M_DEFAULT_TOLERANCE_FIELD;
}  

void ifcs_p3m_get_tolerance_field(void *rd, fcs_float* tolerance_field) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  *tolerance_field = d->tolerance_field;
}

void ifcs_p3m_require_total_energy(void *rd, fcs_int flag) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  d->require_total_energy = flag;
}

FCSResult ifcs_p3m_get_total_energy(void *rd, fcs_float *total_energy) {
  const char* fnc_name = "ifcs_p3m_get_total_energy";
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  if (d->require_total_energy) {
    *total_energy = d->total_energy;
    return NULL;
  } else 
    return 
      fcsResult_create
      (FCS_LOGICAL_ERROR, fnc_name, 
       "Trying to get total energy, but computation was not requested.");
}
