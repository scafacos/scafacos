/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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

void mmm2d_set_far_cutoff(void *rd, fcs_float cutoff) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  if (!fcs_float_is_equal(cutoff, d->far_cut)) {
    d->calculate_far = 0;
    d->far_cut = cutoff;
    d->needs_tuning = 1;
  }
}

void mmm2d_get_far_cutoff(void *rd, fcs_float *cutoff) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  *cutoff = d->far_cut;
}

void mmm2d_set_dielectric_contrasts(void *rd, fcs_float delta_top, fcs_float delta_bot) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  if (!fcs_float_is_equal(delta_top, d->delta_mid_top) || !fcs_float_is_equal(delta_bot, d->delta_mid_bot)) {
    if (delta_top != 0.0 || delta_bot != 0.0) {
      d->dielectric_contrast_on=1;
      d->delta_mid_top= delta_top;
      d->delta_mid_bot= delta_bot;
    } else {
      d->dielectric_contrast_on=0;
      d->delta_mid_top= 0.;
      d->delta_mid_bot= 0.;
    }
  }
}

void mmm2d_get_dielectric_contrasts(void *rd, fcs_float *delta_top, fcs_float *delta_bot) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  *delta_top = d->delta_mid_top;
  *delta_bot = d->delta_mid_bot;
}

void mmm2d_set_maxPWerror(void *rd, fcs_float maxPWerror) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  if (!fcs_float_is_equal(maxPWerror, d->maxPWerror))
    d->needs_tuning = 1;
  d->maxPWerror = maxPWerror;
}

void mmm2d_get_maxPWerror(void *rd, fcs_float *maxPWerror) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  *maxPWerror = d->maxPWerror;
}

void mmm2d_set_box_a(void* rd, fcs_float a) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  if (!fcs_float_is_equal(a, d->box_l[0]))
    d->needs_tuning = 1;
  d->box_l[0] = a;
  d->box_l_i[0] = 1./a;
}

void mmm2d_set_box_b(void* rd, fcs_float b) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  if (!fcs_float_is_equal(b, d->box_l[1]))
    d->needs_tuning = 1;
  d->box_l[1] = b;
  d->box_l_i[1] = 1./b;
}

void mmm2d_set_box_c(void* rd, fcs_float c) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  if (!fcs_float_is_equal(c, d->box_l[2]))
    d->needs_tuning = 1;
  d->box_l[2] = c;
  d->box_l_i[2] = 1./c;
}

void mmm2d_set_layers_per_node(void *rd, fcs_int n_layers) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  if (n_layers!=d->layers_per_node)
    d->needs_tuning = 1;
  d->layers_per_node=n_layers;
  d->n_total_layers=d->comm.size*n_layers+2;
}

void mmm2d_require_total_energy(void *rd, fcs_int flag) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  d->require_total_energy = flag;
}

FCSResult mmm2d_get_total_energy(void *rd, fcs_float *total_energy) {
  const char* fnc_name = "mmm2d_get_total_energy";
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  if (d->require_total_energy) {
    *total_energy = d->total_energy;
    return NULL;
  } else 
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Trying to get total energy, but computation was not requested.");
}

void mmm2d_get_layers_per_node(void *rd, fcs_int *n_layers) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  *n_layers=d->layers_per_node;
}

void mmm2d_set_skin(void *rd, fcs_float skin) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  if (!fcs_float_is_equal(skin, d->skin))
    d->needs_tuning = 1;
  d->skin=skin;
}

void mmm2d_get_skin(void *rd, fcs_float *skin) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  *skin=d->skin;
}
