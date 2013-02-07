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

#ifndef _MMM2D_PARAMETERS_H
#define _MMM2D_PARAMETERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "FCSResult.h"

void mmm2d_set_far_cutoff(void *rd, fcs_float cutoff);
void mmm2d_get_far_cutoff(void *rd, fcs_float *cutoff);

void mmm2d_set_dielectric_contrasts(void *rd, fcs_float delta_top, fcs_float delta_bot);
void mmm2d_get_dielectric_contrasts(void *rd, fcs_float *delta_top, fcs_float *delta_bot);

void mmm2d_set_maxPWerror(void *rd, fcs_float maxerr);
void mmm2d_get_maxPWerror(void *rd, fcs_float *maxerr);

void mmm2d_set_box_a(void *rd, fcs_float a);
void mmm2d_set_box_b(void *rd, fcs_float b);
void mmm2d_set_box_c(void *rd, fcs_float c);

void mmm2d_set_layers_per_node(void *rd, fcs_int n_layers);
void mmm2d_get_layers_per_node(void *rd, fcs_int *n_layers);

void mmm2d_require_total_energy(void *rd, fcs_int flag);
FCSResult mmm2d_get_total_energy(void *rd, fcs_float *total_energy);

void mmm2d_set_skin(void *rd, fcs_float skin);
void mmm2d_get_skin(void *rd, fcs_float *skin);
#endif
