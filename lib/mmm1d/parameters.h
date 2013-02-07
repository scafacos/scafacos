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

#ifndef _MMM1D_PARAMETERS_H
#define _MMM1D_PARAMETERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "FCSResult.h"


/** parameters for MMM1D. Most of the parameters can also be tuned automatically. Unlike
    P3M, this tuning is redone automatically whenever parameters change, but not immediately
    if you set this parameters.
    @param switch_rad at which xy-distance the calculation switches from the far to the
                      near formula. If -1, this parameter will be tuned automatically.
    @param bessel_cutoff the cutoff for the bessel sum, aka far formula. Normally set this
                         to -1, then the cutoff is automatically determined using the error formula.
    @param maxPWerror the maximal allowed error for the potential and the forces without the
                      prefactors, i. e. for the pure lattice 1/r-sum. */
int mmm1d_set_params(double switch_rad, int bessel_cutoff, double maxPWerror);

void mmm1d_set_far_switch_radius_2(void *rd, fcs_float rad2);
void mmm1d_get_far_switch_radius_2(void *rd, fcs_float *rad2);

void mmm1d_set_bessel_cutoff(void *rd, fcs_int cutoff);
void mmm1d_get_bessel_cutoff(void *rd, fcs_int *cutoff);

void mmm1d_set_maxPWerror(void *rd, fcs_float maxerr);
void mmm1d_get_maxPWerror(void *rd, fcs_float *maxerr);

void mmm1d_set_box_a(void *rd, fcs_float a);
void mmm1d_set_box_b(void *rd, fcs_float b);
void mmm1d_set_box_c(void *rd, fcs_float c);
#endif
