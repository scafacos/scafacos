/*
 * Copyright (C) 2011-2013 Michael Pippig
 *
 * This file is part of ScaFaCoS.
 * 
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *	
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nearfield.h"
#include "kernels.h"
#include "types.h"

#define FCS_P2NFFT_NEAR_BASISPOLY 0

static fcs_float compute_self_potential_periodic(
    const void* param);
static fcs_float compute_self_potential_nonperiodic(
    const void* param);
static fcs_float compute_near_potential_periodic(
    const void* param, fcs_float dist);
static fcs_float compute_near_potential_nonperiodic(
    const void* param, fcs_float dist);
static fcs_float compute_near_field_periodic(
    const void* param, fcs_float dist);
static fcs_float compute_near_field_nonperiodic(
    const void* param, fcs_float dist);

fcs_float ifcs_p2nfft_compute_self_potential(
    const void* param
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) param;
  
  if(d->use_ewald)
    return compute_self_potential_periodic(param);
  else
    return compute_self_potential_nonperiodic(param);
}

static fcs_float compute_self_potential_periodic(
    const void* param
    )
{
  fcs_float dist = 0.0;
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) param;

  if(d->interpolation_order >= 0)
    return ifcs_p2nfft_interpolation(
        dist, d->one_over_r_cut,
        d->interpolation_order, d->interpolation_num_nodes,
        d->near_interpolation_table_potential);
  else
    return 2 * d->alpha * FCS_P2NFFT_1_SQRTPI;
}

static fcs_float compute_self_potential_nonperiodic(
    const void* param
    )
{
  fcs_float dist=0.0;
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) param;

  if(d->interpolation_order >= 0)
    return ifcs_p2nfft_interpolation(
        dist, d->one_over_r_cut,
        d->interpolation_order, d->interpolation_num_nodes,
        d->near_interpolation_table_potential
      );
  else {
#if FCS_P2NFFT_NEAR_BASISPOLY
    return ifcs_p2nfft_regkernel(ifcs_p2nfft_one_over_modulus, dist/d->box_scales[0], d->p, 0, d->epsI, d->epsB) / d->box_scales[0];
#else
    return ifcs_p2nfft_nearfield_correction_taylor2p(dist/d->box_scales[0], d->p, d->taylor2p_coeff) / d->box_scales[0];
#endif
  }
}



fcs_float ifcs_p2nfft_compute_near_potential(
    const void* param, fcs_float dist
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) param;
  
  if(d->use_ewald)
    return compute_near_potential_periodic(param, dist);
  else
    return compute_near_potential_nonperiodic(param, dist);
}

static fcs_float compute_near_potential_periodic(
    const void* param, fcs_float dist
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) param;
  if(d->interpolation_order >= 0){
    return 1.0/dist
      - ifcs_p2nfft_interpolation(
          dist, d->one_over_r_cut,
          d->interpolation_order, d->interpolation_num_nodes,
          d->near_interpolation_table_potential);
  } else {
    return (1-erf(d->alpha * dist))/dist; /* use erf instead of erfc to fix ICC performance problems */
  }
}

static fcs_float compute_near_potential_nonperiodic(
    const void* param, fcs_float dist
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) param;

  if(d->interpolation_order >= 0){
    return 1.0/dist
      - ifcs_p2nfft_interpolation(
          dist, d->one_over_r_cut,
          d->interpolation_order, d->interpolation_num_nodes,
          d->near_interpolation_table_potential
        );
  } else {
#if FCS_P2NFFT_NEAR_BASISPOLY
    return 1.0/dist - ifcs_p2nfft_regkernel(ifcs_p2nfft_one_over_modulus, dist/d->box_scales[0], d->p, 0, d->epsI, d->epsB) / d->box_scales[0];
#else
    return 1.0/dist - ifcs_p2nfft_nearfield_correction_taylor2p(dist/d->box_scales[0], d->p, d->taylor2p_coeff) / d->box_scales[0];
#endif
  }
}

/* f(r) = 1/r 
 * grad f(r) = 1/r * (x,y,z)^T * f'(r) 
 * Do not include the factor 1/r * (x,y,z)^T, it will be calculated later. */
fcs_float ifcs_p2nfft_compute_near_field(
    const void* param, fcs_float dist
    )
{
  ifcs_p2nfft_data_struct* d = (ifcs_p2nfft_data_struct*) param;

  if(d->use_ewald)
    return compute_near_field_periodic(param, dist);
  else
    return compute_near_field_nonperiodic(param, dist);
}

static fcs_float compute_near_field_periodic(
    const void* param, fcs_float dist
    )
{
  ifcs_p2nfft_data_struct* d = (ifcs_p2nfft_data_struct*) param;
  
  if(d->interpolation_order >=0 )
    return -1.0/(dist*dist)
      - ifcs_p2nfft_interpolation(
          dist, d->one_over_r_cut,
          d->interpolation_order, d->interpolation_num_nodes,
          d->near_interpolation_table_force);
  else
    return -((1-erf(d->alpha * dist))/dist
          + 2.0*d->alpha*FCS_P2NFFT_1_SQRTPI * exp(- d->alpha*d->alpha * dist*dist)
        ) / dist;
}

static fcs_float compute_near_field_nonperiodic(
    const void* param, fcs_float dist
    )
{
  ifcs_p2nfft_data_struct* d = (ifcs_p2nfft_data_struct*) param;
 
  if(d->interpolation_order >= 0)
    return -1.0/(dist*dist)
      - ifcs_p2nfft_interpolation(
          dist, d->one_over_r_cut,
          d->interpolation_order, d->interpolation_num_nodes,
          d->near_interpolation_table_force
        );
  else
    return -1.0/(dist*dist)
        - ifcs_p2nfft_nearfield_correction_taylor2p_derive(
          dist/d->box_scales[0], d->p, d->taylor2p_derive_coeff
        ) / (d->box_scales[0]*d->box_scales[0]);
}


void ifcs_p2nfft_compute_near(
    const void* param, fcs_float dist,
    fcs_float *f, fcs_float *p
    )
{
  *p = ifcs_p2nfft_compute_near_potential(param, dist);
  *f = ifcs_p2nfft_compute_near_field(param, dist);
}

