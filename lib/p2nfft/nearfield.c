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
static fcs_float evaluate_cos_polynomial_1d(
   fcs_float x, fcs_int N, const fcs_float *coeff);
static fcs_float evaluate_sin_polynomial_1d(
   fcs_float x, fcs_int N, const fcs_float *coeff);


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
        d->interpolation_order, d->near_interpolation_num_nodes,
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

  if(d->interpolation_order >= 0){
    return ifcs_p2nfft_interpolation(
        dist, d->one_over_r_cut,
        d->interpolation_order, d->near_interpolation_num_nodes,
        d->near_interpolation_table_potential
      );
  } else if (d->cg_cos_coeff != NULL){
    return d->epsI/d->r_cut * evaluate_cos_polynomial_1d(dist * d->epsI/d->r_cut, d->N_cg_cos, d->cg_cos_coeff);
  } else if(d->taylor2p_coeff != NULL) {
    return ifcs_p2nfft_nearfield_correction_taylor2p(dist*d->one_over_r_cut, d->p, d->taylor2p_coeff)*d->one_over_r_cut;
  } else
    return 0.0;
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
          d->interpolation_order, d->near_interpolation_num_nodes,
          d->near_interpolation_table_potential);
  } else {
    return (1-fcs_erf(d->alpha * dist))/dist; /* use erf instead of erfc to fix ICC performance problems */
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
          d->interpolation_order, d->near_interpolation_num_nodes,
          d->near_interpolation_table_potential
        );
  } else if (d->cg_cos_coeff != NULL){
    return 1.0/dist - d->epsI/d->r_cut * evaluate_cos_polynomial_1d(dist * d->epsI/d->r_cut, d->N_cg_cos, d->cg_cos_coeff);
  } else if(d->taylor2p_coeff != NULL) {
//     fprintf(stderr, "compute near potential directly with taylor2p\n");
    return 1.0/dist - ifcs_p2nfft_nearfield_correction_taylor2p(dist*d->one_over_r_cut, d->p, d->taylor2p_coeff) * d->one_over_r_cut;
  } else
    return 0.0;
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
          d->interpolation_order, d->near_interpolation_num_nodes,
          d->near_interpolation_table_force);
  else
    return -((1-fcs_erf(d->alpha * dist))/dist
          + 2.0*d->alpha*FCS_P2NFFT_1_SQRTPI * exp(- d->alpha*d->alpha * dist*dist)
        ) / dist;
}

static fcs_float compute_near_field_nonperiodic(
    const void* param, fcs_float dist
    )
{
  ifcs_p2nfft_data_struct* d = (ifcs_p2nfft_data_struct*) param;
 
  if(d->interpolation_order >= 0){
    return -1.0/(dist*dist)
      - ifcs_p2nfft_interpolation(
          dist, d->one_over_r_cut,
          d->interpolation_order, d->near_interpolation_num_nodes,
          d->near_interpolation_table_force
        );
  } else if (d->cg_sin_coeff != NULL){
    fcs_float scale = d->epsI/d->r_cut;
    return -1.0/(dist*dist)
        - evaluate_sin_polynomial_1d(dist*scale, d->N_cg_cos, d->cg_sin_coeff) * scale * scale;
  } else if(d->taylor2p_coeff != NULL){
    return -1.0/(dist*dist)
        - ifcs_p2nfft_nearfield_correction_taylor2p_derive(dist*d->one_over_r_cut, d->p, d->taylor2p_derive_coeff) * d->one_over_r_cut * d->one_over_r_cut;
  } else
    return 0.0;
}


void ifcs_p2nfft_compute_near(
    const void* param, fcs_float dist,
    fcs_float *f, fcs_float *p
    )
{
  *p = ifcs_p2nfft_compute_near_potential(param, dist);
  *f = ifcs_p2nfft_compute_near_field(param, dist);
}

static fcs_float evaluate_cos_polynomial_1d(
   fcs_float x, fcs_int N, const fcs_float *coeff
   )
{
  fcs_float value=0;

  for(int k=0; k<N; k++)
    value += coeff[k] * cos(2*FCS_P2NFFT_PI*k*x);

  return value;
}

static fcs_float evaluate_sin_polynomial_1d(
   fcs_float x, fcs_int N, const fcs_float *coeff
   )
{
  fcs_float value=0;

  for(int k=0; k<N; k++)
    value += coeff[k] * sin(2*FCS_P2NFFT_PI*k*x);

  return value;
}

