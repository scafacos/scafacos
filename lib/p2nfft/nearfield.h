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

#ifndef _P2NFFT_NEARFIELD_H
#define _P2NFFT_NEARFIELD_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "FCSCommon.h"
#include "constants.h"
#include "regularization.h"


typedef struct {
  fcs_int interpolation_order;      /* interpolation order */
  fcs_int interpolation_num_nodes;  /* number of sampled points */
  fcs_float *near_interpolation_table_potential; /* nearfield potential values */
  fcs_float *near_interpolation_table_force;     /* nearfield force values */
  fcs_float one_over_r_cut;
} ifcs_p2nfft_near_params;


fcs_float ifcs_p2nfft_compute_self_potential(
    const void* param);
fcs_float ifcs_p2nfft_compute_near_potential(
    const void* param, fcs_float dist);
fcs_float ifcs_p2nfft_compute_near_field(
    const void* param, fcs_float dist);
void ifcs_p2nfft_compute_near(
    const void* param, fcs_float dist, 
    fcs_float *potential, fcs_float *field);


/* callback functions for near field computations */
static inline fcs_float
ifcs_p2nfft_compute_near_potential_periodic_erfc(
    const void *param, fcs_float dist
    )
{
  fcs_float alpha = *((fcs_float *) param);
  fcs_float adist = alpha * dist;
  return (1.0 - erf(adist)) / dist; /* use erf instead of erfc to fix ICC performance problems */
}

static inline fcs_float
ifcs_p2nfft_compute_near_field_periodic_erfc(
    const void *param, fcs_float dist
    )
{
  fcs_float alpha = *((fcs_float *) param);
  fcs_float inv_dist = 1.0/dist;
  fcs_float adist = alpha * dist;
  fcs_float erfc_part_ri = (1.0 - erf(adist)) * inv_dist; /* use erf instead of erfc to fix ICC performance problems */
  return -(erfc_part_ri + 2.0*alpha*FCS_P2NFFT_1_SQRTPI*exp(-adist*adist)) * inv_dist;
}


static inline void
ifcs_p2nfft_compute_near_periodic_erfc(
    const void *param, fcs_float dist,
    fcs_float *f, fcs_float *p
    )
{
  fcs_float alpha = *((fcs_float *) param);
  fcs_float inv_dist = 1.0/dist;
  fcs_float adist = alpha * dist;
  fcs_float erfc_part_ri = (1.0 - erf(adist)) * inv_dist; /* use erf instead of erfc to fix ICC performance problems */
  *p = erfc_part_ri;
  *f = -(erfc_part_ri + 2.0*alpha*FCS_P2NFFT_1_SQRTPI*exp(-adist*adist)) * inv_dist;
}

static inline fcs_float
ifcs_p2nfft_approx_erfc(
    fcs_float adist
    )
{
  /* approximate \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying formula 7.1.26 from:
     Abramowitz/Stegun: Handbook of Mathematical Functions, Dover (9. ed.), chapter 7.
     Error <= 1.5e-7 */
  fcs_float t = 1.0 / (1.0 + 0.3275911 * adist);
  return fcs_exp(-adist*adist) * 
    (t * (0.254829592 + 
	  t * (-0.284496736 + 
	       t * (1.421413741 + 
		    t * (-1.453152027 + 
			 t * 1.061405429)))));
} 

static inline fcs_float
ifcs_p2nfft_compute_near_potential_periodic_approx_erfc(
    const void *param, fcs_float dist
    )
{
  fcs_float alpha = *((fcs_float *) param);
  fcs_float adist = alpha * dist;
  return ifcs_p2nfft_approx_erfc(adist) / dist;
}

static inline fcs_float
ifcs_p2nfft_compute_near_field_periodic_approx_erfc(
    const void *param, fcs_float dist
    )
{
  fcs_float alpha = *((fcs_float *) param);
  fcs_float inv_dist = 1.0/dist;
  fcs_float adist = alpha * dist;
  fcs_float erfc_part_ri = ifcs_p2nfft_approx_erfc(adist) *inv_dist;
  return -(erfc_part_ri + 2.0*alpha*FCS_P2NFFT_1_SQRTPI*exp(-adist*adist)) * inv_dist;
}

static inline void
ifcs_p2nfft_compute_near_periodic_approx_erfc(
    const void *param, fcs_float dist,
    fcs_float *f, fcs_float *p
    )
{
  fcs_float alpha = *((fcs_float *) param);
  fcs_float inv_dist = 1.0/dist;
  fcs_float adist = alpha * dist;
  fcs_float erfc_part_ri = ifcs_p2nfft_approx_erfc(adist) *inv_dist;

  *p = erfc_part_ri;
  *f = -(erfc_part_ri + 2.0*alpha*FCS_P2NFFT_1_SQRTPI*exp(-adist*adist)) * inv_dist;
}


static inline fcs_float
ifcs_p2nfft_compute_near_potential_interpolation(
    const void *param, fcs_float dist
    )
{
  ifcs_p2nfft_near_params *d = (ifcs_p2nfft_near_params*) param;

  return 1.0/dist - ifcs_p2nfft_interpolation(
      dist, d->one_over_r_cut,
      d->interpolation_order, d->interpolation_num_nodes,
      d->near_interpolation_table_potential);
}

static inline fcs_float
ifcs_p2nfft_compute_near_field_interpolation(
    const void *param, fcs_float dist
    )
{
  ifcs_p2nfft_near_params *d = (ifcs_p2nfft_near_params*) param;
  fcs_float inv_dist = 1.0/dist;

  return -inv_dist*inv_dist - ifcs_p2nfft_interpolation(
        dist, d->one_over_r_cut,
        d->interpolation_order, d->interpolation_num_nodes,
        d->near_interpolation_table_force);
}

static inline void
ifcs_p2nfft_compute_near_interpolation(
    const void *param, fcs_float dist,
    fcs_float *f, fcs_float *p
    )
{
  ifcs_p2nfft_near_params *d = (ifcs_p2nfft_near_params*) param;
  fcs_float inv_dist = 1.0/dist;

  *p = inv_dist - ifcs_p2nfft_interpolation(
      dist, d->one_over_r_cut,
      d->interpolation_order, d->interpolation_num_nodes,
      d->near_interpolation_table_potential);
  *f = -inv_dist*inv_dist - ifcs_p2nfft_interpolation(
        dist, d->one_over_r_cut,
        d->interpolation_order, d->interpolation_num_nodes,
        d->near_interpolation_table_force);
}



/** linear spline interpolation in near field with even kernels */
static inline fcs_float ifcs_p2nfft_intpol_even_const(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c;
  fcs_int r;
  fcs_float f1;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  f1=table[r];
  return f1;
}

/** linear spline interpolation in near field with even kernels */
static inline fcs_float ifcs_p2nfft_intpol_even_lin(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c,c1,c3;
  fcs_int r;
  fcs_float f1,f2;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  f1=table[r];f2=table[r+1];
  c1=c-r;
  c3=c1-1.0;
  return (-f1*c3+f2*c1);
}

/** quadratic spline interpolation in near field with even kernels */
static inline fcs_float ifcs_p2nfft_intpol_even_quad(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c,c1,c2,c3;
  fcs_int r;
  fcs_float f0,f1,f2;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  if (r==0) {f0=table[r+1];f1=table[r];f2=table[r+1];}
  else { f0=table[r-1];f1=table[r];f2=table[r+1];}
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  return (0.5*f0*c1*c3-f1*c2*c3+0.5*f2*c2*c1);
}

/** cubic spline interpolation in near field with even kernels */
static inline fcs_float ifcs_p2nfft_intpol_even_cub(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c,c1,c2,c3,c4;
  fcs_int r;
  fcs_float f0,f1,f2,f3;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  if (r==0) {f0=table[r+1];f1=table[r];f2=table[r+1];f3=table[r+2];}
  else { f0=table[r-1];f1=table[r];f2=table[r+1];f3=table[r+2];}
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  c4=c1-2.0;
  return(-f0*c1*c3*c4+3.0*f1*c2*c3*c4-3.0*f2*c2*c1*c4+f3*c2*c1*c3)/6.0;
}

/* define one inline near field evaluation function for each interpolation order */
#define FCS_P2NFFT_NEAR_INTPOL_FUNC(_suffix_) \
  static inline void \
  ifcs_p2nfft_compute_near_interpolation_ ## _suffix_( \
      const void *param, fcs_float dist, \
      fcs_float *f, fcs_float *p \
      ) \
  { \
    ifcs_p2nfft_near_params *d = (ifcs_p2nfft_near_params*) param; \
    fcs_float inv_dist = 1.0/dist; \
    *p = inv_dist - ifcs_p2nfft_intpol_even_ ## _suffix_( \
        dist, d->near_interpolation_table_potential, \
        d->interpolation_num_nodes, d->one_over_r_cut); \
    *f = -inv_dist*inv_dist - ifcs_p2nfft_intpol_even_ ## _suffix_( \
        dist, d->near_interpolation_table_force, \
        d->interpolation_num_nodes, d->one_over_r_cut); \
  }

FCS_P2NFFT_NEAR_INTPOL_FUNC(const)
FCS_P2NFFT_NEAR_INTPOL_FUNC(lin)
FCS_P2NFFT_NEAR_INTPOL_FUNC(quad)
FCS_P2NFFT_NEAR_INTPOL_FUNC(cub)

#endif
