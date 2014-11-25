/*
  Copyright (C) 2011 Olaf Lenz
  
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

#include "tune.h"
#include <stdlib.h>
// #include <stdio.h>

#include "types.h"

#include <math.h>

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/
static fcs_float mmm1d_far_error(mmm1d_data_struct *d, fcs_int P, fcs_float minrad);
static fcs_int mmm1d_determine_minrad(mmm1d_data_struct *d, fcs_int P);
static void mmm1d_determine_bessel_radii(mmm1d_data_struct *d);
static void mmm1d_prepare_polygamma_series(mmm1d_data_struct *d);
static FCSResult mmm1d_check_system_charges(mmm1d_data_struct *d,
                       fcs_int num_particles, fcs_float *charges);

/***************************************************/
/* IMPLEMENTATION */
/***************************************************/
FCSResult mmm1d_tune(void* rd,
        fcs_int num_particles,
        fcs_float *positions,
        fcs_float *charges) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  const char* fnc_name = "mmm1d_tune";
  
  /* check whether the input parameters are sane */
  /* set switch radius to usually good value */
  if (d->far_switch_radius_2 < 0) {
     fcs_float switch_rad = MMM1D_DEFAULT_FAR_SWITCH_RADIUS*d->box_l[2];
     d->far_switch_radius_2 = switch_rad*switch_rad;
  }

  /* check whether cutoff too large */
  fcs_float rboxz=d->box_l[2]*d->box_l[2];
  
  if (d->far_switch_radius_2 >= rboxz)
    d->far_switch_radius_2 = 0.8*rboxz;

  d->uz  = 1/d->box_l[2];
  d->L2  = d->box_l[2]*d->box_l[2];
  d->uz2 = d->uz*d->uz;
  d->L3_i = d->uz2*d->uz;
  
  mmm1d_determine_bessel_radii(d);
  mmm1d_prepare_polygamma_series(d);
  
  if (d->far_switch_radius_2 <=
      d->bessel_radii[d->bessel_cutoff - 1]*
      d->bessel_radii[d->bessel_cutoff - 1]) {
    // this switching radius is too small for our Bessel series
    return fcs_result_create(FCS_ERROR_LOGICAL_ERROR, fnc_name, "could not tune far formula to require accuracy. Increase far switching radius.");
  }
  
  return NULL;

}

static fcs_float mmm1d_far_error(mmm1d_data_struct *d, fcs_int P, fcs_float minrad)
{
  // this uses an upper bound to all force components and the potential
  fcs_float rhores = 2*M_PI*d->uz*minrad;
  fcs_float pref = 4*d->uz*mmm_dmax(1, 2*M_PI*d->uz);

  return pref*mmm_K1(rhores*P)*exp(rhores)/rhores*(P - 1 + 1/rhores);
}

static fcs_int mmm1d_determine_minrad(mmm1d_data_struct *d, fcs_int P)
{
  // bisection to search for where the error is maxPWerror
  fcs_float rgranularity = 0.01*d->box_l[2];
  fcs_float rmin = rgranularity;
  fcs_float rmax = mmm_dmin(d->box_l[0], d->box_l[1]);
  fcs_float errmin = mmm1d_far_error(d, P, rmin);
  fcs_float errmax = mmm1d_far_error(d, P, rmax);
  if (errmin < d->maxPWerror) {
    // we can do almost all radii with this P
    return rmin;
  }
  if (errmax > d->maxPWerror) {
    // make sure that this switching radius cannot be reached
    return 2*mmm_dmax(d->box_l[0], d->box_l[1]);
  }

  while (rmax - rmin > rgranularity) {
    fcs_float c = 0.5*(rmin + rmax);
    fcs_float errc = mmm1d_far_error(d, P, c);
    if (errc > d->maxPWerror) {
      rmin = c;
    } else {
      rmax = c;
    }
  }
  return 0.5*(rmin + rmax);
}

static void mmm1d_determine_bessel_radii(mmm1d_data_struct *d)
{
  d->bessel_radii = realloc(d->bessel_radii, sizeof(fcs_float)*d->bessel_cutoff);
  for (fcs_int P = 1; P <= d->bessel_cutoff; ++P) {
    d->bessel_radii[P-1] = mmm1d_determine_minrad(d, P);
    // printf("bessel cutoff %d %f\n", P, d->bessel_radii[P-1]);
  }
}

static void mmm1d_prepare_polygamma_series(mmm1d_data_struct *d)
{
  /* polygamma, determine order */
  fcs_int n;
  fcs_float err;
  fcs_float rhomax2nm2, rhomax2 = d->uz2*d->far_switch_radius_2;
  n = 1;
  rhomax2nm2 = 1.0;
  do {
    mmm_create_mod_psi_up_to(d->polTaylor, n+1);
    err = 2*n*fabs(mmm_mod_psi_even(d->polTaylor, n, 0.5))*rhomax2nm2;
    //printf("rank %d, recalc 3 %d, err: %e\n", d->comm.rank, n, err);
    rhomax2nm2 *= rhomax2;
    n++;
  }
  while (err > 0.1*d->maxPWerror);
  // printf("rank %d, rhomax2: %e, recalc_tables: n_modPsi: %d\n", d->comm.rank, rhomax2, (d->polTaylor)->n_modPsi);
}
