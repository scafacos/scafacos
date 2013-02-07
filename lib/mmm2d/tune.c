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

#include "tune.h"
#include <stdlib.h>
#include <stdio.h>

#include "types.h"

#include <math.h>

FCSResult mmm2d_setup_constants(mmm2d_data_struct *d);
FCSResult mmm2d_tune_near(mmm2d_data_struct *d);
FCSResult mmm2d_tune_far(mmm2d_data_struct *d);
static void mmm2d_prepareBernoulliNumbers(mmm2d_data_struct *d);
static FCSResult mmm2d_check_system_charges(mmm2d_data_struct *d,
                fcs_int num_particles, fcs_float *charges);

FCSResult mmm2d_tune(void* rd,
        fcs_int num_particles,
        fcs_float *positions,
        fcs_float *charges) {
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  const char* fnc_name = "mmm2d_tune";
  FCSResult res;
  
  /* Check for charge existence and neutrality */
  mmm2d_check_system_charges(d, num_particles, charges);
  if (!fcs_float_is_zero(d->total_charge))
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "MMM2D requires a zero net charge.");
  
  d->my_bottom = d->comm.rank*d->box_l[2]/(fcs_float)(d->comm.size);
  
  /* Broadcast charges */
  /*
  MPI_Bcast(&num_particles, 1, FCS_MPI_INT, 0, d->comm.mpicomm);
  MPI_Bcast(charges, num_particles, FCS_MPI_FLOAT, 0, d->comm.mpicomm);
  MPI_Bcast(positions, 3*num_particles, FCS_MPI_FLOAT, 0, d->comm.mpicomm);
  */
  
  ///@TODO: skip next steps if d->needs_tuning==0. Check if this is flawless
  if(d->needs_tuning==0)
    return NULL;
  
  /* precalculate some constants */
  mmm2d_setup_constants(d);
  
  /* tune near formula */
  res=mmm2d_tune_near(d);
  if (res) return res;
  
  /* tune far formula */
  if (d->comm.size*d->layers_per_node < 3) {
    d->far_cut = 0.0;
    if (d->dielectric_contrast_on)
      return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Definition of dielectric contrasts requires more than 3 layers");
  } else {
      res=mmm2d_tune_far(d);
      if (res) return res;
  }
  
  d->needs_tuning=0;
  
  /*
  printf("node %d, mmm2d_tune\n", d->comm.rank);
  printf("node %d: ux: %f, uy: %f\n", d->comm.rank, d->ux, d->uy);
  printf("node %d: min_far: %f, far_cut: %f\n", d->comm.rank, d->min_far, d->far_cut);
  printf("node %d: layer_h %f\n", d->comm.rank, d->layer_h);
  printf("node %d: n_layers %d\n", d->comm.rank, d->n_layers);
  */
  
  return NULL;
}

FCSResult mmm2d_setup_constants(mmm2d_data_struct *d) {
  d->ux  = 1./d->box_l[0];
  d->ux2 = d->ux*d->ux;
  d->uy  = 1./d->box_l[1];
  d->uy2 = d->uy*d->uy;  
  d->uz  = 1./d->box_l[2];
  d->layer_h=d->box_l[2]/(d->comm.size*d->layers_per_node);
  d->min_far=d->layer_h - d->skin;
  d->max_near = 2.*d->layer_h + d->skin;
  return NULL;
}

FCSResult mmm2d_tune_near(mmm2d_data_struct *d) {
  const char* fnc_name = "mmm2d_tune_near";
  fcs_int P, n, i, p;
  fcs_float uxrho2m2max, uxrhomax2, T, pref, err, exponent, L, sum;
   
  //printf("node %d, mmm2d_tune_near: max_near: %f, >? box_l[1]/2: %f\n", d->comm.rank, d->max_near, d->box_l[1]/2);
  
  /* yes, it's y only... */
  if (d->max_near > d->box_l[1]/2.)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Layer height too large for MMM2D near formula, increase nodes or layers_per_node.");
  if (d->min_far < 0)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Layer height too small for MMM2D far formula, decrease nodes, layers_per_node or skin.");
  if (d->ux*d->box_l[1] >= 3./M_SQRT2 )
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "box_l[1]/box_l[0] too large for MMM2D near formula, please exchange x and y");

  /* error is split into three parts:
     one part for bessel, one for complex
     and one for polygamma cutoff */
  d->part_error = d->maxPWerror/3.;
  
  /* Bessel sum, determine cutoff */
  P = 2;
  exponent = M_PI*d->ux*d->box_l[1];
  T  = exp(exponent)/exponent;
  pref = 8.*d->ux*mmm_dmax(MMM_COMMON_C_2PI*d->ux, 1.);
  do {
    L = M_PI*d->ux*(P - 1.);
    sum = 0.;
    for (p = 1; p <= P; p++)
      sum += p*exp(-exponent*p);
    err = pref*mmm_K1(d->box_l[1]*L)*(T*((L + d->uy)/M_PI*d->box_l[0] - 1) + sum);
    P++;
  }
  while (err > d->part_error && (P - 1) < MMM2D_MAXIMAL_B_CUT);
  P--;
  if (P == MMM2D_MAXIMAL_B_CUT)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Could not find reasonable Bessel cutoff. Please decrease n_layers or the error bound.");
  
  mmm_realloc_intlist(&(d->besselCutoff), d->besselCutoff.n = P);
  for (p = 1; p < P; p++)
    d->besselCutoff.e[p-1] = (fcs_int)floor(((fcs_float)P)/(2*p)) + 1;
  
  //printf("bessel cutoff %d %g\n", P, err);

  /* complex sum, determine cutoffs (dist dependent) */
  T = log(d->part_error/(16*M_SQRT2)*d->box_l[0]*d->box_l[1]);
  // for 0, the sum is exactly zero, so do not calculate anything
  d->complexCutoff[0] = 0;
  for (i = 1; i <= MMM2D_COMPLEX_STEP; i++)
    d->complexCutoff[i] = (fcs_int)ceil(T/log(i/MMM2D_COMPLEX_FAC));
  mmm2d_prepareBernoulliNumbers(d);

  /* polygamma, determine order */
  n = 1;
  uxrhomax2 = d->ux*d->box_l[1]; //SQR(d->ux*d->box_l[1])/2;
  uxrhomax2 = (uxrhomax2*uxrhomax2)/2.;
  uxrho2m2max = 1.0;
  do {
    mmm_create_mod_psi_up_to(d->polTaylor, n+1);

    err = 2*n*fabs(mmm_mod_psi_even(d->polTaylor, n, 0.5))*uxrho2m2max;
    uxrho2m2max *= uxrhomax2;
    n++;
  }
  while (err > 0.1*d->part_error && n < MMM2D_MAXIMAL_POLYGAMMA);
  if (n == MMM2D_MAXIMAL_POLYGAMMA)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Could find not reasonable Polygamma cutoff. Consider exchanging x and y");
  //printf("polygamma cutoff %d %g\n", n, err);
  
  return NULL;
}

FCSResult mmm2d_tune_far(mmm2d_data_struct *d) {
  const char* fnc_name = "mmm2d_tune_far";
  fcs_float err;
  fcs_float min_inv_boxl = mmm_dmin(d->ux, d->uy);
  d->far_cut = min_inv_boxl;
  do {
    err = exp(-MMM_COMMON_C_2PI*d->far_cut*d->min_far)/d->min_far*
      (MMM_COMMON_C_2PI*d->far_cut + 2.*(d->ux + d->uy) + 1./d->min_far);
      //printf("far tuning: %f -> %f, %f\n", d->far_cut, err, d->min_far);
    d->far_cut += min_inv_boxl;
  }
  while (err > d->maxPWerror && d->far_cut*d->layer_h < MMM2D_MAXIMAL_FAR_CUT);
  if (d->far_cut*d->layer_h >= MMM2D_MAXIMAL_FAR_CUT)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Far cutoff too large, decrease the error bound.");
  d->far_cut -= min_inv_boxl;
  d->far_cut2 = d->far_cut*d->far_cut;
  return NULL;
}

static void mmm2d_prepareBernoulliNumbers(mmm2d_data_struct *d)
{
  fcs_int l, bon_order=d->complexCutoff[MMM2D_COMPLEX_STEP];
  /* BernoulliB[2 n]/(2 n)!(2 Pi)^(2n) up to order 33 */
  static fcs_float bon_table[34] = {
     1.0000000000000000000, 3.2898681336964528729,
    -2.1646464674222763830, 2.0346861239688982794,
    -2.0081547123958886788, 2.0019891502556361707,
    -2.0004921731066160966, 2.0001224962701174097,
    -2.0000305645188173037, 2.0000076345865299997,
    -2.0000019079240677456, 2.0000004769010054555,
    -2.0000001192163781025, 2.0000000298031096567,
    -2.0000000074506680496, 2.0000000018626548648,
    -2.0000000004656623667, 2.0000000001164154418,
    -2.0000000000291038438, 2.0000000000072759591,
    -2.0000000000018189896, 2.0000000000004547474,
    -2.0000000000001136868, 2.0000000000000284217,
    -2.0000000000000071054, 2.0000000000000017764,
    -2.0000000000000004441, 2.0000000000000001110,
    -2.0000000000000000278, 2.0000000000000000069,
    -2.0000000000000000017, 2.0000000000000000004,
    -2.0000000000000000001, 2.0000000000000000000
  };

  if (bon_order < 2)
    bon_order = 2;

  mmm_realloc_doublelist(&(d->bon), d->bon.n = bon_order);

  /* the ux is multiplied in to bessel, complex and psi at once, not here,
     and we use uy*(z + iy), so the uy is also treated below */
  for(l = 1; (l <= bon_order) && (l < 34); l++)
    d->bon.e[l-1] = 2*d->uy*bon_table[l];

  for (; l <= bon_order; l++) {
    if (l & 1)
      d->bon.e[l-1] =  4.0*d->uy;
    else
      d->bon.e[l-1] = -4.0*d->uy;
  }
}

/* compute the net charge of the system */
FCSResult mmm2d_check_system_charges(mmm2d_data_struct *d,
                fcs_int num_particles, fcs_float *charges) {
   /*
  fcs_int i;
  *totcharge=0.0;
  *anycharge=0;
  for (i = 0; i < num_particles; i++) {
    if (!fcs_float_is_zero(charges[i])) {
      *totcharge += charges[i];
      *anycharge=1;
    }
  }
  return NULL;
  */
  fcs_int i;
  fcs_float local_charge=0., total_charge=0.;
  for (i = 0; i < num_particles; i++) {
    if (!fcs_float_is_zero(charges[i])) {
      local_charge += charges[i];
    }
  }
  fprintf(stderr,"check charges: rank %d\n",d->comm.rank);
  fprintf(stderr, "local_charges %d, net charge %e\n", num_particles, local_charge);
  MPI_Allreduce(&local_charge, &total_charge, 1, FCS_MPI_FLOAT, MPI_SUM, d->comm.mpicomm);
  d->total_charge=total_charge;
  fprintf(stderr, "total charge %e\n", d->total_charge);
  return NULL;
}
