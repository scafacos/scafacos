/*
  Copyright (C) 2011,2012 Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.
  
  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/* Implementation of the Ewald summation

   This implementation uses the same notation and variable names as
   Deserno, Holm. "How to mesh up Ewald sums". J. Chem. Phys. 109(18), 1998.
*/

#include "fcs_ewald.h"
#include "FCSResult.h"
#include "FCSInterface.h"
#include "FCSCommon.h"
#include "common/near/near.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FCS_EWALD_USE_ERFC_APPROXIMATION 0

/* Mathematical constants, from gcc's math.h */
#ifndef M_PI
#define M_E             2.7182818284590452353602874713526625L  /* e */
#define M_LOG2E         1.4426950408889634073599246810018921L  /* log_2 e */
#define M_LOG10E        0.4342944819032518276511289189166051L  /* log_10 e */
#define M_LN2           0.6931471805599453094172321214581766L  /* log_e 2 */
#define M_LN10          2.3025850929940456840179914546843642L  /* log_e 10 */
#define M_PI            3.1415926535897932384626433832795029L  /* pi */
#define M_PI_2          1.5707963267948966192313216916397514L  /* pi/2 */
#define M_PI_4          0.7853981633974483096156608458198757L  /* pi/4 */
#define M_1_PI          0.3183098861837906715377675267450287L  /* 1/pi */
#define M_2_PI          0.6366197723675813430755350534900574L  /* 2/pi */
#define M_2_SQRTPI      1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#define M_SQRT2         1.4142135623730950488016887242096981L  /* sqrt(2) */
#define M_SQRT1_2       0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif

#define EWALD_DEFAULT_TOLERANCE_FIELD 1.0e-6

/* precision of alpha */
#ifdef FCS_FLOAT_IS_FLOAT
#define ALPHA_OPT_PREC 1.0e-5
#else
#define ALPHA_OPT_PREC 1.0e-10
#endif

/* Debug macros */
#ifdef FCS_ENABLE_DEBUG
#define ADDITIONAL_CHECKS
#define FCS_TRACE(cmd) if (on_root()) cmd
#define FCS_DEBUG(cmd) if (on_root()) cmd
#define FCS_DEBUG_ALL(cmd) cmd
#else
#define FCS_TRACE(cmd)
#define FCS_DEBUG(cmd) 
#define FCS_DEBUG_ALL(cmd)
#endif

#ifdef FCS_ENABLE_INFO
#define FCS_INFO(cmd) if (on_root()) cmd
#else
#define FCS_INFO(cmd)
#endif

#define MAXKMAX_DEFAULT 100

static inline int on_root()
{
  static int rank = -1;
  if (rank < 0)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == 0;
}

inline static fcs_int 
linindex(fcs_int ix, fcs_int iy, fcs_int iz, 
	 fcs_int kmax) {
  fcs_int km = kmax+1;
  return ix*km*km + iy*km + iz;
}

inline static fcs_float
SQR(fcs_float x) { return x*x; }

typedef struct {
  /* MPI communicator */
  MPI_Comm comm;
  int comm_size;
  int comm_rank;

  /* Cartesian communicator (for near field solver) */
  MPI_Comm comm_cart;

  /* System size in x,y,z */
  fcs_float box_l[3];

  /** Required accuracy. */
  fcs_float tolerance_field;

  /** Splitting parameter */
  fcs_float alpha;

  /** Near-field cutoff */
  fcs_float r_cut;

  /** Kspace cutoff */
  fcs_int kmax;

  /** maximal Kspace cutoff used by tuning */
  fcs_int maxkmax;

  /* influence function */
  fcs_float* G;

  /* Whether or not it is necessary to retune the method before running it. */
  int needs_retune;
  /** Whether or not rcut is to be automatically tuned. */
  int tune_r_cut;
  /** Whether or not kmax is to be automatically tuned. */
  int tune_kmax;
  /** Whether or not alpha is to be automatically tuned. */
  int tune_alpha;

  /* The components of the fields and potentials */
  fcs_float* far_fields;
  fcs_float* near_fields;
  fcs_float* far_potentials;
  fcs_float* near_potentials;
} ewald_data_struct;

/* method to check whether a given FCS is EWALD */
static FCSResult fcs_ewald_check(FCS handle, const char* fnc_name) {
  if (handle == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_METHOD_EWALD)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"EWALD\".");

  return NULL;
}

FCSResult fcs_ewald_init(FCS handle) {
  const char* fnc_name = "ewald_init";
  FCSResult result;

  result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  /* initialize ewald struct */
  ewald_data_struct *d;
  if (handle->method_context == NULL) {
    d = malloc(sizeof(ewald_data_struct));
    handle->method_context = d;
  } else {
    d = (ewald_data_struct*)handle->method_context;
  }

  /* MPI communicator */
  d->comm = handle->communicator;
  MPI_Comm_size(d->comm, &d->comm_size);
  MPI_Comm_rank(d->comm, &d->comm_rank);

  /* set up cartesian communicator */
  int node_grid[3] = {0, 0, 0};
  
  /* compute node grid */
  MPI_Dims_create(d->comm_size, 3, node_grid);

  /* create cartesian communicator */
  int periodicity[3] = {1, 1, 1};
  MPI_Cart_create(d->comm, 3, node_grid, periodicity, 
		  1, &d->comm_cart);
  
  d->box_l[0] = 0.0;
  d->box_l[1] = 0.0;
  d->box_l[2] = 0.0;
  d->tolerance_field = EWALD_DEFAULT_TOLERANCE_FIELD;
  d->alpha = 0.0;
  d->r_cut = 0.0;
  d->kmax = 0;
  d->maxkmax = MAXKMAX_DEFAULT;
  /* d->alpha = 1.0; */
  /* d->r_cut = 3.0; */
  /* d->kmax = 40; */

  d->needs_retune = 1;
  d->tune_alpha = 1;
  d->tune_r_cut = 1;
  d->tune_kmax = 1;
  d->G = NULL;

  d->far_fields = NULL;
  d->near_fields = NULL;
  d->far_potentials = NULL;
  d->near_potentials = NULL;

  return NULL;
}

FCSResult fcs_ewald_set_maxkmax(FCS handle, fcs_int maxkmax) {
  const char *fnc_name = "fcs_ewald_set_maxkmax";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->maxkmax = maxkmax;

  return NULL;
}

FCSResult fcs_ewald_set_maxkmax_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_maxkmax_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->maxkmax = MAXKMAX_DEFAULT;
  
  return NULL;
}

FCSResult fcs_ewald_get_maxkmax(FCS handle, fcs_int *maxkmax) {
  const char *fnc_name = "fcs_ewald_get_kmax";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  *maxkmax = d->maxkmax;

  return NULL;
}

FCSResult fcs_ewald_set_kmax(FCS handle, fcs_int kmax) {
  const char *fnc_name = "fcs_ewald_set_kmax";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  if (kmax != d->kmax)
    d->needs_retune = 1;
  d->kmax = kmax;
  d->tune_kmax = 0;

  return NULL;
}

FCSResult fcs_ewald_set_kmax_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_kmax_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->needs_retune = 1;
  d->tune_kmax = 1;

  return NULL;
}

FCSResult fcs_ewald_get_kmax(FCS handle, fcs_int *kmax) {
  const char *fnc_name = "fcs_ewald_get_kmax";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  *kmax = d->kmax;

  return NULL;
}

FCSResult fcs_ewald_set_r_cut(FCS handle, fcs_float r_cut) {
  const char *fnc_name = "fcs_ewald_set_r_cut";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  if (!fcs_float_is_equal(r_cut, d->r_cut))
    d->needs_retune = 1;
  d->r_cut = r_cut;
  d->tune_r_cut = 0;

  return NULL;
}

FCSResult fcs_ewald_set_r_cut_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_r_cut_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->needs_retune = 1;
  d->tune_r_cut = 1;

  return NULL;
}

FCSResult fcs_ewald_get_r_cut(FCS handle, fcs_float *r_cut) {
  const char *fnc_name = "fcs_ewald_get_r_cut";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  *r_cut = d->r_cut;

  return NULL;
}

FCSResult fcs_ewald_set_alpha(FCS handle, fcs_float alpha) {
  const char *fnc_name = "fcs_ewald_set_alpha";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  if (!fcs_float_is_equal(alpha, d->alpha))
    d->needs_retune = 1;
  d->alpha = alpha;
  d->tune_alpha = 0;

  return NULL;
}

FCSResult fcs_ewald_set_alpha_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_alpha_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  d->needs_retune = 1;
  d->tune_alpha = 1;
  return NULL;
}

FCSResult fcs_ewald_get_alpha(FCS handle, fcs_float *alpha) {
  const char *fnc_name = "fcs_ewald_get_alpha";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  *alpha = d->alpha;
  return NULL;
}

FCSResult fcs_ewald_set_tolerance_field(FCS handle, fcs_float tolerance_field) {
  const char *fnc_name = "fcs_ewald_set_tolerance_field";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  if (!fcs_float_is_equal(tolerance_field, d->tolerance_field))
    d->needs_retune = 1;
  d->tolerance_field = tolerance_field;
  return NULL;
}

FCSResult fcs_ewald_set_tolerance_field_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_tolerance_field_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->needs_retune = 1;
  d->tolerance_field = EWALD_DEFAULT_TOLERANCE_FIELD;
  return NULL;
}  

FCSResult fcs_ewald_get_tolerance_field(FCS handle, fcs_float* tolerance_field) {
  const char *fnc_name = "fcs_ewald_get_tolerance_field";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  *tolerance_field = d->tolerance_field;
  return NULL;
}

FCSResult 
fcs_ewald_get_components(FCS handle, 
			 fcs_float** far_fields, 
			 fcs_float** near_fields,
			 fcs_float** far_potentials, 
			 fcs_float** near_potentials) {
  const char* fnc_name =  "fcs_P3M_get_components";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  if (far_fields != NULL) {
    if (d->far_fields == NULL)
      return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"far_fields wanted but not computed."); 
    else
      *far_fields = d->far_fields;
  }

  if (near_fields != NULL) {
    if (d->near_fields == NULL)
      return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"near_fields wanted but not computed."); 
    else
      *near_fields = d->near_fields;
  }

  if (far_potentials != NULL) {
    if (d->far_potentials == NULL)
      return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"far_potentials wanted but not computed."); 
    else
      *far_potentials = d->far_potentials;
  }

  if (near_potentials != NULL) {
    if (d->near_potentials == NULL)
      return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"near_potentials wanted but not computed."); 
    else
      *near_potentials = d->near_potentials;
  }
  return NULL;
}

/***************************************************
 **** TUNING
 ***************************************************/



/** Computes the real space contribution to the rms error in the
   force (as described by Kolafa and Perram) as well as its
   derivative.
   \param N the number of charged particles in the system.
   \param sum_q2 the sum of square of charges in the system
   \param box_l fcs_float[3] system size
   \param r_cut the real-space cutoff
   \param alpha the Ewald splitting parameter
   \param error (out) the real space error
   \param derivative (out) the derivative of the real space error
*/
static inline void
ewald_real_space_error(fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
		       fcs_float r_cut, fcs_float alpha,
		       fcs_float *error, fcs_float *derivative) {
  const fcs_float V = box_l[0]*box_l[1]*box_l[2];
  *error = 2.0*sum_q2 / sqrt(N*r_cut*V) * exp(-SQR(alpha*r_cut));
  *derivative = -2.0 * alpha * SQR(r_cut) * *error;
}

/** Computes the reciprocal space contribution to the rms error in the
   force (as described by Kolafa and Perram) as well as its derivative.
   \param N the number of charged particles in the system.
   \param sum_q2 the sum of square of charges in the system
   \param box_l fcs_float[3] system size
   \param kmax the maximal k vector
   \param alpha the Ewald splitting parameter
   \param error (out) the recipocal space error
   \param derivative (out) the derivative of the recipocal space error
*/
static inline void
ewald_k_space_error(fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
		    fcs_int kmax, fcs_float alpha,
		    fcs_float *error, fcs_float *derivative) {
  fcs_float Lmax = box_l[0];
  if (box_l[1] > Lmax) Lmax = box_l[1];
  if (box_l[2] > Lmax) Lmax = box_l[2];
  const fcs_float K = 2.0*M_PI*kmax/Lmax;
  const fcs_float V = box_l[0]*box_l[1]*box_l[2];
  const fcs_float fak1 = 2.0*sum_q2 / sqrt(N*M_PI*K*V) * exp(-K*K/(4.0*SQR(alpha)));
  *error = fak1 * alpha;
  *derivative = fak1 * (2.0*SQR(alpha)+SQR(kmax)/(2.0*SQR(alpha)));
}

static inline void
ewald_tune_alpha(fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
                 fcs_float r_cut, fcs_int kmax, 
                 fcs_float *alpha, fcs_float *error, int tune) {
  fcs_float err_r, der_r;
  fcs_float err_k, der_k;
  fcs_float alpha_diff;

  if (tune) {
    /* Newton-Raphson method to find optimal alpha */
    *alpha = 0.1;
    do {
      /* get errors and derivatives */
      ewald_real_space_error(N, sum_q2, box_l, r_cut, *alpha, &err_r, &der_r);
      ewald_k_space_error(N, sum_q2, box_l, kmax, *alpha, &err_k, &der_k);
      
      alpha_diff = (err_r - err_k) / (der_r - der_k);
      *alpha -= alpha_diff;
    } while (fabs(alpha_diff) > ALPHA_OPT_PREC);
  } else {
    ewald_real_space_error(N, sum_q2, box_l, r_cut, *alpha, &err_r, &der_r);
    ewald_k_space_error(N, sum_q2, box_l, kmax, *alpha, &err_k, &der_k);
  }
  *error = sqrt(SQR(err_r) + SQR(err_k));
}


FCSResult fcs_ewald_tune(FCS handle,
			 fcs_int num_particles,
			 fcs_int local_max_particles,
			 fcs_float *positions, 
			 fcs_float *charges) {
  const char *fnc_name = "fcs_ewald_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  FCS_INFO(fprintf(stderr, "fcs_ewald_tune() started...\n"));

  /* Handle periodicity */
  const fcs_int *periodicity = fcs_get_periodicity(handle);
  if (! (periodicity[0] && periodicity[1] && periodicity[2]))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
  			    "EWALD requires periodic boundary conditions.");
  
  /* Handle box size */
  const fcs_float *a = fcs_get_box_a(handle);
  const fcs_float *b = fcs_get_box_b(handle);
  const fcs_float *c = fcs_get_box_c(handle);
  if (!fcs_is_orthogonal(a, b, c))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
  			    "EWALD requires the box to be orthorhombic.");
  
  if (!fcs_uses_principal_axes(a, b, c))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
  			    "EWALD requires the box vectors to be parallel to the principal axes.");
  
  if (!fcs_float_is_equal(a[0], d->box_l[0])) {
    d->box_l[0] = a[0];
    d->needs_retune = 1;
  }
  if (!fcs_float_is_equal(b[1], d->box_l[1])) {
    d->box_l[1] = b[1];
    d->needs_retune = 1;
  }
  if (!fcs_float_is_equal(c[2], d->box_l[2])) {
    d->box_l[2] = c[2];
    d->needs_retune = 1;
  }
  
  /* Check whether retuning is required */
  if (!d->needs_retune) {
    FCS_INFO(fprintf(stderr, "  Retuning is not required.\n"));
    FCS_INFO(fprintf(stderr, "fcs_ewald_tune() finished.\n"));
    return NULL;
  }
  
  FCS_INFO(fprintf(stderr, "  Retuning is required.\n"));
  
  /* Check whether the input parameters are sane */
  if (!d->tune_r_cut) {
    if (d->r_cut < 0.0) {
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "r_cut is negative!");
    }
    
    if (fcs_float_is_zero(d->r_cut)) {
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "r_cut is too small!");
    }
    
    /* check whether cutoff is larger than half a box length */
    if ((d->r_cut > 0.5*d->box_l[0]) ||
  	(d->r_cut > 0.5*d->box_l[1]) ||
  	(d->r_cut > 0.5*d->box_l[2]))
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "r_cut is larger than half a system box length.");
  }

  /* Init tuning */
  /* count the number of charges */
  fcs_int local_num_charges = 0;
  fcs_int num_charges;
  fcs_float local_q2 = 0;
  fcs_float q2;
  for (fcs_int i = 0; i < num_particles; i++) {
    if (!fcs_float_is_zero(charges[i]))
      local_num_charges++;
    local_q2 += charges[i]*charges[i];
  }
  MPI_Reduce(&local_num_charges, &num_charges, 1, FCS_MPI_INT,
  	     MPI_SUM, 0, d->comm);
  MPI_Reduce(&local_q2, &q2, 1, FCS_MPI_FLOAT,
  	     MPI_SUM, 0, d->comm);
  
  if (d->comm_rank == 0) {
    /* Tune r_cut */
    if (d->tune_r_cut) {
      /* compute the average distance between two charges */
      fcs_float avg_dist = pow((d->box_l[0]*d->box_l[1]*d->box_l[2])
  			       / num_charges, 0.33333);
      /* set r_cut to 3 times the average distance between charges */
      d->r_cut = 3.0 * avg_dist;
      /* tune r_cut to half the box length */
      if (0.5*d->box_l[0] < d->r_cut) d->r_cut = 0.5*d->box_l[0];
      if (0.5*d->box_l[1] < d->r_cut) d->r_cut = 0.5*d->box_l[1];
      if (0.5*d->box_l[2] < d->r_cut) d->r_cut = 0.5*d->box_l[2];
      
      FCS_INFO(fprintf(stderr, "    tuned r_cut to %" FCS_LMOD_FLOAT "f\n", d->r_cut));
    }

    fcs_float estimated_accuracy;
    /* Tune kmax and alpha */
    if (d->tune_kmax) {
      estimated_accuracy = 1.0e10;
      d->kmax = 0;
      while (estimated_accuracy > d->tolerance_field) {
  	d->kmax++;
  	if (d->kmax > d->maxkmax) {
  	  static char msg[255];
  	  sprintf(msg, "EWALD cannot achieve required accuracy with given r_cut=%" FCS_LMOD_FLOAT "f with a kmax<=%" FCS_LMOD_INT "d",
  		  d->r_cut, d->maxkmax);
  	  return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, msg);
  	}

	ewald_tune_alpha(num_charges, q2, d->box_l, d->r_cut, d->kmax, 
			 &d->alpha, &estimated_accuracy, d->tune_alpha);
      }
      FCS_INFO(if (d->tune_alpha) fprintf(stderr, "    tuned alpha to %" FCS_LMOD_FLOAT "f\n", d->alpha));
      FCS_INFO(fprintf(stderr, "    tuned kmax to %" FCS_LMOD_INT "d\n", d->kmax));
      FCS_INFO(fprintf(stderr, "    estimated_accuracy=%" FCS_LMOD_FLOAT "e\n", estimated_accuracy));
    } else {
      /* Tune only alpha */
      ewald_tune_alpha(num_charges, q2, d->box_l, d->r_cut, d->kmax, 
		       &d->alpha, &estimated_accuracy, d->tune_alpha);
      FCS_INFO(fprintf(stderr, "    tuned alpha to %" FCS_LMOD_FLOAT "f\n", d->alpha));
      FCS_INFO(fprintf(stderr, "    estimated_accuracy=%" FCS_LMOD_FLOAT "e\n", estimated_accuracy));
      if (estimated_accuracy > d->tolerance_field) {
  	char msg[255];
	if (d->tune_alpha)
	  sprintf(msg, "EWALD cannot achieve required accuracy with given kmax=%" FCS_LMOD_INT "d and r_cut=%" FCS_LMOD_FLOAT "f", d->kmax, d->r_cut);
	else
	  sprintf(msg, "EWALD cannot achieve required accuracy with given kmax=%" FCS_LMOD_INT "d, r_cut=%" FCS_LMOD_FLOAT "f and alpha=%" FCS_LMOD_FLOAT "f", d->kmax, d->r_cut, d->alpha);
  	return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, msg);
      }
    }
  }

  /* broadcast the tuned parameters */
  if (d->tune_r_cut)
    MPI_Bcast(&d->r_cut, 1, FCS_MPI_FLOAT, 0, d->comm);
  if (d->tune_alpha)
    MPI_Bcast(&d->alpha, 1, FCS_MPI_FLOAT, 0, d->comm);
  if (d->tune_kmax)
    MPI_Bcast(&d->kmax, 1, FCS_MPI_INT, 0, d->comm);

  /* Compute influence function */
  /* The influence function is 1/V*g_tilde(k_vec)*gamma_tilde(k_vec) */
  fcs_int size = d->kmax+1;
  size = size*size*size;
  if (d->G == NULL)
    d->G = (fcs_float*)malloc(sizeof(fcs_float)*size);
  else
    d->G = (fcs_float*)realloc(d->G, sizeof(fcs_float)*size);

  for (fcs_int nx=0; nx <= d->kmax; nx++)
    for (fcs_int ny=0; ny <= d->kmax; ny++)
      for (fcs_int nz=0; nz <= d->kmax; nz++) {
	// system length vector L
	const fcs_float Lx = d->box_l[0];
	const fcs_float Ly = d->box_l[1];
	const fcs_float Lz = d->box_l[2];
	// volume V
	const fcs_float V = Lx*Ly*Lz;
	// reciprocal vector k
	const fcs_float kx = 2.0*M_PI*nx / Lx;
	const fcs_float ky = 2.0*M_PI*ny / Ly;
	const fcs_float kz = 2.0*M_PI*nz / Lz;
	// k^2
	const fcs_float ksqr = kx*kx + ky*ky + kz*kz;
	// alpha
	const fcs_float alpha = d->alpha;

  	if ((nx==0 && ny==0 && nz==0) || 
	    nx*nx + ny*ny + nz*nz > d->kmax*d->kmax)
	  // set influence function to 0 for k=0 and outside the sphere
  	  d->G[linindex(nx, ny, nz, d->kmax)] = 0.0;
  	else 
	  // compute the influence function inside the sphere
	  /* Compare to Deserno, Holm (1998) eqs. (5) and (15) */
	  d->G[linindex(nx, ny, nz, d->kmax)]
	    = 1.0/V * 4.0*M_PI/ksqr * exp(-ksqr/(4*SQR(alpha)));
      }

  d->needs_retune = 0;

  FCS_INFO(fprintf(stderr, "fcs_ewald_tune() finished.\n"));
  return NULL;
}

/***************************************************
 **** K-SPACE CONTRIBUTION
 ***************************************************/
void ewald_compute_kspace(ewald_data_struct* d, 
			  fcs_int num_particles,
			  fcs_float *positions, 
			  fcs_float *charges,
			  fcs_float *fields,
			  fcs_float *potentials) {

  FCS_INFO(fprintf(stderr, "ewald_compute_kspace started...\n"));

  /* DISTRIBUTE ALL PARTICLE DATA TO ALL NODES */
  /* Gather all particle numbers */
  int node_num_particles = num_particles;
  int node_particles[d->comm_size]; 
  int node_particles3[d->comm_size]; 
  int displs[d->comm_size];
  int displs3[d->comm_size];
  fcs_int total_particles;

  /* printf("%d: num_particles=%d\n", d->comm_rank, num_particles); */
  MPI_Allgather(&node_num_particles, 1, MPI_INT, node_particles, 1, MPI_INT, d->comm);

  /* compute displacements for MPI_Gatherv */
  total_particles = node_particles[0];
  node_particles3[0] = node_particles[0]*3;
  displs[0] = 0;
  displs3[0] = 0;
  for (fcs_int i=1; i < d->comm_size; i++) {
    total_particles += node_particles[i];
    node_particles3[i] = node_particles[i]*3;
    displs[i] = displs[i-1] + node_particles[i-1];
    displs3[i] = displs3[i-1] + node_particles3[i-1];
  }

  fcs_float all_positions[total_particles*3];
  fcs_float all_charges[total_particles];
  fcs_float node_fields[total_particles*3];
  fcs_float node_potentials[total_particles];

  /* gather all particle data at all nodes */
  MPI_Allgatherv(positions, num_particles*3, FCS_MPI_FLOAT, all_positions, node_particles3, displs3, FCS_MPI_FLOAT, d->comm);
  MPI_Allgatherv(charges, num_particles, FCS_MPI_FLOAT, all_charges, node_particles, displs, FCS_MPI_FLOAT, d->comm);

  /* for (fcs_int i=0; i < total_particles; i++) { */
  /*   printf("%d: all_positions[%d]={%lf, %lf, %lf}\n", d->comm_rank, i, all_positions[i*3], all_positions[3*i+1], all_positions[3*i+2]); */
  /*   printf("%d: all_charges[%d]=%lf\n", d->comm_rank, i, all_charges[i]); */
  /* } */

  /* INIT ALGORITHM */

  /* init fields and potentials */
  if (fields != NULL) {
    for (fcs_int i=0; i < total_particles; i++) {
      node_fields[3*i] = 0.0;
      node_fields[3*i+1] = 0.0;
      node_fields[3*i+2] = 0.0;
    }
  }
  if (potentials != NULL)
    for (fcs_int i=0; i < total_particles; i++)
      node_potentials[i] = 0.0;
  
  /* COMPUTE FAR FIELDS */

  /* evenly distribute the k-vectors onto all tasks */
  fcs_int num_k_per_dir = 2*d->kmax+1;
  fcs_int num_k = num_k_per_dir * num_k_per_dir * num_k_per_dir;
  for (int k_ind = d->comm_rank; k_ind < num_k; k_ind += d->comm_size) {
    /* compute fields and potentials */
    fcs_int nx = 
      k_ind % num_k_per_dir - d->kmax;
    fcs_int ny = 
      k_ind % (num_k_per_dir*num_k_per_dir) / num_k_per_dir - d->kmax;
    fcs_int nz = 
      k_ind / (num_k_per_dir*num_k_per_dir) - d->kmax;
    if (nx*nx + ny*ny + nz*nz <= d->kmax*d->kmax) {
      // system length vector L_vec
      const fcs_float Lx = d->box_l[0];
      const fcs_float Ly = d->box_l[1];
      const fcs_float Lz = d->box_l[2];
      // reciprocal vector k_vec
      const fcs_float kx = 2.0*M_PI*nx / Lx;
      const fcs_float ky = 2.0*M_PI*ny / Ly;
      const fcs_float kz = 2.0*M_PI*nz / Lz;
      /* reciprocal charge density rhohat */
      fcs_float rhohat_re = 0.0;
      fcs_float rhohat_im = 0.0;
      
      /* compute Deserno, Holm (1998) eq. (8) */
      for (fcs_int i=0; i < total_particles; i++) {
	// charge q
	const fcs_float q = all_charges[i];
	if (!fcs_float_is_zero(q)) {
	  // particle position r_vec
	  const fcs_float rx = all_positions[3*i];
	  const fcs_float ry = all_positions[3*i+1];
	  const fcs_float rz = all_positions[3*i+2];
	  // compute k_vec*r_vec
	  fcs_float kr = kx*rx + ky*ry + kz*rz;
	  // rhohat = qi * exp(-i*k_vec*r_vec)
	  rhohat_re += q * cos(kr);
	  rhohat_im += q * -sin(kr);
	}
      }
      
      /* FCS_DEBUG(fprintf(stderr, "  n_vec= (%d, %d, %d) rhohat_re=%e rhohat_im=%e\n", \ */
      /* 		    nx, ny, nz, rhohat_re, rhohat_im)); */
      
      /* fetch influence function */
      fcs_float g = d->G[linindex(abs(nx), abs(ny), abs(nz), d->kmax)];
      for (fcs_int i=0; i < total_particles; i++) {
	// particle position r_vec
	const fcs_float rx = all_positions[3*i];
	const fcs_float ry = all_positions[3*i+1];
	const fcs_float rz = all_positions[3*i+2];
	// compute k_vec*r_vec
	fcs_float kr = kx*rx + ky*ry + kz*rz;
	
	if (fields != NULL) {
	  /* compute field at position of particle i
	     compare to Deserno, Holm (1998) eq. (15) */
	  fcs_float fak1 = g * (rhohat_re*sin(kr) + rhohat_im*cos(kr));
	  node_fields[3*i] += kx * fak1;
	  node_fields[3*i+1] += ky * fak1;
	  node_fields[3*i+2] += kz * fak1;
	}
	if (potentials != NULL) {
	  /* compute potential at position of particle i 
	     compare to Deserno, Holm (1998) eq. (9) */
	  node_potentials[i] += g * (rhohat_re*cos(kr) - rhohat_im*sin(kr));
	}
      }
    }
  }
  
  /* printf("%d: node_fields[0]=%lf\n", d->comm_rank, node_fields[0]); */

  /* REDISTRIBUTE COMPUTED FAR FIELDS AND POTENTIALS  */
  if (fields != NULL) {
    fcs_float all_fields[total_particles*3];
    for (fcs_int pid=0; pid < total_particles; pid++) {
      all_fields[3*pid] = 0.0;
      all_fields[3*pid+1] = 0.0;
      all_fields[3*pid+2] = 0.0;
    }

    /* Combine all fields on master */
    MPI_Reduce(node_fields, all_fields, total_particles*3, 
	       FCS_MPI_FLOAT, MPI_SUM, 0, d->comm);
    /* if (d->comm_rank == 0) */
    /*   printf("%d: all_fields[0]=%lf\n", d->comm_rank, all_fields[0]); */
    /* Scatter the fields to the task that holds the particle */
    MPI_Scatterv(all_fields, node_particles3, displs3, FCS_MPI_FLOAT,
		 fields, num_particles*3, FCS_MPI_FLOAT, 0, d->comm);
  }

  if (potentials != NULL) {
    fcs_float all_potentials[total_particles];
    for (fcs_int pid=0; pid < total_particles; pid++)
      all_potentials[pid] = 0.0;
    /* Combine all potentials on master */
    MPI_Reduce(node_potentials, all_potentials, total_particles, 
	       FCS_MPI_FLOAT, MPI_SUM, 0, d->comm);
    /* Scatter the potentials to the task that holds the particle */
    MPI_Scatterv(all_potentials, node_particles, displs, FCS_MPI_FLOAT,
		 potentials, num_particles, FCS_MPI_FLOAT, 0, d->comm);

    /* subtract self potential */
    FCS_INFO(fprintf(stderr, "  subtracting self potential...\n"));
    for (fcs_int i=0; i < num_particles; i++) {
      potentials[i] -= charges[i] * M_2_SQRTPI * d->alpha;
    }
  }

  /* now each task should have its far field components */
  FCS_INFO(fprintf(stderr, "ewald_compute_kspace finished.\n"));
}

/***************************************************
 **** REAL SPACE CONTRIBUTION
 ***************************************************/
/* callback function for near field computations */
static inline void 
ewald_compute_near(const void *param, fcs_float dist, fcs_float *f, fcs_float *p)
{
  fcs_float alpha = *((fcs_float *) param);
  fcs_float adist = alpha * dist;

#if FCS_EWALD_USE_ERFC_APPROXIMATION

  /* approximate \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying a formula from:
     Abramowitz/Stegun: Handbook of Mathematical Functions, Dover
     (9. ed.), chapter 7 */
  fcs_float t = 1.0 / (1.0 + 0.3275911 * adist);
  fcs_float erfc_part_ri = exp(-adist*dist) * 
    (t * (0.254829592 + 
	  t * (-0.284496736 + 
	       t * (1.421413741 + 
		    t * (-1.453152027 + 
			 t * 1.061405429))))) 
    / dist;

  *p = erfc_part_ri;
  *f = -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498) 
    / dist;
  
#else

  fcs_float erfc_part_ri = erfc(adist) / dist;
  *p = erfc_part_ri;
  *f = -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498*exp(-adist*adist)) 
    / dist;

#endif
}

/* callback function for performing a whole loop of near field computations (using ewald_compute_near) */
FCS_NEAR_LOOP_FP(ewald_compute_near_loop, ewald_compute_near)

void ewald_compute_rspace(ewald_data_struct* d, 
			  fcs_int num_particles,
			  fcs_int max_num_particles,
			  fcs_float *positions, 
			  fcs_float *charges,
			  fcs_float *fields,
			  fcs_float *potentials) {
  FCS_INFO(fprintf(stderr, "ewald_compute_rspace started...\n"));

  fcs_int local_num_particles;
  fcs_int local_num_real_particles;
  fcs_int local_num_ghost_particles;
  fcs_float *local_positions, *local_ghost_positions;
  fcs_float *local_charges, *local_ghost_charges;
  fcs_gridsort_index_t *local_indices, *local_ghost_indices;
  fcs_gridsort_t gridsort;
  const fcs_float box_base[3] = {0.0, 0.0, 0.0 };
  const fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
  const fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
  const fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };

  /* DOMAIN DECOMPOSE */
  fcs_gridsort_create(&gridsort);
  fcs_gridsort_set_system(&gridsort, box_base, box_a, box_b, box_c, NULL);
  fcs_gridsort_set_particles(&gridsort, num_particles, max_num_particles,
			     positions, charges);

  FCS_INFO(fprintf(stderr, "  calling fcs_gridsort_sort_forward()...\n"));
  fcs_gridsort_sort_forward(&gridsort, d->r_cut, d->comm_cart);
  FCS_INFO(fprintf(stderr, "  returning from fcs_gridsort_sort_forward().\n"));

  fcs_gridsort_separate_ghosts(&gridsort, &local_num_real_particles, 
			       &local_num_ghost_particles);
  fcs_gridsort_get_sorted_particles(&gridsort, &local_num_particles, NULL,
				    NULL, NULL, NULL);
  fcs_gridsort_get_real_particles(&gridsort, &local_num_real_particles,
				  &local_positions, &local_charges, 
				  &local_indices);
  fcs_gridsort_get_ghost_particles(&gridsort, 
				   &local_num_ghost_particles, 
				   &local_ghost_positions, 
				   &local_ghost_charges, 
				   &local_ghost_indices);

  FCS_DEBUG_ALL(MPI_Barrier(d->comm_cart));
  FCS_DEBUG(fprintf(stderr,						\
		   "    %d: local_num_particles=%" FCS_LMOD_INT "d local_num_real_particles=%" FCS_LMOD_INT "d local_num_ghost_particles=%" FCS_LMOD_INT "d\n", \
		   d->comm_rank, local_num_particles,			\
		   local_num_real_particles, local_num_ghost_particles));

  /* allocate local fields and potentials */
  fcs_float *local_fields = NULL;
  fcs_float *local_potentials = NULL;
  if (fields != NULL) {
    local_fields = malloc(sizeof(fcs_float)*3*local_num_real_particles);
    for (fcs_int pid=0; pid < local_num_real_particles; pid++) {
      local_fields[3*pid] = 0.0;
      local_fields[3*pid+1] = 0.0;
      local_fields[3*pid+2] = 0.0;
    }
  }
  if (potentials != NULL) {
    local_potentials = malloc(sizeof(fcs_float)*local_num_real_particles);
    for (fcs_int pid=0; pid < local_num_real_particles; pid++)
      local_potentials[pid] = 0.0;
  }

  /* COMPUTE NEAR FIELD */
  fcs_near_t near;
  fcs_near_create(&near);
  fcs_near_set_loop(&near, ewald_compute_near_loop);
  fcs_near_set_system(&near, box_base, box_a, box_b, box_c, NULL);

  fcs_near_set_particles(&near, local_num_real_particles, local_num_real_particles,
			 local_positions, local_charges, 
			 local_indices,
			 (fields != NULL)?local_fields:NULL, 
			 (potentials != NULL)?local_potentials:NULL);
  
  fcs_near_set_ghosts(&near, local_num_ghost_particles,
		      local_ghost_positions, local_ghost_charges, 
		      local_ghost_indices);

  FCS_INFO(fprintf(stderr, "  calling fcs_near_compute()...\n"));
  fcs_near_compute(&near, d->r_cut, &d->alpha, d->comm_cart);
  FCS_INFO(fprintf(stderr, "  returning from fcs_near_compute().\n"));
  fcs_near_destroy(&near);

  /* RECOMPOSE FIELDS AND POTENTIALS */
  FCS_INFO(fprintf(stderr, "  calling fcs_gridsort_sort_backward()...\n"));
  fcs_gridsort_sort_backward(&gridsort,
			     local_fields, local_potentials,
			     fields, potentials, 1,
			     d->comm_cart);
  FCS_INFO(fprintf(stderr, "  returning from fcs_gridsort_sort_backward().\n"));
  
  fcs_gridsort_free(&gridsort);
  fcs_gridsort_destroy(&gridsort);

  FCS_INFO(fprintf(stderr, "ewald_compute_rspace finished.\n"));
}

FCSResult fcs_ewald_run(FCS handle,
			fcs_int num_particles,
			fcs_int local_max_particles,
			fcs_float *positions, 
			fcs_float *charges,
			fcs_float *fields,
			fcs_float *potentials) {
  const char *fnc_name = "fcs_ewald_run";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  /* First run tune */
  fcs_ewald_tune(handle, num_particles, local_max_particles, positions, charges);

  FCS_INFO(fprintf(stderr, "fcs_ewald_run() started...\n"));
  FCS_INFO(fprintf(stderr,						\
		   "    system parameters: box_l=(%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f)\n",	\
		   d->box_l[0], d->box_l[1], d->box_l[2]));
  FCS_INFO(fprintf(stderr,						\
		   "    ewald params: r_cut=%" FCS_LMOD_FLOAT "f, alpha=%" FCS_LMOD_FLOAT "f kmax=%" FCS_LMOD_INT "d\n",		\
		   d->r_cut, d->alpha, d->kmax));

  if (fields != NULL) {
    if (d->far_fields == NULL)
      d->far_fields = malloc(num_particles*sizeof(fcs_float)*3);
    else
      d->far_fields = realloc(d->far_fields, 
			      num_particles*sizeof(fcs_float)*3);
  }
  
  if (potentials != NULL) {
    if (d->far_potentials == NULL)
      d->far_potentials = malloc(num_particles*sizeof(fcs_float));
      else
  	d->far_potentials = realloc(d->far_potentials,
				    num_particles*sizeof(fcs_float));
  }

  /* Compute far field component */
  ewald_compute_kspace(d, num_particles, positions, charges,
  		       fields==NULL ? NULL: d->far_fields,
		       potentials==NULL ? NULL : d->far_potentials
		       );


  if (fields != NULL) {
    if (d->near_fields == NULL)
      d->near_fields = malloc(num_particles*sizeof(fcs_float)*3);
    else
      d->near_fields = realloc(d->near_fields,
			       num_particles*sizeof(fcs_float)*3);
  }
  if (potentials != NULL) {
    if (d->near_potentials == NULL)
      d->near_potentials = malloc(num_particles*sizeof(fcs_float));
    else
      d->near_potentials = realloc(d->near_potentials,
				   num_particles*sizeof(fcs_float));
  }

  /* Compute near field component */
  ewald_compute_rspace(d, num_particles, local_max_particles, positions, charges,
		       fields==NULL ? NULL : d->near_fields, 
		       potentials==NULL ? NULL : d->near_potentials);

  /* Add up components */
  if (fields != NULL) 
    for (fcs_int pid=0; pid < num_particles; pid++) {
      fields[3*pid] = d->far_fields[3*pid] + d->near_fields[3*pid];
      fields[3*pid+1] = d->far_fields[3*pid+1] + d->near_fields[3*pid+1];
      fields[3*pid+2] = d->far_fields[3*pid+2] + d->near_fields[3*pid+2];

      /* printf("%d: far_field=(%e, %e, %e)\n", pid,  */
      /* 	     d->far_fields[3*pid], */
      /* 	     d->far_fields[3*pid+1], */
      /* 	     d->far_fields[3*pid+2]); */
      /* printf("%d: near_field=(%e, %e, %e)\n", pid,  */
      /* 	     d->near_fields[3*pid], */
      /* 	     d->near_fields[3*pid+1], */
      /* 	     d->near_fields[3*pid+2]); */
    }

  if (potentials != NULL)
    for (fcs_int pid=0; pid < num_particles; pid++) {
      potentials[pid] = d->far_potentials[pid] + d->near_potentials[pid];
      /* printf("%d: potential=%e (near=%e, far=%e)\n", pid,  */
      /* 	     potentials[pid], d->near_potentials[pid], d->far_potentials[pid]); */
    }

  FCS_INFO(fprintf(stderr, "fcs_ewald_run() finished.\n"));
  return NULL;
}


FCSResult fcs_ewald_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  const char *fnc_name = "fcs_ewald_set_parameter";

  char *param = *current;
  char *cur = *next;

  *matched = 0;

  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_maxkmax", ewald_set_maxkmax, FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_kmax",    ewald_set_kmax,    FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_r_cut",   ewald_set_r_cut,   FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_alpha",   ewald_set_alpha,   FCS_PARSE_VAL(fcs_float));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}


FCSResult
fcs_ewald_print_parameters(FCS handle) {
  const char *fnc_name = "fcs_ewald_print_parameters";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  printf("ewald tolerance_field=%" FCS_LMOD_FLOAT "f\n", d->tolerance_field);

  if (d->needs_retune && d->tune_kmax)
    printf("ewald kmax=(to be tuned)\n");
  else if (d->tune_kmax)
    printf("ewald kmax=%" FCS_LMOD_INT "d (tuned)\n", d->kmax);
  else
    printf("ewald kmax=%" FCS_LMOD_INT "d\n", d->kmax);

  printf("ewald maxkmax=%" FCS_LMOD_INT "d\n", d->maxkmax);

  if (d->needs_retune && d->tune_r_cut)
    printf("ewald r_cut=(to be tuned)\n");
  else if (d->tune_r_cut)
    printf("ewald r_cut=%" FCS_LMOD_FLOAT "f (tuned)\n", d->r_cut);
  else
    printf("ewald r_cut=%" FCS_LMOD_FLOAT "f\n", d->r_cut);

  if (d->needs_retune && d->tune_alpha)
    printf("ewald alpha=(to be tuned)\n");
  else if (d->tune_alpha)
    printf("ewald alpha=%" FCS_LMOD_FLOAT "f (tuned)\n", d->alpha);
  else
    printf("ewald alpha=%" FCS_LMOD_FLOAT "f\n", d->alpha);
  return NULL;
}

/* safe free */
static void sfree(void* ptr) {
  if (ptr != NULL)
    free(ptr);
}

FCSResult fcs_ewald_destroy(FCS handle) {
  const char *fnc_name = "fcs_ewald_destroy";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  if (handle->method_context != NULL) {
    ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
    sfree(d->G);
    sfree(d->far_fields);
    sfree(d->near_fields);
    sfree(d->far_potentials);
    sfree(d->near_potentials);
    sfree(d);
  }
  return NULL;
}

