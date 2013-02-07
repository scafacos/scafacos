/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (C) 2012 Alexander Köwitsch
 * Copyright (C) 2011 Sebastian Banert
 * Copyright (C) 2011 Olaf Lenz
 * Copyright (C) 2010,2011 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include <stdio.h>
#include <stdlib.h>

#include "tune_0dp_noncubic.h"
#include "types.h"
#include "utils.h"
#include "constants.h"
#include "FCSCommon.h"

#include "kernels.h"
#include "regularization.h"
#include "cg_cos_coeff.h"
#include "cg_cos_err.h"

#define FCS_P2NFFT_DEBUG_TUNING 0
#define FCS_P2NFFT_EXIT_AFTER_TUNING 0
#define FCS_P2NFFT_TEST_GENERAL_ERROR_ESTIMATE 0
#define FCS_P2NFFT_ENABLE_TUNING_BUG 0

/* FORWARD DECLARATIONS OF STATIC FUNCTIONS */
static void init_near_interpolation_table_potential_periodic(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float alpha,
    fcs_float *table);
static void init_near_interpolation_table_force_periodic(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float alpha,
    fcs_float *table);

static fcs_int max_i(fcs_int a, fcs_int b);
static fcs_int calc_interpolation_num_nodes(
    fcs_int interpolation_order, fcs_float eps);

static void init_pnfft(
    FCS_PNFFT(plan) *ths, int dim, const ptrdiff_t *N, const ptrdiff_t *n,
    const fcs_float *x_max, int m,
    unsigned pnfft_flags, fcs_int pnfft_intpol_order, fcs_int pnfft_window,
    unsigned pfft_flags, fcs_int pfft_patience,
    MPI_Comm cart_comm_pnfft);
static int pnfft_is_up_to_date(
    const FCS_PNFFT(plan) ths, int dim, const ptrdiff_t *N, const ptrdiff_t *n,
    const fcs_float *x_max, int m, unsigned pnfft_flags, unsigned pfft_flags);

static void default_tolerance_type(
    fcs_int *periodicity,
    fcs_int *tolerance_type, fcs_float *tolerance);
static FCSResult check_tolerance(
    fcs_int *periodicity, fcs_int tolerance_type, fcs_float tolerance);

static fcs_float p2nfft_real_space_error(
    fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
    fcs_float r_cut, fcs_float alpha);
static fcs_float p2nfft_tune_alpha(
    fcs_int sum_qpart, fcs_float sum_q2, fcs_float box_l[3],
    fcs_float r_cut, ptrdiff_t grid[3], fcs_int cao, fcs_float tolerance_field);
static fcs_float p2nfft_get_accuracy(
    fcs_int sum_qpart, fcs_float sum_q2, fcs_float box_l[3],
    fcs_float r_cut, ptrdiff_t grid[3], fcs_int cao, fcs_float tolerance_field,
    fcs_float alpha,
    fcs_float *rs_err, fcs_float *ks_err);
static fcs_float p2nfft_k_space_error(
    fcs_int N, fcs_float sum_q2, fcs_float box_l[3], ptrdiff_t grid[3],
    fcs_float alpha, fcs_int cao);
static fcs_float p2nfft_k_space_error_sum1(
    fcs_int n, fcs_float grid_i, fcs_int cao);
static void p2nfft_k_space_error_sum2(
    fcs_int nx, fcs_int ny, fcs_int nz,
    ptrdiff_t grid[3], fcs_float grid_i[3],
    fcs_int cao, fcs_float alpha_L_i, 
    fcs_float *alias1, fcs_float *alias2);
static fcs_float p2nfft_k_space_error_approx(
    fcs_int N, fcs_float sum_q2,
    fcs_float box_l[3], ptrdiff_t grid[3],
    fcs_float alpha, fcs_int cao);

static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_periodic(
    const ptrdiff_t *local_N, const ptrdiff_t *local_N_start,
    fcs_float *box_l, fcs_float alpha);


FCSResult ifcs_p2nfft_tune_0dp_noncubic(
    void *rd, fcs_int *periodicity,
    fcs_int local_particles,
    fcs_float *positions, fcs_float *charges,
    fcs_float *box_l, fcs_int short_range_flag
    )
{
  fcs_float box_size = box_l[0];
  int comm_rank;
  const char* fnc_name = "ifcs_p2nfft_tune_3dp_noncubic";
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;
  fcs_int i, num_particles;
  fcs_float sum_q2, sum_q_abs, avg_dist, error=0;
  FCSResult result;

  /* Check consistency of user defined near field radius */
  if( !d->tune_epsI && !d->tune_r_cut )
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Nearfield radius was set in scaled ('epsI') and non-scaled ('r_cut') version. Choose one of them.");

  MPI_Comm_rank(d->cart_comm_3d, &comm_rank);

  d->short_range_flag = short_range_flag;

  /* Retune if total number of particles changed */
  MPI_Allreduce(&local_particles, &num_particles, 1, FCS_MPI_INT, MPI_SUM, d->cart_comm_3d);
  if (num_particles != d->num_nodes) {
    d->num_nodes = num_particles;
    d->needs_retune = 1;
  }

  /* Retune if periodicity changed */
  for(int t=0; t<3; t++){
    if((periodicity[t] == 0) != (d->periodicity[t] == 0))
      d->needs_retune = 1;
    d->periodicity[t] = periodicity[t];
  }
  d->use_ewald = periodicity[0] && periodicity[1] && periodicity[2];  

  /* Now, after the periodicity is clear, we can set the default tolerance type. */
  default_tolerance_type(d->periodicity,
      &d->tolerance_type, &d->tolerance);

  /* Check if P2NFFT can handle the tolerance type */
  result = check_tolerance(d->periodicity, d->tolerance_type, d->tolerance);
  if(result != NULL) return result;

  /* Calculate the sum of the squares of all charges
   * (needed for the error formulae) */
  fcs_float local_sum_q2 = 0;
  for (i = 0; i < local_particles; ++i)
    local_sum_q2 += FCS_P2NFFT_SQR(charges[i]);
  MPI_Allreduce(&local_sum_q2, &sum_q2, 1, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);
  if (!fcs_float_is_equal(sum_q2, d->sum_q2)) {
    d->sum_q2 = sum_q2;
    d->needs_retune = 1;
  }

  /* Calculate the sum of the absolute values of all charges
   * (needed for the error formulae) */
  fcs_float local_sum_q_abs = 0;
  for (i = 0; i < local_particles; ++i)
    local_sum_q_abs += fabs(charges[i]);
  MPI_Allreduce(&local_sum_q_abs, &sum_q_abs, 1, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);

  if (!fcs_float_is_equal(d->box_size, box_size)) {
    d->box_size = box_size;
    d->needs_retune = 1;
  }
  //d->box_l[0] = d->box_l[1] = d->box_l[2] = d->box_size; Non_cubic jetzt
  
  for(fcs_int t=0; t<3; t++)
    d->box_l[t] = box_l[t];

  /* FIXME: number of charged particles may be smaller then number of all particles */
  d->sum_qpart = num_particles;

  /* compute the average distance between two charges  */
  avg_dist = fcs_pow((d->box_l[0]*d->box_l[1]*d->box_l[2]) / d->sum_qpart, 0.33333);

  if (d->needs_retune) {
    fcs_float ks_error, rs_error;
    /* PNFFT calculates with real space cutoff 2*m+2
     * Therefore m is one smaller then the P3M cao. */
    fcs_int cao = 2*d->m;

    /* use full torus for periodic case */
    for(int t=0; t<3; t++)
      d->x_max[t] = 0.5;
   
    if(d->tune_r_cut){
      if(d->tune_epsI){
        /* set r_cut to 3 times the average distance between charges */
        d->r_cut = 3.0 * avg_dist;
        if (d->r_cut > d->box_size) 
          d->r_cut = d->box_size;   //perhaps something to change for noncubic, or better: what to change??
      } else
        d->r_cut = d->epsI * d->box_size;
    }

    /* shift and scale box into [-0.5,0.5)^3 */
    d->box_shift = d->box_size / 2.0;
    //d->box_scale = d->box_size;
    //for noncubic:
    d->box_scale = 1.0;
    for(int t=0; t<3; t++){
      d->box_scales[t] = d->box_l[t];
      d->box_shifts[t] = d->box_l[t] / 2.0;
    }

    /* set normalized near field radius */
    d->epsI = d->r_cut / d->box_scale;  //here the same, i am not sure
    
    /* Tune alpha for fixed N and m. */
    if(!d->tune_N && !d->tune_m){
      if(d->tune_alpha)
        d->alpha = p2nfft_tune_alpha(
            d->sum_qpart, d->sum_q2, d->box_l, d->r_cut, d->N, cao, d->tolerance);
      
      /* User specified N and cao. Therefore we do not necessarily need the error calculation. */
      if(~d->flags & FCS_P2NFFT_IGNORE_TOLERANCE){
        error = p2nfft_get_accuracy(
            d->sum_qpart, d->sum_q2, d->box_l, d->r_cut, d->N, cao,
            d->tolerance, d->alpha,
            &rs_error, &ks_error);

        /* Return error, if tuning of alpha failed. */
        if(error > d->tolerance)
          return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Not able to reach required accuracy (Did not find feasible Ewald splitting parameter 'alpha').");
      } else {
        /* calculation of real space error is fast */
        rs_error = p2nfft_real_space_error(d->sum_qpart, d->sum_q2, d->box_l, d->r_cut, d->alpha);
        /* Set error and ks_error to senseless value to indicate that user did not want to calculate the errors. */
        error = ks_error = -1.0;
      }
    } else {
      if(d->tune_m)
        cao = FCS_P2NFFT_MAXCAO;
   
      /* tune the grid size N */
      if(d->tune_N){
        for (i = FCS_P2NFFT_MINGRID; i <= FCS_P2NFFT_MAXGRID; i *= 2) {
          d->N[0] = i;
	  d->N[1] = 2*ceil(d->N[0]*d->box_l[1]/(d->box_l[0]*2.0));
	  d->N[2] = 2*ceil(d->N[0]*d->box_l[2]/(d->box_l[0]*2.0));
          if(d->tune_alpha)
            d->alpha = p2nfft_tune_alpha(
                d->sum_qpart, d->sum_q2, d->box_l, d->r_cut, d->N, cao, d->tolerance);
          error = p2nfft_get_accuracy(
              d->sum_qpart, d->sum_q2, d->box_l, d->r_cut, d->N, cao,
              d->tolerance, d->alpha,
              &rs_error, &ks_error);
          if (error < d->tolerance) break;
        }

        /* Return error, if tuning of N failed. */
        if(~d->flags & FCS_P2NFFT_IGNORE_TOLERANCE)
          if(error > d->tolerance)
            return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Not able to reach required accuracy (Did not find feasible gridsize 'N').");
      }

      /* tune m, only if we found a feasible gridsize N */
      if(d->tune_m){
        for (cao = 1; cao <= FCS_P2NFFT_MAXCAO; ++cao) {
          if(d->tune_alpha)
            d->alpha = p2nfft_tune_alpha(
                d->sum_qpart, d->sum_q2, d->box_l, d->r_cut, d->N, cao, d->tolerance);
          error = p2nfft_get_accuracy(
              d->sum_qpart, d->sum_q2, d->box_l, d->r_cut, d->N, cao, d->tolerance, d->alpha,
              &rs_error, &ks_error);
          if (error < d->tolerance) break;
        }

        /* Return error, if tuning of m failed. */
        if(~d->flags & FCS_P2NFFT_IGNORE_TOLERANCE)
          if(error > d->tolerance)
            return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Not able to reach required accuracy (Did not find feasible charge assignment order 'm').");
      }
    }
   
    /* Return error, if accuracy tuning failed. */
    if(~d->flags & FCS_P2NFFT_IGNORE_TOLERANCE)
      if(error > d->tolerance)
        return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Not able to reach required accuracy.");

    /* Bugfix: Although we compute with 2*m+2 values, the real space cutoff is still 2*m */
    /* P3M tuning works with cao==2*m */
    /* increase real space cutoff in order to compensate force errors because of analytic differentiation */
    d->m = (cao+1)/2;

    /* default oversampling equals 1 in 3d-periodic case */
    if(d->tune_n){
      for(int t=0; t<3; t++)
        d->n[t] = d->N[t];
    } else {
      /* check user defined oversampling */
      for(int t=0; t<3; t++)
        if(d->N[t] > d->n[t] )
          return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Oversampled gridsize is not allowed to be smaller than normal gridsize.");
    }

    /* Initialize the tables for near field interpolation */
    fcs_int interpolation_order = 3;
    d->interpolation_order = interpolation_order;
    /*   accuracy of 1e-16 needs 14000 interpolation nodes */
    /*   accuracy of 1e-17 needs 24896 interpolation nodes */
    d->interpolation_num_nodes = calc_interpolation_num_nodes(interpolation_order, 1e-16);

    if(d->near_interpolation_table_potential != NULL)
      free(d->near_interpolation_table_potential);
    d->near_interpolation_table_potential = (fcs_float*) malloc(sizeof(fcs_float) * (d->interpolation_num_nodes+3));
    init_near_interpolation_table_potential_periodic(
        d->interpolation_num_nodes,
        d->r_cut, d->alpha, 
        d->near_interpolation_table_potential);

    if(d->near_interpolation_table_force != NULL)
      free(d->near_interpolation_table_force);
    d->near_interpolation_table_force = (fcs_float*) malloc(sizeof(fcs_float) * (d->interpolation_num_nodes+3));
    init_near_interpolation_table_force_periodic(
        d->interpolation_num_nodes,
        d->r_cut, d->alpha,
        d->near_interpolation_table_force);

    /* calculate local data distribution according to PNFFT:
     * local_N, local_N_start, lower_border, upper_border */
    FCS_PNFFT(local_size_guru)(3, d->N, d->n, d->x_max, d->m, d->cart_comm_pnfft,
        d->local_N, d->local_N_start, d->lower_border, d->upper_border);
    
    /* shift and scale decomposition of the torus with box lengths */
    for(int t=0; t<3; t++){
      /* scale from [-x_max,x_max) to [-L/2,L/2) */
      d->lower_border[t] *= d->box_scales[t];
      d->upper_border[t] *= d->box_scales[t];

      /* shift from [-L/2,L/2) to [0,L) */
      d->lower_border[t] += d->box_shifts[t];
      d->upper_border[t] += d->box_shifts[t];

      /* local borders may be out of range because of box shift and scaling */
      if(d->lower_border[t] < 0)
        d->lower_border[t] = 0;
      if(d->upper_border[t] < 0)
        d->upper_border[t] = 0;
      if(d->lower_border[t] > d->box_l[t])
        d->lower_border[t] = d->box_l[t];
      if(d->upper_border[t] > d->box_l[t])
        d->upper_border[t] = d->box_l[t];
    }

    /* precompute Fourier coefficients for convolution */
    d->regkern_hat = malloc_and_precompute_regkern_hat_periodic(
        d->local_N, d->local_N_start, d->box_l, d->alpha);
  }

  /* Initialize the plan for the PNFFT */
  init_pnfft(&d->pnfft, 3, d->N, d->n, d->x_max, d->m, d->pnfft_flags, d->pnfft_interpolation_order, d->pnfft_window,
      d->pfft_flags, d->pfft_patience, d->cart_comm_pnfft);

  d->needs_retune = 0;

  /* Check if gridsize is too small to run with this number of processes. */
  int count=0;
  for(int t=0; t<3; t++)
    if(d->N[t] < d->np[t]) count ++;

  if(count>1)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Procmesh too large for this gridsize.");

  return NULL;
}

static void init_near_interpolation_table_potential_periodic(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float alpha,
    fcs_float *table
    )
{
  fcs_float r;

  for(fcs_int k=0; k<num_nodes+3; k++){
    r = r_cut * (fcs_float) k / num_nodes;
    if (fcs_float_is_zero(r))
      table[k] = 2 * alpha * FCS_P2NFFT_1OVERSQRTPI;
    else
      table[k] = erf(alpha * r)/r;
  }
}

static void init_near_interpolation_table_force_periodic(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float alpha,
    fcs_float *table
    )
{
  fcs_float r;

  for(fcs_int k=0; k<num_nodes+3; k++){
    r = r_cut * (fcs_float) k / num_nodes;
    if (fcs_float_is_zero(r))
      table[k] = 0;
    else
      table[k] = (-erf(alpha * r)/r
          + 2.0*alpha*FCS_P2NFFT_1OVERSQRTPI * fcs_exp(- alpha*alpha * r*r)
        ) / r;
  }
}

static fcs_int max_i(fcs_int a, fcs_int b)
{
  return a >= b ? a : b;
}

/* eps is the relative error */
static fcs_int calc_interpolation_num_nodes(
    fcs_int interpolation_order, fcs_float eps
    )
{
  /* use interpolation order 3 as default case */
  switch(interpolation_order){
    case 1: return (int) ceil(1.7/fcs_pow(eps,1.0/2.0));
    case 2: return (int) ceil(2.2/fcs_pow(eps,1.0/3.0));
    default: return max_i(10, (int) ceil(1.4/fcs_pow(eps,1.0/4.0)));
  }
}

static void init_pnfft(
    FCS_PNFFT(plan) *ths, int dim, const ptrdiff_t *N, const ptrdiff_t *n,
    const fcs_float *x_max, int m,
    unsigned pnfft_flags, fcs_int pnfft_intpol_order, fcs_int pnfft_window,
    unsigned pfft_flags, fcs_int pfft_patience,
    MPI_Comm cart_comm_pnfft
    )
{
  ptrdiff_t no_particles = 0;

#if FCS_ENABLE_INFO
  int myrank;
  MPI_Comm_rank(cart_comm_pnfft, &myrank);
#endif
  
  switch(pnfft_window){
    case 0: pnfft_flags |= PNFFT_WINDOW_KAISER_BESSEL; break;
    case 1: pnfft_flags |= PNFFT_WINDOW_GAUSSIAN; break;
    case 2: pnfft_flags |= PNFFT_WINDOW_BSPLINE; break;
    case 3: pnfft_flags |= PNFFT_WINDOW_SINC_POWER; break;
    case 4: pnfft_flags |= PNFFT_WINDOW_BESSEL_I0; break;
  }

  switch(pnfft_intpol_order){
    case 0: pnfft_flags |= 0; break;
    case 1: pnfft_flags |= PNFFT_PRE_LIN_PSI; break;
    case 2: pnfft_flags |= PNFFT_PRE_QUAD_PSI; break;
    case 3: pnfft_flags |= PNFFT_PRE_KUB_PSI; break;
  }

  switch(pfft_patience){
    case 0: pfft_flags |= FFTW_ESTIMATE; break;
    case 1: pfft_flags |= FFTW_MEASURE; break;
    case 2: pfft_flags |= FFTW_PATIENT; break;
    case 3: pfft_flags |= FFTW_EXHAUSTIVE; break;
    default: pfft_flags |= FFTW_MEASURE; break;
  }

  /* return if nothing to do */
  if( pnfft_is_up_to_date(*ths, dim, N, n, x_max, m, pnfft_flags, pfft_flags) )
    return;
#if FCS_P2NFFT_DEBUG_RETUNE
  else
    fprintf(stderr, "\n!!! pnfft_is_up_to_date fails !!!\n\n");
#endif


#if FCS_ENABLE_INFO
  if(!myrank){
    if(pnfft_flags & PNFFT_WINDOW_GAUSSIAN)
      printf("P2NFFT_INFO: Window function: gaussian (-c pnfft_window,gaussian)\n");
    else if(pnfft_flags & PNFFT_WINDOW_BSPLINE)
      printf("P2NFFT_INFO: Window function: bspline (-c pnfft_window,bspline)\n");
    else if(pnfft_flags & PNFFT_WINDOW_SINC_POWER)
      printf("P2NFFT_INFO: Window function: sinc power (-c pnfft_window,sinc)\n");
    else if(pnfft_flags & PNFFT_WINDOW_BESSEL_I0)
      printf("P2NFFT_INFO: Window function: bessel_i0 (-c pnfft_window,bessel_i0)\n");
    else
      printf("P2NFFT_INFO: Window function: kaiser-bessel (-c pnfft_window,kaiser)\n");

    if(pnfft_flags & PNFFT_FFT_OUT_OF_PLACE)
      printf("P2NFFT_INFO: inplace FFT: no (-c pnfft_fft_in_place,0)\n");
    else
      printf("P2NFFT_INFO: inplace FFT: yes (-c pnfft_fft_in_place,1)\n");
  
    if(pnfft_intpol_order==0)
      printf("P2NFFT_INFO: Interpolation order: direct evaluation of window function (-c pnfft_intpol_order,0)\n");
    else
      printf("P2NFFT_INFO: Interpolation order: %" FCS_LMOD_INT "d (-c pnfft_intpol_order,%" FCS_LMOD_INT "d)\n", pnfft_intpol_order, pnfft_intpol_order);
  }
#endif
#if FCS_ENABLE_INFO
  if(!myrank){
    printf("PNFFT_INIT: N = [%td, %td, %td]\n", N[0], N[1], N[2]);
    printf("PNFFT_INIT: n = [%td, %td, %td]\n", n[0], n[1], n[2]);
    printf("PNFFT_INIT: m = %" FCS_LMOD_INT "d\n", m);
    printf("PNFFT_INIT: pfft_flags = %u\n", pfft_flags);
    printf("PNFFT_INIT: pnfft_flags = %u\n", pnfft_flags);
    printf("PNFFT_INIT: x_max = [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", x_max[0], x_max[1], x_max[2]);
  }
#endif

  /* Start timing of PNFFT tuning */
  FCS_P2NFFT_INIT_TIMING(cart_comm_pnfft);
  FCS_P2NFFT_START_TIMING();
  
  /* finalize, if memory has been already allocated */
  if(*ths != NULL)
    FCS_PNFFT(finalize)(*ths, PNFFT_FREE_F_HAT| PNFFT_FREE_F| PNFFT_FREE_X| PNFFT_FREE_GRAD_F);

  /* call PNFFT planer */
  *ths = FCS_PNFFT(init_guru)(dim, N, n, x_max, no_particles, m, PNFFT_MALLOC_F_HAT| pnfft_flags, pfft_flags, cart_comm_pnfft);
  
  /* Finish timing of PNFFT tuning */
  FCS_P2NFFT_FINISH_TIMING(cart_comm_pnfft, "PNFFT tuning");
}

static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_periodic(
    const ptrdiff_t *local_N, const ptrdiff_t *local_N_start,
    fcs_float *box_l, fcs_float alpha
    )
{
  ptrdiff_t m, k0, k1, k2, alloc_local = 1;
  fcs_float ksqnorm;
  fcs_pnfft_complex *regkern_hat;
 
  for(int t=0; t<3; t++)
    alloc_local *= local_N[t];
  regkern_hat = FCS_PFFT(alloc_complex)(alloc_local);
  
  m = 0;
  for(k1=local_N_start[1]; k1 < local_N_start[1] + local_N[1]; k1++){
    for(k2=local_N_start[2]; k2 < local_N_start[2] + local_N[2]; k2++){
      for(k0=local_N_start[0]; k0 < local_N_start[0] + local_N[0]; k0++, m++){
        ksqnorm = k0*(fcs_float)k0/(box_l[0]*box_l[0]) + k1*(fcs_float)k1/(box_l[1]*box_l[1]) + k2*(fcs_float)k2/(box_l[2]*box_l[2]);
        if ((k0 == 0) && (k1 == 0) && (k2 == 0))
          regkern_hat[m] = 0;
        else
          regkern_hat[m] = fcs_exp(-FCS_P2NFFT_PISQR * ksqnorm/(alpha*alpha))/ksqnorm;
      }
    }
  }
  
  return regkern_hat;
}


static int pnfft_is_up_to_date(
    const FCS_PNFFT(plan) ths, int dim, const ptrdiff_t *N, const ptrdiff_t *n,
    const fcs_float *x_max, int m, unsigned pnfft_flags, unsigned pfft_flags
    )
{
  int plan_d, plan_m;
  ptrdiff_t *plan_N, *plan_n;
  fcs_float *plan_x_max;
  unsigned plan_pnfft_flags, plan_pfft_flags;

  /* plan is uninitialized */
  if( ths == NULL )
    return 0;

  plan_d = FCS_PNFFT(get_d)(ths);
  plan_m = FCS_PNFFT(get_m)(ths);
  plan_N = FCS_PNFFT(get_N)(ths);
  plan_n = FCS_PNFFT(get_n)(ths);
  plan_x_max = FCS_PNFFT(get_x_max)(ths);
  plan_pnfft_flags = FCS_PNFFT(get_pnfft_flags)(ths);
  plan_pfft_flags = pnfft_get_pfft_flags(ths);

  /* check values of d, m */
  if( plan_d != dim
      || plan_m != m
    )
    return 0;

  /* check values of N, n, x_max */
  for(int t=0; t<dim; t++)
    if( plan_N[t] != N[t]
        || plan_n[t] != n[t]
        || !fcs_float_is_zero(plan_x_max[t] - x_max[t])
      )
      return 0;

  /* remove flags that do not imply retuning */
  plan_pnfft_flags &= ~(PNFFT_MALLOC_X| PNFFT_MALLOC_F| PNFFT_MALLOC_GRAD_F);

  /* check PNFFT flags */
  if(plan_pnfft_flags != pnfft_flags)
    return 0;

  /* check PFFT flags */
  if(plan_pfft_flags != pfft_flags)
    return 0;

  return 1;  
}

static void default_tolerance_type(
    fcs_int *periodicity,
    fcs_int *tolerance_type, fcs_float *tolerance
    )
{
  if(periodicity[0] && periodicity[1] && periodicity[2]){
    if(*tolerance < 0.0)
      *tolerance = FCS_P2NFFT_DEFAULT_TOLERANCE;
    if(*tolerance_type == FCS_TOLERANCE_TYPE_UNDEFINED)
      *tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
  } else if(!periodicity[0] && !periodicity[1] && !periodicity[2]){
    if(*tolerance < 0.0)
      *tolerance = FCS_P2NFFT_DEFAULT_TOLERANCE;
    if(*tolerance_type == FCS_TOLERANCE_TYPE_UNDEFINED)
      *tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL;
  }
}

static FCSResult check_tolerance(
    fcs_int *periodicity, fcs_int tolerance_type, fcs_float tolerance
    )
{
  char* fnc_name = "check_tolerance";

  if( (tolerance_type == FCS_TOLERANCE_TYPE_ENERGY)
      || (tolerance_type == FCS_TOLERANCE_TYPE_ENERGY_REL)
      || (tolerance_type == FCS_TOLERANCE_TYPE_POTENTIAL_REL)
      || (tolerance_type == FCS_TOLERANCE_TYPE_FIELD_REL) )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name,"P2NFFT does not support this kind of tolerance.");

  if(tolerance_type == FCS_TOLERANCE_TYPE_POTENTIAL)
    if( periodicity[0] || periodicity[1] || periodicity[2] )
      return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name,"P2NFFT supports FCS_TOLERANCE_POTENTIAL only for non-periodic boundary conditions. Use FCS_TOLERANCE_FIELD instead.");

  if(tolerance_type == FCS_TOLERANCE_TYPE_FIELD)
    if( !periodicity[0] || !periodicity[1] || !periodicity[2] )
      return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name,"P2NFFT supports FCS_TOLERANCE_FIELD only for 3d-periodic boundary conditions. Use FCS_TOLERANCE_POTENTIAL instead.");

  if(tolerance < 0.0)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name,"Tolerance must be non-negative.");

  return NULL;
}









/* Calculates the real space contribution to the rms error in the force
 * (as described by Kolafa and Perram).
 * \param N the number of charged particles in the system.
 * \param sum_q2 the sum of square of charges in the system
 * \param box_l fcs_float[3] system size
 * \param r_cut the P2NFFT cutoff
 * \param alpha the P2NFFT alpha (Ewald splitting parameter)
 * \return the real space error
 */
static fcs_float p2nfft_real_space_error(
    fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
    fcs_float r_cut, fcs_float alpha
    )
{
  return (2.0*sum_q2*fcs_exp(-FCS_P2NFFT_SQR(r_cut*alpha))) / (sqrt((fcs_float)N*r_cut*box_l[0]*box_l[1]*box_l[2]));
}

/* Get the error for this combination of parameters and tune alpha if
 * required. In fact, the real space error is tuned such that it
 * contributes half of the total error, and then the Fourier space
 * error is calculated. Returns the achieved error and the optimal
 * alpha.
 */
// static fcs_float p2nfft_get_accuracy(
//     fcs_int sum_qpart, fcs_float sum_q2, fcs_float box_l[3],
//     fcs_float epsI, ptrdiff_t N[3], fcs_int m, fcs_float tolerance_field,
//     fcs_int tune_alpha, MPI_Comm comm,
//     fcs_float *alpha, fcs_float *rs_error, fcs_float *ks_error
//     )
// {
//   fcs_int grid[3];
//   fcs_float error;
// 
//   for(int t=0; t<3; t++)
//     grid[t] = (fcs_int) N[t];
// 
//   if(tune_alpha)
//     ifcs_fftcommon_determine_good_alpha(
//         sum_qpart, sum_q2, box_l, epsI, grid, m, tolerance_field,
//         alpha, &error, rs_error, ks_error, comm);
//   else
//     ifcs_fftcommon_compute_error_estimate(
//         sum_qpart, sum_q2, box_l, epsI, grid, *alpha, m,
//         &error, rs_error, ks_error, comm);
// 
//   return error;
// }
static fcs_float p2nfft_tune_alpha(
    fcs_int sum_qpart, fcs_float sum_q2, fcs_float box_l[3],
    fcs_float r_cut, ptrdiff_t grid[3], fcs_int cao, fcs_float tolerance_field
    )
{
  fcs_float alpha, rs_err; 
  
  /* calc the maximal real space error for the setting (i.e., set alpha to 0) */
  rs_err = p2nfft_real_space_error(sum_qpart, sum_q2, box_l, r_cut, 0.0);

  if(FCS_P2NFFT_SQRT2 * rs_err > tolerance_field) {
    /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
    alpha = sqrt(log(FCS_P2NFFT_SQRT2 * rs_err/tolerance_field)) / r_cut;
  } else {
    /* even alpha=0 is ok, however, we cannot choose it since it kills the k-space error formula.
     * Anyways, this is very likely NOT the optimal solution */
    alpha = 0.1 / box_l[0];
  }

  return alpha;
}

static fcs_float p2nfft_get_accuracy(
    fcs_int sum_qpart, fcs_float sum_q2, fcs_float box_l[3],
    fcs_float r_cut, ptrdiff_t grid[3], fcs_int cao, fcs_float tolerance_field,
    fcs_float alpha,
    fcs_float *rs_err, fcs_float *ks_err
    )
{
  /* calculate real space and k space error for this alpha */
  *rs_err = p2nfft_real_space_error(sum_qpart, sum_q2, box_l, r_cut, alpha);

  fcs_int full_estimate = 0;
  for(fcs_int t = 0; t < 3 && !full_estimate; t++){
    fcs_float alpha_h = alpha * box_l[t] / grid[t];
    full_estimate = alpha_h > FCS_P2NFFT_FULL_ESTIMATE_ALPHA_H_THRESHOLD;
  }
  
//full_estimate=0;

  if (full_estimate)
    *ks_err = p2nfft_k_space_error(sum_qpart, sum_q2, box_l, grid, alpha, cao);
  else
    *ks_err = p2nfft_k_space_error_approx(sum_qpart, sum_q2, box_l, grid, alpha, cao);

  return sqrt(FCS_P2NFFT_SQR(*rs_err)+FCS_P2NFFT_SQR(*ks_err));
}

/* Calculate the analytic expression of the error estimate for the
 * P2NFFT method in the book of Hockney and Eastwood (Eqn. 8.23) in
 * order to obtain the rms error in the force for a system of N
 * randomly distributed particles in a cubic box (k space part).
 * \param grid     number of grid points in one direction.
 * \param cao      charge assignment order.
 * \param n_c_part number of charged particles in the system.
 * \param sum_q2   sum of square of charges in the system
 * \param alpha    ewald splitting parameter.
 * \return reciprocal (k) space error
*/
static fcs_float p2nfft_k_space_error(
    fcs_int N, fcs_float sum_q2,
    fcs_float box_l[3], ptrdiff_t grid[3],
    fcs_float alpha, fcs_int cao
    )
{
  fcs_float he_q = 0.0;
  fcs_float grid_i[3] = {1.0/grid[0], 1.0/grid[1], 1.0/grid[2]};
  /* @todo Handle non-cubic case  */
  fcs_float alpha_L_i = 1./(alpha*box_l[0]);

#if FCS_P2NFFT_DEBUG_TUNING
  fprintf(stderr, "P2NFFT_DEBUG_TUNING: N = %d, grid_i = [%f, %f, %f], alpha_L_i = %f\n",
      N, grid_i[0], grid_i[1], grid_i[2], alpha_L_i);
#endif

  for (fcs_int nx=-grid[0]/2; nx<grid[0]/2; nx++) {
    fcs_float ctan_x = p2nfft_k_space_error_sum1(nx,grid_i[0],cao);
    for (fcs_int ny=-grid[1]/2; ny<grid[1]/2; ny++) {
      fcs_float ctan_y = ctan_x * p2nfft_k_space_error_sum1(ny,grid_i[1],cao);
      for (fcs_int nz=-grid[2]/2; nz<grid[2]/2; nz++) {
	if((nx!=0) || (ny!=0) || (nz!=0)) {
	  fcs_float n2 = FCS_P2NFFT_SQR(nx) + FCS_P2NFFT_SQR(ny) + FCS_P2NFFT_SQR(nz);
	  fcs_float cs = p2nfft_k_space_error_sum1(nz,grid_i[2],cao)*ctan_y;
	  fcs_float alias1, alias2;
	  p2nfft_k_space_error_sum2(nx,ny,nz,grid,grid_i,cao,alpha_L_i,&alias1,&alias2);
	  fcs_float d = alias1  -  FCS_P2NFFT_SQR(alias2/cs) / n2;
	  /* at high precisions, d can become negative due to extinction;
	     also, don't take values that have no significant digits left*/
	  if (d > 0.0 && (!fcs_float_is_zero(fabs(d/alias1))))
	    he_q += d;
// if( (nx==1) && (ny==0) && (nz==0) )
//   fprintf(stderr, "P3M: k = [%d, %d, %d], N = [%td, %td, %td], alias_k = %e, alias1 = %e, alias2 = %e, cs = %e, n2 = %e, alias_minus = %e, ctan_x = %.16e, ctan_y = %e, ctan_z = %e\n",
//       nx, ny, nz, grid[0], grid[1], grid[2], d, alias1, alias2, cs, n2, FCS_P2NFFT_SQR(alias2/cs) / n2, ctan_x, ctan_y/ctan_x, cs/ctan_y);
	}
      }
    }
  }

#if FCS_P2NFFT_DEBUG_TUNING
  fprintf(stderr, "P2NFFT_DEBUG_TUNING: P3M: he_q = %e\n", he_q);
#endif

#if FCS_P2NFFT_ENABLE_TUNING_BUG
  return 2.0*sum_q2*sqrt(he_q/(fcs_float)N) / (box_l[1]*box_l[2]);
#else
  /* Where does the factor 2.0 come from? */
  return 2.0*sum_q2*sqrt( he_q / (fcs_float)N / (box_l[0]*box_l[1]*box_l[2]) );
#endif
}

/* One of the aliasing sums used by \ref p2nfft_k_space_error.
 * (fortunately the one which is most important (because it converges
 * most slowly, since it is not damped exponentially)) can be
 * calculated analytically. The result (which depends on the order of
 * the spline interpolation) can be written as an even trigonometric
 * polynomial. The results are tabulated here (The employed formula
 * is Eqn. 7.66 in the book of Hockney and Eastwood). */
static fcs_float p2nfft_k_space_error_sum1(
    fcs_int n, fcs_float grid_i, fcs_int cao
    )
{
  fcs_float c, res=0.0;
  c = FCS_P2NFFT_SQR(cos(FCS_P2NFFT_PI*grid_i*(fcs_float)n));
  
  switch (cao) {
    case 1 : { 
      res = 1; 
      break; }
    case 2 : { 
      res = (1.0+c*2.0)/3.0; 
      break; }
    case 3 : { 
      res = (2.0+c*(11.0+c*2.0))/15.0; 
      break; }
    case 4 : { 
      res = (17.0+c*(180.0+c*(114.0+c*4.0)))/315.0; 
      break; }
    case 5 : { 
      res = (62.0+c*(1072.0+c*(1452.0+c*(247.0+c*2.0))))/2835.0; 
      break; }
    case 6 : { 
      res = (1382.0+c*(35396.0+c*(83021.0+c*(34096.0+c*(2026.0+c*4.0)))))/155925.0; 
      break; }
    case 7 : { 
      res = (21844.0+c*(776661.0+c*(2801040.0+c*(2123860.0+c*(349500.0+c*(8166.0+c*4.0))))))/6081075.0; 
      break; }
    default : {
      fprintf(stderr,"INTERNAL_ERROR: The value %" FCS_LMOD_INT "d for the interpolation order should not occur!\n", cao);
      exit(1);
    }
  }
  
  return res;
}


/* aliasing sum used by \ref p2nfft_k_space_error. */
static void p2nfft_k_space_error_sum2(
    fcs_int nx, fcs_int ny, fcs_int nz,
    ptrdiff_t grid[3], fcs_float grid_i[3],
    fcs_int cao, fcs_float alpha_L_i,
    fcs_float *alias1, fcs_float *alias2
    )
{
  fcs_float prefactor = FCS_P2NFFT_SQR(FCS_P2NFFT_PI*alpha_L_i);

  *alias1 = *alias2 = 0.0;
  for (fcs_int mx=-FCS_P2NFFT_BRILLOUIN; mx<=FCS_P2NFFT_BRILLOUIN; mx++) {
    fcs_float nmx = nx + mx * (fcs_float)grid[0];
    fcs_float fnmx = grid_i[0] * nmx;
    for (fcs_int my=-FCS_P2NFFT_BRILLOUIN; my<=FCS_P2NFFT_BRILLOUIN; my++) {
      fcs_float nmy = ny + my * (fcs_float)grid[1];
      fcs_float fnmy = grid_i[1] * nmy;
      for (fcs_int mz=-FCS_P2NFFT_BRILLOUIN; mz<=FCS_P2NFFT_BRILLOUIN; mz++) {
	fcs_float nmz = nz + mz * (fcs_float)grid[2];
	fcs_float fnmz = grid_i[2] * nmz;

	fcs_float nm2 = FCS_P2NFFT_SQR(nmx) + FCS_P2NFFT_SQR(nmy) + FCS_P2NFFT_SQR(nmz);
	fcs_float ex = fcs_exp(-prefactor*nm2);
	
	fcs_float U2 = fcs_pow(sinc(fnmx)*sinc(fnmy)*sinc(fnmz), 2.0*cao);
	
	*alias1 += ex*ex / nm2;
	*alias2 += U2 * ex * (nx*nmx + ny*nmy + nz*nmz) / nm2;
      }
    }
  }
}

/* Calculate the analytical approximation for the k-space part of the
 * error (Eq. 38 in Deserno, Holm; JCP 109,18; 1998). */
static fcs_float p2nfft_k_space_error_approx(
    fcs_int N, fcs_float sum_q2,
    fcs_float box_l[3], ptrdiff_t grid[3],
    fcs_float alpha, fcs_int cao
    )
{
  /* grid spacing */
  /* TODO: non-cubic case*/
  fcs_float h = box_l[0]/grid[0];
  fcs_float ha = h*alpha;

  /* compute the sum in eq. 38 */
  fcs_float sum;
  switch (cao) {
    case 1:
      sum = 2./3.; 
      break;
    case 2:
      sum = 5./294.*fcs_pow(ha,2) + 1./50.; 
      break;
    case 3:
      sum = 
        21./3872.*fcs_pow(ha, 4) + 7./1440.*fcs_pow(ha,2) + 
        1./588.;
      break;
    case 4:
      sum = 
        143./28800.*fcs_pow(ha, 6) + 
        7601./2271360.*fcs_pow(ha, 4) + 3./1936.*fcs_pow(ha,2) + 1./4320.;
      break;
    case 5:
      sum = 
        106640677./11737571328.*fcs_pow(ha,8) + 517231./106536960.*fcs_pow(ha, 6) + 
        143./69120.*fcs_pow(ha, 4) + 7601./13628160.*fcs_pow(ha,2) + 
        1./23232.;
      break;
    case 6:
      sum = 
        326190917./11700633600.*fcs_pow(ha,10) + 
        733191589./59609088000.*fcs_pow(ha,8) + 9694607./2095994880.*fcs_pow(ha, 6) + 
        47021./35512320.*fcs_pow(ha, 4) + 13./57600.*fcs_pow(ha,2) + 
        691./68140800.;
      break;
    case 7:
      sum =
        4887769399./37838389248.*fcs_pow(ha,12) + 1755948832039./36229939200000.*fcs_pow(ha,10) + 
        25091609./1560084480.*fcs_pow(ha,8) + 56399353./12773376000.*fcs_pow(ha, 6) + 
        745739./838397952.*fcs_pow(ha, 4) + 3617./35512320.*fcs_pow(ha,2) + 
        1./345600.;
      break;
    default: 
      fprintf(stderr,"INTERNAL_ERROR: p2nfft_analytical_k_space_error_sum: Charge assignment order of %" FCS_LMOD_INT "d should not occur!\n", cao);
      exit(1);
  };

#if FCS_P2NFFT_ENABLE_TUNING_BUG
  return 
    sum_q2 / (box_l[0]*box_l[0] + box_l[1]*box_l[1] + box_l[2]*box_l[2]) *
    fcs_pow(h*alpha, cao) * 
    sqrt(alpha*box_l[0]/N*sqrt(2.0*FCS_P2NFFT_PI)*sum);
#else
  return 
    sum_q2 / (box_l[0]*box_l[0]) *
    fcs_pow(h*alpha, cao) * 
    sqrt(alpha*box_l[0]/N*sqrt(2.0*FCS_P2NFFT_PI)*sum);
#endif
}




