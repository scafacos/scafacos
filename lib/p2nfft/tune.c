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

#include "tune.h"
#include "types.h"
#include "utils.h"
#include "constants.h"
#include "FCSCommon.h"

#include "kernels.h"
#include "regularization.h"
#include "interpolation.h"
#include "taylor2p.h"
#include "cg_cos_coeff.h"
#include "cg_cos_err.h"
#include "cg_cos_coeff_sym.h"
#include "cg_cos_err_sym.h"

#include "bessel_k.h"
#include "part_derive_one_over_norm_x.h"

#define FCS_P2NFFT_DEBUG_TUNING 0
#define FCS_P2NFFT_EXIT_AFTER_TUNING 0
#define FCS_P2NFFT_TEST_GENERAL_ERROR_ESTIMATE 0
#define FCS_P2NFFT_ENABLE_TUNING_BUG 0
#define FCS_P2NFFT_DEBUG_REGKERN 0

#if FCS_P2NFFT_DEBUG_REGKERN
#  define FCS_P2NFFT_IFDBG_REGKERN(code) code
#else
#  define FCS_P2NFFT_IFDBG_REGKERN(code) 
#endif

/* compute d-th component of A^T * v */
#define At_TIMES_VEC(_A_, _v_, _d_) ( (_v_)[0] * (_A_)[_d_] + (_v_)[1] * (_A_)[_d_ + 3] + (_v_)[2] * (_A_)[_d_ + 6] )

/* FORWARD DECLARATIONS OF STATIC FUNCTIONS */
static void print_command_line_arguments(
    ifcs_p2nfft_data_struct *d, fcs_int verbose);
static void init_near_interpolation_table_potential_3dp(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float alpha,
    fcs_float *table);
static void init_near_interpolation_table_force_3dp(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float alpha,
    fcs_float *table);

static void init_near_interpolation_table_potential_0dp(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float epsI, fcs_int p,
    const fcs_float *taylor2p_coeff,
    fcs_int N_cos, const fcs_float *cos_coeff,
    fcs_float *table);
static void init_far_interpolation_table_potential_0dp(
    fcs_int num_nodes, fcs_int reg_far,
    fcs_float epsB, fcs_int p, fcs_float c,
    fcs_int N_cos, const fcs_float *cos_coeff,
    fcs_float *table);
static void init_far_interpolation_table_potential_0dp_ewald(
    fcs_int num_nodes, fcs_int reg_far,
    fcs_float alpha, fcs_float h,
    fcs_float epsB, fcs_int p, fcs_float c,
    fcs_float *table); 
static void init_far_interpolation_table_potential_1dp(
    fcs_int num_nodes, fcs_int reg_far,
    fcs_int k, fcs_float periodic_box_length, 
    fcs_float alpha, fcs_float h,
    fcs_float epsB, fcs_int p, fcs_float c,
    fcs_float *table );
static void init_near_interpolation_table_force_0dp(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float epsI, fcs_int p,
    const fcs_float *taylor2p_derive_coeff,
    fcs_int N_cos, fcs_float *cos_coeff, fcs_float *sin_coeff,
    fcs_float *table);
static fcs_int max_i(fcs_int a, fcs_int b);
static fcs_int calc_interpolation_num_nodes(
    fcs_int interpolation_order, fcs_float eps);
static fcs_int calc_interpolation_num_nodes_erf(
    fcs_int interpolation_order, fcs_float eps, fcs_float alpha, fcs_float r, unsigned *err);
static fcs_float evaluate_cos_polynomial_1d(
   fcs_float x, fcs_int N, const fcs_float *coeff);
static fcs_float evaluate_sin_polynomial_1d(
   fcs_float x, fcs_int N, const fcs_float *coeff);
static int get_dim_of_smallest_periodic_box_l(
    fcs_int periodicity[3], fcs_float box_l[3]);

static void init_pnfft(
    FCS_PNFFT(plan) *ths, int dim, const ptrdiff_t *N, const ptrdiff_t *n,
    const fcs_float *x_max, int m,
    unsigned pnfft_flags, fcs_int pnfft_intpol_order, fcs_int pnfft_window,
    unsigned pfft_flags, fcs_int pfft_patience,
    MPI_Comm cart_comm_pnfft);
static int pnfft_is_up_to_date(
    const FCS_PNFFT(plan) ths, int dim, const ptrdiff_t *N, const ptrdiff_t *n,
    const fcs_float *x_max, int m, unsigned pnfft_flags, unsigned pfft_flags);
static int reg_far_is_radial(
    fcs_int reg_far);

static void default_tolerance_type(
    fcs_int *periodicity, fcs_int reg_kernel,
    fcs_int *tolerance_type, fcs_float *tolerance);
static FCSResult check_tolerance(
    fcs_int *periodicity, fcs_int reg_kernel, fcs_int tolerance_type, fcs_float tolerance);

static fcs_float p2nfft_real_space_error(
    fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
    fcs_float r_cut, fcs_float alpha);
static fcs_float p2nfft_tune_alpha(
    fcs_int sum_qpart, fcs_float sum_q2, fcs_int dim_tune,
    fcs_float box_l[3], fcs_float r_cut, ptrdiff_t grid[3],
    fcs_int cao, fcs_float tolerance_field, fcs_int plain_ewald_flag);
static fcs_float p2nfft_get_accuracy(
    fcs_int sum_qpart, fcs_float sum_q2, fcs_int dim_tune,
    fcs_float box_l[3], fcs_float r_cut, ptrdiff_t grid[3],
    fcs_int cao, fcs_float tolerance_field,
    fcs_float alpha, fcs_int interlaced, fcs_int plain_ewald_flag,
    fcs_float *rs_err, fcs_float *ks_err);
static fcs_float p2nfft_k_space_error(
    fcs_int N, fcs_float sum_q2, fcs_int dim_tune,
    fcs_float box_l[3], ptrdiff_t grid[3],
    fcs_float alpha, fcs_int cao,
    fcs_int interlaced);
static fcs_float p2nfft_k_space_error_sum1(
    fcs_int n, fcs_float grid_i, fcs_int cao);
static void p2nfft_k_space_error_sum2_ad(
    fcs_int nx, fcs_int ny, fcs_int nz,
    ptrdiff_t grid[3], fcs_float grid_i[3],
    fcs_int cao, fcs_float alpha_L_i, 
    fcs_float *alias1, fcs_float *alias2);
static void p2nfft_k_space_error_sum2_adi(
    fcs_int nx, fcs_int ny, fcs_int nz,
    ptrdiff_t grid[3], fcs_float grid_i[3],
    fcs_int cao, fcs_float alpha_L_i,
    fcs_float *alias1, fcs_float *alias2,
    fcs_float *alias3, fcs_float *alias4,
    fcs_float *alias5, fcs_float *alias6);
static fcs_float p2nfft_k_space_error_approx(
    fcs_int N, fcs_float sum_q2, fcs_int dim_tune,
    fcs_float box_l[3], ptrdiff_t grid[3],
    fcs_float alpha, fcs_int cao);
#if FCS_P2NFFT_TEST_GENERAL_ERROR_ESTIMATE
static fcs_float p2nfft_k_space_error_general_window(
    fcs_int num_part, fcs_float sum_q2,
    const fcs_float box_l[3], const ptrdiff_t N[3],
    fcs_float alpha, fcs_int m,
    unsigned window_flag,
    MPI_Comm comm);
static fcs_float compute_alias_k(
    FCS_PNFFT(plan) *wind_param, fcs_float alpha, const fcs_float box_l[3],
    const ptrdiff_t N[3], const ptrdiff_t k[3]);

#endif

static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_0dp_ewald(
    const ptrdiff_t *N, fcs_float epsB,
//    const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c,
//    const fcs_float *box_inv,
    const fcs_float *box_scales, fcs_float alpha,
    fcs_int p, fcs_float c, fcs_int reg_far,
    fcs_int interpolation_order, fcs_int far_interpolation_num_nodes, const fcs_float *far_interpolation_table_potential,
    MPI_Comm comm_cart, unsigned box_is_cubic
    );
static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_0dp(
    const ptrdiff_t *N, fcs_float r_cut, fcs_float epsI, fcs_float epsB,
    fcs_int p, fcs_float c, fcs_float *box_scales,
    fcs_int reg_near, fcs_int reg_far,
    const fcs_float *taylor2p_coeff, fcs_int N_cos, const fcs_float *cos_coeff,
    fcs_int interpolation_order, fcs_int near_interpolation_num_nodes, fcs_int far_interpolation_num_nodes,
    const fcs_float *near_interpolation_table_potential, const fcs_float *far_interpolation_table_potential,
    MPI_Comm comm_cart, unsigned box_is_cubic);
static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_2dp_and_1dp(
    const ptrdiff_t *N, fcs_float epsB,
    const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c,
    const fcs_float *box_inv, const fcs_float *box_scales, fcs_float alpha, fcs_float kc,
    const fcs_int *periodicity, fcs_int p, fcs_float c, fcs_int reg_far,
    fcs_int interpolation_order, fcs_int far_interpolation_num_nodes, 
    const fcs_float *far_interpolation_table_potential,
    MPI_Comm comm_cart);
static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_3dp(
    const ptrdiff_t *local_N, const ptrdiff_t *local_N_start,
    const fcs_float *ibox, fcs_float alpha, fcs_float kc);

static fcs_int is_cubic(
    fcs_float *box_l);


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
ewald_real_space_error(
    fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
    fcs_float r_cut, fcs_float alpha,
    fcs_float *error, fcs_float *derivative
    )
{
  const fcs_float V = box_l[0]*box_l[1]*box_l[2];
  *error = 2.0*sum_q2 / fcs_sqrt(N*r_cut*V) * fcs_exp(- alpha*alpha * r_cut*r_cut);
  *derivative = -2.0 * alpha * r_cut*r_cut * *error;
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
ewald_k_space_error(
    fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
    fcs_int kmax, fcs_float alpha,
    fcs_float *error, fcs_float *derivative
    )
{
  fcs_float Lmax = box_l[0];
  if (box_l[1] > Lmax) Lmax = box_l[1];
  if (box_l[2] > Lmax) Lmax = box_l[2];
  const fcs_float K = 2.0*M_PI*kmax/Lmax;
  const fcs_float V = box_l[0]*box_l[1]*box_l[2];
  const fcs_float fak1 = 2.0*sum_q2 / fcs_sqrt(N*M_PI*K*V) * fcs_exp(-K*K/(4.0*alpha*alpha));
  *error = fak1 * alpha;
  *derivative = fak1 * (2.0*alpha*alpha + kmax*kmax/(2.0*alpha*alpha));
}

static inline ptrdiff_t
get_kmax(
    ptrdiff_t gridsize
    )
{
  return gridsize/2-1;
}

static inline fcs_float
ewald_tune_alpha(
    fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
    fcs_float r_cut, fcs_int kmax
    )
{
  fcs_float alpha;
  fcs_float err_r, der_r;
  fcs_float err_k, der_k;
  fcs_float alpha_diff;

  /* Newton-Raphson method to find optimal alpha */
  alpha = 0.1;
  do {
    /* get errors and derivatives */
    ewald_real_space_error(N, sum_q2, box_l, r_cut, alpha, &err_r, &der_r);
    ewald_k_space_error(N, sum_q2, box_l, kmax, alpha, &err_k, &der_k);
    
    alpha_diff = (err_r - err_k) / (der_r - der_k);
    alpha -= alpha_diff;
  } while (fcs_fabs(alpha_diff) > FCS_P2NFFT_ALPHA_OPT_PREC);

  return alpha;
}






static fcs_float det_3x3(
    const fcs_float *v0, const fcs_float *v1, const fcs_float *v2
    )
{
  return v0[0] * v1[1] * v2[2] + v1[0] * v2[1] * v0[2] + v2[0] * v0[1] * v1[2] - v2[0] * v1[1] * v0[2] - v1[0] * v0[1] * v2[2] - v0[0] * v2[1] * v1[2];
}

/* get the columns of an 3x3 matrix and output the inverse */
static void invert_3x3(const fcs_float *v0, const fcs_float *v1, const fcs_float *v2, fcs_float *iv)
{
  fcs_float det = det_3x3(v0,v1,v2);

  iv[0] = (v1[1] * v2[2] - v2[1] * v1[2]) / det;
  iv[1] = (v2[0] * v1[2] - v1[0] * v2[2]) / det;
  iv[2] = (v1[0] * v2[1] - v2[0] * v1[1]) / det;

  iv[3] = (v2[1] * v0[2] - v0[1] * v2[2]) / det;
  iv[4] = (v0[0] * v2[2] - v2[0] * v0[2]) / det;
  iv[5] = (v2[0] * v0[1] - v0[0] * v2[1]) / det;

  iv[6] = (v0[1] * v1[2] - v1[1] * v0[2]) / det;
  iv[7] = (v1[0] * v0[2] - v0[0] * v1[2]) / det;
  iv[8] = (v0[0] * v1[1] - v1[0] * v0[1]) / det;
}

static fcs_float volume_of_periodic_dims(
    const fcs_float *a, const fcs_float *b, const fcs_float *c, const fcs_int *periodicity
    )
{
  fcs_float v0[3], v1[3], v2[3];

  v0[0] = periodicity[0] ? a[0] : 1.0;
  v0[1] = periodicity[0] ? a[1] : 0.0;
  v0[2] = periodicity[0] ? a[2] : 0.0;

  v1[0] = periodicity[1] ? b[0] : 0.0;
  v1[1] = periodicity[1] ? b[1] : 1.0;
  v1[2] = periodicity[1] ? b[2] : 0.0;

  v2[0] = periodicity[2] ? c[0] : 0.0;
  v2[1] = periodicity[2] ? c[1] : 0.0;
  v2[2] = periodicity[2] ? c[2] : 1.0;

  return fcs_fabs(det_3x3(v0, v1, v2));
}

static fcs_int box_differs(
    const fcs_float *a1, const fcs_float *b1, const fcs_float *c1,
    const fcs_float *a2, const fcs_float *b2, const fcs_float *c2
    )
{
  for(fcs_int t=0; t<3; t++){
    if(a1[t] != a2[t]) return 1;
    if(b1[t] != b2[t]) return 1;
    if(c1[t] != c2[t]) return 1;
  }
  return 0;
}


static fcs_float sprod_3d(
    const fcs_float *x, const fcs_float *y 
    )
{
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

static int project_vec2_on_vec1(
    const fcs_float *vec1, const fcs_float* vec2
    )
{
  return round(sprod_3d(vec1, vec2) / sprod_3d(vec1, vec1));
}

static void orthogonalize_two_vectors(
    const fcs_float *vec1, fcs_float *vec2
    )
{
  const int proj = project_vec2_on_vec1(vec1, vec2);
  for(int t=0; t<3; ++t){
    vec2[t] -= proj * vec1[t]; 
  }
}

/* shift along periodic dims to maximize angles between box vectors */
static void orthogonalize_periodic_dims(
    const fcs_int *periodicity,
    fcs_float *a, fcs_float *b, fcs_float *c
    )
{
  /* make b and c orthogonal to periodic a */
  if(periodicity[0] == 1){
    orthogonalize_two_vectors(a, b);
    orthogonalize_two_vectors(a, c);
  }

  /* make a and c orthogonal to periodic b */
  if(periodicity[1] == 1){
    orthogonalize_two_vectors(b, a);
    orthogonalize_two_vectors(b, c);
  }

  /* make b and a orthogonal to periodic c */
  if(periodicity[2] == 1){
    orthogonalize_two_vectors(c, b);
    orthogonalize_two_vectors(c, a);
  }
}

FCSResult ifcs_p2nfft_tune(
    void *rd, const fcs_int *periodicity,
    fcs_int local_particles,
    const fcs_float *positions, const fcs_float *charges,
    const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c,
    const fcs_float *offset,
    fcs_int short_range_flag
    )
{
  int comm_rank;
  const char* fnc_name = "ifcs_p2nfft_tune";
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;
  fcs_int local_needs_retune = d->needs_retune;
  fcs_int i, num_particles;
  fcs_float sum_q, sum_q2, sum_q4, sum_q_abs, avg_dist, error=0;
  FCSResult result;

#if FCS_P2NFFT_DEBUG_RETUNE
  if(local_needs_retune)
    fprintf(stderr, "\n!!! begin checks: local_needs_retune = %d !!!\n\n", local_needs_retune);
#endif

  /* Check consistency of user defined near field radius */
  if( !d->tune_epsI && !d->tune_r_cut )
    return fcs_result_create(FCS_ERROR_LOGICAL_ERROR, fnc_name, "Nearfield radius was set in scaled ('epsI') and non-scaled ('r_cut') version. Choose one of them.");

  MPI_Comm_rank(d->cart_comm_3d, &comm_rank);
  FCS_P2NFFT_INIT_TIMING(d->cart_comm_3d);
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);

  d->short_range_flag = short_range_flag;

  /* shift box vectors along periodic dims in order to make box as orthogonal as possible */
  fcs_float box_a_ortho[3] = {box_a[0], box_a[1], box_a[2]};
  fcs_float box_b_ortho[3] = {box_b[0], box_b[1], box_b[2]};
  fcs_float box_c_ortho[3] = {box_c[0], box_c[1], box_c[2]};

  if(d->flags & FCS_P2NFFT_MAX_BOX_ANGLES){
    orthogonalize_periodic_dims(periodicity, box_a_ortho, box_b_ortho, box_c_ortho);
  }

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG  || FCS_P2NFFT_DEBUG_TUNING || FCS_ENABLE_INFO
  printf("input box vectors:   [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f], [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f], [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f]\n",
    box_a[0], box_a[1], box_a[2],
    box_b[0], box_b[1], box_b[2],
    box_c[0], box_c[1], box_c[2]);
  printf("shifted box vectors: [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f], [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f], [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f]\n",
    box_a_ortho[0], box_a_ortho[1], box_a_ortho[2],
    box_b_ortho[0], box_b_ortho[1], box_b_ortho[2],
    box_c_ortho[0], box_c_ortho[1], box_c_ortho[2]);
#endif

  /* Retune if simulation box changed */
  if(box_differs(d->box_a, d->box_b, d->box_c, box_a_ortho, box_b_ortho, box_c_ortho)){
#if FCS_P2NFFT_DEBUG_RETUNE
    fprintf(stderr, "Retune triggered due to changed box shape!\n");
#endif
    local_needs_retune = 1;
  }
  
  for(fcs_int t=0; t<3; t++){
    d->box_a[t] = box_a_ortho[t];
    d->box_b[t] = box_b_ortho[t];
    d->box_c[t] = box_c_ortho[t];
    d->box_base[t] = offset[t];
  }

  /* Get box size for orthorombic boxes */ 
  d->box_l[0] = fcs_norm(d->box_a);
  d->box_l[1] = fcs_norm(d->box_b);
  d->box_l[2] = fcs_norm(d->box_c);

  /* compute dual lattice vectors */
  d->box_V = fcs_fabs(det_3x3(d->box_a, d->box_b, d->box_c));
  invert_3x3(d->box_a, d->box_b, d->box_c, d->box_inv);

  /* Retune if total number of particles changed */
  MPI_Allreduce(&local_particles, &num_particles, 1, FCS_MPI_INT, MPI_SUM, d->cart_comm_3d);
  if (num_particles != d->num_nodes) {
#if FCS_P2NFFT_DEBUG_RETUNE
    fprintf(stderr, "Retune triggered due to changed num_nodes!, num_particles = %d, d->num_nodes = %d\n", num_particles, d->num_nodes);
#endif
    d->num_nodes = num_particles;
    local_needs_retune = 1;
  }

  /* Retune if periodicity changed */
  for(int t=0; t<3; t++){
    if((periodicity[t] == 0) != (d->periodicity[t] == 0)){
#if FCS_P2NFFT_DEBUG_RETUNE
      fprintf(stderr, "Retune triggered due to changed periodicity!, periodicity[%d] = %d, d->periodicity[%d] = %d\n", periodicity[t], t, d->periodicity[t], t);
#endif
      local_needs_retune = 1;
    }
    d->periodicity[t] = periodicity[t];
  }
  d->num_nonperiodic_dims = (periodicity[0]==0) + (periodicity[1]==0) + (periodicity[2]==0);
  d->num_periodic_dims    = (periodicity[0]!=0) + (periodicity[1]!=0) + (periodicity[2]!=0);
  
  /* 1d-periodic case */
  /* check incompl. Bessel function: display tables 1 and 2 from [Slavinsky, Safouhi 2010] */
//   if (d->num_periodic_dims==1)
//     ifcs_p2nfft_plot_slavinsky_safouhi_table();

  fcs_int reg_near=-1, reg_far=-1;
  if(d->reg_near == FCS_P2NFFT_REG_NEAR_DEFAULT)
    reg_near = FCS_P2NFFT_REG_NEAR_CG;
  else
    reg_near = d->reg_near;

  if(d->reg_far == FCS_P2NFFT_REG_FAR_RAD_CG)
    if(d->num_periodic_dims != 0)
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "FCS_P2NFFT_REG_FAR_RAD_CG is only available for 0d-periodicity.");

  if( (d->reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_SYM)
      || (d->reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_EC)
      || (d->reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_IC) )
  {
    if(d->num_periodic_dims != 0)
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Rectangular regularization is only available for 0d-periodicity.");
  }

  if(d->reg_far == FCS_P2NFFT_REG_FAR_DEFAULT){
    if(d->num_periodic_dims == 2)
      reg_far = FCS_P2NFFT_REG_FAR_RAD_T2P_SYM;
    else if(d->num_periodic_dims == 1)
      reg_far = FCS_P2NFFT_REG_FAR_RAD_T2P_IC;
    else if(d->num_periodic_dims == 0)
      reg_far = FCS_P2NFFT_REG_FAR_RAD_T2P_IC;
  } else
    reg_far = d->reg_far;

  if(d->reg_kernel == FCS_P2NFFT_REG_KERNEL_ONE_OVER_ABS_X)
    if(d->num_periodic_dims != 0)
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "FCS_P2NFFT_REG_KERNEL_ONBE_OVER_ABS_X is only available for 0d-periodicity.");


  /* Now, after the periodicity is clear, we can set the default tolerance type. */
  default_tolerance_type(d->periodicity, d->reg_kernel,
      &d->tolerance_type, &d->tolerance);

  /* Check if P2NFFT can handle the tolerance type */
  result = check_tolerance(d->periodicity, d->reg_kernel, d->tolerance_type, d->tolerance);
  if(result != NULL) return result;

  /* Calculate the sum of all charges
   * (pbc need neutral system) */
  fcs_float local_sum_q = 0;
  for (i = 0; i < local_particles; ++i)
    local_sum_q += charges[i];
  MPI_Allreduce(&local_sum_q, &sum_q, 1, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);
//   if (!fcs_float_is_equal(sum_q, d->sum_q)) {
  if (fcs_fabs(sum_q - d->sum_q) > 1e-5) {
#if FCS_P2NFFT_DEBUG_RETUNE
    fprintf(stderr, "Retune triggered due to changed charge sum!, sum_q = %e, d->sum_q = %e\n", sum_q, d->sum_q);
#endif
    d->sum_q = sum_q;
    local_needs_retune = 1;
  }
  
  if (!fcs_float_is_zero(d->sum_q))
    d->bg_charge = d->sum_q / d->num_nodes;
  else
    d->bg_charge = 0.0;

  /* Calculate the sum of the squares of all charges
   * (needed for the error formulae) */
  fcs_float local_sum_q2 = 0;
  for (i = 0; i < local_particles; ++i)
    local_sum_q2 += FCS_P2NFFT_SQR(charges[i] - d->bg_charge);
  MPI_Allreduce(&local_sum_q2, &sum_q2, 1, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);
//   if (!fcs_float_is_equal(sum_q2, d->sum_q2)) {
  if (fcs_fabs(sum_q2 - d->sum_q2 > 1e-5)) {
#if FCS_P2NFFT_DEBUG_RETUNE
    fprintf(stderr, "Retune triggered due to changed square charge sum!, sum_q2 = %e, d->sum_q2 = %e\n", sum_q2, d->sum_q2);
#endif
    d->sum_q2 = sum_q2;
    local_needs_retune = 1;
  }

  fcs_float local_sum_q4 = 0;
  for (i = 0; i < local_particles; ++i)
    local_sum_q4 += fcs_pow(charges[i] - d->bg_charge, 4.0);
  MPI_Allreduce(&local_sum_q4, &sum_q4, 1, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);

  /* Calculate the sum of the absolute values of all charges
   * (needed for the error formulae) */
  fcs_float local_sum_q_abs = 0;
  for (i = 0; i < local_particles; ++i)
    local_sum_q_abs += fabs(charges[i] - d->bg_charge);
  MPI_Allreduce(&local_sum_q_abs, &sum_q_abs, 1, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);
//   if (!fcs_float_is_equal(sum_q_abs, d->sum_q_abs)) {
//     d->sum_q_abs = sum_q_abs;
//     local_needs_retune = 1;
//   }

  /* FIXME: number of charged particles may be less than number of all particles */
  d->sum_qpart = num_particles;

  /* compute the average distance between two charges 
   * (needed for computation of default r_cut) */
  avg_dist = fcs_pow(d->box_V/d->sum_qpart, 0.33333);

#if FCS_P2NFFT_DEBUG_RETUNE
  if(local_needs_retune)
    fprintf(stderr, "\n!!! end checks: local_needs_retune = %d !!!\n\n", local_needs_retune);
#endif

  /* synchronize retuning on all processes */
  MPI_Allreduce(&local_needs_retune, &d->needs_retune, 1, FCS_MPI_INT, MPI_MAX, d->cart_comm_3d);

  if (d->needs_retune) {

    /* determine the dimensions with minimum and maximum box length */
    fcs_int maxdim=0, mindim=0;
    for(fcs_int t=1; t<2; t++){
      if(d->box_l[t] < d->box_l[mindim]) mindim = t;
      if(d->box_l[t] > d->box_l[maxdim]) maxdim = t;
    }

    /* check user defined epsI and epsB */
    if(d->num_nonperiodic_dims > 0)
      if(!d->tune_epsI && !d->tune_epsB)
        if(d->epsI + d->epsB >= 0.5)
          return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Sum of epsI and epsB must be less than 0.5.");

    /* At this point we choose default value 7. */
    if(d->tune_p)
      d->p = 8;

    if((d->num_periodic_dims > 0) || ((d->num_periodic_dims == 0) && (d->reg_kernel == FCS_P2NFFT_REG_KERNEL_EWALD)) ){
      fcs_float ks_error, rs_error;
      /* PNFFT calculates with real space cutoff 2*m+2
       * Therefore m is one less than the P3M cao. */
//      fcs_int cao = d->m + 1;
      fcs_int cao = 2*d->m;

      /* look for a dimension with pbc for tuning alpha and rcut */
      fcs_int dim_tune = get_dim_of_smallest_periodic_box_l(d->periodicity, d->box_l);
  
      if(d->tune_r_cut){
        /* set r_cut to 3 times the average distance between charges */
        d->r_cut = 3.0 * avg_dist;

        /* fulfill minimum image convention for periodic dims */
        for(int t=0; t<3; t++)
          if(periodicity[t])
            if (d->r_cut > d->box_l[t])
              d->r_cut = d->box_l[t];
      }
      d->one_over_r_cut = 1.0/d->r_cut;

      /* set normalized near field radius, relative to minimum box length */
      d->epsI = d->r_cut / d->box_l[mindim];
      
      /* Tune alpha for fixed N and m. */
      if(!d->tune_N && !d->tune_m){
        if(d->tune_alpha)
          d->alpha = p2nfft_tune_alpha(
              d->sum_qpart, d->sum_q2, dim_tune, d->box_l, d->r_cut, d->N, cao, d->tolerance, d->pnfft_direct);
        
        /* User specified N and cao. Therefore we do not necessarily need the error calculation. */
        if(~d->flags & FCS_P2NFFT_IGNORE_TOLERANCE){
          error = p2nfft_get_accuracy(
              d->sum_qpart, d->sum_q2, dim_tune, d->box_l, d->r_cut, d->N, cao,
              d->tolerance, d->alpha, d->pnfft_flags & PNFFT_INTERLACED, d->pnfft_direct,
              &rs_error, &ks_error);

          /* Return error, if tuning of alpha failed. */
          if(error > d->tolerance)
            return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Not able to reach required accuracy (Did not find feasible Ewald splitting parameter 'alpha').");
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
          for (i = FCS_P2NFFT_MINGRID; i <= FCS_P2NFFT_MAXGRID; /* see below */ ) {
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG  || FCS_P2NFFT_DEBUG_TUNING
            if(!comm_rank) printf("P2NFFT_INFO: r_cut = %" FCS_LMOD_FLOAT "f.\n", d->r_cut);
            if(!comm_rank) printf("P2NFFT_INFO: Trying grid size %" FCS_LMOD_INT "d.\n", i);
#endif
            fcs_int t0 = dim_tune, t1 = (t0+1)%3, t2 = (t0+2)%3;
            d->N[t0] = i;
            d->N[t1] = 2*fcs_ceil(d->N[t0]*d->box_l[t1]/(d->box_l[t0]*2.0));
            d->N[t2] = 2*fcs_ceil(d->N[t0]*d->box_l[t2]/(d->box_l[t0]*2.0));
            if(d->tune_alpha)
              d->alpha = p2nfft_tune_alpha(
                  d->sum_qpart, d->sum_q2, dim_tune, d->box_l, d->r_cut, d->N, cao, d->tolerance, d->pnfft_direct);
            error = p2nfft_get_accuracy(
                d->sum_qpart, d->sum_q2, dim_tune, d->box_l, d->r_cut, d->N, cao,
                d->tolerance, d->alpha, d->pnfft_flags & PNFFT_INTERLACED, d->pnfft_direct,
                &rs_error, &ks_error);
#if FCS_P2NFFT_DEBUG_TUNING
            if(!comm_rank) printf("P2NFFT_DEBUG_TUNING: error = %e, rs_error = %e, ks_error = %e, tolerance = %e.\n",
                error, rs_error, ks_error, d->tolerance);
#endif
            if (error < d->tolerance) break;

            if(d->pnfft_direct)
              i += 2;
            else
              i *= 2;
          }

          /* Return error, if tuning of N failed. */
          if(~d->flags & FCS_P2NFFT_IGNORE_TOLERANCE)
            if(error > d->tolerance)
              return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Not able to reach required accuracy (Did not find feasible gridsize 'N').");
        }

        /* tune m, only if we found a feasible gridsize N */
        if(d->tune_m && !d->pnfft_direct){
          for (cao = 1; cao <= FCS_P2NFFT_MAXCAO; ++cao) {
            if(d->tune_alpha)
              d->alpha = p2nfft_tune_alpha(
                  d->sum_qpart, d->sum_q2, dim_tune, d->box_l, d->r_cut, d->N, cao, d->tolerance, d->pnfft_direct);
            error = p2nfft_get_accuracy(
                d->sum_qpart, d->sum_q2, dim_tune, d->box_l, d->r_cut, d->N, cao,
                d->tolerance, d->alpha, d->pnfft_flags & PNFFT_INTERLACED, d->pnfft_direct,
                &rs_error, &ks_error);
            if (error < d->tolerance) break;
          }

          /* Return error, if tuning of m failed. */
          if(~d->flags & FCS_P2NFFT_IGNORE_TOLERANCE)
            if(error > d->tolerance)
              return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Not able to reach required accuracy (Did not find feasible charge assignment order 'm').");
        }
      }
     
      /* Return error, if accuracy tuning failed. */
      if(~d->flags & FCS_P2NFFT_IGNORE_TOLERANCE)
        if(error > d->tolerance)
          return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Not able to reach required accuracy.");

#if FCS_P2NFFT_TEST_GENERAL_ERROR_ESTIMATE
      if(!comm_rank)
        printf("\n\nP2NFFT_INFO: P3M Tuning results in N = %td, cao = %d, alpha = %e\n", d->N[0], cao, d->alpha);

      if(!comm_rank)
        printf("P2NFFT_INFO: P3M: k space error: %e\n", ks_error);

      if(!comm_rank)
        printf("P2NFFT_INFO: P2NFFT: k space error via new error estimate: %e, m = %d\n", 
            p2nfft_k_space_error_general_window(d->sum_qpart, d->sum_q2, d->box_l, d->N, d->alpha, (cao+1)/2, d->pnfft_window_flag, d->cart_comm_3d), (cao+1)/2);
#endif

      /* P3M tuning works with cao==2*m */
      d->m = (cao+1)/2;

      /* default oversampling equals 1 in 3d-periodic case */
      if(d->tune_n){
        for(int t=0; t<3; t++){
          if(d->periodicity[t])
            d->n[t] = d->N[t];
          else
            d->n[t] = 2.0*d->N[t];
        }
      } else {
        /* check user defined oversampling */
        for(int t=0; t<3; t++)
          if(d->N[t] > d->n[t] )
            return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Oversampled gridsize is not allowed to be less than normal gridsize.");
      }

      if(d->tune_epsB){
        /* only one dimension has npbc, choose this one to set epsB */
        for(int t=0; t<3; t++)
          if(!periodicity[t])
            d->epsB = (fcs_float)d->p/d->n[t];
        if(d->epsB > 0.125)
          d->epsB = 0.125;
      }
      
      /* Compute diagonal of box in non periodic dims. Works for orthogonal box vectors, only. */
      fcs_float diagonal = 0.0;
      for(int t=0; t<3; t++){
        if(!d->periodicity[t]) diagonal += d->box_l[t]*d->box_l[t];
      }
      diagonal = fcs_sqrt( diagonal ); /* for cubic boxes this is the same as sqrt(num_nonperiodic_dims)*L */

      for(int t=0; t<3; t++)
      {
        /* shift and scale coordinates into [-0.5,0.5) */
        if(d->periodicity[t]) d->box_expand[t] = 1.0;
        /* shift and scale coordinates into sphere with radius (0.5-epsB) */
        else {
//           d->box_expand[t] = 1.0 / (0.5 - d->epsB) * ( reg_far_is_radial(reg_far) ? fcs_sqrt(d->num_nonperiodic_dims) : 1.0 );
          if(reg_far_is_radial(reg_far)) d->box_expand[t] = 1.0 / (0.5 - d->epsB) * diagonal / d->box_l[t] ;
          else d->box_expand[t] = 1.0 / (0.5 - d->epsB) ;
        }

        /* TODO: remove deprecated variable box_scales */
        d->box_scales[t] = d->box_l[t] * d->box_expand[t]; /* for radial regularization this is equal to 1.0 / (0.5 - d->epsB) * diagonal */

        /* use full torus for periodic boundary conditions, otherwise use appropriate scaling */
        d->x_max[t] = 0.5 / d->box_expand[t];
      }
      
#if FCS_ENABLE_INFO
      if(!comm_rank)
        printf("P2NFFT_INFO: Approximate real space error: %" FCS_LMOD_FLOAT "e, k space error: %" FCS_LMOD_FLOAT "e\n",
            rs_error, ks_error);
      if(!comm_rank)
        printf("P2NFFT_INFO: Tuned alpha: %" FCS_LMOD_FLOAT "f, leading to an error of %" FCS_LMOD_FLOAT "e.\n",
            d->alpha, error);
      if(!comm_rank)
        printf("P2NFFT_INFO: Tuned N: %td, n = [%td, %td, %td], tuned m: %" FCS_LMOD_INT "d, r_cut: %" FCS_LMOD_FLOAT "f, epsI: %" FCS_LMOD_FLOAT "f, epsB: %" FCS_LMOD_FLOAT "f.\n",
            d->N[0], d->n[0], d->n[1], d->n[2], d->m, d->r_cut, d->epsI, d->epsB);
#endif

      /* Initialize the tables for near field interpolation */
      /*   accuracy of 1e-16 needs 14000 interpolation nodes */
      /*   accuracy of 1e-17 needs 24896 interpolation nodes */
//       d->near_interpolation_num_nodes = calc_interpolation_num_nodes(d->interpolation_order, 1e-16);
      unsigned err=0;
      d->near_interpolation_num_nodes = calc_interpolation_num_nodes_erf(d->interpolation_order, 0.1*d->tolerance, d->alpha, d->r_cut, &err);
      if(err)
        return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Number of nodes needed for interpolation exeedes 1e7. Try to use a higher interpolation order or less accuracy.");

      if(d->near_interpolation_table_potential != NULL){
        free(d->near_interpolation_table_potential);
        d->near_interpolation_table_potential = NULL;
      }
      if(d->near_interpolation_table_force != NULL){
        free(d->near_interpolation_table_force);
        d->near_interpolation_table_force = NULL;
      }
      if(d->far_interpolation_table_potential != NULL){
        free(d->far_interpolation_table_potential);
        d->far_interpolation_table_potential = NULL;
      }

      if(d->near_interpolation_num_nodes){
        d->near_interpolation_table_potential = (fcs_float*) malloc(sizeof(fcs_float) * (d->near_interpolation_num_nodes+3));
        init_near_interpolation_table_potential_3dp(
            d->near_interpolation_num_nodes,
            d->r_cut, d->alpha, 
            d->near_interpolation_table_potential);

        d->near_interpolation_table_force = (fcs_float*) malloc(sizeof(fcs_float) * (d->near_interpolation_num_nodes+3));
        init_near_interpolation_table_force_3dp(
            d->near_interpolation_num_nodes,
            d->r_cut, d->alpha,
            d->near_interpolation_table_force);
      }

      /* far field interpolation only works for equal box lengths in all non-periodic dimensions,
       * and only for radial far field regularization */
      d->far_interpolation_num_nodes = 0;
      if( reg_far_is_radial(reg_far) ){
        fcs_int num_intpol = d->near_interpolation_num_nodes;
        fcs_float epsB = d->epsB;
        fcs_float N_total=1.0;
        for(int t=0; t<3; ++t)
          if(!d->periodicity[t])
            N_total *= d->N[t];
        
        /* only use interpolation if it reduces the total number of function evaluations */
        fcs_int num_grid_reg=-1;
        if(d->num_periodic_dims==0)
          num_grid_reg = N_total * (1.0 - 4.0/3.0 * FCS_PI * fcs_pow(0.5-epsB, 3));
        else if(d->num_periodic_dims==1)
          num_grid_reg = N_total * (1.0 - FCS_PI * fcs_pow(0.5-epsB, 2));
        
        if(num_intpol < num_grid_reg)
          d->far_interpolation_num_nodes = 2 * d->near_interpolation_num_nodes;
        FCS_P2NFFT_IFDBG(if(comm_rank==0) fprintf(stderr, "num_intpol = %d, num_grid_reg = %d\n", num_intpol, num_grid_reg));
      }

      if(d->far_interpolation_num_nodes > 0){
        if( (d->num_periodic_dims==0) ){
          d->far_interpolation_table_potential = (fcs_float*) malloc(sizeof(fcs_float) * (d->far_interpolation_num_nodes+4));
          FCS_P2NFFT_IFDBG(double timer=-MPI_Wtime());
          init_far_interpolation_table_potential_0dp_ewald(
              d->far_interpolation_num_nodes, reg_far,
              d->alpha, d->box_scales[0], d->epsB, d->p, d->c,
              d->far_interpolation_table_potential);
          FCS_P2NFFT_IFDBG(timer+=MPI_Wtime(); if(!comm_rank) fprintf(stderr, "init 0dp far interpolation table takes %e s\n", timer));
        }

        if (d->num_periodic_dims==1){
          ptrdiff_t howmany = 1;
          ptrdiff_t local_Ni[3], local_Ni_start[3], local_No[3], local_No_start[3];
      
          /* 1dp: interpolation tables depend on local_Ni and local_No */
          FCS_PFFT(local_size_many_dft)(3, d->N, d->N, d->N, howmany,
              PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, d->cart_comm_pnfft, PFFT_TRANSPOSED_OUT| PFFT_SHIFTED_IN| PFFT_SHIFTED_OUT,
              local_Ni, local_Ni_start, local_No, local_No_start);
      
          fcs_int local_Ni_total=1;
          int pdim = -1, npdim = -1;
      
          for(int t=0; t<3; ++t){
            if(d->periodicity[t])  pdim = t;
            else                   npdim = t;
            if(!d->periodicity[t]) local_Ni_total *= local_Ni[t];
          }
      
          fcs_int offset = (d->far_interpolation_num_nodes+4);
          fcs_int mem = offset * local_Ni_total;
          d->far_interpolation_table_potential = (mem>0) ? (fcs_float*) malloc(sizeof(fcs_float) * mem) : NULL;
      
          if(local_Ni_total > 0){
            FCS_P2NFFT_IFDBG(double timer=-MPI_Wtime());
            fcs_int m=0;
            for(fcs_int k=local_Ni_start[pdim]; k<local_Ni_start[pdim]+local_Ni[pdim]; ++k, ++m)
              init_far_interpolation_table_potential_1dp(
                  d->far_interpolation_num_nodes, reg_far, k, d->box_l[pdim],
                  d->alpha, d->box_scales[npdim], d->epsB, d->p, d->c,
                  d->far_interpolation_table_potential + m * offset);
            FCS_P2NFFT_IFDBG(timer+=MPI_Wtime(); if(!comm_rank) fprintf(stderr, "init 1dp far interpolation tables takes %e s\n", timer));
          }
        }
      }
    } else {
      /********************/
      /* nonperiodic case */
      /********************/

      if(reg_far == FCS_P2NFFT_REG_FAR_RAD_CG)
        if(reg_near != FCS_P2NFFT_REG_NEAR_CG)
          return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Far field regularization FCS_P2NFFT_REG_FAR_RAD_CG is only available in combiniation with FCS_P2NFFT_REG_NEAR_CG.");

      if(!is_cubic(d->box_l))
        if(   reg_far != FCS_P2NFFT_REG_FAR_RAD_T2P_EC
           && reg_far != FCS_P2NFFT_REG_FAR_REC_T2P_SYM
           && reg_far != FCS_P2NFFT_REG_FAR_REC_T2P_IC
           && reg_far != FCS_P2NFFT_REG_FAR_REC_T2P_EC
          )
          return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Noncubic boxes require far field regularization FCS_P2NFFT_REG_FAR_RAD_T2P_EC, FCS_P2NFFT_REG_FAR_REC_T2P_SYM, FCS_P2NFFT_REG_FAR_REC_T2P_IC, or FCS_P2NFFT_REG_FAR_REC_T2P_EC.");

      if(reg_near == FCS_P2NFFT_REG_NEAR_T2P){
        /* TODO: implement parameter tuning for 2-point-Taylor regularization 
         * Question: Do we need it? Afaik, CG-approximation is better for all cases. */
  
        /*****************************************************************/
        /* example result of 2-point-Taylor parameter estimator:         */
        /* % eps_I = 0.062500;  N = 64;                                  */
        /* error = 2.782434e-03                                          */      
        /*****************************************************************/
        if (d->tune_N){
          d->N[mindim] = 64;
          /* Set N for all dimensions - suitable for noncubic geometry */
          for (fcs_int t = 0; t < 3; ++t)
            d->N[t] = 2*fcs_ceil(d->N[mindim]*d->box_l[t]/(d->box_l[mindim]*2.0));
            // This construction guarantees that N[t] is even and rounded up
        }

        if(d->tune_m) d->m = 4;
        if(d->tune_p) d->p = 7;
  
        if(d->tune_epsI){
          d->log2epsI = 5;  //  1.0/16 == 4.0/64
          d->epsI = fcs_pow(0.5, d->log2epsI);
        }
  
        if(d->tune_epsB){
          d->log2epsB = d->log2epsI;
          d->epsB = d->epsI;
        }

        /* shift and scale box with boxlength L/2 into 3d-ball with radius (0.25-epsB/2) */
        for(int t=0; t<3; t++){
          d->box_expand[t] = 1.0 / (0.5 - d->epsB) * ( reg_far_is_radial(reg_far) ? fcs_sqrt(3) : 1.0 );
          d->box_scales[t] = d->box_l[t] * d->box_expand[t];
        }
  
        /* initialize coefficients of 2-point Taylor polynomials */
        if(d->taylor2p_coeff != NULL){
          free(d->taylor2p_coeff);
          d->taylor2p_coeff = NULL;
        }
        d->taylor2p_coeff = (fcs_float*) malloc(sizeof(fcs_float)*(d->p));
        ifcs_p2nfft_load_taylor2p_coefficients(d->p, d->taylor2p_coeff);
  
        if(d->taylor2p_derive_coeff != NULL){
          free(d->taylor2p_derive_coeff);
          d->taylor2p_derive_coeff = NULL;
        }
        d->taylor2p_derive_coeff = (fcs_float*) malloc(sizeof(fcs_float)*(d->p-1));
        ifcs_p2nfft_load_taylor2p_derive_coefficients(d->p, d->taylor2p_derive_coeff);

      } else if(reg_near == FCS_P2NFFT_REG_NEAR_CG) {
        /*********************/
        /* CG regularization */
        /*********************/

        /* Update these values, if new CG-coefficients have been added */
        fcs_int N_min=8, N_max=4096;
        fcs_int log2epsI_min=3, log2epsI_max=11; /* 1/8=0.125, ..., 1/128=0.0078125 */
       
        if(d->tune_epsI){
          if(d->tune_r_cut){
            /* set epsI to 2 times the average distance between charges */
            d->epsI = 2.0 * avg_dist/d->box_l[mindim];
  
            /* good choice for reqiured_accuracy == 1e3 with hammersley_ball_pos_1e4 */
//             d->epsI = 0.5 * avg_dist/d->box_l[mindim];
          } else { /* user defined r_cut, now scale it into unit cube */
            /* invert r_cut = box_scale * epsI, where box_scale = box_l*sqrt(3)/(0.5-epsB) depends on epsI (=epsB) */
            if(reg_far_is_radial(reg_far))
              d->epsI = 0.5 / (d->box_l[mindim] / d->r_cut * fcs_sqrt(3) + 1.0);
            else
              d->epsI = 0.5 / (d->box_l[mindim] / d->r_cut + 1.0);
          }
        }
  
        /* compute the biggest value of type 2^-log2epsI below epsI */
        d->log2epsI = lrint(floor( -log(d->epsI)/log(2.0) ));
  
        /* use maximum/minimum possible epsI */
        if(d->log2epsI < log2epsI_min)
          d->log2epsI = log2epsI_min;
        if(d->log2epsI > log2epsI_max)
          d->log2epsI = log2epsI_max;
  
        d->epsI = fcs_pow(0.5, d->log2epsI);
  
        /* CG-approximation always uses epsI == epsB */
        d->epsB = d->epsI;
        
        /* shift and scale box with boxlength L/2 into 3d-ball with radius (0.25-epsB/2) */
        for(int t=0; t<3; t++){
          d->box_expand[t] = 1.0 / (0.5 - d->epsB) * ( reg_far_is_radial(reg_far) ? fcs_sqrt(3) : 1.0 );
          d->box_scales[t] = d->box_l[t] * d->box_expand[t];
        }

        if(d->tune_N){
          fcs_int N, m, p;
          for(N=N_min; N<=N_max; N*=2){
            d->N_cg_cos = N/2;
            error = ifcs_p2nfft_get_cg_cos_err(d->N_cg_cos, d->log2epsI,
              &m, &p);

#if FCS_P2NFFT_DEBUG_TUNING
            if(!comm_rank){
              printf("P2NFFT_INFO: Trying grid size %d.\n", N);
              if(error < 0.0) 
                printf("P2NFFT_DEBUG_TUNING: No CG approximation available.\n");
              else
                printf("P2NFFT_DEBUG_TUNING: error = %e, sum_q_abs = %e, accuracy = %e, tolerance = %e.\n",
                    error, sum_q_abs, error * sum_q_abs * fcs_sqrt(2.0 / (d->box_scales[0]*d->box_scales[1]*d->box_scales[2])), d->tolerance);
            }
#endif

            /* TODO: find better error estimate */
            if( !(error < 0.0) )
              if (error * sum_q_abs * fcs_sqrt(2.0 / (d->box_scales[0]*d->box_scales[1]*d->box_scales[2])) < d->tolerance) break;
          }
    
          /* Return error, if accuracy tuning failed. */
          if(N==2*N_max)
            if(~d->flags & FCS_P2NFFT_IGNORE_TOLERANCE)
              return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Not able to reach required accuracy.");
  
          /* We found the minimal N. This N corresponds to the minimal box length */
          d->N[mindim] = N;

          /* Set N for all dimensions - suitable for noncubic geometry */
          for (fcs_int t = 0; t < 3; ++t)
            d->N[t] = 2*fcs_ceil(d->N[mindim]*d->box_l[t]/(d->box_l[mindim]*2.0));
            // This construction guarantees that N[t] is even and rounded up
        }
  
        /* initialize cg-regularization */
        d->N_cg_cos = d->N[mindim]/2;
        if(d->cg_cos_coeff != NULL){
          free(d->cg_cos_coeff);
          d->cg_cos_coeff = NULL;
        }
       
        fcs_int m, p; 
        d->cg_cos_coeff = (fcs_float*) malloc(sizeof(fcs_float)*d->N_cg_cos);
        int missed_coeff = 0;
        if(reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_SYM)
          missed_coeff = ifcs_p2nfft_load_cg_cos_coeff_sym(d->N_cg_cos, d->log2epsI,
              &m, &p, d->cg_cos_coeff);
        else
          missed_coeff = ifcs_p2nfft_load_cg_cos_coeff(d->N_cg_cos, d->log2epsI,
            &m, &p, d->cg_cos_coeff);

        if(missed_coeff)
          return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Did not find an appropriate CG approximation.");

        if(d->cg_sin_coeff != NULL){
          free(d->cg_sin_coeff);
          d->cg_sin_coeff = NULL;
        }
        d->cg_sin_coeff = (fcs_float*) malloc(sizeof(fcs_float)*d->N_cg_cos);
        for(fcs_int k=0; k<d->N_cg_cos; k++)
          d->cg_sin_coeff[k] = -2 * FCS_P2NFFT_PI * k * d->cg_cos_coeff[k];

        if(d->tune_m) d->m = m;
        if(d->tune_p) d->p = p;
      } else
        return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Unknown near field regularization flag. Choose either FCS_P2NFFT_REG_NEAR_T2P or FCS_P2NFFT_REG_NEAR_CG.");

      if (d->tune_c){
        d->c = 0.0; /* continue regularization with 0.0 */
//         d->c = 2.0 / d->box_scales[maxdim]; /* corresponds to scaled kernel function at 0.5 */
      }

      /* set unscaled near field radius, relative to minimum box length */
      d->r_cut = d->epsI * d->box_scales[mindim];
      d->one_over_r_cut = 1.0/d->r_cut;

      /* default oversampling equals 2 in nonperiodic case */
      if(d->tune_n){
        for(int t=0; t<3; t++)
          d->n[t] = 2.0 * d->N[t];
      } else {
        /* check user defined oversampling */
        for(int t=0; t<3; t++)
          if(d->N[t] > d->n[t] )
            return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Oversampled gridsize 'n' is not allowed to be less than normal gridsize 'N'.");
      }

      /* TODO: Optimize m for the field error. */
      /* Force calculation needs increased real space cutoff m in comparison to potential calculation */
//       d->m += 1;

      /* Initialize the tables for near field interpolation */
      /*   accuracy of 1e-16 needs 14000 interpolation nodes */
      /*   accuracy of 1e-17 needs 24896 interpolation nodes */
      d->near_interpolation_num_nodes = calc_interpolation_num_nodes(d->interpolation_order, 1e-16);
//       d->near_interpolation_num_nodes = calc_interpolation_num_nodes(d->interpolation_order, d->tolerance);

      if(d->near_interpolation_table_potential != NULL){
        free(d->near_interpolation_table_potential);
        d->near_interpolation_table_potential = NULL;
      }
      if(d->near_interpolation_table_force != NULL){
        free(d->near_interpolation_table_force);
        d->near_interpolation_table_force = NULL;
      }
      if(d->far_interpolation_table_potential != NULL){
        free(d->far_interpolation_table_potential);
        d->far_interpolation_table_potential = NULL;
      }

      if(d->near_interpolation_num_nodes > 0){
        d->near_interpolation_table_potential = (fcs_float*) malloc(sizeof(fcs_float) * (d->near_interpolation_num_nodes+3));
        init_near_interpolation_table_potential_0dp(
            d->near_interpolation_num_nodes,
            d->r_cut, d->epsI, d->p,
            d->taylor2p_coeff,
            d->N_cg_cos, d->cg_cos_coeff,
            d->near_interpolation_table_potential);

        d->near_interpolation_table_force = (fcs_float*) malloc(sizeof(fcs_float) * (d->near_interpolation_num_nodes+3));
        init_near_interpolation_table_force_0dp(
            d->near_interpolation_num_nodes,
            d->r_cut, d->epsI, d->p,
            d->taylor2p_derive_coeff,
            d->N_cg_cos, d->cg_cos_coeff, d->cg_sin_coeff,
            d->near_interpolation_table_force);
      }

      /* far field interpolation only works for cubic boxes and radial far field regularization */
      if( reg_far_is_radial(reg_far)  && is_cubic(d->box_l) ){
        d->far_interpolation_num_nodes = (int) d->near_interpolation_num_nodes * d->epsB / d->epsI;
        FCS_P2NFFT_IFDBG(if(!comm_rank) fprintf(stderr, "interpolation nodes: near = %d, far = %d\n", d->near_interpolation_num_nodes, d->far_interpolation_num_nodes));
      } else
        d->far_interpolation_num_nodes = 0;

      if(d->far_interpolation_num_nodes > 0){
        /* far field regularization table needs one extra point before 0.5-epsB for cubic interpolation */
        d->far_interpolation_table_potential = (fcs_float*) malloc(sizeof(fcs_float) * (d->far_interpolation_num_nodes+4));
        init_far_interpolation_table_potential_0dp(
            d->far_interpolation_num_nodes, reg_far,
            d->epsB, d->p, d->c,
            d->N_cg_cos, d->cg_cos_coeff,
            d->far_interpolation_table_potential);
      }


      for(int t=0; t<3; t++)
        d->x_max[t] = 0.5 / d->box_expand[t];
     
#if FCS_ENABLE_INFO 
      if(!comm_rank){
        printf("P2NFFT_INFO: Tuned parameters: N = [%td, %td, %td], n = [%td, %td, %td], m = %" FCS_LMOD_INT "d, p = %" FCS_LMOD_INT "d, rcut = %" FCS_LMOD_FLOAT "f, epsI = %" FCS_LMOD_FLOAT "f, near field interpolation nodes = %" FCS_LMOD_INT "d, far field interpolation nodes = %" FCS_LMOD_INT "d, cg_err_3d = %" FCS_LMOD_FLOAT "e, tolerance = %" FCS_LMOD_FLOAT "e\n",
            d->N[0], d->N[1], d->N[2], d->n[0], d->n[1], d->n[2], d->m, d->p, d->r_cut, d->epsI, d->near_interpolation_num_nodes, d->far_interpolation_num_nodes, error, d->tolerance);
        printf("P2NFFT_INFO: General error bound (best for non-alternating charges): cg_err_3d * sum_q_abs = %" FCS_LMOD_FLOAT "e\n",
            error*sum_q_abs);
        printf("P2NFFT_INFO: Stochastical error bound for alternating charges: cg_err_3d * sum_q_abs / fcs_sqrt(M) = %" FCS_LMOD_FLOAT "e\n",
            error*sum_q_abs/fcs_sqrt(d->num_nodes));
        printf("P2NFFT_INFO: Test of new Stochastical error bound for near field error: cg_err_3d * sum_q2 / fcs_sqrt(M*V) = %" FCS_LMOD_FLOAT "e\n",
            error * d->sum_q2/fcs_sqrt(d->num_nodes*d->box_l[0]*d->box_l[1]*d->box_l[2]));
        printf("P2NFFT_INFO: Test of new Stochastical error bound for near plus far field error: cg_err_3d * fcs_sqrt( sum_q2*sum_q2 / (M*V) + 1 ) = %" FCS_LMOD_FLOAT "e\n",
            error * fcs_sqrt(d->sum_q2*d->sum_q2/(d->num_nodes*d->box_l[0]*d->box_l[1]*d->box_l[2]) + 1));
        printf("P2NFFT_INFO: Test of new Stochastical error bound for near plus far field error: cg_err_3d * sum_q2 * fcs_sqrt(3 / (M*V)) = %" FCS_LMOD_FLOAT "e\n",
            error * d->sum_q2 * fcs_sqrt(3.0 / (d->num_nodes*d->box_l[0]*d->box_l[1]*d->box_l[2])) );
        printf("P2NFFT_INFO: General error bound (depending on box scale=%" FCS_LMOD_FLOAT "e): cg_err_3d * sum_q_abs * fcs_sqrt(2.0/V_scale) = %" FCS_LMOD_FLOAT "e\n",
            d->box_scales[0], error*sum_q_abs * fcs_sqrt(2.0 / (d->box_scales[0]*d->box_scales[1]*d->box_scales[2]) ));
      }
#endif
    }


    /* compute expanded box vectors */
    d->ebox_a[0] = d->box_a[0] * ( d->periodicity[0] ? 1.0 : d->box_expand[0] ); 
    d->ebox_a[1] = d->box_a[1] * ( d->periodicity[0] ? 1.0 : d->box_expand[0] ); 
    d->ebox_a[2] = d->box_a[2] * ( d->periodicity[0] ? 1.0 : d->box_expand[0] ); 

    d->ebox_b[0] = d->box_b[0] * ( d->periodicity[1] ? 1.0 : d->box_expand[1] ); 
    d->ebox_b[1] = d->box_b[1] * ( d->periodicity[1] ? 1.0 : d->box_expand[1] ); 
    d->ebox_b[2] = d->box_b[2] * ( d->periodicity[1] ? 1.0 : d->box_expand[1] ); 

    d->ebox_c[0] = d->box_c[0] * ( d->periodicity[2] ? 1.0 : d->box_expand[2] ); 
    d->ebox_c[1] = d->box_c[1] * ( d->periodicity[2] ? 1.0 : d->box_expand[2] ); 
    d->ebox_c[2] = d->box_c[2] * ( d->periodicity[2] ? 1.0 : d->box_expand[2] ); 

    /* compute volume and inverse box vectors for expanded box */
    d->ebox_V = fcs_fabs(det_3x3(d->ebox_a, d->ebox_b, d->ebox_c));
    invert_3x3(d->ebox_a, d->ebox_b, d->ebox_c, d->ebox_inv);

#if FCS_P2NFFT_EXIT_AFTER_TUNING
    return fcs_result_create(FCS_ERROR_LOGICAL_ERROR, fnc_name, "End of tuning.");
#endif

    /* calculate local data distribution according to PNFFT:
     * local_N, local_N_start, lower_border, upper_border */
    FCS_PNFFT(local_size_guru)(3, d->N, d->n, d->x_max, d->m, d->cart_comm_pnfft, d->pnfft_flags,
        d->local_N, d->local_N_start, d->lower_border, d->upper_border);

    /* shift decomposition of the torus from [-1/2,1/2)^3 to [0,1)^3 */
    for(int t=0; t<3; t++)
    {
      d->lower_border[t] += 0.5 / d->box_expand[t];
      d->upper_border[t] += 0.5 / d->box_expand[t];
    }
  }
  /* Finish timing of parameter tuning */
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Parameter tuning");

  /* Start timing of precomputation */
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);
  if (d->needs_retune) {
    if(d->regkern_hat != NULL){
      FCS_PFFT(free)(d->regkern_hat);
      d->regkern_hat = NULL;
    }

    /* precompute Fourier coefficients for convolution */
    if (d->num_periodic_dims == 3)
      d->regkern_hat = malloc_and_precompute_regkern_hat_3dp(
          d->local_N, d->local_N_start, d->box_inv, d->alpha, d->k_cut);
    else if ( (d->num_periodic_dims == 2) || (d->num_periodic_dims == 1) )
      d->regkern_hat = malloc_and_precompute_regkern_hat_2dp_and_1dp(
          d->N, d->epsB, d->box_a, d->box_b, d->box_c, d->box_inv, d->box_scales, d->alpha, d->k_cut, d->periodicity, d->p, d->c, reg_far,
          d->interpolation_order, d->far_interpolation_num_nodes, d->far_interpolation_table_potential,
          d->cart_comm_pnfft);
    else if (d->num_periodic_dims == 0) {
      if (d->reg_kernel == FCS_P2NFFT_REG_KERNEL_EWALD) {
        d->regkern_hat = malloc_and_precompute_regkern_hat_0dp_ewald(
            d->N, d->epsB, d->box_scales, d->alpha, d->p, d->c, reg_far,
            d->interpolation_order, d->far_interpolation_num_nodes, d->far_interpolation_table_potential,
            d->cart_comm_pnfft, is_cubic(d->box_l));
      } else if (d->reg_kernel == FCS_P2NFFT_REG_KERNEL_ONE_OVER_ABS_X) {
        d->regkern_hat = malloc_and_precompute_regkern_hat_0dp(
            d->N, d->r_cut, d->epsI, d->epsB, d->p, d->c, d->box_scales, reg_near, reg_far,
            d->taylor2p_coeff, d->N_cg_cos, d->cg_cos_coeff,
            d->interpolation_order, d->near_interpolation_num_nodes, d->far_interpolation_num_nodes,
            d->near_interpolation_table_potential, d->far_interpolation_table_potential,
            d->cart_comm_pnfft, is_cubic(d->box_l));
      }
    }
  }
  /* Finish timing of of precomputation */
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Precomputation of regularization");

  /* Initialize the plan for the PNFFT */
  init_pnfft(&d->pnfft, 3, d->N, d->n, d->x_max, d->m, d->pnfft_flags, d->pnfft_interpolation_order, d->pnfft_window,
      d->pfft_flags, d->pfft_patience, d->cart_comm_pnfft);

  if(d->tune_b){
    FCS_PNFFT(get_b)(d->pnfft, &d->b[0], &d->b[1], &d->b[2]);
  } else {
    /* skip retune if b is the same in P2NFFT and PNFFT plans
     * This is important, since we call fcs_tune as part of fcs_run. */
    fcs_float b[3] = {0,0,0};
    FCS_PNFFT(get_b)(d->pnfft, &b[0], &b[1], &b[2]);

    fcs_int b_is_out_of_sync=0;
    for(int t=0; t<3; t++)
      if(!fcs_float_is_zero(d->b[t]-b[t]))
        b_is_out_of_sync = 1;

    if(b_is_out_of_sync)
      FCS_PNFFT(set_b)(d->b[0], d->b[1], d->b[2], d->pnfft);
  }

  /* Print the command line arguments that recreate this plan. */
  if(d->needs_retune){
    if(d->flags & FCS_P2NFFT_VERBOSE_TUNING){
      if(!comm_rank) printf("P2NFFT_INFO: CMD ARGS: -c ");
      print_command_line_arguments(d, 0);
      if(!comm_rank) printf("P2NFFT_INFO: ALL CMD ARGS: -c ");
      print_command_line_arguments(d, 1);
      if(!comm_rank) printf("P2NFFT_INFO: Q^2 = %.16e\n", d->sum_q2);
      if(!comm_rank) printf("P2NFFT_INFO: Q^4 = %.16e\n", sum_q4);
    }
  }

  d->needs_retune = 0;

  /* Check if gridsize is too small to run with this number of processes. */
  int count=0;
  for(int t=0; t<3; t++)
    if(d->N[t] < d->np[t]) count++;

  if(count>0){
#if FCS_P2NFFT_DEBUG_RETUNE
    printf("P2NFFT_DEBUG: d->N = [%td, %td, %td], d->np = [%d, %d, %d]\n", d->N[0], d->N[1], d->N[2], d->np[0], d->np[1], d->np[2]);
#endif
    return fcs_result_create(FCS_ERROR_LOGICAL_ERROR, fnc_name, "Procmesh too large for this gridsize.");
  }

  return NULL;
}

/* For verbose == 0 only print the user defined parameters. 
 * For verbose != 0 print also the parameters which where determined by the tuning. */
static void print_command_line_arguments(
    ifcs_p2nfft_data_struct *d, fcs_int verbose
    )
{
  int comm_rank;
  MPI_Comm_rank(d->cart_comm_3d, &comm_rank);
 
  /* print full set of command line arguments */ 
  if(!comm_rank){
    if(verbose || !fcs_float_is_equal(d->tolerance, FCS_P2NFFT_DEFAULT_TOLERANCE)){
      switch(d->tolerance_type){
        case FCS_TOLERANCE_TYPE_ENERGY:
          printf("tolerance_energy,"); break;
        case FCS_TOLERANCE_TYPE_POTENTIAL:
          printf("tolerance_potential,"); break;
        case FCS_TOLERANCE_TYPE_FIELD:
          printf("tolerance_field,"); break;
        case FCS_TOLERANCE_TYPE_ENERGY_REL:
          printf("tolerance_energy_rel,"); break;
        case FCS_TOLERANCE_TYPE_POTENTIAL_REL:
          printf("tolerance_potential_rel,"); break;
        case FCS_TOLERANCE_TYPE_FIELD_REL:
          printf("tolerance_field_rel,"); break;
      }
      printf("%e,", d->tolerance);
    }

    /* print P2NFFT specific parameters */
    if(verbose || !d->tune_epsI || !d->tune_r_cut)
      printf("p2nfft_r_cut,%" FCS_LMOD_FLOAT "f,", d->r_cut);
    if(verbose || !d->tune_epsI || !d->tune_r_cut)
      printf("p2nfft_epsI,%" FCS_LMOD_FLOAT "f,", d->epsI);
    if(verbose || !d->tune_epsB)
      printf("p2nfft_epsB,%" FCS_LMOD_FLOAT "f,", d->epsB);
    if(verbose || !d->tune_k_cut)
      printf("p2nfft_k_cut,%" FCS_LMOD_FLOAT "f,", d->k_cut);
    if(verbose || !d->tune_c)
      printf("p2nfft_c,%" FCS_LMOD_FLOAT "f,", d->c);
    if(verbose || !d->tune_alpha)
      printf("p2nfft_alpha,%" FCS_LMOD_FLOAT "f,", d->alpha);
    if(verbose || d->interpolation_order != 3)
      printf("p2nfft_intpol_order,%" FCS_LMOD_INT "d,", d->interpolation_order);
    if(verbose || (d->reg_near != FCS_P2NFFT_REG_NEAR_DEFAULT) ){
      printf("p2nfft_reg_near_name,");
      if(d->reg_near == FCS_P2NFFT_REG_NEAR_DEFAULT)
        printf("default,");
      else if(d->reg_near == FCS_P2NFFT_REG_NEAR_CG)
        printf("cg,");
      else if(d->reg_near == FCS_P2NFFT_REG_NEAR_T2P)
        printf("t2p,");
    }
    if(verbose || (d->reg_far != FCS_P2NFFT_REG_FAR_DEFAULT) ){
      printf("p2nfft_reg_far_name,");
      if(d->reg_far == FCS_P2NFFT_REG_FAR_DEFAULT)
        printf("default,");
      else if(d->reg_far == FCS_P2NFFT_REG_FAR_RAD_CG)
        printf("rad_cg,");
      else if(d->reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM)
        printf("rad_t2p_sym,");
      else if(d->reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC)
        printf("rad_t2p_ec,");
      else if(d->reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC)
        printf("rad_t2p_ic,");
      else if(d->reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_SYM)
        printf("rec_t2p_sym,");
      else if(d->reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_EC)
        printf("rec_t2p_ec,");
      else if(d->reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_IC)
        printf("rec_t2p_ic,");
    }
    if(verbose || (d->reg_kernel != FCS_P2NFFT_REG_KERNEL_DEFAULT) ){
      printf("p2nfft_reg_kernel_name,");
      if(d->reg_kernel == FCS_P2NFFT_REG_KERNEL_EWALD)
        printf("ewald,");
      else if(d->reg_kernel == FCS_P2NFFT_REG_KERNEL_ONE_OVER_ABS_X)
        printf("one_over_abs_x,");
    }
    if(verbose || !d->tune_p)
      printf("p2nfft_p,%" FCS_LMOD_INT "d,", d->p);
    if(verbose || (d->flags & FCS_P2NFFT_IGNORE_TOLERANCE) )
      printf("p2nfft_ignore_tolerance,%d,", (d->flags & FCS_P2NFFT_IGNORE_TOLERANCE) ? 1 : 0);
    if(verbose || (d->flags & FCS_P2NFFT_IGNORE_POTENTIAL) )
      printf("p2nfft_ignore_potential,%d,", (d->flags & FCS_P2NFFT_IGNORE_POTENTIAL) ? 1 : 0);
    if(verbose || (d->flags & FCS_P2NFFT_IGNORE_FIELD) )
      printf("p2nfft_ignore_field,%d,", (d->flags & FCS_P2NFFT_IGNORE_FIELD) ? 1 : 0);
    if(verbose || (d->flags & FCS_P2NFFT_VERBOSE_TUNING) )
      printf("p2nfft_verbose_tuning,%d,", (d->flags & FCS_P2NFFT_VERBOSE_TUNING) ? 1 : 0);
    if(verbose || (d->flags & FCS_P2NFFT_MAX_BOX_ANGLES) )
      printf("p2nfft_max_box_angles,%d,", (d->flags & FCS_P2NFFT_MAX_BOX_ANGLES) ? 1 : 0);
    
    if(verbose || (d->virial != NULL) )
      printf("p2nfft_require_virial,%d,", (d->virial != NULL) ? 1 : 0);

    /* print PNFFT specific parameters */
    if(verbose || !d->tune_N)
      printf("pnfft_N,%td,%td,%td,", d->N[0], d->N[1], d->N[2]);
    if(verbose || !d->tune_n)
      printf("pnfft_n,%td,%td,%td,", d->n[0], d->n[1], d->n[2]);
    if(verbose || (d->pnfft_window != FCS_P2NFFT_DEFAULT_PNFFT_WINDOW) ){
      printf("pnfft_window_name,");
      switch(d->pnfft_window){
        case 0: printf("gaussian,"); break;
        case 1: printf("bspline,"); break;
        case 2: printf("sinc,"); break;
        case 3: printf("kaiser,"); break;
        case 4: printf("bessel_i0,"); break;
        case 5: printf("gaussian_t,"); break;
        default: printf("failure,");
      }
    }
    if(verbose || !d->tune_b)
      printf("pnfft_b,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,", d->b[0], d->b[1], d->b[2]);
    if(verbose || !d->tune_m)
      printf("pnfft_m,%" FCS_LMOD_INT "d,", d->m);
    if(verbose || d->pnfft_interpolation_order != 3)
      printf("pnfft_intpol_order,%" FCS_LMOD_INT "d,", d->pnfft_interpolation_order);
    if(verbose || !d->pnfft_direct)
      printf("pnfft_direct,%" FCS_LMOD_INT "d,", d->pnfft_direct);
    if(verbose || (d->pnfft_flags & PNFFT_PRE_PHI_HAT) )
      printf("pnfft_pre_phi_hat,%d,", (d->pnfft_flags & PNFFT_PRE_PHI_HAT) ? 1 : 0);
    if(verbose || (d->pnfft_flags & PNFFT_FFT_IN_PLACE) )
      printf("pnfft_fft_in_place,%d,", (d->pnfft_flags & PNFFT_FFT_IN_PLACE) ? 1 : 0);
    if(verbose || (d->pnfft_flags & PNFFT_SORT_NODES) )
      printf("pnfft_sort_nodes,%d,", (d->pnfft_flags & PNFFT_SORT_NODES) ? 1 : 0);
    if(verbose || (d->pnfft_flags & PNFFT_INTERLACED) )
      printf("pnfft_interlaced,%d,", (d->pnfft_flags & PNFFT_INTERLACED) ? 1 : 0);
    if(verbose || (d->pnfft_flags & PNFFT_REAL_F) )
      printf("pnfft_real_f,%d,", (d->pnfft_flags & PNFFT_REAL_F) ? 1 : 0);

    /* TODO: upadte to diff_ik */
//     if(d->compute_flags & PNFFT_GRAD_F) ...
//     if(d->pnfft_flags & PNFFT_DIFF_IK)
//       printf("pnfft_diff_ik,%d,", (d->pnfft_flags & PNFFT_DIFF_IK) ? 1 : 0);
//
//     if(d->pnfft_flags & PNFFT_GRAD_IK)
//       printf("pnfft_grad_ik,%d,", (d->pnfft_flags & PNFFT_GRAD_IK) ? 1 : 0);
//     else if(d->pnfft_flags & PNFFT_GRAD_NONE)
//       printf("pnfft_grad_none,%d,", (d->pnfft_flags & PNFFT_GRAD_NONE) ? 1 : 0);
//     else if(verbose)
//       printf("pnfft_grad_ik,0,");

    /* print PFFT specific parameters */
    if(verbose || (d->pfft_patience != FCS_P2NFFT_DEFAULT_PFFT_PATIENCE) ){
      printf("pfft_patience_name,");
      switch(d->pfft_patience){
        case 0 : printf("estimate,"); break;
        case 1 : printf("measure,"); break;
        case 2 : printf("patient,"); break;
        case 3 : printf("exhaustive,"); break;
        default: printf("failure,");
      }
    }
    if(verbose || (d->pfft_flags & PFFT_TUNE) )
      printf("pfft_tune,%d,", (d->pfft_flags & PFFT_TUNE) ? 1 : 0);
    if(verbose || (d->pfft_flags & PFFT_PRESERVE_INPUT) )
      printf("pfft_preserve_input,%d,", (d->pfft_flags & PFFT_PRESERVE_INPUT) ? 1 : 0);

    printf("\n");
  }
}

static void init_near_interpolation_table_potential_3dp(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float alpha,
    fcs_float *table
    )
{
  fcs_float r;

  for(fcs_int k=0; k<num_nodes+3; k++){
    r = r_cut * (fcs_float) k / num_nodes;
    if (fcs_float_is_zero(r))
      table[k] = 2 * alpha * FCS_P2NFFT_1_SQRTPI;
    else
      table[k] = erf(alpha * r)/r;
  }
}

static void init_near_interpolation_table_potential_0dp(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float epsI, fcs_int p,
    const fcs_float *taylor2p_coeff,
    fcs_int N_cos, const fcs_float *cos_coeff,
    fcs_float *table
    )
{
  if(cos_coeff != NULL){
    for(fcs_int k=0; k<num_nodes+3; k++)
      table[k] = evaluate_cos_polynomial_1d(
          epsI * (fcs_float) k / num_nodes, N_cos, cos_coeff) * epsI/r_cut;
  } else if(taylor2p_coeff != NULL) {
    for(fcs_int k=0; k<num_nodes+3; k++) {
      table[k] = ifcs_p2nfft_nearfield_correction_taylor2p(
          (fcs_float) k / num_nodes, p, taylor2p_coeff) / r_cut;
    }
  }
}

/* This function is only used for cubic boxes. */
static void init_far_interpolation_table_potential_0dp(
    fcs_int num_nodes, fcs_int reg_far,
    fcs_float epsB, fcs_int p, fcs_float c,
    fcs_int N_cos, const fcs_float *cos_coeff,
    fcs_float *table 
    )
{
  fcs_float epsI = 0.0;
  fcs_float *param = NULL;

  if(reg_far == FCS_P2NFFT_REG_FAR_RAD_CG){
    /* use CG approximiation */
    for(fcs_int k=-1; k<num_nodes+3; k++)
      table[k+1] = evaluate_cos_polynomial_1d(
          0.5 - epsB + epsB * (fcs_float) k / num_nodes, N_cos, cos_coeff);
  } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
    /* use symmetric Taylor2p, not smooth at r=0.5 in 3d */
    for(fcs_int k=-1; k<num_nodes+3; k++)
      table[k+1] = ifcs_p2nfft_reg_far_rad_sym(
          ifcs_p2nfft_one_over_modulus, param, 0.5 - epsB + epsB * (fcs_float) k / num_nodes, p, epsI, epsB);
  } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
    /* use normal basis polynomials and set constant continuation value 'c' explicitly */
    /* Since we only evaluate regkernel in [0.5-epsB, 0.5] we can savely set epsI = epsB */
    for(fcs_int k=-1; k<num_nodes+3; k++)
      table[k+1] = ifcs_p2nfft_reg_far_rad_ec(
          ifcs_p2nfft_one_over_modulus, param, 0.5 - epsB + epsB * (fcs_float) k / num_nodes, p, epsI, epsB, c);
  } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
    /* use integrated basis polynomials, sets continuation value 'c' implicitily */
    for(fcs_int k=-1; k<num_nodes+3; k++)
      table[k+1] = ifcs_p2nfft_reg_far_rad_ic(
          ifcs_p2nfft_one_over_modulus, param, 0.5 - epsB + epsB * (fcs_float) k / num_nodes, p, epsI, epsB);
  } else
    for(fcs_int k=-1; k<num_nodes+3; k++)
      table[k+1] = 0.0;
}

static void init_far_interpolation_table_potential_0dp_ewald(
    fcs_int num_nodes, fcs_int reg_far,
    fcs_float alpha, fcs_float h,
    fcs_float epsB, fcs_int p, fcs_float c,
    fcs_float *table 
    )
{
  fcs_float param[3] = {alpha, 0, h}; 

  if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
    /* use symmetric Taylor2p, not smooth at r=0.5*h in 3d */
    for(fcs_int k=-1; k<num_nodes+3; k++)
      table[k+1] = ifcs_p2nfft_reg_far_rad_sym_no_singularity(
          ifcs_p2nfft_erfx_over_x, param, h*(0.5 - epsB + epsB * (fcs_float) k / num_nodes), p, epsB);
  } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
    /* use normal basis polynomials and set constant continuation value 'c' explicitly */
    /* Since we only evaluate regkernel in h*[0.5-epsB, 0.5] we can savely set epsI = epsB */
    for(fcs_int k=-1; k<num_nodes+3; k++)
      table[k+1] = ifcs_p2nfft_reg_far_rad_ec_no_singularity(
          ifcs_p2nfft_erfx_over_x, param, h*(0.5 - epsB + epsB * (fcs_float) k / num_nodes), p, epsB, c);
  } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
    /* use integrated basis polynomials, sets continuation value 'c' implicitily */
    for(fcs_int k=-1; k<num_nodes+3; k++)
      table[k+1] = ifcs_p2nfft_reg_far_rad_ic_no_singularity(
          ifcs_p2nfft_erfx_over_x, param, h*(0.5 - epsB + epsB * (fcs_float) k / num_nodes), p, epsB);
  } else
    for(fcs_int k=-1; k<num_nodes+3; k++)
      table[k+1] = 0.0;
}

static void init_far_interpolation_table_potential_1dp(
    fcs_int num_nodes, fcs_int reg_far,
    fcs_int k, fcs_float periodic_box_length, 
    fcs_float alpha, fcs_float h,
    fcs_float epsB, fcs_int p, fcs_float c,
    fcs_float *table 
    )
{

  fcs_float lknorm = fcs_fabs(k / periodic_box_length);
  fcs_float param[3] = {alpha, lknorm, h}; 
  fcs_float scale2 = 1.0/periodic_box_length;
  fcs_float xroot = FCS_PI * lknorm / alpha;
  fcs_float x = xroot * xroot;

  if(k==0){
    if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
      /* use symmetric Taylor2p, not smooth at r=0.5*h in 3d */
      for(fcs_int l=-1; l<num_nodes+3; l++)
        table[l+1] = -ifcs_p2nfft_reg_far_rad_sym_no_singularity(
            ifcs_p2nfft_ewald_1dp_keq0, param, h*(0.5 - epsB + epsB * (fcs_float) l / num_nodes), p, epsB);
    } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
      /* use normal basis polynomials and set constant continuation value 'c' explicitly */
      /* Since we only evaluate regkernel in h*[0.5-epsB, 0.5] we can savely set epsI = epsB */
      for(fcs_int l=-1; l<num_nodes+3; l++)
        table[l+1] = -ifcs_p2nfft_reg_far_rad_ec_no_singularity(
            ifcs_p2nfft_ewald_1dp_keq0, param, h*(0.5 - epsB + epsB * (fcs_float) l / num_nodes), p, epsB, c);
    } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
      /* use integrated basis polynomials, sets continuation value 'c' implicitily */
      for(fcs_int l=-1; l<num_nodes+3; l++)
        table[l+1] = -ifcs_p2nfft_reg_far_rad_ic_no_singularity(
            ifcs_p2nfft_ewald_1dp_keq0, param, h*(0.5 - epsB + epsB * (fcs_float) l / num_nodes), p, epsB);
    } else
      for(fcs_int l=-1; l<num_nodes+3; l++) table[l+1] = 0.0;
  } else {
    if(2.0 * scale2 * ifcs_p2nfft_ewald_1dp_kneq0(h*(0.5+epsB),0,param) < 1e-16){
      /* use periodization (analytical Fourier coefficients, Poisson summation) instead of regularization,
       * For simplicity, we compute 2d-FFTs on zeros and overwrite the values afterward. */
      for(fcs_int l=-1; l<num_nodes+3; l++) table[l+1] = 0.0;
    } else if( x > 34 ){
      for(fcs_int l=-1; l<num_nodes+3; l++) table[l+1] = 0.0;
    } else if( (FCS_PI * lknorm * h*(0.5-epsB) > 19) ){
      for(fcs_int l=-1; l<num_nodes+3; l++) table[l+1] = 0.0;
    } else if(ifcs_p2nfft_ewald_1dp_kneq0(h*(0.5-epsB), 0, param) < 1e-16){
      for(fcs_int l=-1; l<num_nodes+3; l++) table[l+1] = 0.0;
    } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
      /* use symmetric Taylor2p, not smooth at r=0.5*h in 3d */
      for(fcs_int l=-1; l<num_nodes+3; l++)
        table[l+1] = 2.0 * ifcs_p2nfft_reg_far_rad_sym_no_singularity(
            ifcs_p2nfft_ewald_1dp_kneq0, param, h*(0.5 - epsB + epsB * (fcs_float) l / num_nodes), p, epsB);
    } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
      /* use normal basis polynomials and set constant continuation value 'c' explicitly */
      /* Since we only evaluate regkernel in h*[0.5-epsB, 0.5] we can savely set epsI = epsB */
      for(fcs_int l=-1; l<num_nodes+3; l++){
        table[l+1] = 2.0 * ifcs_p2nfft_reg_far_rad_ec_no_singularity(
            ifcs_p2nfft_ewald_1dp_kneq0, param, h*(0.5 - epsB + epsB * (fcs_float) l / num_nodes), p, epsB, c);

      }
    } else if(reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
      /* use integrated basis polynomials, sets continuation value 'c' implicitily */
      for(fcs_int l=-1; l<num_nodes+3; l++){
        table[l+1] = 2.0 * ifcs_p2nfft_reg_far_rad_ic_no_singularity(
            ifcs_p2nfft_ewald_1dp_kneq0, param, h*(0.5 - epsB + epsB * (fcs_float) l / num_nodes), p, epsB);
    }
    } else
      for(fcs_int l=-1; l<num_nodes+3; l++) table[l+1] = 0.0;
  }
}

static void init_near_interpolation_table_force_3dp(
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
          + 2.0*alpha*FCS_P2NFFT_1_SQRTPI * fcs_exp(- alpha*alpha * r*r)
        ) / r;
  }
}

static void init_near_interpolation_table_force_0dp(
    fcs_int num_nodes,
    fcs_float r_cut, fcs_float epsI, fcs_int p,
    const fcs_float *taylor2p_derive_coeff,
    fcs_int N_cos, fcs_float *cos_coeff, fcs_float *sin_coeff,
    fcs_float *table
    )
{
  if(cos_coeff != NULL){
    fcs_float scale = r_cut/epsI;
    for(fcs_int k=0; k<num_nodes+3; k++)
      table[k] = evaluate_sin_polynomial_1d(
          epsI * (fcs_float) k / num_nodes, N_cos, sin_coeff) / (scale*scale);
  } else if(taylor2p_derive_coeff != NULL) {
    for(fcs_int k=0; k<num_nodes+3; k++)
      table[k] = ifcs_p2nfft_nearfield_correction_taylor2p_derive(
          (fcs_float) k / num_nodes, p, taylor2p_derive_coeff) / (r_cut * r_cut);
  }
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

static fcs_int max_i(fcs_int a, fcs_int b)
{
  return a >= b ? a : b;
}

/* Gives the maximum absolute value of the derivative D[Erf[alpha*x]/x, {x,order}] */
static fcs_float get_derivative_bound_erf(
    fcs_int order, fcs_float alpha
    )
{
  fcs_float val=0.0;

  /* We computed the maximum absolute value numerically via Mathematica. */
  switch(order){
    case 0: val = 1.13; break; /* FindMaxValue[Evaluate[D[Erf[x]/x, {x, 0}]], {x, 1}] */
    case 1: val = 0.43; break; /* FindMaxValue[Evaluate[-D[Erf[x]/x, {x, 1}]], {x, 1}] */
    case 2: val = 0.76; break; /* FindMaxValue[Evaluate[-D[Erf[x]/x, {x, 2}]], {x, 1}] */
    case 3: val = 1.06; break; /* FindMaxValue[Evaluate[D[Erf[x]/x, {x, 3}]], {x, 1}] */
    case 4: val = 2.69; break; /* FindMaxValue[Evaluate[D[Erf[x]/x, {x, 4}]], {x, 1}] */
    case 5: val = 6.01; break; /* FindMaxValue[Evaluate[-D[Erf[x]/x, {x, 5}]], {x, 1}] */
  }

  return alpha * val;
}

static fcs_int calc_interpolation_num_nodes(
    fcs_int interpolation_order, fcs_float eps
    )
{
  switch(interpolation_order){
    case 0: return (fcs_int) 100000; /* TODO: Compute correct number of nodes by error estimate */
    case 1: return (fcs_int) fcs_ceil(1.7/fcs_pow(eps,1.0/2.0));
    case 2: return (fcs_int) fcs_ceil(2.2/fcs_pow(eps,1.0/3.0));
    case 3: return max_i(10, (fcs_int) fcs_ceil(1.4/fcs_pow(eps,1.0/4.0)));
    default: return 0; /* no interpolation */
  }
}


/* Here, we calculate the minimum number of nodes that are necessary to guarantee 
 * the relative error 'eps' during the interpolation on the interval [0,r],
 * see 'Methods of Shape-Preserving Spline-Approximation' by B.I.Kvasov, 2000. */
static fcs_int calc_interpolation_num_nodes_erf(
    fcs_int interpolation_order, fcs_float eps, fcs_float alpha, fcs_float r, unsigned *err
    )
{
  fcs_float N, c, M, M_pot, M_force;
  *err=0;

  /* define constants from Taylor expansion */
  switch(interpolation_order){
    case 0: c = 1.0; break;
    case 1: c = 1.0/8.0; break;
    case 2: c = fcs_sqrt(3)/9.0; break;
    case 3: c = 3.0/128.0; break;
    default: return 0; /* no interpolation */
  }
   
  /* Compute the max. value of the regularization derivative one order higher
   * than the interpolation order. This gives the rest term of the taylor expansion. */ 
  M_pot   = get_derivative_bound_erf(interpolation_order+1, alpha);
  M_force = get_derivative_bound_erf(interpolation_order+2, alpha);

  /* We use the same number of interpolation nodes for potentials and forces.
   * Be sure, that accuracy is fulfilled for both. */
  M = (M_force > M_pot) ? M_force : M_pot;

  N = r * fcs_pow(c*M/fcs_fabs(eps-FCS_P2NFFT_EPS) , 1.0 / (1.0 + interpolation_order) ); 

  /* At least use 16 interpolation points. */
  if(N<16) N = 16.0;

  /* Set maximum number of nodes to avoid overflows during conversion to int. */
  if(N>1e7){
    *err=1;
    N = 1e7;
  }

  return (fcs_int) fcs_ceil(N);
}

static void init_pnfft(
    FCS_PNFFT(plan) *ths, int dim, const ptrdiff_t *N, const ptrdiff_t *n,
    const fcs_float *x_max, int m,
    unsigned pnfft_flags, fcs_int pnfft_intpol_order, fcs_int pnfft_window,
    unsigned pfft_flags, fcs_int pfft_patience,
    MPI_Comm cart_comm_pnfft
    )
{
#if FCS_ENABLE_INFO
  int myrank;
  MPI_Comm_rank(cart_comm_pnfft, &myrank);
#endif
  
  switch(pnfft_window){
    case 0: pnfft_flags |= PNFFT_WINDOW_GAUSSIAN; break;
    case 1: pnfft_flags |= PNFFT_WINDOW_BSPLINE; break;
    case 2: pnfft_flags |= PNFFT_WINDOW_SINC_POWER; break;
    case 3: pnfft_flags |= PNFFT_WINDOW_KAISER_BESSEL; break;
    case 4: pnfft_flags |= PNFFT_WINDOW_BESSEL_I0; break;
    case 5: pnfft_flags |= PNFFT_WINDOW_GAUSSIAN_T; break;
  }

  switch(pnfft_intpol_order){
    case 0: pnfft_flags |= PNFFT_PRE_CONST_PSI; break;
    case 1: pnfft_flags |= PNFFT_PRE_LIN_PSI; break;
    case 2: pnfft_flags |= PNFFT_PRE_QUAD_PSI; break;
    case 3: pnfft_flags |= PNFFT_PRE_CUB_PSI; break;
  }

  switch(pfft_patience){
    case 0: pfft_flags |= PFFT_ESTIMATE; break;
    case 1: pfft_flags |= PFFT_MEASURE; break;
    case 2: pfft_flags |= PFFT_PATIENT; break;
    case 3: pfft_flags |= PFFT_EXHAUSTIVE; break;
    default: pfft_flags |= PFFT_MEASURE; break;
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
    if((pnfft_flags & PNFFT_WINDOW_GAUSSIAN) && (pnfft_flags & PNFFT_USE_FK_GAUSSIAN_T))
      printf("P2NFFT_INFO: Window function: truncated gaussian (-c pnfft_window_name,gaussian_t)\n");
    else if(pnfft_flags & PNFFT_WINDOW_GAUSSIAN)
      printf("P2NFFT_INFO: Window function: gaussian (-c pnfft_window_name,gaussian)\n");
    else if(pnfft_flags & PNFFT_WINDOW_BSPLINE)
      printf("P2NFFT_INFO: Window function: bspline (-c pnfft_window_name,bspline)\n");
    else if(pnfft_flags & PNFFT_WINDOW_SINC_POWER)
      printf("P2NFFT_INFO: Window function: sinc power (-c pnfft_window_name,sinc)\n");
    else if(pnfft_flags & PNFFT_WINDOW_BESSEL_I0)
      printf("P2NFFT_INFO: Window function: bessel_i0 (-c pnfft_window_name,bessel_i0)\n");
    else
      printf("P2NFFT_INFO: Window function: kaiser-bessel (-c pnfft_window_name,kaiser)\n");

    if(pnfft_flags & PNFFT_FFT_IN_PLACE)
      printf("P2NFFT_INFO: inplace FFT: yes (-c pnfft_fft_in_place,1)\n");
    else
      printf("P2NFFT_INFO: inplace FFT: no (-c pnfft_fft_in_place,0)\n");
  
    if(pnfft_intpol_order < 0)
      printf("P2NFFT_INFO: Interpolation order: direct evaluation of window function (-c pnfft_intpol_order,-1)\n");
    else
      printf("P2NFFT_INFO: Interpolation order: %" FCS_LMOD_INT "d (-c pnfft_intpol_order,%" FCS_LMOD_INT "d)\n", pnfft_intpol_order, pnfft_intpol_order);
  }
#endif
#if FCS_ENABLE_DEBUG
  if(!myrank){
    printf("PNFFT_INIT: N = [%td, %td, %td]\n", N[0], N[1], N[2]);
    printf("PNFFT_INIT: n = [%td, %td, %td]\n", n[0], n[1], n[2]);
    printf("PNFFT_INIT: m = %d\n", m);
    printf("PNFFT_INIT: pfft_flags = %u\n", pfft_flags);
    printf("PNFFT_INIT: pnfft_flags = %u\n", pnfft_flags);
    printf("PNFFT_INIT: x_max = [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", x_max[0], x_max[1], x_max[2]);
  }
#endif

  /* Start timing of PNFFT tuning */
  FCS_P2NFFT_INIT_TIMING(cart_comm_pnfft);
  FCS_P2NFFT_START_TIMING(cart_comm_pnfft);
  
  /* finalize, if memory has been already allocated */
  if(*ths != NULL)
    FCS_PNFFT(finalize)(*ths, PNFFT_FREE_F_HAT);

  /* call PNFFT planer */
  *ths = FCS_PNFFT(init_guru)(dim, N, n, x_max, m, pnfft_flags, pfft_flags, cart_comm_pnfft);
  
  /* Finish timing of PNFFT tuning */
  FCS_P2NFFT_FINISH_TIMING(cart_comm_pnfft, "PNFFT tuning");
}
  

/* scale epsI and epsB according to box_size == 1 */
static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_0dp(
    const ptrdiff_t *N, fcs_float r_cut, fcs_float epsI, fcs_float epsB,
    fcs_int p, fcs_float c, fcs_float *box_scales,
    fcs_int reg_near, fcs_int reg_far,
    const fcs_float *taylor2p_coeff, fcs_int N_cg_cos, const fcs_float *cg_cos_coeff,
    fcs_int interpolation_order, fcs_int near_interpolation_num_nodes, fcs_int far_interpolation_num_nodes,
    const fcs_float *near_interpolation_table_potential, const fcs_float *far_interpolation_table_potential,
    MPI_Comm comm_cart, unsigned box_is_cubic
    )
{
  ptrdiff_t howmany = 1, alloc_local, m;
  ptrdiff_t local_Ni[3], local_Ni_start[3], local_No[3], local_No_start[3];
  fcs_float x0, x1, x2, x2norm, xsnorm, scale = 1.0;
  FCS_PFFT(plan) pfft;
  fcs_pnfft_complex *regkern_hat;

  alloc_local = FCS_PFFT(local_size_many_dft)(3, N, N, N, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, PFFT_TRANSPOSED_OUT| PFFT_SHIFTED_IN| PFFT_SHIFTED_OUT,
      local_Ni, local_Ni_start, local_No, local_No_start);

  regkern_hat = FCS_PFFT(alloc_complex)(alloc_local);
  
  pfft = FCS_PFFT(plan_many_dft)(3, N, N, N, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, regkern_hat, regkern_hat, comm_cart,
      PFFT_FORWARD, PFFT_TRANSPOSED_OUT| PFFT_SHIFTED_IN| PFFT_SHIFTED_OUT| PFFT_ESTIMATE);

  for(int t=0; t<3; t++)
    scale *= 1.0 / N[t];

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG || FCS_P2NFFT_DEBUG_REGKERN
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  C csum;
  C csum_global;
  csum = 0.0;
#endif
 
  /* shift FFT output via twiddle factors on the input */ 
  m=0;
  for(ptrdiff_t k0 = local_Ni_start[0]; k0 < local_Ni_start[0] + local_Ni[0]; k0++){
    x0 = (fcs_float) k0 / N[0];
    for(ptrdiff_t k1 = local_Ni_start[1]; k1 < local_Ni_start[1] + local_Ni[1]; k1++){
      x1 = (fcs_float) k1 / N[1];
      for(ptrdiff_t k2 = local_Ni_start[2]; k2 < local_Ni_start[2] + local_Ni[2]; k2++, m++){
        x2 = (fcs_float) k2 / N[2];
        xsnorm = fcs_sqrt(x0*x0+x1*x1+x2*x2);

        /* constant continuation outside the ball with radius 0.5 */
        x2norm = fcs_sqrt(x0*x0*box_scales[0]*box_scales[0]
            + x1*x1*box_scales[1]*box_scales[1]
            + x2*x2*box_scales[2]*box_scales[2]);

        /* constant continuation for radii > 0.5 */
        if (xsnorm > 0.5) xsnorm = 0.5;

        /* calculate near and farfield regularization via interpolation */
        if(x2norm < r_cut){
          if(near_interpolation_num_nodes > 0)
            regkern_hat[m] = ifcs_p2nfft_interpolation_near(
                x2norm, 1.0/r_cut, interpolation_order, near_interpolation_num_nodes, near_interpolation_table_potential);
          else if (reg_near == FCS_P2NFFT_REG_NEAR_CG)
            regkern_hat[m] = epsI/r_cut * evaluate_cos_polynomial_1d(x2norm * epsI/r_cut, N_cg_cos, cg_cos_coeff);
          else if (reg_near == FCS_P2NFFT_REG_NEAR_T2P)
            regkern_hat[m] = ifcs_p2nfft_nearfield_correction_taylor2p(x2norm/r_cut, p, taylor2p_coeff) / r_cut;
        } else if(reg_far_is_radial(reg_far)){
          if(xsnorm < 0.5-epsB) {
            regkern_hat[m] = 1.0/x2norm;
          } else {
            if(!box_is_cubic) {
              /* Noncubic regularization works with unscaled coordinates. */
              regkern_hat[m] = ifcs_p2nfft_reg_far_rad_ec_noncubic(
                  ifcs_p2nfft_one_over_modulus, NULL, x2norm, xsnorm, p, r_cut, epsB, c);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "ifcs_p2nfft_reg_far_rad_ec_noncubic: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            } else {
              /* Cubic regularization works with coordinates that are scaled into unit cube.
               * Therefore, use xsnorm, scale the continuation value 'c', and rescale after evaluation.
               * Note that box_scales are the same in every direction for cubic boxes. */
              if(far_interpolation_num_nodes){
                regkern_hat[m] = ifcs_p2nfft_interpolation_far(
                    xsnorm - 0.5 + epsB, 1.0/epsB, interpolation_order, far_interpolation_num_nodes, far_interpolation_table_potential) / box_scales[0];
                FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "ifcs_p2nfft_interpolation_far: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
              } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_CG){
                regkern_hat[m] = evaluate_cos_polynomial_1d(xsnorm, N_cg_cos, cg_cos_coeff) / box_scales[0];
                FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "evaluate_cos_polynomial_1d: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
              } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
                regkern_hat[m] = ifcs_p2nfft_reg_far_rad_sym(
                    ifcs_p2nfft_one_over_modulus, NULL, xsnorm, p, epsI, epsB) / box_scales[0];
                FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "ifcs_p2nfft_reg_far_rad_sym: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
              } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
                regkern_hat[m] = ifcs_p2nfft_reg_far_rad_ec(
                    ifcs_p2nfft_one_over_modulus, NULL, xsnorm, p, epsI, epsB, c*box_scales[0]) / box_scales[0];
                FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "ifcs_p2nfft_reg_far_rad_ec: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
              } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
                regkern_hat[m] = ifcs_p2nfft_reg_far_rad_ic(
                    ifcs_p2nfft_one_over_modulus, NULL, xsnorm, p, epsI, epsB) / box_scales[0];
                FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "ifcs_p2nfft_reg_far_rad_ic: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
              }
            }
          }
        } else {
          fcs_float x[3] = {x0,x1,x2}, h[3], xi[3], xo[3];
          for(int t=0; t<3; t++){
            x[t] *= box_scales[t];
            h[t] = box_scales[t]; //* (0.5-epsB); /* TODO use d->box_l ? */
          }
          if(reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_SYM){
            for(int t=0; t<3; t++){
              xi[t] = (0.5-epsB)*h[t];
              xo[t] = (0.5+epsB)*h[t];
            }
            regkern_hat[m] = ifcs_p2nfft_interpolate_cuboid_symmetric(p, xi, xo, x);
            FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "ifcs_p2nfft_reg_far_rect_sym: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
          } else if(reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_EC){
            for(int t=0; t<3; t++){
              xi[t] = (0.5-epsB)*h[t];
              xo[t] = 0.5*h[t];
            }
            regkern_hat[m] = ifcs_p2nfft_interpolate_cuboid_explicit_continuation(c, p, xi, xo, x);
            FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "ifcs_p2nfft_reg_far_rect_ec: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
          } else if(reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_IC){
            for(int t=0; t<3; t++){
              xi[t] = (0.5-epsB)*h[t];
              xo[t] = 0.5*h[t];
            }
            regkern_hat[m] = ifcs_p2nfft_interpolate_cuboid_implicit_continuation(p, xi, xo, x);
            FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "ifcs_p2nfft_reg_far_rect_ic: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
          }
        }
        
        regkern_hat[m] *= scale;

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  csum += fabs(creal(regkern_hat[m])) + fabs(cimag(regkern_hat[m]));
#endif
      
      }
    }
  }

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
    MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0) fprintf(stderr, "sum of regkernel: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif
    
  FCS_PFFT(execute)(pfft);

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  /* take care of transposed order N1 x N2 x N0 */
  /* shift FFT input via twiddle factors on the output */
  csum = 0.0;
  m=0;
  for(ptrdiff_t k1 = 0; k1 < local_No[1]; k1++){
    for(ptrdiff_t k2 = 0; k2 < local_No[2]; k2++){
      for(ptrdiff_t k0 = 0; k0 < local_No[0]; k0++, m++){
        csum += fabs(creal(regkern_hat[m])) + fabs(cimag(regkern_hat[m]));
      }
    }
  }

  MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myrank == 0) fprintf(stderr, "sum of regkernel_hat: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif
  
  FCS_PFFT(destroy_plan)(pfft);

  return regkern_hat;
}

/* scale epsI and epsB according to box_size == 1 */
static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_0dp_ewald(
    const ptrdiff_t *N, fcs_float epsB,
//    const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c,
//    const fcs_float *box_inv,
    const fcs_float *box_scales, fcs_float alpha,
    fcs_int p, fcs_float c, fcs_int reg_far,
    fcs_int interpolation_order, fcs_int far_interpolation_num_nodes, const fcs_float *far_interpolation_table_potential,
    MPI_Comm comm_cart, unsigned box_is_cubic
    )
{
  ptrdiff_t howmany = 1, alloc_local, m;
  ptrdiff_t local_Ni[3], local_Ni_start[3], local_No[3], local_No_start[3];
  fcs_float x0, x1, x2, x2norm, xsnorm, scale = 1.0;
  FCS_PFFT(plan) pfft;
  fcs_pnfft_complex *regkern_hat;

  alloc_local = FCS_PFFT(local_size_many_dft)(3, N, N, N, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, PFFT_TRANSPOSED_OUT| PFFT_SHIFTED_IN| PFFT_SHIFTED_OUT,
      local_Ni, local_Ni_start, local_No, local_No_start);

  regkern_hat = FCS_PFFT(alloc_complex)(alloc_local);
  
  pfft = FCS_PFFT(plan_many_dft)(3, N, N, N, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, regkern_hat, regkern_hat, comm_cart,
      PFFT_FORWARD, PFFT_TRANSPOSED_OUT| PFFT_SHIFTED_IN| PFFT_SHIFTED_OUT| PFFT_ESTIMATE);

  //TODO make sure if that makes sense wrt (non-rectangular boxes?)
  for(int t=0; t<3; t++)
    scale /= N[t];

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  C csum;
  C csum_global;
  csum = 0.0;
#endif
 
//  fprintf(stderr, "box_scales: %f %f %f\n", box_scales[0], box_scales[1], box_scales[2]);
  /* shift FFT output via twiddle factors on the input */ 
  m=0;
  for(ptrdiff_t k0 = local_Ni_start[0]; k0 < local_Ni_start[0] + local_Ni[0]; k0++){
    x0 = (fcs_float) k0 / N[0];
    for(ptrdiff_t k1 = local_Ni_start[1]; k1 < local_Ni_start[1] + local_Ni[1]; k1++){
      x1 = (fcs_float) k1 / N[1];
      for(ptrdiff_t k2 = local_Ni_start[2]; k2 < local_Ni_start[2] + local_Ni[2]; k2++, m++){
        x2 = (fcs_float) k2 / N[2];
        xsnorm = fcs_sqrt(x0*x0+x1*x1+x2*x2);

        /* constant continuation outside the ball with radius 0.5 */
        x2norm = fcs_sqrt(x0*x0*box_scales[0]*box_scales[0]
            + x1*x1*box_scales[1]*box_scales[1]
            + x2*x2*box_scales[2]*box_scales[2]);

        //fprintf(stderr, "x: %f %f %f  x2norm: %f  xsnorm: %f \n", x0, x1, x2, x2norm, xsnorm);
        /* constant continuation for radii > 0.5 */
        //TODO different lengths?
//        if (x2norm > (0.5-epsB) * box_scales[0])
//          x2norm = (0.5-epsB) * box_scales[0];
//        fprintf(stderr, "x: %f %f %f  x2norm: %f  xsnorm: %f \n", x0, x1, x2, x2norm, xsnorm);


//        /* calculate farfield regularization via interpolation */
//        if(reg_far_is_radial(reg_far)){
//          if(xsnorm < 0.5-epsB) {
//            regkern_hat[m] = 1.0/x2norm;
//          } else {
//            if(!box_is_cubic) {
//              /* Noncubic regularization works with unscaled coordinates. */
//              regkern_hat[m] = ifcs_p2nfft_reg_far_rad_ec_noncubic(
//                  ifcs_p2nfft_one_over_modulus, NULL, x2norm, xsnorm, p, r_cut, epsB, c);
//            } else {
//              /* Cubic regularization works with coordinates that are scaled into unit cube.
//               * Therefore, use xsnorm, scale the continuation value 'c', and rescale after evaluation.
//               * Note that box_scales are the same in every direction for cubic boxes. */
//              if(far_interpolation_num_nodes){
//                regkern_hat[m] = ifcs_p2nfft_interpolation_far(
//                    xsnorm - 0.5 + epsB, 1.0/epsB, interpolation_order, far_interpolation_num_nodes, far_interpolation_table_potential) / box_scales[0];
//              } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_CG){
//                regkern_hat[m] = evaluate_cos_polynomial_1d(xsnorm, N_cg_cos, cg_cos_coeff) / box_scales[0];
//              } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
//                regkern_hat[m] = ifcs_p2nfft_reg_far_rad_sym(
//                    ifcs_p2nfft_one_over_modulus, NULL, xsnorm, p, epsI, epsB) / box_scales[0];
//              } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
//                regkern_hat[m] = ifcs_p2nfft_reg_far_rad_ec(
//                    ifcs_p2nfft_one_over_modulus, NULL, xsnorm, p, epsI, epsB, c*box_scales[0]) / box_scales[0];
//              } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
//                regkern_hat[m] = ifcs_p2nfft_reg_far_rad_ic(
//                    ifcs_p2nfft_one_over_modulus, NULL, xsnorm, p, epsI, epsB) / box_scales[0];
//              }
//            }
//          }
//        } else {
//          fcs_float x[3] = {x0,x1,x2}, h[3];
//          for(int t=0; t<3; t++){
//            x[t] *= box_scales[t];
//            h[t] = box_scales[t]; //* (0.5-epsB); /* TODO use d->box_l ? */
//          }
//          if(reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_SYM)
//            regkern_hat[m] = ifcs_p2nfft_reg_far_rect_sym(x, h, p, epsB);
//          else if(reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_EC)
//            regkern_hat[m] = ifcs_p2nfft_reg_far_rect_ec(x, h, p, epsB, c);
//          else if(reg_far == FCS_P2NFFT_REG_FAR_REC_T2P_IC)
//            regkern_hat[m] = ifcs_p2nfft_reg_far_rect_ic(x, h, p, epsB);
//        }

          fcs_float param[3];
          param[0] = alpha;
          param[1] = x2norm; //TODO revisit this. Probably wrong?
          param[2] = box_scales[0]; //TODO assumes cubic box

          fcs_float h = box_scales[0];

          if( (far_interpolation_num_nodes > 0) && (h*(0.5-epsB) < x2norm) ){
            /* interpolation requires constant continuation */
            fcs_float x2 = (x2norm > h*0.5) ? h*0.5 : x2norm;
            fcs_float xs = (x2norm > h*0.5) ? 0.5   : xsnorm;
            param[1] = x2;

            regkern_hat[m] = ifcs_p2nfft_interpolation_far(
                xs - 0.5 + epsB, 1.0/epsB, interpolation_order, far_interpolation_num_nodes, far_interpolation_table_potential);
          } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
            regkern_hat[m] = ifcs_p2nfft_reg_far_rad_ic_no_singularity(ifcs_p2nfft_erfx_over_x, param,
                x2norm, p, epsB);
          } else {
            fprintf(stderr, "WRONG REGULARIZATION\n");
	        }

//        regkern_hat[m] = ifcs_p2nfft_reg_far_ec_noncubic(
//             ifcs_p2nfft_one_over_modulus, x2norm, xsnorm, p, NULL, r_cut, epsB, c);
//         regkern_hat[m] = ifcs_p2nfft_reg_far_ec(
//             ifcs_p2nfft_one_over_modulus, xsnorm, p, NULL,  epsI,  epsB, c*box_scales[0]) / box_scales[0];
        
        regkern_hat[m] *= scale;

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  csum += fabs(creal(regkern_hat[m])) + fabs(cimag(regkern_hat[m]));
#endif
      
      }
    }
  }

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
    MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0) fprintf(stderr, "sum of regkernel: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif
    
  FCS_PFFT(execute)(pfft);

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
    csum = 0.0;
#endif

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  /* take care of transposed order N1 x N2 x N0 */
  /* shift FFT input via twiddle factors on the output */
  m=0;
  for(ptrdiff_t k1 = 0; k1 < local_No[1]; k1++){
    for(ptrdiff_t k2 = 0; k2 < local_No[2]; k2++){
      for(ptrdiff_t k0 = 0; k0 < local_No[0]; k0++, m++){
        csum += fabs(creal(regkern_hat[m])) + fabs(cimag(regkern_hat[m]));
      }
    }
  }
#endif

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myrank == 0) fprintf(stderr, "sum of regkernel_hat: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif
  
  FCS_PFFT(destroy_plan)(pfft);

  return regkern_hat;
}

static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_3dp(
    const ptrdiff_t *local_N, const ptrdiff_t *local_N_start,
    const fcs_float *ibox, fcs_float alpha, fcs_float kc
    )
{
  ptrdiff_t m, k[3], alloc_local = 1;
  fcs_float ivol = fcs_fabs(det_3x3(ibox+0, ibox+3, ibox+6));
  fcs_pnfft_complex *regkern_hat;
 
  for(int t=0; t<3; t++)
    alloc_local *= local_N[t];
  regkern_hat = FCS_PFFT(alloc_complex)(alloc_local);
  
  m = 0;
  for(k[1]=local_N_start[1]; k[1] < local_N_start[1] + local_N[1]; k[1]++){
    for(k[2]=local_N_start[2]; k[2] < local_N_start[2] + local_N[2]; k[2]++){
      for(k[0]=local_N_start[0]; k[0] < local_N_start[0] + local_N[0]; k[0]++, m++){

        fcs_float kf0 = (fcs_float) k[0];
        fcs_float kf1 = (fcs_float) k[1];
        fcs_float kf2 = (fcs_float) k[2];
        fcs_float ksqnorm = kf0*kf0 + kf1*kf1 + kf2*kf2;

        if ((k[0] == 0) && (k[1] == 0) && (k[2] == 0))
          regkern_hat[m] = 0;
        else if (kc > 0.0 && ksqnorm > kc*kc)
          regkern_hat[m] = 0; /* Apply spherical cutoff for kc > 0 */
        else {
          fcs_float lk0 = At_TIMES_VEC(ibox, k, 0);
          fcs_float lk1 = At_TIMES_VEC(ibox, k, 1);
          fcs_float lk2 = At_TIMES_VEC(ibox, k, 2);
          fcs_float lksqnorm = lk0*lk0 + lk1*lk1 + lk2*lk2;

          regkern_hat[m] = ivol * fcs_exp(-FCS_P2NFFT_PISQR * lksqnorm/(alpha*alpha))/lksqnorm/FCS_P2NFFT_PI;
        }
      }
    }
  }
  
  return regkern_hat;
}

/* scale epsI and epsB according to box_size == 1 */
static fcs_pnfft_complex* malloc_and_precompute_regkern_hat_2dp_and_1dp(
    const ptrdiff_t *N, fcs_float epsB,
    const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c,
    const fcs_float *box_inv, const fcs_float *box_scales, fcs_float alpha, fcs_float kc,
    const fcs_int *periodicity, fcs_int p, fcs_float c, fcs_int reg_far,
    fcs_int interpolation_order, fcs_int far_interpolation_num_nodes, 
    const fcs_float *far_interpolation_table_potential,
    MPI_Comm comm_cart
    )
{
  ptrdiff_t howmany = 1, alloc_local, m;
  ptrdiff_t local_Ni[3], local_Ni_start[3], local_No[3], local_No_start[3], k[3];
  fcs_float lk[3];
  fcs_int num_periodic_dims = (periodicity[0]!=0) + (periodicity[1]!=0) + (periodicity[2]!=0);
  fcs_float scale = 1.0;
  fcs_float scale2 = 1.0;
  fcs_float xs[3];
  FCS_PFFT(plan) pfft;
  fcs_pnfft_complex *regkern_hat;

  alloc_local = FCS_PFFT(local_size_many_dft)(3, N, N, N, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, PFFT_TRANSPOSED_OUT| PFFT_SHIFTED_IN| PFFT_SHIFTED_OUT,
      local_Ni, local_Ni_start, local_No, local_No_start);

  regkern_hat = FCS_PFFT(alloc_complex)(alloc_local);

  int skipped_dims[3];
  for(int t=0; t<3; t++)
    skipped_dims[t] = periodicity[t];

  pfft = FCS_PFFT(plan_many_dft_skipped)(3, N, N, N, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, skipped_dims, regkern_hat, regkern_hat, comm_cart,
      PFFT_FORWARD, PFFT_TRANSPOSED_OUT| PFFT_SHIFTED_IN| PFFT_SHIFTED_OUT| PFFT_ESTIMATE);

  /* compute the volume element corresponding to periodic dims */
  scale /= volume_of_periodic_dims(box_a, box_b, box_c, periodicity);
  scale2 = scale;

  /* this works only for equal nonperiodic box lengths */
  fcs_float h=0.0;
  for(fcs_int t=0; t<3; t++)
    if(!periodicity[t])
      h = box_scales[t];

  /* divide by FFT normalization factor for every nonperiodic dim */
  for(int t=0; t<3; t++) if(!periodicity[t]) scale /= N[t];

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG || FCS_P2NFFT_DEBUG_REGKERN
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  C csum, csum_global;
  csum = 0.0;
  double timer_keq0=0, timer_kne0=0, timer_analytic=0;
  int count_Req0=0, count_Rne0=0, count_all=0, count_analytic=0, count_interpolate=0, count_set_to_zero=0, count_far_interpolate=0;
#endif

// double _time_=-MPI_Wtime();
// fcs_float _alpha_=1.0, _k_ = fcs_sqrt(1.0)/FCS_PI, _r_ = fcs_sqrt(1.1);
// fcs_int _der_=0;
// fcs_float _param_[2] = {_alpha_, _k_};
// fcs_float kern = ifcs_p2nfft_ewald_1dp_kneq0(_r_, _der_, _param_);
// _time_+=MPI_Wtime();
// fprintf(stderr, "\n ifcs_p2nfft_ewald_1dp_kneq0 = %e, time for evaluation = %e\n\n", kern, _time_);
// return regkern_hat;

  /* shift FFT output via twiddle factors on the input */
  /* twiddle only the non-periodic dims, since there we need to calculate 1d-FFTs */ 
  m=0;
  for(ptrdiff_t l0 = local_Ni_start[0]; l0 < local_Ni_start[0] + local_Ni[0]; l0++){
    xs[0] = (periodicity[0]) ? 0.0 : (fcs_float) l0 / N[0];
    k[0] = (periodicity[0]) ? l0 : 0;  
    for(ptrdiff_t l1 = local_Ni_start[1]; l1 < local_Ni_start[1] + local_Ni[1]; l1++){
      xs[1] = (periodicity[1]) ? 0.0 : (fcs_float) l1 / N[1];
      k[1] = (periodicity[1]) ? l1 : 0;  
      for(ptrdiff_t l2 = local_Ni_start[2]; l2 < local_Ni_start[2] + local_Ni[2]; l2++, m++){
        xs[2] = (periodicity[2]) ? 0.0 : (fcs_float) l2 / N[2];
        k[2] = (periodicity[2]) ? l2 : 0;  

        /* New regularization for mixed boundary conditions */
        fcs_float lknorm = 0.0, x2norm = 0.0, xsnorm = 0.0;
	fcs_float kxsqnorm = 0.0;

        lk[0] = At_TIMES_VEC(box_inv, k, 0);
        lk[1] = At_TIMES_VEC(box_inv, k, 1);
        lk[2] = At_TIMES_VEC(box_inv, k, 2);

        fcs_float kf0 = (fcs_float) k[0];
        fcs_float kf1 = (fcs_float) k[1];
        fcs_float kf2 = (fcs_float) k[2];
        fcs_float ksqnorm = kf0*kf0 + kf1*kf1 + kf2*kf2;

        for(fcs_int t=0; t<3; t++){
          lknorm += lk[t] * lk[t];
          x2norm += xs[t] * xs[t] * box_scales[t] * box_scales[t];
          xsnorm += xs[t] * xs[t];
	  kxsqnorm += lk[t] * lk[t] + xs[t]*xs[t]/(box_scales[t]*box_scales[t]*N[t]*N[t]);
        }
        lknorm = fcs_sqrt(lknorm);
        x2norm = fcs_sqrt(x2norm);
        xsnorm = fcs_sqrt(xsnorm);

        fcs_float param[3] = {alpha, lknorm, h};

        /* Check if all indices corresponding to periodic dims are 0.
         * Index correspoding to non-periodic dim will be 0 per default. */
        if ((k[0] == 0) && (k[1] == 0) && (k[2] == 0)){
          FCS_P2NFFT_IFDBG(double tmp = -MPI_Wtime());
          if(num_periodic_dims == 2){
            if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
              regkern_hat[m] = -2.0 * FCS_SQRTPI * ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_2dp_keq0, param, x2norm, p, epsB);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k==0, ifcs_p2nfft_reg_far_rad_sym_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
              regkern_hat[m] = -2.0 * FCS_SQRTPI * ifcs_p2nfft_reg_far_rad_ec_no_singularity(ifcs_p2nfft_ewald_2dp_keq0, param, x2norm, p, epsB, c);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k==0, ifcs_p2nfft_reg_far_rad_ec_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
              regkern_hat[m] = -2.0 * FCS_SQRTPI * ifcs_p2nfft_reg_far_rad_ic_no_singularity(ifcs_p2nfft_ewald_2dp_keq0, param, x2norm, p, epsB);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k==0, ifcs_p2nfft_reg_far_rad_ic_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            }
          } else {
            if( far_interpolation_table_potential && (h*(0.5-epsB) < x2norm) ){
              /* interpolation requires constant continuation */
              fcs_float x2 = (x2norm > h*0.5) ? h*0.5 : x2norm;
              fcs_float xs = (x2norm > h*0.5) ? 0.5   : xsnorm;
              param[1] = x2;

              fcs_int pdim = -1;
              for(int t=0; t<3; t++) if(periodicity[t]) pdim=t;

              fcs_int ind = k[pdim] - local_Ni_start[pdim];
              fcs_int offset = far_interpolation_num_nodes + 4;

              regkern_hat[m] = ifcs_p2nfft_interpolation_far(
                  xs - 0.5 + epsB, 1.0/epsB, interpolation_order, far_interpolation_num_nodes, far_interpolation_table_potential + ind * offset);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k==0, ifcs_p2nfft_interpolation_far: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
              regkern_hat[m] = -ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_1dp_keq0, param, x2norm, p, epsB);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k==0, ifcs_p2nfft_reg_far_rad_sym_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
              regkern_hat[m] = -ifcs_p2nfft_reg_far_rad_ec_no_singularity(ifcs_p2nfft_ewald_1dp_keq0, param, x2norm, p, epsB, c);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k==0, ifcs_p2nfft_reg_far_rad_ec_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
              regkern_hat[m] = -ifcs_p2nfft_reg_far_rad_ic_no_singularity(ifcs_p2nfft_ewald_1dp_keq0, param, x2norm, p, epsB);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k==0, ifcs_p2nfft_reg_far_rad_ic_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            }

            if(isnan(creal(regkern_hat[m]))){
              fprintf(stderr, "keq0: k = [%td, %td, %td], x = [%e, %e, %e], p = %" FCS_LMOD_INT "d, epsB = %e, x2norm = %e, lknorm = %e, alpha = %e\n",
                  k[0], k[1], k[2], xs[0], xs[1], xs[2], p, epsB, x2norm, lknorm, alpha);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          }
          FCS_P2NFFT_IFDBG(++count_interpolate);
          FCS_P2NFFT_IFDBG(tmp += MPI_Wtime(); timer_keq0 += tmp);
        }
        else if (kc > 0.0 && ksqnorm > kc*kc){
          regkern_hat[m] = 0; /* Apply spherical cutoff for kc > 0 */
          FCS_P2NFFT_IFDBG(++count_set_to_zero);
	}
        else { /* k != 0 */
          FCS_P2NFFT_IFDBG(double tmp = -MPI_Wtime());

          if(num_periodic_dims == 2){
	    if(0.5 * scale2 * ifcs_p2nfft_ewald_2dp_kneq0(h*(0.5+epsB),0,param) < 1e-16){
	      /* use periodization (analytical Fourier coefficients, Poisson summation) instead of regularization,
               * For simpilicity, we compute 1d-FFTs on zeros and overwrite the values afterward. */
	      regkern_hat[m] = 0.0;
              FCS_P2NFFT_IFDBG(++count_analytic);
	      // kern_hat_ana[m] = 0.5*2.0/(FCS_PI*kxsqnorm) * fcs_exp(-FCS_P2NFFT_PISQR * kxsqnorm/(alpha*alpha));
	    } else if( (x2norm > h*(0.5-epsB)) && (FCS_PI * lknorm * h*(0.5-epsB) > 19) ){
              /* avoid regularization of functions that are numerically equal to zero (less than 1e-16) */
              regkern_hat[m] = 0.0;
              FCS_P2NFFT_IFDBG(++count_set_to_zero);
            } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
              regkern_hat[m] = 0.5 * ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_2dp_kneq0, param, x2norm, p, epsB) / lknorm;
	      FCS_P2NFFT_IFDBG(++count_interpolate);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k!=0, ifcs_p2nfft_reg_far_rad_sym_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
              regkern_hat[m] = 0.5 * ifcs_p2nfft_reg_far_rad_ec_no_singularity(ifcs_p2nfft_ewald_2dp_kneq0, param, x2norm, p, epsB, c) / lknorm;
	      FCS_P2NFFT_IFDBG(++count_interpolate);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k!=0, ifcs_p2nfft_reg_far_rad_ec_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
              regkern_hat[m] = 0.5 * ifcs_p2nfft_reg_far_rad_ic_no_singularity(ifcs_p2nfft_ewald_2dp_kneq0, param, x2norm, p, epsB) / lknorm;
	      FCS_P2NFFT_IFDBG(++count_interpolate);
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k!=0, ifcs_p2nfft_reg_far_rad_ic_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            }

//             if((k[0] == 0) && (k[1] == 1) && (k[2] == 0))
//               fprintf(stderr, "regkern = %e, x2norm = %e, h = %f\n", creal(2/box_l[1]*regkern_hat[m]), x2norm, box_scales[0]);

//             fprintf(stderr, "regkern = %e, x2norm = %e\n", creal(regkern_hat[m]), x2norm);
//             fprintf(stderr, "kne0: k = [%td, %td, %td], x = [%e, %e, %e], p = %d, epsB = %e, xsnorm = %e, lknorm = %e, alpha = %e\n",
//                 k[0], k[1], k[2], xs[0], xs[1], xs[2], p, epsB, xsnorm, lknorm, alpha);
          } else { /* num_periodic_dims == 1 */
            fcs_float xroot = FCS_PI * lknorm / alpha;
            fcs_float x = xroot * xroot;

//             if(k[2] >= 10) fprintf(stderr, "k = [%td, %td, %td], lknorm = %e, alpha = %e, Pi*k/(alpha*B) = %e\n", k[0], k[1], k[2], lknorm, alpha, FCS_PI*lknorm/alpha);

	    if(2.0 * scale2 * ifcs_p2nfft_ewald_1dp_kneq0(h*(0.5+epsB),0,param) < 1e-16){
	      /* use periodization (analytical Fourier coefficients, Poisson summation) instead of regularization,
               * For simplicity, we compute 2d-FFTs on zeros and overwrite the values afterward. */
	      regkern_hat[m] = 0.0;
              FCS_P2NFFT_IFDBG(++count_analytic);
	      // kern_hat_ana[m] = 2.0/(2.0*FCS_PI*kxsqnorm) * fcs_exp(-FCS_P2NFFT_PISQR * kxsqnorm/(alpha*alpha));
	    } else if( x > 34 ){
              /* avoid regularization of functions that are numerically equal to zero (less than 1e-16) */
              regkern_hat[m] = 0.0;
	      FCS_P2NFFT_IFDBG(++count_set_to_zero);
            } else if( (x2norm > h*(0.5-epsB)) && ( FCS_PI*lknorm * h*(0.5-epsB) > 19) ){
              /* avoid regularization of functions that are numerically equal to zero (less than 1e-16) */
              regkern_hat[m] = 0.0;
	      FCS_P2NFFT_IFDBG(++count_set_to_zero);
            } else if( (x2norm > h*(0.5-epsB)) && (ifcs_p2nfft_ewald_1dp_kneq0(h*(0.5-epsB), 0, param) < 1e-16) ){
              /* avoid regularization of functions that are numerically equal to zero (less than 1e-16) */
              regkern_hat[m] = 0.0;
	      FCS_P2NFFT_IFDBG(++count_set_to_zero);
            } else if( far_interpolation_table_potential && (h*(0.5-epsB) < x2norm) ){
                  /* interpolation requires constant continuation */
                  fcs_float x2 = (x2norm > h*0.5) ? h*0.5 : x2norm;
                  fcs_float xs = (x2norm > h*0.5) ? 0.5   : xsnorm;
                  param[1] = x2;

                  fcs_int pdim = -1;
                  for(int t=0; t<3; t++) if(periodicity[t]) pdim=t;

                  fcs_int ind = k[pdim] - local_Ni_start[pdim];
                  fcs_int offset = far_interpolation_num_nodes + 4;

                  regkern_hat[m] = ifcs_p2nfft_interpolation_far(
                      xs - 0.5 + epsB, 1.0/epsB, interpolation_order, far_interpolation_num_nodes, far_interpolation_table_potential + ind * offset);
                  FCS_P2NFFT_IFDBG(++count_far_interpolate);
                  FCS_P2NFFT_IFDBG(++count_interpolate);
                  FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k!=0, ifcs_p2nfft_interpolation_far: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            } else {
              /* The function evaluations 'in the middle' of the far field are much more expensive than the rest.
               * Here, we use symmetry to reduce the number of expensive function evaluations. 
               * At least in the serial case every regkern value is computed 8 times on one process. 
               * Let's compute it only once and copy afterward. */
              fcs_int found_symmetric_value = 0;
              
              if( (l0>0 || l1>0 || l2>0) ){
                ptrdiff_t k0 = -abs(l0) - local_Ni_start[0];
                ptrdiff_t k1 = -abs(l1) - local_Ni_start[1];
                ptrdiff_t k2 = -abs(l2) - local_Ni_start[2];

                if(0 <= k0 && k0 < local_Ni[0]){
                  if(0 <= k1 && k1 < local_Ni[1]){
                    if(0 <= k2 && k2 < local_Ni[2]){
                      ptrdiff_t ind = k2 + local_Ni[2]* (k1 + local_Ni[1]*k0);
                      regkern_hat[m] = regkern_hat[ind]/scale;
                      found_symmetric_value = 1;
                    }
                  }
                }
              }

              if(!found_symmetric_value){
                if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM){
                  regkern_hat[m] = 2.0 * ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_1dp_kneq0, param, x2norm, p, epsB);
                  FCS_P2NFFT_IFDBG(++count_interpolate);
                  FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k!=0, ifcs_p2nfft_reg_far_rad_sym_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
                } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC){
                  regkern_hat[m] = 2.0 * ifcs_p2nfft_reg_far_rad_ec_no_singularity(ifcs_p2nfft_ewald_1dp_kneq0, param, x2norm, p, epsB, c);
                  FCS_P2NFFT_IFDBG(++count_interpolate);
                  FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k!=0, ifcs_p2nfft_reg_far_rad_ec_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
                } else if (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC){
                  regkern_hat[m] = 2.0 * ifcs_p2nfft_reg_far_rad_ic_no_singularity(ifcs_p2nfft_ewald_1dp_kneq0, param, x2norm, p, epsB);
                  FCS_P2NFFT_IFDBG(++count_interpolate);
                  FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k!=0, ifcs_p2nfft_reg_far_rad_ic_no_singularity: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
                }
              }
            }

            if(isnan(creal(regkern_hat[m]))){
              fprintf(stderr, "kne0: k = [%td, %td, %td], x = [%e, %e, %e], p = %" FCS_LMOD_INT "d, epsB = %e, xsnorm = %e, lknorm = %e, alpha = %e\n",
                  k[0], k[1], k[2], xs[0], xs[1], xs[2], p, epsB, xsnorm, lknorm, alpha);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          }
          FCS_P2NFFT_IFDBG(tmp += MPI_Wtime(); timer_kne0 += tmp);
        }

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
              if(creal(regkern_hat[m]) < 1e-15) ++count_Req0;
              else                              ++count_Rne0;
#endif

        regkern_hat[m] *= scale;

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
        csum += fabs(creal(regkern_hat[m])) + fabs(cimag(regkern_hat[m]));
        ++count_all;
#endif
      }
    }
  }

//   {
//     FILE *file = fopen("reg_keq0.txt", "w");
//     fcs_float h = box_scales[0];
//     fcs_float k = 2.0/box_l[0];
//     fcs_float param[3];
//     param[0] = alpha;
//     param[1] = k;
//     param[2] = h;
//     for(int t=0; t<20; t++){
//       fcs_float x = (0.5-epsB)*t/20.0*h;
//       fprintf(file, "%e %e\n", x, -ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_1dp_keq0, param, x, p, epsB));
//     }
//     for(int t=0; t<=100; t++){
//       fcs_float x = (0.5-epsB + epsB*t/100.0)*h;
//       fprintf(file, "%e %e\n", x, -ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_1dp_keq0, param, x, p, epsB));
//     }
//     for(int t=1; t<=5; t++){
//       fcs_float x = (0.5 + 0.1*t/5.0)*h;
//       fprintf(file, "%e %e\n", x, -ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_1dp_keq0, param, x, p, epsB));
//     }
//     fclose(file);
//   }
//   {
//     FILE *file = fopen("reg_kneq0.txt", "w");
//     fcs_float h = box_scales[0];
//     fcs_float k = 3.0/box_l[0];
//     fcs_float param[3];
//     param[0] = alpha;
//     param[1] = k;
//     param[2] = h;
//     for(int t=0; t<20; t++){
//       fcs_float x = (0.5-epsB)*t/20.0*h;
//       fprintf(file, "%e %e\n", x, -ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_1dp_kneq0, param, x, p, epsB));
//     }
//     for(int t=0; t<=100; t++){
//       fcs_float x = (0.5-epsB + epsB*t/100.0)*h;
//       fprintf(file, "%e %e\n", x, -ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_1dp_kneq0, param, x, p, epsB));
//     }
//     for(int t=1; t<=5; t++){
//       fcs_float x = (0.5 + 0.1*t/5.0)*h;
//       fprintf(file, "%e %e\n", x, -ifcs_p2nfft_reg_far_rad_sym_no_singularity(ifcs_p2nfft_ewald_1dp_kneq0, param, x, p, epsB));
//     }
//     fclose(file);
//   }

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
    MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0) fprintf(stderr, "sum of regkernel: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif
    
  FCS_PFFT(execute)(pfft);
  
  FCS_PFFT(destroy_plan)(pfft);
  
  /* Poisson summation, take analytical Fourier coefficients */
  m=0;
  for(ptrdiff_t l1 = local_No_start[1]; l1 < local_No_start[1] + local_No[1]; ++l1){
    xs[1] = (periodicity[1]) ? 0.0 : (fcs_float) l1;
    k[1] = (periodicity[1]) ? l1 : 0;
    for(ptrdiff_t l2 = local_No_start[2]; l2 < local_No_start[2] + local_No[2]; ++l2){
      xs[2] = (periodicity[2]) ? 0.0 : (fcs_float) l2;
      k[2] = (periodicity[2]) ? l2 : 0;
      for(ptrdiff_t l0 = local_No_start[0]; l0 < local_No_start[0] + local_No[0]; ++l0, ++m){
	xs[0] = (periodicity[0]) ? 0.0 : (fcs_float) l0;
	k[0] = (periodicity[0]) ? l0 : 0;
	
        /* New regularization for mixed boundary conditions */
        fcs_float lknorm = 0.0, kxsqnorm = 0.0;
	
        lk[0] = At_TIMES_VEC(box_inv, k, 0);
        lk[1] = At_TIMES_VEC(box_inv, k, 1);
        lk[2] = At_TIMES_VEC(box_inv, k, 2);

        fcs_float kf0 = (fcs_float) k[0];
        fcs_float kf1 = (fcs_float) k[1];
        fcs_float kf2 = (fcs_float) k[2];
        fcs_float ksqnorm = kf0*kf0 + kf1*kf1 + kf2*kf2;

        for(fcs_int t=0; t<3; t++){
          lknorm   += lk[t] * lk[t];
          kxsqnorm += lk[t] * lk[t] + xs[t]*xs[t] / (h*h); 
        }
        lknorm = fcs_sqrt(lknorm);

        fcs_float param[3];
        param[0] = alpha;
        param[1] = lknorm;
        param[2] = h;

        if ((k[0] == 0) && (k[1] == 0) && (k[2] == 0)){
          ; /* do not use analytic formula for k==0 */
        }
        else if (kc > 0.0 && ksqnorm > kc*kc){
          ; /* apply spherical cutoff in kspace */
        } else { /* no spherical cutoff, take full grid */
          FCS_P2NFFT_IFDBG(double tmp = -MPI_Wtime());
          if(num_periodic_dims == 2){
            if(0.5 * scale2 * ifcs_p2nfft_ewald_2dp_kneq0(h*(0.5+epsB),0,param) < 1e-16){
              regkern_hat[m] = scale2/(h*FCS_PI*kxsqnorm) * fcs_exp(-FCS_P2NFFT_PISQR * kxsqnorm/(alpha*alpha));
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k!=0, analytic Fourier transform: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            }
          } else { /* num_periodic_dims == 1 */
            if(2.0 * scale2 * ifcs_p2nfft_ewald_1dp_kneq0(h*(0.5+epsB),0,param) < 1e-16){
              regkern_hat[m] = scale2/(h*h*FCS_PI*kxsqnorm) * fcs_exp(-FCS_P2NFFT_PISQR * kxsqnorm/(alpha*alpha));
              FCS_P2NFFT_IFDBG_REGKERN(if(myrank==0) fprintf(stderr, "k!=0, analytic Fourier transform: regkern[%td] = %e + I * %e\n", m, creal(regkern_hat[m]), cimag(regkern_hat[m])));
            }
          }
          FCS_P2NFFT_IFDBG(tmp += MPI_Wtime(); timer_analytic += tmp);
        }

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
        csum += fabs(creal(regkern_hat[m])) + fabs(cimag(regkern_hat[m]));
#endif
      }
    }
  }
  
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  if(myrank==0){
    fprintf(stderr, "\nRegularization for k==0 takes %e\n", timer_keq0);
    fprintf(stderr,   "Regularization for k!=0 takes %e\n", timer_kne0);
    fprintf(stderr,   "Analytic Reg.  for k!=0 takes %e\n\n", timer_analytic);

    fprintf(stderr, "# grid points in the regularization domain = %d\n", count_all);  
    fprintf(stderr, "# regularization values that are set to zero = %d\n", count_set_to_zero); 
    fprintf(stderr, "# regularization values that come from extension with two-point Taylor = %d\n", count_interpolate); 
    fprintf(stderr, "# far field regularizations evaluted with interpolation table = %d\n", count_far_interpolate); 
    fprintf(stderr, "# regularization values that come from analytical Fourier transform = %d\n\n", count_analytic);  

    fprintf(stderr, "Regularization numerically     equal to 0 for %d times\n", count_Req0);
    fprintf(stderr, "Regularization numerically NOT equal to 0 for %d times\n\n", count_Rne0);
  }
#endif

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  /* take care of transposed order N1 x N2 x N0 */
  /* shift FFT input via twiddle factors on the output */
  csum = 0.0;
  m=0;
  for(ptrdiff_t k1 = 0; k1 < local_No[1]; k1++){
    for(ptrdiff_t k2 = 0; k2 < local_No[2]; k2++){
      for(ptrdiff_t k0 = 0; k0 < local_No[0]; k0++, m++){
        csum += fabs(creal(regkern_hat[m])) + fabs(cimag(regkern_hat[m]));
      }
    }
  }

  MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myrank == 0) fprintf(stderr, "sum of regkernel_hat: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif
  
  return regkern_hat;
}

static int reg_far_is_radial(
    fcs_int reg_far
    )
{
  return (reg_far == FCS_P2NFFT_REG_FAR_RAD_CG)
    || (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_SYM)
    || (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_EC)
    || (reg_far == FCS_P2NFFT_REG_FAR_RAD_T2P_IC);
}

static int pnfft_is_up_to_date(
    const FCS_PNFFT(plan) ths, int dim, const ptrdiff_t *N, const ptrdiff_t *n,
    const fcs_float *x_max, int m, unsigned pnfft_flags, unsigned pfft_flags
    )
{
  int plan_d, plan_m;
  ptrdiff_t plan_N[3], plan_n[3];
  fcs_float plan_x_max[3];
  unsigned plan_pnfft_flags, plan_pfft_flags;

#if FCS_P2NFFT_DEBUG_RETUNE
    fprintf(stderr, "P2NFFT_DEBUG: pnfft_is_up_to_date starts with: ths==%p\n", ths);
#endif

  /* plan is uninitialized */
  if( ths == NULL )
    return 0;

  plan_d = FCS_PNFFT(get_d)(ths);
  plan_m = FCS_PNFFT(get_m)(ths);
  FCS_PNFFT(get_N)(ths, plan_N);
  FCS_PNFFT(get_n)(ths, plan_n);
  FCS_PNFFT(get_x_max)(ths, plan_x_max);
  plan_pnfft_flags = FCS_PNFFT(get_pnfft_flags)(ths);
  plan_pfft_flags = FCS_PNFFT(get_pfft_flags)(ths);

#if FCS_P2NFFT_DEBUG_RETUNE
  fprintf(stderr, "P2NFFT_DEBUG: pnfft_is_up_to_date: plan_d = %" FCS_LMOD_INT "d, dim = %" FCS_LMOD_INT "d, plan_m = %" FCS_LMOD_INT "d, m = %" FCS_LMOD_INT "d\n", plan_d, dim, plan_m, m);
#endif

  /* check values of d, m */
  if( plan_d != dim
      || plan_m != m
    )
    return 0;

#if FCS_P2NFFT_DEBUG_RETUNE
  fprintf(stderr, "P2NFFT_DEBUG: pnfft_is_up_to_date: plan_N = [%td, %td, %td], N = [%td, %td, %td]\n", plan_N[0], plan_N[1], plan_N[2], N[0], N[1], N[2]);
  fprintf(stderr, "P2NFFT_DEBUG: pnfft_is_up_to_date: plan_n = [%td, %td, %td], n = [%td, %td, %td]\n", plan_n[0], plan_n[1], plan_n[2], n[0], n[1], n[2]);
  fprintf(stderr, "P2NFFT_DEBUG: pnfft_is_up_to_date: plan_x_max = [%f, %f, %f], x_max = [%f, %f, %f]\n", plan_x_max[0], plan_x_max[1], plan_x_max[2], x_max[0], x_max[1], x_max[2]);
#endif

  /* check values of N, n, x_max */
  for(int t=0; t<dim; t++)
    if( plan_N[t] != N[t]
        || plan_n[t] != n[t]
        || !fcs_float_is_zero(plan_x_max[t] - x_max[t])
      )
      return 0;

#if FCS_P2NFFT_DEBUG_RETUNE
  fprintf(stderr, "P2NFFT_DEBUG: pnfft_is_up_to_date: plan_pnfft_flags = %u, pnfft_flags = %u\n", plan_pnfft_flags, pnfft_flags);
#endif

  /* check PNFFT flags */
  if(plan_pnfft_flags != pnfft_flags)
    return 0;

#if FCS_P2NFFT_DEBUG_RETUNE
  fprintf(stderr, "P2NFFT_DEBUG: pnfft_is_up_to_date: plan_pfft_flags = %u, pfft_flags = %u\n", plan_pfft_flags, pfft_flags);
#endif

  /* check PFFT flags */
  if(plan_pfft_flags != pfft_flags)
    return 0;

  return 1;  
}

static void default_tolerance_type(
    fcs_int *periodicity, fcs_int reg_kernel,
    fcs_int *tolerance_type, fcs_float *tolerance
    )
{
  if(!periodicity[0] && !periodicity[1] && !periodicity[2] && reg_kernel == FCS_P2NFFT_REG_KERNEL_ONE_OVER_ABS_X){
    if(*tolerance < 0.0)
      *tolerance = FCS_P2NFFT_DEFAULT_TOLERANCE;
    if(*tolerance_type == FCS_TOLERANCE_TYPE_UNDEFINED)
      *tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL;
  } else {
    if(*tolerance < 0.0)
      *tolerance = FCS_P2NFFT_DEFAULT_TOLERANCE;
    if(*tolerance_type == FCS_TOLERANCE_TYPE_UNDEFINED)
      *tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
  }
}

static FCSResult check_tolerance(
    fcs_int *periodicity, fcs_int reg_kernel, fcs_int tolerance_type, fcs_float tolerance
    )
{
  char* fnc_name = "check_tolerance";

  if( (tolerance_type == FCS_TOLERANCE_TYPE_ENERGY)
      || (tolerance_type == FCS_TOLERANCE_TYPE_ENERGY_REL)
      || (tolerance_type == FCS_TOLERANCE_TYPE_POTENTIAL_REL)
      || (tolerance_type == FCS_TOLERANCE_TYPE_FIELD_REL) )
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,"P2NFFT does not support this kind of tolerance.");

  if(tolerance_type == FCS_TOLERANCE_TYPE_POTENTIAL)
    if( periodicity[0] || periodicity[1] || periodicity[2] || reg_kernel == FCS_P2NFFT_REG_KERNEL_EWALD)
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,"P2NFFT supports FCS_TOLERANCE_POTENTIAL only for non-periodic boundary conditions and 1/x regularized kernel. Use FCS_TOLERANCE_FIELD instead.");

  if(tolerance_type == FCS_TOLERANCE_TYPE_FIELD)
    if( !periodicity[0] && !periodicity[1] && !periodicity[2] && reg_kernel != FCS_P2NFFT_REG_KERNEL_EWALD )
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,"P2NFFT supports FCS_TOLERANCE_FIELD only for 1d-, 2d- and 3d-periodic boundary conditions and non-periodic boundary conditions with Ewald regularized kernel. Use FCS_TOLERANCE_POTENTIAL instead.");

  if(tolerance < 0.0)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,"Tolerance must be non-negative.");

  return NULL;
}

static fcs_int is_cubic(
    fcs_float *box_l
    )
{
  if( fcs_float_is_equal(box_l[0],box_l[1])
      && fcs_float_is_equal(box_l[0],box_l[2]))
    return 1;

  return 0;
}

/* look for the smallest box length with pbc */
static int get_dim_of_smallest_periodic_box_l(
    fcs_int periodicity[3], fcs_float box_l[3]
    )
{
  int tmin = 0;
  for(int t=0; t<3; t++){
    if(periodicity[t]){
      if (box_l[t] < box_l[tmin])
        tmin = t;
    }
  }

  return tmin;
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
  return (2.0*sum_q2*fcs_exp(-FCS_P2NFFT_SQR(r_cut*alpha))) / (fcs_sqrt((fcs_float)N*r_cut*box_l[0]*box_l[1]*box_l[2]));
}

/* Get the error for this combination of parameters and tune alpha if
 * required. In fact, the real space error is tuned such that it
 * contributes half of the total error, and then the Fourier space
 * error is calculated. Returns the achieved error and the optimal
 * alpha.
 */
static fcs_float p2nfft_tune_alpha(
    fcs_int sum_qpart, fcs_float sum_q2, fcs_int dim_tune,
    fcs_float box_l[3], fcs_float r_cut, ptrdiff_t grid[3],
    fcs_int cao, fcs_float tolerance_field, fcs_int plain_ewald_flag
    )
{
  fcs_float alpha, rs_err; 

  /* use plain Ewald error formulae */
  if (plain_ewald_flag)
    return ewald_tune_alpha(sum_qpart, sum_q2, box_l, r_cut, get_kmax(grid[dim_tune]));
  
  /* calc the maximal real space error for the setting (i.e., set alpha to 0) */
  rs_err = p2nfft_real_space_error(sum_qpart, sum_q2, box_l, r_cut, 0.0);

  if(FCS_P2NFFT_SQRT2 * rs_err > tolerance_field) {
    /* assume rs_err = ks_err -> rs_err = accuracy/fcs_sqrt(2.0) -> alpha */
    alpha = fcs_sqrt(fcs_log(FCS_P2NFFT_SQRT2 * rs_err/tolerance_field)) / r_cut;
  } else {
    /* even alpha=0 is ok, however, we cannot choose it since it kills the k-space error formula.
     * Anyways, this is very likely NOT the optimal solution */
    alpha = 0.1 / box_l[dim_tune];
  }

  return alpha;
}

static fcs_float p2nfft_get_accuracy(
    fcs_int sum_qpart, fcs_float sum_q2, fcs_int dim_tune,
    fcs_float box_l[3], fcs_float r_cut, ptrdiff_t grid[3],
    fcs_int cao, fcs_float tolerance_field,
    fcs_float alpha, fcs_int interlaced, fcs_int plain_ewald_flag,
    fcs_float *rs_err, fcs_float *ks_err
    )
{
  if(plain_ewald_flag){
    /* use plain Ewald error formulae */
    fcs_float dummy_der_r, dummy_der_k;

    ewald_real_space_error(sum_qpart, sum_q2, box_l, r_cut, alpha, rs_err, &dummy_der_r);
    ewald_k_space_error(sum_qpart, sum_q2, box_l, get_kmax(grid[dim_tune]), alpha, ks_err, &dummy_der_k);

    return fcs_sqrt( FCS_P2NFFT_SQR(*rs_err) + FCS_P2NFFT_SQR(*ks_err) );
  }

  /* calculate real space and k space error for this alpha */
  *rs_err = p2nfft_real_space_error(sum_qpart, sum_q2, box_l, r_cut, alpha);

  fcs_int full_estimate = 0;
  for(fcs_int t = 0; t < 3 && !full_estimate; t++){
    fcs_float alpha_h = alpha * box_l[t] / grid[t];
    full_estimate = alpha_h > FCS_P2NFFT_FULL_ESTIMATE_ALPHA_H_THRESHOLD;
  }
  
  if (full_estimate)
    *ks_err = p2nfft_k_space_error(sum_qpart, sum_q2, dim_tune, box_l, grid, alpha, cao, interlaced);
  else
    *ks_err = p2nfft_k_space_error_approx(sum_qpart, sum_q2, dim_tune, box_l, grid, alpha, cao);

  return fcs_sqrt( FCS_P2NFFT_SQR(*rs_err) + FCS_P2NFFT_SQR(*ks_err) );
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
    fcs_int N, fcs_float sum_q2, fcs_int dim_tune,
    fcs_float box_l[3], ptrdiff_t grid[3],
    fcs_float alpha, fcs_int cao,
    fcs_int interlaced
    )
{
  fcs_float he_q = 0.0;
  fcs_float grid_i[3] = {1.0/grid[0], 1.0/grid[1], 1.0/grid[2]};
  /* @todo Handle non-cubic case  */
  fcs_float alpha_L_i = 1./(alpha*box_l[dim_tune]);

#if FCS_P2NFFT_DEBUG_TUNING
  fprintf(stderr, "P2NFFT_DEBUG_TUNING: N = %d, grid_i = [%f, %f, %f], alpha_L_i = %f\n",
      N, grid_i[0], grid_i[1], grid_i[2], alpha_L_i);
#endif

  if(interlaced){
    for (fcs_int nx=-grid[0]/2; nx<grid[0]/2; nx++) {
      for (fcs_int ny=-grid[1]/2; ny<grid[1]/2; ny++) {
        for (fcs_int nz=-grid[2]/2; nz<grid[2]/2; nz++) {
          if((nx!=0) || (ny!=0) || (nz!=0)) {

            fcs_float alias1, alias2, alias3, alias4, alias5, alias6;
            fcs_float d;
    
            p2nfft_k_space_error_sum2_adi(nx,ny,nz,grid,grid_i,cao,alpha_L_i,&alias1,&alias2,&alias3,&alias4,&alias5,&alias6);
            d = alias1 - FCS_P2NFFT_SQR(alias2) / (0.5*(alias3*alias4 + alias5*alias6));
    
            if(d > 0.0)
              he_q += d;
          }
        }
      }
    }
  } else {
    for (fcs_int nx=-grid[0]/2; nx<grid[0]/2; nx++) {
      fcs_float ctan_x = p2nfft_k_space_error_sum1(nx,grid_i[0],cao);
      for (fcs_int ny=-grid[1]/2; ny<grid[1]/2; ny++) {
        fcs_float ctan_y = ctan_x * p2nfft_k_space_error_sum1(ny,grid_i[1],cao);
        for (fcs_int nz=-grid[2]/2; nz<grid[2]/2; nz++) {
          if((nx!=0) || (ny!=0) || (nz!=0)) {
            fcs_float n2 = FCS_P2NFFT_SQR(nx) + FCS_P2NFFT_SQR(ny) + FCS_P2NFFT_SQR(nz);
            fcs_float cs = p2nfft_k_space_error_sum1(nz,grid_i[2],cao)*ctan_y;
            fcs_float alias1, alias2;
            p2nfft_k_space_error_sum2_ad(nx,ny,nz,grid,grid_i,cao,alpha_L_i,&alias1,&alias2);
            fcs_float d = alias1  -  FCS_P2NFFT_SQR(alias2/cs) / n2;
            /* at high precisions, d can become negative due to extinction;
               also, don't take values that have no significant digits left*/
            if (d > 0.0 && (!fcs_float_is_zero(fcs_fabs(d/alias1))))
              he_q += d;
  // if( (nx==1) && (ny==0) && (nz==0) )
  //   fprintf(stderr, "P3M: k = [%d, %d, %d], N = [%td, %td, %td], alias_k = %e, alias1 = %e, alias2 = %e, cs = %e, n2 = %e, alias_minus = %e, ctan_x = %.16e, ctan_y = %e, ctan_z = %e\n",
  //       nx, ny, nz, grid[0], grid[1], grid[2], d, alias1, alias2, cs, n2, FCS_P2NFFT_SQR(alias2/cs) / n2, ctan_x, ctan_y/ctan_x, cs/ctan_y);
          }
        }
      }
    }
  }

#if FCS_P2NFFT_DEBUG_TUNING
  fprintf(stderr, "P2NFFT_DEBUG_TUNING: P3M: he_q = %e\n", he_q);
#endif

#if FCS_P2NFFT_ENABLE_TUNING_BUG
  return 2.0*sum_q2*fcs_sqrt(he_q/(fcs_float)N) / (box_l[1]*box_l[2]);
#else
  if(interlaced){
    return 2.0*sum_q2*fcs_sqrt( he_q / (fcs_float)N ) / fcs_pow(box_l[0]*box_l[1]*box_l[2], 2.0/3.0);
  } else {
    /* Where does the factor 2.0 come from? */
    return 2.0*sum_q2*fcs_sqrt( he_q / (fcs_float)N / (box_l[0]*box_l[1]*box_l[2]) );
  }
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
static void p2nfft_k_space_error_sum2_ad(
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

#undef FCS_P2NFFT_BRILLOUIN
#define FCS_P2NFFT_BRILLOUIN 1 

/** aliasing sum used by \ref ifcs_p3m_k_space_error. */
static void p2nfft_k_space_error_sum2_adi(
    fcs_int nx, fcs_int ny, fcs_int nz,
    ptrdiff_t grid[3], fcs_float grid_i[3],
    fcs_int cao, fcs_float alpha_L_i,
    fcs_float *alias1, fcs_float *alias2,
    fcs_float *alias3, fcs_float *alias4,
    fcs_float *alias5, fcs_float *alias6
    )
{
  fcs_float prefactor = FCS_P2NFFT_SQR(FCS_P2NFFT_PI*alpha_L_i);

  *alias1 = *alias2 = *alias3 = *alias4 = *alias5 = *alias6 = 0.0;
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
	*alias2 += U2 * ex;
	*alias3 += U2 * nm2;
	*alias4 += U2;

        if (((mx+my+mz)%2)==0) {					//even term
	  *alias5 += U2 * nm2;
	  *alias6 += U2;
	} else {						//odd term: minus sign!
	  *alias5 -= U2 * nm2;
	  *alias6 -= U2;
	}
      }
    }
  }
}

#undef FCS_P2NFFT_BRILLOUIN
#define FCS_P2NFFT_BRILLOUIN 0 


/* Calculate the analytical approximation for the k-space part of the
 * error (Eq. 38 in Deserno, Holm; JCP 109,18; 1998). */
static fcs_float p2nfft_k_space_error_approx(
    fcs_int N, fcs_float sum_q2, fcs_int dim_tune,
    fcs_float box_l[3], ptrdiff_t grid[3],
    fcs_float alpha, fcs_int cao
    )
{
  /* grid spacing */
  fcs_float h = box_l[dim_tune]/grid[dim_tune];
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
    fcs_sqrt(alpha*box_l[dim_tune]/N*fcs_sqrt(2.0*FCS_P2NFFT_PI)*sum);
#else
  return 
    sum_q2 / (box_l[dim_tune]*box_l[dim_tune]) *
    fcs_pow(h*alpha, cao) * 
    fcs_sqrt(alpha*box_l[dim_tune]/N*fcs_sqrt(2.0*FCS_P2NFFT_PI)*sum);
#endif
}




#if FCS_P2NFFT_TEST_GENERAL_ERROR_ESTIMATE //hier wird nur mit der 3. Boxlänge gearbeitet?
static fcs_float p2nfft_k_space_error_general_window(
    fcs_int num_part, fcs_float sum_q2,
    const fcs_float box_l[3], const ptrdiff_t N[3],
    fcs_float alpha, fcs_int m,
    unsigned window_flag,
    MPI_Comm comm
    )
{
  /* some dummy variables, since we want to use PNFFT data decomposition */
  ptrdiff_t local_N[3], local_N_start[3];
  fcs_float lower_border[3], upper_border[3];
  fcs_float x_max[3] = {0.5, 0.5, 0.5};
  /* PNFFT calculates with real space cutoff 2*m+2
   * Therefore m is one less than the P3M cao. */
  FCS_PNFFT(plan) pnfft;

  fcs_float alias_k, alias_sum, local_alias_sum = 0.0;
  ptrdiff_t k[3];

  FCS_PNFFT(local_size_guru)(3, N, N, x_max, m, comm, PNFFT_TRANSPOSED_F_HAT,
      local_N, local_N_start, lower_border, upper_border);

  /* Fast initialize of PNFFT, since we want to get the inverse Fourier coefficients.
   * Do not allocate arrays and do not tune PFFT */
  FCS_PNFFT(init_guru)(&pnfft, 3, N, N, x_max, m,
      PNFFT_TRANSPOSED_F_HAT | window_flag, PFFT_ESTIMATE, comm);

  for(k[0]=local_N_start[0]; k[0]<local_N_start[0]+local_N[0]; k[0]++){
    for(k[1]=local_N_start[1]; k[1]<local_N_start[1]+local_N[1]; k[1]++){
      for(k[2]=local_N_start[2]; k[2]<local_N_start[2]+local_N[2]; k[2]++){
        alias_k = compute_alias_k(&pnfft, alpha, box_l, N, k);
        /* at high precisions, d can become negative due to extinction;
           also, don't take values that have no significant digits left*/
        if(alias_k > 0.0)
          local_alias_sum += alias_k;
// if( (k[0]==1) && (k[1]==0) && (k[2]==0) )
//   fprintf(stderr, "P2NFFT: k = [%td, %td, %td], N= [%td, %td, %td], alias_k = %e\n", k[0], k[1], k[2], N[0], N[1], N[2], alias_k);
      }
    }
  }

  FCS_PNFFT(finalize)(&pnfft, 0);

  MPI_Allreduce(&local_alias_sum, &alias_sum, 1, FCS_MPI_FLOAT, MPI_SUM, comm);

  fprintf(stderr, "P2NFFT_DEBUG_TUNING: local_alias_sum = %e\n", local_alias_sum);
  fprintf(stderr, "P2NFFT_DEBUG_TUNING: alias_sum = %e\n", alias_sum);

#if FCS_P2NFFT_ENABLE_TUNING_BUG
  return 2.0*sum_q2*fcs_sqrt(alias_sum/(fcs_float)num_part) / (box_l[1]*box_l[2]);
#else
  /* Where does the factor 2.0 come from? */
  return 2.0*sum_q2*fcs_sqrt( alias_sum / (fcs_float)num_part / (box_l[0]*box_l[1]*box_l[2]) );
#endif
}

static fcs_float compute_alias_k(
    FCS_PNFFT(plan) *wind_param, fcs_float alpha, const fcs_float box_l[3],
    const ptrdiff_t N[3], const ptrdiff_t k[3]
    )
{
  const fcs_int r = 10;
  fcs_int s[3];

  for(fcs_int t=0; t<3; t++)
    s[t] = 4;

  fcs_float alias_k, alias_D2, alias_R2, alias_U2, alias_DU2R;
  fcs_float km0, km1, km2, kk, kkm0, kkm1, kkm2, kmkm0, kmkm1, kmkm2;
  fcs_float Wkm0, Wkm1, Wkm2, Wkm, Ekm0, Ekm1, Ekm2;
  fcs_float sum_Wkm[3], km;
  fcs_float p[3];
 
  if((k[0]==0) && (k[1]==0) && (k[2]==0))
    return 0;

  for(fcs_int t=0; t<3; t++)
    p[t] = FCS_P2NFFT_SQR( FCS_P2NFFT_PI/(alpha*box_l[t]) );

  alias_R2 = alias_DU2R = 0.0;
  for(fcs_int m0=-r; m0<=r; m0++){
    km0 = (fcs_float) k[0]+m0*N[0];
    kkm0 = k[0]*km0;
    kmkm0 = km0*km0;
    Wkm0 = FCS_PNFFT(phi_hat)(wind_param, 0, km0);
    Ekm0 = fcs_exp(-p[0]*km0*km0);
    for(fcs_int m1=-r; m1<=r; m1++){
      km1 = (fcs_float) k[1]+m1*N[1];
      kkm1 = kkm0 + k[1]*km1;
      kmkm1 = kmkm0 + km1*km1;
      Wkm1 = Wkm0*FCS_PNFFT(phi_hat)(wind_param, 1, km1);
      Ekm1 = Ekm0*fcs_exp(-p[1]*km1*km1);
      for(fcs_int m2=-r; m2<=r; m2++){
        km2 = (fcs_float) k[2]+m2*N[2];
        kkm2 = kkm1 + k[2]*km2;
        kmkm2 = (fcs_float) k[2]+m2*N[2]; 
        kmkm2 = kmkm1 + km2*km2;
        Wkm2 = Wkm1*FCS_PNFFT(phi_hat)(wind_param, 2, km2);
        Ekm2 = Ekm1*fcs_exp(-p[2]*km2*km2);

        alias_DU2R += Wkm2*Wkm2 * Ekm2 * kkm2/kmkm2;
        alias_R2   += Ekm2*Ekm2 / kmkm2;
      }
    }
  } 

  for(fcs_int t=0; t<3; t++){
    sum_Wkm[t] = 0.0;
    for(fcs_int m=-s[t]; m<=s[t]; m++){
      km = (fcs_float) k[t]+m*N[t]; 
      Wkm = FCS_PNFFT(phi_hat)(wind_param, t, km);
      sum_Wkm[t] += Wkm*Wkm;
    }
  }

  alias_U2 = 1.0;
  for(fcs_int t=0; t<3; t++)
    alias_U2 *= sum_Wkm[t];  

 if( (k[0]==1) && (k[1]==0) && (k[2]==0) ){
   fcs_int t=0;
   for(fcs_int m=-s[t]; m<=s[t]; m++){
       fcs_int kmi = k[t]+m*N[t];
       Wkm = FCS_PNFFT(phi_hat)(wind_param, t, kmi);
       // fprintf(stderr, "P2NFFT: Sinc^(%d)(Pi*(%f)/(%td)) = %e\n", wind_param->m, km, N[0], Wkm);
       fprintf(stderr, "P2NFFT: Phi_hat(m=%d, N=%td, k=%d) = %e\n", wind_param->m, N[0], kmi, Wkm);
   }
 
   fprintf(stderr, "P2NFFT: k = [%td, %td, %td], sum_Wkm = [%e, %e, %e]\n",
       k[0], k[1], k[2], sum_Wkm[0], sum_Wkm[1], sum_Wkm[2]);
 }

 if( (k[0]==1) && (k[1]==0) && (k[2]==0) ){
   fcs_int t=0;
   for(fcs_int m=-N[t]-5; m<=N[t]+5; m++){
       fcs_int kmi = k[t]+m;
       Wkm = FCS_PNFFT(phi_hat)(wind_param, t, kmi);
       // fprintf(stderr, "P2NFFT: Sinc^(%d)(Pi*(%f)/(%td)) = %e\n", wind_param->m, km, N[0], Wkm);
       fprintf(stderr, "P2NFFT: Check Phi_hat(m=%d, N=%td, k=%d) = %e\n", wind_param->m, N[0], kmi, Wkm);
   }

   for(fcs_int m=-N[t]-5; m<=N[t]+5; m++){
       fcs_int kmi = k[t]+m;
       fprintf(stderr, "P2NFFT: Check DPsi(m=%d, N=%td, x=%f) = %e\n",
           wind_param->m, N[0], (fcs_float)kmi/N[0], FCS_PNFFT(dpsi)(wind_param, t, (fcs_float)kmi/N[0]));
   }
 
   for(fcs_int m=-N[t]-5; m<=N[t]+5; m++){
       fcs_int kmi = k[t]+m;
       fprintf(stderr, "P2NFFT: Check Psi(m=%d, N=%td, x=%f) = %e\n",
           wind_param->m, N[0], (fcs_float)kmi/N[0], FCS_PNFFT(psi)(wind_param, t, (fcs_float)kmi/N[0]));
   }
 fprintf(stderr, "P2NFFT: Check DPsi(m=%d, N=%td, x=%f) = %e\n",
     wind_param->m, N[0], 0.3, FCS_PNFFT(dpsi)(wind_param, 0, 0.3));
 fprintf(stderr, "P2NFFT: Check DPsi(m=%d, N=%td, x=%f) = %e\n",
     wind_param->m, N[0], 0.33, FCS_PNFFT(dpsi)(wind_param, 0, 0.33));
 fprintf(stderr, "P2NFFT: Check DPsi(m=%d, N=%td, x=%f) = %e\n",
     wind_param->m, N[0], 0.374, FCS_PNFFT(dpsi)(wind_param, 0, 0.374));
 }
  
  kk = k[0]*(fcs_float)k[0] + k[1]*(fcs_float)k[1] + k[2]*(fcs_float)k[2];
  alias_D2 = kk;
  
  alias_k = alias_R2 - alias_DU2R*alias_DU2R / (alias_D2*alias_U2*alias_U2);

//   if(alias_k > 1e-2){
//     fprintf(stderr, "k = [%td, %td, %td], alias_R2 = %e, alias_D2 = %e, alias_U2 = %e, alias_DU2R = %e\n",
//         k[0], k[1], k[2], alias_R2, alias_D2, alias_U2, alias_DU2R);
//     fprintf(stderr, "alias_k = %e\n", alias_k);
//   }

// if( (k[1]==0) && (k[2]==0) )
//   fprintf(stderr, "P2NFFT: k = [%td, %td, %td], alias_DU2R = %e, alias_U2 = %.16e, alias_D2 = %e, alias1 = %e, alias_minus = %e\n",
//       k[0], k[1], k[2], alias_DU2R, alias_U2, alias_D2, alias_R2, alias_DU2R*alias_DU2R / (alias_D2*alias_U2*alias_U2));

  return alias_k;
}
#endif


