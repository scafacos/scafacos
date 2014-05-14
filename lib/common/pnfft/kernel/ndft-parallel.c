/*
 * Copyright (c) 2011-2013 Michael Pippig
 *
 * This file is part of PNFFT.
 *
 * PNFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PNFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PNFFT.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <complex.h>
#include "pnfft.h"
#include "ipnfft.h"
#include "bessel_i0.h"
#include "bessel_i1.h"
#include "bspline.h"
#include "sinc.h"
#include "matrix_D.h"

#define PNFFT_ENABLE_CALC_INTPOL_NODES 0
#define USE_EWALD_SPLITTING_FUNCTION_AS_WINDOW 0
#define TUNE_B_FOR_EWALD_SPLITTING 0
#define PNFFT_TUNE_LOOP_ADJ_B 0

static void loop_over_particles_adj(
    PNX(plan) ths, INT *local_no_start, INT *local_ngc, INT *gcells_below,
    INT *sorted_index);
static void loop_over_particles_adj_interlaced_0(
    PNX(plan) ths, INT *local_no_start, INT *local_ngc, INT *gcells_below,
    INT *sorted_index);
static void loop_over_particles_adj_interlaced_1(
    PNX(plan) ths, INT *local_no_start, INT *local_ngc, INT *gcells_below,
    INT *sorted_index);

static PNX(plan) mkplan(
    void);

static void local_size_B(
    const PNX(plan) ths,
    INT *local_no, INT *local_no_start);
static void get_size_gcells(
    int m, int cutoff, unsigned pnfft_flags,
    INT *gcells_below, INT *gcells_above);
static void lowest_summation_index(
    const INT *n, int m, const R *x,
    const INT *local_no_start, const INT *gcells_below,
    R *floor_nx_j, INT *u_j);
static void local_array_size(
    const INT *local_n, const INT *gcells_below, const INT *gcells_above,
    INT *local_ngc);

static void pre_psi_tensor(
    const INT *n, const R *b, int m, int cutoff, const R *x, const R *floor_nx,
    const R *exp_const, R *spline_coeffs, unsigned pnfft_flags,
    int intpol_order, INT intpol_num_nodes, R **intpol_tables_psi,
    R *pre_psi);
static void pre_psi_tensor_direct(
    const INT *n, const R *b, int m, int cutoff, const R *x, const R *floor_nx,
    const R *exp_const, R *spline_coeffs,
    unsigned pnfft_flags,
    R *pre_psi);
static void pre_psi_tensor_gaussian(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx,
    R *pre_psi);
static void pre_psi_tensor_fast_gaussian(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, const R *exp_const,
    R *fg_psi);
static void pre_psi_tensor_bspline(
    const INT *n, int m, int cutoff,
    const R *x, const R *floor_nx, R *spline_coeffs, 
    R *pre_psi);
static void pre_psi_tensor_sinc_power(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, 
    R *pre_psi);
static void pre_psi_tensor_bessel_i0(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, 
    R *pre_psi);
static void pre_psi_tensor_kaiser_bessel(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, 
    R *pre_psi);

static void pre_dpsi_tensor(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, R *spline_coeffs,
    int intpol_order, INT intpol_num_nodes, R **intpol_tables_dpsi,
    const R *pre_psi, unsigned pnfft_flags,
    R *pre_dpsi);
static void pre_dpsi_tensor_direct(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, R *spline_coeffs, const R *pre_psi,
    unsigned pnfft_flags,
    R *pre_dpsi);
static void pre_dpsi_tensor_gaussian(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, const R *fg_psi,
    R *fg_dpsi);
static void pre_dpsi_tensor_bspline(
    const INT *n, int m, int cutoff,
    const R *x, const R *floor_nx, R *spline_coeffs,
    R *pre_dpsi);
static void pre_dpsi_tensor_sinc_power(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, const R *pre_psi,
    R *pre_dpsi);
static void pre_dpsi_tensor_bessel_i0(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, 
    R *pre_dpsi);
static void pre_dpsi_tensor_kaiser_bessel(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, const R *pre_psi,
    R *pre_dpsi);
static void sort_nodes_for_better_cache_handle(
    int d, const INT *n, int m, INT local_x_num, const R *local_x,
    INT *ar_x);
static void project_node_to_grid(
    const INT *n, int m, const R *x,
    R *floor_nx_j, INT *u_j);
static void get_mpi_cart_dims_3d(
    MPI_Comm comm_cart,
    int *rnk_pm, int *dims, int *coords);
static int compare_INT(
    const void* a, const void* b);


static R window_bessel_i0_1d(
    R x, INT n, R b, int m);
static R window_bessel_i0_derivative_1d(
    R x, INT n, R b, int m);

static R kaiser_bessel_1d(
    R x, INT n, R b, int m);
static R kaiser_bessel_derivative_1d(
    R x, INT n, R b, int m, R psi);


static void init_intpol_table_psi(
    INT num_nodes_per_interval, int intpol_order, int cutoff,
    INT n, int m, int dim,
    const PNX(plan) wind_param,
    R *table);
static void init_intpol_table_dpsi(
    INT num_nodes_per_interval, int intpol_order, int cutoff,
    INT n, int m, int dim,
    const PNX(plan) wind_param,
    R *table);
// static void init_intpol_table_psi(
//     INT num_nodes, INT n, int m, int dim,
//     const PNX(plan) wind_param,
//     R *table);
// static void init_intpol_table_dpsi(
//     INT num_nodes, INT n, int m, int dim,
//     const PNX(plan) wind_param,
//     R *table);
static R psi_gaussian(
    R x, INT n, R b);
static R dpsi_gaussian(
    R x, INT n, R b);
static R psi_bspline(
    R x, INT n, int m, R *spline_coeffs);
static R dpsi_bspline(
    R x, INT n, int m, R *spline_coeffs);
static R psi_sinc_power(
    R x, INT n, R b, int m);
static R dpsi_sinc_power(
    R x, INT n, R b, int m);
static R psi_bessel_i0(
    R x, INT n, R b, int m);
static R dpsi_bessel_i0(
    R x, INT n, R b, int m);
static R psi_kaiser(
    R x, INT n, R b, int m);
static R dpsi_kaiser(
    R x, INT n, R b, int m);

static void precompute_psi(
    PNX(plan) ths, INT ind, R* x, R* buffer_psi, R* buffer_dpsi, 
    int compute_grad_ad,
    R* pre_psi, R* pre_dpsi);

/* TODO: This function calculates the number of minimum samples of the 2-point-Taylor regularized
 * kernel function 1/x to reach a certain relative error 'eps'. Our windows are likely to be nicer,
 * so this is a good starting point. Think about the optimal number for every window. */
#if PNFFT_ENABLE_CALC_INTPOL_NODES
static INT max_I(INT a, INT b)
{
  return a >= b ? a : b;
}

static R derivative_bound_guess(
    int intpol_order
    )
{
  /* use interpolation order 3 as default case */
  switch(intpol_order){
    case 0: return 2;
    case 1: return 1.7;
    case 2: return 2.2;
    default: return 1.4;
  }
}


static R get_derivative_bound(
    PNX(plan) ths, intpol_order
    )
{
  /* TODO: implement bound for sinc, kaiser, bessel_i0 and gaussian */
  if(ths->pnfft_flags & PNFFT_WINDOW_BSPLINE)
    return PNX(derivative_bound_bspline)(intpol_order, 2*ths->m);
  else
    return derivative_bound_guess(intpol_order);
}

static INT calc_intpol_num_nodes(
    PNX(plan) ths, int intpol_order, R eps, R r, unsigned *err
    )
{
  R N, c, M, M_pot, M_force;
  *err=0;

  /* define constants from Taylor expansion */
  switch(intpol_order){
    case 0: c = 1.0; break;
    case 1: c = 1.0/8.0; break;
    case 2: c = fcs_sqrt(3)/9.0; break;
    case 3: c = 3.0/128.0; break;
    default: return 0; /* no interpolation */
  }
   
  /* Compute the max. value of the regularization derivative one order higher
   * than the interpolation order. This gives the rest term of the taylor expansion. */ 
  M_pot   = get_derivative_bound(ths, interpolation_order+1);
  M_force = get_derivative_bound(ths, interpolation_order+2);

  /* We use the same number of interpolation nodes for potentials and forces.
   * Be sure, that accuracy is fulfilled for both. */
  M = (M_force > M_pot) ? M_force : M_pot;

  N = r * pnfft_pow(c*M/pnfft_fabs(eps-PNFFT_EPSILON) , 1.0 / (1.0 + intpol_order) ); 

  /* Set maximum number of nodes to avoid overflows during conversion to int. */
  if(N>1e7){
    *err=1;
    N = 1e7;
  }

  /* Return the number of nodes needed for interpolation of the interval [0,1/m] */
  N /= 2*ths->m;
  
  /* At least use 4 interpolation points per interval. */
  if(N<2) N = 2.0;

  /* Compute next power of two >= N to optimize memory access */
  N = pnfft_pow( 2.0, pnfft_ceil(pnfft_log(N)/pnfft_log(2)) );

  return (fcs_int) N;
}
#endif

static void init_intpol_table_psi(
    INT num_nodes_per_interval, int intpol_order, int cutoff,
    INT n, int m, int dim,
    const PNX(plan) wind_param,
    R *table
    )
{
  /* interpolation of "f" at grid point "r" of order
   * 0: uses f[r]
   * 1: uses f[r], f[r+1]
   * 2: uses f[r-1], f[r], f[r+1]
   * 3: uses f[r-1], f[r], f[r+1], f[r+2]
   * This equivalent to f[-order/2], ... , f[(order+1)/2]
   * with integer division. */
  INT ind=0;
  for(INT k=0; k<num_nodes_per_interval; k++)
    for(INT c=0; c<cutoff; c++)
      for(INT i=-intpol_order/2; i<=(intpol_order+1)/2; i++)
        table[ind++] = PNX(psi)(wind_param, dim, (-m + (R)(-k+i)/num_nodes_per_interval + c)/n);
}

static void init_intpol_table_dpsi(
    INT num_nodes_per_interval, int intpol_order, int cutoff,
    INT n, int m, int dim,
    const PNX(plan) wind_param,
    R *table
    )
{
  /* interpolation of "f" at grid point "r" of order
   * 0: uses f[r]
   * 1: uses f[r], f[r+1]
   * 2: uses f[r-1], f[r], f[r+1]
   * 3: uses f[r-1], f[r], f[r+1], f[r+2]
   * This equivalent to f[-order/2], ... , f[(order+1)/2]
   * with integer division. */
  INT ind=0;
  for(INT k=0; k<num_nodes_per_interval; k++)
    for(INT c=0; c<cutoff; c++)
      for(INT i=-intpol_order/2; i<=(intpol_order+1)/2; i++)
        table[ind++] = PNX(dpsi)(wind_param, dim, (-m   + (R)(-k+i)/num_nodes_per_interval + c)/n);
}



void PNX(trafo_A)(
    PNX(plan) ths
    )
{
  int np_total, myrnk;
  INT local_Np[3], local_Np_start[3]; 
  C *buffer;

  MPI_Comm_size(ths->comm_cart, &np_total);
  MPI_Comm_rank(ths->comm_cart, &myrnk);

  if (ths->trafo_flag & PNFFTI_TRAFO_C2R) {
    for(INT j=0; j<ths->local_M; j++)  ths->f[j] = 0;
    if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F)
      for(INT j=0; j<3*ths->local_M; j++)  ths->grad_f[j] = 0;
  } else if (ths->trafo_flag & PNFFTI_TRAFO_C2C) {
    for(INT j=0; j<ths->local_M; j++)  ((C*)ths->f)[j] = 0;
    if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F)
      for(INT j=0; j<3*ths->local_M; j++)  ((C*)ths->grad_f)[j] = 0;
  }

  for(int pid=0; pid<np_total; pid++){
    /* compute local_Np, local_Np_start of proc. with rank pid */
    PNX(local_block_internal)(ths->N, ths->no, ths->comm_cart, pid, ths->pnfft_flags, ths->trafo_flag,
        local_Np, local_Np_start);

    INT local_Np_total = PNX(prod_INT)(3, local_Np);

    /* Avoid errors for empty blocks */
    if(local_Np_total == 0) continue;

    buffer = (myrnk == pid) ? ths->f_hat : PNX(malloc_C)(local_Np_total);

    /* broadcast block of Fourier coefficients from p to all procs */
    MPI_Bcast(buffer, 2*local_Np_total, PNFFT_MPI_REAL_TYPE, pid, ths->comm_cart);

    INT t0 = (ths->pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? 1 : 0;
    INT t1 = (ths->pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? 2 : 1;
    INT t2 = (ths->pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? 0 : 2;

    for(INT j=0; j<ths->local_M; j++){
      C exp_x0 = pnfft_cexp(-2.0 * PNFFT_PI * ths->x[3*j+t0] * I);
      C exp_x1 = pnfft_cexp(-2.0 * PNFFT_PI * ths->x[3*j+t1] * I);
      C exp_x2 = pnfft_cexp(-2.0 * PNFFT_PI * ths->x[3*j+t2] * I);

      C exp_kx0_start = pnfft_cexp(-2.0 * PNFFT_PI * local_Np_start[t0] * ths->x[3*j+t0] * I);
      C exp_kx1_start = pnfft_cexp(-2.0 * PNFFT_PI * local_Np_start[t1] * ths->x[3*j+t1] * I);
      C exp_kx2_start = pnfft_cexp(-2.0 * PNFFT_PI * local_Np_start[t2] * ths->x[3*j+t2] * I);

      if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F){
        R grad_f[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        INT m=0;
        C exp_kx0 = exp_kx0_start;
        for(INT k0 = local_Np_start[t0]; k0 < local_Np_start[t0] + local_Np[t0]; k0++){
          C exp_kx1 = exp_kx0 * exp_kx1_start;
          for(INT k1 = local_Np_start[t1]; k1 < local_Np_start[t1] + local_Np[t1]; k1++){
            C exp_kx2 = exp_kx1 * exp_kx2_start;
            for(INT k2 = local_Np_start[t2]; k2 < local_Np_start[t2] + local_Np[t2]; k2++, m++){
              C bufferTimesExp = buffer[m] * exp_kx2;

              if (ths->trafo_flag & PNFFTI_TRAFO_C2R) {
                if (k0 == 0 && k1 == 0 && k2 == 0)
                  ths->f[j] += pnfft_creal(buffer[m]);
                else if ( ! (// these have to be zero
                             ( (k0 == 0 || k0 ==  -ths->N[t0]/2) &&
                               (k1 == 0 || k1 ==  -ths->N[t1]/2) &&
                               (k2 == 0 || k2 ==  -ths->N[t2]/2)    ) ||
                             // these are redundant. we have to skip them because we always add the two
                             // hermitean coefficients at once and we would otherwise add them twice
                             ( k1 > 0 && (k2 == 0 || k2 == -ths->N[t2]/2) ) ||
                             ( k0 > 0 && (k1 == 0 || k1 == -ths->N[t1]/2) && (k2 == 0 || k2 == -ths->N[t2]/2) ))
                        ) {
                  ths->f[j] += 2 * pnfft_creal(bufferTimesExp);
                  grad_f[0] += k0 * pnfft_cimag(bufferTimesExp);
                  grad_f[2] += k1 * pnfft_cimag(bufferTimesExp);
                  grad_f[4] += k2 * pnfft_cimag(bufferTimesExp);
                }
              } else {
                ((C*)ths->f)[j] += bufferTimesExp;
                ((C*)grad_f)[0] += k0 * bufferTimesExp;
                ((C*)grad_f)[1] += k1 * bufferTimesExp;
                ((C*)grad_f)[2] += k2 * bufferTimesExp;
              }

              exp_kx2 *= exp_x2;
            }
            exp_kx1 *= exp_x1;
          }
          exp_kx0 *= exp_x0;
        }

        if (ths->trafo_flag & PNFFTI_TRAFO_C2R) {
          ths->grad_f[3*j+t0] += 2 * grad_f[0];
          ths->grad_f[3*j+t1] += 2 * grad_f[2];
          ths->grad_f[3*j+t2] += 2 * grad_f[4];
        } else if (ths->trafo_flag & PNFFTI_TRAFO_C2C) {
          ((C*)ths->grad_f)[3*j+t0] += grad_f[0];
          ((C*)ths->grad_f)[3*j+t1] += grad_f[1];
          ((C*)ths->grad_f)[3*j+t2] += grad_f[2];
        }

      } else {
        INT m=0;
        C exp_kx0 = exp_kx0_start;
        for(INT k0 = local_Np_start[t0]; k0 < local_Np_start[t0] + local_Np[t0]; k0++){
          C exp_kx1 = exp_kx0 * exp_kx1_start;
          for(INT k1 = local_Np_start[t1]; k1 < local_Np_start[t1] + local_Np[t1]; k1++){
            C exp_kx2 = exp_kx1 * exp_kx2_start;
            for(INT k2 = local_Np_start[t2]; k2 < local_Np_start[t2] + local_Np[t2]; k2++, m++){
              if (ths->trafo_flag & PNFFTI_TRAFO_C2R) {
                if (k0 == 0 && k1 == 0 && k2 == 0)
                  ths->f[j] += pnfft_creal(buffer[m]);
                else if ( ! (// these have to be zero
                             ( (k0 == 0 || k0 ==  -ths->N[t0]/2) &&
                               (k1 == 0 || k1 ==  -ths->N[t1]/2) &&
                               (k2 == 0 || k2 ==  -ths->N[t2]/2)    ) ||
                             // these are redundant. we have to skip them because we always add the two
                             // hermitean coefficients at once and we would otherwise add them twice
                             ( k1 > 0 && (k2 == 0 || k2 == -ths->N[t2]/2) ) ||
                             ( k0 > 0 && (k1 == 0 || k1 == -ths->N[t1]/2) && (k2 == 0 || k2 == -ths->N[t2]/2) ))
                        )
                  ths->f[j] += 2 * pnfft_creal(buffer[m] * exp_kx2);
              } else if (ths->trafo_flag & PNFFTI_TRAFO_C2C)
                ((C*)ths->f)[j] += buffer[m] * exp_kx2;

              exp_kx2 *= exp_x2;
            }
            exp_kx1 *= exp_x1;
          }
          exp_kx0 *= exp_x0;
        }
      }
    }

    if(myrnk != pid) PNX(free)(buffer);
  }

  C minusTwoPi  = -2.0 * PNFFT_PI;
  C minusTwoPiI = minusTwoPi * I;
  if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F) {
    if (ths->trafo_flag & PNFFTI_TRAFO_C2R)
      for(INT j=0; j<3*ths->local_M; j++)
        ths->grad_f[j] *= minusTwoPi;
    else if (ths->trafo_flag & PNFFTI_TRAFO_C2C)
      for(INT j=0; j<3*ths->local_M; j++)
        ((C*)ths->grad_f)[j] *= minusTwoPiI;
  }
}



void PNX(adj_A)(
    PNX(plan) ths
    )
{
  int np_total, myrnk;
  INT local_Np[3], local_Np_start[3]; 
  C *buffer;

  MPI_Comm_size(ths->comm_cart, &np_total);
  MPI_Comm_rank(ths->comm_cart, &myrnk);

  for(int pid=0; pid<np_total; pid++){
    /* compute local_Np, local_Np_start of proc. with rank pid */
    PNX(local_block_internal)(ths->N, ths->no, ths->comm_cart, pid, ths->pnfft_flags, ths->trafo_flag,
        local_Np, local_Np_start);

    INT local_Np_total = PNX(prod_INT)(3, local_Np);

    /* Avoid errors for empty blocks */
    if(local_Np_total == 0) continue;

    buffer = malloc(sizeof(C) * local_Np_total);
    for(INT k=0; k<local_Np_total; k++) buffer[k] = 0;

    INT t0 = (ths->pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? 1 : 0;
    INT t1 = (ths->pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? 2 : 1;
    INT t2 = (ths->pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? 0 : 2;

    for(INT j=0; j<ths->local_M; j++){
      C exp_x0 = pnfft_cexp(+2.0 * PNFFT_PI * ths->x[3*j+t0] * I);
      C exp_x1 = pnfft_cexp(+2.0 * PNFFT_PI * ths->x[3*j+t1] * I);
      C exp_x2 = pnfft_cexp(+2.0 * PNFFT_PI * ths->x[3*j+t2] * I);

      C exp_kx0_start = pnfft_cexp(+2.0 * PNFFT_PI * local_Np_start[t0] * ths->x[3*j+t0] * I);
      C exp_kx1_start = pnfft_cexp(+2.0 * PNFFT_PI * local_Np_start[t1] * ths->x[3*j+t1] * I);
      C exp_kx2_start = pnfft_cexp(+2.0 * PNFFT_PI * local_Np_start[t2] * ths->x[3*j+t2] * I);

      INT m=0;
      C exp_kx0;
      if (ths->trafo_flag & PNFFTI_TRAFO_C2R)
        exp_kx0 = exp_kx0_start * ths->f[j];
      else if (ths->trafo_flag & PNFFTI_TRAFO_C2C)
        exp_kx0 = exp_kx0_start * ((C*)ths->f)[j];
      for(INT k0 = local_Np_start[t0]; k0 < local_Np_start[t0] + local_Np[t0]; k0++){
        C exp_kx1 = exp_kx0 * exp_kx1_start;
        for(INT k1 = local_Np_start[t1]; k1 < local_Np_start[t1] + local_Np[t1]; k1++){
          C exp_kx2 = exp_kx1 * exp_kx2_start;
          for(INT k2 = local_Np_start[t2]; k2 < local_Np_start[t2] + local_Np[t2]; k2++, m++){
            buffer[m] += exp_kx2;

            exp_kx2 *= exp_x2;
          }
          exp_kx1 *= exp_x1;
        }
        exp_kx0 *= exp_x0;
      }
    }

    /* reduce block of Fourier coefficients from all procs to p */
    MPI_Reduce(buffer, ths->f_hat, 2*local_Np_total, PNFFT_MPI_REAL_TYPE, MPI_SUM, pid, ths->comm_cart);
    PNX(free)(buffer);
  }
}



/* The local borders of the FFT output imply the box of nodes
 * that can be handled by the current process. If x_max < 0.5,
 * the FFT output size 'no' includes some extra elements at the
 * beginning and the end (the number of extra elements depends on
 * the interpolation order 'm').
 * Therefore, we must explicitly check, if the borders are not smaller
 * then -x_max and not bigger then x_max.
 */
void PNX(node_borders)(
    const INT *n,
    const INT* local_no, const INT *local_no_start,
    const R* x_max,
    R *lo, R *up
    )
{
  for(int t=0; t<3; t++){
    lo[t] = local_no_start[t] / ((R) n[t]);    
    up[t] = (local_no_start[t] + local_no[t]) / ((R) n[t]);

//     if(pnfft_flags & PNFFT_SHIFTED_OUT){
      /* assure, that borders lie in [-x_max,x_max] */
      if(lo[t] < -x_max[t])
        lo[t] = -x_max[t];
      if(lo[t] > x_max[t])
        lo[t] = x_max[t];
      if(up[t] < -x_max[t])
        up[t] = -x_max[t];
      if(up[t] > x_max[t])
        up[t] = x_max[t];
//     } else {
//       /* assure, that borders lie in [0, x_max] or [1-x_max,1] */
//       if(lo[t] < 0.5){
//         if(lo[t] > x_max[t])
//           lo[t] = x_max[t];
// 
//       if(lo[t] < -x_max[t])
//         lo[t] = -x_max[t];
//       if(lo[t] > x_max[t])
//         lo[t] = x_max[t];
//       if(up[t] < -x_max[t])
//         up[t] = -x_max[t];
//       if(up[t] > x_max[t])
//         up[t] = x_max[t];
//     }
// 
// 
// 
// 
  }
}

static void local_size_B(
    const PNX(plan) ths,
    INT *local_no, INT *local_no_start
    )
{
  for(int t=0; t<ths->d; t++){
    local_no[t] = ths->local_no[t];
    local_no_start[t] = ths->local_no_start[t];
  }
}

INT PNX(local_size_internal)(
    const INT *N, const INT *n, const INT *no,
    MPI_Comm comm_cart,
    unsigned trafo_flag, unsigned pnfft_flags,
    INT *local_N, INT *local_N_start,
    INT *local_no, INT *local_no_start
    )
{
  INT howmany = 1;
  unsigned pfft_flags;

  if (trafo_flag & PNFFTI_TRAFO_C2R) {
    INT alloc_local_data_forw, alloc_local_data_back;
    pfft_flags = (pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? PFFT_TRANSPOSED_IN : 0;

    alloc_local_data_forw = PX(local_size_many_dft_c2r)(3, n, N, no, howmany,
        PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pfft_flags | PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT,
        local_N, local_N_start, local_no, local_no_start);

    pfft_flags = (pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? PFFT_TRANSPOSED_OUT : 0;

    alloc_local_data_back = PX(local_size_many_dft_r2c)(3, n, no, N, howmany,
        PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pfft_flags | PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT,
        local_no, local_no_start, local_N, local_N_start);

    return (alloc_local_data_forw > alloc_local_data_back) ?
        alloc_local_data_forw : alloc_local_data_back;
  } else if (trafo_flag & PNFFTI_TRAFO_C2C) {
    pfft_flags = (pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? PFFT_TRANSPOSED_IN : 0;

    return PX(local_size_many_dft)(3, n, N, no, howmany,
        PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pfft_flags | PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT,
        local_N, local_N_start, local_no, local_no_start);
  }
}

void PNX(local_block_internal)(
    const INT *N, const INT *no,
    MPI_Comm comm_cart, int pid,
    unsigned pnfft_flags, unsigned trafo_flag,
    INT *local_N, INT *local_N_start
    )
{
  unsigned pfft_flags = (pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? PFFT_TRANSPOSED_IN : 0;
  INT dummy_lno[3], dummy_los[3];

  if (trafo_flag & PNFFTI_TRAFO_C2R) {
    PX(local_block_many_dft_c2r)(3, N, no,
        PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pid, pfft_flags | PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT,
        local_N, local_N_start, dummy_lno, dummy_los);
  } else if (trafo_flag & PNFFTI_TRAFO_C2C) {
    PX(local_block_many_dft)(3, N, no,
        PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pid, pfft_flags | PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT,
        local_N, local_N_start, dummy_lno, dummy_los);
  }
}


/* N - size of NFFT
 * n - oversampled FFT size
 * no - FFT output size (if nodes are only in a subset the array) */
PNX(plan) PNX(init_internal)(
    int d, const INT *N, const INT *n, const INT *no,
    INT local_M, int m,
    unsigned trafo_flag, unsigned pnfft_flags, unsigned pfft_opt_flags,
    MPI_Comm comm_cart
    )
{
  unsigned pfft_flags=0;
  INT howmany = 1;
  INT alloc_local_data_in, alloc_local_gc;
  INT gcells_below[3], gcells_above[3];
  INT local_ngc[3], local_gc_start[3];
  PNX(plan) ths;

  /* TODO: apply some parameter checks */
  
  if(pnfft_flags & PNFFT_GRAD_IK && pnfft_flags & PNFFT_GRAD_NONE){
    PX(printf)(comm_cart, "!!! Warning: GRAD_IK and GRAD_NONE can not be used together. Using GRAD_IK for this plan !!!\n");
    pnfft_flags &=  (~PNFFT_GRAD_NONE);
  }

  if(pnfft_flags & PNFFT_PRE_PSI && pnfft_flags & PNFFT_PRE_FULL_PSI){
    PX(printf)(comm_cart, "!!! Warning: PRE_PSI and PRE_FULL_PSI can not be used together. Using PRE_PSI for this plan !!!\n");
    pnfft_flags &= (~PNFFT_PRE_FULL_PSI); /* needed for correct pnfft_finalize */
  }

  if((pnfft_flags & PNFFT_TRANSPOSED_F_HAT) && (trafo_flag & PNFFTI_TRAFO_C2R)) {
    PX(printf)(comm_cart, "!!! Error: transposed input and c2r transform not yet implemented !!!\n");
    return NULL;
  }

  ths = mkplan();

  ths->d = d;
  ths->m= m;
  ths->local_M = local_M;

  ths->N = (INT*) PNX(malloc)(sizeof(INT) * (size_t) d);
  ths->n = (INT*) PNX(malloc)(sizeof(INT) * (size_t) d);
  ths->no= (INT*) PNX(malloc)(sizeof(INT) * (size_t) d);
  for(int t=0; t<d; t++){
    ths->N[t]= N[t];
    ths->n[t]= n[t];
    ths->no[t]= no[t];
  }

  ths->local_N        = (INT*) PNX(malloc)(sizeof(INT) * (size_t) d);
  ths->local_N_start  = (INT*) PNX(malloc)(sizeof(INT) * (size_t) d);
  ths->local_no       = (INT*) PNX(malloc)(sizeof(INT) * (size_t) d);
  ths->local_no_start = (INT*) PNX(malloc)(sizeof(INT) * (size_t) d);

  ths->pnfft_flags = pnfft_flags;
  ths->pfft_opt_flags = pfft_opt_flags;
  ths->trafo_flag = trafo_flag;

  MPI_Comm_dup(comm_cart, &(ths->comm_cart));
  get_mpi_cart_dims_3d(comm_cart, &ths->rnk_pm, ths->np, ths->coords);
  
  ths->cutoff = 2*m+1;
  ths->N_total = ths->n_total = 1;
  for(int t=0; t<d; t++){
    ths->N_total *= N[t];
    ths->n_total *= n[t];
  }
  /* x_max is filled in init_guru */
  ths->x_max = (R*) PNX(malloc)(sizeof(R) * (size_t) d);
  ths->sigma = (R*) PNX(malloc)(sizeof(R) * (size_t) d);
  for(int t = 0;t < d; t++)
    ths->sigma[t] = ((R)n[t])/N[t];

  /* init window specific parameters */
  ths->b = (R*) PNX(malloc)(sizeof(R) * (size_t) d);
  if(pnfft_flags & PNFFT_WINDOW_GAUSSIAN){
    for(int t=0; t<d; t++)
      ths->b[t]= ((R)m / PNFFT_PI) * K(2.0)*ths->sigma[t] / (K(2.0)*ths->sigma[t]-K(1.0));
#if TUNE_B_FOR_EWALD_SPLITTING
    for(int t=0; t<d; t++)
      ths->b[t]= 0.715303;
#endif
#if USE_EWALD_SPLITTING_FUNCTION_AS_WINDOW
    for(int t=0; t<d; t++){
      R C=0.976; /* strange constant from paper by D. Lindbo */
//       R B=10.0; /* box length */
//       R alpha=0.573; /* tuned Ewald splitting parameter */
//       R tmp = (R) ths->n[t] / alpha / B;
//       ths->b[t] -= tmp*tmp;
//       R tmp = (R) ths->n[t];
      ths->b[t] = 2.0 * (R)ths->m / (PNFFT_PI * C*C);
    }
//     for(int t=0; t<d; t++){
//       ths->b[t] = 1.2732;
//     }
#endif
    if(pnfft_flags & PNFFT_FG_PSI){
      ths->exp_const = (R*) PNX(malloc)(sizeof(R) * (size_t) d * ths->cutoff);
      for(int t=0; t<d; t++)
        for(int s=0; s<ths->cutoff; s++)
          ths->exp_const[ths->cutoff*t+s] = pnfft_exp(-s*s/ths->b[t])/(pnfft_sqrt(PNFFT_PI*ths->b[t]));
    }
  } else if(pnfft_flags & PNFFT_WINDOW_BSPLINE){
    /* malloc array for scratch values of de Boor algorithm, no need to initialize */
    ths->spline_coeffs= (R*) PNX(malloc)(sizeof(R)*2*ths->m);
  } else if(pnfft_flags & PNFFT_WINDOW_SINC_POWER){
    /* malloc array for scratch values of de Boor algorithm, no need to initialize */
    ths->spline_coeffs= (R*) PNX(malloc)(sizeof(R)*2*ths->m);
    for(int t=0; t<d; t++)
      ths->b[t]= (R)m * (K(2.0)*ths->sigma[t]) / (K(2.0)*ths->sigma[t]-K(1.0));
#if TUNE_B_FOR_EWALD_SPLITTING
//     fprintf(stderr, "Sinc-Power: old b = %.4e\n", ths->b[0]);
    for(int t=0; t<d; t++)
      ths->b[t]= 2.25;
//     fprintf(stderr, "Sinc-Power: new b = %.4e\n", ths->b[0]);
#endif
  } else if(pnfft_flags & PNFFT_WINDOW_BESSEL_I0){
    for(int t=0; t<d; t++)
      ths->b[t] = (R) PNFFT_PI * (K(2.0) - K(1.0)/ths->sigma[t]);
#if TUNE_B_FOR_EWALD_SPLITTING
    for(int t=0; t<d; t++)
      ths->b[t]= 5.45066;
#endif
  } else { /* default window function is Kaiser-Bessel */
    for(int t=0; t<d; t++)
      ths->b[t] = (R) PNFFT_PI * (K(2.0) - K(1.0)/ths->sigma[t]);
#if TUNE_B_FOR_EWALD_SPLITTING
    for(int t=0; t<d; t++)
//       ths->b[t]= 5.35;
      ths->b[t]= 5.7177;
#endif
  }

  get_size_gcells(m, ths->cutoff, pnfft_flags,
      gcells_below, gcells_above);

  alloc_local_data_in = PNX(local_size_internal)(N, n, no, comm_cart, ths->trafo_flag, ths->pnfft_flags,
      ths->local_N, ths->local_N_start, ths->local_no, ths->local_no_start);

  if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
    alloc_local_gc = PX(local_size_many_gc)(3, ths->local_no, ths->local_no_start,
        2 * alloc_local_data_in, howmany, gcells_below, gcells_above,
        local_ngc, local_gc_start);
  else
    alloc_local_gc = 2 * PX(local_size_many_gc)(3, ths->local_no, ths->local_no_start,
        alloc_local_data_in, howmany, gcells_below, gcells_above,
        local_ngc, local_gc_start);

  ths->local_N_total  = PNX(prod_INT)(d, ths->local_N);
  ths->local_no_total = PNX(prod_INT)(d, ths->local_no);

  if(pnfft_flags & PNFFT_MALLOC_F_HAT)
    ths->f_hat = (ths->local_N_total) ? (C*) PNX(malloc)(sizeof(C) * (size_t) ths->local_N_total) : NULL;

  /* init PFFT all the time (do not use the PNFFT_INIT_FFT flag anymore since
   * the init of parallel FFT is far too complicated for any user) */
  ths->g2 = (alloc_local_gc) ? (R*) PNX(malloc)(sizeof(R) * (size_t) alloc_local_gc) : NULL;
  if(pnfft_flags & PNFFT_FFT_IN_PLACE)
    ths->g1 = ths->g2;
  else
    ths->g1 = (alloc_local_data_in) ? (R*) PNX(malloc)(sizeof(R) * 2 * (size_t) alloc_local_data_in) : NULL;

  /* For derivative in Fourier space we need an extra buffer
   * (since we need to scale the output of the forward FFT with three different factors) */
  if(ths->pnfft_flags & PNFFT_GRAD_IK)
    ths->g1_buffer = (ths->local_N_total) ? (R*) PNX(malloc)(sizeof(R) * 2 * (size_t) ths->local_N_total) : NULL;
  else
    ths->g1_buffer = NULL;

  /* plan PFFT */
  pfft_flags = pfft_opt_flags | PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT;
  if(ths->pnfft_flags & PNFFT_TRANSPOSED_F_HAT)
    pfft_flags |= PFFT_TRANSPOSED_IN;
  if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
    ths->pfft_forw = PX(plan_many_dft_c2r)(3, n, N, no, howmany,
        PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, ths->g1, ths->g2, comm_cart,
        PFFT_FORWARD, pfft_flags);
  else
    ths->pfft_forw = PX(plan_many_dft)(3, n, N, no, howmany,
        PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, ths->g1, ths->g2, comm_cart,
        PFFT_FORWARD, pfft_flags);
  
  pfft_flags = pfft_opt_flags | PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT;
  if(ths->pnfft_flags & PNFFT_TRANSPOSED_F_HAT) 
    pfft_flags |= PFFT_TRANSPOSED_OUT;
  if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
    ths->pfft_back = PX(plan_many_dft_r2c)(3, n, no, N, howmany,
        PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, ths->g2, ths->g1, comm_cart,
        PFFT_BACKWARD, pfft_flags);
  else
    ths->pfft_back = PX(plan_many_dft)(3, n, no, N, howmany,
        PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, ths->g2, ths->g1, comm_cart,
        PFFT_BACKWARD, pfft_flags);

  /* plan ghost cell send and receive */
  if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
    ths->gcplan = PX(plan_many_rgc)(3, no, howmany, PFFT_DEFAULT_BLOCKS,
        gcells_below, gcells_above, ths->g2, comm_cart, 0);
  else
    ths->gcplan = PX(plan_many_cgc)(3, no, howmany, PFFT_DEFAULT_BLOCKS,
        gcells_below, gcells_above, ths->g2, comm_cart, 0);

  /* precompute deconvultion in Fourier space */
  if(pnfft_flags & PNFFT_PRE_PHI_HAT){
    ths->pre_inv_phi_hat_trafo = (C*) malloc(sizeof(C) * PNX(sum_INT)(d, ths->local_N));
    ths->pre_inv_phi_hat_adj   = (C*) malloc(sizeof(C) * PNX(sum_INT)(d, ths->local_N));
    
    PNX(precompute_inv_phi_hat_trafo)(ths,
        ths->pre_inv_phi_hat_trafo);
    PNX(precompute_inv_phi_hat_adj)(ths,
        ths->pre_inv_phi_hat_adj);
  }

  /* init interpolation of window function */
  if(pnfft_flags & PNFFT_PRE_CONST_PSI)
    ths->intpol_order = 0;
  else if(pnfft_flags & PNFFT_PRE_LIN_PSI)
    ths->intpol_order = 1;
  else if(pnfft_flags & PNFFT_PRE_QUAD_PSI)
    ths->intpol_order = 2;
  else if(pnfft_flags & PNFFT_PRE_KUB_PSI)
    ths->intpol_order = 3;
  else
    ths->intpol_order = -1;

  if(pnfft_flags & PNFFT_PRE_INTPOL_PSI){
#if PNFFT_ENABLE_CALC_INTPOL_NODES
    ths->intpol_num_nodes = calc_intpol_num_nodes(ths->intpol_order, 1e-16);
#else
    ths->intpol_num_nodes = 32768;
#endif
    ths->intpol_tables_psi = (R**) PNX(malloc)(sizeof(R*) * (size_t) d);
    for(int t=0; t<d; t++){
      ths->intpol_tables_psi[t] = (R*) PNX(malloc)(sizeof(R)*(ths->intpol_num_nodes * ths->cutoff * (ths->intpol_order+1)));
      init_intpol_table_psi(ths->intpol_num_nodes, ths->intpol_order, ths->cutoff, ths->n[t], ths->m, t, ths,
          ths->intpol_tables_psi[t]);
    }
    if( !(pnfft_flags & (PNFFT_GRAD_IK | PNFFT_GRAD_NONE) ) ){
      ths->intpol_tables_dpsi = (R**) PNX(malloc)(sizeof(R*) * (size_t) ths->d);
      for(int t=0; t<d; t++){
        ths->intpol_tables_dpsi[t] = (R*) PNX(malloc)(sizeof(R)*(ths->intpol_num_nodes * ths->cutoff * (ths->intpol_order+1)));
        init_intpol_table_dpsi(ths->intpol_num_nodes, ths->intpol_order, ths->cutoff, ths->n[t], ths->m, t, ths,
            ths->intpol_tables_dpsi[t]);
      }
    }
  }

  /* init compute flags before malloc f and f_grad */
  ths->compute_flags = 0;

  /* allocate mem and adjust pnfft_flags, compute_flags */
  PNX(malloc_x)(ths, pnfft_flags);
  PNX(malloc_f)(ths, pnfft_flags);
  PNX(malloc_grad_f)(ths, pnfft_flags);

  return ths;
}

/* x and local_M must be initialized */
void PNX(precompute_psi)(
    PNX(plan) ths
    )
{
  INT size;
  INT *sorted_index = NULL;
  R *buffer_psi=NULL, *buffer_dpsi=NULL;
  R x[3];
  int compute_grad_ad;
 
  compute_grad_ad = !(ths->pnfft_flags & (PNFFT_GRAD_NONE | PNFFT_GRAD_IK))
                    && (ths->compute_flags & PNFFT_COMPUTE_GRAD_F);

  /* cleanup old precomputations */
  if(ths->pre_psi != NULL)
    PNX(free)(ths->pre_psi);
  if(ths->pre_dpsi != NULL)
    PNX(free)(ths->pre_dpsi);
  if(ths->pre_psi_il != NULL)
    PNX(free)(ths->pre_psi_il);
  if(ths->pre_dpsi_il != NULL)
    PNX(free)(ths->pre_dpsi_il);

  if(!(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) )
    return;

  /* allocate memory */
  if(ths->pnfft_flags & PNFFT_PRE_PSI){
    size = 3 * ths->cutoff * ths->local_M;
    ths->pre_psi = (size) ? (R*) PNX(malloc)(sizeof(R) * size) : NULL;
    if(ths->pnfft_flags & PNFFT_INTERLACED)
      ths->pre_psi_il = (size) ? (R*) PNX(malloc)(sizeof(R) * size) : NULL;
    if(compute_grad_ad){
      ths->pre_dpsi = (size) ? (R*) PNX(malloc)(sizeof(R) * size) : NULL;
      if(ths->pnfft_flags & PNFFT_INTERLACED)
        ths->pre_dpsi_il = (size) ? (R*) PNX(malloc)(sizeof(R) * size) : NULL;
    }
  }

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI){
    size = PNFFT_POW3(ths->cutoff) * ths->local_M;
    ths->pre_psi = (size) ? (R*) PNX(malloc)(sizeof(R) * size) : NULL;
    if(ths->pnfft_flags & PNFFT_INTERLACED)
      ths->pre_psi_il = (size) ? (R*) PNX(malloc)(sizeof(R) * size) : NULL;
    if(compute_grad_ad){
      ths->pre_dpsi = (size) ? (R*) PNX(malloc)(sizeof(R) * 3 * size) : NULL;
      if(ths->pnfft_flags & PNFFT_INTERLACED)
        ths->pre_dpsi_il = (size) ? (R*) PNX(malloc)(sizeof(R) * 3 * size) : NULL;
    }
  }

  /* save precomputations in the same order as need in matrix B */
  if(ths->pnfft_flags & PNFFT_SORT_NODES){
    sorted_index = (INT*) PNX(malloc)(sizeof(INT) * (size_t) 2*ths->local_M);
    sort_nodes_for_better_cache_handle(
        ths->d, ths->n, ths->m, ths->local_M, ths->x,
        sorted_index);
  }

  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI){
    buffer_psi = (R*) PNX(malloc)(sizeof(R) * (size_t) ths->cutoff*3);
    if(compute_grad_ad)
      buffer_dpsi = (R*) PNX(malloc)(sizeof(R) * (size_t) ths->cutoff*3);
  }
  for(INT p=0; p<ths->local_M; p++){
    INT j = (ths->pnfft_flags & PNFFT_SORT_NODES) ? sorted_index[2*p+1] : p;

    for(int t=0; t<3; t++)
      x[t] = ths->x[ths->d*j+t];
    precompute_psi(ths, p, x, buffer_psi, buffer_dpsi, compute_grad_ad,
        ths->pre_psi, ths->pre_dpsi);

    if(ths->pnfft_flags & PNFFT_INTERLACED){
      /* shift x by half the mesh width */
      for(int t=0; t<3; t++){
        x[t] = ths->x[ths->d*j+t] + 0.5/ths->n[t];
        if(x[t] >= 0.5)
          x[t] -= 1.0;
      }
      precompute_psi(ths, p, x, buffer_psi, buffer_dpsi, compute_grad_ad,
          ths->pre_psi_il, ths->pre_dpsi_il);
    }
  }

  if(sorted_index != NULL)
    PNX(free)(sorted_index);
  if(buffer_psi != NULL)
    PNX(free)(buffer_psi);
  if(buffer_dpsi != NULL)
    PNX(free)(buffer_dpsi);
}

static void precompute_psi(
    PNX(plan) ths, INT ind, R* x, R* buffer_psi, R* buffer_dpsi,
    int compute_grad_ad,
    R* pre_psi, R* pre_dpsi
    )
{
  int cutoff = ths->cutoff;
  R floor_nx[3];
  for(int t=0; t<3; t++)
    floor_nx[t] = pnfft_floor(ths->n[t]*x[t]);

  if(ths->pnfft_flags & PNFFT_PRE_PSI){
    /* shift index to current particle */
    pre_psi += ind * 3 * cutoff;
    if(compute_grad_ad)
      pre_dpsi += ind * 3 * cutoff;

    pre_psi_tensor(
        ths->n, ths->b, ths->m, ths->cutoff, x, floor_nx,
        ths->exp_const, ths->spline_coeffs, ths->pnfft_flags,
        ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_psi,
        pre_psi);

    if(compute_grad_ad)
      pre_dpsi_tensor(
          ths->n, ths->b, ths->m, ths->cutoff, x, floor_nx, ths->spline_coeffs,
          ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_dpsi,
          pre_psi, ths->pnfft_flags,
          pre_dpsi);
  }
  else if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI){
    /* shift index to current particle */
    pre_psi += ind * PNFFT_POW3(cutoff);
    if(compute_grad_ad)
      pre_dpsi += 3 * ind * PNFFT_POW3(cutoff);

    pre_psi_tensor(
        ths->n, ths->b, ths->m, ths->cutoff, x, floor_nx,
        ths->exp_const, ths->spline_coeffs, ths->pnfft_flags,
        ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_psi,
        buffer_psi);

    if(compute_grad_ad)
        pre_dpsi_tensor(
            ths->n, ths->b, ths->m, ths->cutoff, x, floor_nx, ths->spline_coeffs,
            ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_dpsi,
            buffer_psi, ths->pnfft_flags,
            buffer_dpsi);

    if(compute_grad_ad){
      INT m=0;
      INT md=0;
      for(INT l0=0; l0<cutoff; l0++){
        for(INT l1=cutoff; l1<2*cutoff; l1++){
          R psi_xy  = buffer_psi[l0]   * buffer_psi[l1];
          R psi_dxy = buffer_dpsi[l0]  * buffer_psi[l1];
          R psi_xdy = buffer_psi[l0]   * buffer_dpsi[l1];
          for(INT l2=2*cutoff; l2<3*cutoff; l2++, md+=3){
            pre_psi[m++]   = psi_xy  * buffer_psi[l2];
            pre_dpsi[md+0] = psi_dxy * buffer_psi[l2];
            pre_dpsi[md+1] = psi_xdy * buffer_psi[l2];
            pre_dpsi[md+2] = psi_xy  * buffer_dpsi[l2];
          }
        }
      }
    } else {
      INT m=0;
      for(INT l0=0; l0<cutoff; l0++)
        for(INT l1=cutoff; l1<2*cutoff; l1++)
          for(INT l2=2*cutoff; l2<3*cutoff; l2++)
            pre_psi[m++] = buffer_psi[l0] * buffer_psi[l1] * buffer_psi[l2];
    }
  }
}

static PNX(plan) mkplan(
    void
    )
{
  PNX(plan) ths = (plan_s*) malloc(sizeof(plan_s));

  ths->f_hat  = NULL;
  ths->f      = NULL;
  ths->grad_f = NULL;
  ths->x      = NULL;
  ths->N      = NULL;
  ths->sigma  = NULL;
  ths->n      = NULL;
  ths->no     = NULL;
  ths->x_max  = NULL;

  ths->local_N        = NULL;
  ths->local_N_start  = NULL;
  ths->local_N_total  = 0;
  ths->local_no       = NULL;
  ths->local_no_start = NULL;
  ths->local_no_total = 0;

  ths->b              = NULL;
  ths->exp_const      = NULL;
  ths->spline_coeffs  = NULL;

  ths->pre_inv_phi_hat_trafo  = NULL;
  ths->pre_inv_phi_hat_adj    = NULL;

  ths->pre_psi  = NULL;
  ths->pre_dpsi = NULL;
  ths->pre_psi_il  = NULL;
  ths->pre_dpsi_il = NULL;

  ths->g1 = NULL;
  ths->g2 = NULL;
  ths->g1_buffer = NULL;
  
  ths->pfft_forw = NULL;
  ths->pfft_back = NULL;
  ths->gcplan = NULL;

  ths->intpol_num_nodes = 0;
  ths->intpol_tables_psi  = NULL;
  ths->intpol_tables_dpsi = NULL;

  ths->timer_trafo = PNX(mktimer)();
  ths->timer_adj   = PNX(mktimer)();

  return ths;
}


void PNX(rmplan)(
    PNX(plan) ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;

  PNX(rmtimer)(ths->timer_trafo);
  PNX(rmtimer)(ths->timer_adj);

  /* free memory */
  free(ths);
  /* ths=NULL; would be senseless, since we can not change the pointer itself */
}





void PNX(malloc_x)(
    PNX(plan) ths, unsigned pnfft_flags
    )
{
  if( ~pnfft_flags & PNFFT_MALLOC_X )
    return;

  ths->x = (ths->local_M>0) ? (R*) PNX(malloc)(sizeof(R) * (size_t) ths->d*ths->local_M) : NULL;
  ths->pnfft_flags |= PNFFT_MALLOC_X;
}

void PNX(free_x)(
    PNX(plan) ths, unsigned pnfft_finalize_flags
    )
{
  if( ~pnfft_finalize_flags & PNFFT_FREE_X )
    return;

  if(ths->x != NULL)
    PNX(free)(ths->x);
  ths->pnfft_flags &= ~PNFFT_MALLOC_X;
}

void PNX(malloc_f)(
    PNX(plan) ths, unsigned pnfft_flags
    )
{
  if( ~pnfft_flags & PNFFT_MALLOC_F )
    return;

  ths->f = (ths->local_M>0) ? (R*) PNX(malloc)(sizeof(R) * 2 * (size_t) ths->local_M) : NULL;
  ths->pnfft_flags |= PNFFT_COMPUTE_F;
  ths->compute_flags |= PNFFT_COMPUTE_F;
}

void PNX(free_f)(
    PNX(plan) ths, unsigned pnfft_finalize_flags
    )
{
  if( ~pnfft_finalize_flags & PNFFT_FREE_F )
    return;

  if(ths->f != NULL)
    PNX(free)(ths->f);
  ths->pnfft_flags &= ~PNFFT_MALLOC_F;
  ths->compute_flags &= ~PNFFT_COMPUTE_F;
}

void PNX(malloc_grad_f)(
    PNX(plan) ths, unsigned pnfft_flags
    )
{
  if( ~pnfft_flags & PNFFT_MALLOC_GRAD_F )
    return;

  ths->grad_f = (ths->local_M>0) ? (R*) PNX(malloc)(sizeof(R) * 2 * (size_t) ths->d*ths->local_M) : NULL;
  ths->pnfft_flags |= PNFFT_MALLOC_GRAD_F;
  ths->compute_flags |= PNFFT_COMPUTE_GRAD_F;
}

void PNX(free_grad_f)(
    PNX(plan) ths, unsigned pnfft_finalize_flags
    )
{
  if( ~pnfft_finalize_flags & PNFFT_FREE_GRAD_F )
   return;

  if(ths->grad_f != NULL)
    PNX(free)(ths->grad_f);
  ths->pnfft_flags &= ~PNFFT_MALLOC_GRAD_F;
  ths->compute_flags &= ~PNFFT_COMPUTE_GRAD_F;
}


void PNX(trafo_F)(
    PNX(plan) ths
    )
{
#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->g1, ths->local_N[0]*ths->local_N[1]*ths->local_N[2], 1,
      "PNFFT: Sum of Fourier coefficients before FFT");
#endif

  PX(execute)(ths->pfft_forw);
}

void PNX(adjoint_F)(
    PNX(plan) ths
    )
{
  PX(execute)(ths->pfft_back);

#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->g1, ths->local_N[0]*ths->local_N[1]*ths->local_N[2], 1,
      "PNFFT^H: Sum of Fourier coefficients after FFT");
#endif
}


/* Implement ghostcell send for all dimensions */
static void get_size_gcells(
    int m, int cutoff, unsigned pnfft_flags,
    INT *gcells_below, INT *gcells_above
    )
{
  for(int t=0; t<3; t++){
    gcells_below[t] = m;
    gcells_above[t] = cutoff - gcells_below[t] - 1;
    if(pnfft_flags & PNFFT_INTERLACED)
      gcells_above[t] += 1;
  }
}

static void lowest_summation_index(
    const INT *n, int m, const R *x,
    const INT *local_no_start, const INT *gcells_below,
    R *floor_nx_j, INT *u_j
    )
{
  project_node_to_grid(n, m, x, floor_nx_j, u_j);
  for(int t=0; t<3; t++)
    u_j[t] = u_j[t] - local_no_start[t] + gcells_below[t];
}


static void local_array_size(
    const INT *local_n, const INT *gcells_below, const INT *gcells_above,
    INT *local_ngc
    )
{
  for(int t=0; t<3; t++)
    local_ngc[t] = local_n[t] + gcells_below[t] + gcells_above[t];
}



static void pre_tensor_intpol(
    const INT *n, int cutoff, const R *x, const R *floor_nx,
    int intpol_order, INT intpol_num_nodes, R **intpol_tables_psi,
    R *pre_psi
    )
{
  const int d=3;

  for(int t=0; t<d; t++){
    R dist = n[t]*x[t] - floor_nx[t];
    INT k = (INT) pnfft_floor(dist * intpol_num_nodes) * cutoff * (intpol_order+1);
    R dist_k = 1.0 + pnfft_floor(intpol_num_nodes*dist) - intpol_num_nodes*dist;
    switch(intpol_order){
      case 0 :
        for(int s=0; s<cutoff; s++, k++)
          pre_psi[cutoff*t+s] = pnfft_intpol_const(k, intpol_tables_psi[t]);
        break;
      case 1 :
        for(int s=0; s<cutoff; s++, k+=2)
          pre_psi[cutoff*t+s] = pnfft_intpol_lin(k, dist_k, intpol_tables_psi[t]);
        break;
      case 2 :
        for(int s=0; s<cutoff; s++, k+=3)
          pre_psi[cutoff*t+s] = pnfft_intpol_quad(k, dist_k, intpol_tables_psi[t]);
        break;
      default:
        for(int s=0; s<cutoff; s++, k+=4)
          pre_psi[cutoff*t+s] = pnfft_intpol_kub(k, dist_k, intpol_tables_psi[t]);
    }
  }
}


/* switch between direct evaluation and interpolation */
static void pre_psi_tensor(
    const INT *n, const R *b, int m, int cutoff, const R *x, const R *floor_nx,
    const R *exp_const, R *spline_coeffs, unsigned pnfft_flags,
    int intpol_order, INT intpol_num_nodes, R **intpol_tables_psi,
    R *pre_psi
    )
{
  if(pnfft_flags & PNFFT_PRE_INTPOL_PSI)
    pre_tensor_intpol(
        n, cutoff, x, floor_nx,
        intpol_order, intpol_num_nodes, intpol_tables_psi,
        pre_psi);
  else
    pre_psi_tensor_direct(
        n, b, m, cutoff, x, floor_nx,
        exp_const, spline_coeffs, pnfft_flags,
        pre_psi);
}


/* calculate window function */
static void pre_psi_tensor_direct(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx,
    const R *exp_const, R *spline_coeffs,
    unsigned pnfft_flags,
    R *pre_psi
    )
{
  if(pnfft_flags & PNFFT_WINDOW_GAUSSIAN){
    if(pnfft_flags & PNFFT_FG_PSI)
      pre_psi_tensor_fast_gaussian(
          n, b, m, cutoff, x, floor_nx, exp_const,
          pre_psi);
    else
      pre_psi_tensor_gaussian(
          n, b, m, cutoff, x, floor_nx,
          pre_psi);
  } 
  else if(pnfft_flags & PNFFT_WINDOW_BSPLINE)
    pre_psi_tensor_bspline(
        n, m, cutoff, x, floor_nx, spline_coeffs,
        pre_psi);
  else if(pnfft_flags & PNFFT_WINDOW_SINC_POWER)
    pre_psi_tensor_sinc_power(
        n, b, m, cutoff, x, floor_nx,
        pre_psi);
  else if(pnfft_flags & PNFFT_WINDOW_BESSEL_I0)
    pre_psi_tensor_bessel_i0(
        n, b, m, cutoff, x, floor_nx,
        pre_psi);
  else
    pre_psi_tensor_kaiser_bessel(
        n, b, m, cutoff, x, floor_nx,
        pre_psi);
}

static void pre_psi_tensor_gaussian(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx,
    R *pre_psi
    )
{
  const int d=3;
  R u_j;

  for(int t=0; t<d; t++){
    u_j = floor_nx[t] - n[t]*x[t] - m;
    for(int s=0; s<cutoff; s++)
      pre_psi[cutoff*t+s] = pnfft_exp(-PNFFT_SQR(u_j + s) / b[t]) / pnfft_sqrt(PNFFT_PI*b[t]);
  }
}

static void pre_psi_tensor_fast_gaussian(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, const R *exp_const,
    R *fg_psi
    )
{
  const int d=3;
  R u_j, exp_sqr, exp_lin, tmp;

  for(int t=0; t<d; t++){
    u_j = floor_nx[t] - m;
    exp_sqr = pnfft_exp( -PNFFT_SQR( n[t]*x[t]-u_j ) / b[t] );
    exp_lin = pnfft_exp( 2*( n[t]*x[t]-u_j ) / b[t] );

    tmp = exp_sqr;
    for(int s=0; s<cutoff; s++){
      fg_psi[cutoff*t+s] = tmp * exp_const[cutoff*t+s];
      tmp *= exp_lin;
    }
  }
}

static void pre_psi_tensor_bspline(
    const INT *n, int m, int cutoff,
    const R *x, const R *floor_nx, R *spline_coeffs, 
    R *pre_psi
    )
{
  const int d=3;

  if(m<9){
    for(int t=0; t<d; t++){
      /* Bspline is shifted by m */
      R dist = floor_nx[t]  - n[t]*x[t] + 0.5; 
      for(int s=0; s<cutoff; s++)
        pre_psi[cutoff*t+s] = PNX(fast_bspline)(
            s-1, dist, 2*m);
    }
    return;
  }

  for(int t=0; t<d; t++){
    /* Bspline is shifted by m */
    R u_j = floor_nx[t] - n[t]*x[t]; 
    for(int s=0; s<cutoff; s++)
      pre_psi[cutoff*t+s] = PNX(bspline)(
          2*m, u_j + (R)s, spline_coeffs);
  }
}

/* The factor n of the window cancels with the factor 1/n from matrix D. */
static void pre_psi_tensor_sinc_power(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx,
    R *pre_psi
    )
{
  const int d=3;
  R u_j;

  for(int t=0; t<d; t++){
    u_j = floor_nx[t] - n[t]*x[t] - m; 
    for(int s=0; s<cutoff; s++)
      pre_psi[cutoff*t+s] =
        pnfft_pow(
            PNX(sinc)( PNFFT_PI * (u_j + s) / b[t]),
            K(2.0)*(R)m
        ) / b[t];
  }
}

static void pre_psi_tensor_bessel_i0(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx,
    R *pre_psi
    )
{
  const int d=3;
  R u_j;

  for(int t=0; t<d; t++){
    u_j = floor_nx[t] - n[t]*x[t] - m;
    for(int s=0; s<cutoff; s++)
      pre_psi[cutoff*t+s] = window_bessel_i0_1d(
          (u_j + s) / n[t], n[t], b[t], m);
  }
}

static void pre_psi_tensor_kaiser_bessel(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx,
    R *pre_psi
    )
{
  const int d=3;
  R u_j;

  for(int t=0; t<d; t++){
    u_j = floor_nx[t] - n[t]*x[t] - m;
    for(int s=0; s<cutoff; s++)
      pre_psi[cutoff*t+s] = kaiser_bessel_1d(
          (u_j + s) / n[t], n[t], b[t], m);
  }
}

/* switch between direct evaluation and interpolation */
static void pre_dpsi_tensor(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, R *spline_coeffs,
    int intpol_order, INT intpol_num_nodes, R **intpol_tables_dpsi,
    const R *pre_psi, unsigned pnfft_flags,
    R *pre_dpsi
    )
{
  if(pnfft_flags & PNFFT_PRE_INTPOL_PSI)
    pre_tensor_intpol(
        n, cutoff, x, floor_nx,
        intpol_order, intpol_num_nodes, intpol_tables_dpsi,
        pre_dpsi);
  else
    pre_dpsi_tensor_direct(
        n, b, m, cutoff, x, floor_nx, spline_coeffs,
        pre_psi, pnfft_flags,
        pre_dpsi);
}

/* calculate window derivative */
static void pre_dpsi_tensor_direct(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, R *spline_coeffs, const R *pre_psi,
    unsigned pnfft_flags,
    R *pre_dpsi
    )
{
  if(pnfft_flags & PNFFT_WINDOW_GAUSSIAN)
    pre_dpsi_tensor_gaussian(
        n, b, m, cutoff, x, floor_nx, pre_psi,
        pre_dpsi);
  else if(pnfft_flags & PNFFT_WINDOW_BSPLINE)
    pre_dpsi_tensor_bspline(
        n, m, cutoff, x, floor_nx, spline_coeffs,
        pre_dpsi);
  else if(pnfft_flags & PNFFT_WINDOW_SINC_POWER)
    pre_dpsi_tensor_sinc_power(
        n, b, m, cutoff, x, floor_nx, pre_psi,
        pre_dpsi);
  else if(pnfft_flags & PNFFT_WINDOW_BESSEL_I0)
    pre_dpsi_tensor_bessel_i0(
        n, b, m, cutoff, x, floor_nx,
        pre_dpsi);
  else
    pre_dpsi_tensor_kaiser_bessel(
        n, b, m, cutoff, x, floor_nx, pre_psi,
        pre_dpsi);
}

static void pre_dpsi_tensor_gaussian(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, const R *fg_psi,
    R *fg_dpsi
    )
{
  const int d=3;
  R u_j;

  for(int t=0; t<d; t++){
    u_j = n[t]*x[t] - floor_nx[t] + m;
    for(int s=0; s<cutoff; s++)
      fg_dpsi[cutoff*t+s] = -2.0*n[t]/b[t] * (u_j-s) * fg_psi[cutoff*t+s];
  }
}

static void pre_dpsi_tensor_bspline(
    const INT *n, int m, int cutoff,
    const R *x, const R *floor_nx, R *spline_coeffs,
    R *pre_dpsi
    )
{
  const int d=3;

  if(m<9){
    for(int t=0; t<d; t++){
      /* Bspline is shifted by m */
      R dist = floor_nx[t] - n[t]*x[t] + 0.5;
      for(int s=0; s<cutoff; s++)
        pre_dpsi[cutoff*t+s] = -(R)n[t] * PNX(fast_bspline_d)(
            s-1, dist, 2*m);
    }
    return;
  }

  for(int t=0; t<d; t++){
    /* Bspline is shifted by m */
    R u_j = floor_nx[t] - n[t]*x[t];
    for(int s=0; s<cutoff; s++)
      pre_dpsi[cutoff*t+s] = (R)n[t] *
        ( PNX(bspline)(2*m-1, u_j + (R)s - 1, spline_coeffs)
          - PNX(bspline)(2*m-1, u_j + (R)s, spline_coeffs) );
  }
}


/* The factor n of the window cancels with the factor 1/n from matrix D. */
static void pre_dpsi_tensor_sinc_power(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, const R *pre_psi,
    R *pre_dpsi
    )
{
  const int d=3;
  R u_j, y;

  for(int t=0; t<d; t++){
    u_j =  floor_nx[t] - n[t]*x[t] - m;
    for(int s=0; s<cutoff; s++){
      y =  PNFFT_PI * (u_j + s) / b[t];
      if(pnfft_fabs(y) > PNFFT_EPSILON)
        pre_dpsi[cutoff*t+s] =
          2.0 * (R)m * PNFFT_PI * (R)n[t] / b[t] * ( 1.0/y - 1.0/pnfft_tan(y) ) * pre_psi[cutoff*t+s];
      else
        pre_dpsi[cutoff*t+s] = K(0.0);
    }
  }
}

static void pre_dpsi_tensor_bessel_i0(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx,
    R *pre_dpsi
    )
{
  const int d=3;
  R u_j;

  for(int t=0; t<d; t++){
    u_j = floor_nx[t] - n[t]*x[t] - m;
    for(int s=0; s<cutoff; s++)
      pre_dpsi[cutoff*t+s] = window_bessel_i0_derivative_1d(
          (u_j + s) / n[t], n[t], b[t], m);
  }
}


static void pre_dpsi_tensor_kaiser_bessel(
    const INT *n, const R *b, int m, int cutoff,
    const R *x, const R *floor_nx, const R *pre_psi,
    R *pre_dpsi
    )
{
  const int d=3;
  R u_j;

  for(int t=0; t<d; t++){
    u_j =  floor_nx[t] - n[t]*x[t] - m;
    for(int s=0; s<cutoff; s++)
      pre_dpsi[cutoff*t+s] = kaiser_bessel_derivative_1d(
          (u_j + s) / n[t],
          n[t], b[t], m, pre_psi[cutoff*t+s]);
  }
}






/**
 * Sort nodes (index) to get better cache utilization during multiplication
 * with matrix B.
 * The resulting index set is written to ar[2*j+1], the nodes array remains
 * unchanged.
 *
 * \arg n FFTW length (number of oversampled grid points in each dimension)
 * \arg m window length
 * \arg local_x_num number of nodes
 * \arg local_x nodes array
 * \arg ar_x resulting index array
 *
 * \author Toni Volkmer
 */
static void sort_nodes_for_better_cache_handle(
    int d, const INT *n, int m, INT local_x_num, const R *local_x,
    INT *ar_x
    )
{
#if PNFFT_SORT_RADIX
  INT u_j[d], i, j, help, rhigh;
  INT *ar_x_temp;
  R nprod;

  for(i = 0; i < local_x_num; i++) {
    ar_x[2*i] = 0;
    ar_x[2*i+1] = i;
    for(j = 0; j < d; j++) {
      help = pnfft_floor( n[j]*local_x[d*i+j] - m);
      u_j[j] = (help%n[j]+n[j])%n[j];

      ar_x[2*i] += u_j[j];
      if (j+1 < d)
        ar_x[2*i] *= n[j+1];
    }
  }

  for (j = 0, nprod = 1.0; j < d; j++)
    nprod *= n[j];

  rhigh = pnfft_ceil(pnfft_log2(nprod)) - 1;

  ar_x_temp = (INT*) PNX(malloc)(2*local_x_num*sizeof(INT));
  PNX(sort_node_indices_radix_lsdf)(local_x_num, ar_x, ar_x_temp, rhigh);
#  ifdef OMP_ASSERT
  for (i = 1; i < local_x_num; i++)
    assert(ar_x[2*(i-1)] <= ar_x[2*i]);
#  endif
  PNX(free)(ar_x_temp);
#else
  PNX(sort_nodes_indices_qsort_3d)(n, m, local_x_num, local_x, ar_x);
#endif
}





static void project_node_to_grid(
    const INT *n, int m, const R *x,
    R *floor_nx_j, INT *u_j
    )
{
  for(int t=0; t<3; t++){
    floor_nx_j[t] = pnfft_floor(n[t]*x[t]);
    u_j[t] = (INT) floor_nx_j[t] - m;
  }
}

static void get_mpi_cart_dims_3d(
    MPI_Comm comm_cart,
    int *rnk_pm, int *dims, int *coords
    )
{
  int *periods;

  /* takes care of the case rnk_pm < 3 */  
  for(int t=0; t<3; t++) dims[t] = 1;
  for(int t=0; t<3; t++) coords[t] = 0;

  MPI_Cartdim_get(comm_cart, rnk_pm);

  periods = (int*) malloc(sizeof(int) * (size_t) *rnk_pm);
  coords  = (int*) malloc(sizeof(int) * (size_t) *rnk_pm);
  MPI_Cart_get(comm_cart, *rnk_pm, dims, periods, coords);
  free(periods);
}

/* compact support in real space */
static R window_bessel_i0_1d(
    R x, INT n, R b, int m
    )
{
  R d = PNFFT_SQR( (R)m ) - PNFFT_SQR( x*n );

  /* Compact support in real space */
  return (d<0) ? 0.0 : 0.5 * PNX(bessel_i0)(b*pnfft_sqrt(d));
}

static R window_bessel_i0_derivative_1d(
    R x, INT n, R b, int m
    )
{
  R d = PNFFT_SQR( (R)m ) - PNFFT_SQR( x*n );
  R r = (d<0) ? pnfft_sqrt(-d) : pnfft_sqrt(d);

  /* Compact support in real space */
  if(d<0)
   return 0.0;

  /* avoid division by zero */
  return (d>0) ? 0.5 * b * (R)n * (R)n * x * PNX(bessel_i1)(b*r) / r : PNFFT_SQR(0.5*b*n) * x;
}

static R kaiser_bessel_1d(
    R x, INT n, R b, int m
    )
{
  /* TODO: try to avoid case d<0, since in theory d >= 0 */
  R d = PNFFT_SQR( (R)m ) - PNFFT_SQR( x*n );
  R r = (d<0) ? pnfft_sqrt(-d) : pnfft_sqrt(d);

  if(d < 0)
    return pnfft_sin(b*r) / ((R)PNFFT_PI*r);

  /* regular case d >= 0 */
  return (d>0) ? pnfft_sinh(b*r) / ((R)PNFFT_PI*r) : b/PNFFT_PI;
}

static R kaiser_bessel_derivative_1d(
    R x, INT n, R b, int m, R psi
    )
{
  /* TODO: try to avoid case d<0, since in theory d >= 0 */
  R d = PNFFT_SQR( (R)m ) - PNFFT_SQR( x*n );
  R r = (d<0) ? pnfft_sqrt(-d) : pnfft_sqrt(d);

  if(d < 0) /* use of -d results in change of sign within the brackets */
    return n*n*x/d * (b*pnfft_cos(b*r)/PNFFT_PI - psi);

  /* regular case d >= 0 */
  return (d>0) ? n*n*x/d * (b*pnfft_cosh(b*r)/PNFFT_PI - psi) : -n*m*b*b*b/(3.0*PNFFT_PI);
}




R PNX(psi)(
    const PNX(plan) ths, int dim, R x
    )
{
  if(ths->pnfft_flags & PNFFT_WINDOW_GAUSSIAN)
    return psi_gaussian(x, ths->n[dim], ths->b[dim]);
  else if(ths->pnfft_flags & PNFFT_WINDOW_BSPLINE)
    return psi_bspline(x, ths->n[dim], ths->m, ths->spline_coeffs);
  else if(ths->pnfft_flags & PNFFT_WINDOW_SINC_POWER)
    return psi_sinc_power(x, ths->n[dim], ths->b[dim], ths->m);
  else if(ths->pnfft_flags & PNFFT_WINDOW_BESSEL_I0)
    return psi_bessel_i0(x, ths->n[dim], ths->b[dim], ths->m);
  else
    return psi_kaiser(x, ths->n[dim], ths->b[dim], ths->m);
}

R PNX(dpsi)(
    const PNX(plan) ths, int dim, R x
    )
{
  if(ths->pnfft_flags & PNFFT_WINDOW_GAUSSIAN)
    return dpsi_gaussian(x, ths->n[dim], ths->b[dim]);
  else if(ths->pnfft_flags & PNFFT_WINDOW_BSPLINE)
    return dpsi_bspline(x, ths->n[dim], ths->m, ths->spline_coeffs);
  else if(ths->pnfft_flags & PNFFT_WINDOW_SINC_POWER)
    return dpsi_sinc_power(x, ths->n[dim], ths->b[dim], ths->m);
  else if(ths->pnfft_flags & PNFFT_WINDOW_BESSEL_I0)
    return dpsi_bessel_i0(x, ths->n[dim], ths->b[dim], ths->m);
  else
    return dpsi_kaiser(x, ths->n[dim], ths->b[dim], ths->m);
}



static R psi_gaussian(
    R x, INT n, R b
    )
{
  return pnfft_exp( -PNFFT_SQR(n*x)/b ) / pnfft_sqrt(PNFFT_PI*b);
}

static R dpsi_gaussian(
    R x, INT n, R b
    )
{
  return 2.0 * n/b * n*x * pnfft_exp( -PNFFT_SQR(n*x)/b ) / pnfft_sqrt(PNFFT_PI*b);
}

static R psi_bspline(
    R x, INT n, int m, R *spline_coeffs
    )
{
  /* Bspline is shifted by m */
  return PNX(bspline)(2*m, n*x+m, spline_coeffs);
}

static R dpsi_bspline(
    R x, INT n, int m, R *spline_coeffs
    )
{
  /* Bspline is shifted by m */
  return (R)n *( PNX(bspline)(2*m-1, n*x + m - 1, spline_coeffs)
      - PNX(bspline)(2*m-1, n*x + m, spline_coeffs) );
}

static R psi_sinc_power(
    R x, INT n, R b, int m
    )
{
  return pnfft_pow
    (
     PNX(sinc)( PNFFT_PI * (n*x) / b),
     K(2.0)*(R)m
    ) / b;
}

static R dpsi_sinc_power(
    R x, INT n, R b, int m
    )
{
  R y =  PNFFT_PI * n * x / b;
  if(pnfft_fabs(y) > PNFFT_EPSILON)
    return 2.0 * (R)m * PNFFT_PI * (R)n / b * ( 1.0/y - 1.0/pnfft_tan(y) ) * psi_sinc_power(x, n, b, m);
  else
    return 0.0;
}

static R psi_bessel_i0(
    R x, INT n, R b, int m
    )
{
  return window_bessel_i0_1d(x, n, b, m);
}

static R dpsi_bessel_i0(
    R x, INT n, R b, int m
    )
{
  return window_bessel_i0_derivative_1d(x, n, b, m);
}

static R psi_kaiser(
    R x, INT n, R b, int m
    )
{
  return kaiser_bessel_1d(x, n, b, m);
}

static R dpsi_kaiser(
    R x, INT n, R b, int m
    )
{
  return kaiser_bessel_derivative_1d(
      x, n, b, m, kaiser_bessel_1d(x, n, b, m));
}


/* Alternative sorting based on C standard qsort */
void PNX(sort_nodes_indices_qsort_3d)(
    const INT *n, int m, INT local_M, const R *x,
    INT *sort
    )
{ /* sort must be of length 2*local_M */
  int d = 3;
  INT u_j[3];
  R floor_nx_j[3];

  if(local_M == 0)
    return;

  for(INT k = 0; k < local_M; k++){
    project_node_to_grid(n, m, x+d*k, floor_nx_j, u_j);
    sort[2*k]   = PNFFT_PLAIN_INDEX_3D(u_j, n);
    sort[2*k+1] = k;
  }
  qsort(sort, (size_t) local_M, 2*sizeof(INT), compare_INT);
}

static int compare_INT(
    const void* a, const void* b
    )
{
  return (*(const INT*)a <= *(const INT*)b) ? -1 : 1;
}


// static void ik_diff_real_input(
//     PNX(plan) ths
//     )
// {
//   /* duplicate g1 since we have to scale it several times for computing the gradient */
//   for(INT k=0; k<ths->local_N_total; k++)
//     ths->g1_buffer[k] = ths->g1[k];
// 
//   /* calculate potential and 1st component of gradient */
//   scale_ik_diff_r2ci(g1_buffer, ths->local_N_start, ths->local_N, dim=0,
//     ths->g1);
//   PNX(trafo_F)(ths);
//   PNX(new_trafo_B)(ths, ths->grad_f, dim);
//   
//   /* calculate 2nd and 3rd component of gradient */
//   scale_ik_diff_r2cr(g1_buffer, ths->local_N_start, ths->local_N, dim=1,
//     ths->g1);
//   scale_ik_diff_r2ci(g1_buffer, ths->local_N_start, ths->local_N, dim=2,
//     ths->g1);
//   PNX(trafo_F)(ths);
//   PNX(new_trafo_B)(ths, ths->grad_f, dim);
// }


void PNX(trafo_B_grad_ik)(
    PNX(plan) ths, R *f, INT offset, INT stride
    )
{
  const int cutoff = ths->cutoff;
  INT j, m0;
  INT *sorted_index = NULL;
  INT local_no[3], local_no_start[3];
  INT gcells_below[3], gcells_above[3];
  INT u_j[3], local_ngc[3];
  R floor_nx_j[3];
  R *pre_psi = NULL;
 
  local_size_B(ths,
      local_no, local_no_start);

  /* perform fftshift */
  ths->timer_trafo[PNFFT_TIMER_SHIFT_INPUT] -= MPI_Wtime();
  ths->timer_trafo[PNFFT_TIMER_SHIFT_INPUT] += MPI_Wtime();

  get_size_gcells(ths->m, ths->cutoff, ths->pnfft_flags,
      gcells_below, gcells_above);
  local_array_size(local_no, gcells_below, gcells_above,
      local_ngc);

  /* send ghost cells in ring */
  ths->timer_trafo[PNFFT_TIMER_GCELLS] -= MPI_Wtime();
  PX(exchange)(ths->gcplan);
  ths->timer_trafo[PNFFT_TIMER_GCELLS] += MPI_Wtime();

  /* sort indices for better cache handling */
  if(ths->pnfft_flags & PNFFT_SORT_NODES){
    ths->timer_trafo[PNFFT_TIMER_SORT_NODES] -= MPI_Wtime();
    sorted_index = (INT*) PNX(malloc)(sizeof(INT) * (size_t) 2*ths->local_M);
    sort_nodes_for_better_cache_handle(
        ths->d, ths->n, ths->m, ths->local_M, ths->x,
        sorted_index);
    ths->timer_trafo[PNFFT_TIMER_SORT_NODES] += MPI_Wtime();
  }

  ths->timer_trafo[PNFFT_TIMER_LOOP_B] -= MPI_Wtime();
  if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) )
    pre_psi = (R*) PNX(malloc)(sizeof(R) * (size_t) cutoff*3);
  for(INT p=0; p<ths->local_M; p++){
    j = (ths->pnfft_flags & PNFFT_SORT_NODES) ? sorted_index[2*p+1] : p;

    lowest_summation_index(
        ths->n, ths->m, &(ths->x[ths->d*j]), local_no_start, gcells_below,
        floor_nx_j, u_j);
    
    /* evaluate window on axes */
    if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) )
      pre_psi_tensor(
          ths->n, ths->b, ths->m, ths->cutoff, &(ths->x[ths->d*j]), floor_nx_j,
          ths->exp_const, ths->spline_coeffs, ths->pnfft_flags,
          ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_psi,
          pre_psi);

    /* compute f */
    INT ind = j*stride + offset;
    m0 = PNFFT_PLAIN_INDEX_3D(u_j, local_ngc);
    if(ths->pnfft_flags & PNFFT_REAL_F)
      PNX(assign_f_r2r)(
          ths, p, ths->g2, pre_psi, 2*m0, local_ngc, cutoff, 2, 0,
          f + 2*ind);
    else if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
      PNX(assign_f_r2r)(
          ths, p, ths->g2, pre_psi, m0, local_ngc, cutoff, 1, 0,
          f + ind);
    else
      PNX(assign_f_c2c)(
          ths, p, ths->g2, pre_psi, m0, local_ngc, cutoff, 0,
          (C*)f + ind);
  }

  ths->timer_trafo[PNFFT_TIMER_LOOP_B] += MPI_Wtime();
  
  if(sorted_index != NULL)
    PNX(free)(sorted_index);
  if(pre_psi != NULL)
    PNX(free)(pre_psi);
}









void PNX(trafo_B_grad_ad)(
    PNX(plan) ths, int interlaced
    )
{
  const int cutoff = ths->cutoff;
  INT j, m0;
  INT *sorted_index = NULL;
  INT local_no[3], local_no_start[3];
  INT gcells_below[3], gcells_above[3];
  INT u_j[3], local_ngc[3];
  R floor_nx_j[3];
  R *pre_psi = NULL;
  R *pre_dpsi = NULL;
 
  local_size_B(ths,
      local_no, local_no_start);

#if PNFFT_ENABLE_DEBUG
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  if(!myrank) fprintf(stderr, "no = [%td, %td, %td]\n", ths->no[0], ths->no[1], ths->no[2]);
//  if(!myrank) fprintf(stderr, "local_no = [%td, %td, %td], local_ngc = [%td, %td, %td], local_ngc_total = %td\n",
//      local_no[0], local_no[1], local_no[2], local_ngc[0], local_ngc[1], local_ngc[2], local_ngc_total);
      
  PNX(debug_sum_print)(ths->g2, local_no[0]*local_no[1]*local_no[2],
      !(ths->trafo_flag & PNFFTI_TRAFO_C2R),
      "PNFFT: Sum of Fourier coefficients before twiddles");
#endif

  /* perform fftshift */
  ths->timer_trafo[PNFFT_TIMER_SHIFT_INPUT] -= MPI_Wtime();
  ths->timer_trafo[PNFFT_TIMER_SHIFT_INPUT] += MPI_Wtime();

  get_size_gcells(ths->m, ths->cutoff, ths->pnfft_flags,
      gcells_below, gcells_above);
  local_array_size(local_no, gcells_below, gcells_above,
      local_ngc);

#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->g2, local_no[0]*local_no[1]*local_no[2],
      !(ths->trafo_flag & PNFFTI_TRAFO_C2R),
      "PNFFT: Sum of Fourier coefficients before ghostcell send");
#endif

  /* send ghost cells in ring */
  ths->timer_trafo[PNFFT_TIMER_GCELLS] -= MPI_Wtime();
  PX(exchange)(ths->gcplan);
  ths->timer_trafo[PNFFT_TIMER_GCELLS] += MPI_Wtime();

#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->g2, PNX(prod_INT)(3, local_ngc),
      !(ths->trafo_flag & PNFFTI_TRAFO_C2R),
      "PNFFT: Sum of Fourier coefficients after ghostcell send");
#endif  

#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->x, 3*ths->local_M, 0,
      "PNFFT: Sum of x before sort");
#endif

  /* sort indices for better cache handling */
  if(ths->pnfft_flags & PNFFT_SORT_NODES){
    ths->timer_trafo[PNFFT_TIMER_SORT_NODES] -= MPI_Wtime();
    sorted_index = (INT*) PNX(malloc)(sizeof(INT) * (size_t) 2*ths->local_M);
    sort_nodes_for_better_cache_handle(
        ths->d, ths->n, ths->m, ths->local_M, ths->x,
        sorted_index);
    ths->timer_trafo[PNFFT_TIMER_SORT_NODES] += MPI_Wtime();
  }

#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->x, 3*ths->local_M, 0,
      "PNFFT: Sum of x after sort");
#endif

#if PNFFT_ENABLE_DEBUG
  R rsum, rsum_derive, grsum, grsum_derive;
  rsum = rsum_derive = 0.0;
#endif

  ths->timer_trafo[PNFFT_TIMER_LOOP_B] -= MPI_Wtime();
  if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) ){
    pre_psi = (R*) PNX(malloc)(sizeof(R) * (size_t) cutoff*3);
    if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F)
      pre_dpsi = (R*) PNX(malloc)(sizeof(R) * (size_t) cutoff*3);
  }
  for(INT p=0; p<ths->local_M; p++){
    j = (ths->pnfft_flags & PNFFT_SORT_NODES) ? sorted_index[2*p+1] : p;
    
    lowest_summation_index(
        ths->n, ths->m, &(ths->x[ths->d*j]), local_no_start, gcells_below,
        floor_nx_j, u_j);
    
    if(ths->compute_flags & PNFFT_COMPUTE_F) {
      if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
        ths->f[j] = 0;
      else
        ((C*)ths->f)[j] = 0;
    }
    if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F)
      if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
        for(int t=0; t<ths->d; t++)
          ths->grad_f[ths->d*j+t] = 0;
      else
        for(int t=0; t<ths->d; t++)
          ((C*)ths->grad_f)[ths->d*j+t] = 0;

    /* evaluate window on axes */
    if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) ){
      pre_psi_tensor(
          ths->n, ths->b, ths->m, ths->cutoff, &(ths->x[ths->d*j]), floor_nx_j,
          ths->exp_const, ths->spline_coeffs, ths->pnfft_flags,
          ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_psi,
          pre_psi);
  
#if PNFFT_ENABLE_DEBUG
      /* Don't want to use PNX(debug_sum_print) because we are in a loop */
      for(int t=0; t<3*cutoff; t++)
        rsum += pnfft_fabs(pre_psi[t]);
#endif
 
      if(ths->compute_flags & (PNFFT_COMPUTE_GRAD_F)){
        pre_dpsi_tensor(
            ths->n, ths->b, ths->m, ths->cutoff, &(ths->x[ths->d*j]), floor_nx_j, ths->spline_coeffs,
            ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_dpsi,
            pre_psi, ths->pnfft_flags,
            pre_dpsi);

#if PNFFT_ENABLE_DEBUG
        /* Don't want to use PNX(debug_sum_print) because we are in a loop */
        if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F)
          for(int t=0; t<3*cutoff; t++)
            rsum_derive += pnfft_fabs(pre_dpsi[t]);
#endif
      }
    }

    m0 = PNFFT_PLAIN_INDEX_3D(u_j, local_ngc);
    if(interlaced){
      if(ths->compute_flags & PNFFT_COMPUTE_F 
         && ths->compute_flags & PNFFT_COMPUTE_GRAD_F){
        /* compute f and grad_f at once */
        PNX(assign_f_and_grad_f_r2r)(
            ths, p, ths->g2, pre_psi, pre_dpsi, m0, local_ngc, cutoff, 1, 1, 0,
            ths->f + j, ths->grad_f + 3*j);
      }
      else if(ths->compute_flags & PNFFT_COMPUTE_F){
        /* compute f */
        PNX(assign_f_r2r)(
            ths, p, ths->g2, pre_psi, m0, local_ngc, cutoff, 1, 0,
            ths->f + j);
      } 
      else if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F){
        /* compute grad_f */
        PNX(assign_grad_f_r2r)(
            ths, p, ths->g2, pre_psi, pre_dpsi, m0, local_ngc, cutoff, 1, 1, 0,
            ths->grad_f + 3*j);
      }
    } else {
      if(ths->compute_flags & PNFFT_COMPUTE_F 
         && ths->compute_flags & PNFFT_COMPUTE_GRAD_F){
        /* compute f and grad_f at once */
        if(ths->pnfft_flags & PNFFT_REAL_F)
          PNX(assign_f_and_grad_f_r2r)(
              ths, p, ths->g2, pre_psi, pre_dpsi,
              2*m0, local_ngc, cutoff, 2, 2, 0,
              ths->f + 2*j, ths->grad_f + 2*3*j);
        else if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
          PNX(assign_f_and_grad_f_r2r)(
              ths, p, ths->g2, pre_psi, pre_dpsi,
              m0, local_ngc, cutoff, 1, 1, 0,
              ths->f + j, ths->grad_f + 3*j);
        else
          PNX(assign_f_and_grad_f_c2c)(
              ths, p, ths->g2, pre_psi, pre_dpsi,
              m0, local_ngc, cutoff, 0,
              (C*)ths->f + j, (C*)ths->grad_f + 3*j);
      } else if(ths->compute_flags & PNFFT_COMPUTE_F){
        /* compute f */
        if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
          PNX(assign_f_r2r)(
            ths, p, ths->g2, pre_psi, m0, local_ngc, cutoff, 1, 0,
            ths->f + j);
        else
          PNX(assign_f_c2c)(
              ths, p, ths->g2, pre_psi, m0, local_ngc, cutoff, 0,
              (C*)ths->f + j);
      } else if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F){
        /* compute grad_f */
        PNX(assign_grad_f_c2c)(
            ths, p, ths->g2, ths->pre_psi + p*PNFFT_POW3(cutoff), ths->pre_dpsi + 3*p*PNFFT_POW3(cutoff),
            m0, local_ngc, cutoff, 0,
            (C*)ths->grad_f + 3*j);
      }
    }

    if(interlaced){
      R x[3];
      
      /* shift x by half the mesh width */
      for(int t=0; t<3; t++)
        x[t] = ths->x[ths->d*j+t] + 0.5/ths->n[t];

      /* We need to compute the lowest summation index before we fold x back into [-0.5,0.5).
       * Otherwise u_j may be also folded and gets less than the local offset local_no_start. */
      lowest_summation_index(
          ths->n, ths->m, x, local_no_start, gcells_below,
          floor_nx_j, u_j);
      
      /* assure -0.5 <= x < 0.5 */
      for(int t=0; t<3; t++){
        if(x[t] >= 0.5){
          x[t] -= 1.0;
          floor_nx_j[t] -= ths->n[t];
        }
      }

      /* evaluate window on axes */
      if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) ){
        pre_psi_tensor(
            ths->n, ths->b, ths->m, ths->cutoff, x, floor_nx_j,
            ths->exp_const, ths->spline_coeffs, ths->pnfft_flags,
            ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_psi,
            pre_psi);

        if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F)
          pre_dpsi_tensor(
              ths->n, ths->b, ths->m, ths->cutoff, x, floor_nx_j, ths->spline_coeffs,
              ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_dpsi,
              pre_psi, ths->pnfft_flags,
              pre_dpsi);
      }
      
      m0 = PNFFT_PLAIN_INDEX_3D(u_j, local_ngc);
      if(ths->compute_flags & PNFFT_COMPUTE_F 
         && ths->compute_flags & PNFFT_COMPUTE_GRAD_F){
        /* compute f and grad_f at once */
        PNX(assign_f_and_grad_f_r2r)(
            ths, p, ths->g2, pre_psi, pre_dpsi, 2*m0+1, local_ngc, cutoff, 2, 2, 1, 
            ths->f + j*2, ths->grad_f + 3*j*2);
      }
      else if(ths->compute_flags & PNFFT_COMPUTE_F){
        /* compute f */
        PNX(assign_f_r2r)(
            ths, p, ths->g2, pre_psi, 2*m0+1, local_ngc, cutoff, 2, 1,
            ths->f + j*2);
      } 
      else if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F){
        /* compute grad_f */
        PNX(assign_grad_f_r2r)(
            ths, p, ths->g2, pre_psi, pre_dpsi, 2*m0+1, local_ngc, cutoff, 2, 2, 1,
            ths->grad_f + 3*j*2);
      }

      if(ths->compute_flags & PNFFT_COMPUTE_F) {
        if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
          ths->f[j] /= 2.0;
        else
          ((C*)ths->f)[j] /= 2.0;
      }
      if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F) {
        if(ths->trafo_flag & PNFFTI_TRAFO_C2R)
          for(int t=0; t<ths->d; t++)
            ths->grad_f[ths->d*j+t] /= 2.0;
        else
          for(int t=0; t<ths->d; t++)
            ((C*)ths->grad_f)[ths->d*j+t] /= 2.0;
      }
    }
  }

  ths->timer_trafo[PNFFT_TIMER_LOOP_B] += MPI_Wtime();

#if PNFFT_ENABLE_DEBUG
  MPI_Reduce(&rsum, &grsum, 1, PNFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!myrank) fprintf(stderr, "PNFFT: Sum of pre_psi: %e\n", grsum);

  if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F){
    MPI_Reduce(&rsum_derive, &grsum_derive, 1, PNFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(!myrank) fprintf(stderr, "PNFFT: Sum of pre_dpsi: %e\n", grsum_derive);
  }
#endif

#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->f, ths->local_M,
      !(ths->trafo_flag & PNFFTI_TRAFO_C2R),
      "PNFFT: Sum of f");

  if(ths->compute_flags & PNFFT_COMPUTE_GRAD_F){
    PNX(debug_sum_print_strides)(ths->grad_f, ths->local_M, 3,
        !(ths->trafo_flag & PNFFTI_TRAFO_C2R),
        "PNFFT: Sum of %dst component of grad_f");
  }
#endif
  
  if(sorted_index != NULL)
    PNX(free)(sorted_index);
  if(pre_psi != NULL)
    PNX(free)(pre_psi);
  if(pre_dpsi != NULL)
    PNX(free)(pre_dpsi);
}




void PNX(adjoint_B)(
    PNX(plan) ths, int interlaced
    )
{
  INT *sorted_index = NULL;
  INT local_no[3], local_no_start[3];
  INT gcells_below[3], gcells_above[3];
  INT local_ngc[3], local_ngc_total;

  local_size_B(ths,
      local_no, local_no_start);

  get_size_gcells(ths->m, ths->cutoff, ths->pnfft_flags,
      gcells_below, gcells_above);
  local_array_size(local_no, gcells_below, gcells_above,
      local_ngc);

  local_ngc_total = PNX(prod_INT)(3, local_ngc);
  if (ths->trafo_flag & PNFFTI_TRAFO_C2R)
    for(INT k=0; k<local_ngc_total; k++)
      ths->g2[k] = 0;
  else
    for(INT k=0; k<local_ngc_total; k++)
      ((C*)ths->g2)[k] = 0;

#if PNFFT_ENABLE_DEBUG
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  if(!myrank) fprintf(stderr, "no = [%td, %td, %td]\n", ths->no[0], ths->no[1], ths->no[2]);
//  if(!myrank) fprintf(stderr, "local_no = [%td, %td, %td], local_ngc = [%td, %td, %td], local_ngc_total = %td\n",
//      local_no[0], local_no[1], local_no[2], local_ngc[0], local_ngc[1], local_ngc[2], local_ngc_total);
      
  PNX(debug_sum_print)(ths->x, 3*ths->local_M, 0,
      "PNFFT^H: Sum of x before sort");
#endif

  /* sort indices for better cache handling */
  if(ths->pnfft_flags & PNFFT_SORT_NODES){
    ths->timer_adj[PNFFT_TIMER_SORT_NODES] -= MPI_Wtime();
    sorted_index = (INT*) PNX(malloc)(sizeof(INT) * (size_t) 2*ths->local_M);
    sort_nodes_for_better_cache_handle(
        ths->d, ths->n, ths->m, ths->local_M, ths->x,
        sorted_index);
    ths->timer_adj[PNFFT_TIMER_SORT_NODES] += MPI_Wtime();
  }
  
#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->x, 3*ths->local_M, 0,
      "PNFFT^H: Sum of x after sort");
  
  PNX(debug_sum_print)(ths->f, ths->local_M,
      !(ths->trafo_flag & PNFFTI_TRAFO_C2R),
      "PNFFT^H: Sum of f");
#endif
  
  ths->timer_adj[PNFFT_TIMER_LOOP_B] -= MPI_Wtime();
  if(interlaced){
    loop_over_particles_adj_interlaced_0(
        ths, local_no_start, local_ngc, gcells_below, sorted_index);
    loop_over_particles_adj_interlaced_1(
        ths, local_no_start, local_ngc, gcells_below, sorted_index);
  } else {
    loop_over_particles_adj(
        ths, local_no_start, local_ngc, gcells_below, sorted_index);
  }
  ths->timer_adj[PNFFT_TIMER_LOOP_B] += MPI_Wtime();

#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->g2, local_ngc_total,
      !(ths->trafo_flag & PNFFTI_TRAFO_C2R),
      "PNFFT^H: Sum of Fourier coefficients before ghostcell reduce");
#endif  

  /* reduce ghost cells in ring */
  ths->timer_adj[PNFFT_TIMER_GCELLS] -= MPI_Wtime();
  PX(reduce)(ths->gcplan);
  ths->timer_adj[PNFFT_TIMER_GCELLS] += MPI_Wtime();

#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->g2, local_no[0]*local_no[1]*local_no[2],
      !(ths->trafo_flag & PNFFTI_TRAFO_C2R),
      "PNFFT^H: Sum of Fourier coefficients after ghostcell reduce");
#endif

  /* perform fftshift */
  ths->timer_adj[PNFFT_TIMER_SHIFT_INPUT] -= MPI_Wtime();
//   if(ths->pnfft_flags & PNFFT_SHIFTED_IN)
  ths->timer_adj[PNFFT_TIMER_SHIFT_INPUT] += MPI_Wtime();

#if PNFFT_ENABLE_DEBUG
  PNX(debug_sum_print)(ths->g2, local_no[0]*local_no[1]*local_no[2],
      !(ths->trafo_flag & PNFFTI_TRAFO_C2R),
      "PNFFT^H: Sum of Fourier coefficients after twiddles");
#endif
  
  if(sorted_index != NULL)
    PNX(free)(sorted_index);
}


static void loop_over_particles_adj(
    PNX(plan) ths, INT *local_no_start, INT *local_ngc, INT *gcells_below,
    INT *sorted_index
    )
{
  const int cutoff = ths->cutoff;
  INT j, m0, u_j[3];
  R floor_nx_j[3];
  R *pre_psi = NULL;
#if PNFFT_ENABLE_DEBUG
  R rsum = 0.0, grsum;
#endif

  if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) )
    pre_psi = (R*) PNX(malloc)(sizeof(R) * (size_t) cutoff*3);
  for(INT p=0; p<ths->local_M; p++){
    j = (sorted_index) ? sorted_index[2*p+1] : p;

    lowest_summation_index(
        ths->n, ths->m, &(ths->x[ths->d*j]), local_no_start, gcells_below,
        floor_nx_j, u_j);

    /* evaluate window on axes */
    if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) ){
      pre_psi_tensor(
          ths->n, ths->b, ths->m, cutoff, &(ths->x[ths->d*j]), floor_nx_j,
          ths->exp_const, ths->spline_coeffs, ths->pnfft_flags,
          ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_psi,
          pre_psi);

#if PNFFT_ENABLE_DEBUG
      /* Don't want to use PNX(debug_sum_print) because we are in a loop */
      for(int t=0; t<3*cutoff; t++)
        rsum += pnfft_fabs(pre_psi[t]);
#endif
    }

    m0 = PNFFT_PLAIN_INDEX_3D(u_j, local_ngc);
    if (ths->trafo_flag & PNFFTI_TRAFO_C2R)
      PNX(spread_f_r2r)(
          ths, p, ths->f[j], pre_psi, m0, local_ngc, cutoff, 1, 0,
          ths->g2);
    else
      PNX(spread_f_c2c)(
          ths, p, ((C*)ths->f)[j], pre_psi, m0, local_ngc, cutoff, 0,
          ths->g2);
  }

#if PNFFT_ENABLE_DEBUG
  MPI_Reduce(&rsum, &grsum, 1, PNFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  PX(fprintf)(MPI_COMM_WORLD, stderr, "PNFFT^H: Sum of pre_psi: %e\n", grsum);
#endif

  if(pre_psi != NULL)
    PNX(free)(pre_psi);
}

static void loop_over_particles_adj_interlaced_0(
    PNX(plan) ths, INT *local_no_start, INT *local_ngc, INT *gcells_below,
    INT *sorted_index
    )
{
  const int cutoff = ths->cutoff;
  INT j, m0, u_j[3];
  R floor_nx_j[3];
  R *pre_psi = NULL;

  if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) )
    pre_psi = (R*) PNX(malloc)(sizeof(R) * (size_t) cutoff*3);
  for(INT p=0; p<ths->local_M; p++){
    j = (sorted_index) ? sorted_index[2*p+1] : p;

    lowest_summation_index(
        ths->n, ths->m, &(ths->x[ths->d*j]), local_no_start, gcells_below,
        floor_nx_j, u_j);

    /* evaluate window on axes */
    if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) ){
      pre_psi_tensor(
          ths->n, ths->b, ths->m, cutoff, &(ths->x[ths->d*j]), floor_nx_j,
          ths->exp_const, ths->spline_coeffs, ths->pnfft_flags,
          ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_psi,
          pre_psi);
    }

    m0 = PNFFT_PLAIN_INDEX_3D(u_j, local_ngc);
    PNX(spread_f_r2r)(
        ths, p, (R)ths->f[j*2], pre_psi, 2*m0, local_ngc, cutoff, 2, 0,
        ths->g2);
  }

  if(pre_psi != NULL)
    PNX(free)(pre_psi);
}


static void loop_over_particles_adj_interlaced_1(
    PNX(plan) ths, INT *local_no_start, INT *local_ngc, INT *gcells_below,
    INT *sorted_index
    )
{
  const int cutoff = ths->cutoff;
  INT j, m0, u_j[3];
  R floor_nx_j[3];
  R *pre_psi = NULL;
  R x[3];

  if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) )
    pre_psi = (R*) PNX(malloc)(sizeof(R) * (size_t) cutoff*3);
  for(INT p=0; p<ths->local_M; p++){
    j = (sorted_index) ? sorted_index[2*p+1] : p;

    /* shift x by half the mesh width */
    for(int t=0; t<3; t++)
      x[t] = ths->x[ths->d*j+t] + 0.5/ths->n[t];

    /* We need to compute the lowest summation index before we fold x back into [-0.5,0.5).
     * Otherwise u_j may be also folded and gets less than the local offset local_no_start. */
    lowest_summation_index(
        ths->n, ths->m, x, local_no_start, gcells_below,
        floor_nx_j, u_j);

    /* assure -0.5 <= x < 0.5 */
    for(int t=0; t<3; t++){
      if(x[t] >= 0.5){
        x[t] -= 1.0;
        floor_nx_j[t] -= ths->n[t];
      }
    }

    /* evaluate window on axes */
    if( !(ths->pnfft_flags & (PNFFT_PRE_PSI | PNFFT_PRE_FULL_PSI)) ){
      pre_psi_tensor(
          ths->n, ths->b, ths->m, cutoff, x, floor_nx_j,
          ths->exp_const, ths->spline_coeffs, ths->pnfft_flags,
          ths->intpol_order, ths->intpol_num_nodes, ths->intpol_tables_psi,
          pre_psi);
    }

    m0 = PNFFT_PLAIN_INDEX_3D(u_j, local_ngc);
    PNX(spread_f_r2r)(
        ths, p, (R)ths->f[j*2], pre_psi, 2*m0+1, local_ngc, cutoff, 2, 1,
        ths->g2);
  }

  if(pre_psi != NULL)
    PNX(free)(pre_psi);
}


void PNX(scale_ik_diff_c2c)(
    const C* g1_buffer, INT *local_N_start, INT *local_N, int dim, unsigned pnfft_flags,
    C* g1
    )
{
  INT k[3], m=0;

  if(pnfft_flags & PNFFT_TRANSPOSED_F_HAT){
    /* g_hat is transposed N1 x N2 x N0 */
    for(k[1]=local_N_start[1]; k[1]<local_N_start[1] + local_N[1]; k[1]++)
      for(k[2]=local_N_start[2]; k[2]<local_N_start[2] + local_N[2]; k[2]++)
        for(k[0]=local_N_start[0]; k[0]<local_N_start[0] + local_N[0]; k[0]++, m++)
          g1[m] = -2*PNFFT_PI * I * k[dim] * g1_buffer[m];
  } else {
    /* g_hat is non-transposed N0 x N1 x N2 */
    for(k[0]=local_N_start[0]; k[0]<local_N_start[0] + local_N[0]; k[0]++)
      for(k[1]=local_N_start[1]; k[1]<local_N_start[1] + local_N[1]; k[1]++)
        for(k[2]=local_N_start[2]; k[2]<local_N_start[2] + local_N[2]; k[2]++, m++)
          g1[m] = -2*PNFFT_PI * I * k[dim] * g1_buffer[m];
  }
}

