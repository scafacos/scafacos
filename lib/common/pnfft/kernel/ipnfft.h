/*
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
 * Copyright (c) 2010-2013 Michael Pippig
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


/* PNFFT internal header file */
#ifndef __IPNFFT_H__
#define __IPNFFT_H__

#include "config.h"

#define PNFFT_ENABLE_DEBUG 0

#define PNFFT_DEBUG_USE_KAISER_BESSEL 0
#define PNFFT_DEBUG_USE_GAUSSIAN      0
#define PNFFT_DEBUG_USE_BSPLINE       0
#define PNFFT_DEBUG_USE_SINC_POWER    0

#define PNFFT_SORT_RADIX 1

/* Begin: This part is based on ifftw3.h */
#include <stdlib.h>             /* size_t */
#include <stdarg.h>             /* va_list */
#include <stddef.h>             /* ptrdiff_t */
#include <stdio.h>              /* fprintf */
#include <memory.h>             /* memset */
#include <float.h>              /* DBL_EPSILON, ... */

#include <pfft.h>

#define IPNFFT_EXTERN extern

typedef ptrdiff_t INT;
/*
  integral type large enough to contain a stride (what ``int'' should
  have been in the first place.
*/

#define IF(x,a,b) ((x)?(a):(b))

#define PNFFT_MIN(a,b) (((a)<(b))?(a):(b))
#define PNFFT_MAX(a,b) (((a)>(b))?(a):(b))
#define PNFFT_ABS(x) (((x)>0)?(x):(-(x)))
#define PNFFT_SIGN(a) (((a)>=0)?1:-1)
#define PNFFT_SQR(a) ((a)*(a))
#define PNFFT_POW3(a) ((a)*(a)*(a))

#define PNFFT_PLAIN_INDEX_3D(k, n)      ( k[2] + n[2]*(k[1] + n[1]*k[0]) )
#define PNFFT_FFTSHIFT(k, N)            ( ((N)/2-1 < (k)) ? (k)-(N) : (k) )

/* all pnfft identifiers start with pnfft (or pnfftf etc.) */
#define CONCAT(prefix, name) prefix ## name

/* define function names according to used precision */
#if defined(PNFFT_SINGLE)
  typedef float R;
  typedef pnfftf_complex C;
#  define PNFFT_MPI_REAL_TYPE MPI_FLOAT
#  define PNX(name) CONCAT(pnfftf_, name)
#  define PX(name) PFFT_MANGLE_FLOAT(name)
#  define X(name)  FFTW_MANGLE_FLOAT(name)
#  define PNFFT_EPSILON LDBL_EPSILON//4.0E-31L
#  define PNFFT_MATH(name) CONCAT(name, f)
#elif defined(PNFFT_LDOUBLE)
  typedef long double R;
  typedef pnfftl_complex C;
#  define PNFFT_MPI_REAL_TYPE MPI_LONG_DOUBLE
#  define PNX(name) CONCAT(pnfftl_, name)
#  define PX(name) PFFT_MANGLE_LONG_DOUBLE(name)
#  define X(name)  FFTW_MANGLE_LONG_DOUBLE(name)
#  define PNFFT_EPSILON FLT_EPSILON
#  define PNFFT_MATH(name) CONCAT(name, l)
#else
  typedef double R;
  typedef pnfft_complex C;
#  define PNFFT_MPI_REAL_TYPE MPI_DOUBLE
#  define PNX(name) CONCAT(pnfft_, name)
#  define PX(name) PFFT_MANGLE_DOUBLE(name)
#  define X(name)  FFTW_MANGLE_DOUBLE(name)
#  define PNFFT_EPSILON DBL_EPSILON
#  define PNFFT_MATH(name) name
#endif

/* macros for mathematical functions corresponding to the PNFFT float data type */
#  define pnfft_pow(_x_, _y_)   PNFFT_MATH(pow)(_x_, _y_)
#  define pnfft_exp(_x_)   PNFFT_MATH(exp)(_x_)
#  define pnfft_sqrt(_x_)  PNFFT_MATH(sqrt)(_x_)
#  define pnfft_cexp(_x_)  PNFFT_MATH(cexp)(_x_) 
#  define pnfft_creal(_x_) PNFFT_MATH(creal)(_x_) 
#  define pnfft_cimag(_x_) PNFFT_MATH(cimag)(_x_) 
#  define pnfft_floor(_x_) PNFFT_MATH(floor)(_x_) 
#  define pnfft_ceil(_x_)  PNFFT_MATH(ceil)(_x_) 
#  define pnfft_lrint(_x_) PNFFT_MATH(lrint)(_x_) 
#  define pnfft_sin(_x_)   PNFFT_MATH(sin)(_x_) 
#  define pnfft_sinh(_x_)  PNFFT_MATH(sinh)(_x_) 
#  define pnfft_cos(_x_)   PNFFT_MATH(cos)(_x_) 
#  define pnfft_cosh(_x_)  PNFFT_MATH(cosh)(_x_) 
#  define pnfft_tan(_x_)   PNFFT_MATH(tan)(_x_) 
#  define pnfft_fabs(_x_)  PNFFT_MATH(fabs)(_x_) 
#  define pnfft_log2(_x_)  PNFFT_MATH(log2)(_x_) 

/* define MPI-FFTW name mangeling  */
#define XM(name)  X(CONCAT(mpi_, name))

/* Need these defines for compilation of bessel_i0.c, sinc.c (copied from NFFT) */
#ifdef PNFFT_LDOUBLE
#  define K(x) ((R) x##L)
#else
#  define K(x) ((R) x)
#endif

#define A(ex) /* nothing */

#define PNFFT_PRINT_TIMER_BASIC    (1U<<0)
#define PNFFT_PRINT_TIMER_ADV      (1U<<1)


struct PNX(plan_s){                                                                      
  INT N_total;                /**< Total number of Fourier coefficients            */
  C *f_hat;                   /**< Vector of Fourier coefficients                  */
  C *f;                       /**< Vector of samples                               */
  C *grad_f;                  /**< Vector of gradients                             */
  R *x;                       /**< Nodes in time/spatial domain                    */
                                                                                     
  int d;                      /**< Dimension, rank                                 */
  INT *N;                     /**< Multi bandwidth                                 */
  R *sigma;                   /**< Oversampling-factor                             */
  INT *n;                     /**< FFT length, equal to sigma*N                    */
  INT *no;                    /**< FFT output length                               */
  INT n_total;                /**< Total size of FFTW                              */
  int m;                      /**< Cut-off parameter of the window function        */
  R *x_max;                   /**< Upper border for nodes in time/spatial domain   */
  INT *local_N;               /**< Local multi bandwidth                           */
  INT *local_N_start;         /**< Offset of local multi bandwidth                 */
  INT local_N_total;          /**< Total number of local Fourier coefficients      */
  INT *local_no;              /**< Local FFT output length                         */
  INT *local_no_start;        /**< Offset FFT output offset                        */
  INT local_no_total;         /**< Total number of local FFT outputs               */
                                                                                     
  R *b;                       /**< Shape parameter of Gaussian window function     */
  R *exp_const;               /**< Precomputed values for Fast Gaussian window     */
                                                                                     
  R *spline_coeffs;           /**< Input for de Boor algorithm, if B_SPLINE or       
                                   SINC_POWER is used                              */
                                                                                     
  C *pre_inv_phi_hat_trafo;   /**< Precomputed inverse window Fourier coefficients   
                                   for NFFT trafo                                  */
  C *pre_inv_phi_hat_adj;     /**< Precomputed inverse window Fourier coefficients   
                                   for adjoint NFFT                                */

  R *pre_psi;                 /**< Precomputed window function values              */
  R *pre_dpsi;                /**< Precomputed window function derivatives         */
  R *pre_psi_il;              /**< Precomputed window function values, interlaced  */
  R *pre_dpsi_il;             /**< Precomputed window function derivatives, interlaced */
                                                                                     
  unsigned pnfft_flags;        /**< Flags for precomputation, (de)allocation,        
                                   and FFTW usage                                  */
  unsigned compute_flags;     /**< Flags for choice of NFFT results                */
  unsigned pfft_opt_flags;    /**< Flags for PFFT optimization                     */
                                                                                     
  /* internal*/                                                                      
  PX(plan)   pfft_forw;       /**< Forward PFFT plan                               */
  PX(plan)   pfft_back;       /**< Backward PFFT plan                              */
  PX(gcplan) gcplan;          /**< PFFT Ghostcell plan                             */
                                                                                     
  C *g1;                      /**< Input of PFFT                                   */
  C *g2;                      /**< Output of PFFT                                  */
  C *g1_buffer;               /**< Buffer for computing Fourier-space derivatives  */
                                                                                     
  int cutoff;                 /**< cutoff range                                    */
  INT local_M;                /**< Number of local nodes                           */
                                                                                     
  /* parameters for window interpolation table */                                    
  int intpol_order;           /**< order of window interpolation                   */
  INT intpol_num_nodes;       /**< number of sampled points for interpolation      */
  R **intpol_tables_psi;      /**< sampled values of window functions              */
  R **intpol_tables_dpsi;     /**< sampled values of window function derivatives   */
                                                                                     
  MPI_Comm comm_cart;         /**< 2d or 3d Cartesian communicator                 */
  int np[3];                  /**< Size of Cartesian communicator                  */
  int rnk_pm;                 /**< rank of Cartesian communicator                  */
  int coords[3];              /**< 3d coordinates within Cartesian communicator    */
                                                                                     
  double* timer_trafo;        /**< Saves time measurements during PNFFT            */
  double* timer_adj;          /**< Saves time measurements during adjoint PNFFT    */
};
typedef struct PNX(plan_s) plan_s;

#ifndef PNFFT_H
typedef struct PNX(plan_s) *PNX(plan);
#endif /* !PNFFT_H */



/* util.c */
void PNX(message)(
    const char *string, MPI_Comm comm);
INT PNX(prod_INT)(
    int d, const INT *vec);
INT PNX(sum_INT)(
    int d, const INT *vec);
int PNX(equal_INT)(
    int d, const INT *vec1, const INT *vec2);
void PNX(vcopy_INT)(
    int d, const INT *vec1,
    INT *vec2);
void PNX(vadd_INT)(
    int d, const INT *vec1, const INT *vec2,
    INT *sum);
void PNX(vsub_INT)(
    int d, const INT *vec1, const INT *vec2,
    INT *sum);

INT* PNX(malloc_INT)(
    size_t size);
int* PNX(malloc_int)(
    size_t size);
unsigned* PNX(malloc_unsigned)(
    size_t size);
R* PNX(malloc_R)(
    size_t size);
C* PNX(malloc_C)(
    size_t size);

void PNX(sort_node_indices_radix_lsdf)(
    INT n, INT *keys0, INT *keys1, INT rhigh);
void PNX(sort_node_indices_radix_msdf)(
    INT n, INT *keys0, INT *keys1, INT rhigh);
void PNX(sort_nodes_indices_qsort_3d)(
    const INT *n, int m, INT local_M, const R *x,
    INT *sort);


/* malloc.c */
void PNX(die)(
    const char *s, MPI_Comm comm);

/* timer.c */
double* PNX(mktimer)(
    void);
void PNX(rmtimer)(
    double* timer);

/* ndft-parallel.c */
void PNX(rmplan)(
    PNX(plan) ths);
INT PNX(local_size_internal)(
    const INT *N, const INT *n, const INT *no,
    MPI_Comm comm_cart_2d, unsigned pnfft_flags,
    INT *local_N, INT *local_N_start,
    INT *local_no, INT *local_no_start);
PNX(plan) PNX(init_internal)(
    int d, const INT *N, const INT *n, const INT *no,
    INT local_M, int m,
    unsigned pnfft_flags, unsigned pfft_opt_flags,
    MPI_Comm comm_cart_2d);
void PNX(trafo_F)(
    PNX(plan) ths);
void PNX(adjoint_F)(
    PNX(plan) ths);
void PNX(trafo_B_grad_ad)(
    PNX(plan) ths, int interlaced);
void PNX(trafo_B_grad_ik)(
    PNX(plan) ths, C *f, INT offset, INT stride);
void PNX(adjoint_B)(
    PNX(plan) ths, int interlaced);
void PNX(malloc_x)(
    PNX(plan) ths, unsigned pnfft_flags);
void PNX(free_x)(
    PNX(plan) ths, unsigned pnfft_finalize_flags);
void PNX(malloc_f)(
    PNX(plan) ths, unsigned pnfft_flags);
void PNX(free_f)(
    PNX(plan) ths, unsigned pnfft_finalize_flags);
void PNX(malloc_grad_f)(
    PNX(plan) ths, unsigned pnfft_flags);
void PNX(free_grad_f)(
    PNX(plan) ths, unsigned pnfft_finalize_flags);
void PNX(node_borders)(
    const INT *n,
    const INT *local_no, const INT *local_no_start,
    const R* x_max,
    R *lo, R *up);
void PNX(scale_ik_diff_c2c)(
    const C* g1_buffer, INT *local_N_start, INT *local_N, int dim, unsigned pnfft_flags,
    C* g1);

/* assgin.c */
void PNX(spread_f_c2c)(
    PNX(plan) ths, INT ind,
    C f, R *pre_psi, INT m0, INT *local_ngc, int cutoff,
    int interlaced,
    C *grid);
void PNX(spread_f_r2r)(
    PNX(plan) ths, INT ind,
    R f, R *pre_psi, INT m0, INT *local_ngc, int cutoff, INT ostride,
    int interlaced,
    R *grid);
void PNX(assign_f_c2c)(
    PNX(plan) ths, INT ind,
    C *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff, int interlaced,
    C *f);
void PNX(assign_f_r2r)(
    PNX(plan) ths, INT ind,
    R *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff, INT istride, int interlaced,
    R *f);
void PNX(assign_grad_f_c2c)(
    PNX(plan) ths, INT ind,
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, int interlaced,
    C *grad_f);
void PNX(assign_grad_f_r2r)(
    PNX(plan) ths, INT ind,
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    INT istride, INT ostride, int interlaced,
    R *grad_f);
void PNX(assign_f_and_grad_f_c2c)(
    PNX(plan) ths, INT ind,
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, int interlaced,
    C *f, C *grad_f);
void PNX(assign_f_and_grad_f_r2r)(
    PNX(plan) ths, INT ind,
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    INT istride, INT ostride, int interlaced,
    R *f, R *grad_f);



void PNX(spread_f_c2c_pre_psi)(
    C f, R *pre_psi, INT m0, INT *grid_size, int cutoff,
    C *grid);
void PNX(spread_f_c2c_pre_full_psi)(
    C f, R *pre_psi, INT m0, INT *grid_size, int cutoff,
    C *grid);
void PNX(spread_f_r2r_pre_psi)(
    R f, R *pre_psi, INT m0, INT *grid_size, int cutoff, INT ostride,
    R *grid);
void PNX(spread_f_r2r_pre_full_psi)(
    R f, R *pre_psi, INT m0, INT *grid_size, int cutoff, INT ostride,
    R *grid);

void PNX(assign_f_c2c_pre_psi)(
    C *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff,
    C *fv);
void PNX(assign_f_c2c_pre_full_psi)(
    C *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff,
    C *fv);
void PNX(assign_f_r2r_pre_psi)(
    R *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff, INT istride,
    R *fv);
void PNX(assign_f_r2r_pre_full_psi)(
    R *grid, R *pre_psi, INT m0, INT *grid_size, int cutoff, int istride,
    R *fv);

void PNX(assign_grad_f_c2c_pre_psi)(
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    C *grad_f);
void PNX(assign_grad_f_c2c_pre_full_psi)(
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    C *grad_f);
void PNX(assign_grad_f_r2r_pre_psi)(
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, INT istride, INT ostride,
    R *grad_f);
void PNX(assign_grad_f_r2r_pre_full_psi)(
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, INT istride, INT ostride,
    R *grad_f);

void PNX(assign_f_and_grad_f_c2c_pre_psi)(
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    C *fv, C *grad_f);
void PNX(assign_f_and_grad_f_c2c_pre_full_psi)(
    C *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff,
    C *fv, C *grad_f);
void PNX(assign_f_and_grad_f_r2r_pre_psi)(
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, INT istride, INT ostride,
    R *fv, R *grad_f);
void PNX(assign_f_and_grad_f_r2r_pre_full_psi)(
    R *grid, R *pre_psi, R *pre_dpsi, INT m0, INT *grid_size, int cutoff, INT istride, INT ostride,
    R *fv, R *grad_f);




/** constant spline interpolation */
static inline R pnfft_intpol_const(
    int k, const R *table
    )
{
  return table[k];
}

/** linear spline interpolation */
static inline R pnfft_intpol_lin(
    int k, R dist, const R *table
    )
{
  R f1,f2;
  f1=table[k]; f2=table[k+1];
  return f1*(1.0-dist) + f2*dist;
}

/** quadratic spline interpolation */
static inline R pnfft_intpol_quad(
    int k, R dist, const R *table
    )
{
  R c0,c1,c2,f0,f1,f2;
  f0=table[k]; f1=table[k+1]; f2=table[k+2];
  c0=dist+1.0;
  c1=dist;
  c2=dist-1.0;
  return (0.5*f0*c1*c2-c0*f1*c2+0.5*c0*c1*f2);
}

/** cubic spline interpolation */
static inline R pnfft_intpol_kub(
    int k, R dist, const R *table
    )
{
  R c0,c1,c2,c3;
  R f0,f1,f2,f3;
  f0=table[k]; f1=table[k+1]; f2=table[k+2]; f3=table[k+3];
  c0=dist+1.0;
  c1=dist;
  c2=dist-1.0;
  c3=dist-2.0;
  return(-f0*c1*c2*c3+3.0*c0*f1*c2*c3-3.0*c0*c1*f2*c3+c0*c1*c2*f3)/6.0;
}

#endif /* __IPNFFT_H__ */
