/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (C) 2011 Sebastian Banert
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

#ifndef _P2NFFT_TYPES_H
#define _P2NFFT_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <complex.h>
#include "pnfft.h"
#include <float.h>

#include "common/gridsort/gridsort_resort.h"


#define FCS_P2NFFT_DEBUG 0
#define FCS_P2NFFT_DEBUG_RETUNE 0
#define FCS_P2NFFT_TIMING 0

#define FCS_P2NFFT_NORMALIZED_2DP_EWALD 0

#define CONCAT2(prefix, name)  prefix ## name
#define CONCAT(prefix, name)  CONCAT2(prefix, name)

// typedef CONCAT(FFTW_PREFIX, fftw_complex) C;
// #define FFTW(name) FFTW_MANGLE_DOUBLE(name)

#if defined(FCS_FLOAT_IS_DOUBLE)
typedef pnfft_complex fcs_pnfft_complex;
# define FCS_PNFFT(name)  PNFFT_MANGLE_DOUBLE(name)
# define FCS_PFFT(name)   PFFT_MANGLE_DOUBLE(name)
# define FCS_P2NFFT_EPS   DBL_EPSILON
#elif defined(FCS_FLOAT_IS_FLOAT)
typedef pnfftf_complex fcs_pnfft_complex;
# define FCS_PNFFT(name)  PNFFT_MANGLE_FLOAT(name)
# define FCS_PFFT(name)   PFFT_MANGLE_FLOAT(name)
# define FCS_P2NFFT_EPS   FLT_EPSILON
#elif defined(FCS_FLOAT_IS_LONG_DOUBLE)
typedef pnfftl_complex fcs_pnfft_complex;
# define FCS_PNFFT(name)  PNFFT_MANGLE_LONG_DOUBLE(name)
# define FCS_PFFT(name)   PFFT_MANGLE_LONG_DOUBLE(name)
# define FCS_P2NFFT_EPS   LDBL_EPSILON
#else
# error "fcs_float is neither double, float nor long double"
#endif

#define FCS_P2NFFT_SQR(x)  ((x) * (x))

typedef pfft_complex C;
typedef ptrdiff_t INT;

/* We always check for the value of macros. Avoid warnings about empty macro by setting it to 0. */
#ifndef FCS_ENABLE_DEBUG
#  define FCS_ENABLE_DEBUG 0
#endif

#define FCS_P2NFFT_REG_NEAR_DEFAULT     (-1)
#define FCS_P2NFFT_REG_NEAR_CG           0
#define FCS_P2NFFT_REG_NEAR_T2P          1

#define FCS_P2NFFT_REG_FAR_DEFAULT      (-1)
#define FCS_P2NFFT_REG_FAR_RAD_CG        0
#define FCS_P2NFFT_REG_FAR_RAD_T2P_SYM   1
#define FCS_P2NFFT_REG_FAR_RAD_T2P_EC    2
#define FCS_P2NFFT_REG_FAR_RAD_T2P_IC    3
#define FCS_P2NFFT_REG_FAR_REC_T2P_SYM   4
#define FCS_P2NFFT_REG_FAR_REC_T2P_EC    5
#define FCS_P2NFFT_REG_FAR_REC_T2P_IC    6

/* p2nfft_flags */
#define FCS_P2NFFT_CHECK_TOLERANCE           (0U)
#define FCS_P2NFFT_IGNORE_TOLERANCE          (1U << 0)

#define FCS_P2NFFT_DEFAULT_TOLERANCE         0.01

#define FCS_P2NFFT_DEFAULT_PNFFT_WINDOW  1 /* Bspline */
#define FCS_P2NFFT_DEFAULT_PFFT_PATIENCE 1 /* Measure */

#if FCS_ENABLE_TIMING || FCS_P2NFFT_TIMING
#define FCS_P2NFFT_INIT_TIMING(comm) \
  int tm_rank; \
  MPI_Comm_rank(comm, &tm_rank); \
  double tm_timer, tm_global_timer;
#define FCS_P2NFFT_START_TIMING() \
  tm_timer = -MPI_Wtime();
#define FCS_P2NFFT_FINISH_TIMING(comm, str) \
  tm_timer += MPI_Wtime(); \
  MPI_Reduce(&tm_timer, &tm_global_timer, 1, MPI_DOUBLE, MPI_MAX, 0, comm); \
  if(!tm_rank) printf("P2NFFT_TIMING: %s takes %e s\n", str, tm_global_timer);
#else
#define FCS_P2NFFT_INIT_TIMING(comm)
#define FCS_P2NFFT_START_TIMING()
#define FCS_P2NFFT_FINISH_TIMING(comm, str)
#endif

typedef struct {
  MPI_Comm cart_comm_pnfft;   /**< @brief Cartesian communicator used as a
                                processor grid for PNFFT. */
  MPI_Comm cart_comm_3d;      /**< @brief 3d Cartesian communicator used as a
                                processor grid for near field. */
  int np[3];                  /**< @brief Procmesh size */
  FCS_PNFFT(plan) pnfft;      /**< @brief Pointer to the PNFFT plan used
                                for all NDFTs. */
  ptrdiff_t N[3];             /**< @brief The dimensions of the frequency grid. */
  ptrdiff_t local_N[3];       /**< @brief Local dimensions of the grid. */
  ptrdiff_t local_N_start[3]; /**< @brief Local offsets of the grid. */
  fcs_float lower_border[3];
  fcs_float upper_border[3];
  fcs_int tolerance_type;
  fcs_float tolerance;
  fcs_int needs_retune;
  fcs_int num_nodes;
  fcs_int tune_alpha;
  fcs_int tune_r_cut;
  fcs_int tune_epsI;
  fcs_int tune_epsB;
  fcs_int tune_N;
  fcs_int tune_n;
  fcs_int tune_m;
  fcs_int tune_p;
  fcs_int tune_c;
  fcs_int sum_qpart;
  fcs_float sum_q;
  fcs_float sum_q2;
  fcs_float bg_charge;
  fcs_float box_l[3];

  /* introduce box shifts for non-cubic boxes */
  fcs_float box_scales[3];
  fcs_float box_shifts[3];
  
  /* P2NFFT flags */
  unsigned flags;

  fcs_int m;

  unsigned pnfft_flags;
  fcs_int  pnfft_interpolation_order;
  fcs_int  pnfft_window;
  unsigned pfft_flags;
  fcs_int  pfft_patience;

  fcs_int short_range_flag;
  fcs_int reg_near;
  fcs_int reg_far;

  fcs_int periodicity[3];

  /* NFFT specific parameters */
  ptrdiff_t n[3];
  fcs_float x_max[3];
  
  /* Near field parameters */
  fcs_int use_ewald;  /* switch between fully periodic and non-periodic case */
  fcs_float r_cut;    /* near field radius (unscaled) */
  fcs_float one_over_r_cut;    /* inverse near field radius (unscaled) */
  fcs_int num_nonperiodic_dims; /* number of dimensions with nonperiodic boundary conditions */
  fcs_int num_periodic_dims;    /* number of dimensions with periodic boundary conditions */

  /* parameters for periodic case */
  fcs_float alpha;  /* Ewald splitting parameter */

  /* parameters for non-periodic case */
  fcs_float epsI;  /* near field size scaled into unit cube */
  fcs_int log2epsI;
  fcs_float epsB;  /* size of regualization border scaled into unit cube */
  fcs_int log2epsB;
  fcs_int p;
  fcs_float c;
  fcs_float *taylor2p_coeff;
  fcs_float *taylor2p_derive_coeff;

  /* parameters for interpolation table */
  fcs_int interpolation_order;      /* interpolation order */
  fcs_int near_interpolation_num_nodes;  /* number of sampled points in near field */
  fcs_int far_interpolation_num_nodes;   /* number of sampled points in far field */
  fcs_float *near_interpolation_table_potential; /* nearfield potential values */
  fcs_float *near_interpolation_table_force;     /* nearfield force values */
  fcs_float *far_interpolation_table_potential;  /* potential values between far field border and 0.5 */

  /* Cosine coefficients of CG-optimized regularization */
  fcs_int N_cg_cos;
  fcs_float *cg_cos_coeff;
  fcs_float *cg_sin_coeff;
  
  /* Fourier coefficients of regularized kernel approximation */
  fcs_pnfft_complex *regkern_hat;

  /* array to store the virial matrix */
  fcs_float *virial;

  /* resort parameters */
  fcs_float max_particle_move;
  fcs_int resort, local_num_particles;
  fcs_gridsort_resort_t gridsort_resort;

  /* gridsort cache */
  fcs_gridsort_cache_t gridsort_cache;

} ifcs_p2nfft_data_struct;

#endif
