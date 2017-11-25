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
#include <pnfft.h>
#include <float.h>

#include "common/gridsort/gridsort.h"
#include "common/gridsort/gridsort_resort.h"


#define FCS_P2NFFT_DEBUG 0
#define FCS_P2NFFT_DEBUG_RETUNE 0
#define FCS_P2NFFT_TIMING 0

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
#  define FCS_P2NFFT_IFDBG(code) code
#else
#  define FCS_P2NFFT_IFDBG(code) 
#endif

#if FCS_P2NFFT_DEBUG_RETUNE
#  define FCS_P2NFFT_IFDBG_RETUNE(code) code
#else
#  define FCS_P2NFFT_IFDBG_RETUNE(code)
#endif


#define CONCAT2(prefix, name)  prefix ## name
#define CONCAT(prefix, name)  CONCAT2(prefix, name)

// typedef CONCAT(FFTW_PREFIX, fftw_complex) C;
// #define FFTW(name) FFTW_MANGLE_DOUBLE(name)

#if defined(FCS_FLOAT_IS_DOUBLE)
typedef pnfft_complex fcs_pnfft_complex;
# define FCS_PNFFT(name)  PNFFT_MANGLE_DOUBLE(name)
# define FCS_PFFT(name)   PFFT_MANGLE_DOUBLE(name)
# define FCS_P2NFFT_EPS   DBL_EPSILON
# define FCS_P2NFFT_ALPHA_OPT_PREC 1.0e-10
#elif defined(FCS_FLOAT_IS_FLOAT)
typedef pnfftf_complex fcs_pnfft_complex;
# define FCS_PNFFT(name)  PNFFT_MANGLE_FLOAT(name)
# define FCS_PFFT(name)   PFFT_MANGLE_FLOAT(name)
# define FCS_P2NFFT_EPS   FLT_EPSILON
# define FCS_P2NFFT_ALPHA_OPT_PREC 1.0e-5
#elif defined(FCS_FLOAT_IS_LONG_DOUBLE)
typedef pnfftl_complex fcs_pnfft_complex;
# define FCS_PNFFT(name)  PNFFT_MANGLE_LONG_DOUBLE(name)
# define FCS_PFFT(name)   PFFT_MANGLE_LONG_DOUBLE(name)
# define FCS_P2NFFT_EPS   LDBL_EPSILON
# define FCS_P2NFFT_ALPHA_OPT_PREC 1.0e-15
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
#define FCS_P2NFFT_IGNORE_POTENTIAL          (1U << 1)
#define FCS_P2NFFT_IGNORE_FIELD              (1U << 2)
#define FCS_P2NFFT_VERBOSE_TUNING            (1U << 3)
#define FCS_P2NFFT_MAX_BOX_ANGLES            (1U << 4)

/* p2nfft_reg_kernels */
#define FCS_P2NFFT_REG_KERNEL_EWALD             0
#define FCS_P2NFFT_REG_KERNEL_ONE_OVER_ABS_X    1
#define FCS_P2NFFT_REG_KERNEL_DEFAULT           FCS_P2NFFT_REG_KERNEL_EWALD

#define FCS_P2NFFT_DEFAULT_TOLERANCE         0.01

#define FCS_P2NFFT_DEFAULT_PNFFT_WINDOW  1 /* Bspline */
#define FCS_P2NFFT_DEFAULT_PFFT_PATIENCE 1 /* Measure */

#if FCS_ENABLE_TIMING || FCS_P2NFFT_TIMING
#define FCS_P2NFFT_INIT_TIMING(comm) \
  int tm_rank; \
  MPI_Comm_rank(comm, &tm_rank); \
  double tm_timer, tm_global_timer_max, tm_global_timer_min;
#define FCS_P2NFFT_START_TIMING(comm) \
  MPI_Barrier(comm); \
  tm_timer = -MPI_Wtime();
#define FCS_P2NFFT_FINISH_TIMING(comm, str) \
  tm_timer += MPI_Wtime(); \
  MPI_Reduce(&tm_timer, &tm_global_timer_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm); \
  MPI_Reduce(&tm_timer, &tm_global_timer_min, 1, MPI_DOUBLE, MPI_MIN, 0, comm); \
  if(!tm_rank) printf("P2NFFT_TIMING: %s takes %.2e s\n", str, tm_global_timer_max); \
  if(!tm_rank) printf("               (minimum %.2e s, load imbalance %.2f %%)\n", tm_global_timer_min, 100.0*(tm_global_timer_max-tm_global_timer_min)/tm_global_timer_min);
#else
#define FCS_P2NFFT_INIT_TIMING(comm)
#define FCS_P2NFFT_START_TIMING(comm)
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
  FCS_PNFFT(nodes) charges;   /**< @brief Pointer to the PNFFT nodes of charges. */
  FCS_PNFFT(nodes) dipoles;   /**< @brief Pointer to the PNFFT nodes of charges. */

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
  fcs_int tune_k_cut;
  fcs_int tune_N;
  fcs_int tune_n;
  fcs_int tune_m;
  fcs_int tune_p;
  fcs_int tune_b;
  fcs_int tune_c;
  fcs_int sum_qpart;
  fcs_float sum_q;
  fcs_float sum_q2;
  fcs_float bg_charge;

  /* TODO: deprecated, only valid for orthorombic boxes */
  fcs_float box_l[3];


  /* pricipal axes of triclinic box */
  fcs_float box_a[3];
  fcs_float box_b[3];
  fcs_float box_c[3];
  fcs_float box_base[3];    /* box offset */
  fcs_float box_V;          /* volume of unit cell */

  /* pricipal axes of dual lattice (inverse box) */
  fcs_float box_inv[9];

  /* principal axes of the extended and regularized triclinic box */
  fcs_float ebox_a[3];
  fcs_float ebox_b[3];
  fcs_float ebox_c[3];
  fcs_float ebox_V;

  /* pricipal axes of inverse extended box */
  fcs_float ebox_inv[9];

  /* box to ebox expansion factors */
  fcs_float box_expand[3];

  /* box_scales = box_l * box_expand (unit box to ebox expansion factors) */
  fcs_float box_scales[3];
  
  /* P2NFFT flags */
  unsigned flags;

  fcs_int m;

  unsigned pnfft_flags;
  fcs_int  pnfft_interpolation_order;
  fcs_int  pnfft_window;
  fcs_int  pnfft_direct;
  unsigned pfft_flags;
  fcs_int  pfft_patience;

  fcs_int short_range_flag;
  fcs_int reg_near;
  fcs_int reg_far;
  fcs_int reg_kernel;

  fcs_int periodicity[3];

  /* NFFT specific parameters */
  ptrdiff_t n[3];
  fcs_float x_max[3];
  fcs_float b[3];
  
  /* Near field parameters */
  fcs_float r_cut;    /* near field radius (unscaled) */
  fcs_float one_over_r_cut;    /* inverse near field radius (unscaled) */
  fcs_int num_nonperiodic_dims; /* number of dimensions with nonperiodic boundary conditions */
  fcs_int num_periodic_dims;    /* number of dimensions with periodic boundary conditions */

  /* parameters for periodic case */
  fcs_float alpha;  /* Ewald splitting parameter */
  fcs_float k_cut;    /* near field radius (unscaled) */

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
  fcs_int resort, local_num_particles, local_num_dipole_particles;
  fcs_gridsort_resort_t gridsort_resort;

  /* gridsort cache */
  fcs_gridsort_cache_t gridsort_cache;

} ifcs_p2nfft_data_struct;

#endif
