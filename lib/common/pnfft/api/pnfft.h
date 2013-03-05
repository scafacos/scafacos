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

#ifndef PNFFT_H
#define PNFFT_H 1

#include <math.h>
#include <mpi.h>
#include <pfft.h>

#ifdef __cplusplus
  extern "C"
  {
#endif /* __cplusplus */


#define PNFFT_EXTERN extern
#define PNFFT_CONCAT(prefix, name) prefix ## name


  /*
    huge second-order macro that defines prototypes for all API
    functions.  We expand this macro for each supported precision

    PNX: PNFFT name-mangling macro (parallel NFFT)
    PX: PFFT name-mangling macro (parallel FFT)
    X : FFTW name-mangling macro (serial FFT)
    R : real data type
    C : complex data type
  */
#define PNFFT_DEFINE_API(PNX, PX, X, R, C, INT)                                         \
                                                                                        \
  typedef struct PNX(plan_s) *PNX(plan);                                                \
                                                                                        \
  PNFFT_EXTERN int PNX(create_procmesh_2d)(                                             \
      MPI_Comm comm, int np0, int np1, MPI_Comm *comm_cart_2d);                         \
  PNFFT_EXTERN int PNX(create_procmesh)(                                                \
      int rnk, MPI_Comm comm, const int *np, MPI_Comm *comm_cart);                      \
                                                                                        \
  PNFFT_EXTERN void PNX(local_size_3d)(                                                 \
      const INT *N, MPI_Comm comm_cart,                                                 \
      INT *local_N, INT *local_N_start,                                                 \
      R *lower_border, R *upper_border);                                                \
  PNFFT_EXTERN void PNX(local_size_adv)(                                                \
      int d, const INT *N, MPI_Comm comm_cart,                                          \
      INT *local_N, INT *local_N_start,                                                 \
      R *lower_border, R *upper_border);                                                \
  PNFFT_EXTERN void PNX(local_size_guru)(                                               \
      int d, const INT *N, const INT *n, const R *x_max, int m,                         \
      MPI_Comm comm_cart,                                                               \
      INT *local_N, INT *local_N_start,                                                 \
      R *lower_border, R *upper_border);                                                \
                                                                                        \
  PNFFT_EXTERN PNX(plan) PNX(init_3d)(                                                  \
      const INT *N, INT local_M,                                                        \
      MPI_Comm comm_cart);                                                              \
  PNFFT_EXTERN PNX(plan) PNX(init_adv)(                                                 \
      int d, const INT *N, INT local_M,                                                 \
      unsigned pnfft_flags, unsigned fftw_flags, MPI_Comm comm_cart);                   \
  PNFFT_EXTERN PNX(plan) PNX(init_guru)(                                                \
        int d,                                                                          \
        const INT *N, const INT *n, const R *x_max,                                     \
        INT local_M, int m,                                                             \
        unsigned pnfft_flags, unsigned fftw_flags,                                      \
        MPI_Comm comm_cart);                                                            \
                                                                                        \
  PNFFT_EXTERN void PNX(init_nodes)(                                                    \
      PNX(plan) ths, INT local_M,                                                       \
      unsigned pnfft_flags, unsigned pnfft_finalize_flags);                             \
                                                                                        \
  PNFFT_EXTERN void PNX(precompute_psi)(                                                \
      PNX(plan) ths);                                                                   \
                                                                                        \
  PNFFT_EXTERN void PNX(set_f_hat)(                                                     \
      C *f_hat, PNX(plan) ths);                                                         \
  PNFFT_EXTERN void PNX(set_f)(                                                         \
      C *f, PNX(plan) ths);                                                             \
  PNFFT_EXTERN void PNX(set_grad_f)(                                                    \
      C *grad_f, PNX(plan) ths);                                                        \
  PNFFT_EXTERN void PNX(set_x)(                                                         \
      R *x, PNX(plan) ths);                                                             \
                                                                                        \
  PNFFT_EXTERN C* PNX(get_f_hat)(                                                       \
      const PNX(plan) ths);                                                             \
  PNFFT_EXTERN C* PNX(get_f)(                                                           \
      const PNX(plan) ths);                                                             \
  PNFFT_EXTERN C* PNX(get_grad_f)(                                                      \
      const PNX(plan) ths);                                                             \
  PNFFT_EXTERN R* PNX(get_x)(                                                           \
      const PNX(plan) ths);                                                             \
                                                                                        \
  PNFFT_EXTERN int PNX(get_d)(                                                          \
      const PNX(plan) ths);                                                             \
  PNFFT_EXTERN int PNX(get_m)(                                                          \
      const PNX(plan) ths);                                                             \
  PNFFT_EXTERN R* PNX(get_x_max)(                                                       \
      const PNX(plan) ths);                                                             \
  PNFFT_EXTERN INT* PNX(get_N)(                                                         \
      const PNX(plan) ths);                                                             \
  PNFFT_EXTERN INT* PNX(get_n)(                                                         \
      const PNX(plan) ths);                                                             \
  PNFFT_EXTERN unsigned PNX(get_pnfft_flags)(                                           \
      const PNX(plan) ths);                                                             \
  PNFFT_EXTERN unsigned PNX(get_pfft_flags)(                                            \
      const PNX(plan) ths);                                                             \
                                                                                        \
  PNFFT_EXTERN void PNX(finalize)(                                                      \
        PNX(plan) ths, unsigned pnfft_finalize_flags);                                  \
                                                                                        \
  PNFFT_EXTERN void PNX(trafo)(                                                         \
      PNX(plan) ths);                                                                   \
  PNFFT_EXTERN void PNX(adj)(                                                           \
      PNX(plan) ths);                                                                   \
                                                                                        \
  PNFFT_EXTERN void PNX(init)(                                                          \
      void);                                                                            \
  PNFFT_EXTERN void PNX(cleanup)(                                                       \
      void);                                                                            \
                                                                                        \
  PNFFT_EXTERN void *PNX(malloc)(size_t n);					        \
  PNFFT_EXTERN R *PNX(alloc_real)(size_t n);					        \
  PNFFT_EXTERN C *PNX(alloc_complex)(size_t n);				                \
  PNFFT_EXTERN void PNX(free)(void *p);					                \
                                                                                        \
  PNFFT_EXTERN void PNX(init_f_hat_3d)(                                                 \
      const INT *N, const INT *local_N,                                                 \
      const INT *local_N_start,                                                         \
      C *data);                                                                         \
  PNFFT_EXTERN void PNX(init_x_3d)(                                                     \
      const R *lo, const R *up, INT loc_M,                                              \
      R *x);                                                                            \
  PNFFT_EXTERN void PNX(init_x_3d_adv)(                                                 \
      const R *lo, const R *up, const R *x_max, INT loc_M,                              \
      R *x);                                                                            \
                                                                                        \
  PNFFT_EXTERN R PNX(inv_phi_hat)(                                                      \
      const PNX(plan) ths, int dim, INT k);                                             \
  PNFFT_EXTERN R PNX(phi_hat)(                                                          \
      const PNX(plan) ths, int dim, INT k);                                             \
  PNFFT_EXTERN R PNX(psi)(                                                              \
     const PNX(plan) ths, int dim, R x);                                                \
  PNFFT_EXTERN R PNX(dpsi)(                                                             \
      const PNX(plan) ths, int dim, R x);                                               \
                                                                                        \
  PNFFT_EXTERN void PNX(vpr_complex)(                                                   \
      C *data, INT N, const char *name, MPI_Comm comm);                                 \
  PNFFT_EXTERN void PNX(vpr_real)(                                                      \
      R *data, INT N, const char *name, MPI_Comm comm);                                 \
  PNFFT_EXTERN void PNX(apr_complex_3d)(                                                \
      C *data, INT *local_N, INT *local_N_start,                                        \
      const char *name, MPI_Comm comm);                                                 \
                                                                                        \
  PNFFT_EXTERN double* PNX(get_timer_trafo)(                                            \
      PNX(plan) ths);                                                                   \
  PNFFT_EXTERN double* PNX(get_timer_adj)(                                              \
      PNX(plan) ths);                                                                   \
  PNFFT_EXTERN void PNX(timer_average)(                                                 \
      double *timer);                                                                   \
  PNFFT_EXTERN double* PNX(timer_copy)(                                                 \
      const double *orig);                                                              \
  PNFFT_EXTERN double* PNX(timer_reduce_max)(                                           \
      MPI_Comm comm, double* timer);                                                    \
  PNFFT_EXTERN double* PNX(timer_add)(                                                  \
      const double* sum1, const double* sum2);                                          \
  PNFFT_EXTERN void PNX(timer_free)(                                                    \
      double* ths);                                                                     \
                                                                                        \
  PNFFT_EXTERN void PNX(reset_timer)(                                                   \
      PNX(plan) ths);                                                                   \
  PNFFT_EXTERN void PNX(print_average_timer)(                                           \
      const PNX(plan) ths, MPI_Comm comm);                                              \
  PNFFT_EXTERN void PNX(print_average_timer_adv)(                                       \
      const PNX(plan) ths, MPI_Comm comm);                                              \
  PNFFT_EXTERN void PNX(write_average_timer)(                                           \
      const PNX(plan) ths, const char *name, MPI_Comm comm);                            \
  PNFFT_EXTERN void PNX(write_average_timer_adv)(                                       \
      const PNX(plan) ths, const char *name, MPI_Comm comm);                            \
                                                                                        \
  PNFFT_EXTERN void PNX(get_args)(                                                      \
      int argc, char **argv, const char *name,                                          \
      int neededArgs, unsigned type,                                                    \
      void *parameter);
    

typedef pfft_complex pnfft_complex;
typedef pfftf_complex pnfftf_complex;
typedef pfftl_complex pnfftl_complex;

#define PNFFT_MANGLE_DOUBLE(name) PNFFT_CONCAT(pnfft_, name)
#define PNFFT_MANGLE_FLOAT(name) PNFFT_CONCAT(pnfftf_, name)
#define PNFFT_MANGLE_LONG_DOUBLE(name) PNFFT_CONCAT(pnfftl_, name)

PNFFT_DEFINE_API(
    PNFFT_MANGLE_DOUBLE, PFFT_MANGLE_DOUBLE, FFTW_MANGLE_DOUBLE,
    double, pnfft_complex, ptrdiff_t)
PNFFT_DEFINE_API(
    PNFFT_MANGLE_FLOAT, PFFT_MANGLE_FLOAT, FFTW_MANGLE_FLOAT,
    float, pnfftf_complex, ptrdiff_t)
PNFFT_DEFINE_API(
    PNFFT_MANGLE_LONG_DOUBLE, PFFT_MANGLE_LONG_DOUBLE, FFTW_MANGLE_LONG_DOUBLE,
    long double, pnfftl_complex, ptrdiff_t)




#ifndef PNFFT_PI
#define PNFFT_PI           3.14159265358979323846
#endif

/* Use double braces for simple generation of Fortran wrappers */

/* This flags is equal to PRE_PHI_HAT. We introduce it for
 * compliance with the serial NFFT interface. */
#define PNFFT_PRE_PHI_HUT      (1U<< 0)
#define PNFFT_PRE_PHI_HAT      (1U<< 0)
#define PNFFT_FG_PSI           (1U<< 1)
#define PNFFT_PRE_CONST_PSI    (1U<< 2)
#define PNFFT_PRE_LIN_PSI      (1U<< 3)
#define PNFFT_PRE_QUAD_PSI     (1U<< 4)
#define PNFFT_PRE_KUB_PSI      (1U<< 5)
#define PNFFT_PRE_PSI          (1U<< 6)
#define PNFFT_PRE_FG_PSI       ((PNFFT_PRE_PSI | PNFFT_FG_PSI))
#define PNFFT_PRE_FULL_PSI     (1U<< 7)
#define PNFFT_PRE_FULL_FG_PSI  ((PNFFT_PRE_FULL_PSI | PNFFT_FG_PSI))

#define PNFFT_MALLOC_X         (1U<< 8)
#define PNFFT_MALLOC_F_HAT     (1U<< 9)
#define PNFFT_MALLOC_F         (1U<< 10)
#define PNFFT_MALLOC_GRAD_F    (1U<< 11)

#define PNFFT_FFT_OUT_OF_PLACE (0U)
#define PNFFT_FFT_IN_PLACE     (1U<< 12)
#define PNFFT_SORT_NODES       (1U<< 13)
#define PNFFT_INTERLACED       (1U<< 14)
#define PNFFT_SHIFTED_IN       (1U<< 15)
#define PNFFT_SHIFTED_OUT      (1U<< 16)

#define PNFFT_GRAD_AD          (0U)
#define PNFFT_GRAD_IK          (1U<< 17)
#define PNFFT_GRAD_NONE        (1U<< 18) /* turn off gradient NFFT and save memory for buffers */

/* enable some optimizations for real inputs */
#define PNFFT_REAL_F           (1U<< 19)

/* default window function is Kaiser-Bessel */
#define PNFFT_WINDOW_KAISER_BESSEL  (0U)
#define PNFFT_WINDOW_GAUSSIAN       (1U<< 20)
#define PNFFT_WINDOW_BSPLINE        (1U<< 21)
#define PNFFT_WINDOW_SINC_POWER     (1U<< 22)
#define PNFFT_WINDOW_BESSEL_I0      (1U<< 23)

#define PNFFT_PRE_INTPOL_PSI ((PNFFT_PRE_CONST_PSI| PNFFT_PRE_LIN_PSI| PNFFT_PRE_QUAD_PSI| PNFFT_PRE_KUB_PSI))
#define PNFFT_PRE_ONE_PSI    ((PNFFT_PRE_INTPOL_PSI| PNFFT_PRE_FG_PSI| PNFFT_PRE_PSI| PNFFT_PRE_FULL_PSI))

#define PNFFT_FREE_X           ((PNFFT_MALLOC_X))
#define PNFFT_FREE_F_HAT       ((PNFFT_MALLOC_F_HAT))
#define PNFFT_FREE_F           ((PNFFT_MALLOC_F))
#define PNFFT_FREE_GRAD_F      ((PNFFT_MALLOC_GRAD_F))

#define PNFFT_COMPUTE_F        ((PNFFT_MALLOC_F))
#define PNFFT_COMPUTE_GRAD_F   ((PNFFT_MALLOC_GRAD_F))

#define PNFFT_INT            ((PFFT_INT))
#define PNFFT_PTRDIFF_T      ((PFFT_PTRDIFF_T))
#define PNFFT_FLOAT          ((PFFT_FLOAT))
#define PNFFT_DOUBLE         ((PFFT_DOUBLE))
#define PNFFT_UNSIGNED       ((PFFT_UNSIGNED))

/* Make length of Timer public available */
#define PNFFT_TIMER_LENGTH          (10)

#define PNFFT_TIMER_ITER            (0)
#define PNFFT_TIMER_WHOLE           (1)
#define PNFFT_TIMER_LOOP_B          (2)
#define PNFFT_TIMER_SORT_NODES      (3)
#define PNFFT_TIMER_GCELLS          (4)
#define PNFFT_TIMER_MATRIX_B        (5)
#define PNFFT_TIMER_MATRIX_F        (6)
#define PNFFT_TIMER_MATRIX_D        (7)
#define PNFFT_TIMER_SHIFT_INPUT     (8)
#define PNFFT_TIMER_SHIFT_OUTPUT    (9)





#ifdef __cplusplus
  }  /* extern "C" */
#endif /* __cplusplus */

#endif /* PNFFT_H */

