/*
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
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

/* Functions in the PNFFT Fortran API, mangled according to the
   FORT(...) macro.  This file is designed to be #included by
   fortran_api.c, possibly multiple times in order to support multiple
   compiler manglings (via redefinition of FORT). */

/* Generate prototypes from fortran-wrappers.h using the vim search/replace
 * command :%s#)\n{[A-Za-z0-9_,\n ();*=&/\[\]-]*}#);#gc
 */

/* include prototypes of all fortran wrappers to avoid pedantic gcc warnings  */
#include <fortran-prototypes.h>

PNFFT_VOIDFUNC FORT(create_procmesh_2d, CREATE_PROCMESH_2D)(
    int *ierror,
    MPI_Fint *comm, int *np2, int *np3,
    MPI_Fint *comm_cart_2d
    )
{
  MPI_Comm comm_cart_2d_clike;

  *ierror = PX(create_procmesh_2d)(MPI_Comm_f2c(*comm), *np3, *np2, &comm_cart_2d_clike);
  
  /* C to Fortran type conversion */
  *comm_cart_2d = MPI_Comm_c2f(comm_cart_2d_clike);
}

PNFFT_VOIDFUNC FORT(create_procmesh, CREATE_PROCMESH)(
    int *ierror,
    int *rnk, MPI_Fint *comm, int *np,
    MPI_Fint *comm_cart
    )
{
  MPI_Comm comm_cart_clike;
  int *np_rev = malloc_and_revert_int(*rnk, np);

  *ierror = PX(create_procmesh)(*rnk, MPI_Comm_f2c(*comm), np_rev, &comm_cart_clike);
  
  /* C to Fortran type conversion */
  *comm_cart = MPI_Comm_c2f(comm_cart_clike);

  free(np_rev);
}

PNFFT_VOIDFUNC FORT(local_size_3d, LOCAL_SIZE_3D)(
    const INT *N, MPI_Fint *comm_cart,
    INT *local_N, INT *local_N_start,
    R *lower_border, R *upper_border
    )
{
  INT *N_rev = malloc_and_revert_INT(3, N);
  INT *local_N_rev = PNX(malloc_INT)(3);
  INT *local_N_start_rev = PNX(malloc_INT)(3);
  R *lower_border_rev = PNX(malloc_R)(3);
  R *upper_border_rev = PNX(malloc_R)(3);

  PNX(local_size_3d)(
      N_rev, MPI_Comm_f2c(*comm_cart),
      local_N_rev, local_N_start_rev,
      lower_border_rev, upper_border_rev);

  revert_INT(3, local_N_rev, local_N);
  revert_and_add_ones_INT(3, local_N_start_rev, local_N_start);
  revert_R(3, lower_border_rev, lower_border);
  revert_R(3, upper_border_rev, upper_border);

  PNX(free)(local_N_rev); PNX(free)(local_N_start_rev);
  PNX(free)(lower_border_rev); PNX(free)(upper_border_rev);
}

PNFFT_VOIDFUNC FORT(local_size_adv, LOCAL_SIZE_ADV)(
    int *d, const INT *N, MPI_Fint *comm_cart,
    INT *local_N, INT *local_N_start,
    R *lower_border, R *upper_border
    )
{
  INT *N_rev = malloc_and_revert_INT(*d, N);
  INT *local_N_rev = PNX(malloc_INT)(*d);
  INT *local_N_start_rev = PNX(malloc_INT)(*d);
  R *lower_border_rev = PNX(malloc_R)(*d);
  R *upper_border_rev = PNX(malloc_R)(*d);

  PNX(local_size_adv)(
      *d, N_rev, MPI_Comm_f2c(*comm_cart),
      local_N_rev, local_N_start_rev,
      lower_border_rev, upper_border_rev);

  revert_INT(*d, local_N_rev, local_N);
  revert_and_add_ones_INT(*d, local_N_start_rev, local_N_start);
  revert_R(*d, lower_border_rev, lower_border);
  revert_R(*d, upper_border_rev, upper_border);

  PNX(free)(local_N_rev); PNX(free)(local_N_start_rev);
  PNX(free)(lower_border_rev); PNX(free)(upper_border_rev);
}

PNFFT_VOIDFUNC FORT(local_size_guru, LOCAL_SIZE_GURU)(
    int *d, const INT *N, const INT *n, const R *x_max, int *m,
    MPI_Fint *comm_cart,
    INT *local_N, INT *local_N_start,
    R *lower_border, R *upper_border
    )
{
  INT *N_rev = malloc_and_revert_INT(*d, N);
  INT *n_rev = malloc_and_revert_INT(*d, n);
  INT *local_N_rev = PNX(malloc_INT)(*d);
  INT *local_N_start_rev = PNX(malloc_INT)(*d);
  R *lo_rev = PNX(malloc_R)(*d);
  R *up_rev = PNX(malloc_R)(*d);
  R *x_max_rev = malloc_and_revert_R(*d, x_max);

  PNX(local_size_guru)(
      *d, N_rev, n_rev, x_max_rev, *m, MPI_Comm_f2c(*comm_cart),
      local_N_rev, local_N_start_rev,
      lo_rev, up_rev);

  revert_INT(*d, local_N_rev, local_N);
  revert_and_add_ones_INT(*d, local_N_start_rev, local_N_start);
  revert_R(*d, lo_rev, lower_border);
  revert_R(*d, up_rev, upper_border);

  PNX(free)(local_N_rev); PNX(free)(local_N_start_rev);
  PNX(free)(lo_rev); PNX(free)(up_rev); PNX(free)(x_max_rev);
}

PNFFT_VOIDFUNC FORT(init_3d, INIT_3D)(
    PNX(plan) *p,
    const INT* N, INT *local_M, MPI_Fint *comm_cart
    )
{
  INT *N_rev = malloc_and_revert_INT(3, N);
  *p = PNX(init_3d)(N_rev, *local_M, MPI_Comm_f2c(*comm_cart));
  PNX(free)(N_rev);
}

PNFFT_VOIDFUNC FORT(init_adv, INIT_ADV)(
    PNX(plan) *p,
    int *d, const INT* N, INT *local_M,
    int *pnfft_flags, int *fftw_flags, 
    MPI_Fint *comm_cart
    )
{
  INT *N_rev = malloc_and_revert_INT(*d, N);
  *p = PNX(init_adv)(*d, N_rev, *local_M, (unsigned) *pnfft_flags, (unsigned) *fftw_flags, MPI_Comm_f2c(*comm_cart));
  PNX(free)(N_rev);
}

PNFFT_VOIDFUNC FORT(init_guru, INIT_GURU)(
    PNX(plan) *p,
    int *d, const INT* N, const INT *n, const R *x_max,
    INT *local_M, int *m,
    int *pnfft_flags, int *fftw_flags, 
    MPI_Fint *comm_cart
    )
{
  INT *N_rev = malloc_and_revert_INT(*d, N);
  INT *n_rev = malloc_and_revert_INT(*d, n);
  R *x_max_rev = malloc_and_revert_R(*d, x_max);

  *p = PNX(init_guru)(*d, N_rev, n_rev, x_max_rev, *local_M, *m,
      (unsigned) *pnfft_flags, (unsigned) *fftw_flags, MPI_Comm_f2c(*comm_cart));
  
  PNX(free)(N_rev); PNX(free)(n_rev);
  PNX(free)(x_max_rev);
}

PNFFT_VOIDFUNC FORT(init_nodes, INIT_NODES)(
    PNX(plan) *ths, INT *local_M,
    int *pnfft_flags, int *pnfft_finalize_flags
    )
{
  PNX(init_nodes)(
      *ths, *local_M, (unsigned) *pnfft_flags, (unsigned) *pnfft_finalize_flags);
}

PNFFT_VOIDFUNC FORT(precompute_psi, PRECOMPUTE_PSI)(
    PNX(plan) *ths
    )
{
  PNX(precompute_psi)(*ths);
}

PNFFT_VOIDFUNC FORT(set_f_hat, SET_F_HAT)(
    C *f_hat, PNX(plan) *ths
    )
{
  PNX(set_f_hat)(f_hat, *ths);
}
PNFFT_VOIDFUNC FORT(set_f, SET_F)(
    C *f, PNX(plan) *ths
    )
{
  PNX(set_f)(f, *ths);
}
PNFFT_VOIDFUNC FORT(set_grad_f, SET_GRAD_F)(
    C *grad_f, PNX(plan) *ths
    )
{
  PNX(set_grad_f)(grad_f, *ths);
}
PNFFT_VOIDFUNC FORT(set_x, SET_X)(
    R *x, PNX(plan) *ths
    )
{
  PNX(set_x)(x, *ths);
}
                                     
PNFFT_VOIDFUNC FORT(get_f_hat, GET_F_HAT)(
    C **f_hat, PNX(plan) * const ths
    )
{
  *f_hat = PNX(get_f_hat)(*ths);
}
PNFFT_VOIDFUNC FORT(get_f, GET_F)(
    C **f, PNX(plan) * const ths
    )
{
  *f = PNX(get_f)(*ths);
}
PNFFT_VOIDFUNC FORT(get_grad_f, GET_GRAD_F)(
    C **grad_f, PNX(plan) * const ths
    )
{
  *grad_f = PNX(get_grad_f)(*ths);
}
PNFFT_VOIDFUNC FORT(get_x, GET_X)(
    R **x, PNX(plan) * const ths
    )
{
  *x = PNX(get_x)(*ths);
}

PNFFT_VOIDFUNC FORT(get_d, GET_D)(
    int *d, PNX(plan) * const ths
    )
{
  *d = PNX(get_d)(*ths);
}
PNFFT_VOIDFUNC FORT(get_m, GET_M)(
    int *m, PNX(plan) * const ths
    )
{
  *m = PNX(get_m)(*ths);
}
PNFFT_VOIDFUNC FORT(get_x_max, GET_X_MAX)(
    R **x_max, PNX(plan) * const ths
    )
{
  *x_max = PNX(get_x_max)(*ths);
}
PNFFT_VOIDFUNC FORT(get_N, GET_N)(
    INT **N, PNX(plan) * const ths
    )
{
  *N = PNX(get_N)(*ths);
}
PNFFT_VOIDFUNC FORT(get_nos, GET_NOS)(
    INT **n, PNX(plan) * const ths
    )
{
  *n = PNX(get_n)(*ths);
}
PNFFT_VOIDFUNC FORT(get_pnfft_flags, GET_PNFFT_FLAGS)(
    int *pnfft_flags, PNX(plan) * const ths
    )
{
  *pnfft_flags = PNX(get_pnfft_flags)(*ths);
}
PNFFT_VOIDFUNC FORT(get_pfft_flags, GET_PFFT_FLAGS)(
    int *pfft_flags, PNX(plan) * const ths
    )
{
  *pfft_flags = PNX(get_pfft_flags)(*ths);
}

PNFFT_VOIDFUNC FORT(finalize, FINALIZE)(
    PNX(plan) * const ths, int *pnfft_finalize_flags
    )
{
  PNX(finalize)(*ths, *pnfft_finalize_flags);
}

PNFFT_VOIDFUNC FORT(trafo, TRAFO)(
    PNX(plan) * const ths
    )
{
  PNX(trafo)(*ths);
}
PNFFT_VOIDFUNC FORT(adj, ADJ)(
    PNX(plan) * const ths
    )
{
  PNX(adj)(*ths);
}


PNFFT_VOIDFUNC FORT(init, INIT)(void)
{
  PNX(init)();
}

PNFFT_VOIDFUNC FORT(cleanup, CLEANUP)(void)
{
  PNX(cleanup)();
}

PNFFT_VOIDFUNC FORT(malloc, MALLOC)(
    void **ptr, size_t n
    )
{
  *ptr = PNX(malloc)(n);
}
PNFFT_VOIDFUNC FORT(alloc_real, ALLOC_REAL)(
    R **ptr, size_t n
    )
{
  *ptr = PNX(alloc_real)(n);
}
PNFFT_VOIDFUNC FORT(alloc_complex, ALLOC_COMPLEX)(
    C **ptr, size_t n
    )
{
  *ptr = PNX(alloc_complex)(n);
}
PNFFT_VOIDFUNC FORT(free, FREE)(
    void *ptr
    )
{
  PNX(free)(ptr);
}

PNFFT_VOIDFUNC FORT(init_f_hat_3d, INIT_F_HAT_3D)(
      const INT *N, const INT *local_N, const INT *local_N_start,
      C *data
      )
{
  INT *N_rev = malloc_and_revert_INT(3, N);
  INT *local_N_rev = malloc_and_revert_INT(3, local_N);
  INT *local_N_start_rev = malloc_and_revert_INT(3, local_N_start);

  PNX(init_f_hat_3d)(
      N_rev, local_N_rev, local_N_start_rev,
      data);

  PNX(free)(N_rev); PNX(free)(local_N_rev); PNX(free)(local_N_start_rev);
}
PNFFT_VOIDFUNC FORT(init_x_3d, INIT_X_3D)(
      const R *lo, const R *up, const INT *local_M,
      R *x
      )
{
  R *lo_rev = malloc_and_revert_R(3, lo);
  R *up_rev = malloc_and_revert_R(3, up);

  PNX(init_x_3d)(
      lo_rev, up_rev, *local_M,
      x);

  PNX(free)(lo_rev); PNX(free)(up_rev);
}
PNFFT_VOIDFUNC FORT(init_x_3d_adv, INIT_X_3D_ADV)(
      const R *lo, const R *up, const R *x_max, const INT *local_M,
      R *x
      )
{
  R *lo_rev = malloc_and_revert_R(3, lo);
  R *up_rev = malloc_and_revert_R(3, up);
  R *x_max_rev = malloc_and_revert_R(3, x_max);

  PNX(init_x_3d_adv)(
      lo_rev, up_rev, x_max, *local_M,
      x);

  PNX(free)(lo_rev); PNX(free)(up_rev); PNX(free)(x_max_rev);
}

PNFFT_VOIDFUNC FORT(inv_phi_hat, INV_PHI_HAT)(
    R *res,
    const PNX(plan) *p, int *dim, INT *k
    )
{
  *res = PNX(inv_phi_hat)(*p, *dim, *k);
}
PNFFT_VOIDFUNC FORT(phi_hat, PHI_HAT)(
    R *res,
    const PNX(plan) *p, int *dim, INT *k
    )
{
  *res = PNX(phi_hat)(*p, *dim, *k);
}
PNFFT_VOIDFUNC FORT(psi, PSI)(
    R *res,
    const PNX(plan) *p, int *dim, R *x
    )
{
  *res = PNX(psi)(*p, *dim, *x);
}
PNFFT_VOIDFUNC FORT(dpsi, DPSI)(
    R *res,
    const PNX(plan) *p, int *dim, R *x
    )
{
  *res = PNX(dpsi)(*p, *dim, *x);
}

PNFFT_VOIDFUNC FORT(vpr_complex, VPR_COMPLEX)(
    C *data, INT *N, const char *name, MPI_Fint *comm_cart
    )
{
  PNX(vpr_complex)(data, *N, name, MPI_Comm_f2c(*comm_cart));
}
PNFFT_VOIDFUNC FORT(vpr_real, VPR_REAL)(
    R *data, INT *N, const char *name, MPI_Fint *comm_cart
    )
{
  PNX(vpr_real)(data, *N, name, MPI_Comm_f2c(*comm_cart));
}
PNFFT_VOIDFUNC FORT(apr_complex_3d, APR_COMPLEX_3D)(
    C *data, const INT* local_N, const INT *local_N_start,
    const char *name, MPI_Fint *comm_cart
    )
{
  INT *local_N_rev = malloc_and_revert_INT(3, local_N);
  INT *local_N_start_rev = malloc_and_revert_INT(3, local_N_start);
  PNX(apr_complex_3d)(data, local_N_rev, local_N_start_rev, name, MPI_Comm_f2c(*comm_cart));
  PNX(free)(local_N_rev); PNX(free)(local_N_start_rev);
}


PNFFT_VOIDFUNC FORT(get_timer_trafo, GET_TIMER_TRAFO)(
    double **timer,
    const PNX(plan) *p
    )
{
  *timer = PNX(get_timer_trafo)(*p);
}
PNFFT_VOIDFUNC FORT(get_timer_adj, GET_TIMER_ADJ)(
    double **timer,
    const PNX(plan) *p
    )
{
  *timer = PNX(get_timer_adj)(*p);
}
PNFFT_VOIDFUNC FORT(timer_average, TIMER_AVERAGE)(
    double **timer
    )
{
  PNX(timer_average)(*timer);
}
PNFFT_VOIDFUNC FORT(timer_copy, TIMER_COPY)(
    double **timer,
    const double *orig
    )
{
  *timer = PNX(timer_copy)(orig);
}
PNFFT_VOIDFUNC FORT(timer_reduce_max, TIMER_REDUCE_MAX)(
    MPI_Fint *comm, double *timer
    )
{
  PNX(timer_reduce_max)(MPI_Comm_f2c(*comm), timer);
}
PNFFT_VOIDFUNC FORT(timer_add, TIMER_ADD)(
    double **timer,
    const double *sum1, const double *sum2
    )
{
  *timer = PNX(timer_add)(sum1, sum2);
}
PNFFT_VOIDFUNC FORT(timer_free, TIMER_FREE)(
    double **timer
    )
{
  PNX(timer_free)(*timer);
}

PNFFT_VOIDFUNC FORT(reset_timer, RESET_TIMER)(
    PNX(plan) *ths
    )
{
  PNX(reset_timer)(*ths);
}
PNFFT_VOIDFUNC FORT(print_average_timer, PRINT_AVERAGE_TIMER)(
    const PNX(plan) *ths, MPI_Fint *comm
    )
{
  PNX(print_average_timer)(*ths, MPI_Comm_f2c(*comm));
}
PNFFT_VOIDFUNC FORT(print_average_timer_adv, PRINT_AVERAGE_TIMER_ADV)(
    const PNX(plan) *ths, MPI_Fint *comm
    )
{
  PNX(print_average_timer_adv)(*ths, MPI_Comm_f2c(*comm));
}
PNFFT_VOIDFUNC FORT(write_average_timer, PRINT_AVERAGE_TIMER)(
    const PNX(plan) *ths, const char *name, MPI_Fint *comm
    )
{
  PNX(write_average_timer)(*ths, name, MPI_Comm_f2c(*comm));
}
PNFFT_VOIDFUNC FORT(write_average_timer_adv, PRINT_AVERAGE_TIMER_ADV)(
    const PNX(plan) *ths, const char *name, MPI_Fint *comm
    )
{
  PNX(write_average_timer_adv)(*ths, name, MPI_Comm_f2c(*comm));
}
