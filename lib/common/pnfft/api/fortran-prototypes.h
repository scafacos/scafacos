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

PNFFT_VOIDFUNC FORT(create_procmesh_2d, CREATE_PROCMESH_2D)(
    int *ierror,
    MPI_Fint *comm, int *np2, int *np3,
    MPI_Fint *comm_cart_2d
    );

PNFFT_VOIDFUNC FORT(create_procmesh, CREATE_PROCMESH)(
    int *ierror,
    int *rnk, MPI_Fint *comm, int *np,
    MPI_Fint *comm_cart
    );

PNFFT_VOIDFUNC FORT(local_size_3d, LOCAL_SIZE_3D)(
    const INT *N, MPI_Fint *comm_cart,
    INT *local_N, INT *local_N_start,
    R *lower_border, R *upper_border
    );

PNFFT_VOIDFUNC FORT(local_size_adv, LOCAL_SIZE_ADV)(
    int *d, const INT *N, MPI_Fint *comm_cart,
    INT *local_N, INT *local_N_start,
    R *lower_border, R *upper_border
    );

PNFFT_VOIDFUNC FORT(local_size_guru, LOCAL_SIZE_GURU)(
    int *d, const INT *N, const INT *n, const R *x_max, int *m,
    MPI_Fint *comm_cart,
    INT *local_N, INT *local_N_start,
    R *lower_border, R *upper_border
    );

PNFFT_VOIDFUNC FORT(init_3d, INIT_3D)(
    PNX(plan) *p,
    const INT* N, INT *local_M, MPI_Fint *comm_cart
    );

PNFFT_VOIDFUNC FORT(init_adv, INIT_ADV)(
    PNX(plan) *p,
    int *d, const INT* N, INT *local_M,
    int *pnfft_flags, int *fftw_flags, 
    MPI_Fint *comm_cart
    );

PNFFT_VOIDFUNC FORT(init_guru, INIT_GURU)(
    PNX(plan) *p,
    int *d, const INT* N, const INT *n, const R *x_max,
    INT *local_M, int *m,
    int *pnfft_flags, int *fftw_flags, 
    MPI_Fint *comm_cart
    );

PNFFT_VOIDFUNC FORT(init_nodes, INIT_NODES)(
    PNX(plan) *ths, INT *local_M,
    int *pnfft_flags, int *pnfft_finalize_flags
    );

PNFFT_VOIDFUNC FORT(precompute_psi, PRECOMPUTE_PSI)(
    PNX(plan) *ths
    );

PNFFT_VOIDFUNC FORT(set_f_hat, SET_F_HAT)(
    C *f_hat, PNX(plan) *ths
    );
PNFFT_VOIDFUNC FORT(set_f, SET_F)(
    C *f, PNX(plan) *ths
    );
PNFFT_VOIDFUNC FORT(set_grad_f, SET_GRAD_F)(
    C *grad_f, PNX(plan) *ths
    );
PNFFT_VOIDFUNC FORT(set_x, SET_X)(
    R *x, PNX(plan) *ths
    );
                                     
PNFFT_VOIDFUNC FORT(get_f_hat, GET_F_HAT)(
    C **f_hat, PNX(plan) * const ths
    );
PNFFT_VOIDFUNC FORT(get_f, GET_F)(
    C **f, PNX(plan) * const ths
    );
PNFFT_VOIDFUNC FORT(get_grad_f, GET_GRAD_F)(
    C **grad_f, PNX(plan) * const ths
    );
PNFFT_VOIDFUNC FORT(get_x, GET_X)(
    R **x, PNX(plan) * const ths
    );

PNFFT_VOIDFUNC FORT(get_d, GET_D)(
    int *d, PNX(plan) * const ths
    );
PNFFT_VOIDFUNC FORT(get_m, GET_M)(
    int *m, PNX(plan) * const ths
    );
PNFFT_VOIDFUNC FORT(get_x_max, GET_X_MAX)(
    PNX(plan) * const ths, R *x_max
    );
PNFFT_VOIDFUNC FORT(get_N, GET_N)(
    PNX(plan) * const ths, INT *N
    );
PNFFT_VOIDFUNC FORT(get_nos, GET_NOS)(
    PNX(plan) * const ths, INT *n
    );
PNFFT_VOIDFUNC FORT(get_pnfft_flags, GET_PNFFT_FLAGS)(
    int *pnfft_flags, PNX(plan) * const ths
    );
PNFFT_VOIDFUNC FORT(get_pfft_flags, GET_PFFT_FLAGS)(
    int *pfft_flags, PNX(plan) * const ths
    );

PNFFT_VOIDFUNC FORT(finalize, FINALIZE)(
    PNX(plan) * const ths, int *pnfft_finalize_flags
    );

PNFFT_VOIDFUNC FORT(trafo, TRAFO)(
    PNX(plan) * const ths
    );
PNFFT_VOIDFUNC FORT(adj, ADJ)(
    PNX(plan) * const ths
    );


PNFFT_VOIDFUNC FORT(init, INIT)(void);

PNFFT_VOIDFUNC FORT(cleanup, CLEANUP)(void);

PNFFT_VOIDFUNC FORT(malloc, MALLOC)(
    void **ptr, size_t n
    );
PNFFT_VOIDFUNC FORT(alloc_real, ALLOC_REAL)(
    R **ptr, size_t n
    );
PNFFT_VOIDFUNC FORT(alloc_complex, ALLOC_COMPLEX)(
    C **ptr, size_t n
    );
PNFFT_VOIDFUNC FORT(free, FREE)(
    void *ptr
    );

PNFFT_VOIDFUNC FORT(init_f_hat_3d, INIT_F_HAT_3D)(
      const INT *N, const INT *local_N, const INT *local_N_start,
      C *data
      );
PNFFT_VOIDFUNC FORT(init_x_3d, INIT_X_3D)(
      const R *lo, const R *up, const INT *local_M,
      R *x
      );
PNFFT_VOIDFUNC FORT(init_x_3d_adv, INIT_X_3D_ADV)(
      const R *lo, const R *up, const R *x_max, const INT *local_M,
      R *x
      );

PNFFT_VOIDFUNC FORT(inv_phi_hat, INV_PHI_HAT)(
    R *res,
    const PNX(plan) *p, int *dim, INT *k
    );
PNFFT_VOIDFUNC FORT(phi_hat, PHI_HAT)(
    R *res,
    const PNX(plan) *p, int *dim, INT *k
    );
PNFFT_VOIDFUNC FORT(psi, PSI)(
    R *res,
    const PNX(plan) *p, int *dim, R *x
    );
PNFFT_VOIDFUNC FORT(dpsi, DPSI)(
    R *res,
    const PNX(plan) *p, int *dim, R *x
    );

PNFFT_VOIDFUNC FORT(vpr_complex, VPR_COMPLEX)(
    C *data, INT *N, const char *name, MPI_Fint *comm_cart
    );
PNFFT_VOIDFUNC FORT(vpr_real, VPR_REAL)(
    R *data, INT *N, const char *name, MPI_Fint *comm_cart
    );
PNFFT_VOIDFUNC FORT(apr_complex_3d, APR_COMPLEX_3D)(
    C *data, const INT* local_N, const INT *local_N_start,
    const char *name, MPI_Fint *comm_cart
    );


PNFFT_VOIDFUNC FORT(get_timer_trafo, GET_TIMER_TRAFO)(
    double **timer,
    const PNX(plan) *p
    );
PNFFT_VOIDFUNC FORT(get_timer_adj, GET_TIMER_ADJ)(
    double **timer,
    const PNX(plan) *p
    );
PNFFT_VOIDFUNC FORT(timer_average, TIMER_AVERAGE)(
    double **timer
    );
PNFFT_VOIDFUNC FORT(timer_copy, TIMER_COPY)(
    double **timer,
    const double *orig
    );
PNFFT_VOIDFUNC FORT(timer_reduce_max, TIMER_REDUCE_MAX)(
    MPI_Fint *comm, double *timer
    );
PNFFT_VOIDFUNC FORT(timer_add, TIMER_ADD)(
    double **timer,
    const double *sum1, const double *sum2
    );
PNFFT_VOIDFUNC FORT(timer_free, TIMER_FREE)(
    double **timer
    );

PNFFT_VOIDFUNC FORT(reset_timer, RESET_TIMER)(
    PNX(plan) *ths
    );
PNFFT_VOIDFUNC FORT(print_average_timer, PRINT_AVERAGE_TIMER)(
    const PNX(plan) *ths, MPI_Fint *comm
    );
PNFFT_VOIDFUNC FORT(print_average_timer_adv, PRINT_AVERAGE_TIMER_ADV)(
    const PNX(plan) *ths, MPI_Fint *comm
    );
PNFFT_VOIDFUNC FORT(write_average_timer, PRINT_AVERAGE_TIMER)(
    const PNX(plan) *ths, const char *name, MPI_Fint *comm
    );
PNFFT_VOIDFUNC FORT(write_average_timer_adv, PRINT_AVERAGE_TIMER_ADV)(
    const PNX(plan) *ths, const char *name, MPI_Fint *comm
    );
