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
#include "matrix_D.h"

/* wrappers for pfft init and cleanup */
void PNX(init) (void){
  PX(init)();
}


void PNX(cleanup) (void){
  PX(cleanup)();
}


void PNX(local_size_3d)(
    const INT *N, MPI_Comm comm_cart,
    unsigned pnfft_flags,
    INT *local_N, INT *local_N_start,
    R *lower_border, R *upper_border
    )
{
  const int d=3;

  PNX(local_size_adv)(
      d, N, comm_cart, pnfft_flags,
      local_N, local_N_start, lower_border, upper_border);
}


PNX(plan) PNX(init_3d)(
    const INT *N, 
    INT local_M,
    MPI_Comm comm_cart
    )
{
  const int d=3;
  unsigned pnfft_flags, pfft_flags;

  pnfft_flags = PNFFT_MALLOC_X | PNFFT_MALLOC_F_HAT | PNFFT_MALLOC_F;
  pfft_flags = PFFT_MEASURE| PFFT_DESTROY_INPUT;

  return PNX(init_adv)(
      d, N, local_M,
      pnfft_flags, pfft_flags, comm_cart);
}


void PNX(init_nodes)(
    PNX(plan) ths, INT local_M,
    unsigned pnfft_flags, unsigned pnfft_finalize_flags
    )
{
  /* free mem and adjust pnfft_flags, compute_flags */
  PNX(free_x)(ths, pnfft_finalize_flags);
  PNX(free_f)(ths, pnfft_finalize_flags);
  PNX(free_grad_f)(ths, pnfft_finalize_flags);

  /* allocate mem and adjust pnfft_flags, compute_flags */
  ths->local_M = local_M;
  PNX(malloc_x)(ths, pnfft_flags);
  PNX(malloc_f)(ths, pnfft_flags);
  PNX(malloc_grad_f)(ths, pnfft_flags);
}

static void grad_ik_complex_input(
    PNX(plan) ths
    )
{
  /* duplicate g1 since we have to scale it several times for computing the gradient */
  ths->timer_trafo[PNFFT_TIMER_MATRIX_D] -= MPI_Wtime();
  for(INT k=0; k<ths->local_N_total; k++)
    ths->g1_buffer[k] = ths->g1[k];
  ths->timer_trafo[PNFFT_TIMER_MATRIX_D] += MPI_Wtime();

  /* calculate potentials */
  ths->timer_trafo[PNFFT_TIMER_MATRIX_F] -= MPI_Wtime();
  PNX(trafo_F)(ths);
  ths->timer_trafo[PNFFT_TIMER_MATRIX_F] += MPI_Wtime();

  ths->timer_trafo[PNFFT_TIMER_MATRIX_B] -= MPI_Wtime();
  for(INT j=0; j<ths->local_M; j++) ths->f[j] = 0;
  PNX(trafo_B_grad_ik)(ths, ths->f, 0, 1);
  ths->timer_trafo[PNFFT_TIMER_MATRIX_B] += MPI_Wtime();

  /* calculate gradient component wise */
  for(int dim =0; dim<3; dim++){
    ths->timer_trafo[PNFFT_TIMER_MATRIX_D] -= MPI_Wtime();
    PNX(scale_ik_diff_c2c)(ths->g1_buffer, ths->local_N_start, ths->local_N, dim, ths->pnfft_flags,
      ths->g1);
    ths->timer_trafo[PNFFT_TIMER_MATRIX_D] += MPI_Wtime();
    
    ths->timer_trafo[PNFFT_TIMER_MATRIX_F] -= MPI_Wtime();
    PNX(trafo_F)(ths);
    ths->timer_trafo[PNFFT_TIMER_MATRIX_F] += MPI_Wtime();

    ths->timer_trafo[PNFFT_TIMER_MATRIX_B] -= MPI_Wtime();
    for(INT j=0; j<ths->local_M; j++) ths->grad_f[3*j+dim] = 0;
    PNX(trafo_B_grad_ik)(ths, ths->grad_f, dim, 3);
    ths->timer_trafo[PNFFT_TIMER_MATRIX_B] += MPI_Wtime();
  }
}

/* parallel 3dNFFT with different window functions */
void PNX(trafo)(
    PNX(plan) ths
    )
{
  if(ths==NULL){
    PX(fprintf)(MPI_COMM_WORLD, stderr, "!!! Error: Can not execute PNFFT Plan == NULL !!!\n");
    return;
  }

  ths->timer_trafo[PNFFT_TIMER_WHOLE] -= MPI_Wtime();

  /* multiplication with matrix D */
  ths->timer_trafo[PNFFT_TIMER_MATRIX_D] -= MPI_Wtime();
  PNX(trafo_D)(ths);
  ths->timer_trafo[PNFFT_TIMER_MATRIX_D] += MPI_Wtime();
 
  if(ths->pnfft_flags & PNFFT_GRAD_IK){
    grad_ik_complex_input(ths);
  } else {
    /* multiplication with matrix F */
    ths->timer_trafo[PNFFT_TIMER_MATRIX_F] -= MPI_Wtime();
    PNX(trafo_F)(ths);
    ths->timer_trafo[PNFFT_TIMER_MATRIX_F] += MPI_Wtime();

    /* multiplication with matrix B */
    ths->timer_trafo[PNFFT_TIMER_MATRIX_B] -= MPI_Wtime();
    if(ths->pnfft_flags & PNFFT_INTERLACED)
      PNX(trafo_B_grad_ad)(ths, 1);
    else
      PNX(trafo_B_grad_ad)(ths, 0);
    ths->timer_trafo[PNFFT_TIMER_MATRIX_B] += MPI_Wtime();
  }
 
  ths->timer_trafo[PNFFT_TIMER_ITER]++;
  ths->timer_trafo[PNFFT_TIMER_WHOLE] += MPI_Wtime();
}

void PNX(adj)(
    PNX(plan) ths
    )
{
  if(ths==NULL){
    PX(fprintf)(MPI_COMM_WORLD, stderr, "!!! Error: Can not execute PNFFT Plan == NULL !!!\n");
    return;
  }

  ths->timer_adj[PNFFT_TIMER_WHOLE] -= MPI_Wtime();

  /* multiplication with matrix B^T */
  ths->timer_adj[PNFFT_TIMER_MATRIX_B] -= MPI_Wtime();
  if(ths->pnfft_flags & PNFFT_INTERLACED)
    PNX(adjoint_B)(ths, 1);
  else
    PNX(adjoint_B)(ths, 0);
  ths->timer_adj[PNFFT_TIMER_MATRIX_B] += MPI_Wtime();

  /* multiplication with matrix F^H */
  ths->timer_adj[PNFFT_TIMER_MATRIX_F] -= MPI_Wtime();
  PNX(adjoint_F)(ths);
  ths->timer_adj[PNFFT_TIMER_MATRIX_F] += MPI_Wtime();

  /* multiplication with matrix D */
  ths->timer_adj[PNFFT_TIMER_MATRIX_D] -= MPI_Wtime();
  PNX(adjoint_D)(ths);
  ths->timer_adj[PNFFT_TIMER_MATRIX_D] += MPI_Wtime();

  ths->timer_adj[PNFFT_TIMER_ITER]++;
  ths->timer_adj[PNFFT_TIMER_WHOLE] += MPI_Wtime();
}


void PNX(finalize)(
    PNX(plan) ths, unsigned pnfft_finalize_flags
    )
{
  if((pnfft_finalize_flags & PNFFT_FREE_F_HAT) && (ths->f_hat != NULL))
    PNX(free)(ths->f_hat);
  if((pnfft_finalize_flags & PNFFT_FREE_GRAD_F) && (ths->grad_f != NULL))
    PNX(free)(ths->grad_f);
  if((pnfft_finalize_flags & PNFFT_FREE_F) && (ths->f != NULL))
    PNX(free)(ths->f);
  if((pnfft_finalize_flags & PNFFT_FREE_X) && (ths->x != NULL))
    PNX(free)(ths->x);

  PX(destroy_plan)(ths->pfft_forw);
  PX(destroy_plan)(ths->pfft_back);
  PX(destroy_gcplan)(ths->gcplan);

  if(ths->g2 != ths->g1){
    if(ths->g2 != NULL)
      PNX(free)(ths->g2);
  }
  if(ths->g1 != NULL)
    PNX(free)(ths->g1);
  if(ths->g1_buffer != NULL)
    PNX(free)(ths->g1_buffer);

  if(ths->intpol_tables_psi != NULL){
    for(int t=0;t<ths->d; t++)
      if(ths->intpol_tables_psi[t] != NULL)
        PNX(free)(ths->intpol_tables_psi[t]);
    PNX(free)(ths->intpol_tables_psi);
  }

  if(ths->intpol_tables_dpsi != NULL){
    for(int t=0;t<ths->d; t++)
      if(ths->intpol_tables_dpsi[t] != NULL)
        PNX(free)(ths->intpol_tables_dpsi[t]);
    PNX(free)(ths->intpol_tables_dpsi);
  }

  PNX(free)(ths->x_max);
  PNX(free)(ths->sigma);
  PNX(free)(ths->n);
  PNX(free)(ths->no);
  PNX(free)(ths->N);
  PNX(free)(ths->local_N);
  PNX(free)(ths->local_N_start);
  PNX(free)(ths->local_no);
  PNX(free)(ths->local_no_start);

  MPI_Comm_free(&(ths->comm_cart));

  /* finalize window specific parameters */
  if(ths->b != NULL)
    PNX(free)(ths->b);
  if(ths->exp_const != NULL)
    PNX(free)(ths->exp_const);
  if(ths->spline_coeffs != NULL)
    PNX(free)(ths->spline_coeffs);
  if(ths->pre_inv_phi_hat_trafo != NULL)
    PNX(free)(ths->pre_inv_phi_hat_trafo);
  if(ths->pre_inv_phi_hat_adj != NULL)
    PNX(free)(ths->pre_inv_phi_hat_adj);

  if(ths->pre_psi != NULL)     PNX(free)(ths->pre_psi);
  if(ths->pre_dpsi != NULL)    PNX(free)(ths->pre_dpsi);
  if(ths->pre_psi_il != NULL)  PNX(free)(ths->pre_psi_il);
  if(ths->pre_dpsi_il != NULL) PNX(free)(ths->pre_dpsi_il);

  /* free mem of struct */
  PNX(rmplan)(ths);
}

void PNX(set_f_hat)(
    C *f_hat, PNX(plan) ths
    )
{
  ths->f_hat = f_hat;
}

C* PNX(get_f_hat)(
    const PNX(plan) ths
    )
{
  return ths->f_hat;
}

void PNX(set_f)(
    C *f, PNX(plan) ths
    )
{
  ths->f = f;
}

C* PNX(get_f)(
    const PNX(plan) ths
    )
{
  return ths->f;
}

void PNX(set_grad_f)(
    C* grad_f, PNX(plan) ths
    )
{
  ths->grad_f = grad_f;
}

C* PNX(get_grad_f)(
    const PNX(plan) ths
    )
{
  return ths->grad_f;
}

void PNX(set_x)(
    R *x, PNX(plan) ths
    )
{
  ths->x = x;
}

R* PNX(get_x)(
    const PNX(plan) ths
    )
{
  return ths->x;
}


/* getters for PNFFT internal parameters
 * No setters are implemented for these parameters.
 * Use finalize and init_guru instead. */
int PNX(get_d)(
    const PNX(plan) ths
    )
{
  return ths->d;
}

int PNX(get_m)(
    const PNX(plan) ths
    )
{
  return ths->m;
}

void PNX(get_x_max)(
    const PNX(plan) ths,
    R *x_max
    )
{
  for(int t=0; t<ths->d; t++)
    x_max[t] = ths->x_max[t];
}

void PNX(get_N)(
    const PNX(plan) ths,
    INT *N
    )
{
  for(int t=0; t<ths->d; t++)
    N[t] = ths->N[t];
}

void PNX(get_n)(
    const PNX(plan) ths,
    INT *n
    )
{
  for(int t=0; t<ths->d; t++)
    n[t] = ths->n[t];
}

unsigned PNX(get_pnfft_flags)(
    const PNX(plan) ths
    )
{
  return ths->pnfft_flags;
}

unsigned PNX(get_pfft_flags)(
    const PNX(plan) ths
    )
{
  return ths->pfft_opt_flags;
}

void PNX(init_f_hat_3d)(
    const INT *N, const INT *local_N, const INT *local_N_start,
    unsigned pnfft_flags,
    C *data
    )
{
  INT local_Nt[3], local_Nt_start[3];
  int shift = (pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? 1 : 0;

  for(int t=0; t<3; t++){
    local_Nt[t] = local_N[(t + shift) % 3];
    local_Nt_start[t] = local_N_start[(t + shift) % 3];
  }

  PX(init_input_complex_3d)(N, local_Nt, local_Nt_start,
      data);
}

void PNX(init_f)(
    INT local_M,
    C *data
    )
{
  for (INT j=0; j<local_M; j++){
    R real = 100.0 * (R) rand() / RAND_MAX;
    R imag = 100.0 * (R) rand() / RAND_MAX;
    data[j] = real + imag * I;
  }
}

void PNX(init_x_3d)(
    const R *lo, const R *up, INT loc_M,
    R *x
    )
{
  R x_max[3] = {0.5,0.5,0.5};

  PNX(init_x_3d_adv)(lo, up, x_max, loc_M,
      x);
}


static void print_complex_vector(
    R *data, INT N
    )
{
  for(INT l=0; l<N; l++){
    if(l%4 == 0)
      printf("\n%4td.", l/4);
    printf(" %.2e+%.2ei,", data[2*l], data[2*l+1]);
  }
  printf("\n");
}


void PNX(vpr_complex)(
    C *data, INT N,
    const char *name, MPI_Comm comm
    )
{
  int size, myrank;

  if(N < 1)
    return;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &myrank);
  
  fflush(stdout);
  MPI_Barrier(comm);
  for(int t=0; t<size; t++){
    if(t==myrank){
      printf("\nRank %d, %s", myrank, name);
      print_complex_vector((R*) data, N);
      fflush(stdout);
    }
    MPI_Barrier(comm);
  }
}


void PNX(vpr_real)(
    R *data, INT N,
    const char *name, MPI_Comm comm
    )
{
  int size, myrank;

  if(N < 1)
    return;
  
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &myrank);
  
  fflush(stdout);
  MPI_Barrier(comm);
  for(int t=0; t<size; t++){
    if(t==myrank){
      printf("\nRank %d, %s", myrank, name);
      for(INT l=0; l<N; l++){
        if(l%8 == 0)
          printf("\n%4td.", l/8);
        printf(" %e,", data[l]);
      }
      printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(comm);
  }
}


void PNX(apr_complex_3d)(
     C *data, INT *local_N, INT *local_N_start, unsigned pnfft_flags,
     const char *name, MPI_Comm comm
     )
{
  INT local_Nt[3], local_Nt_start[3];
  int shift = (pnfft_flags & PNFFT_TRANSPOSED_F_HAT) ? 1 : 0;

  for(int t=0; t<3; t++){
    local_Nt[t] = local_N[(t + shift) % 3];
    local_Nt_start[t] = local_N_start[(t + shift) % 3];
  }

  PX(apr_complex_3d)(
     data, local_Nt, local_Nt_start, name, comm);
}
