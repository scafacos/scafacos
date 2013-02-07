/*
 * Copyright (C) 2012,2013 Michael Pippig
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
#include <math.h>

#include "kernels.h"
#include "regularization.h"
#include "cg_approx_1d.h"

/* Optimize the approximation in the intervall [eps_I, 0.5-eps_B] only.
 * The result may not be smooth to a constant continuation.
 * This is for comparison. It shows how much we loose because of the far field regularization. */
#define ONLY_FARFIELD 0

/* If this flag is set, the approximation error is calculated for all x in [eps_I, 0.5-eps_B],
 * otherwise for all x in [eps_I, 0.5]. Especially the part [0.5-eps_B,0.5] can be very expensive
 * for high smoothness p. */
#define NO_FULL_ERROR_CALC 1

/* wrap the kernel function to avoid the nasty param variable */
static double kernel(
    double r, int p, double eps_I, double eps_B);

/* wrapper for generating different node distributions */
static int init_nodes_1d(
    int M, double eps_I, double eps_B, double *x);

/* wrapper for counting the node of a distribution */
static int count_nodes_1d(
    int M, double eps_I, double eps_B);
/* count the nodes of an equispaced distribution */
static int count_equispaced_nodes_1d(
    int M, double eps_I, double eps_B);

/* generate different node distributions */
static int init_chebyshev_nodes_1d(
    int M, double eps_I, double eps_B,
    double* x);

/* evaluate kernel at given nodes y = K(x) */
static void init_samples_1d(
    int num_nodes, double* x, double* y,
    int p, double eps_I, double eps_B);

/* calculate error of approximation, max{f - K(x)} */
static double calculate_error_1d(
    int num_nodes, double* x, double* f,
    int p, double eps_I, double eps_B,
    int verbose);




/* This wrapper is used to hide the parameters of the CG algorithm.
 * Here N is the degree of the cosine polynomial. Therefore, it is only one half
 * of the exponential degree. */
double optimize_cos_coeff_via_cg(
    int N, double eps_I, double eps_B, int p,
    double *f_hat_opt
    )
{
  int M = 400;
  int M_check = 10001; 
  int m = 8;
  int iter = 500;
  int verbose = 0;

  return optimize_cos_coeff_via_cg_guru(N, M, M_check, eps_I, eps_B, p, m, iter, verbose,
      f_hat_opt);
}


double optimize_cos_coeff_via_cg_guru(
    int N, int M, double eps_I, double eps_B, int p, int iter, int verbose,
    double *f_hat_opt
    )
{
  int k,l,n;                        /**< index for nodes, freqencies, iter */
//   nfct_plan plan, plan_check;       /**< plan for the nfcts                */
//   solver_plan_double ip;            /**< plan for the inverse nfct         */
 
  double error;
  int num_nodes;
  
  /** initialize one dimensional nfct plan for the CG method */
  num_nodes = count_nodes_1d(M, eps_I, eps_B);
  nfct_init_guru(&plan, 1, &N, num_nodes, &n, m,
                 PRE_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|
                 FFTW_INIT| FFT_OUT_OF_PLACE,
                 FFTW_ESTIMATE| FFTW_DESTROY_INPUT);
  init_nodes_1d(M, eps_I, eps_B, plan.x);
  
  /* precompute psi, the entries of the matrix B */
//   if(plan.nfct_flags & PRE_PSI)
//     nfct_precompute_psi(&plan);

  /** initialise CG method (inverse plan) */
//  solver_init_double(&ip,(mv_plan_double*)(&plan));
  fcs_float *y, *r_iter, *f_hat_iter, *p_hat_iter, *z_hat_iter, *v_iter;

  y          = (double*)nfft_malloc(ths->mv->M_total*sizeof(double));
  r_iter     = (double*)nfft_malloc(ths->mv->M_total*sizeof(double));
  f_hat_iter = (double*)nfft_malloc(ths->mv->N_total*sizeof(double));
  p_hat_iter = (double*)nfft_malloc(ths->mv->N_total*sizeof(double));

  z_hat_iter = (fcs_float*) malloc(sizeof(fcs_float) * N[0]);
  v_iter     = (fcs_float*) malloc(sizeof(fcs_float) * M);




  /** init samples one which we want to minimize the error with the CG method */
  /* init_samples_1d is a wrapper to different node distributions */
  init_samples_1d(plan.M_total, plan.x, ip.y, p, eps_I, eps_B);
//   nfft_vpr_double(ip.y, plan.M_total, "Given data, vector y");

  fcs_float *f_hat_iter;
  f_hat_iter = (fcs_float) malloc(sizeof(fcs_float)*N);

  /** initialise some guess f_hat_0 */
  for(k=0; k<N; k++)
    f_hat_iter[k]=0;
  
  /** execute the CG method */
  double min_error = 1e5;
  int min_iter=0;

  /* inverse trafo */
//  solver_before_loop_double(&ip);

//   nfft_cp_double(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);
//   NFFT_SWAP_double(ths->r_iter, ths->mv->f);
//   ths->mv->mv_trafo(ths->mv);
//   NFFT_SWAP_double(ths->r_iter, ths->mv->f);
  ndct_trafo_1d(N[0], M, f_hat_iter, x,
     r_iter);

//   nfft_upd_axpy_double(ths->r_iter, -1.0, ths->y, ths->mv->M_total);
  for(fcs_int k=0;k<M;k++)
    r_iter[k] = y[k] - r_iter[k];

//   nfft_cp_double(ths->mv->f, ths->r_iter, ths->mv->M_total);
//   NFFT_SWAP_double(ths->z_hat_iter, ths->mv->f_hat);
//   ths->mv->mv_adjoint(ths->mv);
//   NFFT_SWAP_double(ths->z_hat_iter, ths->mv->f_hat);
  ndct_adj_1d(N[0], M, x, r_iter,
     z_hat_iter);

//   if(ths->flags & CGNR)
//     nfft_cp_double(ths->p_hat_iter, ths->z_hat_iter, ths->mv->N_total);
  for(fcs_int k=0;k<N[0];k++)
    p_hat_iter[k] = z_hat_iter[k];



  for(l=0;l<iter;l++)
  {
    /* check error after several steps */
    NFFT_SWAP_double(ip.f_hat_iter,plan_check.f_hat);
    nfct_trafo(&plan_check);
    NFFT_SWAP_double(ip.f_hat_iter,plan_check.f_hat);
  
    error = calculate_error_1d(plan_check.M_total, plan_check.x, plan_check.f, p, eps_I, eps_B, 0);
    if(error < min_error){
      min_error = error;
      min_iter = l;
      for(int k=0; k<N; k++)
        f_hat_opt[k] = ip.f_hat_iter[k];
    }

    /* give some feedback about convergence */
    if(l%10==0 && verbose){
      fprintf(stderr,"%d. Iteration, Residual ||r||=%e, max. double error = %.6e\n",l,sqrt(ip.dot_r_iter), error);
      fprintf(stderr,"%d. Iteration, beta = %.6e, ||grad|| = %.6e\n", l, ip.beta_iter, sqrt(ip.dot_z_hat_iter));
    }
    
    /* execute one CG step */
//     solver_loop_one_step_double(&ip);

//   nfft_cp_double(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);
//   NFFT_SWAP_double(ths->v_iter,ths->mv->f);
//   ths->mv->mv_trafo(ths->mv);
//   NFFT_SWAP_double(ths->v_iter,ths->mv->f);
  ndct_trafo_1d(N[0], M, p_hat_iter,
      v_iter);

  dot_v_iter = scalar_product(M, v_iter, v_iter);


  dot_v_iter = 0;
  for(fcs_int j=0; j<M; j++)
    dot_v_iter += v_iter[k] * v_iter[k];

  /*-----------------*/
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /*-----------------*/
  nfft_upd_xpay_double(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  nfft_upd_xpay_double(ths->r_iter, -ths->alpha_iter, ths->v_iter,
			ths->mv->M_total);

  ths->dot_r_iter = nfft_dot_double(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  nfft_cp_double(ths->mv->f, ths->r_iter, ths->mv->M_total);

  NFFT_SWAP_double(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  NFFT_SWAP_double(ths->z_hat_iter,ths->mv->f_hat);

  ths->dot_z_hat_iter_old = ths->dot_z_hat_iter;
  ths->dot_z_hat_iter = nfft_dot_double(ths->z_hat_iter, ths->mv->N_total);

  /*-----------------*/
  ths->beta_iter = ths->dot_z_hat_iter / ths->dot_z_hat_iter_old;

  /*-----------------*/
  nfft_upd_axpy_double(ths->p_hat_iter, ths->beta_iter, ths->z_hat_iter,
			ths->mv->N_total);












  }

  for(int k=0; k<N; k++)
    f_hat_opt[k] = f_hat_iter[k];

  /* free mem */
//   solver_finalize_double(&ip);
  free(z_hat_iter);
  free(v_iter);
//   nfct_finalize(&plan);

  return min_error;
}





static double kernel(double r, int p, double eps_I, double eps_B)
{
  double c=2.063296; /* optimal for a=b=4 with p=5 */
//   double c=2.1311445; /* optimal for a=b=8 with p=7 */
  return ifcs_p2nfft_regkernel(ifcs_p2nfft_one_over_modulus, r, p, &c, eps_I, eps_B);
}

static int init_nodes_1d(int M, double eps_I, double eps_B, double *x){
  return init_chebyshev_nodes_1d(M, eps_I, eps_B, x);
}

static int count_nodes_1d(int M, double eps_I, double eps_B){
  return init_nodes_1d(M, eps_I, eps_B, NULL);
}

static int count_equispaced_nodes_1d(int M, double eps_I, double eps_B){
  return init_equispaced_nodes_1d(M, eps_I, eps_B, NULL);
}


static int init_chebyshev_nodes_1d(int M, double eps_I, double eps_B, double* x){
  int counter=0;
  double x0, a, b;
#if ONLY_FARFIELD
    int only_farfield = 1;
#else
    int only_farfield = 0;
#endif
  
  /** Chebyshev nodes */
//   X = (1:n+1)';
//   X = cos((2*X-1)./(2*n+2)*pi); => X \in (-1,1)
//   X = 0.5*(a+b)+0.5*(b-a)*X;    => X \in (a,b)
  for(int k=0; k<M; k++){
    x0 = cos((2.0*k+1)/(2*M)*PI); /* x0 \in (-1,1) */
    a = eps_I;
    b = only_farfield ? 0.5-eps_B : 0.5;
    x0 = 0.5*(a+b)+0.5*(b-a)*x0; /* x0 \in (a,b) */
    
    if( x0 <= b){
      if( x != NULL )
        x[counter] = x0;
      counter++;
    }
  }
  return counter;
}


static void init_samples_1d(
    int num_nodes, double* x, double* y,
    int p, double eps_I, double eps_B
    )
{
  double r;
  
  for(int k=0; k<num_nodes; k++){
    r = sqrt(x[k]*x[k]);
    y[k] = kernel(r, p, eps_I, eps_B);
  }
}

static double calculate_error_1d(
    int num_nodes, double* x, double* f,
    int p, double eps_I, double eps_B,
    int verbose
  )
{
  double r, error, maxerror=0, error_sqr=0;
  
  for(int k=0; k<num_nodes; k++){
    r = sqrt(x[k]*x[k]);
#if NO_FULL_ERROR_CALC
    if(r < 0.5-eps_B)
#endif
    {
      error = fabs(kernel(r, p, eps_I, eps_B) - f[k]);
      if( error > maxerror )
        maxerror = error;
      error_sqr += (fabs(1/r-f[k])) * (fabs(1/r-f[k]));
    }
  }
  error_sqr = sqrt(error_sqr);
  if(verbose)
    fprintf(stderr, "square_error = %.6e\n", error_sqr);
  return maxerror;
}


