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
#include "bessel_i0.h"
#include "bspline.h"
#include "sinc.h"

//static R PNX(fmax)(R arg1, R arg2){ return (arg1>arg2) ? arg1 : arg2; }

/* The factor 1/n from matrix D cancels with the factor n of the inverse Fourier coefficients. */
#define PNFFT_INV_PHI_HAT_GAUSS(k,n,b) \
  (pnfft_exp( PNFFT_PI*PNFFT_PI*(k)*(k)*1.0/((n)*(n))*(b) ))

#define PNFFT_PHI_HAT_GAUSS(k,n,b) \
  (pnfft_exp(-PNFFT_PI*PNFFT_PI*(k)*(k)*1.0/((n)*(n))*(b) ))

/* define this macro to keep the notation short for Kaiser-Bessel and BesselI0 windows */
#define PNFFT_KAISER_SQRT1(k,n,b) ( pnfft_sqrt( PNFFT_SQR(b) - PNFFT_SQR(2.0*PNFFT_PI*(k)/(n)) ) )
#define PNFFT_KAISER_SQRT2(k,n,b) ( pnfft_sqrt( PNFFT_SQR(2.0*PNFFT_PI*(k)/(n)) - PNFFT_SQR(b) ) )

/* For oversampling factor sigma==1 we have to assure that the argument of the square root is not less than zero.
 * This is theoretically true (since k runs between -n/2 and n/2-1), but can be numerically wrong for k=-n/2. */
/* The factor 1/n from matrix D cancels with the factor n of the inverse Fourier coefficients. */
#define PNFFT_INV_PHI_HAT_KAISER(k,N,n,b,m) \
  ( (PNFFT_ABS(k) >= (n) - (N)/2) ? 0.0 : 1.0 / PNX(bessel_i0)( (m) * PNFFT_KAISER_SQRT1(k,n,b) ) )

#define PNFFT_PHI_HAT_KAISER(k,N,n,b,m) \
  ( (PNFFT_ABS(k) >= (n) - (N)/2) ? 0.0 : PNX(bessel_i0)( (m) * PNFFT_KAISER_SQRT1(k,n,b) ) )


/* #define PNFFT_INV_PHI_HAT_BESSEL_I0(k,N,n,b,m) \
   ( (PNFFT_ABS(k) == (n) - (N)/2) ? 1.0/(m) : PNFFT_KAISER_SQRT1(k,n,b) / pnfft_sinh( (m) * PNFFT_KAISER_SQRT1(k,n,b) ) ) */

#define PNFFT_INV_PHI_HAT_BESSEL_I0(k,N,n,b,m) \
  ( (PNFFT_ABS(k) < (n) - (N)/2) \
    ? PNFFT_KAISER_SQRT1(k,n,b) / pnfft_sinh( (m) * PNFFT_KAISER_SQRT1(k,n,b) ) \
    : (PNFFT_ABS(k) > (n) - (N)/2) \
    ? PNFFT_KAISER_SQRT2(k,n,b) / pnfft_sin( (m) * PNFFT_KAISER_SQRT2(k,n,b) ) \
    : 1.0 / (m) \
  )

#define PNFFT_PHI_HAT_BESSEL_I0(k,N,n,b,m) \
  ( (PNFFT_ABS(k) < (n) - (N)/2) \
    ? pnfft_sinh( (m) * PNFFT_KAISER_SQRT1(k,n,b) ) / PNFFT_KAISER_SQRT1(k,n,b) \
    : (PNFFT_ABS(k) > (n) - (N)/2) \
    ? pnfft_sin( (m) * PNFFT_KAISER_SQRT2(k,n,b) ) / PNFFT_KAISER_SQRT2(k,n,b) \
    : (m) \
  )

/* Usage of MAX:
 * For oversampling factor sigma==1 we have to assure that the argument of the square root is not less than zero.
 * This is theoretically true (since k runs between -n/2 and n/2-1), but can be numerically wrong for k=-n/2.
 * Attention: Usage of PNFFT_MAX in this macro results in NANs on Jugene with -O3 optimization  */
/*
#define PNFFT_INV_PHI_HAT_KAISER(k,n,b,m) \
  ( 1.0 / PNX(bessel_i0)( (m) * pnfft_sqrt( PNX(fmax)( 0.0, PNFFT_SQR(b) - PNFFT_SQR(2.0*PNFFT_PI*(k)/(n))))))
*/

/* The factor 1/n from matrix D cancels with the factor n of the inverse Fourier coefficients. */
#define PNFFT_INV_PHI_HAT_BSPLINE(k,n,m) \
  ( pnfft_pow( PNX(sinc)((k) * PNFFT_PI / (n)), K(-2.0) * (m)) )

#define PNFFT_PHI_HAT_BSPLINE(k,n,m) \
  ( pnfft_pow( PNX(sinc)((k) * PNFFT_PI / (n)), K(2.0) * (m)) )

/* For oversampling factor sigma==1 avoid division by zero. */
/* The factor 1/n from matrix D is computed in matrix B (There it cancels with the factor N of the window). */
#define PNFFT_INV_PHI_HAT_SINC_POWER(k,n,b,m,spline_coeffs) \
  ( ((k) == -(n)/2) ? 0.0 : 1.0 / PNX(bspline)(2 * (m), (R)(k) * (b) / ((R) n) + (R)(m), (spline_coeffs)) )

#define PNFFT_PHI_HAT_SINC_POWER(k,n,b,m,spline_coeffs) \
  ( (abs(k) >= (n)/2) ? 0.0 : PNX(bspline)(2 * (m), (R)(k) * (b) / ((R) n) + (R)(m), (spline_coeffs)) )


static void convolution_with_general_window(
    const C *in,
    const INT *n, const INT *no,
    const INT *local_N, const INT *local_N_start,
    const PNX(plan) window_param, int sign,
    C *out);
static void convolution_with_pre_inv_phi_hat(
    const C *in,
    const INT *local_N,
    const C *pre_inv_phi_hat,
    C *out);
static void precompute_inv_phi_hat_general_window(
    const INT *local_N, const INT *local_N_start,
    const PNX(plan) window_param,
    C *pre_inv_phi_hat);
static void twiddle_phi_hat(
    const INT *n, const INT *no,
    const INT *local_N, const INT *local_N_start,
    int sign,
    C *pre_phi_hat);

/* Return the inverse window Fourier coefficients.
 * Since we use tensor product structure, only the return the factor that belongs to dimension 'dim'. */
R PNX(inv_phi_hat)(
    const PNX(plan) ths, int dim, INT k
    )
{
  if(ths->pnfft_flags & PNFFT_WINDOW_GAUSSIAN)
    return PNFFT_INV_PHI_HAT_GAUSS(k, ths->n[dim], ths->b[dim]);
  else if(ths->pnfft_flags & PNFFT_WINDOW_BSPLINE)
    return PNFFT_INV_PHI_HAT_BSPLINE(k, ths->n[dim], ths->m);
  else if(ths->pnfft_flags & PNFFT_WINDOW_SINC_POWER)
    return PNFFT_INV_PHI_HAT_SINC_POWER(k, ths->n[dim], ths->b[dim], ths->m, ths->spline_coeffs);
  else if(ths->pnfft_flags & PNFFT_WINDOW_BESSEL_I0)
    return PNFFT_INV_PHI_HAT_BESSEL_I0(k, ths->N[dim], ths->n[dim], ths->b[dim], ths->m);
  else
    return PNFFT_INV_PHI_HAT_KAISER(k, ths->N[dim], ths->n[dim], ths->b[dim], ths->m);
}

R PNX(phi_hat)(
    const PNX(plan) ths, int dim, INT k
    )
{
  if(ths->pnfft_flags & PNFFT_WINDOW_GAUSSIAN)
    return PNFFT_PHI_HAT_GAUSS(k, ths->n[dim], ths->b[dim]);
  else if(ths->pnfft_flags & PNFFT_WINDOW_BSPLINE)
    return PNFFT_PHI_HAT_BSPLINE(k, ths->n[dim], ths->m);
  else if(ths->pnfft_flags & PNFFT_WINDOW_SINC_POWER)
    return PNFFT_PHI_HAT_SINC_POWER(k, ths->n[dim], ths->b[dim], ths->m, ths->spline_coeffs);
  else if(ths->pnfft_flags & PNFFT_WINDOW_BESSEL_I0)
    return PNFFT_PHI_HAT_BESSEL_I0(k, ths->N[dim], ths->n[dim], ths->b[dim], ths->m);
  else
    return PNFFT_PHI_HAT_KAISER(k, ths->N[dim], ths->n[dim], ths->b[dim], ths->m);
}

void PNX(trafo_D)(
    PNX(plan) ths
    )
{
#if PNFFT_ENABLE_DEBUG
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  C csum, gcsum;

  csum = 0.0;
  for(INT t=0; t<ths->local_N[0]*ths->local_N[1]*ths->local_N[2]; t++)
    csum += pnfft_fabs(pnfft_creal(ths->f_hat[t])) + _Complex_I * pnfft_fabs(pnfft_cimag(ths->f_hat[t])) ;
  MPI_Reduce(&csum, &gcsum, 2, PNFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!myrank) fprintf(stderr, "PNFFT: Sum of Fourier coefficients before deconvolution: %e + I* %e\n", pnfft_creal(gcsum), pnfft_cimag(gcsum));
#endif

  /* use precomputed window Fourier coefficients if possible */
  if(ths->pnfft_flags & PNFFT_PRE_PHI_HAT){
    convolution_with_pre_inv_phi_hat(
        ths->f_hat, ths->local_N, ths->pre_inv_phi_hat_trafo,
        ths->g1);
  } else {
    convolution_with_general_window(
        ths->f_hat, ths->n, ths->no, ths->local_N, ths->local_N_start, ths, FFTW_FORWARD,
        ths->g1);
  }
}


void PNX(adjoint_D)(
    PNX(plan) ths
    )
{
  /* use precomputed window Fourier coefficients if possible */
  if(ths->pnfft_flags & PNFFT_PRE_PHI_HAT){
    convolution_with_pre_inv_phi_hat(
        ths->g1, ths->local_N, ths->pre_inv_phi_hat_adj,
        ths->f_hat);
  } else {
    convolution_with_general_window(
        ths->g1, ths->n, ths->no, ths->local_N, ths->local_N_start, ths, FFTW_BACKWARD,
        ths->f_hat);
  }

#if PNFFT_ENABLE_DEBUG
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  C csum, gcsum;

  csum = 0.0;
  for(INT t=0; t<ths->local_N[0]*ths->local_N[1]*ths->local_N[2]; t++)
    csum += pnfft_fabs(pnfft_creal(ths->f_hat[t])) + _Complex_I * pnfft_fabs(pnfft_cimag(ths->f_hat[t])) ;
  MPI_Reduce(&csum, &gcsum, 2, PNFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!myrank) fprintf(stderr, "PNFFT^H: Sum of Fourier coefficients after deconvolution: %e + I* %e\n", pnfft_creal(gcsum), pnfft_cimag(gcsum));
#endif
}

static void convolution_with_general_window(
    const C *in,
    const INT *n, const INT *no,
    const INT *local_N, const INT *local_N_start,
    const PNX(plan) window_param, int sign,
    C *out
    )
{
  INT k0, k1, k2, k=0;
  C inv_phi_k1, inv_phi_k1k2, inv_phi_k1k2k0;

  /* g_hat is transposed N1 x N2 x N0 */
  for(k1=local_N_start[1]; k1<local_N_start[1] + local_N[1]; k1++){
    inv_phi_k1 = PNX(inv_phi_hat)(window_param, 1, k1) * pnfft_cexp(-sign * PNFFT_PI * _Complex_I * ( k1*no[1]/((R)n[1]) ) );
    for(k2=local_N_start[2]; k2<local_N_start[2] + local_N[2]; k2++){
      inv_phi_k1k2 = inv_phi_k1 * PNX(inv_phi_hat)(window_param, 2, k2) * pnfft_cexp(-sign * PNFFT_PI * _Complex_I * ( k1*no[2]/((R)n[2]) ) );
      for(k0=local_N_start[0]; k0<local_N_start[0] + local_N[0]; k0++, k++){
        inv_phi_k1k2k0 = inv_phi_k1k2 * PNX(inv_phi_hat)(window_param, 0, k0) * pnfft_cexp(-sign * PNFFT_PI * _Complex_I * ( k1*no[0]/((R)n[0]) ) );
        out[k] = in[k] * inv_phi_k1k2k0;
      }
    }
  }
}

static void convolution_with_pre_inv_phi_hat(
    const C *in,
    const INT *local_N,
    const C *pre_inv_phi_hat,
    C *out
    )
{
  INT k0, k1, k2, k=0;
  const C *inv_phi_hat0 = pre_inv_phi_hat;
  const C *inv_phi_hat1 = inv_phi_hat0 + local_N[0];
  const C *inv_phi_hat2 = inv_phi_hat1 + local_N[1];

  /* g_hat is transposed N1 x N2 x N0 */
  for(k1=0; k1<local_N[1]; k1++)
    for(k2=0; k2<local_N[2]; k2++)
      for(k0=0; k0<local_N[0]; k0++, k++)
        out[k] = in[k] * inv_phi_hat0[k0] * inv_phi_hat1[k1] * inv_phi_hat2[k2];
}

void PNX(precompute_inv_phi_hat_trafo)(
    PNX(plan) ths,
    C *pre_phi_hat_trafo
    )
{
  precompute_inv_phi_hat_general_window(ths->local_N, ths->local_N_start, ths,
      pre_phi_hat_trafo);
  twiddle_phi_hat(ths->n, ths->no, ths->local_N, ths->local_N_start, FFTW_FORWARD,
      pre_phi_hat_trafo);
}

void PNX(precompute_inv_phi_hat_adj)(
    PNX(plan) ths,
    C *pre_phi_hat_adj
    )
{
  precompute_inv_phi_hat_general_window(ths->local_N, ths->local_N_start, ths,
      pre_phi_hat_adj);
  twiddle_phi_hat(ths->n, ths->no, ths->local_N, ths->local_N_start, FFTW_BACKWARD,
      pre_phi_hat_adj);
}

static void precompute_inv_phi_hat_general_window(
    const INT *local_N, const INT *local_N_start,
    const PNX(plan) window_param,
    C *pre_inv_phi_hat
    )
{
  INT l=0;

  for(INT t=0; t<3; t++)
    for(INT k=local_N_start[t]; k<local_N_start[t] + local_N[t]; k++, l++)
      pre_inv_phi_hat[l] = PNX(inv_phi_hat)(window_param, t, k);
}

static void twiddle_phi_hat(
    const INT *n, const INT *no,
    const INT *local_N, const INT *local_N_start,
    int sign,
    C *pre_inv_phi_hat
    )
{
  INT l=0;

  for(INT t=0; t<3; t++)
    for(INT k=local_N_start[t]; k<local_N_start[t] + local_N[t]; k++, l++)
      pre_inv_phi_hat[l] *= pnfft_cexp(-sign * PNFFT_PI * _Complex_I * ( k*no[t]/((R)n[t]) ) );
}

