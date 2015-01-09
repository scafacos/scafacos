/*
 * Copyright (C) 2011-2014 Michael Pippig
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include "taylor2p.h"
#include "part_derive_one_over_norm_x.h"

/* Switch to use Lagrange instead of Newton basis polynomials for interpolation.
 * Default are Newton polynomials, since they are faster. */
#define FCS_P2NFFT_INTPOL_LAGRANGE 1



#if FCS_P2NFFT_INTPOL_LAGRANGE

/*****************************************************
 * Interpolation polynomials based on Lagrange basis *
 *****************************************************/

/** factorial */
static fcs_float fak(fcs_int n)
{
  if (n<=1) return 1.0;
  else return (fcs_float)n*fak(n-1);
}

/** binomial coefficient */
static fcs_float binom(fcs_int n, fcs_int m)
{
  return fak(n)/fak(m)/fak(n-m);
}

/** basis polynomial for regularized kernel */
static fcs_float BasisPoly(fcs_int m, fcs_int r, fcs_float xx)
{
  fcs_int k;
  fcs_float sum=0.0;

//   if(fcs_float_is_zero(xx+1.0))
//     return 1.0;

  for (k=0; k<=m-r; k++) {
    sum+=binom(m+k,k)*fcs_pow((xx+1.0)/2.0,(fcs_float)k);
  }
  return sum*fcs_pow((xx+1.0),(fcs_float)r)*fcs_pow(1.0-xx,(fcs_float)(m+1))/(1<<(m+1))/fak(r); /* 1<<(m+1) = 2^(m+1) */
}

/** integrated basis polynomial for regularized kernel */
static fcs_float IntBasisPoly(fcs_int p, fcs_int j, fcs_float y)
{
  fcs_int k,l;
  fcs_float sum1=0.0, sum2=0.0;
  
  for (l=0; l<=p; l++) {
    sum1 = 0.0;
    for (k=0; k<=p-j-1; k++) {
      sum1 += binom(p+k-1,k)*fak(j+k)/fak(j+k+1+l)*fcs_pow((1.0+y)/2.0,(fcs_float)k);
    }
    sum2 += fcs_pow(1.0+y,(fcs_float)l)*fcs_pow(1.0-y,(fcs_float)(p-l))/fak(p-l)*sum1;
  }
  return sum2 * fak(p)/fak(j)/(1<<p)*fcs_pow(1.0+y,(fcs_float)(j+1)); /* 1<<p = 2^p */
}

/** Use Lagrange basis polynomials for evaluating the interpolation polynomial P(x) that satisfies
 *  the odd symmetric Hermite interpolation conditions in two nodes up to the p-th derivative, i.e.,
 *      P^k(x0) = (-1)^k P^k(x1) = kernel^k(x0), for all k=0,...,p-1
 */
static fcs_float interpolate_lagrange_sym(
    const fcs_float *kernel, fcs_int p,
    fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float sum=0.0;

  /* canonicalize x to y \in [-1,1] */
  fcs_float r = 0.5*(x1-x0);
  fcs_float m = 0.5*(x1+x0);
  fcs_float y = (x-m)/r;

  for (fcs_int i=0; i<p; i++)
    sum += (BasisPoly(p-1,i,y) + BasisPoly(p-1,i,-y)) * fcs_pow(r,(fcs_float)i) * kernel[i];

  return sum;
}

/** Use Lagrange basis polynomials for evaluating the interpolation polynomial P(x) that satisfies
 *  the non-symmetric Hermite interpolation conditions in two nodes up to the p-th derivative, i.e.,
 *      P^k(x0) = kernel^k(x0), for all k=0,...,p-1,
 *      P(x1)   = c,
 *      p^k(x1) = 0,            for all k=1,...,p-1
 */
static fcs_float interpolate_lagrange_ec(
    const fcs_float *kernel, fcs_int p,
    fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float sum=0.0;

  /* canonicalize x to y \in [-1,1] */
  fcs_float r = 0.5*(x1-x0);
  fcs_float m = 0.5*(x1+x0);
  fcs_float y = (x-m)/r;

  for (fcs_int i=0; i<p; i++)
    sum += BasisPoly(p-1,i,y) * fcs_pow(r,(fcs_float)i) * kernel[i];

  return sum + kernel[p] * BasisPoly(p-1,0,-y);
}

/** Use Lagrangian basis polynomials for evaluating the interpolation polynomial P(x) that satisfies
 *  the non-symmetric Hermite interpolation conditions in two nodes up to the p-th derivative, i.e.,
 *      P^k(x0) = kernel^k(x0), for all k=0,...,p-1,
 *      P(x1) is not determined
 *      P^k(x1) = 0,            for all k=1,...,p-1,
 */
static fcs_float interpolate_lagrange_ic(
    const fcs_float *kernel, fcs_int p,
    fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float sum=0.0;

  /* canonicalize x to y \in [-1,1] */
  fcs_float r = 0.5*(x1-x0);
  fcs_float m = 0.5*(x1+x0);
  fcs_float y = (x-m)/r;

  for (fcs_int j=0; j<=p-2; j++)
    sum += fcs_pow(r,j+1.0) * kernel[j+1]
      * (IntBasisPoly(p-1,j,y) - IntBasisPoly(p-1,j,-1));

  return sum + kernel[0];
}

#else

/***************************************************
 * Interpolation polynomials based on Newton basis *
 ***************************************************/

/* Evaluate polynomial of the type \sum_{j=0}^{2p} w_j(x), where 
 * w_0(x) := 1
 * w_j(x) := \prod_{i=0}^{j-1} x-z_j, for p=1,..,2p-1
 * and
 * z_j := x_0, for j=0,..,p
 * z_j := x_1, for j=p+1,..,2p */
static fcs_float eval_newton_horner(
    fcs_int p, const fcs_float *coeff, fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float res=coeff[2*p+1];

  for(int i=2*p; i>p; --i)
    res = res * (x-x1) + coeff[i];
  for(int i=p; i>=0; --i)
    res = res * (x-x0) + coeff[i];

  return res;
}

/* Evaluate polynomial of the type \sum_{j=0}^{2p} w_j(x), where 
 * w_0(x) := 1
 * w_j(x) := \prod_{i=0}^{j-1} x-z_j, for p=1,..,2p
 * and
 * z_j := x_0, for j=0,..,p+1
 * z_j := x_1, for j=p+2,..,2p */
fcs_float eval_newton_horner_ic(
    fcs_int p, const fcs_float *coeff, fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float res=coeff[2*p];

  for(int i=2*p-1; i>p; --i)
    res = res * (x-x1) + coeff[i];
  for(int i=p; i>=0; --i)
    res = res * (x-x0) + coeff[i];

  return res;
}

/** Compute the coefficients of the Newton basis polynomials corresponding to 
 *  the interpolation polynomial P(x) that satisfies the Hermite interpolation
 *  conditions in two nodes up to the p-th derivative, i.e.,
 *      P^k(x0) = kernel^k(x0), P^k(x1) = kernel^k(x1), for all k=0,...,p
 */
static void newton_coeff(
    const fcs_float *kernel, fcs_int p,
    fcs_float x0, fcs_float x1, 
    fcs_float *coeff
    )
{
  fcs_float *buf = (fcs_float*) malloc(sizeof(fcs_float) * (p+2)); 

  /* init phase 1
   * Start iteration with functions values at x0 and x1. */
  coeff[0] = buf[p] = kernel[0];
  buf[p+1] = kernel[p+1];

  fcs_float one_over_x01 = 1.0 / (x1-x0);
  fcs_float one_over_fak_k = 1.0;

  /* phase 1 (expand coefficients from 2 to p+2):
   * Start with the 2 function values at x0 and x1 in step k=0.
   * In step k we use a Neville-Hermite recursion to compute the next k+2 from the k+1 preceding values (builds up a triangle structure).
   * The values of the k-th derivatives enter the recursion in step k. */
  for(int k=1; k<=p; k++){
    one_over_fak_k /= k;

    /* index i=p-k introduces k-th derivative at x0 */
    buf[p-k] = one_over_fak_k * kernel[k];

    /* use recursion for all values in the middle */
    for(int i=p-k+1; i<=p; i++)
      buf[i] = (buf[i+1] - buf[i]) * one_over_x01;

    /* index i=p+1 introduces k-th derivative at x1 */
    buf[p+1] = one_over_fak_k * kernel[k+p+1];

    /* save coefficients of interpolation polynomial in Newton basis */
    coeff[k] = buf[p-k];
  }

  /* phase 2 (reduce coefficients from p+2 to 1)
   * After phase 1 all the derivatives have entered the Neville-table.
   * Now use the recursion to reduce the 2p values by 1 in every step. */
  for(int k=p+1; k<=2*p+1; k++){
    for(int i=0; i<=2*p+1-k; i++)
      buf[i] = (buf[i+1] - buf[i]) * one_over_x01;
    coeff[k] = buf[0];
  }

  free(buf);
}

/** Compute the coefficients of the Newton basis polynomials corresponding to 
 *  the interpolation polynomial P(x) that satisfies the following Hermite interpolation
 *  conditions in two nodes up to the p-th derivative, i.e.,
 *      P^k(x0) = kernel^k(x0), for all k=0,...,p
 *      P(x1) is not determined
 *      P^k(x1) = kernel^k(x1), for all k=1,...,p
 */
static void newton_coeff_ic(
    const fcs_float *kernel, fcs_int p,
    fcs_float x0, fcs_float x1, 
    fcs_float *coeff
    )
{
  fcs_float x01 = x1-x0;

  /* At first, compute the interpolating derivative of Q(x).
   * Therefore, we compute the Newton basis expansion of the unique
   * polynomial P(x) that fulfills
   * P^(j)(x_0) = a_{j+1}, P^(j)(x_1) = b_{j+1}, p=0,..,p-1. */
  newton_coeff(kernel+1, p-1, x0, x1,
      coeff);

  /* Compute the coefficients of Q(x) from the coefficients of its derivative P(x).
   * This 
   * ) in the Newton basis. */
  if(p>0)
    coeff[2*p] = coeff[2*p-1] / (2*p);
  for(int j=2*p-1; j>p; --j)
    coeff[j] = ( coeff[j-1] - (j-p)*x01*coeff[j+1] ) / j;
  for(int j=p; j>0; --j)
    coeff[j] = coeff[j-1] / j;
    
  /* Finally, set the function value at x0. */
  coeff[0] = kernel[0];
}

static fcs_float interpolate_newton(
    const fcs_float *kernel, fcs_int p,
    fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float Px;
  fcs_float *coeff = (fcs_float*) malloc(sizeof(fcs_float) * (2*p+2)); 

  /* compute the coefficients of P(x) in Newton basis */
  newton_coeff(kernel, p, x0, x1,
      coeff);

  /* evaluate P(x) via Horner scheme */
  Px = eval_newton_horner(p, coeff, x0, x1, x);

  free(coeff);
  return Px;
}

static fcs_float interpolate_newton_ic(
    const fcs_float *kernel, fcs_int p,
    fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float Px;
  fcs_float *coeff = (fcs_float*) malloc(sizeof(fcs_float) * (2*p+1)); 

  /* compute the coefficients of P(x) in Newton basis */
  newton_coeff_ic(kernel, p, x0, x1,
      coeff);

  /* evaluate P(x) via Horner scheme */
  Px = eval_newton_horner_ic(p, coeff, x0, x1, x);

  free(coeff);
  return Px;
}

#endif


/****************************************************
 * Switch between Lagrange and Newton interpolation *
 ****************************************************/

static fcs_float interpolate_sym(
    const fcs_float *kernel, fcs_int p,
    fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float retval;
  fcs_float *buf = (fcs_float*) malloc(sizeof(fcs_float) * (2*p));

  fcs_float minus_one_to_i = 1.0;
  for(int i=0; i<p; i++){
    buf[i]   = kernel[i];  
    buf[i+p] = kernel[i] * minus_one_to_i;
    minus_one_to_i *= -1.0;
  }

#if FCS_P2NFFT_INTPOL_LAGRANGE
  retval = interpolate_lagrange_sym(buf, p, x0, x1, x);
#else
  retval = interpolate_newton(buf, p-1, x0, x1, x);
#endif

  free(buf);
  return retval;
}

/** Use Lagrange or Newton basis polynomials for evaluating the interpolation polynomial P(x) that satisfies
 *  the odd symmetric Hermite interpolation conditions in two nodes up to the p-th derivative, i.e.,
 *      P^k(x0) = (-1)^k P^k(x1) = kernel^k(x0), for all k=0,...,p-1
 */
fcs_float ifcs_p2nfft_interpolate_symmetric(
    ifcs_p2nfft_kernel k, const fcs_float *param,
    fcs_int p, fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float retval;
  fcs_float *kernel = (fcs_float*) malloc(sizeof(fcs_float) * (2*p));

  fcs_float minus_one_to_i = 1.0;
  for(int i=0; i<p; i++){
    kernel[i]   = k(x0, i, param);
    kernel[i+p] = kernel[i] * minus_one_to_i;
    minus_one_to_i *= -1.0;
  }

#if FCS_P2NFFT_INTPOL_LAGRANGE
  retval = interpolate_lagrange_sym(kernel, p, x0, x1, x);
#else
  retval = interpolate_newton(kernel, p-1, x0, x1, x);
#endif

  free(kernel);
  return retval;
}

static fcs_float interpolate_ec(
    const fcs_float *kernel, fcs_float c, fcs_int p,
    fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float retval;
  fcs_float *buf = (fcs_float*) malloc(sizeof(fcs_float) * (2*p));

  for(int i=0; i<p; i++){
    buf[i]     = kernel[i];
    buf[i+p] = 0;
  }
  buf[p] = c;

#if FCS_P2NFFT_INTPOL_LAGRANGE
  retval = interpolate_lagrange_ec(buf, p, x0, x1, x);
#else
  retval = interpolate_newton(buf, p-1, x0, x1, x);
#endif

  free(buf);
  return retval;
}

/** Use Lagrange or Newton basis polynomials for evaluating the interpolation polynomial P(x) that satisfies
 *  the non-symmetric Hermite interpolation conditions in two nodes up to the p-th derivative, i.e.,
 *      P^k(x0) = kernel^k(x0), for all k=0,...,p-1,
 *      P(x1)   = c,
 *      p^k(x1) = 0,            for all k=1,...,p-1
 */
fcs_float ifcs_p2nfft_interpolate_explicit_continuation(
    ifcs_p2nfft_kernel k, const fcs_float *param, fcs_float c,
    fcs_int p, fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float retval;
  fcs_float *kernel = (fcs_float*) malloc(sizeof(fcs_float) * (2*p));

  for(int i=0; i<p; i++){
    kernel[i]     = k(x0, i, param);
    kernel[i+p] = 0;
  }
  kernel[p] = c;

#if FCS_P2NFFT_INTPOL_LAGRANGE
  retval = interpolate_lagrange_ec(kernel, p, x0, x1, x);
#else
  retval = interpolate_newton(kernel, p-1, x0, x1, x);
#endif

  free(kernel);
  return retval;
}

static fcs_float interpolate_ic(
    const fcs_float *kernel, fcs_int p,
    fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float retval;
  fcs_float *buf = (fcs_float*) malloc(sizeof(fcs_float) * (2*p-1));

  for(int i=0; i<p; i++)
    buf[i]   = kernel[i];
  for(int i=0; i<p-1; i++)
    buf[i+p] = 0;

#if FCS_P2NFFT_INTPOL_LAGRANGE
  retval = interpolate_lagrange_ic(buf, p, x0, x1, x);
#else
  retval = interpolate_newton_ic(buf, p-1, x0, x1, x);
#endif

  free(buf);
  return retval;
}

/** Use Lagrange or Newton basis polynomials for evaluating the interpolation polynomial P(x) that satisfies
 *  the non-symmetric Hermite interpolation conditions in two nodes up to the p-th derivative, i.e.,
 *      P^k(x0) = kernel^k(x0), for all k=0,...,p-1,
 *      P(x1) is not determined
 *      P^k(x1) = 0,            for all k=1,...,p-1,
 */
fcs_float ifcs_p2nfft_interpolate_implicit_continuation(
    ifcs_p2nfft_kernel k, const fcs_float *param,
    fcs_int p, fcs_float x0, fcs_float x1, fcs_float x
    )
{
  fcs_float retval;
  fcs_float *kernel = (fcs_float*) malloc(sizeof(fcs_float) * (2*p-1));

  for(int i=0; i<p; i++)
    kernel[i]     = k(x0, i, param);
  for(int i=0; i<p-1; i++)
    kernel[i+p] = 0;

#if FCS_P2NFFT_INTPOL_LAGRANGE
  retval = interpolate_lagrange_ic(kernel, p, x0, x1, x);
#else
  retval = interpolate_newton_ic(kernel, p-1, x0, x1, x);
#endif

  free(kernel);
  return retval;
}

/* Regularize the kernel function 1/norm(x) in every dimension i,
 * where xi exceeds (0.5-epsB) * hi, based on one-dimensional symmetric two-point Taylor polynomials, i.e.
 *   d^j/dx^j K(mi-ri) = P^(j)(mi-ri) = (-1)^j P^(j)(mi+ri)
 * We nest these one-dimensional two-point Taylor polynomials in order to construct multi-dimensional ones.
 * The number of variables of P is equal to the number of dimensions with xi > (0.5-epsB) * hi.
 */
/** Construct 
 *  1d interpolation polynomials at the sides,
 *  2d interpolation polynomials on the edges, and
 *  3d interpolation polynomials on the corners
 *  of a cuboid domain. */
fcs_float ifcs_p2nfft_interpolate_cuboid_symmetric(
    fcs_int p, const fcs_float *x0, const fcs_float *x1, const fcs_float *x_signed
    )
{
  fcs_int reg_dims[3];
  fcs_float z[3], x[3];

  for(fcs_int t=0; t<3; t++){
    /* Due to symmetry we can use fabs(x) and regularize only the positive octant */
    x[t] = fcs_fabs(x_signed[t]);

    /* decide which coordinates are located in the regularization domain */
    reg_dims[t] = (x[t] > x0[t]);

    /* Evaluation point of the kernel function */
    z[t] = reg_dims[t] ? x0[t] : x[t];
  }

  /* Try to implement this with Newton interpolation */
  fcs_float *buf = (fcs_float*) malloc(sizeof(fcs_float) * p*p*p);

  fcs_int l = 0;
  for (fcs_int i0=0; i0<p; ++i0) {
    for (fcs_int i1=0; i1<p; ++i1) {
      for (fcs_int i2=0; i2<p; ++i2) {
        buf[l++] = ifcs_p2nfft_part_derive_one_over_norm_x(i0,i1,i2,z[0],z[1],z[2]);
        if(!reg_dims[2]) break;
      }

      if(reg_dims[2]){
        l -= p;
        buf[l] = interpolate_sym(buf+l, p, x0[2], x1[2], x[2]);
        ++l;
      }

      if(!reg_dims[1]) break;
    }

    if(reg_dims[1]){
      l -= p;
      buf[l] = interpolate_sym(buf+l, p, x0[1], x1[1], x[1]);
      ++l;
    }

    if(!reg_dims[0]) break;
  }

  if(reg_dims[0]){
    l -= p;
    buf[l] = interpolate_sym(buf+l, p, x0[0], x1[0], x[0]);
  }

  fcs_float retval = buf[0];

  free(buf);
  return retval;
}

fcs_float ifcs_p2nfft_interpolate_cuboid_explicit_continuation(
    fcs_float c, fcs_int p, const fcs_float *x0, const fcs_float *x1, const fcs_float *x_signed
    )
{
  fcs_int reg_dims[3];
  fcs_float z[3], x[3];

  for(fcs_int t=0; t<3; t++){
    /* Due to symmetry we can use fabs(x) and regularize only the positive octant */
    x[t] = fcs_fabs(x_signed[t]);

    /* decide which coordinates are located in the regularization domain */
    reg_dims[t] = (x[t] > x0[t]);

    /* Evaluation point of the kernel function */
    z[t] = reg_dims[t] ? x0[t] : x[t];
  }

  /* Try to implement this with Newton interpolation */
  fcs_float *buf = (fcs_float*) malloc(sizeof(fcs_float) * p*p*p);

  fcs_int l = 0;
  for (fcs_int i0=0; i0<p; ++i0) {
    for (fcs_int i1=0; i1<p; ++i1) {
      for (fcs_int i2=0; i2<p; ++i2) {
        buf[l++] = ifcs_p2nfft_part_derive_one_over_norm_x(i0,i1,i2,z[0],z[1],z[2]);
        if(!reg_dims[2]) break;
      }

      if(reg_dims[2]){
        l -= p;
        buf[l] = interpolate_ec(buf+l, c, p, x0[2], x1[2], x[2]);
        ++l;
      }

      if(!reg_dims[1]) break;
    }

    if(reg_dims[1]){
      l -= p;
      buf[l] = interpolate_ec(buf+l, c, p, x0[1], x1[1], x[1]);
      ++l;
    }

    if(!reg_dims[0]) break;
  }

  if(reg_dims[0]){
    l -= p;
    buf[l] = interpolate_ec(buf+l, c, p, x0[0], x1[0], x[0]);
  }

  fcs_float retval = buf[0];

  free(buf);
  return retval;
}

fcs_float ifcs_p2nfft_interpolate_cuboid_implicit_continuation(
    fcs_int p, const fcs_float *x0, const fcs_float *x1, const fcs_float *x_signed
    )
{
  fcs_int reg_dims[3];
  fcs_float z[3], x[3];

  for(fcs_int t=0; t<3; t++){
    /* Due to symmetry we can use fabs(x) and regularize only the positive octant */
    x[t] = fcs_fabs(x_signed[t]);

    /* decide which coordinates are located in the regularization domain */
    reg_dims[t] = (x[t] > x0[t]);

    /* Evaluation point of the kernel function */
    z[t] = reg_dims[t] ? x0[t] : x[t];
  }

  /* Try to implement this with Newton interpolation */
  fcs_float *buf = (fcs_float*) malloc(sizeof(fcs_float) * p*p*p);

  fcs_int l = 0;
  for (fcs_int i0=0; i0<p; ++i0) {
    for (fcs_int i1=0; i1<p; ++i1) {
      for (fcs_int i2=0; i2<p; ++i2) {
        buf[l++] = ifcs_p2nfft_part_derive_one_over_norm_x(i0,i1,i2,z[0],z[1],z[2]);
        if(!reg_dims[2]) break;
      }

      if(reg_dims[2]){
        l -= p;
        buf[l] = interpolate_ic(buf+l, p, x0[2], x1[2], x[2]);
        ++l;
      }

      if(!reg_dims[1]) break;
    }

    if(reg_dims[1]){
      l -= p;
      buf[l] = interpolate_ic(buf+l, p, x0[1], x1[1], x[1]);
      ++l;
    }

    if(!reg_dims[0]) break;
  }

  if(reg_dims[0]){
    l -= p;
    buf[l] = interpolate_ic(buf+l, p, x0[0], x1[0], x[0]);
  }

  fcs_float retval = buf[0];

  free(buf);
  return retval;
}




/**********************************************************
 * Precomputed coefficients of 2-point Taylor polynomials *
 * for 1/abs(x) and its derivative                        *
 **********************************************************/

/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth into x=1/2 (used in dD, d>1)
 */
fcs_float ifcs_p2nfft_nearfield_correction_taylor2p(
    fcs_float xx, fcs_int p, const fcs_float *param
    )
{
  fcs_float horner = param[0];

  xx=fabs(xx);

  /* use horner scheme for next p-2 coefficients */
  for(fcs_int t=1; t<p; t++)
    horner = xx * xx * horner + param[t];

  return horner;
}

/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth into x=1/2 (used in dD, d>1)
 */
fcs_float ifcs_p2nfft_nearfield_correction_taylor2p_derive(
    fcs_float xx, fcs_int p, const fcs_float *param
    )
{
  fcs_float horner = param[0];

  xx=fabs(xx);

  /* use horner schema for next p-2 coefficients */
  for(fcs_int t=1; t<p-1; t++)
    horner = xx * xx * horner + param[t];

  /* multiply last time by xx, since we have an odd polynomial */
  return horner * xx;
}


fcs_int ifcs_p2nfft_load_taylor2p_coefficients(
   fcs_int p,
   fcs_float *param
   )
{
  /* Case p<2 is senseless, for p>16 we need to add the coefficients. */
  /* These coefficients are calculated analytically with a Maple script. */
  switch(p){
    case(2):
      param[0] = -1.0/2.0;
      param[1] = 3.0/2.0;
      break;
    case(3):
      param[0] = 3.0/8.0;
      param[1] = -5.0/4.0;
      param[2] = 15.0/8.0;
      break;
    case(4):
      param[0] = -5.0/16.0;
      param[1] = 21.0/16.0;
      param[2] = -35.0/16.0;
      param[3] = 35.0/16.0;
      break;
    case(5):
      param[0] = 35.0/128.0;
      param[1] = -45.0/32.0;
      param[2] = 189.0/64.0;
      param[3] = -105.0/32.0;
      param[4] = 315.0/128.0;
      break;
    case(6):
      param[0] = -63.0/256.0;
      param[1] = 385.0/256.0;
      param[2] = -495.0/128.0;
      param[3] = 693.0/128.0;
      param[4] = -1155.0/256.0;
      param[5] = 693.0/256.0;
      break;
    case(7):
      param[0] = 231.0/1024.0;
      param[1] = -819.0/512.0;
      param[2] = 5005.0/1024.0;
      param[3] = -2145.0/256.0;
      param[4] = 9009.0/1024.0;
      param[5] = -3003.0/512.0;
      param[6] = 3003.0/1024.0;
      break;
    case(8):
      param[0] = -429.0/2048.0;
      param[1] = 3465.0/2048.0;
      param[2] = -12285.0/2048.0;
      param[3] = 25025.0/2048.0;
      param[4] = -32175.0/2048.0;
      param[5] = 27027.0/2048.0;
      param[6] = -15015.0/2048.0;
      param[7] = 6435.0/2048.0;
      break;
    case(9):
      param[0] = 6435.0/32768.0;
      param[1] = -7293.0/4096.0;
      param[2] = 58905.0/8192.0;
      param[3] = -69615.0/4096.0;
      param[4] = 425425.0/16384.0;
      param[5] = -109395.0/4096.0;
      param[6] = 153153.0/8192.0;
      param[7] = -36465.0/4096.0;
      param[8] = 109395.0/32768.0;
      break;
    case(10):
      param[0] = -12155.0/65536.0;
      param[1] = 122265.0/65536.0;
      param[2] = -138567.0/16384.0;
      param[3] = 373065.0/16384.0;
      param[4] = -1322685.0/32768.0;
      param[5] = 1616615.0/32768.0;
      param[6] = -692835.0/16384.0;
      param[7] = 415701.0/16384.0;
      param[8] = -692835.0/65536.0;
      param[9] = 230945.0/65536.0;
      break;
    case(11):
      param[0] = 46189.0/262144.0;
      param[1] = -255255.0/131072.0;
      param[2] = 2567565.0/262144.0;
      param[3] = -969969.0/32768.0;
      param[4] = 7834365.0/131072.0;
      param[5] = -5555277.0/65536.0;
      param[6] = 11316305.0/131072.0;
      param[7] = -2078505.0/32768.0;
      param[8] = 8729721.0/262144.0;
      param[9] = -1616615.0/131072.0;
      param[10] = 969969.0/262144.0;
      break;
    case(12):
      param[0] = -88179.0/524288.0;
      param[1] = 1062347.0/524288.0;
      param[2] = -5870865.0/524288.0;
      param[3] = 19684665.0/524288.0;
      param[4] = -22309287.0/262144.0;
      param[5] = 36038079.0/262144.0;
      param[6] = -42590457.0/262144.0;
      param[7] = 37182145.0/262144.0;
      param[8] = -47805615.0/524288.0;
      param[9] = 22309287.0/524288.0;
      param[10] = -7436429.0/524288.0;
      param[11] = 2028117.0/524288.0;
      break;
    case(13):
      param[0] = 676039.0/4194304.0;
      param[1] = -2204475.0/1048576.0;
      param[2] = 26558675.0/2097152.0;
      param[3] = -48923875.0/1048576.0;
      param[4] = 492116625.0/4194304.0;
      param[5] = -111546435.0/524288.0;
      param[6] = 300317325.0/1048576.0;
      param[7] = -152108775.0/524288.0;
      param[8] = 929553625.0/4194304.0;
      param[9] = -132793375.0/1048576.0;
      param[10] = 111546435.0/2097152.0;
      param[11] = -16900975.0/1048576.0;
      param[12] = 16900975.0/4194304.0;
      break;
    case(14):
      param[0] = -1300075.0/8388608.0;
      param[1] = 18253053.0/8388608.0;
      param[2] = -59520825.0/4194304.0;
      param[3] = 239028075.0/4194304.0;
      param[4] = -1320944625.0/8388608.0;
      param[5] = 2657429775.0/8388608.0;
      param[6] = -1003917915.0/2097152.0;
      param[7] = 1158366825.0/2097152.0;
      param[8] = -4106936925.0/8388608.0;
      param[9] = 2788660875.0/8388608.0;
      param[10] = -717084225.0/4194304.0;
      param[11] = 273795795.0/4194304.0;
      param[12] = -152108775.0/8388608.0;
      param[13] = 35102025.0/8388608.0;
      break;
    case(15):
      param[0] = 5014575.0/33554432.0;
      param[1] = -37702175.0/16777216.0;
      param[2] = 529338537.0/33554432.0;
      param[3] = -575367975.0/8388608.0;
      param[4] = 6931814175.0/33554432.0;
      param[5] = -7661478825.0/16777216.0;
      param[6] = 25688487825.0/33554432.0;
      param[7] = -4159088505.0/4194304.0;
      param[8] = 33592637925.0/33554432.0;
      param[9] = -13233463425.0/16777216.0;
      param[10] = 16174233075.0/33554432.0;
      param[11] = -1890494775.0/8388608.0;
      param[12] = 2646692685.0/33554432.0;
      param[13] = -339319575.0/16777216.0;
      param[14] = 145422675.0/33554432.0;
      break;
    case(16):
      param[0] = -9694845.0/67108864.0;
      param[1] = 155451825.0/67108864.0;
      param[2] = -1168767425.0/67108864.0;
      param[3] = 5469831549.0/67108864.0;
      param[4] = -17836407225.0/67108864.0;
      param[5] = 42977247885.0/67108864.0;
      param[6] = -79168614525.0/67108864.0;
      param[7] = 113763303225.0/67108864.0;
      param[8] = -128931743655.0/67108864.0;
      param[9] = 115707975075.0/67108864.0;
      param[10] = -82047473235.0/67108864.0;
      param[11] = 45581929575.0/67108864.0;
      param[12] = -19535112675.0/67108864.0;
      param[13] = 6311344095.0/67108864.0;
      param[14] = -1502700975.0/67108864.0;
      param[15] = 300540195.0/67108864.0;
      break;
    default:
      return 1;
  }
  return 0;
}




fcs_int ifcs_p2nfft_load_taylor2p_derive_coefficients(
   fcs_int p,
   fcs_float *param
   )
{
  /* Case p<2 is senseless, for p>16 we need to add the coefficients. */
  /* These coefficients are calculated analytically with a Maple script. */
  switch(p){
    case(2):
      param[0] = -1.0;
      break;
    case(3):
      param[0] = 3.0/2.0;
      param[1] = -5.0/2.0;
      break;
    case(4):
      param[0] = -15.0/8.0;
      param[1] = 21.0/4.0;
      param[2] = -35.0/8.0;
      break;
    case(5):
      param[0] = 35.0/16.0;
      param[1] = -135.0/16.0;
      param[2] = 189.0/16.0;
      param[3] = -105.0/16.0;
      break;
    case(6):
      param[0] = -315.0/128.0;
      param[1] = 385.0/32.0;
      param[2] = -1485.0/64.0;
      param[3] = 693.0/32.0;
      param[4] = -1155.0/128.0;
      break;
    case(7):
      param[0] = 693.0/256.0;
      param[1] = -4095.0/256.0;
      param[2] = 5005.0/128.0;
      param[3] = -6435.0/128.0;
      param[4] = 9009.0/256.0;
      param[5] = -3003.0/256.0;
      break;
    case(8):
      param[0] = -3003.0/1024.0;
      param[1] = 10395.0/512.0;
      param[2] = -61425.0/1024.0;
      param[3] = 25025.0/256.0;
      param[4] = -96525.0/1024.0;
      param[5] = 27027.0/512.0;
      param[6] = -15015.0/1024.0;
      break;
    case(9):
      param[0] = 6435.0/2048.0;
      param[1] = -51051.0/2048.0;
      param[2] = 176715.0/2048.0;
      param[3] = -348075.0/2048.0;
      param[4] = 425425.0/2048.0;
      param[5] = -328185.0/2048.0;
      param[6] = 153153.0/2048.0;
      param[7] = -36465.0/2048.0;
      break;
    case(10):
      param[0] = -109395.0/32768.0;
      param[1] = 122265.0/4096.0;
      param[2] = -969969.0/8192.0;
      param[3] = 1119195.0/4096.0;
      param[4] = -6613425.0/16384.0;
      param[5] = 1616615.0/4096.0;
      param[6] = -2078505.0/8192.0;
      param[7] = 415701.0/4096.0;
      param[8] = -692835.0/32768.0;
      break;
    case(11):
      param[0] = 230945.0/65536.0;
      param[1] = -2297295.0/65536.0;
      param[2] = 2567565.0/16384.0;
      param[3] = -6789783.0/16384.0;
      param[4] = 23503095.0/32768.0;
      param[5] = -27776385.0/32768.0;
      param[6] = 11316305.0/16384.0;
      param[7] = -6235515.0/16384.0;
      param[8] = 8729721.0/65536.0;
      param[9] = -1616615.0/65536.0;
      break;
    case(12):
      param[0] = -969969.0/262144.0;
      param[1] = 5311735.0/131072.0;
      param[2] = -52837785.0/262144.0;
      param[3] = 19684665.0/32768.0;
      param[4] = -156165009.0/131072.0;
      param[5] = 108114237.0/65536.0;
      param[6] = -212952285.0/131072.0;
      param[7] = 37182145.0/32768.0;
      param[8] = -143416845.0/262144.0;
      param[9] = 22309287.0/131072.0;
      param[10] = -7436429.0/262144.0;
      break;
    case(13):
      param[0] = 2028117.0/524288.0;
      param[1] = -24249225.0/524288.0;
      param[2] = 132793375.0/524288.0;
      param[3] = -440314875.0/524288.0;
      param[4] = 492116625.0/262144.0;
      param[5] = -780825045.0/262144.0;
      param[6] = 900951975.0/262144.0;
      param[7] = -760543875.0/262144.0;
      param[8] = 929553625.0/524288.0;
      param[9] = -398380125.0/524288.0;
      param[10] = 111546435.0/524288.0;
      param[11] = -16900975.0/524288.0;
      break;
    case(14):
      param[0] = -16900975.0/4194304.0;
      param[1] = 54759159.0/1048576.0;
      param[2] = -654729075.0/2097152.0;
      param[3] = 1195140375.0/1048576.0;
      param[4] = -11888501625.0/4194304.0;
      param[5] = 2657429775.0/524288.0;
      param[6] = -7027425405.0/1048576.0;
      param[7] = 3475100475.0/524288.0;
      param[8] = -20534684625.0/4194304.0;
      param[9] = 2788660875.0/1048576.0;
      param[10] = -2151252675.0/2097152.0;
      param[11] = 273795795.0/1048576.0;
      param[12] = -152108775.0/4194304.0;
      break;
    case(15):
      param[0] = 35102025.0/8388608.0;
      param[1] = -490128275.0/8388608.0;
      param[2] = 1588015611.0/4194304.0;
      param[3] = -6329047725.0/4194304.0;
      param[4] = 34659070875.0/8388608.0;
      param[5] = -68953309425.0/8388608.0;
      param[6] = 25688487825.0/2097152.0;
      param[7] = -29113619535.0/2097152.0;
      param[8] = 100777913775.0/8388608.0;
      param[9] = -66167317125.0/8388608.0;
      param[10] = 16174233075.0/4194304.0;
      param[11] = -5671484325.0/4194304.0;
      param[12] = 2646692685.0/8388608.0;
      param[13] = -339319575.0/8388608.0;
      break;
    case(16):
      param[0] = -145422675.0/33554432.0;
      param[1] = 1088162775.0/16777216.0;
      param[2] = -15193976525.0/33554432.0;
      param[3] = 16409494647.0/8388608.0;
      param[4] = -196200479475.0/33554432.0;
      param[5] = 214886239425.0/16777216.0;
      param[6] = -712517530725.0/33554432.0;
      param[7] = 113763303225.0/4194304.0;
      param[8] = -902522205585.0/33554432.0;
      param[9] = 347123925225.0/16777216.0;
      param[10] = -410237366175.0/33554432.0;
      param[11] = 45581929575.0/8388608.0;
      param[12] = -58605338025.0/33554432.0;
      param[13] = 6311344095.0/16777216.0;
      param[14] = -1502700975.0/33554432.0;
      break;
    default:
      return 1;
  }
  return 0;

}



