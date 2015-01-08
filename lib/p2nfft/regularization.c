/*
 * Copyright (C) 2011-2013 Michael Pippig
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
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "regularization.h"
#include "taylor2p.h"
#include "part_derive_one_over_norm_x.h"


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

static fcs_float IntBasisPoly2(fcs_int p, fcs_int j, fcs_float y)
{
  return (j<0) ? 1.0 : IntBasisPoly(p-1,j,y) - IntBasisPoly(p-1,j,-1);
}

/*********************************************************
 * Implementation of different regularization approaches *
 *********************************************************/

/** regularized kernel for even kernels with K_I even
 *  and K_B symmetric to K(1/2) (used in 1D)
 */
fcs_float ifcs_p2nfft_reg_far_rad_sym(
    ifcs_p2nfft_kernel k, fcs_float xsnorm,
    fcs_int p, const fcs_float *param,
    fcs_float epsI, fcs_float epsB
    )
{
  xsnorm = fcs_fabs(xsnorm);

  /* constant continuation for radii > 0.5 */
  if (xsnorm > 0.5)
    xsnorm = 0.5;

  /* regularization at farfield border */
  if ( xsnorm > 0.5-epsB )
    return ifcs_p2nfft_interpolate_symmetric(k, param, p, 0.5-epsB, 0.5+epsB, xsnorm);
  
  /* nearfield regularization */
  if ( xsnorm < epsI )
    return ifcs_p2nfft_interpolate_symmetric(k, param, p, -epsI, epsI, xsnorm);

  /* farfield: original kernel function */ 
  return k(xsnorm,0,param);
}


/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth to K(1/2) (used in dD, d>1)
 */
fcs_float ifcs_p2nfft_reg_far_rad_expl_cont(
    ifcs_p2nfft_kernel k, fcs_float xsnorm,
    fcs_int p, const fcs_float *param,
    fcs_float epsI, fcs_float epsB, fcs_float c
    )
{
  xsnorm = fcs_fabs(xsnorm);

  /* constant continuation for radii > 0.5 */
  if (xsnorm > 0.5)
    xsnorm = 0.5;

  /* regularization at farfield border */
  if ( xsnorm > 0.5-epsB )
    return ifcs_p2nfft_interpolate_explicit_continuation(k, param, c, p, 0.5-epsB, 0.5, xsnorm);
  
  /* nearfield regularization */
  if ( xsnorm < epsI )
    return ifcs_p2nfft_interpolate_symmetric(k, param, p, -epsI, epsI, xsnorm);

  /* farfield: original kernel function */ 
  return k(xsnorm,0,param);
}



/* regularized kernel for even kernels
 * and in a noncubic box
 * The value xsnorm is assumed to be transformed
 * to [-0.5, 0.5)^3. 
 *
 * Parameters:
 * k: Function handle of the kernel,
 * x2norm: 2-norm of the node,
 * xsnorm: 2-norm of the scaled node into unit cube
 * p: interpolation order,
 * param: parameters for the kernel function,
 * r_cut: inner cutoff radius (without transformation),
 * epsB: outer cutoff radius (with transformation),
 * c: constant continuation value for far field regularization,
 * */
fcs_float ifcs_p2nfft_reg_far_rad_expl_cont_noncubic(
    ifcs_p2nfft_kernel k, fcs_float x2norm, fcs_float xsnorm,
    fcs_int p, const fcs_float *param,
    fcs_float r_cut, fcs_float epsB, fcs_float c
    )
{
  xsnorm = fcs_fabs(xsnorm);

  /* constant continuation for radii > 0.5 */
  if (xsnorm > 0.5)
    xsnorm = 0.5;

  /* regulariziton at far field border */
  if (xsnorm > 0.5-epsB)
    return ifcs_p2nfft_interpolate_explicit_continuation(k, param, c, p, (0.5-epsB)*x2norm/xsnorm, 0.5*x2norm/xsnorm, x2norm);

  /* nearfield regularization */
  if ( x2norm < r_cut )
    return ifcs_p2nfft_interpolate_symmetric(k, param, p, -r_cut, r_cut, x2norm);

  /* farfield: original kernel function */
  return k(x2norm, 0, param);
}

/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth into x=1/2 (used in dD, d>1)
 */
fcs_float ifcs_p2nfft_reg_far_rad_impl_cont(
    ifcs_p2nfft_kernel k, fcs_float xsnorm,
    fcs_int p, const fcs_float *param,
    fcs_float epsI, fcs_float epsB
    )
{
  xsnorm = fcs_fabs(xsnorm);

  /* constant continuation for radii > 0.5 */
  if (xsnorm > 0.5)
    xsnorm = 0.5;

  /* regularization at farfield border */
  if ( xsnorm > 0.5-epsB )
    return ifcs_p2nfft_interpolate_implicit_continuation(k, param, p, 0.5-epsB, 0.5, xsnorm);
 
  /* nearfield regularization */
  if ( xsnorm < epsI )
    return ifcs_p2nfft_interpolate_symmetric(k, param, p, -epsI, epsI, xsnorm);

  /* farfield: original kernel function */ 
  return k(xsnorm,0,param);
} 

/** regularized kernel for even kernels without singularity, i.e. no K_I needed,
 *  and K_B in [1/2-epsB,1/2+epsB) (used in 1D)
 */
fcs_float ifcs_p2nfft_reg_far_rad_sym_no_singularity(
    ifcs_p2nfft_kernel k, fcs_float x2norm, fcs_int p, const fcs_float *param, fcs_float epsB
    )
{
  fcs_float h = param[2];

  x2norm = fcs_fabs(x2norm);

  /* inner and outer border of regularization area */
  fcs_float xi = h * (0.5 - epsB);
  fcs_float xo = h * (0.5 + epsB);

  /* constant continuation for radii > h/2 */
  if (x2norm > xo)
    x2norm = xo;

  /* regularization at farfield border */
  if ( xi < x2norm )
    return ifcs_p2nfft_interpolate_symmetric(k, param, p, xi, xo, x2norm);
 
  /* near- and farfield (no singularity): original kernel function */ 
  return k(x2norm,0,param);
} 

/** regularized kernel for even kernels without singularity, i.e. no K_I needed,
 *  and K_B mirrored smooth into x=1/2 with explicit continuation value 'c' at 1/2 (used in dD, d>1)
 */
fcs_float ifcs_p2nfft_reg_far_rad_ec_no_singularity(
    ifcs_p2nfft_kernel k, fcs_float x2norm, fcs_int p, const fcs_float *param, fcs_float epsB, fcs_float c
    )
{
  fcs_float h = param[2];

  x2norm = fcs_fabs(x2norm);

  /* inner and outer border of regularization area */
  fcs_float xi = h * (0.5 - epsB);
  fcs_float xo = h * 0.5;

  /* constant continuation for radii > xo */
  if (x2norm > xo)
    x2norm = xo;

  /* regularization at farfield border */
  if ( xi < x2norm )
    return ifcs_p2nfft_interpolate_explicit_continuation(k, param, c, p, xi, xo, x2norm);
 
  /* near- and farfield (no singularity): original kernel function */ 
  return k(x2norm,0,param);
} 

/** regularized kernel for even kernels without singularity, i.e. no K_I needed,
 *  and K_B mirrored smooth into x=1/2 with implicit continuation value at 1/2 (used in dD, d>1)
 */
fcs_float ifcs_p2nfft_reg_far_rad_ic_no_singularity(
    ifcs_p2nfft_kernel k, fcs_float x2norm, fcs_int p, const fcs_float *param, fcs_float epsB
    )
{
  fcs_float h = param[2];

  x2norm = fcs_fabs(x2norm);

  /* inner and outer border of regularization area */
  fcs_float xi = h * (0.5 - epsB);
  fcs_float xo = h * 0.5;

  /* constant continuation for radii > xo */
  if (x2norm > xo)
    x2norm = xo;

  /* regularization at farfield border */
  if ( xi < x2norm )
    ifcs_p2nfft_interpolate_implicit_continuation(k, param, p, xi, xo, x2norm);
 
  /* near- and farfield (no singularity): original kernel function */ 
  return k(x2norm,0,param);
} 


/* Regularize the kernel function 1/norm(x) in every dimension i,
 * where xi exceeds (0.5-epsB) * hi, based on one-dimensional symmetric two-point Taylor polynomials, i.e.
 *   d^j/dx^j K(mi-ri) = P^(j)(mi-ri) = (-1)^j P^(j)(mi+ri)
 * We nest these one-dimensional two-point Taylor polynomials in order to construct multi-dimensional ones.
 * The number of variables of P is equal to the number of dimensions with xi > (0.5-epsB) * hi.
 */
fcs_float ifcs_p2nfft_reg_far_rect_sym(
    const fcs_float *x, const fcs_float *h,
    fcs_int p, fcs_float epsB
    )
{
  fcs_int reg_dims[3];
  fcs_float tmp0, tmp1, tmp2, sum=0.0;
  fcs_float m[3], r[3], y[3], z[3];

  /* Due to symmetry we can use fabs(x) and regularize only the positive octant */
  for(int t=0; t<3; t++){
    reg_dims[t] = (fcs_fabs(x[t]) > (0.5-epsB) * h[t]);
    /* Only needed for regularized dimensions */
    if(reg_dims[t]){
      m[t] = 0.5 * h[t];
      r[t] = epsB * h[t];
      y[t] = (fcs_fabs(x[t])-m[t])/r[t];
    }
    /* Evaluation point of the kernel function */
    z[t] = reg_dims[t] ? m[t]-r[t] : fcs_fabs(x[t]);
  }

  for (fcs_int i0=0; i0<p; i0++) {
    tmp0 = reg_dims[0] ? fcs_pow(r[0],(fcs_float)i0) * (BasisPoly(p-1,i0,y[0])+BasisPoly(p-1,i0,-y[0])) : 1.0;
    for (fcs_int i1=0; i1<p; i1++) {
      tmp1 = tmp0 * (reg_dims[1] ? fcs_pow(r[1],(fcs_float)i1) * (BasisPoly(p-1,i1,y[1])+BasisPoly(p-1,i1,-y[1])) : 1.0);
      for (fcs_int i2=0; i2<p; i2++) {
        tmp2 = tmp1 * (reg_dims[2] ? fcs_pow(r[2],(fcs_float)i2) * (BasisPoly(p-1,i2,y[2])+BasisPoly(p-1,i2,-y[2])) : 1.0);
        sum += tmp2 * ifcs_p2nfft_part_derive_one_over_norm_x(i0,i1,i2,z[0],z[1],z[2]);
        if(!reg_dims[2]) break;
      }
      if(!reg_dims[1]) break;
    }
    if(!reg_dims[0]) break;
  }

  return sum;
}

static fcs_float px(
    fcs_float *x, fcs_float *m, fcs_float *r, fcs_int p
    )
{
  fcs_float sum = 0;
  for(fcs_int i=0; i<p; i++){
      fcs_float bp0 = (BasisPoly(p-1,i,(x[0]-m[0])/r[0]) + BasisPoly(p-1,i,-(x[0]-m[0])/r[0])) * fcs_pow(r[0], (fcs_float)i);
      sum += bp0 * ifcs_p2nfft_part_derive_one_over_norm_x(i,0,0,m[0]-r[0],x[1],x[2]);
  }
  return sum;
}

static fcs_float pxy(
    fcs_float *x, fcs_float *m, fcs_float *r, fcs_int p
    )
{
  fcs_float sum = 0;
  for(fcs_int j=0; j<p; j++){
    fcs_float bp1 = (BasisPoly(p-1,j,(x[1]-m[1])/r[1]) + BasisPoly(p-1,j,-(x[1]-m[1])/r[1])) * fcs_pow(r[1],(fcs_float)j);
    for(fcs_int i=0; i<p; i++){
      fcs_float bp0 = (BasisPoly(p-1,i,(x[0]-m[0])/r[0]) + BasisPoly(p-1,i,-(x[0]-m[0])/r[0])) * fcs_pow(r[0], (fcs_float)i);
      sum += bp1 * bp0 * ifcs_p2nfft_part_derive_one_over_norm_x(i,j,0,m[0]-r[1],m[1]-r[1],x[2]);
    }
  }
  return sum;
}

static fcs_float pxyz(
    fcs_float *x, fcs_float *m, fcs_float *r, fcs_int p
    )
{
  fcs_float sum = 0;
  for(fcs_int k=0; k<p; k++){
    fcs_float bp2 = (BasisPoly(p-1,k,(x[2]-m[2])/r[2]) + BasisPoly(p-1,k,-(x[2]-m[2])/r[2])) * fcs_pow(r[2],(fcs_float)k);
    for(fcs_int j=0; j<p; j++){
      fcs_float bp1 = (BasisPoly(p-1,j,(x[1]-m[1])/r[1]) + BasisPoly(p-1,j,-(x[1]-m[1])/r[1])) * fcs_pow(r[1],(fcs_float)j);
      for(fcs_int i=0; i<p; i++){
        fcs_float bp0 = (BasisPoly(p-1,i,(x[0]-m[0])/r[0]) + BasisPoly(p-1,i,-(x[0]-m[0])/r[0])) * fcs_pow(r[0], (fcs_float)i);
        sum += bp2 * bp1 * bp0 * ifcs_p2nfft_part_derive_one_over_norm_x(i,j,k,m[0]-r[0],m[1]-r[1],m[2]-r[2]);
      }
    }
  }
  return sum;
}

static void sort(
    fcs_float *x, fcs_float *y
    )
{
  if(*x < *y){
    fcs_float tmpx = *x; *x = *y; *y = tmpx;
  }
}

/* alternative implementation of reg_far_rect_sym */
fcs_float ifcs_p2nfft_reg_far_rect_sym_version2(
    const fcs_float *x, const fcs_float *h,
    fcs_int p, fcs_float epsB
    )
{
  fcs_float y[3], m[3], r[3];

  for(fcs_int t=0; t<3; t++)
    y[t] = fcs_fabs(x[t]);

  /* sort coordinates in descending order */
  sort(&y[0], &y[1]); sort(&y[1], &y[2]); sort(&y[0], &y[1]);

  for(fcs_int t=0; t<3; t++){
    m[t] = 0.5 * h[t];
    r[t] = epsB * h[t];
  }

  if( y[2] > m[2]-r[2] )
    return pxyz(y, m, r, p);

  if( y[1] > m[1]-r[1] )
    return pxy(y, m, r, p);

  if( y[0] > m[0]-r[0] )
    return px(y, m, r, p);

  return  ifcs_p2nfft_part_derive_one_over_norm_x(0,0,0,x[0],x[1],x[2]);
}



fcs_float ifcs_p2nfft_reg_far_rect_expl_cont(
    const fcs_float *x, const fcs_float *h,
    fcs_int p, fcs_float epsB, fcs_float c
    )
{
  fcs_int reg_dims[3];
  fcs_float tmp0, tmp1, tmp2, sum=0.0;
  fcs_float m[3], r[3], y[3], z[3];

  /* Due to symmetry we can use fabs(x) and regularize only the positive octant */
  for(int t=0; t<3; t++){
    reg_dims[t] = (fcs_fabs(x[t]) > (0.5-epsB) * h[t]);
    /* Only needed for regularized dimensions */
    if(reg_dims[t]){
      m[t] = (0.5 - epsB/2.0) * h[t];
      r[t] = epsB/2.0 * h[t];
      y[t] = (fcs_fabs(x[t])-m[t])/r[t];
    }
    /* Evaluation point of the kernel function */
    z[t] = reg_dims[t] ? m[t]-r[t] : fcs_fabs(x[t]);
  }

  sum = 0;
  for(int t=0; t<3; t++)
    if(reg_dims[t]) sum = sum * BasisPoly(p-1,0,y[t]) + BasisPoly(p-1,0,-y[t]) * c;

  for (fcs_int i0=0; i0<p; i0++) {
    tmp0 = reg_dims[0] ? fcs_pow(r[0],(fcs_float)i0) * BasisPoly(p-1,i0,y[0]) : 1.0;
    for (fcs_int i1=0; i1<p; i1++) {
      tmp1 = tmp0 * (reg_dims[1] ? fcs_pow(r[1],(fcs_float)i1) * BasisPoly(p-1,i1,y[1]) : 1.0);
      for (fcs_int i2=0; i2<p; i2++) {
        tmp2 = tmp1 * (reg_dims[2] ? fcs_pow(r[2],(fcs_float)i2) * BasisPoly(p-1,i2,y[2]) : 1.0);
        sum += tmp2 * ifcs_p2nfft_part_derive_one_over_norm_x(i0,i1,i2,z[0],z[1],z[2]);
        if(!reg_dims[2]) break;
      }
      if(!reg_dims[1]) break;
    }
    if(!reg_dims[0]) break;
  }

  return sum;
}

fcs_float ifcs_p2nfft_reg_far_rect_impl_cont(
    const fcs_float *x, const fcs_float *h,
    fcs_int p, fcs_float epsB
    )
{
  fcs_int reg_dims[3];
  fcs_float tmp0, tmp1, tmp2, sum=0.0;
  fcs_float m[3], r[3], y[3], z[3];

  /* Due to symmetry we can use fabs(x) and regularize only the positive octant */
  for(int t=0; t<3; t++){
    reg_dims[t] = (fcs_fabs(x[t]) > (0.5-epsB) * h[t]);
    /* Only needed for regularized dimensions */
    if(reg_dims[t]){
      m[t] = (0.5 - epsB/2.0) * h[t];
      r[t] = epsB/2.0 * h[t];
      y[t] = (fcs_fabs(x[t])-m[t])/r[t];
    }
    /* Evaluation point of the kernel function */
    z[t] = reg_dims[t] ? m[t]-r[t] : fcs_fabs(x[t]);
  }

  for (fcs_int i0=-1; i0<p-1; i0++) {
    tmp0 = reg_dims[0] ? fcs_pow(r[0],(fcs_float)i0+1) * IntBasisPoly2(p,i0,y[0]) : 1.0;
    for (fcs_int i1=-1; i1<p-1; i1++) {
      tmp1 = tmp0 * (reg_dims[1] ? fcs_pow(r[1],(fcs_float)i1+1) * IntBasisPoly2(p,i1,y[1]) : 1.0);
      for (fcs_int i2=-1; i2<p-1; i2++) {
        tmp2 = tmp1 * (reg_dims[2] ? fcs_pow(r[2],(fcs_float)i2+1) * IntBasisPoly2(p,i2,y[2]) : 1.0);
        sum += tmp2 * ifcs_p2nfft_part_derive_one_over_norm_x(i0+1,i1+1,i2+1,z[0],z[1],z[2]);
        if(!reg_dims[2]) break;
      }
      if(!reg_dims[1]) break;
    }
    if(!reg_dims[0]) break;
  }

  return sum;
}


