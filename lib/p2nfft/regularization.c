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
#include "part_derive_one_over_norm_x.h"

static fcs_float intpol_even_const(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I);
static fcs_float intpol_even_lin(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I);
static fcs_float intpol_even_quad(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I);
static fcs_float intpol_even_cub(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I);




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
  return sum*pow((xx+1.0),(fcs_float)r)*fcs_pow(1.0-xx,(fcs_float)(m+1))/(1<<(m+1))/fak(r); /* 1<<(m+1) = 2^(m+1) */
}

/** integrated basis polynomial for regularized kernel */
static fcs_float IntBasisPoly(fcs_int p, fcs_int j, fcs_float y)
{
  fcs_int k,l;
  fcs_float sum1=0.0, sum2=0.0;
  
  for (l=0; l<=p; l++) {
    sum1 = 0.0;
    for (k=0; k<=p-j-1; k++) {
      sum1 += binom(p+k-1,k)*fak(j+k)/fak(j+k+1+l)*pow((1.0+y)/2.0,(fcs_float)k);
    }
    sum2 += pow(1.0+y,(fcs_float)l)*pow(1.0-y,(fcs_float)(p-l))/fak(p-l)*sum1;
  }
  return sum2 * fak(p)/fak(j)/(1<<p)*pow(1.0+y,(fcs_float)(j+1)); /* 1<<p = 2^p */
}

static fcs_float IntBasisPoly2(fcs_int p, fcs_int j, fcs_float y)
{
  return (j<0) ? 1.0 : IntBasisPoly(p-1,j,y) - IntBasisPoly(p-1,j,-1);
}

/** regularized kernel for even kernels with K_I even
 *  and K_B symmetric to K(1/2) (used in 1D)
 */
fcs_float ifcs_p2nfft_reg_far_rad_sym(
    ifcs_p2nfft_kernel k, fcs_float xsnorm,
    fcs_int p, const fcs_float *param,
    fcs_float epsI, fcs_float epsB
    )
{
  fcs_float sum=0.0;

  xsnorm = fcs_fabs(xsnorm);

  /* constant continuation for radii > 0.5 */
  if (xsnorm > 0.5)
    xsnorm = 0.5;

  /* regularization at farfield border */
  if ( xsnorm > 0.5-epsB ) {
    fcs_float r = epsB;
    fcs_float m = 0.5;
    fcs_float y = (xsnorm-m)/r;
    for (fcs_int i=0; i<p; i++)
      sum += (BasisPoly(p-1,i,y) + BasisPoly(p-1,i,-y)) * fcs_pow(r,(fcs_float)i) * k(m-r,i,param);
    return sum;
  }
  
  /* nearfield regularization */
  if ( xsnorm < epsI ){
    for (fcs_int i=0; i<p; i++)
      sum += fcs_pow(-epsI,(fcs_float)i) * k(epsI,i,param)
          * (BasisPoly(p-1,i,xsnorm/epsI)+BasisPoly(p-1,i,-xsnorm/epsI));
    return sum;
  }

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
  fcs_int r;
  fcs_float sum=0.0;

  xsnorm = fcs_fabs(xsnorm);

  /* constant continuation for radii > 0.5 */
  if (xsnorm > 0.5)
    xsnorm = 0.5;

  /* regularization at farfield border */
  if ( xsnorm > 0.5-epsB ) {
    sum = BasisPoly(p-1,0,-2.0*xsnorm/epsB+(1.0-epsB)/epsB) * c;
    for (r=0; r<p; r++) {
      sum += fcs_pow(epsB/2.0,(fcs_float)r) * k(0.5-epsB,r,param)
          * BasisPoly(p-1,r,2.0*xsnorm/epsB-(1.0-epsB)/epsB);
    }
    return sum;
  }
  
  /* nearfield regularization */
  if ( xsnorm < epsI ){
    for (r=0; r<p; r++) {
      sum += fcs_pow(-epsI,(fcs_float)r) * k(epsI,r,param)
          * (BasisPoly(p-1,r,xsnorm/epsI)+BasisPoly(p-1,r,-xsnorm/epsI));
    }
    return sum;
  }

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
    fcs_float r_cut, fcs_float eps_B, fcs_float c
    )
{
  fcs_int j;

  xsnorm = fcs_fabs(xsnorm);

  /* constant continuation for radii > 0.5 */
  if (xsnorm > 0.5)
    xsnorm = 0.5;

  /* regulariziton at far field border */
  if (xsnorm > 0.5-eps_B) {
    fcs_float m = (0.5 - eps_B/2.0)/xsnorm;
    fcs_float r = eps_B/(2.0*xsnorm);
    fcs_float y = (1.0-m)/r;
    fcs_float result = BasisPoly(p, 0, -y) * c;
    for (j = 0; j < p; ++j) {
      result += BasisPoly(p, j, y) * pow(r, j) * k((0.5-eps_B)/xsnorm*x2norm, j, param) * pow(x2norm, j);
    }
    return result;
  };

  /* nearfield regularization */
  if (x2norm < r_cut) {
    fcs_float result = 0.0;
    for (j=0; j<p; j++) {
      result+=pow(-r_cut,(fcs_float)j)*k(r_cut,j,param)
        *(BasisPoly(p-1,j,x2norm/r_cut)+BasisPoly(p-1,j,-x2norm/r_cut));
    };
    return result;
  };

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
  fcs_int j;
  fcs_float sum=0.0;

  xsnorm = fcs_fabs(xsnorm);

  /* constant continuation for radii > 0.5 */
  if (xsnorm > 0.5)
    xsnorm = 0.5;

  /* regularization at farfield border */
  if ( xsnorm > 0.5-epsB ) {
    for (j=0; j<=p-2; j++) {
      sum += fcs_pow(epsB/2.0,(fcs_float)j+1) * k(0.5-epsB,j+1,param)
          * (IntBasisPoly(p-1,j,2.0*xsnorm/epsB-(1.0-epsB)/epsB) - IntBasisPoly(p-1,j,-1)
       );
    }
    return sum + k(0.5-epsB,0,param);
  }
 
  /* nearfield regularization */
  if ( xsnorm < epsI ){
    for (j=0; j<p; j++) {
      sum += fcs_pow(-epsI,(fcs_float)j) * k(epsI,j,param)
          * (BasisPoly(p-1,j,xsnorm/epsI)+BasisPoly(p-1,j,-xsnorm/epsI));
    }
    return sum;
  }

  /* farfield: original kernel function */ 
  return k(xsnorm,0,param);
} 



/** regularized kernel for even kernels without singularity, e.g. no K_I needed,
 *  and K_B mirrored smooth into x=1/2 (used in dD, d>1)
 */
fcs_float ifcs_p2nfft_reg_far_no_singularity(
    ifcs_p2nfft_kernel k, fcs_float x2norm, fcs_int p, const fcs_float *param, fcs_float epsB
    )
{
  fcs_int j;
  fcs_float sum=0.0;
  fcs_float h = param[2];

  x2norm = fcs_fabs(x2norm);

  /* inner and outer border of regularization area */
  fcs_float xi = h * (0.5 - epsB);
  fcs_float xo = h * 0.5;

  /* constant continuation for radii > xo */
  if (x2norm > xo)
    x2norm = xo;

  /* canonicalize x to y \in [-1,1] */
  fcs_float r = 0.5 * (xo - xi);
  fcs_float m = 0.5 * (xo + xi);
  fcs_float y = (x2norm-m)/r;

  /* regularization at farfield border */
  if ( xi < x2norm ) {
    for (j=0; j<=p-2; j++) {
      sum += 
        fcs_pow(r,j+1.0) * k(xi,j+1,param)
        * (IntBasisPoly(p-1,j,y) - IntBasisPoly(p-1,j,-1));
    }
    return sum + k(xi,0,param);
  }
 
  /* near- and farfield (no singularity): original kernel function */ 
  return k(x2norm,0,param);
} 


fcs_float ifcs_p2nfft_interpolation(
    fcs_float x, fcs_float one_over_epsI,
    fcs_int order, fcs_int num_nodes,
    const fcs_float *table
    )
{
  switch(order){
    case 0: return intpol_even_const(x, table, num_nodes, one_over_epsI);
    case 1: return intpol_even_lin(x, table, num_nodes, one_over_epsI);
    case 2: return intpol_even_quad(x, table, num_nodes, one_over_epsI);
    default: return intpol_even_cub(x, table, num_nodes, one_over_epsI);
  }
}

/* Regularize the kernel function 1/norm(x) in every dimension i,
 * where xi exceeds (0.5-epsB) * hi, based on one-dimensional symmetric two-point Taylor polynomials, i.e.
 *   d^j/dx^j K(mi-ri) = P^(j)(mi-ri) = (-1)^j P^(j)(mi+ri)
 * We nest these one-dimensional two-point Taylor polynomials in order to construct multi-dimensional ones.
 * The number of variables of P is equal to the number of dimensions with xi > (0.5-epsB) * hi.
 */
fcs_float ifcs_p2nfft_reg_far_rect_sym(
    fcs_float *x, fcs_float *h,
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
    fcs_float *x, fcs_float *h,
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
    fcs_float *x, fcs_float *h,
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
    fcs_float *x, fcs_float *h,
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



fcs_int ifcs_p2nfft_load_taylor2p_coefficients_old(
   fcs_float epsI, fcs_int p,
   fcs_float *param
   )
{
  fcs_float a = epsI;

  /* Case p<2 is senseless, for p>16 we need to add the coefficients. */
  /* These coefficients are calculated analytically with a Maple script. */ 
  switch(p){
    case 2:
      param[0] = -1/(a*a*a)/2.0;
      param[1] = 3.0/2.0/a;
      break;
    case 3:
      param[0] = 3.0/8.0/(a*a*a*a*a);
      param[1] = -5.0/4.0/(a*a*a);
      param[2] = 15.0/8.0/a;
      break;
    case 4:
      param[0] = -5.0/16.0/(a*a*a*a*a*a*a);
      param[1] = 21.0/16.0/(a*a*a*a*a);
      param[2] = -35.0/16.0/(a*a*a);
      param[3] = 35.0/16.0/a;
      break;
    case 5:
      param[0] = 35.0/128.0/(a*a*a*a*a*a*a*a*a);
      param[1] = -45.0/32.0/(a*a*a*a*a*a*a);
      param[2] = 189.0/64.0/(a*a*a*a*a);
      param[3] = -105.0/32.0/(a*a*a);
      param[4] = 315.0/128.0/a;
      break;
    case 6:
      param[0] = -63.0/256.0/pow(a,11.0);
      param[1] = 385.0/256.0/(a*a*a*a*a*a*a*a*a);
      param[2] = -495.0/128.0/(a*a*a*a*a*a*a);
      param[3] = 693.0/128.0/(a*a*a*a*a);
      param[4] = -1155.0/256.0/(a*a*a);
      param[5] = 693.0/256.0/a;
      break;
    case 7:
      param[0] = 231.0/1024.0/pow(a,13.0);
      param[1] = -819.0/512.0/pow(a,11.0);
      param[2] = 5005.0/1024.0/(a*a*a*a*a*a*a*a*a);
      param[3] = -2145.0/256.0/(a*a*a*a*a*a*a);
      param[4] = 9009.0/1024.0/(a*a*a*a*a);
      param[5] = -3003.0/512.0/(a*a*a);
      param[6] = 3003.0/1024.0/a;
      break;
    case 8:
      param[0] = -429.0/2048.0/pow(a,15.0);
      param[1] = 3465.0/2048.0/pow(a,13.0);
      param[2] = -12285.0/2048.0/pow(a,11.0);
      param[3] = 25025.0/2048.0/(a*a*a*a*a*a*a*a*a);
      param[4] = -32175.0/2048.0/(a*a*a*a*a*a*a);
      param[5] = 27027.0/2048.0/(a*a*a*a*a);
      param[6] = -15015.0/2048.0/(a*a*a);
      param[7] = 6435.0/2048.0/a;
      break;
    case 9:
      param[0] = 6435.0/32768.0/pow(a,17.0);
      param[1] = -7293.0/4096.0/pow(a,15.0);
      param[2] = 58905.0/8192.0/pow(a,13.0);
      param[3] = -69615.0/4096.0/pow(a,11.0);
      param[4] = 425425.0/16384.0/(a*a*a*a*a*a*a*a*a);
      param[5] = -109395.0/4096.0/(a*a*a*a*a*a*a);
      param[6] = 153153.0/8192.0/(a*a*a*a*a);
      param[7] = -36465.0/4096.0/(a*a*a);
      param[8] = 109395.0/32768.0/a;
      break;
    case 10:
      param[0] = -12155.0/65536.0/pow(a,19.0);
      param[1] = 122265.0/65536.0/pow(a,17.0);
      param[2] = -138567.0/16384.0/pow(a,15.0);
      param[3] = 373065.0/16384.0/pow(a,13.0);
      param[4] = -1322685.0/32768.0/pow(a,11.0);
      param[5] = 1616615.0/32768.0/(a*a*a*a*a*a*a*a*a);
      param[6] = -692835.0/16384.0/(a*a*a*a*a*a*a);
      param[7] = 415701.0/16384.0/(a*a*a*a*a);
      param[8] = -692835.0/65536.0/(a*a*a);
      param[9] = 230945.0/65536.0/a;
      break;
    case 11:
      param[0] = 46189.0/262144.0/pow(a,21.0);
      param[1] = -255255.0/131072.0/pow(a,19.0);
      param[2] = 2567565.0/262144.0/pow(a,17.0);
      param[3] = -969969.0/32768.0/pow(a,15.0);
      param[4] = 7834365.0/131072.0/pow(a,13.0);
      param[5] = -5555277.0/65536.0/pow(a,11.0);
      param[6] = 11316305.0/131072.0/(a*a*a*a*a*a*a*a*a);
      param[7] = -2078505.0/32768.0/(a*a*a*a*a*a*a);
      param[8] = 8729721.0/262144.0/(a*a*a*a*a);
      param[9] = -1616615.0/131072.0/(a*a*a);
      param[10] = 969969.0/262144.0/a;
      break;
    case 12:
      param[0] = -88179.0/524288.0/pow(a,23.0);
      param[1] = 1062347.0/524288.0/pow(a,21.0);
      param[2] = -5870865.0/524288.0/pow(a,19.0);
      param[3] = 19684665.0/524288.0/pow(a,17.0);
      param[4] = -22309287.0/262144.0/pow(a,15.0);
      param[5] = 36038079.0/262144.0/pow(a,13.0);
      param[6] = -42590457.0/262144.0/pow(a,11.0);
      param[7] = 37182145.0/262144.0/(a*a*a*a*a*a*a*a*a);
      param[8] = -47805615.0/524288.0/(a*a*a*a*a*a*a);
      param[9] = 22309287.0/524288.0/(a*a*a*a*a);
      param[10] = -7436429.0/524288.0/(a*a*a);
      param[11] = 2028117.0/524288.0/a;
      break;
    case 13:
      param[0] = 676039.0/4194304.0/pow(a,25.0);
      param[1] = -2204475.0/1048576.0/pow(a,23.0);
      param[2] = 26558675.0/2097152.0/pow(a,21.0);
      param[3] = -48923875.0/1048576.0/pow(a,19.0);
      param[4] = 492116625.0/4194304.0/pow(a,17.0);
      param[5] = -111546435.0/524288.0/pow(a,15.0);
      param[6] = 300317325.0/1048576.0/pow(a,13.0);
      param[7] = -152108775.0/524288.0/pow(a,11.0);
      param[8] = 929553625.0/4194304.0/(a*a*a*a*a*a*a*a*a);
      param[9] = -132793375.0/1048576.0/(a*a*a*a*a*a*a);
      param[10] = 111546435.0/2097152.0/(a*a*a*a*a);
      param[11] = -16900975.0/1048576.0/(a*a*a);
      param[12] = 16900975.0/4194304.0/a;
      break;
    case 14:
      param[0] = -1300075.0/8388608.0/pow(a,27.0);
      param[1] = 18253053.0/8388608.0/pow(a,25.0);
      param[2] = -59520825.0/4194304.0/pow(a,23.0);
      param[3] = 239028075.0/4194304.0/pow(a,21.0);
      param[4] = -1320944625.0/8388608.0/pow(a,19.0);
      param[5] = 2657429775.0/8388608.0/pow(a,17.0);
      param[6] = -1003917915.0/2097152.0/pow(a,15.0);
      param[7] = 1158366825.0/2097152.0/pow(a,13.0);
      param[8] = -4106936925.0/8388608.0/pow(a,11.0);
      param[9] = 2788660875.0/8388608.0/(a*a*a*a*a*a*a*a*a);
      param[10] = -717084225.0/4194304.0/(a*a*a*a*a*a*a);
      param[11] = 273795795.0/4194304.0/(a*a*a*a*a);
      param[12] = -152108775.0/8388608.0/(a*a*a);
      param[13] = 35102025.0/8388608.0/a;
      break;
    case 15:
      param[0] = 5014575.0/33554432.0/pow(a,29.0);
      param[1] = -37702175.0/16777216.0/pow(a,27.0);
      param[2] = 529338537.0/33554432.0/pow(a,25.0);
      param[3] = -575367975.0/8388608.0/pow(a,23.0);
      param[4] = 6931814175.0/33554432.0/pow(a,21.0);
      param[5] = -7661478825.0/16777216.0/pow(a,19.0);
      param[6] = 25688487825.0/33554432.0/pow(a,17.0);
      param[7] = -4159088505.0/4194304.0/pow(a,15.0);
      param[8] = 33592637925.0/33554432.0/pow(a,13.0);
      param[9] = -13233463425.0/16777216.0/pow(a,11.0);
      param[10] = 16174233075.0/33554432.0/(a*a*a*a*a*a*a*a*a);
      param[11] = -1890494775.0/8388608.0/(a*a*a*a*a*a*a);
      param[12] = 2646692685.0/33554432.0/(a*a*a*a*a);
      param[13] = -339319575.0/16777216.0/(a*a*a);
      param[14] = 145422675.0/33554432.0/a;
      break;
    case 16:
      param[0] = -9694845.0/67108864.0/pow(a,31.0);
      param[1] = 155451825.0/67108864.0/pow(a,29.0);
      param[2] = -1168767425.0/67108864.0/pow(a,27.0);
      param[3] = 5469831549.0/67108864.0/pow(a,25.0);
      param[4] = -17836407225.0/67108864.0/pow(a,23.0);
      param[5] = 42977247885.0/67108864.0/pow(a,21.0);
      param[6] = -79168614525.0/67108864.0/pow(a,19.0);
      param[7] = 113763303225.0/67108864.0/pow(a,17.0);
      param[8] = -128931743655.0/67108864.0/pow(a,15.0);
      param[9] = 115707975075.0/67108864.0/pow(a,13.0);
      param[10] = -82047473235.0/67108864.0/pow(a,11.0);
      param[11] = 45581929575.0/67108864.0/(a*a*a*a*a*a*a*a*a);
      param[12] = -19535112675.0/67108864.0/(a*a*a*a*a*a*a);
      param[13] = 6311344095.0/67108864.0/(a*a*a*a*a);
      param[14] = -1502700975.0/67108864.0/(a*a*a);
      param[15] = 300540195.0/67108864.0/a;
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

fcs_int ifcs_p2nfft_load_taylor2p_derive_coefficients_old(
   fcs_float epsI, fcs_int p,
   fcs_float *param
   )
{
  fcs_float a = epsI;

  /* Case p<2 is senseless, for p>16 we need to add the coefficients. */
  /* These coefficients are calculated analytically with a Maple script. */ 
  switch(p){
    case 2:
      param[0] = -1/(a*a*a);
      break;
    case 3:
      param[0] = 3.0/2.0/(a*a*a*a*a);
      param[1] = -5.0/2.0/(a*a*a);
      break;
    case 4:
      param[0] = -15.0/8.0/(a*a*a*a*a*a*a);
      param[1] = 21.0/4.0/(a*a*a*a*a);
      param[2] = -35.0/8.0/(a*a*a);
      break;
    case 5:
      param[0] = 35.0/16.0/(a*a*a*a*a*a*a*a*a);
      param[1] = -135.0/16.0/(a*a*a*a*a*a*a);
      param[2] = 189.0/16.0/(a*a*a*a*a);
      param[3] = -105.0/16.0/(a*a*a);
      break;
    case 6:
      param[0] = -315.0/128.0/pow(a,11.0);
      param[1] = 385.0/32.0/(a*a*a*a*a*a*a*a*a);
      param[2] = -1485.0/64.0/(a*a*a*a*a*a*a);
      param[3] = 693.0/32.0/(a*a*a*a*a);
      param[4] = -1155.0/128.0/(a*a*a);
      break;
    case 7:
      param[0] = 693.0/256.0/pow(a,13.0);
      param[1] = -4095.0/256.0/pow(a,11.0);
      param[2] = 5005.0/128.0/(a*a*a*a*a*a*a*a*a);
      param[3] = -6435.0/128.0/(a*a*a*a*a*a*a);
      param[4] = 9009.0/256.0/(a*a*a*a*a);
      param[5] = -3003.0/256.0/(a*a*a);
      break;
    case 8:
      param[0] = -3003.0/1024.0/pow(a,15.0);
      param[1] = 10395.0/512.0/pow(a,13.0);
      param[2] = -61425.0/1024.0/pow(a,11.0);
      param[3] = 25025.0/256.0/(a*a*a*a*a*a*a*a*a);
      param[4] = -96525.0/1024.0/(a*a*a*a*a*a*a);
      param[5] = 27027.0/512.0/(a*a*a*a*a);
      param[6] = -15015.0/1024.0/(a*a*a);
      break;
    case 9:
      param[0] = 6435.0/2048.0/pow(a,17.0);
      param[1] = -51051.0/2048.0/pow(a,15.0);
      param[2] = 176715.0/2048.0/pow(a,13.0);
      param[3] = -348075.0/2048.0/pow(a,11.0);
      param[4] = 425425.0/2048.0/(a*a*a*a*a*a*a*a*a);
      param[5] = -328185.0/2048.0/(a*a*a*a*a*a*a);
      param[6] = 153153.0/2048.0/(a*a*a*a*a);
      param[7] = -36465.0/2048.0/(a*a*a);
      break;
    case 10:
      param[0] = -109395.0/32768.0/pow(a,19.0);
      param[1] = 122265.0/4096.0/pow(a,17.0);
      param[2] = -969969.0/8192.0/pow(a,15.0);
      param[3] = 1119195.0/4096.0/pow(a,13.0);
      param[4] = -6613425.0/16384.0/pow(a,11.0);
      param[5] = 1616615.0/4096.0/(a*a*a*a*a*a*a*a*a);
      param[6] = -2078505.0/8192.0/(a*a*a*a*a*a*a);
      param[7] = 415701.0/4096.0/(a*a*a*a*a);
      param[8] = -692835.0/32768.0/(a*a*a);
      break;
    case 11:
      param[0] = 230945.0/65536.0/pow(a,21.0);
      param[1] = -2297295.0/65536.0/pow(a,19.0);
      param[2] = 2567565.0/16384.0/pow(a,17.0);
      param[3] = -6789783.0/16384.0/pow(a,15.0);
      param[4] = 23503095.0/32768.0/pow(a,13.0);
      param[5] = -27776385.0/32768.0/pow(a,11.0);
      param[6] = 11316305.0/16384.0/(a*a*a*a*a*a*a*a*a);
      param[7] = -6235515.0/16384.0/(a*a*a*a*a*a*a);
      param[8] = 8729721.0/65536.0/(a*a*a*a*a);
      param[9] = -1616615.0/65536.0/(a*a*a);
      break;
    case 12:
      param[0] = -969969.0/262144.0/pow(a,23.0);
      param[1] = 5311735.0/131072.0/pow(a,21.0);
      param[2] = -52837785.0/262144.0/pow(a,19.0);
      param[3] = 19684665.0/32768.0/pow(a,17.0);
      param[4] = -156165009.0/131072.0/pow(a,15.0);
      param[5] = 108114237.0/65536.0/pow(a,13.0);
      param[6] = -212952285.0/131072.0/pow(a,11.0);
      param[7] = 37182145.0/32768.0/(a*a*a*a*a*a*a*a*a);
      param[8] = -143416845.0/262144.0/(a*a*a*a*a*a*a);
      param[9] = 22309287.0/131072.0/(a*a*a*a*a);
      param[10] = -7436429.0/262144.0/(a*a*a);
      break;
    case 13:
      param[0] = 2028117.0/524288.0/pow(a,25.0);
      param[1] = -24249225.0/524288.0/pow(a,23.0);
      param[2] = 132793375.0/524288.0/pow(a,21.0);
      param[3] = -440314875.0/524288.0/pow(a,19.0);
      param[4] = 492116625.0/262144.0/pow(a,17.0);
      param[5] = -780825045.0/262144.0/pow(a,15.0);
      param[6] = 900951975.0/262144.0/pow(a,13.0);
      param[7] = -760543875.0/262144.0/pow(a,11.0);
      param[8] = 929553625.0/524288.0/(a*a*a*a*a*a*a*a*a);
      param[9] = -398380125.0/524288.0/(a*a*a*a*a*a*a);
      param[10] = 111546435.0/524288.0/(a*a*a*a*a);
      param[11] = -16900975.0/524288.0/(a*a*a);
      break;
    case 14:
      param[0] = -16900975.0/4194304.0/pow(a,27.0);
      param[1] = 54759159.0/1048576.0/pow(a,25.0);
      param[2] = -654729075.0/2097152.0/pow(a,23.0);
      param[3] = 1195140375.0/1048576.0/pow(a,21.0);
      param[4] = -11888501625.0/4194304.0/pow(a,19.0);
      param[5] = 2657429775.0/524288.0/pow(a,17.0);
      param[6] = -7027425405.0/1048576.0/pow(a,15.0);
      param[7] = 3475100475.0/524288.0/pow(a,13.0);
      param[8] = -20534684625.0/4194304.0/pow(a,11.0);
      param[9] = 2788660875.0/1048576.0/(a*a*a*a*a*a*a*a*a);
      param[10] = -2151252675.0/2097152.0/(a*a*a*a*a*a*a);
      param[11] = 273795795.0/1048576.0/(a*a*a*a*a);
      param[12] = -152108775.0/4194304.0/(a*a*a);
      break;
    case 15:
      param[0] = 35102025.0/8388608.0/pow(a,29.0);
      param[1] = -490128275.0/8388608.0/pow(a,27.0);
      param[2] = 1588015611.0/4194304.0/pow(a,25.0);
      param[3] = -6329047725.0/4194304.0/pow(a,23.0);
      param[4] = 34659070875.0/8388608.0/pow(a,21.0);
      param[5] = -68953309425.0/8388608.0/pow(a,19.0);
      param[6] = 25688487825.0/2097152.0/pow(a,17.0);
      param[7] = -29113619535.0/2097152.0/pow(a,15.0);
      param[8] = 100777913775.0/8388608.0/pow(a,13.0);
      param[9] = -66167317125.0/8388608.0/pow(a,11.0);
      param[10] = 16174233075.0/4194304.0/(a*a*a*a*a*a*a*a*a);
      param[11] = -5671484325.0/4194304.0/(a*a*a*a*a*a*a);
      param[12] = 2646692685.0/8388608.0/(a*a*a*a*a);
      param[13] = -339319575.0/8388608.0/(a*a*a);
      break;
    case 16:
      param[0] = -145422675.0/33554432.0/pow(a,31.0);
      param[1] = 1088162775.0/16777216.0/pow(a,29.0);
      param[2] = -15193976525.0/33554432.0/pow(a,27.0);
      param[3] = 16409494647.0/8388608.0/pow(a,25.0);
      param[4] = -196200479475.0/33554432.0/pow(a,23.0);
      param[5] = 214886239425.0/16777216.0/pow(a,21.0);
      param[6] = -712517530725.0/33554432.0/pow(a,19.0);
      param[7] = 113763303225.0/4194304.0/pow(a,17.0);
      param[8] = -902522205585.0/33554432.0/pow(a,15.0);
      param[9] = 347123925225.0/16777216.0/pow(a,13.0);
      param[10] = -410237366175.0/33554432.0/pow(a,11.0);
      param[11] = 45581929575.0/8388608.0/(a*a*a*a*a*a*a*a*a);
      param[12] = -58605338025.0/33554432.0/(a*a*a*a*a*a*a);
      param[13] = 6311344095.0/16777216.0/(a*a*a*a*a);
      param[14] = -1502700975.0/33554432.0/(a*a*a);
      break;
    default: 
      return 1;
  }

  return 0;
}



/** linear spline interpolation in near field with even kernels */
static fcs_float intpol_even_const(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c;
  fcs_int r;
  fcs_float f1;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  f1=table[r];
  return f1;
}

/** linear spline interpolation in near field with even kernels */
static fcs_float intpol_even_lin(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c,c1,c3;
  fcs_int r;
  fcs_float f1,f2;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  f1=table[r];f2=table[r+1];
  c1=c-r;
  c3=c1-1.0;
  return (-f1*c3+f2*c1);
}

/** quadratic spline interpolation in near field with even kernels */
static fcs_float intpol_even_quad(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c,c1,c2,c3;
  fcs_int r;
  fcs_float f0,f1,f2;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  if (r==0) {f0=table[r+1];f1=table[r];f2=table[r+1];}
  else { f0=table[r-1];f1=table[r];f2=table[r+1];}
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  return (0.5*f0*c1*c3-f1*c2*c3+0.5*f2*c2*c1);
}

/** cubic spline interpolation in near field with even kernels */
static fcs_float intpol_even_cub(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_I
    )
{
  fcs_float c,c1,c2,c3,c4;
  fcs_int r;
  fcs_float f0,f1,f2,f3;

  c=fcs_fabs(x*num_nodes*one_over_eps_I);
  r=(fcs_int)c;
  if (r==0) {f0=table[r+1];f1=table[r];f2=table[r+1];f3=table[r+2];}
  else { f0=table[r-1];f1=table[r];f2=table[r+1];f3=table[r+2];}
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  c4=c1-2.0;
  return(-f0*c1*c3*c4+3.0*f1*c2*c3*c4-3.0*f2*c2*c1*c4+f3*c2*c1*c3)/6.0;
}



