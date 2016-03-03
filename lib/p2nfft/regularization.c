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


/*********************************************************
 * Implementation of different regularization approaches *
 *********************************************************/

/** regularized kernel for even kernels with K_I even
 *  and K_B symmetric to K(1/2) (used in 1D)
 */
fcs_float ifcs_p2nfft_reg_far_rad_sym(
    ifcs_p2nfft_kernel k,  const fcs_float *param,
    fcs_float xsnorm, fcs_int p,
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
fcs_float ifcs_p2nfft_reg_far_rad_ec(
    ifcs_p2nfft_kernel k,  const fcs_float *param,
    fcs_float xsnorm, fcs_int p,
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
fcs_float ifcs_p2nfft_reg_far_rad_ec_noncubic(
    ifcs_p2nfft_kernel k,  const fcs_float *param,
    fcs_float x2norm, fcs_float xsnorm, fcs_int p,
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
fcs_float ifcs_p2nfft_reg_far_rad_ic(
    ifcs_p2nfft_kernel k,  const fcs_float *param,
    fcs_float xsnorm, fcs_int p,
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
    ifcs_p2nfft_kernel k, const fcs_float *param,
    fcs_float x2norm, fcs_int p, fcs_float epsB
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
    ifcs_p2nfft_kernel k, const fcs_float *param,
    fcs_float x2norm, fcs_int p, fcs_float epsB, fcs_float c
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
    ifcs_p2nfft_kernel k, const fcs_float *param,
    fcs_float x2norm, fcs_int p, fcs_float epsB
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
    return ifcs_p2nfft_interpolate_implicit_continuation(k, param, p, xi, xo, x2norm);
 
  /* near- and farfield (no singularity): original kernel function */ 
  return k(x2norm,0,param);
} 


