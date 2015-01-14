/*
 * Copyright (C) 2011-2013 Michael Pippig
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

#ifndef _P2NFFT_INTERPOLATION_H
#define _P2NFFT_INTERPOLATION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "FCSDefinitions.h"

/** constant spline interpolation in near field with even kernels */
static inline fcs_float ifcs_p2nfft_intpol_even_const(
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

/** constant spline interpolation in far field without symmetry */
static inline fcs_float ifcs_p2nfft_intpol_const(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_B
    )
{
  fcs_float c;
  fcs_int r;
  fcs_float f1;

  c=fcs_fabs(x*num_nodes*one_over_eps_B);
  r=(fcs_int)c;
  f1=table[r+1];
  return f1;
}

/** linear spline interpolation in near field with even kernels */
static inline fcs_float ifcs_p2nfft_intpol_even_lin(
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

/** linear spline interpolation in far field without symmetry */
static inline fcs_float ifcs_p2nfft_intpol_lin(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_B
    )
{
  fcs_float c,c1,c3;
  fcs_int r;
  fcs_float f1,f2;

  c=fcs_fabs(x*num_nodes*one_over_eps_B);
  r=(fcs_int)c;
  f1=table[r+1];f2=table[r+2];
  c1=c-r;
  c3=c1-1.0;
  return (-f1*c3+f2*c1);
}

/** quadratic spline interpolation in near field with even kernels */
static inline fcs_float ifcs_p2nfft_intpol_even_quad(
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

/** quadratic spline interpolation in far field without symmetry */
static inline fcs_float ifcs_p2nfft_intpol_quad(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_B
    )
{
  fcs_float c,c1,c2,c3;
  fcs_int r;
  fcs_float f0,f1,f2;

  c=fcs_fabs(x*num_nodes*one_over_eps_B);
  r=(fcs_int)c;
  f0=table[r];f1=table[r+1];f2=table[r+2];
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  return (0.5*f0*c1*c3-f1*c2*c3+0.5*f2*c2*c1);
}

/** cubic spline interpolation in near field with even kernels */
static inline fcs_float ifcs_p2nfft_intpol_even_cub(
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

/** cubic spline interpolation in far field without symmetry */
static inline fcs_float ifcs_p2nfft_intpol_cub(
    const fcs_float x, const fcs_float *table,
    const fcs_int num_nodes, const fcs_float one_over_eps_B
    )
{
  fcs_float c,c1,c2,c3,c4;
  fcs_int r;
  fcs_float f0,f1,f2,f3;

  c=fcs_fabs(x*num_nodes*one_over_eps_B);
  r=(fcs_int)c;
  f0=table[r];f1=table[r+1];f2=table[r+2];f3=table[r+3];
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  c4=c1-2.0;
  return(-f0*c1*c3*c4+3.0*f1*c2*c3*c4-3.0*f2*c2*c1*c4+f3*c2*c1*c3)/6.0;
}

fcs_float ifcs_p2nfft_interpolation_near(
    fcs_float x, fcs_float one_over_epsI,
    fcs_int order, fcs_int num_nodes,
    const fcs_float *table);
fcs_float ifcs_p2nfft_interpolation_far(
    fcs_float x, fcs_float one_over_epsB,
    fcs_int order, fcs_int num_nodes,
    const fcs_float *table);

#endif

