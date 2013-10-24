/*
  Copyright (C) 2011 Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _P3M_UTILS_H
#define _P3M_UTILS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>
#include <stddef.h>
#include <mpi.h>
#include <stdio.h>

/* Debug macros */
#ifdef FCS_ENABLE_DEBUG
#define ADDITIONAL_CHECKS
#define P3M_TRACE(cmd) if (on_root()) cmd
#define P3M_DEBUG(cmd) if (on_root()) cmd
#define P3M_TRACE_LOCAL(cmd) cmd
#define P3M_DEBUG_LOCAL(cmd) cmd
#define P3M_ENABLE_DEBUG 1
#else
#define P3M_TRACE(cmd)
#define P3M_TRACE_LOCAL(cmd)
#define P3M_DEBUG(cmd) 
#define P3M_DEBUG_LOCAL(cmd) 
#endif

#ifdef FCS_ENABLE_INFO
#define P3M_INFO(cmd) if (on_root()) cmd
#define P3M_INFO_LOCAL(cmd) cmd
#define P3M_ENABLE_INFO 1
#else
#define P3M_INFO(cmd)
#define P3M_INFO_LOCAL(cmd)
#endif

/** maximal precision */
#ifdef FCS_FLOAT_IS_DOUBLE
#define ROUND_ERROR_PREC 1.0e-14
#else
#define ROUND_ERROR_PREC 1.0e-6
#endif

/* Wrap math functions */
#ifdef FCS_FLOAT_IS_FLOAT
#define fabs(x) fabsf(x)
#endif

#ifdef FCS_FLOAT_IS_LONG_DOUBLE
#define fabs(x) fabsl(x)
#endif

/** get the linear index from the position (a,b,c) in a 3D grid
 *  of dimensions adim[]. returns linear index.
 *
 * @return        the linear index
 * @param a       x position 
 * @param b       y position 
 * @param c       z position 
 * @param adim    dimensions of the underlying grid  
 */
static inline int get_linear_index(int a, int b, int c, int adim[3])
{
  return a + adim[0]*(b + adim[1]*c);
}

/** get the position a[] from the linear index in a 3D grid
 *  of dimensions adim[].
 *
 * @param i       linear index
 * @param a       x position (return value) 
 * @param b       y position (return value) 
 * @param c       z position (return value) 
 * @param adim    dimensions of the underlying grid  
 */
static inline void get_grid_pos(int i, int *a, int *b, int *c, int adim[3])
{
  *a = i % adim[0];
  i /= adim[0];
  *b = i % adim[1];
  i /= adim[1];
  *c = i;
}

/** Calculates the maximum of 'int'-typed a and b, returning 'int'. */
static inline int imax(int a, int b) { return (a>b) ? a : b; }

/** Calculates the minimum of 'int'-typed a and b, returning 'int'. */
static inline int imin(int a, int b) { return (a<b) ? a : b; }

/** permute an interger array field of size size about permute positions. */
static inline void permute_ifield(int *field, int size, int permute) {
  int i,tmp;

  if(permute==0) return;
  if(permute<0) permute = (size + permute);
  while(permute>0) {
    tmp=field[0];
    for(i=1;i<size;i++) field[i-1] = field[i];
    field[size-1]=tmp;
    permute--;
  }
}

/** Calculates the SQuaRe of 'fcs_float' x, returning 'fcs_float'. */
static inline fcs_float SQR(fcs_float x) { return x*x; }

/** Calculates the pow(x,3) returning 'fcs_float'. */
static inline fcs_float pow3(fcs_float x) { return x*x*x; }
static inline fcs_int pow3i(fcs_int x) { return x*x*x; }

/** Calculates the sinc-function as sin(PI*x)/(PI*x).
 *
 * (same convention as in Hockney/Eastwood). In order to avoid
 * divisions by 0, arguments, whose modulus is smaller than epsi, will
 * be evaluated by an 8th order Taylor expansion of the sinc
 * function. Note that the difference between sinc(x) and this
 * expansion is smaller than 0.235e-12, if x is smaller than 0.1. (The
 * next term in the expansion is the 10th order contribution
 * PI^10/39916800 * x^10 = 0.2346...*x^12).  This expansion should
 * also save time, since it reduces the number of function calls to
 * sin().  
*/
static inline fcs_float sinc(fcs_float d)
{
  const fcs_float epsi = 0.1;
  const fcs_float c2 = -0.1666666666667e-0;
  const fcs_float c4 = 0.8333333333333e-2;
  const fcs_float c6 = -0.1984126984127e-3;
  const fcs_float c8 = 0.2755731922399e-5;

  fcs_float PId = FCS_PI*d, PId2;

  if (fabs(d)>epsi)
    return sin(PId)/PId;
  else {
    PId2 = SQR(PId);
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
}

static inline int on_root()
{
  static int rank = -1;
  if (rank < 0)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == 0;
}

void errexit();

#endif
