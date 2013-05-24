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

//#define P3M_TIMINGS

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

#if defined(P3M_TIMINGS)
#define FCS_P3M_INIT_TIMING(comm)               			\
  int tm_rank;                                                          \
  MPI_Comm_rank(comm, &tm_rank);                                        \
  double tm_timer, tm_global_timer;
#define FCS_P3M_START_TIMING()                  \
  tm_timer = -MPI_Wtime();
#define FCS_P3M_FINISH_TIMING(comm, str)                                \
  tm_timer += MPI_Wtime();                                              \
  MPI_Reduce(&tm_timer, &tm_global_timer, 1,                            \
             MPI_DOUBLE, MPI_MAX, 0, comm);                             \
  if(!tm_rank) printf("P3M_TIMING: %s takes %e s\n", str, tm_global_timer);
#else
#define FCS_P3M_INIT_TIMING(comm)
#define FCS_P3M_START_TIMING()
#define FCS_P3M_FINISH_TIMING(comm, str)
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


/* Mathematical constants, from gcc's math.h */
#ifndef M_PI
#define M_E             2.7182818284590452353602874713526625L  /* e */
#define M_LOG2E         1.4426950408889634073599246810018921L  /* log_2 e */
#define M_LOG10E        0.4342944819032518276511289189166051L  /* log_10 e */
#define M_LN2           0.6931471805599453094172321214581766L  /* log_e 2 */
#define M_LN10          2.3025850929940456840179914546843642L  /* log_e 10 */
#define M_PI            3.1415926535897932384626433832795029L  /* pi */
#define M_PI_2          1.5707963267948966192313216916397514L  /* pi/2 */
#define M_PI_4          0.7853981633974483096156608458198757L  /* pi/4 */
#define M_1_PI          0.3183098861837906715377675267450287L  /* 1/pi */
#define M_2_PI          0.6366197723675813430755350534900574L  /* 2/pi */
#define M_2_SQRTPI      1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#define M_SQRT2         1.4142135623730950488016887242096981L  /* sqrt(2) */
#define M_SQRT1_2       0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
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

  fcs_float PId = M_PI*d, PId2;

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
