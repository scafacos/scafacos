/*
  Copyright (C) 2011, 2012, 2013 Rene Halver, Michael Hofmann

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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/*#include "FCSDefinitions.h"*/
#include "FCSCommon.h"


/** maximal precision */
#ifdef FCS_FLOAT_IS_DOUBLE
#define FCS_FLOAT_PREC 1.0e-14
#else
#define FCS_FLOAT_PREC 1.0e-6
#endif


fcs_int fcs_float_is_equal(fcs_float x, fcs_float y)
{
  return (fcs_fabs(x-y) < FCS_FLOAT_PREC);
}


fcs_int fcs_float_is_zero(fcs_float x)
{
  return (fcs_fabs(x) < FCS_FLOAT_PREC);
}


fcs_int fcs_is_power_of_two(fcs_int x)
{
  return ((x & (x - 1)) == 0);
}


fcs_float fcs_norm(const fcs_float *x)
{
  return fcs_sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}


fcs_int fcs_two_are_orthogonal(const fcs_float *a, const fcs_float *b)
{
  return fcs_float_is_zero(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}


fcs_int fcs_three_are_orthogonal(const fcs_float *a, const fcs_float *b, const fcs_float *c)
{
  return (fcs_two_are_orthogonal(a, b) && fcs_two_are_orthogonal(b, c) && fcs_two_are_orthogonal(c, a));
}


fcs_int fcs_is_orthogonal(const fcs_float *a, const fcs_float *b, const fcs_float *c)
{
  return fcs_three_are_orthogonal(a, b, c);
}


fcs_int fcs_is_cubic(const fcs_float *a, const fcs_float *b, const fcs_float *c)
{
  return fcs_is_orthogonal(a,b,c) && 
    fcs_float_is_equal(fcs_norm(a), fcs_norm(b)) && 
    fcs_float_is_equal(fcs_norm(a), fcs_norm(c));
}


fcs_int fcs_uses_principal_axes(const fcs_float *a, const fcs_float *b, const fcs_float *c)
{
  return
    !fcs_float_is_zero(a[0]) &&  fcs_float_is_zero(a[1]) &&  fcs_float_is_zero(a[2]) &&
     fcs_float_is_zero(b[0]) && !fcs_float_is_zero(b[1]) &&  fcs_float_is_zero(b[2]) &&
     fcs_float_is_zero(c[0]) &&  fcs_float_is_zero(c[1]) && !fcs_float_is_zero(c[2]);
}


static void invert_3x3(const fcs_float *v0, const fcs_float *v1, const fcs_float *v2, fcs_float *iv)
{
  fcs_float det;


  det = v0[0] * v1[1] * v2[2] + v1[0] * v2[1] * v0[2] + v2[0] * v0[1] * v1[2] - v2[0] * v1[1] * v0[2] - v1[0] * v0[1] * v2[2] - v0[0] * v2[1] * v1[2];

  iv[0] = (v1[1] * v2[2] - v2[1] * v1[2]) / det;
  iv[1] = (v2[1] * v0[2] - v0[1] * v2[2]) / det;
  iv[2] = (v0[1] * v1[2] - v1[1] * v0[2]) / det;

  iv[3] = (v2[0] * v1[2] - v1[0] * v2[2]) / det;
  iv[4] = (v0[0] * v2[2] - v2[0] * v0[2]) / det;
  iv[5] = (v1[0] * v0[2] - v0[0] * v1[2]) / det;

  iv[6] = (v1[0] * v2[1] - v2[0] * v1[1]) / det;
  iv[7] = (v2[0] * v0[1] - v0[0] * v2[1]) / det;
  iv[8] = (v0[0] * v1[1] - v1[0] * v0[1]) / det;
}


void fcs_wrap_positions(fcs_int nparticles, fcs_float *positions, const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c, const fcs_float *offset, const fcs_int *periodicity)
{
  fcs_int i;
  fcs_float ibox[9], x[3], y[3];


  if (periodicity == NULL) return;

  invert_3x3(box_a, box_b, box_c, ibox);

  for (i = 0; i < nparticles; ++i)
  {
    x[0] = positions[3 * i + 0] - offset[0];
    x[1] = positions[3 * i + 1] - offset[0];
    x[2] = positions[3 * i + 2] - offset[0];

    y[0] = periodicity[0]?fcs_floor(ibox[0] * x[0] + ibox[3] * x[1] + ibox[6] * x[2]):0;
    y[1] = periodicity[1]?fcs_floor(ibox[1] * x[0] + ibox[4] * x[1] + ibox[7] * x[2]):0;
    y[2] = periodicity[2]?fcs_floor(ibox[2] * x[0] + ibox[5] * x[1] + ibox[8] * x[2]):0;

    positions[3 * i + 0] -= box_a[0] * y[0] + box_b[0] * y[1] + box_c[0] * y[2];
    positions[3 * i + 1] -= box_a[1] * y[0] + box_b[1] * y[1] + box_c[1] * y[2];
    positions[3 * i + 2] -= box_a[2] * y[0] + box_b[2] * y[1] + box_c[2] * y[2];
  }
}


static void fcs_determine_position_extents(fcs_int nparticles, fcs_float *positions, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_float *offset, fcs_float *min, fcs_float *max)
{
  fcs_int i;
  fcs_float ibox[9], x[3], y[3];


  invert_3x3(box_a, box_b, box_c, ibox);

  for (i = 0; i < nparticles; ++i)
  {
    x[0] = positions[3 * i + 0] - offset[0];
    x[1] = positions[3 * i + 1] - offset[1];
    x[2] = positions[3 * i + 2] - offset[2];

    y[0] = ibox[0] * x[0] + ibox[3] * x[1] + ibox[6] * x[2];
    y[1] = ibox[1] * x[0] + ibox[4] * x[1] + ibox[7] * x[2];
    y[2] = ibox[2] * x[0] + ibox[5] * x[1] + ibox[8] * x[2];

    if (i == 0)
    {
      min[0] = max[0] = y[0];
      min[1] = max[1] = y[1];
      min[2] = max[2] = y[2];

    } else
    {
      if (y[0] < min[0]) min[0] = y[0];
      if (y[0] > max[0]) max[0] = y[0];
      if (y[1] < min[1]) min[1] = y[1];
      if (y[1] > max[1]) max[1] = y[1];
      if (y[2] < min[2]) min[2] = y[2];
      if (y[2] > max[2]) max[2] = y[2];
    }
  }
}


void fcs_expand_system_box(fcs_int nparticles, fcs_float *positions, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_float *offset, fcs_int *periodicity)
{
  fcs_float s, min[3], max[3];

  const fcs_float eps1 =
#if defined(FCS_FLOAT_IS_FLOAT)
    1.0f + FLT_EPSILON
#elif defined(FCS_FLOAT_IS_DOUBLE)
    1.0 + DBL_EPSILON
#elif defined(FCS_FLOAT_IS_LONG_DOUBLE)
    1.0l + LDBL_EPSILON
#else
# error FCS float data type is unknown
#endif
    ;


  fcs_determine_position_extents(nparticles, positions, box_a, box_b, box_c, offset, min, max);

/*  printf("min-max: [%f:%f]  [%f:%f]  [%f:%f]\n", min[0], max[0], min[1], max[1], min[2], max[2]);

  printf("old-system: offset: [%.16e,%.16e,%.16e], size: [%.16e,%.16e,%.16e]x[%.16e,%.16e,%.16e]x[%.16e,%.16e,%.16e]\n",
    offset[0], offset[1], offset[2], box_a[0], box_a[1], box_a[2], box_b[0], box_b[1], box_b[2], box_c[0], box_c[1], box_c[2]);*/

  if (periodicity == NULL || !periodicity[0])
  {
    s = 0;
    if (min[0] < 0)
    {
      offset[0] += min[0] * box_a[0];
      offset[1] += min[0] * box_a[1];
      offset[2] += min[0] * box_a[2];
      s -= min[0];
    }
    s += (max[0] >= 1)?max[0]*eps1:1;
    box_a[0] *= s;
    box_a[1] *= s;
    box_a[2] *= s;
  }

  if (periodicity == NULL || !periodicity[1])
  {
    s = 0;
    if (min[1] < 0)
    {
      offset[0] += min[1] * box_b[0];
      offset[1] += min[1] * box_b[1];
      offset[2] += min[1] * box_b[2];
      s -= min[1];
    }
    s += (max[1] >= 1)?max[1]*eps1:1;
    box_b[0] *= s;
    box_b[1] *= s;
    box_b[2] *= s;
  }

  if (periodicity == NULL || !periodicity[2])
  {
    s = 0;
    if (min[2] < 0)
    {
      offset[0] += min[2] * box_c[0];
      offset[1] += min[2] * box_c[1];
      offset[2] += min[2] * box_c[2];
      s -= min[2];
    }
    s += (max[2] >= 1)?max[2]*eps1:1;
    box_c[0] *= s;
    box_c[1] *= s;
    box_c[2] *= s;
  }

/*  printf("new-system: offset: [%.16e,%.16e,%.16e], size: [%.16e,%.16e,%.16e]x[%.16e,%.16e,%.16e]x[%.16e,%.16e,%.16e]\n",
    offset[0], offset[1], offset[2], box_a[0], box_a[1], box_a[2], box_b[0], box_b[1], box_b[2], box_c[0], box_c[1], box_c[2]);*/
}


/******************************/
/* box transforming functions */
/******************************/
void fcs_ftransform_positions(fcs_float *positions, fcs_float *offset, fcs_int local_particles)
{
  int i;
  if (offset != NULL && !(fcs_float_is_zero(offset[0]) && fcs_float_is_zero(offset[1]) && fcs_float_is_zero(offset[2])) )
    for (i = 0; i < local_particles; ++i)
    {
      positions[3*i]   -= offset[0];
      positions[3*i+1] -= offset[1];
      positions[3*i+2] -= offset[2];
    }
}

void fcs_btransform_positions(fcs_float *positions, fcs_float *offset, fcs_int local_particles)
{
  int i;
  if (offset != NULL && !(fcs_float_is_zero(offset[0]) && fcs_float_is_zero(offset[1]) && fcs_float_is_zero(offset[2])) )
    for (i = 0; i < local_particles; ++i)
    {
      positions[3*i]   += offset[0];
      positions[3*i+1] += offset[1];
      positions[3*i+2] += offset[2];
    }
}
