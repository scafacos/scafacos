/*
  Copyright (C) 2011
  
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

/** \file mmm-common.c
    Common parts of the MMM family of methods for the
    electrostatic interaction, MMM1D, MMM2D.  This file contains the code for the polygamma
    expansions used for the near formulas of MMM1D and MMM2D.

    The expansion of the polygamma functions is fairly easy and follows directly from Abramowitz and Stegun.
    For details, see Axel Arnold and Christian Holm, "MMM2D: A fast and accurate summation method for
    electrostatic interactions in 2D slab geometries", Comp. Phys. Comm., 148/3(2002),327-348.
    */

#include <stdlib.h>
#include <math.h>
#include "mmm-common.h"

#include <stdio.h>

fcs_float mmm_dmax(fcs_float x, fcs_float y) {
  if(x >= y)
    return x;
  else
    return y;
}

fcs_float mmm_dmin(fcs_float x, fcs_float y) {
  if(x <= y)
    return x;
  else
    return y;
}

void mmm_distance2vec(fcs_float x1, fcs_float y1, fcs_float z1, fcs_float x2, fcs_float y2, fcs_float z2, fcs_float vec[3]) {
  vec[0] = x1-x2;
  vec[1] = y1-y2;
  vec[2] = z1-z2;
}

/** Allocate an \ref SizedList of size size. If you need an \ref SizedList
    with variable size better use \ref realloc_doublelist */
void mmm_alloc_doublelist(SizedList *dl, fcs_int size)
{
  dl->max = size;
  dl->e = (fcs_float *) malloc(sizeof(fcs_float)*dl->max);
}

void mmm_realloc_intlist(IntList *il, fcs_int size)
{
  if(size != il->max) {
    il->max = size;
    il->e = (fcs_int *) realloc(il->e, sizeof(fcs_int)*il->max);
  }
}

void mmm_realloc_doublelist(SizedList *dl, fcs_int size)
{
  if(size != dl->max) {
    dl->max = size;
    dl->e = (fcs_float *) realloc(dl->e, sizeof(fcs_float)*dl->max);
  }
}

/** Initialize an \ref DoubleList.  */
void mmm_init_doublelist(SizedList *il)
{
  il->n   = 0;
  il->max = 0;
  il->e   = NULL;
}
