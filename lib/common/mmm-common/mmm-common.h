/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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

/** \file mmm-common.h
    modified polygamma functions. See Arnold,Holm 2002
*/
#ifndef MMM_COMMON_H
#define MMM_COMMON_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

/** \name Math Constants */
/* Mathematical constants, from gcc's math.h */
#ifndef M_PI
#define M_E             2.7182818284590452353602874713526625L  /* e */
#define M_LOG2E         1.4426950408889634073599246810018921L  /* log_2 e */
#define M_LOG10E        0.4342944819032518276511289189166051L  /* log_10 e */
#define M_LN2           0.6931471805599453094172321214581766L  /* log_e 2 */
#define M_LN10  2.3025850929940456840179914546843642L  /* log_e 10 */
#define M_PI            3.1415926535897932384626433832795029L  /* pi */
#define M_PI_2  1.5707963267948966192313216916397514L  /* pi/2 */
#define M_PI_4  0.7853981633974483096156608458198757L  /* pi/4 */
#define M_1_PI  0.3183098861837906715377675267450287L  /* 1/pi */
#define M_2_PI  0.6366197723675813430755350534900574L  /* 2/pi */
#define M_2_SQRTPI      1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#define M_SQRT2         1.4142135623730950488016887242096981L  /* sqrt(2) */
#define M_SQRT1_2       0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif
#define MMM_COMMON_C_2PI     (2*M_PI)
#define MMM_COMMON_C_GAMMA   0.57721566490153286060651209008
#define MMM_COMMON_C_2LOG4PI -5.0620484939385815859557831885
#define MMM_COMMON_C_2PISQR  MMM_COMMON_C_2PI*MMM_COMMON_C_2PI

typedef struct {
  /** Dynamically allocated double field. */
  fcs_float *e;
  /** number of used elements in the double field. */
  fcs_int n;
  /** allocated size of the double field. This value is ONLY changed
      in the routines specified in list operations ! */
  fcs_int max;
} SizedList;

/** Integer list. 
    Use the functions specified in list operations. */
typedef struct {
  /** Dynamically allocated integer field. */
  fcs_int *e;
  /** number of used elements in the integer field. */
  fcs_int n;
  /** allocated size of the integer field. This value is ONLY changed
      in the routines specified in list operations ! */
  fcs_int max;
} IntList;

/** table of the Taylor expansions of the modified polygamma functions */
typedef struct {
  SizedList *modPsi;
  fcs_int   n_modPsi;
} mmm_data_struct;

fcs_float mmm_dmax(fcs_float x, fcs_float y);
fcs_float mmm_dmin(fcs_float x, fcs_float y);
void mmm_distance2vec(fcs_float x1, fcs_float y1, fcs_float z1, fcs_float x2, fcs_float y2, fcs_float z2, fcs_float vec[3]);

/* CUSTOM TYPE FUNCTIONS */
/** Reallocate an \ref IntList */
void mmm_realloc_intlist(IntList *il, fcs_int size);
/** Reallocate an \ref DoubleList */
void mmm_init_doublelist(SizedList *il);
void mmm_alloc_doublelist(SizedList *dl, fcs_int size);
void mmm_realloc_doublelist(SizedList *dl, fcs_int size);

#endif
