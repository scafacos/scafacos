/*
  Copyright (C) 2011-2012 Rene Halver

  This file is part of ScaFaCoS.

  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.

  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/



#ifndef FCS_DEFINITIONS_INCLUDED
#define FCS_DEFINITIONS_INCLUDED

/**
 * @file FCSInterface.h
 * @brief definitions for the ScaFaCoS library interface
 * @author Rene Halver
 */


#if defined(FCS_FLOAT_IS_FLOAT)
# define FCS_CONST(_c_)  _c_ ## F
#elif defined(FCS_FLOAT_IS_LONG_DOUBLE)
# define FCS_CONST(_c_)  _c_ ## L
#else
# define FCS_CONST(_c_)  _c_
#endif

/**
 * @brief definitions of mathematical constants with FCS namespace (from gcc's math.h)
 */
#define FCS_E         FCS_CONST(2.7182818284590452353602874713526625)  /* e */
#define FCS_LOG2E     FCS_CONST(1.4426950408889634073599246810018921)  /* log_2 e */
#define FCS_LOG10E    FCS_CONST(0.4342944819032518276511289189166051)  /* log_10 e */
#define FCS_LN2       FCS_CONST(0.6931471805599453094172321214581766)  /* log_e 2 */
#define FCS_LN10      FCS_CONST(2.3025850929940456840179914546843642)  /* log_e 10 */
#define FCS_PI        FCS_CONST(3.1415926535897932384626433832795029)  /* pi */
#define FCS_PI_2      FCS_CONST(1.5707963267948966192313216916397514)  /* pi/2 */
#define FCS_PI_4      FCS_CONST(0.7853981633974483096156608458198757)  /* pi/4 */
#define FCS_1_PI      FCS_CONST(0.3183098861837906715377675267450287)  /* 1/pi */
#define FCS_2_PI      FCS_CONST(0.6366197723675813430755350534900574)  /* 2/pi */
#define FCS_2_SQRTPI  FCS_CONST(1.1283791670955125738961589031215452)  /* 2/sqrt(pi) */
#define FCS_SQRT2     FCS_CONST(1.4142135623730950488016887242096981)  /* sqrt(2) */
#define FCS_SQRT1_2   FCS_CONST(0.7071067811865475244008443621048490)  /* 1/sqrt(2) */

/**
 * @brief definitions of mathematical constants with FCS namespace (not included in gcc's math.h)
 */
#define FCS_1_SQRTPI  FCS_CONST(0.5641895835477562869480794515607726)  /* 1/sqrt(pi) */
#define FCS_SQRTPI    FCS_CONST(1.7724538509055160272981674833411452)  /* sqrt(pi) */
#define FCS_PISQR     FCS_CONST(9.8696044010893586188344909998761511)  /* pi^2 */
#define FCS_EULER     FCS_CONST(0.5772156649015328606065120900824024)  /* Euler-Mascheroni constant */

/**
 * @brief definition of a boolean data type
 */
typedef fcs_int fcs_bool;
#define FCS_TRUE           1
#define FCS_FALSE          0
#define FCS_IS_TRUE(_b_)   ((_b_))
#define FCS_IS_FALSE(_b_)  (!(_b_))


/**
 * @brief definitions of return codes
 */
#define FCS_SUCCESS 0
#define FCS_NULL_ARGUMENT 1
#define FCS_ALLOC_FAILED 2
#define FCS_WRONG_ARGUMENT 3
#define FCS_MISSING_ELEMENT 4
#define FCS_LOGICAL_ERROR 5
#define FCS_INCOMPATIBLE_METHOD 6
#define FCS_MPI_ERROR 7
#define FCS_FORTRAN_CALL_ERROR 8
/**
 * @brief definitions of method flags
 */
#define FCS_FMM 32
#define FCS_P2NFFT 33
#define FCS_PEPC 34
#define FCS_P3M 35
#define FCS_PP3MG 36
#define FCS_VMG 37
#define FCS_DIRECT 38
#define FCS_MEMD 39
#define FCS_NO_METHOD_CHOSEN 40
#define FCS_MMM1D 41
#define FCS_EWALD 42
#define FCS_MMM2D 43
#define FCS_WOLF 44

/**
 * @brief definitions of tolerance types
 */
#define FCS_TOLERANCE_TYPE_UNDEFINED      0
#define FCS_TOLERANCE_TYPE_ENERGY         1
#define FCS_TOLERANCE_TYPE_ENERGY_REL     2
#define FCS_TOLERANCE_TYPE_POTENTIAL      3
#define FCS_TOLERANCE_TYPE_POTENTIAL_REL  4
#define FCS_TOLERANCE_TYPE_FIELD          5
#define FCS_TOLERANCE_TYPE_FIELD_REL      6

#endif
