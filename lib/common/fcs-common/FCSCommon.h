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


#ifndef FCS_COMMON_INCLUDED
#define FCS_COMMON_INCLUDED


#ifdef __cplusplus
extern "C" {
#endif


/* fallback definition, see "6.47 Function Names as Strings" in gcc-4.9 doc */
#if __STDC_VERSION__ < 199901L && !defined(__func__)
# if __GNUC__ >= 2
#  define __func__ __FUNCTION__
# else
#  define __func__ "<unknown>"
# endif
#endif


/**
 * @brief macros for mathematical functions corresponding to the FCS float data type
 */
#if defined(FCS_FLOAT_IS_FLOAT)
# define FCS_MATH(_f_)  _f_##f
#elif defined(FCS_FLOAT_IS_DOUBLE)
# define FCS_MATH(_f_)  _f_
#elif defined(FCS_FLOAT_IS_LONG_DOUBLE)
# define FCS_MATH(_f_)  _f_##l
#else
# error FCS float data type is unknown
#endif

#include <math.h>

#define fcs_sqrt(_x_)      FCS_MATH(sqrt)(_x_)
#define fcs_fabs(_x_)      FCS_MATH(fabs)(_x_)
#define fcs_floor(_x_)     FCS_MATH(floor)(_x_)
#define fcs_ceil(_x_)      FCS_MATH(ceil)(_x_)
#define fcs_exp(_x_)       FCS_MATH(exp)(_x_)
#define fcs_sin(_x_)       FCS_MATH(sin)(_x_)
#define fcs_cos(_x_)       FCS_MATH(cos)(_x_)
#define fcs_sinh(_x_)      FCS_MATH(sinh)(_x_)
#define fcs_cosh(_x_)      FCS_MATH(cosh)(_x_)
#define fcs_log(_x_)       FCS_MATH(log)(_x_)
#define fcs_erf(_x_)       FCS_MATH(erf)(_x_)
#define fcs_erfc(_x_)      FCS_MATH(erfc)(_x_)
#define fcs_pow(_x_, _y_)  FCS_MATH(pow)(_x_, _y_)
#define fcs_creal(_x_)     FCS_MATH(creal)(_x_)
#define fcs_cimag(_x_)     FCS_MATH(cimag)(_x_)

#define fcs_isnan(_x_)     isnan(_x_)
#define fcs_isinf(_x_)     isinf(_x_)

#define fcs_xabs(_x_)      (((_x_) >= 0)?(_x_):(-_x_))

/**
 * @brief function to determine if two float values are equal
 * @param x - fcs_float first float value
 * @param y - fcs_float second float value
 * @return fcs_int 1 if x and y are equal, 0 otherwise
 */ 
fcs_int fcs_float_is_equal(fcs_float x, fcs_float y);

/** 
 * @brief function to determine if a float value is zero
 * @param x - fcs_float float value
 * @return fcs_int 1 if x is equal to 0.0, 0 otherwise
 */ 
fcs_int fcs_float_is_zero(fcs_float x);

/** 
 * @brief function to determine if an integer value is a power of two 
 * @param x - fcs_int integer value
 * @return fcs_int 1 if x is a power of two, 0 otherwise
 */ 
fcs_int fcs_is_power_of_two(fcs_int x);

/**
 * @brief function to calculate the Euclidean norm of a given (3D)-vector
 * @param x - fcs_float* vector
 * @return fcs_float Euclidean norm of given vector x
 */
fcs_float fcs_norm(const fcs_float* x);

/**
 * @brief function to determine if two (3D)-vectors are orthogonal
 * @param a - fcs_float* first vector
 * @param b - fcs_float* second vector
 * @return fcs_int 1 if the vectors are orthogonal, 0 otherwise
 */
fcs_int fcs_two_are_orthogonal(const fcs_float *a, const fcs_float *b);

/**
 * @brief function to determine if three (3D)-vectors are mutually orthogonal
 * @param a - fcs_float* first vector
 * @param b - fcs_float* second vector
 * @param c - fcs_float* third vector
 * @return fcs_int 1 if the vectors are mutually orthogonal, 0 otherwise
 */
fcs_int fcs_three_are_orthogonal(const fcs_float *a, const fcs_float *b, const fcs_float *c);

/**
 * @brief function to check if a (3D)-system of base vectors is orthogonal
 * @param a - fcs_float* first base vector
 * @param b - fcs_float* second base vector
 * @param c - fcs_float* third base vector
 * @return fcs_int 1 if the system is orthogonal, 0 otherwise
 */
fcs_int fcs_is_orthogonal(const fcs_float* a, const fcs_float* b, const fcs_float* c);

/**
 * @brief function to check if a (3D)-system of base vectors is cubic
 * @param a - fcs_float* first base vector
 * @param b - fcs_float* second base vector
 * @param c - fcs_float* third base vector
 * @return fcs_int 1 if the system is cubic, 0 otherwise
 */
fcs_int fcs_is_cubic(const fcs_float *a, const fcs_float *b, const fcs_float *c);

/**
 * @brief function to check if the base vectors of a (3D)-system are parallel to the principal axes.
 * @param a - fcs_float* first base vector
 * @param b - fcs_float* second base vector
 * @param c - fcs_float* third base vector
 * @return fcs_int 1 if the system uses the principal axes, 0 otherwise
 */
fcs_int fcs_uses_principal_axes(const fcs_float *a, const fcs_float *b, const fcs_float *c);

/**
 * @brief wrap particle positions according to periodic dimensions
 * @param nparticles fcs_int number of particles
 * @param positions fcs_float* particle positions
 * @param box_a fcs_float* first base vector
 * @param box_b fcs_float* second base vector
 * @param box_c fcs_float* third base vector
 * @param offset fcs_float* offet vector of system box
 * @param periodicity fcs_int* periodic boundaries
 */
void fcs_wrap_positions(fcs_int nparticles, fcs_float *positions, const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c, const fcs_float *offset, const fcs_int *periodicity);

/**
 * @brief expand particle system in open dimensions to enclose the given particles
 * @param nparticles fcs_int number of particles
 * @param positions fcs_float* particle positions
 * @param box_a fcs_float* first base vector
 * @param box_b fcs_float* second base vector
 * @param box_c fcs_float* third base vector
 * @param offset fcs_float* offet vector of system box
 * @param periodicity fcs_int* periodic boundaries
 */
void fcs_expand_system_box(fcs_int nparticles, fcs_float *positions, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_float *offset, fcs_int *periodicity);

/**
 * @brief shift particle positions (i.e. subtracting the offset value)
 * @param nparticles fcs_int number of particles
 * @param positions fcs_float* particle positions
 * @param shift fcs_float* offset of the shift
 */
void fcs_shift_positions(fcs_int nparticles, fcs_float *positions, const fcs_float *offset);

/**
 * @brief undo shift of particle positions (i.e. adding the offset value)
 * @param nparticles fcs_int number of particles
 * @param positions fcs_float* particle positions
 * @param shift fcs_float* offset of the shift
 */
void fcs_unshift_positions(fcs_int nparticles, fcs_float *positions, const fcs_float *offset);


#ifdef __cplusplus
}
#endif


#endif
