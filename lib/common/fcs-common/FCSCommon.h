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
fcs_float fcs_norm(fcs_float* x);

/**
 * @brief function to determine if two (3D)-vectors are orthogonal
 * @param a - fcs_float* first vector
 * @param b - fcs_float* second vector
 * @return fcs_int 1 if the vectors are orthogonal, 0 otherwise
 */
fcs_int fcs_two_are_orthogonal(fcs_float *a, fcs_float *b);

/**
 * @brief function to determine if three (3D)-vectors are mutually orthogonal
 * @param a - fcs_float* first vector
 * @param b - fcs_float* second vector
 * @param c - fcs_float* third vector
 * @return fcs_int 1 if the vectors are mutually orthogonal, 0 otherwise
 */
fcs_int fcs_three_are_orthogonal(fcs_float *a, fcs_float *b, fcs_float *c);

/**
 * @brief function to check if a (3D)-system of base vectors is orthogonal
 * @param a - fcs_float* first base vector
 * @param b - fcs_float* second base vector
 * @param c - fcs_float* third base vector
 * @return fcs_int 1 if the system is orthogonal, 0 otherwise
 */
fcs_int fcs_is_orthogonal(fcs_float* a, fcs_float* b, fcs_float* c);

/**
 * @brief function to check if a (3D)-system of base vectors is cubic
 * @param a - fcs_float* first base vector
 * @param b - fcs_float* second base vector
 * @param c - fcs_float* third base vector
 * @return fcs_int 1 if the system is cubic, 0 otherwise
 */
fcs_int fcs_is_cubic(fcs_float *a, fcs_float *b, fcs_float *c);

/**
 * @brief function to check if the base vectors of a (3D)-system are parallel to the principal axes.
 * @param a - fcs_float* first base vector
 * @param b - fcs_float* second base vector
 * @param c - fcs_float* third base vector
 * @return fcs_int 1 if the system uses the principal axes, 0 otherwise
 */
fcs_int fcs_uses_principal_axes(fcs_float *a, fcs_float *b, fcs_float *c);

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
void fcs_wrap_positions(fcs_int nparticles, fcs_float *positions, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_float *offset, fcs_int *periodicity);

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
 * @brief transform the positions of particles to a [0:L][0:L][0:L] box
 * @param positions not transformed positions of values_changed
 * @param offset offset particles need to be moved
 * @param local_particles amount of local particles
 */
void fcs_ftransform_positions(fcs_float *positions, fcs_float *offset, fcs_int local_particles);

/**
 * @brief retransform the positions of particles to original box ([-L/2:L/2][-L/2:L/2][-L/2:L/2])
 * @param positions not transformed positions of values_changed
 * @param offset offset particles need to be moved
 * @param local_particles amount of local particles
 */
void fcs_btransform_positions(fcs_float *positions, fcs_float *offset, fcs_int local_particles);


#ifdef __cplusplus
}
#endif


#endif
