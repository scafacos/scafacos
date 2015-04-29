/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS.
 *  
 *  ScaFaCoS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  

 *  
 *  SL - Sorting Library, michael <dot> hofmann <at> informatik <dot> tu-chemnitz <dot> de
 */


#ifndef __ZMPIL_SIMPLE_H__
#define __ZMPIL_SIMPLE_H__


/* quick'n'dirty solution (limited to non-sparse MPI data types at least) */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include <mpi.h>


/*#define ZMPIL_TRACE*/


typedef struct _zmpil_simple_t
{
/*  MPI_Datatype type;*/
  MPI_Aint true_lb, true_extent;

} zmpil_simple_t, *zmpil_simple_p;


int zmpil_simple_create_derived(zmpil_simple_t *mpil, MPI_Datatype type, int count);
void zmpil_simple_destroy(zmpil_simple_t *mpil);

/*void zmpil_simple_memcpy(void *dest, const void *src, size_t n, zmpil_t *mpil);
void zmpil_simple_memmove(void *dest, const void *src, size_t n, zmpil_t *mpil);
void zmpil_simple_memcpy_conv(void *dest, const void *src, size_t n, zmpil_t *zmpil_dest, zmpil_t *zmpil_src);
void zmpil_simple_memmove_conv(void *dest, const void *src, size_t n, zmpil_t *zmpil_dest, zmpil_t *zmpil_src);*/

/* FIXME: extent for offset and true_extent for alloc? */

#define ZMPIL_SIMPLE_TYPE_EXTENT_A(_t_) \
  ((_t_) == MPI_CHAR)?sizeof(char): \
  ((_t_) == MPI_SHORT)?sizeof(signed short int): \
  ((_t_) == MPI_INT)?sizeof(signed int): \
  ((_t_) == MPI_LONG)?sizeof(signed long int): \
  ((_t_) == MPI_LONG_LONG_INT || (_t_) == MPI_LONG_LONG)?sizeof(signed long long int): \
  ((_t_) == MPI_SIGNED_CHAR)?sizeof(signed char): \
  ((_t_) == MPI_UNSIGNED_CHAR)?sizeof(unsigned char): \
  ((_t_) == MPI_UNSIGNED_SHORT)?sizeof(unsigned short int): \
  ((_t_) == MPI_UNSIGNED)?sizeof(unsigned int): \
  ((_t_) == MPI_UNSIGNED_LONG)?sizeof(unsigned long int): \
  ((_t_) == MPI_UNSIGNED_LONG_LONG)?sizeof(unsigned long long int): \
  ((_t_) == MPI_FLOAT)?sizeof(float): \
  ((_t_) == MPI_DOUBLE)?sizeof(double): \
  ((_t_) == MPI_LONG_DOUBLE)?sizeof(long double): \
  ((_t_) == MPI_BYTE)?1:

#define ZMPIL_SIMPLE_TYPE_EXTENT_B(_t_) \
  ((_t_) == MPI_WCHAR)?sizeof(wchar_t):

#define ZMPIL_SIMPLE_TYPE_EXTENT_C(_t_) \
  ((_t_) == MPI_INT8_T)?sizeof(int8_t): \
  ((_t_) == MPI_INT16_T)?sizeof(int16_t): \
  ((_t_) == MPI_INT32_T)?sizeof(int32_t): \
  ((_t_) == MPI_INT64_T)?sizeof(int64_t): \
  ((_t_) == MPI_UINT8_T)?sizeof(uint8_t): \
  ((_t_) == MPI_UINT16_T)?sizeof(uint16_t): \
  ((_t_) == MPI_UINT32_T)?sizeof(uint32_t): \
  ((_t_) == MPI_UINT64_T)?sizeof(uint64_t):

#define ZMPIL_SIMPLE_TYPE_EXTENT_D(_t_) \
  ((_t_) == MPI_C_BOOL)?sizeof(_Bool): \
  ((_t_) == MPI_C_COMPLEX || (_t_) == MPI_C_FLOAT_COMPLEX)?sizeof(float _Complex): \
  ((_t_) == MPI_C_DOUBLE_COMPLEX)?sizeof(double _Complex): \
  ((_t_) == MPI_C_LONG_DOUBLE_COMPLEX)?sizeof(long double _Complex):

#define zmpil_simple_create(_m_, _t_, _n_)  ( \
  (_m_)->true_lb = 0, \
  (_m_)->true_extent = ( \
  ZMPIL_SIMPLE_TYPE_EXTENT_A(_t_) \
  ZMPIL_SIMPLE_TYPE_EXTENT_B(_t_) \
  ZMPIL_SIMPLE_TYPE_EXTENT_C(_t_) \
  ZMPIL_SIMPLE_TYPE_EXTENT_D(_t_) \
  0) * (_n_), \
  ((_m_)->true_extent)?0:zmpil_simple_create_derived((_m_), (_t_), (_n_)) )

#define zmpil_simple_at(_buf_, _n_, _m_)        ((void *) (((char *) _buf_) + ((_n_) * (_m_)->true_extent)))

#define zmpil_simple_diff(_bfrom_, _bto_, _m_)  ((int) ((((char *) _bto_) - ((char *) _bfrom_)) / (_m_)->true_extent))

#define zmpil_simple_extent(_m_)                ((_m_)->true_extent)
#define zmpil_simple_nextent(_n_, _m_)          ((_n_) * (_m_)->true_extent)

#define zmpil_simple_sizefor(_s_, _m_)          ((size_t) (_s_) / (_m_)->true_extent)

#define zmpil_simple_copy(_s_, _d_)             ((_d_)->true_lb = (_s_)->true_lb, (_d_)->true_extent = (_s_)->true_extent)

#ifdef ZMPIL_TRACE
# define zmpil_simple_memcpy(_d_, _s_, _n_, _m_)                               Z_MOP(printf("mpil: copy %d from %p to %p\n", (_n_), (_s_), (_d_)); memcpy(_d_, _s_, (_n_) * (_m_)->true_extent);)
# define zmpil_simple_memmove(_d_, _s_, _n_, _m_)                              Z_MOP(printf("mpil: move %d from %p to %p\n", (_n_), (_s_), (_d_)); memmove(_d_, _s_, (_n_) * (_m_)->true_extent);)
#else
# define zmpil_simple_memcpy(_d_, _s_, _n_, _m_)                               memcpy(_d_, _s_, (_n_) * (_m_)->true_extent)
# define zmpil_simple_memmove(_d_, _s_, _n_, _m_)                              memmove(_d_, _s_, (_n_) * (_m_)->true_extent)
#endif

#define zmpil_simple_memcpy_at(_d_, _dat_, _s_, _sat_, _n_, _m_)               zmpil_simple_memcpy(zmpil_simple_at(_d_, _dat_, _m_), zmpil_simple_at(_s_, _sat_, _m_), _n_, _m_)
#define zmpil_simple_memmove_at(_d_, _dat_, _s_, _sat_, _n_, _m_)              zmpil_simple_memmove(zmpil_simple_at(_d_, _dat_, _m_), zmpil_simple_at(_s_, _sat_, _m_), _n_, _m_)

#define zmpil_simple_memcpy_conv(_d_, _s_, _n_, _md_, _ms_)                    zmpil_simple_memcpy(_d_, _s_, _n_, _md_)
#define zmpil_simple_memmove_conv(_d_, _s_, _n_, _md_, _ms_)                   zmpil_simple_memmove(_d_, _s_, _n_, _md_)

#define zmpil_simple_memcpy_conv_at(_d_, _dat_, _s_, _sat_, _n_, _md_, _ms_)   zmpil_simple_memcpy_conv(zmpil_simple_at(_d_, _dat_, _md_), zmpil_simple_at(_s_, _sat_, _ms_), _n_, _md_, _ms_)
#define zmpil_simple_memmove_conv_at(_d_, _dat_, _s_, _sat_, _n_, _md_, _ms_)  zmpil_simple_memmove_conv(zmpil_simple_at(_d_, _dat_, _md_), zmpil_simple_at(_s_, _sat_, _ms_), _n_, _md_, _ms_)


#endif /* __ZMPIL_SIMPLE_H__ */
