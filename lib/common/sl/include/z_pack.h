/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Michael Hofmann
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


#ifndef __Z_PACK_H__
#define __Z_PACK_H__


#ifdef __cplusplus
# define Z_DECLARE_FUNCTION_BEGIN  extern "C" {
# define Z_DECLARE_FUNCTION_END    }
#else
# define Z_DECLARE_FUNCTION_BEGIN
# define Z_DECLARE_FUNCTION_END
#endif


#ifndef Z_IGNORE_CONFIG_H
# ifdef HAVE_CONFIG_H
#  include <config.h>
# endif
#endif

#ifndef Z_IGNORE_Z_CONFIG_H
# ifdef HAVE_Z_CONFIG_H
#  include "z_config.h"
# endif
#endif


#include "z_pack_conf.h"

#ifdef Z_PACK_RENAME
# include "z_pack_rename.h"
#endif


#define Z_MOP(_mop_)  do { _mop_ } while (0)
#define Z_NOP()       Z_MOP()

#define Z_CONCAT(_a_, _b_)           Z_CONCAT_(_a_, _b_)
#define Z_CONCAT_(_a_, _b_)          _a_##_b_

#define Z_CONCONCAT(_a_, _b_, _c_)   Z_CONCONCAT_(_a_, _b_, _c_)
#define Z_CONCONCAT_(_a_, _b_, _c_)  _a_##_b_##_c_

#define Z_STRINGIFY(_a_)   Z_STRINGIFY_(_a_)
#define Z_STRINGIFY_(_a_)  #_a_


#ifdef Z_PACK_MPI

# include <mpi.h>

# ifdef Z_MPI_COMM_PTR
#  define zcomm_fmt       "p"
#  define zcomm_val(_c_)  (_c_)
# else
#  define zcomm_fmt       "s"
#  define zcomm_val(_c_)  (((_c_) != MPI_COMM_NULL)?"valid":"null")
# endif

Z_DECLARE_FUNCTION_BEGIN
void z_mpi_remap_cart_topology(int from_ndims, int *from_dims, int *from_torus, int *from_pos, int to_ndims, int *to_dims, int *to_torus, int *to_pos);
void z_mpi_get_cart_topology(int *ndims, int *dims, int *torus, int *pos);
void z_mpi_get_grid4d(int *dims, int *pos);

#endif /* Z_PACK_MPI */


#ifdef Z_PACK_NUMERIC

#define z_max(_a_, _b_)           (((_a_)>(_b_))?(_a_):(_b_))
#define z_min(_a_, _b_)           (((_a_)<(_b_))?(_a_):(_b_))
#define z_max3(_a_, _b_, _c_)     z_max(_a_,z_max(_b_,_c_))
#define z_min3(_a_, _b_, _c_)     z_min(_a_,z_min(_b_,_c_))
#define z_minmax(_a_, _b_, _c_)   (((_b_)<(_a_))?(_a_):(((_b_)>(_c_))?(_c_):(_b_)))
#define z_abs(_a_)                (((_a_) >= 0)?(_a_):-(_a_))
#define z_swap(_a_, _b_, _t_)     do { (_t_) = (_a_); (_a_) = (_b_); (_b_) = (_t_); } while (0)

#if HAVE_MATH_H
# ifdef __cplusplus
#  include <cmath>
# else
#  include <math.h>
# endif
#endif

#if defined(HAVE_ROUND) || (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L)
# define z_round(_a_)             round(_a_)
#else
# define z_round(_a_)             (((_a_) >= 0)?floor((_a_) + 0.5):ceil((_a_) - 0.5))
#endif

#define z_powof2_typed(_a_, _t_)  (((_t_) 1) << (_a_))
#define z_powof2(_a_)             z_powof2_typed(_a_, z_int_t)

#define z_ispowof2(_a_)           ((_a_ & (_a_ - 1)) == 0)

#define z_get1d(_x0_)                                       (_x0_)
#define z_get2d(_x1_, _d0_, _x0_)                          ((_x0_) + (_d0_) *  (_x1_))
#define z_get3d(_x2_, _d1_, _x1_, _d0_, _x0_)              ((_x0_) + (_d0_) * ((_x1_) + (_d1_) *  (_x2_)))
#define z_get4d(_x3_, _d2_, _x2_, _d1_, _x1_, _d0_, _x0_)  ((_x0_) + (_d0_) * ((_x1_) + (_d1_) * ((_x2_) + (_d2_) * (_x3_))))

#define z_switch(_x_, _y_, _t_)  ((_t_) = (_x_), (_x_) = (_y_), (_y_) = (_t_))

#define Z_E        2.7182818284590452354   /* e */
#define Z_LOG2E    1.4426950408889634074   /* log_2 e */
#define Z_LOG10E   0.43429448190325182765  /* log_10 e */
#define Z_LN2      0.69314718055994530942  /* log_e 2 */
#define Z_LN10     2.30258509299404568402  /* log_e 10 */
#define Z_PI       3.14159265358979323846  /* pi */
#define Z_PI_2     1.57079632679489661923  /* pi/2 */
#define Z_PI_4     0.78539816339744830962  /* pi/4 */
#define Z_1_PI     0.31830988618379067154  /* 1/pi */
#define Z_2_PI     0.63661977236758134308  /* 2/pi */
#define Z_2_SQRTPI 1.12837916709551257390  /* 2/sqrt(pi) */
#define Z_SQRT2    1.41421356237309504880  /* sqrt(2) */
#define Z_SQRT1_2  0.70710678118654752440  /* 1/sqrt(2) */

#endif /* Z_PACK_NUMERIC */


#ifdef Z_PACK_DEBUG

#if HAVE_STDIO_H || STDC_HEADERS
# ifdef __cplusplus
#  include <cstdio>
# else
#  include <stdio.h>
# endif
#endif

extern FILE *z_notice_fstream, *z_error_fstream, *z_debug_fstream;

#define Z_NOTICE_FSTREAM  (z_notice_fstream?z_notice_fstream:stdout)
#define Z_ERROR_FSTREAM   (z_error_fstream?z_error_fstream:stderr)
#ifdef Z_DEBUG_FSTREAM_STDERR
# define Z_DEBUG_FSTREAM  (z_debug_fstream?z_debug_fstream:stderr)
#else
# define Z_DEBUG_FSTREAM  (z_debug_fstream?z_debug_fstream:stdout)
#endif

#if !defined(Z_DEBUG_MESG_STR) || !defined(Z_DEBUG_MESG_ARG)
# undef Z_DEBUG_MESG_STR
# undef Z_DEBUG_MESG_ARG
# ifdef Z_PACK_MPI
#  define Z_DEBUG_MESG_STR  "%d: "
#  define Z_DEBUG_MESG_ARG  Z_PACK_MPI_RANK
# else
#  define Z_DEBUG_MESG_STR  "%s"
#  define Z_DEBUG_MESG_ARG  ""
# endif
#endif

#if !defined(Z_DEBUG_CODE_STR) || !defined(Z_DEBUG_CODE_ARG)
# undef Z_DEBUG_CODE_STR
# undef Z_DEBUG_CODE_ARG
# define Z_DEBUG_CODE_STR  "%s:%i:%s: "
# define Z_DEBUG_CODE_ARG(_fi_, _li_, _fu_)  _fi_, _li_, _fu_
#endif

#define Z_FPRINTF(_stream_, _format_, ...)                          fprintf(_stream_, _format_ "%s", __VA_ARGS__)
#define Z_FPRINTF_MESG(_stream_, _format_, ...)                     fprintf(_stream_, Z_DEBUG_MESG_STR _format_ "%s", Z_DEBUG_MESG_ARG, __VA_ARGS__)
#define Z_FPRINTF_CODE(_stream_, _format_, ...)                     fprintf(_stream_, Z_DEBUG_MESG_STR Z_DEBUG_CODE_STR _format_ "%s", Z_DEBUG_MESG_ARG, Z_DEBUG_CODE_ARG(__FILE__, __LINE__, __func__), __VA_ARGS__)
#define Z_FPRINTF_CODE_(_fi_, _li_, _fu_, _stream_, _format_, ...)  fprintf(_stream_, Z_DEBUG_MESG_STR Z_DEBUG_CODE_STR _format_ "%s", Z_DEBUG_MESG_ARG, Z_DEBUG_CODE_ARG(_fi_, _li_, _fu_), __VA_ARGS__)

#define Z_NOTICE(...)                              Z_FPRINTF_MESG(Z_NOTICE_FSTREAM, __VA_ARGS__, "\n")
#define Z_NOTICE_IF(_if_, ...)                     Z_MOP(if (_if_) Z_NOTICE(__VA_ARGS__);)
#define Z_NOTICE_(_fi_, _li_, _fu_, ...)           Z_FPRINTF_MESG_(_fi_, _li_, _fu_, Z_NOTICE_FSTREAM, __VA_ARGS__, "\n")
#define Z_NOTICE_IF_(_fi_, _li_, _fu_, _if_, ...)  Z_MOP(if (_if_) Z_NOTICE_(_fi_, _li_, _fu_, __VA_ARGS__);)

#define Z_ERROR(...)                              Z_FPRINTF_MESG(Z_ERROR_FSTREAM, __VA_ARGS__, "\n")
#define Z_ERROR_IF(_if_, ...)                     Z_MOP(if (_if_) Z_ERROR(__VA_ARGS__);)
#define Z_ERROR_(_fi_, _li_, _fu_, ...)           Z_FPRINTF_MESG_(_fi_, _li_, _fu_, Z_ERROR_FSTREAM, __VA_ARGS__, "\n")
#define Z_ERROR_IF_(_fi_, _li_, _fu_, _if_, ...)  Z_MOP(if (_if_) Z_ERROR_(_fi_, _li_, _fu_, __VA_ARGS__);)

#ifdef Z_DEBUG_LEVEL
# define Z_DEBUG_INTRO(_level_, ...)                        Z_MOP(if ((_level_) <= (Z_DEBUG_LEVEL)) Z_FPRINTF_CODE(Z_DEBUG_FSTREAM, __VA_ARGS__, "");)
# define Z_DEBUG_CORE(_level_, ...)                         Z_MOP(if ((_level_) <= (Z_DEBUG_LEVEL)) Z_FPRINTF(Z_DEBUG_FSTREAM, __VA_ARGS__, "");)
# define Z_DEBUG_OUTRO(_level_, ...)                        Z_MOP(if ((_level_) <= (Z_DEBUG_LEVEL)) Z_FPRINTF(Z_DEBUG_FSTREAM, __VA_ARGS__, "\n");)
# define Z_DEBUG(_level_, ...)                              (((_level_) <= (Z_DEBUG_LEVEL))?(Z_FPRINTF_CODE(Z_DEBUG_FSTREAM, __VA_ARGS__, "\n")):0)
# define Z_DEBUG_IF(_if_, _level_, ...)                     (((_if_) && ((_level_) <= (Z_DEBUG_LEVEL)))?(Z_FPRINTF_CODE(Z_DEBUG_FSTREAM, __VA_ARGS__, "\n")):0)
# define Z_DEBUG_(_fi_, _li_, _fu_, _level_, ...)           Z_MOP(if ((_level_) <= (Z_DEBUG_LEVEL)) Z_FPRINTF_CODE_(_fi_, _li_, _fu_, Z_DEBUG_FSTREAM, __VA_ARGS__, "\n");)
# define Z_DEBUG_IF_(_fi_, _li_, _fu_, _if_, _level_, ...)  Z_MOP(if ((_if_) && ((_level_) <= (Z_DEBUG_LEVEL))) Z_FPRINTF_CODE_(_fi_, _li_, _fu_, Z_DEBUG_FSTREAM, __VA_ARGS__, "\n");)
#else
# define Z_DEBUG_INTRO(...)                                 Z_NOP()
# define Z_DEBUG_CORE(...)                                  Z_NOP()
# define Z_DEBUG_OUTRO(...)                                 Z_NOP()
# define Z_DEBUG(...)                                       Z_NOP()
# define Z_DEBUG_IF(...)                                    Z_NOP()
# define Z_DEBUG_(...)                                      Z_NOP()
# define Z_DEBUG_IF_(...)                                   Z_NOP()
#endif

#define Z_TRACE_LEVEL  3

#define Z_TRACE_DECLARE(_c_)                      _c_
#define Z_TRACE_COMMAND(_c_)                      Z_MOP(_c_)

#define Z_TRACE(...)                              Z_DEBUG(Z_TRACE_LEVEL, __VA_ARGS__)
#define Z_TRACE_IF(_if_, ...)                     Z_DEBUG_IF(_if_, Z_TRACE_LEVEL, __VA_ARGS__)
#define Z_TRACE_(_fi_, _li_, _fu_, ...)           Z_DEBUG_(_fi_, _li_, _fu_, Z_TRACE_LEVEL, __VA_ARGS__)
#define Z_TRACE_IF_(_fi_, _li_, _fu_, _if_, ...)  Z_DEBUG_IF_(_fi_, _li_, _fu_, _if_, Z_TRACE_LEVEL, __VA_ARGS__)

#define Z_TRACE_ARRAY(_i_, _n_, _ef_, _e_, ...) \
  Z_MOP(Z_DEBUG_INTRO(Z_TRACE_LEVEL, __VA_ARGS__); \
        for (_i_ = 0; _i_ < (_n_); ++_i_) Z_DEBUG_CORE(Z_TRACE_LEVEL, _ef_, _e_); \
        Z_DEBUG_OUTRO(Z_TRACE_LEVEL, "");)
#define Z_TRACE_ARRAY_IF(_if_, ...)  Z_MOP(if (_if_) Z_TRACE_ARRAY(__VA_ARGS__);)

#define Z_ASSERT_LEVEL  0

#define Z_ASSERT(_x_)           Z_MOP(if (_x_); else Z_DEBUG(Z_ASSERT_LEVEL, "ASSERT: '" #_x_ "' failed."); )
#define Z_ASSERT_IF(_if_, _x_)  Z_MOP(if (_if_) Z_ASSERT(_x_); )

#elif defined(Z_PACK_DEBUG_OFF)

#define Z_NOTICE(...)           Z_NOP()
#define Z_NOTICE_IF(...)        Z_NOP()
#define Z_ERROR(...)            Z_NOP()
#define Z_ERROR_IF(...)         Z_NOP()
#define Z_DEBUG(...)            Z_NOP()
#define Z_DEBUG_IF(...)         Z_NOP()
#define Z_TRACE_DECLARE(_c_)
#define Z_TRACE_COMMAND(_c_)    Z_NOP()
#define Z_TRACE(...)            Z_NOP()
#define Z_TRACE_IF(...)         Z_NOP()
#define Z_TRACE_ARRAY(...)      Z_NOP()
#define Z_TRACE_ARRAY_IF(...)   Z_NOP()
#define Z_ASSERT(_x_)           Z_NOP()
#define Z_ASSERT_IF(_if_, _x_)  Z_NOP()

#endif /* Z_PACK_DEBUG_OFF */


#ifdef Z_PACK_ALLOC

#if HAVE_STDLIB_H || STDC_HEADERS
# ifdef __cplusplus
#  include <cstdlib>
# else
#  include <stdlib.h>
# endif
#endif

#ifndef z_alloc_hook
# define z_alloc_hook_func(_n_, _s_, _r_, _fi_, _li_, _fu_)  (_r_)
#else
Z_DECLARE_FUNCTION_BEGIN
inline static void *z_alloc_hook_func(z_int_t n, z_int_t s, void *r, const char *fi, int li, const char *fu)
{
  z_alloc_hook(n, s, r, fi, li, fu);
  return r;
}
Z_DECLARE_FUNCTION_END
#endif
#ifndef z_realloc_hook_func
# define z_realloc_hook_func(_p_, _n_, _s_, _r_, _fi_, _li_, _fu_)  (_r_)
#else
Z_DECLARE_FUNCTION_BEGIN
inline static void *z_realloc_hook_func(void *p, z_int_t n, z_int_t s, void *r, const char *fi, int li, const char *fu)
{
  z_realloc_hook(p, n, s, r, fi, li, fu);
  return r;
}
Z_DECLARE_FUNCTION_END
#endif
#ifndef z_free_hook
# define z_free_hook(_p_)
#endif

#ifdef Z_ALLOC_DEBUG
# define z_alloc(_n_, _s_)         z_alloc_hook_func((z_int_t) _n_, (z_int_t) _s_, calloc((_n_), (_s_)), __FILE__, __LINE__, __func__)
# define z_realloc(_p_, _n_, _s_)  z_realloc_hook_func(_p_, (z_int_t) _n_, (z_int_t) _s_, realloc(_p_, (_n_) * (_s_)), __FILE__, __LINE__, __func__)
# define z_free(_p_)               Z_MOP(z_free_hook(_p_); free(_p_); _p_ = NULL;)
#else
# define z_alloc(_n_, _s_)         z_alloc_hook_func((z_int_t) _n_, (z_int_t) _s_, malloc((_n_) * (_s_)), __FILE__, __LINE__, __func__)
# define z_realloc(_p_, _n_, _s_)  z_realloc_hook_func(_p_, (z_int_t) _n_, (z_int_t) _s_, realloc(_p_, (_n_) * (_s_)), __FILE__, __LINE__, __func__)
# define z_free(_p_)               Z_MOP(z_free_hook(_p_); free(_p_);)
#endif

#ifndef z_alloca_hook
# define z_alloca_hook_func(_n_, _s_, _p_, _fi_, _li_, _fu_)  (_p_)
#else
Z_DECLARE_FUNCTION_BEGIN
inline static void *z_alloca_hook_func(z_int_t n, z_int_t s, void *p, const char *fi, int li, const char *fu)
{
  z_alloca_hook(n, s, p, fi, li, fu);
  return p;
}
Z_DECLARE_FUNCTION_END
#endif
#ifndef z_freea_hook
# define z_freea_hook(_p_)
#endif

#ifdef HAVE_ALLOCA_H
# include <alloca.h>
#endif

#define z_alloca(_n_, _s_)  z_alloca_hook_func((z_int_t) _n_, (z_int_t) _s_, alloca((_n_) * (_s_)), __FILE__, __LINE__, __func__)
#ifdef Z_ALLOC_DEBUG
# define z_freea(_p_)       Z_MOP(z_freea_hook(_p_); _p_ = NULL;)
#else
# define z_freea(_p_)       Z_MOP(z_freea_hook(_p_);)
#endif

#endif /* Z_PACK_ALLOC */


#if defined(Z_PACK_TIME) || defined(Z_PACK_TIMING)

#ifdef Z_PACK_MPI
typedef double z_time_t;
# define z_time_save(_t_)            (_t_ = MPI_Wtime())
# define z_time_diff_s(_f_, _t_)     ((_t_) - (_f_))
# define z_time_get_s()              (MPI_Wtime())
Z_DECLARE_FUNCTION_BEGIN
inline static double z_time_wtime() { return MPI_Wtime(); }
Z_DECLARE_FUNCTION_END
#else
# if HAVE_STDDEF_H
#  ifdef __cplusplus
#   include <cstddef>
#  else
#   include <stddef.h>
#  endif
# endif
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# endif
typedef struct timeval z_time_t;
# define z_time_save(_t_)            (gettimeofday(&(_t_), NULL))
# define z_time_diff_s(_f_, _t_)     ((double) (((_t_).tv_sec - (_f_).tv_sec) + ((_t_).tv_usec - (_f_).tv_usec) / 1000000.0))
# define z_time_get_s()              (z_time_wtime())
Z_DECLARE_FUNCTION_BEGIN
inline static double z_time_wtime()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double) tv.tv_sec + (tv.tv_usec / 1000000.0);
}
Z_DECLARE_FUNCTION_END
#endif

#endif /* Z_PACK_TIME */


#ifdef Z_PACK_TIMING

#define Z_TIMING
#define Z_TIMING_DECL(_decl_)                     _decl_
#define Z_TIMING_CMD(_cmd_)                       Z_MOP(_cmd_)
#if defined(Z_TIMING_DO_SYNC) && defined(Z_PACK_MPI)
# define Z_TIMING_SYNC(_c_)                       MPI_Barrier(_c_)
#else
# define Z_TIMING_SYNC(_c_)                       Z_NOP()
#endif
#define Z_TIMING_START(_t_)                       ((_t_) = z_time_wtime())
#define Z_TIMING_STOP(_t_)                        ((_t_) = z_time_wtime() - (_t_))
#define Z_TIMING_STOP_ADD(_t_, _r_)               ((_r_) += z_time_wtime() - (_t_))
#ifndef Z_TIMING_PRINT
# define Z_TIMING_PRINT(_i_, _s_, _n_, _v_, _r_)  Z_NOP()
#endif
#ifndef Z_TIMING_PRINT_PREFIX
# define Z_TIMING_PRINT_PREFIX                    "TIMING: "
#endif

#if HAVE_STDIO_H || STDC_HEADERS
# ifdef __cplusplus
#  include <cstdio>
# else
#  include <stdio.h>
# endif
#endif

Z_DECLARE_FUNCTION_BEGIN
inline static void z_timing_print_default(int id, const char *s, z_int_t n, double *v, int rank)
{
  int i;
  printf(Z_TIMING_PRINT_PREFIX "%d: %s:", rank, s);
  for (i = 0; i < n; ++i) printf(" %f", v[i]);
  printf("\n");
}
Z_DECLARE_FUNCTION_END

#else

#define Z_TIMING_DECL(_decl_)
#define Z_TIMING_CMD(_cmd_)                      Z_NOP()
#define Z_TIMING_SYNC(_c_)                       Z_NOP()
#define Z_TIMING_START(_t_)                      Z_NOP()
#define Z_TIMING_STOP(_t_)                       Z_NOP()
#define Z_TIMING_STOP_ADD(_t_, _r_)              Z_NOP()
#undef Z_TIMING_PRINT
#define Z_TIMING_PRINT(_i_, _s_, _n_, _v_, _r_)  Z_NOP()

#endif /* Z_PACK_TIMING */


#ifdef Z_PACK_RANDOM

#ifndef z_rand
# undef z_srand
# undef Z_RAND_MIN
# undef Z_RAND_MAX
# if defined(HAVE_RANDOM) && !defined(__STRICT_ANSI__)
#  if HAVE_STDLIB_H || STDC_HEADERS
#   ifdef __cplusplus
#    include <cstdlib>
#   else
#    include <stdlib.h>
#   endif
#  endif
#  define Z_RAND_MIN    0
#  define Z_RAND_MAX    RAND_MAX
#  define z_srand(_s_)  srandom(_s_)
#  define z_rand()      random()
# elif defined(HAVE_RAND)
#  if HAVE_STDLIB_H || STDC_HEADERS
#   ifdef __cplusplus
#    include <cstdlib>
#   else
#    include <stdlib.h>
#   endif
#  endif
#  define Z_RAND_MIN    0
#  define Z_RAND_MAX    RAND_MAX
#  define z_srand(_s_)  srand(_s_)
#  define z_rand()      rand()
# else
#  define Z_RANDOM_REQUIRED
Z_DECLARE_FUNCTION_BEGIN
void z_srandom(unsigned long seed);
unsigned long z_random();
Z_DECLARE_FUNCTION_END
#  define Z_RAND_MIN    0
#  define Z_RAND_MAX    0xFFFFFFFF
#  define z_srand(_s_)  z_srandom(_s_)
#  define z_rand()      z_random()
# endif
#endif

#ifndef z_rand_minmax
# define z_rand_minmax(_min_, _max_)  ((_min_) + ((double) ((_max_) - (_min_)) * (double) (z_rand() - (Z_RAND_MIN)) / (double) ((Z_RAND_MAX) - (Z_RAND_MIN))))
#endif

#ifndef z_rand_minmax_lt
# define z_rand_minmax_lt(_min_, _max_)  ((_min_) + ((double) ((_max_) - (_min_)) * (double) (z_rand() - (Z_RAND_MIN)) / (double) ((Z_RAND_MAX) - (Z_RAND_MIN) + 1.0)))
#endif

#ifndef z_rand_01
# define z_rand_01()  ((double) (z_rand() - (Z_RAND_MIN)) / (double) ((Z_RAND_MAX) - (Z_RAND_MIN)))
#endif

#ifndef z_rand_01_lt
# define z_rand_01_lt()  ((double) (z_rand() - (Z_RAND_MIN)) / (double) ((Z_RAND_MAX) - (Z_RAND_MIN) + 1.0))
#endif

Z_DECLARE_FUNCTION_BEGIN
void z_srandom64(unsigned long seed);
long long z_random64();
long long z_random64_minmax(long long min, long long max);
unsigned long long z_random64u();
unsigned long long z_random64u_minmax(unsigned long long min, unsigned long long max);

void z_nrandom_seed(unsigned long s);
double z_nrandom();
void z_urandom_seed(unsigned long s);
double z_urandom();
Z_DECLARE_FUNCTION_END

#endif /* Z_PACK_RANDOM */


#ifdef Z_PACK_DIGEST

Z_DECLARE_FUNCTION_BEGIN
z_int_t z_digest_sum_buffer(const void *buffer, z_int_t length, void *sum);
#if HAVE_GCRYPT_H
extern int z_digest_hash_gcrypt_algo;
#endif
void z_digest_hash_open(void **hdl);
void z_digest_hash_close(void *hdl);
void z_digest_hash_write(void *hdl, const void *buffer, z_int_t length);
z_int_t z_digest_hash_read(void *hdl, void *hash);
Z_DECLARE_FUNCTION_END

#endif /* Z_PACK_DIGEST */


#ifdef Z_PACK_CRC32

extern const z_int_t z_crc32_table_size;
extern const z_crc32_t z_crc32_table[];

Z_DECLARE_FUNCTION_BEGIN
void z_crc32_make_table(z_crc32_t *tbl);
void z_crc32_print_table(z_crc32_t *tbl);

z_crc32_t z_crc32_buffer_update(z_crc32_t crc, const void *buffer, z_int_t length);
z_crc32_t z_crc32_buffer(const void *buffer, z_int_t length);
Z_DECLARE_FUNCTION_END

#endif /* Z_CRC32 */


#if defined(Z_PACK_GMP) && HAVE_GMP_H

#if HAVE_GMP_H
# include <gmp.h>
#endif

Z_DECLARE_FUNCTION_BEGIN
void z_gmp_mpz_set_ull(mpz_t z, unsigned long long v);
void z_gmp_mpz_set_sll(mpz_t z, long long v);
unsigned long long z_gmp_mpz_get_ull(mpz_t z);
long long z_gmp_mpz_get_sll(mpz_t z);
Z_DECLARE_FUNCTION_END

#endif /* Z_PACK_GMP && HAVE_GMP_H */


#ifdef Z_PACK_FS

typedef struct
{
  z_int_t is_directory, is_file, is_link, file_size;

} z_fs_stat_t;

Z_DECLARE_FUNCTION_BEGIN
z_int_t z_fs_exists(const char *pathname);
z_int_t z_fs_is_directory(const char *pathname);
z_int_t z_fs_is_file(const char *pathname);
z_int_t z_fs_is_link(const char *pathname);
z_int_t z_fs_get_file_size(const char *pathname);
z_int_t z_fs_get_link_target(const char *pathname, char *target, z_int_t size);
z_int_t z_fs_stat(const char *pathname, z_fs_stat_t *stat);
z_int_t z_fs_mkdir(const char *pathname);
z_int_t z_fs_mkdir_p(const char *pathname);
z_int_t z_fs_rm(const char *pathname);
z_int_t z_fs_rm_r(const char *pathname);
Z_DECLARE_FUNCTION_END

#endif /* Z_PACK_FS */


#ifdef Z_PACK_STDIO

#if HAVE_STDIO_H || STDC_HEADERS
#  ifdef __cplusplus
#   include <cstdio>
#  else
#   include <stdio.h>
#  endif
#endif

int z_snscanf(const char *str, size_t size, const char *format, ...);

#endif /* Z_PACK_STDIO */


#endif /* __Z_PACK_H__ */
