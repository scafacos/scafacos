/*
 *  Copyright (C) 2011, 2012, 2013 Michael Hofmann
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


#ifndef __SL_PEPCKEYS_H__
#define __SL_PEPCKEYS_H__

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif /* SL_USE_MPI */

#define SL_PROTO(_f_)  _f_

#include "fortran2c_types.h"


/* enable runtime_informations */
/*#define pepckeys_SL_USE_RTI
#define pepckeys_SL_USE_RTI_TIM*/


/* standard (SL) integer data type */
#define pepckeys_sl_int_type_c          long long
#define pepckeys_sl_int_type_mpi        MPI_LONG_LONG
#define pepckeys_sl_int_size_mpi        1
#define pepckeys_sl_int_type_fmt        "lld"


/* index data type */
#define pepckeys_sl_index_type_c        FINT_TYPE_C
#define pepckeys_sl_index_type_mpi      FINT_TYPE_MPI
#define pepckeys_sl_index_size_mpi      1
#define pepckeys_sl_index_type_fmt      FINT_TYPE_FMT

/* use indices */
#define pepckeys_SL_INDEX


/* keys */
#define pepckeys_sl_key_type_c          FINT8_TYPE_C
#define pepckeys_sl_key_type_mpi        FINT8_TYPE_MPI
#define pepckeys_sl_key_size_mpi        1
#define pepckeys_sl_key_type_fmt        FINT8_TYPE_FMT
#define pepckeys_sl_key_integer

/* data0: work loads */
#define pepckeys_SL_DATA0
#define pepckeys_sl_data0_type_c        FREAL8_TYPE_C
#define pepckeys_sl_data0_size_c        1
#define pepckeys_sl_data0_type_mpi      FREAL8_TYPE_MPI
#define pepckeys_sl_data0_size_mpi      1

/* weighted elements */
#define pepckeys_sl_elem_weight(e, at)  ((e)->data0[at])

#define pepckeys_sl_data0_weight

#define pepckeys_MSEG_BORDER_UPDATE_REDUCTION

#define pepckeys_MSEG_INFO

/*#define pepckeys_MSS_ROOT*/


/* do reduce+bcast instead of allreduce on jugene */
#ifdef JUGENE
# define GLOBAL_REDUCEBCAST_THRESHOLD  0
#endif




#if defined(MSEG_ROOT) && !defined(pepckeys_MSEG_ROOT)
# define pepckeys_MSEG_ROOT  MSEG_ROOT
#endif

#if defined(MSEG_BORDER_UPDATE_REDUCTION) && !defined(pepckeys_MSEG_BORDER_UPDATE_REDUCTION)
# define pepckeys_MSEG_BORDER_UPDATE_REDUCTION  MSEG_BORDER_UPDATE_REDUCTION
#endif

#if defined(MSEG_DISABLE_BEST_CHOICE) && !defined(pepckeys_MSEG_DISABLE_BEST_CHOICE)
# define pepckeys_MSEG_DISABLE_BEST_CHOICE  MSEG_DISABLE_BEST_CHOICE
#endif

#if defined(MSEG_DISABLE_MINMAX) && !defined(pepckeys_MSEG_DISABLE_MINMAX)
# define pepckeys_MSEG_DISABLE_MINMAX  MSEG_DISABLE_MINMAX
#endif

#if defined(MSEG_ENABLE_OPTIMZED_LOWHIGH) && !defined(pepckeys_MSEG_ENABLE_OPTIMZED_LOWHIGH)
# define pepckeys_MSEG_ENABLE_OPTIMZED_LOWHIGH  MSEG_ENABLE_OPTIMZED_LOWHIGH
#endif

#if defined(MSEG_FORWARD_ONLY) && !defined(pepckeys_MSEG_FORWARD_ONLY)
# define pepckeys_MSEG_FORWARD_ONLY  MSEG_FORWARD_ONLY
#endif

#if defined(MSEG_INFO) && !defined(pepckeys_MSEG_INFO)
# define pepckeys_MSEG_INFO  MSEG_INFO
#endif

#if defined(MSEG_TRACE_IF) && !defined(pepckeys_MSEG_TRACE_IF)
# define pepckeys_MSEG_TRACE_IF  MSEG_TRACE_IF
#endif






/* override SL_USE_MPI from sl_config.h */
#ifdef SL_USE_MPI_IGNORE
# undef SL_USE_MPI
#endif

#ifdef SL_USE_MPI_FORCE
# ifndef SL_USE_MPI
#  define SL_USE_MPI
# endif
#endif

/* override SL_USE_OMP from sl_config.h */
#ifdef SL_USE_OMP_IGNORE
# undef SL_USE_OMP
#endif

#ifdef SL_USE_OMP_FORCE
# ifndef SL_USE_OMP
#  define SL_USE_OMP
# endif
#endif


#ifndef pepckeys_SL_INDEX
# undef pepckeys_SL_PACKED_INDEX
#endif


/* if no special datatype for (sl default) integer ... */
#ifndef pepckeys_sl_int_type_c
  /* ... use a default one */
# define pepckeys_sl_int_type_c               long      /* sl_macro */
# undef pepckeys_sl_int_type_mpi
# define pepckeys_sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef pepckeys_sl_int_size_mpi
# define pepckeys_sl_int_size_mpi             1         /* sl_macro */
# undef pepckeys_sl_int_type_fmt
# define pepckeys_sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepckeys_sl_int_type_mpi) || !defined(pepckeys_sl_int_size_mpi)
#   error "pepckeys_sl_int_type_mpi and/or pepckeys_sl_int_size_mpi missing"
#  endif
# endif
# ifndef pepckeys_sl_int_type_fmt
#  error "pepckeys_sl_int_type_fmt macro is missing, using d as default"
#  define pepckeys_sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for (intern) weight ... */
#ifndef pepckeys_sl_weight_type_c
 /* ... use (sl default) integer */
# define pepckeys_sl_weight_type_c             pepckeys_sl_int_type_c    /* sl_macro */
# undef pepckeys_sl_weight_type_mpi
# define pepckeys_sl_weight_type_mpi           pepckeys_sl_int_type_mpi  /* sl_macro */
# undef pepckeys_sl_weight_size_mpi
# define pepckeys_sl_weight_size_mpi           pepckeys_sl_int_size_mpi  /* sl_macro */
# undef pepckeys_sl_weight_type_fmt
# define pepckeys_sl_weight_type_fmt           pepckeys_sl_int_type_fmt  /* sl_macro */
# undef pepckeys_sl_weight_intequiv
# define pepckeys_sl_weight_intequiv                            /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepckeys_sl_weight_type_mpi) || !defined(pepckeys_sl_weight_size_mpi)
#   error "pepckeys_sl_weight_type_mpi and/or pepckeys_sl_weight_size_mpi missing"
#  endif
# endif
# ifndef pepckeys_sl_weight_type_fmt
#  error "pepckeys_sl_weight_type_fmt macro is missing, using f as default"
#  define pepckeys_sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef pepckeys_sl_index_type_c
 /* ... use the primary integer type */
# define pepckeys_sl_index_type_c             pepckeys_sl_int_type_c
# undef pepckeys_sl_index_type_mpi
# define pepckeys_sl_index_type_mpi           pepckeys_sl_int_type_mpi
# undef pepckeys_sl_index_size_mpi
# define pepckeys_sl_index_size_mpi           pepckeys_sl_int_size_mpi
# undef pepckeys_sl_index_type_fmt
# define pepckeys_sl_index_type_fmt           pepckeys_sl_int_type_fmt
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepckeys_sl_index_type_mpi) || !defined(pepckeys_sl_index_size_mpi)
#   error "pepckeys_sl_index_type_mpi and/or pepckeys_sl_index_size_mpi missing"
#  endif
# endif
# ifndef pepckeys_sl_index_type_fmt
#  error "pepckeys_sl_index_type_fmt macro is missing, using d as default"
#  define pepckeys_sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef pepckeys_sl_key_pure_type_c
# define pepckeys_sl_key_pure_type_c          pepckeys_sl_key_type_c  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_type_mpi
# define pepckeys_sl_key_pure_type_mpi        pepckeys_sl_key_type_mpi  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_size_mpi
# define pepckeys_sl_key_pure_size_mpi        pepckeys_sl_key_size_mpi  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_type_fmt
# ifdef pepckeys_sl_key_type_fmt
#  define pepckeys_sl_key_pure_type_fmt       pepckeys_sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef pepckeys_sl_key_purify
 /* key val -> key val */
 #define pepckeys_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define pepckeys_sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_set_pure
 /* key component pointer and key val */
 #define pepckeys_sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef pepckeys_sl_key_pure_cmp_eq
 #define pepckeys_sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_ne
 #define pepckeys_sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_lt
 #define pepckeys_sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_le
 #define pepckeys_sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_gt
 #define pepckeys_sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_ge
 #define pepckeys_sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef pepckeys_sl_key_cmp_eq
 #define pepckeys_sl_key_cmp_eq(k0, k1)       (pepckeys_sl_key_pure_cmp_eq(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_ne
 #define pepckeys_sl_key_cmp_ne(k0, k1)       (pepckeys_sl_key_pure_cmp_ne(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_lt
 #define pepckeys_sl_key_cmp_lt(k0, k1)       (pepckeys_sl_key_pure_cmp_lt(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_le
 #define pepckeys_sl_key_cmp_le(k0, k1)       (pepckeys_sl_key_pure_cmp_le(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_gt
 #define pepckeys_sl_key_cmp_gt(k0, k1)       (pepckeys_sl_key_pure_cmp_gt(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_ge
 #define pepckeys_sl_key_cmp_ge(k0, k1)       (pepckeys_sl_key_pure_cmp_ge(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif


/* default random key */
#ifdef pepckeys_sl_key_integer
# if !defined(pepckeys_sl_key_val_srand) || !defined(pepckeys_sl_key_val_rand) || !defined(pepckeys_sl_key_val_rand_minmax)
#  undef pepckeys_sl_key_val_srand
#  undef pepckeys_sl_key_val_rand
#  undef pepckeys_sl_key_val_rand_minmax
#  define pepckeys_sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define pepckeys_sl_key_val_rand()                     ((pepckeys_sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define pepckeys_sl_key_val_rand_minmax(_min_, _max_)  ((pepckeys_sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef pepckeys_SL_DATA0_IGNORE
# undef pepckeys_SL_DATA0
#endif
#ifdef pepckeys_SL_DATA1_IGNORE
# undef pepckeys_SL_DATA1
#endif
#ifdef pepckeys_SL_DATA2_IGNORE
# undef pepckeys_SL_DATA2
#endif
#ifdef pepckeys_SL_DATA3_IGNORE
# undef pepckeys_SL_DATA3
#endif
#ifdef pepckeys_SL_DATA4_IGNORE
# undef pepckeys_SL_DATA4
#endif
#ifdef pepckeys_SL_DATA5_IGNORE
# undef pepckeys_SL_DATA5
#endif
#ifdef pepckeys_SL_DATA6_IGNORE
# undef pepckeys_SL_DATA6
#endif
#ifdef pepckeys_SL_DATA7_IGNORE
# undef pepckeys_SL_DATA7
#endif
#ifdef pepckeys_SL_DATA8_IGNORE
# undef pepckeys_SL_DATA8
#endif
#ifdef pepckeys_SL_DATA9_IGNORE
# undef pepckeys_SL_DATA9
#endif
#ifdef pepckeys_SL_DATA10_IGNORE
# undef pepckeys_SL_DATA10
#endif
#ifdef pepckeys_SL_DATA11_IGNORE
# undef pepckeys_SL_DATA11
#endif
#ifdef pepckeys_SL_DATA12_IGNORE
# undef pepckeys_SL_DATA12
#endif
#ifdef pepckeys_SL_DATA13_IGNORE
# undef pepckeys_SL_DATA13
#endif
#ifdef pepckeys_SL_DATA14_IGNORE
# undef pepckeys_SL_DATA14
#endif
#ifdef pepckeys_SL_DATA15_IGNORE
# undef pepckeys_SL_DATA15
#endif
#ifdef pepckeys_SL_DATA16_IGNORE
# undef pepckeys_SL_DATA16
#endif
#ifdef pepckeys_SL_DATA17_IGNORE
# undef pepckeys_SL_DATA17
#endif
#ifdef pepckeys_SL_DATA18_IGNORE
# undef pepckeys_SL_DATA18
#endif
#ifdef pepckeys_SL_DATA19_IGNORE
# undef pepckeys_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro pepckeys_sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef pepckeys_sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef pepckeys_sl_data0_weight
# undef pepckeys_sl_data1_weight
# undef pepckeys_sl_data2_weight
# undef pepckeys_sl_data3_weight
# undef pepckeys_sl_data4_weight
# undef pepckeys_sl_data5_weight
# undef pepckeys_sl_data6_weight
# undef pepckeys_sl_data7_weight
# undef pepckeys_sl_data8_weight
# undef pepckeys_sl_data9_weight
# undef pepckeys_sl_data10_weight
# undef pepckeys_sl_data11_weight
# undef pepckeys_sl_data12_weight
# undef pepckeys_sl_data13_weight
# undef pepckeys_sl_data14_weight
# undef pepckeys_sl_data15_weight
# undef pepckeys_sl_data16_weight
# undef pepckeys_sl_data17_weight
# undef pepckeys_sl_data18_weight
# undef pepckeys_sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable pepckeys_sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(pepckeys_sl_data0_weight) && !defined(pepckeys_SL_DATA0)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data1_weight) && !defined(pepckeys_SL_DATA1)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data2_weight) && !defined(pepckeys_SL_DATA2)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data3_weight) && !defined(pepckeys_SL_DATA3)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data4_weight) && !defined(pepckeys_SL_DATA4)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data5_weight) && !defined(pepckeys_SL_DATA5)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data6_weight) && !defined(pepckeys_SL_DATA6)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data7_weight) && !defined(pepckeys_SL_DATA7)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data8_weight) && !defined(pepckeys_SL_DATA8)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data9_weight) && !defined(pepckeys_SL_DATA9)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data10_weight) && !defined(pepckeys_SL_DATA10)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data11_weight) && !defined(pepckeys_SL_DATA11)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data12_weight) && !defined(pepckeys_SL_DATA12)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data13_weight) && !defined(pepckeys_SL_DATA13)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data14_weight) && !defined(pepckeys_SL_DATA14)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data15_weight) && !defined(pepckeys_SL_DATA15)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data16_weight) && !defined(pepckeys_SL_DATA16)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data17_weight) && !defined(pepckeys_SL_DATA17)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data18_weight) && !defined(pepckeys_SL_DATA18)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data19_weight) && !defined(pepckeys_SL_DATA19)
# undef pepckeys_sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro pepckeys_FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef pepckeys_SL_DATA0
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data0_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA1
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data1_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA2
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data2_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA3
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data3_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA4
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data4_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA5
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data5_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA6
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data6_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA7
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data7_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA8
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data8_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA9
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data9_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA10
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data10_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA11
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data11_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA12
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data12_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA13
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data13_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA14
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data14_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA15
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data15_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA16
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data16_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA17
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data17_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA18
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data18_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA19
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data19_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
/* DATAX_TEMPLATE_END */






#ifdef sort_radix_width_default

 #ifndef sort_radix_width_max
  #define sort_radix_width_max       sort_radix_width_default
 #endif

#else /* sort_radix_width_default */

  #ifndef sort_radix_width_max
   #define sort_radix_width_max      12
  #endif

  #define sort_radix_width_default   sort_radix_width_max

#endif /* sort_radix_width_default */

#ifndef sort_radix_db_width_default
# define sort_radix_db_width_default  sort_radix_width_default
#endif

#ifndef sort_radix_ip_width_default
# define sort_radix_ip_width_default  sort_radix_width_default
#endif

#ifndef sort_radix_iter_width_default
# define sort_radix_iter_width_default  sort_radix_width_default
#endif

/* radix-sort thresholds */
#ifndef sort_radix_threshold
# define sort_radix_threshold  256
#endif

#ifndef sort_radix_ip_threshold
# define sort_radix_ip_threshold  sort_radix_threshold
#endif

#ifndef sort_radix_db_threshold
# define sort_radix_db_threshold  sort_radix_threshold
#endif

#ifndef sort_radix_iter_threshold
# define sort_radix_iter_threshold  sort_radix_threshold
#endif


#ifndef ncopy_auto_loop_border_o4o
 #define ncopy_auto_loop_border_o4o  2
#endif

#ifndef ncopy_auto_loop_border_a4o
 #define ncopy_auto_loop_border_a4o  2
#endif

#ifndef nmove_auto_loop_border_o4o
 #define nmove_auto_loop_border_o4o  2
#endif

#ifndef nmove_auto_loop_border_a4o
 #define nmove_auto_loop_border_a4o  2
#endif


/* configure tuneable */
#ifdef SL_TUNEABLE

#if 0
 /* sort_radix_threshold_rec */
 extern int tuneable_sort_radix_threshold_rec;

 /* sort_radix_threshold_iter */
 extern int tuneable_sort_radix_threshold_iter;

#endif
#endif






#ifndef SL_RTI_TID_DECLARED
# define SL_RTI_TID_DECLARED

enum rti_tid
{
  /* src/base_mpi/base_mpi.c */
  rti_tid_mpi_merge2,
  rti_tid_mpi_merge2_find,
  rti_tid_mpi_merge2_moveright,
  rti_tid_mpi_merge2_exchange,
  rti_tid_mpi_merge2_moveleft,
  rti_tid_mpi_merge2_local,
  rti_tid_mpi_mergek_equal,
  rti_tid_mpi_mergek_equal_while,
  rti_tid_mpi_mergek_equal_while_merge2,
  rti_tid_mpi_mergek_sorted,
  rti_tid_mpi_mergek_sorted_while,
  rti_tid_mpi_mergek_sorted_while_check,
  rti_tid_mpi_mergek_sorted_while_oddeven,
  rti_tid_mpi_mergek_sorted2,
  rti_tid_mpi_mergek_sorted2_while,
  rti_tid_mpi_mergek_sorted2_while_check,
  rti_tid_mpi_mergek_sorted2_while_oddeven,
  rti_tid_mpi_mergek,
  rti_tid_mpi_mergek_equalike,
  rti_tid_mpi_mergek_while,
  rti_tid_mpi_mergek_while_check,
  rti_tid_mpi_mergek_while_oddeven,
  rti_tid_mpi_partition_exact_generic,
  rti_tid_mpi_partition_exact_generic_select,
  rti_tid_mpi_partition_exact_generic_rcounts,
  rti_tid_mpi_partition_exact_radix_ngroups,
  rti_tid_mpi_partition_exact_radix_ngroups_pconds,
  rti_tid_mpi_partition_exact_radix_ngroups_idxin,
  rti_tid_mpi_partition_exact_radix_ngroups_up,
  rti_tid_mpi_partition_exact_radix_ngroups_down,
  rti_tid_mpi_partition_exact_radix_ngroups_down_select,
  rti_tid_mpi_partition_exact_radix_ngroups_down_alltoall,
  rti_tid_mpi_partition_exact_radix_ngroups_down_x2suby,
  rti_tid_mpi_partition_exact_radix_ngroups_down_merge,
  rti_tid_mpi_partition_exact_radix_ngroups_idxout,
  rti_tid_mpi_partition_exact_radix_ngroups_idxout_loop,
  rti_tid_mpi_partition_exact_radix_ngroups_idxout_alltoall,
  rti_tid_mpi_partition_exact_radix_2groups,
  rti_tid_mpi_partition_exact_radix_2groups_pconds,
  rti_tid_mpi_partition_exact_radix_2groups_select1st,
  rti_tid_mpi_partition_exact_radix_2groups_x2suby,
  rti_tid_mpi_partition_exact_radix_2groups_alltoall,
  rti_tid_mpi_partition_exact_radix_2groups_select2nd,
  rti_tid_mpi_partition_exact_radix_2groups_subx2y,
  rti_tid_mpi_partition_sample,
  rti_tid_mpi_partition_sample_select,
  rti_tid_mpi_partition_sample_rcounts,
  rti_tid_mpi_select_exact_generic,
  rti_tid_mpi_select_exact_generic_sync_init,
  rti_tid_mpi_select_exact_generic_sync_exit,
  rti_tid_mpi_select_exact_generic_while,
  rti_tid_mpi_select_exact_generic_while_check,
  rti_tid_mpi_select_exact_generic_while_check_bins,
  rti_tid_mpi_select_exact_generic_while_check_bins_local,
  rti_tid_mpi_select_exact_generic_while_check_bins_global,
  rti_tid_mpi_select_exact_generic_while_check_round1,
  rti_tid_mpi_select_exact_generic_while_check_pre,
  rti_tid_mpi_select_exact_generic_while_check_part,
  rti_tid_mpi_select_exact_generic_while_check_part_root,
  rti_tid_mpi_select_exact_generic_while_check_final,
  rti_tid_mpi_select_exact_generic_while_check_final_root,
  rti_tid_mpi_select_exact_generic_while_check_post,


  rti_tid_all,
  rti_tid_sort_insert,
  rti_tid_sort_quick,
  rti_tid_sort_radix,
  rti_tid_sort_radix_iter,
  rti_tid_sort_permute_forward,
  rti_tid_sort_permute_backward,

  rti_tid_mpi_all,

  rti_tid_mpi_merge2_fe,
  rti_tid_mpi_merge2_xchg,

  rti_tid_mpi_mergek_merge2,

  rti_tid_mpi_splitk_exact,
  rti_tid_mpi_splitk_exact_init,
  rti_tid_mpi_splitk_exact_loop,
  rti_tid_mpi_splitk_exact_loop_walk,
  rti_tid_mpi_splitk_exact_loop_flow,
  rti_tid_mpi_splitk_exact_loop_flow_gather,
  rti_tid_mpi_splitk_exact_loop_flow_create,
  rti_tid_mpi_splitk_exact_loop_flow_reduce,
  rti_tid_mpi_splitk_exact_loop_flow_unbalance,
  rti_tid_mpi_splitk_exact_loop_dist,
  rti_tid_mpi_splitk_exact_loop_dist_pre,
  rti_tid_mpi_splitk_exact_loop_dist_a2av,

  rti_tid_mpi_splitk_dummy,
  rti_tid_mpi_splitk_dummy_init,
  rti_tid_mpi_splitk_dummy_loop,

  rti_tid_mpi_partition_joink,
  rti_tid_mpi_partition_joink_init,
  rti_tid_mpi_partition_joink_loop,
  rti_tid_mpi_partition_joink_loop_flow,
  rti_tid_mpi_partition_joink_loop_dist,

  rti_tid_mpi_select_radix_final,

  rti_tid_mpi_partition_radix2,
  rti_tid_mpi_partition_radix2_sync,
  rti_tid_mpi_partition_radix2_sync_init,
  rti_tid_mpi_partition_radix2_sync_exit,
/*  rti_tid_mpi_partition_radix2_while,
  rti_tid_mpi_partition_radix2_while_count,
  rti_tid_mpi_partition_radix2_while_allreduce,
  rti_tid_mpi_partition_radix2_while_round1,
  rti_tid_mpi_partition_radix2_while_round1_allgather,
  rti_tid_mpi_partition_radix2_while_exscan,
  rti_tid_mpi_partition_radix2_while_check,
  rti_tid_mpi_partition_radix2_while_check_pre,
  rti_tid_mpi_partition_radix2_while_check_classes,
  rti_tid_mpi_partition_radix2_while_check_final,
  rti_tid_mpi_partition_radix2_while_check_post,*/
  rti_tid_mpi_partition_radix2_final,

  rti_tid_mpi_sample_complete,
  rti_tid_mpi_sample_complete_gather,
  rti_tid_mpi_sample_complete_detect,
  rti_tid_mpi_sample_complete_bcast,

  rti_tid_mpi_select_qs,
  rti_tid_mpi_select_qs_pre,
  rti_tid_mpi_select_qs_loop,
  rti_tid_mpi_select_qs_part,
  rti_tid_mpi_select_qs_reduce_sizes,
  rti_tid_mpi_select_qs_area,
  rti_tid_mpi_select_qs_pivot_new,
  rti_tid_mpi_select_qs_pivot_gather,
  rti_tid_mpi_select_qs_pivot_detect,

  rti_tid_mpi_sample_select_qs,
  rti_tid_mpi_sample_select_qs_pre,
  rti_tid_mpi_sample_select_qs_select,

  rti_tid_mpi_sample_precise,
  rti_tid_mpi_sample_precise_llec,
  rti_tid_mpi_sample_precise_gather,
  rti_tid_mpi_sample_precise_detect,

  rti_tid_mpi_sample_permutation,

  rti_tid_mpi_sm_simple,
  rti_tid_mpi_sm_simple_sort,
  rti_tid_mpi_sm_simple_merge,

  rti_tids
};

#endif /* SL_RTI_TID_DECLARED */






#define pepckeys_SPEC_TLOC

typedef pepckeys_sl_int_type_c pepckeys_spec_int_t;

typedef int pepckeys_spec_proc_t;

#define pepckeys_SPEC_LOC_NONE   -1
#ifdef SL_USE_MPI
# define pepckeys_SPEC_PROC_NONE  MPI_PROC_NULL
#else
# define pepckeys_SPEC_PROC_NONE  -1
#endif

typedef void *pepckeys_spec_tloc_data_t;
typedef void *pepckeys_spec_tproc_data_t;

struct pepckeys__elements_t;

typedef struct pepckeys__elements_t *pepckeys_spec_elem_buf_t;

typedef struct pepckeys__elements_t pepckeys_spec_elem_t;

typedef pepckeys_sl_int_type_c pepckeys_spec_elem_index_t;

#define pepckeys_spec_elem_set_n(_e_, _n_)     pepckeys_elem_set_size((_e_), (_n_))
#define pepckeys_spec_elem_get_n(_e_)          pepckeys_elem_get_size((_e_))
#define pepckeys_spec_elem_set_nmax(_e_, _n_)  pepckeys_elem_set_max_size((_e_), (_n_))
#define pepckeys_spec_elem_get_nmax(_e_)       pepckeys_elem_get_max_size((_e_))

#define pepckeys_spec_elem_set_buf(_e_, _b_)   *(_e_) = *(_b_)
#define pepckeys_spec_elem_get_buf(_e_)        (_e_)

#define pepckeys_spec_elem_copy_at(_se_, _sat_, _de_, _dat_) \
  elem_copy_at((_se_), (_sat_), (_de_), (_dat_))

#define pepckeys_spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) \
  elem_xchange_at((_s0_), (_s0at_), (_s1_), (_s1at_), (_t_))






/* tproc count */

/* sp_macro pepckeys_SPEC_DECLARE_TPROC_COUNT_DB */
#define pepckeys_SPEC_DECLARE_TPROC_COUNT_DB \
  struct { pepckeys_spec_elem_index_t i; pepckeys_spec_proc_t p; } spec0cd;

/* sp_macro pepckeys_SPEC_DO_TPROC_COUNT_DB */
#define pepckeys_SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec0cd.i = 0; spec0cd.i < pepckeys_spec_elem_get_n(_b_); ++spec0cd.i) { \
    spec0cd.p = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == pepckeys_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROC_COUNT_DB */
#define pepckeys_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts) \
{ \
  pepckeys_SPEC_DECLARE_TPROC_COUNT_DB \
  pepckeys_SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TPROC_COUNT_IP */
#define pepckeys_SPEC_DECLARE_TPROC_COUNT_IP \
  struct { pepckeys_spec_elem_index_t i, t; pepckeys_spec_proc_t p; } spec0ci;

/* sp_macro pepckeys_SPEC_DO_TPROC_COUNT_IP */
#define pepckeys_SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < pepckeys_spec_elem_get_n(_b_); ++spec0ci.i) { \
    spec0ci.p = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == pepckeys_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) pepckeys_spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  pepckeys_spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROC_COUNT_IP */
#define pepckeys_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts) \
{ \
  pepckeys_SPEC_DECLARE_TPROC_COUNT_IP \
  pepckeys_SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro pepckeys_SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define pepckeys_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { pepckeys_spec_elem_index_t i; pepckeys_spec_proc_t p; } spec1cd;

/* sp_macro pepckeys_SPEC_DO_TPROC_MOD_COUNT_DB */
#define pepckeys_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec1cd.i = 0; spec1cd.i < pepckeys_spec_elem_get_n(_b_); ++spec1cd.i) { \
    spec1cd.p = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == pepckeys_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define pepckeys_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts) \
{ \
  pepckeys_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  pepckeys_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define pepckeys_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { pepckeys_spec_elem_index_t i, t; pepckeys_spec_proc_t p; } spec1ci;

/* sp_macro pepckeys_SPEC_DO_TPROC_MOD_COUNT_IP */
#define pepckeys_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < pepckeys_spec_elem_get_n(_b_); ++spec1ci.i) { \
    spec1ci.p = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == pepckeys_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) pepckeys_spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  pepckeys_spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define pepckeys_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts) \
{ \
  pepckeys_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  pepckeys_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro pepckeys_SPEC_DECLARE_TPROCS_COUNT_DB */
#define pepckeys_SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { pepckeys_spec_elem_index_t i; pepckeys_spec_int_t j, n; } spec2cd;

/* sp_macro pepckeys_SPEC_DO_TPROCS_COUNT_DB */
#define pepckeys_SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec2cd.i = 0; spec2cd.i < pepckeys_spec_elem_get_n(_b_); ++spec2cd.i) { \
    spec2cd.n = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.n; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROCS_COUNT_DB */
#define pepckeys_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts, pepckeys_spec_proc_t *procs) \
{ \
  pepckeys_SPEC_DECLARE_TPROCS_COUNT_DB \
  pepckeys_SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TPROCS_COUNT_IP */
#define pepckeys_SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { pepckeys_spec_elem_index_t i, t; pepckeys_spec_int_t j, n; } spec2ci;

/* sp_macro pepckeys_SPEC_DO_TPROCS_COUNT_IP */
#define pepckeys_SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < pepckeys_spec_elem_get_n(_b_); ++spec2ci.i) { \
    spec2ci.n = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), (_ps_)); \
    if (spec2ci.n <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.n; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) pepckeys_spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  pepckeys_spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROCS_COUNT_IP */
#define pepckeys_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts, pepckeys_spec_proc_t *procs) \
{ \
  pepckeys_SPEC_DECLARE_TPROCS_COUNT_IP \
  pepckeys_SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro pepckeys_SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define pepckeys_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { pepckeys_spec_elem_index_t i; pepckeys_spec_int_t j, n; } spec3cd;

/* sp_macro pepckeys_SPEC_DO_TPROCS_MOD_COUNT_DB */
#define pepckeys_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec3cd.i = 0; spec3cd.i < pepckeys_spec_elem_get_n(_b_); ++spec3cd.i) \
  { \
    spec3cd.n = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.n; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define pepckeys_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts, pepckeys_spec_proc_t *procs) \
{ \
  pepckeys_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  pepckeys_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define pepckeys_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { pepckeys_spec_elem_index_t i, t; pepckeys_spec_int_t j, n; } spec3ci;

/* sp_macro pepckeys_SPEC_DO_TPROCS_MOD_COUNT_IP */
#define pepckeys_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < pepckeys_spec_elem_get_n(_b_); ++spec3ci.i) { \
    spec3ci.n = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), (_ps_), NULL); \
    if (spec3ci.n <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.n; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) pepckeys_spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  pepckeys_spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define pepckeys_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts, pepckeys_spec_proc_t *procs) \
{ \
  pepckeys_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  pepckeys_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tproc rearrange */

/* sp_macro pepckeys_SPEC_DECLARE_TPROC_REARRANGE_DB */
#define pepckeys_SPEC_DECLARE_TPROC_REARRANGE_DB \
  struct { pepckeys_spec_elem_index_t i; pepckeys_spec_proc_t p; } spec0d;

/* sp_macro pepckeys_SPEC_DO_TPROC_REARRANGE_DB */
#define pepckeys_SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  for (spec0d.i = 0; spec0d.i < pepckeys_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tp_)(pepckeys_spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == pepckeys_SPEC_PROC_NONE) continue; \
    pepckeys_spec_elem_copy_at((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p]); \
    ++(_ds_)[spec0d.p]; \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROC_REARRANGE_DB */
#define pepckeys_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tproc_data_t tproc_data, int *displs) \
{ \
  pepckeys_SPEC_DECLARE_TPROC_REARRANGE_DB \
  pepckeys_SPEC_DO_TPROC_REARRANGE_DB(_tp_, tproc_data, s, d, displs); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TPROC_REARRANGE_IP */
#define pepckeys_SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { pepckeys_spec_elem_index_t e, i, j; pepckeys_spec_proc_t p, np; } spec0i;

/* sp_macro pepckeys_SPEC_DO_TPROC_REARRANGE_IP */
#define pepckeys_SPEC_DO_TPROC_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_)  do { \
  for (spec0i.e = 0, spec0i.i = 0; spec0i.i < (_n_); ++spec0i.i) { \
    spec0i.e += (_cs_)[spec0i.i]; \
    spec0i.j = (_ds_)[spec0i.i]; \
    while (spec0i.j < spec0i.e) { \
      spec0i.p = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec0i.j, _tpd_); \
      while (spec0i.p != spec0i.i) { \
        spec0i.np = (_tp_)(pepckeys_spec_elem_get_buf(_b_), (_ds_)[spec0i.p], _tpd_); \
        if (spec0i.np != spec0i.p) pepckeys_spec_elem_exchange_at((_b_), (_ds_)[spec0i.p], (_b_), spec0i.j, (_xb_)); \
        ++(_ds_)[spec0i.p]; \
        spec0i.p = spec0i.np; \
      } \
      ++spec0i.j; \
    } \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROC_REARRANGE_IP */
#define pepckeys_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_rearrange_ip(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tproc_data_t tproc_data, int *displs, int *counts, pepckeys_spec_int_t n) \
{ \
  pepckeys_SPEC_DECLARE_TPROC_REARRANGE_IP \
  pepckeys_SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro pepckeys_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define pepckeys_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  struct { pepckeys_spec_elem_index_t i; pepckeys_spec_proc_t p; } spec1d;

/* sp_macro pepckeys_SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define pepckeys_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < pepckeys_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(pepckeys_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, pepckeys_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == pepckeys_SPEC_PROC_NONE) continue; \
      pepckeys_spec_elem_copy_at((_ib_), 0, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < pepckeys_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(pepckeys_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, NULL); \
      if (spec1d.p == pepckeys_SPEC_PROC_NONE) continue; \
      pepckeys_spec_elem_copy_at((_sb_), spec1d.i, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define pepckeys_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tproc_data_t tproc_data, int *displs, pepckeys_spec_elem_t *mod) \
{ \
  pepckeys_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  pepckeys_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define pepckeys_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { pepckeys_spec_elem_index_t e, i, j; pepckeys_spec_proc_t p, np; } spec1i;

/* sp_macro pepckeys_SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define pepckeys_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
      spec1i.e += (_cs_)[spec1i.i]; \
      spec1i.j = (_ds_)[spec1i.i]; \
      while (spec1i.j < spec1i.e) { \
        spec1i.p = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec1i.j, _tpd_, pepckeys_spec_elem_get_buf(_ib_)); \
        pepckeys_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(pepckeys_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, pepckeys_spec_elem_get_buf(_ib_)); \
          if (spec1i.np != spec1i.p) { \
            pepckeys_spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
            pepckeys_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
          } else pepckeys_spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
          ++(_ds_)[spec1i.p]; \
          spec1i.p = spec1i.np; \
        } \
        ++spec1i.j; \
      } \
    } \
  } else { \
    for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
      spec1i.e += (_cs_)[spec1i.i]; \
      spec1i.j = (_ds_)[spec1i.i]; \
      while (spec1i.j < spec1i.e) { \
        spec1i.p = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec1i.j, _tpd_, NULL); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(pepckeys_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, NULL); \
          if (spec1i.np != spec1i.p) pepckeys_spec_elem_exchange_at((_b_), (_ds_)[spec1i.p], (_b_), spec1i.j, (_xb_)); \
          ++(_ds_)[spec1i.p]; \
          spec1i.p = spec1i.np; \
        } \
        ++spec1i.j; \
      } \
    } \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define pepckeys_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tproc_data_t tproc_data, int *displs, int *counts, pepckeys_spec_int_t n, pepckeys_spec_elem_t *mod) \
{ \
  pepckeys_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  pepckeys_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro pepckeys_SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define pepckeys_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  struct { pepckeys_spec_elem_index_t i; pepckeys_spec_int_t j, n; } spec2d;

/* sp_macro pepckeys_SPEC_DO_TPROCS_REARRANGE_DB */
#define pepckeys_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  for (spec2d.i = 0; spec2d.i < pepckeys_spec_elem_get_n(_sb_); ++spec2d.i) { \
    spec2d.n = (_tp_)(pepckeys_spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.n; ++spec2d.j) { \
      pepckeys_spec_elem_copy_at((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]]); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROCS_REARRANGE_DB */
#define pepckeys_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tproc_data_t tproc_data, int *displs, pepckeys_spec_proc_t *procs) \
{ \
  pepckeys_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  pepckeys_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define pepckeys_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { pepckeys_spec_elem_index_t e, j, fe, fc, le, lc; pepckeys_spec_int_t i, n, f, l, o; } spec2i;

/* sp_macro pepckeys_SPEC_DO_TPROCS_REARRANGE_IP */
#define pepckeys_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_)  do { \
  spec2i.f = 0; spec2i.fe = (_cs_)[0]; spec2i.fc = pepckeys_spec_elem_get_n(_b_); \
  while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; } \
  spec2i.l = 0; spec2i.le = (_cs_)[0]; spec2i.lc = pepckeys_spec_elem_get_n(_b_) - 1; \
  while (spec2i.lc >= spec2i.le) { ++spec2i.l; spec2i.le += (_cs_)[spec2i.l]; } \
  for (spec2i.e = 0, spec2i.i = 0; spec2i.i < (_n_); ++spec2i.i) { \
    spec2i.e += (_cs_)[spec2i.i]; \
    spec2i.j = (_ds_)[spec2i.i]; \
    while (spec2i.j < spec2i.e) { \
      spec2i.n = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec2i.j, (_tpd_), (_ps_)); \
      spec2i.o = -1; \
      while (spec2i.n > 0) { \
        --spec2i.n; \
        if ((_ps_)[spec2i.n] == spec2i.i && spec2i.o < 0) spec2i.o = spec2i.n; \
        else if ((_ds_)[(_ps_)[spec2i.n]] < spec2i.fc) { \
          spec2i.l = spec2i.f; spec2i.le = spec2i.fe; spec2i.lc = spec2i.fc; \
          if (spec2i.fc < spec2i.fe) { \
            pepckeys_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_b_), spec2i.fc); \
            ++spec2i.fc; \
          } else pepckeys_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec2i.n]] == spec2i.fc) ++spec2i.fc; \
        if (spec2i.j != (_ds_)[(_ps_)[spec2i.n]]) pepckeys_spec_elem_copy_at((_b_), spec2i.j, (_b_), (_ds_)[(_ps_)[spec2i.n]]); \
        ++(_ds_)[(_ps_)[spec2i.n]]; \
        while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; spec2i.fc = (_ds_)[spec2i.f]; } \
      } \
      if (spec2i.o < 0) { \
        if (spec2i.lc < spec2i.le) {  \
          pepckeys_spec_elem_copy_at((_b_), spec2i.lc, (_b_), spec2i.j); \
          spec2i.f = spec2i.l; spec2i.fe = spec2i.le; spec2i.fc = spec2i.lc; \
          --spec2i.lc; \
          while (spec2i.l > 0 && spec2i.lc < (_ds_)[spec2i.l]) { spec2i.le -= (_cs_)[spec2i.l]; spec2i.lc = spec2i.le - 1; --spec2i.l; } \
        } else pepckeys_spec_elem_copy_at((_xb_), 0, (_b_), spec2i.j); \
      } \
      spec2i.j = (_ds_)[spec2i.i]; \
    } \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROCS_REARRANGE_IP */
#define pepckeys_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_ip(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tproc_data_t tproc_data, int *displs, int *counts, pepckeys_spec_int_t n, pepckeys_spec_proc_t *procs) \
{ \
  pepckeys_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  pepckeys_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro pepckeys_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define pepckeys_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  struct { pepckeys_spec_elem_index_t i; pepckeys_spec_int_t j, n; } spec3d;

/* sp_macro pepckeys_SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define pepckeys_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  if (_ib_) { \
    for (spec3d.i = 0; spec3d.i < pepckeys_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(pepckeys_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), pepckeys_spec_elem_get_buf(_ib_)); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        pepckeys_spec_elem_copy_at((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } else { \
    for (spec3d.i = 0; spec3d.i < pepckeys_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(pepckeys_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), NULL); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        pepckeys_spec_elem_copy_at((_sb_), spec3d.i, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define pepckeys_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tproc_data_t tproc_data, int *displs, pepckeys_spec_proc_t *procs, pepckeys_spec_elem_t *mod) \
{ \
  pepckeys_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  pepckeys_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define pepckeys_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { pepckeys_spec_elem_index_t e, j, fe, fc, le, lc; pepckeys_spec_int_t i, n, f, l, o; } spec3i;

/* sp_macro pepckeys_SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define pepckeys_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  if (_ib_) { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = pepckeys_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = pepckeys_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), pepckeys_spec_elem_get_buf(_ib_)); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              pepckeys_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else pepckeys_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          pepckeys_spec_elem_copy_at((_ib_), spec3i.n, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            pepckeys_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else pepckeys_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } else { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = pepckeys_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = pepckeys_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(pepckeys_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), NULL); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              pepckeys_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else pepckeys_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          if (spec3i.j != (_ds_)[(_ps_)[spec3i.n]]) pepckeys_spec_elem_copy_at((_b_), spec3i.j, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            pepckeys_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else pepckeys_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define pepckeys_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tproc_data_t tproc_data, int *displs, int *counts, pepckeys_spec_int_t n, pepckeys_spec_proc_t *procs, pepckeys_spec_elem_t *mod) \
{ \
  pepckeys_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  pepckeys_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}

/* sp_macro pepckeys_SPEC_DEFINE_TPROC */
#define pepckeys_SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  pepckeys_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro pepckeys_SPEC_DEFINE_TPROC_MOD */
#define pepckeys_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  pepckeys_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro pepckeys_SPEC_DEFINE_TPROCS */
#define pepckeys_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  pepckeys_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro pepckeys_SPEC_DEFINE_TPROCS_MOD */
#define pepckeys_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  pepckeys_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  pepckeys_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro pepckeys_SPEC_EXT_PARAM_TPROC pepckeys_SPEC_EXT_PARAM_TPROC_NULL pepckeys_SPEC_EXT_PARAM_TPROC_MOD pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL pepckeys_SPEC_EXT_PARAM_TPROCS pepckeys_SPEC_EXT_PARAM_TPROCS_NULL pepckeys_SPEC_EXT_PARAM_TPROCS_MOD pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define pepckeys_SPEC_EXT_PARAM_TPROC(_name_)       _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip
#define pepckeys_SPEC_EXT_PARAM_TPROC_NULL          NULL, NULL, NULL, NULL
#define pepckeys_SPEC_EXT_PARAM_TPROC_MOD(_name_)   _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip
#define pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL      NULL, NULL, NULL, NULL
#define pepckeys_SPEC_EXT_PARAM_TPROCS(_name_)      _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip
#define pepckeys_SPEC_EXT_PARAM_TPROCS_NULL         NULL, NULL, NULL, NULL
#define pepckeys_SPEC_EXT_PARAM_TPROCS_MOD(_name_)  _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip
#define pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL     NULL, NULL, NULL, NULL


/* sp_type pepckeys_spec_tproc_f pepckeys_spec_tproc_count_f pepckeys_spec_tproc_rearrange_db_f pepckeys_spec_tproc_rearrange_ip_f */
typedef pepckeys_spec_proc_t pepckeys_spec_tproc_f(pepckeys_spec_elem_buf_t b, pepckeys_spec_elem_index_t x, pepckeys_spec_tproc_data_t tproc_data);
typedef void pepckeys_spec_tproc_count_f(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts);
typedef void pepckeys_spec_tproc_rearrange_db_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tproc_data_t tproc_data, int *displs);
typedef void pepckeys_spec_tproc_rearrange_ip_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tproc_data_t tproc_data, int *displs, int *counts, pepckeys_spec_int_t n);

/* sp_type pepckeys_spec_tproc_mod_f pepckeys_spec_tproc_mod_count_f pepckeys_spec_tproc_mod_rearrange_db_f pepckeys_spec_tproc_mod_rearrange_ip_f */
typedef pepckeys_spec_proc_t pepckeys_spec_tproc_mod_f(pepckeys_spec_elem_buf_t b, pepckeys_spec_elem_index_t x, pepckeys_spec_tproc_data_t tproc_data, pepckeys_spec_elem_buf_t mod);
typedef void pepckeys_spec_tproc_mod_count_f(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts);
typedef void pepckeys_spec_tproc_mod_rearrange_db_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tproc_data_t tproc_data, int *displs, pepckeys_spec_elem_t *mod);
typedef void pepckeys_spec_tproc_mod_rearrange_ip_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tproc_data_t tproc_data, int *displs, int *counts, pepckeys_spec_int_t n, pepckeys_spec_elem_t *mod);

/* sp_type pepckeys_spec_tprocs_f pepckeys_spec_tprocs_count_f pepckeys_spec_tprocs_rearrange_db_f pepckeys_spec_tprocs_rearrange_ip_f */
typedef pepckeys_spec_int_t pepckeys_spec_tprocs_f(pepckeys_spec_elem_buf_t b, pepckeys_spec_elem_index_t x, pepckeys_spec_tproc_data_t tproc_data, pepckeys_spec_proc_t *procs);
typedef void pepckeys_spec_tprocs_count_f(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts, pepckeys_spec_proc_t *procs);
typedef void pepckeys_spec_tprocs_rearrange_db_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tproc_data_t tproc_data, int *displs, pepckeys_spec_proc_t *procs);
typedef void pepckeys_spec_tprocs_rearrange_ip_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tproc_data_t tproc_data, int *displs, int *counts, pepckeys_spec_int_t n, pepckeys_spec_proc_t *procs);

/* sp_type pepckeys_spec_tprocs_mod_f pepckeys_spec_tprocs_mod_count_f pepckeys_spec_tprocs_mod_rearrange_db_f pepckeys_spec_tprocs_mod_rearrange_ip_f */
typedef pepckeys_spec_int_t pepckeys_spec_tprocs_mod_f(pepckeys_spec_elem_buf_t b, pepckeys_spec_elem_index_t x, pepckeys_spec_tproc_data_t tproc_data, pepckeys_spec_proc_t *procs, pepckeys_spec_elem_buf_t mod);
typedef void pepckeys_spec_tprocs_mod_count_f(pepckeys_spec_elem_t *s, pepckeys_spec_tproc_data_t tproc_data, int *counts, pepckeys_spec_proc_t *procs);
typedef void pepckeys_spec_tprocs_mod_rearrange_db_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tproc_data_t tproc_data, int *displs, pepckeys_spec_proc_t *procs, pepckeys_spec_elem_t *mod);
typedef void pepckeys_spec_tprocs_mod_rearrange_ip_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tproc_data_t tproc_data, int *displs, int *counts, pepckeys_spec_int_t n, pepckeys_spec_proc_t *procs, pepckeys_spec_elem_t *mod);

/* sp_type pepckeys_spec_tproc_reset_f */
typedef void pepckeys_spec_tproc_reset_f(pepckeys_spec_tproc_data_t tproc_data);


/* enable tloc features */
#ifdef pepckeys_SPEC_TLOC

/* sp_macro pepckeys_SPEC_TLOC pepckeys_SPEC_LOC_NONE */


/* tloc rearrange */

/* sp_macro pepckeys_SPEC_DECLARE_TLOC_REARRANGE_DB */
#define pepckeys_SPEC_DECLARE_TLOC_REARRANGE_DB \
  struct { pepckeys_spec_int_t i, p; } spec0d;

/* sp_macro pepckeys_SPEC_DO_TLOC_REARRANGE_DB */
#define pepckeys_SPEC_DO_TLOC_REARRANGE_DB(_tl_, _tld_, _sb_, _db_)  do { \
  for (spec0d.i = 0; spec0d.i < pepckeys_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tl_)(pepckeys_spec_elem_get_buf(_sb_), spec0d.i, _tld_); \
    if (spec0d.p == pepckeys_SPEC_LOC_NONE) continue; \
    pepckeys_spec_elem_copy_at((_sb_), spec0d.i, (_db_), spec0d.p); \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TLOC_REARRANGE_DB */
#define pepckeys_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_rearrange_db(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tloc_data_t tloc_data) \
{ \
  pepckeys_SPEC_DECLARE_TLOC_REARRANGE_DB \
  pepckeys_SPEC_DO_TLOC_REARRANGE_DB(_tl_, tloc_data, s, d); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TLOC_REARRANGE_IP */
#define pepckeys_SPEC_DECLARE_TLOC_REARRANGE_IP \
  struct { pepckeys_spec_int_t i, p, np; } spec0i;

/* sp_macro pepckeys_SPEC_DO_TLOC_REARRANGE_IP */
#define pepckeys_SPEC_DO_TLOC_REARRANGE_IP(_tl_, _tld_, _b_, _xb_)  do { \
  for (spec0i.i = 0; spec0i.i < pepckeys_spec_elem_get_n(_b_); ++spec0i.i) { \
    spec0i.p = (_tl_)(pepckeys_spec_elem_get_buf(_b_), spec0i.i, _tld_); \
    if (spec0i.p == pepckeys_SPEC_LOC_NONE) continue; \
    while (spec0i.i != spec0i.p) { \
      spec0i.np = (_tl_)(pepckeys_spec_elem_get_buf(_b_), spec0i.p, _tld_); \
      if (spec0i.np == pepckeys_SPEC_LOC_NONE) { pepckeys_spec_elem_copy_at((_b_), spec0i.i, (_b_), spec0i.p); break; } \
      pepckeys_spec_elem_exchange_at((_b_), spec0i.i, (_b_), spec0i.p, (_xb_)); \
      spec0i.p = spec0i.np; \
    } \
  } } while (0)

/* sp_macro pepckeys_SPEC_FUNC_TLOC_REARRANGE_IP */
#define pepckeys_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_rearrange_ip(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tloc_data_t tloc_data) \
{ \
  pepckeys_SPEC_DECLARE_TLOC_REARRANGE_IP \
  pepckeys_SPEC_DO_TLOC_REARRANGE_IP(_tl_, tloc_data, s, x); \
}


/* tloc_mod_mod rearrange */

/* sp_macro pepckeys_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB */
#define pepckeys_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  struct { pepckeys_spec_int_t i, p; } spec1d;

/* sp_macro pepckeys_SPEC_DO_TLOC_MOD_REARRANGE_DB */
#define pepckeys_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, _tld_, _sb_, _db_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < pepckeys_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(pepckeys_spec_elem_get_buf(_sb_), spec1d.i, _tld_, pepckeys_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == pepckeys_SPEC_LOC_NONE) continue; \
      pepckeys_spec_elem_copy_at((_ib_), 0, (_db_), spec1d.p); \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < pepckeys_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(pepckeys_spec_elem_get_buf(_sb_), spec1d.i, _tld_, NULL); \
      if (spec1d.p == pepckeys_SPEC_LOC_NONE) continue; \
      pepckeys_spec_elem_copy_at((_sb_), spec1d.i, (_db_), spec1d.p); \
    } \
  } } while (0) 

/* sp_macro pepckeys_SPEC_FUNC_TLOC_MOD_REARRANGE_DB */
#define pepckeys_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_mod_rearrange_db(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tloc_data_t tloc_data, pepckeys_spec_elem_t *mod) \
{ \
  pepckeys_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  pepckeys_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, tloc_data, s, d, mod); \
}

/* sp_macro pepckeys_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP */
#define pepckeys_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  struct { pepckeys_spec_int_t i, p, np; } spec1i;

/* sp_macro pepckeys_SPEC_DO_TLOC_MOD_REARRANGE_IP */
#define pepckeys_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, _tld_, _b_, _xb_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.i = 0; spec1i.i < pepckeys_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(pepckeys_spec_elem_get_buf(_b_), spec1i.i, _tld_, pepckeys_spec_elem_get_buf(_ib_)); \
      if (spec1i.p == pepckeys_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(pepckeys_spec_elem_get_buf(_b_), spec1i.p, _tld_, pepckeys_spec_elem_get_buf(_xb_)); \
        if (spec1i.np == pepckeys_SPEC_LOC_NONE) break; \
        pepckeys_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.p); \
        pepckeys_spec_elem_copy_at((_xb_), 0, (_ib_), 0); \
        spec1i.p = spec1i.np; \
      } \
      pepckeys_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.i); \
    } \
  } else { \
    for (spec1i.i = 0; spec1i.i < pepckeys_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(pepckeys_spec_elem_get_buf(_b_), spec1i.i, _tld_, NULL); \
      if (spec1i.p == pepckeys_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(pepckeys_spec_elem_get_buf(_b_), spec1i.p, _tld_, NULL); \
        if (spec1i.np == pepckeys_SPEC_LOC_NONE) { pepckeys_spec_elem_copy_at((_b_), spec1i.i, (_b_), spec1i.p); break; } \
        pepckeys_spec_elem_exchange_at((_b_), spec1i.i, (_b_), spec1i.p, (_xb_)); \
        spec1i.p = spec1i.np; \
      } \
    } \
 } } while (0) 

/* sp_macro pepckeys_SPEC_FUNC_TLOC_MOD_REARRANGE_IP */
#define pepckeys_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_mod_rearrange_ip(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tloc_data_t tloc_data, pepckeys_spec_elem_t *mod) \
{ \
  pepckeys_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  pepckeys_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, tloc_data, s, x, mod); \
}

/* sp_macro pepckeys_SPEC_DEFINE_TLOC */
#define pepckeys_SPEC_DEFINE_TLOC(_name_, _tl_, _s_...) \
  pepckeys_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_) \
  pepckeys_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro pepckeys_SPEC_DEFINE_TLOC_MOD */
#define pepckeys_SPEC_DEFINE_TLOC_MOD(_name_, _tl_, _s_...) \
  pepckeys_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_) \
  pepckeys_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro pepckeys_SPEC_EXT_PARAM_TLOC pepckeys_SPEC_EXT_PARAM_TLOC_NULL pepckeys_SPEC_EXT_PARAM_TLOC_MOD pepckeys_SPEC_EXT_PARAM_TLOC_MOD_NULL */
#define pepckeys_SPEC_EXT_PARAM_TLOC(_name_)      _name_##_tloc_rearrange_db, _name_##_tloc_rearrange_ip
#define pepckeys_SPEC_EXT_PARAM_TLOC_NULL         NULL, NULL
#define pepckeys_SPEC_EXT_PARAM_TLOC_MOD(_name_)  _name_##_tloc_mod_rearrange_db, _name_##_tloc_mod_rearrange_ip
#define pepckeys_SPEC_EXT_PARAM_TLOC_MOD_NULL     NULL, NULL


/* sp_type pepckeys_spec_tloc_f pepckeys_spec_tloc_rearrange_db_f pepckeys_spec_tloc_rearrange_ip_f */
typedef pepckeys_spec_elem_index_t pepckeys_spec_tloc_f(pepckeys_spec_elem_buf_t b, pepckeys_spec_elem_index_t x, pepckeys_spec_tloc_data_t tloc_data);
typedef void pepckeys_spec_tloc_rearrange_db_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tloc_data_t tloc_data);
typedef void pepckeys_spec_tloc_rearrange_ip_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tloc_data_t tloc_data);

/* sp_type pepckeys_spec_tloc_mod_f pepckeys_spec_tloc_mod_rearrange_db_f pepckeys_spec_tloc_mod_rearrange_ip_f */
typedef pepckeys_spec_elem_index_t pepckeys_spec_tloc_mod_f(pepckeys_spec_elem_buf_t b, pepckeys_spec_elem_index_t x, pepckeys_spec_tloc_data_t tloc_data, pepckeys_spec_elem_buf_t mod);
typedef void pepckeys_spec_tloc_mod_rearrange_db_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *d, pepckeys_spec_tloc_data_t tloc_data, pepckeys_spec_elem_t *mod);
typedef void pepckeys_spec_tloc_mod_rearrange_ip_f(pepckeys_spec_elem_t *s, pepckeys_spec_elem_t *x, pepckeys_spec_tloc_data_t tloc_data, pepckeys_spec_elem_t *mod);


#endif /* pepckeys_SPEC_TLOC */






#ifdef SL_USE_MPI
# include <mpi.h>
#endif


/* sl_type pepckeys_slint_t pepckeys_slint */
typedef pepckeys_sl_int_type_c pepckeys_slint_t, pepckeys_slint;  /* deprecated 'pepckeys_slint' */

#define pepckeys_slint_fmt   pepckeys_sl_int_type_fmt    /* sl_macro */

/* sl_type pepckeys_slindex_t */
typedef pepckeys_sl_index_type_c pepckeys_slindex_t;

#define pepckeys_sindex_fmt  pepckeys_sl_index_type_fmt  /* sl_macro */

/* sl_type pepckeys_slkey_t */
typedef pepckeys_sl_key_type_c pepckeys_slkey_t;

/* sl_type pepckeys_slkey_pure_t pepckeys_slpkey_t */
typedef pepckeys_sl_key_pure_type_c pepckeys_slkey_pure_t, pepckeys_slpkey_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type pepckeys_sldata0_t */
#ifdef pepckeys_sl_data0_type_c
typedef pepckeys_sl_data0_type_c pepckeys_sldata0_t;
#endif
/* sl_type pepckeys_sldata1_t */
#ifdef pepckeys_sl_data1_type_c
typedef pepckeys_sl_data1_type_c pepckeys_sldata1_t;
#endif
/* sl_type pepckeys_sldata2_t */
#ifdef pepckeys_sl_data2_type_c
typedef pepckeys_sl_data2_type_c pepckeys_sldata2_t;
#endif
/* sl_type pepckeys_sldata3_t */
#ifdef pepckeys_sl_data3_type_c
typedef pepckeys_sl_data3_type_c pepckeys_sldata3_t;
#endif
/* sl_type pepckeys_sldata4_t */
#ifdef pepckeys_sl_data4_type_c
typedef pepckeys_sl_data4_type_c pepckeys_sldata4_t;
#endif
/* sl_type pepckeys_sldata5_t */
#ifdef pepckeys_sl_data5_type_c
typedef pepckeys_sl_data5_type_c pepckeys_sldata5_t;
#endif
/* sl_type pepckeys_sldata6_t */
#ifdef pepckeys_sl_data6_type_c
typedef pepckeys_sl_data6_type_c pepckeys_sldata6_t;
#endif
/* sl_type pepckeys_sldata7_t */
#ifdef pepckeys_sl_data7_type_c
typedef pepckeys_sl_data7_type_c pepckeys_sldata7_t;
#endif
/* sl_type pepckeys_sldata8_t */
#ifdef pepckeys_sl_data8_type_c
typedef pepckeys_sl_data8_type_c pepckeys_sldata8_t;
#endif
/* sl_type pepckeys_sldata9_t */
#ifdef pepckeys_sl_data9_type_c
typedef pepckeys_sl_data9_type_c pepckeys_sldata9_t;
#endif
/* sl_type pepckeys_sldata10_t */
#ifdef pepckeys_sl_data10_type_c
typedef pepckeys_sl_data10_type_c pepckeys_sldata10_t;
#endif
/* sl_type pepckeys_sldata11_t */
#ifdef pepckeys_sl_data11_type_c
typedef pepckeys_sl_data11_type_c pepckeys_sldata11_t;
#endif
/* sl_type pepckeys_sldata12_t */
#ifdef pepckeys_sl_data12_type_c
typedef pepckeys_sl_data12_type_c pepckeys_sldata12_t;
#endif
/* sl_type pepckeys_sldata13_t */
#ifdef pepckeys_sl_data13_type_c
typedef pepckeys_sl_data13_type_c pepckeys_sldata13_t;
#endif
/* sl_type pepckeys_sldata14_t */
#ifdef pepckeys_sl_data14_type_c
typedef pepckeys_sl_data14_type_c pepckeys_sldata14_t;
#endif
/* sl_type pepckeys_sldata15_t */
#ifdef pepckeys_sl_data15_type_c
typedef pepckeys_sl_data15_type_c pepckeys_sldata15_t;
#endif
/* sl_type pepckeys_sldata16_t */
#ifdef pepckeys_sl_data16_type_c
typedef pepckeys_sl_data16_type_c pepckeys_sldata16_t;
#endif
/* sl_type pepckeys_sldata17_t */
#ifdef pepckeys_sl_data17_type_c
typedef pepckeys_sl_data17_type_c pepckeys_sldata17_t;
#endif
/* sl_type pepckeys_sldata18_t */
#ifdef pepckeys_sl_data18_type_c
typedef pepckeys_sl_data18_type_c pepckeys_sldata18_t;
#endif
/* sl_type pepckeys_sldata19_t */
#ifdef pepckeys_sl_data19_type_c
typedef pepckeys_sl_data19_type_c pepckeys_sldata19_t;
#endif
/* DATAX_TEMPLATE_END */

#define SL_DATA_NMAX (0 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
 + 1 \
)

/* sl_type pepckeys_slweight_t */
typedef pepckeys_sl_weight_type_c pepckeys_slweight_t;

#define pepckeys_slweight_fmt  pepckeys_sl_weight_type_fmt  /* sl_macro */

#if defined(pepckeys_sl_elem_weight) && defined(pepckeys_sl_weight_intequiv)
typedef pepckeys_sl_weight_type_c pepckeys_slcount_t;       /* sl_type pepckeys_slcount_t */
# define pepckeys_slcount_fmt  pepckeys_sl_weight_type_fmt  /* sl_macro */
#else
typedef pepckeys_sl_int_type_c pepckeys_slcount_t;
# define pepckeys_slcount_fmt  pepckeys_sl_int_type_fmt
#endif


/* sl_type pepckeys__slpwkey_t pepckeys_slpwkey_t */
typedef struct pepckeys__slpwkey_t
{
  pepckeys_slpkey_t pkey;
  pepckeys_slweight_t weight;

} pepckeys_slpwkey_t;


/* sl_type pepckeys__elements_t pepckeys_elements_t */
typedef struct pepckeys__elements_t
{
  pepckeys_slint_t size, max_size;
  pepckeys_slkey_t *keys;

#ifdef pepckeys_SL_INDEX
  pepckeys_slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepckeys_SL_DATA0
  pepckeys_sldata0_t *data0;
#endif
#ifdef pepckeys_SL_DATA1
  pepckeys_sldata1_t *data1;
#endif
#ifdef pepckeys_SL_DATA2
  pepckeys_sldata2_t *data2;
#endif
#ifdef pepckeys_SL_DATA3
  pepckeys_sldata3_t *data3;
#endif
#ifdef pepckeys_SL_DATA4
  pepckeys_sldata4_t *data4;
#endif
#ifdef pepckeys_SL_DATA5
  pepckeys_sldata5_t *data5;
#endif
#ifdef pepckeys_SL_DATA6
  pepckeys_sldata6_t *data6;
#endif
#ifdef pepckeys_SL_DATA7
  pepckeys_sldata7_t *data7;
#endif
#ifdef pepckeys_SL_DATA8
  pepckeys_sldata8_t *data8;
#endif
#ifdef pepckeys_SL_DATA9
  pepckeys_sldata9_t *data9;
#endif
#ifdef pepckeys_SL_DATA10
  pepckeys_sldata10_t *data10;
#endif
#ifdef pepckeys_SL_DATA11
  pepckeys_sldata11_t *data11;
#endif
#ifdef pepckeys_SL_DATA12
  pepckeys_sldata12_t *data12;
#endif
#ifdef pepckeys_SL_DATA13
  pepckeys_sldata13_t *data13;
#endif
#ifdef pepckeys_SL_DATA14
  pepckeys_sldata14_t *data14;
#endif
#ifdef pepckeys_SL_DATA15
  pepckeys_sldata15_t *data15;
#endif
#ifdef pepckeys_SL_DATA16
  pepckeys_sldata16_t *data16;
#endif
#ifdef pepckeys_SL_DATA17
  pepckeys_sldata17_t *data17;
#endif
#ifdef pepckeys_SL_DATA18
  pepckeys_sldata18_t *data18;
#endif
#ifdef pepckeys_SL_DATA19
  pepckeys_sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} pepckeys_elements_t;


/* sl_type pepckeys__packed_element_t pepckeys_packed_element_t */
typedef struct pepckeys__packed_element_t
{
  pepckeys_slkey_t key;

#ifdef pepckeys_SL_PACKED_INDEX
  pepckeys_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepckeys_SL_DATA0
# ifdef pepckeys_sl_data0_flex
  pepckeys_sldata0_t data0[];
# else
  pepckeys_sldata0_t data0[pepckeys_sl_data0_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA1
# ifdef pepckeys_sl_data1_flex
  pepckeys_sldata1_t data1[];
# else
  pepckeys_sldata1_t data1[pepckeys_sl_data1_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA2
# ifdef pepckeys_sl_data2_flex
  pepckeys_sldata2_t data2[];
# else
  pepckeys_sldata2_t data2[pepckeys_sl_data2_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA3
# ifdef pepckeys_sl_data3_flex
  pepckeys_sldata3_t data3[];
# else
  pepckeys_sldata3_t data3[pepckeys_sl_data3_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA4
# ifdef pepckeys_sl_data4_flex
  pepckeys_sldata4_t data4[];
# else
  pepckeys_sldata4_t data4[pepckeys_sl_data4_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA5
# ifdef pepckeys_sl_data5_flex
  pepckeys_sldata5_t data5[];
# else
  pepckeys_sldata5_t data5[pepckeys_sl_data5_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA6
# ifdef pepckeys_sl_data6_flex
  pepckeys_sldata6_t data6[];
# else
  pepckeys_sldata6_t data6[pepckeys_sl_data6_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA7
# ifdef pepckeys_sl_data7_flex
  pepckeys_sldata7_t data7[];
# else
  pepckeys_sldata7_t data7[pepckeys_sl_data7_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA8
# ifdef pepckeys_sl_data8_flex
  pepckeys_sldata8_t data8[];
# else
  pepckeys_sldata8_t data8[pepckeys_sl_data8_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA9
# ifdef pepckeys_sl_data9_flex
  pepckeys_sldata9_t data9[];
# else
  pepckeys_sldata9_t data9[pepckeys_sl_data9_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA10
# ifdef pepckeys_sl_data10_flex
  pepckeys_sldata10_t data10[];
# else
  pepckeys_sldata10_t data10[pepckeys_sl_data10_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA11
# ifdef pepckeys_sl_data11_flex
  pepckeys_sldata11_t data11[];
# else
  pepckeys_sldata11_t data11[pepckeys_sl_data11_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA12
# ifdef pepckeys_sl_data12_flex
  pepckeys_sldata12_t data12[];
# else
  pepckeys_sldata12_t data12[pepckeys_sl_data12_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA13
# ifdef pepckeys_sl_data13_flex
  pepckeys_sldata13_t data13[];
# else
  pepckeys_sldata13_t data13[pepckeys_sl_data13_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA14
# ifdef pepckeys_sl_data14_flex
  pepckeys_sldata14_t data14[];
# else
  pepckeys_sldata14_t data14[pepckeys_sl_data14_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA15
# ifdef pepckeys_sl_data15_flex
  pepckeys_sldata15_t data15[];
# else
  pepckeys_sldata15_t data15[pepckeys_sl_data15_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA16
# ifdef pepckeys_sl_data16_flex
  pepckeys_sldata16_t data16[];
# else
  pepckeys_sldata16_t data16[pepckeys_sl_data16_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA17
# ifdef pepckeys_sl_data17_flex
  pepckeys_sldata17_t data17[];
# else
  pepckeys_sldata17_t data17[pepckeys_sl_data17_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA18
# ifdef pepckeys_sl_data18_flex
  pepckeys_sldata18_t data18[];
# else
  pepckeys_sldata18_t data18[pepckeys_sl_data18_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA19
# ifdef pepckeys_sl_data19_flex
  pepckeys_sldata19_t data19[];
# else
  pepckeys_sldata19_t data19[pepckeys_sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} pepckeys_packed_element_t;


/* sl_type pepckeys__packed_elements_t pepckeys_packed_elements_t */
typedef struct pepckeys__packed_elements_t
{
  pepckeys_slint_t size, max_size;
  
  pepckeys_packed_element_t *elements;
  
} pepckeys_packed_elements_t;


#ifndef SLCINT_T
#define SLCINT_T
typedef long long int slcint_t;
#define slcint_fmt  "ll"
/*#define slcint_sfx  LL*/
#endif


#define SLCM_KEYS     (((slcint_t) 1) << 0)
#define SLCM_INDICES  (((slcint_t) 1) << 1)
#define SLCM_WEIGHTS  (((slcint_t) 1) << 2)

/* DATAX_TEMPLATE_BEGIN */
#define SLCM_DATA0    (((slcint_t) 1) << (3+0))
#define SLCM_DATA1    (((slcint_t) 1) << (3+1))
#define SLCM_DATA2    (((slcint_t) 1) << (3+2))
#define SLCM_DATA3    (((slcint_t) 1) << (3+3))
#define SLCM_DATA4    (((slcint_t) 1) << (3+4))
#define SLCM_DATA5    (((slcint_t) 1) << (3+5))
#define SLCM_DATA6    (((slcint_t) 1) << (3+6))
#define SLCM_DATA7    (((slcint_t) 1) << (3+7))
#define SLCM_DATA8    (((slcint_t) 1) << (3+8))
#define SLCM_DATA9    (((slcint_t) 1) << (3+9))
#define SLCM_DATA10    (((slcint_t) 1) << (3+10))
#define SLCM_DATA11    (((slcint_t) 1) << (3+11))
#define SLCM_DATA12    (((slcint_t) 1) << (3+12))
#define SLCM_DATA13    (((slcint_t) 1) << (3+13))
#define SLCM_DATA14    (((slcint_t) 1) << (3+14))
#define SLCM_DATA15    (((slcint_t) 1) << (3+15))
#define SLCM_DATA16    (((slcint_t) 1) << (3+16))
#define SLCM_DATA17    (((slcint_t) 1) << (3+17))
#define SLCM_DATA18    (((slcint_t) 1) << (3+18))
#define SLCM_DATA19    (((slcint_t) 1) << (3+19))
/* DATAX_TEMPLATE_END */

#define SLCM_DATA     (((slcint_t) 0) \
  |SLCM_DATA0 \
  |SLCM_DATA1 \
  |SLCM_DATA2 \
  |SLCM_DATA3 \
  |SLCM_DATA4 \
  |SLCM_DATA5 \
  |SLCM_DATA6 \
  |SLCM_DATA7 \
  |SLCM_DATA8 \
  |SLCM_DATA9 \
  |SLCM_DATA10 \
  |SLCM_DATA11 \
  |SLCM_DATA12 \
  |SLCM_DATA13 \
  |SLCM_DATA14 \
  |SLCM_DATA15 \
  |SLCM_DATA16 \
  |SLCM_DATA17 \
  |SLCM_DATA18 \
  |SLCM_DATA19 \
  )

#define SLCM_ALL      (~((slcint_t) 0))


/* sl_type pepckeys__classification_info_t pepckeys_classification_info_t pepckeys_classification_info */
typedef struct pepckeys__classification_info_t
{
  pepckeys_slint_t nclasses;
  pepckeys_slkey_pure_t *keys;
  pepckeys_slint_t *counts;
  pepckeys_slint_t *masks;

  /* */
  pepckeys_slint_t *all_local_sizes;
  pepckeys_slint_t *local_lt_eq_counts;
  pepckeys_slint_t *all_local_lt_eq_counts;

} pepckeys_classification_info_t, pepckeys_classification_info;  /* deprecated 'pepckeys_classification_info' */


/* key2class, sl_type pepckeys_key2class_f */
typedef pepckeys_slint_t (*pepckeys_key2class_f)(pepckeys_slkey_t *, pepckeys_slint, void *);

/* pivot-element, sl_type pepckeys_pivot_f */
typedef pepckeys_slint_t (*pepckeys_pivot_f)(pepckeys_elements_t *);

/* sorting-network, sl_type pepckeys_sortnet_f pepckeys_sortnet_data_t */
typedef void *pepckeys_sortnet_data_t;
typedef pepckeys_slint_t (*pepckeys_sortnet_f)(pepckeys_slint_t size, pepckeys_slint_t rank, pepckeys_slint_t stage, pepckeys_sortnet_data_t snd, pepckeys_slint_t *up);

/* merge2, sl_type pepckeys_merge2x_f pepckeys_merge2X_f */
typedef pepckeys_slint_t (*pepckeys_merge2x_f)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
typedef pepckeys_slint_t (*pepckeys_merge2X_f)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);

/* sl_type pepckeys__permute_generic_t pepckeys_permute_generic_t */
typedef struct pepckeys__permute_generic_t
{
  int type;

  pepckeys_spec_tloc_f *tloc;
  pepckeys_spec_tloc_rearrange_db_f *tloc_rearrange_db;
  pepckeys_spec_tloc_rearrange_ip_f *tloc_rearrange_ip;

  pepckeys_spec_tloc_mod_f *tloc_mod;
  pepckeys_spec_tloc_mod_rearrange_db_f *tloc_mod_rearrange_db;
  pepckeys_spec_tloc_mod_rearrange_ip_f *tloc_mod_rearrange_ip;

} pepckeys_permute_generic_t;

/* sl_macro pepckeys_PERMUTE_GENERIC_DEFINE_TLOC pepckeys_PERMUTE_GENERIC_INIT_TLOC pepckeys_PERMUTE_GENERIC_INIT_EXT_TLOC */
#define pepckeys_PERMUTE_GENERIC_DEFINE_TLOC(_tl_, _s_...)      pepckeys_SPEC_DEFINE_TLOC(_tl_, _tl_, _s_)
#define pepckeys_PERMUTE_GENERIC_INIT_TLOC(_tl_)                { 1, _tl_, pepckeys_SPEC_EXT_PARAM_TLOC_NULL,  NULL, pepckeys_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define pepckeys_PERMUTE_GENERIC_INIT_EXT_TLOC(_tl_)            { 1, _tl_, pepckeys_SPEC_EXT_PARAM_TLOC(_tl_), NULL, pepckeys_SPEC_EXT_PARAM_TLOC_MOD_NULL }

/* sl_macro pepckeys_PERMUTE_GENERIC_DEFINE_TLOC_MOD pepckeys_PERMUTE_GENERIC_INIT_TLOC_MOD pepckeys_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD */
#define pepckeys_PERMUTE_GENERIC_DEFINE_TLOC_MOD(_tl_, _s_...)  pepckeys_SPEC_DEFINE_TLOC_MOD(_tl_, _tl_, _s_)
#define pepckeys_PERMUTE_GENERIC_INIT_TLOC_MOD(_tl_)            { 2, NULL, pepckeys_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, pepckeys_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define pepckeys_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD(_tl_)        { 2, NULL, pepckeys_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, pepckeys_SPEC_EXT_PARAM_TLOC_MOD(_tl_) }

/* sl_type pepckeys__split_generic_t pepckeys_split_generic_t */
typedef struct pepckeys__split_generic_t
{
  int type;

  pepckeys_spec_tproc_f *tproc;
  pepckeys_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  pepckeys_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  pepckeys_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  pepckeys_spec_tproc_mod_f *tproc_mod;
  pepckeys_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  pepckeys_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  pepckeys_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  pepckeys_spec_tprocs_f *tprocs;
  pepckeys_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  pepckeys_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  pepckeys_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  pepckeys_spec_tprocs_mod_f *tprocs_mod;
  pepckeys_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  pepckeys_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  pepckeys_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

  pepckeys_spec_tproc_reset_f *reset;

} pepckeys_split_generic_t;

/* sl_macro pepckeys_SPLIT_GENERIC_DEFINE_TPROC pepckeys_SPLIT_GENERIC_INIT_TPROC pepckeys_SPLIT_GENERIC_INIT_EXT_TPROC */
#define pepckeys_SPLIT_GENERIC_DEFINE_TPROC(_tp_, _s_...)         pepckeys_SPEC_DEFINE_TPROC(_tp_, _tp_, _s_)
#define pepckeys_SPLIT_GENERIC_INIT_TPROC(_tp_, _r_...)           { 1, _tp_, pepckeys_SPEC_EXT_PARAM_TPROC_NULL,  NULL, pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define pepckeys_SPLIT_GENERIC_INIT_EXT_TPROC(_tp_, _r_...)       { 1, _tp_, pepckeys_SPEC_EXT_PARAM_TPROC(_tp_), NULL, pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro pepckeys_SPLIT_GENERIC_DEFINE_TPROC_MOD pepckeys_SPLIT_GENERIC_INIT_TPROC_MOD pepckeys_SPLIT_GENERIC_INIT_EXT_TPROC_MOD */
#define pepckeys_SPLIT_GENERIC_DEFINE_TPROC_MOD(_tp_, _s_...)     pepckeys_SPEC_DEFINE_TPROC_MOD(_tp_, _tp_, _s_)
#define pepckeys_SPLIT_GENERIC_INIT_TPROC_MOD(_tp_, _r_...)       { 2, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_NULL, _tp_, pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL,  NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define pepckeys_SPLIT_GENERIC_INIT_EXT_TPROC_MOD(_tp_, _r_...)   { 2, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_NULL, _tp_, pepckeys_SPEC_EXT_PARAM_TPROC_MOD(_tp_), NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro pepckeys_SPLIT_GENERIC_DEFINE_TPROCS pepckeys_SPLIT_GENERIC_INIT_TPROCS pepckeys_SPLIT_GENERIC_INIT_EXT_TPROCS */
#define pepckeys_SPLIT_GENERIC_DEFINE_TPROCS(_tp_, _s_...)        pepckeys_SPEC_DEFINE_TPROCS(_tp_, _tp_, _s_)
#define pepckeys_SPLIT_GENERIC_INIT_TPROCS(_tp_, _r_...)          { 3, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, pepckeys_SPEC_EXT_PARAM_TPROCS_NULL,  NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define pepckeys_SPLIT_GENERIC_INIT_EXT_TPROCS(_tp_, _r_...)      { 3, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, pepckeys_SPEC_EXT_PARAM_TPROCS(_tp_), NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro pepckeys_SPLIT_GENERIC_DEFINE_TPROCS_MOD pepckeys_SPLIT_GENERIC_INIT_TPROCS_MOD pepckeys_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD */
#define pepckeys_SPLIT_GENERIC_DEFINE_TPROCS_MOD(_tp_, _s_...)    pepckeys_SPEC_DEFINE_TPROCS_MOD(_tp_, _tp_, _s_)
#define pepckeys_SPLIT_GENERIC_INIT_TPROCS_MOD(_tp_, _r_...)      { 4, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define pepckeys_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD(_tp_, _r_...)  { 4, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD(_tp_), _r_ }

/* sl_type pepckeys_tloc_f pepckeys_tloc_mod_f */
typedef pepckeys_slint_t pepckeys_tloc_f(pepckeys_elements_t *b, pepckeys_slint_t x, void *tloc_data);
typedef pepckeys_slint_t pepckeys_tloc_mod_f(pepckeys_elements_t *b, pepckeys_slint_t x, void *tloc_data, pepckeys_elements_t *mod);

/* sl_type pepckeys_tproc_f pepckeys_tproc_mod_f pepckeys_tprocs_f pepckeys_tprocs_mod_f */
typedef int pepckeys_tproc_f(pepckeys_elements_t *b, pepckeys_slint_t x, void *tproc_data);
typedef int pepckeys_tproc_mod_f(pepckeys_elements_t *b, pepckeys_slint_t x, void *tproc_data, pepckeys_elements_t *mod);
typedef pepckeys_slint_t pepckeys_tprocs_f(pepckeys_elements_t *b, pepckeys_slint_t x, void *tproc_data, int *procs);
typedef pepckeys_slint_t pepckeys_tprocs_mod_f(pepckeys_elements_t *b, pepckeys_slint_t x, void *tproc_data, int *procs, pepckeys_elements_t *mod);

/* sl_type pepckeys_tproc_reset_f */
typedef void pepckeys_tproc_reset_f(void *tproc_data);

/* sl_macro pepckeys_TPROC_RESET_NULL */
#define pepckeys_TPROC_RESET_NULL  NULL

/* sl_type pepckeys__tproc_t pepckeys_tproc_t */
typedef struct pepckeys__tproc_t *pepckeys_tproc_t;

/* sl_type pepckeys__tproc_exdef pepckeys_tproc_exdef */
typedef struct pepckeys__tproc_exdef {
  int type;

  pepckeys_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  pepckeys_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  pepckeys_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  pepckeys_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  pepckeys_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  pepckeys_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  pepckeys_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  pepckeys_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  pepckeys_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  pepckeys_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  pepckeys_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  pepckeys_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

} const *pepckeys_tproc_exdef;

/* sl_macro pepckeys_TPROC_EXDEF_NULL */
#define pepckeys_TPROC_EXDEF_NULL  NULL

/* sl_macro pepckeys_TPROC_EXDEF_DEFINE_TPROC pepckeys_TPROC_EXDEF_DEFINE_TPROC_MOD pepckeys_TPROC_EXDEF_DEFINE_TPROCS pepckeys_TPROC_EXDEF_DEFINE_TPROCS_MOD */
#define pepckeys_TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  pepckeys_SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct pepckeys__tproc_exdef _##_name_ = { 1, pepckeys_SPEC_EXT_PARAM_TPROC(_name_), pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define pepckeys_TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  pepckeys_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct pepckeys__tproc_exdef _##_name_ = { 2, pepckeys_SPEC_EXT_PARAM_TPROC_NULL, pepckeys_SPEC_EXT_PARAM_TPROC_MOD(_name_), pepckeys_SPEC_EXT_PARAM_TPROCS_NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define pepckeys_TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  pepckeys_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct pepckeys__tproc_exdef _##_name_ = { 3, pepckeys_SPEC_EXT_PARAM_TPROC_NULL, pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL, pepckeys_SPEC_EXT_PARAM_TPROCS(_name_), pepckeys_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define pepckeys_TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  pepckeys_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct pepckeys__tproc_exdef _##_name_ = { 4, pepckeys_SPEC_EXT_PARAM_TPROC_NULL, pepckeys_SPEC_EXT_PARAM_TPROC_MOD_NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_NULL, pepckeys_SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;


/* deprecated, sl_type pepckeys_k2c_func pepckeys_pivot_func pepckeys_sn_func pepckeys_m2x_func pepckeys_m2X_func */
typedef pepckeys_key2class_f pepckeys_k2c_func;
typedef pepckeys_pivot_f pepckeys_pivot_func;
typedef pepckeys_sortnet_f pepckeys_sn_func;
typedef pepckeys_merge2x_f pepckeys_m2x_func;
typedef pepckeys_merge2X_f pepckeys_m2X_func;


/* sl_type pepckeys__mergek_t pepckeys_mergek_t */
typedef struct pepckeys__mergek_t
{
  pepckeys_sortnet_f sn;
  pepckeys_sortnet_data_t snd;

  pepckeys_merge2x_f m2x;
  pepckeys_elements_t *sx;

} pepckeys_mergek_t;


/* sl_type pepckeys_keys_init_type_t pepckeys_keys_init_data_t */
typedef pepckeys_slint_t pepckeys_keys_init_type_t;
typedef void *pepckeys_keys_init_data_t;

/* sl_type pepckeys_key_set_data_t pepckeys_key_set_f */
typedef void *pepckeys_key_set_data_t;
typedef void (*pepckeys_key_set_f)(pepckeys_slkey_pure_t *k, pepckeys_key_set_data_t d);


#undef SL_EKIT_SET
#define SL_EKIT_SET         1
#undef SL_EKIT_SET_FUNC
#define SL_EKIT_SET_FUNC    2
#undef SL_EKIT_RAND
#define SL_EKIT_RAND        3
#undef SL_EKIT_RAND_QUAD
#define SL_EKIT_RAND_QUAD   4
#undef SL_EKIT_RAND_AND
#define SL_EKIT_RAND_AND    5
#undef SL_EKIT_URAND
#define SL_EKIT_URAND       6
#undef SL_EKIT_NRAND
#define SL_EKIT_NRAND       7


#ifndef SL_EIK_OFFSET
# define SL_EIK_OFFSET     65536LL
#endif

#ifndef SL_EIK_SET
# define SL_EIK_SET        SL_EIK_OFFSET*1
#endif

#ifndef SL_EIK_RAND
# define SL_EIK_RAND       SL_EIK_OFFSET*2
#endif

#ifndef SL_EIK_RAND_QUAD
# define SL_EIK_RAND_QUAD  SL_EIK_OFFSET*3
#endif

#ifndef SL_EIK_RAND_AND
# define SL_EIK_RAND_AND   SL_EIK_OFFSET*4
#endif

#ifndef SL_EIK_RAND_NORM
# define SL_EIK_RAND_NORM  SL_EIK_OFFSET*5
#endif


/* pepckeys_elements_keys_stats */
#ifndef SL_EKS_MIN
# define SL_EKS_MIN   0
#endif

#ifndef SL_EKS_MAX
# define SL_EKS_MAX   1
#endif

#ifndef SL_EKS_SUM
# define SL_EKS_SUM   2
#endif

#ifndef SL_EKS_AVG
# define SL_EKS_AVG   3
#endif

#ifndef SL_EKS_STD
# define SL_EKS_STD   4
#endif

#ifndef SL_EKS_SIZE
# define SL_EKS_SIZE  5
#endif


#ifndef SL_SORTED_IN
# define SL_SORTED_IN   0x1LL
#endif

#ifndef SL_SORTED_OUT
# define SL_SORTED_OUT  0x2LL
#endif


#ifndef SL_MSEG_FM_EXACT
# define SL_MSEG_FM_EXACT         0
#endif
#ifndef SL_MSEG_FM_ALLORNOTHING
# define SL_MSEG_FM_ALLORNOTHING  1
#endif
#ifndef SL_MSEG_FM_MIDDLE
# define SL_MSEG_FM_MIDDLE        2
#endif


/* partition conditions, sl_type pepckeys__partcond2_t pepckeys_partcond2_t */
typedef struct pepckeys__partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} pepckeys_partcond2_t;


#ifndef SLPC_COUNTS_MM
# define SLPC_COUNTS_MM   0x1
#endif
#ifndef SLPC_COUNTS_LH
# define SLPC_COUNTS_LH   0x2
#endif
#ifndef SLPC_WEIGHTS_MM
# define SLPC_WEIGHTS_MM  0x4
#endif
#ifndef SLPC_WEIGHTS_LH
# define SLPC_WEIGHTS_LH  0x8
#endif

/* partition conditions, sl_type pepckeys__partcond_t pepckeys_partcond_t pepckeys_partcond_p */
typedef struct pepckeys__partcond_t
{
  pepckeys_slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} pepckeys_partcond_t, *pepckeys_partcond_p;


/* internal partition conditions, sl_type pepckeys__partcond_intern_t pepckeys_partcond_intern_t pepckeys_partcond_intern_p */
typedef struct pepckeys__partcond_intern_t
{
  pepckeys_slint_t pcm;
  pepckeys_slint_t count_min, count_max;
  pepckeys_slint_t count_low, count_high;
#ifdef elem_weight
  pepckeys_slweight_t weight_min, weight_max;
  pepckeys_slweight_t weight_low, weight_high;
#endif

} pepckeys_partcond_intern_t, *pepckeys_partcond_intern_p;


/* sl_type pepckeys__parttype_t pepckeys_parttype_t pepckeys_parttype_p */
typedef struct pepckeys__parttype_t
{
  pepckeys_slint_t type;

} pepckeys_parttype_t, *pepckeys_parttype_p;


/* generic binning method */

/* sl_type pepckeys__bin_t pepckeys_bin_t */
typedef struct pepckeys__bin_t
{
  pepckeys_elements_t s;

#ifdef elem_weight
  pepckeys_slweight_t weight;
#endif

} pepckeys_bin_t;


/* sl_type pepckeys__splitter_t pepckeys_splitter_t */
typedef struct pepckeys__splitter_t
{
  pepckeys_slint_t n;

  int *displs;
  pepckeys_slkey_pure_t *s;
  pepckeys_slint_t *sn;

} pepckeys_splitter_t;


struct pepckeys__binning_t;

/* sl_type pepckeys_binning_pre_f pepckeys_binning_exec_f pepckeys_binning_refine_f pepckeys_binning_hit_f pepckeys_binning_finalize_f pepckeys_binning_post_f */
typedef pepckeys_slint_t (*pepckeys_binning_pre_f)(struct pepckeys__binning_t *bm);
typedef pepckeys_slint_t (*pepckeys_binning_exec_f)(struct pepckeys__binning_t *bm, pepckeys_bin_t *bin, pepckeys_slcount_t *counts, pepckeys_slweight_t *weights);
typedef pepckeys_slint_t (*pepckeys_binning_refine_f)(struct pepckeys__binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t k, pepckeys_slcount_t *counts, pepckeys_slweight_t *weights, pepckeys_splitter_t *sp, pepckeys_slint_t s, pepckeys_bin_t *new_bin);
typedef pepckeys_slint_t (*pepckeys_binning_hit_f)(struct pepckeys__binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t k, pepckeys_slcount_t *counts, pepckeys_splitter_t *sp, pepckeys_slint_t s);
typedef pepckeys_slint_t (*pepckeys_binning_finalize_f)(struct pepckeys__binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t dc, pepckeys_slweight_t dw, pepckeys_slint_t lc_min, pepckeys_slint_t lc_max, pepckeys_slcount_t *lcs, pepckeys_slweight_t *lws, pepckeys_splitter_t *sp, pepckeys_slint_t s);
typedef pepckeys_slint_t (*pepckeys_binning_post_f)(struct pepckeys__binning_t *bm);


/* sl_type pepckeys__binning_data_t pepckeys_binning_data_t */
typedef union pepckeys__binning_data_t
{
  struct
  {
    pepckeys_slint_t rhigh, rlow, rwidth;
    pepckeys_slint_t rcurrent;
    pepckeys_slkey_pure_t bit_mask;

    pepckeys_elements_t sx;

  } radix;

} pepckeys_binning_data_t;


/* sl_type pepckeys__binning_t pepckeys_binning_t */
typedef struct pepckeys__binning_t
{
  pepckeys_slint_t nbins, max_nbins;
  
  pepckeys_binning_pre_f pre;
  pepckeys_binning_exec_f exec;
  pepckeys_binning_refine_f refine;
  pepckeys_binning_hit_f hit;
  pepckeys_binning_finalize_f finalize;
  pepckeys_binning_post_f post;

  pepckeys_slint_t sorted;

  pepckeys_slint_t docounts;
#ifdef elem_weight
  pepckeys_slint_t doweights;
#endif

  pepckeys_binning_data_t bd;

} pepckeys_binning_t;


/* sl_type pepckeys__local_bins_t pepckeys_local_bins_t */
typedef struct pepckeys__local_bins_t
{
  pepckeys_binning_t *bm;

  pepckeys_slint_t nbins, max_nbins;
  pepckeys_slint_t nelements;

  pepckeys_slint_t docounts;
#ifdef elem_weight
  pepckeys_slint_t doweights;
#endif

  pepckeys_slint_t nbinnings, max_nbinnings;

  pepckeys_slint_t nbins_new, last_new_b, last_new_k;
  pepckeys_bin_t *bins, *bins_new;
  pepckeys_bin_t *bins0, *bins1;

  pepckeys_slint_t *bcws;

#if defined(elem_weight) && defined(pepckeys_sl_weight_intequiv)
  pepckeys_slint_t cw_factor, w_index, bin_cw_factor;
  pepckeys_slweight_t *cws, *bin_cws;
  pepckeys_slweight_t *prefix_cws;
#else
  pepckeys_slint_t *cs, *bin_cs;
  pepckeys_slint_t *prefix_cs;
# ifdef elem_weight
  pepckeys_slweight_t *ws, *bin_ws;
  pepckeys_slweight_t *prefix_ws;
# endif
#endif

  pepckeys_slint_t last_exec_b;

} pepckeys_local_bins_t;


/* sl_type pepckeys__global_bins_t pepckeys_global_bins_t */
typedef struct pepckeys__global_bins_t
{
  pepckeys_binning_t *bm;
  
  pepckeys_local_bins_t lb;

  pepckeys_slint_t *bcws;

#if defined(elem_weight) && defined(pepckeys_sl_weight_intequiv)
  pepckeys_slweight_t *cws;
  pepckeys_slweight_t *prefix_cws;
#else
  pepckeys_slint_t *cs;
  pepckeys_slint_t *prefix_cs;
# ifdef elem_weight
  pepckeys_slweight_t *ws;
  pepckeys_slweight_t *prefix_ws;
# endif
#endif

} pepckeys_global_bins_t;


/* sl_type pepckeys_rti_cmc_t */
typedef struct
{
  pepckeys_slint_t cmp, movek, moved;

} pepckeys_rti_cmc_t;

#ifndef my_rti_ccmp
# define my_rti_ccmp(m)    m.cmc.cmp
# define my_rti_cmovek(m)  m.cmc.movek
# define my_rti_cmoved(m)  m.cmc.moved
#endif


/* sl_type pepckeys_rti_tim_t */
typedef struct
{
  double start, stop;
  double last, cumu;

  pepckeys_slint_t num;

} pepckeys_rti_tim_t[rti_tids];

#ifndef my_rti_tlast
# define my_rti_tlast(m, t)  m.tim[t].last
# define my_rti_tcumu(m, t)  m.tim[t].cumu
# define my_rti_tnum(m, t)   m.tim[t].num
#endif


/* sl_type pepckeys_rti_mem_t */
typedef struct
{
  pepckeys_slint_t nalloc, nfree;
  pepckeys_slint_t max, cur, cur_max;

} pepckeys_rti_mem_t;


/* sl_type pepckeys_rti_t */
typedef struct
{
  /* compare-move-counter */
  pepckeys_rti_cmc_t cmc;
  /* timer */
  pepckeys_rti_tim_t tim;
  /* memory */
  pepckeys_rti_mem_t mem;

} pepckeys_rti_t;

#ifndef my_rti_reset
# define my_rti_reset(m)  memset((void *) &m, 0, sizeof(m))
#endif


/* sl_type pepckeys__sl_context_t pepckeys_sl_context_t */
typedef struct pepckeys__sl_context_t
{

/* src/base/base.c */
  struct {
int dummy_rank;
  } sl;
#ifdef pepckeys_SL_USE_RTI
pepckeys_rti_t rti;
#endif
  struct {
pepckeys_slint_t ip_threshold;
pepckeys_slint_t db_threshold;
pepckeys_slint_t ma_threshold;
  } sr;
  struct {
pepckeys_slint_t threshold;
  } sri;
/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
  struct {
MPI_Datatype int_datatype;
MPI_Datatype key_datatype;
MPI_Datatype pkey_datatype;
MPI_Datatype pwkey_datatype;
MPI_Datatype index_datatype;
MPI_Datatype weight_datatype;
MPI_Datatype data_datatype[SL_DATA_NMAX + 1];
int rank;
  } mpi;
#endif
#ifdef SL_USE_MPI
  struct {
void *sendrecv_replace_mem;
pepckeys_slint_t sendrecv_replace_memsize;
pepckeys_slint_t sendrecv_replace_mpi_maxsize;
  } me;
#endif
#ifdef SL_USE_MPI
  struct {
double t[2];
pepckeys_slint_t max_nprocs;
pepckeys_slint_t packed;
pepckeys_slint_t minalloc;
double overalloc;
  } meas;
#endif
#ifdef SL_USE_MPI
  struct {
pepckeys_slint_t packed;
pepckeys_slint_t db_packed;
pepckeys_slint_t ip_packed;
  } mea;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef pepckeys_MSEG_ROOT
int root;
#endif
#ifdef pepckeys_MSEG_BORDER_UPDATE_REDUCTION
double border_update_count_reduction;
double border_update_weight_reduction;
#endif
#ifdef pepckeys_MSEG_FORWARD_ONLY
pepckeys_slint_t forward_only;
#endif
#ifdef pepckeys_MSEG_INFO
pepckeys_slint_t info_rounds;
pepckeys_slint_t *info_finish_rounds;
double info_finish_rounds_avg;
pepckeys_slweight_t info_total_weights;
#endif
pepckeys_slint_t binnings;
pepckeys_slint_t finalize_mode;
  } mseg;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef pepckeys_MSS_ROOT
int root;
#endif
  } mss;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
pepckeys_slint_t sync;
  } msm;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
pepckeys_slint_t sync;
pepckeys_partcond_t *r_pc;
  } msp;
#endif
#ifdef SL_USE_MPI
  struct {
double i_t[3];
double p_t[3];
double b_t[3];
pepckeys_slint_t sync;
pepckeys_slint_t i_sync;
pepckeys_slint_t p_sync;
pepckeys_slint_t b_sync;
pepckeys_slint_t back_packed;
  } mssp;
#endif
} pepckeys_sl_context_t;






/* sl_macro pepckeys_elem_set_size pepckeys_elem_set_max_size pepckeys_elem_set_keys pepckeys_elem_set_indices */
#define pepckeys_elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pepckeys_elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pepckeys_elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define pepckeys_elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))
/* DATAX_TEMPLATE_BEGIN */
#define pepckeys_elem_set_data0(_e_, _b_)     ((_e_)->data0 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data1(_e_, _b_)     ((_e_)->data1 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data2(_e_, _b_)     ((_e_)->data2 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data3(_e_, _b_)     ((_e_)->data3 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data4(_e_, _b_)     ((_e_)->data4 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data5(_e_, _b_)     ((_e_)->data5 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data6(_e_, _b_)     ((_e_)->data6 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data7(_e_, _b_)     ((_e_)->data7 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data8(_e_, _b_)     ((_e_)->data8 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data9(_e_, _b_)     ((_e_)->data9 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data10(_e_, _b_)     ((_e_)->data10 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data11(_e_, _b_)     ((_e_)->data11 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data12(_e_, _b_)     ((_e_)->data12 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data13(_e_, _b_)     ((_e_)->data13 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data14(_e_, _b_)     ((_e_)->data14 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data15(_e_, _b_)     ((_e_)->data15 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data16(_e_, _b_)     ((_e_)->data16 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data17(_e_, _b_)     ((_e_)->data17 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data18(_e_, _b_)     ((_e_)->data18 = (_b_))  /* sl_macro */
#define pepckeys_elem_set_data19(_e_, _b_)     ((_e_)->data19 = (_b_))  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro pepckeys_elem_get_size pepckeys_elem_get_max_size pepckeys_elem_get_keys pepckeys_elem_get_indices */
#define pepckeys_elem_get_size(_e_)           (_e_)->size
#define pepckeys_elem_get_max_size(_e_)       (_e_)->max_size
#define pepckeys_elem_get_keys(_e_)           (_e_)->keys
#define pepckeys_elem_get_indices(_e_)        (_e_)->indices
/* DATAX_TEMPLATE_BEGIN */
#define pepckeys_elem_get_data0(_e_)          (_e_)->data0  /* sl_macro */
#define pepckeys_elem_get_data1(_e_)          (_e_)->data1  /* sl_macro */
#define pepckeys_elem_get_data2(_e_)          (_e_)->data2  /* sl_macro */
#define pepckeys_elem_get_data3(_e_)          (_e_)->data3  /* sl_macro */
#define pepckeys_elem_get_data4(_e_)          (_e_)->data4  /* sl_macro */
#define pepckeys_elem_get_data5(_e_)          (_e_)->data5  /* sl_macro */
#define pepckeys_elem_get_data6(_e_)          (_e_)->data6  /* sl_macro */
#define pepckeys_elem_get_data7(_e_)          (_e_)->data7  /* sl_macro */
#define pepckeys_elem_get_data8(_e_)          (_e_)->data8  /* sl_macro */
#define pepckeys_elem_get_data9(_e_)          (_e_)->data9  /* sl_macro */
#define pepckeys_elem_get_data10(_e_)          (_e_)->data10  /* sl_macro */
#define pepckeys_elem_get_data11(_e_)          (_e_)->data11  /* sl_macro */
#define pepckeys_elem_get_data12(_e_)          (_e_)->data12  /* sl_macro */
#define pepckeys_elem_get_data13(_e_)          (_e_)->data13  /* sl_macro */
#define pepckeys_elem_get_data14(_e_)          (_e_)->data14  /* sl_macro */
#define pepckeys_elem_get_data15(_e_)          (_e_)->data15  /* sl_macro */
#define pepckeys_elem_get_data16(_e_)          (_e_)->data16  /* sl_macro */
#define pepckeys_elem_get_data17(_e_)          (_e_)->data17  /* sl_macro */
#define pepckeys_elem_get_data18(_e_)          (_e_)->data18  /* sl_macro */
#define pepckeys_elem_get_data19(_e_)          (_e_)->data19  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro pepckeys_elem_set_block pepckeys_elem_set_block_size pepckeys_elem_get_block pepckeys_elem_get_block_size */
#define pepckeys_elem_set_block(_e_, _b_)       ((_e_)->keys = (_b_), (_e_)->max_size = -1)
#define pepckeys_elem_set_block_size(_e_, _s_)  ((_e_)->size = (_s_))
#define pepckeys_elem_get_block(_e_)            ((void *) (((_e_)->max_size < 0)?(_e_)->keys:NULL))
#define pepckeys_elem_get_block_size(_e_)       (_e_)->size

/* sl_macro pepckeys_pelem_set_size pepckeys_pelem_set_max_size pepckeys_pelem_set_elements */
#define pepckeys_pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pepckeys_pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pepckeys_pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))

/* sl_macro pepckeys_pelem_get_size pepckeys_pelem_get_max_size pepckeys_pelem_get_elements */
#define pepckeys_pelem_get_size(_e_)           (_e_)->size
#define pepckeys_pelem_get_max_size(_e_)       (_e_)->max_size
#define pepckeys_pelem_get_elements(_e_)       (_e_)->elements

/* sl_macro pepckeys_SL_DEFCON */
#define pepckeys_SL_DEFCON(_v_)  (pepckeys_sl_default_context._v_)






/* src/base/base.c */
extern pepckeys_sl_context_t pepckeys_sl_default_context;
extern const int pepckeys_default_sl_dummy_rank;
#ifdef pepckeys_SL_USE_RTI
extern const pepckeys_rti_t pepckeys_default_rti;
#endif
extern const pepckeys_slint_t pepckeys_default_sr_ip_threshold;
extern const pepckeys_slint_t pepckeys_default_sr_db_threshold;
extern const pepckeys_slint_t pepckeys_default_sr_ma_threshold;
extern const pepckeys_slint_t pepckeys_default_sri_threshold;

/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
extern const MPI_Datatype pepckeys_default_mpi_int_datatype;
extern const MPI_Datatype pepckeys_default_mpi_key_datatype;
extern const MPI_Datatype pepckeys_default_mpi_pkey_datatype;
extern const MPI_Datatype pepckeys_default_mpi_pwkey_datatype;
extern const MPI_Datatype pepckeys_default_mpi_index_datatype;
extern const MPI_Datatype pepckeys_default_mpi_weight_datatype;
extern const MPI_Datatype pepckeys_default_mpi_data_datatype[];
extern const int pepckeys_default_mpi_rank;
#endif
extern const void *pepckeys_default_me_sendrecv_replace_mem;
extern const pepckeys_slint_t pepckeys_default_me_sendrecv_replace_memsize;
extern const pepckeys_slint_t pepckeys_default_me_sendrecv_replace_mpi_maxsize;
extern const double pepckeys_default_meas_t[];
extern const pepckeys_slint_t pepckeys_default_meas_max_nprocs;
extern const pepckeys_slint_t pepckeys_default_meas_packed;
extern const pepckeys_slint_t pepckeys_default_meas_minalloc;
extern const double pepckeys_default_meas_overalloc;
extern const pepckeys_slint_t pepckeys_default_mea_packed;
extern const pepckeys_slint_t pepckeys_default_mea_db_packed;
extern const pepckeys_slint_t pepckeys_default_mea_ip_packed;
#ifdef pepckeys_MSEG_ROOT
extern const int pepckeys_default_mseg_root;
#endif
#ifdef pepckeys_MSEG_BORDER_UPDATE_REDUCTION
extern const double pepckeys_default_mseg_border_update_count_reduction;
extern const double pepckeys_default_mseg_border_update_weight_reduction;
#endif
#ifdef pepckeys_MSEG_FORWARD_ONLY
extern const pepckeys_slint_t pepckeys_default_mseg_forward_only;
#endif
#ifdef pepckeys_MSEG_INFO
extern const pepckeys_slint_t pepckeys_default_mseg_info_rounds;
extern const pepckeys_slint_t *pepckeys_default_mseg_info_finish_rounds;
extern const double pepckeys_default_mseg_info_finish_rounds_avg;
extern const pepckeys_slweight_t pepckeys_default_mseg_info_total_weights;
#endif
extern const pepckeys_slint_t pepckeys_default_mseg_binnings;
extern const pepckeys_slint_t pepckeys_default_mseg_finalize_mode;
#ifdef pepckeys_MSS_ROOT
extern const int pepckeys_default_mss_root;
#endif
extern const double pepckeys_default_msm_t[];
extern const pepckeys_slint_t pepckeys_default_msm_sync;
extern const double pepckeys_default_msp_t[];
extern const pepckeys_slint_t pepckeys_default_msp_sync;
extern const pepckeys_partcond_t *pepckeys_default_msp_r_pc;
extern const double pepckeys_default_mssp_i_t[];
extern const double pepckeys_default_mssp_p_t[];
extern const double pepckeys_default_mssp_b_t[];
extern const pepckeys_slint_t pepckeys_default_mssp_sync;
extern const pepckeys_slint_t pepckeys_default_mssp_i_sync;
extern const pepckeys_slint_t pepckeys_default_mssp_p_sync;
extern const pepckeys_slint_t pepckeys_default_mssp_b_sync;
extern const pepckeys_slint_t pepckeys_default_mssp_back_packed;






/* src/base/base.c */
pepckeys_slint_t SL_PROTO(pepckeys_binning_create)(pepckeys_local_bins_t *lb, pepckeys_slint_t max_nbins, pepckeys_slint_t max_nbinnings, pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t docounts, pepckeys_slint_t doweights, pepckeys_binning_t *bm);
pepckeys_slint_t SL_PROTO(pepckeys_binning_destroy)(pepckeys_local_bins_t *lb);
pepckeys_slint_t SL_PROTO(pepckeys_binning_pre)(pepckeys_local_bins_t *lb);
pepckeys_slint_t SL_PROTO(pepckeys_binning_exec_reset)(pepckeys_local_bins_t *lb, pepckeys_slint_t do_bins, pepckeys_slint_t do_prefixes);
pepckeys_slint_t SL_PROTO(pepckeys_binning_exec)(pepckeys_local_bins_t *lb, pepckeys_slint_t b, pepckeys_slint_t do_bins, pepckeys_slint_t do_prefixes);
pepckeys_slint_t SL_PROTO(pepckeys_binning_refine)(pepckeys_local_bins_t *lb, pepckeys_slint_t b, pepckeys_slint_t k, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_hit)(pepckeys_local_bins_t *lb, pepckeys_slint_t b, pepckeys_slint_t k, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_finalize)(pepckeys_local_bins_t *lb, pepckeys_slint_t b, pepckeys_slint_t dc, pepckeys_slweight_t dw, pepckeys_slint_t lc_min, pepckeys_slint_t lc_max, pepckeys_slcount_t *lcs, pepckeys_slweight_t *lws, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_post)(pepckeys_local_bins_t *lb);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_create)(pepckeys_binning_t *bm, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, pepckeys_slint_t sorted);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_destroy)(pepckeys_binning_t *bm);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_pre)(pepckeys_binning_t *bm);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_exec)(pepckeys_binning_t *bm, pepckeys_bin_t *bin, pepckeys_slcount_t *counts, pepckeys_slweight_t *weights);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_refine)(pepckeys_binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t k, pepckeys_slcount_t *counts, pepckeys_slweight_t *weights, pepckeys_splitter_t *sp, pepckeys_slint_t s, pepckeys_bin_t *new_bin);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_hit)(pepckeys_binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t k, pepckeys_slcount_t *counts, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_finalize)(pepckeys_binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t dc, pepckeys_slweight_t dw, pepckeys_slint_t lc_min, pepckeys_slint_t lc_max, pepckeys_slcount_t *lcs, pepckeys_slweight_t *lws, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_post)(pepckeys_binning_t *bm);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components);
pepckeys_slint_t SL_PROTO(pepckeys_elements_free)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_realloc)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloca)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components);
pepckeys_slint_t SL_PROTO(pepckeys_elements_freea)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_from_blocks)(pepckeys_elements_t *s, pepckeys_slint_t nblocks, void **blocks, pepckeys_slint_t *blocksizes, pepckeys_slint_t alignment, pepckeys_slint_t nmax, slcint_t components);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_from_block)(pepckeys_elements_t *s, void *block, pepckeys_slint_t blocksize, pepckeys_slint_t alignment, pepckeys_slint_t nmax, slcint_t components);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_block)(pepckeys_elements_t *s, void **block, pepckeys_slint_t *blocksize, pepckeys_slint_t alignment, pepckeys_slint_t maxblocksize);
pepckeys_slint_t SL_PROTO(pepckeys_elements_copy)(pepckeys_elements_t *s, pepckeys_elements_t *d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_copy_at)(pepckeys_elements_t *s, pepckeys_slint_t sat, pepckeys_elements_t *d, pepckeys_slint_t dat);
pepckeys_slint_t SL_PROTO(pepckeys_elements_ncopy)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_elements_nmove)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_elements_printf)(pepckeys_elements_t *s, const char *prefix);
pepckeys_slint_t SL_PROTO(pepckeys_elements_extract)(pepckeys_elements_t *src, pepckeys_slint_t nelements, pepckeys_elements_t *dst0, pepckeys_elements_t *dst1);
pepckeys_slint_t SL_PROTO(pepckeys_elements_touch)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_digest_sum)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(pepckeys_elements_crc32)(pepckeys_elements_t *s, pepckeys_slint nelements, pepckeys_slint_t keys, pepckeys_slint_t data);
pepckeys_slint_t SL_PROTO(pepckeys_elements_digest_hash)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components, void *hash);
pepckeys_slint_t SL_PROTO(pepckeys_elements_random_exchange)(pepckeys_elements_t *s, pepckeys_slint_t rounds, pepckeys_elements_t *xs);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_init_seed)(unsigned long s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_init)(pepckeys_elements_t *s, pepckeys_keys_init_type_t t, pepckeys_keys_init_data_t d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_init_randomized)(pepckeys_elements_t *s, pepckeys_slint_t nkeys, pepckeys_keys_init_type_t t, pepckeys_keys_init_data_t d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_init_from_file)(pepckeys_elements_t *s, pepckeys_slint_t data, char *filename, pepckeys_slint_t from, pepckeys_slint_t to, pepckeys_slint_t const_bytes_per_line);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_save_to_file)(pepckeys_elements_t *s, char *filename);
pepckeys_slint_t SL_PROTO(pepckeys_elements_validate_order)(pepckeys_elements_t *s, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_elements_validate_order_bmask)(pepckeys_elements_t *s, pepckeys_slint_t n, pepckeys_slkey_pure_t bmask);
pepckeys_slint_t SL_PROTO(pepckeys_elements_validate_order_weight)(pepckeys_elements_t *s, pepckeys_slint_t n, pepckeys_slkey_pure_t weight);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_stats)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *stats);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_stats_print)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_print_keys)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_print_all)(pepckeys_elements_t *s);
pepckeys_slweight_t SL_PROTO(pepckeys_elements_get_weight)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_get_minmax_keys)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slkey_pure_t *minmaxkeys);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_packed)(pepckeys_packed_elements_t *s, pepckeys_slint_t nelements);
pepckeys_slint_t SL_PROTO(pepckeys_elements_free_packed)(pepckeys_packed_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_packed_from_block)(pepckeys_packed_elements_t *s, void *block, pepckeys_slint_t blocksize, pepckeys_slint_t alignment, pepckeys_slint_t nmax);
pepckeys_slint_t SL_PROTO(pepckeys_elements_pack_indexed)(pepckeys_elements_t *s, pepckeys_packed_elements_t *d, pepckeys_slindex_t *rindx, pepckeys_slindex_t *windx);
pepckeys_slint_t SL_PROTO(pepckeys_elements_pack)(pepckeys_elements_t *s, pepckeys_packed_elements_t *d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_pack_at)(pepckeys_elements_t *s, pepckeys_slint_t sat, pepckeys_packed_elements_t *d, pepckeys_slint_t dat);
pepckeys_slint_t SL_PROTO(pepckeys_elements_unpack_indexed)(pepckeys_packed_elements_t *s, pepckeys_elements_t *d, pepckeys_slindex_t *rindx, pepckeys_slindex_t *windx);
pepckeys_slint_t SL_PROTO(pepckeys_elements_unpack)(pepckeys_packed_elements_t *s, pepckeys_elements_t *d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_unpack_at)(pepckeys_packed_elements_t *s, pepckeys_slint_t sat, pepckeys_elements_t *d, pepckeys_slint_t dat);
pepckeys_slint_t SL_PROTO(pepckeys_elements_unpack_keys)(pepckeys_packed_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_auto_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_m2x_func _x0_1, pepckeys_m2x_func _0x_1);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_01_X)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t, pepckeys_m2X_func _X0_1, pepckeys_m2X_func _0X_1);
pepckeys_slint SL_PROTO(pepckeys_merge2_simplify_s1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_slint s1elements);
pepckeys_slint SL_PROTO(pepckeys_merge2_memory_adaptive)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_compo_hula)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_sseq_x0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_sseq_0x_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_sseq_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_sseq_01)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_sbin_x0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_sbin_0x_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_sbin_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_sbin_01)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_shyb_x0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_shyb_0x_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_shyb_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_shyb_01)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_straight_x0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_straight_0x_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_straight_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_straight_x_0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_straight_X0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_straight_0X_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_straight_01_X)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_basic_straight_X0_1u)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_merge2_compo_tridgell)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_2way_ip_int)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t p, int *displs, pepckeys_merge2x_f m2x);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_2way_ip_int_rec)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t p, int *displs, pepckeys_merge2x_f m2x);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_int)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, int *displs, int *counts);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_int_idx)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, int *displs, int *counts);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_idx)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, pepckeys_slindex_t *displs, pepckeys_slindex_t *counts);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_unpack_idx)(pepckeys_packed_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, pepckeys_slindex_t *displs, pepckeys_slindex_t *counts);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_unpack_idxonly)(pepckeys_packed_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, pepckeys_slindex_t *displs, pepckeys_slindex_t *counts);
pepckeys_slint_t SL_PROTO(pepckeys_permute_generic_db)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_permute_generic_t *pg, void *pg_data);
pepckeys_slint_t SL_PROTO(pepckeys_permute_generic_ip)(pepckeys_elements_t *s, pepckeys_elements_t *x, pepckeys_permute_generic_t *pg, void *pg_data);
pepckeys_slint SL_PROTO(pepckeys_sl_search_sequential_lt)(pepckeys_elements_t *s, pepckeys_slpkey_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_sequential_le)(pepckeys_elements_t *s, pepckeys_slpkey_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_sequential_gt)(pepckeys_elements_t *s, pepckeys_slpkey_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_sequential_ge)(pepckeys_elements_t *s, pepckeys_slpkey_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_sequential_lt)(pepckeys_elements_t *s, pepckeys_slpkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_sequential_le)(pepckeys_elements_t *s, pepckeys_slpkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_sequential_gt)(pepckeys_elements_t *s, pepckeys_slpkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_sequential_ge)(pepckeys_elements_t *s, pepckeys_slpkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_lt)(pepckeys_elements_t *s, pepckeys_slpkey_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_le)(pepckeys_elements_t *s, pepckeys_slpkey_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_gt)(pepckeys_elements_t *s, pepckeys_slpkey_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_ge)(pepckeys_elements_t *s, pepckeys_slpkey_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_binary_lt)(pepckeys_elements_t *s, pepckeys_slpkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_binary_le)(pepckeys_elements_t *s, pepckeys_slpkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_binary_gt)(pepckeys_elements_t *s, pepckeys_slpkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_binary_ge)(pepckeys_elements_t *s, pepckeys_slpkey_t *k);
pepckeys_slint_t SL_PROTO(pepckeys_sl_search_binary_lt_bmask)(pepckeys_elements_t *s, pepckeys_slpkey_t k, pepckeys_slpkey_t bmask);
pepckeys_slint_t SL_PROTO(pepckeys_sl_search_binary_le_bmask)(pepckeys_elements_t *s, pepckeys_slpkey_t k, pepckeys_slpkey_t bmask);
pepckeys_slint_t SL_PROTO(pepckeys_sl_search_binary_sign_switch)(pepckeys_elements_t *s);
pepckeys_slint SL_PROTO(pepckeys_sl_search_hybrid_lt)(pepckeys_elements_t *s, pepckeys_slpkey_t k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_hybrid_le)(pepckeys_elements_t *s, pepckeys_slpkey_t k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_hybrid_gt)(pepckeys_elements_t *s, pepckeys_slpkey_t k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_hybrid_ge)(pepckeys_elements_t *s, pepckeys_slpkey_t k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_hybrid_lt)(pepckeys_elements_t *s, pepckeys_slpkey_t *k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_hybrid_le)(pepckeys_elements_t *s, pepckeys_slpkey_t *k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_hybrid_gt)(pepckeys_elements_t *s, pepckeys_slpkey_t *k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_p_hybrid_ge)(pepckeys_elements_t *s, pepckeys_slpkey_t *k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_ilog2c)(pepckeys_slint x);
pepckeys_slint SL_PROTO(pepckeys_ilog2f)(pepckeys_slint x);
pepckeys_slint SL_PROTO(pepckeys_print_bits)(pepckeys_slint v);
pepckeys_slint SL_PROTO(pepckeys_pivot_random)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_counts2displs)(pepckeys_slint_t n, int *counts, int *displs);
pepckeys_slint_t SL_PROTO(pepckeys_displs2counts)(pepckeys_slint_t n, int *displs, int *counts, pepckeys_slint_t total_counts);
void SL_PROTO(pepckeys_get_displcounts_extent)(pepckeys_slint_t n, int *displs, int *counts, pepckeys_slint_t *lb, pepckeys_slint_t *extent);
void SL_PROTO(pepckeys_elem_set_data)(pepckeys_elements_t *e, ...);
pepckeys_slint_t SL_PROTO(pepckeys_elem_get_max_byte)();
pepckeys_slint_t SL_PROTO(pepckeys_elem_reverse)(pepckeys_elements_t *e, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_nxchange_at)(pepckeys_elements_t *e0, pepckeys_slint_t at0, pepckeys_elements_t *e1, pepckeys_slint_t at1, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_nxchange)(pepckeys_elements_t *e0, pepckeys_elements_t *e1, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_nxchange_ro0)(pepckeys_elements_t *e0, pepckeys_elements_t *e1, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_rotate)(pepckeys_elements_t *e, pepckeys_slint_t m, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_rotate_ro0)(pepckeys_elements_t *e, pepckeys_slint_t m, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_rotate_ro1)(pepckeys_elements_t *e, pepckeys_slint_t m, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_sort_counting_use_displs)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ndispls, pepckeys_slint_t *displs);
pepckeys_slint_t SL_PROTO(pepckeys_sort_counting_use_counts)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ncounts, pepckeys_slint_t *counts);
pepckeys_slint_t SL_PROTO(pepckeys_sort_counting_get_counts)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ncounts, pepckeys_slint_t *counts);
pepckeys_slint_t SL_PROTO(pepckeys_sort_counting)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ncounts);
pepckeys_slint SL_PROTO(pepckeys_sort_heap)(pepckeys_elements_t *s, pepckeys_elements_t *xs);
pepckeys_slint_t SL_PROTO(pepckeys_sort_insert_bmask_kernel)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slkey_pure_t bmask);
pepckeys_slint_t SL_PROTO(pepckeys_sort_insert)(pepckeys_elements_t *s, pepckeys_elements_t *sx);
pepckeys_slint_t SL_PROTO(pepckeys_sort_permute_forward)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t *perm, pepckeys_slint_t offset, pepckeys_slint_t mask_bit);
pepckeys_slint_t SL_PROTO(pepckeys_sort_permute_backward)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t *perm, pepckeys_slint_t offset, pepckeys_slint_t mask_bit);
pepckeys_slint SL_PROTO(pepckeys_sort_quick)(pepckeys_elements_t *s, pepckeys_elements_t *xs);
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix_ip)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth);
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix_db)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth);
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix_ma)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth);
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth);
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix_1bit_kernel)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow);
pepckeys_slint SL_PROTO(pepckeys_sort_radix_1bit)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow);
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix_iter)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t presorted, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth);
pepckeys_slint SL_PROTO(pepckeys_sn_hypercube_lh)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_hypercube_hl)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_odd_even_trans)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_odd)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_even)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_batcher)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_bitonic)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_connected)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint_t SL_PROTO(pepckeys_split_generic_db)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_split_generic_t *sg, void *sg_data, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_split_generic_ip)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_split_generic_t *sg, void *sg_data, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_split_generic_count_db)(pepckeys_elements_t *s, pepckeys_split_generic_t *sg, void *sg_data, int *counts, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_split_generic_count_ip)(pepckeys_elements_t *s, pepckeys_split_generic_t *sg, void *sg_data, int *counts, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_split_generic_rearrange_db)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_split_generic_t *sg, void *sg_data, int *counts, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_split_generic_rearrange_ip)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_split_generic_t *sg, void *sg_data, int *counts, int *displs, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_splitter_reset)(pepckeys_splitter_t *sp);
pepckeys_slint_t SL_PROTO(pepckeys_splitx_radix)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t nclasses, pepckeys_slint_t shl, pepckeys_slint_t *counts);
pepckeys_slint SL_PROTO(pepckeys_split2_lt_ge)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t);
pepckeys_slint SL_PROTO(pepckeys_split2_le_gt)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t);
pepckeys_slint SL_PROTO(pepckeys_split3_lt_eq_gt)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t, pepckeys_slint *nlt, pepckeys_slint *nle);
pepckeys_slint SL_PROTO(pepckeys_split3_lt_eq_gt_old)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t, pepckeys_slint *nlt, pepckeys_slint *nle);
pepckeys_slint SL_PROTO(pepckeys_split2_b)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slkey_pure_t bmask);
pepckeys_slint SL_PROTO(pepckeys_splitk_k2c_af)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint k, pepckeys_slint *c, pepckeys_k2c_func k2c, void *k2c_data);
pepckeys_slint SL_PROTO(pepckeys_splitk_k2c)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint k, pepckeys_slint *c, pepckeys_k2c_func k2c, void *k2c_data);
pepckeys_slint SL_PROTO(pepckeys_splitk_k2c_count)(pepckeys_elements_t *s, pepckeys_slint k, pepckeys_slint *c, pepckeys_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI





/* src/base_mpi/base_mpi.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_create)(pepckeys_global_bins_t *gb, pepckeys_slint_t max_nbins, pepckeys_slint_t max_nbinnings, pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t docounts, pepckeys_slint_t doweights, pepckeys_binning_t *bm, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_destroy)(pepckeys_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_pre)(pepckeys_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_exec_reset)(pepckeys_global_bins_t *gb, pepckeys_slint_t do_bins, pepckeys_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_exec_local)(pepckeys_global_bins_t *gb, pepckeys_slint_t b, pepckeys_slint_t do_bins, pepckeys_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_exec_global)(pepckeys_global_bins_t *gb, pepckeys_slint_t do_bins, pepckeys_slint_t do_prefixes, pepckeys_slint_t root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_refine)(pepckeys_global_bins_t *gb, pepckeys_slint_t b, pepckeys_slint_t k, pepckeys_splitter_t *sp, pepckeys_slint_t s, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_hit)(pepckeys_global_bins_t *gb, pepckeys_slint_t b, pepckeys_slint_t k, pepckeys_splitter_t *sp, pepckeys_slint_t s, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_finalize)(pepckeys_global_bins_t *gb, pepckeys_slint_t b, pepckeys_slint_t dc, pepckeys_slweight_t dw, pepckeys_slint_t lc_min, pepckeys_slint_t lc_max, pepckeys_slcount_t *lcs, pepckeys_slweight_t *lws, pepckeys_splitter_t *sp, pepckeys_slint_t s, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_post)(pepckeys_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_datatypes_init)();
pepckeys_slint_t SL_PROTO(pepckeys_mpi_datatypes_release)();
pepckeys_slint_t SL_PROTO(pepckeys_mpi_get_grid_properties)(pepckeys_slint_t ndims, pepckeys_slint_t *dims, pepckeys_slint_t *pos, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_subgroups_create)(pepckeys_slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_subgroups_delete)(pepckeys_slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);
int SL_PROTO(pepckeys_sl_MPI_Allreduce)(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank);
int SL_PROTO(pepckeys_sl_MPI_Alltoall_int)(void *sendbuf, int sendcount, void *recvbuf, int recvcount, MPI_Comm comm, int size, int rank);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_keys_init_from_file)(pepckeys_elements_t *s, char *filename, pepckeys_slint from, pepckeys_slint to, pepckeys_slint const_bytes_per_line, pepckeys_slint root, int size, int rank, MPI_Comm comm);
pepckeys_slint SL_PROTO(pepckeys_mpi_elements_validate_order)(pepckeys_elements_t *s, pepckeys_slint n, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_linear_exchange_pure_keys)(pepckeys_slkey_pure_t *in, pepckeys_slkey_pure_t *out, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_check_order)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t *orders, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_check_global_order)(pepckeys_slkey_pure_t local_min, pepckeys_slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_digest_sum)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(pepckeys_mpi_elements_crc32)(pepckeys_elements_t *s, pepckeys_slint_t n, pepckeys_slint_t keys, pepckeys_slint_t data, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_digest_hash)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_get_counts)(pepckeys_elements_t *s, pepckeys_slint_t *clocal, pepckeys_slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
pepckeys_slweight_t SL_PROTO(pepckeys_mpi_elements_get_weights)(pepckeys_elements_t *s, pepckeys_slweight_t *wlocal, pepckeys_slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_get_counts_and_weights)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t *counts, pepckeys_slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_sendrecv_replace)(pepckeys_elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_tproc_create_tproc)(pepckeys_tproc_t *tproc, pepckeys_tproc_f *tfn, pepckeys_tproc_reset_f *rfn, pepckeys_tproc_exdef exdef);
pepckeys_slint_t SL_PROTO(pepckeys_tproc_create_tproc_mod)(pepckeys_tproc_t *tproc, pepckeys_tproc_mod_f *tfn, pepckeys_tproc_reset_f *rfn, pepckeys_tproc_exdef exdef);
pepckeys_slint_t SL_PROTO(pepckeys_tproc_create_tprocs)(pepckeys_tproc_t *tproc, pepckeys_tprocs_f *tfn, pepckeys_tproc_reset_f *rfn, pepckeys_tproc_exdef exdef);
pepckeys_slint_t SL_PROTO(pepckeys_tproc_create_tprocs_mod)(pepckeys_tproc_t *tproc, pepckeys_tprocs_mod_f *tfn, pepckeys_tproc_reset_f *rfn, pepckeys_tproc_exdef exdef);
pepckeys_slint_t SL_PROTO(pepckeys_tproc_free)(pepckeys_tproc_t *tproc);
pepckeys_slint_t SL_PROTO(pepckeys_tproc_set_proclists)(pepckeys_tproc_t *tproc, pepckeys_slint_t nsend_procs, int *send_procs, pepckeys_slint_t nrecv_procs, int *recv_procs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_tproc_verify)(pepckeys_tproc_t tproc, void *data, pepckeys_elements_t *s, int proc);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoall_specific)(pepckeys_elements_t *sin, pepckeys_elements_t *sout, pepckeys_elements_t *xs, pepckeys_tproc_t tproc, void *data, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_db_packed)(pepckeys_elements_t *sbuf, int *scounts, int *sdispls, pepckeys_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_db)(pepckeys_elements_t *sbuf, int *scounts, int *sdispls, pepckeys_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_ip_packed)(pepckeys_elements_t *s, pepckeys_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_ip_double)(pepckeys_elements_t *s, pepckeys_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_ip_mpi)(pepckeys_elements_t *s, pepckeys_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_ip_dash)(pepckeys_elements_t *s, pepckeys_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_ip)(pepckeys_elements_t *s, pepckeys_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_proclists_db)(pepckeys_elements_t *sbuf, int *scounts, int *sdispls, int nsendprocs, int *sendprocs, pepckeys_elements_t *rbuf, int *rcounts, int *rdispls, int nrecvprocs, int *recvprocs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, pepckeys_slint_t structured);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_find_exact_equal)(pepckeys_elements_t *s, pepckeys_slint_t other_rank, pepckeys_slint_t high_rank, pepckeys_slint_t *ex_start, pepckeys_slint_t *ex_size, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_find_exact)(pepckeys_elements_t *s, pepckeys_slint_t other_rank, pepckeys_slint_t high_rank, pepckeys_slint_t *dst_size, pepckeys_slint_t *ex_start, pepckeys_slint_t *ex_sizes, pepckeys_slint_t *nx_move, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_merge2)(pepckeys_elements_t *s, pepckeys_slint_t other_rank, pepckeys_slint_t high_rank, pepckeys_slint_t *dst_size, pepckeys_merge2x_f m2, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_mergek_equal)(pepckeys_elements_t *s, pepckeys_sortnet_f sn, pepckeys_sortnet_data_t snd, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_mergek_sorted)(pepckeys_elements_t *s, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_mergek_sorted2)(pepckeys_elements_t *s, pepckeys_sortnet_f sn, pepckeys_sortnet_data_t snd, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_mergek)(pepckeys_elements_t *s, pepckeys_sortnet_f sn, pepckeys_sortnet_data_t snd, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_mergek_equal2)(pepckeys_elements_t *s, pepckeys_sortnet_f sn, pepckeys_sortnet_data_t snd, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_exact_generic)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, pepckeys_binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_exact_radix)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, pepckeys_slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_exact_radix_ngroups)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, pepckeys_slint_t ngroups, MPI_Comm *group_comms, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_exact_radix_2groups)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, MPI_Comm group_comm, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_sample_regular)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_rebalance)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_slint_t stable, pepckeys_slint_t *dst_size, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_rebalance_alltoallv)(pepckeys_elements_t *sbuf, int *scounts, int *sdispls, pepckeys_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
void SL_PROTO(pepckeys_mpi_partcond_set_even)(pepckeys_partcond_t *pcond, pepckeys_slint_t pcm, pepckeys_slint_t ntotal, double nimba, double wtotal, double wimba, int size, int rank);
pepckeys_slint_t SL_PROTO(pepckeys_init_partconds)(pepckeys_slint_t npconds, pepckeys_partcond_t *pconds, pepckeys_slint_t nparts, pepckeys_slint_t total_count, pepckeys_slweight_t total_weight);
pepckeys_slint_t SL_PROTO(pepckeys_init_partconds_intern)(pepckeys_slint_t npconds, pepckeys_partcond_intern_t *pci, pepckeys_partcond_t *pc, pepckeys_slint_t nparts, pepckeys_slint_t total_count, pepckeys_slweight_t total_weight);
pepckeys_slint_t SL_PROTO(pepckeys_merge_partconds)(pepckeys_partcond_t *pconds_in, pepckeys_slint_t npconds_in, pepckeys_partcond_t *pcond_out);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_gather_partconds_grouped)(pepckeys_partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, pepckeys_partcond_t *pconds_out, pepckeys_slint_t *npconds_out, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_gather_partconds)(pepckeys_partcond_t *pcond_in, pepckeys_partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_allgather_partconds)(pepckeys_partcond_t *pcond_in, pepckeys_partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_bcast_partconds)(pepckeys_slint_t npconds, pepckeys_partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_post_check_partconds)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_post_check_partconds_intern)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_stats)(pepckeys_elements_t *s, pepckeys_slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_generic_bulk)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, pepckeys_binning_t *bm, pepckeys_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_generic_grouped)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, pepckeys_binning_t *bm, pepckeys_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_generic)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, pepckeys_binning_t *bm, pepckeys_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_radix)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, pepckeys_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_radix_grouped)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, pepckeys_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_sample_regular)(pepckeys_elements_t *s, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, pepckeys_slint_t nsamples, pepckeys_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_merge)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_merge2)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t merge_type, pepckeys_slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_merge_radix)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t merge_type, pepckeys_slint_t sort_type, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t part_type, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition_radix)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t part_type, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition_exact_radix)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_partcond_t *pcond, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition_exact_radix_ngroups)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_partcond_t *pcond, pepckeys_slint_t ngroups, MPI_Comm *group_comms, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition_exact_radix_2groups)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_partcond_t *pcond, MPI_Comm group_comm, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_insert_radix)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slpkey_t *mmkeys, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_presorted_radix)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t merge_type, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_back)(pepckeys_elements_t *sin, pepckeys_elements_t *sout, pepckeys_elements_t *sx, pepckeys_slpkey_t *lh, pepckeys_slint_t ntotal, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, pepckeys_slint_t ytotal, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_xcounts2ycounts_grouped)(int *xcounts, pepckeys_slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_subxdispls2ycounts)(pepckeys_slint_t nsubs, int *sub_xdispls, pepckeys_slint_t *sub_sources, pepckeys_slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#undef SL_PROTO
#endif /* __SL_PEPCKEYS_H__ */
