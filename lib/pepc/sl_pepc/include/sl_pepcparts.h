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


#ifndef __SL_PEPCPARTS_H__
#define __SL_PEPCPARTS_H__

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif /* SL_USE_MPI */

#define SL_PROTO(_f_)  _f_

#include "fortran2c_types.h"


/* enable runtime_informations */
/*#define pepcparts_SL_USE_RTI
#define pepcparts_SL_USE_RTI_TIM*/


/* standard (SL) integer data type */
#define pepcparts_sl_int_type_c          long long
#define pepcparts_sl_int_type_mpi        MPI_LONG_LONG
#define pepcparts_sl_int_size_mpi        1
#define pepcparts_sl_int_type_fmt        "lld"


/* index data type */
#define pepcparts_sl_index_type_c        FINT_TYPE_C
#define pepcparts_sl_index_type_mpi      FINT_TYPE_MPI
#define pepcparts_sl_index_size_mpi      1
#define pepcparts_sl_index_type_fmt      FINT_TYPE_FMT

/* use indices */
#define pepcparts_SL_INDEX


/* keys */
#define pepcparts_sl_key_type_c          FINT8_TYPE_C
#define pepcparts_sl_key_type_mpi        FINT8_TYPE_MPI
#define pepcparts_sl_key_size_mpi        1
#define pepcparts_sl_key_type_fmt        FINT8_TYPE_FMT
#define pepcparts_sl_key_integer

/* data0: x */
#define pepcparts_SL_DATA0
#define pepcparts_sl_data0_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data0_size_c        1
#define pepcparts_sl_data0_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data0_size_mpi      1

/* data1: y */
#define pepcparts_SL_DATA1
#define pepcparts_sl_data1_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data1_size_c        1
#define pepcparts_sl_data1_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data1_size_mpi      1

/* data2: z */
#define pepcparts_SL_DATA2
#define pepcparts_sl_data2_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data2_size_c        1
#define pepcparts_sl_data2_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data2_size_mpi      1

/* data3: ux */
#define pepcparts_SL_DATA3
#define pepcparts_sl_data3_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data3_size_c        1
#define pepcparts_sl_data3_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data3_size_mpi      1

/* data4: uy */
#define pepcparts_SL_DATA4
#define pepcparts_sl_data4_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data4_size_c        1
#define pepcparts_sl_data4_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data4_size_mpi      1

/* data5: uz */
#define pepcparts_SL_DATA5
#define pepcparts_sl_data5_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data5_size_c        1
#define pepcparts_sl_data5_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data5_size_mpi      1

/* data6: q */
#define pepcparts_SL_DATA6
#define pepcparts_sl_data6_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data6_size_c        1
#define pepcparts_sl_data6_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data6_size_mpi      1

/* data7: m */
#define pepcparts_SL_DATA7
#define pepcparts_sl_data7_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data7_size_c        1
#define pepcparts_sl_data7_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data7_size_mpi      1

/* data8: work */
#define pepcparts_SL_DATA8
#define pepcparts_sl_data8_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data8_size_c        1
#define pepcparts_sl_data8_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data8_size_mpi      1

/* data9: ex */
#define pepcparts_SL_DATA9
#define pepcparts_sl_data9_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data9_size_c        1
#define pepcparts_sl_data9_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data9_size_mpi      1

/* data10: ey */
#define pepcparts_SL_DATA10
#define pepcparts_sl_data10_type_c       FREAL8_TYPE_C
#define pepcparts_sl_data10_size_c       1
#define pepcparts_sl_data10_type_mpi     FREAL8_TYPE_MPI
#define pepcparts_sl_data10_size_mpi     1

/* data11: ez */
#define pepcparts_SL_DATA11
#define pepcparts_sl_data11_type_c       FREAL8_TYPE_C
#define pepcparts_sl_data11_size_c       1
#define pepcparts_sl_data11_type_mpi     FREAL8_TYPE_MPI
#define pepcparts_sl_data11_size_mpi     1

/* data12: pelabel */
#define pepcparts_SL_DATA12
#define pepcparts_sl_data12_type_c       FINT_TYPE_C
#define pepcparts_sl_data12_size_c       1
#define pepcparts_sl_data12_type_mpi     FINT_TYPE_MPI
#define pepcparts_sl_data12_size_mpi     1


/* weighted elements */
#define pepcparts_sl_elem_weight(e, at)  ((e)->data8[at])

#define pepcparts_sl_data8_weight

#define pepcparts_MSEG_BORDER_UPDATE_REDUCTION

#define pepcparts_MSEG_INFO

/*#define pepcparts_MSS_ROOT*/


/* do reduce+bcast instead of allreduce on jugene */
#ifdef JUGENE
# define GLOBAL_REDUCEBCAST_THRESHOLD  0
#endif




#if defined(MSEG_ROOT) && !defined(pepcparts_MSEG_ROOT)
# define pepcparts_MSEG_ROOT  MSEG_ROOT
#endif

#if defined(MSEG_BORDER_UPDATE_REDUCTION) && !defined(pepcparts_MSEG_BORDER_UPDATE_REDUCTION)
# define pepcparts_MSEG_BORDER_UPDATE_REDUCTION  MSEG_BORDER_UPDATE_REDUCTION
#endif

#if defined(MSEG_DISABLE_BEST_CHOICE) && !defined(pepcparts_MSEG_DISABLE_BEST_CHOICE)
# define pepcparts_MSEG_DISABLE_BEST_CHOICE  MSEG_DISABLE_BEST_CHOICE
#endif

#if defined(MSEG_DISABLE_MINMAX) && !defined(pepcparts_MSEG_DISABLE_MINMAX)
# define pepcparts_MSEG_DISABLE_MINMAX  MSEG_DISABLE_MINMAX
#endif

#if defined(MSEG_ENABLE_OPTIMZED_LOWHIGH) && !defined(pepcparts_MSEG_ENABLE_OPTIMZED_LOWHIGH)
# define pepcparts_MSEG_ENABLE_OPTIMZED_LOWHIGH  MSEG_ENABLE_OPTIMZED_LOWHIGH
#endif

#if defined(MSEG_FORWARD_ONLY) && !defined(pepcparts_MSEG_FORWARD_ONLY)
# define pepcparts_MSEG_FORWARD_ONLY  MSEG_FORWARD_ONLY
#endif

#if defined(MSEG_INFO) && !defined(pepcparts_MSEG_INFO)
# define pepcparts_MSEG_INFO  MSEG_INFO
#endif

#if defined(MSEG_TRACE_IF) && !defined(pepcparts_MSEG_TRACE_IF)
# define pepcparts_MSEG_TRACE_IF  MSEG_TRACE_IF
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


#ifndef pepcparts_SL_INDEX
# undef pepcparts_SL_PACKED_INDEX
#endif


/* if no special datatype for (sl default) integer ... */
#ifndef pepcparts_sl_int_type_c
  /* ... use a default one */
# define pepcparts_sl_int_type_c               long      /* sl_macro */
# undef pepcparts_sl_int_type_mpi
# define pepcparts_sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef pepcparts_sl_int_size_mpi
# define pepcparts_sl_int_size_mpi             1         /* sl_macro */
# undef pepcparts_sl_int_type_fmt
# define pepcparts_sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepcparts_sl_int_type_mpi) || !defined(pepcparts_sl_int_size_mpi)
#   error "pepcparts_sl_int_type_mpi and/or pepcparts_sl_int_size_mpi missing"
#  endif
# endif
# ifndef pepcparts_sl_int_type_fmt
#  error "pepcparts_sl_int_type_fmt macro is missing, using d as default"
#  define pepcparts_sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for (intern) weight ... */
#ifndef pepcparts_sl_weight_type_c
 /* ... use (sl default) integer */
# define pepcparts_sl_weight_type_c             pepcparts_sl_int_type_c    /* sl_macro */
# undef pepcparts_sl_weight_type_mpi
# define pepcparts_sl_weight_type_mpi           pepcparts_sl_int_type_mpi  /* sl_macro */
# undef pepcparts_sl_weight_size_mpi
# define pepcparts_sl_weight_size_mpi           pepcparts_sl_int_size_mpi  /* sl_macro */
# undef pepcparts_sl_weight_type_fmt
# define pepcparts_sl_weight_type_fmt           pepcparts_sl_int_type_fmt  /* sl_macro */
# undef pepcparts_sl_weight_intequiv
# define pepcparts_sl_weight_intequiv                            /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepcparts_sl_weight_type_mpi) || !defined(pepcparts_sl_weight_size_mpi)
#   error "pepcparts_sl_weight_type_mpi and/or pepcparts_sl_weight_size_mpi missing"
#  endif
# endif
# ifndef pepcparts_sl_weight_type_fmt
#  error "pepcparts_sl_weight_type_fmt macro is missing, using f as default"
#  define pepcparts_sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef pepcparts_sl_index_type_c
 /* ... use the primary integer type */
# define pepcparts_sl_index_type_c             pepcparts_sl_int_type_c
# undef pepcparts_sl_index_type_mpi
# define pepcparts_sl_index_type_mpi           pepcparts_sl_int_type_mpi
# undef pepcparts_sl_index_size_mpi
# define pepcparts_sl_index_size_mpi           pepcparts_sl_int_size_mpi
# undef pepcparts_sl_index_type_fmt
# define pepcparts_sl_index_type_fmt           pepcparts_sl_int_type_fmt
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepcparts_sl_index_type_mpi) || !defined(pepcparts_sl_index_size_mpi)
#   error "pepcparts_sl_index_type_mpi and/or pepcparts_sl_index_size_mpi missing"
#  endif
# endif
# ifndef pepcparts_sl_index_type_fmt
#  error "pepcparts_sl_index_type_fmt macro is missing, using d as default"
#  define pepcparts_sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef pepcparts_sl_key_pure_type_c
# define pepcparts_sl_key_pure_type_c          pepcparts_sl_key_type_c  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_type_mpi
# define pepcparts_sl_key_pure_type_mpi        pepcparts_sl_key_type_mpi  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_size_mpi
# define pepcparts_sl_key_pure_size_mpi        pepcparts_sl_key_size_mpi  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_type_fmt
# ifdef pepcparts_sl_key_type_fmt
#  define pepcparts_sl_key_pure_type_fmt       pepcparts_sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef pepcparts_sl_key_purify
 /* key val -> key val */
 #define pepcparts_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define pepcparts_sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_set_pure
 /* key component pointer and key val */
 #define pepcparts_sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef pepcparts_sl_key_pure_cmp_eq
 #define pepcparts_sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_ne
 #define pepcparts_sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_lt
 #define pepcparts_sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_le
 #define pepcparts_sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_gt
 #define pepcparts_sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_ge
 #define pepcparts_sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef pepcparts_sl_key_cmp_eq
 #define pepcparts_sl_key_cmp_eq(k0, k1)       (pepcparts_sl_key_pure_cmp_eq(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_ne
 #define pepcparts_sl_key_cmp_ne(k0, k1)       (pepcparts_sl_key_pure_cmp_ne(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_lt
 #define pepcparts_sl_key_cmp_lt(k0, k1)       (pepcparts_sl_key_pure_cmp_lt(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_le
 #define pepcparts_sl_key_cmp_le(k0, k1)       (pepcparts_sl_key_pure_cmp_le(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_gt
 #define pepcparts_sl_key_cmp_gt(k0, k1)       (pepcparts_sl_key_pure_cmp_gt(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_ge
 #define pepcparts_sl_key_cmp_ge(k0, k1)       (pepcparts_sl_key_pure_cmp_ge(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif


/* default random key */
#ifdef pepcparts_sl_key_integer
# if !defined(pepcparts_sl_key_val_srand) || !defined(pepcparts_sl_key_val_rand) || !defined(pepcparts_sl_key_val_rand_minmax)
#  undef pepcparts_sl_key_val_srand
#  undef pepcparts_sl_key_val_rand
#  undef pepcparts_sl_key_val_rand_minmax
#  define pepcparts_sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define pepcparts_sl_key_val_rand()                     ((pepcparts_sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define pepcparts_sl_key_val_rand_minmax(_min_, _max_)  ((pepcparts_sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef pepcparts_SL_DATA0_IGNORE
# undef pepcparts_SL_DATA0
#endif
#ifdef pepcparts_SL_DATA1_IGNORE
# undef pepcparts_SL_DATA1
#endif
#ifdef pepcparts_SL_DATA2_IGNORE
# undef pepcparts_SL_DATA2
#endif
#ifdef pepcparts_SL_DATA3_IGNORE
# undef pepcparts_SL_DATA3
#endif
#ifdef pepcparts_SL_DATA4_IGNORE
# undef pepcparts_SL_DATA4
#endif
#ifdef pepcparts_SL_DATA5_IGNORE
# undef pepcparts_SL_DATA5
#endif
#ifdef pepcparts_SL_DATA6_IGNORE
# undef pepcparts_SL_DATA6
#endif
#ifdef pepcparts_SL_DATA7_IGNORE
# undef pepcparts_SL_DATA7
#endif
#ifdef pepcparts_SL_DATA8_IGNORE
# undef pepcparts_SL_DATA8
#endif
#ifdef pepcparts_SL_DATA9_IGNORE
# undef pepcparts_SL_DATA9
#endif
#ifdef pepcparts_SL_DATA10_IGNORE
# undef pepcparts_SL_DATA10
#endif
#ifdef pepcparts_SL_DATA11_IGNORE
# undef pepcparts_SL_DATA11
#endif
#ifdef pepcparts_SL_DATA12_IGNORE
# undef pepcparts_SL_DATA12
#endif
#ifdef pepcparts_SL_DATA13_IGNORE
# undef pepcparts_SL_DATA13
#endif
#ifdef pepcparts_SL_DATA14_IGNORE
# undef pepcparts_SL_DATA14
#endif
#ifdef pepcparts_SL_DATA15_IGNORE
# undef pepcparts_SL_DATA15
#endif
#ifdef pepcparts_SL_DATA16_IGNORE
# undef pepcparts_SL_DATA16
#endif
#ifdef pepcparts_SL_DATA17_IGNORE
# undef pepcparts_SL_DATA17
#endif
#ifdef pepcparts_SL_DATA18_IGNORE
# undef pepcparts_SL_DATA18
#endif
#ifdef pepcparts_SL_DATA19_IGNORE
# undef pepcparts_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro pepcparts_sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef pepcparts_sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef pepcparts_sl_data0_weight
# undef pepcparts_sl_data1_weight
# undef pepcparts_sl_data2_weight
# undef pepcparts_sl_data3_weight
# undef pepcparts_sl_data4_weight
# undef pepcparts_sl_data5_weight
# undef pepcparts_sl_data6_weight
# undef pepcparts_sl_data7_weight
# undef pepcparts_sl_data8_weight
# undef pepcparts_sl_data9_weight
# undef pepcparts_sl_data10_weight
# undef pepcparts_sl_data11_weight
# undef pepcparts_sl_data12_weight
# undef pepcparts_sl_data13_weight
# undef pepcparts_sl_data14_weight
# undef pepcparts_sl_data15_weight
# undef pepcparts_sl_data16_weight
# undef pepcparts_sl_data17_weight
# undef pepcparts_sl_data18_weight
# undef pepcparts_sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable pepcparts_sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(pepcparts_sl_data0_weight) && !defined(pepcparts_SL_DATA0)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data1_weight) && !defined(pepcparts_SL_DATA1)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data2_weight) && !defined(pepcparts_SL_DATA2)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data3_weight) && !defined(pepcparts_SL_DATA3)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data4_weight) && !defined(pepcparts_SL_DATA4)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data5_weight) && !defined(pepcparts_SL_DATA5)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data6_weight) && !defined(pepcparts_SL_DATA6)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data7_weight) && !defined(pepcparts_SL_DATA7)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data8_weight) && !defined(pepcparts_SL_DATA8)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data9_weight) && !defined(pepcparts_SL_DATA9)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data10_weight) && !defined(pepcparts_SL_DATA10)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data11_weight) && !defined(pepcparts_SL_DATA11)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data12_weight) && !defined(pepcparts_SL_DATA12)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data13_weight) && !defined(pepcparts_SL_DATA13)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data14_weight) && !defined(pepcparts_SL_DATA14)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data15_weight) && !defined(pepcparts_SL_DATA15)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data16_weight) && !defined(pepcparts_SL_DATA16)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data17_weight) && !defined(pepcparts_SL_DATA17)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data18_weight) && !defined(pepcparts_SL_DATA18)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data19_weight) && !defined(pepcparts_SL_DATA19)
# undef pepcparts_sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro pepcparts_FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef pepcparts_SL_DATA0
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data0_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA1
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data1_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA2
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data2_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA3
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data3_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA4
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data4_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA5
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data5_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA6
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data6_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA7
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data7_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA8
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data8_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA9
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data9_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA10
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data10_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA11
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data11_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA12
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data12_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA13
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data13_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA14
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data14_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA15
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data15_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA16
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data16_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA17
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data17_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA18
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data18_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA19
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data19_flex
#   define pepcparts_FLECKS_GUARD
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






#define pepcparts_SPEC_TLOC

typedef pepcparts_sl_int_type_c pepcparts_spec_int_t;

typedef int pepcparts_spec_proc_t;

#define pepcparts_SPEC_LOC_NONE   -1
#define pepcparts_SPEC_PROC_NONE  MPI_PROC_NULL

typedef void *pepcparts_spec_tloc_data_t;
typedef void *pepcparts_spec_tproc_data_t;

struct pepcparts__elements_t;

typedef struct pepcparts__elements_t *pepcparts_spec_elem_buf_t;

typedef struct pepcparts__elements_t pepcparts_spec_elem_t;

typedef pepcparts_sl_int_type_c pepcparts_spec_elem_index_t;

#define pepcparts_spec_elem_set_n(_e_, _n_)     pepcparts_elem_set_size((_e_), (_n_))
#define pepcparts_spec_elem_get_n(_e_)          pepcparts_elem_get_size((_e_))
#define pepcparts_spec_elem_set_nmax(_e_, _n_)  pepcparts_elem_set_max_size((_e_), (_n_))
#define pepcparts_spec_elem_get_nmax(_e_)       pepcparts_elem_get_max_size((_e_))

#define pepcparts_spec_elem_set_buf(_e_, _b_)   *(_e_) = *(_b_)
#define pepcparts_spec_elem_get_buf(_e_)        (_e_)

#define pepcparts_spec_elem_copy_at(_se_, _sat_, _de_, _dat_) \
  elem_copy_at((_se_), (_sat_), (_de_), (_dat_))

#define pepcparts_spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) \
  elem_xchange_at((_s0_), (_s0at_), (_s1_), (_s1at_), (_t_))






/* tproc count */

/* sp_macro pepcparts_SPEC_DECLARE_TPROC_COUNT_DB */
#define pepcparts_SPEC_DECLARE_TPROC_COUNT_DB \
  struct { pepcparts_spec_elem_index_t i; pepcparts_spec_proc_t p; } spec0cd;

/* sp_macro pepcparts_SPEC_DO_TPROC_COUNT_DB */
#define pepcparts_SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec0cd.i = 0; spec0cd.i < pepcparts_spec_elem_get_n(_b_); ++spec0cd.i) { \
    spec0cd.p = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == pepcparts_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROC_COUNT_DB */
#define pepcparts_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts) \
{ \
  pepcparts_SPEC_DECLARE_TPROC_COUNT_DB \
  pepcparts_SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TPROC_COUNT_IP */
#define pepcparts_SPEC_DECLARE_TPROC_COUNT_IP \
  struct { pepcparts_spec_elem_index_t i, t; pepcparts_spec_proc_t p; } spec0ci;

/* sp_macro pepcparts_SPEC_DO_TPROC_COUNT_IP */
#define pepcparts_SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < pepcparts_spec_elem_get_n(_b_); ++spec0ci.i) { \
    spec0ci.p = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == pepcparts_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) pepcparts_spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  pepcparts_spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROC_COUNT_IP */
#define pepcparts_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts) \
{ \
  pepcparts_SPEC_DECLARE_TPROC_COUNT_IP \
  pepcparts_SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro pepcparts_SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define pepcparts_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { pepcparts_spec_elem_index_t i; pepcparts_spec_proc_t p; } spec1cd;

/* sp_macro pepcparts_SPEC_DO_TPROC_MOD_COUNT_DB */
#define pepcparts_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec1cd.i = 0; spec1cd.i < pepcparts_spec_elem_get_n(_b_); ++spec1cd.i) { \
    spec1cd.p = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == pepcparts_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define pepcparts_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts) \
{ \
  pepcparts_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  pepcparts_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define pepcparts_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { pepcparts_spec_elem_index_t i, t; pepcparts_spec_proc_t p; } spec1ci;

/* sp_macro pepcparts_SPEC_DO_TPROC_MOD_COUNT_IP */
#define pepcparts_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < pepcparts_spec_elem_get_n(_b_); ++spec1ci.i) { \
    spec1ci.p = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == pepcparts_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) pepcparts_spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  pepcparts_spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define pepcparts_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts) \
{ \
  pepcparts_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  pepcparts_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro pepcparts_SPEC_DECLARE_TPROCS_COUNT_DB */
#define pepcparts_SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { pepcparts_spec_elem_index_t i; pepcparts_spec_int_t j, n; } spec2cd;

/* sp_macro pepcparts_SPEC_DO_TPROCS_COUNT_DB */
#define pepcparts_SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec2cd.i = 0; spec2cd.i < pepcparts_spec_elem_get_n(_b_); ++spec2cd.i) { \
    spec2cd.n = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.n; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROCS_COUNT_DB */
#define pepcparts_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts, pepcparts_spec_proc_t *procs) \
{ \
  pepcparts_SPEC_DECLARE_TPROCS_COUNT_DB \
  pepcparts_SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TPROCS_COUNT_IP */
#define pepcparts_SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { pepcparts_spec_elem_index_t i, t; pepcparts_spec_int_t j, n; } spec2ci;

/* sp_macro pepcparts_SPEC_DO_TPROCS_COUNT_IP */
#define pepcparts_SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < pepcparts_spec_elem_get_n(_b_); ++spec2ci.i) { \
    spec2ci.n = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), (_ps_)); \
    if (spec2ci.n <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.n; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) pepcparts_spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  pepcparts_spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROCS_COUNT_IP */
#define pepcparts_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts, pepcparts_spec_proc_t *procs) \
{ \
  pepcparts_SPEC_DECLARE_TPROCS_COUNT_IP \
  pepcparts_SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro pepcparts_SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define pepcparts_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { pepcparts_spec_elem_index_t i; pepcparts_spec_int_t j, n; } spec3cd;

/* sp_macro pepcparts_SPEC_DO_TPROCS_MOD_COUNT_DB */
#define pepcparts_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec3cd.i = 0; spec3cd.i < pepcparts_spec_elem_get_n(_b_); ++spec3cd.i) \
  { \
    spec3cd.n = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.n; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define pepcparts_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts, pepcparts_spec_proc_t *procs) \
{ \
  pepcparts_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  pepcparts_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define pepcparts_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { pepcparts_spec_elem_index_t i, t; pepcparts_spec_int_t j, n; } spec3ci;

/* sp_macro pepcparts_SPEC_DO_TPROCS_MOD_COUNT_IP */
#define pepcparts_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < pepcparts_spec_elem_get_n(_b_); ++spec3ci.i) { \
    spec3ci.n = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), (_ps_), NULL); \
    if (spec3ci.n <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.n; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) pepcparts_spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  pepcparts_spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define pepcparts_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts, pepcparts_spec_proc_t *procs) \
{ \
  pepcparts_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  pepcparts_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tproc rearrange */

/* sp_macro pepcparts_SPEC_DECLARE_TPROC_REARRANGE_DB */
#define pepcparts_SPEC_DECLARE_TPROC_REARRANGE_DB \
  struct { pepcparts_spec_elem_index_t i; pepcparts_spec_proc_t p; } spec0d;

/* sp_macro pepcparts_SPEC_DO_TPROC_REARRANGE_DB */
#define pepcparts_SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  for (spec0d.i = 0; spec0d.i < pepcparts_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tp_)(pepcparts_spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == pepcparts_SPEC_PROC_NONE) continue; \
    pepcparts_spec_elem_copy_at((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p]); \
    ++(_ds_)[spec0d.p]; \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROC_REARRANGE_DB */
#define pepcparts_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tproc_data_t tproc_data, int *displs) \
{ \
  pepcparts_SPEC_DECLARE_TPROC_REARRANGE_DB \
  pepcparts_SPEC_DO_TPROC_REARRANGE_DB(_tp_, tproc_data, s, d, displs); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TPROC_REARRANGE_IP */
#define pepcparts_SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { pepcparts_spec_elem_index_t e, i, j; pepcparts_spec_proc_t p, np; } spec0i;

/* sp_macro pepcparts_SPEC_DO_TPROC_REARRANGE_IP */
#define pepcparts_SPEC_DO_TPROC_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_)  do { \
  for (spec0i.e = 0, spec0i.i = 0; spec0i.i < (_n_); ++spec0i.i) { \
    spec0i.e += (_cs_)[spec0i.i]; \
    spec0i.j = (_ds_)[spec0i.i]; \
    while (spec0i.j < spec0i.e) { \
      spec0i.p = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec0i.j, _tpd_); \
      while (spec0i.p != spec0i.i) { \
        spec0i.np = (_tp_)(pepcparts_spec_elem_get_buf(_b_), (_ds_)[spec0i.p], _tpd_); \
        if (spec0i.np != spec0i.p) pepcparts_spec_elem_exchange_at((_b_), (_ds_)[spec0i.p], (_b_), spec0i.j, (_xb_)); \
        ++(_ds_)[spec0i.p]; \
        spec0i.p = spec0i.np; \
      } \
      ++spec0i.j; \
    } \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROC_REARRANGE_IP */
#define pepcparts_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_rearrange_ip(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tproc_data_t tproc_data, int *displs, int *counts, pepcparts_spec_int_t n) \
{ \
  pepcparts_SPEC_DECLARE_TPROC_REARRANGE_IP \
  pepcparts_SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro pepcparts_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define pepcparts_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  struct { pepcparts_spec_elem_index_t i; pepcparts_spec_proc_t p; } spec1d;

/* sp_macro pepcparts_SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define pepcparts_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < pepcparts_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(pepcparts_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, pepcparts_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == pepcparts_SPEC_PROC_NONE) continue; \
      pepcparts_spec_elem_copy_at((_ib_), 0, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < pepcparts_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(pepcparts_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, NULL); \
      if (spec1d.p == pepcparts_SPEC_PROC_NONE) continue; \
      pepcparts_spec_elem_copy_at((_sb_), spec1d.i, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define pepcparts_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tproc_data_t tproc_data, int *displs, pepcparts_spec_elem_t *mod) \
{ \
  pepcparts_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  pepcparts_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define pepcparts_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { pepcparts_spec_elem_index_t e, i, j; pepcparts_spec_proc_t p, np; } spec1i;

/* sp_macro pepcparts_SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define pepcparts_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
      spec1i.e += (_cs_)[spec1i.i]; \
      spec1i.j = (_ds_)[spec1i.i]; \
      while (spec1i.j < spec1i.e) { \
        spec1i.p = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec1i.j, _tpd_, pepcparts_spec_elem_get_buf(_ib_)); \
        pepcparts_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(pepcparts_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, pepcparts_spec_elem_get_buf(_ib_)); \
          if (spec1i.np != spec1i.p) { \
            pepcparts_spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
            pepcparts_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
          } else pepcparts_spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
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
        spec1i.p = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec1i.j, _tpd_, NULL); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(pepcparts_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, NULL); \
          if (spec1i.np != spec1i.p) pepcparts_spec_elem_exchange_at((_b_), (_ds_)[spec1i.p], (_b_), spec1i.j, (_xb_)); \
          ++(_ds_)[spec1i.p]; \
          spec1i.p = spec1i.np; \
        } \
        ++spec1i.j; \
      } \
    } \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define pepcparts_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tproc_data_t tproc_data, int *displs, int *counts, pepcparts_spec_int_t n, pepcparts_spec_elem_t *mod) \
{ \
  pepcparts_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  pepcparts_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro pepcparts_SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define pepcparts_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  struct { pepcparts_spec_elem_index_t i; pepcparts_spec_int_t j, n; } spec2d;

/* sp_macro pepcparts_SPEC_DO_TPROCS_REARRANGE_DB */
#define pepcparts_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  for (spec2d.i = 0; spec2d.i < pepcparts_spec_elem_get_n(_sb_); ++spec2d.i) { \
    spec2d.n = (_tp_)(pepcparts_spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.n; ++spec2d.j) { \
      pepcparts_spec_elem_copy_at((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]]); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROCS_REARRANGE_DB */
#define pepcparts_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tproc_data_t tproc_data, int *displs, pepcparts_spec_proc_t *procs) \
{ \
  pepcparts_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  pepcparts_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define pepcparts_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { pepcparts_spec_elem_index_t e, j, fe, fc, le, lc; pepcparts_spec_int_t i, n, f, l, o; } spec2i;

/* sp_macro pepcparts_SPEC_DO_TPROCS_REARRANGE_IP */
#define pepcparts_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_)  do { \
  spec2i.f = 0; spec2i.fe = (_cs_)[0]; spec2i.fc = pepcparts_spec_elem_get_n(_b_); \
  while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; } \
  spec2i.l = 0; spec2i.le = (_cs_)[0]; spec2i.lc = pepcparts_spec_elem_get_n(_b_) - 1; \
  while (spec2i.lc >= spec2i.le) { ++spec2i.l; spec2i.le += (_cs_)[spec2i.l]; } \
  for (spec2i.e = 0, spec2i.i = 0; spec2i.i < (_n_); ++spec2i.i) { \
    spec2i.e += (_cs_)[spec2i.i]; \
    spec2i.j = (_ds_)[spec2i.i]; \
    while (spec2i.j < spec2i.e) { \
      spec2i.n = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec2i.j, (_tpd_), (_ps_)); \
      spec2i.o = -1; \
      while (spec2i.n > 0) { \
        --spec2i.n; \
        if ((_ps_)[spec2i.n] == spec2i.i && spec2i.o < 0) spec2i.o = spec2i.n; \
        else if ((_ds_)[(_ps_)[spec2i.n]] < spec2i.fc) { \
          spec2i.l = spec2i.f; spec2i.le = spec2i.fe; spec2i.lc = spec2i.fc; \
          if (spec2i.fc < spec2i.fe) { \
            pepcparts_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_b_), spec2i.fc); \
            ++spec2i.fc; \
          } else pepcparts_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec2i.n]] == spec2i.fc) ++spec2i.fc; \
        if (spec2i.j != (_ds_)[(_ps_)[spec2i.n]]) pepcparts_spec_elem_copy_at((_b_), spec2i.j, (_b_), (_ds_)[(_ps_)[spec2i.n]]); \
        ++(_ds_)[(_ps_)[spec2i.n]]; \
        while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; spec2i.fc = (_ds_)[spec2i.f]; } \
      } \
      if (spec2i.o < 0) { \
        if (spec2i.lc < spec2i.le) {  \
          pepcparts_spec_elem_copy_at((_b_), spec2i.lc, (_b_), spec2i.j); \
          spec2i.f = spec2i.l; spec2i.fe = spec2i.le; spec2i.fc = spec2i.lc; \
          --spec2i.lc; \
          while (spec2i.l > 0 && spec2i.lc < (_ds_)[spec2i.l]) { spec2i.le -= (_cs_)[spec2i.l]; spec2i.lc = spec2i.le - 1; --spec2i.l; } \
        } else pepcparts_spec_elem_copy_at((_xb_), 0, (_b_), spec2i.j); \
      } \
      spec2i.j = (_ds_)[spec2i.i]; \
    } \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROCS_REARRANGE_IP */
#define pepcparts_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_ip(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tproc_data_t tproc_data, int *displs, int *counts, pepcparts_spec_int_t n, pepcparts_spec_proc_t *procs) \
{ \
  pepcparts_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  pepcparts_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro pepcparts_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define pepcparts_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  struct { pepcparts_spec_elem_index_t i; pepcparts_spec_int_t j, n; } spec3d;

/* sp_macro pepcparts_SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define pepcparts_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  if (_ib_) { \
    for (spec3d.i = 0; spec3d.i < pepcparts_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(pepcparts_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), pepcparts_spec_elem_get_buf(_ib_)); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        pepcparts_spec_elem_copy_at((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } else { \
    for (spec3d.i = 0; spec3d.i < pepcparts_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(pepcparts_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), NULL); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        pepcparts_spec_elem_copy_at((_sb_), spec3d.i, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define pepcparts_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tproc_data_t tproc_data, int *displs, pepcparts_spec_proc_t *procs, pepcparts_spec_elem_t *mod) \
{ \
  pepcparts_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  pepcparts_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define pepcparts_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { pepcparts_spec_elem_index_t e, j, fe, fc, le, lc; pepcparts_spec_int_t i, n, f, l, o; } spec3i;

/* sp_macro pepcparts_SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define pepcparts_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  if (_ib_) { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = pepcparts_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = pepcparts_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), pepcparts_spec_elem_get_buf(_ib_)); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              pepcparts_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else pepcparts_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          pepcparts_spec_elem_copy_at((_ib_), spec3i.n, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            pepcparts_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else pepcparts_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } else { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = pepcparts_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = pepcparts_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(pepcparts_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), NULL); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              pepcparts_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else pepcparts_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          if (spec3i.j != (_ds_)[(_ps_)[spec3i.n]]) pepcparts_spec_elem_copy_at((_b_), spec3i.j, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            pepcparts_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else pepcparts_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define pepcparts_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tproc_data_t tproc_data, int *displs, int *counts, pepcparts_spec_int_t n, pepcparts_spec_proc_t *procs, pepcparts_spec_elem_t *mod) \
{ \
  pepcparts_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  pepcparts_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}

/* sp_macro pepcparts_SPEC_DEFINE_TPROC */
#define pepcparts_SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  pepcparts_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro pepcparts_SPEC_DEFINE_TPROC_MOD */
#define pepcparts_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  pepcparts_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro pepcparts_SPEC_DEFINE_TPROCS */
#define pepcparts_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  pepcparts_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro pepcparts_SPEC_DEFINE_TPROCS_MOD */
#define pepcparts_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  pepcparts_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  pepcparts_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro pepcparts_SPEC_EXT_PARAM_TPROC pepcparts_SPEC_EXT_PARAM_TPROC_NULL pepcparts_SPEC_EXT_PARAM_TPROC_MOD pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL pepcparts_SPEC_EXT_PARAM_TPROCS pepcparts_SPEC_EXT_PARAM_TPROCS_NULL pepcparts_SPEC_EXT_PARAM_TPROCS_MOD pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define pepcparts_SPEC_EXT_PARAM_TPROC(_name_)       _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip
#define pepcparts_SPEC_EXT_PARAM_TPROC_NULL          NULL, NULL, NULL, NULL
#define pepcparts_SPEC_EXT_PARAM_TPROC_MOD(_name_)   _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip
#define pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL      NULL, NULL, NULL, NULL
#define pepcparts_SPEC_EXT_PARAM_TPROCS(_name_)      _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip
#define pepcparts_SPEC_EXT_PARAM_TPROCS_NULL         NULL, NULL, NULL, NULL
#define pepcparts_SPEC_EXT_PARAM_TPROCS_MOD(_name_)  _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip
#define pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL     NULL, NULL, NULL, NULL


/* sp_type pepcparts_spec_tproc_f pepcparts_spec_tproc_count_f pepcparts_spec_tproc_rearrange_db_f pepcparts_spec_tproc_rearrange_ip_f */
typedef pepcparts_spec_proc_t pepcparts_spec_tproc_f(pepcparts_spec_elem_buf_t b, pepcparts_spec_elem_index_t x, pepcparts_spec_tproc_data_t tproc_data);
typedef void pepcparts_spec_tproc_count_f(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts);
typedef void pepcparts_spec_tproc_rearrange_db_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tproc_data_t tproc_data, int *displs);
typedef void pepcparts_spec_tproc_rearrange_ip_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tproc_data_t tproc_data, int *displs, int *counts, pepcparts_spec_int_t n);

/* sp_type pepcparts_spec_tproc_mod_f pepcparts_spec_tproc_mod_count_f pepcparts_spec_tproc_mod_rearrange_db_f pepcparts_spec_tproc_mod_rearrange_ip_f */
typedef pepcparts_spec_proc_t pepcparts_spec_tproc_mod_f(pepcparts_spec_elem_buf_t b, pepcparts_spec_elem_index_t x, pepcparts_spec_tproc_data_t tproc_data, pepcparts_spec_elem_buf_t mod);
typedef void pepcparts_spec_tproc_mod_count_f(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts);
typedef void pepcparts_spec_tproc_mod_rearrange_db_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tproc_data_t tproc_data, int *displs, pepcparts_spec_elem_t *mod);
typedef void pepcparts_spec_tproc_mod_rearrange_ip_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tproc_data_t tproc_data, int *displs, int *counts, pepcparts_spec_int_t n, pepcparts_spec_elem_t *mod);

/* sp_type pepcparts_spec_tprocs_f pepcparts_spec_tprocs_count_f pepcparts_spec_tprocs_rearrange_db_f pepcparts_spec_tprocs_rearrange_ip_f */
typedef pepcparts_spec_int_t pepcparts_spec_tprocs_f(pepcparts_spec_elem_buf_t b, pepcparts_spec_elem_index_t x, pepcparts_spec_tproc_data_t tproc_data, pepcparts_spec_proc_t *procs);
typedef void pepcparts_spec_tprocs_count_f(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts, pepcparts_spec_proc_t *procs);
typedef void pepcparts_spec_tprocs_rearrange_db_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tproc_data_t tproc_data, int *displs, pepcparts_spec_proc_t *procs);
typedef void pepcparts_spec_tprocs_rearrange_ip_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tproc_data_t tproc_data, int *displs, int *counts, pepcparts_spec_int_t n, pepcparts_spec_proc_t *procs);

/* sp_type pepcparts_spec_tprocs_mod_f pepcparts_spec_tprocs_mod_count_f pepcparts_spec_tprocs_mod_rearrange_db_f pepcparts_spec_tprocs_mod_rearrange_ip_f */
typedef pepcparts_spec_int_t pepcparts_spec_tprocs_mod_f(pepcparts_spec_elem_buf_t b, pepcparts_spec_elem_index_t x, pepcparts_spec_tproc_data_t tproc_data, pepcparts_spec_proc_t *procs, pepcparts_spec_elem_buf_t mod);
typedef void pepcparts_spec_tprocs_mod_count_f(pepcparts_spec_elem_t *s, pepcparts_spec_tproc_data_t tproc_data, int *counts, pepcparts_spec_proc_t *procs);
typedef void pepcparts_spec_tprocs_mod_rearrange_db_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tproc_data_t tproc_data, int *displs, pepcparts_spec_proc_t *procs, pepcparts_spec_elem_t *mod);
typedef void pepcparts_spec_tprocs_mod_rearrange_ip_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tproc_data_t tproc_data, int *displs, int *counts, pepcparts_spec_int_t n, pepcparts_spec_proc_t *procs, pepcparts_spec_elem_t *mod);

/* sp_type pepcparts_spec_tproc_reset_f */
typedef void pepcparts_spec_tproc_reset_f(pepcparts_spec_tproc_data_t tproc_data);


/* enable tloc features */
#ifdef pepcparts_SPEC_TLOC

/* sp_macro pepcparts_SPEC_TLOC pepcparts_SPEC_LOC_NONE */


/* tloc rearrange */

/* sp_macro pepcparts_SPEC_DECLARE_TLOC_REARRANGE_DB */
#define pepcparts_SPEC_DECLARE_TLOC_REARRANGE_DB \
  struct { pepcparts_spec_int_t i, p; } spec0d;

/* sp_macro pepcparts_SPEC_DO_TLOC_REARRANGE_DB */
#define pepcparts_SPEC_DO_TLOC_REARRANGE_DB(_tl_, _tld_, _sb_, _db_)  do { \
  for (spec0d.i = 0; spec0d.i < pepcparts_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tl_)(pepcparts_spec_elem_get_buf(_sb_), spec0d.i, _tld_); \
    if (spec0d.p == pepcparts_SPEC_LOC_NONE) continue; \
    pepcparts_spec_elem_copy_at((_sb_), spec0d.i, (_db_), spec0d.p); \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TLOC_REARRANGE_DB */
#define pepcparts_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_rearrange_db(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tloc_data_t tloc_data) \
{ \
  pepcparts_SPEC_DECLARE_TLOC_REARRANGE_DB \
  pepcparts_SPEC_DO_TLOC_REARRANGE_DB(_tl_, tloc_data, s, d); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TLOC_REARRANGE_IP */
#define pepcparts_SPEC_DECLARE_TLOC_REARRANGE_IP \
  struct { pepcparts_spec_int_t i, p, np; } spec0i;

/* sp_macro pepcparts_SPEC_DO_TLOC_REARRANGE_IP */
#define pepcparts_SPEC_DO_TLOC_REARRANGE_IP(_tl_, _tld_, _b_, _xb_)  do { \
  for (spec0i.i = 0; spec0i.i < pepcparts_spec_elem_get_n(_b_); ++spec0i.i) { \
    spec0i.p = (_tl_)(pepcparts_spec_elem_get_buf(_b_), spec0i.i, _tld_); \
    if (spec0i.p == pepcparts_SPEC_LOC_NONE) continue; \
    while (spec0i.i != spec0i.p) { \
      spec0i.np = (_tl_)(pepcparts_spec_elem_get_buf(_b_), spec0i.p, _tld_); \
      if (spec0i.np == pepcparts_SPEC_LOC_NONE) { pepcparts_spec_elem_copy_at((_b_), spec0i.i, (_b_), spec0i.p); break; } \
      pepcparts_spec_elem_exchange_at((_b_), spec0i.i, (_b_), spec0i.p, (_xb_)); \
      spec0i.p = spec0i.np; \
    } \
  } } while (0)

/* sp_macro pepcparts_SPEC_FUNC_TLOC_REARRANGE_IP */
#define pepcparts_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_rearrange_ip(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tloc_data_t tloc_data) \
{ \
  pepcparts_SPEC_DECLARE_TLOC_REARRANGE_IP \
  pepcparts_SPEC_DO_TLOC_REARRANGE_IP(_tl_, tloc_data, s, x); \
}


/* tloc_mod_mod rearrange */

/* sp_macro pepcparts_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB */
#define pepcparts_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  struct { pepcparts_spec_int_t i, p; } spec1d;

/* sp_macro pepcparts_SPEC_DO_TLOC_MOD_REARRANGE_DB */
#define pepcparts_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, _tld_, _sb_, _db_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < pepcparts_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(pepcparts_spec_elem_get_buf(_sb_), spec1d.i, _tld_, pepcparts_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == pepcparts_SPEC_LOC_NONE) continue; \
      pepcparts_spec_elem_copy_at((_ib_), 0, (_db_), spec1d.p); \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < pepcparts_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(pepcparts_spec_elem_get_buf(_sb_), spec1d.i, _tld_, NULL); \
      if (spec1d.p == pepcparts_SPEC_LOC_NONE) continue; \
      pepcparts_spec_elem_copy_at((_sb_), spec1d.i, (_db_), spec1d.p); \
    } \
  } } while (0) 

/* sp_macro pepcparts_SPEC_FUNC_TLOC_MOD_REARRANGE_DB */
#define pepcparts_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_mod_rearrange_db(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tloc_data_t tloc_data, pepcparts_spec_elem_t *mod) \
{ \
  pepcparts_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  pepcparts_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, tloc_data, s, d, mod); \
}

/* sp_macro pepcparts_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP */
#define pepcparts_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  struct { pepcparts_spec_int_t i, p, np; } spec1i;

/* sp_macro pepcparts_SPEC_DO_TLOC_MOD_REARRANGE_IP */
#define pepcparts_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, _tld_, _b_, _xb_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.i = 0; spec1i.i < pepcparts_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(pepcparts_spec_elem_get_buf(_b_), spec1i.i, _tld_, pepcparts_spec_elem_get_buf(_ib_)); \
      if (spec1i.p == pepcparts_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(pepcparts_spec_elem_get_buf(_b_), spec1i.p, _tld_, pepcparts_spec_elem_get_buf(_xb_)); \
        if (spec1i.np == pepcparts_SPEC_LOC_NONE) break; \
        pepcparts_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.p); \
        pepcparts_spec_elem_copy_at((_xb_), 0, (_ib_), 0); \
        spec1i.p = spec1i.np; \
      } \
      pepcparts_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.i); \
    } \
  } else { \
    for (spec1i.i = 0; spec1i.i < pepcparts_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(pepcparts_spec_elem_get_buf(_b_), spec1i.i, _tld_, NULL); \
      if (spec1i.p == pepcparts_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(pepcparts_spec_elem_get_buf(_b_), spec1i.p, _tld_, NULL); \
        if (spec1i.np == pepcparts_SPEC_LOC_NONE) { pepcparts_spec_elem_copy_at((_b_), spec1i.i, (_b_), spec1i.p); break; } \
        pepcparts_spec_elem_exchange_at((_b_), spec1i.i, (_b_), spec1i.p, (_xb_)); \
        spec1i.p = spec1i.np; \
      } \
    } \
 } } while (0) 

/* sp_macro pepcparts_SPEC_FUNC_TLOC_MOD_REARRANGE_IP */
#define pepcparts_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_mod_rearrange_ip(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tloc_data_t tloc_data, pepcparts_spec_elem_t *mod) \
{ \
  pepcparts_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  pepcparts_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, tloc_data, s, x, mod); \
}

/* sp_macro pepcparts_SPEC_DEFINE_TLOC */
#define pepcparts_SPEC_DEFINE_TLOC(_name_, _tl_, _s_...) \
  pepcparts_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_) \
  pepcparts_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro pepcparts_SPEC_DEFINE_TLOC_MOD */
#define pepcparts_SPEC_DEFINE_TLOC_MOD(_name_, _tl_, _s_...) \
  pepcparts_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_) \
  pepcparts_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro pepcparts_SPEC_EXT_PARAM_TLOC pepcparts_SPEC_EXT_PARAM_TLOC_NULL pepcparts_SPEC_EXT_PARAM_TLOC_MOD pepcparts_SPEC_EXT_PARAM_TLOC_MOD_NULL */
#define pepcparts_SPEC_EXT_PARAM_TLOC(_name_)      _name_##_tloc_rearrange_db, _name_##_tloc_rearrange_ip
#define pepcparts_SPEC_EXT_PARAM_TLOC_NULL         NULL, NULL
#define pepcparts_SPEC_EXT_PARAM_TLOC_MOD(_name_)  _name_##_tloc_mod_rearrange_db, _name_##_tloc_mod_rearrange_ip
#define pepcparts_SPEC_EXT_PARAM_TLOC_MOD_NULL     NULL, NULL


/* sp_type pepcparts_spec_tloc_f pepcparts_spec_tloc_rearrange_db_f pepcparts_spec_tloc_rearrange_ip_f */
typedef pepcparts_spec_elem_index_t pepcparts_spec_tloc_f(pepcparts_spec_elem_buf_t b, pepcparts_spec_elem_index_t x, pepcparts_spec_tloc_data_t tloc_data);
typedef void pepcparts_spec_tloc_rearrange_db_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tloc_data_t tloc_data);
typedef void pepcparts_spec_tloc_rearrange_ip_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tloc_data_t tloc_data);

/* sp_type pepcparts_spec_tloc_mod_f pepcparts_spec_tloc_mod_rearrange_db_f pepcparts_spec_tloc_mod_rearrange_ip_f */
typedef pepcparts_spec_elem_index_t pepcparts_spec_tloc_mod_f(pepcparts_spec_elem_buf_t b, pepcparts_spec_elem_index_t x, pepcparts_spec_tloc_data_t tloc_data, pepcparts_spec_elem_buf_t mod);
typedef void pepcparts_spec_tloc_mod_rearrange_db_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *d, pepcparts_spec_tloc_data_t tloc_data, pepcparts_spec_elem_t *mod);
typedef void pepcparts_spec_tloc_mod_rearrange_ip_f(pepcparts_spec_elem_t *s, pepcparts_spec_elem_t *x, pepcparts_spec_tloc_data_t tloc_data, pepcparts_spec_elem_t *mod);


#endif /* pepcparts_SPEC_TLOC */






#ifdef SL_USE_MPI
# include <mpi.h>
#endif


/* sl_type pepcparts_slint_t pepcparts_slint */
typedef pepcparts_sl_int_type_c pepcparts_slint_t, pepcparts_slint;  /* deprecated 'pepcparts_slint' */

#define pepcparts_slint_fmt   pepcparts_sl_int_type_fmt    /* sl_macro */

/* sl_type pepcparts_slindex_t */
typedef pepcparts_sl_index_type_c pepcparts_slindex_t;

#define pepcparts_sindex_fmt  pepcparts_sl_index_type_fmt  /* sl_macro */

/* sl_type pepcparts_slkey_t */
typedef pepcparts_sl_key_type_c pepcparts_slkey_t;

/* sl_type pepcparts_slkey_pure_t pepcparts_slpkey_t */
typedef pepcparts_sl_key_pure_type_c pepcparts_slkey_pure_t, pepcparts_slpkey_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type pepcparts_sldata0_t */
#ifdef pepcparts_sl_data0_type_c
typedef pepcparts_sl_data0_type_c pepcparts_sldata0_t;
#endif
/* sl_type pepcparts_sldata1_t */
#ifdef pepcparts_sl_data1_type_c
typedef pepcparts_sl_data1_type_c pepcparts_sldata1_t;
#endif
/* sl_type pepcparts_sldata2_t */
#ifdef pepcparts_sl_data2_type_c
typedef pepcparts_sl_data2_type_c pepcparts_sldata2_t;
#endif
/* sl_type pepcparts_sldata3_t */
#ifdef pepcparts_sl_data3_type_c
typedef pepcparts_sl_data3_type_c pepcparts_sldata3_t;
#endif
/* sl_type pepcparts_sldata4_t */
#ifdef pepcparts_sl_data4_type_c
typedef pepcparts_sl_data4_type_c pepcparts_sldata4_t;
#endif
/* sl_type pepcparts_sldata5_t */
#ifdef pepcparts_sl_data5_type_c
typedef pepcparts_sl_data5_type_c pepcparts_sldata5_t;
#endif
/* sl_type pepcparts_sldata6_t */
#ifdef pepcparts_sl_data6_type_c
typedef pepcparts_sl_data6_type_c pepcparts_sldata6_t;
#endif
/* sl_type pepcparts_sldata7_t */
#ifdef pepcparts_sl_data7_type_c
typedef pepcparts_sl_data7_type_c pepcparts_sldata7_t;
#endif
/* sl_type pepcparts_sldata8_t */
#ifdef pepcparts_sl_data8_type_c
typedef pepcparts_sl_data8_type_c pepcparts_sldata8_t;
#endif
/* sl_type pepcparts_sldata9_t */
#ifdef pepcparts_sl_data9_type_c
typedef pepcparts_sl_data9_type_c pepcparts_sldata9_t;
#endif
/* sl_type pepcparts_sldata10_t */
#ifdef pepcparts_sl_data10_type_c
typedef pepcparts_sl_data10_type_c pepcparts_sldata10_t;
#endif
/* sl_type pepcparts_sldata11_t */
#ifdef pepcparts_sl_data11_type_c
typedef pepcparts_sl_data11_type_c pepcparts_sldata11_t;
#endif
/* sl_type pepcparts_sldata12_t */
#ifdef pepcparts_sl_data12_type_c
typedef pepcparts_sl_data12_type_c pepcparts_sldata12_t;
#endif
/* sl_type pepcparts_sldata13_t */
#ifdef pepcparts_sl_data13_type_c
typedef pepcparts_sl_data13_type_c pepcparts_sldata13_t;
#endif
/* sl_type pepcparts_sldata14_t */
#ifdef pepcparts_sl_data14_type_c
typedef pepcparts_sl_data14_type_c pepcparts_sldata14_t;
#endif
/* sl_type pepcparts_sldata15_t */
#ifdef pepcparts_sl_data15_type_c
typedef pepcparts_sl_data15_type_c pepcparts_sldata15_t;
#endif
/* sl_type pepcparts_sldata16_t */
#ifdef pepcparts_sl_data16_type_c
typedef pepcparts_sl_data16_type_c pepcparts_sldata16_t;
#endif
/* sl_type pepcparts_sldata17_t */
#ifdef pepcparts_sl_data17_type_c
typedef pepcparts_sl_data17_type_c pepcparts_sldata17_t;
#endif
/* sl_type pepcparts_sldata18_t */
#ifdef pepcparts_sl_data18_type_c
typedef pepcparts_sl_data18_type_c pepcparts_sldata18_t;
#endif
/* sl_type pepcparts_sldata19_t */
#ifdef pepcparts_sl_data19_type_c
typedef pepcparts_sl_data19_type_c pepcparts_sldata19_t;
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

/* sl_type pepcparts_slweight_t */
typedef pepcparts_sl_weight_type_c pepcparts_slweight_t;

#define pepcparts_slweight_fmt  pepcparts_sl_weight_type_fmt  /* sl_macro */

#if defined(pepcparts_sl_elem_weight) && defined(pepcparts_sl_weight_intequiv)
typedef pepcparts_sl_weight_type_c pepcparts_slcount_t;       /* sl_type pepcparts_slcount_t */
# define pepcparts_slcount_fmt  pepcparts_sl_weight_type_fmt  /* sl_macro */
#else
typedef pepcparts_sl_int_type_c pepcparts_slcount_t;
# define pepcparts_slcount_fmt  pepcparts_sl_int_type_fmt
#endif


/* sl_type pepcparts__slpwkey_t pepcparts_slpwkey_t */
typedef struct pepcparts__slpwkey_t
{
  pepcparts_slpkey_t pkey;
  pepcparts_slweight_t weight;

} pepcparts_slpwkey_t;


/* sl_type pepcparts__elements_t pepcparts_elements_t */
typedef struct pepcparts__elements_t
{
  pepcparts_slint_t size, max_size;
  pepcparts_slkey_t *keys;

#ifdef pepcparts_SL_INDEX
  pepcparts_slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepcparts_SL_DATA0
  pepcparts_sldata0_t *data0;
#endif
#ifdef pepcparts_SL_DATA1
  pepcparts_sldata1_t *data1;
#endif
#ifdef pepcparts_SL_DATA2
  pepcparts_sldata2_t *data2;
#endif
#ifdef pepcparts_SL_DATA3
  pepcparts_sldata3_t *data3;
#endif
#ifdef pepcparts_SL_DATA4
  pepcparts_sldata4_t *data4;
#endif
#ifdef pepcparts_SL_DATA5
  pepcparts_sldata5_t *data5;
#endif
#ifdef pepcparts_SL_DATA6
  pepcparts_sldata6_t *data6;
#endif
#ifdef pepcparts_SL_DATA7
  pepcparts_sldata7_t *data7;
#endif
#ifdef pepcparts_SL_DATA8
  pepcparts_sldata8_t *data8;
#endif
#ifdef pepcparts_SL_DATA9
  pepcparts_sldata9_t *data9;
#endif
#ifdef pepcparts_SL_DATA10
  pepcparts_sldata10_t *data10;
#endif
#ifdef pepcparts_SL_DATA11
  pepcparts_sldata11_t *data11;
#endif
#ifdef pepcparts_SL_DATA12
  pepcparts_sldata12_t *data12;
#endif
#ifdef pepcparts_SL_DATA13
  pepcparts_sldata13_t *data13;
#endif
#ifdef pepcparts_SL_DATA14
  pepcparts_sldata14_t *data14;
#endif
#ifdef pepcparts_SL_DATA15
  pepcparts_sldata15_t *data15;
#endif
#ifdef pepcparts_SL_DATA16
  pepcparts_sldata16_t *data16;
#endif
#ifdef pepcparts_SL_DATA17
  pepcparts_sldata17_t *data17;
#endif
#ifdef pepcparts_SL_DATA18
  pepcparts_sldata18_t *data18;
#endif
#ifdef pepcparts_SL_DATA19
  pepcparts_sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} pepcparts_elements_t;


/* sl_type pepcparts__packed_element_t pepcparts_packed_element_t */
typedef struct pepcparts__packed_element_t
{
  pepcparts_slkey_t key;

#ifdef pepcparts_SL_PACKED_INDEX
  pepcparts_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepcparts_SL_DATA0
# ifdef pepcparts_sl_data0_flex
  pepcparts_sldata0_t data0[];
# else
  pepcparts_sldata0_t data0[pepcparts_sl_data0_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA1
# ifdef pepcparts_sl_data1_flex
  pepcparts_sldata1_t data1[];
# else
  pepcparts_sldata1_t data1[pepcparts_sl_data1_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA2
# ifdef pepcparts_sl_data2_flex
  pepcparts_sldata2_t data2[];
# else
  pepcparts_sldata2_t data2[pepcparts_sl_data2_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA3
# ifdef pepcparts_sl_data3_flex
  pepcparts_sldata3_t data3[];
# else
  pepcparts_sldata3_t data3[pepcparts_sl_data3_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA4
# ifdef pepcparts_sl_data4_flex
  pepcparts_sldata4_t data4[];
# else
  pepcparts_sldata4_t data4[pepcparts_sl_data4_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA5
# ifdef pepcparts_sl_data5_flex
  pepcparts_sldata5_t data5[];
# else
  pepcparts_sldata5_t data5[pepcparts_sl_data5_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA6
# ifdef pepcparts_sl_data6_flex
  pepcparts_sldata6_t data6[];
# else
  pepcparts_sldata6_t data6[pepcparts_sl_data6_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA7
# ifdef pepcparts_sl_data7_flex
  pepcparts_sldata7_t data7[];
# else
  pepcparts_sldata7_t data7[pepcparts_sl_data7_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA8
# ifdef pepcparts_sl_data8_flex
  pepcparts_sldata8_t data8[];
# else
  pepcparts_sldata8_t data8[pepcparts_sl_data8_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA9
# ifdef pepcparts_sl_data9_flex
  pepcparts_sldata9_t data9[];
# else
  pepcparts_sldata9_t data9[pepcparts_sl_data9_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA10
# ifdef pepcparts_sl_data10_flex
  pepcparts_sldata10_t data10[];
# else
  pepcparts_sldata10_t data10[pepcparts_sl_data10_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA11
# ifdef pepcparts_sl_data11_flex
  pepcparts_sldata11_t data11[];
# else
  pepcparts_sldata11_t data11[pepcparts_sl_data11_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA12
# ifdef pepcparts_sl_data12_flex
  pepcparts_sldata12_t data12[];
# else
  pepcparts_sldata12_t data12[pepcparts_sl_data12_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA13
# ifdef pepcparts_sl_data13_flex
  pepcparts_sldata13_t data13[];
# else
  pepcparts_sldata13_t data13[pepcparts_sl_data13_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA14
# ifdef pepcparts_sl_data14_flex
  pepcparts_sldata14_t data14[];
# else
  pepcparts_sldata14_t data14[pepcparts_sl_data14_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA15
# ifdef pepcparts_sl_data15_flex
  pepcparts_sldata15_t data15[];
# else
  pepcparts_sldata15_t data15[pepcparts_sl_data15_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA16
# ifdef pepcparts_sl_data16_flex
  pepcparts_sldata16_t data16[];
# else
  pepcparts_sldata16_t data16[pepcparts_sl_data16_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA17
# ifdef pepcparts_sl_data17_flex
  pepcparts_sldata17_t data17[];
# else
  pepcparts_sldata17_t data17[pepcparts_sl_data17_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA18
# ifdef pepcparts_sl_data18_flex
  pepcparts_sldata18_t data18[];
# else
  pepcparts_sldata18_t data18[pepcparts_sl_data18_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA19
# ifdef pepcparts_sl_data19_flex
  pepcparts_sldata19_t data19[];
# else
  pepcparts_sldata19_t data19[pepcparts_sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} pepcparts_packed_element_t;


/* sl_type pepcparts__packed_elements_t pepcparts_packed_elements_t */
typedef struct pepcparts__packed_elements_t
{
  pepcparts_slint_t size, max_size;
  
  pepcparts_packed_element_t *elements;
  
} pepcparts_packed_elements_t;


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


/* sl_type pepcparts__classification_info_t pepcparts_classification_info_t pepcparts_classification_info */
typedef struct pepcparts__classification_info_t
{
  pepcparts_slint_t nclasses;
  pepcparts_slkey_pure_t *keys;
  pepcparts_slint_t *counts;
  pepcparts_slint_t *masks;

  /* */
  pepcparts_slint_t *all_local_sizes;
  pepcparts_slint_t *local_lt_eq_counts;
  pepcparts_slint_t *all_local_lt_eq_counts;

} pepcparts_classification_info_t, pepcparts_classification_info;  /* deprecated 'pepcparts_classification_info' */


/* key2class, sl_type pepcparts_key2class_f */
typedef pepcparts_slint_t (*pepcparts_key2class_f)(pepcparts_slkey_t *, pepcparts_slint, void *);

/* pivot-element, sl_type pepcparts_pivot_f */
typedef pepcparts_slint_t (*pepcparts_pivot_f)(pepcparts_elements_t *);

/* sorting-network, sl_type pepcparts_sortnet_f pepcparts_sortnet_data_t */
typedef void *pepcparts_sortnet_data_t;
typedef pepcparts_slint_t (*pepcparts_sortnet_f)(pepcparts_slint_t size, pepcparts_slint_t rank, pepcparts_slint_t stage, pepcparts_sortnet_data_t snd, pepcparts_slint_t *up);

/* merge2, sl_type pepcparts_merge2x_f pepcparts_merge2X_f */
typedef pepcparts_slint_t (*pepcparts_merge2x_f)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
typedef pepcparts_slint_t (*pepcparts_merge2X_f)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);

/* sl_type pepcparts__permute_generic_t pepcparts_permute_generic_t */
typedef struct pepcparts__permute_generic_t
{
  int type;

  pepcparts_spec_tloc_f *tloc;
  pepcparts_spec_tloc_rearrange_db_f *tloc_rearrange_db;
  pepcparts_spec_tloc_rearrange_ip_f *tloc_rearrange_ip;

  pepcparts_spec_tloc_mod_f *tloc_mod;
  pepcparts_spec_tloc_mod_rearrange_db_f *tloc_mod_rearrange_db;
  pepcparts_spec_tloc_mod_rearrange_ip_f *tloc_mod_rearrange_ip;

} pepcparts_permute_generic_t;

/* sl_macro pepcparts_PERMUTE_GENERIC_DEFINE_TLOC pepcparts_PERMUTE_GENERIC_INIT_TLOC pepcparts_PERMUTE_GENERIC_INIT_EXT_TLOC */
#define pepcparts_PERMUTE_GENERIC_DEFINE_TLOC(_tl_, _s_...)      pepcparts_SPEC_DEFINE_TLOC(_tl_, _tl_, _s_)
#define pepcparts_PERMUTE_GENERIC_INIT_TLOC(_tl_)                { 1, _tl_, pepcparts_SPEC_EXT_PARAM_TLOC_NULL,  NULL, pepcparts_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define pepcparts_PERMUTE_GENERIC_INIT_EXT_TLOC(_tl_)            { 1, _tl_, pepcparts_SPEC_EXT_PARAM_TLOC(_tl_), NULL, pepcparts_SPEC_EXT_PARAM_TLOC_MOD_NULL }

/* sl_macro pepcparts_PERMUTE_GENERIC_DEFINE_TLOC_MOD pepcparts_PERMUTE_GENERIC_INIT_TLOC_MOD pepcparts_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD */
#define pepcparts_PERMUTE_GENERIC_DEFINE_TLOC_MOD(_tl_, _s_...)  pepcparts_SPEC_DEFINE_TLOC_MOD(_tl_, _tl_, _s_)
#define pepcparts_PERMUTE_GENERIC_INIT_TLOC_MOD(_tl_)            { 2, NULL, pepcparts_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, pepcparts_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define pepcparts_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD(_tl_)        { 2, NULL, pepcparts_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, pepcparts_SPEC_EXT_PARAM_TLOC_MOD(_tl_) }

/* sl_type pepcparts__split_generic_t pepcparts_split_generic_t */
typedef struct pepcparts__split_generic_t
{
  int type;

  pepcparts_spec_tproc_f *tproc;
  pepcparts_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  pepcparts_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  pepcparts_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  pepcparts_spec_tproc_mod_f *tproc_mod;
  pepcparts_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  pepcparts_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  pepcparts_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  pepcparts_spec_tprocs_f *tprocs;
  pepcparts_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  pepcparts_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  pepcparts_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  pepcparts_spec_tprocs_mod_f *tprocs_mod;
  pepcparts_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  pepcparts_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  pepcparts_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

  pepcparts_spec_tproc_reset_f *reset;

} pepcparts_split_generic_t;

/* sl_macro pepcparts_SPLIT_GENERIC_DEFINE_TPROC pepcparts_SPLIT_GENERIC_INIT_TPROC pepcparts_SPLIT_GENERIC_INIT_EXT_TPROC */
#define pepcparts_SPLIT_GENERIC_DEFINE_TPROC(_tp_, _s_...)         pepcparts_SPEC_DEFINE_TPROC(_tp_, _tp_, _s_)
#define pepcparts_SPLIT_GENERIC_INIT_TPROC(_tp_, _r_...)           { 1, _tp_, pepcparts_SPEC_EXT_PARAM_TPROC_NULL,  NULL, pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define pepcparts_SPLIT_GENERIC_INIT_EXT_TPROC(_tp_, _r_...)       { 1, _tp_, pepcparts_SPEC_EXT_PARAM_TPROC(_tp_), NULL, pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro pepcparts_SPLIT_GENERIC_DEFINE_TPROC_MOD pepcparts_SPLIT_GENERIC_INIT_TPROC_MOD pepcparts_SPLIT_GENERIC_INIT_EXT_TPROC_MOD */
#define pepcparts_SPLIT_GENERIC_DEFINE_TPROC_MOD(_tp_, _s_...)     pepcparts_SPEC_DEFINE_TPROC_MOD(_tp_, _tp_, _s_)
#define pepcparts_SPLIT_GENERIC_INIT_TPROC_MOD(_tp_, _r_...)       { 2, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_NULL, _tp_, pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL,  NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define pepcparts_SPLIT_GENERIC_INIT_EXT_TPROC_MOD(_tp_, _r_...)   { 2, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_NULL, _tp_, pepcparts_SPEC_EXT_PARAM_TPROC_MOD(_tp_), NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro pepcparts_SPLIT_GENERIC_DEFINE_TPROCS pepcparts_SPLIT_GENERIC_INIT_TPROCS pepcparts_SPLIT_GENERIC_INIT_EXT_TPROCS */
#define pepcparts_SPLIT_GENERIC_DEFINE_TPROCS(_tp_, _s_...)        pepcparts_SPEC_DEFINE_TPROCS(_tp_, _tp_, _s_)
#define pepcparts_SPLIT_GENERIC_INIT_TPROCS(_tp_, _r_...)          { 3, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, pepcparts_SPEC_EXT_PARAM_TPROCS_NULL,  NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define pepcparts_SPLIT_GENERIC_INIT_EXT_TPROCS(_tp_, _r_...)      { 3, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, pepcparts_SPEC_EXT_PARAM_TPROCS(_tp_), NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro pepcparts_SPLIT_GENERIC_DEFINE_TPROCS_MOD pepcparts_SPLIT_GENERIC_INIT_TPROCS_MOD pepcparts_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD */
#define pepcparts_SPLIT_GENERIC_DEFINE_TPROCS_MOD(_tp_, _s_...)    pepcparts_SPEC_DEFINE_TPROCS_MOD(_tp_, _tp_, _s_)
#define pepcparts_SPLIT_GENERIC_INIT_TPROCS_MOD(_tp_, _r_...)      { 4, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define pepcparts_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD(_tp_, _r_...)  { 4, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD(_tp_), _r_ }

/* sl_type pepcparts_tloc_f pepcparts_tloc_mod_f */
typedef pepcparts_slint_t pepcparts_tloc_f(pepcparts_elements_t *b, pepcparts_slint_t x, void *tloc_data);
typedef pepcparts_slint_t pepcparts_tloc_mod_f(pepcparts_elements_t *b, pepcparts_slint_t x, void *tloc_data, pepcparts_elements_t *mod);

/* sl_type pepcparts_tproc_f pepcparts_tproc_mod_f pepcparts_tprocs_f pepcparts_tprocs_mod_f */
typedef int pepcparts_tproc_f(pepcparts_elements_t *b, pepcparts_slint_t x, void *tproc_data);
typedef int pepcparts_tproc_mod_f(pepcparts_elements_t *b, pepcparts_slint_t x, void *tproc_data, pepcparts_elements_t *mod);
typedef pepcparts_slint_t pepcparts_tprocs_f(pepcparts_elements_t *b, pepcparts_slint_t x, void *tproc_data, int *procs);
typedef pepcparts_slint_t pepcparts_tprocs_mod_f(pepcparts_elements_t *b, pepcparts_slint_t x, void *tproc_data, int *procs, pepcparts_elements_t *mod);

/* sl_type pepcparts_tproc_reset_f */
typedef void pepcparts_tproc_reset_f(void *tproc_data);

/* sl_macro pepcparts_TPROC_RESET_NULL */
#define pepcparts_TPROC_RESET_NULL  NULL

/* sl_type pepcparts__tproc_t pepcparts_tproc_t */
typedef struct pepcparts__tproc_t *pepcparts_tproc_t;

/* sl_type pepcparts__tproc_exdef pepcparts_tproc_exdef */
typedef struct pepcparts__tproc_exdef {
  int type;

  pepcparts_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  pepcparts_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  pepcparts_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  pepcparts_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  pepcparts_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  pepcparts_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  pepcparts_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  pepcparts_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  pepcparts_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  pepcparts_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  pepcparts_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  pepcparts_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

} const *pepcparts_tproc_exdef;

/* sl_macro pepcparts_TPROC_EXDEF_NULL */
#define pepcparts_TPROC_EXDEF_NULL  NULL

/* sl_macro pepcparts_TPROC_EXDEF_DEFINE_TPROC pepcparts_TPROC_EXDEF_DEFINE_TPROC_MOD pepcparts_TPROC_EXDEF_DEFINE_TPROCS pepcparts_TPROC_EXDEF_DEFINE_TPROCS_MOD */
#define pepcparts_TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  pepcparts_SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct pepcparts__tproc_exdef _##_name_ = { 1, pepcparts_SPEC_EXT_PARAM_TPROC(_name_), pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define pepcparts_TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  pepcparts_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct pepcparts__tproc_exdef _##_name_ = { 2, pepcparts_SPEC_EXT_PARAM_TPROC_NULL, pepcparts_SPEC_EXT_PARAM_TPROC_MOD(_name_), pepcparts_SPEC_EXT_PARAM_TPROCS_NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define pepcparts_TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  pepcparts_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct pepcparts__tproc_exdef _##_name_ = { 3, pepcparts_SPEC_EXT_PARAM_TPROC_NULL, pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL, pepcparts_SPEC_EXT_PARAM_TPROCS(_name_), pepcparts_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define pepcparts_TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  pepcparts_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct pepcparts__tproc_exdef _##_name_ = { 4, pepcparts_SPEC_EXT_PARAM_TPROC_NULL, pepcparts_SPEC_EXT_PARAM_TPROC_MOD_NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_NULL, pepcparts_SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;


/* deprecated, sl_type pepcparts_k2c_func pepcparts_pivot_func pepcparts_sn_func pepcparts_m2x_func pepcparts_m2X_func */
typedef pepcparts_key2class_f pepcparts_k2c_func;
typedef pepcparts_pivot_f pepcparts_pivot_func;
typedef pepcparts_sortnet_f pepcparts_sn_func;
typedef pepcparts_merge2x_f pepcparts_m2x_func;
typedef pepcparts_merge2X_f pepcparts_m2X_func;


/* sl_type pepcparts__mergek_t pepcparts_mergek_t */
typedef struct pepcparts__mergek_t
{
  pepcparts_sortnet_f sn;
  pepcparts_sortnet_data_t snd;

  pepcparts_merge2x_f m2x;
  pepcparts_elements_t *sx;

} pepcparts_mergek_t;


/* sl_type pepcparts_keys_init_type_t pepcparts_keys_init_data_t */
typedef pepcparts_slint_t pepcparts_keys_init_type_t;
typedef void *pepcparts_keys_init_data_t;

/* sl_type pepcparts_key_set_data_t pepcparts_key_set_f */
typedef void *pepcparts_key_set_data_t;
typedef void (*pepcparts_key_set_f)(pepcparts_slkey_pure_t *k, pepcparts_key_set_data_t d);


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


/* pepcparts_elements_keys_stats */
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


/* partition conditions, sl_type pepcparts__partcond2_t pepcparts_partcond2_t */
typedef struct pepcparts__partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} pepcparts_partcond2_t;


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

/* partition conditions, sl_type pepcparts__partcond_t pepcparts_partcond_t pepcparts_partcond_p */
typedef struct pepcparts__partcond_t
{
  pepcparts_slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} pepcparts_partcond_t, *pepcparts_partcond_p;


/* internal partition conditions, sl_type pepcparts__partcond_intern_t pepcparts_partcond_intern_t pepcparts_partcond_intern_p */
typedef struct pepcparts__partcond_intern_t
{
  pepcparts_slint_t pcm;
  pepcparts_slint_t count_min, count_max;
  pepcparts_slint_t count_low, count_high;
#ifdef elem_weight
  pepcparts_slweight_t weight_min, weight_max;
  pepcparts_slweight_t weight_low, weight_high;
#endif

} pepcparts_partcond_intern_t, *pepcparts_partcond_intern_p;


/* sl_type pepcparts__parttype_t pepcparts_parttype_t pepcparts_parttype_p */
typedef struct pepcparts__parttype_t
{
  pepcparts_slint_t type;

} pepcparts_parttype_t, *pepcparts_parttype_p;


/* generic binning method */

/* sl_type pepcparts__bin_t pepcparts_bin_t */
typedef struct pepcparts__bin_t
{
  pepcparts_elements_t s;

#ifdef elem_weight
  pepcparts_slweight_t weight;
#endif

} pepcparts_bin_t;


/* sl_type pepcparts__splitter_t pepcparts_splitter_t */
typedef struct pepcparts__splitter_t
{
  pepcparts_slint_t n;

  int *displs;
  pepcparts_slkey_pure_t *s;
  pepcparts_slint_t *sn;

} pepcparts_splitter_t;


struct pepcparts__binning_t;

/* sl_type pepcparts_binning_pre_f pepcparts_binning_exec_f pepcparts_binning_refine_f pepcparts_binning_hit_f pepcparts_binning_finalize_f pepcparts_binning_post_f */
typedef pepcparts_slint_t (*pepcparts_binning_pre_f)(struct pepcparts__binning_t *bm);
typedef pepcparts_slint_t (*pepcparts_binning_exec_f)(struct pepcparts__binning_t *bm, pepcparts_bin_t *bin, pepcparts_slcount_t *counts, pepcparts_slweight_t *weights);
typedef pepcparts_slint_t (*pepcparts_binning_refine_f)(struct pepcparts__binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, pepcparts_slcount_t *counts, pepcparts_slweight_t *weights, pepcparts_splitter_t *sp, pepcparts_slint_t s, pepcparts_bin_t *new_bin);
typedef pepcparts_slint_t (*pepcparts_binning_hit_f)(struct pepcparts__binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, pepcparts_slcount_t *counts, pepcparts_splitter_t *sp, pepcparts_slint_t s);
typedef pepcparts_slint_t (*pepcparts_binning_finalize_f)(struct pepcparts__binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t dc, pepcparts_slweight_t dw, pepcparts_slint_t lc_min, pepcparts_slint_t lc_max, pepcparts_slcount_t *lcs, pepcparts_slweight_t *lws, pepcparts_splitter_t *sp, pepcparts_slint_t s);
typedef pepcparts_slint_t (*pepcparts_binning_post_f)(struct pepcparts__binning_t *bm);


/* sl_type pepcparts__binning_data_t pepcparts_binning_data_t */
typedef union pepcparts__binning_data_t
{
  struct
  {
    pepcparts_slint_t rhigh, rlow, rwidth;
    pepcparts_slint_t rcurrent;
    pepcparts_slkey_pure_t bit_mask;

    pepcparts_elements_t sx;

  } radix;

} pepcparts_binning_data_t;


/* sl_type pepcparts__binning_t pepcparts_binning_t */
typedef struct pepcparts__binning_t
{
  pepcparts_slint_t nbins, max_nbins;
  
  pepcparts_binning_pre_f pre;
  pepcparts_binning_exec_f exec;
  pepcparts_binning_refine_f refine;
  pepcparts_binning_hit_f hit;
  pepcparts_binning_finalize_f finalize;
  pepcparts_binning_post_f post;

  pepcparts_slint_t sorted;

  pepcparts_slint_t docounts;
#ifdef elem_weight
  pepcparts_slint_t doweights;
#endif

  pepcparts_binning_data_t bd;

} pepcparts_binning_t;


/* sl_type pepcparts__local_bins_t pepcparts_local_bins_t */
typedef struct pepcparts__local_bins_t
{
  pepcparts_binning_t *bm;

  pepcparts_slint_t nbins, max_nbins;
  pepcparts_slint_t nelements;

  pepcparts_slint_t docounts;
#ifdef elem_weight
  pepcparts_slint_t doweights;
#endif

  pepcparts_slint_t nbinnings, max_nbinnings;

  pepcparts_slint_t nbins_new, last_new_b, last_new_k;
  pepcparts_bin_t *bins, *bins_new;
  pepcparts_bin_t *bins0, *bins1;

  pepcparts_slint_t *bcws;

#if defined(elem_weight) && defined(pepcparts_sl_weight_intequiv)
  pepcparts_slint_t cw_factor, w_index, bin_cw_factor;
  pepcparts_slweight_t *cws, *bin_cws;
  pepcparts_slweight_t *prefix_cws;
#else
  pepcparts_slint_t *cs, *bin_cs;
  pepcparts_slint_t *prefix_cs;
# ifdef elem_weight
  pepcparts_slweight_t *ws, *bin_ws;
  pepcparts_slweight_t *prefix_ws;
# endif
#endif

  pepcparts_slint_t last_exec_b;

} pepcparts_local_bins_t;


/* sl_type pepcparts__global_bins_t pepcparts_global_bins_t */
typedef struct pepcparts__global_bins_t
{
  pepcparts_binning_t *bm;
  
  pepcparts_local_bins_t lb;

  pepcparts_slint_t *bcws;

#if defined(elem_weight) && defined(pepcparts_sl_weight_intequiv)
  pepcparts_slweight_t *cws;
  pepcparts_slweight_t *prefix_cws;
#else
  pepcparts_slint_t *cs;
  pepcparts_slint_t *prefix_cs;
# ifdef elem_weight
  pepcparts_slweight_t *ws;
  pepcparts_slweight_t *prefix_ws;
# endif
#endif

} pepcparts_global_bins_t;


/* sl_type pepcparts_rti_cmc_t */
typedef struct
{
  pepcparts_slint_t cmp, movek, moved;

} pepcparts_rti_cmc_t;

#ifndef my_rti_ccmp
# define my_rti_ccmp(m)    m.cmc.cmp
# define my_rti_cmovek(m)  m.cmc.movek
# define my_rti_cmoved(m)  m.cmc.moved
#endif


/* sl_type pepcparts_rti_tim_t */
typedef struct
{
  double start, stop;
  double last, cumu;

  pepcparts_slint_t num;

} pepcparts_rti_tim_t[rti_tids];

#ifndef my_rti_tlast
# define my_rti_tlast(m, t)  m.tim[t].last
# define my_rti_tcumu(m, t)  m.tim[t].cumu
# define my_rti_tnum(m, t)   m.tim[t].num
#endif


/* sl_type pepcparts_rti_mem_t */
typedef struct
{
  pepcparts_slint_t nalloc, nfree;
  pepcparts_slint_t max, cur, cur_max;

} pepcparts_rti_mem_t;


/* sl_type pepcparts_rti_t */
typedef struct
{
  /* compare-move-counter */
  pepcparts_rti_cmc_t cmc;
  /* timer */
  pepcparts_rti_tim_t tim;
  /* memory */
  pepcparts_rti_mem_t mem;

} pepcparts_rti_t;

#ifndef my_rti_reset
# define my_rti_reset(m)  memset((void *) &m, 0, sizeof(m))
#endif


/* sl_type pepcparts__sl_context_t pepcparts_sl_context_t */
typedef struct pepcparts__sl_context_t
{

/* src/base/base.c */
  struct {
int dummy_rank;
  } sl;
#ifdef pepcparts_SL_USE_RTI
pepcparts_rti_t rti;
#endif
  struct {
pepcparts_slint_t ip_threshold;
pepcparts_slint_t db_threshold;
pepcparts_slint_t ma_threshold;
  } sr;
  struct {
pepcparts_slint_t threshold;
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
pepcparts_slint_t sendrecv_replace_memsize;
pepcparts_slint_t sendrecv_replace_mpi_maxsize;
  } me;
#endif
#ifdef SL_USE_MPI
  struct {
double t[2];
pepcparts_slint_t max_nprocs;
pepcparts_slint_t packed;
pepcparts_slint_t minalloc;
double overalloc;
  } meas;
#endif
#ifdef SL_USE_MPI
  struct {
pepcparts_slint_t packed;
pepcparts_slint_t db_packed;
pepcparts_slint_t ip_packed;
  } mea;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef pepcparts_MSEG_ROOT
int root;
#endif
#ifdef pepcparts_MSEG_BORDER_UPDATE_REDUCTION
double border_update_count_reduction;
double border_update_weight_reduction;
#endif
#ifdef pepcparts_MSEG_FORWARD_ONLY
pepcparts_slint_t forward_only;
#endif
#ifdef pepcparts_MSEG_INFO
pepcparts_slint_t info_rounds;
pepcparts_slint_t *info_finish_rounds;
double info_finish_rounds_avg;
pepcparts_slweight_t info_total_weights;
#endif
pepcparts_slint_t binnings;
pepcparts_slint_t finalize_mode;
  } mseg;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef pepcparts_MSS_ROOT
int root;
#endif
  } mss;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
pepcparts_slint_t sync;
  } msm;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
pepcparts_slint_t sync;
pepcparts_partcond_t *r_pc;
  } msp;
#endif
#ifdef SL_USE_MPI
  struct {
double i_t[3];
double p_t[3];
double b_t[3];
pepcparts_slint_t sync;
pepcparts_slint_t i_sync;
pepcparts_slint_t p_sync;
pepcparts_slint_t b_sync;
pepcparts_slint_t back_packed;
  } mssp;
#endif
} pepcparts_sl_context_t;






/* sl_macro pepcparts_elem_set_size pepcparts_elem_set_max_size pepcparts_elem_set_keys pepcparts_elem_set_indices */
#define pepcparts_elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pepcparts_elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pepcparts_elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define pepcparts_elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))
/* DATAX_TEMPLATE_BEGIN */
#define pepcparts_elem_set_data0(_e_, _b_)     ((_e_)->data0 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data1(_e_, _b_)     ((_e_)->data1 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data2(_e_, _b_)     ((_e_)->data2 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data3(_e_, _b_)     ((_e_)->data3 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data4(_e_, _b_)     ((_e_)->data4 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data5(_e_, _b_)     ((_e_)->data5 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data6(_e_, _b_)     ((_e_)->data6 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data7(_e_, _b_)     ((_e_)->data7 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data8(_e_, _b_)     ((_e_)->data8 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data9(_e_, _b_)     ((_e_)->data9 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data10(_e_, _b_)     ((_e_)->data10 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data11(_e_, _b_)     ((_e_)->data11 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data12(_e_, _b_)     ((_e_)->data12 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data13(_e_, _b_)     ((_e_)->data13 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data14(_e_, _b_)     ((_e_)->data14 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data15(_e_, _b_)     ((_e_)->data15 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data16(_e_, _b_)     ((_e_)->data16 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data17(_e_, _b_)     ((_e_)->data17 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data18(_e_, _b_)     ((_e_)->data18 = (_b_))  /* sl_macro */
#define pepcparts_elem_set_data19(_e_, _b_)     ((_e_)->data19 = (_b_))  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro pepcparts_elem_get_size pepcparts_elem_get_max_size pepcparts_elem_get_keys pepcparts_elem_get_indices */
#define pepcparts_elem_get_size(_e_)           (_e_)->size
#define pepcparts_elem_get_max_size(_e_)       (_e_)->max_size
#define pepcparts_elem_get_keys(_e_)           (_e_)->keys
#define pepcparts_elem_get_indices(_e_)        (_e_)->indices
/* DATAX_TEMPLATE_BEGIN */
#define pepcparts_elem_get_data0(_e_)          (_e_)->data0  /* sl_macro */
#define pepcparts_elem_get_data1(_e_)          (_e_)->data1  /* sl_macro */
#define pepcparts_elem_get_data2(_e_)          (_e_)->data2  /* sl_macro */
#define pepcparts_elem_get_data3(_e_)          (_e_)->data3  /* sl_macro */
#define pepcparts_elem_get_data4(_e_)          (_e_)->data4  /* sl_macro */
#define pepcparts_elem_get_data5(_e_)          (_e_)->data5  /* sl_macro */
#define pepcparts_elem_get_data6(_e_)          (_e_)->data6  /* sl_macro */
#define pepcparts_elem_get_data7(_e_)          (_e_)->data7  /* sl_macro */
#define pepcparts_elem_get_data8(_e_)          (_e_)->data8  /* sl_macro */
#define pepcparts_elem_get_data9(_e_)          (_e_)->data9  /* sl_macro */
#define pepcparts_elem_get_data10(_e_)          (_e_)->data10  /* sl_macro */
#define pepcparts_elem_get_data11(_e_)          (_e_)->data11  /* sl_macro */
#define pepcparts_elem_get_data12(_e_)          (_e_)->data12  /* sl_macro */
#define pepcparts_elem_get_data13(_e_)          (_e_)->data13  /* sl_macro */
#define pepcparts_elem_get_data14(_e_)          (_e_)->data14  /* sl_macro */
#define pepcparts_elem_get_data15(_e_)          (_e_)->data15  /* sl_macro */
#define pepcparts_elem_get_data16(_e_)          (_e_)->data16  /* sl_macro */
#define pepcparts_elem_get_data17(_e_)          (_e_)->data17  /* sl_macro */
#define pepcparts_elem_get_data18(_e_)          (_e_)->data18  /* sl_macro */
#define pepcparts_elem_get_data19(_e_)          (_e_)->data19  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro pepcparts_elem_set_block pepcparts_elem_set_block_size pepcparts_elem_get_block pepcparts_elem_get_block_size */
#define pepcparts_elem_set_block(_e_, _b_)       ((_e_)->keys = (_b_), (_e_)->max_size = -1)
#define pepcparts_elem_set_block_size(_e_, _s_)  ((_e_)->size = (_s_))
#define pepcparts_elem_get_block(_e_)            ((void *) (((_e_)->max_size < 0)?(_e_)->keys:NULL))
#define pepcparts_elem_get_block_size(_e_)       (_e_)->size

/* sl_macro pepcparts_pelem_set_size pepcparts_pelem_set_max_size pepcparts_pelem_set_elements */
#define pepcparts_pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pepcparts_pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pepcparts_pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))

/* sl_macro pepcparts_pelem_get_size pepcparts_pelem_get_max_size pepcparts_pelem_get_elements */
#define pepcparts_pelem_get_size(_e_)           (_e_)->size
#define pepcparts_pelem_get_max_size(_e_)       (_e_)->max_size
#define pepcparts_pelem_get_elements(_e_)       (_e_)->elements

/* sl_macro pepcparts_SL_DEFCON */
#define pepcparts_SL_DEFCON(_v_)  (pepcparts_sl_default_context._v_)






/* src/base/base.c */
extern pepcparts_sl_context_t pepcparts_sl_default_context;
extern const int pepcparts_default_sl_dummy_rank;
#ifdef pepcparts_SL_USE_RTI
extern const pepcparts_rti_t pepcparts_default_rti;
#endif
extern const pepcparts_slint_t pepcparts_default_sr_ip_threshold;
extern const pepcparts_slint_t pepcparts_default_sr_db_threshold;
extern const pepcparts_slint_t pepcparts_default_sr_ma_threshold;
extern const pepcparts_slint_t pepcparts_default_sri_threshold;

/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
extern const MPI_Datatype pepcparts_default_mpi_int_datatype;
extern const MPI_Datatype pepcparts_default_mpi_key_datatype;
extern const MPI_Datatype pepcparts_default_mpi_pkey_datatype;
extern const MPI_Datatype pepcparts_default_mpi_pwkey_datatype;
extern const MPI_Datatype pepcparts_default_mpi_index_datatype;
extern const MPI_Datatype pepcparts_default_mpi_weight_datatype;
extern const MPI_Datatype pepcparts_default_mpi_data_datatype[];
extern const int pepcparts_default_mpi_rank;
#endif
extern const void *pepcparts_default_me_sendrecv_replace_mem;
extern const pepcparts_slint_t pepcparts_default_me_sendrecv_replace_memsize;
extern const pepcparts_slint_t pepcparts_default_me_sendrecv_replace_mpi_maxsize;
extern const double pepcparts_default_meas_t[];
extern const pepcparts_slint_t pepcparts_default_meas_max_nprocs;
extern const pepcparts_slint_t pepcparts_default_meas_packed;
extern const pepcparts_slint_t pepcparts_default_meas_minalloc;
extern const double pepcparts_default_meas_overalloc;
extern const pepcparts_slint_t pepcparts_default_mea_packed;
extern const pepcparts_slint_t pepcparts_default_mea_db_packed;
extern const pepcparts_slint_t pepcparts_default_mea_ip_packed;
#ifdef pepcparts_MSEG_ROOT
extern const int pepcparts_default_mseg_root;
#endif
#ifdef pepcparts_MSEG_BORDER_UPDATE_REDUCTION
extern const double pepcparts_default_mseg_border_update_count_reduction;
extern const double pepcparts_default_mseg_border_update_weight_reduction;
#endif
#ifdef pepcparts_MSEG_FORWARD_ONLY
extern const pepcparts_slint_t pepcparts_default_mseg_forward_only;
#endif
#ifdef pepcparts_MSEG_INFO
extern const pepcparts_slint_t pepcparts_default_mseg_info_rounds;
extern const pepcparts_slint_t *pepcparts_default_mseg_info_finish_rounds;
extern const double pepcparts_default_mseg_info_finish_rounds_avg;
extern const pepcparts_slweight_t pepcparts_default_mseg_info_total_weights;
#endif
extern const pepcparts_slint_t pepcparts_default_mseg_binnings;
extern const pepcparts_slint_t pepcparts_default_mseg_finalize_mode;
#ifdef pepcparts_MSS_ROOT
extern const int pepcparts_default_mss_root;
#endif
extern const double pepcparts_default_msm_t[];
extern const pepcparts_slint_t pepcparts_default_msm_sync;
extern const double pepcparts_default_msp_t[];
extern const pepcparts_slint_t pepcparts_default_msp_sync;
extern const pepcparts_partcond_t *pepcparts_default_msp_r_pc;
extern const double pepcparts_default_mssp_i_t[];
extern const double pepcparts_default_mssp_p_t[];
extern const double pepcparts_default_mssp_b_t[];
extern const pepcparts_slint_t pepcparts_default_mssp_sync;
extern const pepcparts_slint_t pepcparts_default_mssp_i_sync;
extern const pepcparts_slint_t pepcparts_default_mssp_p_sync;
extern const pepcparts_slint_t pepcparts_default_mssp_b_sync;
extern const pepcparts_slint_t pepcparts_default_mssp_back_packed;






/* src/base/base.c */
pepcparts_slint_t SL_PROTO(pepcparts_binning_create)(pepcparts_local_bins_t *lb, pepcparts_slint_t max_nbins, pepcparts_slint_t max_nbinnings, pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t docounts, pepcparts_slint_t doweights, pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_destroy)(pepcparts_local_bins_t *lb);
pepcparts_slint_t SL_PROTO(pepcparts_binning_pre)(pepcparts_local_bins_t *lb);
pepcparts_slint_t SL_PROTO(pepcparts_binning_exec_reset)(pepcparts_local_bins_t *lb, pepcparts_slint_t do_bins, pepcparts_slint_t do_prefixes);
pepcparts_slint_t SL_PROTO(pepcparts_binning_exec)(pepcparts_local_bins_t *lb, pepcparts_slint_t b, pepcparts_slint_t do_bins, pepcparts_slint_t do_prefixes);
pepcparts_slint_t SL_PROTO(pepcparts_binning_refine)(pepcparts_local_bins_t *lb, pepcparts_slint_t b, pepcparts_slint_t k, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_hit)(pepcparts_local_bins_t *lb, pepcparts_slint_t b, pepcparts_slint_t k, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_finalize)(pepcparts_local_bins_t *lb, pepcparts_slint_t b, pepcparts_slint_t dc, pepcparts_slweight_t dw, pepcparts_slint_t lc_min, pepcparts_slint_t lc_max, pepcparts_slcount_t *lcs, pepcparts_slweight_t *lws, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_post)(pepcparts_local_bins_t *lb);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_create)(pepcparts_binning_t *bm, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, pepcparts_slint_t sorted);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_destroy)(pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_pre)(pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_exec)(pepcparts_binning_t *bm, pepcparts_bin_t *bin, pepcparts_slcount_t *counts, pepcparts_slweight_t *weights);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_refine)(pepcparts_binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, pepcparts_slcount_t *counts, pepcparts_slweight_t *weights, pepcparts_splitter_t *sp, pepcparts_slint_t s, pepcparts_bin_t *new_bin);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_hit)(pepcparts_binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, pepcparts_slcount_t *counts, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_finalize)(pepcparts_binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t dc, pepcparts_slweight_t dw, pepcparts_slint_t lc_min, pepcparts_slint_t lc_max, pepcparts_slcount_t *lcs, pepcparts_slweight_t *lws, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_post)(pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components);
pepcparts_slint_t SL_PROTO(pepcparts_elements_free)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_realloc)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloca)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components);
pepcparts_slint_t SL_PROTO(pepcparts_elements_freea)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_from_blocks)(pepcparts_elements_t *s, pepcparts_slint_t nblocks, void **blocks, pepcparts_slint_t *blocksizes, pepcparts_slint_t alignment, pepcparts_slint_t nmax, slcint_t components);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_from_block)(pepcparts_elements_t *s, void *block, pepcparts_slint_t blocksize, pepcparts_slint_t alignment, pepcparts_slint_t nmax, slcint_t components);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_block)(pepcparts_elements_t *s, void **block, pepcparts_slint_t *blocksize, pepcparts_slint_t alignment, pepcparts_slint_t maxblocksize);
pepcparts_slint_t SL_PROTO(pepcparts_elements_copy)(pepcparts_elements_t *s, pepcparts_elements_t *d);
pepcparts_slint_t SL_PROTO(pepcparts_elements_copy_at)(pepcparts_elements_t *s, pepcparts_slint_t sat, pepcparts_elements_t *d, pepcparts_slint_t dat);
pepcparts_slint_t SL_PROTO(pepcparts_elements_ncopy)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_elements_nmove)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_elements_printf)(pepcparts_elements_t *s, const char *prefix);
pepcparts_slint_t SL_PROTO(pepcparts_elements_extract)(pepcparts_elements_t *src, pepcparts_slint_t nelements, pepcparts_elements_t *dst0, pepcparts_elements_t *dst1);
pepcparts_slint_t SL_PROTO(pepcparts_elements_touch)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_digest_sum)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(pepcparts_elements_crc32)(pepcparts_elements_t *s, pepcparts_slint nelements, pepcparts_slint_t keys, pepcparts_slint_t data);
pepcparts_slint_t SL_PROTO(pepcparts_elements_digest_hash)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components, void *hash);
pepcparts_slint_t SL_PROTO(pepcparts_elements_random_exchange)(pepcparts_elements_t *s, pepcparts_slint_t rounds, pepcparts_elements_t *xs);
pepcparts_slint_t SL_PROTO(pepcparts_elements_keys_init_seed)(unsigned long s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_keys_init)(pepcparts_elements_t *s, pepcparts_keys_init_type_t t, pepcparts_keys_init_data_t d);
pepcparts_slint_t SL_PROTO(pepcparts_elements_keys_init_randomized)(pepcparts_elements_t *s, pepcparts_slint_t nkeys, pepcparts_keys_init_type_t t, pepcparts_keys_init_data_t d);
pepcparts_slint_t SL_PROTO(pepcparts_elements_keys_init_from_file)(pepcparts_elements_t *s, pepcparts_slint_t data, char *filename, pepcparts_slint_t from, pepcparts_slint_t to, pepcparts_slint_t const_bytes_per_line);
pepcparts_slint_t SL_PROTO(pepcparts_elements_keys_save_to_file)(pepcparts_elements_t *s, char *filename);
pepcparts_slint_t SL_PROTO(pepcparts_elements_validate_order)(pepcparts_elements_t *s, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_elements_validate_order_bmask)(pepcparts_elements_t *s, pepcparts_slint_t n, pepcparts_slkey_pure_t bmask);
pepcparts_slint_t SL_PROTO(pepcparts_elements_validate_order_weight)(pepcparts_elements_t *s, pepcparts_slint_t n, pepcparts_slkey_pure_t weight);
pepcparts_slint_t SL_PROTO(pepcparts_elements_keys_stats)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *stats);
pepcparts_slint_t SL_PROTO(pepcparts_elements_keys_stats_print)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_print_keys)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_print_all)(pepcparts_elements_t *s);
pepcparts_slweight_t SL_PROTO(pepcparts_elements_get_weight)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_get_minmax_keys)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slkey_pure_t *minmaxkeys);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_packed)(pepcparts_packed_elements_t *s, pepcparts_slint_t nelements);
pepcparts_slint_t SL_PROTO(pepcparts_elements_free_packed)(pepcparts_packed_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_packed_from_block)(pepcparts_packed_elements_t *s, void *block, pepcparts_slint_t blocksize, pepcparts_slint_t alignment, pepcparts_slint_t nmax);
pepcparts_slint_t SL_PROTO(pepcparts_elements_pack_indexed)(pepcparts_elements_t *s, pepcparts_packed_elements_t *d, pepcparts_slindex_t *rindx, pepcparts_slindex_t *windx);
pepcparts_slint_t SL_PROTO(pepcparts_elements_pack)(pepcparts_elements_t *s, pepcparts_packed_elements_t *d);
pepcparts_slint_t SL_PROTO(pepcparts_elements_pack_at)(pepcparts_elements_t *s, pepcparts_slint_t sat, pepcparts_packed_elements_t *d, pepcparts_slint_t dat);
pepcparts_slint_t SL_PROTO(pepcparts_elements_unpack_indexed)(pepcparts_packed_elements_t *s, pepcparts_elements_t *d, pepcparts_slindex_t *rindx, pepcparts_slindex_t *windx);
pepcparts_slint_t SL_PROTO(pepcparts_elements_unpack)(pepcparts_packed_elements_t *s, pepcparts_elements_t *d);
pepcparts_slint_t SL_PROTO(pepcparts_elements_unpack_at)(pepcparts_packed_elements_t *s, pepcparts_slint_t sat, pepcparts_elements_t *d, pepcparts_slint_t dat);
pepcparts_slint_t SL_PROTO(pepcparts_elements_unpack_keys)(pepcparts_packed_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_auto_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_m2x_func _x0_1, pepcparts_m2x_func _0x_1);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_01_X)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t, pepcparts_m2X_func _X0_1, pepcparts_m2X_func _0X_1);
pepcparts_slint SL_PROTO(pepcparts_merge2_simplify_s1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_slint s1elements);
pepcparts_slint SL_PROTO(pepcparts_merge2_memory_adaptive)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_compo_hula)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_sseq_x0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_sseq_0x_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_sseq_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_sseq_01)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_sbin_x0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_sbin_0x_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_sbin_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_sbin_01)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_shyb_x0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_shyb_0x_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_shyb_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_shyb_01)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_straight_x0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_straight_0x_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_straight_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_straight_x_0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_straight_X0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_straight_0X_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_straight_01_X)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_basic_straight_X0_1u)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_merge2_compo_tridgell)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_2way_ip_int)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t p, int *displs, pepcparts_merge2x_f m2x);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_2way_ip_int_rec)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t p, int *displs, pepcparts_merge2x_f m2x);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_int)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, int *displs, int *counts);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_int_idx)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, int *displs, int *counts);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_idx)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, pepcparts_slindex_t *displs, pepcparts_slindex_t *counts);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_unpack_idx)(pepcparts_packed_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, pepcparts_slindex_t *displs, pepcparts_slindex_t *counts);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_unpack_idxonly)(pepcparts_packed_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, pepcparts_slindex_t *displs, pepcparts_slindex_t *counts);
pepcparts_slint_t SL_PROTO(pepcparts_permute_generic_db)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_permute_generic_t *pg, void *pg_data);
pepcparts_slint_t SL_PROTO(pepcparts_permute_generic_ip)(pepcparts_elements_t *s, pepcparts_elements_t *x, pepcparts_permute_generic_t *pg, void *pg_data);
pepcparts_slint SL_PROTO(pepcparts_sl_search_sequential_lt)(pepcparts_elements_t *s, pepcparts_slpkey_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_sequential_le)(pepcparts_elements_t *s, pepcparts_slpkey_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_sequential_gt)(pepcparts_elements_t *s, pepcparts_slpkey_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_sequential_ge)(pepcparts_elements_t *s, pepcparts_slpkey_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_sequential_lt)(pepcparts_elements_t *s, pepcparts_slpkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_sequential_le)(pepcparts_elements_t *s, pepcparts_slpkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_sequential_gt)(pepcparts_elements_t *s, pepcparts_slpkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_sequential_ge)(pepcparts_elements_t *s, pepcparts_slpkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_lt)(pepcparts_elements_t *s, pepcparts_slpkey_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_le)(pepcparts_elements_t *s, pepcparts_slpkey_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_gt)(pepcparts_elements_t *s, pepcparts_slpkey_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_ge)(pepcparts_elements_t *s, pepcparts_slpkey_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_binary_lt)(pepcparts_elements_t *s, pepcparts_slpkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_binary_le)(pepcparts_elements_t *s, pepcparts_slpkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_binary_gt)(pepcparts_elements_t *s, pepcparts_slpkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_binary_ge)(pepcparts_elements_t *s, pepcparts_slpkey_t *k);
pepcparts_slint_t SL_PROTO(pepcparts_sl_search_binary_lt_bmask)(pepcparts_elements_t *s, pepcparts_slpkey_t k, pepcparts_slpkey_t bmask);
pepcparts_slint_t SL_PROTO(pepcparts_sl_search_binary_le_bmask)(pepcparts_elements_t *s, pepcparts_slpkey_t k, pepcparts_slpkey_t bmask);
pepcparts_slint_t SL_PROTO(pepcparts_sl_search_binary_sign_switch)(pepcparts_elements_t *s);
pepcparts_slint SL_PROTO(pepcparts_sl_search_hybrid_lt)(pepcparts_elements_t *s, pepcparts_slpkey_t k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_hybrid_le)(pepcparts_elements_t *s, pepcparts_slpkey_t k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_hybrid_gt)(pepcparts_elements_t *s, pepcparts_slpkey_t k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_hybrid_ge)(pepcparts_elements_t *s, pepcparts_slpkey_t k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_hybrid_lt)(pepcparts_elements_t *s, pepcparts_slpkey_t *k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_hybrid_le)(pepcparts_elements_t *s, pepcparts_slpkey_t *k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_hybrid_gt)(pepcparts_elements_t *s, pepcparts_slpkey_t *k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_p_hybrid_ge)(pepcparts_elements_t *s, pepcparts_slpkey_t *k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_ilog2c)(pepcparts_slint x);
pepcparts_slint SL_PROTO(pepcparts_ilog2f)(pepcparts_slint x);
pepcparts_slint SL_PROTO(pepcparts_print_bits)(pepcparts_slint v);
pepcparts_slint SL_PROTO(pepcparts_pivot_random)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_counts2displs)(pepcparts_slint_t n, int *counts, int *displs);
pepcparts_slint_t SL_PROTO(pepcparts_displs2counts)(pepcparts_slint_t n, int *displs, int *counts, pepcparts_slint_t total_counts);
void SL_PROTO(pepcparts_get_displcounts_extent)(pepcparts_slint_t n, int *displs, int *counts, pepcparts_slint_t *lb, pepcparts_slint_t *extent);
void SL_PROTO(pepcparts_elem_set_data)(pepcparts_elements_t *e, ...);
pepcparts_slint_t SL_PROTO(pepcparts_elem_get_max_byte)();
pepcparts_slint_t SL_PROTO(pepcparts_elem_reverse)(pepcparts_elements_t *e, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_nxchange_at)(pepcparts_elements_t *e0, pepcparts_slint_t at0, pepcparts_elements_t *e1, pepcparts_slint_t at1, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_nxchange)(pepcparts_elements_t *e0, pepcparts_elements_t *e1, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_nxchange_ro0)(pepcparts_elements_t *e0, pepcparts_elements_t *e1, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_rotate)(pepcparts_elements_t *e, pepcparts_slint_t m, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_rotate_ro0)(pepcparts_elements_t *e, pepcparts_slint_t m, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_rotate_ro1)(pepcparts_elements_t *e, pepcparts_slint_t m, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_sort_counting_use_displs)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ndispls, pepcparts_slint_t *displs);
pepcparts_slint_t SL_PROTO(pepcparts_sort_counting_use_counts)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ncounts, pepcparts_slint_t *counts);
pepcparts_slint_t SL_PROTO(pepcparts_sort_counting_get_counts)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ncounts, pepcparts_slint_t *counts);
pepcparts_slint_t SL_PROTO(pepcparts_sort_counting)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ncounts);
pepcparts_slint SL_PROTO(pepcparts_sort_heap)(pepcparts_elements_t *s, pepcparts_elements_t *xs);
pepcparts_slint_t SL_PROTO(pepcparts_sort_insert_bmask_kernel)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slkey_pure_t bmask);
pepcparts_slint_t SL_PROTO(pepcparts_sort_insert)(pepcparts_elements_t *s, pepcparts_elements_t *sx);
pepcparts_slint_t SL_PROTO(pepcparts_sort_permute_forward)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t *perm, pepcparts_slint_t offset, pepcparts_slint_t mask_bit);
pepcparts_slint_t SL_PROTO(pepcparts_sort_permute_backward)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t *perm, pepcparts_slint_t offset, pepcparts_slint_t mask_bit);
pepcparts_slint SL_PROTO(pepcparts_sort_quick)(pepcparts_elements_t *s, pepcparts_elements_t *xs);
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix_ip)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth);
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix_db)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth);
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix_ma)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth);
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth);
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix_1bit_kernel)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow);
pepcparts_slint SL_PROTO(pepcparts_sort_radix_1bit)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow);
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix_iter)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t presorted, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth);
pepcparts_slint SL_PROTO(pepcparts_sn_hypercube_lh)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_hypercube_hl)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_odd_even_trans)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_odd)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_even)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_batcher)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_bitonic)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_connected)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint_t SL_PROTO(pepcparts_split_generic_db)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_split_generic_t *sg, void *sg_data, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_split_generic_ip)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_split_generic_t *sg, void *sg_data, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_split_generic_count_db)(pepcparts_elements_t *s, pepcparts_split_generic_t *sg, void *sg_data, int *counts, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_split_generic_count_ip)(pepcparts_elements_t *s, pepcparts_split_generic_t *sg, void *sg_data, int *counts, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_split_generic_rearrange_db)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_split_generic_t *sg, void *sg_data, int *counts, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_split_generic_rearrange_ip)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_split_generic_t *sg, void *sg_data, int *counts, int *displs, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_splitter_reset)(pepcparts_splitter_t *sp);
pepcparts_slint_t SL_PROTO(pepcparts_splitx_radix)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t nclasses, pepcparts_slint_t shl, pepcparts_slint_t *counts);
pepcparts_slint SL_PROTO(pepcparts_split2_lt_ge)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t);
pepcparts_slint SL_PROTO(pepcparts_split2_le_gt)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t);
pepcparts_slint SL_PROTO(pepcparts_split3_lt_eq_gt)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t, pepcparts_slint *nlt, pepcparts_slint *nle);
pepcparts_slint SL_PROTO(pepcparts_split3_lt_eq_gt_old)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t, pepcparts_slint *nlt, pepcparts_slint *nle);
pepcparts_slint SL_PROTO(pepcparts_split2_b)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slkey_pure_t bmask);
pepcparts_slint SL_PROTO(pepcparts_splitk_k2c_af)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint k, pepcparts_slint *c, pepcparts_k2c_func k2c, void *k2c_data);
pepcparts_slint SL_PROTO(pepcparts_splitk_k2c)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint k, pepcparts_slint *c, pepcparts_k2c_func k2c, void *k2c_data);
pepcparts_slint SL_PROTO(pepcparts_splitk_k2c_count)(pepcparts_elements_t *s, pepcparts_slint k, pepcparts_slint *c, pepcparts_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI





/* src/base_mpi/base_mpi.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_create)(pepcparts_global_bins_t *gb, pepcparts_slint_t max_nbins, pepcparts_slint_t max_nbinnings, pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t docounts, pepcparts_slint_t doweights, pepcparts_binning_t *bm, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_destroy)(pepcparts_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_pre)(pepcparts_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_exec_reset)(pepcparts_global_bins_t *gb, pepcparts_slint_t do_bins, pepcparts_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_exec_local)(pepcparts_global_bins_t *gb, pepcparts_slint_t b, pepcparts_slint_t do_bins, pepcparts_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_exec_global)(pepcparts_global_bins_t *gb, pepcparts_slint_t do_bins, pepcparts_slint_t do_prefixes, pepcparts_slint_t root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_refine)(pepcparts_global_bins_t *gb, pepcparts_slint_t b, pepcparts_slint_t k, pepcparts_splitter_t *sp, pepcparts_slint_t s, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_hit)(pepcparts_global_bins_t *gb, pepcparts_slint_t b, pepcparts_slint_t k, pepcparts_splitter_t *sp, pepcparts_slint_t s, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_finalize)(pepcparts_global_bins_t *gb, pepcparts_slint_t b, pepcparts_slint_t dc, pepcparts_slweight_t dw, pepcparts_slint_t lc_min, pepcparts_slint_t lc_max, pepcparts_slcount_t *lcs, pepcparts_slweight_t *lws, pepcparts_splitter_t *sp, pepcparts_slint_t s, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_post)(pepcparts_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_datatypes_init)();
pepcparts_slint_t SL_PROTO(pepcparts_mpi_datatypes_release)();
pepcparts_slint_t SL_PROTO(pepcparts_mpi_get_grid_properties)(pepcparts_slint_t ndims, pepcparts_slint_t *dims, pepcparts_slint_t *pos, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_subgroups_create)(pepcparts_slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_subgroups_delete)(pepcparts_slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);
int SL_PROTO(pepcparts_sl_MPI_Allreduce)(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank);
int SL_PROTO(pepcparts_sl_MPI_Alltoall_int)(void *sendbuf, int sendcount, void *recvbuf, int recvcount, MPI_Comm comm, int size, int rank);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_keys_init_from_file)(pepcparts_elements_t *s, char *filename, pepcparts_slint from, pepcparts_slint to, pepcparts_slint const_bytes_per_line, pepcparts_slint root, int size, int rank, MPI_Comm comm);
pepcparts_slint SL_PROTO(pepcparts_mpi_elements_validate_order)(pepcparts_elements_t *s, pepcparts_slint n, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_linear_exchange_pure_keys)(pepcparts_slkey_pure_t *in, pepcparts_slkey_pure_t *out, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_check_order)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t *orders, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_check_global_order)(pepcparts_slkey_pure_t local_min, pepcparts_slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_digest_sum)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(pepcparts_mpi_elements_crc32)(pepcparts_elements_t *s, pepcparts_slint_t n, pepcparts_slint_t keys, pepcparts_slint_t data, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_digest_hash)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_get_counts)(pepcparts_elements_t *s, pepcparts_slint_t *clocal, pepcparts_slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
pepcparts_slweight_t SL_PROTO(pepcparts_mpi_elements_get_weights)(pepcparts_elements_t *s, pepcparts_slweight_t *wlocal, pepcparts_slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_get_counts_and_weights)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t *counts, pepcparts_slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_sendrecv_replace)(pepcparts_elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_tproc_create_tproc)(pepcparts_tproc_t *tproc, pepcparts_tproc_f *tfn, pepcparts_tproc_reset_f *rfn, pepcparts_tproc_exdef exdef);
pepcparts_slint_t SL_PROTO(pepcparts_tproc_create_tproc_mod)(pepcparts_tproc_t *tproc, pepcparts_tproc_mod_f *tfn, pepcparts_tproc_reset_f *rfn, pepcparts_tproc_exdef exdef);
pepcparts_slint_t SL_PROTO(pepcparts_tproc_create_tprocs)(pepcparts_tproc_t *tproc, pepcparts_tprocs_f *tfn, pepcparts_tproc_reset_f *rfn, pepcparts_tproc_exdef exdef);
pepcparts_slint_t SL_PROTO(pepcparts_tproc_create_tprocs_mod)(pepcparts_tproc_t *tproc, pepcparts_tprocs_mod_f *tfn, pepcparts_tproc_reset_f *rfn, pepcparts_tproc_exdef exdef);
pepcparts_slint_t SL_PROTO(pepcparts_tproc_free)(pepcparts_tproc_t *tproc);
pepcparts_slint_t SL_PROTO(pepcparts_tproc_set_proclist)(pepcparts_tproc_t *tproc, pepcparts_slint_t nsend_procs, int *send_procs, pepcparts_slint_t nrecv_procs, int *recv_procs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_tproc_verify)(pepcparts_tproc_t tproc, void *data, pepcparts_elements_t *s, int proc);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoall_specific)(pepcparts_elements_t *sin, pepcparts_elements_t *sout, pepcparts_elements_t *xs, pepcparts_tproc_t tproc, void *data, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoallv_db_packed)(pepcparts_elements_t *sbuf, int *scounts, int *sdispls, pepcparts_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoallv_db)(pepcparts_elements_t *sbuf, int *scounts, int *sdispls, pepcparts_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoallv_ip_packed)(pepcparts_elements_t *s, pepcparts_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoallv_ip_double)(pepcparts_elements_t *s, pepcparts_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoallv_ip_mpi)(pepcparts_elements_t *s, pepcparts_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoallv_ip_dash)(pepcparts_elements_t *s, pepcparts_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoallv_ip)(pepcparts_elements_t *s, pepcparts_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, pepcparts_slint_t structured);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_find_exact_equal)(pepcparts_elements_t *s, pepcparts_slint_t other_rank, pepcparts_slint_t high_rank, pepcparts_slint_t *ex_start, pepcparts_slint_t *ex_size, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_find_exact)(pepcparts_elements_t *s, pepcparts_slint_t other_rank, pepcparts_slint_t high_rank, pepcparts_slint_t *dst_size, pepcparts_slint_t *ex_start, pepcparts_slint_t *ex_sizes, pepcparts_slint_t *nx_move, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_merge2)(pepcparts_elements_t *s, pepcparts_slint_t other_rank, pepcparts_slint_t high_rank, pepcparts_slint_t *dst_size, pepcparts_merge2x_f m2, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_mergek_equal)(pepcparts_elements_t *s, pepcparts_sortnet_f sn, pepcparts_sortnet_data_t snd, pepcparts_merge2x_f m2x, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_mergek_sorted)(pepcparts_elements_t *s, pepcparts_merge2x_f m2x, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_mergek)(pepcparts_elements_t *s, pepcparts_sortnet_f sn, pepcparts_sortnet_data_t snd, pepcparts_merge2x_f m2x, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_mergek_equal2)(pepcparts_elements_t *s, pepcparts_sortnet_f sn, pepcparts_sortnet_data_t snd, pepcparts_merge2x_f m2x, pepcparts_elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_exact_generic)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, pepcparts_binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_exact_radix)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, pepcparts_slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_exact_radix_ngroups)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, pepcparts_slint_t ngroups, MPI_Comm *group_comms, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_exact_radix_2groups)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, MPI_Comm group_comm, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_sample_regular)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_rebalance)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_slint_t stable, pepcparts_slint_t *dst_size, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_rebalance_alltoallv)(pepcparts_elements_t *sbuf, int *scounts, int *sdispls, pepcparts_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
void SL_PROTO(pepcparts_mpi_partcond_set_even)(pepcparts_partcond_t *pcond, pepcparts_slint_t pcm, pepcparts_slint_t ntotal, double nimba, double wtotal, double wimba, int size, int rank);
pepcparts_slint_t SL_PROTO(pepcparts_init_partconds)(pepcparts_slint_t npconds, pepcparts_partcond_t *pconds, pepcparts_slint_t nparts, pepcparts_slint_t total_count, pepcparts_slweight_t total_weight);
pepcparts_slint_t SL_PROTO(pepcparts_init_partconds_intern)(pepcparts_slint_t npconds, pepcparts_partcond_intern_t *pci, pepcparts_partcond_t *pc, pepcparts_slint_t nparts, pepcparts_slint_t total_count, pepcparts_slweight_t total_weight);
pepcparts_slint_t SL_PROTO(pepcparts_merge_partconds)(pepcparts_partcond_t *pconds_in, pepcparts_slint_t npconds_in, pepcparts_partcond_t *pcond_out);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_gather_partconds_grouped)(pepcparts_partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, pepcparts_partcond_t *pconds_out, pepcparts_slint_t *npconds_out, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_gather_partconds)(pepcparts_partcond_t *pcond_in, pepcparts_partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_allgather_partconds)(pepcparts_partcond_t *pcond_in, pepcparts_partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_bcast_partconds)(pepcparts_slint_t npconds, pepcparts_partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_post_check_partconds)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_post_check_partconds_intern)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_stats)(pepcparts_elements_t *s, pepcparts_slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_generic_bulk)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, pepcparts_binning_t *bm, pepcparts_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_generic_grouped)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, pepcparts_binning_t *bm, pepcparts_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_generic)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, pepcparts_binning_t *bm, pepcparts_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_radix)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, pepcparts_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_radix_grouped)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, pepcparts_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_sample_regular)(pepcparts_elements_t *s, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, pepcparts_slint_t nsamples, pepcparts_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_merge)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_merge2)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, pepcparts_slint_t merge_type, pepcparts_slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_merge_radix)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, pepcparts_slint_t merge_type, pepcparts_slint_t sort_type, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, pepcparts_slint_t part_type, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition_radix)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, pepcparts_slint_t part_type, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition_exact_radix)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_partcond_t *pcond, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition_exact_radix_ngroups)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_partcond_t *pcond, pepcparts_slint_t ngroups, MPI_Comm *group_comms, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition_exact_radix_2groups)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_partcond_t *pcond, MPI_Comm group_comm, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_insert_radix)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, pepcparts_slpkey_t *mmkeys, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_presorted_radix)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, pepcparts_slint_t merge_type, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_back)(pepcparts_elements_t *sin, pepcparts_elements_t *sout, pepcparts_elements_t *sx, pepcparts_slpkey_t *lh, pepcparts_slint_t ntotal, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, pepcparts_slint_t ytotal, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_xcounts2ycounts_grouped)(int *xcounts, pepcparts_slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_subxdispls2ycounts)(pepcparts_slint_t nsubs, int *sub_xdispls, pepcparts_slint_t *sub_sources, pepcparts_slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#undef SL_PROTO
#endif /* __SL_PEPCPARTS_H__ */
