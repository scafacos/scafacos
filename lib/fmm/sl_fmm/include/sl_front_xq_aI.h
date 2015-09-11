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


#ifndef __SL_FRONT_XQ_AI_H__
#define __SL_FRONT_XQ_AI_H__

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif /* SL_USE_MPI */

#define SL_PROTO(_f_)  _f_

#include "config_fmm_sort.h"


/* standard (SL) integer data type */
#define fcs_front_xq_aI_sl_int_type_c             SL_INTEGER_C
#define fcs_front_xq_aI_sl_int_type_mpi           SL_INTEGER_MPI
#define fcs_front_xq_aI_sl_int_size_mpi           1
#define fcs_front_xq_aI_sl_int_type_fmt           SL_INTEGER_FMT


/* key section */
#define fcs_front_xq_aI_sl_key_type_c             INTEGER_C
#define fcs_front_xq_aI_sl_key_type_mpi           INTEGER_MPI
#define fcs_front_xq_aI_sl_key_size_mpi           1

#define fcs_front_xq_aI_sl_key_integer
#define fcs_front_xq_aI_sl_key_type_fmt           INTEGER_FMT


/* data section */
#define fcs_front_xq_aI_SL_DATA0                  /* xyz */
#define fcs_front_xq_aI_sl_data0_type_c           REAL_C
#define fcs_front_xq_aI_sl_data0_size_c           3
#define fcs_front_xq_aI_sl_data0_type_mpi         REAL_MPI
#define fcs_front_xq_aI_sl_data0_size_mpi         3

#define fcs_front_xq_aI_SL_DATA1                  /* q */
#define fcs_front_xq_aI_sl_data1_type_c           REAL_C
#define fcs_front_xq_aI_sl_data1_size_c           1
#define fcs_front_xq_aI_sl_data1_type_mpi         REAL_MPI
#define fcs_front_xq_aI_sl_data1_size_mpi         1

#undef fcs_front_xq_aI_SL_DATA2                   /* scr */
#define fcs_front_xq_aI_sl_data2_type_c           INTEGER_C
#define fcs_front_xq_aI_sl_data2_size_c           1
#define fcs_front_xq_aI_sl_data2_type_mpi         INTEGER_MPI
#define fcs_front_xq_aI_sl_data2_size_mpi         1

#define fcs_front_xq_aI_SL_DATA3                  /* addr */
#define fcs_front_xq_aI_sl_data3_type_c           INTEGER_C
#define fcs_front_xq_aI_sl_data3_size_c           1
#define fcs_front_xq_aI_sl_data3_type_mpi         INTEGER_MPI
#define fcs_front_xq_aI_sl_data3_size_mpi         1




#if defined(MSEG_ROOT) && !defined(fcs_front_xq_aI_MSEG_ROOT)
# define fcs_front_xq_aI_MSEG_ROOT  MSEG_ROOT
#endif

#if defined(MSEG_BORDER_UPDATE_REDUCTION) && !defined(fcs_front_xq_aI_MSEG_BORDER_UPDATE_REDUCTION)
# define fcs_front_xq_aI_MSEG_BORDER_UPDATE_REDUCTION  MSEG_BORDER_UPDATE_REDUCTION
#endif

#if defined(MSEG_DISABLE_BEST_CHOICE) && !defined(fcs_front_xq_aI_MSEG_DISABLE_BEST_CHOICE)
# define fcs_front_xq_aI_MSEG_DISABLE_BEST_CHOICE  MSEG_DISABLE_BEST_CHOICE
#endif

#if defined(MSEG_DISABLE_MINMAX) && !defined(fcs_front_xq_aI_MSEG_DISABLE_MINMAX)
# define fcs_front_xq_aI_MSEG_DISABLE_MINMAX  MSEG_DISABLE_MINMAX
#endif

#if defined(MSEG_ENABLE_OPTIMZED_LOWHIGH) && !defined(fcs_front_xq_aI_MSEG_ENABLE_OPTIMZED_LOWHIGH)
# define fcs_front_xq_aI_MSEG_ENABLE_OPTIMZED_LOWHIGH  MSEG_ENABLE_OPTIMZED_LOWHIGH
#endif

#if defined(MSEG_FORWARD_ONLY) && !defined(fcs_front_xq_aI_MSEG_FORWARD_ONLY)
# define fcs_front_xq_aI_MSEG_FORWARD_ONLY  MSEG_FORWARD_ONLY
#endif

#if defined(MSEG_INFO) && !defined(fcs_front_xq_aI_MSEG_INFO)
# define fcs_front_xq_aI_MSEG_INFO  MSEG_INFO
#endif

#if defined(MSEG_TRACE_IF) && !defined(fcs_front_xq_aI_MSEG_TRACE_IF)
# define fcs_front_xq_aI_MSEG_TRACE_IF  MSEG_TRACE_IF
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


#ifndef fcs_front_xq_aI_SL_INDEX
# undef fcs_front_xq_aI_SL_PACKED_INDEX
#endif


/* if no special datatype for (sl default) integer ... */
#ifndef fcs_front_xq_aI_sl_int_type_c
  /* ... use a default one */
# define fcs_front_xq_aI_sl_int_type_c               long      /* sl_macro */
# undef fcs_front_xq_aI_sl_int_type_mpi
# define fcs_front_xq_aI_sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef fcs_front_xq_aI_sl_int_size_mpi
# define fcs_front_xq_aI_sl_int_size_mpi             1         /* sl_macro */
# undef fcs_front_xq_aI_sl_int_type_fmt
# define fcs_front_xq_aI_sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(fcs_front_xq_aI_sl_int_type_mpi) || !defined(fcs_front_xq_aI_sl_int_size_mpi)
#   error "fcs_front_xq_aI_sl_int_type_mpi and/or fcs_front_xq_aI_sl_int_size_mpi missing"
#  endif
# endif
# ifndef fcs_front_xq_aI_sl_int_type_fmt
#  error "fcs_front_xq_aI_sl_int_type_fmt macro is missing, using d as default"
#  define fcs_front_xq_aI_sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for (intern) weight ... */
#ifndef fcs_front_xq_aI_sl_weight_type_c
 /* ... use (sl default) integer */
# define fcs_front_xq_aI_sl_weight_type_c             fcs_front_xq_aI_sl_int_type_c    /* sl_macro */
# undef fcs_front_xq_aI_sl_weight_type_mpi
# define fcs_front_xq_aI_sl_weight_type_mpi           fcs_front_xq_aI_sl_int_type_mpi  /* sl_macro */
# undef fcs_front_xq_aI_sl_weight_size_mpi
# define fcs_front_xq_aI_sl_weight_size_mpi           fcs_front_xq_aI_sl_int_size_mpi  /* sl_macro */
# undef fcs_front_xq_aI_sl_weight_type_fmt
# define fcs_front_xq_aI_sl_weight_type_fmt           fcs_front_xq_aI_sl_int_type_fmt  /* sl_macro */
# undef fcs_front_xq_aI_sl_weight_intequiv
# define fcs_front_xq_aI_sl_weight_intequiv                            /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(fcs_front_xq_aI_sl_weight_type_mpi) || !defined(fcs_front_xq_aI_sl_weight_size_mpi)
#   error "fcs_front_xq_aI_sl_weight_type_mpi and/or fcs_front_xq_aI_sl_weight_size_mpi missing"
#  endif
# endif
# ifndef fcs_front_xq_aI_sl_weight_type_fmt
#  error "fcs_front_xq_aI_sl_weight_type_fmt macro is missing, using f as default"
#  define fcs_front_xq_aI_sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef fcs_front_xq_aI_sl_index_type_c
 /* ... use the primary integer type */
# define fcs_front_xq_aI_sl_index_type_c             fcs_front_xq_aI_sl_int_type_c
# undef fcs_front_xq_aI_sl_index_type_mpi
# define fcs_front_xq_aI_sl_index_type_mpi           fcs_front_xq_aI_sl_int_type_mpi
# undef fcs_front_xq_aI_sl_index_size_mpi
# define fcs_front_xq_aI_sl_index_size_mpi           fcs_front_xq_aI_sl_int_size_mpi
# undef fcs_front_xq_aI_sl_index_type_fmt
# define fcs_front_xq_aI_sl_index_type_fmt           fcs_front_xq_aI_sl_int_type_fmt
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(fcs_front_xq_aI_sl_index_type_mpi) || !defined(fcs_front_xq_aI_sl_index_size_mpi)
#   error "fcs_front_xq_aI_sl_index_type_mpi and/or fcs_front_xq_aI_sl_index_size_mpi missing"
#  endif
# endif
# ifndef fcs_front_xq_aI_sl_index_type_fmt
#  error "fcs_front_xq_aI_sl_index_type_fmt macro is missing, using d as default"
#  define fcs_front_xq_aI_sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef fcs_front_xq_aI_sl_key_pure_type_c
# define fcs_front_xq_aI_sl_key_pure_type_c          fcs_front_xq_aI_sl_key_type_c  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_pure_type_mpi
# define fcs_front_xq_aI_sl_key_pure_type_mpi        fcs_front_xq_aI_sl_key_type_mpi  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_pure_size_mpi
# define fcs_front_xq_aI_sl_key_pure_size_mpi        fcs_front_xq_aI_sl_key_size_mpi  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_pure_type_fmt
# ifdef fcs_front_xq_aI_sl_key_type_fmt
#  define fcs_front_xq_aI_sl_key_pure_type_fmt       fcs_front_xq_aI_sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef fcs_front_xq_aI_sl_key_purify
 /* key val -> key val */
 #define fcs_front_xq_aI_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define fcs_front_xq_aI_sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_set_pure
 /* key component pointer and key val */
 #define fcs_front_xq_aI_sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef fcs_front_xq_aI_sl_key_pure_cmp_eq
 #define fcs_front_xq_aI_sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_pure_cmp_ne
 #define fcs_front_xq_aI_sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_pure_cmp_lt
 #define fcs_front_xq_aI_sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_pure_cmp_le
 #define fcs_front_xq_aI_sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_pure_cmp_gt
 #define fcs_front_xq_aI_sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_pure_cmp_ge
 #define fcs_front_xq_aI_sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef fcs_front_xq_aI_sl_key_cmp_eq
 #define fcs_front_xq_aI_sl_key_cmp_eq(k0, k1)       (fcs_front_xq_aI_sl_key_pure_cmp_eq(fcs_front_xq_aI_sl_key_purify(k0), fcs_front_xq_aI_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_cmp_ne
 #define fcs_front_xq_aI_sl_key_cmp_ne(k0, k1)       (fcs_front_xq_aI_sl_key_pure_cmp_ne(fcs_front_xq_aI_sl_key_purify(k0), fcs_front_xq_aI_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_cmp_lt
 #define fcs_front_xq_aI_sl_key_cmp_lt(k0, k1)       (fcs_front_xq_aI_sl_key_pure_cmp_lt(fcs_front_xq_aI_sl_key_purify(k0), fcs_front_xq_aI_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_cmp_le
 #define fcs_front_xq_aI_sl_key_cmp_le(k0, k1)       (fcs_front_xq_aI_sl_key_pure_cmp_le(fcs_front_xq_aI_sl_key_purify(k0), fcs_front_xq_aI_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_cmp_gt
 #define fcs_front_xq_aI_sl_key_cmp_gt(k0, k1)       (fcs_front_xq_aI_sl_key_pure_cmp_gt(fcs_front_xq_aI_sl_key_purify(k0), fcs_front_xq_aI_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_cmp_ge
 #define fcs_front_xq_aI_sl_key_cmp_ge(k0, k1)       (fcs_front_xq_aI_sl_key_pure_cmp_ge(fcs_front_xq_aI_sl_key_purify(k0), fcs_front_xq_aI_sl_key_purify(k1)))  /* sl_macro */
#endif


/* default random key */
#ifdef fcs_front_xq_aI_sl_key_integer
# if !defined(fcs_front_xq_aI_sl_key_val_srand) || !defined(fcs_front_xq_aI_sl_key_val_rand) || !defined(fcs_front_xq_aI_sl_key_val_rand_minmax)
#  undef fcs_front_xq_aI_sl_key_val_srand
#  undef fcs_front_xq_aI_sl_key_val_rand
#  undef fcs_front_xq_aI_sl_key_val_rand_minmax
#  define fcs_front_xq_aI_sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define fcs_front_xq_aI_sl_key_val_rand()                     ((fcs_front_xq_aI_sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define fcs_front_xq_aI_sl_key_val_rand_minmax(_min_, _max_)  ((fcs_front_xq_aI_sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef fcs_front_xq_aI_SL_DATA0_IGNORE
# undef fcs_front_xq_aI_SL_DATA0
#endif
#ifdef fcs_front_xq_aI_SL_DATA1_IGNORE
# undef fcs_front_xq_aI_SL_DATA1
#endif
#ifdef fcs_front_xq_aI_SL_DATA2_IGNORE
# undef fcs_front_xq_aI_SL_DATA2
#endif
#ifdef fcs_front_xq_aI_SL_DATA3_IGNORE
# undef fcs_front_xq_aI_SL_DATA3
#endif
#ifdef fcs_front_xq_aI_SL_DATA4_IGNORE
# undef fcs_front_xq_aI_SL_DATA4
#endif
#ifdef fcs_front_xq_aI_SL_DATA5_IGNORE
# undef fcs_front_xq_aI_SL_DATA5
#endif
#ifdef fcs_front_xq_aI_SL_DATA6_IGNORE
# undef fcs_front_xq_aI_SL_DATA6
#endif
#ifdef fcs_front_xq_aI_SL_DATA7_IGNORE
# undef fcs_front_xq_aI_SL_DATA7
#endif
#ifdef fcs_front_xq_aI_SL_DATA8_IGNORE
# undef fcs_front_xq_aI_SL_DATA8
#endif
#ifdef fcs_front_xq_aI_SL_DATA9_IGNORE
# undef fcs_front_xq_aI_SL_DATA9
#endif
#ifdef fcs_front_xq_aI_SL_DATA10_IGNORE
# undef fcs_front_xq_aI_SL_DATA10
#endif
#ifdef fcs_front_xq_aI_SL_DATA11_IGNORE
# undef fcs_front_xq_aI_SL_DATA11
#endif
#ifdef fcs_front_xq_aI_SL_DATA12_IGNORE
# undef fcs_front_xq_aI_SL_DATA12
#endif
#ifdef fcs_front_xq_aI_SL_DATA13_IGNORE
# undef fcs_front_xq_aI_SL_DATA13
#endif
#ifdef fcs_front_xq_aI_SL_DATA14_IGNORE
# undef fcs_front_xq_aI_SL_DATA14
#endif
#ifdef fcs_front_xq_aI_SL_DATA15_IGNORE
# undef fcs_front_xq_aI_SL_DATA15
#endif
#ifdef fcs_front_xq_aI_SL_DATA16_IGNORE
# undef fcs_front_xq_aI_SL_DATA16
#endif
#ifdef fcs_front_xq_aI_SL_DATA17_IGNORE
# undef fcs_front_xq_aI_SL_DATA17
#endif
#ifdef fcs_front_xq_aI_SL_DATA18_IGNORE
# undef fcs_front_xq_aI_SL_DATA18
#endif
#ifdef fcs_front_xq_aI_SL_DATA19_IGNORE
# undef fcs_front_xq_aI_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro fcs_front_xq_aI_sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef fcs_front_xq_aI_sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef fcs_front_xq_aI_sl_data0_weight
# undef fcs_front_xq_aI_sl_data1_weight
# undef fcs_front_xq_aI_sl_data2_weight
# undef fcs_front_xq_aI_sl_data3_weight
# undef fcs_front_xq_aI_sl_data4_weight
# undef fcs_front_xq_aI_sl_data5_weight
# undef fcs_front_xq_aI_sl_data6_weight
# undef fcs_front_xq_aI_sl_data7_weight
# undef fcs_front_xq_aI_sl_data8_weight
# undef fcs_front_xq_aI_sl_data9_weight
# undef fcs_front_xq_aI_sl_data10_weight
# undef fcs_front_xq_aI_sl_data11_weight
# undef fcs_front_xq_aI_sl_data12_weight
# undef fcs_front_xq_aI_sl_data13_weight
# undef fcs_front_xq_aI_sl_data14_weight
# undef fcs_front_xq_aI_sl_data15_weight
# undef fcs_front_xq_aI_sl_data16_weight
# undef fcs_front_xq_aI_sl_data17_weight
# undef fcs_front_xq_aI_sl_data18_weight
# undef fcs_front_xq_aI_sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable fcs_front_xq_aI_sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(fcs_front_xq_aI_sl_data0_weight) && !defined(fcs_front_xq_aI_SL_DATA0)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data1_weight) && !defined(fcs_front_xq_aI_SL_DATA1)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data2_weight) && !defined(fcs_front_xq_aI_SL_DATA2)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data3_weight) && !defined(fcs_front_xq_aI_SL_DATA3)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data4_weight) && !defined(fcs_front_xq_aI_SL_DATA4)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data5_weight) && !defined(fcs_front_xq_aI_SL_DATA5)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data6_weight) && !defined(fcs_front_xq_aI_SL_DATA6)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data7_weight) && !defined(fcs_front_xq_aI_SL_DATA7)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data8_weight) && !defined(fcs_front_xq_aI_SL_DATA8)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data9_weight) && !defined(fcs_front_xq_aI_SL_DATA9)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data10_weight) && !defined(fcs_front_xq_aI_SL_DATA10)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data11_weight) && !defined(fcs_front_xq_aI_SL_DATA11)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data12_weight) && !defined(fcs_front_xq_aI_SL_DATA12)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data13_weight) && !defined(fcs_front_xq_aI_SL_DATA13)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data14_weight) && !defined(fcs_front_xq_aI_SL_DATA14)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data15_weight) && !defined(fcs_front_xq_aI_SL_DATA15)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data16_weight) && !defined(fcs_front_xq_aI_SL_DATA16)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data17_weight) && !defined(fcs_front_xq_aI_SL_DATA17)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data18_weight) && !defined(fcs_front_xq_aI_SL_DATA18)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
#if defined(fcs_front_xq_aI_sl_data19_weight) && !defined(fcs_front_xq_aI_SL_DATA19)
# undef fcs_front_xq_aI_sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro fcs_front_xq_aI_FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef fcs_front_xq_aI_SL_DATA0
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data0_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA1
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data1_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA2
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data2_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA3
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data3_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA4
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data4_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA5
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data5_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA6
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data6_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA7
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data7_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA8
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data8_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA9
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data9_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA10
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data10_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA11
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data11_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA12
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data12_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA13
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data13_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA14
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data14_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA15
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data15_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA16
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data16_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA17
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data17_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA18
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data18_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA19
# ifdef fcs_front_xq_aI_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_front_xq_aI_sl_data19_flex
#   define fcs_front_xq_aI_FLECKS_GUARD
#  endif
# endif
#endif
/* DATAX_TEMPLATE_END */



#define fcs_front_xq_aI_MC_ALLTOALL_INT_2STEP_THRESHOLD  1024




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



/* SL_DATA_IGNORE -> prevents that the present data is processed */




/* prevent the following functions from processing the DATAs (even though they may exist) */
#ifdef SL_DATA_IGNORE
 /* disable the single DATAs */
/* DATAX_TEMPLATE_BEGIN */
 #undef fcs_front_xq_aI_SL_DATA0
 #undef fcs_front_xq_aI_SL_DATA1
 #undef fcs_front_xq_aI_SL_DATA2
 #undef fcs_front_xq_aI_SL_DATA3
 #undef fcs_front_xq_aI_SL_DATA4
 #undef fcs_front_xq_aI_SL_DATA5
 #undef fcs_front_xq_aI_SL_DATA6
 #undef fcs_front_xq_aI_SL_DATA7
 #undef fcs_front_xq_aI_SL_DATA8
 #undef fcs_front_xq_aI_SL_DATA9
 #undef fcs_front_xq_aI_SL_DATA10
 #undef fcs_front_xq_aI_SL_DATA11
 #undef fcs_front_xq_aI_SL_DATA12
 #undef fcs_front_xq_aI_SL_DATA13
 #undef fcs_front_xq_aI_SL_DATA14
 #undef fcs_front_xq_aI_SL_DATA15
 #undef fcs_front_xq_aI_SL_DATA16
 #undef fcs_front_xq_aI_SL_DATA17
 #undef fcs_front_xq_aI_SL_DATA18
 #undef fcs_front_xq_aI_SL_DATA19
/* DATAX_TEMPLATE_END */
#endif /* SL_DATA_IGNORE */


#ifndef SL_ARRAYX_COPY
# define SL_ARRAYX_COPY
# define SL_ARRAY1_COPY(_s_, _d_)  ((_d_)[0] = (_s_)[0])
# define SL_ARRAY2_COPY(_s_, _d_)  SL_ARRAY1_COPY(_s_, _d_), ((_d_)[1] = (_s_)[1])
# define SL_ARRAY3_COPY(_s_, _d_)  SL_ARRAY2_COPY(_s_, _d_), ((_d_)[2] = (_s_)[2])
# define SL_ARRAY4_COPY(_s_, _d_)  SL_ARRAY3_COPY(_s_, _d_), ((_d_)[3] = (_s_)[3])
# define SL_ARRAY5_COPY(_s_, _d_)  SL_ARRAY4_COPY(_s_, _d_), ((_d_)[4] = (_s_)[4])
# define SL_ARRAY6_COPY(_s_, _d_)  SL_ARRAY5_COPY(_s_, _d_), ((_d_)[5] = (_s_)[5])
# define SL_ARRAY7_COPY(_s_, _d_)  SL_ARRAY6_COPY(_s_, _d_), ((_d_)[6] = (_s_)[6])
# define SL_ARRAY8_COPY(_s_, _d_)  SL_ARRAY7_COPY(_s_, _d_), ((_d_)[7] = (_s_)[7])
# define SL_ARRAY9_COPY(_s_, _d_)  SL_ARRAY8_COPY(_s_, _d_), ((_d_)[8] = (_s_)[8])
#endif





/* sl_macro fcs_front_xq_aI_sl_key_type_c fcs_front_xq_aI_sl_key_type_mpi fcs_front_xq_aI_sl_key_size_mpi fcs_front_xq_aI_sl_key_type_fmt fcs_front_xq_aI_sl_key_integer fcs_front_xq_aI_sl_key_memcpy */


#define fcs_front_xq_aI_sl_key_byte                          ((fcs_front_xq_aI_slint_t) sizeof(fcs_front_xq_aI_sl_key_type_c))  /* sl_macro */

#ifndef fcs_front_xq_aI_sl_key_copy
 #ifndef fcs_front_xq_aI_sl_key_memcpy
  #define fcs_front_xq_aI_sl_key_copy(src, dst)              SL_ARRAY1_COPY(src, dst)
 #else
  #define fcs_front_xq_aI_sl_key_copy(src, dst)              memcpy(dst, src, fcs_front_xq_aI_sl_key_byte)  /* sl_macro */
 #endif
#endif
#ifndef fcs_front_xq_aI_sl_key_ncopy
 #define fcs_front_xq_aI_sl_key_ncopy(src, dst, n)           memcpy(dst, src, (n) * fcs_front_xq_aI_sl_key_byte)  /* sl_macro */
#endif
#ifndef fcs_front_xq_aI_sl_key_nmove
 #define fcs_front_xq_aI_sl_key_nmove(src, dst, n)           memmove(dst, src, (n) * fcs_front_xq_aI_sl_key_byte)  /* sl_macro */
#endif


#define fcs_front_xq_aI_key_type_c                           fcs_front_xq_aI_sl_key_type_c  /* sl_macro */
#define fcs_front_xq_aI_key_type_mpi                         (fcs_front_xq_aI_sl_key_type_mpi)  /* sl_macro */
#define fcs_front_xq_aI_key_size_mpi                         (fcs_front_xq_aI_sl_key_size_mpi)  /* sl_macro */
#ifdef fcs_front_xq_aI_sl_key_type_fmt
# define fcs_front_xq_aI_key_type_fmt                        fcs_front_xq_aI_sl_key_type_fmt  /* sl_macro */
#endif
#define fcs_front_xq_aI_key_integer                          fcs_front_xq_aI_sl_key_integer  /* sl_macro */

#define fcs_front_xq_aI_key_pure_type_c                      fcs_front_xq_aI_sl_key_pure_type_c  /* sl_macro */
#define fcs_front_xq_aI_key_pure_type_mpi                    (fcs_front_xq_aI_sl_key_pure_type_mpi)  /* sl_macro */
#define fcs_front_xq_aI_key_pure_size_mpi                    (fcs_front_xq_aI_sl_key_pure_size_mpi)  /* sl_macro */
#ifdef fcs_front_xq_aI_sl_key_pure_type_fmt
# define fcs_front_xq_aI_key_pure_type_fmt                   fcs_front_xq_aI_sl_key_pure_type_fmt  /* sl_macro */
#endif

#define fcs_front_xq_aI_key_purify(k)                        (fcs_front_xq_aI_sl_key_purify(k))  /* sl_macro */
#define fcs_front_xq_aI_key_get_pure(k)                      (fcs_front_xq_aI_sl_key_get_pure(k))  /* sl_macro */
#define fcs_front_xq_aI_key_set_pure(k, p)                   (fcs_front_xq_aI_sl_key_set_pure(k, p))  /* sl_macro */

#ifdef fcs_front_xq_aI_key_integer
# define fcs_front_xq_aI_key_integer_unsigned                (((fcs_front_xq_aI_key_pure_type_c) ~((fcs_front_xq_aI_key_pure_type_c) 0)) >= ((fcs_front_xq_aI_key_pure_type_c) 0))  /* sl_macro */
#endif

#define fcs_front_xq_aI_key_n                                1  /* sl_macro */
#define fcs_front_xq_aI_key_byte                             (fcs_front_xq_aI_sl_key_byte)  /* sl_macro */

#define fcs_front_xq_aI_key_cmp_eq(k0, k1)                   (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_cmp_eq((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_cmp_ne(k0, k1)                   (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_cmp_ne((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_cmp_lt(k0, k1)                   (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_cmp_lt((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_cmp_le(k0, k1)                   (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_cmp_le((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_cmp_gt(k0, k1)                   (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_cmp_gt((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_cmp_ge(k0, k1)                   (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_cmp_ge((k0), (k1)))  /* sl_macro */

#define fcs_front_xq_aI_key_pure_cmp_eq(k0, k1)              (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_pure_cmp_eq((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_pure_cmp_ne(k0, k1)              (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_pure_cmp_ne((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_pure_cmp_lt(k0, k1)              (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_pure_cmp_lt((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_pure_cmp_le(k0, k1)              (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_pure_cmp_le((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_pure_cmp_gt(k0, k1)              (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_pure_cmp_gt((k0), (k1)))  /* sl_macro */
#define fcs_front_xq_aI_key_pure_cmp_ge(k0, k1)              (fcs_front_xq_aI_cc_rti_cadd_cmp(1) fcs_front_xq_aI_sl_key_pure_cmp_ge((k0), (k1)))  /* sl_macro */

#ifdef fcs_front_xq_aI_sl_key_val_srand
# define fcs_front_xq_aI_key_val_srand(_s_)                  fcs_front_xq_aI_sl_key_val_srand(_s_)  /* sl_macro */
#endif
#ifdef fcs_front_xq_aI_sl_key_val_rand
# define fcs_front_xq_aI_key_val_rand()                      fcs_front_xq_aI_sl_key_val_rand()  /* sl_macro */
# define fcs_front_xq_aI_have_key_val_rand                   1  /* sl_macro */
#else
# define fcs_front_xq_aI_key_val_rand()                      Z_NOP()
# define fcs_front_xq_aI_have_key_val_rand                   0
#endif
#ifdef fcs_front_xq_aI_sl_key_val_rand_minmax
# define fcs_front_xq_aI_key_val_rand_minmax(_min_, _max_)   fcs_front_xq_aI_sl_key_val_rand_minmax(_min_, _max_)  /* sl_macro */
# define fcs_front_xq_aI_have_key_val_rand_minmax            1  /* sl_macro */
#else
# define fcs_front_xq_aI_key_val_rand_minmax(_min_, _max_)   Z_NOP()
# define fcs_front_xq_aI_have_key_val_rand_minmax            0
#endif


#define fcs_front_xq_aI_key_at(_s_, _sat_)                   ((_s_) + (_sat_))  /* sl_macro */

#define fcs_front_xq_aI_key_assign(src, dst)                 (dst = src)  /* sl_macro */
#define fcs_front_xq_aI_key_assign_at(src, sat, dst)         (dst = &src[sat])  /* sl_macro */
#define fcs_front_xq_aI_key_null(k)                          (k = NULL)  /* sl_macro */
#define fcs_front_xq_aI_key_inc(k)                           (++k)  /* sl_macro */
#define fcs_front_xq_aI_key_dec(k)                           (--k)  /* sl_macro */
#define fcs_front_xq_aI_key_add(k, n)                        (k += n)  /* sl_macro */
#define fcs_front_xq_aI_key_sub(k, n)                        (k -= n)  /* sl_macro */

#define fcs_front_xq_aI_key_copy(src, dst)                   (fcs_front_xq_aI_cc_rti_cadd_movek(1) fcs_front_xq_aI_sl_key_copy(src, dst))  /* sl_macro */
#define fcs_front_xq_aI_key_ncopy(src, dst, n)               (fcs_front_xq_aI_cc_rti_cadd_movek(n) fcs_front_xq_aI_sl_key_ncopy(src, dst, n))  /* sl_macro */
#define fcs_front_xq_aI_key_nmove(src, dst, n)               (fcs_front_xq_aI_cc_rti_cadd_movek(n) fcs_front_xq_aI_sl_key_nmove(src, dst, n))  /* sl_macro */

#define fcs_front_xq_aI_key_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_key_copy(&(src)[sat], &(dst)[dat])  /* sl_macro */
#define fcs_front_xq_aI_key_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_key_ncopy(&(src)[sat], &(dst)[dat], n)  /* sl_macro */
#define fcs_front_xq_aI_key_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_key_nmove(&(src)[sat], &(dst)[dat], n)  /* sl_macro */

#define fcs_front_xq_aI_key_xchange(k0, k1, t)               (fcs_front_xq_aI_key_copy(k0, t), fcs_front_xq_aI_key_copy(k1, k0), fcs_front_xq_aI_key_copy(t, k1))  /* sl_macro */
#define fcs_front_xq_aI_key_xchange_at(k0, at0, k1, at1, t)  (fcs_front_xq_aI_key_copy_at(k0, at0, t, 0), fcs_front_xq_aI_key_copy_at(k1, at1, k0, at0), fcs_front_xq_aI_key_copy_at(t, 0, k1, at1))  /* sl_macro */

#define fcs_front_xq_aI_key_cm                               SLCM_KEYS  /* sl_macro */

#ifdef fcs_front_xq_aI_key_integer
# define fcs_front_xq_aI_key_radix_low                       ((fcs_front_xq_aI_slint_t) 0)  /* sl_macro */
# define fcs_front_xq_aI_key_radix_high                      ((fcs_front_xq_aI_slint_t) (sizeof(fcs_front_xq_aI_key_pure_type_c) * 8 - 1))  /* sl_macro */
# define fcs_front_xq_aI_key_radix_key2class(_k_, _x_, _y_)  (((_k_) >> (_x_)) & (_y_))  /* sl_macro */
#endif







/* sl_macro fcs_front_xq_aI_SL_INDEX fcs_front_xq_aI_SL_PACKED_INDEX fcs_front_xq_aI_sl_index_type_c fcs_front_xq_aI_sl_index_type_mpi fcs_front_xq_aI_sl_index_size_mpi fcs_front_xq_aI_sl_index_type_fmt fcs_front_xq_aI_sl_index_integer fcs_front_xq_aI_sl_index_memcpy */

#ifdef fcs_front_xq_aI_SL_INDEX

# define fcs_front_xq_aI_sl_index_byte                                       ((fcs_front_xq_aI_slint_t) sizeof(fcs_front_xq_aI_sl_index_type_c))  /* sl_macro */

# ifndef fcs_front_xq_aI_sl_index_copy
#  ifndef fcs_front_xq_aI_sl_index_memcpy
#   define fcs_front_xq_aI_sl_index_copy(_s_, _d_)                           SL_ARRAY1_COPY(_s_, _d_)
#  else
#   define fcs_front_xq_aI_sl_index_copy(_s_, _d_)                           memcpy(_d_, _s_, fcs_front_xq_aI_sl_index_byte)  /* sl_macro */
#  endif
# endif
# ifndef fcs_front_xq_aI_sl_index_ncopy
#  define fcs_front_xq_aI_sl_index_ncopy(_s_, _d_, _n_)                      memcpy(_d_, _s_, (_n_) * fcs_front_xq_aI_sl_index_byte)  /* sl_macro */
# endif
# ifndef fcs_front_xq_aI_sl_index_nmove
#  define fcs_front_xq_aI_sl_index_nmove(_s_, _d_, _n_)                      memmove(_d_, _s_, (_n_) * fcs_front_xq_aI_sl_index_byte)  /* sl_macro */
# endif


# define fcs_front_xq_aI_index_type_c                                        fcs_front_xq_aI_sl_index_type_c  /* sl_macro */
# define fcs_front_xq_aI_index_type_mpi                                      (fcs_front_xq_aI_sl_index_type_mpi)  /* sl_macro */
# define fcs_front_xq_aI_index_size_mpi                                      (fcs_front_xq_aI_sl_index_size_mpi)  /* sl_macro */
# define fcs_front_xq_aI_index_type_fmt                                      fcs_front_xq_aI_sl_index_type_fmt  /* sl_macro */

# define fcs_front_xq_aI_index_n                                             1  /* sl_macro */
# define fcs_front_xq_aI_index_byte                                          (fcs_front_xq_aI_sl_index_byte)  /* sl_macro */

/* commands for regular use */
# define fcs_front_xq_aI_index_assign(_s_, _d_)                              (_d_ = _s_)  /* sl_macro */
# define fcs_front_xq_aI_index_assign_at(_s_, _sat_, _d_)                    (_d_ = &_s_[_sat_])  /* sl_macro */
# define fcs_front_xq_aI_index_null(_i_)                                     (_i_ = NULL)  /* sl_macro */
# define fcs_front_xq_aI_index_inc(_i_)                                      (++_i_)  /* sl_macro */
# define fcs_front_xq_aI_index_dec(_i_)                                      (--_i_)  /* sl_macro */
# define fcs_front_xq_aI_index_add(_i_, _n_)                                 (_i_ += _n_)  /* sl_macro */
# define fcs_front_xq_aI_index_sub(_i_, _n_)                                 (_i_ -= _n_)  /* sl_macro */

# define fcs_front_xq_aI_index_copy(_s_, _d_)                                fcs_front_xq_aI_sl_index_copy(_s_, _d_)  /* sl_macro */
# define fcs_front_xq_aI_index_ncopy(_s_, _d_, _n_)                          fcs_front_xq_aI_sl_index_ncopy(_s_, _d_, _n_)  /* sl_macro */
# define fcs_front_xq_aI_index_nmove(_s_, _d_, _n_)                          fcs_front_xq_aI_sl_index_nmove(_s_, _d_, _n_)  /* sl_macro */

# define fcs_front_xq_aI_index_copy_at(_s_, _sat_, _d_, _dat_)               fcs_front_xq_aI_index_copy(&(_s_)[_sat_], &(_d_)[_dat_])  /* sl_macro */
# define fcs_front_xq_aI_index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)         fcs_front_xq_aI_index_ncopy(&(_s_)[_sat_], &(_d_)[_dat_], _n_)  /* sl_macro */
# define fcs_front_xq_aI_index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)         fcs_front_xq_aI_index_nmove(&(_s_)[_sat_], &(_d_)[_dat_], _n_)  /* sl_macro */

# define fcs_front_xq_aI_index_xchange(_i0_, _i1_, _t_)                      (fcs_front_xq_aI_index_copy(_i0_, _t_), fcs_front_xq_aI_index_copy(_i1_, _i0_), fcs_front_xq_aI_index_copy(_t_, _i1_))  /* sl_macro */
# define fcs_front_xq_aI_index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)     (fcs_front_xq_aI_index_copy_at(_i0_, _at0_, _t_, 0), fcs_front_xq_aI_index_copy_at(_i1_, _at1_, _i0_, _at0_), fcs_front_xq_aI_index_copy_at(_t_, 0, _i1_, _at1_))  /* sl_macro */

/* chained command versions */
# define fcs_front_xq_aI_cc_index_assign(_s_, _d_)                           , fcs_front_xq_aI_index_assign(_s_, _d_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_assign_at(_s_, _sat_, _d_)                 , fcs_front_xq_aI_index_assign_at(_s_, _sat_, _d_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_null(_i_)                                  , fcs_front_xq_aI_index_null(_i_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_inc(_i_)                                   , fcs_front_xq_aI_index_inc(_i_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_dec(_i_)                                   , fcs_front_xq_aI_index_dec(_i_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_add(_i_, _n_)                              , fcs_front_xq_aI_index_add(_i_, _n_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_sub(_i_, _n_)                              , fcs_front_xq_aI_index_sub(_i_, _n_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_copy(_s_, _d_)                             , fcs_front_xq_aI_index_copy(_s_, _d_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_ncopy(_s_, _d_, _n_)                       , fcs_front_xq_aI_index_ncopy(_s_, _d_, _n_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_nmove(_s_, _d_, _n_)                       , fcs_front_xq_aI_index_nmove(_s_, _d_, _n_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_copy_at(_s_, _sat_, _d_, _dat_)            , fcs_front_xq_aI_index_copy_at(_s_, _sat_, _d_, _dat_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      , fcs_front_xq_aI_index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      , fcs_front_xq_aI_index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_xchange(_i0_, _i1_, _t_)                   , fcs_front_xq_aI_index_xchange(_i0_, _i1_, _t_)  /* sl_macro */
# define fcs_front_xq_aI_cc_index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)  , fcs_front_xq_aI_index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)  /* sl_macro */

#else /* fcs_front_xq_aI_SL_INDEX */

# define fcs_front_xq_aI_index_n                                             0
# define fcs_front_xq_aI_index_byte                                          0

/* commands for regular use */
# define fcs_front_xq_aI_index_assign(_s_, _d_)                              Z_NOP()
# define fcs_front_xq_aI_index_assign_at(_s_, _sat_, _d_)                    Z_NOP()
# define fcs_front_xq_aI_index_null(_i_)                                     Z_NOP()
# define fcs_front_xq_aI_index_inc(_i_)                                      Z_NOP()
# define fcs_front_xq_aI_index_dec(_i_)                                      Z_NOP()
# define fcs_front_xq_aI_index_add(_i_, _n_)                                 Z_NOP()
# define fcs_front_xq_aI_index_sub(_i_, _n_)                                 Z_NOP()
# define fcs_front_xq_aI_index_copy(_s_, _d_)                                Z_NOP()
# define fcs_front_xq_aI_index_ncopy(_s_, _d_, _n_)                          Z_NOP()
# define fcs_front_xq_aI_index_nmove(_s_, _d_, _n_)                          Z_NOP()
# define fcs_front_xq_aI_index_copy_at(_s_, _sat_, _d_, _dat_)               Z_NOP()
# define fcs_front_xq_aI_index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)         Z_NOP()
# define fcs_front_xq_aI_index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)         Z_NOP()
# define fcs_front_xq_aI_index_xchange(_i0_, _i1_, _t_)                      Z_NOP()
# define fcs_front_xq_aI_index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)     Z_NOP()

/* chained command versions */
# define fcs_front_xq_aI_cc_index_assign(_s_, _d_)
# define fcs_front_xq_aI_cc_index_assign_at(_s_, _sat_, _d_)
# define fcs_front_xq_aI_cc_index_null(_i_)
# define fcs_front_xq_aI_cc_index_inc(_i_)
# define fcs_front_xq_aI_cc_index_dec(_i_)
# define fcs_front_xq_aI_cc_index_add(_i_, _n_)
# define fcs_front_xq_aI_cc_index_sub(_i_, _n_)
# define fcs_front_xq_aI_cc_index_copy(_s_, _d_)
# define fcs_front_xq_aI_cc_index_ncopy(_s_, _d_, _n_)
# define fcs_front_xq_aI_cc_index_nmove(_s_, _d_, _n_)
# define fcs_front_xq_aI_cc_index_copy_at(_s_, _sat_, _d_, _dat_)
# define fcs_front_xq_aI_cc_index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
# define fcs_front_xq_aI_cc_index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
# define fcs_front_xq_aI_cc_index_xchange(_i0_, _i1_, _t_)
# define fcs_front_xq_aI_cc_index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)

#endif /* fcs_front_xq_aI_SL_INDEX */

#define fcs_front_xq_aI_index_cm                                             SLCM_INDICES  /* sl_macro */



#undef SL_DATA





/* DATAX_TEMPLATE_BEGIN */

/* sl_macro fcs_front_xq_aI_SL_DATA0 fcs_front_xq_aI_SL_DATA0_IGNORE fcs_front_xq_aI_sl_data0_type_c fcs_front_xq_aI_sl_data0_size_c fcs_front_xq_aI_sl_data0_type_mpi fcs_front_xq_aI_sl_data0_size_mpi fcs_front_xq_aI_sl_data0_memcpy fcs_front_xq_aI_sl_data0_weight fcs_front_xq_aI_sl_data0_flex */

#ifdef fcs_front_xq_aI_SL_DATA0

 #define fcs_front_xq_aI_sl_data0_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data0_size_c) * sizeof(fcs_front_xq_aI_sl_data0_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data0_copy
  #if fcs_front_xq_aI_sl_data0_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data0_memcpy)
   #if fcs_front_xq_aI_sl_data0_size_c == 1
    #define fcs_front_xq_aI_sl_data0_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data0_size_c == 2
    #define fcs_front_xq_aI_sl_data0_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data0_size_c == 3
    #define fcs_front_xq_aI_sl_data0_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data0_size_c == 4
    #define fcs_front_xq_aI_sl_data0_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data0_size_c == 5
    #define fcs_front_xq_aI_sl_data0_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data0_size_c == 6
    #define fcs_front_xq_aI_sl_data0_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data0_size_c == 7
    #define fcs_front_xq_aI_sl_data0_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data0_size_c == 8
    #define fcs_front_xq_aI_sl_data0_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data0_size_c == 9
    #define fcs_front_xq_aI_sl_data0_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data0_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data0_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data0_ncopy
  #define fcs_front_xq_aI_sl_data0_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data0_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data0_nmove
  #define fcs_front_xq_aI_sl_data0_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data0_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data0_type_c                              fcs_front_xq_aI_sl_data0_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data0_size_c                              (fcs_front_xq_aI_sl_data0_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data0_type_mpi                            (fcs_front_xq_aI_sl_data0_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data0_size_mpi                            (fcs_front_xq_aI_sl_data0_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data0_idx                                 0  /* sl_macro */

 #define fcs_front_xq_aI_data0_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data0_byte                                (fcs_front_xq_aI_sl_data0_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data0_ptr(e)                              (e)->data0  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data0_flex
 # define fcs_front_xq_aI_data0_byte_flex                          (fcs_front_xq_aI_sl_data0_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data0_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data0_weight
 # define fcs_front_xq_aI_data0_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data0_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data0_assign(src, dst)                    ((dst)->data0 = (src)->data0)  /* sl_macro */
 #define fcs_front_xq_aI_data0_assign_at(src, sat, dst)            ((dst)->data0 = &(src)->data0[(sat) * fcs_front_xq_aI_data0_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data0_null(e)                             ((e)->data0 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data0_inc(e)                              ((e)->data0 += fcs_front_xq_aI_data0_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data0_dec(e)                              ((e)->data0 -= fcs_front_xq_aI_data0_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data0_add(e, n)                           ((e)->data0 += (n) * fcs_front_xq_aI_data0_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data0_sub(e, n)                           ((e)->data0 -= (n) * fcs_front_xq_aI_data0_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data0_copy(src, dst)                      fcs_front_xq_aI_sl_data0_copy((src)->data0, (dst)->data0)  /* sl_macro */
 #define fcs_front_xq_aI_data0_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data0_ncopy((src)->data0, (dst)->data0, n)  /* sl_macro */
 #define fcs_front_xq_aI_data0_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data0_nmove((src)->data0, (dst)->data0, n)  /* sl_macro */

 #define fcs_front_xq_aI_data0_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data0_copy(&(src)->data0[(sat) * fcs_front_xq_aI_data0_size_c], &(dst)->data0[(dat) * fcs_front_xq_aI_data0_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data0_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data0_ncopy(&(src)->data0[(sat) * fcs_front_xq_aI_data0_size_c], &(dst)->data0[(dat) * fcs_front_xq_aI_data0_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data0_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data0_nmove(&(src)->data0[(sat) * fcs_front_xq_aI_data0_size_c], &(dst)->data0[(dat) * fcs_front_xq_aI_data0_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data0_xchange(e0, e1, t)                  (fcs_front_xq_aI_data0_copy(e0, t), fcs_front_xq_aI_data0_copy(e1, e0), fcs_front_xq_aI_data0_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data0_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data0_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data0_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data0_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data0_assign(src, dst)                 , fcs_front_xq_aI_data0_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_assign_at(src, sat, dst)         , fcs_front_xq_aI_data0_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_null(e)                          , fcs_front_xq_aI_data0_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_inc(e)                           , fcs_front_xq_aI_data0_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_dec(e)                           , fcs_front_xq_aI_data0_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_add(e, n)                        , fcs_front_xq_aI_data0_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_sub(e, n)                        , fcs_front_xq_aI_data0_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_copy(src, dst)                   , fcs_front_xq_aI_data0_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_ncopy(src, dst, n)               , fcs_front_xq_aI_data0_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_nmove(src, dst, n)               , fcs_front_xq_aI_data0_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data0_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data0_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data0_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_xchange(e0, e1, t)               , fcs_front_xq_aI_data0_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data0_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data0_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data0_assign(src, dst)                 fcs_front_xq_aI_data0_assign(src, dst)
 #define fcs_front_xq_aI_cc_data0_assign_at(src, sat, dst)         fcs_front_xq_aI_data0_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data0_null(e)                          fcs_front_xq_aI_data0_null(e)
 #define fcs_front_xq_aI_cc_data0_inc(e)                           fcs_front_xq_aI_data0_inc(e)
 #define fcs_front_xq_aI_cc_data0_dec(e)                           fcs_front_xq_aI_data0_dec(e)
 #define fcs_front_xq_aI_cc_data0_add(e, n)                        fcs_front_xq_aI_data0_add(e, n)
 #define fcs_front_xq_aI_cc_data0_sub(e, n)                        fcs_front_xq_aI_data0_sub(e, n)
 #define fcs_front_xq_aI_cc_data0_copy(src, dst)                   fcs_front_xq_aI_data0_copy(src, dst)
 #define fcs_front_xq_aI_cc_data0_ncopy(src, dst, n)               fcs_front_xq_aI_data0_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data0_nmove(src, dst, n)               fcs_front_xq_aI_data0_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data0_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data0_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data0_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data0_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data0_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data0_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data0_xchange(e0, e1, t)               fcs_front_xq_aI_data0_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data0_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data0_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA0 */

 #define fcs_front_xq_aI_data0_n                                   0
 #define fcs_front_xq_aI_data0_byte                                0
/* #define fcs_front_xq_aI_data0_ptr(e)*/

 #define fcs_front_xq_aI_data0_byte_flex                           0
 #define fcs_front_xq_aI_data0_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data0_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data0_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data0_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data0_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data0_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data0_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data0_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data0_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data0_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data0_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data0_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data0_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data0_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data0_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data0_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data0_assign(src, dst)
 #define fcs_front_xq_aI_cc_data0_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data0_null(e)
 #define fcs_front_xq_aI_cc_data0_inc(e)
 #define fcs_front_xq_aI_cc_data0_dec(e)
 #define fcs_front_xq_aI_cc_data0_add(e, n)
 #define fcs_front_xq_aI_cc_data0_sub(e, n)
 #define fcs_front_xq_aI_cc_data0_copy(src, dst)
 #define fcs_front_xq_aI_cc_data0_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data0_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data0_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data0_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data0_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data0_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data0_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA0 */

#define fcs_front_xq_aI_data0_cm                                   SLCM_DATA0  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA1 fcs_front_xq_aI_SL_DATA1_IGNORE fcs_front_xq_aI_sl_data1_type_c fcs_front_xq_aI_sl_data1_size_c fcs_front_xq_aI_sl_data1_type_mpi fcs_front_xq_aI_sl_data1_size_mpi fcs_front_xq_aI_sl_data1_memcpy fcs_front_xq_aI_sl_data1_weight fcs_front_xq_aI_sl_data1_flex */

#ifdef fcs_front_xq_aI_SL_DATA1

 #define fcs_front_xq_aI_sl_data1_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data1_size_c) * sizeof(fcs_front_xq_aI_sl_data1_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data1_copy
  #if fcs_front_xq_aI_sl_data1_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data1_memcpy)
   #if fcs_front_xq_aI_sl_data1_size_c == 1
    #define fcs_front_xq_aI_sl_data1_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data1_size_c == 2
    #define fcs_front_xq_aI_sl_data1_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data1_size_c == 3
    #define fcs_front_xq_aI_sl_data1_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data1_size_c == 4
    #define fcs_front_xq_aI_sl_data1_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data1_size_c == 5
    #define fcs_front_xq_aI_sl_data1_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data1_size_c == 6
    #define fcs_front_xq_aI_sl_data1_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data1_size_c == 7
    #define fcs_front_xq_aI_sl_data1_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data1_size_c == 8
    #define fcs_front_xq_aI_sl_data1_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data1_size_c == 9
    #define fcs_front_xq_aI_sl_data1_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data1_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data1_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data1_ncopy
  #define fcs_front_xq_aI_sl_data1_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data1_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data1_nmove
  #define fcs_front_xq_aI_sl_data1_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data1_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data1_type_c                              fcs_front_xq_aI_sl_data1_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data1_size_c                              (fcs_front_xq_aI_sl_data1_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data1_type_mpi                            (fcs_front_xq_aI_sl_data1_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data1_size_mpi                            (fcs_front_xq_aI_sl_data1_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data1_idx                                 1  /* sl_macro */

 #define fcs_front_xq_aI_data1_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data1_byte                                (fcs_front_xq_aI_sl_data1_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data1_ptr(e)                              (e)->data1  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data1_flex
 # define fcs_front_xq_aI_data1_byte_flex                          (fcs_front_xq_aI_sl_data1_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data1_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data1_weight
 # define fcs_front_xq_aI_data1_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data1_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data1_assign(src, dst)                    ((dst)->data1 = (src)->data1)  /* sl_macro */
 #define fcs_front_xq_aI_data1_assign_at(src, sat, dst)            ((dst)->data1 = &(src)->data1[(sat) * fcs_front_xq_aI_data1_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data1_null(e)                             ((e)->data1 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data1_inc(e)                              ((e)->data1 += fcs_front_xq_aI_data1_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data1_dec(e)                              ((e)->data1 -= fcs_front_xq_aI_data1_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data1_add(e, n)                           ((e)->data1 += (n) * fcs_front_xq_aI_data1_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data1_sub(e, n)                           ((e)->data1 -= (n) * fcs_front_xq_aI_data1_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data1_copy(src, dst)                      fcs_front_xq_aI_sl_data1_copy((src)->data1, (dst)->data1)  /* sl_macro */
 #define fcs_front_xq_aI_data1_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data1_ncopy((src)->data1, (dst)->data1, n)  /* sl_macro */
 #define fcs_front_xq_aI_data1_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data1_nmove((src)->data1, (dst)->data1, n)  /* sl_macro */

 #define fcs_front_xq_aI_data1_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data1_copy(&(src)->data1[(sat) * fcs_front_xq_aI_data1_size_c], &(dst)->data1[(dat) * fcs_front_xq_aI_data1_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data1_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data1_ncopy(&(src)->data1[(sat) * fcs_front_xq_aI_data1_size_c], &(dst)->data1[(dat) * fcs_front_xq_aI_data1_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data1_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data1_nmove(&(src)->data1[(sat) * fcs_front_xq_aI_data1_size_c], &(dst)->data1[(dat) * fcs_front_xq_aI_data1_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data1_xchange(e0, e1, t)                  (fcs_front_xq_aI_data1_copy(e0, t), fcs_front_xq_aI_data1_copy(e1, e0), fcs_front_xq_aI_data1_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data1_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data1_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data1_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data1_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data1_assign(src, dst)                 , fcs_front_xq_aI_data1_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_assign_at(src, sat, dst)         , fcs_front_xq_aI_data1_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_null(e)                          , fcs_front_xq_aI_data1_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_inc(e)                           , fcs_front_xq_aI_data1_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_dec(e)                           , fcs_front_xq_aI_data1_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_add(e, n)                        , fcs_front_xq_aI_data1_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_sub(e, n)                        , fcs_front_xq_aI_data1_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_copy(src, dst)                   , fcs_front_xq_aI_data1_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_ncopy(src, dst, n)               , fcs_front_xq_aI_data1_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_nmove(src, dst, n)               , fcs_front_xq_aI_data1_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data1_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data1_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data1_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_xchange(e0, e1, t)               , fcs_front_xq_aI_data1_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data1_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data1_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data1_assign(src, dst)                 fcs_front_xq_aI_data1_assign(src, dst)
 #define fcs_front_xq_aI_cc_data1_assign_at(src, sat, dst)         fcs_front_xq_aI_data1_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data1_null(e)                          fcs_front_xq_aI_data1_null(e)
 #define fcs_front_xq_aI_cc_data1_inc(e)                           fcs_front_xq_aI_data1_inc(e)
 #define fcs_front_xq_aI_cc_data1_dec(e)                           fcs_front_xq_aI_data1_dec(e)
 #define fcs_front_xq_aI_cc_data1_add(e, n)                        fcs_front_xq_aI_data1_add(e, n)
 #define fcs_front_xq_aI_cc_data1_sub(e, n)                        fcs_front_xq_aI_data1_sub(e, n)
 #define fcs_front_xq_aI_cc_data1_copy(src, dst)                   fcs_front_xq_aI_data1_copy(src, dst)
 #define fcs_front_xq_aI_cc_data1_ncopy(src, dst, n)               fcs_front_xq_aI_data1_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data1_nmove(src, dst, n)               fcs_front_xq_aI_data1_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data1_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data1_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data1_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data1_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data1_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data1_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data1_xchange(e0, e1, t)               fcs_front_xq_aI_data1_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data1_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data1_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA1 */

 #define fcs_front_xq_aI_data1_n                                   0
 #define fcs_front_xq_aI_data1_byte                                0
/* #define fcs_front_xq_aI_data1_ptr(e)*/

 #define fcs_front_xq_aI_data1_byte_flex                           0
 #define fcs_front_xq_aI_data1_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data1_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data1_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data1_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data1_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data1_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data1_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data1_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data1_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data1_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data1_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data1_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data1_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data1_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data1_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data1_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data1_assign(src, dst)
 #define fcs_front_xq_aI_cc_data1_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data1_null(e)
 #define fcs_front_xq_aI_cc_data1_inc(e)
 #define fcs_front_xq_aI_cc_data1_dec(e)
 #define fcs_front_xq_aI_cc_data1_add(e, n)
 #define fcs_front_xq_aI_cc_data1_sub(e, n)
 #define fcs_front_xq_aI_cc_data1_copy(src, dst)
 #define fcs_front_xq_aI_cc_data1_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data1_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data1_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data1_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data1_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data1_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data1_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA1 */

#define fcs_front_xq_aI_data1_cm                                   SLCM_DATA1  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA2 fcs_front_xq_aI_SL_DATA2_IGNORE fcs_front_xq_aI_sl_data2_type_c fcs_front_xq_aI_sl_data2_size_c fcs_front_xq_aI_sl_data2_type_mpi fcs_front_xq_aI_sl_data2_size_mpi fcs_front_xq_aI_sl_data2_memcpy fcs_front_xq_aI_sl_data2_weight fcs_front_xq_aI_sl_data2_flex */

#ifdef fcs_front_xq_aI_SL_DATA2

 #define fcs_front_xq_aI_sl_data2_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data2_size_c) * sizeof(fcs_front_xq_aI_sl_data2_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data2_copy
  #if fcs_front_xq_aI_sl_data2_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data2_memcpy)
   #if fcs_front_xq_aI_sl_data2_size_c == 1
    #define fcs_front_xq_aI_sl_data2_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data2_size_c == 2
    #define fcs_front_xq_aI_sl_data2_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data2_size_c == 3
    #define fcs_front_xq_aI_sl_data2_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data2_size_c == 4
    #define fcs_front_xq_aI_sl_data2_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data2_size_c == 5
    #define fcs_front_xq_aI_sl_data2_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data2_size_c == 6
    #define fcs_front_xq_aI_sl_data2_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data2_size_c == 7
    #define fcs_front_xq_aI_sl_data2_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data2_size_c == 8
    #define fcs_front_xq_aI_sl_data2_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data2_size_c == 9
    #define fcs_front_xq_aI_sl_data2_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data2_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data2_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data2_ncopy
  #define fcs_front_xq_aI_sl_data2_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data2_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data2_nmove
  #define fcs_front_xq_aI_sl_data2_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data2_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data2_type_c                              fcs_front_xq_aI_sl_data2_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data2_size_c                              (fcs_front_xq_aI_sl_data2_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data2_type_mpi                            (fcs_front_xq_aI_sl_data2_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data2_size_mpi                            (fcs_front_xq_aI_sl_data2_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data2_idx                                 2  /* sl_macro */

 #define fcs_front_xq_aI_data2_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data2_byte                                (fcs_front_xq_aI_sl_data2_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data2_ptr(e)                              (e)->data2  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data2_flex
 # define fcs_front_xq_aI_data2_byte_flex                          (fcs_front_xq_aI_sl_data2_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data2_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data2_weight
 # define fcs_front_xq_aI_data2_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data2_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data2_assign(src, dst)                    ((dst)->data2 = (src)->data2)  /* sl_macro */
 #define fcs_front_xq_aI_data2_assign_at(src, sat, dst)            ((dst)->data2 = &(src)->data2[(sat) * fcs_front_xq_aI_data2_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data2_null(e)                             ((e)->data2 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data2_inc(e)                              ((e)->data2 += fcs_front_xq_aI_data2_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data2_dec(e)                              ((e)->data2 -= fcs_front_xq_aI_data2_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data2_add(e, n)                           ((e)->data2 += (n) * fcs_front_xq_aI_data2_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data2_sub(e, n)                           ((e)->data2 -= (n) * fcs_front_xq_aI_data2_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data2_copy(src, dst)                      fcs_front_xq_aI_sl_data2_copy((src)->data2, (dst)->data2)  /* sl_macro */
 #define fcs_front_xq_aI_data2_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data2_ncopy((src)->data2, (dst)->data2, n)  /* sl_macro */
 #define fcs_front_xq_aI_data2_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data2_nmove((src)->data2, (dst)->data2, n)  /* sl_macro */

 #define fcs_front_xq_aI_data2_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data2_copy(&(src)->data2[(sat) * fcs_front_xq_aI_data2_size_c], &(dst)->data2[(dat) * fcs_front_xq_aI_data2_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data2_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data2_ncopy(&(src)->data2[(sat) * fcs_front_xq_aI_data2_size_c], &(dst)->data2[(dat) * fcs_front_xq_aI_data2_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data2_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data2_nmove(&(src)->data2[(sat) * fcs_front_xq_aI_data2_size_c], &(dst)->data2[(dat) * fcs_front_xq_aI_data2_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data2_xchange(e0, e1, t)                  (fcs_front_xq_aI_data2_copy(e0, t), fcs_front_xq_aI_data2_copy(e1, e0), fcs_front_xq_aI_data2_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data2_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data2_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data2_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data2_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data2_assign(src, dst)                 , fcs_front_xq_aI_data2_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_assign_at(src, sat, dst)         , fcs_front_xq_aI_data2_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_null(e)                          , fcs_front_xq_aI_data2_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_inc(e)                           , fcs_front_xq_aI_data2_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_dec(e)                           , fcs_front_xq_aI_data2_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_add(e, n)                        , fcs_front_xq_aI_data2_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_sub(e, n)                        , fcs_front_xq_aI_data2_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_copy(src, dst)                   , fcs_front_xq_aI_data2_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_ncopy(src, dst, n)               , fcs_front_xq_aI_data2_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_nmove(src, dst, n)               , fcs_front_xq_aI_data2_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data2_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data2_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data2_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_xchange(e0, e1, t)               , fcs_front_xq_aI_data2_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data2_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data2_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data2_assign(src, dst)                 fcs_front_xq_aI_data2_assign(src, dst)
 #define fcs_front_xq_aI_cc_data2_assign_at(src, sat, dst)         fcs_front_xq_aI_data2_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data2_null(e)                          fcs_front_xq_aI_data2_null(e)
 #define fcs_front_xq_aI_cc_data2_inc(e)                           fcs_front_xq_aI_data2_inc(e)
 #define fcs_front_xq_aI_cc_data2_dec(e)                           fcs_front_xq_aI_data2_dec(e)
 #define fcs_front_xq_aI_cc_data2_add(e, n)                        fcs_front_xq_aI_data2_add(e, n)
 #define fcs_front_xq_aI_cc_data2_sub(e, n)                        fcs_front_xq_aI_data2_sub(e, n)
 #define fcs_front_xq_aI_cc_data2_copy(src, dst)                   fcs_front_xq_aI_data2_copy(src, dst)
 #define fcs_front_xq_aI_cc_data2_ncopy(src, dst, n)               fcs_front_xq_aI_data2_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data2_nmove(src, dst, n)               fcs_front_xq_aI_data2_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data2_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data2_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data2_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data2_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data2_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data2_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data2_xchange(e0, e1, t)               fcs_front_xq_aI_data2_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data2_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data2_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA2 */

 #define fcs_front_xq_aI_data2_n                                   0
 #define fcs_front_xq_aI_data2_byte                                0
/* #define fcs_front_xq_aI_data2_ptr(e)*/

 #define fcs_front_xq_aI_data2_byte_flex                           0
 #define fcs_front_xq_aI_data2_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data2_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data2_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data2_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data2_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data2_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data2_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data2_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data2_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data2_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data2_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data2_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data2_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data2_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data2_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data2_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data2_assign(src, dst)
 #define fcs_front_xq_aI_cc_data2_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data2_null(e)
 #define fcs_front_xq_aI_cc_data2_inc(e)
 #define fcs_front_xq_aI_cc_data2_dec(e)
 #define fcs_front_xq_aI_cc_data2_add(e, n)
 #define fcs_front_xq_aI_cc_data2_sub(e, n)
 #define fcs_front_xq_aI_cc_data2_copy(src, dst)
 #define fcs_front_xq_aI_cc_data2_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data2_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data2_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data2_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data2_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data2_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data2_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA2 */

#define fcs_front_xq_aI_data2_cm                                   SLCM_DATA2  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA3 fcs_front_xq_aI_SL_DATA3_IGNORE fcs_front_xq_aI_sl_data3_type_c fcs_front_xq_aI_sl_data3_size_c fcs_front_xq_aI_sl_data3_type_mpi fcs_front_xq_aI_sl_data3_size_mpi fcs_front_xq_aI_sl_data3_memcpy fcs_front_xq_aI_sl_data3_weight fcs_front_xq_aI_sl_data3_flex */

#ifdef fcs_front_xq_aI_SL_DATA3

 #define fcs_front_xq_aI_sl_data3_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data3_size_c) * sizeof(fcs_front_xq_aI_sl_data3_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data3_copy
  #if fcs_front_xq_aI_sl_data3_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data3_memcpy)
   #if fcs_front_xq_aI_sl_data3_size_c == 1
    #define fcs_front_xq_aI_sl_data3_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data3_size_c == 2
    #define fcs_front_xq_aI_sl_data3_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data3_size_c == 3
    #define fcs_front_xq_aI_sl_data3_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data3_size_c == 4
    #define fcs_front_xq_aI_sl_data3_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data3_size_c == 5
    #define fcs_front_xq_aI_sl_data3_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data3_size_c == 6
    #define fcs_front_xq_aI_sl_data3_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data3_size_c == 7
    #define fcs_front_xq_aI_sl_data3_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data3_size_c == 8
    #define fcs_front_xq_aI_sl_data3_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data3_size_c == 9
    #define fcs_front_xq_aI_sl_data3_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data3_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data3_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data3_ncopy
  #define fcs_front_xq_aI_sl_data3_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data3_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data3_nmove
  #define fcs_front_xq_aI_sl_data3_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data3_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data3_type_c                              fcs_front_xq_aI_sl_data3_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data3_size_c                              (fcs_front_xq_aI_sl_data3_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data3_type_mpi                            (fcs_front_xq_aI_sl_data3_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data3_size_mpi                            (fcs_front_xq_aI_sl_data3_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data3_idx                                 3  /* sl_macro */

 #define fcs_front_xq_aI_data3_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data3_byte                                (fcs_front_xq_aI_sl_data3_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data3_ptr(e)                              (e)->data3  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data3_flex
 # define fcs_front_xq_aI_data3_byte_flex                          (fcs_front_xq_aI_sl_data3_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data3_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data3_weight
 # define fcs_front_xq_aI_data3_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data3_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data3_assign(src, dst)                    ((dst)->data3 = (src)->data3)  /* sl_macro */
 #define fcs_front_xq_aI_data3_assign_at(src, sat, dst)            ((dst)->data3 = &(src)->data3[(sat) * fcs_front_xq_aI_data3_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data3_null(e)                             ((e)->data3 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data3_inc(e)                              ((e)->data3 += fcs_front_xq_aI_data3_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data3_dec(e)                              ((e)->data3 -= fcs_front_xq_aI_data3_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data3_add(e, n)                           ((e)->data3 += (n) * fcs_front_xq_aI_data3_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data3_sub(e, n)                           ((e)->data3 -= (n) * fcs_front_xq_aI_data3_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data3_copy(src, dst)                      fcs_front_xq_aI_sl_data3_copy((src)->data3, (dst)->data3)  /* sl_macro */
 #define fcs_front_xq_aI_data3_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data3_ncopy((src)->data3, (dst)->data3, n)  /* sl_macro */
 #define fcs_front_xq_aI_data3_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data3_nmove((src)->data3, (dst)->data3, n)  /* sl_macro */

 #define fcs_front_xq_aI_data3_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data3_copy(&(src)->data3[(sat) * fcs_front_xq_aI_data3_size_c], &(dst)->data3[(dat) * fcs_front_xq_aI_data3_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data3_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data3_ncopy(&(src)->data3[(sat) * fcs_front_xq_aI_data3_size_c], &(dst)->data3[(dat) * fcs_front_xq_aI_data3_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data3_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data3_nmove(&(src)->data3[(sat) * fcs_front_xq_aI_data3_size_c], &(dst)->data3[(dat) * fcs_front_xq_aI_data3_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data3_xchange(e0, e1, t)                  (fcs_front_xq_aI_data3_copy(e0, t), fcs_front_xq_aI_data3_copy(e1, e0), fcs_front_xq_aI_data3_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data3_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data3_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data3_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data3_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data3_assign(src, dst)                 , fcs_front_xq_aI_data3_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_assign_at(src, sat, dst)         , fcs_front_xq_aI_data3_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_null(e)                          , fcs_front_xq_aI_data3_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_inc(e)                           , fcs_front_xq_aI_data3_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_dec(e)                           , fcs_front_xq_aI_data3_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_add(e, n)                        , fcs_front_xq_aI_data3_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_sub(e, n)                        , fcs_front_xq_aI_data3_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_copy(src, dst)                   , fcs_front_xq_aI_data3_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_ncopy(src, dst, n)               , fcs_front_xq_aI_data3_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_nmove(src, dst, n)               , fcs_front_xq_aI_data3_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data3_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data3_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data3_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_xchange(e0, e1, t)               , fcs_front_xq_aI_data3_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data3_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data3_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data3_assign(src, dst)                 fcs_front_xq_aI_data3_assign(src, dst)
 #define fcs_front_xq_aI_cc_data3_assign_at(src, sat, dst)         fcs_front_xq_aI_data3_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data3_null(e)                          fcs_front_xq_aI_data3_null(e)
 #define fcs_front_xq_aI_cc_data3_inc(e)                           fcs_front_xq_aI_data3_inc(e)
 #define fcs_front_xq_aI_cc_data3_dec(e)                           fcs_front_xq_aI_data3_dec(e)
 #define fcs_front_xq_aI_cc_data3_add(e, n)                        fcs_front_xq_aI_data3_add(e, n)
 #define fcs_front_xq_aI_cc_data3_sub(e, n)                        fcs_front_xq_aI_data3_sub(e, n)
 #define fcs_front_xq_aI_cc_data3_copy(src, dst)                   fcs_front_xq_aI_data3_copy(src, dst)
 #define fcs_front_xq_aI_cc_data3_ncopy(src, dst, n)               fcs_front_xq_aI_data3_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data3_nmove(src, dst, n)               fcs_front_xq_aI_data3_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data3_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data3_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data3_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data3_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data3_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data3_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data3_xchange(e0, e1, t)               fcs_front_xq_aI_data3_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data3_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data3_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA3 */

 #define fcs_front_xq_aI_data3_n                                   0
 #define fcs_front_xq_aI_data3_byte                                0
/* #define fcs_front_xq_aI_data3_ptr(e)*/

 #define fcs_front_xq_aI_data3_byte_flex                           0
 #define fcs_front_xq_aI_data3_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data3_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data3_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data3_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data3_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data3_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data3_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data3_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data3_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data3_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data3_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data3_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data3_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data3_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data3_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data3_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data3_assign(src, dst)
 #define fcs_front_xq_aI_cc_data3_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data3_null(e)
 #define fcs_front_xq_aI_cc_data3_inc(e)
 #define fcs_front_xq_aI_cc_data3_dec(e)
 #define fcs_front_xq_aI_cc_data3_add(e, n)
 #define fcs_front_xq_aI_cc_data3_sub(e, n)
 #define fcs_front_xq_aI_cc_data3_copy(src, dst)
 #define fcs_front_xq_aI_cc_data3_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data3_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data3_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data3_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data3_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data3_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data3_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA3 */

#define fcs_front_xq_aI_data3_cm                                   SLCM_DATA3  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA4 fcs_front_xq_aI_SL_DATA4_IGNORE fcs_front_xq_aI_sl_data4_type_c fcs_front_xq_aI_sl_data4_size_c fcs_front_xq_aI_sl_data4_type_mpi fcs_front_xq_aI_sl_data4_size_mpi fcs_front_xq_aI_sl_data4_memcpy fcs_front_xq_aI_sl_data4_weight fcs_front_xq_aI_sl_data4_flex */

#ifdef fcs_front_xq_aI_SL_DATA4

 #define fcs_front_xq_aI_sl_data4_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data4_size_c) * sizeof(fcs_front_xq_aI_sl_data4_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data4_copy
  #if fcs_front_xq_aI_sl_data4_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data4_memcpy)
   #if fcs_front_xq_aI_sl_data4_size_c == 1
    #define fcs_front_xq_aI_sl_data4_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data4_size_c == 2
    #define fcs_front_xq_aI_sl_data4_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data4_size_c == 3
    #define fcs_front_xq_aI_sl_data4_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data4_size_c == 4
    #define fcs_front_xq_aI_sl_data4_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data4_size_c == 5
    #define fcs_front_xq_aI_sl_data4_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data4_size_c == 6
    #define fcs_front_xq_aI_sl_data4_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data4_size_c == 7
    #define fcs_front_xq_aI_sl_data4_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data4_size_c == 8
    #define fcs_front_xq_aI_sl_data4_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data4_size_c == 9
    #define fcs_front_xq_aI_sl_data4_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data4_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data4_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data4_ncopy
  #define fcs_front_xq_aI_sl_data4_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data4_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data4_nmove
  #define fcs_front_xq_aI_sl_data4_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data4_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data4_type_c                              fcs_front_xq_aI_sl_data4_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data4_size_c                              (fcs_front_xq_aI_sl_data4_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data4_type_mpi                            (fcs_front_xq_aI_sl_data4_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data4_size_mpi                            (fcs_front_xq_aI_sl_data4_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data4_idx                                 4  /* sl_macro */

 #define fcs_front_xq_aI_data4_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data4_byte                                (fcs_front_xq_aI_sl_data4_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data4_ptr(e)                              (e)->data4  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data4_flex
 # define fcs_front_xq_aI_data4_byte_flex                          (fcs_front_xq_aI_sl_data4_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data4_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data4_weight
 # define fcs_front_xq_aI_data4_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data4_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data4_assign(src, dst)                    ((dst)->data4 = (src)->data4)  /* sl_macro */
 #define fcs_front_xq_aI_data4_assign_at(src, sat, dst)            ((dst)->data4 = &(src)->data4[(sat) * fcs_front_xq_aI_data4_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data4_null(e)                             ((e)->data4 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data4_inc(e)                              ((e)->data4 += fcs_front_xq_aI_data4_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data4_dec(e)                              ((e)->data4 -= fcs_front_xq_aI_data4_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data4_add(e, n)                           ((e)->data4 += (n) * fcs_front_xq_aI_data4_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data4_sub(e, n)                           ((e)->data4 -= (n) * fcs_front_xq_aI_data4_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data4_copy(src, dst)                      fcs_front_xq_aI_sl_data4_copy((src)->data4, (dst)->data4)  /* sl_macro */
 #define fcs_front_xq_aI_data4_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data4_ncopy((src)->data4, (dst)->data4, n)  /* sl_macro */
 #define fcs_front_xq_aI_data4_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data4_nmove((src)->data4, (dst)->data4, n)  /* sl_macro */

 #define fcs_front_xq_aI_data4_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data4_copy(&(src)->data4[(sat) * fcs_front_xq_aI_data4_size_c], &(dst)->data4[(dat) * fcs_front_xq_aI_data4_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data4_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data4_ncopy(&(src)->data4[(sat) * fcs_front_xq_aI_data4_size_c], &(dst)->data4[(dat) * fcs_front_xq_aI_data4_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data4_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data4_nmove(&(src)->data4[(sat) * fcs_front_xq_aI_data4_size_c], &(dst)->data4[(dat) * fcs_front_xq_aI_data4_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data4_xchange(e0, e1, t)                  (fcs_front_xq_aI_data4_copy(e0, t), fcs_front_xq_aI_data4_copy(e1, e0), fcs_front_xq_aI_data4_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data4_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data4_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data4_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data4_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data4_assign(src, dst)                 , fcs_front_xq_aI_data4_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_assign_at(src, sat, dst)         , fcs_front_xq_aI_data4_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_null(e)                          , fcs_front_xq_aI_data4_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_inc(e)                           , fcs_front_xq_aI_data4_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_dec(e)                           , fcs_front_xq_aI_data4_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_add(e, n)                        , fcs_front_xq_aI_data4_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_sub(e, n)                        , fcs_front_xq_aI_data4_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_copy(src, dst)                   , fcs_front_xq_aI_data4_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_ncopy(src, dst, n)               , fcs_front_xq_aI_data4_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_nmove(src, dst, n)               , fcs_front_xq_aI_data4_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data4_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data4_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data4_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_xchange(e0, e1, t)               , fcs_front_xq_aI_data4_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data4_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data4_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data4_assign(src, dst)                 fcs_front_xq_aI_data4_assign(src, dst)
 #define fcs_front_xq_aI_cc_data4_assign_at(src, sat, dst)         fcs_front_xq_aI_data4_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data4_null(e)                          fcs_front_xq_aI_data4_null(e)
 #define fcs_front_xq_aI_cc_data4_inc(e)                           fcs_front_xq_aI_data4_inc(e)
 #define fcs_front_xq_aI_cc_data4_dec(e)                           fcs_front_xq_aI_data4_dec(e)
 #define fcs_front_xq_aI_cc_data4_add(e, n)                        fcs_front_xq_aI_data4_add(e, n)
 #define fcs_front_xq_aI_cc_data4_sub(e, n)                        fcs_front_xq_aI_data4_sub(e, n)
 #define fcs_front_xq_aI_cc_data4_copy(src, dst)                   fcs_front_xq_aI_data4_copy(src, dst)
 #define fcs_front_xq_aI_cc_data4_ncopy(src, dst, n)               fcs_front_xq_aI_data4_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data4_nmove(src, dst, n)               fcs_front_xq_aI_data4_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data4_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data4_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data4_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data4_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data4_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data4_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data4_xchange(e0, e1, t)               fcs_front_xq_aI_data4_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data4_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data4_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA4 */

 #define fcs_front_xq_aI_data4_n                                   0
 #define fcs_front_xq_aI_data4_byte                                0
/* #define fcs_front_xq_aI_data4_ptr(e)*/

 #define fcs_front_xq_aI_data4_byte_flex                           0
 #define fcs_front_xq_aI_data4_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data4_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data4_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data4_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data4_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data4_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data4_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data4_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data4_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data4_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data4_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data4_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data4_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data4_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data4_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data4_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data4_assign(src, dst)
 #define fcs_front_xq_aI_cc_data4_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data4_null(e)
 #define fcs_front_xq_aI_cc_data4_inc(e)
 #define fcs_front_xq_aI_cc_data4_dec(e)
 #define fcs_front_xq_aI_cc_data4_add(e, n)
 #define fcs_front_xq_aI_cc_data4_sub(e, n)
 #define fcs_front_xq_aI_cc_data4_copy(src, dst)
 #define fcs_front_xq_aI_cc_data4_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data4_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data4_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data4_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data4_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data4_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data4_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA4 */

#define fcs_front_xq_aI_data4_cm                                   SLCM_DATA4  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA5 fcs_front_xq_aI_SL_DATA5_IGNORE fcs_front_xq_aI_sl_data5_type_c fcs_front_xq_aI_sl_data5_size_c fcs_front_xq_aI_sl_data5_type_mpi fcs_front_xq_aI_sl_data5_size_mpi fcs_front_xq_aI_sl_data5_memcpy fcs_front_xq_aI_sl_data5_weight fcs_front_xq_aI_sl_data5_flex */

#ifdef fcs_front_xq_aI_SL_DATA5

 #define fcs_front_xq_aI_sl_data5_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data5_size_c) * sizeof(fcs_front_xq_aI_sl_data5_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data5_copy
  #if fcs_front_xq_aI_sl_data5_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data5_memcpy)
   #if fcs_front_xq_aI_sl_data5_size_c == 1
    #define fcs_front_xq_aI_sl_data5_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data5_size_c == 2
    #define fcs_front_xq_aI_sl_data5_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data5_size_c == 3
    #define fcs_front_xq_aI_sl_data5_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data5_size_c == 4
    #define fcs_front_xq_aI_sl_data5_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data5_size_c == 5
    #define fcs_front_xq_aI_sl_data5_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data5_size_c == 6
    #define fcs_front_xq_aI_sl_data5_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data5_size_c == 7
    #define fcs_front_xq_aI_sl_data5_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data5_size_c == 8
    #define fcs_front_xq_aI_sl_data5_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data5_size_c == 9
    #define fcs_front_xq_aI_sl_data5_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data5_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data5_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data5_ncopy
  #define fcs_front_xq_aI_sl_data5_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data5_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data5_nmove
  #define fcs_front_xq_aI_sl_data5_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data5_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data5_type_c                              fcs_front_xq_aI_sl_data5_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data5_size_c                              (fcs_front_xq_aI_sl_data5_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data5_type_mpi                            (fcs_front_xq_aI_sl_data5_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data5_size_mpi                            (fcs_front_xq_aI_sl_data5_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data5_idx                                 5  /* sl_macro */

 #define fcs_front_xq_aI_data5_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data5_byte                                (fcs_front_xq_aI_sl_data5_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data5_ptr(e)                              (e)->data5  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data5_flex
 # define fcs_front_xq_aI_data5_byte_flex                          (fcs_front_xq_aI_sl_data5_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data5_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data5_weight
 # define fcs_front_xq_aI_data5_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data5_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data5_assign(src, dst)                    ((dst)->data5 = (src)->data5)  /* sl_macro */
 #define fcs_front_xq_aI_data5_assign_at(src, sat, dst)            ((dst)->data5 = &(src)->data5[(sat) * fcs_front_xq_aI_data5_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data5_null(e)                             ((e)->data5 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data5_inc(e)                              ((e)->data5 += fcs_front_xq_aI_data5_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data5_dec(e)                              ((e)->data5 -= fcs_front_xq_aI_data5_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data5_add(e, n)                           ((e)->data5 += (n) * fcs_front_xq_aI_data5_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data5_sub(e, n)                           ((e)->data5 -= (n) * fcs_front_xq_aI_data5_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data5_copy(src, dst)                      fcs_front_xq_aI_sl_data5_copy((src)->data5, (dst)->data5)  /* sl_macro */
 #define fcs_front_xq_aI_data5_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data5_ncopy((src)->data5, (dst)->data5, n)  /* sl_macro */
 #define fcs_front_xq_aI_data5_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data5_nmove((src)->data5, (dst)->data5, n)  /* sl_macro */

 #define fcs_front_xq_aI_data5_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data5_copy(&(src)->data5[(sat) * fcs_front_xq_aI_data5_size_c], &(dst)->data5[(dat) * fcs_front_xq_aI_data5_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data5_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data5_ncopy(&(src)->data5[(sat) * fcs_front_xq_aI_data5_size_c], &(dst)->data5[(dat) * fcs_front_xq_aI_data5_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data5_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data5_nmove(&(src)->data5[(sat) * fcs_front_xq_aI_data5_size_c], &(dst)->data5[(dat) * fcs_front_xq_aI_data5_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data5_xchange(e0, e1, t)                  (fcs_front_xq_aI_data5_copy(e0, t), fcs_front_xq_aI_data5_copy(e1, e0), fcs_front_xq_aI_data5_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data5_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data5_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data5_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data5_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data5_assign(src, dst)                 , fcs_front_xq_aI_data5_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_assign_at(src, sat, dst)         , fcs_front_xq_aI_data5_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_null(e)                          , fcs_front_xq_aI_data5_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_inc(e)                           , fcs_front_xq_aI_data5_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_dec(e)                           , fcs_front_xq_aI_data5_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_add(e, n)                        , fcs_front_xq_aI_data5_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_sub(e, n)                        , fcs_front_xq_aI_data5_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_copy(src, dst)                   , fcs_front_xq_aI_data5_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_ncopy(src, dst, n)               , fcs_front_xq_aI_data5_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_nmove(src, dst, n)               , fcs_front_xq_aI_data5_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data5_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data5_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data5_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_xchange(e0, e1, t)               , fcs_front_xq_aI_data5_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data5_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data5_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data5_assign(src, dst)                 fcs_front_xq_aI_data5_assign(src, dst)
 #define fcs_front_xq_aI_cc_data5_assign_at(src, sat, dst)         fcs_front_xq_aI_data5_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data5_null(e)                          fcs_front_xq_aI_data5_null(e)
 #define fcs_front_xq_aI_cc_data5_inc(e)                           fcs_front_xq_aI_data5_inc(e)
 #define fcs_front_xq_aI_cc_data5_dec(e)                           fcs_front_xq_aI_data5_dec(e)
 #define fcs_front_xq_aI_cc_data5_add(e, n)                        fcs_front_xq_aI_data5_add(e, n)
 #define fcs_front_xq_aI_cc_data5_sub(e, n)                        fcs_front_xq_aI_data5_sub(e, n)
 #define fcs_front_xq_aI_cc_data5_copy(src, dst)                   fcs_front_xq_aI_data5_copy(src, dst)
 #define fcs_front_xq_aI_cc_data5_ncopy(src, dst, n)               fcs_front_xq_aI_data5_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data5_nmove(src, dst, n)               fcs_front_xq_aI_data5_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data5_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data5_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data5_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data5_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data5_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data5_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data5_xchange(e0, e1, t)               fcs_front_xq_aI_data5_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data5_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data5_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA5 */

 #define fcs_front_xq_aI_data5_n                                   0
 #define fcs_front_xq_aI_data5_byte                                0
/* #define fcs_front_xq_aI_data5_ptr(e)*/

 #define fcs_front_xq_aI_data5_byte_flex                           0
 #define fcs_front_xq_aI_data5_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data5_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data5_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data5_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data5_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data5_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data5_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data5_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data5_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data5_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data5_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data5_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data5_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data5_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data5_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data5_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data5_assign(src, dst)
 #define fcs_front_xq_aI_cc_data5_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data5_null(e)
 #define fcs_front_xq_aI_cc_data5_inc(e)
 #define fcs_front_xq_aI_cc_data5_dec(e)
 #define fcs_front_xq_aI_cc_data5_add(e, n)
 #define fcs_front_xq_aI_cc_data5_sub(e, n)
 #define fcs_front_xq_aI_cc_data5_copy(src, dst)
 #define fcs_front_xq_aI_cc_data5_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data5_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data5_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data5_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data5_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data5_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data5_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA5 */

#define fcs_front_xq_aI_data5_cm                                   SLCM_DATA5  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA6 fcs_front_xq_aI_SL_DATA6_IGNORE fcs_front_xq_aI_sl_data6_type_c fcs_front_xq_aI_sl_data6_size_c fcs_front_xq_aI_sl_data6_type_mpi fcs_front_xq_aI_sl_data6_size_mpi fcs_front_xq_aI_sl_data6_memcpy fcs_front_xq_aI_sl_data6_weight fcs_front_xq_aI_sl_data6_flex */

#ifdef fcs_front_xq_aI_SL_DATA6

 #define fcs_front_xq_aI_sl_data6_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data6_size_c) * sizeof(fcs_front_xq_aI_sl_data6_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data6_copy
  #if fcs_front_xq_aI_sl_data6_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data6_memcpy)
   #if fcs_front_xq_aI_sl_data6_size_c == 1
    #define fcs_front_xq_aI_sl_data6_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data6_size_c == 2
    #define fcs_front_xq_aI_sl_data6_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data6_size_c == 3
    #define fcs_front_xq_aI_sl_data6_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data6_size_c == 4
    #define fcs_front_xq_aI_sl_data6_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data6_size_c == 5
    #define fcs_front_xq_aI_sl_data6_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data6_size_c == 6
    #define fcs_front_xq_aI_sl_data6_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data6_size_c == 7
    #define fcs_front_xq_aI_sl_data6_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data6_size_c == 8
    #define fcs_front_xq_aI_sl_data6_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data6_size_c == 9
    #define fcs_front_xq_aI_sl_data6_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data6_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data6_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data6_ncopy
  #define fcs_front_xq_aI_sl_data6_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data6_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data6_nmove
  #define fcs_front_xq_aI_sl_data6_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data6_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data6_type_c                              fcs_front_xq_aI_sl_data6_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data6_size_c                              (fcs_front_xq_aI_sl_data6_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data6_type_mpi                            (fcs_front_xq_aI_sl_data6_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data6_size_mpi                            (fcs_front_xq_aI_sl_data6_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data6_idx                                 6  /* sl_macro */

 #define fcs_front_xq_aI_data6_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data6_byte                                (fcs_front_xq_aI_sl_data6_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data6_ptr(e)                              (e)->data6  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data6_flex
 # define fcs_front_xq_aI_data6_byte_flex                          (fcs_front_xq_aI_sl_data6_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data6_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data6_weight
 # define fcs_front_xq_aI_data6_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data6_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data6_assign(src, dst)                    ((dst)->data6 = (src)->data6)  /* sl_macro */
 #define fcs_front_xq_aI_data6_assign_at(src, sat, dst)            ((dst)->data6 = &(src)->data6[(sat) * fcs_front_xq_aI_data6_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data6_null(e)                             ((e)->data6 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data6_inc(e)                              ((e)->data6 += fcs_front_xq_aI_data6_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data6_dec(e)                              ((e)->data6 -= fcs_front_xq_aI_data6_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data6_add(e, n)                           ((e)->data6 += (n) * fcs_front_xq_aI_data6_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data6_sub(e, n)                           ((e)->data6 -= (n) * fcs_front_xq_aI_data6_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data6_copy(src, dst)                      fcs_front_xq_aI_sl_data6_copy((src)->data6, (dst)->data6)  /* sl_macro */
 #define fcs_front_xq_aI_data6_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data6_ncopy((src)->data6, (dst)->data6, n)  /* sl_macro */
 #define fcs_front_xq_aI_data6_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data6_nmove((src)->data6, (dst)->data6, n)  /* sl_macro */

 #define fcs_front_xq_aI_data6_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data6_copy(&(src)->data6[(sat) * fcs_front_xq_aI_data6_size_c], &(dst)->data6[(dat) * fcs_front_xq_aI_data6_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data6_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data6_ncopy(&(src)->data6[(sat) * fcs_front_xq_aI_data6_size_c], &(dst)->data6[(dat) * fcs_front_xq_aI_data6_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data6_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data6_nmove(&(src)->data6[(sat) * fcs_front_xq_aI_data6_size_c], &(dst)->data6[(dat) * fcs_front_xq_aI_data6_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data6_xchange(e0, e1, t)                  (fcs_front_xq_aI_data6_copy(e0, t), fcs_front_xq_aI_data6_copy(e1, e0), fcs_front_xq_aI_data6_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data6_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data6_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data6_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data6_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data6_assign(src, dst)                 , fcs_front_xq_aI_data6_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_assign_at(src, sat, dst)         , fcs_front_xq_aI_data6_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_null(e)                          , fcs_front_xq_aI_data6_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_inc(e)                           , fcs_front_xq_aI_data6_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_dec(e)                           , fcs_front_xq_aI_data6_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_add(e, n)                        , fcs_front_xq_aI_data6_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_sub(e, n)                        , fcs_front_xq_aI_data6_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_copy(src, dst)                   , fcs_front_xq_aI_data6_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_ncopy(src, dst, n)               , fcs_front_xq_aI_data6_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_nmove(src, dst, n)               , fcs_front_xq_aI_data6_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data6_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data6_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data6_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_xchange(e0, e1, t)               , fcs_front_xq_aI_data6_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data6_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data6_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data6_assign(src, dst)                 fcs_front_xq_aI_data6_assign(src, dst)
 #define fcs_front_xq_aI_cc_data6_assign_at(src, sat, dst)         fcs_front_xq_aI_data6_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data6_null(e)                          fcs_front_xq_aI_data6_null(e)
 #define fcs_front_xq_aI_cc_data6_inc(e)                           fcs_front_xq_aI_data6_inc(e)
 #define fcs_front_xq_aI_cc_data6_dec(e)                           fcs_front_xq_aI_data6_dec(e)
 #define fcs_front_xq_aI_cc_data6_add(e, n)                        fcs_front_xq_aI_data6_add(e, n)
 #define fcs_front_xq_aI_cc_data6_sub(e, n)                        fcs_front_xq_aI_data6_sub(e, n)
 #define fcs_front_xq_aI_cc_data6_copy(src, dst)                   fcs_front_xq_aI_data6_copy(src, dst)
 #define fcs_front_xq_aI_cc_data6_ncopy(src, dst, n)               fcs_front_xq_aI_data6_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data6_nmove(src, dst, n)               fcs_front_xq_aI_data6_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data6_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data6_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data6_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data6_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data6_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data6_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data6_xchange(e0, e1, t)               fcs_front_xq_aI_data6_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data6_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data6_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA6 */

 #define fcs_front_xq_aI_data6_n                                   0
 #define fcs_front_xq_aI_data6_byte                                0
/* #define fcs_front_xq_aI_data6_ptr(e)*/

 #define fcs_front_xq_aI_data6_byte_flex                           0
 #define fcs_front_xq_aI_data6_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data6_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data6_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data6_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data6_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data6_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data6_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data6_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data6_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data6_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data6_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data6_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data6_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data6_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data6_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data6_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data6_assign(src, dst)
 #define fcs_front_xq_aI_cc_data6_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data6_null(e)
 #define fcs_front_xq_aI_cc_data6_inc(e)
 #define fcs_front_xq_aI_cc_data6_dec(e)
 #define fcs_front_xq_aI_cc_data6_add(e, n)
 #define fcs_front_xq_aI_cc_data6_sub(e, n)
 #define fcs_front_xq_aI_cc_data6_copy(src, dst)
 #define fcs_front_xq_aI_cc_data6_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data6_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data6_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data6_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data6_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data6_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data6_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA6 */

#define fcs_front_xq_aI_data6_cm                                   SLCM_DATA6  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA7 fcs_front_xq_aI_SL_DATA7_IGNORE fcs_front_xq_aI_sl_data7_type_c fcs_front_xq_aI_sl_data7_size_c fcs_front_xq_aI_sl_data7_type_mpi fcs_front_xq_aI_sl_data7_size_mpi fcs_front_xq_aI_sl_data7_memcpy fcs_front_xq_aI_sl_data7_weight fcs_front_xq_aI_sl_data7_flex */

#ifdef fcs_front_xq_aI_SL_DATA7

 #define fcs_front_xq_aI_sl_data7_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data7_size_c) * sizeof(fcs_front_xq_aI_sl_data7_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data7_copy
  #if fcs_front_xq_aI_sl_data7_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data7_memcpy)
   #if fcs_front_xq_aI_sl_data7_size_c == 1
    #define fcs_front_xq_aI_sl_data7_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data7_size_c == 2
    #define fcs_front_xq_aI_sl_data7_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data7_size_c == 3
    #define fcs_front_xq_aI_sl_data7_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data7_size_c == 4
    #define fcs_front_xq_aI_sl_data7_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data7_size_c == 5
    #define fcs_front_xq_aI_sl_data7_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data7_size_c == 6
    #define fcs_front_xq_aI_sl_data7_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data7_size_c == 7
    #define fcs_front_xq_aI_sl_data7_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data7_size_c == 8
    #define fcs_front_xq_aI_sl_data7_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data7_size_c == 9
    #define fcs_front_xq_aI_sl_data7_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data7_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data7_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data7_ncopy
  #define fcs_front_xq_aI_sl_data7_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data7_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data7_nmove
  #define fcs_front_xq_aI_sl_data7_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data7_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data7_type_c                              fcs_front_xq_aI_sl_data7_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data7_size_c                              (fcs_front_xq_aI_sl_data7_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data7_type_mpi                            (fcs_front_xq_aI_sl_data7_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data7_size_mpi                            (fcs_front_xq_aI_sl_data7_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data7_idx                                 7  /* sl_macro */

 #define fcs_front_xq_aI_data7_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data7_byte                                (fcs_front_xq_aI_sl_data7_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data7_ptr(e)                              (e)->data7  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data7_flex
 # define fcs_front_xq_aI_data7_byte_flex                          (fcs_front_xq_aI_sl_data7_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data7_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data7_weight
 # define fcs_front_xq_aI_data7_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data7_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data7_assign(src, dst)                    ((dst)->data7 = (src)->data7)  /* sl_macro */
 #define fcs_front_xq_aI_data7_assign_at(src, sat, dst)            ((dst)->data7 = &(src)->data7[(sat) * fcs_front_xq_aI_data7_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data7_null(e)                             ((e)->data7 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data7_inc(e)                              ((e)->data7 += fcs_front_xq_aI_data7_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data7_dec(e)                              ((e)->data7 -= fcs_front_xq_aI_data7_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data7_add(e, n)                           ((e)->data7 += (n) * fcs_front_xq_aI_data7_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data7_sub(e, n)                           ((e)->data7 -= (n) * fcs_front_xq_aI_data7_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data7_copy(src, dst)                      fcs_front_xq_aI_sl_data7_copy((src)->data7, (dst)->data7)  /* sl_macro */
 #define fcs_front_xq_aI_data7_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data7_ncopy((src)->data7, (dst)->data7, n)  /* sl_macro */
 #define fcs_front_xq_aI_data7_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data7_nmove((src)->data7, (dst)->data7, n)  /* sl_macro */

 #define fcs_front_xq_aI_data7_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data7_copy(&(src)->data7[(sat) * fcs_front_xq_aI_data7_size_c], &(dst)->data7[(dat) * fcs_front_xq_aI_data7_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data7_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data7_ncopy(&(src)->data7[(sat) * fcs_front_xq_aI_data7_size_c], &(dst)->data7[(dat) * fcs_front_xq_aI_data7_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data7_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data7_nmove(&(src)->data7[(sat) * fcs_front_xq_aI_data7_size_c], &(dst)->data7[(dat) * fcs_front_xq_aI_data7_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data7_xchange(e0, e1, t)                  (fcs_front_xq_aI_data7_copy(e0, t), fcs_front_xq_aI_data7_copy(e1, e0), fcs_front_xq_aI_data7_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data7_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data7_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data7_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data7_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data7_assign(src, dst)                 , fcs_front_xq_aI_data7_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_assign_at(src, sat, dst)         , fcs_front_xq_aI_data7_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_null(e)                          , fcs_front_xq_aI_data7_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_inc(e)                           , fcs_front_xq_aI_data7_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_dec(e)                           , fcs_front_xq_aI_data7_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_add(e, n)                        , fcs_front_xq_aI_data7_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_sub(e, n)                        , fcs_front_xq_aI_data7_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_copy(src, dst)                   , fcs_front_xq_aI_data7_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_ncopy(src, dst, n)               , fcs_front_xq_aI_data7_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_nmove(src, dst, n)               , fcs_front_xq_aI_data7_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data7_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data7_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data7_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_xchange(e0, e1, t)               , fcs_front_xq_aI_data7_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data7_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data7_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data7_assign(src, dst)                 fcs_front_xq_aI_data7_assign(src, dst)
 #define fcs_front_xq_aI_cc_data7_assign_at(src, sat, dst)         fcs_front_xq_aI_data7_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data7_null(e)                          fcs_front_xq_aI_data7_null(e)
 #define fcs_front_xq_aI_cc_data7_inc(e)                           fcs_front_xq_aI_data7_inc(e)
 #define fcs_front_xq_aI_cc_data7_dec(e)                           fcs_front_xq_aI_data7_dec(e)
 #define fcs_front_xq_aI_cc_data7_add(e, n)                        fcs_front_xq_aI_data7_add(e, n)
 #define fcs_front_xq_aI_cc_data7_sub(e, n)                        fcs_front_xq_aI_data7_sub(e, n)
 #define fcs_front_xq_aI_cc_data7_copy(src, dst)                   fcs_front_xq_aI_data7_copy(src, dst)
 #define fcs_front_xq_aI_cc_data7_ncopy(src, dst, n)               fcs_front_xq_aI_data7_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data7_nmove(src, dst, n)               fcs_front_xq_aI_data7_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data7_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data7_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data7_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data7_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data7_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data7_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data7_xchange(e0, e1, t)               fcs_front_xq_aI_data7_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data7_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data7_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA7 */

 #define fcs_front_xq_aI_data7_n                                   0
 #define fcs_front_xq_aI_data7_byte                                0
/* #define fcs_front_xq_aI_data7_ptr(e)*/

 #define fcs_front_xq_aI_data7_byte_flex                           0
 #define fcs_front_xq_aI_data7_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data7_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data7_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data7_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data7_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data7_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data7_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data7_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data7_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data7_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data7_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data7_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data7_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data7_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data7_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data7_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data7_assign(src, dst)
 #define fcs_front_xq_aI_cc_data7_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data7_null(e)
 #define fcs_front_xq_aI_cc_data7_inc(e)
 #define fcs_front_xq_aI_cc_data7_dec(e)
 #define fcs_front_xq_aI_cc_data7_add(e, n)
 #define fcs_front_xq_aI_cc_data7_sub(e, n)
 #define fcs_front_xq_aI_cc_data7_copy(src, dst)
 #define fcs_front_xq_aI_cc_data7_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data7_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data7_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data7_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data7_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data7_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data7_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA7 */

#define fcs_front_xq_aI_data7_cm                                   SLCM_DATA7  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA8 fcs_front_xq_aI_SL_DATA8_IGNORE fcs_front_xq_aI_sl_data8_type_c fcs_front_xq_aI_sl_data8_size_c fcs_front_xq_aI_sl_data8_type_mpi fcs_front_xq_aI_sl_data8_size_mpi fcs_front_xq_aI_sl_data8_memcpy fcs_front_xq_aI_sl_data8_weight fcs_front_xq_aI_sl_data8_flex */

#ifdef fcs_front_xq_aI_SL_DATA8

 #define fcs_front_xq_aI_sl_data8_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data8_size_c) * sizeof(fcs_front_xq_aI_sl_data8_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data8_copy
  #if fcs_front_xq_aI_sl_data8_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data8_memcpy)
   #if fcs_front_xq_aI_sl_data8_size_c == 1
    #define fcs_front_xq_aI_sl_data8_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data8_size_c == 2
    #define fcs_front_xq_aI_sl_data8_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data8_size_c == 3
    #define fcs_front_xq_aI_sl_data8_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data8_size_c == 4
    #define fcs_front_xq_aI_sl_data8_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data8_size_c == 5
    #define fcs_front_xq_aI_sl_data8_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data8_size_c == 6
    #define fcs_front_xq_aI_sl_data8_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data8_size_c == 7
    #define fcs_front_xq_aI_sl_data8_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data8_size_c == 8
    #define fcs_front_xq_aI_sl_data8_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data8_size_c == 9
    #define fcs_front_xq_aI_sl_data8_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data8_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data8_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data8_ncopy
  #define fcs_front_xq_aI_sl_data8_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data8_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data8_nmove
  #define fcs_front_xq_aI_sl_data8_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data8_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data8_type_c                              fcs_front_xq_aI_sl_data8_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data8_size_c                              (fcs_front_xq_aI_sl_data8_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data8_type_mpi                            (fcs_front_xq_aI_sl_data8_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data8_size_mpi                            (fcs_front_xq_aI_sl_data8_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data8_idx                                 8  /* sl_macro */

 #define fcs_front_xq_aI_data8_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data8_byte                                (fcs_front_xq_aI_sl_data8_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data8_ptr(e)                              (e)->data8  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data8_flex
 # define fcs_front_xq_aI_data8_byte_flex                          (fcs_front_xq_aI_sl_data8_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data8_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data8_weight
 # define fcs_front_xq_aI_data8_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data8_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data8_assign(src, dst)                    ((dst)->data8 = (src)->data8)  /* sl_macro */
 #define fcs_front_xq_aI_data8_assign_at(src, sat, dst)            ((dst)->data8 = &(src)->data8[(sat) * fcs_front_xq_aI_data8_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data8_null(e)                             ((e)->data8 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data8_inc(e)                              ((e)->data8 += fcs_front_xq_aI_data8_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data8_dec(e)                              ((e)->data8 -= fcs_front_xq_aI_data8_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data8_add(e, n)                           ((e)->data8 += (n) * fcs_front_xq_aI_data8_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data8_sub(e, n)                           ((e)->data8 -= (n) * fcs_front_xq_aI_data8_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data8_copy(src, dst)                      fcs_front_xq_aI_sl_data8_copy((src)->data8, (dst)->data8)  /* sl_macro */
 #define fcs_front_xq_aI_data8_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data8_ncopy((src)->data8, (dst)->data8, n)  /* sl_macro */
 #define fcs_front_xq_aI_data8_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data8_nmove((src)->data8, (dst)->data8, n)  /* sl_macro */

 #define fcs_front_xq_aI_data8_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data8_copy(&(src)->data8[(sat) * fcs_front_xq_aI_data8_size_c], &(dst)->data8[(dat) * fcs_front_xq_aI_data8_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data8_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data8_ncopy(&(src)->data8[(sat) * fcs_front_xq_aI_data8_size_c], &(dst)->data8[(dat) * fcs_front_xq_aI_data8_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data8_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data8_nmove(&(src)->data8[(sat) * fcs_front_xq_aI_data8_size_c], &(dst)->data8[(dat) * fcs_front_xq_aI_data8_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data8_xchange(e0, e1, t)                  (fcs_front_xq_aI_data8_copy(e0, t), fcs_front_xq_aI_data8_copy(e1, e0), fcs_front_xq_aI_data8_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data8_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data8_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data8_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data8_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data8_assign(src, dst)                 , fcs_front_xq_aI_data8_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_assign_at(src, sat, dst)         , fcs_front_xq_aI_data8_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_null(e)                          , fcs_front_xq_aI_data8_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_inc(e)                           , fcs_front_xq_aI_data8_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_dec(e)                           , fcs_front_xq_aI_data8_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_add(e, n)                        , fcs_front_xq_aI_data8_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_sub(e, n)                        , fcs_front_xq_aI_data8_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_copy(src, dst)                   , fcs_front_xq_aI_data8_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_ncopy(src, dst, n)               , fcs_front_xq_aI_data8_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_nmove(src, dst, n)               , fcs_front_xq_aI_data8_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data8_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data8_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data8_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_xchange(e0, e1, t)               , fcs_front_xq_aI_data8_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data8_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data8_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data8_assign(src, dst)                 fcs_front_xq_aI_data8_assign(src, dst)
 #define fcs_front_xq_aI_cc_data8_assign_at(src, sat, dst)         fcs_front_xq_aI_data8_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data8_null(e)                          fcs_front_xq_aI_data8_null(e)
 #define fcs_front_xq_aI_cc_data8_inc(e)                           fcs_front_xq_aI_data8_inc(e)
 #define fcs_front_xq_aI_cc_data8_dec(e)                           fcs_front_xq_aI_data8_dec(e)
 #define fcs_front_xq_aI_cc_data8_add(e, n)                        fcs_front_xq_aI_data8_add(e, n)
 #define fcs_front_xq_aI_cc_data8_sub(e, n)                        fcs_front_xq_aI_data8_sub(e, n)
 #define fcs_front_xq_aI_cc_data8_copy(src, dst)                   fcs_front_xq_aI_data8_copy(src, dst)
 #define fcs_front_xq_aI_cc_data8_ncopy(src, dst, n)               fcs_front_xq_aI_data8_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data8_nmove(src, dst, n)               fcs_front_xq_aI_data8_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data8_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data8_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data8_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data8_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data8_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data8_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data8_xchange(e0, e1, t)               fcs_front_xq_aI_data8_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data8_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data8_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA8 */

 #define fcs_front_xq_aI_data8_n                                   0
 #define fcs_front_xq_aI_data8_byte                                0
/* #define fcs_front_xq_aI_data8_ptr(e)*/

 #define fcs_front_xq_aI_data8_byte_flex                           0
 #define fcs_front_xq_aI_data8_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data8_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data8_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data8_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data8_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data8_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data8_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data8_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data8_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data8_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data8_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data8_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data8_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data8_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data8_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data8_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data8_assign(src, dst)
 #define fcs_front_xq_aI_cc_data8_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data8_null(e)
 #define fcs_front_xq_aI_cc_data8_inc(e)
 #define fcs_front_xq_aI_cc_data8_dec(e)
 #define fcs_front_xq_aI_cc_data8_add(e, n)
 #define fcs_front_xq_aI_cc_data8_sub(e, n)
 #define fcs_front_xq_aI_cc_data8_copy(src, dst)
 #define fcs_front_xq_aI_cc_data8_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data8_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data8_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data8_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data8_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data8_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data8_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA8 */

#define fcs_front_xq_aI_data8_cm                                   SLCM_DATA8  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA9 fcs_front_xq_aI_SL_DATA9_IGNORE fcs_front_xq_aI_sl_data9_type_c fcs_front_xq_aI_sl_data9_size_c fcs_front_xq_aI_sl_data9_type_mpi fcs_front_xq_aI_sl_data9_size_mpi fcs_front_xq_aI_sl_data9_memcpy fcs_front_xq_aI_sl_data9_weight fcs_front_xq_aI_sl_data9_flex */

#ifdef fcs_front_xq_aI_SL_DATA9

 #define fcs_front_xq_aI_sl_data9_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data9_size_c) * sizeof(fcs_front_xq_aI_sl_data9_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data9_copy
  #if fcs_front_xq_aI_sl_data9_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data9_memcpy)
   #if fcs_front_xq_aI_sl_data9_size_c == 1
    #define fcs_front_xq_aI_sl_data9_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data9_size_c == 2
    #define fcs_front_xq_aI_sl_data9_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data9_size_c == 3
    #define fcs_front_xq_aI_sl_data9_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data9_size_c == 4
    #define fcs_front_xq_aI_sl_data9_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data9_size_c == 5
    #define fcs_front_xq_aI_sl_data9_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data9_size_c == 6
    #define fcs_front_xq_aI_sl_data9_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data9_size_c == 7
    #define fcs_front_xq_aI_sl_data9_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data9_size_c == 8
    #define fcs_front_xq_aI_sl_data9_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data9_size_c == 9
    #define fcs_front_xq_aI_sl_data9_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data9_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data9_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data9_ncopy
  #define fcs_front_xq_aI_sl_data9_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data9_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data9_nmove
  #define fcs_front_xq_aI_sl_data9_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data9_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data9_type_c                              fcs_front_xq_aI_sl_data9_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data9_size_c                              (fcs_front_xq_aI_sl_data9_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data9_type_mpi                            (fcs_front_xq_aI_sl_data9_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data9_size_mpi                            (fcs_front_xq_aI_sl_data9_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data9_idx                                 9  /* sl_macro */

 #define fcs_front_xq_aI_data9_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data9_byte                                (fcs_front_xq_aI_sl_data9_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data9_ptr(e)                              (e)->data9  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data9_flex
 # define fcs_front_xq_aI_data9_byte_flex                          (fcs_front_xq_aI_sl_data9_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data9_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data9_weight
 # define fcs_front_xq_aI_data9_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data9_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data9_assign(src, dst)                    ((dst)->data9 = (src)->data9)  /* sl_macro */
 #define fcs_front_xq_aI_data9_assign_at(src, sat, dst)            ((dst)->data9 = &(src)->data9[(sat) * fcs_front_xq_aI_data9_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data9_null(e)                             ((e)->data9 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data9_inc(e)                              ((e)->data9 += fcs_front_xq_aI_data9_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data9_dec(e)                              ((e)->data9 -= fcs_front_xq_aI_data9_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data9_add(e, n)                           ((e)->data9 += (n) * fcs_front_xq_aI_data9_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data9_sub(e, n)                           ((e)->data9 -= (n) * fcs_front_xq_aI_data9_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data9_copy(src, dst)                      fcs_front_xq_aI_sl_data9_copy((src)->data9, (dst)->data9)  /* sl_macro */
 #define fcs_front_xq_aI_data9_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data9_ncopy((src)->data9, (dst)->data9, n)  /* sl_macro */
 #define fcs_front_xq_aI_data9_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data9_nmove((src)->data9, (dst)->data9, n)  /* sl_macro */

 #define fcs_front_xq_aI_data9_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data9_copy(&(src)->data9[(sat) * fcs_front_xq_aI_data9_size_c], &(dst)->data9[(dat) * fcs_front_xq_aI_data9_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data9_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data9_ncopy(&(src)->data9[(sat) * fcs_front_xq_aI_data9_size_c], &(dst)->data9[(dat) * fcs_front_xq_aI_data9_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data9_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data9_nmove(&(src)->data9[(sat) * fcs_front_xq_aI_data9_size_c], &(dst)->data9[(dat) * fcs_front_xq_aI_data9_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data9_xchange(e0, e1, t)                  (fcs_front_xq_aI_data9_copy(e0, t), fcs_front_xq_aI_data9_copy(e1, e0), fcs_front_xq_aI_data9_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data9_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data9_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data9_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data9_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data9_assign(src, dst)                 , fcs_front_xq_aI_data9_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_assign_at(src, sat, dst)         , fcs_front_xq_aI_data9_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_null(e)                          , fcs_front_xq_aI_data9_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_inc(e)                           , fcs_front_xq_aI_data9_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_dec(e)                           , fcs_front_xq_aI_data9_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_add(e, n)                        , fcs_front_xq_aI_data9_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_sub(e, n)                        , fcs_front_xq_aI_data9_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_copy(src, dst)                   , fcs_front_xq_aI_data9_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_ncopy(src, dst, n)               , fcs_front_xq_aI_data9_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_nmove(src, dst, n)               , fcs_front_xq_aI_data9_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data9_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data9_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data9_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_xchange(e0, e1, t)               , fcs_front_xq_aI_data9_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data9_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data9_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data9_assign(src, dst)                 fcs_front_xq_aI_data9_assign(src, dst)
 #define fcs_front_xq_aI_cc_data9_assign_at(src, sat, dst)         fcs_front_xq_aI_data9_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data9_null(e)                          fcs_front_xq_aI_data9_null(e)
 #define fcs_front_xq_aI_cc_data9_inc(e)                           fcs_front_xq_aI_data9_inc(e)
 #define fcs_front_xq_aI_cc_data9_dec(e)                           fcs_front_xq_aI_data9_dec(e)
 #define fcs_front_xq_aI_cc_data9_add(e, n)                        fcs_front_xq_aI_data9_add(e, n)
 #define fcs_front_xq_aI_cc_data9_sub(e, n)                        fcs_front_xq_aI_data9_sub(e, n)
 #define fcs_front_xq_aI_cc_data9_copy(src, dst)                   fcs_front_xq_aI_data9_copy(src, dst)
 #define fcs_front_xq_aI_cc_data9_ncopy(src, dst, n)               fcs_front_xq_aI_data9_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data9_nmove(src, dst, n)               fcs_front_xq_aI_data9_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data9_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data9_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data9_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data9_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data9_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data9_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data9_xchange(e0, e1, t)               fcs_front_xq_aI_data9_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data9_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data9_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA9 */

 #define fcs_front_xq_aI_data9_n                                   0
 #define fcs_front_xq_aI_data9_byte                                0
/* #define fcs_front_xq_aI_data9_ptr(e)*/

 #define fcs_front_xq_aI_data9_byte_flex                           0
 #define fcs_front_xq_aI_data9_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data9_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data9_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data9_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data9_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data9_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data9_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data9_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data9_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data9_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data9_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data9_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data9_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data9_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data9_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data9_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data9_assign(src, dst)
 #define fcs_front_xq_aI_cc_data9_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data9_null(e)
 #define fcs_front_xq_aI_cc_data9_inc(e)
 #define fcs_front_xq_aI_cc_data9_dec(e)
 #define fcs_front_xq_aI_cc_data9_add(e, n)
 #define fcs_front_xq_aI_cc_data9_sub(e, n)
 #define fcs_front_xq_aI_cc_data9_copy(src, dst)
 #define fcs_front_xq_aI_cc_data9_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data9_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data9_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data9_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data9_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data9_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data9_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA9 */

#define fcs_front_xq_aI_data9_cm                                   SLCM_DATA9  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA10 fcs_front_xq_aI_SL_DATA10_IGNORE fcs_front_xq_aI_sl_data10_type_c fcs_front_xq_aI_sl_data10_size_c fcs_front_xq_aI_sl_data10_type_mpi fcs_front_xq_aI_sl_data10_size_mpi fcs_front_xq_aI_sl_data10_memcpy fcs_front_xq_aI_sl_data10_weight fcs_front_xq_aI_sl_data10_flex */

#ifdef fcs_front_xq_aI_SL_DATA10

 #define fcs_front_xq_aI_sl_data10_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data10_size_c) * sizeof(fcs_front_xq_aI_sl_data10_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data10_copy
  #if fcs_front_xq_aI_sl_data10_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data10_memcpy)
   #if fcs_front_xq_aI_sl_data10_size_c == 1
    #define fcs_front_xq_aI_sl_data10_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data10_size_c == 2
    #define fcs_front_xq_aI_sl_data10_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data10_size_c == 3
    #define fcs_front_xq_aI_sl_data10_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data10_size_c == 4
    #define fcs_front_xq_aI_sl_data10_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data10_size_c == 5
    #define fcs_front_xq_aI_sl_data10_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data10_size_c == 6
    #define fcs_front_xq_aI_sl_data10_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data10_size_c == 7
    #define fcs_front_xq_aI_sl_data10_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data10_size_c == 8
    #define fcs_front_xq_aI_sl_data10_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data10_size_c == 9
    #define fcs_front_xq_aI_sl_data10_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data10_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data10_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data10_ncopy
  #define fcs_front_xq_aI_sl_data10_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data10_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data10_nmove
  #define fcs_front_xq_aI_sl_data10_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data10_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data10_type_c                              fcs_front_xq_aI_sl_data10_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data10_size_c                              (fcs_front_xq_aI_sl_data10_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data10_type_mpi                            (fcs_front_xq_aI_sl_data10_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data10_size_mpi                            (fcs_front_xq_aI_sl_data10_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data10_idx                                 10  /* sl_macro */

 #define fcs_front_xq_aI_data10_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data10_byte                                (fcs_front_xq_aI_sl_data10_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data10_ptr(e)                              (e)->data10  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data10_flex
 # define fcs_front_xq_aI_data10_byte_flex                          (fcs_front_xq_aI_sl_data10_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data10_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data10_weight
 # define fcs_front_xq_aI_data10_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data10_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data10_assign(src, dst)                    ((dst)->data10 = (src)->data10)  /* sl_macro */
 #define fcs_front_xq_aI_data10_assign_at(src, sat, dst)            ((dst)->data10 = &(src)->data10[(sat) * fcs_front_xq_aI_data10_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data10_null(e)                             ((e)->data10 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data10_inc(e)                              ((e)->data10 += fcs_front_xq_aI_data10_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data10_dec(e)                              ((e)->data10 -= fcs_front_xq_aI_data10_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data10_add(e, n)                           ((e)->data10 += (n) * fcs_front_xq_aI_data10_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data10_sub(e, n)                           ((e)->data10 -= (n) * fcs_front_xq_aI_data10_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data10_copy(src, dst)                      fcs_front_xq_aI_sl_data10_copy((src)->data10, (dst)->data10)  /* sl_macro */
 #define fcs_front_xq_aI_data10_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data10_ncopy((src)->data10, (dst)->data10, n)  /* sl_macro */
 #define fcs_front_xq_aI_data10_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data10_nmove((src)->data10, (dst)->data10, n)  /* sl_macro */

 #define fcs_front_xq_aI_data10_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data10_copy(&(src)->data10[(sat) * fcs_front_xq_aI_data10_size_c], &(dst)->data10[(dat) * fcs_front_xq_aI_data10_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data10_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data10_ncopy(&(src)->data10[(sat) * fcs_front_xq_aI_data10_size_c], &(dst)->data10[(dat) * fcs_front_xq_aI_data10_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data10_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data10_nmove(&(src)->data10[(sat) * fcs_front_xq_aI_data10_size_c], &(dst)->data10[(dat) * fcs_front_xq_aI_data10_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data10_xchange(e0, e1, t)                  (fcs_front_xq_aI_data10_copy(e0, t), fcs_front_xq_aI_data10_copy(e1, e0), fcs_front_xq_aI_data10_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data10_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data10_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data10_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data10_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data10_assign(src, dst)                 , fcs_front_xq_aI_data10_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_assign_at(src, sat, dst)         , fcs_front_xq_aI_data10_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_null(e)                          , fcs_front_xq_aI_data10_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_inc(e)                           , fcs_front_xq_aI_data10_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_dec(e)                           , fcs_front_xq_aI_data10_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_add(e, n)                        , fcs_front_xq_aI_data10_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_sub(e, n)                        , fcs_front_xq_aI_data10_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_copy(src, dst)                   , fcs_front_xq_aI_data10_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_ncopy(src, dst, n)               , fcs_front_xq_aI_data10_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_nmove(src, dst, n)               , fcs_front_xq_aI_data10_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data10_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data10_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data10_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_xchange(e0, e1, t)               , fcs_front_xq_aI_data10_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data10_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data10_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data10_assign(src, dst)                 fcs_front_xq_aI_data10_assign(src, dst)
 #define fcs_front_xq_aI_cc_data10_assign_at(src, sat, dst)         fcs_front_xq_aI_data10_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data10_null(e)                          fcs_front_xq_aI_data10_null(e)
 #define fcs_front_xq_aI_cc_data10_inc(e)                           fcs_front_xq_aI_data10_inc(e)
 #define fcs_front_xq_aI_cc_data10_dec(e)                           fcs_front_xq_aI_data10_dec(e)
 #define fcs_front_xq_aI_cc_data10_add(e, n)                        fcs_front_xq_aI_data10_add(e, n)
 #define fcs_front_xq_aI_cc_data10_sub(e, n)                        fcs_front_xq_aI_data10_sub(e, n)
 #define fcs_front_xq_aI_cc_data10_copy(src, dst)                   fcs_front_xq_aI_data10_copy(src, dst)
 #define fcs_front_xq_aI_cc_data10_ncopy(src, dst, n)               fcs_front_xq_aI_data10_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data10_nmove(src, dst, n)               fcs_front_xq_aI_data10_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data10_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data10_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data10_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data10_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data10_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data10_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data10_xchange(e0, e1, t)               fcs_front_xq_aI_data10_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data10_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data10_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA10 */

 #define fcs_front_xq_aI_data10_n                                   0
 #define fcs_front_xq_aI_data10_byte                                0
/* #define fcs_front_xq_aI_data10_ptr(e)*/

 #define fcs_front_xq_aI_data10_byte_flex                           0
 #define fcs_front_xq_aI_data10_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data10_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data10_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data10_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data10_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data10_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data10_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data10_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data10_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data10_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data10_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data10_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data10_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data10_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data10_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data10_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data10_assign(src, dst)
 #define fcs_front_xq_aI_cc_data10_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data10_null(e)
 #define fcs_front_xq_aI_cc_data10_inc(e)
 #define fcs_front_xq_aI_cc_data10_dec(e)
 #define fcs_front_xq_aI_cc_data10_add(e, n)
 #define fcs_front_xq_aI_cc_data10_sub(e, n)
 #define fcs_front_xq_aI_cc_data10_copy(src, dst)
 #define fcs_front_xq_aI_cc_data10_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data10_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data10_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data10_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data10_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data10_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data10_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA10 */

#define fcs_front_xq_aI_data10_cm                                   SLCM_DATA10  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA11 fcs_front_xq_aI_SL_DATA11_IGNORE fcs_front_xq_aI_sl_data11_type_c fcs_front_xq_aI_sl_data11_size_c fcs_front_xq_aI_sl_data11_type_mpi fcs_front_xq_aI_sl_data11_size_mpi fcs_front_xq_aI_sl_data11_memcpy fcs_front_xq_aI_sl_data11_weight fcs_front_xq_aI_sl_data11_flex */

#ifdef fcs_front_xq_aI_SL_DATA11

 #define fcs_front_xq_aI_sl_data11_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data11_size_c) * sizeof(fcs_front_xq_aI_sl_data11_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data11_copy
  #if fcs_front_xq_aI_sl_data11_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data11_memcpy)
   #if fcs_front_xq_aI_sl_data11_size_c == 1
    #define fcs_front_xq_aI_sl_data11_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data11_size_c == 2
    #define fcs_front_xq_aI_sl_data11_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data11_size_c == 3
    #define fcs_front_xq_aI_sl_data11_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data11_size_c == 4
    #define fcs_front_xq_aI_sl_data11_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data11_size_c == 5
    #define fcs_front_xq_aI_sl_data11_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data11_size_c == 6
    #define fcs_front_xq_aI_sl_data11_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data11_size_c == 7
    #define fcs_front_xq_aI_sl_data11_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data11_size_c == 8
    #define fcs_front_xq_aI_sl_data11_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data11_size_c == 9
    #define fcs_front_xq_aI_sl_data11_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data11_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data11_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data11_ncopy
  #define fcs_front_xq_aI_sl_data11_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data11_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data11_nmove
  #define fcs_front_xq_aI_sl_data11_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data11_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data11_type_c                              fcs_front_xq_aI_sl_data11_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data11_size_c                              (fcs_front_xq_aI_sl_data11_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data11_type_mpi                            (fcs_front_xq_aI_sl_data11_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data11_size_mpi                            (fcs_front_xq_aI_sl_data11_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data11_idx                                 11  /* sl_macro */

 #define fcs_front_xq_aI_data11_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data11_byte                                (fcs_front_xq_aI_sl_data11_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data11_ptr(e)                              (e)->data11  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data11_flex
 # define fcs_front_xq_aI_data11_byte_flex                          (fcs_front_xq_aI_sl_data11_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data11_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data11_weight
 # define fcs_front_xq_aI_data11_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data11_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data11_assign(src, dst)                    ((dst)->data11 = (src)->data11)  /* sl_macro */
 #define fcs_front_xq_aI_data11_assign_at(src, sat, dst)            ((dst)->data11 = &(src)->data11[(sat) * fcs_front_xq_aI_data11_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data11_null(e)                             ((e)->data11 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data11_inc(e)                              ((e)->data11 += fcs_front_xq_aI_data11_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data11_dec(e)                              ((e)->data11 -= fcs_front_xq_aI_data11_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data11_add(e, n)                           ((e)->data11 += (n) * fcs_front_xq_aI_data11_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data11_sub(e, n)                           ((e)->data11 -= (n) * fcs_front_xq_aI_data11_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data11_copy(src, dst)                      fcs_front_xq_aI_sl_data11_copy((src)->data11, (dst)->data11)  /* sl_macro */
 #define fcs_front_xq_aI_data11_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data11_ncopy((src)->data11, (dst)->data11, n)  /* sl_macro */
 #define fcs_front_xq_aI_data11_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data11_nmove((src)->data11, (dst)->data11, n)  /* sl_macro */

 #define fcs_front_xq_aI_data11_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data11_copy(&(src)->data11[(sat) * fcs_front_xq_aI_data11_size_c], &(dst)->data11[(dat) * fcs_front_xq_aI_data11_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data11_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data11_ncopy(&(src)->data11[(sat) * fcs_front_xq_aI_data11_size_c], &(dst)->data11[(dat) * fcs_front_xq_aI_data11_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data11_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data11_nmove(&(src)->data11[(sat) * fcs_front_xq_aI_data11_size_c], &(dst)->data11[(dat) * fcs_front_xq_aI_data11_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data11_xchange(e0, e1, t)                  (fcs_front_xq_aI_data11_copy(e0, t), fcs_front_xq_aI_data11_copy(e1, e0), fcs_front_xq_aI_data11_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data11_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data11_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data11_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data11_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data11_assign(src, dst)                 , fcs_front_xq_aI_data11_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_assign_at(src, sat, dst)         , fcs_front_xq_aI_data11_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_null(e)                          , fcs_front_xq_aI_data11_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_inc(e)                           , fcs_front_xq_aI_data11_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_dec(e)                           , fcs_front_xq_aI_data11_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_add(e, n)                        , fcs_front_xq_aI_data11_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_sub(e, n)                        , fcs_front_xq_aI_data11_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_copy(src, dst)                   , fcs_front_xq_aI_data11_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_ncopy(src, dst, n)               , fcs_front_xq_aI_data11_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_nmove(src, dst, n)               , fcs_front_xq_aI_data11_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data11_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data11_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data11_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_xchange(e0, e1, t)               , fcs_front_xq_aI_data11_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data11_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data11_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data11_assign(src, dst)                 fcs_front_xq_aI_data11_assign(src, dst)
 #define fcs_front_xq_aI_cc_data11_assign_at(src, sat, dst)         fcs_front_xq_aI_data11_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data11_null(e)                          fcs_front_xq_aI_data11_null(e)
 #define fcs_front_xq_aI_cc_data11_inc(e)                           fcs_front_xq_aI_data11_inc(e)
 #define fcs_front_xq_aI_cc_data11_dec(e)                           fcs_front_xq_aI_data11_dec(e)
 #define fcs_front_xq_aI_cc_data11_add(e, n)                        fcs_front_xq_aI_data11_add(e, n)
 #define fcs_front_xq_aI_cc_data11_sub(e, n)                        fcs_front_xq_aI_data11_sub(e, n)
 #define fcs_front_xq_aI_cc_data11_copy(src, dst)                   fcs_front_xq_aI_data11_copy(src, dst)
 #define fcs_front_xq_aI_cc_data11_ncopy(src, dst, n)               fcs_front_xq_aI_data11_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data11_nmove(src, dst, n)               fcs_front_xq_aI_data11_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data11_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data11_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data11_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data11_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data11_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data11_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data11_xchange(e0, e1, t)               fcs_front_xq_aI_data11_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data11_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data11_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA11 */

 #define fcs_front_xq_aI_data11_n                                   0
 #define fcs_front_xq_aI_data11_byte                                0
/* #define fcs_front_xq_aI_data11_ptr(e)*/

 #define fcs_front_xq_aI_data11_byte_flex                           0
 #define fcs_front_xq_aI_data11_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data11_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data11_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data11_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data11_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data11_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data11_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data11_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data11_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data11_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data11_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data11_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data11_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data11_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data11_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data11_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data11_assign(src, dst)
 #define fcs_front_xq_aI_cc_data11_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data11_null(e)
 #define fcs_front_xq_aI_cc_data11_inc(e)
 #define fcs_front_xq_aI_cc_data11_dec(e)
 #define fcs_front_xq_aI_cc_data11_add(e, n)
 #define fcs_front_xq_aI_cc_data11_sub(e, n)
 #define fcs_front_xq_aI_cc_data11_copy(src, dst)
 #define fcs_front_xq_aI_cc_data11_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data11_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data11_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data11_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data11_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data11_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data11_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA11 */

#define fcs_front_xq_aI_data11_cm                                   SLCM_DATA11  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA12 fcs_front_xq_aI_SL_DATA12_IGNORE fcs_front_xq_aI_sl_data12_type_c fcs_front_xq_aI_sl_data12_size_c fcs_front_xq_aI_sl_data12_type_mpi fcs_front_xq_aI_sl_data12_size_mpi fcs_front_xq_aI_sl_data12_memcpy fcs_front_xq_aI_sl_data12_weight fcs_front_xq_aI_sl_data12_flex */

#ifdef fcs_front_xq_aI_SL_DATA12

 #define fcs_front_xq_aI_sl_data12_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data12_size_c) * sizeof(fcs_front_xq_aI_sl_data12_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data12_copy
  #if fcs_front_xq_aI_sl_data12_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data12_memcpy)
   #if fcs_front_xq_aI_sl_data12_size_c == 1
    #define fcs_front_xq_aI_sl_data12_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data12_size_c == 2
    #define fcs_front_xq_aI_sl_data12_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data12_size_c == 3
    #define fcs_front_xq_aI_sl_data12_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data12_size_c == 4
    #define fcs_front_xq_aI_sl_data12_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data12_size_c == 5
    #define fcs_front_xq_aI_sl_data12_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data12_size_c == 6
    #define fcs_front_xq_aI_sl_data12_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data12_size_c == 7
    #define fcs_front_xq_aI_sl_data12_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data12_size_c == 8
    #define fcs_front_xq_aI_sl_data12_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data12_size_c == 9
    #define fcs_front_xq_aI_sl_data12_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data12_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data12_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data12_ncopy
  #define fcs_front_xq_aI_sl_data12_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data12_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data12_nmove
  #define fcs_front_xq_aI_sl_data12_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data12_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data12_type_c                              fcs_front_xq_aI_sl_data12_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data12_size_c                              (fcs_front_xq_aI_sl_data12_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data12_type_mpi                            (fcs_front_xq_aI_sl_data12_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data12_size_mpi                            (fcs_front_xq_aI_sl_data12_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data12_idx                                 12  /* sl_macro */

 #define fcs_front_xq_aI_data12_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data12_byte                                (fcs_front_xq_aI_sl_data12_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data12_ptr(e)                              (e)->data12  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data12_flex
 # define fcs_front_xq_aI_data12_byte_flex                          (fcs_front_xq_aI_sl_data12_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data12_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data12_weight
 # define fcs_front_xq_aI_data12_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data12_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data12_assign(src, dst)                    ((dst)->data12 = (src)->data12)  /* sl_macro */
 #define fcs_front_xq_aI_data12_assign_at(src, sat, dst)            ((dst)->data12 = &(src)->data12[(sat) * fcs_front_xq_aI_data12_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data12_null(e)                             ((e)->data12 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data12_inc(e)                              ((e)->data12 += fcs_front_xq_aI_data12_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data12_dec(e)                              ((e)->data12 -= fcs_front_xq_aI_data12_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data12_add(e, n)                           ((e)->data12 += (n) * fcs_front_xq_aI_data12_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data12_sub(e, n)                           ((e)->data12 -= (n) * fcs_front_xq_aI_data12_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data12_copy(src, dst)                      fcs_front_xq_aI_sl_data12_copy((src)->data12, (dst)->data12)  /* sl_macro */
 #define fcs_front_xq_aI_data12_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data12_ncopy((src)->data12, (dst)->data12, n)  /* sl_macro */
 #define fcs_front_xq_aI_data12_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data12_nmove((src)->data12, (dst)->data12, n)  /* sl_macro */

 #define fcs_front_xq_aI_data12_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data12_copy(&(src)->data12[(sat) * fcs_front_xq_aI_data12_size_c], &(dst)->data12[(dat) * fcs_front_xq_aI_data12_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data12_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data12_ncopy(&(src)->data12[(sat) * fcs_front_xq_aI_data12_size_c], &(dst)->data12[(dat) * fcs_front_xq_aI_data12_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data12_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data12_nmove(&(src)->data12[(sat) * fcs_front_xq_aI_data12_size_c], &(dst)->data12[(dat) * fcs_front_xq_aI_data12_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data12_xchange(e0, e1, t)                  (fcs_front_xq_aI_data12_copy(e0, t), fcs_front_xq_aI_data12_copy(e1, e0), fcs_front_xq_aI_data12_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data12_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data12_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data12_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data12_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data12_assign(src, dst)                 , fcs_front_xq_aI_data12_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_assign_at(src, sat, dst)         , fcs_front_xq_aI_data12_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_null(e)                          , fcs_front_xq_aI_data12_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_inc(e)                           , fcs_front_xq_aI_data12_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_dec(e)                           , fcs_front_xq_aI_data12_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_add(e, n)                        , fcs_front_xq_aI_data12_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_sub(e, n)                        , fcs_front_xq_aI_data12_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_copy(src, dst)                   , fcs_front_xq_aI_data12_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_ncopy(src, dst, n)               , fcs_front_xq_aI_data12_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_nmove(src, dst, n)               , fcs_front_xq_aI_data12_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data12_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data12_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data12_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_xchange(e0, e1, t)               , fcs_front_xq_aI_data12_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data12_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data12_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data12_assign(src, dst)                 fcs_front_xq_aI_data12_assign(src, dst)
 #define fcs_front_xq_aI_cc_data12_assign_at(src, sat, dst)         fcs_front_xq_aI_data12_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data12_null(e)                          fcs_front_xq_aI_data12_null(e)
 #define fcs_front_xq_aI_cc_data12_inc(e)                           fcs_front_xq_aI_data12_inc(e)
 #define fcs_front_xq_aI_cc_data12_dec(e)                           fcs_front_xq_aI_data12_dec(e)
 #define fcs_front_xq_aI_cc_data12_add(e, n)                        fcs_front_xq_aI_data12_add(e, n)
 #define fcs_front_xq_aI_cc_data12_sub(e, n)                        fcs_front_xq_aI_data12_sub(e, n)
 #define fcs_front_xq_aI_cc_data12_copy(src, dst)                   fcs_front_xq_aI_data12_copy(src, dst)
 #define fcs_front_xq_aI_cc_data12_ncopy(src, dst, n)               fcs_front_xq_aI_data12_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data12_nmove(src, dst, n)               fcs_front_xq_aI_data12_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data12_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data12_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data12_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data12_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data12_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data12_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data12_xchange(e0, e1, t)               fcs_front_xq_aI_data12_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data12_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data12_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA12 */

 #define fcs_front_xq_aI_data12_n                                   0
 #define fcs_front_xq_aI_data12_byte                                0
/* #define fcs_front_xq_aI_data12_ptr(e)*/

 #define fcs_front_xq_aI_data12_byte_flex                           0
 #define fcs_front_xq_aI_data12_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data12_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data12_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data12_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data12_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data12_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data12_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data12_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data12_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data12_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data12_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data12_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data12_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data12_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data12_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data12_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data12_assign(src, dst)
 #define fcs_front_xq_aI_cc_data12_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data12_null(e)
 #define fcs_front_xq_aI_cc_data12_inc(e)
 #define fcs_front_xq_aI_cc_data12_dec(e)
 #define fcs_front_xq_aI_cc_data12_add(e, n)
 #define fcs_front_xq_aI_cc_data12_sub(e, n)
 #define fcs_front_xq_aI_cc_data12_copy(src, dst)
 #define fcs_front_xq_aI_cc_data12_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data12_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data12_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data12_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data12_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data12_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data12_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA12 */

#define fcs_front_xq_aI_data12_cm                                   SLCM_DATA12  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA13 fcs_front_xq_aI_SL_DATA13_IGNORE fcs_front_xq_aI_sl_data13_type_c fcs_front_xq_aI_sl_data13_size_c fcs_front_xq_aI_sl_data13_type_mpi fcs_front_xq_aI_sl_data13_size_mpi fcs_front_xq_aI_sl_data13_memcpy fcs_front_xq_aI_sl_data13_weight fcs_front_xq_aI_sl_data13_flex */

#ifdef fcs_front_xq_aI_SL_DATA13

 #define fcs_front_xq_aI_sl_data13_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data13_size_c) * sizeof(fcs_front_xq_aI_sl_data13_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data13_copy
  #if fcs_front_xq_aI_sl_data13_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data13_memcpy)
   #if fcs_front_xq_aI_sl_data13_size_c == 1
    #define fcs_front_xq_aI_sl_data13_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data13_size_c == 2
    #define fcs_front_xq_aI_sl_data13_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data13_size_c == 3
    #define fcs_front_xq_aI_sl_data13_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data13_size_c == 4
    #define fcs_front_xq_aI_sl_data13_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data13_size_c == 5
    #define fcs_front_xq_aI_sl_data13_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data13_size_c == 6
    #define fcs_front_xq_aI_sl_data13_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data13_size_c == 7
    #define fcs_front_xq_aI_sl_data13_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data13_size_c == 8
    #define fcs_front_xq_aI_sl_data13_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data13_size_c == 9
    #define fcs_front_xq_aI_sl_data13_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data13_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data13_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data13_ncopy
  #define fcs_front_xq_aI_sl_data13_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data13_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data13_nmove
  #define fcs_front_xq_aI_sl_data13_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data13_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data13_type_c                              fcs_front_xq_aI_sl_data13_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data13_size_c                              (fcs_front_xq_aI_sl_data13_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data13_type_mpi                            (fcs_front_xq_aI_sl_data13_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data13_size_mpi                            (fcs_front_xq_aI_sl_data13_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data13_idx                                 13  /* sl_macro */

 #define fcs_front_xq_aI_data13_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data13_byte                                (fcs_front_xq_aI_sl_data13_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data13_ptr(e)                              (e)->data13  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data13_flex
 # define fcs_front_xq_aI_data13_byte_flex                          (fcs_front_xq_aI_sl_data13_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data13_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data13_weight
 # define fcs_front_xq_aI_data13_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data13_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data13_assign(src, dst)                    ((dst)->data13 = (src)->data13)  /* sl_macro */
 #define fcs_front_xq_aI_data13_assign_at(src, sat, dst)            ((dst)->data13 = &(src)->data13[(sat) * fcs_front_xq_aI_data13_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data13_null(e)                             ((e)->data13 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data13_inc(e)                              ((e)->data13 += fcs_front_xq_aI_data13_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data13_dec(e)                              ((e)->data13 -= fcs_front_xq_aI_data13_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data13_add(e, n)                           ((e)->data13 += (n) * fcs_front_xq_aI_data13_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data13_sub(e, n)                           ((e)->data13 -= (n) * fcs_front_xq_aI_data13_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data13_copy(src, dst)                      fcs_front_xq_aI_sl_data13_copy((src)->data13, (dst)->data13)  /* sl_macro */
 #define fcs_front_xq_aI_data13_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data13_ncopy((src)->data13, (dst)->data13, n)  /* sl_macro */
 #define fcs_front_xq_aI_data13_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data13_nmove((src)->data13, (dst)->data13, n)  /* sl_macro */

 #define fcs_front_xq_aI_data13_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data13_copy(&(src)->data13[(sat) * fcs_front_xq_aI_data13_size_c], &(dst)->data13[(dat) * fcs_front_xq_aI_data13_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data13_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data13_ncopy(&(src)->data13[(sat) * fcs_front_xq_aI_data13_size_c], &(dst)->data13[(dat) * fcs_front_xq_aI_data13_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data13_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data13_nmove(&(src)->data13[(sat) * fcs_front_xq_aI_data13_size_c], &(dst)->data13[(dat) * fcs_front_xq_aI_data13_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data13_xchange(e0, e1, t)                  (fcs_front_xq_aI_data13_copy(e0, t), fcs_front_xq_aI_data13_copy(e1, e0), fcs_front_xq_aI_data13_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data13_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data13_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data13_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data13_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data13_assign(src, dst)                 , fcs_front_xq_aI_data13_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_assign_at(src, sat, dst)         , fcs_front_xq_aI_data13_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_null(e)                          , fcs_front_xq_aI_data13_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_inc(e)                           , fcs_front_xq_aI_data13_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_dec(e)                           , fcs_front_xq_aI_data13_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_add(e, n)                        , fcs_front_xq_aI_data13_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_sub(e, n)                        , fcs_front_xq_aI_data13_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_copy(src, dst)                   , fcs_front_xq_aI_data13_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_ncopy(src, dst, n)               , fcs_front_xq_aI_data13_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_nmove(src, dst, n)               , fcs_front_xq_aI_data13_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data13_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data13_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data13_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_xchange(e0, e1, t)               , fcs_front_xq_aI_data13_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data13_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data13_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data13_assign(src, dst)                 fcs_front_xq_aI_data13_assign(src, dst)
 #define fcs_front_xq_aI_cc_data13_assign_at(src, sat, dst)         fcs_front_xq_aI_data13_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data13_null(e)                          fcs_front_xq_aI_data13_null(e)
 #define fcs_front_xq_aI_cc_data13_inc(e)                           fcs_front_xq_aI_data13_inc(e)
 #define fcs_front_xq_aI_cc_data13_dec(e)                           fcs_front_xq_aI_data13_dec(e)
 #define fcs_front_xq_aI_cc_data13_add(e, n)                        fcs_front_xq_aI_data13_add(e, n)
 #define fcs_front_xq_aI_cc_data13_sub(e, n)                        fcs_front_xq_aI_data13_sub(e, n)
 #define fcs_front_xq_aI_cc_data13_copy(src, dst)                   fcs_front_xq_aI_data13_copy(src, dst)
 #define fcs_front_xq_aI_cc_data13_ncopy(src, dst, n)               fcs_front_xq_aI_data13_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data13_nmove(src, dst, n)               fcs_front_xq_aI_data13_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data13_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data13_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data13_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data13_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data13_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data13_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data13_xchange(e0, e1, t)               fcs_front_xq_aI_data13_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data13_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data13_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA13 */

 #define fcs_front_xq_aI_data13_n                                   0
 #define fcs_front_xq_aI_data13_byte                                0
/* #define fcs_front_xq_aI_data13_ptr(e)*/

 #define fcs_front_xq_aI_data13_byte_flex                           0
 #define fcs_front_xq_aI_data13_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data13_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data13_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data13_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data13_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data13_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data13_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data13_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data13_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data13_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data13_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data13_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data13_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data13_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data13_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data13_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data13_assign(src, dst)
 #define fcs_front_xq_aI_cc_data13_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data13_null(e)
 #define fcs_front_xq_aI_cc_data13_inc(e)
 #define fcs_front_xq_aI_cc_data13_dec(e)
 #define fcs_front_xq_aI_cc_data13_add(e, n)
 #define fcs_front_xq_aI_cc_data13_sub(e, n)
 #define fcs_front_xq_aI_cc_data13_copy(src, dst)
 #define fcs_front_xq_aI_cc_data13_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data13_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data13_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data13_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data13_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data13_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data13_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA13 */

#define fcs_front_xq_aI_data13_cm                                   SLCM_DATA13  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA14 fcs_front_xq_aI_SL_DATA14_IGNORE fcs_front_xq_aI_sl_data14_type_c fcs_front_xq_aI_sl_data14_size_c fcs_front_xq_aI_sl_data14_type_mpi fcs_front_xq_aI_sl_data14_size_mpi fcs_front_xq_aI_sl_data14_memcpy fcs_front_xq_aI_sl_data14_weight fcs_front_xq_aI_sl_data14_flex */

#ifdef fcs_front_xq_aI_SL_DATA14

 #define fcs_front_xq_aI_sl_data14_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data14_size_c) * sizeof(fcs_front_xq_aI_sl_data14_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data14_copy
  #if fcs_front_xq_aI_sl_data14_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data14_memcpy)
   #if fcs_front_xq_aI_sl_data14_size_c == 1
    #define fcs_front_xq_aI_sl_data14_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data14_size_c == 2
    #define fcs_front_xq_aI_sl_data14_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data14_size_c == 3
    #define fcs_front_xq_aI_sl_data14_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data14_size_c == 4
    #define fcs_front_xq_aI_sl_data14_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data14_size_c == 5
    #define fcs_front_xq_aI_sl_data14_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data14_size_c == 6
    #define fcs_front_xq_aI_sl_data14_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data14_size_c == 7
    #define fcs_front_xq_aI_sl_data14_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data14_size_c == 8
    #define fcs_front_xq_aI_sl_data14_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data14_size_c == 9
    #define fcs_front_xq_aI_sl_data14_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data14_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data14_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data14_ncopy
  #define fcs_front_xq_aI_sl_data14_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data14_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data14_nmove
  #define fcs_front_xq_aI_sl_data14_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data14_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data14_type_c                              fcs_front_xq_aI_sl_data14_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data14_size_c                              (fcs_front_xq_aI_sl_data14_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data14_type_mpi                            (fcs_front_xq_aI_sl_data14_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data14_size_mpi                            (fcs_front_xq_aI_sl_data14_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data14_idx                                 14  /* sl_macro */

 #define fcs_front_xq_aI_data14_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data14_byte                                (fcs_front_xq_aI_sl_data14_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data14_ptr(e)                              (e)->data14  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data14_flex
 # define fcs_front_xq_aI_data14_byte_flex                          (fcs_front_xq_aI_sl_data14_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data14_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data14_weight
 # define fcs_front_xq_aI_data14_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data14_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data14_assign(src, dst)                    ((dst)->data14 = (src)->data14)  /* sl_macro */
 #define fcs_front_xq_aI_data14_assign_at(src, sat, dst)            ((dst)->data14 = &(src)->data14[(sat) * fcs_front_xq_aI_data14_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data14_null(e)                             ((e)->data14 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data14_inc(e)                              ((e)->data14 += fcs_front_xq_aI_data14_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data14_dec(e)                              ((e)->data14 -= fcs_front_xq_aI_data14_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data14_add(e, n)                           ((e)->data14 += (n) * fcs_front_xq_aI_data14_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data14_sub(e, n)                           ((e)->data14 -= (n) * fcs_front_xq_aI_data14_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data14_copy(src, dst)                      fcs_front_xq_aI_sl_data14_copy((src)->data14, (dst)->data14)  /* sl_macro */
 #define fcs_front_xq_aI_data14_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data14_ncopy((src)->data14, (dst)->data14, n)  /* sl_macro */
 #define fcs_front_xq_aI_data14_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data14_nmove((src)->data14, (dst)->data14, n)  /* sl_macro */

 #define fcs_front_xq_aI_data14_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data14_copy(&(src)->data14[(sat) * fcs_front_xq_aI_data14_size_c], &(dst)->data14[(dat) * fcs_front_xq_aI_data14_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data14_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data14_ncopy(&(src)->data14[(sat) * fcs_front_xq_aI_data14_size_c], &(dst)->data14[(dat) * fcs_front_xq_aI_data14_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data14_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data14_nmove(&(src)->data14[(sat) * fcs_front_xq_aI_data14_size_c], &(dst)->data14[(dat) * fcs_front_xq_aI_data14_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data14_xchange(e0, e1, t)                  (fcs_front_xq_aI_data14_copy(e0, t), fcs_front_xq_aI_data14_copy(e1, e0), fcs_front_xq_aI_data14_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data14_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data14_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data14_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data14_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data14_assign(src, dst)                 , fcs_front_xq_aI_data14_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_assign_at(src, sat, dst)         , fcs_front_xq_aI_data14_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_null(e)                          , fcs_front_xq_aI_data14_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_inc(e)                           , fcs_front_xq_aI_data14_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_dec(e)                           , fcs_front_xq_aI_data14_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_add(e, n)                        , fcs_front_xq_aI_data14_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_sub(e, n)                        , fcs_front_xq_aI_data14_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_copy(src, dst)                   , fcs_front_xq_aI_data14_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_ncopy(src, dst, n)               , fcs_front_xq_aI_data14_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_nmove(src, dst, n)               , fcs_front_xq_aI_data14_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data14_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data14_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data14_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_xchange(e0, e1, t)               , fcs_front_xq_aI_data14_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data14_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data14_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data14_assign(src, dst)                 fcs_front_xq_aI_data14_assign(src, dst)
 #define fcs_front_xq_aI_cc_data14_assign_at(src, sat, dst)         fcs_front_xq_aI_data14_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data14_null(e)                          fcs_front_xq_aI_data14_null(e)
 #define fcs_front_xq_aI_cc_data14_inc(e)                           fcs_front_xq_aI_data14_inc(e)
 #define fcs_front_xq_aI_cc_data14_dec(e)                           fcs_front_xq_aI_data14_dec(e)
 #define fcs_front_xq_aI_cc_data14_add(e, n)                        fcs_front_xq_aI_data14_add(e, n)
 #define fcs_front_xq_aI_cc_data14_sub(e, n)                        fcs_front_xq_aI_data14_sub(e, n)
 #define fcs_front_xq_aI_cc_data14_copy(src, dst)                   fcs_front_xq_aI_data14_copy(src, dst)
 #define fcs_front_xq_aI_cc_data14_ncopy(src, dst, n)               fcs_front_xq_aI_data14_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data14_nmove(src, dst, n)               fcs_front_xq_aI_data14_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data14_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data14_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data14_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data14_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data14_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data14_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data14_xchange(e0, e1, t)               fcs_front_xq_aI_data14_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data14_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data14_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA14 */

 #define fcs_front_xq_aI_data14_n                                   0
 #define fcs_front_xq_aI_data14_byte                                0
/* #define fcs_front_xq_aI_data14_ptr(e)*/

 #define fcs_front_xq_aI_data14_byte_flex                           0
 #define fcs_front_xq_aI_data14_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data14_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data14_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data14_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data14_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data14_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data14_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data14_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data14_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data14_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data14_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data14_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data14_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data14_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data14_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data14_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data14_assign(src, dst)
 #define fcs_front_xq_aI_cc_data14_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data14_null(e)
 #define fcs_front_xq_aI_cc_data14_inc(e)
 #define fcs_front_xq_aI_cc_data14_dec(e)
 #define fcs_front_xq_aI_cc_data14_add(e, n)
 #define fcs_front_xq_aI_cc_data14_sub(e, n)
 #define fcs_front_xq_aI_cc_data14_copy(src, dst)
 #define fcs_front_xq_aI_cc_data14_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data14_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data14_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data14_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data14_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data14_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data14_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA14 */

#define fcs_front_xq_aI_data14_cm                                   SLCM_DATA14  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA15 fcs_front_xq_aI_SL_DATA15_IGNORE fcs_front_xq_aI_sl_data15_type_c fcs_front_xq_aI_sl_data15_size_c fcs_front_xq_aI_sl_data15_type_mpi fcs_front_xq_aI_sl_data15_size_mpi fcs_front_xq_aI_sl_data15_memcpy fcs_front_xq_aI_sl_data15_weight fcs_front_xq_aI_sl_data15_flex */

#ifdef fcs_front_xq_aI_SL_DATA15

 #define fcs_front_xq_aI_sl_data15_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data15_size_c) * sizeof(fcs_front_xq_aI_sl_data15_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data15_copy
  #if fcs_front_xq_aI_sl_data15_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data15_memcpy)
   #if fcs_front_xq_aI_sl_data15_size_c == 1
    #define fcs_front_xq_aI_sl_data15_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data15_size_c == 2
    #define fcs_front_xq_aI_sl_data15_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data15_size_c == 3
    #define fcs_front_xq_aI_sl_data15_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data15_size_c == 4
    #define fcs_front_xq_aI_sl_data15_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data15_size_c == 5
    #define fcs_front_xq_aI_sl_data15_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data15_size_c == 6
    #define fcs_front_xq_aI_sl_data15_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data15_size_c == 7
    #define fcs_front_xq_aI_sl_data15_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data15_size_c == 8
    #define fcs_front_xq_aI_sl_data15_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data15_size_c == 9
    #define fcs_front_xq_aI_sl_data15_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data15_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data15_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data15_ncopy
  #define fcs_front_xq_aI_sl_data15_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data15_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data15_nmove
  #define fcs_front_xq_aI_sl_data15_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data15_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data15_type_c                              fcs_front_xq_aI_sl_data15_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data15_size_c                              (fcs_front_xq_aI_sl_data15_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data15_type_mpi                            (fcs_front_xq_aI_sl_data15_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data15_size_mpi                            (fcs_front_xq_aI_sl_data15_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data15_idx                                 15  /* sl_macro */

 #define fcs_front_xq_aI_data15_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data15_byte                                (fcs_front_xq_aI_sl_data15_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data15_ptr(e)                              (e)->data15  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data15_flex
 # define fcs_front_xq_aI_data15_byte_flex                          (fcs_front_xq_aI_sl_data15_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data15_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data15_weight
 # define fcs_front_xq_aI_data15_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data15_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data15_assign(src, dst)                    ((dst)->data15 = (src)->data15)  /* sl_macro */
 #define fcs_front_xq_aI_data15_assign_at(src, sat, dst)            ((dst)->data15 = &(src)->data15[(sat) * fcs_front_xq_aI_data15_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data15_null(e)                             ((e)->data15 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data15_inc(e)                              ((e)->data15 += fcs_front_xq_aI_data15_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data15_dec(e)                              ((e)->data15 -= fcs_front_xq_aI_data15_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data15_add(e, n)                           ((e)->data15 += (n) * fcs_front_xq_aI_data15_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data15_sub(e, n)                           ((e)->data15 -= (n) * fcs_front_xq_aI_data15_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data15_copy(src, dst)                      fcs_front_xq_aI_sl_data15_copy((src)->data15, (dst)->data15)  /* sl_macro */
 #define fcs_front_xq_aI_data15_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data15_ncopy((src)->data15, (dst)->data15, n)  /* sl_macro */
 #define fcs_front_xq_aI_data15_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data15_nmove((src)->data15, (dst)->data15, n)  /* sl_macro */

 #define fcs_front_xq_aI_data15_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data15_copy(&(src)->data15[(sat) * fcs_front_xq_aI_data15_size_c], &(dst)->data15[(dat) * fcs_front_xq_aI_data15_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data15_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data15_ncopy(&(src)->data15[(sat) * fcs_front_xq_aI_data15_size_c], &(dst)->data15[(dat) * fcs_front_xq_aI_data15_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data15_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data15_nmove(&(src)->data15[(sat) * fcs_front_xq_aI_data15_size_c], &(dst)->data15[(dat) * fcs_front_xq_aI_data15_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data15_xchange(e0, e1, t)                  (fcs_front_xq_aI_data15_copy(e0, t), fcs_front_xq_aI_data15_copy(e1, e0), fcs_front_xq_aI_data15_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data15_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data15_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data15_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data15_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data15_assign(src, dst)                 , fcs_front_xq_aI_data15_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_assign_at(src, sat, dst)         , fcs_front_xq_aI_data15_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_null(e)                          , fcs_front_xq_aI_data15_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_inc(e)                           , fcs_front_xq_aI_data15_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_dec(e)                           , fcs_front_xq_aI_data15_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_add(e, n)                        , fcs_front_xq_aI_data15_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_sub(e, n)                        , fcs_front_xq_aI_data15_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_copy(src, dst)                   , fcs_front_xq_aI_data15_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_ncopy(src, dst, n)               , fcs_front_xq_aI_data15_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_nmove(src, dst, n)               , fcs_front_xq_aI_data15_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data15_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data15_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data15_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_xchange(e0, e1, t)               , fcs_front_xq_aI_data15_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data15_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data15_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data15_assign(src, dst)                 fcs_front_xq_aI_data15_assign(src, dst)
 #define fcs_front_xq_aI_cc_data15_assign_at(src, sat, dst)         fcs_front_xq_aI_data15_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data15_null(e)                          fcs_front_xq_aI_data15_null(e)
 #define fcs_front_xq_aI_cc_data15_inc(e)                           fcs_front_xq_aI_data15_inc(e)
 #define fcs_front_xq_aI_cc_data15_dec(e)                           fcs_front_xq_aI_data15_dec(e)
 #define fcs_front_xq_aI_cc_data15_add(e, n)                        fcs_front_xq_aI_data15_add(e, n)
 #define fcs_front_xq_aI_cc_data15_sub(e, n)                        fcs_front_xq_aI_data15_sub(e, n)
 #define fcs_front_xq_aI_cc_data15_copy(src, dst)                   fcs_front_xq_aI_data15_copy(src, dst)
 #define fcs_front_xq_aI_cc_data15_ncopy(src, dst, n)               fcs_front_xq_aI_data15_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data15_nmove(src, dst, n)               fcs_front_xq_aI_data15_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data15_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data15_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data15_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data15_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data15_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data15_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data15_xchange(e0, e1, t)               fcs_front_xq_aI_data15_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data15_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data15_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA15 */

 #define fcs_front_xq_aI_data15_n                                   0
 #define fcs_front_xq_aI_data15_byte                                0
/* #define fcs_front_xq_aI_data15_ptr(e)*/

 #define fcs_front_xq_aI_data15_byte_flex                           0
 #define fcs_front_xq_aI_data15_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data15_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data15_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data15_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data15_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data15_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data15_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data15_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data15_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data15_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data15_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data15_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data15_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data15_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data15_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data15_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data15_assign(src, dst)
 #define fcs_front_xq_aI_cc_data15_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data15_null(e)
 #define fcs_front_xq_aI_cc_data15_inc(e)
 #define fcs_front_xq_aI_cc_data15_dec(e)
 #define fcs_front_xq_aI_cc_data15_add(e, n)
 #define fcs_front_xq_aI_cc_data15_sub(e, n)
 #define fcs_front_xq_aI_cc_data15_copy(src, dst)
 #define fcs_front_xq_aI_cc_data15_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data15_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data15_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data15_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data15_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data15_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data15_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA15 */

#define fcs_front_xq_aI_data15_cm                                   SLCM_DATA15  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA16 fcs_front_xq_aI_SL_DATA16_IGNORE fcs_front_xq_aI_sl_data16_type_c fcs_front_xq_aI_sl_data16_size_c fcs_front_xq_aI_sl_data16_type_mpi fcs_front_xq_aI_sl_data16_size_mpi fcs_front_xq_aI_sl_data16_memcpy fcs_front_xq_aI_sl_data16_weight fcs_front_xq_aI_sl_data16_flex */

#ifdef fcs_front_xq_aI_SL_DATA16

 #define fcs_front_xq_aI_sl_data16_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data16_size_c) * sizeof(fcs_front_xq_aI_sl_data16_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data16_copy
  #if fcs_front_xq_aI_sl_data16_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data16_memcpy)
   #if fcs_front_xq_aI_sl_data16_size_c == 1
    #define fcs_front_xq_aI_sl_data16_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data16_size_c == 2
    #define fcs_front_xq_aI_sl_data16_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data16_size_c == 3
    #define fcs_front_xq_aI_sl_data16_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data16_size_c == 4
    #define fcs_front_xq_aI_sl_data16_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data16_size_c == 5
    #define fcs_front_xq_aI_sl_data16_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data16_size_c == 6
    #define fcs_front_xq_aI_sl_data16_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data16_size_c == 7
    #define fcs_front_xq_aI_sl_data16_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data16_size_c == 8
    #define fcs_front_xq_aI_sl_data16_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data16_size_c == 9
    #define fcs_front_xq_aI_sl_data16_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data16_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data16_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data16_ncopy
  #define fcs_front_xq_aI_sl_data16_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data16_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data16_nmove
  #define fcs_front_xq_aI_sl_data16_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data16_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data16_type_c                              fcs_front_xq_aI_sl_data16_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data16_size_c                              (fcs_front_xq_aI_sl_data16_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data16_type_mpi                            (fcs_front_xq_aI_sl_data16_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data16_size_mpi                            (fcs_front_xq_aI_sl_data16_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data16_idx                                 16  /* sl_macro */

 #define fcs_front_xq_aI_data16_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data16_byte                                (fcs_front_xq_aI_sl_data16_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data16_ptr(e)                              (e)->data16  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data16_flex
 # define fcs_front_xq_aI_data16_byte_flex                          (fcs_front_xq_aI_sl_data16_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data16_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data16_weight
 # define fcs_front_xq_aI_data16_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data16_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data16_assign(src, dst)                    ((dst)->data16 = (src)->data16)  /* sl_macro */
 #define fcs_front_xq_aI_data16_assign_at(src, sat, dst)            ((dst)->data16 = &(src)->data16[(sat) * fcs_front_xq_aI_data16_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data16_null(e)                             ((e)->data16 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data16_inc(e)                              ((e)->data16 += fcs_front_xq_aI_data16_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data16_dec(e)                              ((e)->data16 -= fcs_front_xq_aI_data16_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data16_add(e, n)                           ((e)->data16 += (n) * fcs_front_xq_aI_data16_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data16_sub(e, n)                           ((e)->data16 -= (n) * fcs_front_xq_aI_data16_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data16_copy(src, dst)                      fcs_front_xq_aI_sl_data16_copy((src)->data16, (dst)->data16)  /* sl_macro */
 #define fcs_front_xq_aI_data16_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data16_ncopy((src)->data16, (dst)->data16, n)  /* sl_macro */
 #define fcs_front_xq_aI_data16_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data16_nmove((src)->data16, (dst)->data16, n)  /* sl_macro */

 #define fcs_front_xq_aI_data16_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data16_copy(&(src)->data16[(sat) * fcs_front_xq_aI_data16_size_c], &(dst)->data16[(dat) * fcs_front_xq_aI_data16_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data16_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data16_ncopy(&(src)->data16[(sat) * fcs_front_xq_aI_data16_size_c], &(dst)->data16[(dat) * fcs_front_xq_aI_data16_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data16_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data16_nmove(&(src)->data16[(sat) * fcs_front_xq_aI_data16_size_c], &(dst)->data16[(dat) * fcs_front_xq_aI_data16_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data16_xchange(e0, e1, t)                  (fcs_front_xq_aI_data16_copy(e0, t), fcs_front_xq_aI_data16_copy(e1, e0), fcs_front_xq_aI_data16_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data16_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data16_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data16_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data16_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data16_assign(src, dst)                 , fcs_front_xq_aI_data16_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_assign_at(src, sat, dst)         , fcs_front_xq_aI_data16_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_null(e)                          , fcs_front_xq_aI_data16_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_inc(e)                           , fcs_front_xq_aI_data16_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_dec(e)                           , fcs_front_xq_aI_data16_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_add(e, n)                        , fcs_front_xq_aI_data16_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_sub(e, n)                        , fcs_front_xq_aI_data16_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_copy(src, dst)                   , fcs_front_xq_aI_data16_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_ncopy(src, dst, n)               , fcs_front_xq_aI_data16_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_nmove(src, dst, n)               , fcs_front_xq_aI_data16_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data16_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data16_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data16_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_xchange(e0, e1, t)               , fcs_front_xq_aI_data16_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data16_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data16_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data16_assign(src, dst)                 fcs_front_xq_aI_data16_assign(src, dst)
 #define fcs_front_xq_aI_cc_data16_assign_at(src, sat, dst)         fcs_front_xq_aI_data16_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data16_null(e)                          fcs_front_xq_aI_data16_null(e)
 #define fcs_front_xq_aI_cc_data16_inc(e)                           fcs_front_xq_aI_data16_inc(e)
 #define fcs_front_xq_aI_cc_data16_dec(e)                           fcs_front_xq_aI_data16_dec(e)
 #define fcs_front_xq_aI_cc_data16_add(e, n)                        fcs_front_xq_aI_data16_add(e, n)
 #define fcs_front_xq_aI_cc_data16_sub(e, n)                        fcs_front_xq_aI_data16_sub(e, n)
 #define fcs_front_xq_aI_cc_data16_copy(src, dst)                   fcs_front_xq_aI_data16_copy(src, dst)
 #define fcs_front_xq_aI_cc_data16_ncopy(src, dst, n)               fcs_front_xq_aI_data16_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data16_nmove(src, dst, n)               fcs_front_xq_aI_data16_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data16_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data16_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data16_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data16_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data16_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data16_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data16_xchange(e0, e1, t)               fcs_front_xq_aI_data16_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data16_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data16_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA16 */

 #define fcs_front_xq_aI_data16_n                                   0
 #define fcs_front_xq_aI_data16_byte                                0
/* #define fcs_front_xq_aI_data16_ptr(e)*/

 #define fcs_front_xq_aI_data16_byte_flex                           0
 #define fcs_front_xq_aI_data16_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data16_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data16_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data16_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data16_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data16_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data16_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data16_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data16_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data16_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data16_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data16_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data16_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data16_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data16_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data16_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data16_assign(src, dst)
 #define fcs_front_xq_aI_cc_data16_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data16_null(e)
 #define fcs_front_xq_aI_cc_data16_inc(e)
 #define fcs_front_xq_aI_cc_data16_dec(e)
 #define fcs_front_xq_aI_cc_data16_add(e, n)
 #define fcs_front_xq_aI_cc_data16_sub(e, n)
 #define fcs_front_xq_aI_cc_data16_copy(src, dst)
 #define fcs_front_xq_aI_cc_data16_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data16_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data16_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data16_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data16_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data16_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data16_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA16 */

#define fcs_front_xq_aI_data16_cm                                   SLCM_DATA16  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA17 fcs_front_xq_aI_SL_DATA17_IGNORE fcs_front_xq_aI_sl_data17_type_c fcs_front_xq_aI_sl_data17_size_c fcs_front_xq_aI_sl_data17_type_mpi fcs_front_xq_aI_sl_data17_size_mpi fcs_front_xq_aI_sl_data17_memcpy fcs_front_xq_aI_sl_data17_weight fcs_front_xq_aI_sl_data17_flex */

#ifdef fcs_front_xq_aI_SL_DATA17

 #define fcs_front_xq_aI_sl_data17_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data17_size_c) * sizeof(fcs_front_xq_aI_sl_data17_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data17_copy
  #if fcs_front_xq_aI_sl_data17_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data17_memcpy)
   #if fcs_front_xq_aI_sl_data17_size_c == 1
    #define fcs_front_xq_aI_sl_data17_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data17_size_c == 2
    #define fcs_front_xq_aI_sl_data17_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data17_size_c == 3
    #define fcs_front_xq_aI_sl_data17_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data17_size_c == 4
    #define fcs_front_xq_aI_sl_data17_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data17_size_c == 5
    #define fcs_front_xq_aI_sl_data17_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data17_size_c == 6
    #define fcs_front_xq_aI_sl_data17_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data17_size_c == 7
    #define fcs_front_xq_aI_sl_data17_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data17_size_c == 8
    #define fcs_front_xq_aI_sl_data17_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data17_size_c == 9
    #define fcs_front_xq_aI_sl_data17_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data17_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data17_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data17_ncopy
  #define fcs_front_xq_aI_sl_data17_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data17_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data17_nmove
  #define fcs_front_xq_aI_sl_data17_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data17_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data17_type_c                              fcs_front_xq_aI_sl_data17_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data17_size_c                              (fcs_front_xq_aI_sl_data17_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data17_type_mpi                            (fcs_front_xq_aI_sl_data17_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data17_size_mpi                            (fcs_front_xq_aI_sl_data17_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data17_idx                                 17  /* sl_macro */

 #define fcs_front_xq_aI_data17_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data17_byte                                (fcs_front_xq_aI_sl_data17_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data17_ptr(e)                              (e)->data17  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data17_flex
 # define fcs_front_xq_aI_data17_byte_flex                          (fcs_front_xq_aI_sl_data17_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data17_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data17_weight
 # define fcs_front_xq_aI_data17_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data17_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data17_assign(src, dst)                    ((dst)->data17 = (src)->data17)  /* sl_macro */
 #define fcs_front_xq_aI_data17_assign_at(src, sat, dst)            ((dst)->data17 = &(src)->data17[(sat) * fcs_front_xq_aI_data17_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data17_null(e)                             ((e)->data17 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data17_inc(e)                              ((e)->data17 += fcs_front_xq_aI_data17_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data17_dec(e)                              ((e)->data17 -= fcs_front_xq_aI_data17_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data17_add(e, n)                           ((e)->data17 += (n) * fcs_front_xq_aI_data17_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data17_sub(e, n)                           ((e)->data17 -= (n) * fcs_front_xq_aI_data17_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data17_copy(src, dst)                      fcs_front_xq_aI_sl_data17_copy((src)->data17, (dst)->data17)  /* sl_macro */
 #define fcs_front_xq_aI_data17_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data17_ncopy((src)->data17, (dst)->data17, n)  /* sl_macro */
 #define fcs_front_xq_aI_data17_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data17_nmove((src)->data17, (dst)->data17, n)  /* sl_macro */

 #define fcs_front_xq_aI_data17_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data17_copy(&(src)->data17[(sat) * fcs_front_xq_aI_data17_size_c], &(dst)->data17[(dat) * fcs_front_xq_aI_data17_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data17_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data17_ncopy(&(src)->data17[(sat) * fcs_front_xq_aI_data17_size_c], &(dst)->data17[(dat) * fcs_front_xq_aI_data17_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data17_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data17_nmove(&(src)->data17[(sat) * fcs_front_xq_aI_data17_size_c], &(dst)->data17[(dat) * fcs_front_xq_aI_data17_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data17_xchange(e0, e1, t)                  (fcs_front_xq_aI_data17_copy(e0, t), fcs_front_xq_aI_data17_copy(e1, e0), fcs_front_xq_aI_data17_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data17_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data17_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data17_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data17_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data17_assign(src, dst)                 , fcs_front_xq_aI_data17_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_assign_at(src, sat, dst)         , fcs_front_xq_aI_data17_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_null(e)                          , fcs_front_xq_aI_data17_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_inc(e)                           , fcs_front_xq_aI_data17_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_dec(e)                           , fcs_front_xq_aI_data17_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_add(e, n)                        , fcs_front_xq_aI_data17_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_sub(e, n)                        , fcs_front_xq_aI_data17_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_copy(src, dst)                   , fcs_front_xq_aI_data17_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_ncopy(src, dst, n)               , fcs_front_xq_aI_data17_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_nmove(src, dst, n)               , fcs_front_xq_aI_data17_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data17_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data17_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data17_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_xchange(e0, e1, t)               , fcs_front_xq_aI_data17_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data17_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data17_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data17_assign(src, dst)                 fcs_front_xq_aI_data17_assign(src, dst)
 #define fcs_front_xq_aI_cc_data17_assign_at(src, sat, dst)         fcs_front_xq_aI_data17_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data17_null(e)                          fcs_front_xq_aI_data17_null(e)
 #define fcs_front_xq_aI_cc_data17_inc(e)                           fcs_front_xq_aI_data17_inc(e)
 #define fcs_front_xq_aI_cc_data17_dec(e)                           fcs_front_xq_aI_data17_dec(e)
 #define fcs_front_xq_aI_cc_data17_add(e, n)                        fcs_front_xq_aI_data17_add(e, n)
 #define fcs_front_xq_aI_cc_data17_sub(e, n)                        fcs_front_xq_aI_data17_sub(e, n)
 #define fcs_front_xq_aI_cc_data17_copy(src, dst)                   fcs_front_xq_aI_data17_copy(src, dst)
 #define fcs_front_xq_aI_cc_data17_ncopy(src, dst, n)               fcs_front_xq_aI_data17_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data17_nmove(src, dst, n)               fcs_front_xq_aI_data17_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data17_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data17_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data17_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data17_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data17_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data17_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data17_xchange(e0, e1, t)               fcs_front_xq_aI_data17_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data17_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data17_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA17 */

 #define fcs_front_xq_aI_data17_n                                   0
 #define fcs_front_xq_aI_data17_byte                                0
/* #define fcs_front_xq_aI_data17_ptr(e)*/

 #define fcs_front_xq_aI_data17_byte_flex                           0
 #define fcs_front_xq_aI_data17_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data17_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data17_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data17_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data17_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data17_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data17_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data17_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data17_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data17_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data17_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data17_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data17_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data17_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data17_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data17_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data17_assign(src, dst)
 #define fcs_front_xq_aI_cc_data17_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data17_null(e)
 #define fcs_front_xq_aI_cc_data17_inc(e)
 #define fcs_front_xq_aI_cc_data17_dec(e)
 #define fcs_front_xq_aI_cc_data17_add(e, n)
 #define fcs_front_xq_aI_cc_data17_sub(e, n)
 #define fcs_front_xq_aI_cc_data17_copy(src, dst)
 #define fcs_front_xq_aI_cc_data17_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data17_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data17_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data17_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data17_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data17_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data17_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA17 */

#define fcs_front_xq_aI_data17_cm                                   SLCM_DATA17  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA18 fcs_front_xq_aI_SL_DATA18_IGNORE fcs_front_xq_aI_sl_data18_type_c fcs_front_xq_aI_sl_data18_size_c fcs_front_xq_aI_sl_data18_type_mpi fcs_front_xq_aI_sl_data18_size_mpi fcs_front_xq_aI_sl_data18_memcpy fcs_front_xq_aI_sl_data18_weight fcs_front_xq_aI_sl_data18_flex */

#ifdef fcs_front_xq_aI_SL_DATA18

 #define fcs_front_xq_aI_sl_data18_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data18_size_c) * sizeof(fcs_front_xq_aI_sl_data18_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data18_copy
  #if fcs_front_xq_aI_sl_data18_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data18_memcpy)
   #if fcs_front_xq_aI_sl_data18_size_c == 1
    #define fcs_front_xq_aI_sl_data18_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data18_size_c == 2
    #define fcs_front_xq_aI_sl_data18_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data18_size_c == 3
    #define fcs_front_xq_aI_sl_data18_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data18_size_c == 4
    #define fcs_front_xq_aI_sl_data18_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data18_size_c == 5
    #define fcs_front_xq_aI_sl_data18_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data18_size_c == 6
    #define fcs_front_xq_aI_sl_data18_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data18_size_c == 7
    #define fcs_front_xq_aI_sl_data18_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data18_size_c == 8
    #define fcs_front_xq_aI_sl_data18_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data18_size_c == 9
    #define fcs_front_xq_aI_sl_data18_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data18_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data18_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data18_ncopy
  #define fcs_front_xq_aI_sl_data18_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data18_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data18_nmove
  #define fcs_front_xq_aI_sl_data18_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data18_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data18_type_c                              fcs_front_xq_aI_sl_data18_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data18_size_c                              (fcs_front_xq_aI_sl_data18_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data18_type_mpi                            (fcs_front_xq_aI_sl_data18_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data18_size_mpi                            (fcs_front_xq_aI_sl_data18_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data18_idx                                 18  /* sl_macro */

 #define fcs_front_xq_aI_data18_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data18_byte                                (fcs_front_xq_aI_sl_data18_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data18_ptr(e)                              (e)->data18  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data18_flex
 # define fcs_front_xq_aI_data18_byte_flex                          (fcs_front_xq_aI_sl_data18_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data18_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data18_weight
 # define fcs_front_xq_aI_data18_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data18_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data18_assign(src, dst)                    ((dst)->data18 = (src)->data18)  /* sl_macro */
 #define fcs_front_xq_aI_data18_assign_at(src, sat, dst)            ((dst)->data18 = &(src)->data18[(sat) * fcs_front_xq_aI_data18_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data18_null(e)                             ((e)->data18 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data18_inc(e)                              ((e)->data18 += fcs_front_xq_aI_data18_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data18_dec(e)                              ((e)->data18 -= fcs_front_xq_aI_data18_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data18_add(e, n)                           ((e)->data18 += (n) * fcs_front_xq_aI_data18_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data18_sub(e, n)                           ((e)->data18 -= (n) * fcs_front_xq_aI_data18_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data18_copy(src, dst)                      fcs_front_xq_aI_sl_data18_copy((src)->data18, (dst)->data18)  /* sl_macro */
 #define fcs_front_xq_aI_data18_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data18_ncopy((src)->data18, (dst)->data18, n)  /* sl_macro */
 #define fcs_front_xq_aI_data18_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data18_nmove((src)->data18, (dst)->data18, n)  /* sl_macro */

 #define fcs_front_xq_aI_data18_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data18_copy(&(src)->data18[(sat) * fcs_front_xq_aI_data18_size_c], &(dst)->data18[(dat) * fcs_front_xq_aI_data18_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data18_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data18_ncopy(&(src)->data18[(sat) * fcs_front_xq_aI_data18_size_c], &(dst)->data18[(dat) * fcs_front_xq_aI_data18_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data18_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data18_nmove(&(src)->data18[(sat) * fcs_front_xq_aI_data18_size_c], &(dst)->data18[(dat) * fcs_front_xq_aI_data18_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data18_xchange(e0, e1, t)                  (fcs_front_xq_aI_data18_copy(e0, t), fcs_front_xq_aI_data18_copy(e1, e0), fcs_front_xq_aI_data18_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data18_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data18_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data18_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data18_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data18_assign(src, dst)                 , fcs_front_xq_aI_data18_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_assign_at(src, sat, dst)         , fcs_front_xq_aI_data18_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_null(e)                          , fcs_front_xq_aI_data18_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_inc(e)                           , fcs_front_xq_aI_data18_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_dec(e)                           , fcs_front_xq_aI_data18_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_add(e, n)                        , fcs_front_xq_aI_data18_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_sub(e, n)                        , fcs_front_xq_aI_data18_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_copy(src, dst)                   , fcs_front_xq_aI_data18_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_ncopy(src, dst, n)               , fcs_front_xq_aI_data18_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_nmove(src, dst, n)               , fcs_front_xq_aI_data18_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data18_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data18_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data18_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_xchange(e0, e1, t)               , fcs_front_xq_aI_data18_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data18_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data18_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data18_assign(src, dst)                 fcs_front_xq_aI_data18_assign(src, dst)
 #define fcs_front_xq_aI_cc_data18_assign_at(src, sat, dst)         fcs_front_xq_aI_data18_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data18_null(e)                          fcs_front_xq_aI_data18_null(e)
 #define fcs_front_xq_aI_cc_data18_inc(e)                           fcs_front_xq_aI_data18_inc(e)
 #define fcs_front_xq_aI_cc_data18_dec(e)                           fcs_front_xq_aI_data18_dec(e)
 #define fcs_front_xq_aI_cc_data18_add(e, n)                        fcs_front_xq_aI_data18_add(e, n)
 #define fcs_front_xq_aI_cc_data18_sub(e, n)                        fcs_front_xq_aI_data18_sub(e, n)
 #define fcs_front_xq_aI_cc_data18_copy(src, dst)                   fcs_front_xq_aI_data18_copy(src, dst)
 #define fcs_front_xq_aI_cc_data18_ncopy(src, dst, n)               fcs_front_xq_aI_data18_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data18_nmove(src, dst, n)               fcs_front_xq_aI_data18_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data18_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data18_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data18_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data18_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data18_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data18_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data18_xchange(e0, e1, t)               fcs_front_xq_aI_data18_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data18_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data18_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA18 */

 #define fcs_front_xq_aI_data18_n                                   0
 #define fcs_front_xq_aI_data18_byte                                0
/* #define fcs_front_xq_aI_data18_ptr(e)*/

 #define fcs_front_xq_aI_data18_byte_flex                           0
 #define fcs_front_xq_aI_data18_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data18_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data18_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data18_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data18_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data18_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data18_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data18_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data18_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data18_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data18_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data18_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data18_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data18_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data18_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data18_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data18_assign(src, dst)
 #define fcs_front_xq_aI_cc_data18_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data18_null(e)
 #define fcs_front_xq_aI_cc_data18_inc(e)
 #define fcs_front_xq_aI_cc_data18_dec(e)
 #define fcs_front_xq_aI_cc_data18_add(e, n)
 #define fcs_front_xq_aI_cc_data18_sub(e, n)
 #define fcs_front_xq_aI_cc_data18_copy(src, dst)
 #define fcs_front_xq_aI_cc_data18_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data18_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data18_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data18_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data18_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data18_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data18_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA18 */

#define fcs_front_xq_aI_data18_cm                                   SLCM_DATA18  /* sl_macro */


/* sl_macro fcs_front_xq_aI_SL_DATA19 fcs_front_xq_aI_SL_DATA19_IGNORE fcs_front_xq_aI_sl_data19_type_c fcs_front_xq_aI_sl_data19_size_c fcs_front_xq_aI_sl_data19_type_mpi fcs_front_xq_aI_sl_data19_size_mpi fcs_front_xq_aI_sl_data19_memcpy fcs_front_xq_aI_sl_data19_weight fcs_front_xq_aI_sl_data19_flex */

#ifdef fcs_front_xq_aI_SL_DATA19

 #define fcs_front_xq_aI_sl_data19_byte                             ((fcs_front_xq_aI_slint_t) (fcs_front_xq_aI_sl_data19_size_c) * sizeof(fcs_front_xq_aI_sl_data19_type_c))  /* sl_macro */

 #ifndef fcs_front_xq_aI_sl_data19_copy
  #if fcs_front_xq_aI_sl_data19_size_c <= 9 && !defined(fcs_front_xq_aI_sl_data19_memcpy)
   #if fcs_front_xq_aI_sl_data19_size_c == 1
    #define fcs_front_xq_aI_sl_data19_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data19_size_c == 2
    #define fcs_front_xq_aI_sl_data19_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data19_size_c == 3
    #define fcs_front_xq_aI_sl_data19_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data19_size_c == 4
    #define fcs_front_xq_aI_sl_data19_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data19_size_c == 5
    #define fcs_front_xq_aI_sl_data19_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data19_size_c == 6
    #define fcs_front_xq_aI_sl_data19_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data19_size_c == 7
    #define fcs_front_xq_aI_sl_data19_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data19_size_c == 8
    #define fcs_front_xq_aI_sl_data19_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif fcs_front_xq_aI_sl_data19_size_c == 9
    #define fcs_front_xq_aI_sl_data19_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define fcs_front_xq_aI_sl_data19_copy(src, dst)                 memcpy(dst, src, fcs_front_xq_aI_sl_data19_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef fcs_front_xq_aI_sl_data19_ncopy
  #define fcs_front_xq_aI_sl_data19_ncopy(src, dst, n)              memcpy(dst, src, (n) * fcs_front_xq_aI_sl_data19_byte)  /* sl_macro */
 #endif
 #ifndef fcs_front_xq_aI_sl_data19_nmove
  #define fcs_front_xq_aI_sl_data19_nmove(src, dst, n)              memmove(dst, src, (n) * fcs_front_xq_aI_sl_data19_byte)  /* sl_macro */
 #endif

 #define fcs_front_xq_aI_data19_type_c                              fcs_front_xq_aI_sl_data19_type_c  /* sl_macro */
 #define fcs_front_xq_aI_data19_size_c                              (fcs_front_xq_aI_sl_data19_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data19_type_mpi                            (fcs_front_xq_aI_sl_data19_type_mpi)  /* sl_macro */
 #define fcs_front_xq_aI_data19_size_mpi                            (fcs_front_xq_aI_sl_data19_size_mpi)  /* sl_macro */

 #define fcs_front_xq_aI_data19_idx                                 19  /* sl_macro */

 #define fcs_front_xq_aI_data19_n                                   1  /* sl_macro */
 #define fcs_front_xq_aI_data19_byte                                (fcs_front_xq_aI_sl_data19_byte)  /* sl_macro */
 #define fcs_front_xq_aI_data19_ptr(e)                              (e)->data19  /* sl_macro */

 #ifdef fcs_front_xq_aI_sl_data19_flex
 # define fcs_front_xq_aI_data19_byte_flex                          (fcs_front_xq_aI_sl_data19_byte)  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data19_byte_flex                          0
 #endif

 #ifdef fcs_front_xq_aI_sl_data19_weight
 # define fcs_front_xq_aI_data19_weight                             1  /* sl_macro */
 #else
 # define fcs_front_xq_aI_data19_weight                             0
 #endif

 /* commands for regular use */
 #define fcs_front_xq_aI_data19_assign(src, dst)                    ((dst)->data19 = (src)->data19)  /* sl_macro */
 #define fcs_front_xq_aI_data19_assign_at(src, sat, dst)            ((dst)->data19 = &(src)->data19[(sat) * fcs_front_xq_aI_data19_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data19_null(e)                             ((e)->data19 = NULL)  /* sl_macro */
 #define fcs_front_xq_aI_data19_inc(e)                              ((e)->data19 += fcs_front_xq_aI_data19_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data19_dec(e)                              ((e)->data19 -= fcs_front_xq_aI_data19_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data19_add(e, n)                           ((e)->data19 += (n) * fcs_front_xq_aI_data19_size_c)  /* sl_macro */
 #define fcs_front_xq_aI_data19_sub(e, n)                           ((e)->data19 -= (n) * fcs_front_xq_aI_data19_size_c)  /* sl_macro */

 #define fcs_front_xq_aI_data19_copy(src, dst)                      fcs_front_xq_aI_sl_data19_copy((src)->data19, (dst)->data19)  /* sl_macro */
 #define fcs_front_xq_aI_data19_ncopy(src, dst, n)                  fcs_front_xq_aI_sl_data19_ncopy((src)->data19, (dst)->data19, n)  /* sl_macro */
 #define fcs_front_xq_aI_data19_nmove(src, dst, n)                  fcs_front_xq_aI_sl_data19_nmove((src)->data19, (dst)->data19, n)  /* sl_macro */

 #define fcs_front_xq_aI_data19_copy_at(src, sat, dst, dat)         fcs_front_xq_aI_sl_data19_copy(&(src)->data19[(sat) * fcs_front_xq_aI_data19_size_c], &(dst)->data19[(dat) * fcs_front_xq_aI_data19_size_c])  /* sl_macro */
 #define fcs_front_xq_aI_data19_ncopy_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data19_ncopy(&(src)->data19[(sat) * fcs_front_xq_aI_data19_size_c], &(dst)->data19[(dat) * fcs_front_xq_aI_data19_size_c], n)  /* sl_macro */
 #define fcs_front_xq_aI_data19_nmove_at(src, sat, dst, dat, n)     fcs_front_xq_aI_sl_data19_nmove(&(src)->data19[(sat) * fcs_front_xq_aI_data19_size_c], &(dst)->data19[(dat) * fcs_front_xq_aI_data19_size_c], n)  /* sl_macro */

 #define fcs_front_xq_aI_data19_xchange(e0, e1, t)                  (fcs_front_xq_aI_data19_copy(e0, t), fcs_front_xq_aI_data19_copy(e1, e0), fcs_front_xq_aI_data19_copy(t, e1))  /* sl_macro */
 #define fcs_front_xq_aI_data19_xchange_at(e0, at0, e1, at1, t)     (fcs_front_xq_aI_data19_copy_at(e0, at0, t, 0), fcs_front_xq_aI_data19_copy_at(e1, at1, e0, at0), fcs_front_xq_aI_data19_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define fcs_front_xq_aI_cc_data19_assign(src, dst)                 , fcs_front_xq_aI_data19_assign(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_assign_at(src, sat, dst)         , fcs_front_xq_aI_data19_assign_at(src, sat, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_null(e)                          , fcs_front_xq_aI_data19_null(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_inc(e)                           , fcs_front_xq_aI_data19_inc(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_dec(e)                           , fcs_front_xq_aI_data19_dec(e)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_add(e, n)                        , fcs_front_xq_aI_data19_add(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_sub(e, n)                        , fcs_front_xq_aI_data19_sub(e, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_copy(src, dst)                   , fcs_front_xq_aI_data19_copy(src, dst)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_ncopy(src, dst, n)               , fcs_front_xq_aI_data19_ncopy(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_nmove(src, dst, n)               , fcs_front_xq_aI_data19_nmove(src, dst, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_copy_at(src, sat, dst, dat)      , fcs_front_xq_aI_data19_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_ncopy_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data19_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_nmove_at(src, sat, dst, dat, n)  , fcs_front_xq_aI_data19_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_xchange(e0, e1, t)               , fcs_front_xq_aI_data19_xchange(e0, e1, t)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data19_xchange_at(e0, at0, e1, at1, t)  , fcs_front_xq_aI_data19_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define fcs_front_xq_aI_cc_data19_assign(src, dst)                 fcs_front_xq_aI_data19_assign(src, dst)
 #define fcs_front_xq_aI_cc_data19_assign_at(src, sat, dst)         fcs_front_xq_aI_data19_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data19_null(e)                          fcs_front_xq_aI_data19_null(e)
 #define fcs_front_xq_aI_cc_data19_inc(e)                           fcs_front_xq_aI_data19_inc(e)
 #define fcs_front_xq_aI_cc_data19_dec(e)                           fcs_front_xq_aI_data19_dec(e)
 #define fcs_front_xq_aI_cc_data19_add(e, n)                        fcs_front_xq_aI_data19_add(e, n)
 #define fcs_front_xq_aI_cc_data19_sub(e, n)                        fcs_front_xq_aI_data19_sub(e, n)
 #define fcs_front_xq_aI_cc_data19_copy(src, dst)                   fcs_front_xq_aI_data19_copy(src, dst)
 #define fcs_front_xq_aI_cc_data19_ncopy(src, dst, n)               fcs_front_xq_aI_data19_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data19_nmove(src, dst, n)               fcs_front_xq_aI_data19_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data19_copy_at(src, sat, dst, dat)      fcs_front_xq_aI_data19_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data19_ncopy_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data19_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data19_nmove_at(src, sat, dst, dat, n)  fcs_front_xq_aI_data19_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data19_xchange(e0, e1, t)               fcs_front_xq_aI_data19_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data19_xchange_at(e0, at0, e1, at1, t)  fcs_front_xq_aI_data19_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* fcs_front_xq_aI_SL_DATA19 */

 #define fcs_front_xq_aI_data19_n                                   0
 #define fcs_front_xq_aI_data19_byte                                0
/* #define fcs_front_xq_aI_data19_ptr(e)*/

 #define fcs_front_xq_aI_data19_byte_flex                           0
 #define fcs_front_xq_aI_data19_weight                              0

 /* commands for regular use */
 #define fcs_front_xq_aI_data19_assign(src, dst)                    Z_NOP()
 #define fcs_front_xq_aI_data19_assign_at(src, sat, dst)            Z_NOP()
 #define fcs_front_xq_aI_data19_null(e)                             Z_NOP()
 #define fcs_front_xq_aI_data19_inc(e)                              Z_NOP()
 #define fcs_front_xq_aI_data19_dec(e)                              Z_NOP()
 #define fcs_front_xq_aI_data19_add(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data19_sub(e, n)                           Z_NOP()
 #define fcs_front_xq_aI_data19_copy(src, dst)                      Z_NOP()
 #define fcs_front_xq_aI_data19_ncopy(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data19_nmove(src, dst, n)                  Z_NOP()
 #define fcs_front_xq_aI_data19_copy_at(src, sat, dst, dat)         Z_NOP()
 #define fcs_front_xq_aI_data19_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data19_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define fcs_front_xq_aI_data19_xchange(e0, e1, t)                  Z_NOP()
 #define fcs_front_xq_aI_data19_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define fcs_front_xq_aI_cc_data19_assign(src, dst)
 #define fcs_front_xq_aI_cc_data19_assign_at(src, sat, dst)
 #define fcs_front_xq_aI_cc_data19_null(e)
 #define fcs_front_xq_aI_cc_data19_inc(e)
 #define fcs_front_xq_aI_cc_data19_dec(e)
 #define fcs_front_xq_aI_cc_data19_add(e, n)
 #define fcs_front_xq_aI_cc_data19_sub(e, n)
 #define fcs_front_xq_aI_cc_data19_copy(src, dst)
 #define fcs_front_xq_aI_cc_data19_ncopy(src, dst, n)
 #define fcs_front_xq_aI_cc_data19_nmove(src, dst, n)
 #define fcs_front_xq_aI_cc_data19_copy_at(src, sat, dst, dat)
 #define fcs_front_xq_aI_cc_data19_ncopy_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data19_nmove_at(src, sat, dst, dat, n)
 #define fcs_front_xq_aI_cc_data19_xchange(e0, e1, t)
 #define fcs_front_xq_aI_cc_data19_xchange_at(e0, at0, e1, at1, t)

#endif /* fcs_front_xq_aI_SL_DATA19 */

#define fcs_front_xq_aI_data19_cm                                   SLCM_DATA19  /* sl_macro */

/* DATAX_TEMPLATE_END */







/* sl_macro fcs_front_xq_aI_data_nmax */
#define fcs_front_xq_aI_data_nmax (0 \
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

/* sl_macro fcs_front_xq_aI_data_n */
#define fcs_front_xq_aI_data_n (0 \
 + (fcs_front_xq_aI_data0_n) \
 + (fcs_front_xq_aI_data1_n) \
 + (fcs_front_xq_aI_data2_n) \
 + (fcs_front_xq_aI_data3_n) \
 + (fcs_front_xq_aI_data4_n) \
 + (fcs_front_xq_aI_data5_n) \
 + (fcs_front_xq_aI_data6_n) \
 + (fcs_front_xq_aI_data7_n) \
 + (fcs_front_xq_aI_data8_n) \
 + (fcs_front_xq_aI_data9_n) \
 + (fcs_front_xq_aI_data10_n) \
 + (fcs_front_xq_aI_data11_n) \
 + (fcs_front_xq_aI_data12_n) \
 + (fcs_front_xq_aI_data13_n) \
 + (fcs_front_xq_aI_data14_n) \
 + (fcs_front_xq_aI_data15_n) \
 + (fcs_front_xq_aI_data16_n) \
 + (fcs_front_xq_aI_data17_n) \
 + (fcs_front_xq_aI_data18_n) \
 + (fcs_front_xq_aI_data19_n) \
)

/* sl_macro fcs_front_xq_aI_data_byte */
#define fcs_front_xq_aI_data_byte (0 \
 + (fcs_front_xq_aI_data0_byte) \
 + (fcs_front_xq_aI_data1_byte) \
 + (fcs_front_xq_aI_data2_byte) \
 + (fcs_front_xq_aI_data3_byte) \
 + (fcs_front_xq_aI_data4_byte) \
 + (fcs_front_xq_aI_data5_byte) \
 + (fcs_front_xq_aI_data6_byte) \
 + (fcs_front_xq_aI_data7_byte) \
 + (fcs_front_xq_aI_data8_byte) \
 + (fcs_front_xq_aI_data9_byte) \
 + (fcs_front_xq_aI_data10_byte) \
 + (fcs_front_xq_aI_data11_byte) \
 + (fcs_front_xq_aI_data12_byte) \
 + (fcs_front_xq_aI_data13_byte) \
 + (fcs_front_xq_aI_data14_byte) \
 + (fcs_front_xq_aI_data15_byte) \
 + (fcs_front_xq_aI_data16_byte) \
 + (fcs_front_xq_aI_data17_byte) \
 + (fcs_front_xq_aI_data18_byte) \
 + (fcs_front_xq_aI_data19_byte) \
)

/* sl_macro fcs_front_xq_aI_data_byte_flex */
#define fcs_front_xq_aI_data_byte_flex (0 \
 + (fcs_front_xq_aI_data0_byte_flex) \
 + (fcs_front_xq_aI_data1_byte_flex) \
 + (fcs_front_xq_aI_data2_byte_flex) \
 + (fcs_front_xq_aI_data3_byte_flex) \
 + (fcs_front_xq_aI_data4_byte_flex) \
 + (fcs_front_xq_aI_data5_byte_flex) \
 + (fcs_front_xq_aI_data6_byte_flex) \
 + (fcs_front_xq_aI_data7_byte_flex) \
 + (fcs_front_xq_aI_data8_byte_flex) \
 + (fcs_front_xq_aI_data9_byte_flex) \
 + (fcs_front_xq_aI_data10_byte_flex) \
 + (fcs_front_xq_aI_data11_byte_flex) \
 + (fcs_front_xq_aI_data12_byte_flex) \
 + (fcs_front_xq_aI_data13_byte_flex) \
 + (fcs_front_xq_aI_data14_byte_flex) \
 + (fcs_front_xq_aI_data15_byte_flex) \
 + (fcs_front_xq_aI_data16_byte_flex) \
 + (fcs_front_xq_aI_data17_byte_flex) \
 + (fcs_front_xq_aI_data18_byte_flex) \
 + (fcs_front_xq_aI_data19_byte_flex) \
)

/* sl_macro fcs_front_xq_aI_data_assign */
#define fcs_front_xq_aI_data_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data0_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data1_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data2_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data3_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data4_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data5_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data6_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data7_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data8_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data9_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data10_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data11_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data12_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data13_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data14_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data15_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data16_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data17_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data18_assign(_s_, _d_) \
 fcs_front_xq_aI_cc_data19_assign(_s_, _d_) \

/* sl_macro fcs_front_xq_aI_data_assign_at */
#define fcs_front_xq_aI_data_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data0_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data1_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data2_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data3_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data4_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data5_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data6_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data7_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data8_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data9_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data10_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data11_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data12_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data13_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data14_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data15_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data16_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data17_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data18_assign_at(_s_, _sat_, _d_) \
 fcs_front_xq_aI_cc_data19_assign_at(_s_, _sat_, _d_) \

/* sl_macro fcs_front_xq_aI_data_null */
#define fcs_front_xq_aI_data_null(_e_) \
 fcs_front_xq_aI_cc_data0_null(_e_) \
 fcs_front_xq_aI_cc_data1_null(_e_) \
 fcs_front_xq_aI_cc_data2_null(_e_) \
 fcs_front_xq_aI_cc_data3_null(_e_) \
 fcs_front_xq_aI_cc_data4_null(_e_) \
 fcs_front_xq_aI_cc_data5_null(_e_) \
 fcs_front_xq_aI_cc_data6_null(_e_) \
 fcs_front_xq_aI_cc_data7_null(_e_) \
 fcs_front_xq_aI_cc_data8_null(_e_) \
 fcs_front_xq_aI_cc_data9_null(_e_) \
 fcs_front_xq_aI_cc_data10_null(_e_) \
 fcs_front_xq_aI_cc_data11_null(_e_) \
 fcs_front_xq_aI_cc_data12_null(_e_) \
 fcs_front_xq_aI_cc_data13_null(_e_) \
 fcs_front_xq_aI_cc_data14_null(_e_) \
 fcs_front_xq_aI_cc_data15_null(_e_) \
 fcs_front_xq_aI_cc_data16_null(_e_) \
 fcs_front_xq_aI_cc_data17_null(_e_) \
 fcs_front_xq_aI_cc_data18_null(_e_) \
 fcs_front_xq_aI_cc_data19_null(_e_) \

/* sl_macro fcs_front_xq_aI_data_inc */
#define fcs_front_xq_aI_data_inc(_e_) \
 fcs_front_xq_aI_cc_data0_inc(_e_) \
 fcs_front_xq_aI_cc_data1_inc(_e_) \
 fcs_front_xq_aI_cc_data2_inc(_e_) \
 fcs_front_xq_aI_cc_data3_inc(_e_) \
 fcs_front_xq_aI_cc_data4_inc(_e_) \
 fcs_front_xq_aI_cc_data5_inc(_e_) \
 fcs_front_xq_aI_cc_data6_inc(_e_) \
 fcs_front_xq_aI_cc_data7_inc(_e_) \
 fcs_front_xq_aI_cc_data8_inc(_e_) \
 fcs_front_xq_aI_cc_data9_inc(_e_) \
 fcs_front_xq_aI_cc_data10_inc(_e_) \
 fcs_front_xq_aI_cc_data11_inc(_e_) \
 fcs_front_xq_aI_cc_data12_inc(_e_) \
 fcs_front_xq_aI_cc_data13_inc(_e_) \
 fcs_front_xq_aI_cc_data14_inc(_e_) \
 fcs_front_xq_aI_cc_data15_inc(_e_) \
 fcs_front_xq_aI_cc_data16_inc(_e_) \
 fcs_front_xq_aI_cc_data17_inc(_e_) \
 fcs_front_xq_aI_cc_data18_inc(_e_) \
 fcs_front_xq_aI_cc_data19_inc(_e_) \

/* sl_macro fcs_front_xq_aI_data_dec */
#define fcs_front_xq_aI_data_dec(_e_) \
 fcs_front_xq_aI_cc_data0_dec(_e_) \
 fcs_front_xq_aI_cc_data1_dec(_e_) \
 fcs_front_xq_aI_cc_data2_dec(_e_) \
 fcs_front_xq_aI_cc_data3_dec(_e_) \
 fcs_front_xq_aI_cc_data4_dec(_e_) \
 fcs_front_xq_aI_cc_data5_dec(_e_) \
 fcs_front_xq_aI_cc_data6_dec(_e_) \
 fcs_front_xq_aI_cc_data7_dec(_e_) \
 fcs_front_xq_aI_cc_data8_dec(_e_) \
 fcs_front_xq_aI_cc_data9_dec(_e_) \
 fcs_front_xq_aI_cc_data10_dec(_e_) \
 fcs_front_xq_aI_cc_data11_dec(_e_) \
 fcs_front_xq_aI_cc_data12_dec(_e_) \
 fcs_front_xq_aI_cc_data13_dec(_e_) \
 fcs_front_xq_aI_cc_data14_dec(_e_) \
 fcs_front_xq_aI_cc_data15_dec(_e_) \
 fcs_front_xq_aI_cc_data16_dec(_e_) \
 fcs_front_xq_aI_cc_data17_dec(_e_) \
 fcs_front_xq_aI_cc_data18_dec(_e_) \
 fcs_front_xq_aI_cc_data19_dec(_e_) \

/* sl_macro fcs_front_xq_aI_data_add */
#define fcs_front_xq_aI_data_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data0_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data1_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data2_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data3_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data4_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data5_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data6_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data7_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data8_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data9_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data10_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data11_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data12_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data13_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data14_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data15_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data16_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data17_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data18_add(_e_, _n_) \
 fcs_front_xq_aI_cc_data19_add(_e_, _n_) \

/* sl_macro fcs_front_xq_aI_data_sub */
#define fcs_front_xq_aI_data_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data0_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data1_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data2_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data3_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data4_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data5_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data6_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data7_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data8_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data9_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data10_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data11_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data12_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data13_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data14_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data15_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data16_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data17_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data18_sub(_e_, _n_) \
 fcs_front_xq_aI_cc_data19_sub(_e_, _n_) \

/* FIXME: add fcs_front_xq_aI_rti_cadd_moved(_n_,cmd) -> only ifdef SL_DATA else empty (like dataX) */

/* sl_macro fcs_front_xq_aI_data_copy */
#define fcs_front_xq_aI_data_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data0_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data1_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data2_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data3_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data4_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data5_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data6_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data7_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data8_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data9_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data10_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data11_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data12_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data13_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data14_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data15_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data16_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data17_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data18_copy(_s_, _d_) \
 fcs_front_xq_aI_cc_data19_copy(_s_, _d_) \

/* sl_macro fcs_front_xq_aI_data_ncopy */
#define fcs_front_xq_aI_data_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data0_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data1_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data2_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data3_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data4_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data5_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data6_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data7_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data8_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data9_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data10_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data11_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data12_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data13_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data14_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data15_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data16_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data17_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data18_ncopy(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data19_ncopy(_s_, _d_, _n_) \

/* sl_macro fcs_front_xq_aI_data_nmove */
#define fcs_front_xq_aI_data_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data0_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data1_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data2_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data3_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data4_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data5_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data6_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data7_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data8_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data9_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data10_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data11_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data12_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data13_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data14_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data15_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data16_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data17_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data18_nmove(_s_, _d_, _n_) \
 fcs_front_xq_aI_cc_data19_nmove(_s_, _d_, _n_) \

/* sl_macro fcs_front_xq_aI_data_copy_at */
#define fcs_front_xq_aI_data_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data0_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data1_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data2_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data3_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data4_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data5_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data6_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data7_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data8_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data9_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data10_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data11_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data12_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data13_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data14_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data15_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data16_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data17_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data18_copy_at(_s_, _sat_, _d_, _dat_) \
 fcs_front_xq_aI_cc_data19_copy_at(_s_, _sat_, _d_, _dat_) \

/* sl_macro fcs_front_xq_aI_data_ncopy_at */
#define fcs_front_xq_aI_data_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data0_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data1_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data2_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data3_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data4_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data5_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data6_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data7_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data8_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data9_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data10_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data11_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data12_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data13_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data14_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data15_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data16_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data17_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data18_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data19_ncopy_at(_s_, _sat_, _d_, _dat_, _n_) \

/* sl_macro fcs_front_xq_aI_data_nmove_at */
#define fcs_front_xq_aI_data_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data0_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data1_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data2_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data3_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data4_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data5_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data6_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data7_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data8_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data9_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data10_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data11_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data12_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data13_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data14_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data15_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data16_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data17_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data18_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \
 fcs_front_xq_aI_cc_data19_nmove_at(_s_, _sat_, _d_, _dat_, _n_) \

/* sl_macro fcs_front_xq_aI_data_xchange */
#define fcs_front_xq_aI_data_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data0_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data1_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data2_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data3_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data4_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data5_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data6_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data7_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data8_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data9_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data10_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data11_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data12_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data13_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data14_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data15_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data16_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data17_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data18_xchange(_e0_, _e1_, _t_) \
 fcs_front_xq_aI_cc_data19_xchange(_e0_, _e1_, _t_) \

/* sl_macro fcs_front_xq_aI_data_xchange_at */
#define fcs_front_xq_aI_data_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data0_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data1_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data2_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data3_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data4_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data5_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data6_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data7_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data8_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data9_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data10_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data11_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data12_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data13_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data14_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data15_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data16_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data17_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data18_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \
 fcs_front_xq_aI_cc_data19_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_) \

/* chained versions */
#ifdef SL_DATA

 #define fcs_front_xq_aI_cc_data_assign(_s_, _d_)                           , fcs_front_xq_aI_data_assign(_s_, _d_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_assign_at(_s_, _sat_, _d_)                 , fcs_front_xq_aI_data_assign_at(_s_, _sat_, _d_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_null(_e_)                                  , fcs_front_xq_aI_data_null(_e_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_inc(_e_)                                   , fcs_front_xq_aI_data_inc(_e_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_dec(_e_)                                   , fcs_front_xq_aI_data_dec(_e_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_add(_e_, _n_)                              , fcs_front_xq_aI_data_add(_e_, _n_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_sub(_e_, _n_)                              , fcs_front_xq_aI_data_sub(_e_, _n_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_copy(_s_, _d_)                             , fcs_front_xq_aI_data_copy(_s_, _d_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_ncopy(_s_, _d_, _n_)                       , fcs_front_xq_aI_data_ncopy(_s_, _d_, _n_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_nmove(_s_, _d_, _n_)                       , fcs_front_xq_aI_data_nmove(_s_, _d_, _n_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_copy_at(_s_, _sat_, _d_, _dat_)            , fcs_front_xq_aI_data_copy_at(_s_, _sat_, _d_, _dat_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      , fcs_front_xq_aI_data_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      , fcs_front_xq_aI_data_nmove_at(_s_, _sat_, _d_, _dat_, _n_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_xchange(_e0_, _e1_, _t_)                   , fcs_front_xq_aI_data_xchange(_e0_, _e1_, _t_)  /* sl_macro */
 #define fcs_front_xq_aI_cc_data_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  , fcs_front_xq_aI_data_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  /* sl_macro */

#else /* SL_DATA */

/* #define SL_DATA*/
 #define fcs_front_xq_aI_cc_data_assign(_s_, _d_)
 #define fcs_front_xq_aI_cc_data_assign_at(_s_, _sat_, _d_)
 #define fcs_front_xq_aI_cc_data_null(_e_)
 #define fcs_front_xq_aI_cc_data_inc(_e_)
 #define fcs_front_xq_aI_cc_data_dec(_e_)
 #define fcs_front_xq_aI_cc_data_add(_e_, _n_)
 #define fcs_front_xq_aI_cc_data_sub(_e_, _n_)
 #define fcs_front_xq_aI_cc_data_copy(_s_, _d_)
 #define fcs_front_xq_aI_cc_data_ncopy(_s_, _d_, _n_)
 #define fcs_front_xq_aI_cc_data_nmove(_s_, _d_, _n_)
 #define fcs_front_xq_aI_cc_data_copy_at(_s_, _sat_, _d_, _dat_)
 #define fcs_front_xq_aI_cc_data_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
 #define fcs_front_xq_aI_cc_data_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
 #define fcs_front_xq_aI_cc_data_xchange(_e0_, _e1_, _t_)
 #define fcs_front_xq_aI_cc_data_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)

#endif /* SL_DATA */



#undef SL_DATA


#define fcs_front_xq_aI_elem_n                                          (fcs_front_xq_aI_key_n + fcs_front_xq_aI_index_n + fcs_front_xq_aI_data_n)  /* sl_macro */
#define fcs_front_xq_aI_elem_byte                                       (fcs_front_xq_aI_key_byte + fcs_front_xq_aI_index_byte + fcs_front_xq_aI_data_byte)  /* sl_macro */

#define fcs_front_xq_aI_elem_key_at(_s_, _sat_)                         (fcs_front_xq_aI_key_at((_s_)->keys, _sat_))  /* sl_macro */

#define fcs_front_xq_aI_elem_assign(_s_, _d_)                           (*(_d_) = *(_s_))  /* sl_macro */
#define fcs_front_xq_aI_elem_assign_at(_s_, _sat_, _d_)                 ((_d_)->size = (_s_)->size - (_sat_), (_d_)->max_size = (_s_)->max_size - (_sat_), \
                                                         fcs_front_xq_aI_key_assign_at((_s_)->keys, _sat_, (_d_)->keys) \
                                                         fcs_front_xq_aI_cc_index_assign_at((_s_)->indices, _sat_, (_d_)->indices) \
                                                         fcs_front_xq_aI_cc_data_assign_at(_s_, _sat_, _d_))  /* sl_macro fcs_front_xq_aI_elem_assign_at */
#define fcs_front_xq_aI_elem_null(_e_)                                  ((_e_)->size = 0, (_e_)->max_size = 0, \
                                                         fcs_front_xq_aI_key_null((_e_)->keys) \
                                                         fcs_front_xq_aI_cc_index_null((_e_)->indices) \
                                                         fcs_front_xq_aI_cc_data_null(_e_))  /* sl_macro fcs_front_xq_aI_elem_null */
#define fcs_front_xq_aI_elem_inc(_e_)                                   (fcs_front_xq_aI_key_inc((_e_)->keys) \
                                                         fcs_front_xq_aI_cc_index_inc((_e_)->indices) \
                                                         fcs_front_xq_aI_cc_data_inc(_e_))  /* sl_macro fcs_front_xq_aI_elem_inc */
#define fcs_front_xq_aI_elem_dec(_e_)                                   (fcs_front_xq_aI_key_dec((_e_)->keys) \
                                                         fcs_front_xq_aI_cc_index_dec((_e_)->indices) \
                                                         fcs_front_xq_aI_cc_data_dec(_e_))  /* sl_macro fcs_front_xq_aI_elem_dec */
#define fcs_front_xq_aI_elem_add(_e_, _n_)                              (fcs_front_xq_aI_key_add((_e_)->keys, _n_) \
                                                         fcs_front_xq_aI_cc_index_add((_e_)->indices, _n_) \
                                                         fcs_front_xq_aI_cc_data_add(_e_, _n_))  /* sl_macro fcs_front_xq_aI_elem_add */
#define fcs_front_xq_aI_elem_sub(_e_, _n_)                              (fcs_front_xq_aI_key_sub((_e_)->keys, _n_) \
                                                         fcs_front_xq_aI_cc_index_sub((_e_)->indices, _n_) \
                                                         fcs_front_xq_aI_cc_data_sub(_e_, _n_))  /* sl_macro fcs_front_xq_aI_elem_sub */

#define fcs_front_xq_aI_elem_copy(_s_, _d_)                             (fcs_front_xq_aI_key_copy((_s_)->keys, (_d_)->keys) \
                                                         fcs_front_xq_aI_cc_index_copy((_s_)->indices, (_d_)->indices) \
                                                         fcs_front_xq_aI_cc_data_copy(_s_, _d_))  /* sl_macro fcs_front_xq_aI_elem_copy */
#define fcs_front_xq_aI_elem_ncopy(_s_, _d_, _n_)                       (fcs_front_xq_aI_key_ncopy((_s_)->keys, (_d_)->keys, _n_) \
                                                         fcs_front_xq_aI_cc_index_ncopy((_s_)->indices, (_d_)->indices, _n_) \
                                                         fcs_front_xq_aI_cc_data_ncopy(_s_, _d_, _n_))  /* sl_macro fcs_front_xq_aI_elem_ncopy */
#define fcs_front_xq_aI_elem_nmove(_s_, _d_, _n_)                       (fcs_front_xq_aI_key_nmove((_s_)->keys, (_d_)->keys, _n_) \
                                                         fcs_front_xq_aI_cc_index_nmove((_s_)->indices, (_d_)->indices, _n_) \
                                                         fcs_front_xq_aI_cc_data_nmove(_s_, _d_, _n_))  /* sl_macro fcs_front_xq_aI_elem_nmove */

#define fcs_front_xq_aI_elem_copy_at(_s_, _sat_, _d_, _dat_)            (fcs_front_xq_aI_key_copy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_) \
                                                         fcs_front_xq_aI_cc_index_copy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_) \
                                                         fcs_front_xq_aI_cc_data_copy_at(_s_, _sat_, _d_, _dat_))  /* sl_macro fcs_front_xq_aI_elem_copy_at */
#define fcs_front_xq_aI_elem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      (fcs_front_xq_aI_key_ncopy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_) \
                                                         fcs_front_xq_aI_cc_index_ncopy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_) \
                                                         fcs_front_xq_aI_cc_data_ncopy_at(_s_, _sat_, _d_, _dat_, _n_))  /* sl_macro fcs_front_xq_aI_elem_ncopy_at */
#define fcs_front_xq_aI_elem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      (fcs_front_xq_aI_key_nmove_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_) \
                                                         fcs_front_xq_aI_cc_index_ncopy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_) \
                                                         fcs_front_xq_aI_cc_data_nmove_at(_s_, _sat_, _d_, _dat_, _n_))  /* sl_macro fcs_front_xq_aI_elem_nmove_at */

#define fcs_front_xq_aI_elem_xchange(_e0_, _e1_, _t_)                   (fcs_front_xq_aI_key_xchange((_e0_)->keys, (_e1_)->keys, (_t_)->keys) \
                                                         fcs_front_xq_aI_cc_index_xchange((_e0_)->indices, (_e1_)->indices, (_t_)->indices) \
                                                         fcs_front_xq_aI_cc_data_xchange(_e0_, _e1_, _t_))  /* sl_macro fcs_front_xq_aI_elem_xchange */
#define fcs_front_xq_aI_elem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  (fcs_front_xq_aI_key_xchange_at((_e0_)->keys, _at0_, (_e1_)->keys, _at1_, (_t_)->keys) \
                                                         fcs_front_xq_aI_cc_index_xchange_at((_e0_)->indices, _at0_, (_e1_)->indices, _at1_, (_t_)->indices) \
                                                         fcs_front_xq_aI_cc_data_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_))  /* sl_macro fcs_front_xq_aI_elem_xchange_at */

#ifdef fcs_front_xq_aI_sl_elem_weight
# define fcs_front_xq_aI_elem_has_weight                                1  /* sl_macro */
# define fcs_front_xq_aI_elem_weight_ifelse(_if_, _el_)                 (_if_)  /* sl_macro */
# define fcs_front_xq_aI_elem_weight(_e_, _at_)                         ((fcs_front_xq_aI_slweight_t) fcs_front_xq_aI_sl_elem_weight((_e_), (_at_)))  /* sl_macro */
# define fcs_front_xq_aI_elem_weight_one(_e_, _at_)                     ((fcs_front_xq_aI_slweight_t) fcs_front_xq_aI_sl_elem_weight((_e_), (_at_)))  /* sl_macro */
# ifdef sl_elem_weight_set
#  define fcs_front_xq_aI_elem_weight_set(_e_, _at_, _w_)               sl_elem_weight_set((_e_), (_at_), (_w_))  /* sl_macro */
# else
#  define fcs_front_xq_aI_elem_weight_set(_e_, _at_, _w_)               fcs_front_xq_aI_sl_elem_weight((_e_), (_at_)) = (_w_)
# endif
#else
# define fcs_front_xq_aI_elem_has_weight                                0
# define fcs_front_xq_aI_elem_weight_ifelse(_if_, _el_)                 (_el_)
# define fcs_front_xq_aI_elem_weight_one(_e_, _at_)                     ((fcs_front_xq_aI_slweight_t) 1)
# define fcs_front_xq_aI_elem_weight_set(_e_, _at_, _w_)                Z_NOP()
#endif

#define fcs_front_xq_aI_elem_pack(_s_, _d_)                             (fcs_front_xq_aI_key_copy((_s_)->keys, &(_d_)->elements[0].key) fcs_front_xq_aI_cc_data_copy(_d_, &(_s_)->elements[0]))  /* sl_macro */
#define fcs_front_xq_aI_elem_pack_at(_s_, _sat_, _d_, _dat_)            (fcs_front_xq_aI_key_copy_at((_s_)->keys, _sat_, &(_d_)->elements[_dat_].key, 0) fcs_front_xq_aI_cc_data_copy_at(_s_, _sat_, &(_d_)->elements[_dat_], 0))  /* sl_macro */

#define fcs_front_xq_aI_elem_npack(_s_, _d_, _n_)                       elem_npack_at(_s_, 0, _d_, 0, _n_)  /* sl_macro */






#ifndef SL_RTI_TID_DECLARED
# define SL_RTI_TID_DECLARED

enum rti_tid
{
  /* base_mpi/base_mpi.c */
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






/* sl_macro fcs_front_xq_aI_SL_USE_RTI fcs_front_xq_aI_SL_USE_RTI_CMC fcs_front_xq_aI_SL_USE_RTI_TIM fcs_front_xq_aI_SL_USE_RTI_MEM */

#ifndef fcs_front_xq_aI_SL_USE_RTI

 #undef fcs_front_xq_aI_SL_USE_RTI_CMC  /* compare-move-counter */
 #undef fcs_front_xq_aI_SL_USE_RTI_TIM  /* timing */
 #undef fcs_front_xq_aI_SL_USE_RTI_MEM  /* memory */

#endif

#ifdef fcs_front_xq_aI_SL_USE_RTI_CMC

 /* regular commands */
 #define fcs_front_xq_aI_rti_cadd_cmp(n)          (fcs_front_xq_aI_SL_DEFCON(rti).cmc.cmp += n)  /* sl_macro */
 #define fcs_front_xq_aI_rti_cadd_movek(n)        (fcs_front_xq_aI_SL_DEFCON(rti).cmc.movek += n)  /* sl_macro */
 #define fcs_front_xq_aI_rti_cadd_moved(n)        (fcs_front_xq_aI_SL_DEFCON(rti).cmc.moved += n)  /* sl_macro */
 #define fcs_front_xq_aI_rti_cclear_cmp()         (fcs_front_xq_aI_SL_DEFCON(rti).cmc.cmp = 0)  /* sl_macro */
 #define fcs_front_xq_aI_rti_cclear_movek()       (fcs_front_xq_aI_SL_DEFCON(rti).cmc.movek = 0)  /* sl_macro */
 #define fcs_front_xq_aI_rti_cclear_moved()       (fcs_front_xq_aI_SL_DEFCON(rti).cmc.moved = 0)  /* sl_macro */
 #define fcs_front_xq_aI_rti_cclear_all()         (fcs_front_xq_aI_SL_DEFCON(rti).cmc.cmp = fcs_front_xq_aI_SL_DEFCON(rti).cmc.movek = fcs_front_xq_aI_SL_DEFCON(rti).cmc.moved = 0)  /* sl_macro */
 #define fcs_front_xq_aI_rti_ccmp()               my_rti_ccmp(fcs_front_xq_aI_SL_DEFCON(rti))  /* sl_macro */
 #define fcs_front_xq_aI_rti_cmovek()             my_rti_cmovek(fcs_front_xq_aI_SL_DEFCON(rti))  /* sl_macro */
 #define fcs_front_xq_aI_rti_cmoved()             my_rti_cmoved(fcs_front_xq_aI_SL_DEFCON(rti))  /* sl_macro */

 /* chained commands */
 #define fcs_front_xq_aI_cc_rti_cadd_cmp(n)       fcs_front_xq_aI_rti_cadd_cmp(n),  /* sl_macro */
 #define fcs_front_xq_aI_cc_rti_cadd_movek(n)     fcs_front_xq_aI_rti_cadd_movek(n),  /* sl_macro */
 #define fcs_front_xq_aI_cc_rti_cadd_moved(n)     fcs_front_xq_aI_rti_cadd_moved(n),  /* sl_macro */

#else /* fcs_front_xq_aI_SL_USE_RTI_CMC */

 /* regular commands */
 #define fcs_front_xq_aI_rti_cadd_cmp(n)
 #define fcs_front_xq_aI_rti_cadd_movek(n)
 #define fcs_front_xq_aI_rti_cadd_moved(n)
 #define fcs_front_xq_aI_rti_cclear_cmp()
 #define fcs_front_xq_aI_rti_cclear_movek()
 #define fcs_front_xq_aI_rti_cclear_moved()
 #define fcs_front_xq_aI_rti_cclear_all()
 #define fcs_front_xq_aI_rti_ccmp()               0
 #define fcs_front_xq_aI_rti_cmovek()             0
 #define fcs_front_xq_aI_rti_cmoved()             0

 /* chained commands */
 #define fcs_front_xq_aI_cc_rti_cadd_cmp(n)
 #define fcs_front_xq_aI_cc_rti_cadd_movek(n)
 #define fcs_front_xq_aI_cc_rti_cadd_moved(n)

#endif /* fcs_front_xq_aI_SL_USE_RTI_CMC */


#ifdef fcs_front_xq_aI_SL_USE_RTI_TIM

 #define fcs_front_xq_aI_rti_tstart(t)            (fcs_front_xq_aI_SL_DEFCON(rti).tim[t].start = z_time_get_s(), ++fcs_front_xq_aI_SL_DEFCON(rti).tim[t].num)  /* sl_macro */
 #define fcs_front_xq_aI_rti_tstop(t)             (fcs_front_xq_aI_SL_DEFCON(rti).tim[t].stop = z_time_get_s(), fcs_front_xq_aI_SL_DEFCON(rti).tim[t].cumu += (fcs_front_xq_aI_SL_DEFCON(rti).tim[t].last = fcs_front_xq_aI_SL_DEFCON(rti).tim[t].stop - fcs_front_xq_aI_SL_DEFCON(rti).tim[t].start))  /* sl_macro */
 #define fcs_front_xq_aI_rti_tclear(t)            (fcs_front_xq_aI_SL_DEFCON(rti).tim[t].last = 0)  /* sl_macro */
 #define fcs_front_xq_aI_rti_treset(t)            (fcs_front_xq_aI_SL_DEFCON(rti).tim[t].last = fcs_front_xq_aI_SL_DEFCON(rti).tim[t].cumu = 0, fcs_front_xq_aI_SL_DEFCON(rti).tim[t].num = 0)  /* sl_macro */
 #define fcs_front_xq_aI_rti_tlast(t)             my_rti_tlast(fcs_front_xq_aI_SL_DEFCON(rti), t)  /* sl_macro */
 #define fcs_front_xq_aI_rti_tcumu(t)             my_rti_tcumu(fcs_front_xq_aI_SL_DEFCON(rti), t)  /* sl_macro */
 #define fcs_front_xq_aI_rti_tnum(t)              my_rti_tnum(fcs_front_xq_aI_SL_DEFCON(rti), t)  /* sl_macro */

#else

 #define fcs_front_xq_aI_rti_tstart(t)            Z_NOP()
 #define fcs_front_xq_aI_rti_tstop(t)             Z_NOP()
 #define fcs_front_xq_aI_rti_tclear(t)            Z_NOP()
 #define fcs_front_xq_aI_rti_treset(t)            Z_NOP()
 #define fcs_front_xq_aI_rti_tlast(t)             0
 #define fcs_front_xq_aI_rti_tcumu(t)             0
 #define fcs_front_xq_aI_rti_tnum(t)              0

#endif


#ifdef fcs_front_xq_aI_SL_USE_RTI_MEM

 #define fcs_front_xq_aI_rti_minc_alloc()         ++fcs_front_xq_aI_SL_DEFCON(rti).mem.nalloc  /* sl_macro */
 #define fcs_front_xq_aI_rti_minc_free()          ++fcs_front_xq_aI_SL_DEFCON(rti).mem.nfree  /* sl_macro */
 #define fcs_front_xq_aI_rti_malloc(_s_)          (fcs_front_xq_aI_SL_DEFCON(rti).mem.max = z_max(_s_, fcs_front_xq_aI_SL_DEFCON(rti).mem.max), fcs_front_xq_aI_SL_DEFCON(rti).mem.cur += _s_, fcs_front_xq_aI_SL_DEFCON(rti).mem.cur_max = z_max(fcs_front_xq_aI_SL_DEFCON(rti).mem.cur, fcs_front_xq_aI_SL_DEFCON(rti).mem.cur_max))  /* sl_macro */
 #define fcs_front_xq_aI_rti_mfree(_s_)           (fcs_front_xq_aI_SL_DEFCON(rti).mem.cur -= _s_)  /* sl_macro */

 #define fcs_front_xq_aI_cc_rti_minc_alloc()      fcs_front_xq_aI_rti_minc_alloc(),  /* sl_macro */
 #define fcs_front_xq_aI_cc_rti_minc_free()       fcs_front_xq_aI_rti_minc_free(),  /* sl_macro */
 #define fcs_front_xq_aI_cc_rti_malloc(_s_)       fcs_front_xq_aI_rti_malloc(_s_),  /* sl_macro */
 #define fcs_front_xq_aI_cc_rti_mfree(_s_)        fcs_front_xq_aI_rti_mfree(_s_),  /* sl_macro */

#else

 #define fcs_front_xq_aI_rti_minc_alloc()         Z_NOP()
 #define fcs_front_xq_aI_rti_minc_free()          Z_NOP()
 #define fcs_front_xq_aI_rti_malloc(_s_)          Z_NOP()
 #define fcs_front_xq_aI_rti_mfree(_s_)           Z_NOP()

 #define fcs_front_xq_aI_cc_rti_minc_alloc()
 #define fcs_front_xq_aI_cc_rti_minc_free()
 #define fcs_front_xq_aI_cc_rti_malloc(_s_)
 #define fcs_front_xq_aI_cc_rti_mfree(_s_)

#endif


#ifdef fcs_front_xq_aI_SL_USE_RTI
 #define fcs_front_xq_aI_rti_reset()              my_rti_reset(fcs_front_xq_aI_SL_DEFCON(rti))  /* sl_macro */
#else
 #define fcs_front_xq_aI_rti_reset()              Z_NOP()
#endif






#define fcs_front_xq_aI_SPEC_TLOC

typedef fcs_front_xq_aI_sl_int_type_c fcs_front_xq_aI_spec_int_t;

typedef int fcs_front_xq_aI_spec_proc_t;

#define fcs_front_xq_aI_SPEC_LOC_NONE   -1
#ifdef SL_USE_MPI
# define fcs_front_xq_aI_SPEC_PROC_NONE  MPI_PROC_NULL
#else
# define fcs_front_xq_aI_SPEC_PROC_NONE  -1
#endif

typedef void *fcs_front_xq_aI_spec_tloc_data_t;
typedef void *fcs_front_xq_aI_spec_tproc_data_t;

struct fcs_front_xq_aI__elements_t;

typedef struct fcs_front_xq_aI__elements_t *fcs_front_xq_aI_spec_elem_buf_t;

typedef struct fcs_front_xq_aI__elements_t fcs_front_xq_aI_spec_elem_t;

typedef fcs_front_xq_aI_sl_int_type_c fcs_front_xq_aI_spec_elem_index_t;

#define fcs_front_xq_aI_spec_elem_unset(_e_)          fcs_front_xq_aI_elem_null((_e_))

#define fcs_front_xq_aI_spec_elem_set_n(_e_, _n_)     fcs_front_xq_aI_elem_set_size((_e_), (_n_))
#define fcs_front_xq_aI_spec_elem_get_n(_e_)          fcs_front_xq_aI_elem_get_size((_e_))
#define fcs_front_xq_aI_spec_elem_set_nmax(_e_, _n_)  fcs_front_xq_aI_elem_set_max_size((_e_), (_n_))
#define fcs_front_xq_aI_spec_elem_get_nmax(_e_)       fcs_front_xq_aI_elem_get_max_size((_e_))

#define fcs_front_xq_aI_spec_elem_set_buf(_e_, _b_)   *(_e_) = *(_b_)
#define fcs_front_xq_aI_spec_elem_get_buf(_e_)        (_e_)

#define fcs_front_xq_aI_spec_elem_copy_at(_se_, _sat_, _de_, _dat_) \
  fcs_front_xq_aI_elem_copy_at((_se_), (_sat_), (_de_), (_dat_))

#define fcs_front_xq_aI_spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) \
  fcs_front_xq_aI_elem_xchange_at((_s0_), (_s0at_), (_s1_), (_s1at_), (_t_))






/* tproc count */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_COUNT_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_COUNT_DB \
  struct { fcs_front_xq_aI_spec_elem_index_t n, i; fcs_front_xq_aI_spec_proc_t p; } spec0cd;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_COUNT_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0cd.n = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  for (spec0cd.i = 0; spec0cd.i < spec0cd.n; ++spec0cd.i) { \
    spec0cd.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == fcs_front_xq_aI_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_COUNT_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROC_COUNT_DB \
  fcs_front_xq_aI_SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_COUNT_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_COUNT_IP \
  struct { fcs_front_xq_aI_spec_elem_index_t n, t, i; fcs_front_xq_aI_spec_proc_t p; } spec0ci;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_COUNT_IP */
#define fcs_front_xq_aI_SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.n = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < spec0ci.n; ++spec0ci.i) { \
    spec0ci.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == fcs_front_xq_aI_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) fcs_front_xq_aI_spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  fcs_front_xq_aI_spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_COUNT_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROC_COUNT_IP \
  fcs_front_xq_aI_SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { fcs_front_xq_aI_spec_elem_index_t n, i; fcs_front_xq_aI_spec_proc_t p; } spec1cd;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_MOD_COUNT_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1cd.n = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  for (spec1cd.i = 0; spec1cd.i < spec1cd.n; ++spec1cd.i) { \
    spec1cd.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == fcs_front_xq_aI_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  fcs_front_xq_aI_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { fcs_front_xq_aI_spec_elem_index_t n, t, i; fcs_front_xq_aI_spec_proc_t p; } spec1ci;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_MOD_COUNT_IP */
#define fcs_front_xq_aI_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.n = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < spec1ci.n; ++spec1ci.i) { \
    spec1ci.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == fcs_front_xq_aI_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) fcs_front_xq_aI_spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  fcs_front_xq_aI_spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  fcs_front_xq_aI_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_COUNT_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { fcs_front_xq_aI_spec_elem_index_t n, i; fcs_front_xq_aI_spec_int_t j, m; } spec2cd;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_COUNT_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2cd.n = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  for (spec2cd.i = 0; spec2cd.i < spec2cd.n; ++spec2cd.i) { \
    (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), &spec2cd.m, (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.m; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_COUNT_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts, fcs_front_xq_aI_spec_proc_t *procs) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROCS_COUNT_DB \
  fcs_front_xq_aI_SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_COUNT_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { fcs_front_xq_aI_spec_elem_index_t n, t, i; fcs_front_xq_aI_spec_int_t j, m; } spec2ci;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_COUNT_IP */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.n = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < spec2ci.n; ++spec2ci.i) { \
    (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), &spec2ci.m, (_ps_)); \
    if (spec2ci.m <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.m; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) fcs_front_xq_aI_spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  fcs_front_xq_aI_spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_COUNT_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts, fcs_front_xq_aI_spec_proc_t *procs) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROCS_COUNT_IP \
  fcs_front_xq_aI_SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { fcs_front_xq_aI_spec_elem_index_t n, i; fcs_front_xq_aI_spec_int_t j, m; } spec3cd;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_COUNT_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3cd.n = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  for (spec3cd.i = 0; spec3cd.i < spec3cd.n; ++spec3cd.i) { \
    (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), &spec3cd.m, (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.m; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts, fcs_front_xq_aI_spec_proc_t *procs) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { fcs_front_xq_aI_spec_elem_index_t n, t, i; fcs_front_xq_aI_spec_int_t j, m; } spec3ci;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_COUNT_IP */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.n = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < spec3ci.n; ++spec3ci.i) { \
    (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), &spec3ci.m, (_ps_), NULL); \
    if (spec3ci.m <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.m; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) fcs_front_xq_aI_spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  fcs_front_xq_aI_spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts, fcs_front_xq_aI_spec_proc_t *procs) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* un-fixed macros, sp_macro fcs_front_xq_aI_spec_fixed_default_declare fcs_front_xq_aI_spec_fixed_default_create fcs_front_xq_aI_spec_fixed_default_copy_at fcs_front_xq_aI_spec_fixed_default_exchange_at fcs_front_xq_aI_spec_fixed_default_destroy */
#define fcs_front_xq_aI_spec_fixed_default_declare(_fx_, _fxp_)
#define fcs_front_xq_aI_spec_fixed_default_create(_fx_, _fxp_)
#define fcs_front_xq_aI_spec_fixed_default_copy_at(_se_, _sat_, _de_, _dat_, _fx_, _fxp_)             fcs_front_xq_aI_spec_elem_copy_at(_se_, _sat_, _de_, _dat_)
#define fcs_front_xq_aI_spec_fixed_default_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_, _fx_, _fxp_)  fcs_front_xq_aI_spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_)
#define fcs_front_xq_aI_spec_fixed_default_destroy(_fx_, _fxp_)


/* tproc rearrange */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { fcs_front_xq_aI_spec_elem_index_t n, i; fcs_front_xq_aI_spec_proc_t p; _fxdcl_(fx, _fxp_) } spec0d;

/* sp_macro fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  _fxc_(spec0d.fx, _fxp_); \
  spec0d.n = fcs_front_xq_aI_spec_elem_get_n(_sb_); \
  for (spec0d.i = 0; spec0d.i < spec0d.n; ++spec0d.i) { \
    spec0d.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == fcs_front_xq_aI_SPEC_PROC_NONE) continue; \
    _fxca_((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p], spec0d.fx, _fxp_); \
    ++(_ds_)[spec0d.p]; \
  } \
  _fxd_(spec0d.fx, _fxp_); \
  } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, int *displs) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB(_fxdcl_, _fxp_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_REARRANGE_DB \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_REARRANGE_DB(fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM)

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_REARRANGE_DB(NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_REARRANGE_DB(_name_, fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { fcs_front_xq_aI_spec_elem_index_t e, i, j; fcs_front_xq_aI_spec_proc_t p, np; } spec0i;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DO_TPROC_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_)  do { \
  for (spec0i.e = 0, spec0i.i = 0; spec0i.i < (_n_); ++spec0i.i) { \
    spec0i.e += (_cs_)[spec0i.i]; \
    spec0i.j = (_ds_)[spec0i.i]; \
    while (spec0i.j < spec0i.e) { \
      spec0i.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec0i.j, _tpd_); \
      while (spec0i.p != spec0i.i) { \
        spec0i.np = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), (_ds_)[spec0i.p], _tpd_); \
        if (spec0i.np != spec0i.p) fcs_front_xq_aI_spec_elem_exchange_at((_b_), (_ds_)[spec0i.p], (_b_), spec0i.j, (_xb_)); \
        ++(_ds_)[spec0i.p]; \
        spec0i.p = spec0i.np; \
      } \
      ++spec0i.j; \
    } \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_rearrange_ip(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, int *displs, int *counts, fcs_front_xq_aI_spec_int_t n) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROC_REARRANGE_IP \
  fcs_front_xq_aI_SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { fcs_front_xq_aI_spec_elem_index_t n, i; fcs_front_xq_aI_spec_proc_t p; _fxdcl_(fx, _fxp_) } spec1d;

/* sp_macro fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  spec1d.n = fcs_front_xq_aI_spec_elem_get_n(_sb_); \
  _fxc_(spec0d.fx, _fxp_); \
  for (spec1d.i = 0; spec1d.i < spec1d.n; ++spec1d.i) { \
    spec1d.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, fcs_front_xq_aI_spec_elem_get_buf(_ib_)); \
    if (spec1d.p == fcs_front_xq_aI_SPEC_PROC_NONE) continue; \
    _fxca_((_ib_), 0, (_db_), (_ds_)[spec1d.p], spec1d.fx, _fxp_); \
    ++(_ds_)[spec1d.p]; \
  } \
  _fxd_(spec0d.fx, _fxp_); \
  } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, int *displs, fcs_front_xq_aI_spec_elem_t *mod) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_MOD_REARRANGE_DB(fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM)

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_MOD_REARRANGE_DB(NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_, _ib_)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB(_name_, fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { fcs_front_xq_aI_spec_elem_index_t e, i, j; fcs_front_xq_aI_spec_proc_t p, np; } spec1i;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
    spec1i.e += (_cs_)[spec1i.i]; \
    spec1i.j = (_ds_)[spec1i.i]; \
    while (spec1i.j < spec1i.e) { \
      spec1i.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec1i.j, _tpd_, fcs_front_xq_aI_spec_elem_get_buf(_ib_)); \
      fcs_front_xq_aI_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
      while (spec1i.p != spec1i.i) { \
        spec1i.np = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, fcs_front_xq_aI_spec_elem_get_buf(_ib_)); \
        if (spec1i.np != spec1i.p) { \
          fcs_front_xq_aI_spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
          fcs_front_xq_aI_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        } else fcs_front_xq_aI_spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
        ++(_ds_)[spec1i.p]; \
        spec1i.p = spec1i.np; \
      } \
      ++spec1i.j; \
    } \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, int *displs, int *counts, fcs_front_xq_aI_spec_int_t n, fcs_front_xq_aI_spec_elem_t *mod) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  fcs_front_xq_aI_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { fcs_front_xq_aI_spec_elem_index_t n, i; fcs_front_xq_aI_spec_int_t j, m; _fxdcl_(fx, _fxp_) } spec2d;

/* sp_macro fcs_front_xq_aI_SPEC_DO_FIXED_TPROCS_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_FIXED_TPROCS_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  _fxc_(spec2d.fx, _fxp_); \
  spec2d.n = fcs_front_xq_aI_spec_elem_get_n(_sb_); \
  for (spec2d.i = 0; spec2d.i < spec2d.n; ++spec2d.i) { \
    (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), &spec2d.m, (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.m; ++spec2d.j) { \
      _fxca_((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]], spec2d.fx, _fxp_); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } \
  _fxd_(spec2d.fx, _fxp_); \
  } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, int *displs, fcs_front_xq_aI_spec_proc_t *procs) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB(_fxdcl_, _fxp_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROCS_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROCS_REARRANGE_DB(fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM)

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROCS_REARRANGE_DB(NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB(_name_, fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { fcs_front_xq_aI_spec_elem_index_t e, j, fe, fc, le, lc; fcs_front_xq_aI_spec_int_t i, n, f, l, o; } spec2i;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_)  do { \
  spec2i.f = 0; spec2i.fe = (_cs_)[0]; spec2i.fc = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; } \
  spec2i.l = 0; spec2i.le = (_cs_)[0]; spec2i.lc = fcs_front_xq_aI_spec_elem_get_n(_b_) - 1; \
  while (spec2i.lc >= spec2i.le) { ++spec2i.l; spec2i.le += (_cs_)[spec2i.l]; } \
  for (spec2i.e = 0, spec2i.i = 0; spec2i.i < (_n_); ++spec2i.i) { \
    spec2i.e += (_cs_)[spec2i.i]; \
    spec2i.j = (_ds_)[spec2i.i]; \
    while (spec2i.j < spec2i.e) { \
      (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec2i.j, (_tpd_), &spec2i.n, (_ps_)); \
      spec2i.o = -1; \
      while (spec2i.n > 0) { \
        --spec2i.n; \
        if ((_ps_)[spec2i.n] == spec2i.i && spec2i.o < 0) spec2i.o = spec2i.n; \
        else if ((_ds_)[(_ps_)[spec2i.n]] < spec2i.fc) { \
          spec2i.l = spec2i.f; spec2i.le = spec2i.fe; spec2i.lc = spec2i.fc; \
          if (spec2i.fc < spec2i.fe) { \
            fcs_front_xq_aI_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_b_), spec2i.fc); \
            ++spec2i.fc; \
          } else fcs_front_xq_aI_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec2i.n]] == spec2i.fc) ++spec2i.fc; \
        if (spec2i.j != (_ds_)[(_ps_)[spec2i.n]]) fcs_front_xq_aI_spec_elem_copy_at((_b_), spec2i.j, (_b_), (_ds_)[(_ps_)[spec2i.n]]); \
        ++(_ds_)[(_ps_)[spec2i.n]]; \
        while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; spec2i.fc = (_ds_)[spec2i.f]; } \
      } \
      if (spec2i.o < 0) { \
        if (spec2i.lc < spec2i.le) {  \
          fcs_front_xq_aI_spec_elem_copy_at((_b_), spec2i.lc, (_b_), spec2i.j); \
          spec2i.f = spec2i.l; spec2i.fe = spec2i.le; spec2i.fc = spec2i.lc; \
          --spec2i.lc; \
          while (spec2i.l > 0 && spec2i.lc < (_ds_)[spec2i.l]) { spec2i.le -= (_cs_)[spec2i.l]; spec2i.lc = spec2i.le - 1; --spec2i.l; } \
        } else fcs_front_xq_aI_spec_elem_copy_at((_xb_), 0, (_b_), spec2i.j); \
      } \
      spec2i.j = (_ds_)[spec2i.i]; \
    } \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_ip(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, int *displs, int *counts, fcs_front_xq_aI_spec_int_t n, fcs_front_xq_aI_spec_proc_t *procs) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  fcs_front_xq_aI_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  struct { fcs_front_xq_aI_spec_elem_index_t n, i; fcs_front_xq_aI_spec_int_t j, m; _fxdcl_(fx, _fxp_) } spec3d;

/* sp_macro fcs_front_xq_aI_SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  _fxc_(spec3d.fx, _fxp_); \
  spec3d.n = fcs_front_xq_aI_spec_elem_get_n(_sb_); \
  for (spec3d.i = 0; spec3d.i < spec3d.n; ++spec3d.i) { \
    (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), &spec3d.m, (_ps_), fcs_front_xq_aI_spec_elem_get_buf(_ib_)); \
    for (spec3d.j = 0; spec3d.j < spec3d.m; ++spec3d.j) { \
      _fxca_((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]], spec3d.fx, _fxp_); \
      ++(_ds_)[(_ps_)[spec3d.j]]; \
    } \
  } \
  _fxd_(spec3d.fx, _fxp_); \
  } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, int *displs, fcs_front_xq_aI_spec_proc_t *procs, fcs_front_xq_aI_spec_elem_t *mod) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB(_fxdcl_, _fxp_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROCS_MOD_REARRANGE_DB(fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM)

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROCS_MOD_REARRANGE_DB(NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB(_name_, fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { fcs_front_xq_aI_spec_elem_index_t e, j, fe, fc, le, lc; fcs_front_xq_aI_spec_int_t i, n, f, l, o; } spec3i;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = fcs_front_xq_aI_spec_elem_get_n(_b_); \
  while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
  spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = fcs_front_xq_aI_spec_elem_get_n(_b_) - 1; \
  while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
  for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
    spec3i.e += (_cs_)[spec3i.i]; \
    spec3i.j = (_ds_)[spec3i.i]; \
    while (spec3i.j < spec3i.e) { \
      (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), &spec3i.n, (_ps_), fcs_front_xq_aI_spec_elem_get_buf(_ib_)); \
      spec3i.o = -1; \
      while (spec3i.n > 0) { \
        --spec3i.n; \
        if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
        else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
          spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
          if (spec3i.fc < spec3i.fe) { \
            fcs_front_xq_aI_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
            ++spec3i.fc; \
          } else fcs_front_xq_aI_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
        fcs_front_xq_aI_spec_elem_copy_at((_ib_), spec3i.n, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
        ++(_ds_)[(_ps_)[spec3i.n]]; \
        while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
      } \
      if (spec3i.o < 0) { \
        if (spec3i.lc < spec3i.le) {  \
          fcs_front_xq_aI_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
          spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
          --spec3i.lc; \
          while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
        } else fcs_front_xq_aI_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
      } \
      spec3i.j = (_ds_)[spec3i.i]; \
    } \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, int *displs, int *counts, fcs_front_xq_aI_spec_int_t n, fcs_front_xq_aI_spec_proc_t *procs, fcs_front_xq_aI_spec_elem_t *mod) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}


/* tproc indices */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_INDICES_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_INDICES_DB \
  struct { fcs_front_xq_aI_spec_elem_index_t i; fcs_front_xq_aI_spec_proc_t p; } spec0xd;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_INDICES_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROC_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_)  do { \
  for (spec0xd.i = 0; spec0xd.i < fcs_front_xq_aI_spec_elem_get_n(_b_); ++spec0xd.i) { \
    spec0xd.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec0xd.i, (_tpd_)); \
    if (spec0xd.p == fcs_front_xq_aI_SPEC_PROC_NONE) continue; \
    (_ix_)[(_id_)[spec0xd.p]] = spec0xd.i; \
    ++(_id_)[spec0xd.p]; \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_INDICES_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_indices_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *indices, int *idispls) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROC_INDICES_DB \
  fcs_front_xq_aI_SPEC_DO_TPROC_INDICES_DB(_tp_, tproc_data, s, indices, idispls); \
}


/* tproc_mod indices */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_INDICES_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_INDICES_DB \
  struct { fcs_front_xq_aI_spec_elem_index_t i, k; fcs_front_xq_aI_spec_proc_t p; } spec1xd;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_MOD_INDICES_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROC_MOD_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_, _ib_, _d_)  do { \
  spec1xd.k = 0; \
  for (spec1xd.i = 0; spec1xd.i < fcs_front_xq_aI_spec_elem_get_n(_b_); ++spec1xd.i) { \
    spec1xd.p = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec1xd.i, (_tpd_), fcs_front_xq_aI_spec_elem_get_buf(_ib_)); \
    if (spec1xd.p == fcs_front_xq_aI_SPEC_PROC_NONE) continue; \
    fcs_front_xq_aI_spec_elem_copy_at((_ib_), 0, (_d_), spec1xd.k); \
    (_ix_)[(_id_)[spec1xd.p]] = spec1xd.k; \
    ++spec1xd.k; \
    ++(_id_)[spec1xd.p]; \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_INDICES_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_indices_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *indices, int *idispls, fcs_front_xq_aI_spec_elem_t *mod, fcs_front_xq_aI_spec_elem_t *d) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROC_MOD_INDICES_DB \
  fcs_front_xq_aI_SPEC_DO_TPROC_MOD_INDICES_DB(_tp_, tproc_data, s, indices, idispls, mod, d); \
}


/* tprocs indices */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_INDICES_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_INDICES_DB \
  struct { fcs_front_xq_aI_spec_elem_index_t i; fcs_front_xq_aI_spec_int_t j, n; } spec2xd;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_INDICES_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_, _ps_)  do { \
  for (spec2xd.i = 0; spec2xd.i < fcs_front_xq_aI_spec_elem_get_n(_b_); ++spec2xd.i) { \
    (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec2xd.i, (_tpd_), &spec2xd.n, (_ps_)); \
    for (spec2xd.j = 0; spec2xd.j < spec2xd.n; ++spec2xd.j) { \
      (_ix_)[(_id_)[(_ps_)[spec2xd.j]]] = spec2xd.i; \
      ++(_id_)[(_ps_)[spec2xd.j]]; \
    } \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_INDICES_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_indices_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *indices, int *idispls, fcs_front_xq_aI_spec_proc_t *procs) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROCS_INDICES_DB \
  fcs_front_xq_aI_SPEC_DO_TPROCS_INDICES_DB(_tp_, tproc_data, s, indices, idispls, procs); \
}


/* tprocs_mod indices */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_INDICES_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_INDICES_DB \
  struct { fcs_front_xq_aI_spec_elem_index_t i, k; fcs_front_xq_aI_spec_int_t j, n; } spec3xd;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_INDICES_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_INDICES_DB(_tp_, _tpd_, _b_, _ix_, _id_, _ps_, _ib_, _d_)  do { \
  spec3xd.k = 0; \
  for (spec3xd.i = 0; spec3xd.i < fcs_front_xq_aI_spec_elem_get_n(_b_); ++spec3xd.i) { \
    (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec3xd.i, (_tpd_), &spec3xd.n, (_ps_), fcs_front_xq_aI_spec_elem_get_buf(_ib_)); \
    for (spec3xd.j = 0; spec3xd.j < spec3xd.n; ++spec3xd.j) { \
      fcs_front_xq_aI_spec_elem_copy_at((_ib_), spec3xd.j, (_d_), spec3xd.k); \
      (_ix_)[(_id_)[(_ps_)[spec3xd.j]]] = spec3xd.k; \
      ++spec3xd.k; \
      ++(_id_)[(_ps_)[spec3xd.j]]; \
    } \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_INDICES_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_INDICES_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_indices_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *indices, int *idispls, fcs_front_xq_aI_spec_proc_t *procs, fcs_front_xq_aI_spec_elem_t *mod, fcs_front_xq_aI_spec_elem_t *d) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TPROCS_MOD_INDICES_DB \
  fcs_front_xq_aI_SPEC_DO_TPROCS_MOD_INDICES_DB(_tp_, tproc_data, s, indices, idispls, procs, mod, d); \
}


/* tproc sendrecv */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(_fxdcl_, _fxp_)
/*#define fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(_fxdcl_, _fxp_) \
  struct { _fxdcl_(fx, _fxp_) } spec0srd;*/

/* sp_macro fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_SENDRECV_DB */
#define fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_SENDRECV_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _tpd_, _sb_, _rb_, _sc_, _sd_, _rd_, _ab_, _ad_, _as_, _aq_, _aqn_, _aqs_, _r_, _p_)  do { \
  _fxc_(spec0srd.fx, _fxp_); \
  while (*(_sd_) < (_sc_)) { \
    if ((_p_) == fcs_front_xq_aI_SPEC_PROC_NONE) (_p_) = (_tp_)(fcs_front_xq_aI_spec_elem_get_buf(_sb_), *(_sd_), (_tpd_)); \
    if ((_p_) != fcs_front_xq_aI_SPEC_PROC_NONE) { \
      if ((_p_) == (_r_)) { \
        _fxca_((_sb_), *(_sd_), (_rb_), *(_rd_), spec0srd.fx, _fxp_); \
        ++(*(_rd_)); \
      } else { \
        if ((_ad_)[_p_] >= ((_p_) + 1) * (_as_)) break; \
        _fxca_((_sb_), *(_sd_), (_ab_), (_ad_)[_p_], spec0srd.fx, _fxp_); \
        ++(_ad_)[_p_]; \
        if ((_ad_)[_p_] >= ((_p_) + 1) * (_as_)) { \
          (_aq_)[*(_aqn_)] = (_p_); ++(*(_aqn_)); *(_aqn_) %= (_aqs_); \
        } \
      } \
    } \
    (_p_) = fcs_front_xq_aI_SPEC_PROC_NONE; \
    ++(*(_sd_)); \
 } \
 _fxd_(spec0srd.fx, _fxp_); \
 } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_SENDRECV_DB */
#define fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_SENDRECV_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
_s_ fcs_front_xq_aI_spec_proc_t _name_##_tproc_sendrecv_db(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *sb, fcs_front_xq_aI_spec_elem_t *rb, fcs_front_xq_aI_spec_int_t scount, fcs_front_xq_aI_spec_int_t *sdispl, fcs_front_xq_aI_spec_int_t *rdispl, fcs_front_xq_aI_spec_elem_t *ax, fcs_front_xq_aI_spec_int_t *aux_displs, fcs_front_xq_aI_spec_int_t aux_size_max, fcs_front_xq_aI_spec_int_t *aux_queue, fcs_front_xq_aI_spec_int_t *aux_queue_next, fcs_front_xq_aI_spec_int_t aux_queue_size, fcs_front_xq_aI_spec_proc_t rank, fcs_front_xq_aI_spec_proc_t p) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(_fxdcl_, _fxp_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_SENDRECV_DB(_fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, tproc_data, sb, rb, scount, sdispl, rdispl, ax, aux_displs, aux_size_max, aux_queue, aux_queue_next, aux_queue_size, rank, p); \
  return p; \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TPROC_SENDRECV_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TPROC_SENDRECV_DB \
  fcs_front_xq_aI_SPEC_DECLARE_FIXED_TPROC_SENDRECV_DB(fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM)

/* sp_macro fcs_front_xq_aI_SPEC_DO_TPROC_SENDRECV_DB */
#define fcs_front_xq_aI_SPEC_DO_TPROC_SENDRECV_DB(_tp_, _tpd_, _sb_, _rb_, _sc_, _sd_, _rd_, _ab_, _ad_, _as_, _aq_, _aqn_, _aqs_, _r_, _p_) \
  fcs_front_xq_aI_SPEC_DO_FIXED_TPROC_SENDRECV_DB(NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _tpd_, _sb_, _rb_, _sc_, _sd_, _rd_, _ab_, _ad_, _as_, _aq_, _aqn_, _aqs_, _r_, _p_)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TPROC_SENDRECV_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TPROC_SENDRECV_DB(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_SENDRECV_DB(_name_, fcs_front_xq_aI_spec_fixed_default_declare, NOPARAM, fcs_front_xq_aI_spec_fixed_default_create, fcs_front_xq_aI_spec_fixed_default_copy_at, fcs_front_xq_aI_spec_fixed_default_exchange_at, fcs_front_xq_aI_spec_fixed_default_destroy, _tp_, _s_)


/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_TPROC */
#define fcs_front_xq_aI_SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_INDICES_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_SENDRECV_DB(_name_, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_TPROC_MOD */
#define fcs_front_xq_aI_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_TPROCS */
#define fcs_front_xq_aI_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_TPROCS_MOD */
#define fcs_front_xq_aI_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_FIXED_TPROC */
#define fcs_front_xq_aI_SPEC_DEFINE_FIXED_TPROC(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_INDICES_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_SENDRECV_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_FIXED_TPROC_MOD */
#define fcs_front_xq_aI_SPEC_DEFINE_FIXED_TPROC_MOD(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROC_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROC_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_FIXED_TPROCS */
#define fcs_front_xq_aI_SPEC_DEFINE_FIXED_TPROCS(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROCS_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_INDICES_DB(_name_, _tp_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_FIXED_TPROCS_MOD */
#define fcs_front_xq_aI_SPEC_DEFINE_FIXED_TPROCS_MOD(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_FIXED_TPROCS_MOD_REARRANGE_DB(_name_, _fxdcl_, _fxp_, _fxc_, _fxca_, _fxxa_, _fxd_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TPROCS_MOD_INDICES_DB(_name_, _tp_, _s_)

/* sp_type fcs_front_xq_aI_spec_tproc_f fcs_front_xq_aI_spec_tproc_count_f fcs_front_xq_aI_spec_tproc_rearrange_db_f fcs_front_xq_aI_spec_tproc_rearrange_ip_f fcs_front_xq_aI_spec_tproc_indices_db_f fcs_front_xq_aI_spec_tproc_sendrecv_db_f */
typedef fcs_front_xq_aI_spec_proc_t fcs_front_xq_aI_spec_tproc_f(fcs_front_xq_aI_spec_elem_buf_t b, fcs_front_xq_aI_spec_elem_index_t x, fcs_front_xq_aI_spec_tproc_data_t tproc_data);
typedef void fcs_front_xq_aI_spec_tproc_count_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts);
typedef void fcs_front_xq_aI_spec_tproc_rearrange_db_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, int *displs);
typedef void fcs_front_xq_aI_spec_tproc_rearrange_ip_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, int *displs, int *counts, fcs_front_xq_aI_spec_int_t n);
typedef void fcs_front_xq_aI_spec_tproc_indices_db_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *indices, int *idispls);
typedef fcs_front_xq_aI_spec_proc_t fcs_front_xq_aI_spec_tproc_sendrecv_db_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *sb, fcs_front_xq_aI_spec_elem_t *rb, fcs_front_xq_aI_spec_int_t scount, fcs_front_xq_aI_spec_int_t *sdispl, fcs_front_xq_aI_spec_int_t *rdispl, fcs_front_xq_aI_spec_elem_t *ax, fcs_front_xq_aI_spec_int_t *aux_displs, fcs_front_xq_aI_spec_int_t aux_size_max, fcs_front_xq_aI_spec_int_t *aux_queue, fcs_front_xq_aI_spec_int_t *aux_queue_next, fcs_front_xq_aI_spec_int_t aux_queue_size, fcs_front_xq_aI_spec_proc_t rank, fcs_front_xq_aI_spec_proc_t p);

/* sp_type fcs_front_xq_aI_spec_tproc_mod_f fcs_front_xq_aI_spec_tproc_mod_count_f fcs_front_xq_aI_spec_tproc_mod_rearrange_db_f fcs_front_xq_aI_spec_tproc_mod_rearrange_ip_f fcs_front_xq_aI_spec_tproc_mod_indices_db_f */
typedef fcs_front_xq_aI_spec_proc_t fcs_front_xq_aI_spec_tproc_mod_f(fcs_front_xq_aI_spec_elem_buf_t b, fcs_front_xq_aI_spec_elem_index_t x, fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_buf_t mod);
typedef void fcs_front_xq_aI_spec_tproc_mod_count_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts);
typedef void fcs_front_xq_aI_spec_tproc_mod_rearrange_db_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, int *displs, fcs_front_xq_aI_spec_elem_t *mod);
typedef void fcs_front_xq_aI_spec_tproc_mod_rearrange_ip_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, int *displs, int *counts, fcs_front_xq_aI_spec_int_t n, fcs_front_xq_aI_spec_elem_t *mod);
typedef void fcs_front_xq_aI_spec_tproc_mod_indices_db_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *indices, int *idispls, fcs_front_xq_aI_spec_elem_t *mod, fcs_front_xq_aI_spec_elem_t *d);

/* sp_type fcs_front_xq_aI_spec_tprocs_f fcs_front_xq_aI_spec_tprocs_count_f fcs_front_xq_aI_spec_tprocs_rearrange_db_f fcs_front_xq_aI_spec_tprocs_rearrange_ip_f fcs_front_xq_aI_spec_tprocs_indices_db_f */
typedef void fcs_front_xq_aI_spec_tprocs_f(fcs_front_xq_aI_spec_elem_buf_t b, fcs_front_xq_aI_spec_elem_index_t x, fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_int_t *nprocs, fcs_front_xq_aI_spec_proc_t *procs);
typedef void fcs_front_xq_aI_spec_tprocs_count_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts, fcs_front_xq_aI_spec_proc_t *procs);
typedef void fcs_front_xq_aI_spec_tprocs_rearrange_db_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, int *displs, fcs_front_xq_aI_spec_proc_t *procs);
typedef void fcs_front_xq_aI_spec_tprocs_rearrange_ip_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, int *displs, int *counts, fcs_front_xq_aI_spec_int_t n, fcs_front_xq_aI_spec_proc_t *procs);
typedef void fcs_front_xq_aI_spec_tprocs_indices_db_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *indices, int *idispls, fcs_front_xq_aI_spec_proc_t *procs);

/* sp_type fcs_front_xq_aI_spec_tprocs_mod_f fcs_front_xq_aI_spec_tprocs_mod_count_f fcs_front_xq_aI_spec_tprocs_mod_rearrange_db_f fcs_front_xq_aI_spec_tprocs_mod_rearrange_ip_f fcs_front_xq_aI_spec_tprocs_mod_indices_db_f */
typedef void fcs_front_xq_aI_spec_tprocs_mod_f(fcs_front_xq_aI_spec_elem_buf_t b, fcs_front_xq_aI_spec_elem_index_t x, fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_int_t *nprocs, fcs_front_xq_aI_spec_proc_t *procs, fcs_front_xq_aI_spec_elem_buf_t mod);
typedef void fcs_front_xq_aI_spec_tprocs_mod_count_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *counts, fcs_front_xq_aI_spec_proc_t *procs);
typedef void fcs_front_xq_aI_spec_tprocs_mod_rearrange_db_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, int *displs, fcs_front_xq_aI_spec_proc_t *procs, fcs_front_xq_aI_spec_elem_t *mod);
typedef void fcs_front_xq_aI_spec_tprocs_mod_rearrange_ip_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, int *displs, int *counts, fcs_front_xq_aI_spec_int_t n, fcs_front_xq_aI_spec_proc_t *procs, fcs_front_xq_aI_spec_elem_t *mod);
typedef void fcs_front_xq_aI_spec_tprocs_mod_indices_db_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data, fcs_front_xq_aI_spec_elem_t *s, int *indices, int *idispls, fcs_front_xq_aI_spec_proc_t *procs, fcs_front_xq_aI_spec_elem_t *mod, fcs_front_xq_aI_spec_elem_t *d);

/* sp_type fcs_front_xq_aI_spec_tproc_reset_f */
typedef void fcs_front_xq_aI_spec_tproc_reset_f(fcs_front_xq_aI_spec_tproc_data_t tproc_data);


/* sp_type fcs_front_xq_aI__spec_tproc_ext_t fcs_front_xq_aI_spec_tproc_ext_t */
typedef struct fcs_front_xq_aI__spec_tproc_ext_t
{
  fcs_front_xq_aI_spec_tproc_count_f *count_db, *count_ip;
  fcs_front_xq_aI_spec_tproc_rearrange_db_f *rearrange_db;
  fcs_front_xq_aI_spec_tproc_rearrange_ip_f *rearrange_ip;
  fcs_front_xq_aI_spec_tproc_indices_db_f *indices_db;
  fcs_front_xq_aI_spec_tproc_sendrecv_db_f *sendrecv_db;

} fcs_front_xq_aI_spec_tproc_ext_t;

/* sp_type fcs_front_xq_aI__spec_tproc_mod_ext_tproc_t fcs_front_xq_aI_spec_tproc_mod_ext_t */
typedef struct fcs_front_xq_aI__spec_tproc_mod_ext_tproc_t
{
  fcs_front_xq_aI_spec_tproc_mod_count_f *count_db, *count_ip;
  fcs_front_xq_aI_spec_tproc_mod_rearrange_db_f *rearrange_db;
  fcs_front_xq_aI_spec_tproc_mod_rearrange_ip_f *rearrange_ip;
  fcs_front_xq_aI_spec_tproc_mod_indices_db_f *indices_db;

} fcs_front_xq_aI_spec_tproc_mod_ext_t;

/* sp_type fcs_front_xq_aI__spec_tprocs_ext_t fcs_front_xq_aI_spec_tprocs_ext_t */
typedef struct fcs_front_xq_aI__spec_tprocs_ext_t
{
  fcs_front_xq_aI_spec_tprocs_count_f *count_db, *count_ip;
  fcs_front_xq_aI_spec_tprocs_rearrange_db_f *rearrange_db;
  fcs_front_xq_aI_spec_tprocs_rearrange_ip_f *rearrange_ip;
  fcs_front_xq_aI_spec_tprocs_indices_db_f *indices_db;

} fcs_front_xq_aI_spec_tprocs_ext_t;

/* sp_type fcs_front_xq_aI__spec_tprocs_mod_ext_t fcs_front_xq_aI_spec_tprocs_mod_ext_t */
typedef struct fcs_front_xq_aI__spec_tprocs_mod_ext_t
{
  fcs_front_xq_aI_spec_tprocs_mod_count_f *count_db, *count_ip;
  fcs_front_xq_aI_spec_tprocs_mod_rearrange_db_f *rearrange_db;
  fcs_front_xq_aI_spec_tprocs_mod_rearrange_ip_f *rearrange_ip;
  fcs_front_xq_aI_spec_tprocs_mod_indices_db_f *indices_db;

} fcs_front_xq_aI_spec_tprocs_mod_ext_t;

/* sp_macro fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC(_name_)       { _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip, _name_##_tproc_indices_db, _name_##_tproc_sendrecv_db }
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL          { NULL, NULL, NULL, NULL, NULL, NULL }
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD(_name_)   { _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip, _name_##_tproc_mod_indices_db }
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL      { NULL, NULL, NULL, NULL, NULL }
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS(_name_)      { _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip, _name_##_tprocs_indices_db }
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL         { NULL, NULL, NULL, NULL, NULL }
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD(_name_)  { _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip, _name_##_tprocs_mod_indices_db }
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL     { NULL, NULL, NULL, NULL, NULL }


/* enable tloc features */
#ifdef fcs_front_xq_aI_SPEC_TLOC

/* sp_macro fcs_front_xq_aI_SPEC_TLOC fcs_front_xq_aI_SPEC_LOC_NONE */


/* tloc rearrange */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TLOC_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TLOC_REARRANGE_DB \
  struct { fcs_front_xq_aI_spec_int_t i, p; } spec0d;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TLOC_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_TLOC_REARRANGE_DB(_tl_, _tld_, _sb_, _db_)  do { \
  for (spec0d.i = 0; spec0d.i < fcs_front_xq_aI_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tl_)(fcs_front_xq_aI_spec_elem_get_buf(_sb_), spec0d.i, _tld_); \
    if (spec0d.p == fcs_front_xq_aI_SPEC_LOC_NONE) continue; \
    fcs_front_xq_aI_spec_elem_copy_at((_sb_), spec0d.i, (_db_), spec0d.p); \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TLOC_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_rearrange_db(fcs_front_xq_aI_spec_tloc_data_t tloc_data, fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TLOC_REARRANGE_DB \
  fcs_front_xq_aI_SPEC_DO_TLOC_REARRANGE_DB(_tl_, tloc_data, s, d); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TLOC_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TLOC_REARRANGE_IP \
  struct { fcs_front_xq_aI_spec_int_t i, p, np; } spec0i;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TLOC_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DO_TLOC_REARRANGE_IP(_tl_, _tld_, _b_, _xb_)  do { \
  for (spec0i.i = 0; spec0i.i < fcs_front_xq_aI_spec_elem_get_n(_b_); ++spec0i.i) { \
    spec0i.p = (_tl_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec0i.i, _tld_); \
    if (spec0i.p == fcs_front_xq_aI_SPEC_LOC_NONE) continue; \
    while (spec0i.i != spec0i.p) { \
      spec0i.np = (_tl_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec0i.p, _tld_); \
      if (spec0i.np == fcs_front_xq_aI_SPEC_LOC_NONE) { fcs_front_xq_aI_spec_elem_copy_at((_b_), spec0i.i, (_b_), spec0i.p); break; } \
      fcs_front_xq_aI_spec_elem_exchange_at((_b_), spec0i.i, (_b_), spec0i.p, (_xb_)); \
      spec0i.p = spec0i.np; \
    } \
  } } while (0)

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TLOC_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_rearrange_ip(fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, fcs_front_xq_aI_spec_tloc_data_t tloc_data) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TLOC_REARRANGE_IP \
  fcs_front_xq_aI_SPEC_DO_TLOC_REARRANGE_IP(_tl_, tloc_data, s, x); \
}


/* tloc_mod_mod rearrange */

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  struct { fcs_front_xq_aI_spec_int_t i, p; } spec1d;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TLOC_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, _tld_, _sb_, _db_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < fcs_front_xq_aI_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(fcs_front_xq_aI_spec_elem_get_buf(_sb_), spec1d.i, _tld_, fcs_front_xq_aI_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == fcs_front_xq_aI_SPEC_LOC_NONE) continue; \
      fcs_front_xq_aI_spec_elem_copy_at((_ib_), 0, (_db_), spec1d.p); \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < fcs_front_xq_aI_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(fcs_front_xq_aI_spec_elem_get_buf(_sb_), spec1d.i, _tld_, NULL); \
      if (spec1d.p == fcs_front_xq_aI_SPEC_LOC_NONE) continue; \
      fcs_front_xq_aI_spec_elem_copy_at((_sb_), spec1d.i, (_db_), spec1d.p); \
    } \
  } } while (0) 

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TLOC_MOD_REARRANGE_DB */
#define fcs_front_xq_aI_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_mod_rearrange_db(fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, fcs_front_xq_aI_spec_tloc_data_t tloc_data, fcs_front_xq_aI_spec_elem_t *mod) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  fcs_front_xq_aI_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, tloc_data, s, d, mod); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  struct { fcs_front_xq_aI_spec_int_t i, p, np; } spec1i;

/* sp_macro fcs_front_xq_aI_SPEC_DO_TLOC_MOD_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, _tld_, _b_, _xb_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.i = 0; spec1i.i < fcs_front_xq_aI_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec1i.i, _tld_, fcs_front_xq_aI_spec_elem_get_buf(_ib_)); \
      if (spec1i.p == fcs_front_xq_aI_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec1i.p, _tld_, fcs_front_xq_aI_spec_elem_get_buf(_xb_)); \
        if (spec1i.np == fcs_front_xq_aI_SPEC_LOC_NONE) break; \
        fcs_front_xq_aI_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.p); \
        fcs_front_xq_aI_spec_elem_copy_at((_xb_), 0, (_ib_), 0); \
        spec1i.p = spec1i.np; \
      } \
      fcs_front_xq_aI_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.i); \
    } \
  } else { \
    for (spec1i.i = 0; spec1i.i < fcs_front_xq_aI_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec1i.i, _tld_, NULL); \
      if (spec1i.p == fcs_front_xq_aI_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(fcs_front_xq_aI_spec_elem_get_buf(_b_), spec1i.p, _tld_, NULL); \
        if (spec1i.np == fcs_front_xq_aI_SPEC_LOC_NONE) { fcs_front_xq_aI_spec_elem_copy_at((_b_), spec1i.i, (_b_), spec1i.p); break; } \
        fcs_front_xq_aI_spec_elem_exchange_at((_b_), spec1i.i, (_b_), spec1i.p, (_xb_)); \
        spec1i.p = spec1i.np; \
      } \
    } \
 } } while (0) 

/* sp_macro fcs_front_xq_aI_SPEC_FUNC_TLOC_MOD_REARRANGE_IP */
#define fcs_front_xq_aI_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_mod_rearrange_ip(fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, fcs_front_xq_aI_spec_tloc_data_t tloc_data, fcs_front_xq_aI_spec_elem_t *mod) \
{ \
  fcs_front_xq_aI_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  fcs_front_xq_aI_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, tloc_data, s, x, mod); \
}

/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_TLOC */
#define fcs_front_xq_aI_SPEC_DEFINE_TLOC(_name_, _tl_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_DEFINE_TLOC_MOD */
#define fcs_front_xq_aI_SPEC_DEFINE_TLOC_MOD(_name_, _tl_, _s_...) \
  fcs_front_xq_aI_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_) \
  fcs_front_xq_aI_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_NULL fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD_NULL */
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC(_name_)      _name_##_tloc_rearrange_db, _name_##_tloc_rearrange_ip
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_NULL         NULL, NULL
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD(_name_)  _name_##_tloc_mod_rearrange_db, _name_##_tloc_mod_rearrange_ip
#define fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD_NULL     NULL, NULL


/* sp_type fcs_front_xq_aI_spec_tloc_f fcs_front_xq_aI_spec_tloc_rearrange_db_f fcs_front_xq_aI_spec_tloc_rearrange_ip_f */
typedef fcs_front_xq_aI_spec_elem_index_t fcs_front_xq_aI_spec_tloc_f(fcs_front_xq_aI_spec_elem_buf_t b, fcs_front_xq_aI_spec_elem_index_t x, fcs_front_xq_aI_spec_tloc_data_t tloc_data);
typedef void fcs_front_xq_aI_spec_tloc_rearrange_db_f(fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, fcs_front_xq_aI_spec_tloc_data_t tloc_data);
typedef void fcs_front_xq_aI_spec_tloc_rearrange_ip_f(fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, fcs_front_xq_aI_spec_tloc_data_t tloc_data);

/* sp_type fcs_front_xq_aI_spec_tloc_mod_f fcs_front_xq_aI_spec_tloc_mod_rearrange_db_f fcs_front_xq_aI_spec_tloc_mod_rearrange_ip_f */
typedef fcs_front_xq_aI_spec_elem_index_t fcs_front_xq_aI_spec_tloc_mod_f(fcs_front_xq_aI_spec_elem_buf_t b, fcs_front_xq_aI_spec_elem_index_t x, fcs_front_xq_aI_spec_tloc_data_t tloc_data, fcs_front_xq_aI_spec_elem_buf_t mod);
typedef void fcs_front_xq_aI_spec_tloc_mod_rearrange_db_f(fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *d, fcs_front_xq_aI_spec_tloc_data_t tloc_data, fcs_front_xq_aI_spec_elem_t *mod);
typedef void fcs_front_xq_aI_spec_tloc_mod_rearrange_ip_f(fcs_front_xq_aI_spec_elem_t *s, fcs_front_xq_aI_spec_elem_t *x, fcs_front_xq_aI_spec_tloc_data_t tloc_data, fcs_front_xq_aI_spec_elem_t *mod);


#endif /* fcs_front_xq_aI_SPEC_TLOC */






#ifdef SL_USE_MPI
# include <mpi.h>
#endif


/* sl_type fcs_front_xq_aI_slint_t fcs_front_xq_aI_slint */
typedef fcs_front_xq_aI_sl_int_type_c fcs_front_xq_aI_slint_t, fcs_front_xq_aI_slint;  /* deprecated 'fcs_front_xq_aI_slint' */

#define fcs_front_xq_aI_slint_fmt   fcs_front_xq_aI_sl_int_type_fmt    /* sl_macro */

/* sl_type fcs_front_xq_aI_slindex_t */
typedef fcs_front_xq_aI_sl_index_type_c fcs_front_xq_aI_slindex_t;

#define fcs_front_xq_aI_sindex_fmt  fcs_front_xq_aI_sl_index_type_fmt  /* sl_macro */

/* sl_type fcs_front_xq_aI_slkey_t */
typedef fcs_front_xq_aI_sl_key_type_c fcs_front_xq_aI_slkey_t;

/* sl_type fcs_front_xq_aI_slkey_pure_t fcs_front_xq_aI_slpkey_t */
typedef fcs_front_xq_aI_sl_key_pure_type_c fcs_front_xq_aI_slkey_pure_t, fcs_front_xq_aI_slpkey_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type fcs_front_xq_aI_sldata0_t */
#ifdef fcs_front_xq_aI_sl_data0_type_c
typedef fcs_front_xq_aI_sl_data0_type_c fcs_front_xq_aI_sldata0_t;
#endif
/* sl_type fcs_front_xq_aI_sldata1_t */
#ifdef fcs_front_xq_aI_sl_data1_type_c
typedef fcs_front_xq_aI_sl_data1_type_c fcs_front_xq_aI_sldata1_t;
#endif
/* sl_type fcs_front_xq_aI_sldata2_t */
#ifdef fcs_front_xq_aI_sl_data2_type_c
typedef fcs_front_xq_aI_sl_data2_type_c fcs_front_xq_aI_sldata2_t;
#endif
/* sl_type fcs_front_xq_aI_sldata3_t */
#ifdef fcs_front_xq_aI_sl_data3_type_c
typedef fcs_front_xq_aI_sl_data3_type_c fcs_front_xq_aI_sldata3_t;
#endif
/* sl_type fcs_front_xq_aI_sldata4_t */
#ifdef fcs_front_xq_aI_sl_data4_type_c
typedef fcs_front_xq_aI_sl_data4_type_c fcs_front_xq_aI_sldata4_t;
#endif
/* sl_type fcs_front_xq_aI_sldata5_t */
#ifdef fcs_front_xq_aI_sl_data5_type_c
typedef fcs_front_xq_aI_sl_data5_type_c fcs_front_xq_aI_sldata5_t;
#endif
/* sl_type fcs_front_xq_aI_sldata6_t */
#ifdef fcs_front_xq_aI_sl_data6_type_c
typedef fcs_front_xq_aI_sl_data6_type_c fcs_front_xq_aI_sldata6_t;
#endif
/* sl_type fcs_front_xq_aI_sldata7_t */
#ifdef fcs_front_xq_aI_sl_data7_type_c
typedef fcs_front_xq_aI_sl_data7_type_c fcs_front_xq_aI_sldata7_t;
#endif
/* sl_type fcs_front_xq_aI_sldata8_t */
#ifdef fcs_front_xq_aI_sl_data8_type_c
typedef fcs_front_xq_aI_sl_data8_type_c fcs_front_xq_aI_sldata8_t;
#endif
/* sl_type fcs_front_xq_aI_sldata9_t */
#ifdef fcs_front_xq_aI_sl_data9_type_c
typedef fcs_front_xq_aI_sl_data9_type_c fcs_front_xq_aI_sldata9_t;
#endif
/* sl_type fcs_front_xq_aI_sldata10_t */
#ifdef fcs_front_xq_aI_sl_data10_type_c
typedef fcs_front_xq_aI_sl_data10_type_c fcs_front_xq_aI_sldata10_t;
#endif
/* sl_type fcs_front_xq_aI_sldata11_t */
#ifdef fcs_front_xq_aI_sl_data11_type_c
typedef fcs_front_xq_aI_sl_data11_type_c fcs_front_xq_aI_sldata11_t;
#endif
/* sl_type fcs_front_xq_aI_sldata12_t */
#ifdef fcs_front_xq_aI_sl_data12_type_c
typedef fcs_front_xq_aI_sl_data12_type_c fcs_front_xq_aI_sldata12_t;
#endif
/* sl_type fcs_front_xq_aI_sldata13_t */
#ifdef fcs_front_xq_aI_sl_data13_type_c
typedef fcs_front_xq_aI_sl_data13_type_c fcs_front_xq_aI_sldata13_t;
#endif
/* sl_type fcs_front_xq_aI_sldata14_t */
#ifdef fcs_front_xq_aI_sl_data14_type_c
typedef fcs_front_xq_aI_sl_data14_type_c fcs_front_xq_aI_sldata14_t;
#endif
/* sl_type fcs_front_xq_aI_sldata15_t */
#ifdef fcs_front_xq_aI_sl_data15_type_c
typedef fcs_front_xq_aI_sl_data15_type_c fcs_front_xq_aI_sldata15_t;
#endif
/* sl_type fcs_front_xq_aI_sldata16_t */
#ifdef fcs_front_xq_aI_sl_data16_type_c
typedef fcs_front_xq_aI_sl_data16_type_c fcs_front_xq_aI_sldata16_t;
#endif
/* sl_type fcs_front_xq_aI_sldata17_t */
#ifdef fcs_front_xq_aI_sl_data17_type_c
typedef fcs_front_xq_aI_sl_data17_type_c fcs_front_xq_aI_sldata17_t;
#endif
/* sl_type fcs_front_xq_aI_sldata18_t */
#ifdef fcs_front_xq_aI_sl_data18_type_c
typedef fcs_front_xq_aI_sl_data18_type_c fcs_front_xq_aI_sldata18_t;
#endif
/* sl_type fcs_front_xq_aI_sldata19_t */
#ifdef fcs_front_xq_aI_sl_data19_type_c
typedef fcs_front_xq_aI_sl_data19_type_c fcs_front_xq_aI_sldata19_t;
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

/* sl_type fcs_front_xq_aI_slweight_t */
typedef fcs_front_xq_aI_sl_weight_type_c fcs_front_xq_aI_slweight_t;

#define fcs_front_xq_aI_slweight_fmt  fcs_front_xq_aI_sl_weight_type_fmt  /* sl_macro */

#if defined(fcs_front_xq_aI_sl_elem_weight) && defined(fcs_front_xq_aI_sl_weight_intequiv)
typedef fcs_front_xq_aI_sl_weight_type_c fcs_front_xq_aI_slcount_t;       /* sl_type fcs_front_xq_aI_slcount_t */
# define fcs_front_xq_aI_slcount_fmt  fcs_front_xq_aI_sl_weight_type_fmt  /* sl_macro */
#else
typedef fcs_front_xq_aI_sl_int_type_c fcs_front_xq_aI_slcount_t;
# define fcs_front_xq_aI_slcount_fmt  fcs_front_xq_aI_sl_int_type_fmt
#endif


/* sl_type fcs_front_xq_aI__slpwkey_t fcs_front_xq_aI_slpwkey_t */
typedef struct fcs_front_xq_aI__slpwkey_t
{
  fcs_front_xq_aI_slpkey_t pkey;
  fcs_front_xq_aI_slweight_t weight;

} fcs_front_xq_aI_slpwkey_t;


/* sl_type fcs_front_xq_aI__elements_t fcs_front_xq_aI_elements_t */
typedef struct fcs_front_xq_aI__elements_t
{
  fcs_front_xq_aI_slint_t size, max_size;
  fcs_front_xq_aI_slkey_t *keys;

#ifdef fcs_front_xq_aI_SL_INDEX
  fcs_front_xq_aI_slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef fcs_front_xq_aI_SL_DATA0
  fcs_front_xq_aI_sldata0_t *data0;
#endif
#ifdef fcs_front_xq_aI_SL_DATA1
  fcs_front_xq_aI_sldata1_t *data1;
#endif
#ifdef fcs_front_xq_aI_SL_DATA2
  fcs_front_xq_aI_sldata2_t *data2;
#endif
#ifdef fcs_front_xq_aI_SL_DATA3
  fcs_front_xq_aI_sldata3_t *data3;
#endif
#ifdef fcs_front_xq_aI_SL_DATA4
  fcs_front_xq_aI_sldata4_t *data4;
#endif
#ifdef fcs_front_xq_aI_SL_DATA5
  fcs_front_xq_aI_sldata5_t *data5;
#endif
#ifdef fcs_front_xq_aI_SL_DATA6
  fcs_front_xq_aI_sldata6_t *data6;
#endif
#ifdef fcs_front_xq_aI_SL_DATA7
  fcs_front_xq_aI_sldata7_t *data7;
#endif
#ifdef fcs_front_xq_aI_SL_DATA8
  fcs_front_xq_aI_sldata8_t *data8;
#endif
#ifdef fcs_front_xq_aI_SL_DATA9
  fcs_front_xq_aI_sldata9_t *data9;
#endif
#ifdef fcs_front_xq_aI_SL_DATA10
  fcs_front_xq_aI_sldata10_t *data10;
#endif
#ifdef fcs_front_xq_aI_SL_DATA11
  fcs_front_xq_aI_sldata11_t *data11;
#endif
#ifdef fcs_front_xq_aI_SL_DATA12
  fcs_front_xq_aI_sldata12_t *data12;
#endif
#ifdef fcs_front_xq_aI_SL_DATA13
  fcs_front_xq_aI_sldata13_t *data13;
#endif
#ifdef fcs_front_xq_aI_SL_DATA14
  fcs_front_xq_aI_sldata14_t *data14;
#endif
#ifdef fcs_front_xq_aI_SL_DATA15
  fcs_front_xq_aI_sldata15_t *data15;
#endif
#ifdef fcs_front_xq_aI_SL_DATA16
  fcs_front_xq_aI_sldata16_t *data16;
#endif
#ifdef fcs_front_xq_aI_SL_DATA17
  fcs_front_xq_aI_sldata17_t *data17;
#endif
#ifdef fcs_front_xq_aI_SL_DATA18
  fcs_front_xq_aI_sldata18_t *data18;
#endif
#ifdef fcs_front_xq_aI_SL_DATA19
  fcs_front_xq_aI_sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} fcs_front_xq_aI_elements_t;


/* sl_type fcs_front_xq_aI__packed_element_t fcs_front_xq_aI_packed_element_t */
typedef struct fcs_front_xq_aI__packed_element_t
{
  fcs_front_xq_aI_slkey_t key;

#ifdef fcs_front_xq_aI_SL_PACKED_INDEX
  fcs_front_xq_aI_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef fcs_front_xq_aI_SL_DATA0
# ifdef fcs_front_xq_aI_sl_data0_flex
  fcs_front_xq_aI_sldata0_t data0[];
# else
  fcs_front_xq_aI_sldata0_t data0[fcs_front_xq_aI_sl_data0_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA1
# ifdef fcs_front_xq_aI_sl_data1_flex
  fcs_front_xq_aI_sldata1_t data1[];
# else
  fcs_front_xq_aI_sldata1_t data1[fcs_front_xq_aI_sl_data1_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA2
# ifdef fcs_front_xq_aI_sl_data2_flex
  fcs_front_xq_aI_sldata2_t data2[];
# else
  fcs_front_xq_aI_sldata2_t data2[fcs_front_xq_aI_sl_data2_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA3
# ifdef fcs_front_xq_aI_sl_data3_flex
  fcs_front_xq_aI_sldata3_t data3[];
# else
  fcs_front_xq_aI_sldata3_t data3[fcs_front_xq_aI_sl_data3_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA4
# ifdef fcs_front_xq_aI_sl_data4_flex
  fcs_front_xq_aI_sldata4_t data4[];
# else
  fcs_front_xq_aI_sldata4_t data4[fcs_front_xq_aI_sl_data4_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA5
# ifdef fcs_front_xq_aI_sl_data5_flex
  fcs_front_xq_aI_sldata5_t data5[];
# else
  fcs_front_xq_aI_sldata5_t data5[fcs_front_xq_aI_sl_data5_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA6
# ifdef fcs_front_xq_aI_sl_data6_flex
  fcs_front_xq_aI_sldata6_t data6[];
# else
  fcs_front_xq_aI_sldata6_t data6[fcs_front_xq_aI_sl_data6_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA7
# ifdef fcs_front_xq_aI_sl_data7_flex
  fcs_front_xq_aI_sldata7_t data7[];
# else
  fcs_front_xq_aI_sldata7_t data7[fcs_front_xq_aI_sl_data7_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA8
# ifdef fcs_front_xq_aI_sl_data8_flex
  fcs_front_xq_aI_sldata8_t data8[];
# else
  fcs_front_xq_aI_sldata8_t data8[fcs_front_xq_aI_sl_data8_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA9
# ifdef fcs_front_xq_aI_sl_data9_flex
  fcs_front_xq_aI_sldata9_t data9[];
# else
  fcs_front_xq_aI_sldata9_t data9[fcs_front_xq_aI_sl_data9_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA10
# ifdef fcs_front_xq_aI_sl_data10_flex
  fcs_front_xq_aI_sldata10_t data10[];
# else
  fcs_front_xq_aI_sldata10_t data10[fcs_front_xq_aI_sl_data10_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA11
# ifdef fcs_front_xq_aI_sl_data11_flex
  fcs_front_xq_aI_sldata11_t data11[];
# else
  fcs_front_xq_aI_sldata11_t data11[fcs_front_xq_aI_sl_data11_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA12
# ifdef fcs_front_xq_aI_sl_data12_flex
  fcs_front_xq_aI_sldata12_t data12[];
# else
  fcs_front_xq_aI_sldata12_t data12[fcs_front_xq_aI_sl_data12_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA13
# ifdef fcs_front_xq_aI_sl_data13_flex
  fcs_front_xq_aI_sldata13_t data13[];
# else
  fcs_front_xq_aI_sldata13_t data13[fcs_front_xq_aI_sl_data13_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA14
# ifdef fcs_front_xq_aI_sl_data14_flex
  fcs_front_xq_aI_sldata14_t data14[];
# else
  fcs_front_xq_aI_sldata14_t data14[fcs_front_xq_aI_sl_data14_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA15
# ifdef fcs_front_xq_aI_sl_data15_flex
  fcs_front_xq_aI_sldata15_t data15[];
# else
  fcs_front_xq_aI_sldata15_t data15[fcs_front_xq_aI_sl_data15_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA16
# ifdef fcs_front_xq_aI_sl_data16_flex
  fcs_front_xq_aI_sldata16_t data16[];
# else
  fcs_front_xq_aI_sldata16_t data16[fcs_front_xq_aI_sl_data16_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA17
# ifdef fcs_front_xq_aI_sl_data17_flex
  fcs_front_xq_aI_sldata17_t data17[];
# else
  fcs_front_xq_aI_sldata17_t data17[fcs_front_xq_aI_sl_data17_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA18
# ifdef fcs_front_xq_aI_sl_data18_flex
  fcs_front_xq_aI_sldata18_t data18[];
# else
  fcs_front_xq_aI_sldata18_t data18[fcs_front_xq_aI_sl_data18_size_c];
# endif
#endif
#ifdef fcs_front_xq_aI_SL_DATA19
# ifdef fcs_front_xq_aI_sl_data19_flex
  fcs_front_xq_aI_sldata19_t data19[];
# else
  fcs_front_xq_aI_sldata19_t data19[fcs_front_xq_aI_sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} fcs_front_xq_aI_packed_element_t;


/* sl_type fcs_front_xq_aI__packed_elements_t fcs_front_xq_aI_packed_elements_t */
typedef struct fcs_front_xq_aI__packed_elements_t
{
  fcs_front_xq_aI_slint_t size, max_size;
  
  fcs_front_xq_aI_packed_element_t *elements;
  
} fcs_front_xq_aI_packed_elements_t;


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


/* sl_type fcs_front_xq_aI__classification_info_t fcs_front_xq_aI_classification_info_t fcs_front_xq_aI_classification_info */
typedef struct fcs_front_xq_aI__classification_info_t
{
  fcs_front_xq_aI_slint_t nclasses;
  fcs_front_xq_aI_slkey_pure_t *keys;
  fcs_front_xq_aI_slint_t *counts;
  fcs_front_xq_aI_slint_t *masks;

  /* */
  fcs_front_xq_aI_slint_t *all_local_sizes;
  fcs_front_xq_aI_slint_t *local_lt_eq_counts;
  fcs_front_xq_aI_slint_t *all_local_lt_eq_counts;

} fcs_front_xq_aI_classification_info_t, fcs_front_xq_aI_classification_info;  /* deprecated 'fcs_front_xq_aI_classification_info' */


/* key2class, sl_type fcs_front_xq_aI_key2class_f */
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_key2class_f)(fcs_front_xq_aI_slkey_t *, fcs_front_xq_aI_slint, void *);

/* pivot-element, sl_type fcs_front_xq_aI_pivot_f */
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_pivot_f)(fcs_front_xq_aI_elements_t *);

/* sorting-network, sl_type fcs_front_xq_aI_sortnet_f fcs_front_xq_aI_sortnet_data_t */
typedef void *fcs_front_xq_aI_sortnet_data_t;
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_sortnet_f)(fcs_front_xq_aI_slint_t size, fcs_front_xq_aI_slint_t rank, fcs_front_xq_aI_slint_t stage, fcs_front_xq_aI_sortnet_data_t snd, fcs_front_xq_aI_slint_t *up);

/* merge2, sl_type fcs_front_xq_aI_merge2x_f fcs_front_xq_aI_merge2X_f */
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_merge2x_f)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_merge2X_f)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_elements_t *t);

/* sl_type fcs_front_xq_aI__permute_generic_t fcs_front_xq_aI_permute_generic_t */
typedef struct fcs_front_xq_aI__permute_generic_t
{
  int type;

  fcs_front_xq_aI_spec_tloc_f *tloc;
  fcs_front_xq_aI_spec_tloc_rearrange_db_f *tloc_rearrange_db;
  fcs_front_xq_aI_spec_tloc_rearrange_ip_f *tloc_rearrange_ip;

  fcs_front_xq_aI_spec_tloc_mod_f *tloc_mod;
  fcs_front_xq_aI_spec_tloc_mod_rearrange_db_f *tloc_mod_rearrange_db;
  fcs_front_xq_aI_spec_tloc_mod_rearrange_ip_f *tloc_mod_rearrange_ip;

} fcs_front_xq_aI_permute_generic_t;

/* sl_macro fcs_front_xq_aI_PERMUTE_GENERIC_DEFINE_TLOC fcs_front_xq_aI_PERMUTE_GENERIC_INIT_TLOC fcs_front_xq_aI_PERMUTE_GENERIC_INIT_EXT_TLOC */
#define fcs_front_xq_aI_PERMUTE_GENERIC_DEFINE_TLOC(_tl_, _s_...)      fcs_front_xq_aI_SPEC_DEFINE_TLOC(_tl_, _tl_, _s_)
#define fcs_front_xq_aI_PERMUTE_GENERIC_INIT_TLOC(_tl_)                { 1, _tl_, fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_NULL,  NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define fcs_front_xq_aI_PERMUTE_GENERIC_INIT_EXT_TLOC(_tl_)            { 1, _tl_, fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC(_tl_), NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD_NULL }

/* sl_macro fcs_front_xq_aI_PERMUTE_GENERIC_DEFINE_TLOC_MOD fcs_front_xq_aI_PERMUTE_GENERIC_INIT_TLOC_MOD fcs_front_xq_aI_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD */
#define fcs_front_xq_aI_PERMUTE_GENERIC_DEFINE_TLOC_MOD(_tl_, _s_...)  fcs_front_xq_aI_SPEC_DEFINE_TLOC_MOD(_tl_, _tl_, _s_)
#define fcs_front_xq_aI_PERMUTE_GENERIC_INIT_TLOC_MOD(_tl_)            { 2, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define fcs_front_xq_aI_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD(_tl_)        { 2, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, fcs_front_xq_aI_SPEC_EXT_PARAM_TLOC_MOD(_tl_) }

/* sl_type fcs_front_xq_aI__split_generic_t fcs_front_xq_aI_split_generic_t */
typedef struct fcs_front_xq_aI__split_generic_t
{
  int type;

  fcs_front_xq_aI_slint_t max_tprocs;

  fcs_front_xq_aI_spec_tproc_f *tproc;
  fcs_front_xq_aI_spec_tproc_ext_t tproc_ext;

  fcs_front_xq_aI_spec_tproc_mod_f *tproc_mod;
  fcs_front_xq_aI_spec_tproc_mod_ext_t tproc_mod_ext;

  fcs_front_xq_aI_spec_tprocs_f *tprocs;
  fcs_front_xq_aI_spec_tprocs_ext_t tprocs_ext;

  fcs_front_xq_aI_spec_tprocs_mod_f *tprocs_mod;
  fcs_front_xq_aI_spec_tprocs_mod_ext_t tprocs_mod_ext;

  fcs_front_xq_aI_spec_tproc_reset_f *reset;

} fcs_front_xq_aI_split_generic_t;

/* sl_macro fcs_front_xq_aI_SPLIT_GENERIC_DEFINE_TPROC fcs_front_xq_aI_SPLIT_GENERIC_INIT_TPROC fcs_front_xq_aI_SPLIT_GENERIC_INIT_EXT_TPROC */
#define fcs_front_xq_aI_SPLIT_GENERIC_DEFINE_TPROC(_tp_, _s_...)                fcs_front_xq_aI_SPEC_DEFINE_TPROC(_tp_, _tp_, _s_)
#define fcs_front_xq_aI_SPLIT_GENERIC_INIT_TPROC(_tp_, _r_...)                  { 1, 0, _tp_, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL,  NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }
#define fcs_front_xq_aI_SPLIT_GENERIC_INIT_EXT_TPROC(_tp_, _r_...)              { 1, 0, _tp_, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC(_tp_), NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }

/* sl_macro fcs_front_xq_aI_SPLIT_GENERIC_DEFINE_TPROC_MOD fcs_front_xq_aI_SPLIT_GENERIC_INIT_TPROC_MOD fcs_front_xq_aI_SPLIT_GENERIC_INIT_EXT_TPROC_MOD */
#define fcs_front_xq_aI_SPLIT_GENERIC_DEFINE_TPROC_MOD(_tp_, _s_...)            fcs_front_xq_aI_SPEC_DEFINE_TPROC_MOD(_tp_, _tp_, _s_)
#define fcs_front_xq_aI_SPLIT_GENERIC_INIT_TPROC_MOD(_tp_, _r_...)              { 2, 0, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL, _tp_, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL,  NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }
#define fcs_front_xq_aI_SPLIT_GENERIC_INIT_EXT_TPROC_MOD(_tp_, _r_...)          { 2, 0, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL, _tp_, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD(_tp_), NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }

/* sl_macro fcs_front_xq_aI_SPLIT_GENERIC_DEFINE_TPROCS fcs_front_xq_aI_SPLIT_GENERIC_INIT_TPROCS fcs_front_xq_aI_SPLIT_GENERIC_INIT_EXT_TPROCS */
#define fcs_front_xq_aI_SPLIT_GENERIC_DEFINE_TPROCS(_tp_, _s_...)               fcs_front_xq_aI_SPEC_DEFINE_TPROCS(_tp_, _tp_, _s_)
#define fcs_front_xq_aI_SPLIT_GENERIC_INIT_TPROCS(_tp_, _xtp_, _r_...)          { 3, (_xtp_), NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL,  NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }
#define fcs_front_xq_aI_SPLIT_GENERIC_INIT_EXT_TPROCS(_tp_, _xtp_, _r_...)      { 3, (_xtp_), NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS(_tp_), NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL, _r_ }

/* sl_macro fcs_front_xq_aI_SPLIT_GENERIC_DEFINE_TPROCS_MOD fcs_front_xq_aI_SPLIT_GENERIC_INIT_TPROCS_MOD fcs_front_xq_aI_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD */
#define fcs_front_xq_aI_SPLIT_GENERIC_DEFINE_TPROCS_MOD(_tp_, _s_...)           fcs_front_xq_aI_SPEC_DEFINE_TPROCS_MOD(_tp_, _tp_, _s_)
#define fcs_front_xq_aI_SPLIT_GENERIC_INIT_TPROCS_MOD(_tp_, _xtp_, _r_...)      { 4, (_xtp_), NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL, _tp_, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define fcs_front_xq_aI_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD(_tp_, _xtp_, _r_...)  { 4, (_xtp_), NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL, _tp_, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD(_tp_), _r_ }

/* sl_type fcs_front_xq_aI_tloc_f fcs_front_xq_aI_tloc_mod_f */
typedef fcs_front_xq_aI_slint_t fcs_front_xq_aI_tloc_f(fcs_front_xq_aI_elements_t *b, fcs_front_xq_aI_slint_t x, void *tloc_data);
typedef fcs_front_xq_aI_slint_t fcs_front_xq_aI_tloc_mod_f(fcs_front_xq_aI_elements_t *b, fcs_front_xq_aI_slint_t x, void *tloc_data, fcs_front_xq_aI_elements_t *mod);

/* sl_type fcs_front_xq_aI_tproc_f fcs_front_xq_aI_tproc_mod_f fcs_front_xq_aI_tprocs_f fcs_front_xq_aI_tprocs_mod_f */
typedef int fcs_front_xq_aI_tproc_f(fcs_front_xq_aI_elements_t *b, fcs_front_xq_aI_slint_t x, void *tproc_data);
typedef int fcs_front_xq_aI_tproc_mod_f(fcs_front_xq_aI_elements_t *b, fcs_front_xq_aI_slint_t x, void *tproc_data, fcs_front_xq_aI_elements_t *mod);
typedef void fcs_front_xq_aI_tprocs_f(fcs_front_xq_aI_elements_t *b, fcs_front_xq_aI_slint_t x, void *tproc_data, fcs_front_xq_aI_slint_t *nprocs, int *procs);
typedef void fcs_front_xq_aI_tprocs_mod_f(fcs_front_xq_aI_elements_t *b, fcs_front_xq_aI_slint_t x, void *tproc_data, fcs_front_xq_aI_slint_t *nprocs, int *procs, fcs_front_xq_aI_elements_t *mod);

/* sl_type fcs_front_xq_aI_tproc_reset_f */
typedef void fcs_front_xq_aI_tproc_reset_f(void *tproc_data);

/* sl_macro fcs_front_xq_aI_TPROC_RESET_NULL */
#define fcs_front_xq_aI_TPROC_RESET_NULL  NULL

/* sl_type fcs_front_xq_aI_tproc_t */
typedef struct fcs_front_xq_aI__spec_tproc_t *fcs_front_xq_aI_tproc_t;

/* sl_type fcs_front_xq_aI__tproc_exdef fcs_front_xq_aI_tproc_exdef */
typedef struct fcs_front_xq_aI__tproc_exdef
{
  int type;

  fcs_front_xq_aI_spec_tproc_ext_t tproc_ext;
  fcs_front_xq_aI_spec_tproc_mod_ext_t tproc_mod_ext;
  fcs_front_xq_aI_spec_tprocs_ext_t tprocs_ext;
  fcs_front_xq_aI_spec_tprocs_mod_ext_t tprocs_mod_ext;

} const *fcs_front_xq_aI_tproc_exdef;

/* sl_macro fcs_front_xq_aI_TPROC_EXDEF_NULL */
#define fcs_front_xq_aI_TPROC_EXDEF_NULL  NULL

/* sl_macro fcs_front_xq_aI_TPROC_EXDEF_DEFINE_TPROC fcs_front_xq_aI_TPROC_EXDEF_DEFINE_TPROC_MOD fcs_front_xq_aI_TPROC_EXDEF_DEFINE_TPROCS fcs_front_xq_aI_TPROC_EXDEF_DEFINE_TPROCS_MOD */
#define fcs_front_xq_aI_TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct fcs_front_xq_aI__tproc_exdef _##_name_ = { 1, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC(_name_), fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define fcs_front_xq_aI_TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct fcs_front_xq_aI__tproc_exdef _##_name_ = { 2, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD(_name_), fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define fcs_front_xq_aI_TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct fcs_front_xq_aI__tproc_exdef _##_name_ = { 3, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS(_name_), fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define fcs_front_xq_aI_TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  fcs_front_xq_aI_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct fcs_front_xq_aI__tproc_exdef _##_name_ = { 4, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROC_MOD_NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_NULL, fcs_front_xq_aI_SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;


/* fcs_front_xq_aI_mpi_elements_alltoall_specific */
#ifndef SL_MEAS_TYPE_ALLTOALLV
# define SL_MEAS_TYPE_ALLTOALLV  0
#endif

#ifndef SL_MEAS_TYPE_SENDRECV
# define SL_MEAS_TYPE_SENDRECV   1
#endif


/* deprecated, sl_type fcs_front_xq_aI_k2c_func fcs_front_xq_aI_pivot_func fcs_front_xq_aI_sn_func fcs_front_xq_aI_m2x_func fcs_front_xq_aI_m2X_func */
typedef fcs_front_xq_aI_key2class_f fcs_front_xq_aI_k2c_func;
typedef fcs_front_xq_aI_pivot_f fcs_front_xq_aI_pivot_func;
typedef fcs_front_xq_aI_sortnet_f fcs_front_xq_aI_sn_func;
typedef fcs_front_xq_aI_merge2x_f fcs_front_xq_aI_m2x_func;
typedef fcs_front_xq_aI_merge2X_f fcs_front_xq_aI_m2X_func;


/* sl_type fcs_front_xq_aI__mergek_t fcs_front_xq_aI_mergek_t */
typedef struct fcs_front_xq_aI__mergek_t
{
  fcs_front_xq_aI_sortnet_f sn;
  fcs_front_xq_aI_sortnet_data_t snd;

  fcs_front_xq_aI_merge2x_f m2x;
  fcs_front_xq_aI_elements_t *sx;

} fcs_front_xq_aI_mergek_t;


/* sl_type fcs_front_xq_aI_keys_init_type_t fcs_front_xq_aI_keys_init_data_t */
typedef fcs_front_xq_aI_slint_t fcs_front_xq_aI_keys_init_type_t;
typedef void *fcs_front_xq_aI_keys_init_data_t;

/* sl_type fcs_front_xq_aI_key_set_data_t fcs_front_xq_aI_key_set_f */
typedef void *fcs_front_xq_aI_key_set_data_t;
typedef void (*fcs_front_xq_aI_key_set_f)(fcs_front_xq_aI_slkey_pure_t *k, fcs_front_xq_aI_key_set_data_t d);


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


/* fcs_front_xq_aI_elements_keys_stats */
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


/* partition conditions, sl_type fcs_front_xq_aI__partcond2_t fcs_front_xq_aI_partcond2_t */
typedef struct fcs_front_xq_aI__partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} fcs_front_xq_aI_partcond2_t;


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

/* partition conditions, sl_type fcs_front_xq_aI__partcond_t fcs_front_xq_aI_partcond_t fcs_front_xq_aI_partcond_p */
typedef struct fcs_front_xq_aI__partcond_t
{
  fcs_front_xq_aI_slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} fcs_front_xq_aI_partcond_t, *fcs_front_xq_aI_partcond_p;


/* internal partition conditions, sl_type fcs_front_xq_aI__partcond_intern_t fcs_front_xq_aI_partcond_intern_t fcs_front_xq_aI_partcond_intern_p */
typedef struct fcs_front_xq_aI__partcond_intern_t
{
  fcs_front_xq_aI_slint_t pcm;
  fcs_front_xq_aI_slint_t count_min, count_max;
  fcs_front_xq_aI_slint_t count_low, count_high;
#ifdef fcs_front_xq_aI_elem_weight
  fcs_front_xq_aI_slweight_t weight_min, weight_max;
  fcs_front_xq_aI_slweight_t weight_low, weight_high;
#endif

} fcs_front_xq_aI_partcond_intern_t, *fcs_front_xq_aI_partcond_intern_p;


/* sl_type fcs_front_xq_aI__parttype_t fcs_front_xq_aI_parttype_t fcs_front_xq_aI_parttype_p */
typedef struct fcs_front_xq_aI__parttype_t
{
  fcs_front_xq_aI_slint_t type;

} fcs_front_xq_aI_parttype_t, *fcs_front_xq_aI_parttype_p;


/* generic binning method */

/* sl_type fcs_front_xq_aI__bin_t fcs_front_xq_aI_bin_t */
typedef struct fcs_front_xq_aI__bin_t
{
  fcs_front_xq_aI_elements_t s;

#ifdef fcs_front_xq_aI_elem_weight
  fcs_front_xq_aI_slweight_t weight;
#endif

} fcs_front_xq_aI_bin_t;


/* sl_type fcs_front_xq_aI__splitter_t fcs_front_xq_aI_splitter_t */
typedef struct fcs_front_xq_aI__splitter_t
{
  fcs_front_xq_aI_slint_t n;

  int *displs;
  fcs_front_xq_aI_slkey_pure_t *s;
  fcs_front_xq_aI_slint_t *sn;

} fcs_front_xq_aI_splitter_t;


struct fcs_front_xq_aI__binning_t;

/* sl_type fcs_front_xq_aI_binning_pre_f fcs_front_xq_aI_binning_exec_f fcs_front_xq_aI_binning_refine_f fcs_front_xq_aI_binning_hit_f fcs_front_xq_aI_binning_finalize_f fcs_front_xq_aI_binning_post_f */
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_binning_pre_f)(struct fcs_front_xq_aI__binning_t *bm);
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_binning_exec_f)(struct fcs_front_xq_aI__binning_t *bm, fcs_front_xq_aI_bin_t *bin, fcs_front_xq_aI_slcount_t *counts, fcs_front_xq_aI_slweight_t *weights);
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_binning_refine_f)(struct fcs_front_xq_aI__binning_t *bm, fcs_front_xq_aI_bin_t *bin, fcs_front_xq_aI_slint_t k, fcs_front_xq_aI_slcount_t *counts, fcs_front_xq_aI_slweight_t *weights, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s, fcs_front_xq_aI_bin_t *new_bin);
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_binning_hit_f)(struct fcs_front_xq_aI__binning_t *bm, fcs_front_xq_aI_bin_t *bin, fcs_front_xq_aI_slint_t k, fcs_front_xq_aI_slcount_t *counts, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s);
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_binning_finalize_f)(struct fcs_front_xq_aI__binning_t *bm, fcs_front_xq_aI_bin_t *bin, fcs_front_xq_aI_slint_t dc, fcs_front_xq_aI_slweight_t dw, fcs_front_xq_aI_slint_t lc_min, fcs_front_xq_aI_slint_t lc_max, fcs_front_xq_aI_slcount_t *lcs, fcs_front_xq_aI_slweight_t *lws, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s);
typedef fcs_front_xq_aI_slint_t (*fcs_front_xq_aI_binning_post_f)(struct fcs_front_xq_aI__binning_t *bm);


/* sl_type fcs_front_xq_aI__binning_data_t fcs_front_xq_aI_binning_data_t */
typedef union fcs_front_xq_aI__binning_data_t
{
  struct
  {
    fcs_front_xq_aI_slint_t rhigh, rlow, rwidth;
    fcs_front_xq_aI_slint_t rcurrent;
    fcs_front_xq_aI_slkey_pure_t bit_mask;

    fcs_front_xq_aI_elements_t sx;

  } radix;

} fcs_front_xq_aI_binning_data_t;


/* sl_type fcs_front_xq_aI__binning_t fcs_front_xq_aI_binning_t */
typedef struct fcs_front_xq_aI__binning_t
{
  fcs_front_xq_aI_slint_t nbins, max_nbins;
  
  fcs_front_xq_aI_binning_pre_f pre;
  fcs_front_xq_aI_binning_exec_f exec;
  fcs_front_xq_aI_binning_refine_f refine;
  fcs_front_xq_aI_binning_hit_f hit;
  fcs_front_xq_aI_binning_finalize_f finalize;
  fcs_front_xq_aI_binning_post_f post;

  fcs_front_xq_aI_slint_t sorted;

  fcs_front_xq_aI_slint_t docounts;
#ifdef fcs_front_xq_aI_elem_weight
  fcs_front_xq_aI_slint_t doweights;
#endif

  fcs_front_xq_aI_binning_data_t bd;

} fcs_front_xq_aI_binning_t;


/* sl_type fcs_front_xq_aI__local_bins_t fcs_front_xq_aI_local_bins_t */
typedef struct fcs_front_xq_aI__local_bins_t
{
  fcs_front_xq_aI_binning_t *bm;

  fcs_front_xq_aI_slint_t nbins, max_nbins;
  fcs_front_xq_aI_slint_t nelements;

  fcs_front_xq_aI_slint_t docounts;
#ifdef fcs_front_xq_aI_elem_weight
  fcs_front_xq_aI_slint_t doweights;
#endif

  fcs_front_xq_aI_slint_t nbinnings, max_nbinnings;

  fcs_front_xq_aI_slint_t nbins_new, last_new_b, last_new_k;
  fcs_front_xq_aI_bin_t *bins, *bins_new;
  fcs_front_xq_aI_bin_t *bins0, *bins1;

  fcs_front_xq_aI_slint_t *bcws;

#if defined(fcs_front_xq_aI_elem_weight) && defined(fcs_front_xq_aI_sl_weight_intequiv)
  fcs_front_xq_aI_slint_t cw_factor, w_index, bin_cw_factor;
  fcs_front_xq_aI_slweight_t *cws, *bin_cws;
  fcs_front_xq_aI_slweight_t *prefix_cws;
#else
  fcs_front_xq_aI_slint_t *cs, *bin_cs;
  fcs_front_xq_aI_slint_t *prefix_cs;
# ifdef fcs_front_xq_aI_elem_weight
  fcs_front_xq_aI_slweight_t *ws, *bin_ws;
  fcs_front_xq_aI_slweight_t *prefix_ws;
# endif
#endif

  fcs_front_xq_aI_slint_t last_exec_b;

} fcs_front_xq_aI_local_bins_t;


/* sl_type fcs_front_xq_aI__global_bins_t fcs_front_xq_aI_global_bins_t */
typedef struct fcs_front_xq_aI__global_bins_t
{
  fcs_front_xq_aI_binning_t *bm;
  
  fcs_front_xq_aI_local_bins_t lb;

  fcs_front_xq_aI_slint_t *bcws;

#if defined(fcs_front_xq_aI_elem_weight) && defined(fcs_front_xq_aI_sl_weight_intequiv)
  fcs_front_xq_aI_slweight_t *cws;
  fcs_front_xq_aI_slweight_t *prefix_cws;
#else
  fcs_front_xq_aI_slint_t *cs;
  fcs_front_xq_aI_slint_t *prefix_cs;
# ifdef fcs_front_xq_aI_elem_weight
  fcs_front_xq_aI_slweight_t *ws;
  fcs_front_xq_aI_slweight_t *prefix_ws;
# endif
#endif

} fcs_front_xq_aI_global_bins_t;


/* sl_type fcs_front_xq_aI_rti_cmc_t */
typedef struct
{
  fcs_front_xq_aI_slint_t cmp, movek, moved;

} fcs_front_xq_aI_rti_cmc_t;

#ifndef my_rti_ccmp
# define my_rti_ccmp(m)    m.cmc.cmp
# define my_rti_cmovek(m)  m.cmc.movek
# define my_rti_cmoved(m)  m.cmc.moved
#endif


/* sl_type fcs_front_xq_aI_rti_tim_t */
typedef struct
{
  double start, stop;
  double last, cumu;

  fcs_front_xq_aI_slint_t num;

} fcs_front_xq_aI_rti_tim_t[rti_tids];

#ifndef my_rti_tlast
# define my_rti_tlast(m, t)  m.tim[t].last
# define my_rti_tcumu(m, t)  m.tim[t].cumu
# define my_rti_tnum(m, t)   m.tim[t].num
#endif


/* sl_type fcs_front_xq_aI_rti_mem_t */
typedef struct
{
  fcs_front_xq_aI_slint_t nalloc, nfree;
  fcs_front_xq_aI_slint_t max, cur, cur_max;

} fcs_front_xq_aI_rti_mem_t;


/* sl_type fcs_front_xq_aI_rti_t */
typedef struct
{
  /* compare-move-counter */
  fcs_front_xq_aI_rti_cmc_t cmc;
  /* timer */
  fcs_front_xq_aI_rti_tim_t tim;
  /* memory */
  fcs_front_xq_aI_rti_mem_t mem;

} fcs_front_xq_aI_rti_t;

#ifndef my_rti_reset
# define my_rti_reset(m)  memset((void *) &m, 0, sizeof(m))
#endif


#ifdef SL_USE_MPI
/* sl_type fcs_front_xq_aI__sl_mpi_context_t fcs_front_xq_aI_sl_mpi_context_t */
typedef struct fcs_front_xq_aI__sl_mpi_context_t
{
  int size, rank;
  MPI_Comm comm;

} fcs_front_xq_aI_sl_mpi_context_t;
#endif


#ifdef SL_USE_OMP
/* sl_type fcs_front_xq_aI__sl_omp_context_t fcs_front_xq_aI_sl_omp_context_t */
typedef struct fcs_front_xq_aI__sl_omp_context_t
{
  int thread_num, num_threads, *coop_thread_nums;

} fcs_front_xq_aI_sl_omp_context_t;
#endif


/* sl_type fcs_front_xq_aI__sl_context_t fcs_front_xq_aI_sl_context_t */
typedef struct fcs_front_xq_aI__sl_context_t
{

/* base/base.c */
  struct {
int dummy_rank;
  } sl;
#ifdef fcs_front_xq_aI_SL_USE_RTI
fcs_front_xq_aI_rti_t rti;
#endif
  struct {
fcs_front_xq_aI_slint_t ip_threshold;
fcs_front_xq_aI_slint_t db_threshold;
fcs_front_xq_aI_slint_t ma_threshold;
  } sr;
  struct {
fcs_front_xq_aI_slint_t threshold;
  } sri;
/* base_mpi/base_mpi.c */
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
fcs_front_xq_aI_slint_t sendrecv_replace_memsize;
fcs_front_xq_aI_slint_t sendrecv_replace_mpi_maxsize;
  } me;
#endif
#ifdef SL_USE_MPI
  struct {
double t[2];
fcs_front_xq_aI_slint_t packed;
fcs_front_xq_aI_slint_t minalloc;
double overalloc;
fcs_front_xq_aI_slint_t type;
void *sendrecv_aux;
fcs_front_xq_aI_slint_t sendrecv_aux_size;
fcs_front_xq_aI_slint_t sendrecv_requests;
  } meas;
#endif
#ifdef SL_USE_MPI
  struct {
fcs_front_xq_aI_slint_t packed;
fcs_front_xq_aI_slint_t db_packed;
fcs_front_xq_aI_slint_t ip_packed;
  } mea;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef fcs_front_xq_aI_MSEG_ROOT
int root;
#endif
#ifdef fcs_front_xq_aI_MSEG_BORDER_UPDATE_REDUCTION
double border_update_count_reduction;
double border_update_weight_reduction;
#endif
#ifdef fcs_front_xq_aI_MSEG_FORWARD_ONLY
fcs_front_xq_aI_slint_t forward_only;
#endif
#ifdef fcs_front_xq_aI_MSEG_INFO
fcs_front_xq_aI_slint_t info_rounds;
fcs_front_xq_aI_slint_t *info_finish_rounds;
double info_finish_rounds_avg;
fcs_front_xq_aI_slweight_t info_total_weights;
#endif
fcs_front_xq_aI_slint_t binnings;
fcs_front_xq_aI_slint_t finalize_mode;
  } mseg;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef fcs_front_xq_aI_MSS_ROOT
int root;
#endif
  } mss;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
fcs_front_xq_aI_slint_t sync;
  } msm;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
fcs_front_xq_aI_slint_t sync;
fcs_front_xq_aI_partcond_t *r_pc;
  } msp;
#endif
#ifdef SL_USE_MPI
  struct {
double i_t[3];
double p_t[3];
double b_t[3];
fcs_front_xq_aI_slint_t sync;
fcs_front_xq_aI_slint_t i_sync;
fcs_front_xq_aI_slint_t p_sync;
fcs_front_xq_aI_slint_t b_sync;
fcs_front_xq_aI_slint_t back_packed;
  } mssp;
#endif
} fcs_front_xq_aI_sl_context_t;






/* sl_macro fcs_front_xq_aI_elem_set_size fcs_front_xq_aI_elem_set_max_size fcs_front_xq_aI_elem_set_keys fcs_front_xq_aI_elem_set_indices */
#define fcs_front_xq_aI_elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define fcs_front_xq_aI_elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define fcs_front_xq_aI_elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define fcs_front_xq_aI_elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))
/* DATAX_TEMPLATE_BEGIN */
#define fcs_front_xq_aI_elem_set_data0(_e_, _b_)     ((_e_)->data0 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data1(_e_, _b_)     ((_e_)->data1 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data2(_e_, _b_)     ((_e_)->data2 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data3(_e_, _b_)     ((_e_)->data3 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data4(_e_, _b_)     ((_e_)->data4 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data5(_e_, _b_)     ((_e_)->data5 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data6(_e_, _b_)     ((_e_)->data6 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data7(_e_, _b_)     ((_e_)->data7 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data8(_e_, _b_)     ((_e_)->data8 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data9(_e_, _b_)     ((_e_)->data9 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data10(_e_, _b_)     ((_e_)->data10 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data11(_e_, _b_)     ((_e_)->data11 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data12(_e_, _b_)     ((_e_)->data12 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data13(_e_, _b_)     ((_e_)->data13 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data14(_e_, _b_)     ((_e_)->data14 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data15(_e_, _b_)     ((_e_)->data15 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data16(_e_, _b_)     ((_e_)->data16 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data17(_e_, _b_)     ((_e_)->data17 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data18(_e_, _b_)     ((_e_)->data18 = (_b_))  /* sl_macro */
#define fcs_front_xq_aI_elem_set_data19(_e_, _b_)     ((_e_)->data19 = (_b_))  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro fcs_front_xq_aI_elem_get_size fcs_front_xq_aI_elem_get_max_size fcs_front_xq_aI_elem_get_keys fcs_front_xq_aI_elem_get_indices */
#define fcs_front_xq_aI_elem_get_size(_e_)           (_e_)->size
#define fcs_front_xq_aI_elem_get_max_size(_e_)       (_e_)->max_size
#define fcs_front_xq_aI_elem_get_keys(_e_)           (_e_)->keys
#define fcs_front_xq_aI_elem_get_indices(_e_)        (_e_)->indices
/* DATAX_TEMPLATE_BEGIN */
#define fcs_front_xq_aI_elem_get_data0(_e_)          (_e_)->data0  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data1(_e_)          (_e_)->data1  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data2(_e_)          (_e_)->data2  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data3(_e_)          (_e_)->data3  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data4(_e_)          (_e_)->data4  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data5(_e_)          (_e_)->data5  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data6(_e_)          (_e_)->data6  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data7(_e_)          (_e_)->data7  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data8(_e_)          (_e_)->data8  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data9(_e_)          (_e_)->data9  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data10(_e_)          (_e_)->data10  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data11(_e_)          (_e_)->data11  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data12(_e_)          (_e_)->data12  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data13(_e_)          (_e_)->data13  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data14(_e_)          (_e_)->data14  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data15(_e_)          (_e_)->data15  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data16(_e_)          (_e_)->data16  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data17(_e_)          (_e_)->data17  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data18(_e_)          (_e_)->data18  /* sl_macro */
#define fcs_front_xq_aI_elem_get_data19(_e_)          (_e_)->data19  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro fcs_front_xq_aI_elem_set_block fcs_front_xq_aI_elem_set_block_size fcs_front_xq_aI_elem_get_block fcs_front_xq_aI_elem_get_block_size */
#define fcs_front_xq_aI_elem_set_block(_e_, _b_)       ((_e_)->keys = (_b_), (_e_)->max_size = -1)
#define fcs_front_xq_aI_elem_set_block_size(_e_, _s_)  ((_e_)->size = (_s_))
#define fcs_front_xq_aI_elem_get_block(_e_)            ((void *) (((_e_)->max_size < 0)?(_e_)->keys:NULL))
#define fcs_front_xq_aI_elem_get_block_size(_e_)       (_e_)->size

/* sl_macro fcs_front_xq_aI_pelem_set_size fcs_front_xq_aI_pelem_set_max_size fcs_front_xq_aI_pelem_set_elements */
#define fcs_front_xq_aI_pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define fcs_front_xq_aI_pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define fcs_front_xq_aI_pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))

/* sl_macro fcs_front_xq_aI_pelem_get_size fcs_front_xq_aI_pelem_get_max_size fcs_front_xq_aI_pelem_get_elements */
#define fcs_front_xq_aI_pelem_get_size(_e_)           (_e_)->size
#define fcs_front_xq_aI_pelem_get_max_size(_e_)       (_e_)->max_size
#define fcs_front_xq_aI_pelem_get_elements(_e_)       (_e_)->elements

/* sl_macro fcs_front_xq_aI_SL_DEFCON */
#define fcs_front_xq_aI_SL_DEFCON(_v_)  (fcs_front_xq_aI_sl_default_context._v_)






/* base/base.c */
extern fcs_front_xq_aI_sl_context_t fcs_front_xq_aI_sl_default_context;
extern const int fcs_front_xq_aI_default_sl_dummy_rank;
#ifdef fcs_front_xq_aI_SL_USE_RTI
extern const fcs_front_xq_aI_rti_t fcs_front_xq_aI_default_rti;
#endif
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_sr_ip_threshold;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_sr_db_threshold;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_sr_ma_threshold;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_sri_threshold;

/* base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
extern const MPI_Datatype fcs_front_xq_aI_default_mpi_int_datatype;
extern const MPI_Datatype fcs_front_xq_aI_default_mpi_key_datatype;
extern const MPI_Datatype fcs_front_xq_aI_default_mpi_pkey_datatype;
extern const MPI_Datatype fcs_front_xq_aI_default_mpi_pwkey_datatype;
extern const MPI_Datatype fcs_front_xq_aI_default_mpi_index_datatype;
extern const MPI_Datatype fcs_front_xq_aI_default_mpi_weight_datatype;
extern const MPI_Datatype fcs_front_xq_aI_default_mpi_data_datatype[];
extern const int fcs_front_xq_aI_default_mpi_rank;
#endif
extern const void *fcs_front_xq_aI_default_me_sendrecv_replace_mem;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_me_sendrecv_replace_memsize;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_me_sendrecv_replace_mpi_maxsize;
extern const double fcs_front_xq_aI_default_meas_t[];
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_meas_packed;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_meas_minalloc;
extern const double fcs_front_xq_aI_default_meas_overalloc;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_meas_type;
extern const void *fcs_front_xq_aI_default_meas_sendrecv_aux;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_meas_sendrecv_aux_size;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_meas_sendrecv_requests;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mea_packed;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mea_db_packed;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mea_ip_packed;
#ifdef fcs_front_xq_aI_MSEG_ROOT
extern const int fcs_front_xq_aI_default_mseg_root;
#endif
#ifdef fcs_front_xq_aI_MSEG_BORDER_UPDATE_REDUCTION
extern const double fcs_front_xq_aI_default_mseg_border_update_count_reduction;
extern const double fcs_front_xq_aI_default_mseg_border_update_weight_reduction;
#endif
#ifdef fcs_front_xq_aI_MSEG_FORWARD_ONLY
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mseg_forward_only;
#endif
#ifdef fcs_front_xq_aI_MSEG_INFO
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mseg_info_rounds;
extern const fcs_front_xq_aI_slint_t *fcs_front_xq_aI_default_mseg_info_finish_rounds;
extern const double fcs_front_xq_aI_default_mseg_info_finish_rounds_avg;
extern const fcs_front_xq_aI_slweight_t fcs_front_xq_aI_default_mseg_info_total_weights;
#endif
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mseg_binnings;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mseg_finalize_mode;
#ifdef fcs_front_xq_aI_MSS_ROOT
extern const int fcs_front_xq_aI_default_mss_root;
#endif
extern const double fcs_front_xq_aI_default_msm_t[];
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_msm_sync;
extern const double fcs_front_xq_aI_default_msp_t[];
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_msp_sync;
extern const fcs_front_xq_aI_partcond_t *fcs_front_xq_aI_default_msp_r_pc;
extern const double fcs_front_xq_aI_default_mssp_i_t[];
extern const double fcs_front_xq_aI_default_mssp_p_t[];
extern const double fcs_front_xq_aI_default_mssp_b_t[];
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mssp_sync;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mssp_i_sync;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mssp_p_sync;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mssp_b_sync;
extern const fcs_front_xq_aI_slint_t fcs_front_xq_aI_default_mssp_back_packed;






/* base/base.c */
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_create)(fcs_front_xq_aI_local_bins_t *lb, fcs_front_xq_aI_slint_t max_nbins, fcs_front_xq_aI_slint_t max_nbinnings, fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slint_t docounts, fcs_front_xq_aI_slint_t doweights, fcs_front_xq_aI_binning_t *bm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_destroy)(fcs_front_xq_aI_local_bins_t *lb);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_pre)(fcs_front_xq_aI_local_bins_t *lb);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_exec_reset)(fcs_front_xq_aI_local_bins_t *lb, fcs_front_xq_aI_slint_t do_bins, fcs_front_xq_aI_slint_t do_prefixes);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_exec)(fcs_front_xq_aI_local_bins_t *lb, fcs_front_xq_aI_slint_t b, fcs_front_xq_aI_slint_t do_bins, fcs_front_xq_aI_slint_t do_prefixes);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_refine)(fcs_front_xq_aI_local_bins_t *lb, fcs_front_xq_aI_slint_t b, fcs_front_xq_aI_slint_t k, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_hit)(fcs_front_xq_aI_local_bins_t *lb, fcs_front_xq_aI_slint_t b, fcs_front_xq_aI_slint_t k, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_finalize)(fcs_front_xq_aI_local_bins_t *lb, fcs_front_xq_aI_slint_t b, fcs_front_xq_aI_slint_t dc, fcs_front_xq_aI_slweight_t dw, fcs_front_xq_aI_slint_t lc_min, fcs_front_xq_aI_slint_t lc_max, fcs_front_xq_aI_slcount_t *lcs, fcs_front_xq_aI_slweight_t *lws, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_post)(fcs_front_xq_aI_local_bins_t *lb);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_radix_create)(fcs_front_xq_aI_binning_t *bm, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, fcs_front_xq_aI_slint_t sorted);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_radix_destroy)(fcs_front_xq_aI_binning_t *bm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_radix_pre)(fcs_front_xq_aI_binning_t *bm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_radix_exec)(fcs_front_xq_aI_binning_t *bm, fcs_front_xq_aI_bin_t *bin, fcs_front_xq_aI_slcount_t *counts, fcs_front_xq_aI_slweight_t *weights);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_radix_refine)(fcs_front_xq_aI_binning_t *bm, fcs_front_xq_aI_bin_t *bin, fcs_front_xq_aI_slint_t k, fcs_front_xq_aI_slcount_t *counts, fcs_front_xq_aI_slweight_t *weights, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s, fcs_front_xq_aI_bin_t *new_bin);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_radix_hit)(fcs_front_xq_aI_binning_t *bm, fcs_front_xq_aI_bin_t *bin, fcs_front_xq_aI_slint_t k, fcs_front_xq_aI_slcount_t *counts, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_radix_finalize)(fcs_front_xq_aI_binning_t *bm, fcs_front_xq_aI_bin_t *bin, fcs_front_xq_aI_slint_t dc, fcs_front_xq_aI_slweight_t dw, fcs_front_xq_aI_slint_t lc_min, fcs_front_xq_aI_slint_t lc_max, fcs_front_xq_aI_slcount_t *lcs, fcs_front_xq_aI_slweight_t *lws, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_binning_radix_post)(fcs_front_xq_aI_binning_t *bm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_alloc)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, slcint_t components);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_free)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_realloc)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, slcint_t components);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_alloca)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, slcint_t components);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_freea)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_alloc_from_blocks)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nblocks, void **blocks, fcs_front_xq_aI_slint_t *blocksizes, fcs_front_xq_aI_slint_t alignment, fcs_front_xq_aI_slint_t nmax, slcint_t components);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_alloc_from_block)(fcs_front_xq_aI_elements_t *s, void *block, fcs_front_xq_aI_slint_t blocksize, fcs_front_xq_aI_slint_t alignment, fcs_front_xq_aI_slint_t nmax, slcint_t components);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_block_alloc)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, slcint_t components);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_block_free)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_alloc_block)(fcs_front_xq_aI_elements_t *s, void **block, fcs_front_xq_aI_slint_t *blocksize, fcs_front_xq_aI_slint_t alignment, fcs_front_xq_aI_slint_t maxblocksize);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_copy)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_copy_at)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t sat, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t dat);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_ncopy)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t n);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_nmove)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t n);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_printf)(fcs_front_xq_aI_elements_t *s, const char *prefix);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_extract)(fcs_front_xq_aI_elements_t *src, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_elements_t *dst0, fcs_front_xq_aI_elements_t *dst1);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_touch)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_digest_sum)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(fcs_front_xq_aI_elements_crc32)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint nelements, fcs_front_xq_aI_slint_t keys, fcs_front_xq_aI_slint_t data);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_digest_hash)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, slcint_t components, void *hash);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_random_exchange)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t rounds, fcs_front_xq_aI_elements_t *xs);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_keys_init_seed)(unsigned long s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_keys_init)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_keys_init_type_t t, fcs_front_xq_aI_keys_init_data_t d);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_keys_init_randomized)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nkeys, fcs_front_xq_aI_keys_init_type_t t, fcs_front_xq_aI_keys_init_data_t d);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_keys_init_from_file)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t data, char *filename, fcs_front_xq_aI_slint_t from, fcs_front_xq_aI_slint_t to, fcs_front_xq_aI_slint_t const_bytes_per_line);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_keys_save_to_file)(fcs_front_xq_aI_elements_t *s, char *filename);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_validate_order)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t n);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_validate_order_bmask)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t n, fcs_front_xq_aI_slkey_pure_t bmask);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_validate_order_weight)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t n, fcs_front_xq_aI_slkey_pure_t weight);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_keys_stats)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slkey_pure_t *stats);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_keys_stats_print)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_print_keys)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_print_all)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slweight_t SL_PROTO(fcs_front_xq_aI_elements_get_weight)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_get_minmax_keys)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slkey_pure_t *minmaxkeys);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_alloc_packed)(fcs_front_xq_aI_packed_elements_t *s, fcs_front_xq_aI_slint_t nelements);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_free_packed)(fcs_front_xq_aI_packed_elements_t *s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_alloc_packed_from_block)(fcs_front_xq_aI_packed_elements_t *s, void *block, fcs_front_xq_aI_slint_t blocksize, fcs_front_xq_aI_slint_t alignment, fcs_front_xq_aI_slint_t nmax);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_pack_indexed)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_packed_elements_t *d, fcs_front_xq_aI_slindex_t *rindx, fcs_front_xq_aI_slindex_t *windx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_pack)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_packed_elements_t *d);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_pack_at)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t sat, fcs_front_xq_aI_packed_elements_t *d, fcs_front_xq_aI_slint_t dat);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_unpack_indexed)(fcs_front_xq_aI_packed_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slindex_t *rindx, fcs_front_xq_aI_slindex_t *windx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_unpack)(fcs_front_xq_aI_packed_elements_t *s, fcs_front_xq_aI_elements_t *d);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_unpack_at)(fcs_front_xq_aI_packed_elements_t *s, fcs_front_xq_aI_slint_t sat, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t dat);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elements_unpack_keys)(fcs_front_xq_aI_packed_elements_t *s, fcs_front_xq_aI_slkey_t *k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_merge2_basic_auto_01_x)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_merge2_basic_01_x)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_m2x_func _x0_1, fcs_front_xq_aI_m2x_func _0x_1);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_merge2_basic_01_X)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_elements_t *t, fcs_front_xq_aI_m2X_func _X0_1, fcs_front_xq_aI_m2X_func _0X_1);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_merge2_simplify_s1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint s1elements);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_merge2_memory_adaptive)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_compo_hula)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *xs);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_sseq_x0_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_sseq_0x_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_sseq_01_x)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_sseq_01)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_sbin_x0_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_sbin_0x_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_sbin_01_x)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_sbin_01)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_shyb_x0_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_shyb_0x_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_shyb_01_x)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_shyb_01)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_straight_x0_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_straight_0x_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_straight_01_x)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_straight_x_0_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_straight_X0_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_straight_0X_1)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_straight_01_X)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_basic_straight_X0_1u)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge2_compo_tridgell)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mergep_2way_ip_int)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t p, int *displs, fcs_front_xq_aI_merge2x_f m2x);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mergep_2way_ip_int_rec)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t p, int *displs, fcs_front_xq_aI_merge2x_f m2x);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mergep_heap_int)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t p, int *displs, int *counts);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mergep_heap_int_idx)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t p, int *displs, int *counts);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mergep_heap_idx)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t p, fcs_front_xq_aI_slindex_t *displs, fcs_front_xq_aI_slindex_t *counts);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mergep_heap_unpack_idx)(fcs_front_xq_aI_packed_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t p, fcs_front_xq_aI_slindex_t *displs, fcs_front_xq_aI_slindex_t *counts);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mergep_heap_unpack_idxonly)(fcs_front_xq_aI_packed_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t p, fcs_front_xq_aI_slindex_t *displs, fcs_front_xq_aI_slindex_t *counts);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_permute_generic_db)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_permute_generic_t *pg, void *pg_data);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_permute_generic_ip)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *x, fcs_front_xq_aI_permute_generic_t *pg, void *pg_data);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_sequential_lt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_sequential_le)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_sequential_gt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_sequential_ge)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_sequential_lt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_sequential_le)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_sequential_gt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_sequential_ge)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_binary_lt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_binary_le)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_binary_gt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_binary_ge)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_binary_lt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_binary_le)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_binary_gt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_binary_ge)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sl_search_binary_lt_bmask)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k, fcs_front_xq_aI_slpkey_t bmask);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sl_search_binary_le_bmask)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k, fcs_front_xq_aI_slpkey_t bmask);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sl_search_binary_sign_switch)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_hybrid_lt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k, fcs_front_xq_aI_slint t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_hybrid_le)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k, fcs_front_xq_aI_slint t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_hybrid_gt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k, fcs_front_xq_aI_slint t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_hybrid_ge)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t k, fcs_front_xq_aI_slint t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_hybrid_lt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k, fcs_front_xq_aI_slint t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_hybrid_le)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k, fcs_front_xq_aI_slint t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_hybrid_gt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k, fcs_front_xq_aI_slint t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sl_search_p_hybrid_ge)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slpkey_t *k, fcs_front_xq_aI_slint t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_ilog2c)(fcs_front_xq_aI_slint x);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_ilog2f)(fcs_front_xq_aI_slint x);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_print_bits)(fcs_front_xq_aI_slint v);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_pivot_random)(fcs_front_xq_aI_elements_t *s);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_counts2displs)(fcs_front_xq_aI_slint_t n, int *counts, int *displs);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_displs2counts)(fcs_front_xq_aI_slint_t n, int *displs, int *counts, fcs_front_xq_aI_slint_t total_counts);
void SL_PROTO(fcs_front_xq_aI_get_displcounts_extent)(fcs_front_xq_aI_slint_t n, int *displs, int *counts, fcs_front_xq_aI_slint_t *lb, fcs_front_xq_aI_slint_t *extent);
void SL_PROTO(fcs_front_xq_aI_elem_set_data)(fcs_front_xq_aI_elements_t *e, ...);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elem_get_max_byte)();
void SL_PROTO(fcs_front_xq_aI_elem_set)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d);
void SL_PROTO(fcs_front_xq_aI_elem_set_at)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t sat, fcs_front_xq_aI_elements_t *d);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elem_reverse)(fcs_front_xq_aI_elements_t *e, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elem_nxchange_at)(fcs_front_xq_aI_elements_t *e0, fcs_front_xq_aI_slint_t at0, fcs_front_xq_aI_elements_t *e1, fcs_front_xq_aI_slint_t at1, fcs_front_xq_aI_slint_t n, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elem_nxchange)(fcs_front_xq_aI_elements_t *e0, fcs_front_xq_aI_elements_t *e1, fcs_front_xq_aI_slint_t n, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elem_nxchange_ro0)(fcs_front_xq_aI_elements_t *e0, fcs_front_xq_aI_elements_t *e1, fcs_front_xq_aI_slint_t n, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elem_rotate)(fcs_front_xq_aI_elements_t *e, fcs_front_xq_aI_slint_t m, fcs_front_xq_aI_slint_t n, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elem_rotate_ro0)(fcs_front_xq_aI_elements_t *e, fcs_front_xq_aI_slint_t m, fcs_front_xq_aI_slint_t n, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_elem_rotate_ro1)(fcs_front_xq_aI_elements_t *e, fcs_front_xq_aI_slint_t m, fcs_front_xq_aI_slint_t n, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_counting_use_displs)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t ndispls, fcs_front_xq_aI_slint_t *displs);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_counting_use_counts)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t ncounts, fcs_front_xq_aI_slint_t *counts);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_counting_get_counts)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t ncounts, fcs_front_xq_aI_slint_t *counts);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_counting)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_slint_t ncounts);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sort_heap)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *xs);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_insert_bmask_kernel)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slkey_pure_t bmask);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_insert_kernel)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_insert)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_permute_forward)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t *perm, fcs_front_xq_aI_slint_t offset, fcs_front_xq_aI_slint_t mask_bit);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_permute_backward)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t *perm, fcs_front_xq_aI_slint_t offset, fcs_front_xq_aI_slint_t mask_bit);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sort_quick)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *xs);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_radix_ip)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_radix_db)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_radix_ma)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_radix)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_radix_1bit_kernel)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sort_radix_1bit)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_sort_radix_iter)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t presorted, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sn_hypercube_lh)(fcs_front_xq_aI_slint size, fcs_front_xq_aI_slint rank, fcs_front_xq_aI_slint stage, void *snp, fcs_front_xq_aI_slint *up);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sn_hypercube_hl)(fcs_front_xq_aI_slint size, fcs_front_xq_aI_slint rank, fcs_front_xq_aI_slint stage, void *snp, fcs_front_xq_aI_slint *up);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sn_odd_even_trans)(fcs_front_xq_aI_slint size, fcs_front_xq_aI_slint rank, fcs_front_xq_aI_slint stage, void *snp, fcs_front_xq_aI_slint *up);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sn_odd)(fcs_front_xq_aI_slint size, fcs_front_xq_aI_slint rank, fcs_front_xq_aI_slint stage, void *snp, fcs_front_xq_aI_slint *up);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sn_even)(fcs_front_xq_aI_slint size, fcs_front_xq_aI_slint rank, fcs_front_xq_aI_slint stage, void *snp, fcs_front_xq_aI_slint *up);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sn_batcher)(fcs_front_xq_aI_slint size, fcs_front_xq_aI_slint rank, fcs_front_xq_aI_slint stage, void *snp, fcs_front_xq_aI_slint *up);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sn_bitonic)(fcs_front_xq_aI_slint size, fcs_front_xq_aI_slint rank, fcs_front_xq_aI_slint stage, void *snp, fcs_front_xq_aI_slint *up);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_sn_connected)(fcs_front_xq_aI_slint size, fcs_front_xq_aI_slint rank, fcs_front_xq_aI_slint stage, void *snp, fcs_front_xq_aI_slint *up);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_split_generic_db)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_split_generic_t *sg, void *sg_data, fcs_front_xq_aI_slint_t n);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_split_generic_ip)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_split_generic_t *sg, void *sg_data, fcs_front_xq_aI_slint_t n);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_split_generic_count_db)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_split_generic_t *sg, void *sg_data, int *counts, fcs_front_xq_aI_slint_t n);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_split_generic_count_ip)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_split_generic_t *sg, void *sg_data, int *counts, fcs_front_xq_aI_slint_t n);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_split_generic_rearrange_db)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_split_generic_t *sg, void *sg_data, int *counts, fcs_front_xq_aI_slint_t n);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_split_generic_rearrange_ip)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *d, fcs_front_xq_aI_split_generic_t *sg, void *sg_data, int *counts, int *displs, fcs_front_xq_aI_slint_t n);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_splitter_reset)(fcs_front_xq_aI_splitter_t *sp);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_splitx_radix)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t nclasses, fcs_front_xq_aI_slint_t shl, fcs_front_xq_aI_slint_t *counts);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_split2_lt_ge)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slkey_pure_t *k, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_split2_le_gt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slkey_pure_t *k, fcs_front_xq_aI_elements_t *t);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_split3_lt_eq_gt)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slkey_pure_t *k, fcs_front_xq_aI_elements_t *t, fcs_front_xq_aI_slint *nlt, fcs_front_xq_aI_slint *nle);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_split3_lt_eq_gt_old)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slkey_pure_t *k, fcs_front_xq_aI_elements_t *t, fcs_front_xq_aI_slint *nlt, fcs_front_xq_aI_slint *nle);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_split2_b)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slkey_pure_t bmask);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_splitk_k2c_af)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint k, fcs_front_xq_aI_slint *c, fcs_front_xq_aI_k2c_func k2c, void *k2c_data);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_splitk_k2c)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint k, fcs_front_xq_aI_slint *c, fcs_front_xq_aI_k2c_func k2c, void *k2c_data);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_splitk_k2c_count)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint k, fcs_front_xq_aI_slint *c, fcs_front_xq_aI_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI





/* base_mpi/base_mpi.c */
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_create)(fcs_front_xq_aI_global_bins_t *gb, fcs_front_xq_aI_slint_t max_nbins, fcs_front_xq_aI_slint_t max_nbinnings, fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slint_t docounts, fcs_front_xq_aI_slint_t doweights, fcs_front_xq_aI_binning_t *bm, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_destroy)(fcs_front_xq_aI_global_bins_t *gb, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_pre)(fcs_front_xq_aI_global_bins_t *gb, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_exec_reset)(fcs_front_xq_aI_global_bins_t *gb, fcs_front_xq_aI_slint_t do_bins, fcs_front_xq_aI_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_exec_local)(fcs_front_xq_aI_global_bins_t *gb, fcs_front_xq_aI_slint_t b, fcs_front_xq_aI_slint_t do_bins, fcs_front_xq_aI_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_exec_global)(fcs_front_xq_aI_global_bins_t *gb, fcs_front_xq_aI_slint_t do_bins, fcs_front_xq_aI_slint_t do_prefixes, fcs_front_xq_aI_slint_t root, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_refine)(fcs_front_xq_aI_global_bins_t *gb, fcs_front_xq_aI_slint_t b, fcs_front_xq_aI_slint_t k, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_hit)(fcs_front_xq_aI_global_bins_t *gb, fcs_front_xq_aI_slint_t b, fcs_front_xq_aI_slint_t k, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_finalize)(fcs_front_xq_aI_global_bins_t *gb, fcs_front_xq_aI_slint_t b, fcs_front_xq_aI_slint_t dc, fcs_front_xq_aI_slweight_t dw, fcs_front_xq_aI_slint_t lc_min, fcs_front_xq_aI_slint_t lc_max, fcs_front_xq_aI_slcount_t *lcs, fcs_front_xq_aI_slweight_t *lws, fcs_front_xq_aI_splitter_t *sp, fcs_front_xq_aI_slint_t s, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_binning_post)(fcs_front_xq_aI_global_bins_t *gb, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_datatypes_init)();
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_datatypes_release)();
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_get_grid_properties)(fcs_front_xq_aI_slint_t ndims, fcs_front_xq_aI_slint_t *dims, fcs_front_xq_aI_slint_t *pos, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_subgroups_create)(fcs_front_xq_aI_slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_subgroups_delete)(fcs_front_xq_aI_slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);
int SL_PROTO(fcs_front_xq_aI_sl_MPI_Allreduce)(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank);
int SL_PROTO(fcs_front_xq_aI_sl_MPI_Alltoall_int)(void *sendbuf, int sendcount, void *recvbuf, int recvcount, MPI_Comm comm, int size, int rank);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_keys_init_from_file)(fcs_front_xq_aI_elements_t *s, char *filename, fcs_front_xq_aI_slint from, fcs_front_xq_aI_slint to, fcs_front_xq_aI_slint const_bytes_per_line, fcs_front_xq_aI_slint root, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint SL_PROTO(fcs_front_xq_aI_mpi_elements_validate_order)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint n, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_linear_exchange_pure_keys)(fcs_front_xq_aI_slkey_pure_t *in, fcs_front_xq_aI_slkey_pure_t *out, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_check_order)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slint_t *orders, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_check_global_order)(fcs_front_xq_aI_slkey_pure_t local_min, fcs_front_xq_aI_slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_digest_sum)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(fcs_front_xq_aI_mpi_elements_crc32)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t n, fcs_front_xq_aI_slint_t keys, fcs_front_xq_aI_slint_t data, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_digest_hash)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_get_counts)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t *clocal, fcs_front_xq_aI_slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slweight_t SL_PROTO(fcs_front_xq_aI_mpi_elements_get_weights)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slweight_t *wlocal, fcs_front_xq_aI_slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_get_counts_and_weights)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slint_t *counts, fcs_front_xq_aI_slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_sendrecv)(fcs_front_xq_aI_elements_t *sb, int sendcount, int dest, int sendtag, fcs_front_xq_aI_elements_t *rb, int recvcount, int recvtag, int source, fcs_front_xq_aI_slint_t *received, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_sendrecv_replace)(fcs_front_xq_aI_elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_isend_components)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t at, int count, int dest, int tag, MPI_Request *reqs, slcint_t components, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_irecv_components)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t at, int count, int source, int tag, MPI_Request *reqs, slcint_t components, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_tproc_create_tproc)(fcs_front_xq_aI_tproc_t *tproc, fcs_front_xq_aI_tproc_f *tfn, fcs_front_xq_aI_tproc_reset_f *rfn, fcs_front_xq_aI_tproc_exdef exdef);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_tproc_create_tproc_mod)(fcs_front_xq_aI_tproc_t *tproc, fcs_front_xq_aI_tproc_mod_f *tfn, fcs_front_xq_aI_tproc_reset_f *rfn, fcs_front_xq_aI_tproc_exdef exdef);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_tproc_create_tprocs)(fcs_front_xq_aI_tproc_t *tproc, fcs_front_xq_aI_slint_t max_tprocs, fcs_front_xq_aI_tprocs_f *tfn, fcs_front_xq_aI_tproc_reset_f *rfn, fcs_front_xq_aI_tproc_exdef exdef);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_tproc_create_tprocs_mod)(fcs_front_xq_aI_tproc_t *tproc, fcs_front_xq_aI_slint_t max_tprocs, fcs_front_xq_aI_tprocs_mod_f *tfn, fcs_front_xq_aI_tproc_reset_f *rfn, fcs_front_xq_aI_tproc_exdef exdef);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_tproc_free)(fcs_front_xq_aI_tproc_t *tproc);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_tproc_set_proclists)(fcs_front_xq_aI_tproc_t tproc, fcs_front_xq_aI_slint_t nsend_procs, int *send_procs, fcs_front_xq_aI_slint_t nrecv_procs, int *recv_procs, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_tproc_verify)(fcs_front_xq_aI_tproc_t tproc, void *data, fcs_front_xq_aI_elements_t *s, int proc);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_alltoall_specific)(fcs_front_xq_aI_elements_t *sin, fcs_front_xq_aI_elements_t *sout, fcs_front_xq_aI_elements_t *xs, fcs_front_xq_aI_tproc_t tproc, void *data, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_alltoallv_db_packed)(fcs_front_xq_aI_elements_t *sbuf, int *scounts, int *sdispls, fcs_front_xq_aI_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_alltoallv_db)(fcs_front_xq_aI_elements_t *sbuf, int *scounts, int *sdispls, fcs_front_xq_aI_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_alltoallv_ip_packed)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_alltoallv_ip_double)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_alltoallv_ip_mpi)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_alltoallv_ip_dash)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_alltoallv_ip)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_alltoallv_proclists_db)(fcs_front_xq_aI_elements_t *sbuf, int *scounts, int *sdispls, int nsendprocs, int *sendprocs, fcs_front_xq_aI_elements_t *rbuf, int *rcounts, int *rdispls, int nrecvprocs, int *recvprocs, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, fcs_front_xq_aI_slint_t structured);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_find_exact_equal)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t other_rank, fcs_front_xq_aI_slint_t high_rank, fcs_front_xq_aI_slint_t *ex_start, fcs_front_xq_aI_slint_t *ex_size, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_find_exact)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t other_rank, fcs_front_xq_aI_slint_t high_rank, fcs_front_xq_aI_slint_t *dst_size, fcs_front_xq_aI_slint_t *ex_start, fcs_front_xq_aI_slint_t *ex_sizes, fcs_front_xq_aI_slint_t *nx_move, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_merge2)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t other_rank, fcs_front_xq_aI_slint_t high_rank, fcs_front_xq_aI_slint_t *dst_size, fcs_front_xq_aI_merge2x_f m2, fcs_front_xq_aI_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_mergek_equal)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_sortnet_f sn, fcs_front_xq_aI_sortnet_data_t snd, fcs_front_xq_aI_merge2x_f m2x, fcs_front_xq_aI_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_mergek_sorted)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_merge2x_f m2x, fcs_front_xq_aI_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_mergek_sorted2)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_sortnet_f sn, fcs_front_xq_aI_sortnet_data_t snd, fcs_front_xq_aI_merge2x_f m2x, fcs_front_xq_aI_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_mergek)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_sortnet_f sn, fcs_front_xq_aI_sortnet_data_t snd, fcs_front_xq_aI_merge2x_f m2x, fcs_front_xq_aI_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_mergek_equal2)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_sortnet_f sn, fcs_front_xq_aI_sortnet_data_t snd, fcs_front_xq_aI_merge2x_f m2x, fcs_front_xq_aI_elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_partition_exact_generic)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_partcond_t *pcond, fcs_front_xq_aI_binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_partition_exact_radix)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_partcond_t *pcond, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, fcs_front_xq_aI_slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_partition_exact_radix_ngroups)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_partcond_t *pcond, fcs_front_xq_aI_slint_t ngroups, MPI_Comm *group_comms, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_partition_exact_radix_2groups)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_partcond_t *pcond, MPI_Comm group_comm, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_partition_sample_regular)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_rebalance)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_slint_t stable, fcs_front_xq_aI_slint_t *dst_size, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_rebalance_alltoallv)(fcs_front_xq_aI_elements_t *sbuf, int *scounts, int *sdispls, fcs_front_xq_aI_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
void SL_PROTO(fcs_front_xq_aI_mpi_partcond_set_even)(fcs_front_xq_aI_partcond_t *pcond, fcs_front_xq_aI_slint_t pcm, fcs_front_xq_aI_slint_t ntotal, double nimba, double wtotal, double wimba, int size, int rank);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_init_partconds)(fcs_front_xq_aI_slint_t npconds, fcs_front_xq_aI_partcond_t *pconds, fcs_front_xq_aI_slint_t nparts, fcs_front_xq_aI_slint_t total_count, fcs_front_xq_aI_slweight_t total_weight);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_init_partconds_intern)(fcs_front_xq_aI_slint_t npconds, fcs_front_xq_aI_partcond_intern_t *pci, fcs_front_xq_aI_partcond_t *pc, fcs_front_xq_aI_slint_t nparts, fcs_front_xq_aI_slint_t total_count, fcs_front_xq_aI_slweight_t total_weight);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_merge_partconds)(fcs_front_xq_aI_partcond_t *pconds_in, fcs_front_xq_aI_slint_t npconds_in, fcs_front_xq_aI_partcond_t *pcond_out);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_gather_partconds_grouped)(fcs_front_xq_aI_partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, fcs_front_xq_aI_partcond_t *pconds_out, fcs_front_xq_aI_slint_t *npconds_out, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_gather_partconds)(fcs_front_xq_aI_partcond_t *pcond_in, fcs_front_xq_aI_partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_allgather_partconds)(fcs_front_xq_aI_partcond_t *pcond_in, fcs_front_xq_aI_partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_bcast_partconds)(fcs_front_xq_aI_slint_t npconds, fcs_front_xq_aI_partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_post_check_partconds)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slint_t nparts, fcs_front_xq_aI_partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_post_check_partconds_intern)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slint_t nparts, fcs_front_xq_aI_partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_select_stats)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_select_exact_generic_bulk)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slint_t nparts, fcs_front_xq_aI_partcond_t *pconds, fcs_front_xq_aI_binning_t *bm, fcs_front_xq_aI_splitter_t *sp, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_select_exact_generic_grouped)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, fcs_front_xq_aI_binning_t *bm, fcs_front_xq_aI_splitter_t *sp, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_select_exact_generic)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slint_t nparts, fcs_front_xq_aI_partcond_t *pconds, fcs_front_xq_aI_binning_t *bm, fcs_front_xq_aI_splitter_t *sp, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_select_exact_radix)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_slint_t nparts, fcs_front_xq_aI_partcond_t *pconds, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, fcs_front_xq_aI_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_select_exact_radix_grouped)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nelements, fcs_front_xq_aI_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, fcs_front_xq_aI_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_select_sample_regular)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_slint_t nparts, fcs_front_xq_aI_partcond_t *pconds, fcs_front_xq_aI_slint_t nsamples, fcs_front_xq_aI_splitter_t *sp, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_merge)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_merge2)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *xs, fcs_front_xq_aI_slint_t merge_type, fcs_front_xq_aI_slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_merge_radix)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *xs, fcs_front_xq_aI_slint_t merge_type, fcs_front_xq_aI_slint_t sort_type, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_partition)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *xs, fcs_front_xq_aI_slint_t part_type, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_partition_radix)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *xs, fcs_front_xq_aI_slint_t part_type, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_partition_exact_radix)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_partcond_t *pcond, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_partition_exact_radix_ngroups)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_partcond_t *pcond, fcs_front_xq_aI_slint_t ngroups, MPI_Comm *group_comms, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_partition_exact_radix_2groups)(fcs_front_xq_aI_elements_t *s, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_partcond_t *pcond, MPI_Comm group_comm, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_insert_radix)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *xs, fcs_front_xq_aI_slpkey_t *mmkeys, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_presorted_radix)(fcs_front_xq_aI_elements_t *s0, fcs_front_xq_aI_elements_t *s1, fcs_front_xq_aI_elements_t *xs, fcs_front_xq_aI_slint_t merge_type, fcs_front_xq_aI_slint_t rhigh, fcs_front_xq_aI_slint_t rlow, fcs_front_xq_aI_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_sort_back)(fcs_front_xq_aI_elements_t *sin, fcs_front_xq_aI_elements_t *sout, fcs_front_xq_aI_elements_t *sx, fcs_front_xq_aI_slpkey_t *lh, fcs_front_xq_aI_slint_t ntotal, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, fcs_front_xq_aI_slint_t ytotal, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_xcounts2ycounts_grouped)(int *xcounts, fcs_front_xq_aI_slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
fcs_front_xq_aI_slint_t SL_PROTO(fcs_front_xq_aI_mpi_subxdispls2ycounts)(fcs_front_xq_aI_slint_t nsubs, int *sub_xdispls, fcs_front_xq_aI_slint_t *sub_sources, fcs_front_xq_aI_slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#undef SL_PROTO
#endif /* __SL_FRONT_XQ_AI_H__ */
