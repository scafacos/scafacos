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


#ifndef __SL_NEAR_FP_H__
#define __SL_NEAR_FP_H__

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif /* SL_USE_MPI */

#define SL_PROTO(_f_)  _f_

#define fcs_near_fp_sl_key_type_c        long long
#define fcs_near_fp_sl_key_type_mpi      MPI_LONG_LONG
#define fcs_near_fp_sl_key_size_mpi      1
#define fcs_near_fp_sl_key_type_fmt      "lld"

#define fcs_near_fp_sl_key_integer

#define fcs_near_fp_SL_DATA0             /* positions */
#define fcs_near_fp_sl_data0_type_c      fcs_float
#define fcs_near_fp_sl_data0_size_c      3
#define fcs_near_fp_sl_data0_type_mpi    FCS_MPI_FLOAT
#define fcs_near_fp_sl_data0_size_mpi    3

#define fcs_near_fp_SL_DATA1             /* charges */
#define fcs_near_fp_sl_data1_type_c      fcs_float
#define fcs_near_fp_sl_data1_size_c      1
#define fcs_near_fp_sl_data1_type_mpi    FCS_MPI_FLOAT
#define fcs_near_fp_sl_data1_size_mpi    1

#define fcs_near_fp_SL_DATA2             /* indices */
#define fcs_near_fp_sl_data2_type_c      long long
#define fcs_near_fp_sl_data2_size_c      1
#define fcs_near_fp_sl_data2_type_mpi    MPI_LONG_LONG
#define fcs_near_fp_sl_data2_size_mpi    1

#define fcs_near_fp_SL_DATA3             /* field */
#define fcs_near_fp_sl_data3_type_c      fcs_float
#define fcs_near_fp_sl_data3_size_c      3
#define fcs_near_fp_sl_data3_type_mpi    FCS_MPI_FLOAT
#define fcs_near_fp_sl_data3_size_mpi    3

#define fcs_near_fp_SL_DATA4             /* potentials */
#define fcs_near_fp_sl_data4_type_c      fcs_float
#define fcs_near_fp_sl_data4_size_c      1
#define fcs_near_fp_sl_data4_type_mpi    FCS_MPI_FLOAT
#define fcs_near_fp_sl_data4_size_mpi    1




#if defined(MSEG_ROOT) && !defined(fcs_near_fp_MSEG_ROOT)
# define fcs_near_fp_MSEG_ROOT  MSEG_ROOT
#endif

#if defined(MSEG_BORDER_UPDATE_REDUCTION) && !defined(fcs_near_fp_MSEG_BORDER_UPDATE_REDUCTION)
# define fcs_near_fp_MSEG_BORDER_UPDATE_REDUCTION  MSEG_BORDER_UPDATE_REDUCTION
#endif

#if defined(MSEG_DISABLE_BEST_CHOICE) && !defined(fcs_near_fp_MSEG_DISABLE_BEST_CHOICE)
# define fcs_near_fp_MSEG_DISABLE_BEST_CHOICE  MSEG_DISABLE_BEST_CHOICE
#endif

#if defined(MSEG_DISABLE_MINMAX) && !defined(fcs_near_fp_MSEG_DISABLE_MINMAX)
# define fcs_near_fp_MSEG_DISABLE_MINMAX  MSEG_DISABLE_MINMAX
#endif

#if defined(MSEG_ENABLE_OPTIMZED_LOWHIGH) && !defined(fcs_near_fp_MSEG_ENABLE_OPTIMZED_LOWHIGH)
# define fcs_near_fp_MSEG_ENABLE_OPTIMZED_LOWHIGH  MSEG_ENABLE_OPTIMZED_LOWHIGH
#endif

#if defined(MSEG_FORWARD_ONLY) && !defined(fcs_near_fp_MSEG_FORWARD_ONLY)
# define fcs_near_fp_MSEG_FORWARD_ONLY  MSEG_FORWARD_ONLY
#endif

#if defined(MSEG_INFO) && !defined(fcs_near_fp_MSEG_INFO)
# define fcs_near_fp_MSEG_INFO  MSEG_INFO
#endif

#if defined(MSEG_TRACE_IF) && !defined(fcs_near_fp_MSEG_TRACE_IF)
# define fcs_near_fp_MSEG_TRACE_IF  MSEG_TRACE_IF
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


#ifndef fcs_near_fp_SL_INDEX
# undef fcs_near_fp_SL_PACKED_INDEX
#endif


/* if no special datatype for (sl default) integer ... */
#ifndef fcs_near_fp_sl_int_type_c
  /* ... use a default one */
# define fcs_near_fp_sl_int_type_c               long      /* sl_macro */
# undef fcs_near_fp_sl_int_type_mpi
# define fcs_near_fp_sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef fcs_near_fp_sl_int_size_mpi
# define fcs_near_fp_sl_int_size_mpi             1         /* sl_macro */
# undef fcs_near_fp_sl_int_type_fmt
# define fcs_near_fp_sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(fcs_near_fp_sl_int_type_mpi) || !defined(fcs_near_fp_sl_int_size_mpi)
#   error "fcs_near_fp_sl_int_type_mpi and/or fcs_near_fp_sl_int_size_mpi missing"
#  endif
# endif
# ifndef fcs_near_fp_sl_int_type_fmt
#  error "fcs_near_fp_sl_int_type_fmt macro is missing, using d as default"
#  define fcs_near_fp_sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for (intern) weight ... */
#ifndef fcs_near_fp_sl_weight_type_c
 /* ... use (sl default) integer */
# define fcs_near_fp_sl_weight_type_c             fcs_near_fp_sl_int_type_c    /* sl_macro */
# undef fcs_near_fp_sl_weight_type_mpi
# define fcs_near_fp_sl_weight_type_mpi           fcs_near_fp_sl_int_type_mpi  /* sl_macro */
# undef fcs_near_fp_sl_weight_size_mpi
# define fcs_near_fp_sl_weight_size_mpi           fcs_near_fp_sl_int_size_mpi  /* sl_macro */
# undef fcs_near_fp_sl_weight_type_fmt
# define fcs_near_fp_sl_weight_type_fmt           fcs_near_fp_sl_int_type_fmt  /* sl_macro */
# undef fcs_near_fp_sl_weight_intequiv
# define fcs_near_fp_sl_weight_intequiv                            /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(fcs_near_fp_sl_weight_type_mpi) || !defined(fcs_near_fp_sl_weight_size_mpi)
#   error "fcs_near_fp_sl_weight_type_mpi and/or fcs_near_fp_sl_weight_size_mpi missing"
#  endif
# endif
# ifndef fcs_near_fp_sl_weight_type_fmt
#  error "fcs_near_fp_sl_weight_type_fmt macro is missing, using f as default"
#  define fcs_near_fp_sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef fcs_near_fp_sl_index_type_c
 /* ... use the primary integer type */
# define fcs_near_fp_sl_index_type_c             fcs_near_fp_sl_int_type_c
# undef fcs_near_fp_sl_index_type_mpi
# define fcs_near_fp_sl_index_type_mpi           fcs_near_fp_sl_int_type_mpi
# undef fcs_near_fp_sl_index_size_mpi
# define fcs_near_fp_sl_index_size_mpi           fcs_near_fp_sl_int_size_mpi
# undef fcs_near_fp_sl_index_type_fmt
# define fcs_near_fp_sl_index_type_fmt           fcs_near_fp_sl_int_type_fmt
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(fcs_near_fp_sl_index_type_mpi) || !defined(fcs_near_fp_sl_index_size_mpi)
#   error "fcs_near_fp_sl_index_type_mpi and/or fcs_near_fp_sl_index_size_mpi missing"
#  endif
# endif
# ifndef fcs_near_fp_sl_index_type_fmt
#  error "fcs_near_fp_sl_index_type_fmt macro is missing, using d as default"
#  define fcs_near_fp_sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef fcs_near_fp_sl_key_pure_type_c
# define fcs_near_fp_sl_key_pure_type_c          fcs_near_fp_sl_key_type_c  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_pure_type_mpi
# define fcs_near_fp_sl_key_pure_type_mpi        fcs_near_fp_sl_key_type_mpi  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_pure_size_mpi
# define fcs_near_fp_sl_key_pure_size_mpi        fcs_near_fp_sl_key_size_mpi  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_pure_type_fmt
# ifdef fcs_near_fp_sl_key_type_fmt
#  define fcs_near_fp_sl_key_pure_type_fmt       fcs_near_fp_sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef fcs_near_fp_sl_key_purify
 /* key val -> key val */
 #define fcs_near_fp_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define fcs_near_fp_sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_set_pure
 /* key component pointer and key val */
 #define fcs_near_fp_sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef fcs_near_fp_sl_key_pure_cmp_eq
 #define fcs_near_fp_sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_pure_cmp_ne
 #define fcs_near_fp_sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_pure_cmp_lt
 #define fcs_near_fp_sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_pure_cmp_le
 #define fcs_near_fp_sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_pure_cmp_gt
 #define fcs_near_fp_sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_pure_cmp_ge
 #define fcs_near_fp_sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef fcs_near_fp_sl_key_cmp_eq
 #define fcs_near_fp_sl_key_cmp_eq(k0, k1)       (fcs_near_fp_sl_key_pure_cmp_eq(fcs_near_fp_sl_key_purify(k0), fcs_near_fp_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_cmp_ne
 #define fcs_near_fp_sl_key_cmp_ne(k0, k1)       (fcs_near_fp_sl_key_pure_cmp_ne(fcs_near_fp_sl_key_purify(k0), fcs_near_fp_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_cmp_lt
 #define fcs_near_fp_sl_key_cmp_lt(k0, k1)       (fcs_near_fp_sl_key_pure_cmp_lt(fcs_near_fp_sl_key_purify(k0), fcs_near_fp_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_cmp_le
 #define fcs_near_fp_sl_key_cmp_le(k0, k1)       (fcs_near_fp_sl_key_pure_cmp_le(fcs_near_fp_sl_key_purify(k0), fcs_near_fp_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_cmp_gt
 #define fcs_near_fp_sl_key_cmp_gt(k0, k1)       (fcs_near_fp_sl_key_pure_cmp_gt(fcs_near_fp_sl_key_purify(k0), fcs_near_fp_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef fcs_near_fp_sl_key_cmp_ge
 #define fcs_near_fp_sl_key_cmp_ge(k0, k1)       (fcs_near_fp_sl_key_pure_cmp_ge(fcs_near_fp_sl_key_purify(k0), fcs_near_fp_sl_key_purify(k1)))  /* sl_macro */
#endif


/* default random key */
#ifdef fcs_near_fp_sl_key_integer
# if !defined(fcs_near_fp_sl_key_val_srand) || !defined(fcs_near_fp_sl_key_val_rand) || !defined(fcs_near_fp_sl_key_val_rand_minmax)
#  undef fcs_near_fp_sl_key_val_srand
#  undef fcs_near_fp_sl_key_val_rand
#  undef fcs_near_fp_sl_key_val_rand_minmax
#  define fcs_near_fp_sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define fcs_near_fp_sl_key_val_rand()                     ((fcs_near_fp_sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define fcs_near_fp_sl_key_val_rand_minmax(_min_, _max_)  ((fcs_near_fp_sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef fcs_near_fp_SL_DATA0_IGNORE
# undef fcs_near_fp_SL_DATA0
#endif
#ifdef fcs_near_fp_SL_DATA1_IGNORE
# undef fcs_near_fp_SL_DATA1
#endif
#ifdef fcs_near_fp_SL_DATA2_IGNORE
# undef fcs_near_fp_SL_DATA2
#endif
#ifdef fcs_near_fp_SL_DATA3_IGNORE
# undef fcs_near_fp_SL_DATA3
#endif
#ifdef fcs_near_fp_SL_DATA4_IGNORE
# undef fcs_near_fp_SL_DATA4
#endif
#ifdef fcs_near_fp_SL_DATA5_IGNORE
# undef fcs_near_fp_SL_DATA5
#endif
#ifdef fcs_near_fp_SL_DATA6_IGNORE
# undef fcs_near_fp_SL_DATA6
#endif
#ifdef fcs_near_fp_SL_DATA7_IGNORE
# undef fcs_near_fp_SL_DATA7
#endif
#ifdef fcs_near_fp_SL_DATA8_IGNORE
# undef fcs_near_fp_SL_DATA8
#endif
#ifdef fcs_near_fp_SL_DATA9_IGNORE
# undef fcs_near_fp_SL_DATA9
#endif
#ifdef fcs_near_fp_SL_DATA10_IGNORE
# undef fcs_near_fp_SL_DATA10
#endif
#ifdef fcs_near_fp_SL_DATA11_IGNORE
# undef fcs_near_fp_SL_DATA11
#endif
#ifdef fcs_near_fp_SL_DATA12_IGNORE
# undef fcs_near_fp_SL_DATA12
#endif
#ifdef fcs_near_fp_SL_DATA13_IGNORE
# undef fcs_near_fp_SL_DATA13
#endif
#ifdef fcs_near_fp_SL_DATA14_IGNORE
# undef fcs_near_fp_SL_DATA14
#endif
#ifdef fcs_near_fp_SL_DATA15_IGNORE
# undef fcs_near_fp_SL_DATA15
#endif
#ifdef fcs_near_fp_SL_DATA16_IGNORE
# undef fcs_near_fp_SL_DATA16
#endif
#ifdef fcs_near_fp_SL_DATA17_IGNORE
# undef fcs_near_fp_SL_DATA17
#endif
#ifdef fcs_near_fp_SL_DATA18_IGNORE
# undef fcs_near_fp_SL_DATA18
#endif
#ifdef fcs_near_fp_SL_DATA19_IGNORE
# undef fcs_near_fp_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro fcs_near_fp_sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef fcs_near_fp_sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef fcs_near_fp_sl_data0_weight
# undef fcs_near_fp_sl_data1_weight
# undef fcs_near_fp_sl_data2_weight
# undef fcs_near_fp_sl_data3_weight
# undef fcs_near_fp_sl_data4_weight
# undef fcs_near_fp_sl_data5_weight
# undef fcs_near_fp_sl_data6_weight
# undef fcs_near_fp_sl_data7_weight
# undef fcs_near_fp_sl_data8_weight
# undef fcs_near_fp_sl_data9_weight
# undef fcs_near_fp_sl_data10_weight
# undef fcs_near_fp_sl_data11_weight
# undef fcs_near_fp_sl_data12_weight
# undef fcs_near_fp_sl_data13_weight
# undef fcs_near_fp_sl_data14_weight
# undef fcs_near_fp_sl_data15_weight
# undef fcs_near_fp_sl_data16_weight
# undef fcs_near_fp_sl_data17_weight
# undef fcs_near_fp_sl_data18_weight
# undef fcs_near_fp_sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable fcs_near_fp_sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(fcs_near_fp_sl_data0_weight) && !defined(fcs_near_fp_SL_DATA0)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data1_weight) && !defined(fcs_near_fp_SL_DATA1)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data2_weight) && !defined(fcs_near_fp_SL_DATA2)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data3_weight) && !defined(fcs_near_fp_SL_DATA3)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data4_weight) && !defined(fcs_near_fp_SL_DATA4)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data5_weight) && !defined(fcs_near_fp_SL_DATA5)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data6_weight) && !defined(fcs_near_fp_SL_DATA6)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data7_weight) && !defined(fcs_near_fp_SL_DATA7)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data8_weight) && !defined(fcs_near_fp_SL_DATA8)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data9_weight) && !defined(fcs_near_fp_SL_DATA9)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data10_weight) && !defined(fcs_near_fp_SL_DATA10)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data11_weight) && !defined(fcs_near_fp_SL_DATA11)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data12_weight) && !defined(fcs_near_fp_SL_DATA12)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data13_weight) && !defined(fcs_near_fp_SL_DATA13)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data14_weight) && !defined(fcs_near_fp_SL_DATA14)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data15_weight) && !defined(fcs_near_fp_SL_DATA15)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data16_weight) && !defined(fcs_near_fp_SL_DATA16)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data17_weight) && !defined(fcs_near_fp_SL_DATA17)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data18_weight) && !defined(fcs_near_fp_SL_DATA18)
# undef fcs_near_fp_sl_elem_weight
#endif
#if defined(fcs_near_fp_sl_data19_weight) && !defined(fcs_near_fp_SL_DATA19)
# undef fcs_near_fp_sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro fcs_near_fp_FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef fcs_near_fp_SL_DATA0
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data0_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA1
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data1_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA2
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data2_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA3
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data3_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA4
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data4_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA5
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data5_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA6
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data6_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA7
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data7_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA8
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data8_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA9
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data9_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA10
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data10_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA11
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data11_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA12
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data12_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA13
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data13_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA14
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data14_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA15
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data15_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA16
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data16_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA17
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data17_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA18
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data18_flex
#   define fcs_near_fp_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef fcs_near_fp_SL_DATA19
# ifdef fcs_near_fp_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef fcs_near_fp_sl_data19_flex
#   define fcs_near_fp_FLECKS_GUARD
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






#define fcs_near_fp_SPEC_TLOC

typedef fcs_near_fp_sl_int_type_c fcs_near_fp_spec_int_t;

typedef int fcs_near_fp_spec_proc_t;

#define fcs_near_fp_SPEC_LOC_NONE   -1
#ifdef SL_USE_MPI
# define fcs_near_fp_SPEC_PROC_NONE  MPI_PROC_NULL
#else
# define fcs_near_fp_SPEC_PROC_NONE  -1
#endif

typedef void *fcs_near_fp_spec_tloc_data_t;
typedef void *fcs_near_fp_spec_tproc_data_t;

struct fcs_near_fp__elements_t;

typedef struct fcs_near_fp__elements_t *fcs_near_fp_spec_elem_buf_t;

typedef struct fcs_near_fp__elements_t fcs_near_fp_spec_elem_t;

typedef fcs_near_fp_sl_int_type_c fcs_near_fp_spec_elem_index_t;

#define fcs_near_fp_spec_elem_set_n(_e_, _n_)     fcs_near_fp_elem_set_size((_e_), (_n_))
#define fcs_near_fp_spec_elem_get_n(_e_)          fcs_near_fp_elem_get_size((_e_))
#define fcs_near_fp_spec_elem_set_nmax(_e_, _n_)  fcs_near_fp_elem_set_max_size((_e_), (_n_))
#define fcs_near_fp_spec_elem_get_nmax(_e_)       fcs_near_fp_elem_get_max_size((_e_))

#define fcs_near_fp_spec_elem_set_buf(_e_, _b_)   *(_e_) = *(_b_)
#define fcs_near_fp_spec_elem_get_buf(_e_)        (_e_)

#define fcs_near_fp_spec_elem_copy_at(_se_, _sat_, _de_, _dat_) \
  elem_copy_at((_se_), (_sat_), (_de_), (_dat_))

#define fcs_near_fp_spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) \
  elem_xchange_at((_s0_), (_s0at_), (_s1_), (_s1at_), (_t_))






/* tproc count */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROC_COUNT_DB */
#define fcs_near_fp_SPEC_DECLARE_TPROC_COUNT_DB \
  struct { fcs_near_fp_spec_elem_index_t i; fcs_near_fp_spec_proc_t p; } spec0cd;

/* sp_macro fcs_near_fp_SPEC_DO_TPROC_COUNT_DB */
#define fcs_near_fp_SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec0cd.i = 0; spec0cd.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec0cd.i) { \
    spec0cd.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == fcs_near_fp_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROC_COUNT_DB */
#define fcs_near_fp_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROC_COUNT_DB \
  fcs_near_fp_SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROC_COUNT_IP */
#define fcs_near_fp_SPEC_DECLARE_TPROC_COUNT_IP \
  struct { fcs_near_fp_spec_elem_index_t i, t; fcs_near_fp_spec_proc_t p; } spec0ci;

/* sp_macro fcs_near_fp_SPEC_DO_TPROC_COUNT_IP */
#define fcs_near_fp_SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec0ci.i) { \
    spec0ci.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == fcs_near_fp_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) fcs_near_fp_spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  fcs_near_fp_spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROC_COUNT_IP */
#define fcs_near_fp_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROC_COUNT_IP \
  fcs_near_fp_SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define fcs_near_fp_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { fcs_near_fp_spec_elem_index_t i; fcs_near_fp_spec_proc_t p; } spec1cd;

/* sp_macro fcs_near_fp_SPEC_DO_TPROC_MOD_COUNT_DB */
#define fcs_near_fp_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec1cd.i = 0; spec1cd.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec1cd.i) { \
    spec1cd.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == fcs_near_fp_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define fcs_near_fp_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  fcs_near_fp_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define fcs_near_fp_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { fcs_near_fp_spec_elem_index_t i, t; fcs_near_fp_spec_proc_t p; } spec1ci;

/* sp_macro fcs_near_fp_SPEC_DO_TPROC_MOD_COUNT_IP */
#define fcs_near_fp_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec1ci.i) { \
    spec1ci.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == fcs_near_fp_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) fcs_near_fp_spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  fcs_near_fp_spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define fcs_near_fp_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  fcs_near_fp_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROCS_COUNT_DB */
#define fcs_near_fp_SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { fcs_near_fp_spec_elem_index_t i; fcs_near_fp_spec_int_t j, n; } spec2cd;

/* sp_macro fcs_near_fp_SPEC_DO_TPROCS_COUNT_DB */
#define fcs_near_fp_SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec2cd.i = 0; spec2cd.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec2cd.i) { \
    spec2cd.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.n; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROCS_COUNT_DB */
#define fcs_near_fp_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts, fcs_near_fp_spec_proc_t *procs) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROCS_COUNT_DB \
  fcs_near_fp_SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROCS_COUNT_IP */
#define fcs_near_fp_SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { fcs_near_fp_spec_elem_index_t i, t; fcs_near_fp_spec_int_t j, n; } spec2ci;

/* sp_macro fcs_near_fp_SPEC_DO_TPROCS_COUNT_IP */
#define fcs_near_fp_SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec2ci.i) { \
    spec2ci.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), (_ps_)); \
    if (spec2ci.n <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.n; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) fcs_near_fp_spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  fcs_near_fp_spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROCS_COUNT_IP */
#define fcs_near_fp_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts, fcs_near_fp_spec_proc_t *procs) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROCS_COUNT_IP \
  fcs_near_fp_SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { fcs_near_fp_spec_elem_index_t i; fcs_near_fp_spec_int_t j, n; } spec3cd;

/* sp_macro fcs_near_fp_SPEC_DO_TPROCS_MOD_COUNT_DB */
#define fcs_near_fp_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec3cd.i = 0; spec3cd.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec3cd.i) \
  { \
    spec3cd.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.n; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define fcs_near_fp_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts, fcs_near_fp_spec_proc_t *procs) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  fcs_near_fp_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { fcs_near_fp_spec_elem_index_t i, t; fcs_near_fp_spec_int_t j, n; } spec3ci;

/* sp_macro fcs_near_fp_SPEC_DO_TPROCS_MOD_COUNT_IP */
#define fcs_near_fp_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec3ci.i) { \
    spec3ci.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), (_ps_), NULL); \
    if (spec3ci.n <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.n; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) fcs_near_fp_spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  fcs_near_fp_spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define fcs_near_fp_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts, fcs_near_fp_spec_proc_t *procs) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  fcs_near_fp_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tproc rearrange */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROC_REARRANGE_DB */
#define fcs_near_fp_SPEC_DECLARE_TPROC_REARRANGE_DB \
  struct { fcs_near_fp_spec_elem_index_t i; fcs_near_fp_spec_proc_t p; } spec0d;

/* sp_macro fcs_near_fp_SPEC_DO_TPROC_REARRANGE_DB */
#define fcs_near_fp_SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  for (spec0d.i = 0; spec0d.i < fcs_near_fp_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == fcs_near_fp_SPEC_PROC_NONE) continue; \
    fcs_near_fp_spec_elem_copy_at((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p]); \
    ++(_ds_)[spec0d.p]; \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROC_REARRANGE_DB */
#define fcs_near_fp_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROC_REARRANGE_DB \
  fcs_near_fp_SPEC_DO_TPROC_REARRANGE_DB(_tp_, tproc_data, s, d, displs); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROC_REARRANGE_IP */
#define fcs_near_fp_SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { fcs_near_fp_spec_elem_index_t e, i, j; fcs_near_fp_spec_proc_t p, np; } spec0i;

/* sp_macro fcs_near_fp_SPEC_DO_TPROC_REARRANGE_IP */
#define fcs_near_fp_SPEC_DO_TPROC_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_)  do { \
  for (spec0i.e = 0, spec0i.i = 0; spec0i.i < (_n_); ++spec0i.i) { \
    spec0i.e += (_cs_)[spec0i.i]; \
    spec0i.j = (_ds_)[spec0i.i]; \
    while (spec0i.j < spec0i.e) { \
      spec0i.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec0i.j, _tpd_); \
      while (spec0i.p != spec0i.i) { \
        spec0i.np = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), (_ds_)[spec0i.p], _tpd_); \
        if (spec0i.np != spec0i.p) fcs_near_fp_spec_elem_exchange_at((_b_), (_ds_)[spec0i.p], (_b_), spec0i.j, (_xb_)); \
        ++(_ds_)[spec0i.p]; \
        spec0i.p = spec0i.np; \
      } \
      ++spec0i.j; \
    } \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROC_REARRANGE_IP */
#define fcs_near_fp_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_rearrange_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, int *counts, fcs_near_fp_spec_int_t n) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROC_REARRANGE_IP \
  fcs_near_fp_SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define fcs_near_fp_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  struct { fcs_near_fp_spec_elem_index_t i; fcs_near_fp_spec_proc_t p; } spec1d;

/* sp_macro fcs_near_fp_SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define fcs_near_fp_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < fcs_near_fp_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, fcs_near_fp_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == fcs_near_fp_SPEC_PROC_NONE) continue; \
      fcs_near_fp_spec_elem_copy_at((_ib_), 0, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < fcs_near_fp_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, NULL); \
      if (spec1d.p == fcs_near_fp_SPEC_PROC_NONE) continue; \
      fcs_near_fp_spec_elem_copy_at((_sb_), spec1d.i, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define fcs_near_fp_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, fcs_near_fp_spec_elem_t *mod) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  fcs_near_fp_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define fcs_near_fp_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { fcs_near_fp_spec_elem_index_t e, i, j; fcs_near_fp_spec_proc_t p, np; } spec1i;

/* sp_macro fcs_near_fp_SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define fcs_near_fp_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
      spec1i.e += (_cs_)[spec1i.i]; \
      spec1i.j = (_ds_)[spec1i.i]; \
      while (spec1i.j < spec1i.e) { \
        spec1i.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec1i.j, _tpd_, fcs_near_fp_spec_elem_get_buf(_ib_)); \
        fcs_near_fp_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, fcs_near_fp_spec_elem_get_buf(_ib_)); \
          if (spec1i.np != spec1i.p) { \
            fcs_near_fp_spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
            fcs_near_fp_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
          } else fcs_near_fp_spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
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
        spec1i.p = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec1i.j, _tpd_, NULL); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, NULL); \
          if (spec1i.np != spec1i.p) fcs_near_fp_spec_elem_exchange_at((_b_), (_ds_)[spec1i.p], (_b_), spec1i.j, (_xb_)); \
          ++(_ds_)[spec1i.p]; \
          spec1i.p = spec1i.np; \
        } \
        ++spec1i.j; \
      } \
    } \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define fcs_near_fp_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, int *counts, fcs_near_fp_spec_int_t n, fcs_near_fp_spec_elem_t *mod) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  fcs_near_fp_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define fcs_near_fp_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  struct { fcs_near_fp_spec_elem_index_t i; fcs_near_fp_spec_int_t j, n; } spec2d;

/* sp_macro fcs_near_fp_SPEC_DO_TPROCS_REARRANGE_DB */
#define fcs_near_fp_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  for (spec2d.i = 0; spec2d.i < fcs_near_fp_spec_elem_get_n(_sb_); ++spec2d.i) { \
    spec2d.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.n; ++spec2d.j) { \
      fcs_near_fp_spec_elem_copy_at((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]]); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROCS_REARRANGE_DB */
#define fcs_near_fp_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, fcs_near_fp_spec_proc_t *procs) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  fcs_near_fp_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define fcs_near_fp_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { fcs_near_fp_spec_elem_index_t e, j, fe, fc, le, lc; fcs_near_fp_spec_int_t i, n, f, l, o; } spec2i;

/* sp_macro fcs_near_fp_SPEC_DO_TPROCS_REARRANGE_IP */
#define fcs_near_fp_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_)  do { \
  spec2i.f = 0; spec2i.fe = (_cs_)[0]; spec2i.fc = fcs_near_fp_spec_elem_get_n(_b_); \
  while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; } \
  spec2i.l = 0; spec2i.le = (_cs_)[0]; spec2i.lc = fcs_near_fp_spec_elem_get_n(_b_) - 1; \
  while (spec2i.lc >= spec2i.le) { ++spec2i.l; spec2i.le += (_cs_)[spec2i.l]; } \
  for (spec2i.e = 0, spec2i.i = 0; spec2i.i < (_n_); ++spec2i.i) { \
    spec2i.e += (_cs_)[spec2i.i]; \
    spec2i.j = (_ds_)[spec2i.i]; \
    while (spec2i.j < spec2i.e) { \
      spec2i.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec2i.j, (_tpd_), (_ps_)); \
      spec2i.o = -1; \
      while (spec2i.n > 0) { \
        --spec2i.n; \
        if ((_ps_)[spec2i.n] == spec2i.i && spec2i.o < 0) spec2i.o = spec2i.n; \
        else if ((_ds_)[(_ps_)[spec2i.n]] < spec2i.fc) { \
          spec2i.l = spec2i.f; spec2i.le = spec2i.fe; spec2i.lc = spec2i.fc; \
          if (spec2i.fc < spec2i.fe) { \
            fcs_near_fp_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_b_), spec2i.fc); \
            ++spec2i.fc; \
          } else fcs_near_fp_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec2i.n]] == spec2i.fc) ++spec2i.fc; \
        if (spec2i.j != (_ds_)[(_ps_)[spec2i.n]]) fcs_near_fp_spec_elem_copy_at((_b_), spec2i.j, (_b_), (_ds_)[(_ps_)[spec2i.n]]); \
        ++(_ds_)[(_ps_)[spec2i.n]]; \
        while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; spec2i.fc = (_ds_)[spec2i.f]; } \
      } \
      if (spec2i.o < 0) { \
        if (spec2i.lc < spec2i.le) {  \
          fcs_near_fp_spec_elem_copy_at((_b_), spec2i.lc, (_b_), spec2i.j); \
          spec2i.f = spec2i.l; spec2i.fe = spec2i.le; spec2i.fc = spec2i.lc; \
          --spec2i.lc; \
          while (spec2i.l > 0 && spec2i.lc < (_ds_)[spec2i.l]) { spec2i.le -= (_cs_)[spec2i.l]; spec2i.lc = spec2i.le - 1; --spec2i.l; } \
        } else fcs_near_fp_spec_elem_copy_at((_xb_), 0, (_b_), spec2i.j); \
      } \
      spec2i.j = (_ds_)[spec2i.i]; \
    } \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROCS_REARRANGE_IP */
#define fcs_near_fp_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, int *counts, fcs_near_fp_spec_int_t n, fcs_near_fp_spec_proc_t *procs) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  fcs_near_fp_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  struct { fcs_near_fp_spec_elem_index_t i; fcs_near_fp_spec_int_t j, n; } spec3d;

/* sp_macro fcs_near_fp_SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define fcs_near_fp_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  if (_ib_) { \
    for (spec3d.i = 0; spec3d.i < fcs_near_fp_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), fcs_near_fp_spec_elem_get_buf(_ib_)); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        fcs_near_fp_spec_elem_copy_at((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } else { \
    for (spec3d.i = 0; spec3d.i < fcs_near_fp_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), NULL); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        fcs_near_fp_spec_elem_copy_at((_sb_), spec3d.i, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define fcs_near_fp_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, fcs_near_fp_spec_proc_t *procs, fcs_near_fp_spec_elem_t *mod) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  fcs_near_fp_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { fcs_near_fp_spec_elem_index_t e, j, fe, fc, le, lc; fcs_near_fp_spec_int_t i, n, f, l, o; } spec3i;

/* sp_macro fcs_near_fp_SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define fcs_near_fp_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  if (_ib_) { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = fcs_near_fp_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = fcs_near_fp_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), fcs_near_fp_spec_elem_get_buf(_ib_)); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              fcs_near_fp_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else fcs_near_fp_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          fcs_near_fp_spec_elem_copy_at((_ib_), spec3i.n, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            fcs_near_fp_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else fcs_near_fp_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } else { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = fcs_near_fp_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = fcs_near_fp_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(fcs_near_fp_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), NULL); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              fcs_near_fp_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else fcs_near_fp_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          if (spec3i.j != (_ds_)[(_ps_)[spec3i.n]]) fcs_near_fp_spec_elem_copy_at((_b_), spec3i.j, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            fcs_near_fp_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else fcs_near_fp_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define fcs_near_fp_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, int *counts, fcs_near_fp_spec_int_t n, fcs_near_fp_spec_proc_t *procs, fcs_near_fp_spec_elem_t *mod) \
{ \
  fcs_near_fp_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  fcs_near_fp_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}

/* sp_macro fcs_near_fp_SPEC_DEFINE_TPROC */
#define fcs_near_fp_SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  fcs_near_fp_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro fcs_near_fp_SPEC_DEFINE_TPROC_MOD */
#define fcs_near_fp_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  fcs_near_fp_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro fcs_near_fp_SPEC_DEFINE_TPROCS */
#define fcs_near_fp_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  fcs_near_fp_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro fcs_near_fp_SPEC_DEFINE_TPROCS_MOD */
#define fcs_near_fp_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  fcs_near_fp_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  fcs_near_fp_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro fcs_near_fp_SPEC_EXT_PARAM_TPROC fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL fcs_near_fp_SPEC_EXT_PARAM_TPROCS fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define fcs_near_fp_SPEC_EXT_PARAM_TPROC(_name_)       _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip
#define fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL          NULL, NULL, NULL, NULL
#define fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD(_name_)   _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip
#define fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL      NULL, NULL, NULL, NULL
#define fcs_near_fp_SPEC_EXT_PARAM_TPROCS(_name_)      _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip
#define fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL         NULL, NULL, NULL, NULL
#define fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD(_name_)  _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip
#define fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL     NULL, NULL, NULL, NULL


/* sp_type fcs_near_fp_spec_tproc_f fcs_near_fp_spec_tproc_count_f fcs_near_fp_spec_tproc_rearrange_db_f fcs_near_fp_spec_tproc_rearrange_ip_f */
typedef fcs_near_fp_spec_proc_t fcs_near_fp_spec_tproc_f(fcs_near_fp_spec_elem_buf_t b, fcs_near_fp_spec_elem_index_t x, fcs_near_fp_spec_tproc_data_t tproc_data);
typedef void fcs_near_fp_spec_tproc_count_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts);
typedef void fcs_near_fp_spec_tproc_rearrange_db_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs);
typedef void fcs_near_fp_spec_tproc_rearrange_ip_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, int *counts, fcs_near_fp_spec_int_t n);

/* sp_type fcs_near_fp_spec_tproc_mod_f fcs_near_fp_spec_tproc_mod_count_f fcs_near_fp_spec_tproc_mod_rearrange_db_f fcs_near_fp_spec_tproc_mod_rearrange_ip_f */
typedef fcs_near_fp_spec_proc_t fcs_near_fp_spec_tproc_mod_f(fcs_near_fp_spec_elem_buf_t b, fcs_near_fp_spec_elem_index_t x, fcs_near_fp_spec_tproc_data_t tproc_data, fcs_near_fp_spec_elem_buf_t mod);
typedef void fcs_near_fp_spec_tproc_mod_count_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts);
typedef void fcs_near_fp_spec_tproc_mod_rearrange_db_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, fcs_near_fp_spec_elem_t *mod);
typedef void fcs_near_fp_spec_tproc_mod_rearrange_ip_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, int *counts, fcs_near_fp_spec_int_t n, fcs_near_fp_spec_elem_t *mod);

/* sp_type fcs_near_fp_spec_tprocs_f fcs_near_fp_spec_tprocs_count_f fcs_near_fp_spec_tprocs_rearrange_db_f fcs_near_fp_spec_tprocs_rearrange_ip_f */
typedef fcs_near_fp_spec_int_t fcs_near_fp_spec_tprocs_f(fcs_near_fp_spec_elem_buf_t b, fcs_near_fp_spec_elem_index_t x, fcs_near_fp_spec_tproc_data_t tproc_data, fcs_near_fp_spec_proc_t *procs);
typedef void fcs_near_fp_spec_tprocs_count_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts, fcs_near_fp_spec_proc_t *procs);
typedef void fcs_near_fp_spec_tprocs_rearrange_db_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, fcs_near_fp_spec_proc_t *procs);
typedef void fcs_near_fp_spec_tprocs_rearrange_ip_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, int *counts, fcs_near_fp_spec_int_t n, fcs_near_fp_spec_proc_t *procs);

/* sp_type fcs_near_fp_spec_tprocs_mod_f fcs_near_fp_spec_tprocs_mod_count_f fcs_near_fp_spec_tprocs_mod_rearrange_db_f fcs_near_fp_spec_tprocs_mod_rearrange_ip_f */
typedef fcs_near_fp_spec_int_t fcs_near_fp_spec_tprocs_mod_f(fcs_near_fp_spec_elem_buf_t b, fcs_near_fp_spec_elem_index_t x, fcs_near_fp_spec_tproc_data_t tproc_data, fcs_near_fp_spec_proc_t *procs, fcs_near_fp_spec_elem_buf_t mod);
typedef void fcs_near_fp_spec_tprocs_mod_count_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_tproc_data_t tproc_data, int *counts, fcs_near_fp_spec_proc_t *procs);
typedef void fcs_near_fp_spec_tprocs_mod_rearrange_db_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, fcs_near_fp_spec_proc_t *procs, fcs_near_fp_spec_elem_t *mod);
typedef void fcs_near_fp_spec_tprocs_mod_rearrange_ip_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tproc_data_t tproc_data, int *displs, int *counts, fcs_near_fp_spec_int_t n, fcs_near_fp_spec_proc_t *procs, fcs_near_fp_spec_elem_t *mod);

/* sp_type fcs_near_fp_spec_tproc_reset_f */
typedef void fcs_near_fp_spec_tproc_reset_f(fcs_near_fp_spec_tproc_data_t tproc_data);


/* enable tloc features */
#ifdef fcs_near_fp_SPEC_TLOC

/* sp_macro fcs_near_fp_SPEC_TLOC fcs_near_fp_SPEC_LOC_NONE */


/* tloc rearrange */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TLOC_REARRANGE_DB */
#define fcs_near_fp_SPEC_DECLARE_TLOC_REARRANGE_DB \
  struct { fcs_near_fp_spec_int_t i, p; } spec0d;

/* sp_macro fcs_near_fp_SPEC_DO_TLOC_REARRANGE_DB */
#define fcs_near_fp_SPEC_DO_TLOC_REARRANGE_DB(_tl_, _tld_, _sb_, _db_)  do { \
  for (spec0d.i = 0; spec0d.i < fcs_near_fp_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tl_)(fcs_near_fp_spec_elem_get_buf(_sb_), spec0d.i, _tld_); \
    if (spec0d.p == fcs_near_fp_SPEC_LOC_NONE) continue; \
    fcs_near_fp_spec_elem_copy_at((_sb_), spec0d.i, (_db_), spec0d.p); \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TLOC_REARRANGE_DB */
#define fcs_near_fp_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_rearrange_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tloc_data_t tloc_data) \
{ \
  fcs_near_fp_SPEC_DECLARE_TLOC_REARRANGE_DB \
  fcs_near_fp_SPEC_DO_TLOC_REARRANGE_DB(_tl_, tloc_data, s, d); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TLOC_REARRANGE_IP */
#define fcs_near_fp_SPEC_DECLARE_TLOC_REARRANGE_IP \
  struct { fcs_near_fp_spec_int_t i, p, np; } spec0i;

/* sp_macro fcs_near_fp_SPEC_DO_TLOC_REARRANGE_IP */
#define fcs_near_fp_SPEC_DO_TLOC_REARRANGE_IP(_tl_, _tld_, _b_, _xb_)  do { \
  for (spec0i.i = 0; spec0i.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec0i.i) { \
    spec0i.p = (_tl_)(fcs_near_fp_spec_elem_get_buf(_b_), spec0i.i, _tld_); \
    if (spec0i.p == fcs_near_fp_SPEC_LOC_NONE) continue; \
    while (spec0i.i != spec0i.p) { \
      spec0i.np = (_tl_)(fcs_near_fp_spec_elem_get_buf(_b_), spec0i.p, _tld_); \
      if (spec0i.np == fcs_near_fp_SPEC_LOC_NONE) { fcs_near_fp_spec_elem_copy_at((_b_), spec0i.i, (_b_), spec0i.p); break; } \
      fcs_near_fp_spec_elem_exchange_at((_b_), spec0i.i, (_b_), spec0i.p, (_xb_)); \
      spec0i.p = spec0i.np; \
    } \
  } } while (0)

/* sp_macro fcs_near_fp_SPEC_FUNC_TLOC_REARRANGE_IP */
#define fcs_near_fp_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_rearrange_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tloc_data_t tloc_data) \
{ \
  fcs_near_fp_SPEC_DECLARE_TLOC_REARRANGE_IP \
  fcs_near_fp_SPEC_DO_TLOC_REARRANGE_IP(_tl_, tloc_data, s, x); \
}


/* tloc_mod_mod rearrange */

/* sp_macro fcs_near_fp_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB */
#define fcs_near_fp_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  struct { fcs_near_fp_spec_int_t i, p; } spec1d;

/* sp_macro fcs_near_fp_SPEC_DO_TLOC_MOD_REARRANGE_DB */
#define fcs_near_fp_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, _tld_, _sb_, _db_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < fcs_near_fp_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(fcs_near_fp_spec_elem_get_buf(_sb_), spec1d.i, _tld_, fcs_near_fp_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == fcs_near_fp_SPEC_LOC_NONE) continue; \
      fcs_near_fp_spec_elem_copy_at((_ib_), 0, (_db_), spec1d.p); \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < fcs_near_fp_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(fcs_near_fp_spec_elem_get_buf(_sb_), spec1d.i, _tld_, NULL); \
      if (spec1d.p == fcs_near_fp_SPEC_LOC_NONE) continue; \
      fcs_near_fp_spec_elem_copy_at((_sb_), spec1d.i, (_db_), spec1d.p); \
    } \
  } } while (0) 

/* sp_macro fcs_near_fp_SPEC_FUNC_TLOC_MOD_REARRANGE_DB */
#define fcs_near_fp_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_mod_rearrange_db(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tloc_data_t tloc_data, fcs_near_fp_spec_elem_t *mod) \
{ \
  fcs_near_fp_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  fcs_near_fp_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, tloc_data, s, d, mod); \
}

/* sp_macro fcs_near_fp_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP */
#define fcs_near_fp_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  struct { fcs_near_fp_spec_int_t i, p, np; } spec1i;

/* sp_macro fcs_near_fp_SPEC_DO_TLOC_MOD_REARRANGE_IP */
#define fcs_near_fp_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, _tld_, _b_, _xb_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.i = 0; spec1i.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(fcs_near_fp_spec_elem_get_buf(_b_), spec1i.i, _tld_, fcs_near_fp_spec_elem_get_buf(_ib_)); \
      if (spec1i.p == fcs_near_fp_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(fcs_near_fp_spec_elem_get_buf(_b_), spec1i.p, _tld_, fcs_near_fp_spec_elem_get_buf(_xb_)); \
        if (spec1i.np == fcs_near_fp_SPEC_LOC_NONE) break; \
        fcs_near_fp_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.p); \
        fcs_near_fp_spec_elem_copy_at((_xb_), 0, (_ib_), 0); \
        spec1i.p = spec1i.np; \
      } \
      fcs_near_fp_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.i); \
    } \
  } else { \
    for (spec1i.i = 0; spec1i.i < fcs_near_fp_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(fcs_near_fp_spec_elem_get_buf(_b_), spec1i.i, _tld_, NULL); \
      if (spec1i.p == fcs_near_fp_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(fcs_near_fp_spec_elem_get_buf(_b_), spec1i.p, _tld_, NULL); \
        if (spec1i.np == fcs_near_fp_SPEC_LOC_NONE) { fcs_near_fp_spec_elem_copy_at((_b_), spec1i.i, (_b_), spec1i.p); break; } \
        fcs_near_fp_spec_elem_exchange_at((_b_), spec1i.i, (_b_), spec1i.p, (_xb_)); \
        spec1i.p = spec1i.np; \
      } \
    } \
 } } while (0) 

/* sp_macro fcs_near_fp_SPEC_FUNC_TLOC_MOD_REARRANGE_IP */
#define fcs_near_fp_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_mod_rearrange_ip(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tloc_data_t tloc_data, fcs_near_fp_spec_elem_t *mod) \
{ \
  fcs_near_fp_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  fcs_near_fp_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, tloc_data, s, x, mod); \
}

/* sp_macro fcs_near_fp_SPEC_DEFINE_TLOC */
#define fcs_near_fp_SPEC_DEFINE_TLOC(_name_, _tl_, _s_...) \
  fcs_near_fp_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_) \
  fcs_near_fp_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro fcs_near_fp_SPEC_DEFINE_TLOC_MOD */
#define fcs_near_fp_SPEC_DEFINE_TLOC_MOD(_name_, _tl_, _s_...) \
  fcs_near_fp_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_) \
  fcs_near_fp_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro fcs_near_fp_SPEC_EXT_PARAM_TLOC fcs_near_fp_SPEC_EXT_PARAM_TLOC_NULL fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD_NULL */
#define fcs_near_fp_SPEC_EXT_PARAM_TLOC(_name_)      _name_##_tloc_rearrange_db, _name_##_tloc_rearrange_ip
#define fcs_near_fp_SPEC_EXT_PARAM_TLOC_NULL         NULL, NULL
#define fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD(_name_)  _name_##_tloc_mod_rearrange_db, _name_##_tloc_mod_rearrange_ip
#define fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD_NULL     NULL, NULL


/* sp_type fcs_near_fp_spec_tloc_f fcs_near_fp_spec_tloc_rearrange_db_f fcs_near_fp_spec_tloc_rearrange_ip_f */
typedef fcs_near_fp_spec_elem_index_t fcs_near_fp_spec_tloc_f(fcs_near_fp_spec_elem_buf_t b, fcs_near_fp_spec_elem_index_t x, fcs_near_fp_spec_tloc_data_t tloc_data);
typedef void fcs_near_fp_spec_tloc_rearrange_db_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tloc_data_t tloc_data);
typedef void fcs_near_fp_spec_tloc_rearrange_ip_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tloc_data_t tloc_data);

/* sp_type fcs_near_fp_spec_tloc_mod_f fcs_near_fp_spec_tloc_mod_rearrange_db_f fcs_near_fp_spec_tloc_mod_rearrange_ip_f */
typedef fcs_near_fp_spec_elem_index_t fcs_near_fp_spec_tloc_mod_f(fcs_near_fp_spec_elem_buf_t b, fcs_near_fp_spec_elem_index_t x, fcs_near_fp_spec_tloc_data_t tloc_data, fcs_near_fp_spec_elem_buf_t mod);
typedef void fcs_near_fp_spec_tloc_mod_rearrange_db_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *d, fcs_near_fp_spec_tloc_data_t tloc_data, fcs_near_fp_spec_elem_t *mod);
typedef void fcs_near_fp_spec_tloc_mod_rearrange_ip_f(fcs_near_fp_spec_elem_t *s, fcs_near_fp_spec_elem_t *x, fcs_near_fp_spec_tloc_data_t tloc_data, fcs_near_fp_spec_elem_t *mod);


#endif /* fcs_near_fp_SPEC_TLOC */






#ifdef SL_USE_MPI
# include <mpi.h>
#endif


/* sl_type fcs_near_fp_slint_t fcs_near_fp_slint */
typedef fcs_near_fp_sl_int_type_c fcs_near_fp_slint_t, fcs_near_fp_slint;  /* deprecated 'fcs_near_fp_slint' */

#define fcs_near_fp_slint_fmt   fcs_near_fp_sl_int_type_fmt    /* sl_macro */

/* sl_type fcs_near_fp_slindex_t */
typedef fcs_near_fp_sl_index_type_c fcs_near_fp_slindex_t;

#define fcs_near_fp_sindex_fmt  fcs_near_fp_sl_index_type_fmt  /* sl_macro */

/* sl_type fcs_near_fp_slkey_t */
typedef fcs_near_fp_sl_key_type_c fcs_near_fp_slkey_t;

/* sl_type fcs_near_fp_slkey_pure_t fcs_near_fp_slpkey_t */
typedef fcs_near_fp_sl_key_pure_type_c fcs_near_fp_slkey_pure_t, fcs_near_fp_slpkey_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type fcs_near_fp_sldata0_t */
#ifdef fcs_near_fp_sl_data0_type_c
typedef fcs_near_fp_sl_data0_type_c fcs_near_fp_sldata0_t;
#endif
/* sl_type fcs_near_fp_sldata1_t */
#ifdef fcs_near_fp_sl_data1_type_c
typedef fcs_near_fp_sl_data1_type_c fcs_near_fp_sldata1_t;
#endif
/* sl_type fcs_near_fp_sldata2_t */
#ifdef fcs_near_fp_sl_data2_type_c
typedef fcs_near_fp_sl_data2_type_c fcs_near_fp_sldata2_t;
#endif
/* sl_type fcs_near_fp_sldata3_t */
#ifdef fcs_near_fp_sl_data3_type_c
typedef fcs_near_fp_sl_data3_type_c fcs_near_fp_sldata3_t;
#endif
/* sl_type fcs_near_fp_sldata4_t */
#ifdef fcs_near_fp_sl_data4_type_c
typedef fcs_near_fp_sl_data4_type_c fcs_near_fp_sldata4_t;
#endif
/* sl_type fcs_near_fp_sldata5_t */
#ifdef fcs_near_fp_sl_data5_type_c
typedef fcs_near_fp_sl_data5_type_c fcs_near_fp_sldata5_t;
#endif
/* sl_type fcs_near_fp_sldata6_t */
#ifdef fcs_near_fp_sl_data6_type_c
typedef fcs_near_fp_sl_data6_type_c fcs_near_fp_sldata6_t;
#endif
/* sl_type fcs_near_fp_sldata7_t */
#ifdef fcs_near_fp_sl_data7_type_c
typedef fcs_near_fp_sl_data7_type_c fcs_near_fp_sldata7_t;
#endif
/* sl_type fcs_near_fp_sldata8_t */
#ifdef fcs_near_fp_sl_data8_type_c
typedef fcs_near_fp_sl_data8_type_c fcs_near_fp_sldata8_t;
#endif
/* sl_type fcs_near_fp_sldata9_t */
#ifdef fcs_near_fp_sl_data9_type_c
typedef fcs_near_fp_sl_data9_type_c fcs_near_fp_sldata9_t;
#endif
/* sl_type fcs_near_fp_sldata10_t */
#ifdef fcs_near_fp_sl_data10_type_c
typedef fcs_near_fp_sl_data10_type_c fcs_near_fp_sldata10_t;
#endif
/* sl_type fcs_near_fp_sldata11_t */
#ifdef fcs_near_fp_sl_data11_type_c
typedef fcs_near_fp_sl_data11_type_c fcs_near_fp_sldata11_t;
#endif
/* sl_type fcs_near_fp_sldata12_t */
#ifdef fcs_near_fp_sl_data12_type_c
typedef fcs_near_fp_sl_data12_type_c fcs_near_fp_sldata12_t;
#endif
/* sl_type fcs_near_fp_sldata13_t */
#ifdef fcs_near_fp_sl_data13_type_c
typedef fcs_near_fp_sl_data13_type_c fcs_near_fp_sldata13_t;
#endif
/* sl_type fcs_near_fp_sldata14_t */
#ifdef fcs_near_fp_sl_data14_type_c
typedef fcs_near_fp_sl_data14_type_c fcs_near_fp_sldata14_t;
#endif
/* sl_type fcs_near_fp_sldata15_t */
#ifdef fcs_near_fp_sl_data15_type_c
typedef fcs_near_fp_sl_data15_type_c fcs_near_fp_sldata15_t;
#endif
/* sl_type fcs_near_fp_sldata16_t */
#ifdef fcs_near_fp_sl_data16_type_c
typedef fcs_near_fp_sl_data16_type_c fcs_near_fp_sldata16_t;
#endif
/* sl_type fcs_near_fp_sldata17_t */
#ifdef fcs_near_fp_sl_data17_type_c
typedef fcs_near_fp_sl_data17_type_c fcs_near_fp_sldata17_t;
#endif
/* sl_type fcs_near_fp_sldata18_t */
#ifdef fcs_near_fp_sl_data18_type_c
typedef fcs_near_fp_sl_data18_type_c fcs_near_fp_sldata18_t;
#endif
/* sl_type fcs_near_fp_sldata19_t */
#ifdef fcs_near_fp_sl_data19_type_c
typedef fcs_near_fp_sl_data19_type_c fcs_near_fp_sldata19_t;
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

/* sl_type fcs_near_fp_slweight_t */
typedef fcs_near_fp_sl_weight_type_c fcs_near_fp_slweight_t;

#define fcs_near_fp_slweight_fmt  fcs_near_fp_sl_weight_type_fmt  /* sl_macro */

#if defined(fcs_near_fp_sl_elem_weight) && defined(fcs_near_fp_sl_weight_intequiv)
typedef fcs_near_fp_sl_weight_type_c fcs_near_fp_slcount_t;       /* sl_type fcs_near_fp_slcount_t */
# define fcs_near_fp_slcount_fmt  fcs_near_fp_sl_weight_type_fmt  /* sl_macro */
#else
typedef fcs_near_fp_sl_int_type_c fcs_near_fp_slcount_t;
# define fcs_near_fp_slcount_fmt  fcs_near_fp_sl_int_type_fmt
#endif


/* sl_type fcs_near_fp__slpwkey_t fcs_near_fp_slpwkey_t */
typedef struct fcs_near_fp__slpwkey_t
{
  fcs_near_fp_slpkey_t pkey;
  fcs_near_fp_slweight_t weight;

} fcs_near_fp_slpwkey_t;


/* sl_type fcs_near_fp__elements_t fcs_near_fp_elements_t */
typedef struct fcs_near_fp__elements_t
{
  fcs_near_fp_slint_t size, max_size;
  fcs_near_fp_slkey_t *keys;

#ifdef fcs_near_fp_SL_INDEX
  fcs_near_fp_slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef fcs_near_fp_SL_DATA0
  fcs_near_fp_sldata0_t *data0;
#endif
#ifdef fcs_near_fp_SL_DATA1
  fcs_near_fp_sldata1_t *data1;
#endif
#ifdef fcs_near_fp_SL_DATA2
  fcs_near_fp_sldata2_t *data2;
#endif
#ifdef fcs_near_fp_SL_DATA3
  fcs_near_fp_sldata3_t *data3;
#endif
#ifdef fcs_near_fp_SL_DATA4
  fcs_near_fp_sldata4_t *data4;
#endif
#ifdef fcs_near_fp_SL_DATA5
  fcs_near_fp_sldata5_t *data5;
#endif
#ifdef fcs_near_fp_SL_DATA6
  fcs_near_fp_sldata6_t *data6;
#endif
#ifdef fcs_near_fp_SL_DATA7
  fcs_near_fp_sldata7_t *data7;
#endif
#ifdef fcs_near_fp_SL_DATA8
  fcs_near_fp_sldata8_t *data8;
#endif
#ifdef fcs_near_fp_SL_DATA9
  fcs_near_fp_sldata9_t *data9;
#endif
#ifdef fcs_near_fp_SL_DATA10
  fcs_near_fp_sldata10_t *data10;
#endif
#ifdef fcs_near_fp_SL_DATA11
  fcs_near_fp_sldata11_t *data11;
#endif
#ifdef fcs_near_fp_SL_DATA12
  fcs_near_fp_sldata12_t *data12;
#endif
#ifdef fcs_near_fp_SL_DATA13
  fcs_near_fp_sldata13_t *data13;
#endif
#ifdef fcs_near_fp_SL_DATA14
  fcs_near_fp_sldata14_t *data14;
#endif
#ifdef fcs_near_fp_SL_DATA15
  fcs_near_fp_sldata15_t *data15;
#endif
#ifdef fcs_near_fp_SL_DATA16
  fcs_near_fp_sldata16_t *data16;
#endif
#ifdef fcs_near_fp_SL_DATA17
  fcs_near_fp_sldata17_t *data17;
#endif
#ifdef fcs_near_fp_SL_DATA18
  fcs_near_fp_sldata18_t *data18;
#endif
#ifdef fcs_near_fp_SL_DATA19
  fcs_near_fp_sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} fcs_near_fp_elements_t;


/* sl_type fcs_near_fp__packed_element_t fcs_near_fp_packed_element_t */
typedef struct fcs_near_fp__packed_element_t
{
  fcs_near_fp_slkey_t key;

#ifdef fcs_near_fp_SL_PACKED_INDEX
  fcs_near_fp_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef fcs_near_fp_SL_DATA0
# ifdef fcs_near_fp_sl_data0_flex
  fcs_near_fp_sldata0_t data0[];
# else
  fcs_near_fp_sldata0_t data0[fcs_near_fp_sl_data0_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA1
# ifdef fcs_near_fp_sl_data1_flex
  fcs_near_fp_sldata1_t data1[];
# else
  fcs_near_fp_sldata1_t data1[fcs_near_fp_sl_data1_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA2
# ifdef fcs_near_fp_sl_data2_flex
  fcs_near_fp_sldata2_t data2[];
# else
  fcs_near_fp_sldata2_t data2[fcs_near_fp_sl_data2_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA3
# ifdef fcs_near_fp_sl_data3_flex
  fcs_near_fp_sldata3_t data3[];
# else
  fcs_near_fp_sldata3_t data3[fcs_near_fp_sl_data3_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA4
# ifdef fcs_near_fp_sl_data4_flex
  fcs_near_fp_sldata4_t data4[];
# else
  fcs_near_fp_sldata4_t data4[fcs_near_fp_sl_data4_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA5
# ifdef fcs_near_fp_sl_data5_flex
  fcs_near_fp_sldata5_t data5[];
# else
  fcs_near_fp_sldata5_t data5[fcs_near_fp_sl_data5_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA6
# ifdef fcs_near_fp_sl_data6_flex
  fcs_near_fp_sldata6_t data6[];
# else
  fcs_near_fp_sldata6_t data6[fcs_near_fp_sl_data6_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA7
# ifdef fcs_near_fp_sl_data7_flex
  fcs_near_fp_sldata7_t data7[];
# else
  fcs_near_fp_sldata7_t data7[fcs_near_fp_sl_data7_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA8
# ifdef fcs_near_fp_sl_data8_flex
  fcs_near_fp_sldata8_t data8[];
# else
  fcs_near_fp_sldata8_t data8[fcs_near_fp_sl_data8_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA9
# ifdef fcs_near_fp_sl_data9_flex
  fcs_near_fp_sldata9_t data9[];
# else
  fcs_near_fp_sldata9_t data9[fcs_near_fp_sl_data9_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA10
# ifdef fcs_near_fp_sl_data10_flex
  fcs_near_fp_sldata10_t data10[];
# else
  fcs_near_fp_sldata10_t data10[fcs_near_fp_sl_data10_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA11
# ifdef fcs_near_fp_sl_data11_flex
  fcs_near_fp_sldata11_t data11[];
# else
  fcs_near_fp_sldata11_t data11[fcs_near_fp_sl_data11_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA12
# ifdef fcs_near_fp_sl_data12_flex
  fcs_near_fp_sldata12_t data12[];
# else
  fcs_near_fp_sldata12_t data12[fcs_near_fp_sl_data12_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA13
# ifdef fcs_near_fp_sl_data13_flex
  fcs_near_fp_sldata13_t data13[];
# else
  fcs_near_fp_sldata13_t data13[fcs_near_fp_sl_data13_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA14
# ifdef fcs_near_fp_sl_data14_flex
  fcs_near_fp_sldata14_t data14[];
# else
  fcs_near_fp_sldata14_t data14[fcs_near_fp_sl_data14_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA15
# ifdef fcs_near_fp_sl_data15_flex
  fcs_near_fp_sldata15_t data15[];
# else
  fcs_near_fp_sldata15_t data15[fcs_near_fp_sl_data15_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA16
# ifdef fcs_near_fp_sl_data16_flex
  fcs_near_fp_sldata16_t data16[];
# else
  fcs_near_fp_sldata16_t data16[fcs_near_fp_sl_data16_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA17
# ifdef fcs_near_fp_sl_data17_flex
  fcs_near_fp_sldata17_t data17[];
# else
  fcs_near_fp_sldata17_t data17[fcs_near_fp_sl_data17_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA18
# ifdef fcs_near_fp_sl_data18_flex
  fcs_near_fp_sldata18_t data18[];
# else
  fcs_near_fp_sldata18_t data18[fcs_near_fp_sl_data18_size_c];
# endif
#endif
#ifdef fcs_near_fp_SL_DATA19
# ifdef fcs_near_fp_sl_data19_flex
  fcs_near_fp_sldata19_t data19[];
# else
  fcs_near_fp_sldata19_t data19[fcs_near_fp_sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} fcs_near_fp_packed_element_t;


/* sl_type fcs_near_fp__packed_elements_t fcs_near_fp_packed_elements_t */
typedef struct fcs_near_fp__packed_elements_t
{
  fcs_near_fp_slint_t size, max_size;
  
  fcs_near_fp_packed_element_t *elements;
  
} fcs_near_fp_packed_elements_t;


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


/* sl_type fcs_near_fp__classification_info_t fcs_near_fp_classification_info_t fcs_near_fp_classification_info */
typedef struct fcs_near_fp__classification_info_t
{
  fcs_near_fp_slint_t nclasses;
  fcs_near_fp_slkey_pure_t *keys;
  fcs_near_fp_slint_t *counts;
  fcs_near_fp_slint_t *masks;

  /* */
  fcs_near_fp_slint_t *all_local_sizes;
  fcs_near_fp_slint_t *local_lt_eq_counts;
  fcs_near_fp_slint_t *all_local_lt_eq_counts;

} fcs_near_fp_classification_info_t, fcs_near_fp_classification_info;  /* deprecated 'fcs_near_fp_classification_info' */


/* key2class, sl_type fcs_near_fp_key2class_f */
typedef fcs_near_fp_slint_t (*fcs_near_fp_key2class_f)(fcs_near_fp_slkey_t *, fcs_near_fp_slint, void *);

/* pivot-element, sl_type fcs_near_fp_pivot_f */
typedef fcs_near_fp_slint_t (*fcs_near_fp_pivot_f)(fcs_near_fp_elements_t *);

/* sorting-network, sl_type fcs_near_fp_sortnet_f fcs_near_fp_sortnet_data_t */
typedef void *fcs_near_fp_sortnet_data_t;
typedef fcs_near_fp_slint_t (*fcs_near_fp_sortnet_f)(fcs_near_fp_slint_t size, fcs_near_fp_slint_t rank, fcs_near_fp_slint_t stage, fcs_near_fp_sortnet_data_t snd, fcs_near_fp_slint_t *up);

/* merge2, sl_type fcs_near_fp_merge2x_f fcs_near_fp_merge2X_f */
typedef fcs_near_fp_slint_t (*fcs_near_fp_merge2x_f)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
typedef fcs_near_fp_slint_t (*fcs_near_fp_merge2X_f)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx, fcs_near_fp_elements_t *t);

/* sl_type fcs_near_fp__permute_generic_t fcs_near_fp_permute_generic_t */
typedef struct fcs_near_fp__permute_generic_t
{
  int type;

  fcs_near_fp_spec_tloc_f *tloc;
  fcs_near_fp_spec_tloc_rearrange_db_f *tloc_rearrange_db;
  fcs_near_fp_spec_tloc_rearrange_ip_f *tloc_rearrange_ip;

  fcs_near_fp_spec_tloc_mod_f *tloc_mod;
  fcs_near_fp_spec_tloc_mod_rearrange_db_f *tloc_mod_rearrange_db;
  fcs_near_fp_spec_tloc_mod_rearrange_ip_f *tloc_mod_rearrange_ip;

} fcs_near_fp_permute_generic_t;

/* sl_macro fcs_near_fp_PERMUTE_GENERIC_DEFINE_TLOC fcs_near_fp_PERMUTE_GENERIC_INIT_TLOC fcs_near_fp_PERMUTE_GENERIC_INIT_EXT_TLOC */
#define fcs_near_fp_PERMUTE_GENERIC_DEFINE_TLOC(_tl_, _s_...)      fcs_near_fp_SPEC_DEFINE_TLOC(_tl_, _tl_, _s_)
#define fcs_near_fp_PERMUTE_GENERIC_INIT_TLOC(_tl_)                { 1, _tl_, fcs_near_fp_SPEC_EXT_PARAM_TLOC_NULL,  NULL, fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define fcs_near_fp_PERMUTE_GENERIC_INIT_EXT_TLOC(_tl_)            { 1, _tl_, fcs_near_fp_SPEC_EXT_PARAM_TLOC(_tl_), NULL, fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD_NULL }

/* sl_macro fcs_near_fp_PERMUTE_GENERIC_DEFINE_TLOC_MOD fcs_near_fp_PERMUTE_GENERIC_INIT_TLOC_MOD fcs_near_fp_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD */
#define fcs_near_fp_PERMUTE_GENERIC_DEFINE_TLOC_MOD(_tl_, _s_...)  fcs_near_fp_SPEC_DEFINE_TLOC_MOD(_tl_, _tl_, _s_)
#define fcs_near_fp_PERMUTE_GENERIC_INIT_TLOC_MOD(_tl_)            { 2, NULL, fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define fcs_near_fp_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD(_tl_)        { 2, NULL, fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, fcs_near_fp_SPEC_EXT_PARAM_TLOC_MOD(_tl_) }

/* sl_type fcs_near_fp__split_generic_t fcs_near_fp_split_generic_t */
typedef struct fcs_near_fp__split_generic_t
{
  int type;

  fcs_near_fp_spec_tproc_f *tproc;
  fcs_near_fp_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  fcs_near_fp_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  fcs_near_fp_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  fcs_near_fp_spec_tproc_mod_f *tproc_mod;
  fcs_near_fp_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  fcs_near_fp_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  fcs_near_fp_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  fcs_near_fp_spec_tprocs_f *tprocs;
  fcs_near_fp_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  fcs_near_fp_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  fcs_near_fp_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  fcs_near_fp_spec_tprocs_mod_f *tprocs_mod;
  fcs_near_fp_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  fcs_near_fp_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  fcs_near_fp_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

  fcs_near_fp_spec_tproc_reset_f *reset;

} fcs_near_fp_split_generic_t;

/* sl_macro fcs_near_fp_SPLIT_GENERIC_DEFINE_TPROC fcs_near_fp_SPLIT_GENERIC_INIT_TPROC fcs_near_fp_SPLIT_GENERIC_INIT_EXT_TPROC */
#define fcs_near_fp_SPLIT_GENERIC_DEFINE_TPROC(_tp_, _s_...)         fcs_near_fp_SPEC_DEFINE_TPROC(_tp_, _tp_, _s_)
#define fcs_near_fp_SPLIT_GENERIC_INIT_TPROC(_tp_, _r_...)           { 1, _tp_, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL,  NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define fcs_near_fp_SPLIT_GENERIC_INIT_EXT_TPROC(_tp_, _r_...)       { 1, _tp_, fcs_near_fp_SPEC_EXT_PARAM_TPROC(_tp_), NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro fcs_near_fp_SPLIT_GENERIC_DEFINE_TPROC_MOD fcs_near_fp_SPLIT_GENERIC_INIT_TPROC_MOD fcs_near_fp_SPLIT_GENERIC_INIT_EXT_TPROC_MOD */
#define fcs_near_fp_SPLIT_GENERIC_DEFINE_TPROC_MOD(_tp_, _s_...)     fcs_near_fp_SPEC_DEFINE_TPROC_MOD(_tp_, _tp_, _s_)
#define fcs_near_fp_SPLIT_GENERIC_INIT_TPROC_MOD(_tp_, _r_...)       { 2, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL, _tp_, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL,  NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define fcs_near_fp_SPLIT_GENERIC_INIT_EXT_TPROC_MOD(_tp_, _r_...)   { 2, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL, _tp_, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD(_tp_), NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro fcs_near_fp_SPLIT_GENERIC_DEFINE_TPROCS fcs_near_fp_SPLIT_GENERIC_INIT_TPROCS fcs_near_fp_SPLIT_GENERIC_INIT_EXT_TPROCS */
#define fcs_near_fp_SPLIT_GENERIC_DEFINE_TPROCS(_tp_, _s_...)        fcs_near_fp_SPEC_DEFINE_TPROCS(_tp_, _tp_, _s_)
#define fcs_near_fp_SPLIT_GENERIC_INIT_TPROCS(_tp_, _r_...)          { 3, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL,  NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define fcs_near_fp_SPLIT_GENERIC_INIT_EXT_TPROCS(_tp_, _r_...)      { 3, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, fcs_near_fp_SPEC_EXT_PARAM_TPROCS(_tp_), NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro fcs_near_fp_SPLIT_GENERIC_DEFINE_TPROCS_MOD fcs_near_fp_SPLIT_GENERIC_INIT_TPROCS_MOD fcs_near_fp_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD */
#define fcs_near_fp_SPLIT_GENERIC_DEFINE_TPROCS_MOD(_tp_, _s_...)    fcs_near_fp_SPEC_DEFINE_TPROCS_MOD(_tp_, _tp_, _s_)
#define fcs_near_fp_SPLIT_GENERIC_INIT_TPROCS_MOD(_tp_, _r_...)      { 4, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define fcs_near_fp_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD(_tp_, _r_...)  { 4, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD(_tp_), _r_ }

/* sl_type fcs_near_fp_tloc_f fcs_near_fp_tloc_mod_f */
typedef fcs_near_fp_slint_t fcs_near_fp_tloc_f(fcs_near_fp_elements_t *b, fcs_near_fp_slint_t x, void *tloc_data);
typedef fcs_near_fp_slint_t fcs_near_fp_tloc_mod_f(fcs_near_fp_elements_t *b, fcs_near_fp_slint_t x, void *tloc_data, fcs_near_fp_elements_t *mod);

/* sl_type fcs_near_fp_tproc_f fcs_near_fp_tproc_mod_f fcs_near_fp_tprocs_f fcs_near_fp_tprocs_mod_f */
typedef int fcs_near_fp_tproc_f(fcs_near_fp_elements_t *b, fcs_near_fp_slint_t x, void *tproc_data);
typedef int fcs_near_fp_tproc_mod_f(fcs_near_fp_elements_t *b, fcs_near_fp_slint_t x, void *tproc_data, fcs_near_fp_elements_t *mod);
typedef fcs_near_fp_slint_t fcs_near_fp_tprocs_f(fcs_near_fp_elements_t *b, fcs_near_fp_slint_t x, void *tproc_data, int *procs);
typedef fcs_near_fp_slint_t fcs_near_fp_tprocs_mod_f(fcs_near_fp_elements_t *b, fcs_near_fp_slint_t x, void *tproc_data, int *procs, fcs_near_fp_elements_t *mod);

/* sl_type fcs_near_fp_tproc_reset_f */
typedef void fcs_near_fp_tproc_reset_f(void *tproc_data);

/* sl_macro fcs_near_fp_TPROC_RESET_NULL */
#define fcs_near_fp_TPROC_RESET_NULL  NULL

/* sl_type fcs_near_fp__tproc_t fcs_near_fp_tproc_t */
typedef struct fcs_near_fp__tproc_t *fcs_near_fp_tproc_t;

/* sl_type fcs_near_fp__tproc_exdef fcs_near_fp_tproc_exdef */
typedef struct fcs_near_fp__tproc_exdef {
  int type;

  fcs_near_fp_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  fcs_near_fp_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  fcs_near_fp_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  fcs_near_fp_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  fcs_near_fp_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  fcs_near_fp_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  fcs_near_fp_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  fcs_near_fp_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  fcs_near_fp_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  fcs_near_fp_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  fcs_near_fp_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  fcs_near_fp_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

} const *fcs_near_fp_tproc_exdef;

/* sl_macro fcs_near_fp_TPROC_EXDEF_NULL */
#define fcs_near_fp_TPROC_EXDEF_NULL  NULL

/* sl_macro fcs_near_fp_TPROC_EXDEF_DEFINE_TPROC fcs_near_fp_TPROC_EXDEF_DEFINE_TPROC_MOD fcs_near_fp_TPROC_EXDEF_DEFINE_TPROCS fcs_near_fp_TPROC_EXDEF_DEFINE_TPROCS_MOD */
#define fcs_near_fp_TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  fcs_near_fp_SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct fcs_near_fp__tproc_exdef _##_name_ = { 1, fcs_near_fp_SPEC_EXT_PARAM_TPROC(_name_), fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define fcs_near_fp_TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  fcs_near_fp_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct fcs_near_fp__tproc_exdef _##_name_ = { 2, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD(_name_), fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define fcs_near_fp_TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  fcs_near_fp_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct fcs_near_fp__tproc_exdef _##_name_ = { 3, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS(_name_), fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define fcs_near_fp_TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  fcs_near_fp_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct fcs_near_fp__tproc_exdef _##_name_ = { 4, fcs_near_fp_SPEC_EXT_PARAM_TPROC_NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROC_MOD_NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_NULL, fcs_near_fp_SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;


/* deprecated, sl_type fcs_near_fp_k2c_func fcs_near_fp_pivot_func fcs_near_fp_sn_func fcs_near_fp_m2x_func fcs_near_fp_m2X_func */
typedef fcs_near_fp_key2class_f fcs_near_fp_k2c_func;
typedef fcs_near_fp_pivot_f fcs_near_fp_pivot_func;
typedef fcs_near_fp_sortnet_f fcs_near_fp_sn_func;
typedef fcs_near_fp_merge2x_f fcs_near_fp_m2x_func;
typedef fcs_near_fp_merge2X_f fcs_near_fp_m2X_func;


/* sl_type fcs_near_fp__mergek_t fcs_near_fp_mergek_t */
typedef struct fcs_near_fp__mergek_t
{
  fcs_near_fp_sortnet_f sn;
  fcs_near_fp_sortnet_data_t snd;

  fcs_near_fp_merge2x_f m2x;
  fcs_near_fp_elements_t *sx;

} fcs_near_fp_mergek_t;


/* sl_type fcs_near_fp_keys_init_type_t fcs_near_fp_keys_init_data_t */
typedef fcs_near_fp_slint_t fcs_near_fp_keys_init_type_t;
typedef void *fcs_near_fp_keys_init_data_t;

/* sl_type fcs_near_fp_key_set_data_t fcs_near_fp_key_set_f */
typedef void *fcs_near_fp_key_set_data_t;
typedef void (*fcs_near_fp_key_set_f)(fcs_near_fp_slkey_pure_t *k, fcs_near_fp_key_set_data_t d);


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


/* fcs_near_fp_elements_keys_stats */
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


/* partition conditions, sl_type fcs_near_fp__partcond2_t fcs_near_fp_partcond2_t */
typedef struct fcs_near_fp__partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} fcs_near_fp_partcond2_t;


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

/* partition conditions, sl_type fcs_near_fp__partcond_t fcs_near_fp_partcond_t fcs_near_fp_partcond_p */
typedef struct fcs_near_fp__partcond_t
{
  fcs_near_fp_slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} fcs_near_fp_partcond_t, *fcs_near_fp_partcond_p;


/* internal partition conditions, sl_type fcs_near_fp__partcond_intern_t fcs_near_fp_partcond_intern_t fcs_near_fp_partcond_intern_p */
typedef struct fcs_near_fp__partcond_intern_t
{
  fcs_near_fp_slint_t pcm;
  fcs_near_fp_slint_t count_min, count_max;
  fcs_near_fp_slint_t count_low, count_high;
#ifdef elem_weight
  fcs_near_fp_slweight_t weight_min, weight_max;
  fcs_near_fp_slweight_t weight_low, weight_high;
#endif

} fcs_near_fp_partcond_intern_t, *fcs_near_fp_partcond_intern_p;


/* sl_type fcs_near_fp__parttype_t fcs_near_fp_parttype_t fcs_near_fp_parttype_p */
typedef struct fcs_near_fp__parttype_t
{
  fcs_near_fp_slint_t type;

} fcs_near_fp_parttype_t, *fcs_near_fp_parttype_p;


/* generic binning method */

/* sl_type fcs_near_fp__bin_t fcs_near_fp_bin_t */
typedef struct fcs_near_fp__bin_t
{
  fcs_near_fp_elements_t s;

#ifdef elem_weight
  fcs_near_fp_slweight_t weight;
#endif

} fcs_near_fp_bin_t;


/* sl_type fcs_near_fp__splitter_t fcs_near_fp_splitter_t */
typedef struct fcs_near_fp__splitter_t
{
  fcs_near_fp_slint_t n;

  int *displs;
  fcs_near_fp_slkey_pure_t *s;
  fcs_near_fp_slint_t *sn;

} fcs_near_fp_splitter_t;


struct fcs_near_fp__binning_t;

/* sl_type fcs_near_fp_binning_pre_f fcs_near_fp_binning_exec_f fcs_near_fp_binning_refine_f fcs_near_fp_binning_hit_f fcs_near_fp_binning_finalize_f fcs_near_fp_binning_post_f */
typedef fcs_near_fp_slint_t (*fcs_near_fp_binning_pre_f)(struct fcs_near_fp__binning_t *bm);
typedef fcs_near_fp_slint_t (*fcs_near_fp_binning_exec_f)(struct fcs_near_fp__binning_t *bm, fcs_near_fp_bin_t *bin, fcs_near_fp_slcount_t *counts, fcs_near_fp_slweight_t *weights);
typedef fcs_near_fp_slint_t (*fcs_near_fp_binning_refine_f)(struct fcs_near_fp__binning_t *bm, fcs_near_fp_bin_t *bin, fcs_near_fp_slint_t k, fcs_near_fp_slcount_t *counts, fcs_near_fp_slweight_t *weights, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s, fcs_near_fp_bin_t *new_bin);
typedef fcs_near_fp_slint_t (*fcs_near_fp_binning_hit_f)(struct fcs_near_fp__binning_t *bm, fcs_near_fp_bin_t *bin, fcs_near_fp_slint_t k, fcs_near_fp_slcount_t *counts, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s);
typedef fcs_near_fp_slint_t (*fcs_near_fp_binning_finalize_f)(struct fcs_near_fp__binning_t *bm, fcs_near_fp_bin_t *bin, fcs_near_fp_slint_t dc, fcs_near_fp_slweight_t dw, fcs_near_fp_slint_t lc_min, fcs_near_fp_slint_t lc_max, fcs_near_fp_slcount_t *lcs, fcs_near_fp_slweight_t *lws, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s);
typedef fcs_near_fp_slint_t (*fcs_near_fp_binning_post_f)(struct fcs_near_fp__binning_t *bm);


/* sl_type fcs_near_fp__binning_data_t fcs_near_fp_binning_data_t */
typedef union fcs_near_fp__binning_data_t
{
  struct
  {
    fcs_near_fp_slint_t rhigh, rlow, rwidth;
    fcs_near_fp_slint_t rcurrent;
    fcs_near_fp_slkey_pure_t bit_mask;

    fcs_near_fp_elements_t sx;

  } radix;

} fcs_near_fp_binning_data_t;


/* sl_type fcs_near_fp__binning_t fcs_near_fp_binning_t */
typedef struct fcs_near_fp__binning_t
{
  fcs_near_fp_slint_t nbins, max_nbins;
  
  fcs_near_fp_binning_pre_f pre;
  fcs_near_fp_binning_exec_f exec;
  fcs_near_fp_binning_refine_f refine;
  fcs_near_fp_binning_hit_f hit;
  fcs_near_fp_binning_finalize_f finalize;
  fcs_near_fp_binning_post_f post;

  fcs_near_fp_slint_t sorted;

  fcs_near_fp_slint_t docounts;
#ifdef elem_weight
  fcs_near_fp_slint_t doweights;
#endif

  fcs_near_fp_binning_data_t bd;

} fcs_near_fp_binning_t;


/* sl_type fcs_near_fp__local_bins_t fcs_near_fp_local_bins_t */
typedef struct fcs_near_fp__local_bins_t
{
  fcs_near_fp_binning_t *bm;

  fcs_near_fp_slint_t nbins, max_nbins;
  fcs_near_fp_slint_t nelements;

  fcs_near_fp_slint_t docounts;
#ifdef elem_weight
  fcs_near_fp_slint_t doweights;
#endif

  fcs_near_fp_slint_t nbinnings, max_nbinnings;

  fcs_near_fp_slint_t nbins_new, last_new_b, last_new_k;
  fcs_near_fp_bin_t *bins, *bins_new;
  fcs_near_fp_bin_t *bins0, *bins1;

  fcs_near_fp_slint_t *bcws;

#if defined(elem_weight) && defined(fcs_near_fp_sl_weight_intequiv)
  fcs_near_fp_slint_t cw_factor, w_index, bin_cw_factor;
  fcs_near_fp_slweight_t *cws, *bin_cws;
  fcs_near_fp_slweight_t *prefix_cws;
#else
  fcs_near_fp_slint_t *cs, *bin_cs;
  fcs_near_fp_slint_t *prefix_cs;
# ifdef elem_weight
  fcs_near_fp_slweight_t *ws, *bin_ws;
  fcs_near_fp_slweight_t *prefix_ws;
# endif
#endif

  fcs_near_fp_slint_t last_exec_b;

} fcs_near_fp_local_bins_t;


/* sl_type fcs_near_fp__global_bins_t fcs_near_fp_global_bins_t */
typedef struct fcs_near_fp__global_bins_t
{
  fcs_near_fp_binning_t *bm;
  
  fcs_near_fp_local_bins_t lb;

  fcs_near_fp_slint_t *bcws;

#if defined(elem_weight) && defined(fcs_near_fp_sl_weight_intequiv)
  fcs_near_fp_slweight_t *cws;
  fcs_near_fp_slweight_t *prefix_cws;
#else
  fcs_near_fp_slint_t *cs;
  fcs_near_fp_slint_t *prefix_cs;
# ifdef elem_weight
  fcs_near_fp_slweight_t *ws;
  fcs_near_fp_slweight_t *prefix_ws;
# endif
#endif

} fcs_near_fp_global_bins_t;


/* sl_type fcs_near_fp_rti_cmc_t */
typedef struct
{
  fcs_near_fp_slint_t cmp, movek, moved;

} fcs_near_fp_rti_cmc_t;

#ifndef my_rti_ccmp
# define my_rti_ccmp(m)    m.cmc.cmp
# define my_rti_cmovek(m)  m.cmc.movek
# define my_rti_cmoved(m)  m.cmc.moved
#endif


/* sl_type fcs_near_fp_rti_tim_t */
typedef struct
{
  double start, stop;
  double last, cumu;

  fcs_near_fp_slint_t num;

} fcs_near_fp_rti_tim_t[rti_tids];

#ifndef my_rti_tlast
# define my_rti_tlast(m, t)  m.tim[t].last
# define my_rti_tcumu(m, t)  m.tim[t].cumu
# define my_rti_tnum(m, t)   m.tim[t].num
#endif


/* sl_type fcs_near_fp_rti_mem_t */
typedef struct
{
  fcs_near_fp_slint_t nalloc, nfree;
  fcs_near_fp_slint_t max, cur, cur_max;

} fcs_near_fp_rti_mem_t;


/* sl_type fcs_near_fp_rti_t */
typedef struct
{
  /* compare-move-counter */
  fcs_near_fp_rti_cmc_t cmc;
  /* timer */
  fcs_near_fp_rti_tim_t tim;
  /* memory */
  fcs_near_fp_rti_mem_t mem;

} fcs_near_fp_rti_t;

#ifndef my_rti_reset
# define my_rti_reset(m)  memset((void *) &m, 0, sizeof(m))
#endif


/* sl_type fcs_near_fp__sl_context_t fcs_near_fp_sl_context_t */
typedef struct fcs_near_fp__sl_context_t
{

/* src/base/base.c */
  struct {
int dummy_rank;
  } sl;
#ifdef fcs_near_fp_SL_USE_RTI
fcs_near_fp_rti_t rti;
#endif
  struct {
fcs_near_fp_slint_t ip_threshold;
fcs_near_fp_slint_t db_threshold;
fcs_near_fp_slint_t ma_threshold;
  } sr;
  struct {
fcs_near_fp_slint_t threshold;
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
fcs_near_fp_slint_t sendrecv_replace_memsize;
fcs_near_fp_slint_t sendrecv_replace_mpi_maxsize;
  } me;
#endif
#ifdef SL_USE_MPI
  struct {
double t[2];
fcs_near_fp_slint_t max_nprocs;
fcs_near_fp_slint_t packed;
fcs_near_fp_slint_t minalloc;
double overalloc;
  } meas;
#endif
#ifdef SL_USE_MPI
  struct {
fcs_near_fp_slint_t packed;
fcs_near_fp_slint_t db_packed;
fcs_near_fp_slint_t ip_packed;
  } mea;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef fcs_near_fp_MSEG_ROOT
int root;
#endif
#ifdef fcs_near_fp_MSEG_BORDER_UPDATE_REDUCTION
double border_update_count_reduction;
double border_update_weight_reduction;
#endif
#ifdef fcs_near_fp_MSEG_FORWARD_ONLY
fcs_near_fp_slint_t forward_only;
#endif
#ifdef fcs_near_fp_MSEG_INFO
fcs_near_fp_slint_t info_rounds;
fcs_near_fp_slint_t *info_finish_rounds;
double info_finish_rounds_avg;
fcs_near_fp_slweight_t info_total_weights;
#endif
fcs_near_fp_slint_t binnings;
fcs_near_fp_slint_t finalize_mode;
  } mseg;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef fcs_near_fp_MSS_ROOT
int root;
#endif
  } mss;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
fcs_near_fp_slint_t sync;
  } msm;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
fcs_near_fp_slint_t sync;
fcs_near_fp_partcond_t *r_pc;
  } msp;
#endif
#ifdef SL_USE_MPI
  struct {
double i_t[3];
double p_t[3];
double b_t[3];
fcs_near_fp_slint_t sync;
fcs_near_fp_slint_t i_sync;
fcs_near_fp_slint_t p_sync;
fcs_near_fp_slint_t b_sync;
fcs_near_fp_slint_t back_packed;
  } mssp;
#endif
} fcs_near_fp_sl_context_t;






/* sl_macro fcs_near_fp_elem_set_size fcs_near_fp_elem_set_max_size fcs_near_fp_elem_set_keys fcs_near_fp_elem_set_indices */
#define fcs_near_fp_elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define fcs_near_fp_elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define fcs_near_fp_elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define fcs_near_fp_elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))
/* DATAX_TEMPLATE_BEGIN */
#define fcs_near_fp_elem_set_data0(_e_, _b_)     ((_e_)->data0 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data1(_e_, _b_)     ((_e_)->data1 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data2(_e_, _b_)     ((_e_)->data2 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data3(_e_, _b_)     ((_e_)->data3 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data4(_e_, _b_)     ((_e_)->data4 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data5(_e_, _b_)     ((_e_)->data5 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data6(_e_, _b_)     ((_e_)->data6 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data7(_e_, _b_)     ((_e_)->data7 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data8(_e_, _b_)     ((_e_)->data8 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data9(_e_, _b_)     ((_e_)->data9 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data10(_e_, _b_)     ((_e_)->data10 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data11(_e_, _b_)     ((_e_)->data11 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data12(_e_, _b_)     ((_e_)->data12 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data13(_e_, _b_)     ((_e_)->data13 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data14(_e_, _b_)     ((_e_)->data14 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data15(_e_, _b_)     ((_e_)->data15 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data16(_e_, _b_)     ((_e_)->data16 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data17(_e_, _b_)     ((_e_)->data17 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data18(_e_, _b_)     ((_e_)->data18 = (_b_))  /* sl_macro */
#define fcs_near_fp_elem_set_data19(_e_, _b_)     ((_e_)->data19 = (_b_))  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro fcs_near_fp_elem_get_size fcs_near_fp_elem_get_max_size fcs_near_fp_elem_get_keys fcs_near_fp_elem_get_indices */
#define fcs_near_fp_elem_get_size(_e_)           (_e_)->size
#define fcs_near_fp_elem_get_max_size(_e_)       (_e_)->max_size
#define fcs_near_fp_elem_get_keys(_e_)           (_e_)->keys
#define fcs_near_fp_elem_get_indices(_e_)        (_e_)->indices
/* DATAX_TEMPLATE_BEGIN */
#define fcs_near_fp_elem_get_data0(_e_)          (_e_)->data0  /* sl_macro */
#define fcs_near_fp_elem_get_data1(_e_)          (_e_)->data1  /* sl_macro */
#define fcs_near_fp_elem_get_data2(_e_)          (_e_)->data2  /* sl_macro */
#define fcs_near_fp_elem_get_data3(_e_)          (_e_)->data3  /* sl_macro */
#define fcs_near_fp_elem_get_data4(_e_)          (_e_)->data4  /* sl_macro */
#define fcs_near_fp_elem_get_data5(_e_)          (_e_)->data5  /* sl_macro */
#define fcs_near_fp_elem_get_data6(_e_)          (_e_)->data6  /* sl_macro */
#define fcs_near_fp_elem_get_data7(_e_)          (_e_)->data7  /* sl_macro */
#define fcs_near_fp_elem_get_data8(_e_)          (_e_)->data8  /* sl_macro */
#define fcs_near_fp_elem_get_data9(_e_)          (_e_)->data9  /* sl_macro */
#define fcs_near_fp_elem_get_data10(_e_)          (_e_)->data10  /* sl_macro */
#define fcs_near_fp_elem_get_data11(_e_)          (_e_)->data11  /* sl_macro */
#define fcs_near_fp_elem_get_data12(_e_)          (_e_)->data12  /* sl_macro */
#define fcs_near_fp_elem_get_data13(_e_)          (_e_)->data13  /* sl_macro */
#define fcs_near_fp_elem_get_data14(_e_)          (_e_)->data14  /* sl_macro */
#define fcs_near_fp_elem_get_data15(_e_)          (_e_)->data15  /* sl_macro */
#define fcs_near_fp_elem_get_data16(_e_)          (_e_)->data16  /* sl_macro */
#define fcs_near_fp_elem_get_data17(_e_)          (_e_)->data17  /* sl_macro */
#define fcs_near_fp_elem_get_data18(_e_)          (_e_)->data18  /* sl_macro */
#define fcs_near_fp_elem_get_data19(_e_)          (_e_)->data19  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro fcs_near_fp_elem_set_block fcs_near_fp_elem_set_block_size fcs_near_fp_elem_get_block fcs_near_fp_elem_get_block_size */
#define fcs_near_fp_elem_set_block(_e_, _b_)       ((_e_)->keys = (_b_), (_e_)->max_size = -1)
#define fcs_near_fp_elem_set_block_size(_e_, _s_)  ((_e_)->size = (_s_))
#define fcs_near_fp_elem_get_block(_e_)            ((void *) (((_e_)->max_size < 0)?(_e_)->keys:NULL))
#define fcs_near_fp_elem_get_block_size(_e_)       (_e_)->size

/* sl_macro fcs_near_fp_pelem_set_size fcs_near_fp_pelem_set_max_size fcs_near_fp_pelem_set_elements */
#define fcs_near_fp_pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define fcs_near_fp_pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define fcs_near_fp_pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))

/* sl_macro fcs_near_fp_pelem_get_size fcs_near_fp_pelem_get_max_size fcs_near_fp_pelem_get_elements */
#define fcs_near_fp_pelem_get_size(_e_)           (_e_)->size
#define fcs_near_fp_pelem_get_max_size(_e_)       (_e_)->max_size
#define fcs_near_fp_pelem_get_elements(_e_)       (_e_)->elements

/* sl_macro fcs_near_fp_SL_DEFCON */
#define fcs_near_fp_SL_DEFCON(_v_)  (fcs_near_fp_sl_default_context._v_)






/* src/base/base.c */
extern fcs_near_fp_sl_context_t fcs_near_fp_sl_default_context;
extern const int fcs_near_fp_default_sl_dummy_rank;
#ifdef fcs_near_fp_SL_USE_RTI
extern const fcs_near_fp_rti_t fcs_near_fp_default_rti;
#endif
extern const fcs_near_fp_slint_t fcs_near_fp_default_sr_ip_threshold;
extern const fcs_near_fp_slint_t fcs_near_fp_default_sr_db_threshold;
extern const fcs_near_fp_slint_t fcs_near_fp_default_sr_ma_threshold;
extern const fcs_near_fp_slint_t fcs_near_fp_default_sri_threshold;

/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
extern const MPI_Datatype fcs_near_fp_default_mpi_int_datatype;
extern const MPI_Datatype fcs_near_fp_default_mpi_key_datatype;
extern const MPI_Datatype fcs_near_fp_default_mpi_pkey_datatype;
extern const MPI_Datatype fcs_near_fp_default_mpi_pwkey_datatype;
extern const MPI_Datatype fcs_near_fp_default_mpi_index_datatype;
extern const MPI_Datatype fcs_near_fp_default_mpi_weight_datatype;
extern const MPI_Datatype fcs_near_fp_default_mpi_data_datatype[];
extern const int fcs_near_fp_default_mpi_rank;
#endif
extern const void *fcs_near_fp_default_me_sendrecv_replace_mem;
extern const fcs_near_fp_slint_t fcs_near_fp_default_me_sendrecv_replace_memsize;
extern const fcs_near_fp_slint_t fcs_near_fp_default_me_sendrecv_replace_mpi_maxsize;
extern const double fcs_near_fp_default_meas_t[];
extern const fcs_near_fp_slint_t fcs_near_fp_default_meas_max_nprocs;
extern const fcs_near_fp_slint_t fcs_near_fp_default_meas_packed;
extern const fcs_near_fp_slint_t fcs_near_fp_default_meas_minalloc;
extern const double fcs_near_fp_default_meas_overalloc;
extern const fcs_near_fp_slint_t fcs_near_fp_default_mea_packed;
extern const fcs_near_fp_slint_t fcs_near_fp_default_mea_db_packed;
extern const fcs_near_fp_slint_t fcs_near_fp_default_mea_ip_packed;
#ifdef fcs_near_fp_MSEG_ROOT
extern const int fcs_near_fp_default_mseg_root;
#endif
#ifdef fcs_near_fp_MSEG_BORDER_UPDATE_REDUCTION
extern const double fcs_near_fp_default_mseg_border_update_count_reduction;
extern const double fcs_near_fp_default_mseg_border_update_weight_reduction;
#endif
#ifdef fcs_near_fp_MSEG_FORWARD_ONLY
extern const fcs_near_fp_slint_t fcs_near_fp_default_mseg_forward_only;
#endif
#ifdef fcs_near_fp_MSEG_INFO
extern const fcs_near_fp_slint_t fcs_near_fp_default_mseg_info_rounds;
extern const fcs_near_fp_slint_t *fcs_near_fp_default_mseg_info_finish_rounds;
extern const double fcs_near_fp_default_mseg_info_finish_rounds_avg;
extern const fcs_near_fp_slweight_t fcs_near_fp_default_mseg_info_total_weights;
#endif
extern const fcs_near_fp_slint_t fcs_near_fp_default_mseg_binnings;
extern const fcs_near_fp_slint_t fcs_near_fp_default_mseg_finalize_mode;
#ifdef fcs_near_fp_MSS_ROOT
extern const int fcs_near_fp_default_mss_root;
#endif
extern const double fcs_near_fp_default_msm_t[];
extern const fcs_near_fp_slint_t fcs_near_fp_default_msm_sync;
extern const double fcs_near_fp_default_msp_t[];
extern const fcs_near_fp_slint_t fcs_near_fp_default_msp_sync;
extern const fcs_near_fp_partcond_t *fcs_near_fp_default_msp_r_pc;
extern const double fcs_near_fp_default_mssp_i_t[];
extern const double fcs_near_fp_default_mssp_p_t[];
extern const double fcs_near_fp_default_mssp_b_t[];
extern const fcs_near_fp_slint_t fcs_near_fp_default_mssp_sync;
extern const fcs_near_fp_slint_t fcs_near_fp_default_mssp_i_sync;
extern const fcs_near_fp_slint_t fcs_near_fp_default_mssp_p_sync;
extern const fcs_near_fp_slint_t fcs_near_fp_default_mssp_b_sync;
extern const fcs_near_fp_slint_t fcs_near_fp_default_mssp_back_packed;






/* src/base/base.c */
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_create)(fcs_near_fp_local_bins_t *lb, fcs_near_fp_slint_t max_nbins, fcs_near_fp_slint_t max_nbinnings, fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slint_t docounts, fcs_near_fp_slint_t doweights, fcs_near_fp_binning_t *bm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_destroy)(fcs_near_fp_local_bins_t *lb);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_pre)(fcs_near_fp_local_bins_t *lb);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_exec_reset)(fcs_near_fp_local_bins_t *lb, fcs_near_fp_slint_t do_bins, fcs_near_fp_slint_t do_prefixes);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_exec)(fcs_near_fp_local_bins_t *lb, fcs_near_fp_slint_t b, fcs_near_fp_slint_t do_bins, fcs_near_fp_slint_t do_prefixes);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_refine)(fcs_near_fp_local_bins_t *lb, fcs_near_fp_slint_t b, fcs_near_fp_slint_t k, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_hit)(fcs_near_fp_local_bins_t *lb, fcs_near_fp_slint_t b, fcs_near_fp_slint_t k, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_finalize)(fcs_near_fp_local_bins_t *lb, fcs_near_fp_slint_t b, fcs_near_fp_slint_t dc, fcs_near_fp_slweight_t dw, fcs_near_fp_slint_t lc_min, fcs_near_fp_slint_t lc_max, fcs_near_fp_slcount_t *lcs, fcs_near_fp_slweight_t *lws, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_post)(fcs_near_fp_local_bins_t *lb);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_radix_create)(fcs_near_fp_binning_t *bm, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, fcs_near_fp_slint_t sorted);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_radix_destroy)(fcs_near_fp_binning_t *bm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_radix_pre)(fcs_near_fp_binning_t *bm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_radix_exec)(fcs_near_fp_binning_t *bm, fcs_near_fp_bin_t *bin, fcs_near_fp_slcount_t *counts, fcs_near_fp_slweight_t *weights);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_radix_refine)(fcs_near_fp_binning_t *bm, fcs_near_fp_bin_t *bin, fcs_near_fp_slint_t k, fcs_near_fp_slcount_t *counts, fcs_near_fp_slweight_t *weights, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s, fcs_near_fp_bin_t *new_bin);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_radix_hit)(fcs_near_fp_binning_t *bm, fcs_near_fp_bin_t *bin, fcs_near_fp_slint_t k, fcs_near_fp_slcount_t *counts, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_radix_finalize)(fcs_near_fp_binning_t *bm, fcs_near_fp_bin_t *bin, fcs_near_fp_slint_t dc, fcs_near_fp_slweight_t dw, fcs_near_fp_slint_t lc_min, fcs_near_fp_slint_t lc_max, fcs_near_fp_slcount_t *lcs, fcs_near_fp_slweight_t *lws, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_binning_radix_post)(fcs_near_fp_binning_t *bm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_alloc)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, slcint_t components);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_free)(fcs_near_fp_elements_t *s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_realloc)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, slcint_t components);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_alloca)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, slcint_t components);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_freea)(fcs_near_fp_elements_t *s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_alloc_from_blocks)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nblocks, void **blocks, fcs_near_fp_slint_t *blocksizes, fcs_near_fp_slint_t alignment, fcs_near_fp_slint_t nmax, slcint_t components);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_alloc_from_block)(fcs_near_fp_elements_t *s, void *block, fcs_near_fp_slint_t blocksize, fcs_near_fp_slint_t alignment, fcs_near_fp_slint_t nmax, slcint_t components);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_alloc_block)(fcs_near_fp_elements_t *s, void **block, fcs_near_fp_slint_t *blocksize, fcs_near_fp_slint_t alignment, fcs_near_fp_slint_t maxblocksize);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_copy)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_copy_at)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t sat, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t dat);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_ncopy)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t n);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_nmove)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t n);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_printf)(fcs_near_fp_elements_t *s, const char *prefix);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_extract)(fcs_near_fp_elements_t *src, fcs_near_fp_slint_t nelements, fcs_near_fp_elements_t *dst0, fcs_near_fp_elements_t *dst1);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_touch)(fcs_near_fp_elements_t *s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_digest_sum)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(fcs_near_fp_elements_crc32)(fcs_near_fp_elements_t *s, fcs_near_fp_slint nelements, fcs_near_fp_slint_t keys, fcs_near_fp_slint_t data);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_digest_hash)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, slcint_t components, void *hash);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_random_exchange)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t rounds, fcs_near_fp_elements_t *xs);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_keys_init_seed)(unsigned long s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_keys_init)(fcs_near_fp_elements_t *s, fcs_near_fp_keys_init_type_t t, fcs_near_fp_keys_init_data_t d);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_keys_init_randomized)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nkeys, fcs_near_fp_keys_init_type_t t, fcs_near_fp_keys_init_data_t d);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_keys_init_from_file)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t data, char *filename, fcs_near_fp_slint_t from, fcs_near_fp_slint_t to, fcs_near_fp_slint_t const_bytes_per_line);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_keys_save_to_file)(fcs_near_fp_elements_t *s, char *filename);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_validate_order)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t n);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_validate_order_bmask)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t n, fcs_near_fp_slkey_pure_t bmask);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_validate_order_weight)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t n, fcs_near_fp_slkey_pure_t weight);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_keys_stats)(fcs_near_fp_elements_t *s, fcs_near_fp_slkey_pure_t *stats);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_keys_stats_print)(fcs_near_fp_elements_t *s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_print_keys)(fcs_near_fp_elements_t *s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_print_all)(fcs_near_fp_elements_t *s);
fcs_near_fp_slweight_t SL_PROTO(fcs_near_fp_elements_get_weight)(fcs_near_fp_elements_t *s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_get_minmax_keys)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slkey_pure_t *minmaxkeys);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_alloc_packed)(fcs_near_fp_packed_elements_t *s, fcs_near_fp_slint_t nelements);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_free_packed)(fcs_near_fp_packed_elements_t *s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_alloc_packed_from_block)(fcs_near_fp_packed_elements_t *s, void *block, fcs_near_fp_slint_t blocksize, fcs_near_fp_slint_t alignment, fcs_near_fp_slint_t nmax);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_pack_indexed)(fcs_near_fp_elements_t *s, fcs_near_fp_packed_elements_t *d, fcs_near_fp_slindex_t *rindx, fcs_near_fp_slindex_t *windx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_pack)(fcs_near_fp_elements_t *s, fcs_near_fp_packed_elements_t *d);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_pack_at)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t sat, fcs_near_fp_packed_elements_t *d, fcs_near_fp_slint_t dat);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_unpack_indexed)(fcs_near_fp_packed_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slindex_t *rindx, fcs_near_fp_slindex_t *windx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_unpack)(fcs_near_fp_packed_elements_t *s, fcs_near_fp_elements_t *d);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_unpack_at)(fcs_near_fp_packed_elements_t *s, fcs_near_fp_slint_t sat, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t dat);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elements_unpack_keys)(fcs_near_fp_packed_elements_t *s, fcs_near_fp_slkey_t *k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_merge2_basic_auto_01_x)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_merge2_basic_01_x)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx, fcs_near_fp_m2x_func _x0_1, fcs_near_fp_m2x_func _0x_1);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_merge2_basic_01_X)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx, fcs_near_fp_elements_t *t, fcs_near_fp_m2X_func _X0_1, fcs_near_fp_m2X_func _0X_1);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_merge2_simplify_s1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx, fcs_near_fp_slint s1elements);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_merge2_memory_adaptive)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_compo_hula)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *xs);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_sseq_x0_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_sseq_0x_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_sseq_01_x)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_sseq_01)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_sbin_x0_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_sbin_0x_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_sbin_01_x)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_sbin_01)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_shyb_x0_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_shyb_0x_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_shyb_01_x)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_shyb_01)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_straight_x0_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_straight_0x_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_straight_01_x)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_straight_x_0_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_straight_X0_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_straight_0X_1)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_straight_01_X)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_basic_straight_X0_1u)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge2_compo_tridgell)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mergep_2way_ip_int)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t p, int *displs, fcs_near_fp_merge2x_f m2x);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mergep_2way_ip_int_rec)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t p, int *displs, fcs_near_fp_merge2x_f m2x);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mergep_heap_int)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t p, int *displs, int *counts);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mergep_heap_int_idx)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t p, int *displs, int *counts);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mergep_heap_idx)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t p, fcs_near_fp_slindex_t *displs, fcs_near_fp_slindex_t *counts);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mergep_heap_unpack_idx)(fcs_near_fp_packed_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t p, fcs_near_fp_slindex_t *displs, fcs_near_fp_slindex_t *counts);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mergep_heap_unpack_idxonly)(fcs_near_fp_packed_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t p, fcs_near_fp_slindex_t *displs, fcs_near_fp_slindex_t *counts);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_permute_generic_db)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_permute_generic_t *pg, void *pg_data);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_permute_generic_ip)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *x, fcs_near_fp_permute_generic_t *pg, void *pg_data);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_sequential_lt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_sequential_le)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_sequential_gt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_sequential_ge)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_sequential_lt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_sequential_le)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_sequential_gt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_sequential_ge)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_binary_lt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_binary_le)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_binary_gt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_binary_ge)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_binary_lt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_binary_le)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_binary_gt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_binary_ge)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sl_search_binary_lt_bmask)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k, fcs_near_fp_slpkey_t bmask);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sl_search_binary_le_bmask)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k, fcs_near_fp_slpkey_t bmask);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sl_search_binary_sign_switch)(fcs_near_fp_elements_t *s);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_hybrid_lt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k, fcs_near_fp_slint t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_hybrid_le)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k, fcs_near_fp_slint t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_hybrid_gt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k, fcs_near_fp_slint t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_hybrid_ge)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t k, fcs_near_fp_slint t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_hybrid_lt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k, fcs_near_fp_slint t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_hybrid_le)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k, fcs_near_fp_slint t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_hybrid_gt)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k, fcs_near_fp_slint t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sl_search_p_hybrid_ge)(fcs_near_fp_elements_t *s, fcs_near_fp_slpkey_t *k, fcs_near_fp_slint t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_ilog2c)(fcs_near_fp_slint x);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_ilog2f)(fcs_near_fp_slint x);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_print_bits)(fcs_near_fp_slint v);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_pivot_random)(fcs_near_fp_elements_t *s);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_counts2displs)(fcs_near_fp_slint_t n, int *counts, int *displs);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_displs2counts)(fcs_near_fp_slint_t n, int *displs, int *counts, fcs_near_fp_slint_t total_counts);
void SL_PROTO(fcs_near_fp_get_displcounts_extent)(fcs_near_fp_slint_t n, int *displs, int *counts, fcs_near_fp_slint_t *lb, fcs_near_fp_slint_t *extent);
void SL_PROTO(fcs_near_fp_elem_set_data)(fcs_near_fp_elements_t *e, ...);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elem_get_max_byte)();
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elem_reverse)(fcs_near_fp_elements_t *e, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elem_nxchange_at)(fcs_near_fp_elements_t *e0, fcs_near_fp_slint_t at0, fcs_near_fp_elements_t *e1, fcs_near_fp_slint_t at1, fcs_near_fp_slint_t n, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elem_nxchange)(fcs_near_fp_elements_t *e0, fcs_near_fp_elements_t *e1, fcs_near_fp_slint_t n, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elem_nxchange_ro0)(fcs_near_fp_elements_t *e0, fcs_near_fp_elements_t *e1, fcs_near_fp_slint_t n, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elem_rotate)(fcs_near_fp_elements_t *e, fcs_near_fp_slint_t m, fcs_near_fp_slint_t n, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elem_rotate_ro0)(fcs_near_fp_elements_t *e, fcs_near_fp_slint_t m, fcs_near_fp_slint_t n, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_elem_rotate_ro1)(fcs_near_fp_elements_t *e, fcs_near_fp_slint_t m, fcs_near_fp_slint_t n, fcs_near_fp_elements_t *t);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_counting_use_displs)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t ndispls, fcs_near_fp_slint_t *displs);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_counting_use_counts)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t ncounts, fcs_near_fp_slint_t *counts);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_counting_get_counts)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t ncounts, fcs_near_fp_slint_t *counts);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_counting)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_slint_t ncounts);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sort_heap)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *xs);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_insert_bmask_kernel)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slkey_pure_t bmask);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_insert)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_permute_forward)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t *perm, fcs_near_fp_slint_t offset, fcs_near_fp_slint_t mask_bit);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_permute_backward)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t *perm, fcs_near_fp_slint_t offset, fcs_near_fp_slint_t mask_bit);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sort_quick)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *xs);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_radix_ip)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_radix_db)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_radix_ma)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_radix)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_radix_1bit_kernel)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sort_radix_1bit)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_sort_radix_iter)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t presorted, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sn_hypercube_lh)(fcs_near_fp_slint size, fcs_near_fp_slint rank, fcs_near_fp_slint stage, void *snp, fcs_near_fp_slint *up);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sn_hypercube_hl)(fcs_near_fp_slint size, fcs_near_fp_slint rank, fcs_near_fp_slint stage, void *snp, fcs_near_fp_slint *up);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sn_odd_even_trans)(fcs_near_fp_slint size, fcs_near_fp_slint rank, fcs_near_fp_slint stage, void *snp, fcs_near_fp_slint *up);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sn_odd)(fcs_near_fp_slint size, fcs_near_fp_slint rank, fcs_near_fp_slint stage, void *snp, fcs_near_fp_slint *up);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sn_even)(fcs_near_fp_slint size, fcs_near_fp_slint rank, fcs_near_fp_slint stage, void *snp, fcs_near_fp_slint *up);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sn_batcher)(fcs_near_fp_slint size, fcs_near_fp_slint rank, fcs_near_fp_slint stage, void *snp, fcs_near_fp_slint *up);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sn_bitonic)(fcs_near_fp_slint size, fcs_near_fp_slint rank, fcs_near_fp_slint stage, void *snp, fcs_near_fp_slint *up);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_sn_connected)(fcs_near_fp_slint size, fcs_near_fp_slint rank, fcs_near_fp_slint stage, void *snp, fcs_near_fp_slint *up);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_split_generic_db)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_split_generic_t *sg, void *sg_data, fcs_near_fp_slint_t n);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_split_generic_ip)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_split_generic_t *sg, void *sg_data, fcs_near_fp_slint_t n);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_split_generic_count_db)(fcs_near_fp_elements_t *s, fcs_near_fp_split_generic_t *sg, void *sg_data, int *counts, fcs_near_fp_slint_t n);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_split_generic_count_ip)(fcs_near_fp_elements_t *s, fcs_near_fp_split_generic_t *sg, void *sg_data, int *counts, fcs_near_fp_slint_t n);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_split_generic_rearrange_db)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_split_generic_t *sg, void *sg_data, int *counts, fcs_near_fp_slint_t n);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_split_generic_rearrange_ip)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *d, fcs_near_fp_split_generic_t *sg, void *sg_data, int *counts, int *displs, fcs_near_fp_slint_t n);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_splitter_reset)(fcs_near_fp_splitter_t *sp);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_splitx_radix)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t nclasses, fcs_near_fp_slint_t shl, fcs_near_fp_slint_t *counts);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_split2_lt_ge)(fcs_near_fp_elements_t *s, fcs_near_fp_slkey_pure_t *k, fcs_near_fp_elements_t *t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_split2_le_gt)(fcs_near_fp_elements_t *s, fcs_near_fp_slkey_pure_t *k, fcs_near_fp_elements_t *t);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_split3_lt_eq_gt)(fcs_near_fp_elements_t *s, fcs_near_fp_slkey_pure_t *k, fcs_near_fp_elements_t *t, fcs_near_fp_slint *nlt, fcs_near_fp_slint *nle);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_split3_lt_eq_gt_old)(fcs_near_fp_elements_t *s, fcs_near_fp_slkey_pure_t *k, fcs_near_fp_elements_t *t, fcs_near_fp_slint *nlt, fcs_near_fp_slint *nle);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_split2_b)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slkey_pure_t bmask);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_splitk_k2c_af)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint k, fcs_near_fp_slint *c, fcs_near_fp_k2c_func k2c, void *k2c_data);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_splitk_k2c)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_slint k, fcs_near_fp_slint *c, fcs_near_fp_k2c_func k2c, void *k2c_data);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_splitk_k2c_count)(fcs_near_fp_elements_t *s, fcs_near_fp_slint k, fcs_near_fp_slint *c, fcs_near_fp_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI





/* src/base_mpi/base_mpi.c */
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_create)(fcs_near_fp_global_bins_t *gb, fcs_near_fp_slint_t max_nbins, fcs_near_fp_slint_t max_nbinnings, fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slint_t docounts, fcs_near_fp_slint_t doweights, fcs_near_fp_binning_t *bm, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_destroy)(fcs_near_fp_global_bins_t *gb, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_pre)(fcs_near_fp_global_bins_t *gb, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_exec_reset)(fcs_near_fp_global_bins_t *gb, fcs_near_fp_slint_t do_bins, fcs_near_fp_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_exec_local)(fcs_near_fp_global_bins_t *gb, fcs_near_fp_slint_t b, fcs_near_fp_slint_t do_bins, fcs_near_fp_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_exec_global)(fcs_near_fp_global_bins_t *gb, fcs_near_fp_slint_t do_bins, fcs_near_fp_slint_t do_prefixes, fcs_near_fp_slint_t root, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_refine)(fcs_near_fp_global_bins_t *gb, fcs_near_fp_slint_t b, fcs_near_fp_slint_t k, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_hit)(fcs_near_fp_global_bins_t *gb, fcs_near_fp_slint_t b, fcs_near_fp_slint_t k, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_finalize)(fcs_near_fp_global_bins_t *gb, fcs_near_fp_slint_t b, fcs_near_fp_slint_t dc, fcs_near_fp_slweight_t dw, fcs_near_fp_slint_t lc_min, fcs_near_fp_slint_t lc_max, fcs_near_fp_slcount_t *lcs, fcs_near_fp_slweight_t *lws, fcs_near_fp_splitter_t *sp, fcs_near_fp_slint_t s, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_binning_post)(fcs_near_fp_global_bins_t *gb, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_datatypes_init)();
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_datatypes_release)();
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_get_grid_properties)(fcs_near_fp_slint_t ndims, fcs_near_fp_slint_t *dims, fcs_near_fp_slint_t *pos, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_subgroups_create)(fcs_near_fp_slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_subgroups_delete)(fcs_near_fp_slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);
int SL_PROTO(fcs_near_fp_sl_MPI_Allreduce)(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank);
int SL_PROTO(fcs_near_fp_sl_MPI_Alltoall_int)(void *sendbuf, int sendcount, void *recvbuf, int recvcount, MPI_Comm comm, int size, int rank);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_keys_init_from_file)(fcs_near_fp_elements_t *s, char *filename, fcs_near_fp_slint from, fcs_near_fp_slint to, fcs_near_fp_slint const_bytes_per_line, fcs_near_fp_slint root, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint SL_PROTO(fcs_near_fp_mpi_elements_validate_order)(fcs_near_fp_elements_t *s, fcs_near_fp_slint n, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_linear_exchange_pure_keys)(fcs_near_fp_slkey_pure_t *in, fcs_near_fp_slkey_pure_t *out, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_check_order)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slint_t *orders, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_check_global_order)(fcs_near_fp_slkey_pure_t local_min, fcs_near_fp_slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_digest_sum)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(fcs_near_fp_mpi_elements_crc32)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t n, fcs_near_fp_slint_t keys, fcs_near_fp_slint_t data, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_digest_hash)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_get_counts)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t *clocal, fcs_near_fp_slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
fcs_near_fp_slweight_t SL_PROTO(fcs_near_fp_mpi_elements_get_weights)(fcs_near_fp_elements_t *s, fcs_near_fp_slweight_t *wlocal, fcs_near_fp_slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_get_counts_and_weights)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slint_t *counts, fcs_near_fp_slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_sendrecv_replace)(fcs_near_fp_elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_tproc_create_tproc)(fcs_near_fp_tproc_t *tproc, fcs_near_fp_tproc_f *tfn, fcs_near_fp_tproc_reset_f *rfn, fcs_near_fp_tproc_exdef exdef);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_tproc_create_tproc_mod)(fcs_near_fp_tproc_t *tproc, fcs_near_fp_tproc_mod_f *tfn, fcs_near_fp_tproc_reset_f *rfn, fcs_near_fp_tproc_exdef exdef);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_tproc_create_tprocs)(fcs_near_fp_tproc_t *tproc, fcs_near_fp_tprocs_f *tfn, fcs_near_fp_tproc_reset_f *rfn, fcs_near_fp_tproc_exdef exdef);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_tproc_create_tprocs_mod)(fcs_near_fp_tproc_t *tproc, fcs_near_fp_tprocs_mod_f *tfn, fcs_near_fp_tproc_reset_f *rfn, fcs_near_fp_tproc_exdef exdef);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_tproc_free)(fcs_near_fp_tproc_t *tproc);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_tproc_set_proclists)(fcs_near_fp_tproc_t *tproc, fcs_near_fp_slint_t nsend_procs, int *send_procs, fcs_near_fp_slint_t nrecv_procs, int *recv_procs, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_tproc_verify)(fcs_near_fp_tproc_t tproc, void *data, fcs_near_fp_elements_t *s, int proc);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_alltoall_specific)(fcs_near_fp_elements_t *sin, fcs_near_fp_elements_t *sout, fcs_near_fp_elements_t *xs, fcs_near_fp_tproc_t tproc, void *data, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_alltoallv_db_packed)(fcs_near_fp_elements_t *sbuf, int *scounts, int *sdispls, fcs_near_fp_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_alltoallv_db)(fcs_near_fp_elements_t *sbuf, int *scounts, int *sdispls, fcs_near_fp_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_alltoallv_ip_packed)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_alltoallv_ip_double)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_alltoallv_ip_mpi)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_alltoallv_ip_dash)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_alltoallv_ip)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_alltoallv_proclists_db)(fcs_near_fp_elements_t *sbuf, int *scounts, int *sdispls, int nsendprocs, int *sendprocs, fcs_near_fp_elements_t *rbuf, int *rcounts, int *rdispls, int nrecvprocs, int *recvprocs, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, fcs_near_fp_slint_t structured);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_find_exact_equal)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t other_rank, fcs_near_fp_slint_t high_rank, fcs_near_fp_slint_t *ex_start, fcs_near_fp_slint_t *ex_size, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_find_exact)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t other_rank, fcs_near_fp_slint_t high_rank, fcs_near_fp_slint_t *dst_size, fcs_near_fp_slint_t *ex_start, fcs_near_fp_slint_t *ex_sizes, fcs_near_fp_slint_t *nx_move, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_merge2)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t other_rank, fcs_near_fp_slint_t high_rank, fcs_near_fp_slint_t *dst_size, fcs_near_fp_merge2x_f m2, fcs_near_fp_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_mergek_equal)(fcs_near_fp_elements_t *s, fcs_near_fp_sortnet_f sn, fcs_near_fp_sortnet_data_t snd, fcs_near_fp_merge2x_f m2x, fcs_near_fp_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_mergek_sorted)(fcs_near_fp_elements_t *s, fcs_near_fp_merge2x_f m2x, fcs_near_fp_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_mergek_sorted2)(fcs_near_fp_elements_t *s, fcs_near_fp_sortnet_f sn, fcs_near_fp_sortnet_data_t snd, fcs_near_fp_merge2x_f m2x, fcs_near_fp_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_mergek)(fcs_near_fp_elements_t *s, fcs_near_fp_sortnet_f sn, fcs_near_fp_sortnet_data_t snd, fcs_near_fp_merge2x_f m2x, fcs_near_fp_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_mergek_equal2)(fcs_near_fp_elements_t *s, fcs_near_fp_sortnet_f sn, fcs_near_fp_sortnet_data_t snd, fcs_near_fp_merge2x_f m2x, fcs_near_fp_elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_partition_exact_generic)(fcs_near_fp_elements_t *s, fcs_near_fp_partcond_t *pcond, fcs_near_fp_binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_partition_exact_radix)(fcs_near_fp_elements_t *s, fcs_near_fp_partcond_t *pcond, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, fcs_near_fp_slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_partition_exact_radix_ngroups)(fcs_near_fp_elements_t *s, fcs_near_fp_partcond_t *pcond, fcs_near_fp_slint_t ngroups, MPI_Comm *group_comms, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_partition_exact_radix_2groups)(fcs_near_fp_elements_t *s, fcs_near_fp_partcond_t *pcond, MPI_Comm group_comm, fcs_near_fp_elements_t *sx, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_partition_sample_regular)(fcs_near_fp_elements_t *s, fcs_near_fp_partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_rebalance)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_slint_t stable, fcs_near_fp_slint_t *dst_size, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_rebalance_alltoallv)(fcs_near_fp_elements_t *sbuf, int *scounts, int *sdispls, fcs_near_fp_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
void SL_PROTO(fcs_near_fp_mpi_partcond_set_even)(fcs_near_fp_partcond_t *pcond, fcs_near_fp_slint_t pcm, fcs_near_fp_slint_t ntotal, double nimba, double wtotal, double wimba, int size, int rank);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_init_partconds)(fcs_near_fp_slint_t npconds, fcs_near_fp_partcond_t *pconds, fcs_near_fp_slint_t nparts, fcs_near_fp_slint_t total_count, fcs_near_fp_slweight_t total_weight);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_init_partconds_intern)(fcs_near_fp_slint_t npconds, fcs_near_fp_partcond_intern_t *pci, fcs_near_fp_partcond_t *pc, fcs_near_fp_slint_t nparts, fcs_near_fp_slint_t total_count, fcs_near_fp_slweight_t total_weight);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_merge_partconds)(fcs_near_fp_partcond_t *pconds_in, fcs_near_fp_slint_t npconds_in, fcs_near_fp_partcond_t *pcond_out);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_gather_partconds_grouped)(fcs_near_fp_partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, fcs_near_fp_partcond_t *pconds_out, fcs_near_fp_slint_t *npconds_out, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_gather_partconds)(fcs_near_fp_partcond_t *pcond_in, fcs_near_fp_partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_allgather_partconds)(fcs_near_fp_partcond_t *pcond_in, fcs_near_fp_partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_bcast_partconds)(fcs_near_fp_slint_t npconds, fcs_near_fp_partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_post_check_partconds)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slint_t nparts, fcs_near_fp_partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_post_check_partconds_intern)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slint_t nparts, fcs_near_fp_partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_select_stats)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_select_exact_generic_bulk)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slint_t nparts, fcs_near_fp_partcond_t *pconds, fcs_near_fp_binning_t *bm, fcs_near_fp_splitter_t *sp, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_select_exact_generic_grouped)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, fcs_near_fp_binning_t *bm, fcs_near_fp_splitter_t *sp, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_select_exact_generic)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slint_t nparts, fcs_near_fp_partcond_t *pconds, fcs_near_fp_binning_t *bm, fcs_near_fp_splitter_t *sp, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_select_exact_radix)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_slint_t nparts, fcs_near_fp_partcond_t *pconds, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, fcs_near_fp_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_select_exact_radix_grouped)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nelements, fcs_near_fp_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, fcs_near_fp_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_select_sample_regular)(fcs_near_fp_elements_t *s, fcs_near_fp_slint_t nparts, fcs_near_fp_partcond_t *pconds, fcs_near_fp_slint_t nsamples, fcs_near_fp_splitter_t *sp, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_merge)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *xs, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_merge2)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *xs, fcs_near_fp_slint_t merge_type, fcs_near_fp_slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_merge_radix)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *xs, fcs_near_fp_slint_t merge_type, fcs_near_fp_slint_t sort_type, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_partition)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *xs, fcs_near_fp_slint_t part_type, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_partition_radix)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *xs, fcs_near_fp_slint_t part_type, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_partition_exact_radix)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_partcond_t *pcond, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_partition_exact_radix_ngroups)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_partcond_t *pcond, fcs_near_fp_slint_t ngroups, MPI_Comm *group_comms, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_partition_exact_radix_2groups)(fcs_near_fp_elements_t *s, fcs_near_fp_elements_t *sx, fcs_near_fp_partcond_t *pcond, MPI_Comm group_comm, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_insert_radix)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *xs, fcs_near_fp_slpkey_t *mmkeys, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_presorted_radix)(fcs_near_fp_elements_t *s0, fcs_near_fp_elements_t *s1, fcs_near_fp_elements_t *xs, fcs_near_fp_slint_t merge_type, fcs_near_fp_slint_t rhigh, fcs_near_fp_slint_t rlow, fcs_near_fp_slint_t rwidth, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_sort_back)(fcs_near_fp_elements_t *sin, fcs_near_fp_elements_t *sout, fcs_near_fp_elements_t *sx, fcs_near_fp_slpkey_t *lh, fcs_near_fp_slint_t ntotal, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, fcs_near_fp_slint_t ytotal, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_xcounts2ycounts_grouped)(int *xcounts, fcs_near_fp_slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
fcs_near_fp_slint_t SL_PROTO(fcs_near_fp_mpi_subxdispls2ycounts)(fcs_near_fp_slint_t nsubs, int *sub_xdispls, fcs_near_fp_slint_t *sub_sources, fcs_near_fp_slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#undef SL_PROTO
#endif /* __SL_NEAR_FP_H__ */
