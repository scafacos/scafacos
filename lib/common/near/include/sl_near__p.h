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


#ifndef __SL_NEAR__P_H__
#define __SL_NEAR__P_H__

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif /* SL_USE_MPI */

#define SL_PROTO(_f_)  _f_

#define near__p_sl_key_type_c        long long
#define near__p_sl_key_type_mpi      MPI_LONG_LONG
#define near__p_sl_key_size_mpi      1
#define near__p_sl_key_type_fmt      "lld"

#define near__p_sl_key_integer

#define near__p_SL_DATA0             /* positions */
#define near__p_sl_data0_type_c      fcs_float
#define near__p_sl_data0_size_c      3
#define near__p_sl_data0_type_mpi    FCS_MPI_FLOAT
#define near__p_sl_data0_size_mpi    3

#define near__p_SL_DATA1             /* charges */
#define near__p_sl_data1_type_c      fcs_float
#define near__p_sl_data1_size_c      1
#define near__p_sl_data1_type_mpi    FCS_MPI_FLOAT
#define near__p_sl_data1_size_mpi    1

#define near__p_SL_DATA2             /* indices */
#define near__p_sl_data2_type_c      long long
#define near__p_sl_data2_size_c      1
#define near__p_sl_data2_type_mpi    MPI_LONG_LONG
#define near__p_sl_data2_size_mpi    1

#undef near__p_SL_DATA3              /* field */
#define near__p_sl_data3_type_c      fcs_float
#define near__p_sl_data3_size_c      3
#define near__p_sl_data3_type_mpi    FCS_MPI_FLOAT
#define near__p_sl_data3_size_mpi    3

#define near__p_SL_DATA4             /* potentials */
#define near__p_sl_data4_type_c      fcs_float
#define near__p_sl_data4_size_c      1
#define near__p_sl_data4_type_mpi    FCS_MPI_FLOAT
#define near__p_sl_data4_size_mpi    1




#if defined(MSEG_ROOT) && !defined(near__p_MSEG_ROOT)
# define near__p_MSEG_ROOT  MSEG_ROOT
#endif

#if defined(MSEG_BORDER_UPDATE_REDUCTION) && !defined(near__p_MSEG_BORDER_UPDATE_REDUCTION)
# define near__p_MSEG_BORDER_UPDATE_REDUCTION  MSEG_BORDER_UPDATE_REDUCTION
#endif

#if defined(MSEG_DISABLE_BEST_CHOICE) && !defined(near__p_MSEG_DISABLE_BEST_CHOICE)
# define near__p_MSEG_DISABLE_BEST_CHOICE  MSEG_DISABLE_BEST_CHOICE
#endif

#if defined(MSEG_DISABLE_MINMAX) && !defined(near__p_MSEG_DISABLE_MINMAX)
# define near__p_MSEG_DISABLE_MINMAX  MSEG_DISABLE_MINMAX
#endif

#if defined(MSEG_ENABLE_OPTIMZED_LOWHIGH) && !defined(near__p_MSEG_ENABLE_OPTIMZED_LOWHIGH)
# define near__p_MSEG_ENABLE_OPTIMZED_LOWHIGH  MSEG_ENABLE_OPTIMZED_LOWHIGH
#endif

#if defined(MSEG_FORWARD_ONLY) && !defined(near__p_MSEG_FORWARD_ONLY)
# define near__p_MSEG_FORWARD_ONLY  MSEG_FORWARD_ONLY
#endif

#if defined(MSEG_INFO) && !defined(near__p_MSEG_INFO)
# define near__p_MSEG_INFO  MSEG_INFO
#endif

#if defined(MSEG_TRACE_IF) && !defined(near__p_MSEG_TRACE_IF)
# define near__p_MSEG_TRACE_IF  MSEG_TRACE_IF
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


#ifndef near__p_SL_INDEX
# undef near__p_SL_PACKED_INDEX
#endif


/* if no special datatype for (sl default) integer ... */
#ifndef near__p_sl_int_type_c
  /* ... use a default one */
# define near__p_sl_int_type_c               long      /* sl_macro */
# undef near__p_sl_int_type_mpi
# define near__p_sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef near__p_sl_int_size_mpi
# define near__p_sl_int_size_mpi             1         /* sl_macro */
# undef near__p_sl_int_type_fmt
# define near__p_sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(near__p_sl_int_type_mpi) || !defined(near__p_sl_int_size_mpi)
#   error "near__p_sl_int_type_mpi and/or near__p_sl_int_size_mpi missing"
#  endif
# endif
# ifndef near__p_sl_int_type_fmt
#  error "near__p_sl_int_type_fmt macro is missing, using d as default"
#  define near__p_sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for (intern) weight ... */
#ifndef near__p_sl_weight_type_c
 /* ... use (sl default) integer */
# define near__p_sl_weight_type_c             near__p_sl_int_type_c    /* sl_macro */
# undef near__p_sl_weight_type_mpi
# define near__p_sl_weight_type_mpi           near__p_sl_int_type_mpi  /* sl_macro */
# undef near__p_sl_weight_size_mpi
# define near__p_sl_weight_size_mpi           near__p_sl_int_size_mpi  /* sl_macro */
# undef near__p_sl_weight_type_fmt
# define near__p_sl_weight_type_fmt           near__p_sl_int_type_fmt  /* sl_macro */
# undef near__p_sl_weight_intequiv
# define near__p_sl_weight_intequiv                            /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(near__p_sl_weight_type_mpi) || !defined(near__p_sl_weight_size_mpi)
#   error "near__p_sl_weight_type_mpi and/or near__p_sl_weight_size_mpi missing"
#  endif
# endif
# ifndef near__p_sl_weight_type_fmt
#  error "near__p_sl_weight_type_fmt macro is missing, using f as default"
#  define near__p_sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef near__p_sl_index_type_c
 /* ... use the primary integer type */
# define near__p_sl_index_type_c             near__p_sl_int_type_c
# undef near__p_sl_index_type_mpi
# define near__p_sl_index_type_mpi           near__p_sl_int_type_mpi
# undef near__p_sl_index_size_mpi
# define near__p_sl_index_size_mpi           near__p_sl_int_size_mpi
# undef near__p_sl_index_type_fmt
# define near__p_sl_index_type_fmt           near__p_sl_int_type_fmt
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(near__p_sl_index_type_mpi) || !defined(near__p_sl_index_size_mpi)
#   error "near__p_sl_index_type_mpi and/or near__p_sl_index_size_mpi missing"
#  endif
# endif
# ifndef near__p_sl_index_type_fmt
#  error "near__p_sl_index_type_fmt macro is missing, using d as default"
#  define near__p_sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef near__p_sl_key_pure_type_c
# define near__p_sl_key_pure_type_c          near__p_sl_key_type_c  /* sl_macro */
#endif
#ifndef near__p_sl_key_pure_type_mpi
# define near__p_sl_key_pure_type_mpi        near__p_sl_key_type_mpi  /* sl_macro */
#endif
#ifndef near__p_sl_key_pure_size_mpi
# define near__p_sl_key_pure_size_mpi        near__p_sl_key_size_mpi  /* sl_macro */
#endif
#ifndef near__p_sl_key_pure_type_fmt
# ifdef near__p_sl_key_type_fmt
#  define near__p_sl_key_pure_type_fmt       near__p_sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef near__p_sl_key_purify
 /* key val -> key val */
 #define near__p_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef near__p_sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define near__p_sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef near__p_sl_key_set_pure
 /* key component pointer and key val */
 #define near__p_sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef near__p_sl_key_pure_cmp_eq
 #define near__p_sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef near__p_sl_key_pure_cmp_ne
 #define near__p_sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef near__p_sl_key_pure_cmp_lt
 #define near__p_sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef near__p_sl_key_pure_cmp_le
 #define near__p_sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef near__p_sl_key_pure_cmp_gt
 #define near__p_sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef near__p_sl_key_pure_cmp_ge
 #define near__p_sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef near__p_sl_key_cmp_eq
 #define near__p_sl_key_cmp_eq(k0, k1)       (near__p_sl_key_pure_cmp_eq(near__p_sl_key_purify(k0), near__p_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near__p_sl_key_cmp_ne
 #define near__p_sl_key_cmp_ne(k0, k1)       (near__p_sl_key_pure_cmp_ne(near__p_sl_key_purify(k0), near__p_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near__p_sl_key_cmp_lt
 #define near__p_sl_key_cmp_lt(k0, k1)       (near__p_sl_key_pure_cmp_lt(near__p_sl_key_purify(k0), near__p_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near__p_sl_key_cmp_le
 #define near__p_sl_key_cmp_le(k0, k1)       (near__p_sl_key_pure_cmp_le(near__p_sl_key_purify(k0), near__p_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near__p_sl_key_cmp_gt
 #define near__p_sl_key_cmp_gt(k0, k1)       (near__p_sl_key_pure_cmp_gt(near__p_sl_key_purify(k0), near__p_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near__p_sl_key_cmp_ge
 #define near__p_sl_key_cmp_ge(k0, k1)       (near__p_sl_key_pure_cmp_ge(near__p_sl_key_purify(k0), near__p_sl_key_purify(k1)))  /* sl_macro */
#endif


/* default random key */
#ifdef near__p_sl_key_integer
# if !defined(near__p_sl_key_val_srand) || !defined(near__p_sl_key_val_rand) || !defined(near__p_sl_key_val_rand_minmax)
#  undef near__p_sl_key_val_srand
#  undef near__p_sl_key_val_rand
#  undef near__p_sl_key_val_rand_minmax
#  define near__p_sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define near__p_sl_key_val_rand()                     ((near__p_sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define near__p_sl_key_val_rand_minmax(_min_, _max_)  ((near__p_sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef near__p_SL_DATA0_IGNORE
# undef near__p_SL_DATA0
#endif
#ifdef near__p_SL_DATA1_IGNORE
# undef near__p_SL_DATA1
#endif
#ifdef near__p_SL_DATA2_IGNORE
# undef near__p_SL_DATA2
#endif
#ifdef near__p_SL_DATA3_IGNORE
# undef near__p_SL_DATA3
#endif
#ifdef near__p_SL_DATA4_IGNORE
# undef near__p_SL_DATA4
#endif
#ifdef near__p_SL_DATA5_IGNORE
# undef near__p_SL_DATA5
#endif
#ifdef near__p_SL_DATA6_IGNORE
# undef near__p_SL_DATA6
#endif
#ifdef near__p_SL_DATA7_IGNORE
# undef near__p_SL_DATA7
#endif
#ifdef near__p_SL_DATA8_IGNORE
# undef near__p_SL_DATA8
#endif
#ifdef near__p_SL_DATA9_IGNORE
# undef near__p_SL_DATA9
#endif
#ifdef near__p_SL_DATA10_IGNORE
# undef near__p_SL_DATA10
#endif
#ifdef near__p_SL_DATA11_IGNORE
# undef near__p_SL_DATA11
#endif
#ifdef near__p_SL_DATA12_IGNORE
# undef near__p_SL_DATA12
#endif
#ifdef near__p_SL_DATA13_IGNORE
# undef near__p_SL_DATA13
#endif
#ifdef near__p_SL_DATA14_IGNORE
# undef near__p_SL_DATA14
#endif
#ifdef near__p_SL_DATA15_IGNORE
# undef near__p_SL_DATA15
#endif
#ifdef near__p_SL_DATA16_IGNORE
# undef near__p_SL_DATA16
#endif
#ifdef near__p_SL_DATA17_IGNORE
# undef near__p_SL_DATA17
#endif
#ifdef near__p_SL_DATA18_IGNORE
# undef near__p_SL_DATA18
#endif
#ifdef near__p_SL_DATA19_IGNORE
# undef near__p_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro near__p_sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef near__p_sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef near__p_sl_data0_weight
# undef near__p_sl_data1_weight
# undef near__p_sl_data2_weight
# undef near__p_sl_data3_weight
# undef near__p_sl_data4_weight
# undef near__p_sl_data5_weight
# undef near__p_sl_data6_weight
# undef near__p_sl_data7_weight
# undef near__p_sl_data8_weight
# undef near__p_sl_data9_weight
# undef near__p_sl_data10_weight
# undef near__p_sl_data11_weight
# undef near__p_sl_data12_weight
# undef near__p_sl_data13_weight
# undef near__p_sl_data14_weight
# undef near__p_sl_data15_weight
# undef near__p_sl_data16_weight
# undef near__p_sl_data17_weight
# undef near__p_sl_data18_weight
# undef near__p_sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable near__p_sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(near__p_sl_data0_weight) && !defined(near__p_SL_DATA0)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data1_weight) && !defined(near__p_SL_DATA1)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data2_weight) && !defined(near__p_SL_DATA2)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data3_weight) && !defined(near__p_SL_DATA3)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data4_weight) && !defined(near__p_SL_DATA4)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data5_weight) && !defined(near__p_SL_DATA5)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data6_weight) && !defined(near__p_SL_DATA6)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data7_weight) && !defined(near__p_SL_DATA7)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data8_weight) && !defined(near__p_SL_DATA8)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data9_weight) && !defined(near__p_SL_DATA9)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data10_weight) && !defined(near__p_SL_DATA10)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data11_weight) && !defined(near__p_SL_DATA11)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data12_weight) && !defined(near__p_SL_DATA12)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data13_weight) && !defined(near__p_SL_DATA13)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data14_weight) && !defined(near__p_SL_DATA14)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data15_weight) && !defined(near__p_SL_DATA15)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data16_weight) && !defined(near__p_SL_DATA16)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data17_weight) && !defined(near__p_SL_DATA17)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data18_weight) && !defined(near__p_SL_DATA18)
# undef near__p_sl_elem_weight
#endif
#if defined(near__p_sl_data19_weight) && !defined(near__p_SL_DATA19)
# undef near__p_sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro near__p_FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef near__p_SL_DATA0
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data0_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA1
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data1_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA2
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data2_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA3
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data3_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA4
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data4_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA5
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data5_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA6
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data6_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA7
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data7_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA8
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data8_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA9
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data9_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA10
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data10_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA11
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data11_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA12
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data12_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA13
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data13_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA14
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data14_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA15
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data15_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA16
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data16_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA17
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data17_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA18
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data18_flex
#   define near__p_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near__p_SL_DATA19
# ifdef near__p_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near__p_sl_data19_flex
#   define near__p_FLECKS_GUARD
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






#define near__p_SPEC_TLOC

typedef near__p_sl_int_type_c near__p_spec_int_t;

typedef int near__p_spec_proc_t;

#define near__p_SPEC_LOC_NONE   -1
#define near__p_SPEC_PROC_NONE  MPI_PROC_NULL

typedef void *near__p_spec_tloc_data_t;
typedef void *near__p_spec_tproc_data_t;

struct near__p__elements_t;

typedef struct near__p__elements_t *near__p_spec_elem_buf_t;

typedef struct near__p__elements_t near__p_spec_elem_t;

typedef near__p_sl_int_type_c near__p_spec_elem_index_t;

#define near__p_spec_elem_set_n(_e_, _n_)     near__p_elem_set_size((_e_), (_n_))
#define near__p_spec_elem_get_n(_e_)          near__p_elem_get_size((_e_))
#define near__p_spec_elem_set_nmax(_e_, _n_)  near__p_elem_set_max_size((_e_), (_n_))
#define near__p_spec_elem_get_nmax(_e_)       near__p_elem_get_max_size((_e_))

#define near__p_spec_elem_set_buf(_e_, _b_)   *(_e_) = *(_b_)
#define near__p_spec_elem_get_buf(_e_)        (_e_)

#define near__p_spec_elem_copy_at(_se_, _sat_, _de_, _dat_) \
  elem_copy_at((_se_), (_sat_), (_de_), (_dat_))

#define near__p_spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) \
  elem_xchange_at((_s0_), (_s0at_), (_s1_), (_s1at_), (_t_))






/* tproc count */

/* sp_macro near__p_SPEC_DECLARE_TPROC_COUNT_DB */
#define near__p_SPEC_DECLARE_TPROC_COUNT_DB \
  struct { near__p_spec_elem_index_t i; near__p_spec_proc_t p; } spec0cd;

/* sp_macro near__p_SPEC_DO_TPROC_COUNT_DB */
#define near__p_SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec0cd.i = 0; spec0cd.i < near__p_spec_elem_get_n(_b_); ++spec0cd.i) { \
    spec0cd.p = (_tp_)(near__p_spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == near__p_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROC_COUNT_DB */
#define near__p_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts) \
{ \
  near__p_SPEC_DECLARE_TPROC_COUNT_DB \
  near__p_SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro near__p_SPEC_DECLARE_TPROC_COUNT_IP */
#define near__p_SPEC_DECLARE_TPROC_COUNT_IP \
  struct { near__p_spec_elem_index_t i, t; near__p_spec_proc_t p; } spec0ci;

/* sp_macro near__p_SPEC_DO_TPROC_COUNT_IP */
#define near__p_SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < near__p_spec_elem_get_n(_b_); ++spec0ci.i) { \
    spec0ci.p = (_tp_)(near__p_spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == near__p_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) near__p_spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  near__p_spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro near__p_SPEC_FUNC_TPROC_COUNT_IP */
#define near__p_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts) \
{ \
  near__p_SPEC_DECLARE_TPROC_COUNT_IP \
  near__p_SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro near__p_SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define near__p_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { near__p_spec_elem_index_t i; near__p_spec_proc_t p; } spec1cd;

/* sp_macro near__p_SPEC_DO_TPROC_MOD_COUNT_DB */
#define near__p_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec1cd.i = 0; spec1cd.i < near__p_spec_elem_get_n(_b_); ++spec1cd.i) { \
    spec1cd.p = (_tp_)(near__p_spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == near__p_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define near__p_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts) \
{ \
  near__p_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  near__p_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro near__p_SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define near__p_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { near__p_spec_elem_index_t i, t; near__p_spec_proc_t p; } spec1ci;

/* sp_macro near__p_SPEC_DO_TPROC_MOD_COUNT_IP */
#define near__p_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < near__p_spec_elem_get_n(_b_); ++spec1ci.i) { \
    spec1ci.p = (_tp_)(near__p_spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == near__p_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) near__p_spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  near__p_spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro near__p_SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define near__p_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts) \
{ \
  near__p_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  near__p_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro near__p_SPEC_DECLARE_TPROCS_COUNT_DB */
#define near__p_SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { near__p_spec_elem_index_t i; near__p_spec_int_t j, n; } spec2cd;

/* sp_macro near__p_SPEC_DO_TPROCS_COUNT_DB */
#define near__p_SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec2cd.i = 0; spec2cd.i < near__p_spec_elem_get_n(_b_); ++spec2cd.i) { \
    spec2cd.n = (_tp_)(near__p_spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.n; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROCS_COUNT_DB */
#define near__p_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts, near__p_spec_proc_t *procs) \
{ \
  near__p_SPEC_DECLARE_TPROCS_COUNT_DB \
  near__p_SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro near__p_SPEC_DECLARE_TPROCS_COUNT_IP */
#define near__p_SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { near__p_spec_elem_index_t i, t; near__p_spec_int_t j, n; } spec2ci;

/* sp_macro near__p_SPEC_DO_TPROCS_COUNT_IP */
#define near__p_SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < near__p_spec_elem_get_n(_b_); ++spec2ci.i) { \
    spec2ci.n = (_tp_)(near__p_spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), (_ps_)); \
    if (spec2ci.n <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.n; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) near__p_spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  near__p_spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro near__p_SPEC_FUNC_TPROCS_COUNT_IP */
#define near__p_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts, near__p_spec_proc_t *procs) \
{ \
  near__p_SPEC_DECLARE_TPROCS_COUNT_IP \
  near__p_SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro near__p_SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define near__p_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { near__p_spec_elem_index_t i; near__p_spec_int_t j, n; } spec3cd;

/* sp_macro near__p_SPEC_DO_TPROCS_MOD_COUNT_DB */
#define near__p_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec3cd.i = 0; spec3cd.i < near__p_spec_elem_get_n(_b_); ++spec3cd.i) \
  { \
    spec3cd.n = (_tp_)(near__p_spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.n; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define near__p_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts, near__p_spec_proc_t *procs) \
{ \
  near__p_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  near__p_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro near__p_SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define near__p_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { near__p_spec_elem_index_t i, t; near__p_spec_int_t j, n; } spec3ci;

/* sp_macro near__p_SPEC_DO_TPROCS_MOD_COUNT_IP */
#define near__p_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < near__p_spec_elem_get_n(_b_); ++spec3ci.i) { \
    spec3ci.n = (_tp_)(near__p_spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), (_ps_), NULL); \
    if (spec3ci.n <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.n; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) near__p_spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  near__p_spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro near__p_SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define near__p_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts, near__p_spec_proc_t *procs) \
{ \
  near__p_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  near__p_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tproc rearrange */

/* sp_macro near__p_SPEC_DECLARE_TPROC_REARRANGE_DB */
#define near__p_SPEC_DECLARE_TPROC_REARRANGE_DB \
  struct { near__p_spec_elem_index_t i; near__p_spec_proc_t p; } spec0d;

/* sp_macro near__p_SPEC_DO_TPROC_REARRANGE_DB */
#define near__p_SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  for (spec0d.i = 0; spec0d.i < near__p_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tp_)(near__p_spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == near__p_SPEC_PROC_NONE) continue; \
    near__p_spec_elem_copy_at((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p]); \
    ++(_ds_)[spec0d.p]; \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROC_REARRANGE_DB */
#define near__p_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tproc_data_t tproc_data, int *displs) \
{ \
  near__p_SPEC_DECLARE_TPROC_REARRANGE_DB \
  near__p_SPEC_DO_TPROC_REARRANGE_DB(_tp_, tproc_data, s, d, displs); \
}

/* sp_macro near__p_SPEC_DECLARE_TPROC_REARRANGE_IP */
#define near__p_SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { near__p_spec_elem_index_t e, i, j; near__p_spec_proc_t p, np; } spec0i;

/* sp_macro near__p_SPEC_DO_TPROC_REARRANGE_IP */
#define near__p_SPEC_DO_TPROC_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_)  do { \
  for (spec0i.e = 0, spec0i.i = 0; spec0i.i < (_n_); ++spec0i.i) { \
    spec0i.e += (_cs_)[spec0i.i]; \
    spec0i.j = (_ds_)[spec0i.i]; \
    while (spec0i.j < spec0i.e) { \
      spec0i.p = (_tp_)(near__p_spec_elem_get_buf(_b_), spec0i.j, _tpd_); \
      while (spec0i.p != spec0i.i) { \
        spec0i.np = (_tp_)(near__p_spec_elem_get_buf(_b_), (_ds_)[spec0i.p], _tpd_); \
        if (spec0i.np != spec0i.p) near__p_spec_elem_exchange_at((_b_), (_ds_)[spec0i.p], (_b_), spec0i.j, (_xb_)); \
        ++(_ds_)[spec0i.p]; \
        spec0i.p = spec0i.np; \
      } \
      ++spec0i.j; \
    } \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROC_REARRANGE_IP */
#define near__p_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_rearrange_ip(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tproc_data_t tproc_data, int *displs, int *counts, near__p_spec_int_t n) \
{ \
  near__p_SPEC_DECLARE_TPROC_REARRANGE_IP \
  near__p_SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro near__p_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define near__p_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  struct { near__p_spec_elem_index_t i; near__p_spec_proc_t p; } spec1d;

/* sp_macro near__p_SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define near__p_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < near__p_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(near__p_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, near__p_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == near__p_SPEC_PROC_NONE) continue; \
      near__p_spec_elem_copy_at((_ib_), 0, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < near__p_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(near__p_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, NULL); \
      if (spec1d.p == near__p_SPEC_PROC_NONE) continue; \
      near__p_spec_elem_copy_at((_sb_), spec1d.i, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define near__p_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tproc_data_t tproc_data, int *displs, near__p_spec_elem_t *mod) \
{ \
  near__p_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  near__p_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro near__p_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define near__p_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { near__p_spec_elem_index_t e, i, j; near__p_spec_proc_t p, np; } spec1i;

/* sp_macro near__p_SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define near__p_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
      spec1i.e += (_cs_)[spec1i.i]; \
      spec1i.j = (_ds_)[spec1i.i]; \
      while (spec1i.j < spec1i.e) { \
        spec1i.p = (_tp_)(near__p_spec_elem_get_buf(_b_), spec1i.j, _tpd_, near__p_spec_elem_get_buf(_ib_)); \
        near__p_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(near__p_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, near__p_spec_elem_get_buf(_ib_)); \
          if (spec1i.np != spec1i.p) { \
            near__p_spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
            near__p_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
          } else near__p_spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
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
        spec1i.p = (_tp_)(near__p_spec_elem_get_buf(_b_), spec1i.j, _tpd_, NULL); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(near__p_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, NULL); \
          if (spec1i.np != spec1i.p) near__p_spec_elem_exchange_at((_b_), (_ds_)[spec1i.p], (_b_), spec1i.j, (_xb_)); \
          ++(_ds_)[spec1i.p]; \
          spec1i.p = spec1i.np; \
        } \
        ++spec1i.j; \
      } \
    } \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define near__p_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tproc_data_t tproc_data, int *displs, int *counts, near__p_spec_int_t n, near__p_spec_elem_t *mod) \
{ \
  near__p_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  near__p_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro near__p_SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define near__p_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  struct { near__p_spec_elem_index_t i; near__p_spec_int_t j, n; } spec2d;

/* sp_macro near__p_SPEC_DO_TPROCS_REARRANGE_DB */
#define near__p_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  for (spec2d.i = 0; spec2d.i < near__p_spec_elem_get_n(_sb_); ++spec2d.i) { \
    spec2d.n = (_tp_)(near__p_spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.n; ++spec2d.j) { \
      near__p_spec_elem_copy_at((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]]); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROCS_REARRANGE_DB */
#define near__p_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tproc_data_t tproc_data, int *displs, near__p_spec_proc_t *procs) \
{ \
  near__p_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  near__p_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro near__p_SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define near__p_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { near__p_spec_elem_index_t e, j, fe, fc, le, lc; near__p_spec_int_t i, n, f, l, o; } spec2i;

/* sp_macro near__p_SPEC_DO_TPROCS_REARRANGE_IP */
#define near__p_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_)  do { \
  spec2i.f = 0; spec2i.fe = (_cs_)[0]; spec2i.fc = near__p_spec_elem_get_n(_b_); \
  while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; } \
  spec2i.l = 0; spec2i.le = (_cs_)[0]; spec2i.lc = near__p_spec_elem_get_n(_b_) - 1; \
  while (spec2i.lc >= spec2i.le) { ++spec2i.l; spec2i.le += (_cs_)[spec2i.l]; } \
  for (spec2i.e = 0, spec2i.i = 0; spec2i.i < (_n_); ++spec2i.i) { \
    spec2i.e += (_cs_)[spec2i.i]; \
    spec2i.j = (_ds_)[spec2i.i]; \
    while (spec2i.j < spec2i.e) { \
      spec2i.n = (_tp_)(near__p_spec_elem_get_buf(_b_), spec2i.j, (_tpd_), (_ps_)); \
      spec2i.o = -1; \
      while (spec2i.n > 0) { \
        --spec2i.n; \
        if ((_ps_)[spec2i.n] == spec2i.i && spec2i.o < 0) spec2i.o = spec2i.n; \
        else if ((_ds_)[(_ps_)[spec2i.n]] < spec2i.fc) { \
          spec2i.l = spec2i.f; spec2i.le = spec2i.fe; spec2i.lc = spec2i.fc; \
          if (spec2i.fc < spec2i.fe) { \
            near__p_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_b_), spec2i.fc); \
            ++spec2i.fc; \
          } else near__p_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec2i.n]] == spec2i.fc) ++spec2i.fc; \
        if (spec2i.j != (_ds_)[(_ps_)[spec2i.n]]) near__p_spec_elem_copy_at((_b_), spec2i.j, (_b_), (_ds_)[(_ps_)[spec2i.n]]); \
        ++(_ds_)[(_ps_)[spec2i.n]]; \
        while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; spec2i.fc = (_ds_)[spec2i.f]; } \
      } \
      if (spec2i.o < 0) { \
        if (spec2i.lc < spec2i.le) {  \
          near__p_spec_elem_copy_at((_b_), spec2i.lc, (_b_), spec2i.j); \
          spec2i.f = spec2i.l; spec2i.fe = spec2i.le; spec2i.fc = spec2i.lc; \
          --spec2i.lc; \
          while (spec2i.l > 0 && spec2i.lc < (_ds_)[spec2i.l]) { spec2i.le -= (_cs_)[spec2i.l]; spec2i.lc = spec2i.le - 1; --spec2i.l; } \
        } else near__p_spec_elem_copy_at((_xb_), 0, (_b_), spec2i.j); \
      } \
      spec2i.j = (_ds_)[spec2i.i]; \
    } \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROCS_REARRANGE_IP */
#define near__p_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_ip(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tproc_data_t tproc_data, int *displs, int *counts, near__p_spec_int_t n, near__p_spec_proc_t *procs) \
{ \
  near__p_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  near__p_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro near__p_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define near__p_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  struct { near__p_spec_elem_index_t i; near__p_spec_int_t j, n; } spec3d;

/* sp_macro near__p_SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define near__p_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  if (_ib_) { \
    for (spec3d.i = 0; spec3d.i < near__p_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(near__p_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), near__p_spec_elem_get_buf(_ib_)); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        near__p_spec_elem_copy_at((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } else { \
    for (spec3d.i = 0; spec3d.i < near__p_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(near__p_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), NULL); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        near__p_spec_elem_copy_at((_sb_), spec3d.i, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define near__p_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tproc_data_t tproc_data, int *displs, near__p_spec_proc_t *procs, near__p_spec_elem_t *mod) \
{ \
  near__p_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  near__p_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro near__p_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define near__p_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { near__p_spec_elem_index_t e, j, fe, fc, le, lc; near__p_spec_int_t i, n, f, l, o; } spec3i;

/* sp_macro near__p_SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define near__p_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  if (_ib_) { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = near__p_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = near__p_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(near__p_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), near__p_spec_elem_get_buf(_ib_)); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              near__p_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else near__p_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          near__p_spec_elem_copy_at((_ib_), spec3i.n, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            near__p_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else near__p_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } else { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = near__p_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = near__p_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(near__p_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), NULL); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              near__p_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else near__p_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          if (spec3i.j != (_ds_)[(_ps_)[spec3i.n]]) near__p_spec_elem_copy_at((_b_), spec3i.j, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            near__p_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else near__p_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define near__p_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tproc_data_t tproc_data, int *displs, int *counts, near__p_spec_int_t n, near__p_spec_proc_t *procs, near__p_spec_elem_t *mod) \
{ \
  near__p_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  near__p_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}

/* sp_macro near__p_SPEC_DEFINE_TPROC */
#define near__p_SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  near__p_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro near__p_SPEC_DEFINE_TPROC_MOD */
#define near__p_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  near__p_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro near__p_SPEC_DEFINE_TPROCS */
#define near__p_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  near__p_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro near__p_SPEC_DEFINE_TPROCS_MOD */
#define near__p_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  near__p_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  near__p_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro near__p_SPEC_EXT_PARAM_TPROC near__p_SPEC_EXT_PARAM_TPROC_NULL near__p_SPEC_EXT_PARAM_TPROC_MOD near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL near__p_SPEC_EXT_PARAM_TPROCS near__p_SPEC_EXT_PARAM_TPROCS_NULL near__p_SPEC_EXT_PARAM_TPROCS_MOD near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define near__p_SPEC_EXT_PARAM_TPROC(_name_)       _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip
#define near__p_SPEC_EXT_PARAM_TPROC_NULL          NULL, NULL, NULL, NULL
#define near__p_SPEC_EXT_PARAM_TPROC_MOD(_name_)   _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip
#define near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL      NULL, NULL, NULL, NULL
#define near__p_SPEC_EXT_PARAM_TPROCS(_name_)      _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip
#define near__p_SPEC_EXT_PARAM_TPROCS_NULL         NULL, NULL, NULL, NULL
#define near__p_SPEC_EXT_PARAM_TPROCS_MOD(_name_)  _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip
#define near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL     NULL, NULL, NULL, NULL


/* sp_type near__p_spec_tproc_f near__p_spec_tproc_count_f near__p_spec_tproc_rearrange_db_f near__p_spec_tproc_rearrange_ip_f */
typedef near__p_spec_proc_t near__p_spec_tproc_f(near__p_spec_elem_buf_t b, near__p_spec_elem_index_t x, near__p_spec_tproc_data_t tproc_data);
typedef void near__p_spec_tproc_count_f(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts);
typedef void near__p_spec_tproc_rearrange_db_f(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tproc_data_t tproc_data, int *displs);
typedef void near__p_spec_tproc_rearrange_ip_f(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tproc_data_t tproc_data, int *displs, int *counts, near__p_spec_int_t n);

/* sp_type near__p_spec_tproc_mod_f near__p_spec_tproc_mod_count_f near__p_spec_tproc_mod_rearrange_db_f near__p_spec_tproc_mod_rearrange_ip_f */
typedef near__p_spec_proc_t near__p_spec_tproc_mod_f(near__p_spec_elem_buf_t b, near__p_spec_elem_index_t x, near__p_spec_tproc_data_t tproc_data, near__p_spec_elem_buf_t mod);
typedef void near__p_spec_tproc_mod_count_f(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts);
typedef void near__p_spec_tproc_mod_rearrange_db_f(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tproc_data_t tproc_data, int *displs, near__p_spec_elem_t *mod);
typedef void near__p_spec_tproc_mod_rearrange_ip_f(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tproc_data_t tproc_data, int *displs, int *counts, near__p_spec_int_t n, near__p_spec_elem_t *mod);

/* sp_type near__p_spec_tprocs_f near__p_spec_tprocs_count_f near__p_spec_tprocs_rearrange_db_f near__p_spec_tprocs_rearrange_ip_f */
typedef near__p_spec_int_t near__p_spec_tprocs_f(near__p_spec_elem_buf_t b, near__p_spec_elem_index_t x, near__p_spec_tproc_data_t tproc_data, near__p_spec_proc_t *procs);
typedef void near__p_spec_tprocs_count_f(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts, near__p_spec_proc_t *procs);
typedef void near__p_spec_tprocs_rearrange_db_f(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tproc_data_t tproc_data, int *displs, near__p_spec_proc_t *procs);
typedef void near__p_spec_tprocs_rearrange_ip_f(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tproc_data_t tproc_data, int *displs, int *counts, near__p_spec_int_t n, near__p_spec_proc_t *procs);

/* sp_type near__p_spec_tprocs_mod_f near__p_spec_tprocs_mod_count_f near__p_spec_tprocs_mod_rearrange_db_f near__p_spec_tprocs_mod_rearrange_ip_f */
typedef near__p_spec_int_t near__p_spec_tprocs_mod_f(near__p_spec_elem_buf_t b, near__p_spec_elem_index_t x, near__p_spec_tproc_data_t tproc_data, near__p_spec_proc_t *procs, near__p_spec_elem_buf_t mod);
typedef void near__p_spec_tprocs_mod_count_f(near__p_spec_elem_t *s, near__p_spec_tproc_data_t tproc_data, int *counts, near__p_spec_proc_t *procs);
typedef void near__p_spec_tprocs_mod_rearrange_db_f(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tproc_data_t tproc_data, int *displs, near__p_spec_proc_t *procs, near__p_spec_elem_t *mod);
typedef void near__p_spec_tprocs_mod_rearrange_ip_f(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tproc_data_t tproc_data, int *displs, int *counts, near__p_spec_int_t n, near__p_spec_proc_t *procs, near__p_spec_elem_t *mod);

/* sp_type near__p_spec_tproc_reset_f */
typedef void near__p_spec_tproc_reset_f(near__p_spec_tproc_data_t tproc_data);


/* enable tloc features */
#ifdef near__p_SPEC_TLOC

/* sp_macro near__p_SPEC_TLOC near__p_SPEC_LOC_NONE */


/* tloc rearrange */

/* sp_macro near__p_SPEC_DECLARE_TLOC_REARRANGE_DB */
#define near__p_SPEC_DECLARE_TLOC_REARRANGE_DB \
  struct { near__p_spec_int_t i, p; } spec0d;

/* sp_macro near__p_SPEC_DO_TLOC_REARRANGE_DB */
#define near__p_SPEC_DO_TLOC_REARRANGE_DB(_tl_, _tld_, _sb_, _db_)  do { \
  for (spec0d.i = 0; spec0d.i < near__p_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tl_)(near__p_spec_elem_get_buf(_sb_), spec0d.i, _tld_); \
    if (spec0d.p == near__p_SPEC_LOC_NONE) continue; \
    near__p_spec_elem_copy_at((_sb_), spec0d.i, (_db_), spec0d.p); \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TLOC_REARRANGE_DB */
#define near__p_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_rearrange_db(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tloc_data_t tloc_data) \
{ \
  near__p_SPEC_DECLARE_TLOC_REARRANGE_DB \
  near__p_SPEC_DO_TLOC_REARRANGE_DB(_tl_, tloc_data, s, d); \
}

/* sp_macro near__p_SPEC_DECLARE_TLOC_REARRANGE_IP */
#define near__p_SPEC_DECLARE_TLOC_REARRANGE_IP \
  struct { near__p_spec_int_t i, p, np; } spec0i;

/* sp_macro near__p_SPEC_DO_TLOC_REARRANGE_IP */
#define near__p_SPEC_DO_TLOC_REARRANGE_IP(_tl_, _tld_, _b_, _xb_)  do { \
  for (spec0i.i = 0; spec0i.i < near__p_spec_elem_get_n(_b_); ++spec0i.i) { \
    spec0i.p = (_tl_)(near__p_spec_elem_get_buf(_b_), spec0i.i, _tld_); \
    if (spec0i.p == near__p_SPEC_LOC_NONE) continue; \
    while (spec0i.i != spec0i.p) { \
      spec0i.np = (_tl_)(near__p_spec_elem_get_buf(_b_), spec0i.p, _tld_); \
      if (spec0i.np == near__p_SPEC_LOC_NONE) { near__p_spec_elem_copy_at((_b_), spec0i.i, (_b_), spec0i.p); break; } \
      near__p_spec_elem_exchange_at((_b_), spec0i.i, (_b_), spec0i.p, (_xb_)); \
      spec0i.p = spec0i.np; \
    } \
  } } while (0)

/* sp_macro near__p_SPEC_FUNC_TLOC_REARRANGE_IP */
#define near__p_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_rearrange_ip(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tloc_data_t tloc_data) \
{ \
  near__p_SPEC_DECLARE_TLOC_REARRANGE_IP \
  near__p_SPEC_DO_TLOC_REARRANGE_IP(_tl_, tloc_data, s, x); \
}


/* tloc_mod_mod rearrange */

/* sp_macro near__p_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB */
#define near__p_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  struct { near__p_spec_int_t i, p; } spec1d;

/* sp_macro near__p_SPEC_DO_TLOC_MOD_REARRANGE_DB */
#define near__p_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, _tld_, _sb_, _db_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < near__p_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(near__p_spec_elem_get_buf(_sb_), spec1d.i, _tld_, near__p_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == near__p_SPEC_LOC_NONE) continue; \
      near__p_spec_elem_copy_at((_ib_), 0, (_db_), spec1d.p); \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < near__p_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(near__p_spec_elem_get_buf(_sb_), spec1d.i, _tld_, NULL); \
      if (spec1d.p == near__p_SPEC_LOC_NONE) continue; \
      near__p_spec_elem_copy_at((_sb_), spec1d.i, (_db_), spec1d.p); \
    } \
  } } while (0) 

/* sp_macro near__p_SPEC_FUNC_TLOC_MOD_REARRANGE_DB */
#define near__p_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_mod_rearrange_db(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tloc_data_t tloc_data, near__p_spec_elem_t *mod) \
{ \
  near__p_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  near__p_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, tloc_data, s, d, mod); \
}

/* sp_macro near__p_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP */
#define near__p_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  struct { near__p_spec_int_t i, p, np; } spec1i;

/* sp_macro near__p_SPEC_DO_TLOC_MOD_REARRANGE_IP */
#define near__p_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, _tld_, _b_, _xb_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.i = 0; spec1i.i < near__p_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(near__p_spec_elem_get_buf(_b_), spec1i.i, _tld_, near__p_spec_elem_get_buf(_ib_)); \
      if (spec1i.p == near__p_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(near__p_spec_elem_get_buf(_b_), spec1i.p, _tld_, near__p_spec_elem_get_buf(_xb_)); \
        if (spec1i.np == near__p_SPEC_LOC_NONE) break; \
        near__p_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.p); \
        near__p_spec_elem_copy_at((_xb_), 0, (_ib_), 0); \
        spec1i.p = spec1i.np; \
      } \
      near__p_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.i); \
    } \
  } else { \
    for (spec1i.i = 0; spec1i.i < near__p_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(near__p_spec_elem_get_buf(_b_), spec1i.i, _tld_, NULL); \
      if (spec1i.p == near__p_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(near__p_spec_elem_get_buf(_b_), spec1i.p, _tld_, NULL); \
        if (spec1i.np == near__p_SPEC_LOC_NONE) { near__p_spec_elem_copy_at((_b_), spec1i.i, (_b_), spec1i.p); break; } \
        near__p_spec_elem_exchange_at((_b_), spec1i.i, (_b_), spec1i.p, (_xb_)); \
        spec1i.p = spec1i.np; \
      } \
    } \
 } } while (0) 

/* sp_macro near__p_SPEC_FUNC_TLOC_MOD_REARRANGE_IP */
#define near__p_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_mod_rearrange_ip(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tloc_data_t tloc_data, near__p_spec_elem_t *mod) \
{ \
  near__p_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  near__p_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, tloc_data, s, x, mod); \
}

/* sp_macro near__p_SPEC_DEFINE_TLOC */
#define near__p_SPEC_DEFINE_TLOC(_name_, _tl_, _s_...) \
  near__p_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_) \
  near__p_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro near__p_SPEC_DEFINE_TLOC_MOD */
#define near__p_SPEC_DEFINE_TLOC_MOD(_name_, _tl_, _s_...) \
  near__p_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_) \
  near__p_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro near__p_SPEC_EXT_PARAM_TLOC near__p_SPEC_EXT_PARAM_TLOC_NULL near__p_SPEC_EXT_PARAM_TLOC_MOD near__p_SPEC_EXT_PARAM_TLOC_MOD_NULL */
#define near__p_SPEC_EXT_PARAM_TLOC(_name_)      _name_##_tloc_rearrange_db, _name_##_tloc_rearrange_ip
#define near__p_SPEC_EXT_PARAM_TLOC_NULL         NULL, NULL
#define near__p_SPEC_EXT_PARAM_TLOC_MOD(_name_)  _name_##_tloc_mod_rearrange_db, _name_##_tloc_mod_rearrange_ip
#define near__p_SPEC_EXT_PARAM_TLOC_MOD_NULL     NULL, NULL


/* sp_type near__p_spec_tloc_f near__p_spec_tloc_rearrange_db_f near__p_spec_tloc_rearrange_ip_f */
typedef near__p_spec_elem_index_t near__p_spec_tloc_f(near__p_spec_elem_buf_t b, near__p_spec_elem_index_t x, near__p_spec_tloc_data_t tloc_data);
typedef void near__p_spec_tloc_rearrange_db_f(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tloc_data_t tloc_data);
typedef void near__p_spec_tloc_rearrange_ip_f(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tloc_data_t tloc_data);

/* sp_type near__p_spec_tloc_mod_f near__p_spec_tloc_mod_rearrange_db_f near__p_spec_tloc_mod_rearrange_ip_f */
typedef near__p_spec_elem_index_t near__p_spec_tloc_mod_f(near__p_spec_elem_buf_t b, near__p_spec_elem_index_t x, near__p_spec_tloc_data_t tloc_data, near__p_spec_elem_buf_t mod);
typedef void near__p_spec_tloc_mod_rearrange_db_f(near__p_spec_elem_t *s, near__p_spec_elem_t *d, near__p_spec_tloc_data_t tloc_data, near__p_spec_elem_t *mod);
typedef void near__p_spec_tloc_mod_rearrange_ip_f(near__p_spec_elem_t *s, near__p_spec_elem_t *x, near__p_spec_tloc_data_t tloc_data, near__p_spec_elem_t *mod);


#endif /* near__p_SPEC_TLOC */






#ifdef SL_USE_MPI
# include <mpi.h>
#endif


/* sl_type near__p_slint_t near__p_slint */
typedef near__p_sl_int_type_c near__p_slint_t, near__p_slint;  /* deprecated 'near__p_slint' */

#define near__p_slint_fmt   near__p_sl_int_type_fmt    /* sl_macro */

/* sl_type near__p_slindex_t */
typedef near__p_sl_index_type_c near__p_slindex_t;

#define near__p_sindex_fmt  near__p_sl_index_type_fmt  /* sl_macro */

/* sl_type near__p_slkey_t */
typedef near__p_sl_key_type_c near__p_slkey_t;

/* sl_type near__p_slkey_pure_t near__p_slpkey_t */
typedef near__p_sl_key_pure_type_c near__p_slkey_pure_t, near__p_slpkey_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type near__p_sldata0_t */
#ifdef near__p_sl_data0_type_c
typedef near__p_sl_data0_type_c near__p_sldata0_t;
#endif
/* sl_type near__p_sldata1_t */
#ifdef near__p_sl_data1_type_c
typedef near__p_sl_data1_type_c near__p_sldata1_t;
#endif
/* sl_type near__p_sldata2_t */
#ifdef near__p_sl_data2_type_c
typedef near__p_sl_data2_type_c near__p_sldata2_t;
#endif
/* sl_type near__p_sldata3_t */
#ifdef near__p_sl_data3_type_c
typedef near__p_sl_data3_type_c near__p_sldata3_t;
#endif
/* sl_type near__p_sldata4_t */
#ifdef near__p_sl_data4_type_c
typedef near__p_sl_data4_type_c near__p_sldata4_t;
#endif
/* sl_type near__p_sldata5_t */
#ifdef near__p_sl_data5_type_c
typedef near__p_sl_data5_type_c near__p_sldata5_t;
#endif
/* sl_type near__p_sldata6_t */
#ifdef near__p_sl_data6_type_c
typedef near__p_sl_data6_type_c near__p_sldata6_t;
#endif
/* sl_type near__p_sldata7_t */
#ifdef near__p_sl_data7_type_c
typedef near__p_sl_data7_type_c near__p_sldata7_t;
#endif
/* sl_type near__p_sldata8_t */
#ifdef near__p_sl_data8_type_c
typedef near__p_sl_data8_type_c near__p_sldata8_t;
#endif
/* sl_type near__p_sldata9_t */
#ifdef near__p_sl_data9_type_c
typedef near__p_sl_data9_type_c near__p_sldata9_t;
#endif
/* sl_type near__p_sldata10_t */
#ifdef near__p_sl_data10_type_c
typedef near__p_sl_data10_type_c near__p_sldata10_t;
#endif
/* sl_type near__p_sldata11_t */
#ifdef near__p_sl_data11_type_c
typedef near__p_sl_data11_type_c near__p_sldata11_t;
#endif
/* sl_type near__p_sldata12_t */
#ifdef near__p_sl_data12_type_c
typedef near__p_sl_data12_type_c near__p_sldata12_t;
#endif
/* sl_type near__p_sldata13_t */
#ifdef near__p_sl_data13_type_c
typedef near__p_sl_data13_type_c near__p_sldata13_t;
#endif
/* sl_type near__p_sldata14_t */
#ifdef near__p_sl_data14_type_c
typedef near__p_sl_data14_type_c near__p_sldata14_t;
#endif
/* sl_type near__p_sldata15_t */
#ifdef near__p_sl_data15_type_c
typedef near__p_sl_data15_type_c near__p_sldata15_t;
#endif
/* sl_type near__p_sldata16_t */
#ifdef near__p_sl_data16_type_c
typedef near__p_sl_data16_type_c near__p_sldata16_t;
#endif
/* sl_type near__p_sldata17_t */
#ifdef near__p_sl_data17_type_c
typedef near__p_sl_data17_type_c near__p_sldata17_t;
#endif
/* sl_type near__p_sldata18_t */
#ifdef near__p_sl_data18_type_c
typedef near__p_sl_data18_type_c near__p_sldata18_t;
#endif
/* sl_type near__p_sldata19_t */
#ifdef near__p_sl_data19_type_c
typedef near__p_sl_data19_type_c near__p_sldata19_t;
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

/* sl_type near__p_slweight_t */
typedef near__p_sl_weight_type_c near__p_slweight_t;

#define near__p_slweight_fmt  near__p_sl_weight_type_fmt  /* sl_macro */

#if defined(near__p_sl_elem_weight) && defined(near__p_sl_weight_intequiv)
typedef near__p_sl_weight_type_c near__p_slcount_t;       /* sl_type near__p_slcount_t */
# define near__p_slcount_fmt  near__p_sl_weight_type_fmt  /* sl_macro */
#else
typedef near__p_sl_int_type_c near__p_slcount_t;
# define near__p_slcount_fmt  near__p_sl_int_type_fmt
#endif


/* sl_type near__p__slpwkey_t near__p_slpwkey_t */
typedef struct near__p__slpwkey_t
{
  near__p_slpkey_t pkey;
  near__p_slweight_t weight;

} near__p_slpwkey_t;


/* sl_type near__p__elements_t near__p_elements_t */
typedef struct near__p__elements_t
{
  near__p_slint_t size, max_size;
  near__p_slkey_t *keys;

#ifdef near__p_SL_INDEX
  near__p_slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef near__p_SL_DATA0
  near__p_sldata0_t *data0;
#endif
#ifdef near__p_SL_DATA1
  near__p_sldata1_t *data1;
#endif
#ifdef near__p_SL_DATA2
  near__p_sldata2_t *data2;
#endif
#ifdef near__p_SL_DATA3
  near__p_sldata3_t *data3;
#endif
#ifdef near__p_SL_DATA4
  near__p_sldata4_t *data4;
#endif
#ifdef near__p_SL_DATA5
  near__p_sldata5_t *data5;
#endif
#ifdef near__p_SL_DATA6
  near__p_sldata6_t *data6;
#endif
#ifdef near__p_SL_DATA7
  near__p_sldata7_t *data7;
#endif
#ifdef near__p_SL_DATA8
  near__p_sldata8_t *data8;
#endif
#ifdef near__p_SL_DATA9
  near__p_sldata9_t *data9;
#endif
#ifdef near__p_SL_DATA10
  near__p_sldata10_t *data10;
#endif
#ifdef near__p_SL_DATA11
  near__p_sldata11_t *data11;
#endif
#ifdef near__p_SL_DATA12
  near__p_sldata12_t *data12;
#endif
#ifdef near__p_SL_DATA13
  near__p_sldata13_t *data13;
#endif
#ifdef near__p_SL_DATA14
  near__p_sldata14_t *data14;
#endif
#ifdef near__p_SL_DATA15
  near__p_sldata15_t *data15;
#endif
#ifdef near__p_SL_DATA16
  near__p_sldata16_t *data16;
#endif
#ifdef near__p_SL_DATA17
  near__p_sldata17_t *data17;
#endif
#ifdef near__p_SL_DATA18
  near__p_sldata18_t *data18;
#endif
#ifdef near__p_SL_DATA19
  near__p_sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} near__p_elements_t;


/* sl_type near__p__packed_element_t near__p_packed_element_t */
typedef struct near__p__packed_element_t
{
  near__p_slkey_t key;

#ifdef near__p_SL_PACKED_INDEX
  near__p_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef near__p_SL_DATA0
# ifdef near__p_sl_data0_flex
  near__p_sldata0_t data0[];
# else
  near__p_sldata0_t data0[near__p_sl_data0_size_c];
# endif
#endif
#ifdef near__p_SL_DATA1
# ifdef near__p_sl_data1_flex
  near__p_sldata1_t data1[];
# else
  near__p_sldata1_t data1[near__p_sl_data1_size_c];
# endif
#endif
#ifdef near__p_SL_DATA2
# ifdef near__p_sl_data2_flex
  near__p_sldata2_t data2[];
# else
  near__p_sldata2_t data2[near__p_sl_data2_size_c];
# endif
#endif
#ifdef near__p_SL_DATA3
# ifdef near__p_sl_data3_flex
  near__p_sldata3_t data3[];
# else
  near__p_sldata3_t data3[near__p_sl_data3_size_c];
# endif
#endif
#ifdef near__p_SL_DATA4
# ifdef near__p_sl_data4_flex
  near__p_sldata4_t data4[];
# else
  near__p_sldata4_t data4[near__p_sl_data4_size_c];
# endif
#endif
#ifdef near__p_SL_DATA5
# ifdef near__p_sl_data5_flex
  near__p_sldata5_t data5[];
# else
  near__p_sldata5_t data5[near__p_sl_data5_size_c];
# endif
#endif
#ifdef near__p_SL_DATA6
# ifdef near__p_sl_data6_flex
  near__p_sldata6_t data6[];
# else
  near__p_sldata6_t data6[near__p_sl_data6_size_c];
# endif
#endif
#ifdef near__p_SL_DATA7
# ifdef near__p_sl_data7_flex
  near__p_sldata7_t data7[];
# else
  near__p_sldata7_t data7[near__p_sl_data7_size_c];
# endif
#endif
#ifdef near__p_SL_DATA8
# ifdef near__p_sl_data8_flex
  near__p_sldata8_t data8[];
# else
  near__p_sldata8_t data8[near__p_sl_data8_size_c];
# endif
#endif
#ifdef near__p_SL_DATA9
# ifdef near__p_sl_data9_flex
  near__p_sldata9_t data9[];
# else
  near__p_sldata9_t data9[near__p_sl_data9_size_c];
# endif
#endif
#ifdef near__p_SL_DATA10
# ifdef near__p_sl_data10_flex
  near__p_sldata10_t data10[];
# else
  near__p_sldata10_t data10[near__p_sl_data10_size_c];
# endif
#endif
#ifdef near__p_SL_DATA11
# ifdef near__p_sl_data11_flex
  near__p_sldata11_t data11[];
# else
  near__p_sldata11_t data11[near__p_sl_data11_size_c];
# endif
#endif
#ifdef near__p_SL_DATA12
# ifdef near__p_sl_data12_flex
  near__p_sldata12_t data12[];
# else
  near__p_sldata12_t data12[near__p_sl_data12_size_c];
# endif
#endif
#ifdef near__p_SL_DATA13
# ifdef near__p_sl_data13_flex
  near__p_sldata13_t data13[];
# else
  near__p_sldata13_t data13[near__p_sl_data13_size_c];
# endif
#endif
#ifdef near__p_SL_DATA14
# ifdef near__p_sl_data14_flex
  near__p_sldata14_t data14[];
# else
  near__p_sldata14_t data14[near__p_sl_data14_size_c];
# endif
#endif
#ifdef near__p_SL_DATA15
# ifdef near__p_sl_data15_flex
  near__p_sldata15_t data15[];
# else
  near__p_sldata15_t data15[near__p_sl_data15_size_c];
# endif
#endif
#ifdef near__p_SL_DATA16
# ifdef near__p_sl_data16_flex
  near__p_sldata16_t data16[];
# else
  near__p_sldata16_t data16[near__p_sl_data16_size_c];
# endif
#endif
#ifdef near__p_SL_DATA17
# ifdef near__p_sl_data17_flex
  near__p_sldata17_t data17[];
# else
  near__p_sldata17_t data17[near__p_sl_data17_size_c];
# endif
#endif
#ifdef near__p_SL_DATA18
# ifdef near__p_sl_data18_flex
  near__p_sldata18_t data18[];
# else
  near__p_sldata18_t data18[near__p_sl_data18_size_c];
# endif
#endif
#ifdef near__p_SL_DATA19
# ifdef near__p_sl_data19_flex
  near__p_sldata19_t data19[];
# else
  near__p_sldata19_t data19[near__p_sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} near__p_packed_element_t;


/* sl_type near__p__packed_elements_t near__p_packed_elements_t */
typedef struct near__p__packed_elements_t
{
  near__p_slint_t size, max_size;
  
  near__p_packed_element_t *elements;
  
} near__p_packed_elements_t;


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


/* sl_type near__p__classification_info_t near__p_classification_info_t near__p_classification_info */
typedef struct near__p__classification_info_t
{
  near__p_slint_t nclasses;
  near__p_slkey_pure_t *keys;
  near__p_slint_t *counts;
  near__p_slint_t *masks;

  /* */
  near__p_slint_t *all_local_sizes;
  near__p_slint_t *local_lt_eq_counts;
  near__p_slint_t *all_local_lt_eq_counts;

} near__p_classification_info_t, near__p_classification_info;  /* deprecated 'near__p_classification_info' */


/* key2class, sl_type near__p_key2class_f */
typedef near__p_slint_t (*near__p_key2class_f)(near__p_slkey_t *, near__p_slint, void *);

/* pivot-element, sl_type near__p_pivot_f */
typedef near__p_slint_t (*near__p_pivot_f)(near__p_elements_t *);

/* sorting-network, sl_type near__p_sortnet_f near__p_sortnet_data_t */
typedef void *near__p_sortnet_data_t;
typedef near__p_slint_t (*near__p_sortnet_f)(near__p_slint_t size, near__p_slint_t rank, near__p_slint_t stage, near__p_sortnet_data_t snd, near__p_slint_t *up);

/* merge2, sl_type near__p_merge2x_f near__p_merge2X_f */
typedef near__p_slint_t (*near__p_merge2x_f)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
typedef near__p_slint_t (*near__p_merge2X_f)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx, near__p_elements_t *t);

/* sl_type near__p__permute_generic_t near__p_permute_generic_t */
typedef struct near__p__permute_generic_t
{
  int type;

  near__p_spec_tloc_f *tloc;
  near__p_spec_tloc_rearrange_db_f *tloc_rearrange_db;
  near__p_spec_tloc_rearrange_ip_f *tloc_rearrange_ip;

  near__p_spec_tloc_mod_f *tloc_mod;
  near__p_spec_tloc_mod_rearrange_db_f *tloc_mod_rearrange_db;
  near__p_spec_tloc_mod_rearrange_ip_f *tloc_mod_rearrange_ip;

} near__p_permute_generic_t;

/* sl_macro near__p_PERMUTE_GENERIC_DEFINE_TLOC near__p_PERMUTE_GENERIC_INIT_TLOC near__p_PERMUTE_GENERIC_INIT_EXT_TLOC */
#define near__p_PERMUTE_GENERIC_DEFINE_TLOC(_tl_, _s_...)      near__p_SPEC_DEFINE_TLOC(_tl_, _tl_, _s_)
#define near__p_PERMUTE_GENERIC_INIT_TLOC(_tl_)                { 1, _tl_, near__p_SPEC_EXT_PARAM_TLOC_NULL,  NULL, near__p_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define near__p_PERMUTE_GENERIC_INIT_EXT_TLOC(_tl_)            { 1, _tl_, near__p_SPEC_EXT_PARAM_TLOC(_tl_), NULL, near__p_SPEC_EXT_PARAM_TLOC_MOD_NULL }

/* sl_macro near__p_PERMUTE_GENERIC_DEFINE_TLOC_MOD near__p_PERMUTE_GENERIC_INIT_TLOC_MOD near__p_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD */
#define near__p_PERMUTE_GENERIC_DEFINE_TLOC_MOD(_tl_, _s_...)  near__p_SPEC_DEFINE_TLOC_MOD(_tl_, _tl_, _s_)
#define near__p_PERMUTE_GENERIC_INIT_TLOC_MOD(_tl_)            { 2, NULL, near__p_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, near__p_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define near__p_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD(_tl_)        { 2, NULL, near__p_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, near__p_SPEC_EXT_PARAM_TLOC_MOD(_tl_) }

/* sl_type near__p__split_generic_t near__p_split_generic_t */
typedef struct near__p__split_generic_t
{
  int type;

  near__p_spec_tproc_f *tproc;
  near__p_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  near__p_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  near__p_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  near__p_spec_tproc_mod_f *tproc_mod;
  near__p_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  near__p_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  near__p_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  near__p_spec_tprocs_f *tprocs;
  near__p_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  near__p_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  near__p_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  near__p_spec_tprocs_mod_f *tprocs_mod;
  near__p_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  near__p_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  near__p_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

  near__p_spec_tproc_reset_f *reset;

} near__p_split_generic_t;

/* sl_macro near__p_SPLIT_GENERIC_DEFINE_TPROC near__p_SPLIT_GENERIC_INIT_TPROC near__p_SPLIT_GENERIC_INIT_EXT_TPROC */
#define near__p_SPLIT_GENERIC_DEFINE_TPROC(_tp_, _s_...)         near__p_SPEC_DEFINE_TPROC(_tp_, _tp_, _s_)
#define near__p_SPLIT_GENERIC_INIT_TPROC(_tp_, _r_...)           { 1, _tp_, near__p_SPEC_EXT_PARAM_TPROC_NULL,  NULL, near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROCS_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define near__p_SPLIT_GENERIC_INIT_EXT_TPROC(_tp_, _r_...)       { 1, _tp_, near__p_SPEC_EXT_PARAM_TPROC(_tp_), NULL, near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROCS_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro near__p_SPLIT_GENERIC_DEFINE_TPROC_MOD near__p_SPLIT_GENERIC_INIT_TPROC_MOD near__p_SPLIT_GENERIC_INIT_EXT_TPROC_MOD */
#define near__p_SPLIT_GENERIC_DEFINE_TPROC_MOD(_tp_, _s_...)     near__p_SPEC_DEFINE_TPROC_MOD(_tp_, _tp_, _s_)
#define near__p_SPLIT_GENERIC_INIT_TPROC_MOD(_tp_, _r_...)       { 2, NULL, near__p_SPEC_EXT_PARAM_TPROC_NULL, _tp_, near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL,  NULL, near__p_SPEC_EXT_PARAM_TPROCS_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define near__p_SPLIT_GENERIC_INIT_EXT_TPROC_MOD(_tp_, _r_...)   { 2, NULL, near__p_SPEC_EXT_PARAM_TPROC_NULL, _tp_, near__p_SPEC_EXT_PARAM_TPROC_MOD(_tp_), NULL, near__p_SPEC_EXT_PARAM_TPROCS_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro near__p_SPLIT_GENERIC_DEFINE_TPROCS near__p_SPLIT_GENERIC_INIT_TPROCS near__p_SPLIT_GENERIC_INIT_EXT_TPROCS */
#define near__p_SPLIT_GENERIC_DEFINE_TPROCS(_tp_, _s_...)        near__p_SPEC_DEFINE_TPROCS(_tp_, _tp_, _s_)
#define near__p_SPLIT_GENERIC_INIT_TPROCS(_tp_, _r_...)          { 3, NULL, near__p_SPEC_EXT_PARAM_TPROC_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, near__p_SPEC_EXT_PARAM_TPROCS_NULL,  NULL, near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define near__p_SPLIT_GENERIC_INIT_EXT_TPROCS(_tp_, _r_...)      { 3, NULL, near__p_SPEC_EXT_PARAM_TPROC_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, near__p_SPEC_EXT_PARAM_TPROCS(_tp_), NULL, near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro near__p_SPLIT_GENERIC_DEFINE_TPROCS_MOD near__p_SPLIT_GENERIC_INIT_TPROCS_MOD near__p_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD */
#define near__p_SPLIT_GENERIC_DEFINE_TPROCS_MOD(_tp_, _s_...)    near__p_SPEC_DEFINE_TPROCS_MOD(_tp_, _tp_, _s_)
#define near__p_SPLIT_GENERIC_INIT_TPROCS_MOD(_tp_, _r_...)      { 4, NULL, near__p_SPEC_EXT_PARAM_TPROC_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define near__p_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD(_tp_, _r_...)  { 4, NULL, near__p_SPEC_EXT_PARAM_TPROC_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, near__p_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, near__p_SPEC_EXT_PARAM_TPROCS_MOD(_tp_), _r_ }

/* sl_type near__p_tloc_f near__p_tloc_mod_f */
typedef near__p_slint_t near__p_tloc_f(near__p_elements_t *b, near__p_slint_t x, void *tloc_data);
typedef near__p_slint_t near__p_tloc_mod_f(near__p_elements_t *b, near__p_slint_t x, void *tloc_data, near__p_elements_t *mod);

/* sl_type near__p_tproc_f near__p_tproc_mod_f near__p_tprocs_f near__p_tprocs_mod_f */
typedef int near__p_tproc_f(near__p_elements_t *b, near__p_slint_t x, void *tproc_data);
typedef int near__p_tproc_mod_f(near__p_elements_t *b, near__p_slint_t x, void *tproc_data, near__p_elements_t *mod);
typedef near__p_slint_t near__p_tprocs_f(near__p_elements_t *b, near__p_slint_t x, void *tproc_data, int *procs);
typedef near__p_slint_t near__p_tprocs_mod_f(near__p_elements_t *b, near__p_slint_t x, void *tproc_data, int *procs, near__p_elements_t *mod);

/* sl_type near__p_tproc_reset_f */
typedef void near__p_tproc_reset_f(void *tproc_data);

/* sl_macro near__p_TPROC_RESET_NULL */
#define near__p_TPROC_RESET_NULL  NULL

/* sl_type near__p__tproc_t near__p_tproc_t */
typedef struct near__p__tproc_t *near__p_tproc_t;

/* sl_type near__p__tproc_exdef near__p_tproc_exdef */
typedef struct near__p__tproc_exdef {
  int type;

  near__p_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  near__p_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  near__p_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  near__p_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  near__p_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  near__p_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  near__p_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  near__p_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  near__p_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  near__p_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  near__p_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  near__p_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

} const *near__p_tproc_exdef;

/* sl_macro near__p_TPROC_EXDEF_NULL */
#define near__p_TPROC_EXDEF_NULL  NULL

/* sl_macro near__p_TPROC_EXDEF_DEFINE_TPROC near__p_TPROC_EXDEF_DEFINE_TPROC_MOD near__p_TPROC_EXDEF_DEFINE_TPROCS near__p_TPROC_EXDEF_DEFINE_TPROCS_MOD */
#define near__p_TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  near__p_SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct near__p__tproc_exdef _##_name_ = { 1, near__p_SPEC_EXT_PARAM_TPROC(_name_), near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL, near__p_SPEC_EXT_PARAM_TPROCS_NULL, near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define near__p_TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  near__p_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct near__p__tproc_exdef _##_name_ = { 2, near__p_SPEC_EXT_PARAM_TPROC_NULL, near__p_SPEC_EXT_PARAM_TPROC_MOD(_name_), near__p_SPEC_EXT_PARAM_TPROCS_NULL, near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define near__p_TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  near__p_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct near__p__tproc_exdef _##_name_ = { 3, near__p_SPEC_EXT_PARAM_TPROC_NULL, near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL, near__p_SPEC_EXT_PARAM_TPROCS(_name_), near__p_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define near__p_TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  near__p_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct near__p__tproc_exdef _##_name_ = { 4, near__p_SPEC_EXT_PARAM_TPROC_NULL, near__p_SPEC_EXT_PARAM_TPROC_MOD_NULL, near__p_SPEC_EXT_PARAM_TPROCS_NULL, near__p_SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;


/* deprecated, sl_type near__p_k2c_func near__p_pivot_func near__p_sn_func near__p_m2x_func near__p_m2X_func */
typedef near__p_key2class_f near__p_k2c_func;
typedef near__p_pivot_f near__p_pivot_func;
typedef near__p_sortnet_f near__p_sn_func;
typedef near__p_merge2x_f near__p_m2x_func;
typedef near__p_merge2X_f near__p_m2X_func;


/* sl_type near__p__mergek_t near__p_mergek_t */
typedef struct near__p__mergek_t
{
  near__p_sortnet_f sn;
  near__p_sortnet_data_t snd;

  near__p_merge2x_f m2x;
  near__p_elements_t *sx;

} near__p_mergek_t;


/* sl_type near__p_keys_init_type_t near__p_keys_init_data_t */
typedef near__p_slint_t near__p_keys_init_type_t;
typedef void *near__p_keys_init_data_t;

/* sl_type near__p_key_set_data_t near__p_key_set_f */
typedef void *near__p_key_set_data_t;
typedef void (*near__p_key_set_f)(near__p_slkey_pure_t *k, near__p_key_set_data_t d);


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


/* near__p_elements_keys_stats */
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


/* partition conditions, sl_type near__p__partcond2_t near__p_partcond2_t */
typedef struct near__p__partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} near__p_partcond2_t;


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

/* partition conditions, sl_type near__p__partcond_t near__p_partcond_t near__p_partcond_p */
typedef struct near__p__partcond_t
{
  near__p_slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} near__p_partcond_t, *near__p_partcond_p;


/* internal partition conditions, sl_type near__p__partcond_intern_t near__p_partcond_intern_t near__p_partcond_intern_p */
typedef struct near__p__partcond_intern_t
{
  near__p_slint_t pcm;
  near__p_slint_t count_min, count_max;
  near__p_slint_t count_low, count_high;
#ifdef elem_weight
  near__p_slweight_t weight_min, weight_max;
  near__p_slweight_t weight_low, weight_high;
#endif

} near__p_partcond_intern_t, *near__p_partcond_intern_p;


/* sl_type near__p__parttype_t near__p_parttype_t near__p_parttype_p */
typedef struct near__p__parttype_t
{
  near__p_slint_t type;

} near__p_parttype_t, *near__p_parttype_p;


/* generic binning method */

/* sl_type near__p__bin_t near__p_bin_t */
typedef struct near__p__bin_t
{
  near__p_elements_t s;

#ifdef elem_weight
  near__p_slweight_t weight;
#endif

} near__p_bin_t;


/* sl_type near__p__splitter_t near__p_splitter_t */
typedef struct near__p__splitter_t
{
  near__p_slint_t n;

  int *displs;
  near__p_slkey_pure_t *s;
  near__p_slint_t *sn;

} near__p_splitter_t;


struct near__p__binning_t;

/* sl_type near__p_binning_pre_f near__p_binning_exec_f near__p_binning_refine_f near__p_binning_hit_f near__p_binning_finalize_f near__p_binning_post_f */
typedef near__p_slint_t (*near__p_binning_pre_f)(struct near__p__binning_t *bm);
typedef near__p_slint_t (*near__p_binning_exec_f)(struct near__p__binning_t *bm, near__p_bin_t *bin, near__p_slcount_t *counts, near__p_slweight_t *weights);
typedef near__p_slint_t (*near__p_binning_refine_f)(struct near__p__binning_t *bm, near__p_bin_t *bin, near__p_slint_t k, near__p_slcount_t *counts, near__p_slweight_t *weights, near__p_splitter_t *sp, near__p_slint_t s, near__p_bin_t *new_bin);
typedef near__p_slint_t (*near__p_binning_hit_f)(struct near__p__binning_t *bm, near__p_bin_t *bin, near__p_slint_t k, near__p_slcount_t *counts, near__p_splitter_t *sp, near__p_slint_t s);
typedef near__p_slint_t (*near__p_binning_finalize_f)(struct near__p__binning_t *bm, near__p_bin_t *bin, near__p_slint_t dc, near__p_slweight_t dw, near__p_slint_t lc_min, near__p_slint_t lc_max, near__p_slcount_t *lcs, near__p_slweight_t *lws, near__p_splitter_t *sp, near__p_slint_t s);
typedef near__p_slint_t (*near__p_binning_post_f)(struct near__p__binning_t *bm);


/* sl_type near__p__binning_data_t near__p_binning_data_t */
typedef union near__p__binning_data_t
{
  struct
  {
    near__p_slint_t rhigh, rlow, rwidth;
    near__p_slint_t rcurrent;
    near__p_slkey_pure_t bit_mask;

    near__p_elements_t sx;

  } radix;

} near__p_binning_data_t;


/* sl_type near__p__binning_t near__p_binning_t */
typedef struct near__p__binning_t
{
  near__p_slint_t nbins, max_nbins;
  
  near__p_binning_pre_f pre;
  near__p_binning_exec_f exec;
  near__p_binning_refine_f refine;
  near__p_binning_hit_f hit;
  near__p_binning_finalize_f finalize;
  near__p_binning_post_f post;

  near__p_slint_t sorted;

  near__p_slint_t docounts;
#ifdef elem_weight
  near__p_slint_t doweights;
#endif

  near__p_binning_data_t bd;

} near__p_binning_t;


/* sl_type near__p__local_bins_t near__p_local_bins_t */
typedef struct near__p__local_bins_t
{
  near__p_binning_t *bm;

  near__p_slint_t nbins, max_nbins;
  near__p_slint_t nelements;

  near__p_slint_t docounts;
#ifdef elem_weight
  near__p_slint_t doweights;
#endif

  near__p_slint_t nbinnings, max_nbinnings;

  near__p_slint_t nbins_new, last_new_b, last_new_k;
  near__p_bin_t *bins, *bins_new;
  near__p_bin_t *bins0, *bins1;

  near__p_slint_t *bcws;

#if defined(elem_weight) && defined(near__p_sl_weight_intequiv)
  near__p_slint_t cw_factor, w_index, bin_cw_factor;
  near__p_slweight_t *cws, *bin_cws;
  near__p_slweight_t *prefix_cws;
#else
  near__p_slint_t *cs, *bin_cs;
  near__p_slint_t *prefix_cs;
# ifdef elem_weight
  near__p_slweight_t *ws, *bin_ws;
  near__p_slweight_t *prefix_ws;
# endif
#endif

  near__p_slint_t last_exec_b;

} near__p_local_bins_t;


/* sl_type near__p__global_bins_t near__p_global_bins_t */
typedef struct near__p__global_bins_t
{
  near__p_binning_t *bm;
  
  near__p_local_bins_t lb;

  near__p_slint_t *bcws;

#if defined(elem_weight) && defined(near__p_sl_weight_intequiv)
  near__p_slweight_t *cws;
  near__p_slweight_t *prefix_cws;
#else
  near__p_slint_t *cs;
  near__p_slint_t *prefix_cs;
# ifdef elem_weight
  near__p_slweight_t *ws;
  near__p_slweight_t *prefix_ws;
# endif
#endif

} near__p_global_bins_t;


/* sl_type near__p_rti_cmc_t */
typedef struct
{
  near__p_slint_t cmp, movek, moved;

} near__p_rti_cmc_t;

#ifndef my_rti_ccmp
# define my_rti_ccmp(m)    m.cmc.cmp
# define my_rti_cmovek(m)  m.cmc.movek
# define my_rti_cmoved(m)  m.cmc.moved
#endif


/* sl_type near__p_rti_tim_t */
typedef struct
{
  double start, stop;
  double last, cumu;

  near__p_slint_t num;

} near__p_rti_tim_t[rti_tids];

#ifndef my_rti_tlast
# define my_rti_tlast(m, t)  m.tim[t].last
# define my_rti_tcumu(m, t)  m.tim[t].cumu
# define my_rti_tnum(m, t)   m.tim[t].num
#endif


/* sl_type near__p_rti_mem_t */
typedef struct
{
  near__p_slint_t nalloc, nfree;
  near__p_slint_t max, cur, cur_max;

} near__p_rti_mem_t;


/* sl_type near__p_rti_t */
typedef struct
{
  /* compare-move-counter */
  near__p_rti_cmc_t cmc;
  /* timer */
  near__p_rti_tim_t tim;
  /* memory */
  near__p_rti_mem_t mem;

} near__p_rti_t;

#ifndef my_rti_reset
# define my_rti_reset(m)  memset((void *) &m, 0, sizeof(m))
#endif


/* sl_type near__p__sl_context_t near__p_sl_context_t */
typedef struct near__p__sl_context_t
{

/* src/base/base.c */
  struct {
int dummy_rank;
  } sl;
#ifdef near__p_SL_USE_RTI
near__p_rti_t rti;
#endif
  struct {
near__p_slint_t ip_threshold;
near__p_slint_t db_threshold;
near__p_slint_t ma_threshold;
  } sr;
  struct {
near__p_slint_t threshold;
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
near__p_slint_t sendrecv_replace_memsize;
near__p_slint_t sendrecv_replace_mpi_maxsize;
  } me;
#endif
#ifdef SL_USE_MPI
  struct {
double t[2];
near__p_slint_t max_nprocs;
near__p_slint_t packed;
near__p_slint_t minalloc;
double overalloc;
  } meas;
#endif
#ifdef SL_USE_MPI
  struct {
near__p_slint_t packed;
near__p_slint_t db_packed;
near__p_slint_t ip_packed;
  } mea;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef near__p_MSEG_ROOT
int root;
#endif
#ifdef near__p_MSEG_BORDER_UPDATE_REDUCTION
double border_update_count_reduction;
double border_update_weight_reduction;
#endif
#ifdef near__p_MSEG_FORWARD_ONLY
near__p_slint_t forward_only;
#endif
#ifdef near__p_MSEG_INFO
near__p_slint_t info_rounds;
near__p_slint_t *info_finish_rounds;
double info_finish_rounds_avg;
near__p_slweight_t info_total_weights;
#endif
near__p_slint_t binnings;
near__p_slint_t finalize_mode;
  } mseg;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef near__p_MSS_ROOT
int root;
#endif
  } mss;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
near__p_slint_t sync;
  } msm;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
near__p_slint_t sync;
near__p_partcond_t *r_pc;
  } msp;
#endif
#ifdef SL_USE_MPI
  struct {
double i_t[3];
double p_t[3];
double b_t[3];
near__p_slint_t sync;
near__p_slint_t i_sync;
near__p_slint_t p_sync;
near__p_slint_t b_sync;
near__p_slint_t back_packed;
  } mssp;
#endif
} near__p_sl_context_t;






/* sl_macro near__p_elem_set_size near__p_elem_set_max_size near__p_elem_set_keys near__p_elem_set_indices */
#define near__p_elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define near__p_elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define near__p_elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define near__p_elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))
/* DATAX_TEMPLATE_BEGIN */
#define near__p_elem_set_data0(_e_, _b_)     ((_e_)->data0 = (_b_))  /* sl_macro */
#define near__p_elem_set_data1(_e_, _b_)     ((_e_)->data1 = (_b_))  /* sl_macro */
#define near__p_elem_set_data2(_e_, _b_)     ((_e_)->data2 = (_b_))  /* sl_macro */
#define near__p_elem_set_data3(_e_, _b_)     ((_e_)->data3 = (_b_))  /* sl_macro */
#define near__p_elem_set_data4(_e_, _b_)     ((_e_)->data4 = (_b_))  /* sl_macro */
#define near__p_elem_set_data5(_e_, _b_)     ((_e_)->data5 = (_b_))  /* sl_macro */
#define near__p_elem_set_data6(_e_, _b_)     ((_e_)->data6 = (_b_))  /* sl_macro */
#define near__p_elem_set_data7(_e_, _b_)     ((_e_)->data7 = (_b_))  /* sl_macro */
#define near__p_elem_set_data8(_e_, _b_)     ((_e_)->data8 = (_b_))  /* sl_macro */
#define near__p_elem_set_data9(_e_, _b_)     ((_e_)->data9 = (_b_))  /* sl_macro */
#define near__p_elem_set_data10(_e_, _b_)     ((_e_)->data10 = (_b_))  /* sl_macro */
#define near__p_elem_set_data11(_e_, _b_)     ((_e_)->data11 = (_b_))  /* sl_macro */
#define near__p_elem_set_data12(_e_, _b_)     ((_e_)->data12 = (_b_))  /* sl_macro */
#define near__p_elem_set_data13(_e_, _b_)     ((_e_)->data13 = (_b_))  /* sl_macro */
#define near__p_elem_set_data14(_e_, _b_)     ((_e_)->data14 = (_b_))  /* sl_macro */
#define near__p_elem_set_data15(_e_, _b_)     ((_e_)->data15 = (_b_))  /* sl_macro */
#define near__p_elem_set_data16(_e_, _b_)     ((_e_)->data16 = (_b_))  /* sl_macro */
#define near__p_elem_set_data17(_e_, _b_)     ((_e_)->data17 = (_b_))  /* sl_macro */
#define near__p_elem_set_data18(_e_, _b_)     ((_e_)->data18 = (_b_))  /* sl_macro */
#define near__p_elem_set_data19(_e_, _b_)     ((_e_)->data19 = (_b_))  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro near__p_elem_get_size near__p_elem_get_max_size near__p_elem_get_keys near__p_elem_get_indices */
#define near__p_elem_get_size(_e_)           (_e_)->size
#define near__p_elem_get_max_size(_e_)       (_e_)->max_size
#define near__p_elem_get_keys(_e_)           (_e_)->keys
#define near__p_elem_get_indices(_e_)        (_e_)->indices
/* DATAX_TEMPLATE_BEGIN */
#define near__p_elem_get_data0(_e_)          (_e_)->data0  /* sl_macro */
#define near__p_elem_get_data1(_e_)          (_e_)->data1  /* sl_macro */
#define near__p_elem_get_data2(_e_)          (_e_)->data2  /* sl_macro */
#define near__p_elem_get_data3(_e_)          (_e_)->data3  /* sl_macro */
#define near__p_elem_get_data4(_e_)          (_e_)->data4  /* sl_macro */
#define near__p_elem_get_data5(_e_)          (_e_)->data5  /* sl_macro */
#define near__p_elem_get_data6(_e_)          (_e_)->data6  /* sl_macro */
#define near__p_elem_get_data7(_e_)          (_e_)->data7  /* sl_macro */
#define near__p_elem_get_data8(_e_)          (_e_)->data8  /* sl_macro */
#define near__p_elem_get_data9(_e_)          (_e_)->data9  /* sl_macro */
#define near__p_elem_get_data10(_e_)          (_e_)->data10  /* sl_macro */
#define near__p_elem_get_data11(_e_)          (_e_)->data11  /* sl_macro */
#define near__p_elem_get_data12(_e_)          (_e_)->data12  /* sl_macro */
#define near__p_elem_get_data13(_e_)          (_e_)->data13  /* sl_macro */
#define near__p_elem_get_data14(_e_)          (_e_)->data14  /* sl_macro */
#define near__p_elem_get_data15(_e_)          (_e_)->data15  /* sl_macro */
#define near__p_elem_get_data16(_e_)          (_e_)->data16  /* sl_macro */
#define near__p_elem_get_data17(_e_)          (_e_)->data17  /* sl_macro */
#define near__p_elem_get_data18(_e_)          (_e_)->data18  /* sl_macro */
#define near__p_elem_get_data19(_e_)          (_e_)->data19  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro near__p_elem_set_block near__p_elem_set_block_size near__p_elem_get_block near__p_elem_get_block_size */
#define near__p_elem_set_block(_e_, _b_)       ((_e_)->keys = (_b_), (_e_)->max_size = -1)
#define near__p_elem_set_block_size(_e_, _s_)  ((_e_)->size = (_s_))
#define near__p_elem_get_block(_e_)            ((void *) (((_e_)->max_size < 0)?(_e_)->keys:NULL))
#define near__p_elem_get_block_size(_e_)       (_e_)->size

/* sl_macro near__p_pelem_set_size near__p_pelem_set_max_size near__p_pelem_set_elements */
#define near__p_pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define near__p_pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define near__p_pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))

/* sl_macro near__p_pelem_get_size near__p_pelem_get_max_size near__p_pelem_get_elements */
#define near__p_pelem_get_size(_e_)           (_e_)->size
#define near__p_pelem_get_max_size(_e_)       (_e_)->max_size
#define near__p_pelem_get_elements(_e_)       (_e_)->elements

/* sl_macro near__p_SL_DEFCON */
#define near__p_SL_DEFCON(_v_)  (near__p_sl_default_context._v_)






/* src/base/base.c */
extern near__p_sl_context_t near__p_sl_default_context;
extern const int near__p_default_sl_dummy_rank;
#ifdef near__p_SL_USE_RTI
extern const near__p_rti_t near__p_default_rti;
#endif
extern const near__p_slint_t near__p_default_sr_ip_threshold;
extern const near__p_slint_t near__p_default_sr_db_threshold;
extern const near__p_slint_t near__p_default_sr_ma_threshold;
extern const near__p_slint_t near__p_default_sri_threshold;

/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
extern const MPI_Datatype near__p_default_mpi_int_datatype;
extern const MPI_Datatype near__p_default_mpi_key_datatype;
extern const MPI_Datatype near__p_default_mpi_pkey_datatype;
extern const MPI_Datatype near__p_default_mpi_pwkey_datatype;
extern const MPI_Datatype near__p_default_mpi_index_datatype;
extern const MPI_Datatype near__p_default_mpi_weight_datatype;
extern const MPI_Datatype near__p_default_mpi_data_datatype[];
extern const int near__p_default_mpi_rank;
#endif
extern const void *near__p_default_me_sendrecv_replace_mem;
extern const near__p_slint_t near__p_default_me_sendrecv_replace_memsize;
extern const near__p_slint_t near__p_default_me_sendrecv_replace_mpi_maxsize;
extern const double near__p_default_meas_t[];
extern const near__p_slint_t near__p_default_meas_max_nprocs;
extern const near__p_slint_t near__p_default_meas_packed;
extern const near__p_slint_t near__p_default_meas_minalloc;
extern const double near__p_default_meas_overalloc;
extern const near__p_slint_t near__p_default_mea_packed;
extern const near__p_slint_t near__p_default_mea_db_packed;
extern const near__p_slint_t near__p_default_mea_ip_packed;
#ifdef near__p_MSEG_ROOT
extern const int near__p_default_mseg_root;
#endif
#ifdef near__p_MSEG_BORDER_UPDATE_REDUCTION
extern const double near__p_default_mseg_border_update_count_reduction;
extern const double near__p_default_mseg_border_update_weight_reduction;
#endif
#ifdef near__p_MSEG_FORWARD_ONLY
extern const near__p_slint_t near__p_default_mseg_forward_only;
#endif
#ifdef near__p_MSEG_INFO
extern const near__p_slint_t near__p_default_mseg_info_rounds;
extern const near__p_slint_t *near__p_default_mseg_info_finish_rounds;
extern const double near__p_default_mseg_info_finish_rounds_avg;
extern const near__p_slweight_t near__p_default_mseg_info_total_weights;
#endif
extern const near__p_slint_t near__p_default_mseg_binnings;
extern const near__p_slint_t near__p_default_mseg_finalize_mode;
#ifdef near__p_MSS_ROOT
extern const int near__p_default_mss_root;
#endif
extern const double near__p_default_msm_t[];
extern const near__p_slint_t near__p_default_msm_sync;
extern const double near__p_default_msp_t[];
extern const near__p_slint_t near__p_default_msp_sync;
extern const near__p_partcond_t *near__p_default_msp_r_pc;
extern const double near__p_default_mssp_i_t[];
extern const double near__p_default_mssp_p_t[];
extern const double near__p_default_mssp_b_t[];
extern const near__p_slint_t near__p_default_mssp_sync;
extern const near__p_slint_t near__p_default_mssp_i_sync;
extern const near__p_slint_t near__p_default_mssp_p_sync;
extern const near__p_slint_t near__p_default_mssp_b_sync;
extern const near__p_slint_t near__p_default_mssp_back_packed;






/* src/base/base.c */
near__p_slint_t SL_PROTO(near__p_binning_create)(near__p_local_bins_t *lb, near__p_slint_t max_nbins, near__p_slint_t max_nbinnings, near__p_elements_t *s, near__p_slint_t nelements, near__p_slint_t docounts, near__p_slint_t doweights, near__p_binning_t *bm);
near__p_slint_t SL_PROTO(near__p_binning_destroy)(near__p_local_bins_t *lb);
near__p_slint_t SL_PROTO(near__p_binning_pre)(near__p_local_bins_t *lb);
near__p_slint_t SL_PROTO(near__p_binning_exec_reset)(near__p_local_bins_t *lb, near__p_slint_t do_bins, near__p_slint_t do_prefixes);
near__p_slint_t SL_PROTO(near__p_binning_exec)(near__p_local_bins_t *lb, near__p_slint_t b, near__p_slint_t do_bins, near__p_slint_t do_prefixes);
near__p_slint_t SL_PROTO(near__p_binning_refine)(near__p_local_bins_t *lb, near__p_slint_t b, near__p_slint_t k, near__p_splitter_t *sp, near__p_slint_t s);
near__p_slint_t SL_PROTO(near__p_binning_hit)(near__p_local_bins_t *lb, near__p_slint_t b, near__p_slint_t k, near__p_splitter_t *sp, near__p_slint_t s);
near__p_slint_t SL_PROTO(near__p_binning_finalize)(near__p_local_bins_t *lb, near__p_slint_t b, near__p_slint_t dc, near__p_slweight_t dw, near__p_slint_t lc_min, near__p_slint_t lc_max, near__p_slcount_t *lcs, near__p_slweight_t *lws, near__p_splitter_t *sp, near__p_slint_t s);
near__p_slint_t SL_PROTO(near__p_binning_post)(near__p_local_bins_t *lb);
near__p_slint_t SL_PROTO(near__p_binning_radix_create)(near__p_binning_t *bm, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, near__p_slint_t sorted);
near__p_slint_t SL_PROTO(near__p_binning_radix_destroy)(near__p_binning_t *bm);
near__p_slint_t SL_PROTO(near__p_binning_radix_pre)(near__p_binning_t *bm);
near__p_slint_t SL_PROTO(near__p_binning_radix_exec)(near__p_binning_t *bm, near__p_bin_t *bin, near__p_slcount_t *counts, near__p_slweight_t *weights);
near__p_slint_t SL_PROTO(near__p_binning_radix_refine)(near__p_binning_t *bm, near__p_bin_t *bin, near__p_slint_t k, near__p_slcount_t *counts, near__p_slweight_t *weights, near__p_splitter_t *sp, near__p_slint_t s, near__p_bin_t *new_bin);
near__p_slint_t SL_PROTO(near__p_binning_radix_hit)(near__p_binning_t *bm, near__p_bin_t *bin, near__p_slint_t k, near__p_slcount_t *counts, near__p_splitter_t *sp, near__p_slint_t s);
near__p_slint_t SL_PROTO(near__p_binning_radix_finalize)(near__p_binning_t *bm, near__p_bin_t *bin, near__p_slint_t dc, near__p_slweight_t dw, near__p_slint_t lc_min, near__p_slint_t lc_max, near__p_slcount_t *lcs, near__p_slweight_t *lws, near__p_splitter_t *sp, near__p_slint_t s);
near__p_slint_t SL_PROTO(near__p_binning_radix_post)(near__p_binning_t *bm);
near__p_slint_t SL_PROTO(near__p_elements_alloc)(near__p_elements_t *s, near__p_slint_t nelements, slcint_t components);
near__p_slint_t SL_PROTO(near__p_elements_free)(near__p_elements_t *s);
near__p_slint_t SL_PROTO(near__p_elements_realloc)(near__p_elements_t *s, near__p_slint_t nelements, slcint_t components);
near__p_slint_t SL_PROTO(near__p_elements_alloca)(near__p_elements_t *s, near__p_slint_t nelements, slcint_t components);
near__p_slint_t SL_PROTO(near__p_elements_freea)(near__p_elements_t *s);
near__p_slint_t SL_PROTO(near__p_elements_alloc_from_blocks)(near__p_elements_t *s, near__p_slint_t nblocks, void **blocks, near__p_slint_t *blocksizes, near__p_slint_t alignment, near__p_slint_t nmax, slcint_t components);
near__p_slint_t SL_PROTO(near__p_elements_alloc_from_block)(near__p_elements_t *s, void *block, near__p_slint_t blocksize, near__p_slint_t alignment, near__p_slint_t nmax, slcint_t components);
near__p_slint_t SL_PROTO(near__p_elements_alloc_block)(near__p_elements_t *s, void **block, near__p_slint_t *blocksize, near__p_slint_t alignment, near__p_slint_t maxblocksize);
near__p_slint_t SL_PROTO(near__p_elements_copy)(near__p_elements_t *s, near__p_elements_t *d);
near__p_slint_t SL_PROTO(near__p_elements_copy_at)(near__p_elements_t *s, near__p_slint_t sat, near__p_elements_t *d, near__p_slint_t dat);
near__p_slint_t SL_PROTO(near__p_elements_ncopy)(near__p_elements_t *s, near__p_elements_t *d, near__p_slint_t n);
near__p_slint_t SL_PROTO(near__p_elements_nmove)(near__p_elements_t *s, near__p_elements_t *d, near__p_slint_t n);
near__p_slint_t SL_PROTO(near__p_elements_printf)(near__p_elements_t *s, const char *prefix);
near__p_slint_t SL_PROTO(near__p_elements_extract)(near__p_elements_t *src, near__p_slint_t nelements, near__p_elements_t *dst0, near__p_elements_t *dst1);
near__p_slint_t SL_PROTO(near__p_elements_touch)(near__p_elements_t *s);
near__p_slint_t SL_PROTO(near__p_elements_digest_sum)(near__p_elements_t *s, near__p_slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(near__p_elements_crc32)(near__p_elements_t *s, near__p_slint nelements, near__p_slint_t keys, near__p_slint_t data);
near__p_slint_t SL_PROTO(near__p_elements_digest_hash)(near__p_elements_t *s, near__p_slint_t nelements, slcint_t components, void *hash);
near__p_slint_t SL_PROTO(near__p_elements_random_exchange)(near__p_elements_t *s, near__p_slint_t rounds, near__p_elements_t *xs);
near__p_slint_t SL_PROTO(near__p_elements_keys_init_seed)(unsigned long s);
near__p_slint_t SL_PROTO(near__p_elements_keys_init)(near__p_elements_t *s, near__p_keys_init_type_t t, near__p_keys_init_data_t d);
near__p_slint_t SL_PROTO(near__p_elements_keys_init_randomized)(near__p_elements_t *s, near__p_slint_t nkeys, near__p_keys_init_type_t t, near__p_keys_init_data_t d);
near__p_slint_t SL_PROTO(near__p_elements_keys_init_from_file)(near__p_elements_t *s, near__p_slint_t data, char *filename, near__p_slint_t from, near__p_slint_t to, near__p_slint_t const_bytes_per_line);
near__p_slint_t SL_PROTO(near__p_elements_keys_save_to_file)(near__p_elements_t *s, char *filename);
near__p_slint_t SL_PROTO(near__p_elements_validate_order)(near__p_elements_t *s, near__p_slint_t n);
near__p_slint_t SL_PROTO(near__p_elements_validate_order_bmask)(near__p_elements_t *s, near__p_slint_t n, near__p_slkey_pure_t bmask);
near__p_slint_t SL_PROTO(near__p_elements_validate_order_weight)(near__p_elements_t *s, near__p_slint_t n, near__p_slkey_pure_t weight);
near__p_slint_t SL_PROTO(near__p_elements_keys_stats)(near__p_elements_t *s, near__p_slkey_pure_t *stats);
near__p_slint_t SL_PROTO(near__p_elements_keys_stats_print)(near__p_elements_t *s);
near__p_slint_t SL_PROTO(near__p_elements_print_keys)(near__p_elements_t *s);
near__p_slint_t SL_PROTO(near__p_elements_print_all)(near__p_elements_t *s);
near__p_slweight_t SL_PROTO(near__p_elements_get_weight)(near__p_elements_t *s);
near__p_slint_t SL_PROTO(near__p_elements_get_minmax_keys)(near__p_elements_t *s, near__p_slint_t nelements, near__p_slkey_pure_t *minmaxkeys);
near__p_slint_t SL_PROTO(near__p_elements_alloc_packed)(near__p_packed_elements_t *s, near__p_slint_t nelements);
near__p_slint_t SL_PROTO(near__p_elements_free_packed)(near__p_packed_elements_t *s);
near__p_slint_t SL_PROTO(near__p_elements_alloc_packed_from_block)(near__p_packed_elements_t *s, void *block, near__p_slint_t blocksize, near__p_slint_t alignment, near__p_slint_t nmax);
near__p_slint_t SL_PROTO(near__p_elements_pack_indexed)(near__p_elements_t *s, near__p_packed_elements_t *d, near__p_slindex_t *rindx, near__p_slindex_t *windx);
near__p_slint_t SL_PROTO(near__p_elements_pack)(near__p_elements_t *s, near__p_packed_elements_t *d);
near__p_slint_t SL_PROTO(near__p_elements_pack_at)(near__p_elements_t *s, near__p_slint_t sat, near__p_packed_elements_t *d, near__p_slint_t dat);
near__p_slint_t SL_PROTO(near__p_elements_unpack_indexed)(near__p_packed_elements_t *s, near__p_elements_t *d, near__p_slindex_t *rindx, near__p_slindex_t *windx);
near__p_slint_t SL_PROTO(near__p_elements_unpack)(near__p_packed_elements_t *s, near__p_elements_t *d);
near__p_slint_t SL_PROTO(near__p_elements_unpack_at)(near__p_packed_elements_t *s, near__p_slint_t sat, near__p_elements_t *d, near__p_slint_t dat);
near__p_slint_t SL_PROTO(near__p_elements_unpack_keys)(near__p_packed_elements_t *s, near__p_slkey_t *k);
near__p_slint SL_PROTO(near__p_merge2_basic_auto_01_x)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint SL_PROTO(near__p_merge2_basic_01_x)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx, near__p_m2x_func _x0_1, near__p_m2x_func _0x_1);
near__p_slint SL_PROTO(near__p_merge2_basic_01_X)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx, near__p_elements_t *t, near__p_m2X_func _X0_1, near__p_m2X_func _0X_1);
near__p_slint SL_PROTO(near__p_merge2_simplify_s1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx, near__p_slint s1elements);
near__p_slint SL_PROTO(near__p_merge2_memory_adaptive)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_compo_hula)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *xs);
near__p_slint_t SL_PROTO(near__p_merge2_basic_sseq_x0_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_sseq_0x_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_sseq_01_x)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_sseq_01)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_merge2_basic_sbin_x0_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_sbin_0x_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_sbin_01_x)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_sbin_01)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_merge2_basic_shyb_x0_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_shyb_0x_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_shyb_01_x)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_shyb_01)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_merge2_basic_straight_x0_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_straight_0x_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_straight_01_x)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_straight_x_0_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_merge2_basic_straight_X0_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_merge2_basic_straight_0X_1)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_merge2_basic_straight_01_X)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_merge2_basic_straight_X0_1u)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_merge2_compo_tridgell)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_mergep_2way_ip_int)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t p, int *displs, near__p_merge2x_f m2x);
near__p_slint_t SL_PROTO(near__p_mergep_2way_ip_int_rec)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t p, int *displs, near__p_merge2x_f m2x);
near__p_slint_t SL_PROTO(near__p_mergep_heap_int)(near__p_elements_t *s, near__p_elements_t *d, near__p_slint_t p, int *displs, int *counts);
near__p_slint_t SL_PROTO(near__p_mergep_heap_int_idx)(near__p_elements_t *s, near__p_elements_t *d, near__p_slint_t p, int *displs, int *counts);
near__p_slint_t SL_PROTO(near__p_mergep_heap_idx)(near__p_elements_t *s, near__p_elements_t *d, near__p_slint_t p, near__p_slindex_t *displs, near__p_slindex_t *counts);
near__p_slint_t SL_PROTO(near__p_mergep_heap_unpack_idx)(near__p_packed_elements_t *s, near__p_elements_t *d, near__p_slint_t p, near__p_slindex_t *displs, near__p_slindex_t *counts);
near__p_slint_t SL_PROTO(near__p_mergep_heap_unpack_idxonly)(near__p_packed_elements_t *s, near__p_elements_t *d, near__p_slint_t p, near__p_slindex_t *displs, near__p_slindex_t *counts);
near__p_slint_t SL_PROTO(near__p_permute_generic_db)(near__p_elements_t *s, near__p_elements_t *d, near__p_permute_generic_t *pg, void *pg_data);
near__p_slint_t SL_PROTO(near__p_permute_generic_ip)(near__p_elements_t *s, near__p_elements_t *x, near__p_permute_generic_t *pg, void *pg_data);
near__p_slint SL_PROTO(near__p_sl_search_sequential_lt)(near__p_elements_t *s, near__p_slpkey_t k);
near__p_slint SL_PROTO(near__p_sl_search_sequential_le)(near__p_elements_t *s, near__p_slpkey_t k);
near__p_slint SL_PROTO(near__p_sl_search_sequential_gt)(near__p_elements_t *s, near__p_slpkey_t k);
near__p_slint SL_PROTO(near__p_sl_search_sequential_ge)(near__p_elements_t *s, near__p_slpkey_t k);
near__p_slint SL_PROTO(near__p_sl_search_p_sequential_lt)(near__p_elements_t *s, near__p_slpkey_t *k);
near__p_slint SL_PROTO(near__p_sl_search_p_sequential_le)(near__p_elements_t *s, near__p_slpkey_t *k);
near__p_slint SL_PROTO(near__p_sl_search_p_sequential_gt)(near__p_elements_t *s, near__p_slpkey_t *k);
near__p_slint SL_PROTO(near__p_sl_search_p_sequential_ge)(near__p_elements_t *s, near__p_slpkey_t *k);
near__p_slint SL_PROTO(near__p_sl_search_binary_lt)(near__p_elements_t *s, near__p_slpkey_t k);
near__p_slint SL_PROTO(near__p_sl_search_binary_le)(near__p_elements_t *s, near__p_slpkey_t k);
near__p_slint SL_PROTO(near__p_sl_search_binary_gt)(near__p_elements_t *s, near__p_slpkey_t k);
near__p_slint SL_PROTO(near__p_sl_search_binary_ge)(near__p_elements_t *s, near__p_slpkey_t k);
near__p_slint SL_PROTO(near__p_sl_search_p_binary_lt)(near__p_elements_t *s, near__p_slpkey_t *k);
near__p_slint SL_PROTO(near__p_sl_search_p_binary_le)(near__p_elements_t *s, near__p_slpkey_t *k);
near__p_slint SL_PROTO(near__p_sl_search_p_binary_gt)(near__p_elements_t *s, near__p_slpkey_t *k);
near__p_slint SL_PROTO(near__p_sl_search_p_binary_ge)(near__p_elements_t *s, near__p_slpkey_t *k);
near__p_slint_t SL_PROTO(near__p_sl_search_binary_lt_bmask)(near__p_elements_t *s, near__p_slpkey_t k, near__p_slpkey_t bmask);
near__p_slint_t SL_PROTO(near__p_sl_search_binary_le_bmask)(near__p_elements_t *s, near__p_slpkey_t k, near__p_slpkey_t bmask);
near__p_slint_t SL_PROTO(near__p_sl_search_binary_sign_switch)(near__p_elements_t *s);
near__p_slint SL_PROTO(near__p_sl_search_hybrid_lt)(near__p_elements_t *s, near__p_slpkey_t k, near__p_slint t);
near__p_slint SL_PROTO(near__p_sl_search_hybrid_le)(near__p_elements_t *s, near__p_slpkey_t k, near__p_slint t);
near__p_slint SL_PROTO(near__p_sl_search_hybrid_gt)(near__p_elements_t *s, near__p_slpkey_t k, near__p_slint t);
near__p_slint SL_PROTO(near__p_sl_search_hybrid_ge)(near__p_elements_t *s, near__p_slpkey_t k, near__p_slint t);
near__p_slint SL_PROTO(near__p_sl_search_p_hybrid_lt)(near__p_elements_t *s, near__p_slpkey_t *k, near__p_slint t);
near__p_slint SL_PROTO(near__p_sl_search_p_hybrid_le)(near__p_elements_t *s, near__p_slpkey_t *k, near__p_slint t);
near__p_slint SL_PROTO(near__p_sl_search_p_hybrid_gt)(near__p_elements_t *s, near__p_slpkey_t *k, near__p_slint t);
near__p_slint SL_PROTO(near__p_sl_search_p_hybrid_ge)(near__p_elements_t *s, near__p_slpkey_t *k, near__p_slint t);
near__p_slint SL_PROTO(near__p_ilog2c)(near__p_slint x);
near__p_slint SL_PROTO(near__p_ilog2f)(near__p_slint x);
near__p_slint SL_PROTO(near__p_print_bits)(near__p_slint v);
near__p_slint SL_PROTO(near__p_pivot_random)(near__p_elements_t *s);
near__p_slint_t SL_PROTO(near__p_counts2displs)(near__p_slint_t n, int *counts, int *displs);
near__p_slint_t SL_PROTO(near__p_displs2counts)(near__p_slint_t n, int *displs, int *counts, near__p_slint_t total_counts);
void SL_PROTO(near__p_get_displcounts_extent)(near__p_slint_t n, int *displs, int *counts, near__p_slint_t *lb, near__p_slint_t *extent);
void SL_PROTO(near__p_elem_set_data)(near__p_elements_t *e, ...);
near__p_slint_t SL_PROTO(near__p_elem_get_max_byte)();
near__p_slint_t SL_PROTO(near__p_elem_reverse)(near__p_elements_t *e, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_elem_nxchange_at)(near__p_elements_t *e0, near__p_slint_t at0, near__p_elements_t *e1, near__p_slint_t at1, near__p_slint_t n, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_elem_nxchange)(near__p_elements_t *e0, near__p_elements_t *e1, near__p_slint_t n, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_elem_nxchange_ro0)(near__p_elements_t *e0, near__p_elements_t *e1, near__p_slint_t n, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_elem_rotate)(near__p_elements_t *e, near__p_slint_t m, near__p_slint_t n, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_elem_rotate_ro0)(near__p_elements_t *e, near__p_slint_t m, near__p_slint_t n, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_elem_rotate_ro1)(near__p_elements_t *e, near__p_slint_t m, near__p_slint_t n, near__p_elements_t *t);
near__p_slint_t SL_PROTO(near__p_sort_counting_use_displs)(near__p_elements_t *s, near__p_elements_t *d, near__p_slint_t ndispls, near__p_slint_t *displs);
near__p_slint_t SL_PROTO(near__p_sort_counting_use_counts)(near__p_elements_t *s, near__p_elements_t *d, near__p_slint_t ncounts, near__p_slint_t *counts);
near__p_slint_t SL_PROTO(near__p_sort_counting_get_counts)(near__p_elements_t *s, near__p_elements_t *d, near__p_slint_t ncounts, near__p_slint_t *counts);
near__p_slint_t SL_PROTO(near__p_sort_counting)(near__p_elements_t *s, near__p_elements_t *d, near__p_slint_t ncounts);
near__p_slint SL_PROTO(near__p_sort_heap)(near__p_elements_t *s, near__p_elements_t *xs);
near__p_slint_t SL_PROTO(near__p_sort_insert_bmask_kernel)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slkey_pure_t bmask);
near__p_slint_t SL_PROTO(near__p_sort_insert)(near__p_elements_t *s, near__p_elements_t *sx);
near__p_slint_t SL_PROTO(near__p_sort_permute_forward)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t *perm, near__p_slint_t offset, near__p_slint_t mask_bit);
near__p_slint_t SL_PROTO(near__p_sort_permute_backward)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t *perm, near__p_slint_t offset, near__p_slint_t mask_bit);
near__p_slint SL_PROTO(near__p_sort_quick)(near__p_elements_t *s, near__p_elements_t *xs);
near__p_slint_t SL_PROTO(near__p_sort_radix_ip)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth);
near__p_slint_t SL_PROTO(near__p_sort_radix_db)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth);
near__p_slint_t SL_PROTO(near__p_sort_radix_ma)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth);
near__p_slint_t SL_PROTO(near__p_sort_radix)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth);
near__p_slint_t SL_PROTO(near__p_sort_radix_1bit_kernel)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t rhigh, near__p_slint_t rlow);
near__p_slint SL_PROTO(near__p_sort_radix_1bit)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t rhigh, near__p_slint_t rlow);
near__p_slint_t SL_PROTO(near__p_sort_radix_iter)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t presorted, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth);
near__p_slint SL_PROTO(near__p_sn_hypercube_lh)(near__p_slint size, near__p_slint rank, near__p_slint stage, void *snp, near__p_slint *up);
near__p_slint SL_PROTO(near__p_sn_hypercube_hl)(near__p_slint size, near__p_slint rank, near__p_slint stage, void *snp, near__p_slint *up);
near__p_slint SL_PROTO(near__p_sn_odd_even_trans)(near__p_slint size, near__p_slint rank, near__p_slint stage, void *snp, near__p_slint *up);
near__p_slint SL_PROTO(near__p_sn_odd)(near__p_slint size, near__p_slint rank, near__p_slint stage, void *snp, near__p_slint *up);
near__p_slint SL_PROTO(near__p_sn_even)(near__p_slint size, near__p_slint rank, near__p_slint stage, void *snp, near__p_slint *up);
near__p_slint SL_PROTO(near__p_sn_batcher)(near__p_slint size, near__p_slint rank, near__p_slint stage, void *snp, near__p_slint *up);
near__p_slint SL_PROTO(near__p_sn_bitonic)(near__p_slint size, near__p_slint rank, near__p_slint stage, void *snp, near__p_slint *up);
near__p_slint SL_PROTO(near__p_sn_connected)(near__p_slint size, near__p_slint rank, near__p_slint stage, void *snp, near__p_slint *up);
near__p_slint_t SL_PROTO(near__p_split_generic_db)(near__p_elements_t *s, near__p_elements_t *d, near__p_split_generic_t *sg, void *sg_data, near__p_slint_t n);
near__p_slint_t SL_PROTO(near__p_split_generic_ip)(near__p_elements_t *s, near__p_elements_t *d, near__p_split_generic_t *sg, void *sg_data, near__p_slint_t n);
near__p_slint_t SL_PROTO(near__p_split_generic_count_db)(near__p_elements_t *s, near__p_split_generic_t *sg, void *sg_data, int *counts, near__p_slint_t n);
near__p_slint_t SL_PROTO(near__p_split_generic_count_ip)(near__p_elements_t *s, near__p_split_generic_t *sg, void *sg_data, int *counts, near__p_slint_t n);
near__p_slint_t SL_PROTO(near__p_split_generic_rearrange_db)(near__p_elements_t *s, near__p_elements_t *d, near__p_split_generic_t *sg, void *sg_data, int *counts, near__p_slint_t n);
near__p_slint_t SL_PROTO(near__p_split_generic_rearrange_ip)(near__p_elements_t *s, near__p_elements_t *d, near__p_split_generic_t *sg, void *sg_data, int *counts, int *displs, near__p_slint_t n);
near__p_slint_t SL_PROTO(near__p_splitter_reset)(near__p_splitter_t *sp);
near__p_slint_t SL_PROTO(near__p_splitx_radix)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint_t nclasses, near__p_slint_t shl, near__p_slint_t *counts);
near__p_slint SL_PROTO(near__p_split2_lt_ge)(near__p_elements_t *s, near__p_slkey_pure_t *k, near__p_elements_t *t);
near__p_slint SL_PROTO(near__p_split2_le_gt)(near__p_elements_t *s, near__p_slkey_pure_t *k, near__p_elements_t *t);
near__p_slint SL_PROTO(near__p_split3_lt_eq_gt)(near__p_elements_t *s, near__p_slkey_pure_t *k, near__p_elements_t *t, near__p_slint *nlt, near__p_slint *nle);
near__p_slint SL_PROTO(near__p_split3_lt_eq_gt_old)(near__p_elements_t *s, near__p_slkey_pure_t *k, near__p_elements_t *t, near__p_slint *nlt, near__p_slint *nle);
near__p_slint SL_PROTO(near__p_split2_b)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slkey_pure_t bmask);
near__p_slint SL_PROTO(near__p_splitk_k2c_af)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint k, near__p_slint *c, near__p_k2c_func k2c, void *k2c_data);
near__p_slint SL_PROTO(near__p_splitk_k2c)(near__p_elements_t *s, near__p_elements_t *sx, near__p_slint k, near__p_slint *c, near__p_k2c_func k2c, void *k2c_data);
near__p_slint SL_PROTO(near__p_splitk_k2c_count)(near__p_elements_t *s, near__p_slint k, near__p_slint *c, near__p_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI





/* src/base_mpi/base_mpi.c */
near__p_slint_t SL_PROTO(near__p_mpi_binning_create)(near__p_global_bins_t *gb, near__p_slint_t max_nbins, near__p_slint_t max_nbinnings, near__p_elements_t *s, near__p_slint_t nelements, near__p_slint_t docounts, near__p_slint_t doweights, near__p_binning_t *bm, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_binning_destroy)(near__p_global_bins_t *gb, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_binning_pre)(near__p_global_bins_t *gb, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_binning_exec_reset)(near__p_global_bins_t *gb, near__p_slint_t do_bins, near__p_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_binning_exec_local)(near__p_global_bins_t *gb, near__p_slint_t b, near__p_slint_t do_bins, near__p_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_binning_exec_global)(near__p_global_bins_t *gb, near__p_slint_t do_bins, near__p_slint_t do_prefixes, near__p_slint_t root, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_binning_refine)(near__p_global_bins_t *gb, near__p_slint_t b, near__p_slint_t k, near__p_splitter_t *sp, near__p_slint_t s, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_binning_hit)(near__p_global_bins_t *gb, near__p_slint_t b, near__p_slint_t k, near__p_splitter_t *sp, near__p_slint_t s, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_binning_finalize)(near__p_global_bins_t *gb, near__p_slint_t b, near__p_slint_t dc, near__p_slweight_t dw, near__p_slint_t lc_min, near__p_slint_t lc_max, near__p_slcount_t *lcs, near__p_slweight_t *lws, near__p_splitter_t *sp, near__p_slint_t s, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_binning_post)(near__p_global_bins_t *gb, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_datatypes_init)();
near__p_slint_t SL_PROTO(near__p_mpi_datatypes_release)();
near__p_slint_t SL_PROTO(near__p_mpi_get_grid_properties)(near__p_slint_t ndims, near__p_slint_t *dims, near__p_slint_t *pos, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_subgroups_create)(near__p_slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_subgroups_delete)(near__p_slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);
int SL_PROTO(near__p_sl_MPI_Allreduce)(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank);
int SL_PROTO(near__p_sl_MPI_Alltoall_int)(void *sendbuf, int sendcount, void *recvbuf, int recvcount, MPI_Comm comm, int size, int rank);
near__p_slint_t SL_PROTO(near__p_mpi_elements_keys_init_from_file)(near__p_elements_t *s, char *filename, near__p_slint from, near__p_slint to, near__p_slint const_bytes_per_line, near__p_slint root, int size, int rank, MPI_Comm comm);
near__p_slint SL_PROTO(near__p_mpi_elements_validate_order)(near__p_elements_t *s, near__p_slint n, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_linear_exchange_pure_keys)(near__p_slkey_pure_t *in, near__p_slkey_pure_t *out, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_check_order)(near__p_elements_t *s, near__p_slint_t nelements, near__p_slint_t *orders, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_check_global_order)(near__p_slkey_pure_t local_min, near__p_slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_digest_sum)(near__p_elements_t *s, near__p_slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(near__p_mpi_elements_crc32)(near__p_elements_t *s, near__p_slint_t n, near__p_slint_t keys, near__p_slint_t data, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_digest_hash)(near__p_elements_t *s, near__p_slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_get_counts)(near__p_elements_t *s, near__p_slint_t *clocal, near__p_slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
near__p_slweight_t SL_PROTO(near__p_mpi_elements_get_weights)(near__p_elements_t *s, near__p_slweight_t *wlocal, near__p_slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_get_counts_and_weights)(near__p_elements_t *s, near__p_slint_t nelements, near__p_slint_t *counts, near__p_slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_sendrecv_replace)(near__p_elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_tproc_create_tproc)(near__p_tproc_t *tproc, near__p_tproc_f *tfn, near__p_tproc_reset_f *rfn, near__p_tproc_exdef exdef);
near__p_slint_t SL_PROTO(near__p_tproc_create_tproc_mod)(near__p_tproc_t *tproc, near__p_tproc_mod_f *tfn, near__p_tproc_reset_f *rfn, near__p_tproc_exdef exdef);
near__p_slint_t SL_PROTO(near__p_tproc_create_tprocs)(near__p_tproc_t *tproc, near__p_tprocs_f *tfn, near__p_tproc_reset_f *rfn, near__p_tproc_exdef exdef);
near__p_slint_t SL_PROTO(near__p_tproc_create_tprocs_mod)(near__p_tproc_t *tproc, near__p_tprocs_mod_f *tfn, near__p_tproc_reset_f *rfn, near__p_tproc_exdef exdef);
near__p_slint_t SL_PROTO(near__p_tproc_free)(near__p_tproc_t *tproc);
near__p_slint_t SL_PROTO(near__p_tproc_set_proclists)(near__p_tproc_t *tproc, near__p_slint_t nsend_procs, int *send_procs, near__p_slint_t nrecv_procs, int *recv_procs, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_tproc_verify)(near__p_tproc_t tproc, void *data, near__p_elements_t *s, int proc);
near__p_slint_t SL_PROTO(near__p_mpi_elements_alltoall_specific)(near__p_elements_t *sin, near__p_elements_t *sout, near__p_elements_t *xs, near__p_tproc_t tproc, void *data, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_alltoallv_db_packed)(near__p_elements_t *sbuf, int *scounts, int *sdispls, near__p_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_alltoallv_db)(near__p_elements_t *sbuf, int *scounts, int *sdispls, near__p_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_alltoallv_ip_packed)(near__p_elements_t *s, near__p_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_alltoallv_ip_double)(near__p_elements_t *s, near__p_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_alltoallv_ip_mpi)(near__p_elements_t *s, near__p_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_alltoallv_ip_dash)(near__p_elements_t *s, near__p_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_alltoallv_ip)(near__p_elements_t *s, near__p_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, near__p_slint_t structured);
near__p_slint_t SL_PROTO(near__p_mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);
near__p_slint_t SL_PROTO(near__p_mpi_find_exact_equal)(near__p_elements_t *s, near__p_slint_t other_rank, near__p_slint_t high_rank, near__p_slint_t *ex_start, near__p_slint_t *ex_size, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_find_exact)(near__p_elements_t *s, near__p_slint_t other_rank, near__p_slint_t high_rank, near__p_slint_t *dst_size, near__p_slint_t *ex_start, near__p_slint_t *ex_sizes, near__p_slint_t *nx_move, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_merge2)(near__p_elements_t *s, near__p_slint_t other_rank, near__p_slint_t high_rank, near__p_slint_t *dst_size, near__p_merge2x_f m2, near__p_elements_t *xs, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_mergek_equal)(near__p_elements_t *s, near__p_sortnet_f sn, near__p_sortnet_data_t snd, near__p_merge2x_f m2x, near__p_elements_t *xs, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_mergek_sorted)(near__p_elements_t *s, near__p_merge2x_f m2x, near__p_elements_t *xs, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_mergek)(near__p_elements_t *s, near__p_sortnet_f sn, near__p_sortnet_data_t snd, near__p_merge2x_f m2x, near__p_elements_t *xs, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_mergek_equal2)(near__p_elements_t *s, near__p_sortnet_f sn, near__p_sortnet_data_t snd, near__p_merge2x_f m2x, near__p_elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);
near__p_slint_t SL_PROTO(near__p_mpi_partition_exact_generic)(near__p_elements_t *s, near__p_partcond_t *pcond, near__p_binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_partition_exact_radix)(near__p_elements_t *s, near__p_partcond_t *pcond, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, near__p_slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_partition_exact_radix_ngroups)(near__p_elements_t *s, near__p_partcond_t *pcond, near__p_slint_t ngroups, MPI_Comm *group_comms, near__p_elements_t *sx, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_partition_exact_radix_2groups)(near__p_elements_t *s, near__p_partcond_t *pcond, MPI_Comm group_comm, near__p_elements_t *sx, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_partition_sample_regular)(near__p_elements_t *s, near__p_partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_rebalance)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_slint_t stable, near__p_slint_t *dst_size, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_rebalance_alltoallv)(near__p_elements_t *sbuf, int *scounts, int *sdispls, near__p_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
void SL_PROTO(near__p_mpi_partcond_set_even)(near__p_partcond_t *pcond, near__p_slint_t pcm, near__p_slint_t ntotal, double nimba, double wtotal, double wimba, int size, int rank);
near__p_slint_t SL_PROTO(near__p_init_partconds)(near__p_slint_t npconds, near__p_partcond_t *pconds, near__p_slint_t nparts, near__p_slint_t total_count, near__p_slweight_t total_weight);
near__p_slint_t SL_PROTO(near__p_init_partconds_intern)(near__p_slint_t npconds, near__p_partcond_intern_t *pci, near__p_partcond_t *pc, near__p_slint_t nparts, near__p_slint_t total_count, near__p_slweight_t total_weight);
near__p_slint_t SL_PROTO(near__p_merge_partconds)(near__p_partcond_t *pconds_in, near__p_slint_t npconds_in, near__p_partcond_t *pcond_out);
near__p_slint_t SL_PROTO(near__p_mpi_gather_partconds_grouped)(near__p_partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, near__p_partcond_t *pconds_out, near__p_slint_t *npconds_out, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_gather_partconds)(near__p_partcond_t *pcond_in, near__p_partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_allgather_partconds)(near__p_partcond_t *pcond_in, near__p_partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_bcast_partconds)(near__p_slint_t npconds, near__p_partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_post_check_partconds)(near__p_elements_t *s, near__p_slint_t nelements, near__p_slint_t nparts, near__p_partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_post_check_partconds_intern)(near__p_elements_t *s, near__p_slint_t nelements, near__p_slint_t nparts, near__p_partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_select_stats)(near__p_elements_t *s, near__p_slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_select_exact_generic_bulk)(near__p_elements_t *s, near__p_slint_t nelements, near__p_slint_t nparts, near__p_partcond_t *pconds, near__p_binning_t *bm, near__p_splitter_t *sp, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_select_exact_generic_grouped)(near__p_elements_t *s, near__p_slint_t nelements, near__p_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, near__p_binning_t *bm, near__p_splitter_t *sp, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_select_exact_generic)(near__p_elements_t *s, near__p_slint_t nelements, near__p_slint_t nparts, near__p_partcond_t *pconds, near__p_binning_t *bm, near__p_splitter_t *sp, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_select_exact_radix)(near__p_elements_t *s, near__p_slint_t nelements, near__p_slint_t nparts, near__p_partcond_t *pconds, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, near__p_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_select_exact_radix_grouped)(near__p_elements_t *s, near__p_slint_t nelements, near__p_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, near__p_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_select_sample_regular)(near__p_elements_t *s, near__p_slint_t nparts, near__p_partcond_t *pconds, near__p_slint_t nsamples, near__p_splitter_t *sp, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_merge)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *xs, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_merge2)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *xs, near__p_slint_t merge_type, near__p_slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_merge_radix)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *xs, near__p_slint_t merge_type, near__p_slint_t sort_type, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_partition)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *xs, near__p_slint_t part_type, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_partition_radix)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *xs, near__p_slint_t part_type, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_partition_exact_radix)(near__p_elements_t *s, near__p_elements_t *sx, near__p_partcond_t *pcond, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_partition_exact_radix_ngroups)(near__p_elements_t *s, near__p_elements_t *sx, near__p_partcond_t *pcond, near__p_slint_t ngroups, MPI_Comm *group_comms, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_partition_exact_radix_2groups)(near__p_elements_t *s, near__p_elements_t *sx, near__p_partcond_t *pcond, MPI_Comm group_comm, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_insert_radix)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *xs, near__p_slpkey_t *mmkeys, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_presorted_radix)(near__p_elements_t *s0, near__p_elements_t *s1, near__p_elements_t *xs, near__p_slint_t merge_type, near__p_slint_t rhigh, near__p_slint_t rlow, near__p_slint_t rwidth, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_sort_back)(near__p_elements_t *sin, near__p_elements_t *sout, near__p_elements_t *sx, near__p_slpkey_t *lh, near__p_slint_t ntotal, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, near__p_slint_t ytotal, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_xcounts2ycounts_grouped)(int *xcounts, near__p_slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
near__p_slint_t SL_PROTO(near__p_mpi_subxdispls2ycounts)(near__p_slint_t nsubs, int *sub_xdispls, near__p_slint_t *sub_sources, near__p_slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#undef SL_PROTO
#endif /* __SL_NEAR__P_H__ */
