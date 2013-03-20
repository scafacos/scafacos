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


#ifndef __SL_NEAR____H__
#define __SL_NEAR____H__

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif /* SL_USE_MPI */

#define SL_PROTO(_f_)  _f_

#define near____sl_key_type_c        long long
#define near____sl_key_type_mpi      MPI_LONG_LONG
#define near____sl_key_size_mpi      1
#define near____sl_key_type_fmt      "lld"

#define near____sl_key_integer

#define near____SL_DATA0             /* positions */
#define near____sl_data0_type_c      fcs_float
#define near____sl_data0_size_c      3
#define near____sl_data0_type_mpi    FCS_MPI_FLOAT
#define near____sl_data0_size_mpi    3

#define near____SL_DATA1             /* charges */
#define near____sl_data1_type_c      fcs_float
#define near____sl_data1_size_c      1
#define near____sl_data1_type_mpi    FCS_MPI_FLOAT
#define near____sl_data1_size_mpi    1

#define near____SL_DATA2             /* indices */
#define near____sl_data2_type_c      long long
#define near____sl_data2_size_c      1
#define near____sl_data2_type_mpi    MPI_LONG_LONG
#define near____sl_data2_size_mpi    1

#undef near____SL_DATA3              /* field */
#define near____sl_data3_type_c      fcs_float
#define near____sl_data3_size_c      3
#define near____sl_data3_type_mpi    FCS_MPI_FLOAT
#define near____sl_data3_size_mpi    3

#undef near____SL_DATA4              /* potentials */
#define near____sl_data4_type_c      fcs_float
#define near____sl_data4_size_c      1
#define near____sl_data4_type_mpi    FCS_MPI_FLOAT
#define near____sl_data4_size_mpi    1




#if defined(MSEG_ROOT) && !defined(near____MSEG_ROOT)
# define near____MSEG_ROOT  MSEG_ROOT
#endif

#if defined(MSEG_BORDER_UPDATE_REDUCTION) && !defined(near____MSEG_BORDER_UPDATE_REDUCTION)
# define near____MSEG_BORDER_UPDATE_REDUCTION  MSEG_BORDER_UPDATE_REDUCTION
#endif

#if defined(MSEG_DISABLE_BEST_CHOICE) && !defined(near____MSEG_DISABLE_BEST_CHOICE)
# define near____MSEG_DISABLE_BEST_CHOICE  MSEG_DISABLE_BEST_CHOICE
#endif

#if defined(MSEG_DISABLE_MINMAX) && !defined(near____MSEG_DISABLE_MINMAX)
# define near____MSEG_DISABLE_MINMAX  MSEG_DISABLE_MINMAX
#endif

#if defined(MSEG_ENABLE_OPTIMZED_LOWHIGH) && !defined(near____MSEG_ENABLE_OPTIMZED_LOWHIGH)
# define near____MSEG_ENABLE_OPTIMZED_LOWHIGH  MSEG_ENABLE_OPTIMZED_LOWHIGH
#endif

#if defined(MSEG_FORWARD_ONLY) && !defined(near____MSEG_FORWARD_ONLY)
# define near____MSEG_FORWARD_ONLY  MSEG_FORWARD_ONLY
#endif

#if defined(MSEG_INFO) && !defined(near____MSEG_INFO)
# define near____MSEG_INFO  MSEG_INFO
#endif

#if defined(MSEG_TRACE_IF) && !defined(near____MSEG_TRACE_IF)
# define near____MSEG_TRACE_IF  MSEG_TRACE_IF
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


#ifndef near____SL_INDEX
# undef near____SL_PACKED_INDEX
#endif


/* if no special datatype for (sl default) integer ... */
#ifndef near____sl_int_type_c
  /* ... use a default one */
# define near____sl_int_type_c               long      /* sl_macro */
# undef near____sl_int_type_mpi
# define near____sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef near____sl_int_size_mpi
# define near____sl_int_size_mpi             1         /* sl_macro */
# undef near____sl_int_type_fmt
# define near____sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(near____sl_int_type_mpi) || !defined(near____sl_int_size_mpi)
#   error "near____sl_int_type_mpi and/or near____sl_int_size_mpi missing"
#  endif
# endif
# ifndef near____sl_int_type_fmt
#  error "near____sl_int_type_fmt macro is missing, using d as default"
#  define near____sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for (intern) weight ... */
#ifndef near____sl_weight_type_c
 /* ... use (sl default) integer */
# define near____sl_weight_type_c             near____sl_int_type_c    /* sl_macro */
# undef near____sl_weight_type_mpi
# define near____sl_weight_type_mpi           near____sl_int_type_mpi  /* sl_macro */
# undef near____sl_weight_size_mpi
# define near____sl_weight_size_mpi           near____sl_int_size_mpi  /* sl_macro */
# undef near____sl_weight_type_fmt
# define near____sl_weight_type_fmt           near____sl_int_type_fmt  /* sl_macro */
# undef near____sl_weight_intequiv
# define near____sl_weight_intequiv                            /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(near____sl_weight_type_mpi) || !defined(near____sl_weight_size_mpi)
#   error "near____sl_weight_type_mpi and/or near____sl_weight_size_mpi missing"
#  endif
# endif
# ifndef near____sl_weight_type_fmt
#  error "near____sl_weight_type_fmt macro is missing, using f as default"
#  define near____sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef near____sl_index_type_c
 /* ... use the primary integer type */
# define near____sl_index_type_c             near____sl_int_type_c
# undef near____sl_index_type_mpi
# define near____sl_index_type_mpi           near____sl_int_type_mpi
# undef near____sl_index_size_mpi
# define near____sl_index_size_mpi           near____sl_int_size_mpi
# undef near____sl_index_type_fmt
# define near____sl_index_type_fmt           near____sl_int_type_fmt
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(near____sl_index_type_mpi) || !defined(near____sl_index_size_mpi)
#   error "near____sl_index_type_mpi and/or near____sl_index_size_mpi missing"
#  endif
# endif
# ifndef near____sl_index_type_fmt
#  error "near____sl_index_type_fmt macro is missing, using d as default"
#  define near____sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef near____sl_key_pure_type_c
# define near____sl_key_pure_type_c          near____sl_key_type_c  /* sl_macro */
#endif
#ifndef near____sl_key_pure_type_mpi
# define near____sl_key_pure_type_mpi        near____sl_key_type_mpi  /* sl_macro */
#endif
#ifndef near____sl_key_pure_size_mpi
# define near____sl_key_pure_size_mpi        near____sl_key_size_mpi  /* sl_macro */
#endif
#ifndef near____sl_key_pure_type_fmt
# ifdef near____sl_key_type_fmt
#  define near____sl_key_pure_type_fmt       near____sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef near____sl_key_purify
 /* key val -> key val */
 #define near____sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef near____sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define near____sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef near____sl_key_set_pure
 /* key component pointer and key val */
 #define near____sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef near____sl_key_pure_cmp_eq
 #define near____sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef near____sl_key_pure_cmp_ne
 #define near____sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef near____sl_key_pure_cmp_lt
 #define near____sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef near____sl_key_pure_cmp_le
 #define near____sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef near____sl_key_pure_cmp_gt
 #define near____sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef near____sl_key_pure_cmp_ge
 #define near____sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef near____sl_key_cmp_eq
 #define near____sl_key_cmp_eq(k0, k1)       (near____sl_key_pure_cmp_eq(near____sl_key_purify(k0), near____sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near____sl_key_cmp_ne
 #define near____sl_key_cmp_ne(k0, k1)       (near____sl_key_pure_cmp_ne(near____sl_key_purify(k0), near____sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near____sl_key_cmp_lt
 #define near____sl_key_cmp_lt(k0, k1)       (near____sl_key_pure_cmp_lt(near____sl_key_purify(k0), near____sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near____sl_key_cmp_le
 #define near____sl_key_cmp_le(k0, k1)       (near____sl_key_pure_cmp_le(near____sl_key_purify(k0), near____sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near____sl_key_cmp_gt
 #define near____sl_key_cmp_gt(k0, k1)       (near____sl_key_pure_cmp_gt(near____sl_key_purify(k0), near____sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef near____sl_key_cmp_ge
 #define near____sl_key_cmp_ge(k0, k1)       (near____sl_key_pure_cmp_ge(near____sl_key_purify(k0), near____sl_key_purify(k1)))  /* sl_macro */
#endif


/* default random key */
#ifdef near____sl_key_integer
# if !defined(near____sl_key_val_srand) || !defined(near____sl_key_val_rand) || !defined(near____sl_key_val_rand_minmax)
#  undef near____sl_key_val_srand
#  undef near____sl_key_val_rand
#  undef near____sl_key_val_rand_minmax
#  define near____sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define near____sl_key_val_rand()                     ((near____sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define near____sl_key_val_rand_minmax(_min_, _max_)  ((near____sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef near____SL_DATA0_IGNORE
# undef near____SL_DATA0
#endif
#ifdef near____SL_DATA1_IGNORE
# undef near____SL_DATA1
#endif
#ifdef near____SL_DATA2_IGNORE
# undef near____SL_DATA2
#endif
#ifdef near____SL_DATA3_IGNORE
# undef near____SL_DATA3
#endif
#ifdef near____SL_DATA4_IGNORE
# undef near____SL_DATA4
#endif
#ifdef near____SL_DATA5_IGNORE
# undef near____SL_DATA5
#endif
#ifdef near____SL_DATA6_IGNORE
# undef near____SL_DATA6
#endif
#ifdef near____SL_DATA7_IGNORE
# undef near____SL_DATA7
#endif
#ifdef near____SL_DATA8_IGNORE
# undef near____SL_DATA8
#endif
#ifdef near____SL_DATA9_IGNORE
# undef near____SL_DATA9
#endif
#ifdef near____SL_DATA10_IGNORE
# undef near____SL_DATA10
#endif
#ifdef near____SL_DATA11_IGNORE
# undef near____SL_DATA11
#endif
#ifdef near____SL_DATA12_IGNORE
# undef near____SL_DATA12
#endif
#ifdef near____SL_DATA13_IGNORE
# undef near____SL_DATA13
#endif
#ifdef near____SL_DATA14_IGNORE
# undef near____SL_DATA14
#endif
#ifdef near____SL_DATA15_IGNORE
# undef near____SL_DATA15
#endif
#ifdef near____SL_DATA16_IGNORE
# undef near____SL_DATA16
#endif
#ifdef near____SL_DATA17_IGNORE
# undef near____SL_DATA17
#endif
#ifdef near____SL_DATA18_IGNORE
# undef near____SL_DATA18
#endif
#ifdef near____SL_DATA19_IGNORE
# undef near____SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro near____sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef near____sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef near____sl_data0_weight
# undef near____sl_data1_weight
# undef near____sl_data2_weight
# undef near____sl_data3_weight
# undef near____sl_data4_weight
# undef near____sl_data5_weight
# undef near____sl_data6_weight
# undef near____sl_data7_weight
# undef near____sl_data8_weight
# undef near____sl_data9_weight
# undef near____sl_data10_weight
# undef near____sl_data11_weight
# undef near____sl_data12_weight
# undef near____sl_data13_weight
# undef near____sl_data14_weight
# undef near____sl_data15_weight
# undef near____sl_data16_weight
# undef near____sl_data17_weight
# undef near____sl_data18_weight
# undef near____sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable near____sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(near____sl_data0_weight) && !defined(near____SL_DATA0)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data1_weight) && !defined(near____SL_DATA1)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data2_weight) && !defined(near____SL_DATA2)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data3_weight) && !defined(near____SL_DATA3)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data4_weight) && !defined(near____SL_DATA4)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data5_weight) && !defined(near____SL_DATA5)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data6_weight) && !defined(near____SL_DATA6)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data7_weight) && !defined(near____SL_DATA7)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data8_weight) && !defined(near____SL_DATA8)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data9_weight) && !defined(near____SL_DATA9)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data10_weight) && !defined(near____SL_DATA10)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data11_weight) && !defined(near____SL_DATA11)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data12_weight) && !defined(near____SL_DATA12)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data13_weight) && !defined(near____SL_DATA13)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data14_weight) && !defined(near____SL_DATA14)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data15_weight) && !defined(near____SL_DATA15)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data16_weight) && !defined(near____SL_DATA16)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data17_weight) && !defined(near____SL_DATA17)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data18_weight) && !defined(near____SL_DATA18)
# undef near____sl_elem_weight
#endif
#if defined(near____sl_data19_weight) && !defined(near____SL_DATA19)
# undef near____sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro near____FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef near____SL_DATA0
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data0_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA1
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data1_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA2
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data2_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA3
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data3_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA4
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data4_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA5
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data5_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA6
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data6_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA7
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data7_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA8
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data8_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA9
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data9_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA10
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data10_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA11
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data11_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA12
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data12_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA13
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data13_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA14
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data14_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA15
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data15_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA16
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data16_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA17
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data17_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA18
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data18_flex
#   define near____FLECKS_GUARD
#  endif
# endif
#endif
#ifdef near____SL_DATA19
# ifdef near____FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef near____sl_data19_flex
#   define near____FLECKS_GUARD
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






#define near____SPEC_TLOC

typedef near____sl_int_type_c near____spec_int_t;

typedef int near____spec_proc_t;

#define near____SPEC_LOC_NONE   -1
#define near____SPEC_PROC_NONE  MPI_PROC_NULL

typedef void *near____spec_tloc_data_t;
typedef void *near____spec_tproc_data_t;

struct near_____elements_t;

typedef struct near_____elements_t *near____spec_elem_buf_t;

typedef struct near_____elements_t near____spec_elem_t;

typedef near____sl_int_type_c near____spec_elem_index_t;

#define near____spec_elem_set_n(_e_, _n_)     near____elem_set_size((_e_), (_n_))
#define near____spec_elem_get_n(_e_)          near____elem_get_size((_e_))
#define near____spec_elem_set_nmax(_e_, _n_)  near____elem_set_max_size((_e_), (_n_))
#define near____spec_elem_get_nmax(_e_)       near____elem_get_max_size((_e_))

#define near____spec_elem_set_buf(_e_, _b_)   *(_e_) = *(_b_)
#define near____spec_elem_get_buf(_e_)        (_e_)

#define near____spec_elem_copy_at(_se_, _sat_, _de_, _dat_) \
  elem_copy_at((_se_), (_sat_), (_de_), (_dat_))

#define near____spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) \
  elem_xchange_at((_s0_), (_s0at_), (_s1_), (_s1at_), (_t_))






/* tproc count */

/* sp_macro near____SPEC_DECLARE_TPROC_COUNT_DB */
#define near____SPEC_DECLARE_TPROC_COUNT_DB \
  struct { near____spec_elem_index_t i; near____spec_proc_t p; } spec0cd;

/* sp_macro near____SPEC_DO_TPROC_COUNT_DB */
#define near____SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec0cd.i = 0; spec0cd.i < near____spec_elem_get_n(_b_); ++spec0cd.i) { \
    spec0cd.p = (_tp_)(near____spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == near____SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROC_COUNT_DB */
#define near____SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts) \
{ \
  near____SPEC_DECLARE_TPROC_COUNT_DB \
  near____SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro near____SPEC_DECLARE_TPROC_COUNT_IP */
#define near____SPEC_DECLARE_TPROC_COUNT_IP \
  struct { near____spec_elem_index_t i, t; near____spec_proc_t p; } spec0ci;

/* sp_macro near____SPEC_DO_TPROC_COUNT_IP */
#define near____SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < near____spec_elem_get_n(_b_); ++spec0ci.i) { \
    spec0ci.p = (_tp_)(near____spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == near____SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) near____spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  near____spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro near____SPEC_FUNC_TPROC_COUNT_IP */
#define near____SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts) \
{ \
  near____SPEC_DECLARE_TPROC_COUNT_IP \
  near____SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro near____SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define near____SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { near____spec_elem_index_t i; near____spec_proc_t p; } spec1cd;

/* sp_macro near____SPEC_DO_TPROC_MOD_COUNT_DB */
#define near____SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec1cd.i = 0; spec1cd.i < near____spec_elem_get_n(_b_); ++spec1cd.i) { \
    spec1cd.p = (_tp_)(near____spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == near____SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define near____SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts) \
{ \
  near____SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  near____SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro near____SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define near____SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { near____spec_elem_index_t i, t; near____spec_proc_t p; } spec1ci;

/* sp_macro near____SPEC_DO_TPROC_MOD_COUNT_IP */
#define near____SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < near____spec_elem_get_n(_b_); ++spec1ci.i) { \
    spec1ci.p = (_tp_)(near____spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == near____SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) near____spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  near____spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro near____SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define near____SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts) \
{ \
  near____SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  near____SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro near____SPEC_DECLARE_TPROCS_COUNT_DB */
#define near____SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { near____spec_elem_index_t i; near____spec_int_t j, n; } spec2cd;

/* sp_macro near____SPEC_DO_TPROCS_COUNT_DB */
#define near____SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec2cd.i = 0; spec2cd.i < near____spec_elem_get_n(_b_); ++spec2cd.i) { \
    spec2cd.n = (_tp_)(near____spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.n; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROCS_COUNT_DB */
#define near____SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts, near____spec_proc_t *procs) \
{ \
  near____SPEC_DECLARE_TPROCS_COUNT_DB \
  near____SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro near____SPEC_DECLARE_TPROCS_COUNT_IP */
#define near____SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { near____spec_elem_index_t i, t; near____spec_int_t j, n; } spec2ci;

/* sp_macro near____SPEC_DO_TPROCS_COUNT_IP */
#define near____SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < near____spec_elem_get_n(_b_); ++spec2ci.i) { \
    spec2ci.n = (_tp_)(near____spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), (_ps_)); \
    if (spec2ci.n <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.n; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) near____spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  near____spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro near____SPEC_FUNC_TPROCS_COUNT_IP */
#define near____SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts, near____spec_proc_t *procs) \
{ \
  near____SPEC_DECLARE_TPROCS_COUNT_IP \
  near____SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro near____SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define near____SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { near____spec_elem_index_t i; near____spec_int_t j, n; } spec3cd;

/* sp_macro near____SPEC_DO_TPROCS_MOD_COUNT_DB */
#define near____SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec3cd.i = 0; spec3cd.i < near____spec_elem_get_n(_b_); ++spec3cd.i) \
  { \
    spec3cd.n = (_tp_)(near____spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.n; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define near____SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts, near____spec_proc_t *procs) \
{ \
  near____SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  near____SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro near____SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define near____SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { near____spec_elem_index_t i, t; near____spec_int_t j, n; } spec3ci;

/* sp_macro near____SPEC_DO_TPROCS_MOD_COUNT_IP */
#define near____SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < near____spec_elem_get_n(_b_); ++spec3ci.i) { \
    spec3ci.n = (_tp_)(near____spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), (_ps_), NULL); \
    if (spec3ci.n <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.n; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) near____spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  near____spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro near____SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define near____SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts, near____spec_proc_t *procs) \
{ \
  near____SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  near____SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tproc rearrange */

/* sp_macro near____SPEC_DECLARE_TPROC_REARRANGE_DB */
#define near____SPEC_DECLARE_TPROC_REARRANGE_DB \
  struct { near____spec_elem_index_t i; near____spec_proc_t p; } spec0d;

/* sp_macro near____SPEC_DO_TPROC_REARRANGE_DB */
#define near____SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  for (spec0d.i = 0; spec0d.i < near____spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tp_)(near____spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == near____SPEC_PROC_NONE) continue; \
    near____spec_elem_copy_at((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p]); \
    ++(_ds_)[spec0d.p]; \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROC_REARRANGE_DB */
#define near____SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tproc_data_t tproc_data, int *displs) \
{ \
  near____SPEC_DECLARE_TPROC_REARRANGE_DB \
  near____SPEC_DO_TPROC_REARRANGE_DB(_tp_, tproc_data, s, d, displs); \
}

/* sp_macro near____SPEC_DECLARE_TPROC_REARRANGE_IP */
#define near____SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { near____spec_elem_index_t e, i, j; near____spec_proc_t p, np; } spec0i;

/* sp_macro near____SPEC_DO_TPROC_REARRANGE_IP */
#define near____SPEC_DO_TPROC_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_)  do { \
  for (spec0i.e = 0, spec0i.i = 0; spec0i.i < (_n_); ++spec0i.i) { \
    spec0i.e += (_cs_)[spec0i.i]; \
    spec0i.j = (_ds_)[spec0i.i]; \
    while (spec0i.j < spec0i.e) { \
      spec0i.p = (_tp_)(near____spec_elem_get_buf(_b_), spec0i.j, _tpd_); \
      while (spec0i.p != spec0i.i) { \
        spec0i.np = (_tp_)(near____spec_elem_get_buf(_b_), (_ds_)[spec0i.p], _tpd_); \
        if (spec0i.np != spec0i.p) near____spec_elem_exchange_at((_b_), (_ds_)[spec0i.p], (_b_), spec0i.j, (_xb_)); \
        ++(_ds_)[spec0i.p]; \
        spec0i.p = spec0i.np; \
      } \
      ++spec0i.j; \
    } \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROC_REARRANGE_IP */
#define near____SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_rearrange_ip(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tproc_data_t tproc_data, int *displs, int *counts, near____spec_int_t n) \
{ \
  near____SPEC_DECLARE_TPROC_REARRANGE_IP \
  near____SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro near____SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define near____SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  struct { near____spec_elem_index_t i; near____spec_proc_t p; } spec1d;

/* sp_macro near____SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define near____SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < near____spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(near____spec_elem_get_buf(_sb_), spec1d.i, _tpd_, near____spec_elem_get_buf(_ib_)); \
      if (spec1d.p == near____SPEC_PROC_NONE) continue; \
      near____spec_elem_copy_at((_ib_), 0, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < near____spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(near____spec_elem_get_buf(_sb_), spec1d.i, _tpd_, NULL); \
      if (spec1d.p == near____SPEC_PROC_NONE) continue; \
      near____spec_elem_copy_at((_sb_), spec1d.i, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define near____SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tproc_data_t tproc_data, int *displs, near____spec_elem_t *mod) \
{ \
  near____SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  near____SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro near____SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define near____SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { near____spec_elem_index_t e, i, j; near____spec_proc_t p, np; } spec1i;

/* sp_macro near____SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define near____SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
      spec1i.e += (_cs_)[spec1i.i]; \
      spec1i.j = (_ds_)[spec1i.i]; \
      while (spec1i.j < spec1i.e) { \
        spec1i.p = (_tp_)(near____spec_elem_get_buf(_b_), spec1i.j, _tpd_, near____spec_elem_get_buf(_ib_)); \
        near____spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(near____spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, near____spec_elem_get_buf(_ib_)); \
          if (spec1i.np != spec1i.p) { \
            near____spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
            near____spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
          } else near____spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
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
        spec1i.p = (_tp_)(near____spec_elem_get_buf(_b_), spec1i.j, _tpd_, NULL); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(near____spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, NULL); \
          if (spec1i.np != spec1i.p) near____spec_elem_exchange_at((_b_), (_ds_)[spec1i.p], (_b_), spec1i.j, (_xb_)); \
          ++(_ds_)[spec1i.p]; \
          spec1i.p = spec1i.np; \
        } \
        ++spec1i.j; \
      } \
    } \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define near____SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tproc_data_t tproc_data, int *displs, int *counts, near____spec_int_t n, near____spec_elem_t *mod) \
{ \
  near____SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  near____SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro near____SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define near____SPEC_DECLARE_TPROCS_REARRANGE_DB \
  struct { near____spec_elem_index_t i; near____spec_int_t j, n; } spec2d;

/* sp_macro near____SPEC_DO_TPROCS_REARRANGE_DB */
#define near____SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  for (spec2d.i = 0; spec2d.i < near____spec_elem_get_n(_sb_); ++spec2d.i) { \
    spec2d.n = (_tp_)(near____spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.n; ++spec2d.j) { \
      near____spec_elem_copy_at((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]]); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROCS_REARRANGE_DB */
#define near____SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tproc_data_t tproc_data, int *displs, near____spec_proc_t *procs) \
{ \
  near____SPEC_DECLARE_TPROCS_REARRANGE_DB \
  near____SPEC_DO_TPROCS_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro near____SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define near____SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { near____spec_elem_index_t e, j, fe, fc, le, lc; near____spec_int_t i, n, f, l, o; } spec2i;

/* sp_macro near____SPEC_DO_TPROCS_REARRANGE_IP */
#define near____SPEC_DO_TPROCS_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_)  do { \
  spec2i.f = 0; spec2i.fe = (_cs_)[0]; spec2i.fc = near____spec_elem_get_n(_b_); \
  while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; } \
  spec2i.l = 0; spec2i.le = (_cs_)[0]; spec2i.lc = near____spec_elem_get_n(_b_) - 1; \
  while (spec2i.lc >= spec2i.le) { ++spec2i.l; spec2i.le += (_cs_)[spec2i.l]; } \
  for (spec2i.e = 0, spec2i.i = 0; spec2i.i < (_n_); ++spec2i.i) { \
    spec2i.e += (_cs_)[spec2i.i]; \
    spec2i.j = (_ds_)[spec2i.i]; \
    while (spec2i.j < spec2i.e) { \
      spec2i.n = (_tp_)(near____spec_elem_get_buf(_b_), spec2i.j, (_tpd_), (_ps_)); \
      spec2i.o = -1; \
      while (spec2i.n > 0) { \
        --spec2i.n; \
        if ((_ps_)[spec2i.n] == spec2i.i && spec2i.o < 0) spec2i.o = spec2i.n; \
        else if ((_ds_)[(_ps_)[spec2i.n]] < spec2i.fc) { \
          spec2i.l = spec2i.f; spec2i.le = spec2i.fe; spec2i.lc = spec2i.fc; \
          if (spec2i.fc < spec2i.fe) { \
            near____spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_b_), spec2i.fc); \
            ++spec2i.fc; \
          } else near____spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec2i.n]] == spec2i.fc) ++spec2i.fc; \
        if (spec2i.j != (_ds_)[(_ps_)[spec2i.n]]) near____spec_elem_copy_at((_b_), spec2i.j, (_b_), (_ds_)[(_ps_)[spec2i.n]]); \
        ++(_ds_)[(_ps_)[spec2i.n]]; \
        while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; spec2i.fc = (_ds_)[spec2i.f]; } \
      } \
      if (spec2i.o < 0) { \
        if (spec2i.lc < spec2i.le) {  \
          near____spec_elem_copy_at((_b_), spec2i.lc, (_b_), spec2i.j); \
          spec2i.f = spec2i.l; spec2i.fe = spec2i.le; spec2i.fc = spec2i.lc; \
          --spec2i.lc; \
          while (spec2i.l > 0 && spec2i.lc < (_ds_)[spec2i.l]) { spec2i.le -= (_cs_)[spec2i.l]; spec2i.lc = spec2i.le - 1; --spec2i.l; } \
        } else near____spec_elem_copy_at((_xb_), 0, (_b_), spec2i.j); \
      } \
      spec2i.j = (_ds_)[spec2i.i]; \
    } \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROCS_REARRANGE_IP */
#define near____SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_ip(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tproc_data_t tproc_data, int *displs, int *counts, near____spec_int_t n, near____spec_proc_t *procs) \
{ \
  near____SPEC_DECLARE_TPROCS_REARRANGE_IP \
  near____SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro near____SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define near____SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  struct { near____spec_elem_index_t i; near____spec_int_t j, n; } spec3d;

/* sp_macro near____SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define near____SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  if (_ib_) { \
    for (spec3d.i = 0; spec3d.i < near____spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(near____spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), near____spec_elem_get_buf(_ib_)); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        near____spec_elem_copy_at((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } else { \
    for (spec3d.i = 0; spec3d.i < near____spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(near____spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), NULL); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        near____spec_elem_copy_at((_sb_), spec3d.i, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define near____SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tproc_data_t tproc_data, int *displs, near____spec_proc_t *procs, near____spec_elem_t *mod) \
{ \
  near____SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  near____SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro near____SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define near____SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { near____spec_elem_index_t e, j, fe, fc, le, lc; near____spec_int_t i, n, f, l, o; } spec3i;

/* sp_macro near____SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define near____SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  if (_ib_) { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = near____spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = near____spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(near____spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), near____spec_elem_get_buf(_ib_)); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              near____spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else near____spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          near____spec_elem_copy_at((_ib_), spec3i.n, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            near____spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else near____spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } else { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = near____spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = near____spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(near____spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), NULL); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              near____spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else near____spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          if (spec3i.j != (_ds_)[(_ps_)[spec3i.n]]) near____spec_elem_copy_at((_b_), spec3i.j, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            near____spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else near____spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define near____SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tproc_data_t tproc_data, int *displs, int *counts, near____spec_int_t n, near____spec_proc_t *procs, near____spec_elem_t *mod) \
{ \
  near____SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  near____SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}

/* sp_macro near____SPEC_DEFINE_TPROC */
#define near____SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  near____SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro near____SPEC_DEFINE_TPROC_MOD */
#define near____SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  near____SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro near____SPEC_DEFINE_TPROCS */
#define near____SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  near____SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro near____SPEC_DEFINE_TPROCS_MOD */
#define near____SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  near____SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  near____SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro near____SPEC_EXT_PARAM_TPROC near____SPEC_EXT_PARAM_TPROC_NULL near____SPEC_EXT_PARAM_TPROC_MOD near____SPEC_EXT_PARAM_TPROC_MOD_NULL near____SPEC_EXT_PARAM_TPROCS near____SPEC_EXT_PARAM_TPROCS_NULL near____SPEC_EXT_PARAM_TPROCS_MOD near____SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define near____SPEC_EXT_PARAM_TPROC(_name_)       _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip
#define near____SPEC_EXT_PARAM_TPROC_NULL          NULL, NULL, NULL, NULL
#define near____SPEC_EXT_PARAM_TPROC_MOD(_name_)   _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip
#define near____SPEC_EXT_PARAM_TPROC_MOD_NULL      NULL, NULL, NULL, NULL
#define near____SPEC_EXT_PARAM_TPROCS(_name_)      _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip
#define near____SPEC_EXT_PARAM_TPROCS_NULL         NULL, NULL, NULL, NULL
#define near____SPEC_EXT_PARAM_TPROCS_MOD(_name_)  _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip
#define near____SPEC_EXT_PARAM_TPROCS_MOD_NULL     NULL, NULL, NULL, NULL


/* sp_type near____spec_tproc_f near____spec_tproc_count_f near____spec_tproc_rearrange_db_f near____spec_tproc_rearrange_ip_f */
typedef near____spec_proc_t near____spec_tproc_f(near____spec_elem_buf_t b, near____spec_elem_index_t x, near____spec_tproc_data_t tproc_data);
typedef void near____spec_tproc_count_f(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts);
typedef void near____spec_tproc_rearrange_db_f(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tproc_data_t tproc_data, int *displs);
typedef void near____spec_tproc_rearrange_ip_f(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tproc_data_t tproc_data, int *displs, int *counts, near____spec_int_t n);

/* sp_type near____spec_tproc_mod_f near____spec_tproc_mod_count_f near____spec_tproc_mod_rearrange_db_f near____spec_tproc_mod_rearrange_ip_f */
typedef near____spec_proc_t near____spec_tproc_mod_f(near____spec_elem_buf_t b, near____spec_elem_index_t x, near____spec_tproc_data_t tproc_data, near____spec_elem_buf_t mod);
typedef void near____spec_tproc_mod_count_f(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts);
typedef void near____spec_tproc_mod_rearrange_db_f(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tproc_data_t tproc_data, int *displs, near____spec_elem_t *mod);
typedef void near____spec_tproc_mod_rearrange_ip_f(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tproc_data_t tproc_data, int *displs, int *counts, near____spec_int_t n, near____spec_elem_t *mod);

/* sp_type near____spec_tprocs_f near____spec_tprocs_count_f near____spec_tprocs_rearrange_db_f near____spec_tprocs_rearrange_ip_f */
typedef near____spec_int_t near____spec_tprocs_f(near____spec_elem_buf_t b, near____spec_elem_index_t x, near____spec_tproc_data_t tproc_data, near____spec_proc_t *procs);
typedef void near____spec_tprocs_count_f(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts, near____spec_proc_t *procs);
typedef void near____spec_tprocs_rearrange_db_f(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tproc_data_t tproc_data, int *displs, near____spec_proc_t *procs);
typedef void near____spec_tprocs_rearrange_ip_f(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tproc_data_t tproc_data, int *displs, int *counts, near____spec_int_t n, near____spec_proc_t *procs);

/* sp_type near____spec_tprocs_mod_f near____spec_tprocs_mod_count_f near____spec_tprocs_mod_rearrange_db_f near____spec_tprocs_mod_rearrange_ip_f */
typedef near____spec_int_t near____spec_tprocs_mod_f(near____spec_elem_buf_t b, near____spec_elem_index_t x, near____spec_tproc_data_t tproc_data, near____spec_proc_t *procs, near____spec_elem_buf_t mod);
typedef void near____spec_tprocs_mod_count_f(near____spec_elem_t *s, near____spec_tproc_data_t tproc_data, int *counts, near____spec_proc_t *procs);
typedef void near____spec_tprocs_mod_rearrange_db_f(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tproc_data_t tproc_data, int *displs, near____spec_proc_t *procs, near____spec_elem_t *mod);
typedef void near____spec_tprocs_mod_rearrange_ip_f(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tproc_data_t tproc_data, int *displs, int *counts, near____spec_int_t n, near____spec_proc_t *procs, near____spec_elem_t *mod);

/* sp_type near____spec_tproc_reset_f */
typedef void near____spec_tproc_reset_f(near____spec_tproc_data_t tproc_data);


/* enable tloc features */
#ifdef near____SPEC_TLOC

/* sp_macro near____SPEC_TLOC near____SPEC_LOC_NONE */


/* tloc rearrange */

/* sp_macro near____SPEC_DECLARE_TLOC_REARRANGE_DB */
#define near____SPEC_DECLARE_TLOC_REARRANGE_DB \
  struct { near____spec_int_t i, p; } spec0d;

/* sp_macro near____SPEC_DO_TLOC_REARRANGE_DB */
#define near____SPEC_DO_TLOC_REARRANGE_DB(_tl_, _tld_, _sb_, _db_)  do { \
  for (spec0d.i = 0; spec0d.i < near____spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tl_)(near____spec_elem_get_buf(_sb_), spec0d.i, _tld_); \
    if (spec0d.p == near____SPEC_LOC_NONE) continue; \
    near____spec_elem_copy_at((_sb_), spec0d.i, (_db_), spec0d.p); \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TLOC_REARRANGE_DB */
#define near____SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_rearrange_db(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tloc_data_t tloc_data) \
{ \
  near____SPEC_DECLARE_TLOC_REARRANGE_DB \
  near____SPEC_DO_TLOC_REARRANGE_DB(_tl_, tloc_data, s, d); \
}

/* sp_macro near____SPEC_DECLARE_TLOC_REARRANGE_IP */
#define near____SPEC_DECLARE_TLOC_REARRANGE_IP \
  struct { near____spec_int_t i, p, np; } spec0i;

/* sp_macro near____SPEC_DO_TLOC_REARRANGE_IP */
#define near____SPEC_DO_TLOC_REARRANGE_IP(_tl_, _tld_, _b_, _xb_)  do { \
  for (spec0i.i = 0; spec0i.i < near____spec_elem_get_n(_b_); ++spec0i.i) { \
    spec0i.p = (_tl_)(near____spec_elem_get_buf(_b_), spec0i.i, _tld_); \
    if (spec0i.p == near____SPEC_LOC_NONE) continue; \
    while (spec0i.i != spec0i.p) { \
      spec0i.np = (_tl_)(near____spec_elem_get_buf(_b_), spec0i.p, _tld_); \
      if (spec0i.np == near____SPEC_LOC_NONE) { near____spec_elem_copy_at((_b_), spec0i.i, (_b_), spec0i.p); break; } \
      near____spec_elem_exchange_at((_b_), spec0i.i, (_b_), spec0i.p, (_xb_)); \
      spec0i.p = spec0i.np; \
    } \
  } } while (0)

/* sp_macro near____SPEC_FUNC_TLOC_REARRANGE_IP */
#define near____SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_rearrange_ip(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tloc_data_t tloc_data) \
{ \
  near____SPEC_DECLARE_TLOC_REARRANGE_IP \
  near____SPEC_DO_TLOC_REARRANGE_IP(_tl_, tloc_data, s, x); \
}


/* tloc_mod_mod rearrange */

/* sp_macro near____SPEC_DECLARE_TLOC_MOD_REARRANGE_DB */
#define near____SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  struct { near____spec_int_t i, p; } spec1d;

/* sp_macro near____SPEC_DO_TLOC_MOD_REARRANGE_DB */
#define near____SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, _tld_, _sb_, _db_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < near____spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(near____spec_elem_get_buf(_sb_), spec1d.i, _tld_, near____spec_elem_get_buf(_ib_)); \
      if (spec1d.p == near____SPEC_LOC_NONE) continue; \
      near____spec_elem_copy_at((_ib_), 0, (_db_), spec1d.p); \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < near____spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(near____spec_elem_get_buf(_sb_), spec1d.i, _tld_, NULL); \
      if (spec1d.p == near____SPEC_LOC_NONE) continue; \
      near____spec_elem_copy_at((_sb_), spec1d.i, (_db_), spec1d.p); \
    } \
  } } while (0) 

/* sp_macro near____SPEC_FUNC_TLOC_MOD_REARRANGE_DB */
#define near____SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_mod_rearrange_db(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tloc_data_t tloc_data, near____spec_elem_t *mod) \
{ \
  near____SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  near____SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, tloc_data, s, d, mod); \
}

/* sp_macro near____SPEC_DECLARE_TLOC_MOD_REARRANGE_IP */
#define near____SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  struct { near____spec_int_t i, p, np; } spec1i;

/* sp_macro near____SPEC_DO_TLOC_MOD_REARRANGE_IP */
#define near____SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, _tld_, _b_, _xb_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.i = 0; spec1i.i < near____spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(near____spec_elem_get_buf(_b_), spec1i.i, _tld_, near____spec_elem_get_buf(_ib_)); \
      if (spec1i.p == near____SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(near____spec_elem_get_buf(_b_), spec1i.p, _tld_, near____spec_elem_get_buf(_xb_)); \
        if (spec1i.np == near____SPEC_LOC_NONE) break; \
        near____spec_elem_copy_at((_ib_), 0, (_b_), spec1i.p); \
        near____spec_elem_copy_at((_xb_), 0, (_ib_), 0); \
        spec1i.p = spec1i.np; \
      } \
      near____spec_elem_copy_at((_ib_), 0, (_b_), spec1i.i); \
    } \
  } else { \
    for (spec1i.i = 0; spec1i.i < near____spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(near____spec_elem_get_buf(_b_), spec1i.i, _tld_, NULL); \
      if (spec1i.p == near____SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(near____spec_elem_get_buf(_b_), spec1i.p, _tld_, NULL); \
        if (spec1i.np == near____SPEC_LOC_NONE) { near____spec_elem_copy_at((_b_), spec1i.i, (_b_), spec1i.p); break; } \
        near____spec_elem_exchange_at((_b_), spec1i.i, (_b_), spec1i.p, (_xb_)); \
        spec1i.p = spec1i.np; \
      } \
    } \
 } } while (0) 

/* sp_macro near____SPEC_FUNC_TLOC_MOD_REARRANGE_IP */
#define near____SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_mod_rearrange_ip(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tloc_data_t tloc_data, near____spec_elem_t *mod) \
{ \
  near____SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  near____SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, tloc_data, s, x, mod); \
}

/* sp_macro near____SPEC_DEFINE_TLOC */
#define near____SPEC_DEFINE_TLOC(_name_, _tl_, _s_...) \
  near____SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_) \
  near____SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro near____SPEC_DEFINE_TLOC_MOD */
#define near____SPEC_DEFINE_TLOC_MOD(_name_, _tl_, _s_...) \
  near____SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_) \
  near____SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro near____SPEC_EXT_PARAM_TLOC near____SPEC_EXT_PARAM_TLOC_NULL near____SPEC_EXT_PARAM_TLOC_MOD near____SPEC_EXT_PARAM_TLOC_MOD_NULL */
#define near____SPEC_EXT_PARAM_TLOC(_name_)      _name_##_tloc_rearrange_db, _name_##_tloc_rearrange_ip
#define near____SPEC_EXT_PARAM_TLOC_NULL         NULL, NULL
#define near____SPEC_EXT_PARAM_TLOC_MOD(_name_)  _name_##_tloc_mod_rearrange_db, _name_##_tloc_mod_rearrange_ip
#define near____SPEC_EXT_PARAM_TLOC_MOD_NULL     NULL, NULL


/* sp_type near____spec_tloc_f near____spec_tloc_rearrange_db_f near____spec_tloc_rearrange_ip_f */
typedef near____spec_elem_index_t near____spec_tloc_f(near____spec_elem_buf_t b, near____spec_elem_index_t x, near____spec_tloc_data_t tloc_data);
typedef void near____spec_tloc_rearrange_db_f(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tloc_data_t tloc_data);
typedef void near____spec_tloc_rearrange_ip_f(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tloc_data_t tloc_data);

/* sp_type near____spec_tloc_mod_f near____spec_tloc_mod_rearrange_db_f near____spec_tloc_mod_rearrange_ip_f */
typedef near____spec_elem_index_t near____spec_tloc_mod_f(near____spec_elem_buf_t b, near____spec_elem_index_t x, near____spec_tloc_data_t tloc_data, near____spec_elem_buf_t mod);
typedef void near____spec_tloc_mod_rearrange_db_f(near____spec_elem_t *s, near____spec_elem_t *d, near____spec_tloc_data_t tloc_data, near____spec_elem_t *mod);
typedef void near____spec_tloc_mod_rearrange_ip_f(near____spec_elem_t *s, near____spec_elem_t *x, near____spec_tloc_data_t tloc_data, near____spec_elem_t *mod);


#endif /* near____SPEC_TLOC */






#ifdef SL_USE_MPI
# include <mpi.h>
#endif


/* sl_type near____slint_t near____slint */
typedef near____sl_int_type_c near____slint_t, near____slint;  /* deprecated 'near____slint' */

#define near____slint_fmt   near____sl_int_type_fmt    /* sl_macro */

/* sl_type near____slindex_t */
typedef near____sl_index_type_c near____slindex_t;

#define near____sindex_fmt  near____sl_index_type_fmt  /* sl_macro */

/* sl_type near____slkey_t */
typedef near____sl_key_type_c near____slkey_t;

/* sl_type near____slkey_pure_t near____slpkey_t */
typedef near____sl_key_pure_type_c near____slkey_pure_t, near____slpkey_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type near____sldata0_t */
#ifdef near____sl_data0_type_c
typedef near____sl_data0_type_c near____sldata0_t;
#endif
/* sl_type near____sldata1_t */
#ifdef near____sl_data1_type_c
typedef near____sl_data1_type_c near____sldata1_t;
#endif
/* sl_type near____sldata2_t */
#ifdef near____sl_data2_type_c
typedef near____sl_data2_type_c near____sldata2_t;
#endif
/* sl_type near____sldata3_t */
#ifdef near____sl_data3_type_c
typedef near____sl_data3_type_c near____sldata3_t;
#endif
/* sl_type near____sldata4_t */
#ifdef near____sl_data4_type_c
typedef near____sl_data4_type_c near____sldata4_t;
#endif
/* sl_type near____sldata5_t */
#ifdef near____sl_data5_type_c
typedef near____sl_data5_type_c near____sldata5_t;
#endif
/* sl_type near____sldata6_t */
#ifdef near____sl_data6_type_c
typedef near____sl_data6_type_c near____sldata6_t;
#endif
/* sl_type near____sldata7_t */
#ifdef near____sl_data7_type_c
typedef near____sl_data7_type_c near____sldata7_t;
#endif
/* sl_type near____sldata8_t */
#ifdef near____sl_data8_type_c
typedef near____sl_data8_type_c near____sldata8_t;
#endif
/* sl_type near____sldata9_t */
#ifdef near____sl_data9_type_c
typedef near____sl_data9_type_c near____sldata9_t;
#endif
/* sl_type near____sldata10_t */
#ifdef near____sl_data10_type_c
typedef near____sl_data10_type_c near____sldata10_t;
#endif
/* sl_type near____sldata11_t */
#ifdef near____sl_data11_type_c
typedef near____sl_data11_type_c near____sldata11_t;
#endif
/* sl_type near____sldata12_t */
#ifdef near____sl_data12_type_c
typedef near____sl_data12_type_c near____sldata12_t;
#endif
/* sl_type near____sldata13_t */
#ifdef near____sl_data13_type_c
typedef near____sl_data13_type_c near____sldata13_t;
#endif
/* sl_type near____sldata14_t */
#ifdef near____sl_data14_type_c
typedef near____sl_data14_type_c near____sldata14_t;
#endif
/* sl_type near____sldata15_t */
#ifdef near____sl_data15_type_c
typedef near____sl_data15_type_c near____sldata15_t;
#endif
/* sl_type near____sldata16_t */
#ifdef near____sl_data16_type_c
typedef near____sl_data16_type_c near____sldata16_t;
#endif
/* sl_type near____sldata17_t */
#ifdef near____sl_data17_type_c
typedef near____sl_data17_type_c near____sldata17_t;
#endif
/* sl_type near____sldata18_t */
#ifdef near____sl_data18_type_c
typedef near____sl_data18_type_c near____sldata18_t;
#endif
/* sl_type near____sldata19_t */
#ifdef near____sl_data19_type_c
typedef near____sl_data19_type_c near____sldata19_t;
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

/* sl_type near____slweight_t */
typedef near____sl_weight_type_c near____slweight_t;

#define near____slweight_fmt  near____sl_weight_type_fmt  /* sl_macro */

#if defined(near____sl_elem_weight) && defined(near____sl_weight_intequiv)
typedef near____sl_weight_type_c near____slcount_t;       /* sl_type near____slcount_t */
# define near____slcount_fmt  near____sl_weight_type_fmt  /* sl_macro */
#else
typedef near____sl_int_type_c near____slcount_t;
# define near____slcount_fmt  near____sl_int_type_fmt
#endif


/* sl_type near_____slpwkey_t near____slpwkey_t */
typedef struct near_____slpwkey_t
{
  near____slpkey_t pkey;
  near____slweight_t weight;

} near____slpwkey_t;


/* sl_type near_____elements_t near____elements_t */
typedef struct near_____elements_t
{
  near____slint_t size, max_size;
  near____slkey_t *keys;

#ifdef near____SL_INDEX
  near____slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef near____SL_DATA0
  near____sldata0_t *data0;
#endif
#ifdef near____SL_DATA1
  near____sldata1_t *data1;
#endif
#ifdef near____SL_DATA2
  near____sldata2_t *data2;
#endif
#ifdef near____SL_DATA3
  near____sldata3_t *data3;
#endif
#ifdef near____SL_DATA4
  near____sldata4_t *data4;
#endif
#ifdef near____SL_DATA5
  near____sldata5_t *data5;
#endif
#ifdef near____SL_DATA6
  near____sldata6_t *data6;
#endif
#ifdef near____SL_DATA7
  near____sldata7_t *data7;
#endif
#ifdef near____SL_DATA8
  near____sldata8_t *data8;
#endif
#ifdef near____SL_DATA9
  near____sldata9_t *data9;
#endif
#ifdef near____SL_DATA10
  near____sldata10_t *data10;
#endif
#ifdef near____SL_DATA11
  near____sldata11_t *data11;
#endif
#ifdef near____SL_DATA12
  near____sldata12_t *data12;
#endif
#ifdef near____SL_DATA13
  near____sldata13_t *data13;
#endif
#ifdef near____SL_DATA14
  near____sldata14_t *data14;
#endif
#ifdef near____SL_DATA15
  near____sldata15_t *data15;
#endif
#ifdef near____SL_DATA16
  near____sldata16_t *data16;
#endif
#ifdef near____SL_DATA17
  near____sldata17_t *data17;
#endif
#ifdef near____SL_DATA18
  near____sldata18_t *data18;
#endif
#ifdef near____SL_DATA19
  near____sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} near____elements_t;


/* sl_type near_____packed_element_t near____packed_element_t */
typedef struct near_____packed_element_t
{
  near____slkey_t key;

#ifdef near____SL_PACKED_INDEX
  near____slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef near____SL_DATA0
# ifdef near____sl_data0_flex
  near____sldata0_t data0[];
# else
  near____sldata0_t data0[near____sl_data0_size_c];
# endif
#endif
#ifdef near____SL_DATA1
# ifdef near____sl_data1_flex
  near____sldata1_t data1[];
# else
  near____sldata1_t data1[near____sl_data1_size_c];
# endif
#endif
#ifdef near____SL_DATA2
# ifdef near____sl_data2_flex
  near____sldata2_t data2[];
# else
  near____sldata2_t data2[near____sl_data2_size_c];
# endif
#endif
#ifdef near____SL_DATA3
# ifdef near____sl_data3_flex
  near____sldata3_t data3[];
# else
  near____sldata3_t data3[near____sl_data3_size_c];
# endif
#endif
#ifdef near____SL_DATA4
# ifdef near____sl_data4_flex
  near____sldata4_t data4[];
# else
  near____sldata4_t data4[near____sl_data4_size_c];
# endif
#endif
#ifdef near____SL_DATA5
# ifdef near____sl_data5_flex
  near____sldata5_t data5[];
# else
  near____sldata5_t data5[near____sl_data5_size_c];
# endif
#endif
#ifdef near____SL_DATA6
# ifdef near____sl_data6_flex
  near____sldata6_t data6[];
# else
  near____sldata6_t data6[near____sl_data6_size_c];
# endif
#endif
#ifdef near____SL_DATA7
# ifdef near____sl_data7_flex
  near____sldata7_t data7[];
# else
  near____sldata7_t data7[near____sl_data7_size_c];
# endif
#endif
#ifdef near____SL_DATA8
# ifdef near____sl_data8_flex
  near____sldata8_t data8[];
# else
  near____sldata8_t data8[near____sl_data8_size_c];
# endif
#endif
#ifdef near____SL_DATA9
# ifdef near____sl_data9_flex
  near____sldata9_t data9[];
# else
  near____sldata9_t data9[near____sl_data9_size_c];
# endif
#endif
#ifdef near____SL_DATA10
# ifdef near____sl_data10_flex
  near____sldata10_t data10[];
# else
  near____sldata10_t data10[near____sl_data10_size_c];
# endif
#endif
#ifdef near____SL_DATA11
# ifdef near____sl_data11_flex
  near____sldata11_t data11[];
# else
  near____sldata11_t data11[near____sl_data11_size_c];
# endif
#endif
#ifdef near____SL_DATA12
# ifdef near____sl_data12_flex
  near____sldata12_t data12[];
# else
  near____sldata12_t data12[near____sl_data12_size_c];
# endif
#endif
#ifdef near____SL_DATA13
# ifdef near____sl_data13_flex
  near____sldata13_t data13[];
# else
  near____sldata13_t data13[near____sl_data13_size_c];
# endif
#endif
#ifdef near____SL_DATA14
# ifdef near____sl_data14_flex
  near____sldata14_t data14[];
# else
  near____sldata14_t data14[near____sl_data14_size_c];
# endif
#endif
#ifdef near____SL_DATA15
# ifdef near____sl_data15_flex
  near____sldata15_t data15[];
# else
  near____sldata15_t data15[near____sl_data15_size_c];
# endif
#endif
#ifdef near____SL_DATA16
# ifdef near____sl_data16_flex
  near____sldata16_t data16[];
# else
  near____sldata16_t data16[near____sl_data16_size_c];
# endif
#endif
#ifdef near____SL_DATA17
# ifdef near____sl_data17_flex
  near____sldata17_t data17[];
# else
  near____sldata17_t data17[near____sl_data17_size_c];
# endif
#endif
#ifdef near____SL_DATA18
# ifdef near____sl_data18_flex
  near____sldata18_t data18[];
# else
  near____sldata18_t data18[near____sl_data18_size_c];
# endif
#endif
#ifdef near____SL_DATA19
# ifdef near____sl_data19_flex
  near____sldata19_t data19[];
# else
  near____sldata19_t data19[near____sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} near____packed_element_t;


/* sl_type near_____packed_elements_t near____packed_elements_t */
typedef struct near_____packed_elements_t
{
  near____slint_t size, max_size;
  
  near____packed_element_t *elements;
  
} near____packed_elements_t;


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


/* sl_type near_____classification_info_t near____classification_info_t near____classification_info */
typedef struct near_____classification_info_t
{
  near____slint_t nclasses;
  near____slkey_pure_t *keys;
  near____slint_t *counts;
  near____slint_t *masks;

  /* */
  near____slint_t *all_local_sizes;
  near____slint_t *local_lt_eq_counts;
  near____slint_t *all_local_lt_eq_counts;

} near____classification_info_t, near____classification_info;  /* deprecated 'near____classification_info' */


/* key2class, sl_type near____key2class_f */
typedef near____slint_t (*near____key2class_f)(near____slkey_t *, near____slint, void *);

/* pivot-element, sl_type near____pivot_f */
typedef near____slint_t (*near____pivot_f)(near____elements_t *);

/* sorting-network, sl_type near____sortnet_f near____sortnet_data_t */
typedef void *near____sortnet_data_t;
typedef near____slint_t (*near____sortnet_f)(near____slint_t size, near____slint_t rank, near____slint_t stage, near____sortnet_data_t snd, near____slint_t *up);

/* merge2, sl_type near____merge2x_f near____merge2X_f */
typedef near____slint_t (*near____merge2x_f)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
typedef near____slint_t (*near____merge2X_f)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx, near____elements_t *t);

/* sl_type near_____permute_generic_t near____permute_generic_t */
typedef struct near_____permute_generic_t
{
  int type;

  near____spec_tloc_f *tloc;
  near____spec_tloc_rearrange_db_f *tloc_rearrange_db;
  near____spec_tloc_rearrange_ip_f *tloc_rearrange_ip;

  near____spec_tloc_mod_f *tloc_mod;
  near____spec_tloc_mod_rearrange_db_f *tloc_mod_rearrange_db;
  near____spec_tloc_mod_rearrange_ip_f *tloc_mod_rearrange_ip;

} near____permute_generic_t;

/* sl_macro near____PERMUTE_GENERIC_DEFINE_TLOC near____PERMUTE_GENERIC_INIT_TLOC near____PERMUTE_GENERIC_INIT_EXT_TLOC */
#define near____PERMUTE_GENERIC_DEFINE_TLOC(_tl_, _s_...)      near____SPEC_DEFINE_TLOC(_tl_, _tl_, _s_)
#define near____PERMUTE_GENERIC_INIT_TLOC(_tl_)                { 1, _tl_, near____SPEC_EXT_PARAM_TLOC_NULL,  NULL, near____SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define near____PERMUTE_GENERIC_INIT_EXT_TLOC(_tl_)            { 1, _tl_, near____SPEC_EXT_PARAM_TLOC(_tl_), NULL, near____SPEC_EXT_PARAM_TLOC_MOD_NULL }

/* sl_macro near____PERMUTE_GENERIC_DEFINE_TLOC_MOD near____PERMUTE_GENERIC_INIT_TLOC_MOD near____PERMUTE_GENERIC_INIT_EXT_TLOC_MOD */
#define near____PERMUTE_GENERIC_DEFINE_TLOC_MOD(_tl_, _s_...)  near____SPEC_DEFINE_TLOC_MOD(_tl_, _tl_, _s_)
#define near____PERMUTE_GENERIC_INIT_TLOC_MOD(_tl_)            { 2, NULL, near____SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, near____SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define near____PERMUTE_GENERIC_INIT_EXT_TLOC_MOD(_tl_)        { 2, NULL, near____SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, near____SPEC_EXT_PARAM_TLOC_MOD(_tl_) }

/* sl_type near_____split_generic_t near____split_generic_t */
typedef struct near_____split_generic_t
{
  int type;

  near____spec_tproc_f *tproc;
  near____spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  near____spec_tproc_rearrange_db_f *tproc_rearrange_db;
  near____spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  near____spec_tproc_mod_f *tproc_mod;
  near____spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  near____spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  near____spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  near____spec_tprocs_f *tprocs;
  near____spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  near____spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  near____spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  near____spec_tprocs_mod_f *tprocs_mod;
  near____spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  near____spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  near____spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

  near____spec_tproc_reset_f *reset;

} near____split_generic_t;

/* sl_macro near____SPLIT_GENERIC_DEFINE_TPROC near____SPLIT_GENERIC_INIT_TPROC near____SPLIT_GENERIC_INIT_EXT_TPROC */
#define near____SPLIT_GENERIC_DEFINE_TPROC(_tp_, _s_...)         near____SPEC_DEFINE_TPROC(_tp_, _tp_, _s_)
#define near____SPLIT_GENERIC_INIT_TPROC(_tp_, _r_...)           { 1, _tp_, near____SPEC_EXT_PARAM_TPROC_NULL,  NULL, near____SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, near____SPEC_EXT_PARAM_TPROCS_NULL, NULL, near____SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define near____SPLIT_GENERIC_INIT_EXT_TPROC(_tp_, _r_...)       { 1, _tp_, near____SPEC_EXT_PARAM_TPROC(_tp_), NULL, near____SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, near____SPEC_EXT_PARAM_TPROCS_NULL, NULL, near____SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro near____SPLIT_GENERIC_DEFINE_TPROC_MOD near____SPLIT_GENERIC_INIT_TPROC_MOD near____SPLIT_GENERIC_INIT_EXT_TPROC_MOD */
#define near____SPLIT_GENERIC_DEFINE_TPROC_MOD(_tp_, _s_...)     near____SPEC_DEFINE_TPROC_MOD(_tp_, _tp_, _s_)
#define near____SPLIT_GENERIC_INIT_TPROC_MOD(_tp_, _r_...)       { 2, NULL, near____SPEC_EXT_PARAM_TPROC_NULL, _tp_, near____SPEC_EXT_PARAM_TPROC_MOD_NULL,  NULL, near____SPEC_EXT_PARAM_TPROCS_NULL, NULL, near____SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define near____SPLIT_GENERIC_INIT_EXT_TPROC_MOD(_tp_, _r_...)   { 2, NULL, near____SPEC_EXT_PARAM_TPROC_NULL, _tp_, near____SPEC_EXT_PARAM_TPROC_MOD(_tp_), NULL, near____SPEC_EXT_PARAM_TPROCS_NULL, NULL, near____SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro near____SPLIT_GENERIC_DEFINE_TPROCS near____SPLIT_GENERIC_INIT_TPROCS near____SPLIT_GENERIC_INIT_EXT_TPROCS */
#define near____SPLIT_GENERIC_DEFINE_TPROCS(_tp_, _s_...)        near____SPEC_DEFINE_TPROCS(_tp_, _tp_, _s_)
#define near____SPLIT_GENERIC_INIT_TPROCS(_tp_, _r_...)          { 3, NULL, near____SPEC_EXT_PARAM_TPROC_NULL, NULL, near____SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, near____SPEC_EXT_PARAM_TPROCS_NULL,  NULL, near____SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define near____SPLIT_GENERIC_INIT_EXT_TPROCS(_tp_, _r_...)      { 3, NULL, near____SPEC_EXT_PARAM_TPROC_NULL, NULL, near____SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, near____SPEC_EXT_PARAM_TPROCS(_tp_), NULL, near____SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro near____SPLIT_GENERIC_DEFINE_TPROCS_MOD near____SPLIT_GENERIC_INIT_TPROCS_MOD near____SPLIT_GENERIC_INIT_EXT_TPROCS_MOD */
#define near____SPLIT_GENERIC_DEFINE_TPROCS_MOD(_tp_, _s_...)    near____SPEC_DEFINE_TPROCS_MOD(_tp_, _tp_, _s_)
#define near____SPLIT_GENERIC_INIT_TPROCS_MOD(_tp_, _r_...)      { 4, NULL, near____SPEC_EXT_PARAM_TPROC_NULL, NULL, near____SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, near____SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, near____SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define near____SPLIT_GENERIC_INIT_EXT_TPROCS_MOD(_tp_, _r_...)  { 4, NULL, near____SPEC_EXT_PARAM_TPROC_NULL, NULL, near____SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, near____SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, near____SPEC_EXT_PARAM_TPROCS_MOD(_tp_), _r_ }

/* sl_type near____tloc_f near____tloc_mod_f */
typedef near____slint_t near____tloc_f(near____elements_t *b, near____slint_t x, void *tloc_data);
typedef near____slint_t near____tloc_mod_f(near____elements_t *b, near____slint_t x, void *tloc_data, near____elements_t *mod);

/* sl_type near____tproc_f near____tproc_mod_f near____tprocs_f near____tprocs_mod_f */
typedef int near____tproc_f(near____elements_t *b, near____slint_t x, void *tproc_data);
typedef int near____tproc_mod_f(near____elements_t *b, near____slint_t x, void *tproc_data, near____elements_t *mod);
typedef near____slint_t near____tprocs_f(near____elements_t *b, near____slint_t x, void *tproc_data, int *procs);
typedef near____slint_t near____tprocs_mod_f(near____elements_t *b, near____slint_t x, void *tproc_data, int *procs, near____elements_t *mod);

/* sl_type near____tproc_reset_f */
typedef void near____tproc_reset_f(void *tproc_data);

/* sl_macro near____TPROC_RESET_NULL */
#define near____TPROC_RESET_NULL  NULL

/* sl_type near_____tproc_t near____tproc_t */
typedef struct near_____tproc_t *near____tproc_t;

/* sl_type near_____tproc_exdef near____tproc_exdef */
typedef struct near_____tproc_exdef {
  int type;

  near____spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  near____spec_tproc_rearrange_db_f *tproc_rearrange_db;
  near____spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  near____spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  near____spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  near____spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  near____spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  near____spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  near____spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  near____spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  near____spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  near____spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

} const *near____tproc_exdef;

/* sl_macro near____TPROC_EXDEF_NULL */
#define near____TPROC_EXDEF_NULL  NULL

/* sl_macro near____TPROC_EXDEF_DEFINE_TPROC near____TPROC_EXDEF_DEFINE_TPROC_MOD near____TPROC_EXDEF_DEFINE_TPROCS near____TPROC_EXDEF_DEFINE_TPROCS_MOD */
#define near____TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  near____SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct near_____tproc_exdef _##_name_ = { 1, near____SPEC_EXT_PARAM_TPROC(_name_), near____SPEC_EXT_PARAM_TPROC_MOD_NULL, near____SPEC_EXT_PARAM_TPROCS_NULL, near____SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define near____TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  near____SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct near_____tproc_exdef _##_name_ = { 2, near____SPEC_EXT_PARAM_TPROC_NULL, near____SPEC_EXT_PARAM_TPROC_MOD(_name_), near____SPEC_EXT_PARAM_TPROCS_NULL, near____SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define near____TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  near____SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct near_____tproc_exdef _##_name_ = { 3, near____SPEC_EXT_PARAM_TPROC_NULL, near____SPEC_EXT_PARAM_TPROC_MOD_NULL, near____SPEC_EXT_PARAM_TPROCS(_name_), near____SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define near____TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  near____SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct near_____tproc_exdef _##_name_ = { 4, near____SPEC_EXT_PARAM_TPROC_NULL, near____SPEC_EXT_PARAM_TPROC_MOD_NULL, near____SPEC_EXT_PARAM_TPROCS_NULL, near____SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;


/* deprecated, sl_type near____k2c_func near____pivot_func near____sn_func near____m2x_func near____m2X_func */
typedef near____key2class_f near____k2c_func;
typedef near____pivot_f near____pivot_func;
typedef near____sortnet_f near____sn_func;
typedef near____merge2x_f near____m2x_func;
typedef near____merge2X_f near____m2X_func;


/* sl_type near_____mergek_t near____mergek_t */
typedef struct near_____mergek_t
{
  near____sortnet_f sn;
  near____sortnet_data_t snd;

  near____merge2x_f m2x;
  near____elements_t *sx;

} near____mergek_t;


/* sl_type near____keys_init_type_t near____keys_init_data_t */
typedef near____slint_t near____keys_init_type_t;
typedef void *near____keys_init_data_t;

/* sl_type near____key_set_data_t near____key_set_f */
typedef void *near____key_set_data_t;
typedef void (*near____key_set_f)(near____slkey_pure_t *k, near____key_set_data_t d);


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


/* near____elements_keys_stats */
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


/* partition conditions, sl_type near_____partcond2_t near____partcond2_t */
typedef struct near_____partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} near____partcond2_t;


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

/* partition conditions, sl_type near_____partcond_t near____partcond_t near____partcond_p */
typedef struct near_____partcond_t
{
  near____slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} near____partcond_t, *near____partcond_p;


/* internal partition conditions, sl_type near_____partcond_intern_t near____partcond_intern_t near____partcond_intern_p */
typedef struct near_____partcond_intern_t
{
  near____slint_t pcm;
  near____slint_t count_min, count_max;
  near____slint_t count_low, count_high;
#ifdef elem_weight
  near____slweight_t weight_min, weight_max;
  near____slweight_t weight_low, weight_high;
#endif

} near____partcond_intern_t, *near____partcond_intern_p;


/* sl_type near_____parttype_t near____parttype_t near____parttype_p */
typedef struct near_____parttype_t
{
  near____slint_t type;

} near____parttype_t, *near____parttype_p;


/* generic binning method */

/* sl_type near_____bin_t near____bin_t */
typedef struct near_____bin_t
{
  near____elements_t s;

#ifdef elem_weight
  near____slweight_t weight;
#endif

} near____bin_t;


/* sl_type near_____splitter_t near____splitter_t */
typedef struct near_____splitter_t
{
  near____slint_t n;

  int *displs;
  near____slkey_pure_t *s;
  near____slint_t *sn;

} near____splitter_t;


struct near_____binning_t;

/* sl_type near____binning_pre_f near____binning_exec_f near____binning_refine_f near____binning_hit_f near____binning_finalize_f near____binning_post_f */
typedef near____slint_t (*near____binning_pre_f)(struct near_____binning_t *bm);
typedef near____slint_t (*near____binning_exec_f)(struct near_____binning_t *bm, near____bin_t *bin, near____slcount_t *counts, near____slweight_t *weights);
typedef near____slint_t (*near____binning_refine_f)(struct near_____binning_t *bm, near____bin_t *bin, near____slint_t k, near____slcount_t *counts, near____slweight_t *weights, near____splitter_t *sp, near____slint_t s, near____bin_t *new_bin);
typedef near____slint_t (*near____binning_hit_f)(struct near_____binning_t *bm, near____bin_t *bin, near____slint_t k, near____slcount_t *counts, near____splitter_t *sp, near____slint_t s);
typedef near____slint_t (*near____binning_finalize_f)(struct near_____binning_t *bm, near____bin_t *bin, near____slint_t dc, near____slweight_t dw, near____slint_t lc_min, near____slint_t lc_max, near____slcount_t *lcs, near____slweight_t *lws, near____splitter_t *sp, near____slint_t s);
typedef near____slint_t (*near____binning_post_f)(struct near_____binning_t *bm);


/* sl_type near_____binning_data_t near____binning_data_t */
typedef union near_____binning_data_t
{
  struct
  {
    near____slint_t rhigh, rlow, rwidth;
    near____slint_t rcurrent;
    near____slkey_pure_t bit_mask;

    near____elements_t sx;

  } radix;

} near____binning_data_t;


/* sl_type near_____binning_t near____binning_t */
typedef struct near_____binning_t
{
  near____slint_t nbins, max_nbins;
  
  near____binning_pre_f pre;
  near____binning_exec_f exec;
  near____binning_refine_f refine;
  near____binning_hit_f hit;
  near____binning_finalize_f finalize;
  near____binning_post_f post;

  near____slint_t sorted;

  near____slint_t docounts;
#ifdef elem_weight
  near____slint_t doweights;
#endif

  near____binning_data_t bd;

} near____binning_t;


/* sl_type near_____local_bins_t near____local_bins_t */
typedef struct near_____local_bins_t
{
  near____binning_t *bm;

  near____slint_t nbins, max_nbins;
  near____slint_t nelements;

  near____slint_t docounts;
#ifdef elem_weight
  near____slint_t doweights;
#endif

  near____slint_t nbinnings, max_nbinnings;

  near____slint_t nbins_new, last_new_b, last_new_k;
  near____bin_t *bins, *bins_new;
  near____bin_t *bins0, *bins1;

  near____slint_t *bcws;

#if defined(elem_weight) && defined(near____sl_weight_intequiv)
  near____slint_t cw_factor, w_index, bin_cw_factor;
  near____slweight_t *cws, *bin_cws;
  near____slweight_t *prefix_cws;
#else
  near____slint_t *cs, *bin_cs;
  near____slint_t *prefix_cs;
# ifdef elem_weight
  near____slweight_t *ws, *bin_ws;
  near____slweight_t *prefix_ws;
# endif
#endif

  near____slint_t last_exec_b;

} near____local_bins_t;


/* sl_type near_____global_bins_t near____global_bins_t */
typedef struct near_____global_bins_t
{
  near____binning_t *bm;
  
  near____local_bins_t lb;

  near____slint_t *bcws;

#if defined(elem_weight) && defined(near____sl_weight_intequiv)
  near____slweight_t *cws;
  near____slweight_t *prefix_cws;
#else
  near____slint_t *cs;
  near____slint_t *prefix_cs;
# ifdef elem_weight
  near____slweight_t *ws;
  near____slweight_t *prefix_ws;
# endif
#endif

} near____global_bins_t;


/* sl_type near____rti_cmc_t */
typedef struct
{
  near____slint_t cmp, movek, moved;

} near____rti_cmc_t;

#ifndef my_rti_ccmp
# define my_rti_ccmp(m)    m.cmc.cmp
# define my_rti_cmovek(m)  m.cmc.movek
# define my_rti_cmoved(m)  m.cmc.moved
#endif


/* sl_type near____rti_tim_t */
typedef struct
{
  double start, stop;
  double last, cumu;

  near____slint_t num;

} near____rti_tim_t[rti_tids];

#ifndef my_rti_tlast
# define my_rti_tlast(m, t)  m.tim[t].last
# define my_rti_tcumu(m, t)  m.tim[t].cumu
# define my_rti_tnum(m, t)   m.tim[t].num
#endif


/* sl_type near____rti_mem_t */
typedef struct
{
  near____slint_t nalloc, nfree;
  near____slint_t max, cur, cur_max;

} near____rti_mem_t;


/* sl_type near____rti_t */
typedef struct
{
  /* compare-move-counter */
  near____rti_cmc_t cmc;
  /* timer */
  near____rti_tim_t tim;
  /* memory */
  near____rti_mem_t mem;

} near____rti_t;

#ifndef my_rti_reset
# define my_rti_reset(m)  memset((void *) &m, 0, sizeof(m))
#endif


/* sl_type near_____sl_context_t near____sl_context_t */
typedef struct near_____sl_context_t
{

/* src/base/base.c */
  struct {
int dummy_rank;
  } sl;
#ifdef near____SL_USE_RTI
near____rti_t rti;
#endif
  struct {
near____slint_t ip_threshold;
near____slint_t db_threshold;
near____slint_t ma_threshold;
  } sr;
  struct {
near____slint_t threshold;
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
near____slint_t sendrecv_replace_memsize;
near____slint_t sendrecv_replace_mpi_maxsize;
  } me;
#endif
#ifdef SL_USE_MPI
  struct {
double t[2];
near____slint_t max_nprocs;
near____slint_t packed;
near____slint_t minalloc;
double overalloc;
  } meas;
#endif
#ifdef SL_USE_MPI
  struct {
near____slint_t packed;
near____slint_t db_packed;
near____slint_t ip_packed;
  } mea;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef near____MSEG_ROOT
int root;
#endif
#ifdef near____MSEG_BORDER_UPDATE_REDUCTION
double border_update_count_reduction;
double border_update_weight_reduction;
#endif
#ifdef near____MSEG_FORWARD_ONLY
near____slint_t forward_only;
#endif
#ifdef near____MSEG_INFO
near____slint_t info_rounds;
near____slint_t *info_finish_rounds;
double info_finish_rounds_avg;
near____slweight_t info_total_weights;
#endif
near____slint_t binnings;
near____slint_t finalize_mode;
  } mseg;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef near____MSS_ROOT
int root;
#endif
  } mss;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
near____slint_t sync;
  } msm;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
near____slint_t sync;
near____partcond_t *r_pc;
  } msp;
#endif
#ifdef SL_USE_MPI
  struct {
double i_t[3];
double p_t[3];
double b_t[3];
near____slint_t sync;
near____slint_t i_sync;
near____slint_t p_sync;
near____slint_t b_sync;
near____slint_t back_packed;
  } mssp;
#endif
} near____sl_context_t;






/* sl_macro near____elem_set_size near____elem_set_max_size near____elem_set_keys near____elem_set_indices */
#define near____elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define near____elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define near____elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define near____elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))
/* DATAX_TEMPLATE_BEGIN */
#define near____elem_set_data0(_e_, _b_)     ((_e_)->data0 = (_b_))  /* sl_macro */
#define near____elem_set_data1(_e_, _b_)     ((_e_)->data1 = (_b_))  /* sl_macro */
#define near____elem_set_data2(_e_, _b_)     ((_e_)->data2 = (_b_))  /* sl_macro */
#define near____elem_set_data3(_e_, _b_)     ((_e_)->data3 = (_b_))  /* sl_macro */
#define near____elem_set_data4(_e_, _b_)     ((_e_)->data4 = (_b_))  /* sl_macro */
#define near____elem_set_data5(_e_, _b_)     ((_e_)->data5 = (_b_))  /* sl_macro */
#define near____elem_set_data6(_e_, _b_)     ((_e_)->data6 = (_b_))  /* sl_macro */
#define near____elem_set_data7(_e_, _b_)     ((_e_)->data7 = (_b_))  /* sl_macro */
#define near____elem_set_data8(_e_, _b_)     ((_e_)->data8 = (_b_))  /* sl_macro */
#define near____elem_set_data9(_e_, _b_)     ((_e_)->data9 = (_b_))  /* sl_macro */
#define near____elem_set_data10(_e_, _b_)     ((_e_)->data10 = (_b_))  /* sl_macro */
#define near____elem_set_data11(_e_, _b_)     ((_e_)->data11 = (_b_))  /* sl_macro */
#define near____elem_set_data12(_e_, _b_)     ((_e_)->data12 = (_b_))  /* sl_macro */
#define near____elem_set_data13(_e_, _b_)     ((_e_)->data13 = (_b_))  /* sl_macro */
#define near____elem_set_data14(_e_, _b_)     ((_e_)->data14 = (_b_))  /* sl_macro */
#define near____elem_set_data15(_e_, _b_)     ((_e_)->data15 = (_b_))  /* sl_macro */
#define near____elem_set_data16(_e_, _b_)     ((_e_)->data16 = (_b_))  /* sl_macro */
#define near____elem_set_data17(_e_, _b_)     ((_e_)->data17 = (_b_))  /* sl_macro */
#define near____elem_set_data18(_e_, _b_)     ((_e_)->data18 = (_b_))  /* sl_macro */
#define near____elem_set_data19(_e_, _b_)     ((_e_)->data19 = (_b_))  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro near____elem_get_size near____elem_get_max_size near____elem_get_keys near____elem_get_indices */
#define near____elem_get_size(_e_)           (_e_)->size
#define near____elem_get_max_size(_e_)       (_e_)->max_size
#define near____elem_get_keys(_e_)           (_e_)->keys
#define near____elem_get_indices(_e_)        (_e_)->indices
/* DATAX_TEMPLATE_BEGIN */
#define near____elem_get_data0(_e_)          (_e_)->data0  /* sl_macro */
#define near____elem_get_data1(_e_)          (_e_)->data1  /* sl_macro */
#define near____elem_get_data2(_e_)          (_e_)->data2  /* sl_macro */
#define near____elem_get_data3(_e_)          (_e_)->data3  /* sl_macro */
#define near____elem_get_data4(_e_)          (_e_)->data4  /* sl_macro */
#define near____elem_get_data5(_e_)          (_e_)->data5  /* sl_macro */
#define near____elem_get_data6(_e_)          (_e_)->data6  /* sl_macro */
#define near____elem_get_data7(_e_)          (_e_)->data7  /* sl_macro */
#define near____elem_get_data8(_e_)          (_e_)->data8  /* sl_macro */
#define near____elem_get_data9(_e_)          (_e_)->data9  /* sl_macro */
#define near____elem_get_data10(_e_)          (_e_)->data10  /* sl_macro */
#define near____elem_get_data11(_e_)          (_e_)->data11  /* sl_macro */
#define near____elem_get_data12(_e_)          (_e_)->data12  /* sl_macro */
#define near____elem_get_data13(_e_)          (_e_)->data13  /* sl_macro */
#define near____elem_get_data14(_e_)          (_e_)->data14  /* sl_macro */
#define near____elem_get_data15(_e_)          (_e_)->data15  /* sl_macro */
#define near____elem_get_data16(_e_)          (_e_)->data16  /* sl_macro */
#define near____elem_get_data17(_e_)          (_e_)->data17  /* sl_macro */
#define near____elem_get_data18(_e_)          (_e_)->data18  /* sl_macro */
#define near____elem_get_data19(_e_)          (_e_)->data19  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro near____elem_set_block near____elem_set_block_size near____elem_get_block near____elem_get_block_size */
#define near____elem_set_block(_e_, _b_)       ((_e_)->keys = (_b_), (_e_)->max_size = -1)
#define near____elem_set_block_size(_e_, _s_)  ((_e_)->size = (_s_))
#define near____elem_get_block(_e_)            ((void *) (((_e_)->max_size < 0)?(_e_)->keys:NULL))
#define near____elem_get_block_size(_e_)       (_e_)->size

/* sl_macro near____pelem_set_size near____pelem_set_max_size near____pelem_set_elements */
#define near____pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define near____pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define near____pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))

/* sl_macro near____pelem_get_size near____pelem_get_max_size near____pelem_get_elements */
#define near____pelem_get_size(_e_)           (_e_)->size
#define near____pelem_get_max_size(_e_)       (_e_)->max_size
#define near____pelem_get_elements(_e_)       (_e_)->elements

/* sl_macro near____SL_DEFCON */
#define near____SL_DEFCON(_v_)  (near____sl_default_context._v_)






/* src/base/base.c */
extern near____sl_context_t near____sl_default_context;
extern const int near____default_sl_dummy_rank;
#ifdef near____SL_USE_RTI
extern const near____rti_t near____default_rti;
#endif
extern const near____slint_t near____default_sr_ip_threshold;
extern const near____slint_t near____default_sr_db_threshold;
extern const near____slint_t near____default_sr_ma_threshold;
extern const near____slint_t near____default_sri_threshold;

/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
extern const MPI_Datatype near____default_mpi_int_datatype;
extern const MPI_Datatype near____default_mpi_key_datatype;
extern const MPI_Datatype near____default_mpi_pkey_datatype;
extern const MPI_Datatype near____default_mpi_pwkey_datatype;
extern const MPI_Datatype near____default_mpi_index_datatype;
extern const MPI_Datatype near____default_mpi_weight_datatype;
extern const MPI_Datatype near____default_mpi_data_datatype[];
extern const int near____default_mpi_rank;
#endif
extern const void *near____default_me_sendrecv_replace_mem;
extern const near____slint_t near____default_me_sendrecv_replace_memsize;
extern const near____slint_t near____default_me_sendrecv_replace_mpi_maxsize;
extern const double near____default_meas_t[];
extern const near____slint_t near____default_meas_max_nprocs;
extern const near____slint_t near____default_meas_packed;
extern const near____slint_t near____default_meas_minalloc;
extern const double near____default_meas_overalloc;
extern const near____slint_t near____default_mea_packed;
extern const near____slint_t near____default_mea_db_packed;
extern const near____slint_t near____default_mea_ip_packed;
#ifdef near____MSEG_ROOT
extern const int near____default_mseg_root;
#endif
#ifdef near____MSEG_BORDER_UPDATE_REDUCTION
extern const double near____default_mseg_border_update_count_reduction;
extern const double near____default_mseg_border_update_weight_reduction;
#endif
#ifdef near____MSEG_FORWARD_ONLY
extern const near____slint_t near____default_mseg_forward_only;
#endif
#ifdef near____MSEG_INFO
extern const near____slint_t near____default_mseg_info_rounds;
extern const near____slint_t *near____default_mseg_info_finish_rounds;
extern const double near____default_mseg_info_finish_rounds_avg;
extern const near____slweight_t near____default_mseg_info_total_weights;
#endif
extern const near____slint_t near____default_mseg_binnings;
extern const near____slint_t near____default_mseg_finalize_mode;
#ifdef near____MSS_ROOT
extern const int near____default_mss_root;
#endif
extern const double near____default_msm_t[];
extern const near____slint_t near____default_msm_sync;
extern const double near____default_msp_t[];
extern const near____slint_t near____default_msp_sync;
extern const near____partcond_t *near____default_msp_r_pc;
extern const double near____default_mssp_i_t[];
extern const double near____default_mssp_p_t[];
extern const double near____default_mssp_b_t[];
extern const near____slint_t near____default_mssp_sync;
extern const near____slint_t near____default_mssp_i_sync;
extern const near____slint_t near____default_mssp_p_sync;
extern const near____slint_t near____default_mssp_b_sync;
extern const near____slint_t near____default_mssp_back_packed;






/* src/base/base.c */
near____slint_t SL_PROTO(near____binning_create)(near____local_bins_t *lb, near____slint_t max_nbins, near____slint_t max_nbinnings, near____elements_t *s, near____slint_t nelements, near____slint_t docounts, near____slint_t doweights, near____binning_t *bm);
near____slint_t SL_PROTO(near____binning_destroy)(near____local_bins_t *lb);
near____slint_t SL_PROTO(near____binning_pre)(near____local_bins_t *lb);
near____slint_t SL_PROTO(near____binning_exec_reset)(near____local_bins_t *lb, near____slint_t do_bins, near____slint_t do_prefixes);
near____slint_t SL_PROTO(near____binning_exec)(near____local_bins_t *lb, near____slint_t b, near____slint_t do_bins, near____slint_t do_prefixes);
near____slint_t SL_PROTO(near____binning_refine)(near____local_bins_t *lb, near____slint_t b, near____slint_t k, near____splitter_t *sp, near____slint_t s);
near____slint_t SL_PROTO(near____binning_hit)(near____local_bins_t *lb, near____slint_t b, near____slint_t k, near____splitter_t *sp, near____slint_t s);
near____slint_t SL_PROTO(near____binning_finalize)(near____local_bins_t *lb, near____slint_t b, near____slint_t dc, near____slweight_t dw, near____slint_t lc_min, near____slint_t lc_max, near____slcount_t *lcs, near____slweight_t *lws, near____splitter_t *sp, near____slint_t s);
near____slint_t SL_PROTO(near____binning_post)(near____local_bins_t *lb);
near____slint_t SL_PROTO(near____binning_radix_create)(near____binning_t *bm, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, near____slint_t sorted);
near____slint_t SL_PROTO(near____binning_radix_destroy)(near____binning_t *bm);
near____slint_t SL_PROTO(near____binning_radix_pre)(near____binning_t *bm);
near____slint_t SL_PROTO(near____binning_radix_exec)(near____binning_t *bm, near____bin_t *bin, near____slcount_t *counts, near____slweight_t *weights);
near____slint_t SL_PROTO(near____binning_radix_refine)(near____binning_t *bm, near____bin_t *bin, near____slint_t k, near____slcount_t *counts, near____slweight_t *weights, near____splitter_t *sp, near____slint_t s, near____bin_t *new_bin);
near____slint_t SL_PROTO(near____binning_radix_hit)(near____binning_t *bm, near____bin_t *bin, near____slint_t k, near____slcount_t *counts, near____splitter_t *sp, near____slint_t s);
near____slint_t SL_PROTO(near____binning_radix_finalize)(near____binning_t *bm, near____bin_t *bin, near____slint_t dc, near____slweight_t dw, near____slint_t lc_min, near____slint_t lc_max, near____slcount_t *lcs, near____slweight_t *lws, near____splitter_t *sp, near____slint_t s);
near____slint_t SL_PROTO(near____binning_radix_post)(near____binning_t *bm);
near____slint_t SL_PROTO(near____elements_alloc)(near____elements_t *s, near____slint_t nelements, slcint_t components);
near____slint_t SL_PROTO(near____elements_free)(near____elements_t *s);
near____slint_t SL_PROTO(near____elements_realloc)(near____elements_t *s, near____slint_t nelements, slcint_t components);
near____slint_t SL_PROTO(near____elements_alloca)(near____elements_t *s, near____slint_t nelements, slcint_t components);
near____slint_t SL_PROTO(near____elements_freea)(near____elements_t *s);
near____slint_t SL_PROTO(near____elements_alloc_from_blocks)(near____elements_t *s, near____slint_t nblocks, void **blocks, near____slint_t *blocksizes, near____slint_t alignment, near____slint_t nmax, slcint_t components);
near____slint_t SL_PROTO(near____elements_alloc_from_block)(near____elements_t *s, void *block, near____slint_t blocksize, near____slint_t alignment, near____slint_t nmax, slcint_t components);
near____slint_t SL_PROTO(near____elements_alloc_block)(near____elements_t *s, void **block, near____slint_t *blocksize, near____slint_t alignment, near____slint_t maxblocksize);
near____slint_t SL_PROTO(near____elements_copy)(near____elements_t *s, near____elements_t *d);
near____slint_t SL_PROTO(near____elements_copy_at)(near____elements_t *s, near____slint_t sat, near____elements_t *d, near____slint_t dat);
near____slint_t SL_PROTO(near____elements_ncopy)(near____elements_t *s, near____elements_t *d, near____slint_t n);
near____slint_t SL_PROTO(near____elements_nmove)(near____elements_t *s, near____elements_t *d, near____slint_t n);
near____slint_t SL_PROTO(near____elements_printf)(near____elements_t *s, const char *prefix);
near____slint_t SL_PROTO(near____elements_extract)(near____elements_t *src, near____slint_t nelements, near____elements_t *dst0, near____elements_t *dst1);
near____slint_t SL_PROTO(near____elements_touch)(near____elements_t *s);
near____slint_t SL_PROTO(near____elements_digest_sum)(near____elements_t *s, near____slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(near____elements_crc32)(near____elements_t *s, near____slint nelements, near____slint_t keys, near____slint_t data);
near____slint_t SL_PROTO(near____elements_digest_hash)(near____elements_t *s, near____slint_t nelements, slcint_t components, void *hash);
near____slint_t SL_PROTO(near____elements_random_exchange)(near____elements_t *s, near____slint_t rounds, near____elements_t *xs);
near____slint_t SL_PROTO(near____elements_keys_init_seed)(unsigned long s);
near____slint_t SL_PROTO(near____elements_keys_init)(near____elements_t *s, near____keys_init_type_t t, near____keys_init_data_t d);
near____slint_t SL_PROTO(near____elements_keys_init_randomized)(near____elements_t *s, near____slint_t nkeys, near____keys_init_type_t t, near____keys_init_data_t d);
near____slint_t SL_PROTO(near____elements_keys_init_from_file)(near____elements_t *s, near____slint_t data, char *filename, near____slint_t from, near____slint_t to, near____slint_t const_bytes_per_line);
near____slint_t SL_PROTO(near____elements_keys_save_to_file)(near____elements_t *s, char *filename);
near____slint_t SL_PROTO(near____elements_validate_order)(near____elements_t *s, near____slint_t n);
near____slint_t SL_PROTO(near____elements_validate_order_bmask)(near____elements_t *s, near____slint_t n, near____slkey_pure_t bmask);
near____slint_t SL_PROTO(near____elements_validate_order_weight)(near____elements_t *s, near____slint_t n, near____slkey_pure_t weight);
near____slint_t SL_PROTO(near____elements_keys_stats)(near____elements_t *s, near____slkey_pure_t *stats);
near____slint_t SL_PROTO(near____elements_keys_stats_print)(near____elements_t *s);
near____slint_t SL_PROTO(near____elements_print_keys)(near____elements_t *s);
near____slint_t SL_PROTO(near____elements_print_all)(near____elements_t *s);
near____slweight_t SL_PROTO(near____elements_get_weight)(near____elements_t *s);
near____slint_t SL_PROTO(near____elements_get_minmax_keys)(near____elements_t *s, near____slint_t nelements, near____slkey_pure_t *minmaxkeys);
near____slint_t SL_PROTO(near____elements_alloc_packed)(near____packed_elements_t *s, near____slint_t nelements);
near____slint_t SL_PROTO(near____elements_free_packed)(near____packed_elements_t *s);
near____slint_t SL_PROTO(near____elements_alloc_packed_from_block)(near____packed_elements_t *s, void *block, near____slint_t blocksize, near____slint_t alignment, near____slint_t nmax);
near____slint_t SL_PROTO(near____elements_pack_indexed)(near____elements_t *s, near____packed_elements_t *d, near____slindex_t *rindx, near____slindex_t *windx);
near____slint_t SL_PROTO(near____elements_pack)(near____elements_t *s, near____packed_elements_t *d);
near____slint_t SL_PROTO(near____elements_pack_at)(near____elements_t *s, near____slint_t sat, near____packed_elements_t *d, near____slint_t dat);
near____slint_t SL_PROTO(near____elements_unpack_indexed)(near____packed_elements_t *s, near____elements_t *d, near____slindex_t *rindx, near____slindex_t *windx);
near____slint_t SL_PROTO(near____elements_unpack)(near____packed_elements_t *s, near____elements_t *d);
near____slint_t SL_PROTO(near____elements_unpack_at)(near____packed_elements_t *s, near____slint_t sat, near____elements_t *d, near____slint_t dat);
near____slint_t SL_PROTO(near____elements_unpack_keys)(near____packed_elements_t *s, near____slkey_t *k);
near____slint SL_PROTO(near____merge2_basic_auto_01_x)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint SL_PROTO(near____merge2_basic_01_x)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx, near____m2x_func _x0_1, near____m2x_func _0x_1);
near____slint SL_PROTO(near____merge2_basic_01_X)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx, near____elements_t *t, near____m2X_func _X0_1, near____m2X_func _0X_1);
near____slint SL_PROTO(near____merge2_simplify_s1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx, near____slint s1elements);
near____slint SL_PROTO(near____merge2_memory_adaptive)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_compo_hula)(near____elements_t *s0, near____elements_t *s1, near____elements_t *xs);
near____slint_t SL_PROTO(near____merge2_basic_sseq_x0_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_sseq_0x_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_sseq_01_x)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_sseq_01)(near____elements_t *s0, near____elements_t *s1, near____elements_t *t);
near____slint_t SL_PROTO(near____merge2_basic_sbin_x0_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_sbin_0x_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_sbin_01_x)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_sbin_01)(near____elements_t *s0, near____elements_t *s1, near____elements_t *t);
near____slint_t SL_PROTO(near____merge2_basic_shyb_x0_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_shyb_0x_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_shyb_01_x)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_shyb_01)(near____elements_t *s0, near____elements_t *s1, near____elements_t *t);
near____slint_t SL_PROTO(near____merge2_basic_straight_x0_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_straight_0x_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_straight_01_x)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_straight_x_0_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____merge2_basic_straight_X0_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx, near____elements_t *t);
near____slint_t SL_PROTO(near____merge2_basic_straight_0X_1)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx, near____elements_t *t);
near____slint_t SL_PROTO(near____merge2_basic_straight_01_X)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx, near____elements_t *t);
near____slint_t SL_PROTO(near____merge2_basic_straight_X0_1u)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx, near____elements_t *t);
near____slint_t SL_PROTO(near____merge2_compo_tridgell)(near____elements_t *s0, near____elements_t *s1, near____elements_t *sx);
near____slint_t SL_PROTO(near____mergep_2way_ip_int)(near____elements_t *s, near____elements_t *sx, near____slint_t p, int *displs, near____merge2x_f m2x);
near____slint_t SL_PROTO(near____mergep_2way_ip_int_rec)(near____elements_t *s, near____elements_t *sx, near____slint_t p, int *displs, near____merge2x_f m2x);
near____slint_t SL_PROTO(near____mergep_heap_int)(near____elements_t *s, near____elements_t *d, near____slint_t p, int *displs, int *counts);
near____slint_t SL_PROTO(near____mergep_heap_int_idx)(near____elements_t *s, near____elements_t *d, near____slint_t p, int *displs, int *counts);
near____slint_t SL_PROTO(near____mergep_heap_idx)(near____elements_t *s, near____elements_t *d, near____slint_t p, near____slindex_t *displs, near____slindex_t *counts);
near____slint_t SL_PROTO(near____mergep_heap_unpack_idx)(near____packed_elements_t *s, near____elements_t *d, near____slint_t p, near____slindex_t *displs, near____slindex_t *counts);
near____slint_t SL_PROTO(near____mergep_heap_unpack_idxonly)(near____packed_elements_t *s, near____elements_t *d, near____slint_t p, near____slindex_t *displs, near____slindex_t *counts);
near____slint_t SL_PROTO(near____permute_generic_db)(near____elements_t *s, near____elements_t *d, near____permute_generic_t *pg, void *pg_data);
near____slint_t SL_PROTO(near____permute_generic_ip)(near____elements_t *s, near____elements_t *x, near____permute_generic_t *pg, void *pg_data);
near____slint SL_PROTO(near____sl_search_sequential_lt)(near____elements_t *s, near____slpkey_t k);
near____slint SL_PROTO(near____sl_search_sequential_le)(near____elements_t *s, near____slpkey_t k);
near____slint SL_PROTO(near____sl_search_sequential_gt)(near____elements_t *s, near____slpkey_t k);
near____slint SL_PROTO(near____sl_search_sequential_ge)(near____elements_t *s, near____slpkey_t k);
near____slint SL_PROTO(near____sl_search_p_sequential_lt)(near____elements_t *s, near____slpkey_t *k);
near____slint SL_PROTO(near____sl_search_p_sequential_le)(near____elements_t *s, near____slpkey_t *k);
near____slint SL_PROTO(near____sl_search_p_sequential_gt)(near____elements_t *s, near____slpkey_t *k);
near____slint SL_PROTO(near____sl_search_p_sequential_ge)(near____elements_t *s, near____slpkey_t *k);
near____slint SL_PROTO(near____sl_search_binary_lt)(near____elements_t *s, near____slpkey_t k);
near____slint SL_PROTO(near____sl_search_binary_le)(near____elements_t *s, near____slpkey_t k);
near____slint SL_PROTO(near____sl_search_binary_gt)(near____elements_t *s, near____slpkey_t k);
near____slint SL_PROTO(near____sl_search_binary_ge)(near____elements_t *s, near____slpkey_t k);
near____slint SL_PROTO(near____sl_search_p_binary_lt)(near____elements_t *s, near____slpkey_t *k);
near____slint SL_PROTO(near____sl_search_p_binary_le)(near____elements_t *s, near____slpkey_t *k);
near____slint SL_PROTO(near____sl_search_p_binary_gt)(near____elements_t *s, near____slpkey_t *k);
near____slint SL_PROTO(near____sl_search_p_binary_ge)(near____elements_t *s, near____slpkey_t *k);
near____slint_t SL_PROTO(near____sl_search_binary_lt_bmask)(near____elements_t *s, near____slpkey_t k, near____slpkey_t bmask);
near____slint_t SL_PROTO(near____sl_search_binary_le_bmask)(near____elements_t *s, near____slpkey_t k, near____slpkey_t bmask);
near____slint_t SL_PROTO(near____sl_search_binary_sign_switch)(near____elements_t *s);
near____slint SL_PROTO(near____sl_search_hybrid_lt)(near____elements_t *s, near____slpkey_t k, near____slint t);
near____slint SL_PROTO(near____sl_search_hybrid_le)(near____elements_t *s, near____slpkey_t k, near____slint t);
near____slint SL_PROTO(near____sl_search_hybrid_gt)(near____elements_t *s, near____slpkey_t k, near____slint t);
near____slint SL_PROTO(near____sl_search_hybrid_ge)(near____elements_t *s, near____slpkey_t k, near____slint t);
near____slint SL_PROTO(near____sl_search_p_hybrid_lt)(near____elements_t *s, near____slpkey_t *k, near____slint t);
near____slint SL_PROTO(near____sl_search_p_hybrid_le)(near____elements_t *s, near____slpkey_t *k, near____slint t);
near____slint SL_PROTO(near____sl_search_p_hybrid_gt)(near____elements_t *s, near____slpkey_t *k, near____slint t);
near____slint SL_PROTO(near____sl_search_p_hybrid_ge)(near____elements_t *s, near____slpkey_t *k, near____slint t);
near____slint SL_PROTO(near____ilog2c)(near____slint x);
near____slint SL_PROTO(near____ilog2f)(near____slint x);
near____slint SL_PROTO(near____print_bits)(near____slint v);
near____slint SL_PROTO(near____pivot_random)(near____elements_t *s);
near____slint_t SL_PROTO(near____counts2displs)(near____slint_t n, int *counts, int *displs);
near____slint_t SL_PROTO(near____displs2counts)(near____slint_t n, int *displs, int *counts, near____slint_t total_counts);
void SL_PROTO(near____get_displcounts_extent)(near____slint_t n, int *displs, int *counts, near____slint_t *lb, near____slint_t *extent);
void SL_PROTO(near____elem_set_data)(near____elements_t *e, ...);
near____slint_t SL_PROTO(near____elem_get_max_byte)();
near____slint_t SL_PROTO(near____elem_reverse)(near____elements_t *e, near____elements_t *t);
near____slint_t SL_PROTO(near____elem_nxchange_at)(near____elements_t *e0, near____slint_t at0, near____elements_t *e1, near____slint_t at1, near____slint_t n, near____elements_t *t);
near____slint_t SL_PROTO(near____elem_nxchange)(near____elements_t *e0, near____elements_t *e1, near____slint_t n, near____elements_t *t);
near____slint_t SL_PROTO(near____elem_nxchange_ro0)(near____elements_t *e0, near____elements_t *e1, near____slint_t n, near____elements_t *t);
near____slint_t SL_PROTO(near____elem_rotate)(near____elements_t *e, near____slint_t m, near____slint_t n, near____elements_t *t);
near____slint_t SL_PROTO(near____elem_rotate_ro0)(near____elements_t *e, near____slint_t m, near____slint_t n, near____elements_t *t);
near____slint_t SL_PROTO(near____elem_rotate_ro1)(near____elements_t *e, near____slint_t m, near____slint_t n, near____elements_t *t);
near____slint_t SL_PROTO(near____sort_counting_use_displs)(near____elements_t *s, near____elements_t *d, near____slint_t ndispls, near____slint_t *displs);
near____slint_t SL_PROTO(near____sort_counting_use_counts)(near____elements_t *s, near____elements_t *d, near____slint_t ncounts, near____slint_t *counts);
near____slint_t SL_PROTO(near____sort_counting_get_counts)(near____elements_t *s, near____elements_t *d, near____slint_t ncounts, near____slint_t *counts);
near____slint_t SL_PROTO(near____sort_counting)(near____elements_t *s, near____elements_t *d, near____slint_t ncounts);
near____slint SL_PROTO(near____sort_heap)(near____elements_t *s, near____elements_t *xs);
near____slint_t SL_PROTO(near____sort_insert_bmask_kernel)(near____elements_t *s, near____elements_t *sx, near____slkey_pure_t bmask);
near____slint_t SL_PROTO(near____sort_insert)(near____elements_t *s, near____elements_t *sx);
near____slint_t SL_PROTO(near____sort_permute_forward)(near____elements_t *s, near____elements_t *sx, near____slint_t *perm, near____slint_t offset, near____slint_t mask_bit);
near____slint_t SL_PROTO(near____sort_permute_backward)(near____elements_t *s, near____elements_t *sx, near____slint_t *perm, near____slint_t offset, near____slint_t mask_bit);
near____slint SL_PROTO(near____sort_quick)(near____elements_t *s, near____elements_t *xs);
near____slint_t SL_PROTO(near____sort_radix_ip)(near____elements_t *s, near____elements_t *sx, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth);
near____slint_t SL_PROTO(near____sort_radix_db)(near____elements_t *s, near____elements_t *sx, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth);
near____slint_t SL_PROTO(near____sort_radix_ma)(near____elements_t *s, near____elements_t *sx, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth);
near____slint_t SL_PROTO(near____sort_radix)(near____elements_t *s, near____elements_t *sx, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth);
near____slint_t SL_PROTO(near____sort_radix_1bit_kernel)(near____elements_t *s, near____elements_t *sx, near____slint_t rhigh, near____slint_t rlow);
near____slint SL_PROTO(near____sort_radix_1bit)(near____elements_t *s, near____elements_t *sx, near____slint_t rhigh, near____slint_t rlow);
near____slint_t SL_PROTO(near____sort_radix_iter)(near____elements_t *s, near____elements_t *sx, near____slint_t presorted, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth);
near____slint SL_PROTO(near____sn_hypercube_lh)(near____slint size, near____slint rank, near____slint stage, void *snp, near____slint *up);
near____slint SL_PROTO(near____sn_hypercube_hl)(near____slint size, near____slint rank, near____slint stage, void *snp, near____slint *up);
near____slint SL_PROTO(near____sn_odd_even_trans)(near____slint size, near____slint rank, near____slint stage, void *snp, near____slint *up);
near____slint SL_PROTO(near____sn_odd)(near____slint size, near____slint rank, near____slint stage, void *snp, near____slint *up);
near____slint SL_PROTO(near____sn_even)(near____slint size, near____slint rank, near____slint stage, void *snp, near____slint *up);
near____slint SL_PROTO(near____sn_batcher)(near____slint size, near____slint rank, near____slint stage, void *snp, near____slint *up);
near____slint SL_PROTO(near____sn_bitonic)(near____slint size, near____slint rank, near____slint stage, void *snp, near____slint *up);
near____slint SL_PROTO(near____sn_connected)(near____slint size, near____slint rank, near____slint stage, void *snp, near____slint *up);
near____slint_t SL_PROTO(near____split_generic_db)(near____elements_t *s, near____elements_t *d, near____split_generic_t *sg, void *sg_data, near____slint_t n);
near____slint_t SL_PROTO(near____split_generic_ip)(near____elements_t *s, near____elements_t *d, near____split_generic_t *sg, void *sg_data, near____slint_t n);
near____slint_t SL_PROTO(near____split_generic_count_db)(near____elements_t *s, near____split_generic_t *sg, void *sg_data, int *counts, near____slint_t n);
near____slint_t SL_PROTO(near____split_generic_count_ip)(near____elements_t *s, near____split_generic_t *sg, void *sg_data, int *counts, near____slint_t n);
near____slint_t SL_PROTO(near____split_generic_rearrange_db)(near____elements_t *s, near____elements_t *d, near____split_generic_t *sg, void *sg_data, int *counts, near____slint_t n);
near____slint_t SL_PROTO(near____split_generic_rearrange_ip)(near____elements_t *s, near____elements_t *d, near____split_generic_t *sg, void *sg_data, int *counts, int *displs, near____slint_t n);
near____slint_t SL_PROTO(near____splitter_reset)(near____splitter_t *sp);
near____slint_t SL_PROTO(near____splitx_radix)(near____elements_t *s, near____elements_t *sx, near____slint_t nclasses, near____slint_t shl, near____slint_t *counts);
near____slint SL_PROTO(near____split2_lt_ge)(near____elements_t *s, near____slkey_pure_t *k, near____elements_t *t);
near____slint SL_PROTO(near____split2_le_gt)(near____elements_t *s, near____slkey_pure_t *k, near____elements_t *t);
near____slint SL_PROTO(near____split3_lt_eq_gt)(near____elements_t *s, near____slkey_pure_t *k, near____elements_t *t, near____slint *nlt, near____slint *nle);
near____slint SL_PROTO(near____split3_lt_eq_gt_old)(near____elements_t *s, near____slkey_pure_t *k, near____elements_t *t, near____slint *nlt, near____slint *nle);
near____slint SL_PROTO(near____split2_b)(near____elements_t *s, near____elements_t *sx, near____slkey_pure_t bmask);
near____slint SL_PROTO(near____splitk_k2c_af)(near____elements_t *s, near____elements_t *sx, near____slint k, near____slint *c, near____k2c_func k2c, void *k2c_data);
near____slint SL_PROTO(near____splitk_k2c)(near____elements_t *s, near____elements_t *sx, near____slint k, near____slint *c, near____k2c_func k2c, void *k2c_data);
near____slint SL_PROTO(near____splitk_k2c_count)(near____elements_t *s, near____slint k, near____slint *c, near____k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI





/* src/base_mpi/base_mpi.c */
near____slint_t SL_PROTO(near____mpi_binning_create)(near____global_bins_t *gb, near____slint_t max_nbins, near____slint_t max_nbinnings, near____elements_t *s, near____slint_t nelements, near____slint_t docounts, near____slint_t doweights, near____binning_t *bm, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_binning_destroy)(near____global_bins_t *gb, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_binning_pre)(near____global_bins_t *gb, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_binning_exec_reset)(near____global_bins_t *gb, near____slint_t do_bins, near____slint_t do_prefixes, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_binning_exec_local)(near____global_bins_t *gb, near____slint_t b, near____slint_t do_bins, near____slint_t do_prefixes, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_binning_exec_global)(near____global_bins_t *gb, near____slint_t do_bins, near____slint_t do_prefixes, near____slint_t root, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_binning_refine)(near____global_bins_t *gb, near____slint_t b, near____slint_t k, near____splitter_t *sp, near____slint_t s, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_binning_hit)(near____global_bins_t *gb, near____slint_t b, near____slint_t k, near____splitter_t *sp, near____slint_t s, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_binning_finalize)(near____global_bins_t *gb, near____slint_t b, near____slint_t dc, near____slweight_t dw, near____slint_t lc_min, near____slint_t lc_max, near____slcount_t *lcs, near____slweight_t *lws, near____splitter_t *sp, near____slint_t s, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_binning_post)(near____global_bins_t *gb, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_datatypes_init)();
near____slint_t SL_PROTO(near____mpi_datatypes_release)();
near____slint_t SL_PROTO(near____mpi_get_grid_properties)(near____slint_t ndims, near____slint_t *dims, near____slint_t *pos, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_subgroups_create)(near____slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_subgroups_delete)(near____slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);
int SL_PROTO(near____sl_MPI_Allreduce)(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank);
int SL_PROTO(near____sl_MPI_Alltoall_int)(void *sendbuf, int sendcount, void *recvbuf, int recvcount, MPI_Comm comm, int size, int rank);
near____slint_t SL_PROTO(near____mpi_elements_keys_init_from_file)(near____elements_t *s, char *filename, near____slint from, near____slint to, near____slint const_bytes_per_line, near____slint root, int size, int rank, MPI_Comm comm);
near____slint SL_PROTO(near____mpi_elements_validate_order)(near____elements_t *s, near____slint n, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_linear_exchange_pure_keys)(near____slkey_pure_t *in, near____slkey_pure_t *out, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_check_order)(near____elements_t *s, near____slint_t nelements, near____slint_t *orders, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_check_global_order)(near____slkey_pure_t local_min, near____slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_digest_sum)(near____elements_t *s, near____slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(near____mpi_elements_crc32)(near____elements_t *s, near____slint_t n, near____slint_t keys, near____slint_t data, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_digest_hash)(near____elements_t *s, near____slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_get_counts)(near____elements_t *s, near____slint_t *clocal, near____slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
near____slweight_t SL_PROTO(near____mpi_elements_get_weights)(near____elements_t *s, near____slweight_t *wlocal, near____slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_get_counts_and_weights)(near____elements_t *s, near____slint_t nelements, near____slint_t *counts, near____slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_sendrecv_replace)(near____elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____tproc_create_tproc)(near____tproc_t *tproc, near____tproc_f *tfn, near____tproc_reset_f *rfn, near____tproc_exdef exdef);
near____slint_t SL_PROTO(near____tproc_create_tproc_mod)(near____tproc_t *tproc, near____tproc_mod_f *tfn, near____tproc_reset_f *rfn, near____tproc_exdef exdef);
near____slint_t SL_PROTO(near____tproc_create_tprocs)(near____tproc_t *tproc, near____tprocs_f *tfn, near____tproc_reset_f *rfn, near____tproc_exdef exdef);
near____slint_t SL_PROTO(near____tproc_create_tprocs_mod)(near____tproc_t *tproc, near____tprocs_mod_f *tfn, near____tproc_reset_f *rfn, near____tproc_exdef exdef);
near____slint_t SL_PROTO(near____tproc_free)(near____tproc_t *tproc);
near____slint_t SL_PROTO(near____tproc_set_proclists)(near____tproc_t *tproc, near____slint_t nsend_procs, int *send_procs, near____slint_t nrecv_procs, int *recv_procs, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____tproc_verify)(near____tproc_t tproc, void *data, near____elements_t *s, int proc);
near____slint_t SL_PROTO(near____mpi_elements_alltoall_specific)(near____elements_t *sin, near____elements_t *sout, near____elements_t *xs, near____tproc_t tproc, void *data, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_alltoallv_db_packed)(near____elements_t *sbuf, int *scounts, int *sdispls, near____elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_alltoallv_db)(near____elements_t *sbuf, int *scounts, int *sdispls, near____elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_alltoallv_ip_packed)(near____elements_t *s, near____elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_alltoallv_ip_double)(near____elements_t *s, near____elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_alltoallv_ip_mpi)(near____elements_t *s, near____elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_alltoallv_ip_dash)(near____elements_t *s, near____elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_alltoallv_ip)(near____elements_t *s, near____elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_alltoallv_proclists_db)(near____elements_t *sbuf, int *scounts, int *sdispls, int nsendprocs, int *sendprocs, near____elements_t *rbuf, int *rcounts, int *rdispls, int nrecvprocs, int *recvprocs, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, near____slint_t structured);
near____slint_t SL_PROTO(near____mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);
near____slint_t SL_PROTO(near____mpi_find_exact_equal)(near____elements_t *s, near____slint_t other_rank, near____slint_t high_rank, near____slint_t *ex_start, near____slint_t *ex_size, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_find_exact)(near____elements_t *s, near____slint_t other_rank, near____slint_t high_rank, near____slint_t *dst_size, near____slint_t *ex_start, near____slint_t *ex_sizes, near____slint_t *nx_move, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_merge2)(near____elements_t *s, near____slint_t other_rank, near____slint_t high_rank, near____slint_t *dst_size, near____merge2x_f m2, near____elements_t *xs, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_mergek_equal)(near____elements_t *s, near____sortnet_f sn, near____sortnet_data_t snd, near____merge2x_f m2x, near____elements_t *xs, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_mergek_sorted)(near____elements_t *s, near____merge2x_f m2x, near____elements_t *xs, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_mergek_sorted2)(near____elements_t *s, near____sortnet_f sn, near____sortnet_data_t snd, near____merge2x_f m2x, near____elements_t *xs, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_mergek)(near____elements_t *s, near____sortnet_f sn, near____sortnet_data_t snd, near____merge2x_f m2x, near____elements_t *xs, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_mergek_equal2)(near____elements_t *s, near____sortnet_f sn, near____sortnet_data_t snd, near____merge2x_f m2x, near____elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);
near____slint_t SL_PROTO(near____mpi_partition_exact_generic)(near____elements_t *s, near____partcond_t *pcond, near____binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_partition_exact_radix)(near____elements_t *s, near____partcond_t *pcond, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, near____slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_partition_exact_radix_ngroups)(near____elements_t *s, near____partcond_t *pcond, near____slint_t ngroups, MPI_Comm *group_comms, near____elements_t *sx, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_partition_exact_radix_2groups)(near____elements_t *s, near____partcond_t *pcond, MPI_Comm group_comm, near____elements_t *sx, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_partition_sample_regular)(near____elements_t *s, near____partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_rebalance)(near____elements_t *s0, near____elements_t *s1, near____slint_t stable, near____slint_t *dst_size, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_rebalance_alltoallv)(near____elements_t *sbuf, int *scounts, int *sdispls, near____elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
void SL_PROTO(near____mpi_partcond_set_even)(near____partcond_t *pcond, near____slint_t pcm, near____slint_t ntotal, double nimba, double wtotal, double wimba, int size, int rank);
near____slint_t SL_PROTO(near____init_partconds)(near____slint_t npconds, near____partcond_t *pconds, near____slint_t nparts, near____slint_t total_count, near____slweight_t total_weight);
near____slint_t SL_PROTO(near____init_partconds_intern)(near____slint_t npconds, near____partcond_intern_t *pci, near____partcond_t *pc, near____slint_t nparts, near____slint_t total_count, near____slweight_t total_weight);
near____slint_t SL_PROTO(near____merge_partconds)(near____partcond_t *pconds_in, near____slint_t npconds_in, near____partcond_t *pcond_out);
near____slint_t SL_PROTO(near____mpi_gather_partconds_grouped)(near____partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, near____partcond_t *pconds_out, near____slint_t *npconds_out, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_gather_partconds)(near____partcond_t *pcond_in, near____partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_allgather_partconds)(near____partcond_t *pcond_in, near____partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_bcast_partconds)(near____slint_t npconds, near____partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_post_check_partconds)(near____elements_t *s, near____slint_t nelements, near____slint_t nparts, near____partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_post_check_partconds_intern)(near____elements_t *s, near____slint_t nelements, near____slint_t nparts, near____partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_select_stats)(near____elements_t *s, near____slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_select_exact_generic_bulk)(near____elements_t *s, near____slint_t nelements, near____slint_t nparts, near____partcond_t *pconds, near____binning_t *bm, near____splitter_t *sp, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_select_exact_generic_grouped)(near____elements_t *s, near____slint_t nelements, near____partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, near____binning_t *bm, near____splitter_t *sp, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_select_exact_generic)(near____elements_t *s, near____slint_t nelements, near____slint_t nparts, near____partcond_t *pconds, near____binning_t *bm, near____splitter_t *sp, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_select_exact_radix)(near____elements_t *s, near____slint_t nelements, near____slint_t nparts, near____partcond_t *pconds, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, near____slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_select_exact_radix_grouped)(near____elements_t *s, near____slint_t nelements, near____partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, near____slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_select_sample_regular)(near____elements_t *s, near____slint_t nparts, near____partcond_t *pconds, near____slint_t nsamples, near____splitter_t *sp, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_merge)(near____elements_t *s0, near____elements_t *s1, near____elements_t *xs, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_merge2)(near____elements_t *s0, near____elements_t *s1, near____elements_t *xs, near____slint_t merge_type, near____slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_merge_radix)(near____elements_t *s0, near____elements_t *s1, near____elements_t *xs, near____slint_t merge_type, near____slint_t sort_type, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_partition)(near____elements_t *s0, near____elements_t *s1, near____elements_t *xs, near____slint_t part_type, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_partition_radix)(near____elements_t *s0, near____elements_t *s1, near____elements_t *xs, near____slint_t part_type, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_partition_exact_radix)(near____elements_t *s, near____elements_t *sx, near____partcond_t *pcond, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_partition_exact_radix_ngroups)(near____elements_t *s, near____elements_t *sx, near____partcond_t *pcond, near____slint_t ngroups, MPI_Comm *group_comms, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_partition_exact_radix_2groups)(near____elements_t *s, near____elements_t *sx, near____partcond_t *pcond, MPI_Comm group_comm, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_insert_radix)(near____elements_t *s0, near____elements_t *s1, near____elements_t *xs, near____slpkey_t *mmkeys, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_presorted_radix)(near____elements_t *s0, near____elements_t *s1, near____elements_t *xs, near____slint_t merge_type, near____slint_t rhigh, near____slint_t rlow, near____slint_t rwidth, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_sort_back)(near____elements_t *sin, near____elements_t *sout, near____elements_t *sx, near____slpkey_t *lh, near____slint_t ntotal, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, near____slint_t ytotal, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_xcounts2ycounts_grouped)(int *xcounts, near____slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
near____slint_t SL_PROTO(near____mpi_subxdispls2ycounts)(near____slint_t nsubs, int *sub_xdispls, near____slint_t *sub_sources, near____slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#undef SL_PROTO
#endif /* __SL_NEAR____H__ */
