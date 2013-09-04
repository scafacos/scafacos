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


#ifndef __SL_FORW_H__
#define __SL_FORW_H__

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif /* SL_USE_MPI */

#define SL_PROTO(_f_)  _f_

#define forw_sl_key_type_c              long long
#define forw_sl_key_type_mpi            MPI_LONG_LONG
#define forw_sl_key_size_mpi            1
#define forw_sl_key_type_fmt            "lld"

#define forw_sl_key_integer

#define forw_SL_DATA0
#define forw_sl_data0_type_c            fcs_float
#define forw_sl_data0_size_c            3
#define forw_sl_data0_type_mpi          FCS_MPI_FLOAT
#define forw_sl_data0_size_mpi          3

#define forw_SL_DATA1
#define forw_sl_data1_type_c            fcs_float
#define forw_sl_data1_size_c            1
#define forw_sl_data1_type_mpi          FCS_MPI_FLOAT
#define forw_sl_data1_size_mpi          1

#define forw_MC_ALLTOALL_INT_2STEP_THRESHOLD  1024




#if defined(MSEG_ROOT) && !defined(forw_MSEG_ROOT)
# define forw_MSEG_ROOT  MSEG_ROOT
#endif

#if defined(MSEG_BORDER_UPDATE_REDUCTION) && !defined(forw_MSEG_BORDER_UPDATE_REDUCTION)
# define forw_MSEG_BORDER_UPDATE_REDUCTION  MSEG_BORDER_UPDATE_REDUCTION
#endif

#if defined(MSEG_DISABLE_BEST_CHOICE) && !defined(forw_MSEG_DISABLE_BEST_CHOICE)
# define forw_MSEG_DISABLE_BEST_CHOICE  MSEG_DISABLE_BEST_CHOICE
#endif

#if defined(MSEG_DISABLE_MINMAX) && !defined(forw_MSEG_DISABLE_MINMAX)
# define forw_MSEG_DISABLE_MINMAX  MSEG_DISABLE_MINMAX
#endif

#if defined(MSEG_ENABLE_OPTIMZED_LOWHIGH) && !defined(forw_MSEG_ENABLE_OPTIMZED_LOWHIGH)
# define forw_MSEG_ENABLE_OPTIMZED_LOWHIGH  MSEG_ENABLE_OPTIMZED_LOWHIGH
#endif

#if defined(MSEG_FORWARD_ONLY) && !defined(forw_MSEG_FORWARD_ONLY)
# define forw_MSEG_FORWARD_ONLY  MSEG_FORWARD_ONLY
#endif

#if defined(MSEG_INFO) && !defined(forw_MSEG_INFO)
# define forw_MSEG_INFO  MSEG_INFO
#endif

#if defined(MSEG_TRACE_IF) && !defined(forw_MSEG_TRACE_IF)
# define forw_MSEG_TRACE_IF  MSEG_TRACE_IF
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


#ifndef forw_SL_INDEX
# undef forw_SL_PACKED_INDEX
#endif


/* if no special datatype for (sl default) integer ... */
#ifndef forw_sl_int_type_c
  /* ... use a default one */
# define forw_sl_int_type_c               long      /* sl_macro */
# undef forw_sl_int_type_mpi
# define forw_sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef forw_sl_int_size_mpi
# define forw_sl_int_size_mpi             1         /* sl_macro */
# undef forw_sl_int_type_fmt
# define forw_sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(forw_sl_int_type_mpi) || !defined(forw_sl_int_size_mpi)
#   error "forw_sl_int_type_mpi and/or forw_sl_int_size_mpi missing"
#  endif
# endif
# ifndef forw_sl_int_type_fmt
#  error "forw_sl_int_type_fmt macro is missing, using d as default"
#  define forw_sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for (intern) weight ... */
#ifndef forw_sl_weight_type_c
 /* ... use (sl default) integer */
# define forw_sl_weight_type_c             forw_sl_int_type_c    /* sl_macro */
# undef forw_sl_weight_type_mpi
# define forw_sl_weight_type_mpi           forw_sl_int_type_mpi  /* sl_macro */
# undef forw_sl_weight_size_mpi
# define forw_sl_weight_size_mpi           forw_sl_int_size_mpi  /* sl_macro */
# undef forw_sl_weight_type_fmt
# define forw_sl_weight_type_fmt           forw_sl_int_type_fmt  /* sl_macro */
# undef forw_sl_weight_intequiv
# define forw_sl_weight_intequiv                            /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(forw_sl_weight_type_mpi) || !defined(forw_sl_weight_size_mpi)
#   error "forw_sl_weight_type_mpi and/or forw_sl_weight_size_mpi missing"
#  endif
# endif
# ifndef forw_sl_weight_type_fmt
#  error "forw_sl_weight_type_fmt macro is missing, using f as default"
#  define forw_sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef forw_sl_index_type_c
 /* ... use the primary integer type */
# define forw_sl_index_type_c             forw_sl_int_type_c
# undef forw_sl_index_type_mpi
# define forw_sl_index_type_mpi           forw_sl_int_type_mpi
# undef forw_sl_index_size_mpi
# define forw_sl_index_size_mpi           forw_sl_int_size_mpi
# undef forw_sl_index_type_fmt
# define forw_sl_index_type_fmt           forw_sl_int_type_fmt
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(forw_sl_index_type_mpi) || !defined(forw_sl_index_size_mpi)
#   error "forw_sl_index_type_mpi and/or forw_sl_index_size_mpi missing"
#  endif
# endif
# ifndef forw_sl_index_type_fmt
#  error "forw_sl_index_type_fmt macro is missing, using d as default"
#  define forw_sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef forw_sl_key_pure_type_c
# define forw_sl_key_pure_type_c          forw_sl_key_type_c  /* sl_macro */
#endif
#ifndef forw_sl_key_pure_type_mpi
# define forw_sl_key_pure_type_mpi        forw_sl_key_type_mpi  /* sl_macro */
#endif
#ifndef forw_sl_key_pure_size_mpi
# define forw_sl_key_pure_size_mpi        forw_sl_key_size_mpi  /* sl_macro */
#endif
#ifndef forw_sl_key_pure_type_fmt
# ifdef forw_sl_key_type_fmt
#  define forw_sl_key_pure_type_fmt       forw_sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef forw_sl_key_purify
 /* key val -> key val */
 #define forw_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef forw_sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define forw_sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef forw_sl_key_set_pure
 /* key component pointer and key val */
 #define forw_sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef forw_sl_key_pure_cmp_eq
 #define forw_sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef forw_sl_key_pure_cmp_ne
 #define forw_sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef forw_sl_key_pure_cmp_lt
 #define forw_sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef forw_sl_key_pure_cmp_le
 #define forw_sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef forw_sl_key_pure_cmp_gt
 #define forw_sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef forw_sl_key_pure_cmp_ge
 #define forw_sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef forw_sl_key_cmp_eq
 #define forw_sl_key_cmp_eq(k0, k1)       (forw_sl_key_pure_cmp_eq(forw_sl_key_purify(k0), forw_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef forw_sl_key_cmp_ne
 #define forw_sl_key_cmp_ne(k0, k1)       (forw_sl_key_pure_cmp_ne(forw_sl_key_purify(k0), forw_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef forw_sl_key_cmp_lt
 #define forw_sl_key_cmp_lt(k0, k1)       (forw_sl_key_pure_cmp_lt(forw_sl_key_purify(k0), forw_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef forw_sl_key_cmp_le
 #define forw_sl_key_cmp_le(k0, k1)       (forw_sl_key_pure_cmp_le(forw_sl_key_purify(k0), forw_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef forw_sl_key_cmp_gt
 #define forw_sl_key_cmp_gt(k0, k1)       (forw_sl_key_pure_cmp_gt(forw_sl_key_purify(k0), forw_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef forw_sl_key_cmp_ge
 #define forw_sl_key_cmp_ge(k0, k1)       (forw_sl_key_pure_cmp_ge(forw_sl_key_purify(k0), forw_sl_key_purify(k1)))  /* sl_macro */
#endif


/* default random key */
#ifdef forw_sl_key_integer
# if !defined(forw_sl_key_val_srand) || !defined(forw_sl_key_val_rand) || !defined(forw_sl_key_val_rand_minmax)
#  undef forw_sl_key_val_srand
#  undef forw_sl_key_val_rand
#  undef forw_sl_key_val_rand_minmax
#  define forw_sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define forw_sl_key_val_rand()                     ((forw_sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define forw_sl_key_val_rand_minmax(_min_, _max_)  ((forw_sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef forw_SL_DATA0_IGNORE
# undef forw_SL_DATA0
#endif
#ifdef forw_SL_DATA1_IGNORE
# undef forw_SL_DATA1
#endif
#ifdef forw_SL_DATA2_IGNORE
# undef forw_SL_DATA2
#endif
#ifdef forw_SL_DATA3_IGNORE
# undef forw_SL_DATA3
#endif
#ifdef forw_SL_DATA4_IGNORE
# undef forw_SL_DATA4
#endif
#ifdef forw_SL_DATA5_IGNORE
# undef forw_SL_DATA5
#endif
#ifdef forw_SL_DATA6_IGNORE
# undef forw_SL_DATA6
#endif
#ifdef forw_SL_DATA7_IGNORE
# undef forw_SL_DATA7
#endif
#ifdef forw_SL_DATA8_IGNORE
# undef forw_SL_DATA8
#endif
#ifdef forw_SL_DATA9_IGNORE
# undef forw_SL_DATA9
#endif
#ifdef forw_SL_DATA10_IGNORE
# undef forw_SL_DATA10
#endif
#ifdef forw_SL_DATA11_IGNORE
# undef forw_SL_DATA11
#endif
#ifdef forw_SL_DATA12_IGNORE
# undef forw_SL_DATA12
#endif
#ifdef forw_SL_DATA13_IGNORE
# undef forw_SL_DATA13
#endif
#ifdef forw_SL_DATA14_IGNORE
# undef forw_SL_DATA14
#endif
#ifdef forw_SL_DATA15_IGNORE
# undef forw_SL_DATA15
#endif
#ifdef forw_SL_DATA16_IGNORE
# undef forw_SL_DATA16
#endif
#ifdef forw_SL_DATA17_IGNORE
# undef forw_SL_DATA17
#endif
#ifdef forw_SL_DATA18_IGNORE
# undef forw_SL_DATA18
#endif
#ifdef forw_SL_DATA19_IGNORE
# undef forw_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro forw_sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef forw_sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef forw_sl_data0_weight
# undef forw_sl_data1_weight
# undef forw_sl_data2_weight
# undef forw_sl_data3_weight
# undef forw_sl_data4_weight
# undef forw_sl_data5_weight
# undef forw_sl_data6_weight
# undef forw_sl_data7_weight
# undef forw_sl_data8_weight
# undef forw_sl_data9_weight
# undef forw_sl_data10_weight
# undef forw_sl_data11_weight
# undef forw_sl_data12_weight
# undef forw_sl_data13_weight
# undef forw_sl_data14_weight
# undef forw_sl_data15_weight
# undef forw_sl_data16_weight
# undef forw_sl_data17_weight
# undef forw_sl_data18_weight
# undef forw_sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable forw_sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(forw_sl_data0_weight) && !defined(forw_SL_DATA0)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data1_weight) && !defined(forw_SL_DATA1)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data2_weight) && !defined(forw_SL_DATA2)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data3_weight) && !defined(forw_SL_DATA3)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data4_weight) && !defined(forw_SL_DATA4)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data5_weight) && !defined(forw_SL_DATA5)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data6_weight) && !defined(forw_SL_DATA6)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data7_weight) && !defined(forw_SL_DATA7)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data8_weight) && !defined(forw_SL_DATA8)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data9_weight) && !defined(forw_SL_DATA9)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data10_weight) && !defined(forw_SL_DATA10)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data11_weight) && !defined(forw_SL_DATA11)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data12_weight) && !defined(forw_SL_DATA12)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data13_weight) && !defined(forw_SL_DATA13)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data14_weight) && !defined(forw_SL_DATA14)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data15_weight) && !defined(forw_SL_DATA15)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data16_weight) && !defined(forw_SL_DATA16)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data17_weight) && !defined(forw_SL_DATA17)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data18_weight) && !defined(forw_SL_DATA18)
# undef forw_sl_elem_weight
#endif
#if defined(forw_sl_data19_weight) && !defined(forw_SL_DATA19)
# undef forw_sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro forw_FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef forw_SL_DATA0
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data0_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA1
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data1_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA2
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data2_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA3
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data3_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA4
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data4_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA5
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data5_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA6
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data6_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA7
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data7_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA8
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data8_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA9
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data9_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA10
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data10_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA11
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data11_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA12
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data12_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA13
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data13_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA14
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data14_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA15
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data15_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA16
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data16_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA17
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data17_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA18
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data18_flex
#   define forw_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef forw_SL_DATA19
# ifdef forw_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef forw_sl_data19_flex
#   define forw_FLECKS_GUARD
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






#define forw_SPEC_TLOC

typedef forw_sl_int_type_c forw_spec_int_t;

typedef int forw_spec_proc_t;

#define forw_SPEC_LOC_NONE   -1
#ifdef SL_USE_MPI
# define forw_SPEC_PROC_NONE  MPI_PROC_NULL
#else
# define forw_SPEC_PROC_NONE  -1
#endif

typedef void *forw_spec_tloc_data_t;
typedef void *forw_spec_tproc_data_t;

struct forw__elements_t;

typedef struct forw__elements_t *forw_spec_elem_buf_t;

typedef struct forw__elements_t forw_spec_elem_t;

typedef forw_sl_int_type_c forw_spec_elem_index_t;

#define forw_spec_elem_set_n(_e_, _n_)     forw_elem_set_size((_e_), (_n_))
#define forw_spec_elem_get_n(_e_)          forw_elem_get_size((_e_))
#define forw_spec_elem_set_nmax(_e_, _n_)  forw_elem_set_max_size((_e_), (_n_))
#define forw_spec_elem_get_nmax(_e_)       forw_elem_get_max_size((_e_))

#define forw_spec_elem_set_buf(_e_, _b_)   *(_e_) = *(_b_)
#define forw_spec_elem_get_buf(_e_)        (_e_)

#define forw_spec_elem_copy_at(_se_, _sat_, _de_, _dat_) \
  elem_copy_at((_se_), (_sat_), (_de_), (_dat_))

#define forw_spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) \
  elem_xchange_at((_s0_), (_s0at_), (_s1_), (_s1at_), (_t_))






/* tproc count */

/* sp_macro forw_SPEC_DECLARE_TPROC_COUNT_DB */
#define forw_SPEC_DECLARE_TPROC_COUNT_DB \
  struct { forw_spec_elem_index_t i; forw_spec_proc_t p; } spec0cd;

/* sp_macro forw_SPEC_DO_TPROC_COUNT_DB */
#define forw_SPEC_DO_TPROC_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec0cd.i = 0; spec0cd.i < forw_spec_elem_get_n(_b_); ++spec0cd.i) { \
    spec0cd.p = (_tp_)(forw_spec_elem_get_buf(_b_), spec0cd.i, _tpd_); \
    if (spec0cd.p == forw_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0cd.p]; \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROC_COUNT_DB */
#define forw_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_db(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts) \
{ \
  forw_SPEC_DECLARE_TPROC_COUNT_DB \
  forw_SPEC_DO_TPROC_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro forw_SPEC_DECLARE_TPROC_COUNT_IP */
#define forw_SPEC_DECLARE_TPROC_COUNT_IP \
  struct { forw_spec_elem_index_t i, t; forw_spec_proc_t p; } spec0ci;

/* sp_macro forw_SPEC_DO_TPROC_COUNT_IP */
#define forw_SPEC_DO_TPROC_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec0ci.t = 0; \
  for (spec0ci.i = 0; spec0ci.i < forw_spec_elem_get_n(_b_); ++spec0ci.i) { \
    spec0ci.p = (_tp_)(forw_spec_elem_get_buf(_b_), spec0ci.i, _tpd_); \
    if (spec0ci.p == forw_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec0ci.p]; \
    if (spec0ci.t < spec0ci.i) forw_spec_elem_copy_at((_b_), spec0ci.i, (_b_), spec0ci.t); \
    ++spec0ci.t; \
  } \
  forw_spec_elem_set_n(_b_, spec0ci.t); \
} while (0)

/* sp_macro forw_SPEC_FUNC_TPROC_COUNT_IP */
#define forw_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_count_ip(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts) \
{ \
  forw_SPEC_DECLARE_TPROC_COUNT_IP \
  forw_SPEC_DO_TPROC_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tproc_mod count */

/* sp_macro forw_SPEC_DECLARE_TPROC_MOD_COUNT_DB */
#define forw_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  struct { forw_spec_elem_index_t i; forw_spec_proc_t p; } spec1cd;

/* sp_macro forw_SPEC_DO_TPROC_MOD_COUNT_DB */
#define forw_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_)  do { \
  for (spec1cd.i = 0; spec1cd.i < forw_spec_elem_get_n(_b_); ++spec1cd.i) { \
    spec1cd.p = (_tp_)(forw_spec_elem_get_buf(_b_), spec1cd.i, _tpd_, NULL); \
    if (spec1cd.p == forw_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1cd.p]; \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROC_MOD_COUNT_DB */
#define forw_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_db(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts) \
{ \
  forw_SPEC_DECLARE_TPROC_MOD_COUNT_DB \
  forw_SPEC_DO_TPROC_MOD_COUNT_DB(_tp_, tproc_data, s, counts); \
}

/* sp_macro forw_SPEC_DECLARE_TPROC_MOD_COUNT_IP */
#define forw_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  struct { forw_spec_elem_index_t i, t; forw_spec_proc_t p; } spec1ci;

/* sp_macro forw_SPEC_DO_TPROC_MOD_COUNT_IP */
#define forw_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_)  do { \
  spec1ci.t = 0; \
  for (spec1ci.i = 0; spec1ci.i < forw_spec_elem_get_n(_b_); ++spec1ci.i) { \
    spec1ci.p = (_tp_)(forw_spec_elem_get_buf(_b_), spec1ci.i, _tpd_, NULL); \
    if (spec1ci.p == forw_SPEC_PROC_NONE) continue; \
    ++(_cs_)[spec1ci.p]; \
    if (spec1ci.t < spec1ci.i) forw_spec_elem_copy_at((_b_), spec1ci.i, (_b_), spec1ci.t); \
    ++spec1ci.t; \
  } \
  forw_spec_elem_set_n(_b_, spec1ci.t); \
} while (0)

/* sp_macro forw_SPEC_FUNC_TPROC_MOD_COUNT_IP */
#define forw_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_count_ip(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts) \
{ \
  forw_SPEC_DECLARE_TPROC_MOD_COUNT_IP \
  forw_SPEC_DO_TPROC_MOD_COUNT_IP(_tp_, tproc_data, s, counts); \
}


/* tprocs count */

/* sp_macro forw_SPEC_DECLARE_TPROCS_COUNT_DB */
#define forw_SPEC_DECLARE_TPROCS_COUNT_DB \
  struct { forw_spec_elem_index_t i; forw_spec_int_t j, n; } spec2cd;

/* sp_macro forw_SPEC_DO_TPROCS_COUNT_DB */
#define forw_SPEC_DO_TPROCS_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec2cd.i = 0; spec2cd.i < forw_spec_elem_get_n(_b_); ++spec2cd.i) { \
    spec2cd.n = (_tp_)(forw_spec_elem_get_buf(_b_), spec2cd.i, (_tpd_), (_ps_)); \
    for (spec2cd.j = 0; spec2cd.j < spec2cd.n; ++spec2cd.j) ++(_cs_)[(_ps_)[spec2cd.j]]; \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROCS_COUNT_DB */
#define forw_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_db(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts, forw_spec_proc_t *procs) \
{ \
  forw_SPEC_DECLARE_TPROCS_COUNT_DB \
  forw_SPEC_DO_TPROCS_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro forw_SPEC_DECLARE_TPROCS_COUNT_IP */
#define forw_SPEC_DECLARE_TPROCS_COUNT_IP \
  struct { forw_spec_elem_index_t i, t; forw_spec_int_t j, n; } spec2ci;

/* sp_macro forw_SPEC_DO_TPROCS_COUNT_IP */
#define forw_SPEC_DO_TPROCS_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec2ci.t = 0; \
  for (spec2ci.i = 0; spec2ci.i < forw_spec_elem_get_n(_b_); ++spec2ci.i) { \
    spec2ci.n = (_tp_)(forw_spec_elem_get_buf(_b_), spec2ci.i, (_tpd_), (_ps_)); \
    if (spec2ci.n <= 0) continue; \
    for (spec2ci.j = 0; spec2ci.j < spec2ci.n; ++spec2ci.j) ++(_cs_)[(_ps_)[spec2ci.j]]; \
    if (spec2ci.t < spec2ci.i) forw_spec_elem_copy_at((_b_), spec2ci.i, (_b_), spec2ci.t); \
    ++spec2ci.t; \
  } \
  forw_spec_elem_set_n(_b_, spec2ci.t); \
} while (0)

/* sp_macro forw_SPEC_FUNC_TPROCS_COUNT_IP */
#define forw_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_count_ip(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts, forw_spec_proc_t *procs) \
{ \
  forw_SPEC_DECLARE_TPROCS_COUNT_IP \
  forw_SPEC_DO_TPROCS_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tprocs_mod count */

/* sp_macro forw_SPEC_DECLARE_TPROCS_MOD_COUNT_DB */
#define forw_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  struct { forw_spec_elem_index_t i; forw_spec_int_t j, n; } spec3cd;

/* sp_macro forw_SPEC_DO_TPROCS_MOD_COUNT_DB */
#define forw_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  for (spec3cd.i = 0; spec3cd.i < forw_spec_elem_get_n(_b_); ++spec3cd.i) \
  { \
    spec3cd.n = (_tp_)(forw_spec_elem_get_buf(_b_), spec3cd.i, (_tpd_), (_ps_), NULL); \
    for (spec3cd.j = 0; spec3cd.j < spec3cd.n; ++spec3cd.j) ++(_cs_)[(_ps_)[spec3cd.j]]; \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROCS_MOD_COUNT_DB */
#define forw_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_db(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts, forw_spec_proc_t *procs) \
{ \
  forw_SPEC_DECLARE_TPROCS_MOD_COUNT_DB \
  forw_SPEC_DO_TPROCS_MOD_COUNT_DB(_tp_, tproc_data, s, counts, procs); \
}

/* sp_macro forw_SPEC_DECLARE_TPROCS_MOD_COUNT_IP */
#define forw_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  struct { forw_spec_elem_index_t i, t; forw_spec_int_t j, n; } spec3ci;

/* sp_macro forw_SPEC_DO_TPROCS_MOD_COUNT_IP */
#define forw_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, _tpd_, _b_, _cs_, _ps_)  do { \
  spec3ci.t = 0; \
  for (spec3ci.i = 0; spec3ci.i < forw_spec_elem_get_n(_b_); ++spec3ci.i) { \
    spec3ci.n = (_tp_)(forw_spec_elem_get_buf(_b_), spec3ci.i, (_tpd_), (_ps_), NULL); \
    if (spec3ci.n <= 0) continue; \
    for (spec3ci.j = 0; spec3ci.j < spec3ci.n; ++spec3ci.j) ++(_cs_)[(_ps_)[spec3ci.j]]; \
    if (spec3ci.t < spec3ci.i) forw_spec_elem_copy_at((_b_), spec3ci.i, (_b_), spec3ci.t); \
    ++spec3ci.t; \
  } \
  forw_spec_elem_set_n(_b_, spec3ci.t); \
} while (0)

/* sp_macro forw_SPEC_FUNC_TPROCS_MOD_COUNT_IP */
#define forw_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_count_ip(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts, forw_spec_proc_t *procs) \
{ \
  forw_SPEC_DECLARE_TPROCS_MOD_COUNT_IP \
  forw_SPEC_DO_TPROCS_MOD_COUNT_IP(_tp_, tproc_data, s, counts, procs); \
}


/* tproc rearrange */

/* sp_macro forw_SPEC_DECLARE_TPROC_REARRANGE_DB */
#define forw_SPEC_DECLARE_TPROC_REARRANGE_DB \
  struct { forw_spec_elem_index_t i; forw_spec_proc_t p; } spec0d;

/* sp_macro forw_SPEC_DO_TPROC_REARRANGE_DB */
#define forw_SPEC_DO_TPROC_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_)  do { \
  for (spec0d.i = 0; spec0d.i < forw_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tp_)(forw_spec_elem_get_buf(_sb_), spec0d.i, _tpd_); \
    if (spec0d.p == forw_SPEC_PROC_NONE) continue; \
    forw_spec_elem_copy_at((_sb_), spec0d.i, (_db_), (_ds_)[spec0d.p]); \
    ++(_ds_)[spec0d.p]; \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROC_REARRANGE_DB */
#define forw_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_rearrange_db(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tproc_data_t tproc_data, int *displs) \
{ \
  forw_SPEC_DECLARE_TPROC_REARRANGE_DB \
  forw_SPEC_DO_TPROC_REARRANGE_DB(_tp_, tproc_data, s, d, displs); \
}

/* sp_macro forw_SPEC_DECLARE_TPROC_REARRANGE_IP */
#define forw_SPEC_DECLARE_TPROC_REARRANGE_IP \
  struct { forw_spec_elem_index_t e, i, j; forw_spec_proc_t p, np; } spec0i;

/* sp_macro forw_SPEC_DO_TPROC_REARRANGE_IP */
#define forw_SPEC_DO_TPROC_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_)  do { \
  for (spec0i.e = 0, spec0i.i = 0; spec0i.i < (_n_); ++spec0i.i) { \
    spec0i.e += (_cs_)[spec0i.i]; \
    spec0i.j = (_ds_)[spec0i.i]; \
    while (spec0i.j < spec0i.e) { \
      spec0i.p = (_tp_)(forw_spec_elem_get_buf(_b_), spec0i.j, _tpd_); \
      while (spec0i.p != spec0i.i) { \
        spec0i.np = (_tp_)(forw_spec_elem_get_buf(_b_), (_ds_)[spec0i.p], _tpd_); \
        if (spec0i.np != spec0i.p) forw_spec_elem_exchange_at((_b_), (_ds_)[spec0i.p], (_b_), spec0i.j, (_xb_)); \
        ++(_ds_)[spec0i.p]; \
        spec0i.p = spec0i.np; \
      } \
      ++spec0i.j; \
    } \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROC_REARRANGE_IP */
#define forw_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_rearrange_ip(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tproc_data_t tproc_data, int *displs, int *counts, forw_spec_int_t n) \
{ \
  forw_SPEC_DECLARE_TPROC_REARRANGE_IP \
  forw_SPEC_DO_TPROC_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n); \
}


/* tproc_mod rearrange */

/* sp_macro forw_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB */
#define forw_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  struct { forw_spec_elem_index_t i; forw_spec_proc_t p; } spec1d;

/* sp_macro forw_SPEC_DO_TPROC_MOD_REARRANGE_DB */
#define forw_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < forw_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(forw_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, forw_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == forw_SPEC_PROC_NONE) continue; \
      forw_spec_elem_copy_at((_ib_), 0, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < forw_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tp_)(forw_spec_elem_get_buf(_sb_), spec1d.i, _tpd_, NULL); \
      if (spec1d.p == forw_SPEC_PROC_NONE) continue; \
      forw_spec_elem_copy_at((_sb_), spec1d.i, (_db_), (_ds_)[spec1d.p]); \
      ++(_ds_)[spec1d.p]; \
    } \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROC_MOD_REARRANGE_DB */
#define forw_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tproc_mod_rearrange_db(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tproc_data_t tproc_data, int *displs, forw_spec_elem_t *mod) \
{ \
  forw_SPEC_DECLARE_TPROC_MOD_REARRANGE_DB \
  forw_SPEC_DO_TPROC_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, mod); \
}

/* sp_macro forw_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP */
#define forw_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  struct { forw_spec_elem_index_t e, i, j; forw_spec_proc_t p, np; } spec1i;

/* sp_macro forw_SPEC_DO_TPROC_MOD_REARRANGE_IP */
#define forw_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.e = 0, spec1i.i = 0; spec1i.i < (_n_); ++spec1i.i) { \
      spec1i.e += (_cs_)[spec1i.i]; \
      spec1i.j = (_ds_)[spec1i.i]; \
      while (spec1i.j < spec1i.e) { \
        spec1i.p = (_tp_)(forw_spec_elem_get_buf(_b_), spec1i.j, _tpd_, forw_spec_elem_get_buf(_ib_)); \
        forw_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(forw_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, forw_spec_elem_get_buf(_ib_)); \
          if (spec1i.np != spec1i.p) { \
            forw_spec_elem_copy_at((_b_), spec1i.j, (_b_), (_ds_)[spec1i.p]); \
            forw_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.j); \
          } else forw_spec_elem_copy_at((_ib_), 0, (_b_), (_ds_)[spec1i.p]); \
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
        spec1i.p = (_tp_)(forw_spec_elem_get_buf(_b_), spec1i.j, _tpd_, NULL); \
        while (spec1i.p != spec1i.i) { \
          spec1i.np = (_tp_)(forw_spec_elem_get_buf(_b_), (_ds_)[spec1i.p], _tpd_, NULL); \
          if (spec1i.np != spec1i.p) forw_spec_elem_exchange_at((_b_), (_ds_)[spec1i.p], (_b_), spec1i.j, (_xb_)); \
          ++(_ds_)[spec1i.p]; \
          spec1i.p = spec1i.np; \
        } \
        ++spec1i.j; \
      } \
    } \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROC_MOD_REARRANGE_IP */
#define forw_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_) \
_s_ void _name_##_tproc_mod_rearrange_ip(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tproc_data_t tproc_data, int *displs, int *counts, forw_spec_int_t n, forw_spec_elem_t *mod) \
{ \
  forw_SPEC_DECLARE_TPROC_MOD_REARRANGE_IP \
  forw_SPEC_DO_TPROC_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, mod); \
}


/* tprocs rearrange */

/* sp_macro forw_SPEC_DECLARE_TPROCS_REARRANGE_DB */
#define forw_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  struct { forw_spec_elem_index_t i; forw_spec_int_t j, n; } spec2d;

/* sp_macro forw_SPEC_DO_TPROCS_REARRANGE_DB */
#define forw_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_)  do { \
  for (spec2d.i = 0; spec2d.i < forw_spec_elem_get_n(_sb_); ++spec2d.i) { \
    spec2d.n = (_tp_)(forw_spec_elem_get_buf(_sb_), spec2d.i, (_tpd_), (_ps_)); \
    for (spec2d.j = 0; spec2d.j < spec2d.n; ++spec2d.j) { \
      forw_spec_elem_copy_at((_sb_), spec2d.i, (_db_), (_ds_)[(_ps_)[spec2d.j]]); \
      ++(_ds_)[(_ps_)[spec2d.j]]; \
    } \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROCS_REARRANGE_DB */
#define forw_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_db(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tproc_data_t tproc_data, int *displs, forw_spec_proc_t *procs) \
{ \
  forw_SPEC_DECLARE_TPROCS_REARRANGE_DB \
  forw_SPEC_DO_TPROCS_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs); \
}

/* sp_macro forw_SPEC_DECLARE_TPROCS_REARRANGE_IP */
#define forw_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  struct { forw_spec_elem_index_t e, j, fe, fc, le, lc; forw_spec_int_t i, n, f, l, o; } spec2i;

/* sp_macro forw_SPEC_DO_TPROCS_REARRANGE_IP */
#define forw_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_)  do { \
  spec2i.f = 0; spec2i.fe = (_cs_)[0]; spec2i.fc = forw_spec_elem_get_n(_b_); \
  while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; } \
  spec2i.l = 0; spec2i.le = (_cs_)[0]; spec2i.lc = forw_spec_elem_get_n(_b_) - 1; \
  while (spec2i.lc >= spec2i.le) { ++spec2i.l; spec2i.le += (_cs_)[spec2i.l]; } \
  for (spec2i.e = 0, spec2i.i = 0; spec2i.i < (_n_); ++spec2i.i) { \
    spec2i.e += (_cs_)[spec2i.i]; \
    spec2i.j = (_ds_)[spec2i.i]; \
    while (spec2i.j < spec2i.e) { \
      spec2i.n = (_tp_)(forw_spec_elem_get_buf(_b_), spec2i.j, (_tpd_), (_ps_)); \
      spec2i.o = -1; \
      while (spec2i.n > 0) { \
        --spec2i.n; \
        if ((_ps_)[spec2i.n] == spec2i.i && spec2i.o < 0) spec2i.o = spec2i.n; \
        else if ((_ds_)[(_ps_)[spec2i.n]] < spec2i.fc) { \
          spec2i.l = spec2i.f; spec2i.le = spec2i.fe; spec2i.lc = spec2i.fc; \
          if (spec2i.fc < spec2i.fe) { \
            forw_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_b_), spec2i.fc); \
            ++spec2i.fc; \
          } else forw_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec2i.n]], (_xb_), 0); \
        } else if ((_ds_)[(_ps_)[spec2i.n]] == spec2i.fc) ++spec2i.fc; \
        if (spec2i.j != (_ds_)[(_ps_)[spec2i.n]]) forw_spec_elem_copy_at((_b_), spec2i.j, (_b_), (_ds_)[(_ps_)[spec2i.n]]); \
        ++(_ds_)[(_ps_)[spec2i.n]]; \
        while (spec2i.f + 1 < (_n_) && spec2i.fc >= spec2i.fe) { ++spec2i.f; spec2i.fe += (_cs_)[spec2i.f]; spec2i.fc = (_ds_)[spec2i.f]; } \
      } \
      if (spec2i.o < 0) { \
        if (spec2i.lc < spec2i.le) {  \
          forw_spec_elem_copy_at((_b_), spec2i.lc, (_b_), spec2i.j); \
          spec2i.f = spec2i.l; spec2i.fe = spec2i.le; spec2i.fc = spec2i.lc; \
          --spec2i.lc; \
          while (spec2i.l > 0 && spec2i.lc < (_ds_)[spec2i.l]) { spec2i.le -= (_cs_)[spec2i.l]; spec2i.lc = spec2i.le - 1; --spec2i.l; } \
        } else forw_spec_elem_copy_at((_xb_), 0, (_b_), spec2i.j); \
      } \
      spec2i.j = (_ds_)[spec2i.i]; \
    } \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROCS_REARRANGE_IP */
#define forw_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_rearrange_ip(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tproc_data_t tproc_data, int *displs, int *counts, forw_spec_int_t n, forw_spec_proc_t *procs) \
{ \
  forw_SPEC_DECLARE_TPROCS_REARRANGE_IP \
  forw_SPEC_DO_TPROCS_REARRANGE_IP(_tp_, tproc_data, s, d, displs, counts, n, procs); \
}


/* tprocs_mod rearrange */

/* sp_macro forw_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB */
#define forw_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  struct { forw_spec_elem_index_t i; forw_spec_int_t j, n; } spec3d;

/* sp_macro forw_SPEC_DO_TPROCS_MOD_REARRANGE_DB */
#define forw_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, _tpd_, _sb_, _db_, _ds_, _ps_, _ib_)  do { \
  if (_ib_) { \
    for (spec3d.i = 0; spec3d.i < forw_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(forw_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), forw_spec_elem_get_buf(_ib_)); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        forw_spec_elem_copy_at((_ib_), spec3d.j, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } else { \
    for (spec3d.i = 0; spec3d.i < forw_spec_elem_get_n(_sb_); ++spec3d.i) { \
      spec3d.n = (_tp_)(forw_spec_elem_get_buf(_sb_), spec3d.i, (_tpd_), (_ps_), NULL); \
      for (spec3d.j = 0; spec3d.j < spec3d.n; ++spec3d.j) { \
        forw_spec_elem_copy_at((_sb_), spec3d.i, (_db_), (_ds_)[(_ps_)[spec3d.j]]); \
        ++(_ds_)[(_ps_)[spec3d.j]]; \
      } \
    } \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB */
#define forw_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_db(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tproc_data_t tproc_data, int *displs, forw_spec_proc_t *procs, forw_spec_elem_t *mod) \
{ \
  forw_SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB \
  forw_SPEC_DO_TPROCS_MOD_REARRANGE_DB(_tp_, tproc_data, s, d, displs, procs, mod); \
}

/* sp_macro forw_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP */
#define forw_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  struct { forw_spec_elem_index_t e, j, fe, fc, le, lc; forw_spec_int_t i, n, f, l, o; } spec3i;

/* sp_macro forw_SPEC_DO_TPROCS_MOD_REARRANGE_IP */
#define forw_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, _tpd_, _b_, _xb_, _ds_, _cs_, _n_, _ps_, _ib_)  do { \
  if (_ib_) { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = forw_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = forw_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(forw_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), forw_spec_elem_get_buf(_ib_)); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              forw_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else forw_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          forw_spec_elem_copy_at((_ib_), spec3i.n, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            forw_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else forw_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } else { \
    spec3i.f = 0; spec3i.fe = (_cs_)[0]; spec3i.fc = forw_spec_elem_get_n(_b_); \
    while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; } \
    spec3i.l = 0; spec3i.le = (_cs_)[0]; spec3i.lc = forw_spec_elem_get_n(_b_) - 1; \
    while (spec3i.lc >= spec3i.le) { ++spec3i.l; spec3i.le += (_cs_)[spec3i.l]; } \
    for (spec3i.e = 0, spec3i.i = 0; spec3i.i < (_n_); ++spec3i.i) { \
      spec3i.e += (_cs_)[spec3i.i]; \
      spec3i.j = (_ds_)[spec3i.i]; \
      while (spec3i.j < spec3i.e) { \
        spec3i.n = (_tp_)(forw_spec_elem_get_buf(_b_), spec3i.j, (_tpd_), (_ps_), NULL); \
        spec3i.o = -1; \
        while (spec3i.n > 0) { \
          --spec3i.n; \
          if ((_ps_)[spec3i.n] == spec3i.i && spec3i.o < 0) spec3i.o = spec3i.n; \
          else if ((_ds_)[(_ps_)[spec3i.n]] < spec3i.fc) { \
            spec3i.l = spec3i.f; spec3i.le = spec3i.fe; spec3i.lc = spec3i.fc; \
            if (spec3i.fc < spec3i.fe) { \
              forw_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_b_), spec3i.fc); \
              ++spec3i.fc; \
            } else forw_spec_elem_copy_at((_b_), (_ds_)[(_ps_)[spec3i.n]], (_xb_), 0); \
          } else if ((_ds_)[(_ps_)[spec3i.n]] == spec3i.fc) ++spec3i.fc; \
          if (spec3i.j != (_ds_)[(_ps_)[spec3i.n]]) forw_spec_elem_copy_at((_b_), spec3i.j, (_b_), (_ds_)[(_ps_)[spec3i.n]]); \
          ++(_ds_)[(_ps_)[spec3i.n]]; \
          while (spec3i.f + 1 < (_n_) && spec3i.fc >= spec3i.fe) { ++spec3i.f; spec3i.fe += (_cs_)[spec3i.f]; spec3i.fc = (_ds_)[spec3i.f]; } \
        } \
        if (spec3i.o < 0) { \
          if (spec3i.lc < spec3i.le) {  \
            forw_spec_elem_copy_at((_b_), spec3i.lc, (_b_), spec3i.j); \
            spec3i.f = spec3i.l; spec3i.fe = spec3i.le; spec3i.fc = spec3i.lc; \
            --spec3i.lc; \
            while (spec3i.l > 0 && spec3i.lc < (_ds_)[spec3i.l]) { spec3i.le -= (_cs_)[spec3i.l]; spec3i.lc = spec3i.le - 1; --spec3i.l; } \
          } else forw_spec_elem_copy_at((_xb_), 0, (_b_), spec3i.j); \
        } \
        spec3i.j = (_ds_)[spec3i.i]; \
      } \
    } \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP */
#define forw_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_...) \
_s_ void _name_##_tprocs_mod_rearrange_ip(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tproc_data_t tproc_data, int *displs, int *counts, forw_spec_int_t n, forw_spec_proc_t *procs, forw_spec_elem_t *mod) \
{ \
  forw_SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP \
  forw_SPEC_DO_TPROCS_MOD_REARRANGE_IP(_tp_, tproc_data, s, x, displs, counts, n, procs, mod); \
}

/* sp_macro forw_SPEC_DEFINE_TPROC */
#define forw_SPEC_DEFINE_TPROC(_name_, _tp_, _s_...) \
  forw_SPEC_FUNC_TPROC_COUNT_DB(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROC_COUNT_IP(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROC_REARRANGE_DB(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROC_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro forw_SPEC_DEFINE_TPROC_MOD */
#define forw_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  forw_SPEC_FUNC_TPROC_MOD_COUNT_DB(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROC_MOD_COUNT_IP(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROC_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROC_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro forw_SPEC_DEFINE_TPROCS */
#define forw_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  forw_SPEC_FUNC_TPROCS_COUNT_DB(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROCS_COUNT_IP(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROCS_REARRANGE_DB(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROCS_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro forw_SPEC_DEFINE_TPROCS_MOD */
#define forw_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  forw_SPEC_FUNC_TPROCS_MOD_COUNT_DB(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROCS_MOD_COUNT_IP(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROCS_MOD_REARRANGE_DB(_name_, _tp_, _s_) \
  forw_SPEC_FUNC_TPROCS_MOD_REARRANGE_IP(_name_, _tp_, _s_)

/* sp_macro forw_SPEC_EXT_PARAM_TPROC forw_SPEC_EXT_PARAM_TPROC_NULL forw_SPEC_EXT_PARAM_TPROC_MOD forw_SPEC_EXT_PARAM_TPROC_MOD_NULL forw_SPEC_EXT_PARAM_TPROCS forw_SPEC_EXT_PARAM_TPROCS_NULL forw_SPEC_EXT_PARAM_TPROCS_MOD forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL */
#define forw_SPEC_EXT_PARAM_TPROC(_name_)       _name_##_tproc_count_db, _name_##_tproc_count_ip, _name_##_tproc_rearrange_db, _name_##_tproc_rearrange_ip
#define forw_SPEC_EXT_PARAM_TPROC_NULL          NULL, NULL, NULL, NULL
#define forw_SPEC_EXT_PARAM_TPROC_MOD(_name_)   _name_##_tproc_mod_count_db, _name_##_tproc_mod_count_ip, _name_##_tproc_mod_rearrange_db, _name_##_tproc_mod_rearrange_ip
#define forw_SPEC_EXT_PARAM_TPROC_MOD_NULL      NULL, NULL, NULL, NULL
#define forw_SPEC_EXT_PARAM_TPROCS(_name_)      _name_##_tprocs_count_db, _name_##_tprocs_count_ip, _name_##_tprocs_rearrange_db, _name_##_tprocs_rearrange_ip
#define forw_SPEC_EXT_PARAM_TPROCS_NULL         NULL, NULL, NULL, NULL
#define forw_SPEC_EXT_PARAM_TPROCS_MOD(_name_)  _name_##_tprocs_mod_count_db, _name_##_tprocs_mod_count_ip, _name_##_tprocs_mod_rearrange_db, _name_##_tprocs_mod_rearrange_ip
#define forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL     NULL, NULL, NULL, NULL


/* sp_type forw_spec_tproc_f forw_spec_tproc_count_f forw_spec_tproc_rearrange_db_f forw_spec_tproc_rearrange_ip_f */
typedef forw_spec_proc_t forw_spec_tproc_f(forw_spec_elem_buf_t b, forw_spec_elem_index_t x, forw_spec_tproc_data_t tproc_data);
typedef void forw_spec_tproc_count_f(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts);
typedef void forw_spec_tproc_rearrange_db_f(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tproc_data_t tproc_data, int *displs);
typedef void forw_spec_tproc_rearrange_ip_f(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tproc_data_t tproc_data, int *displs, int *counts, forw_spec_int_t n);

/* sp_type forw_spec_tproc_mod_f forw_spec_tproc_mod_count_f forw_spec_tproc_mod_rearrange_db_f forw_spec_tproc_mod_rearrange_ip_f */
typedef forw_spec_proc_t forw_spec_tproc_mod_f(forw_spec_elem_buf_t b, forw_spec_elem_index_t x, forw_spec_tproc_data_t tproc_data, forw_spec_elem_buf_t mod);
typedef void forw_spec_tproc_mod_count_f(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts);
typedef void forw_spec_tproc_mod_rearrange_db_f(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tproc_data_t tproc_data, int *displs, forw_spec_elem_t *mod);
typedef void forw_spec_tproc_mod_rearrange_ip_f(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tproc_data_t tproc_data, int *displs, int *counts, forw_spec_int_t n, forw_spec_elem_t *mod);

/* sp_type forw_spec_tprocs_f forw_spec_tprocs_count_f forw_spec_tprocs_rearrange_db_f forw_spec_tprocs_rearrange_ip_f */
typedef forw_spec_int_t forw_spec_tprocs_f(forw_spec_elem_buf_t b, forw_spec_elem_index_t x, forw_spec_tproc_data_t tproc_data, forw_spec_proc_t *procs);
typedef void forw_spec_tprocs_count_f(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts, forw_spec_proc_t *procs);
typedef void forw_spec_tprocs_rearrange_db_f(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tproc_data_t tproc_data, int *displs, forw_spec_proc_t *procs);
typedef void forw_spec_tprocs_rearrange_ip_f(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tproc_data_t tproc_data, int *displs, int *counts, forw_spec_int_t n, forw_spec_proc_t *procs);

/* sp_type forw_spec_tprocs_mod_f forw_spec_tprocs_mod_count_f forw_spec_tprocs_mod_rearrange_db_f forw_spec_tprocs_mod_rearrange_ip_f */
typedef forw_spec_int_t forw_spec_tprocs_mod_f(forw_spec_elem_buf_t b, forw_spec_elem_index_t x, forw_spec_tproc_data_t tproc_data, forw_spec_proc_t *procs, forw_spec_elem_buf_t mod);
typedef void forw_spec_tprocs_mod_count_f(forw_spec_elem_t *s, forw_spec_tproc_data_t tproc_data, int *counts, forw_spec_proc_t *procs);
typedef void forw_spec_tprocs_mod_rearrange_db_f(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tproc_data_t tproc_data, int *displs, forw_spec_proc_t *procs, forw_spec_elem_t *mod);
typedef void forw_spec_tprocs_mod_rearrange_ip_f(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tproc_data_t tproc_data, int *displs, int *counts, forw_spec_int_t n, forw_spec_proc_t *procs, forw_spec_elem_t *mod);

/* sp_type forw_spec_tproc_reset_f */
typedef void forw_spec_tproc_reset_f(forw_spec_tproc_data_t tproc_data);


/* enable tloc features */
#ifdef forw_SPEC_TLOC

/* sp_macro forw_SPEC_TLOC forw_SPEC_LOC_NONE */


/* tloc rearrange */

/* sp_macro forw_SPEC_DECLARE_TLOC_REARRANGE_DB */
#define forw_SPEC_DECLARE_TLOC_REARRANGE_DB \
  struct { forw_spec_int_t i, p; } spec0d;

/* sp_macro forw_SPEC_DO_TLOC_REARRANGE_DB */
#define forw_SPEC_DO_TLOC_REARRANGE_DB(_tl_, _tld_, _sb_, _db_)  do { \
  for (spec0d.i = 0; spec0d.i < forw_spec_elem_get_n(_sb_); ++spec0d.i) { \
    spec0d.p = (_tl_)(forw_spec_elem_get_buf(_sb_), spec0d.i, _tld_); \
    if (spec0d.p == forw_SPEC_LOC_NONE) continue; \
    forw_spec_elem_copy_at((_sb_), spec0d.i, (_db_), spec0d.p); \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TLOC_REARRANGE_DB */
#define forw_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_rearrange_db(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tloc_data_t tloc_data) \
{ \
  forw_SPEC_DECLARE_TLOC_REARRANGE_DB \
  forw_SPEC_DO_TLOC_REARRANGE_DB(_tl_, tloc_data, s, d); \
}

/* sp_macro forw_SPEC_DECLARE_TLOC_REARRANGE_IP */
#define forw_SPEC_DECLARE_TLOC_REARRANGE_IP \
  struct { forw_spec_int_t i, p, np; } spec0i;

/* sp_macro forw_SPEC_DO_TLOC_REARRANGE_IP */
#define forw_SPEC_DO_TLOC_REARRANGE_IP(_tl_, _tld_, _b_, _xb_)  do { \
  for (spec0i.i = 0; spec0i.i < forw_spec_elem_get_n(_b_); ++spec0i.i) { \
    spec0i.p = (_tl_)(forw_spec_elem_get_buf(_b_), spec0i.i, _tld_); \
    if (spec0i.p == forw_SPEC_LOC_NONE) continue; \
    while (spec0i.i != spec0i.p) { \
      spec0i.np = (_tl_)(forw_spec_elem_get_buf(_b_), spec0i.p, _tld_); \
      if (spec0i.np == forw_SPEC_LOC_NONE) { forw_spec_elem_copy_at((_b_), spec0i.i, (_b_), spec0i.p); break; } \
      forw_spec_elem_exchange_at((_b_), spec0i.i, (_b_), spec0i.p, (_xb_)); \
      spec0i.p = spec0i.np; \
    } \
  } } while (0)

/* sp_macro forw_SPEC_FUNC_TLOC_REARRANGE_IP */
#define forw_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_rearrange_ip(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tloc_data_t tloc_data) \
{ \
  forw_SPEC_DECLARE_TLOC_REARRANGE_IP \
  forw_SPEC_DO_TLOC_REARRANGE_IP(_tl_, tloc_data, s, x); \
}


/* tloc_mod_mod rearrange */

/* sp_macro forw_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB */
#define forw_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  struct { forw_spec_int_t i, p; } spec1d;

/* sp_macro forw_SPEC_DO_TLOC_MOD_REARRANGE_DB */
#define forw_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, _tld_, _sb_, _db_, _ib_)  do { \
  if (_ib_) { \
    for (spec1d.i = 0; spec1d.i < forw_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(forw_spec_elem_get_buf(_sb_), spec1d.i, _tld_, forw_spec_elem_get_buf(_ib_)); \
      if (spec1d.p == forw_SPEC_LOC_NONE) continue; \
      forw_spec_elem_copy_at((_ib_), 0, (_db_), spec1d.p); \
    } \
  } else { \
    for (spec1d.i = 0; spec1d.i < forw_spec_elem_get_n(_sb_); ++spec1d.i) { \
      spec1d.p = (_tl_)(forw_spec_elem_get_buf(_sb_), spec1d.i, _tld_, NULL); \
      if (spec1d.p == forw_SPEC_LOC_NONE) continue; \
      forw_spec_elem_copy_at((_sb_), spec1d.i, (_db_), spec1d.p); \
    } \
  } } while (0) 

/* sp_macro forw_SPEC_FUNC_TLOC_MOD_REARRANGE_DB */
#define forw_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_...) \
_s_ void _name_##_tloc_mod_rearrange_db(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tloc_data_t tloc_data, forw_spec_elem_t *mod) \
{ \
  forw_SPEC_DECLARE_TLOC_MOD_REARRANGE_DB \
  forw_SPEC_DO_TLOC_MOD_REARRANGE_DB(_tl_, tloc_data, s, d, mod); \
}

/* sp_macro forw_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP */
#define forw_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  struct { forw_spec_int_t i, p, np; } spec1i;

/* sp_macro forw_SPEC_DO_TLOC_MOD_REARRANGE_IP */
#define forw_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, _tld_, _b_, _xb_, _ib_)  do { \
  if (_ib_) { \
    for (spec1i.i = 0; spec1i.i < forw_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(forw_spec_elem_get_buf(_b_), spec1i.i, _tld_, forw_spec_elem_get_buf(_ib_)); \
      if (spec1i.p == forw_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(forw_spec_elem_get_buf(_b_), spec1i.p, _tld_, forw_spec_elem_get_buf(_xb_)); \
        if (spec1i.np == forw_SPEC_LOC_NONE) break; \
        forw_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.p); \
        forw_spec_elem_copy_at((_xb_), 0, (_ib_), 0); \
        spec1i.p = spec1i.np; \
      } \
      forw_spec_elem_copy_at((_ib_), 0, (_b_), spec1i.i); \
    } \
  } else { \
    for (spec1i.i = 0; spec1i.i < forw_spec_elem_get_n(_b_); ++spec1i.i) { \
      spec1i.p = (_tl_)(forw_spec_elem_get_buf(_b_), spec1i.i, _tld_, NULL); \
      if (spec1i.p == forw_SPEC_LOC_NONE) continue; \
      while (spec1i.i != spec1i.p) { \
        spec1i.np = (_tl_)(forw_spec_elem_get_buf(_b_), spec1i.p, _tld_, NULL); \
        if (spec1i.np == forw_SPEC_LOC_NONE) { forw_spec_elem_copy_at((_b_), spec1i.i, (_b_), spec1i.p); break; } \
        forw_spec_elem_exchange_at((_b_), spec1i.i, (_b_), spec1i.p, (_xb_)); \
        spec1i.p = spec1i.np; \
      } \
    } \
 } } while (0) 

/* sp_macro forw_SPEC_FUNC_TLOC_MOD_REARRANGE_IP */
#define forw_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_) \
_s_ void _name_##_tloc_mod_rearrange_ip(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tloc_data_t tloc_data, forw_spec_elem_t *mod) \
{ \
  forw_SPEC_DECLARE_TLOC_MOD_REARRANGE_IP \
  forw_SPEC_DO_TLOC_MOD_REARRANGE_IP(_tl_, tloc_data, s, x, mod); \
}

/* sp_macro forw_SPEC_DEFINE_TLOC */
#define forw_SPEC_DEFINE_TLOC(_name_, _tl_, _s_...) \
  forw_SPEC_FUNC_TLOC_REARRANGE_DB(_name_, _tl_, _s_) \
  forw_SPEC_FUNC_TLOC_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro forw_SPEC_DEFINE_TLOC_MOD */
#define forw_SPEC_DEFINE_TLOC_MOD(_name_, _tl_, _s_...) \
  forw_SPEC_FUNC_TLOC_MOD_REARRANGE_DB(_name_, _tl_, _s_) \
  forw_SPEC_FUNC_TLOC_MOD_REARRANGE_IP(_name_, _tl_, _s_)

/* sp_macro forw_SPEC_EXT_PARAM_TLOC forw_SPEC_EXT_PARAM_TLOC_NULL forw_SPEC_EXT_PARAM_TLOC_MOD forw_SPEC_EXT_PARAM_TLOC_MOD_NULL */
#define forw_SPEC_EXT_PARAM_TLOC(_name_)      _name_##_tloc_rearrange_db, _name_##_tloc_rearrange_ip
#define forw_SPEC_EXT_PARAM_TLOC_NULL         NULL, NULL
#define forw_SPEC_EXT_PARAM_TLOC_MOD(_name_)  _name_##_tloc_mod_rearrange_db, _name_##_tloc_mod_rearrange_ip
#define forw_SPEC_EXT_PARAM_TLOC_MOD_NULL     NULL, NULL


/* sp_type forw_spec_tloc_f forw_spec_tloc_rearrange_db_f forw_spec_tloc_rearrange_ip_f */
typedef forw_spec_elem_index_t forw_spec_tloc_f(forw_spec_elem_buf_t b, forw_spec_elem_index_t x, forw_spec_tloc_data_t tloc_data);
typedef void forw_spec_tloc_rearrange_db_f(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tloc_data_t tloc_data);
typedef void forw_spec_tloc_rearrange_ip_f(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tloc_data_t tloc_data);

/* sp_type forw_spec_tloc_mod_f forw_spec_tloc_mod_rearrange_db_f forw_spec_tloc_mod_rearrange_ip_f */
typedef forw_spec_elem_index_t forw_spec_tloc_mod_f(forw_spec_elem_buf_t b, forw_spec_elem_index_t x, forw_spec_tloc_data_t tloc_data, forw_spec_elem_buf_t mod);
typedef void forw_spec_tloc_mod_rearrange_db_f(forw_spec_elem_t *s, forw_spec_elem_t *d, forw_spec_tloc_data_t tloc_data, forw_spec_elem_t *mod);
typedef void forw_spec_tloc_mod_rearrange_ip_f(forw_spec_elem_t *s, forw_spec_elem_t *x, forw_spec_tloc_data_t tloc_data, forw_spec_elem_t *mod);


#endif /* forw_SPEC_TLOC */






#ifdef SL_USE_MPI
# include <mpi.h>
#endif


/* sl_type forw_slint_t forw_slint */
typedef forw_sl_int_type_c forw_slint_t, forw_slint;  /* deprecated 'forw_slint' */

#define forw_slint_fmt   forw_sl_int_type_fmt    /* sl_macro */

/* sl_type forw_slindex_t */
typedef forw_sl_index_type_c forw_slindex_t;

#define forw_sindex_fmt  forw_sl_index_type_fmt  /* sl_macro */

/* sl_type forw_slkey_t */
typedef forw_sl_key_type_c forw_slkey_t;

/* sl_type forw_slkey_pure_t forw_slpkey_t */
typedef forw_sl_key_pure_type_c forw_slkey_pure_t, forw_slpkey_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type forw_sldata0_t */
#ifdef forw_sl_data0_type_c
typedef forw_sl_data0_type_c forw_sldata0_t;
#endif
/* sl_type forw_sldata1_t */
#ifdef forw_sl_data1_type_c
typedef forw_sl_data1_type_c forw_sldata1_t;
#endif
/* sl_type forw_sldata2_t */
#ifdef forw_sl_data2_type_c
typedef forw_sl_data2_type_c forw_sldata2_t;
#endif
/* sl_type forw_sldata3_t */
#ifdef forw_sl_data3_type_c
typedef forw_sl_data3_type_c forw_sldata3_t;
#endif
/* sl_type forw_sldata4_t */
#ifdef forw_sl_data4_type_c
typedef forw_sl_data4_type_c forw_sldata4_t;
#endif
/* sl_type forw_sldata5_t */
#ifdef forw_sl_data5_type_c
typedef forw_sl_data5_type_c forw_sldata5_t;
#endif
/* sl_type forw_sldata6_t */
#ifdef forw_sl_data6_type_c
typedef forw_sl_data6_type_c forw_sldata6_t;
#endif
/* sl_type forw_sldata7_t */
#ifdef forw_sl_data7_type_c
typedef forw_sl_data7_type_c forw_sldata7_t;
#endif
/* sl_type forw_sldata8_t */
#ifdef forw_sl_data8_type_c
typedef forw_sl_data8_type_c forw_sldata8_t;
#endif
/* sl_type forw_sldata9_t */
#ifdef forw_sl_data9_type_c
typedef forw_sl_data9_type_c forw_sldata9_t;
#endif
/* sl_type forw_sldata10_t */
#ifdef forw_sl_data10_type_c
typedef forw_sl_data10_type_c forw_sldata10_t;
#endif
/* sl_type forw_sldata11_t */
#ifdef forw_sl_data11_type_c
typedef forw_sl_data11_type_c forw_sldata11_t;
#endif
/* sl_type forw_sldata12_t */
#ifdef forw_sl_data12_type_c
typedef forw_sl_data12_type_c forw_sldata12_t;
#endif
/* sl_type forw_sldata13_t */
#ifdef forw_sl_data13_type_c
typedef forw_sl_data13_type_c forw_sldata13_t;
#endif
/* sl_type forw_sldata14_t */
#ifdef forw_sl_data14_type_c
typedef forw_sl_data14_type_c forw_sldata14_t;
#endif
/* sl_type forw_sldata15_t */
#ifdef forw_sl_data15_type_c
typedef forw_sl_data15_type_c forw_sldata15_t;
#endif
/* sl_type forw_sldata16_t */
#ifdef forw_sl_data16_type_c
typedef forw_sl_data16_type_c forw_sldata16_t;
#endif
/* sl_type forw_sldata17_t */
#ifdef forw_sl_data17_type_c
typedef forw_sl_data17_type_c forw_sldata17_t;
#endif
/* sl_type forw_sldata18_t */
#ifdef forw_sl_data18_type_c
typedef forw_sl_data18_type_c forw_sldata18_t;
#endif
/* sl_type forw_sldata19_t */
#ifdef forw_sl_data19_type_c
typedef forw_sl_data19_type_c forw_sldata19_t;
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

/* sl_type forw_slweight_t */
typedef forw_sl_weight_type_c forw_slweight_t;

#define forw_slweight_fmt  forw_sl_weight_type_fmt  /* sl_macro */

#if defined(forw_sl_elem_weight) && defined(forw_sl_weight_intequiv)
typedef forw_sl_weight_type_c forw_slcount_t;       /* sl_type forw_slcount_t */
# define forw_slcount_fmt  forw_sl_weight_type_fmt  /* sl_macro */
#else
typedef forw_sl_int_type_c forw_slcount_t;
# define forw_slcount_fmt  forw_sl_int_type_fmt
#endif


/* sl_type forw__slpwkey_t forw_slpwkey_t */
typedef struct forw__slpwkey_t
{
  forw_slpkey_t pkey;
  forw_slweight_t weight;

} forw_slpwkey_t;


/* sl_type forw__elements_t forw_elements_t */
typedef struct forw__elements_t
{
  forw_slint_t size, max_size;
  forw_slkey_t *keys;

#ifdef forw_SL_INDEX
  forw_slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef forw_SL_DATA0
  forw_sldata0_t *data0;
#endif
#ifdef forw_SL_DATA1
  forw_sldata1_t *data1;
#endif
#ifdef forw_SL_DATA2
  forw_sldata2_t *data2;
#endif
#ifdef forw_SL_DATA3
  forw_sldata3_t *data3;
#endif
#ifdef forw_SL_DATA4
  forw_sldata4_t *data4;
#endif
#ifdef forw_SL_DATA5
  forw_sldata5_t *data5;
#endif
#ifdef forw_SL_DATA6
  forw_sldata6_t *data6;
#endif
#ifdef forw_SL_DATA7
  forw_sldata7_t *data7;
#endif
#ifdef forw_SL_DATA8
  forw_sldata8_t *data8;
#endif
#ifdef forw_SL_DATA9
  forw_sldata9_t *data9;
#endif
#ifdef forw_SL_DATA10
  forw_sldata10_t *data10;
#endif
#ifdef forw_SL_DATA11
  forw_sldata11_t *data11;
#endif
#ifdef forw_SL_DATA12
  forw_sldata12_t *data12;
#endif
#ifdef forw_SL_DATA13
  forw_sldata13_t *data13;
#endif
#ifdef forw_SL_DATA14
  forw_sldata14_t *data14;
#endif
#ifdef forw_SL_DATA15
  forw_sldata15_t *data15;
#endif
#ifdef forw_SL_DATA16
  forw_sldata16_t *data16;
#endif
#ifdef forw_SL_DATA17
  forw_sldata17_t *data17;
#endif
#ifdef forw_SL_DATA18
  forw_sldata18_t *data18;
#endif
#ifdef forw_SL_DATA19
  forw_sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} forw_elements_t;


/* sl_type forw__packed_element_t forw_packed_element_t */
typedef struct forw__packed_element_t
{
  forw_slkey_t key;

#ifdef forw_SL_PACKED_INDEX
  forw_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef forw_SL_DATA0
# ifdef forw_sl_data0_flex
  forw_sldata0_t data0[];
# else
  forw_sldata0_t data0[forw_sl_data0_size_c];
# endif
#endif
#ifdef forw_SL_DATA1
# ifdef forw_sl_data1_flex
  forw_sldata1_t data1[];
# else
  forw_sldata1_t data1[forw_sl_data1_size_c];
# endif
#endif
#ifdef forw_SL_DATA2
# ifdef forw_sl_data2_flex
  forw_sldata2_t data2[];
# else
  forw_sldata2_t data2[forw_sl_data2_size_c];
# endif
#endif
#ifdef forw_SL_DATA3
# ifdef forw_sl_data3_flex
  forw_sldata3_t data3[];
# else
  forw_sldata3_t data3[forw_sl_data3_size_c];
# endif
#endif
#ifdef forw_SL_DATA4
# ifdef forw_sl_data4_flex
  forw_sldata4_t data4[];
# else
  forw_sldata4_t data4[forw_sl_data4_size_c];
# endif
#endif
#ifdef forw_SL_DATA5
# ifdef forw_sl_data5_flex
  forw_sldata5_t data5[];
# else
  forw_sldata5_t data5[forw_sl_data5_size_c];
# endif
#endif
#ifdef forw_SL_DATA6
# ifdef forw_sl_data6_flex
  forw_sldata6_t data6[];
# else
  forw_sldata6_t data6[forw_sl_data6_size_c];
# endif
#endif
#ifdef forw_SL_DATA7
# ifdef forw_sl_data7_flex
  forw_sldata7_t data7[];
# else
  forw_sldata7_t data7[forw_sl_data7_size_c];
# endif
#endif
#ifdef forw_SL_DATA8
# ifdef forw_sl_data8_flex
  forw_sldata8_t data8[];
# else
  forw_sldata8_t data8[forw_sl_data8_size_c];
# endif
#endif
#ifdef forw_SL_DATA9
# ifdef forw_sl_data9_flex
  forw_sldata9_t data9[];
# else
  forw_sldata9_t data9[forw_sl_data9_size_c];
# endif
#endif
#ifdef forw_SL_DATA10
# ifdef forw_sl_data10_flex
  forw_sldata10_t data10[];
# else
  forw_sldata10_t data10[forw_sl_data10_size_c];
# endif
#endif
#ifdef forw_SL_DATA11
# ifdef forw_sl_data11_flex
  forw_sldata11_t data11[];
# else
  forw_sldata11_t data11[forw_sl_data11_size_c];
# endif
#endif
#ifdef forw_SL_DATA12
# ifdef forw_sl_data12_flex
  forw_sldata12_t data12[];
# else
  forw_sldata12_t data12[forw_sl_data12_size_c];
# endif
#endif
#ifdef forw_SL_DATA13
# ifdef forw_sl_data13_flex
  forw_sldata13_t data13[];
# else
  forw_sldata13_t data13[forw_sl_data13_size_c];
# endif
#endif
#ifdef forw_SL_DATA14
# ifdef forw_sl_data14_flex
  forw_sldata14_t data14[];
# else
  forw_sldata14_t data14[forw_sl_data14_size_c];
# endif
#endif
#ifdef forw_SL_DATA15
# ifdef forw_sl_data15_flex
  forw_sldata15_t data15[];
# else
  forw_sldata15_t data15[forw_sl_data15_size_c];
# endif
#endif
#ifdef forw_SL_DATA16
# ifdef forw_sl_data16_flex
  forw_sldata16_t data16[];
# else
  forw_sldata16_t data16[forw_sl_data16_size_c];
# endif
#endif
#ifdef forw_SL_DATA17
# ifdef forw_sl_data17_flex
  forw_sldata17_t data17[];
# else
  forw_sldata17_t data17[forw_sl_data17_size_c];
# endif
#endif
#ifdef forw_SL_DATA18
# ifdef forw_sl_data18_flex
  forw_sldata18_t data18[];
# else
  forw_sldata18_t data18[forw_sl_data18_size_c];
# endif
#endif
#ifdef forw_SL_DATA19
# ifdef forw_sl_data19_flex
  forw_sldata19_t data19[];
# else
  forw_sldata19_t data19[forw_sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} forw_packed_element_t;


/* sl_type forw__packed_elements_t forw_packed_elements_t */
typedef struct forw__packed_elements_t
{
  forw_slint_t size, max_size;
  
  forw_packed_element_t *elements;
  
} forw_packed_elements_t;


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


/* sl_type forw__classification_info_t forw_classification_info_t forw_classification_info */
typedef struct forw__classification_info_t
{
  forw_slint_t nclasses;
  forw_slkey_pure_t *keys;
  forw_slint_t *counts;
  forw_slint_t *masks;

  /* */
  forw_slint_t *all_local_sizes;
  forw_slint_t *local_lt_eq_counts;
  forw_slint_t *all_local_lt_eq_counts;

} forw_classification_info_t, forw_classification_info;  /* deprecated 'forw_classification_info' */


/* key2class, sl_type forw_key2class_f */
typedef forw_slint_t (*forw_key2class_f)(forw_slkey_t *, forw_slint, void *);

/* pivot-element, sl_type forw_pivot_f */
typedef forw_slint_t (*forw_pivot_f)(forw_elements_t *);

/* sorting-network, sl_type forw_sortnet_f forw_sortnet_data_t */
typedef void *forw_sortnet_data_t;
typedef forw_slint_t (*forw_sortnet_f)(forw_slint_t size, forw_slint_t rank, forw_slint_t stage, forw_sortnet_data_t snd, forw_slint_t *up);

/* merge2, sl_type forw_merge2x_f forw_merge2X_f */
typedef forw_slint_t (*forw_merge2x_f)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
typedef forw_slint_t (*forw_merge2X_f)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx, forw_elements_t *t);

/* sl_type forw__permute_generic_t forw_permute_generic_t */
typedef struct forw__permute_generic_t
{
  int type;

  forw_spec_tloc_f *tloc;
  forw_spec_tloc_rearrange_db_f *tloc_rearrange_db;
  forw_spec_tloc_rearrange_ip_f *tloc_rearrange_ip;

  forw_spec_tloc_mod_f *tloc_mod;
  forw_spec_tloc_mod_rearrange_db_f *tloc_mod_rearrange_db;
  forw_spec_tloc_mod_rearrange_ip_f *tloc_mod_rearrange_ip;

} forw_permute_generic_t;

/* sl_macro forw_PERMUTE_GENERIC_DEFINE_TLOC forw_PERMUTE_GENERIC_INIT_TLOC forw_PERMUTE_GENERIC_INIT_EXT_TLOC */
#define forw_PERMUTE_GENERIC_DEFINE_TLOC(_tl_, _s_...)      forw_SPEC_DEFINE_TLOC(_tl_, _tl_, _s_)
#define forw_PERMUTE_GENERIC_INIT_TLOC(_tl_)                { 1, _tl_, forw_SPEC_EXT_PARAM_TLOC_NULL,  NULL, forw_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define forw_PERMUTE_GENERIC_INIT_EXT_TLOC(_tl_)            { 1, _tl_, forw_SPEC_EXT_PARAM_TLOC(_tl_), NULL, forw_SPEC_EXT_PARAM_TLOC_MOD_NULL }

/* sl_macro forw_PERMUTE_GENERIC_DEFINE_TLOC_MOD forw_PERMUTE_GENERIC_INIT_TLOC_MOD forw_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD */
#define forw_PERMUTE_GENERIC_DEFINE_TLOC_MOD(_tl_, _s_...)  forw_SPEC_DEFINE_TLOC_MOD(_tl_, _tl_, _s_)
#define forw_PERMUTE_GENERIC_INIT_TLOC_MOD(_tl_)            { 2, NULL, forw_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, forw_SPEC_EXT_PARAM_TLOC_MOD_NULL }
#define forw_PERMUTE_GENERIC_INIT_EXT_TLOC_MOD(_tl_)        { 2, NULL, forw_SPEC_EXT_PARAM_TLOC_MOD_NULL, _tl_, forw_SPEC_EXT_PARAM_TLOC_MOD(_tl_) }

/* sl_type forw__split_generic_t forw_split_generic_t */
typedef struct forw__split_generic_t
{
  int type;

  forw_spec_tproc_f *tproc;
  forw_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  forw_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  forw_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  forw_spec_tproc_mod_f *tproc_mod;
  forw_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  forw_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  forw_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  forw_spec_tprocs_f *tprocs;
  forw_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  forw_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  forw_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  forw_spec_tprocs_mod_f *tprocs_mod;
  forw_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  forw_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  forw_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

  forw_spec_tproc_reset_f *reset;

} forw_split_generic_t;

/* sl_macro forw_SPLIT_GENERIC_DEFINE_TPROC forw_SPLIT_GENERIC_INIT_TPROC forw_SPLIT_GENERIC_INIT_EXT_TPROC */
#define forw_SPLIT_GENERIC_DEFINE_TPROC(_tp_, _s_...)         forw_SPEC_DEFINE_TPROC(_tp_, _tp_, _s_)
#define forw_SPLIT_GENERIC_INIT_TPROC(_tp_, _r_...)           { 1, _tp_, forw_SPEC_EXT_PARAM_TPROC_NULL,  NULL, forw_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, forw_SPEC_EXT_PARAM_TPROCS_NULL, NULL, forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define forw_SPLIT_GENERIC_INIT_EXT_TPROC(_tp_, _r_...)       { 1, _tp_, forw_SPEC_EXT_PARAM_TPROC(_tp_), NULL, forw_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, forw_SPEC_EXT_PARAM_TPROCS_NULL, NULL, forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro forw_SPLIT_GENERIC_DEFINE_TPROC_MOD forw_SPLIT_GENERIC_INIT_TPROC_MOD forw_SPLIT_GENERIC_INIT_EXT_TPROC_MOD */
#define forw_SPLIT_GENERIC_DEFINE_TPROC_MOD(_tp_, _s_...)     forw_SPEC_DEFINE_TPROC_MOD(_tp_, _tp_, _s_)
#define forw_SPLIT_GENERIC_INIT_TPROC_MOD(_tp_, _r_...)       { 2, NULL, forw_SPEC_EXT_PARAM_TPROC_NULL, _tp_, forw_SPEC_EXT_PARAM_TPROC_MOD_NULL,  NULL, forw_SPEC_EXT_PARAM_TPROCS_NULL, NULL, forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define forw_SPLIT_GENERIC_INIT_EXT_TPROC_MOD(_tp_, _r_...)   { 2, NULL, forw_SPEC_EXT_PARAM_TPROC_NULL, _tp_, forw_SPEC_EXT_PARAM_TPROC_MOD(_tp_), NULL, forw_SPEC_EXT_PARAM_TPROCS_NULL, NULL, forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro forw_SPLIT_GENERIC_DEFINE_TPROCS forw_SPLIT_GENERIC_INIT_TPROCS forw_SPLIT_GENERIC_INIT_EXT_TPROCS */
#define forw_SPLIT_GENERIC_DEFINE_TPROCS(_tp_, _s_...)        forw_SPEC_DEFINE_TPROCS(_tp_, _tp_, _s_)
#define forw_SPLIT_GENERIC_INIT_TPROCS(_tp_, _r_...)          { 3, NULL, forw_SPEC_EXT_PARAM_TPROC_NULL, NULL, forw_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, forw_SPEC_EXT_PARAM_TPROCS_NULL,  NULL, forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define forw_SPLIT_GENERIC_INIT_EXT_TPROCS(_tp_, _r_...)      { 3, NULL, forw_SPEC_EXT_PARAM_TPROC_NULL, NULL, forw_SPEC_EXT_PARAM_TPROC_MOD_NULL, _tp_, forw_SPEC_EXT_PARAM_TPROCS(_tp_), NULL, forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }

/* sl_macro forw_SPLIT_GENERIC_DEFINE_TPROCS_MOD forw_SPLIT_GENERIC_INIT_TPROCS_MOD forw_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD */
#define forw_SPLIT_GENERIC_DEFINE_TPROCS_MOD(_tp_, _s_...)    forw_SPEC_DEFINE_TPROCS_MOD(_tp_, _tp_, _s_)
#define forw_SPLIT_GENERIC_INIT_TPROCS_MOD(_tp_, _r_...)      { 4, NULL, forw_SPEC_EXT_PARAM_TPROC_NULL, NULL, forw_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, forw_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL,  _r_ }
#define forw_SPLIT_GENERIC_INIT_EXT_TPROCS_MOD(_tp_, _r_...)  { 4, NULL, forw_SPEC_EXT_PARAM_TPROC_NULL, NULL, forw_SPEC_EXT_PARAM_TPROC_MOD_NULL, NULL, forw_SPEC_EXT_PARAM_TPROCS_NULL,  _tp_, forw_SPEC_EXT_PARAM_TPROCS_MOD(_tp_), _r_ }

/* sl_type forw_tloc_f forw_tloc_mod_f */
typedef forw_slint_t forw_tloc_f(forw_elements_t *b, forw_slint_t x, void *tloc_data);
typedef forw_slint_t forw_tloc_mod_f(forw_elements_t *b, forw_slint_t x, void *tloc_data, forw_elements_t *mod);

/* sl_type forw_tproc_f forw_tproc_mod_f forw_tprocs_f forw_tprocs_mod_f */
typedef int forw_tproc_f(forw_elements_t *b, forw_slint_t x, void *tproc_data);
typedef int forw_tproc_mod_f(forw_elements_t *b, forw_slint_t x, void *tproc_data, forw_elements_t *mod);
typedef forw_slint_t forw_tprocs_f(forw_elements_t *b, forw_slint_t x, void *tproc_data, int *procs);
typedef forw_slint_t forw_tprocs_mod_f(forw_elements_t *b, forw_slint_t x, void *tproc_data, int *procs, forw_elements_t *mod);

/* sl_type forw_tproc_reset_f */
typedef void forw_tproc_reset_f(void *tproc_data);

/* sl_macro forw_TPROC_RESET_NULL */
#define forw_TPROC_RESET_NULL  NULL

/* sl_type forw__tproc_t forw_tproc_t */
typedef struct forw__tproc_t *forw_tproc_t;

/* sl_type forw__tproc_exdef forw_tproc_exdef */
typedef struct forw__tproc_exdef {
  int type;

  forw_spec_tproc_count_f *tproc_count_db, *tproc_count_ip;
  forw_spec_tproc_rearrange_db_f *tproc_rearrange_db;
  forw_spec_tproc_rearrange_ip_f *tproc_rearrange_ip;

  forw_spec_tproc_mod_count_f *tproc_mod_count_db, *tproc_mod_count_ip;
  forw_spec_tproc_mod_rearrange_db_f *tproc_mod_rearrange_db;
  forw_spec_tproc_mod_rearrange_ip_f *tproc_mod_rearrange_ip;

  forw_spec_tprocs_count_f *tprocs_count_db, *tprocs_count_ip;
  forw_spec_tprocs_rearrange_db_f *tprocs_rearrange_db;
  forw_spec_tprocs_rearrange_ip_f *tprocs_rearrange_ip;

  forw_spec_tprocs_mod_count_f *tprocs_mod_count_db, *tprocs_mod_count_ip;
  forw_spec_tprocs_mod_rearrange_db_f *tprocs_mod_rearrange_db;
  forw_spec_tprocs_mod_rearrange_ip_f *tprocs_mod_rearrange_ip;

} const *forw_tproc_exdef;

/* sl_macro forw_TPROC_EXDEF_NULL */
#define forw_TPROC_EXDEF_NULL  NULL

/* sl_macro forw_TPROC_EXDEF_DEFINE_TPROC forw_TPROC_EXDEF_DEFINE_TPROC_MOD forw_TPROC_EXDEF_DEFINE_TPROCS forw_TPROC_EXDEF_DEFINE_TPROCS_MOD */
#define forw_TPROC_EXDEF_DEFINE_TPROC(_name_, _tp_, _s_...) \
  forw_SPEC_DEFINE_TPROC(_name_, _tp_, _s_) \
  _s_ const struct forw__tproc_exdef _##_name_ = { 1, forw_SPEC_EXT_PARAM_TPROC(_name_), forw_SPEC_EXT_PARAM_TPROC_MOD_NULL, forw_SPEC_EXT_PARAM_TPROCS_NULL, forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define forw_TPROC_EXDEF_DEFINE_TPROC_MOD(_name_, _tp_, _s_...) \
  forw_SPEC_DEFINE_TPROC_MOD(_name_, _tp_, _s_) \
  _s_ const struct forw__tproc_exdef _##_name_ = { 2, forw_SPEC_EXT_PARAM_TPROC_NULL, forw_SPEC_EXT_PARAM_TPROC_MOD(_name_), forw_SPEC_EXT_PARAM_TPROCS_NULL, forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define forw_TPROC_EXDEF_DEFINE_TPROCS(_name_, _tp_, _s_...) \
  forw_SPEC_DEFINE_TPROCS(_name_, _tp_, _s_) \
  _s_ const struct forw__tproc_exdef _##_name_ = { 3, forw_SPEC_EXT_PARAM_TPROC_NULL, forw_SPEC_EXT_PARAM_TPROC_MOD_NULL, forw_SPEC_EXT_PARAM_TPROCS(_name_), forw_SPEC_EXT_PARAM_TPROCS_MOD_NULL }, *_name_ = &_##_name_;

#define forw_TPROC_EXDEF_DEFINE_TPROCS_MOD(_name_, _tp_, _s_...) \
  forw_SPEC_DEFINE_TPROCS_MOD(_name_, _tp_, _s_) \
  _s_ const struct forw__tproc_exdef _##_name_ = { 4, forw_SPEC_EXT_PARAM_TPROC_NULL, forw_SPEC_EXT_PARAM_TPROC_MOD_NULL, forw_SPEC_EXT_PARAM_TPROCS_NULL, forw_SPEC_EXT_PARAM_TPROCS_MOD(_name_) }, *_name_ = &_##_name_;


/* deprecated, sl_type forw_k2c_func forw_pivot_func forw_sn_func forw_m2x_func forw_m2X_func */
typedef forw_key2class_f forw_k2c_func;
typedef forw_pivot_f forw_pivot_func;
typedef forw_sortnet_f forw_sn_func;
typedef forw_merge2x_f forw_m2x_func;
typedef forw_merge2X_f forw_m2X_func;


/* sl_type forw__mergek_t forw_mergek_t */
typedef struct forw__mergek_t
{
  forw_sortnet_f sn;
  forw_sortnet_data_t snd;

  forw_merge2x_f m2x;
  forw_elements_t *sx;

} forw_mergek_t;


/* sl_type forw_keys_init_type_t forw_keys_init_data_t */
typedef forw_slint_t forw_keys_init_type_t;
typedef void *forw_keys_init_data_t;

/* sl_type forw_key_set_data_t forw_key_set_f */
typedef void *forw_key_set_data_t;
typedef void (*forw_key_set_f)(forw_slkey_pure_t *k, forw_key_set_data_t d);


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


/* forw_elements_keys_stats */
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


/* partition conditions, sl_type forw__partcond2_t forw_partcond2_t */
typedef struct forw__partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} forw_partcond2_t;


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

/* partition conditions, sl_type forw__partcond_t forw_partcond_t forw_partcond_p */
typedef struct forw__partcond_t
{
  forw_slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} forw_partcond_t, *forw_partcond_p;


/* internal partition conditions, sl_type forw__partcond_intern_t forw_partcond_intern_t forw_partcond_intern_p */
typedef struct forw__partcond_intern_t
{
  forw_slint_t pcm;
  forw_slint_t count_min, count_max;
  forw_slint_t count_low, count_high;
#ifdef elem_weight
  forw_slweight_t weight_min, weight_max;
  forw_slweight_t weight_low, weight_high;
#endif

} forw_partcond_intern_t, *forw_partcond_intern_p;


/* sl_type forw__parttype_t forw_parttype_t forw_parttype_p */
typedef struct forw__parttype_t
{
  forw_slint_t type;

} forw_parttype_t, *forw_parttype_p;


/* generic binning method */

/* sl_type forw__bin_t forw_bin_t */
typedef struct forw__bin_t
{
  forw_elements_t s;

#ifdef elem_weight
  forw_slweight_t weight;
#endif

} forw_bin_t;


/* sl_type forw__splitter_t forw_splitter_t */
typedef struct forw__splitter_t
{
  forw_slint_t n;

  int *displs;
  forw_slkey_pure_t *s;
  forw_slint_t *sn;

} forw_splitter_t;


struct forw__binning_t;

/* sl_type forw_binning_pre_f forw_binning_exec_f forw_binning_refine_f forw_binning_hit_f forw_binning_finalize_f forw_binning_post_f */
typedef forw_slint_t (*forw_binning_pre_f)(struct forw__binning_t *bm);
typedef forw_slint_t (*forw_binning_exec_f)(struct forw__binning_t *bm, forw_bin_t *bin, forw_slcount_t *counts, forw_slweight_t *weights);
typedef forw_slint_t (*forw_binning_refine_f)(struct forw__binning_t *bm, forw_bin_t *bin, forw_slint_t k, forw_slcount_t *counts, forw_slweight_t *weights, forw_splitter_t *sp, forw_slint_t s, forw_bin_t *new_bin);
typedef forw_slint_t (*forw_binning_hit_f)(struct forw__binning_t *bm, forw_bin_t *bin, forw_slint_t k, forw_slcount_t *counts, forw_splitter_t *sp, forw_slint_t s);
typedef forw_slint_t (*forw_binning_finalize_f)(struct forw__binning_t *bm, forw_bin_t *bin, forw_slint_t dc, forw_slweight_t dw, forw_slint_t lc_min, forw_slint_t lc_max, forw_slcount_t *lcs, forw_slweight_t *lws, forw_splitter_t *sp, forw_slint_t s);
typedef forw_slint_t (*forw_binning_post_f)(struct forw__binning_t *bm);


/* sl_type forw__binning_data_t forw_binning_data_t */
typedef union forw__binning_data_t
{
  struct
  {
    forw_slint_t rhigh, rlow, rwidth;
    forw_slint_t rcurrent;
    forw_slkey_pure_t bit_mask;

    forw_elements_t sx;

  } radix;

} forw_binning_data_t;


/* sl_type forw__binning_t forw_binning_t */
typedef struct forw__binning_t
{
  forw_slint_t nbins, max_nbins;
  
  forw_binning_pre_f pre;
  forw_binning_exec_f exec;
  forw_binning_refine_f refine;
  forw_binning_hit_f hit;
  forw_binning_finalize_f finalize;
  forw_binning_post_f post;

  forw_slint_t sorted;

  forw_slint_t docounts;
#ifdef elem_weight
  forw_slint_t doweights;
#endif

  forw_binning_data_t bd;

} forw_binning_t;


/* sl_type forw__local_bins_t forw_local_bins_t */
typedef struct forw__local_bins_t
{
  forw_binning_t *bm;

  forw_slint_t nbins, max_nbins;
  forw_slint_t nelements;

  forw_slint_t docounts;
#ifdef elem_weight
  forw_slint_t doweights;
#endif

  forw_slint_t nbinnings, max_nbinnings;

  forw_slint_t nbins_new, last_new_b, last_new_k;
  forw_bin_t *bins, *bins_new;
  forw_bin_t *bins0, *bins1;

  forw_slint_t *bcws;

#if defined(elem_weight) && defined(forw_sl_weight_intequiv)
  forw_slint_t cw_factor, w_index, bin_cw_factor;
  forw_slweight_t *cws, *bin_cws;
  forw_slweight_t *prefix_cws;
#else
  forw_slint_t *cs, *bin_cs;
  forw_slint_t *prefix_cs;
# ifdef elem_weight
  forw_slweight_t *ws, *bin_ws;
  forw_slweight_t *prefix_ws;
# endif
#endif

  forw_slint_t last_exec_b;

} forw_local_bins_t;


/* sl_type forw__global_bins_t forw_global_bins_t */
typedef struct forw__global_bins_t
{
  forw_binning_t *bm;
  
  forw_local_bins_t lb;

  forw_slint_t *bcws;

#if defined(elem_weight) && defined(forw_sl_weight_intequiv)
  forw_slweight_t *cws;
  forw_slweight_t *prefix_cws;
#else
  forw_slint_t *cs;
  forw_slint_t *prefix_cs;
# ifdef elem_weight
  forw_slweight_t *ws;
  forw_slweight_t *prefix_ws;
# endif
#endif

} forw_global_bins_t;


/* sl_type forw_rti_cmc_t */
typedef struct
{
  forw_slint_t cmp, movek, moved;

} forw_rti_cmc_t;

#ifndef my_rti_ccmp
# define my_rti_ccmp(m)    m.cmc.cmp
# define my_rti_cmovek(m)  m.cmc.movek
# define my_rti_cmoved(m)  m.cmc.moved
#endif


/* sl_type forw_rti_tim_t */
typedef struct
{
  double start, stop;
  double last, cumu;

  forw_slint_t num;

} forw_rti_tim_t[rti_tids];

#ifndef my_rti_tlast
# define my_rti_tlast(m, t)  m.tim[t].last
# define my_rti_tcumu(m, t)  m.tim[t].cumu
# define my_rti_tnum(m, t)   m.tim[t].num
#endif


/* sl_type forw_rti_mem_t */
typedef struct
{
  forw_slint_t nalloc, nfree;
  forw_slint_t max, cur, cur_max;

} forw_rti_mem_t;


/* sl_type forw_rti_t */
typedef struct
{
  /* compare-move-counter */
  forw_rti_cmc_t cmc;
  /* timer */
  forw_rti_tim_t tim;
  /* memory */
  forw_rti_mem_t mem;

} forw_rti_t;

#ifndef my_rti_reset
# define my_rti_reset(m)  memset((void *) &m, 0, sizeof(m))
#endif


/* sl_type forw__sl_context_t forw_sl_context_t */
typedef struct forw__sl_context_t
{

/* src/base/base.c */
  struct {
int dummy_rank;
  } sl;
#ifdef forw_SL_USE_RTI
forw_rti_t rti;
#endif
  struct {
forw_slint_t ip_threshold;
forw_slint_t db_threshold;
forw_slint_t ma_threshold;
  } sr;
  struct {
forw_slint_t threshold;
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
forw_slint_t sendrecv_replace_memsize;
forw_slint_t sendrecv_replace_mpi_maxsize;
  } me;
#endif
#ifdef SL_USE_MPI
  struct {
double t[2];
forw_slint_t max_nprocs;
forw_slint_t packed;
forw_slint_t minalloc;
double overalloc;
  } meas;
#endif
#ifdef SL_USE_MPI
  struct {
forw_slint_t packed;
forw_slint_t db_packed;
forw_slint_t ip_packed;
  } mea;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef forw_MSEG_ROOT
int root;
#endif
#ifdef forw_MSEG_BORDER_UPDATE_REDUCTION
double border_update_count_reduction;
double border_update_weight_reduction;
#endif
#ifdef forw_MSEG_FORWARD_ONLY
forw_slint_t forward_only;
#endif
#ifdef forw_MSEG_INFO
forw_slint_t info_rounds;
forw_slint_t *info_finish_rounds;
double info_finish_rounds_avg;
forw_slweight_t info_total_weights;
#endif
forw_slint_t binnings;
forw_slint_t finalize_mode;
  } mseg;
#endif
#ifdef SL_USE_MPI
  struct {
#ifdef forw_MSS_ROOT
int root;
#endif
  } mss;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
forw_slint_t sync;
  } msm;
#endif
#ifdef SL_USE_MPI
  struct {
double t[4];
forw_slint_t sync;
forw_partcond_t *r_pc;
  } msp;
#endif
#ifdef SL_USE_MPI
  struct {
double i_t[3];
double p_t[3];
double b_t[3];
forw_slint_t sync;
forw_slint_t i_sync;
forw_slint_t p_sync;
forw_slint_t b_sync;
forw_slint_t back_packed;
  } mssp;
#endif
} forw_sl_context_t;






/* sl_macro forw_elem_set_size forw_elem_set_max_size forw_elem_set_keys forw_elem_set_indices */
#define forw_elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define forw_elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define forw_elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define forw_elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))
/* DATAX_TEMPLATE_BEGIN */
#define forw_elem_set_data0(_e_, _b_)     ((_e_)->data0 = (_b_))  /* sl_macro */
#define forw_elem_set_data1(_e_, _b_)     ((_e_)->data1 = (_b_))  /* sl_macro */
#define forw_elem_set_data2(_e_, _b_)     ((_e_)->data2 = (_b_))  /* sl_macro */
#define forw_elem_set_data3(_e_, _b_)     ((_e_)->data3 = (_b_))  /* sl_macro */
#define forw_elem_set_data4(_e_, _b_)     ((_e_)->data4 = (_b_))  /* sl_macro */
#define forw_elem_set_data5(_e_, _b_)     ((_e_)->data5 = (_b_))  /* sl_macro */
#define forw_elem_set_data6(_e_, _b_)     ((_e_)->data6 = (_b_))  /* sl_macro */
#define forw_elem_set_data7(_e_, _b_)     ((_e_)->data7 = (_b_))  /* sl_macro */
#define forw_elem_set_data8(_e_, _b_)     ((_e_)->data8 = (_b_))  /* sl_macro */
#define forw_elem_set_data9(_e_, _b_)     ((_e_)->data9 = (_b_))  /* sl_macro */
#define forw_elem_set_data10(_e_, _b_)     ((_e_)->data10 = (_b_))  /* sl_macro */
#define forw_elem_set_data11(_e_, _b_)     ((_e_)->data11 = (_b_))  /* sl_macro */
#define forw_elem_set_data12(_e_, _b_)     ((_e_)->data12 = (_b_))  /* sl_macro */
#define forw_elem_set_data13(_e_, _b_)     ((_e_)->data13 = (_b_))  /* sl_macro */
#define forw_elem_set_data14(_e_, _b_)     ((_e_)->data14 = (_b_))  /* sl_macro */
#define forw_elem_set_data15(_e_, _b_)     ((_e_)->data15 = (_b_))  /* sl_macro */
#define forw_elem_set_data16(_e_, _b_)     ((_e_)->data16 = (_b_))  /* sl_macro */
#define forw_elem_set_data17(_e_, _b_)     ((_e_)->data17 = (_b_))  /* sl_macro */
#define forw_elem_set_data18(_e_, _b_)     ((_e_)->data18 = (_b_))  /* sl_macro */
#define forw_elem_set_data19(_e_, _b_)     ((_e_)->data19 = (_b_))  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro forw_elem_get_size forw_elem_get_max_size forw_elem_get_keys forw_elem_get_indices */
#define forw_elem_get_size(_e_)           (_e_)->size
#define forw_elem_get_max_size(_e_)       (_e_)->max_size
#define forw_elem_get_keys(_e_)           (_e_)->keys
#define forw_elem_get_indices(_e_)        (_e_)->indices
/* DATAX_TEMPLATE_BEGIN */
#define forw_elem_get_data0(_e_)          (_e_)->data0  /* sl_macro */
#define forw_elem_get_data1(_e_)          (_e_)->data1  /* sl_macro */
#define forw_elem_get_data2(_e_)          (_e_)->data2  /* sl_macro */
#define forw_elem_get_data3(_e_)          (_e_)->data3  /* sl_macro */
#define forw_elem_get_data4(_e_)          (_e_)->data4  /* sl_macro */
#define forw_elem_get_data5(_e_)          (_e_)->data5  /* sl_macro */
#define forw_elem_get_data6(_e_)          (_e_)->data6  /* sl_macro */
#define forw_elem_get_data7(_e_)          (_e_)->data7  /* sl_macro */
#define forw_elem_get_data8(_e_)          (_e_)->data8  /* sl_macro */
#define forw_elem_get_data9(_e_)          (_e_)->data9  /* sl_macro */
#define forw_elem_get_data10(_e_)          (_e_)->data10  /* sl_macro */
#define forw_elem_get_data11(_e_)          (_e_)->data11  /* sl_macro */
#define forw_elem_get_data12(_e_)          (_e_)->data12  /* sl_macro */
#define forw_elem_get_data13(_e_)          (_e_)->data13  /* sl_macro */
#define forw_elem_get_data14(_e_)          (_e_)->data14  /* sl_macro */
#define forw_elem_get_data15(_e_)          (_e_)->data15  /* sl_macro */
#define forw_elem_get_data16(_e_)          (_e_)->data16  /* sl_macro */
#define forw_elem_get_data17(_e_)          (_e_)->data17  /* sl_macro */
#define forw_elem_get_data18(_e_)          (_e_)->data18  /* sl_macro */
#define forw_elem_get_data19(_e_)          (_e_)->data19  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro forw_elem_set_block forw_elem_set_block_size forw_elem_get_block forw_elem_get_block_size */
#define forw_elem_set_block(_e_, _b_)       ((_e_)->keys = (_b_), (_e_)->max_size = -1)
#define forw_elem_set_block_size(_e_, _s_)  ((_e_)->size = (_s_))
#define forw_elem_get_block(_e_)            ((void *) (((_e_)->max_size < 0)?(_e_)->keys:NULL))
#define forw_elem_get_block_size(_e_)       (_e_)->size

/* sl_macro forw_pelem_set_size forw_pelem_set_max_size forw_pelem_set_elements */
#define forw_pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define forw_pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define forw_pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))

/* sl_macro forw_pelem_get_size forw_pelem_get_max_size forw_pelem_get_elements */
#define forw_pelem_get_size(_e_)           (_e_)->size
#define forw_pelem_get_max_size(_e_)       (_e_)->max_size
#define forw_pelem_get_elements(_e_)       (_e_)->elements

/* sl_macro forw_SL_DEFCON */
#define forw_SL_DEFCON(_v_)  (forw_sl_default_context._v_)






/* src/base/base.c */
extern forw_sl_context_t forw_sl_default_context;
extern const int forw_default_sl_dummy_rank;
#ifdef forw_SL_USE_RTI
extern const forw_rti_t forw_default_rti;
#endif
extern const forw_slint_t forw_default_sr_ip_threshold;
extern const forw_slint_t forw_default_sr_db_threshold;
extern const forw_slint_t forw_default_sr_ma_threshold;
extern const forw_slint_t forw_default_sri_threshold;

/* src/base_mpi/base_mpi.c */
#ifdef SL_USE_MPI
extern const MPI_Datatype forw_default_mpi_int_datatype;
extern const MPI_Datatype forw_default_mpi_key_datatype;
extern const MPI_Datatype forw_default_mpi_pkey_datatype;
extern const MPI_Datatype forw_default_mpi_pwkey_datatype;
extern const MPI_Datatype forw_default_mpi_index_datatype;
extern const MPI_Datatype forw_default_mpi_weight_datatype;
extern const MPI_Datatype forw_default_mpi_data_datatype[];
extern const int forw_default_mpi_rank;
#endif
extern const void *forw_default_me_sendrecv_replace_mem;
extern const forw_slint_t forw_default_me_sendrecv_replace_memsize;
extern const forw_slint_t forw_default_me_sendrecv_replace_mpi_maxsize;
extern const double forw_default_meas_t[];
extern const forw_slint_t forw_default_meas_max_nprocs;
extern const forw_slint_t forw_default_meas_packed;
extern const forw_slint_t forw_default_meas_minalloc;
extern const double forw_default_meas_overalloc;
extern const forw_slint_t forw_default_mea_packed;
extern const forw_slint_t forw_default_mea_db_packed;
extern const forw_slint_t forw_default_mea_ip_packed;
#ifdef forw_MSEG_ROOT
extern const int forw_default_mseg_root;
#endif
#ifdef forw_MSEG_BORDER_UPDATE_REDUCTION
extern const double forw_default_mseg_border_update_count_reduction;
extern const double forw_default_mseg_border_update_weight_reduction;
#endif
#ifdef forw_MSEG_FORWARD_ONLY
extern const forw_slint_t forw_default_mseg_forward_only;
#endif
#ifdef forw_MSEG_INFO
extern const forw_slint_t forw_default_mseg_info_rounds;
extern const forw_slint_t *forw_default_mseg_info_finish_rounds;
extern const double forw_default_mseg_info_finish_rounds_avg;
extern const forw_slweight_t forw_default_mseg_info_total_weights;
#endif
extern const forw_slint_t forw_default_mseg_binnings;
extern const forw_slint_t forw_default_mseg_finalize_mode;
#ifdef forw_MSS_ROOT
extern const int forw_default_mss_root;
#endif
extern const double forw_default_msm_t[];
extern const forw_slint_t forw_default_msm_sync;
extern const double forw_default_msp_t[];
extern const forw_slint_t forw_default_msp_sync;
extern const forw_partcond_t *forw_default_msp_r_pc;
extern const double forw_default_mssp_i_t[];
extern const double forw_default_mssp_p_t[];
extern const double forw_default_mssp_b_t[];
extern const forw_slint_t forw_default_mssp_sync;
extern const forw_slint_t forw_default_mssp_i_sync;
extern const forw_slint_t forw_default_mssp_p_sync;
extern const forw_slint_t forw_default_mssp_b_sync;
extern const forw_slint_t forw_default_mssp_back_packed;






/* src/base/base.c */
forw_slint_t SL_PROTO(forw_binning_create)(forw_local_bins_t *lb, forw_slint_t max_nbins, forw_slint_t max_nbinnings, forw_elements_t *s, forw_slint_t nelements, forw_slint_t docounts, forw_slint_t doweights, forw_binning_t *bm);
forw_slint_t SL_PROTO(forw_binning_destroy)(forw_local_bins_t *lb);
forw_slint_t SL_PROTO(forw_binning_pre)(forw_local_bins_t *lb);
forw_slint_t SL_PROTO(forw_binning_exec_reset)(forw_local_bins_t *lb, forw_slint_t do_bins, forw_slint_t do_prefixes);
forw_slint_t SL_PROTO(forw_binning_exec)(forw_local_bins_t *lb, forw_slint_t b, forw_slint_t do_bins, forw_slint_t do_prefixes);
forw_slint_t SL_PROTO(forw_binning_refine)(forw_local_bins_t *lb, forw_slint_t b, forw_slint_t k, forw_splitter_t *sp, forw_slint_t s);
forw_slint_t SL_PROTO(forw_binning_hit)(forw_local_bins_t *lb, forw_slint_t b, forw_slint_t k, forw_splitter_t *sp, forw_slint_t s);
forw_slint_t SL_PROTO(forw_binning_finalize)(forw_local_bins_t *lb, forw_slint_t b, forw_slint_t dc, forw_slweight_t dw, forw_slint_t lc_min, forw_slint_t lc_max, forw_slcount_t *lcs, forw_slweight_t *lws, forw_splitter_t *sp, forw_slint_t s);
forw_slint_t SL_PROTO(forw_binning_post)(forw_local_bins_t *lb);
forw_slint_t SL_PROTO(forw_binning_radix_create)(forw_binning_t *bm, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, forw_slint_t sorted);
forw_slint_t SL_PROTO(forw_binning_radix_destroy)(forw_binning_t *bm);
forw_slint_t SL_PROTO(forw_binning_radix_pre)(forw_binning_t *bm);
forw_slint_t SL_PROTO(forw_binning_radix_exec)(forw_binning_t *bm, forw_bin_t *bin, forw_slcount_t *counts, forw_slweight_t *weights);
forw_slint_t SL_PROTO(forw_binning_radix_refine)(forw_binning_t *bm, forw_bin_t *bin, forw_slint_t k, forw_slcount_t *counts, forw_slweight_t *weights, forw_splitter_t *sp, forw_slint_t s, forw_bin_t *new_bin);
forw_slint_t SL_PROTO(forw_binning_radix_hit)(forw_binning_t *bm, forw_bin_t *bin, forw_slint_t k, forw_slcount_t *counts, forw_splitter_t *sp, forw_slint_t s);
forw_slint_t SL_PROTO(forw_binning_radix_finalize)(forw_binning_t *bm, forw_bin_t *bin, forw_slint_t dc, forw_slweight_t dw, forw_slint_t lc_min, forw_slint_t lc_max, forw_slcount_t *lcs, forw_slweight_t *lws, forw_splitter_t *sp, forw_slint_t s);
forw_slint_t SL_PROTO(forw_binning_radix_post)(forw_binning_t *bm);
forw_slint_t SL_PROTO(forw_elements_alloc)(forw_elements_t *s, forw_slint_t nelements, slcint_t components);
forw_slint_t SL_PROTO(forw_elements_free)(forw_elements_t *s);
forw_slint_t SL_PROTO(forw_elements_realloc)(forw_elements_t *s, forw_slint_t nelements, slcint_t components);
forw_slint_t SL_PROTO(forw_elements_alloca)(forw_elements_t *s, forw_slint_t nelements, slcint_t components);
forw_slint_t SL_PROTO(forw_elements_freea)(forw_elements_t *s);
forw_slint_t SL_PROTO(forw_elements_alloc_from_blocks)(forw_elements_t *s, forw_slint_t nblocks, void **blocks, forw_slint_t *blocksizes, forw_slint_t alignment, forw_slint_t nmax, slcint_t components);
forw_slint_t SL_PROTO(forw_elements_alloc_from_block)(forw_elements_t *s, void *block, forw_slint_t blocksize, forw_slint_t alignment, forw_slint_t nmax, slcint_t components);
forw_slint_t SL_PROTO(forw_elements_alloc_block)(forw_elements_t *s, void **block, forw_slint_t *blocksize, forw_slint_t alignment, forw_slint_t maxblocksize);
forw_slint_t SL_PROTO(forw_elements_copy)(forw_elements_t *s, forw_elements_t *d);
forw_slint_t SL_PROTO(forw_elements_copy_at)(forw_elements_t *s, forw_slint_t sat, forw_elements_t *d, forw_slint_t dat);
forw_slint_t SL_PROTO(forw_elements_ncopy)(forw_elements_t *s, forw_elements_t *d, forw_slint_t n);
forw_slint_t SL_PROTO(forw_elements_nmove)(forw_elements_t *s, forw_elements_t *d, forw_slint_t n);
forw_slint_t SL_PROTO(forw_elements_printf)(forw_elements_t *s, const char *prefix);
forw_slint_t SL_PROTO(forw_elements_extract)(forw_elements_t *src, forw_slint_t nelements, forw_elements_t *dst0, forw_elements_t *dst1);
forw_slint_t SL_PROTO(forw_elements_touch)(forw_elements_t *s);
forw_slint_t SL_PROTO(forw_elements_digest_sum)(forw_elements_t *s, forw_slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(forw_elements_crc32)(forw_elements_t *s, forw_slint nelements, forw_slint_t keys, forw_slint_t data);
forw_slint_t SL_PROTO(forw_elements_digest_hash)(forw_elements_t *s, forw_slint_t nelements, slcint_t components, void *hash);
forw_slint_t SL_PROTO(forw_elements_random_exchange)(forw_elements_t *s, forw_slint_t rounds, forw_elements_t *xs);
forw_slint_t SL_PROTO(forw_elements_keys_init_seed)(unsigned long s);
forw_slint_t SL_PROTO(forw_elements_keys_init)(forw_elements_t *s, forw_keys_init_type_t t, forw_keys_init_data_t d);
forw_slint_t SL_PROTO(forw_elements_keys_init_randomized)(forw_elements_t *s, forw_slint_t nkeys, forw_keys_init_type_t t, forw_keys_init_data_t d);
forw_slint_t SL_PROTO(forw_elements_keys_init_from_file)(forw_elements_t *s, forw_slint_t data, char *filename, forw_slint_t from, forw_slint_t to, forw_slint_t const_bytes_per_line);
forw_slint_t SL_PROTO(forw_elements_keys_save_to_file)(forw_elements_t *s, char *filename);
forw_slint_t SL_PROTO(forw_elements_validate_order)(forw_elements_t *s, forw_slint_t n);
forw_slint_t SL_PROTO(forw_elements_validate_order_bmask)(forw_elements_t *s, forw_slint_t n, forw_slkey_pure_t bmask);
forw_slint_t SL_PROTO(forw_elements_validate_order_weight)(forw_elements_t *s, forw_slint_t n, forw_slkey_pure_t weight);
forw_slint_t SL_PROTO(forw_elements_keys_stats)(forw_elements_t *s, forw_slkey_pure_t *stats);
forw_slint_t SL_PROTO(forw_elements_keys_stats_print)(forw_elements_t *s);
forw_slint_t SL_PROTO(forw_elements_print_keys)(forw_elements_t *s);
forw_slint_t SL_PROTO(forw_elements_print_all)(forw_elements_t *s);
forw_slweight_t SL_PROTO(forw_elements_get_weight)(forw_elements_t *s);
forw_slint_t SL_PROTO(forw_elements_get_minmax_keys)(forw_elements_t *s, forw_slint_t nelements, forw_slkey_pure_t *minmaxkeys);
forw_slint_t SL_PROTO(forw_elements_alloc_packed)(forw_packed_elements_t *s, forw_slint_t nelements);
forw_slint_t SL_PROTO(forw_elements_free_packed)(forw_packed_elements_t *s);
forw_slint_t SL_PROTO(forw_elements_alloc_packed_from_block)(forw_packed_elements_t *s, void *block, forw_slint_t blocksize, forw_slint_t alignment, forw_slint_t nmax);
forw_slint_t SL_PROTO(forw_elements_pack_indexed)(forw_elements_t *s, forw_packed_elements_t *d, forw_slindex_t *rindx, forw_slindex_t *windx);
forw_slint_t SL_PROTO(forw_elements_pack)(forw_elements_t *s, forw_packed_elements_t *d);
forw_slint_t SL_PROTO(forw_elements_pack_at)(forw_elements_t *s, forw_slint_t sat, forw_packed_elements_t *d, forw_slint_t dat);
forw_slint_t SL_PROTO(forw_elements_unpack_indexed)(forw_packed_elements_t *s, forw_elements_t *d, forw_slindex_t *rindx, forw_slindex_t *windx);
forw_slint_t SL_PROTO(forw_elements_unpack)(forw_packed_elements_t *s, forw_elements_t *d);
forw_slint_t SL_PROTO(forw_elements_unpack_at)(forw_packed_elements_t *s, forw_slint_t sat, forw_elements_t *d, forw_slint_t dat);
forw_slint_t SL_PROTO(forw_elements_unpack_keys)(forw_packed_elements_t *s, forw_slkey_t *k);
forw_slint SL_PROTO(forw_merge2_basic_auto_01_x)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint SL_PROTO(forw_merge2_basic_01_x)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx, forw_m2x_func _x0_1, forw_m2x_func _0x_1);
forw_slint SL_PROTO(forw_merge2_basic_01_X)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx, forw_elements_t *t, forw_m2X_func _X0_1, forw_m2X_func _0X_1);
forw_slint SL_PROTO(forw_merge2_simplify_s1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx, forw_slint s1elements);
forw_slint SL_PROTO(forw_merge2_memory_adaptive)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_compo_hula)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *xs);
forw_slint_t SL_PROTO(forw_merge2_basic_sseq_x0_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_sseq_0x_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_sseq_01_x)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_sseq_01)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_merge2_basic_sbin_x0_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_sbin_0x_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_sbin_01_x)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_sbin_01)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_merge2_basic_shyb_x0_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_shyb_0x_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_shyb_01_x)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_shyb_01)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_merge2_basic_straight_x0_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_straight_0x_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_straight_01_x)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_straight_x_0_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_merge2_basic_straight_X0_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_merge2_basic_straight_0X_1)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_merge2_basic_straight_01_X)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_merge2_basic_straight_X0_1u)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_merge2_compo_tridgell)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_mergep_2way_ip_int)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t p, int *displs, forw_merge2x_f m2x);
forw_slint_t SL_PROTO(forw_mergep_2way_ip_int_rec)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t p, int *displs, forw_merge2x_f m2x);
forw_slint_t SL_PROTO(forw_mergep_heap_int)(forw_elements_t *s, forw_elements_t *d, forw_slint_t p, int *displs, int *counts);
forw_slint_t SL_PROTO(forw_mergep_heap_int_idx)(forw_elements_t *s, forw_elements_t *d, forw_slint_t p, int *displs, int *counts);
forw_slint_t SL_PROTO(forw_mergep_heap_idx)(forw_elements_t *s, forw_elements_t *d, forw_slint_t p, forw_slindex_t *displs, forw_slindex_t *counts);
forw_slint_t SL_PROTO(forw_mergep_heap_unpack_idx)(forw_packed_elements_t *s, forw_elements_t *d, forw_slint_t p, forw_slindex_t *displs, forw_slindex_t *counts);
forw_slint_t SL_PROTO(forw_mergep_heap_unpack_idxonly)(forw_packed_elements_t *s, forw_elements_t *d, forw_slint_t p, forw_slindex_t *displs, forw_slindex_t *counts);
forw_slint_t SL_PROTO(forw_permute_generic_db)(forw_elements_t *s, forw_elements_t *d, forw_permute_generic_t *pg, void *pg_data);
forw_slint_t SL_PROTO(forw_permute_generic_ip)(forw_elements_t *s, forw_elements_t *x, forw_permute_generic_t *pg, void *pg_data);
forw_slint SL_PROTO(forw_sl_search_sequential_lt)(forw_elements_t *s, forw_slpkey_t k);
forw_slint SL_PROTO(forw_sl_search_sequential_le)(forw_elements_t *s, forw_slpkey_t k);
forw_slint SL_PROTO(forw_sl_search_sequential_gt)(forw_elements_t *s, forw_slpkey_t k);
forw_slint SL_PROTO(forw_sl_search_sequential_ge)(forw_elements_t *s, forw_slpkey_t k);
forw_slint SL_PROTO(forw_sl_search_p_sequential_lt)(forw_elements_t *s, forw_slpkey_t *k);
forw_slint SL_PROTO(forw_sl_search_p_sequential_le)(forw_elements_t *s, forw_slpkey_t *k);
forw_slint SL_PROTO(forw_sl_search_p_sequential_gt)(forw_elements_t *s, forw_slpkey_t *k);
forw_slint SL_PROTO(forw_sl_search_p_sequential_ge)(forw_elements_t *s, forw_slpkey_t *k);
forw_slint SL_PROTO(forw_sl_search_binary_lt)(forw_elements_t *s, forw_slpkey_t k);
forw_slint SL_PROTO(forw_sl_search_binary_le)(forw_elements_t *s, forw_slpkey_t k);
forw_slint SL_PROTO(forw_sl_search_binary_gt)(forw_elements_t *s, forw_slpkey_t k);
forw_slint SL_PROTO(forw_sl_search_binary_ge)(forw_elements_t *s, forw_slpkey_t k);
forw_slint SL_PROTO(forw_sl_search_p_binary_lt)(forw_elements_t *s, forw_slpkey_t *k);
forw_slint SL_PROTO(forw_sl_search_p_binary_le)(forw_elements_t *s, forw_slpkey_t *k);
forw_slint SL_PROTO(forw_sl_search_p_binary_gt)(forw_elements_t *s, forw_slpkey_t *k);
forw_slint SL_PROTO(forw_sl_search_p_binary_ge)(forw_elements_t *s, forw_slpkey_t *k);
forw_slint_t SL_PROTO(forw_sl_search_binary_lt_bmask)(forw_elements_t *s, forw_slpkey_t k, forw_slpkey_t bmask);
forw_slint_t SL_PROTO(forw_sl_search_binary_le_bmask)(forw_elements_t *s, forw_slpkey_t k, forw_slpkey_t bmask);
forw_slint_t SL_PROTO(forw_sl_search_binary_sign_switch)(forw_elements_t *s);
forw_slint SL_PROTO(forw_sl_search_hybrid_lt)(forw_elements_t *s, forw_slpkey_t k, forw_slint t);
forw_slint SL_PROTO(forw_sl_search_hybrid_le)(forw_elements_t *s, forw_slpkey_t k, forw_slint t);
forw_slint SL_PROTO(forw_sl_search_hybrid_gt)(forw_elements_t *s, forw_slpkey_t k, forw_slint t);
forw_slint SL_PROTO(forw_sl_search_hybrid_ge)(forw_elements_t *s, forw_slpkey_t k, forw_slint t);
forw_slint SL_PROTO(forw_sl_search_p_hybrid_lt)(forw_elements_t *s, forw_slpkey_t *k, forw_slint t);
forw_slint SL_PROTO(forw_sl_search_p_hybrid_le)(forw_elements_t *s, forw_slpkey_t *k, forw_slint t);
forw_slint SL_PROTO(forw_sl_search_p_hybrid_gt)(forw_elements_t *s, forw_slpkey_t *k, forw_slint t);
forw_slint SL_PROTO(forw_sl_search_p_hybrid_ge)(forw_elements_t *s, forw_slpkey_t *k, forw_slint t);
forw_slint SL_PROTO(forw_ilog2c)(forw_slint x);
forw_slint SL_PROTO(forw_ilog2f)(forw_slint x);
forw_slint SL_PROTO(forw_print_bits)(forw_slint v);
forw_slint SL_PROTO(forw_pivot_random)(forw_elements_t *s);
forw_slint_t SL_PROTO(forw_counts2displs)(forw_slint_t n, int *counts, int *displs);
forw_slint_t SL_PROTO(forw_displs2counts)(forw_slint_t n, int *displs, int *counts, forw_slint_t total_counts);
void SL_PROTO(forw_get_displcounts_extent)(forw_slint_t n, int *displs, int *counts, forw_slint_t *lb, forw_slint_t *extent);
void SL_PROTO(forw_elem_set_data)(forw_elements_t *e, ...);
forw_slint_t SL_PROTO(forw_elem_get_max_byte)();
forw_slint_t SL_PROTO(forw_elem_reverse)(forw_elements_t *e, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_elem_nxchange_at)(forw_elements_t *e0, forw_slint_t at0, forw_elements_t *e1, forw_slint_t at1, forw_slint_t n, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_elem_nxchange)(forw_elements_t *e0, forw_elements_t *e1, forw_slint_t n, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_elem_nxchange_ro0)(forw_elements_t *e0, forw_elements_t *e1, forw_slint_t n, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_elem_rotate)(forw_elements_t *e, forw_slint_t m, forw_slint_t n, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_elem_rotate_ro0)(forw_elements_t *e, forw_slint_t m, forw_slint_t n, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_elem_rotate_ro1)(forw_elements_t *e, forw_slint_t m, forw_slint_t n, forw_elements_t *t);
forw_slint_t SL_PROTO(forw_sort_counting_use_displs)(forw_elements_t *s, forw_elements_t *d, forw_slint_t ndispls, forw_slint_t *displs);
forw_slint_t SL_PROTO(forw_sort_counting_use_counts)(forw_elements_t *s, forw_elements_t *d, forw_slint_t ncounts, forw_slint_t *counts);
forw_slint_t SL_PROTO(forw_sort_counting_get_counts)(forw_elements_t *s, forw_elements_t *d, forw_slint_t ncounts, forw_slint_t *counts);
forw_slint_t SL_PROTO(forw_sort_counting)(forw_elements_t *s, forw_elements_t *d, forw_slint_t ncounts);
forw_slint SL_PROTO(forw_sort_heap)(forw_elements_t *s, forw_elements_t *xs);
forw_slint_t SL_PROTO(forw_sort_insert_bmask_kernel)(forw_elements_t *s, forw_elements_t *sx, forw_slkey_pure_t bmask);
forw_slint_t SL_PROTO(forw_sort_insert)(forw_elements_t *s, forw_elements_t *sx);
forw_slint_t SL_PROTO(forw_sort_permute_forward)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t *perm, forw_slint_t offset, forw_slint_t mask_bit);
forw_slint_t SL_PROTO(forw_sort_permute_backward)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t *perm, forw_slint_t offset, forw_slint_t mask_bit);
forw_slint SL_PROTO(forw_sort_quick)(forw_elements_t *s, forw_elements_t *xs);
forw_slint_t SL_PROTO(forw_sort_radix_ip)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth);
forw_slint_t SL_PROTO(forw_sort_radix_db)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth);
forw_slint_t SL_PROTO(forw_sort_radix_ma)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth);
forw_slint_t SL_PROTO(forw_sort_radix)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth);
forw_slint_t SL_PROTO(forw_sort_radix_1bit_kernel)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t rhigh, forw_slint_t rlow);
forw_slint SL_PROTO(forw_sort_radix_1bit)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t rhigh, forw_slint_t rlow);
forw_slint_t SL_PROTO(forw_sort_radix_iter)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t presorted, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth);
forw_slint SL_PROTO(forw_sn_hypercube_lh)(forw_slint size, forw_slint rank, forw_slint stage, void *snp, forw_slint *up);
forw_slint SL_PROTO(forw_sn_hypercube_hl)(forw_slint size, forw_slint rank, forw_slint stage, void *snp, forw_slint *up);
forw_slint SL_PROTO(forw_sn_odd_even_trans)(forw_slint size, forw_slint rank, forw_slint stage, void *snp, forw_slint *up);
forw_slint SL_PROTO(forw_sn_odd)(forw_slint size, forw_slint rank, forw_slint stage, void *snp, forw_slint *up);
forw_slint SL_PROTO(forw_sn_even)(forw_slint size, forw_slint rank, forw_slint stage, void *snp, forw_slint *up);
forw_slint SL_PROTO(forw_sn_batcher)(forw_slint size, forw_slint rank, forw_slint stage, void *snp, forw_slint *up);
forw_slint SL_PROTO(forw_sn_bitonic)(forw_slint size, forw_slint rank, forw_slint stage, void *snp, forw_slint *up);
forw_slint SL_PROTO(forw_sn_connected)(forw_slint size, forw_slint rank, forw_slint stage, void *snp, forw_slint *up);
forw_slint_t SL_PROTO(forw_split_generic_db)(forw_elements_t *s, forw_elements_t *d, forw_split_generic_t *sg, void *sg_data, forw_slint_t n);
forw_slint_t SL_PROTO(forw_split_generic_ip)(forw_elements_t *s, forw_elements_t *d, forw_split_generic_t *sg, void *sg_data, forw_slint_t n);
forw_slint_t SL_PROTO(forw_split_generic_count_db)(forw_elements_t *s, forw_split_generic_t *sg, void *sg_data, int *counts, forw_slint_t n);
forw_slint_t SL_PROTO(forw_split_generic_count_ip)(forw_elements_t *s, forw_split_generic_t *sg, void *sg_data, int *counts, forw_slint_t n);
forw_slint_t SL_PROTO(forw_split_generic_rearrange_db)(forw_elements_t *s, forw_elements_t *d, forw_split_generic_t *sg, void *sg_data, int *counts, forw_slint_t n);
forw_slint_t SL_PROTO(forw_split_generic_rearrange_ip)(forw_elements_t *s, forw_elements_t *d, forw_split_generic_t *sg, void *sg_data, int *counts, int *displs, forw_slint_t n);
forw_slint_t SL_PROTO(forw_splitter_reset)(forw_splitter_t *sp);
forw_slint_t SL_PROTO(forw_splitx_radix)(forw_elements_t *s, forw_elements_t *sx, forw_slint_t nclasses, forw_slint_t shl, forw_slint_t *counts);
forw_slint SL_PROTO(forw_split2_lt_ge)(forw_elements_t *s, forw_slkey_pure_t *k, forw_elements_t *t);
forw_slint SL_PROTO(forw_split2_le_gt)(forw_elements_t *s, forw_slkey_pure_t *k, forw_elements_t *t);
forw_slint SL_PROTO(forw_split3_lt_eq_gt)(forw_elements_t *s, forw_slkey_pure_t *k, forw_elements_t *t, forw_slint *nlt, forw_slint *nle);
forw_slint SL_PROTO(forw_split3_lt_eq_gt_old)(forw_elements_t *s, forw_slkey_pure_t *k, forw_elements_t *t, forw_slint *nlt, forw_slint *nle);
forw_slint SL_PROTO(forw_split2_b)(forw_elements_t *s, forw_elements_t *sx, forw_slkey_pure_t bmask);
forw_slint SL_PROTO(forw_splitk_k2c_af)(forw_elements_t *s, forw_elements_t *sx, forw_slint k, forw_slint *c, forw_k2c_func k2c, void *k2c_data);
forw_slint SL_PROTO(forw_splitk_k2c)(forw_elements_t *s, forw_elements_t *sx, forw_slint k, forw_slint *c, forw_k2c_func k2c, void *k2c_data);
forw_slint SL_PROTO(forw_splitk_k2c_count)(forw_elements_t *s, forw_slint k, forw_slint *c, forw_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI





/* src/base_mpi/base_mpi.c */
forw_slint_t SL_PROTO(forw_mpi_binning_create)(forw_global_bins_t *gb, forw_slint_t max_nbins, forw_slint_t max_nbinnings, forw_elements_t *s, forw_slint_t nelements, forw_slint_t docounts, forw_slint_t doweights, forw_binning_t *bm, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_binning_destroy)(forw_global_bins_t *gb, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_binning_pre)(forw_global_bins_t *gb, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_binning_exec_reset)(forw_global_bins_t *gb, forw_slint_t do_bins, forw_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_binning_exec_local)(forw_global_bins_t *gb, forw_slint_t b, forw_slint_t do_bins, forw_slint_t do_prefixes, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_binning_exec_global)(forw_global_bins_t *gb, forw_slint_t do_bins, forw_slint_t do_prefixes, forw_slint_t root, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_binning_refine)(forw_global_bins_t *gb, forw_slint_t b, forw_slint_t k, forw_splitter_t *sp, forw_slint_t s, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_binning_hit)(forw_global_bins_t *gb, forw_slint_t b, forw_slint_t k, forw_splitter_t *sp, forw_slint_t s, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_binning_finalize)(forw_global_bins_t *gb, forw_slint_t b, forw_slint_t dc, forw_slweight_t dw, forw_slint_t lc_min, forw_slint_t lc_max, forw_slcount_t *lcs, forw_slweight_t *lws, forw_splitter_t *sp, forw_slint_t s, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_binning_post)(forw_global_bins_t *gb, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_datatypes_init)();
forw_slint_t SL_PROTO(forw_mpi_datatypes_release)();
forw_slint_t SL_PROTO(forw_mpi_get_grid_properties)(forw_slint_t ndims, forw_slint_t *dims, forw_slint_t *pos, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_subgroups_create)(forw_slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_subgroups_delete)(forw_slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);
int SL_PROTO(forw_sl_MPI_Allreduce)(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank);
int SL_PROTO(forw_sl_MPI_Alltoall_int)(void *sendbuf, int sendcount, void *recvbuf, int recvcount, MPI_Comm comm, int size, int rank);
forw_slint_t SL_PROTO(forw_mpi_elements_keys_init_from_file)(forw_elements_t *s, char *filename, forw_slint from, forw_slint to, forw_slint const_bytes_per_line, forw_slint root, int size, int rank, MPI_Comm comm);
forw_slint SL_PROTO(forw_mpi_elements_validate_order)(forw_elements_t *s, forw_slint n, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_linear_exchange_pure_keys)(forw_slkey_pure_t *in, forw_slkey_pure_t *out, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_check_order)(forw_elements_t *s, forw_slint_t nelements, forw_slint_t *orders, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_check_global_order)(forw_slkey_pure_t local_min, forw_slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_digest_sum)(forw_elements_t *s, forw_slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(forw_mpi_elements_crc32)(forw_elements_t *s, forw_slint_t n, forw_slint_t keys, forw_slint_t data, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_digest_hash)(forw_elements_t *s, forw_slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_get_counts)(forw_elements_t *s, forw_slint_t *clocal, forw_slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
forw_slweight_t SL_PROTO(forw_mpi_elements_get_weights)(forw_elements_t *s, forw_slweight_t *wlocal, forw_slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_get_counts_and_weights)(forw_elements_t *s, forw_slint_t nelements, forw_slint_t *counts, forw_slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_sendrecv_replace)(forw_elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_tproc_create_tproc)(forw_tproc_t *tproc, forw_tproc_f *tfn, forw_tproc_reset_f *rfn, forw_tproc_exdef exdef);
forw_slint_t SL_PROTO(forw_tproc_create_tproc_mod)(forw_tproc_t *tproc, forw_tproc_mod_f *tfn, forw_tproc_reset_f *rfn, forw_tproc_exdef exdef);
forw_slint_t SL_PROTO(forw_tproc_create_tprocs)(forw_tproc_t *tproc, forw_tprocs_f *tfn, forw_tproc_reset_f *rfn, forw_tproc_exdef exdef);
forw_slint_t SL_PROTO(forw_tproc_create_tprocs_mod)(forw_tproc_t *tproc, forw_tprocs_mod_f *tfn, forw_tproc_reset_f *rfn, forw_tproc_exdef exdef);
forw_slint_t SL_PROTO(forw_tproc_free)(forw_tproc_t *tproc);
forw_slint_t SL_PROTO(forw_tproc_set_proclists)(forw_tproc_t *tproc, forw_slint_t nsend_procs, int *send_procs, forw_slint_t nrecv_procs, int *recv_procs, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_tproc_verify)(forw_tproc_t tproc, void *data, forw_elements_t *s, int proc);
forw_slint_t SL_PROTO(forw_mpi_elements_alltoall_specific)(forw_elements_t *sin, forw_elements_t *sout, forw_elements_t *xs, forw_tproc_t tproc, void *data, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_alltoallv_db_packed)(forw_elements_t *sbuf, int *scounts, int *sdispls, forw_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_alltoallv_db)(forw_elements_t *sbuf, int *scounts, int *sdispls, forw_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_alltoallv_ip_packed)(forw_elements_t *s, forw_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_alltoallv_ip_double)(forw_elements_t *s, forw_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_alltoallv_ip_mpi)(forw_elements_t *s, forw_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_alltoallv_ip_dash)(forw_elements_t *s, forw_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_alltoallv_ip)(forw_elements_t *s, forw_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_alltoallv_proclists_db)(forw_elements_t *sbuf, int *scounts, int *sdispls, int nsendprocs, int *sendprocs, forw_elements_t *rbuf, int *rcounts, int *rdispls, int nrecvprocs, int *recvprocs, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, forw_slint_t structured);
forw_slint_t SL_PROTO(forw_mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);
forw_slint_t SL_PROTO(forw_mpi_find_exact_equal)(forw_elements_t *s, forw_slint_t other_rank, forw_slint_t high_rank, forw_slint_t *ex_start, forw_slint_t *ex_size, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_find_exact)(forw_elements_t *s, forw_slint_t other_rank, forw_slint_t high_rank, forw_slint_t *dst_size, forw_slint_t *ex_start, forw_slint_t *ex_sizes, forw_slint_t *nx_move, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_merge2)(forw_elements_t *s, forw_slint_t other_rank, forw_slint_t high_rank, forw_slint_t *dst_size, forw_merge2x_f m2, forw_elements_t *xs, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_mergek_equal)(forw_elements_t *s, forw_sortnet_f sn, forw_sortnet_data_t snd, forw_merge2x_f m2x, forw_elements_t *xs, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_mergek_sorted)(forw_elements_t *s, forw_merge2x_f m2x, forw_elements_t *xs, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_mergek_sorted2)(forw_elements_t *s, forw_sortnet_f sn, forw_sortnet_data_t snd, forw_merge2x_f m2x, forw_elements_t *xs, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_mergek)(forw_elements_t *s, forw_sortnet_f sn, forw_sortnet_data_t snd, forw_merge2x_f m2x, forw_elements_t *xs, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_mergek_equal2)(forw_elements_t *s, forw_sortnet_f sn, forw_sortnet_data_t snd, forw_merge2x_f m2x, forw_elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);
forw_slint_t SL_PROTO(forw_mpi_partition_exact_generic)(forw_elements_t *s, forw_partcond_t *pcond, forw_binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_partition_exact_radix)(forw_elements_t *s, forw_partcond_t *pcond, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, forw_slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_partition_exact_radix_ngroups)(forw_elements_t *s, forw_partcond_t *pcond, forw_slint_t ngroups, MPI_Comm *group_comms, forw_elements_t *sx, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_partition_exact_radix_2groups)(forw_elements_t *s, forw_partcond_t *pcond, MPI_Comm group_comm, forw_elements_t *sx, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_partition_sample_regular)(forw_elements_t *s, forw_partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_rebalance)(forw_elements_t *s0, forw_elements_t *s1, forw_slint_t stable, forw_slint_t *dst_size, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_rebalance_alltoallv)(forw_elements_t *sbuf, int *scounts, int *sdispls, forw_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
void SL_PROTO(forw_mpi_partcond_set_even)(forw_partcond_t *pcond, forw_slint_t pcm, forw_slint_t ntotal, double nimba, double wtotal, double wimba, int size, int rank);
forw_slint_t SL_PROTO(forw_init_partconds)(forw_slint_t npconds, forw_partcond_t *pconds, forw_slint_t nparts, forw_slint_t total_count, forw_slweight_t total_weight);
forw_slint_t SL_PROTO(forw_init_partconds_intern)(forw_slint_t npconds, forw_partcond_intern_t *pci, forw_partcond_t *pc, forw_slint_t nparts, forw_slint_t total_count, forw_slweight_t total_weight);
forw_slint_t SL_PROTO(forw_merge_partconds)(forw_partcond_t *pconds_in, forw_slint_t npconds_in, forw_partcond_t *pcond_out);
forw_slint_t SL_PROTO(forw_mpi_gather_partconds_grouped)(forw_partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, forw_partcond_t *pconds_out, forw_slint_t *npconds_out, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_gather_partconds)(forw_partcond_t *pcond_in, forw_partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_allgather_partconds)(forw_partcond_t *pcond_in, forw_partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_bcast_partconds)(forw_slint_t npconds, forw_partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_post_check_partconds)(forw_elements_t *s, forw_slint_t nelements, forw_slint_t nparts, forw_partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_post_check_partconds_intern)(forw_elements_t *s, forw_slint_t nelements, forw_slint_t nparts, forw_partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_select_stats)(forw_elements_t *s, forw_slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_select_exact_generic_bulk)(forw_elements_t *s, forw_slint_t nelements, forw_slint_t nparts, forw_partcond_t *pconds, forw_binning_t *bm, forw_splitter_t *sp, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_select_exact_generic_grouped)(forw_elements_t *s, forw_slint_t nelements, forw_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, forw_binning_t *bm, forw_splitter_t *sp, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_select_exact_generic)(forw_elements_t *s, forw_slint_t nelements, forw_slint_t nparts, forw_partcond_t *pconds, forw_binning_t *bm, forw_splitter_t *sp, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_select_exact_radix)(forw_elements_t *s, forw_slint_t nelements, forw_slint_t nparts, forw_partcond_t *pconds, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, forw_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_select_exact_radix_grouped)(forw_elements_t *s, forw_slint_t nelements, forw_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, forw_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_select_sample_regular)(forw_elements_t *s, forw_slint_t nparts, forw_partcond_t *pconds, forw_slint_t nsamples, forw_splitter_t *sp, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_merge)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *xs, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_merge2)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *xs, forw_slint_t merge_type, forw_slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_merge_radix)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *xs, forw_slint_t merge_type, forw_slint_t sort_type, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_partition)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *xs, forw_slint_t part_type, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_partition_radix)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *xs, forw_slint_t part_type, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_partition_exact_radix)(forw_elements_t *s, forw_elements_t *sx, forw_partcond_t *pcond, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_partition_exact_radix_ngroups)(forw_elements_t *s, forw_elements_t *sx, forw_partcond_t *pcond, forw_slint_t ngroups, MPI_Comm *group_comms, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_partition_exact_radix_2groups)(forw_elements_t *s, forw_elements_t *sx, forw_partcond_t *pcond, MPI_Comm group_comm, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_insert_radix)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *xs, forw_slpkey_t *mmkeys, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_presorted_radix)(forw_elements_t *s0, forw_elements_t *s1, forw_elements_t *xs, forw_slint_t merge_type, forw_slint_t rhigh, forw_slint_t rlow, forw_slint_t rwidth, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_sort_back)(forw_elements_t *sin, forw_elements_t *sout, forw_elements_t *sx, forw_slpkey_t *lh, forw_slint_t ntotal, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, forw_slint_t ytotal, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_xcounts2ycounts_grouped)(int *xcounts, forw_slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
forw_slint_t SL_PROTO(forw_mpi_subxdispls2ycounts)(forw_slint_t nsubs, int *sub_xdispls, forw_slint_t *sub_sources, forw_slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#undef SL_PROTO
#endif /* __SL_FORW_H__ */
