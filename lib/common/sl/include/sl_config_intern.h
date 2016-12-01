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


#ifndef __SL_CONFIG_INTERN_H__
#define __SL_CONFIG_INTERN_H__


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


#ifndef SL_INDEX
# undef SL_PACKED_INDEX
#endif


/* if no special datatype for (sl default) integer ... */
#ifndef sl_int_type_c
  /* ... use a default one */
# define sl_int_type_c               long      /* sl_macro */
# undef sl_int_type_mpi
# define sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef sl_int_size_mpi
# define sl_int_size_mpi             1         /* sl_macro */
# undef sl_int_type_fmt
# define sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(sl_int_type_mpi) || !defined(sl_int_size_mpi)
#   error "sl_int_type_mpi and/or sl_int_size_mpi missing"
#  endif
# endif
# ifndef sl_int_type_fmt
#  error "sl_int_type_fmt macro is missing, using d as default"
#  define sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for (intern) weight ... */
#ifndef sl_weight_type_c
 /* ... use (sl default) integer */
# define sl_weight_type_c             sl_int_type_c    /* sl_macro */
# undef sl_weight_type_mpi
# define sl_weight_type_mpi           sl_int_type_mpi  /* sl_macro */
# undef sl_weight_size_mpi
# define sl_weight_size_mpi           sl_int_size_mpi  /* sl_macro */
# undef sl_weight_type_fmt
# define sl_weight_type_fmt           sl_int_type_fmt  /* sl_macro */
# undef sl_weight_intequiv
# define sl_weight_intequiv                            /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(sl_weight_type_mpi) || !defined(sl_weight_size_mpi)
#   error "sl_weight_type_mpi and/or sl_weight_size_mpi missing"
#  endif
# endif
# ifndef sl_weight_type_fmt
#  error "sl_weight_type_fmt macro is missing, using f as default"
#  define sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef sl_index_type_c
 /* ... use the primary integer type */
# define sl_index_type_c             sl_int_type_c
# undef sl_index_type_mpi
# define sl_index_type_mpi           sl_int_type_mpi
# undef sl_index_size_mpi
# define sl_index_size_mpi           sl_int_size_mpi
# undef sl_index_type_fmt
# define sl_index_type_fmt           sl_int_type_fmt
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(sl_index_type_mpi) || !defined(sl_index_size_mpi)
#   error "sl_index_type_mpi and/or sl_index_size_mpi missing"
#  endif
# endif
# ifndef sl_index_type_fmt
#  error "sl_index_type_fmt macro is missing, using d as default"
#  define sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef sl_key_pure_type_c
# define sl_key_pure_type_c          sl_key_type_c  /* sl_macro */
#endif
#ifndef sl_key_pure_type_mpi
# define sl_key_pure_type_mpi        sl_key_type_mpi  /* sl_macro */
#endif
#ifndef sl_key_pure_size_mpi
# define sl_key_pure_size_mpi        sl_key_size_mpi  /* sl_macro */
#endif
#ifndef sl_key_pure_type_fmt
# ifdef sl_key_type_fmt
#  define sl_key_pure_type_fmt       sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef sl_key_purify
 /* key val -> key val */
 #define sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef sl_key_set_pure
 /* key component pointer and key val */
 #define sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef sl_key_pure_cmp_eq
 #define sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_ne
 #define sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_lt
 #define sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_le
 #define sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_gt
 #define sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_ge
 #define sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef sl_key_cmp_eq
 #define sl_key_cmp_eq(k0, k1)       (sl_key_pure_cmp_eq(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_ne
 #define sl_key_cmp_ne(k0, k1)       (sl_key_pure_cmp_ne(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_lt
 #define sl_key_cmp_lt(k0, k1)       (sl_key_pure_cmp_lt(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_le
 #define sl_key_cmp_le(k0, k1)       (sl_key_pure_cmp_le(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_gt
 #define sl_key_cmp_gt(k0, k1)       (sl_key_pure_cmp_gt(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_ge
 #define sl_key_cmp_ge(k0, k1)       (sl_key_pure_cmp_ge(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif


/* default random key */
#ifdef sl_key_integer
# if !defined(sl_key_val_srand) || !defined(sl_key_val_rand) || !defined(sl_key_val_rand_minmax)
#  undef sl_key_val_srand
#  undef sl_key_val_rand
#  undef sl_key_val_rand_minmax
#  define sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define sl_key_val_rand()                     ((sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define sl_key_val_rand_minmax(_min_, _max_)  ((sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef SL_DATA0_IGNORE
# undef SL_DATA0
#endif
#ifdef SL_DATA1_IGNORE
# undef SL_DATA1
#endif
#ifdef SL_DATA2_IGNORE
# undef SL_DATA2
#endif
#ifdef SL_DATA3_IGNORE
# undef SL_DATA3
#endif
#ifdef SL_DATA4_IGNORE
# undef SL_DATA4
#endif
#ifdef SL_DATA5_IGNORE
# undef SL_DATA5
#endif
#ifdef SL_DATA6_IGNORE
# undef SL_DATA6
#endif
#ifdef SL_DATA7_IGNORE
# undef SL_DATA7
#endif
#ifdef SL_DATA8_IGNORE
# undef SL_DATA8
#endif
#ifdef SL_DATA9_IGNORE
# undef SL_DATA9
#endif
#ifdef SL_DATA10_IGNORE
# undef SL_DATA10
#endif
#ifdef SL_DATA11_IGNORE
# undef SL_DATA11
#endif
#ifdef SL_DATA12_IGNORE
# undef SL_DATA12
#endif
#ifdef SL_DATA13_IGNORE
# undef SL_DATA13
#endif
#ifdef SL_DATA14_IGNORE
# undef SL_DATA14
#endif
#ifdef SL_DATA15_IGNORE
# undef SL_DATA15
#endif
#ifdef SL_DATA16_IGNORE
# undef SL_DATA16
#endif
#ifdef SL_DATA17_IGNORE
# undef SL_DATA17
#endif
#ifdef SL_DATA18_IGNORE
# undef SL_DATA18
#endif
#ifdef SL_DATA19_IGNORE
# undef SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef sl_data0_weight
# undef sl_data1_weight
# undef sl_data2_weight
# undef sl_data3_weight
# undef sl_data4_weight
# undef sl_data5_weight
# undef sl_data6_weight
# undef sl_data7_weight
# undef sl_data8_weight
# undef sl_data9_weight
# undef sl_data10_weight
# undef sl_data11_weight
# undef sl_data12_weight
# undef sl_data13_weight
# undef sl_data14_weight
# undef sl_data15_weight
# undef sl_data16_weight
# undef sl_data17_weight
# undef sl_data18_weight
# undef sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(sl_data0_weight) && !defined(SL_DATA0)
# undef sl_elem_weight
#endif
#if defined(sl_data1_weight) && !defined(SL_DATA1)
# undef sl_elem_weight
#endif
#if defined(sl_data2_weight) && !defined(SL_DATA2)
# undef sl_elem_weight
#endif
#if defined(sl_data3_weight) && !defined(SL_DATA3)
# undef sl_elem_weight
#endif
#if defined(sl_data4_weight) && !defined(SL_DATA4)
# undef sl_elem_weight
#endif
#if defined(sl_data5_weight) && !defined(SL_DATA5)
# undef sl_elem_weight
#endif
#if defined(sl_data6_weight) && !defined(SL_DATA6)
# undef sl_elem_weight
#endif
#if defined(sl_data7_weight) && !defined(SL_DATA7)
# undef sl_elem_weight
#endif
#if defined(sl_data8_weight) && !defined(SL_DATA8)
# undef sl_elem_weight
#endif
#if defined(sl_data9_weight) && !defined(SL_DATA9)
# undef sl_elem_weight
#endif
#if defined(sl_data10_weight) && !defined(SL_DATA10)
# undef sl_elem_weight
#endif
#if defined(sl_data11_weight) && !defined(SL_DATA11)
# undef sl_elem_weight
#endif
#if defined(sl_data12_weight) && !defined(SL_DATA12)
# undef sl_elem_weight
#endif
#if defined(sl_data13_weight) && !defined(SL_DATA13)
# undef sl_elem_weight
#endif
#if defined(sl_data14_weight) && !defined(SL_DATA14)
# undef sl_elem_weight
#endif
#if defined(sl_data15_weight) && !defined(SL_DATA15)
# undef sl_elem_weight
#endif
#if defined(sl_data16_weight) && !defined(SL_DATA16)
# undef sl_elem_weight
#endif
#if defined(sl_data17_weight) && !defined(SL_DATA17)
# undef sl_elem_weight
#endif
#if defined(sl_data18_weight) && !defined(SL_DATA18)
# undef sl_elem_weight
#endif
#if defined(sl_data19_weight) && !defined(SL_DATA19)
# undef sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef SL_DATA0
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data0_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA1
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data1_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA2
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data2_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA3
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data3_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA4
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data4_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA5
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data5_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA6
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data6_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA7
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data7_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA8
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data8_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA9
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data9_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA10
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data10_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA11
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data11_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA12
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data12_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA13
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data13_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA14
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data14_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA15
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data15_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA16
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data16_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA17
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data17_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA18
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data18_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
#ifdef SL_DATA19
# ifdef FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef sl_data19_flex
#   define FLECKS_GUARD
#  endif
# endif
#endif
/* DATAX_TEMPLATE_END */


#endif /* __SL_CONFIG_INTERN_H__ */
