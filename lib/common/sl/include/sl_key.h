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


#ifndef __SL_KEY_H__
#define __SL_KEY_H__


/* sl_macro sl_key_type_c sl_key_type_mpi sl_key_size_mpi sl_key_type_fmt sl_key_integer sl_key_memcpy */


#define sl_key_byte                          ((slint_t) sizeof(sl_key_type_c))  /* sl_macro */

#ifndef sl_key_copy
 #ifndef sl_key_memcpy
  #define sl_key_copy(src, dst)              SL_ARRAY1_COPY(src, dst)
 #else
  #define sl_key_copy(src, dst)              memcpy(dst, src, sl_key_byte)  /* sl_macro */
 #endif
#endif
#ifndef sl_key_ncopy
 #define sl_key_ncopy(src, dst, n)           memcpy(dst, src, (n) * sl_key_byte)  /* sl_macro */
#endif
#ifndef sl_key_nmove
 #define sl_key_nmove(src, dst, n)           memmove(dst, src, (n) * sl_key_byte)  /* sl_macro */
#endif


#define key_type_c                           sl_key_type_c  /* sl_macro */
#define key_type_mpi                         (sl_key_type_mpi)  /* sl_macro */
#define key_size_mpi                         (sl_key_size_mpi)  /* sl_macro */
#ifdef sl_key_type_fmt
# define key_type_fmt                        sl_key_type_fmt  /* sl_macro */
#endif
#define key_integer                          sl_key_integer  /* sl_macro */

#define key_pure_type_c                      sl_key_pure_type_c  /* sl_macro */
#define key_pure_type_mpi                    (sl_key_pure_type_mpi)  /* sl_macro */
#define key_pure_size_mpi                    (sl_key_pure_size_mpi)  /* sl_macro */
#ifdef sl_key_pure_type_fmt
# define key_pure_type_fmt                   sl_key_pure_type_fmt  /* sl_macro */
#endif

#define key_purify(k)                        (sl_key_purify(k))  /* sl_macro */
#define key_get_pure(k)                      (sl_key_get_pure(k))  /* sl_macro */
#define key_set_pure(k, p)                   (sl_key_set_pure(k, p))  /* sl_macro */

#ifdef key_integer
# define key_integer_unsigned                (((key_pure_type_c) ~((key_pure_type_c) 0)) >= ((key_pure_type_c) 0))  /* sl_macro */
#endif

#define key_n                                1  /* sl_macro */
#define key_byte                             (sl_key_byte)  /* sl_macro */

#define key_cmp_eq(k0, k1)                   (cc_rti_cadd_cmp(1) sl_key_cmp_eq((k0), (k1)))  /* sl_macro */
#define key_cmp_ne(k0, k1)                   (cc_rti_cadd_cmp(1) sl_key_cmp_ne((k0), (k1)))  /* sl_macro */
#define key_cmp_lt(k0, k1)                   (cc_rti_cadd_cmp(1) sl_key_cmp_lt((k0), (k1)))  /* sl_macro */
#define key_cmp_le(k0, k1)                   (cc_rti_cadd_cmp(1) sl_key_cmp_le((k0), (k1)))  /* sl_macro */
#define key_cmp_gt(k0, k1)                   (cc_rti_cadd_cmp(1) sl_key_cmp_gt((k0), (k1)))  /* sl_macro */
#define key_cmp_ge(k0, k1)                   (cc_rti_cadd_cmp(1) sl_key_cmp_ge((k0), (k1)))  /* sl_macro */

#define key_pure_cmp_eq(k0, k1)              (cc_rti_cadd_cmp(1) sl_key_pure_cmp_eq((k0), (k1)))  /* sl_macro */
#define key_pure_cmp_ne(k0, k1)              (cc_rti_cadd_cmp(1) sl_key_pure_cmp_ne((k0), (k1)))  /* sl_macro */
#define key_pure_cmp_lt(k0, k1)              (cc_rti_cadd_cmp(1) sl_key_pure_cmp_lt((k0), (k1)))  /* sl_macro */
#define key_pure_cmp_le(k0, k1)              (cc_rti_cadd_cmp(1) sl_key_pure_cmp_le((k0), (k1)))  /* sl_macro */
#define key_pure_cmp_gt(k0, k1)              (cc_rti_cadd_cmp(1) sl_key_pure_cmp_gt((k0), (k1)))  /* sl_macro */
#define key_pure_cmp_ge(k0, k1)              (cc_rti_cadd_cmp(1) sl_key_pure_cmp_ge((k0), (k1)))  /* sl_macro */

#ifdef sl_key_val_srand
# define key_val_srand(_s_)                  sl_key_val_srand(_s_)  /* sl_macro */
#endif
#ifdef sl_key_val_rand
# define key_val_rand()                      sl_key_val_rand()  /* sl_macro */
# define have_key_val_rand                   1  /* sl_macro */
#else
# define key_val_rand()                      Z_NOP()
# define have_key_val_rand                   0
#endif
#ifdef sl_key_val_rand_minmax
# define key_val_rand_minmax(_min_, _max_)   sl_key_val_rand_minmax(_min_, _max_)  /* sl_macro */
# define have_key_val_rand_minmax            1  /* sl_macro */
#else
# define key_val_rand_minmax(_min_, _max_)   Z_NOP()
# define have_key_val_rand_minmax            0
#endif


#define key_at(_s_, _sat_)                   ((_s_) + (_sat_))  /* sl_macro */

#define key_assign(src, dst)                 (dst = src)  /* sl_macro */
#define key_assign_at(src, sat, dst)         (dst = &src[sat])  /* sl_macro */
#define key_null(k)                          (k = NULL)  /* sl_macro */
#define key_inc(k)                           (++k)  /* sl_macro */
#define key_dec(k)                           (--k)  /* sl_macro */
#define key_add(k, n)                        (k += n)  /* sl_macro */
#define key_sub(k, n)                        (k -= n)  /* sl_macro */

#define key_copy(src, dst)                   (cc_rti_cadd_movek(1) sl_key_copy(src, dst))  /* sl_macro */
#define key_ncopy(src, dst, n)               (cc_rti_cadd_movek(n) sl_key_ncopy(src, dst, n))  /* sl_macro */
#define key_nmove(src, dst, n)               (cc_rti_cadd_movek(n) sl_key_nmove(src, dst, n))  /* sl_macro */

#define key_copy_at(src, sat, dst, dat)      key_copy(&(src)[sat], &(dst)[dat])  /* sl_macro */
#define key_ncopy_at(src, sat, dst, dat, n)  key_ncopy(&(src)[sat], &(dst)[dat], n)  /* sl_macro */
#define key_nmove_at(src, sat, dst, dat, n)  key_nmove(&(src)[sat], &(dst)[dat], n)  /* sl_macro */

#define key_xchange(k0, k1, t)               (key_copy(k0, t), key_copy(k1, k0), key_copy(t, k1))  /* sl_macro */
#define key_xchange_at(k0, at0, k1, at1, t)  (key_copy_at(k0, at0, t, 0), key_copy_at(k1, at1, k0, at0), key_copy_at(t, 0, k1, at1))  /* sl_macro */

#define key_cm                               SLCM_KEYS  /* sl_macro */

#ifdef key_integer
# define key_radix_low                       ((slint_t) 0)  /* sl_macro */
# define key_radix_high                      ((slint_t) (sizeof(key_pure_type_c) * 8 - 1))  /* sl_macro */
# define key_radix_key2class(_k_, _x_, _y_)  (((_k_) >> (_x_)) & (_y_))  /* sl_macro */
#endif


#endif /* __SL_KEY_H__ */
