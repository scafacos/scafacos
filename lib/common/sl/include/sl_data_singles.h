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


#ifndef __SL_DATA_SINGLES_H__
#define __SL_DATA_SINGLES_H__


/* DATAX_TEMPLATE_BEGIN */

/* sl_macro SL_DATA0 SL_DATA0_IGNORE sl_data0_type_c sl_data0_size_c sl_data0_type_mpi sl_data0_size_mpi sl_data0_memcpy sl_data0_weight sl_data0_flex */

#ifdef SL_DATA0

 #define sl_data0_byte                             ((slint_t) (sl_data0_size_c) * sizeof(sl_data0_type_c))  /* sl_macro */

 #ifndef sl_data0_copy
  #if sl_data0_size_c <= 9 && !defined(sl_data0_memcpy)
   #if sl_data0_size_c == 1
    #define sl_data0_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data0_size_c == 2
    #define sl_data0_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data0_size_c == 3
    #define sl_data0_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data0_size_c == 4
    #define sl_data0_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data0_size_c == 5
    #define sl_data0_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data0_size_c == 6
    #define sl_data0_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data0_size_c == 7
    #define sl_data0_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data0_size_c == 8
    #define sl_data0_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data0_size_c == 9
    #define sl_data0_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data0_copy(src, dst)                 memcpy(dst, src, sl_data0_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data0_ncopy
  #define sl_data0_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data0_byte)  /* sl_macro */
 #endif
 #ifndef sl_data0_nmove
  #define sl_data0_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data0_byte)  /* sl_macro */
 #endif

 #define data0_type_c                              sl_data0_type_c  /* sl_macro */
 #define data0_size_c                              (sl_data0_size_c)  /* sl_macro */
 #define data0_type_mpi                            (sl_data0_type_mpi)  /* sl_macro */
 #define data0_size_mpi                            (sl_data0_size_mpi)  /* sl_macro */

 #define data0_idx                                 0  /* sl_macro */

 #define data0_n                                   1  /* sl_macro */
 #define data0_byte                                (sl_data0_byte)  /* sl_macro */
 #define data0_ptr(e)                              (e)->data0  /* sl_macro */

 #ifdef sl_data0_flex
 # define data0_byte_flex                          (sl_data0_byte)  /* sl_macro */
 #else
 # define data0_byte_flex                          0
 #endif

 #ifdef sl_data0_weight
 # define data0_weight                             1  /* sl_macro */
 #else
 # define data0_weight                             0
 #endif

 /* commands for regular use */
 #define data0_assign(src, dst)                    ((dst)->data0 = (src)->data0)  /* sl_macro */
 #define data0_assign_at(src, sat, dst)            ((dst)->data0 = &(src)->data0[(sat) * data0_size_c])  /* sl_macro */
 #define data0_null(e)                             ((e)->data0 = NULL)  /* sl_macro */
 #define data0_inc(e)                              ((e)->data0 += data0_size_c)  /* sl_macro */
 #define data0_dec(e)                              ((e)->data0 -= data0_size_c)  /* sl_macro */
 #define data0_add(e, n)                           ((e)->data0 += (n) * data0_size_c)  /* sl_macro */
 #define data0_sub(e, n)                           ((e)->data0 -= (n) * data0_size_c)  /* sl_macro */

 #define data0_copy(src, dst)                      sl_data0_copy((src)->data0, (dst)->data0)  /* sl_macro */
 #define data0_ncopy(src, dst, n)                  sl_data0_ncopy((src)->data0, (dst)->data0, n)  /* sl_macro */
 #define data0_nmove(src, dst, n)                  sl_data0_nmove((src)->data0, (dst)->data0, n)  /* sl_macro */

 #define data0_copy_at(src, sat, dst, dat)         sl_data0_copy(&(src)->data0[(sat) * data0_size_c], &(dst)->data0[(dat) * data0_size_c])  /* sl_macro */
 #define data0_ncopy_at(src, sat, dst, dat, n)     sl_data0_ncopy(&(src)->data0[(sat) * data0_size_c], &(dst)->data0[(dat) * data0_size_c], n)  /* sl_macro */
 #define data0_nmove_at(src, sat, dst, dat, n)     sl_data0_nmove(&(src)->data0[(sat) * data0_size_c], &(dst)->data0[(dat) * data0_size_c], n)  /* sl_macro */

 #define data0_xchange(e0, e1, t)                  (data0_copy(e0, t), data0_copy(e1, e0), data0_copy(t, e1))  /* sl_macro */
 #define data0_xchange_at(e0, at0, e1, at1, t)     (data0_copy_at(e0, at0, t, 0), data0_copy_at(e1, at1, e0, at0), data0_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data0_assign(src, dst)                 , data0_assign(src, dst)  /* sl_macro */
 #define cc_data0_assign_at(src, sat, dst)         , data0_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data0_null(e)                          , data0_null(e)  /* sl_macro */
 #define cc_data0_inc(e)                           , data0_inc(e)  /* sl_macro */
 #define cc_data0_dec(e)                           , data0_dec(e)  /* sl_macro */
 #define cc_data0_add(e, n)                        , data0_add(e, n)  /* sl_macro */
 #define cc_data0_sub(e, n)                        , data0_sub(e, n)  /* sl_macro */
 #define cc_data0_copy(src, dst)                   , data0_copy(src, dst)  /* sl_macro */
 #define cc_data0_ncopy(src, dst, n)               , data0_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data0_nmove(src, dst, n)               , data0_nmove(src, dst, n)  /* sl_macro */
 #define cc_data0_copy_at(src, sat, dst, dat)      , data0_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data0_ncopy_at(src, sat, dst, dat, n)  , data0_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data0_nmove_at(src, sat, dst, dat, n)  , data0_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data0_xchange(e0, e1, t)               , data0_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data0_xchange_at(e0, at0, e1, at1, t)  , data0_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data0_assign(src, dst)                 data0_assign(src, dst)
 #define cc_data0_assign_at(src, sat, dst)         data0_assign_at(src, sat, dst)
 #define cc_data0_null(e)                          data0_null(e)
 #define cc_data0_inc(e)                           data0_inc(e)
 #define cc_data0_dec(e)                           data0_dec(e)
 #define cc_data0_add(e, n)                        data0_add(e, n)
 #define cc_data0_sub(e, n)                        data0_sub(e, n)
 #define cc_data0_copy(src, dst)                   data0_copy(src, dst)
 #define cc_data0_ncopy(src, dst, n)               data0_ncopy(src, dst, n)
 #define cc_data0_nmove(src, dst, n)               data0_nmove(src, dst, n)
 #define cc_data0_copy_at(src, sat, dst, dat)      data0_copy_at(src, sat, dst, dat)
 #define cc_data0_ncopy_at(src, sat, dst, dat, n)  data0_ncopy_at(src, sat, dst, dat, n)
 #define cc_data0_nmove_at(src, sat, dst, dat, n)  data0_nmove_at(src, sat, dst, dat, n)
 #define cc_data0_xchange(e0, e1, t)               data0_xchange(e0, e1, t)
 #define cc_data0_xchange_at(e0, at0, e1, at1, t)  data0_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA0 */

 #define data0_n                                   0
 #define data0_byte                                0
/* #define data0_ptr(e)*/

 #define data0_byte_flex                           0
 #define data0_weight                              0

 /* commands for regular use */
 #define data0_assign(src, dst)                    Z_NOP()
 #define data0_assign_at(src, sat, dst)            Z_NOP()
 #define data0_null(e)                             Z_NOP()
 #define data0_inc(e)                              Z_NOP()
 #define data0_dec(e)                              Z_NOP()
 #define data0_add(e, n)                           Z_NOP()
 #define data0_sub(e, n)                           Z_NOP()
 #define data0_copy(src, dst)                      Z_NOP()
 #define data0_ncopy(src, dst, n)                  Z_NOP()
 #define data0_nmove(src, dst, n)                  Z_NOP()
 #define data0_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data0_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data0_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data0_xchange(e0, e1, t)                  Z_NOP()
 #define data0_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data0_assign(src, dst)
 #define cc_data0_assign_at(src, sat, dst)
 #define cc_data0_null(e)
 #define cc_data0_inc(e)
 #define cc_data0_dec(e)
 #define cc_data0_add(e, n)
 #define cc_data0_sub(e, n)
 #define cc_data0_copy(src, dst)
 #define cc_data0_ncopy(src, dst, n)
 #define cc_data0_nmove(src, dst, n)
 #define cc_data0_copy_at(src, sat, dst, dat)
 #define cc_data0_ncopy_at(src, sat, dst, dat, n)
 #define cc_data0_nmove_at(src, sat, dst, dat, n)
 #define cc_data0_xchange(e0, e1, t)
 #define cc_data0_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA0 */

#define data0_cm                                   SLCM_DATA0  /* sl_macro */


/* sl_macro SL_DATA1 SL_DATA1_IGNORE sl_data1_type_c sl_data1_size_c sl_data1_type_mpi sl_data1_size_mpi sl_data1_memcpy sl_data1_weight sl_data1_flex */

#ifdef SL_DATA1

 #define sl_data1_byte                             ((slint_t) (sl_data1_size_c) * sizeof(sl_data1_type_c))  /* sl_macro */

 #ifndef sl_data1_copy
  #if sl_data1_size_c <= 9 && !defined(sl_data1_memcpy)
   #if sl_data1_size_c == 1
    #define sl_data1_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data1_size_c == 2
    #define sl_data1_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data1_size_c == 3
    #define sl_data1_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data1_size_c == 4
    #define sl_data1_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data1_size_c == 5
    #define sl_data1_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data1_size_c == 6
    #define sl_data1_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data1_size_c == 7
    #define sl_data1_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data1_size_c == 8
    #define sl_data1_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data1_size_c == 9
    #define sl_data1_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data1_copy(src, dst)                 memcpy(dst, src, sl_data1_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data1_ncopy
  #define sl_data1_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data1_byte)  /* sl_macro */
 #endif
 #ifndef sl_data1_nmove
  #define sl_data1_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data1_byte)  /* sl_macro */
 #endif

 #define data1_type_c                              sl_data1_type_c  /* sl_macro */
 #define data1_size_c                              (sl_data1_size_c)  /* sl_macro */
 #define data1_type_mpi                            (sl_data1_type_mpi)  /* sl_macro */
 #define data1_size_mpi                            (sl_data1_size_mpi)  /* sl_macro */

 #define data1_idx                                 1  /* sl_macro */

 #define data1_n                                   1  /* sl_macro */
 #define data1_byte                                (sl_data1_byte)  /* sl_macro */
 #define data1_ptr(e)                              (e)->data1  /* sl_macro */

 #ifdef sl_data1_flex
 # define data1_byte_flex                          (sl_data1_byte)  /* sl_macro */
 #else
 # define data1_byte_flex                          0
 #endif

 #ifdef sl_data1_weight
 # define data1_weight                             1  /* sl_macro */
 #else
 # define data1_weight                             0
 #endif

 /* commands for regular use */
 #define data1_assign(src, dst)                    ((dst)->data1 = (src)->data1)  /* sl_macro */
 #define data1_assign_at(src, sat, dst)            ((dst)->data1 = &(src)->data1[(sat) * data1_size_c])  /* sl_macro */
 #define data1_null(e)                             ((e)->data1 = NULL)  /* sl_macro */
 #define data1_inc(e)                              ((e)->data1 += data1_size_c)  /* sl_macro */
 #define data1_dec(e)                              ((e)->data1 -= data1_size_c)  /* sl_macro */
 #define data1_add(e, n)                           ((e)->data1 += (n) * data1_size_c)  /* sl_macro */
 #define data1_sub(e, n)                           ((e)->data1 -= (n) * data1_size_c)  /* sl_macro */

 #define data1_copy(src, dst)                      sl_data1_copy((src)->data1, (dst)->data1)  /* sl_macro */
 #define data1_ncopy(src, dst, n)                  sl_data1_ncopy((src)->data1, (dst)->data1, n)  /* sl_macro */
 #define data1_nmove(src, dst, n)                  sl_data1_nmove((src)->data1, (dst)->data1, n)  /* sl_macro */

 #define data1_copy_at(src, sat, dst, dat)         sl_data1_copy(&(src)->data1[(sat) * data1_size_c], &(dst)->data1[(dat) * data1_size_c])  /* sl_macro */
 #define data1_ncopy_at(src, sat, dst, dat, n)     sl_data1_ncopy(&(src)->data1[(sat) * data1_size_c], &(dst)->data1[(dat) * data1_size_c], n)  /* sl_macro */
 #define data1_nmove_at(src, sat, dst, dat, n)     sl_data1_nmove(&(src)->data1[(sat) * data1_size_c], &(dst)->data1[(dat) * data1_size_c], n)  /* sl_macro */

 #define data1_xchange(e0, e1, t)                  (data1_copy(e0, t), data1_copy(e1, e0), data1_copy(t, e1))  /* sl_macro */
 #define data1_xchange_at(e0, at0, e1, at1, t)     (data1_copy_at(e0, at0, t, 0), data1_copy_at(e1, at1, e0, at0), data1_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data1_assign(src, dst)                 , data1_assign(src, dst)  /* sl_macro */
 #define cc_data1_assign_at(src, sat, dst)         , data1_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data1_null(e)                          , data1_null(e)  /* sl_macro */
 #define cc_data1_inc(e)                           , data1_inc(e)  /* sl_macro */
 #define cc_data1_dec(e)                           , data1_dec(e)  /* sl_macro */
 #define cc_data1_add(e, n)                        , data1_add(e, n)  /* sl_macro */
 #define cc_data1_sub(e, n)                        , data1_sub(e, n)  /* sl_macro */
 #define cc_data1_copy(src, dst)                   , data1_copy(src, dst)  /* sl_macro */
 #define cc_data1_ncopy(src, dst, n)               , data1_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data1_nmove(src, dst, n)               , data1_nmove(src, dst, n)  /* sl_macro */
 #define cc_data1_copy_at(src, sat, dst, dat)      , data1_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data1_ncopy_at(src, sat, dst, dat, n)  , data1_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data1_nmove_at(src, sat, dst, dat, n)  , data1_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data1_xchange(e0, e1, t)               , data1_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data1_xchange_at(e0, at0, e1, at1, t)  , data1_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data1_assign(src, dst)                 data1_assign(src, dst)
 #define cc_data1_assign_at(src, sat, dst)         data1_assign_at(src, sat, dst)
 #define cc_data1_null(e)                          data1_null(e)
 #define cc_data1_inc(e)                           data1_inc(e)
 #define cc_data1_dec(e)                           data1_dec(e)
 #define cc_data1_add(e, n)                        data1_add(e, n)
 #define cc_data1_sub(e, n)                        data1_sub(e, n)
 #define cc_data1_copy(src, dst)                   data1_copy(src, dst)
 #define cc_data1_ncopy(src, dst, n)               data1_ncopy(src, dst, n)
 #define cc_data1_nmove(src, dst, n)               data1_nmove(src, dst, n)
 #define cc_data1_copy_at(src, sat, dst, dat)      data1_copy_at(src, sat, dst, dat)
 #define cc_data1_ncopy_at(src, sat, dst, dat, n)  data1_ncopy_at(src, sat, dst, dat, n)
 #define cc_data1_nmove_at(src, sat, dst, dat, n)  data1_nmove_at(src, sat, dst, dat, n)
 #define cc_data1_xchange(e0, e1, t)               data1_xchange(e0, e1, t)
 #define cc_data1_xchange_at(e0, at0, e1, at1, t)  data1_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA1 */

 #define data1_n                                   0
 #define data1_byte                                0
/* #define data1_ptr(e)*/

 #define data1_byte_flex                           0
 #define data1_weight                              0

 /* commands for regular use */
 #define data1_assign(src, dst)                    Z_NOP()
 #define data1_assign_at(src, sat, dst)            Z_NOP()
 #define data1_null(e)                             Z_NOP()
 #define data1_inc(e)                              Z_NOP()
 #define data1_dec(e)                              Z_NOP()
 #define data1_add(e, n)                           Z_NOP()
 #define data1_sub(e, n)                           Z_NOP()
 #define data1_copy(src, dst)                      Z_NOP()
 #define data1_ncopy(src, dst, n)                  Z_NOP()
 #define data1_nmove(src, dst, n)                  Z_NOP()
 #define data1_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data1_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data1_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data1_xchange(e0, e1, t)                  Z_NOP()
 #define data1_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data1_assign(src, dst)
 #define cc_data1_assign_at(src, sat, dst)
 #define cc_data1_null(e)
 #define cc_data1_inc(e)
 #define cc_data1_dec(e)
 #define cc_data1_add(e, n)
 #define cc_data1_sub(e, n)
 #define cc_data1_copy(src, dst)
 #define cc_data1_ncopy(src, dst, n)
 #define cc_data1_nmove(src, dst, n)
 #define cc_data1_copy_at(src, sat, dst, dat)
 #define cc_data1_ncopy_at(src, sat, dst, dat, n)
 #define cc_data1_nmove_at(src, sat, dst, dat, n)
 #define cc_data1_xchange(e0, e1, t)
 #define cc_data1_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA1 */

#define data1_cm                                   SLCM_DATA1  /* sl_macro */


/* sl_macro SL_DATA2 SL_DATA2_IGNORE sl_data2_type_c sl_data2_size_c sl_data2_type_mpi sl_data2_size_mpi sl_data2_memcpy sl_data2_weight sl_data2_flex */

#ifdef SL_DATA2

 #define sl_data2_byte                             ((slint_t) (sl_data2_size_c) * sizeof(sl_data2_type_c))  /* sl_macro */

 #ifndef sl_data2_copy
  #if sl_data2_size_c <= 9 && !defined(sl_data2_memcpy)
   #if sl_data2_size_c == 1
    #define sl_data2_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data2_size_c == 2
    #define sl_data2_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data2_size_c == 3
    #define sl_data2_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data2_size_c == 4
    #define sl_data2_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data2_size_c == 5
    #define sl_data2_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data2_size_c == 6
    #define sl_data2_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data2_size_c == 7
    #define sl_data2_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data2_size_c == 8
    #define sl_data2_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data2_size_c == 9
    #define sl_data2_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data2_copy(src, dst)                 memcpy(dst, src, sl_data2_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data2_ncopy
  #define sl_data2_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data2_byte)  /* sl_macro */
 #endif
 #ifndef sl_data2_nmove
  #define sl_data2_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data2_byte)  /* sl_macro */
 #endif

 #define data2_type_c                              sl_data2_type_c  /* sl_macro */
 #define data2_size_c                              (sl_data2_size_c)  /* sl_macro */
 #define data2_type_mpi                            (sl_data2_type_mpi)  /* sl_macro */
 #define data2_size_mpi                            (sl_data2_size_mpi)  /* sl_macro */

 #define data2_idx                                 2  /* sl_macro */

 #define data2_n                                   1  /* sl_macro */
 #define data2_byte                                (sl_data2_byte)  /* sl_macro */
 #define data2_ptr(e)                              (e)->data2  /* sl_macro */

 #ifdef sl_data2_flex
 # define data2_byte_flex                          (sl_data2_byte)  /* sl_macro */
 #else
 # define data2_byte_flex                          0
 #endif

 #ifdef sl_data2_weight
 # define data2_weight                             1  /* sl_macro */
 #else
 # define data2_weight                             0
 #endif

 /* commands for regular use */
 #define data2_assign(src, dst)                    ((dst)->data2 = (src)->data2)  /* sl_macro */
 #define data2_assign_at(src, sat, dst)            ((dst)->data2 = &(src)->data2[(sat) * data2_size_c])  /* sl_macro */
 #define data2_null(e)                             ((e)->data2 = NULL)  /* sl_macro */
 #define data2_inc(e)                              ((e)->data2 += data2_size_c)  /* sl_macro */
 #define data2_dec(e)                              ((e)->data2 -= data2_size_c)  /* sl_macro */
 #define data2_add(e, n)                           ((e)->data2 += (n) * data2_size_c)  /* sl_macro */
 #define data2_sub(e, n)                           ((e)->data2 -= (n) * data2_size_c)  /* sl_macro */

 #define data2_copy(src, dst)                      sl_data2_copy((src)->data2, (dst)->data2)  /* sl_macro */
 #define data2_ncopy(src, dst, n)                  sl_data2_ncopy((src)->data2, (dst)->data2, n)  /* sl_macro */
 #define data2_nmove(src, dst, n)                  sl_data2_nmove((src)->data2, (dst)->data2, n)  /* sl_macro */

 #define data2_copy_at(src, sat, dst, dat)         sl_data2_copy(&(src)->data2[(sat) * data2_size_c], &(dst)->data2[(dat) * data2_size_c])  /* sl_macro */
 #define data2_ncopy_at(src, sat, dst, dat, n)     sl_data2_ncopy(&(src)->data2[(sat) * data2_size_c], &(dst)->data2[(dat) * data2_size_c], n)  /* sl_macro */
 #define data2_nmove_at(src, sat, dst, dat, n)     sl_data2_nmove(&(src)->data2[(sat) * data2_size_c], &(dst)->data2[(dat) * data2_size_c], n)  /* sl_macro */

 #define data2_xchange(e0, e1, t)                  (data2_copy(e0, t), data2_copy(e1, e0), data2_copy(t, e1))  /* sl_macro */
 #define data2_xchange_at(e0, at0, e1, at1, t)     (data2_copy_at(e0, at0, t, 0), data2_copy_at(e1, at1, e0, at0), data2_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data2_assign(src, dst)                 , data2_assign(src, dst)  /* sl_macro */
 #define cc_data2_assign_at(src, sat, dst)         , data2_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data2_null(e)                          , data2_null(e)  /* sl_macro */
 #define cc_data2_inc(e)                           , data2_inc(e)  /* sl_macro */
 #define cc_data2_dec(e)                           , data2_dec(e)  /* sl_macro */
 #define cc_data2_add(e, n)                        , data2_add(e, n)  /* sl_macro */
 #define cc_data2_sub(e, n)                        , data2_sub(e, n)  /* sl_macro */
 #define cc_data2_copy(src, dst)                   , data2_copy(src, dst)  /* sl_macro */
 #define cc_data2_ncopy(src, dst, n)               , data2_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data2_nmove(src, dst, n)               , data2_nmove(src, dst, n)  /* sl_macro */
 #define cc_data2_copy_at(src, sat, dst, dat)      , data2_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data2_ncopy_at(src, sat, dst, dat, n)  , data2_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data2_nmove_at(src, sat, dst, dat, n)  , data2_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data2_xchange(e0, e1, t)               , data2_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data2_xchange_at(e0, at0, e1, at1, t)  , data2_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data2_assign(src, dst)                 data2_assign(src, dst)
 #define cc_data2_assign_at(src, sat, dst)         data2_assign_at(src, sat, dst)
 #define cc_data2_null(e)                          data2_null(e)
 #define cc_data2_inc(e)                           data2_inc(e)
 #define cc_data2_dec(e)                           data2_dec(e)
 #define cc_data2_add(e, n)                        data2_add(e, n)
 #define cc_data2_sub(e, n)                        data2_sub(e, n)
 #define cc_data2_copy(src, dst)                   data2_copy(src, dst)
 #define cc_data2_ncopy(src, dst, n)               data2_ncopy(src, dst, n)
 #define cc_data2_nmove(src, dst, n)               data2_nmove(src, dst, n)
 #define cc_data2_copy_at(src, sat, dst, dat)      data2_copy_at(src, sat, dst, dat)
 #define cc_data2_ncopy_at(src, sat, dst, dat, n)  data2_ncopy_at(src, sat, dst, dat, n)
 #define cc_data2_nmove_at(src, sat, dst, dat, n)  data2_nmove_at(src, sat, dst, dat, n)
 #define cc_data2_xchange(e0, e1, t)               data2_xchange(e0, e1, t)
 #define cc_data2_xchange_at(e0, at0, e1, at1, t)  data2_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA2 */

 #define data2_n                                   0
 #define data2_byte                                0
/* #define data2_ptr(e)*/

 #define data2_byte_flex                           0
 #define data2_weight                              0

 /* commands for regular use */
 #define data2_assign(src, dst)                    Z_NOP()
 #define data2_assign_at(src, sat, dst)            Z_NOP()
 #define data2_null(e)                             Z_NOP()
 #define data2_inc(e)                              Z_NOP()
 #define data2_dec(e)                              Z_NOP()
 #define data2_add(e, n)                           Z_NOP()
 #define data2_sub(e, n)                           Z_NOP()
 #define data2_copy(src, dst)                      Z_NOP()
 #define data2_ncopy(src, dst, n)                  Z_NOP()
 #define data2_nmove(src, dst, n)                  Z_NOP()
 #define data2_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data2_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data2_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data2_xchange(e0, e1, t)                  Z_NOP()
 #define data2_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data2_assign(src, dst)
 #define cc_data2_assign_at(src, sat, dst)
 #define cc_data2_null(e)
 #define cc_data2_inc(e)
 #define cc_data2_dec(e)
 #define cc_data2_add(e, n)
 #define cc_data2_sub(e, n)
 #define cc_data2_copy(src, dst)
 #define cc_data2_ncopy(src, dst, n)
 #define cc_data2_nmove(src, dst, n)
 #define cc_data2_copy_at(src, sat, dst, dat)
 #define cc_data2_ncopy_at(src, sat, dst, dat, n)
 #define cc_data2_nmove_at(src, sat, dst, dat, n)
 #define cc_data2_xchange(e0, e1, t)
 #define cc_data2_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA2 */

#define data2_cm                                   SLCM_DATA2  /* sl_macro */


/* sl_macro SL_DATA3 SL_DATA3_IGNORE sl_data3_type_c sl_data3_size_c sl_data3_type_mpi sl_data3_size_mpi sl_data3_memcpy sl_data3_weight sl_data3_flex */

#ifdef SL_DATA3

 #define sl_data3_byte                             ((slint_t) (sl_data3_size_c) * sizeof(sl_data3_type_c))  /* sl_macro */

 #ifndef sl_data3_copy
  #if sl_data3_size_c <= 9 && !defined(sl_data3_memcpy)
   #if sl_data3_size_c == 1
    #define sl_data3_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data3_size_c == 2
    #define sl_data3_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data3_size_c == 3
    #define sl_data3_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data3_size_c == 4
    #define sl_data3_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data3_size_c == 5
    #define sl_data3_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data3_size_c == 6
    #define sl_data3_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data3_size_c == 7
    #define sl_data3_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data3_size_c == 8
    #define sl_data3_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data3_size_c == 9
    #define sl_data3_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data3_copy(src, dst)                 memcpy(dst, src, sl_data3_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data3_ncopy
  #define sl_data3_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data3_byte)  /* sl_macro */
 #endif
 #ifndef sl_data3_nmove
  #define sl_data3_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data3_byte)  /* sl_macro */
 #endif

 #define data3_type_c                              sl_data3_type_c  /* sl_macro */
 #define data3_size_c                              (sl_data3_size_c)  /* sl_macro */
 #define data3_type_mpi                            (sl_data3_type_mpi)  /* sl_macro */
 #define data3_size_mpi                            (sl_data3_size_mpi)  /* sl_macro */

 #define data3_idx                                 3  /* sl_macro */

 #define data3_n                                   1  /* sl_macro */
 #define data3_byte                                (sl_data3_byte)  /* sl_macro */
 #define data3_ptr(e)                              (e)->data3  /* sl_macro */

 #ifdef sl_data3_flex
 # define data3_byte_flex                          (sl_data3_byte)  /* sl_macro */
 #else
 # define data3_byte_flex                          0
 #endif

 #ifdef sl_data3_weight
 # define data3_weight                             1  /* sl_macro */
 #else
 # define data3_weight                             0
 #endif

 /* commands for regular use */
 #define data3_assign(src, dst)                    ((dst)->data3 = (src)->data3)  /* sl_macro */
 #define data3_assign_at(src, sat, dst)            ((dst)->data3 = &(src)->data3[(sat) * data3_size_c])  /* sl_macro */
 #define data3_null(e)                             ((e)->data3 = NULL)  /* sl_macro */
 #define data3_inc(e)                              ((e)->data3 += data3_size_c)  /* sl_macro */
 #define data3_dec(e)                              ((e)->data3 -= data3_size_c)  /* sl_macro */
 #define data3_add(e, n)                           ((e)->data3 += (n) * data3_size_c)  /* sl_macro */
 #define data3_sub(e, n)                           ((e)->data3 -= (n) * data3_size_c)  /* sl_macro */

 #define data3_copy(src, dst)                      sl_data3_copy((src)->data3, (dst)->data3)  /* sl_macro */
 #define data3_ncopy(src, dst, n)                  sl_data3_ncopy((src)->data3, (dst)->data3, n)  /* sl_macro */
 #define data3_nmove(src, dst, n)                  sl_data3_nmove((src)->data3, (dst)->data3, n)  /* sl_macro */

 #define data3_copy_at(src, sat, dst, dat)         sl_data3_copy(&(src)->data3[(sat) * data3_size_c], &(dst)->data3[(dat) * data3_size_c])  /* sl_macro */
 #define data3_ncopy_at(src, sat, dst, dat, n)     sl_data3_ncopy(&(src)->data3[(sat) * data3_size_c], &(dst)->data3[(dat) * data3_size_c], n)  /* sl_macro */
 #define data3_nmove_at(src, sat, dst, dat, n)     sl_data3_nmove(&(src)->data3[(sat) * data3_size_c], &(dst)->data3[(dat) * data3_size_c], n)  /* sl_macro */

 #define data3_xchange(e0, e1, t)                  (data3_copy(e0, t), data3_copy(e1, e0), data3_copy(t, e1))  /* sl_macro */
 #define data3_xchange_at(e0, at0, e1, at1, t)     (data3_copy_at(e0, at0, t, 0), data3_copy_at(e1, at1, e0, at0), data3_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data3_assign(src, dst)                 , data3_assign(src, dst)  /* sl_macro */
 #define cc_data3_assign_at(src, sat, dst)         , data3_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data3_null(e)                          , data3_null(e)  /* sl_macro */
 #define cc_data3_inc(e)                           , data3_inc(e)  /* sl_macro */
 #define cc_data3_dec(e)                           , data3_dec(e)  /* sl_macro */
 #define cc_data3_add(e, n)                        , data3_add(e, n)  /* sl_macro */
 #define cc_data3_sub(e, n)                        , data3_sub(e, n)  /* sl_macro */
 #define cc_data3_copy(src, dst)                   , data3_copy(src, dst)  /* sl_macro */
 #define cc_data3_ncopy(src, dst, n)               , data3_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data3_nmove(src, dst, n)               , data3_nmove(src, dst, n)  /* sl_macro */
 #define cc_data3_copy_at(src, sat, dst, dat)      , data3_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data3_ncopy_at(src, sat, dst, dat, n)  , data3_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data3_nmove_at(src, sat, dst, dat, n)  , data3_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data3_xchange(e0, e1, t)               , data3_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data3_xchange_at(e0, at0, e1, at1, t)  , data3_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data3_assign(src, dst)                 data3_assign(src, dst)
 #define cc_data3_assign_at(src, sat, dst)         data3_assign_at(src, sat, dst)
 #define cc_data3_null(e)                          data3_null(e)
 #define cc_data3_inc(e)                           data3_inc(e)
 #define cc_data3_dec(e)                           data3_dec(e)
 #define cc_data3_add(e, n)                        data3_add(e, n)
 #define cc_data3_sub(e, n)                        data3_sub(e, n)
 #define cc_data3_copy(src, dst)                   data3_copy(src, dst)
 #define cc_data3_ncopy(src, dst, n)               data3_ncopy(src, dst, n)
 #define cc_data3_nmove(src, dst, n)               data3_nmove(src, dst, n)
 #define cc_data3_copy_at(src, sat, dst, dat)      data3_copy_at(src, sat, dst, dat)
 #define cc_data3_ncopy_at(src, sat, dst, dat, n)  data3_ncopy_at(src, sat, dst, dat, n)
 #define cc_data3_nmove_at(src, sat, dst, dat, n)  data3_nmove_at(src, sat, dst, dat, n)
 #define cc_data3_xchange(e0, e1, t)               data3_xchange(e0, e1, t)
 #define cc_data3_xchange_at(e0, at0, e1, at1, t)  data3_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA3 */

 #define data3_n                                   0
 #define data3_byte                                0
/* #define data3_ptr(e)*/

 #define data3_byte_flex                           0
 #define data3_weight                              0

 /* commands for regular use */
 #define data3_assign(src, dst)                    Z_NOP()
 #define data3_assign_at(src, sat, dst)            Z_NOP()
 #define data3_null(e)                             Z_NOP()
 #define data3_inc(e)                              Z_NOP()
 #define data3_dec(e)                              Z_NOP()
 #define data3_add(e, n)                           Z_NOP()
 #define data3_sub(e, n)                           Z_NOP()
 #define data3_copy(src, dst)                      Z_NOP()
 #define data3_ncopy(src, dst, n)                  Z_NOP()
 #define data3_nmove(src, dst, n)                  Z_NOP()
 #define data3_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data3_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data3_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data3_xchange(e0, e1, t)                  Z_NOP()
 #define data3_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data3_assign(src, dst)
 #define cc_data3_assign_at(src, sat, dst)
 #define cc_data3_null(e)
 #define cc_data3_inc(e)
 #define cc_data3_dec(e)
 #define cc_data3_add(e, n)
 #define cc_data3_sub(e, n)
 #define cc_data3_copy(src, dst)
 #define cc_data3_ncopy(src, dst, n)
 #define cc_data3_nmove(src, dst, n)
 #define cc_data3_copy_at(src, sat, dst, dat)
 #define cc_data3_ncopy_at(src, sat, dst, dat, n)
 #define cc_data3_nmove_at(src, sat, dst, dat, n)
 #define cc_data3_xchange(e0, e1, t)
 #define cc_data3_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA3 */

#define data3_cm                                   SLCM_DATA3  /* sl_macro */


/* sl_macro SL_DATA4 SL_DATA4_IGNORE sl_data4_type_c sl_data4_size_c sl_data4_type_mpi sl_data4_size_mpi sl_data4_memcpy sl_data4_weight sl_data4_flex */

#ifdef SL_DATA4

 #define sl_data4_byte                             ((slint_t) (sl_data4_size_c) * sizeof(sl_data4_type_c))  /* sl_macro */

 #ifndef sl_data4_copy
  #if sl_data4_size_c <= 9 && !defined(sl_data4_memcpy)
   #if sl_data4_size_c == 1
    #define sl_data4_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data4_size_c == 2
    #define sl_data4_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data4_size_c == 3
    #define sl_data4_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data4_size_c == 4
    #define sl_data4_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data4_size_c == 5
    #define sl_data4_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data4_size_c == 6
    #define sl_data4_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data4_size_c == 7
    #define sl_data4_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data4_size_c == 8
    #define sl_data4_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data4_size_c == 9
    #define sl_data4_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data4_copy(src, dst)                 memcpy(dst, src, sl_data4_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data4_ncopy
  #define sl_data4_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data4_byte)  /* sl_macro */
 #endif
 #ifndef sl_data4_nmove
  #define sl_data4_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data4_byte)  /* sl_macro */
 #endif

 #define data4_type_c                              sl_data4_type_c  /* sl_macro */
 #define data4_size_c                              (sl_data4_size_c)  /* sl_macro */
 #define data4_type_mpi                            (sl_data4_type_mpi)  /* sl_macro */
 #define data4_size_mpi                            (sl_data4_size_mpi)  /* sl_macro */

 #define data4_idx                                 4  /* sl_macro */

 #define data4_n                                   1  /* sl_macro */
 #define data4_byte                                (sl_data4_byte)  /* sl_macro */
 #define data4_ptr(e)                              (e)->data4  /* sl_macro */

 #ifdef sl_data4_flex
 # define data4_byte_flex                          (sl_data4_byte)  /* sl_macro */
 #else
 # define data4_byte_flex                          0
 #endif

 #ifdef sl_data4_weight
 # define data4_weight                             1  /* sl_macro */
 #else
 # define data4_weight                             0
 #endif

 /* commands for regular use */
 #define data4_assign(src, dst)                    ((dst)->data4 = (src)->data4)  /* sl_macro */
 #define data4_assign_at(src, sat, dst)            ((dst)->data4 = &(src)->data4[(sat) * data4_size_c])  /* sl_macro */
 #define data4_null(e)                             ((e)->data4 = NULL)  /* sl_macro */
 #define data4_inc(e)                              ((e)->data4 += data4_size_c)  /* sl_macro */
 #define data4_dec(e)                              ((e)->data4 -= data4_size_c)  /* sl_macro */
 #define data4_add(e, n)                           ((e)->data4 += (n) * data4_size_c)  /* sl_macro */
 #define data4_sub(e, n)                           ((e)->data4 -= (n) * data4_size_c)  /* sl_macro */

 #define data4_copy(src, dst)                      sl_data4_copy((src)->data4, (dst)->data4)  /* sl_macro */
 #define data4_ncopy(src, dst, n)                  sl_data4_ncopy((src)->data4, (dst)->data4, n)  /* sl_macro */
 #define data4_nmove(src, dst, n)                  sl_data4_nmove((src)->data4, (dst)->data4, n)  /* sl_macro */

 #define data4_copy_at(src, sat, dst, dat)         sl_data4_copy(&(src)->data4[(sat) * data4_size_c], &(dst)->data4[(dat) * data4_size_c])  /* sl_macro */
 #define data4_ncopy_at(src, sat, dst, dat, n)     sl_data4_ncopy(&(src)->data4[(sat) * data4_size_c], &(dst)->data4[(dat) * data4_size_c], n)  /* sl_macro */
 #define data4_nmove_at(src, sat, dst, dat, n)     sl_data4_nmove(&(src)->data4[(sat) * data4_size_c], &(dst)->data4[(dat) * data4_size_c], n)  /* sl_macro */

 #define data4_xchange(e0, e1, t)                  (data4_copy(e0, t), data4_copy(e1, e0), data4_copy(t, e1))  /* sl_macro */
 #define data4_xchange_at(e0, at0, e1, at1, t)     (data4_copy_at(e0, at0, t, 0), data4_copy_at(e1, at1, e0, at0), data4_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data4_assign(src, dst)                 , data4_assign(src, dst)  /* sl_macro */
 #define cc_data4_assign_at(src, sat, dst)         , data4_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data4_null(e)                          , data4_null(e)  /* sl_macro */
 #define cc_data4_inc(e)                           , data4_inc(e)  /* sl_macro */
 #define cc_data4_dec(e)                           , data4_dec(e)  /* sl_macro */
 #define cc_data4_add(e, n)                        , data4_add(e, n)  /* sl_macro */
 #define cc_data4_sub(e, n)                        , data4_sub(e, n)  /* sl_macro */
 #define cc_data4_copy(src, dst)                   , data4_copy(src, dst)  /* sl_macro */
 #define cc_data4_ncopy(src, dst, n)               , data4_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data4_nmove(src, dst, n)               , data4_nmove(src, dst, n)  /* sl_macro */
 #define cc_data4_copy_at(src, sat, dst, dat)      , data4_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data4_ncopy_at(src, sat, dst, dat, n)  , data4_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data4_nmove_at(src, sat, dst, dat, n)  , data4_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data4_xchange(e0, e1, t)               , data4_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data4_xchange_at(e0, at0, e1, at1, t)  , data4_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data4_assign(src, dst)                 data4_assign(src, dst)
 #define cc_data4_assign_at(src, sat, dst)         data4_assign_at(src, sat, dst)
 #define cc_data4_null(e)                          data4_null(e)
 #define cc_data4_inc(e)                           data4_inc(e)
 #define cc_data4_dec(e)                           data4_dec(e)
 #define cc_data4_add(e, n)                        data4_add(e, n)
 #define cc_data4_sub(e, n)                        data4_sub(e, n)
 #define cc_data4_copy(src, dst)                   data4_copy(src, dst)
 #define cc_data4_ncopy(src, dst, n)               data4_ncopy(src, dst, n)
 #define cc_data4_nmove(src, dst, n)               data4_nmove(src, dst, n)
 #define cc_data4_copy_at(src, sat, dst, dat)      data4_copy_at(src, sat, dst, dat)
 #define cc_data4_ncopy_at(src, sat, dst, dat, n)  data4_ncopy_at(src, sat, dst, dat, n)
 #define cc_data4_nmove_at(src, sat, dst, dat, n)  data4_nmove_at(src, sat, dst, dat, n)
 #define cc_data4_xchange(e0, e1, t)               data4_xchange(e0, e1, t)
 #define cc_data4_xchange_at(e0, at0, e1, at1, t)  data4_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA4 */

 #define data4_n                                   0
 #define data4_byte                                0
/* #define data4_ptr(e)*/

 #define data4_byte_flex                           0
 #define data4_weight                              0

 /* commands for regular use */
 #define data4_assign(src, dst)                    Z_NOP()
 #define data4_assign_at(src, sat, dst)            Z_NOP()
 #define data4_null(e)                             Z_NOP()
 #define data4_inc(e)                              Z_NOP()
 #define data4_dec(e)                              Z_NOP()
 #define data4_add(e, n)                           Z_NOP()
 #define data4_sub(e, n)                           Z_NOP()
 #define data4_copy(src, dst)                      Z_NOP()
 #define data4_ncopy(src, dst, n)                  Z_NOP()
 #define data4_nmove(src, dst, n)                  Z_NOP()
 #define data4_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data4_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data4_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data4_xchange(e0, e1, t)                  Z_NOP()
 #define data4_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data4_assign(src, dst)
 #define cc_data4_assign_at(src, sat, dst)
 #define cc_data4_null(e)
 #define cc_data4_inc(e)
 #define cc_data4_dec(e)
 #define cc_data4_add(e, n)
 #define cc_data4_sub(e, n)
 #define cc_data4_copy(src, dst)
 #define cc_data4_ncopy(src, dst, n)
 #define cc_data4_nmove(src, dst, n)
 #define cc_data4_copy_at(src, sat, dst, dat)
 #define cc_data4_ncopy_at(src, sat, dst, dat, n)
 #define cc_data4_nmove_at(src, sat, dst, dat, n)
 #define cc_data4_xchange(e0, e1, t)
 #define cc_data4_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA4 */

#define data4_cm                                   SLCM_DATA4  /* sl_macro */


/* sl_macro SL_DATA5 SL_DATA5_IGNORE sl_data5_type_c sl_data5_size_c sl_data5_type_mpi sl_data5_size_mpi sl_data5_memcpy sl_data5_weight sl_data5_flex */

#ifdef SL_DATA5

 #define sl_data5_byte                             ((slint_t) (sl_data5_size_c) * sizeof(sl_data5_type_c))  /* sl_macro */

 #ifndef sl_data5_copy
  #if sl_data5_size_c <= 9 && !defined(sl_data5_memcpy)
   #if sl_data5_size_c == 1
    #define sl_data5_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data5_size_c == 2
    #define sl_data5_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data5_size_c == 3
    #define sl_data5_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data5_size_c == 4
    #define sl_data5_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data5_size_c == 5
    #define sl_data5_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data5_size_c == 6
    #define sl_data5_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data5_size_c == 7
    #define sl_data5_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data5_size_c == 8
    #define sl_data5_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data5_size_c == 9
    #define sl_data5_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data5_copy(src, dst)                 memcpy(dst, src, sl_data5_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data5_ncopy
  #define sl_data5_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data5_byte)  /* sl_macro */
 #endif
 #ifndef sl_data5_nmove
  #define sl_data5_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data5_byte)  /* sl_macro */
 #endif

 #define data5_type_c                              sl_data5_type_c  /* sl_macro */
 #define data5_size_c                              (sl_data5_size_c)  /* sl_macro */
 #define data5_type_mpi                            (sl_data5_type_mpi)  /* sl_macro */
 #define data5_size_mpi                            (sl_data5_size_mpi)  /* sl_macro */

 #define data5_idx                                 5  /* sl_macro */

 #define data5_n                                   1  /* sl_macro */
 #define data5_byte                                (sl_data5_byte)  /* sl_macro */
 #define data5_ptr(e)                              (e)->data5  /* sl_macro */

 #ifdef sl_data5_flex
 # define data5_byte_flex                          (sl_data5_byte)  /* sl_macro */
 #else
 # define data5_byte_flex                          0
 #endif

 #ifdef sl_data5_weight
 # define data5_weight                             1  /* sl_macro */
 #else
 # define data5_weight                             0
 #endif

 /* commands for regular use */
 #define data5_assign(src, dst)                    ((dst)->data5 = (src)->data5)  /* sl_macro */
 #define data5_assign_at(src, sat, dst)            ((dst)->data5 = &(src)->data5[(sat) * data5_size_c])  /* sl_macro */
 #define data5_null(e)                             ((e)->data5 = NULL)  /* sl_macro */
 #define data5_inc(e)                              ((e)->data5 += data5_size_c)  /* sl_macro */
 #define data5_dec(e)                              ((e)->data5 -= data5_size_c)  /* sl_macro */
 #define data5_add(e, n)                           ((e)->data5 += (n) * data5_size_c)  /* sl_macro */
 #define data5_sub(e, n)                           ((e)->data5 -= (n) * data5_size_c)  /* sl_macro */

 #define data5_copy(src, dst)                      sl_data5_copy((src)->data5, (dst)->data5)  /* sl_macro */
 #define data5_ncopy(src, dst, n)                  sl_data5_ncopy((src)->data5, (dst)->data5, n)  /* sl_macro */
 #define data5_nmove(src, dst, n)                  sl_data5_nmove((src)->data5, (dst)->data5, n)  /* sl_macro */

 #define data5_copy_at(src, sat, dst, dat)         sl_data5_copy(&(src)->data5[(sat) * data5_size_c], &(dst)->data5[(dat) * data5_size_c])  /* sl_macro */
 #define data5_ncopy_at(src, sat, dst, dat, n)     sl_data5_ncopy(&(src)->data5[(sat) * data5_size_c], &(dst)->data5[(dat) * data5_size_c], n)  /* sl_macro */
 #define data5_nmove_at(src, sat, dst, dat, n)     sl_data5_nmove(&(src)->data5[(sat) * data5_size_c], &(dst)->data5[(dat) * data5_size_c], n)  /* sl_macro */

 #define data5_xchange(e0, e1, t)                  (data5_copy(e0, t), data5_copy(e1, e0), data5_copy(t, e1))  /* sl_macro */
 #define data5_xchange_at(e0, at0, e1, at1, t)     (data5_copy_at(e0, at0, t, 0), data5_copy_at(e1, at1, e0, at0), data5_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data5_assign(src, dst)                 , data5_assign(src, dst)  /* sl_macro */
 #define cc_data5_assign_at(src, sat, dst)         , data5_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data5_null(e)                          , data5_null(e)  /* sl_macro */
 #define cc_data5_inc(e)                           , data5_inc(e)  /* sl_macro */
 #define cc_data5_dec(e)                           , data5_dec(e)  /* sl_macro */
 #define cc_data5_add(e, n)                        , data5_add(e, n)  /* sl_macro */
 #define cc_data5_sub(e, n)                        , data5_sub(e, n)  /* sl_macro */
 #define cc_data5_copy(src, dst)                   , data5_copy(src, dst)  /* sl_macro */
 #define cc_data5_ncopy(src, dst, n)               , data5_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data5_nmove(src, dst, n)               , data5_nmove(src, dst, n)  /* sl_macro */
 #define cc_data5_copy_at(src, sat, dst, dat)      , data5_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data5_ncopy_at(src, sat, dst, dat, n)  , data5_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data5_nmove_at(src, sat, dst, dat, n)  , data5_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data5_xchange(e0, e1, t)               , data5_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data5_xchange_at(e0, at0, e1, at1, t)  , data5_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data5_assign(src, dst)                 data5_assign(src, dst)
 #define cc_data5_assign_at(src, sat, dst)         data5_assign_at(src, sat, dst)
 #define cc_data5_null(e)                          data5_null(e)
 #define cc_data5_inc(e)                           data5_inc(e)
 #define cc_data5_dec(e)                           data5_dec(e)
 #define cc_data5_add(e, n)                        data5_add(e, n)
 #define cc_data5_sub(e, n)                        data5_sub(e, n)
 #define cc_data5_copy(src, dst)                   data5_copy(src, dst)
 #define cc_data5_ncopy(src, dst, n)               data5_ncopy(src, dst, n)
 #define cc_data5_nmove(src, dst, n)               data5_nmove(src, dst, n)
 #define cc_data5_copy_at(src, sat, dst, dat)      data5_copy_at(src, sat, dst, dat)
 #define cc_data5_ncopy_at(src, sat, dst, dat, n)  data5_ncopy_at(src, sat, dst, dat, n)
 #define cc_data5_nmove_at(src, sat, dst, dat, n)  data5_nmove_at(src, sat, dst, dat, n)
 #define cc_data5_xchange(e0, e1, t)               data5_xchange(e0, e1, t)
 #define cc_data5_xchange_at(e0, at0, e1, at1, t)  data5_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA5 */

 #define data5_n                                   0
 #define data5_byte                                0
/* #define data5_ptr(e)*/

 #define data5_byte_flex                           0
 #define data5_weight                              0

 /* commands for regular use */
 #define data5_assign(src, dst)                    Z_NOP()
 #define data5_assign_at(src, sat, dst)            Z_NOP()
 #define data5_null(e)                             Z_NOP()
 #define data5_inc(e)                              Z_NOP()
 #define data5_dec(e)                              Z_NOP()
 #define data5_add(e, n)                           Z_NOP()
 #define data5_sub(e, n)                           Z_NOP()
 #define data5_copy(src, dst)                      Z_NOP()
 #define data5_ncopy(src, dst, n)                  Z_NOP()
 #define data5_nmove(src, dst, n)                  Z_NOP()
 #define data5_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data5_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data5_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data5_xchange(e0, e1, t)                  Z_NOP()
 #define data5_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data5_assign(src, dst)
 #define cc_data5_assign_at(src, sat, dst)
 #define cc_data5_null(e)
 #define cc_data5_inc(e)
 #define cc_data5_dec(e)
 #define cc_data5_add(e, n)
 #define cc_data5_sub(e, n)
 #define cc_data5_copy(src, dst)
 #define cc_data5_ncopy(src, dst, n)
 #define cc_data5_nmove(src, dst, n)
 #define cc_data5_copy_at(src, sat, dst, dat)
 #define cc_data5_ncopy_at(src, sat, dst, dat, n)
 #define cc_data5_nmove_at(src, sat, dst, dat, n)
 #define cc_data5_xchange(e0, e1, t)
 #define cc_data5_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA5 */

#define data5_cm                                   SLCM_DATA5  /* sl_macro */


/* sl_macro SL_DATA6 SL_DATA6_IGNORE sl_data6_type_c sl_data6_size_c sl_data6_type_mpi sl_data6_size_mpi sl_data6_memcpy sl_data6_weight sl_data6_flex */

#ifdef SL_DATA6

 #define sl_data6_byte                             ((slint_t) (sl_data6_size_c) * sizeof(sl_data6_type_c))  /* sl_macro */

 #ifndef sl_data6_copy
  #if sl_data6_size_c <= 9 && !defined(sl_data6_memcpy)
   #if sl_data6_size_c == 1
    #define sl_data6_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data6_size_c == 2
    #define sl_data6_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data6_size_c == 3
    #define sl_data6_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data6_size_c == 4
    #define sl_data6_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data6_size_c == 5
    #define sl_data6_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data6_size_c == 6
    #define sl_data6_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data6_size_c == 7
    #define sl_data6_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data6_size_c == 8
    #define sl_data6_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data6_size_c == 9
    #define sl_data6_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data6_copy(src, dst)                 memcpy(dst, src, sl_data6_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data6_ncopy
  #define sl_data6_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data6_byte)  /* sl_macro */
 #endif
 #ifndef sl_data6_nmove
  #define sl_data6_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data6_byte)  /* sl_macro */
 #endif

 #define data6_type_c                              sl_data6_type_c  /* sl_macro */
 #define data6_size_c                              (sl_data6_size_c)  /* sl_macro */
 #define data6_type_mpi                            (sl_data6_type_mpi)  /* sl_macro */
 #define data6_size_mpi                            (sl_data6_size_mpi)  /* sl_macro */

 #define data6_idx                                 6  /* sl_macro */

 #define data6_n                                   1  /* sl_macro */
 #define data6_byte                                (sl_data6_byte)  /* sl_macro */
 #define data6_ptr(e)                              (e)->data6  /* sl_macro */

 #ifdef sl_data6_flex
 # define data6_byte_flex                          (sl_data6_byte)  /* sl_macro */
 #else
 # define data6_byte_flex                          0
 #endif

 #ifdef sl_data6_weight
 # define data6_weight                             1  /* sl_macro */
 #else
 # define data6_weight                             0
 #endif

 /* commands for regular use */
 #define data6_assign(src, dst)                    ((dst)->data6 = (src)->data6)  /* sl_macro */
 #define data6_assign_at(src, sat, dst)            ((dst)->data6 = &(src)->data6[(sat) * data6_size_c])  /* sl_macro */
 #define data6_null(e)                             ((e)->data6 = NULL)  /* sl_macro */
 #define data6_inc(e)                              ((e)->data6 += data6_size_c)  /* sl_macro */
 #define data6_dec(e)                              ((e)->data6 -= data6_size_c)  /* sl_macro */
 #define data6_add(e, n)                           ((e)->data6 += (n) * data6_size_c)  /* sl_macro */
 #define data6_sub(e, n)                           ((e)->data6 -= (n) * data6_size_c)  /* sl_macro */

 #define data6_copy(src, dst)                      sl_data6_copy((src)->data6, (dst)->data6)  /* sl_macro */
 #define data6_ncopy(src, dst, n)                  sl_data6_ncopy((src)->data6, (dst)->data6, n)  /* sl_macro */
 #define data6_nmove(src, dst, n)                  sl_data6_nmove((src)->data6, (dst)->data6, n)  /* sl_macro */

 #define data6_copy_at(src, sat, dst, dat)         sl_data6_copy(&(src)->data6[(sat) * data6_size_c], &(dst)->data6[(dat) * data6_size_c])  /* sl_macro */
 #define data6_ncopy_at(src, sat, dst, dat, n)     sl_data6_ncopy(&(src)->data6[(sat) * data6_size_c], &(dst)->data6[(dat) * data6_size_c], n)  /* sl_macro */
 #define data6_nmove_at(src, sat, dst, dat, n)     sl_data6_nmove(&(src)->data6[(sat) * data6_size_c], &(dst)->data6[(dat) * data6_size_c], n)  /* sl_macro */

 #define data6_xchange(e0, e1, t)                  (data6_copy(e0, t), data6_copy(e1, e0), data6_copy(t, e1))  /* sl_macro */
 #define data6_xchange_at(e0, at0, e1, at1, t)     (data6_copy_at(e0, at0, t, 0), data6_copy_at(e1, at1, e0, at0), data6_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data6_assign(src, dst)                 , data6_assign(src, dst)  /* sl_macro */
 #define cc_data6_assign_at(src, sat, dst)         , data6_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data6_null(e)                          , data6_null(e)  /* sl_macro */
 #define cc_data6_inc(e)                           , data6_inc(e)  /* sl_macro */
 #define cc_data6_dec(e)                           , data6_dec(e)  /* sl_macro */
 #define cc_data6_add(e, n)                        , data6_add(e, n)  /* sl_macro */
 #define cc_data6_sub(e, n)                        , data6_sub(e, n)  /* sl_macro */
 #define cc_data6_copy(src, dst)                   , data6_copy(src, dst)  /* sl_macro */
 #define cc_data6_ncopy(src, dst, n)               , data6_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data6_nmove(src, dst, n)               , data6_nmove(src, dst, n)  /* sl_macro */
 #define cc_data6_copy_at(src, sat, dst, dat)      , data6_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data6_ncopy_at(src, sat, dst, dat, n)  , data6_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data6_nmove_at(src, sat, dst, dat, n)  , data6_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data6_xchange(e0, e1, t)               , data6_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data6_xchange_at(e0, at0, e1, at1, t)  , data6_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data6_assign(src, dst)                 data6_assign(src, dst)
 #define cc_data6_assign_at(src, sat, dst)         data6_assign_at(src, sat, dst)
 #define cc_data6_null(e)                          data6_null(e)
 #define cc_data6_inc(e)                           data6_inc(e)
 #define cc_data6_dec(e)                           data6_dec(e)
 #define cc_data6_add(e, n)                        data6_add(e, n)
 #define cc_data6_sub(e, n)                        data6_sub(e, n)
 #define cc_data6_copy(src, dst)                   data6_copy(src, dst)
 #define cc_data6_ncopy(src, dst, n)               data6_ncopy(src, dst, n)
 #define cc_data6_nmove(src, dst, n)               data6_nmove(src, dst, n)
 #define cc_data6_copy_at(src, sat, dst, dat)      data6_copy_at(src, sat, dst, dat)
 #define cc_data6_ncopy_at(src, sat, dst, dat, n)  data6_ncopy_at(src, sat, dst, dat, n)
 #define cc_data6_nmove_at(src, sat, dst, dat, n)  data6_nmove_at(src, sat, dst, dat, n)
 #define cc_data6_xchange(e0, e1, t)               data6_xchange(e0, e1, t)
 #define cc_data6_xchange_at(e0, at0, e1, at1, t)  data6_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA6 */

 #define data6_n                                   0
 #define data6_byte                                0
/* #define data6_ptr(e)*/

 #define data6_byte_flex                           0
 #define data6_weight                              0

 /* commands for regular use */
 #define data6_assign(src, dst)                    Z_NOP()
 #define data6_assign_at(src, sat, dst)            Z_NOP()
 #define data6_null(e)                             Z_NOP()
 #define data6_inc(e)                              Z_NOP()
 #define data6_dec(e)                              Z_NOP()
 #define data6_add(e, n)                           Z_NOP()
 #define data6_sub(e, n)                           Z_NOP()
 #define data6_copy(src, dst)                      Z_NOP()
 #define data6_ncopy(src, dst, n)                  Z_NOP()
 #define data6_nmove(src, dst, n)                  Z_NOP()
 #define data6_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data6_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data6_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data6_xchange(e0, e1, t)                  Z_NOP()
 #define data6_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data6_assign(src, dst)
 #define cc_data6_assign_at(src, sat, dst)
 #define cc_data6_null(e)
 #define cc_data6_inc(e)
 #define cc_data6_dec(e)
 #define cc_data6_add(e, n)
 #define cc_data6_sub(e, n)
 #define cc_data6_copy(src, dst)
 #define cc_data6_ncopy(src, dst, n)
 #define cc_data6_nmove(src, dst, n)
 #define cc_data6_copy_at(src, sat, dst, dat)
 #define cc_data6_ncopy_at(src, sat, dst, dat, n)
 #define cc_data6_nmove_at(src, sat, dst, dat, n)
 #define cc_data6_xchange(e0, e1, t)
 #define cc_data6_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA6 */

#define data6_cm                                   SLCM_DATA6  /* sl_macro */


/* sl_macro SL_DATA7 SL_DATA7_IGNORE sl_data7_type_c sl_data7_size_c sl_data7_type_mpi sl_data7_size_mpi sl_data7_memcpy sl_data7_weight sl_data7_flex */

#ifdef SL_DATA7

 #define sl_data7_byte                             ((slint_t) (sl_data7_size_c) * sizeof(sl_data7_type_c))  /* sl_macro */

 #ifndef sl_data7_copy
  #if sl_data7_size_c <= 9 && !defined(sl_data7_memcpy)
   #if sl_data7_size_c == 1
    #define sl_data7_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data7_size_c == 2
    #define sl_data7_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data7_size_c == 3
    #define sl_data7_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data7_size_c == 4
    #define sl_data7_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data7_size_c == 5
    #define sl_data7_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data7_size_c == 6
    #define sl_data7_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data7_size_c == 7
    #define sl_data7_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data7_size_c == 8
    #define sl_data7_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data7_size_c == 9
    #define sl_data7_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data7_copy(src, dst)                 memcpy(dst, src, sl_data7_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data7_ncopy
  #define sl_data7_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data7_byte)  /* sl_macro */
 #endif
 #ifndef sl_data7_nmove
  #define sl_data7_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data7_byte)  /* sl_macro */
 #endif

 #define data7_type_c                              sl_data7_type_c  /* sl_macro */
 #define data7_size_c                              (sl_data7_size_c)  /* sl_macro */
 #define data7_type_mpi                            (sl_data7_type_mpi)  /* sl_macro */
 #define data7_size_mpi                            (sl_data7_size_mpi)  /* sl_macro */

 #define data7_idx                                 7  /* sl_macro */

 #define data7_n                                   1  /* sl_macro */
 #define data7_byte                                (sl_data7_byte)  /* sl_macro */
 #define data7_ptr(e)                              (e)->data7  /* sl_macro */

 #ifdef sl_data7_flex
 # define data7_byte_flex                          (sl_data7_byte)  /* sl_macro */
 #else
 # define data7_byte_flex                          0
 #endif

 #ifdef sl_data7_weight
 # define data7_weight                             1  /* sl_macro */
 #else
 # define data7_weight                             0
 #endif

 /* commands for regular use */
 #define data7_assign(src, dst)                    ((dst)->data7 = (src)->data7)  /* sl_macro */
 #define data7_assign_at(src, sat, dst)            ((dst)->data7 = &(src)->data7[(sat) * data7_size_c])  /* sl_macro */
 #define data7_null(e)                             ((e)->data7 = NULL)  /* sl_macro */
 #define data7_inc(e)                              ((e)->data7 += data7_size_c)  /* sl_macro */
 #define data7_dec(e)                              ((e)->data7 -= data7_size_c)  /* sl_macro */
 #define data7_add(e, n)                           ((e)->data7 += (n) * data7_size_c)  /* sl_macro */
 #define data7_sub(e, n)                           ((e)->data7 -= (n) * data7_size_c)  /* sl_macro */

 #define data7_copy(src, dst)                      sl_data7_copy((src)->data7, (dst)->data7)  /* sl_macro */
 #define data7_ncopy(src, dst, n)                  sl_data7_ncopy((src)->data7, (dst)->data7, n)  /* sl_macro */
 #define data7_nmove(src, dst, n)                  sl_data7_nmove((src)->data7, (dst)->data7, n)  /* sl_macro */

 #define data7_copy_at(src, sat, dst, dat)         sl_data7_copy(&(src)->data7[(sat) * data7_size_c], &(dst)->data7[(dat) * data7_size_c])  /* sl_macro */
 #define data7_ncopy_at(src, sat, dst, dat, n)     sl_data7_ncopy(&(src)->data7[(sat) * data7_size_c], &(dst)->data7[(dat) * data7_size_c], n)  /* sl_macro */
 #define data7_nmove_at(src, sat, dst, dat, n)     sl_data7_nmove(&(src)->data7[(sat) * data7_size_c], &(dst)->data7[(dat) * data7_size_c], n)  /* sl_macro */

 #define data7_xchange(e0, e1, t)                  (data7_copy(e0, t), data7_copy(e1, e0), data7_copy(t, e1))  /* sl_macro */
 #define data7_xchange_at(e0, at0, e1, at1, t)     (data7_copy_at(e0, at0, t, 0), data7_copy_at(e1, at1, e0, at0), data7_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data7_assign(src, dst)                 , data7_assign(src, dst)  /* sl_macro */
 #define cc_data7_assign_at(src, sat, dst)         , data7_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data7_null(e)                          , data7_null(e)  /* sl_macro */
 #define cc_data7_inc(e)                           , data7_inc(e)  /* sl_macro */
 #define cc_data7_dec(e)                           , data7_dec(e)  /* sl_macro */
 #define cc_data7_add(e, n)                        , data7_add(e, n)  /* sl_macro */
 #define cc_data7_sub(e, n)                        , data7_sub(e, n)  /* sl_macro */
 #define cc_data7_copy(src, dst)                   , data7_copy(src, dst)  /* sl_macro */
 #define cc_data7_ncopy(src, dst, n)               , data7_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data7_nmove(src, dst, n)               , data7_nmove(src, dst, n)  /* sl_macro */
 #define cc_data7_copy_at(src, sat, dst, dat)      , data7_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data7_ncopy_at(src, sat, dst, dat, n)  , data7_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data7_nmove_at(src, sat, dst, dat, n)  , data7_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data7_xchange(e0, e1, t)               , data7_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data7_xchange_at(e0, at0, e1, at1, t)  , data7_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data7_assign(src, dst)                 data7_assign(src, dst)
 #define cc_data7_assign_at(src, sat, dst)         data7_assign_at(src, sat, dst)
 #define cc_data7_null(e)                          data7_null(e)
 #define cc_data7_inc(e)                           data7_inc(e)
 #define cc_data7_dec(e)                           data7_dec(e)
 #define cc_data7_add(e, n)                        data7_add(e, n)
 #define cc_data7_sub(e, n)                        data7_sub(e, n)
 #define cc_data7_copy(src, dst)                   data7_copy(src, dst)
 #define cc_data7_ncopy(src, dst, n)               data7_ncopy(src, dst, n)
 #define cc_data7_nmove(src, dst, n)               data7_nmove(src, dst, n)
 #define cc_data7_copy_at(src, sat, dst, dat)      data7_copy_at(src, sat, dst, dat)
 #define cc_data7_ncopy_at(src, sat, dst, dat, n)  data7_ncopy_at(src, sat, dst, dat, n)
 #define cc_data7_nmove_at(src, sat, dst, dat, n)  data7_nmove_at(src, sat, dst, dat, n)
 #define cc_data7_xchange(e0, e1, t)               data7_xchange(e0, e1, t)
 #define cc_data7_xchange_at(e0, at0, e1, at1, t)  data7_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA7 */

 #define data7_n                                   0
 #define data7_byte                                0
/* #define data7_ptr(e)*/

 #define data7_byte_flex                           0
 #define data7_weight                              0

 /* commands for regular use */
 #define data7_assign(src, dst)                    Z_NOP()
 #define data7_assign_at(src, sat, dst)            Z_NOP()
 #define data7_null(e)                             Z_NOP()
 #define data7_inc(e)                              Z_NOP()
 #define data7_dec(e)                              Z_NOP()
 #define data7_add(e, n)                           Z_NOP()
 #define data7_sub(e, n)                           Z_NOP()
 #define data7_copy(src, dst)                      Z_NOP()
 #define data7_ncopy(src, dst, n)                  Z_NOP()
 #define data7_nmove(src, dst, n)                  Z_NOP()
 #define data7_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data7_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data7_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data7_xchange(e0, e1, t)                  Z_NOP()
 #define data7_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data7_assign(src, dst)
 #define cc_data7_assign_at(src, sat, dst)
 #define cc_data7_null(e)
 #define cc_data7_inc(e)
 #define cc_data7_dec(e)
 #define cc_data7_add(e, n)
 #define cc_data7_sub(e, n)
 #define cc_data7_copy(src, dst)
 #define cc_data7_ncopy(src, dst, n)
 #define cc_data7_nmove(src, dst, n)
 #define cc_data7_copy_at(src, sat, dst, dat)
 #define cc_data7_ncopy_at(src, sat, dst, dat, n)
 #define cc_data7_nmove_at(src, sat, dst, dat, n)
 #define cc_data7_xchange(e0, e1, t)
 #define cc_data7_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA7 */

#define data7_cm                                   SLCM_DATA7  /* sl_macro */


/* sl_macro SL_DATA8 SL_DATA8_IGNORE sl_data8_type_c sl_data8_size_c sl_data8_type_mpi sl_data8_size_mpi sl_data8_memcpy sl_data8_weight sl_data8_flex */

#ifdef SL_DATA8

 #define sl_data8_byte                             ((slint_t) (sl_data8_size_c) * sizeof(sl_data8_type_c))  /* sl_macro */

 #ifndef sl_data8_copy
  #if sl_data8_size_c <= 9 && !defined(sl_data8_memcpy)
   #if sl_data8_size_c == 1
    #define sl_data8_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data8_size_c == 2
    #define sl_data8_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data8_size_c == 3
    #define sl_data8_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data8_size_c == 4
    #define sl_data8_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data8_size_c == 5
    #define sl_data8_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data8_size_c == 6
    #define sl_data8_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data8_size_c == 7
    #define sl_data8_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data8_size_c == 8
    #define sl_data8_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data8_size_c == 9
    #define sl_data8_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data8_copy(src, dst)                 memcpy(dst, src, sl_data8_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data8_ncopy
  #define sl_data8_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data8_byte)  /* sl_macro */
 #endif
 #ifndef sl_data8_nmove
  #define sl_data8_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data8_byte)  /* sl_macro */
 #endif

 #define data8_type_c                              sl_data8_type_c  /* sl_macro */
 #define data8_size_c                              (sl_data8_size_c)  /* sl_macro */
 #define data8_type_mpi                            (sl_data8_type_mpi)  /* sl_macro */
 #define data8_size_mpi                            (sl_data8_size_mpi)  /* sl_macro */

 #define data8_idx                                 8  /* sl_macro */

 #define data8_n                                   1  /* sl_macro */
 #define data8_byte                                (sl_data8_byte)  /* sl_macro */
 #define data8_ptr(e)                              (e)->data8  /* sl_macro */

 #ifdef sl_data8_flex
 # define data8_byte_flex                          (sl_data8_byte)  /* sl_macro */
 #else
 # define data8_byte_flex                          0
 #endif

 #ifdef sl_data8_weight
 # define data8_weight                             1  /* sl_macro */
 #else
 # define data8_weight                             0
 #endif

 /* commands for regular use */
 #define data8_assign(src, dst)                    ((dst)->data8 = (src)->data8)  /* sl_macro */
 #define data8_assign_at(src, sat, dst)            ((dst)->data8 = &(src)->data8[(sat) * data8_size_c])  /* sl_macro */
 #define data8_null(e)                             ((e)->data8 = NULL)  /* sl_macro */
 #define data8_inc(e)                              ((e)->data8 += data8_size_c)  /* sl_macro */
 #define data8_dec(e)                              ((e)->data8 -= data8_size_c)  /* sl_macro */
 #define data8_add(e, n)                           ((e)->data8 += (n) * data8_size_c)  /* sl_macro */
 #define data8_sub(e, n)                           ((e)->data8 -= (n) * data8_size_c)  /* sl_macro */

 #define data8_copy(src, dst)                      sl_data8_copy((src)->data8, (dst)->data8)  /* sl_macro */
 #define data8_ncopy(src, dst, n)                  sl_data8_ncopy((src)->data8, (dst)->data8, n)  /* sl_macro */
 #define data8_nmove(src, dst, n)                  sl_data8_nmove((src)->data8, (dst)->data8, n)  /* sl_macro */

 #define data8_copy_at(src, sat, dst, dat)         sl_data8_copy(&(src)->data8[(sat) * data8_size_c], &(dst)->data8[(dat) * data8_size_c])  /* sl_macro */
 #define data8_ncopy_at(src, sat, dst, dat, n)     sl_data8_ncopy(&(src)->data8[(sat) * data8_size_c], &(dst)->data8[(dat) * data8_size_c], n)  /* sl_macro */
 #define data8_nmove_at(src, sat, dst, dat, n)     sl_data8_nmove(&(src)->data8[(sat) * data8_size_c], &(dst)->data8[(dat) * data8_size_c], n)  /* sl_macro */

 #define data8_xchange(e0, e1, t)                  (data8_copy(e0, t), data8_copy(e1, e0), data8_copy(t, e1))  /* sl_macro */
 #define data8_xchange_at(e0, at0, e1, at1, t)     (data8_copy_at(e0, at0, t, 0), data8_copy_at(e1, at1, e0, at0), data8_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data8_assign(src, dst)                 , data8_assign(src, dst)  /* sl_macro */
 #define cc_data8_assign_at(src, sat, dst)         , data8_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data8_null(e)                          , data8_null(e)  /* sl_macro */
 #define cc_data8_inc(e)                           , data8_inc(e)  /* sl_macro */
 #define cc_data8_dec(e)                           , data8_dec(e)  /* sl_macro */
 #define cc_data8_add(e, n)                        , data8_add(e, n)  /* sl_macro */
 #define cc_data8_sub(e, n)                        , data8_sub(e, n)  /* sl_macro */
 #define cc_data8_copy(src, dst)                   , data8_copy(src, dst)  /* sl_macro */
 #define cc_data8_ncopy(src, dst, n)               , data8_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data8_nmove(src, dst, n)               , data8_nmove(src, dst, n)  /* sl_macro */
 #define cc_data8_copy_at(src, sat, dst, dat)      , data8_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data8_ncopy_at(src, sat, dst, dat, n)  , data8_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data8_nmove_at(src, sat, dst, dat, n)  , data8_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data8_xchange(e0, e1, t)               , data8_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data8_xchange_at(e0, at0, e1, at1, t)  , data8_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data8_assign(src, dst)                 data8_assign(src, dst)
 #define cc_data8_assign_at(src, sat, dst)         data8_assign_at(src, sat, dst)
 #define cc_data8_null(e)                          data8_null(e)
 #define cc_data8_inc(e)                           data8_inc(e)
 #define cc_data8_dec(e)                           data8_dec(e)
 #define cc_data8_add(e, n)                        data8_add(e, n)
 #define cc_data8_sub(e, n)                        data8_sub(e, n)
 #define cc_data8_copy(src, dst)                   data8_copy(src, dst)
 #define cc_data8_ncopy(src, dst, n)               data8_ncopy(src, dst, n)
 #define cc_data8_nmove(src, dst, n)               data8_nmove(src, dst, n)
 #define cc_data8_copy_at(src, sat, dst, dat)      data8_copy_at(src, sat, dst, dat)
 #define cc_data8_ncopy_at(src, sat, dst, dat, n)  data8_ncopy_at(src, sat, dst, dat, n)
 #define cc_data8_nmove_at(src, sat, dst, dat, n)  data8_nmove_at(src, sat, dst, dat, n)
 #define cc_data8_xchange(e0, e1, t)               data8_xchange(e0, e1, t)
 #define cc_data8_xchange_at(e0, at0, e1, at1, t)  data8_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA8 */

 #define data8_n                                   0
 #define data8_byte                                0
/* #define data8_ptr(e)*/

 #define data8_byte_flex                           0
 #define data8_weight                              0

 /* commands for regular use */
 #define data8_assign(src, dst)                    Z_NOP()
 #define data8_assign_at(src, sat, dst)            Z_NOP()
 #define data8_null(e)                             Z_NOP()
 #define data8_inc(e)                              Z_NOP()
 #define data8_dec(e)                              Z_NOP()
 #define data8_add(e, n)                           Z_NOP()
 #define data8_sub(e, n)                           Z_NOP()
 #define data8_copy(src, dst)                      Z_NOP()
 #define data8_ncopy(src, dst, n)                  Z_NOP()
 #define data8_nmove(src, dst, n)                  Z_NOP()
 #define data8_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data8_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data8_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data8_xchange(e0, e1, t)                  Z_NOP()
 #define data8_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data8_assign(src, dst)
 #define cc_data8_assign_at(src, sat, dst)
 #define cc_data8_null(e)
 #define cc_data8_inc(e)
 #define cc_data8_dec(e)
 #define cc_data8_add(e, n)
 #define cc_data8_sub(e, n)
 #define cc_data8_copy(src, dst)
 #define cc_data8_ncopy(src, dst, n)
 #define cc_data8_nmove(src, dst, n)
 #define cc_data8_copy_at(src, sat, dst, dat)
 #define cc_data8_ncopy_at(src, sat, dst, dat, n)
 #define cc_data8_nmove_at(src, sat, dst, dat, n)
 #define cc_data8_xchange(e0, e1, t)
 #define cc_data8_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA8 */

#define data8_cm                                   SLCM_DATA8  /* sl_macro */


/* sl_macro SL_DATA9 SL_DATA9_IGNORE sl_data9_type_c sl_data9_size_c sl_data9_type_mpi sl_data9_size_mpi sl_data9_memcpy sl_data9_weight sl_data9_flex */

#ifdef SL_DATA9

 #define sl_data9_byte                             ((slint_t) (sl_data9_size_c) * sizeof(sl_data9_type_c))  /* sl_macro */

 #ifndef sl_data9_copy
  #if sl_data9_size_c <= 9 && !defined(sl_data9_memcpy)
   #if sl_data9_size_c == 1
    #define sl_data9_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data9_size_c == 2
    #define sl_data9_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data9_size_c == 3
    #define sl_data9_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data9_size_c == 4
    #define sl_data9_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data9_size_c == 5
    #define sl_data9_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data9_size_c == 6
    #define sl_data9_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data9_size_c == 7
    #define sl_data9_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data9_size_c == 8
    #define sl_data9_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data9_size_c == 9
    #define sl_data9_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data9_copy(src, dst)                 memcpy(dst, src, sl_data9_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data9_ncopy
  #define sl_data9_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data9_byte)  /* sl_macro */
 #endif
 #ifndef sl_data9_nmove
  #define sl_data9_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data9_byte)  /* sl_macro */
 #endif

 #define data9_type_c                              sl_data9_type_c  /* sl_macro */
 #define data9_size_c                              (sl_data9_size_c)  /* sl_macro */
 #define data9_type_mpi                            (sl_data9_type_mpi)  /* sl_macro */
 #define data9_size_mpi                            (sl_data9_size_mpi)  /* sl_macro */

 #define data9_idx                                 9  /* sl_macro */

 #define data9_n                                   1  /* sl_macro */
 #define data9_byte                                (sl_data9_byte)  /* sl_macro */
 #define data9_ptr(e)                              (e)->data9  /* sl_macro */

 #ifdef sl_data9_flex
 # define data9_byte_flex                          (sl_data9_byte)  /* sl_macro */
 #else
 # define data9_byte_flex                          0
 #endif

 #ifdef sl_data9_weight
 # define data9_weight                             1  /* sl_macro */
 #else
 # define data9_weight                             0
 #endif

 /* commands for regular use */
 #define data9_assign(src, dst)                    ((dst)->data9 = (src)->data9)  /* sl_macro */
 #define data9_assign_at(src, sat, dst)            ((dst)->data9 = &(src)->data9[(sat) * data9_size_c])  /* sl_macro */
 #define data9_null(e)                             ((e)->data9 = NULL)  /* sl_macro */
 #define data9_inc(e)                              ((e)->data9 += data9_size_c)  /* sl_macro */
 #define data9_dec(e)                              ((e)->data9 -= data9_size_c)  /* sl_macro */
 #define data9_add(e, n)                           ((e)->data9 += (n) * data9_size_c)  /* sl_macro */
 #define data9_sub(e, n)                           ((e)->data9 -= (n) * data9_size_c)  /* sl_macro */

 #define data9_copy(src, dst)                      sl_data9_copy((src)->data9, (dst)->data9)  /* sl_macro */
 #define data9_ncopy(src, dst, n)                  sl_data9_ncopy((src)->data9, (dst)->data9, n)  /* sl_macro */
 #define data9_nmove(src, dst, n)                  sl_data9_nmove((src)->data9, (dst)->data9, n)  /* sl_macro */

 #define data9_copy_at(src, sat, dst, dat)         sl_data9_copy(&(src)->data9[(sat) * data9_size_c], &(dst)->data9[(dat) * data9_size_c])  /* sl_macro */
 #define data9_ncopy_at(src, sat, dst, dat, n)     sl_data9_ncopy(&(src)->data9[(sat) * data9_size_c], &(dst)->data9[(dat) * data9_size_c], n)  /* sl_macro */
 #define data9_nmove_at(src, sat, dst, dat, n)     sl_data9_nmove(&(src)->data9[(sat) * data9_size_c], &(dst)->data9[(dat) * data9_size_c], n)  /* sl_macro */

 #define data9_xchange(e0, e1, t)                  (data9_copy(e0, t), data9_copy(e1, e0), data9_copy(t, e1))  /* sl_macro */
 #define data9_xchange_at(e0, at0, e1, at1, t)     (data9_copy_at(e0, at0, t, 0), data9_copy_at(e1, at1, e0, at0), data9_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data9_assign(src, dst)                 , data9_assign(src, dst)  /* sl_macro */
 #define cc_data9_assign_at(src, sat, dst)         , data9_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data9_null(e)                          , data9_null(e)  /* sl_macro */
 #define cc_data9_inc(e)                           , data9_inc(e)  /* sl_macro */
 #define cc_data9_dec(e)                           , data9_dec(e)  /* sl_macro */
 #define cc_data9_add(e, n)                        , data9_add(e, n)  /* sl_macro */
 #define cc_data9_sub(e, n)                        , data9_sub(e, n)  /* sl_macro */
 #define cc_data9_copy(src, dst)                   , data9_copy(src, dst)  /* sl_macro */
 #define cc_data9_ncopy(src, dst, n)               , data9_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data9_nmove(src, dst, n)               , data9_nmove(src, dst, n)  /* sl_macro */
 #define cc_data9_copy_at(src, sat, dst, dat)      , data9_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data9_ncopy_at(src, sat, dst, dat, n)  , data9_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data9_nmove_at(src, sat, dst, dat, n)  , data9_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data9_xchange(e0, e1, t)               , data9_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data9_xchange_at(e0, at0, e1, at1, t)  , data9_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data9_assign(src, dst)                 data9_assign(src, dst)
 #define cc_data9_assign_at(src, sat, dst)         data9_assign_at(src, sat, dst)
 #define cc_data9_null(e)                          data9_null(e)
 #define cc_data9_inc(e)                           data9_inc(e)
 #define cc_data9_dec(e)                           data9_dec(e)
 #define cc_data9_add(e, n)                        data9_add(e, n)
 #define cc_data9_sub(e, n)                        data9_sub(e, n)
 #define cc_data9_copy(src, dst)                   data9_copy(src, dst)
 #define cc_data9_ncopy(src, dst, n)               data9_ncopy(src, dst, n)
 #define cc_data9_nmove(src, dst, n)               data9_nmove(src, dst, n)
 #define cc_data9_copy_at(src, sat, dst, dat)      data9_copy_at(src, sat, dst, dat)
 #define cc_data9_ncopy_at(src, sat, dst, dat, n)  data9_ncopy_at(src, sat, dst, dat, n)
 #define cc_data9_nmove_at(src, sat, dst, dat, n)  data9_nmove_at(src, sat, dst, dat, n)
 #define cc_data9_xchange(e0, e1, t)               data9_xchange(e0, e1, t)
 #define cc_data9_xchange_at(e0, at0, e1, at1, t)  data9_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA9 */

 #define data9_n                                   0
 #define data9_byte                                0
/* #define data9_ptr(e)*/

 #define data9_byte_flex                           0
 #define data9_weight                              0

 /* commands for regular use */
 #define data9_assign(src, dst)                    Z_NOP()
 #define data9_assign_at(src, sat, dst)            Z_NOP()
 #define data9_null(e)                             Z_NOP()
 #define data9_inc(e)                              Z_NOP()
 #define data9_dec(e)                              Z_NOP()
 #define data9_add(e, n)                           Z_NOP()
 #define data9_sub(e, n)                           Z_NOP()
 #define data9_copy(src, dst)                      Z_NOP()
 #define data9_ncopy(src, dst, n)                  Z_NOP()
 #define data9_nmove(src, dst, n)                  Z_NOP()
 #define data9_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data9_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data9_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data9_xchange(e0, e1, t)                  Z_NOP()
 #define data9_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data9_assign(src, dst)
 #define cc_data9_assign_at(src, sat, dst)
 #define cc_data9_null(e)
 #define cc_data9_inc(e)
 #define cc_data9_dec(e)
 #define cc_data9_add(e, n)
 #define cc_data9_sub(e, n)
 #define cc_data9_copy(src, dst)
 #define cc_data9_ncopy(src, dst, n)
 #define cc_data9_nmove(src, dst, n)
 #define cc_data9_copy_at(src, sat, dst, dat)
 #define cc_data9_ncopy_at(src, sat, dst, dat, n)
 #define cc_data9_nmove_at(src, sat, dst, dat, n)
 #define cc_data9_xchange(e0, e1, t)
 #define cc_data9_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA9 */

#define data9_cm                                   SLCM_DATA9  /* sl_macro */


/* sl_macro SL_DATA10 SL_DATA10_IGNORE sl_data10_type_c sl_data10_size_c sl_data10_type_mpi sl_data10_size_mpi sl_data10_memcpy sl_data10_weight sl_data10_flex */

#ifdef SL_DATA10

 #define sl_data10_byte                             ((slint_t) (sl_data10_size_c) * sizeof(sl_data10_type_c))  /* sl_macro */

 #ifndef sl_data10_copy
  #if sl_data10_size_c <= 9 && !defined(sl_data10_memcpy)
   #if sl_data10_size_c == 1
    #define sl_data10_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data10_size_c == 2
    #define sl_data10_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data10_size_c == 3
    #define sl_data10_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data10_size_c == 4
    #define sl_data10_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data10_size_c == 5
    #define sl_data10_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data10_size_c == 6
    #define sl_data10_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data10_size_c == 7
    #define sl_data10_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data10_size_c == 8
    #define sl_data10_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data10_size_c == 9
    #define sl_data10_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data10_copy(src, dst)                 memcpy(dst, src, sl_data10_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data10_ncopy
  #define sl_data10_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data10_byte)  /* sl_macro */
 #endif
 #ifndef sl_data10_nmove
  #define sl_data10_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data10_byte)  /* sl_macro */
 #endif

 #define data10_type_c                              sl_data10_type_c  /* sl_macro */
 #define data10_size_c                              (sl_data10_size_c)  /* sl_macro */
 #define data10_type_mpi                            (sl_data10_type_mpi)  /* sl_macro */
 #define data10_size_mpi                            (sl_data10_size_mpi)  /* sl_macro */

 #define data10_idx                                 10  /* sl_macro */

 #define data10_n                                   1  /* sl_macro */
 #define data10_byte                                (sl_data10_byte)  /* sl_macro */
 #define data10_ptr(e)                              (e)->data10  /* sl_macro */

 #ifdef sl_data10_flex
 # define data10_byte_flex                          (sl_data10_byte)  /* sl_macro */
 #else
 # define data10_byte_flex                          0
 #endif

 #ifdef sl_data10_weight
 # define data10_weight                             1  /* sl_macro */
 #else
 # define data10_weight                             0
 #endif

 /* commands for regular use */
 #define data10_assign(src, dst)                    ((dst)->data10 = (src)->data10)  /* sl_macro */
 #define data10_assign_at(src, sat, dst)            ((dst)->data10 = &(src)->data10[(sat) * data10_size_c])  /* sl_macro */
 #define data10_null(e)                             ((e)->data10 = NULL)  /* sl_macro */
 #define data10_inc(e)                              ((e)->data10 += data10_size_c)  /* sl_macro */
 #define data10_dec(e)                              ((e)->data10 -= data10_size_c)  /* sl_macro */
 #define data10_add(e, n)                           ((e)->data10 += (n) * data10_size_c)  /* sl_macro */
 #define data10_sub(e, n)                           ((e)->data10 -= (n) * data10_size_c)  /* sl_macro */

 #define data10_copy(src, dst)                      sl_data10_copy((src)->data10, (dst)->data10)  /* sl_macro */
 #define data10_ncopy(src, dst, n)                  sl_data10_ncopy((src)->data10, (dst)->data10, n)  /* sl_macro */
 #define data10_nmove(src, dst, n)                  sl_data10_nmove((src)->data10, (dst)->data10, n)  /* sl_macro */

 #define data10_copy_at(src, sat, dst, dat)         sl_data10_copy(&(src)->data10[(sat) * data10_size_c], &(dst)->data10[(dat) * data10_size_c])  /* sl_macro */
 #define data10_ncopy_at(src, sat, dst, dat, n)     sl_data10_ncopy(&(src)->data10[(sat) * data10_size_c], &(dst)->data10[(dat) * data10_size_c], n)  /* sl_macro */
 #define data10_nmove_at(src, sat, dst, dat, n)     sl_data10_nmove(&(src)->data10[(sat) * data10_size_c], &(dst)->data10[(dat) * data10_size_c], n)  /* sl_macro */

 #define data10_xchange(e0, e1, t)                  (data10_copy(e0, t), data10_copy(e1, e0), data10_copy(t, e1))  /* sl_macro */
 #define data10_xchange_at(e0, at0, e1, at1, t)     (data10_copy_at(e0, at0, t, 0), data10_copy_at(e1, at1, e0, at0), data10_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data10_assign(src, dst)                 , data10_assign(src, dst)  /* sl_macro */
 #define cc_data10_assign_at(src, sat, dst)         , data10_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data10_null(e)                          , data10_null(e)  /* sl_macro */
 #define cc_data10_inc(e)                           , data10_inc(e)  /* sl_macro */
 #define cc_data10_dec(e)                           , data10_dec(e)  /* sl_macro */
 #define cc_data10_add(e, n)                        , data10_add(e, n)  /* sl_macro */
 #define cc_data10_sub(e, n)                        , data10_sub(e, n)  /* sl_macro */
 #define cc_data10_copy(src, dst)                   , data10_copy(src, dst)  /* sl_macro */
 #define cc_data10_ncopy(src, dst, n)               , data10_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data10_nmove(src, dst, n)               , data10_nmove(src, dst, n)  /* sl_macro */
 #define cc_data10_copy_at(src, sat, dst, dat)      , data10_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data10_ncopy_at(src, sat, dst, dat, n)  , data10_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data10_nmove_at(src, sat, dst, dat, n)  , data10_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data10_xchange(e0, e1, t)               , data10_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data10_xchange_at(e0, at0, e1, at1, t)  , data10_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data10_assign(src, dst)                 data10_assign(src, dst)
 #define cc_data10_assign_at(src, sat, dst)         data10_assign_at(src, sat, dst)
 #define cc_data10_null(e)                          data10_null(e)
 #define cc_data10_inc(e)                           data10_inc(e)
 #define cc_data10_dec(e)                           data10_dec(e)
 #define cc_data10_add(e, n)                        data10_add(e, n)
 #define cc_data10_sub(e, n)                        data10_sub(e, n)
 #define cc_data10_copy(src, dst)                   data10_copy(src, dst)
 #define cc_data10_ncopy(src, dst, n)               data10_ncopy(src, dst, n)
 #define cc_data10_nmove(src, dst, n)               data10_nmove(src, dst, n)
 #define cc_data10_copy_at(src, sat, dst, dat)      data10_copy_at(src, sat, dst, dat)
 #define cc_data10_ncopy_at(src, sat, dst, dat, n)  data10_ncopy_at(src, sat, dst, dat, n)
 #define cc_data10_nmove_at(src, sat, dst, dat, n)  data10_nmove_at(src, sat, dst, dat, n)
 #define cc_data10_xchange(e0, e1, t)               data10_xchange(e0, e1, t)
 #define cc_data10_xchange_at(e0, at0, e1, at1, t)  data10_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA10 */

 #define data10_n                                   0
 #define data10_byte                                0
/* #define data10_ptr(e)*/

 #define data10_byte_flex                           0
 #define data10_weight                              0

 /* commands for regular use */
 #define data10_assign(src, dst)                    Z_NOP()
 #define data10_assign_at(src, sat, dst)            Z_NOP()
 #define data10_null(e)                             Z_NOP()
 #define data10_inc(e)                              Z_NOP()
 #define data10_dec(e)                              Z_NOP()
 #define data10_add(e, n)                           Z_NOP()
 #define data10_sub(e, n)                           Z_NOP()
 #define data10_copy(src, dst)                      Z_NOP()
 #define data10_ncopy(src, dst, n)                  Z_NOP()
 #define data10_nmove(src, dst, n)                  Z_NOP()
 #define data10_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data10_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data10_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data10_xchange(e0, e1, t)                  Z_NOP()
 #define data10_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data10_assign(src, dst)
 #define cc_data10_assign_at(src, sat, dst)
 #define cc_data10_null(e)
 #define cc_data10_inc(e)
 #define cc_data10_dec(e)
 #define cc_data10_add(e, n)
 #define cc_data10_sub(e, n)
 #define cc_data10_copy(src, dst)
 #define cc_data10_ncopy(src, dst, n)
 #define cc_data10_nmove(src, dst, n)
 #define cc_data10_copy_at(src, sat, dst, dat)
 #define cc_data10_ncopy_at(src, sat, dst, dat, n)
 #define cc_data10_nmove_at(src, sat, dst, dat, n)
 #define cc_data10_xchange(e0, e1, t)
 #define cc_data10_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA10 */

#define data10_cm                                   SLCM_DATA10  /* sl_macro */


/* sl_macro SL_DATA11 SL_DATA11_IGNORE sl_data11_type_c sl_data11_size_c sl_data11_type_mpi sl_data11_size_mpi sl_data11_memcpy sl_data11_weight sl_data11_flex */

#ifdef SL_DATA11

 #define sl_data11_byte                             ((slint_t) (sl_data11_size_c) * sizeof(sl_data11_type_c))  /* sl_macro */

 #ifndef sl_data11_copy
  #if sl_data11_size_c <= 9 && !defined(sl_data11_memcpy)
   #if sl_data11_size_c == 1
    #define sl_data11_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data11_size_c == 2
    #define sl_data11_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data11_size_c == 3
    #define sl_data11_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data11_size_c == 4
    #define sl_data11_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data11_size_c == 5
    #define sl_data11_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data11_size_c == 6
    #define sl_data11_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data11_size_c == 7
    #define sl_data11_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data11_size_c == 8
    #define sl_data11_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data11_size_c == 9
    #define sl_data11_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data11_copy(src, dst)                 memcpy(dst, src, sl_data11_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data11_ncopy
  #define sl_data11_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data11_byte)  /* sl_macro */
 #endif
 #ifndef sl_data11_nmove
  #define sl_data11_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data11_byte)  /* sl_macro */
 #endif

 #define data11_type_c                              sl_data11_type_c  /* sl_macro */
 #define data11_size_c                              (sl_data11_size_c)  /* sl_macro */
 #define data11_type_mpi                            (sl_data11_type_mpi)  /* sl_macro */
 #define data11_size_mpi                            (sl_data11_size_mpi)  /* sl_macro */

 #define data11_idx                                 11  /* sl_macro */

 #define data11_n                                   1  /* sl_macro */
 #define data11_byte                                (sl_data11_byte)  /* sl_macro */
 #define data11_ptr(e)                              (e)->data11  /* sl_macro */

 #ifdef sl_data11_flex
 # define data11_byte_flex                          (sl_data11_byte)  /* sl_macro */
 #else
 # define data11_byte_flex                          0
 #endif

 #ifdef sl_data11_weight
 # define data11_weight                             1  /* sl_macro */
 #else
 # define data11_weight                             0
 #endif

 /* commands for regular use */
 #define data11_assign(src, dst)                    ((dst)->data11 = (src)->data11)  /* sl_macro */
 #define data11_assign_at(src, sat, dst)            ((dst)->data11 = &(src)->data11[(sat) * data11_size_c])  /* sl_macro */
 #define data11_null(e)                             ((e)->data11 = NULL)  /* sl_macro */
 #define data11_inc(e)                              ((e)->data11 += data11_size_c)  /* sl_macro */
 #define data11_dec(e)                              ((e)->data11 -= data11_size_c)  /* sl_macro */
 #define data11_add(e, n)                           ((e)->data11 += (n) * data11_size_c)  /* sl_macro */
 #define data11_sub(e, n)                           ((e)->data11 -= (n) * data11_size_c)  /* sl_macro */

 #define data11_copy(src, dst)                      sl_data11_copy((src)->data11, (dst)->data11)  /* sl_macro */
 #define data11_ncopy(src, dst, n)                  sl_data11_ncopy((src)->data11, (dst)->data11, n)  /* sl_macro */
 #define data11_nmove(src, dst, n)                  sl_data11_nmove((src)->data11, (dst)->data11, n)  /* sl_macro */

 #define data11_copy_at(src, sat, dst, dat)         sl_data11_copy(&(src)->data11[(sat) * data11_size_c], &(dst)->data11[(dat) * data11_size_c])  /* sl_macro */
 #define data11_ncopy_at(src, sat, dst, dat, n)     sl_data11_ncopy(&(src)->data11[(sat) * data11_size_c], &(dst)->data11[(dat) * data11_size_c], n)  /* sl_macro */
 #define data11_nmove_at(src, sat, dst, dat, n)     sl_data11_nmove(&(src)->data11[(sat) * data11_size_c], &(dst)->data11[(dat) * data11_size_c], n)  /* sl_macro */

 #define data11_xchange(e0, e1, t)                  (data11_copy(e0, t), data11_copy(e1, e0), data11_copy(t, e1))  /* sl_macro */
 #define data11_xchange_at(e0, at0, e1, at1, t)     (data11_copy_at(e0, at0, t, 0), data11_copy_at(e1, at1, e0, at0), data11_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data11_assign(src, dst)                 , data11_assign(src, dst)  /* sl_macro */
 #define cc_data11_assign_at(src, sat, dst)         , data11_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data11_null(e)                          , data11_null(e)  /* sl_macro */
 #define cc_data11_inc(e)                           , data11_inc(e)  /* sl_macro */
 #define cc_data11_dec(e)                           , data11_dec(e)  /* sl_macro */
 #define cc_data11_add(e, n)                        , data11_add(e, n)  /* sl_macro */
 #define cc_data11_sub(e, n)                        , data11_sub(e, n)  /* sl_macro */
 #define cc_data11_copy(src, dst)                   , data11_copy(src, dst)  /* sl_macro */
 #define cc_data11_ncopy(src, dst, n)               , data11_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data11_nmove(src, dst, n)               , data11_nmove(src, dst, n)  /* sl_macro */
 #define cc_data11_copy_at(src, sat, dst, dat)      , data11_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data11_ncopy_at(src, sat, dst, dat, n)  , data11_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data11_nmove_at(src, sat, dst, dat, n)  , data11_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data11_xchange(e0, e1, t)               , data11_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data11_xchange_at(e0, at0, e1, at1, t)  , data11_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data11_assign(src, dst)                 data11_assign(src, dst)
 #define cc_data11_assign_at(src, sat, dst)         data11_assign_at(src, sat, dst)
 #define cc_data11_null(e)                          data11_null(e)
 #define cc_data11_inc(e)                           data11_inc(e)
 #define cc_data11_dec(e)                           data11_dec(e)
 #define cc_data11_add(e, n)                        data11_add(e, n)
 #define cc_data11_sub(e, n)                        data11_sub(e, n)
 #define cc_data11_copy(src, dst)                   data11_copy(src, dst)
 #define cc_data11_ncopy(src, dst, n)               data11_ncopy(src, dst, n)
 #define cc_data11_nmove(src, dst, n)               data11_nmove(src, dst, n)
 #define cc_data11_copy_at(src, sat, dst, dat)      data11_copy_at(src, sat, dst, dat)
 #define cc_data11_ncopy_at(src, sat, dst, dat, n)  data11_ncopy_at(src, sat, dst, dat, n)
 #define cc_data11_nmove_at(src, sat, dst, dat, n)  data11_nmove_at(src, sat, dst, dat, n)
 #define cc_data11_xchange(e0, e1, t)               data11_xchange(e0, e1, t)
 #define cc_data11_xchange_at(e0, at0, e1, at1, t)  data11_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA11 */

 #define data11_n                                   0
 #define data11_byte                                0
/* #define data11_ptr(e)*/

 #define data11_byte_flex                           0
 #define data11_weight                              0

 /* commands for regular use */
 #define data11_assign(src, dst)                    Z_NOP()
 #define data11_assign_at(src, sat, dst)            Z_NOP()
 #define data11_null(e)                             Z_NOP()
 #define data11_inc(e)                              Z_NOP()
 #define data11_dec(e)                              Z_NOP()
 #define data11_add(e, n)                           Z_NOP()
 #define data11_sub(e, n)                           Z_NOP()
 #define data11_copy(src, dst)                      Z_NOP()
 #define data11_ncopy(src, dst, n)                  Z_NOP()
 #define data11_nmove(src, dst, n)                  Z_NOP()
 #define data11_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data11_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data11_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data11_xchange(e0, e1, t)                  Z_NOP()
 #define data11_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data11_assign(src, dst)
 #define cc_data11_assign_at(src, sat, dst)
 #define cc_data11_null(e)
 #define cc_data11_inc(e)
 #define cc_data11_dec(e)
 #define cc_data11_add(e, n)
 #define cc_data11_sub(e, n)
 #define cc_data11_copy(src, dst)
 #define cc_data11_ncopy(src, dst, n)
 #define cc_data11_nmove(src, dst, n)
 #define cc_data11_copy_at(src, sat, dst, dat)
 #define cc_data11_ncopy_at(src, sat, dst, dat, n)
 #define cc_data11_nmove_at(src, sat, dst, dat, n)
 #define cc_data11_xchange(e0, e1, t)
 #define cc_data11_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA11 */

#define data11_cm                                   SLCM_DATA11  /* sl_macro */


/* sl_macro SL_DATA12 SL_DATA12_IGNORE sl_data12_type_c sl_data12_size_c sl_data12_type_mpi sl_data12_size_mpi sl_data12_memcpy sl_data12_weight sl_data12_flex */

#ifdef SL_DATA12

 #define sl_data12_byte                             ((slint_t) (sl_data12_size_c) * sizeof(sl_data12_type_c))  /* sl_macro */

 #ifndef sl_data12_copy
  #if sl_data12_size_c <= 9 && !defined(sl_data12_memcpy)
   #if sl_data12_size_c == 1
    #define sl_data12_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data12_size_c == 2
    #define sl_data12_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data12_size_c == 3
    #define sl_data12_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data12_size_c == 4
    #define sl_data12_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data12_size_c == 5
    #define sl_data12_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data12_size_c == 6
    #define sl_data12_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data12_size_c == 7
    #define sl_data12_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data12_size_c == 8
    #define sl_data12_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data12_size_c == 9
    #define sl_data12_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data12_copy(src, dst)                 memcpy(dst, src, sl_data12_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data12_ncopy
  #define sl_data12_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data12_byte)  /* sl_macro */
 #endif
 #ifndef sl_data12_nmove
  #define sl_data12_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data12_byte)  /* sl_macro */
 #endif

 #define data12_type_c                              sl_data12_type_c  /* sl_macro */
 #define data12_size_c                              (sl_data12_size_c)  /* sl_macro */
 #define data12_type_mpi                            (sl_data12_type_mpi)  /* sl_macro */
 #define data12_size_mpi                            (sl_data12_size_mpi)  /* sl_macro */

 #define data12_idx                                 12  /* sl_macro */

 #define data12_n                                   1  /* sl_macro */
 #define data12_byte                                (sl_data12_byte)  /* sl_macro */
 #define data12_ptr(e)                              (e)->data12  /* sl_macro */

 #ifdef sl_data12_flex
 # define data12_byte_flex                          (sl_data12_byte)  /* sl_macro */
 #else
 # define data12_byte_flex                          0
 #endif

 #ifdef sl_data12_weight
 # define data12_weight                             1  /* sl_macro */
 #else
 # define data12_weight                             0
 #endif

 /* commands for regular use */
 #define data12_assign(src, dst)                    ((dst)->data12 = (src)->data12)  /* sl_macro */
 #define data12_assign_at(src, sat, dst)            ((dst)->data12 = &(src)->data12[(sat) * data12_size_c])  /* sl_macro */
 #define data12_null(e)                             ((e)->data12 = NULL)  /* sl_macro */
 #define data12_inc(e)                              ((e)->data12 += data12_size_c)  /* sl_macro */
 #define data12_dec(e)                              ((e)->data12 -= data12_size_c)  /* sl_macro */
 #define data12_add(e, n)                           ((e)->data12 += (n) * data12_size_c)  /* sl_macro */
 #define data12_sub(e, n)                           ((e)->data12 -= (n) * data12_size_c)  /* sl_macro */

 #define data12_copy(src, dst)                      sl_data12_copy((src)->data12, (dst)->data12)  /* sl_macro */
 #define data12_ncopy(src, dst, n)                  sl_data12_ncopy((src)->data12, (dst)->data12, n)  /* sl_macro */
 #define data12_nmove(src, dst, n)                  sl_data12_nmove((src)->data12, (dst)->data12, n)  /* sl_macro */

 #define data12_copy_at(src, sat, dst, dat)         sl_data12_copy(&(src)->data12[(sat) * data12_size_c], &(dst)->data12[(dat) * data12_size_c])  /* sl_macro */
 #define data12_ncopy_at(src, sat, dst, dat, n)     sl_data12_ncopy(&(src)->data12[(sat) * data12_size_c], &(dst)->data12[(dat) * data12_size_c], n)  /* sl_macro */
 #define data12_nmove_at(src, sat, dst, dat, n)     sl_data12_nmove(&(src)->data12[(sat) * data12_size_c], &(dst)->data12[(dat) * data12_size_c], n)  /* sl_macro */

 #define data12_xchange(e0, e1, t)                  (data12_copy(e0, t), data12_copy(e1, e0), data12_copy(t, e1))  /* sl_macro */
 #define data12_xchange_at(e0, at0, e1, at1, t)     (data12_copy_at(e0, at0, t, 0), data12_copy_at(e1, at1, e0, at0), data12_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data12_assign(src, dst)                 , data12_assign(src, dst)  /* sl_macro */
 #define cc_data12_assign_at(src, sat, dst)         , data12_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data12_null(e)                          , data12_null(e)  /* sl_macro */
 #define cc_data12_inc(e)                           , data12_inc(e)  /* sl_macro */
 #define cc_data12_dec(e)                           , data12_dec(e)  /* sl_macro */
 #define cc_data12_add(e, n)                        , data12_add(e, n)  /* sl_macro */
 #define cc_data12_sub(e, n)                        , data12_sub(e, n)  /* sl_macro */
 #define cc_data12_copy(src, dst)                   , data12_copy(src, dst)  /* sl_macro */
 #define cc_data12_ncopy(src, dst, n)               , data12_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data12_nmove(src, dst, n)               , data12_nmove(src, dst, n)  /* sl_macro */
 #define cc_data12_copy_at(src, sat, dst, dat)      , data12_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data12_ncopy_at(src, sat, dst, dat, n)  , data12_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data12_nmove_at(src, sat, dst, dat, n)  , data12_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data12_xchange(e0, e1, t)               , data12_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data12_xchange_at(e0, at0, e1, at1, t)  , data12_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data12_assign(src, dst)                 data12_assign(src, dst)
 #define cc_data12_assign_at(src, sat, dst)         data12_assign_at(src, sat, dst)
 #define cc_data12_null(e)                          data12_null(e)
 #define cc_data12_inc(e)                           data12_inc(e)
 #define cc_data12_dec(e)                           data12_dec(e)
 #define cc_data12_add(e, n)                        data12_add(e, n)
 #define cc_data12_sub(e, n)                        data12_sub(e, n)
 #define cc_data12_copy(src, dst)                   data12_copy(src, dst)
 #define cc_data12_ncopy(src, dst, n)               data12_ncopy(src, dst, n)
 #define cc_data12_nmove(src, dst, n)               data12_nmove(src, dst, n)
 #define cc_data12_copy_at(src, sat, dst, dat)      data12_copy_at(src, sat, dst, dat)
 #define cc_data12_ncopy_at(src, sat, dst, dat, n)  data12_ncopy_at(src, sat, dst, dat, n)
 #define cc_data12_nmove_at(src, sat, dst, dat, n)  data12_nmove_at(src, sat, dst, dat, n)
 #define cc_data12_xchange(e0, e1, t)               data12_xchange(e0, e1, t)
 #define cc_data12_xchange_at(e0, at0, e1, at1, t)  data12_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA12 */

 #define data12_n                                   0
 #define data12_byte                                0
/* #define data12_ptr(e)*/

 #define data12_byte_flex                           0
 #define data12_weight                              0

 /* commands for regular use */
 #define data12_assign(src, dst)                    Z_NOP()
 #define data12_assign_at(src, sat, dst)            Z_NOP()
 #define data12_null(e)                             Z_NOP()
 #define data12_inc(e)                              Z_NOP()
 #define data12_dec(e)                              Z_NOP()
 #define data12_add(e, n)                           Z_NOP()
 #define data12_sub(e, n)                           Z_NOP()
 #define data12_copy(src, dst)                      Z_NOP()
 #define data12_ncopy(src, dst, n)                  Z_NOP()
 #define data12_nmove(src, dst, n)                  Z_NOP()
 #define data12_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data12_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data12_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data12_xchange(e0, e1, t)                  Z_NOP()
 #define data12_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data12_assign(src, dst)
 #define cc_data12_assign_at(src, sat, dst)
 #define cc_data12_null(e)
 #define cc_data12_inc(e)
 #define cc_data12_dec(e)
 #define cc_data12_add(e, n)
 #define cc_data12_sub(e, n)
 #define cc_data12_copy(src, dst)
 #define cc_data12_ncopy(src, dst, n)
 #define cc_data12_nmove(src, dst, n)
 #define cc_data12_copy_at(src, sat, dst, dat)
 #define cc_data12_ncopy_at(src, sat, dst, dat, n)
 #define cc_data12_nmove_at(src, sat, dst, dat, n)
 #define cc_data12_xchange(e0, e1, t)
 #define cc_data12_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA12 */

#define data12_cm                                   SLCM_DATA12  /* sl_macro */


/* sl_macro SL_DATA13 SL_DATA13_IGNORE sl_data13_type_c sl_data13_size_c sl_data13_type_mpi sl_data13_size_mpi sl_data13_memcpy sl_data13_weight sl_data13_flex */

#ifdef SL_DATA13

 #define sl_data13_byte                             ((slint_t) (sl_data13_size_c) * sizeof(sl_data13_type_c))  /* sl_macro */

 #ifndef sl_data13_copy
  #if sl_data13_size_c <= 9 && !defined(sl_data13_memcpy)
   #if sl_data13_size_c == 1
    #define sl_data13_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data13_size_c == 2
    #define sl_data13_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data13_size_c == 3
    #define sl_data13_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data13_size_c == 4
    #define sl_data13_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data13_size_c == 5
    #define sl_data13_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data13_size_c == 6
    #define sl_data13_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data13_size_c == 7
    #define sl_data13_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data13_size_c == 8
    #define sl_data13_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data13_size_c == 9
    #define sl_data13_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data13_copy(src, dst)                 memcpy(dst, src, sl_data13_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data13_ncopy
  #define sl_data13_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data13_byte)  /* sl_macro */
 #endif
 #ifndef sl_data13_nmove
  #define sl_data13_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data13_byte)  /* sl_macro */
 #endif

 #define data13_type_c                              sl_data13_type_c  /* sl_macro */
 #define data13_size_c                              (sl_data13_size_c)  /* sl_macro */
 #define data13_type_mpi                            (sl_data13_type_mpi)  /* sl_macro */
 #define data13_size_mpi                            (sl_data13_size_mpi)  /* sl_macro */

 #define data13_idx                                 13  /* sl_macro */

 #define data13_n                                   1  /* sl_macro */
 #define data13_byte                                (sl_data13_byte)  /* sl_macro */
 #define data13_ptr(e)                              (e)->data13  /* sl_macro */

 #ifdef sl_data13_flex
 # define data13_byte_flex                          (sl_data13_byte)  /* sl_macro */
 #else
 # define data13_byte_flex                          0
 #endif

 #ifdef sl_data13_weight
 # define data13_weight                             1  /* sl_macro */
 #else
 # define data13_weight                             0
 #endif

 /* commands for regular use */
 #define data13_assign(src, dst)                    ((dst)->data13 = (src)->data13)  /* sl_macro */
 #define data13_assign_at(src, sat, dst)            ((dst)->data13 = &(src)->data13[(sat) * data13_size_c])  /* sl_macro */
 #define data13_null(e)                             ((e)->data13 = NULL)  /* sl_macro */
 #define data13_inc(e)                              ((e)->data13 += data13_size_c)  /* sl_macro */
 #define data13_dec(e)                              ((e)->data13 -= data13_size_c)  /* sl_macro */
 #define data13_add(e, n)                           ((e)->data13 += (n) * data13_size_c)  /* sl_macro */
 #define data13_sub(e, n)                           ((e)->data13 -= (n) * data13_size_c)  /* sl_macro */

 #define data13_copy(src, dst)                      sl_data13_copy((src)->data13, (dst)->data13)  /* sl_macro */
 #define data13_ncopy(src, dst, n)                  sl_data13_ncopy((src)->data13, (dst)->data13, n)  /* sl_macro */
 #define data13_nmove(src, dst, n)                  sl_data13_nmove((src)->data13, (dst)->data13, n)  /* sl_macro */

 #define data13_copy_at(src, sat, dst, dat)         sl_data13_copy(&(src)->data13[(sat) * data13_size_c], &(dst)->data13[(dat) * data13_size_c])  /* sl_macro */
 #define data13_ncopy_at(src, sat, dst, dat, n)     sl_data13_ncopy(&(src)->data13[(sat) * data13_size_c], &(dst)->data13[(dat) * data13_size_c], n)  /* sl_macro */
 #define data13_nmove_at(src, sat, dst, dat, n)     sl_data13_nmove(&(src)->data13[(sat) * data13_size_c], &(dst)->data13[(dat) * data13_size_c], n)  /* sl_macro */

 #define data13_xchange(e0, e1, t)                  (data13_copy(e0, t), data13_copy(e1, e0), data13_copy(t, e1))  /* sl_macro */
 #define data13_xchange_at(e0, at0, e1, at1, t)     (data13_copy_at(e0, at0, t, 0), data13_copy_at(e1, at1, e0, at0), data13_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data13_assign(src, dst)                 , data13_assign(src, dst)  /* sl_macro */
 #define cc_data13_assign_at(src, sat, dst)         , data13_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data13_null(e)                          , data13_null(e)  /* sl_macro */
 #define cc_data13_inc(e)                           , data13_inc(e)  /* sl_macro */
 #define cc_data13_dec(e)                           , data13_dec(e)  /* sl_macro */
 #define cc_data13_add(e, n)                        , data13_add(e, n)  /* sl_macro */
 #define cc_data13_sub(e, n)                        , data13_sub(e, n)  /* sl_macro */
 #define cc_data13_copy(src, dst)                   , data13_copy(src, dst)  /* sl_macro */
 #define cc_data13_ncopy(src, dst, n)               , data13_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data13_nmove(src, dst, n)               , data13_nmove(src, dst, n)  /* sl_macro */
 #define cc_data13_copy_at(src, sat, dst, dat)      , data13_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data13_ncopy_at(src, sat, dst, dat, n)  , data13_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data13_nmove_at(src, sat, dst, dat, n)  , data13_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data13_xchange(e0, e1, t)               , data13_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data13_xchange_at(e0, at0, e1, at1, t)  , data13_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data13_assign(src, dst)                 data13_assign(src, dst)
 #define cc_data13_assign_at(src, sat, dst)         data13_assign_at(src, sat, dst)
 #define cc_data13_null(e)                          data13_null(e)
 #define cc_data13_inc(e)                           data13_inc(e)
 #define cc_data13_dec(e)                           data13_dec(e)
 #define cc_data13_add(e, n)                        data13_add(e, n)
 #define cc_data13_sub(e, n)                        data13_sub(e, n)
 #define cc_data13_copy(src, dst)                   data13_copy(src, dst)
 #define cc_data13_ncopy(src, dst, n)               data13_ncopy(src, dst, n)
 #define cc_data13_nmove(src, dst, n)               data13_nmove(src, dst, n)
 #define cc_data13_copy_at(src, sat, dst, dat)      data13_copy_at(src, sat, dst, dat)
 #define cc_data13_ncopy_at(src, sat, dst, dat, n)  data13_ncopy_at(src, sat, dst, dat, n)
 #define cc_data13_nmove_at(src, sat, dst, dat, n)  data13_nmove_at(src, sat, dst, dat, n)
 #define cc_data13_xchange(e0, e1, t)               data13_xchange(e0, e1, t)
 #define cc_data13_xchange_at(e0, at0, e1, at1, t)  data13_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA13 */

 #define data13_n                                   0
 #define data13_byte                                0
/* #define data13_ptr(e)*/

 #define data13_byte_flex                           0
 #define data13_weight                              0

 /* commands for regular use */
 #define data13_assign(src, dst)                    Z_NOP()
 #define data13_assign_at(src, sat, dst)            Z_NOP()
 #define data13_null(e)                             Z_NOP()
 #define data13_inc(e)                              Z_NOP()
 #define data13_dec(e)                              Z_NOP()
 #define data13_add(e, n)                           Z_NOP()
 #define data13_sub(e, n)                           Z_NOP()
 #define data13_copy(src, dst)                      Z_NOP()
 #define data13_ncopy(src, dst, n)                  Z_NOP()
 #define data13_nmove(src, dst, n)                  Z_NOP()
 #define data13_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data13_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data13_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data13_xchange(e0, e1, t)                  Z_NOP()
 #define data13_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data13_assign(src, dst)
 #define cc_data13_assign_at(src, sat, dst)
 #define cc_data13_null(e)
 #define cc_data13_inc(e)
 #define cc_data13_dec(e)
 #define cc_data13_add(e, n)
 #define cc_data13_sub(e, n)
 #define cc_data13_copy(src, dst)
 #define cc_data13_ncopy(src, dst, n)
 #define cc_data13_nmove(src, dst, n)
 #define cc_data13_copy_at(src, sat, dst, dat)
 #define cc_data13_ncopy_at(src, sat, dst, dat, n)
 #define cc_data13_nmove_at(src, sat, dst, dat, n)
 #define cc_data13_xchange(e0, e1, t)
 #define cc_data13_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA13 */

#define data13_cm                                   SLCM_DATA13  /* sl_macro */


/* sl_macro SL_DATA14 SL_DATA14_IGNORE sl_data14_type_c sl_data14_size_c sl_data14_type_mpi sl_data14_size_mpi sl_data14_memcpy sl_data14_weight sl_data14_flex */

#ifdef SL_DATA14

 #define sl_data14_byte                             ((slint_t) (sl_data14_size_c) * sizeof(sl_data14_type_c))  /* sl_macro */

 #ifndef sl_data14_copy
  #if sl_data14_size_c <= 9 && !defined(sl_data14_memcpy)
   #if sl_data14_size_c == 1
    #define sl_data14_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data14_size_c == 2
    #define sl_data14_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data14_size_c == 3
    #define sl_data14_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data14_size_c == 4
    #define sl_data14_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data14_size_c == 5
    #define sl_data14_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data14_size_c == 6
    #define sl_data14_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data14_size_c == 7
    #define sl_data14_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data14_size_c == 8
    #define sl_data14_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data14_size_c == 9
    #define sl_data14_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data14_copy(src, dst)                 memcpy(dst, src, sl_data14_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data14_ncopy
  #define sl_data14_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data14_byte)  /* sl_macro */
 #endif
 #ifndef sl_data14_nmove
  #define sl_data14_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data14_byte)  /* sl_macro */
 #endif

 #define data14_type_c                              sl_data14_type_c  /* sl_macro */
 #define data14_size_c                              (sl_data14_size_c)  /* sl_macro */
 #define data14_type_mpi                            (sl_data14_type_mpi)  /* sl_macro */
 #define data14_size_mpi                            (sl_data14_size_mpi)  /* sl_macro */

 #define data14_idx                                 14  /* sl_macro */

 #define data14_n                                   1  /* sl_macro */
 #define data14_byte                                (sl_data14_byte)  /* sl_macro */
 #define data14_ptr(e)                              (e)->data14  /* sl_macro */

 #ifdef sl_data14_flex
 # define data14_byte_flex                          (sl_data14_byte)  /* sl_macro */
 #else
 # define data14_byte_flex                          0
 #endif

 #ifdef sl_data14_weight
 # define data14_weight                             1  /* sl_macro */
 #else
 # define data14_weight                             0
 #endif

 /* commands for regular use */
 #define data14_assign(src, dst)                    ((dst)->data14 = (src)->data14)  /* sl_macro */
 #define data14_assign_at(src, sat, dst)            ((dst)->data14 = &(src)->data14[(sat) * data14_size_c])  /* sl_macro */
 #define data14_null(e)                             ((e)->data14 = NULL)  /* sl_macro */
 #define data14_inc(e)                              ((e)->data14 += data14_size_c)  /* sl_macro */
 #define data14_dec(e)                              ((e)->data14 -= data14_size_c)  /* sl_macro */
 #define data14_add(e, n)                           ((e)->data14 += (n) * data14_size_c)  /* sl_macro */
 #define data14_sub(e, n)                           ((e)->data14 -= (n) * data14_size_c)  /* sl_macro */

 #define data14_copy(src, dst)                      sl_data14_copy((src)->data14, (dst)->data14)  /* sl_macro */
 #define data14_ncopy(src, dst, n)                  sl_data14_ncopy((src)->data14, (dst)->data14, n)  /* sl_macro */
 #define data14_nmove(src, dst, n)                  sl_data14_nmove((src)->data14, (dst)->data14, n)  /* sl_macro */

 #define data14_copy_at(src, sat, dst, dat)         sl_data14_copy(&(src)->data14[(sat) * data14_size_c], &(dst)->data14[(dat) * data14_size_c])  /* sl_macro */
 #define data14_ncopy_at(src, sat, dst, dat, n)     sl_data14_ncopy(&(src)->data14[(sat) * data14_size_c], &(dst)->data14[(dat) * data14_size_c], n)  /* sl_macro */
 #define data14_nmove_at(src, sat, dst, dat, n)     sl_data14_nmove(&(src)->data14[(sat) * data14_size_c], &(dst)->data14[(dat) * data14_size_c], n)  /* sl_macro */

 #define data14_xchange(e0, e1, t)                  (data14_copy(e0, t), data14_copy(e1, e0), data14_copy(t, e1))  /* sl_macro */
 #define data14_xchange_at(e0, at0, e1, at1, t)     (data14_copy_at(e0, at0, t, 0), data14_copy_at(e1, at1, e0, at0), data14_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data14_assign(src, dst)                 , data14_assign(src, dst)  /* sl_macro */
 #define cc_data14_assign_at(src, sat, dst)         , data14_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data14_null(e)                          , data14_null(e)  /* sl_macro */
 #define cc_data14_inc(e)                           , data14_inc(e)  /* sl_macro */
 #define cc_data14_dec(e)                           , data14_dec(e)  /* sl_macro */
 #define cc_data14_add(e, n)                        , data14_add(e, n)  /* sl_macro */
 #define cc_data14_sub(e, n)                        , data14_sub(e, n)  /* sl_macro */
 #define cc_data14_copy(src, dst)                   , data14_copy(src, dst)  /* sl_macro */
 #define cc_data14_ncopy(src, dst, n)               , data14_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data14_nmove(src, dst, n)               , data14_nmove(src, dst, n)  /* sl_macro */
 #define cc_data14_copy_at(src, sat, dst, dat)      , data14_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data14_ncopy_at(src, sat, dst, dat, n)  , data14_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data14_nmove_at(src, sat, dst, dat, n)  , data14_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data14_xchange(e0, e1, t)               , data14_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data14_xchange_at(e0, at0, e1, at1, t)  , data14_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data14_assign(src, dst)                 data14_assign(src, dst)
 #define cc_data14_assign_at(src, sat, dst)         data14_assign_at(src, sat, dst)
 #define cc_data14_null(e)                          data14_null(e)
 #define cc_data14_inc(e)                           data14_inc(e)
 #define cc_data14_dec(e)                           data14_dec(e)
 #define cc_data14_add(e, n)                        data14_add(e, n)
 #define cc_data14_sub(e, n)                        data14_sub(e, n)
 #define cc_data14_copy(src, dst)                   data14_copy(src, dst)
 #define cc_data14_ncopy(src, dst, n)               data14_ncopy(src, dst, n)
 #define cc_data14_nmove(src, dst, n)               data14_nmove(src, dst, n)
 #define cc_data14_copy_at(src, sat, dst, dat)      data14_copy_at(src, sat, dst, dat)
 #define cc_data14_ncopy_at(src, sat, dst, dat, n)  data14_ncopy_at(src, sat, dst, dat, n)
 #define cc_data14_nmove_at(src, sat, dst, dat, n)  data14_nmove_at(src, sat, dst, dat, n)
 #define cc_data14_xchange(e0, e1, t)               data14_xchange(e0, e1, t)
 #define cc_data14_xchange_at(e0, at0, e1, at1, t)  data14_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA14 */

 #define data14_n                                   0
 #define data14_byte                                0
/* #define data14_ptr(e)*/

 #define data14_byte_flex                           0
 #define data14_weight                              0

 /* commands for regular use */
 #define data14_assign(src, dst)                    Z_NOP()
 #define data14_assign_at(src, sat, dst)            Z_NOP()
 #define data14_null(e)                             Z_NOP()
 #define data14_inc(e)                              Z_NOP()
 #define data14_dec(e)                              Z_NOP()
 #define data14_add(e, n)                           Z_NOP()
 #define data14_sub(e, n)                           Z_NOP()
 #define data14_copy(src, dst)                      Z_NOP()
 #define data14_ncopy(src, dst, n)                  Z_NOP()
 #define data14_nmove(src, dst, n)                  Z_NOP()
 #define data14_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data14_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data14_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data14_xchange(e0, e1, t)                  Z_NOP()
 #define data14_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data14_assign(src, dst)
 #define cc_data14_assign_at(src, sat, dst)
 #define cc_data14_null(e)
 #define cc_data14_inc(e)
 #define cc_data14_dec(e)
 #define cc_data14_add(e, n)
 #define cc_data14_sub(e, n)
 #define cc_data14_copy(src, dst)
 #define cc_data14_ncopy(src, dst, n)
 #define cc_data14_nmove(src, dst, n)
 #define cc_data14_copy_at(src, sat, dst, dat)
 #define cc_data14_ncopy_at(src, sat, dst, dat, n)
 #define cc_data14_nmove_at(src, sat, dst, dat, n)
 #define cc_data14_xchange(e0, e1, t)
 #define cc_data14_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA14 */

#define data14_cm                                   SLCM_DATA14  /* sl_macro */


/* sl_macro SL_DATA15 SL_DATA15_IGNORE sl_data15_type_c sl_data15_size_c sl_data15_type_mpi sl_data15_size_mpi sl_data15_memcpy sl_data15_weight sl_data15_flex */

#ifdef SL_DATA15

 #define sl_data15_byte                             ((slint_t) (sl_data15_size_c) * sizeof(sl_data15_type_c))  /* sl_macro */

 #ifndef sl_data15_copy
  #if sl_data15_size_c <= 9 && !defined(sl_data15_memcpy)
   #if sl_data15_size_c == 1
    #define sl_data15_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data15_size_c == 2
    #define sl_data15_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data15_size_c == 3
    #define sl_data15_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data15_size_c == 4
    #define sl_data15_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data15_size_c == 5
    #define sl_data15_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data15_size_c == 6
    #define sl_data15_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data15_size_c == 7
    #define sl_data15_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data15_size_c == 8
    #define sl_data15_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data15_size_c == 9
    #define sl_data15_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data15_copy(src, dst)                 memcpy(dst, src, sl_data15_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data15_ncopy
  #define sl_data15_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data15_byte)  /* sl_macro */
 #endif
 #ifndef sl_data15_nmove
  #define sl_data15_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data15_byte)  /* sl_macro */
 #endif

 #define data15_type_c                              sl_data15_type_c  /* sl_macro */
 #define data15_size_c                              (sl_data15_size_c)  /* sl_macro */
 #define data15_type_mpi                            (sl_data15_type_mpi)  /* sl_macro */
 #define data15_size_mpi                            (sl_data15_size_mpi)  /* sl_macro */

 #define data15_idx                                 15  /* sl_macro */

 #define data15_n                                   1  /* sl_macro */
 #define data15_byte                                (sl_data15_byte)  /* sl_macro */
 #define data15_ptr(e)                              (e)->data15  /* sl_macro */

 #ifdef sl_data15_flex
 # define data15_byte_flex                          (sl_data15_byte)  /* sl_macro */
 #else
 # define data15_byte_flex                          0
 #endif

 #ifdef sl_data15_weight
 # define data15_weight                             1  /* sl_macro */
 #else
 # define data15_weight                             0
 #endif

 /* commands for regular use */
 #define data15_assign(src, dst)                    ((dst)->data15 = (src)->data15)  /* sl_macro */
 #define data15_assign_at(src, sat, dst)            ((dst)->data15 = &(src)->data15[(sat) * data15_size_c])  /* sl_macro */
 #define data15_null(e)                             ((e)->data15 = NULL)  /* sl_macro */
 #define data15_inc(e)                              ((e)->data15 += data15_size_c)  /* sl_macro */
 #define data15_dec(e)                              ((e)->data15 -= data15_size_c)  /* sl_macro */
 #define data15_add(e, n)                           ((e)->data15 += (n) * data15_size_c)  /* sl_macro */
 #define data15_sub(e, n)                           ((e)->data15 -= (n) * data15_size_c)  /* sl_macro */

 #define data15_copy(src, dst)                      sl_data15_copy((src)->data15, (dst)->data15)  /* sl_macro */
 #define data15_ncopy(src, dst, n)                  sl_data15_ncopy((src)->data15, (dst)->data15, n)  /* sl_macro */
 #define data15_nmove(src, dst, n)                  sl_data15_nmove((src)->data15, (dst)->data15, n)  /* sl_macro */

 #define data15_copy_at(src, sat, dst, dat)         sl_data15_copy(&(src)->data15[(sat) * data15_size_c], &(dst)->data15[(dat) * data15_size_c])  /* sl_macro */
 #define data15_ncopy_at(src, sat, dst, dat, n)     sl_data15_ncopy(&(src)->data15[(sat) * data15_size_c], &(dst)->data15[(dat) * data15_size_c], n)  /* sl_macro */
 #define data15_nmove_at(src, sat, dst, dat, n)     sl_data15_nmove(&(src)->data15[(sat) * data15_size_c], &(dst)->data15[(dat) * data15_size_c], n)  /* sl_macro */

 #define data15_xchange(e0, e1, t)                  (data15_copy(e0, t), data15_copy(e1, e0), data15_copy(t, e1))  /* sl_macro */
 #define data15_xchange_at(e0, at0, e1, at1, t)     (data15_copy_at(e0, at0, t, 0), data15_copy_at(e1, at1, e0, at0), data15_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data15_assign(src, dst)                 , data15_assign(src, dst)  /* sl_macro */
 #define cc_data15_assign_at(src, sat, dst)         , data15_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data15_null(e)                          , data15_null(e)  /* sl_macro */
 #define cc_data15_inc(e)                           , data15_inc(e)  /* sl_macro */
 #define cc_data15_dec(e)                           , data15_dec(e)  /* sl_macro */
 #define cc_data15_add(e, n)                        , data15_add(e, n)  /* sl_macro */
 #define cc_data15_sub(e, n)                        , data15_sub(e, n)  /* sl_macro */
 #define cc_data15_copy(src, dst)                   , data15_copy(src, dst)  /* sl_macro */
 #define cc_data15_ncopy(src, dst, n)               , data15_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data15_nmove(src, dst, n)               , data15_nmove(src, dst, n)  /* sl_macro */
 #define cc_data15_copy_at(src, sat, dst, dat)      , data15_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data15_ncopy_at(src, sat, dst, dat, n)  , data15_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data15_nmove_at(src, sat, dst, dat, n)  , data15_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data15_xchange(e0, e1, t)               , data15_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data15_xchange_at(e0, at0, e1, at1, t)  , data15_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data15_assign(src, dst)                 data15_assign(src, dst)
 #define cc_data15_assign_at(src, sat, dst)         data15_assign_at(src, sat, dst)
 #define cc_data15_null(e)                          data15_null(e)
 #define cc_data15_inc(e)                           data15_inc(e)
 #define cc_data15_dec(e)                           data15_dec(e)
 #define cc_data15_add(e, n)                        data15_add(e, n)
 #define cc_data15_sub(e, n)                        data15_sub(e, n)
 #define cc_data15_copy(src, dst)                   data15_copy(src, dst)
 #define cc_data15_ncopy(src, dst, n)               data15_ncopy(src, dst, n)
 #define cc_data15_nmove(src, dst, n)               data15_nmove(src, dst, n)
 #define cc_data15_copy_at(src, sat, dst, dat)      data15_copy_at(src, sat, dst, dat)
 #define cc_data15_ncopy_at(src, sat, dst, dat, n)  data15_ncopy_at(src, sat, dst, dat, n)
 #define cc_data15_nmove_at(src, sat, dst, dat, n)  data15_nmove_at(src, sat, dst, dat, n)
 #define cc_data15_xchange(e0, e1, t)               data15_xchange(e0, e1, t)
 #define cc_data15_xchange_at(e0, at0, e1, at1, t)  data15_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA15 */

 #define data15_n                                   0
 #define data15_byte                                0
/* #define data15_ptr(e)*/

 #define data15_byte_flex                           0
 #define data15_weight                              0

 /* commands for regular use */
 #define data15_assign(src, dst)                    Z_NOP()
 #define data15_assign_at(src, sat, dst)            Z_NOP()
 #define data15_null(e)                             Z_NOP()
 #define data15_inc(e)                              Z_NOP()
 #define data15_dec(e)                              Z_NOP()
 #define data15_add(e, n)                           Z_NOP()
 #define data15_sub(e, n)                           Z_NOP()
 #define data15_copy(src, dst)                      Z_NOP()
 #define data15_ncopy(src, dst, n)                  Z_NOP()
 #define data15_nmove(src, dst, n)                  Z_NOP()
 #define data15_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data15_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data15_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data15_xchange(e0, e1, t)                  Z_NOP()
 #define data15_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data15_assign(src, dst)
 #define cc_data15_assign_at(src, sat, dst)
 #define cc_data15_null(e)
 #define cc_data15_inc(e)
 #define cc_data15_dec(e)
 #define cc_data15_add(e, n)
 #define cc_data15_sub(e, n)
 #define cc_data15_copy(src, dst)
 #define cc_data15_ncopy(src, dst, n)
 #define cc_data15_nmove(src, dst, n)
 #define cc_data15_copy_at(src, sat, dst, dat)
 #define cc_data15_ncopy_at(src, sat, dst, dat, n)
 #define cc_data15_nmove_at(src, sat, dst, dat, n)
 #define cc_data15_xchange(e0, e1, t)
 #define cc_data15_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA15 */

#define data15_cm                                   SLCM_DATA15  /* sl_macro */


/* sl_macro SL_DATA16 SL_DATA16_IGNORE sl_data16_type_c sl_data16_size_c sl_data16_type_mpi sl_data16_size_mpi sl_data16_memcpy sl_data16_weight sl_data16_flex */

#ifdef SL_DATA16

 #define sl_data16_byte                             ((slint_t) (sl_data16_size_c) * sizeof(sl_data16_type_c))  /* sl_macro */

 #ifndef sl_data16_copy
  #if sl_data16_size_c <= 9 && !defined(sl_data16_memcpy)
   #if sl_data16_size_c == 1
    #define sl_data16_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data16_size_c == 2
    #define sl_data16_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data16_size_c == 3
    #define sl_data16_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data16_size_c == 4
    #define sl_data16_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data16_size_c == 5
    #define sl_data16_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data16_size_c == 6
    #define sl_data16_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data16_size_c == 7
    #define sl_data16_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data16_size_c == 8
    #define sl_data16_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data16_size_c == 9
    #define sl_data16_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data16_copy(src, dst)                 memcpy(dst, src, sl_data16_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data16_ncopy
  #define sl_data16_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data16_byte)  /* sl_macro */
 #endif
 #ifndef sl_data16_nmove
  #define sl_data16_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data16_byte)  /* sl_macro */
 #endif

 #define data16_type_c                              sl_data16_type_c  /* sl_macro */
 #define data16_size_c                              (sl_data16_size_c)  /* sl_macro */
 #define data16_type_mpi                            (sl_data16_type_mpi)  /* sl_macro */
 #define data16_size_mpi                            (sl_data16_size_mpi)  /* sl_macro */

 #define data16_idx                                 16  /* sl_macro */

 #define data16_n                                   1  /* sl_macro */
 #define data16_byte                                (sl_data16_byte)  /* sl_macro */
 #define data16_ptr(e)                              (e)->data16  /* sl_macro */

 #ifdef sl_data16_flex
 # define data16_byte_flex                          (sl_data16_byte)  /* sl_macro */
 #else
 # define data16_byte_flex                          0
 #endif

 #ifdef sl_data16_weight
 # define data16_weight                             1  /* sl_macro */
 #else
 # define data16_weight                             0
 #endif

 /* commands for regular use */
 #define data16_assign(src, dst)                    ((dst)->data16 = (src)->data16)  /* sl_macro */
 #define data16_assign_at(src, sat, dst)            ((dst)->data16 = &(src)->data16[(sat) * data16_size_c])  /* sl_macro */
 #define data16_null(e)                             ((e)->data16 = NULL)  /* sl_macro */
 #define data16_inc(e)                              ((e)->data16 += data16_size_c)  /* sl_macro */
 #define data16_dec(e)                              ((e)->data16 -= data16_size_c)  /* sl_macro */
 #define data16_add(e, n)                           ((e)->data16 += (n) * data16_size_c)  /* sl_macro */
 #define data16_sub(e, n)                           ((e)->data16 -= (n) * data16_size_c)  /* sl_macro */

 #define data16_copy(src, dst)                      sl_data16_copy((src)->data16, (dst)->data16)  /* sl_macro */
 #define data16_ncopy(src, dst, n)                  sl_data16_ncopy((src)->data16, (dst)->data16, n)  /* sl_macro */
 #define data16_nmove(src, dst, n)                  sl_data16_nmove((src)->data16, (dst)->data16, n)  /* sl_macro */

 #define data16_copy_at(src, sat, dst, dat)         sl_data16_copy(&(src)->data16[(sat) * data16_size_c], &(dst)->data16[(dat) * data16_size_c])  /* sl_macro */
 #define data16_ncopy_at(src, sat, dst, dat, n)     sl_data16_ncopy(&(src)->data16[(sat) * data16_size_c], &(dst)->data16[(dat) * data16_size_c], n)  /* sl_macro */
 #define data16_nmove_at(src, sat, dst, dat, n)     sl_data16_nmove(&(src)->data16[(sat) * data16_size_c], &(dst)->data16[(dat) * data16_size_c], n)  /* sl_macro */

 #define data16_xchange(e0, e1, t)                  (data16_copy(e0, t), data16_copy(e1, e0), data16_copy(t, e1))  /* sl_macro */
 #define data16_xchange_at(e0, at0, e1, at1, t)     (data16_copy_at(e0, at0, t, 0), data16_copy_at(e1, at1, e0, at0), data16_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data16_assign(src, dst)                 , data16_assign(src, dst)  /* sl_macro */
 #define cc_data16_assign_at(src, sat, dst)         , data16_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data16_null(e)                          , data16_null(e)  /* sl_macro */
 #define cc_data16_inc(e)                           , data16_inc(e)  /* sl_macro */
 #define cc_data16_dec(e)                           , data16_dec(e)  /* sl_macro */
 #define cc_data16_add(e, n)                        , data16_add(e, n)  /* sl_macro */
 #define cc_data16_sub(e, n)                        , data16_sub(e, n)  /* sl_macro */
 #define cc_data16_copy(src, dst)                   , data16_copy(src, dst)  /* sl_macro */
 #define cc_data16_ncopy(src, dst, n)               , data16_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data16_nmove(src, dst, n)               , data16_nmove(src, dst, n)  /* sl_macro */
 #define cc_data16_copy_at(src, sat, dst, dat)      , data16_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data16_ncopy_at(src, sat, dst, dat, n)  , data16_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data16_nmove_at(src, sat, dst, dat, n)  , data16_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data16_xchange(e0, e1, t)               , data16_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data16_xchange_at(e0, at0, e1, at1, t)  , data16_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data16_assign(src, dst)                 data16_assign(src, dst)
 #define cc_data16_assign_at(src, sat, dst)         data16_assign_at(src, sat, dst)
 #define cc_data16_null(e)                          data16_null(e)
 #define cc_data16_inc(e)                           data16_inc(e)
 #define cc_data16_dec(e)                           data16_dec(e)
 #define cc_data16_add(e, n)                        data16_add(e, n)
 #define cc_data16_sub(e, n)                        data16_sub(e, n)
 #define cc_data16_copy(src, dst)                   data16_copy(src, dst)
 #define cc_data16_ncopy(src, dst, n)               data16_ncopy(src, dst, n)
 #define cc_data16_nmove(src, dst, n)               data16_nmove(src, dst, n)
 #define cc_data16_copy_at(src, sat, dst, dat)      data16_copy_at(src, sat, dst, dat)
 #define cc_data16_ncopy_at(src, sat, dst, dat, n)  data16_ncopy_at(src, sat, dst, dat, n)
 #define cc_data16_nmove_at(src, sat, dst, dat, n)  data16_nmove_at(src, sat, dst, dat, n)
 #define cc_data16_xchange(e0, e1, t)               data16_xchange(e0, e1, t)
 #define cc_data16_xchange_at(e0, at0, e1, at1, t)  data16_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA16 */

 #define data16_n                                   0
 #define data16_byte                                0
/* #define data16_ptr(e)*/

 #define data16_byte_flex                           0
 #define data16_weight                              0

 /* commands for regular use */
 #define data16_assign(src, dst)                    Z_NOP()
 #define data16_assign_at(src, sat, dst)            Z_NOP()
 #define data16_null(e)                             Z_NOP()
 #define data16_inc(e)                              Z_NOP()
 #define data16_dec(e)                              Z_NOP()
 #define data16_add(e, n)                           Z_NOP()
 #define data16_sub(e, n)                           Z_NOP()
 #define data16_copy(src, dst)                      Z_NOP()
 #define data16_ncopy(src, dst, n)                  Z_NOP()
 #define data16_nmove(src, dst, n)                  Z_NOP()
 #define data16_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data16_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data16_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data16_xchange(e0, e1, t)                  Z_NOP()
 #define data16_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data16_assign(src, dst)
 #define cc_data16_assign_at(src, sat, dst)
 #define cc_data16_null(e)
 #define cc_data16_inc(e)
 #define cc_data16_dec(e)
 #define cc_data16_add(e, n)
 #define cc_data16_sub(e, n)
 #define cc_data16_copy(src, dst)
 #define cc_data16_ncopy(src, dst, n)
 #define cc_data16_nmove(src, dst, n)
 #define cc_data16_copy_at(src, sat, dst, dat)
 #define cc_data16_ncopy_at(src, sat, dst, dat, n)
 #define cc_data16_nmove_at(src, sat, dst, dat, n)
 #define cc_data16_xchange(e0, e1, t)
 #define cc_data16_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA16 */

#define data16_cm                                   SLCM_DATA16  /* sl_macro */


/* sl_macro SL_DATA17 SL_DATA17_IGNORE sl_data17_type_c sl_data17_size_c sl_data17_type_mpi sl_data17_size_mpi sl_data17_memcpy sl_data17_weight sl_data17_flex */

#ifdef SL_DATA17

 #define sl_data17_byte                             ((slint_t) (sl_data17_size_c) * sizeof(sl_data17_type_c))  /* sl_macro */

 #ifndef sl_data17_copy
  #if sl_data17_size_c <= 9 && !defined(sl_data17_memcpy)
   #if sl_data17_size_c == 1
    #define sl_data17_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data17_size_c == 2
    #define sl_data17_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data17_size_c == 3
    #define sl_data17_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data17_size_c == 4
    #define sl_data17_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data17_size_c == 5
    #define sl_data17_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data17_size_c == 6
    #define sl_data17_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data17_size_c == 7
    #define sl_data17_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data17_size_c == 8
    #define sl_data17_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data17_size_c == 9
    #define sl_data17_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data17_copy(src, dst)                 memcpy(dst, src, sl_data17_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data17_ncopy
  #define sl_data17_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data17_byte)  /* sl_macro */
 #endif
 #ifndef sl_data17_nmove
  #define sl_data17_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data17_byte)  /* sl_macro */
 #endif

 #define data17_type_c                              sl_data17_type_c  /* sl_macro */
 #define data17_size_c                              (sl_data17_size_c)  /* sl_macro */
 #define data17_type_mpi                            (sl_data17_type_mpi)  /* sl_macro */
 #define data17_size_mpi                            (sl_data17_size_mpi)  /* sl_macro */

 #define data17_idx                                 17  /* sl_macro */

 #define data17_n                                   1  /* sl_macro */
 #define data17_byte                                (sl_data17_byte)  /* sl_macro */
 #define data17_ptr(e)                              (e)->data17  /* sl_macro */

 #ifdef sl_data17_flex
 # define data17_byte_flex                          (sl_data17_byte)  /* sl_macro */
 #else
 # define data17_byte_flex                          0
 #endif

 #ifdef sl_data17_weight
 # define data17_weight                             1  /* sl_macro */
 #else
 # define data17_weight                             0
 #endif

 /* commands for regular use */
 #define data17_assign(src, dst)                    ((dst)->data17 = (src)->data17)  /* sl_macro */
 #define data17_assign_at(src, sat, dst)            ((dst)->data17 = &(src)->data17[(sat) * data17_size_c])  /* sl_macro */
 #define data17_null(e)                             ((e)->data17 = NULL)  /* sl_macro */
 #define data17_inc(e)                              ((e)->data17 += data17_size_c)  /* sl_macro */
 #define data17_dec(e)                              ((e)->data17 -= data17_size_c)  /* sl_macro */
 #define data17_add(e, n)                           ((e)->data17 += (n) * data17_size_c)  /* sl_macro */
 #define data17_sub(e, n)                           ((e)->data17 -= (n) * data17_size_c)  /* sl_macro */

 #define data17_copy(src, dst)                      sl_data17_copy((src)->data17, (dst)->data17)  /* sl_macro */
 #define data17_ncopy(src, dst, n)                  sl_data17_ncopy((src)->data17, (dst)->data17, n)  /* sl_macro */
 #define data17_nmove(src, dst, n)                  sl_data17_nmove((src)->data17, (dst)->data17, n)  /* sl_macro */

 #define data17_copy_at(src, sat, dst, dat)         sl_data17_copy(&(src)->data17[(sat) * data17_size_c], &(dst)->data17[(dat) * data17_size_c])  /* sl_macro */
 #define data17_ncopy_at(src, sat, dst, dat, n)     sl_data17_ncopy(&(src)->data17[(sat) * data17_size_c], &(dst)->data17[(dat) * data17_size_c], n)  /* sl_macro */
 #define data17_nmove_at(src, sat, dst, dat, n)     sl_data17_nmove(&(src)->data17[(sat) * data17_size_c], &(dst)->data17[(dat) * data17_size_c], n)  /* sl_macro */

 #define data17_xchange(e0, e1, t)                  (data17_copy(e0, t), data17_copy(e1, e0), data17_copy(t, e1))  /* sl_macro */
 #define data17_xchange_at(e0, at0, e1, at1, t)     (data17_copy_at(e0, at0, t, 0), data17_copy_at(e1, at1, e0, at0), data17_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data17_assign(src, dst)                 , data17_assign(src, dst)  /* sl_macro */
 #define cc_data17_assign_at(src, sat, dst)         , data17_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data17_null(e)                          , data17_null(e)  /* sl_macro */
 #define cc_data17_inc(e)                           , data17_inc(e)  /* sl_macro */
 #define cc_data17_dec(e)                           , data17_dec(e)  /* sl_macro */
 #define cc_data17_add(e, n)                        , data17_add(e, n)  /* sl_macro */
 #define cc_data17_sub(e, n)                        , data17_sub(e, n)  /* sl_macro */
 #define cc_data17_copy(src, dst)                   , data17_copy(src, dst)  /* sl_macro */
 #define cc_data17_ncopy(src, dst, n)               , data17_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data17_nmove(src, dst, n)               , data17_nmove(src, dst, n)  /* sl_macro */
 #define cc_data17_copy_at(src, sat, dst, dat)      , data17_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data17_ncopy_at(src, sat, dst, dat, n)  , data17_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data17_nmove_at(src, sat, dst, dat, n)  , data17_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data17_xchange(e0, e1, t)               , data17_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data17_xchange_at(e0, at0, e1, at1, t)  , data17_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data17_assign(src, dst)                 data17_assign(src, dst)
 #define cc_data17_assign_at(src, sat, dst)         data17_assign_at(src, sat, dst)
 #define cc_data17_null(e)                          data17_null(e)
 #define cc_data17_inc(e)                           data17_inc(e)
 #define cc_data17_dec(e)                           data17_dec(e)
 #define cc_data17_add(e, n)                        data17_add(e, n)
 #define cc_data17_sub(e, n)                        data17_sub(e, n)
 #define cc_data17_copy(src, dst)                   data17_copy(src, dst)
 #define cc_data17_ncopy(src, dst, n)               data17_ncopy(src, dst, n)
 #define cc_data17_nmove(src, dst, n)               data17_nmove(src, dst, n)
 #define cc_data17_copy_at(src, sat, dst, dat)      data17_copy_at(src, sat, dst, dat)
 #define cc_data17_ncopy_at(src, sat, dst, dat, n)  data17_ncopy_at(src, sat, dst, dat, n)
 #define cc_data17_nmove_at(src, sat, dst, dat, n)  data17_nmove_at(src, sat, dst, dat, n)
 #define cc_data17_xchange(e0, e1, t)               data17_xchange(e0, e1, t)
 #define cc_data17_xchange_at(e0, at0, e1, at1, t)  data17_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA17 */

 #define data17_n                                   0
 #define data17_byte                                0
/* #define data17_ptr(e)*/

 #define data17_byte_flex                           0
 #define data17_weight                              0

 /* commands for regular use */
 #define data17_assign(src, dst)                    Z_NOP()
 #define data17_assign_at(src, sat, dst)            Z_NOP()
 #define data17_null(e)                             Z_NOP()
 #define data17_inc(e)                              Z_NOP()
 #define data17_dec(e)                              Z_NOP()
 #define data17_add(e, n)                           Z_NOP()
 #define data17_sub(e, n)                           Z_NOP()
 #define data17_copy(src, dst)                      Z_NOP()
 #define data17_ncopy(src, dst, n)                  Z_NOP()
 #define data17_nmove(src, dst, n)                  Z_NOP()
 #define data17_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data17_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data17_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data17_xchange(e0, e1, t)                  Z_NOP()
 #define data17_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data17_assign(src, dst)
 #define cc_data17_assign_at(src, sat, dst)
 #define cc_data17_null(e)
 #define cc_data17_inc(e)
 #define cc_data17_dec(e)
 #define cc_data17_add(e, n)
 #define cc_data17_sub(e, n)
 #define cc_data17_copy(src, dst)
 #define cc_data17_ncopy(src, dst, n)
 #define cc_data17_nmove(src, dst, n)
 #define cc_data17_copy_at(src, sat, dst, dat)
 #define cc_data17_ncopy_at(src, sat, dst, dat, n)
 #define cc_data17_nmove_at(src, sat, dst, dat, n)
 #define cc_data17_xchange(e0, e1, t)
 #define cc_data17_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA17 */

#define data17_cm                                   SLCM_DATA17  /* sl_macro */


/* sl_macro SL_DATA18 SL_DATA18_IGNORE sl_data18_type_c sl_data18_size_c sl_data18_type_mpi sl_data18_size_mpi sl_data18_memcpy sl_data18_weight sl_data18_flex */

#ifdef SL_DATA18

 #define sl_data18_byte                             ((slint_t) (sl_data18_size_c) * sizeof(sl_data18_type_c))  /* sl_macro */

 #ifndef sl_data18_copy
  #if sl_data18_size_c <= 9 && !defined(sl_data18_memcpy)
   #if sl_data18_size_c == 1
    #define sl_data18_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data18_size_c == 2
    #define sl_data18_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data18_size_c == 3
    #define sl_data18_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data18_size_c == 4
    #define sl_data18_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data18_size_c == 5
    #define sl_data18_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data18_size_c == 6
    #define sl_data18_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data18_size_c == 7
    #define sl_data18_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data18_size_c == 8
    #define sl_data18_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data18_size_c == 9
    #define sl_data18_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data18_copy(src, dst)                 memcpy(dst, src, sl_data18_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data18_ncopy
  #define sl_data18_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data18_byte)  /* sl_macro */
 #endif
 #ifndef sl_data18_nmove
  #define sl_data18_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data18_byte)  /* sl_macro */
 #endif

 #define data18_type_c                              sl_data18_type_c  /* sl_macro */
 #define data18_size_c                              (sl_data18_size_c)  /* sl_macro */
 #define data18_type_mpi                            (sl_data18_type_mpi)  /* sl_macro */
 #define data18_size_mpi                            (sl_data18_size_mpi)  /* sl_macro */

 #define data18_idx                                 18  /* sl_macro */

 #define data18_n                                   1  /* sl_macro */
 #define data18_byte                                (sl_data18_byte)  /* sl_macro */
 #define data18_ptr(e)                              (e)->data18  /* sl_macro */

 #ifdef sl_data18_flex
 # define data18_byte_flex                          (sl_data18_byte)  /* sl_macro */
 #else
 # define data18_byte_flex                          0
 #endif

 #ifdef sl_data18_weight
 # define data18_weight                             1  /* sl_macro */
 #else
 # define data18_weight                             0
 #endif

 /* commands for regular use */
 #define data18_assign(src, dst)                    ((dst)->data18 = (src)->data18)  /* sl_macro */
 #define data18_assign_at(src, sat, dst)            ((dst)->data18 = &(src)->data18[(sat) * data18_size_c])  /* sl_macro */
 #define data18_null(e)                             ((e)->data18 = NULL)  /* sl_macro */
 #define data18_inc(e)                              ((e)->data18 += data18_size_c)  /* sl_macro */
 #define data18_dec(e)                              ((e)->data18 -= data18_size_c)  /* sl_macro */
 #define data18_add(e, n)                           ((e)->data18 += (n) * data18_size_c)  /* sl_macro */
 #define data18_sub(e, n)                           ((e)->data18 -= (n) * data18_size_c)  /* sl_macro */

 #define data18_copy(src, dst)                      sl_data18_copy((src)->data18, (dst)->data18)  /* sl_macro */
 #define data18_ncopy(src, dst, n)                  sl_data18_ncopy((src)->data18, (dst)->data18, n)  /* sl_macro */
 #define data18_nmove(src, dst, n)                  sl_data18_nmove((src)->data18, (dst)->data18, n)  /* sl_macro */

 #define data18_copy_at(src, sat, dst, dat)         sl_data18_copy(&(src)->data18[(sat) * data18_size_c], &(dst)->data18[(dat) * data18_size_c])  /* sl_macro */
 #define data18_ncopy_at(src, sat, dst, dat, n)     sl_data18_ncopy(&(src)->data18[(sat) * data18_size_c], &(dst)->data18[(dat) * data18_size_c], n)  /* sl_macro */
 #define data18_nmove_at(src, sat, dst, dat, n)     sl_data18_nmove(&(src)->data18[(sat) * data18_size_c], &(dst)->data18[(dat) * data18_size_c], n)  /* sl_macro */

 #define data18_xchange(e0, e1, t)                  (data18_copy(e0, t), data18_copy(e1, e0), data18_copy(t, e1))  /* sl_macro */
 #define data18_xchange_at(e0, at0, e1, at1, t)     (data18_copy_at(e0, at0, t, 0), data18_copy_at(e1, at1, e0, at0), data18_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data18_assign(src, dst)                 , data18_assign(src, dst)  /* sl_macro */
 #define cc_data18_assign_at(src, sat, dst)         , data18_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data18_null(e)                          , data18_null(e)  /* sl_macro */
 #define cc_data18_inc(e)                           , data18_inc(e)  /* sl_macro */
 #define cc_data18_dec(e)                           , data18_dec(e)  /* sl_macro */
 #define cc_data18_add(e, n)                        , data18_add(e, n)  /* sl_macro */
 #define cc_data18_sub(e, n)                        , data18_sub(e, n)  /* sl_macro */
 #define cc_data18_copy(src, dst)                   , data18_copy(src, dst)  /* sl_macro */
 #define cc_data18_ncopy(src, dst, n)               , data18_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data18_nmove(src, dst, n)               , data18_nmove(src, dst, n)  /* sl_macro */
 #define cc_data18_copy_at(src, sat, dst, dat)      , data18_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data18_ncopy_at(src, sat, dst, dat, n)  , data18_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data18_nmove_at(src, sat, dst, dat, n)  , data18_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data18_xchange(e0, e1, t)               , data18_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data18_xchange_at(e0, at0, e1, at1, t)  , data18_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data18_assign(src, dst)                 data18_assign(src, dst)
 #define cc_data18_assign_at(src, sat, dst)         data18_assign_at(src, sat, dst)
 #define cc_data18_null(e)                          data18_null(e)
 #define cc_data18_inc(e)                           data18_inc(e)
 #define cc_data18_dec(e)                           data18_dec(e)
 #define cc_data18_add(e, n)                        data18_add(e, n)
 #define cc_data18_sub(e, n)                        data18_sub(e, n)
 #define cc_data18_copy(src, dst)                   data18_copy(src, dst)
 #define cc_data18_ncopy(src, dst, n)               data18_ncopy(src, dst, n)
 #define cc_data18_nmove(src, dst, n)               data18_nmove(src, dst, n)
 #define cc_data18_copy_at(src, sat, dst, dat)      data18_copy_at(src, sat, dst, dat)
 #define cc_data18_ncopy_at(src, sat, dst, dat, n)  data18_ncopy_at(src, sat, dst, dat, n)
 #define cc_data18_nmove_at(src, sat, dst, dat, n)  data18_nmove_at(src, sat, dst, dat, n)
 #define cc_data18_xchange(e0, e1, t)               data18_xchange(e0, e1, t)
 #define cc_data18_xchange_at(e0, at0, e1, at1, t)  data18_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA18 */

 #define data18_n                                   0
 #define data18_byte                                0
/* #define data18_ptr(e)*/

 #define data18_byte_flex                           0
 #define data18_weight                              0

 /* commands for regular use */
 #define data18_assign(src, dst)                    Z_NOP()
 #define data18_assign_at(src, sat, dst)            Z_NOP()
 #define data18_null(e)                             Z_NOP()
 #define data18_inc(e)                              Z_NOP()
 #define data18_dec(e)                              Z_NOP()
 #define data18_add(e, n)                           Z_NOP()
 #define data18_sub(e, n)                           Z_NOP()
 #define data18_copy(src, dst)                      Z_NOP()
 #define data18_ncopy(src, dst, n)                  Z_NOP()
 #define data18_nmove(src, dst, n)                  Z_NOP()
 #define data18_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data18_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data18_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data18_xchange(e0, e1, t)                  Z_NOP()
 #define data18_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data18_assign(src, dst)
 #define cc_data18_assign_at(src, sat, dst)
 #define cc_data18_null(e)
 #define cc_data18_inc(e)
 #define cc_data18_dec(e)
 #define cc_data18_add(e, n)
 #define cc_data18_sub(e, n)
 #define cc_data18_copy(src, dst)
 #define cc_data18_ncopy(src, dst, n)
 #define cc_data18_nmove(src, dst, n)
 #define cc_data18_copy_at(src, sat, dst, dat)
 #define cc_data18_ncopy_at(src, sat, dst, dat, n)
 #define cc_data18_nmove_at(src, sat, dst, dat, n)
 #define cc_data18_xchange(e0, e1, t)
 #define cc_data18_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA18 */

#define data18_cm                                   SLCM_DATA18  /* sl_macro */


/* sl_macro SL_DATA19 SL_DATA19_IGNORE sl_data19_type_c sl_data19_size_c sl_data19_type_mpi sl_data19_size_mpi sl_data19_memcpy sl_data19_weight sl_data19_flex */

#ifdef SL_DATA19

 #define sl_data19_byte                             ((slint_t) (sl_data19_size_c) * sizeof(sl_data19_type_c))  /* sl_macro */

 #ifndef sl_data19_copy
  #if sl_data19_size_c <= 9 && !defined(sl_data19_memcpy)
   #if sl_data19_size_c == 1
    #define sl_data19_copy(src, dst)                (SL_ARRAY1_COPY(src, dst))
   #elif sl_data19_size_c == 2
    #define sl_data19_copy(src, dst)                (SL_ARRAY2_COPY(src, dst))
   #elif sl_data19_size_c == 3
    #define sl_data19_copy(src, dst)                (SL_ARRAY3_COPY(src, dst))
   #elif sl_data19_size_c == 4
    #define sl_data19_copy(src, dst)                (SL_ARRAY4_COPY(src, dst))
   #elif sl_data19_size_c == 5
    #define sl_data19_copy(src, dst)                (SL_ARRAY5_COPY(src, dst))
   #elif sl_data19_size_c == 6
    #define sl_data19_copy(src, dst)                (SL_ARRAY6_COPY(src, dst))
   #elif sl_data19_size_c == 7
    #define sl_data19_copy(src, dst)                (SL_ARRAY7_COPY(src, dst))
   #elif sl_data19_size_c == 8
    #define sl_data19_copy(src, dst)                (SL_ARRAY8_COPY(src, dst))
   #elif sl_data19_size_c == 9
    #define sl_data19_copy(src, dst)                (SL_ARRAY9_COPY(src, dst))
   #endif
  #else
   #define sl_data19_copy(src, dst)                 memcpy(dst, src, sl_data19_byte)  /* sl_macro */
  #endif
 #endif
 #ifndef sl_data19_ncopy
  #define sl_data19_ncopy(src, dst, n)              memcpy(dst, src, (n) * sl_data19_byte)  /* sl_macro */
 #endif
 #ifndef sl_data19_nmove
  #define sl_data19_nmove(src, dst, n)              memmove(dst, src, (n) * sl_data19_byte)  /* sl_macro */
 #endif

 #define data19_type_c                              sl_data19_type_c  /* sl_macro */
 #define data19_size_c                              (sl_data19_size_c)  /* sl_macro */
 #define data19_type_mpi                            (sl_data19_type_mpi)  /* sl_macro */
 #define data19_size_mpi                            (sl_data19_size_mpi)  /* sl_macro */

 #define data19_idx                                 19  /* sl_macro */

 #define data19_n                                   1  /* sl_macro */
 #define data19_byte                                (sl_data19_byte)  /* sl_macro */
 #define data19_ptr(e)                              (e)->data19  /* sl_macro */

 #ifdef sl_data19_flex
 # define data19_byte_flex                          (sl_data19_byte)  /* sl_macro */
 #else
 # define data19_byte_flex                          0
 #endif

 #ifdef sl_data19_weight
 # define data19_weight                             1  /* sl_macro */
 #else
 # define data19_weight                             0
 #endif

 /* commands for regular use */
 #define data19_assign(src, dst)                    ((dst)->data19 = (src)->data19)  /* sl_macro */
 #define data19_assign_at(src, sat, dst)            ((dst)->data19 = &(src)->data19[(sat) * data19_size_c])  /* sl_macro */
 #define data19_null(e)                             ((e)->data19 = NULL)  /* sl_macro */
 #define data19_inc(e)                              ((e)->data19 += data19_size_c)  /* sl_macro */
 #define data19_dec(e)                              ((e)->data19 -= data19_size_c)  /* sl_macro */
 #define data19_add(e, n)                           ((e)->data19 += (n) * data19_size_c)  /* sl_macro */
 #define data19_sub(e, n)                           ((e)->data19 -= (n) * data19_size_c)  /* sl_macro */

 #define data19_copy(src, dst)                      sl_data19_copy((src)->data19, (dst)->data19)  /* sl_macro */
 #define data19_ncopy(src, dst, n)                  sl_data19_ncopy((src)->data19, (dst)->data19, n)  /* sl_macro */
 #define data19_nmove(src, dst, n)                  sl_data19_nmove((src)->data19, (dst)->data19, n)  /* sl_macro */

 #define data19_copy_at(src, sat, dst, dat)         sl_data19_copy(&(src)->data19[(sat) * data19_size_c], &(dst)->data19[(dat) * data19_size_c])  /* sl_macro */
 #define data19_ncopy_at(src, sat, dst, dat, n)     sl_data19_ncopy(&(src)->data19[(sat) * data19_size_c], &(dst)->data19[(dat) * data19_size_c], n)  /* sl_macro */
 #define data19_nmove_at(src, sat, dst, dat, n)     sl_data19_nmove(&(src)->data19[(sat) * data19_size_c], &(dst)->data19[(dat) * data19_size_c], n)  /* sl_macro */

 #define data19_xchange(e0, e1, t)                  (data19_copy(e0, t), data19_copy(e1, e0), data19_copy(t, e1))  /* sl_macro */
 #define data19_xchange_at(e0, at0, e1, at1, t)     (data19_copy_at(e0, at0, t, 0), data19_copy_at(e1, at1, e0, at0), data19_copy_at(t, 0, e1, at1))  /* sl_macro */

 /* chained command versions */
#ifdef SL_DATA
 #define cc_data19_assign(src, dst)                 , data19_assign(src, dst)  /* sl_macro */
 #define cc_data19_assign_at(src, sat, dst)         , data19_assign_at(src, sat, dst)  /* sl_macro */
 #define cc_data19_null(e)                          , data19_null(e)  /* sl_macro */
 #define cc_data19_inc(e)                           , data19_inc(e)  /* sl_macro */
 #define cc_data19_dec(e)                           , data19_dec(e)  /* sl_macro */
 #define cc_data19_add(e, n)                        , data19_add(e, n)  /* sl_macro */
 #define cc_data19_sub(e, n)                        , data19_sub(e, n)  /* sl_macro */
 #define cc_data19_copy(src, dst)                   , data19_copy(src, dst)  /* sl_macro */
 #define cc_data19_ncopy(src, dst, n)               , data19_ncopy(src, dst, n)  /* sl_macro */
 #define cc_data19_nmove(src, dst, n)               , data19_nmove(src, dst, n)  /* sl_macro */
 #define cc_data19_copy_at(src, sat, dst, dat)      , data19_copy_at(src, sat, dst, dat)  /* sl_macro */
 #define cc_data19_ncopy_at(src, sat, dst, dat, n)  , data19_ncopy_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data19_nmove_at(src, sat, dst, dat, n)  , data19_nmove_at(src, sat, dst, dat, n)  /* sl_macro */
 #define cc_data19_xchange(e0, e1, t)               , data19_xchange(e0, e1, t)  /* sl_macro */
 #define cc_data19_xchange_at(e0, at0, e1, at1, t)  , data19_xchange_at(e0, at0, e1, at1, t)  /* sl_macro */
#else /* SL_DATA */
 #define SL_DATA
 #define cc_data19_assign(src, dst)                 data19_assign(src, dst)
 #define cc_data19_assign_at(src, sat, dst)         data19_assign_at(src, sat, dst)
 #define cc_data19_null(e)                          data19_null(e)
 #define cc_data19_inc(e)                           data19_inc(e)
 #define cc_data19_dec(e)                           data19_dec(e)
 #define cc_data19_add(e, n)                        data19_add(e, n)
 #define cc_data19_sub(e, n)                        data19_sub(e, n)
 #define cc_data19_copy(src, dst)                   data19_copy(src, dst)
 #define cc_data19_ncopy(src, dst, n)               data19_ncopy(src, dst, n)
 #define cc_data19_nmove(src, dst, n)               data19_nmove(src, dst, n)
 #define cc_data19_copy_at(src, sat, dst, dat)      data19_copy_at(src, sat, dst, dat)
 #define cc_data19_ncopy_at(src, sat, dst, dat, n)  data19_ncopy_at(src, sat, dst, dat, n)
 #define cc_data19_nmove_at(src, sat, dst, dat, n)  data19_nmove_at(src, sat, dst, dat, n)
 #define cc_data19_xchange(e0, e1, t)               data19_xchange(e0, e1, t)
 #define cc_data19_xchange_at(e0, at0, e1, at1, t)  data19_xchange_at(e0, at0, e1, at1, t)
#endif /* SL_DATA */

#else /* SL_DATA19 */

 #define data19_n                                   0
 #define data19_byte                                0
/* #define data19_ptr(e)*/

 #define data19_byte_flex                           0
 #define data19_weight                              0

 /* commands for regular use */
 #define data19_assign(src, dst)                    Z_NOP()
 #define data19_assign_at(src, sat, dst)            Z_NOP()
 #define data19_null(e)                             Z_NOP()
 #define data19_inc(e)                              Z_NOP()
 #define data19_dec(e)                              Z_NOP()
 #define data19_add(e, n)                           Z_NOP()
 #define data19_sub(e, n)                           Z_NOP()
 #define data19_copy(src, dst)                      Z_NOP()
 #define data19_ncopy(src, dst, n)                  Z_NOP()
 #define data19_nmove(src, dst, n)                  Z_NOP()
 #define data19_copy_at(src, sat, dst, dat)         Z_NOP()
 #define data19_ncopy_at(src, sat, dst, dat, n)     Z_NOP()
 #define data19_nmove_at(src, sat, dst, dat, n)     Z_NOP()
 #define data19_xchange(e0, e1, t)                  Z_NOP()
 #define data19_xchange_at(e0, at0, e1, at1, t)     Z_NOP()

 /* chained command versions */
 #define cc_data19_assign(src, dst)
 #define cc_data19_assign_at(src, sat, dst)
 #define cc_data19_null(e)
 #define cc_data19_inc(e)
 #define cc_data19_dec(e)
 #define cc_data19_add(e, n)
 #define cc_data19_sub(e, n)
 #define cc_data19_copy(src, dst)
 #define cc_data19_ncopy(src, dst, n)
 #define cc_data19_nmove(src, dst, n)
 #define cc_data19_copy_at(src, sat, dst, dat)
 #define cc_data19_ncopy_at(src, sat, dst, dat, n)
 #define cc_data19_nmove_at(src, sat, dst, dat, n)
 #define cc_data19_xchange(e0, e1, t)
 #define cc_data19_xchange_at(e0, at0, e1, at1, t)

#endif /* SL_DATA19 */

#define data19_cm                                   SLCM_DATA19  /* sl_macro */

/* DATAX_TEMPLATE_END */


#endif /* __SL_DATA_SINGLES_H__ */
