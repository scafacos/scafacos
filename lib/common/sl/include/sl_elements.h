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

/* SL_DATA_IGNORE -> prevents that the present data is processed */


#ifndef __SL_ELEMENTS_H__
#define __SL_ELEMENTS_H__


/* prevent the following functions from processing the DATAs (even though they may exist) */
#ifdef SL_DATA_IGNORE
 /* disable the single DATAs */
/* DATAX_TEMPLATE_BEGIN */
 #undef SL_DATA0
 #undef SL_DATA1
 #undef SL_DATA2
 #undef SL_DATA3
 #undef SL_DATA4
 #undef SL_DATA5
 #undef SL_DATA6
 #undef SL_DATA7
 #undef SL_DATA8
 #undef SL_DATA9
 #undef SL_DATA10
 #undef SL_DATA11
 #undef SL_DATA12
 #undef SL_DATA13
 #undef SL_DATA14
 #undef SL_DATA15
 #undef SL_DATA16
 #undef SL_DATA17
 #undef SL_DATA18
 #undef SL_DATA19
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

#include "sl_key.h"

#include "sl_index.h"

#undef SL_DATA

#include "sl_data_singles.h"

#include "sl_data.h"

#undef SL_DATA


#define elem_n                                          (key_n + index_n + data_n)  /* sl_macro */
#define elem_byte                                       (key_byte + index_byte + data_byte)  /* sl_macro */

#define elem_key_at(_s_, _sat_)                         (key_at((_s_)->keys, _sat_))  /* sl_macro */

#define elem_assign(_s_, _d_)                           (*(_d_) = *(_s_))  /* sl_macro */
#define elem_assign_at(_s_, _sat_, _d_)                 ((_d_)->size = (_s_)->size - (_sat_), (_d_)->max_size = (_s_)->max_size - (_sat_), \
                                                         key_assign_at((_s_)->keys, _sat_, (_d_)->keys) \
                                                         cc_index_assign_at((_s_)->indices, _sat_, (_d_)->indices) \
                                                         cc_data_assign_at(_s_, _sat_, _d_))  /* sl_macro elem_assign_at */
#define elem_null(_e_)                                  ((_e_)->size = 0, (_e_)->max_size = 0, \
                                                         key_null((_e_)->keys) \
                                                         cc_index_null((_e_)->indices) \
                                                         cc_data_null(_e_))  /* sl_macro elem_null */
#define elem_inc(_e_)                                   (key_inc((_e_)->keys) \
                                                         cc_index_inc((_e_)->indices) \
                                                         cc_data_inc(_e_))  /* sl_macro elem_inc */
#define elem_dec(_e_)                                   (key_dec((_e_)->keys) \
                                                         cc_index_dec((_e_)->indices) \
                                                         cc_data_dec(_e_))  /* sl_macro elem_dec */
#define elem_add(_e_, _n_)                              (key_add((_e_)->keys, _n_) \
                                                         cc_index_add((_e_)->indices, _n_) \
                                                         cc_data_add(_e_, _n_))  /* sl_macro elem_add */
#define elem_sub(_e_, _n_)                              (key_sub((_e_)->keys, _n_) \
                                                         cc_index_sub((_e_)->indices, _n_) \
                                                         cc_data_sub(_e_, _n_))  /* sl_macro elem_sub */

#define elem_copy(_s_, _d_)                             (key_copy((_s_)->keys, (_d_)->keys) \
                                                         cc_index_copy((_s_)->indices, (_d_)->indices) \
                                                         cc_data_copy(_s_, _d_))  /* sl_macro elem_copy */
#define elem_ncopy(_s_, _d_, _n_)                       (key_ncopy((_s_)->keys, (_d_)->keys, _n_) \
                                                         cc_index_ncopy((_s_)->indices, (_d_)->indices, _n_) \
                                                         cc_data_ncopy(_s_, _d_, _n_))  /* sl_macro elem_ncopy */
#define elem_nmove(_s_, _d_, _n_)                       (key_nmove((_s_)->keys, (_d_)->keys, _n_) \
                                                         cc_index_nmove((_s_)->indices, (_d_)->indices, _n_) \
                                                         cc_data_nmove(_s_, _d_, _n_))  /* sl_macro elem_nmove */

#define elem_copy_at(_s_, _sat_, _d_, _dat_)            (key_copy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_) \
                                                         cc_index_copy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_) \
                                                         cc_data_copy_at(_s_, _sat_, _d_, _dat_))  /* sl_macro elem_copy_at */
#define elem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      (key_ncopy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_) \
                                                         cc_index_ncopy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_) \
                                                         cc_data_ncopy_at(_s_, _sat_, _d_, _dat_, _n_))  /* sl_macro elem_ncopy_at */
#define elem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      (key_nmove_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_) \
                                                         cc_index_ncopy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_) \
                                                         cc_data_nmove_at(_s_, _sat_, _d_, _dat_, _n_))  /* sl_macro elem_nmove_at */

#define elem_xchange(_e0_, _e1_, _t_)                   (key_xchange((_e0_)->keys, (_e1_)->keys, (_t_)->keys) \
                                                         cc_index_xchange((_e0_)->indices, (_e1_)->indices, (_t_)->indices) \
                                                         cc_data_xchange(_e0_, _e1_, _t_))  /* sl_macro elem_xchange */
#define elem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  (key_xchange_at((_e0_)->keys, _at0_, (_e1_)->keys, _at1_, (_t_)->keys) \
                                                         cc_index_xchange_at((_e0_)->indices, _at0_, (_e1_)->indices, _at1_, (_t_)->indices) \
                                                         cc_data_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_))  /* sl_macro elem_xchange_at */

#ifdef sl_elem_weight
# define elem_has_weight                                1  /* sl_macro */
# define elem_weight_ifelse(_if_, _el_)                 (_if_)  /* sl_macro */
# define elem_weight(_e_, _at_)                         ((slweight_t) sl_elem_weight((_e_), (_at_)))  /* sl_macro */
# define elem_weight_one(_e_, _at_)                     ((slweight_t) sl_elem_weight((_e_), (_at_)))  /* sl_macro */
# ifdef sl_elem_weight_set
#  define elem_weight_set(_e_, _at_, _w_)               sl_elem_weight_set((_e_), (_at_), (_w_))  /* sl_macro */
# else
#  define elem_weight_set(_e_, _at_, _w_)               sl_elem_weight((_e_), (_at_)) = (_w_)
# endif
#else
# define elem_has_weight                                0
# define elem_weight_ifelse(_if_, _el_)                 (_el_)
# define elem_weight_one(_e_, _at_)                     ((slweight_t) 1)
# define elem_weight_set(_e_, _at_, _w_)                Z_NOP()
#endif

#define elem_pack(_s_, _d_)                             (key_copy((_s_)->keys, &(_d_)->elements[0].key) cc_data_copy(_d_, &(_s_)->elements[0]))  /* sl_macro */
#define elem_pack_at(_s_, _sat_, _d_, _dat_)            (key_copy_at((_s_)->keys, _sat_, &(_d_)->elements[_dat_].key, 0) cc_data_copy_at(_s_, _sat_, &(_d_)->elements[_dat_], 0))  /* sl_macro */

#define elem_npack(_s_, _d_, _n_)                       elem_npack_at(_s_, 0, _d_, 0, _n_)  /* sl_macro */


#endif /* __SL_ELEMENTS_H__ */
