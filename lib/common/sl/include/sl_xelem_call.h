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


#define xelem_name                                       keys
#define xelem_name_packed                                key
#define xelem_name_str                                   "keys"
#define xelem_type_c                                     key_type_c
#define xelem_size_c                                     1
#define xelem_type_mpi                                   key_type_mpi
#define xelem_size_mpi                                   key_size_mpi
#define xelem_sltype_t                                   slkey_t
#define xelem_byte                                       key_byte
#define xelem_weight                                     0
#define xelem_assign(_s_, _d_)                           key_assign((_s_)->keys, (_d_)->keys)
#define xelem_assign_at(_s_, _sat_, _d_)                 key_assign_at((_s_)->keys, _sat_, (_d_)->keys)
#define xelem_null(_e_)                                  key_null((_e_)->keys)
#define xelem_inc(_e_)                                   key_inc((_e_)->keys)
#define xelem_dec(_e_)                                   key_dec((_e_)->keys)
#define xelem_add(_e_, _n_)                              key_add((_e_)->keys, _n_)
#define xelem_sub(_e_, _n_)                              key_sub((_e_)->keys, _n_)
#define xelem_copy(_s_, _d_)                             key_copy((_s_)->keys, (_d_)->keys)
#define xelem_ncopy(_s_, _d_, _n_)                       key_ncopy((_s_)->keys, (_d_)->keys, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       key_nmove((_s_)->keys, (_d_)->keys, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            key_copy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      key_ncopy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      key_nmove_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   key_xchange((_e0_)->keys, (_e1_)->keys, (_t_)->keys)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  key_xchange_at((_e0_)->keys, _at0_, (_e1_)->keys, _at1_, (_t_)->keys)
#define xelem_cm                                         key_cm
#define xelem_set(_e_, _k_)                              elem_set_keys(_e_, _k_)
#define xelem_get(_e_)                                   elem_get_keys(_e_)
#define xelem_buf(_e_)                                   (_e_)->keys
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->keys[_at_])
#define xelem_mpi_datatype                               key_mpi_datatype

#if defined(xelem_key_if)
if (xelem_key_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if !defined(xelem_key_not)
# ifdef xelem_call_key
xelem_call_key
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_key_not
xelem_call_key_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_key_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       indices
#define xelem_name_packed                                index
#define xelem_name_str                                   "indices"
#define xelem_type_c                                     index_type_c
#define xelem_size_c                                     1
#define xelem_type_mpi                                   index_type_mpi
#define xelem_size_mpi                                   index_size_mpi
#define xelem_sltype_t                                   slindex_t
#define xelem_byte                                       index_byte
#define xelem_weight                                     0
#define xelem_assign(_s_, _d_)                           index_assign((_s_)->indices, (_d_)->indices)
#define xelem_assign_at(_s_, _sat_, _d_)                 index_assign_at((_s_)->indices, _sat_, (_d_)->indices)
#define xelem_null(_e_)                                  index_null((_e_)->indices)
#define xelem_inc(_e_)                                   index_inc((_e_)->indices)
#define xelem_dec(_e_)                                   index_dec((_e_)->indices)
#define xelem_add(_e_, _n_)                              index_add((_e_)->indices, _n_)
#define xelem_sub(_e_, _n_)                              index_sub((_e_)->indices, _n_)
#define xelem_copy(_s_, _d_)                             index_copy((_s_)->indices, (_d_)->indices)
#define xelem_ncopy(_s_, _d_, _n_)                       index_ncopy((_s_)->indices, (_d_)->indices, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       index_nmove((_s_)->indices, (_d_)->indices, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            index_copy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      index_ncopy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      index_nmove_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   index_xchange((_e0_)->indices, (_e1_)->indices, (_t_)->indices)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  index_xchange_at((_e0_)->indices, _at0_, (_e1_)->indices, _at1_, (_t_)->indices)
#define xelem_cm                                         index_cm
#define xelem_set(_e_, _i_)                              elem_set_indices(_e_, _i_)
#define xelem_get(_e_)                                   elem_get_indices(_e_)
#define xelem_buf(_e_)                                   (_e_)->indices
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->indices[_at_])
#define xelem_mpi_datatype                               index_mpi_datatype

#if defined(xelem_index_if)
if (xelem_index_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_INDEX) && !defined(xelem_index_not)
# ifdef xelem_call_index
xelem_call_index
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_index_not
xelem_call_index_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_index_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


/* DATAX_TEMPLATE_BEGIN */

#define xelem_name                                       data0
#define xelem_name_packed                                data0
#define xelem_name_str                                   "data0"
#define xelem_type_c                                     data0_type_c
#define xelem_size_c                                     data0_size_c
#define xelem_type_mpi                                   data0_type_mpi
#define xelem_size_mpi                                   data0_size_mpi
#define xelem_sltype_t                                   sldata0_t
#define xelem_byte                                       data0_byte
#define xelem_weight                                     data0_weight
#define xelem_assign(_s_, _d_)                           data0_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data0_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data0_null(_e_)
#define xelem_inc(_e_)                                   data0_inc(_e_)
#define xelem_dec(_e_)                                   data0_dec(_e_)
#define xelem_add(_e_, _n_)                              data0_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data0_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data0_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data0_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data0_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data0_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data0_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data0_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data0_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data0_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data0_cm
#define xelem_set(_e_, _b_)                              elem_set_data0(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data0(_e_)
#define xelem_buf(_e_)                                   (_e_)->data0
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data0[(_at_) * data0_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[0]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA0) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data0_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data1
#define xelem_name_packed                                data1
#define xelem_name_str                                   "data1"
#define xelem_type_c                                     data1_type_c
#define xelem_size_c                                     data1_size_c
#define xelem_type_mpi                                   data1_type_mpi
#define xelem_size_mpi                                   data1_size_mpi
#define xelem_sltype_t                                   sldata1_t
#define xelem_byte                                       data1_byte
#define xelem_weight                                     data1_weight
#define xelem_assign(_s_, _d_)                           data1_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data1_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data1_null(_e_)
#define xelem_inc(_e_)                                   data1_inc(_e_)
#define xelem_dec(_e_)                                   data1_dec(_e_)
#define xelem_add(_e_, _n_)                              data1_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data1_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data1_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data1_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data1_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data1_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data1_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data1_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data1_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data1_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data1_cm
#define xelem_set(_e_, _b_)                              elem_set_data1(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data1(_e_)
#define xelem_buf(_e_)                                   (_e_)->data1
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data1[(_at_) * data1_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[1]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA1) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data1_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data2
#define xelem_name_packed                                data2
#define xelem_name_str                                   "data2"
#define xelem_type_c                                     data2_type_c
#define xelem_size_c                                     data2_size_c
#define xelem_type_mpi                                   data2_type_mpi
#define xelem_size_mpi                                   data2_size_mpi
#define xelem_sltype_t                                   sldata2_t
#define xelem_byte                                       data2_byte
#define xelem_weight                                     data2_weight
#define xelem_assign(_s_, _d_)                           data2_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data2_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data2_null(_e_)
#define xelem_inc(_e_)                                   data2_inc(_e_)
#define xelem_dec(_e_)                                   data2_dec(_e_)
#define xelem_add(_e_, _n_)                              data2_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data2_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data2_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data2_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data2_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data2_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data2_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data2_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data2_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data2_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data2_cm
#define xelem_set(_e_, _b_)                              elem_set_data2(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data2(_e_)
#define xelem_buf(_e_)                                   (_e_)->data2
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data2[(_at_) * data2_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[2]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA2) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data2_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data3
#define xelem_name_packed                                data3
#define xelem_name_str                                   "data3"
#define xelem_type_c                                     data3_type_c
#define xelem_size_c                                     data3_size_c
#define xelem_type_mpi                                   data3_type_mpi
#define xelem_size_mpi                                   data3_size_mpi
#define xelem_sltype_t                                   sldata3_t
#define xelem_byte                                       data3_byte
#define xelem_weight                                     data3_weight
#define xelem_assign(_s_, _d_)                           data3_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data3_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data3_null(_e_)
#define xelem_inc(_e_)                                   data3_inc(_e_)
#define xelem_dec(_e_)                                   data3_dec(_e_)
#define xelem_add(_e_, _n_)                              data3_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data3_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data3_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data3_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data3_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data3_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data3_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data3_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data3_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data3_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data3_cm
#define xelem_set(_e_, _b_)                              elem_set_data3(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data3(_e_)
#define xelem_buf(_e_)                                   (_e_)->data3
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data3[(_at_) * data3_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[3]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA3) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data3_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data4
#define xelem_name_packed                                data4
#define xelem_name_str                                   "data4"
#define xelem_type_c                                     data4_type_c
#define xelem_size_c                                     data4_size_c
#define xelem_type_mpi                                   data4_type_mpi
#define xelem_size_mpi                                   data4_size_mpi
#define xelem_sltype_t                                   sldata4_t
#define xelem_byte                                       data4_byte
#define xelem_weight                                     data4_weight
#define xelem_assign(_s_, _d_)                           data4_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data4_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data4_null(_e_)
#define xelem_inc(_e_)                                   data4_inc(_e_)
#define xelem_dec(_e_)                                   data4_dec(_e_)
#define xelem_add(_e_, _n_)                              data4_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data4_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data4_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data4_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data4_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data4_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data4_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data4_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data4_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data4_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data4_cm
#define xelem_set(_e_, _b_)                              elem_set_data4(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data4(_e_)
#define xelem_buf(_e_)                                   (_e_)->data4
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data4[(_at_) * data4_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[4]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA4) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data4_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data5
#define xelem_name_packed                                data5
#define xelem_name_str                                   "data5"
#define xelem_type_c                                     data5_type_c
#define xelem_size_c                                     data5_size_c
#define xelem_type_mpi                                   data5_type_mpi
#define xelem_size_mpi                                   data5_size_mpi
#define xelem_sltype_t                                   sldata5_t
#define xelem_byte                                       data5_byte
#define xelem_weight                                     data5_weight
#define xelem_assign(_s_, _d_)                           data5_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data5_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data5_null(_e_)
#define xelem_inc(_e_)                                   data5_inc(_e_)
#define xelem_dec(_e_)                                   data5_dec(_e_)
#define xelem_add(_e_, _n_)                              data5_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data5_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data5_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data5_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data5_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data5_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data5_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data5_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data5_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data5_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data5_cm
#define xelem_set(_e_, _b_)                              elem_set_data5(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data5(_e_)
#define xelem_buf(_e_)                                   (_e_)->data5
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data5[(_at_) * data5_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[5]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA5) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data5_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data6
#define xelem_name_packed                                data6
#define xelem_name_str                                   "data6"
#define xelem_type_c                                     data6_type_c
#define xelem_size_c                                     data6_size_c
#define xelem_type_mpi                                   data6_type_mpi
#define xelem_size_mpi                                   data6_size_mpi
#define xelem_sltype_t                                   sldata6_t
#define xelem_byte                                       data6_byte
#define xelem_weight                                     data6_weight
#define xelem_assign(_s_, _d_)                           data6_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data6_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data6_null(_e_)
#define xelem_inc(_e_)                                   data6_inc(_e_)
#define xelem_dec(_e_)                                   data6_dec(_e_)
#define xelem_add(_e_, _n_)                              data6_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data6_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data6_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data6_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data6_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data6_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data6_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data6_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data6_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data6_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data6_cm
#define xelem_set(_e_, _b_)                              elem_set_data6(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data6(_e_)
#define xelem_buf(_e_)                                   (_e_)->data6
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data6[(_at_) * data6_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[6]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA6) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data6_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data7
#define xelem_name_packed                                data7
#define xelem_name_str                                   "data7"
#define xelem_type_c                                     data7_type_c
#define xelem_size_c                                     data7_size_c
#define xelem_type_mpi                                   data7_type_mpi
#define xelem_size_mpi                                   data7_size_mpi
#define xelem_sltype_t                                   sldata7_t
#define xelem_byte                                       data7_byte
#define xelem_weight                                     data7_weight
#define xelem_assign(_s_, _d_)                           data7_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data7_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data7_null(_e_)
#define xelem_inc(_e_)                                   data7_inc(_e_)
#define xelem_dec(_e_)                                   data7_dec(_e_)
#define xelem_add(_e_, _n_)                              data7_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data7_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data7_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data7_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data7_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data7_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data7_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data7_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data7_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data7_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data7_cm
#define xelem_set(_e_, _b_)                              elem_set_data7(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data7(_e_)
#define xelem_buf(_e_)                                   (_e_)->data7
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data7[(_at_) * data7_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[7]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA7) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data7_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data8
#define xelem_name_packed                                data8
#define xelem_name_str                                   "data8"
#define xelem_type_c                                     data8_type_c
#define xelem_size_c                                     data8_size_c
#define xelem_type_mpi                                   data8_type_mpi
#define xelem_size_mpi                                   data8_size_mpi
#define xelem_sltype_t                                   sldata8_t
#define xelem_byte                                       data8_byte
#define xelem_weight                                     data8_weight
#define xelem_assign(_s_, _d_)                           data8_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data8_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data8_null(_e_)
#define xelem_inc(_e_)                                   data8_inc(_e_)
#define xelem_dec(_e_)                                   data8_dec(_e_)
#define xelem_add(_e_, _n_)                              data8_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data8_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data8_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data8_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data8_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data8_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data8_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data8_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data8_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data8_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data8_cm
#define xelem_set(_e_, _b_)                              elem_set_data8(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data8(_e_)
#define xelem_buf(_e_)                                   (_e_)->data8
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data8[(_at_) * data8_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[8]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA8) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data8_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data9
#define xelem_name_packed                                data9
#define xelem_name_str                                   "data9"
#define xelem_type_c                                     data9_type_c
#define xelem_size_c                                     data9_size_c
#define xelem_type_mpi                                   data9_type_mpi
#define xelem_size_mpi                                   data9_size_mpi
#define xelem_sltype_t                                   sldata9_t
#define xelem_byte                                       data9_byte
#define xelem_weight                                     data9_weight
#define xelem_assign(_s_, _d_)                           data9_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data9_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data9_null(_e_)
#define xelem_inc(_e_)                                   data9_inc(_e_)
#define xelem_dec(_e_)                                   data9_dec(_e_)
#define xelem_add(_e_, _n_)                              data9_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data9_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data9_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data9_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data9_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data9_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data9_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data9_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data9_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data9_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data9_cm
#define xelem_set(_e_, _b_)                              elem_set_data9(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data9(_e_)
#define xelem_buf(_e_)                                   (_e_)->data9
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data9[(_at_) * data9_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[9]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA9) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data9_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data10
#define xelem_name_packed                                data10
#define xelem_name_str                                   "data10"
#define xelem_type_c                                     data10_type_c
#define xelem_size_c                                     data10_size_c
#define xelem_type_mpi                                   data10_type_mpi
#define xelem_size_mpi                                   data10_size_mpi
#define xelem_sltype_t                                   sldata10_t
#define xelem_byte                                       data10_byte
#define xelem_weight                                     data10_weight
#define xelem_assign(_s_, _d_)                           data10_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data10_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data10_null(_e_)
#define xelem_inc(_e_)                                   data10_inc(_e_)
#define xelem_dec(_e_)                                   data10_dec(_e_)
#define xelem_add(_e_, _n_)                              data10_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data10_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data10_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data10_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data10_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data10_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data10_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data10_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data10_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data10_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data10_cm
#define xelem_set(_e_, _b_)                              elem_set_data10(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data10(_e_)
#define xelem_buf(_e_)                                   (_e_)->data10
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data10[(_at_) * data10_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[10]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA10) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data10_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data11
#define xelem_name_packed                                data11
#define xelem_name_str                                   "data11"
#define xelem_type_c                                     data11_type_c
#define xelem_size_c                                     data11_size_c
#define xelem_type_mpi                                   data11_type_mpi
#define xelem_size_mpi                                   data11_size_mpi
#define xelem_sltype_t                                   sldata11_t
#define xelem_byte                                       data11_byte
#define xelem_weight                                     data11_weight
#define xelem_assign(_s_, _d_)                           data11_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data11_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data11_null(_e_)
#define xelem_inc(_e_)                                   data11_inc(_e_)
#define xelem_dec(_e_)                                   data11_dec(_e_)
#define xelem_add(_e_, _n_)                              data11_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data11_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data11_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data11_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data11_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data11_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data11_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data11_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data11_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data11_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data11_cm
#define xelem_set(_e_, _b_)                              elem_set_data11(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data11(_e_)
#define xelem_buf(_e_)                                   (_e_)->data11
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data11[(_at_) * data11_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[11]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA11) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data11_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data12
#define xelem_name_packed                                data12
#define xelem_name_str                                   "data12"
#define xelem_type_c                                     data12_type_c
#define xelem_size_c                                     data12_size_c
#define xelem_type_mpi                                   data12_type_mpi
#define xelem_size_mpi                                   data12_size_mpi
#define xelem_sltype_t                                   sldata12_t
#define xelem_byte                                       data12_byte
#define xelem_weight                                     data12_weight
#define xelem_assign(_s_, _d_)                           data12_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data12_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data12_null(_e_)
#define xelem_inc(_e_)                                   data12_inc(_e_)
#define xelem_dec(_e_)                                   data12_dec(_e_)
#define xelem_add(_e_, _n_)                              data12_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data12_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data12_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data12_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data12_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data12_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data12_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data12_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data12_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data12_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data12_cm
#define xelem_set(_e_, _b_)                              elem_set_data12(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data12(_e_)
#define xelem_buf(_e_)                                   (_e_)->data12
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data12[(_at_) * data12_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[12]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA12) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data12_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data13
#define xelem_name_packed                                data13
#define xelem_name_str                                   "data13"
#define xelem_type_c                                     data13_type_c
#define xelem_size_c                                     data13_size_c
#define xelem_type_mpi                                   data13_type_mpi
#define xelem_size_mpi                                   data13_size_mpi
#define xelem_sltype_t                                   sldata13_t
#define xelem_byte                                       data13_byte
#define xelem_weight                                     data13_weight
#define xelem_assign(_s_, _d_)                           data13_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data13_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data13_null(_e_)
#define xelem_inc(_e_)                                   data13_inc(_e_)
#define xelem_dec(_e_)                                   data13_dec(_e_)
#define xelem_add(_e_, _n_)                              data13_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data13_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data13_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data13_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data13_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data13_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data13_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data13_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data13_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data13_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data13_cm
#define xelem_set(_e_, _b_)                              elem_set_data13(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data13(_e_)
#define xelem_buf(_e_)                                   (_e_)->data13
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data13[(_at_) * data13_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[13]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA13) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data13_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data14
#define xelem_name_packed                                data14
#define xelem_name_str                                   "data14"
#define xelem_type_c                                     data14_type_c
#define xelem_size_c                                     data14_size_c
#define xelem_type_mpi                                   data14_type_mpi
#define xelem_size_mpi                                   data14_size_mpi
#define xelem_sltype_t                                   sldata14_t
#define xelem_byte                                       data14_byte
#define xelem_weight                                     data14_weight
#define xelem_assign(_s_, _d_)                           data14_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data14_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data14_null(_e_)
#define xelem_inc(_e_)                                   data14_inc(_e_)
#define xelem_dec(_e_)                                   data14_dec(_e_)
#define xelem_add(_e_, _n_)                              data14_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data14_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data14_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data14_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data14_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data14_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data14_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data14_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data14_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data14_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data14_cm
#define xelem_set(_e_, _b_)                              elem_set_data14(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data14(_e_)
#define xelem_buf(_e_)                                   (_e_)->data14
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data14[(_at_) * data14_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[14]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA14) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data14_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data15
#define xelem_name_packed                                data15
#define xelem_name_str                                   "data15"
#define xelem_type_c                                     data15_type_c
#define xelem_size_c                                     data15_size_c
#define xelem_type_mpi                                   data15_type_mpi
#define xelem_size_mpi                                   data15_size_mpi
#define xelem_sltype_t                                   sldata15_t
#define xelem_byte                                       data15_byte
#define xelem_weight                                     data15_weight
#define xelem_assign(_s_, _d_)                           data15_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data15_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data15_null(_e_)
#define xelem_inc(_e_)                                   data15_inc(_e_)
#define xelem_dec(_e_)                                   data15_dec(_e_)
#define xelem_add(_e_, _n_)                              data15_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data15_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data15_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data15_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data15_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data15_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data15_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data15_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data15_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data15_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data15_cm
#define xelem_set(_e_, _b_)                              elem_set_data15(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data15(_e_)
#define xelem_buf(_e_)                                   (_e_)->data15
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data15[(_at_) * data15_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[15]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA15) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data15_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data16
#define xelem_name_packed                                data16
#define xelem_name_str                                   "data16"
#define xelem_type_c                                     data16_type_c
#define xelem_size_c                                     data16_size_c
#define xelem_type_mpi                                   data16_type_mpi
#define xelem_size_mpi                                   data16_size_mpi
#define xelem_sltype_t                                   sldata16_t
#define xelem_byte                                       data16_byte
#define xelem_weight                                     data16_weight
#define xelem_assign(_s_, _d_)                           data16_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data16_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data16_null(_e_)
#define xelem_inc(_e_)                                   data16_inc(_e_)
#define xelem_dec(_e_)                                   data16_dec(_e_)
#define xelem_add(_e_, _n_)                              data16_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data16_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data16_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data16_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data16_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data16_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data16_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data16_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data16_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data16_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data16_cm
#define xelem_set(_e_, _b_)                              elem_set_data16(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data16(_e_)
#define xelem_buf(_e_)                                   (_e_)->data16
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data16[(_at_) * data16_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[16]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA16) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data16_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data17
#define xelem_name_packed                                data17
#define xelem_name_str                                   "data17"
#define xelem_type_c                                     data17_type_c
#define xelem_size_c                                     data17_size_c
#define xelem_type_mpi                                   data17_type_mpi
#define xelem_size_mpi                                   data17_size_mpi
#define xelem_sltype_t                                   sldata17_t
#define xelem_byte                                       data17_byte
#define xelem_weight                                     data17_weight
#define xelem_assign(_s_, _d_)                           data17_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data17_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data17_null(_e_)
#define xelem_inc(_e_)                                   data17_inc(_e_)
#define xelem_dec(_e_)                                   data17_dec(_e_)
#define xelem_add(_e_, _n_)                              data17_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data17_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data17_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data17_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data17_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data17_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data17_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data17_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data17_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data17_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data17_cm
#define xelem_set(_e_, _b_)                              elem_set_data17(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data17(_e_)
#define xelem_buf(_e_)                                   (_e_)->data17
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data17[(_at_) * data17_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[17]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA17) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data17_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data18
#define xelem_name_packed                                data18
#define xelem_name_str                                   "data18"
#define xelem_type_c                                     data18_type_c
#define xelem_size_c                                     data18_size_c
#define xelem_type_mpi                                   data18_type_mpi
#define xelem_size_mpi                                   data18_size_mpi
#define xelem_sltype_t                                   sldata18_t
#define xelem_byte                                       data18_byte
#define xelem_weight                                     data18_weight
#define xelem_assign(_s_, _d_)                           data18_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data18_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data18_null(_e_)
#define xelem_inc(_e_)                                   data18_inc(_e_)
#define xelem_dec(_e_)                                   data18_dec(_e_)
#define xelem_add(_e_, _n_)                              data18_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data18_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data18_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data18_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data18_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data18_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data18_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data18_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data18_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data18_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data18_cm
#define xelem_set(_e_, _b_)                              elem_set_data18(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data18(_e_)
#define xelem_buf(_e_)                                   (_e_)->data18
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data18[(_at_) * data18_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[18]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA18) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data18_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype


#define xelem_name                                       data19
#define xelem_name_packed                                data19
#define xelem_name_str                                   "data19"
#define xelem_type_c                                     data19_type_c
#define xelem_size_c                                     data19_size_c
#define xelem_type_mpi                                   data19_type_mpi
#define xelem_size_mpi                                   data19_size_mpi
#define xelem_sltype_t                                   sldata19_t
#define xelem_byte                                       data19_byte
#define xelem_weight                                     data19_weight
#define xelem_assign(_s_, _d_)                           data19_assign(_s_, _d_)
#define xelem_assign_at(_s_, _sat_, _d_)                 data19_assign_at(_s_, _sat_, _d_)
#define xelem_null(_e_)                                  data19_null(_e_)
#define xelem_inc(_e_)                                   data19_inc(_e_)
#define xelem_dec(_e_)                                   data19_dec(_e_)
#define xelem_add(_e_, _n_)                              data19_add(_e_, _n_)
#define xelem_sub(_e_, _n_)                              data19_sub(_e_, _n_)
#define xelem_copy(_s_, _d_)                             data19_copy(_s_, _d_)
#define xelem_ncopy(_s_, _d_, _n_)                       data19_ncopy(_s_, _d_, _n_)
#define xelem_nmove(_s_, _d_, _n_)                       data19_nmove(_s_, _d_, _n_)
#define xelem_copy_at(_s_, _sat_, _d_, _dat_)            data19_copy_at(_s_, _sat_, _d_, _dat_)
#define xelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      data19_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      data19_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
#define xelem_xchange(_e0_, _e1_, _t_)                   data19_xchange(_e0_, _e1_, _t_)
#define xelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  data19_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)
#define xelem_cm                                         data19_cm
#define xelem_set(_e_, _b_)                              elem_set_data19(_e_, _b_)
#define xelem_get(_e_)                                   elem_get_data19(_e_)
#define xelem_buf(_e_)                                   (_e_)->data19
#define xelem_buf_at(_e_, _at_)                          (&(_e_)->data19[(_at_) * data19_size_c])
#define xelem_mpi_datatype                               data_mpi_datatype[19]

#if defined(xelem_data_if)
if (xelem_data_if) {
#elif defined(xelem_if)
if (xelem_if) {
#endif

#if defined(SL_DATA19) && (!defined(xelem_data_not) || (defined(xelem_data_weight) && defined(sl_data19_weight)))
# ifdef xelem_call_data
xelem_call_data
# elif defined(xelem_call)
xelem_call
# endif
#else
# ifdef xelem_call_data_not
xelem_call_data_not
# elif defined(xelem_call_not)
xelem_call_not
# endif
#endif

#if defined(xelem_data_if) || defined(xelem_if)
} else;
#endif

#undef xelem_name
#undef xelem_name_packed
#undef xelem_name_str
#undef xelem_type_c
#undef xelem_size_c
#undef xelem_type_mpi
#undef xelem_size_mpi
#undef xelem_sltype_t
#undef xelem_byte
#undef xelem_weight
#undef xelem_assign
#undef xelem_assign_at
#undef xelem_null
#undef xelem_inc
#undef xelem_dec
#undef xelem_add
#undef xelem_sub
#undef xelem_copy
#undef xelem_ncopy
#undef xelem_nmove
#undef xelem_copy_at
#undef xelem_ncopy_at
#undef xelem_nmove_at
#undef xelem_xchange
#undef xelem_xchange_at
#undef xelem_cm
#undef xelem_set
#undef xelem_get
#undef xelem_buf
#undef xelem_buf_at
#undef xelem_mpi_datatype

/* DATAX_TEMPLATE_END */


#undef xelem_if
#undef xelem_key_if
#undef xelem_index_if
#undef xelem_data_if

#undef xelem_key_not
#undef xelem_index_not
#undef xelem_data_not

#undef xelem_data_weight

#undef xelem_call_key
#undef xelem_call_key_not
#undef xelem_call_index
#undef xelem_call_index_not
#undef xelem_call_data
#undef xelem_call_data_not
#undef xelem_call
#undef xelem_call_not
