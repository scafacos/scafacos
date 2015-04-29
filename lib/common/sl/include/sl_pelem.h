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


#ifndef __SL_PELEM_H__
#define __SL_PELEM_H__


#define pelem_byte                                       (sizeof(packed_element_t) + data_byte_flex)

#define pelem_at(_s_, _sat_)                             ((_s_)->elements + (_sat_))
#define pelem_key_at(_s_, _sat_)                         (&pelem_at((_s_), (_sat_))->key)

#define pelem_assign(_s_, _d_)                           (*(_d_) = *(_s_))
#define pelem_assign_at(_s_, _sat_, _d_)                 ((_d_)->size = (_s_)->size - (_sat_), (_d_)->max_size = (_s_)->max_size - (_sat_), (_d_)->elements = (_s_)->elements + (_sat_))
#define pelem_null(_e_)                                  ((_e_)->size = 0, (_e_)->max_size = 0, (_e_)->elements = NULL)
#define pelem_inc(_e_)                                   (++(_e_)->elements)
#define pelem_dec(_e_)                                   (--(_e_)->elements)
#define pelem_add(_e_, _n_)                              ((_e_)->elements += (_n_))
#define pelem_sub(_e_, _n_)                              ((_e_)->elements -= (_n_))

#define pelem_copy(_s_, _d_)                             (*(_d_)->elements = *(_s_)->elements)
#define pelem_ncopy(_s_, _d_, _n_)                       (memcpy((_d_)->elements, (_s_)->elements, (_n_) * pelem_byte))
#define pelem_nmove(_s_, _d_, _n_)                       (memmove((_d_)->elements, (_s_)->elements, (_n_) * pelem_byte))

#define pelem_copy_at(_s_, _sat_, _d_, _dat_)            ((_d_)->elements[_dat_] = (_s_)->elements[_sat_])
#define pelem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      (memcpy(&(_d_)->elements[_dat_], &(_s_)->elements[_sat_], (_n_) * pelem_byte))
#define pelem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      (memmove(&(_d_)->elements[_dat_], &(_s_)->elements[_sat_], (_n_) * pelem_byte))

#define pelem_xchange(_e0_, _e1_, _t_)                   (pelem_copy(_e0_, _t_), pelem_copy(_e1_, _e0_), pelem_copy(_t_, _e1_))
#define pelem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  (pelem_copy_at(_e0_, _at0_, _t_, 0), pelem_copy_at(_e1_, _at1_, _e0_, _at0_), pelem_copy_at(_t_, 0, _e1_, _at1_))

#define pelem_weight(_e_, _at_)                          sl_elem_weight((_e_), (_at_))

#define pelem_unpack(_s_, _d_)                           (key_copy(&(_s_)->elements[0].key, (_d_)->keys) cc_data_copy(&(_s_)->elements[0], _d_))
#define pelem_unpack_at(_s_, _sat_, _d_, _dat_)          (key_copy_at(&(_s_)->elements[_sat_].key, 0, (_d_)->keys, _dat_) cc_data_copy_at(&(_s_)->elements[_sat_], 0, _d_, _dat_))

#define pelem_nunpack(_s_, _d_, _n_)                     pelem_nunpack_at(_s_, 0, _d_, 0, _n_)


#endif /* __SL_ELEMENTS_H__ */
