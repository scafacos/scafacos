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


#ifndef __SL_ADDS_H__
#define __SL_ADDS_H__


/* sl_macro elem_set_size elem_set_max_size elem_set_keys elem_set_indices */
#define elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))
/* DATAX_TEMPLATE_BEGIN */
#define elem_set_data0(_e_, _b_)     ((_e_)->data0 = (_b_))  /* sl_macro */
#define elem_set_data1(_e_, _b_)     ((_e_)->data1 = (_b_))  /* sl_macro */
#define elem_set_data2(_e_, _b_)     ((_e_)->data2 = (_b_))  /* sl_macro */
#define elem_set_data3(_e_, _b_)     ((_e_)->data3 = (_b_))  /* sl_macro */
#define elem_set_data4(_e_, _b_)     ((_e_)->data4 = (_b_))  /* sl_macro */
#define elem_set_data5(_e_, _b_)     ((_e_)->data5 = (_b_))  /* sl_macro */
#define elem_set_data6(_e_, _b_)     ((_e_)->data6 = (_b_))  /* sl_macro */
#define elem_set_data7(_e_, _b_)     ((_e_)->data7 = (_b_))  /* sl_macro */
#define elem_set_data8(_e_, _b_)     ((_e_)->data8 = (_b_))  /* sl_macro */
#define elem_set_data9(_e_, _b_)     ((_e_)->data9 = (_b_))  /* sl_macro */
#define elem_set_data10(_e_, _b_)     ((_e_)->data10 = (_b_))  /* sl_macro */
#define elem_set_data11(_e_, _b_)     ((_e_)->data11 = (_b_))  /* sl_macro */
#define elem_set_data12(_e_, _b_)     ((_e_)->data12 = (_b_))  /* sl_macro */
#define elem_set_data13(_e_, _b_)     ((_e_)->data13 = (_b_))  /* sl_macro */
#define elem_set_data14(_e_, _b_)     ((_e_)->data14 = (_b_))  /* sl_macro */
#define elem_set_data15(_e_, _b_)     ((_e_)->data15 = (_b_))  /* sl_macro */
#define elem_set_data16(_e_, _b_)     ((_e_)->data16 = (_b_))  /* sl_macro */
#define elem_set_data17(_e_, _b_)     ((_e_)->data17 = (_b_))  /* sl_macro */
#define elem_set_data18(_e_, _b_)     ((_e_)->data18 = (_b_))  /* sl_macro */
#define elem_set_data19(_e_, _b_)     ((_e_)->data19 = (_b_))  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro elem_get_size elem_get_max_size elem_get_keys elem_get_indices */
#define elem_get_size(_e_)           (_e_)->size
#define elem_get_max_size(_e_)       (_e_)->max_size
#define elem_get_keys(_e_)           (_e_)->keys
#define elem_get_indices(_e_)        (_e_)->indices
/* DATAX_TEMPLATE_BEGIN */
#define elem_get_data0(_e_)          (_e_)->data0  /* sl_macro */
#define elem_get_data1(_e_)          (_e_)->data1  /* sl_macro */
#define elem_get_data2(_e_)          (_e_)->data2  /* sl_macro */
#define elem_get_data3(_e_)          (_e_)->data3  /* sl_macro */
#define elem_get_data4(_e_)          (_e_)->data4  /* sl_macro */
#define elem_get_data5(_e_)          (_e_)->data5  /* sl_macro */
#define elem_get_data6(_e_)          (_e_)->data6  /* sl_macro */
#define elem_get_data7(_e_)          (_e_)->data7  /* sl_macro */
#define elem_get_data8(_e_)          (_e_)->data8  /* sl_macro */
#define elem_get_data9(_e_)          (_e_)->data9  /* sl_macro */
#define elem_get_data10(_e_)          (_e_)->data10  /* sl_macro */
#define elem_get_data11(_e_)          (_e_)->data11  /* sl_macro */
#define elem_get_data12(_e_)          (_e_)->data12  /* sl_macro */
#define elem_get_data13(_e_)          (_e_)->data13  /* sl_macro */
#define elem_get_data14(_e_)          (_e_)->data14  /* sl_macro */
#define elem_get_data15(_e_)          (_e_)->data15  /* sl_macro */
#define elem_get_data16(_e_)          (_e_)->data16  /* sl_macro */
#define elem_get_data17(_e_)          (_e_)->data17  /* sl_macro */
#define elem_get_data18(_e_)          (_e_)->data18  /* sl_macro */
#define elem_get_data19(_e_)          (_e_)->data19  /* sl_macro */
/* DATAX_TEMPLATE_END */

/* sl_macro elem_set_block elem_set_block_size elem_get_block elem_get_block_size */
#define elem_set_block(_e_, _b_)       ((_e_)->keys = (_b_), (_e_)->max_size = -1)
#define elem_set_block_size(_e_, _s_)  ((_e_)->size = (_s_))
#define elem_get_block(_e_)            ((void *) (((_e_)->max_size < 0)?(_e_)->keys:NULL))
#define elem_get_block_size(_e_)       (_e_)->size

/* sl_macro pelem_set_size pelem_set_max_size pelem_set_elements */
#define pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))

/* sl_macro pelem_get_size pelem_get_max_size pelem_get_elements */
#define pelem_get_size(_e_)           (_e_)->size
#define pelem_get_max_size(_e_)       (_e_)->max_size
#define pelem_get_elements(_e_)       (_e_)->elements

/* sl_macro SL_DEFCON */
#define SL_DEFCON(_v_)  (sl_default_context._v_)


#endif /* __SL_ADDS_H__ */
