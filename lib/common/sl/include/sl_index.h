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


#ifndef __SL_INDEX_H__
#define __SL_INDEX_H__


/* sl_macro SL_INDEX SL_PACKED_INDEX sl_index_type_c sl_index_type_mpi sl_index_size_mpi sl_index_type_fmt sl_index_integer sl_index_memcpy */

#ifdef SL_INDEX

# define sl_index_byte                                       ((slint_t) sizeof(sl_index_type_c))  /* sl_macro */

# ifndef sl_index_copy
#  ifndef sl_index_memcpy
#   define sl_index_copy(_s_, _d_)                           SL_ARRAY1_COPY(_s_, _d_)
#  else
#   define sl_index_copy(_s_, _d_)                           memcpy(_d_, _s_, sl_index_byte)  /* sl_macro */
#  endif
# endif
# ifndef sl_index_ncopy
#  define sl_index_ncopy(_s_, _d_, _n_)                      memcpy(_d_, _s_, (_n_) * sl_index_byte)  /* sl_macro */
# endif
# ifndef sl_index_nmove
#  define sl_index_nmove(_s_, _d_, _n_)                      memmove(_d_, _s_, (_n_) * sl_index_byte)  /* sl_macro */
# endif


# define index_type_c                                        sl_index_type_c  /* sl_macro */
# define index_type_mpi                                      (sl_index_type_mpi)  /* sl_macro */
# define index_size_mpi                                      (sl_index_size_mpi)  /* sl_macro */
# define index_type_fmt                                      sl_index_type_fmt  /* sl_macro */

# define index_n                                             1  /* sl_macro */
# define index_byte                                          (sl_index_byte)  /* sl_macro */

/* commands for regular use */
# define index_assign(_s_, _d_)                              (_d_ = _s_)  /* sl_macro */
# define index_assign_at(_s_, _sat_, _d_)                    (_d_ = &_s_[_sat_])  /* sl_macro */
# define index_null(_i_)                                     (_i_ = NULL)  /* sl_macro */
# define index_inc(_i_)                                      (++_i_)  /* sl_macro */
# define index_dec(_i_)                                      (--_i_)  /* sl_macro */
# define index_add(_i_, _n_)                                 (_i_ += _n_)  /* sl_macro */
# define index_sub(_i_, _n_)                                 (_i_ -= _n_)  /* sl_macro */

# define index_copy(_s_, _d_)                                sl_index_copy(_s_, _d_)  /* sl_macro */
# define index_ncopy(_s_, _d_, _n_)                          sl_index_ncopy(_s_, _d_, _n_)  /* sl_macro */
# define index_nmove(_s_, _d_, _n_)                          sl_index_nmove(_s_, _d_, _n_)  /* sl_macro */

# define index_copy_at(_s_, _sat_, _d_, _dat_)               index_copy(&(_s_)[_sat_], &(_d_)[_dat_])  /* sl_macro */
# define index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)         index_ncopy(&(_s_)[_sat_], &(_d_)[_dat_], _n_)  /* sl_macro */
# define index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)         index_nmove(&(_s_)[_sat_], &(_d_)[_dat_], _n_)  /* sl_macro */

# define index_xchange(_i0_, _i1_, _t_)                      (index_copy(_i0_, _t_), index_copy(_i1_, _i0_), index_copy(_t_, _i1_))  /* sl_macro */
# define index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)     (index_copy_at(_i0_, _at0_, _t_, 0), index_copy_at(_i1_, _at1_, _i0_, _at0_), index_copy_at(_t_, 0, _i1_, _at1_))  /* sl_macro */

/* chained command versions */
# define cc_index_assign(_s_, _d_)                           , index_assign(_s_, _d_)  /* sl_macro */
# define cc_index_assign_at(_s_, _sat_, _d_)                 , index_assign_at(_s_, _sat_, _d_)  /* sl_macro */
# define cc_index_null(_i_)                                  , index_null(_i_)  /* sl_macro */
# define cc_index_inc(_i_)                                   , index_inc(_i_)  /* sl_macro */
# define cc_index_dec(_i_)                                   , index_dec(_i_)  /* sl_macro */
# define cc_index_add(_i_, _n_)                              , index_add(_i_, _n_)  /* sl_macro */
# define cc_index_sub(_i_, _n_)                              , index_sub(_i_, _n_)  /* sl_macro */
# define cc_index_copy(_s_, _d_)                             , index_copy(_s_, _d_)  /* sl_macro */
# define cc_index_ncopy(_s_, _d_, _n_)                       , index_ncopy(_s_, _d_, _n_)  /* sl_macro */
# define cc_index_nmove(_s_, _d_, _n_)                       , index_nmove(_s_, _d_, _n_)  /* sl_macro */
# define cc_index_copy_at(_s_, _sat_, _d_, _dat_)            , index_copy_at(_s_, _sat_, _d_, _dat_)  /* sl_macro */
# define cc_index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      , index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)  /* sl_macro */
# define cc_index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      , index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)  /* sl_macro */
# define cc_index_xchange(_i0_, _i1_, _t_)                   , index_xchange(_i0_, _i1_, _t_)  /* sl_macro */
# define cc_index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)  , index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)  /* sl_macro */

#else /* SL_INDEX */

# define index_n                                             0
# define index_byte                                          0

/* commands for regular use */
# define index_assign(_s_, _d_)                              Z_NOP()
# define index_assign_at(_s_, _sat_, _d_)                    Z_NOP()
# define index_null(_i_)                                     Z_NOP()
# define index_inc(_i_)                                      Z_NOP()
# define index_dec(_i_)                                      Z_NOP()
# define index_add(_i_, _n_)                                 Z_NOP()
# define index_sub(_i_, _n_)                                 Z_NOP()
# define index_copy(_s_, _d_)                                Z_NOP()
# define index_ncopy(_s_, _d_, _n_)                          Z_NOP()
# define index_nmove(_s_, _d_, _n_)                          Z_NOP()
# define index_copy_at(_s_, _sat_, _d_, _dat_)               Z_NOP()
# define index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)         Z_NOP()
# define index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)         Z_NOP()
# define index_xchange(_i0_, _i1_, _t_)                      Z_NOP()
# define index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)     Z_NOP()

/* chained command versions */
# define cc_index_assign(_s_, _d_)
# define cc_index_assign_at(_s_, _sat_, _d_)
# define cc_index_null(_i_)
# define cc_index_inc(_i_)
# define cc_index_dec(_i_)
# define cc_index_add(_i_, _n_)
# define cc_index_sub(_i_, _n_)
# define cc_index_copy(_s_, _d_)
# define cc_index_ncopy(_s_, _d_, _n_)
# define cc_index_nmove(_s_, _d_, _n_)
# define cc_index_copy_at(_s_, _sat_, _d_, _dat_)
# define cc_index_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)
# define cc_index_nmove_at(_s_, _sat_, _d_, _dat_, _n_)
# define cc_index_xchange(_i0_, _i1_, _t_)
# define cc_index_xchange_at(_i0_, _at0_, _i1_, _at1_, _t_)

#endif /* SL_INDEX */

#define index_cm                                             SLCM_INDICES  /* sl_macro */


#endif /* __SL_INDEX_H__ */
