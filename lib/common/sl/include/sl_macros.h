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


#ifndef __SL_MACROS_H__
#define __SL_MACROS_H__


/* count */

/* sl_macro SL_KEY2CLASS_DECL_COUNT_DB */
#define SL_KEY2CLASS_DECL_COUNT_DB \
  struct { slkey_t *ki, *kend; slint_t p; } k2c0cd;

/* sl_macro SL_KEY2CLASS_EXEC_COUNT_DB */
#define SL_KEY2CLASS_EXEC_COUNT_DB(_k2c_, _s_, _nc_, _cs_)  do { \
  for (k2c0cd.p = 0; k2c0cd.p < (_nc_); k2c0cd.p++) (_cs_)[k2c0cd.p] = 0; \
  key_assign_at((_s_)->keys, (_s_)->size, k2c0cd.kend); \
  for (key_assign((_s_)->keys, k2c0cd.ki); k2c0cd.ki < k2c0cd.kend; key_inc(k2c0cd.ki)) { \
    k2c0cd.p = (_k2c_)(*k2c0cd.ki); \
    ++(_cs_)[k2c0cd.p]; \
  } } while (0)


#endif /* __SL_MACROS_H__ */
