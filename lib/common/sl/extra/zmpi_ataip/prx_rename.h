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


#ifndef __PRX_RENAME_H__
#define __PRX_RENAME_H__


#define PRX_CONCAT(_a_, _b_)           PRX_CONCAT_(_a_, _b_)
#define PRX_CONCAT_(_a_, _b_)          _a_##_b_

#define PRX_CONCONCAT(_a_, _b_, _c_)   PRX_CONCONCAT_(_a_, _b_, _c_)
#define PRX_CONCONCAT_(_a_, _b_, _c_)  _a_##_b_##_c_

#ifdef PRX_PREFIX
# define PRX_VAR(_v_)   PRX_CONCAT(PRX_PREFIX, _v_)
# define PRX_FUNC(_f_)  PRX_CONCAT(PRX_PREFIX, _f_)
#else
# define PRX_VAR(_v_)   _v_
# define PRX_FUNC(_f_)  _f_
#endif


/* prx.c */
#define prx_seed  PRX_FUNC(prx_seed)
#define prx_permutation  PRX_FUNC(prx_permutation)
#define prx_enumerate_create  PRX_FUNC(prx_enumerate_create)
#define prx_enumerate_destroy  PRX_FUNC(prx_enumerate_destroy)
#define prx_enumerate_print  PRX_FUNC(prx_enumerate_print)
#define prx_enumerate  PRX_FUNC(prx_enumerate)


#endif /* __PRX_RENAME_H__ */
