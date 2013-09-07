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


#ifndef __SPEC_PUBLIC_CONF_H__
#define __SPEC_PUBLIC_CONF_H__


#define SPEC_TLOC

typedef sl_int_type_c spec_int_t;

typedef int spec_proc_t;

#define SPEC_LOC_NONE   -1
#ifdef SL_USE_MPI
# define SPEC_PROC_NONE  MPI_PROC_NULL
#else
# define SPEC_PROC_NONE  -1
#endif

typedef void *spec_tloc_data_t;
typedef void *spec_tproc_data_t;

struct _elements_t;

typedef struct _elements_t *spec_elem_buf_t;

typedef struct _elements_t spec_elem_t;

typedef sl_int_type_c spec_elem_index_t;

#define spec_elem_set_n(_e_, _n_)     elem_set_size((_e_), (_n_))
#define spec_elem_get_n(_e_)          elem_get_size((_e_))
#define spec_elem_set_nmax(_e_, _n_)  elem_set_max_size((_e_), (_n_))
#define spec_elem_get_nmax(_e_)       elem_get_max_size((_e_))

#define spec_elem_set_buf(_e_, _b_)   *(_e_) = *(_b_)
#define spec_elem_get_buf(_e_)        (_e_)

#define spec_elem_copy_at(_se_, _sat_, _de_, _dat_) \
  elem_copy_at((_se_), (_sat_), (_de_), (_dat_))

#define spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) \
  elem_xchange_at((_s0_), (_s0at_), (_s1_), (_s1at_), (_t_))


#endif /* __SPEC_PUBLIC_CONF_H__ */
