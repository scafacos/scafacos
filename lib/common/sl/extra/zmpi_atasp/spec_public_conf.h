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


#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef HAVE_ZMPI_LOCAL_H
# include "zmpi_local.h"
#endif


typedef int spec_int_t;

typedef int spec_proc_t;

#define SPEC_PROC_NONE  MPI_PROC_NULL

typedef void *spec_tproc_data_t;

typedef void *spec_elem_buf_t;

typedef struct
{
  void *buf;
  int count, max_count;
  MPI_Datatype mpi_type;

#ifdef HAVE_ZMPI_LOCAL_H
  zmpil_t zmpil_type;
#endif

} spec_elem_t;

typedef long spec_elem_index_t;

#define spec_elem_set_n(_e_, _n_)     (_e_)->count = (_n_)
#define spec_elem_get_n(_e_)          (_e_)->count
#define spec_elem_set_nmax(_e_, _n_)  (_e_)->max_count = (_n_)
#define spec_elem_get_nmax(_e_)       (_e_)->max_count

#define spec_elem_set_buf(_e_, _b_)   (_e_)->buf = (_b_)
#define spec_elem_get_buf(_e_)        (_e_)->buf

#ifdef HAVE_ZMPI_LOCAL_H
# define spec_elem_copy_at(_se_, _sat_, _de_, _dat_)  zmpil_memcpy_at((_de_)->buf, (_dat_), (_se_)->buf, (_sat_), 1, &(_se_)->zmpil_type)
#endif

#define spec_elem_exchange_at(_s0_, _s0at_, _s1_, _s1at_, _t_) do { \
    spec_elem_copy_at((_s0_), (_s0at_), (_t_), 0); \
    spec_elem_copy_at((_s1_), (_s1at_), (_s0_), (_s0at_)); \
    spec_elem_copy_at((_t_), 0, (_s1_), (_s1at_)); \
  } while (0)


#endif /* __SPEC_PUBLIC_CONF_H__ */
