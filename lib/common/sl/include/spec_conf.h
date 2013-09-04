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


#ifndef __SPEC_CONF_H__
#define __SPEC_CONF_H__


#include "sl_common.h"

#ifdef HAVE_ZMPI_ATAIP_H
# include "zmpi_ataip.h"
#endif

#ifdef HAVE_ZMPI_TOOLS_H
# include "zmpi_tools.h"
#endif

#include "spec_public_conf.h"


#define z_mpi_rank  SL_PROC_RANK

#ifdef SL_PREFIX
# define SP_RENAME
# define SP_PREFIX  SL_PREFIX
#endif


#define spec_int_fmt  sl_int_type_fmt
#define spec_int_mpi  sl_int_type_mpi

#define spec_proc_fmt  "d"

#define spec_elem_index_fmt  sl_int_type_fmt

slint_t mpi_elements_alltoall_specific_alloc_size(slint_t n);

#define spec_elem_alloc_buf(_e_, _n_)  elements_alloc((_e_), mpi_elements_alltoall_specific_alloc_size(_n_), SLCM_ALL)
#define spec_elem_free_buf(_e_)        elements_free((_e_))

#define spec_elem_copy_type(_s_, _d_)  Z_NOP()

#define spec_elem_ncopy_at(_se_, _sat_, _de_, _dat_, _n_) \
  elem_ncopy_at((_se_), (_sat_), (_de_), (_dat_), (_n_))

#define spec_elem_alltoallv_db(_sb_, _sc_, _sd_, _rb_, _rc_, _rd_, _s_, _r_, _c_) \
  mpi_elements_alltoallv_db((_sb_), (_sc_), (_sd_), (_rb_), (_rc_), (_rd_), (_s_), (_r_), (_c_))

#define spec_elem_alltoallv_proclists_db(_sb_, _sc_, _sd_, _nsp_, _sp_, _rb_, _rc_, _rd_, _nrp_, _rp_, _s_, _r_, _c_) \
  mpi_elements_alltoallv_proclists_db((_sb_), (_sc_), (_sd_), (_nsp_), (_sp_), (_rb_), (_rc_), (_rd_), (_nrp_), (_rp_), (_s_), (_r_), (_c_))

#define spec_elem_alltoallv_ip(_b_, _xb_, _sc_, _sd_, _rc_, _rd_, _s_, _r_, _c_) \
  mpi_elements_alltoallv_ip((_b_), (_xb_), (_sc_), (_sd_), (_rc_), (_rd_), (_s_), (_r_), (_c_))

#define spec_elem_alloc_rbuf(_e_)  ((_e_)->max_size <= 0)


/*#define SPEC_GLOBAL_EXIT_ON_ERROR*/

#define SPEC_PROCLIST

/*#define SPEC_TIMING*/

/*#define SPEC_ERROR_FILE*/

#define SPEC_EXIT_SUCCESS  MPI_SUCCESS
#define SPEC_EXIT_FAILED   1

#ifdef MC_ALLTOALL_INT_2STEP_THRESHOLD
# define SPEC_MPI_ALLTOALL_2STEP_THRESHOLD  MC_ALLTOALL_INT_2STEP_THRESHOLD
#endif


#endif /* __SPEC_CONF_H__ */
