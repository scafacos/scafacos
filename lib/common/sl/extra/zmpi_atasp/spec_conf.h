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


#ifndef __SPEC_CONF_H__
#define __SPEC_CONF_H__


#include <mpi.h>

#ifdef HAVE_ZMPI_ATAIP_H
# include "zmpi_ataip.h"
#endif

#ifdef HAVE_ZMPI_TOOLS_H
# include "zmpi_tools.h"
#endif

#include "spec_public_conf.h"


#define SP_RENAME

#ifdef ZMPI_PREFIX
# define SP_PREFIX  SP_CONCAT(ZMPI_PREFIX, zmpi_)
#else
# define SP_PREFIX  zmpi_
#endif


#define spec_int_fmt  "d"
#define spec_int_mpi  MPI_INT

#define spec_proc_fmt  "d"

#if MPI_VERSION >= 3
# define spec_elem_index_fmt  ((sizeof(MPI_Count) == sizeof(short))?"h":(sizeof(MPI_Count) == sizeof(long))?"ld":(sizeof(MPI_Count) == sizeof(long long))?"ld":"d")
#else
# define spec_elem_index_fmt  "ld"
#endif

#define spec_elem_alloc_buf(_e_, _n_)  Z_MOP((_e_)->buf = z_alloc(zmpil_nextent((_n_), &(_e_)->zmpil_type), 1); (_e_)->count = (_e_)->max_count = (_n_);)
#define spec_elem_free_buf(_e_)        z_free((_e_)->buf)

/*#define spec_elem_alloc_rbuf(_e_)*/

#define spec_elem_alloc_tmp_from_block(_e_, _b_, _s_)  Z_MOP((_e_)->buf = _b_; (_e_)->count = (_e_)->max_count = zmpil_sizefor(_s_, &(_e_)->zmpil_type);)

#define spec_elem_sizefor(_e_, _s_)    zmpil_sizefor(_s_, &(_e_)->zmpil_type)

#define spec_elem_copy_type(_s_, _d_)  Z_MOP((_d_)->mpi_type = (_s_)->mpi_type; zmpil_copy(&(_s_)->zmpil_type, &(_d_)->zmpil_type);)

#define spec_elem_ncopy_at(_se_, _sat_, _de_, _dat_, _n_) \
  zmpil_memcpy_at((_de_)->buf, (_dat_), (_se_)->buf, (_sat_), (_n_), &(_se_)->zmpil_type)

#define spec_elem_alltoallv_db(_sb_, _sc_, _sd_, _rb_, _rc_, _rd_, _s_, _r_, _c_) \
  MPI_Alltoallv((_sb_)->buf, (_sc_), (_sd_), (_sb_)->mpi_type, (_rb_)->buf, (_rc_), (_rd_), (_rb_)->mpi_type, (_c_))

#ifdef HAVE_ZMPI_ALLTOALLV_PROCLISTS
#define spec_elem_alltoallv_proclists_db(_sb_, _sc_, _sd_, _nsp_, _sp_, _rb_, _rc_, _rd_, _nrp_, _rp_, _s_, _r_, _c_) \
  ZMPI_Alltoallv_proclists((_sb_)->buf, (_sc_), (_sd_), (_sb_)->mpi_type, (_nsp_), (_sp_), (_rb_)->buf, (_rc_), (_rd_), (_rb_)->mpi_type, (_nrp_), (_rp_), (_c_))
#endif

#ifdef HAVE_ZMPI_ATAIP_H
#define spec_elem_alltoallv_ip(_b_, _xb_, _sc_, _sd_, _rc_, _rd_, _s_, _r_, _c_) do { \
    ZMPI_Alltoallv_inplace_aux = (_xb_)->buf; \
    ZMPI_Alltoallv_inplace_aux_size = zmpil_nextent((_xb_)->count, &(_xb_)->zmpil_type); \
    ZMPI_Alltoallv_inplace(MPI_IN_PLACE, (_sc_), (_sd_), (_b_)->mpi_type, (_b_)->buf, (_rc_), (_rd_), (_b_)->mpi_type, (_c_)); \
  } while (0)
#endif

#define spec_elem_x                         1
#define spec_elem_x_type_mpi(_e_, _x_)      ((_e_)->mpi_type)
#define spec_elem_x_size_mpi(_e_, _x_)      1
#define spec_elem_x_at(_e_, _x_, _at_)      zmpil_simple_at((_e_)->buf, _at_, &(_e_)->zmpil_type)
#define spec_elem_x_extent(_e_, _x_)        zmpil_simple_extent(&(_e_)->zmpil_type)
#define spec_elem_x_nextent(_e_, _x_, _n_)  zmpil_simple_nextent((_n_), &(_e_)->zmpil_type)

#define spec_elem_isend(_e_, _at_, _n_, _p_, _t_, _rq_, _s_, _r_, _c_) \
  MPI_Isend(zmpil_simple_at((_e_)->buf, _at_, &(_e_)->zmpil_type), _n_, (_e_)->mpi_type, _p_, _t_, _c_, _rq_)

#define spec_elem_isend_first(_e_, _at_, _n_, _p_, _t_, _rq_, _s_, _r_, _c_) \
  MPI_Isend(zmpil_simple_at((_e_)->buf, _at_, &(_e_)->zmpil_type), _n_, (_e_)->mpi_type, _p_, _t_, _c_, _rq_)

#define spec_elem_isend_next(_e_, _at_, _n_, _p_, _t_, _rq_, _s_, _r_, _c_) Z_NOP()

#define spec_elem_irecv_first(_e_, _at_, _n_, _p_, _t_, _rq_, _s_, _r_, _c_) \
  MPI_Irecv(zmpil_simple_at((_e_)->buf, _at_, &(_e_)->zmpil_type), _n_, (_e_)->mpi_type, _p_, _t_, _c_, _rq_); \

#define spec_elem_irecv_next(_e_, _at_, _n_, _p_, _t_, _rq_, _s_, _r_, _c_) Z_NOP()

#define spec_elem_get_recv_count(_e_, _s_, _rc_)  MPI_Get_count(_s_, (_e_)->mpi_type, _rc_);

#define SPEC_GLOBAL_EXIT_ON_ERROR

#define SPEC_PROCLISTS

/*#define SPEC_TIMING
#define SPEC_TIMING_PRINT*/

#ifdef SPEC_TIMING
# define Z_PACK_TIMING
#endif
#ifdef SPEC_TIMING_PRINT
# define Z_TIMING_PRINT_PREFIX  "ZMPI_ATASP_TIMING: "
# define Z_TIMING_PRINT(_i_, _s_, _n_, _v_, _r_)  Z_MOP(if ((_i_) == (_r_)) z_timing_print_default(_i_, _s_, _n_, _v_, _r_);)
#endif

#define SPEC_ERROR_FILE    stderr

#define SPEC_EXIT_SUCCESS  MPI_SUCCESS
#define SPEC_EXIT_FAILED   1

#define SPEC_ALLTOALLV
#define SPEC_SENDRECV

#define SPEC_REDISTRIBUTE_COUNTS_2STEP_THRESHOLD  1024


#endif /* __SPEC_CONF_H__ */
