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


#ifndef __SPEC_CORE_H__
#define __SPEC_CORE_H__


#include "spec_conf.h"

#ifdef SP_RENAME
# include "spec_rename.h"
#endif

#include "spec_public.h"

#include "z_pack.h"


/* sp_type spec_int_t */
/* sp_macro spec_int_fmt */
/* sp_type spec_proc_t */
/* sp_macro spec_proc_fmt */
/* sp_macro SPEC_LOC_NONE SPEC_PROC_NONE */
/* sp_type spec_tloc_data_t spec_tproc_data_t */
/* sp_type spec_elem_buf_t spec_elem_t */
/* sp_type spec_elem_index_t */
/* sp_macro spec_elem_index_fmt spidx_fmt */
/* sp_macro spec_elem_alloc_buf spec_elem_free_buf spec_elem_alloc_tmp spec_elem_free_tmp spec_elem_sizefor spec_elem_copy_type */
/* sp_macro spec_elem_set_n spec_elem_get_n spec_elem_set_nmax spec_elem_get_nmax spec_elem_set_buf spec_elem_get_buf spec_elem_x_get_buf */
/* sp_macro spec_elem_copy_at spec_elem_ncopy_at spec_elem_exchange_at */
/* sp_macro spec_elem_alltoallv_db spec_elem_alltoallv_ip spec_elem_alltoallw_indexed_db */
/* sp_macro spec_elem_alltoallv_proclists_db spec_elem_alltoallv_proclists_ip spec_elem_alltoallw_indexed_proclists_db */


typedef spec_int_t spint_t;
#define spint_fmt  spec_int_fmt
#define MPI_SPINT  spec_int_mpi

typedef spec_proc_t sproc_t;
#define sproc_fmt  spec_proc_fmt

typedef spec_elem_index_t spidx_t;
#define spidx_fmt  spec_elem_index_fmt


#ifndef spec_elem_alloc_tmp
# define spec_elem_alloc_tmp(_e_, _n_)  spec_elem_alloc_buf(_e_, _n_)
#endif
#ifndef spec_elem_free_tmp
# define spec_elem_free_tmp(_e_)        spec_elem_free_buf(_e_)
#endif


#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
extern double *spec_timing;
#endif


/* sp_type _spec_tproc_t spec_tproc_t */
typedef struct _spec_tproc_t
{
  spint_t max_tprocs;

  spec_tproc_f *tproc;
  spec_tproc_ext_t tproc_ext;

  spec_tproc_mod_f *tproc_mod;
  spec_tproc_mod_ext_t tproc_mod_ext;

  spec_tprocs_f *tprocs;
  spec_tprocs_ext_t tprocs_ext;

  spec_tprocs_mod_f *tprocs_mod;
  spec_tprocs_mod_ext_t tprocs_mod_ext;

  spec_tproc_reset_f *reset;

#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, nrecv_procs;
  sproc_t *send_procs, *recv_procs;
#endif

} *spec_tproc_t;

/* sp_macro SPEC_TPROC_NULL */
#define SPEC_TPROC_NULL  NULL


spint_t spec_tproc_create(spec_tproc_t *tproc, spec_tproc_f *func, spec_tproc_mod_f *func_mod, spec_tprocs_f *func_s, spec_tprocs_mod_f *func_s_mod, spint_t max_tprocs);
spint_t spec_tproc_destroy(spec_tproc_t *tproc);
spint_t spec_tproc_duplicate(spec_tproc_t tproc, spec_tproc_t *newtproc);

spint_t spec_tproc_set_tproc(spec_tproc_t tproc, spec_tproc_f *func);
spint_t spec_tproc_set_ext_tproc(spec_tproc_t tproc, const spec_tproc_ext_t *tproc_ext);

spint_t spec_tproc_set_tproc_mod(spec_tproc_t tproc, spec_tproc_mod_f *func_mod);
spint_t spec_tproc_set_ext_tproc_mod(spec_tproc_t tproc, const spec_tproc_mod_ext_t *tproc_mod_ext);

spint_t spec_tproc_set_tprocs(spec_tproc_t tproc, spec_tprocs_f *func_s, spint_t max_tprocs);
spint_t spec_tproc_set_ext_tprocs(spec_tproc_t tproc, const spec_tprocs_ext_t *tprocs_ext);

spint_t spec_tproc_set_tprocs_mod(spec_tproc_t tproc, spec_tprocs_mod_f *func_s_mod, spint_t max_tprocs);
spint_t spec_tproc_set_ext_tprocs_mod(spec_tproc_t tproc, const spec_tprocs_mod_ext_t *tprocs_mod_ext);

spint_t spec_tproc_set_reset(spec_tproc_t tproc, spec_tproc_reset_f *reset);

#ifdef SPEC_PROCLISTS
void spec_make_recv_proclist(spint_t nsend_procs, sproc_t *send_procs, spint_t *nrecv_procs, sproc_t **recv_procs, int size, int rank, MPI_Comm comm);
spint_t spec_tproc_set_proclists(spec_tproc_t tproc, spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs, int size, int rank, MPI_Comm comm);
#endif

spint_t spec_print(spec_tproc_t tproc, spec_tproc_data_t tproc_data, spec_elem_t *b);

#ifdef SPEC_ALLTOALLV
spint_t spec_alltoallv_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm);
spint_t spec_alltoallv_ip(spec_elem_t *b, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm);
#endif

#ifdef SPEC_ALLTOALLW
spint_t spec_alltoallw_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm);
#endif

#ifdef SPEC_PUT
spint_t spec_put_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm);
spint_t spec_put_2phases_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm);
#endif

#ifdef SPEC_SENDRECV
extern void *spec_sendrecv_aux;
extern spint_t spec_sendrecv_aux_size, spec_sendrecv_send_requests, spec_sendrecv_receive_requests;
spint_t spec_sendrecv_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm);
spint_t spec_sendrecv_db2(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm);
#endif


#endif /* __SPEC_CORE_H__ */
