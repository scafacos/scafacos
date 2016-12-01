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


#ifndef __SPEC_COMMON_H__
#define __SPEC_COMMON_H__


#include "spec_conf.h"

#ifdef SP_RENAME
# include "spec_rename.h"
#endif

#include "spec_public.h"

#include "z_pack.h"


spint_t spec_check_tproc_support(spec_tproc_t tproc, spint_t have_tproc, spint_t have_tproc_mod, spint_t have_tprocs, spint_t have_tprocs_mod, int rank, const char *name);
spint_t spec_check_buffer_size(spec_elem_t *b, spint_t min_size, spint_t allocatable, MPI_Comm comm, int rank, const char *name, const char *buf_name);

void spec_tproc_setup(spec_tproc_t tproc, spec_elem_t *b, spec_proc_t **procs, spec_elem_t **mods);
void spec_tproc_release(spec_proc_t **procs, spec_elem_t **mods);

spint_t spec_make_counts(spec_tproc_t tproc, spec_tproc_data_t tproc_data, spec_elem_t *b, int ip, int size, int *counts, spec_proc_t *procs);

extern spint_t spec_redistribute_counts_type, spec_redistribute_counts_proclists_type;

#define SPEC_REDISTRIBUTE_COUNTS_DEFAULT            0
#define SPEC_REDISTRIBUTE_COUNTS_ALLTOALL           1
#define SPEC_REDISTRIBUTE_COUNTS_2STEP              2
#define SPEC_REDISTRIBUTE_COUNTS_PUT                3
#define SPEC_REDISTRIBUTE_COUNTS_PUT_ALLOC          4
#define SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES        5
#define SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES_ALLOC  6
#define SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES        7
#define SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES_ALLOC  8

#define SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_DEFAULT            0
#define SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_IGNORE             1
#define SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ISENDIRECV         2
#define SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ALLTOALLV          3
#define SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT                4
#define SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_ALLOC          5
#define SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES        6
#define SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES_ALLOC  7

spint_t spec_redistribute_counts(int *scounts, int *rcounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm);

extern spint_t spec_reduce_scatter_counts_type, spec_reduce_scatter_counts_proclists_type;

#define SPEC_REDUCE_SCATTER_COUNTS_DEFAULT     0
#define SPEC_REDUCE_SCATTER_COUNTS_REDSCAT     1
#define SPEC_REDUCE_SCATTER_COUNTS_ACCUMULATE  2

#define SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_DEFAULT     0
#define SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_IGNORE      1
#define SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ISENDIRECV  2
#define SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ALLTOALLV   3
#define SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ACCUMULATE  4

spint_t spec_reduce_scatter_counts(int *scounts, int *rcounts, int ncounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm);

extern spint_t spec_prefix_counts_type, spec_prefix_counts_proclists_type;

#define SPEC_PREFIX_COUNTS_SCAN        0
/*#define SPEC_PREFIX_COUNTS_ACCUMULATE  1*/
#define SPEC_PREFIX_COUNTS_DEFAULT     SPEC_PREFIX_COUNTS_SCAN

#define SPEC_PREFIX_COUNTS_PROCLISTS_IGNORE      0
/*#define SPEC_PREFIX_COUNTS_PROCLISTS_ISENDIRECV  1
#define SPEC_PREFIX_COUNTS_PROCLISTS_ALLTOALLV   2
#define SPEC_PREFIX_COUNTS_PROCLISTS_ACCUMULATE  3*/
#define SPEC_PREFIX_COUNTS_PROCLISTS_DEFAULT     SPEC_PREFIX_COUNTS_PROCLISTS_IGNORE

spint_t spec_prefix_counts(int *scounts, int *rcounts, int ncounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm);


#endif /* __SPEC_COMMON_H__ */
