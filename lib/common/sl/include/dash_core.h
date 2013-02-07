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


#ifndef __DASH_CORE_H__
#define __DASH_CORE_H__

#include <mpi.h>

#include "dash_conf.h"

#ifdef DS_RENAME
# include "dash_rename.h"
#endif

#ifdef HAVE_ZMPI_TOOLS_H
# include "zmpi_tools.h"
#endif

#include "z_pack.h"


#ifdef DS_TIMING
# define DS_TIMING_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define DS_TIMING_CMD(_cmd_)  Z_NOP()
#endif


extern dsint_t ds_core_run_sync;

#define DS_TIMES_RUN                 0
#define DS_TIMES_RUN_PRE             1
#define DS_TIMES_RUN_WHILE           2
#define DS_TIMES_RUN_WHILE_PRE       3
#define DS_TIMES_RUN_WHILE_MAKE      4
#define DS_TIMES_RUN_WHILE_POST      5
#define DS_TIMES_RUN_WHILE_FINISHED  6
#define DS_TIMES_RUN_POST            7
#define DS_TIMES_RUN_NTIMES          8

extern double ds_times[DS_TIMES_RUN_NTIMES];


struct _ds_t;
struct _ds_sched_t;
struct _ds_exec_t;


typedef dsint_t (*ds_sched_std_f)(struct _ds_sched_t *sched);
typedef dsint_t (*ds_sched_max_n_f)(struct _ds_sched_t *sched, dsint_t *max_nsends, dsint_t *max_nrecvs, dsint_t *max_nlocals, dsint_t *max_nsyms);


typedef struct _ds_sched_t
{
  struct _ds_t *ds;

  dsint_t round;

  dsint_t nbufs;
  struct {
    dsint_t addr_id;
  } bufs[DS_MAX_NBUFFERS];

  struct {
    dsint_t buf_id, exec_id, n;
    dspint_t *counts, *displs, *ranks;
  } send, recv;
  
  struct {
    dsint_t buf_id, buf_size;
  } aux;

  void *cxt;
  ds_sched_max_n_f max_n;
  ds_sched_std_f pre_run, post_run, finished, pre, post;

} ds_sched_t, *ds_sched_p;


dsint_t ds_sched_add_buffer(ds_sched_t *sched, dsint_t addr_id);
dsint_t ds_sched_set_send(ds_sched_t *sched, dsint_t buf_id, dsint_t n, dspint_t *scounts, dspint_t *sdispls, dspint_t *sranks, dsint_t exec_id);
dsint_t ds_sched_set_recv(ds_sched_t *sched, dsint_t buf_id, dsint_t n, dspint_t *rcounts, dspint_t *rdispls, dspint_t *rranks, dsint_t exec_id);
dsint_t ds_sched_set_aux(ds_sched_t *sched, dsint_t buf_id, dsint_t buf_size);


typedef dsint_t (*ds_exec_std_f)(struct _ds_exec_t *exec);
typedef void (*ds_exec_move_f)(struct _ds_exec_t *exec, dsint_t exec_id, dsint_t src_buf_id, dsint_t src_displs, dsint_t dst_buf_id, dsint_t dst_displs, dsint_t count);
typedef void (*ds_exec_sendrecv_replace_f)(struct _ds_exec_t *exec, int proc, dsint_t exec_id, dsint_t buf_id, dspint_t displs, dspint_t count);
/*typedef dsint_t (*ds_exec_set_max)(struct _ds_exec_t *exec, dsint_t max_nsends, dsint_t max_nrecvs, dsint_t max_nlocals);*/

typedef struct _ds_exec_t
{
  struct _ds_t *ds;

  dsint_t naddrs;
  void *addrs[DS_MAX_NBUFFERS];

  void *cxt;
  ds_exec_std_f pre_run, post_run, make;
  ds_exec_move_f move;
#ifdef WITH_SYMMETRIC
  ds_exec_sendrecv_replace_f sendrecv_replace;
#endif

} ds_exec_t, *_ds_exec_p;


dsint_t ds_exec_add_address(ds_exec_t *exec, void *addr);


typedef struct _ds_t
{
  ds_sched_t *sched;
  ds_exec_t *exec;
  
  MPI_Comm comm;
  int comm_size, comm_rank;

  struct {
    dsint_t nmax, n, execute, *buf_ids, exec_id;
    dspint_t *counts, *displs, *dsts;
  } sends;

  struct {
    dsint_t nmax, n, execute, *buf_ids, exec_id;
    dspint_t *counts, *displs, *srcs;
  } recvs;

  struct {
    dsint_t nmax, n, execute;
    dsint_t *src_buf_ids, *dst_buf_ids, *src_exec_ids, *dst_exec_ids;
    dspint_t *src_counts, *src_displs, *dst_counts, *dst_displs;
  } locals;

#ifdef WITH_SYMMETRIC
  struct {
    dsint_t nmax, n, execute, *buf_ids, exec_id;
    dspint_t *counts, *displs;
  } syms;
#endif

} ds_t, *ds_p;


dsint_t ds_create(ds_t *ds, ds_sched_t *sched, ds_exec_t *exec);
dsint_t ds_destroy(ds_t *ds);
dsint_t ds_run(ds_t *ds);


#endif /* __DASH_CORE_H__ */
