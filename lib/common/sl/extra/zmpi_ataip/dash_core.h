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

#include "dash_common.h"
#include "dash_exec.h"
#include "dash_sched.h"


extern dsint_t ds_core_run_sync;

#define DS_TIMES_RUN                 0
#define DS_TIMES_RUN_PRE             1
#define DS_TIMES_RUN_WHILE           2
#define DS_TIMES_RUN_WHILE_PRE       3
#define DS_TIMES_RUN_WHILE_MAKE      4
#define DS_TIMES_RUN_WHILE_POST      5
#define DS_TIMES_RUN_WHILE_FINISHED  6
#define DS_TIMES_RUN_WHILE_ROUNDS    7
#define DS_TIMES_RUN_POST            8
#define DS_TIMES_RUN_NTIMES          9

extern double ds_times[DS_TIMES_RUN_NTIMES];


typedef struct _ds_t
{
  ds_sched_t *sched;
  ds_exec_t *exec;
  
  MPI_Comm comm;
  int comm_size, comm_rank;

  struct {
    dsint_t nmax, n, execute, *buf_ids, exec_id;
    dspint_t *counts, *displs, *proc_ids;
    int *ranks;
  } sends;

  struct {
    dsint_t nmax, n, execute, *buf_ids, exec_id;
    dspint_t *counts, *displs, *proc_ids;
    int *ranks;
  } recvs;

  struct {
    dsint_t nmax, n, execute;
    dsint_t *src_buf_ids, *dst_buf_ids, *src_exec_ids, *dst_exec_ids;
    dspint_t *src_counts, *src_displs, *dst_counts, *dst_displs;
  } locals;

#ifdef DASH_SYMMETRIC
  struct {
    dsint_t nmax, n, execute, *buf_ids, exec_id;
    dspint_t *counts, *displs, *proc_ids;
    int *ranks;
  } syms;
#endif

} ds_t, *ds_p;


dsint_t ds_create(ds_t *ds, ds_sched_t *sched, ds_exec_t *exec);
dsint_t ds_destroy(ds_t *ds);
dsint_t ds_run(ds_t *ds);


#endif /* __DASH_CORE_H__ */
