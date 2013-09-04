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


#ifndef __DASH_SCHED_H__
#define __DASH_SCHED_H__


#include <mpi.h>

#include "dash_conf.h"

#ifdef DS_RENAME
# include "dash_rename.h"
#endif

#include "z_pack.h"


struct _ds_t;
struct _ds_sched_t;

typedef dsint_t (*ds_sched_std_f)(struct _ds_sched_t *sched);
typedef dsint_t (*ds_sched_max_n_f)(struct _ds_sched_t *sched, dsint_t *max_nsends, dsint_t *max_nrecvs, dsint_t *max_nlocals, dsint_t *max_nsyms);

typedef struct _ds_sched_t
{
  struct _ds_t *ds;

  dsint_t round;

  dsint_t nbufs;
  struct {
    dsint_t addr_id;
  } bufs[DASH_MAX_NBUFFERS];

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


dsint_t ds_sched_create(ds_sched_t *sched);
dsint_t ds_sched_destroy(ds_sched_t *sched);

dsint_t ds_sched_add_buffer(ds_sched_t *sched, dsint_t addr_id);
dsint_t ds_sched_set_send(ds_sched_t *sched, dsint_t buf_id, dsint_t n, dspint_t *scounts, dspint_t *sdispls, dspint_t *sranks, dsint_t exec_id);
dsint_t ds_sched_set_recv(ds_sched_t *sched, dsint_t buf_id, dsint_t n, dspint_t *rcounts, dspint_t *rdispls, dspint_t *rranks, dsint_t exec_id);
dsint_t ds_sched_set_aux(ds_sched_t *sched, dsint_t buf_id, dsint_t buf_size);


#endif /* __DASH_SCHED_H__ */
