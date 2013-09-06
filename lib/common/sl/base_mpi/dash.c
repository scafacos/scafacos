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


#include "dash_sched_a2av_aux.h"
#include "local_generic_heap.h"


#ifndef DS_TRACE_IF
# define DS_TRACE_IF  (z_mpi_rank == -1)
#endif

#define MAP_PROC_TO_BLOCK(_p_, _n_)  (_p_ % _n_)


typedef struct _aux_heap_t
{
  lgh_t lgh;

  dspint_t nblocks;
  struct {
    lgh_segment_t *seg;
    dsint_t begin, end;

  } *blocks;

} aux_heap_t;


static dsint_t ds_aux_heap_acquire(ds_sched_a2av_aux_t *aux, dsint_t count, dspint_t proc_id);
static dsint_t ds_aux_heap_get_count(ds_sched_a2av_aux_t *aux, dspint_t proc_id);
static dsint_t ds_aux_heap_get_displ(ds_sched_a2av_aux_t *aux, dspint_t proc_id);
static void ds_aux_heap_vacate(ds_sched_a2av_aux_t *aux, void *schedptr);
static void ds_aux_heap_accept_recv(ds_sched_a2av_aux_t *aux, dspint_t proc_id, dsint_t count);


void ds_aux_heap_create(ds_sched_a2av_aux_t *aux, dsint_t size, dsint_t nblocks) /* ds_func ds_aux_heap_create */
{
  aux_heap_t *d;
  dsint_t i;


  ds_aux_create(aux, size);

  aux->acquire = ds_aux_heap_acquire;

  aux->get_count = ds_aux_heap_get_count;
  aux->get_displ = ds_aux_heap_get_displ;

  aux->vacate = ds_aux_heap_vacate;

  aux->accept_recv = ds_aux_heap_accept_recv;

  d = z_alloc(1, sizeof(aux_heap_t));

  aux->data = d;

  lgh_create(&d->lgh, size);

  d->nblocks = nblocks;

  Z_TRACE_IF(DS_TRACE_IF, "nblocks: %" dsint_fmt, nblocks);

  d->blocks = z_alloc(d->nblocks, sizeof(*d->blocks));

  for (i = 0; i < d->nblocks; ++i)
  {
    d->blocks[i].seg = NULL;
    d->blocks[i].begin = d->blocks[i].end = -1;
  }
}


void ds_aux_heap_destroy(ds_sched_a2av_aux_t *aux) /* ds_func ds_aux_heap_destroy */
{
  aux_heap_t *d = (aux_heap_t *) aux->data;

  lgh_destroy(&d->lgh);

  z_free(d->blocks);

  z_free(aux->data);

  ds_aux_destroy(aux);
}


static dsint_t ds_aux_heap_acquire(ds_sched_a2av_aux_t *aux, dsint_t count, dspint_t proc_id)
{
  aux_heap_t *d = (aux_heap_t *) aux->data;
  dsint_t i;


  i = MAP_PROC_TO_BLOCK(proc_id, d->nblocks);

  if (d->blocks[i].begin < d->blocks[i].end) return 0;

  d->blocks[i].seg = lgh_alloc_minmax(&d->lgh, 1, count);

  if (d->blocks[i].seg == NULL) return 0;

  d->blocks[i].begin = d->blocks[i].seg->offset;
  d->blocks[i].end = d->blocks[i].seg->offset + d->blocks[i].seg->size;

  Z_TRACE_IF(DS_TRACE_IF, "acquired: %" dsint_fmt " (of %" dsint_fmt ") from %" dspint_fmt, d->blocks[i].end - d->blocks[i].begin, count, proc_id);

  return d->blocks[i].end - d->blocks[i].begin;
}


static dsint_t ds_aux_heap_get_count(ds_sched_a2av_aux_t *aux, dspint_t proc_id)
{
  aux_heap_t *d = (aux_heap_t *) aux->data;
  dsint_t i;


  i = MAP_PROC_TO_BLOCK(proc_id, d->nblocks);

  Z_TRACE_IF(DS_TRACE_IF, "count: %" dsint_fmt, d->blocks[i].end - d->blocks[i].begin);

  return d->blocks[i].end - d->blocks[i].begin;
}


static dsint_t ds_aux_heap_get_displ(ds_sched_a2av_aux_t *aux, dspint_t proc_id)
{
  aux_heap_t *d = (aux_heap_t *) aux->data;
  dsint_t i;


  i = MAP_PROC_TO_BLOCK(proc_id, d->nblocks);

  Z_TRACE_IF(DS_TRACE_IF, "displ: %" dsint_fmt, d->blocks[i].begin);

  return d->blocks[i].begin;
}


static void ds_aux_heap_vacate(ds_sched_a2av_aux_t *aux, void *schedptr)
{
  aux_heap_t *d = (aux_heap_t *) aux->data;
  dsint_t i, count;


  for (i = 0; i < d->nblocks; ++i)
  if (d->blocks[i].begin < d->blocks[i].end)
  {
    count = aux->vacate_aux(schedptr, i, d->blocks[i].end - d->blocks[i].begin, d->blocks[i].begin);

    d->blocks[i].begin += count;

    if (d->blocks[i].begin >= d->blocks[i].end)
    {
      lgh_free(&d->lgh, d->blocks[i].seg);
      d->blocks[i].seg = NULL;
    }

    Z_TRACE_IF(DS_TRACE_IF, "moved %" dsint_fmt " and %" dsint_fmt " are left%s", count, d->blocks[i].end - d->blocks[i].begin, (d->blocks[i].end >= d->blocks[i].begin)?" -> empty":"");
  }
}


static void ds_aux_heap_accept_recv(ds_sched_a2av_aux_t *aux, dspint_t proc_id, dsint_t count)
{
  aux_heap_t *d = (aux_heap_t *) aux->data;
  dsint_t i;


  i = MAP_PROC_TO_BLOCK(proc_id, d->nblocks);

  d->blocks[i].end = d->blocks[i].begin + count;

  if (count <= 0)
  {
    lgh_free(&d->lgh, d->blocks[i].seg);
    d->blocks[i].seg = NULL;
  }
}


#undef DS_TRACE_IF

#undef MAP_PROC_TO_BLOCK



#include "dash_sched_a2av_aux.h"


#ifndef DS_TRACE_IF
# define DS_TRACE_IF  (z_mpi_rank == -1)
#endif

#define MAP_PROC_TO_BLOCK(_p_, _n_)  (_p_ % _n_)


typedef struct _aux_static_t
{
  dspint_t nblocks;
  struct {
    dspint_t proc_id;
    dsint_t base, size, new;
    dsint_t begin, end;

  } *blocks;

} aux_static_t;


static dsint_t ds_aux_static_acquire(ds_sched_a2av_aux_t *aux, dsint_t count, dspint_t proc_id);
static dsint_t ds_aux_static_get_count(ds_sched_a2av_aux_t *aux, dspint_t proc_id);
static dsint_t ds_aux_static_get_displ(ds_sched_a2av_aux_t *aux, dspint_t proc_id);
static void ds_aux_static_vacate(ds_sched_a2av_aux_t *aux, void *schedptr);
static void ds_aux_static_accept_recv(ds_sched_a2av_aux_t *aux, dspint_t proc_id, dsint_t count);


void ds_aux_static_create(ds_sched_a2av_aux_t *aux, dsint_t size, dsint_t nblocks) /* ds_func ds_aux_static_create */
{
  aux_static_t *d;
  dsint_t i, base;


  ds_aux_create(aux, size);

  aux->acquire = ds_aux_static_acquire;

  aux->get_count = ds_aux_static_get_count;
  aux->get_displ = ds_aux_static_get_displ;

  aux->vacate = ds_aux_static_vacate;

  aux->accept_recv = ds_aux_static_accept_recv;

  d = z_alloc(1, sizeof(aux_static_t));

  aux->data = d;

  d->nblocks = z_minmax(1, nblocks, size);

  Z_TRACE_IF(DS_TRACE_IF, "nblocks: %" dsint_fmt, nblocks);

  d->blocks = z_alloc(d->nblocks, sizeof(*d->blocks));

  base = 0;
  for (i = 0; i < d->nblocks; ++i)
  {
    d->blocks[i].proc_id = -1;
    d->blocks[i].new = 0;
    d->blocks[i].base = base;
    d->blocks[i].size = (dsint_t) z_round((double) (size - base) / (double) (d->nblocks - i));
    d->blocks[i].begin = d->blocks[i].end = base;

    base += d->blocks[i].size;
  }
}


void ds_aux_static_destroy(ds_sched_a2av_aux_t *aux) /* ds_func ds_aux_static_destroy */
{
  aux_static_t *d = (aux_static_t *) aux->data;


  z_free(d->blocks);

  z_free(aux->data);

  ds_aux_destroy(aux);
}


static dsint_t ds_aux_static_acquire(ds_sched_a2av_aux_t *aux, dsint_t count, dspint_t proc_id)
{
  aux_static_t *d = (aux_static_t *) aux->data;
  dsint_t i;


  i = MAP_PROC_TO_BLOCK(proc_id, d->nblocks);

  if (d->blocks[i].begin < d->blocks[i].end) return 0;

  d->blocks[i].proc_id = proc_id;
  d->blocks[i].new = 1;
  d->blocks[i].begin = d->blocks[i].base;
  d->blocks[i].end = d->blocks[i].begin + z_min(d->blocks[i].size, count);

  Z_TRACE_IF(DS_TRACE_IF, "acquired: %" dsint_fmt " (of %" dsint_fmt ") from %" dspint_fmt, d->blocks[i].end - d->blocks[i].begin, count, d->blocks[i].proc_id);

  return d->blocks[i].end - d->blocks[i].begin;
}


static dsint_t ds_aux_static_get_count(ds_sched_a2av_aux_t *aux, dspint_t proc_id)
{
  aux_static_t *d = (aux_static_t *) aux->data;
  dsint_t i;


  i = MAP_PROC_TO_BLOCK(proc_id, d->nblocks);

  if (d->blocks[i].proc_id != proc_id) return 0;

  Z_TRACE_IF(DS_TRACE_IF, "count: %" dsint_fmt, d->blocks[i].end - d->blocks[i].begin);

  return d->blocks[i].end - d->blocks[i].begin;
}


static dsint_t ds_aux_static_get_displ(ds_sched_a2av_aux_t *aux, dspint_t proc_id)
{
  aux_static_t *d = (aux_static_t *) aux->data;
  dsint_t i;


  i = MAP_PROC_TO_BLOCK(proc_id, d->nblocks);

  Z_ASSERT(d->blocks[i].new);
  Z_ASSERT(d->blocks[i].proc_id == proc_id);

  Z_TRACE_IF(DS_TRACE_IF, "displ: %" dsint_fmt, d->blocks[i].begin);

  return d->blocks[i].begin;
}


static void ds_aux_static_vacate(ds_sched_a2av_aux_t *aux, void *schedptr)
{
  aux_static_t *d = (aux_static_t *) aux->data;
  dsint_t i, count;


  for (i = 0; i < d->nblocks; ++i)
  if (d->blocks[i].begin < d->blocks[i].end)
  {
    count = aux->vacate_aux(schedptr, d->blocks[i].proc_id, d->blocks[i].end - d->blocks[i].begin, d->blocks[i].begin);

    d->blocks[i].begin += count;

    Z_TRACE_IF(DS_TRACE_IF, "moved %" dsint_fmt " and %" dsint_fmt " are left%s", count, d->blocks[i].end - d->blocks[i].begin, (d->blocks[i].end >= d->blocks[i].begin)?" -> empty":"");
  }
}


static void ds_aux_static_accept_recv(ds_sched_a2av_aux_t *aux, dspint_t proc_id, dsint_t count)
{
  aux_static_t *d = (aux_static_t *) aux->data;
  dsint_t i;


  i = MAP_PROC_TO_BLOCK(proc_id, d->nblocks);

  Z_ASSERT(d->blocks[i].new);
  Z_ASSERT(d->blocks[i].proc_id == proc_id);


  d->blocks[i].new = 0;
  d->blocks[i].end = d->blocks[i].begin + count;
}


#undef DS_TRACE_IF

#undef MAP_PROC_TO_BLOCK



#include "dash_core.h"
#include "prx.h"


/*#define QSORT*/

#ifdef QSORT
static int compare_dsints(const void *b0, const void *b1)
{
  return *((dsint_t *) b0) - *((dsint_t *) b1);
}
#endif


void ds_sort_dsints(dsint_t *ints, dsint_t n, dsint_t x) /* ds_func ds_sort_dsints */
{
#ifdef QSORT

  qsort(ints, n, x * sizeof(dsint_t), compare_dsints);

#else

  dsint_t k, t;
#define GET(_i_)        (ints[(_i_) * x])
#define XCHG(_i_, _j_)  Z_MOP(for (k = 0; k < x; ++k) { t = ints[(_i_) * x + k]; ints[(_i_) * x + k] = ints[(_j_) * x + k]; ints[(_j_) * x + k] = t; })

  DASH_RADIXSORT_SORT_DECLARE(dsint_t)

  DASH_RADIXSORT_SORT(n, 0, sizeof(dsint_t) * 8 - 1, GET, XCHG);

#undef GET
#undef XCHG

#endif
}

#undef QSORT



#include "dash_core.h"


#ifndef DS_TRACE_IF
# define DS_TRACE_IF  (z_mpi_rank == -1)
#endif


dsint_t ds_core_run_sync = 0; /* ds_var ds_core_run_sync */

double ds_times[DS_TIMES_RUN_NTIMES];  /* ds_var ds_times */


dsint_t ds_create(ds_t *ds, ds_sched_t *sched, ds_exec_t *exec) /* ds_func ds_create */
{
  dsint_t max_nsends, max_nrecvs, max_nlocals, max_nsyms;

  Z_TRACE_IF(DS_TRACE_IF, "START");

  ds->sched = sched; sched->ds = ds;
  ds->exec = exec; exec->ds = ds;

  sched->max_n(sched, &max_nsends, &max_nrecvs, &max_nlocals, &max_nsyms);

  ds->sends.nmax = max_nsends;
  ds->sends.n = 0;
  ds->sends.execute = 0;
  ds->sends.buf_ids = z_alloc(max_nsends, sizeof(dsint_t));
  ds->sends.counts   = z_alloc(max_nsends * 3, sizeof(dspint_t));
  ds->sends.displs   = ds->sends.counts + 1 * max_nsends;
  ds->sends.proc_ids = ds->sends.counts + 2 * max_nsends;
  ds->sends.ranks = NULL;

  ds->recvs.nmax = max_nrecvs;
  ds->recvs.n = 0;
  ds->recvs.execute = 0;
  ds->recvs.buf_ids = z_alloc(max_nrecvs, sizeof(dsint_t));
  ds->recvs.counts   = z_alloc(max_nrecvs * 3, sizeof(dspint_t));
  ds->recvs.displs   = ds->recvs.counts + 1 * max_nrecvs;
  ds->recvs.proc_ids = ds->recvs.counts + 2 * max_nrecvs;
  ds->recvs.ranks = NULL;

  ds->locals.nmax = max_nlocals;
  ds->locals.n = 0;
  ds->locals.execute = 0;
  ds->locals.src_buf_ids  = z_alloc(max_nlocals * 4, sizeof(dsint_t));
  ds->locals.dst_buf_ids  = ds->locals.src_buf_ids + 1 * max_nlocals;
  ds->locals.src_exec_ids = ds->locals.src_buf_ids + 2 * max_nlocals;
  ds->locals.dst_exec_ids = ds->locals.src_buf_ids + 3 * max_nlocals;
  ds->locals.src_counts = z_alloc(max_nlocals * 4, sizeof(dspint_t));
  ds->locals.src_displs = ds->locals.src_counts + 1 * max_nlocals;
  ds->locals.dst_counts = ds->locals.src_counts + 2 * max_nlocals;
  ds->locals.dst_displs = ds->locals.src_counts + 3 * max_nlocals;

#ifdef DASH_SYMMETRIC
  ds->syms.nmax = max_nsyms;
  ds->syms.n = 0;
  ds->syms.execute = 0;
  ds->syms.buf_ids = z_alloc(max_nsyms, sizeof(dsint_t));
  ds->syms.counts   = z_alloc(max_nsyms * 3, sizeof(dspint_t));
  ds->syms.displs   = ds->syms.counts + 1 * max_nsyms;
  ds->syms.proc_ids = ds->syms.counts + 2 * max_nsyms;
  ds->syms.ranks = NULL;
#endif

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


dsint_t ds_destroy(ds_t *ds) /* ds_func ds_destroy */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  ds->sched = NULL;
  ds->exec = NULL;

  z_free(ds->sends.buf_ids);
  ds->sends.buf_ids = NULL;
  z_free(ds->sends.counts);
  ds->sends.counts = ds->sends.displs = ds->sends.proc_ids = NULL;
  ds->sends.ranks = NULL;

  z_free(ds->recvs.buf_ids);
  ds->recvs.buf_ids = NULL;
  z_free(ds->recvs.counts);
  ds->recvs.counts = ds->recvs.displs = ds->recvs.proc_ids = NULL;
  ds->recvs.ranks = NULL;

  z_free(ds->locals.src_buf_ids);
  ds->locals.src_buf_ids = ds->locals.dst_buf_ids = ds->locals.src_exec_ids = ds->locals.dst_exec_ids = NULL;
  z_free(ds->locals.src_counts);
  ds->locals.src_counts = ds->locals.src_displs = ds->locals.dst_counts = ds->locals.dst_displs = NULL;

#ifdef DASH_SYMMETRIC
  z_free(ds->syms.buf_ids);
  ds->syms.buf_ids = NULL;
  z_free(ds->syms.counts);
  ds->syms.counts = ds->syms.displs = ds->syms.proc_ids = NULL;
  ds->syms.ranks = NULL;
#endif

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


static void ds_print(ds_t *ds)
{
  Z_TRACE_DECLARE(dsint_t i;)


  Z_TRACE_IF(DS_TRACE_IF, "sends: execute: %" dsint_fmt ", n: %" dsint_fmt " / %" dsint_fmt, ds->sends.execute, ds->sends.n, ds->sends.nmax);

  if (ds->sends.execute)
  {
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->sends.n, " %" dspint_fmt, ds->sends.counts[i], "send counts:");
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->sends.n, " %" dspint_fmt, ds->sends.displs[i], "send displs:");
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->sends.n, " %" dspint_fmt, ds->sends.proc_ids[i], "send proc_ids:");
  }

  Z_TRACE_IF(DS_TRACE_IF, "receives: execute: %" dsint_fmt ", n: %" dsint_fmt " / %" dsint_fmt, ds->recvs.execute, ds->recvs.n, ds->recvs.nmax);

  if (ds->recvs.execute)
  {
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->recvs.n, " %" dspint_fmt, ds->recvs.counts[i], "recv counts:");
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->recvs.n, " %" dspint_fmt, ds->recvs.displs[i], "recv displs:");
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->recvs.n, " %" dspint_fmt, ds->recvs.proc_ids[i], "recv proc_ids:");
  }

  Z_TRACE_IF(DS_TRACE_IF, "locals: execute: %" dsint_fmt ", n: %" dsint_fmt " / %" dsint_fmt, ds->locals.execute, ds->locals.n, ds->locals.nmax);

  if (ds->locals.execute)
  {
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->locals.n, " %" dspint_fmt, ds->locals.src_counts[i], "local src counts:");
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->locals.n, " %" dspint_fmt, ds->locals.src_displs[i], "local src displs:");
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->locals.n, " %" dspint_fmt, ds->locals.dst_counts[i], "local dst counts:");
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->locals.n, " %" dspint_fmt, ds->locals.dst_displs[i], "local dst displs:");
  }

#ifdef DASH_SYMMETRIC
  Z_TRACE_IF(DS_TRACE_IF, "syms: execute: %" dsint_fmt ", n: %" dsint_fmt " / %" dsint_fmt, ds->syms.execute, ds->syms.n, ds->syms.nmax);

  if (ds->syms.execute)
  {
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->syms.nmax, " %" dspint_fmt, ds->syms.counts[i], "sym counts:");
    Z_TRACE_ARRAY_IF(DS_TRACE_IF, i, ds->syms.nmax, " %" dspint_fmt, ds->syms.displs[i], "sym displs:");
  }
#endif
}


dsint_t ds_run(ds_t *ds) /* ds_func ds_run */
{
#ifdef DASH_TIMING
  double t;
#endif


  Z_TRACE_IF(DS_TRACE_IF, "START");

  if (ds_core_run_sync) MPI_Barrier(ds->comm);
  DS_TIMING_CMD(ds_times[DS_TIMES_RUN] = z_time_get_s(););
  DS_TIMING_CMD(ds_times[DS_TIMES_RUN_PRE] = z_time_get_s(););

  if (ds->sched->pre_run) ds->sched->pre_run(ds->sched);
  ds->exec->pre_run(ds->exec);

  if (ds_core_run_sync) MPI_Barrier(ds->comm);
  DS_TIMING_CMD(ds_times[DS_TIMES_RUN_PRE] = z_time_get_s() - ds_times[DS_TIMES_RUN_PRE];);
  DS_TIMING_CMD(ds_times[DS_TIMES_RUN_WHILE] = z_time_get_s(););

  DS_TIMING_CMD(ds_times[DS_TIMES_RUN_WHILE_PRE] = ds_times[DS_TIMES_RUN_WHILE_MAKE] = ds_times[DS_TIMES_RUN_WHILE_POST] = ds_times[DS_TIMES_RUN_WHILE_FINISHED] = 0;);

  DS_TIMING_CMD(t = z_time_get_s(););

  while (!ds->sched->finished(ds->sched))
  {
    if (ds_core_run_sync) MPI_Barrier(ds->comm);
    DS_TIMING_CMD(ds_times[DS_TIMES_RUN_WHILE_FINISHED] += z_time_get_s() - t; t = z_time_get_s(););

    if (ds->sched->pre) ds->sched->pre(ds->sched);

    if (ds_core_run_sync) MPI_Barrier(ds->comm);
    DS_TIMING_CMD(ds_times[DS_TIMES_RUN_WHILE_PRE] += z_time_get_s() - t; t = z_time_get_s(););

    ds_print(ds);

    ds_exec_make(ds->exec);

    if (ds_core_run_sync) MPI_Barrier(ds->comm);
    DS_TIMING_CMD(ds_times[DS_TIMES_RUN_WHILE_MAKE] += z_time_get_s() - t; t = z_time_get_s(););

    if (ds->sched->post) ds->sched->post(ds->sched);

    if (ds_core_run_sync) MPI_Barrier(ds->comm);
    DS_TIMING_CMD(ds_times[DS_TIMES_RUN_WHILE_POST] += z_time_get_s() - t; t = z_time_get_s(););
  }

  DS_TIMING_CMD(
    ds_times[DS_TIMES_RUN_WHILE_FINISHED] += z_time_get_s() - t;

    ds_times[DS_TIMES_RUN_WHILE_ROUNDS] = ds->sched->round;

    ds_times[DS_TIMES_RUN_WHILE] = z_time_get_s() - ds_times[DS_TIMES_RUN_WHILE];
    ds_times[DS_TIMES_RUN_POST] = z_time_get_s();
  );

  if (ds->sched->post_run) ds->sched->post_run(ds->sched);
  ds->exec->post_run(ds->exec);

  if (ds_core_run_sync) MPI_Barrier(ds->comm);
  DS_TIMING_CMD(
    ds_times[DS_TIMES_RUN_POST] = z_time_get_s() - ds_times[DS_TIMES_RUN_POST];
    ds_times[DS_TIMES_RUN] = z_time_get_s() - ds_times[DS_TIMES_RUN];
  );

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


#undef DS_TRACE_IF



#include "dash_exec.h"
#include "dash_core.h"
#include "prx.h"


#ifndef DS_TRACE_IF
# define DS_TRACE_IF  (z_mpi_rank == -1)
#endif


#define PROC_ID2RANK(_id_, _rs_)  ((int) ((_rs_)?(_rs_)[_id_]:(_id_)))


dsint_t ds_exec_create(ds_exec_t *exec) /* ds_func ds_exec_create */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  exec->naddrs = 0;

  exec->cxt = NULL;

  exec->pre_run = NULL;
  exec->post_run = NULL;
  exec->make = NULL;
  exec->move = NULL;

#ifdef DASH_SYMMETRIC
  exec->make_sym = DASH_EXEC_MAKE_SYM_DEFAULT;
  exec->sendrecv_replace = NULL;
# ifdef DASH_SYMMETRIC_AUX
  exec->sendrecv_aux_setup = NULL;
  exec->sendrecv_aux = NULL;
  exec->sendrecv_aux_intermediate = NULL;
  exec->sendrecv_aux_finish = NULL;
# endif
#endif

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


dsint_t ds_exec_destroy(ds_exec_t *exec) /* ds_func ds_exec_destroy */
{
  return 0;
}


dsint_t ds_exec_add_address(ds_exec_t *exec, void *addr) /* ds_func ds_exec_add_address */
{
  dsint_t addr_id = -1;

  Z_TRACE_IF(DS_TRACE_IF, "START");

  if (exec->naddrs < DASH_MAX_NBUFFERS)
  {
    addr_id = exec->naddrs;
    ++exec->naddrs;

    exec->addrs[addr_id] = addr;
  }

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return addr_id;
}


#ifdef DASH_SYMMETRIC

#ifdef DASH_SYMMETRIC_AUX


typedef struct _sym_aux_t
{
  dsint_t offset;

} sym_aux_t;


static void sym_aux_setup(ds_exec_t *exec, sym_aux_t *sym_aux)
{
  if (exec->sendrecv_aux_setup) exec->sendrecv_aux_setup(exec, exec->ds->comm_size);

  sym_aux->offset = 0;
}


static void sym_aux_intermediate(ds_exec_t *exec, sym_aux_t *sym_aux)
{
  if (exec->sendrecv_aux_intermediate) exec->sendrecv_aux_intermediate(exec);

  sym_aux->offset = 0;
}


static dsint_t sym_aux_sendrecv(ds_exec_t *exec, sym_aux_t *sym_aux, int proc)
{
  if (sym_aux->offset + exec->ds->syms.counts[proc] > exec->ds->sched->aux.buf_size) return 0;

  if (!exec->sendrecv_aux(exec, proc, exec->ds->syms.exec_id, exec->ds->syms.buf_ids[proc], exec->ds->syms.displs[proc], exec->ds->syms.counts[proc], exec->ds->syms.exec_id, exec->ds->sched->aux.buf_id, sym_aux->offset))
  {
    sym_aux_intermediate(exec, sym_aux);

    return sym_aux_sendrecv(exec, sym_aux, proc);
  }

  sym_aux->offset += exec->ds->syms.counts[proc];

  return 1;
}


static void sym_aux_finish(ds_exec_t *exec, sym_aux_t *sym_aux)
{
  sym_aux_intermediate(exec, sym_aux);

  if (exec->sendrecv_aux_finish) exec->sendrecv_aux_finish(exec);
}


#endif /* DASH_SYMMETRIC_AUX */


#define SYM_SPARSE_THRESHOLD  0.5

static void sym_prepare_full(ds_exec_t *exec)
{
  dsint_t i, j, t_dsint;
  dspint_t t_dspint;


  Z_TRACE_IF(DS_TRACE_IF, "original syms (%" dsint_fmt ")", exec->ds->syms.n);
  for (i = 0; i < exec->ds->syms.n; ++i)
    Z_TRACE_IF(DS_TRACE_IF, "  %" dsint_fmt ": rank: %d, count: %" dspint_fmt ", displ: %" dspint_fmt, i, (exec->ds->syms.proc_ids[i] >= 0)?PROC_ID2RANK(exec->ds->syms.proc_ids[i], exec->ds->syms.ranks):-1, exec->ds->syms.counts[i], exec->ds->syms.displs[i]);

  Z_TRACE_IF(DS_TRACE_IF, "mode: %s", (exec->ds->syms.nmax < exec->ds->comm_size || exec->ds->syms.n < SYM_SPARSE_THRESHOLD * exec->ds->comm_size)?"sparse":"full");

  /* skip full mode if nmax too small (i.e., sparse!!!) or only "few" sym-exchanges in comparison to comm_size (i.e., sorting the few sym-exchanges is faster than full iteration) */
  if (exec->ds->syms.nmax < exec->ds->comm_size || exec->ds->syms.n < SYM_SPARSE_THRESHOLD * exec->ds->comm_size) return;

  /* fill up with empty sym-exchanges */
  while (exec->ds->syms.n < exec->ds->comm_size)
  {
    exec->ds->syms.buf_ids[exec->ds->syms.n] = 0;
    exec->ds->syms.counts[exec->ds->syms.n] = 0;
    exec->ds->syms.displs[exec->ds->syms.n] = 0;
    exec->ds->syms.proc_ids[exec->ds->syms.n] = -1;

    ++exec->ds->syms.n;
  }

  Z_TRACE_IF(DS_TRACE_IF, "expanded syms (%" dsint_fmt ")", exec->ds->syms.n);
  for (i = 0; i < exec->ds->syms.n; ++i)
    Z_TRACE_IF(DS_TRACE_IF, "  %" dsint_fmt ": rank: %d, count: %" dspint_fmt ", displ: %" dspint_fmt, i, (exec->ds->syms.proc_ids[i] >= 0)?PROC_ID2RANK(exec->ds->syms.proc_ids[i], exec->ds->syms.ranks):-1, exec->ds->syms.counts[i], exec->ds->syms.displs[i]);

  /* permute sym-exchanges */
  if (exec->ds->syms.ranks)
  {
    for (i = 0; i < exec->ds->syms.n; ++i)
    {
      j = exec->ds->syms.proc_ids[i];

      if (j >= 0)
      {
        j = exec->ds->syms.ranks[j];

        while (j != i)
        {
          z_swap(exec->ds->syms.buf_ids[i], exec->ds->syms.buf_ids[j], t_dsint);
          z_swap(exec->ds->syms.counts[i], exec->ds->syms.counts[j], t_dspint);
          z_swap(exec->ds->syms.displs[i], exec->ds->syms.displs[j], t_dspint);
          z_swap(exec->ds->syms.proc_ids[i], exec->ds->syms.proc_ids[j], t_dspint);

          /* break if j < i since swapping entries i and j with j < i always makes entry i empty (because all entries less than i are either done (i.e., no swap) or empty) */
          if (j < i) break;

          j = exec->ds->syms.proc_ids[i];

          /* break on empty entry */
          if (j < 0) break;

          j = exec->ds->syms.ranks[j];
        }
      }

      /* override the -1 empty markers */
      exec->ds->syms.proc_ids[i] = i;
    }

    exec->ds->syms.ranks = NULL;

  } else
  {
    for (i = 0; i < exec->ds->syms.n; ++i)
    {
      j = exec->ds->syms.proc_ids[i];

      if (j >= 0)
      {
        while (j != i)
        {
          z_swap(exec->ds->syms.buf_ids[i], exec->ds->syms.buf_ids[j], t_dsint);
          z_swap(exec->ds->syms.counts[i], exec->ds->syms.counts[j], t_dspint);
          z_swap(exec->ds->syms.displs[i], exec->ds->syms.displs[j], t_dspint);
          z_swap(exec->ds->syms.proc_ids[i], exec->ds->syms.proc_ids[j], t_dspint);

          /* break if j < i since swapping entries i and j with j < i always makes entry i empty (because all entries less than i are either done (i.e., no swap) or empty) */
          if (j < i) break;

          j = exec->ds->syms.proc_ids[i];

          /* break on empty entry */
          if (j < 0) break;
        }
      }

      /* override the -1 empty markers */
      exec->ds->syms.proc_ids[i] = i;
    }
  }

  Z_TRACE_IF(DS_TRACE_IF, "permuted syms (%" dsint_fmt ")", exec->ds->syms.n);
  for (i = 0; i < exec->ds->syms.n; ++i)
    Z_TRACE_IF(DS_TRACE_IF, "  %" dsint_fmt ": rank: %d, count: %" dspint_fmt ", displ: %" dspint_fmt, i, (exec->ds->syms.proc_ids[i] >= 0)?PROC_ID2RANK(exec->ds->syms.proc_ids[i], exec->ds->syms.ranks):-1, exec->ds->syms.counts[i], exec->ds->syms.displs[i]);
}


static void sym_order_destroy(dsint_t *proc_ids)
{
  z_free(proc_ids);
}


static dsint_t *sym_order_linear_create(ds_exec_t *exec, dsint_t randomize_ranks)
{
  dsint_t i, t;
  dsint_t *order;
  prx_enumerate_t random_enumerate = PRX_ENUMERATE_NULL;
  DASH_RADIXSORT_SORT_DECLARE(dsint_t)
  dsint_t rhigh;


  order = z_alloc(exec->ds->syms.n, sizeof(dsint_t));

  if (exec->ds->syms.n >= exec->ds->comm_size)
  {
    /* full mode (i.e., stored rank-wise) */

    if (randomize_ranks)
    {
      prx_seed(2501);
      prx_enumerate_create(&random_enumerate, exec->ds->comm_size, PRX_FISHER_YATES_SHUFFLE);

      for (i = 0; i < exec->ds->syms.n; ++i)
      {
#define GET(_i_)  ((int) prx_enumerate(random_enumerate, (exec->ds->comm_rank + (_i_)) % exec->ds->comm_size))
        order[GET(i)] = i;
#undef GET
      }

      prx_enumerate_destroy(&random_enumerate);

    } else
    {
      for (i = 0; i < exec->ds->syms.n; ++i)
      {
#define GET(_i_)  ((int) ((exec->ds->comm_rank + (_i_)) % exec->ds->comm_size))
        order[GET(i)] = i;
#undef GET
      }
    }

  } else
  {
    /* sparse mode (i.e., need to sort) */

    for (i = 0; i < exec->ds->syms.n; ++i) order[i] = i;

#define XCHG(_i_, _j_)  Z_MOP(t = order[_i_]; order[_i_] = order[_j_]; order[_j_] = t;)

    if (randomize_ranks)
    {
      prx_seed(2501);
      prx_enumerate_create(&random_enumerate, exec->ds->comm_size, PRX_FISHER_YATES_SHUFFLE);

      rhigh = (dsint_t) floor(log(exec->ds->comm_size) / log(2));

      if (exec->ds->syms.ranks)
      {
#define GET(_i_)  ((int) prx_enumerate(random_enumerate, (exec->ds->comm_rank + exec->ds->syms.ranks[exec->ds->syms.proc_ids[order[_i_]]]) % exec->ds->comm_size))
        DASH_RADIXSORT_SORT(exec->ds->syms.n, 0, rhigh, GET, XCHG);
#undef GET
      } else
      {
#define GET(_i_)  ((int) prx_enumerate(random_enumerate, (exec->ds->comm_rank + exec->ds->syms.proc_ids[order[_i_]]) % exec->ds->comm_size))
        DASH_RADIXSORT_SORT(exec->ds->syms.n, 0, rhigh, GET, XCHG);
#undef GET
      }

      prx_enumerate_destroy(&random_enumerate);

    } else
    {
      rhigh = (dsint_t) floor(log(exec->ds->comm_size) / log(2));

      if (exec->ds->syms.ranks)
      {
#define GET(_i_)  ((int) ((exec->ds->comm_rank + exec->ds->syms.ranks[exec->ds->syms.proc_ids[order[_i_]]]) % exec->ds->comm_size))
        DASH_RADIXSORT_SORT(exec->ds->syms.n, 0, rhigh, GET, XCHG);
#undef GET
      } else
      {
#define GET(_i_)  ((int) ((exec->ds->comm_rank + exec->ds->syms.proc_ids[order[_i_]]) % exec->ds->comm_size))
        DASH_RADIXSORT_SORT(exec->ds->syms.n, 0, rhigh, GET, XCHG);
#undef GET
      }
    }

#undef XCHG
  }

  Z_TRACE_IF(DS_TRACE_IF, "linear order:");
  for (i = 0; i < exec->ds->syms.n; ++i)
    Z_TRACE_IF(DS_TRACE_IF, "  %" dsint_fmt " -> %" dsint_fmt " -> %d", i, order[i], (order[i] >= 0)?PROC_ID2RANK(exec->ds->syms.proc_ids[order[i]], exec->ds->syms.ranks):-1);

  return order;
}


static void make_sym_linear(ds_exec_t *exec, dsint_t randomize)
{
  dsint_t i, j;
  int r;
  dsint_t *order;
#ifdef DASH_SYMMETRIC_AUX
  sym_aux_t sym_aux;
  dsint_t sym_aux_success;
#endif


  order = sym_order_linear_create(exec, randomize);

#ifdef DASH_SYMMETRIC_AUX
  sym_aux_setup(exec, &sym_aux);
#endif

  for (i = 0; i < exec->ds->syms.n; ++i)
  {
    j = order[i];

    if (j < 0) continue;

    r = PROC_ID2RANK(exec->ds->syms.proc_ids[j], exec->ds->syms.ranks);

    Z_TRACE_IF(DS_TRACE_IF, "sendrecv_replace with %d", r);

    if (exec->ds->comm_rank != r && exec->ds->syms.counts[j] > 0)
    {
#ifdef DASH_SYMMETRIC_AUX
      sym_aux_success = sym_aux_sendrecv(exec, &sym_aux, j);
      if (!sym_aux_success)
      {
#endif
        exec->sendrecv_replace(exec, r, exec->ds->syms.exec_id, exec->ds->syms.buf_ids[j], exec->ds->syms.displs[j], exec->ds->syms.counts[j]);
#ifdef DASH_SYMMETRIC_AUX
# ifdef DASH_SYMMETRIC_AUX_IMMEDIATELY
        sym_aux_intermediate(exec, &sym_aux);
# endif
      }
#endif
    }
  }

#ifdef DASH_SYMMETRIC_AUX
  sym_aux_finish(exec, &sym_aux);
#endif

  sym_order_destroy(order);
}


#ifdef OLD_HIERARCHIC


static void make_sym_hierarchic(ds_exec_t *exec, dsint_t randomize)
{
  int low, high, mid, skip, i, j, n;
#ifdef DASH_SYMMETRIC_AUX
  sym_aux_t sym_aux;
  dsint_t sym_aux_success;
#endif
  prx_enumerate_t random_enumerate = PRX_ENUMERATE_NULL;


  low = 0;
  high = exec->ds->comm_size;

#ifdef DASH_SYMMETRIC_AUX
  sym_aux_setup(exec, &sym_aux);
#endif

  while (low + 1 < high)
  {
    mid = (low + high) / 2;

#if 0

    if (exec->ds->comm_rank < mid)
    {
      n = high - mid;

      for (i = 0; i < n; ++i)
      {
        skip = z_min(exec->ds->comm_rank - low, n - 1);
        j = mid + (i + skip) % n;

        if (exec->ds->syms.counts[j] > 0)
        {
#ifdef DASH_SYMMETRIC_AUX
          sym_aux_success = sym_aux_sendrecv(exec, &sym_aux, j);
          if (!sym_aux_success)
          {
#endif
            exec->sendrecv_replace(exec, j, exec->ds->syms.exec_id, exec->ds->syms.buf_ids[j], exec->ds->syms.displs[j], exec->ds->syms.counts[j]);
#ifdef DASH_SYMMETRIC_AUX
# ifdef DASH_SYMMETRIC_AUX_IMMEDIATELY
            sym_aux_intermediate(exec, &sym_aux);
# endif
          }
#endif
        }
      }

      high = mid;

    } else
    {
      n = mid - low;

      for (i = 0; i < n; ++i)
      {
        skip = z_min(exec->ds->comm_rank - mid, n - 1);
        j = low + (skip - i + n) % n;

        if (exec->ds->syms.counts[j] > 0)
        {
#ifdef DASH_SYMMETRIC_AUX
          sym_aux_success = sym_aux_sendrecv(exec, &sym_aux, j);
          if (!sym_aux_success)
          {
#endif
            exec->sendrecv_replace(exec, j, exec->ds->syms.exec_id, exec->ds->syms.buf_ids[j], exec->ds->syms.displs[j], exec->ds->syms.counts[j]);
#ifdef DASH_SYMMETRIC_AUX
# ifdef DASH_SYMMETRIC_AUX_IMMEDIATELY
            sym_aux_intermediate(exec, &sym_aux);
# endif
          }
#endif
        }
      }

      low = mid;
    }

#else

    n = z_max(mid - low, high - mid);

    if (randomize)
    {
      prx_seed(2501);
      prx_enumerate_create(&random_enumerate, n, PRX_FISHER_YATES_SHUFFLE);
    }

    if (exec->ds->comm_rank < mid)
    {
      for (i = 0; i < n; ++i)
      {
        skip = exec->ds->comm_rank - low;
        j = mid + ((randomize?prx_enumerate(random_enumerate, i):i) + skip) % n;

        if (j >= mid && exec->ds->syms.counts[j] > 0)
        {
#ifdef DASH_SYMMETRIC_AUX
          sym_aux_success = sym_aux_sendrecv(exec, &sym_aux, j);
          if (!sym_aux_success)
          {
#endif
            exec->sendrecv_replace(exec, j, exec->ds->syms.exec_id, exec->ds->syms.buf_ids[j], exec->ds->syms.displs[j], exec->ds->syms.counts[j]);
#ifdef DASH_SYMMETRIC_AUX
# ifdef DASH_SYMMETRIC_AUX_IMMEDIATELY
            sym_aux_intermediate(exec, &sym_aux);
# endif
          }
#endif
        }
      }

      high = mid;

    } else
    {
      for (i = 0; i < n; ++i)
      {
        skip = exec->ds->comm_rank - mid;
        j = low + (skip - (randomize?prx_enumerate(random_enumerate, i):i) + n) % n;

        if (j < mid && exec->ds->syms.counts[j] > 0)
        {
#ifdef DASH_SYMMETRIC_AUX
          sym_aux_success = sym_aux_sendrecv(exec, &sym_aux, j);
          if (!sym_aux_success)
          {
#endif
            exec->sendrecv_replace(exec, j, exec->ds->syms.exec_id, exec->ds->syms.buf_ids[j], exec->ds->syms.displs[j], exec->ds->syms.counts[j]);
#ifdef DASH_SYMMETRIC_AUX
# ifdef DASH_SYMMETRIC_AUX_IMMEDIATELY
            sym_aux_intermediate(exec, &sym_aux);
# endif
          }
#endif
        }
      }

      low = mid;
    }

    if (randomize) prx_enumerate_destroy(&random_enumerate);
#endif
  }

#ifdef DASH_SYMMETRIC_AUX
  sym_aux_finish(exec, &sym_aux);
#endif
}


#else /* OLD_HIERARCHIC */


static dsint_t *sym_order_hierarchic_create(ds_exec_t *exec, dsint_t randomize_ranks)
{
  dsint_t i, t;
  dsint_t *order;
  DASH_RADIXSORT_SORT_DECLARE(dsint_t)
  dsint_t rhigh;


  if (exec->ds->syms.n >= exec->ds->comm_size)
  {
    /* full mode (i.e., stored rank-wise, no need for separate rank order) */

    order = z_alloc(exec->ds->syms.n, sizeof(dsint_t));

  } else
  {
    /* sparse mode (i.e., need to sort to create rank order) */

    order = z_alloc(2 * exec->ds->syms.n + 2, sizeof(dsint_t));

    for (i = 0; i < exec->ds->syms.n; ++i) order[i] = i;

#define XCHG(_i_, _j_)  Z_MOP(t = order[_i_]; order[_i_] = order[_j_]; order[_j_] = t;)

    rhigh = (dsint_t) floor(log(exec->ds->comm_size) / log(2));

    if (exec->ds->syms.ranks)
    {
#define GET(_i_)  ((int) exec->ds->syms.ranks[exec->ds->syms.proc_ids[order[_i_]]])
      DASH_RADIXSORT_SORT(exec->ds->syms.n, 0, rhigh, GET, XCHG);
#undef GET
    } else
    {
#define GET(_i_)  ((int) exec->ds->syms.proc_ids[order[_i_]])
      DASH_RADIXSORT_SORT(exec->ds->syms.n, 0, rhigh, GET, XCHG);
#undef GET
    }

#undef XCHG

    /* store range for next search in select */
    order[2 * exec->ds->syms.n + 0] = 0;
    order[2 * exec->ds->syms.n + 1] = exec->ds->syms.n;
  }

  return order;
}


static dsint_t sym_order_hierarchic_search_rank_ge(ds_exec_t *exec, dsint_t *order, int rank, dsint_t range_low, dsint_t range_high)
{
  DASH_BINARY_SEARCH_DECLARE()
  dsint_t i;


  Z_TRACE_IF(DS_TRACE_IF, "search rank ge %d in [%" dsint_fmt ",%" dsint_fmt "]", rank, range_low, range_high);
  for (i = range_low; i < range_high; ++i)
    Z_TRACE_IF(DS_TRACE_IF, "  %" dsint_fmt ": %" dsint_fmt " -> %" dspint_fmt, i, order[i], exec->ds->syms.proc_ids[order[i]]);

  --range_high;

  if (exec->ds->syms.ranks)
  {
#define GET(_i_)  exec->ds->syms.ranks[exec->ds->syms.proc_ids[order[_i_]]]
    DASH_BINARY_SEARCH_GE(range_low, range_high, rank, GET, i);
#undef GET

  } else
  {
#define GET(_i_)  exec->ds->syms.proc_ids[order[_i_]]
    DASH_BINARY_SEARCH_GE(range_low, range_high, rank, GET, i);
#undef GET
  }

  Z_TRACE_IF(DS_TRACE_IF, "found %" dsint_fmt, i);

  return i;
}


static dsint_t sym_order_hierarchic_search_rank_le(ds_exec_t *exec, dsint_t *order, int rank, dsint_t range_low, dsint_t range_high)
{
  DASH_BINARY_SEARCH_DECLARE()
  dsint_t i;


  Z_TRACE_IF(DS_TRACE_IF, "search rank le %d in [%" dsint_fmt ",%" dsint_fmt "]", rank, range_low, range_high);
  for (i = range_low; i < range_high; ++i)
    Z_TRACE_IF(DS_TRACE_IF, "  %" dsint_fmt ": %" dsint_fmt " -> %" dspint_fmt, i, order[i], exec->ds->syms.proc_ids[order[i]]);

  --range_high;

  if (exec->ds->syms.ranks)
  {
#define GET(_i_)  exec->ds->syms.ranks[exec->ds->syms.proc_ids[order[_i_]]]
    DASH_BINARY_SEARCH_LE(range_low, range_high, rank, GET, i);
#undef GET

  } else
  {
#define GET(_i_)  exec->ds->syms.proc_ids[order[_i_]]
    DASH_BINARY_SEARCH_LE(range_low, range_high, rank, GET, i);
#undef GET
  }

  Z_TRACE_IF(DS_TRACE_IF, "found %" dsint_fmt, i);

  return i;
}


static dsint_t *sym_order_hierarchic_select(ds_exec_t *exec, dsint_t randomize_ranks, dsint_t *order, int low, int mid, int high, dsint_t *n)
{
  dsint_t i, t, l, h, xsize, obase, osize, mskip, f;
  dsint_t *select;
  prx_enumerate_t random_enumerate = PRX_ENUMERATE_NULL;
  DASH_RADIXSORT_SORT_DECLARE(dsint_t)
  dsint_t rhigh;


  Z_TRACE_IF(DS_TRACE_IF, "low = %d, mid = %d, high = %d", low, mid, high);

  select = order;

  if (exec->ds->syms.n >= exec->ds->comm_size)
  {
    /* full mode (i.e., stored rank-wise) */

    if (exec->ds->comm_rank < mid)
    {
      i = low;
      l = mid;
      h = high;

    } else
    {
      l = low;
      h = mid;
      i = high;
    }

  } else
  {
    /* sparse mode (i.e., have rank order) */

    select += exec->ds->syms.n;

    if (exec->ds->comm_rank < mid)
    {
      i = sym_order_hierarchic_search_rank_ge(exec, order, low, order[2 * exec->ds->syms.n + 0], order[2 * exec->ds->syms.n + 1]);
      l = sym_order_hierarchic_search_rank_ge(exec, order, mid, order[2 * exec->ds->syms.n + 0], order[2 * exec->ds->syms.n + 1]);
      h = sym_order_hierarchic_search_rank_ge(exec, order, high, order[2 * exec->ds->syms.n + 0], order[2 * exec->ds->syms.n + 1]);

      order[2 * exec->ds->syms.n + 0] = i;
      order[2 * exec->ds->syms.n + 1] = l;

    } else
    {
      l = sym_order_hierarchic_search_rank_ge(exec, order, low, order[2 * exec->ds->syms.n + 0], order[2 * exec->ds->syms.n + 1]);
      h = sym_order_hierarchic_search_rank_ge(exec, order, mid, order[2 * exec->ds->syms.n + 0], order[2 * exec->ds->syms.n + 1]);
      i = sym_order_hierarchic_search_rank_ge(exec, order, high, order[2 * exec->ds->syms.n + 0], order[2 * exec->ds->syms.n + 1]);

      order[2 * exec->ds->syms.n + 0] = h;
      order[2 * exec->ds->syms.n + 1] = i;
    }
  }

  xsize = z_max(mid - low, high - mid);

  Z_TRACE_IF(DS_TRACE_IF, "xsize = %" dsint_fmt ", l = %" dsint_fmt ", h = %" dsint_fmt ", i = %" dsint_fmt, xsize, l, h, i);

  if (exec->ds->comm_rank < mid)
  {
    mskip = exec->ds->comm_rank - low;
    obase = mid;
    osize = high - mid;
    f = 1;

    Z_TRACE_IF(DS_TRACE_IF, "i'm %d in lower -> exchange with %" dsint_fmt " higher!", exec->ds->comm_rank - low, osize);

  } else
  {
    mskip = exec->ds->comm_rank - mid;
    obase = low;
    osize = mid - low;
    f = -1;

    Z_TRACE_IF(DS_TRACE_IF, "i'm %d in higher -> exchange with %" dsint_fmt " lower!", exec->ds->comm_rank - mid, osize);
  }

  Z_TRACE_IF(DS_TRACE_IF, "mskip = %" dsint_fmt ", obase = %" dsint_fmt ", osize = %" dsint_fmt ", f = %" dsint_fmt, mskip, obase, osize, f);

  if (exec->ds->syms.n >= exec->ds->comm_size)
  {
    /* full mode (i.e., stored rank-wise) */

    *n = xsize;

    if (randomize_ranks)
    {
      prx_seed(2501);
      prx_enumerate_create(&random_enumerate, xsize, PRX_FISHER_YATES_SHUFFLE);

      for (i = 0; i < xsize; ++i)
      {
#define GET(_i_)  ((int) prx_enumerate(random_enumerate, (xsize + f * ((_i_) - mskip)) % xsize))
        select[GET(i)] = (i < osize)?(obase + i):-1;
#undef GET
      }

      prx_enumerate_destroy(&random_enumerate);

    } else
    {
      for (i = 0; i < xsize; ++i)
      {
#define GET(_i_)  ((int) ((xsize + f * ((_i_) - mskip)) % xsize))
        select[GET(i)] = (i < osize)?(obase + i):-1;
#undef GET
      }
    }

  } else
  {
    /* sparse mode (i.e., have rank order) */

    *n = h - l;

    Z_TRACE_IF(DS_TRACE_IF, "using order [%" dsint_fmt ",%" dsint_fmt "]", l, h);
    for (i = 0; i < exec->ds->syms.n; ++i)
      Z_TRACE_IF(DS_TRACE_IF, "  %" dsint_fmt ": %" dsint_fmt " -> %" dspint_fmt "", i, order[i], exec->ds->syms.proc_ids[order[i]]);

    if (randomize_ranks)
    {
      for (i = 0; i < *n; ++i) select[i] = order[l + i];

#define XCHG(_i_, _j_)  Z_MOP(t = select[_i_]; select[_i_] = select[_j_]; select[_j_] = t;)

      prx_seed(2501);
      prx_enumerate_create(&random_enumerate, xsize, PRX_FISHER_YATES_SHUFFLE);

      rhigh = (dsint_t) floor(log(xsize) / log(2));

      if (exec->ds->syms.ranks)
      {
#define GET(_i_)  ((int) prx_enumerate(random_enumerate, (xsize + f * (exec->ds->syms.ranks[exec->ds->syms.proc_ids[select[_i_]]] - obase - mskip)) % xsize))
        DASH_RADIXSORT_SORT(*n, 0, rhigh, GET, XCHG);
#undef GET
      } else
      {
#define GET(_i_)  ((int) prx_enumerate(random_enumerate, (xsize + f * (exec->ds->syms.proc_ids[select[_i_]] - obase - mskip)) % xsize))
        DASH_RADIXSORT_SORT(*n, 0, rhigh, GET, XCHG);
#undef GET
      }

      prx_enumerate_destroy(&random_enumerate);

#undef XCHG

    } else
    {
      if (f > 0) t = sym_order_hierarchic_search_rank_ge(exec, order, obase + mskip, l, h);
      else t = sym_order_hierarchic_search_rank_le(exec, order, obase + mskip, l, h);

      for (i = 0; i < *n; ++i) select[i] = order[((*n + t + i * f - l) % *n) + l];
    }
  }

  Z_TRACE_IF(DS_TRACE_IF, "hierarchic select:");
  for (i = 0; i < *n; ++i)
    Z_TRACE_IF(DS_TRACE_IF, "  %" dsint_fmt " -> %" dsint_fmt " -> %d", i, select[i], (select[i] >= 0)?PROC_ID2RANK(exec->ds->syms.proc_ids[select[i]], exec->ds->syms.ranks):-1);

  return select;
}


static void make_sym_hierarchic(ds_exec_t *exec, dsint_t randomize)
{
  dsint_t i, j, n;
  int r, low, high, mid;
  dsint_t *order, *select;
#ifdef DASH_SYMMETRIC_AUX
  sym_aux_t sym_aux;
  dsint_t sym_aux_success;
#endif


  order = sym_order_hierarchic_create(exec, randomize);

  low = 0;
  high = exec->ds->comm_size;

#ifdef DASH_SYMMETRIC_AUX
  sym_aux_setup(exec, &sym_aux);
#endif

  while (low + 1 < high)
  {
    mid = (low + high) / 2;

    select = sym_order_hierarchic_select(exec, randomize, order, low, mid, high, &n);

    for (i = 0; i < n; ++i)
    {
      j = select[i];

      if (j < 0) continue;

      r = PROC_ID2RANK(exec->ds->syms.proc_ids[select[i]], exec->ds->syms.ranks);

      Z_TRACE_IF(DS_TRACE_IF, "sendrecv_replace with %d", r);

      if (exec->ds->syms.counts[j] > 0)
      {
#ifdef DASH_SYMMETRIC_AUX
        sym_aux_success = sym_aux_sendrecv(exec, &sym_aux, j);
        if (!sym_aux_success)
        {
#endif
          exec->sendrecv_replace(exec, r, exec->ds->syms.exec_id, exec->ds->syms.buf_ids[j], exec->ds->syms.displs[j], exec->ds->syms.counts[j]);
#ifdef DASH_SYMMETRIC_AUX
# ifdef DASH_SYMMETRIC_AUX_IMMEDIATELY
          sym_aux_intermediate(exec, &sym_aux);
# endif
        }
#endif
      }
    }

    if (exec->ds->comm_rank < mid)
    {
      high = mid;

    } else
    {
      low = mid;
    }
  }

#ifdef DASH_SYMMETRIC_AUX
  sym_aux_finish(exec, &sym_aux);
#endif

  sym_order_destroy(order);
}


#endif /* OLD_HIERARCHIC */


static void make_sym(ds_exec_t *exec)
{
  dsint_t i;


  Z_TRACE_IF(DS_TRACE_IF, "START");

  Z_TRACE_IF(DS_TRACE_IF, "comm_size: %d, syms: max: %" dsint_fmt ", n: %" dsint_fmt, exec->ds->comm_size, exec->ds->syms.nmax, exec->ds->syms.n);
  for (i = 0; i < exec->ds->syms.n; ++i)
    Z_TRACE_IF(DS_TRACE_IF, "  %" dsint_fmt ": count: %" dspint_fmt ", displ: %" dspint_fmt ", proc_id: %" dspint_fmt, i, exec->ds->syms.counts[i], exec->ds->syms.displs[i], exec->ds->syms.proc_ids[i]);

  sym_prepare_full(exec);

  i = exec->make_sym;

switch_sym:

  switch (exec->make_sym)
  {
    case DASH_EXEC_MAKE_SYM_LINEAR:
      Z_TRACE_IF(DS_TRACE_IF, "make_sym_linear");
      make_sym_linear(exec, 0);
      break;
    case DASH_EXEC_MAKE_SYM_LINEAR_RANDOM:
      Z_TRACE_IF(DS_TRACE_IF, "make_sym_linear_random");
      make_sym_linear(exec, 1);
      break;
    case DASH_EXEC_MAKE_SYM_HIERARCHIC:
      Z_TRACE_IF(DS_TRACE_IF, "make_sym_hierarchic");
      make_sym_hierarchic(exec, 0);
      break;
    case DASH_EXEC_MAKE_SYM_HIERARCHIC_RANDOM:
      Z_TRACE_IF(DS_TRACE_IF, "make_sym_hierarchic_random");
      make_sym_hierarchic(exec, 1);
      break;
    default:
      Z_TRACE_IF(DS_TRACE_IF, "make_sym default");
      i = DASH_EXEC_MAKE_SYM_DEFAULT;
      goto switch_sym;
  }

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


#endif /* DASH_SYMMETRIC */


dsint_t ds_exec_make(ds_exec_t *exec) /* ds_func ds_exec_make */
{
#ifdef DASH_SYMMETRIC
  if (exec->ds->syms.execute)
  {
    make_sym(exec);
    exec->ds->syms.execute = 0;
  }
#endif

  if (exec->ds->sends.execute || exec->ds->recvs.execute || exec->ds->locals.execute)
  {
    exec->make(exec);
    exec->ds->sends.execute = exec->ds->recvs.execute = exec->ds->locals.execute = 0;
  }

  return 0;
}


#undef PROC_ID2RANK

#undef DS_TRACE_IF



#include "dash_exec.h"
#include "dash_exec_sl.h"


#ifndef DS_TRACE_IF
# define DS_TRACE_IF  (z_mpi_rank == -1)
#endif


static dsint_t ds_exec_sl_pre_run(ds_exec_t *exec);
static dsint_t ds_exec_sl_post_run(ds_exec_t *exec);
static dsint_t ds_exec_sl_make(ds_exec_t *exec);

static void ds_exec_sl_move(ds_exec_t *exec, dsint_t exec_id, dsint_t src_buf_id, dsint_t src_displs, dsint_t dst_buf_id, dsint_t dst_displs, dsint_t count);
#ifdef DASH_SYMMETRIC
static void ds_exec_sl_sendrecv_replace(ds_exec_t *exec, int proc, dsint_t exec_id, dsint_t buf_id, dspint_t displ, dspint_t count);
#endif

dsint_t ds_exec_sl_create(ds_exec_t *exec) /* sl_func ds_exec_sl_create */
{
  ds_exec_sl_t *exec_sl;


  ds_exec_create(exec);

  exec->cxt = exec_sl = z_alloc(1, sizeof(ds_exec_sl_t));

  exec->pre_run = ds_exec_sl_pre_run;
  exec->post_run = ds_exec_sl_post_run;
  exec->make = ds_exec_sl_make;
  exec->move = ds_exec_sl_move;

#ifdef DASH_SYMMETRIC
  exec->sendrecv_replace = ds_exec_sl_sendrecv_replace;
#endif

  exec_sl->nmax = 0;
  exec_sl->n = 0;
  exec_sl->reqs = NULL;
  exec_sl->stats = NULL;

  return 0;
}


dsint_t ds_exec_sl_destroy(ds_exec_t *exec) /* sl_func ds_exec_sl_destroy */
{
  z_free(exec->cxt);

  ds_exec_destroy(exec);

  return 0;
}


static dsint_t ds_exec_sl_pre_run(ds_exec_t *exec)
{
  DEFINE_EXEC_SL(exec, exec_sl);


  exec_sl->nmax = exec->ds->sends.nmax + exec->ds->recvs.nmax;
  exec_sl->n = 0;
  exec_sl->reqs = z_alloc(exec_sl->nmax, sizeof(MPI_Request));
  exec_sl->stats = z_alloc(exec_sl->nmax, sizeof(MPI_Status));

  return 0;
}


static dsint_t ds_exec_sl_post_run(ds_exec_t *exec)
{
  DEFINE_EXEC_SL(exec, exec_sl);


  z_free(exec_sl->reqs); exec_sl->reqs = NULL;
  z_free(exec_sl->stats); exec_sl->stats = NULL;

  return 0;
}


static dsint_t ds_exec_sl_make(ds_exec_t *exec)
{
  DEFINE_EXEC_SL(exec, exec_sl);
  dsint_t i, req_id, j;
  dsint_t src_buf_id, dst_buf_id;
  elements_t *src_addr, *dst_addr;

  dsint_t first = 1;
  int rcount;


#define xelem_call \
  req_id = 0; \
  for (i = 0; i < exec->ds->recvs.n; ++i) \
  { \
    dst_buf_id = exec->ds->recvs.buf_ids[i]; \
    dst_addr = exec->addrs[exec->ds->sched->bufs[dst_buf_id].addr_id]; \
    j = exec->ds->recvs.proc_ids[i]; \
    if (exec->ds->recvs.ranks) j = exec->ds->recvs.ranks[j]; \
    MPI_Irecv(xelem_buf_at(dst_addr, exec->ds->recvs.displs[i]), exec->ds->recvs.counts[i], xelem_mpi_datatype, j, DS_EXEC_SL_ISENDRECV_TAG, exec->ds->comm, &exec_sl->reqs[req_id]); \
    ++req_id; \
  } \
  for (i = 0; i < exec->ds->sends.n; ++i) \
  { \
    src_buf_id = exec->ds->sends.buf_ids[i]; \
    src_addr = exec->addrs[exec->ds->sched->bufs[src_buf_id].addr_id]; \
    j = exec->ds->sends.proc_ids[i]; \
    if (exec->ds->sends.ranks) j = exec->ds->sends.ranks[j]; \
    MPI_Isend(xelem_buf_at(src_addr, exec->ds->sends.displs[i]), exec->ds->sends.counts[i], xelem_mpi_datatype, j, DS_EXEC_SL_ISENDRECV_TAG, exec->ds->comm, &exec_sl->reqs[req_id]); \
    ++req_id; \
  } \
  for (i = 0; i < exec->ds->locals.n; ++i) \
  { \
    src_buf_id = exec->ds->locals.src_buf_ids[i]; \
    dst_buf_id = exec->ds->locals.dst_buf_ids[i]; \
    src_addr = exec->addrs[exec->ds->sched->bufs[src_buf_id].addr_id]; \
    dst_addr = exec->addrs[exec->ds->sched->bufs[dst_buf_id].addr_id]; \
    if (src_buf_id == dst_buf_id && z_max(exec->ds->locals.src_displs[i], exec->ds->locals.src_displs[i]) < z_min(exec->ds->locals.src_displs[i] + exec->ds->locals.src_counts[i], exec->ds->locals.src_displs[i] + exec->ds->locals.src_counts[i])) \
    { \
      xelem_nmove_at(src_addr, exec->ds->locals.src_displs[i], dst_addr, exec->ds->locals.dst_displs[i], exec->ds->locals.src_counts[i]); \
    } else \
    { \
      xelem_ncopy_at(src_addr, exec->ds->locals.src_displs[i], dst_addr, exec->ds->locals.dst_displs[i], exec->ds->locals.src_counts[i]); \
    } \
  } \
  MPI_Waitall(req_id, exec_sl->reqs, exec_sl->stats); \
  if (first) \
  { \
    for (i = 0; i < exec->ds->recvs.n; ++i) \
    { \
      MPI_Get_count(&exec_sl->stats[i], xelem_mpi_datatype, &rcount); \
      exec->ds->recvs.counts[i] = rcount; \
    } \
    first = 0; \
  }
#include "sl_xelem_call.h"

  return 0;
}


dsint_t ds_exec_sl_add_address(ds_exec_t *exec, elements_t *s) /* sl_func ds_exec_sl_add_address */
{
  return ds_exec_add_address(exec, s);
}


static void ds_exec_sl_move(ds_exec_t *exec, dsint_t exec_id, dsint_t src_buf_id, dsint_t src_displs, dsint_t dst_buf_id, dsint_t dst_displs, dsint_t count)
{
  elements_t *src_addr, *dst_addr;


  src_addr = exec->addrs[exec->ds->sched->bufs[src_buf_id].addr_id];
  dst_addr = exec->addrs[exec->ds->sched->bufs[dst_buf_id].addr_id];

  if (src_buf_id != dst_buf_id || src_displs + count < dst_displs || dst_displs + count < src_displs)
  {
#define xelem_call \
    xelem_ncopy_at(src_addr, src_displs, dst_addr, dst_displs, count);
#include "sl_xelem_call.h"
  } else
  {
#define xelem_call \
    xelem_nmove_at(src_addr, src_displs, dst_addr, dst_displs, count);
#include "sl_xelem_call.h"
  }
}


#ifdef DASH_SYMMETRIC


static void ds_exec_sl_sendrecv_replace(ds_exec_t *exec, int proc, dsint_t exec_id, dsint_t buf_id, dspint_t displ, dspint_t count)
{
  elements_t *addr;
  MPI_Status status;


#define xelem_call \
  addr = exec->addrs[exec->ds->sched->bufs[buf_id].addr_id]; \
  MPI_Sendrecv_replace(addr, count, xelem_mpi_datatype, proc, DS_EXEC_SL_SENDRECV_REPLACE_TAG, proc, DS_EXEC_SL_SENDRECV_REPLACE_TAG, exec->ds->comm, &status);
#include "sl_xelem_call.h"
}


#endif /* DASH_SYMMETRIC */


#undef DS_TRACE_IF



#include "dash_sched.h"
#include "dash_core.h"


#ifndef DS_TRACE_IF
# define DS_TRACE_IF  (z_mpi_rank == -1)
#endif


dsint_t ds_sched_create(ds_sched_t *sched) /* ds_func ds_sched_create */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  sched->round = 0;

  sched->nbufs = 0;

  sched->max_n = NULL;
  sched->pre_run = NULL;
  sched->post_run = NULL;
  sched->finished = NULL;
  sched->pre = NULL;
  sched->post = NULL;

  sched->cxt = NULL;

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


dsint_t ds_sched_destroy(ds_sched_t *sched) /* ds_func ds_sched_destroy */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


dsint_t ds_sched_add_buffer(ds_sched_t *sched, dsint_t addr_id) /* ds_func ds_sched_add_buffer */
{
  dsint_t buf_id = -1;

  if (sched->nbufs < DASH_MAX_NBUFFERS)
  {
    buf_id = sched->nbufs;
    ++sched->nbufs;

    sched->bufs[buf_id].addr_id = addr_id;
  }

  return buf_id;
}


dsint_t ds_sched_set_send(ds_sched_t *sched, dsint_t buf_id, dsint_t n, dspint_t *scounts, dspint_t *sdispls, dspint_t *sranks, dsint_t exec_id) /* ds_func ds_sched_set_send */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  sched->send.buf_id = buf_id;
  sched->send.exec_id = exec_id;

  sched->send.n = n;
  sched->send.counts = scounts;
  sched->send.displs = sdispls;
  sched->send.ranks = sranks;

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return buf_id;
}


dsint_t ds_sched_set_recv(ds_sched_t *sched, dsint_t buf_id, dsint_t n, dspint_t *rcounts, dspint_t *rdispls, dspint_t *rranks, dsint_t exec_id) /* ds_func ds_sched_set_recv */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  sched->recv.buf_id = buf_id;
  sched->recv.exec_id = exec_id;

  sched->recv.n = n;
  sched->recv.counts = rcounts;
  sched->recv.displs = rdispls;
  sched->recv.ranks = rranks;

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return buf_id;
}


dsint_t ds_sched_set_aux(ds_sched_t *sched, dsint_t buf_id, dsint_t buf_size) /* ds_func ds_sched_set_aux */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  sched->aux.buf_id = buf_id;
  sched->aux.buf_size = buf_size;

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return buf_id;
}


#undef DS_TRACE_IF



#include "dash_sched_a2av.h"


#ifndef DS_TRACE_IF
# define DS_TRACE_IF  (z_mpi_rank == -1)
#endif

#define DS_PRINT_BLOCKS

dsint_t ds_sched_a2av_aux_blocks = -8; /* ds_var ds_sched_a2av_aux_blocks */


static dsint_t ds_sched_a2av_max_n(ds_sched_t *sched, dsint_t *max_nsends, dsint_t *max_n, dsint_t *max_nlocals, dsint_t *max_nsyms);
static dsint_t ds_sched_a2av_pre_run(ds_sched_t *sched);
static dsint_t ds_sched_a2av_post_run(ds_sched_t *sched);
static dsint_t ds_sched_a2av_finished(ds_sched_t *sched);
static dsint_t ds_sched_a2av_pre(ds_sched_t *sched);
static dsint_t ds_sched_a2av_post(ds_sched_t *sched);

/* block management */
static void ds_init_blocks(ds_sched_t *sched);
static void ds_free_blocks(ds_sched_t *sched);
static void ds_print_blocks(ds_sched_t *sched);
static void ds_insert_blocks(ds_sched_t *sched);
#ifdef DASH_SCHED_A2AV_OVERLAP
static void ds_resort_sblocks(ds_sched_a2av_t *sched_a2av, dsint_t b);
static void ds_recreate_matchings(ds_sched_a2av_t *sched_a2av, dsint_t b);
#endif
static void ds_compactify_blocks(ds_sched_t *sched, dsint_t b);

static block_t *ds_add_block(ds_sched_a2av_t *sched_a2av, dspint_t displ, dspint_t count, dspint_t proc_id, block_t *prev, dsint_t buf_id, dsint_t send);
static void ds_del_block(ds_sched_a2av_t *sched_a2av, block_t *blk, dsint_t buf_id, dsint_t send);
static block_t *ds_split_block(ds_sched_a2av_t *sched_a2av, block_t *blk, dsint_t p);
static void ds_check_block(ds_sched_a2av_t *sched_a2av, block_t *blk, dsint_t buf_id, dsint_t send);

/* misc */
static void ds_sort_displs(dsint_t ndispls, dspint_t *displs, dspint_t *counts, dsint_t *order);
#ifdef DASH_SCHED_A2AV_OVERLAP
static void ds_sort_blocks_begin(dsint_t n, block_t **blks);
#endif

/* block memory pooling */
static void ds_block_pool_add_buffer(block_pool_t *bp, void *buf, dsint_t bufsize);
static void ds_block_pool_init(block_pool_t *bp, dsint_t nblocks);
static void ds_block_pool_free(block_pool_t *bp);
static block_t *ds_block_pool_get(block_pool_t *bp);
static void ds_block_pool_put(block_pool_t *bp, block_t *blk);

/* send/receive request handling */
static void ds_req_init(ds_sched_t *sched);
static void ds_req_free(ds_sched_t *sched);
static void ds_req_reset_recvs(ds_sched_t *sched);
static void ds_req_exchange(ds_sched_t *sched);
static void ds_req_set_recv(ds_sched_t *sched, dsint_t proc_id, dsint_t req_id, dsint_t n);
static dsint_t ds_req_get_recv(ds_sched_t *sched, dsint_t proc_id, dsint_t req_id);
static dsint_t ds_req_get_send(ds_sched_t *sched, dsint_t proc_id, dsint_t req_id);


dsint_t ds_sched_a2av_create(ds_sched_t *sched) /* ds_func ds_sched_a2av_create */
{
  ds_sched_a2av_t *sched_a2av;


  Z_TRACE_IF(DS_TRACE_IF, "START");

  ds_sched_create(sched);

  sched->max_n = ds_sched_a2av_max_n;
  sched->pre_run = ds_sched_a2av_pre_run;
  sched->post_run = ds_sched_a2av_post_run;
  sched->finished = ds_sched_a2av_finished;
  sched->pre = ds_sched_a2av_pre;
  sched->post = ds_sched_a2av_post;

  sched_a2av = z_alloc(1, sizeof(ds_sched_a2av_t));
  sched_a2av->aux = NULL;
  sched_a2av->skip_sym = 0;

  sched->cxt = sched_a2av;

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


dsint_t ds_sched_a2av_destroy(ds_sched_t *sched) /* ds_func ds_sched_a2av_destroy */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  z_free(sched->cxt);

  ds_sched_destroy(sched);

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


void ds_sched_a2av_set_aux(ds_sched_t *sched, ds_sched_a2av_aux_t *aux) /* ds_func ds_sched_a2av_set_aux */
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);


  sched_a2av->aux = aux;
}


void ds_sched_a2av_skip_sym(ds_sched_t *sched, dsint_t skip_sym) /* ds_func ds_sched_a2av_skip_sym */
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);


  sched_a2av->skip_sym = skip_sym;
}


static dsint_t ds_sched_a2av_max_n(ds_sched_t *sched, dsint_t *max_nsends, dsint_t *max_nrecvs, dsint_t *max_nlocals, dsint_t *max_nsyms)
{
  *max_nsends = NREQ_PER_PROC * sched->send.n;
  *max_nrecvs = NREQ_PER_PROC * sched->recv.n;
  *max_nlocals = NREQ_PER_PROC * 1;

  *max_nsyms = z_min(sched->send.n, sched->recv.n);

  return 0;
}


static dsint_t ds_sched_a2av_pre_run(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);
  dsint_t i;


  Z_TRACE_IF(DS_TRACE_IF, "START");

  if (sched->send.ranks)
  {
    sched_a2av->local_send_proc_id = -1;
    for (i = 0; i < sched->send.n; ++i)
    if (sched->send.ranks[i] == sched->ds->comm_rank) { sched_a2av->local_send_proc_id = i; break; }

  } else sched_a2av->local_send_proc_id = sched->ds->comm_rank;

  if (sched->recv.ranks)
  {
    sched_a2av->local_recv_proc_id = -1;
    for (i = 0; i < sched->recv.n; ++i)
    if (sched->recv.ranks[i] == sched->ds->comm_rank) { sched_a2av->local_recv_proc_id = i; break; }

  } else sched_a2av->local_recv_proc_id = sched->ds->comm_rank;

  ds_block_pool_init(&sched_a2av->bp, 3 * sched->send.n + sched->recv.n);

  ds_init_blocks(sched);

  ds_req_init(sched);

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


static dsint_t ds_sched_a2av_post_run(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  Z_TRACE_IF(DS_TRACE_IF, "START");

  ds_free_blocks(sched);

  ds_block_pool_free(&sched_a2av->bp);

  ds_req_free(sched);

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


static dsint_t ds_sched_a2av_finished(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);
  dsint_t local, global;
/*  const dsint_t allreduce_interval = 10;*/

  Z_TRACE_IF(DS_TRACE_IF, "START: round = %" dsint_fmt ", stotal = %" dsint_fmt ", rtotal = %" dsint_fmt, sched->round, sched_a2av->stotal, sched_a2av->rtotal);

#ifdef DS_PRINT_BLOCKS
  ds_print_blocks(sched);
#endif

/*  if ((sched->round % allreduce_interval) + 1 == allreduce_interval)
  {*/
    local = (sched_a2av->stotal != 0 || sched_a2av->rtotal != 0);

    MPI_Allreduce(&local, &global, 1, MPI_DSINT, MPI_SUM, sched->ds->comm);

/*  } else global = 1;*/

#if 0
  if (sched->round > 4) global = 0;
#endif

  Z_TRACE_IF(DS_TRACE_IF, "END: %s", (global == 0)?"finished":"continue");

  return (global == 0);
}


#ifdef DASH_SYMMETRIC

static void ds_sched_a2av_pre_symmetric(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  dsint_t b, i, j, naux;

  dsint_t match_group_end, sblk_begin, sblk_end, sblk_begin_extra;
  block_t *sblk, *rblk, **sblks;
  dsint_t *sblk_skips;
  dsint_t my_sym[3], ot_sym[3], delta_s, delta_r, length, delta;


  Z_TRACE_IF(DS_TRACE_IF, "START: round = %" dsint_fmt, sched->round);

  sblks = z_alloc(sched_a2av->bufs[0].nrblocks, sizeof(block_t *));
  sblk_skips = z_alloc(sched_a2av->bufs[0].nrblocks, sizeof(dsint_t));

  /* reset request handling */
  ds_req_reset_recvs(sched);

  /* check receive blocks */
  for (b = 0; b < sched->nbufs; ++b)
  {
    rblk = sched_a2av->bufs[b].rblock_first;

    sblk = NULL;

    if (rblk != NULL) sblk_begin_extra = rblk->begin;

    while (rblk != NULL)
    {
      j = rblk->proc_id;

      sblks[j] = NULL;
      sblk_skips[j] = 0;

      Z_TRACE_IF(DS_TRACE_IF, "proc: %" dsint_fmt ", check rblock: %p: " PRINT_BLOCK_STR, j, rblk, PRINT_BLOCK_PARAMS(*rblk));

      if (rblk->match == NULL) sblk = NULL;

      if (sblk != NULL && sblk->end <= rblk->begin) sblk = NULL;

      if (sblk == NULL) sblk = rblk->match;

      if (sblk != NULL)
      {
        /* account for the items that are already in aux */
        naux = ds_aux_get_count(sched_a2av->aux, j);

        match_group_end = sblk->begin - 1;

        while (1)
        {
          Z_TRACE_IF(DS_TRACE_IF, "match_group_end: %" dsint_fmt ", sblk: %p", match_group_end, sblk);

          /* current send-block and next send-block NOT connected or overlapping? */
          if (match_group_end < sblk->begin)
          {
            Z_TRACE_IF(DS_TRACE_IF, "start new match-group with %p", sblk);

            /* start new match-group */
            match_group_end = sblk->end;

            sblk_begin = z_max(rblk->begin, sblk_begin_extra);
            sblk_end = rblk->end;

          } else
          {
            Z_TRACE_IF(DS_TRACE_IF, "extent match-group with %p", sblk);

            /* extent current match-group */
            match_group_end = z_max(match_group_end, sblk->end);
          }

          if (sblk->proc_id == j)
          {
            sblks[j] = sblk;

            sblk_begin = z_max(sblk_begin, sblk->begin + naux);
            sblk_end = z_min(sblk_end, sblk->end);

            sblk_begin_extra = z_max(sblk_begin_extra, sblk->end);

          } else
          {
            if (sblks[j] == NULL) sblk_begin = z_max(sblk_begin, sblk->end);
            else sblk_end = z_min(sblk_end, sblk->begin);
          }

          Z_TRACE_IF(DS_TRACE_IF, "updated sblk: [%" dsint_fmt ",%" dsint_fmt "]", sblk_begin, sblk_end);

          if (sblk->proc_id != j && rblk->end <= sblk->end) break;

          Z_TRACE_IF(DS_TRACE_IF, "after break 1");

          if (sblk->next == NULL || sblk->next->begin >= rblk->end) break;

          Z_TRACE_IF(DS_TRACE_IF, "after break 2");

          sblk = sblk->next;
        }

        if (sblks[j] != NULL)
        {
          Z_TRACE_IF(DS_TRACE_IF, "sblk_begin: %" dsint_fmt ", sblk_end: %" dsint_fmt, sblk_begin, sblk_end);

          sblk_skips[j] = sblk_begin - sblks[j]->begin;

          Z_TRACE_IF(DS_TRACE_IF, "%" dsint_fmt " - %" dspint_fmt " - %" dsint_fmt "", sblk_begin, sched->send.displs[j], sblks[j]->moved);
          Z_TRACE_IF(DS_TRACE_IF, "%" dsint_fmt " - %" dspint_fmt "", sblk_begin, sched->recv.displs[j]);

          my_sym[0] = sblk_begin - sched->send.displs[j] - sblks[j]->moved;
          my_sym[1] = sblk_begin - sched->recv.displs[j];
          my_sym[2] = sblk_end - sblk_begin;

          Z_TRACE_IF(DS_TRACE_IF, "proc %" dsint_fmt ": sym request: %" dsint_fmt ", %" dsint_fmt ", %" dsint_fmt, j, my_sym[0], my_sym[1], my_sym[2]);

          ds_req_set_recv(sched, j, 0, my_sym[0]);
          ds_req_set_recv(sched, j, 1, my_sym[1]);
          ds_req_set_recv(sched, j, 2, my_sym[2]);
        }
      }

      rblk = rblk->next;
    }
  }

  /* exchange requests */
  ds_req_exchange(sched);

  /* prepare symmetric exchanges */
  sched->ds->syms.n = 0;
/*  sched->ds->syms.execute = 0;*/
  sched->ds->syms.exec_id = sched->send.exec_id;

  /* check receive requests */
  for (b = 0; b < sched->nbufs; ++b)
  for (i = 0; i < sched_a2av->bufs[b].nrblocks; ++i)
  {
    j = (i + z_max(sched_a2av->local_recv_proc_id, 0) + sched_a2av->bufs[b].nrblocks) % sched_a2av->bufs[b].nrblocks;

    my_sym[0] = ds_req_get_recv(sched, j, 0);
    my_sym[1] = ds_req_get_recv(sched, j, 1);
    my_sym[2] = ds_req_get_recv(sched, j, 2);

    ot_sym[0] = ds_req_get_send(sched, j, 0);
    ot_sym[1] = ds_req_get_send(sched, j, 1);
    ot_sym[2] = ds_req_get_send(sched, j, 2);

    Z_TRACE_IF(DS_TRACE_IF, "proc %" dsint_fmt ": sym: %" dsint_fmt ", %" dsint_fmt ", %" dsint_fmt " vs. %" dsint_fmt ", %" dsint_fmt ", %" dsint_fmt,
      j, my_sym[0], my_sym[1], my_sym[2], ot_sym[0], ot_sym[1], ot_sym[2]);

    if (my_sym[2] > 0 && ot_sym[2] > 0)
    {
      delta_s = ot_sym[1] - my_sym[0];
      delta_r = ot_sym[0] - my_sym[1];

      Z_TRACE_IF(DS_TRACE_IF, "proc %" dsint_fmt ": delta_s = %" dsint_fmt ", delta_r = %" dsint_fmt, j, delta_s, delta_r);

      if (delta_s == delta_r)
      {
        delta = 0;
        if (delta_s > 0) { my_sym[2] -= delta_s; delta = delta_s; }
        if (delta_s < 0) { ot_sym[2] += delta_s; }
        length = z_min(my_sym[2], ot_sym[2]);

        Z_TRACE_IF(DS_TRACE_IF, "proc %" dsint_fmt ": sendrecv_replace: " dsint_fmt ": offset = %" dsint_fmt ", size = %" dsint_fmt, j, sblks[j]->begin + sblk_skips[j] + delta, length);

        if (length > 0)
        {
          sched->ds->syms.buf_ids[sched->ds->syms.n] = b;
          sched->ds->syms.counts[sched->ds->syms.n] = length;
          sched->ds->syms.displs[sched->ds->syms.n] = sblks[j]->begin + sblk_skips[j] + delta;
          sched->ds->syms.proc_ids[sched->ds->syms.n] = j;

          ++sched->ds->syms.n;
        }
      }
    }
  }

  sched->ds->syms.execute = (sched->ds->syms.n > 0);
  sched->ds->syms.ranks = sched->recv.ranks;

  z_free(sblks);
  z_free(sblk_skips);

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


static void ds_accept_sym(ds_sched_a2av_t *sched_a2av, dsint_t buf_id, dsint_t proc_id, dsint_t displ, dsint_t count)
{
  block_t *sblk, *rblk;


  rblk = sched_a2av->bufs[buf_id].rblocks[proc_id];

  /* select symblock directly and go to its second if necessary */
  sblk = sched_a2av->bufs[buf_id].sblocks[proc_id];
  if (sblk->end <= rblk->begin || rblk->end <= sblk->begin) sblk = sblk->second;

  /* accept sym in sblock */
  sblk->sym_count = count;
  sblk->sym_displ = displ;
  sblk->sym_offset = displ - sblk->begin;

  Z_TRACE_IF(DS_TRACE_IF, "send-block: %p: " PRINT_BLOCK_STR ", sym_count: %" dsint_fmt ", sym_displ: %" dsint_fmt ", sym_offset: %" dsint_fmt,
    sblk, PRINT_BLOCK_PARAMS(*sblk), sblk->sym_count, sblk->sym_displ, sblk->sym_offset);

  ds_check_block(sched_a2av, sblk, buf_id, 1);

  /* accept sym in rblock */
  rblk->sym_count = count;
  rblk->sym_displ = displ;
  rblk->sym_offset = displ - rblk->begin;

  Z_TRACE_IF(DS_TRACE_IF, "receive-block: %p " PRINT_BLOCK_STR ", sym_count: %" dsint_fmt ", sym_displ: %" dsint_fmt ", sym_offset: %" dsint_fmt,
    rblk, PRINT_BLOCK_PARAMS(*rblk), rblk->sym_count, rblk->sym_displ, rblk->sym_offset);

  ds_check_block(sched_a2av, rblk, buf_id, 0);
}


static void ds_sched_a2av_post_symmetric(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  dsint_t i;


  Z_TRACE_IF(DS_TRACE_IF, "START: round = %" dsint_fmt, sched->round);

  /* accept symmetrics */
  for (i = 0; i < sched->ds->syms.n; ++i)
  if (sched->ds->syms.counts[i])
  {
    Z_TRACE_IF(DS_TRACE_IF, "accept symmetric: %" dspint_fmt " @ %" dspint_fmt " to %" dspint_fmt, sched->ds->syms.counts[i], sched->ds->syms.displs[i], sched->ds->syms.proc_ids[i]);

    ds_accept_sym(sched_a2av, sched->ds->syms.buf_ids[i], sched->ds->syms.proc_ids[i], sched->ds->syms.displs[i], sched->ds->syms.counts[i]);
  }

  Z_TRACE_IF(DS_TRACE_IF, "END");
}

#endif


static dsint_t ds_get_connected(dsint_t max_count, block_t *blk, dsint_t blk_skip, block_t **next, dsint_t *next_skip)
{
  dsint_t count, l, end;
#ifdef DASH_SYMMETRIC
  dsint_t sym_cross;
#endif

  Z_TRACE_IF(DS_TRACE_IF, "get %" dsint_fmt " of block %p with skip %" dsint_fmt, max_count, blk, blk_skip);

  count = 0;
  l = 0;
  end = (blk)?blk->begin:-1;

  while (max_count > 0 && blk != NULL && blk->begin == end)
  {
    l = z_min(max_count, blk->end - blk->begin - blk_skip);

#ifdef DASH_SYMMETRIC
    if (blk->sym_count != 0)
    {
      /* if skipped begin is left of sym */
      if (blk->begin + blk_skip < blk->sym_displ) l = z_min(l, blk->sym_displ - blk->begin - blk_skip);  /* ... then limit to start of sym */
      else                                                                                               /* ... else skipped begin is right of sym */
      {
        /* skip crosses sym? */
        sym_cross = (blk->begin < blk->sym_displ && blk->sym_displ < blk->begin + blk_skip)?blk->sym_count:0;

        if (blk_skip - sym_cross < blk->sym_offset) l = z_min(l, blk->sym_offset - blk_skip + sym_cross);
      }
    }
#endif

    count += l;
    max_count -= l;

    blk_skip += l;

#ifdef DASH_SYMMETRIC
    /* if any sym and current skip hits this sym */
    if (blk->sym_count != 0 && blk->begin + blk_skip == blk->sym_displ)
    {
      /* skip the sym and break (cannot go on since blk is interrupted by the sym!) */
      blk_skip += z_abs(blk->sym_count);

      break;
    }
#endif

    if (blk->begin + blk_skip >= blk->end)
    {
      end = blk->end;

#ifdef DASH_SCHED_A2AV_OVERLAP
      /* skip to second part */
      blk = blk->second;
#else
      /* skip to next part */
      blk = blk->next;
#endif
      blk_skip = 0;
    }
  }

  *next = blk;
  *next_skip = blk_skip;

  Z_ASSERT(!(max_count > 0 && blk == NULL));

  return count;
}


static dsint_t ds_add_aux_recv_request(void *schedptr, dsint_t src, dsint_t n)
{
  ds_sched_t *sched = schedptr;

  if (n > 0)
  {
    ds_req_set_recv(sched, src, 1, -n);
    Z_TRACE_IF(DS_TRACE_IF, "recv-req-aux: #%d: recv %" dsint_fmt " from %" dsint_fmt, 1, n, src);
  }

  return (n > 0);
}


static void ds_sched_a2av_pre_individual(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  dsint_t b, i, j, n, isnotaux, r;
  block_t *sblk, *rblk;
  dsint_t sblk_skip;


  Z_TRACE_IF(DS_TRACE_IF, "START: round = %" dsint_fmt, sched->round);

  /* reset request handling */
  ds_req_reset_recvs(sched);

  /* pre acquire aux */
  ds_aux_pre_acquire(sched_a2av->aux);

  /* make recv-requests */
  for (b = 0; b < sched->nbufs; ++b)
  for (i = 0; i < sched_a2av->bufs[b].nrblocks; ++i)
  {
    j = (i + z_max(sched_a2av->local_recv_proc_id, 0) + sched_a2av->bufs[b].nrblocks) % sched_a2av->bufs[b].nrblocks;

    rblk = sched_a2av->bufs[b].rblocks[j];

    if (rblk == NULL) continue;

    Z_TRACE_IF(DS_TRACE_IF, "proc: %" dsint_fmt ", check rblock: %p: " PRINT_BLOCK_STR, j, rblk, PRINT_BLOCK_PARAMS(*rblk));

    /* have a match? */
    if (rblk->match)
    {
      Z_TRACE_IF(DS_TRACE_IF, "have match: %p [%" dsint_fmt ",%" dsint_fmt "]", rblk->match, rblk->match->begin, rblk->match->end);

      if (sched_a2av->local_recv_proc_id == j)
      {
        /* determine local receive*/

        /* if match is local send-block AND no local items stored in aux then optimized local move is possible */
        if (rblk->match == sched_a2av->bufs[b].sblocks[j] && ds_aux_get_count(sched_a2av->aux, j) <= 0)
        {
          if (rblk->match->next && rblk->match->next->begin < rblk->end)
          {
            n = z_max(0, rblk->match->next->begin - rblk->begin);

            Z_TRACE_IF(DS_TRACE_IF, "limited optimized local receive, n = %" dsint_fmt, n);

          } else
          {
            n = z_max(0, rblk->end - rblk->begin);

            Z_TRACE_IF(DS_TRACE_IF, "full optimized local receive, n = %" dsint_fmt, n);
          }

        } else
        {
          n = z_max(0, rblk->match->begin - rblk->begin);

          Z_TRACE_IF(DS_TRACE_IF, "normal local receive, n = %" dsint_fmt, n);
        }

      } else
      {
        /* normal receive up to the begin of the matching send block (difference is negative, if match is in front of the recv-block) */
        n = z_max(0, rblk->match->begin - rblk->begin);

        Z_TRACE_IF(DS_TRACE_IF, "normal receive, n = %" dsint_fmt, n);
      }

    } else
    {
      /* have no match, can receive all */
      n = rblk->end - rblk->begin;

      Z_TRACE_IF(DS_TRACE_IF, "no match, n = %" dsint_fmt, n);
    }

#ifdef DASH_SYMMETRIC
    /* if receive block has an active sym then limit number of items to receive to the sym offset */
    if (rblk->sym_count != 0)
    {
      n = z_min(n, rblk->sym_offset);

      Z_TRACE_IF(DS_TRACE_IF, "have sym, n = %" dsint_fmt, n);
    }
#endif

    ds_req_set_recv(sched, j, 0, n);
    Z_TRACE_IF(DS_TRACE_IF, "recv-req: %" dsint_fmt " #%d: recv %" dsint_fmt " from %" dsint_fmt, i, 0, n, j);

    /* what's left if the request is fullfilled? */
    n = (rblk->end - rblk->begin) - n;

#ifdef DASH_SYMMETRIC
    if (rblk->sym_count != 0) n -= rblk->sym_count;
#endif

    /* go to next if no items are left */
    if (n <= 0) continue;

    /* account for the items that are already in aux */
    n -= ds_aux_get_count(sched_a2av->aux, j);

    /* acquire aux */
    n = ds_aux_acquire(sched_a2av->aux, n, j);

    if (n > 0) ds_add_aux_recv_request(sched, j, n);
  }

  /* post acquire aux */
  sched_a2av->aux->add_aux_recv_request = ds_add_aux_recv_request;
  ds_aux_post_acquire(sched_a2av->aux);

  /* exchange requests */
  ds_req_exchange(sched);

  /* prepare non-local receives */
  sched->ds->recvs.n = 0;
  sched->ds->recvs.execute = 1;
  sched->ds->recvs.exec_id = sched->recv.exec_id;
  sched->ds->recvs.ranks = sched->recv.ranks;

  for (b = 0; b < sched->nbufs; ++b)
  for (i = 0; i < sched_a2av->bufs[b].nrblocks; ++i)
  {
    j = (i + z_max(sched_a2av->local_recv_proc_id, 0) + sched_a2av->bufs[b].nrblocks) % sched_a2av->bufs[b].nrblocks;

    rblk = sched_a2av->bufs[b].rblocks[j];

    if (rblk == NULL) continue;

    Z_TRACE_IF(DS_TRACE_IF, "proc: %" dsint_fmt ", rblock: %p: " PRINT_BLOCK_STR, j, rblk, PRINT_BLOCK_PARAMS(*rblk));

    for (r = 0; r < NREQ_PER_PROC; ++r)
    {
      n = ds_req_get_recv(sched, j, r);

      Z_TRACE_IF(DS_TRACE_IF, "receive request: %" dsint_fmt ", %" dsint_fmt ", n: %" dsint_fmt, j, r, n);

      /* skip if nothing or local */
      if (n == 0 || sched_a2av->local_recv_proc_id == j) continue;

      Z_ASSERT(sched->recv.buf_id == b);

      isnotaux = (n > 0);
      if (!isnotaux) n *= -1;

      if (isnotaux)
      {
        /* prepare normal receive */
        sched->ds->recvs.buf_ids[sched->ds->recvs.n] = b;
        sched->ds->recvs.counts[sched->ds->recvs.n] = n;
        sched->ds->recvs.displs[sched->ds->recvs.n] = rblk->begin;
        sched->ds->recvs.proc_ids[sched->ds->recvs.n] = j;

        Z_TRACE_IF(DS_TRACE_IF, "prepare receive: #%" dsint_fmt ": %" dspint_fmt " @ %" dspint_fmt " from %" dspint_fmt, \
          sched->ds->recvs.n, sched->ds->recvs.counts[sched->ds->recvs.n], sched->ds->recvs.displs[sched->ds->recvs.n], sched->ds->recvs.proc_ids[sched->ds->recvs.n]);

      } else
      {
        /* prepare aux receive */
        sched->ds->recvs.buf_ids[sched->ds->recvs.n] = sched->aux.buf_id;
        sched->ds->recvs.counts[sched->ds->recvs.n] = n;
        sched->ds->recvs.displs[sched->ds->recvs.n] = ds_aux_get_displ(sched_a2av->aux, j);
        sched->ds->recvs.proc_ids[sched->ds->recvs.n] = j;

        Z_TRACE_IF(DS_TRACE_IF, "prepare aux receive: #%" dsint_fmt ": %" dspint_fmt " @ %" dspint_fmt " from %" dspint_fmt, \
          sched->ds->recvs.n, sched->ds->recvs.counts[sched->ds->recvs.n], sched->ds->recvs.displs[sched->ds->recvs.n], sched->ds->recvs.proc_ids[sched->ds->recvs.n]);
      }

      ++sched->ds->recvs.n;
    }
  }

  /* make exec from send-requests */
  sched->ds->sends.n = 0;
  sched->ds->sends.execute = 1;
  sched->ds->sends.exec_id = sched->send.exec_id;
  sched->ds->sends.ranks = sched->send.ranks;

  sched->ds->locals.n = 0;
  sched->ds->locals.execute = 1;

  for (b = 0; b < sched->nbufs; ++b)
  for (i = 0; i < sched_a2av->bufs[b].nsblocks; ++i)
  {
    j = (i + z_max(sched_a2av->local_send_proc_id, 0) + sched_a2av->bufs[b].nsblocks) % sched_a2av->bufs[b].nsblocks;

    sblk = sched_a2av->bufs[b].sblocks[j];
    sblk_skip = 0;

    if (sblk == NULL) continue;

    Z_TRACE_IF(DS_TRACE_IF, "proc: %" dsint_fmt ", sblock (with skip %" dsint_fmt "): %p: " PRINT_BLOCK_STR, j, sblk_skip, sblk, PRINT_BLOCK_PARAMS(*sblk));

    for (r = 0; r < NREQ_PER_PROC; ++r)
    {
      n = ds_req_get_send(sched, j, r);

      Z_TRACE_IF(DS_TRACE_IF, "send request: %" dsint_fmt ", %" dsint_fmt ", n: %" dsint_fmt, j, r, n);

      /* skip if nothing */
      if (n == 0) continue;

      Z_ASSERT(sched->send.buf_id == b);

      isnotaux = (n > 0);
      if (!isnotaux) n *= -1;

      Z_TRACE_IF(DS_TRACE_IF, "sending: %" dsint_fmt " #%" dsint_fmt": send %" dsint_fmt " to %" dsint_fmt "%s", i, r, n, j, (isnotaux)?"":" (aux)");

      if (sched_a2av->local_send_proc_id != j)
      {
        /* prepare send */
        sched->ds->sends.buf_ids[sched->ds->sends.n] = b;
        sched->ds->sends.displs[sched->ds->sends.n] = sblk->begin + sblk_skip;
        sched->ds->sends.counts[sched->ds->sends.n] = ds_get_connected(n, sblk, sblk_skip, &sblk, &sblk_skip);
        sched->ds->sends.proc_ids[sched->ds->sends.n] = j;

        Z_TRACE_IF(DS_TRACE_IF, "prepare send: #%" dsint_fmt ": %" dspint_fmt " @ %" dspint_fmt " to %" dspint_fmt, \
          sched->ds->sends.n, sched->ds->sends.counts[sched->ds->sends.n], sched->ds->sends.displs[sched->ds->sends.n], sched->ds->sends.proc_ids[sched->ds->sends.n]);

        ++sched->ds->sends.n;

      } else
      {
        /* prepare local */
        sched->ds->locals.src_buf_ids[sched->ds->locals.n] = b;
        sched->ds->locals.src_exec_ids[sched->ds->locals.n] = sched->send.exec_id;

        sched->ds->locals.src_displs[sched->ds->locals.n] = sblk->begin + sblk_skip;
        sched->ds->locals.src_counts[sched->ds->locals.n] = ds_get_connected(n, sblk, sblk_skip, &sblk, &sblk_skip);

        sched->ds->locals.dst_exec_ids[sched->ds->locals.n] = sched->recv.exec_id;
        sched->ds->locals.dst_counts[sched->ds->locals.n] = sched->ds->locals.src_counts[sched->ds->locals.n];

        if (isnotaux)
        {
          sched->ds->locals.dst_buf_ids[sched->ds->locals.n] = sched->recv.buf_id;
          sched->ds->locals.dst_displs[sched->ds->locals.n] = sched_a2av->bufs[sched->recv.buf_id].rblocks[j]->begin;

        } else
        {
          sched->ds->locals.dst_buf_ids[sched->ds->locals.n] = sched->aux.buf_id;
          sched->ds->locals.dst_displs[sched->ds->locals.n] = ds_aux_get_displ(sched_a2av->aux, j);
        }

        Z_TRACE_IF(DS_TRACE_IF, "prepare local: #%" dsint_fmt ": %" dspint_fmt " @ %" dspint_fmt " to%s %" dspint_fmt " @ %" dspint_fmt, \
          sched->ds->locals.n, \
          sched->ds->locals.src_counts[sched->ds->locals.n], sched->ds->locals.src_displs[sched->ds->locals.n],
          (isnotaux?"":" aux"),
          sched->ds->locals.dst_counts[sched->ds->locals.n], sched->ds->locals.dst_displs[sched->ds->locals.n]);

        ++sched->ds->locals.n;
      }
    }
  }

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


static void ds_accept_send(ds_sched_a2av_t *sched_a2av, dsint_t buf_id, dsint_t dst_proc_id, dsint_t count)
{
  dsint_t l;
  block_t *blk;


  sched_a2av->stotal -= count;

  while (count > 0)
  {
    blk = sched_a2av->bufs[buf_id].sblocks[dst_proc_id];

    l = z_min(count, blk->end - blk->begin);

#ifdef DASH_SYMMETRIC
    /* if active sym */
    if (blk->sym_count > 0)
    {
      /* if begin is left of active sym */
      if (blk->begin < blk->sym_displ) l = z_min(l, blk->sym_displ - blk->begin);  /* ... then limit to start of active sym */
      else l = z_min(l, blk->sym_offset);                                          /* ... else limit to sym_offset items */
    }

    /* if passive sym then limit to sym_offset items */
    if (blk->sym_count < 0) l = z_min(l, blk->sym_offset);
#endif

    blk->begin += l;

#ifdef DASH_SYMMETRIC
    if (blk->sym_count != 0) blk->sym_offset -= l;
#endif

    Z_ASSERT(blk->begin <= blk->end);

    ds_check_block(sched_a2av, blk, buf_id, 1);

    count -= l;
  }
}


static void ds_accept_recv(ds_sched_a2av_t *sched_a2av, dsint_t buf_id, dsint_t src_proc_id, dsint_t count)
{
  block_t *blk;


  sched_a2av->rtotal -= count;

  blk = sched_a2av->bufs[buf_id].rblocks[src_proc_id];

  blk->begin += count;

#ifdef DASH_SYMMETRIC
  if (blk->sym_count != 0) blk->sym_offset -= count;
#endif

  Z_ASSERT(blk->begin <= blk->end);

  ds_check_block(sched_a2av, blk, buf_id, 0);
}


static dsint_t ds_vacate_aux(void *schedptr, dspint_t proc_id, dsint_t count, dsint_t displ)
{
  ds_sched_t *sched = schedptr;
  DEFINE_SCHED_A2AV(sched, sched_a2av);
  dsint_t n;
  block_t *rblk;


  Z_TRACE_IF(DS_TRACE_IF, "vacate %" dsint_fmt " @ %" dsint_fmt " to %" dspint_fmt, count, displ, proc_id);

  rblk = sched_a2av->bufs[sched->recv.buf_id].rblocks[proc_id];

  if (rblk->match)
  {
    Z_TRACE_IF(DS_TRACE_IF, "rblk->match->begin = %" dsint_fmt ", rblk->begin = %" dsint_fmt, rblk->match->begin, rblk->begin);

    /* normal receive up to the begin of the matching send block (difference is negative, if match is in front of the recv-block) */
    n = z_max(0, rblk->match->begin - rblk->begin);

    count = z_min(n, count);

    Z_TRACE_IF(DS_TRACE_IF, "have match, limited to count = %" dsint_fmt, count);
  }

#ifdef DASH_SYMMETRIC
  /* if receive block has an active sym then limit number of items to receive to the sym offset */
  if (rblk->sym_count != 0)
  {
    count = z_min(count, rblk->sym_offset);

    Z_TRACE_IF(DS_TRACE_IF, "have sym, limited to count = %" dsint_fmt, count);
  }
#endif

  /* move content */
  sched->ds->exec->move(sched->ds->exec, sched->recv.exec_id, sched->aux.buf_id, displ, sched->recv.buf_id, sched_a2av->bufs[sched->recv.buf_id].rblocks[proc_id]->begin, count);

  ds_accept_recv(sched_a2av, sched->recv.buf_id, proc_id, count);

  return count;
}


static void ds_sched_a2av_post_individual(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  dsint_t i, b;


  Z_TRACE_IF(DS_TRACE_IF, "START: round = %" dsint_fmt, sched->round);

  /* accept locals */
  for (i = 0; i < sched->ds->locals.n; ++i)
  {
    Z_TRACE_IF(DS_TRACE_IF, "accept local send of %" dspint_fmt " to me @ *", sched->ds->locals.src_counts[i]);

    ds_accept_send(sched_a2av, sched->ds->locals.src_buf_ids[i], sched_a2av->local_send_proc_id, sched->ds->locals.src_counts[i]);

    if (sched->ds->locals.dst_buf_ids[i] == sched->recv.buf_id)
    {
      Z_TRACE_IF(DS_TRACE_IF, "accept local receive of %" dspint_fmt " from me @ *", sched->ds->locals.dst_counts[i]);

      ds_accept_recv(sched_a2av, sched->ds->locals.dst_buf_ids[i], sched_a2av->local_recv_proc_id, sched->ds->locals.dst_counts[i]);

    } else
    {
      Z_TRACE_IF(DS_TRACE_IF, "accept local aux receive of %" dspint_fmt " from me @ *", sched->ds->locals.dst_counts[i]);

      ds_aux_accept_recv(sched_a2av->aux, sched_a2av->local_recv_proc_id, sched->ds->locals.dst_counts[i]);
    }
  }

  /* accept sends */
  for (i = 0; i < sched->ds->sends.n; ++i)
  {
    Z_TRACE_IF(DS_TRACE_IF, "accept send of %" dspint_fmt " to %" dspint_fmt " @ *", sched->ds->sends.counts[i], sched->ds->sends.proc_ids[i]);

    ds_accept_send(sched_a2av, sched->ds->sends.buf_ids[i], sched->ds->sends.proc_ids[i], sched->ds->sends.counts[i]);
  }

  /* accept receives */
  for (i = 0; i < sched->ds->recvs.n; ++i)
  {
    if (sched->ds->recvs.buf_ids[i] == sched->recv.buf_id)
    {
      Z_TRACE_IF(DS_TRACE_IF, "accept receive of %" dspint_fmt " from % "dspint_fmt " @ *", sched->ds->recvs.counts[i], sched->ds->recvs.proc_ids[i]);

      ds_accept_recv(sched_a2av, sched->ds->recvs.buf_ids[i], sched->ds->recvs.proc_ids[i], sched->ds->recvs.counts[i]);

    } else
    {
      Z_TRACE_IF(DS_TRACE_IF, "accept aux receive of %" dspint_fmt " from %" dspint_fmt " @ *", sched->ds->recvs.counts[i], sched->ds->recvs.proc_ids[i]);

      ds_aux_accept_recv(sched_a2av->aux, sched->ds->recvs.proc_ids[i], sched->ds->recvs.counts[i]);
    }
  }

  for (b = 0; b < sched->nbufs; ++b)
  {
#ifdef DASH_SCHED_A2AV_OVERLAP
    /* resort the splitted send-blocks */
    ds_resort_sblocks(sched_a2av, b);

    /* create matchings */
    ds_recreate_matchings(sched_a2av, b);
#endif

    /* compactify */
    ds_compactify_blocks(sched, b);
  }

#ifdef DS_PRINT_BLOCKS
  ds_print_blocks(sched);
#endif

  /* move from aux to the receive buffer */
  sched_a2av->aux->vacate_aux = ds_vacate_aux;
  ds_aux_vacate(sched_a2av->aux, sched);

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


static dsint_t ds_sched_a2av_pre(ds_sched_t *sched)
{
#ifdef DASH_SYMMETRIC
  DEFINE_SCHED_A2AV(sched, sched_a2av);
#endif


  Z_TRACE_IF(DS_TRACE_IF, "START: round = %" dsint_fmt, sched->round);

#ifdef DS_PRINT_BLOCKS
  ds_print_blocks(sched);
#endif

#ifdef DASH_SYMMETRIC
  if (sched->round == 0 && !sched_a2av->skip_sym)
  {
    Z_TRACE_IF(DS_TRACE_IF, "symmetric");
    ds_sched_a2av_pre_symmetric(sched);

  } else
#endif
  {
    Z_TRACE_IF(DS_TRACE_IF, "individual");
    ds_sched_a2av_pre_individual(sched);
  }

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


static dsint_t ds_sched_a2av_post(ds_sched_t *sched)
{
#ifdef DASH_SYMMETRIC
  DEFINE_SCHED_A2AV(sched, sched_a2av);
#endif


  Z_TRACE_IF(DS_TRACE_IF, "START: round = %" dsint_fmt, sched->round);

#ifdef DASH_SYMMETRIC
  if (sched->round == 0 && !sched_a2av->skip_sym)
  {
    Z_TRACE_IF(DS_TRACE_IF, "symmetric");
    ds_sched_a2av_post_symmetric(sched);

  } else
#endif
  {
    Z_TRACE_IF(DS_TRACE_IF, "individual");
    ds_sched_a2av_post_individual(sched);
  }

  ++sched->round;

  return 0;
}

#undef SKIP_SYMMETRIC


static void ds_init_blocks(ds_sched_t *sched)
{
#ifdef DASH_SCHED_A2AV_OVERLAP
  DEFINE_SCHED_A2AV(sched, sched_a2av);
#endif

  dsint_t b;


  Z_TRACE_IF(DS_TRACE_IF, "START");

  ds_insert_blocks(sched);

  for (b = 0; b < sched->nbufs; ++b)
  {
#ifdef DASH_SCHED_A2AV_OVERLAP
    /* resort the splitted send-blocks */
    ds_resort_sblocks(sched_a2av, b);

    /* create matchings */
    ds_recreate_matchings(sched_a2av, b);
#endif

#ifndef DASH_SYMMETRIC
    /* omit initial compactify, because it may destroy existing initial symmetric opportunities */
    ds_compactify_blocks(sched, b);
#endif
  }

#ifdef DS_PRINT_BLOCKS
  ds_print_blocks(sched);
#endif

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


static void ds_free_blocks(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  dsint_t b;


  Z_TRACE_IF(DS_TRACE_IF, "START");

  for (b = 0; b < sched->nbufs; ++b)
  {
    if (sched_a2av->bufs[b].sblocks) z_free(sched_a2av->bufs[b].sblocks);
    if (sched_a2av->bufs[b].rblocks) z_free(sched_a2av->bufs[b].rblocks);
  }

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


static void ds_print_blocks(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  dsint_t b;
  block_t *blk;

  Z_TRACE_IF(DS_TRACE_IF, "START");

  for (b = 0; b < sched->nbufs; ++b)
  {
    Z_TRACE_IF(DS_TRACE_IF, "buf: %" dsint_fmt, b);

    blk = sched_a2av->bufs[b].sblock_first;
    while (blk != NULL)
    {
      Z_TRACE_IF(DS_TRACE_IF, " s %" dsint_fmt ": %p, " PRINT_BLOCK_STR, blk->proc_id, blk, PRINT_BLOCK_PARAMS(*blk));
      blk = blk->next;
    }

    blk = sched_a2av->bufs[b].rblock_first;
    while (blk != NULL)
    {
      Z_TRACE_IF(DS_TRACE_IF, " r %" dsint_fmt ": %p, " PRINT_BLOCK_STR, blk->proc_id, blk, PRINT_BLOCK_PARAMS(*blk));
      blk = blk->next;
    }
  }

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


#ifdef DASH_SCHED_A2AV_OVERLAP

static void ds_insert_blocks(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  dsint_t b, i, j;
  dsint_t *order;
  block_t *blk, *rblk;


  Z_TRACE_IF(DS_TRACE_IF, "START");

  for (b = 0; b < sched->nbufs; ++b)
  {
    sched_a2av->bufs[b].nsblocks = sched_a2av->bufs[b].nrblocks = 0;
    sched_a2av->bufs[b].sblocks = sched_a2av->bufs[b].rblocks = NULL;
    sched_a2av->bufs[b].sblock_first = sched_a2av->bufs[b].rblock_first = NULL;
  }

  order = z_alloca(z_max(sched->send.n, sched->recv.n), sizeof(dsint_t));

  sched_a2av->stotal = sched_a2av->rtotal = 0;

  /* create recv-blocks */
  b = sched->recv.buf_id;

  ds_sort_displs(sched->recv.n, sched->recv.displs, NULL, order);

  Z_TRACE_IF(DS_TRACE_IF, "adding %" dsint_fmt " rblocks", sched->recv.n);
  sched_a2av->bufs[b].nrblocks = sched->recv.n;
  sched_a2av->bufs[b].rblocks = z_alloc(sched->recv.n, sizeof(block_t *));

  blk = NULL;
  for (i = 0; i < sched->recv.n; ++i)
  {
    j = order[i];

    if (sched->recv.counts[j] <= 0)
    {
      sched_a2av->bufs[b].rblocks[j] = NULL;
      continue;
    }

    sched_a2av->bufs[b].rblocks[j] = blk = ds_add_block(sched_a2av, sched->recv.displs[j], sched->recv.counts[j], j, blk, b, 0);

    Z_TRACE_IF(DS_TRACE_IF, "add new recv-block %p [%" dsint_fmt ",%" dsint_fmt "]", blk, blk->begin, blk->end);

    sched_a2av->rtotal += sched->recv.counts[j];
  }

  /* create send-blocks */
  b = sched->send.buf_id;

  ds_sort_displs(sched->send.n, sched->send.displs, sched->send.counts, order);

  Z_TRACE_IF(DS_TRACE_IF, "adding %" dsint_fmt " sblocks", sched->send.n);
  sched_a2av->bufs[b].nsblocks = sched->send.n;
  sched_a2av->bufs[b].sblocks = z_alloc(sched->send.n, sizeof(block_t *));

  rblk = sched_a2av->bufs[b].rblock_first;

  blk = NULL;
  for (i = 0; i < sched->send.n; ++i)
  {
    j = order[i];

    if (sched->send.counts[j] <= 0)
    {
      sched_a2av->bufs[b].sblocks[j] = NULL;
      continue;
    }

    sched_a2av->bufs[b].sblocks[j] = blk = ds_add_block(sched_a2av, sched->send.displs[j], sched->send.counts[j], j, blk, b, 1);

    sched_a2av->stotal += sched->send.counts[j];

    /* skip recv-blocks that end before the current send-block */
    while (rblk != NULL && rblk->end <= blk->end) rblk = rblk->next;

    Z_TRACE_IF(DS_TRACE_IF, "try split new send-block %p [%" dsint_fmt ",%" dsint_fmt "] at current recv-block %p",
      blk, blk->begin, blk->end, rblk);

    /* split send-block that encloses the begin of the current recv-block */
    if (rblk != NULL && blk->begin < rblk->begin && rblk->begin < blk->end) blk = ds_split_block(sched_a2av, blk, rblk->begin);

    Z_TRACE_IF(DS_TRACE_IF, "current send-block after split: %p [%" dsint_fmt ",%" dsint_fmt "]",
      blk, blk->begin, blk->end);
  }

  z_freea(order);

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


static void ds_resort_sblocks(ds_sched_a2av_t *sched_a2av, dsint_t b)
{
  block_t **blks, *sb;
  dsint_t i, n;


  blks = z_alloca(sched_a2av->bufs[b].nsblocks * 2, sizeof(block_t *));

  sb = sched_a2av->bufs[b].sblock_first;
  n = 0;

  while (sb != NULL)
  {
    blks[n] = sb;
    ++n;
    sb = sb->next;
  }

  if (n <= 0) return;

  ds_sort_blocks_begin(n, blks);

  sched_a2av->bufs[b].sblock_first = blks[0];

  sb = NULL;
  for (i = 0; i < n - 1; ++i)
  {
    blks[i]->prev = sb;
    blks[i]->next = blks[i + 1];

    sb = blks[i];
  }
  blks[n - 1]->prev = sb;
  blks[n - 1]->next = NULL;
}


static void ds_recreate_matchings(ds_sched_a2av_t *sched_a2av, dsint_t b)
{
  block_t *sb, *rb;


  sb = sched_a2av->bufs[b].sblock_first;
  if (sb) sb->match = NULL;

  rb = sched_a2av->bufs[b].rblock_first;
  if (rb) rb->match = NULL;

  while (sb != NULL && rb != NULL)
  {
    Z_TRACE_IF(DS_TRACE_IF, "send-block: %p, " PRINT_BLOCK_STR, sb, PRINT_BLOCK_PARAMS(*sb));
    Z_TRACE_IF(DS_TRACE_IF, "recv-block: %p, " PRINT_BLOCK_STR, rb, PRINT_BLOCK_PARAMS(*rb));

    if (sb->end <= rb->begin)
    {
      Z_TRACE_IF(DS_TRACE_IF, "skip send-block %p", sb);
      sb = sb->next;
      if (sb) sb->match = NULL;
      continue;
    }

    if (rb->end <= sb->begin)
    {
      Z_TRACE_IF(DS_TRACE_IF, "skip recv-block %p", rb);
      rb = rb->next;
      if (rb) rb->match = NULL;
      continue;
    }

    if (sb->match == NULL) sb->match = rb;
    if (rb->match == NULL) rb->match = sb;

    if (sb->end <= rb->end)
    {
      sb = sb->next;
      if (sb) sb->match = NULL;

    } else
    {
      rb = rb->next;
      if (rb) rb->match = NULL;
    }
  }

  /* unset send-block remaining matchings */
  while (sb != NULL)
  {
    sb = sb->next;
    if (sb) sb->match = NULL;
  }

  /* unset recv-block remaining matchings */
  while (rb != NULL)
  {
    rb = rb->next;
    if (rb) rb->match = NULL;
  }
}

#else

static void ds_split_sblocks(ds_sched_a2av_t *sched_a2av, dsint_t b)
{
  block_t *sb, *rb, *nb;


  sb = sched_a2av->bufs[b].sblock_first;
  rb = sched_a2av->bufs[b].rblock_first;

  while (sb != NULL && rb != NULL)
  {
    Z_TRACE_IF(DS_TRACE_IF, "comparing: %p vs. %p", sb, rb);

    Z_TRACE_IF(DS_TRACE_IF, "send-block: %p, " PRINT_BLOCK_STR, sb, PRINT_BLOCK_PARAMS(*sb));
    Z_TRACE_IF(DS_TRACE_IF, "recv-block: %p, " PRINT_BLOCK_STR, rb, PRINT_BLOCK_PARAMS(*rb));

    if (sb->end <= rb->begin)
    {
      Z_TRACE_IF(DS_TRACE_IF, "skip send-block %p", sb);
      sb = sb->next;
      continue;
    }

    if (rb->end <= sb->begin)
    {
      Z_TRACE_IF(DS_TRACE_IF, "skip recv-block %p", rb);
      rb = rb->next;
      continue;
    }

    Z_TRACE_IF(DS_TRACE_IF, "intersect: send-block %p with recv-block %p", sb, rb);

    Z_TRACE_IF(DS_TRACE_IF, "split send-block %p at begin %" dsint_fmt " of %p", sb, rb->begin, rb);
    nb = ds_split_block(sched_a2av, sb, rb->begin);
    if (nb)
    {
      Z_TRACE_IF(DS_TRACE_IF, "split done, nb = %p becomes new sb", nb);
      sb = nb;
    }

    Z_TRACE_IF(DS_TRACE_IF, "set match: send-block %p -> recv-block %p", sb, rb);
    sb->match = rb;

    Z_TRACE_IF(DS_TRACE_IF, "set match: recv-block %p -> send-block %p %s", rb, sb, (!rb->match)?"YES":"NO");
    if (!rb->match) rb->match = sb;

    Z_TRACE_IF(DS_TRACE_IF, "split send-block %p at end %" dsint_fmt " of %p", sb, rb->end, rb);
    nb = ds_split_block(sched_a2av, sb, rb->end);
    if (nb)
    {
      Z_TRACE_IF(DS_TRACE_IF, "split done, nb = %p becomes new sb", nb);
      sb = nb;

    } else
    {
      Z_TRACE_IF(DS_TRACE_IF, "split failed, next %p becomes new sb", sb->next);
      sb = sb->next;
    }
  }
}


static void ds_insert_blocks(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  dsint_t b, i, j;
  dsint_t *order;
  block_t *blk;


  Z_TRACE_IF(DS_TRACE_IF, "START");

  for (b = 0; b < sched->nbufs; ++b)
  {
    sched_a2av->bufs[b].nsblocks = sched_a2av->bufs[b].nrblocks = 0;
    sched_a2av->bufs[b].sblocks = sched_a2av->bufs[b].rblocks = NULL;
    sched_a2av->bufs[b].sblock_first = sched_a2av->bufs[b].rblock_first = NULL;
  }

  order = z_alloca(z_max(sched->send.n, sched->recv.n), sizeof(dsint_t));

  sched_a2av->stotal = sched_a2av->rtotal = 0;

  /* create send-blocks */
  b = sched->send.buf_id;

  ds_sort_displs(sched->send.n, sched->send.displs, NULL, order);

  Z_TRACE_IF(DS_TRACE_IF, "adding %" dsint_fmt " sblocks", sched->send.n);
  sched_a2av->bufs[b].nsblocks = sched->send.n;
  sched_a2av->bufs[b].sblocks = z_alloc(sched->send.n, sizeof(block_t *));

  blk = NULL;
  for (i = 0; i < sched->send.n; ++i)
  {
    j = order[i];
    if (sched->send.counts[j] > 0) sched_a2av->bufs[b].sblocks[j] = blk = ds_add_block(sched_a2av, sched->send.displs[j], sched->send.counts[j], j, blk, b, 1);
    else sched_a2av->bufs[b].sblocks[j] = NULL;

    sched_a2av->stotal += sched->send.counts[j];
  }

  /* create recv-blocks */
  b = sched->recv.buf_id;

  ds_sort_displs(sched->recv.n, sched->recv.displs, NULL, order);

  Z_TRACE_IF(DS_TRACE_IF, "adding %" dsint_fmt " rblocks", sched->recv.n);
  sched_a2av->bufs[b].nrblocks = sched->recv.n;
  sched_a2av->bufs[b].rblocks = z_alloc(sched->recv.n, sizeof(block_t *));

  blk = NULL;
  for (i = 0; i < sched->recv.n; ++i)
  {
    j = order[i];
    if (sched->recv.counts[j] > 0) sched_a2av->bufs[b].rblocks[j] = blk = ds_add_block(sched_a2av, sched->recv.displs[j], sched->recv.counts[j], j, blk, b, 0);
    else sched_a2av->bufs[b].rblocks[j] = NULL;

    sched_a2av->rtotal += sched->recv.counts[j];
  }

  z_freea(order);

#ifdef DS_PRINT_BLOCKS
  ds_print_blocks(sched);
#endif

  for (b = 0; b < sched->nbufs; ++b)
  {
    ds_split_sblocks(sched_a2av, b);
  }

  Z_TRACE_IF(DS_TRACE_IF, "END");
}

#endif


#ifdef DASH_SYMMETRIC

static void ds_intersect_block_at(dsint_t begin, dsint_t end, dsint_t at_begin, dsint_t at_end, dsint_t *l_begin, dsint_t *l_end, dsint_t *m_begin, dsint_t *m_end, dsint_t *h_begin, dsint_t *h_end)
{
  if (l_begin) *l_begin = begin;
  if (l_end)   *l_end = z_min(end, at_begin);

  if (m_begin) *m_begin = z_max(begin, at_begin);
  if (m_end)   *m_end = z_min(end, at_end);

  if (h_begin) *h_begin = z_max(begin, at_end);
  if (h_end)   *h_end = end;
}

#endif


#ifndef DASH_SCHED_A2AV_OVERLAP

static void ds_move_sblock(ds_sched_t *sched, block_t *sblk, dsint_t move)
{
#ifdef DASH_SYMMETRIC
  block_t *rblk;

  dsint_t move_skip, old_sblk_begin, sblk_unset;
  dsint_t S_begin, S_end, l_begin, l_end, l_mid, m_begin, m_end, h_begin, h_end;
#endif


  Z_TRACE_IF(DS_TRACE_IF, "move sblock %p", sblk);

#ifdef DASH_SYMMETRIC
  rblk = sblk->match;

  if (sblk->sym_count > 0)
  {
    S_begin = sblk->sym_displ;
    S_end = S_begin + sblk->sym_count;

  } else if (rblk->sym_count > 0)
  {
    S_begin = rblk->sym_displ;
    S_end = S_begin + rblk->sym_count;

  } else S_begin = S_end = sblk->begin;

  if (sblk->sym_count > 0) Z_ASSERT(S_begin == rblk->sym_displ && S_end == S_begin + rblk->sym_count);

  /* current end is not inside sym */
  Z_ASSERT(sblk->end <= S_begin || S_end <= sblk->end);

  /* moved end is not inside sym */
  Z_ASSERT(sblk->end + move <= S_begin || S_end <= sblk->end + move);

  /* reduce move if current end is left and moved end is right (of active sym) */
  if (sblk->end <= S_begin && sblk->end + move >= S_end) move -= S_end - S_begin;

  move_skip = move + (S_end - S_begin);

  old_sblk_begin = sblk->begin;

  if (move == 0)
  {
    Z_TRACE_IF(DS_TRACE_IF, "move = 0 after skip over sym");
    return;
  }

  ds_intersect_block_at(sblk->begin, sblk->end, S_begin, S_end, &l_begin, &l_end, &m_begin, &m_end, &h_begin, &h_end);

  Z_TRACE_IF(DS_TRACE_IF, "intersection at [%" dsint_fmt ",%" dsint_fmt "]: low: [%" dsint_fmt ",%" dsint_fmt "], mid: [%" dsint_fmt ",%" dsint_fmt "], high: [%" dsint_fmt ",%" dsint_fmt "]: ",
    S_begin, S_end, l_begin, l_end, m_begin, m_end, h_begin, h_end);

  sblk_unset = 1;

  if (h_begin < h_end)
  {
    Z_TRACE_IF(DS_TRACE_IF, "high: move from %" dsint_fmt " to %" dsint_fmt, h_begin, h_begin + move);

    sched->ds->exec->move(sched->ds->exec, sched->recv.exec_id, sched->recv.buf_id, h_begin, sched->recv.buf_id, h_begin + move, h_end - h_begin);

    sblk->begin = h_begin + move;
    sblk->end = h_end + move;

    sblk_unset = 0;

    Z_TRACE_IF(DS_TRACE_IF, " new sblk: [%" dsint_fmt ",%" dsint_fmt "]", sblk->begin, sblk->end);
  }

  if (l_begin < l_end)
  {
    l_mid = z_minmax(l_begin, S_begin - move, l_end);

    Z_TRACE_IF(DS_TRACE_IF, "low: [%" dsint_fmt ",%" dsint_fmt ",%" dsint_fmt "]", l_begin, l_mid, l_end);

    if (l_mid < l_end)
    {
      Z_TRACE_IF(DS_TRACE_IF, " move 1st from %" dsint_fmt " to %" dsint_fmt, l_mid, l_mid + move_skip);

      sched->ds->exec->move(sched->ds->exec, sched->recv.exec_id, sched->recv.buf_id, l_mid, sched->recv.buf_id, l_mid + move_skip, l_end - l_mid);

      if (sblk_unset || l_mid + move_skip < sblk->begin) sblk->begin = l_mid + move_skip;
      if (sblk_unset || l_end + move_skip > sblk->end) sblk->end = l_end + move_skip;

      sblk_unset = 0;

      Z_TRACE_IF(DS_TRACE_IF, " new sblk: [%" dsint_fmt ",%" dsint_fmt "]", sblk->begin, sblk->end);
    }

    if (l_begin < l_mid)
    {
      Z_TRACE_IF(DS_TRACE_IF, " move 2nd from %" dsint_fmt " to %" dsint_fmt, l_begin, l_begin + move);

      sched->ds->exec->move(sched->ds->exec, sched->recv.exec_id, sched->recv.buf_id, l_begin, sched->recv.buf_id, l_begin + move, l_mid - l_begin);

      if (sblk_unset || l_begin + move < sblk->begin) sblk->begin = l_begin + move;
      if (sblk_unset || l_mid + move > sblk->end) sblk->end = l_mid + move;

      sblk_unset = 0;

      Z_TRACE_IF(DS_TRACE_IF, " new sblk: [%" dsint_fmt ",%" dsint_fmt "]", sblk->begin, sblk->end);
    }
  }

  Z_ASSERT(!sblk_unset);

  /* how far has the begin moved this time? */
  sblk->moved += sblk->begin - old_sblk_begin;

  /* if no active sym */
  if (sblk->sym_count <= 0)
  {
    /* if matching rblock has an active sym and it is enclosed by sblock */
    if (rblk->sym_count > 0 && sblk->begin < S_begin && S_end <= sblk->end)
    {
      /* set passive sym */
      sblk->sym_count = -rblk->sym_count;
      sblk->sym_displ = rblk->sym_displ;
      sblk->sym_offset = sblk->sym_displ - sblk->begin;

    } else sblk->sym_count = sblk->sym_displ = sblk->sym_offset = 0; /* unset passive sym */
  }

#else
  sched->ds->exec->move(sched->ds->exec, sched->recv.exec_id, sched->recv.buf_id, sblk->begin, sched->recv.buf_id, sblk->begin + move, sblk->end - sblk->begin);

  sblk->begin += move;
  sblk->end += move;
#endif
}

#else

static block_t *ds_move_sblock_group(ds_sched_t *sched, block_t *last, dsint_t group_size, dsint_t group_begin, dsint_t group_end, dsint_t move
#ifdef DASH_SYMMETRIC
  , dsint_t S_begin, dsint_t S_end
#endif
  )
{
  dsint_t i, sblk_begin, sblk_end;
#ifdef DASH_SYMMETRIC
  dsint_t sblk_begin_old, sblk_unset, move_straight, move_skip;
  dsint_t l_begin, l_end, l_mid, m_begin, m_end, h_begin, h_end;
#endif


  Z_TRACE_IF(DS_TRACE_IF, "move-group: last: %p, group-size: %" dsint_fmt ", [%" dsint_fmt ",%" dsint_fmt "] by %" dsint_fmt, last, group_size, group_begin, group_end, move);
#ifdef DASH_SYMMETRIC
  Z_TRACE_IF(DS_TRACE_IF, "  with sym: [%" dsint_fmt ",%" dsint_fmt "]", S_begin, S_end);
#endif

  if (move <= 0)
  {
    for (i = 0; i < group_size; ++i) last = last->prev;

    return last;
  }

  for (i = -1; i < group_size; ++i)
  {
    if (i >= 0)
    {
      Z_TRACE_IF(DS_TRACE_IF, "move sblock %p by %" dsint_fmt, last, move);

      sblk_begin = last->begin;
      sblk_end = last->end;

    } else
    {
      Z_TRACE_IF(DS_TRACE_IF, "move content [%" dsint_fmt ",%" dsint_fmt "] by %" dsint_fmt, group_begin, group_end, move);

      sblk_begin = group_begin;
      sblk_end = group_end;
    }

#ifdef DASH_SYMMETRIC
    /* current end is not inside sym */
    Z_ASSERT(sblk_end <= S_begin || S_end <= sblk_end);

    move_straight = move;

#if 0
    /* moved end is not inside sym */
    Z_ASSERT(sblk_end + move <= S_begin || S_end <= sblk_end + move);
#endif

    /* reduce move if current end is left and moved end is right (of active sym) */
/*    if (sblk_end <= S_begin && sblk_end + move >= S_end) move_straight -= S_end - S_begin;*/
    if (group_end <= S_begin && group_end + move >= S_end) move_straight -= S_end - S_begin;

    if (move_straight == 0)
    {
      Z_TRACE_IF(DS_TRACE_IF, "move straight equal 0 after skip over sym");
      goto move_next_block;
    }

    sblk_begin_old = sblk_begin;

    sblk_unset = 1;

    move_skip = move_straight + (S_end - S_begin);

    Z_TRACE_IF(DS_TRACE_IF, "straight move: %" dsint_fmt ", skip move: %" dsint_fmt, move_straight, move_skip);

    ds_intersect_block_at(sblk_begin, sblk_end, S_begin, S_end, &l_begin, &l_end, &m_begin, &m_end, &h_begin, &h_end);

    Z_TRACE_IF(DS_TRACE_IF, "intersection at [%" dsint_fmt ",%" dsint_fmt "]: low: [%" dsint_fmt ",%" dsint_fmt "], mid: [%" dsint_fmt ",%" dsint_fmt "], high: [%" dsint_fmt ",%" dsint_fmt "]: ",
      S_begin, S_end, l_begin, l_end, m_begin, m_end, h_begin, h_end);

    if (h_begin < h_end)
    {
      Z_TRACE_IF(DS_TRACE_IF, "high: move from %" dsint_fmt " to %" dsint_fmt, h_begin, h_begin + move_straight);

      if (i < 0) sched->ds->exec->move(sched->ds->exec, sched->recv.exec_id, sched->recv.buf_id, h_begin, sched->recv.buf_id, h_begin + move_straight, h_end - h_begin);

      sblk_begin = h_begin + move_straight;
      sblk_end = h_end + move_straight;

      sblk_unset = 0;

      Z_TRACE_IF(DS_TRACE_IF, " new sblk: [%" dsint_fmt ",%" dsint_fmt "]", sblk_begin, sblk_end);
    }

    if (l_begin < l_end)
    {
      l_mid = z_minmax(l_begin, S_begin - move_straight, l_end);

      Z_TRACE_IF(DS_TRACE_IF, "low: [%" dsint_fmt ",%" dsint_fmt ",%" dsint_fmt "]", l_begin, l_mid, l_end);

      if (l_mid < l_end)
      {
        Z_TRACE_IF(DS_TRACE_IF, " move 1st from %" dsint_fmt " to %" dsint_fmt, l_mid, l_mid + move_skip);

        if (i < 0) sched->ds->exec->move(sched->ds->exec, sched->recv.exec_id, sched->recv.buf_id, l_mid, sched->recv.buf_id, l_mid + move_skip, l_end - l_mid);

        if (sblk_unset || l_mid + move_skip < sblk_begin) sblk_begin = l_mid + move_skip;
        if (sblk_unset || l_end + move_skip > sblk_end) sblk_end = l_end + move_skip;

        sblk_unset = 0;

        Z_TRACE_IF(DS_TRACE_IF, " new sblk: [%" dsint_fmt ",%" dsint_fmt "]", sblk_begin, sblk_end);
      }

      if (l_begin < l_mid)
      {
        Z_TRACE_IF(DS_TRACE_IF, " move 2nd from %" dsint_fmt " to %" dsint_fmt, l_begin, l_begin + move_straight);

        if (i < 0) sched->ds->exec->move(sched->ds->exec, sched->recv.exec_id, sched->recv.buf_id, l_begin, sched->recv.buf_id, l_begin + move_straight, l_mid - l_begin);

        if (sblk_unset || l_begin + move_straight < sblk_begin) sblk_begin = l_begin + move_straight;
        if (sblk_unset || l_mid + move_straight > sblk_end) sblk_end = l_mid + move_straight;

        sblk_unset = 0;

        Z_TRACE_IF(DS_TRACE_IF, " new sblk: [%" dsint_fmt ",%" dsint_fmt "]", sblk_begin, sblk_end);
      }
    }

    Z_ASSERT(!sblk_unset);

    if (i >= 0)
    {
      /* set moved extent */
      last->begin = sblk_begin;
      last->end = sblk_end;

      /* how far has the begin moved this time? */
      last->moved += sblk_begin - sblk_begin_old;

      /* if no active sym */
      if (last->sym_count <= 0)
      {
        /* if matching rblock has an active sym and it is enclosed by sblock */
        if (last->match->sym_count > 0 && last->begin < S_begin && S_end <= last->end)
        {
          /* set passive sym */
          last->sym_count = -last->match->sym_count;
          last->sym_displ = last->match->sym_displ;
          last->sym_offset = last->sym_displ - last->begin;

        } else last->sym_count = last->sym_displ = last->sym_offset = 0; /* unset passive sym */
      }
    }

move_next_block:

#else
    if (i >= 0)
    {
      last->begin += move;
      last->end += move;

    } else
    {
      sched->ds->exec->move(sched->ds->exec, sched->recv.exec_id, sched->recv.buf_id, sblk_begin, sched->recv.buf_id, sblk_begin + move, sblk_end - sblk_begin);
    }
#endif

    if (i >= 0) last = last->prev;
  }

  return last;
}

#endif


#ifdef DASH_SCHED_A2AV_OVERLAP

static void ds_compactify_blocks(ds_sched_t *sched, dsint_t b)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  block_t *sblk_last, *rblk;
  dsint_t end, move;
  dsint_t nmove_groups, *mg_sizes, *mg_begins, *mg_ends;


  Z_TRACE_IF(DS_TRACE_IF, "START");

  mg_sizes  = z_alloca(3 * sched->send.n, sizeof(dsint_t));
  mg_begins = mg_sizes + 1 * sched->send.n;
  mg_ends   = mg_sizes + 2 * sched->send.n;

  rblk = sched_a2av->bufs[b].rblock_first;

  sblk_last = NULL;

  while (rblk != NULL)
  {
    Z_TRACE_IF(DS_TRACE_IF, "proc: %" dsint_fmt ", recv-block: %p: " PRINT_BLOCK_STR, rblk->proc_id, rblk, PRINT_BLOCK_PARAMS(*rblk));

    if (rblk->match == NULL) sblk_last = NULL;

    if (sblk_last != NULL && sblk_last->end <= rblk->begin) sblk_last = NULL;

    if (sblk_last == NULL) sblk_last = rblk->match;

    if (sblk_last != NULL)
    {
      nmove_groups = 0;

      end = sblk_last->begin - 1;

      while (1)
      {
        /* current last and next last NOT connected or overlapping? */
        if (end < sblk_last->begin)
        {
          /* start new move-group */
          ++nmove_groups;
          mg_sizes[nmove_groups - 1]  = 1;
          mg_begins[nmove_groups - 1] = sblk_last->begin;
          mg_ends[nmove_groups - 1]   = sblk_last->end;

        } else
        {
          /* extent current move-group */
          ++mg_sizes[nmove_groups - 1];
          mg_begins[nmove_groups - 1] = z_min(mg_begins[nmove_groups - 1], sblk_last->begin);
          mg_ends[nmove_groups - 1]   = z_max(mg_ends[nmove_groups - 1], sblk_last->end);
        }

        if (rblk->end <= sblk_last->end) break;

        if (sblk_last->next == NULL || sblk_last->next->begin >= rblk->end) break;

        end = mg_ends[nmove_groups - 1];

        sblk_last = sblk_last->next;
      }

      end = rblk->end;

      Z_TRACE_IF(DS_TRACE_IF, "compactify: match: %p, last: %p, end: %" dsint_fmt ", move-groups: %" dsint_fmt, rblk->match, sblk_last, end, nmove_groups);

      /* go back to the first block */
      while (nmove_groups > 0)
      {
        --nmove_groups;

        /* distance of current move-group end to the current end (if move-group end spans across the blk end (thus, move is less than zero), then nothing should be moved) */
        move = z_max(0, end - mg_ends[nmove_groups]);

        Z_TRACE_IF(DS_TRACE_IF, "move-group %" dsint_fmt ": %" dsint_fmt ", %" dsint_fmt ", %" dsint_fmt ", end: %" dsint_fmt ", move: %" dsint_fmt,
          nmove_groups, mg_sizes[nmove_groups], mg_begins[nmove_groups], mg_ends[nmove_groups], end, move);

        sblk_last = ds_move_sblock_group(sched, sblk_last, mg_sizes[nmove_groups], mg_begins[nmove_groups], mg_ends[nmove_groups], move
#ifdef DASH_SYMMETRIC
          , sblk_last->match->sym_displ, sblk_last->match->sym_displ + sblk_last->match->sym_count
#endif
          );

        end = mg_begins[nmove_groups];
      }
    }

    rblk = rblk->next;
  }

  z_freea(mg_sizes);

  Z_TRACE_IF(DS_TRACE_IF, "END");
}

#else

static void ds_compactify(ds_sched_t *sched, block_t *blk)
{
  block_t *first, *last, *current;
  dsint_t end, move;

  first = last = blk->match;

  if (!first) return;

  /* go to the last block */
  while (last->next && last->next->match == blk) last = last->next;

  end = blk->end;

  Z_TRACE_IF(DS_TRACE_IF, "compactify: first: %p, last: %p, end: %" dsint_fmt, first, last, end);

  /* go back to the first block */
  do
  {
    current = last;
    last = last->prev;

    Z_TRACE_IF(DS_TRACE_IF, "current: %p [%" dsint_fmt ",%" dsint_fmt "]", current, current->begin, current->end);

    /* current block can not be empty! */
    Z_ASSERT(current->begin < current->end);

    /* distance of current block to the current end */
    move = end - current->end;

    if (move > 0)
    {
      Z_TRACE_IF(DS_TRACE_IF, "moving %p from %" dsint_fmt " by %" dsint_fmt " (count %" dsint_fmt ")", current, current->begin, move, current->end - current->begin);

      /* move content */
      ds_move_sblock(sched, current, move);
    }

    end = current->begin;

  } while (current != first);
}

static void ds_compactify_blocks(ds_sched_t *sched, dsint_t b)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  block_t *rblk;

  Z_TRACE_IF(DS_TRACE_IF, "START");

  /* compactify */
#if 0
  dsint_t i;
  for (i = 0; i < sched_a2av->bufs[b].nrblocks; ++i)
  if (sched_a2av->bufs[b].rblocks[i]) ds_compactify(sched, sched_a2av->bufs[b].rblocks[i]);
#else
  rblk = sched_a2av->bufs[b].rblock_first;
  while (rblk != NULL)
  {
    ds_compactify(sched, rblk);

    rblk = rblk->next;
  }
#endif

  Z_TRACE_IF(DS_TRACE_IF, "END");
}

#endif


static block_t *ds_add_block(ds_sched_a2av_t *sched_a2av, dspint_t displ, dspint_t count, dspint_t proc_id, block_t *prev, dsint_t buf_id, dsint_t send)
{
  block_t *blk = NULL;

  blk = ds_block_pool_get(&sched_a2av->bp);

  blk->begin = displ;
  blk->end = displ + count;

  blk->proc_id = proc_id;

#ifdef DASH_SYMMETRIC
  blk->moved = 0;
  blk->sym_count = blk->sym_displ = blk->sym_offset = 0;
#endif

  blk->match = NULL;
#ifdef DASH_SCHED_A2AV_OVERLAP
  blk->second = NULL;
#endif

  if (prev)
  {
    prev->next = blk;
    blk->prev = prev;

  } else
  {
    if (send) sched_a2av->bufs[buf_id].sblock_first = blk;
    else sched_a2av->bufs[buf_id].rblock_first = blk;
    blk->prev = NULL;
  }
  blk->next = NULL;

  return blk;
}


static void ds_del_block(ds_sched_a2av_t *sched_a2av, block_t *blk, dsint_t buf_id, dsint_t send)
{
  if (blk->prev) blk->prev->next = blk->next;
  if (blk->next) blk->next->prev = blk->prev;

  /* if there is a matching rblock that matches this sblock */
  if (blk->match && blk->match->match == blk)
  {
    /* update the matching rblock */
    if (blk->next && blk->next->match == blk->match) blk->match->match = blk->next;
    else blk->match->match = NULL;
  }

  if (send)
  {
    if (sched_a2av->bufs[buf_id].sblock_first == blk) sched_a2av->bufs[buf_id].sblock_first = blk->next;

    if (sched_a2av->bufs[buf_id].sblocks[blk->proc_id] == blk)
    {
#ifdef DASH_SCHED_A2AV_OVERLAP
      sched_a2av->bufs[buf_id].sblocks[blk->proc_id] = blk->second;
#else
      if (blk->next && blk->next->proc_id == blk->proc_id) sched_a2av->bufs[buf_id].sblocks[blk->proc_id] = blk->next;
      else sched_a2av->bufs[buf_id].sblocks[blk->proc_id] = NULL;
#endif
    }

  } else
  {
    if (sched_a2av->bufs[buf_id].rblock_first == blk) sched_a2av->bufs[buf_id].rblock_first = blk->next;

    sched_a2av->bufs[buf_id].rblocks[blk->proc_id] = NULL;
  }

  ds_block_pool_put(&sched_a2av->bp, blk);
}


static block_t *ds_split_block(ds_sched_a2av_t *sched_a2av, block_t *blk, dsint_t p)
{
  block_t *new_blk;


  /* don't split if 'p' is out of range */
  if (p <= blk->begin || blk->end <= p) return NULL;

  new_blk = ds_block_pool_get(&sched_a2av->bp);

  new_blk->begin = p;
  new_blk->end = blk->end;

#ifdef DASH_SYMMETRIC
  new_blk->moved = blk->moved;
  new_blk->sym_count = new_blk->sym_displ = new_blk->sym_offset = 0;
#endif

  blk->end = p;

  new_blk->proc_id = blk->proc_id;

  new_blk->match = NULL;
#ifdef DASH_SCHED_A2AV_OVERLAP
  new_blk->second = NULL;
#endif
  new_blk->prev = blk;
  new_blk->next = blk->next;

  if (new_blk->next) new_blk->next->prev = new_blk;

  blk->next = new_blk;

#ifdef DASH_SCHED_A2AV_OVERLAP
  blk->second = new_blk;
#endif

  return new_blk;
}


static void ds_check_block(ds_sched_a2av_t *sched_a2av, block_t *blk, dsint_t buf_id, dsint_t send)
{
#ifdef DASH_SYMMETRIC
  Z_TRACE_IF(DS_TRACE_IF, "blk: %p" PRINT_BLOCK_STR, blk, PRINT_BLOCK_PARAMS(*blk));

  if (blk->sym_count != 0 && blk->begin == blk->sym_displ)
  {
    Z_TRACE_IF(DS_TRACE_IF, "skip sym of |%" dsint_fmt "|", blk->sym_count);

    blk->begin += z_abs(blk->sym_count);
  }

  Z_ASSERT(blk->sym_offset >= 0);

  if (blk->sym_offset == 0)
  {
    Z_TRACE_IF(DS_TRACE_IF, "unset and free sym");

    /* take only active syms into account */
    if (blk->sym_count > 0)
    {
      if (send) sched_a2av->stotal -= blk->sym_count;
      else sched_a2av->rtotal -= blk->sym_count;
    }

    blk->sym_count = blk->sym_displ = blk->sym_offset = 0;
  }
#endif

  if (blk->begin >= blk->end)
  {
    Z_TRACE_IF(DS_TRACE_IF, "delete block %p, proc: %" dsint_fmt, blk, blk->proc_id);

    ds_del_block(sched_a2av, blk, buf_id, send);
  }
}


static void ds_sort_displs(dsint_t ndispls, dspint_t *displs, dspint_t *counts, dsint_t *order)
{
  dsint_t i, s[2 * ndispls];

  Z_TRACE_IF(DS_TRACE_IF, "START: ndispls = %" dsint_fmt, ndispls);

  if (counts)
  {
    for (i = 0; i < ndispls; ++i)
    {
      s[2 * i + 0] = displs[i] + counts[i];
      s[2 * i + 1] = i;
    }

  } else
  {
    for (i = 0; i < ndispls; ++i)
    {
      s[2 * i + 0] = displs[i];
      s[2 * i + 1] = i;
    }
  }

  ds_sort_dsints(s, ndispls, 2);

  for (i = 0; i < ndispls; ++i) order[i] = s[2 * i + 1];

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


#ifdef DASH_SCHED_A2AV_OVERLAP

/*#define QSORT*/

#ifdef QSORT
static int ds_compare_begins(const void *b0, const void *b1)
{
  return (*((block_t **) b0))->begin - (*((block_t **) b1))->begin;
}
#endif

static void ds_sort_blocks_begin(dsint_t n, block_t **blks)
{
  Z_TRACE_IF(DS_TRACE_IF, "START: n = %" dsint_fmt, n);

#ifdef QSORT

  qsort(blks, n, sizeof(block_t *), ds_compare_begins);

#else

  block_t *t;
#define GET(_i_)        blks[_i_]->begin
#define XCHG(_i_, _j_)  Z_MOP(t = blks[_i_]; blks[_i_] = blks[_j_]; blks[_j_] = t;)

  DASH_RADIXSORT_SORT_DECLARE(dsint_t)

  DASH_RADIXSORT_SORT(n, 0, sizeof(dsint_t) * 8 - 1, GET, XCHG);

#undef GET
#undef XCHG

#endif

  Z_TRACE_IF(DS_TRACE_IF, "END");
}

#undef QSORT

#endif


static void ds_block_pool_add_buffer(block_pool_t *bp, void *buf, dsint_t bufsize)
{
  dsint_t nblocks;
  block_t *blocks;

  Z_TRACE_IF(DS_TRACE_IF, "START");

  nblocks = bufsize / sizeof(block_t);
  blocks = buf;

  if (!bp->free && nblocks > 0)
  {
    blocks->prev = blocks->next = NULL;
    bp->free = blocks;
    ++bp->nfree;

    --nblocks;
    ++blocks;
  }

  while (nblocks > 0)
  {
    blocks->prev = NULL;
    blocks->next = bp->free;
    bp->free = blocks;
    ++bp->nfree;

    --nblocks;
    ++blocks;
  }

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


static void ds_block_pool_init(block_pool_t *bp, dsint_t nblocks)
{
  Z_TRACE_IF(DS_TRACE_IF, "START: nblocks = %" dsint_fmt, nblocks);

  bp->nused = 0;
  bp->nfree = 0;
  bp->free = NULL;

  if (nblocks > 0)
  {
    bp->mem = z_alloc(nblocks, sizeof(block_t));

    ds_block_pool_add_buffer(bp, bp->mem, nblocks * sizeof(block_t));

  } else bp->mem = NULL;

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


static void ds_block_pool_free(block_pool_t *bp)
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  Z_ASSERT(bp->nused == 0);

  z_free(bp->mem);

  Z_TRACE_IF(DS_TRACE_IF, "END");
}


static block_t *ds_block_pool_get(block_pool_t *bp)
{
  block_t *blk = NULL;

  if (bp->nfree > 0)
  {
    blk = bp->free;
    bp->free = blk->next;
    ++bp->nused;

    blk->next = NULL;
  }

  return blk;
}


static void ds_block_pool_put(block_pool_t *bp, block_t *blk)
{
  Z_ASSERT(bp->nused > 0);

  blk->next = bp->free;
  bp->free = blk;
  --bp->nused;
}


static void ds_req_init(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  /* FIXME: choose full_reqs collectively depending on sparse */
  sched_a2av->full_reqs = 1;

  if (sched_a2av->full_reqs) sched_a2av->max_nrecv_reqs = sched_a2av->max_nsend_reqs = sched->ds->comm_size * NREQ_MAX;
  else
  {
    sched_a2av->max_nrecv_reqs = sched->recv.n * NREQ_MAX;
    sched_a2av->max_nsend_reqs = sched->send.n * NREQ_MAX;
  }

  sched_a2av->recv_reqs = z_alloc(sched_a2av->max_nrecv_reqs, sizeof(dsint_t));
  sched_a2av->send_reqs = z_alloc(sched_a2av->max_nsend_reqs, sizeof(dsint_t));
}


static void ds_req_free(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  z_free(sched_a2av->recv_reqs);
  z_free(sched_a2av->send_reqs);
}


static void ds_req_reset_recvs(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  dsint_t i;

  for (i = 0; i < sched_a2av->max_nrecv_reqs; ++i) sched_a2av->recv_reqs[i] = 0;
}


static void ds_req_exchange(ds_sched_t *sched)
{
  DEFINE_SCHED_A2AV(sched, sched_a2av);

  if (sched_a2av->full_reqs)
  {
    MPI_Alltoall(sched_a2av->recv_reqs, NREQ_MAX, MPI_DSINT, sched_a2av->send_reqs, NREQ_MAX, MPI_DSINT, sched->ds->comm);
    return;
  }

  /* FIXME: implement sparse request exchange */
}


static void ds_req_set_recv(ds_sched_t *sched, dsint_t proc_id, dsint_t req_id, dsint_t n)
{
  dsint_t i;


  DEFINE_SCHED_A2AV(sched, sched_a2av);

  Z_ASSERT(req_id < NREQ_MAX);

  i = (sched_a2av->full_reqs && sched->recv.ranks)?sched->recv.ranks[proc_id]:proc_id;

  Z_ASSERT(i * NREQ_MAX + req_id < sched_a2av->max_nrecv_reqs);

  sched_a2av->recv_reqs[i * NREQ_MAX + req_id] = n;
}


static dsint_t ds_req_get_recv(ds_sched_t *sched, dsint_t proc_id, dsint_t req_id)
{
  dsint_t i;


  DEFINE_SCHED_A2AV(sched, sched_a2av);

  Z_ASSERT(req_id < NREQ_MAX);

  i = (sched_a2av->full_reqs && sched->recv.ranks)?sched->recv.ranks[proc_id]:proc_id;

  Z_ASSERT(i * NREQ_MAX + req_id < sched_a2av->max_nrecv_reqs);

  return sched_a2av->recv_reqs[i * NREQ_MAX + req_id];
}


static dsint_t ds_req_get_send(ds_sched_t *sched, dsint_t proc_id, dsint_t req_id)
{
  dsint_t i;


  DEFINE_SCHED_A2AV(sched, sched_a2av);

  Z_ASSERT(req_id < NREQ_MAX);

  i = (sched_a2av->full_reqs && sched->send.ranks)?sched->send.ranks[proc_id]:proc_id;

  Z_ASSERT(i * NREQ_MAX + req_id < sched_a2av->max_nsend_reqs);

  return sched_a2av->send_reqs[i * NREQ_MAX + req_id];
}


#undef DS_TRACE_IF

#undef DS_PRINT_BLOCKS



#include "dash_sched_a2av_aux.h"


#ifndef DS_TRACE_IF
# define DS_TRACE_IF  (z_mpi_rank == -1)
#endif


void ds_aux_create(ds_sched_a2av_aux_t *aux, dsint_t size) /* ds_func ds_aux_create */
{
  aux->size = size;

  Z_TRACE_IF(DS_TRACE_IF, "size: %" dsint_fmt, size);

  aux->pre_acquire = NULL;
  aux->acquire = NULL;
  aux->post_acquire = NULL;

  aux->get_count = NULL;
  aux->get_displ = NULL;

  aux->vacate = NULL;

  aux->accept_recv = NULL;
}


void ds_aux_destroy(ds_sched_a2av_aux_t *aux) /* ds_func ds_aux_destroy */
{
}


void ds_aux_pre_acquire(ds_sched_a2av_aux_t *aux) /* ds_func ds_aux_pre_acquire */
{
  if (aux->pre_acquire) aux->pre_acquire(aux);
}


dsint_t ds_aux_acquire(ds_sched_a2av_aux_t *aux, dsint_t count, dspint_t proc_id) /* ds_func ds_aux_acquire */
{
  return aux->acquire(aux, count, proc_id);
}


void ds_aux_post_acquire(ds_sched_a2av_aux_t *aux) /* ds_func ds_aux_post_acquire */
{
  if (aux->post_acquire) aux->post_acquire(aux);
}


dsint_t ds_aux_get_count(ds_sched_a2av_aux_t *aux, dspint_t proc_id) /* ds_func ds_aux_get_count */
{
  return aux->get_count(aux, proc_id);
}


dsint_t ds_aux_get_displ(ds_sched_a2av_aux_t *aux, dspint_t proc_id) /* ds_func ds_aux_get_displ */
{
  return aux->get_displ(aux, proc_id);
}


void ds_aux_vacate(ds_sched_a2av_aux_t *aux, void *schedptr) /* ds_func ds_aux_vacate */
{
  aux->vacate(aux, schedptr);
}


void ds_aux_accept_recv(ds_sched_a2av_aux_t *aux, dspint_t proc_id, dsint_t count) /* ds_func ds_aux_accept_recv */
{
  aux->accept_recv(aux, proc_id, count);
}


#undef DS_TRACE_IF



#include "dash_sched_a2av_sym.h"


#ifndef DS_TRACE_IF
# define DS_TRACE_IF  (z_mpi_rank == -1)
#endif


#ifdef DASH_SYMMETRIC


static dsint_t ds_sched_a2av_sym_max_n(ds_sched_t *sched, dsint_t *max_nsends, dsint_t *max_n, dsint_t *max_nlocals, dsint_t *max_nsyms);
/*static dsint_t ds_sched_a2av_sym_pre_run(ds_sched_t *sched);
static dsint_t ds_sched_a2av_sym_post_run(ds_sched_t *sched);*/
static dsint_t ds_sched_a2av_sym_finished(ds_sched_t *sched);
static dsint_t ds_sched_a2av_sym_pre(ds_sched_t *sched);
/*static dsint_t ds_sched_a2av_sym_post(ds_sched_t *sched);*/


dsint_t ds_sched_a2av_sym_create(ds_sched_t *sched) /* ds_func ds_sched_a2av_sym_create */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  ds_sched_create(sched);

  sched->max_n = ds_sched_a2av_sym_max_n;
  sched->finished = ds_sched_a2av_sym_finished;
  sched->pre = ds_sched_a2av_sym_pre;

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


dsint_t ds_sched_a2av_sym_destroy(ds_sched_t *sched) /* ds_func ds_sched_a2av_sym_destroy */
{
  Z_TRACE_IF(DS_TRACE_IF, "START");

  ds_sched_destroy(sched);

  Z_TRACE_IF(DS_TRACE_IF, "END");

  return 0;
}


static dsint_t ds_sched_a2av_sym_max_n(ds_sched_t *sched, dsint_t *max_nsends, dsint_t *max_nrecvs, dsint_t *max_nlocals, dsint_t *max_nsyms)
{
#ifdef DASH_SYMMETRIC_AUX
  *max_nsends = 0;
  *max_nrecvs = 2 * sched->recv.n;
#else
  *max_nsends = 0;
  *max_nrecvs = 0;
#endif
  *max_nlocals = 0;

  *max_nsyms = sched->recv.n;

  return 0;
}


static dsint_t ds_sched_a2av_sym_finished(ds_sched_t *sched)
{
  if (sched->cxt) return 1;

  sched->cxt = (void *) 1;

  return 0;
}


static dsint_t ds_sched_a2av_sym_pre(ds_sched_t *sched)
{
  dsint_t i;


  sched->ds->syms.n = sched->recv.n;
  sched->ds->syms.execute = 0;
  sched->ds->syms.exec_id = sched->send.exec_id;

  for (i = 0; i < sched->recv.n; ++i)
  {
    sched->ds->syms.buf_ids[i] = 0;
    sched->ds->syms.counts[i] = sched->recv.counts[i];
    sched->ds->syms.displs[i] = sched->recv.displs[i];
    sched->ds->syms.proc_ids[i] = i;

    if (sched->recv.counts[i]) sched->ds->syms.execute = 1;
  }

  return 0;
}

#endif /* DASH_SYMMETRIC */


#undef DS_TRACE_IF
