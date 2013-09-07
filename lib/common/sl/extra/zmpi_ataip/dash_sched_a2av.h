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


#ifndef __DASH_SCHED_A2AV_H__
#define __DASH_SCHED_A2AV_H__


#include "dash_sched.h"
#include "dash_sched_a2av_aux.h"


typedef struct _block_t
{
  dsint_t begin, end;

  dsint_t proc_id;
  
#ifdef DASH_SYMMETRIC
  dsint_t moved;
  dsint_t sym_count,   /* > 0 -> number of items in active sym (sblock and rblock), < 0 -> neg. number of items in passive sym (sblock only) */
          sym_displ,   /* original (un-moved) position where the current (active or passive) sym begins */
          sym_offset;  /* number of items between begin and active sym (sblock-own sym -> number of pre-sym items) */
#endif

  struct _block_t *match, *prev, *next;
#ifdef DASH_SCHED_A2AV_OVERLAP
  struct _block_t *second;
#endif

} block_t, *block_p;

#ifdef DASH_SYMMETRIC
# define PRINT_BLOCK_SYM_STR          ", moved: %" dsint_fmt ", sym_count: %" dsint_fmt ", sym_displ: %" dsint_fmt ", sym_offset: %" dsint_fmt
# define PRINT_BLOCK_SYM_PARAMS(_b_)  (_b_).moved, (_b_).sym_count, (_b_).sym_displ, (_b_).sym_offset, 
#else
# define PRINT_BLOCK_SYM_STR          ""
# define PRINT_BLOCK_SYM_PARAMS(_b_)
#endif

#define PRINT_BLOCK_STR          "[%" dsint_fmt ",%" dsint_fmt "], proc: %" dsint_fmt PRINT_BLOCK_SYM_STR ", match: %p, prev: %p, next: %p"
#define PRINT_BLOCK_PARAMS(_b_)  (_b_).begin, (_b_).end, (_b_).proc_id, PRINT_BLOCK_SYM_PARAMS(_b_) (_b_).match, (_b_).prev, (_b_).next

typedef struct _block_pool_t
{
  void *mem;

  dsint_t nused, nfree;
  block_t *free;

} block_pool_t;


#define NREQ_PER_PROC  2
#define NREQ_SYM       3
#define NREQ_MAX       z_max(NREQ_PER_PROC, NREQ_SYM)

typedef struct _ds_sched_a2av_t
{
  dsint_t local_send_proc_id, local_recv_proc_id;

  block_pool_t bp;

  ds_sched_a2av_aux_t *aux;

  struct {
    dsint_t nsblocks, nrblocks;
    block_t **sblocks, **rblocks;

    block_t *sblock_first, *rblock_first;

  } bufs[DASH_MAX_NBUFFERS];

  dsint_t full_reqs;
  dsint_t max_nrecv_reqs, max_nsend_reqs;
  dsint_t *recv_reqs, *send_reqs;

  dsint_t stotal, rtotal;

  dsint_t skip_sym;

} ds_sched_a2av_t, *ds_sched_a2av_p;

#define DEFINE_SCHED_A2AV(_s_, _v_)  ds_sched_a2av_t *_v_ = (_s_)->cxt


extern dsint_t ds_sched_a2av_aux_blocks;

dsint_t ds_sched_a2av_create(ds_sched_t *sched);
dsint_t ds_sched_a2av_destroy(ds_sched_t *sched);

void ds_sched_a2av_set_aux(ds_sched_t *sched, ds_sched_a2av_aux_t *aux);
void ds_sched_a2av_skip_sym(ds_sched_t *sched, dsint_t skip_sym);


#endif /* __DASH_SCHED_A2AV_H__ */
