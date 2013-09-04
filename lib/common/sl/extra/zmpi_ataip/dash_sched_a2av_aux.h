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


#ifndef __DASH_SCHED_A2AV_AUX_H__
#define __DASH_SCHED_A2AV_AUX_H__


#include "dash_core.h"


struct _ds_sched_a2av_aux_t;

typedef void ds_aux_pre_acquire_f(struct _ds_sched_a2av_aux_t *aux);
typedef dsint_t ds_aux_acquire_f(struct _ds_sched_a2av_aux_t *aux, dsint_t count, dspint_t proc_id);
typedef void ds_aux_post_acquire_f(struct _ds_sched_a2av_aux_t *aux);

typedef dsint_t ds_aux_get_count_f(struct _ds_sched_a2av_aux_t *aux, dspint_t proc_id);
typedef dsint_t ds_aux_get_displ_f(struct _ds_sched_a2av_aux_t *aux, dspint_t proc_id);

typedef void ds_aux_vacate_f(struct _ds_sched_a2av_aux_t *aux, void *schedptr);

typedef void ds_aux_accept_recv_f(struct _ds_sched_a2av_aux_t *aux, dspint_t proc_id, dsint_t count);

typedef struct _ds_sched_a2av_aux_t
{
  dsint_t size;

  dsint_t (*add_aux_recv_request)(void *, dsint_t, dsint_t);
  dsint_t (*vacate_aux)(void *, dspint_t, dsint_t, dsint_t);

  ds_aux_pre_acquire_f *pre_acquire;
  ds_aux_acquire_f *acquire;
  ds_aux_post_acquire_f *post_acquire;

  ds_aux_get_count_f *get_count;
  ds_aux_get_displ_f *get_displ;

  ds_aux_vacate_f *vacate;

  ds_aux_accept_recv_f *accept_recv;

  void *data;

} ds_sched_a2av_aux_t;


void ds_aux_create(ds_sched_a2av_aux_t *aux, dsint_t size);
void ds_aux_destroy(ds_sched_a2av_aux_t *aux);

void ds_aux_pre_acquire(ds_sched_a2av_aux_t *aux);
dsint_t ds_aux_acquire(ds_sched_a2av_aux_t *aux, dsint_t count, dspint_t proc_id);
void ds_aux_post_acquire(ds_sched_a2av_aux_t *aux);

dsint_t ds_aux_get_count(ds_sched_a2av_aux_t *aux, dspint_t proc_id);
dsint_t ds_aux_get_displ(ds_sched_a2av_aux_t *aux, dspint_t proc_id);

void ds_aux_vacate(ds_sched_a2av_aux_t *aux, void *schedptr);

void ds_aux_accept_recv(ds_sched_a2av_aux_t *aux, dspint_t proc_id, dsint_t count);


#ifdef DASH_SCHED_A2AV_AUX_STATIC
# include "dash_aux_static.h"
#endif
#ifdef DASH_SCHED_A2AV_AUX_HEAP
# include "dash_aux_heap.h"
#endif


#endif /* __DASH_SCHED_A2AV_AUX_H__ */
