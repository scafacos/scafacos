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


#ifndef __DASH_EXEC_H__
#define __DASH_EXEC_H__


#include <mpi.h>

#include "dash_conf.h"

#ifdef DS_RENAME
# include "dash_rename.h"
#endif

#include "z_pack.h"


#define DASH_EXEC_MAKE_SYM_LINEAR              0
#define DASH_EXEC_MAKE_SYM_LINEAR_RANDOM       1
#define DASH_EXEC_MAKE_SYM_HIERARCHIC          2
#define DASH_EXEC_MAKE_SYM_HIERARCHIC_RANDOM   3

#define DASH_EXEC_MAKE_SYM_DEFAULT  DASH_EXEC_MAKE_SYM_HIERARCHIC


struct _ds_t;
struct _ds_exec_t;

typedef dsint_t (*ds_exec_std_f)(struct _ds_exec_t *exec);
typedef void (*ds_exec_move_f)(struct _ds_exec_t *exec, dsint_t exec_id, dsint_t src_buf_id, dsint_t src_displs, dsint_t dst_buf_id, dsint_t dst_displs, dsint_t count);
#ifdef DASH_SYMMETRIC
typedef void (*ds_exec_sendrecv_replace_f)(struct _ds_exec_t *exec, int proc, dsint_t exec_id, dsint_t buf_id, dspint_t displ, dspint_t count);
# ifdef DASH_SYMMETRIC_AUX
typedef void (*ds_exec_sendrecv_aux_setup_f)(struct _ds_exec_t *exec, dsint_t max_nsyms);
typedef dsint_t (*ds_exec_sendrecv_aux_f)(struct _ds_exec_t *exec, int proc, dsint_t exec_id, dsint_t buf_id, dspint_t displ, dspint_t count, dsint_t aux_exec_id, dsint_t aux_buf_id, dspint_t aux_displ);
typedef void (*ds_exec_sendrecv_aux_intermediate_f)(struct _ds_exec_t *exec);
typedef void (*ds_exec_sendrecv_aux_finish_f)(struct _ds_exec_t *exec);
# endif
#endif

typedef struct _ds_exec_t
{
  struct _ds_t *ds;

  dsint_t naddrs;
  void *addrs[DASH_MAX_NBUFFERS];

  void *cxt;
  ds_exec_std_f pre_run, post_run, make;
  ds_exec_move_f move;
#ifdef DASH_SYMMETRIC
  dsint_t make_sym;
  ds_exec_sendrecv_replace_f sendrecv_replace;
# ifdef DASH_SYMMETRIC_AUX
  ds_exec_sendrecv_aux_setup_f sendrecv_aux_setup;
  ds_exec_sendrecv_aux_f sendrecv_aux;
  ds_exec_sendrecv_aux_intermediate_f sendrecv_aux_intermediate;
  ds_exec_sendrecv_aux_finish_f sendrecv_aux_finish;
# endif
#endif

} ds_exec_t, *_ds_exec_p;


dsint_t ds_exec_create(ds_exec_t *exec);
dsint_t ds_exec_destroy(ds_exec_t *exec);

dsint_t ds_exec_add_address(ds_exec_t *exec, void *addr);

dsint_t ds_exec_make(ds_exec_t *exec);


#endif /* __DASH_EXEC_H__ */
