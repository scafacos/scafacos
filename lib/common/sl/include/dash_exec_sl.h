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


#ifndef __DS_EXEC_SL_H__
#define __DS_EXEC_SL_H__

#include "sl_common.h"

#include "dash_core.h"


typedef struct _ds_exec_sl_t
{
  dsint_t nmax, n;
  MPI_Request *reqs;
  MPI_Status *stats;

} ds_exec_sl_t, *ds_exec_sl_p;

#define DEFINE_EXEC_SL(_s_, _v_)  ds_exec_sl_t *_v_ = (_s_)->cxt

#define DS_EXEC_SL_DATA_TAG  0


dsint_t ds_exec_sl_create(ds_exec_t *exec);
dsint_t ds_exec_sl_destroy(ds_exec_t *exec);

dsint_t ds_exec_sl_pre_run(ds_exec_t *exec);
dsint_t ds_exec_sl_post_run(ds_exec_t *exec);
dsint_t ds_exec_sl_make(ds_exec_t *exec);

dsint_t ds_exec_sl_add_address(ds_exec_t *exec, elements_t *s);

void ds_exec_sl_move(ds_exec_t *exec, dsint_t exec_id, dsint_t src_buf_id, dsint_t src_displs, dsint_t dst_buf_id, dsint_t dst_displs, dsint_t count);


#endif /* __DS_EXEC_SL_H__ */
