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


#ifndef __DASH_EXEC_MPI_H__
#define __DASH_EXEC_MPI_H__


#include <mpi.h>

#include "zmpi_local.h"

#include "dash_core.h"


typedef struct _ds_exec_mpi_t
{
  dsint_t ntypes;
  MPI_Datatype mpi_types[DS_MAX_NBUFFERS];
  zmpil_t zmpil_types[DS_MAX_NBUFFERS];

  dsint_t nmax, n;
  MPI_Request *reqs;
  MPI_Status *stats;

  int *scounts, *sdispls, *rcounts, *rdispls;
  MPI_Datatype *stypes, *rtypes;
  MPI_Datatype *addr_types;

} ds_exec_mpi_t, *ds_exec_mpi_p;

#define DEFINE_EXEC_MPI(_s_, _v_)  ds_exec_mpi_t *_v_ = (_s_)->cxt

#define DS_EXEC_MPI_ISENDRECV_TAG         0
#define DS_EXEC_MPI_SENDRECV_REPLACE_TAG  0


dsint_t ds_exec_mpi_create(ds_exec_t *exec);
dsint_t ds_exec_mpi_destroy(ds_exec_t *exec);

dsint_t ds_exec_mpi_pre_run(ds_exec_t *exec);
dsint_t ds_exec_mpi_post_run(ds_exec_t *exec);
dsint_t ds_exec_mpi_make_isendrecv(ds_exec_t *exec);
dsint_t ds_exec_mpi_make_alltoallw(ds_exec_t *exec);

dsint_t ds_exec_mpi_add_address(ds_exec_t *exec, void *addr);
void ds_exec_mpi_move(ds_exec_t *exec, dsint_t exec_id, dsint_t src_buf_id, dsint_t src_displs, dsint_t dst_buf_id, dsint_t dst_displs, dsint_t count);
void ds_exec_mpi_sendrecv_replace(ds_exec_t *exec, int proc, dsint_t exec_id, dsint_t buf_id, dspint_t displ, dspint_t count);

dsint_t ds_exec_mpi_add_type(ds_exec_t *exec, MPI_Datatype type);
dsint_t ds_exec_mpi_sizefor(ds_exec_t *exec, dsint_t exec_id, dsint_t size);
dsint_t ds_exec_mpi_extent(ds_exec_t *exec, dsint_t exec_id);


#endif /* __DASH_EXEC_MPI_H__ */
