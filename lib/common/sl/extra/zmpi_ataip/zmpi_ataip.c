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


#include <mpi.h>

#include "dash_core.h"
#include "dash_sched_a2av.h"
#include "dash_sched_a2av_aux.h"
#include "dash_exec.h"
#include "dash_exec_mpi.h"
#include "zmpi_ataip.h"


#define AD_TRACE_IF  (z_mpi_rank == -1)


/* zmpi_var ZMPI_Alltoallv_inplace_aux_type */
int ZMPI_Alltoallv_inplace_aux_type = DS_AUX_TYPE_HEAP;

/* zmpi_var ZMPI_Alltoallv_inplace_aux ZMPI_Alltoallv_inplace_aux_size ZMPI_Alltoallv_inplace_aux_blocks */
void *ZMPI_Alltoallv_inplace_aux = NULL;
int ZMPI_Alltoallv_inplace_aux_size = 1000000;  /* default: 1 MB */
int ZMPI_Alltoallv_inplace_aux_blocks = -8;

/* zmpi_var ZMPI_Alltoallv_inplace_sync_on_init ZMPI_Alltoallv_inplace_sync_on_exit ZMPI_Alltoallv_inplace_sync_run */
int ZMPI_Alltoallv_inplace_sync_on_init = 0;
int ZMPI_Alltoallv_inplace_sync_on_exit = 0;
int ZMPI_Alltoallv_inplace_sync_run = 0;

/* zmpi_var ZMPI_Alltoallv_inplace_times */
double ZMPI_Alltoallv_inplace_times[9];


int ZMPI_Alltoallv_inplace(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallv_inplace */
{
  int exit_code = MPI_SUCCESS;

  ds_t ds;
  ds_sched_t sched;
  ds_sched_a2av_aux_t sched_a2av_aux;
  ds_exec_t exec;

  void *aux;
  dsint_t aux_size, aux_blocks, aux_buf_size;

  dsint_t stype_id, rtype_id;
  dsint_t rbuf_addr_id, aux_addr_id;
  dsint_t rbuf_buf_id, aux_buf_id;

  int comm_size, comm_rank;


  if (ZMPI_Alltoallv_inplace_sync_on_init) MPI_Barrier(comm);

  if (rbuf == NULL) return 0;

  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[0] = z_time_get_s(););
  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[1] = z_time_get_s(););

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  ds_exec_mpi_create(&exec);

  stype_id = ds_exec_mpi_add_type(&exec, stype);
  rtype_id = ds_exec_mpi_add_type(&exec, rtype);

  if (ZMPI_Alltoallv_inplace_aux_blocks < 0) aux_blocks = comm_size / -ZMPI_Alltoallv_inplace_aux_blocks;
  else aux_blocks = ZMPI_Alltoallv_inplace_aux_blocks;
  aux_blocks = z_minmax(1, aux_blocks, comm_size);

  aux_size = ZMPI_Alltoallv_inplace_aux_size;

  if (ds_exec_mpi_sizefor(&exec, rtype_id, aux_size) < aux_blocks)
  {
    aux_size = ds_exec_mpi_extent(&exec, rtype_id) * aux_blocks;
    ZMPI_Alltoallv_inplace_aux = NULL;
  }

  if (ZMPI_Alltoallv_inplace_aux == NULL) aux = z_alloc(aux_size, 1);
  else aux = ZMPI_Alltoallv_inplace_aux;

  if (aux == NULL)
  {
    exit_code = 1;
    goto exit_on_aux;
  }

  aux_buf_size = ds_exec_mpi_sizefor(&exec, rtype_id, aux_size);

  rbuf_addr_id = ds_exec_mpi_add_address(&exec, rbuf);
  aux_addr_id = ds_exec_mpi_add_address(&exec, aux);

  ds_sched_a2av_create(&sched);

  rbuf_buf_id = ds_sched_add_buffer(&sched, rbuf_addr_id);
  ds_sched_set_send(&sched, rbuf_buf_id, comm_size, scounts, sdispls, NULL, stype_id);
  ds_sched_set_recv(&sched, rbuf_buf_id, comm_size, rcounts, rdispls, NULL, rtype_id);

  aux_buf_id = ds_sched_add_buffer(&sched, aux_addr_id);
  ds_sched_set_aux(&sched, aux_buf_id, aux_buf_size);

  Z_TRACE_IF(AD_TRACE_IF, "aux type = %d", ZMPI_Alltoallv_inplace_aux_type);

  switch (ZMPI_Alltoallv_inplace_aux_type)
  {
    case DS_AUX_TYPE_STATIC:
      ds_aux_static_create(&sched_a2av_aux, aux_buf_size, aux_blocks);
      break;
    case DS_AUX_TYPE_HEAP:
      ds_aux_heap_create(&sched_a2av_aux, aux_buf_size, comm_size);
      break;
  }

  ds_sched_a2av_set_aux(&sched, &sched_a2av_aux);

  ds_create(&ds, &sched, &exec);
  
  ds.comm = comm;
  ds.comm_size = comm_size;
  ds.comm_rank = comm_rank;

  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[1] = z_time_get_s() - ZMPI_Alltoallv_inplace_times[1];);
  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[2] = z_time_get_s(););

  ds_core_run_sync = ZMPI_Alltoallv_inplace_sync_run;

  ds_run(&ds);

  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[2] = z_time_get_s() - ZMPI_Alltoallv_inplace_times[2];);
  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[3] = z_time_get_s(););

  ds_destroy(&ds);

  switch (ZMPI_Alltoallv_inplace_aux_type)
  {
    case DS_AUX_TYPE_STATIC:
      ds_aux_static_destroy(&sched_a2av_aux);
      break;
    case DS_AUX_TYPE_HEAP:
      ds_aux_heap_destroy(&sched_a2av_aux);
      break;
  }

  ds_sched_a2av_destroy(&sched);

  if (ZMPI_Alltoallv_inplace_aux == NULL) z_free(aux);
  else aux = NULL;

exit_on_aux:
  ds_exec_mpi_destroy(&exec);

  if (ZMPI_Alltoallv_inplace_sync_on_exit) MPI_Barrier(comm);

  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[3] = z_time_get_s() - ZMPI_Alltoallv_inplace_times[3];);
  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[0] = z_time_get_s() - ZMPI_Alltoallv_inplace_times[0];);

  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[4] = ds_times[DS_TIMES_RUN_WHILE];);
  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[5] = ds_times[DS_TIMES_RUN_WHILE_PRE];);
  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[6] = ds_times[DS_TIMES_RUN_WHILE_MAKE];);
  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[7] = ds_times[DS_TIMES_RUN_WHILE_POST];);
  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[8] = ds_times[DS_TIMES_RUN_WHILE_FINISHED];);

  return exit_code;
}


int ZMPI_Alltoallv_inplace_static(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallv_inplace_static */
{
  int ret, old_type = ZMPI_Alltoallv_inplace_aux_type;

  ZMPI_Alltoallv_inplace_aux_type = DS_AUX_TYPE_STATIC;

  ret = ZMPI_Alltoallv_inplace(sbuf, scounts, sdispls, stype, rbuf, rcounts, rdispls, rtype, comm);

  ZMPI_Alltoallv_inplace_aux_type = old_type;

  return ret;
}


int ZMPI_Alltoallv_inplace_heap(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallv_inplace_heap */
{
  int ret, old_type = ZMPI_Alltoallv_inplace_aux_type;

  ZMPI_Alltoallv_inplace_aux_type = DS_AUX_TYPE_HEAP;

  ret = ZMPI_Alltoallv_inplace(sbuf, scounts, sdispls, stype, rbuf, rcounts, rdispls, rtype, comm);

  ZMPI_Alltoallv_inplace_aux_type = old_type;

  return ret;
}
