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
#include "dash_sched_a2a_sym.h"
#include "dash_sched_a2av_sym.h"
#include "dash_exec_mpi.h"
#include "zmpi_ataip.h"


#ifndef AD_TRACE_IF
# define AD_TRACE_IF  (z_mpi_rank == -1)
#endif


/* zmpi_var ZMPI_Alltoallv_inplace_sym_type */
int ZMPI_Alltoallv_inplace_sym_type = ZMPI_ALLTOALLV_INPLACE_SYM_TYPE_DEFAULT;

/* zmpi_var ZMPI_Alltoallv_inplace_aux_type ZMPI_Alltoallv_inplace_aux ZMPI_Alltoallv_inplace_aux_size ZMPI_Alltoallv_inplace_aux_blocks */
int ZMPI_Alltoallv_inplace_aux_type = ZMPI_ALLTOALLV_INPLACE_AUX_TYPE_HEAP;
void *ZMPI_Alltoallv_inplace_aux = NULL;
int ZMPI_Alltoallv_inplace_aux_size = 1000000;  /* default: 1 MB */
int ZMPI_Alltoallv_inplace_aux_static_blocks = -8;

/* zmpi_var ZMPI_Alltoallv_inplace_sync_on_init ZMPI_Alltoallv_inplace_sync_on_exit ZMPI_Alltoallv_inplace_sync_run */
int ZMPI_Alltoallv_inplace_sync_on_init = 0;
int ZMPI_Alltoallv_inplace_sync_on_exit = 0;
int ZMPI_Alltoallv_inplace_sync_run = 0;

/* zmpi_var ZMPI_Alltoall_inplace_symmetric_sym_type */
int ZMPI_Alltoall_inplace_symmetric_sym_type = ZMPI_ALLTOALL_INPLACE_SYMMETRIC_SYM_TYPE_DEFAULT;

/* zmpi_var ZMPI_Alltoallv_inplace_symmetric_sym_type */
int ZMPI_Alltoallv_inplace_symmetric_sym_type = ZMPI_ALLTOALLV_INPLACE_SYMMETRIC_SYM_TYPE_DEFAULT;

/* zmpi_var ZMPI_Alltoallv_inplace_symmetric_aux ZMPI_Alltoallv_inplace_symmetric_aux_size */
void *ZMPI_Alltoallv_inplace_symmetric_aux = NULL;
int ZMPI_Alltoallv_inplace_symmetric_aux_size = 0;  /* default: 0 */

/* zmpi_var ZMPI_Alltoallv_inplace_times */
double ZMPI_Alltoallv_inplace_times[10];


static int create_proclists_from_mpi_comm_topology(MPI_Comm comm, int rank, int *nsend_procs, int **send_procs, int *nrecv_procs, int **recv_procs)
{
  int status;

  int n, i;


  MPI_Topo_test(comm, &status);

  switch (status)
  {
    case MPI_CART:
      MPI_Cartdim_get(comm, &n);
      *nsend_procs = *nrecv_procs = 2 * n;
      *send_procs = *recv_procs = z_alloc(2 * n, sizeof(int));
      for (i = 0; i < n; ++i) MPI_Cart_shift(comm, i, -1, &(*send_procs)[2 * i + 1], &(*send_procs)[2 * i + 0]);
      break;
    case MPI_GRAPH:
      MPI_Graph_neighbors_count(comm, rank, &n);
      *nsend_procs = *nrecv_procs = n;
      *send_procs = *recv_procs = z_alloc(n, sizeof(int));
      MPI_Graph_neighbors(comm, rank, n, *send_procs);
      break;
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 2
    case MPI_DIST_GRAPH:
      MPI_Dist_graph_neighbors_count(comm, nrecv_procs, nsend_procs, &n);
      *send_procs = z_alloc(*nsend_procs + *nrecv_procs, sizeof(int));
      *recv_procs = *send_procs + *nsend_procs;
      MPI_Dist_graph_neighbors(comm, *nrecv_procs, *recv_procs, MPI_UNWEIGHTED, *nsend_procs, *send_procs, MPI_UNWEIGHTED);
      break;
#endif
    default:
      return 0;
  }

  return 1;
}


static void destroy_proclists(int **send_procs, int **recv_procs)
{
  z_free(*send_procs);
}


static int _ZMPI_Alltoallv_inplace(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, int nsends, int *sranks, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, int nrecvs, int *rranks, MPI_Comm comm)
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

#ifdef DASH_SYMMETRIC
  switch (ZMPI_Alltoallv_inplace_sym_type)
  {
    case ZMPI_SYM_TYPE_LINEAR:
      exec.make_sym = DASH_EXEC_MAKE_SYM_LINEAR;
      break;
    case ZMPI_SYM_TYPE_LINEAR_RANDOM:
      exec.make_sym = DASH_EXEC_MAKE_SYM_LINEAR_RANDOM;
      break;
    case ZMPI_SYM_TYPE_HIERARCHIC:
      exec.make_sym = DASH_EXEC_MAKE_SYM_HIERARCHIC;
      break;
    case ZMPI_SYM_TYPE_HIERARCHIC_RANDOM:
      exec.make_sym = DASH_EXEC_MAKE_SYM_HIERARCHIC_RANDOM;
      break;
  }
#endif

  stype_id = ds_exec_mpi_add_type(&exec, stype);
  rtype_id = ds_exec_mpi_add_type(&exec, rtype);

  if (ZMPI_Alltoallv_inplace_aux_static_blocks < 0) aux_blocks = comm_size / -ZMPI_Alltoallv_inplace_aux_static_blocks;
  else aux_blocks = ZMPI_Alltoallv_inplace_aux_static_blocks;
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

  ds_sched_a2av_skip_sym(&sched, (ZMPI_Alltoallv_inplace_sym_type == ZMPI_SYM_TYPE_SKIP)?1:0);

  rbuf_buf_id = ds_sched_add_buffer(&sched, rbuf_addr_id);
  ds_sched_set_send(&sched, rbuf_buf_id, (nsends < 0)?comm_size:nsends, scounts, sdispls, sranks, stype_id);
  ds_sched_set_recv(&sched, rbuf_buf_id, (nrecvs < 0)?comm_size:nrecvs, rcounts, rdispls, rranks, rtype_id);

  aux_buf_id = ds_sched_add_buffer(&sched, aux_addr_id);
  ds_sched_set_aux(&sched, aux_buf_id, aux_buf_size);

  switch (ZMPI_Alltoallv_inplace_aux_type)
  {
    case ZMPI_ALLTOALLV_INPLACE_AUX_TYPE_STATIC:
      ds_aux_static_create(&sched_a2av_aux, aux_buf_size, aux_blocks);
      break;
    case ZMPI_ALLTOALLV_INPLACE_AUX_TYPE_HEAP:
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
    case ZMPI_ALLTOALLV_INPLACE_AUX_TYPE_STATIC:
      ds_aux_static_destroy(&sched_a2av_aux);
      break;
    case ZMPI_ALLTOALLV_INPLACE_AUX_TYPE_HEAP:
      ds_aux_heap_destroy(&sched_a2av_aux);
      break;
  }

  ds_sched_a2av_destroy(&sched);

  if (ZMPI_Alltoallv_inplace_aux == NULL) z_free(aux);

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
  DS_TIMING_CMD(ZMPI_Alltoallv_inplace_times[9] = ds_times[DS_TIMES_RUN_WHILE_ROUNDS];);

  DS_TIMING_CMD(
    if (comm_rank == 0)
      printf("%d: _ZMPI_Alltoallv_inplace: %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", comm_rank,
        ZMPI_Alltoallv_inplace_times[0], ZMPI_Alltoallv_inplace_times[1], ZMPI_Alltoallv_inplace_times[2], ZMPI_Alltoallv_inplace_times[3], ZMPI_Alltoallv_inplace_times[4],
        ZMPI_Alltoallv_inplace_times[5], ZMPI_Alltoallv_inplace_times[6], ZMPI_Alltoallv_inplace_times[7], ZMPI_Alltoallv_inplace_times[8], ZMPI_Alltoallv_inplace_times[9]);
  );

  return exit_code;
}


int ZMPI_Alltoallv_inplace(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallv_inplace */
{
  return _ZMPI_Alltoallv_inplace(sbuf, scounts, sdispls, stype, -1, NULL, rbuf, rcounts, rdispls, rtype, -1, NULL, comm);
}


int ZMPI_Neighbor_alltoallv_inplace(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm) /* zmpi_func ZMPI_Neighbor_alltoallv_inplace */
{
  int exit_code;
  int comm_size, comm_rank;

  int nsend_procs, *send_procs, nrecv_procs, *recv_procs;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if (!create_proclists_from_mpi_comm_topology(comm, comm_rank, &nsend_procs, &send_procs, &nrecv_procs, &recv_procs))
  {
#ifdef ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE
    fprintf(ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE, "%d: ZMPI_Neighbor_alltoallv_inplace: error: communicator has no (supported) topology!\n", comm_rank);
#endif
    exit_code = 1;
    goto exit;
  }

  exit_code = _ZMPI_Alltoallv_inplace(sbuf, scounts, sdispls, stype, nsend_procs, send_procs, rbuf, rcounts, rdispls, rtype, nrecv_procs, recv_procs, comm);

  destroy_proclists(&send_procs, &recv_procs);

exit:

  return exit_code;
}


#ifdef DASH_SYMMETRIC


static int _ZMPI_Alltoall_inplace_symmetric(void *sbuf, int scount, MPI_Datatype stype, int nsends, int *sranks, void *rbuf, int rcount, MPI_Datatype rtype, int nrecvs, int *rranks, MPI_Comm comm){
  int exit_code = MPI_SUCCESS;

  ds_t ds;
  ds_sched_t sched;
  ds_exec_t exec;

  dsint_t rtype_id;
  dsint_t rbuf_addr_id, rbuf_buf_id;

#ifdef DASH_SYMMETRIC_AUX
  void *aux;
  dsint_t aux_size, aux_buf_size;

  dsint_t aux_addr_id, aux_buf_id;
#endif

  int comm_size, comm_rank;


  if (ZMPI_Alltoallv_inplace_sync_on_init) MPI_Barrier(comm);

  if (rbuf == NULL) return 0;

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  ds_exec_mpi_create(&exec);

  switch (ZMPI_Alltoall_inplace_symmetric_sym_type)
  {
    case ZMPI_SYM_TYPE_LINEAR:
      exec.make_sym = DASH_EXEC_MAKE_SYM_LINEAR;
      break;
    case ZMPI_SYM_TYPE_LINEAR_RANDOM:
      exec.make_sym = DASH_EXEC_MAKE_SYM_LINEAR_RANDOM;
      break;
    case ZMPI_SYM_TYPE_HIERARCHIC:
      exec.make_sym = DASH_EXEC_MAKE_SYM_HIERARCHIC;
      break;
    case ZMPI_SYM_TYPE_HIERARCHIC_RANDOM:
      exec.make_sym = DASH_EXEC_MAKE_SYM_HIERARCHIC_RANDOM;
      break;
  }

  rtype_id = ds_exec_mpi_add_type(&exec, rtype);

#ifdef DASH_SYMMETRIC_AUX
  aux_size = ZMPI_Alltoallv_inplace_symmetric_aux_size;

  if (ZMPI_Alltoallv_inplace_symmetric_aux == NULL && aux_size > 0) aux = z_alloc(aux_size, 1);
  else aux = ZMPI_Alltoallv_inplace_symmetric_aux;

  aux_buf_size = ds_exec_mpi_sizefor(&exec, rtype_id, aux_size);
#endif

  rbuf_addr_id = ds_exec_mpi_add_address(&exec, rbuf);
#ifdef DASH_SYMMETRIC_AUX
  aux_addr_id = ds_exec_mpi_add_address(&exec, aux);
#endif

  ds_sched_a2a_sym_create(&sched);

  rbuf_buf_id = ds_sched_add_buffer(&sched, rbuf_addr_id);
  ds_sched_set_send(&sched, rbuf_buf_id, (nsends < 0)?comm_size:nsends, &scount, NULL, sranks, rtype_id);
  ds_sched_set_recv(&sched, rbuf_buf_id, (nrecvs < 0)?comm_size:nrecvs, &rcount, NULL, rranks, rtype_id);

#ifdef DASH_SYMMETRIC_AUX
  aux_buf_id = ds_sched_add_buffer(&sched, aux_addr_id);
  ds_sched_set_aux(&sched, aux_buf_id, aux_buf_size);
#endif

  ds_create(&ds, &sched, &exec);

  ds.comm = comm;
  ds.comm_size = comm_size;
  ds.comm_rank = comm_rank;

  ds_core_run_sync = ZMPI_Alltoallv_inplace_sync_run;

  ds_run(&ds);

  ds_destroy(&ds);

  ds_sched_a2a_sym_destroy(&sched);

#ifdef DASH_SYMMETRIC_AUX
  if (ZMPI_Alltoallv_inplace_symmetric_aux == NULL && aux_size > 0) z_free(aux);
#endif

  ds_exec_mpi_destroy(&exec);

  if (ZMPI_Alltoallv_inplace_sync_on_exit) MPI_Barrier(comm);

  DS_TIMING_CMD(
    if (comm_rank == 0)
      printf("%d: _ZMPI_Alltoall_inplace_symmetric: %f  %f  %f  %f  %f  %f  %f  %f  %f\n", comm_rank,
        ds_times[DS_TIMES_RUN], ds_times[DS_TIMES_RUN_PRE], ds_times[DS_TIMES_RUN_WHILE], ds_times[DS_TIMES_RUN_WHILE_PRE], ds_times[DS_TIMES_RUN_WHILE_MAKE],
        ds_times[DS_TIMES_RUN_WHILE_POST], ds_times[DS_TIMES_RUN_WHILE_FINISHED], ds_times[DS_TIMES_RUN_WHILE_ROUNDS], ds_times[DS_TIMES_RUN_POST]);
  );

  return exit_code;
}


int ZMPI_Alltoall_inplace_symmetric(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_inplace_symmetric */
{
  return _ZMPI_Alltoall_inplace_symmetric(sbuf, scount, stype, -1, NULL, rbuf, rcount, rtype, -1, NULL, comm);
}


int ZMPI_Neighbor_alltoall_inplace_symmetric(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, MPI_Comm comm) /* zmpi_func ZMPI_Neighbor_alltoall_inplace_symmetric */
{
  int exit_code;
  int comm_size, comm_rank;

  int nsend_procs, *send_procs, nrecv_procs, *recv_procs;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if (!create_proclists_from_mpi_comm_topology(comm, comm_rank, &nsend_procs, &send_procs, &nrecv_procs, &recv_procs))
  {
#ifdef ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE
    fprintf(ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE, "%d: ZMPI_Neighbor_alltoall_inplace_symmetric: error: communicator has no (supported) topology!\n", comm_rank);
#endif
    exit_code = 1;
    goto exit;
  }

  exit_code = _ZMPI_Alltoall_inplace_symmetric(sbuf, scount, stype, nsend_procs, send_procs, rbuf, rcount, rtype, nrecv_procs, recv_procs, comm);

  destroy_proclists(&send_procs, &recv_procs);

exit:

  return exit_code;
}


static int _ZMPI_Alltoallv_inplace_symmetric(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, int nsends, int *sranks, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, int nrecvs, int *rranks, MPI_Comm comm)
{
  int exit_code = MPI_SUCCESS;

  ds_t ds;
  ds_sched_t sched;
  ds_exec_t exec;

  dsint_t rtype_id;
  dsint_t rbuf_addr_id, rbuf_buf_id;

#ifdef DASH_SYMMETRIC_AUX
  void *aux;
  dsint_t aux_size, aux_buf_size;

  dsint_t aux_addr_id, aux_buf_id;
#endif

  int comm_size, comm_rank;


  if (ZMPI_Alltoallv_inplace_sync_on_init) MPI_Barrier(comm);

  if (rbuf == NULL) return 0;

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  ds_exec_mpi_create(&exec);

  switch (ZMPI_Alltoallv_inplace_symmetric_sym_type)
  {
    case ZMPI_SYM_TYPE_LINEAR:
      exec.make_sym = DASH_EXEC_MAKE_SYM_LINEAR;
      break;
    case ZMPI_SYM_TYPE_LINEAR_RANDOM:
      exec.make_sym = DASH_EXEC_MAKE_SYM_LINEAR_RANDOM;
      break;
    case ZMPI_SYM_TYPE_HIERARCHIC:
      exec.make_sym = DASH_EXEC_MAKE_SYM_HIERARCHIC;
      break;
    case ZMPI_SYM_TYPE_HIERARCHIC_RANDOM:
      exec.make_sym = DASH_EXEC_MAKE_SYM_HIERARCHIC_RANDOM;
      break;
  }

  rtype_id = ds_exec_mpi_add_type(&exec, rtype);

#ifdef DASH_SYMMETRIC_AUX
  aux_size = ZMPI_Alltoallv_inplace_symmetric_aux_size;

  if (ZMPI_Alltoallv_inplace_symmetric_aux == NULL && aux_size > 0) aux = z_alloc(aux_size, 1);
  else aux = ZMPI_Alltoallv_inplace_symmetric_aux;

  aux_buf_size = ds_exec_mpi_sizefor(&exec, rtype_id, aux_size);
#endif

  rbuf_addr_id = ds_exec_mpi_add_address(&exec, rbuf);
#ifdef DASH_SYMMETRIC_AUX
  aux_addr_id = ds_exec_mpi_add_address(&exec, aux);
#endif

  ds_sched_a2av_sym_create(&sched);

  rbuf_buf_id = ds_sched_add_buffer(&sched, rbuf_addr_id);
  ds_sched_set_send(&sched, rbuf_buf_id, (nsends < 0)?comm_size:nsends, scounts, sdispls, sranks, rtype_id);
  ds_sched_set_recv(&sched, rbuf_buf_id, (nrecvs < 0)?comm_size:nrecvs, rcounts, rdispls, rranks, rtype_id);

#ifdef DASH_SYMMETRIC_AUX
  aux_buf_id = ds_sched_add_buffer(&sched, aux_addr_id);
  ds_sched_set_aux(&sched, aux_buf_id, aux_buf_size);
#endif

  ds_create(&ds, &sched, &exec);

  ds.comm = comm;
  ds.comm_size = comm_size;
  ds.comm_rank = comm_rank;

  ds_core_run_sync = ZMPI_Alltoallv_inplace_sync_run;

  ds_run(&ds);

  ds_destroy(&ds);

  ds_sched_a2av_sym_destroy(&sched);

#ifdef DASH_SYMMETRIC_AUX
  if (ZMPI_Alltoallv_inplace_symmetric_aux == NULL && aux_size > 0) z_free(aux);
#endif

  ds_exec_mpi_destroy(&exec);

  if (ZMPI_Alltoallv_inplace_sync_on_exit) MPI_Barrier(comm);

  return exit_code;
}


int ZMPI_Alltoallv_inplace_symmetric(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallv_inplace_symmetric */
{
  return _ZMPI_Alltoallv_inplace_symmetric(sbuf, scounts, sdispls, stype, -1, NULL, rbuf, rcounts, rdispls, rtype, -1, NULL, comm);
}


int ZMPI_Neighbor_alltoallv_inplace_symmetric(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm) /* zmpi_func ZMPI_Neighbor_alltoallv_inplace_symmetric */
{
  int exit_code;
  int comm_size, comm_rank;

  int nsend_procs, *send_procs, nrecv_procs, *recv_procs;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if (!create_proclists_from_mpi_comm_topology(comm, comm_rank, &nsend_procs, &send_procs, &nrecv_procs, &recv_procs))
  {
#ifdef ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE
    fprintf(ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE, "%d: ZMPI_Neighbor_alltoallv_inplace_symmetric: error: communicator has no (supported) topology!\n", comm_rank);
#endif
    exit_code = 1;
    goto exit;
  }

  exit_code = _ZMPI_Alltoallv_inplace_symmetric(sbuf, scounts, sdispls, stype, nsend_procs, send_procs, rbuf, rcounts, rdispls, rtype, nrecv_procs, recv_procs, comm);

  destroy_proclists(&send_procs, &recv_procs);

exit:

  return exit_code;
}


#endif /* DASH_SYMMETRIC */


#undef AD_TRACE_IF
