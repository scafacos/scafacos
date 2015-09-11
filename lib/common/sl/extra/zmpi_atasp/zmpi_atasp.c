/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015 Michael Hofmann
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


#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "spec_core.h"
#include "zmpi_atasp.h"


#define ZMPI_FAILURE  1


/* zmpi_var ZMPI_Alltoall_specific_type */
int ZMPI_Alltoall_specific_type = ZMPI_ALLTOALL_SPECIFIC_TYPE_DEFAULT;

/* zmpi_var ZMPI_Neighbor_alltoall_specific_type */
int ZMPI_Neighbor_alltoall_specific_type = ZMPI_NEIGHBOR_ALLTOALL_SPECIFIC_TYPE_DEFAULT;


int ZMPI_Tproc_create_tproc(ZMPI_Tproc *tproc, ZMPI_TPROC_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef) /* zmpi_func ZMPI_Tproc_create_tproc */
{
  if (exdef != ZMPI_TPROC_EXDEF_NULL && exdef->type != 1) return ZMPI_FAILURE;

  spec_tproc_create(tproc, tfn, NULL, NULL, NULL, 0);

  spec_tproc_set_reset(*tproc, rfn);

  if (exdef != ZMPI_TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tproc(*tproc, &exdef->tproc_ext);

  return MPI_SUCCESS;
}


int ZMPI_Tproc_create_tproc_mod(ZMPI_Tproc *tproc, ZMPI_TPROC_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef) /* zmpi_func ZMPI_Tproc_create_tproc_mod */
{
  if (exdef != ZMPI_TPROC_EXDEF_NULL && exdef->type != 2) return ZMPI_FAILURE;

  spec_tproc_create(tproc, NULL, tfn, NULL, NULL, 0);

  spec_tproc_set_reset(*tproc, rfn);

  if (exdef != ZMPI_TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tproc_mod(*tproc, &exdef->tproc_mod_ext);

  return MPI_SUCCESS;
}


int ZMPI_Tproc_create_tprocs(ZMPI_Tproc *tproc, int max_tprocs, ZMPI_TPROCS_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef) /* zmpi_func ZMPI_Tproc_create_tprocs */
{
  if (exdef != ZMPI_TPROC_EXDEF_NULL && exdef->type != 3) return ZMPI_FAILURE;

  spec_tproc_create(tproc, NULL, NULL, tfn, NULL, max_tprocs);

  spec_tproc_set_reset(*tproc, rfn);

  if (exdef != ZMPI_TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tprocs(*tproc, &exdef->tprocs_ext);

  return MPI_SUCCESS;
}


int ZMPI_Tproc_create_tprocs_mod(ZMPI_Tproc *tproc, int max_tprocs, ZMPI_TPROCS_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, const ZMPI_Tproc_exdef exdef) /* zmpi_func ZMPI_Tproc_create_tprocs_mod */
{
  if (exdef != ZMPI_TPROC_EXDEF_NULL && exdef->type != 4) return ZMPI_FAILURE;

  spec_tproc_create(tproc, NULL, NULL, NULL, tfn, max_tprocs);

  spec_tproc_set_reset(*tproc, rfn);

  if (exdef != ZMPI_TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tprocs_mod(*tproc, &exdef->tprocs_mod_ext);

  return MPI_SUCCESS;
}


int ZMPI_Tproc_free(ZMPI_Tproc *tproc) /* zmpi_func ZMPI_Tproc_free */
{
  spec_tproc_destroy(tproc);

  return MPI_SUCCESS;
}


int ZMPI_Tproc_set_neighbors(ZMPI_Tproc tproc, int nneighbors, int *neighbors, MPI_Comm comm) /* zmpi_func ZMPI_Tproc_set_neighbors */
{
  int comm_size, comm_rank;

#ifdef SPEC_PROCLISTS
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  spec_tproc_set_proclists(tproc, nneighbors, neighbors, -1, NULL, comm_size, comm_rank, comm);

  return MPI_SUCCESS;
#else
  return 1;
#endif
}


int ZMPI_Tproc_set_proclists(ZMPI_Tproc tproc, int ndstprocs, int *dstprocs, int nsrcprocs, int *srcprocs, MPI_Comm comm) /* zmpi_func ZMPI_Tproc_set_proclists */
{
  int comm_size, comm_rank;

#ifdef SPEC_PROCLISTS
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  spec_tproc_set_proclists(tproc, ndstprocs, dstprocs, nsrcprocs, srcprocs, comm_size, comm_rank, comm);

  return MPI_SUCCESS;
#else
  return 1;
#endif
}


#if MPI_VERSION < 3

int ZMPI_Get_elements(const ZMPI_Status *status, MPI_Datatype datatype, int *count)
{
  if (status == ZMPI_STATUS_IGNORE) return MPI_ERR_ARG;

  *count = *status;

  return MPI_SUCCESS;
}

#endif


static int _ZMPI_Alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, int comm_size, int comm_rank, MPI_Comm comm,
#if MPI_VERSION >= 3
  MPI_Status
#else
  ZMPI_Status
#endif
  *status, int type, const char *name)
{
  int exit_code, received = 0;

  spec_elem_t sb, rb, *b, *xb;


  if (sbuf == MPI_IN_PLACE && rbuf == MPI_IN_PLACE)
  {
#ifdef ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE
    fprintf(ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE, "%d: %s: error: either sbuf or rbuf can be MPI_IN_PLACE!\n", comm_rank, name);
#endif
    exit_code = ZMPI_FAILURE;
    goto exit;
  }

  if (tproc == ZMPI_TPROC_NULL)
  {
    exit_code = MPI_SUCCESS;
    goto exit;
  }

  if (sbuf != MPI_IN_PLACE)
  {
    sb.buf = sbuf;
    sb.count = sb.max_count = scount;
    sb.mpi_type = stype;
    zmpil_create(&sb.zmpil_type, stype, 1);

    b = &sb;
  }

  if (rbuf != MPI_IN_PLACE)
  {
    rb.buf = rbuf;
    rb.count = rb.max_count = rcount;
    rb.mpi_type = rtype;
    zmpil_create(&rb.zmpil_type, rtype, 1);

    b = &rb;
  }

  xb = NULL;

  if (sbuf == MPI_IN_PLACE || rbuf == MPI_IN_PLACE)
  {
    b->count = scount;
    b->max_count = rcount;

    switch (type)
    {
#ifdef SPEC_ALLTOALLV
      case ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLV:
        exit_code = spec_alltoallv_ip(b, xb, tproc, tproc_data, comm_size, comm_rank, comm);
        break;
#endif
      default:
#ifdef ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE
        fprintf(ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE, "%d: %s: error: in-place implementation of type %d is not available!\n", comm_rank, name);
#endif
        exit_code = ZMPI_FAILURE;
        break;
    }
    received = b->count;

  } else {

    switch (type)
    {
#ifdef SPEC_ALLTOALLV
      case ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLV:
        exit_code = spec_alltoallv_db(&sb, &rb, xb, tproc, tproc_data, comm_size, comm_rank, comm);
        break;
#endif
#ifdef SPEC_ALLTOALLW
      case ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLW:
        exit_code = spec_alltoallw_db(&sb, &rb, xb, tproc, tproc_data, comm_size, comm_rank, comm);
        break;
#endif
#ifdef SPEC_PUT
      case ZMPI_ALLTOALL_SPECIFIC_TYPE_PUT:
        exit_code = spec_put_db(&sb, &rb, xb, tproc, tproc_data, comm_size, comm_rank, comm);
        break;
      case ZMPI_ALLTOALL_SPECIFIC_TYPE_PUT_2PHASES:
        exit_code = spec_put_2phases_db(&sb, &rb, xb, tproc, tproc_data, comm_size, comm_rank, comm);
        break;
#endif
#ifdef SPEC_SENDRECV
      case ZMPI_ALLTOALL_SPECIFIC_TYPE_SENDRECV:
        exit_code = spec_sendrecv_db(&sb, &rb, xb, tproc, tproc_data, comm_size, comm_rank, comm);
        break;
#endif
      default:
#ifdef ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE
        fprintf(ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE, "%d: %s: error: implementation of type %d is not available!\n", comm_rank, name);
#endif
        exit_code = ZMPI_FAILURE;
        break;
    }
    received = rb.count;
  }

  zmpil_destroy(&sb.zmpil_type);

  zmpil_destroy(&rb.zmpil_type);

exit:

#if MPI_VERSION >= 3
  if (status != MPI_STATUS_IGNORE) MPI_Status_set_elements(status, rtype, received);
#else
  if (status != ZMPI_STATUS_IGNORE) *status = received;
#endif

  return exit_code;
}


int ZMPI_Alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, MPI_Comm comm,
#if MPI_VERSION >= 3
  MPI_Status
#else
  ZMPI_Status
#endif
  *status) /* zmpi_func ZMPI_Alltoall_specific */
{
  int comm_size, comm_rank;

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  return _ZMPI_Alltoall_specific(sbuf, scount, stype, rbuf, rcount, rtype, tproc, tproc_data, comm_size, comm_rank, comm, status, ZMPI_Alltoall_specific_type, "ZMPI_Alltoall_specific");
}


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


int ZMPI_Neighbor_alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, MPI_Comm comm,
#if MPI_VERSION >= 3
  MPI_Status
#else
  ZMPI_Status
#endif
  *status) /* zmpi_func ZMPI_Neighbor_alltoall_specific */
{
  int exit_code;
  int comm_size, comm_rank;

  int nsend_procs, *send_procs, nrecv_procs, *recv_procs;

  ZMPI_Tproc newtproc;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if (!create_proclists_from_mpi_comm_topology(comm, comm_rank, &nsend_procs, &send_procs, &nrecv_procs, &recv_procs))
  {
#ifdef ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE
    fprintf(ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE, "%d: ZMPI_Neighbor_alltoall_specific: error: communicator has no (supported) topology!\n", comm_rank);
#endif
    exit_code = ZMPI_FAILURE;
    goto exit;
  }

  spec_tproc_duplicate(tproc, &newtproc);

  ZMPI_Tproc_set_proclists(newtproc, nsend_procs, send_procs, nrecv_procs, recv_procs, comm);

  exit_code = _ZMPI_Alltoall_specific(sbuf, scount, stype, rbuf, rcount, rtype, newtproc, tproc_data, comm_size, comm_rank, comm, status, ZMPI_Neighbor_alltoall_specific_type, "ZMPI_Neighbor_alltoall_specific");

  destroy_proclists(&send_procs, &recv_procs);

  spec_tproc_destroy(&newtproc);

exit:

  return exit_code;
}
