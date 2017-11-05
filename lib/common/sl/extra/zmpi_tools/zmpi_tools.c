/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Michael Hofmann
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
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "z_pack.h"

#include "zmpi_tools.h"


/*#define WITH_TIMING*/

#ifdef WITH_TIMING
# define TSTART(_t_, _c_)  Z_MOP(MPI_Barrier(_c_); Z_TIMING_START(_t_);)
# define TSTOP(_t_, _c_)   Z_MOP(MPI_Barrier(_c_); Z_TIMING_STOP(_t_);)
#else
# define TSTART(_t_, _c_)  Z_NOP()
# define TSTOP(_t_, _c_)   Z_NOP()
#endif

#define int_copy_at(_s_, _sat_, _d_, _dat_, _c_, _ex_)  do { \
  for (_i = 0; _i < (_c_); ++_i) { \
    *(((int *) (_d_)) + (_dat_) + _i) = *(((int *) (_s_)) + (_sat_) + _i); \
  } \
} while (0)

/*    printf("dat: %d, sat: %d\n", (_dat_), (_sat_)); \*/


int ZMPI_Alltoall_2step_int(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_2step_int */
{
  int i, j, k, l;
  int _i;

  int size, rank;
  int sq, size_div, size_mod, size_fac, rank_div, rank_mod;

  int *scounts, *sdispls, *rcounts, *rdispls;

  char *tempbuf0, *tempbuf1;

  MPI_Aint stype_lb, stype_extent, rtype_lb, rtype_extent;  

#ifdef WITH_TIMING
  double t[6];
#endif


  TSTART(t[0], comm);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  sq = (int) ceil(sqrt(size));

  size_div = size / sq;
  size_mod = size % sq;
  size_fac = size_div + (size_mod > 0);

  rank_mod = rank % sq;
  rank_div = rank / sq;

  scounts = alloca(4 * size * sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

  MPI_Type_get_extent(sendtype, &stype_lb, &stype_extent);
  MPI_Type_get_extent(recvtype, &rtype_lb, &rtype_extent);

  tempbuf0 = alloca((size_fac + 1) * sq * sendcount * stype_extent);
  tempbuf1 = tempbuf0 + size_fac * sq * sendcount * stype_extent;

#ifdef PRINT_BUF
  memset(recvbuf, 0, size * recvcount * rtype_extent);
  memset(tempbuf0, 0, (size_fac + 1) * sq * sendcount * stype_extent);
#endif

#ifdef PRINT_BUF
  if (rank == 0) printf("input:\n");
  testset_print(size, sendcount, (int *) sendbuf, size, rank, comm);
#endif

  TSTART(t[1], comm);

  memset(scounts, 0, 4 * size * sizeof(int));
  if (rank_div < size_div)
  {
    k = 0;
    for (i = 0; i < sq; ++i)
    {
      l = (i < size_mod)?size_fac:size_div;
      scounts[rank_div * sq + i] = l * sendcount;
      rcounts[rank_div * sq + i] = ((rank_mod < size_mod)?size_fac:size_div) * sendcount;
      sdispls[rank_div * sq + i] = k;
      rdispls[rank_div * sq + i] = i * size_fac * recvcount;
      for (j = 0; j < l; ++j) { int_copy_at(sendbuf, (j * sq + i) * sendcount, recvbuf, k, sendcount, stype_extent); k += sendcount; }
    }

  } else
  {
    k = 0;
    for (i = 0; i < size_mod; ++i)
    {
      scounts[rank_div * sq + i] = rcounts[rank_div * sq + i] = size_fac * sendcount;
      sdispls[rank_div * sq + i] = rdispls[rank_div * sq + i] = k;
      for (j = 0; j < size_fac; ++j) { int_copy_at(sendbuf, (j * sq + i) * sendcount, recvbuf, k, sendcount, stype_extent); k += sendcount; }
    }
  }

#ifdef PRINT_BUF
  if (rank == 0) printf("packed:\n");
  testset_print(size, sendcount, (int *) recvbuf, size, rank, comm);
#endif

  TSTOP(t[1], comm);

  TSTART(t[2], comm);
  MPI_Alltoallv(recvbuf, scounts, sdispls, sendtype, tempbuf0, rcounts, rdispls, recvtype, comm);
  TSTOP(t[2], comm);

#ifdef PRINT_BUF
  if (rank == 0) printf("local transpose packed:\n");
  testset_print(size_fac * sq, recvcount, (int *) tempbuf0, size, rank, comm);
#endif

  TSTART(t[3], comm);

  memset(scounts, 0, 4 * size * sizeof(int));
  if (rank_div < size_div)
  {
    for (i = 0; i < size_div; ++i)
    {
      scounts[i * sq + rank_mod] = rcounts[i * sq + rank_mod] = sq * recvcount;
      sdispls[i * sq + rank_mod] = rdispls[i * sq + rank_mod] = i * sq * recvcount;
      
      for (j = 0; j < sq; ++j) int_copy_at(tempbuf0, (j * size_fac + i) * recvcount, recvbuf, (i * sq + j) * recvcount, recvcount, rtype_extent);
    }

    if (size_mod > 0)
    {
      if (rank_mod < size_mod)
      {
        scounts[size_div * sq + rank_mod] = sq * recvcount;
        rcounts[size_div * sq + rank_mod] = size_mod * recvcount;

        sdispls[size_div * sq + rank_mod] = rdispls[size_div * sq + rank_mod] = size_div * sq * recvcount;

        for (j = 0; j < sq; ++j) int_copy_at(tempbuf0, (j * size_fac + size_div) * recvcount, tempbuf1, j * recvcount, recvcount, rtype_extent);

      } else
      {
        for (j = 0; j < size_mod; ++j)
        {
          rcounts[size_div * sq + j] = recvcount;
          rdispls[size_div * sq + j] = (size_div * sq + j) * recvcount;
        }
      }

      int_copy_at(tempbuf1, 0, tempbuf0, size_div * sq * recvcount, sq * recvcount, rtype_extent);
    }

    int_copy_at(recvbuf, 0, tempbuf0, 0, size_div * sq * recvcount, rtype_extent);

  } else
  {
    for (i = 0; i < size_div; ++i)
    {
      scounts[i * sq + rank_mod] = size_mod * recvcount;
      rcounts[i * sq + rank_mod] = sq * recvcount;

      sdispls[i * sq + rank_mod] = rdispls[i * sq + rank_mod] = i * sq * recvcount;

      for (j = 0; j < size_mod; ++j) int_copy_at(tempbuf0, (j * size_fac + i) * recvcount, recvbuf, (i * sq + j) * recvcount, recvcount, rtype_extent);

      for (j = size_mod; j < sq; ++j)
      {
        scounts[i * sq + j] = recvcount;
        sdispls[i * sq + j] = (i * sq + j) * recvcount;
        
        int_copy_at(sendbuf, (i * sq + j) * recvcount, recvbuf, (i * sq + j) * recvcount, recvcount, rtype_extent);
      }
    }

    scounts[size_div * sq + rank_mod] = rcounts[size_div * sq + rank_mod] = size_mod * recvcount;
    sdispls[size_div * sq + rank_mod] = rdispls[size_div * sq + rank_mod] = size_div * sq * recvcount;
    
    for (j = 0; j < size_mod; ++j) int_copy_at(tempbuf0, (j * size_fac + size_div) * recvcount, recvbuf, (size_div * sq + j) * recvcount, recvcount, rtype_extent);

    int_copy_at(recvbuf, 0, tempbuf0, 0, size * recvcount, rtype_extent);
  }

  TSTOP(t[3], comm);

#ifdef PRINT_BUF
  if (rank == 0) printf("local transpose:\n");
  testset_print(size_fac * sq, recvcount, (int *) tempbuf0, size, rank, comm);
#endif

  TSTART(t[4], comm);
  MPI_Alltoallv(tempbuf0, scounts, sdispls, recvtype, recvbuf, rcounts, rdispls, recvtype, comm);
  TSTOP(t[4], comm);

#ifdef PRINT_BUF
  if (rank == 0) printf("global transpose:\n");
  testset_print(size, recvcount, recvbuf, size, rank, comm);
#endif

  TSTOP(t[0], comm);

#ifdef WITH_TIMING
  if (rank == 0) printf("%d: ZMPI_Alltoall_2step_int: %d  %f  %f  %f  %f  %f\n", rank, sendcount, t[0], t[1], t[2], t[3], t[4]);
#endif

  return MPI_SUCCESS;
}

#undef TSTART
#undef TSTOP
#undef int_copy_at
#undef WITH_TIMING



#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "z_pack.h"

#include "zmpi_tools.h"


#define DEFAULT_INT  0


int ZMPI_Alltoall_int_alltoall(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_alltoall */
{
  return MPI_Alltoall(sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, comm);
}


int ZMPI_Alltoall_int_2step(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_2step */
{
  return ZMPI_Alltoall_2step_int(sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, comm);
}


static int _ZMPI_Alltoall_int_proclists_put(int alloc_mem, int nphases, int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm)
{
  int i, p, size, rank, *rcounts_put;

  MPI_Win win;


  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (alloc_mem) MPI_Alloc_mem(size * sizeof(int), MPI_INFO_NULL, &rcounts_put);
  else rcounts_put = recvbuf;

  if (nrprocs >= 0)
    for (i = 0; i < nrprocs; ++i) rcounts_put[rprocs[i]] = DEFAULT_INT;
  else
    for (i = 0; i < size; ++i) rcounts_put[i] = DEFAULT_INT;

  MPI_Win_create(rcounts_put, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &win);
  MPI_Win_fence(MPI_MODE_NOSTORE|MPI_MODE_NOPRECEDE, win);

  for (p = 0; p < nphases; ++p)
  {
/*    printf("%d: phase = %d of %d\n", rank, p, nphases);*/
  
    if (rank % nphases == p)
    {
      if (nsprocs >= 0)
      {
        for (i = 0; i < nsprocs; ++i)
          if (sendbuf[sprocs[i]] != DEFAULT_INT) MPI_Put(&sendbuf[sprocs[i]], 1, MPI_INT, sprocs[i], rank, 1, MPI_INT, win);

      } else
      {
        for (i = 0; i < size; ++i)
          if (sendbuf[i] != DEFAULT_INT) MPI_Put(&sendbuf[i], 1, MPI_INT, i, rank, 1, MPI_INT, win);
      }
    }

    if (p < nphases - 1) MPI_Win_fence(0, win);
  }

  MPI_Win_fence(MPI_MODE_NOPUT|MPI_MODE_NOSUCCEED, win);
  MPI_Win_free(&win);

  if (alloc_mem)
  {
    if (nrprocs >= 0)
      for (i = 0; i < nrprocs; ++i) recvbuf[rprocs[i]] = rcounts_put[rprocs[i]];
    else
      for (i = 0; i < size; ++i) recvbuf[i] = rcounts_put[i];

    MPI_Free_mem(rcounts_put);    
  }

  return MPI_SUCCESS;
}


int ZMPI_Alltoall_int_put(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 1, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_alloc(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 1, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_2phases(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_2phases */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 2, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_2phases_alloc(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_2phases_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 2, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_3phases(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_3phases */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 3, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_3phases_alloc(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_3phases_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 3, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_proclists_isendirecv(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_isendirecv */
{
  return ZMPI_Alltoall_proclists_isendirecv(sendbuf, 1, MPI_INT, nsprocs, sprocs, recvbuf, 1, MPI_INT, nrprocs, rprocs, comm);
}


int ZMPI_Alltoall_int_proclists_alltoallv(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_alltoallv */
{
  int i, size;

  int *scounts2, *sdispls2, *rcounts2, *rdispls2;


  MPI_Comm_size(comm, &size);

  scounts2 = z_alloc(4 * size, sizeof(int));
  sdispls2 = scounts2 + 1 * size;
  rcounts2 = scounts2 + 2 * size;
  rdispls2 = scounts2 + 3 * size;

  for (i = 0; i < size; ++i)
  {
    scounts2[i] = rcounts2[i] = DEFAULT_INT;
    sdispls2[i] = rdispls2[i] = i;
    recvbuf[i] = 0;
  }

  for (i = 0; i < nsprocs; ++i) scounts2[sprocs[i]] = 1;
  for (i = 0; i < nrprocs; ++i) rcounts2[rprocs[i]] = 1;

  MPI_Alltoallv(sendbuf, scounts2, sdispls2, MPI_INT, recvbuf, rcounts2, rdispls2, MPI_INT, comm);

  z_free(scounts2);

  return MPI_SUCCESS;
}


int ZMPI_Alltoall_int_proclists_put(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_put */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 1, sendbuf, nsprocs, sprocs, recvbuf, nrprocs, rprocs, comm);
}


int ZMPI_Alltoall_int_proclists_put_alloc(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_put_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 1, sendbuf, nsprocs, sprocs, recvbuf, nrprocs, rprocs, comm);
}


int ZMPI_Alltoall_int_proclists_put_2phases(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_put_2phases */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 2, sendbuf, nsprocs, sprocs, recvbuf, nrprocs, rprocs, comm);
}


int ZMPI_Alltoall_int_proclists_put_2phases_alloc(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_put_2phases_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 2, sendbuf, nsprocs, sprocs, recvbuf, nrprocs, rprocs, comm);
}



#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "z_pack.h"

#include "zmpi_tools.h"


int ZMPI_Alltoall_proclists_isendirecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int recvcount, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_proclists_isendirecv */
{
  int i, j;

  const int tag = 0;

  int nreqs;
  MPI_Request *reqs;
  MPI_Status *stats;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;


  reqs = z_alloc(nrecvprocs + nsendprocs, sizeof(MPI_Request) + sizeof(MPI_Status));
  stats = (MPI_Status *) (reqs + nrecvprocs + nsendprocs);

  MPI_Type_get_extent(sendtype, &sendtype_lb, &sendtype_extent);
  MPI_Type_get_extent(recvtype, &recvtype_lb, &recvtype_extent);

  nreqs = 0;

  for (i = 0; i < nrecvprocs; ++i)
  {
    j = recvprocs[i];
    MPI_Irecv(((char *) recvbuf) + (j * recvcount * recvtype_extent), recvcount, recvtype, j, tag, comm, &reqs[nreqs]);
    ++nreqs;
  }

  for (i = 0; i < nsendprocs; ++i)
  {
    j = sendprocs[i];
    MPI_Isend(((char *) sendbuf) + (j * sendcount * sendtype_extent), sendcount, sendtype, j, tag, comm, &reqs[nreqs]);
    ++nreqs;
  }

  MPI_Waitall(nreqs, reqs, stats);

  z_free(reqs);

  return MPI_SUCCESS;
}


int ZMPI_Alltoall_proclists(void *sendbuf, int sendcount, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int recvcount, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_proclists */
{
  return ZMPI_Alltoall_proclists_isendirecv(sendbuf, sendcount, sendtype, nsendprocs, sendprocs, recvbuf, recvcount, recvtype, nrecvprocs, recvprocs, comm);
}


int ZMPI_Alltoallv_proclists_isendirecv(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallv_proclists_isendirecv */
{
  int i, j;

  const int tag = 0;

  int nreqs;
  MPI_Request *reqs;
  MPI_Status *stats;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;


  reqs = z_alloc(nrecvprocs + nsendprocs, sizeof(MPI_Request) + sizeof(MPI_Status));
  stats = (MPI_Status *) (reqs + nrecvprocs + nsendprocs);

  MPI_Type_get_extent(sendtype, &sendtype_lb, &sendtype_extent);
  MPI_Type_get_extent(recvtype, &recvtype_lb, &recvtype_extent);

  nreqs = 0;

  for (i = 0; i < nrecvprocs; ++i)
  {
    j = recvprocs[i];
    if (recvcounts[j] > 0)
    {
      MPI_Irecv(((char *) recvbuf) + (recvdispls[j] * recvtype_extent), recvcounts[j], recvtype, j, tag, comm, &reqs[nreqs]);
      ++nreqs;
    }
  }

  for (i = 0; i < nsendprocs; ++i)
  {
    j = sendprocs[i];
    if (sendcounts[j] > 0)
    {
      MPI_Isend(((char *) sendbuf) + (senddispls[j] * sendtype_extent), sendcounts[j], sendtype, j, tag, comm, &reqs[nreqs]);
      ++nreqs;
    }
  }

  MPI_Waitall(nreqs, reqs, stats);

  z_free(reqs);

  return MPI_SUCCESS;
}


int ZMPI_Alltoallv_proclists(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallv_proclists */
{
  return ZMPI_Alltoallv_proclists_isendirecv(sendbuf, sendcounts, senddispls, sendtype, nsendprocs, sendprocs, recvbuf, recvcounts, recvdispls, recvtype, nrecvprocs, recvprocs, comm);
}


int ZMPI_Alltoallw_proclists_isendirecv(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype *sendtypes, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype *recvtypes, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallw_proclists_isendirecv */
{
  int i, j;

  const int tag = 0;

  int nreqs;
  MPI_Request *reqs;
  MPI_Status *stats;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;


  reqs = z_alloc(nrecvprocs + nsendprocs, sizeof(MPI_Request) + sizeof(MPI_Status));
  stats = (MPI_Status *) (reqs + nrecvprocs + nsendprocs);

  nreqs = 0;

  for (i = 0; i < nrecvprocs; ++i)
  {
    j = recvprocs[i];
    if (recvcounts[j] > 0)
    {
      MPI_Type_get_extent(recvtypes[j], &recvtype_lb, &recvtype_extent);
      MPI_Irecv(((char *) recvbuf) + recvdispls[j], recvcounts[j], recvtypes[j], j, tag, comm, &reqs[nreqs]);
      ++nreqs;
    }
  }

  for (i = 0; i < nsendprocs; ++i)
  {
    j = sendprocs[i];
    if (sendcounts[j] > 0)
    {
      MPI_Type_get_extent(sendtypes[j], &sendtype_lb, &sendtype_extent);
      MPI_Isend(((char *) sendbuf) + senddispls[j], sendcounts[j], sendtypes[j], j, tag, comm, &reqs[nreqs]);
      ++nreqs;
    }
  }

  MPI_Waitall(nreqs, reqs, stats);

  z_free(reqs);

  return MPI_SUCCESS;
}


int ZMPI_Alltoallw_proclists(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype *sendtypes, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype *recvtypes, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallw_proclists */
{
  return ZMPI_Alltoallw_proclists_isendirecv(sendbuf, sendcounts, senddispls, sendtypes, nsendprocs, sendprocs, recvbuf, recvcounts, recvdispls, recvtypes, nrecvprocs, recvprocs, comm);
}



#include <mpi.h>

#include "z_pack.h"

#include "zmpi_tools.h"


int ZMPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block */
{
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 2

  return MPI_Reduce_scatter_block((void *) sendbuf, recvbuf, recvcount, datatype,op, comm);

#else

  int comm_size, *recvcounts, i, exit_code;


  MPI_Comm_size(comm, &comm_size);

  recvcounts = z_alloc(comm_size, sizeof(int));

  for (i = 0; i < comm_size; ++i) recvcounts[i] = recvcount;

  exit_code = MPI_Reduce_scatter((void *) sendbuf, recvbuf, recvcounts, datatype, op, comm);

  z_free(recvcounts);

  return exit_code;
#endif
}



#include <string.h>
#include <mpi.h>

#include "z_pack.h"

#include "zmpi_tools.h"


#define DEFAULT_INT  0


static int _ZMPI_Reduce_scatter_block_intsum_accumulate(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm)
{
  int i, j, size, rank;

  MPI_Win win;


  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  for (i = 0; i < recvcount; ++i) recvbuf[i] = DEFAULT_INT;

  MPI_Win_create(recvbuf, recvcount * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &win);
  MPI_Win_fence(MPI_MODE_NOSTORE|MPI_MODE_NOPRECEDE, win);

  if (nsendprocs >= 0)
  {
    for (j = 0; j < nsendprocs; ++j)
    {
      for (i = 0; i < recvcount; ++i) if (sendbuf[sendprocs[j] * recvcount + i] != DEFAULT_INT) break;

      if (i < recvcount) MPI_Accumulate((void *) &sendbuf[sendprocs[j] * recvcount], recvcount, MPI_INT, sendprocs[j], 0, recvcount, MPI_INT, MPI_SUM, win);
    }

  } else
  {
    for (j = 0; j < size; ++j)
    {
      for (i = 0; i < recvcount; ++i) if (sendbuf[j * recvcount + i] != DEFAULT_INT) break;

      if (i < recvcount) MPI_Accumulate((void *) &sendbuf[j * recvcount], recvcount, MPI_INT, j, 0, recvcount, MPI_INT, MPI_SUM, win);
    }
  }

  MPI_Win_fence(MPI_MODE_NOPUT|MPI_MODE_NOSUCCEED, win);
  MPI_Win_free(&win);

  return MPI_SUCCESS;
}


int ZMPI_Reduce_scatter_block_intsum_accumulate(const int *sendbuf, int *recvbuf, int recvcount, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block_intsum_accumulate */
{
  return _ZMPI_Reduce_scatter_block_intsum_accumulate(sendbuf, -1, NULL, recvbuf, recvcount, -1, NULL, comm);
}


int ZMPI_Reduce_scatter_block_intsum_proclists_isendirecv(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block_intsum_proclists_isendirecv */
{
  int i, j;

  int *recvbuf_full;
  MPI_Request *reqs;


  recvbuf_full = z_alloc(nrecvprocs * recvcount, sizeof(int));

  reqs = z_alloc(nsendprocs + nrecvprocs, sizeof(MPI_Request));

  for (j = 0; j < nrecvprocs; ++j) MPI_Irecv(&recvbuf_full[j * recvcount], recvcount, MPI_INT, recvprocs[j], 0, comm, &reqs[j]);

  for (j = 0; j < nsendprocs; ++j) MPI_Isend((void *) &sendbuf[sendprocs[j] * recvcount], recvcount, MPI_INT, sendprocs[j], 0, comm, &reqs[nrecvprocs + j]);

  MPI_Waitall(nsendprocs + nrecvprocs, reqs, MPI_STATUSES_IGNORE);

  for (i = 0; i < recvcount; ++i) recvbuf[i] = DEFAULT_INT;

  for (j = 0; j < nrecvprocs; ++j)
    for (i = 0; i < recvcount; ++i) recvbuf[i] += recvbuf_full[j * recvcount + i];

  z_free(reqs);

  z_free(recvbuf_full);

  return MPI_SUCCESS;
}


int ZMPI_Reduce_scatter_block_intsum_proclists_alltoallv(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block_intsum_proclists_alltoallv */
{
  int i, j, size, rank;

  int *recvbuf_full;
  int *scounts, *sdispls, *rcounts, *rdispls;


  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  recvbuf_full = z_alloc(nrecvprocs * recvcount, sizeof(int));

  scounts = z_alloc(4 * size, sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

  memset(scounts, 0, 4 * size * sizeof(int));

  for (j = 0; j < nrecvprocs; ++j)
  {
    rcounts[recvprocs[j]] = recvcount;
    rdispls[recvprocs[j]] = j * recvcount;
  }

  for (j = 0; j < nsendprocs; ++j)
  {
    scounts[sendprocs[j]] = recvcount;
    sdispls[sendprocs[j]] = sendprocs[j] * recvcount;
  }

  MPI_Alltoallv((void *) sendbuf, scounts, sdispls, MPI_INT, recvbuf_full, rcounts, rdispls, MPI_INT, comm);

  for (i = 0; i < recvcount; ++i) recvbuf[i] = DEFAULT_INT;

  for (j = 0; j < nrecvprocs; ++j)
    for (i = 0; i < recvcount; ++i) recvbuf[i] += recvbuf_full[j * recvcount + i];

  z_free(scounts);

  z_free(recvbuf_full);

  return MPI_SUCCESS;
}


int ZMPI_Reduce_scatter_block_intsum_proclists_accumulate(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block_intsum_proclists_accumulate */
{
  return _ZMPI_Reduce_scatter_block_intsum_accumulate(sendbuf, nsendprocs, sendprocs, recvbuf, recvcount, nrecvprocs, recvprocs, comm);
}
