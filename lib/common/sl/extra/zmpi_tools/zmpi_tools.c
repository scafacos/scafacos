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
