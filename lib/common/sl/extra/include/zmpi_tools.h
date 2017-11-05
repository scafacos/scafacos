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


#ifndef __ZMPI_TOOLS_H__
#define __ZMPI_TOOLS_H__





#ifdef HAVE_CONFIG_H
# include "config.h"
#endif


#define HAVE_ZMPI_ALLTOALL_2STEP
#define HAVE_ZMPI_ALLTOALLX_PROCLISTS
#define HAVE_ZMPI_ALLTOALL_INT
#define HAVE_ZMPI_REDUCE_SCATTER_BLOCK
#define HAVE_ZMPI_REDUCE_SCATTER_BLOCK_INTSUM



#ifdef ZMPI_PREFIX
# include "zmpi_tools_rename.h"
#endif


#ifdef HAVE_ZMPI_ALLTOALL_2STEP
int ZMPI_Alltoall_2step_int(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
#endif

#ifdef HAVE_ZMPI_ALLTOALLX_ISENDIRECV
extern int ZMPI_Alltoallx_isendirecv_max_procs;
int ZMPI_Alltoall_isendirecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int ZMPI_Alltoallv_isendirecv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);
int ZMPI_Alltoallw_isendirecv(void* sendbuf, int sendcounts[], int sdispls[], MPI_Datatype sendtypes[], void *recvbuf, int recvcounts[], int rdispls[], MPI_Datatype recvtypes[], MPI_Comm comm);
#endif

#ifdef HAVE_ZMPI_ALLTOALLX_SENDRECV
int ZMPI_Alltoall_sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int ZMPI_Alltoallv_sendrecv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);
int ZMPI_Alltoallw_sendrecv(void* sendbuf, int sendcounts[], int sdispls[], MPI_Datatype sendtypes[], void *recvbuf, int recvcounts[], int rdispls[], MPI_Datatype recvtypes[], MPI_Comm comm);
#endif

#ifdef HAVE_ZMPI_ALLTOALLX_PROCLISTS
int ZMPI_Alltoall_proclists_isendirecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int recvcount, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm);
int ZMPI_Alltoall_proclists(void *sendbuf, int sendcount, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int recvcount, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm);
int ZMPI_Alltoallv_proclists_isendirecv(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm);
int ZMPI_Alltoallv_proclists(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm);
int ZMPI_Alltoallw_proclists_isendirecv(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype *sendtypes, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype *recvtypes, int nrecvprocs, int *recvprocs, MPI_Comm comm);
int ZMPI_Alltoallw_proclists(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype *sendtypes, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype *recvtypes, int nrecvprocs, int *recvprocs, MPI_Comm comm);
#endif

#ifdef HAVE_ZMPI_ALLTOALL_INT
int ZMPI_Alltoall_int_alltoall(int *sendbuf, int *recvbuf, MPI_Comm comm);
int ZMPI_Alltoall_int_2step(int *sendbuf, int *recvbuf, MPI_Comm comm);
int ZMPI_Alltoall_int_put(int *sendbuf, int *recvbuf, MPI_Comm comm);
int ZMPI_Alltoall_int_put_alloc(int *sendbuf, int *recvbuf, MPI_Comm comm);
int ZMPI_Alltoall_int_put_2phases(int *sendbuf, int *recvbuf, MPI_Comm comm);
int ZMPI_Alltoall_int_put_2phases_alloc(int *sendbuf, int *recvbuf, MPI_Comm comm);
int ZMPI_Alltoall_int_put_3phases(int *sendbuf, int *recvbuf, MPI_Comm comm);
int ZMPI_Alltoall_int_put_3phases_alloc(int *sendbuf, int *recvbuf, MPI_Comm comm);
int ZMPI_Alltoall_int_proclists_isendirecv(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm);
int ZMPI_Alltoall_int_proclists_alltoallv(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm);
int ZMPI_Alltoall_int_proclists_put(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm);
int ZMPI_Alltoall_int_proclists_put_alloc(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm);
int ZMPI_Alltoall_int_proclists_put_2phases(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm);
int ZMPI_Alltoall_int_proclists_put_2phases_alloc(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm);
#endif

#ifdef HAVE_ZMPI_REDUCE_SCATTER_BLOCK
int ZMPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
#endif

#ifdef HAVE_ZMPI_REDUCE_SCATTER_BLOCK_INTSUM
int ZMPI_Reduce_scatter_block_intsum_accumulate(const int *sendbuf, int *recvbuf, int recvcount, MPI_Comm comm);
int ZMPI_Reduce_scatter_block_intsum_proclists_isendirecv(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm);
int ZMPI_Reduce_scatter_block_intsum_proclists_alltoallv(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm);
int ZMPI_Reduce_scatter_block_intsum_proclists_accumulate(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm);
#endif


#endif /* __ZMPI_TOOLS_H__ */
