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


#ifndef __ZMPI_TOOLS_H__
#define __ZMPI_TOOLS_H__





#ifdef HAVE_CONFIG_H
# include "config.h"
#endif


#define HAVE_ZMPI_ALLTOALL_2STEP
#define HAVE_ZMPI_ALLTOALLV_PROCLISTS



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

#ifdef HAVE_ZMPI_ALLTOALLV_PROCLISTS
int ZMPI_Alltoallv_proclists(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm);
#endif


#endif /* __ZMPI_TOOLS_H__ */
