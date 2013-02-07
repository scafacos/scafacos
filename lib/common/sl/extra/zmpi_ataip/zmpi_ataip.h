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


#ifndef __ZMPI_ATAIP_H__
#define __ZMPI_ATAIP_H__


#include <mpi.h>


#ifdef ZMPI_PREFIX
# include "zmpi_ataip_rename.h"
#endif


extern void *ZMPI_Alltoallv_inplace_aux;
extern int ZMPI_Alltoallv_inplace_aux_size, ZMPI_Alltoallv_inplace_aux_blocks;
extern int ZMPI_Alltoallv_inplace_sync_on_init, ZMPI_Alltoallv_inplace_sync_on_exit, ZMPI_Alltoallv_inplace_sync_run;
extern double ZMPI_Alltoallv_inplace_times[9];

int ZMPI_Alltoallv_inplace(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm);
int ZMPI_Alltoallv_inplace_static(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm);
int ZMPI_Alltoallv_inplace_heap(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm);


#endif /* __ZMPI_ATAIP_H__ */
