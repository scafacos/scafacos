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


#define ZMPI_ALLTOALLV_INPLACE_AUX_TYPE_STATIC  0
#define ZMPI_ALLTOALLV_INPLACE_AUX_TYPE_HEAP    1

#define ZMPI_SYM_TYPE_SKIP               -1
#define ZMPI_SYM_TYPE_LINEAR              0
#define ZMPI_SYM_TYPE_LINEAR_RANDOM       1
#define ZMPI_SYM_TYPE_HIERARCHIC          2
#define ZMPI_SYM_TYPE_HIERARCHIC_RANDOM   3

#define ZMPI_ALLTOALLV_INPLACE_SYM_TYPE_DEFAULT  ZMPI_SYM_TYPE_HIERARCHIC

#define ZMPI_ALLTOALL_INPLACE_SYMMETRIC_SYM_TYPE_DEFAULT  ZMPI_SYM_TYPE_HIERARCHIC

#define ZMPI_ALLTOALLV_INPLACE_SYMMETRIC_SYM_TYPE_DEFAULT  ZMPI_SYM_TYPE_HIERARCHIC


/* ZMPI_Alltoallv_inplace and ZMPI_Neighbor_alltoallv_inplace */
extern int ZMPI_Alltoallv_inplace_sym_type;
extern int ZMPI_Alltoallv_inplace_aux_type;
extern void *ZMPI_Alltoallv_inplace_aux;
extern int ZMPI_Alltoallv_inplace_aux_size, ZMPI_Alltoallv_inplace_aux_static_blocks;
extern int ZMPI_Alltoallv_inplace_sync_on_init, ZMPI_Alltoallv_inplace_sync_on_exit, ZMPI_Alltoallv_inplace_sync_run;

extern double ZMPI_Alltoallv_inplace_times[10];

int ZMPI_Alltoallv_inplace(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm);

int ZMPI_Neighbor_alltoallv_inplace(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm);

/* ZMPI_Alltoall_inplace_symmetric and ZMPI_Neighbor_alltoall_inplace_symmetric */
extern int ZMPI_Alltoall_inplace_symmetric_sym_type;

int ZMPI_Alltoall_inplace_symmetric(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, MPI_Comm comm);

int ZMPI_Neighbor_alltoall_inplace_symmetric(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, MPI_Comm comm);

/* ZMPI_Alltoallv_inplace_symmetric and ZMPI_Neighbor_alltoallv_inplace_symmetric */
extern int ZMPI_Alltoallv_inplace_symmetric_sym_type;
extern void *ZMPI_Alltoallv_inplace_symmetric_aux;
extern int ZMPI_Alltoallv_inplace_symmetric_aux_size;

int ZMPI_Alltoallv_inplace_symmetric(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm);

int ZMPI_Neighbor_alltoallv_inplace_symmetric(void *sbuf, int *scounts, int *sdispls, MPI_Datatype stype, void *rbuf, int *rcounts, int *rdispls, MPI_Datatype rtype, MPI_Comm comm);


#endif /* __ZMPI_ATAIP_H__ */
