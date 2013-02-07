/********************************************************************
 * The following is a notice of limited availability of the code, and disclaimer
 * which must be included in the prologue of the code and in all source listings
 * of the code.
 *
 * Copyright (c) 2010 Argonne Leadership Computing Facility, Argonne National Laboratory
 *
 * Permission is hereby granted to use, reproduce, prepare derivative works, and
 * to redistribute to others.
 *
 *                 LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer listed
 *    in this license in the documentation and/or other materials
 *    provided with the distribution.
 *
 *  - Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * The copyright holders provide no reassurances that the source code
 * provided does not infringe any patent, copyright, or any other
 * intellectual property rights of third parties.  The copyright holders
 * disclaim any liability to any recipient for claims brought against
 * recipient by any third party for infringement of that parties
 * intellectual property rights.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************/

#include "a1d_coll.h"

/*********************************************************************/

int A1D_Barrier(MPI_Comm comm)
{
    int mpi_status = MPI_SUCCESS;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Barrier() \n");
#endif

    mpi_status = MPI_Barrier(comm);
    assert(mpi_status==MPI_SUCCESS);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Barrier() \n");
#endif

    return(0);
}

int A1D_Allgather(MPI_Comm comm, void * local, void * gout, int local_bytes )
{
    int mpi_status = MPI_SUCCESS;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Allgather(void * local, void * gout, int local_bytes ) \n");
#endif

    mpi_status = MPI_Allgather( local, local_bytes, MPI_BYTE, gout, local_bytes, MPI_BYTE, comm );
    assert(mpi_status==MPI_SUCCESS);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Allgather(void * local, void * gout, int local_bytes ) \n");
#endif

    return(0);
}

int A1D_Allreduce_max32(MPI_Comm comm, int32_t in, int32_t * out)
{
    int mpi_status = MPI_SUCCESS;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Allreduce_max32(int32_t in, int32_t * out) \n");
#endif

#ifdef __bgp__
    mpi_status = MPI_Allreduce( &in, out, 1, MPI_INT, MPI_MAX, comm );
#else
    mpi_status = MPI_Allreduce( &in, out, 1, MPI_INT32_T, MPI_MAX, comm );
#endif
    assert(mpi_status==MPI_SUCCESS);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Allreduce_max32(int32_t in, int32_t * out) \n");
#endif

    return(0);
}

int A1D_Allreduce_issame32(MPI_Comm comm, int32_t value, int * flag)
{
    int mpi_status = MPI_SUCCESS;
    int32_t in[2], out[2];

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Allreduce_issame32(int32_t value, int * flag) \n");
#endif

    *flag = 0;

    in[0]  = value;
    in[1]  = -value;
    out[0] = 0;
    out[1] = 0;

#ifdef __bgp__
    mpi_status = MPI_Allreduce( in, out, 2, MPI_INT, MPI_MAX, comm );
#else
    mpi_status = MPI_Allreduce( in, out, 2, MPI_INT32_T, MPI_MAX, comm );
#endif
    assert(mpi_status==MPI_SUCCESS);

    if ( (out[0] == value) && (out[1] = -value) ) 
        (*flag)=1;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Allreduce_issame32(int32_t value, int * flag) \n");
#endif

    return(0);
}

int A1D_Allreduce_issame64(MPI_Comm comm, int64_t value, int * flag)
{
    int mpi_status = MPI_SUCCESS;
    int64_t in[2], out[2];

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Allreduce_issame64(int64_t value, int * flag) \n");
#endif

    *flag = 0;

    in[0]  = value;
    in[1]  = -value;
    out[0] = 0;
    out[1] = 0;

#ifdef __bgp__
    mpi_status = MPI_Allreduce( in, out, 2, MPI_LONG, MPI_MAX, comm );
#else
    mpi_status = MPI_Allreduce( in, out, 2, MPI_INT64_T, MPI_MAX, comm );
#endif
    assert(mpi_status==MPI_SUCCESS);

    if ( (out[0] == value) && (out[1] = -value) ) 
        (*flag)=1;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Allreduce_issame64(int64_t value, int * flag) \n");
#endif

    return(0);
}
