/* 
 * The following is a notice of limited availability of the code, and disclaimer
 * which must be included in the prologue of the code and in all source listings
 * of the code.
 *
 * Copyright (c) 2010  Argonne Leadership Computing Facility, Argonne National
 * Laboratory
 *
 * Permission is hereby granted to use, reproduce, prepare derivative works, and
 * to redistribute to others.
 *
 *
 *                          LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer listed
 *   in this license in the documentation and/or other materials
 *   provided with the distribution.
 *
 * - Neither the name of the copyright holders nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include "parmci.h"

#if defined(__bgp__)
typedef int32_t atomic_int_t;
#elif defined(__CRAYXE)
typedef int64_t atomic_int_t;
#else
typedef int64_t atomic_int_t;
#endif

int main(int argc, char *argv[])
{
    int rank, size;
    int provided;

#if defined(__bgp__)
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    assert(provided==MPI_THREAD_MULTIPLE);
#else
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    //assert(provided>MPI_THREAD_SINGLE);
#endif
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert( size > 1 );

    PARMCI_Init_args(&argc, &argv);

    int w = ( argc > 1 ? atoi(argv[1]) : 10 );

    if ( rank == 0 ) printf( "testing RMW: %d ranks, %d integers\n", size, w );

    atomic_int_t ** window;
    window  = (atomic_int_t **) PARMCI_Malloc_local( size * sizeof(atomic_int_t *) );
    PARMCI_Malloc( (void **) window, w * sizeof(atomic_int_t) );
    for (int i = 0; i < w; i++) window[rank][i] = rank*1000+i;

    atomic_int_t * buffer;
    buffer = (atomic_int_t *) PARMCI_Malloc_local(  w * sizeof(atomic_int_t) );
    for (int i = 0; i < w; i++) buffer[i] = -17037L;

    PARMCI_Barrier();

    if (rank == 0)
        for (int t=1; t<size; t++)
        {
            int bytes = w * sizeof(atomic_int_t);

            for (int i = 0; i < w; i++) buffer[i] = (atomic_int_t)(1000000+t);

            for (int i = 0; i < w; i++) 
                printf("rank %d (before put+fence) buffer[%d] = %ld \n", rank, i, (long)buffer[i] );

            PARMCI_Put( buffer, window[t], bytes, t );
            PARMCI_Fence( t );

            for (int i = 0; i < w; i++) buffer[i] = 11111;

            PARMCI_Get( window[t], buffer, bytes, t );

            for (int i = 0; i < w; i++) 
                printf("rank %d (after get 1) buffer[%d] = %d \n", rank, i, buffer[i] );

            for (int i = 0; i < w; i++) 
#if defined(__bgp__)
                PARMCI_Rmw( ARMCI_FETCH_AND_ADD,      &buffer[i], &window[t][i], (atomic_int_t)1000, t );
#elif defined(__CRAYXE)
                PARMCI_Rmw( ARMCI_FETCH_AND_ADD_LONG, &buffer[i], &window[t][i], (atomic_int_t)1000, t );
#endif

            for (int i = 0; i < w; i++) 
                printf("rank %d (after atomic fetch-and-add-1000) buffer[%d] = %ld \n", rank, i, (long)buffer[i] );

            for (int i = 0; i < w; i++) buffer[i] = 55555;

            PARMCI_Get( window[t], buffer, bytes, t );

            for (int i = 0; i < w; i++) 
                printf("rank %d (after get 2) buffer[%d] = %ld \n", rank, i, (long)buffer[i] );
        }

    PARMCI_Barrier();

    PARMCI_Free_local( (void *) buffer );

    PARMCI_Free( (void *) window[rank] );
    PARMCI_Free_local( (void *) window );

    PARMCI_Finalize();

    printf("%d: all done \n", rank );
    fflush(stdout);


    MPI_Finalize();

    return 0;
}
