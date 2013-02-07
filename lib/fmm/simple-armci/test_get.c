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

    int w, maxwinsize = ( argc > 1 ? atoi(argv[1]) : 1000000 );

    if ( rank == 0 ) printf( "size = %d maxwinsize = %d floats\n", size, maxwinsize );

    for ( w=1 ; w<maxwinsize ; w*=2 )
    {
        float ** window = (float **) PARMCI_Malloc_local( size * sizeof(float *) );
        PARMCI_Malloc( (void **) window, w * sizeof(float) );
        for (int i = 0; i < w; i++) window[rank][i] = (float)(rank);
        //for (i = 0; i < w; i++) printf("window[%d][%d] = %lf \n", rank, i, window[rank][i] );

        float * buffer = (float *) PARMCI_Malloc_local(  w * sizeof(float) );
        for (int i = 0; i < w; i++) buffer[i] = (float)(-rank);
        //for (i = 0; i < w; i++) printf("%d: buffer[%d] = %lf \n", rank, i, buffer[i] );

        PARMCI_Barrier();

        if (rank == 0)
            for (int t=1; t<size; t++)
            {
                int bytes, repeat;
                double t0, t1, dt, bw;

                bytes  = w * sizeof(float);
                repeat = 1; //10 * (maxwinsize/w);

                PARMCI_Get( window[t], buffer, bytes, t );
                for (int i = (w-1); i >=0 ; i--) assert( buffer[i] == (float)t );
                //for (int i = 0; i < w; i++) 
                //    if ( buffer[i] != (float)i ) printf("rank %d buffer[%d] = %lf \n", rank, i, buffer[t] );

                t0 = MPI_Wtime();
                //for (int r=0; r<repeat; r++)
                PARMCI_Get( window[t], buffer, bytes, t );
                t1 = MPI_Wtime();

                dt  = (t1-t0) / repeat;
                bw  = bytes / dt;
                bw /= 1000000.0;
                printf("PARMCI_Get of from rank %4d to rank %4d of %9d bytes took %lf seconds (%lf MB/s)\n", t, 0, bytes, dt, bw);
                fflush(stdout);
            }

        PARMCI_Barrier();

        PARMCI_Free_local( (void *) buffer );

        PARMCI_Free( (void *) window[rank] );
        PARMCI_Free_local( (void *) window );
    }

    PARMCI_Finalize();

    printf("%d: all done \n", rank );
    fflush(stdout);

    MPI_Finalize();

    return 0;
}
