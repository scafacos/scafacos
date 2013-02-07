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
#include <parmci.h>

int main(int argc, char *argv[])
{
    int i;
    int rank, size;
    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    PARMCI_Init_args(&argc, &argv);

    int winsize = ( argc > 1 ? atoi(argv[1]) : 100 );
    if ( rank == 0 ) printf( "winsize = %d doubles\n", winsize );

    double ** window = (double **) PARMCI_Malloc_local( size * sizeof(double*) );
    PARMCI_Malloc((void **) window, winsize * sizeof(double) );
    for (i = 0; i < winsize; i++) 
    {
        printf("%d: i=%d \n",rank,i);
        window[rank][i] = (double)rank;
    }

    double * buffer = PARMCI_Malloc_local(winsize);
    for (i = 0; i < winsize; i++) buffer[i] = -1.0;

    PARMCI_Barrier();

    if (rank == 0)
        for (i=1; i<size; i++)
        {
            double t0 = MPI_Wtime();
            PARMCI_Get( window[i], buffer, winsize * sizeof(double), i );
            double t1 = MPI_Wtime();

            for (i = 0; i < winsize; i++) assert( buffer[i] == (double)i );

            double bw = winsize * sizeof(double) / (t1 - t0);
            printf("PARMCI_Get from rank %d to rank %d of %d bytes took %lf seconds (%lf MB/s)\n",i,0,winsize,t1-t0,bw);
            fflush(stdout);
        }

    PARMCI_Barrier();

    PARMCI_Free_local(buffer);

    PARMCI_Free(window[rank]);

    PARMCI_Finalize();

    MPI_Finalize();

    return 0;
}
