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
#include <mpi.h>
#include <armci.h>

#define MAX_MSG_SIZE 1024*1024
#define ITERATIONS 100
#define SKIP 10

int main(int argc, char **argv)
{
    int provided;
    int i, j, rank, nranks, msgsize;
    long bufsize;
    double scaling;
    double **buffer;
    double t_start, t_stop, t_latency = 0;
    int count[2], src_stride, trg_stride, stride_level;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    ARMCI_Init_args(&argc, &argv);

    buffer = (double **) malloc(sizeof(double *) * nranks);

    bufsize = MAX_MSG_SIZE * (ITERATIONS + SKIP);
    ARMCI_Malloc((void **) buffer, bufsize);

    scaling = 2.0;
    for (i = 0; i < bufsize / sizeof(double); i++)
    {
        *(buffer[rank] + i) = 1.0 + rank;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("TESTING ALL-FLUSH-ALL\n");
        printf("ARMCI_Put + ARMCI_Barrier Latency - in usec \n");
        printf("%20s %22s\n", "Message Size", "Latency");
        fflush(stdout);
    }

    for (msgsize = sizeof(double); msgsize < MAX_MSG_SIZE; msgsize *= 2)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        for (i = 0; i < ITERATIONS + SKIP; i++)
        {
            for (j = 0; j < nranks; j++)
            {
                ARMCI_Put((void *) ((size_t) buffer[rank] + (size_t)(i * msgsize)),
                          (void *) ((size_t) buffer[j] + (size_t)(i * msgsize)),
                          msgsize,
                          j);
            }
            t_start = MPI_Wtime();
            ARMCI_AllFence();
            MPI_Barrier(MPI_COMM_WORLD);
            t_stop = MPI_Wtime();
            if (i >= SKIP) t_latency = t_latency + (t_stop - t_start);
        }
        printf("%20d %20.2f \n", msgsize, ((t_latency) * 1000000) / ITERATIONS);
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (0 == rank)
    {
        printf("\n");
        printf("ARMCI_Acc + ARMCI_Barrier - in usec \n");
        printf("%20s %22s\n", "Message Size", "Latency");
        fflush(stdout);
    }

    stride_level = 0;
    for (msgsize = sizeof(double); msgsize < MAX_MSG_SIZE; msgsize *= 2)
    {
        src_stride = msgsize * sizeof(double);
        trg_stride = msgsize * sizeof(double);
        count[0] = msgsize * sizeof(double);
        count[1] = msgsize;
        MPI_Barrier(MPI_COMM_WORLD);
        for (i = 0; i < ITERATIONS + SKIP; i++)
        {
            for (j = 0; j < nranks; j++)
            {
                ARMCI_AccS(ARMCI_ACC_DBL,
                           (void *) &scaling,
                           (void *) ((size_t) buffer[rank] + (size_t)(i * msgsize)),
                           &src_stride,
                           (void *) ((size_t) buffer[j] + (size_t)(i * msgsize)),
                           &trg_stride,
                           count,
                           stride_level,
                           j);
            }
            t_start = MPI_Wtime();
            ARMCI_AllFence();
            MPI_Barrier(MPI_COMM_WORLD);
            t_stop = MPI_Wtime();
            if (i >= SKIP) t_latency = t_latency + (t_stop - t_start);
        }
        printf("%20d %20.2f \n", msgsize, ((t_latency) * 1000000) / ITERATIONS);
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    ARMCI_Free(buffer[rank]);
    ARMCI_Finalize();
    MPI_Finalize();

    return 0;
}
