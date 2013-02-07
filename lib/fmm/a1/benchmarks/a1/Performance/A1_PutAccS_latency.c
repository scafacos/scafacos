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
#include <a1.h>

#define MAX_DIM 1024 
#define ITERATIONS 100
#define SKIP 10

int main()
{

    size_t i, j, rank, nranks, msgsize, peer;
    size_t dim;
    unsigned long bufsize;
    double **buffer;
    double t_start, t_stop, t_latency;
    int count[2], src_stride, trg_stride, stride_level;
    double scaling;

    A1_Initialize(A1_THREAD_SINGLE);

    rank = A1_Process_id(A1_GROUP_WORLD);
    nranks = A1_Process_total(A1_GROUP_WORLD);

    buffer = (double **) malloc(sizeof(int32_t *) * nranks);

    bufsize = MAX_DIM * MAX_DIM * sizeof(double);
    A1_Alloc_segment((void **) &(buffer[rank]), bufsize);
    A1_Exchange_segments(A1_GROUP_WORLD, (void **) buffer);

    if (rank == 0)
    {
        printf("A1_PutAccS Latency - local and remote completions - in usec \n");
        printf("%30s %22s %22s\n",
               "Dimensions(array of double)",
               "Local Completion",
               "Remote completion");
        fflush(stdout);
    }

    for (i = 0; i < bufsize / sizeof(double); i++)
    {
        *(buffer[rank] + i) = 1.0 + rank;
    }
    scaling = 2.0;

    src_stride = MAX_DIM * sizeof(double);
    trg_stride = MAX_DIM * sizeof(double);
    stride_level = 1;

    A1_Barrier_group(A1_GROUP_WORLD);

    for (dim = 1; dim <= MAX_DIM; dim *= 2)
    {

        count[0] = dim*sizeof(double);
        count[1] = dim;

            if (rank == 0)
            {

                peer = 1;

                for (i = 0; i < ITERATIONS + SKIP; i++)
                {

                    if (i == SKIP) t_start = A1_Time_seconds();

                    A1_PutAccS(1,
                               stride_level,
                               count,
                               (void *) buffer[rank],
                               &src_stride,
                               (void *) buffer[peer],
                               &trg_stride,
                               A1_DOUBLE,
                               (void *) &scaling);

                }
                t_stop = A1_Time_seconds();
                A1_Flush(1);

                char temp[10];
                sprintf(temp, "%dX%d", dim, dim);
                printf("%30s %20.2f ", temp, ((t_stop - t_start) * 1000000)
                        / ITERATIONS);
                fflush(stdout);

                A1_Barrier_group(A1_GROUP_WORLD);

                A1_Barrier_group(A1_GROUP_WORLD);

                for (i = 0; i < ITERATIONS + SKIP; i++)
                {

                    if (i == SKIP) t_start = A1_Time_seconds();

                    A1_PutAccS(1,
                               stride_level,
                               count,
                               (void *) buffer[rank],
                               &src_stride,
                               (void *) buffer[peer],
                               &trg_stride,
                               A1_DOUBLE,
                               (void *) &scaling);
                    A1_Flush(1);

                }
                t_stop = A1_Time_seconds();
                printf("%20.2f \n", ((t_stop - t_start) * 1000000) / ITERATIONS);
                fflush(stdout);

                A1_Barrier_group(A1_GROUP_WORLD);

                A1_Barrier_group(A1_GROUP_WORLD);

            }
            else
            {

                peer = 0;

                A1_Barrier_group(A1_GROUP_WORLD);

                for (i = 0; i < dim; i++)
                {
                    for (j = 0; j < dim; j++)
                    {
                        if (*(buffer[rank] + i * MAX_DIM + j) != ((1.0 + rank)
                                + scaling * (1.0 + peer) * (ITERATIONS + SKIP)))
                        {
                            printf("Data validation failed at X: %d Y: %d Expected : %f Actual : %f \n",
                                   i,
                                   j,
                                   ((1.0 + rank) + scaling * (1.0 + peer)),
                                   *(buffer[rank] + i * MAX_DIM + j));
                            fflush(stdout);
                            return -1;
                        }
                    }
                }

                for (i = 0; i < bufsize / sizeof(double); i++)
                {
                    *(buffer[rank] + i) = 1.0 + rank;
                }

                A1_Barrier_group(A1_GROUP_WORLD);

                A1_Barrier_group(A1_GROUP_WORLD);

                for (i = 0; i < dim; i++)
                {
                    for (j = 0; j < dim; j++)
                    {
                        if (*(buffer[rank] + i * MAX_DIM + j) != ((1.0 + rank)
                                + scaling * (1.0 + peer) * (ITERATIONS + SKIP)))
                        {
                            printf("Data validation failed at X: %d Y: %d Expected : %f Actual : %f \n",
                                   i,
                                   j,
                                   ((1.0 + rank) + scaling * (1.0 + peer)),
                                   *(buffer[rank] + i * MAX_DIM + j));
                            fflush(stdout);
                            return -1;
                        }
                    }
                }

                for (i = 0; i < bufsize / sizeof(double); i++)
                {
                    *(buffer[rank] + i) = 1.0 + rank;
                }

                A1_Barrier_group(A1_GROUP_WORLD);

            }

    }

    A1_Barrier_group(A1_GROUP_WORLD);

    A1_Release_segments(A1_GROUP_WORLD, (void *) buffer[rank]);
    A1_Free_segment((void *) buffer[rank]);

    A1_Finalize();

    return 0;
}
