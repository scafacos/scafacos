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

#define MAX_MSGSIZE 2*1024*1024
#define ITERATIONS_VERYSMALL 20000 
#define ITERATIONS_SMALL 4000 
#define ITERATIONS_MEDIUM 1000 
#define ITERATIONS_LARGE 500 

int main()
{

    size_t i, rank, nranks, msgsize, dest;
    size_t iterations, max_msgsize;
    int bufsize;
    double **buffer;
    double t_start, t_stop, t_total, d_total;
    double expected, bandwidth;
    A1_handle_t a1_handle;

    max_msgsize = MAX_MSGSIZE;

    A1_Initialize(A1_THREAD_SINGLE);

    rank = A1_Process_id(A1_GROUP_WORLD);
    nranks = A1_Process_total(A1_GROUP_WORLD);

    bufsize = max_msgsize * ITERATIONS_LARGE;
    buffer = (double **) malloc(sizeof(double *) * nranks);
    A1_Alloc_segment((void **) &(buffer[rank]), bufsize);
    A1_Exchange_segments(A1_GROUP_WORLD, (void **) buffer);

    for (i = 0; i < bufsize / sizeof(double); i++)
    {
        *(buffer[rank] + i) = 1.0 + rank;
    }

    A1_Allocate_handle(&a1_handle);

    A1_Barrier_group(A1_GROUP_WORLD);

    if (rank == 0)
    {

        printf("A1_Put Bandwidth in MBPS \n");
        printf("%20s %22s \n", "Message Size", "Bandwidth");
        fflush(stdout);

        dest = 1;
        expected = 1 + dest;

        for (msgsize = sizeof(double); msgsize <= max_msgsize; msgsize *= 2)
        {

            if (msgsize <= 16 * 1024) iterations = ITERATIONS_VERYSMALL;
            else if (msgsize <= 64 * 1024) iterations = ITERATIONS_SMALL;
            else if (msgsize <= 512 * 1024) iterations = ITERATIONS_MEDIUM;
            else iterations = ITERATIONS_LARGE;

            t_start = A1_Time_seconds();

            for (i = 0; i < iterations; i++)
            {

                A1_NbPut(dest, (void *) ((size_t) buffer[dest] + (size_t)(i
                        * msgsize)), (void *) ((size_t) buffer[rank]
                        + (size_t)(i * msgsize)), msgsize, a1_handle);

            }

            A1_Wait_handle(a1_handle);

            t_stop = A1_Time_seconds();
            d_total = (iterations * msgsize) / (1024 * 1024);
            t_total = t_stop - t_start;
            bandwidth = d_total / t_total;
            printf("%20d %20.4lf \n", msgsize, bandwidth);
            fflush(stdout);
           
            A1_Flush(dest);
        }

    }

    A1_Barrier_group(A1_GROUP_WORLD);

    A1_Release_handle(a1_handle);

    A1_Release_segments(A1_GROUP_WORLD, buffer[rank]);

    A1_Finalize();

    return 0;
}
