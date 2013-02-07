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
#include <armci.h>

#define MAX_MSG_SIZE 4*1024*1024

int main()
{

    int i, rank, nranks, msgsize, peer;
    long bufsize;
    double *buffer;
    double expected;

    A1_Initialize(A1_THREAD_SINGLE);

    rank = A1_Process_id(A1_GROUP_WORLD);
    nranks = A1_Process_total(A1_GROUP_WORLD);

    bufsize = MAX_MSG_SIZE;
    buffer = (double *) malloc(bufsize);

    if (rank == 0)
    {
        printf("A1_Allreduce Test - in usec \n");
        fflush(stdout);
    }

    A1_Barrier_group(A1_GROUP_WORLD);

    for (msgsize = sizeof(double); msgsize < MAX_MSG_SIZE; msgsize *= 2)
    {

            for (i = 0; i < bufsize/sizeof(double); i++)
            {
                 buffer[i] = rank + 1.0;
            }

            A1_Allreduce_group(A1_GROUP_WORLD,
                               msgsize/sizeof(double),
                               A1_SUM,
                               A1_DOUBLE,
                               (void *) buffer,
                               (void *) buffer);

            expected = (nranks)*(nranks+1)/2;
            for (i = 0; i < msgsize/sizeof(double); i++)
            {
               if(buffer[i] - expected != 0)
               {
                   printf("[%d] Validation has failed Expected: %lf, Actual: %lf, i: %d \n",
                               rank, expected, buffer[i], i);
                   fflush(stdout);
                   exit(-1);
               }
            }

            printf("[%d] %d message sum reduce successful\n", rank, msgsize);
            fflush(stdout);

            for (i = 0; i < bufsize/sizeof(double); i++)
            {
                  buffer[i] = 1.0;
            }

            A1_Allreduce_group(A1_GROUP_WORLD,
                               msgsize/sizeof(int),
                               A1_PROD,
                               A1_DOUBLE,
                               (void *) buffer,
                               (void *) buffer);

            expected = 1.0;
            for (i = 0; i < msgsize/sizeof(double); i++)
            {
               if(buffer[i] - expected != 0)
               {
                   printf("[%d] Validation has failed Expected: %lf, Actual: %lf, i: %d \n",
                               rank, expected, buffer[i], i);
                   fflush(stdout);
                   exit(-1);
               }
            }

            printf("[%d] %d message product reduce successful\n", rank, msgsize);
            fflush(stdout);

    }

    free(buffer);
    A1_Finalize();

    return 0;
}
