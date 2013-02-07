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

#define COUNT 1024*1024

int main(int argc, char* argv[])
{
    int provided;
    int i, rank, nranks, msgsize, target;
    long bufsize;
    int **counter;
    int *complete;
    int increment;
    int counter_fetch;
    int counters_received;
    int t_start, t_stop, t_latency;
    int expected;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    ARMCI_Init_args(&argc, &argv);

    complete = (int *) malloc(sizeof(int) * COUNT);

    counter = (int**) ARMCI_Malloc_local( nranks * sizeof(int*) );
    ARMCI_Malloc((void **) counter[rank], sizeof(int));

    if (rank == 0)
    {
        printf("ARMCI_RMW Test - in usec \n");
        fflush(stdout);
    }

    target = 0; 

    for(i=0; i<COUNT; i++)
    {
       complete[i] = 0;
    } 
    if(rank == target) 
    { 
       *(counter[rank]) = 0;
    }
    increment = 1;
    counter_fetch = 0;
    counters_received = 0;

    MPI_Barrier(MPI_COMM_WORLD);
 
    while(counter_fetch < COUNT)
    {  
        ARMCI_Rmw(ARMCI_FETCH_AND_ADD,
                  (void *) &counter_fetch,
                  (void *) counter[target],
                  increment,
                  target);

        /* s/1/rank/ means we will know who got the counter */
        if (counter_fetch < COUNT) complete[counter_fetch] = rank;
        counters_received++;
    }

    MPI_Allreduce(MPI_IN_PLACE,complete,COUNT,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    for(i=0; i<COUNT; i++)
    {
       if (complete[i] == 0)
       {
           printf("[%d] The RMW update failed at index: %d \n", rank, i);
           fflush(stdout);
           exit(-1);
       }   
    }
    printf("[%d] The RMW update completed successfully \n", rank);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (0==rank)
    {
        printf("Checking for fairness...\n", rank);
        fflush(stdout);
        for(i=0; i<COUNT; i++)
        {
           printf("counter value %d was received by process %d\n", i, complete[i]);
        }
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    printf("process %d received %d counters\n", rank, counters_received);
    fflush(stdout);

    ARMCI_Free(counter[rank]);
    ARMCI_Free_local(counter);

    ARMCI_Finalize();

    MPI_Finalize();

    return 0;
}
