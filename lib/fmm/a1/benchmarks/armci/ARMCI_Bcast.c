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

int main(int argc, char **argv)
{
    int i, j, rank, nranks, msgsize;
    int *buffer;
    int provided;
    char op = '+';

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    ARMCI_Init_args(&argc, &argv);

    ARMCI_Barrier();

    buffer = (int *) malloc(MAX_MSG_SIZE);

    for(i=0; i<MAX_MSG_SIZE/sizeof(int); i++)
    {
       if(rank == 0) 
          buffer[i] = (2<<20 - 1);
       else
          buffer[i] = 0;
    }

    if(rank == 0)
    {
      printf("Testing functionality of ARMCI_Bcast \n");
      fflush(stdout);
    } 

    for(msgsize=sizeof(int); msgsize<=MAX_MSG_SIZE; msgsize*=2)
    {
       armci_msg_bcast(buffer, msgsize, 0); 

       for(i=0; i<msgsize/sizeof(int); i++) 
       {
          if(buffer[i] != (2<<20 - 1))
          {
             printf("[%d] Validation failed for msg size: %d at index: %d expected: %d actual: %d \n",
                     rank, msgsize, i, (2<<20 - 1), buffer[i]);
             fflush(stdout);
             exit(-1);
          }  
       }

       for(i=0; i<MAX_MSG_SIZE/sizeof(int); i++)
       {
          if(rank == 0)
             buffer[i] = (2<<20 - 1);
          else
             buffer[i] = 0;
       }

       ARMCI_Barrier();

       if(rank == 0)
       {
         printf("Validation successful for msg size: %d\n", msgsize);
         fflush(stdout);
       }
    }

    free(buffer);

    ARMCI_Finalize();

    MPI_Finalize();

    return 0;
}
