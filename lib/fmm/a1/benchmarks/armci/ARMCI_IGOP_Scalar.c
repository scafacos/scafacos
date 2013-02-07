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

void checkanswer(int num, int buffer, int answer)
{
    if ( buffer != answer ) printf("Test %d failed expected: %d actual: %d \n", num, answer, buffer);
    else                    printf("Test %d passed expected: %d actual: %d \n", num, answer, buffer);
    fflush(stdout);
}

int main(int argc, char **argv)
{
    int provided;

    int rank, size;

    char* sum    = "+";
    char* prod   = "*";
    char* max    = "max";
    char* min    = "min";
    char* absmax = "absmax";
    char* absmin = "absmin";
    char* or     = "or";

    int count = 0;
    int buffer;
    int answer;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ARMCI_Init_args(&argc, &argv);

    ARMCI_Barrier();

    /* TESTING SUM OPERATOR */

    if (rank==0) printf("Test %d: sum A\n",count);
    buffer = ( rank==1 ? 1 : 0 );
    armci_msg_igop(&buffer, sizeof(int), sum); 
    answer = 1;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    if (rank==0) printf("Test %d: sum B\n",count);
    buffer = 1;
    armci_msg_igop(&buffer, sizeof(int), sum);
    answer = size;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    /* TESTING PROD OPERATOR */

    if (rank==0) printf("Test %d: prod A\n",count);
    buffer = ( rank==0 ? 0 : 1 );
    armci_msg_igop(&buffer, sizeof(int), prod);
    answer = 0;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    if (rank==0) printf("Test %d: prod B\n",count);
    buffer = 2;
    armci_msg_igop(&buffer, sizeof(int), prod);
    answer = 1 << size;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    /* TEST OF MAX OPERATOR */

    if (rank==0) printf("Test %d: max A\n",count);
    buffer = rank;
    armci_msg_igop(&buffer, sizeof(int), max);
    answer = (size-1);
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    if (rank==0) printf("Test %d: max B\n",count);
    buffer = -rank;
    armci_msg_igop(&buffer, sizeof(int), min);
    answer = 0;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    /* TEST OF MIN OPERATOR */

    if (rank==0) printf("Test %d: min A\n",count);
    buffer = rank;
    armci_msg_igop(&buffer, sizeof(int), min);
    answer = 0;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    if (rank==0) printf("Test %d: min B\n",count);
    buffer = -rank;
    armci_msg_igop(&buffer, sizeof(int), min);
    answer = -(size-1);
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    /* TEST OF ABSMAX OPERATOR */

    if (rank==0) printf("Test %d: absmax A\n",count);
    buffer = rank;
    armci_msg_igop(&buffer, sizeof(int), absmax);
    answer = (size-1);
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    if (rank==0) printf("Test %d: absmax B\n",count);
    buffer = -rank;
    armci_msg_igop(&buffer, sizeof(int), absmax);
    answer = -(size-1);
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    /* TEST OF ABSMIN OPERATOR */

    if (rank==0) printf("Test %d: absmin A\n",count);
    buffer = rank;
    armci_msg_igop(&buffer, sizeof(int), absmin);
    answer = 0;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    if (rank==0) printf("Test %d: absmin B\n",count);
    buffer = -rank;
    armci_msg_igop(&buffer, sizeof(int), absmin);
    answer = 0;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    /* TESTING OR OPERATOR */

    if (rank==0) printf("Test %d: or A\n",count);
    buffer = ( rank==0 ? 1 : 0 );
    armci_msg_igop(&buffer, sizeof(int), or);
    answer = 1;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    if (rank==0) printf("Test %d: or B\n",count);
    buffer = 0;
    armci_msg_igop(&buffer, sizeof(int), or);
    answer = 0;
    checkanswer(count,buffer,answer);
    count++;
    ARMCI_Barrier();

    /* END OF TESTS */

    ARMCI_Finalize();

    MPI_Finalize();

    return 0;
}
