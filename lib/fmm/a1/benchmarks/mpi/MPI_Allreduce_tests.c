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

int main(int argc, char **argv)
{
    int provided;

    int rank, size;

    int count;

    int in, out, ans;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    /* TESTING SUM OPERATOR */

    if (rank==0) printf("Test %d: sum 1 out-of-place\n",count++);
    in  = ( rank==0 ? 1 : 0 );
    out = 0;
    ans = 1;
    MPI_Allreduce(&in,&out,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (out==ans) printf("%4d: in=%d out=%d ans=%d PASSED \n",rank,in,out,ans);
    else          printf("%4d: in=%d out=%d ans=%d FAILED \n",rank,in,out,ans);
    fflush(stdout);

    if (rank==0) printf("Test %d: sum 1 in-place\n",count++);
    in  = ( rank==0 ? 1 : 0 );
    out = in;
    ans = 1;
    MPI_Allreduce(MPI_IN_PLACE,&out,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (out==ans) printf("%4d: in=%d out=%d ans=%d PASSED \n",rank,in,out,ans);
    else          printf("%4d: in=%d out=%d ans=%d FAILED \n",rank,in,out,ans);
    fflush(stdout);

    if (rank==0) printf("Test %d: sum 2 out-of-place\n",count++);
    in  = 1;
    out = 0;
    ans = size;
    MPI_Allreduce(&in,&out,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (out==ans) printf("%4d: in=%d out=%d ans=%d PASSED \n",rank,in,out,ans);
    else          printf("%4d: in=%d out=%d ans=%d FAILED \n",rank,in,out,ans);
    fflush(stdout);

    if (rank==0) printf("Test %d: sum 2 in-place\n",count++);
    in  = 1;
    out = in;
    ans = size;
    MPI_Allreduce(MPI_IN_PLACE,&out,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (out==ans) printf("%4d: in=%d out=%d ans=%d PASSED \n",rank,in,out,ans);
    else          printf("%4d: in=%d out=%d ans=%d FAILED \n",rank,in,out,ans);
    fflush(stdout);

    /* TESTING PROD OPERATOR */

    if (rank==0) printf("Test %d: prod 1 out-of-place\n",count++);
    in  = ( rank==0 ? 1 : 0 );
    out = 0;
    ans = 0;
    MPI_Allreduce(&in,&out,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
    if (out==ans) printf("%4d: in=%d out=%d ans=%d PASSED \n",rank,in,out,ans);
    else          printf("%4d: in=%d out=%d ans=%d FAILED \n",rank,in,out,ans);
    fflush(stdout);

    if (rank==0) printf("Test %d: prod 1 in-place\n",count++);
    in  = ( rank==0 ? 1 : 0 );
    out = in;
    ans = 0;
    MPI_Allreduce(MPI_IN_PLACE,&out,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
    if (out==ans) printf("%4d: in=%d out=%d ans=%d PASSED \n",rank,in,out,ans);
    else          printf("%4d: in=%d out=%d ans=%d FAILED \n",rank,in,out,ans);;
    fflush(stdout);

    if (rank==0) printf("Test %d: prod 2 out-of-place\n",count++);
    in  = 1;
    out = 0;
    ans = 1;
    MPI_Allreduce(&in,&out,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
    if (out==ans) printf("%4d: in=%d out=%d ans=%d PASSED \n",rank,in,out,ans);
    else          printf("%4d: in=%d out=%d ans=%d FAILED \n",rank,in,out,ans);
    fflush(stdout);

    if (rank==0) printf("Test %d: prod 2 in-place\n",count++);
    in  = 1;
    out = in;
    ans = 1;
    MPI_Allreduce(MPI_IN_PLACE,&out,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
    if (out==ans) printf("%4d: in=%d out=%d ans=%d PASSED \n",rank,in,out,ans);
    else          printf("%4d: in=%d out=%d ans=%d FAILED \n",rank,in,out,ans);
    fflush(stdout);

    /* END OF TESTS */

    MPI_Finalize();

    return 0;
}
