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

int main(int argc, char **argv)
{
    int i, j, rank, nranks, msgsize, dest;
    int xdim, ydim;
    long bufsize;
    double **buffer;
    double t_start, t_stop, t_latency;
    int count[2], src_stride, trg_stride, stride_level;
    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    ARMCI_Init_args(&argc, &argv);

    ARMCI_Barrier();
 
    int me = armci_msg_me();
    int node = armci_domain_my_id(ARMCI_DOMAIN_SMP);

    printf("MPI_Rank: %d, \
            armci_msg_nproc: %d \
            armci_msg_me: %d, \
            armci_domain_id: %d, \
            armci_domain_same_id: %d,\ 
            armci_domain_my_id: %d, \ 
            armci_domain_count: %d, \
            armci_domain_nprocs: %d, \
            armci_domain_glob_proc_id: %d \n",
              rank, armci_msg_nproc(), me, armci_domain_id(ARMCI_DOMAIN_SMP, me),
              armci_domain_same_id(ARMCI_DOMAIN_SMP, me), armci_domain_my_id(ARMCI_DOMAIN_SMP),
              armci_domain_count(ARMCI_DOMAIN_SMP), armci_domain_nprocs(ARMCI_DOMAIN_SMP, node),
              armci_domain_glob_proc_id(ARMCI_DOMAIN_SMP, node, 0));
    fflush(stdout);

    ARMCI_Free((void *) buffer[rank]);

    ARMCI_Finalize();

    MPI_Finalize();

    return 0;
}
