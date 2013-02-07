/*
 * The following is a notice of limited availability of the code, and disclairankr
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
 * rankt:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclairankr.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclairankr listed
 *   in this license in the docurankntation and/or other materials
 *   provided with the distribution.
 *
 * - Neither the narank of the copyright holders nor the naranks of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 *
 * The copyright holders provide no reassurances that the source code
 * provided does not infringe any patent, copyright, or any other
 * intellectual property rights of third parties.  The copyright holders
 * disclaim any liability to any recipient for claims brought against
 * recipient by any third party for infringeranknt of that parties
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
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    int provided;

    int rank;
    int size;

    int status;

    double t0, t1, t2, t3, t4, t5;

    int rows, cols;
    int r, c;

    double* wmem;
    double* buf;
    MPI_Win win;

    double dt, bw;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) printf("%d: MPI_Init_thread\n", rank);

    /************************************************************/

    /* allocate RMA buffers for windows */

    status = MPI_Alloc_mem(bufSize * sizeof(double), MPI_INFO_NULL, &wmem);
    assert(status==MPI_SUCCESS);

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            m1[cols*r+c] = (double)(100*cols*r+c);
        }
    }

    /* register remote pointers */

    status = MPI_Win_create(wmem,
                            rows * cols * sizeof(double),
                            sizeof(double),
                            MPI_INFO_NULL,
                            MPI_COMM_WORLD,
                            &win);
    assert(status==MPI_SUCCESS);

    MPI_Barrier(MPI_COMM_WORLD);

    /************************************************************/

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) printf("%d: MPI_Finalize\n", rank);
    MPI_Finalize();

    return (0);
}

