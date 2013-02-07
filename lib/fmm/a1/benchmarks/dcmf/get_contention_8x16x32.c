/* 
 * The following is a notice of limited availability of the code, and disclaimer
 * which must be included in the prologue of the code and in all source listings
 * of the code.
 *
 * Copyright (c) 2010  Argonne Leadership Comgeting Facility, Argonne National
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

#include "bench.h"
#include "float.h"

void get_contention()
{

    unsigned int iter, size, dst;
    unsigned int i, j, k, s;
    unsigned int xdim, ydim, zdim;
    unsigned int xdisp, ydisp, zdisp;
    DCMF_Request_t get_req[ITERATIONS];
    DCMF_Callback_t get_done;
    unsigned int done_count;
    DCMF_NetworkCoord_t myaddr, dstaddr;
    DCMF_Network ntwk;
    char buf[50];

    get_done.function = done;
    get_done.clientdata = (void *) &done_count;

    DCMF_Messager_rank2network(nranks - 1, DCMF_TORUS_NETWORK, &dstaddr);
    xdim = dstaddr.torus.x + 1;
    ydim = dstaddr.torus.y + 1;
    zdim = dstaddr.torus.z + 1;

    if (myrank == 0)
    {
        printf("Dimensions of Torus : %d, %d, %d \n", xdim, ydim, zdim);
        fflush(stdout);
    }

    DCMF_Messager_rank2network(myrank, DCMF_TORUS_NETWORK, &myaddr);
    dstaddr.network = myaddr.network;
    dstaddr.torus.t = myaddr.torus.t;

    int size_array[] = { 8, 64, 512, 4096, 32768, 262144, 1048576 };
    int size_count = sizeof(size_array) / sizeof(int);

    int disp_array[][3] = { { 0, 0, 1 }, { 0, 0, 16 }, { 0, 8, 16 }, { 4, 8,
                                                                        16 },
                             { 0, 1, 16 }, { 1, 1, 16 }, { 0, 8, 16 }, { 1, 8,
                                                                         16 },
                             { 8, 8, 16 }, { 1, 16, 16 }, { 8, 16, 16 } };
    int disp_count = sizeof(disp_array) / (sizeof(int) * 3);

    for (s = 0; s < size_count; s++)
    {
        size = size_array[s];

        if (myrank == 0)
        {
            printf("Message Size : %20d \n", size);
            printf("%30s  %20s \n",
                   "Displacement b/w Pairs",
                   "Avg Bandwidth (Mbps)");
            fflush(stdout);
        }

        /*Assumes all dimensions are equal*/
        for (i = 0; i < disp_count; i++)
        {
            xdisp = disp_array[i][0];
            ydisp = disp_array[i][1];
            zdisp = disp_array[i][2];

            dstaddr.torus.x = (myaddr.torus.x + xdisp) % xdim;
            dstaddr.torus.y = (myaddr.torus.y + ydisp) % ydim;
            dstaddr.torus.z = (myaddr.torus.z + zdisp) % zdim;

            DCMF_Messager_network2rank(&dstaddr, &dst, &ntwk);

            barrier();

            /***********************
             * start timer          *
             ***********************/
            t_start = DCMF_Timebase();

            done_count = ITERATIONS;
            for (iter = 0; iter < ITERATIONS; iter++)
            {
                DCMF_Get(&get_reg,
                         &get_req[iter],
                         get_done,
                         DCMF_SEQUENTIAL_CONSISTENCY,
                         dst,
                         size,
                         memregion[myrank],
                         memregion[dst],
                         iter * size,
                         MAX_MSG_SIZE * ITERATIONS + iter * size);
            }
            while (done_count)
                DCMF_Messager_advance();

            t_stop = DCMF_Timebase();
            /***********************
             * stop timer          *
             ***********************/
            t_sec = (t_stop - t_start) / (clockMHz * 1000000);
            bw = (ITERATIONS * size) / (t_sec * 1024 * 1024);

            barrier();
            allreduce(-1,
                      (char *) &bw,
                      (char *) &bw_avg,
                      1,
                      DCMF_DOUBLE,
                      DCMF_SUM);

            if (myrank == 0)
            {
                bw_avg = bw_avg / nranks;
                sprintf(buf, "(%d)(%d)(%d)", xdisp, ydisp, zdisp);
                printf("%30s %20.0f \n", buf, bw_avg);
                fflush(stdout);
            }
        }

    }
}

int main()
{
    DCMF_Messager_initialize();

    init();

    barrier_init(DCMF_DEFAULT_GLOBALBARRIER_PROTOCOL);

    allreduce_init(DCMF_DEFAULT_GLOBALALLREDUCE_PROTOCOL);

    control_init(DCMF_DEFAULT_CONTROL_PROTOCOL, DCMF_DEFAULT_NETWORK);

    memregion_init(MAX_MSG_SIZE * ITERATIONS * 2);

    get_init(DCMF_DEFAULT_PUT_PROTOCOL, DCMF_TORUS_NETWORK);

    if (myrank == 0)
    {
        printf("Get Bandwidth - All processes communication in pairs \n");
        fflush(stdout);
    }
    get_contention();

    if (myrank == 0)
    {
        printf("Benchmark Complete \n");
        fflush(stdout);
    }

    memregion_finalize();

    DCMF_Messager_finalize();

    return 0;
}
