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

#include "bench.h"

void put_localvsremote()
{

    DCMF_Request_t skip_req[SKIP];
    DCMF_Request_t local_req1[ITERATIONS], local_req2[ITERATIONS];
    DCMF_Request_t flush_req1, flush_req2;
    DCMF_Request_t remote_req1[ITERATIONS], remote_req2[ITERATIONS];
    DCMF_Request_t req1[ITERATIONS], req2[ITERATIONS];
    DCMF_Callback_t put_done, put_ack, put_nocallback;
    int done_count, ack_count;
    unsigned int msgsize, i, dst;
    DCMF_NetworkCoord_t myaddr, dstaddr;
    DCMF_Network ntwk;

    DCMF_Messager_rank2network(myrank, DCMF_TORUS_NETWORK, &myaddr);

    dstaddr.network = myaddr.network;
    dstaddr.torus.x = (myaddr.torus.x + 3) % 8;
    dstaddr.torus.y = (myaddr.torus.y + 3) % 8;
    dstaddr.torus.z = (myaddr.torus.z + 3) % 8;
    dstaddr.torus.t = myaddr.torus.t;

    DCMF_Messager_network2rank(&dstaddr, &dst, &ntwk);

    put_done.function = done;
    put_done.clientdata = (void *) &done_count;
    put_ack.function = done;
    put_ack.clientdata = (void *) &ack_count;
    put_nocallback.function = NULL;
    put_nocallback.clientdata = NULL;

    if (myrank == 0)
    {
        printf("Put call overhead in usec\n");
        fflush(stdout);
    }

    if (myrank == 0)
    {
        char buffer[100];
        sprintf(buffer,
                "%20s  %20s %20s",
                "Msg Size",
                "Farthest pairs",
                "Closest pairs");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    for (msgsize = 1; msgsize < MAX_MSG_SIZE; msgsize *= 2)
    {

        /***********************
         * warmup               *
         ***********************/
        done_count = SKIP;
        ack_count = SKIP;
        for (i = 0; i < SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &skip_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }
        while (done_count || ack_count)
            DCMF_Messager_advance();

        t_avg = 0;
        t_avg1 = 0;
        barrier();

        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        for (i = 0; i < ITERATIONS; i++)
        {
            DCMF_Put(&put_reg,
                     &local_req1[i],
                     put_nocallback,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_nocallback);
        }

        t_stop = DCMF_Timebase();
        t_usec = (t_stop - t_start) / (clockMHz * ITERATIONS);

        done_count = 1;
        ack_count = 1;
        DCMF_Put(&put_reg,
                 &flush_req1,
                 put_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 dst,
                 1,
                 memregion[myrank],
                 memregion[dst],
                 0,
                 0,
                 put_ack);
        while (done_count || ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_usec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

        t_start = DCMF_Timebase();

        for (i = 0; i < ITERATIONS; i++)
        {
            DCMF_Put(&put_reg,
                     &local_req2[i],
                     put_nocallback,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     (myrank + 1) % nranks,
                     msgsize,
                     memregion[myrank],
                     memregion[(myrank + 1) % nranks],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_nocallback);
        }

        t_stop = DCMF_Timebase();
        t_usec1 = (t_stop - t_start) / (clockMHz * ITERATIONS);

        done_count = 1;
        ack_count = 1;
        DCMF_Put(&put_reg,
                 &flush_req2,
                 put_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 (myrank + 1) % nranks,
                 1,
                 memregion[myrank],
                 memregion[(myrank + 1) % nranks],
                 0,
                 0,
                 put_ack);
        while (done_count || ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_usec1,
                  (char *) &t_avg1,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            t_avg1 = t_avg1 / nranks;
            printf("%20d %20.2f %20.2f \n", msgsize, t_avg, t_avg1);
            fflush(stdout);
        }

    }

    if (myrank == 0)
    {
        printf("Put latency in usec with local vs remot completion \n");
        fflush(stdout);
    }

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer,
                "%20s  %20s  %20s  %20s  %20s %20s  %20s",
                "Msg Size",
                "Farthest pairs-local",
                "Farthest pairs-remote",
                "Farthest pairs-both",
                "Closest pairs-local",
                "Closest pairs-remote",
                "Closest pairs-both");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    barrier();

    for (msgsize = 1; msgsize < MAX_MSG_SIZE; msgsize *= 2)
    {

        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        for (i = 0; i < ITERATIONS; i++)
        {
            done_count = 1;
            DCMF_Put(&put_reg,
                     &local_req1[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_nocallback);
            while (done_count)
                DCMF_Messager_advance();
        }

        t_stop = DCMF_Timebase();
        t_usec = (t_stop - t_start) / (clockMHz * ITERATIONS);

        done_count = 1;
        ack_count = 1;
        DCMF_Put(&put_reg,
                 &flush_req1,
                 put_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 dst,
                 1,
                 memregion[myrank],
                 memregion[dst],
                 0,
                 0,
                 put_ack);
        while (done_count || ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_usec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

        t_start = DCMF_Timebase();

        for (i = 0; i < ITERATIONS; i++)
        {
            ack_count = 1;
            DCMF_Put(&put_reg,
                     &remote_req1[i],
                     put_nocallback,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_ack);
            while (ack_count)
                DCMF_Messager_advance();
        }

        t_stop = DCMF_Timebase();
        t_usec1 = (t_stop - t_start) / (clockMHz * ITERATIONS);

        barrier();
        allreduce(-1,
                  (char *) &t_usec1,
                  (char *) &t_avg1,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

        t_start = DCMF_Timebase();

        for (i = 0; i < ITERATIONS; i++)
        {
            done_count = 1;
            ack_count = 1;
            DCMF_Put(&put_reg,
                     &req1[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_ack);
            while (done_count || ack_count)
                DCMF_Messager_advance();
        }

        t_stop = DCMF_Timebase();
        t_usec2 = (t_stop - t_start) / (clockMHz * ITERATIONS);

        /***********************
         * stop timer          *
         ***********************/

        barrier();
        allreduce(-1,
                  (char *) &t_usec2,
                  (char *) &t_avg2,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            t_avg1 = t_avg1 / nranks;
            t_avg2 = t_avg2 / nranks;
            printf("%20d %20.2f %20.2f %20.2f", msgsize, t_avg, t_avg1, t_avg2);
            fflush(stdout);
        }

        t_avg = 0;
        t_avg1 = 0, t_avg2 = 0;
        barrier();
        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        for (i = 0; i < ITERATIONS; i++)
        {
            done_count = 1;
            DCMF_Put(&put_reg,
                     &local_req2[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     (myrank + 1) % nranks,
                     msgsize,
                     memregion[myrank],
                     memregion[(myrank + 1) % nranks],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_nocallback);
            while (done_count)
                DCMF_Messager_advance();
        }

        t_stop = DCMF_Timebase();
        t_usec = (t_stop - t_start) / (clockMHz * ITERATIONS);

        done_count = 1;
        ack_count = 1;
        DCMF_Put(&put_reg,
                 &flush_req2,
                 put_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 (myrank + 1) % nranks,
                 1,
                 memregion[myrank],
                 memregion[(myrank + 1) % nranks],
                 0,
                 0,
                 put_ack);
        while (done_count || ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_usec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

        t_start = DCMF_Timebase();

        for (i = 0; i < ITERATIONS; i++)
        {
            ack_count = 1;
            DCMF_Put(&put_reg,
                     &remote_req2[i],
                     put_nocallback,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     (myrank + 1) % nranks,
                     msgsize,
                     memregion[myrank],
                     memregion[(myrank + 1) % nranks],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_ack);
            while (ack_count)
                DCMF_Messager_advance();
        }

        t_stop = DCMF_Timebase();
        t_usec1 = (t_stop - t_start) / (clockMHz * ITERATIONS);

        barrier();
        allreduce(-1,
                  (char *) &t_usec1,
                  (char *) &t_avg1,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

        t_start = DCMF_Timebase();

        for (i = 0; i < ITERATIONS; i++)
        {
            done_count = 1;
            ack_count = 1;
            DCMF_Put(&put_reg,
                     &req2[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     (myrank + 1) % nranks,
                     msgsize,
                     memregion[myrank],
                     memregion[(myrank + 1) % nranks],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_ack);
            while (done_count || ack_count)
                DCMF_Messager_advance();
        }

        t_stop = DCMF_Timebase();
        t_usec2 = (t_stop - t_start) / (clockMHz * ITERATIONS);

        /***********************
         * stop timer          *
         ***********************/

        barrier();
        allreduce(-1,
                  (char *) &t_usec2,
                  (char *) &t_avg2,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            t_avg1 = t_avg1 / nranks;
            t_avg2 = t_avg2 / nranks;
            printf("%20.2f %20.2f  %20.2f \n", t_avg, t_avg1, t_avg2);
            fflush(stdout);
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

    memregion_init(1 + MAX_BUF_SIZE * 2);

    put_init(DCMF_DEFAULT_PUT_PROTOCOL, DCMF_TORUS_NETWORK);

    barrier();

    put_localvsremote();

    barrier();

    if (myrank == 0)
    {
        printf("[%d] Benchmark Complete \n", myrank);
        fflush(stdout);
    }

    memregion_finalize();

    DCMF_Messager_finalize();

    return 0;
}
