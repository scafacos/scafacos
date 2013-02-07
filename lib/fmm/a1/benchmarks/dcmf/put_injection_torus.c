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

void put_injection_noremotecallback()
{

    DCMF_Request_t put_req[ITERATIONS + SKIP];
    DCMF_Request_t put_flush;
    DCMF_Callback_t put_done, put_ack;
    DCMF_Callback_t put_flush_done, put_flush_ack;
    int flush_ack_count;
    int ack_count, done_count;
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
    put_ack.function = NULL;
    put_ack.clientdata = NULL;

    put_flush_done.function = NULL;
    put_flush_done.clientdata = NULL;
    put_flush_ack.function = done;
    put_flush_ack.clientdata = (void *) &flush_ack_count;

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer,
                "%20s  %30s  %20s",
                "Msg Size",
                "InjectionRt(MBps)-farthest pairs",
                "closest pairs");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    barrier();

    for (msgsize = 1; msgsize < MAX_MSG_SIZE; msgsize *= 2)
    {

        /***********************
         * warmup               *
         ***********************/
        done_count = SKIP;
        for (i = 0; i < SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     1 + i * msgsize,
                     1 + MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }
        while (done_count)
            DCMF_Messager_advance();

        flush_ack_count = 1;
        DCMF_Put(&put_reg,
                 &put_flush,
                 put_flush_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 dst,
                 1,
                 memregion[myrank],
                 memregion[dst],
                 0,
                 0,
                 put_flush_ack);
        while (flush_ack_count)
            DCMF_Messager_advance();

        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        done_count = ITERATIONS;
        for (i = SKIP; i < ITERATIONS + SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     1 + i * msgsize,
                     1 + MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }
        while (done_count)
            DCMF_Messager_advance();
        t_stop = DCMF_Timebase();
        t_sec = ((t_stop - t_start) / clockMHz) / 1000000;

        /***********************
         * stop timer          *
         ***********************/

        flush_ack_count = 1;
        DCMF_Put(&put_reg,
                 &put_flush,
                 put_flush_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 dst,
                 1,
                 memregion[myrank],
                 memregion[dst],
                 0,
                 0,
                 put_flush_ack);
        while (flush_ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_sec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            printf("%20d %30.2f ", msgsize, ((ITERATIONS) * msgsize) / (t_avg
                    * 1024 * 1024));
            fflush(stdout);
        }

        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        done_count = ITERATIONS;

        for (i = SKIP; i < ITERATIONS + SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     (myrank + 1) % nranks,
                     msgsize,
                     memregion[myrank],
                     memregion[(myrank + 1) % nranks],
                     1 + i * msgsize,
                     1 + MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }

        while (done_count)
            DCMF_Messager_advance();
        t_stop = DCMF_Timebase();
        t_sec = ((t_stop - t_start) / clockMHz) / 1000000;

        /***********************
         * stop timer          *
         ***********************/

        flush_ack_count = 1;
        DCMF_Put(&put_reg,
                 &put_flush,
                 put_flush_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 (myrank + 1) % nranks,
                 1,
                 memregion[myrank],
                 memregion[(myrank + 1) % nranks],
                 0,
                 0,
                 put_flush_ack);
        while (flush_ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_sec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            printf("%20.2f \n", ((ITERATIONS) * msgsize)
                    / (t_avg * 1024 * 1024));
            fflush(stdout);
        }

    }
}
void put_injection_nolocalcallback()
{

    DCMF_Request_t put_req[ITERATIONS + SKIP];
    DCMF_Request_t put_flush;
    DCMF_Callback_t put_done, put_ack;
    DCMF_Callback_t put_flush_done, put_flush_ack;
    int flush_done_count, flush_ack_count;
    int ack_count;
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

    put_done.function = NULL;
    put_done.clientdata = NULL;
    put_ack.function = done;
    put_ack.clientdata = (void *) &ack_count;

    put_flush_done.function = done;
    put_flush_done.clientdata = (void *) &flush_done_count;
    put_flush_ack.function = done;
    put_flush_ack.clientdata = (void *) &flush_ack_count;

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer,
                "%20s  %30s  %20s",
                "Msg Size",
                "InjectionRt(MBps)-farthest pairs",
                "closest pairs");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    barrier();

    for (msgsize = 1; msgsize < MAX_MSG_SIZE; msgsize *= 2)
    {

        /***********************
         * warmup               *
         ***********************/
        ack_count = SKIP;
        for (i = 0; i < SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     1 + i * msgsize,
                     1 + MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }

        flush_done_count = 1;
        flush_ack_count = 1;
        DCMF_Put(&put_reg,
                 &put_flush,
                 put_flush_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 dst,
                 1,
                 memregion[myrank],
                 memregion[dst],
                 0,
                 0,
                 put_flush_ack);
        while (flush_done_count)
            DCMF_Messager_advance();

        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        ack_count += ITERATIONS;
        for (i = SKIP; i < ITERATIONS + SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     1 + i * msgsize,
                     1 + MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }

        flush_done_count = 1;
        flush_ack_count += 1;
        DCMF_Put(&put_reg,
                 &put_flush,
                 put_flush_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 dst,
                 1,
                 memregion[myrank],
                 memregion[dst],
                 0,
                 0,
                 put_flush_ack);
        while (flush_done_count)
            DCMF_Messager_advance();
        t_stop = DCMF_Timebase();
        t_sec = ((t_stop - t_start) / clockMHz) / 1000000;

        /***********************
         * stop timer          *
         ***********************/

        while (ack_count || flush_ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_sec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            printf("%20d %30.2f ", msgsize, ((ITERATIONS) * msgsize) / (t_avg
                    * 1024 * 1024));
            fflush(stdout);
        }

        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        ack_count += ITERATIONS;
        for (i = SKIP; i < ITERATIONS + SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     (myrank + 1) % nranks,
                     msgsize,
                     memregion[myrank],
                     memregion[(myrank + 1) % nranks],
                     1 + i * msgsize,
                     1 + MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }

        flush_done_count = 1;
        flush_ack_count += 1;
        DCMF_Put(&put_reg,
                 &put_flush,
                 put_flush_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 (myrank + 1) % nranks,
                 1,
                 memregion[myrank],
                 memregion[(myrank + 1) % nranks],
                 0,
                 0,
                 put_flush_ack);
        while (flush_done_count)
            DCMF_Messager_advance();
        t_stop = DCMF_Timebase();
        t_sec = ((t_stop - t_start) / clockMHz) / 1000000;

        /***********************
         * stop timer          *
         ***********************/

        while (ack_count || flush_ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_sec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            printf("%20.2f \n", ((ITERATIONS) * msgsize)
                    / (t_avg * 1024 * 1024));
            fflush(stdout);
        }

    }
}

void put_injection_nocallback()
{

    DCMF_Request_t put_req[ITERATIONS + SKIP];
    DCMF_Request_t put_flush;
    DCMF_Callback_t put_done, put_ack;
    DCMF_Callback_t put_flush_done, put_flush_ack;
    int flush_done_count, flush_ack_count;
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

    put_done.function = NULL;
    put_done.clientdata = NULL;
    put_ack.function = NULL;
    put_ack.clientdata = NULL;

    put_flush_done.function = done;
    put_flush_done.clientdata = (void *) &flush_done_count;
    put_flush_ack.function = done;
    put_flush_ack.clientdata = (void *) &flush_ack_count;

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer,
                "%20s  %30s  %20s",
                "Msg Size",
                "InjectionRt(MBps)-farthest pairs",
                "closest pairs");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    barrier();

    for (msgsize = 1; msgsize < MAX_MSG_SIZE; msgsize *= 2)
    {

        /***********************
         * warmup               *
         ***********************/
        for (i = 0; i < SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     1 + i * msgsize,
                     1 + MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }

        flush_done_count = 1;
        flush_ack_count = 1;
        DCMF_Put(&put_reg,
                 &put_flush,
                 put_flush_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 dst,
                 1,
                 memregion[myrank],
                 memregion[dst],
                 0,
                 0,
                 put_flush_ack);
        while (flush_done_count)
            DCMF_Messager_advance();

        while (flush_ack_count)
            DCMF_Messager_advance();
        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        for (i = SKIP; i < ITERATIONS + SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     dst,
                     msgsize,
                     memregion[myrank],
                     memregion[dst],
                     1 + i * msgsize,
                     1 + MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }

        flush_done_count = 1;
        flush_ack_count += 1;
        DCMF_Put(&put_reg,
                 &put_flush,
                 put_flush_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 dst,
                 1,
                 memregion[myrank],
                 memregion[dst],
                 0,
                 0,
                 put_flush_ack);
        while (flush_done_count)
            DCMF_Messager_advance();
        t_stop = DCMF_Timebase();
        t_sec = ((t_stop - t_start) / clockMHz) / 1000000;

        /***********************
         * stop timer          *
         ***********************/

        while (flush_ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_sec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            printf("%20d %30.2f ", msgsize, ((ITERATIONS) * msgsize) / (t_avg
                    * 1024 * 1024));
            fflush(stdout);
        }

        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        for (i = SKIP; i < ITERATIONS + SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     (myrank + 1) % nranks,
                     msgsize,
                     memregion[myrank],
                     memregion[(myrank + 1) % nranks],
                     1 + i * msgsize,
                     1 + MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }

        flush_done_count = 1;
        flush_ack_count += 1;
        DCMF_Put(&put_reg,
                 &put_flush,
                 put_flush_done,
                 DCMF_SEQUENTIAL_CONSISTENCY,
                 (myrank + 1) % nranks,
                 1,
                 memregion[myrank],
                 memregion[(myrank + 1) % nranks],
                 0,
                 0,
                 put_flush_ack);
        while (flush_done_count)
            DCMF_Messager_advance();
        t_stop = DCMF_Timebase();
        t_sec = ((t_stop - t_start) / clockMHz) / 1000000;

        /***********************
         * stop timer          *
         ***********************/

        while (flush_ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_sec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            printf("%20.2f \n", ((ITERATIONS) * msgsize)
                    / (t_avg * 1024 * 1024));
            fflush(stdout);
        }

    }
}

void put_injection_callback()
{

    DCMF_Request_t put_req[ITERATIONS + SKIP];
    DCMF_Callback_t put_done, put_ack;
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

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer,
                "%20s  %30s  %20s",
                "Msg Size",
                "InjectionRt(MBps)-farthest pairs",
                "closest pairs");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    barrier();

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
                     &put_req[i],
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
        while (done_count)
            DCMF_Messager_advance();

        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        done_count = ITERATIONS;
        ack_count += ITERATIONS;

        for (i = SKIP; i < ITERATIONS + SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
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

        while (done_count)
            DCMF_Messager_advance();
        t_stop = DCMF_Timebase();
        t_sec = ((t_stop - t_start) / clockMHz) / 1000000;

        /***********************
         * stop timer          *
         ***********************/

        while (ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_sec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            printf("%20d %30.2f ", msgsize, ((ITERATIONS) * msgsize) / (t_avg
                    * 1024 * 1024));
            fflush(stdout);
        }

        /***********************
         * start timer          *
         ***********************/

        t_start = DCMF_Timebase();

        done_count = ITERATIONS;
        ack_count += ITERATIONS;

        for (i = SKIP; i < ITERATIONS + SKIP; i++)
        {
            DCMF_Put(&put_reg,
                     &put_req[i],
                     put_done,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     (myrank + 1) % nranks,
                     msgsize,
                     memregion[myrank],
                     memregion[(myrank + 1) % nranks],
                     i * msgsize,
                     MAX_BUF_SIZE + i * msgsize,
                     put_ack);
        }

        while (done_count)
            DCMF_Messager_advance();
        t_stop = DCMF_Timebase();
        t_sec = ((t_stop - t_start) / clockMHz) / 1000000;

        /***********************
         * stop timer          *
         ***********************/

        while (ack_count)
            DCMF_Messager_advance();

        barrier();
        allreduce(-1,
                  (char *) &t_sec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);

        if (myrank == 0)
        {
            t_avg = t_avg / nranks;
            printf("%20.2f \n", ((ITERATIONS) * msgsize)
                    / (t_avg * 1024 * 1024));
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

    if (myrank == 0)
    {
        printf("Put Injection Rate (MBps) with static routing and callback \n");
        fflush(stdout);
    }
    put_injection_callback();

    if (myrank == 0)
    {
        printf("Put Injection Rate (MBps) with static routing and no local callback \n");
        fflush(stdout);
    }
    put_injection_nolocalcallback();

    if (myrank == 0)
    {
        printf("Put Injection Rate (MBps) with static routing and no remote callback \n");
        fflush(stdout);
    }
    put_injection_noremotecallback();

    if (myrank == 0)
    {
        printf("Put Injection Rate (MBps) with static routing and nocallback \n");
        fflush(stdout);
    }
    put_injection_nocallback();

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
