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

#define ITERATIONS_LOCAL 10000
#define MAX_MSG_SIZE_LOCAL 1*1024*1024

void put_remoteadvance()
{

    DCMF_Request_t *put_req;
    DCMF_Callback_t put_done, put_ack;
    int done_count, ack_count;
    unsigned int msgsize, i, dst;

    put_req = (DCMF_Request_t *) malloc(sizeof(DCMF_Request_t)
            * ITERATIONS_LOCAL);

    put_done.function = done;
    put_done.clientdata = (void *) &done_count;
    put_ack.function = done;
    put_ack.clientdata = (void *) &ack_count;

    if (myrank == 0)
    {
        printf("Put latency in usec\n");
        fflush(stdout);
    }

    if (myrank == 0)
    {
        char buffer[100];
        sprintf(buffer,
                "%20s  %20s %20s",
                "Msg Size",
                "Put-Remote Barrier",
                "Put-Remote Sleep");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    if (myrank == 0)
    {

        for (msgsize = 1; msgsize < MAX_MSG_SIZE_LOCAL; msgsize *= 2)
        {

            /***********************
             * start timer          *
             ***********************/

            t_start = DCMF_Timebase();

            done_count = 10000;
            ack_count = 10000;

            for (i = 0; i < ITERATIONS_LOCAL; i++)
            {
                DCMF_Put(&put_reg,
                         &put_req[i],
                         put_done,
                         DCMF_SEQUENTIAL_CONSISTENCY,
                         (myrank + 1) % nranks,
                         msgsize,
                         memregion[0],
                         memregion[1],
                         0,
                         0,
                         put_ack);
            }

            while (done_count > 0 || ack_count > 0)
                DCMF_Messager_advance();

            t_stop = DCMF_Timebase();
            t_usec = (t_stop - t_start) / (clockMHz * ITERATIONS_LOCAL);

            /***********************
             * stop timer          *
             ***********************/

            if (myrank == 0)
            {
                printf("%20d %20.2f ", msgsize, t_usec);
                fflush(stdout);
            }

            barrier();

            /***********************
             * start timer          *
             ***********************/

            t_start = DCMF_Timebase();

            done_count = 10000;
            ack_count = 10000;

            for (i = 0; i < ITERATIONS_LOCAL; i++)
            {
                DCMF_Put(&put_reg,
                         &put_req[i],
                         put_done,
                         DCMF_SEQUENTIAL_CONSISTENCY,
                         (myrank + 1) % nranks,
                         msgsize,
                         memregion[0],
                         memregion[1],
                         0,
                         0,
                         put_ack);
            }

            while (done_count > 0 || ack_count > 0)
                DCMF_Messager_advance();

            t_stop = DCMF_Timebase();
            t_usec = (t_stop - t_start) / (clockMHz * ITERATIONS_LOCAL);

            /***********************
             * stop timer          *
             ***********************/

            if (myrank == 0)
            {
                printf("%20.2f \n", t_usec);
                fflush(stdout);
            }

            barrier();

        }

    }
    else
    {

        for (msgsize = 1; msgsize < MAX_MSG_SIZE_LOCAL; msgsize *= 2)
        {

            barrier();

            DCMF_CriticalSection_enter(0);

            sleep(10);

            DCMF_CriticalSection_exit(0);

            barrier();
        }

    }

    barrier();
}

int main()
{
    DCMF_Messager_initialize();

    init();

    barrier_init(DCMF_DEFAULT_GLOBALBARRIER_PROTOCOL);

    control_init(DCMF_DEFAULT_CONTROL_PROTOCOL, DCMF_DEFAULT_NETWORK);

    memregion_init(MAX_MSG_SIZE_LOCAL);

    put_init(DCMF_DEFAULT_PUT_PROTOCOL, DCMF_TORUS_NETWORK);

    barrier();

    put_remoteadvance();

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
