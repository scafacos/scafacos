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

void put_static_localcallback()
{

        DCMF_Request_t put_req;
        DCMF_Callback_t put_done, put_ack, no_callback;
        int done_count, ack_count;
        int msgsize, i;

        put_done.function = done;
        put_done.clientdata = (void *) &done_count;
        no_callback.function = NULL;
        no_callback.clientdata = NULL;
        put_ack.function = done;
        put_ack.clientdata = (void *) &ack_count;

        char buffer[50];
        sprintf(buffer, "%20s  %20s", "Msg Size", "Latency (usec)");
        printf("%s \n", buffer);
        fflush(stdout);

        barrier();

        if(myrank == 0)
        {

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
                         &put_req,
                         put_done,
                         DCMF_SEQUENTIAL_CONSISTENCY,
                         1,
                         msgsize,
                         memregion[myrank],
                         memregion[myrank + 1],
                         i * msgsize,
                         i * msgsize,
                         no_callback);

                while(done_count)
                      DCMF_Messager_advance();
            }

            t_stop = DCMF_Timebase();

            /***********************
             * stop timer          *
             ***********************/

            t_usec = ((t_stop - t_start) / clockMHz);
            printf("%20d %20.2f \n", msgsize, t_usec / (ITERATIONS));

            ack_count = 1;

            DCMF_Put(&put_reg,
                     &put_req,
                     no_callback,
                     DCMF_SEQUENTIAL_CONSISTENCY,
                     1,
                     msgsize,
                     memregion[myrank],
                     memregion[myrank + 1],
                     i * msgsize,
                     i * msgsize,
                     put_ack);

            while(ack_count)
                  DCMF_Messager_advance();

        }

    }
    barrier();
}

void put_static_remotecallback()
{

        DCMF_Request_t put_req;
        DCMF_Callback_t put_done, put_ack;
        int done_count, ack_count;
        int msgsize, i;

        put_done.function = done;
        put_done.clientdata = (void *) &done_count;
        put_ack.function = done;
        put_ack.clientdata = (void *) &ack_count;

        char buffer[50];
        sprintf(buffer, "%20s  %20s", "Msg Size", "Latency (usec)");
        printf("%s \n", buffer);
        fflush(stdout);

        barrier();

        if(myrank == 0)
        {

        for (msgsize = 1; msgsize < MAX_MSG_SIZE; msgsize *= 2)
        {

            /***********************
             * start timer          *
             ***********************/

            t_start = DCMF_Timebase();

            for (i = 0; i < ITERATIONS; i++)
            {
                done_count = 1; 
                ack_count = 1; 

                DCMF_Put(&put_reg,
                         &put_req,
                         put_done,
                         DCMF_SEQUENTIAL_CONSISTENCY,
                         1,
                         msgsize,
                         memregion[myrank],
                         memregion[myrank + 1],
                         i * msgsize,
                         i * msgsize,
                         put_ack);

                while(done_count || ack_count)
                      DCMF_Messager_advance();
            }

            t_stop = DCMF_Timebase();

            /***********************
             * stop timer          *
             ***********************/

            t_usec = ((t_stop - t_start) / clockMHz);
            printf("%20d %20.2f \n", msgsize, t_usec / (ITERATIONS));

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

    memregion_init(MAX_BUF_SIZE * nranks);

    put_init(DCMF_DEFAULT_PUT_PROTOCOL, DCMF_TORUS_NETWORK);

    barrier();

    if (myrank == 0)
    {
        printf("Put Latency (usec) with local callback \n");
        fflush(stdout);
    }
    put_static_localcallback();

    if (myrank == 0)
    {
        printf("Put Latency (usec) with remote callback \n");
        fflush(stdout);
    }
    put_static_remotecallback();

    barrier();

    printf("[%d] Benchmark complete\n", myrank);
    fflush(stdout);

    memregion_finalize();

    DCMF_Messager_finalize();

    return 0;
}
