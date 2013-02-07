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

void flush_put()
{

    DCMF_Request_t put_req[nranks];
    DCMF_Callback_t put_done, put_ack;
    int ack_count;
    int dest, i;

    put_done.function = NULL;
    put_done.clientdata = NULL;
    put_ack.function = done;
    put_ack.clientdata = (void *) &ack_count;

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer,
                "%20s %30s",
                "Flush Latency (us)",
                "Flush Restart Latency(us)");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    barrier();

    /***********************
     * warmup               *
     ***********************/
    for (i = 0; i < SKIP; i++)
    {

        ack_count = nranks - 1;
        for (dest = 0; dest < nranks; dest++)
        {
            if (dest != myrank)
            {
                DCMF_Put(&put_reg,
                         &put_req[dest],
                         put_done,
                         DCMF_SEQUENTIAL_CONSISTENCY,
                         dest,
                         1,
                         memregion[myrank],
                         memregion[dest],
                         0,
                         1,
                         put_ack);

            }
        }
        while (ack_count)
            DCMF_Messager_advance();

    }

    /***********************
     * start timer          *
     ***********************/
    t_start = DCMF_Timebase();

    for (i = SKIP; i < ITERATIONS + SKIP; i++)
    {

        ack_count = nranks - 1;
        for (dest = 0; dest < nranks; dest++)
        {
            if (dest != myrank)
            {
                DCMF_Put(&put_reg,
                         &put_req[dest],
                         put_done,
                         DCMF_SEQUENTIAL_CONSISTENCY,
                         dest,
                         1,
                         memregion[myrank],
                         memregion[dest],
                         0,
                         1,
                         put_ack);

            }
        }
        while (ack_count)
            DCMF_Messager_advance();

    }

    t_stop = DCMF_Timebase();
    /***********************
     * stop timer          *
     ***********************/
    t_usec = ((t_stop - t_start) / clockMHz);
    t_usec = t_usec / (ITERATIONS);

    barrier();

    allreduce(-1, (char *) &t_usec, (char *) &t_max, 1, DCMF_DOUBLE, DCMF_MAX);

    barrier();

    if (myrank == 0)
    {
        printf("%20.0f ", t_max);
        fflush(stdout);
    }

    /***********************
     * start timer          *
     ***********************/
    t_start = DCMF_Timebase();

    for (i = SKIP; i < ITERATIONS + SKIP; i++)
    {

        ack_count = nranks - 1;
        for (dest = 0; dest < nranks; dest++)
        {
            if (dest != myrank)
            {
                DCMF_Restart(&put_req[dest]);
            }
        }
        while (ack_count)
            DCMF_Messager_advance();

    }

    t_stop = DCMF_Timebase();
    /***********************
     * stop timer          *
     ***********************/
    t_usec = ((t_stop - t_start) / clockMHz);
    t_usec = t_usec / ITERATIONS;

    barrier();

    allreduce(-1, (char *) &t_usec, (char *) &t_max, 1, DCMF_DOUBLE, DCMF_MAX);

    barrier();

    if (myrank == 0)
    {
        printf("%20.0f \n", t_max);
        fflush(stdout);
    }

}

void flush_multicast()
{

    int i;
    barrier();

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer, "%20s", "Flush Latency (us)");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    for (i = 0; i < SKIP; i++)
    {

        mc_active = 1;
        mc_rcv_active += (nranks - 1);
        DCMF_Multicast(&mc_info);
        while (mc_active || mc_rcv_active)
            DCMF_Messager_advance();

        mc_active = 1;
        mc_rcv_active += (nranks - 1);
        DCMF_Multicast(&mc_info);
        while (mc_active || mc_rcv_active)
            DCMF_Messager_advance();

    }

    t_start = DCMF_Timebase();
    /***********************
     * start timer          *
     ***********************/
    for (i = SKIP; i < ITERATIONS + SKIP; i++)
    {

        mc_active = 1;
        mc_rcv_active += (nranks - 1);
        DCMF_Multicast(&mc_info);
        while (mc_active || mc_rcv_active)
            DCMF_Messager_advance();

        mc_active = 1;
        mc_rcv_active += (nranks - 1);
        DCMF_Multicast(&mc_info);
        while (mc_active || mc_rcv_active)
            DCMF_Messager_advance();

    }
    t_stop = DCMF_Timebase();
    /***********************
     * stop timer          *
     ***********************/

    t_usec = ((t_stop - t_start) / clockMHz);
    t_usec = t_usec / (ITERATIONS);
    fflush(stdout);

    barrier();

    allreduce(-1, (char *) &t_usec, (char *) &t_max, 1, DCMF_DOUBLE, DCMF_MAX);

    barrier();

    if (myrank == 0)
    {
        printf("%20.0f \n", t_max);
        fflush(stdout);
    }

}

void flush_send(unsigned int size)
{

    barrier();

    int i, j;
    DCMF_Request_t snd_req[nranks];
    DCMF_Callback_t snd_callback;
    char snd_buffer[size];
    unsigned int snd_active;

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer, "%20s", "Flush Latency (us)");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    snd_callback.function = done;
    snd_callback.clientdata = (void *) &snd_active;

    for (i = 0; i < SKIP; i++)
    {
        snd_active = nranks - 1;
        snd_rcv_active += nranks - 1;
        ack_rcv_active += nranks - 1;

        for (j = 0; j < nranks; j++)
        {
            if (j != myrank)
            {
                DCMF_Send(&flush_snd_reg,
                          &snd_req[j],
                          snd_callback,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          j,
                          size,
                          snd_buffer,
                          snd_msginfo,
                          1);
            }
        }

        while (snd_active || snd_rcv_active || ack_rcv_active)
            DCMF_Messager_advance();
    }

    /***********************
     * start timer          *
     ***********************/
    t_start = DCMF_Timebase();
    for (i = 0; i < ITERATIONS; i++)
    {
        snd_active = nranks - 1;
        snd_rcv_active += nranks - 1;
        ack_rcv_active += nranks - 1;

        for (j = 0; j < nranks; j++)
        {
            if (j != myrank)
            {
                DCMF_Send(&flush_snd_reg,
                          &snd_req[j],
                          snd_callback,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          j,
                          size,
                          snd_buffer,
                          snd_msginfo,
                          1);
            }
        }

        while (snd_active || snd_rcv_active || ack_rcv_active)
            DCMF_Messager_advance();
    }
    t_stop = DCMF_Timebase();
    /***********************
     * stop timer          *
     ***********************/
    t_usec = ((t_stop - t_start) / clockMHz);
    t_usec = t_usec / ITERATIONS;
    fflush(stdout);

    barrier();

    allreduce(-1, (char *) &t_usec, (char *) &t_max, 1, DCMF_DOUBLE, DCMF_MAX);

    barrier();

    if (myrank == 0)
    {
        printf("%20.0f \n", t_max);
        fflush(stdout);
    }
}

int main()
{
    DCMF_Messager_initialize();

    init();

    barrier_init(DCMF_DEFAULT_GLOBALBARRIER_PROTOCOL);

    allreduce_init(DCMF_DEFAULT_GLOBALALLREDUCE_PROTOCOL);

    control_init(DCMF_DEFAULT_CONTROL_PROTOCOL, DCMF_DEFAULT_NETWORK);

    memregion_init(2);

    put_init(DCMF_DEFAULT_PUT_PROTOCOL, DCMF_TORUS_NETWORK);

    barrier();

    if (myrank == 0)
    {
        printf("Number of processes : %d\n", nranks);
        printf("Latency (usec) of Flush with Puts\n");
        fflush(stdout);
    }

    flush_put();

    memregion_finalize();

    multicast_init(DCMF_MEMFIFO_DMA_MSEND_PROTOCOL, 240);

    if (myrank == 0)
    {
        printf("Latency (usec) of Flush with Multicast\n");
        fflush(stdout);
    }

    flush_multicast();

    /*
     flush_send_init(DCMF_EAGER_SEND_PROTOCOL, DCMF_TORUS_NETWORK);
     ack_init();
     ack_control_init();

     if(myrank == 0) {
     printf("Latency (usec) of Flush with Send\n");
     fflush(stdout);
     }

     flush_send(1);
     */

    barrier();

    if (myrank == 0)
    {
        printf("Benchmark complete\n");
        fflush(stdout);
    }

    DCMF_Messager_finalize();

    return 0;
}
