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

#define MAX_MSG_SIZE 1024*1024

#include <bench.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define THREAD_NUM 1 
#define LOCAL_MAX_MSG_SIZE 1024
#define LOCAL_ITERATIONS 1000
#define LOCAL_MAX_BUF_SIZE LOCAL_MAX_MSG_SIZE*LOCAL_ITERATIONS

pthread_barrier_t ptbarrier, ptbarrier1;

void *mrate_test(void *threadid)
{

    unsigned i, channel, base, count, msgsize;
    DCMF_Protocol_t protocol;
    DCMF_Callback_t cb_done;
    uint32_t done_count;
    DCMF_Request_t *req;
    DCMF_Messager_advance_options adv_options;
    DCQuad msginfo;

    req = (DCMF_Request_t *) malloc(sizeof(DCMF_Request_t) * LOCAL_ITERATIONS);

    cb_done.function = done;
    cb_done.clientdata = (void *) &done_count;

    DCMF_Channel_info(&protocol, &base, &count, &channel);
    DCMF_Channel_acquire(base + (long) threadid);

    msgsize = 512;

    if ((long) threadid == 0)
    {
        printf("%10s %20s %30s \n",
               "Thread ID",
               "Message Size",
               "Injection Rate (MBps)");
        fflush(stdout);
    }

    pthread_barrier_wait(&ptbarrier);

    t_start = DCMF_Timebase();

    done_count = LOCAL_ITERATIONS;
    for (i = 0; i < LOCAL_ITERATIONS; i++)
    {

        DCMF_Send(&snd_reg,
                  &req[i],
                  cb_done,
                  DCMF_SEQUENTIAL_CONSISTENCY,
                  (long) threadid + 1,
                  msgsize,
                  source,
                  &msginfo,
                  1);

    }

    adv_options.channel = base + (long) threadid;
    while (done_count > 0)
        DCMF_Messager_advance_expert(adv_options);

    t_stop = DCMF_Timebase();
    t_sec = (t_stop - t_start) / (clockMHz * 1000000);
    printf("%10d %20d %26.4f \n",
           (long) threadid,
           msgsize,
           ((double) LOCAL_ITERATIONS * msgsize) / (t_sec * (double) 1024
                   * 1024));
    fflush(stdout);

    pthread_barrier_wait(&ptbarrier1);

    DCMF_Channel_release();

    return;
}

int main(int argc, void* argv[])
{
    DCMF_Configure_t config;

    config.thread_level = DCMF_THREAD_MULTIPLE;

    DCMF_Messager_initialize();

    DCMF_Messager_configure(&config, &config);

    init();

    if (nranks != (THREAD_NUM + 1))
    {
        printf("This test requires only %d processes \n", (THREAD_NUM + 1));
        fflush(stdout);
        return -1;
    }

    barrier_init(DCMF_DEFAULT_GLOBALBARRIER_PROTOCOL);

    control_init(DCMF_DEFAULT_CONTROL_PROTOCOL, DCMF_DEFAULT_NETWORK);

    memregion_init(LOCAL_MAX_BUF_SIZE * THREAD_NUM);

    get_init(DCMF_DEFAULT_PUT_PROTOCOL, DCMF_TORUS_NETWORK);

    source = (char *) malloc(LOCAL_MAX_BUF_SIZE * THREAD_NUM);
    target = (char *) malloc(LOCAL_MAX_BUF_SIZE * THREAD_NUM);

    send_init(DCMF_DEFAULT_SEND_PROTOCOL, DCMF_TORUS_NETWORK);

    int status;
    long i;

    if (myrank == 0)
    {

        pthread_t threads[THREAD_NUM];
        pthread_barrier_init(&ptbarrier, NULL, THREAD_NUM);
        pthread_barrier_init(&ptbarrier1, NULL, THREAD_NUM);

        for (i = 0; i < THREAD_NUM; i++)
        {
            pthread_create(&threads[i], NULL, mrate_test, (void *) i);
        }

        for (i = 0; i < THREAD_NUM; i++)
        {
            pthread_join(threads[i], (void *) &status);
        }
    }
    else
    {

        snd_rcv_active += LOCAL_ITERATIONS;
        while (snd_rcv_active > 0)
            DCMF_Messager_advance();

    }

    barrier();

    DCMF_Messager_finalize();

    if (myrank == 0)
    {
        printf("Benchmark Complete \n");
        fflush(stdout);
    }

    return (0);
}
;
