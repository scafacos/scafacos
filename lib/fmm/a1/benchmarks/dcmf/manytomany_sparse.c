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

#define MAXNODES   8192
#define MAXMSGSIZE 1024

int send_active[2];
int recv_active[2];

DCMF_Request_t send_req[2] __attribute ((__aligned__(16)));
DCMF_Request_t recv_req[2] __attribute ((__aligned__(16)));

char recv_buf[MAXMSGSIZE + 1] __attribute ((__aligned__(16)));
unsigned recv_len[2];
unsigned recv_displ[2];
unsigned recv_counter[2][2];

char send_buf[MAXMSGSIZE] __attribute ((__aligned__(16)));
unsigned send_rank;
unsigned send_len;
unsigned send_displ;
unsigned send_counter[2];
unsigned permutation;

unsigned rank, size;

void cb_Manytomany_senddone(void * clientdata, DCMF_Error_t *error)
{
    send_active[(unsigned) clientdata]--;
}

void cb_Manytomany_recvdone(void * clientdata, DCMF_Error_t *error)
{
    recv_active[(unsigned) clientdata]--;
}

DCMF_Request_t* cb_recv_Manytomany(unsigned connid,
                                   void * arg,
                                   char ** recvbuf,
                                   unsigned ** recvlens,
                                   unsigned ** recvdispls,
                                   unsigned ** recvcounters,
                                   unsigned * numranks,
                                   unsigned * ridx,
                                   DCMF_Callback_t* const cb_info)
{
    *recvbuf = (char *) recv_buf;
    *recvlens = recv_len;
    *recvdispls = recv_displ;
    *recvcounters = recv_counter[connid];

    *numranks = 2;
    *ridx = 1;

    cb_info->function = cb_Manytomany_recvdone;
    cb_info->clientdata = (void *) connid;

    return &recv_req[connid];
}

int main()
{
    unsigned sendlen;
    unsigned i;

    DCMF_Messager_initialize();

    init();

    barrier_init(DCMF_DEFAULT_GLOBALBARRIER_PROTOCOL);

    DCMF_Protocol_t protocol;
    DCMF_Manytomany_Configuration_t mconfig;
    DCMF_Callback_t cb_done[2];

    mconfig.protocol = DCMF_MEMFIFO_DMA_M2M_PROTOCOL;
    mconfig.cb_recv = cb_recv_Manytomany;
    mconfig.arg = NULL;
    mconfig.nconnections = 2;

    DCMF_Manytomany_register(&protocol, &mconfig);

    cb_done[0].function = cb_Manytomany_senddone;
    cb_done[0].clientdata = (void *) 0;
    cb_done[1].function = cb_Manytomany_senddone;
    cb_done[1].clientdata = (void *) 1;

    for (i = 0; i < MAXMSGSIZE; i++)
    {
        send_buf[i] = '*';
    }

    send_rank = (myrank + 1) % nranks;
    permutation = 0;
    send_len = 0;
    send_displ = 0;

    recv_len[0] = 0;
    recv_len[1] = 1;
    recv_displ[0] = 0;
    recv_displ[1] = MAXMSGSIZE;

    for (sendlen = 1; sendlen <= MAXMSGSIZE; sendlen *= 2)
    {

        send_len = sendlen;
        recv_len[0] = sendlen;

        for (i = 0; i < MAXMSGSIZE + 1; i++)
        {
            recv_buf[i] = '_';
        }

        send_active[0] = 1;
        recv_active[0] = 1;

        DCMF_Manytomany(&protocol,
                        &send_req[0],
                        cb_done[0],
                        DCMF_MATCH_CONSISTENCY,
                        0,
                        0,
                        0,
                        NULL,
                        (char *) send_buf,
                        &send_len,
                        &send_displ,
                        &send_counter[0],
                        &send_rank,
                        &permutation,
                        1);

        while (send_active[0] > 0 || recv_active[0] > 0)
            DCMF_Messager_advance();

        printf("[%d] After manytomany 1 \n", myrank);
        fflush(stdout);

        barrier();

        send_active[1] = 1;
        recv_active[1] = 1;

        DCMF_Manytomany(&protocol,
                        &send_req[1],
                        cb_done[1],
                        DCMF_MATCH_CONSISTENCY,
                        1,
                        0,
                        0,
                        NULL,
                        (char *) send_buf,
                        &send_len,
                        &send_displ,
                        &send_counter[1],
                        &send_rank,
                        &permutation,
                        1);

        while (send_active[1] > 0 || recv_active[1] > 0)
            DCMF_Messager_advance();

        printf("[%d] After manytomany 2 \n", myrank);
        fflush(stdout);

        barrier();

        int niter = 2;
        unsigned long long start, time;

        start = DCMF_Timebase();
        for (i = 0; i < niter; i++)
        {
            int conn = i % 2;
            send_active[conn] = 1;
            recv_active[conn] = 1;

            DCMF_Manytomany(&protocol,
                            &send_req[conn],
                            cb_done[conn],
                            DCMF_MATCH_CONSISTENCY,
                            conn,
                            0,
                            0,
                            NULL,
                            (char *) send_buf,
                            &send_len,
                            &send_displ,
                            &send_counter[conn],
                            &send_rank,
                            &permutation,
                            1);

            while (send_active[conn] > 0 || recv_active[conn] > 0)
                DCMF_Messager_advance();
        }

        time = DCMF_Timebase() - start;
        if (myrank == 0)
        {
            printf("Time For Many-to-many with size %d = %d cycles, buffer = %s \n",
                   sendlen,
                   (unsigned) (time / niter),
                   recv_buf);
            fflush(stdout);
        }

    }

    return 0;
}
