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

void put_direct(int dim1, int dim2)
{

    memregion_init(MAX_DIM * MAX_DIM * sizeof(double));

    if (myrank == 0)
    {

        DCMF_Request_t put_req[dim1];
        DCMF_Callback_t put_done, put_ack;
        int done_count, ack_count;
        int j, k, dst;

        put_done.function = done;
        put_done.clientdata = (void *) &done_count;
        put_ack.function = done;
        put_ack.clientdata = (void *) &ack_count;

        t_usec = 0;
        for (dst = 0; dst < nranks; dst++)
        {

            if (dst != myrank)
            {

                t_start = DCMF_Timebase();
                for (j = 0; j < ITERS; j++)
                {

                    ack_count = dim1;
                    for (k = 0; k < dim1; k++)
                    {
                        DCMF_Put(&put_reg,
                                 &put_req[k],
                                 put_done,
                                 DCMF_SEQUENTIAL_CONSISTENCY,
                                 dst,
                                 dim2 * sizeof(double),
                                 memregion[0],
                                 memregion[dst],
                                 k * MAX_DIM * sizeof(double),
                                 k * MAX_DIM * sizeof(double),
                                 put_ack);
                    }
                    while (ack_count)
                        DCMF_Messager_advance();

                }
                t_stop = DCMF_Timebase();
                t_usec += ((t_stop - t_start) / clockMHz);
            }

        }
        printf("%20.0f ", t_usec / ((nranks - 1) * ITERS));
        fflush(stdout);

    }

    barrier();

    memregion_finalize();

}

void send_pack(int dim1, int dim2)
{

    memregion_init(MAX_DIM * MAX_DIM * sizeof(double));

    if (myrank == 0)
    {

        DCMF_Request_t snd_req;
        DCMF_Callback_t snd_done;
        DCQuad msginfo;
        int done_count;
        unsigned int j, k, dst, size;
        char *pack_buffer;
        struct noncontig_header pack_header;

        pack_buffer = (char *) malloc(sizeof(struct noncontig_header) + dim1
                * dim2 * sizeof(double));

        snd_done.function = done;
        snd_done.clientdata = (void *) &done_count;

        for (dst = 0; dst < nranks; dst++)
        {

            barrier();

            if (dst != myrank)
            {

                for (j = 0; j < ITERS; j++)
                {

                    done_count = 1;

                    DCMF_Memregion_query(memregion[dst],
                                         &size,
                                         &(pack_header.vaddress));
                    pack_header.stride = MAX_DIM * sizeof(double);
                    pack_header.d1 = dim1;
                    pack_header.d2 = dim2 * sizeof(double);
                    memcpy(pack_buffer, &pack_header, sizeof(pack_header));
                    int msgsize = sizeof(pack_header);
                    for (k = 0; k < dim1; k++)
                    {
                        memcpy(pack_buffer + msgsize, window + k * MAX_DIM
                                * sizeof(double), dim2 * sizeof(double));
                        msgsize += dim2 * sizeof(double);
                    }

                    DCMF_Send(&snd_noncontig_reg,
                              &snd_req,
                              snd_done,
                              DCMF_SEQUENTIAL_CONSISTENCY,
                              dst,
                              msgsize,
                              pack_buffer,
                              &msginfo,
                              1);

                    while (done_count)
                        DCMF_Messager_advance();

                }

            }

        }

        free(pack_buffer);

        t_usec = 0;
        barrier();
        allreduce(-1,
                  (char *) &t_usec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

        t_avg = t_avg / (nranks - 1);
        printf("%20.0f \n", t_avg);
        fflush(stdout);

    }
    else
    {

        int dst;

        t_usec = 0;
        for (dst = 0; dst < nranks; dst++)
        {

            barrier();

            if (dst == myrank)
            {

                t_start = DCMF_Timebase();

                snd_rcv_noncontig_active = ITERS;
                while (snd_rcv_noncontig_active)
                    DCMF_Messager_advance();

                t_stop = DCMF_Timebase();
                t_usec = ((t_stop - t_start) / clockMHz);
                t_usec = t_usec / ITERS;

            }
        }

        barrier();
        allreduce(-1,
                  (char *) &t_usec,
                  (char *) &t_avg,
                  1,
                  DCMF_DOUBLE,
                  DCMF_SUM);
        barrier();

    }

    memregion_finalize();

}

void manytomany()
{

    memregion_init(MAX_DIM * MAX_DIM * sizeof(double));

    DCMF_Request_t snd_req;
    DCMF_Callback_t snd_done, m2m_done;
    DCQuad msginfo;
    int done_count;
    unsigned int dst, size;
    struct noncontig_header pack_header;
    unsigned sndlens, snddispls, sndcounters, permutation;
    unsigned ranks;

    dst = (myrank + 1) % nranks;

    snd_done.function = done;
    snd_done.clientdata = (void *) &done_count;

    m2m_done.function = done;
    m2m_done.clientdata = (void *) &m2m_active;

    barrier();

    done_count = 1;
    snd_rcv_manytomany_active += 1;

    DCMF_Memregion_query(memregion[dst], &size, &(pack_header.vaddress));
    pack_header.stride = MAX_DIM * sizeof(double);
    pack_header.d1 = 1;
    pack_header.d2 = 64 * sizeof(double);

    DCMF_Send(&snd_manytomany_reg,
              &snd_req,
              snd_done,
              DCMF_SEQUENTIAL_CONSISTENCY,
              dst,
              sizeof(struct noncontig_header),
              (char *) &pack_header,
              &msginfo,
              1);

    while (done_count || snd_rcv_manytomany_active)
        DCMF_Messager_advance();

    printf("[%d] Done sending and receiving header information \n", myrank);
    fflush(stdout);

    sndlens = pack_header.d2;
    snddispls = 0;
    ranks = dst;
    permutation = dst;

    m2m_active = 1;
    m2m_rcv_active = 1;

    DCMF_Manytomany(&m2m_reg,
                    &m2m_req,
                    m2m_done,
                    DCMF_SEQUENTIAL_CONSISTENCY,
                    0,
                    myrank,
                    window,
                    &sndlens,
                    &snddispls,
                    &sndcounters,
                    &ranks,
                    &permutation,
                    1);

    printf("[%d] posted manytomany \n", myrank);
    fflush(stdout);

    while (m2m_active)
        DCMF_Messager_advance();

    printf("[%d] done sending manytomany \n", myrank);
    fflush(stdout);

    while (m2m_rcv_active)
        DCMF_Messager_advance();

    printf("[%d] done sending manytomany \n", myrank);
    fflush(stdout);

    memregion_finalize();
}

int main()
{
    DCMF_Messager_initialize();
    /* int dim1, dim2; */
    char buffer[50];

    init();

    barrier_init(DCMF_DEFAULT_GLOBALBARRIER_PROTOCOL);

    allreduce_init(DCMF_DEFAULT_GLOBALALLREDUCE_PROTOCOL);

    control_init(DCMF_DEFAULT_CONTROL_PROTOCOL, DCMF_DEFAULT_NETWORK);

    /* put_init(DCMF_DEFAULT_PUT_PROTOCOL, DCMF_TORUS_NETWORK);

     send_noncontig_init(DCMF_DEFAULT_SEND_PROTOCOL, DCMF_TORUS_NETWORK); */

    send_manytomany_init(DCMF_DEFAULT_SEND_PROTOCOL, DCMF_TORUS_NETWORK);

    manytomany_init(DCMF_MEMFIFO_DMA_M2M_PROTOCOL);

    printf("[%d] Registration complete \n", myrank);
    fflush(stdout);

    if (myrank == 0)
    {
        sprintf(buffer,
                "%20s  %20s  %20s",
                "Dimensions",
                "DirectPut Latency (us)",
                "PackSend Latency (us)");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    barrier();

    /* for(dim1=2; dim1<=512; dim1*=2) {

     for(dim2=1; dim2<=512; dim2*=2) {

     if(myrank == 0) {
     sprintf(buffer, "%dX%d", dim1, dim2);
     printf("%20s  ", buffer);
     fflush(stdout);
     }

     barrier();

     put_direct(dim1, dim2);

     send_pack(dim1, dim2);

     manytomany(dim1, dim2);

     }

     } */

    manytomany();

    printf("[%d] Benchmark complete\n", myrank);
    fflush(stdout);

    DCMF_Messager_finalize();

    return 0;
}
