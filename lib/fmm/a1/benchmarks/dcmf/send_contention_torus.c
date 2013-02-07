#include "bench.h"
#include "float.h"

void send_static(int pick)
{

    unsigned int size, i, dst, src;
    DCMF_Request_t snd_req[ITERATIONS];
    DCMF_Callback_t snd_done;
    DCQuad msginfo[ITERATIONS];
    DCMF_NetworkCoord_t myaddr, dstaddr;
    DCMF_Network ntwk;

    DCMF_Request_t* rcv_req2 = (DCMF_Request_t *) malloc(sizeof(DCMF_Request_t) * ITERATIONS);
    rcv_req_index = 0;

    DCMF_Messager_rank2network(myrank, DCMF_TORUS_NETWORK, &myaddr);

    dstaddr.network = myaddr.network;
    dstaddr.torus.x = (myaddr.torus.x + 4) % 8;
    dstaddr.torus.y = (myaddr.torus.y + 4) % 8;
    dstaddr.torus.z = (myaddr.torus.z + 4) % 8;
    dstaddr.torus.t = myaddr.torus.t;

    DCMF_Messager_network2rank(&dstaddr, &dst, &ntwk);
    src = dst;

    snd_done.function = done;
    snd_done.clientdata = (void *) &snd_active;

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer,
                "%20s  %20s  %20s  %20s  %20s",
                "Msg Size",
                "Max Latency (us)",
                "Min Latency (us)",
                "Avg Latency (us)",
                "Latency at 0 (usec)");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    for (size = 1; size <= MAX_MSG_SIZE; size *= 2)
    {

        snd_active = 0;
        snd_rcv_active = 0;
        ack_rcv_active = 0;

        barrier();

        if ((myaddr.torus.x) % pick == 0)
        {

            snd_active += ITERATIONS;
            snd_rcv_active += ITERATIONS;

            t_start = DCMF_Timebase();

            for (i = 0; i < ITERATIONS; i++)
            {
                DCMF_Send(&snd_reg,
                          &snd_req[i],
                          snd_done,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          dst,
                          size,
                          source + i * size,
                          &msginfo[i],
                          1);
            }

            while (snd_active > 0 || snd_rcv_active > 0)
                DCMF_Messager_advance();

            snd_active += 1;
            ack_rcv_active += 1;
            DCMF_Send(&ack_reg,
                      snd_req,
                      snd_done,
                      DCMF_SEQUENTIAL_CONSISTENCY,
                      src,
                      1,
                      source,
                      msginfo,
                      1);

            while (snd_active > 0 || ack_rcv_active > 0)
                DCMF_Messager_advance();

            t_stop = DCMF_Timebase();
            t_usec = ((t_stop - t_start) / clockMHz);
            t_usec = t_usec / ITERATIONS;

            barrier();

            allreduce(-1,
                      (char *) &t_usec,
                      (char *) &t_max,
                      1,
                      DCMF_DOUBLE,
                      DCMF_MAX);
            allreduce(-1,
                      (char *) &t_usec,
                      (char *) &t_min,
                      1,
                      DCMF_DOUBLE,
                      DCMF_MIN);
            allreduce(-1,
                      (char *) &t_usec,
                      (char *) &t_avg,
                      1,
                      DCMF_DOUBLE,
                      DCMF_SUM);

            barrier();

            if (myrank == 0)
            {
                t_avg = t_avg / (nranks / pick);
                printf("%20d %20.0f  %20.0f  %20.0f %20.0f\n",
                       size,
                       t_max,
                       t_min,
                       t_avg,
                       t_usec);
                fflush(stdout);
            }

        }
        else
        {

            double d_min, d_max, d_avg;
            d_max = DBL_MAX;
            d_min = DBL_MIN;
            d_avg = 0;

            barrier();

            allreduce(-1,
                      (char *) &d_min,
                      (char *) &t_max,
                      1,
                      DCMF_DOUBLE,
                      DCMF_MAX);
            allreduce(-1,
                      (char *) &d_max,
                      (char *) &t_min,
                      1,
                      DCMF_DOUBLE,
                      DCMF_MIN);
            allreduce(-1,
                      (char *) &d_avg,
                      (char *) &t_avg,
                      1,
                      DCMF_DOUBLE,
                      DCMF_SUM);

            barrier();

        }

    }
}

void send_dynamic(int pick)
{

    unsigned int size, i, dst, src;
    DCMF_Request_t snd_req[ITERATIONS];
    DCMF_Callback_t snd_done;
    DCQuad msginfo[ITERATIONS];
    DCMF_NetworkCoord_t myaddr, dstaddr;
    DCMF_Network ntwk;

    DCMF_Request_t* rcv_req2 = (DCMF_Request_t*) malloc(sizeof(DCMF_Request_t) * ITERATIONS);
    rcv_req_index = 0;

    DCMF_Messager_rank2network(myrank, DCMF_TORUS_NETWORK, &myaddr);

    dstaddr.network = myaddr.network;
    dstaddr.torus.x = (myaddr.torus.x + 4) % 8;
    dstaddr.torus.y = (myaddr.torus.y + 4) % 8;
    dstaddr.torus.z = (myaddr.torus.z + 4) % 8;
    dstaddr.torus.t = myaddr.torus.t;

    DCMF_Messager_network2rank(&dstaddr, &dst, &ntwk);
    src = dst;

    snd_done.function = done;
    snd_done.clientdata = (void *) &snd_active;

    if (myrank == 0)
    {
        char buffer[50];
        sprintf(buffer,
                "%20s  %20s  %20s  %20s  %20s",
                "Msg Size",
                "Max Latency (us)",
                "Min Latency (us)",
                "Avg Latency (us)",
                "Latency at 0 (usec)");
        printf("%s \n", buffer);
        fflush(stdout);
    }

    for (size = 1; size <= MAX_MSG_SIZE; size *= 2)
    {

        snd_active = 0;
        snd_rcv_active = 0;
        ack_rcv_active = 0;

        barrier();

        if ((myaddr.torus.x) % pick == 0)
        {

            snd_active += ITERATIONS;
            snd_rcv_active += ITERATIONS;

            t_start = DCMF_Timebase();

            for (i = 0; i < ITERATIONS; i++)
            {
                DCMF_Send(&snd_reg,
                          &snd_req[i],
                          snd_done,
                          DCMF_RELAXED_CONSISTENCY,
                          dst,
                          size,
                          source + i * size,
                          &msginfo[i],
                          1);
            }

            while (snd_active > 0 || snd_rcv_active > 0)
                DCMF_Messager_advance();

            snd_active += 1;
            ack_rcv_active += 1;
            DCMF_Send(&ack_reg,
                      snd_req,
                      snd_done,
                      DCMF_RELAXED_CONSISTENCY,
                      src,
                      1,
                      source,
                      msginfo,
                      1);

            while (snd_active > 0 || ack_rcv_active > 0)
                DCMF_Messager_advance();

            t_stop = DCMF_Timebase();
            t_usec = ((t_stop - t_start) / clockMHz);
            t_usec = t_usec / ITERATIONS;

            barrier();

            allreduce(-1,
                      (char *) &t_usec,
                      (char *) &t_max,
                      1,
                      DCMF_DOUBLE,
                      DCMF_MAX);
            allreduce(-1,
                      (char *) &t_usec,
                      (char *) &t_min,
                      1,
                      DCMF_DOUBLE,
                      DCMF_MIN);
            allreduce(-1,
                      (char *) &t_usec,
                      (char *) &t_avg,
                      1,
                      DCMF_DOUBLE,
                      DCMF_SUM);

            barrier();

            if (myrank == 0)
            {
                t_avg = t_avg / (nranks / pick);
                printf("%20d %20.0f  %20.0f  %20.0f %20.0f\n",
                       size,
                       t_max,
                       t_min,
                       t_avg,
                       t_usec);
                fflush(stdout);
            }

        }
        else
        {

            double d_min, d_max, d_avg;
            d_max = DBL_MAX;
            d_min = DBL_MIN;
            d_avg = 0;

            barrier();

            allreduce(-1,
                      (char *) &d_min,
                      (char *) &t_max,
                      1,
                      DCMF_DOUBLE,
                      DCMF_MAX);
            allreduce(-1,
                      (char *) &d_max,
                      (char *) &t_min,
                      1,
                      DCMF_DOUBLE,
                      DCMF_MIN);
            allreduce(-1,
                      (char *) &d_avg,
                      (char *) &t_avg,
                      1,
                      DCMF_DOUBLE,
                      DCMF_SUM);

            barrier();

        }

    }
}

int main()
{
    DCMF_Messager_initialize();

    init();

    barrier_init(DCMF_DEFAULT_GLOBALBARRIER_PROTOCOL);

    allreduce_init(DCMF_DEFAULT_GLOBALALLREDUCE_PROTOCOL);

    source = (char *) malloc(MAX_MSG_SIZE * ITERATIONS);
    target = (char *) malloc(MAX_MSG_SIZE * ITERATIONS);
    target_index = 0;

    send_init(DCMF_RZV_SEND_PROTOCOL, DCMF_TORUS_NETWORK);

    ack_init();

    barrier();

    if (myrank == 0)
    {
        printf("Send Latency (usec) with static routing - 100p sturation \n");
        fflush(stdout);
    }
    send_static(1);

    if (myrank == 0)
    {
        printf("Send Latency (usec) with static routing - 50p sturation \n");
        fflush(stdout);
    }
    send_static(2);

    if (myrank == 0)
    {
        printf("Send Latency (usec) with static routing - 25p sturation \n");
        fflush(stdout);
    }
    send_static(4);

    if (myrank == 0)
    {
        printf("Send Latency (usec) with dynamic routing - 100p sturation \n");
        fflush(stdout);
    }
    send_dynamic(1);

    if (myrank == 0)
    {
        printf("Send Latency (usec) with dynamic routing - 50p sturation \n");
        fflush(stdout);
    }
    send_dynamic(2);

    if (myrank == 0)
    {
        printf("Send Latency (usec) with dynamic routing - 25p sturation \n");
        fflush(stdout);
    }
    send_dynamic(4);

    barrier();

    if (myrank == 0)
    {
        printf("[%d] Benchmark complete\n", myrank);
        fflush(stdout);
    }

    free((void *) source);
    free((void *) target);

    DCMF_Messager_finalize();

    return 0;
}
