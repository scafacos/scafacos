#include <assert.h>
#include <stdio.h>
#include <spi/bgp_SPI.h>
#include <dcmf.h>
#include <pthread.h>

#define NITERS 100

pthread_t pt[3];
pthread_barrier_t pt_bar;

LockBox_Mutex_t global_mutex;
volatile int shared[4 * NITERS] __attribute__((__aligned__(16)));
volatile int shared_idx __attribute__((__aligned__(16)));

void *execute(void * dummy)
{
    int i, idx, coreid;
    unsigned long long t_start, t_stop;
    coreid = Kernel_PhysicalProcessorID();
    //pthread_barrier_wait(&pt_bar);
    t_start = DCMF_Timebase();
    for (i = 0; i < NITERS; i++)
    {
        LockBox_MutexLock(global_mutex);
        _bgp_dcache_touch_line(&shared_idx);
        idx = shared_idx;
        shared[idx] = coreid;
        idx++;
        shared_idx = idx;
        _bgp_mbar();
        LockBox_MutexUnlock(global_mutex);
    }
    t_stop = DCMF_Timebase();
    printf("Time per interation at id %d is: %lld t_start %lld t_stop %lld\n",
           coreid,
           (t_stop - t_start),
           t_start,
           t_stop);
}

int main()
{
    DCMF_Configure_t conf;

    DCMF_Messager_initialize();

    conf.thread_level = DCMF_THREAD_MULTIPLE;
    conf.interrupts = DCMF_INTERRUPTS_OFF;

    DCMF_Messager_configure(&conf, &conf);

    int i;
    for (i = 100; i < 1024; i++)
    {
        if (!LockBox_AllocateMutex(i, &global_mutex, 0, 1, 0)) break;
    }
    if (i == 1024)
    {
        printf("Lockbox allocation failed \n");
        return -1;
    }

    //pthread_barrier_init(&pt_bar, NULL, 4);
    shared_idx = 0;

    /*
    pthread_create(&pt[0], NULL, &execute, NULL);
    pthread_create(&pt[1], NULL, &execute, NULL);
    pthread_create(&pt[2], NULL, &execute, NULL);
    */

    execute(NULL);

    DCMF_Messager_finalize();

    for (i = 0; i < 4 * NITERS; i++)
    {
        printf("%d \t", shared[i]);
    }
    printf("\n");
    fflush(stdout);

    return (0);
}
