#include <assert.h>
#include <stdio.h>
#include <bpcore/bgp_atomic_ops.h>
#include <dcmf.h>
#include <pthread.h>

#define NITERS 100

pthread_t pt[3];
pthread_barrier_t pt_bar;

_BGP_Atomic global_atomic;
volatile int shared[4 * NITERS] __attribute__((__aligned__(16)));
volatile int shared_idx __attribute__((__aligned__(16)));

#define A1DI_GLOBAL_ATOMIC_ACQUIRE()                    \
 {                                                      \
   volatile int done __attribute__((__aligned__(16)));  \
   done = 0;                                            \
   do {                                                 \
     while(global_atomic.atom);                         \
     done = _bgp_test_and_set(&global_atomic, 1);       \
   } while(!done);                                      \
 }                                                      \

#define A1DI_GLOBAL_ATOMIC_RELEASE() do{ global_atomic.atom = 0; _bgp_mbar(); }while(0)

void *execute(void * dummy)
{
    int i, idx, coreid;
    unsigned long long t_start, t_stop;
    coreid = Kernel_PhysicalProcessorID();
    //pthread_barrier_wait(&pt_bar);
    t_start = DCMF_Timebase();
    for (i = 0; i < NITERS; i++)
    {
        A1DI_GLOBAL_ATOMIC_ACQUIRE();
        _bgp_dcache_touch_line(&shared_idx);
        idx = shared_idx;
        shared[idx] = coreid;
        idx++;
        shared_idx = idx;
        _bgp_mbar();
        A1DI_GLOBAL_ATOMIC_RELEASE();
    }
    t_stop = DCMF_Timebase();
    printf("Time at id %d is: %lld t_start %lld t_stop %lld\n", coreid, (t_stop
            - t_start), t_start, t_stop);
}

int main()
{
    DCMF_Configure_t conf;

    DCMF_Messager_initialize();

    conf.thread_level = DCMF_THREAD_MULTIPLE;
    conf.interrupts = DCMF_INTERRUPTS_OFF;

    DCMF_Messager_configure(&conf, &conf);

    //pthread_barrier_init(&pt_bar, NULL, 4);
    shared_idx = 0;

    /*
    pthread_create(&pt[0], NULL, &execute, NULL);
    pthread_create(&pt[1], NULL, &execute, NULL);
    pthread_create(&pt[2], NULL, &execute, NULL);
    */

    execute(NULL);

    DCMF_Messager_finalize();

    int i;
    for (i = 0; i < 4 * NITERS; i++)
    {
        printf("%d \t", shared[i]);
    }
    printf("\n");
    fflush(stdout);

    return (0);
}
