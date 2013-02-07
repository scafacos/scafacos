#include <stdio.h>

#if defined(__powerpc__)

__inline__ unsigned long long getticks(void)
{
    unsigned long long int result=0;
    unsigned long int upper, lower,tmp;
    __asm__ volatile(
            "0:                  \n"
            "\tmftbu   %0           \n"
            "\tmftb    %1           \n"
            "\tmftbu   %2           \n"
            "\tcmpw    %2,%0        \n"
            "\tbne     0b         \n"
            : "=r"(upper),"=r"(lower),"=r"(tmp)
    );
    result = upper;
    result = result<<32;
    result = result|lower;

    return(result);
}

#else

#error PPC required

#endif

#define SIZE 100*1024*1024

int main(int argc, char **argv)
{

    int i;
    double source[SIZE] __attribute__((__aligned__(16)));
    double target[SIZE] __attribute__((__aligned__(16)));
    double scaling[2] __attribute__((__aligned__(16)));
    double _Complex csource __attribute__((__aligned__(16)));
    double _Complex cscaling __attribute__((__aligned__(16)));
    double _Complex ctarget __attribute__((__aligned__(16)));
    unsigned long long t_start, t_stop;

    scaling[0] = 2.0;
    scaling[1] = 2.0;
    for (i = 0; i < SIZE; i++)
    {
        source[i] = i + 1;
        target[i] = i + 1;
    }

    t_start = getticks();
    cscaling = __lfpd(scaling);
    for (i = 0; i < SIZE; i = i + 2)
    {
        csource = __lfpd(&(source[i]));
        ctarget = __lfpd(&(target[i]));
        ctarget = __fpmadd(ctarget, cscaling, csource);
        __stfpd(&(target[i]), ctarget);
    }
    t_stop = getticks();

    printf("Time (in msec) with double hummmer: %llu target: %f\n", t_stop
            - t_start, target[0]);
    fflush(stdout);

    for (i = 0; i < SIZE; i++)
    {
        target[i] = i + 1;
    }

    t_start = getticks();
    for (i = 0; i < SIZE; i++)
    {
        target[i] = target[i] + source[i] * scaling[0];
    }
    t_stop = getticks();

    printf("Time (in msec) with regular arithmetic: %llu target: %f\n", t_stop
            - t_start, target[0]);
    fflush(stdout);

    return 0;
}
