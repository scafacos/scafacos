#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define likely_if(x) if(__builtin_expect(x,1))
#define unlikely_if(x) if(__builtin_expect(x,0))

static __inline__ unsigned long long getticks(void)
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

static __inline__ void aligned_intrinsics_accumulate(int count, double* target, double* source, double scale)
{
    int i, j, num;

    #pragma disjoint (*target, *source)
    __alignx(16,target);
    __alignx(16,source);

    num = ( count - 16 );
    for (i = 0; i <= num; i += 16)
    {
        {
            double _Complex s0, s2, s4, s6, s8, s10, s12, s14;
            double _Complex t0, t2, t4, t6, t8, t10, t12, t14;

            s0  = __lfpd(&source[i     ]);
            s2  = __lfpd(&source[i +  2]);
            s4  = __lfpd(&source[i +  4]);
            s6  = __lfpd(&source[i +  6]);
            s8  = __lfpd(&source[i +  8]);
            s10 = __lfpd(&source[i + 10]);
            s12 = __lfpd(&source[i + 12]);
            s14 = __lfpd(&source[i + 14]);

            t0  = __lfpd(&target[i     ]);
            t2  = __lfpd(&target[i +  2]);
            t4  = __lfpd(&target[i +  4]);
            t6  = __lfpd(&target[i +  6]);
            t8  = __lfpd(&target[i +  8]);
            t10 = __lfpd(&target[i + 10]);
            t12 = __lfpd(&target[i + 12]);
            t14 = __lfpd(&target[i + 14]);

            t0  = __fxcpmadd( t0,  s0, scale);
            t2  = __fxcpmadd( t2,  s2, scale);
            t4  = __fxcpmadd( t4,  s4, scale);
            t6  = __fxcpmadd( t6,  s6, scale);
            t8  = __fxcpmadd( t8,  s8, scale);
            t10 = __fxcpmadd(t10, s10, scale);
            t12 = __fxcpmadd(t12, s12, scale);
            t14 = __fxcpmadd(t14, s14, scale);

            __stfpd(&target[i     ],  t0);
            __stfpd(&target[i +  2],  t2);
            __stfpd(&target[i +  4],  t4);
            __stfpd(&target[i +  6],  t6);
            __stfpd(&target[i +  8],  t8);
            __stfpd(&target[i + 10], t10);
            __stfpd(&target[i + 12], t12);
            __stfpd(&target[i + 14], t14);
        }
    }
    num = ( count & 15 );
    for ( j = 0; j < num; j++ )
    {
        target[j+i] += ( scale * source[j+i] );
    }
}

static __inline__ void aligned_unrolled_accumulate(int count, double* target, double* source, double scale)
{
    int i;

    #pragma disjoint (*target, *source)
    __alignx(16,target);
    __alignx(16,source);
    for ( i = 0; i < count; i++ ) target[i] += ( scale * source[i] );
}

static __inline__ void unaligned_unrolled_accumulate(int count, double* target, double* source, double scale)
{
    int i;

    #pragma disjoint (*target, *source)
    __alignx( 8,target);
    __alignx(16,source);
    for ( i = 0; i < count; i++ ) target[i] += ( scale * source[i] );
}

static __inline__ void generic_accumulate(int count, double* target, double* source, double scale)
{
    int i;
    for ( i = 0; i < count; i++ ) target[i] += ( scale * source[i] );
}

static __inline__ void optimized_accumulate(int count, double* target, double* source, double scale)
{
    /* are both x & f 16 byte aligned? */
    //if ( ((((int) x) | ((int) f)) & 0xf) == 0){};
    /* alternative implementation: */
    // if (((int) x % 16 == 0) && ((int) f % 16) == 0)

    likely_if (count>1)
    {
          assert( ( ( (int) source ) % 16 ) == 0 );
        if      ( ( ( (int) target ) % 16 ) == 0 ) aligned_intrinsics_accumulate(count,target,source,scale);
        else if ( ( ( (int) target ) %  8 ) == 0 ) unaligned_unrolled_accumulate(count,target,source,scale);
        else { fprintf(stderr,"IMPROPER ALIGNMENT\n"); abort(); }
    }
    else
    {
        if (count==1) target[0] += ( scale * source[0] );
    }
}

int main(int argc, char* argv[])
{
    printf("TEST OF DP ACCUMULATE USING ALIGNMENT\n");

    int dim;
    int max = ( argc>1 ? atoi(argv[1]) : 5000 );

    printf("%6s %44s %33s\n","","-----------aligned target------------","-----unaligned target-----");
    printf("%6s %10s %10s %10s %10s %10s %10s %10s\n","dim","gen","opt","unal-ur","al-ur","gen","opt","unal-ur");

    for (dim=1;dim<max;dim++)
    {
        int i;

        unsigned long long t0, t1, dt0, dt1, dt2, dt3, dt4, dt5, dt6, dt7;

        double  scale = 0.1;

        double* target0;
        double* source0;
        double* target1;
        double* source1;
        double* target2;
        double* source2;
        double* target3;
        double* source3;

        posix_memalign((void**)&target0, 16*sizeof(double), dim*sizeof(double));
        posix_memalign((void**)&source0, 16*sizeof(double), dim*sizeof(double));
        posix_memalign((void**)&target1, 16*sizeof(double), dim*sizeof(double));
        posix_memalign((void**)&source1, 16*sizeof(double), dim*sizeof(double));
        posix_memalign((void**)&target2, 16*sizeof(double), dim*sizeof(double));
        posix_memalign((void**)&source2, 16*sizeof(double), dim*sizeof(double));
        posix_memalign((void**)&target3, 16*sizeof(double), dim*sizeof(double));
        posix_memalign((void**)&source3, 16*sizeof(double), dim*sizeof(double));

        for (i=0;i<dim;i++) target0[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;
        for (i=0;i<dim;i++) source0[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;
        for (i=0;i<dim;i++) target1[i] = target0[i];
        for (i=0;i<dim;i++) source1[i] = source0[i];
        for (i=0;i<dim;i++) target2[i] = target0[i];
        for (i=0;i<dim;i++) source2[i] = source0[i];
        for (i=0;i<dim;i++) target3[i] = target0[i];
        for (i=0;i<dim;i++) source3[i] = source0[i];

        t0 = getticks();
        generic_accumulate(dim, target0, source0, scale);
        t1 = getticks();
        dt0 = t1 - t0;

        t0 = getticks();
        optimized_accumulate(dim, target1, source1, scale);
        t1 = getticks();
        dt1 = t1 - t0;

        t0 = getticks();
        unaligned_unrolled_accumulate(dim, target2, source2, scale);
        t1 = getticks();
        dt2 = t1 - t0;

        t0 = getticks();
        aligned_unrolled_accumulate(dim, target3, source3, scale);
        t1 = getticks();
        dt3 = t1 - t0;

        for (i=0;i<dim;i++) target0[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;
        for (i=0;i<dim;i++) source0[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;
        for (i=0;i<dim;i++) target1[i] = target0[i];
        for (i=0;i<dim;i++) source1[i] = source0[i];

        t0 = getticks();
        generic_accumulate(dim-1, &target0[1], source0, scale);
        t1 = getticks();
        dt4 = t1 - t0;

        t0 = getticks();
        optimized_accumulate(dim-1, &target1[1], source1, scale);
        t1 = getticks();
        dt5 = t1 - t0;

        t0 = getticks();
        unaligned_unrolled_accumulate(dim, target2, source2, scale);
        t1 = getticks();
        dt6 = t1 - t0;

        // VERIFICATION
        for (i=0;i<dim;i++)
        {
            if (target0[i] != target1[i])
            {
                printf("%4d %30.15f %30.15f\n",i,target0[i],target1[i]);
            }
        }

        printf("%6d %10llu %10llu %10llu %10llu %10llu %10llu %10llu\n",dim,dt0,dt1,dt2,dt3,dt4,dt5,dt6);

        free(target0);
        free(source0);
        free(target1);
        free(source1);
        free(target2);
        free(source2);
        free(target3);
        free(source3);
    }

    printf("%6s %10s %10s %10s %10s %10s %10s %10s\n","dim","gen","opt","unal-ur","al-ur","gen","opt","unal-ur");
    printf("%6s %44s %33s\n","","-----------aligned target------------","-----unaligned target-----");

    printf("ALL DONE\n");

    return(0);
}
