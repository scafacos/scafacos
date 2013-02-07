#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

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

int main()
{
    int i,j,count=10,dim=1024;
    unsigned long long t0,t1;
    double* x;
    double* y1;
    double* y2;

    posix_memalign((void**)&x , 16*sizeof(double), dim*sizeof(double));
    posix_memalign((void**)&y1, 16*sizeof(double), dim*sizeof(double));
    posix_memalign((void**)&y2, 16*sizeof(double), dim*sizeof(double));

    printf("TESTING IMPLEMENTATIONS OF VECTOR ABSOLUTE VALUE\n");

    for (i=0;i<dim;i++) x[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;
    //for (i=0;i<dim;i++) printf("x[%4d] = %30.15f\n",i,x[i]);
    for (i=0;i<dim;i++) y1[i] = 0.0;
    for (i=0;i<dim;i++) y2[i] = 0.0;

    printf("if-else VERSION\n");

    // WARM-UP
    for (i=0;i<dim;i++)
    {
        if (x[i]>0.0) y1[i] = x[i];
        else          y1[i] = -x[i];
    }

    // TIMING
    t0 = getticks();
    for (j=0;j<count;j++)
    {
        for (i=0;i<dim;i++)
        {
            if (x[i]>0.0) y1[i] = x[i];
            else          y1[i] = -x[i];
        }
    }
    t1 = getticks();
    printf("time for if-else version = %30llu\n",t1-t0);

    printf("a?:b:c VERSION\n");

    // WARM-UP
    for (i=0;i<dim;i++)
    {
        y1[i] = x[i]>0.0 ? x[i] : -x[i];
    }

    // TIMING
    t0 = getticks();
    for (j=0;j<count;j++)
    {
        for (i=0;i<dim;i++)
        {
            y1[i] = x[i]>0.0 ? x[i] : -x[i];
        }
    }
    t1 = getticks();
    printf("time for  a?b:c  version = %30llu\n",t1-t0);

    printf("FABS VERSION\n");

    // WARM-UP
    for (i=0;i<dim;i+=2)
    {
        y1[i] = fabs(x[i]);
    }

    // TIMING
    t0 = getticks();
    for (j=0;j<count;j++)
    {
        for (i=0;i<dim;i+=2)
        {
            y1[i] = fabs(x[i]);
        }
    }
    t1 = getticks();
    printf("time for  fabs() version = %30llu\n",t1-t0);

    double* a;
    double* b;
    posix_memalign((void**)&a, 16*sizeof(double), dim*sizeof(double));
    posix_memalign((void**)&b, 16*sizeof(double), dim*sizeof(double));

    for (i=0;i<dim;i++) a[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;
    for (i=0;i<dim;i++) b[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;
    for (i=0;i<dim;i++) y1[i] = 0.0;
    for (i=0;i<dim;i++) y2[i] = 0.0;

    // WARM-UP
    for (i=0;i<dim;i++)
    {
        y1[i] = fabs(a[i]);// - b[i];
    }

    printf("ASM1 VERSION\n");

    // TIMING
    t0 = getticks();
    for (j=0;j<count;j++)
    {
        for (i=0;i<dim;i+=2)
        {
            __stfpd(&y2[i], __fpabs( __lfpd(&a[i]) ) );
        }
    }
    t1 = getticks();
    printf("time for   ASM1  version = %30llu\n",t1-t0);

    printf("VERIFYING\n");

    for (i=0;i<dim;i++)
    {
        if (y1[i] != y2[i])
        {
            printf("%4d %30.15f %30.15f\n",i,y1[i],y2[i]);
        }
    }
    for (i=0;i<dim;i++) y2[i] = 0.0;

    printf("ASM2 VERSION\n");

    // TIMING
    t0 = getticks();
    for (j=0;j<count;j++)
    {
        for (i=0;i<dim;i+=4)
        {
            __stfpd(&y2[i  ], __fpabs( __lfpd(&a[i  ]) ) );
            __stfpd(&y2[i+2], __fpabs( __lfpd(&a[i+2]) ) );
        }
    }
    t1 = getticks();
    printf("time for   ASM2  version = %30llu\n",t1-t0);

    printf("VERIFYING\n");

    for (i=0;i<dim;i++)
    {
        if (y1[i] != y2[i])
        {
            printf("%4d %30.15f %30.15f\n",i,y1[i],y2[i]);
        }
    }
    for (i=0;i<dim;i++) y2[i] = 0.0;

    printf("ASM3 VERSION\n");

    // TIMING
    t0 = getticks();
    for (j=0;j<count;j++)
    {
        for (i=0;i<dim;i+=8)
        {
            __stfpd(&y2[i  ], __fpabs( __lfpd(&a[i  ]) ) );
            __stfpd(&y2[i+2], __fpabs( __lfpd(&a[i+2]) ) );
            __stfpd(&y2[i+4], __fpabs( __lfpd(&a[i+4]) ) );
            __stfpd(&y2[i+6], __fpabs( __lfpd(&a[i+6]) ) );
        }
    }
    t1 = getticks();
    printf("time for   ASM3  version = %30llu\n",t1-t0);

    printf("VERIFYING\n");

    for (i=0;i<dim;i++)
    {
        if (y1[i] != y2[i])
        {
            printf("%4d %30.15f %30.15f\n",i,y1[i],y2[i]);
        }
    }
    for (i=0;i<dim;i++) y2[i] = 0.0;

    printf("ASM4 VERSION\n");

    // TIMING
    t0 = getticks();
    for (j=0;j<count;j++)
    {
        for (i=0;i<dim;i+=16)
        {
            __stfpd(&y2[i   ], __fpabs( __lfpd(&a[i   ]) ) );
            __stfpd(&y2[i+ 2], __fpabs( __lfpd(&a[i+ 2]) ) );
            __stfpd(&y2[i+ 4], __fpabs( __lfpd(&a[i+ 4]) ) );
            __stfpd(&y2[i+ 6], __fpabs( __lfpd(&a[i+ 6]) ) );
            __stfpd(&y2[i+ 8], __fpabs( __lfpd(&a[i+ 8]) ) );
            __stfpd(&y2[i+10], __fpabs( __lfpd(&a[i+10]) ) );
            __stfpd(&y2[i+12], __fpabs( __lfpd(&a[i+12]) ) );
            __stfpd(&y2[i+14], __fpabs( __lfpd(&a[i+14]) ) );
        }
    }
    t1 = getticks();
    printf("time for   ASM4  version = %30llu\n",t1-t0);

    printf("VERIFYING\n");

    for (i=0;i<dim;i++)
    {
        if (y1[i] != y2[i])
        {
            printf("%4d %30.15f %30.15f\n",i,y1[i],y2[i]);
        }
    }
    for (i=0;i<dim;i++) y2[i] = 0.0;

    printf("ALL DONE\n");

    return(0);
}
