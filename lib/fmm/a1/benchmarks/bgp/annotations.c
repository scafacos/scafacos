#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int main()
{
    int n=1024;
    int i;
    double* x = malloc(n*sizeof(double));
    double* y = malloc(n*sizeof(double));
    double* a = malloc(n*sizeof(double));
    double* b = malloc(n*sizeof(double));

    printf("START\n");

    for (i=0;i<n;i++) x[i] = 0.0;
    for (i=0;i<n;i++) y[i] = 0.0;
    for (i=0;i<n;i++) a[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;
    for (i=0;i<n;i++) b[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;

    for (i=0;i<n;i++)
    {
        y[i] = a[i] - b[i];
    }

    for (i=0;i<n;i+=4)
    {
        {
            double _Complex t1, t2, t3;
            t1 = __lfpd(&a[i]);
            t2 = __lfpd(&b[i]);
            t3 = __fpsub(t1,t2);
            __stfpd(&x[i],t3);
        }
        {
            double _Complex t1, t2, t3;
            t1 = __lfpd(&a[i+2]);
            t2 = __lfpd(&b[i+2]);
            t3 = __fpsub(t1,t2);
            __stfpd(&x[i+2],t3);
        }
    }

    printf("VERIFYING\n");

    for (i=0;i<n;i++)
    {
        if (x[i] != y[i])
            printf("%4d %30.15f %30.15f\n",i,x[i],y[i]);
    }

    printf("FINISH\n");

    return(0);
}
