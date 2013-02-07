#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

int main(int argc, char* argv[])
{
    fprintf(stderr,"BEGIN TESTING OF DP ACCUMULATE\n");
    printf("%18s %18s %18s %18s %18s %18s\n","dim","basic","1-hummer","2-hummers","4-hummers","8-hummers");

    int k;
    for (k=6;k<20;k++)
    {
        int dim = pow(2,k);

        int count = 10;

        int i,j;

        unsigned long long t0, t1;
        unsigned long long dt0, dt1, dt2, dt3, dt4;

        double* a;
        double* b;
        double* c;

        double  scale = 0.1;

        posix_memalign((void**)&a, 16*sizeof(double), dim*sizeof(double));
        posix_memalign((void**)&b, 16*sizeof(double), dim*sizeof(double));
        posix_memalign((void**)&c, 16*sizeof(double), dim*sizeof(double));

        for (i=0;i<dim;i++) a[i] = 1.0 - 2*(double)rand()/(double)RAND_MAX;

        fprintf(stderr,"BASIC VERSION\n");

        // WARM-UP
        for (i=0;i<dim;i++) b[i] = 0.0;
        for (i=0;i<dim;i++) b[i] += scale*a[i];

        // TIMING
        for (i=0;i<dim;i++) b[i] = 0.0;
        t0 = getticks();
        for (j=0;j<count;j++)
        {
            for (i=0;i<dim;i++) b[i] += scale*a[i];
        }
        t1 = getticks();
        dt0 = t1 - t0;

        fprintf(stderr,"INTRINSICS VERSION 1\n");

        // WARM-UP
        for (i=0;i<dim;i++) c[i] = 0.0;
        for (i=0;i<dim;i+=2)
        {
            __stfpd(&c[i], __fxcpmadd( __lfpd(&c[i]), __lfpd(&a[i]), scale) );
        }

        // TIMING
        for (i=0;i<dim;i++) c[i] = 0.0;
        t0 = getticks();
        for (j=0;j<count;j++)
        {
            for (i=0;i<dim;i+=2)
            {
                __stfpd(&c[i], __fxcpmadd( __lfpd(&c[i]), __lfpd(&a[i]), scale) );
            }
        }
        t1 = getticks();
        dt1 = t1 - t0;

        // VERIFICATION
        for (i=0;i<dim;i++)
        {
            if (b[i] != c[i])
            {
                printf("%4d %30.15f %30.15f\n",i,b[i],c[i]);
            }
        }

        fprintf(stderr,"INTRINSICS VERSION 2\n");

        // WARM-UP
        for (i=0;i<dim;i++) c[i] = 0.0;
        for (i=0;i<dim;i+=4)
        {
            {
                double _Complex a0, a2, c0, c2;
                a0 = __lfpd(&a[i  ]);
                a2 = __lfpd(&a[i+2]);
                c0 = __lfpd(&c[i  ]);
                c2 = __lfpd(&c[i+2]);
                c0 = __fxcpmadd(c0,a0,scale);
                c2 = __fxcpmadd(c2,a2,scale);
                __stfpd(&c[i  ],c0);
                __stfpd(&c[i+2],c2);
            }
        }

        // TIMING
        for (i=0;i<dim;i++) c[i] = 0.0;
        t0 = getticks();
        for (j=0;j<count;j++)
        {
            for (i=0;i<dim;i+=4)
            {
                {
                    double _Complex a0, a2, c0, c2;
                    a0 = __lfpd(&a[i  ]);
                    a2 = __lfpd(&a[i+2]);
                    c0 = __lfpd(&c[i  ]);
                    c2 = __lfpd(&c[i+2]);
                    c0 = __fxcpmadd(c0,a0,scale);
                    c2 = __fxcpmadd(c2,a2,scale);
                    __stfpd(&c[i  ],c0);
                    __stfpd(&c[i+2],c2);
                }
            }
        }
        t1 = getticks();
        dt2 = t1 - t0;

        // VERIFICATION
        for (i=0;i<dim;i++)
        {
            if (b[i] != c[i])
            {
                printf("%4d %30.15f %30.15f\n",i,b[i],c[i]);
            }
        }

        fprintf(stderr,"INTRINSICS VERSION 3\n");

        // WARM-UP
        for (i=0;i<dim;i++) c[i] = 0.0;
        for (i=0;i<dim;i+=8)
        {
            {
                double _Complex a0, a2, a4, a6;
                double _Complex c0, c2, c4, c6;
                a0 = __lfpd(&a[i  ]);
                a2 = __lfpd(&a[i+2]);
                a4 = __lfpd(&a[i+4]);
                a6 = __lfpd(&a[i+6]);
                c0 = __lfpd(&c[i  ]);
                c2 = __lfpd(&c[i+2]);
                c4 = __lfpd(&c[i+4]);
                c6 = __lfpd(&c[i+6]);
                c0 = __fxcpmadd(c0,a0,scale);
                c2 = __fxcpmadd(c2,a2,scale);
                c4 = __fxcpmadd(c4,a4,scale);
                c6 = __fxcpmadd(c6,a6,scale);
                __stfpd(&c[i  ],c0);
                __stfpd(&c[i+2],c2);
                __stfpd(&c[i+4],c4);
                __stfpd(&c[i+6],c6);
            }
        }

        // TIMING
        for (i=0;i<dim;i++) c[i] = 0.0;
        t0 = getticks();
        for (j=0;j<count;j++)
        {
            for (i=0;i<dim;i+=8)
            {
                {
                    double _Complex a0, a2, a4, a6;
                    double _Complex c0, c2, c4, c6;
                    a0 = __lfpd(&a[i  ]);
                    a2 = __lfpd(&a[i+2]);
                    a4 = __lfpd(&a[i+4]);
                    a6 = __lfpd(&a[i+6]);
                    c0 = __lfpd(&c[i  ]);
                    c2 = __lfpd(&c[i+2]);
                    c4 = __lfpd(&c[i+4]);
                    c6 = __lfpd(&c[i+6]);
                    c0 = __fxcpmadd(c0,a0,scale);
                    c2 = __fxcpmadd(c2,a2,scale);
                    c4 = __fxcpmadd(c4,a4,scale);
                    c6 = __fxcpmadd(c6,a6,scale);
                    __stfpd(&c[i  ],c0);
                    __stfpd(&c[i+2],c2);
                    __stfpd(&c[i+4],c4);
                    __stfpd(&c[i+6],c6);
                }
            }
        }
        t1 = getticks();
        dt3 = t1 - t0;

        // VERIFICATION
        for (i=0;i<dim;i++)
        {
            if (b[i] != c[i])
            {
                printf("%4d %30.15f %30.15f\n",i,b[i],c[i]);
            }
        }

        fprintf(stderr,"INTRINSICS VERSION 4\n");

        // WARM-UP
        for (i=0;i<dim;i++) c[i] = 0.0;
        for (i=0;i<dim;i+=16)
        {
            {
                double _Complex a0, a2, a4, a6, a8, a10, a12, a14;
                double _Complex c0, c2, c4, c6, c8, c10, c12, c14;
                a0 = __lfpd(&a[i   ]);
                a2 = __lfpd(&a[i+ 2]);
                a4 = __lfpd(&a[i+ 4]);
                a6 = __lfpd(&a[i+ 6]);
                a4 = __lfpd(&a[i+ 8]);
                a6 = __lfpd(&a[i+10]);
                a4 = __lfpd(&a[i+12]);
                a6 = __lfpd(&a[i+14]);
                c0 = __lfpd(&c[i   ]);
                c2 = __lfpd(&c[i+ 2]);
                c4 = __lfpd(&c[i+ 4]);
                c6 = __lfpd(&c[i+ 6]);
                c4 = __lfpd(&c[i+ 8]);
                c6 = __lfpd(&c[i+10]);
                c4 = __lfpd(&c[i+12]);
                c6 = __lfpd(&c[i+14]);
                c0 = __fxcpmadd( c0, a0,scale);
                c2 = __fxcpmadd( c2, a2,scale);
                c4 = __fxcpmadd( c4, a4,scale);
                c6 = __fxcpmadd( c6, a6,scale);
                c4 = __fxcpmadd( c8, a8,scale);
                c6 = __fxcpmadd(c10,a10,scale);
                c4 = __fxcpmadd(c12,a12,scale);
                c6 = __fxcpmadd(c14,a14,scale);
                __stfpd(&c[i   ],c0);
                __stfpd(&c[i+ 2],c2);
                __stfpd(&c[i+ 4],c4);
                __stfpd(&c[i+ 6],c6);
                __stfpd(&c[i+ 8],c4);
                __stfpd(&c[i+10],c6);
                __stfpd(&c[i+12],c4);
                __stfpd(&c[i+14],c6);
            }
        }

        // TIMING
        for (i=0;i<dim;i++) c[i] = 0.0;
        t0 = getticks();
        for (j=0;j<count;j++)
        {
            for (i=0;i<dim;i+=16)
            {
                {
                    double _Complex a0, a2, a4, a6, a8, a10, a12, a14;
                    double _Complex c0, c2, c4, c6, c8, c10, c12, c14;
                    a0  = __lfpd(&a[i   ]);
                    a2  = __lfpd(&a[i+ 2]);
                    a4  = __lfpd(&a[i+ 4]);
                    a6  = __lfpd(&a[i+ 6]);
                    a8  = __lfpd(&a[i+ 8]);
                    a10 = __lfpd(&a[i+10]);
                    a12 = __lfpd(&a[i+12]);
                    a14 = __lfpd(&a[i+14]);
                    c0  = __lfpd(&c[i   ]);
                    c2  = __lfpd(&c[i+ 2]);
                    c4  = __lfpd(&c[i+ 4]);
                    c6  = __lfpd(&c[i+ 6]);
                    c8  = __lfpd(&c[i+ 8]);
                    c10 = __lfpd(&c[i+10]);
                    c12 = __lfpd(&c[i+12]);
                    c14 = __lfpd(&c[i+14]);
                    c0  = __fxcpmadd( c0, a0,scale);
                    c2  = __fxcpmadd( c2, a2,scale);
                    c4  = __fxcpmadd( c4, a4,scale);
                    c6  = __fxcpmadd( c6, a6,scale);
                    c8  = __fxcpmadd( c8, a8,scale);
                    c10 = __fxcpmadd(c10,a10,scale);
                    c12 = __fxcpmadd(c12,a12,scale);
                    c14 = __fxcpmadd(c14,a14,scale);
                    __stfpd(&c[i   ], c0);
                    __stfpd(&c[i+ 2], c2);
                    __stfpd(&c[i+ 4], c4);
                    __stfpd(&c[i+ 6], c6);
                    __stfpd(&c[i+ 8], c8);
                    __stfpd(&c[i+10],c10);
                    __stfpd(&c[i+12],c12);
                    __stfpd(&c[i+14],c14);
                }
            }
        }
        t1 = getticks();
        dt4 = t1 - t0;

        // VERIFICATION
        for (i=0;i<dim;i++)
        {
            if (b[i] != c[i])
            {
                printf("%4d %30.15f %30.15f\n",i,b[i],c[i]);
            }
        }

        printf("%18d %18llu %18llu %18llu %18llu %18llu\n",dim,dt0,dt1,dt2,dt3,dt4);

        free(a);
        free(b);
        free(c);

    }

    fprintf(stderr,"ALL DONE\n");

    return(0);
}
