/***************************************************************************
 *   Copyright (C) 2009 by Jeff Hammond                                    *
 *   jeff.science@gmail.com                                                *
 *                                                                         *
 * Redistribution and use in source and binary forms, with or without      *
 * modification, are permitted provided that the following conditions      *
 * are met:                                                                *
 * 1. Redistributions of source code must retain the above copyright       *
 *    notice, this list of conditions and the following disclaimer.        *
 * 2. Redistributions in binary form must reproduce the above copyright    *
 *    notice, this list of conditions and the following disclaimer in the  *
 *    documentation and/or other materials provided with the distribution. *
 * 3. The name of the author may not be used to endorse or promote         *
 *    products derived from this software without specific prior written   *
 *    permission.                                                          *
 *                                                                         *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR    *
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED          *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE  *
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,      *
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR      *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)      *
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,     *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING   *
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE      *
 * POSSIBILITY OF SUCH DAMAGE.                                             *
 *                                                                         *
 ***************************************************************************/

#include "driver.h"

/***************************************************************************
 *                                                                         *
 * simple_get:                                                             *
 *       -demonstrates how to allocate some shared segements with ARMCI    *
 *       -demonstrates how to do one-sided point-to-point communication    *
 *       -inspired by armci/examples/features/non-blocking/simple/simple.c *
 *                                                                         *
 ***************************************************************************/

int simple_get(int me, int nproc, int len)
{
    int status;
    int n,i;
    double t0,t1;

    double** addr_vec = (double **) malloc( nproc * sizeof(double *) );
    ARMCI_Malloc((void **) addr_vec, len*sizeof(double));
    MPI_Barrier(MPI_COMM_WORLD);

    /* initialization of local segments */
    for( i=0 ; i<len ; i++ ){
       addr_vec[me][i] = (double) (1000*me+i);    
    }

    /* print before exchange */
    for( n=0 ; n<nproc ; n++){
       MPI_Barrier(MPI_COMM_WORLD);
       if (n==me){
          printf("values before exchange\n");
          for( i=0 ; i<len ; i++ ){
             printf("proc %d: addr_vec[%d][%d] = %f\n", n, n, i, addr_vec[n][i]);
          }
          fflush(stdout);
       }
       MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    /* even processes get from odd right neighbors */
    if (me%2 == 0){
       t0 = MPI_Wtime();
       status = ARMCI_Get(addr_vec[me+1], addr_vec[me], len*sizeof(double), me+1);
       t1 = MPI_Wtime();
       if(status != 0){
    	  if (me == 0) printf("%s: ARMCI_Get failed at line %d\n",__FILE__,__LINE__);
       }
       printf("Proc %d: Get Latency=%lf microseconds\n",me,1e6*(t1-t0)/len);
       fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    /* print after exchange */
    for( n=0 ; n<nproc ; n++){
       MPI_Barrier(MPI_COMM_WORLD);
       if (n==me){
          printf("values after exchange\n");
          for( i=0 ; i<len ; i++ ){
             printf("proc %d: addr_vec[%d][%d] = %f\n", n, n, i, addr_vec[n][i]);
          }
          fflush(stdout);
       }
       MPI_Barrier(MPI_COMM_WORLD);
    }

    return(0);
}
