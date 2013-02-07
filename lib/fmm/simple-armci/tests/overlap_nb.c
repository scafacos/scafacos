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
 * overlap_nb:                                                             *
 *       -test of overlapping communication and computation for Vitali     *
 *                                                                         *
 ***************************************************************************/

#define REPS 28

unsigned long long int DCMF_Timebase();
void delay( unsigned long long delay_ticks );

int overlap_nb(int me, int nproc, int len)
{
    int status;
    int i;
#ifdef DEBUG
    int n;
#endif
    unsigned long long int delays[ REPS ];
    unsigned long long int t0, t1;
    unsigned long long int cp, cm, tt;
    double ov;
    armci_hdl_t nb_handle;

    /* setup delays */
    delays[0] = 0;
    for ( i = 1; i < REPS; i++ )
    {
        delays[i] = pow(2,i) - 1;
    }

    /* register remote pointers */
    double** addr_vec1 = (double **) malloc( nproc * sizeof(double *) );
    double** addr_vec2 = (double **) malloc( nproc * sizeof(double *) );
    ARMCI_Malloc((void **) addr_vec1, len*sizeof(double));
    ARMCI_Malloc((void **) addr_vec2, len*sizeof(double));
    MPI_Barrier(MPI_COMM_WORLD);

    /* initialization of local segments */
    for( i=0 ; i<len ; i++ ){
       addr_vec1[me][i] = (double) +1*(1000*me+i);    
    }
    for( i=0 ; i<len ; i++ ){
       addr_vec2[me][i] = (double) -1*(1000*me+i);    
    }

#ifdef DEBUG
    /* print before exchange */
    for( n=0 ; n<nproc ; n++){
       MPI_Barrier(MPI_COMM_WORLD);
       if (n==me){
          printf("values before exchange\n");
          for( i=0 ; i<len ; i++ ){
             printf("proc %d: addr_vec1[%d][%d] = %f\n", n, n, i, addr_vec1[n][i]);
          }
          for( i=0 ; i<len ; i++ ){
             printf("proc %d: addr_vec2[%d][%d] = %f\n", n, n, i, addr_vec2[n][i]);
          }
          fflush(stdout);
       }
       MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    for ( i = 0; i < REPS; i++ ){

        MPI_Barrier(MPI_COMM_WORLD);
        ARMCI_INIT_HANDLE(&nb_handle);

        t0 = DCMF_Timebase();

        if (me == 0){

           status = ARMCI_NbGet(addr_vec1[1], addr_vec1[0], len*sizeof(double), 1, &nb_handle);
           delay( delays[i] );

        } else if (me == 1){

           status = ARMCI_NbGet(addr_vec2[0], addr_vec2[1], len*sizeof(double), 0, &nb_handle);
           delay( delays[i] );

        }

        if((status != 0) && (me == 0)) printf("%s: ARMCI_Get failed at line %d\n",__FILE__,__LINE__);

        ARMCI_Wait(&nb_handle);
        MPI_Barrier(MPI_COMM_WORLD);

        t1 = DCMF_Timebase();

        if (me == 0){
           //printf("Iter %6d Proc %6d: (t0,t1) = %16lld %16lld\n",i,me,t0,t1);
           tt = t1 - t0;
           cp = delays[i];
           cm = tt - cp;
           ov = (double)cp / (double)(tt);
           printf("NONBLOCK %6d: comp, comm, total, ratio:  %16lld  %16lld  %16lld  %18.8lf\n", me, cp, cm, tt, ov );
        }
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
    /* print after exchange */
    for( n=0 ; n<nproc ; n++){
       MPI_Barrier(MPI_COMM_WORLD);
       if (n==me){
          printf("values after exchange\n");
          for( i=0 ; i<len ; i++ ){
             printf("proc %d: addr_vec1[%d][%d] = %f\n", n, n, i, addr_vec1[n][i]);
          }
          for( i=0 ; i<len ; i++ ){
             printf("proc %d: addr_vec2[%d][%d] = %f\n", n, n, i, addr_vec2[n][i]);
          }
          fflush(stdout);
       }
       MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    return(0);
}
