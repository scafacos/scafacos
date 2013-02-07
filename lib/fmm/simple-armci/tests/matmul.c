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

int main(int argc, char **argv)
{
	int me,nproc;
    int status;
    int rank;

    /* initialization */
    MPI_Init(&argc, &argv);
    ARMCI_Init();

#ifdef HPC_PROFILING
    HPM_Init();
#endif

    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

#ifdef DEBUG
    if(me == 0){
       printf("The result of MPI_Comm_size is %d\n",nproc);
       fflush(stdout);
    }
#endif

    /* get the matrix parameters */
    if (argc > 1){
        rank = atoi(argv[1]);
    } else {
        rank = 8;
    }
    if (me == 0){
        printf("Running matmul.x with rank = %d\n",rank);
        fflush(stdout);
    }

    /* register remote pointers */
    double** addr_A = (double **) ARMCI_Malloc_local(sizeof(double *) * nproc);
    if (addr_A == NULL) ARMCI_Error("malloc A failed at line",0);

    double** addr_B = (double **) ARMCI_Malloc_local(sizeof(double *) * nproc);
    if (addr_B == NULL) ARMCI_Error("malloc B failed at line",0);

    double** addr_C = (double **) ARMCI_Malloc_local(sizeof(double *) * nproc);
    if (addr_C == NULL) ARMCI_Error("malloc C failed at line",0);

#ifdef DEBUG
    if(me == 0) printf("ARMCI_Malloc A requests %lu bytes\n",rank*rank*sizeof(double));
    fflush(stdout);
#endif
    status = ARMCI_Malloc((void **) addr_A, rank*rank*sizeof(double));
    if (status != 0) ARMCI_Error("ARMCI_Malloc A failed",status);

#ifdef DEBUG
    if(me == 0) printf("ARMCI_Malloc B requests %lu bytes\n",rank*rank*sizeof(double));
    fflush(stdout);
#endif
    status = ARMCI_Malloc((void **) addr_B, rank*rank*sizeof(double));
    if (status != 0) ARMCI_Error("ARMCI_Malloc B failed",status);

#ifdef DEBUG
    if(me == 0) printf("ARMCI_Malloc C requests %lu bytes\n",rank*rank*sizeof(double));
    fflush(stdout);
#endif
    status = ARMCI_Malloc((void **) addr_C, rank*rank*sizeof(double));
    if (status != 0) ARMCI_Error("ARMCI_Malloc C failed",status);

    MPI_Barrier(MPI_COMM_WORLD);

    /* free ARMCI pointers */
    ARMCI_Free_local(addr_C);
    ARMCI_Free_local(addr_B);
    ARMCI_Free_local(addr_A);

#ifdef HPC_PROFILING
    HPM_Print();
#endif

    /* the end */
    ARMCI_Finalize();
    MPI_Finalize();

    return(0);
}



