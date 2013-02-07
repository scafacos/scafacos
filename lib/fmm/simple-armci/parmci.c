/********************************************************************
 * The following is a notice of limited availability of the code, and disclaimer
 * which must be included in the prologue of the code and in all source listings
 * of the code.
 *
 * Copyright (c) 2010 Argonne Leadership Computing Facility, Argonne National Laboratory
 *
 * Permission is hereby granted to use, reproduce, prepare derivative works, and
 * to redistribute to others.
 *
 *                 LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer listed
 *    in this license in the documentation and/or other materials
 *    provided with the distribution.
 *
 *  - Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * The copyright holders provide no reassurances that the source code
 * provided does not infringe any patent, copyright, or any other
 * intellectual property rights of third parties.  The copyright holders
 * disclaim any liability to any recipient for claims brought against
 * recipient by any third party for infringement of that parties
 * intellectual property rights.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************/

#include "parmci.h"

/* initialization and termination */

int PARMCI_Init(void)
{
    return A1D_Initialize();
}

int PARMCI_Init_args(int *argc, char ***argv)
{
    //fprintf(stderr,"PARMCI_Init_args: argc/argv may not be setup properly by device \n");
    return A1D_Initialize();
}

void PARMCI_Finalize(void)
{
    A1D_Finalize();
    return;
}

/* memory management */

void * PARMCI_Malloc_local(int bytes)
{
    return A1D_Allocate_local(bytes);
}

int PARMCI_Free_local(void * ptr)
{
    A1D_Free_local(ptr);
    return(0);
}

int PARMCI_Malloc(void * ptr_arr[], int bytes)
{
    return A1D_Allocate_comm(MPI_COMM_WORLD, ptr_arr, bytes);
}

int PARMCIX_Malloc_comm(MPI_Comm comm, void * ptr_arr[], int bytes)
{
    return A1D_Allocate_comm(comm, ptr_arr, bytes);
}


int PARMCI_Free(void * ptr)
{
    A1D_Free_comm(MPI_COMM_WORLD, ptr);
    return(0);
}

int PARMCIX_Free_comm(MPI_Comm comm, void * ptr)
{
    A1D_Free_comm(comm, ptr);
    return(0);
}

/* NOT USED
void * PARMCI_Memat(armci_meminfo_t* meminfo, int memflg)
{
    fprintf(stderr,"PARMCI_Memat: not implemented \n");
    assert(0);
    return (void *) NULL;
}

void PARMCI_Memget(size_t bytes, armci_meminfo_t* meminfo, int memflg)
{
    fprintf(stderr,"PARMCI_Memat: not implemented \n");
    assert(0);
    return;
}
 */

/* synchronization */

void PARMCI_Barrier(void)
{
    A1D_Flush_all();
    MPI_Barrier(MPI_COMM_WORLD);
    return;
}

void PARMCI_Fence(int proc)
{
    A1D_Flush(proc);
    return;
}
void PARMCI_AllFence(void)
{
    A1D_Flush_all();
    return;
}

void PARMCIX_Barrier_comm(MPI_Comm comm)
{
    int comm_result;

    MPI_Comm_compare( MPI_COMM_WORLD, comm, &comm_result );

    if ( comm_result==MPI_IDENT || comm_result==MPI_CONGRUENT )
        A1D_Flush_all();
    else
        A1D_Flush_comm(comm);

    MPI_Barrier(comm);

    return;
}

void PARMCIX_AllFence_comm(MPI_Comm comm)
{
    int comm_result;

    MPI_Comm_compare( MPI_COMM_WORLD, comm, &comm_result );

    if ( comm_result==MPI_IDENT || comm_result==MPI_CONGRUENT )
        A1D_Flush_all();
    else
        A1D_Flush_comm(comm);

    return;
}

int PARMCI_Test(armci_hdl_t * nb_handle)
{
    int status = -1;
    A1D_Test(nb_handle->a1d_handle,&status);
    /* ARMCI: 0=complete, 1=in-progress */
    /* A1D:   1=complete, 0=incomplete */
    if      (status==1) return 0;
    else if (status==0) return 1;
    else                return -1;
}

int PARMCI_Wait(armci_hdl_t * nb_handle)
{
    return A1D_Wait(nb_handle->a1d_handle);
}

int PARMCI_WaitProc(int proc)
{
    fprintf(stderr,"WARNING: PARMCI_WaitProc(int proc) only synchronizes implicit nonblocking operations! \n");
    PARMCI_Fence(proc);
    return(0);
}

int PARMCI_WaitAll(void)
{
    fprintf(stderr,"WARNING: PARMCI_WaitAll() only synchronizes implicit nonblocking operations! \n");
    A1D_Flush_all();
    return(0);
}

/* remote atomic update and mutexes */

#if defined(USE_32B_ATOMICS) && defined(USE_64B_ATOMICS)
#error The call syntax of (P)ARMCI_Rmw is stupid.  Use 32B or 64B atomics but not both.
#endif

#if !defined(USE_32B_ATOMICS) && !defined(USE_64B_ATOMICS)
#error You must define atomics to be 32B or 64B!
#endif

#if defined(USE_32B_ATOMICS)
int32_t PARMCI_Rmw(int optype, int32_t * local, int32_t * remote, int32_t incr, int proc)
#elif defined(USE_64B_ATOMICS)
int64_t PARMCI_Rmw(int optype, int64_t * local, int64_t * remote, int64_t incr, int proc)
#endif
{
#if defined(USE_32B_ATOMICS)
    int32_t copy;
#elif defined(USE_64B_ATOMICS)
    int64_t copy;
#else
#warning no atomics available
#endif

    switch (optype)
    {
#ifdef USE_32B_ATOMICS
        case ARMCI_FETCH:
            A1D_Fetch_and_inc32(proc, remote, local, (int32_t)0 );
            return *local;

        case ARMCI_FETCH_AND_ADD:
            A1D_Fetch_and_inc32(proc, remote, local, incr );
            return *local;
#elif defined(USE_64B_ATOMICS)
        case ARMCI_FETCH_LONG:
            A1D_Fetch_and_inc64(proc, remote, local, (int64_t)0 );
            return *local;

        case ARMCI_FETCH_AND_ADD_LONG:
            A1D_Fetch_and_inc64(proc, remote, local, incr );
            return *local;
#endif
        default:
            fprintf(stderr,"PARMCI_Rmw: unknown operation request! \n");
            assert(0);
            return -1;
    }
    return -1;
}

int PARMCI_Create_mutexes(int num)
{
    fprintf(stderr,"PARMCI_Create_mutexes: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_Destroy_mutexes(void)
{
    fprintf(stderr,"PARMCI_Destroy_mutexes: not implemented \n");
    assert(0);
    return(-1);
}

void PARMCI_Lock(int mutex, int proc)
{
    fprintf(stderr,"PARMCI_Lock: not implemented \n");
    assert(0);
    return;
}

void PARMCI_Unlock(int mutex, int proc)
{
    fprintf(stderr,"PARMCI_Unlock: not implemented \n");
    assert(0);
    return;
}

/* blocking one-sided */

int PARMCI_Get(void *src, void *dst, int bytes, int proc)
{
    return A1D_GetC(proc, bytes, src, dst);
}

int PARMCI_Put(void *src, void *dst, int bytes, int proc)
{
    return A1D_PutC(proc, bytes, src, dst);
}

int PARMCI_Acc(armci_acc_t type, void *scale, void *src, void* dst, int bytes, int proc)
{
    switch(type)
    {
        case ARMCI_ACC_DBL:
            return A1D_AccC(proc, bytes, src, dst, A1D_DOUBLE, scale);
            break;

        case ARMCI_ACC_FLT:
            return A1D_AccC(proc, bytes, src, dst, A1D_SINGLE, scale);
            break;

#ifdef A1D_USE_COMPLEX
        case ARMCI_ACC_DCP:
            return A1D_AccC(proc, bytes, src, dst, A1D_DOUBLE_COMPLEX, scale);
            break;

        case ARMCI_ACC_CPL:
            return A1D_AccC(proc, bytes, src, dst, A1D_SINGLE_COMPLEX, scale);
            break;

#endif
        case ARMCI_ACC_INT:
            return A1D_AccC(proc, bytes, src, dst, A1D_INT32, scale);
            break;

        case ARMCI_ACC_LNG:
            return A1D_AccC(proc, bytes, src, dst, A1D_INT64, scale);
            break;
    }
    return -1;
}

int PARMCI_GetS(void *src_ptr, int *src_stride_arr,
                void *dst_ptr, int *dst_stride_arr,
                int *block_sizes, int stride_levels, int proc)
{
    fprintf(stderr,"PARMCI_GetS: not implemented \n");
    assert(0);
    return(-1);

    //    return A1D_GetS(proc, stride_levels, block_sizes,
    //                    src_ptr, src_stride_arr,
    //                    dst_ptr, dst_stride_arr);
}

int PARMCI_PutS(void *src_ptr, int *src_stride_arr,
                void *dst_ptr, int *dst_stride_arr,
                int *block_sizes, int stride_levels, int proc)
{
    fprintf(stderr,"PARMCI_PutS: not implemented \n");
    assert(0);
    return(-1);

    //    return A1D_PutS(proc, stride_levels, block_sizes,
    //                    src_ptr, src_stride_arr,
    //                    dst_ptr, dst_stride_arr);
}

int PARMCI_AccS(armci_acc_t optype, void *scale,
                void *src_ptr, int *src_stride_arr,
                void *dst_ptr, int *dst_stride_arr,
                int *count, int stride_levels, int proc)
{
    fprintf(stderr,"PARMCI_AccS: not implemented \n");
    assert(0);
    return(-1);

    //    return A1D_AccS(proc, stride_levels, block_sizes,
    //                    src_ptr, src_stride_arr,
    //                    dst_ptr, dst_stride_arr,
    //                    type, scale);
}


int PARMCI_GetV(armci_giov_t * array_descr, int len, int proc)
{
    fprintf(stderr,"PARMCI_GetV: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_PutV(armci_giov_t * array_descr, int len, int proc)
{
    fprintf(stderr,"PARMCI_PutV: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_AccV(armci_acc_t type, void *scale, armci_giov_t * array_descr, int len, int proc)
{
    fprintf(stderr,"PARMCI_AccV: not implemented \n");
    assert(0);
    return(-1);
}

/* non-blocking one-sided */

int PARMCI_NbGet(void *src, void *dst, int bytes, int proc,
                 armci_hdl_t * nb_handle)
{
    return A1D_GetC(proc, bytes, src, dst);
}

int PARMCI_NbPut(void *src, void *dst, int bytes, int proc,
                 armci_hdl_t * nb_handle)
{
    return A1D_PutC(proc, bytes, src, dst);
}

int PARMCI_NbAcc(int type, void *scale, void *src, void* dst, int bytes, int proc,
                 armci_hdl_t * nb_handle)
{
    return A1D_AccC(proc, bytes, src, dst, type, scale);
}

int PARMCI_NbGetS(void *src_ptr, int *src_stride_arr,
                  void *dst_ptr, int *dst_stride_arr,
                  int *block_sizes, int stride_levels, int proc,
                  armci_hdl_t * nb_handle)
{
    fprintf(stderr,"PARMCI_NbGetS: not implemented \n");
    assert(0);
    return(-1);

    //    return A1D_GetS(proc, stride_levels, block_sizes,
    //                    src_ptr, src_stride_arr,
    //                    dst_ptr, dst_stride_arr);
}

int PARMCI_NbPutS(void *src_ptr, int *src_stride_arr,
                  void *dst_ptr, int *dst_stride_arr,
                  int *block_sizes, int stride_levels, int proc,
                  armci_hdl_t * nb_handle)
{
    fprintf(stderr,"PARMCI_NbPutS: not implemented \n");
    assert(0);
    return(-1);

    //    return A1D_PutS(proc, stride_levels, block_sizes,
    //                    src_ptr, src_stride_arr,
    //                    dst_ptr, dst_stride_arr);
}

int PARMCI_NbAccS(int optype, void *scale,
                  void *src_ptr, int *src_stride_arr,
                  void *dst_ptr, int *dst_stride_arr,
                  int *count, int stride_levels, int proc,
                  armci_hdl_t * nb_handle)
{
    fprintf(stderr,"PARMCI_NbAccS: not implemented \n");
    assert(0);
    return(-1);

    //    return A1D_AccS(proc, stride_levels, block_sizes,
    //                    src_ptr, src_stride_arr,
    //                    dst_ptr, dst_stride_arr,
    //                    type, scale);
}


int PARMCI_NbGetV(armci_giov_t * array_descr, int len, int proc,
                  armci_hdl_t * nb_handle)
{
    fprintf(stderr,"PARMCI_NbGetV: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_NbPutV(armci_giov_t * array_descr, int len, int proc,
                  armci_hdl_t * nb_handle)
{
    fprintf(stderr,"PARMCI_NbPutV: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_NbAccV(int type, void *scale, armci_giov_t * array_descr, int len, int proc,
                  armci_hdl_t * nb_handle)
{
    fprintf(stderr,"PARMCI_NbAccV: not implemented \n");
    assert(0);
    return(-1);
}

/* ??? extensions */

int PARMCI_Put_flag(void *src_ptr, void *dst_ptr, int bytes, int *flag, int val, int proc)
{
    fprintf(stderr,"PARMCI_Put_flag: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_PutS_flag(void *src_ptr, int *src_stride_arr,
                     void *dst_ptr, int *dst_stride_arr,
                     int *count, int stride_levels,
                     int *flag, int val, int proc)
{
    fprintf(stderr,"PARMCI_PutS_flag: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_PutS_flag_dir(void *src_ptr, int *src_stride_arr,
                         void *dst_ptr, int *dst_stride_arr,
                         int *count, int stride_levels,
                         int *flag, int val, int proc)
{
    fprintf(stderr,"PARMCI_PutS_flag_dir: not implemented \n");
    assert(0);
    return(-1);
}

/* CAF extensions */

int PARMCI_PutValueInt(int src, void *dst, int proc)
{
    fprintf(stderr,"PARMCI_PutValueInt: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_PutValueLong(long src, void *dst, int proc)
{
    fprintf(stderr,"PARMCI_PutValueLong: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_PutValueFloat(float src, void *dst, int proc)
{
    fprintf(stderr,"PARMCI_PutValueFloat: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_PutValueDouble(double src, void *dst, int proc)
{
    fprintf(stderr,"PARMCI_PutValueDouble: not implemented \n");
    assert(0);
    return(-1);
}

int PARMCI_GetValueInt(void *src, int proc)
{
    fprintf(stderr,"PARMCI_GetValueInt: not implemented \n");
    assert(0);
    return(-1);
}

long PARMCI_GetValueLong(void *src, int proc)
{
    fprintf(stderr,"PARMCI_GetValueLong: not implemented \n");
    assert(0);
    return(-1);
}

float PARMCI_GetValueFloat(void *src, int proc)
{
    fprintf(stderr,"PARMCI_GetValueFloat: not implemented \n");
    assert(0);
    return(-1);
}

double PARMCI_GetValueDouble(void *src, int proc)
{
    fprintf(stderr,"PARMCI_GetValueDouble: not implemented \n");
    assert(0);
    return(-1);
}

/**********************************************/
