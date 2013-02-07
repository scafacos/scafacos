/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "armci.h"
#include "a1.h"
#include "a1d.h"
#include "a1u.h"
#include "assert.h"
#include "mpi.h"

/* #define A1_ARMCI_PROFILING */

#ifdef A1_ARMCI_PROFILING
int __a1_prof_me = -1;
char* __a1_prof_name[32];
double __a1_prof_t0;
double __a1_prof_t1;
double __a1_prof_dt;

#define AAP_INIT()                                    \
        do {                                              \
            __a1_prof_me = A1_Process_id(A1_GROUP_WORLD); \
        } while (0)

#define AAP_START(a)                      \
        do {                                  \
            __a1_prof_name[0] = a;            \
            __a1_prof_t0 = A1_Time_seconds(); \
        } while (0)

#define AAP_STOP()                                                                      \
        do {                                                                                \
            __a1_prof_t1 = A1_Time_seconds();                                               \
            __a1_prof_dt = __a1_prof_t1 - __a1_prof_t0;                                     \
            printf("iam %d: %s took %10.4lf s\n",__a1_prof_me,__a1_prof_name,__a1_prof_dt); \
            fflush(stdout);                                                                 \
        } while (0)

#define AAP_ARGS printf

#else
#define AAP_INIT()
#define AAP_START(a)
#define AAP_STOP(a)
#define AAP_ARGS(...)
#endif

/* TODO We really should have ARMCI-to-A1 type/op conversion functions/macros
 *      instead of repeating that code so many times */

int ARMCI_Init_args(int *argc, char ***argv)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    AAP_INIT();

    status = A1_Initialize(A1_THREAD_SINGLE);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Initialize returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;

}

int ARMCI_Init(void)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    AAP_INIT();

    status = A1_Initialize(A1_THREAD_SINGLE);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Initialize returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;

}

int ARMCI_Finalize(void)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1_Finalize();
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Finalize returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;

}

int ARMCI_Malloc(void* ptr[], armci_size_t bytes)
{
    int status = A1_SUCCESS;
    int my_rank;
    int my_size = sizeof(void*)

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    my_rank = A1_Process_id(A1_GROUP_WORLD);

    if (bytes == 0)
    {
        ptr[my_rank] = NULL;
    }
    else
    {
        status = A1_Alloc_segment(&ptr[my_rank], bytes);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Alloc_segment returned an error\n");
    }

    status = MPI_Allgather(MPI_IN_PLACE,my_size,MPI_BYTE,ptr,my_size,MPI_BYTE,MPI_COMM_WORLD);
    A1U_ERR_POP(status != MPI_SUCCESS,"MPI_Allgather failed");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;

}

void* ARMCI_Malloc_local(armci_size_t bytes)
{
    int status = A1_SUCCESS;
    void *segment_ptr;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1_Alloc_segment(&segment_ptr, bytes);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Alloc_segement returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return segment_ptr;

    fn_fail: goto fn_exit;
}

int ARMCI_Free(void *ptr)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1_Release_segments(A1_GROUP_WORLD, ptr);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Release_segments returned an error\n");

    status = A1_Free_segment(ptr);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Free_segment returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_Free_local(void *ptr)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1_Free_segment(ptr);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Free_segment returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void ARMCI_INIT_HANDLE(armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1_Allocate_handle(&a1_handle);
    A1U_ERR_ABORT(status != A1_SUCCESS,
            "A1_Allocate_handle returned an error\n");

    *handle = a1_handle;

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

int ARMCI_Put(void* src, void* dst, int bytes, int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    AAP_ARGS("iam %d: A1_Put proc = %d, bytes = %d\n",__a1_prof_me,proc,bytes);AAP_START("A1_Put          ");

    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }

    status = A1_Put(proc, src, dst, bytes);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Put returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_NbPut(void* src, void* dst, int bytes, int proc, armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    AAP_ARGS("iam %d: A1_NbPut proc = %d, levels = %d, count[0] = %d, count[1] = %d\n",__a1_prof_me,proc,stride_levels,count[0],count[stride_levels-1]);
    AAP_START("A1_NbPutS          ");
    status = A1_NbPut(proc, src, dst, bytes, a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_NbPut returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_PutS(void* src_ptr,
               int src_stride_ar[],
               void* dst_ptr,
               int dst_stride_ar[],
               int count[],
               int stride_levels,
               int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    AAP_ARGS("iam %d: A1_PutS proc = %d, levels = %d, count[0] = %d, count[1] = %d\n",__a1_prof_me,proc,stride_levels,count[0],count[stride_levels-1]);
    AAP_START("A1_PutS          ");

    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }

    status = A1_PutS(proc,
                     stride_levels,
                     count,
                     src_ptr,
                     src_stride_ar,
                     dst_ptr,
                     dst_stride_ar);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_PutS returned an error\n");AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_NbPutS(void* src_ptr,
                 int src_stride_ar[],
                 void* dst_ptr,
                 int dst_stride_ar[],
                 int count[],
                 int stride_levels,
                 int proc,
                 armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    AAP_ARGS("iam %d: A1_NbPutS proc = %d, levels = %d, count[0] = %d, count[1] = %d\n",__a1_prof_me,proc,stride_levels,count[0],count[stride_levels-1]);
    AAP_START("A1_NbPutS          ");
    status = A1_NbPutS(proc,
                       stride_levels,
                       count,
                       src_ptr,
                       src_stride_ar,
                       dst_ptr,
                       dst_stride_ar,
                       a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_NbPutS returned an error\n");AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_PutV(armci_giov_t *dsrc_arr, int arr_len, int proc)
{
    int status = A1_SUCCESS;
    A1_iov_t *a1_iov_ar;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    /* ARMCI iov and A1 iov are similar structures but follow 
     * different naming conventions. So we make a copy.*/
    /* TODO Why not use A1D_Malloc here? A1D_Malloc is not exposed here, currently*/
    status = posix_memalign((void **) &a1_iov_ar, 16, sizeof(A1_iov_t)
            * arr_len);
    A1U_ERR_POP(status != 0, "posix_memalign returned an error\n");

    memcpy((void *) a1_iov_ar, (void *) dsrc_arr, sizeof(A1_iov_t) * arr_len);

    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }

    status = A1_PutV(proc, a1_iov_ar, arr_len);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_PutV returned an error\n");

    free(a1_iov_ar);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_NbPutV(armci_giov_t *dsrc_arr,
                 int arr_len,
                 int proc,
                 armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_iov_t *a1_iov_ar;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    /* ARMCI iov and A1 iov are similar structures but follow
     * different naming conventions. So we make a copy.*/
    /* TODO Why not use A1D_Malloc here? A1D_Malloc is not exposed here, currently*/
    status = posix_memalign((void **) &a1_iov_ar, 16, sizeof(A1_iov_t)
            * arr_len);
    A1U_ERR_POP(status != 0, "posix_memalign returned an error\n");

    memcpy((void *) a1_iov_ar, (void *) dsrc_arr, sizeof(A1_iov_t) * arr_len);

    status = A1_NbPutV(proc, a1_iov_ar, arr_len, a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_NbPutV returned an error\n");

    /* it is okay to free this as soon as the function returns?  so we copy it upon
     * entry into A1D_NbPutV? */
    free(a1_iov_ar);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_Get(void* src, void* dst, int bytes, int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    AAP_ARGS("iam %d: A1_Get proc = %d, bytes = %d\n",__a1_prof_me,proc,bytes);
    AAP_START("A1_Get         ");

    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }

    status = A1_Get(proc, src, dst, bytes);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Get returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_NbGet(void* src, void* dst, int bytes, int proc, armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    AAP_ARGS("iam %d: A1_NbGet proc = %d, levels = %d, count[0] = %d, count[1] = %d\n",__a1_prof_me,proc,stride_levels,count[0],count[stride_levels-1]);
    AAP_START("A1_NbGet           ");
    status = A1_NbGet(proc, src, dst, bytes, a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_NbGet returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_GetS(void* src_ptr,
               int src_stride_ar[],
               void* dst_ptr,
               int dst_stride_ar[],
               int count[],
               int stride_levels,
               int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    AAP_ARGS("iam %d: A1_GetS proc = %d, levels = %d, count[0] = %d, count[1] = %d\n",__a1_prof_me,proc,stride_levels,count[0],count[stride_levels-1]);
    AAP_START("A1_GetS           ");

    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }

    status = A1_GetS(proc,
                     stride_levels,
                     count,
                     src_ptr,
                     src_stride_ar,
                     dst_ptr,
                     dst_stride_ar);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_GetS returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_NbGetS(void* src_ptr,
                 int src_stride_ar[],
                 void* dst_ptr,
                 int dst_stride_ar[],
                 int count[],
                 int stride_levels,
                 int proc,
                 armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    AAP_ARGS("iam %d: A1_NbGetS proc = %d, levels = %d, count[0] = %d, count[1] = %d\n",__a1_prof_me,proc,stride_levels,count[0],count[stride_levels-1]);
    AAP_START("A1_NbGetS           ");
    status = A1_NbGetS(proc,
                       stride_levels,
                       count,
                       src_ptr,
                       src_stride_ar,
                       dst_ptr,
                       dst_stride_ar,
                       a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_NbPutS returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_GetV(armci_giov_t *dsrc_arr, int arr_len, int proc)
{
    int status = A1_SUCCESS;
    A1_iov_t *a1_iov_ar;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    /* ARMCI iov and A1 iov are similar structures but follow
     * different naming conventions. So we make a copy.*/

    /* TODO Why not use A1D_Malloc here? */
    status = posix_memalign((void **) &a1_iov_ar, 16, sizeof(A1_iov_t)
            * arr_len);
    A1U_ERR_POP(status != 0, "posix_memalign returned an error\n");

    memcpy((void *) a1_iov_ar, (void *) dsrc_arr, sizeof(A1_iov_t) * arr_len);

    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }

    status = A1_GetV(proc, a1_iov_ar, arr_len);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_GetV returned an error\n");

    free(a1_iov_ar);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_NbGetV(armci_giov_t *dsrc_arr,
                 int arr_len,
                 int proc,
                 armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_iov_t *a1_iov_ar;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    /* ARMCI iov and A1 iov are similar structures but follow
     * different naming conventions. So we make a copy.*/

    /* TODO Why not use A1D_Malloc here? */
    status = posix_memalign((void **) &a1_iov_ar, 16, sizeof(A1_iov_t)
            * arr_len);
    A1U_ERR_POP(status != 0, "posix_memalign returned an error\n");

    memcpy((void *) a1_iov_ar, (void *) dsrc_arr, sizeof(A1_iov_t) * arr_len);

    status = A1_NbGetV(proc, a1_iov_ar, arr_len, a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_NbGetV returned an error\n");

    free(a1_iov_ar);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_Acc(int datatype,
              void *scale,
              void* src,
              void* dst,
              int bytes,
              int proc)
{
    int status = A1_SUCCESS;
    A1_datatype_t a1_type;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    switch (datatype)
    {
        case ARMCI_INT:
        case ARMCI_LONG:
        case ARMCI_ACC_INT:
            a1_type = A1_INT32;
            break;
        case ARMCI_LONG_LONG:
            a1_type = A1_INT64;
            break;
        case ARMCI_FLOAT:
        case ARMCI_ACC_FLT:
            a1_type = A1_FLOAT;
            break;
        case ARMCI_DOUBLE:
        case ARMCI_ACC_DBL:
            a1_type = A1_DOUBLE;
            break;
        case ARMCI_ACC_CPL:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_CPL datatype not supported\n");
        case ARMCI_ACC_DCP:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_DCP datatype not supported\n");
        default:
            A1U_ERR_ABORT(status != A1_ERROR, "invalid datatype\n");
    }

    AAP_ARGS("iam %d: A1_PutAcc proc = %d, bytes = %d\n",__a1_prof_me,proc,bytes);
    AAP_START("A1_PutAcc             ");

    /* accumulate flushes puts before and holds the lock throughout
    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }
    */

    status = A1_PutAcc(proc, src, dst, bytes, a1_type, scale);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_PutAcc returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_NbAcc(int datatype,
                void *scale,
                void* src,
                void* dst,
                int bytes,
                int proc,
                armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_datatype_t a1_type;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    switch (datatype)
    {
        case ARMCI_INT:
        case ARMCI_LONG:
        case ARMCI_ACC_INT:
            a1_type = A1_INT32;
            break;
        case ARMCI_LONG_LONG:
            a1_type = A1_INT64;
            break;
        case ARMCI_FLOAT:
        case ARMCI_ACC_FLT:
            a1_type = A1_FLOAT;
            break;
        case ARMCI_DOUBLE:
        case ARMCI_ACC_DBL:
            a1_type = A1_DOUBLE;
            break;
        case ARMCI_ACC_CPL:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_CPL datatype not supported\n");
        case ARMCI_ACC_DCP:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_DCP datatype not supported\n");
        default:
            A1U_ERR_ABORT(status != A1_ERROR, "invalid datatype\n");
    }

    AAP_ARGS("iam %d: A1_NbPutAcc proc = %d, bytes = %d\n",__a1_prof_me,proc,bytes);
    AAP_START("A1_NbPutAcc             ");
    status = A1_NbPutAcc(proc, src, dst, bytes, a1_type, scale, a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_NbPutAcc returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_AccS(int datatype,
               void *scale,
               void* src_ptr,
               int src_stride_ar[],
               void* dst_ptr,
               int dst_stride_ar[],
               int count[],
               int stride_levels,
               int proc)
{
    int status = A1_SUCCESS;
    A1_datatype_t a1_type;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    switch (datatype)
    {
        case ARMCI_INT:
        case ARMCI_LONG:
        case ARMCI_ACC_INT:
        case ARMCI_ACC_LNG:
            a1_type = A1_INT32;
            break;
        case ARMCI_LONG_LONG:
            a1_type = A1_INT64;
            break;
        case ARMCI_FLOAT:
        case ARMCI_ACC_FLT:
            a1_type = A1_FLOAT;
            break;
        case ARMCI_DOUBLE:
        case ARMCI_ACC_DBL:
            a1_type = A1_DOUBLE;
            break;
        case ARMCI_ACC_CPL:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_CPL datatype not supported\n");
        case ARMCI_ACC_DCP:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_DCP datatype not supported\n");
        default:
            A1U_ERR_ABORT(status != A1_ERROR, "invalid datatype %d \n", datatype);
    }

    AAP_ARGS("iam %d: A1_PutAccS proc = %d, levels = %d, count[0] = %d, count[1] = %d\n",__a1_prof_me,proc,stride_levels,count[0],count[stride_levels-1]);
    AAP_START("A1_PutAccS             ");

    /* accumulate flushes puts before and holds the lock throughout
    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }
    */

    status = A1_PutAccS(proc,
                        stride_levels,
                        count,
                        src_ptr,
                        src_stride_ar,
                        dst_ptr,
                        dst_stride_ar,
                        a1_type,
                        scale);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_PutAccS returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_NbAccS(int datatype,
                 void *scale,
                 void* src_ptr,
                 int src_stride_ar[],
                 void* dst_ptr,
                 int dst_stride_ar[],
                 int count[],
                 int stride_levels,
                 int proc,
                 armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_datatype_t a1_type;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    switch (datatype)
    {
        case ARMCI_INT:
        case ARMCI_LONG:
        case ARMCI_ACC_INT:
            a1_type = A1_INT32;
            break;
        case ARMCI_LONG_LONG:
            a1_type = A1_INT64;
            break;
        case ARMCI_FLOAT:
        case ARMCI_ACC_FLT:
            a1_type = A1_FLOAT;
            break;
        case ARMCI_DOUBLE:
        case ARMCI_ACC_DBL:
            a1_type = A1_DOUBLE;
            break;
        case ARMCI_ACC_CPL:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_CPL datatype not supported\n");
        case ARMCI_ACC_DCP:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_DCP datatype not supported\n");
        default:
            A1U_ERR_ABORT(status != A1_ERROR, "invalid datatype\n");
    }

    AAP_ARGS("iam %d: A1_NbPutAccS proc = %d, levels = %d, count[0] = %d, count[1] = %d\n",__a1_prof_me,proc,stride_levels,count[0],count[stride_levels-1]);
    AAP_START("A1_NbPutAccS             ");
    status = A1_NbPutAccS(proc,
                          stride_levels,
                          count,
                          src_ptr,
                          src_stride_ar,
                          dst_ptr,
                          dst_stride_ar,
                          a1_type,
                          scale,
                          a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_NbPutAccS returned an error\n");
    AAP_STOP();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_AccV(int datatype,
               void *scale,
               armci_giov_t *dsrc_arr,
               int arr_len,
               int proc)
{
    int status = A1_SUCCESS;
    A1_iov_t *a1_iov_ar;
    A1_reduce_op_t a1_op; /* only used by PutModV */
    A1_datatype_t a1_type;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    switch (datatype)
    {
        case ARMCI_INT:
        case ARMCI_LONG:
        case ARMCI_ACC_INT:
            a1_type = A1_INT32;
            break;
        case ARMCI_LONG_LONG:
            a1_type = A1_INT64;
            break;
        case ARMCI_FLOAT:
        case ARMCI_ACC_FLT:
            a1_type = A1_FLOAT;
            break;
        case ARMCI_DOUBLE:
        case ARMCI_ACC_DBL:
            a1_type = A1_DOUBLE;
            break;
        case ARMCI_ACC_RA:
            a1_type = A1_INT32;
            a1_op = A1_BXOR;
            status = posix_memalign((void **) &a1_iov_ar, 16, sizeof(A1_iov_t) * arr_len);
            A1U_ERR_POP(status != 0, "posix_memalign returned an error\n");
            memcpy((void *) a1_iov_ar, (void *) dsrc_arr, sizeof(A1_iov_t) * arr_len);
            status = A1_PutModV(proc, a1_iov_ar, arr_len, a1_op, a1_type);
            A1U_ERR_POP(status != A1_SUCCESS, "A1_PutModV returned an error\n");
            free(a1_iov_ar);
            goto fn_exit;
            break;
        case ARMCI_ACC_CPL:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_CPL datatype not supported\n");
        case ARMCI_ACC_DCP:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_DCP datatype not supported\n");
        default:
            A1U_ERR_ABORT(status != A1_ERROR, "invalid datatype\n");
    }

    /* ARMCI iov and A1 iov are similar structures but follow
     * different naming conventions. So we make a copy.*/
    status = posix_memalign((void **) &a1_iov_ar, 16, sizeof(A1_iov_t)
            * arr_len);
    A1U_ERR_POP(status != 0, "posix_memalign returned an error\n");

    memcpy((void *) a1_iov_ar, (void *) dsrc_arr, sizeof(A1_iov_t) * arr_len);

    /* accumulate flushes puts before and holds the lock throughout
    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }
    */

    status = A1_PutAccV(proc, a1_iov_ar, arr_len, a1_type, scale);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_PutAccV returned an error\n");

    free(a1_iov_ar);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_NbAccV(int datatype,
                 void *scale,
                 armci_giov_t *dsrc_arr,
                 int arr_len,
                 int proc,
                 armci_hdl_t* handle)
{
    int status = A1_SUCCESS;
    A1_iov_t *a1_iov_ar;
    A1_datatype_t a1_type;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    switch (datatype)
    {
        case ARMCI_INT:
        case ARMCI_LONG:
        case ARMCI_ACC_INT:
            a1_type = A1_INT32;
            break;
        case ARMCI_LONG_LONG:
            a1_type = A1_INT64;
            break;
        case ARMCI_FLOAT:
        case ARMCI_ACC_FLT:
            a1_type = A1_FLOAT;
            break;
        case ARMCI_DOUBLE:
        case ARMCI_ACC_DBL:
            a1_type = A1_DOUBLE;
            break;
        case ARMCI_ACC_CPL:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_CPL datatype not supported\n");
        case ARMCI_ACC_DCP:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_DCP datatype not supported\n");
        default:
            A1U_ERR_ABORT(status != A1_ERROR, "invalid datatype\n");
    }

    /* ARMCI iov and A1 iov are similar structures but follow
     * different naming conventions. So we make a copy.*/
    status = posix_memalign((void **) &a1_iov_ar, 16, sizeof(A1_iov_t)
            * arr_len);
    A1U_ERR_POP(status != 0, "posix_memalign returned an error\n");

    memcpy((void *) a1_iov_ar, (void *) dsrc_arr, sizeof(A1_iov_t) * arr_len);

    status = A1_NbPutAccV(proc, a1_iov_ar, arr_len, a1_type, scale, a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_NbPutAccV returned an error\n");

    free(a1_iov_ar);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_Rmw(int op, void *ploc, void *prem, int value, int proc)
{

    int status = A1_SUCCESS;
    A1_atomic_op_t a1_op;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    if (op == ARMCI_FETCH_AND_ADD || op == ARMCI_FETCH_AND_ADD_LONG)
    {
        a1_op = A1_FETCH_AND_ADD;
    }
    else if (op == ARMCI_SWAP || op == ARMCI_SWAP_LONG)
    {
        a1_op = ARMCI_SWAP;
    }
    else
    {
        A1U_ERR_POP(status != A1_ERROR, "Unsupported rmw operation : %d \n", op);
    }

    /* accumulate flushes puts before and holds the lock throughout
    if ( 1==a1u_settings.armci_strict_ordering )
    {
        status = A1_Flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1_Flush returned an error\n");
    }
    */

    /*Assuming int and long as 32bit signed integers*/
    status = A1_Rmw(proc, &value, ploc, prem, sizeof(int), a1_op, A1_INT32);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Rmw returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_Wait(armci_hdl_t* handle)
{

    int status = A1_SUCCESS;
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    status = A1_Wait_handle(a1_handle);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Wait_handle returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_Test(armci_hdl_t* handle)
{

    int status = A1_SUCCESS;
    A1_handle_t a1_handle;
    A1_bool_t complete;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    a1_handle = (A1_handle_t) *handle;

    status = A1_Test_handle(a1_handle, &complete);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Test_handle returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return !complete;

    fn_fail: goto fn_exit;
}

int ARMCI_WaitAll(void)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1_Wait_handle_all();
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Wait_handle_all returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void ARMCI_Fence(int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1_Flush(proc);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Flush returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void ARMCI_AllFence(void)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1_Flush_group(A1_GROUP_WORLD);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Flush_group returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

int ARMCI_Barrier(void)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1_Sync_group(A1_GROUP_WORLD);
    A1U_ERR_POP(status != A1_SUCCESS, "A1_Sync_group returned an error\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int armci_msg_nproc(void)
{
    int nproc;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    nproc = A1_Process_total(A1_GROUP_WORLD);

    fn_exit: A1U_FUNC_EXIT();
    return nproc;

    fn_fail: goto fn_exit;
}

int armci_msg_me(void)
{
    int me;

    A1U_FUNC_ENTER();

    me = A1_Process_id(A1_GROUP_WORLD);

    fn_exit: A1U_FUNC_EXIT();
    return me;

    fn_fail: goto fn_exit;
}

int armci_domain_id(armci_domain_t domain, int glob_proc_id)
{
    int domain_id;

    A1U_FUNC_ENTER();

    domain_id = 0;

    fn_exit: A1U_FUNC_EXIT();
    return domain_id;

    fn_fail: goto fn_exit;
}

/**********************************************
 Some Dummy ARMCI Function implemenations
 but have to be implemented soon
 ***********************************************/

int ARMCI_Create_mutexes(int num)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void ARMCI_Lock(int mutex, int proc)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void ARMCI_Unlock(int mutex, int proc)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

int ARMCI_Destroy_mutexes(void)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void ARMCI_Copy(void *src, void *dst, int n)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

/*************************************************
 Some More Dummy ARMCI Function implementations
 which are less important
 **************************************************/

int ARMCI_Absolute_id(ARMCI_Group *group, int group_rank)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void ARMCI_Cleanup(void)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void ARMCI_Error(char *message, int code)
{
    A1U_FUNC_ENTER();

    A1_Abort(code, message);

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

int ARMCI_Uses_shm(void)
{
    A1U_FUNC_ENTER();

    fn_exit: A1U_FUNC_EXIT();
    return 0;

    fn_fail: goto fn_exit;
}

void ARMCI_Set_shm_limit(unsigned long shmemlimit)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

int ARMCI_Uses_shm_grp(ARMCI_Group *group)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void ARMCI_Group_get_world(ARMCI_Group *group_out)
{
    A1U_FUNC_ENTER();

    *group_out = A1_GROUP_WORLD;

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void ARMCI_Group_set_default(ARMCI_Group *group)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void ARMCI_Group_create(int n, int *pid_list, ARMCI_Group *group_out)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1 \n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void ARMCI_Group_free(ARMCI_Group *group)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

int ARMCI_Free_group(void *ptr, ARMCI_Group *group)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_Malloc_group(void *ptr_arr[], armci_size_t bytes, ARMCI_Group *group)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void armci_msg_gop_scope(int scope, void *x, int n, char* op, int type)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;
    A1_datatype_t a1_type;

    A1U_FUNC_ENTER();

    if (scope != SCOPE_ALL)
    {
        A1U_ERR_ABORT(A1_ERROR, "Only SCOPE_ALL is supported in armci_msg_gop_scope");
    }

    switch (type)
    {
        case ARMCI_INT:
        case ARMCI_LONG:
            a1_type = A1_INT32;
            break;
        case ARMCI_LONG_LONG:
            a1_type = A1_INT64;
            break;
        case ARMCI_FLOAT:
            a1_type = A1_FLOAT;
            break;
        case ARMCI_DOUBLE:
            a1_type = A1_DOUBLE;
            break;
        default:
            A1U_ERR_ABORT(A1_ERROR, "Invalid datatype received in armci_msg_group_gop_scope");
            break;
    }

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_ABORT(A1_ERROR, "Invalid op received in armci_msg_group_gop_scope");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                a1_type,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_snd(int tag, void* buffer, int len, int to)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_msg_rcv(int tag, void* buffer, int buflen, int *msglen, int from)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

int armci_msg_rcvany(int tag, void* buffer, int buflen, int *msglen)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void armci_write_strided(void *ptr,
                         int stride_levels,
                         int stride_arr[],
                         int count[],
                         char *buf)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_read_strided(void *ptr,
                        int stride_levels,
                        int stride_arr[],
                        int count[],
                        char *buf)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

int ARMCI_PutS_flag_dir(void *src_ptr,
                        int src_stride_arr[],
                        void* dst_ptr,
                        int dst_stride_arr[],
                        int count[],
                        int stride_levels,
                        int *flag,
                        int val,
                        int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_PutS_flag(void *src_ptr,
                    int src_stride_arr[],
                    void* dst_ptr,
                    int dst_stride_arr[],
                    int count[],
                    int stride_levels,
                    int *flag,
                    int val,
                    int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int ARMCI_Same_node(int proc)
{
    int val;

    A1U_FUNC_ENTER();

    /* each process is its own node */
    val = (proc == A1_Process_id(A1_GROUP_WORLD) ? 1 : 0);

    fn_exit: A1U_FUNC_EXIT();
    return val;

    fn_fail: goto fn_exit;
}

int armci_domain_same_id(armci_domain_t domain, int proc)
{
    int val;

    A1U_FUNC_ENTER();

    /* this function always returns false */
    val = 0;

    fn_exit: A1U_FUNC_EXIT();
    return val;

    fn_fail: goto fn_exit;
}

int armci_domain_my_id(armci_domain_t domain)
{
    int val;

    A1U_FUNC_ENTER();

    /* this function always returns false */
    val = 0;

    fn_exit: A1U_FUNC_EXIT();
    return val;

    fn_fail: goto fn_exit;
}

int armci_domain_count(armci_domain_t domain)
{
    int val;

    A1U_FUNC_ENTER();

    /* this function always returns one */
    val = 1;

    fn_exit: A1U_FUNC_EXIT();
    return val;

    fn_fail: goto fn_exit;
}

int armci_domain_nprocs(armci_domain_t domain, int id)
{
    int nproc;

    A1U_FUNC_ENTER();

    nproc = A1_Process_total(A1_GROUP_WORLD);

    fn_exit: A1U_FUNC_EXIT();
    return nproc;

    fn_fail: goto fn_exit;
}

int armci_domain_glob_proc_id(armci_domain_t domain, int id, int loc_proc_id)
{
    int val;

    A1U_FUNC_ENTER();

    /* this function always returns zero */
    val = 0;

    fn_exit: A1U_FUNC_EXIT();
    return val;

    fn_fail: goto fn_exit;
}

void armci_msg_llgop(long long *x, int n, char* op)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_INT64,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_bcast(void* buffer, int len, int root)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    status = A1_Bcast_group(A1_GROUP_WORLD, root, len, buffer);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Bcast_group returned error\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_msg_barrier(void)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    status = A1_Barrier_group(A1_GROUP_WORLD);

    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Barrier_group returned error\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_msg_dgop(double *x, int n, char* op)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_DOUBLE,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_fgop(float *x, int n, char* op)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_FLOAT,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_igop(int *x, int n, char* op)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_INT32,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_lgop(long *x, int n, char* op)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_INT32,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_reduce(void *x, int n, char* op, int datatype)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;
    A1_datatype_t a1_type;

    A1U_FUNC_ENTER();

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    switch (datatype)
    {
        case ARMCI_INT:
        case ARMCI_LONG:
        case ARMCI_ACC_INT:
            a1_type = A1_INT32;
            break;
        case ARMCI_LONG_LONG:
            a1_type = A1_INT64;
            break;
        case ARMCI_FLOAT:
        case ARMCI_ACC_FLT:
            a1_type = A1_FLOAT;
            break;
        case ARMCI_DOUBLE:
        case ARMCI_ACC_DBL:
            a1_type = A1_DOUBLE;
            break;
        case ARMCI_ACC_CPL:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_CPL datatype not supported\n");
        case ARMCI_ACC_DCP:
            A1U_ERR_ABORT(status != A1_ERROR, "ARMCI_ACC_DCP datatype not supported\n");
        default:
            A1U_ERR_ABORT(status != A1_ERROR, "invalid datatype\n");
    }

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                a1_type,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_bintree(int scope, int* Root, int *Up, int *Left, int *Right)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_exchange_address(void *ptr_ar[], int n)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_msg_group_igop(int *x, int n, char* op, ARMCI_Group *group)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    /*We need a check if it is a world group, we assume it for now*/

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_INT32,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_group_lgop(long *x, int n, char* op, ARMCI_Group *group)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    /*We need a check if it is a world group, we assume it for now*/

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_INT32,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_group_llgop(long long *x, int n, char* op, ARMCI_Group *group)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    /*We need a check if it is a world group, we assume it for now*/

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_INT64,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_group_fgop(float *x, int n, char* op, ARMCI_Group *group)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    /*We need a check if it is a world group, we assume it for now*/

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_FLOAT,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_group_dgop(double *x, int n, char* op, ARMCI_Group *group)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;

    A1U_FUNC_ENTER();

    /*We need a check if it is a world group, we assume it for now*/

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else  A1U_ERR_POP(A1_ERROR, "Invalid op received\n");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                A1_DOUBLE,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_msg_group_bcast_scope(int scope,
                                 void *buf,
                                 int len,
                                 int root,
                                 ARMCI_Group *group)
{
    int status = A1_SUCCESS;
    A1U_FUNC_ENTER();

    if (scope != SCOPE_ALL)
    {
        A1U_ERR_ABORT(A1_ERROR, "Only SCOPE_ALL is supported \n");
    }

    /*We need a check if it is a world group, we assume it for now*/

    status = A1_Bcast_group(A1_GROUP_WORLD, root, len, buf);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Bcast_group returned error\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_msg_group_barrier(ARMCI_Group *group)
{
    A1U_FUNC_ENTER();

    /*We need a check if it is a world group, we assume it for now*/

    A1_Barrier_group(A1_GROUP_WORLD);

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_msg_group_gop_scope(int scope,
                               void *x,
                               int n,
                               char* op,
                               int type,
                               ARMCI_Group *group)
{
    int status = A1_SUCCESS;
    A1_reduce_op_t a1_op;
    A1_datatype_t a1_type;

    A1U_FUNC_ENTER();

    if (scope != SCOPE_ALL)
    {
        A1U_ERR_ABORT(A1_ERROR, "Only SCOPE_ALL is supported in armci_msg_gop_scope");
    }

    /*We need a check if it is a world group, we assume it for now*/

    switch (type)
    {
        case ARMCI_INT:
        case ARMCI_LONG:
            a1_type = A1_INT32;
            break;
        case ARMCI_LONG_LONG:
            a1_type = A1_INT64;
            break;
        case ARMCI_FLOAT:
            a1_type = A1_FLOAT;
            break;
        case ARMCI_DOUBLE:
            a1_type = A1_DOUBLE;
            break;
        default:
            A1U_ERR_ABORT(A1_ERROR, "Invalid datatype received in armci_msg_group_gop_scope");
            break;
    }

    if (strncmp(op, "+", 1) == 0) a1_op = A1_SUM;
    else if (strncmp(op, "*", 1) == 0) a1_op = A1_PROD;
    else if (strncmp(op, "max", 3) == 0) a1_op = A1_MAX;
    else if (strncmp(op, "min", 3) == 0) a1_op = A1_MIN;
    else if (strncmp(op, "absmax", 6) == 0) a1_op = A1_MAXABS;
    else if (strncmp(op, "absmin", 6) == 0) a1_op = A1_MINABS;
    else if (strncmp(op, "or", 2) == 0) a1_op = A1_OR;
    else A1U_ERR_ABORT(A1_ERROR, "Invalid op received in armci_msg_group_gop_scope");

    status = A1_Allreduce_group(A1_GROUP_WORLD,
                                n,
                                a1_op,
                                a1_type,
                                x,
                                x);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Allreduce_group returned error\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return;

    fn_fail:
    goto fn_exit;
}

void armci_grp_clus_brdcst(void *buf,
                           int len,
                           int grp_master,
                           int grp_clus_nproc,
                           ARMCI_Group *mastergroup)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_exchange_address_grp(void *ptr_arr[], int n, ARMCI_Group *group)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_msg_bcast_scope(int scope, void* buffer, int len, int root)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_msg_brdcst(void* buffer, int len, int root)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    status = A1_Bcast_group(A1_GROUP_WORLD, root, len, buffer);
    A1U_ERR_ABORT(status != A1_SUCCESS, "A1_Bcast_group returned error\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void armci_msg_sel_scope(int scope,
                         void *x,
                         int n,
                         char* op,
                         int type,
                         int contribute)
{
    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(A1_ERROR, "This function is not supported in ARMCI-A1\n");

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}
