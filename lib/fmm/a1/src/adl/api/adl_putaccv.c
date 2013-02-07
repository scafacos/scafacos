/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

/* This is here because the build system does not yet have the necessary
 * logic to set these options for each device. */
#define A1D_IMPLEMENTS_PUTACCV

#if defined A1D_IMPLEMENTS_PUTACCV

int A1_PutAccV(int target,
               A1_iov_t *iov_ar,
               int ar_len,
               A1_datatype_t a1_type,
               void* scaling)
{
    int status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
        int i, total_bytes = 0;
        for (i = 0; i < ar_len; i++)
            total_bytes += iov_ar[i].ptr_array_len * iov_ar[i].bytes;
        TAU_TRACE_SENDMSG (A1_TAU_TAG_PUTACCV, target, total_bytes);
    }
#   endif

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_AccV_memcpy(iov_ar, ar_len, a1_type, scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_AccV_memcpy returned an error\n");
    }
    else
    {
        status = A1D_PutAccV(target, iov_ar, ar_len, a1_type, scaling);
        A1U_ERR_POP(status, "A1D_PutAccV returned error\n");
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1_NbPutAccV(int target,
                 A1_iov_t *iov_ar,
                 int ar_len,
                 A1_datatype_t a1_type,
                 void* scaling,
                 A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
        int i, total_bytes = 0;
        for (i = 0; i < ar_len; i++)
            total_bytes += iov_ar[i].ptr_array_len * iov_ar[i].bytes;
        TAU_TRACE_SENDMSG (A1_TAU_TAG_NBPUTACCV, target, total_bytes);
    }
#   endif

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_AccV_memcpy(iov_ar, ar_len, a1_type, scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_AccV_memcpy returned an error\n");
    }
    else
    {
        status = A1D_NbPutAccV(target,
                               iov_ar,
                               ar_len,
                               a1_type,
                               scaling,
                               a1_handle);
        A1U_ERR_POP(status, "A1D_NbPutAccV returned error\n");
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

#else

int A1_PutAccV(int target,
        A1_iov_t *iov_ar,
        int ar_len,
        A1_datatype_t a1_type,
        void* scaling);
{
    int i, j, status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
        int i, total_bytes = 0;
        for (i = 0; i < ar_len; i++)
            total_bytes += iov_ar[i].ptr_array_len * iov_ar[i].bytes;
        TAU_TRACE_SENDMSG (A1_TAU_TAG_PUTACCV, target, total_bytes);
    }
#   endif

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_AccV_memcpy(iov_ar, ar_len, a1_type, scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_AccV_memcpy returned an error\n");
    }
    else
    {
        status = A1D_Allocate_handle(&a1_handle);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Allocate_handle returned error\n");

        for (i=0; i<ar_len; i++)
        {
            for(j=0; j<iov_ar[i].ptr_ar_len; j++)
            {
                status = A1D_NbPutAcc(target,
                                      iov_ar[i].source_ptr_ar[j],
                                      iov_ar[i].target_ptr_ar[j],
                                      iov_ar[i].size,
                                      a1_type,
                                      scaling,
                                      a1_handle);
                A1U_ERR_POP(status != A1_SUCCESS, "A1D_NbPutAcc returned with an error \n");
            }
        }
    }

    status = A1D_Wait_handle(a1_handle);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Wait_handle returned error\n");

    fn_exit:
    /* Could also test for NULL, assuming we set it as such in the declaration. */
    if(target == my_rank && a1u_settings.network_bypass) A1D_Release_handle(a1_handle);
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

int A1_NbPutAccV(int target,
        A1_iov_t *iov_ar,
        int ar_len,
        A1_datatype_t a1_type,
        void* scaling,
        A1_handle_t a1_handle);
{
    int i, j, status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
        int i, total_bytes = 0;
        for (i = 0; i < ar_len; i++)
            total_bytes += iov_ar[i].ptr_array_len * iov_ar[i].bytes;
        TAU_TRACE_SENDMSG (A1_TAU_TAG_NBPUTACCV, target, total_bytes);
    }
#   endif

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_AccV_memcpy(iov_ar, ar_len, a1_type, scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_AccV_memcpy returned an error\n");
    }
    else
    {
        for (i=0; i<ar_len; i++)
        {
            for(j=0; j<iov_ar[i].ptr_ar_len; j++)
            {
                status = A1D_NbPutAcc(target,
                                      iov_ar[i].source_ptr_ar[j],
                                      iov_ar[i].target_ptr_ar[j],
                                      iov_ar[i].size,
                                      a1_type,
                                      scaling,
                                      a1_handle);
                A1U_ERR_POP(status != A1_SUCCESS, "A1D_NbPutAcc returned with an error \n");
            }
        }
    }

    fn_exit:
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

#endif /* A1D_IMPLEMENTS_PUTACCV */
