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
#define A1D_IMPLEMENTS_PUTV

#if defined A1D_IMPLEMENTS_PUTV

int A1_PutV(int target, A1_iov_t *iov_ar, int ar_len)
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
        TAU_TRACE_SENDMSG (A1_TAU_TAG_PUTV, target, total_bytes);
    }
#   endif

    /* It isn't worth trying to optimize for the threshold here because
     * these operations aren't used much in GA. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_PutV_memcpy(iov_ar, ar_len);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1U_PutV_memcpy returned error\n");
    }
    else
    {
        status = A1D_PutV(target, iov_ar, ar_len);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1D_PutV returned error\n");
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1_NbPutV(int target, A1_iov_t *iov_ar, int ar_len, A1_handle_t a1_handle)
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
        TAU_TRACE_SENDMSG (A1_TAU_TAG_NBPUTV, target, total_bytes);
    }
#   endif

    /* It isn't worth trying to optimize for the threshold here because
     * these operations aren't used much in GA. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_PutV_memcpy(iov_ar, ar_len);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1U_PutV_memcpy returned error\n");
    }
    else
    {
        status = A1D_NbPutV(target, iov_ar, ar_len, a1_handle);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1D_NbPutV returned error\n");
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

#else

int A1_PutV(int target,
        A1_iov_t *iov_ar,
        int ar_len)
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
        TAU_TRACE_SENDMSG (A1_TAU_TAG_PUTV, target, total_bytes);
    }
#   endif

    /* It isn't worth trying to optimize for the threshold here because
     * these operations aren't used much in GA. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_PutV_memcpy(iov_ar, ar_len);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1U_PutV_memcpy returned error\n");
    }
    else
    {
        status = A1D_Allocate_handle(&a1_handle);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Allocate_handle returned error\n");

        for (i=0; i<ar_len; i++)
        {
            for(j=0; j<iov_ar[i].ptr_ar_len; j++)
            {
                status = A1D_NbPut(target,
                        iov_ar[i].source_ptr_ar[j],
                        iov_ar[i].target_ptr_ar[j],
                        iov_ar[i].size,
                        a1_handle);
                A1U_ERR_POP(status != A1_SUCCESS, "A1D_NbPut returned with an error \n");
            }
        }

        status = A1D_Wait_handle(a1_handle);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Wait_handle returned error\n");
    }

    fn_exit:
    /* Could also test for NULL, assuming we set it as such in the declaration. */
    if(target == my_rank && a1u_settings.network_bypass) A1D_Release_handle(a1_handle);
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

int A1_NbPutV(int target,
        A1_iov_t *iov_ar,
        int ar_len,
        A1_handle_t a1_handle)
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
        TAU_TRACE_SENDMSG (A1_TAU_TAG_NBPUTV, target, total_bytes);
    }
#   endif

    /* It isn't worth trying to optimize for the threshold here because
     * these operations aren't used much in GA. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_PutV_memcpy(iov_ar, ar_len);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1U_PutV_memcpy returned error\n");
    }
    else
    {
        for (i=0; i<ar_len; i++)
        {
            for(j=0; j<iov_ar[i].ptr_ar_len; j++)
            {
                status = A1D_NbPut(target,
                        iov_ar[i].source_ptr_ar[j],
                        iov_ar[i].target_ptr_ar[j],
                        iov_ar[i].size,
                        a1_handle);
                A1U_ERR_POP(status != A1_SUCCESS, "A1D_NbPut returned with an error \n");
            }
        }
    }

    fn_exit:
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

#endif /* A1D_IMPLEMENTS_PUTV */
