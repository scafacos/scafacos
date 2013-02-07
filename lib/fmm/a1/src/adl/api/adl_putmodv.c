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
#define A1D_IMPLEMENTS_PUTMODV

#if defined A1D_IMPLEMENTS_PUTMODV

int A1_PutModV(int target,
               A1_iov_t *iov_ar,
               int ar_len,
               A1_reduce_op_t a1_op,
               A1_datatype_t a1_type)
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
        TAU_TRACE_SENDMSG (A1_TAU_TAG_PUTMODV, target, total_bytes);
    }
#   endif

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_ModV_memcpy(iov_ar, ar_len, a1_op, a1_type);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_ModV_memcpy returned an error\n");
    }
    else
    {
        status = A1D_PutModV(target, iov_ar, ar_len, a1_op, a1_type);
        A1U_ERR_POP(status, "A1D_PutModV returned error\n");
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

#else

int A1_PutModV(int target,
        A1_iov_t *iov_ar,
        int ar_len,
        A1_reduce_op_t a1_op,
        A1_datatype_t a1_type);
{
    int i, j, status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    A1U_ERR_POP(status!=A1_SUCCESS, "A1_PutModV requires device-level implementation.\n");

    fn_exit:
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

#endif /* A1D_IMPLEMENTS_PUTMODV */
