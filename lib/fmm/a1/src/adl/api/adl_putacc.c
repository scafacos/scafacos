/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

int A1_PutAcc(int target,
              void* source_ptr,
              void* target_ptr,
              int bytes,
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
        TAU_TRACE_SENDMSG (A1_TAU_TAG_PUTACC, target, bytes);
    }
#   endif

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_Acc_memcpy(source_ptr,
                                target_ptr,
                                bytes,
                                a1_type,
                                scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_Acc_memcpy returned an error\n");
    }
    else
    {
        status = A1D_PutAcc(target,
                            source_ptr,
                            target_ptr,
                            bytes,
                            a1_type,
                            scaling);
        A1U_ERR_POP((status!=A1_SUCCESS), "A1D_PutAcc returned error\n");
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1_NbPutAcc(int target,
                void* source_ptr,
                void* target_ptr,
                int bytes,
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
        TAU_TRACE_SENDMSG (A1_TAU_TAG_NBPUTACC, target, bytes);
    }
#   endif

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_Acc_memcpy(source_ptr,
                                target_ptr,
                                bytes,
                                a1_type,
                                scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_Acc_memcpy returned an error\n");
    }
    else
    {
        status = A1D_NbPutAcc(target,
                              source_ptr,
                              target_ptr,
                              bytes,
                              a1_type,
                              scaling,
                              a1_handle);
        A1U_ERR_POP((status!=A1_SUCCESS), "A1D_NbPutAcc returned error\n");
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}
