/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

int A1_Put(int target, void* src, void* dst, int bytes)
{
    int status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
        TAU_TRACE_SENDMSG (A1_TAU_TAG_PUT, target, bytes);
    }
#   endif

    if(target == my_rank && (bytes < a1u_settings.network_bypass_upper_limit_1d) )
    {
       status = A1U_Put_memcpy(src, dst, bytes);
       A1U_ERR_POP(status != A1_SUCCESS, "A1U_Put_memcpy returned an error\n");
    }
    else
    {
        status = A1D_Put(target, src, dst, bytes);
        A1U_ERR_POP(status != A1_SUCCESS, "A1D_Put returned an error\n");
    }

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1_NbPut(int target, void* src, void* dst, int bytes, A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
        TAU_TRACE_SENDMSG (A1_TAU_TAG_NBPUT, target, bytes);
    }
#   endif

    /* Not sure if what is the right strategy for bypass.  A1U_*_memcpy are blocking
     * but the overhead of going into DCMF_Put is likely not worth the savings
     * from said call being non-blocking.  This is especially true under heavy load
     * since we have determined that DMA vs. memcpy turns over when the NIC is getting
     * hammered.
     */
    if(target == my_rank && (bytes < a1u_settings.network_bypass_upper_limit_1d) )
    {
       status = A1U_Put_memcpy(src, dst, bytes);
       A1U_ERR_POP(status != A1_SUCCESS, "A1U_Put_memcpy returned an error\n");
    }
    else
    {
        status = A1D_NbPut(target, src, dst, bytes, a1_handle);
        A1U_ERR_POP(status != A1_SUCCESS, "A1D_NbPut returned an error\n");
    }

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}
