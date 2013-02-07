/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

int A1_Rmw(int target,
           void* source_ptr_in,
           void* source_ptr_out,
           void* target_ptr,
           int bytes,
           A1_atomic_op_t op,
           A1_datatype_t a1_type)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
      TAU_TRACE_SENDMSG (A1_TAU_TAG_RMW, target, bytes);
    }
#   endif

    status = A1D_Rmw(target,
                     source_ptr_in,
                     source_ptr_out,
                     target_ptr,
                     bytes,
                     op,
                     a1_type);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Rmw returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
