/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

int A1_Flush(int proc) 
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
      TAU_TRACE_SENDMSG (A1_TAU_TAG_FLUSH, proc, 8);
    }
#   endif

    status = A1D_Flush(proc); 
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Flush returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
