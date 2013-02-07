/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

int A1D_Release_segments(A1_group_t* group, void *ptr)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    /*This functions does nothing becuase BG does not involve
      any registration. It has to do a barrier syncrhonization
      to ensure everyone is agreeing on the release*/

    if (group == A1_GROUP_WORLD || group == NULL)
    {
        status = A1DI_GlobalBarrier();
        A1U_ERR_ABORT(status != A1_SUCCESS, "DCMF_GlobalBarrier returned with an error");
        goto fn_exit;
    }
    else
    {
        A1U_ERR_POP(1, "A1D_Release_segments not implemented for non-world groups!");
        goto fn_fail;
    }

  fn_exit:
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1D_Free_segment(void *ptr)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    A1DI_Free(ptr);

  fn_exit:
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
