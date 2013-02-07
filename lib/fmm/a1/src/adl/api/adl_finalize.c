/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

int A1_Finalize(void)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

    status = A1D_Finalize();
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Finalize returned error\n");

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

void A1_Abort(int error_code, char error_message[])
{
    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

    A1D_Abort(error_code,  error_message);

  fn_exit: 
    A1U_FUNC_EXIT();
    return;

  fn_fail: 
    goto fn_exit;
}
