/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

int A1_Exchange_segments(A1_group_t* group, void* ptr[])
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1D_Exchange_segments(group, ptr);
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Exchange_segments returned an error\n");

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail:  
    goto fn_exit;
}

int  A1_Alloc_segment(void** pointer, int bytes)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

    if(bytes == 0)
    {
      *pointer = NULL;
      goto fn_exit;
    }

    status = A1D_Alloc_segment(pointer, bytes);
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Alloc_segment returned an error\n");

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

