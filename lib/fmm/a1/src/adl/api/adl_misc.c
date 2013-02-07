/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

int A1_Process_id(A1_group_t* group)
{
    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

  fn_exit: 
    A1U_FUNC_EXIT();
    return A1D_Process_id(group);

  fn_fail: 
    goto fn_exit;
}

int A1_Process_total(A1_group_t* group)
{
    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

  fn_exit: 
    A1U_FUNC_EXIT();
    return A1D_Process_total(group);

  fn_fail: 
    goto fn_exit;
}

int A1_Node_id(A1_group_t* group)
{
    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

  fn_exit: 
    A1U_FUNC_EXIT();
    return A1D_Node_id(group);

  fn_fail: 
    goto fn_exit;
}

int A1_Node_total(A1_group_t* group)
{
    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

  fn_exit: 
    A1U_FUNC_EXIT();
    return A1D_Node_total(group);

  fn_fail: 
    goto fn_exit;
}

double A1_Time_seconds(void)
{
    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

  fn_exit: 
    A1U_FUNC_EXIT();
    return A1D_Time_seconds();

  fn_fail: 
    goto fn_exit;
}

unsigned long long A1_Time_cycles(void)
{
    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

  fn_exit: 
    A1U_FUNC_EXIT();
    return A1D_Time_cycles();

  fn_fail: 
    goto fn_exit;
}
