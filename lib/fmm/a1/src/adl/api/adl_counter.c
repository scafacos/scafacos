/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

int A1_Create_counter(A1_group_t* group,
                      A1_counter_t *counter)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1D_Create_counter(group,
                                counter);
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Create_counter returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_Destroy_counter(A1_group_t* group,
                       A1_counter_t *counter)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1D_Destroy_counter(group,
                                 counter);
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Destroy_counter returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_Incr_counter(A1_counter_t counter,
                    long increment,
                    long* original)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1D_Incr_counter(counter,
                              increment,
                              original);
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Incr_counter returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
