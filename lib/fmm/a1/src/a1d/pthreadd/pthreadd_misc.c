/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "pthreaddimpl.h"

int A1D_Process_id(A1_group_t* group)
{
    int id;

    A1U_FUNC_ENTER();

    if (group == A1_GROUP_WORLD || group == NULL)
    {
        id = A1D_Process_info.my_rank;
        goto fn_exit;
    }
    else
    {
        id = -1;
        A1U_ERR_ABORT(1, "A1D_Process_id not implement for non-world groups!");
        goto fn_fail;
    }

  fn_exit:
    A1U_FUNC_EXIT();
    return id;

  fn_fail:
    goto fn_exit;
}

int A1D_Process_total(A1_group_t* group)
{
    int total;

    A1U_FUNC_ENTER();

    if (group == A1_GROUP_WORLD || group == NULL)
    {
        total = A1D_Process_info.num_ranks;
        goto fn_exit;
    }
    else
    {
        total = -1;
        A1U_ERR_ABORT(1,
                      "A1D_Process_total not implement for non-world groups!");
        goto fn_fail;
    }

  fn_exit:
    A1U_FUNC_EXIT();
    return total;

  fn_fail:
    goto fn_exit;
}

int A1D_Node_id(A1_group_t* group)
{
    int id;

    A1U_FUNC_ENTER();

    if (group == A1_GROUP_WORLD || group == NULL)
    {
        id = A1D_Process_info.my_node;
        goto fn_exit;
    }
    else
    {
        id = -1;
        A1U_ERR_ABORT(1, "A1D_Node_id not implement for non-world groups!");
        goto fn_fail;
    }

  fn_exit:
    A1U_FUNC_EXIT();
    return id;

  fn_fail:
    goto fn_exit;
}

int A1D_Node_total(A1_group_t* group)
{
    int total;

    A1U_FUNC_ENTER();

    if (group == A1_GROUP_WORLD || group == NULL)
    {
        total = A1D_Process_info.num_nodes;
        goto fn_exit;
    }
    else
    {
        total = -1;
        A1U_ERR_ABORT(1, "A1D_Node_total not implement for non-world groups!");
        goto fn_fail;
    }

  fn_exit:
    A1U_FUNC_EXIT();
    return total;

  fn_fail:
    goto fn_exit;
}

double A1D_Time_seconds()
{
    A1U_FUNC_ENTER();

    /* TODO: implement this function */

  fn_exit:
    A1U_FUNC_EXIT();
    return 0.0;

  fn_fail:
    goto fn_exit;
}

unsigned long long A1D_Time_cycles()
{
    A1U_FUNC_ENTER();

    /* TODO: implement this function */

  fn_exit:
    A1U_FUNC_EXIT();
    return 0;

  fn_fail:
    goto fn_exit;
}
