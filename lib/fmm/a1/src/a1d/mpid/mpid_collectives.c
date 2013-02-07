/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpidimpl.h"

void A1D_Barrier_group(A1_group_t* group) {

    A1U_FUNC_ENTER();

    if(group == A1_GROUP_WORLD) {
        MPI_Barrier(MPI_COMM_WORLD);
        goto fn_exit;
    } else {
        A1U_ERR_ABORT(1,"A1D_Barrier_group not implement for non-world groups!");
        goto fn_fail;
    }

  fn_exit:
    A1U_FUNC_EXIT();
    return;

  fn_fail:
    goto fn_exit;

}
