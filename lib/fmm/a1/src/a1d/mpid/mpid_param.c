/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpi2rmadimpl.h"

int A1DI_Read_parameters() {

    int result = A1_SUCCESS;
    char* value = NULL;

    A1U_FUNC_ENTER();

//    if ((value = getenv("A1_ALIGNMENT")) != NULL) {
//        a1_alignment = atoi(value);
//    }

  fn_exit:
    A1U_FUNC_EXIT();
    return result;

  fn_fail:
    goto fn_exit;
}
