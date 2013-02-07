/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpidimpl.h"

int A1D_Initialize(int thread_level)
{

    int result = MPI_SUCCESS;
    int required;
    int provided;
    const int zero = 0;

    A1U_FUNC_ENTER();

    switch (thread_level)
    {
    case A1_THREAD_MULTIPLE:
        required = MPI_THREAD_MULTIPLE;
        break;

    case A1_THREAD_SERIALIZED:
        required = MPI_THREAD_SERIALIZED;
        break;

    case A1_THREAD_FUNNELED:
        required = MPI_THREAD_FUNNELED;
        break;

    case A1_THREAD_SINGLE:
        required = MPI_THREAD_SINGLE;
        break;

    default:
        A1U_ERR_POP(1, "Invalid choice for A1_thread_level\n");

    }

    MPI_Init_thread(&zero, NULL, required, &provided);

    A1U_ERR_POP(required != provided,
                "MPI cannot provide requested thread support.");

    A1D_Messager_info.thread_level = thread_level;

    MPI_Comm_size(MPI_COMM_WORLD, &(A1D_Process_info.num_ranks));
    MPI_Comm_rank(MPI_COMM_WORLD, &(A1D_Process_info.my_rank));

    result = A1DI_Read_parameters();
    A1U_ERR_POP(result != A1_SUCCESS, "A1DI_Read_parameters failed");

    /* FIXME: Need to do stuff here! */

    fn_exit: A1U_FUNC_EXIT();
    return result;

    fn_fail: goto fn_exit;
}

