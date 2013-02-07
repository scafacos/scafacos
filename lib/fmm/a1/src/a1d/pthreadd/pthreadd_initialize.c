/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "pthreaddimpl.h"

pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

int A1D_Initialize(int thread_level)
{

    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

//    status = A1DI_Read_parameters();
//    A1U_ERR_POP(status != A1_SUCCESS,
//                "A1DI_Read_parameters returned with error \n");

    switch (thread_level)
    {
    case A1_THREAD_SINGLE:
        a1d_settings.thread_safe = 0;
        break;
    case A1_THREAD_FUNNELED:
        a1d_settings.thread_safe = 0;
        break;
    case A1_THREAD_SERIALIZED:
        a1d_settings.thread_safe = 0;
        break;
    case A1_THREAD_MULTIPLE:
        a1d_settings.thread_safe = 1;
        break;
    default:
        A1U_ERR_POP(A1_ERROR,
                    "Unsupported thread level provided in A1D_Initialize \n");
        break;
    }

    /* NOTE: Clearly, the Pthread device only supports a single process. */
    A1D_Process_info.my_rank = 0;
    A1D_Process_info.num_ranks = 1;

    A1D_Process_info.my_node = 0;
    A1D_Process_info.num_nodes = 1;

//    status = A1DI_Print_parameters();
//    A1U_ERR_POP(status != A1_SUCCESS,
//                "A1DI_Print_parameters returned with error \n");

    A1DI_CRITICAL_ENTER();

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}
