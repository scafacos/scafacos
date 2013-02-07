/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Protocol_t A1D_Generic_get_protocol;

int A1DI_Get_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Get_Configuration_t conf;

    A1U_FUNC_ENTER();

    conf.protocol = DCMF_DEFAULT_GET_PROTOCOL;
    conf.network = DCMF_TORUS_NETWORK;
    status = DCMF_Get_register(&A1D_Generic_get_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "get registartion returned with error %d \n",
                status);

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1D_Get(int target, void* src, void* dst, int bytes)
{
    int status = A1_SUCCESS;
    DCMF_Request_t request;
    DCMF_Callback_t callback;
    volatile int active;
    unsigned src_disp, dst_disp;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    callback.function = A1DI_Generic_done;
    callback.clientdata = (void *) &active;

    src_disp = (size_t) src - (size_t) A1D_Membase_global[target];
    dst_disp = (size_t) dst - (size_t) A1D_Membase_global[A1D_Process_info.my_rank];

    active = 1;

    status = DCMF_Get(&A1D_Generic_get_protocol,
                      &request,
                      callback,
                      DCMF_RELAXED_CONSISTENCY,
                      target,
                      bytes,
                      &A1D_Memregion_global[target],
                      &A1D_Memregion_global[A1D_Process_info.my_rank],
                      src_disp,
                      dst_disp);
    A1U_ERR_POP(status, "DCMF_Get returned with an error \n");

    A1DI_Conditional_advance(active > 0);

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}


int A1D_NbGet(int target, void* src, void* dst, int bytes, A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    A1D_Handle_t* a1d_handle = NULL;
    A1D_Request_t* a1d_request = NULL;
    DCMF_Callback_t callback;
    unsigned src_disp, dst_disp;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_handle = (A1D_Handle_t *) a1_handle;

    a1d_handle->active++;

    a1d_request = A1DI_Get_request(1);
    A1U_ERR_POP(status = (a1d_request == NULL),
                "A1DI_Get_request returned error.\n");
    A1DI_Set_handle(a1d_request, a1d_handle); 

    callback.function = A1DI_Request_done;
    callback.clientdata = (void *) a1d_request;

    src_disp = (size_t) src - (size_t) A1D_Membase_global[target];
    dst_disp = (size_t) dst - (size_t) A1D_Membase_global[A1D_Process_info.my_rank];

    status = DCMF_Get(&A1D_Generic_get_protocol,
                      &(a1d_request->request),
                      callback,
                      DCMF_RELAXED_CONSISTENCY,
                      target,
                      bytes,
                      &A1D_Memregion_global[target],
                      &A1D_Memregion_global[A1D_Process_info.my_rank],
                      src_disp,
                      dst_disp);
    A1U_ERR_POP(status, "DCMF_Get returned with an error \n");

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}
