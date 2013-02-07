/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Protocol_t A1D_Generic_put_protocol;

int A1DI_Put_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Put_Configuration_t conf;

    A1U_FUNC_ENTER();

    conf.protocol = DCMF_DEFAULT_PUT_PROTOCOL;
    conf.network = DCMF_TORUS_NETWORK;
    status = DCMF_Put_register(&A1D_Generic_put_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "DCMF_Put_register returned with error %d \n",
                status);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_Put(int target,
            void* src,
            void* dst,
            int bytes)
{
    int status = A1_SUCCESS;
    DCMF_Request_t request;
    DCMF_Callback_t done_callback, ack_callback;
    volatile int done_active, ack_active;
    size_t src_disp, dst_disp;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    done_callback.function = A1DI_Generic_done;
    done_callback.clientdata = (void *) &done_active;
    done_active = 1;
    A1D_Connection_put_active[target]++;

    src_disp = (size_t) src - (size_t) A1D_Membase_global[A1D_Process_info.my_rank];
    dst_disp = (size_t) dst - (size_t) A1D_Membase_global[target];

    status = DCMF_Put(&A1D_Generic_put_protocol,
                      &request,
                      done_callback,
                      DCMF_SEQUENTIAL_CONSISTENCY,
                      target,
                      bytes,
                      &A1D_Memregion_global[A1D_Process_info.my_rank],
                      &A1D_Memregion_global[target],
                      src_disp,
                      dst_disp,
                      A1D_Nocallback);
    A1U_ERR_POP(status, "DCMF_Put returned with an error \n");

    A1DI_Conditional_advance(done_active > 0);

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
   goto fn_exit;
}

int A1D_NbPut(int target,
              void* src,
              void* dst,
              int bytes,
              A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle;
    A1D_Request_t *a1d_request;
    DCMF_Callback_t done_callback, ack_callback;
    size_t src_disp, dst_disp;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_handle = (A1D_Handle_t *) a1_handle;

    a1d_request = A1DI_Get_request(1);
    A1U_ERR_POP(status = (a1d_request == NULL),
                "A1DI_Get_request returned error \n");
    A1DI_Set_handle(a1d_request, a1d_handle);

    done_callback.function = A1DI_Request_done;
    done_callback.clientdata = (void *) a1d_request;

    a1d_handle->active++;

    src_disp = (size_t) src - (size_t) A1D_Membase_global[A1D_Process_info.my_rank];
    dst_disp = (size_t) dst - (size_t) A1D_Membase_global[target];

    status = DCMF_Put(&A1D_Generic_put_protocol,
                      &(a1d_request->request),
                      done_callback,
                      DCMF_SEQUENTIAL_CONSISTENCY,
                      target,
                      bytes,
                      &A1D_Memregion_global[A1D_Process_info.my_rank],
                      &A1D_Memregion_global[target],
                      src_disp,
                      dst_disp,
                      A1D_Nocallback);
    A1U_ERR_POP(status != DCMF_SUCCESS, "DCMF_Put returned with an error \n");

    A1D_Connection_put_active[target]++;

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}
