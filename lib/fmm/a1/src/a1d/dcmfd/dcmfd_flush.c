/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Protocol_t A1D_Send_flush_protocol;
DCMF_Protocol_t A1D_Control_flushack_protocol;
void **A1D_Put_Flushcounter_ptr;
volatile int A1D_Control_flushack_active;
volatile int A1D_Put_flushack_active;
volatile int *A1D_Connection_send_active;
volatile int *A1D_Connection_put_active;

/**************************************************************** 
 * Control protocol used to send acknowledgements to send-flush *
 * messages                                                     *
 ****************************************************************/

/* TODO: If this is just decrementing the payload like A1DI_Generic_done, why don't we have
 * another A1D_Generic_control_done?  Better yet, call them A1DI_Send_decrement and
 * A1D_Control_decrement. */
void A1DI_Control_flushack_callback(void *clientdata,
                                    const DCMF_Control_t *info,
                                    size_t peer)
{
    --(*((uint32_t *) clientdata));

}

int A1DI_Control_flushack_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Control_Configuration_t conf;

    A1U_FUNC_ENTER();

    A1D_Control_flushack_active = 0;

    conf.protocol = DCMF_DEFAULT_CONTROL_PROTOCOL;
    conf.network = DCMF_DEFAULT_NETWORK;
    conf.cb_recv = A1DI_Control_flushack_callback;
    conf.cb_recv_clientdata = (void *) &A1D_Control_flushack_active;

    status = DCMF_Control_register(&A1D_Control_flushack_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "Control flushack registartion returned with error %d \n",
                status);

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

/****************************************************************
 * Send protocol used to flush send message                     *
 ****************************************************************/

void A1DI_RecvSendShort_flush_callback(void *clientdata,
                                       const DCQuad *msginfo,
                                       unsigned count,
                                       size_t peer,
                                       const char *src,
                                       size_t bytes)
{
    int status = A1_SUCCESS;
    DCMF_Control_t info;

    status = DCMF_Control(&A1D_Control_flushack_protocol,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          peer,
                          &info);
    A1U_ERR_ABORT(status != DCMF_SUCCESS,
                  "DCMF_Control failed in A1DI_RecvSendShort_flush_callback\n");

}

int A1DI_Send_flush_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Send_Configuration_t conf;

    A1U_FUNC_ENTER();

    /* FIXME: The recv callback should be implemented when Send might be used *
     * with large messages */

    conf.protocol = DCMF_DEFAULT_SEND_PROTOCOL;
    conf.network = DCMF_TORUS_NETWORK;
    conf.cb_recv_short = A1DI_RecvSendShort_flush_callback;
    conf.cb_recv_short_clientdata = NULL;
    conf.cb_recv = NULL;
    conf.cb_recv_clientdata = NULL;

    status = DCMF_Send_register(&A1D_Send_flush_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "DCMF_Send_register failed with error %d \n",
                status);

    /* Allocating memory for vector that tracks connections with active sends */
    status = A1DI_Malloc((void **) &A1D_Connection_send_active,
                                 sizeof(int) * A1D_Process_info.num_ranks);
    A1U_ERR_POP(status != 0,
                "A1DI_Malloc failed \n");
    A1DI_Memset((void *) A1D_Connection_send_active,
                0,
                sizeof(int) * A1D_Process_info.num_ranks);

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

/****************************************************************
 * Put protocol used to flush put messages                      *
 ****************************************************************/

int A1DI_Put_flush_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Control_t info;
    int rank;

    A1U_FUNC_ENTER();

    status = A1DI_Malloc((void **) &A1D_Put_Flushcounter_ptr,
                                 sizeof(void *) * A1D_Process_info.num_ranks);
    A1U_ERR_POP(status != 0, "A1DI_Malloc failed \n");

    status = A1DI_Malloc((void **) &(A1D_Put_Flushcounter_ptr[A1D_Process_info.my_rank]),
                                  2);
    A1U_ERR_POP(status != 0, "A1DI_Malloc failed \n");

    /*TODO: Use DCMF_Send operations instead to exploit TORUS network */
    A1D_Control_xchange_info.xchange_ptr = (void *) A1D_Put_Flushcounter_ptr;
    A1D_Control_xchange_info.xchange_size = sizeof(void *);
    A1D_Control_xchange_info.rcv_active += A1D_Process_info.num_ranks - 1;

    A1DI_GlobalBarrier();

    A1DI_Memcpy((void *) &info,
                (void *) &(A1D_Put_Flushcounter_ptr[A1D_Process_info.my_rank]),
                sizeof(void *));
    for (rank = 0; rank < A1D_Process_info.num_ranks; rank++)
    {
        likely_if (rank != A1D_Process_info.my_rank)
        {
            status = DCMF_Control(&A1D_Control_xchange_info.protocol,
                                  DCMF_SEQUENTIAL_CONSISTENCY,
                                  rank,
                                  &info);
            A1U_ERR_POP(status != DCMF_SUCCESS, "DCMF_Control failed\n");
        }
    }
    A1DI_Conditional_advance(A1D_Control_xchange_info.rcv_active > 0);

    /* Allocating memory for vector thats tracks connections with active puts */
    status = A1DI_Malloc((void **) &A1D_Connection_put_active,
                                 sizeof(int) * A1D_Process_info.num_ranks);
    A1U_ERR_POP(status != 0,
                "Connection put active buffer allocation Failed \n");
    A1DI_Memset((void *) A1D_Connection_put_active,
                0,
                sizeof(int) * A1D_Process_info.num_ranks);

    A1D_Put_flushack_active = 0;

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;

}

int A1DI_Send_flush(int proc)
{
    int status = A1_SUCCESS;
    DCMF_Request_t request;
    DCQuad msginfo;

    A1U_FUNC_ENTER();

    A1D_Control_flushack_active++;

    status = DCMF_Send(&A1D_Send_flush_protocol,
                       &request,
                       A1D_Nocallback,
                       DCMF_SEQUENTIAL_CONSISTENCY,
                       proc,
                       0,
                       NULL,
                       &msginfo,
                       1);
    A1U_ERR_POP(status != DCMF_SUCCESS, "DCMF_Send returned with an error \n");

    A1DI_Conditional_advance(A1D_Control_flushack_active > 0);

    A1D_Connection_send_active[proc] = 0;
    A1D_Connection_put_active[proc] = 0;

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1DI_Send_flush_start(int proc)
{
    int status = A1_SUCCESS;
    DCMF_Request_t request;
    DCMF_Callback_t callback;
    DCQuad msginfo;

    A1U_FUNC_ENTER();

    A1D_Control_flushack_active++;

    status = DCMF_Send(&A1D_Send_flush_protocol,
                       &request,
                       A1D_Nocallback,
                       DCMF_SEQUENTIAL_CONSISTENCY,
                       proc,
                       0,
                       NULL,
                       &msginfo,
                       1);
    A1U_ERR_POP(status != DCMF_SUCCESS, "DCMF_Send returned with an error \n");

    A1D_Connection_send_active[proc] = 0;
    A1D_Connection_put_active[proc] = 0;

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1DI_Put_flush(int proc)
{
    int status = A1_SUCCESS;
    DCMF_Request_t request;
    DCMF_Callback_t callback;
    volatile int active;
    size_t src_disp, dst_disp;

    A1U_FUNC_ENTER();

    callback.function = A1DI_Generic_done;
    callback.clientdata = (void *) &A1D_Put_flushack_active;

    src_disp = (size_t)(A1D_Put_Flushcounter_ptr[A1D_Process_info.my_rank])
             - (size_t) A1D_Membase_global[A1D_Process_info.my_rank];
    dst_disp = (size_t)(A1D_Put_Flushcounter_ptr[proc])
             - (size_t) A1D_Membase_global[proc] + 1;

    A1D_Connection_put_active[proc] = 0;
    A1D_Put_flushack_active++;

    status = DCMF_Put(&A1D_Generic_put_protocol,
                      &request,
                      A1D_Nocallback,
                      DCMF_SEQUENTIAL_CONSISTENCY,
                      proc,
                      1,
                      &A1D_Memregion_global[A1D_Process_info.my_rank],
                      &A1D_Memregion_global[proc],
                      src_disp,
                      dst_disp,
                      callback);
    A1U_ERR_POP(status != DCMF_SUCCESS, "DCMF_Put returned with an error \n");

    A1DI_Conditional_advance(A1D_Put_flushack_active > 0);

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1DI_Put_flush_start(int proc)
{
    int status = A1_SUCCESS;
    DCMF_Request_t request;
    DCMF_Callback_t ack_callback;
    size_t src_disp, dst_disp;

    A1U_FUNC_ENTER();

    ack_callback.function = A1DI_Generic_done;
    ack_callback.clientdata = (void *) &A1D_Put_flushack_active;

    src_disp = (size_t)(A1D_Put_Flushcounter_ptr[A1D_Process_info.my_rank])
             - (size_t) A1D_Membase_global[A1D_Process_info.my_rank];
    dst_disp = (size_t)(A1D_Put_Flushcounter_ptr[proc])
             - (size_t) A1D_Membase_global[proc] + 1;

    A1D_Put_flushack_active++;

    status = DCMF_Put(&A1D_Generic_put_protocol,
                      &request,
                      A1D_Nocallback,
                      DCMF_SEQUENTIAL_CONSISTENCY,
                      proc,
                      1,
                      &A1D_Memregion_global[A1D_Process_info.my_rank],
                      &A1D_Memregion_global[proc],
                      src_disp,
                      dst_disp,
                      ack_callback);
    A1U_ERR_POP(status != DCMF_SUCCESS, "DCMF_Put returned with an error \n");

    A1D_Connection_put_active[proc] = 0;

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1D_Flush(int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    if (A1D_Connection_send_active[proc])
    {
        status = A1DI_Send_flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1DI_Send_flush failed\n");
    }
    else if (A1D_Connection_put_active[proc])
    {
        status = A1DI_Put_flush(proc);
        A1U_ERR_POP(status != A1_SUCCESS, "A1DI_Put_flush failed \n");
    }

    A1DI_Conditional_advance(A1D_Put_flushack_active > 0 || A1D_Control_flushack_active > 0);

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}
