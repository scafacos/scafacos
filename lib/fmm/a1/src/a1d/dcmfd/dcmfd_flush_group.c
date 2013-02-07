/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

void A1DI_Flush_all() 
{
    int status = A1_SUCCESS;
    int dst,request_count;
    DCMF_Request_t *request;
    DCQuad msginfo;
    DCMF_Callback_t ack_callback;
    volatile int pending_count;
    size_t src_disp, dst_disp;

    A1U_FUNC_ENTER();

    request_count = (a1d_settings.flushall_pending_limit < A1D_Process_info.num_ranks) ? 
                      a1d_settings.flushall_pending_limit : A1D_Process_info.num_ranks;
    status = A1DI_Malloc((void **) &request, sizeof(DCMF_Request_t) * request_count); 
    A1U_ERR_POP(status != 0, "A1DI_Malloc failed in A1DI_Flush_all\n"); 

    pending_count = 0;
    ack_callback.function = A1DI_Generic_done;
    ack_callback.clientdata = (void *) &A1D_Put_flushack_active;

    for (dst = 0; dst < A1D_Process_info.num_ranks; dst++)
    {
        likely_if (dst != A1D_Process_info.my_rank)
        {

            if (A1D_Connection_send_active[dst] > 0)
            {

                A1D_Control_flushack_active++;

                status = DCMF_Send(&A1D_Send_flush_protocol,
                                   &request[pending_count],
                                   A1D_Nocallback,
                                   DCMF_SEQUENTIAL_CONSISTENCY,
                                   dst,
                                   0,
                                   NULL,
                                   &msginfo,
                                   1);
                A1U_ERR_POP(status != DCMF_SUCCESS, "DCMF_Send returned with an error\n");
                pending_count++;

            }
            else if (A1D_Connection_put_active[dst] > 0)
            {

                src_disp = (size_t) A1D_Put_Flushcounter_ptr[A1D_Process_info.my_rank]
                         - (size_t) A1D_Membase_global[A1D_Process_info.my_rank];
                dst_disp = (size_t) A1D_Put_Flushcounter_ptr[dst]
                         - (size_t) A1D_Membase_global[dst] + 1;

                A1D_Put_flushack_active++;

                status = DCMF_Put(&A1D_Generic_put_protocol,
                                  &request[pending_count],
                                  A1D_Nocallback,
                                  DCMF_SEQUENTIAL_CONSISTENCY,
                                  dst,
                                  1,
                                  &A1D_Memregion_global[A1D_Process_info.my_rank],
                                  &A1D_Memregion_global[dst],
                                  src_disp,
                                  dst_disp,
                                  ack_callback);
                A1U_ERR_POP(status != DCMF_SUCCESS, "DCMF_Put returned with an error\n");
                pending_count++;

            }

            if (pending_count >= a1d_settings.flushall_pending_limit)
            {
                A1DI_Conditional_advance(A1D_Control_flushack_active > 0 || A1D_Put_flushack_active > 0);
                pending_count = 0;
            }

        }
    }
    A1DI_Conditional_advance(A1D_Control_flushack_active > 0 || A1D_Put_flushack_active > 0);

    A1DI_Memset((void *) A1D_Connection_send_active, 0, sizeof(uint32_t) * A1D_Process_info.num_ranks);
    A1DI_Memset((void *) A1D_Connection_put_active, 0, sizeof(uint32_t) * A1D_Process_info.num_ranks);

  fn_exit: 
    A1DI_Free(request); 
    A1U_FUNC_EXIT();
    return;

  fn_fail: 
    goto fn_exit;
}

int A1D_Flush_group(A1_group_t* group)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    if (group == A1_GROUP_WORLD || group == NULL)
    {
        A1DI_Flush_all();
        goto fn_exit;
    }
    else
    {
        A1U_ERR_POP(1, "A1D_Flush_group not implemented for non-world groups!");
        goto fn_fail;
    }

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}
