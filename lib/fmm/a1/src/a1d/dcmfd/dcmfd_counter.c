/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Protocol_t A1D_Counter_create_protocol;
volatile int counter_create_active;
long* counter_create_response;

DCMF_Protocol_t A1D_Counter_protocol;
volatile int counter_incr_active;
long counter_incr_response;

void A1DI_Counter_create_callback(void *clientdata,
                                  const DCMF_Control_t *info,
                                  size_t peer)
{
    A1DI_Memcpy((void *) &counter_create_response, (void *) info, sizeof(void *));
    counter_create_active--;
}

void A1DI_Counter_callback(void *clientdata,
                           const DCMF_Control_t *info,
                           size_t peer)
{
    int status = A1_SUCCESS;
    A1D_Counter_pkt_t *counter_pkt = (A1D_Counter_pkt_t *) info;

    if (counter_pkt->value_ptr == NULL)
    {
        /*This is a response packet*/
        counter_incr_response = counter_pkt->value;
        counter_incr_active--;
    }
    else
    {
        A1D_Counter_pkt_t response_pkt;
        DCMF_Control_t cmsg;
        long original;

        original = *(counter_pkt->value_ptr);
        *(counter_pkt->value_ptr) += response_pkt.value;

        response_pkt.value_ptr = NULL;
        response_pkt.value = original;

        A1DI_Memcpy(&cmsg, &response_pkt, sizeof(A1D_Counter_pkt_t));

        status = DCMF_Control(&A1D_Counter_protocol,
                              DCMF_SEQUENTIAL_CONSISTENCY,
                              peer,
                              &cmsg);
        A1U_ERR_ABORT(status != DCMF_SUCCESS,
                      "DCMF_Control failed in A1DI_Counter_callback\n");
    }
}

int A1DI_Counter_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Control_Configuration_t conf;

    A1U_FUNC_ENTER();

    /* Protocol used to create counters */
    conf.protocol = DCMF_DEFAULT_CONTROL_PROTOCOL;
    conf.network = DCMF_DEFAULT_NETWORK;
    conf.cb_recv = A1DI_Counter_create_callback;
    conf.cb_recv_clientdata = NULL;

    status = DCMF_Control_register(&A1D_Counter_create_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "DCMF_Control_register returned with error %d \n",
                status);

    /* Protocol used for counter operations */
    conf.protocol = DCMF_DEFAULT_CONTROL_PROTOCOL;
    conf.network = DCMF_DEFAULT_NETWORK;
    conf.cb_recv = A1DI_Counter_callback;
    conf.cb_recv_clientdata = NULL;

    status = DCMF_Control_register(&A1D_Counter_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "Counter protocol registartion returned with error %d \n",
                status);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_Create_counter(A1_group_t* group, A1_counter_t *counter_ptr)
{
    int index, status = A1_SUCCESS;
    DCMF_Control_t cmsg;
    A1D_Counter_t *a1d_counter;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    A1U_ERR_POP(status = (group != A1_GROUP_WORLD),
                "Counters are not implemented for non-world groups!");

    status = A1DI_Malloc((void **) &a1d_counter, sizeof(A1D_Counter_t));
    A1U_ERR_POP(status != 0,
                "A1DI_Malloc returned error in A1D_Create_counter\n");

    /* TODO: We have to find a way to select the location of counters dynamically. Will 
     distributing counters across processes in a round-robin fashion be a good way? */
    a1d_counter->rank = A1D_Process_info.num_ranks - 1;

    unlikely_if (a1d_counter->rank == A1D_Process_info.my_rank)
    {

        a1d_counter->value_ptr = NULL;

        for (index = 0; index < A1D_Process_info.num_ranks; index++)
        {
            if (index != A1D_Process_info.my_rank)
            {
                status = DCMF_Control(&A1D_Counter_create_protocol,
                                      DCMF_SEQUENTIAL_CONSISTENCY,
                                      index,
                                      &cmsg);
                A1U_ERR_POP(status != DCMF_SUCCESS,
                            "DCMF_Control failed in A1D_Alloc_counter\n");
            }
        }

    }
    else
    {
        counter_create_active++;

        A1DI_Conditional_advance(counter_create_active > 0);

        a1d_counter->value_ptr = counter_create_response;
    }

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1D_Destroy_counter(A1_group_t* group, A1_counter_t *counter_ptr)
{
    int index, status = A1_SUCCESS;
    DCMF_Control_t cmsg;
    A1D_Counter_t *a1d_counter;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    A1U_ERR_POP(status = (group != A1_GROUP_WORLD),
                "Counters are not implemented for non-world groups!");

    A1DI_GlobalBarrier();

    a1d_counter = *counter_ptr;

    A1DI_Free(a1d_counter);

    *counter_ptr = NULL;

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_Incr_counter(A1_counter_t counter, long increment, long* original)
{
    int status = A1_SUCCESS;
    A1D_Counter_t *a1d_counter;
    A1D_Counter_pkt_t counter_pkt;
    DCMF_Control_t cmsg;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_counter = (A1D_Counter_t *) counter;

    if (a1d_counter->rank == A1D_Process_info.my_rank)
    {
        *original = a1d_counter->value;
        a1d_counter->value = a1d_counter->value + increment;
    }
    else
    {

        counter_pkt.value_ptr = a1d_counter->value_ptr;
        counter_pkt.value = increment;

        A1DI_Memcpy(&cmsg, &counter_pkt, sizeof(A1D_Counter_pkt_t));

        counter_incr_active = 1;

        status = DCMF_Control(&A1D_Counter_protocol,
                              DCMF_SEQUENTIAL_CONSISTENCY,
                              a1d_counter->rank,
                              &cmsg);
        A1U_ERR_POP(status != DCMF_SUCCESS,
                    "DCMF_Control failed in A1D_Incr_counter\n");

        A1DI_Conditional_advance(counter_incr_active > 0);

        *original = counter_incr_response;

    }

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

