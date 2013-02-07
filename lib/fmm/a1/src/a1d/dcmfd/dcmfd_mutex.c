/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Protocol_t A1D_Mutex_protocol;

int *A1D_Mutexes_count;
A1D_Mutex_t *A1D_Mutexes;

volatile int mutex_request_active = 0;
volatile A1_bool_t mutex_acquired;

void A1DI_Mutex_callback(void *clientdata,
                         const DCMF_Control_t *info,
                         size_t peer)
{
    int status = A1_SUCCESS;
    DCMF_Control_t cmsg;
    A1D_Mutex_pkt_t *mutex_pkt = (A1D_Mutex_pkt_t *) info;
    A1D_Mutex_pkt_t response_pkt;
    A1D_Mutex_request_t *mutex_request;

    if (mutex_pkt->response == -1) 
    {
        /*This is a mutex request packet*/
        if(mutex_pkt->mutex_op == A1D_MUTEX_LOCK)
        {
            if(A1D_Mutexes[mutex_pkt->mutex_idx].mutex == -1) 
            {
                A1D_Mutexes[mutex_pkt->mutex_idx].mutex = peer;
                
                /* other fields do not matter in response packet 
                   as the requester will be waiting for just one
                   mutex response at a time */
                response_pkt.response = A1_TRUE;

                A1DI_Memcpy((void *) &cmsg,(void *) &response_pkt, sizeof(A1D_Mutex_pkt_t));

                status = DCMF_Control(&A1D_Mutex_protocol,
                                      DCMF_SEQUENTIAL_CONSISTENCY,
                                      peer,
                                      &cmsg);
                A1U_ERR_ABORT(status != DCMF_SUCCESS,
                       "DCMF_Control failed in A1DI_Mutex_callback \n"); 
            }
            else
            {
               status = A1DI_Malloc((void **) &mutex_request, sizeof(A1D_Mutex_request_t));
               A1U_ERR_ABORT(status != A1_SUCCESS,
                       "A1DI_Malloc failed in A1DI_Mutex_callback \n");  
               
               if(A1D_Mutexes[mutex_pkt->mutex_idx].tail == NULL)
               {
                  A1D_Mutexes[mutex_pkt->mutex_idx].tail = mutex_request;
                  A1D_Mutexes[mutex_pkt->mutex_idx].head = mutex_request; 
               }
               else
               {
                  A1D_Mutexes[mutex_pkt->mutex_idx].tail->next = mutex_request; 
                  A1D_Mutexes[mutex_pkt->mutex_idx].tail = mutex_request;
               }
            }
        } 
        else if(mutex_pkt->mutex_op == A1D_MUTEX_TRYLOCK)
        {
            /* other fields does not matter in response packet
               as the requestor will be waiting for just one
               mutex response at a time */
            if(A1D_Mutexes[mutex_pkt->mutex_idx].mutex == -1)
            {
                A1D_Mutexes[mutex_pkt->mutex_idx].mutex = peer;
                response_pkt.response = A1_TRUE;
            }
            else
            {
                response_pkt.response = A1_FALSE;
            }

            A1DI_Memcpy((void *) &cmsg,(void *) &response_pkt, sizeof(A1D_Mutex_pkt_t));

            status = DCMF_Control(&A1D_Mutex_protocol,
                                  DCMF_SEQUENTIAL_CONSISTENCY,
                                  peer,
                                  &cmsg);
            A1U_ERR_ABORT(status != DCMF_SUCCESS,
                   "DCMF_Control failed in A1DI_Mutex_callback \n");              
        }
        else if(mutex_pkt->mutex_op == A1D_MUTEX_UNLOCK) 
        {
            if(A1D_Mutexes[mutex_pkt->mutex_idx].mutex != peer)
            {
                A1U_ERR_ABORT(status = A1_ERROR,
                   "Invalid unlock request received \n");     
            }
 
            if(A1D_Mutexes[mutex_pkt->mutex_idx].head != NULL)
            {
                /*retrieve the next request from the queue*/
                mutex_request = A1D_Mutexes[mutex_pkt->mutex_idx].head;
                A1D_Mutexes[mutex_pkt->mutex_idx].head = A1D_Mutexes[mutex_pkt->mutex_idx].head->next;

                A1D_Mutexes[mutex_pkt->mutex_idx].mutex = mutex_request->rank; 

                response_pkt.response = A1_TRUE;

                A1DI_Memcpy((void *) &cmsg,(void *) &response_pkt, sizeof(A1D_Mutex_pkt_t));

                status = DCMF_Control(&A1D_Mutex_protocol,
                                      DCMF_SEQUENTIAL_CONSISTENCY,
                                      mutex_request->rank,
                                      &cmsg);
                A1U_ERR_ABORT(status != DCMF_SUCCESS,
                      "DCMF_Control failed in A1DI_Mutex_callback \n");               

                A1DI_Free(mutex_request); 
            }
            else
            { 
                A1D_Mutexes[mutex_pkt->mutex_idx].mutex = -1;  
            }
        }
        else
        {
            A1U_ERR_ABORT(status = A1_ERROR,
                   "Invalid mutex request received \n");

        }
    }
    else
    {
        /*this is a mutex reponse packet*/
        mutex_acquired = mutex_pkt->response;        

        mutex_request_active--; 
    }
}

int A1DI_Mutex_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Control_Configuration_t conf;

    A1U_FUNC_ENTER();

    /* Protocol used for mutex operations */
    conf.protocol = DCMF_DEFAULT_CONTROL_PROTOCOL;
    conf.network = DCMF_DEFAULT_NETWORK;
    conf.cb_recv = A1DI_Mutex_callback;
    conf.cb_recv_clientdata = NULL;

    status = DCMF_Control_register(&A1D_Mutex_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "Mutexes protocol registration returned with error %d \n",
                status);

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1D_Create_mutexes(A1_group_t* group, int mutex_count, int *mutex_count_ar)
{
    int index, rank, status = A1_SUCCESS;
    DCMF_Control_t cmsg;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    A1U_ERR_POP(status = (group != A1_GROUP_WORLD),
                "Mutexes are not implemented for non-world groups!");

    status = A1DI_Malloc((void **) &A1D_Mutexes, sizeof(A1D_Mutex_t)*mutex_count);
    A1U_ERR_POP(status != 0,
                "A1DI_Malloc returned error in A1D_Create_counter\n");

    for(index=0; index<mutex_count; index++)
    {
       A1D_Mutexes[index].mutex = -1;  
    }

    status = A1DI_Malloc((void **) &A1D_Mutexes_count, sizeof(int)*A1D_Process_info.num_ranks);
    A1U_ERR_POP(status != 0,
                "A1DI_Malloc returned error in A1D_Create_counter\n");   

    A1D_Mutexes_count[A1D_Process_info.my_rank] = mutex_count;

    /*Exchanging mutex count information among processes in the group*/
    A1D_Control_xchange_info.xchange_ptr = (void *) A1D_Mutexes_count;
    A1D_Control_xchange_info.xchange_size = sizeof(int);
    A1D_Control_xchange_info.rcv_active += A1D_Process_info.num_ranks - 1;

    A1DI_GlobalBarrier();

    A1DI_Memcpy((void *) &cmsg,
                (void *) &A1D_Mutexes_count[A1D_Process_info.my_rank],
                sizeof(int));

    for (rank = 0; rank < A1D_Process_info.num_ranks; rank++)
    {
        likely_if (rank != A1D_Process_info.my_rank)
        {
            DCMF_Control(&A1D_Control_xchange_info.protocol,
                         DCMF_SEQUENTIAL_CONSISTENCY,
                         rank,
                         &cmsg);
        }
    }

    A1DI_Conditional_advance(A1D_Control_xchange_info.rcv_active > 0); 

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1D_Destroy_mutexes(A1_group_t* group)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    A1U_ERR_POP(status = (group != A1_GROUP_WORLD),
                "Mutexes are not implemented for non-world groups!");

    A1DI_GlobalBarrier();

    A1DI_Free(A1D_Mutexes);

    A1DI_Free(A1D_Mutexes_count);

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1D_Lock_mutex(A1_group_t* group, int mutex, int proc)
{
    int status = A1_SUCCESS;
    DCMF_Control_t cmsg;
    A1D_Mutex_pkt_t packet;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    A1U_ERR_POP(status = (mutex > A1D_Mutexes_count[proc]),
                "Non-existent mutex being requested!");

    mutex_request_active++;

    packet.mutex_idx = mutex;
    packet.mutex_op = A1D_MUTEX_LOCK;
    packet.response = -1; 

    A1DI_Memcpy((void *) &cmsg,(void *) &packet, sizeof(A1D_Mutex_pkt_t));      

    status = DCMF_Control(&A1D_Mutex_protocol,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          proc,
                          &cmsg);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "DCMF_Control failed in A1D_Lock_mutex \n");

    A1DI_Conditional_advance(mutex_request_active > 0);

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1D_Trylock_mutex(A1_group_t* group, int mutex, int proc, A1_bool_t *acquired)
{
    int status = A1_SUCCESS;
    DCMF_Control_t cmsg;
    A1D_Mutex_pkt_t packet;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    A1U_ERR_POP(status = (mutex > A1D_Mutexes_count[proc]),
                "Non-existent mutex being requested!");

    mutex_request_active++;

    packet.mutex_idx = mutex;
    packet.mutex_op = A1D_MUTEX_TRYLOCK;
    packet.response = -1;

    A1DI_Memcpy((void *) &cmsg,(void *) &packet, sizeof(A1D_Mutex_pkt_t));

    status = DCMF_Control(&A1D_Mutex_protocol,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          proc,
                          &cmsg);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "DCMF_Control failed in A1D_Lock_mutex \n");

    A1DI_Conditional_advance(mutex_request_active > 0);

    *acquired = mutex_acquired;            

  fn_exit:
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1D_Unlock_mutex(A1_group_t* group, int mutex, int proc)
{
    int status = A1_SUCCESS;
    DCMF_Control_t cmsg;
    A1D_Mutex_pkt_t packet;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    A1U_ERR_POP(status = (mutex > A1D_Mutexes_count[proc]),
                "Non-existent mutex being requested!");

    packet.mutex_idx = mutex;
    packet.mutex_op = A1D_MUTEX_UNLOCK;
    packet.response = -1;

    A1DI_Memcpy((void *) &cmsg,(void *) &packet, sizeof(A1D_Mutex_pkt_t));

    status = DCMF_Control(&A1D_Mutex_protocol,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          proc,
                          &cmsg);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "DCMF_Control failed in A1D_Lock_mutex \n");

  fn_exit:
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
