/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Configure_t A1D_Messager_info;
A1D_Process_info_t A1D_Process_info;

DCMF_Callback_t A1D_Nocallback;

int A1D_Initialize(int thread_level)
{

    int status = A1_SUCCESS;
    int count = 0;

    A1U_FUNC_ENTER();

    /* TODO: need a non-DCMF lock here to make this function thread-safe */
    /* TODO: need to set "A1 is alive" global variable */

    count = DCMF_Messager_initialize();
    /* Do not issue this warning if using MPI since in that case we know DCMF
       will be initialized by MPI before A1 (assuming GA->ARMCI->A1 call path). */
    // OMIT THIS WARNING FOR NOW.  WE HAVE BIGGER PROBLEMS AT THE MOMENT.
    //if(!a1d_settings.mpi_active)
    //{
    //    A1U_WARNING(count == 0,
    //                "DCMF_Messager_initialize has been called more than once.");
    //}

    if ( DCMF_Messager_size() > 1 ) DCMF_Collective_initialize();

    A1D_Nocallback.function = NULL;
    A1D_Nocallback.clientdata = NULL;

    status = A1DI_Read_parameters();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Read_parameters returned with error \n");

    if (a1d_settings.enable_cht)
    {
        if (!a1d_settings.mpi_active)
        {
            /* We can use THREAD_SERIALIZED if we are implementing out own locks
             * ~AND~ MPI is not active */
            A1D_Messager_info.thread_level = DCMF_THREAD_SERIALIZED;
        }
        else
        {
            /* If MPI is active, DCMF_Critical_section requires DCMF_THREAD_MULTIPLE
             * to work properly */
            A1D_Messager_info.thread_level = DCMF_THREAD_MULTIPLE;
        }
        A1D_Messager_info.interrupts = DCMF_INTERRUPTS_OFF;
    }
    else
    {
        switch (thread_level)
        {
        case A1_THREAD_SINGLE:
            thread_level = DCMF_THREAD_SINGLE;
            break;
        case A1_THREAD_FUNNELED:
            thread_level = DCMF_THREAD_FUNNELED;
            break;
        case A1_THREAD_SERIALIZED:
            thread_level = DCMF_THREAD_SERIALIZED;
            break;
        case A1_THREAD_MULTIPLE:
            thread_level = DCMF_THREAD_MULTIPLE;
            break;
        default:
            A1U_ERR_POP(A1_ERROR,
                        "Unsupported thread level provided in A1D_Initialize \n");
            break;
        }
        A1D_Messager_info.interrupts = ( a1d_settings.enable_interrupts ? DCMF_INTERRUPTS_ON : DCMF_INTERRUPTS_OFF );
    }

    status = DCMF_Messager_configure(&A1D_Messager_info, &A1D_Messager_info);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "DCMF_Messager_configure returned with error \n");

    A1D_Process_info.my_rank = DCMF_Messager_rank();
    A1D_Process_info.num_ranks = DCMF_Messager_size();

    /* TODO: initialize node rank/size properly on BGP */
    A1D_Process_info.my_node = DCMF_Messager_rank();
    A1D_Process_info.num_nodes = DCMF_Messager_size();

    if (a1d_settings.enable_cht)
    {
        /* Initialize LockBox if it is the locking mechanism used */
        // A1DI_GLOBAL_LBMUTEX_INITIALIZE();

        /* Create CHT */
        status = pthread_create(&A1DI_CHT_pthread,
                                NULL,
                                &A1DI_CHT_advance_lock,
                                NULL);
        A1U_ERR_POP(status != 0, "pthread_create returned with error \n");
    }

    A1DI_CRITICAL_ENTER();

    status = A1DI_Print_parameters();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Print_parameters returned with error \n");

    status = A1DI_Control_xchange_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Control_xchange_initialize returned with error \n");

    status = A1DI_Control_flushack_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Control_flushack_initialize returned with error \n");

    status = A1DI_GlobalBarrier_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_GlobalBarrier_initialize returned with error \n");

    status = A1DI_GlobalAllreduce_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_GlobalAllreduce_initialize returned with error \n");

    status = A1DI_GlobalBcast_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_GlobalBcast_initialize returned with error \n");

    status = A1DI_Put_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Put_initialize returned with error \n");

    status = A1DI_Get_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Get_initialize returned with error \n");

    status = A1DI_Putacc_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Putacc_initialize returned with error \n");

    status = A1DI_Putmod_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Putmod_initialize returned with error \n");

    status = A1DI_Packed_puts_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Packed_puts_initialize returned with error \n");

    status = A1DI_Packed_gets_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Packed_gets_initialize returned with error \n");

    status = A1DI_Packed_putaccs_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Packed_putaccs_initialize returned with error \n");

    status = A1DI_Rmw_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Rmw_initialize returned with error \n");

    status = A1DI_Counter_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Counter_initialize returned with error \n");

    status = A1DI_Send_flush_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Send_flush_initialize returned with error \n");

    status = A1DI_Put_flush_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "Put flush initialize returned with error \n");

    status = A1DI_Request_pool_initialize();
    A1U_ERR_POP(status != A1_SUCCESS, "A1DI_Request_pool_initialize failed \n");

    status = A1DI_Buffer_pool_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Buffer_pool_initialize returned with error \n");

    status = A1DI_Handle_pool_initialize();
    A1U_ERR_POP(status != A1_SUCCESS, "A1DI_Handle_pool_initialize failed \n");

    /* TODO: Do we need to barrier before this call?
     *       Won't this call fail internally if one process is late? */
    /* Resolution: This function has a barrier inside, before the exchange
     *       happens. So a barrier is not required here */
    status = A1DI_Memregion_Global_initialize();
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Memregion_Global_initialize returned with error \n");
  
    /*waiting for everyone*/
    status = A1DI_GlobalBarrier();
    A1U_ERR_POP(status != A1_SUCCESS,
              "A1DI_GlobalBarrier returned with an error");

  fn_exit: 
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

