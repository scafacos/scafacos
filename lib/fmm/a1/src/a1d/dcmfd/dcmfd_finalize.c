/* *- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in toplevel directory.
 */

#include "dcmfdimpl.h"

int A1D_Finalize(void)
{
    int status = A1_SUCCESS;
    int count = 0;

    A1U_FUNC_ENTER();

    /* TODO: need to unset "A1 is alive" global variable */

    A1DI_CRITICAL_ENTER();

    /*waiting for everyone*/
    status = A1DI_GlobalBarrier();
    A1U_ERR_POP(status != A1_SUCCESS, 
              "A1DI_GlobalBarrier returned with an error");

    /* Freeing request pool */
    A1DI_Request_pool_finalize();

    /* Freeing handle pool */
    A1DI_Handle_pool_finalize();

    /* Freeing buffer pool */
    A1DI_Buffer_pool_finalize();

    /* Freeing memory region pointers and local memroy region*/
    A1DI_Free(A1D_Membase_global);
    A1DI_Free(A1D_Memregion_global);

    /* Freeing conenction active counters */
    A1DI_Free((void *) A1D_Connection_send_active);
    A1DI_Free((void *) A1D_Connection_put_active);

    /* Freeing put flush local counters and pointers */
    A1DI_Free(A1D_Put_Flushcounter_ptr[A1D_Process_info.my_rank]);
    A1DI_Free(A1D_Put_Flushcounter_ptr);

    if (a1d_settings.enable_cht)
    {
        status = pthread_cancel(A1DI_CHT_pthread);
    }

    A1DI_CRITICAL_EXIT();

    /* NOTE: exit critical section before finalize since CS may not work after DCMF is terminated */

    count = DCMF_Messager_finalize();
    /* Do not issue this warning if using MPI since in that case we know DCMF
       will be initialized by MPI before A1 (assuming GA->ARMCI->A1 call path). */
    //if(!a1d_settings.mpi_active)
    //{
    //    A1U_WARNING(count == 0,
    //                "DCMF_Messager_finalize has been called more than once.");
    //}

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void A1D_Abort(int error_code, char error_message[])
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1U_ERR_ABORT(status = A1_ERROR,
            "User called A1_ABORT with error code %d, error msg: %s Program terminating abnormally \n",
            error_code,
            error_message);

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

