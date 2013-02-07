/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

A1D_Request_pool_t A1D_Request_pool;

A1D_Request_t* A1DI_Get_request(int wait_and_advance)
{
    int status = A1_SUCCESS;
    A1D_Request_t *a1d_request = NULL;

    A1U_FUNC_ENTER();

    if(!A1D_Request_pool.head)
         A1U_DEBUG_PRINT("Request pool exhausted. Wait and advance is: %d \n", 
                   wait_and_advance);

    if (!wait_and_advance && !A1D_Request_pool.head)
    {
        status = A1DI_Malloc((void **) &a1d_request, sizeof(A1D_Request_t));
        A1U_ERR_ABORT(status != 0,
                   "A1DI_Malloc failed while allocating request pool\
                        in A1DI_Get_request\n");
        a1d_request->in_pool = 0;
    }
    else
    {
        A1DI_Conditional_advance(!A1D_Request_pool.head);

        a1d_request = A1D_Request_pool.head;
        A1D_Request_pool.head = A1D_Request_pool.head->next;
        a1d_request->in_pool = 1;
    }

    a1d_request->next = NULL;
    a1d_request->buffer_ptr = NULL;
    a1d_request->a1d_buffer_ptr = NULL;
    a1d_request->handle_ptr = NULL;

    fn_exit: A1U_FUNC_EXIT();
    return a1d_request;

    fn_fail: goto fn_exit;
}

void A1DI_Release_request(A1D_Request_t *a1d_request)
{
    A1U_FUNC_ENTER();
    A1D_Handle_t *a1d_handle;

    if (a1d_request->a1d_buffer_ptr != NULL)
    {
        A1DI_Release_buffer(a1d_request->a1d_buffer_ptr);
    }

    if ((a1d_request->buffer_ptr) != NULL)
    {
        A1DI_Free((void *) (a1d_request->buffer_ptr));
    }

    if (a1d_request->handle_ptr != NULL)
    {
        a1d_handle = a1d_request->handle_ptr;
        --(a1d_handle->active);
    }

    if (a1d_request->in_pool == 0)
    {
        A1DI_Free(a1d_request);
    }
    else
    {
        a1d_request->next = A1D_Request_pool.head;
        A1D_Request_pool.head = a1d_request;
    }

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

int A1DI_Request_pool_initialize()
{

    int status = A1_SUCCESS;
    int index;
    A1D_Request_t *a1d_request;

    A1U_FUNC_ENTER();

    status = A1DI_Malloc((void **) &(A1D_Request_pool.region_ptr),
                                 sizeof(A1D_Request_t)
                                         * a1d_settings.requestpool_size);
    A1U_ERR_POP(status != A1_SUCCESS,
                "A1DI_Malloc failed while allocating request pool\
                      in A1DI_Request_pool_initialize\n");

    a1d_request = A1D_Request_pool.region_ptr;
    A1D_Request_pool.head = a1d_request;
    for (index = 0; index < a1d_settings.requestpool_size - 1; index++)
    {
        a1d_request[index].next = &a1d_request[index + 1];
    }
    a1d_request[index].next = NULL;

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

void A1DI_Request_pool_finalize()
{
    int i;

    A1U_FUNC_ENTER();

    A1DI_Free((void *) (A1D_Request_pool.region_ptr));

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}
