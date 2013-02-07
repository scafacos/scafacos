/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

int A1D_Wait_handle_all()
{
    int status = A1_SUCCESS;
    int index;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    for (index = 0; index < a1d_settings.handlepool_size; index++)
    {
        if (A1D_Active_handle_list[index] != NULL)
            A1DI_Conditional_advance((A1D_Active_handle_list[index])->active > 0);
    }

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_Wait_handle(A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_handle = (A1D_Handle_t *) a1_handle;

    A1DI_Conditional_advance(a1d_handle->active > 0);

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_Wait_handle_list(int count, A1_handle_t *a1_handle)
{
    int status = A1_SUCCESS;
    int index;
    A1D_Handle_t *a1d_handle;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    for (index = 0; index < count; index++)
    {
        /* TODO: Do we really need to do a copy here?
         *        Isn't it enough to pass the argument as const?
         */
        a1d_handle = (A1D_Handle_t *) a1_handle[index];
        A1DI_Conditional_advance(a1d_handle->active > 0);
    }

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_Test_handle(A1_handle_t a1_handle, A1_bool_t* completed)
{
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    A1DI_Advance();
    /* TODO: Do we really need to do a copy here?
     *        Isn't it enough to pass the argument as const?
     */
    a1d_handle = (A1D_Handle_t *) a1_handle;
    *completed = (a1d_handle->active > 0) ? A1_FALSE : A1_TRUE;

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_Test_handle_list(int count,
                         A1_handle_t *a1_handle,
                         A1_bool_t* *completed)
{
    int i;
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    /* TODO: verify that this is the correct implementation */
    A1DI_Advance();
    for (i = 0; i < count; i++)
    {
        /* TODO: Do we really need to do a copy here?
         *        Isn't it enough to pass the argument as const?
         */
        a1d_handle = (A1D_Handle_t *) a1_handle[i];
        *completed[i] = (a1d_handle->active > 0) ? A1_FALSE : A1_TRUE;
    }

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}
