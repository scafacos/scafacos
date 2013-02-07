/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

int A1DI_Direct_putaccv(int target,
                        A1_iov_t *iov_ar,
                        int ar_len,
                        A1_datatype_t a1_type,
                        void *scaling,
                        A1D_Handle_t *a1d_handle)
{
    int i, j, status = A1_SUCCESS;
    A1D_Putacc_header_t header;
    A1D_Request_t *a1d_request;
    DCMF_Callback_t done_callback;

    A1U_FUNC_ENTER();

    header.datatype = a1_type;
    switch (a1_type)
    {
        case A1_DOUBLE:
            (header.scaling).double_value = *((double *) scaling);
            break;
        case A1_INT32:
            (header.scaling).int32_value = *((int32_t *) scaling);
            break;
        case A1_INT64:
            (header.scaling).int64_value = *((int64_t *) scaling);
            break;
        case A1_UINT32:
            (header.scaling).uint32_value = *((uint32_t *) scaling);
            break;
        case A1_UINT64:
            (header.scaling).uint64_value = *((uint64_t *) scaling);
            break;
        case A1_FLOAT:
            (header.scaling).float_value = *((float *) scaling);
            break;
        default:
            status = A1_ERROR;
            A1U_ERR_POP((status != A1_SUCCESS),"Invalid data type in putacc \n");
            break;
    }

    for (i=0; i<ar_len; i++)
    {
        for(j=0; j<iov_ar[i].ptr_ar_len; j++)
        {

           a1d_request = A1DI_Get_request(1);
           A1U_ERR_POP(status = (a1d_request == NULL),
                "A1DI_Get_request returned error.\n");
           A1DI_Set_handle(a1d_request, a1d_handle);

           done_callback.function = A1DI_Request_done;
           done_callback.clientdata = (void *) a1d_request;
 
           a1d_handle->active++;

           header.target_ptr = iov_ar[i].target_ptr_ar[j];
 
           status = DCMF_Send(&A1D_Generic_putacc_protocol,
                              &(a1d_request->request),
                              done_callback,
                              DCMF_SEQUENTIAL_CONSISTENCY,
                              target,
                              iov_ar[i].size,
                              iov_ar[i].source_ptr_ar[j],
                              (DCQuad *) &header,
                              (unsigned) 2);
           A1U_ERR_POP((status != DCMF_SUCCESS), "Putacc returned with an error \n");
 
           A1D_Connection_send_active[target]++;
        }
    }

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1D_PutAccV(int target,
                A1_iov_t *iov_ar,
                int ar_len,
                A1_datatype_t a1_type,
                void* scaling)
{
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_handle = A1DI_Get_handle();
    A1U_ERR_POP(status = (a1d_handle == NULL),
                "A1DI_Get_handle returned NULL in A1D_PutAccS.\n");

    status = A1DI_Direct_putaccv(target,
                                 iov_ar,
                                 ar_len,
                                 a1_type,
                                 scaling,
                                 a1d_handle);
    A1U_ERR_POP(status, "Direct putaccv function returned with an error \n");

    A1DI_Conditional_advance(a1d_handle->active > 0);

  fn_exit:
    A1DI_Release_handle(a1d_handle);
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1D_NbPutAccV(int target,
                  A1_iov_t *iov_ar,
                  int ar_len,
                  A1_datatype_t a1_type,
                  void* scaling,
                  A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_handle = (A1D_Handle_t *) a1_handle;

    status = A1DI_Direct_putaccv(target,
                                 iov_ar,
                                 ar_len,
                                 a1_type,
                                 scaling,
                                 a1d_handle);
    A1U_ERR_POP(status, "Direct putaccv function returned with an error \n");

  fn_exit:
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
