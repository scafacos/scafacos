/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Protocol_t A1D_Generic_putmod_protocol;

void A1DI_RecvDone_putmod_callback(void *clientdata, DCMF_Error_t *error)
{
    int status = A1_SUCCESS;

    A1D_Request_t *a1d_request = (A1D_Request_t *) clientdata;
    A1D_Putmod_header_t *header = (A1D_Putmod_header_t *) a1d_request->buffer_ptr;
    void* source_ptr = (void *) ((size_t) a1d_request->buffer_ptr + sizeof(A1D_Putmod_header_t));
    int data_size = a1d_request->buffer_size - sizeof(A1D_Putmod_header_t);

    switch (header->op)
    {
        case A1_BXOR:
            switch (header->datatype)
            {
                case A1_INT32:
                    A1DI_MOD_BXOR(int32_t,
                            source_ptr,
                            header->target_ptr,
                            data_size / sizeof(int32_t));
                    break;
                case A1_INT64:
                    A1DI_MOD_BXOR(int64_t,
                            source_ptr,
                            header->target_ptr,
                            data_size / sizeof(int64_t));
                    break;
                case A1_UINT32:
                    A1DI_MOD_BXOR(uint32_t,
                            source_ptr,
                            header->target_ptr,
                            data_size / sizeof(uint32_t));
                    break;
                case A1_UINT64:
                    A1DI_MOD_BXOR(uint64_t,
                            source_ptr,
                            header->target_ptr,
                            data_size / sizeof(uint64_t));
                    break;
                default:
                    A1U_ERR_ABORT(status, "Invalid datatype received in PutdodV\n");
                    break;
            }
            break;
        default:
            status = A1_ERROR;
            A1U_ERR_ABORT((status != A1_SUCCESS), "Invalid op type in PutdodV\n");
            break;
    }

    A1DI_Release_request(a1d_request);
}

DCMF_Request_t* A1DI_RecvSend_putmod_callback(void *clientdata,
                                              const DCQuad *msginfo,
                                              unsigned count,
                                              size_t peer,
                                              size_t sndlen,
                                              size_t *rcvlen,
                                              char **rcvbuf,
                                              DCMF_Callback_t *cb_done)
{
    int status = 0;
    A1D_Request_t *a1d_request;

    a1d_request = A1DI_Get_request(0);
    A1U_ERR_ABORT(status = (a1d_request == NULL),
            "A1DI_Get_request returned NULL in A1DI_RecvSend_putmod_callback\n");

    A1U_ASSERT_ABORT(sizeof(A1D_Putmod_header_t) == count*sizeof(DCQuad),
            "Header of invalid size received in A1DI_RecvSend_putmod_callback\n")

    a1d_request->buffer_size = sndlen + sizeof(A1D_Putmod_header_t);
    status = A1DI_Malloc((void **) &(a1d_request->buffer_ptr),
                         sndlen + sizeof(A1D_Putmod_header_t));
    A1U_ERR_ABORT(status != 0,
            "A1DI_Malloc failed in A1DI_RecvSend_packedputs_callback\n");

    A1DI_Memcpy(a1d_request->buffer_ptr,(void *) msginfo,sizeof(A1D_Putmod_header_t));

    *rcvlen = sndlen;
    *rcvbuf = (char *) ((size_t) a1d_request->buffer_ptr + sizeof(A1D_Putmod_header_t));

    cb_done->function = A1DI_RecvDone_putmod_callback;
    cb_done->clientdata = (void *) a1d_request;

    return &(a1d_request->request);
}

void A1DI_RecvSendShort_putmod_callback(void *clientdata,
                                        const DCQuad *msginfo,
                                        unsigned count,
                                        size_t peer,
                                        const char *src,
                                        size_t bytes)
{
    int status = A1_SUCCESS;
    A1D_Putmod_header_t *header = (A1D_Putmod_header_t *) msginfo;
    void* source_ptr = (void *) src;

    switch (header->op)
    {
        case A1_BXOR:
            switch (header->datatype)
            {
                case A1_INT32:
                    A1DI_MOD_BXOR(int32_t,
                            source_ptr,
                            header->target_ptr,
                            bytes / sizeof(int32_t));
                    break;
                case A1_INT64:
                    A1DI_MOD_BXOR(int64_t,
                            source_ptr,
                            header->target_ptr,
                            bytes / sizeof(int64_t));
                    break;
                case A1_UINT32:
                    A1DI_MOD_BXOR(uint32_t,
                            source_ptr,
                            header->target_ptr,
                            bytes / sizeof(uint32_t));
                    break;
                case A1_UINT64:
                    A1DI_MOD_BXOR(uint64_t,
                            source_ptr,
                            header->target_ptr,
                            bytes / sizeof(uint64_t));
                    break;
                default:
                    A1U_ERR_ABORT(status, "Invalid datatype received in PutdodV\n");
                    break;
            }
            break;
        default:
            status = A1_ERROR;
            A1U_ERR_ABORT((status != A1_SUCCESS), "Invalid op type in PutdodV\n");
            break;
    }
}

int A1DI_Putmod_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Send_Configuration_t conf;

    A1U_FUNC_ENTER();

    conf.protocol = DCMF_DEFAULT_SEND_PROTOCOL;
    conf.network = DCMF_TORUS_NETWORK;
    conf.cb_recv_short = A1DI_RecvSendShort_putmod_callback;
    conf.cb_recv_short_clientdata = NULL;
    conf.cb_recv = A1DI_RecvSend_putmod_callback;
    conf.cb_recv_clientdata = NULL;

    status = DCMF_Send_register(&A1D_Generic_putmod_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
            "DCMF_Send_register registration returned with error %d \n",
            status);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}


int A1DI_Direct_putmodv(int target,
                        A1_iov_t *iov_ar,
                        int ar_len,
                        A1_reduce_op_t a1_op,
                        A1_datatype_t a1_type,
                        A1D_Handle_t *a1d_handle)
{
    int i, j, status = A1_SUCCESS;
    A1D_Putmod_header_t header;
    A1D_Request_t *a1d_request;
    DCMF_Callback_t done_callback;

    A1U_FUNC_ENTER();

    header.datatype = a1_type;
    header.op = a1_op;

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
 
           status = DCMF_Send(&A1D_Generic_putmod_protocol,
                              &(a1d_request->request),
                              done_callback,
                              DCMF_SEQUENTIAL_CONSISTENCY,
                              target,
                              iov_ar[i].size,
                              iov_ar[i].source_ptr_ar[j],
                              (DCQuad *) &header,
                              (unsigned) 2);
           A1U_ERR_POP((status != DCMF_SUCCESS), "DCMF_Send returned with an error \n");
 
           A1D_Connection_send_active[target]++;
        }
    }

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1D_PutModV(int target,
                A1_iov_t *iov_ar,
                int ar_len,
                A1_reduce_op_t a1_op,
                A1_datatype_t a1_type)
{
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_handle = A1DI_Get_handle();
    A1U_ERR_POP(status = (a1d_handle == NULL),
                "A1DI_Get_handle returned NULL in A1D_PutModV.\n");

    status = A1DI_Direct_putmodv(target,
                                 iov_ar,
                                 ar_len,
                                 a1_op,
                                 a1_type,
                                 a1d_handle);
    A1U_ERR_POP(status, "A1DI_Direct_putmodv returned with an error \n");

    A1DI_Conditional_advance(a1d_handle->active > 0);

  fn_exit:
    A1DI_Release_handle(a1d_handle);
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}
