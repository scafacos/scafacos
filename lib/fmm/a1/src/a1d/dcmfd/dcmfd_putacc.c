/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Protocol_t A1D_Generic_putacc_protocol;

void A1DI_RecvDone_putacc_callback(void *clientdata, DCMF_Error_t *error)
{
    int status = A1_SUCCESS;

    A1D_Request_t *a1d_request = (A1D_Request_t *) clientdata;
    A1D_Putacc_header_t *header = (A1D_Putacc_header_t *) a1d_request->buffer_ptr;
    void* source_ptr = (void *) ((size_t) a1d_request->buffer_ptr + sizeof(A1D_Putacc_header_t));
    int data_size = a1d_request->buffer_size - sizeof(A1D_Putacc_header_t);

    switch (header->datatype)
    {
        case A1_DOUBLE:
            A1DI_ACC(double,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).double_value,
                    data_size/sizeof(double));
            break;
        case A1_INT32:
            A1DI_ACC(int32_t,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).int32_value,
                    data_size / sizeof(int32_t));
            break;
        case A1_INT64:
            A1DI_ACC(int64_t,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).int64_value,
                    data_size / sizeof(int64_t));
            break;
        case A1_UINT32:
            A1DI_ACC(uint32_t,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).uint32_value,
                    data_size / sizeof(uint32_t));
            break;
        case A1_UINT64:
            A1DI_ACC(uint64_t,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).uint64_value,
                    data_size / sizeof(uint64_t));
            break;
        case A1_FLOAT:
            A1DI_ACC(float,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).float_value,
                    data_size / sizeof(float));
            break;
        default:
            A1U_ERR_ABORT(status, "Invalid datatype received in Putacc operation\n");
            break;
    }

    A1DI_Release_request(a1d_request);
}

DCMF_Request_t* A1DI_RecvSend_putacc_callback(void *clientdata,
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
            "A1DI_Get_request returned NULL in A1DI_RecvSend_putacc_callback\n");

    A1U_ASSERT_ABORT(sizeof(A1D_Putacc_header_t) == count*sizeof(DCQuad),
            "Header of invalid size received in A1DI_RecvSend_putacc_callback\n")

    a1d_request->buffer_size = sndlen + sizeof(A1D_Putacc_header_t);
    status = A1DI_Malloc((void **) &(a1d_request->buffer_ptr),
                         sndlen + sizeof(A1D_Putacc_header_t));
    A1U_ERR_ABORT(status != 0,
            "A1DI_Malloc failed in A1DI_RecvSend_packedputs_callback\n");

    A1DI_Memcpy(a1d_request->buffer_ptr,(void *) msginfo,sizeof(A1D_Putacc_header_t));

    *rcvlen = sndlen;
    *rcvbuf = (char *) ((size_t) a1d_request->buffer_ptr + sizeof(A1D_Putacc_header_t));

    cb_done->function = A1DI_RecvDone_putacc_callback;
    cb_done->clientdata = (void *) a1d_request;

    return &(a1d_request->request);
}

void A1DI_RecvSendShort_putacc_callback(void *clientdata,
                                        const DCQuad *msginfo,
                                        unsigned count,
                                        size_t peer,
                                        const char *src,
                                        size_t bytes)
{
    int status = A1_SUCCESS;
    A1D_Putacc_header_t *header = (A1D_Putacc_header_t *) msginfo;
    void* source_ptr = (void *) src;

    switch (header->datatype)
    {
        case A1_DOUBLE:
            A1DI_ACC(double,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).double_value,
                    bytes/sizeof(double));
            break;
        case A1_INT32:
            A1DI_ACC(int32_t,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).int32_value,
                    bytes / sizeof(int32_t));
            break;
        case A1_INT64:
            A1DI_ACC(int64_t,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).int64_value,
                    bytes / sizeof(int64_t));
            break;
        case A1_UINT32:
            A1DI_ACC(uint32_t,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).uint32_value,
                    bytes / sizeof(uint32_t));
            break;
        case A1_UINT64:
            A1DI_ACC(uint64_t,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).uint64_value,
                    bytes / sizeof(uint64_t));
            break;
        case A1_FLOAT:
            A1DI_ACC(float,
                    source_ptr,
                    header->target_ptr,
                    (header->scaling).float_value,
                    bytes/sizeof(float));
            break;
        default:
            A1U_ERR_ABORT(status, "Invalid datatype received in Putacc operation \n")
            ;
            break;
    }
}

int A1DI_Putacc_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Send_Configuration_t conf;

    A1U_FUNC_ENTER();

    conf.protocol = DCMF_DEFAULT_SEND_PROTOCOL;
    conf.network = DCMF_TORUS_NETWORK;
    conf.cb_recv_short = A1DI_RecvSendShort_putacc_callback;
    conf.cb_recv_short_clientdata = NULL;
    conf.cb_recv = A1DI_RecvSend_putacc_callback;
    conf.cb_recv_clientdata = NULL;

    status = DCMF_Send_register(&A1D_Generic_putacc_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
            "DCMF_Send_register registration returned with error %d \n",
            status);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_PutAcc(int target,
               void* source_ptr,
               void* target_ptr,
               int bytes,
               A1_datatype_t a1_type,
               void* scaling)
{
    int status = A1_SUCCESS;
    DCMF_Request_t request;
    DCMF_Callback_t done_callback;
    volatile int active;
    A1D_Putacc_header_t header;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    done_callback.function = A1DI_Generic_done;
    done_callback.clientdata = (void *) &active;
    active = 1;

    header.target_ptr = target_ptr;
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
            A1U_ERR_POP((status != A1_SUCCESS), "Invalid data type in putacc \n");
            break;
    }

    status = DCMF_Send(&A1D_Generic_putacc_protocol,
                       &request,
                       done_callback,
                       DCMF_SEQUENTIAL_CONSISTENCY,
                       target,
                       bytes,
                       source_ptr,
                       (DCQuad *) &header,
                       (unsigned) 2);
    A1U_ERR_POP((status != A1_SUCCESS), "DCMF_Send returned with an error \n");

    A1D_Connection_send_active[target]++;
    A1DI_Conditional_advance(active > 0);

    fn_exit:
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_NbPutAcc(int target,
                 void* source_ptr,
                 void* target_ptr,
                 int bytes,
                 A1_datatype_t a1_type,
                 void* scaling,
                 A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle;
    DCMF_Callback_t done_callback;
    A1D_Request_t *a1d_request;
    A1D_Putacc_header_t header;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_handle = (A1D_Handle_t *) a1_handle;

    a1d_handle->active++;

    a1d_request = A1DI_Get_request(1);
    A1U_ERR_POP(status = (a1d_request == NULL),
            "A1DI_Get_request returned NULL request \n");
    A1DI_Set_handle(a1d_request, a1d_handle);

    done_callback.function = A1DI_Request_done;
    done_callback.clientdata = (void *) a1d_request;

    header.target_ptr = target_ptr;
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
            A1U_ERR_POP(status, "Invalid data type in putacc \n");
            break;
    }

    status = DCMF_Send(&A1D_Generic_putacc_protocol,
                       &(a1d_request->request),
                       done_callback,
                       DCMF_SEQUENTIAL_CONSISTENCY,
                       target,
                       bytes,
                       source_ptr,
                       (DCQuad *) &header,
                       (unsigned) 2);
    A1U_ERR_POP((status != DCMF_SUCCESS), "DCMF_Send returned with an error \n");

    A1D_Connection_send_active[target]++;

    fn_exit:
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}
