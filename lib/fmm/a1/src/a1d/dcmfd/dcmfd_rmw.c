/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Protocol_t A1D_Rmw_protocol;
DCMF_Protocol_t A1D_Rmw_response_protocol;

void A1DI_RecvDone_rmw_response_callback(void *clientdata,
                                         DCMF_Error_t *error)
{
    int status = A1_SUCCESS;
    A1D_Request_t *a1d_request =  (A1D_Request_t *) clientdata;
    A1D_Buffer_t *a1d_buffer = a1d_request->a1d_buffer_ptr;
    A1D_Rmw_response_header_t *header = (A1D_Rmw_response_header_t *) a1d_buffer->buffer_ptr;
    A1D_Handle_t *a1d_handle = (A1D_Handle_t *) header->handle_ptr;
    void *reponse_data = (void *) (a1d_buffer->buffer_ptr + sizeof(A1D_Rmw_response_header_t));
    void *source_ptr_out = header->source_ptr_out;

    A1DI_Memcpy(source_ptr_out, reponse_data, header->bytes); 

    (a1d_handle->active)--;

    A1DI_Release_request(a1d_request);
}

DCMF_Request_t* A1DI_RecvSend_rmw_response_callback(void *clientdata,
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
    A1D_Buffer_t *a1d_buffer;

    a1d_request = A1DI_Get_request(0);
    A1U_ERR_ABORT(status = (a1d_request == NULL),
                "A1DI_Get_request returned NULL in A1DI_RecvSend_rmw_response_callback.\n");

    a1d_buffer = A1DI_Get_buffer(sndlen + sizeof(A1D_Rmw_response_header_t), 0);
    A1U_ERR_ABORT(status = (a1d_buffer == NULL),
                  "A1DI_Get_buffer returned NULL.\n");

    A1DI_Memcpy(a1d_buffer->buffer_ptr, msginfo, sizeof(DCQuad)*count);

    *rcvlen = sndlen;
    *rcvbuf = a1d_buffer->buffer_ptr + sizeof(A1D_Rmw_response_header_t);

    a1d_request->a1d_buffer_ptr = a1d_buffer;
    cb_done->function = A1DI_RecvDone_rmw_response_callback;
    cb_done->clientdata = (void *) a1d_request;

    return &(a1d_request->request);
}

void A1DI_RecvSendShort_rmw_response_callback(void *clientdata,
                                              const DCQuad *msginfo,
                                              unsigned count,
                                              size_t peer,
                                              const char *src,
                                              size_t bytes)
{
    int status = A1_SUCCESS;
    A1D_Rmw_response_header_t *header = (A1D_Rmw_response_header_t *) msginfo;
    A1D_Handle_t *a1d_handle = (A1D_Handle_t *) header->handle_ptr;
    void *response_data = (void *) src;
    void *source_ptr_out = header->source_ptr_out;
    
    A1DI_Memcpy(source_ptr_out, response_data, header->bytes);

    (a1d_handle->active)--;

}

void A1DI_RecvDone_rmw_callback(void *clientdata,
                                DCMF_Error_t *error)
{
    int status = A1_SUCCESS;
    DCMF_Callback_t done_callback;
    A1D_Request_t *response_request = NULL;
    A1D_Rmw_response_header_t response_header;
    A1D_Request_t *a1d_request =  (A1D_Request_t *) clientdata;
    A1D_Buffer_t *a1d_buffer = a1d_request->a1d_buffer_ptr;
    A1D_Rmw_header_t *header = (A1D_Rmw_header_t *) a1d_buffer->buffer_ptr;
    A1_datatype_t datatype = header->datatype;
    void *source = (void *) (a1d_buffer->buffer_ptr + sizeof(A1D_Rmw_header_t));
    void *target = header->target_ptr;
    void *original = NULL;

    status = A1DI_Malloc(&original, header->bytes);
    A1U_ERR_ABORT(status,
                "A1DI_Malloc returned error in A1DI_RecvDone_rmw_callback\n"); 
 
    if(header->op == A1_FETCH_AND_ADD)
    {
       likely_if(datatype == A1_DOUBLE) 
       {
           A1DI_FETCHANDADD_EXECUTE(double, source, target, original, header->bytes/sizeof(double));
       }
       else
       {
          switch (datatype)
          {
          case A1_INT32:
              A1DI_FETCHANDADD_EXECUTE(int32_t, source, target, original, header->bytes/sizeof(int32_t));
              break;
          case A1_INT64:
              A1DI_FETCHANDADD_EXECUTE(int64_t, source, target, original, header->bytes/sizeof(int64_t));
              break; 
          case A1_UINT32:
              A1DI_FETCHANDADD_EXECUTE(uint32_t, source, target, original, header->bytes/sizeof(uint32_t));
              break;
          case A1_UINT64:
              A1DI_FETCHANDADD_EXECUTE(uint64_t, source, target, original, header->bytes/sizeof(uint64_t));
              break;
          case A1_FLOAT:
              A1DI_FETCHANDADD_EXECUTE(float, source, target, original, header->bytes/sizeof(float));
              break;
          default:
              status = A1_ERROR;
              A1U_ERR_ABORT((status != A1_SUCCESS), "Invalid data type in rmw \n");
              break;
          }
       }    
    }
    else if(header->op = A1_SWAP)
    {
        A1DI_Memcpy(original, target, header->bytes); 
        A1DI_Memcpy(target, source, header->bytes); 
    }

    response_request = A1DI_Get_request(0);
    A1U_ERR_ABORT(status = (a1d_request == NULL), "A1DI_Get_request returned NULL in RMW callback\n");

    response_request->buffer_ptr = original;
 
    response_header.bytes = header->bytes;
    response_header.source_ptr_out = header->source_ptr_out;
    response_header.handle_ptr = header->handle_ptr;

    done_callback.function = A1DI_Request_done;
    done_callback.clientdata = (void *) response_request;

    status = DCMF_Send(&A1D_Rmw_response_protocol,
                       &(response_request->request),
                       done_callback,
                       DCMF_SEQUENTIAL_CONSISTENCY,
                       header->source,
                       header->bytes,
                       original,
                       (DCQuad *) &response_header,
                       (unsigned) 1);
    A1U_ERR_ABORT((status != DCMF_SUCCESS), "DCMF_Send returned with an error A1D_Rmw\n"); 

    A1DI_Release_request(a1d_request);
}

DCMF_Request_t* A1DI_RecvSend_rmw_callback(void *clientdata,
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
    A1D_Buffer_t *a1d_buffer;

    a1d_request = A1DI_Get_request(0);
    A1U_ERR_ABORT(status = (a1d_request == NULL),
                "A1DI_Get_request returned NULL in A1DI_RecvSend_rmw_callback.\n");

    a1d_buffer = A1DI_Get_buffer(sndlen + sizeof(A1D_Rmw_header_t), 0);
    A1U_ERR_ABORT(status = (a1d_buffer == NULL),
                  "A1DI_Get_buffer returned NULL.\n");

    A1DI_Memcpy(a1d_buffer->buffer_ptr, msginfo, sizeof(DCQuad)*count);

    *rcvlen = sndlen;
    *rcvbuf = a1d_buffer->buffer_ptr + sizeof(A1D_Rmw_header_t);

    a1d_request->a1d_buffer_ptr = a1d_buffer;
    cb_done->function = A1DI_RecvDone_rmw_callback;
    cb_done->clientdata = (void *) a1d_request;

    return &(a1d_request->request);
}

void A1DI_RecvSendShort_rmw_callback(void *clientdata,
                                     const DCQuad *msginfo,
                                     unsigned count,
                                     size_t peer,
                                     const char *src,
                                     size_t bytes)
{
    int status = A1_SUCCESS;
    DCMF_Callback_t done_callback;
    A1D_Rmw_header_t *header = (A1D_Rmw_header_t *) msginfo;
    A1_datatype_t datatype = header->datatype;
    void *source = (void *) src;
    void *target = header->target_ptr;
    void *original = NULL;
    A1D_Request_t *response_request = NULL;
    A1D_Rmw_response_header_t response_header;

    status = A1DI_Malloc(&original, bytes);
    A1U_ERR_ABORT(status,
                "A1DI_Malloc returned error in A1DI_RecvSendShort_rmw_callback\n"); 
 
    if(header->op == A1_FETCH_AND_ADD)
    {
       likely_if(datatype == A1_DOUBLE) 
       {
           A1DI_FETCHANDADD_EXECUTE(double, source, target, original, header->bytes/sizeof(double));
       }
       else
       {
          switch (datatype)
          {
          case A1_INT32:
              A1DI_FETCHANDADD_EXECUTE(int32_t, source, target, original, header->bytes/sizeof(int32_t));
              break;
          case A1_INT64:
              A1DI_FETCHANDADD_EXECUTE(int64_t, source, target, original, header->bytes/sizeof(int64_t));
              break; 
          case A1_UINT32:
              A1DI_FETCHANDADD_EXECUTE(uint32_t, source, target, original, header->bytes/sizeof(uint32_t));
              break;
          case A1_UINT64:
              A1DI_FETCHANDADD_EXECUTE(uint64_t, source, target, original, header->bytes/sizeof(uint64_t));
              break;
          case A1_FLOAT:
              A1DI_FETCHANDADD_EXECUTE(float, source, target, original, header->bytes/sizeof(float));
              break;
          default:
              status = A1_ERROR;
              A1U_ERR_ABORT((status != A1_SUCCESS), "Invalid data type in rmw \n");
              break;
          }
       }    
    }
    else if(header->op = A1_SWAP)
    {
        A1DI_Memcpy(original, target, header->bytes); 
        A1DI_Memcpy(target, source, header->bytes); 
    }

    response_request = A1DI_Get_request(0);
    A1U_ERR_ABORT(status = (response_request == NULL), "A1DI_Get_request returned NULL in RMW callback\n");

    response_request->buffer_ptr = original;
 
    response_header.bytes = header->bytes;
    response_header.source_ptr_out = header->source_ptr_out;
    response_header.handle_ptr = header->handle_ptr;

    done_callback.function = A1DI_Request_done;
    done_callback.clientdata = (void *) response_request;

    status = DCMF_Send(&A1D_Rmw_response_protocol,
                       &(response_request->request),
                       done_callback,
                       DCMF_SEQUENTIAL_CONSISTENCY,
                       header->source,
                       header->bytes,
                       original,
                       (DCQuad *) &response_header,
                       (unsigned) 1);
    A1U_ERR_ABORT((status != DCMF_SUCCESS), "DCMF_Send returned with an error A1D_Rmw\n"); 
    
}

int A1DI_Rmw_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Send_Configuration_t conf;
    
    A1U_FUNC_ENTER();
    
    conf.protocol = DCMF_DEFAULT_SEND_PROTOCOL;
    conf.network = DCMF_TORUS_NETWORK;
    conf.cb_recv_short = A1DI_RecvSendShort_rmw_callback;
    conf.cb_recv_short_clientdata = NULL; 
    conf.cb_recv = A1DI_RecvSend_rmw_callback;
    conf.cb_recv_clientdata = NULL;

    status = DCMF_Send_register(&A1D_Rmw_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "rmw registartion returned with error %d \n",
                status);

    conf.cb_recv_short = A1DI_RecvSendShort_rmw_response_callback;
    conf.cb_recv = A1DI_RecvSend_rmw_response_callback;

    status = DCMF_Send_register(&A1D_Rmw_response_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "rmw_response registartion returned with error %d \n",
                status);

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1D_Rmw(int target,
            void* source_ptr_in,
            void* source_ptr_out,
            void* target_ptr,
            int bytes,
            A1_atomic_op_t op,
            A1_datatype_t a1_type)
{
    int status = A1_SUCCESS;
    DCMF_Callback_t done_callback;
    DCMF_Request_t request;
    A1D_Rmw_header_t header;
    A1D_Handle_t *a1d_handle;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_handle = A1DI_Get_handle();
    A1U_ERR_POP(status = (a1d_handle == NULL),
                "A1DI_Get_handle returned NULL in A1D_Rmw\n");

    header.source = A1D_Process_info.my_rank; 
    header.source_ptr_out = source_ptr_out;
    header.target_ptr = target_ptr;
    header.bytes = bytes;
    header.op = op;
    header.datatype = a1_type;
    header.handle_ptr = a1d_handle;
 
    a1d_handle->active++;

    status = DCMF_Send(&A1D_Rmw_protocol,
                       &request,
                       A1D_Nocallback,
                       DCMF_SEQUENTIAL_CONSISTENCY,
                       target,
                       bytes,
                       source_ptr_in,
                       (DCQuad *) &header,
                       (unsigned) 2);
    A1U_ERR_POP((status != DCMF_SUCCESS), "DCMF_Send returned with an error A1D_Rmw\n");    

    A1DI_Conditional_advance(a1d_handle->active > 0);

  fn_exit:
    A1DI_Release_handle(a1d_handle);
    A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
