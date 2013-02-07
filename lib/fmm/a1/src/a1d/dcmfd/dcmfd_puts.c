/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

DCMF_Protocol_t A1D_Packed_puts_protocol;

void A1DI_RecvDone_packedputs_callback(void *clientdata, DCMF_Error_t *error)
{
    A1D_Request_t *a1d_request = (A1D_Request_t *) clientdata;
    A1D_Buffer_t *a1d_buffer = a1d_request->a1d_buffer_ptr;
    A1D_Packed_puts_header_t *header =
            (A1D_Packed_puts_header_t *) a1d_buffer->buffer_ptr;
    int complete = 0;

    A1DI_Unpack_strided((void *) ((size_t) a1d_buffer->buffer_ptr + sizeof(A1D_Packed_puts_header_t)),
                        header->data_size,
                        header->stride_level,
                        header->block_sizes,
                        header->target_ptr,
                        header->trg_stride_ar,
                        header->block_idx,
                        &complete);

    A1DI_Release_request(a1d_request);
}

DCMF_Request_t* A1DI_RecvSend_packedputs_callback(void *clientdata,
                                                  const DCQuad *msginfo,
                                                  unsigned count, /* TODO: this is not used */
                                                  size_t peer,
                                                  size_t sndlen,
                                                  size_t *rcvlen,
                                                  char **rcvbuf,
                                                  DCMF_Callback_t *cb_done)
{
    int status = 0;
    A1D_Request_t *a1d_request = NULL;
    A1D_Buffer_t *a1d_buffer = NULL;

    a1d_request = A1DI_Get_request(0);
    A1U_ERR_ABORT(status = (a1d_request == NULL),
                  "A1DI_Get_request returned NULL in A1DI_RecvSend_packedputs_callback.\n");

    a1d_buffer = A1DI_Get_buffer(a1d_settings.put_packetsize, 0);
    A1U_ERR_ABORT(status = (a1d_buffer == NULL),
                  "A1DI_Get_buffer returned NULL.\n");

    *rcvlen = sndlen;
    *rcvbuf = a1d_buffer->buffer_ptr;

    a1d_request->a1d_buffer_ptr = a1d_buffer;
    cb_done->function = A1DI_RecvDone_packedputs_callback;
    cb_done->clientdata = (void *) a1d_request;

    return &(a1d_request->request);
}

void A1DI_RecvSendShort_packedputs_callback(void *clientdata,
                                            const DCQuad *msginfo,
                                            unsigned count, /* TODO: this is not used */
                                            size_t peer,
                                            const char *src,
                                            size_t bytes)
{
    A1D_Packed_puts_header_t *header;
    void *packet_ptr = (void *) src;
    int complete = 0;

    A1U_FUNC_ENTER();

    header = (A1D_Packed_puts_header_t *) packet_ptr;

    A1DI_Unpack_strided((void *) ((size_t) packet_ptr + sizeof(A1D_Packed_puts_header_t)),
                        header->data_size,
                        header->stride_level,
                        header->block_sizes,
                        header->target_ptr,
                        header->trg_stride_ar,
                        header->block_idx,
                        &complete);
}

int A1DI_Packed_puts_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Send_Configuration_t conf;

    A1U_FUNC_ENTER();

    conf.protocol = DCMF_DEFAULT_SEND_PROTOCOL;
    conf.network = DCMF_TORUS_NETWORK;
    conf.cb_recv_short = A1DI_RecvSendShort_packedputs_callback;
    conf.cb_recv_short_clientdata = NULL;
    conf.cb_recv = A1DI_RecvSend_packedputs_callback;
    conf.cb_recv_clientdata = NULL;

    status = DCMF_Send_register(&A1D_Packed_puts_protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "DCMF_Send_register returned with error %d \n",
                status);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1DI_Packed_puts(int target,
                     int stride_level,
                     int *block_sizes,
                     void* source_ptr,
                     int *src_stride_ar,
                     void* target_ptr,
                     int *trg_stride_ar)
{
    int status = A1_SUCCESS;
    DCMF_Callback_t done_callback;
    A1D_Request_t *a1d_request = NULL;
    A1D_Buffer_t *a1d_buffer = NULL;
    A1D_Packed_puts_header_t header;
    void *packet_ptr = NULL, *data_ptr = NULL;
    int packet_size = 0, data_size = 0 , data_limit = 0;
    int block_idx[A1C_MAX_STRIDED_DIM];
    int complete = 0;

    A1U_FUNC_ENTER();

    A1DI_Memset(block_idx, 0, (stride_level + 1) * sizeof(int));

    header.stride_level = stride_level;
    A1DI_Memcpy(header.trg_stride_ar, trg_stride_ar, stride_level * sizeof(int));
    A1DI_Memcpy(header.block_sizes, block_sizes, (stride_level + 1) * sizeof(int));

    while (!complete)
    {
        header.target_ptr = target_ptr;
        A1DI_Memcpy(header.block_idx, block_idx, (stride_level + 1) * sizeof(int));

        /*Fetching buffer from the pool*/
        a1d_buffer = A1DI_Get_buffer(a1d_settings.put_packetsize, 1);
        A1U_ERR_ABORT(status = (a1d_buffer == NULL),
                  "A1DI_Get_buffer returned NULL.\n");

        packet_ptr = a1d_buffer->buffer_ptr;

        data_ptr = (void *) ((size_t) packet_ptr
              + sizeof(A1D_Packed_puts_header_t));
        data_limit = a1d_settings.put_packetsize
              - sizeof(A1D_Packed_puts_header_t);
        data_size = 0;

        /*The packing function can modify the source ptr, target ptr, and block index*/
        A1DI_Pack_strided(data_ptr,
                          data_limit,
                          stride_level,
                          block_sizes,
                          &source_ptr,
                          src_stride_ar,
                          &target_ptr,
                          trg_stride_ar,
                          block_idx,
                          &data_size,
                          &complete);

        /*Setting data size information in the header and copying it into the packet*/
        header.data_size = data_size;
        A1DI_Memcpy((void *) packet_ptr, (void *) &header, sizeof(A1D_Packed_puts_header_t));

        packet_size = data_size + sizeof(A1D_Packed_puts_header_t);

        /*Fetching request from the pool*/
        a1d_request = A1DI_Get_request(1);
        A1U_ERR_POP(status = (a1d_request == NULL),
                    "A1DI_Get_request returned error.\n");
        a1d_request->a1d_buffer_ptr = a1d_buffer;

        done_callback.function = A1DI_Request_done;
        done_callback.clientdata = (void *) a1d_request;

        status = DCMF_Send(&A1D_Packed_puts_protocol,
                           &(a1d_request->request),
                           done_callback,
                           DCMF_SEQUENTIAL_CONSISTENCY,
                           target,
                           packet_size,
                           packet_ptr,
                           NULL,
                           0);
        A1U_ERR_POP(status != DCMF_SUCCESS,
                    "DCMF_Send returned with an error \n");

        A1D_Connection_send_active[target]++;
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1DI_Direct_puts(int target,
                     int stride_level,
                     int *block_sizes,
                     void *source_ptr,
                     int *src_stride_ar,
                     void *target_ptr,
                     int *trg_stride_ar,
                     A1D_Handle_t *a1d_handle)
{
    int i, j, status = A1_SUCCESS;
    size_t src_disp, dst_disp;
    DCMF_Callback_t done_callback;
    A1D_Request_t *a1d_request;
    int chunk_count = 1;
    int *block_sizes_w;
    int y = 0;

    A1U_FUNC_ENTER();

    status = A1DI_Malloc((void **) &block_sizes_w, sizeof(int) * (stride_level + 1));
    A1U_ERR_POP(status != 0,
                "A1DI_Malloc returned error in A1DI_Direct_puts");

    A1DI_Memcpy(block_sizes_w, block_sizes, sizeof(int) * (stride_level + 1));

    for (i = 1; i <= stride_level; i++)
        chunk_count = block_sizes[i] * chunk_count;

    for (i = 0; i < chunk_count; i++)
    {
        src_disp = (size_t) source_ptr - (size_t) A1D_Membase_global[A1D_Process_info.my_rank];
        dst_disp = (size_t) target_ptr - (size_t) A1D_Membase_global[target];

        a1d_request = A1DI_Get_request(1);
        A1U_ERR_POP(status = (a1d_request == NULL),
                    "A1DI_Get_request returned error.\n");
        A1DI_Set_handle(a1d_request, a1d_handle);

        done_callback.function = A1DI_Request_done;
        done_callback.clientdata = (void *) a1d_request;

        a1d_handle->active++;

        status = DCMF_Put(&A1D_Generic_put_protocol,
                          &(a1d_request->request),
                          done_callback,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          target,
                          block_sizes[0],
                          &A1D_Memregion_global[A1D_Process_info.my_rank],
                          &A1D_Memregion_global[target],
                          src_disp,
                          dst_disp,
                          A1D_Nocallback);
        A1U_ERR_POP(status != DCMF_SUCCESS,
                    "DCMF_Put returned with an error \n");

        A1D_Connection_put_active[target]++;

        block_sizes_w[1]--;
        if (block_sizes_w[1] == 0)
        {
            y = 1;
            while (block_sizes_w[y] == 0)
            {
                if (y == stride_level)
                {
                    A1U_ASSERT(i == chunk_count - 1, status);
                    return status;
                }
                y++;
            }
            block_sizes_w[y]--;

            /*The strides done on lower dimension should be subtracted as these are
              included in the stride along the current dimension*/
            source_ptr = (void *) ((size_t) source_ptr 
                   + src_stride_ar[y - 1] 
                   - (block_sizes[y-1] - 1) * src_stride_ar[y-2]);
            target_ptr = (void *) ((size_t) target_ptr 
                   + trg_stride_ar[y - 1] 
                   - (block_sizes[y-1] - 1) * trg_stride_ar[y-2]);

            y--;
            while (y >= 1)
            {
                block_sizes_w[y] = block_sizes[y];
                y--;
            }
        }
        else
        {
            source_ptr = (void *) ((size_t) source_ptr + src_stride_ar[0]);
            target_ptr = (void *) ((size_t) target_ptr + trg_stride_ar[0]);
        }
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1DI_Recursive_puts(int target,
                        int stride_level,
                        int *block_sizes,
                        void *source_ptr,
                        int *src_stride_ar,
                        void *target_ptr,
                        int *trg_stride_ar,
                        A1D_Handle_t *a1d_handle)
{
    int i, status = A1_SUCCESS;
    size_t src_disp, dst_disp;
    DCMF_Callback_t done_callback;
    A1D_Request_t *a1d_request;

    A1U_FUNC_ENTER();

    if (stride_level > 0)
    {
        for (i = 0; i < block_sizes[stride_level]; i++)
        {
            status = A1DI_Recursive_puts(target,
                                         stride_level - 1,
                                         block_sizes,
                                         (void *) ((size_t) source_ptr + i * src_stride_ar[stride_level - 1]),
                                         src_stride_ar,
                                         (void *) ((size_t) target_ptr + i * trg_stride_ar[stride_level - 1]),
                                         trg_stride_ar,
                                         a1d_handle);
            A1U_ERR_POP(status != A1_SUCCESS,
                        "A1DI_Recursive_puts returned error in A1DI_Recursive_puts.\n");
        }
    }
    else
    {
        src_disp = (size_t) source_ptr - (size_t) A1D_Membase_global[A1D_Process_info.my_rank];
        dst_disp = (size_t) target_ptr - (size_t) A1D_Membase_global[target];

        a1d_request = A1DI_Get_request(1);
        A1U_ERR_POP(status = (a1d_request == NULL),
                    "A1DI_Get_request returned error.\n");
        A1DI_Set_handle(a1d_request, a1d_handle);

        done_callback.function = A1DI_Request_done;
        done_callback.clientdata = (void *) a1d_request;

        a1d_handle->active++;

        status = DCMF_Put(&A1D_Generic_put_protocol,
                          &(a1d_request->request),
                          done_callback,
                          DCMF_SEQUENTIAL_CONSISTENCY,
                          target,
                          block_sizes[0],
                          &A1D_Memregion_global[A1D_Process_info.my_rank],
                          &A1D_Memregion_global[target],
                          src_disp,
                          dst_disp,
                          A1D_Nocallback);
        A1U_ERR_POP(status != DCMF_SUCCESS,
                    "DCMF_Put returned with an error \n");

        A1D_Connection_put_active[target]++;
    }

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_PutS(int target,
             int stride_level,
             int *block_sizes,
             void* source_ptr,
             int *src_stride_ar,
             void* target_ptr,
             int *trg_stride_ar)
{
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle = NULL;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    if (block_sizes[0] > a1d_settings.put_packing_limit)
    {
        a1d_handle = A1DI_Get_handle();
        A1U_ERR_POP(status = (a1d_handle == NULL),
                    "A1DI_Get_handle returned NULL in A1D_PutS\n");

        status = A1DI_Direct_puts(target,
                                  stride_level,
                                  block_sizes,
                                  source_ptr,
                                  src_stride_ar,
                                  target_ptr,
                                  trg_stride_ar,
                                  a1d_handle);
        A1U_ERR_POP(status, "A1DI_Direct_puts returned with an error \n");

        A1DI_Conditional_advance(a1d_handle->active > 0);

        A1DI_Release_handle(a1d_handle);
    }
    else
    {
        status = A1DI_Packed_puts(target,
                                  stride_level,
                                  block_sizes,
                                  source_ptr,
                                  src_stride_ar,
                                  target_ptr,
                                  trg_stride_ar);
        A1U_ERR_POP(status, "A1DI_Packed_puts returned with an error \n");
    }

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

int A1D_NbPutS(int target,
               int stride_level,
               int *block_sizes,
               void* source_ptr,
               int *src_stride_ar,
               void* target_ptr,
               int *trg_stride_ar,
               A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    A1D_Handle_t *a1d_handle = NULL;

    A1U_FUNC_ENTER();

    A1DI_CRITICAL_ENTER();

    a1d_handle = (A1D_Handle_t *) a1_handle;

    if (block_sizes[0] > a1d_settings.put_packing_limit)
    {

        status = A1DI_Direct_puts(target,
                                  stride_level,
                                  block_sizes,
                                  source_ptr,
                                  src_stride_ar,
                                  target_ptr,
                                  trg_stride_ar,
                                  a1d_handle);
        A1U_ERR_POP(status, "A1DI_Direct_puts returned with an error \n");

    }
    else
    {

        if (a1d_settings.use_handoff)
        {
            A1D_Op_handoff *op_handoff;
            status = A1DI_Malloc((void **) &op_handoff, sizeof(A1D_Op_handoff));
            A1U_ERR_POP(status != 0,
                       "A1DI_Malloc returned with error \n");

            op_handoff->op_type = A1D_PACKED_PUTS;
            op_handoff->op.puts_op.target = target;
            op_handoff->op.puts_op.stride_level = stride_level;
            op_handoff->op.puts_op.block_sizes = block_sizes;
            op_handoff->op.puts_op.source_ptr = source_ptr;
            op_handoff->op.puts_op.src_stride_ar = src_stride_ar;
            op_handoff->op.puts_op.target_ptr = target_ptr;
            op_handoff->op.puts_op.trg_stride_ar = trg_stride_ar;
            op_handoff->op.puts_op.a1d_handle = a1d_handle;

            a1d_handle->active++;

            if (A1D_Op_handoff_queuetail == NULL)
            {
                A1D_Op_handoff_queuehead = op_handoff;
                A1D_Op_handoff_queuetail = op_handoff;
            }
            else
            {
                A1D_Op_handoff_queuetail->next = op_handoff;
                A1D_Op_handoff_queuetail = op_handoff;
            }
            op_handoff->next = NULL;

        }
        else
        {
            status = A1DI_Packed_puts(target,
                                      stride_level,
                                      block_sizes,
                                      source_ptr,
                                      src_stride_ar,
                                      target_ptr,
                                      trg_stride_ar);
            A1U_ERR_POP(status, "A1DI_Packed_puts returned with an error \n");
        }

    }

    fn_exit: A1DI_CRITICAL_EXIT();
    A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}
