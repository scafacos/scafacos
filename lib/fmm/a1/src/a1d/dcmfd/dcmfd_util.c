/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

A1D_Control_xchange_info_t A1D_Control_xchange_info;

/*************************************************************
 Local Completion Callbacks
 **************************************************************/

void A1DI_Generic_done(void *clientdata, DCMF_Error_t *error)
{
    --(*((uint32_t *) clientdata));
}

void A1DI_Request_done(void *clientdata, DCMF_Error_t *error)
{
    A1DI_Release_request((A1D_Request_t *) clientdata);
}

/*************************************************************
 Control Protocol for Information Exchange
 **************************************************************/

void A1DI_Control_xchange_callback(void *clientdata,
                                   const DCMF_Control_t *info,
                                   size_t peer)
{
    A1DI_Memcpy((void *) ((size_t) A1D_Control_xchange_info.xchange_ptr
                   + (size_t)(peer * A1D_Control_xchange_info.xchange_size)),
                (void *) info,
                A1D_Control_xchange_info.xchange_size);

    --(*((uint32_t *) clientdata));
}

int A1DI_Control_xchange_initialize()
{
    int status = A1_SUCCESS;
    DCMF_Control_Configuration_t conf;

    A1U_FUNC_ENTER();

    conf.protocol = DCMF_DEFAULT_CONTROL_PROTOCOL;
    conf.network = DCMF_DEFAULT_NETWORK;
    conf.cb_recv = A1DI_Control_xchange_callback;
    conf.cb_recv_clientdata = (void *) &A1D_Control_xchange_info.rcv_active;

    status = DCMF_Control_register(&A1D_Control_xchange_info.protocol, &conf);
    A1U_ERR_POP(status != DCMF_SUCCESS,
                "Control xchange registartion returned with error %d \n",
                status);

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail: goto fn_exit;
}

/*************************************************************
 Data Packing Code
 **************************************************************/

int A1DI_Pack_strided(void *data_ptr,
                      int data_limit,
                      int stride_level,
                      int *block_sizes,
                      void **source_ptr,
                      int *src_stride_ar,
                      void **target_ptr,
                      int *trg_stride_ar,
                      int *block_idx,
                      int *data_size,
                      int *complete)
{
    int status = A1_SUCCESS;
    int y, index, size_data;
    int block_sizes_w[A1C_MAX_STRIDED_DIM];

    A1U_FUNC_ENTER();

    *complete = 0;
    *data_size = 0;

    while ((*data_size + block_sizes[0]) <= data_limit)
    {

        A1DI_Memcpy(data_ptr, *source_ptr, block_sizes[0]);
        data_ptr = (void *) ((size_t) data_ptr + block_sizes[0]);
        *data_size = *data_size + block_sizes[0];

        block_idx[1]++;
        if (block_idx[1] == block_sizes[1])
        {
            y = 1;
            while (block_idx[y] == block_sizes[y])
            {
                if (y == stride_level)
                {
                    *complete = 1;
                    return status;
                }
                y++;
            }
            block_idx[y]++;

            /*The strides done on lower dimension should be subtracted as these are
              included in the stride along the current dimension*/ 
            *source_ptr = (void *) ((size_t) *source_ptr
                    + src_stride_ar[y - 1] 
                    - (block_sizes[y-1] - 1) * src_stride_ar[y-2]);
            *target_ptr = (void *) ((size_t) *target_ptr
                    + trg_stride_ar[y - 1] 
                    - (block_sizes[y-1] - 1) * trg_stride_ar[y-2]);

            y--;
            while (y >= 1)
            {
                block_idx[y] = 0;
                y--;
            }
        }
        else
        {
            *source_ptr = (void *) ((size_t) *source_ptr + src_stride_ar[0]);
            *target_ptr = (void *) ((size_t) *target_ptr + trg_stride_ar[0]);
        }
    }

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1DI_Unpack_strided(void *data_ptr,
                        int data_size,
                        int stride_level,
                        int *block_sizes,
                        void *target_ptr,
                        int *trg_stride_ar,
                        int *block_idx,
                        int *complete)
{
    int status = A1_SUCCESS;
    int y, index;

    A1U_FUNC_ENTER();

    while (data_size > 0)
    {
        A1DI_Memcpy(target_ptr, data_ptr, block_sizes[0]);

        data_ptr = (void *) ((size_t) data_ptr + block_sizes[0]);
        data_size = data_size - block_sizes[0];

        block_idx[1]++;
        if (block_idx[1] == block_sizes[1])
        {
            y = 1;
            while (block_idx[y] == block_sizes[y])
            {
                if (y == stride_level)
                {
                    *complete = 1;
                    break;
                }
                y++;
            }
            block_idx[y]++;

            /*The strides done on lower dimension should be subtracted as these are
              included in the stride along the current dimension*/
            target_ptr = (void *) ((size_t) target_ptr 
                   + trg_stride_ar[y - 1] 
                   - (block_sizes[y-1] - 1) * trg_stride_ar[y-2]);

            y--;
            while (y >= 1)
            {
                block_idx[y] = 0;
                y--;
            }
        }
        else
        {
            target_ptr = (void *) ((size_t) target_ptr + trg_stride_ar[0]);
        }
    }

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1DI_Unpack_strided_acc(void *data_ptr,
                            int data_size,
                            int stride_level,
                            int *block_sizes,
                            void *target_ptr,
                            int *trg_stride_ar,
                            int *block_idx,
                            A1_datatype_t a1_type,
                            void *scaling,
                            int *complete)
{
    int status = A1_SUCCESS;
    int y;

    A1U_FUNC_ENTER();

    while (data_size > 0)
    {
        likely_if(a1_type == A1_DOUBLE)
        {
            A1DI_ACC(double,
                    data_ptr,
                    target_ptr,
                    *((double *) scaling),
                    block_sizes[0]/sizeof(double));
        }
        else
        {
            switch (a1_type)
            {
                case A1_INT32:
                A1DI_ACC(int32_t,
                        data_ptr,
                        target_ptr,
                        *((int32_t *) scaling),
                        block_sizes[0] / sizeof(int32_t));
                break;
                case A1_INT64:
                A1DI_ACC(int64_t,
                        data_ptr,
                        target_ptr,
                        *((int64_t *) scaling),
                        block_sizes[0] / sizeof(int64_t));
                break;
                case A1_UINT32:
                A1DI_ACC(uint32_t,
                        data_ptr,
                        target_ptr,
                        *((uint32_t *) scaling),
                        block_sizes[0] / sizeof(uint32_t));
                break;
                case A1_UINT64:
                A1DI_ACC(uint64_t,
                        data_ptr,
                        target_ptr,
                        *((uint64_t *) scaling),
                        block_sizes[0] / sizeof(uint64_t));
                break;
                case A1_FLOAT:
                A1DI_ACC(float,
                        data_ptr,
                        target_ptr,
                        *((float *) scaling),
                        block_sizes[0]/sizeof(float));
                break;
                default:
                A1U_ERR_ABORT(A1_ERROR,
                        "Invalid datatype received in Putacc operation \n");
                break;
            }
        }

        data_ptr = (void *) ((size_t) data_ptr + block_sizes[0]);
        data_size = data_size - block_sizes[0];

        block_idx[1]++;
        if (block_idx[1] == block_sizes[1])
        {
            y = 1;
            while (block_idx[y] == block_sizes[y])
            {
                if (y == stride_level)
                {
                    *complete = 1;
                    break;
                }
                y++;
            }
            block_idx[y]++;

            /*The strides done on lower dimension should be subtracted as these are
              included in the stride along the current dimension*/
            target_ptr = (void *) ((size_t) target_ptr 
                   + trg_stride_ar[y - 1] 
                   - (block_sizes[y-1] - 1) * trg_stride_ar[y-2]);

            y--;
            while (y >= 1)
            {
                block_idx[y] = 0;
                y--;
            }
        }
        else
        {
            target_ptr = (void *) ((size_t) target_ptr + trg_stride_ar[0]);
        }
    }

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

/**** Exposing locking to ADL layer ****/

void A1D_Global_lock_acquire()
{
    A1U_FUNC_ENTER();

    A1DI_GLOBAL_LOCK_ACQUIRE();

  fn_exit:
    A1U_FUNC_EXIT();
    return;

  fn_fail:
    goto fn_exit;
}

void A1D_Global_lock_release()
{
    A1U_FUNC_ENTER();

    A1DI_GLOBAL_LOCK_RELEASE();

  fn_exit:
    A1U_FUNC_EXIT();
    return;

  fn_fail:
    goto fn_exit;
}
