/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1u.h"
#include "a1d.h"

#define A1UI_ACC(datatype, source, target, scaling, count)                  \
        do {                                                                     \
            int w;                                                                 \
            datatype *s = (datatype *) source;                                     \
            datatype *t = (datatype *) target;                                     \
            datatype c = (datatype) scaling;                                       \
            for(w=0; w<count; w++)                                                 \
            t[w] += s[w]*c;                                                   \
        } while(0)                                                               \

int A1U_Acc_memcpy(void* source_ptr,
                   void* target_ptr,
                   int bytes,
                   A1_datatype_t a1_type,
                   void* scaling)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1D_Global_lock_acquire();

    switch (a1_type)
    {
        case A1_DOUBLE:
            A1UI_ACC(double,
                     source_ptr,
                     target_ptr,
                     *((double *) scaling),
                     bytes/sizeof(double));
        case A1_INT32:
            A1UI_ACC(int32_t,
                     source_ptr,
                     target_ptr,
                     *((int32_t *) scaling),
                     bytes/sizeof(int32_t));
            break;
        case A1_INT64:
            A1UI_ACC(int64_t,
                     source_ptr,
                     target_ptr,
                     *((int64_t *) scaling),
                     bytes/sizeof(int64_t));
            break;
        case A1_UINT32:
            A1UI_ACC(uint32_t,
                     source_ptr,
                     target_ptr,
                     *((uint32_t *) scaling),
                     bytes/sizeof(uint32_t));
            break;
        case A1_UINT64:
            A1UI_ACC(uint64_t,
                     source_ptr,
                     target_ptr,
                     *((uint64_t *) scaling),
                     bytes/sizeof(uint64_t));
            break;
        case A1_FLOAT:
            A1UI_ACC(float,
                     source_ptr,
                     target_ptr,
                     *((float *) scaling),
                     bytes/sizeof(float));
            break;
        default:
            status = A1_ERROR;
            A1U_ERR_POP((status != A1_SUCCESS), "Invalid datatype in A1U_Acc_memcpy \n");
            break;
    }

    A1D_Global_lock_release();

    fn_exit:
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

int A1U_AccS_memcpy(int stride_level,
                    int *block_sizes,
                    void* source_ptr,
                    int *src_stride_ar,
                    void* target_ptr,
                    int *trg_stride_ar,
                    A1_datatype_t a1_type,
                    void* scaling)
{
    int status = A1_SUCCESS;
    int chunk_count = 1;
    int *block_sizes_w;
    int i, y;

    A1U_FUNC_ENTER();

    A1D_Global_lock_acquire();

    block_sizes_w = malloc(sizeof(int) * (stride_level + 1));
    A1U_ERR_POP((status = (NULL == block_sizes_w)),
                "malloc failed in A1U_PutS_memcpy");

    memcpy(block_sizes_w, block_sizes, sizeof(int) * (stride_level + 1));

    for (i = 1; i <= stride_level; i++)
        chunk_count = block_sizes[i] * chunk_count;

    for (i = 0; i < chunk_count; i++)
    {
        switch (a1_type)
        {
            case A1_DOUBLE:
                A1UI_ACC(double,
                         source_ptr,
                         target_ptr,
                         *((double *) scaling),
                         block_sizes[0]/sizeof(double));
                break;
            case A1_INT32:
                A1UI_ACC(int32_t,
                         source_ptr,
                         target_ptr,
                         *((int32_t *) scaling),
                         block_sizes[0]/sizeof(int32_t));
                break;
            case A1_INT64:
                A1UI_ACC(int64_t,
                         source_ptr,
                         target_ptr,
                         *((int64_t *) scaling),
                         block_sizes[0]/sizeof(int64_t));
                break;
            case A1_UINT32:
                A1UI_ACC(uint32_t,
                         source_ptr,
                         target_ptr,
                         *((uint32_t *) scaling),
                         block_sizes[0]/sizeof(uint32_t));
                break;
            case A1_UINT64:
                A1UI_ACC(uint64_t,
                         source_ptr,
                         target_ptr,
                         *((uint64_t *) scaling),
                         block_sizes[0]/sizeof(uint64_t));
                break;
            case A1_FLOAT:
                A1UI_ACC(float,
                         source_ptr,
                         target_ptr,
                         *((float *) scaling),
                         block_sizes[0]/sizeof(float));
                break;
            default:
                status = A1_ERROR;
                A1U_ERR_POP((status != A1_SUCCESS), "Invalid data type in putacc \n");
                break;
        }

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

            /*The strides done on lower dimensions should be subtracted as these are
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

    A1D_Global_lock_release();

    fn_exit:
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

int A1U_AccV_memcpy(A1_iov_t *iov_ar,
                    int ar_len,
                    A1_datatype_t a1_type,
                    void* scaling)
{
    int i, j, status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    A1D_Global_lock_acquire();

    for (i=0; i<ar_len; i++)
    {
        for(j=0; j<iov_ar[i].ptr_ar_len; j++) 
        {
            switch (a1_type)
            {
                case A1_DOUBLE:
                    A1UI_ACC(double,
                             iov_ar[i].source_ptr_ar[j],
                             iov_ar[i].target_ptr_ar[j],
                             *((double *) scaling),
                             (iov_ar[i].size)/sizeof(double));
                    break;
                case A1_INT32:
                    A1UI_ACC(int32_t,
                             iov_ar[i].source_ptr_ar[j],
                             iov_ar[i].target_ptr_ar[j],
                             *((int32_t *) scaling),
                             (iov_ar[i].size)/sizeof(int32_t));
                    break;
                case A1_INT64:
                    A1UI_ACC(int64_t,
                             iov_ar[i].source_ptr_ar[j],
                             iov_ar[i].target_ptr_ar[j],
                             *((int64_t *) scaling),
                             (iov_ar[i].size)/sizeof(int64_t));
                    break;
                case A1_UINT32:
                    A1UI_ACC(uint32_t,
                             iov_ar[i].source_ptr_ar[j],
                             iov_ar[i].target_ptr_ar[j],
                             *((uint32_t *) scaling),
                             (iov_ar[i].size)/sizeof(uint32_t));
                    break;
                case A1_UINT64:
                    A1UI_ACC(uint64_t,
                             iov_ar[i].source_ptr_ar[j],
                             iov_ar[i].target_ptr_ar[j],
                             *((uint64_t *) scaling),
                             (iov_ar[i].size)/sizeof(uint64_t));
                    break;
                case A1_FLOAT:
                    A1UI_ACC(float,
                             iov_ar[i].source_ptr_ar[j],
                             iov_ar[i].target_ptr_ar[j],
                             *((float *) scaling),
                             (iov_ar[i].size)/sizeof(float));
                    break;
                default:
                    status = A1_ERROR;
                    A1U_ERR_POP((status != A1_SUCCESS), "Invalid data type in putacc \n");
                    break;
            }
        }
    }

    A1D_Global_lock_release();

    fn_exit: A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}
