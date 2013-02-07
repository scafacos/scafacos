/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1u.h"
#include "a1d.h"

#define A1UI_MOD_BXOR(datatype, source, target, count)                             \
        do {                                                                       \
            int w;                                                                 \
            datatype *s = (datatype *) source;                                     \
            datatype *t = (datatype *) target;                                     \
            for(w=0; w<count; w++)                                                 \
            t[w] ^= s[w];                                                          \
        } while(0)                                                                 \

int A1U_ModV_memcpy(A1_iov_t *iov_ar,
                    int ar_len,
                    A1_reduce_op_t a1_op,
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
            switch (a1_op)
            {
                case A1_BXOR:
                    switch (a1_type)
                    {
                        case A1_INT32:
                            A1UI_MOD_BXOR(int32_t,
                                     iov_ar[i].source_ptr_ar[j],
                                     iov_ar[i].target_ptr_ar[j],
                                     (iov_ar[i].size)/sizeof(int32_t));
                            break;
                        case A1_INT64:
                            A1UI_MOD_BXOR(int64_t,
                                     iov_ar[i].source_ptr_ar[j],
                                     iov_ar[i].target_ptr_ar[j],
                                     (iov_ar[i].size)/sizeof(int64_t));
                            break;
                        case A1_UINT32:
                            A1UI_MOD_BXOR(uint32_t,
                                     iov_ar[i].source_ptr_ar[j],
                                     iov_ar[i].target_ptr_ar[j],
                                     (iov_ar[i].size)/sizeof(uint32_t));
                            break;
                        case A1_UINT64:
                            A1UI_MOD_BXOR(uint64_t,
                                     iov_ar[i].source_ptr_ar[j],
                                     iov_ar[i].target_ptr_ar[j],
                                     (iov_ar[i].size)/sizeof(uint64_t));
                            break;
                        default:
                            status = A1_ERROR;
                            A1U_ERR_POP((status != A1_SUCCESS), "Invalid data type in A1U_AccV_memcpy\n");
                            break;
                    }
                    break;
                default:
                    status = A1_ERROR;
                    A1U_ERR_POP((status != A1_SUCCESS), "Invalid op type in A1U_AccV_memcpy\n");
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
