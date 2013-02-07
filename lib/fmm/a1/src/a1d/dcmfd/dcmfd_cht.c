/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "dcmfdimpl.h"

_BGP_Atomic global_atomic;
LockBox_Mutex_t global_lbmutex;

pthread_t A1DI_CHT_pthread;

A1D_Op_handoff *A1D_Op_handoff_queuehead = NULL;
A1D_Op_handoff *A1D_Op_handoff_queuetail = NULL;

volatile int A1D_Inside_handoff;

void A1DI_Handoff_progress()
{
    int status = A1_SUCCESS;
    A1D_Op_handoff *op_handoff;

    A1U_FUNC_ENTER();

    A1D_Inside_handoff = 1;

    if (A1D_Op_handoff_queuehead)
    {
        op_handoff = A1D_Op_handoff_queuehead;

        if (op_handoff->op_type == A1D_PACKED_PUTS)
        {

            status = A1DI_Packed_puts(op_handoff->op.puts_op.target,
                                      op_handoff->op.puts_op.stride_level,
                                      op_handoff->op.puts_op.block_sizes,
                                      op_handoff->op.puts_op.source_ptr,
                                      op_handoff->op.puts_op.src_stride_ar,
                                      op_handoff->op.puts_op.target_ptr,
                                      op_handoff->op.puts_op.trg_stride_ar);
            A1U_ERR_ABORT(status,
                          "A1DI_Packed_puts returned with an error handoff progress\n");

            op_handoff->op.puts_op.a1d_handle->active--;

        }
        else if (op_handoff->op_type == A1D_PACKED_PUTACCS)
        {

            status = A1DI_Packed_putaccs(op_handoff->op.putaccs_op.target,
                                         op_handoff->op.putaccs_op.stride_level,
                                         op_handoff->op.putaccs_op.block_sizes,
                                         op_handoff->op.putaccs_op.source_ptr,
                                         op_handoff->op.putaccs_op.src_stride_ar,
                                         op_handoff->op.putaccs_op.target_ptr,
                                         op_handoff->op.putaccs_op.trg_stride_ar,
                                         op_handoff->op.putaccs_op.datatype,
                                         op_handoff->op.putaccs_op.scaling);
            A1U_ERR_ABORT(status,
                          "A1DI_Packed_putaccs returned with an error handoff progress\n");

            op_handoff->op.putaccs_op.a1d_handle->active--;

        }
        else
        {
            A1U_ERR_ABORT(status = A1_ERROR,
                          "Invalid op encountered in handoff progress. \n");
        }

        if (A1D_Op_handoff_queuehead == A1D_Op_handoff_queuetail)
        {
            A1D_Op_handoff_queuehead = NULL;
            A1D_Op_handoff_queuetail = NULL;
        }
        else
        {
            A1D_Op_handoff_queuehead = A1D_Op_handoff_queuehead->next;
        }

        A1DI_Free(op_handoff);
    }

    A1D_Inside_handoff = 0;

    fn_exit: A1U_FUNC_EXIT();
    return;

    fn_fail: goto fn_exit;
}

void *A1DI_CHT_advance_lock(void * dummy)
{
    A1DI_GLOBAL_LOCK_ACQUIRE();
    while (1)
    {
        DCMF_Messager_advance(0);
        if (a1d_settings.use_handoff && (A1D_Inside_handoff==0))
        {
            A1DI_Handoff_progress();
        }
        A1DI_GLOBAL_LOCK_RELEASE();
        A1DI_Wait_cycles(a1d_settings.cht_pause_cycles);
        A1DI_GLOBAL_LOCK_ACQUIRE();
    }
    A1DI_GLOBAL_LOCK_RELEASE();
}

