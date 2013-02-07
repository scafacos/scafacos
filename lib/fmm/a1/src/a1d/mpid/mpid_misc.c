/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpidimpl.h"

int A1D_Rank()
{
    A1U_FUNC_ENTER();

  fn_exit:
    A1U_FUNC_EXIT();
    return A1D_Process_info.my_rank;

  fn_fail:
    goto fn_exit;
}

int A1D_Size()
{
    A1U_FUNC_ENTER();

  fn_exit:
    A1U_FUNC_EXIT();
    return A1D_Process_info.num_ranks;

  fn_fail:
    goto fn_exit;
}

double A1D_Time_seconds()
{
    double time;

  fn_exit:
    A1U_FUNC_EXIT();
    return MPI_Wtime();

  fn_fail:
    goto fn_exit;
} 


unsigned long long A1D_Time_cycles()
{
    A1U_FUNC_ENTER();

    /* FIXME: implement this function using Kaz's ASM */

    A1U_error_printf("A1D_Time_cycles not implemented.\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return 0;

  fn_fail:
    goto fn_exit;
}
