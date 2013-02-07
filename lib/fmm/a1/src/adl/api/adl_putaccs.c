/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

/* This is here because the build system does not yet have the necessary
 * logic to set these options for each device. */
#define A1D_IMPLEMENTS_PUTACCS

#if defined A1D_IMPLEMENTS_PUTACCS

int A1_PutAccS(int target,
               int stride_level,
               int *block_sizes,
               void* source_ptr,
               int *src_stride_ar,
               void* target_ptr,
               int *trg_stride_ar,
               A1_datatype_t a1_type,
               void* scaling)
{
    int status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
      int i, bytes = 1;
      for (i = 0; i <= stride_levels; i++) total_bytes *= count[i];
      TAU_TRACE_SENDMSG (A1_TAU_TAG_PUTACCS, target, total_bytes);
    }
#   endif

    /*Check if it is a contiguous transfer, issue a contiguous op*/
    if(stride_level == 0)
    {
        /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
        if (target == my_rank && a1u_settings.network_bypass)
        {
           status = A1U_Acc_memcpy(source_ptr,
                                   target_ptr,
                                   block_sizes[0],
                                   a1_type,
                                   scaling);
           A1U_ERR_POP(status != A1_SUCCESS, "A1U_Acc_memcpy returned an error\n");
        }
        else
        { 
           status = A1D_PutAcc(target,
                               source_ptr,
                               target_ptr,
                               block_sizes[0],
                               a1_type,
                               scaling);
           A1U_ERR_POP((status!=A1_SUCCESS), "A1D_PutAcc returned error\n");
        }
        goto fn_exit;
    }

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_AccS_memcpy(stride_level,
                                 block_sizes,
                                 source_ptr,
                                 src_stride_ar,
                                 target_ptr,
                                 trg_stride_ar,
                                 a1_type,
                                 scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_AccS_memcpy returned an error\n");
    }
    else
    {
        status = A1D_PutAccS(target,
                             stride_level,
                             block_sizes,
                             source_ptr,
                             src_stride_ar,
                             target_ptr,
                             trg_stride_ar,
                             a1_type,
                             scaling);
        A1U_ERR_POP(status, "A1D_PutAccS returned error\n");
    }

    fn_exit:
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

int A1_NbPutAccS(int target,
                 int stride_level,
                 int *block_sizes,
                 void* source_ptr,
                 int *src_stride_ar,
                 void* target_ptr,
                 int *trg_stride_ar,
                 A1_datatype_t a1_type,
                 void* scaling,
                 A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
      int i, bytes = 1;
      for (i = 0; i <= stride_levels; i++) total_bytes *= count[i];
      TAU_TRACE_SENDMSG (A1_TAU_TAG_NBPUTACCS, target, total_bytes);
    }
#   endif

    /*Check if it is a contiguous transfer, issue a contiguous op*/
    if(stride_level == 0)
    {
        /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
        if (target == my_rank && a1u_settings.network_bypass)
        {
           status = A1U_Acc_memcpy(source_ptr,
                                   target_ptr,
                                   block_sizes[0],
                                   a1_type,
                                   scaling);
           A1U_ERR_POP(status != A1_SUCCESS, "A1U_Acc_memcpy returned an error\n");
        }
        else
        {
           status = A1D_NbPutAcc(target,
                                 source_ptr,
                                 target_ptr,
                                 block_sizes[0],
                                 a1_type,
                                 scaling,
		   	         a1_handle);
           A1U_ERR_POP((status!=A1_SUCCESS), "A1D_PutAcc returned error\n");
        }
        goto fn_exit;
    }

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_AccS_memcpy(stride_level,
                                 block_sizes,
                                 source_ptr,
                                 src_stride_ar,
                                 target_ptr,
                                 trg_stride_ar,
                                 a1_type,
                                 scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_AccS_memcpy returned an error\n");
    }
    else
    {
        status = A1D_NbPutAccS(target,
                               stride_level,
                               block_sizes,
                               source_ptr,
                               src_stride_ar,
                               target_ptr,
                               trg_stride_ar,
                               a1_type,
                               scaling,
                               a1_handle);
        A1U_ERR_POP(status, "NbPutAccS returned error\n");
    }

    fn_exit:
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

#else

int A1I_Recursive_PutAcc(int target,
                         int stride_level,
                         int *block_sizes,
                         void* source_ptr,
                         int *src_stride_ar,
                         void* target_ptr,
                         int *trg_stride_ar,
                         A1_datatype_t a1_type,
                         void* scaling,
                         A1_handle_t a1_handle)
{
    int i, status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    if (stride_level > 0)
    {
        for (i = 0; i < block_sizes[stride_level]; i++)
        {
            status = A1I_Recursive_PutAcc(target,
                                          stride_level - 1,
                                          block_sizes,
                                          (void *) ((size_t) source_ptr + i * src_stride_ar[stride_level - 1]),
                                          src_stride_ar,
                                          (void *) ((size_t) target_ptr + i * trg_stride_ar[stride_level - 1]),
                                          trg_stride_ar,
                                          a1_type,
                                          scaling,
                                          a1_handle);
            A1U_ERR_POP(status != A1_SUCCESS,
                        "A1I_Recursive_PutAcc returned error in A1I_Recursive_PutAcc.\n");
        }
    }
    else
    {
        status = A1D_NbPutAcc(target,
                              source_ptr,
                              target_ptr,
                              block_sizes[0],
                              a1_type,
                              scaling,
                              a1_handle);
        A1U_ERR_POP(status != A1_SUCCESS, "A1D_NbPutAcc returned with an error \n");
    }

    fn_exit:
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

int A1_PutAccS(int target,
               int stride_level,
               int *block_sizes,
               void* source_ptr,
               int *src_stride_ar,
               void* target_ptr,
               int *trg_stride_ar,
               A1_datatype_t a1_type,
               void* scaling)
{
    int status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);
    A1_handle_t a1_handle;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
      int i, bytes = 1;
      for (i = 0; i <= stride_levels; i++) total_bytes *= count[i];
      TAU_TRACE_SENDMSG (A1_TAU_TAG_PUTACCS, target, total_bytes);
    }
#   endif

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_AccS_memcpy(stride_level,
                                 block_sizes,
                                 source_ptr,
                                 src_stride_ar,
                                 target_ptr,
                                 trg_stride_ar,
                                 a1_type,
                                 scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_AccS_memcpy returned an error\n");
    }
    else
    {
        status = A1D_Allocate_handle(&a1_handle);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Allocate_handle returned error\n");

        status = A1I_Recursive_PutAcc(target,
                                      stride_level,
                                      block_sizes,
                                      source_ptr,
                                      src_stride_ar,
                                      target_ptr,
                                      trg_stride_ar,
                                      a1_type,
                                      scaling,
                                      a1_handle);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1I_Recursive_PutAcc returned error\n");

        status = A1D_Wait_handle(a1_handle);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Wait_handle returned error\n");
    }

    fn_exit:
    /* Could also test for NULL, assuming we set it as such in the declaration. */
    if(target == my_rank && a1u_settings.network_bypass) A1D_Release_handle(a1_handle);
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

int A1_NbPutAccS(int target,
                 int stride_level,
                 int *block_sizes,
                 void* source_ptr,
                 int *src_stride_ar,
                 void* target_ptr,
                 int *trg_stride_ar,
                 A1_datatype_t a1_type,
                 void* scaling,
                 A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;
    int my_rank = A1D_Process_id(A1_GROUP_WORLD);

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

#   ifdef A1_TAU_PROFILING
    {
      int i, bytes = 1;
      for (i = 0; i <= stride_levels; i++) total_bytes *= count[i];
      TAU_TRACE_SENDMSG (A1_TAU_TAG_NBPUTACCS, target, total_bytes);
    }
#   endif

    /* Bypass is ALWAYS better for accumulate; we do not test against threshold. */
    if (target == my_rank && a1u_settings.network_bypass)
    {
        status = A1U_AccS_memcpy(stride_level,
                                 block_sizes,
                                 source_ptr,
                                 src_stride_ar,
                                 target_ptr,
                                 trg_stride_ar,
                                 a1_type,
                                 scaling);
        A1U_ERR_POP(status != A1_SUCCESS, "A1U_AccS_memcpy returned an error\n");
    }
    else
    {
        status = A1I_Recursive_PutAcc(target,
                                      stride_level,
                                      block_sizes,
                                      source_ptr,
                                      src_stride_ar,
                                      target_ptr,
                                      trg_stride_ar,
                                      a1_type,
                                      scaling,
                                      a1_handle);
        A1U_ERR_POP(status!=A1_SUCCESS, "A1I_Recursive_PutAcc returned error\n");
    }

    fn_exit:
    A1U_FUNC_EXIT();
    return status;

    fn_fail:
    goto fn_exit;
}

#endif /* A1D_IMPLEMENTS_PUTACCV */
