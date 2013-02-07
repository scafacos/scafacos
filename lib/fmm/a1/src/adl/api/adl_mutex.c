/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1d.h"
#include "a1u.h"

int A1_Create_mutexes(A1_group_t* group, 
                      int mutex_count, 
                      int *mutex_count_ar)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1D_Create_mutexes(group,  
                                mutex_count, 
                                mutex_count_ar);
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Create_mutexes returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}


int A1_Destroy_mutexes(A1_group_t* group)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1D_Destroy_mutexes(group);
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Destroy_mutexes returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_Lock_mutex(A1_group_t* group, 
                  int mutex, 
                  int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1D_Lock_mutex(group, mutex, proc);
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Lock_mutex returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_Trylock_mutex(A1_group_t* group, 
                     int mutex, 
                     int proc,
                     A1_bool_t *acquired)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1D_Trylock_mutex(group, mutex, proc, acquired);
    A1U_ERR_POP(status, "A1D_Lock_mutex returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_Unlock_mutex(A1_group_t* group, 
                    int mutex, 
                    int proc)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

#   ifdef HAVE_ERROR_CHECKING
#   endif

    status = A1D_Unlock_mutex(group, mutex, proc);
    A1U_ERR_POP(status != A1_SUCCESS, "A1D_Unlock_mutex returned an error\n");

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
