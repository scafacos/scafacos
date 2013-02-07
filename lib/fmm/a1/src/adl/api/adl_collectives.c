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

#define A1_USES_MPI_COLLECTIVES

#ifdef A1_USES_MPI_COLLECTIVES
#include "mpi.h"
#endif

#define ABS(datatype, source, target, count)                                \
   do {                                                                     \
     int i;                                                                 \
     datatype *s = (datatype *) source;                                     \
     datatype *t = (datatype *) target;                                     \
     for(i=0; i<count; i++) t[i] = ( s[i] > 0 ? s[i] : -s[i]);              \
   } while(0)                                                               \

int A1_Barrier_group(A1_group_t* group)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

#ifdef A1_USES_MPI_COLLECTIVES

    if (group == A1_GROUP_WORLD || group == NULL)
    {
        status = MPI_Barrier(MPI_COMM_WORLD);
        switch (status)
        {
            case MPI_ERR_COMM:
                A1U_ERR_POP(1,"MPI_Barrier returned MPI_ERR_COMM.");
                break;
            default:
                status = A1_SUCCESS;
                goto fn_exit;
                break;
        }
    }
    else
    {
        A1U_ERR_POP(1,"A1_Barrier_group not implemented for non-world groups!");
    }

#else

    /* barrier is meaningless with 1 process */
    if (1==A1D_Process_total(A1_GROUP_WORLD)) goto fn_exit;

    status = A1D_Barrier_group(group);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Barrier_group returned an error\n");

#endif

  fn_exit: 
    A1U_FUNC_EXIT();
    return status;

  fn_fail: 
    goto fn_exit;
}

int A1_NbBarrier_group(A1_group_t* group, A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

#ifdef A1_USES_MPI_COLLECTIVES

    A1U_ERR_POP(1,"A1_NbBarrier_group not implemented for when A1_USES_MPI_COLLECTIVES is defined.");

#else

    /* barrier is meaningless with 1 process */
    if ( 1==A1D_Process_total(A1_GROUP_WORLD) ) goto fn_exit;

    status = A1D_NbBarrier_group(group, a1_handle);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_NbBarrier_group returned an error\n");

#endif

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_Sync_group(A1_group_t* group)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

#ifdef A1_USES_MPI_COLLECTIVES

    status = A1D_Flush_group(group);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Flush_group returned an error\n");

    if (group == A1_GROUP_WORLD || group == NULL)
    {
        status = MPI_Barrier(MPI_COMM_WORLD);
        switch (status)
        {
            case MPI_ERR_COMM:
                A1U_ERR_POP(1,"MPI_Barrier returned MPI_ERR_COMM.");
                break;
            default:
                status = A1_SUCCESS;
                goto fn_exit;
                break;
        }
    }
    else
    {
        A1U_ERR_POP(1,"A1_Sync_group not implemented for non-world groups!");
    }

#else

    /* no collective bypass for 1 proc here because we will use DCMF for some operations
     * which need to be completed by flush */

    status = A1D_Sync_group(group);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Sync_group returned an error\n");

#endif


  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_NbSync_group(A1_group_t* group, A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

#ifdef A1_USES_MPI_COLLECTIVES

    A1U_ERR_POP(1,"A1_NbSync_group not implemented for when A1_USES_MPI_COLLECTIVES is defined.");

#else

    /* no collective bypass for 1 proc here because we will use DCMF for some operations
     * which need to be completed by flush */

    status = A1D_NbSync_group(group, a1_handle);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_NbSync_group returned an error\n");

#endif

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_Allreduce_group(A1_group_t* group,
                       int count,
                       A1_reduce_op_t a1_op,
                       A1_datatype_t a1_type,
                       void* in,
                       void* out)
{
#ifdef A1_USES_MPI_COLLECTIVES
    MPI_Datatype mpi_type;
    MPI_Op mpi_oper;
    int bytes;
    void *in_absolute = NULL;
#endif /* A1_USES_MPI_COLLECTIVES */
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

#ifdef A1_USES_MPI_COLLECTIVES
    if (group == A1_GROUP_WORLD || group == NULL)
    {
        switch (a1_type)
        {
            case A1_DOUBLE:
                mpi_type = MPI_DOUBLE;
                break;
            case A1_INT32:
                mpi_type = MPI_LONG;
                break;
            case A1_INT64:
                mpi_type = MPI_LONG_LONG;
                break;
            case A1_UINT32:
                mpi_type = MPI_UNSIGNED_LONG;
                break;
            case A1_UINT64:
                mpi_type = MPI_UNSIGNED_LONG_LONG;
                break;
            case A1_FLOAT:
                mpi_type = MPI_FLOAT;
                break;
            default:
                A1U_ERR_POP(status!=A1_SUCCESS, "Unsupported A1_datatype\n");
                break;
        }

        switch (a1_op)
        {
            case A1_SUM:
                mpi_oper = MPI_SUM;
                break;
            case A1_PROD:
                mpi_oper = MPI_PROD;
                break;
            case A1_MAX:
                mpi_oper = MPI_MAX;
                break;
            case A1_MIN:
                mpi_oper = MPI_MIN;
                break;
            case A1_OR:
                mpi_oper = MPI_LOR;
                break;
            case A1_MAXABS:
                mpi_oper = MPI_MAX;
                break;
            case A1_MINABS:
                mpi_oper = MPI_MIN;
                break;
            case A1_SAME:
                A1U_ERR_POP(1, "A1_SAME is not supported when A1_USES_MPI_COLLECTIVES is defined.\n");
                break;
            default:
                A1U_ERR_POP(status!=A1_SUCCESS, "Unsupported A1_op\n");
                break;
        }
 
        if(a1_op == A1_MAXABS || a1_op == A1_MINABS)
        switch (a1_type)
        {
            case A1_DOUBLE:
                bytes = count * sizeof(double);
                in_absolute = malloc(bytes);
                A1U_ERR_POP(in_absolute == NULL,
                            "malloc returned error in A1_Allreduce_group \n");
                ABS(double, in, in_absolute, count);
                in = in_absolute;
                break;
            case A1_INT32:
                bytes = count * sizeof(int32_t);
                in_absolute = malloc(bytes);
                A1U_ERR_POP(in_absolute == NULL,
                            "malloc returned error in A1_Allreduce_group \n");
                ABS(int32_t, in, in_absolute, count);
                in = in_absolute;
                break;
            case A1_INT64:
                bytes = count * sizeof(int64_t);
                in_absolute = malloc(bytes);
                A1U_ERR_POP(in_absolute == NULL,
                            "malloc returned error in A1_Allreduce_group \n");
                ABS(int64_t, in, in_absolute, count);
                in = in_absolute;
                break;
            case A1_UINT32:
                break;
            case A1_UINT64:
                break;
            case A1_FLOAT:
                bytes = count * sizeof(float);
                in_absolute = malloc(bytes);
                A1U_ERR_POP(in_absolute == NULL,
                            "malloc returned error in A1_Allreduce_group \n");
                ABS(float, in, in_absolute, count);
                in = in_absolute;
                break;
            default:
                status = A1_ERROR;
                A1U_ERR_POP(status != A1_SUCCESS, "Unsupported A1_datatype \n");
                break;
        }

        if (in==out) in = MPI_IN_PLACE;
        status = MPI_Allreduce(in,out,count,mpi_type,mpi_oper,MPI_COMM_WORLD);
        switch (status)
        {
            case MPI_ERR_BUFFER:
                A1U_ERR_POP(1,"MPI_Allreduce returned MPI_ERR_BUFFER.");
                break;
            case MPI_ERR_COUNT:
                A1U_error_printf("count = %d\n",count);
                A1U_ERR_POP(1,"MPI_Allreduce returned MPI_ERR_COUNT.");
                break;
            case MPI_ERR_TYPE:
                A1U_ERR_POP(1,"MPI_Allreduce returned MPI_ERR_TYPE.");
                break;
            case MPI_ERR_OP:
                A1U_ERR_POP(1,"MPI_Allreduce returned MPI_ERR_OP.");
                break;
            case MPI_ERR_COMM:
                A1U_ERR_POP(1,"MPI_Allreduce returned MPI_ERR_COMM.");
                break;
            default:
                status = A1_SUCCESS;
                goto fn_exit;
                break;
        }
    }
    else
    {
        A1U_ERR_POP(1,"A1_Allreduce_group not implemented for non-world groups!");
    }

#else

    if (count <= 0) goto fn_exit;

    /* bypass any sort of network API or communication altogether */
    if ( 1==A1D_Process_total(A1_GROUP_WORLD) )
    {
        switch (a1_type)
        {
            case A1_DOUBLE:
                memcpy(in,out,count*sizeof(double));
                break;
            case A1_INT32:
                memcpy(in,out,count*sizeof(int32_t));
                break;
            case A1_INT64:
                memcpy(in,out,count*sizeof(int64_t));
                break;
            case A1_UINT32:
                memcpy(in,out,count*sizeof(uint32_t));
                break;
            case A1_UINT64:
                memcpy(in,out,count*sizeof(uint64_t));
                break;
            case A1_FLOAT:
                memcpy(in,out,count*sizeof(float));
                break;
            default:
                A1U_ERR_POP(status!=A1_SUCCESS, "Unsupported A1_datatype\n");
                break;
        }
        goto fn_exit;
    }

    status = A1D_Allreduce_group(group,
                                 count,
                                 a1_op,
                                 a1_type,
                                 in,
                                 out);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Allreduce_group returned an error\n");

#endif

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_NbAllreduce_group(A1_group_t* group,
                         int count,
                         A1_reduce_op_t a1_op,
                         A1_datatype_t a1_type,
                         void* in,
                         void* out,
                         A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

#ifdef A1_USES_MPI_COLLECTIVES

    A1U_ERR_POP(1,"A1_NbAllreduce_group not implemented for when A1_USES_MPI_COLLECTIVES is defined.");

#else

    if (count <= 0) goto fn_exit;

    /* bypass any sort of network API or communication altogether */
    if ( 1==A1D_Process_total(A1_GROUP_WORLD) )
    {
        switch (a1_type)
        {
            case A1_DOUBLE:
                memcpy(in,out,count*sizeof(double));
                break;
            case A1_INT32:
                memcpy(in,out,count*sizeof(int32_t));
                break;
            case A1_INT64:
                memcpy(in,out,count*sizeof(int64_t));
                break;
            case A1_UINT32:
                memcpy(in,out,count*sizeof(uint32_t));
                break;
            case A1_UINT64:
                memcpy(in,out,count*sizeof(uint64_t));
                break;
            case A1_FLOAT:
                memcpy(in,out,count*sizeof(float));
                break;
            default:
                A1U_ERR_POP(status!=A1_SUCCESS, "Unsupported A1_datatype\n");
                break;
        }
        goto fn_exit;
    }

    status = A1D_NbAllreduce_group(group,
                                   count,
                                   a1_op,
                                   a1_type,
                                   in,
                                   out,
                                   a1_handle);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_NbAllreduce_group returned an error\n");

#endif

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_Bcast_group(A1_group_t* group,
                   int root,
                   int count,
                   void* buffer)
{
#ifdef A1_USES_MPI_COLLECTIVES
    MPI_Datatype mpi_type = MPI_BYTE;
#endif /* A1_USES_MPI_COLLECTIVES */
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

#ifdef A1_USES_MPI_COLLECTIVES
    if (group == A1_GROUP_WORLD || group == NULL)
    {
        status = MPI_Bcast(buffer,count,mpi_type,root,MPI_COMM_WORLD);
        switch (status)
        {
            case MPI_ERR_BUFFER:
                A1U_ERR_POP(1,"MPI_Bcast returned MPI_ERR_BUFFER.");
                break;
            case MPI_ERR_COUNT:
                A1U_error_printf("count = %d\n",count);
                A1U_ERR_POP(1,"MPI_Bcast returned MPI_ERR_COUNT.");
                break;
            case MPI_ERR_TYPE:
                A1U_ERR_POP(1,"MPI_Bcast returned MPI_ERR_TYPE.");
                break;
            case MPI_ERR_ROOT:
                A1U_ERR_POP(1,"MPI_Bcast returned MPI_ERR_ROOT.");
                break;
            case MPI_ERR_COMM:
                A1U_ERR_POP(1,"MPI_Bcast returned MPI_ERR_COMM.");
                break;
            default:
                status = A1_SUCCESS;
                goto fn_exit;
                break;
        }
    }
    else
    {
        A1U_ERR_POP(1,"A1_Barrier_group not implemented for non-world groups!");
    }

#else

    if (count <= 0) goto fn_exit;

    /* bypass any sort of network API or communication altogether */
    if ( 1==A1D_Process_total(A1_GROUP_WORLD) ) goto fn_exit;

    status = A1D_Bcast_group(group,
                             root,
                             count,
                             buffer);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_Bcast_group returned an error\n");

#endif

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}

int A1_NbBcast_group(A1_group_t* group,
                     int root,
                     int count,
                     void* buffer,
                     A1_handle_t a1_handle)
{
    int status = A1_SUCCESS;

    A1U_FUNC_ENTER();

    /* FIXME: The profiling interface needs to go here */

    /* FIXME: Locking functionality needs to go here */

#   ifdef HAVE_ERROR_CHECKING
#   endif

#ifdef A1_USES_MPI_COLLECTIVES

    A1U_ERR_POP(1,"A1_NbBcast_group not implemented for when A1_USES_MPI_COLLECTIVES is defined.");

#else

    if (count <= 0) goto fn_exit;

    /* bypass any sort of network API or communication altogether */
    if ( 1==A1D_Process_total(A1_GROUP_WORLD) ) goto fn_exit;

    status = A1D_NbBcast_group(group,
                               root,
                               count,
                               buffer,
                               a1_handle);
    A1U_ERR_POP(status!=A1_SUCCESS, "A1D_NbBcast_group returned an error\n");

#endif

  fn_exit:
    A1U_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
