/********************************************************************
 * The following is a notice of limited availability of the code, and disclaimer
 * which must be included in the prologue of the code and in all source listings
 * of the code.
 *
 * Copyright (c) 2010 Argonne Leadership Computing Facility, Argonne National Laboratory
 *
 * Permission is hereby granted to use, reproduce, prepare derivative works, and
 * to redistribute to others.
 *
 *                 LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer listed
 *    in this license in the documentation and/or other materials
 *    provided with the distribution.
 *
 *  - Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * The copyright holders provide no reassurances that the source code
 * provided does not infringe any patent, copyright, or any other
 * intellectual property rights of third parties.  The copyright holders
 * disclaim any liability to any recipient for claims brought against
 * recipient by any third party for infringement of that parties
 * intellectual property rights.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************/

#include "a1d_headers.h"
#include "a1d_globals.h"

#include "a1d_api.h"

int mpi_rank;
int mpi_size;

MPI_Comm A1D_COMM_WORLD;
void ** A1D_Baseptr_list;

#ifdef __bgp__
DCMF_Memregion_t * A1D_Memregion_list;
DCMF_Callback_t A1D_Nocallback;
#endif

#ifdef FLUSH_IMPLEMENTED
int* A1D_Put_flush_list;
int* A1D_Send_flush_list;
#endif


int A1D_Rank()
{
    return mpi_rank;
}

int A1D_Size()
{
    return mpi_size;
}

int A1D_Initialize()
{
    int mpi_initialized, mpi_provided;
    int mpi_status;
    int i;
    size_t bytes_in, bytes_out;
#ifdef __bgp__
    DCMF_Result dcmf_result;
    DCMF_Configure_t dcmf_config;
    DCMF_Memregion_t local_memregion;
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Initialize() \n");
#endif

    /***************************************************
     *
     * configure MPI
     *
     ***************************************************/

    /* MPI has to be Initialized for this implementation to work */
    MPI_Initialized(&mpi_initialized);
    assert(mpi_initialized==1);

    /* MPI has to be thread-safe so that DCMF doesn't explode */
    MPI_Query_thread(&mpi_provided);
    assert(mpi_provided==MPI_THREAD_MULTIPLE);

    /* have to use our own communicator for collectives to be proper */
    mpi_status = MPI_Comm_dup(MPI_COMM_WORLD,&A1D_COMM_WORLD);
    assert(mpi_status==0);

    /* get my MPI rank */
    mpi_status = MPI_Comm_rank(A1D_COMM_WORLD,&mpi_rank);
    assert(mpi_status==0);

    /* get MPI world size */
    mpi_status = MPI_Comm_size(A1D_COMM_WORLD,&mpi_size);
    assert(mpi_status==0);

    /* barrier before DCMF_Messager_configure to make sure MPI is ready everywhere */
    mpi_status = MPI_Barrier(A1D_COMM_WORLD);
    assert(mpi_status==0);

    /***************************************************
     *
     * configure DCMF
     *
     ***************************************************/

#ifdef __bgp__
    DCMF_CriticalSection_enter(0);

    /* make sure MPI and DCMF agree */
    assert(mpi_rank==(int)DCMF_Messager_rank());
    assert(mpi_size==(int)DCMF_Messager_size());

    /* probably not necessary since we require MPI_THREAD_MULTIPLE. */
    dcmf_config.thread_level = DCMF_THREAD_MULTIPLE;
    /* interrupts required for accumulate only, Put/Get use DMA
     * if accumulate not used, MPI will query environment for DCMF_INTERRUPTS */
    dcmf_config.interrupts = DCMF_INTERRUPTS_ON;

    /* reconfigure DCMF with interrupts on */
    dcmf_result = DCMF_Messager_configure(&dcmf_config, &dcmf_config);
    assert(dcmf_result==DCMF_SUCCESS);

    DCMF_CriticalSection_exit(0);
#endif

    /* barrier after DCMF_Messager_configure to make sure everyone has the new DCMF config */
    mpi_status = MPI_Barrier(A1D_COMM_WORLD);
    assert(mpi_status==0);

    /***************************************************
     *
     * setup DCMF memregions
     *
     ***************************************************/

#ifdef __bgp__
    DCMF_CriticalSection_enter(0);

    /* allocate memregion list */
    A1D_Memregion_list = malloc( mpi_size * sizeof(DCMF_Memregion_t) );
    assert(A1D_Memregion_list != NULL);

    /* allocate base pointer list */
    A1D_Baseptr_list = malloc( mpi_size * sizeof(void*) );
    assert(A1D_Memregion_list != NULL);

    /* create memregions */
    bytes_in = -1;
    dcmf_result = DCMF_Memregion_create(&local_memregion,&bytes_out,bytes_in,NULL,0);
    assert(dcmf_result==DCMF_SUCCESS);

    DCMF_CriticalSection_exit(0);

    /* exchange memregions because we don't use symmetry heap */
    mpi_status = MPI_Allgather(&local_memregion,sizeof(DCMF_Memregion_t),MPI_BYTE,
                               A1D_Memregion_list,sizeof(DCMF_Memregion_t),MPI_BYTE,
                               A1D_COMM_WORLD);
    assert(mpi_status==0);

    DCMF_CriticalSection_enter(0);

    /* destroy temporary local memregion */
    dcmf_result = DCMF_Memregion_destroy(&local_memregion);
    assert(dcmf_result==DCMF_SUCCESS);

    /* check for valid memregions */
    for (i = 0; i < mpi_size; i++)
    {
        dcmf_result = DCMF_Memregion_query(&A1D_Memregion_list[i],
                                           &bytes_out,
                                           &A1D_Baseptr_list[i]);
        assert(dcmf_result==DCMF_SUCCESS);
    }

    /***************************************************
     *
     * setup protocols and flush list(s)
     *
     ***************************************************/

    A1DI_Atomic_Initialize();

    A1DI_Get_Initialize();

    A1DI_Put_Initialize();
#  ifdef FLUSH_IMPLEMENTED
    /* allocate Put list */
    A1D_Put_flush_list = malloc( mpi_size * sizeof(int) );
    assert(A1D_Put_flush_list != NULL);
#  endif

    A1DI_Acc_Initialize();
#  ifdef FLUSH_IMPLEMENTED
    /* allocate Acc list */
    A1D_Send_flush_list = malloc( mpi_size * sizeof(int) );
    assert(A1D_Send_flush_list != NULL);
#  endif

    /***************************************************
     *
     * define null callback
     *
     ***************************************************/

    A1D_Nocallback.function = NULL;
    A1D_Nocallback.clientdata = NULL;

    DCMF_CriticalSection_exit(0);

#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Initialize() \n");
#endif

    return(0);
}

int A1D_Finalize()
{
    int mpi_status;
    int i;
#ifdef __bgp__
    DCMF_Result dcmf_result;
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Finalize() \n");
#endif

    A1D_Print_stats();

    /* barrier so that no one is able to access remote memregions after they are destroyed */
    mpi_status = MPI_Barrier(A1D_COMM_WORLD);
    assert(mpi_status==0);

#ifdef __bgp__
    DCMF_CriticalSection_enter(0);

#ifdef FLUSH_IMPLEMENTED
    /* free Put list */
    free(A1D_Put_flush_list);
    /* free Acc list */
    free(A1D_Send_flush_list);
#endif

    /* destroy all memregions - not absolutely unnecessary if memregion creation has no side effects */
    for (i = 0; i < mpi_size; i++)
    {
        dcmf_result = DCMF_Memregion_destroy(&A1D_Memregion_list[i]);
        assert(dcmf_result==DCMF_SUCCESS);
    }

    /* free memregion list */
    free(A1D_Memregion_list);

    /* free base pointer list */
    free(A1D_Baseptr_list);

    DCMF_CriticalSection_exit(0);
#endif

    /* free the A1D communicator */
    mpi_status = MPI_Comm_free(&A1D_COMM_WORLD);
    assert(mpi_status==0);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Finalize() \n");
#endif

    return(0);
}

/***************************************************
 *
 * local memory allocation
 *
 ***************************************************/

void * A1D_Allocate_local(int bytes)
{
    void * tmp;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Allocate_local(void ** ptr, int bytes) \n");
#endif

    if (bytes>0) 
    {
        tmp = calloc(bytes,1);
        assert( tmp != NULL );
    }
    else
    {
        if (bytes<0)
        {
            fprintf(stderr, "You requested %d bytes.  What kind of computer do you think I am? \n",bytes);
            fflush(stderr);
        }
        tmp = NULL;
    }

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Allocate_local(void ** ptr, int bytes) \n");
#endif

    return tmp;
}

void A1D_Free_local(void * ptr)
{
#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Free_local(void* ptr) \n");
#endif

    if (ptr != NULL)
    {
        free(ptr);
    }
    else
    {
        fprintf(stderr, "You tried to free a NULL pointer.  Please check your code. \n");
        fflush(stderr);
    }

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Free_local(void* ptr) \n");
#endif

    return;
}

/***************************************************
 *
 * global shared memory allocation
 *
 ***************************************************/

int A1D_Allocate_shared(void * ptrs[], int bytes)
{
    int mpi_status;
    void * tmp_ptr;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Allocate_shared(void* ptrs[], int bytes) \n");
#endif

    tmp_ptr = A1D_Allocate_local(bytes);

    mpi_status = MPI_Allgather(&tmp_ptr, sizeof(void *), MPI_BYTE,
                               ptrs,     sizeof(void *), MPI_BYTE,
                               A1D_COMM_WORLD);
    assert(mpi_status==0);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Allocate_shared(void* ptrs[], int bytes) \n");
#endif

    return(0);
}


void A1D_Free_shared(void * ptr)
{
    int mpi_status;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Free_shared(void* ptr) \n");
#endif

    /* barrier so that no one tries to access memory which is no longer allocated */
    mpi_status = MPI_Barrier(A1D_COMM_WORLD);
    assert(mpi_status==0);

    A1D_Free_local(ptr);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Free_shared(void* ptr) \n");
#endif

    return;
}

/***************************************************
 *
 * communicator-based shared memory allocation
 *
 ***************************************************/

int A1D_Create_window(const MPI_Comm comm, int bytes, A1D_Window_t* window)
{
    int mpi_status;
    int mpi_size;
    int mpi_rank;
    void * tmp_ptr;
    MPI_Comm newcomm;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Create_window(const MPI_Comm comm, int bytes, A1D_Window_t* window) \n");
#endif

    /* save (dup) the communicator into the window object */
    mpi_status = MPI_Comm_dup(comm,&newcomm);
    assert(mpi_status==0);

    window->comm = newcomm;

    /* need array sizeof(comm) for now */
    mpi_status = MPI_Comm_size(window->comm,&mpi_size);
    assert(mpi_status==0);

    /* my rank in this communicator */
    mpi_status = MPI_Comm_rank(window->comm,&mpi_rank);
    assert(mpi_status==0);

    /* allocate list of base pointers for this window */
    window->addr_list = malloc( mpi_size * sizeof(void *) );
    assert(window->addr_list != NULL);

    /* allocate local memory for this window */
    tmp_ptr = A1D_Allocate_local(bytes);

    /* exchange base pointers */
    mpi_status = MPI_Allgather(&tmp_ptr,          sizeof(void *), MPI_BYTE,
                               window->addr_list, sizeof(void *), MPI_BYTE,
                               window->comm);
    assert(mpi_status==0);

#ifndef NO_WINDOW_BOUNDS_CHECKING
    /* allocate list of sizes */
    window->addr_list = malloc(mpi_size*sizeof(int));
    assert(window->size_list != NULL);

    /* exchange sizes pointers */
    mpi_status = MPI_Allgather(&bytes,            sizeof(int), MPI_BYTE,
                               window->size_list, sizeof(int), MPI_BYTE,
                               window->comm);
    assert(mpi_status==0);

#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Create_window(const MPI_Comm comm, int bytes, A1D_Window_t* window) \n");
#endif

    return(0);
}

int A1D_Destroy_window(A1D_Window_t* window)
{
    int mpi_status;
    int mpi_rank;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Destroy_window(A1D_Window_t* window) \n");
#endif

    /* barrier so that no one is able to access remote window memory after it is free */
    mpi_status = MPI_Barrier(window->comm);
    assert(mpi_status==0);

    /* my rank in this communicator */
    mpi_status = MPI_Comm_rank(window->comm,&mpi_rank);
    assert(mpi_status==0);

    /* free the local memory */
    A1D_Free_local(window->addr_list[mpi_rank]);

    /* free the list of base pointers */
    free(window->addr_list);

#ifndef NO_WINDOW_BOUNDS_CHECKING
    /* free list of sizes */
    free(window->size_list);
#endif

    /* free the communicator */
    mpi_status = MPI_Comm_free(&(window->comm));
    assert(mpi_status==0);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Destroy_window(A1D_Window_t* window) \n");
#endif

    return(0);
}

