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

#ifdef DMAPPD_USES_MPI
  MPI_Comm A1D_COMM_WORLD;
#else
# error DMAPP requires MPI for now.
#endif

#ifdef __CRAYXE
  int A1D_Memdesc_list_size = 0;
  dmapp_seg_desc_t * A1D_Memdesc_list = NULL;
  dmapp_seg_desc_t A1D_Sheap_desc;
#endif

#ifdef FLUSH_IMPLEMENTED
  int32_t *  A1D_Put_flush_list;
#endif

int64_t * A1D_Acc_lock;

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

#ifdef DMAPPD_USES_MPI
    int mpi_initialized, mpi_provided;
    int mpi_status = MPI_SUCCESS;

    int namelen;
    char procname[MPI_MAX_PROCESSOR_NAME];
#endif

#ifdef __CRAYXE
    int                                 pmi_status  = PMI_SUCCESS;
    int                                 nodeid = -1;
    rca_mesh_coord_t                    rca_xyz;

    dmapp_return_t                      dmapp_status = DMAPP_RC_SUCCESS;

    dmapp_rma_attrs_ext_t               dmapp_config_in, dmapp_config_out;

    dmapp_jobinfo_t                     dmapp_info;
    dmapp_pe_t                          dmapp_rank = -1;
    int                                 dmapp_size = -1;
#endif
    int                                 sheapflag = 0;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Initialize() \n");
#endif

#ifdef DMAPPD_USES_MPI

    /***************************************************
     *
     * configure MPI
     *
     ***************************************************/

    /* MPI has to be Initialized for this implementation to work */
    MPI_Initialized(&mpi_initialized);
    assert(mpi_initialized==1);

    /* MPI has to tolerate threads because A1 supports them */
    MPI_Query_thread(&mpi_provided);
    //assert(mpi_provided>MPI_THREAD_SINGLE);

    /* have to use our own communicator for collectives to be proper */
    mpi_status = MPI_Comm_dup(MPI_COMM_WORLD,&A1D_COMM_WORLD);
    assert(mpi_status==0);

    /* get my MPI rank */
    mpi_status = MPI_Comm_rank(A1D_COMM_WORLD,&mpi_rank);
    assert(mpi_status==0);

    /* get MPI world size */
    mpi_status = MPI_Comm_size(A1D_COMM_WORLD,&mpi_size);
    assert(mpi_status==0);

    /* in a perfect world, this would provide topology information like BG */
    MPI_Get_processor_name( procname, &namelen );
    printf( "%d: MPI_Get_processor_name = %s\n" , mpi_rank, procname );
    fflush( stdout );

    /* barrier to make sure MPI is ready everywhere */
    mpi_status = MPI_Barrier(A1D_COMM_WORLD);
    assert(mpi_status==0);

#endif

#ifdef __CRAYXE

    /***************************************************
     *
     * query topology
     *
     ***************************************************/

    PMI_Get_nid( mpi_rank, &nodeid );
    assert(pmi_status==PMI_SUCCESS);

    rca_get_meshcoord((uint16_t)nodeid, &rca_xyz);
    printf("%d: rca_get_meshcoord returns (%2u,%2u,%2u)\n", mpi_rank, rca_xyz.mesh_x, rca_xyz.mesh_y, rca_xyz.mesh_z );

#endif

#ifdef __CRAYXE

    /***************************************************
     *
     * configure DMAPP
     *
     ***************************************************/

    dmapp_config_in.max_outstanding_nb   = DMAPP_DEF_OUTSTANDING_NB; /*  512 */
    dmapp_config_in.offload_threshold    = DMAPP_OFFLOAD_THRESHOLD;  /* 4096 */
#ifdef DETERMINISTIC_ROUTING
    dmapp_config_in.put_relaxed_ordering = DMAPP_ROUTING_DETERMINISTIC;
    dmapp_config_in.get_relaxed_ordering = DMAPP_ROUTING_DETERMINISTIC;
#else
    dmapp_config_in.put_relaxed_ordering = DMAPP_ROUTING_ADAPTIVE;
    dmapp_config_in.get_relaxed_ordering = DMAPP_ROUTING_ADAPTIVE;
#endif
    dmapp_config_in.max_concurrency      = 1; /* not thread-safe */
#ifdef FLUSH_IMPLEMENTED
    dmapp_config_in.PI_ordering          = DMAPP_PI_ORDERING_RELAXED;
#else
    dmapp_config_in.PI_ordering          = DMAPP_PI_ORDERING_STRICT;
#endif

    dmapp_status = dmapp_init_ext( &dmapp_config_in, &dmapp_config_out );
    assert(dmapp_status==DMAPP_RC_SUCCESS);

#ifndef FLUSH_IMPLEMENTED
    /* without strict PI ordering, we have to flush remote stores with a get packet to force global visibility */
    assert( dmapp_config_out.PI_ordering == DMAPP_PI_ORDERING_STRICT);
#endif

    dmapp_status = dmapp_get_jobinfo(&dmapp_info);
    assert(dmapp_status==DMAPP_RC_SUCCESS);

    dmapp_rank     = dmapp_info.pe;
    dmapp_size     = dmapp_info.npes;
    A1D_Sheap_desc = dmapp_info.sheap_seg;

    /* make sure PMI and DMAPP agree */
    assert(mpi_rank==dmapp_rank);
    assert(mpi_size==dmapp_size);

#endif

    /***************************************************
     *
     * setup protocols
     *
     ***************************************************/

#ifdef FLUSH_IMPLEMENTED
    /* allocate Put list */
    A1D_Put_flush_list = malloc( mpi_size * sizeof(int32_t) );
    assert(A1D_Put_flush_list != NULL);
#endif

#ifdef __CRAYXE
    A1D_Acc_lock = dmapp_sheap_malloc( sizeof(int64_t) );
#endif

    A1D_Allreduce_issame64((size_t)A1D_Acc_lock, &sheapflag);
    assert(sheapflag==1);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Initialize() \n");
#endif

    return(0);
}

int A1D_Finalize()
{
#ifdef __CRAYXE
    dmapp_return_t dmapp_status = DMAPP_RC_SUCCESS;
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Finalize() \n");
#endif

    A1D_Print_stats();

#ifdef FLUSH_IMPLEMENTED
    free(A1D_Put_flush_list);
#endif

    /* barrier so that no one is able to access remote memregions after they are destroyed */
    mpi_status = MPI_Barrier(A1D_COMM_WORLD);
    assert(mpi_status==0);

#ifdef __CRAYXE
    /* shut down DMAPP */
    dmapp_status = dmapp_finalize();
    assert(dmapp_status==DMAPP_RC_SUCCESS);
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Finalize() \n");
#endif

    return(0);
}

/***************************************************
 *
 * global shared memory allocation
 *
 ***************************************************/

int A1D_Allocate_shared(void * ptrs[], int bytes)
{
    void *  tmp_ptr       = NULL;
    int     max_bytes     = 0;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Allocate_shared(void* ptrs[], int bytes) \n");
#endif

    mpi_status = MPI_Barrier(A1D_COMM_WORLD);
    assert(mpi_status==0);

#ifdef __CRAYXE
    A1D_Allreduce_max32( bytes, &max_bytes );

    /* allocate memory from symmetric heap */
    tmp_ptr = dmapp_sheap_malloc( (size_t)max_bytes );
    assert(tmp_ptr!=NULL);
#endif

    /* allgather addresses into pointer vector */
    A1D_Allgather( &tmp_ptr, ptrs, sizeof(void*) );

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Allocate_shared(void* ptrs[], int bytes) \n");
#endif

    return(0);
}

void A1D_Free_shared(void * ptr)
{
#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Free_shared(void* ptr) \n");
#endif

    /* barrier so that no one tries to access memory which is no longer allocated
     * and to ensure that the user calls this function collectively */
    A1D_Barrier();

#ifdef __CRAYXE
    if (ptr != NULL)
    {
        dmapp_sheap_free(ptr);
    }
    else
    {
        fprintf(stderr, "You tried to free a NULL pointer.  Please check your code. \n");
        fflush(stderr);
    }
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Free_shared(void* ptr) \n");
#endif

    return;
}

/***************************************************
 *
 * local memory allocation
 *
 ***************************************************/

/* A1DI_Acquire_mem_descriptor and A1DI_Release_mem_descriptor
 * are definitely NOT thread-safe and should be used carefully */

#ifdef __CRAYXE
int A1DI_Acquire_mem_descriptor( dmapp_seg_desc_t * mem_desc )
{
#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Acquire_mem_descriptor( dmapp_seg_desc_t * mem_desc) \n");
#endif

    if (A1D_Memdesc_list_size==0)
    {
        A1D_Memdesc_list = (dmapp_seg_desc_t *) malloc( sizeof(dmapp_seg_desc_t) );
        assert(A1D_Memdesc_list!=NULL);

        A1D_Memdesc_list_size=1;
    }
    else
    {
        A1D_Memdesc_list = (dmapp_seg_desc_t *) realloc( list, (A1D_Memdesc_list_size+1) * sizeof(dmapp_seg_desc_t) );
        assert(A1D_Memdesc_list!=NULL);

        A1D_Memdesc_list_size++;
    }

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Acquire_mem_descriptor( dmapp_seg_desc_t * mem_desc) \n");
#endif

    return(0);
}

int A1DI_Release_mem_descriptor( void * address )
{
    dmapp_return_t     dmapp_status  = DMAPP_RC_SUCCESS;
    int                flag          = 0;
    int                offset        = -1;
    dmapp_seg_desc_t * temp_list     = NULL;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Release_mem_descriptor( dmapp_seg_desc_t * mem_desc) \n");
#endif

    assert(A1D_Memdesc_list_size>0);

    for ( int i=0 ; i<A1D_Memdesc_list_size ; i++)
        if ( A1D_Memdesc_list[i].addr == address )
            offset = i;

    if (offset == -1)
    {
        fprintf(stderr,"A1DI_Release_mem_descriptor: mem_desc not found! \n");
        assert(0); /* aborting isn't actually required but makes sense */
    }
    else
    {
        temp_list = (dmapp_seg_desc_t *) malloc( (A1D_Memdesc_list_size-1) * sizeof(dmapp_seg_desc_t) );
        assert(temp_list!=NULL);

        /* copy A1D_Memdesc_list sans mem_desc into temp_list */
        int j = 0;
        for ( int i=0 ; i<A1D_Memdesc_list_size ; i++)
            if (i!=offset)
                temp_list[j++] = A1D_Memdesc_list[i];

        /* resize primary, copy temp into it and decrement size */
        A1D_Memdesc_list = (dmapp_seg_desc_t *) realloc( list, (A1D_Memdesc_list_size-1) * sizeof(dmapp_seg_desc_t) );
        memcpy( A1D_Memdesc_list, temp_list,  (A1D_Memdesc_list_size-1) * sizeof(dmapp_seg_desc_t) );
        A1D_Memdesc_list_size--;

        free(temp_list);
    }

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Release_mem_descriptor( dmapp_seg_desc_t * mem_desc) \n");
#endif

    return(0);
}

#endif

void * A1D_Allocate_local(int bytes)
{
    void * tmp;
#ifdef __CRAYXE
    dmapp_return_t   dmapp_status = DMAPP_RC_SUCCESS;
    dmapp_seg_desc_t * mem_desc = NULL;
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Allocate_local(void ** ptr, int bytes) \n");
#endif

    if (bytes>0)
    {
        posix_memalign( &tmp, (size_t)DMAPP_ALIGNMENT, (size_t)bytes);
        assert(tmp!=NULL);

#ifdef __CRAYXE
        A1DI_Acquire_mem_descriptor( &mem_desc );

        dmapp_status = dmapp_mem_register( tmp, (uint64_t)bytes, &mem_desc );
        assert(dmapp_status==DMAPP_RC_SUCCESS);

        /* verify that DMAPP has registered at my address exactly */
        assert( (*mem_desc).addr==tmp);
#endif
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
#ifdef __CRAYXE
        A1DI_Release_mem_descriptor(ptr);
#endif
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
