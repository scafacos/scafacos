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

#ifdef __bgp__

DCMF_Protocol_t A1D_PutC_protocol;
DCMF_Protocol_t A1D_GetC_protocol;
DCMF_Protocol_t A1D_AccC_protocol;
DCMF_Protocol_t A1D_RemoteDone_protocol;

void A1DI_Done_cb(void * clientdata, DCMF_Error_t * error)
{
#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Done_cb() \n");
#endif

    int32_t * temp = (int32_t *) clientdata;
    (*temp) = 0; /* TODO use atomic operation */

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Done_cb() \n");
#endif
    return;
}

void A1DI_Control_done_cb(void * clientdata, const DCMF_Control_t * info, size_t peer)
{
#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Control_done_cb() \n");
#endif

    int32_t * temp = (int32_t *) clientdata;
    (*temp) = 0; /* TODO use atomic operation */

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Control_done_cb() \n");
#endif
    return;
}

void A1DI_AccC_short_cb(void *clientdata,
                        const DCQuad *msginfo,
                        unsigned count,
                        size_t peer,
                        const char *src,
                        size_t bytes)
{
    size_t i;
    DCMF_Result dcmf_result;
    DCMF_Control_t info;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_AccC_short_cb() \n");
#endif

    dcmf_result =  DCMF_Control(&A1D_RemoteDone_protocol,
                                DCMF_SEQUENTIAL_CONSISTENCY,
                                peer,
                                &info);
    assert(dcmf_result==DCMF_SUCCESS);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_AccC_short_cb() \n");
#endif
    return;
}

/*********************************************************************/

int A1DI_PutC_Initialize()
{
    DCMF_Result dcmf_result;
    DCMF_Put_Configuration_t conf;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_PutC_Initialize() \n");
#endif

    DCMF_CriticalSection_enter(0);

    conf.protocol = DCMF_DEFAULT_PUT_PROTOCOL;
    conf.network  = DCMF_DEFAULT_NETWORK;

    dcmf_result = DCMF_Put_register(&A1D_PutC_protocol, &conf);
    assert(dcmf_result==DCMF_SUCCESS);

    DCMF_CriticalSection_exit(0);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_PutC_Initialize() \n");
#endif

    return(0);
}

int A1DI_GetC_Initialize()
{
    DCMF_Result dcmf_result;
    DCMF_Get_Configuration_t conf;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_GetC_Initialize() \n");
#endif

    conf.protocol = DCMF_DEFAULT_GET_PROTOCOL;
    conf.network  = DCMF_DEFAULT_NETWORK;

    dcmf_result = DCMF_Get_register(&A1D_GetC_protocol, &conf);
    assert(dcmf_result==DCMF_SUCCESS);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_GetC_Initialize() \n");
#endif

    return(0);
}

int A1DI_AccC_Initialize()
{
    DCMF_Result dcmf_result;
    DCMF_Send_Configuration_t send_conf;
    DCMF_Control_Configuration_t control_conf;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_AccC_Initialize() \n");
#endif

    send_conf.protocol                 = DCMF_DEFAULT_SEND_PROTOCOL;
    send_conf.network                  = DCMF_DEFAULT_NETWORK;
    send_conf.cb_recv_short            = NULL;
    send_conf.cb_recv_short_clientdata = NULL;
    send_conf.cb_recv                  = NULL;
    send_conf.cb_recv_clientdata       = NULL;

    dcmf_result = DCMF_Send_register(&A1D_AccC_protocol, &send_conf);
    assert(dcmf_result==DCMF_SUCCESS);

    control_conf.protocol           = DCMF_DEFAULT_CONTROL_PROTOCOL;
    control_conf.network            = DCMF_DEFAULT_NETWORK;
    control_conf.cb_recv            = A1DI_Control_done_cb;
    control_conf.cb_recv_clientdata = NULL;

    dcmf_result = DCMF_Control_register(&A1D_RemoteDone_protocol, &control_conf);
    assert(dcmf_result==DCMF_SUCCESS);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_AccC_Initialize() \n");
#endif

    return(0);
}

/*********************************************************************/

int A1DI_Put_Initialize()
{

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Put_Initialize() \n");
#endif

    A1DI_PutC_Initialize();

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Put_Initialize() \n");
#endif

    return(0);
}

int A1DI_Get_Initialize()
{

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Get_Initialize() \n");
#endif

    A1DI_GetC_Initialize();

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Get_Initialize() \n");
#endif

    return(0);
}

int A1DI_Acc_Initialize()
{

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Acc_Initialize() \n");
#endif

    A1DI_AccC_Initialize();

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Acc_Initialize() \n");
#endif

    return(0);
}

#endif

/*********************************************************************/

#define DCMF_FLUSH_COUNT_MAX 100

/*********************************************************************/

int A1D_Flush(int target)
{
#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Flush(int target) \n");
#endif

#if defined(FLUSH_IMPLEMENTED)
    /* inject put */
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Flush(int target) \n");
#endif

    return(0);
}

int A1D_Flush_all(void)
{
    int count = 0;
    int temp[DCMF_FLUSH_COUNT_MAX+1];

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Flush(int target) \n");
#endif

#if defined(FLUSH_IMPLEMENTED)
    for ( int i=0 ; i<mpi_size ; i++)
    {
        if ( A1D_Put_flush_list[i] > 0 )
        {
            /* inject puts */

            count++;

            if ( count > DCMF_FLUSH_COUNT_MAX )
            {
                /* wait and advance */
            }
        }
    }
#endif

#ifdef FLUSH_IMPLEMENTED
    for ( int i=0 ; i<mpi_size ; i++) A1D_Put_flush_list[i] = 0;
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Flush(int target) \n");
#endif

    return(0);
}

/*********************************************************************/

int A1D_Reset(a1d_nbhandle_t * nbhandle)
{
    return(0);
}

int A1D_Test(a1d_nbhandle_t * nbhandle, int * status)
{
    return(0);
}

int A1D_Wait(a1d_nbhandle_t * nbhandle)
{
    return(0);
}

int A1D_Reset_list(int count, a1d_nbhandle_t * nbhandle)
{
    return(0);
}

int A1D_Test_list(int count, a1d_nbhandle_t * nbhandle, int * statuses)
{
    return(0);
}

int A1D_Wait_list(int count, a1d_nbhandle_t * nbhandle)
{
    return(0);
}

/*********************************************************************/

int A1D_GetC(int target, int bytes, void* src, void* dst)
{
#ifdef __bgp__
    DCMF_Result dcmf_result;
    DCMF_Request_t request;
    DCMF_Callback_t done_callback;
    volatile int32_t done_active;
    size_t src_disp, dst_disp;
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_GetC(int target, int bytes, void* src, void* dst) \n");
#endif

#ifdef __bgp__

    DCMF_CriticalSection_enter(0);

    done_callback.function = A1DI_Done_cb;
    done_callback.clientdata = (void *) &done_active;

    done_active = 1; /* TODO use atomic operation */

    src_disp = (size_t) src - (size_t) A1D_Baseptr_list[mpi_rank];
    dst_disp = (size_t) dst - (size_t) A1D_Baseptr_list[target];

    dcmf_result = DCMF_Get(&A1D_GetC_protocol,
                           &request,
                           done_callback,
                           //DCMF_RELAXED_CONSISTENCY,
                           DCMF_SEQUENTIAL_CONSISTENCY,
                           target,
                           bytes,
                           &A1D_Memregion_list[target],
                           &A1D_Memregion_list[mpi_rank],
                           src_disp,
                           dst_disp);
    assert(dcmf_result==DCMF_SUCCESS);

    A1DI_Conditional_advance(done_active > 0);

    DCMF_CriticalSection_exit(0);

#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_GetC(int target, int bytes, void* src, void* dst) \n");
#endif

    return(0);
}

int A1D_PutC(int target, int bytes, void* src, void* dst)
{
#ifdef __bgp__
    DCMF_Result dcmf_result;
    DCMF_Request_t request;
    DCMF_Callback_t done_callback;
    volatile int done_active;
    size_t src_disp, dst_disp;
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_PutC(int target, int bytes, void* src, void* dst) \n");
#endif

#ifdef __bgp__
    DCMF_CriticalSection_enter(0);

    src_disp = (size_t) src - (size_t) A1D_Baseptr_list[mpi_rank];
    dst_disp = (size_t) dst - (size_t) A1D_Baseptr_list[target];

    done_callback.function = A1DI_Done_cb;
    done_callback.clientdata = (void *) &done_active;

    done_active = 1; /* TODO use atomic operation */

#ifdef FLUSH_IMPLEMENTED

    A1D_Put_flush_list[target]++;

    /* local completion only - must flush later */
    dcmf_result = DCMF_Put(&A1D_PutC_protocol,
                           &request,
                           done_callback, /* local completion */
                           DCMF_SEQUENTIAL_CONSISTENCY,
                           target,
                           bytes,
                           &A1D_Memregion_list[mpi_rank],
                           &A1D_Memregion_list[target],
                           src_disp,
                           dst_disp,
                           A1D_Nocallback); /* remote completion */
    assert(dcmf_result==DCMF_SUCCESS);

    A1DI_Conditional_advance(done_active > 0);

    DCMF_CriticalSection_exit(0);

#else

    DCMF_CriticalSection_enter(0);

    /* end-to-end completion - no flush required */
    dcmf_result = DCMF_Put(&A1D_PutC_protocol,
                           &request,
                           A1D_Nocallback, /* local completion */
                           DCMF_RELAXED_CONSISTENCY,
                           target,
                           bytes,
                           &A1D_Memregion_list[mpi_rank],
                           &A1D_Memregion_list[target],
                           src_disp,
                           dst_disp,
                           done_callback); /* remote completion */
    assert(dcmf_result==DCMF_SUCCESS);

#endif

    A1DI_Conditional_advance(done_active > 0);

    DCMF_CriticalSection_exit(0);
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_PutC(int target, int bytes, void* src, void* dst) \n");
#endif

    return(0);
}

int A1D_AccC(int proc, int bytes, void* src, void* dst, int type, void* scale)
{
    return 0;
}

/*********************************************************************/

#ifdef STRIDED_IMPLEMENTED

int A1D_GetS(int proc, stride_levels, block_sizes,
             src_ptr, src_stride_arr,
             dst_ptr, dst_stride_arr)
{
    return 0;
}

int A1D_PutS(int proc, stride_levels, block_sizes,
             src_ptr, src_stride_arr,
             dst_ptr, dst_stride_arr)
{
    return 0;
}

int A1D_AccS(int proc, stride_levels, block_sizes,
             src_ptr, src_stride_arr,
             dst_ptr, dst_stride_arr,
             int type, void* scale)
{
    return 0;
}

#endif

/*********************************************************************/
