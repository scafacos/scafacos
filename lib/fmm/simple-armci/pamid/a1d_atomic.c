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
 * int64_tellectual property rights of third parties.  The copyright holders
 * disclaim any liability to any recipient for claims brought against
 * recipient by any third party for infringement of that parties
 * int64_tellectual property rights.
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

DCMF_Protocol_t A1D_Fetch64_protocol;
DCMF_Protocol_t A1D_Inc64_protocol;

/***********************************************************************/

void A1DI_Fetch64_cb(void * clientdata, const DCMF_Control_t * info, size_t peer)
{
    int64_t   value = 0;
    int64_t * return_address = NULL;
    volatile uint64_t * active_address = NULL;
    A1D_Fetch64_t data;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Fetch64_cb \n");
#endif

    memcpy( &data, info, sizeof(A1D_Fetch64_t) );

    value          = data.value;
    return_address = data.return_address;
    active_address = data.active_address;

    fprintf(stderr,"A1D_Fetch64_cb A rank = %d peer = %d, value = %d, return_address = %p, active_address = %p, *active_address = %u \n", A1D_Rank(), peer, value, return_address, active_address, *active_address );

    if ( return_address == NULL )
    {
        fprintf(stderr,"A1DI_Fetch64_cb: return_address is a NULL pointer. This is bad. \n");
        assert( return_address != NULL );
    }
    else
    {
        (*return_address) = value;
    }

    if ( active_address == NULL )
    {
        fprintf(stderr,"A1DI_Fetch64_cb: active_address is a NULL pointer. This is bad. \n");
        assert( active_address != NULL );
    }
    else
    {
        fprintf(stderr,"A1D_Fetch64_cb B rank = %d peer = %d, value = %d, return_address = %p, active_address = %p, *active_address = %u \n", A1D_Rank(), peer, value, return_address, active_address, *active_address );

        (*active_address) = 0;

        fprintf(stderr,"A1D_Fetch64_cb C rank = %d peer = %d, value = %d, return_address = %p, active_address = %p, *active_address = %u \n", A1D_Rank(), peer, value, return_address, active_address, *active_address );
    }

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Fetch64_cb \n");
#endif

    return;
}

void A1DI_Inc64_cb(void * clientdata, const DCMF_Control_t * info, size_t peer)
{
    int64_t   incr = 999999999;
    int64_t * incr_address = NULL;
    int64_t * return_address = NULL;
    volatile uint64_t * active_address = NULL;
    A1D_Inc64_t data;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Inc64_cb \n");
#endif

    memcpy( &data, info, sizeof(A1D_Inc64_t) );
    incr           = data.incr;
    incr_address   = data.incr_address;
    return_address = data.return_address;
    active_address = data.active_address;

    fprintf(stderr,"A1D_Inc64_cb rank = %d peer = %d, incr = %d, incr_address = %p return_address = %p active_address = %p \n", A1D_Rank(), peer, incr, incr_address, return_address, active_address );

    if ( incr_address == NULL )
    {
        fprintf(stderr,"A1DI_Inc64_cb: incr_address is a NULL pointer. This is bad. \n");
        assert( incr_address != NULL );
    }
    else
    {
        if ( return_address != NULL )
        {
            /* if sending message back to source, make sure the active address is valid */
            if ( active_address == NULL )
            {
                fprintf(stderr,"A1DI_Inc64_cb: active_address is a NULL pointer. This is bad. \n");
                assert( active_address != NULL );
            }

            DCMF_Result dcmf_result;
            A1D_Fetch64_t return_data;
            DCMF_Control_t return_payload;

            register int64_t old_val, new_val;
            old_val = (*incr_address);
            new_val = old_val + incr; /* TODO use atomic operation */
            (*incr_address) = new_val;

            return_data.value          = old_val;
            return_data.return_address = return_address;
            return_data.active_address = active_address;

            memcpy(&return_payload, &return_data, sizeof(A1D_Fetch64_t));

            dcmf_result =  DCMF_Control(&A1D_Fetch64_protocol,
                                        DCMF_SEQUENTIAL_CONSISTENCY,
                                        peer,
                                        &return_payload);
            assert(dcmf_result==DCMF_SUCCESS);
        }
        /* else return_address == NULL do not send anything back to source  */
    }

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Inc64_cb \n");
#endif

    return;
}

/***********************************************************************/

void A1DI_Fetch64_Initialize()
{
    DCMF_Result dcmf_result;
    DCMF_Control_Configuration_t conf;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Fetch64_Initialize \n");
#endif

    if ( sizeof(A1D_Fetch64_t) > sizeof(DCMF_Control_t) )
    {
        fprintf(stderr,"A1D_Fetch64_t requires more storage than DCMF_Control_t! \n");
        assert( sizeof(A1D_Fetch64_t) == sizeof(DCMF_Control_t) );
    }

    conf.protocol           = DCMF_DEFAULT_CONTROL_PROTOCOL;
    conf.network            = DCMF_DEFAULT_NETWORK;
    conf.cb_recv            = A1DI_Fetch64_cb;
    conf.cb_recv_clientdata = NULL;

    dcmf_result = DCMF_Control_register(&A1D_Fetch64_protocol, &conf);
    assert(dcmf_result==DCMF_SUCCESS);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Fetch64_Initialize \n");
#endif

    return;
}

void A1DI_Inc64_Initialize()
{
    DCMF_Result dcmf_result;
    DCMF_Control_Configuration_t conf;

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Inc64_Initialize \n");
#endif

    if ( sizeof(A1D_Inc64_t) > sizeof(DCMF_Control_t) )
    {
        fprintf(stderr,"A1D_Inc64_t requires more storage than DCMF_Control_t! \n");
        assert( sizeof(A1D_Inc64_t) == sizeof(DCMF_Control_t) );
    }

    conf.protocol           = DCMF_DEFAULT_CONTROL_PROTOCOL;
    conf.network            = DCMF_DEFAULT_NETWORK;
    conf.cb_recv            = A1DI_Inc64_cb;
    conf.cb_recv_clientdata = NULL;

    dcmf_result = DCMF_Control_register(&A1D_Inc64_protocol, &conf);
    assert(dcmf_result==DCMF_SUCCESS);

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Inc64_Initialize \n");
#endif

    return;
}

void A1DI_Atomic_Initialize()
{
#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1DI_Atomic_Initialize \n");
#endif

    A1DI_Fetch64_Initialize();
    A1DI_Inc64_Initialize();

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1DI_Atomic_Initialize \n");
#endif

    return;
}

/* ifdef __bgp__ */

#endif

/***********************************************************************/

void A1D_Fetch_and_inc64(int proc, int64_t * local, int64_t * remote, int64_t incr)
{
#ifdef __bgp__
    DCMF_Result dcmf_result;
    A1D_Inc64_t data;
    DCMF_Control_t payload;
    volatile uint64_t active = 0;
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"entering A1D_Fetch_and_inc64 \n");
#endif

    if ( ((unsigned long)remote & 0x07uL != 0) || ((unsigned long)local & 0x07uL != 0) )
    {
        printf("%d: remote = %p local = %p \n", A1D_Rank(), remote, local );
        assert( ((unsigned long)remote & 0x07uL != 0) && ((unsigned long)local & 0x07uL != 0) );
    }

#ifdef __bgp__
    DCMF_CriticalSection_enter(0);

    data.incr           = incr;
    data.incr_address   = remote;
    data.return_address = local;
    data.active_address = &active;

    memcpy(&payload, &data, sizeof(A1D_Inc64_t));

    active = 1;

    dcmf_result = DCMF_Control(&A1D_Inc64_protocol,
                               DCMF_SEQUENTIAL_CONSISTENCY,
                               proc,
                               &payload);
    assert(dcmf_result==DCMF_SUCCESS);

    A1DI_Conditional_advance(active);

    DCMF_CriticalSection_exit(0);
#endif

#ifdef DEBUG_FUNCTION_ENTER_EXIT
    fprintf(stderr,"exiting A1D_Fetch_and_inc64 \n");
#endif
    return;
}
