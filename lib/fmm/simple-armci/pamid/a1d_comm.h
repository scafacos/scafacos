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

#ifndef A1D_COMM_H
#define A1D_COMM_H

typedef enum
{
    A1D_DOUBLE,
    A1D_SINGLE,
#ifdef A1D_USE_COMPLEX
    A1D_DOUBLE_COMPLEX,
    A1D_SINGLE_COMPLEX,
#endif
    A1D_INT32,
    A1D_UINT32,
    A1D_INT64,
    A1D_UINT64
}
a1d_datatype_t;

typedef struct
{
    uint32_t aggr_size; /* 0 means not aggregate handle */
    union
    {
        int   handle;
        int * handles;
    }
    nbh;
}
a1d_nbhandle_t;

int A1DI_Put_Initialize();
int A1DI_Get_Initialize();
int A1DI_Acc_Initialize();

int A1D_Flush(int target);
int A1D_Flush_all(void);

int A1D_Reset(a1d_nbhandle_t * nbhandle);
int A1D_Test(a1d_nbhandle_t * nbhandle, int * status);
int A1D_Wait(a1d_nbhandle_t * nbhandle);

int A1D_Reset_list(int count, a1d_nbhandle_t * nbhandle);
int A1D_Test_list(int count, a1d_nbhandle_t * nbhandle, int * statuses);
int A1D_Wait_list(int count, a1d_nbhandle_t * nbhandle);

int A1D_GetC(int proc, int bytes, void* src, void* dst);
int A1D_PutC(int proc, int bytes, void* src, void* dst);
int A1D_AccC(int proc, int bytes, void* src, void* dst, int type, void* scale);

int A1D_iGetC(int proc, int bytes, void * src, void * dst, a1d_nbhandle_t * nbhandle);
int A1D_iPutC(int proc, int bytes, void * src, void * dst, a1d_nbhandle_t * nbhandle);
int A1D_iAccC(int proc, int bytes, void * src, void * dst, a1d_nbhandle_t * nbhandle, int type, void * scale);

/*
int A1D_GetS(int proc, int stride_levels, int block_sizes,
                          src_ptr, src_stride_arr,
                          dst_ptr, dst_stride_arr);
int A1D_PutS(int proc, int stride_levels, int block_sizes,
                          src_ptr, src_stride_arr,
                          dst_ptr, dst_stride_arr);
int A1D_AccS(int proc, stride_levels, block_sizes,
                          src_ptr, src_stride_arr,
                          dst_ptr, dst_stride_arr,
                          int type, void* scale);
*/

#endif
