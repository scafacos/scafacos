/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "a1.h"
#include "a1u.h"

#if !defined A1D_H_INCLUDED
#define A1D_H_INCLUDED

/* ********************************************************************* */
/*                                                                       */
/*               A1 data structures                                      */
/*                                                                       */
/* ********************************************************************* */

/* FIXME: This should be in another header but we're using a1u.h
 *         already exists for other purposes.  I will deal with all of
 *         this after we submit the IPDPS paper and fix the build system.
 */

/**
 * \brief A1 device-independent settings type.
 *
 * \see A1_Initialize
 *
 * \ingroup TYPEDEFS
 *
 */

typedef struct
{
    uint32_t network_bypass;
    uint32_t network_bypass_upper_limit_1d;
    uint32_t network_bypass_upper_limit_Nd;
    uint32_t armci_strict_ordering;
} A1U_Settings_t;

extern A1U_Settings_t a1u_settings;

/* ********************************************************************* */
/*                                                                       */
/*               A1 device-independent internal functions                */
/*                                                                       */
/* ********************************************************************* */

/**
 * \brief Parses the device-indpendent parameter information
 *
 * \param[out] rc               The error code
 *
 * \ingroup MANAGEMENT
 */

int A1U_Read_parameters(void);

/**
 * \brief Prints the device-indpendent parameter information
 *
 * \param[out] rc               The error code
 *
 * \ingroup MANAGEMENT
 */

int A1U_Print_parameters(void);

/* ********************************************************************* */
/*                                                                       */
/*               A1 device-level functions                               */
/*                                                                       */
/* ********************************************************************* */

/**
 * \brief Device level implementation of A1_Initialize.
 *
 * \param[out] rc               The error code from initializing A1
 * \param[in]  A1_thread_level  The type of thread support for A1
 *
 * \ingroup MANAGEMENT
 */
int A1D_Initialize(int A1_thread_level);

/**
 * \brief Device level implementation of A1_Finalize.
 *
 * \param[out] rc  The error code from terminating A1.  
 *
 * \ingroup MANAGEMENT
 */
int A1D_Finalize(void);

/**
 * \brief Device level implementation of A1_Abort.
 *
 * \param[in]  error_code    The error code to be returned to the submitting environment.
 * \param[in]  error_message Text string to print to stderr upon termination.
 *
 * \ingroup MANAGEMENT
 */

void A1D_Abort(int error_code, char error_message[]);

/**
 * \brief Device level implementation of A1D_Alloc_segment.
 *
 * A local operation to allocate memory to be used in context of A1 copy operations.
 *
 * \note Memory allocated with this function will be properly aligned for the architecture.
 *
 * \warning Memory allocated with this function must be freed by A1_Free_segment.
 *
 * \param[out] rc            The error code.
 * \param[out] ptr           Pointer to memory.
 * \param[in]  bytes         The amount of memory requested.
 *
 * \ingroup MEMORY
 */
int A1D_Alloc_segment(void** pointer, int bytes);

/**
 * \brief Device level implementation of A1D_Free_segment.
 *
 * A local operation to free memory allocated by A1_Alloc_segment.
 *
 * \warning It is erroneous to attempt to free memory not allocated by A1_Alloc_segment.
 *
 * \param[out] rc            The error code.
 * \param[in] ptr           Pointer to memory.
 *
 * \ingroup MEMORY
 */

int A1D_Free_segment(void* pointer);

/**
 * \brief Device level implementation of A1_Exchange_segments.
 *
 *  A collective operation to allocate memory to be used in context of A1 copy operations.
 *
 * \param[out] rc         The error code.
 * \param[in]  group      Group of processes within which the pointer list is exchanged.
 * \param[in]  ptr        Pointer array. Each one points to memory allocated at one process, 
 *                        in order of ranks.
 * \param[in]  bytes      The size of memory allocated at each process.
 *
 * \ingroup MEMORY 
 */
int A1D_Exchange_segments(A1_group_t* group, void **ptr);

/**
 * \brief Device level implementation of A1_Release_segments.
 *
 * A collective operation to invalidate and de-register memory segments
 * associated with an A1D_Exchange_segments call. 
 *
 * \param[out] rc          The error code.
 * \param[in]  group       Group of processes within which the pointer list was exchanged.
 * \param[in]  ptr         Pointer to the allocated memory.
 *
 * \ingroup MEMORY 
 */
int A1D_Release_segments(A1_group_t* group, void *ptr);

/**
 * \brief Device level implementation of A1_Allocate_handle.
 *
 * Allocates a non-blocking handle.
 *
 * \param[in] handle      Non-blocking handle upon which to be waited.
 *
 * \see A1_handle_t, A1_Wait_handle_list, A1_Test_handle
 *
 * \ingroup MEMORY
 */

int A1D_Allocate_handle(A1_handle_t *handle);


/**
 * \brief Device level implementation of A1_Release_handle.
 * 
 * Releases a non-blocking handle.
 * 
 * \param[in] handle      Non-blocking handle upon which to be waited.
 *
 * \see A1_handle_t, A1_Wait_handle_list, A1_Test_handle
 *
 * \ingroup MEMORY
 */

int A1D_Release_handle(A1_handle_t handle);

/**
 * \brief Device level implementation of A1_Wait_handle_all.
 *
 * Waits for operations on all handle to complete.
 *
 * \param[in] handle      Non-blocking handle upon which to be waited.
 *
 * \see A1_handle_t, A1_Wait_handle_list, A1_Test_handle
 *
 * \ingroup MEMORY
 */

int A1D_Wait_handle_all(void);

/**
 * \brief Device level implementation of A1_Wait_handle.
 *
 * Waits for operations on a handle to complete.
 *
 * \param[in] handle      Non-blocking handle upon which to be waited.
 *
 * \see A1_handle_t, A1_Wait_handle_list, A1_Test_handle
 *
 * \ingroup MEMORY
 */

int A1D_Wait_handle(A1_handle_t handle);

/**
 * \brief Device level implementation of A1_Wait_handle_list.
 *
 * Waits for operations on a list of handles to complete.
 *
 * \param[in] count          Number of handles
 * \param[in] a1_handle      Non-blocking handles upon which to be waited.
 *
 * \see A1_handle_t, A1_Wait_handle_list, A1_Test_handle
 *
 * \ingroup MEMORY
 */

int A1D_Wait_handle_list(int count, A1_handle_t *a1_handle);

/**
 * \brief Device level implementation of A1_Test_handle.
 *
 * Test for completion of operations on a handle.
 *
 * \param[in] handle      Non-blocking handle upon which to be waited.
 *
 * \see A1_handle_t, A1_Wait_handle_list, A1_Test_handle
 *
 * \ingroup MEMORY
 */

int A1D_Test_handle(A1_handle_t handle, A1_bool_t* completed);

/**
 * \brief Device level implementation of A1_Test_handle_list.
 *
 * Test for completion of operations on a list of handles.
 *
 * \param[in] count          Number of handles
 * \param[in] a1_handle      Non-blocking handles upon which to be tested.
 *
 * \see A1_handle_t, A1_Wait_handle_list, A1_Test_handle
 *
 * \ingroup MEMORY
 */
int A1D_Test_handle_list(int count,
                         A1_handle_t *a1_handle,
                         A1_bool_t* *completed);

/**
 * \brief Device level implementation of A1_Barrier_group.
 *
 * On return, this call ensures that all processes within the entire group
 * have reached this point in the program.
 *
 * \param[in] group          Group of processes to synchronize.
 *
 * \ingroup  SYNCHRONIZATION
 */
int A1D_Barrier_group(A1_group_t* group);

/**
 * \brief Device level implementation of A1_NbBarrier_group.
 *
 * \param[in] group          Group of processes to synchronize.
 *
 * \ingroup  SYNCHRONIZATION
 */
int A1_NbBarrier_group(A1_group_t* group, A1_handle_t handle);

/**
 * \brief Device level implementation of A1_Sync_group.
 *
 * On return, this call ensures that all processes within the entire group
 * have reached this point in the program and that all messages have completed remotely.
 *
 * \param[in] group          Group of processes to synchronize.
 *
 * \ingroup  SYNCHRONIZATION
 */
int A1D_Sync_group(A1_group_t* group);

/**
 * \brief Device level implementation of A1_NbSync_group.
 *
 * \param[in] group          Group of processes to synchronize.
 *
 * \ingroup  SYNCHRONIZATION
 */
int A1D_NbSync_group(A1_group_t* group, A1_handle_t handle);

/**
 * \brief Device level implementation of A1_Put.
 *
 * Blocking copy of contiguous data from local memory to remote memory.
 *
 * \param[out] rc            The error code.
 * \param[in]  target        Rank of the remote process.
 * \param[in]  source_ptr    Starting address in the (local) source memory.
 * \param[in]  target_ptr    Starting address in the (remote) target memory.
 * \param[in]  bytes         Amount of data to transfer, in bytes. 
 *
 * \ingroup  COPY OPERATIONS
 */
int A1D_Put(int target, void* src, void* dst, int bytes);

/**
 * \brief Device level implementation of A1_NbPut.
 *
 * Non-Blocking copy of contiguous data from local memory to remote memory.
 *
 * \param[out] rc            The error code.
 * \param[in]  target        Rank of the remote process.
 * \param[in]  source_ptr    Starting address in the (local) source memory.
 * \param[in]  target_ptr    Starting address in the (remote) target memory.
 * \param[in]  bytes         Amount of data to transfer, in bytes.
 * \param[in]  handle        Opaque A1 handle for request
 *
 * \ingroup  COPY OPERATIONS
 */
int A1D_NbPut(int target, void* src, void* dst, int bytes, A1_handle_t handle);

/**
 * \brief Device level implementation of A1_PutS.
 *
 * Blocking copy of non-contiguous (strided) data from local memory to remote memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  stride_level   The number of levels of stride.
 * \param[in]  block_sizes     Block size in each dimension, in bytes.
 * \param[in]  source_ptr      Starting address in the (local) source memory.
 * \param[in]  src_stride_ar   Array of stride distances at source, in bytes.
 * \param[in]  target_ptr      Starting address in the (remote) target memory.
 * \param[in]  trg_stride_ar   Array of stride distances at target, in bytes.
 *
 * \ingroup COPY OPERATIONS
 */
int A1D_PutS(int target,
             int stride_level,
             int block_sizes[],
             void* source_ptr,
             int src_stride_ar[],
             void* target_ptr,
             int trg_stride_ar[]);

/**
 * \brief Device level implementation of A1_NbPutS.
 * 
 * Non-Blocking copy of non-contiguous (strided) data from local memory to remote memory.
 * 
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  stride_level   The number of levels of stride.
 * \param[in]  block_sizes     Block size in each dimension, in bytes.
 * \param[in]  source_ptr      Starting address in the (local) source memory.
 * \param[in]  src_stride_ar   Array of stride distances at source, in bytes.
 * \param[in]  target_ptr      Starting address in the (remote) target memory.
 * \param[in]  trg_stride_ar   Array of stride distances at target, in bytes.
 * \param[in]  handle        Opaque A1 handle for request
 *
 * \ingroup COPY OPERATIONS
 */

int A1D_NbPutS(int target,
               int stride_level,
               int block_sizes[],
               void* source_ptr,
               int src_stride_ar[],
               void* target_ptr,
               int trg_stride_ar[],
               A1_handle_t handle);

/**
 * \brief Device level implementation of A1_PutV.
 * 
 *  Blocking copy of non-contiguous data from local memory to remote memory.
 * 
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  iov_ar          Array of io vectors. Each vector represents a set of
 *                             chunks of same size.
 * \param[in]  ar_len          Number of elements in the array.
 * 
 * \see A1_NbPut, A1_NbPutV, A1_NbMultiPut, A1_NbMultiPutS, A1_NbMultiPutV
 * 
 * \ingroup DATA_TRANSFER
 */

int A1D_PutV(int target, 
             A1_iov_t *iov_ar,
             int ar_len);

/**
 * \brief Device level implementation of A1_NbPutV.
 * 
 *  Non-Blocking copy of non-contiguous data from local memory to remote memory.
 * 
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  iov_ar          Array of io vectors. Each vector represents a set of
 *                             chunks of same size.
 * \param[in]  ar_len          Number of elements in the array.
 * \param[in]  a1_handle       A1 Opaque handle
 * 
 * \see A1_NbPut, A1_NbPutV, A1_NbMultiPut, A1_NbMultiPutS, A1_NbMultiPutV
 * 
 * \ingroup DATA_TRANSFER
 */

int A1D_NbPutV(int target, 
               A1_iov_t *iov_ar,
               int ar_len,
               A1_handle_t a1_handle);

/**
 * \brief Device level implementation of A1_Get.
 *
 * Blocking copy of contiguous data from remote memory to local memory.
 *
 * \param[out] rc            The error code.
 * \param[in]  target        Rank of the remote process.
 * \param[in]  source_ptr    Starting address in the (remote) source memory.
 * \param[in]  target_ptr    Starting address in the (local) target memory.
 * \param[in]  bytes         Amount of data to transfer, in bytes.
 *
 * \ingroup  COPY OPERATIONS
 */
int A1D_Get(int target, void* src, void* dst, int bytes);

/**
 * \brief Device level implementation of A1_Get.
 *
 * Non-Blocking copy of contiguous data from remote memory to local memory.
 *
 * \param[out] rc            The error code.
 * \param[in]  target        Rank of the remote process.
 * \param[in]  source_ptr    Starting address in the (remote) source memory.
 * \param[in]  target_ptr    Starting address in the (local) target memory.
 * \param[in]  bytes         Amount of data to transfer, in bytes.
 * \param[in]  handle        Opaque A1 handle for request
 *
 * \ingroup  COPY OPERATIONS
 */
int A1D_NbGet(int target, void* src, void* dst, int bytes, A1_handle_t handle);

/**
 * \brief Device level implementation of A1_GetS.
 *
 * Blocking copy of non-contiguous (strided) data from remote memory to local memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  stride_level   The number of levels of stride.
 * \param[in]  block_sizes     Block size in each dimension, in bytes.
 * \param[in]  source_ptr      Starting address in the (remote) source memory.
 * \param[in]  src_stride_ar   Array of stride distances at source, in bytes.
 * \param[in]  target_ptr      Starting address in the (local) target memory.
 * \param[in]  trg_stride_ar   Array of stride distances at target, in bytes.
 *
 * \ingroup COPY OPERATIONS
 */
int A1D_GetS(int target,
             int stride_level,
             int block_sizes[],
             void* source_ptr,
             int src_stride_ar[],
             void* target_ptr,
             int trg_stride_ar[]);

/**
 * \brief Device level implementation of A1_NbGetS.
 *
 * Non-Blocking copy of non-contiguous (strided) data from remote memory to local memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  stride_level    The number of levels of stride.
 * \param[in]  block_sizes     Block size in each dimension, in bytes.
 * \param[in]  source_ptr      Starting address in the (remote) source memory.
 * \param[in]  src_stride_ar   Array of stride distances at source, in bytes.
 * \param[in]  target_ptr      Starting address in the (local) target memory.
 * \param[in]  trg_stride_ar   Array of stride distances at target, in bytes.
 * \param[in]  handle          Opaque A1 handle for request
 *
 * \ingroup COPY OPERATIONS
 */
int A1D_NbGetS(int target,
               int stride_level,
               int block_sizes[],
               void* source_ptr,
               int src_stride_ar[],
               void* target_ptr,
               int trg_stride_ar[],
               A1_handle_t handle);

/**
 * \brief  Device level implementation of A1_GetV.
 *
 * Blocking copy of non-contiguous data from remote memory to local memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  iov_ar          Array of io vectors. Each vector represents a set of
 *                             chunks of same size.
 * \param[in]  ar_len          Number of elements in the array.
 * 
 * \see A1_NbPut, A1_NbPutV, A1_NbMultiPut, A1_NbMultiPutS, A1_NbMultiPutV
 * 
 * \ingroup DATA_TRANSFER
 */

int A1D_GetV(int target,
             A1_iov_t *iov_ar,
             int ar_len);

/**
 * \brief Device level implementation of A1_NbGetV
 *
 * Non-Blocking copy of non-contiguous data from remote memory to local memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  iov_ar          Array of io vectors. Each vector represents a set of
 *                             chunks of same size.
 * \param[in]  ar_len          Number of elements in the array.
 * \param[in]  handle          A1 Opaque handle
 * 
 * \see A1_NbPut, A1_NbPutV, A1_NbMultiPut, A1_NbMultiPutS, A1_NbMultiPutV
 * 
 * \ingroup DATA_TRANSFER
 */

int A1D_NbGetV(int target,
               A1_iov_t *iov_ar,
               int ar_len,
               A1_handle_t handle);

/**
 * \brief Device level implementation of A1_PutAcc
 *
 * Blocking accumulate of contiguous data from local memory onto remote memory.
 *
 * \param[out] rc            The error code.
 * \param[in]  target        Rank of the remote process.
 * \param[in]  source_ptr    Starting address in the (local) source memory.
 * \param[in]  target_ptr    Starting address in the (remote) target memory.
 * \param[in]  bytes         Amount of data to transfer, in bytes.
 * \param[in]  a1_type       Amount of data to transfer, in bytes.
 * \param[in]  scaling       Factor for scaling source
 *
 * \ingroup COPY OPERATIONS
 */
int A1D_PutAcc(int target,
               void* source_ptr,
               void* target_ptr,
               int bytes,
               A1_datatype_t a1_type,
               void* scaling);

/**
 * \brief Device level implementation of A1_NbPutAcc
 * 
 * Non-Blocking accumulate of contiguous data from local memory onto remote memory.
 * 
 * \param[out] rc            The error code.
 * \param[in]  target        Rank of the remote process.
 * \param[in]  source_ptr    Starting address in the (local) source memory.
 * \param[in]  target_ptr    Starting address in the (remote) target memory.
 * \param[in]  bytes         Amount of data to transfer, in bytes.
 * \param[in]  a1_type       Amount of data to transfer, in bytes.
 * \param[in]  scaling       Factor for scaling source
 * \param[in]  handle        Opaque A1 handle
 *
 * \ingroup COPY OPERATIONS
 */
int A1D_NbPutAcc(int target,
                 void* source_ptr,
                 void* target_ptr,
                 int bytes,
                 A1_datatype_t a1_type,
                 void* scaling,
                 A1_handle_t handle);

/**
 * \brief Device level implementation of A1_PutAccS 
 *
 * Blocking accumulate of non-contiguous (strided) data from local memory to remote memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  stride_level   The number of levels of stride.
 * \param[in]  block_sizes     Block size in each dimension, in bytes.
 * \param[in]  source_ptr      Starting address in the (local) source memory.
 * \param[in]  src_stride_ar   Array of stride distances at source, in bytes.
 * \param[in]  target_ptr      Starting address in the (remote) target memory.
 * \param[in]  trg_stride_ar   Array of stride distances at target, in bytes.
 * \param[in]  a1_type         Amount of data to transfer, in bytes.
 * \param[in]  scaling         Factor for scaling source
 *
 * \ingroup COPY OPERATIONS
 */
int A1D_PutAccS(int target,
                int stride_level,
                int block_sizes[],
                void* source_ptr,
                int *src_stride_ar,
                void* target_ptr,
                int *trg_stride_ar,
                A1_datatype_t a1_type,
                void* scaling);

/**
 * \brief Device level implementation of A1_NbPutAccS
 *
 * Non-Blocking accumulate of non-contiguous (strided) data from local memory to remote memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  stride_level   The number of levels of stride.
 * \param[in]  block_sizes     Block size in each dimension, in bytes.
 * \param[in]  source_ptr      Starting address in the (local) source memory.
 * \param[in]  src_stride_ar   Array of stride distances at source, in bytes.
 * \param[in]  target_ptr      Starting address in the (remote) target memory.
 * \param[in]  trg_stride_ar   Array of stride distances at target, in bytes.
 * \param[in]  a1_type         Amount of data to transfer, in bytes.
 * \param[in]  scaling         Factor for scaling source
 * \param[in]  handle          Opaque A1 handle
 *
 * \ingroup COPY OPERATIONS
 */
int A1D_NbPutAccS(int target,
                  int stride_level,
                  int block_sizes[],
                  void* source_ptr,
                  int *src_stride_ar,
                  void* target_ptr,
                  int *trg_stride_ar,
                  A1_datatype_t a1_type,
                  void* scaling,
                  A1_handle_t handle);

/**
 * \brief Device level implementation of A1_PutAccV
 * 
 * Blocking accumulate of non-contiguous data from local memory to remote memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  iov_ar          Array of io vectors. Each vector represents a set of
 *                             chunks of same size.
 * \param[in]  ar_len          Number of elements in the array.
 * \param[in]  a1_type         Type of data and scaling factor
 * \param[in]  *scaling        Scaling factor in the accumulate operation.
 *
 * \see A1_NbPut, A1_NbPutV, A1_NbMultiPut, A1_NbMultiPutS, A1_NbMultiPutV
 *
 * \ingroup DATA_TRANSFER
 */

int A1D_PutAccV(int target,
                A1_iov_t *iov_ar,
                int ar_len,
                A1_datatype_t a1_type,
                void* scaling);

/**
 * \brief Device level implementation of A1_PutAccV
 * 
 * Blocking accumulate of non-contiguous data from local memory to remote memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  iov_ar          Array of io vectors. Each vector represents a set of
 *                             chunks of same size.
 * \param[in]  ar_len          Number of elements in the array.
 * \param[in]  a1_type         Type of data and scaling factor
 * \param[in]  *scaling        Scaling factor in the accumulate operation.
 * \param[in]  a1_handle       A1 opaque handle
 *
 * \see A1_NbPut, A1_NbPutV, A1_NbMultiPut, A1_NbMultiPutS, A1_NbMultiPutV
 *
 * \ingroup DATA_TRANSFER
 */

int A1D_NbPutAccV(int target,
                  A1_iov_t *iov_ar,
                  int ar_len,
                  A1_datatype_t a1_type,
                  void* scaling,
                  A1_handle_t a1_handle);

/**
 * \brief Device level implementation of A1_PutModV
 *
 * \brief Blocking remote modify of non-contiguous data from local memory to remote memory.
 *
 * \param[out] rc              The error code.
 * \param[in]  target          Rank of the remote process.
 * \param[in]  iov_ar          Array of io vectors. Each vector represents a set of
 *                             chunks of same size.
 * \param[in]  ar_len          Number of elements in the array.
 * \param[in]  a1_op           Reduce operation
 * \param[in]  a1_type         Type of data and scaling factor
 *
 * \see A1_NbPut, A1_NbPutV, A1_NbMultiPut, A1_NbMultiPutS, A1_NbMultiPutV
 *
 * \ingroup DATA_TRANSFER
 */

int A1D_PutModV(int target,
                A1_iov_t *iov_ar,
                int ar_len,
                A1_reduce_op_t a1_op,
                A1_datatype_t a1_type);

/**
 * \brief Collective operation to allocate a counter.
 *
 * \param[out] rc            The error code.
 * \param[in]  counter       A1 shared counter.
 * 
 * \ingroup Atomics
 */
int A1D_Create_counter(A1_group_t* group, 
                       A1_counter_t *counter);

/**
 * \brief Collective operation to deallocate and deregister a counter.
 *
 * \param[out]    rc            The error code.
 * \param[in]     group         A1 group over which the counter is shared.
 * \param[inout]  counter       A1 shared counter.
 *
 * \see a1_counter_t, A1_Create_counter, A1_Incr_counter, A1_NbIncr_counter
 *
 * \ingroup Atomics
 */

int A1D_Destroy_counter(A1_group_t* group,
                       A1_counter_t *counter); 


/**
 * \brief Atomically updates a shared counter and returns the current value.
 *
 * \param[out] rc            The error code.
 * \param[in]  counter       A1 shared counter.
 * \param[in]  increment     The value to add to the counter.
 * \param[in]  original      The remote value of the counter prior to the increment.
 *
 * \see a1_counter_t, A1_Create_counter, A1_Destroy_counter, A1_NbIncr_counter
 *
 * \ingroup Atomics
 */

int A1D_Incr_counter(A1_counter_t counter,
                     long increment,
                     long* original);

/**
 * \brief Atomically updates a shared counter and returns the current value.
 *
 * \param[out] rc            The error code.
 * \param[in]  target        Rank of the target process.
 * \param[in]  source_ptr    Pointer of variable at source process.
 * \param[in]  target_ptr    Pointer of variable at target process.
 * \param[in]  op            Operation to be performed.
 * \param[in]  value         Local buffer containing the value and which will contain
 *                           the current value after the operation.
 *
 * \ingroup Atomics
 */
int A1D_Rmw(int target,
           void* source_ptr_in,
           void* source_ptr_out,
           void* target_ptr,
           int bytes,
           A1_atomic_op_t op,
           A1_datatype_t a1_type);


/**
 * \brief Device level implementation of A1_Create_mutexes
 *
 * Collective operation to allocate and register a list of mutexes.
 *
 * \param[out]    rc            The error code.
 * \param[in]     group         A1 group over which the mutexes are shared.
 * \param[in]     count         Number of mutexes to be created.
 * \param[int]    count_ar      An arrays storing the number of mutexes on each process
 *
 * \ingroup Atomics
 */
int A1D_Create_mutexes(A1_group_t* group, 
                       int mutex_count, 
                       int *mutex_count_ar);

/**
 * \brief Device level implementation of A1_Destroy_mutexes
 *
 * Collective operation to unregister and deallocate a list of mutexes.
 *
 * \param[out]    rc            The error code.
 * \param[in]     group         A1 group over which the mutexes are shared.
 *
 * \ingroup Atomics
 */
int A1D_Destroy_mutexes(A1_group_t* group);

/**
 * \brief Device level implementation of A1_Lock_mutex
 * 
 * Operation to lock a mutex. Blocks until lock has been acquired.
 *
 * \param[out]    rc            The error code.
 * \param[in]     group         A1 group over which the mutexes are shared.
 * \param[in]     mutex         A1 mutex.
 * \param[in]     proc          Process on which you want to lock mutex on
 *
 * \ingroup Atomics
 */
int A1D_Lock_mutex(A1_group_t* group, 
                   int mutex, 
                   int proc);

/**
 * \brief Device level implementation of A1_Trylock_mutex
 * 
 * Operation to trylock a mutex.
 *
 * \param[out]    rc            The error code.
 * \param[in]     group         A1 group over which the mutexes are shared.
 * \param[in]     mutex         A1 mutex.
 * \param[in]     proc          Process on which you want to lock mutex on
 * \param[out]    acquired      returns 1 if was acquired
 *
 * \ingroup Atomics
 */

int A1D_Trylock_mutex(A1_group_t* group, 
                      int mutex, 
                      int proc, 
                      A1_bool_t *acquired);

/**
 * \brief Device level implementation of A1_Unlock_mutex
 *  
 * Operation to unlock a mutex.  This call blocks until the mutex has been unlocked.
 *
 * \param[out]    rc            The error code.
 * \param[in]     group         A1 group over which the mutexes are shared.
 * \param[in]     mutex         A1 mutex.
 * \param[in]     proc          Process on which you want to lock mutex on
 *
 * \ingroup Atomics
 */

int A1D_Unlock_mutex(A1_group_t* group, 
                     int mutex, 
                     int proc);

/**
 * \brief Device level implementation of A1_Flush 
 * 
 *  On return, this call ensure that all blocking put or accumulate operations
 *  issued to a particular process are complete remotely.
 *
 * \param[in]  proc          Rank of the remote process.
 *
 * \ingroup COMPLETION
 */
int A1D_Flush(int proc);

/**
 * \brief Device level implementation of A1_Flush_group
 *
 *  On return, this call ensure that all blocking put or accumulate operations
 *  issued to the group of processes are complete remotely.
 *
 * \param[in]  group          Group of the remote processs.
 *
 * \ingroup COMPLETION
 */
int A1D_Flush_group(A1_group_t *group);

/**
 * \brief Reduce data from all processes and broadcast results to all processes.  
 *
 * \param[in] group          Group of processes.
 *
 * \see
 *
 * \ingroup MANYTOMANY
 */
int A1D_Allreduce_group(A1_group_t* group,
                       int count,
                       A1_reduce_op_t a1_op,
                       A1_datatype_t a1_type,
                       void* in,
                       void* out);

/**
 * \brief Reduce data from all processes and broadcast results to all processes.
 *
 * \param[in] group          Group of processes.
 *
 * \see
 *
 * \ingroup MANYTOMANY
 */
int A1D_NbAllreduce_group(A1_group_t* group,
                         int count,
                         A1_reduce_op_t a1_op,
                         A1_datatype_t a1_type,
                         void* in,
                         void* out,
                         A1_handle_t a1_handle);

/**
 * \brief
 *
 * \param[in] group          Group of processes.
 *
 * \see
 *
 * \ingroup MANYTOMANY
 */
int A1D_Bcast_group(A1_group_t* group,
                   int root,
                   int count,
                   void* buffer);

/**
 * \brief
 *
 * \param[in] group          Group of processes.
 *
 * \see
 *
 * \ingroup MANYTOMANY
 */

int A1D_NbBcast_group(A1_group_t* group,
                      int root,
                      int count,
                      void* buffer,
                      A1_handle_t a1_handle);

/**
 * \brief Device level implementation of A1_Process_id 
 *
 * Returns process rank relative to the group base specified.
 *
 * \param[out] rc          Process id in the process group.
 * \param[in]  group       Process group.
 *
 * \ingroup INFORMATION
 */
int A1D_Process_id(A1_group_t* group);

/**
 * \brief Device level implementation of A1_Process_total 
 *
 * Returns the total number of processes in  the group base specified.
 * 
 * \param[out] rc          Total number of processes in the process group.
 * \param[in]  group       Process group.
 *
 * \ingroup INFORMATION
 */
int A1D_Process_total(A1_group_t* group);

/**
 * \brief Device level implementation of A1_Node_id
 *
 * Returns node rank relative to the group base specified.
 *
 * \param[out] rc          Node id in the process group.
 * \param[in]  group       Process group.
 *
 * \ingroup INFORMATION
 */
int A1D_Node_id(A1_group_t* group);

/**
 * \brief Device level implementation of A1_Node_total
 *
 * Returns total number of nodes in the group base specified.
 * 
 * \param[out] rc          Total number of nodes in the process group.
 * \param[in]  group       Process group.
 *
 * \ingroup INFORMATION
 */
int A1D_Node_total(A1_group_t* group);

/**
 * \brief Device level implementation of A1_Time_seconds 
 * 
 * Timer in units of seconds.
 *
 * \param[out] rc          Number of secs from an arbitrary time in the past.
 *
 * \ingroup INFORMATION
 */
double A1D_Time_seconds(void);

/**
 * \brief Device level implementation of A1_Time_cycles 
 *
 * Timer in units of cycles.
 *
 * \param[out] rc          Number of cycles from an arbitrary time in the past.
 *
 * \ingroup INFORMATION
 */
unsigned long long A1D_Time_cycles(void);

void A1D_Global_lock_acquire(void);

void A1D_Global_lock_release(void);

#endif /* A1D_H_INCLUDED */
