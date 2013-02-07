/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

/** @file dcmfdimpl.h */

/*! \addtogroup a1 A1D dcmfd device interface
 * @{
 */

#include "a1.h"
#include "a1u.h"
#include "a1d.h"
#include <dcmf.h>
#include <dcmf_globalcollectives.h>
#include <dcmf_collectives.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <bpcore/bgp_atomic_ops.h>
#include <spi/bgp_SPI.h>

/*************************************************
 *                 Constants                     *
 ************************************************/

#define A1C_ALIGNMENT 16

#define A1C_ENABLE_CHT 1
#define A1C_ENABLE_INTERRUPTS 0
#define A1C_MPI_ACTIVE 1
#define A1C_CHT_PAUSE_CYCLES 200

#define A1C_PUT_PACKING_LIMIT 512
#define A1C_GET_PACKING_LIMIT 512
#define A1C_PUTACC_PACKING_LIMIT 2048

#define A1C_PUT_PACKETSIZE 2048 
#define A1C_GET_PACKETSIZE 2048 
#define A1C_PUTACC_PACKETSIZE 8192 

#define A1C_FLUSHALL_PENDING_LIMIT 512 

#define A1C_REQUEST_POOL_SIZE 100 

#define A1C_HANDLE_POOL_SIZE 100

#define A1C_MAX_STRIDED_DIM 8

#define A1C_USE_HANDOFF 1 

/* Currently we have two sizes of buffers to provide for Put/Get and
 * Acc packets. */
#define A1C_BUFFER_SIZES 3

#define A1C_PUT_BUFFERPOOL_SIZE 200;
#define A1C_GET_BUFFERPOOL_SIZE 200;
#define A1C_PUTACC_BUFFERPOOL_SIZE 200;

/*************************************************
*                  BGP Atomics                   *
*************************************************/

extern _BGP_Atomic global_atomic;

#define A1DI_GLOBAL_ATOMIC_ACQUIRE()                 \
 {                                                   \
   volatile int done=0;                              \
   do {                                              \
     while(global_atomic.atom);                      \
     done = _bgp_test_and_set(&global_atomic, 1);    \
   } while(!done);                                   \
 }                                                   \

#define A1DI_GLOBAL_ATOMIC_RELEASE() do{ global_atomic.atom = 0; _bgp_mbar(); }while(0)

/*************************************************
*                  Lockbox                       *
*************************************************/

extern LockBox_Mutex_t global_lbmutex;

/* Different cores which want to use independent lockbox mutexes should 
 * use different counters. So we try to find a free counter in a non-overlapping 
 * range of 200 counters. Counters range from 0-1023 */
#define A1DI_GLOBAL_LBMUTEX_INITIALIZE()            	     	     \
 do {                                                  	     	     \
   int idx, coreid;                                   	     	     \
   coreid = Kernel_PhysicalProcessorID();                            \
   for(idx=200*coreid; idx<200*(coreid+1); idx++)          	     \
   {                                                         	     \
     if(!LockBox_AllocateMutex(idx, &global_lbmutex, coreid, 1, 0))  \
          break;					     	     \
   }       						             \
   A1U_ERR_POP(idx == 200*(coreid+1),			 	     \
         "LockBox_AllocateMutex did not find a free index \n");      \
 } while(0)	                                                     \
       
#define A1DI_GLOBAL_LBMUTEX_ACQUIRE() LockBox_MutexLock(global_lbmutex);

#define A1DI_GLOBAL_LBMUTEX_RELEASE() LockBox_MutexUnlock(global_lbmutex);

/*************************************************
*           Lock Type Selection                  *
*************************************************/

#define A1DI_GLOBAL_LOCK_ACQUIRE()     \
 do {                                  \
    if(!a1d_settings.mpi_active)       \
    {                                  \
         A1DI_GLOBAL_ATOMIC_ACQUIRE(); \
    }                                  \
    else                               \
    {                                  \
        DCMF_CriticalSection_enter(0); \
    }                                  \
 } while(0);                           \

#define A1DI_GLOBAL_LOCK_RELEASE()     \
 do {                                  \
    if(!a1d_settings.mpi_active)       \
    {                                  \
        A1DI_GLOBAL_ATOMIC_RELEASE();  \
    }                                  \
    else                               \
    {                                  \
        DCMF_CriticalSection_exit(0);  \
    }                                  \
 } while(0);                           \

/*************************************************
*           Likely and Unlikely Ifs              *
*************************************************/

#define likely_if(x) if(__builtin_expect(x,1))
#define unlikely_if(x) if(__builtin_expect(x,0))

/*************************************************
 *          Generic  Macros                      *
 *************************************************/

#define A1DI_Wait_cycles(cycles)                     \
   do {                                              \
      unsigned long long start = DCMF_Timebase();    \
      while((DCMF_Timebase() - start) < cycles);     \
   } while(0)                                        \

#define A1DI_Wait_seconds(seconds)               \
   do {                                           \
      double start = DCMF_Timer();               \
      while((DCMF_Timer() - start) < seconds);   \
   } while(0)                                    \

#define A1DI_Set_handle(request, handle)  \
do {                                      \
    request->handle_ptr = handle;         \
   } while(0)                             \

/*************************************************
 *          Memory Allocation Macros             *
 *************************************************/

#define A1DI_Malloc(ptr, num) posix_memalign(ptr, a1d_settings.alignment, num)
/*
 * I don't know why one would want unaligned memory, but here it is for posterity
 * #define A1DI_Malloc(ptr, num)  ((ptr = malloc(num)) == NULL)
 */

#define A1DI_Free(ptr) free(ptr)

#define A1DI_Memset(ptr, val, num)  memset(ptr, val, num)

#define A1DI_Memcpy(trg, src, num)  memcpy(trg, src, num)

/*************************************************
 *          Critical Section Macros              *
 *************************************************/

/**
 * \brief Handles non-contiguous puts which have been handed-off to the CHT.
 *
 * \see A1D_NbPutS, A1DI_CHT_advance_lock
 *
 * \ingroup CHT
 */
void A1DI_Handoff_progress();

#define A1DI_CRITICAL_ENTER()                                     \
    do {                                                          \
      if(a1d_settings.enable_cht || a1d_settings.mpi_active)      \
      {                                                           \
        A1DI_GLOBAL_LOCK_ACQUIRE();                               \
      }     							  \
    } while (0)                                                   \

#define A1DI_CRITICAL_EXIT()                                      \
    do {                                                          \
      if(a1d_settings.enable_cht || a1d_settings.mpi_active)      \
      {                                                           \
        A1DI_GLOBAL_LOCK_RELEASE();                               \
      }                                                           \
    } while (0)                                                   \

#define A1DI_Advance()                                             \
 do {                                                              \
         DCMF_Messager_advance(0);                                 \
         if (a1d_settings.use_handoff && (A1D_Inside_handoff==0))  \
         {        						   \
             A1DI_Handoff_progress();                              \
         }                                                         \
    } while(0)                                                     \

#define A1DI_Conditional_advance(boolean)                           \
    while(boolean) {                                                \
          DCMF_Messager_advance(0);                                 \
          if (a1d_settings.use_handoff && (A1D_Inside_handoff==0))  \
          {							    \
                A1DI_Handoff_progress();                            \
          }							    \
    }                                                               \

/*************************************************
 *          Computation macros                   *
 *************************************************/

#define A1DI_ACC(datatype, source, target, scaling, count)                  \
   do {                                                                     \
     int i;                                                                 \
     datatype *s = (datatype *) source;                                     \
     datatype *t = (datatype *) target;                                     \
     datatype c = (datatype) scaling;                                       \
     for(i=0; i<count; i++)                                                 \
          t[i] += s[i]*c;                                                   \
   } while(0)                                                               \

#define A1DI_MOD_BXOR(datatype, source, target, count)                      \
   do {                                                                     \
     int i;                                                                 \
     datatype *s = (datatype *) source;                                     \
     datatype *t = (datatype *) target;                                     \
     for(i=0; i<count; i++)                                                 \
          t[i] ^= s[i];                                                     \
   } while(0)                                                               \

 /* TODO probably need to optimize these functions for double-hummer */
#define A1DI_ACC_DOUBLE(source, target, scaling, count)                  \
   do {                                                                     \
     int i;                                                                 \
     double *s = (double *) source;                                     \
     double *t = (double *) target;                                     \
     double c = (double) scaling;                                       \
     for(i=0; i<count; i++)                                                 \
          t[i] += s[i]*c;                                                   \
   } while(0)                                                               \

#define A1DI_ABS(datatype, source, target, count)                           \
   do {                                                                     \
     int i;                                                                 \
     datatype *s = (datatype *) source;                                     \
     datatype *t = (datatype *) target;                                     \
     for(i=0; i<count; i++) t[i] = ( s[i] > 0 ? s[i] : -s[i]);              \
   } while(0)                                                               \

/* NOTE: fabs will compile to the best assembly. */
 #define A1DI_ABS_DOUBLE(source, target, count)                           \
    do {                                                                  \
      int i;                                                              \
      double *s = (double *) source;                                     \
      double *t = (double *) target;                                     \
      for(i=0; i<count; i++) t[i] = fabs(s[i]);                            \
    } while(0)                                                             \

#define A1DI_FETCHANDADD_EXECUTE(datatype, source, target, original, count) \
   do {                                                                     \
     int i;                                                                 \
     datatype *s = (datatype *) source;                                     \
     datatype *t = (datatype *) target;                                     \
     datatype *o = (datatype *) original;                                   \
     for(i=0; i<count; i++)                                                 \
     {                                                                      \
          o[i] = t[i];                                                      \
          t[i] += s[i];                                                     \
     }                  						                            \
   } while(0)                                                              \

/*************************************************
 *             Data Structures                   *
 *************************************************/
typedef enum
{
  A1D_MUTEX_LOCK = 0,
  A1D_MUTEX_TRYLOCK,
  A1D_MUTEX_UNLOCK
} A1D_Mutex_op_type;

typedef struct
{
    volatile uint32_t enable_cht;
    volatile uint32_t mpi_active;
    volatile uint32_t cht_pause_cycles;
    volatile uint32_t enable_interrupts;
    volatile uint32_t put_packing_limit;
    volatile uint32_t get_packing_limit;
    volatile uint32_t putacc_packing_limit;
    volatile uint32_t put_packetsize;
    volatile uint32_t get_packetsize;
    volatile uint32_t putacc_packetsize;
    volatile uint32_t flushall_pending_limit;
    volatile uint32_t alignment;
    volatile uint32_t handlepool_size;
    volatile uint32_t requestpool_size;
    volatile uint32_t put_bufferpool_size;
    volatile uint32_t get_bufferpool_size;
    volatile uint32_t putacc_bufferpool_size;
    volatile uint32_t use_handoff;
} A1D_Settings_t;

typedef struct
{
    size_t my_rank;
    size_t my_node;
    size_t num_ranks;
    size_t num_nodes;
    DCMF_Hardware_t hw;
} A1D_Process_info_t;

typedef struct
{
    int   rank;
    long  *value_ptr; /*This will hold ptr it value if counter
                        is located on a remote process*/
    long  value;      /*This will contain the value if counter
                        is located locally*/
} A1D_Counter_t;

typedef struct A1D_Mutex_request_t
{
  int rank;
  struct A1D_Mutex_request_t *next;
} A1D_Mutex_request_t;

typedef struct
{
  int   mutex;
  A1D_Mutex_request_t *head;
  A1D_Mutex_request_t *tail;
} A1D_Mutex_t;

typedef struct A1D_Buffer_t
{
  void *buffer_ptr;
  int pool_index;
  struct A1D_Buffer_t *next;
} A1D_Buffer_t;

typedef struct
{
  A1D_Buffer_t* pool_heads[A1C_BUFFER_SIZES];
  int limits[A1C_BUFFER_SIZES];
  int sizes[A1C_BUFFER_SIZES];
  void* pool_region_ptrs[A1C_BUFFER_SIZES];
  void* mem_region_ptrs[A1C_BUFFER_SIZES];
} A1D_Buffer_pool_t;

typedef struct A1D_Handle_t
{
    volatile int active;
    volatile int active_list_index;
    struct A1D_Handle_t *next;
} A1D_Handle_t;

typedef struct A1D_Handle_pool_t
{
    A1D_Handle_t *head;
    void *region_ptr;
} A1D_Handle_pool_t;

typedef struct A1D_Request_t
{
    DCMF_Request_t request;
    int in_pool;
    void* buffer_ptr;
    A1D_Buffer_t *a1d_buffer_ptr;
    uint32_t buffer_size;
    A1D_Handle_t *handle_ptr;
    struct A1D_Request_t *next;
} A1D_Request_t;

typedef struct
{
    A1D_Request_t *head;
    A1D_Request_t *region_ptr;
} A1D_Request_pool_t;

typedef struct
{
    DCMF_Protocol_t protocol;
    volatile int rcv_active;
    void **xchange_ptr;
    size_t xchange_size;
} A1D_Control_xchange_info_t;

typedef union
{
    DCQuad info[2];
    struct
    {
        void* target_ptr;
        A1_datatype_t datatype;
        union
        {
            int32_t int32_value;
            int64_t int64_value;
            uint32_t uint32_value;
            uint64_t uint64_value;
            float float_value;
            double double_value;
        } scaling;
    };
} A1D_Putacc_header_t;

typedef union
{
    DCQuad info[2];
    struct
    {
        void* target_ptr;
        A1_reduce_op_t op;
        A1_datatype_t datatype;
    };
} A1D_Putmod_header_t;

typedef struct
{
    int stride_level;
    int block_sizes[A1C_MAX_STRIDED_DIM];
    void *target_ptr;
    int trg_stride_ar[A1C_MAX_STRIDED_DIM-1];
    int block_idx[A1C_MAX_STRIDED_DIM];
    int data_size; 
} A1D_Packed_puts_header_t;

typedef struct
{
    int stride_level;
    int block_sizes[A1C_MAX_STRIDED_DIM];
    void *target_ptr;
    int trg_stride_ar[A1C_MAX_STRIDED_DIM-1];
    int block_idx[A1C_MAX_STRIDED_DIM];
    int data_size;
    A1D_Handle_t *handle_ptr;
} A1D_Packed_gets_response_header_t;

typedef struct
{
    uint32_t target;
    int stride_level;
    int block_sizes[A1C_MAX_STRIDED_DIM];
    void* source_ptr;
    int src_stride_ar[A1C_MAX_STRIDED_DIM-1];
    void* target_ptr;
    int trg_stride_ar[A1C_MAX_STRIDED_DIM-1];
    A1D_Handle_t *handle_ptr;
} A1D_Packed_gets_header_t;

typedef union
{
   DCQuad info[2];
   struct
   {
     int bytes;
     int source;
     void* source_ptr_out; 
     void* target_ptr;
     A1_atomic_op_t op;
     A1_datatype_t datatype;
     A1D_Handle_t* handle_ptr;
   };
} A1D_Rmw_header_t;

typedef struct
{
   int bytes;
   void* source_ptr_out; 
   A1D_Handle_t* handle_ptr;
} A1D_Rmw_response_header_t;

typedef struct
{
    long *value_ptr;
    long value;
} A1D_Counter_pkt_t;

typedef struct
{
    int mutex_idx;
    A1D_Mutex_op_type mutex_op;
    int response;
} A1D_Mutex_pkt_t;

typedef struct
{
    int stride_level;
    int block_sizes[A1C_MAX_STRIDED_DIM];
    void *target_ptr;
    int trg_stride_ar[A1C_MAX_STRIDED_DIM-1];
    int block_idx[A1C_MAX_STRIDED_DIM];
    int data_size;
    A1_datatype_t datatype;
    union
    {
        int32_t int32_value;
        int64_t int64_value;
        uint32_t uint32_value;
        uint64_t uint64_value;
        float float_value;
        double double_value;
    } scaling;
} A1D_Packed_putaccs_header_t;


/*************************************************
 *             Global variables                  *
 ************************************************/

extern pthread_t A1DI_CHT_pthread;

extern A1D_Process_info_t A1D_Process_info;
extern A1D_Control_xchange_info_t A1D_Control_xchange_info;
extern A1D_Request_pool_t A1D_Request_pool;
extern A1D_Handle_pool_t A1D_Handle_pool;
extern A1D_Buffer_pool_t A1D_Buffer_pool;

extern int *A1D_Mutexes_count;
extern A1D_Mutex_t *A1D_Mutexes;

extern DCMF_Configure_t A1D_Messager_info;
extern DCMF_Protocol_t A1D_Control_flushack_protocol;
extern DCMF_Protocol_t A1D_Send_flush_protocol;
extern DCMF_Protocol_t A1D_GlobalBarrier_protocol;
extern DCMF_CollectiveProtocol_t A1D_GlobalAllreduce_protocol;
extern DCMF_Protocol_t A1D_GlobalBcast_protocol;
extern DCMF_Protocol_t A1D_Generic_put_protocol;
extern DCMF_Protocol_t A1D_Generic_get_protocol;
extern DCMF_Protocol_t A1D_Generic_putacc_protocol;
extern DCMF_Protocol_t A1D_Generic_putmod_protocol;
extern DCMF_Protocol_t A1D_Rmw_protocol;
extern DCMF_Protocol_t A1D_Rmw_response_protocol;
extern DCMF_Protocol_t A1D_Packed_puts_protocol;
extern DCMF_Protocol_t A1D_Packed_gets_protocol;
extern DCMF_Protocol_t A1D_Packed_gets_response_protocol;
extern DCMF_Protocol_t A1D_Packed_putaccs_protocol;
extern DCMF_Protocol_t A1D_Counter_create_protocol;
extern DCMF_Protocol_t A1D_Counter_protocol;
extern DCMF_Protocol_t A1D_Mutex_protocol;
extern DCMF_Protocol_t A1D_Control_protocol;
extern DCMF_Callback_t A1D_Nocallback;
extern DCMF_Memregion_t *A1D_Memregion_global;

extern void **A1D_Membase_global;
extern void **A1D_Put_Flushcounter_ptr;
extern A1D_Handle_t **A1D_Active_handle_list;
extern volatile int *A1D_Connection_send_active;
extern volatile int *A1D_Connection_put_active;
extern volatile int A1D_Control_flushack_active;
extern volatile int A1D_Put_flushack_active;
extern volatile int A1D_Inside_handoff;

extern A1D_Settings_t a1d_settings;

/************************************************* 
 *             Function Prototypes               *
 ************************************************/

void *A1DI_CHT_advance_lock(void *);

void A1DI_Global_lock_acquire();

void A1DI_Global_lock_release();

void A1DI_Generic_done(void *, DCMF_Error_t *);

void A1DI_Request_done(void *, DCMF_Error_t *);

int A1DI_Memregion_Global_initialize();

int A1DI_Put_initialize();

int A1DI_Packed_puts_initialize();

int A1DI_Get_initialize();

int A1DI_Rmw_initialize();

int A1DI_Counter_initialize();

int A1DI_Request_pool_initialize();

void A1DI_Request_pool_finalize();

A1D_Request_t* A1DI_Get_request(int);

void A1DI_Release_request(A1D_Request_t *);

int A1DI_Buffer_pool_initialize();

void A1DI_Buffer_pool_finalize();

A1D_Buffer_t* A1DI_Get_buffer(int, int);

void A1DI_Release_buffer(A1D_Buffer_t *);

int A1DI_Handle_pool_initialize();

void A1DI_Handle_pool_finalize();

A1D_Handle_t* A1DI_Get_handle();

void A1DI_Release_handle(A1D_Handle_t *);

int A1DI_Packed_gets_initialize();

int A1DI_Putacc_initialize();

int A1DI_Putmod_initialize();

int A1DI_Packed_putaccs_initialize();

int A1DI_GlobalBarrier_initialize();

int A1DI_GlobalAllreduce_initialize();

int A1DI_GlobalBcast_initialize();

int A1DI_Send_flush_initialize();

int A1DI_Put_flush_initialize();

int A1DI_Control_flushack_initialize();

int A1DI_GlobalBarrier();

int A1DI_Read_parameters();

int A1DI_Send_flush(int proc);

int A1DI_Pack_strided(void *packet_ptr,
                      int packet_limit,
                      int stride_level,
                      int *block_sizes,
                      void **source_ptr,
                      int *src_stride_ar,
                      void **target_ptr,
                      int *trg_stride_ar,
                      int *block_idx,
                      int *data_size,
                      int *complete);

int A1DI_Unpack_strided(void *packet_ptr,
                        int data_size,
                        int stride_level,
                        int *block_sizes,
                        void *target_ptr,
                        int *trg_stride_ar,
                        int *block_idx,
                        int *complete);

int A1DI_Unpack_strided_acc(void *data_ptr,
                            int data_size,
                            int stride_level,
                            int *block_sizes,
                            void *target_ptr,
                            int *trg_stride_ar,
                            int *block_idx,
                            A1_datatype_t a1_type,
                            void *scaling,
                            int *complete);

/*****************************************************
                 Packing Handoff
*****************************************************/

typedef enum
{
  A1D_PACKED_PUTS = 0,
  A1D_PACKED_PUTACCS
} A1D_Op_type;

typedef struct
{
   int target;
   int stride_level;
   int *block_sizes;
   void *source_ptr;
   int *src_stride_ar;
   void *target_ptr;
   int *trg_stride_ar;
   A1D_Handle_t *a1d_handle;
} A1D_Puts_op;

typedef struct
{
   int target;
   int stride_level;
   int *block_sizes;
   void *source_ptr;
   int *src_stride_ar;
   void *target_ptr;
   int *trg_stride_ar;
   A1_datatype_t datatype;
   void *scaling;
   A1D_Handle_t *a1d_handle;
} A1D_Putaccs_op;

typedef struct A1D_Op_handoff
{
   A1D_Op_type op_type;
   union
   { 
     A1D_Puts_op puts_op;
     A1D_Putaccs_op putaccs_op;
   } op;  
   void *op_ptr;
   struct A1D_Op_handoff *next;
} A1D_Op_handoff;

extern A1D_Op_handoff *A1D_Op_handoff_queuehead;
extern A1D_Op_handoff *A1D_Op_handoff_queuetail;

/*! @} */
