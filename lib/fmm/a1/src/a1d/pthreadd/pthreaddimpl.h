/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  Copyright (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <pthread.h>
#include <assert.h>
#include "a1.h"
#include "a1u.h"
#include "a1d.h"

/*************************************************
 *                 Constants                     *
 ************************************************/

#define A1C_ALIGNMENT 16
#define A1C_MAX_STRIDED_DIM 4

/*************************************************
*                  Global Lock                   *
*************************************************/

extern pthread_mutex_t global_mutex;

#define A1DI_GLOBAL_MUTEX_ACQUIRE()                 \
 {                                                   \
    int rc = 0                                      \
    rc = pthread_mutex_lock(&global_mutex);          \
    A1U_ERR_POP(rc != 0,                             \
                "pthread_mutex_lock returned %d\n",  \
                rc);                                 \
 }                                                   \

#define A1DI_GLOBAL_MUTEX_RELEASE()                     \
 {                                                      \
     int rc = 0                                         \
     rc = pthread_mutex_unlock(&global_mutex);          \
     A1U_ERR_POP(rc != 0,                               \
                 "pthread_mutex_unlock returned %d\n",  \
                 rc);                                   \
 }

/*************************************************
 *          Memory Allocation Macros             *
 *************************************************/

//#define A1DI_Malloc_aligned(ptr, num) posix_memalign(ptr, a1d_settings.alignment, num)
#define A1DI_Malloc_aligned(ptr, num) ((*ptr = malloc(num)) == NULL)

#define A1DI_Malloc(ptr, num)  ((ptr = malloc(num)) == NULL)

#define A1DI_Free(ptr) free(ptr)

#define A1DI_Memset(ptr, val, num)  memset(ptr, val, num)

#define A1DI_Memcpy(trg, src, num)  memcpy(trg, src, num)

/*************************************************
 *          Critical Section Macros              *
 *************************************************/

#define A1DI_CRITICAL_ENTER()                                    \
    do {                                                          \
      if(a1d_settings.thread_safe)                                 \
      {                                                           \
          A1DI_GLOBAL_MUTEX_ACQUIRE();                            \
      }                                                           \
    } while (0)                                                  \

#define A1DI_CRITICAL_EXIT()                                     \
    do {                                                          \
      if(a1d_settings.thread_safe)                                 \
      {                                                           \
          A1DI_GLOBAL_MUTEX_RELEASE();                            \
      }                                                           \
    } while (0)                                                  \

/*************************************************
 *          Non-contiguous macros                *
 *************************************************/

#define A1DI_ACC_EXECUTE(datatype, source, target, scaling, count)          \
   do {                                                                     \
     int i;                                                                 \
     datatype *a = (datatype *) source;                                     \
     datatype *b = (datatype *) target;                                     \
     datatype c = (datatype) scaling;                                       \
     for(i=0; i<count; i++)                                                 \
          b[i] = b[i] + a[i]*c;                                             \
   } while(0)                                                               \

/*************************************************
 *             Data Structures                   *
 *************************************************/
typedef struct
{
    volatile uint32_t thread_safe;
    volatile uint32_t alignment;
} A1D_Settings_t;

typedef struct
{
    size_t my_rank;
    size_t my_node;
    size_t num_ranks;
    size_t num_nodes;
} A1D_Process_info_t;

/*************************************************
 *             Global variables                  *
 ************************************************/

extern A1D_Process_info_t A1D_Process_info;
extern A1D_Settings_t a1d_settings;

/*************************************************
 *             Function Prototypes               *
 ************************************************/

