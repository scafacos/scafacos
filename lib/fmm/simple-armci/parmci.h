#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <mpi.h>

#if defined(DCMFD)
#  include "dcmfd/a1d_api.h"
#  define USE_32B_ATOMICS
#elif defined(PAMID)
#  include "pamid/a1d_api.h"
#  define USE_64B_ATOMICS
#elif defined(DMAPPD)
#  include "dmappd/a1d_api.h"
#  define USE_64B_ATOMICS
#else
#  error No device defined!
#endif

/* to make sure I haven't broken the API */
/* #include "parmci.h" */
/* #include "armci.h" */

typedef enum
{
    /* the following 32-but atomics are not supported on XE */
    ARMCI_FETCH,
    ARMCI_ADD,
    ARMCI_FETCH_AND_ADD,
    ARMCI_SWAP,
    /* the following 64-bit atomics are not supported on BGP */
    ARMCI_FETCH_LONG,
    ARMCI_ADD_LONG,
    ARMCI_FETCH_AND_ADD_LONG,
    ARMCI_SWAP_LONG
}
armci_rmw_t;

typedef enum
{
    ARMCI_ACC_DBL,
    ARMCI_ACC_FLT,
    ARMCI_ACC_CPL,
    ARMCI_ACC_DCP,
    ARMCI_ACC_INT,
    ARMCI_ACC_LNG
}
armci_acc_t;

/* for vector calls */
typedef struct
{
    void **src_ptr_array;
    void **dst_ptr_array;
    int  ptr_array_len;
    int bytes;
}
armci_giov_t;

typedef struct
{
    void * a1d_handle;
}
armci_hdl_t;

/* NOT USED
typedef struct armci_meminfo_ds
{
  void     *armci_addr;
  void     *addr;
  size_t    size;
  int       cpid;
  long      idlist[64];
}
armci_meminfo_t;
*/

/* initialization and termination */

int PARMCI_Init();
int PARMCI_Init_args(int *argc, char ***argv);

void PARMCI_Finalize();

/* memory management */

int PARMCI_Malloc(void * ptr_arr[], int bytes);
int PARMCIX_Malloc_comm(MPI_Comm comm, void * ptr_arr[], int bytes);

int PARMCI_Free(void * ptr);
int PARMCIX_Free_comm(MPI_Comm comm, void * ptr);

void * PARMCI_Malloc_local(int bytes);
int PARMCI_Free_local(void * ptr);

/* NOT USED
void *PARMCI_Memat(armci_meminfo_t * meminfo, int memflg);
void PARMCI_Memget(size_t bytes, armci_meminfo_t * meminfo, int memflg);
*/

/* synchronization */

void PARMCI_Barrier();
void PARMCIX_Barrier_comm(MPI_Comm comm);

void PARMCI_Fence(int proc);
void PARMCI_AllFence();
void PARMCIX_AllFence_comm(MPI_Comm comm);

int PARMCI_Test(armci_hdl_t * nb_handle);
int PARMCI_Wait(armci_hdl_t * nb_handle);
int PARMCI_WaitProc(int proc);
int PARMCI_WaitAll();

/* remote atomic update and mutexes */

#if defined(USE_32B_ATOMICS) && defined(USE_64B_ATOMICS)
#error The call syntax of (P)ARMCI_Rmw is stupid.  Use 32B or 64B atomics but not both.
#endif

#ifdef USE_32B_ATOMICS
int32_t PARMCI_Rmw(int optype, int32_t * local, int32_t * remote, int32_t incr, int proc);
#endif
#ifdef USE_64B_ATOMICS
int64_t PARMCI_Rmw(int optype, int64_t * local, int64_t * remote, int64_t incr, int proc);
#endif

int PARMCI_Create_mutexes(int num);
int PARMCI_Destroy_mutexes();
void PARMCI_Lock(int mutex, int proc);
void PARMCI_Unlock(int mutex, int proc);

/* blocking one-sided */

int PARMCI_Acc(armci_acc_t optype, void *scale, void *src, void* dst, int bytes, int proc);
int PARMCI_AccS(armci_acc_t optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc);
int PARMCI_AccV(armci_acc_t optype, void *scale, armci_giov_t * darr, int len, int proc);
int PARMCI_Get(void *src, void *dst, int bytes, int proc);
int PARMCI_GetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc);
int PARMCI_GetV(armci_giov_t * darr, int len, int proc);
int PARMCI_Put(void *src, void *dst, int bytes, int proc);
int PARMCI_PutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc);
int PARMCI_PutV(armci_giov_t * darr, int len, int proc);

/* non-blocking one-sided */

int PARMCI_NbAccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t * nb_handle);
int PARMCI_NbAccV(int op, void *scale, armci_giov_t * darr, int len, int proc, armci_hdl_t * nb_handle);
int PARMCI_NbPut(void *src, void *dst, int bytes, int proc, armci_hdl_t * nb_handle);
int PARMCI_NbPutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t * nb_handle);
int PARMCI_NbPutV(armci_giov_t * darr, int len, int proc, armci_hdl_t * nb_handle);
int PARMCI_NbGet(void *src, void *dst, int bytes, int proc, armci_hdl_t * nb_handle);
int PARMCI_NbGetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t * nb_handle);
int PARMCI_NbGetV(armci_giov_t * darr, int len, int proc, armci_hdl_t * nb_handle);

int PARMCI_Put_flag(void *src, void *dst, int bytes, int *f, int v, int proc);
int PARMCI_PutS_flag(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc);
int PARMCI_PutS_flag_dir(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc);

/* CAF extensions */

int PARMCI_PutValueInt(int src, void *dst, int proc);
int PARMCI_PutValueLong(long src, void *dst, int proc);
int PARMCI_PutValueFloat(float src, void *dst, int proc);
int PARMCI_PutValueDouble(double src, void *dst, int proc);

int PARMCI_GetValueInt(void *src, int proc);
long PARMCI_GetValueLong(void *src, int proc);
float PARMCI_GetValueFloat(void *src, int proc);
double PARMCI_GetValueDouble(void *src, int proc);

/**********************************************/
