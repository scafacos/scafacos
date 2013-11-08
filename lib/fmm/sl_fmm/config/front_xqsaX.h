
#include "config_fmm_sort.h"


/* standard (SL) integer data type */
#define sl_int_type_c             SL_INTEGER_C
#define sl_int_type_mpi           SL_INTEGER_MPI
#define sl_int_size_mpi           1
#define sl_int_type_fmt           SL_INTEGER_FMT


/* key section */
#define sl_key_type_c             INTEGER_C
#define sl_key_type_mpi           INTEGER_MPI
#define sl_key_size_mpi           1

#define sl_key_integer
#define sl_key_type_fmt           INTEGER_FMT

#define sl_key_purify(_k_)        (_k_ & SL_FMM_CONFIG_VAR(fmm_front_key_mask))


/* data section */
#define SL_DATA0                  /* xyz */
#define sl_data0_type_c           REAL_C
#define sl_data0_size_c           3
#define sl_data0_type_mpi         REAL_MPI
#define sl_data0_size_mpi         3

#define SL_DATA1                  /* q */
#define sl_data1_type_c           REAL_C
#define sl_data1_size_c           1
#define sl_data1_type_mpi         REAL_MPI
#define sl_data1_size_mpi         1

#define SL_DATA2                  /* scr */
#define sl_data2_type_c           INTEGER_C
#define sl_data2_size_c           1
#define sl_data2_type_mpi         INTEGER_MPI
#define sl_data2_size_mpi         1

#define SL_DATA3                  /* addr */
#define sl_data3_type_c           char
#define sl_data3_size_c           SL_FMM_CONFIG_VAR(fmm_front_aX)
#define sl_data3_type_mpi         MPI_CHAR
#define sl_data3_size_mpi         SL_FMM_CONFIG_VAR(fmm_front_aX)

#define sl_data3_flex
#define sl_data3_memcpy
