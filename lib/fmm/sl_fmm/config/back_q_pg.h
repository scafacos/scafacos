
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


/* data section */
#define SL_DATA0                  /* q */
#define sl_data0_type_c           REAL_C
#define sl_data0_size_c           1
#define sl_data0_type_mpi         REAL_MPI
#define sl_data0_size_mpi         1

#undef SL_DATA1                   /* xyz */
#define sl_data1_type_c           REAL_C
#define sl_data1_size_c           3
#define sl_data1_type_mpi         REAL_MPI
#define sl_data1_size_mpi         3

#define SL_DATA2                  /* pot */
#define sl_data2_type_c           REAL_C
#define sl_data2_size_c           1
#define sl_data2_type_mpi         REAL_MPI
#define sl_data2_size_mpi         1

#define SL_DATA3                  /* grad */
#define sl_data3_type_c           REAL_C
#define sl_data3_size_c           3
#define sl_data3_type_mpi         REAL_MPI
#define sl_data3_size_mpi         3

#undef SL_DATA4                   /* load */
#define sl_data4_type_c           REAL_C
#define sl_data4_size_c           1
#define sl_data4_type_mpi         REAL_MPI
#define sl_data4_size_mpi         1
