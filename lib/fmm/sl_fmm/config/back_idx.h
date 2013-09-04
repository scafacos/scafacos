
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
#define SL_DATA0
#define sl_data0_type_c           long long
#define sl_data0_size_c           1
#define sl_data0_type_mpi         MPI_LONG_LONG
#define sl_data0_size_mpi         1
