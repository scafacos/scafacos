
#include "fortran2c_types.h"


/* enable runtime_informations */
/*#define SL_USE_RTI
#define SL_USE_RTI_TIM*/


/* standard (SL) integer data type */
#define sl_int_type_c          long long
#define sl_int_type_mpi        MPI_LONG_LONG
#define sl_int_size_mpi        1
#define sl_int_type_fmt        "lld"


/* index data type */
#define sl_index_type_c        FINT_TYPE_C
#define sl_index_type_mpi      FINT_TYPE_MPI
#define sl_index_size_mpi      1
#define sl_index_type_fmt      FINT_TYPE_FMT

/* use indices */
#define SL_INDEX


/* keys */
#define sl_key_type_c          FINT8_TYPE_C
#define sl_key_type_mpi        FINT8_TYPE_MPI
#define sl_key_size_mpi        1
#define sl_key_type_fmt        FINT8_TYPE_FMT
#define sl_key_integer

/* data0: work loads */
#define SL_DATA0
#define sl_data0_type_c        FREAL8_TYPE_C
#define sl_data0_size_c        1
#define sl_data0_type_mpi      FREAL8_TYPE_MPI
#define sl_data0_size_mpi      1

/* weighted elements */
#define sl_elem_weight(e, at)  ((e)->data0[at])

#define sl_data0_weight
