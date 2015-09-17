
#define sl_key_type_c        long long
#define sl_key_type_mpi      MPI_LONG_LONG
#define sl_key_size_mpi      1
#define sl_key_type_fmt      "lld"

#define sl_key_integer

#define SL_DATA0             /* positions */
#define sl_data0_type_c      fcs_float
#define sl_data0_size_c      3
#define sl_data0_type_mpi    FCS_MPI_FLOAT
#define sl_data0_size_mpi    3

#define SL_DATA1             /* moments */
#define sl_data1_type_c      fcs_float
#define sl_data1_size_c      3
#define sl_data1_type_mpi    FCS_MPI_FLOAT
#define sl_data1_size_mpi    3

#define SL_DATA2             /* indices */
#define sl_data2_type_c      long long
#define sl_data2_size_c      1
#define sl_data2_type_mpi    MPI_LONG_LONG
#define sl_data2_size_mpi    1

#undef SL_DATA3              /* field */
#define sl_data3_type_c      fcs_float
#define sl_data3_size_c      6
#define sl_data3_type_mpi    FCS_MPI_FLOAT
#define sl_data3_size_mpi    6

#define SL_DATA4             /* potentials */
#define sl_data4_type_c      fcs_float
#define sl_data4_size_c      3
#define sl_data4_type_mpi    FCS_MPI_FLOAT
#define sl_data4_size_mpi    3
