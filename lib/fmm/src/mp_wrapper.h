! get the defines computed by configure
#include <fconfig.h>

!define message passing libraries variables for comparison
#define FMM_MP_ARMCI 1
#define FMM_MP_A1 2
#define FMM_MP_MPI 3
#define FMM_MP_SIMPLE_ARMCI 4

! enable warnings if a sendbuffer is larger than 2GB
#undef FMM_CHECKSIZE

! enable 16Byte mpi_walltime
#undef FMM_TIMER16
