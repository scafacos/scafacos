#ifndef _MP3_H_
#define _MP3_H_

#include <mpi.h>

#define MP_INIT(ARGC,ARGV)   MPI_Init(&(ARGC),&(ARGV))
#define MP_BARRIER()         MPI_Barrier(MPI_COMM_WORLD)
#define MP_FINALIZE()        MPI_Finalize()

#endif /* _MP3_H_ */
