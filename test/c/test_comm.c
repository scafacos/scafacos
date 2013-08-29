
#include <stdio.h>
#include <stdlib.h>

#include "fcs.h"


int size, rank;

int test_solver(const char *method, MPI_Comm sub_comm)
{
  int sub_size, sub_rank;

  fcs_float offset[] = { 0.0, 0.0, 0.0 };
  fcs_float box_a[] = { 1.0, 0.0, 0.0 };
  fcs_float box_b[] = { 0.0, 1.0, 0.0 };
  fcs_float box_c[] = { 0.0, 0.0, 1.0 };
  fcs_int periodicity[] = { 1, 1, 1 };

  fcs_int local_nparticles = 4;
  fcs_float positions[local_nparticles * 3], charges[local_nparticles];
  fcs_float field[local_nparticles * 3], potentials[local_nparticles];

  FCS fcs;
  FCSResult result;

  MPI_Comm_size(sub_comm, &sub_size);
  MPI_Comm_rank(sub_comm, &sub_rank);

  srandom(sub_rank * 2501);

  for (fcs_int i = 0; i < local_nparticles; ++i)
  {
    positions[3 * i + 0] = (fcs_float) random() / (fcs_float) RAND_MAX;
    positions[3 * i + 1] = (fcs_float) random() / (fcs_float) RAND_MAX;
    positions[3 * i + 2] = (fcs_float) random() / (fcs_float) RAND_MAX;
    charges[i] = 1.0;
  }

  result = fcs_init(&fcs, method, sub_comm);
  if (result != FCS_RESULT_SUCCESS)
  {
    if (fcsResult_getReturnCode(result) == FCS_WRONG_ARGUMENT)
    {
      printf("%d: solver method '%s' not available!\n", rank, method);
      return 10;
    }
    printf("%d: fcs_init failed!\n", rank);
    return 11;
  }

  result = fcs_set_common(fcs, 1, box_a, box_b, box_c, offset, periodicity, local_nparticles * sub_size);
  if (result != FCS_RESULT_SUCCESS) { printf("%d: fcs_set_common failed!\n", rank); return 12; }

  result = fcs_run(fcs, local_nparticles, local_nparticles, positions, charges, field, potentials);
  if (result != FCS_RESULT_SUCCESS) { printf("%d: fcs_run failed!\n", rank); return 13; }

  result = fcs_destroy(fcs);
  if (result != FCS_RESULT_SUCCESS) { printf("%d: fcs_destroy failed!\n", rank); return 14; }

  return 0;
}


int main(int argc, char *argv[])
{
  MPI_Comm comm, dup_comm, sub_comm;
  int err = 0;
  double timeout = 2;
#if FCS_ENABLE_PEPC
  int mpi_thread_requested = MPI_THREAD_MULTIPLE;
  int mpi_thread_provided;

  MPI_Init_thread(&argc, &argv, mpi_thread_requested, &mpi_thread_provided);
#else
  MPI_Init(&argc, &argv);
#endif

  comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (argc < 2) goto exit;

  if (argc > 2) timeout = atof(argv[2]);

  if (rank == 0)
  {
    printf("Communicator test with method: '%s', timeout: %f second(s)\n", argv[1], timeout);
  }

  MPI_Comm_dup(comm, &dup_comm);
  MPI_Comm_split(comm, rank % 2, -rank, &sub_comm);

  if (rank % 2) err = test_solver(argv[1], sub_comm);
  else MPI_Barrier(sub_comm);

/*  if (err) goto exit;*/
  
  int b = 0xdeadbeef;

  if (rank > 0) MPI_Send(&b, 1, MPI_INT, 0, 0, dup_comm);
  else
  {
    MPI_Request req = MPI_REQUEST_NULL;
    MPI_Status status;
    int flag;
    int r = size - 1;
    double t = MPI_Wtime() + timeout;

    while (r > 0)
    {
      if (req == MPI_REQUEST_NULL) MPI_Irecv(&b, 1, MPI_INT, MPI_ANY_SOURCE, 0, dup_comm, &req);

      MPI_Test(&req, &flag, &status);

      if (flag) --r;

      if (MPI_Wtime() > t) break;
    }
    
    if (r > 0) MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Bcast(&err, 1, MPI_INT, 1, comm);

exit:
  MPI_Finalize();

  return err;
}
