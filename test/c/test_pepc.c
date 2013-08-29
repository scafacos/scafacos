/**
 * @file   test_pepc.c
 * @author Lukas Arnold <l.arnold@fz-juelich.de>
 * @date   today
 *
 * @brief  Brief test of the pepc lib and interface.
 *         Most of the code is taken from test_direct_new.c.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _USE_MATH_CONSTANTS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "fcs.h"

#ifndef M_PI
# define M_PI  3.1415926535897932384626433832795029L  /* pi */
#endif


void assert_fcs(FCSResult r)
{
  if(r) {
    fcsResult_printResult(r);
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    exit(-1);
  }
}

int main(int argc, char* argv[])
{
  MPI_Comm comm;
  int comm_size, comm_rank;
  fcs_float *x, *q, *f, *p;
  fcs_int n_axis, n_total, n_local, n_local_max;
  fcs_int i, j, k;
  fcs_int p_c, p_start, p_stop,ip;
  fcs_float e_local, e_total;
  fcs_float madelung_approx;
  const fcs_float madelung = 1.74756459463318219;
  int mpi_thread_requested = MPI_THREAD_MULTIPLE;
  int mpi_thread_provided;

  FCS fcs_handle;
  FCSResult fcs_result;

  char method[] = "pepc";
  fcs_float box_a[] = { 1.0, 0.0, 0.0 };
  fcs_float box_b[] = { 0.0, 1.0, 0.0 };
  fcs_float box_c[] = { 0.0, 0.0, 1.0 };
  fcs_float offset[] = { 0.0, 0.0, 0.0 };
  fcs_int periodic[] = { 1, 1, 1 };
  //fcs_int periodic[] = { 0, 0, 0 };

  fcs_float theta   = 0.2;
  fcs_float epsilon = 1.23e-6;


  MPI_Init_thread(&argc, &argv, mpi_thread_requested, &mpi_thread_provided);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  
  if (mpi_thread_provided < mpi_thread_requested && comm_rank == 0) {
    printf("Call to MPI_INIT_THREAD failed. Requested/provided level of multithreading: %d / %d. Continuing but expect program crash.\n", mpi_thread_requested, mpi_thread_provided);
  }
  

  n_axis = 16;
  n_total = n_axis * n_axis * n_axis;

  n_local = n_total / comm_size;
  if(comm_rank == comm_size-1) n_local += n_total % comm_size;
  n_local_max = n_total / comm_size + n_total % comm_size;


  if (comm_rank == 0) {
    printf("*** RUNNING pepc TEST ***\n");
    printf("  n_total =          %" FCS_LMOD_INT "d\n", n_total);
    printf("  n_procs =          %d\n", comm_size);
    printf("  theta =            %e\n", theta);
    printf("  epsilon =          %e\n", epsilon);
  }

  x = (fcs_float*)malloc(3 * n_local * sizeof(fcs_float));
  q = (fcs_float*)malloc(    n_local * sizeof(fcs_float));
  f = (fcs_float*)malloc(3 * n_local * sizeof(fcs_float));
  p = (fcs_float*)malloc(    n_local * sizeof(fcs_float));

  p_c = 0;
  p_start = comm_rank*(n_total/comm_size);
  p_stop  = p_start + n_local;
  for (ip=p_start; ip<p_stop; ip++, p_c++) {
    
    i = ip % n_axis;
    j = (ip / n_axis) % n_axis;
    k = ip / (n_axis*n_axis);
  
    x[3*p_c  ] = offset[0] + (i * box_a[0]) / n_axis;
    x[3*p_c+1] = offset[1] + (j * box_b[1]) / n_axis;
    x[3*p_c+2] = offset[2] + (k * box_c[2]) / n_axis;
    q[p_c] = ((i+j+k)%2 ? 1.0 : -1.0);
  
    /* printf("init positions (rank %d) for particle id %d: %e %e %e %e\n", */
    /* 	   comm_rank, ip, x[3*p_c  ], x[3*p_c+1], x[3*p_c+2], q[p_c]); */
  }

  fcs_result = fcs_init(&fcs_handle, method, comm);
  assert_fcs(fcs_result);

  fcs_result = fcs_set_common(fcs_handle, 1,
			      box_a, box_b, box_c,
			      offset, periodic, n_total);
  assert_fcs(fcs_result);

  fcs_result = fcs_pepc_setup(fcs_handle, epsilon, theta);
  assert_fcs(fcs_result);

  fcs_result = fcs_tune(fcs_handle, n_local, n_local_max, x, q);
  assert_fcs(fcs_result);

  fcs_result = fcs_run(fcs_handle, n_local, n_local_max, x, q, f, p);
  assert_fcs(fcs_result);

  e_local = 0.0;
  for (i=0; i<n_local; ++i)
    e_local += p[i] * q[i];

  MPI_Reduce(&e_local, &e_total, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  //madelung_approx = 8.0 * M_PI / n_axis * e_total / n_total;
  madelung_approx = 1.0 / n_axis * e_total / n_total;

  if (comm_rank == 0) {
    printf("\n");
    printf("  Results:\n");
    printf("    Energy:            %e\n", e_total);
    printf("    Madelung constant: %e\n", madelung_approx);
    printf("    Relative error:    %e\n", fabs(madelung-fabs(madelung_approx))/madelung);
  }

  p_c = 0;
  p_start = comm_rank*(n_total/comm_size);
  p_stop  = p_start + n_local;
  for (ip=p_start; ip<p_stop; ip++, p_c++) {
    
    /* printf("results (rank %d) for particle id %d: %e %e %e %e\n", */
    /* 	   comm_rank, ip, f[3*p_c  ], f[3*p_c+1], f[3*p_c+2], p[p_c]); */
    /* printf("dataout: %e %e %e %e %e %e %e %e\n", x[3*p_c  ], x[3*p_c+1], x[3*p_c+2], q[p_c],  */
    /* 	   f[3*p_c  ], f[3*p_c+1], f[3*p_c+2], p[p_c]); */
  }


  fcs_destroy(fcs_handle);

  free(x);
  free(q);
  free(f);
  free(p);

  if (comm_rank == 0)
    printf("*** pepc DONE ***\n");

  MPI_Finalize();

  return 0;
}
