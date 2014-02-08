/**
 * @file   test_vmg_open.c
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Fri Apr 15 16:30:26 2011
 *
 * @brief  Brief test of the vmg lib and interface.
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
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif

#include "fcs.h"

#ifndef M_PI
# define M_PI  3.1415926535897932384626433832795029L  /* pi */
#endif


void assert_fcs(FCSResult r)
{
  if(r) {
    fcs_result_print_result(r);
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
  fcs_int p_c;

  FCS fcs_handle;
  FCSResult fcs_result;

  char method[] = "vmg";
  fcs_float box_a[] = { 1.0, 0.0, 0.0 };
  fcs_float box_b[] = { 0.0, 1.0, 0.0 };
  fcs_float box_c[] = { 0.0, 0.0, 1.0 };
  fcs_float offset[] = { 0.0, 0.0, 0.0 };
  fcs_int periodic[] = { 0, 0, 0 };
  fcs_int max_level = 6;
  fcs_int max_iterations = 20;
  fcs_int smoothing_steps = 3;
  fcs_int cycle_type = 1;
  fcs_float precision = 1e-10;
  fcs_int near_field_cells = 6;
  fcs_int interpolation_degree = 4;
  fcs_int discretization_order = 2;

  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  n_axis = 16;
  n_total = n_axis * n_axis * n_axis;

  if (comm_rank == 0) {
    n_local = n_total;
    n_local_max = n_total;
  }else {
    n_local = 0;
    n_local_max = 0;
  }

  if (comm_rank == 0) {
    printf("*** RUNNING vmg TEST ***\n");
    printf("  n_total =              %" FCS_LMOD_INT "d\n", n_total);
    printf("  n_procs =              %d\n", comm_size);
    printf("  periodicity_x =        %d\n", periodic[0]);
    printf("  periodicity_y =        %d\n", periodic[1]);
    printf("  periodicity_z =        %d\n", periodic[2]);
    printf("  max_level =            %d\n", max_level);
    printf("  max_iterations =       %d\n", max_iterations);
    printf("  smoothing_steps =      %d\n", smoothing_steps);
    printf("  cycle_type =                %d\n", cycle_type);
    printf("  precision =            %e\n", precision);
    printf("  near_field_cells =     %d\n", near_field_cells);
    printf("  interpolation_degree = %d\n", interpolation_degree);
    printf("  discretization_order = %d\n", discretization_order);
  }

  x = (fcs_float*)malloc(3 * n_local_max * sizeof(fcs_float));
  q = (fcs_float*)malloc(n_local_max * sizeof(fcs_float));
  f = (fcs_float*)malloc(3 * n_local_max * sizeof(fcs_float));
  p = (fcs_float*)malloc(n_local_max * sizeof(fcs_float));

  if (comm_rank == 0) {

    p_c = 0;

    for (i=0; i<n_axis; ++i)
      for (j=0; j<n_axis; ++j)
	for (k=0; k<n_axis; ++k) {
	  x[3*p_c  ] = offset[0] + (i * box_a[0]) / n_axis;
	  x[3*p_c+1] = offset[1] + (j * box_b[1]) / n_axis;
	  x[3*p_c+2] = offset[2] + (k * box_c[2]) / n_axis;
	  q[p_c] = ((i+j+k)%2 ? 1.0 : -1.0);
	  ++p_c;
	}

  }

  fcs_result = fcs_init(&fcs_handle, method, comm);
  assert_fcs(fcs_result);

  fcs_result = fcs_set_common(fcs_handle, 1,
			      box_a, box_b, box_c,
			      offset, periodic, n_total);
  assert_fcs(fcs_result);

  fcs_result = fcs_vmg_setup(fcs_handle, max_level,
			     max_iterations, smoothing_steps,
			     cycle_type, precision, near_field_cells,
                             interpolation_degree, discretization_order);
  assert_fcs(fcs_result);

  fcs_result = fcs_set_max_local_particles(fcs_handle, n_local_max);
  assert_fcs(fcs_result);

  fcs_result = fcs_tune(fcs_handle, n_local, x, q);
  assert_fcs(fcs_result);

  fcs_result = fcs_run(fcs_handle, n_local, x, q, f, p);
  assert_fcs(fcs_result);

  fcs_destroy(fcs_handle);

  free(x);
  free(q);
  free(f);
  free(p);

  if (comm_rank == 0)
    printf("*** vmg DONE ***\n");

  MPI_Finalize();

  return 0;
}
