#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
/* #include <stdlib.h> */
/* #include <math.h> */

/* #include <mpi.h> */

#include <stdlib.h>
#include "fcs.h"
void assert_fcs(FCSResult r)
{
  if (r) {
    fcs_result_print_result(r);
    MPI_Finalize();
    exit(-1);
  }
}

int main(int argc, char **argv)
{
  const char* method = "mmm1d";
  const char* datafile = "../inp_data/mmm/mmm1d_wall.dat";
  MPI_Comm comm = MPI_COMM_WORLD;
  fcs_int comm_rank, comm_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  
  fcs_int periodicity[3] = { 0, 0, 1 };
  fcs_int node_grid[3] = {1, 1, comm_size};
  fcs_int node[3];
  
  if (comm_rank == 0)
    fprintf(stderr, "Creating cartesian communicator...\n");
  MPI_Dims_create(comm_size, 3, node_grid);
  MPI_Cart_create(MPI_COMM_WORLD, 3, node_grid, periodicity, 1, &comm);
  
  if (comm_rank == 0)
    fprintf(stderr, "  node_grid=(%d, %d, %d)\n", node_grid[0], node_grid[1], node_grid[2]);
  MPI_Cart_coords(comm, comm_rank, 3, node);
  
  fcs_int pid;

  if (comm_rank == 0) {
    printf("------------------------------\n");
    printf("Running mmm1d test on %d nodes\n", comm_size);
    printf("------------------------------\n");
  }
  
  fcs_float box_l[3] = { 10.0, 10.0, 10.0 };
  fcs_float offset[3] = {0.0, 0.0, 0.0};
  fcs_int n_particles=0, total_particles = 4;
  
  fcs_float charges[total_particles];
  fcs_float positions[3*total_particles];
  fcs_float forces[3*total_particles];
  fcs_float potentials[total_particles];
  
  fcs_float charge_sum = 0.0;
  if (comm_rank == 0) {
    printf("Reading %s...\n", datafile);
    FILE *data = fopen(datafile, "r");
    if (!data) {
      fprintf(stderr, "ERROR: Can't read %s!\n", datafile);
      perror("ERROR");
      exit(1);
    }
    for (pid = 0; pid < total_particles; pid++) {
      fscanf(data, "%lf %lf %lf", &positions[3*pid], &positions[3*pid+1], &positions[3*pid+2]);
      fscanf(data, "%lf", &charges[pid]);
      printf("read: %d %le %le %le %le\n", pid, charges[pid], positions[3*pid], positions[3*pid+1], positions[3*pid+2]);
      charge_sum += charges[pid];
    }
    fclose(data);
    n_particles=total_particles;
    /*
    printf("Charge 1 with value %e, position %e %e %e\n",charges[0],positions[0],positions[1],positions[2]);
    printf("Charge %d with value %e, position %e %e %e\n", n_particles,charges[n_particles-1],positions[3*n_particles-3],positions[3*n_particles-2],positions[3*n_particles-1]);
    */
  } else n_particles=0;

  FCS handle = NULL;
  FCSResult result = NULL;

  if (comm_rank == 0) printf("\nInitializing mmm1d...\n");
  result = fcs_init(&handle, method, comm);
  assert_fcs(result);

  if (comm_rank == 0) printf("\nSetting parameters...\n");

  fcs_float box_a[3] = { 0.0, 0.0, 0.0 };
  fcs_float box_b[3] = { 0.0, 0.0, 0.0 };
  fcs_float box_c[3] = { 0.0, 0.0, 0.0 };
  box_a[0] = box_l[0];
  box_b[1] = box_l[1];
  box_c[2] = box_l[2];

  result = fcs_set_common(handle, 1, box_a, box_b, box_c, 
           offset, periodicity, total_particles);
  if (result != NULL) {
    fcs_result_print_result(result);
    MPI_Finalize();
    exit(1);
  }

  assert_fcs(result);
  
  /* Tuning */
  if (comm_rank == 0) printf("\nTesting parameters interface...\n");
  fcs_float val;
  fcs_int vali;
  
  fcs_mmm1d_get_far_switch_radius(handle, &val);
  if (comm_rank == 0) printf("default far switch radius: %f\n",val);
  fcs_mmm1d_set_far_switch_radius(handle, 6.0);
  fcs_mmm1d_get_far_switch_radius(handle, &val);
  if (comm_rank == 0) printf("new far switch radius: %f\n",val);
  
  fcs_mmm1d_get_bessel_cutoff(handle, &vali);
  if (comm_rank == 0) printf("default bessel_cutoff: %d\n",vali);
  fcs_mmm1d_set_bessel_cutoff(handle, 3);
  fcs_mmm1d_get_bessel_cutoff(handle, &vali);
  if (comm_rank == 0) printf("new bessel_cutoff: %d\n",vali);
  
  fcs_mmm1d_get_maxPWerror(handle, &val);
  if (comm_rank == 0) printf("default maxPWerror: %f\n",val);
  fcs_mmm1d_set_maxPWerror(handle, 0.0001);
  fcs_mmm1d_get_maxPWerror(handle, &val);
  if (comm_rank == 0) printf("new maxPWerror: %f\n",val);
  
  if (comm_rank == 0) printf("\nTesting tuning routine...\n");
  result = fcs_tune(handle, n_particles, n_particles, positions, charges);
  if (result != NULL) {
    fcs_result_print_result(result);
    MPI_Finalize();
    exit(1);
  }
  
  //fcs_mmm1d_get_bessel_cutoff(handle, &vali);
  //printf("calculated bessel cutoff: %d\n",vali);
  
  if (comm_rank == 0) printf("\nTesting running routine...\n");
  result = fcs_run(handle, n_particles, n_particles, positions, charges, forces, potentials);
  assert_fcs(result);
  MPI_Barrier(comm);
  if (comm_rank == 0) {
    for(pid=0; pid<total_particles; pid++) {
      printf("Potential on particle %d: %e\n", pid, potentials[pid]);
    //printf("Potential on particle 2: %e\n", potentials[1]);
      printf("Force on particle %d: %e, %e, %e\n", pid, forces[3*pid], forces[3*pid+1], forces[3*pid+2]);
    //printf("Force on particle 2: %e, %e, %e\n", forces[3], forces[4], forces[5]);
    }
    printf("\nFinalizing...\n");
  }
  fcs_destroy(handle);

  MPI_Finalize();
  if (comm_rank == 0) printf("Done.\n");

  return 0;
}
