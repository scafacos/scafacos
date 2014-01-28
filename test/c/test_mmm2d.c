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
   
  const char* method = "mmm2d";
  const char* datafile = "../inp_data/mmm/mmm2d_wall.dat";
  
  /***************************************************/
  /*
  /// SET UP MPI
  MPI_Init(&argc, &argv);

  /// Create cartesian communicator
  MPI_Comm comm;
  fcs_int node_grid[3] = {0, 0, 0};
  fcs_int node[3];
  int comm_rank, comm_size;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size); /// numero de procesos
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank); /// que proceso soy yo
  if (comm_rank == 0) {
    fprintf(stderr, "----------------\n");
    fprintf(stderr, "Running P3M test\n");
    fprintf(stderr, "----------------\n");
    fprintf(stderr, "Setting up MPI...\n");
    fprintf(stderr, "  Using %d tasks.\n", comm_size);
    printf("rank %d\n", comm_rank);
  }

  MPI_Dims_create(comm_size, 3, node_grid);
  /// swap first and last dimension, as P3M currently wants to have the
increasing 
  fcs_int tmp = node_grid[2];
  node_grid[2] = node_grid[0];
  node_grid[0] = tmp;
  
  MPI_Cart_create(MPI_COMM_WORLD, 3, node_grid, periodicity, 1, &comm);
  if (comm_rank == 0)
    fprintf(stderr, "  node_grid=(%d, %d, %d)\n", node_grid[0], node_grid[1],
node_grid[2]);

  MPI_Comm_rank(comm, &comm_rank);
  printf("comm_rank: %d\n",comm_rank);
  MPI_Cart_coords(comm, comm_rank, 3, node);
  */
  /***************************************************/
  
  FCS handle = NULL;
  FCSResult result = NULL;
  
  fcs_float box_l[3] = { 1.0, 1.0, .5 };
  fcs_float offset[3] = {0.0, 0.0, 0.0};
  fcs_int periodicity[3] = { 1, 1, 0 };
  fcs_int total_particles = 2, pid;
  fcs_int n_layers_per_node=1;
  fcs_float charges[total_particles];
  fcs_float positions[3*total_particles];
  fcs_float forces[3*total_particles];
  fcs_float potentials[total_particles];
  fcs_int i;
  for (i=0; i<total_particles; i++) {
    charges[i]=0.;
    potentials[i]=0.;
  }
  for (i=0; i<3*total_particles; i++) {
    positions[i]=0.;
    forces[i]=0.;
  }
  
  printf("Initializing MPI...\n");
  /***************************************************/
  /* SET UP MPI */
  MPI_Comm comm;
  MPI_Init(&argc, &argv);
  
  /* Create cartesian communicator */
  fcs_int node[3];
  fcs_int comm_rank, comm_size;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size); /// numero de procesos
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank); /// que proceso soy yo
  if (comm_rank == 0) {
    printf("----------------------------------------\n");
    printf("Running mmm2d test using Using %d tasks.\n", comm_size);
    printf("----------------------------------------\n");
  }
  printf("rank %d\n", comm_rank);
  
  
  
  ////////////// comm_size sera el numero maximo de layers posibles
  fcs_int node_grid[3] = {1, 1, comm_size};
  
  MPI_Dims_create(comm_size, 3, node_grid);
  MPI_Cart_create(MPI_COMM_WORLD, 3, node_grid, periodicity, 1, &comm);
  
  if (comm_rank == 0)
    fprintf(stderr, "  node_grid=(%d, %d, %d)\n", node_grid[0], node_grid[1], node_grid[2]);
  MPI_Cart_coords(comm, comm_rank, 3, node);
  /***************************************************/
  if (comm_rank == 0)
    printf("\nInitializing mmm2d...\n");
  result = fcs_init(&handle, method, comm);
  assert_fcs(result);
  
  fcs_int n_particles=0;
  if (comm_rank == 0) {
    n_particles=total_particles;
    printf("\nReading test system data...\n");
    printf("Reading %s...\n", datafile);
    FILE *data = fopen(datafile, "r");
    if (!data) {
      fprintf(stderr, "ERROR: Can't read %s!\n", datafile);
      perror("ERROR");
      exit(1);
    }
    
    fcs_float charge_sum = 0.0;
    for (pid = 0; pid < n_particles; pid++) {
      fscanf(data, "%lf %lf %lf", &positions[3*pid], &positions[3*pid+1],
&positions[3*pid+2]);
      fscanf(data, "%lf", &charges[pid]);
      charge_sum += charges[pid];
    }
    fclose(data);
  } else n_particles = 0;
  
  if (comm_rank == 0)
    printf("\nSetting parameters...\n");

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
  printf("common parameters set\n");

  assert_fcs(result);
  
  /* Parameters interface */
  fcs_float val1;
  fcs_int vali;
  
  /**/
  printf("\nTesting parameters interface...\n");
  //defaults
  fcs_mmm2d_get_far_cutoff(handle, &val1);
  printf("default far cutoff set as %f\n",val1);
  fcs_mmm2d_get_maxPWerror(handle, &val1);
  printf("default max error set as %f\n",val1);
  /*
  fcs_mmm2d_get_dielectric_contrasts(handle, &val1, &val2);
  printf("default dielectric contrasts set as top: %f bottom: %f\n",val1, val2);
  fcs_mmm2d_set_far_cutoff(handle, 1.73);
  fcs_mmm2d_get_far_cutoff(handle, &val1);
  printf("far cutoff set as %f\n",val1);
  
  fcs_mmm2d_set_dielectric_contrasts(handle, 3.17, 2.13);
  fcs_mmm2d_get_dielectric_contrasts(handle, &val1, &val2);
  printf("dielectric contrasts set as top: %f bottom: %f\n",val1, val2);
  */
  result=fcs_mmm2d_set_maxPWerror(handle, 1.e-3);
  if (result != NULL) {
    fcs_result_print_result(result);
    MPI_Finalize();
    exit(1);
  }
  result=fcs_mmm2d_get_maxPWerror(handle, &val1);
  if (result != NULL) {
    fcs_result_print_result(result);
    MPI_Finalize();
    exit(1);
  }
  if (comm_rank == 0)
    printf("max error set as %e\n",val1);
  
  
  result=fcs_mmm2d_set_layers_per_node(handle, n_layers_per_node);
  if (result != NULL) {
    fcs_result_print_result(result);
    MPI_Finalize();
    exit(1);
  }
  fcs_mmm2d_get_layers_per_node(handle, &vali);
  if (comm_rank == 0)
    printf("number of layers per node set as %d\n",vali);
  
  /*
  fcs_mmm2d_set_skin(handle, 0.5);
  fcs_mmm2d_get_skin(handle, &val1);
  printf("skin set as %f\n",val1);
  */
  
  /* Tuning */
  if (comm_rank == 0)
    printf("\nTesting tuning routine...\n");
  MPI_Barrier(comm);
  result = fcs_tune(handle, n_particles, n_particles, positions, charges);
  if (result != NULL) {
    fcs_result_print_result(result);
    MPI_Finalize();
    exit(1);
  }
  MPI_Barrier(comm);
  
  if (comm_rank == 0) {
    fcs_mmm2d_get_far_cutoff(handle, &val1);
    printf("far cutoff set at tuning as %f\n",val1);
  }
  
  MPI_Barrier(comm);
  
  result = fcs_mmm2d_require_total_energy(handle, 1);
  
  result = fcs_run(handle, n_particles, n_particles, positions, charges, forces,
NULL);
  fcs_float total_energy;
  result=fcs_mmm2d_get_total_energy(handle, &total_energy);

  assert_fcs(result);
  MPI_Barrier(comm);
  if (comm_rank == 0) {
    /*printf("Potential on particle 1: %e\n", potentials[0]);
    printf("Potential on particle 2: %e\n", potentials[1]);
    if (n_particles>=4) {
      printf("Potential on particle 3: %e\n", potentials[2]);
      printf("Potential on particle 4: %e\n", potentials[3]);
    }
    
    fcs_float toteng=0.;
    for (pid=0; pid<n_particles; pid++) {
      toteng +=potentials[pid];
    }
    printf("total energy: %e\n", toteng);
    */
    for(pid=0; pid<total_particles; pid++) {
      printf("Force on particle %d (at %e %e %e): %e, %e, %e\n", pid, positions[3*pid], positions[3*pid+1], positions[3*pid+2], forces[3*pid], forces[3*pid+1], forces[3*pid+2]);
    }
    fcs_mmm2d_get_total_energy(handle, &total_energy);
    printf("total energy: %e\n", total_energy);
    
    /*printf("Force on particle 2: %e, %e, %e\n", forces[3], forces[4],
forces[5]);
    if (n_particles>=4) {
      printf("Force on particle 3: %e, %e, %e\n", forces[6], forces[7],
forces[8]);
      printf("Force on particle 4: %e, %e, %e\n", forces[9], forces[10],
forces[11]);
    }
    */
  }
  
  fcs_destroy(handle);

  MPI_Finalize();
  if (comm_rank == 0)
    printf("Done.\n");
   
   /*
   const char* method = "mmm2d";
   MPI_Comm comm = MPI_COMM_WORLD;
   
  printf("Initializing MPI...\n");
  MPI_Init(&argc, &argv);
  
  printf("Initializing MPI...2\n");
  
  int comm_rank, comm_size;
  MPI_Comm_size(comm, &comm_size);
  printf("Initializing MPI...3\n");
  MPI_Comm_rank(comm, &comm_rank);
  printf("Initializing MPI...4\n");

  FCS handle = NULL;
  FCSOutput output = NULL;
  FCSResult result = NULL;

  
  */
   
   /*
  
  const char* datafile = "test/inp_data/mmm/mmm2d_wall.dat";
  MPI_Comm comm = MPI_COMM_WORLD;
  fcs_int periodicity[3] = { 0, 0, 1 };
  fcs_int pid;

  printf("------------------\n");
  printf("Running mmm2d test\n");
  printf("------------------\n");
  
  fcs_float box_l[3] = { 10.0, 10.0, 10.0 };
  int n_particles = 2;
  
  fcs_float charges[n_particles];
  fcs_float positions[3*n_particles];
  fcs_float forces[3*n_particles];
  fcs_float potentials[n_particles];
  
  printf("Reading %s...\n", datafile);
  FILE *data = fopen(datafile, "r");
  if (!data) {
    fprintf(stderr, "ERROR: Can't read %s!\n", datafile);
    perror("ERROR");
    exit(1);
  }
  
  fcs_float charge_sum = 0.0;
  for (pid = 0; pid < n_particles; pid++) {
    fscanf(data, "%lf %lf %lf", &positions[3*pid], &positions[3*pid+1],
&positions[3*pid+2]);
    fscanf(data, "%lf", &charges[pid]);
    charge_sum += charges[pid];
  }
  fclose(data);
  
  
  printf("Charge 1 with value %e, position %e %e
%e\n",charges[0],positions[0],positions[1],positions[2]);
  printf("Charge %d with value %e, position %e %e
%e\n",n_particles,charges[n_particles-1],positions[3*n_particles-3],positions[3*
n_particles-2],positions[3*n_particles-1]);
  
  printf("Initializing MPI...\n");
  MPI_Init(&argc, &argv);
  
  printf("Initializing MPI...2\n");
  
  int comm_rank, comm_size;
  MPI_Comm_size(comm, &comm_size);
  printf("Initializing MPI...3\n");
  MPI_Comm_rank(comm, &comm_rank);
  printf("Initializing MPI...4\n");

  FCS handle = NULL;
  FCSOutput output = NULL;
  FCSResult result = NULL;

  printf("\nInitializing MMM1D...\n");
  result = fcs_init(&handle, method, comm);
  assert_fcs(result);

  printf("\nSetting parameters...\n");

  fcs_float box_a[3] = { 0.0, 0.0, 0.0 };
  fcs_float box_b[3] = { 0.0, 0.0, 0.0 };
  fcs_float box_c[3] = { 0.0, 0.0, 0.0 };
  box_a[0] = box_l[0];
  box_b[1] = box_l[1];
  box_c[2] = box_l[2];

  result = fcs_set_common(handle, 1, box_a, box_b, box_c, 
           periodicity, n_particles, n_particles);
  if (result != NULL) {
    fcs_result_print_result(result);
    MPI_Finalize();
    exit(1);
  }

  assert_fcs(result);
  
  /// create fcs output structure
  */
  /*
  result=fcsOutput_create(&output);
  assert_fcs(result);
  */

  /* Tuning */
  /*
  printf("\nTesting parameters interface...\n");
  fcs_float val;
  fcs_int vali;
  
  fcs_MMM1D_get_far_switch_radius(handle, &val);
  printf("default far switch radius: %f\n",val);
  fcs_MMM1D_set_far_switch_radius(handle, 6.0);
  fcs_MMM1D_get_far_switch_radius(handle, &val);
  printf("new far switch radius: %f\n",val);
  
  fcs_MMM1D_get_bessel_cutoff(handle, &vali);
  printf("default bessel_cutoff: %d\n",vali);
  fcs_MMM1D_set_bessel_cutoff(handle, 3);
  fcs_MMM1D_get_bessel_cutoff(handle, &vali);
  printf("new bessel_cutoff: %d\n",vali);
  
  fcs_MMM1D_get_maxPWerror(handle, &val);
  printf("default maxPWerror: %f\n",val);
  fcs_MMM1D_set_maxPWerror(handle, 0.0001);
  fcs_MMM1D_get_maxPWerror(handle, &val);
  printf("new maxPWerror: %f\n",val);
  
  fcs_MMM1D_get_coulomb_prefactor(handle, &val);
  printf("default coulomb prefactor: %f\n",val);
  fcs_MMM1D_set_coulomb_prefactor(handle, 1.0);
  fcs_MMM1D_get_coulomb_prefactor(handle, &val);
  printf("new coulomb prefactor: %f\n",val);
  
  printf("\nTesting tuning routine...\n");
  result = fcs_tune(handle, n_particles, n_particles, positions, charges, 0);
  if (result != NULL) {
    fcs_result_print_result(result);
    MPI_Finalize();
    exit(1);
  }
  
  fcs_MMM1D_get_bessel_cutoff(handle, &vali);
  printf("calculated bessel cutoff: %d\n",vali);
  
  printf("\nTesting running routine...\n");
  result = fcs_run(handle, n_particles, n_particles, positions, charges, forces,
potentials, 0, output);
  assert_fcs(result);
  printf("Potential on particle 1: %e\n", potentials[0]);
  printf("Potential on particle 2: %e\n", potentials[1]);
  printf("Force on particle 1: %e, %e, %e\n", forces[0], forces[1], forces[2]);
  printf("Force on particle 2: %e, %e, %e\n", forces[3], forces[4], forces[5]);
  
  printf("\nFinalizing...\n");
  fcs_destroy(handle);
  fcsOutput_destroy(output);

  MPI_Finalize();
  printf("Done.\n");
   */
  return 0;
}
