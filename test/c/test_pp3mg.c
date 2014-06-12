#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "fcs.h"

#define DIM 3
#define PI 3.14159265

#define ASSERT_FCS(err) \
  do { \
    if(err) { \
      fcs_result_print_result(err); MPI_Finalize(); exit(-1); \
    } \
  } while (0)

/**
 * read input data
 */
int read_config_file(char *conf_file_name,
		     char *input_file_name,
		     int *m,
		     int *n,
		     int *o,
		     int *max_iterations,
		     int *max_particles,
		     int* nu1,
		     int *nu2,
		     int *ghosts,
		     double *err_bound){
  int help_size = 300;
  char help[help_size];
  FILE *conf_file;

  if((conf_file = fopen(conf_file_name, "r")) == NULL) {
    printf("No config file was found!\n");
    return -1;
  }


  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);

  sscanf(help, " %d %d %d", m, n, o);

  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);

  sscanf(help, " %d", max_iterations);


  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);

  sscanf(help, " %d", nu1);
  sscanf(help, " %d", nu2);

  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);

  sscanf(help, " %lf", err_bound);

  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);

  sscanf(help, " %d", ghosts);

  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);

  sscanf(help, " %d", max_particles);

  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);
  fgets(help, help_size, conf_file);

  sscanf(help, " %s", input_file_name);

  fclose(conf_file);

  return 0;

}


int main (int argc, char **argv)
{
  char input_file_name[300], *conf_file_name;
  FILE *input_file, *conf_file;
  int i, j, help, p_local = 0;
  FCSResult fcs_result;
  FCS fcs_handle;
  char parameterstring[200];

  int particles_i;
  fcs_int particles = 0;
  fcs_int local_particles = 0;
  fcs_float *pos = NULL;
  fcs_float *charges = NULL;
  fcs_float *field = NULL;
  fcs_float *potentials = NULL;

  int cells_x = 128, cells_y = 128, cells_z = 128;
  int periodic = 1;
  int ghosts = 4;
  int degree = 4;
  int max_particles_i;
  fcs_int max_particles = 2000000;
  int max_iterations;
  double err_bound = 1.e-3;
  int nu1, nu2;

  fcs_float box[DIM][DIM] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}};
  fcs_float offset[DIM] = {0.0, 0.0, 0.0};
  fcs_int periodic_flags[DIM] = {1,1,1};

  /* Variables for analyzing results */
  int read_ret_value = 0;

  double f_sum_local[DIM];
  double f_sum[DIM];
  double f_sum_squared_local[DIM];
  double f_sum_squared[DIM];
  double f_max_abs_local[DIM];
  double f_max_abs[DIM];

  double e_sum_local = 0.0;
  double e_sum = 0.0;
  int min = 0;

  /* Variables for reading data */
  double my_x, my_y, my_z, my_q;

  /* Size of local domain */
  double x_start, y_start, z_start;
  double x_end, y_end, z_end;

  int my_rank;
  int mpi_size;
  MPI_Comm mpi_comm_cart;
  int mpi_dims_i[DIM];
  fcs_int mpi_dims[DIM];
  int mpi_periods[DIM];
  int mpi_coords[DIM];
  MPI_Status status;

  double starttime = 0.0, endtime = 0.0;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    fprintf(stderr, "----------------\n");
    fprintf(stderr, "Running pp3mg test\n");
    fprintf(stderr, "----------------\n");
    fprintf(stderr, "Setting up MPI...\n");
    fprintf(stderr, "  Using %d tasks.\n", mpi_size);
  }

  if(argc == 1) {
    printf("No config file was specified!\n");
    MPI_Finalize();
    exit(-1);
  }

  /* read config and input files */
  if(argc >= 2) {
    conf_file_name = argv[1];
    read_ret_value = read_config_file(conf_file_name,
      input_file_name,
      &cells_x, &cells_y, &cells_z,
      &max_iterations, &max_particles_i,
      &nu1, &nu2, &ghosts, &err_bound);

    if(read_ret_value) {
      printf("Config file couldn't be read!\n");
      MPI_Finalize();
      exit(-1);
    }
    
    max_particles = max_particles_i;
  }

  if((input_file = fopen(input_file_name, "r")) == NULL) {
    printf("input file name: %s\n" , input_file_name);
    printf("No input file was found!\n");
    MPI_Finalize();
    exit(-1);
  }

  fscanf(input_file, " %d", &particles_i);
  particles = particles_i;

  /* create a cartesian communicator */
  for(i = 0; i < DIM; ++i) {
    mpi_periods[i] = 1;
    mpi_dims_i[i] = 0;
  }

  assert(MPI_Dims_create(mpi_size, 3, mpi_dims_i) == MPI_SUCCESS);
  assert(MPI_Cart_create(MPI_COMM_WORLD, 3, mpi_dims_i, mpi_periods, 1, &mpi_comm_cart) == MPI_SUCCESS);
  assert(MPI_Comm_rank(mpi_comm_cart, &my_rank) == MPI_SUCCESS);
  assert(MPI_Cart_coords(mpi_comm_cart, my_rank, 3, mpi_coords) == MPI_SUCCESS);

  for (i = 0; i < DIM; ++i) mpi_dims[i] = mpi_dims_i[i];

  if (mpi_comm_cart != MPI_COMM_NULL) {

    /* initialize parameters for particle grid */
    x_start = ((double)1.)/mpi_dims[0] * (double)mpi_coords[0];
    x_end = ((double)1.)/mpi_dims[0] * (double)(mpi_coords[0]+1.);
    y_start = ((double)1.)/mpi_dims[1] * (double)mpi_coords[1];
    y_end = ((double)1.)/mpi_dims[1] * (double)(mpi_coords[1]+1.);
    z_start = ((double)1.)/mpi_dims[2] * (double)mpi_coords[2];
    z_end = ((double)1.)/mpi_dims[2] * (double)(mpi_coords[2]+1.);

    local_particles = 0;

    /* count local partciles , domain decomposition */
    for(i = 0; i < particles; ++i) {
      fscanf(input_file, " %lf", &my_x);
      fscanf(input_file, " %lf", &my_y);
      fscanf(input_file, " %lf", &my_z);
      fscanf(input_file, " %lf", &my_q);

      if(x_start <= my_x && my_x < x_end &&
	  y_start <= my_y && my_y < y_end &&
	  z_start <= my_z && my_z < z_end)
	local_particles ++;
    }
  }
  fclose(input_file);

  if(local_particles != 0) {

    pos = malloc(local_particles*DIM*sizeof(fcs_float));
    assert(pos);
    charges = malloc(local_particles*sizeof(fcs_float));
    assert(charges);
    field = malloc(local_particles*DIM*sizeof(fcs_float));
    assert(field);
    potentials = malloc(local_particles*sizeof(fcs_float));
    assert(potentials);

    /* open the input file again and read particle data */
    if((input_file = fopen(input_file_name, "r")) == NULL) {
      printf("No input file found!\n");
      MPI_Finalize();
      exit(-1);
    }
    fscanf(input_file, " %d", &particles);

    p_local = 0;

    for(i = 0; i < particles; ++i ) {
      fscanf(input_file, " %lf", &my_x);
      fscanf(input_file, " %lf", &my_y);
      fscanf(input_file, " %lf", &my_z);
      fscanf(input_file, " %lf", &my_q);

      if(x_start <= my_x && my_x < x_end &&
	  y_start <= my_y && my_y < y_end &&
	  z_start <= my_z && my_z < z_end) {

	pos[p_local] = my_x;
	pos[p_local + 1] = my_y;
	pos[p_local + 2] = my_z;
	charges[p_local/3] = my_q;
	p_local += 3;
      }


    }
  }
  /* create FCSParameter object. IMPORTANT: pp3mg requires a cartesian MPI communicator
     which has to be created by the calling program  */
  fcs_result = fcs_init(&fcs_handle, "pp3mg", mpi_comm_cart);
  ASSERT_FCS(fcs_result);
  fcs_result = fcs_set_dimension (fcs_handle, DIM);
  ASSERT_FCS(fcs_result);
  fcs_result = fcs_set_common(fcs_handle, 1,
                                    box[0], box[1], box[2], offset, periodic_flags,
                                    particles);
  ASSERT_FCS(fcs_result);
  fcs_result = fcs_pp3mg_setup (fcs_handle, mpi_dims, cells_x, cells_y,
                                      cells_z, ghosts, max_iterations, max_particles,
                                      periodic_flags[0] /* FIXME */, degree, err_bound);
  ASSERT_FCS(fcs_result);
  fcs_result = fcs_tune(fcs_handle, local_particles, pos, charges);
  ASSERT_FCS(fcs_result);

  MPI_Barrier(mpi_comm_cart);
  if(my_rank == 0) {
    starttime = MPI_Wtime();
  }

  /* 2. step: run pp3mg. IMPORTANT: domain decomposition, particles must be distributed */
  /*          according to the specified MPI communicator */
  fcs_result = fcs_run(fcs_handle, local_particles, pos, charges, field, potentials);
  ASSERT_FCS(fcs_result);

  /* analyze results */
  /* field = fcsOutput_getField(fcs_output);
  potentials = fcsOutput_getPotentials(fcs_output); */

  MPI_Barrier(mpi_comm_cart);
  if(my_rank == 0) {
    endtime = MPI_Wtime();
    printf("pp3mg runtime with %d processors: %f s\n", mpi_size, endtime - starttime);
  }

  e_sum_local = 0.0;
  for(i = 0; i < local_particles; ++i)
    e_sum_local += potentials[i];

  e_sum = 0.0;
  MPI_Reduce(&e_sum_local,
      &e_sum,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      0,
      mpi_comm_cart);

  if(my_rank == 0) {
    min = cells_x;
    if(cells_y < min || cells_z < min)
      min = (cells_y < cells_z) ? cells_y : cells_z;

    printf("Self energy:  %e\n",
	1.0/(4.0 * PI) * 14.0/(5.0*1./(2.*ghosts*min)));
    printf("Approx. Madelung's constant: %e\n", e_sum/particles *2.0 * PI);
    printf( "Total energy:: %e\n\n", e_sum);
  }

  for(i = 0; i < DIM; ++i) {
    f_sum_local[i] = 0.0;
    f_sum[i] = 0.0;
    f_sum_squared_local[i] = 0.0;
    f_sum_squared[i] = 0.0;
    f_max_abs_local[i] = 0.0;
    f_max_abs[i] = 0.0;
  }



  for(i = 0; i < local_particles; ++i) {
    for(j = 0; j < DIM; ++j) {
      f_sum_local[j]  += charges[i] * field[i*DIM + j];
      f_sum_squared_local[j] += charges[i] * charges[i] * field[i*DIM + j]*field[i*DIM + j];
      if(fabs(field[i*DIM + j]) > f_max_abs_local[j]) {
	f_max_abs_local[j] = fabs(charges[i] * field[i*DIM + j]);
      }/* end if */
    }/* end for-j */
  }/* end for-i */

  for(j = 0; j < DIM; ++j) {
    f_max_abs[j] = f_max_abs_local[j];
  }
  MPI_Reduce(f_sum_local, f_sum, 3, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_cart);
  MPI_Reduce(f_sum_squared_local, f_sum_squared, 3,
      MPI_DOUBLE, MPI_SUM, 0, mpi_comm_cart);
  MPI_Reduce(f_max_abs_local, f_max_abs, 3,
      MPI_DOUBLE, MPI_MAX, 0, mpi_comm_cart);
  if(my_rank == 0) {
    printf("Norm of sum of forces: %e\n",
	sqrt(f_sum[0]*f_sum[0] + f_sum[1]*f_sum[1] + f_sum[2]*f_sum[2]) );

    printf("Sqrt. of sum of squares of all components of forces: %f\n",
	sqrt(f_sum_squared[0] + f_sum_squared[1] + f_sum_squared[2]));

    printf("Maximal absolute force components: \n");
    for(i = 0; i < DIM; ++i) {
      printf("%f ", f_max_abs[i]);
    }
    printf("\n");

  }

  if(local_particles != 0)
    fclose(input_file);

  /* 3. step: deallocate resources for pp3mg */
  fcs_destroy(fcs_handle);
  free (pos);
  free (charges);
  free (field);
  free (potentials);
  MPI_Comm_free(&mpi_comm_cart);

  MPI_Finalize();
  exit(0);
 
}
