/*
  Copyright (C) 2011,2012 René Halver, Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.
  
  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "fcs.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "common/near/near.h"

#define INTERFACE_TEST_BOX_SIZE 1.03125l
#undef INTERFACE_TEST_BOX_SIZE
#define INTERFACE_TEST_BOX_SIZE 1.01l
#define INTERFACE_TEST_XYZ_INTERVAL 1
#define INTERFACE_TEST_OUTPUT_INTERVAL 1
#define FCS_NEAR_TEST_LOOP_COUNT 100

/* coulomb potential / field for testing purposes */
fcs_float
fcs_near_coulomb_potential (fcs_float r_ij, fcs_float q_j, const void *param)
{
  return q_j / r_ij;
}

fcs_float
fcs_near_coulomb_field (fcs_float r_ij, fcs_float q_j, const void *param)
{
  return q_j / (r_ij * r_ij * r_ij);
}

int
main (int argc, char **argv)
{
  char **arguments = argv;
  int arg_count = argc;

  /* time values for measuring near field time */
  /* double t1 = 0.0, t2 = 0.0, tmax = 0.0; */

  char *method;
  char *filename;
  fcs_int number_of_runs;
  MPI_Comm communicator = MPI_COMM_WORLD;

  fcs_float *box_a;
  fcs_float *box_b;
  fcs_float *box_c;
  fcs_float *offset;
  fcs_float *particles;
  fcs_float *charges;
  fcs_float *velocities;
  fcs_float *masses;
  fcs_float *field;
  fcs_float *potentials;
  fcs_int periodicity[3];
  /*
  fcs_int near_field_flag = 1;
  */
#if FCS_ENABLE_DIRECT
  char DIRECT_parameters[] = "DIRECT_cutoff,0.0,DIRECT_periodic_images,1,1,1";
#endif
#if FCS_ENABLE_EWALD
  char EWALD_parameters[] = "EWALD_required_accuracy,1e-6";
#endif
#if FCS_ENABLE_FMM
  /*
  fcs_int FMM_absrel = FCS_METHOD_FMM_CUSTOM_RELATIVE;
  fcs_int FMM_dipole_correction = FCS_METHOD_FMM_STANDARD_DIPOLE_CORRECTION;
  fcs_float FMM_deltaE = 1e-6;
  */
  char FMM_parameters[] = "fmm_absrel,2,fmm_dipole_correction,0,fmm_tolerance_energy,1e-6,fmm_internal_tuning,0ll";
#endif
#if FCS_ENABLE_MMM1D
  char MMM1D_parameters[] = "MMM1D_bessel_cutoff,3,MMM1D_far_switch_radius,6.0,MMM1D_maxPWerror,1e-6";
#endif
#if FCS_ENABLE_MMM2D
  char MMM2D_parameters[] = "MMM2D_maxPWerror,1e-3,MMM2D_far_cutoff,1.73,MMM2D_dielectric_contrasts,3.17,2.13,MMM2D_layers_per_node,100,MMM2D_skin,0.5";
#endif
#if FCS_ENABLE_PEPC
  /*
  fcs_float PEPC_epsilon = 0.5;
  fcs_float PEPC_theta = 0.5;
  fcs_int PEPC_debuglevel = -1;
  */
  char PEPC_parameters[] = "PEPC_debuglevel,-1,PEPC_epsilon,0.5,PEPC_theta,0.5";
#endif
#if FCS_ENABLE_PP3MG
  /*
  int *PP3MG_dims;
  fcs_int PP3MG_cells_x = 128;  
  fcs_int PP3MG_cells_y = 128;  
  fcs_int PP3MG_cells_z = 128;  
  fcs_int PP3MG_ghost_cells = 4;
  fcs_int PP3MG_pol_degree = 4; 
  fcs_float PP3MG_err_bound = 1e-6;
  fcs_int PP3MG_max_iterations = 20;
  fcs_int PP3MG_max_particles = 20; 
  fcs_int PP3MG_periodic = 1;
  */
  char PP3MG_parameters[] = "pp3mg_cells_x,128,pp3mg_cells_y,128,pp3mg_cells_z,128,pp3mg_ghosts,4";
#endif
#if FCS_ENABLE_P2NFFT
  /* fcs_float P2NFFT_accuracy = 1e-6; */
  char P2NFFT_parameters[] = "P2NFFT_required_accuracy,1e-6";
#endif
#if FCS_ENABLE_P3M
  /* fcs_float P3M_accuracy = 1e-6; */
  char P3M_parameters[] = "P3M_tolerance_field_abs,1e-6";
#endif
#if FCS_ENABLE_VMG
  char VMG_parameters[] = "VMG_cycle_type,2,VMG_max_iterations,20,VMG_max_level,6,VMG_near_field_cells,6,VMG_precision,1e-6,VMG_smoothing_steps,3";
#endif
  FILE *input_file;
  FILE *xyz_file;
  fcs_int local_particles;
  char line[400];
  fcs_int number_of_particles;

  fcs_int i, j, k;
  int comm_rank, comm_size;
  int dims[3], pers[3];
  char c_dummy;

  FCS handle = NULL;
  FCSResult result = NULL;

/*  const fcs_float time_step = 1e-9;
  const fcs_float charge_scale = 1.602e-19;
  const fcs_float mass_scale = 1.6e-26; */
  /* scaling factor for MD steps */
  /* const fcs_float scaling_factor = time_step * charge_scale / mass_scale; */
  MPI_Comm cart_comm;

  char common_parameters[] = "box_a,1.00,0.0,0.0,box_b,0.0,1.00,0.0,box_c,0.0,0.0,1.00,periodicity,1,1,1,offset,0.0,0.0,0.0,near_field_flag,0";
  
  
  periodicity[0] = 1;
  periodicity[1] = 1;
  periodicity[2] = 1;

  dims[0] = dims[1] = dims[2] = 0;
  pers[0] = periodicity[0];
  pers[1] = periodicity[1];
  pers[2] = periodicity[2];

  /* command line handling */
  MPI_Init (&arg_count, &arguments);

  MPI_Comm_size (communicator, &comm_size);
  fprintf (stderr, "comm_size: %d \n", comm_size);
  MPI_Dims_create (comm_size, 3, dims);
  fprintf (stderr, "dims: %d %d %d\n", dims[0], dims[1], dims[2]);
  MPI_Cart_create (communicator, 3, dims, pers, 0, &cart_comm);
  communicator = cart_comm;
  fprintf (stderr, "method: %s\n", argv[1]);
  
  if (argc < 3 || argc > 4)
    {
      fprintf (stderr,
	       "Usage: %s <solver> <input_file_name> [<number_of_runs>]\n",
	       argv[0]);
      fprintf (stderr,
	       "  currently only the files in test_env/data are supported\n");
      return (-1);
    }

  method = argv[1];
  filename = argv[2];

  number_of_runs = 1;
  if (argc > 3)
    number_of_runs = atoi (argv[3]);

  MPI_Comm_rank (communicator, &comm_rank);
  MPI_Comm_size (communicator, &comm_size);

  if (comm_rank == 0)
    printf ("Interface test\n");

  input_file = fopen (filename, "r");
  if (input_file == NULL)
    {
      fprintf (stderr, "ERROR: file %s could not be opened!\n", filename);
      return 2;
    }

  for (i = 0; i < 4; ++i)
    {
      fgets (line, 400, input_file);
    }
  fscanf (input_file, "%c%" FCS_CONV_INT "d\n", &c_dummy, &number_of_particles);


  /* dividing particles up */
  if (comm_rank == comm_size - 1)
    local_particles =
      number_of_particles - (comm_size -
			     1) * (number_of_particles / comm_size);
  else
    local_particles = number_of_particles / comm_size;


  particles = (fcs_float *) malloc (3 * sizeof (fcs_float) * local_particles);
  charges = (fcs_float *) malloc (3 * sizeof (fcs_float) * local_particles);
  velocities =
    (fcs_float *) malloc (3 * sizeof (fcs_float) * local_particles);
  masses = (fcs_float *) malloc (sizeof (fcs_float) * local_particles);
  field = (fcs_float *) malloc (3 * sizeof (fcs_float) * local_particles);
  potentials = (fcs_float *) malloc (sizeof (fcs_float) * local_particles);

  /*      printf("%d particles: %d/%d\n", comm_rank, local_particles, number_of_particles); */

  j = 0;
  for (i = 0; i < number_of_particles; ++i)
    {
      if (((i >= comm_rank * local_particles)
	   || i >= number_of_particles - local_particles)
	  && i < (comm_rank + 1) * local_particles)
	{
	  fscanf (input_file,
		  "%" FCS_CONV_FLOAT "f\t %" FCS_CONV_FLOAT "f\t %" FCS_CONV_FLOAT "f\t %" FCS_CONV_FLOAT "f\t %" FCS_CONV_FLOAT "f\t %" FCS_CONV_FLOAT "f\t %" FCS_CONV_FLOAT "f\t %" FCS_CONV_FLOAT "f\n",
		  (particles + j * 3), (particles + j * 3 + 1),
		  (particles + j * 3 + 2), masses + j, (charges + j),
		  (velocities + j * 3), (velocities + j * 3 + 1),
		  (velocities + j * 3 + 2));
	  j++;
	}
      else
	{
	  fgets (line, 400, input_file);
	}
    }

  /*      for (i = 0; i < comm_size; ++i)
     {
     if (comm_rank == i)
     {
     for (j = 0; j < local_particles; ++j)
     printf("%d %e %e %e %e\n", comm_rank, *(particles+j*3), *(particles+j*3+1), *(particles+j*3+2), *(charges+j));
     }
     MPI_Barrier(communicator);
     } */

  box_a = (fcs_float *) malloc (3 * sizeof (fcs_float));
  box_b = (fcs_float *) malloc (3 * sizeof (fcs_float));
  box_c = (fcs_float *) malloc (3 * sizeof (fcs_float));
  offset = (fcs_float *) malloc (3 * sizeof (fcs_float));

  box_a[0] = box_a[1] = box_a[2] = box_b[0] = box_b[1] = box_b[2] = box_c[0] =
  box_c[1] = box_c[2] = 0.0;
  box_a[0] = INTERFACE_TEST_BOX_SIZE;
  box_b[1] = INTERFACE_TEST_BOX_SIZE;
  box_c[2] = INTERFACE_TEST_BOX_SIZE;

  offset[0] = 0.0;
  offset[1] = 0.0;
  offset[2] = 0.0;


  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Initializing %s method...\n", method);
  result = fcs_init (&handle, method, communicator);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);

  MPI_Barrier(communicator);

  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Testing generic near field ...\n");
/*
  for (i = 0; i < FCS_NEAR_TEST_LOOP_COUNT; ++i)
  {  
    t1 += MPI_Wtime();
  	fcs_near_field(communicator, box, periodicity, 0.02,  local_particles,number_of_particles,particles,charges,field,potentials,NULL,fcs_near_coulomb_potential,fcs_near_coulomb_field);
  	t2 += MPI_Wtime();
  }

  t2 -= t1;
  MPI_Reduce (&t2, &tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, communicator);
  if (comm_rank == 0)
    fprintf (stdout, "time spent in near field: %14.8e\n",
	     tmax / ((fcs_float) FCS_NEAR_TEST_LOOP_COUNT));
*/

  MPI_Barrier(communicator);         
         
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Setting common data...\n");
  /*result =
    fcs_set_common (handle, near_field_flag, box_a, box_b, box_c, offset,
		    periodicity, number_of_particles);*/
  result = fcs_set_parameters(handle, common_parameters, FCS_FALSE);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  result = fcs_set_total_particles(handle, number_of_particles);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  
  MPI_Barrier(communicator);         
  
  switch (fcs_get_method (handle))
    {
#if FCS_ENABLE_DIRECT
    case FCS_METHOD_DIRECT:
      if (comm_rank == 0)
    printf ("\n");
      if (comm_rank == 0)
    printf ("Setting up DIRECT... \n");
      result = /*
    fcs_FMM_setup (handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction);
              */
              fcs_set_parameters( handle, DIRECT_parameters, FCS_FALSE);
      if (comm_rank == 0)
    fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif
#if FCS_ENABLE_EWALD
    case FCS_METHOD_EWALD:
      if (comm_rank == 0)
    printf ("\n");
      if (comm_rank == 0)
    printf ("Setting up EWALD... \n");
      result = /*
    fcs_FMM_setup (handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction);
              */
              fcs_set_parameters( handle, EWALD_parameters, FCS_FALSE);
      if (comm_rank == 0)
    fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif
#if FCS_ENABLE_FMM
    case FCS_METHOD_FMM:
      if (comm_rank == 0)
    printf ("\n");
      if (comm_rank == 0)
    printf ("Setting up FMM... \n");
      result = /*
    fcs_FMM_setup (handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction);
              */
              fcs_set_parameters( handle, FMM_parameters, FCS_FALSE);
      if (comm_rank == 0)
    fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif
#if FCS_ENABLE_MMM1D
    case FCS_METHOD_MMM1D:
      if (comm_rank == 0)
    printf ("\n");
      if (comm_rank == 0)
    printf ("Setting up MMM1D... \n");
      result = /*
    fcs_FMM_setup (handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction);
              */
              fcs_set_parameters( handle, MMM1D_parameters, FCS_FALSE);
      if (comm_rank == 0)
    fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif
#if FCS_ENABLE_MMM2D
    case FCS_METHOD_MMM2D:
      if (comm_rank == 0)
    printf ("\n");
      if (comm_rank == 0)
    printf ("Setting up MMM2D... \n");
      result = /*
    fcs_FMM_setup (handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction);
              */
              fcs_set_parameters( handle, MMM2D_parameters, FCS_FALSE);
      if (comm_rank == 0)
    fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif
#if FCS_ENABLE_PEPC
    case FCS_METHOD_PEPC:
      if (comm_rank == 0)
	printf ("\n");
      if (comm_rank == 0)
	printf ("Setting up PEPC\n");
      result = /*
	fcs_PEPC_setup (handle, PEPC_epsilon, PEPC_theta, PEPC_debuglevel);
               */
               fcs_set_parameters(handle, PEPC_parameters, FCS_FALSE);
      if (comm_rank == 0)
	fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif
#if FCS_ENABLE_PP3MG
    case FCS_METHOD_PP3MG:
      if (comm_rank == 0)
    printf ("\n");
      if (comm_rank == 0)
    printf ("Setting up PP3MG\n");
      result = /*
    fcs_PP3MG_setup (handle, PP3MG_dims, PP3MG_cells_x, PP3MG_cells_y,
             PP3MG_cells_z, PP3MG_ghost_cells,
             PP3MG_max_iterations, PP3MG_max_particles,
             PP3MG_periodic, PP3MG_pol_degree, PP3MG_err_bound);
             */
             fcs_set_parameters (handle, PP3MG_parameters, FCS_FALSE);
      if (comm_rank == 0)
    fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif
#if FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      if (comm_rank == 0)
    printf ("\n");
      if (comm_rank == 0)
    printf ("Setting up P2NFFT... \n");
      result = /*
    fcs_FMM_setup (handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction);
              */
              fcs_set_parameters( handle, P2NFFT_parameters, FCS_FALSE);
      if (comm_rank == 0)
    fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif
#if FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      if (comm_rank == 0)
    printf ("\n");
      if (comm_rank == 0)
    printf ("Setting up P3M... \n");
      result = /*
    fcs_FMM_setup (handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction);
              */
              fcs_set_parameters( handle, P3M_parameters, FCS_FALSE);
      if (comm_rank == 0)
    fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif
#if FCS_ENABLE_VMG
    case FCS_METHOD_VMG:
      if (comm_rank == 0)
    printf ("\n");
      if (comm_rank == 0)
    printf ("Setting up VMG... \n");
      result = /*
    fcs_FMM_setup (handle, FMM_absrel, FMM_deltaE, FMM_dipole_correction);
              */
              fcs_set_parameters( handle, VMG_parameters, FCS_FALSE);
      if (comm_rank == 0)
    fcs_result_print_result (result);
      fcs_result_destroy (result);
      break;
#endif

    default:
      break;

    }

  MPI_Barrier(communicator);         

  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters... \n");
  if (comm_rank == 0)
    fcs_print_parameters (handle);

  MPI_Barrier(communicator);         
  
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Tuning... \n");
  result =
    fcs_tune (handle, local_particles, particles, charges);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);

  MPI_Barrier(communicator);         
  
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Running... \n");
  for (i = 0; i < number_of_runs; ++i)
    {
      /*if(comm_rank == 0) printf("\n"); */
      /*if(comm_rank == 0) printf("Run %d... \n", i+1); */
      result = fcs_run (handle, local_particles, particles,
			charges, field, potentials);
      /*if(comm_rank == 0) fcs_result_print_result(result);
         if(comm_rank == 0) printf("\n");
         if(comm_rank == 0) printf("Moving particles: %d... \n", i+1); */
/*      
      for (j = 0; j < local_particles; ++j)
	{
	  for (k = 0; k < 3; ++k)
	    {
	      *(velocities + j * 3 + k) +=
		*(field + j * 3 + k) * *(charges + j) / *(masses +
							  j) * scaling_factor;
	      *(particles + j * 3 + k) +=
		*(velocities + j * 3 + k) * time_step;
	      if (*(periodicity + k) == 1)
		{
		  if (*(particles + j * 3 + k) > INTERFACE_TEST_BOX_SIZE)
		    *(particles + j * 3 + k) -= INTERFACE_TEST_BOX_SIZE;
		  else if (*(particles + j * 3 + k) < 0.0l)
		    *(particles + j * 3 + k) += INTERFACE_TEST_BOX_SIZE;
		}
	    }
	}
*/
      if ((i + 1) % INTERFACE_TEST_XYZ_INTERVAL == 0)
	{
	  for (j = 0; j < comm_size; ++j)
	    {
	      if (j == comm_rank)
		{
		  if ((j == 0) && ((i + 1) == INTERFACE_TEST_XYZ_INTERVAL))
		    xyz_file = fopen ("C_test.xyz", "w");
		  else
		    xyz_file = fopen ("C_test.xyz", "a");
		  if (j == 0)
		    {
		      fprintf (xyz_file, "%" FCS_LMOD_INT "d\n%s %s %d %" FCS_LMOD_INT "d %c %c %c\n",
			       number_of_particles, "C", method, comm_size,
			       number_of_runs,
			       (*periodicity == 1) ? 'T' : 'F',
			       (*(periodicity + 1) == 1) ? 'T' : 'F',
			       (*(periodicity + 2) == 1) ? 'T' : 'F');
		    }
		  for (k = 0; k < local_particles; ++k)
		    {
		      fprintf (xyz_file, "O %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e\n",
			       *(particles + 3 * k), *(particles + 3 * k + 1),
			       *(particles + 3 * k + 2), /*
			       *(charges + k),
			       *(velocities + 3 * k),
			       *(velocities + 3 * k + 1),
			       *(velocities + 3 * k + 2)); */
			       *(potentials + k),
                   *(field + 3 * k),
                   *(field + 3 * k + 1),
                   *(field + 3 * k + 2));
		    }
		  fclose (xyz_file);

		}
	      MPI_Barrier (communicator);
	    }

	}

    }

  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters... \n");
  if (comm_rank == 0)
    fcs_print_parameters (handle);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);

  free (box_a);
  free (box_b);
  free (box_c);
  free (offset);

  fcs_destroy (handle);
  /*  printf("address of handle: %p\n",*fcs_get_box_a(handle)); */

  free (particles);
  free (charges);
  free (velocities);
  free (field);
  free (potentials);
  free (masses);

  MPI_Finalize ();


  return (FCS_SUCCESS);
}
