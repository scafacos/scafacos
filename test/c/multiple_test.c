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

#define INTERFACE_TEST_BOX_SIZE 1.0l
#define INTERFACE_TEST_XYZ_INTERVAL 1
#define INTERFACE_TEST_OUTPUT_INTERVAL 1
#define FCS_NEAR_TEST_LOOP_COUNT 100

int main (int argc, char **argv)
{
  char **arguments = argv;
  int arg_count = argc;

  /* time values for measuring near field time */
  /* double t1 = 0.0, t2 = 0.0, tmax = 0.0; */

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
#if FCS_ENABLE_FMM
  char FMM_parameters[] = "fmm_absrel,2,fmm_dipole_correction,0,fmm_tolerance_energy,1e-6,fmm_internal_tuning,0ll";
#endif
#if FCS_ENABLE_PP3MG
  char PP3MG_parameters[] = "pp3mg_cells_x,128,pp3mg_cells_y,128,pp3mg_cells_z,128,pp3mg_ghosts,4";
#endif
#if FCS_ENABLE_P2NFFT
  char P2NFFT_parameters[] = "P2NFFT_required_accuracy,1e-6";
#endif
#if FCS_ENABLE_P3M
  char P3M_parameters[] = "P3M_tolerance_field_abs,1e-6";
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

  /*different handles for different methods*/
#if FCS_ENABLE_FMM  
  FCS handle_fmm    = NULL;
#endif
#if FCS_ENABLE_MEMD
  FCS handle_memd = NULL;
#endif
#if FCS_ENABLE_PEPC
  FCS handle_pepc = NULL;
#endif
#if FCS_ENABLE_P2NFFT
  FCS handle_p2nfft = NULL;
#endif
#if FCS_ENABLE_P3M
  FCS handle_p3m    = NULL;
#endif
#if FCS_ENABLE_PP3MG
  FCS handle_pp3mg = NULL;
#endif
#if FCS_ENABLE_VMG
  FCS handle_vmg = NULL;
#endif
  FCSResult result = NULL;

  const fcs_float time_step = 1e-7;
  const fcs_float charge_scale = 1.602e-19;
  const fcs_float mass_scale = 1.6e-26; 
  /* scaling factor for MD steps */
  const fcs_float scaling_factor = time_step * charge_scale / mass_scale;
  MPI_Comm cart_comm;

  char common_parameters[] = "box_a,1.0,0.0,0.0,box_b,0.0,1.0,0.0,box_c,0.0,0.0,1.0,periodicity,1,1,1,offset,0.0,0.0,0.0,near_field_flag,1";
  
  
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
  
  if (argc < 2 || argc > 3)
    {
      fprintf (stderr,
         "Usage: %s <solver> <input_file_name> [<number_of_runs>]\n",
         argv[0]);
      fprintf (stderr,
         "  currently only the files in test_env/data are supported\n");
      return (-1);
    }

  filename = argv[1];

  number_of_runs = 1;
  if (argc > 2)
    number_of_runs = atoi (argv[2]);

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
#if FCS_ENABLE_PP3MG
  if (comm_rank == 0)
    printf ("Initializing %s method...\n", "pp3mg");
  result = fcs_init (&handle_pp3mg, "pp3mg", communicator);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_FMM
  if (comm_rank == 0)
    printf ("Initializing %s method...\n", "fmm");
  result = fcs_init (&handle_fmm, "fmm", communicator);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_MEMD
  if (comm_rank == 0)
    printf ("Initializing %s method...\n", "memd");
  result = fcs_init (&handle_memd, "memd", communicator);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_PEPC
  if (comm_rank == 0)
    printf ("Initializing %s method...\n", "pepc");
  result = fcs_init (&handle_pepc, "pepc", communicator);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_P2NFFT
  if (comm_rank == 0)
    printf ("Initializing %s method...\n", "p2nfft");
  result = fcs_init (&handle_p2nfft, "p2nfft", communicator);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_P3M
  if (comm_rank == 0)
    printf ("Initializing %s method...\n", "p3m");
  result = fcs_init (&handle_p3m, "p3m", communicator);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_VMG
  if (comm_rank == 0)
    printf ("Initializing %s method...\n", "vmg");
  result = fcs_init (&handle_vmg, "vmg", communicator);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif

  MPI_Barrier(communicator);

  if (comm_rank == 0)
    printf ("\n");

  MPI_Barrier(communicator);         

#if FCS_ENABLE_PP3MG
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Setting common data for pp3mg...\n");
  result = fcs_set_parameters(handle_pp3mg, common_parameters, FCS_FALSE);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  result = fcs_set_total_particles(handle_pp3mg, number_of_particles);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters of pp3mg... \n");
  if (comm_rank == 0)
    fcs_print_parameters (handle_pp3mg);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Tuning pp3mg... \n");
  result =
    fcs_tune (handle_pp3mg, local_particles, local_particles, particles, charges);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_MEMD
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Setting common data for memd...\n");
  result = fcs_set_parameters(handle_memd, common_parameters, FCS_FALSE);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  result = fcs_set_total_particles(handle_memd, number_of_particles);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters of memd... \n");
  if (comm_rank == 0)
    fcs_print_parameters (handle_memd);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Tuning memd... \n");
  result =
    fcs_tune (handle_memd, local_particles, local_particles, particles, charges);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_PEPC
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Setting common data for pepc...\n");
  result = fcs_set_parameters(handle_pepc, common_parameters, FCS_FALSE);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  result = fcs_set_total_particles(handle_pepc, number_of_particles);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters of pepc... \n");
  if (comm_rank == 0)
    fcs_print_parameters (handle_pepc);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Tuning pepc... \n");
  result =
    fcs_tune (handle_pepc, local_particles, local_particles, particles, charges);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_FMM
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Setting common data for fmm...\n");
  result = fcs_set_parameters(handle_fmm, common_parameters, FCS_FALSE);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  result = fcs_set_total_particles(handle_fmm, number_of_particles);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters of fmm... \n");
  if (comm_rank == 0)
    fcs_print_parameters (handle_fmm);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Tuning fmm... \n");
  result =
    fcs_tune (handle_fmm, local_particles, local_particles, particles, charges);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_P2NFFT
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Setting common data for p2nfft...\n");
  result = fcs_set_parameters(handle_p2nfft, common_parameters, FCS_FALSE);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  result = fcs_set_total_particles(handle_p2nfft, number_of_particles);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters of p2nfft ... \n");
  if (comm_rank == 0)
    fcs_print_parameters (handle_p2nfft);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Tuning p2nfft... \n");
  result =
    fcs_tune (handle_p2nfft, local_particles, local_particles, particles, charges);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_P3M
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Setting common data for p3m...\n");
  result = fcs_set_parameters(handle_p3m, common_parameters, FCS_FALSE);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  result = fcs_set_total_particles(handle_p3m, number_of_particles);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters of p3m ... \n");
  if (comm_rank == 0)
    fcs_print_parameters (handle_p3m);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Tuning p3m... \n");
  result =
    fcs_tune (handle_p3m, local_particles, local_particles, particles, charges);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif
#if FCS_ENABLE_VMG
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Setting common data for vmg...\n");
  result = fcs_set_parameters(handle_vmg, common_parameters, FCS_FALSE);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  result = fcs_set_total_particles(handle_vmg, number_of_particles);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Parameters of vmg... \n");
  if (comm_rank == 0)
    fcs_print_parameters (handle_vmg);
  if (comm_rank == 0)
    printf ("\n");
  if (comm_rank == 0)
    printf ("Tuning vmg... \n");
  result =
    fcs_tune (handle_vmg, local_particles, local_particles, particles, charges);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);
#endif

  MPI_Barrier(communicator);         
  
 for (i = 0; i < number_of_runs; ++i)
    {
      switch(i%6)
      {
#if FCS_ENABLE_PP3MG
        case 0:
          if (comm_rank == 0)
            printf ("\n");
          if (comm_rank == 0)
            printf ("Running pp3mg ... \n");
          result = fcs_run (handle_pp3mg, local_particles, local_particles, particles,
                            charges, field, potentials);
          break;
#endif
#if FCS_ENABLE_FMM
        case 1:
          if (comm_rank == 0)
            printf ("\n");
          if (comm_rank == 0)
            printf ("Running fmm ... \n");
          result = fcs_run (handle_fmm, local_particles, local_particles, particles,
                            charges, field, potentials);
          break;
#endif
#if FCS_ENABLE_MEMD       
        case -1:
          if (comm_rank == 0)
            printf ("\n");
          if (comm_rank == 0)
            printf ("Running memd ... \n");
          result = fcs_run (handle_memd, local_particles, local_particles, particles,
                            charges, field, potentials);
          break;
#endif
#if FCS_ENABLE_P2NFFT
        case 2:
          if (comm_rank == 0)
            printf ("\n");
          if (comm_rank == 0)
            printf ("Running p2nfft ... \n");
          result = fcs_run (handle_p2nfft, local_particles, local_particles, particles,
                            charges, field, potentials);
          break;
#endif
#if FCS_ENABLE_P3M        
        case 3:
          if (comm_rank == 0)
            printf ("\n");
          if (comm_rank == 0)
            printf ("Running p3m ... \n");
          result = fcs_run (handle_p3m, local_particles, local_particles, particles,
                            charges, field, potentials);
          break;
#endif
#if FCS_ENABLE_PEPC
        case 4:
          if (comm_rank == 0)
            printf ("\n");
          if (comm_rank == 0)
            printf ("Running pepc ... \n");
          result = fcs_run (handle_pepc, local_particles, local_particles, particles,
                            charges, field, potentials);
          break;
#endif
#if FCS_ENABLE_VMG        
        case 5:
          if (comm_rank == 0)
            printf ("\n");
          if (comm_rank == 0)
            printf ("Running vmg ... \n");
          result = fcs_run (handle_vmg, local_particles, local_particles, particles,
                            charges, field, potentials);
          break;
#endif
        default:
          printf( "Skipping step due to missing method \n");
          continue;
          break;
      }
      if (comm_rank == 0)
        fcs_result_print_result (result);
      fcs_result_destroy (result);
//       if(comm_rank == 0) printf("Moving particles: %d... \n", i+1);
//       for (j = 0; j < local_particles; ++j)
//       {
//         for (k = 0; k < 3; ++k)
//         {
//           *(velocities + j * 3 + k) +=
//             *(field + j * 3 + k) * *(charges + j) / *(masses + j) * scaling_factor;
//           *(particles + j * 3 + k) +=
//             *(velocities + j * 3 + k) * time_step;
//           if (*(periodicity + k) == 1)
//           {
//             if (*(particles + j * 3 + k) > INTERFACE_TEST_BOX_SIZE)
//               *(particles + j * 3 + k) -= INTERFACE_TEST_BOX_SIZE;
//             else if (*(particles + j * 3 + k) < 0.0l)
//               *(particles + j * 3 + k) += INTERFACE_TEST_BOX_SIZE;
//           }
//         }
//       }

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
              fprintf (xyz_file, "%" FCS_LMOD_INT "d\n%" FCS_LMOD_INT "d %s %s %d %" FCS_LMOD_INT "d %c %c %c\n",
                number_of_particles, i+1, "C", (i%4==0)?"pp3mg":(i%4==1)?"fmm":(i%4==2)?"p2nfft":(i%4==3)?"p3m":"unknown", comm_size,
                number_of_runs,
                (*periodicity == 1) ? 'T' : 'F',
                (*(periodicity + 1) == 1) ? 'T' : 'F',
                (*(periodicity + 2) == 1) ? 'T' : 'F');
            }
            for (k = 0; k < local_particles; ++k)
            {
              fprintf (xyz_file, "O %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e\n",
                *(particles + 3 * k), *(particles + 3 * k + 1),
                *(particles + 3 * k + 2),
                *(charges + k),
                *(velocities + 3 * k),
                *(velocities + 3 * k + 1),
                *(velocities + 3 * k + 2));
                /**(potentials + k),
                    *(field + 3 * k),
                    *(field + 3 * k + 1),
                    *(field + 3 * k + 2));*/
            }
            fclose (xyz_file);

          }
          MPI_Barrier (communicator);
        }

      }

    }

//   if (comm_rank == 0)
//     printf ("\n");
//   if (comm_rank == 0)
//     printf ("Parameters... \n");
//   if (comm_rank == 0)
//     fcs_print_parameters (handle);
  if (comm_rank == 0)
    fcs_result_print_result (result);
  fcs_result_destroy (result);

  free (box_a);
  free (box_b);
  free (box_c);
  free (offset);

#if FCS_ENABLE_PP3MG
  fcs_destroy (handle_pp3mg);
#endif
#if FCS_ENABLE_FMM  
  fcs_destroy (handle_fmm);
#endif
#if FCS_ENABLE_P2NFFT  
  fcs_destroy (handle_p2nfft);
#endif
#if FCS_ENABLE_P3M
  fcs_destroy (handle_p3m);
#endif
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
