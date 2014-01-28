/*
 * Copyright (C) 2011 Olaf Lenz, Sebastian Banert, Michael Pippig
 *
 * This file is part of ScaFaCoS.
 * 
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *	
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fcs.h"

static void assert_fcs(FCSResult r)
{
  if (r) {
    fcs_result_print_result(r);
    MPI_Finalize();
    exit(-1);
  }
}

int main(int argc, char **argv)
{
  int comm_rank, comm_size;
  const char* method = "p2nfft";
/* DEBUG */  
//  const char* datafile = "../inp_data/p2nfft/debug_wall_small.dat";
  const char* datafile = "../inp_data/p3m/p3m_wall.dat";
  MPI_Comm comm = MPI_COMM_WORLD;
  fcs_int periodicity[3] = { 1, 1, 1 };
  fcs_int pid;
  fcs_float tolerance = 1.e-3;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if(!comm_rank){
    printf("-------------------\n");
    printf("Running p2nfft test\n");
    printf("-------------------\n");
  }

  fcs_float box_l[3] = { 10.0, 10.0, 10.0 };
  fcs_float offset[3] = {0.0, 0.0, 0.0};
/* DEBUG */  
//  int n_particles = 4;
  int n_particles = 300;

  fcs_float charges[300];
  fcs_float positions[900];
  fcs_float reference_forces[900];

  fcs_int local_particles = 0;
  fcs_float local_charges[300];
  fcs_float local_positions[900];
//   fcs_int global_particle_indices[300];


  if(!comm_rank)
    printf("Reading %s...\n", datafile);
  FILE *data = fopen(datafile, "r");
  if (!data) {
    fprintf(stderr, "ERROR: Can't read %s!", datafile);
    perror("ERROR");
    exit(1);
  }
  
  fcs_float charge_sum = 0.0;
  for (pid = 0; pid < n_particles; pid++) {
    fscanf(data, "%" FCS_CONV_FLOAT "f %" FCS_CONV_FLOAT "f %" FCS_CONV_FLOAT "f",
        &positions[3*pid], &positions[3*pid+1], &positions[3*pid+2]);
    fscanf(data, "%" FCS_CONV_FLOAT "f",
        &charges[pid]);
    fscanf(data, "%" FCS_CONV_FLOAT "f %" FCS_CONV_FLOAT "f %" FCS_CONV_FLOAT "f",
        &reference_forces[3*pid], &reference_forces[3*pid+1], &reference_forces[3*pid+2]);
    charge_sum += charges[pid];
  }

  fclose(data);


  FCS handle = NULL;
  FCSResult result = NULL;

  MPI_Barrier(MPI_COMM_WORLD);
  if(!comm_rank)
    printf("Initializing p2nfft...\n");
  result = fcs_init(&handle, method, comm);
  assert_fcs(result);

  if(!comm_rank)
    printf("Reading particles ... \n");
  for (pid = 0; pid < n_particles; ++pid) {
    if (pid % comm_size == comm_rank) {
      local_charges[local_particles] = charges[pid];
      local_positions[3*local_particles] = positions[3*pid];
      local_positions[3*local_particles+1] = positions[3*pid+1];
      local_positions[3*local_particles+2] = positions[3*pid+2];
//       global_particle_indices[local_particles] = pid;
      ++local_particles;
    }
  }

  if(!comm_rank)
    printf("Setting parameters...\n");

  fcs_float box_a[3] = { 0.0, 0.0, 0.0 };
  fcs_float box_b[3] = { 0.0, 0.0, 0.0 };
  fcs_float box_c[3] = { 0.0, 0.0, 0.0 };
  box_a[0] = box_l[0];
  box_b[1] = box_l[1];
  box_c[2] = box_l[2];

  result = fcs_set_common(handle, 1, box_a, box_b, box_c,
     offset, periodicity, n_particles);
  assert_fcs(result);

  /* Tuning */
  fcs_set_tolerance(handle, FCS_TOLERANCE_TYPE_FIELD, tolerance);

  fcs_float pot[3], dist[3], rcut, alpha;

  for(fcs_int t=0; t<3; t++){
    rcut = 3.0 + t; 
    fcs_p2nfft_set_r_cut(handle, rcut);

    if(!comm_rank)
      printf("Tuning p2nfft to accuracy %" FCS_LMOD_FLOAT "e...\n", tolerance);
    result = fcs_tune(handle, local_particles, n_particles, local_positions, local_charges);
    assert_fcs(result);
  
    result = fcs_p2nfft_get_alpha(handle, &alpha);
    assert_fcs(result);
  
    dist[0] = 0.0;
    result = fcs_compute_near_potential(handle, dist[0], &pot[0]);
    assert_fcs(result);

    dist[1] = 0.5*rcut;
    result = fcs_compute_near_potential(handle, dist[1], &pot[1]);
    assert_fcs(result);

    dist[2] = rcut;
    result = fcs_compute_near_potential(handle, dist[2], &pot[2]);
    assert_fcs(result);

    printf("rank = %d, alpha = %" FCS_LMOD_FLOAT "e, pot(%" FCS_LMOD_FLOAT "e) = %" FCS_LMOD_FLOAT "e, pot(%" FCS_LMOD_FLOAT "e) = %" FCS_LMOD_FLOAT "e, pot(%" FCS_LMOD_FLOAT "e) = %" FCS_LMOD_FLOAT "e\n",
        comm_rank, alpha, dist[0], pot[0], dist[1], pot[1], dist[2], pot[2]);
  }


  /* Far field computation */
//  if(!comm_rank)
//    printf("Running p2nfft (computing far fields and potentials)...\n");
//  result = fcs_run(handle, local_particles, n_particles,
//      local_positions, local_charges, far_fields, far_potentials);
//  assert_fcs(result);
//
//  /* Add components */
//  for (pid = 0; pid < local_particles; pid++) {
//    forces[3*pid]   = local_charges[pid] * far_fields[3*pid];
//    forces[3*pid+1] = local_charges[pid] * far_fields[3*pid+1];
//    forces[3*pid+2] = local_charges[pid] * far_fields[3*pid+2];
//
//    energies[pid] = 0.5 * local_charges[pid] * far_potentials[pid];
//  }
//  
//  /* Compare forces to reference field */
//  fcs_float sqr_sum = 0.0;
//  fcs_float sum_energy = 0.0;
//  for (pid = 0; pid < local_particles; pid++) {
//    sum_energy += energies[pid];
//
//    fcs_float d0 = forces[3*pid]   - reference_forces[3*global_particle_indices[pid]];
//    fcs_float d1 = forces[3*pid+1] - reference_forces[3*global_particle_indices[pid]+1];
//    fcs_float d2 = forces[3*pid+2] - reference_forces[3*global_particle_indices[pid]+2];
//    sqr_sum += d0*d0+d1*d1+d2*d2;
//  }
//
//  /* Reduce to global values */
//  fcs_float global_total_energy, global_sqr_sum;
//  MPI_Reduce(&sum_energy, &global_total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//  MPI_Reduce(&sqr_sum, &global_sqr_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//  if (!comm_rank) {
//    printf("sum_energy=%lf\n", global_total_energy);
//    printf("rms_error=%e\n", sqrt(global_sqr_sum / (fcs_float)n_particles));
//  }

  if(!comm_rank)
    printf("Finalizing...\n");
  fcs_destroy(handle);

  MPI_Finalize();
  if(!comm_rank)
    printf("Done.\n");

  return 0;
}
