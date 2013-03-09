/*
 * Copyright (C) 2011 Michael Pippig
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

#define TEST_BOX_SIZE 2 
#define TEST_N_PARTICLES (TEST_BOX_SIZE * TEST_BOX_SIZE * TEST_BOX_SIZE)

static void assert_fcs(FCSResult r)
{
  if (r) {
    fcsResult_printResult(r);
    MPI_Finalize();
    exit(-1);
  }
}

int main(int argc, char **argv)
{
  fcs_int total_num_particles = TEST_N_PARTICLES;
  fcs_int num_particles, max_num_particles;
  fcs_float box_size = TEST_BOX_SIZE;
  fcs_int pid, px, py, pz;
  fcs_float positions[3*TEST_N_PARTICLES];
  fcs_float charges[TEST_N_PARTICLES];
  fcs_float direct_fields[3*TEST_N_PARTICLES];
  fcs_float direct_potentials[TEST_N_PARTICLES];
  fcs_float p2nfft_fields[3*TEST_N_PARTICLES];
  fcs_float p2nfft_potentials[TEST_N_PARTICLES];
  fcs_float direct_virial[9], p2nfft_virial[9];
  
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int comm_rank, comm_size;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  pid = num_particles = 0;
  for (px = 0; px < TEST_BOX_SIZE; ++px) {
    for (py = 0; py < TEST_BOX_SIZE; ++py) {
      for (pz = 0; pz < TEST_BOX_SIZE; ++pz, ++pid) {
        if (pid % comm_size == comm_rank) {
          positions[3*num_particles] = px + 0.5;
          positions[3*num_particles + 1] = py + 0.5;
          positions[3*num_particles + 2] = pz + 0.5;
          charges[num_particles] = 1.0-((px + py + pz) % 2)*2;
          ++num_particles;
        }
      }
    }
  }
  max_num_particles = num_particles;

/* Debugging */
for(int t=0; t<6; t++)
  fprintf(stderr, "init positions[%d] = %" FCS_LMOD_FLOAT "f\n", t, positions[t]);

  fcs_float box_a[] = { box_size, 0.0, 0.0 };
  fcs_float box_b[] = { 0.0, box_size, 0.0 };
  fcs_float box_c[] = { 0.0, 0.0, box_size };
  fcs_float offset[] = {0.0, 0.0, 0.0};
  fcs_int periodicity[] = {0, 0, 0};
//  fcs_int periodicity[] = {1, 1, 1};

  FCS fcs_handle = NULL;
  FCSResult fcs_result = NULL;

  /* Calculate this system via FCS direct solver */
  fcs_result = fcs_init(&fcs_handle, "direct", comm);
  assert_fcs(fcs_result);

  fcs_result = fcs_set_common(fcs_handle, 1, box_a, box_b, box_c, offset, periodicity, total_num_particles);
  assert_fcs(fcs_result);

  fcs_result = fcs_tune(fcs_handle, num_particles, max_num_particles, positions, charges);
  assert_fcs(fcs_result);

  fcs_result = fcs_require_virial(fcs_handle, 1);
  assert_fcs(fcs_result);

/* Debugging */
for(int t=0; t<6; t++)
  fprintf(stderr, "before direct run: positions[%d] = %" FCS_LMOD_FLOAT "f\n", t, positions[t]);

  fcs_result = fcs_run(fcs_handle, num_particles, max_num_particles, positions, charges,
      direct_fields, direct_potentials);
  assert_fcs(fcs_result);

/* Debugging */
for(int t=0; t<6; t++)
  fprintf(stderr, "after direct run: positions[%d] = %" FCS_LMOD_FLOAT "f\n", t, positions[t]);

  fcs_result = fcs_get_virial(fcs_handle, direct_virial);
  assert_fcs(fcs_result);

  printf("Virial via FCS direct:\n");
  if(!comm_rank)
    printf("virial tensor = [%" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e; %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e; %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e]\n", direct_virial[0], direct_virial[1], direct_virial[2],
        direct_virial[3], direct_virial[4], direct_virial[5], direct_virial[6], direct_virial[7], direct_virial[8]);

  printf("Potentials via FCS direct:\n");
  for (pid = 0; pid < num_particles; ++pid)
    printf("%" FCS_LMOD_FLOAT "f\n", direct_potentials[pid]);

  printf("Fields via FCS direct:\n");
  for (pid = 0; pid < num_particles; ++pid)
    printf("[%" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f]\n", direct_fields[3*pid+0], direct_fields[3*pid+1], direct_fields[3*pid+2]);

  fcs_destroy(fcs_handle);


  /* set p2nfft specific parameters */
  fcs_float tolerance = 1.e-3;
  
  /* Calculate this system via FCS p2nfft solver */
  fcs_result = fcs_init(&fcs_handle, "p2nfft", comm);
  assert_fcs(fcs_result);

  fcs_result = fcs_set_common(fcs_handle, 1, box_a, box_b, box_c, offset, periodicity, total_num_particles);
  assert_fcs(fcs_result);

  fcs_set_tolerance(fcs_handle, FCS_TOLERANCE_TYPE_POTENTIAL, tolerance);
  fcs_result = fcs_tune(fcs_handle, num_particles, max_num_particles, positions, charges);
  assert_fcs(fcs_result);

  fcs_result = fcs_require_virial(fcs_handle, 1);
  assert_fcs(fcs_result);

/* Debugging */
for(int t=0; t<6; t++)
  fprintf(stderr, "test: positions[%d] = %" FCS_LMOD_FLOAT "f\n", t, positions[t]);

  fcs_result = fcs_run(fcs_handle, num_particles, max_num_particles, positions, charges,
      p2nfft_fields, p2nfft_potentials);
  assert_fcs(fcs_result);

  fcs_result = fcs_get_virial(fcs_handle, p2nfft_virial);
  assert_fcs(fcs_result);

  printf("Virial via FCS p2nfft:\n");
  if(!comm_rank)
    printf("virial tensor = [%" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e; %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e; %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e]\n", p2nfft_virial[0], p2nfft_virial[1], p2nfft_virial[2],
        p2nfft_virial[3], p2nfft_virial[4], p2nfft_virial[5], p2nfft_virial[6], p2nfft_virial[7], p2nfft_virial[8]);

  printf("Potentials via FCS p2nfft:\n");
  for (pid = 0; pid < num_particles; ++pid)
    printf("%" FCS_LMOD_FLOAT "f\n", p2nfft_potentials[pid]);

  printf("Fields via FCS p2nfft:\n");
  for (pid = 0; pid < num_particles; ++pid)
    printf("[%" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f]\n", p2nfft_fields[3*pid+0], p2nfft_fields[3*pid+1], p2nfft_fields[3*pid+2]);

  fcs_destroy(fcs_handle);

  /* Compare results of direct and p2nfft solver */

  fcs_float direct_energy = 0, p2nfft_energy = 0;
  for (pid = 0; pid < num_particles; pid++) {
    direct_energy += 0.5 * charges[pid] * direct_potentials[pid];
    p2nfft_energy += 0.5 * charges[pid] * p2nfft_potentials[pid];
  }

  fcs_float sqr_sum = 0.0;
  for (pid = 0; pid < num_particles; ++pid) {
    fcs_float d0 = p2nfft_fields[3*pid]   - direct_fields[3*pid];
    fcs_float d1 = p2nfft_fields[3*pid+1] - direct_fields[3*pid+1];
    fcs_float d2 = p2nfft_fields[3*pid+2] - direct_fields[3*pid+2];
    sqr_sum += d0*d0+d1*d1+d2*d2;
  }

  /* Reduce to global values */
  fcs_float direct_total_energy, p2nfft_total_energy, total_sqr_sum;
  MPI_Reduce(&direct_energy, &direct_total_energy, 1, FCS_MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&p2nfft_energy, &p2nfft_total_energy, 1, FCS_MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sqr_sum, &total_sqr_sum, 1, FCS_MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (!comm_rank) {
    printf("direct_energy = %" FCS_LMOD_FLOAT "f\n", direct_total_energy);
    printf("p2nfft_energy = %" FCS_LMOD_FLOAT "f\n", p2nfft_total_energy);
    printf("rms_error = %e\n", sqrt(total_sqr_sum / (fcs_float)total_num_particles));
  }

  MPI_Finalize();

  return 0;
}
