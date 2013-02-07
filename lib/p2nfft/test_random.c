/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (C) 2011 Sebastian Banert
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

#include <p2nfft.h>
#include <fcs.h>

int main (
    int argc, char** argv
    )
{
  fcs_int num_particles = 8;
  fcs_int i;
  fcs_float positions[24] = {
    0.25, 0.25, 0.25,
    0.25, 0.25, 0.75,
    0.25, 0.75, 0.25,
    0.75, 0.25, 0.25,
    0.25, 0.75, 0.75,
    0.75, 0.25, 0.75,
    0.75, 0.75, 0.25,
    0.75, 0.75, 0.75};
  fcs_float charges[8] = {1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0};
  fcs_float field[24];
  fcs_float potentials[24];
  MPI_Comm comm = MPI_COMM_WORLD;
  
  /* Calculate this system via FCS direct solver */
  fcs_float box_a[] = { 1.0, 0.0, 0.0 };
  fcs_float box_b[] = { 0.0, 1.0, 0.0 };
  fcs_float box_c[] = { 0.0, 0.0, 1.0 };
  fcs_int periodicity[] = {1, 1, 1};
  
  FCS fcs_handle;
  FCSResult fcs_result;
  
  fcs_result = fcs_init(&fcs_handle, "DIRECT", comm);
  fcs_result = fcs_common_set(fcs_handle, 1, box_a, box_b, box_c, periodicity, num_particles, num_particles);
  fcs_result = fcs_DIRECT_require_virial(fcs_result, 0);
  fcs_result = fcs_tune(fcs_handle, 8, 8, positions, charges);
  fcs_result = fcs_run(fcs_handle, 8, 8, positions, charges, field, potentials);

  printf("Potentials via FCS DIRECT:\n");
  for (i = 0; i < num_partcles; ++i)
    printf("%f\n", potentials[i]);
  printf("\n");
  free_destroy(fcs_handle);

  /* Calculate this system via P2NFFT */
  void* rd = NULL;
  ifcs_p2nfft_init(&rd, comm);
  ifcs_p2nfft_tune(rd, num_particles, positions, charges, 1.0);
  ifcs_p2nfft_run(rd, num_particles, positions, charges, potentials, field);

  printf("Potentials via P2NFFT\n");
  for(i = 0; i<num_particles; ++i)
    printf("%f\n", potentials[i]);
  ifcs_p2nfft_destroy(rd);

  return 0;
}
