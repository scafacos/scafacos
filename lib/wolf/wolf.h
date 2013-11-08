/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann
  
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

#ifndef __WOLF_H__
#define __WOLF_H__


#ifdef __cplusplus
extern "C" {
#endif


#include "common/near/near.h"


typedef struct _ifcs_wolf_t
{
  fcs_float box_base[3], box_a[3], box_b[3], box_c[3];
  fcs_int periodicity[3];

  fcs_int nparticles, max_nparticles;
  fcs_float *positions, *charges, *field, *potentials;

  fcs_float virial[9];

  fcs_float cutoff, alpha;

  fcs_float max_particle_move;

  fcs_int resort;
  fcs_near_resort_t near_resort;

} ifcs_wolf_t;


void ifcs_wolf_create(ifcs_wolf_t *wolf);
void ifcs_wolf_destroy(ifcs_wolf_t *wolf);
void ifcs_wolf_set_system(ifcs_wolf_t *wolf, fcs_float *box_base, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_int *periodicity);
void ifcs_wolf_set_particles(ifcs_wolf_t *wolf, fcs_int nparticles, fcs_int max_nparticles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials);
void ifcs_wolf_set_cutoff(ifcs_wolf_t *wolf, fcs_float cutoff);
void ifcs_wolf_get_cutoff(ifcs_wolf_t *wolf, fcs_float *cutoff);
void ifcs_wolf_set_alpha(ifcs_wolf_t *wolf, fcs_float alpha);
void ifcs_wolf_get_alpha(ifcs_wolf_t *wolf, fcs_float *alpha);
void ifcs_wolf_set_max_particle_move(ifcs_wolf_t *wolf, fcs_float max_particle_move);
void ifcs_wolf_set_resort(ifcs_wolf_t *wolf, fcs_int resort);
void ifcs_wolf_get_resort(ifcs_wolf_t *wolf, fcs_int *resort);
void ifcs_wolf_get_resort_availability(ifcs_wolf_t *wolf, fcs_int *availability);
void ifcs_wolf_get_resort_particles(ifcs_wolf_t *wolf, fcs_int *resort_particles);
void ifcs_wolf_run(ifcs_wolf_t *wolf, MPI_Comm comm);
void ifcs_wolf_resort_ints(ifcs_wolf_t *wolf, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm);
void ifcs_wolf_resort_floats(ifcs_wolf_t *wolf, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm);
void ifcs_wolf_resort_bytes(ifcs_wolf_t *wolf, void *src, void *dst, fcs_int n, MPI_Comm comm);


#ifdef __cplusplus
}
#endif


#endif /* __WOLF_H__ */
