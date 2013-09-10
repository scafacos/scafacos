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

#ifndef __DIRECTC_H__
#define __DIRECTC_H__


#ifdef __cplusplus
extern "C" {
#endif


#include "common/near/near.h"


typedef struct _fcs_directc_t
{
  fcs_float box_base[3], box_a[3], box_b[3], box_c[3];
  fcs_int periodicity[3];

  fcs_int nparticles, max_nparticles;
  fcs_float *positions, *charges, *field, *potentials;

  fcs_int in_nparticles;
  fcs_float *in_positions, *in_charges;

  fcs_int out_nparticles;
  fcs_float *out_positions, *out_field, *out_potentials;

  fcs_float virial[9];

  fcs_int periodic_images[3];
  fcs_float cutoff;
  fcs_int cutoff_with_near;

  fcs_float max_particle_move;

  fcs_int resort;
  fcs_near_resort_t near_resort;

} fcs_directc_t;


void fcs_directc_create(fcs_directc_t *directc);
void fcs_directc_destroy(fcs_directc_t *directc);
void fcs_directc_set_system(fcs_directc_t *directc, fcs_float *box_base, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_int *periodicity);
void fcs_directc_set_particles(fcs_directc_t *directc, fcs_int nparticles, fcs_int max_nparticles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials);
void fcs_directc_set_in_particles(fcs_directc_t *directc, fcs_int in_nparticles, fcs_float *in_positions, fcs_float *in_charges);
void fcs_directc_set_out_particles(fcs_directc_t *directc, fcs_int out_nparticles, fcs_float *out_positions, fcs_float *out_field, fcs_float *out_potentials);
void fcs_directc_set_periodic_images(fcs_directc_t *directc, fcs_int *periodic_images);
void fcs_directc_set_cutoff(fcs_directc_t *directc, fcs_float cutoff);
void fcs_directc_get_cutoff(fcs_directc_t *directc, fcs_float *cutoff);
void fcs_directc_set_cutoff_with_near(fcs_directc_t *directc, fcs_int cutoff_with_near);
void fcs_directc_set_max_particle_move(fcs_directc_t *directc, fcs_float max_particle_move);
void fcs_directc_set_resort(fcs_directc_t *directc, fcs_int resort);
void fcs_directc_get_resort(fcs_directc_t *directc, fcs_int *resort);
void fcs_directc_get_resort_availability(fcs_directc_t *directc, fcs_int *availability);
void fcs_directc_get_resort_particles(fcs_directc_t *directc, fcs_int *resort_particles);
void fcs_directc_run(fcs_directc_t *directc, MPI_Comm comm);
void fcs_directc_resort_ints(fcs_directc_t *directc, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm);
void fcs_directc_resort_floats(fcs_directc_t *directc, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm);
void fcs_directc_resort_bytes(fcs_directc_t *directc, void *src, void *dst, fcs_int n, MPI_Comm comm);


#ifdef __cplusplus
}
#endif


#endif /* __DIRECTC_H__ */
