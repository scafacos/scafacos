/*
  Copyright (C) 2011,2014 Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _P3M_SCAFACOS_H
#define _P3M_SCAFACOS_H

#include <config.h>
#include <mpi.h>
#include "FCSResult.h"

#ifdef __cplusplus
extern "C" {
#endif
    
  /** Initialize all structures, parameters and arrays needed for the 
   *  P3M algorithm and set their default values.
   */
  FCSResult ifcs_p3m_init(void **rd, MPI_Comm communicator);

  /** Clean up P3M memory allocations. */
  void ifcs_p3m_destroy(void *rd);

  fcs_int ifcs_p3m_check_triclinic_box(fcs_float ay, fcs_float az, fcs_float bz);
  
  void ifcs_p3m_set_near_field_flag(void *rd, fcs_int flag);
  FCSResult ifcs_p3m_set_triclinic_flag(void *rd);
  
  void ifcs_p3m_set_box_a(void *rd, fcs_float a);
  void ifcs_p3m_set_box_b(void *rd, fcs_float b);
  void ifcs_p3m_set_box_c(void *rd, fcs_float c);
  
  void ifcs_p3m_set_box_geometry(void *rd, const fcs_float *a,
          const fcs_float *b, const fcs_float *c);

  void ifcs_p3m_set_r_cut(void *rd, fcs_float alpha);
  void ifcs_p3m_set_r_cut_tune(void *rd);
  void ifcs_p3m_get_r_cut(void *rd, fcs_float *r_cut);

  void ifcs_p3m_set_alpha(void *rd, fcs_float alpha);
  void ifcs_p3m_set_alpha_tune(void *rd);
  void ifcs_p3m_get_alpha(void *rd, fcs_float *alpha);
  void ifcs_p3m_get_near_params(void* rd, fcs_float *alpha, fcs_float *offset);

  void ifcs_p3m_set_grid(void *rd, fcs_int mesh);
  void ifcs_p3m_set_grid_tune(void *rd);
  void ifcs_p3m_get_grid(void *rd, fcs_int *mesh);

  void ifcs_p3m_set_cao(void *rd, fcs_int cao);
  void ifcs_p3m_set_cao_tune(void *rd);
  void ifcs_p3m_get_cao(void *rd, fcs_int *cao);

  void ifcs_p3m_set_tolerance_field(void *rd, fcs_float tolerance_field);
  void ifcs_p3m_set_tolerance_field_tune(void *rd);
  void ifcs_p3m_get_tolerance_field(void *rd, fcs_float* tolerance_field);

  void ifcs_p3m_require_total_energy(void *rd, fcs_int flag);
  FCSResult ifcs_p3m_get_total_energy(void *rd, fcs_float *total_energy);
  
  void ifcs_p3m_set_potential_shift(void *rd, fcs_int flag);
  void ifcs_p3m_get_potential_shift(void *rd, fcs_int *flag);

  void ifcs_p3m_require_timings(void *rd, fcs_int flag);
  FCSResult 
  ifcs_p3m_get_timings(void *rd, 
                       double *timing, 
                       double *timing_near_field, 
                       double *timing_far_field);

  FCSResult 
  ifcs_p3m_tune(void* rd, fcs_int num_particles, fcs_int max_particles,
                fcs_float *positions, fcs_float *charges);
  
  FCSResult 
  ifcs_p3m_run(void* rd, fcs_int num_particles, fcs_int max_particles,
               fcs_float *positions, fcs_float *charges,
               fcs_float *fields, fcs_float *potentials);
  
#ifdef __cplusplus
}
#endif

#endif
