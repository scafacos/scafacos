/*
  Copyright (C) 2011,2012 Olaf Lenz, Michael Hofmann
  
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

#ifndef _INTEGRATE_HPP
#define _INTEGRATE_HPP


typedef struct _integration_t
{
  fcs_int time_steps;
  
  fcs_float delta_t;
  fcs_int resort;
  fcs_int max_move;
  fcs_int output_steps;

  fcs_float box_a[3], box_b[3], box_c[3], offset[3];
  fcs_int periodicity[3];

} integration_t;


void integ_setup(integration_t *integ, fcs_int time_steps, fcs_int resort, char *conf);
void integ_system_setup(integration_t *integ, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_float *offset, fcs_int *periodicity);
void integ_print_settings(integration_t *integ, const char *prefix = "");
void integ_init(integration_t *integ, fcs_int nparticles, fcs_float *velocities, fcs_float *field);
void integ_update_velocities(integration_t *integ, fcs_int nparticles, fcs_float *v, fcs_float *f_old, fcs_float *f_cur, fcs_float *q);
void integ_update_positions(integration_t *integ, fcs_int nparticles, fcs_float *pos, fcs_float *pos_old, fcs_float *v_old, fcs_float *f_old, fcs_float *q, fcs_float *max_particle_move);
void integ_correct_positions(integration_t *integ, fcs_int nparticles, fcs_float *pos);


#endif /* _INTEGRATE_HPP */
