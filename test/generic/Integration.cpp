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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include <mpi.h>

#include "fcs.h"
#include "common/fcs-common/FCSCommon.h"

#include "common.hpp"
#include "Integration.hpp"


using namespace std;


static void parse_conf(integration_t *integ, char *conf)
{
  char *c, *param, *cur;

  c = new char[strlen(conf) + 1];

  strcpy(c, conf);

  cur = c;

#define PARSE_NEXT(_c_) do { (_c_) = strchr((_c_), ','); if (_c_) { *(_c_) = 0; ++(_c_); } } while (0)

  while (cur)
  {
    param = cur;

    PARSE_NEXT(cur);

    if (strcmp(param, "delta_t") == 0) { integ->delta_t = atof(cur); PARSE_NEXT(cur); }
    else if (strcmp(param, "max_move") == 0) { integ->max_move = atoi(cur); PARSE_NEXT(cur); }
    else if (strcmp(param, "output_steps") == 0) { integ->output_steps = atoi(cur); PARSE_NEXT(cur); }
  }

#undef PARSE_NEXT

  delete[] c;
}


void integ_setup(integration_t *integ, fcs_int time_steps, fcs_int resort, char *conf)
{
  integ->time_steps = time_steps;
  
  integ->delta_t = 1e-2;

  integ->resort = resort;

  integ->max_move = 1;
  
  integ->output_steps = 0;

  parse_conf(integ, conf);
}


void integ_system_setup(integration_t *integ, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_float *offset, fcs_int *periodicity)
{
  integ->box_a[0] = box_a[0]; integ->box_a[1] = box_a[1]; integ->box_a[2] = box_a[2];
  integ->box_b[0] = box_b[0]; integ->box_b[1] = box_b[1]; integ->box_b[2] = box_b[2];
  integ->box_c[0] = box_c[0]; integ->box_c[1] = box_c[1]; integ->box_c[2] = box_c[2];
  integ->offset[0] = offset[0]; integ->offset[1] = offset[1]; integ->offset[2] = offset[2];
  integ->periodicity[0] = periodicity[0]; integ->periodicity[1] = periodicity[1]; integ->periodicity[2] = periodicity[2];
}


void integ_print_settings(integration_t *integ, const char *prefix)
{
  cout << prefix << "Time steps:   " << integ->time_steps << endl;
  cout << prefix << "delta_t:      " << integ->delta_t << endl;
  cout << prefix << "resort:       " << integ->resort << endl;
  cout << prefix << "max_move:     " << integ->max_move << endl;
  cout << prefix << "output_steps: " << integ->output_steps << endl;
}


void integ_init(integration_t *integ, fcs_int nparticles, fcs_float *velocities, fcs_float *field)
{
  for (fcs_int i = 0; i < nparticles; ++i)
  {
    velocities[3 * i + 0] = velocities[3 * i + 1] = velocities[3 * i + 2] = 0.0;
    field[3 * i + 0] = field[3 * i + 1] = field[3 * i + 2] = 0.0;
  }
}


void integ_update_velocities(integration_t *integ, fcs_int nparticles, fcs_float *v, fcs_float *f_old, fcs_float *f_cur, fcs_float *q)
{
  for (fcs_int i = 0; i < nparticles; ++i)
  {
    v[3 * i + 0] += 0.5 * (f_old[3 * i + 0] + f_cur[3 * i + 0]) * q[i] * integ->delta_t;
    v[3 * i + 1] += 0.5 * (f_old[3 * i + 1] + f_cur[3 * i + 1]) * q[i] * integ->delta_t;
    v[3 * i + 2] += 0.5 * (f_old[3 * i + 2] + f_cur[3 * i + 2]) * q[i] * integ->delta_t;
  }
}


void integ_update_positions(integration_t *integ, fcs_int nparticles, fcs_float *pos, fcs_float *pos_old, fcs_float *v_old, fcs_float *f_old, fcs_float *q, fcs_float *max_particle_move)
{
  fcs_float old_pos[3], d;


  *max_particle_move = 0;

  if (pos_old == NULL) pos_old = pos;

  for (fcs_int i = 0; i < nparticles; ++i)
  {
    old_pos[0] = pos_old[3 * i + 0];
    old_pos[1] = pos_old[3 * i + 1];
    old_pos[2] = pos_old[3 * i + 2];

    pos[3 * i + 0] = old_pos[0] + integ->delta_t * (v_old[3 * i + 0] + 0.5 * f_old[3 * i + 0] * q[i] * integ->delta_t);
    pos[3 * i + 1] = old_pos[1] + integ->delta_t * (v_old[3 * i + 1] + 0.5 * f_old[3 * i + 1] * q[i] * integ->delta_t);
    pos[3 * i + 2] = old_pos[2] + integ->delta_t * (v_old[3 * i + 2] + 0.5 * f_old[3 * i + 2] * q[i] * integ->delta_t);

    d = fcs_sqrt((old_pos[0] - pos[3 * i + 0]) * (old_pos[0] - pos[3 * i + 0]) + (old_pos[1] - pos[3 * i + 1]) * (old_pos[1] - pos[3 * i + 1]) + (old_pos[2] - pos[3 * i + 2]) * (old_pos[2] - pos[3 * i + 2]));

    if (d > *max_particle_move) *max_particle_move = d;
  }
}


void integ_correct_positions(integration_t *integ, fcs_int nparticles, fcs_float *pos)
{
  /* wrap particle positions of periodic dimensions */
  fcs_wrap_positions(nparticles, pos, integ->box_a, integ->box_b, integ->box_c, integ->offset, integ->periodicity);
  
  /* increase particle system in open dimensions to enclose all particles */
  fcs_expand_system_box(nparticles, pos, integ->box_a, integ->box_b, integ->box_c, integ->offset, integ->periodicity);
}
