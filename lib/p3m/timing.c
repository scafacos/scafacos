/*
  Copyright (C) 2013 Olaf Lenz
  
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
#include "tune_broadcast.h"
#include "timing.h"
#include "prepare.h"
#include "run.h"
#include <stdlib.h>

/** Test run the method with the current parameters and return the elapsed CPU time. */
fcs_float
ifcs_p3m_timing
(ifcs_p3m_data_struct *d,
 fcs_int _num_particles, fcs_int _max_num_particles,
 fcs_float *_positions, fcs_float *_charges) {
  if (d->comm.rank == 0)
    ifcs_p3m_tune_broadcast_command(d, TIMING);

  ifcs_p3m_prepare(d, _max_num_particles);

  fcs_float *fields = malloc(_num_particles*3*sizeof(fcs_float));
  fcs_float *potentials = malloc(_num_particles*sizeof(fcs_float));
  fcs_float required_time =
    ifcs_p3m_runit(d, _num_particles, _max_num_particles,
                   _positions, _charges,
                   fields, potentials);

  free(fields);
  free(potentials);
  
  return required_time;
}

