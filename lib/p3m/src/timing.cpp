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
#include "utils.hpp"
#include "p3m.hpp"
#include "tune_broadcast.hpp"
#include "timing.hpp"
#include "prepare.hpp"
#include <cstdlib>

namespace ScaFaCoS {
  namespace P3M {
    void timing(data_struct *d, fcs_int num_particles, 
                fcs_float *positions, fcs_float *charges) {
      if (d->comm.rank == 0)
        tune_broadcast_command(d, CMD_TIMING);

      fcs_float *fields = new fcs_float[3*num_particles];
      fcs_float *potentials = new fcs_float[num_particles];

      prepare(d);

      /* store require_timings */
      timingEnum require_timings_before = d->require_timings;
      if(d->require_timings == NONE || d->require_timings == FULL)
        d->require_timings = ESTIMATE_ALL;
      run(d, num_particles, positions, charges, fields, potentials);
      /* Afterwards, d->timings is set */
      /* restore require_timings */
      d->require_timings = require_timings_before;

      delete[] fields;
      delete[] potentials;
    }
  }
}

