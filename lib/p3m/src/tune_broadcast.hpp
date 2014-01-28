/*
  Copyright (C) 2013,2014 Olaf Lenz
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
#ifndef _P3M_TUNE_BROADCAST_HPP
#define _P3M_TUNE_BROADCAST_HPP
#include <config.h>

#include "types.hpp"

namespace ScaFaCoS {
  namespace P3M {
    /* Events during tuning */
    const int CMD_FAILED = -1;
    const int CMD_FINISHED = 0;
    const int CMD_COMPUTE_ERROR_ESTIMATE = 1;
    const int CMD_TIMING = 2;
    
    void tune_broadcast_command(data_struct *d, fcs_int command);
    
    void 
    tune_broadcast_slave(data_struct *d, fcs_int num_particles,
                         fcs_float *positions, fcs_float *charges);
  }
}
#endif
