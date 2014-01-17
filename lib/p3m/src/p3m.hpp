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
#ifndef _P3M_P3M_HPP
#define _P3M_P3M_HPP

#include "types.hpp"

namespace ScaFaCoS {
  namespace P3M {
    void init(data_struct *d, MPI_Comm communicator);

    void destroy(data_struct *d);

    void tune(data_struct *d,
              fcs_int num_particles,
              fcs_int max_particles,
              fcs_float *positions, 
              fcs_float *charges);

    void run(data_struct* d,
             fcs_int _num_particles,
             fcs_int _max_num_particles,
             fcs_float *_positions, 
             fcs_float *_charges,
             fcs_float *_fields,
             fcs_float *_potentials);

  }
}

#endif
