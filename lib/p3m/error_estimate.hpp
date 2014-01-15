/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
#ifndef _P3M_ERROR_ESTIMATE_H
#define _P3M_ERROR_ESTIMATE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <mpi.h>
#include "types.hpp"

namespace ScaFaCoS {
  namespace P3M {
    /** Determines a value for alpha that achieves the wanted_error, if at
        all possible. Also returns the achieved errors with these
        parameters. Check whether wanted_error > achieved_error to see
        whether the required error can actually be met.
    */
    void determine_good_alpha(data_struct *d);

    /** Calculates the rms error estimate in the force (as described in
        the book of Hockney and Eastwood (Eqn. 8.23) for a system of N
        randomly distributed particles.
    */
    void compute_error_estimate(data_struct *d);

    /** Calculates the reciprocal space contribution to the rms error in the
        force (as described in the book of Hockney and Eastwood
        (Eqn. 8.23) (for a system of N randomly distributed particles in a
        cubic box).
    */
    void k_space_error(data_struct *d);
  }
}
#endif
