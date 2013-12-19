/*
  Copyright (C) 2011 Olaf Lenz
  
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
#ifndef _P3M_CAF_HPP
#define _P3M_CAF_HPP

#include "types.hpp"

void ifcs_p3m_interpolate_charge_assignment_function(ifcs_p3m_data_struct *d);

namespace P3M {
  class ChargeAssignmentFunction {
  public:
    /** Computes the \a i'th point of the CAF centered around \a r0 of
        the charge assignment function of order \a cao. */
    static fcs_float compute(fcs_int i, fcs_float r0, fcs_int cao);
    static fcs_float compute_derivative(fcs_int i, fcs_float r0, fcs_int cao);
    
    // /** Get the \a i'th point of the CAF centered around r0. */
    // fcs_float get(const fcs_float r0, const fcs_int i) = 0;
    // /** Get a pointer to an array that contains \a cao values of the caf. 
    //     Note that the array may be modified when this function is
    //     called a second time.
    //     The array is owned by the object, so do not try to free it.
    //  */
    // virtual fcs_float *get_array(const fcs_float r0) = 0;

    // virtual fcs_float get_diff(const fcs_float r0, const fcs_int i) = 0;
    // virtual fcs_float *get_diff_array(const fcs_float r0) = 0;
  };
}

#endif
