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

/** Computes the charge assignment function of for the \a i'th degree
    at value \a x. */
fcs_float ifcs_p3m_caf(fcs_int i, fcs_float x, fcs_int cao_value);
fcs_float ifcs_p3m_caf_d(fcs_int i, fcs_float x, fcs_int cao_value);
void ifcs_p3m_interpolate_charge_assignment_function(ifcs_p3m_data_struct *d);

#endif
