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

#ifndef _P3M_INFLUENCE_FUNCTION_H
#define _P3M_INFLUENCE_FUNCTION_H

#include "types.hpp"

namespace P3M {
  
  void 
  calc_influence_function_ik(data_struct *d);
  
  void 
  calc_influence_function_iki(data_struct *d);
  
  void 
  calc_influence_function_adi(data_struct *d);
  
}

#endif
