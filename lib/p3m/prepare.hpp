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
#ifndef _P3M_INIT_H
#define _P3M_INIT_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "types.hpp"

/* Call this function when all parameters are set correctly and the
   P3M data structures should be prepared. */
void 
ifcs_p3m_prepare(ifcs_p3m_data_struct *d, 
		 fcs_int max_particles);
#endif
