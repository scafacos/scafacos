/*
 Copyright (C) 2010/2011/2012 Florian Fahrenberger
 
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

#ifndef _MEMD_INITIAL_SOLUTION_H
#define _MEMD_INITIAL_SOLUTION_H

#include "data_types.h"

/** calculates initial electric field configuration.
 currently uses simple and slow method of plaquettes and links.
 energy minimization takes up lots of time. */
void ifcs_memd_calc_init_e_field(memd_struct* memd);

#endif