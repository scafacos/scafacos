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

#ifndef _MEMD_FORCES_H
#define _MEMD_FORCES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "data_types.h"

/** propagate the B-field via \f$\frac{\partial}{\partial t}{B} = \nabla\times D\f$ (and prefactor)
 CAREFUL: Usually this function is called twice, with dt/2 each time
 to ensure a time reversible integration scheme!
 @param dt time step for update. Should be half the MD time step
 */
void fcs_memd_propagate_B_field(memd_struct* memd, fcs_float dt);

/** Public function.
 Calculates the actual force on each particle
 by calling all other needed functions (except
 for fcs_memd_propagate_B_field) */
void fcs_memd_calc_forces(memd_struct* memd);

#endif