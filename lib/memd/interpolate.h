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

#ifndef _MEMD_INTERPOLATE_H
#define _MEMD_INTERPOLATE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "data_types.h"

/** finds current lattice site of each particle.
 calculates charge interpolation on cube. */
void ifcs_memd_distribute_particle_charges(memd_struct* memd);

/** Does the actual calculation of the gradient. Parameters:
 @param rel  3dim-array of relative position in cube,
 @param q    charge,
 @param grad huge gradient array to write into
 */
void ifcs_memd_calc_charge_gradients(fcs_float *rel, fcs_float q, fcs_float *grad);

/** finds correct lattice cube for particles.
 calculates charge gradients on this cube.
 @param grad gradient array to write into
 */
void ifcs_memd_update_charge_gradients(memd_struct* memd, fcs_float *grad);

#endif