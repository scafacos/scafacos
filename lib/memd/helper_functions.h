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

#ifndef _MEMD_HELPER_FUNCTIONS_H
#define _MEMD_HELPER_FUNCTIONS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSResult.h"
#include "data_types.h"

/** Turns the 3D index into a linear index.
 adim contains the dimensions of the lattice in
 the 3 directions: adim[0] = x_max, ..., adim[2] = z_max
 z is the first index!!!
 @return linear index
 @param x           index in x direction
 @param y           index in y direction
 @param z           index in z direction
 @param latticedim  dimensions of the lattice
 */
fcs_int ifcs_memd_get_linear_index(fcs_int x, fcs_int y, fcs_int z, fcs_int latticedim[SPACE_DIM]);

/** Counts the total number of charged particles on 
 all processors.
 @return total number of charged particles in the system
 */
fcs_int ifcs_memd_count_charged_particles(memd_struct* memd);

/** Index shift is calculated for moving in various
 directions with the linear index.
 @return Offset of current node from base
 @param index_shift     amount to move in direction
 @param index_base      index of base on current processor
 @param axes            in which direction to move
 @param adim            dimensions of the local lattice
 */
fcs_int ifcs_memd_get_offset(fcs_int index_shift, fcs_int index_base, fcs_int axes, fcs_int adim[3]);

/** For any direction j, write the two other directions
 into the pointers in circular permutation.
 @param j    given first direction
 @param dir1 write second direction into
 @param dir2 write third direction into
 */
void ifcs_memd_calc_directions(fcs_int j, fcs_int* dir1, fcs_int*dir2);

/** Calculates the finite differences rotation in real space in mue-nue plane:
 \f$\frac{\partial}{\partial t}{D} = \nabla \times B\f$ (and prefactors plus current)
 The given "fcs_float* field" should be a B-Field!!
 @return rotation result
 @param mue       direction 1 of plane to rotate in
 @param nue       direction 2 of plane to rotate in
 @param field     input B-field
 @param Neighbor  neighbor lattice site
 @param index     index of current lattice site
 */
fcs_float ifcs_memd_calc_curl(fcs_int mue, fcs_int nue, fcs_float* field, fcs_int* Neighbor, fcs_int index);

/** Calculates the finite differences rotation in dual space in mue-nue plane:
 \f$\frac{\partial}{\partial t}{B} = - \nabla \times D / (\epsilon)\f$
 The given "fcs_float* field" should be a D-Field!!
 @return rotation result
 @param mue       direction 1 of plane to rotate in
 @param nue       direction 2 of plane to rotate in
 @param field     input D-field
 @param Neighbor  neighbor lattice site
 @param index     index of current lattice site
 */
fcs_float ifcs_memd_calc_dual_curl(fcs_int mue, fcs_int nue, fcs_float* field, fcs_int* Neighbor, fcs_int index);

/** updates all D-fields on links of the plaquette
 and the surroundings. delta was calculated before
 in function "ifcs_memd_perform_rot_move_inplane".
 @param mue        direction 1 of update
 @param nue        direction 2 of update
 @param Neighbor   neighbor lattice site
 @param index      index of current lattice site
 @param delta      by which amount to update field
 */
void ifcs_memd_update_plaquette(fcs_int mue, fcs_int nue, fcs_int* Neighbor, fcs_int index, fcs_float delta);

/** Basic sanity checks to see if the code will run.
 @return zero if everything is fine. -1 otherwise.
 */
FCSResult ifcs_memd_sanity_checks(memd_struct* memd);

/** counts the total number of charges in the system
 @param memd            MEMD struct
 @param num_particles   Number of particles on node
 @return                total number of charges
 */
fcs_int ifcs_memd_count_total_charges(memd_struct *memd, fcs_int num_particles);

/** Calculate closest upper number that is 2^n */
fcs_int ifcs_memd_get_next_higher_power_of_two(fcs_float number);

#endif
