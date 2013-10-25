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

#include <stdio.h>

#include "helper_functions.h"
#include "FCSCommon.h"
#include <mpi.h>

/*************************************/
/****** small helper functions: ******/
/*************************************/

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
fcs_int ifcs_memd_get_linear_index(fcs_int x, fcs_int y, fcs_int z, fcs_int latticedim[SPACE_DIM])
{
    return (z + latticedim[ZPLUS]*(y + latticedim[YPLUS]*x));
}

/** Counts the total number of charged particles on 
 all processors.
 @return total number of charged particles in the system
 */
fcs_int ifcs_memd_count_charged_particles(memd_struct* memd)
{  
    memd_cell *cell;
    memd_particle *part;
    fcs_int i,c,np;
    fcs_float node_sum, tot_sum;
	
    node_sum=0.0; 
    tot_sum =0.0;
	
    for (c = 0; c < memd->local_cells.n; c++) {
        cell = memd->local_cells.cell[c];
        part = cell->part;
        np   = cell->n;
        for(i=0;i<np;i++) 
            if( part[i].q != 0.0 ) node_sum += 1.0;
    }
	
    MPI_Reduce(&node_sum, &tot_sum, 1, FCS_MPI_FLOAT, MPI_SUM, 0, memd->mpiparams.communicator);
	
    return tot_sum;
}


/** Index shift is calculated for moving in various
 directions with the linear index.
 @return Offset of current node from base
 @param index_shift     amount to move in direction
 @param index_base      index of base on current processor
 @param axes            in which direction to move
 @param adim            dimensions of the local lattice
 */
fcs_int ifcs_memd_get_offset(fcs_int index_shift, fcs_int index_base, fcs_int axes, fcs_int adim[3])
{
    fcs_int dif;
    dif = index_shift - index_base;
    if(axes <= 1) dif *= adim[2];
    if(axes == 0) dif *= adim[1]; 
	
    return (dif);
}


/** For any direction j, write the two other directions
 into the pointers in circular permutation.
 @param j    given first direction
 @param dir1 write second direction into
 @param dir2 write third direction into
 */
void ifcs_memd_calc_directions(fcs_int j, fcs_int* dir1, fcs_int*dir2)
{
    *dir1 = *dir2 = -1;
    switch(j) {
        case 0 :
            *dir1 = 2;
            *dir2 = 1;
            break;
        case 1 :
            *dir1 = 2;
            *dir2 = 0;
            break;  
        case 2 :
            *dir1 = 1;
            *dir2 = 0;
            break;
    }  
}

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
fcs_float ifcs_memd_calc_curl(fcs_int mue, fcs_int nue, fcs_float* field, fcs_int* Neighbor, fcs_int index)
{
    fcs_float result;
	
    result = field[index+mue] + field[3*Neighbor[OPP_DIR(mue)]+nue] -
    field[3*Neighbor[OPP_DIR(nue)]+mue] - field[index+nue];
	
    return result;
}

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
fcs_float ifcs_memd_calc_dual_curl(fcs_int mue, fcs_int nue, fcs_float* field, fcs_int* Neighbor, fcs_int index)
{
    fcs_float res;
	
    res = field[index+mue] + field[3*Neighbor[mue]+nue] -
    field[3*Neighbor[nue]+mue] - field[index+nue];
	
    return res;
}


/** updates all D-fields on links of the plaquette
 and the surroundings. delta was calculated before
 in function "ifcs_memd_perform_rot_move_inplane".
 @param mue        direction 1 of update
 @param nue        direction 2 of update
 @param Neighbor   neighbor lattice site
 @param index      index of current lattice site
 @param delta      by which amount to update field
 */
void ifcs_memd_update_plaquette(fcs_int mue, fcs_int nue, fcs_int* Neighbor, fcs_int index, fcs_float delta)
{
    fcs_int i = 3*index;
    memd.Dfield[i+mue]             += delta;
    memd.Dfield[3*Neighbor[mue]+nue] += delta;
    memd.Dfield[3*Neighbor[nue]+mue] -= delta;
    memd.Dfield[i+nue]             -= delta;  
}


/** Basic sanity checks to see if the code will run.
 @return zero if everything is fine. -1 otherwise.
 */
FCSResult ifcs_memd_sanity_checks(memd_struct* memd)
{
    FCSResult result = NULL;
    const char* fnc_name = "ifcs_memd_sanity_checks";
//    fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "You can only set 2 of the parameters bjerrum, temperature and epsilon.");
    
    /*
    char *errtxt;
    fcs_int ret = 0;
    fcs_int d;
    fcs_int max_node_grid = 1;
    
    FOR3D(d) if(node_grid[d] > max_node_grid) max_node_grid = node_grid[d];
	
    if (memd.parameters.bjerrum == 0.) {
        errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{301 MEMD: bjerrum length is zero.} ");
        ret = -1;
    }
    else if ( (box_l[0] != box_l[1]) || (box_l[1] != box_l[2]) ) {
        errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{302 MEMD needs cubic box.} ");
        ret = -1;
    }
    if (!PERIODIC(0) || !PERIODIC(1) || !PERIODIC(2)) {
        errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{303 MEMD requires periodicity 1 1 1} ");
        ret = 1;
    }
    else if ( memd.parameters.mesh%max_node_grid != 0 ) {
        errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{304 MEMD: meshsize is incompatible with number of processes.} ");
        ret = -1;
    }
    else if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
        errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{305 MEMD requires domain-decomposition cellsystem.} ");
        ret = -1;
    }
    else if (dd.use_vList) {
        errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{306 MEMD requires no Verlet Lists.} ");
        ret = -1;
    }
    else if (memd.parameters.f_mass < 2. * time_step * time_step * memd.parameters.a * memd.parameters.a) {
        errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{307 MEMD: Speed of light is set too high. Increase f_mass.} ");
        ret = -1;      
    }
    else if (memd.parameters.a < skin) {
        errtxt = runtime_error(128);
        ERROR_SPRINTF(errtxt, "{308 MEMD: Skin should be smaller than MEMD mesh size.} ");
        ret = -1;
    }
    */
    return result;
}


/** Calculate number of charged particles, the sum of the squared
 charges and the squared sum of the charges. */
fcs_int ifcs_memd_count_total_charges(memd_struct *memd, fcs_int num_particles) {
    int i;
    
    fcs_int node_sum, total_sum;

    node_sum = num_particles;
    total_sum = 0;

    MPI_Allreduce(&node_sum, &total_sum, 1, FCS_MPI_INT, MPI_SUM, memd->mpiparams.communicator);
    
    return total_sum;
}

/** Calculate closest upper number that is 2^n */
fcs_int ifcs_memd_get_next_higher_power_of_two(fcs_float number) {
    return pow(2, ceil(log(number)/log(2)));
}

