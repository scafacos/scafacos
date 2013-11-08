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

#include "interpolate.h"
#include <math.h>
#include "helper_functions.h"
#include "FCSCommon.h"

/*********************************************/
/****** interpolate charges on lattice: ******/
/*********************************************/

/** Interpolation function in one dimension.
 Currently only linear interpolation conserves charge.
 @return interpolation value
 @param x relative position of particle
 */
fcs_float ifcs_memd_interpol1D(fcs_float x)
{	
    return x;
    /* Cosine interpolation would be: */
    /* return sqr(sin(M_PI_2*x)); */
}


/** Does the actual interpolating calculation.
 @param first  3dim-array lattice position,
 @param rel    3dim-array relative position in cube,
 @param q      charge to interpolate
 */
void ifcs_memd_interpolate_charge(memd_struct* memd, fcs_int *first, fcs_float *rel, fcs_float q)
{
    fcs_int i, k, l, m, index, temp_ind;
    fcs_int help_index[3];
    fcs_float temp;
    fcs_float help[SPACE_DIM];
	
    FOR3D(i) help[i] = 1. - rel[i];     /** relative pos. w.r.t. first */
    //  printf("first: %d %d %d\n", first[0], first[1], first[2]);
    /** calculate charges at each vertex */
    index = ifcs_memd_get_linear_index(first[0],first[1],first[2],memd->lparams.dim);
	
    FOR3D(i) {
        temp_ind = memd->neighbor[index][i];
        if(temp_ind == NOWHERE) help_index[i] = memd->lparams.volume; /* force huge index */
        else { /* incr. for x-neighbor */
            help_index[i] = ifcs_memd_get_offset(memd->lattice[memd->neighbor[index][i]].r[i], first[i], i, memd->lparams.dim);
        }
		
    }
	
    for(k=0;k<2;k++){   /* jumps from x- to x+ */
        for(l=0;l<2;l++){  /* jumps from y- to y+ */
            for(m=0;m<2;m++){ /* jumps from z- to z+ */      
                if(index < memd->lparams.volume) {
                    temp = q;
                    FOR3D(i) temp *= ifcs_memd_interpol1D(help[i]);
                    memd->lattice[index].charge += temp;
                }
				
                index+=help_index[2];
                help[2]=1.-help[2];
                help_index[2]=-help_index[2];
            }
            index+=help_index[1];
            help[1]=1.-help[1];
            help_index[1]=-help_index[1];
        }
        index+=help_index[0];
        help[0]=1.-help[0];
        help_index[0]=-help_index[0];
    }		
}

/** add charges from ghost cells to lattice sites. */
void ifcs_memd_accumulate_charge_from_ghosts(memd_struct* memd)
{
    memd_cell *cell;
    memd_particle* p;
    fcs_int i, c, d;
    fcs_int np;
    fcs_int flag_inner=0;
    fcs_int first[SPACE_DIM];
    fcs_float q;
    fcs_float pos[SPACE_DIM], rel[SPACE_DIM];
	
    /** loop over ghost cells */
    for (c = 0; c < memd->ghost_cells.n; c++) {
        cell = memd->ghost_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
            if( (q=p[i].q) != 0.0 ) {
                flag_inner=1;
                FOR3D(d) {
                    if(p[i].r[d]<memd->lparams.left_down_position[d]||p[i].r[d]>=memd->lparams.upper_right_position[d])
                    {flag_inner=0; break;}
                }
            }
            if(flag_inner) {
                FOR3D(d) {
                    pos[d]        = (p[i].r[d] - memd->lparams.left_down_position[d])* memd->parameters.inva;
                    first[d]      = (fcs_int) pos[d];
                    rel[d]        = pos[d] - first[d];
                }
                //      	fprintf(stderr,"pos: %f %f %f\n", p[i].r[0], p[i].r[1], p[i].r[2]);
                ifcs_memd_interpolate_charge(memd, first, rel, q);
            }
        }      
    }  
	
}


/** finds current lattice site of each particle.
 calculates charge interpolation on cube. */
void ifcs_memd_distribute_particle_charges(memd_struct* memd)
{	
    memd_cell *cell;
    memd_particle* p;
    fcs_int i, c, d;
    fcs_int np;
    fcs_int first[SPACE_DIM];
    fcs_float q;
    fcs_float pos[SPACE_DIM], rel[SPACE_DIM];
    
    for(i=0;i<memd->lparams.volume;i++) memd->lattice[i].charge = 0.;
    /** === charge assignment === */ 
    /** loop over inner cells */
    for (c = 0; c < memd->local_cells.n; c++) {
        cell = memd->local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
            if( (q=p[i].q) != 0.0 ) {
                FOR3D(d) {
                    pos[d]        = (p[i].r[d] - memd->lparams.left_down_position[d])* memd->parameters.inva;
                    first[d]      = (fcs_int) pos[d];
                    rel[d]        = pos[d] - first[d];
                }
                ifcs_memd_interpolate_charge(memd, first, rel, q);
            }
        }      
    }
    ifcs_memd_accumulate_charge_from_ghosts(memd);	
}

/** Does the actual calculation of the gradient. Parameters:
 @param rel  3dim-array of relative position in cube,
 @param q    charge,
 @param grad huge gradient array to write into
 */
void ifcs_memd_calc_charge_gradients(fcs_float *rel, fcs_float q, fcs_float *grad)
{
    fcs_int i,l,m,index, d;
    fcs_float help[3];
    fcs_int dir1, dir2;
	
    FOR3D(i) help[i] = 1. - rel[i];     /* relative pos. w.r.t. x_int */
	
    index = 0;
	
    FOR3D(d) {
        ifcs_memd_calc_directions(d, &dir1, &dir2);
        for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
            for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */          
				
                /* with q!!! */
                grad[index] = - q * help[dir1] * help[dir2];
				
                index++;
                help[dir1] = 1.-help[dir1];
            }
            help[dir2] = 1.-help[dir2];
        }
    }
}



/** finds correct lattice cube for particles.
 calculates charge gradients on this cube.
 @param grad gradient array to write into
 */
void ifcs_memd_update_charge_gradients(memd_struct* memd, fcs_float *grad)
{
    memd_cell *cell;
    memd_particle* p;
    fcs_int i, c, d, ip;
    fcs_int np;
    fcs_int first[SPACE_DIM];
    fcs_float q;
    fcs_float pos[SPACE_DIM], rel[SPACE_DIM];
    
    /* === grad assignment for real particles and self-force === */ 
    ip = 0;
    for (c = 0; c < memd->local_cells.n; c++) {
        cell = memd->local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
            if( (q=p[i].q) != 0.0 ) {
                FOR3D(d) {
                    pos[d]        = (p[i].r[d] - memd->lparams.left_down_position[d])* memd->parameters.inva;
                    first[d]      = (fcs_int) pos[d];
                    rel[d]        = pos[d] - first[d];
                }
                ifcs_memd_calc_charge_gradients(rel, q, &grad[ip]);
                ip += 12;
            }
        }      
    }  
    
} 






