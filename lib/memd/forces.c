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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "forces.h"
#include "helper_functions.h"
#include "communication.h"
#include "interpolate.h"
#include "FCSCommon.h"

/*********************************************/
/****** calculate currents and E-fields ******/
/*********************************************/

/** calculate the charge flux in direction "dir".
 @param q     charge of the moving particle
 @param help  the relative position in the cube to the opposite of left down front lattice site.
 @param flux  flux variable to write into
 @param dir   direction in which to calculate the flux
 */
void ifcs_memd_calc_charge_fluxes_1D(memd_struct* memd, fcs_float q, fcs_float *help, fcs_float *flux, fcs_int dir)
{
    /* at the moment works only for linear interpolation */
    fcs_int index, dir1, dir2;
    fcs_int l,m; 
    fcs_float q_scaled;
	
    q_scaled = q * memd->parameters.prefactor*memd->parameters.inva;
    index = 0;
	
    ifcs_memd_calc_directions(dir, &dir1, &dir2);   
	
    for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
        for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */   
			
            flux[index] = q_scaled * help[dir1]*help[dir2];
            index++;
			
            help[dir1] = 1. - help[dir1];
        }
        help[dir2] = 1. - help[dir2]; 
    }
}

/** Extend particle trajectories on the one time step and check if the
 trajectory intersects the cell boundary in direction dirx#
 @return zero if no particle crosses intersection
 @param delta    amount by which the particle moves
 @param r_new    current position of particle
 @param dir      direction in which to check for intersection
 @param first    position of particle within lattice cube
 @param t_step   time step
 @param identity particle ID
 */
fcs_int ifcs_memd_check_intersect_1D(fcs_float delta, fcs_float r_new, fcs_int dir, fcs_int first, fcs_float *t_step, fcs_int identity)
{	
    fcs_int candidateplane = -1; /* force alloc error */
    fcs_int f_crossing; 
    fcs_float r_old, temp;
    fcs_float ZERO = 0.0;
    fcs_float ONE  = 1.0;
	
    f_crossing = 0;
    r_old = r_new - delta;
    f_crossing = f_crossing||(r_old>=ONE||r_old<ZERO);
	
    if(dir==2) temp = 1.;
    else       temp = 0.5;
	
    if(f_crossing) {
		
        if(r_old >= ONE) candidateplane = ONE;
        if(r_old < ZERO) candidateplane = ZERO;
		
        /****** Update time step *********************/
        *t_step = temp * fabs((candidateplane-r_new)/delta);
    } /* end if crossing */
    else *t_step = temp;
    return f_crossing;
}

/** updates field on link coupling it with current.
 Force is multiplied by the time_step
 @param index   index of lattice site
 @param flux    charge flux in lattice site
 @param v       speed of particle
 @param dir     first direction of flux calculation
 */
void ifcs_memd_calc_e_field_on_link_1D(memd_struct* memd, fcs_int index, fcs_float *flux, fcs_float v, fcs_int dir)
{  
    fcs_int l, m, ind_flux, dir1, dir2;
    fcs_int temp_ind;
    fcs_int help_index[2];
    fcs_int* anchor_neighb;
    t_site* anchor_site;
	
    ifcs_memd_calc_directions(dir, &dir1, &dir2);
	
    anchor_neighb = &memd->neighbor[index][0]; 
    anchor_site = &memd->lattice[index];
	
    temp_ind = anchor_neighb[dir1];
    if(temp_ind == NOWHERE) help_index[0] = memd->lparams.volume;
    else
        help_index[0] = ifcs_memd_get_offset(memd->lattice[temp_ind].r[dir1], anchor_site->r[dir1], dir1, memd->lparams.dim);
    temp_ind = anchor_neighb[dir2];
    if(temp_ind == NOWHERE) help_index[1] = memd->lparams.volume;
    else
        help_index[1] = ifcs_memd_get_offset(memd->lattice[temp_ind].r[dir2], anchor_site->r[dir2], dir2, memd->lparams.dim);
	
	
    ind_flux = 0;
    for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
        for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */  
			
            if(index < memd->lparams.volume){
                memd->Dfield[3*index+dir] -= flux[ind_flux] * v;
            }
			
            ind_flux++; 
			
            index+=help_index[0];
            help_index[0]=-help_index[0];	
        }
        index+=help_index[1];
        help_index[1]=-help_index[1];     
    }
}


/** loop over all cells and call charge current functions
 @param p           Particle pointer
 @param ghost_cell  flag if cell is ghost cell
 */
void ifcs_memd_add_current_on_segment(memd_struct* memd, memd_particle *p, fcs_int ghost_cell)
{
    fcs_int d;
    fcs_int icoord, dir;
    fcs_int lat_index = -1; /* force alloc error */
    fcs_int first[SPACE_DIM];
    fcs_int f_crossing, flag_update_flux;
    t_dvector r_temp; 
    fcs_float inva_half;
    fcs_float delta;
    fcs_float flux[4], v;
    fcs_float v_inva[SPACE_DIM], v_invasq[SPACE_DIM];
    fcs_float pos[SPACE_DIM], help[3];
    fcs_float t_step;
	
    inva_half = 0.5*memd->parameters.inva;
	
    FOR3D(d) {
        pos[d]   = (p->r[d] - memd->lparams.left_down_position[d])* memd->parameters.inva;
        first[d] = (fcs_int) floor(pos[d]);
        r_temp[d]   = pos[d] - first[d]; /* it is the updated coord (we have to go back) */
        help[d]     = 1. - r_temp[d];
        v_inva[d]   = memd->parameters.inva * p->v[d];
        v_invasq[d] = memd->parameters.inva * v_inva[d];
    }
	
    flag_update_flux = 1;
    if(ghost_cell) {
        FOR3D(d) if(first[d]<memd->lparams.halo_left_down[d] || first[d]>= memd->lparams.halo_upper_right[d])
        {flag_update_flux = 0;break;}
    }
	
    if(flag_update_flux) {
        lat_index = ifcs_memd_get_linear_index(first[0], first[1], first[2], memd->lparams.dim);
    }
	
    /* loop coordinates in order x->y->z->y->x */
    for(dir=0; dir<5; dir++) {
        icoord = dir;
        if(dir>2) icoord = dir%2;
        if(icoord == 2) delta = v_inva[icoord];
        else            delta = 0.5 * v_inva[icoord];
		
        f_crossing = ifcs_memd_check_intersect_1D(delta, r_temp[icoord], icoord, first[icoord], &t_step, *p->identity);
		
        /* calculate flux */
        if(flag_update_flux) {
            ifcs_memd_calc_charge_fluxes_1D(memd, p->q, help, flux, icoord);
			
            v = t_step * v_invasq[icoord];
            ifcs_memd_calc_e_field_on_link_1D(memd, lat_index, flux, v, icoord);
        }
		
        if(f_crossing) {
            if(delta > 0.) {
                first[icoord]--;
                r_temp[icoord] += 1.;
            }
            else {
                first[icoord]++;
                r_temp[icoord] -= 1.;
            }
            if(icoord == 2) t_step = 1.  - t_step;
            else            t_step = 0.5 - t_step;
			
            if(ghost_cell){
                if(flag_update_flux) {
                    if(first[icoord]<memd->lparams.halo_left_down[icoord] || first[icoord]>= memd->lparams.halo_upper_right[icoord])
                    {flag_update_flux = 0;}
                }
                else {
                    flag_update_flux = 1;
                    FOR3D(d) if(first[d]<memd->lparams.halo_left_down[d] || first[d]>= memd->lparams.halo_upper_right[d])
                    {flag_update_flux = 0;break;}
                    if(flag_update_flux) ifcs_memd_calc_charge_fluxes_1D(memd, p->q, help, flux, icoord); 
                }
            }
			
            if(flag_update_flux) {
                v = t_step * v_invasq[icoord];
                lat_index = ifcs_memd_get_linear_index(first[0], first[1], first[2], memd->lparams.dim);
                ifcs_memd_calc_e_field_on_link_1D(memd, lat_index, flux, v, icoord);
            }
        }
        r_temp[icoord] -= delta;
        help[icoord]    = 1. - r_temp[icoord];
    }
}


/** Calculate fluxes and couple them with fields symplectically.  It
 is assumed that the particle can not cross more than one cell
 boundary per direction */
void ifcs_memd_couple_current_to_Dfield(memd_struct* memd)
{
    memd_cell *cell;
    memd_particle* p;
    fcs_int i, c, d, np;
    fcs_int flag_inner;
    fcs_float q;
    fcs_float r1, r2;
	
    /** loop over real particles */
    for (c = 0; c < memd->local_cells.n; c++) {
        cell = memd->local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
            if((q=p[i].q) != 0.) {
                /*	if(sim_time>49.08&&p[i].identity==231) */
                /*	  fprintf(stderr,"time=%f, v=(%f,%f,%f)\n",sim_time, p[i].v[0], p[i].v[1],p[i].v[2]); */
                ifcs_memd_add_current_on_segment(memd, &p[i], 0);
            }/* if particle.q != ZERO */
        }
    }
	
    /** loop over ghost particles */
    for (c = 0; c < memd->ghost_cells.n; c++) {
        cell = memd->ghost_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
            if((q=p[i].q) != 0.) {
                flag_inner = 1;
                FOR3D(d) {
                    r2 = p[i].r[d];
                    r1 = r2 - p[i].v[d];
                    if(((r2 < memd->lparams.left_down_position[d])&&(r1 < memd->lparams.left_down_position[d]))
                       ||((r2 >= memd->lparams.upper_right_position[d] && r1 >= memd->lparams.upper_right_position[d])))
                    {flag_inner = 0; break;}
                }
                if(flag_inner) {
                    ifcs_memd_add_current_on_segment(memd, &p[i], 1);
                }
            }/* if particle.q != ZERO */
        }
    }
}









/*******************************************/
/****** calculate B-fields and forces ******/
/*******************************************/

/** propagate the B-field via \f$\frac{\partial}{\partial t}{B} = \nabla\times D\f$ (and prefactor)
 CAREFUL: Usually this function is called twice, with dt/2 each time
 to ensure a time reversible integration scheme!
 @param dt time step for update. Should be half the MD time step
 */
void fcs_memd_propagate_B_field(memd_struct* memd, fcs_float dt)
{
    fcs_int x, y, z, i, offset, index;
    fcs_int xoffset, yoffset;
    fcs_float help = dt*memd->parameters.invsqrt_f_mass;
    /* B(t+h/2) = B(t-h/2) + h*curlE(t) */ 
	
    offset = ifcs_memd_get_linear_index(1,1,1, memd->lparams.dim);
    yoffset = memd->lparams.dim[2];
    xoffset = 2*memd->lparams.dim[2];
	
    for(x=0;x<memd->lparams.size[0];x++) {
        for(y=0;y<memd->lparams.size[1];y++) {
            for(z=0;z<memd->lparams.size[2];z++) {
				
                i = offset+z;
                index = 3*i;
                memd->Bfield[index+0] += - help*ifcs_memd_calc_dual_curl(1,2, memd->Dfield, memd->neighbor[i], index); 
                memd->Bfield[index+1] += - help*ifcs_memd_calc_dual_curl(2,0, memd->Dfield, memd->neighbor[i], index); 
                memd->Bfield[index+2] += - help*ifcs_memd_calc_dual_curl(0,1, memd->Dfield, memd->neighbor[i], index);  
            }
            offset += yoffset;
        }
        offset += xoffset;
    }
	
    fcs_memd_exchange_surface_patch(memd, memd->Bfield, 3, 0);
}

/** calculate D-field from B-field according to
 \f$\frac{\partial}{\partial t}{D} = \nabla\times B\f$ (and prefactors)
 @param dt MD time step
 */
void ifcs_memd_add_transverse_field(memd_struct* memd, fcs_float dt)
{
    fcs_int i, index;
    fcs_float invasq; 
    fcs_int x, y, z;
    fcs_int offset, xoffset, yoffset;
    fcs_float help;
	
    invasq = SQR(memd->parameters.inva);
    help = dt * invasq * memd->parameters.invsqrt_f_mass;
	
    /***calculate e-field***/ 
    offset = ifcs_memd_get_linear_index(1,1,1, memd->lparams.dim);
    yoffset = memd->lparams.dim[2];
    xoffset = 2*memd->lparams.dim[2];
    for(x=0;x<memd->lparams.size[0];x++) {
        for(y=0;y<memd->lparams.size[1];y++) {
            for(z=0;z<memd->lparams.size[2];z++) {
                /*  FORALL_INNER_SITES(x, y, z) { */
                /*    i = ifcs_memd_get_linear_index(x, y, z, memd->lparams.dim); */
                i = offset+z;
                index = 3*i;
                memd->Dfield[index  ] += help * ifcs_memd_calc_curl(2, 1, memd->Bfield, memd->neighbor[i], index);
                memd->Dfield[index+1] += help * ifcs_memd_calc_curl(0, 2, memd->Bfield, memd->neighbor[i], index);
                memd->Dfield[index+2] += help * ifcs_memd_calc_curl(1, 0, memd->Bfield, memd->neighbor[i], index);
            }
            offset += yoffset;
        }
        offset += xoffset;
    } 
	
    fcs_memd_exchange_surface_patch(memd, memd->Dfield, 3, 0);
}


/** interpolation function, solely for new self force correction. */
void ifcs_memd_interpolate_part_charge_from_grad(fcs_float rel_x, fcs_float *grad, fcs_float *rho)
{
    fcs_int i, k, l, m, index;
    fcs_int grad_ind;
    
    fcs_int help_index[3];
    fcs_float help_x;
    
    help_x = 1. - rel_x;     /* relative pos. w.r.t. first */  
    
    help_index[0] = 4;
    help_index[1] = 2; 
    help_index[2] = 1;
    
    grad_ind = 0;
    index = 0;
    for(i=0;i<8;i++) rho[i] = 0.;
    
    for(k=0;k<2;k++){   /* jumps from x- to x+ */
        for(l=0;l<2;l++){  /* jumps from y- to y+ */
            for(m=0;m<2;m++){ /* jumps from z- to z+ */
                // without q!!!
                if(k==0) rho[index] += - help_x * grad[grad_ind];
                else {
                    rho[index] += - rel_x * grad[grad_ind%4];
                }
                
                grad_ind ++;
                index+=help_index[2];
                help_index[2]=-help_index[2];
            }
            index+=help_index[1];
            help_index[1]=-help_index[1];
        }
        index+=help_index[0];
        help_index[0]=-help_index[0];
    }
}

/** new self force correction with lattice Green's function.
 @param p Particle pointer
 */
void ifcs_memd_calc_part_self_force(memd_struct* memd, memd_particle *p)
{
    fcs_int i, j, k, ip=0;
    fcs_float self, temp;
    static fcs_int help_index[SPACE_DIM];
    
    fcs_int dir1, dir2, d, grad_ind;
    fcs_int l, m, index, temp_ind;
    fcs_float grad[24];
    fcs_float grad2[24];
    fcs_float rho[8];
    fcs_float *force = p->f;
    fcs_float rel = 0.;
    
    help_index[0] = 12;
    help_index[1] = 6;
    help_index[2] = 3; 
    
    ifcs_memd_calc_charge_gradients(&rel, p->q, &grad[ip]);
    ifcs_memd_interpolate_part_charge_from_grad(rel, grad, rho);
    index = 0;
    grad_ind = 0;
    FOR3D(d) {
        ifcs_memd_calc_directions(d, &dir1, &dir2);
        for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
            for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */          
                
                temp_ind = index + d;
                grad2[temp_ind] = grad[grad_ind];
                grad2[temp_ind + help_index[d]] = -grad[grad_ind];
                
                grad_ind++;
                index+=help_index[dir1];
                help_index[dir1]=-help_index[dir1];
            }
            index+=help_index[dir2];
            help_index[dir2]=-help_index[dir2];
            
        }
    }
    
    FOR3D(k) {
        self = 0.; 
        for(i=0;i<8;i++) {
            
            temp = rho[i]*grad2[i*SPACE_DIM + k];
            self += memd->parameters.alpha[i][i] * temp;
            
            for(j=i+1;j<8;j++) {
                temp = rho[i]*grad2[j*SPACE_DIM + k] + rho[j]*grad2[i*SPACE_DIM + k];
                self += memd->parameters.alpha[i][j] * temp;
            }
        }
        force[k] += 2. * self; 
    }
}


/** For each particle P, calculates self energy influence
 with direct Greens function, assuming constant permittivity
 on each lattice site.
 @param P Particle pointer
 */
void ifcs_memd_calc_self_influence(memd_struct* memd, memd_particle* P)
{
    fcs_int k;
    //	fcs_int ix, iy, iz;
    fcs_float invasq = memd->parameters.inva*memd->parameters.inva;
    fcs_float position[SPACE_DIM], left_down_position[SPACE_DIM], relative_position[SPACE_DIM];
    //	fcs_int index = 0;
    //	fcs_int globalindex = 0;
    fcs_float local_force[SPACE_DIM];
    fcs_float particle_charge = P->q;
	
    /* calculate position in cell, normalized to lattice size: */
    FOR3D(k) {
        position[k]           = (P->r[k] - memd->lparams.left_down_position[k]) * memd->parameters.inva;
        left_down_position[k] = floor(position[k]);
        relative_position[k]  = position[k] - left_down_position[k];
        local_force[k] = 0.0;
    }
	
	
    /* Copy permittivity values to the mini-lattice: */
    /*
     for (iz=0;iz<zmax;iz++) {
     for (iy=0;iy<ymax;iy++) {
     for (ix=0;ix<xmax;ix++) {
     index = (iz + zmax*(iy + ymax*ix));
     globalindex = ifcs_memd_get_linear_index((left_down_position[0]+ix),
     (left_down_position[1]+iy),
     (left_down_position[2]+iz), memd->lparams.dim);
     self_influence_correction[index].permittivity = lattice[globalindex].permittivity;
     self_influence_correction[index].charge = 0.0; //background_charge;
     FOR3D(k) {
     charge_gradient[index][k] = 0.0;
     local_E_field[index][k] = 0.0;
     }
     }
     }
     }
     */
	
	
    /* Calculate self-force directly: */	
    FOR3D(k) {
        local_force[k] = pow(((0.5 - relative_position[k])*SELF_FACTOR_1), 3.0) +
        (0.5 - relative_position[k]) * SELF_FACTOR_2;
        local_force[k] *= memd->parameters.bjerrum;
        local_force[k] *= particle_charge * particle_charge;
        local_force[k] *= invasq;
    }
    // memd->parameters.prefactor = 3.544907702 = sqrt(4. * M_PI * memd->parameters.bjerrum * temperature);
    // Fehlender Faktor: 4*pi=12.5663706
    // a^3/c = 2.58448064965803
	
    /* Correct force: */
    FOR3D(k) {
        P->f[k] -= local_force[k];
    }
}

/** Calculate the actual force from the E-Field
 @param p      Particle pointer
 @param index  index of the lattice site
 @param grad   charge gradient
 */
void ifcs_memd_calc_part_link_forces(memd_struct* memd, memd_particle *p, fcs_int index, fcs_float *grad)
{
    static fcs_int init = 1;
    static fcs_int help_index[SPACE_DIM];
    fcs_int ind_grad, j;
    fcs_int dir1, dir2;
    /*  fcs_int* anchor_neighb; */
    fcs_int l,m;
    fcs_float local_force[SPACE_DIM];
	
    if(init) {
        t_site* anchor_site;
        anchor_site = &memd->lattice[index];
        FOR3D(j) {
            help_index[j] = ifcs_memd_get_offset(memd->lattice[memd->neighbor[index][j]].r[j], anchor_site->r[j], j, memd->lparams.dim);
        }
        init = 0;
    }
	
    FOR3D(j){
        local_force[j] = 0.;
    }
	
    ind_grad = 0; 
	
    FOR3D(j) {
        ifcs_memd_calc_directions(j, &dir1, &dir2);
		
        for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
            for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */   
                local_force[j] += -grad[ind_grad]*memd->Dfield[3*index+j];
                //fprintf(stderr, "charge_gradient %d: %1.9f\n", ind_grad, grad[ind_grad]);
                ind_grad++;
                index += help_index[dir1];
                help_index[dir1] = -help_index[dir1];
            }
            index += help_index[dir2];
            help_index[dir2] = -help_index[dir2];
        }
    }  
    /* Attention! Here, the INTERLACING is done! */
    FOR3D(j){
        p->f[j] += memd->parameters.prefactor * local_force[j];
    }
}


/** Public function.
 Calculates the actual force on each particle
 by calling all other needed functions (except
 for fcs_memd_propagate_B_field) */
void fcs_memd_calc_forces(memd_struct* memd)
{ 
    memd_cell *cell;
    static fcs_int init = 1;
    memd_particle *p;
    fcs_int i, c, np, d, index, Npart, ip; 
    fcs_float q;
    /* position of a particle in local lattice units */
    fcs_float pos[SPACE_DIM];
    /* index of first assignment lattice point */
    fcs_int first[3];
    /* charge gradient (number of neighbor sites X number of dimensions) */
    static fcs_float *grad;
	
    if(init) grad = (fcs_float *) realloc(grad, 12*memd->parameters.n_part*sizeof(fcs_float));
	
    /* Hopefully only needed for Yukawa: */
    ifcs_memd_update_charge_gradients(memd, grad);
	
    if(!init) {
        ifcs_memd_couple_current_to_Dfield(memd);
        ifcs_memd_add_transverse_field(memd, memd->parameters.time_step);  
    }
    else init = 0;
	
    ip = 0;
    for (c = 0; c < memd->local_cells.n; c++) {
        cell = memd->local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i=0; i<np; i++) { 
            q = p[i].q;
            if( abs(q) > 1.0e-5 ) {
                FOR3D(d) {
                    pos[d]   = (p[i].r[d] - memd->lparams.left_down_position[d])* memd->parameters.inva;
                    first[d] = (fcs_int) pos[d];
                }
				
                index = ifcs_memd_get_linear_index(first[0],first[1],first[2],memd->lparams.dim);
                ifcs_memd_calc_part_link_forces(memd, &p[i], index, &grad[ip]);
                ifcs_memd_calc_self_influence(memd, &p[i]);
                
                ip+=12;
            }
        }
    }
}

