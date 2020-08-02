/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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

#include "run.h"
#include "types.h"
#include "communication.h"
#include "tune.h"
#include "common/gridsort/gridsort.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "FCSCommon.h"

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/

/********** near formula ***********/
/* self-energy calculation using the near formula */
static void mmm2d_self_energy(mmm2d_data_struct *d, fcs_float *charges, fcs_int num_particles);

/* pair energy calculation using the near formula */
static fcs_float mmm2d_pair_energy_near(mmm2d_data_struct *d, fcs_float charge_factor, fcs_float dv[3], fcs_float dl);

/* actual near formula energy calculation */
static fcs_float mmm2d_calc_pair_energy_near(mmm2d_data_struct *d, fcs_float disp[3]);

/* near formula force calculation */
static void mmm2d_pair_force_near(mmm2d_data_struct *d, fcs_float charge_factor, fcs_float disp[3], fcs_float dl2, fcs_float dl, fcs_float force[3]);

/********** far formula ***********/
/* force and energy far formula contribution */
static fcs_float mmm2d_pair_interactions_far(mmm2d_data_struct *d, fcs_float *forces);
/* main loop for the force far formula */
static void far_force_contribution(mmm2d_data_struct *d, fcs_int p, fcs_int q, fcs_float *forces);
/* main loop for the energy far formula */
static fcs_float far_energy_contribution(mmm2d_data_struct *d, fcs_int p, fcs_int q);

/* 2 pi |z| code */
static void setup_z_force(mmm2d_data_struct *d);
static void setup_z_energy(mmm2d_data_struct *d);
static void   add_z_force(mmm2d_data_struct *d, fcs_float *forces);
static fcs_float  z_energy(mmm2d_data_struct *d);
/* p=0 per frequency code */
static void setup_P(mmm2d_data_struct *d, fcs_int p, fcs_float omega, fcs_float fac);
static void   add_P_force(mmm2d_data_struct *d, fcs_float *forces);
static fcs_float  P_energy(mmm2d_data_struct *d, fcs_float omega);
/* q=0 per frequency code */
static void setup_Q(mmm2d_data_struct *d, fcs_int q, fcs_float omega, fcs_float fac);
static void   add_Q_force(mmm2d_data_struct *d, fcs_float *forces);
static fcs_float  Q_energy(mmm2d_data_struct *d, fcs_float omega);
/* p,q <> 0 per frequency code */
static void setup_PQ(mmm2d_data_struct *d, fcs_int p, fcs_int q, fcs_float omega, fcs_float fac);
static void   add_PQ_force(mmm2d_data_struct *d, fcs_int p, fcs_int q, fcs_float omega, fcs_float *forces);
static fcs_float  PQ_energy(mmm2d_data_struct *d, fcs_float omega);

/* cache management for far formula precalculated factors */
static void realloc_caches(mmm2d_data_struct *d);
static void prepare_scx_cache(mmm2d_data_struct *d);
static void prepare_scy_cache(mmm2d_data_struct *d);

static fcs_float *block(fcs_float *p, fcs_int index, fcs_int size);
static fcs_float *blwentry(fcs_float *p, fcs_int index, fcs_int e_size);
static fcs_float *abventry(fcs_float *p, fcs_int index, fcs_int e_size);

/********** dielectric layers ***********/
static fcs_float dielectric_layers_energy_contribution(mmm2d_data_struct *d);
static void dielectric_layers_force_contribution(mmm2d_data_struct *d, fcs_float *forces);

/* dealing with the image contributions from far outside the simulation box */
/* gather the informations for the far away image charges */
static void gather_image_contributions(mmm2d_data_struct *d, fcs_int e_size);
/* clear the image contributions if there is no dielectric contrast and no image charges */
static void clear_image_contributions(mmm2d_data_struct *d, fcs_int size);
/* spread the top/bottom sums */
static void distribute(mmm2d_data_struct *d, fcs_int e_size, fcs_float fac);

/* vector operations */
/* pdc = 0 */
static void clear_vec(fcs_float *pdc, fcs_int size);
/* pdc_d = pdc_s */
static void copy_vec(fcs_float *pdc_d, fcs_float *pdc_s, fcs_int size);
/* pdc_d = pdc_s1 + pdc_s2 */
static void add_vec(fcs_float *pdc_d, fcs_float *pdc_s1, fcs_float *pdc_s2, fcs_int size);
/* pdc_d = scale*pdc_s1 + pdc_s2 */
static void addscale_vec(fcs_float *pdc_d, fcs_float scale, fcs_float *pdc_s1, fcs_float *pdc_s2, fcs_int size);
/* pdc_d = scale*pdc */
static void scale_vec(fcs_float scale, fcs_float *pdc, fcs_int size);

/*  */
static void layered_displacement_vector(mmm2d_data_struct *d, fcs_float x1, fcs_float y1, fcs_float z1, fcs_float x2, fcs_float y2, fcs_float z2, fcs_float disp[3]);


/***************************************************/
/* IMPLEMENTATION OF THE PUBLIC FUNCTION */
/***************************************************/
void mmm2d_run(void* rd,
        fcs_int num_particles,
        fcs_int max_num_particles,
        fcs_float *positions,
        fcs_float *charges,
        fcs_float *forces,
        fcs_float *potentials) {
  /* Here we assume, that the method is tuned and that all parameters are valid */
  mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
  
  /* calculate self energy */
  mmm2d_self_energy(d, charges, num_particles);
  //printf("rank %d, self_energy %f\n", d->comm.rank, d->self_energy);
  
  /* decompose system */
  fcs_int local_num_particles;
  fcs_int local_num_ghost_particles;
  fcs_float *local_positions, *local_ghost_positions;
  fcs_float *local_charges, *local_ghost_charges;
  fcs_gridsort_index_t *local_indices, *local_ghost_indices;
  fcs_float box_base[3] = {0.0, 0.0, 0.0 };
  fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
  fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
  fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };
  fcs_gridsort_t gridsort;
  fcs_int zslices_nparticles[d->layers_per_node], zslices_ghost_nparticles[2];
  fcs_int i;
  
  fcs_gridsort_create(&gridsort);
  
  fcs_gridsort_set_system(&gridsort, box_base, box_a, box_b, box_c, NULL);

  fcs_gridsort_set_zslices(&gridsort, d->layers_per_node, 1);

  fcs_gridsort_set_particles(&gridsort, num_particles, max_num_particles, positions, charges);
  
  MPI_Barrier(d->comm.mpicomm);
  fprintf(stderr,"mmm2d_run, %d\n", d->comm.rank);
  //printf("layer_h: %f\n",d->layer_h);
  
  //printf("  calling fcs_gridsort_sort_forward()...\n");
  
  fcs_gridsort_sort_forward(&gridsort, 0.0, d->comm.mpicomm);
  
  MPI_Barrier(d->comm.mpicomm);
  //printf("  returning from fcs_gridsort_sort_forward().\n");
  
  fcs_gridsort_separate_ghosts(&gridsort);
  //printf("  rank %d, separate ghosts\n", d->comm.rank);
  
  //printf("  calling separate zslices\n");
  fcs_gridsort_separate_zslices(&gridsort, &zslices_ghost_nparticles[0], zslices_nparticles, &zslices_ghost_nparticles[1]);
  
  //printf("  rank %d, separate zslices: layers_per_node %d... ", d->comm.rank, d->layers_per_node);
  /*
  for (i=0; i<d->layers_per_node; i++) {
    printf("rank %d, layer %d: %d parts, ", d->comm.rank, i+1, zslices_nparticles[i]);
  }
  printf("\n");
  printf("rank %d, low ghost layer: %d parts\n", d->comm.rank, zslices_ghost_nparticles[0]);
  printf("rank %d, high ghost layer: %d parts\n", d->comm.rank, zslices_ghost_nparticles[1]);
  */
  
  MPI_Barrier(d->comm.mpicomm);
  
  fcs_gridsort_get_real_particles(&gridsort, &local_num_particles, &local_positions, &local_charges, &local_indices);
  
  d->n_localpart=local_num_particles;
  d->local_charges=local_charges;
  d->local_positions=local_positions;
  d->zslices_nparticles=zslices_nparticles;
  
  MPI_Barrier(d->comm.mpicomm);
  fprintf(stderr, "rank %d, local particles: %d\n", d->comm.rank, d->n_localpart);
  /*
  if (d->n_localpart>0) {
   printf("rank %d, test part[0]: local pos %f, global pos %f\n", d->comm.rank, local_positions[0], (d->local_positions)[0]);
  }
  */
  
  fcs_gridsort_get_ghost_particles(&gridsort, &local_num_ghost_particles, &local_ghost_positions, &local_ghost_charges, &local_ghost_indices);
  
  MPI_Barrier(d->comm.mpicomm);
  fprintf(stderr, "rank %d, ghost particles: %d\n", d->comm.rank, zslices_nparticles[0]+zslices_nparticles[1]);
  
  /* allocate local forces */
  fprintf(stderr, "rank %d, allocate local containers\n", d->comm.rank);
  fcs_float *local_forces = NULL;
  //fcs_float *local_potentials = NULL;
  if (forces != NULL) {
    local_forces = malloc(sizeof(fcs_float)*3*d->n_localpart);
    for (i=0; i<3*local_num_particles; i++) local_forces[i]=0.;
  }
  /*fcs_float energy, total_energy=0.;
  if (potentials != NULL || d->require_total_energy) {
    local_potentials = malloc(sizeof(fcs_float)*d->n_localpart);
    for (i=0; i<d->n_localpart; i++) local_potentials[i]=0.;
  }*/
  
  fcs_int j, c, ci, cj, di, dj, np, npb, offset=0, offsetb=0;
  fcs_float disp[3]={0., 0., 0.}, displ;
  
  /* near formula force calculation */
  fprintf(stderr, "rank %d, near formula force calculation\n", d->comm.rank);
  if (forces != NULL) {
    fcs_float force[3]={0., 0., 0.};
    for (c = 0; c < d->layers_per_node; c++) {
      np   = d->zslices_nparticles[c];
      ///@TODO: optimize indexes for minimum calculations
      for(i = 0; i < np; i++) {
        ci=i+offset;
        di=3*ci;
        //particles in the same slice
        for(j = i+1; j < np; j++) {
          cj=j+offset;
          dj=3*cj;
          layered_displacement_vector(d, local_positions[di], local_positions[di+1], local_positions[di+2], local_positions[dj], local_positions[dj+1], local_positions[dj+2], disp);
          displ=disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
          mmm2d_pair_force_near(d, local_charges[ci]*local_charges[cj], disp, displ, sqrt(displ), force);
          local_forces[di]+=force[0];
          local_forces[di+1]+=force[1];
          local_forces[di+2]+=force[2];
          local_forces[dj]-= force[0];
          local_forces[dj+1]-= force[1];
          local_forces[dj+2]-= force[2];
        }
        
        //particles in the bottom neighbor
        if (c>0) {
          ////bottom neighbor is a real slice
          ///@TODO: test with np=2, n_layers=2, max_posz=0.049
          npb=d->zslices_nparticles[c-1];
          for(j = 0; j < npb; j++) {
            cj=j+offsetb;
            dj=3*cj;
            layered_displacement_vector(d, local_positions[di], local_positions[di+1], local_positions[di+2], local_positions[dj], local_positions[dj+1], local_positions[dj+2], disp);
            displ=disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
            mmm2d_pair_force_near(d, local_charges[ci]*local_charges[cj], disp, displ, sqrt(displ), force);
            local_forces[di]+=force[0];
            local_forces[di+1]+=force[1];
            local_forces[di+2]+=force[2];
            local_forces[dj]-= force[0];
            local_forces[dj+1]-= force[1];
            local_forces[dj+2]-= force[2];
          }
        } else if(c==0) {
          ////bottom slice, my bottom is a ghost slice
          for(j = 0; j < zslices_ghost_nparticles[0]; j++) {
            dj=3*j;
            layered_displacement_vector(d, local_positions[di], local_positions[di+1], local_positions[di+2], local_ghost_positions[dj], local_ghost_positions[dj+1], local_ghost_positions[dj+2], disp);
            //local_ghost_positions, local_ghost_charges
            displ=disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
            mmm2d_pair_force_near(d, local_charges[ci]*local_ghost_charges[j], disp, displ, sqrt(displ), force);
            local_forces[di]+=force[0];
            local_forces[di+1]+=force[1];
            local_forces[di+2]+=force[2];
          }
        }
        
        if(c==d->layers_per_node-1) {
          ////bottom slice, my bottom neighbor is a ghost slice
          for(j = 0; j < zslices_ghost_nparticles[1]; j++) {
            cj=j+zslices_ghost_nparticles[0];
            dj=3*cj;
            layered_displacement_vector(d, local_positions[di], local_positions[di+1], local_positions[di+2], local_ghost_positions[dj], local_ghost_positions[dj+1], local_ghost_positions[dj+2], disp);
            //local_ghost_positions, local_ghost_charges
            displ=disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
            mmm2d_pair_force_near(d, local_charges[ci]*local_ghost_charges[cj], disp, displ, sqrt(displ), force);
            local_forces[di]+=force[0];
            local_forces[di+1]+=force[1];
            local_forces[di+2]+=force[2];
          }
        }
      }
      offsetb=offset;
      offset += np;
    }
    
  }
  
  /* near formula energy calculations */
  fprintf(stderr, "rank %d, near formula energy calculation\n", d->comm.rank);
  if (d->require_total_energy) {
    offset=0;
    offsetb=0;
    for (c = 0; c < d->layers_per_node; c++) {
      np   = zslices_nparticles[c];
      for(i = 0; i < np; i++) {
        ci=i+offset;
        di=3*ci;
        /// cell itself
        for(j = i+1; j < np; j++) {
          cj=j+offset;
          dj=3*cj;
          layered_displacement_vector(d, local_positions[di], local_positions[di+1], local_positions[di+2], local_positions[dj], local_positions[dj+1], local_positions[dj+2], disp);
          displ=disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
          d->total_energy+=mmm2d_pair_energy_near(d, local_charges[ci]*local_charges[cj], disp, sqrt(displ));
        }
      }
      /// bottom neighbor
      if (c>0) {
        npb = zslices_nparticles[c-1];
        for(i = 0; i < np; i++) {
          ci=i+offset;
          di=3*ci;
          for(j = 0; j < npb; j++) {
            cj=j+offsetb;
            dj=3*cj;
            layered_displacement_vector(d, local_positions[di], local_positions[di+1], local_positions[di+2], local_positions[dj], local_positions[dj+1], local_positions[dj+2], disp);
            displ=disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
            d->total_energy+=mmm2d_pair_energy_near(d, local_charges[ci]*local_charges[cj], disp, sqrt(displ));
          }
        }
        offsetb=offset;
      } else {
        ////bottom slice, my bottom neighbor is a ghost slice
        npb = zslices_ghost_nparticles[0];
        for(i = 0; i < np; i++) {
          ci=i+offset;
          di=3*ci;
          for(j = 0; j < npb; j++) {
            dj=3*j;
            layered_displacement_vector(d, local_positions[di], local_positions[di+1], local_positions[di+2], local_ghost_positions[dj], local_ghost_positions[dj+1], local_ghost_positions[dj+2], disp);
            displ=disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
            d->total_energy+=mmm2d_pair_energy_near(d, local_charges[ci]*local_ghost_charges[j], disp, sqrt(displ));
          }
        }
      }
      /*
      if(c == d->layers_per_node-1){
        ////top slice, interaction with top ghost slice must be calculated
        npb = zslices_ghost_nparticles[1];
        for(i = 0; i < np; i++) {
          ci=i+offset;
          di=3*ci;
          for(j = 0; j < npb; j++) {
            cj=j+zslices_ghost_nparticles[0];
            dj=3*cj;
            layered_displacement_vector(d, local_positions[di], local_positions[di+1], local_positions[di+2], local_ghost_positions[dj], local_ghost_positions[dj+1], local_ghost_positions[dj+2], disp);
            displ=disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
            energy=mmm2d_pair_energy_near(d, local_charges[ci]*local_ghost_charges[cj], disp, sqrt(displ));
            local_potentials[ci]+=energy;
          //
          }
        }
      }
      */
      offset+=np;
    }
  }
  
  /* far formula force and energy calculations */
  fprintf(stderr, "rank %d, far formula calculations\n", d->comm.rank);
  /* allocate far formmula caches */
  ///@TODO: only if far formula is needed
  realloc_caches(d);
  //printf("rank %d, caches reallocated\n", d->comm.rank);
  if (d->require_total_energy) {
    d->total_energy+=mmm2d_pair_interactions_far(d, local_forces);
  } else {
    mmm2d_pair_interactions_far(d, local_forces);
  }
  
  fprintf(stderr, "rank %d, dielectric contrasts contributions\n", d->comm.rank);
  /* forces from dielectric layers */
  if (forces != NULL && d->dielectric_contrast_on) dielectric_layers_force_contribution(d, local_forces);
  
  /* potentials from dielectric layers */
  if (d->require_total_energy && d->dielectric_contrast_on) d->total_energy+=dielectric_layers_energy_contribution(d);
  
  /* add total forces */
  /*
  if (forces != NULL) {
    for (i=0; i<d->n_localpart; i++) {
      di=3*i;
      //printf("local_forces: %d: %f %f %f\n", i, local_forces[di], local_forces[di+1], local_forces[di+2]);
      forces[di]=local_forces[di]; ///@TODO: add, don't set?
      forces[di+1]=local_forces[di+1];
      forces[di+2]=local_forces[di+2];
    }
  }
  */
  
  /* clean up and finish */
  fcs_gridsort_set_sorted_results(&gridsort, local_num_particles, local_forces, NULL);
  fcs_gridsort_set_results(&gridsort, max_num_particles, forces, NULL);
  fcs_gridsort_sort_backward(&gridsort, d->comm.mpicomm);
  
  fcs_gridsort_free(&gridsort);
  
  fcs_gridsort_destroy(&gridsort);
  
  fcs_float total_energy;
  //MPI_Reduce(&d->total_energy, &total_energy, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->comm.mpicomm);
  if (d->require_total_energy){
    MPI_Allreduce(&d->total_energy, &total_energy, 1, FCS_MPI_FLOAT, MPI_SUM, d->comm.mpicomm);
    d->total_energy=total_energy;
  }
  return;
}

/***************************************************/
/* IMPLEMENTATION OF THE PRIVATE FUNCTIONS */
/***************************************************/

/********** near formula ***********/

/* self-energy calculation using the near formula */
static void mmm2d_self_energy(mmm2d_data_struct *d, fcs_float *charges, fcs_int num_particles)
{
  fcs_int c;
  fcs_float seng, dv[3] = {0, 0, 0};
  seng = 0.;
  /* this one gives twice the real self energy, as it is used
     in the far formula which counts everything twice and in
     the end divides by two*/
  for (c=0; c < num_particles; c++) {
    seng += charges[c]*charges[c];
  }
  d->self_energy = seng * mmm2d_calc_pair_energy_near(d, dv);
}

/* pair energy calculation using the near formula */
static fcs_float mmm2d_pair_energy_near(mmm2d_data_struct *d, fcs_float charge_factor,
                                 fcs_float dv[3], fcs_float dl)
{
  fcs_float eng;
  if (charge_factor != 0.0) {
    eng = mmm2d_calc_pair_energy_near(d, dv);
    return charge_factor*(eng + 1./dl);
  }
  return 0.0;
}

/* actual near formula energy calculation */
fcs_float mmm2d_calc_pair_energy_near(mmm2d_data_struct *d, fcs_float disp[3])
{
  fcs_float eng;
  fcs_float z2     = disp[2]*disp[2];
  fcs_float rho2   = disp[1]*disp[1] + z2;
  
  ///* the ux is multiplied in below
  eng = -2*log(4*M_PI*d->uy*d->box_l[0]);
  
  ///* Bessel sum
  {
    fcs_int p, l;
    fcs_float k0Sum;
    fcs_float freq;
    fcs_float rho_l, ypl;
    fcs_float c;
    
    for (p = 1; p < d->besselCutoff.n; p++) {
      k0Sum  = 0;
      
      freq = MMM_COMMON_C_2PI*d->ux*p;
      
      for (l = 1; l < d->besselCutoff.e[p-1]; l++) {
        ypl   = disp[1] + l*d->box_l[1];
        rho_l = sqrt(ypl*ypl + z2);
        k0Sum  += mmm_K0(freq*rho_l);
        
        ypl   = disp[1] - l*d->box_l[1];
        rho_l = sqrt(ypl*ypl + z2);
        k0Sum  += mmm_K0(freq*rho_l);
      }
      
      ///* the ux is multiplied in to bessel, complex and psi at once, not here
      c = 4*cos(freq*disp[0]);
      eng += c*k0Sum;
    }
  }
  
  ///* complex sum
  {
    fcs_float zeta_r, zeta_i;
    fcs_float zet2_r, zet2_i;
    fcs_float ztn_r,  ztn_i;
    fcs_float tmp_r;
    fcs_int end, n;
    
    ztn_r = zeta_r = d->uy*disp[2];
    ztn_i = zeta_i = d->uy*disp[1];
    
    zet2_r = zeta_r*zeta_r - zeta_i*zeta_i;
    zet2_i = 2*zeta_r*zeta_i;
    
    ztn_r = zet2_r;
    ztn_i = zet2_i;
    
    end = (fcs_int)ceil(MMM2D_COMPLEX_FAC*d->uy2*rho2);
    if (end > MMM2D_COMPLEX_STEP) {
      end = MMM2D_COMPLEX_STEP;
      fprintf(stderr, "MMM2D: some particles left the assumed slab, precision might be lost\n");
    }
    end = d->complexCutoff[end];
    for (n = 1; n <= end; n++) {
      eng -= d->box_l[1]/(2*n)*d->bon.e[n-1]*ztn_r;
      
      tmp_r = ztn_r*zet2_r - ztn_i*zet2_i;
      ztn_i = ztn_r*zet2_i + ztn_i*zet2_r;
      ztn_r = tmp_r;
    }
  }
  
  ///* psi sum
  {
    fcs_int n;
    fcs_float add;
    fcs_float uxx = d->ux*disp[0];
    fcs_float uxrho2 = d->ux2*rho2;
    fcs_float uxrho_2n;
    
    ///* n = 0 inflicts only Fx and pot
    ///* one ux is multiplied in to bessel, complex and psi at once, not here
    eng -= mmm_mod_psi_even(d->polTaylor, 0, uxx);
    
    uxrho_2n = uxrho2;
    for (n = 1; n < (d->polTaylor)->n_modPsi; n++) {
      add = uxrho_2n*mmm_mod_psi_even(d->polTaylor, n, uxx);
      eng -= add;
      if (fabs(add) < d->part_error)
        break;
      uxrho_2n *= uxrho2;
    }
  }
  
  eng *= d->ux;
  
  ///* explicitly added potentials r_{-1,0} and r_{1,0}
  {
    fcs_float cx   = disp[0] + d->box_l[0];
    fcs_float rinv = sqrt(1.0/(cx*cx + rho2));
    eng += rinv;
    
    cx   = disp[0] - d->box_l[0];
    rinv = sqrt(1.0/(cx*cx + rho2));
    eng += rinv;
  }
  return eng;
}

/* near formula force calculation */
static void mmm2d_pair_force_near(mmm2d_data_struct *d, fcs_float charge_factor, fcs_float dist[3], fcs_float dl2, fcs_float dl, fcs_float force[3]) {
  //printf("bessel cutoff %d %f %f\ndist %f %f %f\n", d->besselCutoff.n, MMM_COMMON_C_2PI, d->ux, dist[0], dist[1], dist[2]);
  fcs_float F[3];
  fcs_float z2   = dist[2]*dist[2];
  fcs_float rho2 = dist[1]*dist[1] + z2;
  fcs_int i;
  
  if (charge_factor == 0.0)
      return;
  
  F[0] = F[1] = F[2] = 0;
  
  ///* Bessel sum
  {
    fcs_int p, l;
    fcs_float k0, k1;
    fcs_float k0Sum, k1ySum, k1Sum;
    fcs_float freq;
    fcs_float rho_l, ypl;
    fcs_float c, s;

    for (p = 1; p < d->besselCutoff.n; p++) {
      k0Sum  = 0;
      k1ySum = 0;
      k1Sum  = 0;

      freq = MMM_COMMON_C_2PI*d->ux*p;

      for (l = 1; l < d->besselCutoff.e[p-1]; l++) {
        ypl   = dist[1] + l*d->box_l[1];
        rho_l = sqrt(ypl*ypl + z2);
  #ifdef BESSEL_MACHINE_PREC
        k0 = mmm_K0(freq*rho_l);
        k1 = mmm_K1(freq*rho_l);
  #else
        mmm_LPK01(freq*rho_l, &k0, &k1);
  #endif
        k1 /= rho_l;
        k0Sum  += k0;
        k1Sum  += k1;
        k1ySum += k1*ypl;

        ypl   = dist[1] - l*d->box_l[1];
        rho_l = sqrt(ypl*ypl + z2);
  #ifdef BESSEL_MACHINE_PREC
        k0 = mmm_K0(freq*rho_l);
        k1 = mmm_K1(freq*rho_l);
  #else
        mmm_LPK01(freq*rho_l, &k0, &k1);
  #endif
        k1 /= rho_l;
        k0Sum  += k0;
        k1Sum  += k1;
        k1ySum += k1*ypl;
      }

      ///* the ux is multiplied in to bessel, complex and psi at once, not here 
      c = 4*freq*cos(freq*dist[0]);
      s = 4*freq*sin(freq*dist[0]);
      F[0] +=      s*k0Sum;
      F[1] +=      c*k1ySum;
      F[2] += dist[2]*c*k1Sum;
    }
    //printf(" bessel force %f %f %f\n", F[0], F[1], F[2]);
  }
  
  ///* complex sum
  {
    fcs_float zeta_r, zeta_i;
    fcs_float zet2_r, zet2_i;
    fcs_float ztn_r,  ztn_i;
    fcs_float tmp_r;
    fcs_int end, n;
    
    ztn_r = zeta_r = d->uy*dist[2];
    ztn_i = zeta_i = d->uy*dist[1];
    zet2_r = zeta_r*zeta_r - zeta_i*zeta_i;
    zet2_i = 2*zeta_r*zeta_i;
    
    end = (fcs_int)ceil(MMM2D_COMPLEX_FAC*d->uy2*rho2);
    if (end > MMM2D_COMPLEX_STEP) {
      //printf("end: %d, complex step %d, complex fac: %e, uy2: %e, rho2: %e\n", end, MMM2D_COMPLEX_STEP, MMM2D_COMPLEX_FAC, d->uy2, rho2);
      end = MMM2D_COMPLEX_STEP;
      fprintf(stderr, "MMM2D: some particles left the assumed slab, precision might be lost\n");
    }
    if (end < 0) {
      //char *errtxt = runtime_error(100);
      fprintf(stderr, "MMM2D: distance was negative, coordinates probably out of range!!\n");
      end = 0;
    }
    end = d->complexCutoff[end];
    
    for (n = 0; n < end; n++) {
      F[1] -= d->bon.e[n]*ztn_i;
      F[2] += d->bon.e[n]*ztn_r;
      
      tmp_r = ztn_r*zet2_r - ztn_i*zet2_i;
      ztn_i = ztn_r*zet2_i + ztn_i*zet2_r;
      ztn_r = tmp_r;
    }
    //printf("complex force %f %f %f %d\n", F[0], F[1], F[2], end);
  }
    
  ///* psi sum
  {
    fcs_int n;
    fcs_float uxx = d->ux*dist[0];
    fcs_float uxrho2 = d->ux2*rho2;
    fcs_float uxrho_2n, uxrho_2nm2; ///* rho^{2n-2} 
    fcs_float mpe, mpo;
    
    ///* n = 0 inflicts only Fx and pot 
    ///* one ux is multiplied in to bessel, complex and psi at once, not here 
    F[0] += d->ux*mmm_mod_psi_odd(d->polTaylor, 0, uxx);
    
    uxrho_2nm2 = 1.0;
    for (n = 1;n < (d->polTaylor)->n_modPsi; n++) {
      mpe    = mmm_mod_psi_even(d->polTaylor, n, uxx);
      mpo    = mmm_mod_psi_odd(d->polTaylor, n, uxx);
      uxrho_2n = uxrho_2nm2*uxrho2;
      
      F[0] +=     d->ux *uxrho_2n  *mpo;
      F[1] += 2*n*d->ux2*uxrho_2nm2*mpe*dist[1];
      F[2] += 2*n*d->ux2*uxrho_2nm2*mpe*dist[2];
      
      ///* y < rho => ux2*uxrho_2nm2*dist[1] < ux*uxrho_2n 
      if (fabs(2*n*d->ux*uxrho_2n*mpe) < d->part_error)
        break;
      
      uxrho_2nm2 = uxrho_2n;
    }
  }
  
  for (i = 0; i < 3; i++)
    F[i] *= d->ux;
  
  ///* explicitly added potentials r_{-1,0} and r_{1,0}
  {
    fcs_float cx    = dist[0] + d->box_l[0];
    fcs_float rinv2 = 1.0/(cx*cx + rho2), rinv = sqrt(rinv2);
    fcs_float rinv3 = rinv*rinv2;
    F[0] +=   cx*rinv3;
    F[1] += dist[1]*rinv3;
    F[2] += dist[2]*rinv3;
    
    cx   = dist[0] - d->box_l[0];
    rinv2 = 1.0/(cx*cx + rho2); rinv = sqrt(rinv2);
    rinv3 = rinv*rinv2;
    F[0] +=   cx*rinv3;
    F[1] += dist[1]*rinv3;
    F[2] += dist[2]*rinv3;
    
    rinv3 = 1/(dl2*dl);
    F[0] += dist[0]*rinv3;
    F[1] += dist[1]*rinv3;
    F[2] += dist[2]*rinv3;
    
    //printf("explcit force %f %f %f\n", F[0], F[1], F[2]);
  }
  
  for (i = 0; i < 3; i++)
    force[i] = charge_factor*F[i];
  
    //printf("force final: %f %f %f\n", force[0], force[1], force[2]);
}

/********** far formula ***********/
static fcs_float mmm2d_pair_interactions_far(mmm2d_data_struct *d, fcs_float *forces)
{
  //if(f)  printf("doy fuerza\n");
  //if(e) printf("doy energia\n");
  //printf("far_cut: %f\n", mmm2d_params.far_cut);
  fcs_float eng;
  fcs_int f, e, p, q;
  fcs_float R, dR, q2;
  fcs_int *undone;
  
  if (forces==NULL) f=0; else f=1;
  if (!d->require_total_energy) e=0; else e=1;
  
  // It's not really far...
  eng = e ? d->self_energy : 0.;
  
  if (d->far_cut == 0.0) {
    return 0.5*eng;
  }
  
  undone = malloc((d->n_scxcache + 1)*sizeof(fcs_int));
  
  fprintf(stderr, "rank %d, preparing sc caches\n", d->comm.rank);
  prepare_scx_cache(d);
  prepare_scy_cache(d);
  
  //printf("rank %d, pair_interactions prepare caches done\n", d->comm.rank);
  
  /* complicated loop. We work through the p,q vectors in rings
  from outside to inside to avoid problems with cancellation */
  
  // up to which q vector we have to work
  for (p = 0; p <= d->n_scxcache; p++) {
    if (p == 0)
      q =  d->n_scycache;
    else {
      q2 = d->far_cut2 - (d->ux*(p - 1))*(d->ux*(p - 1));
      if (q2 > 0)
        q = 1 + d->box_l[1]*(int)ceil(sqrt(q2));
      else
        q = 1;
      // just to be on the safe side... 
      if (q > d->n_scycache) q = d->n_scycache;
    }
    undone[p] = q;
  }

  dR = -log(MMM2D_FARRELPREC)/MMM_COMMON_C_2PI*d->uz;
  
  //printf("************ 1\n");
  for(R = d->far_cut; R > 0; R -= dR) {
    for (p = d->n_scxcache; p >= 0; p--) {
      for (q = undone[p]; q >= 0; q--) {
        if (d->ux2*p*p  + d->uy2*q*q < R*R)
          break;
        if (f)
          far_force_contribution(d, p, q, forces);
          //printf("************ 1.1\n");
        if (e)
          eng += far_energy_contribution(d, p, q);
          //printf("************ 1.2\n");
      }
      undone[p] = q;
    }
  }
  //printf("************ 2\n");
  // clean up left overs 
  for (p = d->n_scxcache; p >= 0; p--) {
    q = undone[p];
    // fprintf(stderr, "left over %d\n", q);
    for (; q >= 0; q--) {
      // printf("xxxxx %d %d\n", p, q);
      if (f)
        far_force_contribution(d, p, q, forces);
      if (e)
        eng += far_energy_contribution(d, p, q);
    }
  }
  free(undone);
  return 0.5*eng;
}

static void realloc_caches(mmm2d_data_struct *d)
{
  d->n_scxcache = (fcs_int)(ceil(d->far_cut/d->ux) + 1.);
  d->n_scycache = (fcs_int)(ceil(d->far_cut/d->uy) + 1.);
  d->scxcache = realloc(d->scxcache, d->n_scxcache*d->n_localpart*sizeof(mmm2d_SCCache));
  d->scycache = realloc(d->scycache, d->n_scycache*d->n_localpart*sizeof(mmm2d_SCCache));
  d->partblk   = realloc(d->partblk,  d->n_localpart*8*sizeof(fcs_float));
  d->lclcblk   = realloc(d->lclcblk,  d->n_total_layers*8*sizeof(fcs_float));
  d->gblcblk   = realloc(d->gblcblk,  d->layers_per_node*8*sizeof(fcs_float));
}

static void prepare_scx_cache(mmm2d_data_struct *d)
{
  fcs_int i, ic, freq, o;
  fcs_float pref, arg;
  //printf("prepare_scx_cache: n_sxc: %d, n_local %d\n", d->n_scxcache, d->n_localpart);
  for (freq = 1; freq <= d->n_scxcache; freq++) {
    pref = MMM_COMMON_C_2PI*d->ux*freq;
    o = (freq-1)*d->n_localpart;
    ic = 0;
    for (i = 0; i < d->n_localpart; i++) {
      arg = pref*d->local_positions[i*3];
      d->scxcache[o + ic].s = sin(arg);
      d->scxcache[o + ic].c = cos(arg);
      ic++;
    }
  }
}

static void prepare_scy_cache(mmm2d_data_struct *d)
{
  fcs_int i, ic, freq, o;
  fcs_float pref, arg;

  for (freq = 1; freq <= d->n_scycache; freq++) {
    pref = MMM_COMMON_C_2PI*d->uy*freq;
    o = (freq-1)*d->n_localpart;
    ic = 0;
    for (i = 0; i < d->n_localpart; i++) {
      arg = pref*d->local_positions[i*3+1];
      d->scycache[o + ic].s = sin(arg);
      d->scycache[o + ic].c = cos(arg);
      ic++;
    }
  }
}

/*****************************************************************/
/* far formula main loops */
/*****************************************************************/

static void far_force_contribution(mmm2d_data_struct *d, fcs_int p, fcs_int q, fcs_float *forces)
{
   //printf("rank %d, far force contribution 0, p %d, q %d\n", d->comm.rank, p, q);
   //printf(" rank %d, p: %d, q: %d\n",this_node, p, q);
  fcs_float omega, fac;

  if (q == 0) {
    if (p == 0) {
//printf("rank %d, far force contribution 1*\n", d->comm.rank);
      setup_z_force(d);
//printf("rank %d, far force contribution 2*\n", d->comm.rank);
      if (d->dielectric_contrast_on)
        gather_image_contributions(d, 1);
      else
        clear_image_contributions(d, 1);
//printf("rank %d, far force contribution 3*\n", d->comm.rank);
      distribute(d, 1, 1.);
      add_z_force(d, forces);
//printf("rank %d, far force contribution 4*\n", d->comm.rank);
      //checkpoint("************2piz", 0, 0, 1);
    } else {
      omega = MMM_COMMON_C_2PI*d->ux*p;
      fac = exp(-omega*d->layer_h);
      setup_P(d, p, omega, fac);
      if (d->dielectric_contrast_on)
        gather_image_contributions(d, 2);
      else
        clear_image_contributions(d, 2);
      distribute(d, 2, fac);
      add_P_force(d, forces);
      //checkpoint("************distri p", p, 0, 2);
    }
  } else if (p == 0) {
//printf("rank %d, far force contribution 1**\n", d->comm.rank);
    omega = MMM_COMMON_C_2PI*d->uy*q;
    fac = exp(-omega*d->layer_h);
    setup_Q(d, q, omega, fac);
    if (d->dielectric_contrast_on)
      gather_image_contributions(d, 2);
    else
      clear_image_contributions(d, 2);
    distribute(d, 2, fac);
    add_Q_force(d, forces);
//printf("rank %d, far force contribution 2**\n", d->comm.rank);
    //checkpoint("************distri q", 0, q, 2);
  } else {
//printf("rank %d, far force contribution 1***\n", d->comm.rank);
    omega = MMM_COMMON_C_2PI*sqrt((d->ux*p)*(d->ux*p) + (d->uy*q)*(d->uy*q));
    fac = exp(-omega*d->layer_h);
    setup_PQ(d, p, q, omega, fac);
    if (d->dielectric_contrast_on)
      gather_image_contributions(d, 4);
    else
      clear_image_contributions(d, 4);
//printf("rank %d, far force contribution 2***\n", d->comm.rank);
    distribute(d, 4, fac);
    add_PQ_force(d, p, q, omega, forces);
    //checkpoint("************distri pq", p, q, 4);
  }
}

fcs_float far_energy_contribution(mmm2d_data_struct *d, fcs_int p, fcs_int q)
{
  fcs_float eng;
  //printf("energy contribution\n");
  fcs_float omega, fac;
  
  if (q == 0) {
    if (p == 0) {
      setup_z_energy(d);
      clear_image_contributions(d, 2);
      distribute(d, 2, 1.);
      eng = z_energy(d);
      //checkpoint("E************2piz", 0, 0, 2);
    }
    else {
      omega = MMM_COMMON_C_2PI*d->ux*p;
      fac = exp(-omega*d->layer_h);
      setup_P(d, p, omega, fac);
      if (d->dielectric_contrast_on)
   gather_image_contributions(d, 2);
      else
   clear_image_contributions(d, 2);
      distribute(d, 2, fac);
      eng = P_energy(d, omega);
      //checkpoint("************distri p", p, 0, 2);
    }
  }
  else if (p == 0) {
    omega = MMM_COMMON_C_2PI*d->uy*q;
    fac = exp(-omega*d->layer_h);
    setup_Q(d, q, omega, fac);
    if (d->dielectric_contrast_on)
      gather_image_contributions(d, 2);
    else
      clear_image_contributions(d, 2);
    distribute(d, 2, fac);
    eng = Q_energy(d, omega);
    //checkpoint("************distri q", 0, q, 2);
  }
  else {
    omega = MMM_COMMON_C_2PI*sqrt((d->ux*p)*(d->ux*p) + (d->uy*q)*(d->uy*q));
    fac = exp(-omega*d->layer_h);
    setup_PQ(d, p, q, omega, fac);
    if (d->dielectric_contrast_on)
      gather_image_contributions(d, 4);
    else
      clear_image_contributions(d, 4);
    distribute(d, 4, fac);
    eng = PQ_energy(d, omega);
    //checkpoint("************distri pq", p, q, 4);
  }
  return eng;
}

static void setup_z_force(mmm2d_data_struct *d)
{
  fcs_int np, c, i, ci, offset=0;
  fcs_float pref = MMM_COMMON_C_2PI*d->ux*d->uy;
  fcs_float *lclimgebot=NULL,*lclimgetop=NULL;
  fcs_int e_size=1,size = 2;
  fcs_float e, e_di_l, e_di_h;

  fcs_float fac_imgsum;

  /* in case of metallic boundary conditions on both sides, we get an infinite array,
     which only exists for charge neutral systems. But in this case, we can as well
     not sum up the force array, as the net force per image is 0 */

  if (d->delta_mult != 1.0) {
    fac_imgsum = 1/(1 - d->delta_mult);
  }
  else {
    fac_imgsum = 0;
  }

  if (d->dielectric_contrast_on)
    clear_vec(d->lclimge, size);

  if(d->comm.rank==0) {
    lclimgebot=blwentry(d->lclcblk,0,e_size);
    clear_vec(lclimgebot, e_size);
  }

  if(d->comm.rank==d->comm.size-1) {
    lclimgetop=abventry(d->lclcblk,d->layers_per_node+1,e_size);
    clear_vec(lclimgetop, e_size);
  }

  /// calculate local cellblks. partblks don't make sense
  for (c = 1; c <= d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c-1];
    d->lclcblk[size*c] = 0;
    for (i = 0; i < np; i++) {
      ci=i+offset;
      d->lclcblk[size*c] += d->local_charges[ci];

      if (d->dielectric_contrast_on) {
        e_di_l = (d->delta_mult*d->delta_mid_bot
        + d->delta_mult)*fac_imgsum;
        if (c==1 && d->comm.rank==0) {
          e = d->delta_mid_bot;
           lclimgebot[MMM2D_QQEQQP] += d->local_charges[ci]*e;
        }
        else
          e_di_l += d->delta_mid_bot;

        e_di_h = (d->delta_mult*d->delta_mid_top
          + d->delta_mult)*fac_imgsum;

        if (c==d->layers_per_node && d->comm.rank==d->comm.size-1) {
          e = d->delta_mid_top;
          lclimgetop[MMM2D_QQEQQP] += d->local_charges[ci]*e;
        }
        else
          e_di_h += d->delta_mid_top;

        d->lclimge[MMM2D_QQEQQP] += d->local_charges[ci]*e_di_l;
        d->lclimge[MMM2D_QQEQQM] += d->local_charges[ci]*e_di_h;
      }
    }
    d->lclcblk[size*c] *= pref;
    d->lclcblk[size*c+1] = d->lclcblk[size*c];

    offset += np;
  }

  if (d->dielectric_contrast_on) {
    scale_vec(pref, d->lclimge, size);
    if(d->comm.rank==0)
      scale_vec(pref, blwentry(d->lclcblk, 0, e_size), e_size);
    if(d->comm.rank==d->comm.size-1)
      scale_vec(pref, abventry(d->lclcblk, d->layers_per_node + 1, e_size), e_size);
  }
}

static void add_z_force(mmm2d_data_struct *d, fcs_float *forces)
{
  fcs_int c, i, np, offset=0;
  fcs_float add;
  fcs_float *othcblk;
  fcs_int size = 2;

  for (c = 0; c < d->layers_per_node; c++) {
    othcblk = block(d->gblcblk, c, size);
    add = othcblk[MMM2D_QQEQQP] - othcblk[MMM2D_QQEQQM];
    np   = d->zslices_nparticles[c];
    //printf("add_z_force: add %d -> %f\n", c, add);
    for (i = 0; i < np; i++) {
      forces[3*(i+offset)+2] +=d->local_charges[i]*add;
      /*LOG_FORCES(fprintf(stderr, "%d: part %d force %10.3g %10.3g %10.3g\n",
          this_node, part[i].p.identity, part[i].f.f[0],
          part[i].f.f[1], part[i].f.f[2]));*/
    }
    offset+=np;
  }
}

static void setup_z_energy(mmm2d_data_struct *d)
{
  fcs_int np, c, i, ci, offset=0;
  fcs_float pref = -MMM_COMMON_C_2PI*d->ux*d->uy;
  fcs_int e_size = 2, size = 4;

  if (d->comm.rank == 0)
    /// the lowest lclcblk does not contain anything, since there are no charges below the simulation box, at least for this term.
    clear_vec(blwentry(d->lclcblk, 0, e_size), e_size);

  if (d->comm.rank == d->comm.size - 1)
    /// same for the top node
    clear_vec(abventry(d->lclcblk, d->layers_per_node + 1, e_size), e_size);

  /// calculate local cellblks. partblks don't make sense
  for (c = 1; c <= d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c-1];
    clear_vec(blwentry(d->lclcblk, c, e_size), e_size);
    for (i = 0; i < np; i++) {
      ci=i+offset;
      d->lclcblk[size*c + MMM2D_ABEQQP] += d->local_charges[ci];
      d->lclcblk[size*c + MMM2D_ABEQZP] += d->local_charges[ci]*d->local_positions[3*ci+2];
    }
    scale_vec(pref, blwentry(d->lclcblk, c, e_size), e_size);
    /// just to be able to use the standard distribution. Here below and above terms are the same
    copy_vec(abventry(d->lclcblk, c, e_size), blwentry(d->lclcblk, c, e_size), e_size);
    offset += np;
  }
}

static fcs_float z_energy(mmm2d_data_struct *d)
{
  fcs_float eng = 0.;
  fcs_int np, c, i, ci, offset=0;
  fcs_float *othcblk;
  fcs_int size = 4;
  for (c = 1; c <= d->layers_per_node; c++) {
    othcblk = block(d->gblcblk, c - 1, size);
    np   = d->zslices_nparticles[c-1];
    for (i = 0; i < np; i++) {
      ci=i+offset;
      eng += d->local_charges[ci]*(d->local_positions[3*ci+2]*othcblk[MMM2D_ABEQQP] - othcblk[MMM2D_ABEQZP] -
           d->local_positions[3*ci+2]*othcblk[MMM2D_ABEQQM] + othcblk[MMM2D_ABEQZM]);
    }
    offset += np;
  }
  return eng;
}

/* PoQ exp sum */
/*****************************************************************/
static void setup_P(mmm2d_data_struct *d, fcs_int p, fcs_float omega, fcs_float fac)
{
  fcs_int np, c, i, ic, ci, posid, offset=0, o = (p-1)*d->n_localpart;
  fcs_float pref = 4*M_PI*d->ux*d->uy*fac*fac;
  fcs_float h = d->box_l[2];
  fcs_float fac_imgsum = 1/(1 - d->delta_mult*exp(-omega*2*h));
  fcs_float fac_delta_mid_bot = d->delta_mid_bot*fac_imgsum; 
  fcs_float fac_delta_mid_top = d->delta_mid_top*fac_imgsum;
  fcs_float fac_delta         = d->delta_mult*fac_imgsum;
  fcs_float layer_top;
  fcs_float e, e_di_l, e_di_h;
  fcs_float *llclcblk;
  fcs_float *lclimgebot = NULL, *lclimgetop = NULL;
  fcs_int e_size = 2, size = 4;

  if (d->dielectric_contrast_on)
    clear_vec(d->lclimge, size);

  if(d->comm.rank==0) {
    /* on the lowest node, clear the lclcblk below, which only contains the images of the lowest layer
       if there is dielectric contrast, otherwise it is empty */
    lclimgebot = block(d->lclcblk, 0, size);
    clear_vec(blwentry(d->lclcblk, 0, e_size), e_size);
  }
  if(d->comm.rank==d->comm.size-1) {
    /* same for the top node */
    lclimgetop = block(d->lclcblk, d->layers_per_node + 1, size);
    clear_vec(abventry(d->lclcblk, d->layers_per_node + 1, e_size), e_size);
  }

  layer_top = d->my_bottom + d->layer_h;
  ic = 0;
  
  for (c = 1; c <= d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c-1];
    llclcblk = block(d->lclcblk, c, size);

    clear_vec(llclcblk, size);
    
    for (i = 0; i < np; i++) {
      ci=i+offset;
      posid=3*ci+2;
      e = exp(omega*(d->local_positions[posid] - layer_top));

      d->partblk[size*ic + MMM2D_POQESM] = d->local_charges[ci]*d->scxcache[o + ic].s/e;
      d->partblk[size*ic + MMM2D_POQESP] = d->local_charges[ci]*d->scxcache[o + ic].s*e;
      d->partblk[size*ic + MMM2D_POQECM] = d->local_charges[ci]*d->scxcache[o + ic].c/e;
      d->partblk[size*ic + MMM2D_POQECP] = d->local_charges[ci]*d->scxcache[o + ic].c*e;

      /* take images due to different dielectric constants into account */
      if (d->dielectric_contrast_on) {
        if (c==1 && d->comm.rank==0) {
        /* There are image charges at -(2h+z) and -(2h-z) etc. layer_h included due to the shift
        in z */
          e_di_l = ( exp(omega*(-(d->local_positions[posid]) - 2*h + d->layer_h))*(d->delta_mid_bot) +
           exp(omega*( d->local_positions[posid] - 2*h + d->layer_h))                   )*fac_delta;

          e = exp(omega*(-(d->local_positions[posid])))*d->delta_mid_bot;

          lclimgebot[MMM2D_POQESP] += d->local_charges[ci]*d->scxcache[o + ic].s*e;
          lclimgebot[MMM2D_POQECP] += d->local_charges[ci]*d->scxcache[o + ic].c*e;
        }
        else
        /* There are image charges at -(z) and -(2h-z) etc. layer_h included due to the shift in z */
          e_di_l = ( exp(omega*(-(d->local_positions[posid]) + d->layer_h)) +
           exp(omega*( d->local_positions[posid] - 2*h + d->layer_h))*d->delta_mid_top )*fac_delta_mid_bot;

        if (c==d->layers_per_node && d->comm.rank==d->comm.size-1) {
        /* There are image charges at (3h-z) and (h+z) from the top layer etc. layer_h included
           due to the shift in z */
          e_di_h = (exp(omega*( d->local_positions[posid] - 3*h + 2*d->layer_h))*d->delta_mid_top +
          exp(omega*(-(d->local_positions[posid]) - h + 2*d->layer_h)))*fac_delta;

        /* There are image charges at (h-z) layer_h included due to the shift in z */
          e = exp(omega*(d->local_positions[posid] - h + d->layer_h))*d->delta_mid_top;

          lclimgetop[MMM2D_POQESM]+= d->local_charges[ci]*d->scxcache[o + ic].s*e;
          lclimgetop[MMM2D_POQECM]+= d->local_charges[ci]*d->scxcache[o + ic].c*e;
        }
        else
          /* There are image charges at (h-z) and (h+z) from the top layer etc. layer_h included
          due to the shift in z */
            e_di_h = (exp(omega*( d->local_positions[posid] - h + 2*d->layer_h)) +
            exp(omega*(-(d->local_positions[posid]) - h + 2*d->layer_h))*d->delta_mid_bot )*fac_delta_mid_top;

          d->lclimge[MMM2D_POQESP] += d->local_charges[ci]*d->scxcache[o + ic].s*e_di_l;
          d->lclimge[MMM2D_POQECP] += d->local_charges[ci]*d->scxcache[o + ic].c*e_di_l;
          d->lclimge[MMM2D_POQESM] += d->local_charges[ci]*d->scxcache[o + ic].s*e_di_h;
          d->lclimge[MMM2D_POQECM] += d->local_charges[ci]*d->scxcache[o + ic].c*e_di_h;
        }

        add_vec(llclcblk, llclcblk, block(d->partblk, ic, size), size);
        ic++;
      }
      scale_vec(pref, blwentry(d->lclcblk, c, e_size), e_size);
      scale_vec(pref, abventry(d->lclcblk, c, e_size), e_size);

      layer_top += d->layer_h;

      offset += np;
    }

  if (d->dielectric_contrast_on) {
    scale_vec(pref, d->lclimge, size);
    if(d->comm.rank==0)
      scale_vec(pref, blwentry(d->lclcblk, 0, e_size), e_size);
    if(d->comm.rank==d->comm.size-1)
      scale_vec(pref, abventry(d->lclcblk, d->layers_per_node + 1, e_size), e_size);
  }
}

/* compare setup_P */
static void setup_Q(mmm2d_data_struct *d, fcs_int q, fcs_float omega, fcs_float fac)
{

  fcs_int np, c, i, ic, ci, posid, offset=0, o = (q-1)*d->n_localpart;
  fcs_float pref = 4*M_PI*d->ux*d->uy*fac*fac;
  fcs_float h = d->box_l[2];
  fcs_float fac_imgsum = 1/(1 - d->delta_mult*exp(-omega*2*h));
  fcs_float fac_delta_mid_bot = d->delta_mid_bot*fac_imgsum; 
  fcs_float fac_delta_mid_top = d->delta_mid_top*fac_imgsum;
  fcs_float fac_delta         = d->delta_mult*fac_imgsum;
  fcs_float layer_top;
  fcs_float e, e_di_l, e_di_h;
  fcs_float *llclcblk;
  fcs_float *lclimgebot=NULL, *lclimgetop=NULL;
  fcs_int e_size = 2, size = 4;

  if (d->dielectric_contrast_on)
    clear_vec(d->lclimge, size); 

  if(d->comm.rank==0) {
    lclimgebot = block(d->lclcblk, 0, size);
    clear_vec(blwentry(d->lclcblk, 0, e_size), e_size);
  }

  if(d->comm.rank==d->comm.size-1) {
    lclimgetop = block(d->lclcblk, d->layers_per_node + 1, size);
    clear_vec(abventry(d->lclcblk, d->layers_per_node + 1, e_size), e_size);
  }

  layer_top = d->my_bottom + d->layer_h;
  ic = 0;
  for (c = 1; c <= d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c-1];
    llclcblk = block(d->lclcblk, c, size);

    clear_vec(llclcblk, size);

    for (i = 0; i < np; i++) {
      ci=i+offset;
      posid=3*ci+2;
      e = exp(omega*(d->local_positions[posid] - layer_top));

      d->partblk[size*ic + MMM2D_POQESM] = d->local_charges[ci]*d->scycache[o + ic].s/e;
      d->partblk[size*ic + MMM2D_POQESP] = d->local_charges[ci]*d->scycache[o + ic].s*e;
      d->partblk[size*ic + MMM2D_POQECM] = d->local_charges[ci]*d->scycache[o + ic].c/e;
      d->partblk[size*ic + MMM2D_POQECP] = d->local_charges[ci]*d->scycache[o + ic].c*e;

      if (d->dielectric_contrast_on) {
        if(c==1 && d->comm.rank==0) {
          e_di_l = (exp(omega*(-(d->local_positions[posid]) -2*h + d->layer_h))*d->delta_mid_bot +
           exp(omega*(d->local_positions[posid] - 2*h + d->layer_h)))*fac_delta;

          e = exp(omega*(-(d->local_positions[posid])))*d->delta_mid_bot;

          lclimgebot[MMM2D_POQESP] += d->local_charges[ci]*d->scycache[o + ic].s*e;
          lclimgebot[MMM2D_POQECP] += d->local_charges[ci]*d->scycache[o + ic].c*e;
        }
        else
          e_di_l = ( exp(omega*(-(d->local_positions[posid]) + d->layer_h)) +
           exp(omega*( d->local_positions[posid] - 2*h + d->layer_h))*d->delta_mid_top )*fac_delta_mid_bot;

        if(c==d->layers_per_node && d->comm.rank==d->comm.size-1) {
          e_di_h = (exp(omega*( d->local_positions[posid] -3*h + 2*d->layer_h))*d->delta_mid_top +
          exp(omega*(-(d->local_positions[posid]) - h + 2*d->layer_h))                 )*fac_delta;

          e = exp(omega*(d->local_positions[posid] - h + d->layer_h))*d->delta_mid_top;

          lclimgetop[MMM2D_POQESM] += d->local_charges[ci]*d->scycache[o + ic].s*e;
          lclimgetop[MMM2D_POQECM] += d->local_charges[ci]*d->scycache[o + ic].c*e;
        }
        else
          e_di_h = ( exp(omega*( d->local_positions[posid] - h + 2*d->layer_h)) +
           exp(omega*(-(d->local_positions[posid]) - h + 2*d->layer_h))*d->delta_mid_bot )*fac_delta_mid_top;

        d->lclimge[MMM2D_POQESP] += d->local_charges[ci]*d->scycache[o + ic].s*e_di_l;
        d->lclimge[MMM2D_POQECP] += d->local_charges[ci]*d->scycache[o + ic].c*e_di_l;
        d->lclimge[MMM2D_POQESM] += d->local_charges[ci]*d->scycache[o + ic].s*e_di_h;
        d->lclimge[MMM2D_POQECM] += d->local_charges[ci]*d->scycache[o + ic].c*e_di_h;
      }

      add_vec(llclcblk, llclcblk, block(d->partblk, ic, size), size);
      ic++;
    }
    scale_vec(pref, blwentry(d->lclcblk, c, e_size), e_size);
    scale_vec(pref, abventry(d->lclcblk, c, e_size), e_size);

    layer_top += d->layer_h;

    offset += np;
  }

  if (d->dielectric_contrast_on) {
    scale_vec(pref, d->lclimge, size);
    if(d->comm.rank==0)
      scale_vec(pref, blwentry(d->lclcblk, 0, e_size), e_size);
    if(d->comm.rank==d->comm.size-1)
      scale_vec(pref, abventry(d->lclcblk, d->layers_per_node + 1, e_size), e_size);
  }
}

static void add_P_force(mmm2d_data_struct *d, fcs_float *forces)
{
  //printf("add_P_force\n");
  fcs_int np, c, i, ic, ix, iz, offset=0;
  fcs_float *othcblk;
  fcs_int size = 4;

  ic = 0;
  for (c = 0; c < d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c];
    othcblk = block(d->gblcblk, c, size);
    for (i = 0; i < np; i++) {
      ix=3*(i+offset);
      iz=ix+2;
      forces[ix] +=
   d->partblk[size*ic + MMM2D_POQESM]*othcblk[MMM2D_POQECP] - d->partblk[size*ic + MMM2D_POQECM]*othcblk[MMM2D_POQESP] +
   d->partblk[size*ic + MMM2D_POQESP]*othcblk[MMM2D_POQECM] - d->partblk[size*ic + MMM2D_POQECP]*othcblk[MMM2D_POQESM];
      forces[iz] +=
   d->partblk[size*ic + MMM2D_POQECM]*othcblk[MMM2D_POQECP] + d->partblk[size*ic + MMM2D_POQESM]*othcblk[MMM2D_POQESP] -
   d->partblk[size*ic + MMM2D_POQECP]*othcblk[MMM2D_POQECM] - d->partblk[size*ic + MMM2D_POQESP]*othcblk[MMM2D_POQESM];

      /*LOG_FORCES(fprintf(stderr, "%d: part %d force %10.3g %10.3g %10.3g\n",
          this_node, part[i].p.identity, part[i].f.f[0],
          part[i].f.f[1], part[i].f.f[2]));*/
      ic++;
    }
    offset+=np;
  }
}

static fcs_float P_energy(mmm2d_data_struct *d, fcs_float omega)
{
  fcs_float eng = 0.;

  fcs_int np, c, i, ic;
  fcs_float *othcblk;
  fcs_int size = 4;
  fcs_float pref = 1/omega;

  ic = 0;
  for (c = 1; c <= d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c-1];
    othcblk = block(d->gblcblk, c - 1, size);
    for (i = 0; i < np; i++) {
      eng += pref*(d->partblk[size*ic + MMM2D_POQECM]*othcblk[MMM2D_POQECP] + d->partblk[size*ic + MMM2D_POQESM]*othcblk[MMM2D_POQESP] +
         d->partblk[size*ic + MMM2D_POQECP]*othcblk[MMM2D_POQECM] + d->partblk[size*ic + MMM2D_POQESP]*othcblk[MMM2D_POQESM]);
      ic++;
    }
  }
  return eng;
}

static void add_Q_force(mmm2d_data_struct *d, fcs_float *forces)
{
  //printf("add_Q_force\n");
  fcs_int np, c, i, ic, iy, iz, offset=0;
  fcs_float *othcblk;
  fcs_int size = 4;

  ic = 0;
  for (c = 0; c < d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c];
    othcblk = block(d->gblcblk, c, size);

    for (i = 0; i < np; i++) {
      iy=3*(i+offset)+1;
      iz=iy+1;
      forces[iy] +=
   d->partblk[size*ic + MMM2D_POQESM]*othcblk[MMM2D_POQECP] - d->partblk[size*ic + MMM2D_POQECM]*othcblk[MMM2D_POQESP] +
   d->partblk[size*ic + MMM2D_POQESP]*othcblk[MMM2D_POQECM] - d->partblk[size*ic + MMM2D_POQECP]*othcblk[MMM2D_POQESM];
      forces[iz] +=
   d->partblk[size*ic + MMM2D_POQECM]*othcblk[MMM2D_POQECP] + d->partblk[size*ic + MMM2D_POQESM]*othcblk[MMM2D_POQESP] -
   d->partblk[size*ic + MMM2D_POQECP]*othcblk[MMM2D_POQECM] - d->partblk[size*ic + MMM2D_POQESP]*othcblk[MMM2D_POQESM];
      /*LOG_FORCES(fprintf(stderr, "%d: part %d force %10.3g %10.3g %10.3g\n",
          this_node, part[i].p.identity, part[i].f.f[0],
          part[i].f.f[1], part[i].f.f[2]));*/
      ic++;
    }
    offset+=np;
  }
}

static fcs_float Q_energy(mmm2d_data_struct *d, fcs_float omega)
{
  fcs_float eng = 0.;

  fcs_int np, c, i, ic;
  fcs_float *othcblk;
  fcs_int size = 4;
  fcs_float pref = 1/omega;

  ic = 0;
  for (c = 1; c <= d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c-1];
    othcblk = block(d->gblcblk, c - 1, size);
    for (i = 0; i < np; i++) {
      eng += pref*(d->partblk[size*ic + MMM2D_POQECM]*othcblk[MMM2D_POQECP] + d->partblk[size*ic + MMM2D_POQESM]*othcblk[MMM2D_POQESP] +
         d->partblk[size*ic + MMM2D_POQECP]*othcblk[MMM2D_POQECM] + d->partblk[size*ic + MMM2D_POQESP]*othcblk[MMM2D_POQESM]);
      ic++;
    }
  }

  return eng;
}

/*****************************************************************/
/* PQ particle blocks */
/*****************************************************************/

/* compare setup_P */
static void setup_PQ(mmm2d_data_struct *d, fcs_int p, fcs_int q, fcs_float omega, fcs_float fac)
{
  fcs_int np, c, i, ic, ci, posid, offset=0, ox = (p - 1)*d->n_localpart, oy = (q - 1)*d->n_localpart;
  fcs_float pref = 8*M_PI*d->ux*d->uy*fac*fac;
  fcs_float h = d->box_l[2];
  fcs_float fac_imgsum = 1/(1 - d->delta_mult*exp(-omega*2*h));
  fcs_float fac_delta_mid_bot = d->delta_mid_bot*fac_imgsum; 
  fcs_float fac_delta_mid_top = d->delta_mid_top*fac_imgsum;
  fcs_float fac_delta         = d->delta_mult*fac_imgsum;
  fcs_float layer_top;
  fcs_float e, e_di_l, e_di_h;
  fcs_float *llclcblk;
  fcs_float *lclimgebot=NULL, *lclimgetop=NULL;
  fcs_int e_size = 4, size = 8;

  if (d->dielectric_contrast_on)
    clear_vec(d->lclimge, size); 

  if(d->comm.rank==0) {
    lclimgebot = block(d->lclcblk, 0, size);
    clear_vec(blwentry(d->lclcblk, 0, e_size), e_size);
  }
  
  if(d->comm.rank==d->comm.size-1) {
    lclimgetop = block(d->lclcblk, d->layers_per_node + 1, size);
    clear_vec(abventry(d->lclcblk, d->layers_per_node + 1, e_size), e_size);
  }
  
  layer_top = d->my_bottom + d->layer_h;
  ic = 0;
  for (c = 1; c <= d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c-1];
    
    llclcblk = block(d->lclcblk, c, size);

    clear_vec(llclcblk, size);

    for (i = 0; i < np; i++) {
      ci=i+offset;
      posid=3*ci+2;
      e = exp(omega*(d->local_positions[posid] - layer_top));
      
      d->partblk[size*ic + MMM2D_PQESSM] = d->scxcache[ox + ic].s*d->scycache[oy + ic].s*(d->local_charges[ci])/e;
      d->partblk[size*ic + MMM2D_PQESCM] = d->scxcache[ox + ic].s*d->scycache[oy + ic].c*(d->local_charges[ci])/e;
      d->partblk[size*ic + MMM2D_PQECSM] = d->scxcache[ox + ic].c*d->scycache[oy + ic].s*(d->local_charges[ci])/e;
      d->partblk[size*ic + MMM2D_PQECCM] = d->scxcache[ox + ic].c*d->scycache[oy + ic].c*(d->local_charges[ci])/e;

      d->partblk[size*ic + MMM2D_PQESSP] = d->scxcache[ox + ic].s*d->scycache[oy + ic].s*(d->local_charges[ci])*e;
      d->partblk[size*ic + MMM2D_PQESCP] = d->scxcache[ox + ic].s*d->scycache[oy + ic].c*(d->local_charges[ci])*e;
      d->partblk[size*ic + MMM2D_PQECSP] = d->scxcache[ox + ic].c*d->scycache[oy + ic].s*(d->local_charges[ci])*e;
      d->partblk[size*ic + MMM2D_PQECCP] = d->scxcache[ox + ic].c*d->scycache[oy + ic].c*(d->local_charges[ci])*e;

      if (d->dielectric_contrast_on) {
        if(c==1 && d->comm.rank==0) {
          e_di_l = (exp(omega*(-(d->local_positions[posid])- 2*h + d->layer_h))*d->delta_mid_bot +
          exp(omega*( d->local_positions[posid] - 2*h + d->layer_h)))*fac_delta;

          e = exp(omega*(-(d->local_positions[posid])))*d->delta_mid_bot;

          lclimgebot[MMM2D_PQESSP] += d->scxcache[ox + ic].s*d->scycache[oy + ic].s*(d->local_charges[ci])*e;
          lclimgebot[MMM2D_PQESCP] += d->scxcache[ox + ic].s*d->scycache[oy + ic].c*(d->local_charges[ci])*e;
          lclimgebot[MMM2D_PQECSP] += d->scxcache[ox + ic].c*d->scycache[oy + ic].s*(d->local_charges[ci])*e;
          lclimgebot[MMM2D_PQECCP] += d->scxcache[ox + ic].c*d->scycache[oy + ic].c*(d->local_charges[ci])*e;
        }
        else
          e_di_l = ( exp(omega*(-(d->local_positions[posid]) + d->layer_h)) +
           exp(omega*( d->local_positions[posid] - 2*h + d->layer_h))*d->delta_mid_top )*fac_delta_mid_bot;

        if(c==d->layers_per_node && d->comm.rank==d->comm.size-1) {
          e_di_h = (exp(omega*( d->local_positions[posid]- 3*h + 2*d->layer_h))*d->delta_mid_top +
           exp(omega*(-(d->local_positions[posid]) - h + 2*d->layer_h)) )*fac_delta;

          e = exp(omega*(d->local_positions[posid]-h+d->layer_h))*d->delta_mid_top;

          lclimgetop[MMM2D_PQESSM] += d->scxcache[ox + ic].s*d->scycache[oy + ic].s*(d->local_charges[ci])*e;
          lclimgetop[MMM2D_PQESCM] += d->scxcache[ox + ic].s*d->scycache[oy + ic].c*(d->local_charges[ci])*e;
          lclimgetop[MMM2D_PQECSM] += d->scxcache[ox + ic].c*d->scycache[oy + ic].s*(d->local_charges[ci])*e;
          lclimgetop[MMM2D_PQECCM] += d->scxcache[ox + ic].c*d->scycache[oy + ic].c*(d->local_charges[ci])*e;
        }
        else
          e_di_h = ( exp(omega*( d->local_positions[posid] - h + 2*d->layer_h)) +
           exp(omega*(-(d->local_positions[posid]) - h + 2*d->layer_h))*d->delta_mid_bot)*fac_delta_mid_top;

        d->lclimge[MMM2D_PQESSP] += d->scxcache[ox + ic].s*d->scycache[oy + ic].s*(d->local_charges[ci])*e_di_l;
        d->lclimge[MMM2D_PQESCP] += d->scxcache[ox + ic].s*d->scycache[oy + ic].c*(d->local_charges[ci])*e_di_l;
        d->lclimge[MMM2D_PQECSP] += d->scxcache[ox + ic].c*d->scycache[oy + ic].s*(d->local_charges[ci])*e_di_l;
        d->lclimge[MMM2D_PQECCP] += d->scxcache[ox + ic].c*d->scycache[oy + ic].c*(d->local_charges[ci])*e_di_l;

        d->lclimge[MMM2D_PQESSM] += d->scxcache[ox + ic].s*d->scycache[oy + ic].s*(d->local_charges[ci])*e_di_h;
        d->lclimge[MMM2D_PQESCM] += d->scxcache[ox + ic].s*d->scycache[oy + ic].c*(d->local_charges[ci])*e_di_h;
        d->lclimge[MMM2D_PQECSM] += d->scxcache[ox + ic].c*d->scycache[oy + ic].s*(d->local_charges[ci])*e_di_h;
        d->lclimge[MMM2D_PQECCM] += d->scxcache[ox + ic].c*d->scycache[oy + ic].c*(d->local_charges[ci])*e_di_h;
      }
      
      add_vec(llclcblk, llclcblk, block(d->partblk, ic, size), size);
      ic++;
    }
    offset+=i;
    
    scale_vec(pref, blwentry(d->lclcblk, c, e_size), e_size);
    scale_vec(pref, abventry(d->lclcblk, c, e_size), e_size);
    
    layer_top += d->layer_h;
  }

  if (d->dielectric_contrast_on) {
    scale_vec(pref, d->lclimge, size);

    if(d->comm.rank==0)
      scale_vec(pref, blwentry(d->lclcblk, 0, e_size), e_size);
    if(d->comm.rank==d->comm.size-1)
      scale_vec(pref, abventry(d->lclcblk, d->layers_per_node + 1, e_size), e_size);
  }
}

static void add_PQ_force(mmm2d_data_struct *d, fcs_int p, fcs_int q, fcs_float omega, fcs_float *forces)
{
  //printf("add_PQ_force\n");
  fcs_int np, c, i, ic, ix, iy, iz, offset=0;
  fcs_float pref_x = MMM_COMMON_C_2PI*d->ux*p/omega;
  fcs_float pref_y = MMM_COMMON_C_2PI*d->uy*q/omega;
  fcs_float *othcblk;
  fcs_int size = 8;

  ic = 0;
  for (c = 0; c < d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c];
    othcblk = block(d->gblcblk, c, size);

    for (i = 0; i < np; i++) {
      ix=3*(i+offset);
      iy=ix+1;
      iz=iy+1;
      forces[ix] +=
   pref_x*(d->partblk[size*ic + MMM2D_PQESCM]*othcblk[MMM2D_PQECCP] + d->partblk[size*ic + MMM2D_PQESSM]*othcblk[MMM2D_PQECSP] -
      d->partblk[size*ic + MMM2D_PQECCM]*othcblk[MMM2D_PQESCP] - d->partblk[size*ic + MMM2D_PQECSM]*othcblk[MMM2D_PQESSP] +
      d->partblk[size*ic + MMM2D_PQESCP]*othcblk[MMM2D_PQECCM] + d->partblk[size*ic + MMM2D_PQESSP]*othcblk[MMM2D_PQECSM] -
      d->partblk[size*ic + MMM2D_PQECCP]*othcblk[MMM2D_PQESCM] - d->partblk[size*ic + MMM2D_PQECSP]*othcblk[MMM2D_PQESSM]);
      forces[iy] +=
   pref_y*(d->partblk[size*ic + MMM2D_PQECSM]*othcblk[MMM2D_PQECCP] + d->partblk[size*ic + MMM2D_PQESSM]*othcblk[MMM2D_PQESCP] -
      d->partblk[size*ic + MMM2D_PQECCM]*othcblk[MMM2D_PQECSP] - d->partblk[size*ic + MMM2D_PQESCM]*othcblk[MMM2D_PQESSP] +
      d->partblk[size*ic + MMM2D_PQECSP]*othcblk[MMM2D_PQECCM] + d->partblk[size*ic + MMM2D_PQESSP]*othcblk[MMM2D_PQESCM] -
      d->partblk[size*ic + MMM2D_PQECCP]*othcblk[MMM2D_PQECSM] - d->partblk[size*ic + MMM2D_PQESCP]*othcblk[MMM2D_PQESSM]);
      forces[iz] +=
          (d->partblk[size*ic + MMM2D_PQECCM]*othcblk[MMM2D_PQECCP] + d->partblk[size*ic + MMM2D_PQECSM]*othcblk[MMM2D_PQECSP] +
           d->partblk[size*ic + MMM2D_PQESCM]*othcblk[MMM2D_PQESCP] + d->partblk[size*ic + MMM2D_PQESSM]*othcblk[MMM2D_PQESSP] -
           d->partblk[size*ic + MMM2D_PQECCP]*othcblk[MMM2D_PQECCM] - d->partblk[size*ic + MMM2D_PQECSP]*othcblk[MMM2D_PQECSM] -
           d->partblk[size*ic + MMM2D_PQESCP]*othcblk[MMM2D_PQESCM] - d->partblk[size*ic + MMM2D_PQESSP]*othcblk[MMM2D_PQESSM]);

      /*LOG_FORCES(fprintf(stderr, "%d: part %d force %10.3g %10.3g %10.3g\n",
          this_node, part[i].p.identity, part[i].f.f[0],
          part[i].f.f[1], part[i].f.f[2]));*/
      ic++;
    }
    offset+=np;
  }
}

static fcs_float PQ_energy(mmm2d_data_struct *d, fcs_float omega)
{
  fcs_float eng = 0;
  fcs_int np, c, i, ic;
  fcs_float *othcblk;
  fcs_int size = 8;
  fcs_float pref = 1/omega;

  ic = 0;
  for (c = 1; c <= d->layers_per_node; c++) {
    np   = d->zslices_nparticles[c-1];
    othcblk = block(d->gblcblk, c - 1, size);

    for (i = 0; i < np; i++) {
      eng += pref*(d->partblk[size*ic + MMM2D_PQECCM]*othcblk[MMM2D_PQECCP] + d->partblk[size*ic + MMM2D_PQECSM]*othcblk[MMM2D_PQECSP] +
         d->partblk[size*ic + MMM2D_PQESCM]*othcblk[MMM2D_PQESCP] + d->partblk[size*ic + MMM2D_PQESSM]*othcblk[MMM2D_PQESSP] +
         d->partblk[size*ic + MMM2D_PQECCP]*othcblk[MMM2D_PQECCM] + d->partblk[size*ic + MMM2D_PQECSP]*othcblk[MMM2D_PQECSM] +
         d->partblk[size*ic + MMM2D_PQESCP]*othcblk[MMM2D_PQESCM] + d->partblk[size*ic + MMM2D_PQESSP]*othcblk[MMM2D_PQESSM]);
      ic++;
    }
  }

  return eng;
}

static fcs_float dielectric_layers_energy_contribution(mmm2d_data_struct *d)
{
  ///@TODO: revisar
  fcs_int c, ci, cj, i, j;
  fcs_int      npl;
  fcs_float dist2, dv[3];
  fcs_float charge_factor;
  fcs_float a[3], b[3];
  fcs_float eng=0.0;
  // prefactor for the charged plate interaction removal
  fcs_float corr_pref = MMM_COMMON_C_2PI*d->ux*d->uy;

  if (!d->dielectric_contrast_on) return 0.0;
  if(d->comm.rank==0) {
    npl = d->zslices_nparticles[0];
    for(i = 0; i < npl; i++) {
      c=3*i;
      a[0]=d->local_positions[c]; a[1]=d->local_positions[c+1]; a[2]=d->local_positions[c+2];
      for(j = 0; j < npl; j++) {
        c=3*j;
        b[0]=d->local_positions[c]; b[1]=d->local_positions[c+1]; b[2]=-d->local_positions[c+2];
        layered_displacement_vector(d, a[0], a[1], a[2], b[0], b[1], b[2], dv);
        dist2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
        charge_factor = d->delta_mid_bot*d->local_charges[i]*d->local_charges[j];
        eng+=mmm2d_pair_energy_near(d, charge_factor, dv, sqrt(dist2)) + corr_pref*charge_factor*dv[2];
      }
    }
  }
  
  if(d->comm.rank==d->comm.size-1) {
    npl = d->zslices_nparticles[d->layers_per_node-1];
    c=d->n_localpart-npl;
    for(i = 0; i < npl; i++) {
      ci=3*(c+i);
      a[0]=d->local_positions[ci]; a[1]=d->local_positions[ci+1]; a[2]=d->local_positions[ci+2];
      for(j = 0; j < npl; j++) {
        cj=3*(c+j);
        b[0]=d->local_positions[cj]; b[1]=d->local_positions[cj+1]; b[2]=d->box_l[2]*2.0-d->local_positions[cj+2];
        layered_displacement_vector(d, a[0], a[1], a[2], b[0], b[1], b[2], dv);
        dist2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
        charge_factor=d->delta_mid_top*d->local_charges[c+i]*d->local_charges[c+j];
        eng+=mmm2d_pair_energy_near(d, charge_factor, dv, sqrt(dist2)) - corr_pref*charge_factor*dv[2];
      }
    }
  }
  
  return 0.5*eng;
}

static void dielectric_layers_force_contribution(mmm2d_data_struct *d, fcs_float *forces)
{
  ///@TODO: revisar
  fcs_int c, ci, cj, i, j;
  fcs_int      npl;
  fcs_float dist2, dv[3];
  fcs_float charge_factor;
  fcs_float a[3], b[3];
  fcs_float force[3]={0, 0, 0};
  
  if (!d->dielectric_contrast_on) return;
  
  if(d->comm.rank==0) {
    npl = d->zslices_nparticles[0];
    
    for(i = 0; i < npl; i++) {
      ci=3*i;
      a[0]=d->local_positions[ci]; a[1]=d->local_positions[ci+1]; a[2]=d->local_positions[ci+2];
      for(j = 0; j < npl; j++) {
        cj=3*j;
        b[0]=d->local_positions[cj]; b[1]=d->local_positions[cj+1]; b[2]=-d->local_positions[cj+2];
        layered_displacement_vector(d, a[0], a[1], a[2], b[0], b[1], b[2], dv);
        dist2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
        charge_factor=d->local_charges[j]*d->delta_mid_bot;
        mmm2d_pair_force_near(d, charge_factor, dv, dist2, sqrt(dist2), force);
        forces[ci]+=force[0];
        forces[ci+1]+=force[1];
        forces[ci+2]+=force[2];
      }
    }
  }

  if(d->comm.rank==d->comm.size-1) {
    npl = d->zslices_nparticles[d->layers_per_node-1];
    c=d->n_localpart-npl;
    for(i = 0; i < npl; i++) {
      ci=3*(c+i);
      a[0]=d->local_positions[ci]; a[1]=d->local_positions[ci+1]; a[2]=d->local_positions[ci+2];
      for(j = 0; j < npl; j++) {
        cj=3*(c+j);
        b[0]=d->local_positions[cj]; b[1]=d->local_positions[cj+1]; b[2]=d->box_l[2]*2.0-d->local_positions[cj+2];
        layered_displacement_vector(d, a[0], a[1], a[2], b[0], b[1], b[2], dv);
        dist2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
        charge_factor=d->local_charges[j]*d->delta_mid_top;
        mmm2d_pair_force_near(d, charge_factor, dv, dist2, sqrt(dist2), force);
        forces[ci]+=force[0];
        forces[ci+1]+=force[1];
        forces[ci+2]+=force[2];
      }
    }
  }
}

/* the data transfer routine for the lclcblks itself */
static void distribute(mmm2d_data_struct *d, fcs_int e_size, fcs_float fac)
{
  //printf("distribute\n");
  fcs_int c, node, inv_node;
  fcs_float sendbuf[8];
  fcs_float recvbuf[8];
  MPI_Status status;

  /* send/recv to/from other nodes. Also builds up the gblcblk. */
  for (node = 0; node < d->comm.size; node++) {
    inv_node = d->comm.size - node - 1;
    /* up */
    if (node == d->comm.rank) {
      /* calculate sums of cells below */
      for (c = 1; c < d->layers_per_node; c++)
   addscale_vec(blwentry(d->gblcblk, c, e_size), fac, blwentry(d->gblcblk, c - 1, e_size), blwentry(d->lclcblk, c - 1, e_size), e_size);

      /* calculate my ghost contribution only if a node above exists */
      if (node + 1 < d->comm.size) {
   addscale_vec(sendbuf, fac, blwentry(d->gblcblk, d->layers_per_node - 1, e_size), blwentry(d->lclcblk, d->layers_per_node - 1, e_size), e_size);
   copy_vec(sendbuf + e_size, blwentry(d->lclcblk, d->layers_per_node, e_size), e_size);
   //printf("rank %d, node %d, Send 1!!!!!!! %d, %d\n", d->comm.rank, node, e_size, sizeof(sendbuf));
   MPI_Send(sendbuf, 2*e_size, FCS_MPI_FLOAT, node + 1, 0, d->comm.mpicomm);
      }
    }
    else if (node + 1 == d->comm.rank) {
       //printf("rank %d, node %d, Recv 1!!!!!!! %d, %d\n", d->comm.rank, node, e_size, sizeof(recvbuf));
      MPI_Recv(recvbuf, 2*e_size, FCS_MPI_FLOAT, node, 0, d->comm.mpicomm, &status);
      copy_vec(blwentry(d->gblcblk, 0, e_size), recvbuf, e_size);
      copy_vec(blwentry(d->lclcblk, 0, e_size), recvbuf + e_size, e_size);
    }

    /* down */
    
    if (inv_node == d->comm.rank) {
      /* calculate sums of all cells above */
      for (c = d->layers_per_node + 1; c > 2; c--)
   addscale_vec(abventry(d->gblcblk, c - 3, e_size), fac, abventry(d->gblcblk, c - 2, e_size), abventry(d->lclcblk, c, e_size), e_size);
      
      /* calculate my ghost contribution only if a node below exists */
      if (inv_node -  1 >= 0) {
   addscale_vec(sendbuf, fac, abventry(d->gblcblk, 0, e_size), abventry(d->lclcblk, 2, e_size), e_size);
   copy_vec(sendbuf + e_size, abventry(d->lclcblk, 1, e_size), e_size);
   //printf("rank %d, node %d, Send 2!!!!!!! %d, %d\n", d->comm.rank, node, e_size, sizeof(sendbuf));
   MPI_Send(sendbuf, 2*e_size, FCS_MPI_FLOAT, inv_node - 1, 0, d->comm.mpicomm);
      }
    }
    else if (inv_node - 1 == d->comm.rank) {
       //printf("rank %d, node %d, Recv 2!!!!!!! %d, %d\n", d->comm.rank, node, e_size, sizeof(recvbuf));
      MPI_Recv(recvbuf, 2*e_size, FCS_MPI_FLOAT, inv_node, 0, d->comm.mpicomm, &status);
      copy_vec(abventry(d->gblcblk, d->layers_per_node - 1, e_size), recvbuf, e_size);
      copy_vec(abventry(d->lclcblk, d->layers_per_node + 1, e_size), recvbuf + e_size, e_size);
    }
  }
}

/* dealing with the image contributions from far outside the simulation box */
static void gather_image_contributions(mmm2d_data_struct *d, fcs_int e_size)
{
   //printf("gather image contribution\n");
  fcs_float recvbuf[8];

  //* collect the image charge contributions with at least a layer distance 
  MPI_Allreduce(d->lclimge, recvbuf, 2*e_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (d->comm.rank == 0)
    /* the gblcblk contains all contributions from layers deeper than one layer below our system,
       which is precisely what the gblcblk should contain for the lowest layer. */
    copy_vec(blwentry(d->gblcblk, 0, e_size),recvbuf, e_size);

  if (d->comm.rank == d->comm.size - 1)
    //* same for the top node
    copy_vec(abventry(d->gblcblk, d->layers_per_node - 1, e_size), recvbuf + e_size, e_size);
}

static void clear_image_contributions(mmm2d_data_struct *d, fcs_int e_size)
{
  if (d->comm.rank == 0)
    /* the gblcblk contains all contributions from layers deeper than one layer below our system,
       which is precisely what the gblcblk should contain for the lowest layer. */
    clear_vec(blwentry(d->gblcblk, 0, e_size), e_size);

  if (d->comm.rank == d->comm.size - 1)
    /* same for the top node */
    clear_vec(abventry(d->gblcblk, d->layers_per_node - 1, e_size), e_size);
}

/*****************************************************************/
/* data distribution */
/*****************************************************************/

/* vector operations */

/** pdc = 0 */
static void clear_vec(fcs_float *pdc, fcs_int size)
{
  fcs_int i;
  for (i = 0; i < size; i++)
    pdc[i] = 0;
}

/** pdc_d = pdc_s */
static void copy_vec(fcs_float *pdc_d, fcs_float *pdc_s, fcs_int size)
{
  fcs_int i;
  for (i = 0; i < size; i++)
    pdc_d[i] = pdc_s[i];
}

/** pdc_d = pdc_s1 + pdc_s2 */
static void add_vec(fcs_float *pdc_d, fcs_float *pdc_s1, fcs_float *pdc_s2, fcs_int size)
{
  fcs_int i;
  for (i = 0; i < size; i++)
    pdc_d[i] = pdc_s1[i] + pdc_s2[i];
}

/** pdc_d = scale*pdc_s1 + pdc_s2 */
static void addscale_vec(fcs_float *pdc_d, fcs_float scale, fcs_float *pdc_s1, fcs_float *pdc_s2, fcs_int size)
{
  fcs_int i;
  for (i = 0; i < size; i++)
    pdc_d[i] = scale*pdc_s1[i] + pdc_s2[i];
}

/** pdc_d = scale*pdc */
static void scale_vec(fcs_float scale, fcs_float *pdc, fcs_int size)
{
  fcs_int i;
  for (i = 0; i < size; i++)
    pdc[i] *= scale;
}

/* block indexing - has to fit to the PQ block definitions above.
   size gives the full size of the data block,
   e_size is the size of only the top or bottom half, i.e. half of size.
*/
static fcs_float *block(fcs_float *p, fcs_int index, fcs_int size)
{
  return &p[index*size];
}

static fcs_float *blwentry(fcs_float *p, fcs_int index, fcs_int e_size)
{
  return &p[2*index*e_size];
}

static fcs_float *abventry(fcs_float *p, fcs_int index, fcs_int e_size)
{
  return &p[(2*index + 1)*e_size];
}

/**
Ugly function to get the displacement vector between two folded positions
*/
void layered_displacement_vector(mmm2d_data_struct *d, fcs_float x1, fcs_float y1, fcs_float z1, fcs_float x2, fcs_float y2, fcs_float z2, fcs_float disp[3]) {
//printf("disp para %f %f %f <-> %f %f %f = %f %f %f\n", x1, y1, z1, x2, y2, z2, disp[0], disp[1], disp[2]);
  disp[0] = x1 - x2;
  disp[0] -= floor(disp[0]*d->box_l_i[0] + 0.5)*d->box_l[0];
  disp[1] = y1 - y2;
  disp[1] -= floor(disp[1]*d->box_l_i[1] + 0.5)*d->box_l[1];

  disp[2] = z1 - z2;
  //printf("disp para %f %f %f <-> %f %f %f = %f %f %f\n", x1, y1, z1, x2, y2, z2, disp[0], disp[1], disp[2]);
}
