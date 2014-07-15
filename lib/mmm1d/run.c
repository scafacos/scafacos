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
#include "tune.h"
#include "common/gridsort/gridsort.h"
#include "FCSCommon.h"
#include <stdlib.h>
// #include <stdio.h>
#include <math.h>

static fcs_float mmm1d_coulomb_pair_energy(mmm1d_data_struct *d, fcs_float disp[3]);
static void mmm1d_coulomb_pair_force(mmm1d_data_struct *d, fcs_float disp[3], fcs_float force[3]);

void mmm1d_run(void* rd,
        fcs_int num_particles,
        fcs_int max_num_particles,
        fcs_float *positions,
        fcs_float *charges,
        fcs_float *fields,
        fcs_float *potentials) {
  /* Here we assume, that the method is tuned and that all parameters are valid */
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  
  /* decompose system */
  fcs_int local_num_particles=0;
  fcs_int local_num_real_particles=0;
  fcs_int *local_num_ghost_particles;
  fcs_float *local_positions;
  fcs_float *local_charges;
  fcs_gridsort_index_t *local_indexes;
  fcs_float box_base[3] = {0.0, 0.0, 0.0 };
  fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
  fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
  fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };
  fcs_gridsort_t gridsort;
  fcs_int comm_size;
  MPI_Comm_size(d->comm, &comm_size);
  fcs_int comm_rank;
  MPI_Comm_rank(d->comm, &comm_rank);

  fcs_gridsort_create(&gridsort);
  fcs_gridsort_set_system(&gridsort, box_base, box_a, box_b, box_c, NULL);
  fcs_gridsort_set_particles(&gridsort, num_particles, max_num_particles, positions, charges);
  
  fcs_gridsort_sort_forward(&gridsort, 0.0, d->comm);
  
  fcs_gridsort_separate_ghosts(&gridsort, &local_num_real_particles, NULL);
  fcs_gridsort_get_sorted_particles(&gridsort, &local_num_real_particles, NULL, NULL, NULL, NULL);
  
  fcs_gridsort_get_real_particles(&gridsort, &local_num_real_particles, &local_positions, &local_charges, &local_indexes);
  
  // fprintf(stderr, "rank %d, local particles %d\n", comm_rank, local_num_real_particles);
  
  /* assign workload to nodes and distribute ghosts */
  // fprintf(stderr, "rank %d, assign workload to nodes\n", comm_rank);
  fcs_int n_ghost_neighbors=0, n_clairvoyant_neighbors=0, n_ghosts=0;
  fcs_int ghost_neighbors[(comm_size + 3)/2], clairvoyant_neighbors[(comm_size + 3)/2];
  fcs_float **local_ghosts_positions = (fcs_float **)malloc((comm_size-1)*sizeof(fcs_float *));
  fcs_float **local_ghosts_charges = (fcs_float **)malloc((comm_size-1)*sizeof(fcs_float *));
  local_num_ghost_particles = (fcs_int *)malloc(sizeof(fcs_int)*(comm_size-1));
  MPI_Status status;
  for(fcs_int i=0; i< comm_size; i++) {
    fcs_int p1 = i - comm_rank;
    if (p1==0) continue;
    if ((p1 > 0 && p1 % 2 == 0) || (p1 < 0 && -p1 % 2 == 1)) {
      ghost_neighbors[n_ghost_neighbors]=i;
      //printf("rank %d :receive workload from %d\n", comm_rank, i);
      MPI_Recv(&n_ghosts, 1, FCS_MPI_INT, i, i, d->comm, &status);
      local_num_ghost_particles[n_ghost_neighbors] = n_ghosts;
      local_ghosts_positions[n_ghost_neighbors] = (fcs_float *)malloc(sizeof(fcs_float)*3*n_ghosts);
      local_ghosts_charges[n_ghost_neighbors] = (fcs_float *)malloc(sizeof(fcs_float)*n_ghosts);
      MPI_Recv(local_ghosts_positions[n_ghost_neighbors], 3*n_ghosts, FCS_MPI_FLOAT, i, i, d->comm, &status);
      MPI_Recv(local_ghosts_charges[n_ghost_neighbors], n_ghosts, FCS_MPI_FLOAT, i, i, d->comm, &status);
      n_ghost_neighbors++;
    } else {
      clairvoyant_neighbors[n_clairvoyant_neighbors]=i;
      n_clairvoyant_neighbors++;
      MPI_Send(&local_num_real_particles, 1, FCS_MPI_INT, i, comm_rank, d->comm);
      MPI_Send(local_positions, 3*local_num_real_particles, FCS_MPI_FLOAT, i, comm_rank, d->comm);
      MPI_Send(local_charges, local_num_real_particles, FCS_MPI_FLOAT, i, comm_rank, d->comm);
      //printf("rank %d :send workload to %d\n", comm_rank, i);
    }
  }

  MPI_Barrier(d->comm);
  
  /* allocate local containers */
  // fprintf(stderr, "rank %d, allocate containers\n", comm_rank);
  fcs_float *local_fields = NULL;
  fcs_float *local_potentials = NULL;
  local_fields=malloc(sizeof(fcs_float)*3*local_num_real_particles);
  local_potentials=malloc(sizeof(fcs_float)*local_num_real_particles);
  for(fcs_int p1=0; p1<local_num_real_particles; p1++) {
    fcs_int c1=3*p1;
    local_fields[c1]=0;
    local_fields[c1+1]=0;
    local_fields[c1+2]=0;
    local_potentials[p1]=0;
  }
  fcs_float **local_ghosts_fields = (fcs_float **)malloc(n_ghost_neighbors*sizeof(fcs_float *));
  fcs_float **local_ghosts_potentials = (fcs_float **)malloc(n_ghost_neighbors*sizeof(fcs_float *));
  for(fcs_int i=0; i<n_ghost_neighbors; i++) {
    local_ghosts_fields[i]=(fcs_float *)malloc(3*local_num_ghost_particles[i]*sizeof(fcs_float));
    local_ghosts_potentials[i]=(fcs_float *)malloc(local_num_ghost_particles[i]*sizeof(fcs_float));
    for(fcs_int p1=0; p1<local_num_ghost_particles[i]; p1++) {
      fcs_int c1=3*p1;
      local_ghosts_fields[i][c1]=0;
      local_ghosts_fields[i][c1+1]=0;
      local_ghosts_fields[i][c1+2]=0;
      local_ghosts_potentials[i][p1]=0;
    }
  }
  
  /* calculate local interactions */
  // fprintf(stderr, "rank %d, local interactions\n", comm_rank);
  for(fcs_int p1=0; p1<local_num_real_particles-1; p1++) {
    fcs_int c1=3*p1;
    fcs_float x1=local_positions[c1];
    fcs_float y1=local_positions[c1+1];
    fcs_float z1=local_positions[c1+2];
  
    for(fcs_int p2=p1+1; p2<local_num_real_particles; p2++) {
      fcs_int c2=3*p2;
      fcs_float x2=local_positions[c2];
      fcs_float y2=local_positions[c2+1];
      fcs_float z2=local_positions[c2+2];
      
      fcs_float disp[3];
      mmm_distance2vec(x1,y1,z1,x2,y2,z2,disp);
      // printf("mmm1d_run, selected real-real particles: %d, %d\n",p1,p2);

      if (potentials) {
        fcs_float eng = mmm1d_coulomb_pair_energy(d, disp);

        local_potentials[p1]+=local_charges[p2]*eng;
        local_potentials[p2]+=local_charges[p1]*eng;
      }

      if (fields) {
        fcs_float field[3];
        mmm1d_coulomb_pair_force(d, disp, field);

        local_fields[c1]  += local_charges[p2]*field[0];
        local_fields[c1+1]+= local_charges[p2]*field[1];
        local_fields[c1+2]+= local_charges[p2]*field[2];
        local_fields[c2]  +=-local_charges[p1]*field[0];
        local_fields[c2+1]+=-local_charges[p1]*field[1];
        local_fields[c2+2]+=-local_charges[p1]*field[2];
      }
    }
  }
  
  /* calculate ghost interactions */
  // fprintf(stderr, "rank %d, ghost interactions\n", comm_rank);
  for(fcs_int p1=0; p1<local_num_real_particles; p1++) {
    fcs_int c1=3*p1;
    fcs_float x1=local_positions[c1];
    fcs_float y1=local_positions[c1+1];
    fcs_float z1=local_positions[c1+2];
    for(fcs_int i=0; i<n_ghost_neighbors; i++){
      for(fcs_int p2=0; p2<local_num_ghost_particles[i]; p2++) {
        fcs_int c2=3*p2;
        fcs_float x2=local_ghosts_positions[i][c2];
        fcs_float y2=local_ghosts_positions[i][c2+1];
        fcs_float z2=local_ghosts_positions[i][c2+2];


        fcs_float disp[3];      
        mmm_distance2vec(x1,y1,z1,x2,y2,z2,disp);
        //printf("mmm1d_run, selected real-ghost particles: %d, %d\n",p1,p2);

        if (potentials) {
          fcs_float eng = mmm1d_coulomb_pair_energy(d, disp);

          local_potentials[p1]          +=local_ghosts_charges[i][p2]*eng;
          local_ghosts_potentials[i][p2]+=local_charges[p1]*eng;
        }
        if (fields) {
          fcs_float field[3];
          mmm1d_coulomb_pair_force(d, disp, field);

          local_fields[c1]  +=local_ghosts_charges[i][p2]*field[0];
          local_fields[c1+1]+=local_ghosts_charges[i][p2]*field[1];
          local_fields[c1+2]+=local_ghosts_charges[i][p2]*field[2];
          local_ghosts_fields[i][c2]  +=-local_charges[p1]*field[0];
          local_ghosts_fields[i][c2+1]+=-local_charges[p1]*field[1];
          local_ghosts_fields[i][c2+2]+=-local_charges[p1]*field[2];
        }
      }
    }
  }
  
  //printf("rank %d, locals: potential: %e, field 0: %e %e %e\n", comm_rank, local_potentials[0], local_fields[0], local_fields[1], local_fields[2]);
  //if (comm_rank==1) printf("rank %d, local ghosts: potential: %e, field 0: %e %e %e\n", comm_rank, local_ghosts_potentials[0][0], local_ghosts_fields[0][0], local_ghosts_fields[0][1], local_ghosts_fields[0][2]);
  
  /* send back ghost interactions to neighbors */
  // fprintf(stderr, "rank %d, send back ghosts\n", comm_rank);
  fcs_float foreign_ghost_fields[3*local_num_real_particles], foreign_ghost_potentials[local_num_real_particles];
  
  fcs_int c1=0; fcs_int c2=0;
  for(fcs_int i=0; i< comm_size; i++) {
    fcs_int p1 = i - comm_rank;
    if (p1==0) continue;
    if ((p1 > 0 && p1 % 2 == 0) || (p1 < 0 && -p1 % 2 == 1)) {
      //send:
      if(local_num_ghost_particles[c1]>0) {
        MPI_Send(local_ghosts_fields[c1], 3*local_num_ghost_particles[c1], FCS_MPI_FLOAT, ghost_neighbors[c1], comm_rank, d->comm);
        MPI_Send(local_ghosts_potentials[c1], local_num_ghost_particles[c1], FCS_MPI_FLOAT, ghost_neighbors[c1], comm_rank, d->comm);
      }
      c1++;
    } else {
      //receive
      if (local_num_real_particles>0) {
        MPI_Recv(foreign_ghost_fields, 3*local_num_real_particles, FCS_MPI_FLOAT, clairvoyant_neighbors[c2], clairvoyant_neighbors[c2], d->comm, &status);
        MPI_Recv(foreign_ghost_potentials, local_num_real_particles, FCS_MPI_FLOAT, clairvoyant_neighbors[c2], clairvoyant_neighbors[c2], d->comm, &status);
        for(fcs_int p1=0; p1<local_num_real_particles; p1++) {
          fcs_int p2=3*p1;
          local_fields[p2]+=foreign_ghost_fields[p2];
          local_fields[p2+1]+=foreign_ghost_fields[p2+1];
          local_fields[p2+2]+=foreign_ghost_fields[p2+2];
          local_potentials[p1]+=foreign_ghost_potentials[p2];
        }
      }
      c2++;
    }
  }
  
  //printf("rank %d, local_potential 0: %e\n", comm_rank, local_potentials[0]);
  
  /* sort back, clean up and finish */
  fcs_gridsort_sort_backward(&gridsort,
              local_fields, local_potentials,
              fields, potentials, 1,
              d->comm);
  
  fcs_gridsort_free(&gridsort);
  
  fcs_gridsort_destroy(&gridsort);
}

fcs_float mmm1d_coulomb_pair_energy(mmm1d_data_struct *d, fcs_float disp[3])
{
  fcs_float rxy2, rxy2_d, z_d, E;
  
  rxy2   = disp[0]*disp[0] + disp[1]*disp[1];
  rxy2_d = rxy2*d->uz2;
  z_d    = disp[2]*d->uz;
  
  if (rxy2 <= d->far_switch_radius_2) {
    //* near range formula *
    fcs_float r2n, rt, shift_z, add;
    fcs_int n;

    E = -2*MMM_COMMON_C_GAMMA;

    //* polygamma summation *
    r2n = 1.0;
    for (n = 0; n < (d->polTaylor)->n_modPsi; n++) {
      add = mmm_mod_psi_even(d->polTaylor, n, z_d)*r2n;
      E -= add;
      
      if (fabs(add) < d->maxPWerror) break;
      
      r2n *= rxy2_d;
    }
    E *= d->uz;

    //* real space parts *

    E += 1./sqrt(rxy2 + disp[2]*disp[2]); 

    shift_z = disp[2] + d->box_l[2];
    E += 1./sqrt(rxy2 + shift_z*shift_z);

    shift_z = disp[2] - d->box_l[2];
    E += 1./sqrt(rxy2 + shift_z*shift_z);
  }
  else {
    //* far range formula *
    fcs_float rxy   = sqrt(rxy2);
    fcs_float rxy_d = rxy*d->uz;
    fcs_int bp;
    // The first Bessel term will compensate a little bit the log term, so add them close together
    E = -0.25*log(rxy2_d) + 0.5*(M_LN2 - MMM_COMMON_C_GAMMA);
    for (bp = 1; bp < d->bessel_cutoff; bp++) {
      if (d->bessel_radii[bp-1] < rxy) {
        break;
      }
      fcs_float fq = MMM_COMMON_C_2PI*bp;
      E += mmm_K0(fq*rxy_d)*cos(fq*z_d);
    }
    E *= 4.*d->uz;
  }
  
  return E;
}

void mmm1d_coulomb_pair_force(mmm1d_data_struct *d, fcs_float disp[3], fcs_float force[3])
{
  fcs_int dim;
  fcs_float rxy2, rxy2_d, z_d;
  fcs_float pref;
  fcs_float Fx, Fy, Fz;
  
  rxy2   = disp[0]*disp[0] + disp[1]*disp[1];
  rxy2_d = rxy2*d->uz2;
  z_d    = disp[2]*d->uz;
  
  //printf("rank %d, rxy2: %e, rxy2_d: %e, z_d: %e, far_switch_radius_2: %e\n", comm_rank, rxy2, rxy2_d, z_d, d->far_switch_radius_2);
  
  if (rxy2 <= d->far_switch_radius_2) {
    // near range formula
    fcs_float sr, sz, r2nm1, rt, rt2, shift_z;
    fcs_int n;

    // polygamma summation
    sr = 0;
    sz = mmm_mod_psi_odd(d->polTaylor, 0, z_d);
/*printf("rank %d, sz: %e\n", comm_rank, sz);
for(n=0; n<d->polTaylor->n_modPsi; n++) {
   for(dim=0; dim<d->polTaylor->modPsi[n].n; dim++) {
      printf("rank %d, taylor [%d][%d]: %e\n", comm_rank, n, dim, d->polTaylor->modPsi[n].e[dim]);
  }
}*/
    r2nm1 = 1.0;
    for (n = 1; n < (d->polTaylor)->n_modPsi; n++) {
      fcs_float deriv = 2*n;
      fcs_float mpe   = mmm_mod_psi_even(d->polTaylor, n, z_d);
      fcs_float mpo   = mmm_mod_psi_odd(d->polTaylor, n, z_d);
      fcs_float r2n   = r2nm1*rxy2_d;

      sz +=         r2n*mpo;
      sr += deriv*r2nm1*mpe;

      if (fabs(deriv*r2nm1*mpe) < d->maxPWerror) break;

      r2nm1 = r2n;
    }

    Fx = d->L3_i*sr*disp[0];
    Fy = d->L3_i*sr*disp[1];
    Fz = d->uz2*sz;
    //printf("rank %d, F: %e %e %e\n", comm_rank, Fx, Fx, Fz);

    // real space parts

    shift_z = disp[2];
    rt2 = rxy2 + shift_z*shift_z;
    pref = 1./(rt2*sqrt(rt2));
    Fx += pref*disp[0];
    Fy += pref*disp[1];
    Fz += pref*disp[2];

    shift_z = disp[2] + d->box_l[2];
    rt2 = rxy2 + shift_z*shift_z;
    pref = 1./(rt2*sqrt(rt2));
    Fx += pref*disp[0];
    Fy += pref*disp[1];
    Fz += pref*shift_z;

    shift_z = disp[2] - d->box_l[2];
    rt2 = rxy2 + shift_z*shift_z;
    rt  = sqrt(rt2);
    pref = 1./(rt2*rt); 
    Fx += pref*disp[0];
    Fy += pref*disp[1];
    Fz += pref*shift_z;

    force[0] = Fx;
    force[1] = Fy;
    force[2] = Fz;
  }
  else {
    // far range formula
    fcs_float rxy   = sqrt(rxy2);
    fcs_float rxy_d = rxy*d->uz;
    fcs_float sr = 0, sz = 0;
    fcs_int bp;

    for (bp = 1; bp < d->bessel_cutoff; bp++) {
      if (d->bessel_radii[bp-1] < rxy)
        break;

      fcs_float fq = MMM_COMMON_C_2PI*bp, k0, k1;
      
#ifdef MMM_BESSEL_MACHINE_PREC
      k0 = mmm_K0(fq*rxy_d);
      k1 = mmm_K1(fq*rxy_d);
#else
      mmm_LPK01(fq*rxy_d, &k0, &k1);
#endif
      sr += bp*k1*cos(fq*z_d);
      sz += bp*k0*sin(fq*z_d);
    }
    sr *= d->uz2*4*MMM_COMMON_C_2PI;
    sz *= d->uz2*4*MMM_COMMON_C_2PI;
    
    pref = 1.*(sr/rxy + 2*d->uz/rxy2);

    force[0] = pref*disp[0];
    force[1] = pref*disp[1];
    force[2] = sz;
  }

  //printf("rank: %d, disp: %e %e %e, charge1: %e, charge2: %e, force: %e %e %e\n", comm_rank, disp[0], disp[1], disp[2], charge1, charge2,  chpref * F[0], chpref * F[1], chpref * F[2]);
}
