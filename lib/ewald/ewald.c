/*
  Copyright (C) 2011,2012 Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.
  
  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/* Implementation of the Ewald summation

   This implementation uses the same notation and variable names as
   Deserno, Holm. "How to mesh up Ewald sums". J. Chem. Phys. 109(18), 1998.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FCSResult.h"
#include "FCSInterface.h"
#include "FCSCommon.h"
#include "common/near/near.h"
#include "ewald.h"


/***************************************************
 **** TUNING
 ***************************************************/

/** Computes the real space contribution to the rms error in the
   force (as described by Kolafa and Perram) as well as its
   derivative.
   \param N the number of charged particles in the system.
   \param sum_q2 the sum of square of charges in the system
   \param box_l fcs_float[3] system size
   \param r_cut the real-space cutoff
   \param alpha the Ewald splitting parameter
   \param error (out) the real space error
   \param derivative (out) the derivative of the real space error
*/
static inline void
ewald_real_space_error(fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
		       fcs_float r_cut, fcs_float alpha,
		       fcs_float *error, fcs_float *derivative) {
  const fcs_float V = box_l[0]*box_l[1]*box_l[2];
  *error = 2.0*sum_q2 / sqrt(N*r_cut*V) * exp(-SQR(alpha*r_cut));
  *derivative = -2.0 * alpha * SQR(r_cut) * *error;
}

/** Computes the reciprocal space contribution to the rms error in the
   force (as described by Kolafa and Perram) as well as its derivative.
   \param N the number of charged particles in the system.
   \param sum_q2 the sum of square of charges in the system
   \param box_l fcs_float[3] system size
   \param kmax the maximal k vector
   \param alpha the Ewald splitting parameter
   \param error (out) the recipocal space error
   \param derivative (out) the derivative of the recipocal space error
*/
static inline void
ewald_k_space_error(fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
		    fcs_int kmax, fcs_float alpha,
		    fcs_float *error, fcs_float *derivative) {
  fcs_float Lmax = box_l[0];
  if (box_l[1] > Lmax) Lmax = box_l[1];
  if (box_l[2] > Lmax) Lmax = box_l[2];
  const fcs_float K = 2.0*M_PI*kmax/Lmax;
  const fcs_float V = box_l[0]*box_l[1]*box_l[2];
  const fcs_float fak1 = 2.0*sum_q2 / sqrt(N*M_PI*K*V) * exp(-K*K/(4.0*SQR(alpha)));
  *error = fak1 * alpha;
  *derivative = fak1 * (2.0*SQR(alpha)+SQR(kmax)/(2.0*SQR(alpha)));
}

void ewald_tune_alpha(fcs_int N, fcs_float sum_q2, fcs_float box_l[3],
                 fcs_float r_cut, fcs_int kmax, 
                 fcs_float *alpha, fcs_float *error, int tune) {
  fcs_float err_r, der_r;
  fcs_float err_k, der_k;
  fcs_float alpha_diff;

  if (tune) {
    /* Newton-Raphson method to find optimal alpha */
    *alpha = 0.1;
    do {
      /* get errors and derivatives */
      ewald_real_space_error(N, sum_q2, box_l, r_cut, *alpha, &err_r, &der_r);
      ewald_k_space_error(N, sum_q2, box_l, kmax, *alpha, &err_k, &der_k);
      
      alpha_diff = (err_r - err_k) / (der_r - der_k);
      *alpha -= alpha_diff;
    } while (fabs(alpha_diff) > ALPHA_OPT_PREC);
  } else {
    ewald_real_space_error(N, sum_q2, box_l, r_cut, *alpha, &err_r, &der_r);
    ewald_k_space_error(N, sum_q2, box_l, kmax, *alpha, &err_k, &der_k);
  }
  *error = sqrt(SQR(err_r) + SQR(err_k));
}


/***************************************************
 **** K-SPACE CONTRIBUTION
 ***************************************************/
void ewald_compute_kspace(ewald_data_struct* d, 
			  fcs_int num_particles,
			  fcs_float *positions, 
			  fcs_float *charges,
			  fcs_float *fields,
			  fcs_float *potentials) {

  FCS_INFO(fprintf(stderr, "ewald_compute_kspace started...\n"));

  /* DISTRIBUTE ALL PARTICLE DATA TO ALL NODES */
  /* Gather all particle numbers */
  int node_num_particles = num_particles;
  int node_particles[d->comm_size]; 
  int node_particles3[d->comm_size]; 
  int displs[d->comm_size];
  int displs3[d->comm_size];
  fcs_int total_particles;

  /* printf("%d: num_particles=%d\n", d->comm_rank, num_particles); */
  MPI_Allgather(&node_num_particles, 1, MPI_INT, node_particles, 1, MPI_INT, d->comm);

  /* compute displacements for MPI_Gatherv */
  total_particles = node_particles[0];
  node_particles3[0] = node_particles[0]*3;
  displs[0] = 0;
  displs3[0] = 0;
  for (fcs_int i=1; i < d->comm_size; i++) {
    total_particles += node_particles[i];
    node_particles3[i] = node_particles[i]*3;
    displs[i] = displs[i-1] + node_particles[i-1];
    displs3[i] = displs3[i-1] + node_particles3[i-1];
  }

  fcs_float *all_positions, *all_charges, *node_fields, *node_potentials;

  all_positions = malloc(sizeof(fcs_float) * 3 * total_particles);
  all_charges = malloc(sizeof(fcs_float) * total_particles);
  node_fields = malloc(sizeof(fcs_float) * 3 * total_particles);
  node_potentials = malloc(sizeof(fcs_float) * total_particles);

  /* gather all particle data at all nodes */
  MPI_Allgatherv(positions, num_particles*3, FCS_MPI_FLOAT, all_positions, node_particles3, displs3, FCS_MPI_FLOAT, d->comm);
  MPI_Allgatherv(charges, num_particles, FCS_MPI_FLOAT, all_charges, node_particles, displs, FCS_MPI_FLOAT, d->comm);

  /* for (fcs_int i=0; i < total_particles; i++) { */
  /*   printf("%d: all_positions[%d]={%lf, %lf, %lf}\n", d->comm_rank, i, all_positions[i*3], all_positions[3*i+1], all_positions[3*i+2]); */
  /*   printf("%d: all_charges[%d]=%lf\n", d->comm_rank, i, all_charges[i]); */
  /* } */

  /* INIT ALGORITHM */

  /* init fields and potentials */
  if (fields != NULL) {
    for (fcs_int i=0; i < total_particles; i++) {
      node_fields[3*i] = 0.0;
      node_fields[3*i+1] = 0.0;
      node_fields[3*i+2] = 0.0;
    }
  }
  if (potentials != NULL)
    for (fcs_int i=0; i < total_particles; i++)
      node_potentials[i] = 0.0;
  
  /* COMPUTE FAR FIELDS */

  /* evenly distribute the k-vectors onto all tasks */
  fcs_int num_k_per_dir = 2*d->kmax+1;
  fcs_int num_k = num_k_per_dir * num_k_per_dir * num_k_per_dir;
  for (int k_ind = d->comm_rank; k_ind < num_k; k_ind += d->comm_size) {
    /* compute fields and potentials */
    fcs_int nx = 
      k_ind % num_k_per_dir - d->kmax;
    fcs_int ny = 
      k_ind % (num_k_per_dir*num_k_per_dir) / num_k_per_dir - d->kmax;
    fcs_int nz = 
      k_ind / (num_k_per_dir*num_k_per_dir) - d->kmax;
    if (nx*nx + ny*ny + nz*nz <= d->kmax*d->kmax) {
      // system length vector L_vec
      const fcs_float Lx = d->box_l[0];
      const fcs_float Ly = d->box_l[1];
      const fcs_float Lz = d->box_l[2];
      // reciprocal vector k_vec
      const fcs_float kx = 2.0*M_PI*nx / Lx;
      const fcs_float ky = 2.0*M_PI*ny / Ly;
      const fcs_float kz = 2.0*M_PI*nz / Lz;
      /* reciprocal charge density rhohat */
      fcs_float rhohat_re = 0.0;
      fcs_float rhohat_im = 0.0;
      
      /* compute Deserno, Holm (1998) eq. (8) */
      for (fcs_int i=0; i < total_particles; i++) {
	// charge q
	const fcs_float q = all_charges[i];
	if (!fcs_float_is_zero(q)) {
	  // particle position r_vec
	  const fcs_float rx = all_positions[3*i];
	  const fcs_float ry = all_positions[3*i+1];
	  const fcs_float rz = all_positions[3*i+2];
	  // compute k_vec*r_vec
	  fcs_float kr = kx*rx + ky*ry + kz*rz;
	  // rhohat = qi * exp(-i*k_vec*r_vec)
	  rhohat_re += q * cos(kr);
	  rhohat_im += q * -sin(kr);
	}
      }
      
      /* FCS_DEBUG(fprintf(stderr, "  n_vec= (%d, %d, %d) rhohat_re=%e rhohat_im=%e\n", \ */
      /* 		    nx, ny, nz, rhohat_re, rhohat_im)); */
      
      /* fetch influence function */
      fcs_float g = d->G[linindex(abs(nx), abs(ny), abs(nz), d->kmax)];
      for (fcs_int i=0; i < total_particles; i++) {
	// particle position r_vec
	const fcs_float rx = all_positions[3*i];
	const fcs_float ry = all_positions[3*i+1];
	const fcs_float rz = all_positions[3*i+2];
	// compute k_vec*r_vec
	fcs_float kr = kx*rx + ky*ry + kz*rz;
	
	if (fields != NULL) {
	  /* compute field at position of particle i
	     compare to Deserno, Holm (1998) eq. (15) */
	  fcs_float fak1 = g * (rhohat_re*sin(kr) + rhohat_im*cos(kr));
	  node_fields[3*i] += kx * fak1;
	  node_fields[3*i+1] += ky * fak1;
	  node_fields[3*i+2] += kz * fak1;
	}
	if (potentials != NULL) {
	  /* compute potential at position of particle i 
	     compare to Deserno, Holm (1998) eq. (9) */
	  node_potentials[i] += g * (rhohat_re*cos(kr) - rhohat_im*sin(kr));
	}
      }
    }
  }
  
  /* printf("%d: node_fields[0]=%lf\n", d->comm_rank, node_fields[0]); */

  /* REDISTRIBUTE COMPUTED FAR FIELDS AND POTENTIALS  */
  if (fields != NULL) {
    fcs_float all_fields[total_particles*3];
    for (fcs_int pid=0; pid < total_particles; pid++) {
      all_fields[3*pid] = 0.0;
      all_fields[3*pid+1] = 0.0;
      all_fields[3*pid+2] = 0.0;
    }

    /* Combine all fields on master */
    MPI_Reduce(node_fields, all_fields, total_particles*3, 
	       FCS_MPI_FLOAT, MPI_SUM, 0, d->comm);
    /* if (d->comm_rank == 0) */
    /*   printf("%d: all_fields[0]=%lf\n", d->comm_rank, all_fields[0]); */
    /* Scatter the fields to the task that holds the particle */
    MPI_Scatterv(all_fields, node_particles3, displs3, FCS_MPI_FLOAT,
		 fields, num_particles*3, FCS_MPI_FLOAT, 0, d->comm);
  }

  if (potentials != NULL) {
    fcs_float all_potentials[total_particles];
    for (fcs_int pid=0; pid < total_particles; pid++)
      all_potentials[pid] = 0.0;
    /* Combine all potentials on master */
    MPI_Reduce(node_potentials, all_potentials, total_particles, 
	       FCS_MPI_FLOAT, MPI_SUM, 0, d->comm);
    /* Scatter the potentials to the task that holds the particle */
    MPI_Scatterv(all_potentials, node_particles, displs, FCS_MPI_FLOAT,
		 potentials, num_particles, FCS_MPI_FLOAT, 0, d->comm);

    /* subtract self potential */
    FCS_INFO(fprintf(stderr, "  subtracting self potential...\n"));
    for (fcs_int i=0; i < num_particles; i++) {
      potentials[i] -= charges[i] * M_2_SQRTPI * d->alpha;
    }
  }

  free(all_positions);
  free(all_charges);
  free(node_fields);
  free(node_potentials);

  /* now each task should have its far field components */
  FCS_INFO(fprintf(stderr, "ewald_compute_kspace finished.\n"));
}

/***************************************************
 **** REAL SPACE CONTRIBUTION
 ***************************************************/
/* callback function for near field computations */
static inline void 
ewald_compute_near(const void *param, fcs_float dist, fcs_float *f, fcs_float *p)
{
  fcs_float alpha = *((fcs_float *) param);
  fcs_float adist = alpha * dist;

#if FCS_EWALD_USE_ERFC_APPROXIMATION

  /* approximate \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying a formula from:
     Abramowitz/Stegun: Handbook of Mathematical Functions, Dover
     (9. ed.), chapter 7 */
  fcs_float t = 1.0 / (1.0 + 0.3275911 * adist);
  fcs_float erfc_part_ri = exp(-adist*dist) * 
    (t * (0.254829592 + 
	  t * (-0.284496736 + 
	       t * (1.421413741 + 
		    t * (-1.453152027 + 
			 t * 1.061405429))))) 
    / dist;

  *p = erfc_part_ri;
  *f = -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498) 
    / dist;
  
#else

  fcs_float erfc_part_ri = erfc(adist) / dist;
  *p = erfc_part_ri;
  *f = -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498*exp(-adist*adist)) 
    / dist;

#endif
}

/* callback function for performing a whole loop of near field computations (using ewald_compute_near) */
FCS_NEAR_LOOP_FP(ewald_compute_near_loop, ewald_compute_near)

void ewald_compute_rspace(ewald_data_struct* d, 
			  fcs_int num_particles,
			  fcs_int max_num_particles,
			  fcs_float *positions, 
			  fcs_float *charges,
			  fcs_float *fields,
			  fcs_float *potentials) {
  FCS_INFO(fprintf(stderr, "ewald_compute_rspace started...\n"));

  fcs_int local_num_particles;
  fcs_int local_num_real_particles;
  fcs_int local_num_ghost_particles;
  fcs_float *local_positions, *local_ghost_positions;
  fcs_float *local_charges, *local_ghost_charges;
  fcs_gridsort_index_t *local_indices, *local_ghost_indices;
  fcs_gridsort_t gridsort;
  const fcs_float box_base[3] = {0.0, 0.0, 0.0 };
  const fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
  const fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
  const fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };

  /* DOMAIN DECOMPOSE */
  fcs_gridsort_create(&gridsort);
  fcs_gridsort_set_system(&gridsort, box_base, box_a, box_b, box_c, NULL);
  fcs_gridsort_set_particles(&gridsort, num_particles, max_num_particles,
			     positions, charges);

  FCS_INFO(fprintf(stderr, "  calling fcs_gridsort_sort_forward()...\n"));
  fcs_gridsort_sort_forward(&gridsort, d->r_cut, d->comm_cart);
  FCS_INFO(fprintf(stderr, "  returning from fcs_gridsort_sort_forward().\n"));

  fcs_gridsort_separate_ghosts(&gridsort);
  fcs_gridsort_get_sorted_particles(&gridsort, &local_num_particles, NULL,
				    NULL, NULL, NULL);
  fcs_gridsort_get_real_particles(&gridsort, &local_num_real_particles,
				  &local_positions, &local_charges, 
				  &local_indices);
  fcs_gridsort_get_ghost_particles(&gridsort, 
				   &local_num_ghost_particles, 
				   &local_ghost_positions, 
				   &local_ghost_charges, 
				   &local_ghost_indices);

  FCS_DEBUG_ALL(MPI_Barrier(d->comm_cart));
  FCS_DEBUG(fprintf(stderr,						\
		   "    %d: local_num_particles=%" FCS_LMOD_INT "d local_num_real_particles=%" FCS_LMOD_INT "d local_num_ghost_particles=%" FCS_LMOD_INT "d\n", \
		   d->comm_rank, local_num_particles,			\
		   local_num_real_particles, local_num_ghost_particles));

  /* allocate local fields and potentials */
  fcs_float *local_fields = NULL;
  fcs_float *local_potentials = NULL;
  if (fields != NULL) {
    local_fields = malloc(sizeof(fcs_float)*3*local_num_real_particles);
    for (fcs_int pid=0; pid < local_num_real_particles; pid++) {
      local_fields[3*pid] = 0.0;
      local_fields[3*pid+1] = 0.0;
      local_fields[3*pid+2] = 0.0;
    }
  }
  if (potentials != NULL) {
    local_potentials = malloc(sizeof(fcs_float)*local_num_real_particles);
    for (fcs_int pid=0; pid < local_num_real_particles; pid++)
      local_potentials[pid] = 0.0;
  }

  /* COMPUTE NEAR FIELD */
  fcs_near_t near;
  fcs_near_create(&near);
  fcs_near_set_loop(&near, ewald_compute_near_loop);
  fcs_near_set_system(&near, box_base, box_a, box_b, box_c, NULL);

  fcs_near_set_particles(&near, local_num_real_particles, local_num_real_particles,
			 local_positions, local_charges, 
			 local_indices,
			 (fields != NULL)?local_fields:NULL, 
			 (potentials != NULL)?local_potentials:NULL);
  
  fcs_near_set_ghosts(&near, local_num_ghost_particles,
		      local_ghost_positions, local_ghost_charges, 
		      local_ghost_indices);

  FCS_INFO(fprintf(stderr, "  calling fcs_near_compute()...\n"));
  fcs_near_compute(&near, d->r_cut, &d->alpha, d->comm_cart);
  FCS_INFO(fprintf(stderr, "  returning from fcs_near_compute().\n"));
  fcs_near_destroy(&near);

  /* RECOMPOSE FIELDS AND POTENTIALS */
  fcs_gridsort_set_sorted_results(&gridsort, local_num_real_particles, local_fields, local_potentials);
  fcs_gridsort_set_results(&gridsort, max_num_particles, fields, potentials);

  FCS_INFO(fprintf(stderr, "  calling fcs_gridsort_sort_backward()...\n"));
  fcs_gridsort_sort_backward(&gridsort, d->comm_cart);
  FCS_INFO(fprintf(stderr, "  returning from fcs_gridsort_sort_backward().\n"));
  
  fcs_gridsort_free(&gridsort);
  fcs_gridsort_destroy(&gridsort);

  FCS_INFO(fprintf(stderr, "ewald_compute_rspace finished.\n"));
}
