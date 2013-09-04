/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (C) 2012 Alexander Köwitsch
 * Copyright (C) 2011 Sebastian Banert
 *
 * This file is part of ScaFaCoS.
 * 
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *	
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


#include "run.h"
#include "types.h"
#include "utils.h"
#include "nearfield.h"
#include <common/near/near.h>
#include <fcs.h>
//#include "constants.h"

#define WORKAROUND_GRIDSORT_BUG 1
#define FCS_P2NFFT_DISABLE_PNFFT_INFO 1
#define CREATE_GHOSTS_SEPARATE 0

/* callback functions for performing a whole loop of near field computations (using ifcs_p2nfft_compute_near_...) */
static FCS_NEAR_LOOP_FP(ifcs_p2nfft_compute_near_periodic_erfc_loop, ifcs_p2nfft_compute_near_periodic_erfc)
static FCS_NEAR_LOOP_FP(ifcs_p2nfft_compute_near_periodic_approx_erfc_loop, ifcs_p2nfft_compute_near_periodic_approx_erfc)
// static FCS_NEAR_LOOP_FP(ifcs_p2nfft_compute_near_interpolation_loop, ifcs_p2nfft_compute_near_interpolation)
static FCS_NEAR_LOOP_FP(ifcs_p2nfft_compute_near_interpolation_const_loop, ifcs_p2nfft_compute_near_interpolation_const)
static FCS_NEAR_LOOP_FP(ifcs_p2nfft_compute_near_interpolation_lin_loop, ifcs_p2nfft_compute_near_interpolation_lin)
static FCS_NEAR_LOOP_FP(ifcs_p2nfft_compute_near_interpolation_quad_loop, ifcs_p2nfft_compute_near_interpolation_quad)
static FCS_NEAR_LOOP_FP(ifcs_p2nfft_compute_near_interpolation_cub_loop, ifcs_p2nfft_compute_near_interpolation_cub)

static void convolution(
    const INT *local_N, const fcs_pnfft_complex *regkern_hat,
    fcs_pnfft_complex *f_hat);


FCSResult ifcs_p2nfft_run(
    void *rd, fcs_int local_num_particles, fcs_int max_local_num_particles,
    fcs_float *positions, fcs_float *charges,
    fcs_float *potentials, fcs_float *field
    )
{
  const char* fnc_name = "ifcs_p2nfft_run";
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  C csum;
  C csum_global;
  fcs_float rsum;
  fcs_float rsum_global;
#endif

#if FCS_ENABLE_INFO || FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  int myrank;
  MPI_Comm_rank(d->cart_comm_3d, &myrank);
#endif

  FCS_P2NFFT_INIT_TIMING(d->cart_comm_3d);

  /* Tuning was called in fcs_run */
//  result = ifcs_p2nfft_tune(rd, local_num_particles, positions, charges, d->box_l[0]);

  /* handle particles, that left the box [0,L] */
  /* for non-periodic boundary conditions: user must increase the box */
  for(fcs_int j=0; j<local_num_particles; j++)
    for(fcs_int t=0; t<3; t++)
      if(!d->periodicity[t]) /* for mixed periodicity: only handle the non-periodic dimensions */
        if( (positions[3*j+t] < 0) || (positions[3*j+t] > d->box_l[t]) )
//        if( (positions[3*j+t] < 0) || fcs_float_is_zero(d->box_[t] - positions[3*j+t])  || (positions[3*j+t] > d->box_l[t]) )
          return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Box size does not fit. Some particles left the box or reached the upper border.");
  /* TODO: implement additional scaling of particles to ensure x \in [0,L)
   * Idea: use allreduce to get min and max coordinates, adapt scaling of particles for every time step */

  /* Start forw sort timing */
  FCS_P2NFFT_START_TIMING();
  
#if WORKAROUND_GRIDSORT_BUG
  fcs_float lo[3], up[3];
  for(int t=0; t<3; t++){
    lo[t] = d->lower_border[t];
    up[t] = d->upper_border[t];
    if(!d->periodicity[t]){ /* for mixed periodicity: only handle the non-periodic dimensions */
      if(fcs_float_is_zero(lo[t]-d->box_l[t]))
        lo[t] += 0.1;
      if(fcs_float_is_zero(up[t]-d->box_l[t]))
        up[t] += 0.1;
    }
  }
#endif
  
  /* Compute the near field */
  fcs_float box_a[3] = { d->box_l[0], 0, 0 };
  fcs_float box_b[3] = { 0, d->box_l[1], 0 };
  fcs_float box_c[3] = { 0, 0, d->box_l[2] };
  fcs_float box_base[3] = { 0, 0, 0 };
  fcs_near_t near;

#if WORKAROUND_GRIDSORT_BUG
  /* for mixed periodicity: only handle the non-periodic dimensions */
  if(!d->periodicity[0]) box_a[0] += 0.1;
  if(!d->periodicity[1]) box_b[1] += 0.1;
  if(!d->periodicity[2]) box_c[2] += 0.1;
#endif

  fcs_int sorted_num_particles, ghost_num_particles;
  fcs_float *sorted_positions, *ghost_positions;
  fcs_float *sorted_charges, *ghost_charges;
  fcs_gridsort_index_t *sorted_indices, *ghost_indices;
  fcs_gridsort_t gridsort;

  fcs_gridsort_create(&gridsort);
  
  fcs_gridsort_set_system(&gridsort, box_base, box_a, box_b, box_c, d->periodicity);
  
#if WORKAROUND_GRIDSORT_BUG
  fcs_gridsort_set_bounds(&gridsort, lo, up);
#else
  fcs_gridsort_set_bounds(&gridsort, d->lower_border, d->upper_border);
#endif

  fcs_gridsort_set_particles(&gridsort, local_num_particles, max_local_num_particles, positions, charges);

  fcs_gridsort_set_max_particle_move(&gridsort, d->max_particle_move);

  fcs_gridsort_set_cache(&gridsort, &d->gridsort_cache);

#if CREATE_GHOSTS_SEPARATE
  fcs_gridsort_sort_forward(&gridsort, 0, d->cart_comm_3d);
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Forward grid sort");

  /* Start near sort timing */
  FCS_P2NFFT_START_TIMING();
  if (d->short_range_flag) fcs_gridsort_create_ghosts(&gridsort, d->r_cut, d->cart_comm_3d);
#else
  fcs_gridsort_sort_forward(&gridsort, (d->short_range_flag ? d->r_cut: 0.0), d->cart_comm_3d);
#endif

  fcs_gridsort_separate_ghosts(&gridsort, NULL, NULL);

  fcs_gridsort_get_real_particles(&gridsort, &sorted_num_particles, &sorted_positions, &sorted_charges, &sorted_indices);
  fcs_gridsort_get_ghost_particles(&gridsort, &ghost_num_particles, &ghost_positions, &ghost_charges, &ghost_indices);

  /* Finish forw sort timing */
#if CREATE_GHOSTS_SEPARATE
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Forward near sort");
#else
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Forward grid and near sort");
#endif

  /* Handle particles, that left the box [0,L] */
  /* For periodic boundary conditions: just fold them back.
   * We change sorted_positions (and not positions), since we are allowed to overwrite them. */
  if(d->periodicity[0] || d->periodicity[1] || d->periodicity[2]){
    for(fcs_int j=0; j<sorted_num_particles; j++){
      for(fcs_int t=0; t<3; t++){
        if(d->periodicity[t]){
          while(sorted_positions[3*j+t] < 0)
            sorted_positions[3*j+t] += d->box_l[t];
          while(sorted_positions[3*j+t] >= d->box_l[t])
            sorted_positions[3*j+t] -= d->box_l[t];
        }
      }
    }
  }

  /* Start near field timing */
  FCS_P2NFFT_START_TIMING();

/*  printf("%d: input number = %" FCS_LMOD_INT "d, sorted number = %" FCS_LMOD_INT "d, ghost number = %" FCS_LMOD_INT "d\n",
    myrank, local_num_particles, sorted_num_particles, ghost_num_particles);*/

  fcs_int compute_field = (field != NULL);
  fcs_int compute_potentials = (potentials != NULL);

  fcs_float *sorted_field = (compute_field) ? malloc(sizeof(fcs_float)*3*sorted_num_particles) : NULL;
  fcs_float *sorted_potentials = (compute_potentials) ? malloc(sizeof(fcs_float)*sorted_num_particles) : NULL;

  /* Initialize all the potentials */
  for (fcs_int j = 0; j < sorted_num_particles; ++j)
    sorted_potentials[j] = 0;
  
  /* Initialize all the forces */
  for (fcs_int j = 0; j < 3 * sorted_num_particles; ++j)
    sorted_field[j] = 0;

  if(d->short_range_flag){
    fcs_near_create(&near);
  
    if(d->interpolation_order >= 0){
      switch(d->interpolation_order){
        case 0: fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_interpolation_const_loop); break;
        case 1: fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_interpolation_lin_loop); break;
        case 2: fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_interpolation_quad_loop); break;
        case 3: fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_interpolation_cub_loop); break;
        default: return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name,"P2NFFT interpolation order is too large.");
      } 
    } else if(d->use_ewald){
      if(d->interpolation_order == -1)
        fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_periodic_erfc_loop);
      else
        fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_periodic_approx_erfc_loop);
    } else {
      fcs_near_set_field(&near, ifcs_p2nfft_compute_near_field);
      fcs_near_set_potential(&near, ifcs_p2nfft_compute_near_potential);
    }

    // fcs_int periodicity[3] = { 1, 1, 1 };
    // fcs_int periodicity[3] = { 0, 0, 0 };
    // fcs_int *periodicity = NULL; /* sorter uses periodicity of the communicator */
    fcs_near_set_system(&near, box_base, box_a, box_b, box_c, d->periodicity);
  
    fcs_near_set_particles(&near, sorted_num_particles, sorted_num_particles, sorted_positions, sorted_charges, sorted_indices,
        (compute_field)?sorted_field:NULL, (compute_potentials)?sorted_potentials:NULL);
  
    fcs_near_set_ghosts(&near, ghost_num_particles, ghost_positions, ghost_charges, ghost_indices);
  
    if(d->interpolation_order >= 0){
      ifcs_p2nfft_near_params near_params;
      near_params.interpolation_order = d->interpolation_order;
      near_params.interpolation_num_nodes = d->near_interpolation_num_nodes;
      near_params.near_interpolation_table_potential = d->near_interpolation_table_potential;
      near_params.near_interpolation_table_force = d->near_interpolation_table_force;
      near_params.one_over_r_cut = d->one_over_r_cut;

      fcs_near_compute(&near, d->r_cut, &near_params, d->cart_comm_3d);
    } else if(d->use_ewald)
      fcs_near_compute(&near, d->r_cut, &(d->alpha), d->cart_comm_3d);
    else
      fcs_near_compute(&near, d->r_cut, rd, d->cart_comm_3d);
  
    fcs_near_destroy(&near);

  }

  /* Finish near field timing */
//  tm_timer += MPI_Wtime();
//  printf("P2NFFT_TIMING: rank = %d, Near field computation takes %e s\n", myrank, tm_timer);
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Near field computation");
  
  /* Checksum: global sum of nearfield energy */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  fcs_float near_energy = 0.0;
  fcs_float near_global;
  for(fcs_int j = 0; j < sorted_num_particles; ++j) {
    near_energy += 0.5 * sorted_charges[j] * sorted_potentials[j];
  }
  MPI_Reduce(&near_energy, &near_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: near field energy: %" FCS_LMOD_FLOAT "f\n", near_global);
#endif

  /* Checksum: fields resulting from nearfield interactions */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  for(fcs_int t=0; t<3; t++){
    rsum = 0.0;
    for(fcs_int j = 0; (j < sorted_num_particles); ++j)
      rsum += fabs(sorted_field[3*j+t]);
    MPI_Reduce(&rsum, &rsum_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
    if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: near field %" FCS_LMOD_INT "d. component: %" FCS_LMOD_FLOAT "f\n", t, rsum_global);
  }
  
  if (myrank == 0) fprintf(stderr, "E_NEAR(0) = %" FCS_LMOD_FLOAT "e\n", sorted_field[0]);
#endif
      
  /* Reinit PNFFT corresponding to number of sorted nodes */
  FCS_PNFFT(init_nodes)(d->pnfft, sorted_num_particles,
      PNFFT_MALLOC_X| PNFFT_MALLOC_F| PNFFT_MALLOC_GRAD_F,
      PNFFT_FREE_X|   PNFFT_FREE_F|   PNFFT_FREE_GRAD_F);

  fcs_pnfft_complex *f_hat, *f, *grad_f;
  fcs_float *x;

  f_hat  = FCS_PNFFT(get_f_hat)(d->pnfft);
  f      = FCS_PNFFT(get_f)(d->pnfft);
  grad_f = FCS_PNFFT(get_grad_f)(d->pnfft);
  x      = FCS_PNFFT(get_x)(d->pnfft);

// #if FCS_ENABLE_INFO
//   fcs_float min[3], max[3], gmin[3], gmax[3];
// 
//   /* initialize */
//   for(fcs_int t=0; t<3; t++)
//     min[t] = (sorted_num_particles >0) ? sorted_positions[t] : 1e16;
//   for(fcs_int t=0; t<3; t++)
//     max[t] = (sorted_num_particles >0) ? sorted_positions[t] : -1e16;
// 
//   for (fcs_int j = 1; j < sorted_num_particles; ++j){
//     for(fcs_int t=0; t<3; t++){
//       if(sorted_positions[3*j+t] < min[t])
// 	min[t] = sorted_positions[3*j+t];
//       if(sorted_positions[3*j+t] > max[t])
// 	max[t] = sorted_positions[3*j+t];
//     }
//   }
//       
//   MPI_Reduce(&min, &gmin, 3, FCS_MPI_FLOAT, MPI_MIN, 0, d->cart_comm_3d);
//   MPI_Reduce(&max, &gmax, 3, FCS_MPI_FLOAT, MPI_MAX, 0, d->cart_comm_3d);
//   if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: min range of particles: (%e, %e, %e)\n", gmin[0], gmin[1], gmin[2]);
//   if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: max range of particles: (%e, %e, %e)\n", gmax[0], gmax[1], gmax[2]);
//   fprintf(stderr, "myrank = %d, sorted_num_particles = %d\n", myrank, sorted_num_particles);
// #endif
  
  /* Set the NFFT nodes and values */
  for (fcs_int j = 0; j < sorted_num_particles; ++j)
  {
    x[3 * j + 0] = (sorted_positions[3 * j + 0] - d->box_shifts[0]) / d->box_scales[0];
    x[3 * j + 1] = (sorted_positions[3 * j + 1] - d->box_shifts[1]) / d->box_scales[1];
    x[3 * j + 2] = (sorted_positions[3 * j + 2] - d->box_shifts[2]) / d->box_scales[2];
    
    f[j] = sorted_charges[j];
  }

// #if FCS_ENABLE_INFO
//   fcs_float min[3], max[3], gmin[3], gmax[3];
// 
//   /* initialize */
//   for(fcs_int t=0; t<3; t++)
//     min[t] = (sorted_num_particles >0) ? x[t] : 1e16;
//   for(fcs_int t=0; t<3; t++)
//     max[t] = (sorted_num_particles >0) ? x[t] : -1e16;
// 
//   for (fcs_int j = 1; j < sorted_num_particles; ++j){
//     for(fcs_int t=0; t<3; t++){
//       if(x[3*j+t] < min[t])
// 	min[t] = x[3*j+t];
//       if(x[3*j+t] > max[t])
// 	max[t] = x[3*j+t];
//     }
//   }
//       
//   MPI_Reduce(&min, &gmin, 3, FCS_MPI_FLOAT, MPI_MIN, 0, d->cart_comm_3d);
//   MPI_Reduce(&max, &gmax, 3, FCS_MPI_FLOAT, MPI_MAX, 0, d->cart_comm_3d);
//   if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: min range of particles: (%e, %e, %e)\n", gmin[0], gmin[1], gmin[2]);
//   if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: max range of particles: (%e, %e, %e)\n", gmax[0], gmax[1], gmax[2]);
//   fprintf(stderr, "myrank = %d, sorted_num_particles = %d\n", myrank, sorted_num_particles);
// #endif

  FCS_P2NFFT_START_TIMING();
  FCS_PNFFT(precompute_psi)(d->pnfft);
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "pnfft_precompute_psi");

  /* Reset pnfft timer (delete timings from fcs_init and fcs_tune) */  
#if FCS_ENABLE_INFO && !FCS_P2NFFT_DISABLE_PNFFT_INFO
  FCS_PNFFT(reset_timer)(d->pnfft);
#endif

  /* Checksum: Input of adjoint NFFT */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  rsum = 0.0;
  for (fcs_int j = 0; j < 3*sorted_num_particles; ++j)
    rsum += fabs(x[j]);
  MPI_Reduce(&rsum, &rsum_global, 1, MPI_DOUBLE, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: checksum of x: %" FCS_LMOD_FLOAT "e\n", rsum_global);

  csum = 0.0;
  for (fcs_int j = 0; j < sorted_num_particles; ++j)
    csum += fabs(creal(f[j])) + _Complex_I * fabs(cimag(f[j]));
  MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: checksum of NFFT^H input: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif

  /* Start far field timing */
  FCS_P2NFFT_START_TIMING();

  /* Perform adjoint NFFT */
  FCS_PNFFT(adj)(d->pnfft);

  /* Checksum: Output of adjoint NFFT */  
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  csum = 0.0;
  for(fcs_int k = 0; k < d->local_N[0]*d->local_N[1]*d->local_N[2]; ++k)
     csum += fabs(creal(f_hat[k])) + _Complex_I * fabs(cimag(f_hat[k]));
  MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: checksum of Fourier coefficients before convolution: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif

  /* Checksum: Fourier coefficients of regkernel  */  
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  csum = 0.0;
  for(fcs_int k = 0; k < d->local_N[0]*d->local_N[1]*d->local_N[2]; ++k)
     csum += fabs(creal(d->regkern_hat[k])) + _Complex_I * fabs(cimag(d->regkern_hat[k]));
  MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: checksum of Regkernel Fourier coefficients: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif

  /* Multiply with the analytically given Fourier coefficients */
  convolution(d->local_N, d->regkern_hat,
      f_hat);

  /* Checksum: Input of NFFT */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  csum = 0.0;
  for(fcs_int k = 0; k < d->local_N[0]*d->local_N[1]*d->local_N[2]; ++k)
     csum += fabs(creal(f_hat[k])) + _Complex_I * fabs(cimag(f_hat[k]));
  MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: checksum of Fourier coefficients after convolution: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif
    
  /* Perform NFFT */
  FCS_PNFFT(trafo)(d->pnfft);

  fcs_float box_surf = 1.0;
  for(fcs_int t=0; t<3; t++)
    if(d->periodicity[t])
      box_surf *= d->box_scales[t];

  /* Copy the results to the output vector and rescale */
  for (fcs_int j = 0; j < sorted_num_particles; ++j){
    if(d->use_ewald){
      sorted_potentials[j] += creal(f[j]) / box_surf;
      for(fcs_int t=0; t<3; t++)
        sorted_field[3 * j + t] -= creal(grad_f[3 * j + t]) / (box_surf * d->box_scales[t]);
    } else {
      sorted_potentials[j] += creal(f[j]);
      for(fcs_int t=0; t<3; t++)
        sorted_field[3 * j + t] -= creal(grad_f[3 * j + t]) / d->box_scales[t];
    }
  }


  /* Checksum: fields resulting from farfield interactions */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  for(fcs_int t=0; t<3; t++){
    rsum = 0.0;
    for(fcs_int j = 0; (j < sorted_num_particles); ++j)
      rsum += fabs(sorted_field[3*j+t]);
    MPI_Reduce(&rsum, &rsum_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
    if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: near plus far field %" FCS_LMOD_INT "d. component: %" FCS_LMOD_FLOAT "f\n", t, rsum_global);
  }

  if (myrank == 0) fprintf(stderr, "E_NEAR_FAR(0) = %" FCS_LMOD_FLOAT "e\n", sorted_field[0]);
#endif

  /* Finish far field timing */
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Far field computation");

#if FCS_ENABLE_TIMING
  /* Print pnfft timer */
  FCS_P2NFFT_START_TIMING();
  FCS_PNFFT(print_average_timer_adv)(d->pnfft, d->cart_comm_3d);
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Printing of PNFFT timings");
#endif

  /* Checksum: global sum of farfield energy */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  fcs_float far_energy = 0.0;
  fcs_float far_global;
  for(fcs_int j = 0; j < sorted_num_particles; ++j)
    if(d->use_ewald)
      far_energy += 0.5 * sorted_charges[j] * f[j] / box_surf;
    else
      far_energy += 0.5 * sorted_charges[j] * f[j];

  MPI_Reduce(&far_energy, &far_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: far field energy: %" FCS_LMOD_FLOAT "f\n", far_global);
#endif
  
  /* Start self interaction timing */
  FCS_P2NFFT_START_TIMING();

  /* Calculate self-interactions */
  for (fcs_int j = 0; j < sorted_num_particles; ++j)
    sorted_potentials[j] -= sorted_charges[j] * ifcs_p2nfft_compute_self_potential(rd);

  /* Finish self interaction timing */
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "self interaction calculation");

  /* Calculate virial if needed */
  if(d->virial != NULL){
    if(d->use_ewald){
      fcs_float total_energy = 0.0;
      fcs_float total_global;
      for(fcs_int j = 0; j < sorted_num_particles; ++j)
        total_energy += 0.5 * sorted_charges[j] * sorted_potentials[j];

      MPI_Allreduce(&total_energy, &total_global, 1, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);

      /* Approximate virial in 3d-periodic case:
       * Fill the main diagonal with one third of the total energy */      
      for(fcs_int t=0; t<9; t++)
        d->virial[t] = 0.0;
      d->virial[0] = d->virial[4] = d->virial[8] = total_global/3.0;
    }
    else {
      fcs_float local_virial[9];

      for(fcs_int t=0; t<9; t++)
        local_virial[t] = 0.0;
      
      for(fcs_int j = 0; j < sorted_num_particles; ++j)
        for(fcs_int t0=0; t0<3; t0++)
          for(fcs_int t1=0; t1<3; t1++)
            local_virial[t0*3+t1] += sorted_charges[j] * sorted_field[3*j+t0] * sorted_positions[3*j+t1];
      
      MPI_Allreduce(local_virial, d->virial, 9, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);
    }
  }

  /* Checksum: global sum of self energy */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  fcs_float self_energy = 0.0;
  fcs_float self_global;
  for(fcs_int j = 0; j < sorted_num_particles; ++j)
    self_energy -= 0.5 * sorted_charges[j] * sorted_charges[j] * ifcs_p2nfft_compute_self_potential(rd);

  MPI_Reduce(&self_energy, &self_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: self energy: %" FCS_LMOD_FLOAT "f\n", self_global);
#endif
      
  /* Checksum: global sum of total energy */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  fcs_float total_energy = 0.0;
  fcs_float total_global;
  for(fcs_int j = 0; j < sorted_num_particles; ++j)
    total_energy += 0.5 * sorted_charges[j] * sorted_potentials[j];

  MPI_Reduce(&total_energy, &total_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: total energy: %" FCS_LMOD_FLOAT "f\n", total_global);
#endif

  /* Try: calculate total dipol moment */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  fcs_float total_dipol_local[3] = {0.0, 0.0, 0.0};
  fcs_float total_dipol_global[3];
  for(fcs_int j = 0; j < sorted_num_particles; ++j)
    for(fcs_int t=0; t<3; ++t)
//       total_dipol_local[t] += sorted_charges[j] * sorted_positions[3*j+t];
      total_dipol_local[t] += sorted_charges[j] * (sorted_positions[3*j+t] - d->box_l[t]);

  MPI_Allreduce(&total_dipol_local, &total_dipol_global, 3, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: total dipol: (%" FCS_LMOD_FLOAT "e, %" FCS_LMOD_FLOAT "e, %" FCS_LMOD_FLOAT "e)\n", total_dipol_global[0], total_dipol_global[1], total_dipol_global[2]);

  for(fcs_int t=0; t<3; ++t)
//     total_dipol_global[t] *= 4.0*PNFFT_PI/3.0 / (d->box_l[0]*d->box_l[1]*d->box_l[2]);
//     total_dipol_global[t] *= PNFFT_PI / (d->box_l[0]*d->box_l[1]*d->box_l[2]);
    total_dipol_global[t] *= PNFFT_PI/3.0 *PNFFT_PI / (d->box_l[0]*d->box_l[1]*d->box_l[2]);
//     total_dipol_global[t] *= 2.0*PNFFT_PI / (d->box_l[0]*d->box_l[1]*d->box_l[2]);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: scaled total dipol: (%" FCS_LMOD_FLOAT "e, %" FCS_LMOD_FLOAT "e, %" FCS_LMOD_FLOAT "e)\n", total_dipol_global[0], total_dipol_global[1], total_dipol_global[2]);
#endif
      
  /* Start back sort timing */
  FCS_P2NFFT_START_TIMING();

  fcs_int resort;

  if (d->resort) resort = fcs_gridsort_prepare_resort(&gridsort, sorted_field, sorted_potentials, field, potentials, d->cart_comm_3d);
  else resort = 0;

  /* Backsort data into user given ordering (if resort is disabled) */
  if (!resort) fcs_gridsort_sort_backward(&gridsort, sorted_field, sorted_potentials, field, potentials, 1, d->cart_comm_3d);

  fcs_gridsort_resort_destroy(&d->gridsort_resort);

  if (d->resort) fcs_gridsort_resort_create(&d->gridsort_resort, &gridsort, d->cart_comm_3d);
  
  d->local_num_particles = local_num_particles;

  if (sorted_field) free(sorted_field);
  if (sorted_potentials) free(sorted_potentials);

  fcs_gridsort_free(&gridsort);

  fcs_gridsort_destroy(&gridsort);

  /* Finish back sort timing */
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Backward sort");

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  /* print first value of fields */
  if (myrank == 0) printf("P2NFFT_INFO: E(0) = %e\n", creal(field[0]));
  if (myrank == 0) printf("P2NFFT_INFO: dE(0) = %e\n", creal(field[0])+1.619834832399799e-06);
#endif

  return NULL;
}

static void convolution(
    const INT *local_N, const fcs_pnfft_complex *regkern_hat,
    fcs_pnfft_complex *f_hat
    )
{  
  INT m=0;

  for (INT k0 = 0; k0 < local_N[0]; ++k0)
    for (INT k1 = 0; k1 < local_N[1]; ++k1)
      for (INT k2 = 0; k2 < local_N[2]; ++k2, ++m){
//         fprintf(stderr, "f_hat[%td, %td, %td] = %e + I * %e, regkern_hat[%td, %td, %td] = %e I * %e\n", k0, k1, k2, creal(f_hat[m]), cimag(f_hat[m]), k0, k1, k2, creal(regkern_hat[m]), cimag(regkern_hat[m]));
//         fprintf(stderr, "regkern_hat[%td, %td, %td] = %e I * %e\n", k0, k1, k2, creal(regkern_hat[m]), cimag(regkern_hat[m]));
        f_hat[m] *= regkern_hat[m];
      }
}


void ifcs_p2nfft_set_max_particle_move(void *rd, fcs_float max_particle_move)
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;

  d->max_particle_move = max_particle_move;
}

void ifcs_p2nfft_set_resort(void *rd, fcs_int resort)
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;

  d->resort = resort;
}

void ifcs_p2nfft_get_resort(void *rd, fcs_int *resort)
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;

  *resort = d->resort;
}

void ifcs_p2nfft_get_resort_availability(void *rd, fcs_int *availability)
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;

  *availability = fcs_gridsort_resort_is_available(d->gridsort_resort);
}

void ifcs_p2nfft_get_resort_particles(void *rd, fcs_int *resort_particles)
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;

  if (d->gridsort_resort == FCS_GRIDSORT_RESORT_NULL)
  {
    *resort_particles = d->local_num_particles;
    return;
  }
  
  *resort_particles = fcs_gridsort_resort_get_sorted_particles(d->gridsort_resort);
}

void ifcs_p2nfft_resort_ints(void *rd, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;

  if (d->gridsort_resort == FCS_GRIDSORT_RESORT_NULL) return;
  
  fcs_gridsort_resort_ints(d->gridsort_resort, src, dst, n, comm);
}

void ifcs_p2nfft_resort_floats(void *rd, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;

  if (d->gridsort_resort == FCS_GRIDSORT_RESORT_NULL) return;
  
  fcs_gridsort_resort_floats(d->gridsort_resort, src, dst, n, comm);
}

void ifcs_p2nfft_resort_bytes(void *rd, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;

  if (d->gridsort_resort == FCS_GRIDSORT_RESORT_NULL) return;
  
  fcs_gridsort_resort_bytes(d->gridsort_resort, src, dst, n, comm);
}
