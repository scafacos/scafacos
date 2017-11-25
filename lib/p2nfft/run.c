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
#include "interpolation.h"
//#include "constants.h"
#include <common/near/near.h>

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

#define SPROD3(_u_, _v_) ( (_u_)[0] * (_v_)[0] + (_u_)[1] * (_v_)[1] + (_u_)[2] * (_v_)[2] ) 
#define XYZ2TRI(_d_, _x_, _ib_) ( SPROD3((_ib_) + 3*(_d_), (_x_)) )
/* compute d-th component of A^T * v */
#define At_TIMES_VEC(_A_, _v_, _d_) ( (_v_)[0] * (_A_)[_d_] + (_v_)[1] * (_A_)[_d_ + 3] + (_v_)[2] * (_A_)[_d_ + 6] )

#define At_TIMES_SYMMAT_TIMES_A(A, S, i, j) \
  (   (A)[  i] * ( (S)[0] * (A)[j] + (S)[1] * (A)[3+j] + (S)[2] * (A)[6+j]) \
    + (A)[3+i] * ( (S)[1] * (A)[j] + (S)[3] * (A)[3+j] + (S)[4] * (A)[6+j]) \
    + (A)[6+i] * ( (S)[2] * (A)[j] + (S)[4] * (A)[3+j] + (S)[5] * (A)[6+j]) )

static fcs_int box_not_large_enough(
    fcs_int npart, const fcs_float *pos_with_offset, const fcs_float *box_base, const fcs_float *ibox, const fcs_int *periodicity
    )
{
  fcs_float pos[3];
  for(fcs_int j=0; j<npart; j++)
  {
    pos[0] = pos_with_offset[3*j + 0] - box_base[0];
    pos[1] = pos_with_offset[3*j + 1] - box_base[1];
    pos[2] = pos_with_offset[3*j + 2] - box_base[2];

    if( (periodicity[0] == 0) && (XYZ2TRI(0, pos, ibox) < 0.0) ) return 1;
    if( (periodicity[0] == 0) && (XYZ2TRI(0, pos, ibox) > 1.0) ) return 1;

    if( (periodicity[1] == 0) && (XYZ2TRI(1, pos, ibox) < 0.0) ) return 1;
    if( (periodicity[1] == 0) && (XYZ2TRI(1, pos, ibox) > 1.0) ) return 1;

    if( (periodicity[2] == 0) && (XYZ2TRI(2, pos, ibox) < 0.0) ) return 1;
    if( (periodicity[2] == 0) && (XYZ2TRI(2, pos, ibox) > 1.0) ) return 1;
  }

  return 0;
}




FCSResult ifcs_p2nfft_run(
    void *rd, 
    fcs_int local_num_particles, fcs_int max_local_num_particles,
    fcs_float *positions, fcs_float *charges,
    fcs_float *potential, fcs_float *field,
    fcs_int local_num_dipole_particles, fcs_int max_local_num_dipole_particles,
    fcs_float *dipole_positions, fcs_float *dipole_moments,
    fcs_float *dipole_potential, fcs_float *dipole_field
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


  /* handle particles, that left the box [0,L] */
  /* for non-periodic boundary conditions: user must increase the box */
  if(box_not_large_enough(local_num_particles, positions, d->box_base, d->box_inv, d->periodicity))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Box size does not fit. Some particles left the box.");

  if(box_not_large_enough(local_num_dipole_particles, dipole_positions, d->box_base, d->box_inv, d->periodicity))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Box size does not fit. Some dipole particles left the box.");

  /* TODO: implement additional scaling of particles to ensure x \in [0,L)
   * Idea: use allreduce to get min and max coordinates, adapt scaling of particles for every time step */

  /* Start forw sort timing */
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);
  
  /* Compute the near field */
  fcs_near_t near;

  fcs_int sorted_num_particles, ghost_num_particles;
  fcs_float *sorted_positions, *ghost_positions;
  fcs_float *sorted_charges, *ghost_charges;
  fcs_gridsort_index_t *sorted_indices, *ghost_indices;
  fcs_int sorted_num_dipole_particles, ghost_num_dipole_particles;
  fcs_float *sorted_dipole_positions, *ghost_dipole_positions;
  fcs_float *sorted_dipole_moments, *ghost_dipole_moments;
  fcs_gridsort_index_t *sorted_dipole_indices, *ghost_dipole_indices;
  fcs_gridsort_t gridsort;

  fcs_gridsort_create(&gridsort);
  
  fcs_gridsort_set_system(&gridsort, d->box_base, d->ebox_a, d->ebox_b, d->ebox_c, d->periodicity);

  fcs_gridsort_set_bounds(&gridsort, d->lower_border, d->upper_border);

  fcs_gridsort_set_particles(&gridsort, local_num_particles, max_local_num_particles, positions, charges);

  fcs_gridsort_set_dipole_particles(&gridsort, local_num_dipole_particles, max_local_num_dipole_particles, dipole_positions, dipole_moments);

  fcs_gridsort_set_max_particle_move(&gridsort, d->max_particle_move);

  fcs_gridsort_set_cache(&gridsort, &d->gridsort_cache);

#if CREATE_GHOSTS_SEPARATE
  fcs_gridsort_sort_forward(&gridsort, 0, d->cart_comm_3d);
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Forward grid sort");

  /* Start near sort timing */
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);
  if (d->short_range_flag) fcs_gridsort_create_ghosts(&gridsort, d->r_cut, d->cart_comm_3d);
#else
  fcs_gridsort_sort_forward(&gridsort, (d->short_range_flag ? d->r_cut: 0.0), d->cart_comm_3d);
#endif

  fcs_gridsort_separate_ghosts(&gridsort);

  fcs_gridsort_get_real_particles(&gridsort, &sorted_num_particles, &sorted_positions, &sorted_charges, &sorted_indices);
  fcs_gridsort_get_real_dipole_particles(&gridsort, &sorted_num_dipole_particles, &sorted_dipole_positions, &sorted_dipole_moments, &sorted_dipole_indices);

  fcs_gridsort_get_ghost_particles(&gridsort, &ghost_num_particles, &ghost_positions, &ghost_charges, &ghost_indices);
  fcs_gridsort_get_ghost_dipole_particles(&gridsort, &ghost_num_dipole_particles, &ghost_dipole_positions, &ghost_dipole_moments, &ghost_dipole_indices);

  /* Finish forw sort timing */
#if CREATE_GHOSTS_SEPARATE
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Forward near sort");
#else
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Forward grid and near sort");
#endif

  /* Handle particles, that left the box [0,L] */
  /* For periodic boundary conditions: just fold them back.
   * We change sorted_positions (and not positions), since we are allowed to overwrite them. */
  fcs_wrap_positions(sorted_num_particles, sorted_positions, d->box_a, d->box_b, d->box_c, d->box_base, d->periodicity);
  fcs_wrap_positions(sorted_num_dipole_particles, sorted_dipole_positions, d->box_a, d->box_b, d->box_c, d->box_base, d->periodicity);

  /* Start near field timing */
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);

/*  printf("%d: input number = %" FCS_LMOD_INT "d, sorted number = %" FCS_LMOD_INT "d, ghost number = %" FCS_LMOD_INT "d\n",
    myrank, local_num_particles, sorted_num_particles, ghost_num_particles);*/

  /* additional switchs to turn off computation of field / potential */
  if(d->flags & FCS_P2NFFT_IGNORE_FIELD)     field     = NULL;
  if(d->flags & FCS_P2NFFT_IGNORE_POTENTIAL) potential = NULL;
  if(d->flags & FCS_P2NFFT_IGNORE_FIELD)     dipole_field     = NULL;
  if(d->flags & FCS_P2NFFT_IGNORE_POTENTIAL) dipole_potential = NULL;

  fcs_int compute_field     = (field != NULL);
  fcs_int compute_potential = (potential != NULL);
  fcs_int compute_dipole_field     = (dipole_field != NULL);
  fcs_int compute_dipole_potential = (dipole_potential != NULL);

  fcs_float *sorted_field     = (compute_field)     ? malloc(sizeof(fcs_float)*3*sorted_num_particles) : NULL;
  fcs_float *sorted_potential = (compute_potential) ? malloc(sizeof(fcs_float)*sorted_num_particles) : NULL;
  fcs_float *sorted_dipole_field     = (compute_dipole_field)     ? malloc(sizeof(fcs_float)*6*sorted_num_dipole_particles) : NULL;
  fcs_float *sorted_dipole_potential = (compute_dipole_potential) ? malloc(sizeof(fcs_float)*3*sorted_num_dipole_particles) : NULL;

  /* Initialize all the potential */
  if(compute_potential)
    for (fcs_int j = 0; j < sorted_num_particles; ++j)
      sorted_potential[j] = 0;

  if(compute_dipole_potential)
    for (fcs_int j = 0; j < 3 * sorted_num_dipole_particles; ++j)
      sorted_dipole_potential[j] = 0;
  
  /* Initialize all the forces */
  if(compute_field)
    for (fcs_int j = 0; j < 3 * sorted_num_particles; ++j)
      sorted_field[j] = 0;

  if(compute_dipole_field)
    for (fcs_int j = 0; j < 6 * sorted_num_dipole_particles; ++j)
      sorted_dipole_field[j] = 0;

  if(d->short_range_flag){
    fcs_near_create(&near);

    if(d->interpolation_order >= 0 && !compute_dipole_potential && !compute_dipole_field ){
      /* interpolation loop only works for charge-charge interactions */
      /* TODO: implement interpolation and loops for dipoles */
      switch(d->interpolation_order){
        case 0: fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_interpolation_const_loop); break;
        case 1: fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_interpolation_lin_loop); break;
        case 2: fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_interpolation_quad_loop); break;
        case 3: fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_interpolation_cub_loop); break;
        default: return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,"P2NFFT interpolation order is too large.");
      } 
    } else if(d->reg_kernel == FCS_P2NFFT_REG_KERNEL_EWALD &&  !compute_dipole_potential && !compute_dipole_field ){
      /* near field loop only works for charge-charge interactions */
      /* TODO: implement approx erfc interactions for dipoles */
      if(d->interpolation_order < -1 && !compute_dipole_potential && !compute_dipole_field )
        fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_periodic_approx_erfc_loop);
      else /* d->interpolation_order = -1 */
        fcs_near_set_loop(&near, ifcs_p2nfft_compute_near_periodic_erfc_loop);
    } else {
      /* set scalar functions for dipole interactions */
      fcs_near_set_charge_charge(&near, ifcs_p2nfft_compute_near_charge_charge);
      fcs_near_set_charge_dipole(&near, ifcs_p2nfft_compute_near_charge_dipole);
      fcs_near_set_dipole_dipole(&near, ifcs_p2nfft_compute_near_dipole_dipole);

//       fcs_near_set_field(&near, ifcs_p2nfft_compute_near_field);
//       fcs_near_set_potential(&near, ifcs_p2nfft_compute_near_potential);
//       fcs_near_set_field_potential(&near, ifcs_p2nfft_compute_near_field_and_potential);
    }

    // fcs_int *periodicity = NULL; /* sorter uses periodicity of the communicator */
    fcs_near_set_system(&near, d->box_base, d->box_a, d->box_b, d->box_c, d->periodicity);
  
    fcs_near_set_particles(&near, sorted_num_particles, sorted_num_particles, sorted_positions, sorted_charges, sorted_indices,
        (compute_field)?sorted_field:NULL, (compute_potential)?sorted_potential:NULL);
    fcs_near_set_dipole_particles(&near, sorted_num_dipole_particles, sorted_num_dipole_particles, sorted_dipole_positions, sorted_dipole_moments, sorted_dipole_indices,
        (compute_dipole_field)?sorted_dipole_field:NULL, (compute_dipole_potential)?sorted_dipole_potential:NULL);
  
    fcs_near_set_ghosts(&near, ghost_num_particles, ghost_positions, ghost_charges, ghost_indices);
    fcs_near_set_dipole_ghosts(&near, ghost_num_dipole_particles, ghost_dipole_positions, ghost_dipole_moments, ghost_dipole_indices);
  
    if(d->interpolation_order >= 0 && !compute_dipole_potential && !compute_dipole_field ){
      /* interpolation loop only works for charge-charge interactions */
      /* TODO: implement interpolation for dipole interactions */

      ifcs_p2nfft_near_params near_params;
      near_params.interpolation_order = d->interpolation_order;
      near_params.interpolation_num_nodes = d->near_interpolation_num_nodes;
      near_params.near_interpolation_table_potential = d->near_interpolation_table_potential;
      near_params.near_interpolation_table_force = d->near_interpolation_table_force;
      near_params.one_over_r_cut = d->one_over_r_cut;

      fcs_near_compute(&near, d->r_cut, &near_params, d->cart_comm_3d);
    } else // in all other cases use p2nfft handle for transmission of parameters
      fcs_near_compute(&near, d->r_cut, rd, d->cart_comm_3d);
  
    fcs_near_destroy(&near);
  }

  /* Finish near field timing */
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Near field computation");
  
  /* Checksum: global sum of nearfield energy */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  fcs_float near_energy = 0.0;
  fcs_float near_global;
  if(compute_potential)
    for(fcs_int j = 0; j < sorted_num_particles; ++j)
      near_energy += 0.5 * sorted_charges[j] * sorted_potential[j];
  MPI_Reduce(&near_energy, &near_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: near field energy: %" FCS_LMOD_FLOAT "f\n", near_global);
#endif

  /* Checksum: fields resulting from nearfield interactions */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  for(fcs_int t=0; t<3; t++){
    rsum = 0.0;
    if(compute_field)
      for(fcs_int j = 0; (j < sorted_num_particles); ++j)
        rsum += fabs(sorted_field[3*j+t]);
    MPI_Reduce(&rsum, &rsum_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
    if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: near field %" FCS_LMOD_INT "d. component: %" FCS_LMOD_FLOAT "f\n", t, rsum_global);
  }
  
  if(compute_field)
    if (myrank == 0) fprintf(stderr, "E_NEAR(0) = %" FCS_LMOD_FLOAT "e\n", sorted_field[0]);
#endif
      
  /* Reinit PNFFT nodes: Number and positions of nodes typically change between runs */
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);
  {
    unsigned pnfft_malloc_flags = PNFFT_MALLOC_X | PNFFT_MALLOC_F;
    if(compute_field)     pnfft_malloc_flags |= PNFFT_MALLOC_GRAD_F;
    FCS_PNFFT(free_nodes)(d->charges, PNFFT_FREE_ALL);
    d->charges = FCS_PNFFT(init_nodes)(sorted_num_particles, pnfft_malloc_flags);
  }
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Reinit charges");

  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);
  {
    unsigned pnfft_malloc_flags = PNFFT_MALLOC_X | PNFFT_MALLOC_GRAD_F;
    if(compute_dipole_field)     pnfft_malloc_flags |= PNFFT_MALLOC_HESSIAN_F;
    FCS_PNFFT(free_nodes)(d->dipoles, PNFFT_FREE_ALL);
    d->dipoles = FCS_PNFFT(init_nodes)(sorted_num_dipole_particles, pnfft_malloc_flags);
  }
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Reinit dipoles");

  /* catch data pointers from PNFFT opaque plans */
  fcs_pnfft_complex *f_hat             = FCS_PNFFT(get_f_hat)(d->pnfft);

  fcs_pnfft_complex *charges_f         = FCS_PNFFT(get_f)(d->charges);
  fcs_pnfft_complex *charges_grad_f    = FCS_PNFFT(get_grad_f)(d->charges);
  fcs_float         *charges_x         = FCS_PNFFT(get_x)(d->charges);
  fcs_pnfft_complex *dipoles_grad_f    = FCS_PNFFT(get_grad_f)(d->dipoles);
  fcs_pnfft_complex *dipoles_hessian_f = FCS_PNFFT(get_hessian_f)(d->dipoles);
  fcs_float         *dipoles_x         = FCS_PNFFT(get_x)(d->dipoles);

  /* Set NFFT nodes within [-0.5,0.5]^3 */
  for (fcs_int j = 0; j < sorted_num_particles; ++j)
  {
    fcs_float pos[3];

    pos[0] = sorted_positions[3*j + 0] - d->box_base[0];
    pos[1] = sorted_positions[3*j + 1] - d->box_base[1];
    pos[2] = sorted_positions[3*j + 2] - d->box_base[2];

    charges_x[3 * j + 0] = ( XYZ2TRI(0, pos, d->box_inv) - 0.5 ) / d->box_expand[0];
    charges_x[3 * j + 1] = ( XYZ2TRI(1, pos, d->box_inv) - 0.5 ) / d->box_expand[1];
    charges_x[3 * j + 2] = ( XYZ2TRI(2, pos, d->box_inv) - 0.5 ) / d->box_expand[2];
  }
  for (fcs_int j = 0; j < sorted_num_dipole_particles; ++j)
  {
    fcs_float pos[3];

    pos[0] = sorted_dipole_positions[3*j + 0] - d->box_base[0];
    pos[1] = sorted_dipole_positions[3*j + 1] - d->box_base[1];
    pos[2] = sorted_dipole_positions[3*j + 2] - d->box_base[2];

    dipoles_x[3 * j + 0] = ( XYZ2TRI(0, pos, d->box_inv) - 0.5 ) / d->box_expand[0];
    dipoles_x[3 * j + 1] = ( XYZ2TRI(1, pos, d->box_inv) - 0.5 ) / d->box_expand[1];
    dipoles_x[3 * j + 2] = ( XYZ2TRI(2, pos, d->box_inv) - 0.5 ) / d->box_expand[2];
  }
    
  /* Set NFFT values */
  for (fcs_int j = 0; j < sorted_num_particles; ++j)
    charges_f[j] = sorted_charges[j];
  
//   for (fcs_int j = 0; j < 3*sorted_num_dipole_particles; ++j)
//     dipoles_grad_f[j] = sorted_dipole_moments[j];
  
  for (fcs_int j = 0; j < sorted_num_dipole_particles; ++j){
    dipoles_grad_f[3 * j + 0] = XYZ2TRI(0, sorted_dipole_moments + 3*j, d->box_inv) / d->box_expand[0];
    dipoles_grad_f[3 * j + 1] = XYZ2TRI(1, sorted_dipole_moments + 3*j, d->box_inv) / d->box_expand[1];
    dipoles_grad_f[3 * j + 2] = XYZ2TRI(2, sorted_dipole_moments + 3*j, d->box_inv) / d->box_expand[2];
  }

  /* Reset pnfft timer (delete timings from fcs_init and fcs_tune) */  
#if FCS_ENABLE_INFO && !FCS_P2NFFT_DISABLE_PNFFT_INFO
  FCS_PNFFT(reset_timer)(d->pnfft);
#endif

  /* Checksum: Input of adjoint NFFT */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  rsum = 0.0;
  for (fcs_int j = 0; j < 3*sorted_num_particles; ++j)
    rsum += fabs(charges_x[j]);
  MPI_Reduce(&rsum, &rsum_global, 1, MPI_DOUBLE, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: checksum of x: %" FCS_LMOD_FLOAT "e\n", rsum_global);

  csum = 0.0;
  for (fcs_int j = 0; j < sorted_num_particles; ++j)
    csum += fabs(creal(charges_f[j])) + _Complex_I * fabs(cimag(charges_f[j]));
  MPI_Reduce(&csum, &csum_global, 2, MPI_DOUBLE, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: checksum of NFFT^H input: %e + I* %e\n", creal(csum_global), cimag(csum_global));
#endif

  /* Start far field timing */
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);

  /* Perform adjoint NFFT for charges and dipoles */
  {
    unsigned direct_flag = (d->pnfft_direct) ? PNFFT_COMPUTE_DIRECT : 0;
    FCS_PNFFT(zero_f_hat)(d->pnfft);
    FCS_PNFFT(adj)(d->pnfft, d->charges, direct_flag | PNFFT_COMPUTE_ACCUMULATED | PNFFT_COMPUTE_F);
    FCS_PNFFT(adj)(d->pnfft, d->dipoles, direct_flag | PNFFT_COMPUTE_ACCUMULATED | PNFFT_COMPUTE_GRAD_F);
  }

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

  /* Perform NFFT for charges and dipoles */
  {
    unsigned direct_flag = (d->pnfft_direct) ? PNFFT_COMPUTE_DIRECT : 0;
    unsigned compute_flags_charges = 0;
    if(compute_potential) compute_flags_charges |= PNFFT_COMPUTE_F;
    if(compute_field)     compute_flags_charges |= PNFFT_COMPUTE_GRAD_F;
    FCS_PNFFT(trafo)(d->pnfft, d->charges, direct_flag | compute_flags_charges);

    unsigned compute_flags_dipoles = 0;
    if(compute_dipole_potential) compute_flags_dipoles |= PNFFT_COMPUTE_GRAD_F;
    if(compute_dipole_field)     compute_flags_dipoles |= PNFFT_COMPUTE_HESSIAN_F;
    FCS_PNFFT(trafo)(d->pnfft, d->dipoles, direct_flag | compute_flags_dipoles);
  }

  /* Copy the results to the output vector and rescale with L^{-T} */
  if(compute_potential)
    for (fcs_int j = 0; j < sorted_num_particles; ++j)
      sorted_potential[j] += creal(charges_f[j]);

  /* Rescale all gradients L^{-T} * grad_f */
  if(compute_dipole_potential){

// for (fcs_int j = 0; j < sorted_num_dipole_particles; ++j)
//   for(int t=0; t<3; t++)
//     fprintf(stderr, "near_pot[%d, %d] = %.6e\n", j, t, sorted_dipole_potential[3*j+t]);

    for (fcs_int j = 0; j < sorted_num_dipole_particles; ++j){
      sorted_dipole_potential[3 * j + 0] -= fcs_creal( At_TIMES_VEC(d->ebox_inv, dipoles_grad_f + 3*j, 0) );
      sorted_dipole_potential[3 * j + 1] -= fcs_creal( At_TIMES_VEC(d->ebox_inv, dipoles_grad_f + 3*j, 1) );
      sorted_dipole_potential[3 * j + 2] -= fcs_creal( At_TIMES_VEC(d->ebox_inv, dipoles_grad_f + 3*j, 2) );
    }

  }

  if(compute_field){
    for (fcs_int j = 0; j < sorted_num_particles; ++j){
      sorted_field[3 * j + 0] -= fcs_creal( At_TIMES_VEC(d->ebox_inv, charges_grad_f + 3*j, 0) );
      sorted_field[3 * j + 1] -= fcs_creal( At_TIMES_VEC(d->ebox_inv, charges_grad_f + 3*j, 1) );
      sorted_field[3 * j + 2] -= fcs_creal( At_TIMES_VEC(d->ebox_inv, charges_grad_f + 3*j, 2) );
    }
  }

  /* Rescale all (symmetric) Hessian via L^{-T} * Hf * L^{-1} */
  if(compute_dipole_field){
    for (fcs_int j = 0; j < sorted_num_dipole_particles; ++j){
      sorted_dipole_field[6 * j + 0] -= fcs_creal( At_TIMES_SYMMAT_TIMES_A(d->ebox_inv, dipoles_hessian_f + 6*j, 0, 0) );
      sorted_dipole_field[6 * j + 1] -= fcs_creal( At_TIMES_SYMMAT_TIMES_A(d->ebox_inv, dipoles_hessian_f + 6*j, 0, 1) );
      sorted_dipole_field[6 * j + 2] -= fcs_creal( At_TIMES_SYMMAT_TIMES_A(d->ebox_inv, dipoles_hessian_f + 6*j, 0, 2) );
      sorted_dipole_field[6 * j + 3] -= fcs_creal( At_TIMES_SYMMAT_TIMES_A(d->ebox_inv, dipoles_hessian_f + 6*j, 1, 1) );
      sorted_dipole_field[6 * j + 4] -= fcs_creal( At_TIMES_SYMMAT_TIMES_A(d->ebox_inv, dipoles_hessian_f + 6*j, 1, 2) );
      sorted_dipole_field[6 * j + 5] -= fcs_creal( At_TIMES_SYMMAT_TIMES_A(d->ebox_inv, dipoles_hessian_f + 6*j, 2, 2) );
    }
  }

  /* Checksum: fields resulting from farfield interactions */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  for(fcs_int t=0; t<3; t++){
    rsum = 0.0;
    if(compute_field)
      for(fcs_int j = 0; (j < sorted_num_particles); ++j)
        rsum += fabs(sorted_field[3*j+t]);
    MPI_Reduce(&rsum, &rsum_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
    if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: checksum near plus far field %" FCS_LMOD_INT "d. component: %" FCS_LMOD_FLOAT "f\n", t, rsum_global);
  }

  if(compute_field){
    if (myrank == 0) fprintf(stderr, "E_FAR(0) = %" FCS_LMOD_FLOAT "e\n", -fcs_creal( At_TIMES_VEC(d->ebox_inv, charges_grad_f + 3*0, 0) ));
    if (myrank == 0) fprintf(stderr, "E_NEAR_FAR(0) = %" FCS_LMOD_FLOAT "e\n", sorted_field[0]);
  }
#endif

  /* Finish far field timing */
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Far field computation");

#if FCS_ENABLE_TIMING_PNFFT
  /* Print pnfft timer */
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);
  FCS_PNFFT(print_average_timer_adv)(d->pnfft, d->cart_comm_3d);
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Printing of PNFFT timings");
#endif

  /* Checksum: global sum of farfield energy */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  fcs_float far_energy = 0.0;
  fcs_float far_global;

  for(fcs_int j = 0; j < sorted_num_particles; ++j)
    far_energy += 0.5 * sorted_charges[j] * charges_f[j];

  MPI_Reduce(&far_energy, &far_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: far field energy: %" FCS_LMOD_FLOAT "f\n", far_global);
#endif
  
  /* Start self interaction timing */
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);

  /* Calculate self-interactions */
  if(compute_potential){
    fcs_float self = ifcs_p2nfft_compute_self_potential(rd);
    for (fcs_int j = 0; j < sorted_num_particles; ++j)
      sorted_potential[j] -= sorted_charges[j] * self;
  }

  if(compute_dipole_potential){
    fcs_float self = ifcs_p2nfft_compute_self_dipole_potential(rd);
    for (fcs_int j = 0; j < sorted_num_dipole_particles; ++j){
      sorted_dipole_potential[3 * j + 0] +=  sorted_dipole_moments[3 * j + 0] * self; 
      sorted_dipole_potential[3 * j + 1] +=  sorted_dipole_moments[3 * j + 1] * self;
      sorted_dipole_potential[3 * j + 2] +=  sorted_dipole_moments[3 * j + 2] * self;
    }

  }

  /* Finish self interaction timing */
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "self interaction calculation");

  /* Calculate virial if needed */
  if(d->virial != NULL){
    if (d->num_periodic_dims == 3) {
      fcs_float total_energy = 0.0;
      fcs_float total_global;
      if(compute_potential)
        for(fcs_int j = 0; j < sorted_num_particles; ++j)
          total_energy += 0.5 * sorted_charges[j] * sorted_potential[j];

      MPI_Allreduce(&total_energy, &total_global, 1, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);

      /* Approximate virial in 3d-periodic case:
       * Fill the main diagonal with one third of the total energy */      
      for(fcs_int t=0; t<9; t++)
        d->virial[t] = 0.0;
      d->virial[0] = d->virial[4] = d->virial[8] = total_global/3.0;
    } 
    else if (d->num_periodic_dims == 0) {
      fcs_float local_virial[9];

      for(fcs_int t=0; t<9; t++)
        local_virial[t] = 0.0;
      
      if(compute_field)
        for(fcs_int j = 0; j < sorted_num_particles; ++j)
          for(fcs_int t0=0; t0<3; t0++)
            for(fcs_int t1=0; t1<3; t1++)
              local_virial[t0*3+t1] += sorted_charges[j] * sorted_field[3*j+t0] * sorted_positions[3*j+t1];
      
      MPI_Allreduce(local_virial, d->virial, 9, FCS_MPI_FLOAT, MPI_SUM, d->cart_comm_3d);
    }
    else {
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Virial computation is currently not available for mixed boundary conditions.");
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
  if(compute_potential)
    for(fcs_int j = 0; j < sorted_num_particles; ++j)
      total_energy += 0.5 * sorted_charges[j] * sorted_potential[j];

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
  FCS_P2NFFT_START_TIMING(d->cart_comm_3d);

  fcs_gridsort_set_sorted_results(&gridsort, sorted_num_particles, sorted_field, sorted_potential);
  fcs_gridsort_set_results(&gridsort, max_local_num_particles, field, potential);

  fcs_gridsort_set_sorted_dipole_results(&gridsort, sorted_num_dipole_particles, sorted_dipole_field, sorted_dipole_potential);
  fcs_gridsort_set_dipole_results(&gridsort, max_local_num_dipole_particles, dipole_field, dipole_potential);

  fcs_int resort;

  if (d->resort) resort = fcs_gridsort_prepare_resort(&gridsort, d->cart_comm_3d);
  else resort = 0;

  /* Backsort data into user given ordering (if resort is disabled) */
  if (!resort) fcs_gridsort_sort_backward(&gridsort, d->cart_comm_3d);

  fcs_gridsort_resort_destroy(&d->gridsort_resort);

  if (resort) fcs_gridsort_resort_create(&d->gridsort_resort, &gridsort, d->cart_comm_3d);
  
  d->local_num_particles = local_num_particles;
  d->local_num_dipole_particles = local_num_dipole_particles;

  if (sorted_field) free(sorted_field);
  if (sorted_potential) free(sorted_potential);

  if (sorted_dipole_field) free(sorted_dipole_field);
  if (sorted_dipole_potential) free(sorted_dipole_potential);

  fcs_gridsort_free(&gridsort);

  fcs_gridsort_destroy(&gridsort);

  /* Finish back sort timing */
  FCS_P2NFFT_FINISH_TIMING(d->cart_comm_3d, "Backward sort");

#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  /* print first value of fields */
  if(compute_field){
    if (myrank == 0) printf("P2NFFT_INFO: E(0) = %e\n", creal(field[0]));
    if (myrank == 0) printf("P2NFFT_INFO: dE(0) = %e\n", creal(field[0])+1.619834832399799e-06);
  }
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

