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

#include "run_0dp_noncubic.h"
#include "types.h"
#include "utils.h"
#include "nearfield.h"
#include <common/near/near.h>
#include <math.h>
#include <fcs.h>
//#include "constants.h"

#define WORKAROUND_GRIDSORT_BUG 1
#define FCS_P2NFFT_DISABLE_PNFFT_INFO 0

static void convolution(
    const INT *local_N, const fcs_pnfft_complex *regkern_hat,
    fcs_pnfft_complex *f_hat);


FCSResult ifcs_p2nfft_run_0dp_noncubic(
    void *rd, fcs_int local_num_particles, fcs_int max_local_num_particles,
    fcs_float *positions, fcs_float *charges,
    fcs_float *potentials, fcs_float *field
    )
{
//   const char* fnc_name = "ifcs_p2nfft_run";
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct*) rd;
 
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  C csum;
  C csum_global;
  fcs_float rsum;
  fcs_float rsum_global;
#endif

#if FCS_ENABLE_INFO || FCS_P2NFFT_DEBUG || FCS_ENABLE_DEBUG
  int myrank;
  MPI_Comm_rank(d->cart_comm_3d, &myrank);
#endif

  /* Tuning was called in fcs_run */
//  result = ifcs_p2nfft_tune(rd, local_num_particles, positions, charges, d->box_size);

 /* Compute the near field */
  fcs_float box_a[3] = { d->box_l[0], 0, 0 };
  fcs_float box_b[3] = { 0, d->box_l[1], 0 };
  fcs_float box_c[3] = { 0, 0, d->box_l[2] };
  fcs_float box_base[3] = { 0, 0, 0 };
  fcs_near_t near;

  fcs_int sorted_num_particles, ghost_num_particles;
  fcs_float *sorted_positions, *ghost_positions;
  fcs_float *sorted_charges, *ghost_charges;
  fcs_gridsort_index_t *sorted_indices, *ghost_indices;
  fcs_gridsort_t gridsort;

  fcs_gridsort_create(&gridsort);
  
  fcs_gridsort_set_system(&gridsort, box_base, box_a, box_b, box_c, d->periodicity);
  
  fcs_gridsort_set_bounds(&gridsort, d->lower_border, d->upper_border);

  fcs_gridsort_set_particles(&gridsort, local_num_particles, max_local_num_particles, positions, charges);

  fcs_gridsort_sort_forward(&gridsort, (d->short_range_flag ? d->r_cut: 0.0), d->cart_comm_3d);

  fcs_gridsort_separate_ghosts(&gridsort, NULL, NULL);

  fcs_gridsort_get_real_particles(&gridsort, &sorted_num_particles, &sorted_positions, &sorted_charges, &sorted_indices);
  fcs_gridsort_get_ghost_particles(&gridsort, &ghost_num_particles, &ghost_positions, &ghost_charges, &ghost_indices);

  /* Handle particles, that left the box [0,L] */
  /* For periodic boundary conditions: just fold them back  */
  for(fcs_int j=0; j<sorted_num_particles; j++){
    for(fcs_int t=0; t<3; t++){
      while(sorted_positions[3*j+t] < 0)
        sorted_positions[3*j+t] += d->box_l[t];
      while(sorted_positions[3*j+t] >= d->box_l[t])
        sorted_positions[3*j+t] -= d->box_l[t];
    }
  }

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
  
    fcs_near_set_field(&near, ifcs_p2nfft_compute_near_field);
    fcs_near_set_potential(&near, ifcs_p2nfft_compute_near_potential);

    // fcs_int periodicity[3] = { 1, 1, 1 };
    // fcs_int periodicity[3] = { 0, 0, 0 };
    // fcs_int *periodicity = NULL; /* sorter uses periodicity of the communicator */
    fcs_near_set_system(&near, box_base, box_a, box_b, box_c, d->periodicity);
  
    fcs_near_set_particles(&near, sorted_num_particles, sorted_num_particles, sorted_positions, sorted_charges, sorted_indices,
        (compute_field)?sorted_field:NULL, (compute_potentials)?sorted_potentials:NULL);
  
    fcs_near_set_ghosts(&near, ghost_num_particles, ghost_positions, ghost_charges, ghost_indices);
  
    fcs_near_compute(&near, d->r_cut, rd, d->cart_comm_3d);
  
    fcs_near_destroy(&near);

  }
  
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
  
  /* Set the NFFT nodes and values */
  for (fcs_int j = 0; j < sorted_num_particles; ++j)
  {
    x[3 * j + 0] = (sorted_positions[3 * j + 0] - d->box_shifts[0]) / d->box_scales[0];
    x[3 * j + 1] = (sorted_positions[3 * j + 1] - d->box_shifts[1]) / d->box_scales[1];
    x[3 * j + 2] = (sorted_positions[3 * j + 2] - d->box_shifts[2]) / d->box_scales[2];
    
    f[j] = sorted_charges[j];
  }

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

  /* Copy the results to the output vector and rescale */
  for (fcs_int j = 0; j < sorted_num_particles; ++j){
    fcs_float vol = d->box_l[0]*d->box_l[1]*d->box_l[2];
    sorted_potentials[j] += creal(f[j]) / (FCS_P2NFFT_PI * vol);
    for(fcs_int t=0; t<3; t++)
      sorted_field[3 * j + t] -= creal(grad_f[3 * j + t]) / (FCS_P2NFFT_PI * vol * d->box_l[t]);
  }

  /* Calculate self-interactions */
  for (fcs_int j = 0; j < sorted_num_particles; ++j)
    sorted_potentials[j] -= sorted_charges[j] * ifcs_p2nfft_compute_near_potential(rd, 0.0);

  /* Checksum: global sum of farfield energy */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  fcs_float far_energy = 0.0;
  fcs_float far_global;
  for(fcs_int j = 0; j < sorted_num_particles; ++j)
    if(d->use_ewald)
      far_energy += 0.5 * sorted_charges[j] * f[j] / (FCS_P2NFFT_PI * d->box_scale);
    else
      far_energy += 0.5 * sorted_charges[j] * f[j] / d->box_scale;

  MPI_Reduce(&far_energy, &far_global, 1, FCS_MPI_FLOAT, MPI_SUM, 0, d->cart_comm_3d);
  if (myrank == 0) fprintf(stderr, "P2NFFT_DEBUG: far field energy: %" FCS_LMOD_FLOAT "f\n", far_global);
#endif

  /* Calculate virial if needed */
  if(d->virial != NULL){
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

  /* Checksum: global sum of self energy */
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  fcs_float self_energy = 0.0;
  fcs_float self_global;
  for(fcs_int j = 0; j < sorted_num_particles; ++j)
    self_energy -= 0.5 * sorted_charges[j] * sorted_charges[j] * ifcs_p2nfft_compute_near_potential(rd, 0.0);

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

  /* Backsort data into user given ordering */
  fcs_int set_values = 1; /* set (1) or add (0) the field and potentials */

  fcs_gridsort_sort_backward(&gridsort,
      sorted_field, sorted_potentials,
      field, potentials, set_values, d->cart_comm_3d);

  if (sorted_field) free(sorted_field);
  if (sorted_potentials) free(sorted_potentials);

  fcs_gridsort_free(&gridsort);

  fcs_gridsort_destroy(&gridsort);

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
      for (INT k2 = 0; k2 < local_N[2]; ++k2, ++m)
        f_hat[m] *= regkern_hat[m];
}




