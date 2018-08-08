/*
  Copyright (C) 2018 Michael Hofmann
  
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


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "../gridsort/include/sl_forw.h"

#include "z_tools.h"
#include "redist_common.h"
#include "redist_index.h"
#include "redist.h"


#define REDISTRIBUTE_ENABLED  1
#define PRINT_INDICES  0


void fcs_redist_create(fcs_redist_t *redist, MPI_Comm comm)
{
  *redist = malloc(sizeof(**redist));

  (*redist)->comm = comm;

  (*redist)->noriginal_particles = -1;
  (*redist)->max_noriginal_particles = -1;
  (*redist)->original_positions = NULL;
  (*redist)->original_charges = NULL;
  (*redist)->original_field = NULL;
  (*redist)->original_potentials = NULL;

  (*redist)->nredistributed_particles = -1;
  (*redist)->max_nredistributed_particles = -1;
  (*redist)->redistributed_indices = NULL;
  (*redist)->redistributed_positions = NULL;
  (*redist)->redistributed_charges = NULL;
  (*redist)->redistributed_field = NULL;
  (*redist)->redistributed_potentials = NULL;
}


void fcs_redist_destroy(fcs_redist_t *redist)
{
  if (*redist == FCS_REDIST_NULL) return;

  if ((*redist)->redistributed_indices != NULL) fcs_redist_indices_free((*redist)->redistributed_indices);
  if ((*redist)->redistributed_positions != NULL) free((*redist)->redistributed_positions);
  if ((*redist)->redistributed_charges != NULL) free((*redist)->redistributed_charges);
  if ((*redist)->redistributed_field != NULL) free((*redist)->redistributed_field);
  if ((*redist)->redistributed_potentials != NULL) free((*redist)->redistributed_potentials);

  (*redist)->redistributed_indices = NULL;
  (*redist)->redistributed_positions = NULL;
  (*redist)->redistributed_charges = NULL;
  (*redist)->redistributed_field = NULL;
  (*redist)->redistributed_potentials = NULL;

  free(*redist);

  *redist = FCS_REDIST_NULL;
}


void fcs_redist_print(fcs_redist_t redist)
{
  printf("fcs_redist: %p -> [%" FCS_LMOD_INT "d, %" FCS_LMOD_INT "d, %p]\n", redist, redist->noriginal_particles, redist->nredistributed_particles, redist->redistributed_indices);

  if (redist->redistributed_indices)
    fcs_redist_indices_print(redist->nredistributed_particles, redist->redistributed_indices);
  else
    printf("-> ID\n");
}


void fcs_redist_set_original_particles(fcs_redist_t redist, fcs_int noriginal_particles, fcs_int max_noriginal_particles, fcs_float *original_positions, fcs_float *original_charges, fcs_float *original_field, fcs_float *original_potentials)
{
  redist->noriginal_particles = noriginal_particles;
  redist->max_noriginal_particles = max_noriginal_particles;
  redist->original_positions = original_positions;
  redist->original_charges = original_charges;
  redist->original_field = original_field;
  redist->original_potentials = original_potentials;
}


void fcs_redist_get_redistributed_particles(fcs_redist_t redist, fcs_int *nredistributed_particles, fcs_int *max_nredistributed_particles, fcs_float **redistributed_positions, fcs_float **redistributed_charges, fcs_float **redistributed_field, fcs_float **redistributed_potentials)
{
  *nredistributed_particles = redist->nredistributed_particles;
  *max_nredistributed_particles = redist->max_nredistributed_particles;
  *redistributed_positions = redist->redistributed_positions;
  *redistributed_charges = redist->redistributed_charges;
  if (redistributed_field) *redistributed_field = redist->redistributed_field;
  if (redistributed_potentials) *redistributed_potentials = redist->redistributed_potentials;
}


fcs_int fcs_redist_redistribute_forward_equal(fcs_redist_t redist)
{
  INFO_CMD(printf(INFO_PRINT_PREFIX "forward: original particles: %" FCS_LMOD_INT "d / %" FCS_LMOD_INT "d\n", redist->noriginal_particles, redist->max_noriginal_particles););

#if REDISTRIBUTE_ENABLED

  int comm_size, comm_rank;
  fcs_redist_index_t *original_indices;
  fcs_int total_nparticles;
  fcs_forw_slint_t dst_size;
  fcs_forw_elements_t sin0, sout0;


  MPI_Comm_rank(redist->comm, &comm_rank);
  MPI_Comm_size(redist->comm, &comm_size);

  MPI_Allreduce(&redist->noriginal_particles, &total_nparticles, 1, FCS_MPI_INT, MPI_SUM, redist->comm);

  redist->nredistributed_particles = redist->max_nredistributed_particles = (total_nparticles / comm_size) + ((comm_rank < total_nparticles % comm_size)?1:0);

  redist->redistributed_indices    = fcs_redist_indices_alloc(redist->max_nredistributed_particles);
  redist->redistributed_positions  = malloc(3 * redist->max_nredistributed_particles * sizeof(fcs_float));
  redist->redistributed_charges    = malloc(redist->max_nredistributed_particles * sizeof(fcs_float));
  redist->redistributed_field      = (redist->original_field)?(malloc(3 * redist->max_nredistributed_particles * sizeof(fcs_float))):NULL;
  redist->redistributed_potentials = (redist->original_potentials)?(malloc(redist->max_nredistributed_particles * sizeof(fcs_float))):NULL;

  fcs_forw_SL_DEFCON(mpi.rank) = comm_rank;

  fcs_forw_mpi_datatypes_init();

  original_indices = fcs_redist_indices_alloc(redist->noriginal_particles);
  fcs_redist_indices_init(redist->noriginal_particles, original_indices, comm_rank);

#if PRINT_INDICES
  fcs_redist_indices_print(redist->noriginal_particles, original_indices);
#endif

  fcs_forw_elem_set_size(&sin0, redist->noriginal_particles);
  fcs_forw_elem_set_max_size(&sin0, redist->max_noriginal_particles);
  fcs_forw_elem_set_keys(&sin0, original_indices);
  fcs_forw_elem_set_data(&sin0, redist->original_positions, redist->original_charges);

  fcs_forw_elem_set_size(&sout0, redist->nredistributed_particles);
  fcs_forw_elem_set_max_size(&sout0, redist->max_nredistributed_particles);
  fcs_forw_elem_set_keys(&sout0, redist->redistributed_indices);
  fcs_forw_elem_set_data(&sout0, redist->redistributed_positions, redist->redistributed_charges);

  dst_size = redist->nredistributed_particles;

  fcs_forw_mpi_rebalance(&sin0, &sout0, 0, &dst_size, comm_size, comm_rank, redist->comm);

#if PRINT_INDICES
  fcs_redist_indices_print(redist->nredistributed_particles, redist->redistributed_indices);
#endif

  fcs_redist_indices_free(original_indices);

  fcs_forw_mpi_datatypes_release();

#else /* REDISTRIBUTE_ENABLED */

  redist->nredistributed_particles = redist->noriginal_particles;
  redist->max_nredistributed_particles = redist->max_noriginal_particles;
  redist->redistributed_positions = redist->original_positions;
  redist->redistributed_charges = redist->original_charges;
  redist->redistributed_field = redist->original_field;
  redist->redistributed_potentials = redist->original_potentials;

#endif /* REDISTRIBUTE_ENABLED */

  INFO_CMD(printf(INFO_PRINT_PREFIX "forward: redistributed particles: %" FCS_LMOD_INT "d / %" FCS_LMOD_INT "d\n", redist->nredistributed_particles, redist->max_nredistributed_particles););

  return 1;
}


fcs_int fcs_redist_redistribute_backward(fcs_redist_t redist)
{
  INFO_CMD(printf(INFO_PRINT_PREFIX "backward: redistributed particles: %" FCS_LMOD_INT "d / %" FCS_LMOD_INT "d\n", redist->nredistributed_particles, redist->max_nredistributed_particles););

#if REDISTRIBUTE_ENABLED

  fcs_redist_indices_sort_back_results(redist->nredistributed_particles, redist->redistributed_indices, redist->redistributed_field, redist->redistributed_potentials, redist->noriginal_particles, redist->original_field, redist->original_potentials, redist->comm);

  fcs_redist_indices_free(redist->redistributed_indices);
  free(redist->redistributed_positions);
  free(redist->redistributed_charges);
  if (redist->redistributed_field) free(redist->redistributed_field);
  if (redist->redistributed_potentials) free(redist->redistributed_potentials);

  redist->redistributed_indices = NULL;
  redist->redistributed_positions = NULL;
  redist->redistributed_charges = NULL;
  redist->redistributed_field = NULL;
  redist->redistributed_potentials = NULL;

#endif /* REDISTRIBUTE_ENABLED */

  INFO_CMD(printf(INFO_PRINT_PREFIX "backward: original particles: %" FCS_LMOD_INT "d / %" FCS_LMOD_INT "d\n", redist->noriginal_particles, redist->max_noriginal_particles););

  return 1;
}
