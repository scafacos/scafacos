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


#ifndef __REDIST_H__
#define __REDIST_H__


#ifdef __cplusplus
extern "C" {
#endif


#include <mpi.h>

#include "redist_index.h"


/**
 * @brief redist object structure
 */
typedef struct _fcs_redist_t
{
  MPI_Comm comm;

  fcs_int noriginal_particles, max_noriginal_particles;
  fcs_float *original_positions, *original_charges;
  fcs_float *original_field, *original_potentials;

  fcs_int nredistributed_particles, max_nredistributed_particles;
  fcs_redist_index_t *redistributed_indices;
  fcs_float *redistributed_positions, *redistributed_charges;
  fcs_float *redistributed_field, *redistributed_potentials;

} *fcs_redist_t;

#define FCS_REDIST_NULL  NULL


/**
 * @brief create redist object from given gridsort object
 * @param redist fcs_redist_t* redistribution object
 * @param comm MPI_Comm MPI communicator
 */
void fcs_redist_create(fcs_redist_t *redist, MPI_Comm comm);

/**
 * @brief destroy redist object
 * @param redist fcs_redist_t* redistribution object
 */
void fcs_redist_destroy(fcs_redist_t *redist);

/**
 * @brief print redist object information
 * @param redist fcs_redist_t redistribution object
 */
void fcs_redist_print(fcs_redist_t redist);

/**
 * @brief set the original particles (before redistribution)
 * @param redist fcs_redist_t redistribution object
 */
void fcs_redist_set_original_particles(fcs_redist_t redist, fcs_int noriginal_particles, fcs_int max_noriginal_particles, fcs_float *original_positions, fcs_float *original_charges, fcs_float *original_field, fcs_float *original_potentials);

/**
 * @brief return the redistributed particles (after redistribution)
 * @param redist fcs_redist_t redistribution object
 */
void fcs_redist_get_redistributed_particles(fcs_redist_t redist, fcs_int *nredistributed_particles, fcs_int *max_nredistributed_particles, fcs_float **redistributed_positions, fcs_float **redistributed_charges, fcs_float **redistributed_field, fcs_float **redistributed_potentials);

/**
 * @brief redistribute particles equally among all processes (i.e. positions and charges)
 * @param redist fcs_redist_t redistribution object
 * @return exit code
 */
fcs_int fcs_redist_redistribute_forward_equal(fcs_redist_t redist);

/**
 * @brief undo redistribution of particles (i.e. field and potential results)
 * @param redist fcs_redist_t redistribution object
 * @return exit code
 */
fcs_int fcs_redist_redistribute_backward(fcs_redist_t redist);


#ifdef __cplusplus
}
#endif


#endif /* __REDIST_H__ */
