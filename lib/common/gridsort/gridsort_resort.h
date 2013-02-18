/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann
  
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


#ifndef __GRIDSORT_RESORT_H__
#define __GRIDSORT_RESORT_H__


#ifdef __cplusplus
extern "C" {
#endif


#include <mpi.h>

#include "gridsort.h"


/**
 * @brief gridsort_resort object structure
 */
typedef struct _fcs_gridsort_resort_t
{
  fcs_int noriginal_particles, nsorted_particles;
  fcs_gridsort_index_t *indices;

} *fcs_gridsort_resort_t;

#define FCS_GRIDSORT_RESORT_NULL  NULL


/**
 * @brief create gridsort_resort object from given gridsort object
 * @param gsr fcs_gridsort_resort_t* gridsort_resort object
 * @param gs fcs_gridsort_t* gridsort object
 * @param comm MPI_Comm communicator used for previous gridsort sorting
 */
void fcs_gridsort_resort_create(fcs_gridsort_resort_t *gsr, fcs_gridsort_t *gs, MPI_Comm comm);

/**
 * @brief destroy gridsort_resort object
 * @param gsr fcs_gridsort_resort_t* gridsort_resort object
 */
void fcs_gridsort_resort_destroy(fcs_gridsort_resort_t *gsr);

/**
 * @brief print resort information stored in gridsort_resort object (for DEBUG)
 * @param gsr fcs_gridsort_resort_t* gridsort_resort object
 */
void fcs_gridsort_resort_print(fcs_gridsort_resort_t gsr, MPI_Comm comm);

/**
 * @brief return resorting availability, i.e., return 1 if resorting can be performed, otherwise 0
 * @param gsr fcs_gridsort_resort_t* gridsort_resort object
 */
fcs_int fcs_gridsort_resort_is_available(fcs_gridsort_resort_t gsr);

/**
 * @brief return number of original (unsorted, input) particles belonging to the resorting
 * @param gsr fcs_gridsort_resort_t* gridsort_resort object
 */
fcs_int fcs_gridsort_resort_get_original_particles(fcs_gridsort_resort_t gsr);

/**
 * @brief return number of sorted (output) particles belonging to the resorting
 * @param gsr fcs_gridsort_resort_t* gridsort_resort object
 */
fcs_int fcs_gridsort_resort_get_sorted_particles(fcs_gridsort_resort_t gsr);

/**
 * @brief perform resorting of integer values (i.e., fcs_int)
 * @param gsr fcs_gridsort_resort_t* gridsort_resort object
 * @param src fcs_int* array of integer values in original (unsorted, input) order
 * @param dst fcs_int* array to store resorted integer values
 * @param n fcs_int number of integer values to resort for each particle
 * @param comm MPI_Comm communicator to be used for resorting
 */
void fcs_gridsort_resort_ints(fcs_gridsort_resort_t gsr, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm);

/**
 * @brief perform resorting of float values (i.e., fcs_float)
 * @param gsr fcs_gridsort_resort_t* gridsort_resort object
 * @param src fcs_float* array of float values in original (unsorted, input) order
 * @param dst fcs_float* array to store resorted float values
 * @param n fcs_int number of float values to resort for each particle
 * @param comm MPI_Comm communicator to be used for resorting
 */
void fcs_gridsort_resort_floats(fcs_gridsort_resort_t gsr, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm);

/**
 * @brief perform resorting of byte values
 * @param gsr fcs_gridsort_resort_t* gridsort_resort object
 * @param src void* array of byte values in original (unsorted, input) order
 * @param dst void* array to store resorted byte values
 * @param n fcs_int number of byte values to resort for each particle
 * @param comm MPI_Comm communicator to be used for resorting
 */
void fcs_gridsort_resort_bytes(fcs_gridsort_resort_t gsr, void *src, void *dst, fcs_int n, MPI_Comm comm);


#ifdef __cplusplus
}
#endif


#endif /* __GRIDSORT_RESORT_H__ */
