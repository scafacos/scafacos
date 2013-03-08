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


#ifndef __RESORT_H__
#define __RESORT_H__


#ifdef __cplusplus
extern "C" {
#endif


#include <mpi.h>


typedef long long fcs_resort_index_t;
#define FCS_MPI_RESORT_INDEX  MPI_LONG_LONG


#define FCS_RESORT_INDEX_IS_VALID(_x_)       ((_x_) >= 0)
#define FCS_RESORT_INDEX_VAL_PROC(_proc_)    (((fcs_resort_index_t) (_proc_)) << 32)
#define FCS_RESORT_INDEX_VAL_POS(_pos_)      ((fcs_resort_index_t) (_pos_))
#define FCS_RESORT_INDEX_VAL(_proc_, _pos_)  (FCS_RESORT_INDEX_VAL_PROC(_proc_) + FCS_RESORT_INDEX_VAL_POS(_pos_))
#define FCS_RESORT_INDEX_GET_PROC(_x_)       ((_x_) >> 32)
#define FCS_RESORT_INDEX_GET_POS(_x_)        ((_x_) & 0x00000000FFFFFFFFLL)

#define FCS_RESORT_INDEX_STR         "(%lld,%lld)"
#define FCS_RESORT_INDEX_PARAM(_x_)  (FCS_RESORT_INDEX_IS_VALID(_x_)?FCS_RESORT_INDEX_GET_PROC(_x_):-1), (FCS_RESORT_INDEX_IS_VALID(_x_)?FCS_RESORT_INDEX_GET_POS(_x_):-1)


fcs_resort_index_t *fcs_resort_indices_alloc(fcs_int nindices);
void fcs_resort_indices_free(fcs_resort_index_t *indices);
void fcs_resort_indices_init(fcs_int nindices, fcs_resort_index_t *indices, int rank);
void fcs_resort_indices_print(fcs_int nindices, fcs_resort_index_t *indices);


/**
 * @brief resort object structure
 */
typedef struct _fcs_resort_t
{
  fcs_int noriginal_particles, nsorted_particles;
  fcs_resort_index_t *indices;
  
  fcs_int nprocs;
  int *procs;

} *fcs_resort_t;

#define FCS_RESORT_NULL  NULL


/**
 * @brief create resort object from given gridsort object
 * @param resort fcs_resort_t* resort object
 * @param comm MPI_Comm MPI communicator
 */
void fcs_resort_create(fcs_resort_t *resort);

/**
 * @brief destroy resort object
 * @param resort fcs_resort_t* resort object
 */
void fcs_resort_destroy(fcs_resort_t *resort);

/**
 * @brief print resort object information
 * @param resort fcs_resort_t resort object
 */
void fcs_resort_print(fcs_resort_t resort);

/**
 * @brief return resorting availability, i.e., return 1 if resorting can be performed, otherwise 0
 * @param resort fcs_resort_t resort object
 */
fcs_int fcs_resort_is_available(fcs_resort_t resort);

/**
 * @brief set the number of original (unsorted, input) particles
 * @param resort fcs_resort_t resort object
 */
void fcs_resort_set_original_particles(fcs_resort_t resort, fcs_int original_particles);

/**
 * @brief return number of original (unsorted, input) particles belonging to the resorting
 * @param resort fcs_resort_t resort object
 */
fcs_int fcs_resort_get_original_particles(fcs_resort_t resort);

/**
 * @brief set the number of sorted (output) particles
 * @param resort fcs_resort_t resort object
 */
void fcs_resort_set_sorted_particles(fcs_resort_t resort, fcs_int sorted_particles);

/**
 * @brief return number of sorted (output) particles belonging to the resorting
 * @param resort fcs_resort_t resort object
 */
fcs_int fcs_resort_get_sorted_particles(fcs_resort_t resort);

/**
 * @brief allocate memory for resort indices
 * @param resort fcs_resort_t resort object
 */
void fcs_resort_alloc_indices(fcs_resort_t resort);

/**
 * @brief free allocated memory for resort indices
 * @param resort fcs_resort_t resort object
 */
void fcs_resort_free_indices(fcs_resort_t resort);

/**
 * @brief return pointer to allocated memory for resort indices
 * @param resort fcs_resort_t resort object
 */
fcs_resort_index_t *fcs_resort_get_indices(fcs_resort_t resort);

/**
 * @brief set (optionally) the list of processes that the local process will communicate with for resort
 * @param nprocs fcs_int number of processes in the list
 * @param procs *int list of process ranks
 * @param resort fcs_resort_t resort object
 */
void fcs_resort_set_proclists(fcs_resort_t resort, fcs_int nprocs, int *procs);

/**
 * @brief perform resorting of integer values (i.e., fcs_int)
 * @param resort fcs_resort_t resort object
 * @param src fcs_int* array of integer values in original (unsorted, input) order
 * @param dst fcs_int* array to store resorted integer values
 * @param n fcs_int number of integer values to resort for each particle
 * @param comm MPI_Comm MPI communicator
 */
void fcs_resort_resort_ints(fcs_resort_t resort, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm);

/**
 * @brief perform resorting of float values (i.e., fcs_float)
 * @param resort fcs_resort_t resort object
 * @param src fcs_float* array of float values in original (unsorted, input) order
 * @param dst fcs_float* array to store resorted float values
 * @param n fcs_int number of float values to resort for each particle
 * @param comm MPI_Comm MPI communicator
 */
void fcs_resort_resort_floats(fcs_resort_t resort, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm);

/**
 * @brief perform resorting of byte values
 * @param resort fcs_resort_t resort object
 * @param src void* array of byte values in original (unsorted, input) order
 * @param dst void* array to store resorted byte values
 * @param n fcs_int number of byte values to resort for each particle
 * @param comm MPI_Comm MPI communicator
 */
void fcs_resort_resort_bytes(fcs_resort_t resort, void *src, void *dst, fcs_int n, MPI_Comm comm);


#ifdef __cplusplus
}
#endif


#endif /* __RESORT_H__ */
