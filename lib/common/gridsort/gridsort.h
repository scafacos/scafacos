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


#ifndef __GRIDSORT_H__
#define __GRIDSORT_H__


#ifdef __cplusplus
extern "C" {
#endif


#include <mpi.h>


/**
 * @brief gridsort index type
 */
typedef long long fcs_gridsort_index_t;
#define FCS_MPI_GRIDSORT_INDEX  MPI_LONG_LONG

/**
 * @brief gridsort structure form permanent data
 */
typedef struct _fcs_gridsort_data_t
{
  fcs_float box_base[3], box_a[3], box_b[3], box_c[3];
  fcs_int periodicity[3];

  fcs_float lower_bounds[3], upper_bounds[3];

  fcs_float *bounds;

} *fcs_gridsort_cache_t;

#define FCS_GRIDSORT_CACHE_NULL  NULL


/**
 * @brief gridsort object structure
 */
typedef struct _fcs_gridsort_t
{
  struct _fcs_gridsort_data_t d;

  fcs_int local_nzslices, ghost_nzslices, max_ghost_nzslices;

  fcs_int minalloc;
  fcs_float overalloc;

  fcs_int noriginal_particles, max_noriginal_particles;
  fcs_float *original_positions, *original_charges;

  fcs_int nsorted_particles, max_nsorted_particles;
  fcs_float *sorted_positions, *sorted_charges;
  fcs_gridsort_index_t *sorted_indices;

  fcs_int nsorted_real_particles, nsorted_ghost_particles;

  fcs_int max_nsorted_results;
  fcs_float *sorted_field, *sorted_potentials;

  fcs_int max_noriginal_results;
  fcs_float *original_field, *original_potentials;

  fcs_int nresort_particles;

  fcs_float max_particle_move;
  fcs_int nprocs;
  int *procs;

  fcs_float sub_box_base[3], sub_box_a[3], sub_box_b[3], sub_box_c[3];

  fcs_gridsort_cache_t *cache;

} fcs_gridsort_t;


#define GRIDSORT_GHOST_BASE              0x8000000000000000LL

#define GRIDSORT_IS_GHOST(_x_)           ((_x_) < 0)

#define GRIDSORT_PERIODIC_BITS           10
#define GRIDSORT_PERIODIC_CONST(_v_)     (_v_##LL)
#define GRIDSORT_PERIODIC_MASK           ((GRIDSORT_PERIODIC_CONST(1) << GRIDSORT_PERIODIC_BITS) - GRIDSORT_PERIODIC_CONST(1))
#define GRIDSORT_PERIODIC_SET(_v_, _x_)  (((fcs_gridsort_index_t) (_v_)) << ((_x_) * GRIDSORT_PERIODIC_BITS))
#define GRIDSORT_PERIODIC_GET(_v_, _x_)  (((_v_) >> ((_x_) * GRIDSORT_PERIODIC_BITS)) & GRIDSORT_PERIODIC_MASK)


#include "gridsort_resort.h"


/**
 * @brief create gridsort object
 * @param gs fcs_gridsort_t* gridsort object
 */
void fcs_gridsort_create(fcs_gridsort_t *gs);

/**
 * @brief destroy gridsort object
 * @param gs fcs_gridsort_t* gridsort object
 */
void fcs_gridsort_destroy(fcs_gridsort_t *gs);

/**
 * @brief set particle system properties
 * @param gs fcs_gridsort_t* gridsort object
 * @param box_base fcs_float* origin of the system box
 * @param box_a fcs_float* 1st base vector of the system box
 * @param box_b fcs_float* 2nd base vector of the system box
 * @param box_c fcs_float* 3rd base vector of the system box
 * @param periodicity fcs_int* periodicity of the system, if NULL then periodicities of the (cartesian) communicator are used
 */
void fcs_gridsort_set_system(fcs_gridsort_t *gs, const fcs_float *box_base, const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c, const fcs_int *periodicity);

/**
 * @brief set lower and upper bounds of the local process within in the grid
 * @param gs fcs_gridsort_t* gridsort object
 * @param lower_bounds fcs_float* lower bound of the local process within the grid in each dimension
 * @param upper_bounds fcs_float* upper bound of the local process within the grid in each dimension
 */
void fcs_gridsort_set_bounds(fcs_gridsort_t *gs, fcs_float *lower_bounds, fcs_float *upper_bounds);

/**
 * @brief set number of zslices per process and the number of ghost zslices
 * @param gs fcs_gridsort_t* gridsort object
 * @param local_nzslices fcs_int local number of zslices on each process
 * @param ghost_nzslices fcs_int number of zslices used to create ghost particles
 */
void fcs_gridsort_set_zslices(fcs_gridsort_t *gs, fcs_int local_nzslices, fcs_int ghost_nzslices);

/**
 * @brief set minimum size of memory to be allocated for sorted particle data array
 * @param gs fcs_gridsort_t* gridsort object
 * @param minalloc fcs_int minimum size of memory, default: minalloc = 0 (i.e., no minimum)
 * @param ghost_nzslices fcs_int number of zslices used to create ghost particles
 */
void fcs_gridsort_set_minalloc(fcs_gridsort_t *gs, fcs_int minalloc);

/**
 * @brief set size of extra memory to be allocated for sorted particle data array
 * @param gs fcs_gridsort_t* gridsort object
 * @param overalloc fcs_float if overalloc >= 0 then newly allocated arrays for sorted particle data can store overalloc more particles than actually required,
 *   if overalloc < 0 then newly allocated arrays become a factor of 1-overalloc as large (e.g., overalloc = -0.1 leads to 10% larger arrays),
 *   default: overalloc = 0 (i.e., no extra memory)
 * @param ghost_nzslices fcs_int number of zslices used to create ghost particles
 */
void fcs_gridsort_set_overalloc(fcs_gridsort_t *gs, fcs_float overalloc);

/**
 * @brief set information of particles to sort
 * @param gs fcs_gridsort_t* gridsort object
 * @param nparticles fcs_int local number of particles
 * @param max_nparticles fcs_int max number of particles that can be stored in local particle data arrays
 * @param positions fcs_float* array of local particle positions
 * @param charges fcs_float* array of local particle charges
 */
void fcs_gridsort_set_particles(fcs_gridsort_t *gs, fcs_int nparticles, fcs_int max_nparticles, fcs_float *positions, fcs_float *charges);

/**
 * @brief set max. distance particles are away from the domain of the local process, e.g. because they have moved since the last grid sorting
 * @param gs fcs_gridsort_t* gridsort object
 * @param max_particle_move fcs_float max. distance particles have moved
 */
void fcs_gridsort_set_max_particle_move(fcs_gridsort_t *gs, fcs_float max_particle_move);

/**
 * @brief get information of sorted particles
 * @param gs fcs_gridsort_t* gridsort object
 * @param nparticles fcs_int* pointer to local number of particles (NULL if not required)
 * @param max_nparticles fcs_int* pointer to max number of particles that can be stored in local particle data arrays (NULL if not required)
 * @param positions fcs_float** pointer to local array of particle positions (NULL if not required)
 * @param charges fcs_float** pointer to local array of particle charges (NULL if not required)
 * @param indices fcs_gridsort_index_t** pointer to local array of particle indices (NULL if not required)
 */
void fcs_gridsort_get_sorted_particles(fcs_gridsort_t *gs, fcs_int *nparticles, fcs_int *max_nparticles, fcs_float **positions, fcs_float **charges, fcs_gridsort_index_t **indices);

/**
 * @brief get information of real (non-ghost) particles after sorting and separating (with fcs_gridsort_separate_ghosts)
 * @param gs fcs_gridsort_t* gridsort object
 * @param nparticles fcs_int* pointer to real number of particles (NULL if not required)
 * @param positions fcs_float** pointer to real array of particle positions (NULL if not required)
 * @param charges fcs_float** pointer to real array of particle charges (NULL if not required)
 * @param indices fcs_gridsort_index_t** real pointer to array of prickle indices (NULL if not required)
 */
void fcs_gridsort_get_real_particles(fcs_gridsort_t *gs, fcs_int *nparticles, fcs_float **positions, fcs_float **charges, fcs_gridsort_index_t **indices);

/**
 * @brief get information of ghost particles after sorting and separating (with fcs_gridsort_separate_ghosts)
 * @param gs fcs_gridsort_t* gridsort object
 * @param nparticles fcs_int* pointer to ghost number of particles (NULL if not required)
 * @param positions fcs_float** pointer to ghost array of particle positions (NULL if not required)
 * @param charges fcs_float** pointer to ghost array of particle charges (NULL if not required)
 * @param indices fcs_gridsort_index_t** ghost pointer to array of prickle indices (NULL if not required)
 */
void fcs_gridsort_get_ghost_particles(fcs_gridsort_t *gs, fcs_int *nparticles, fcs_float **positions, fcs_float **charges, fcs_gridsort_index_t **indices);

/**
 * @brief set gridsort cache to store and reuse gridsort data across several gridsort instances
 * @param gs fcs_gridsort_t* gridsort object
 * @param cache fcs_gridsort_cache_t* gridsort cache
 */
void fcs_gridsort_set_cache(fcs_gridsort_t *gs, fcs_gridsort_cache_t *cache);

/**
 * @brief free gridsort cache
 * @param cache fcs_gridsort_cache_t* gridsort cache
 */
void fcs_gridsort_release_cache(fcs_gridsort_cache_t *cache);

/**
 * @brief sort particles to the process grid and create ghost particles
 * @param gs fcs_gridsort_t* gridsort object
 * @param ghost_range fcs_float range around process borders where ghost particles should be created
 * @param comm MPI_Comm Cartesian communicator specifying the process grid
 */
fcs_int fcs_gridsort_sort_forward(fcs_gridsort_t *gs, fcs_float ghost_range, MPI_Comm comm);

/**
 * @brief create ghost particles according to the given ghost range
 * @param gs fcs_gridsort_t* gridsort object
 * @param ghost_range fcs_float range around process borders where ghost particles should be created
 * @param comm MPI_Comm Cartesian communicator specifying the process grid
 */
fcs_int fcs_gridsort_create_ghosts(fcs_gridsort_t *gs, fcs_float ghost_range, MPI_Comm comm);

/**
 * @brief separate real (non-ghost) and ghost particles after sorting
 * @param gs fcs_gridsort_t* gridsort object
 */
void fcs_gridsort_separate_ghosts(fcs_gridsort_t *gs);

/**
 * @brief separate all particles (real and ghosts) according to their z-coordinate into equidistant slices, the slices are ordered according to their z coordinates, the number of local and ghost slices has to be set with fcs_gridsort_set_zslices
 * @param gs fcs_gridsort_t* gridsort object
 * @param zslices_low_ghost_nparticles fcs_int* pointer to an array to return the number of ghost particles in each lower slice (can be NULL, if not NULL, then the array has to be large enough to hold ghost_nzslices entries, see fcs_gridsort_set_zslices)
 * @param zslices_real_nparticles fcs_int* pointer to an array to return the number of real particles in each slice (can be NULL, if not NULL, then the array has to be large enough to hold local_nzslices entries, see fcs_gridsort_set_zslices)
 * @param zslices_high_ghost_nparticles fcs_int* pointer to an array to return the number of ghost particles in each lower slice (can be NULL, if not NULL, then the array has to be large enough to hold ghost_nzslices entries, see fcs_gridsort_set_zslices)
 */
void fcs_gridsort_separate_zslices(fcs_gridsort_t *gs, fcs_int *zslices_low_ghost_nparticles, fcs_int *zslices_real_nparticles, fcs_int *zslices_high_ghost_nparticles);

/**
 * @brief sort particles randomly
 * @param gs fcs_gridsort_t* gridsort object
 * @param comm MPI_Comm communicator
 */
fcs_int fcs_gridsort_sort_random(fcs_gridsort_t *gs, MPI_Comm comm);

/**
 * @brief unfold positions of ghost particles created at periodic boundaries
 * @param nparticles fcs_int local number of particles
 * @param indices fcs_gridsort_index_t* array of local particle indices
 * @param positions fcs_float* array of local particle positions
 * @param box_a fcs_float* 1st base vector of the system box
 * @param box_b fcs_float* 2nd base vector of the system box
 * @param box_c fcs_float* 3rd base vector of the system box
 */
void fcs_gridsort_unfold_periodic_particles(fcs_int nparticles, fcs_gridsort_index_t *indices, fcs_float *positions, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c);

/**
 * @brief set information of results to sort back
 * @param gs fcs_gridsort_t* gridsort object
 * @param max_nparticles fcs_int max number of results that can be stored in local result data arrays
 * @param field fcs_float* array of field values (can be NULL)
 * @param potentials fcs_float* array of potential values (can be NULL)
 */
void fcs_gridsort_set_sorted_results(fcs_gridsort_t *gs, fcs_int max_nresults, fcs_float *field, fcs_float *potentials);

/**
 * @brief set information for storing results of sort back
 * @param gs fcs_gridsort_t* gridsort object
 * @param max_nparticles fcs_int max number of results that can be stored in local result data arrays
 * @param field fcs_float* array for storing field values (can be NULL)
 * @param potentials fcs_float* array for storing potential values (can be NULL)
 */
void fcs_gridsort_set_results(fcs_gridsort_t *gs, fcs_int max_nresults, fcs_float *field, fcs_float *potentials);

/**
 * @brief sort particle information back into their original order
 * @param gs fcs_gridsort_t* gridsort object
 * @param comm MPI_Comm communicator to use
 */
fcs_int fcs_gridsort_sort_backward(fcs_gridsort_t *gs, MPI_Comm comm);

/**
 * @brief prepare resorting of additional particle data, i.e.
 *   (1) enable creation of gridsort_resort object, and
 *   (2) copy sorted potential and field values locally to original potential and field arrays
 * @param gs fcs_gridsort_t* gridsort object
 * @param comm MPI_Comm communicator to use
 */
fcs_int fcs_gridsort_prepare_resort(fcs_gridsort_t *gs, MPI_Comm comm);

/**
 * @brief free arrays allocated by fcs_gridsort_sort_forward (all arrays returned by getters become invalid)
 * @param gs fcs_gridsort_t* gridsort object
 */
void fcs_gridsort_free(fcs_gridsort_t *gs);


#ifdef __cplusplus
}
#endif


#endif /* __GRIDSORT_H__ */
