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


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "common/fcs-common/FCSCommon.h"

#include "sl_forw.h"
#include "sl_back_fp.h"
#include "sl_back_f_.h"
#include "sl_back__p.h"

#include "z_tools.h"
#include "common.h"
#include "gridsort.h"
#include "gridsort_resort.h"


/*#define PRINT_FORWARD_SORTED*/

#define ALLTOALLV_PACKED(_p_, _n_)  ((_p_) >= 1024 && (_n_) <= 200)

/*#define ALLTOALL_SPECIFIC_IN_PLACE*/

#define GRID_DATA_BASE          (0 * 3)
#define GRID_DATA_LOW           (1 * 3)
#define GRID_DATA_HIGH          (2 * 3)
#define GRID_DATA_A             (3 * 3)
#define GRID_DATA_B             (4 * 3)
#define GRID_DATA_C             (5 * 3)
#define GRID_DATA_ZSLICES_LOW   (6 * 3)
#define GRID_DATA_ZSLICES_HIGH  (7 * 3)
#define GRID_DATA_LAST          (8 * 3)


#define BOUNDS_XYZ2COORDS_NAME  bounds_xyz2coords
#include "bounds_xyz2coords.h"

#define BOUNDS_XYZ2COORDS_GHOST
#define BOUNDS_XYZ2COORDS_NAME  bounds_xyz2coords_ghost
#include "bounds_xyz2coords.h"

#define BOUNDS_XYZ2COORDS_PERIODIC
#define BOUNDS_XYZ2COORDS_NAME  bounds_xyz2coords_ghost_periodic
#include "bounds_xyz2coords.h"


/* no ghosts and not periodic */
#undef GRIDSORT_FRONT_TPROC_GHOST
#undef GRIDSORT_FRONT_TPROC_PERIODIC

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_tricl
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_bounds
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#define GRIDSORT_FRONT_TPROC_ZSLICES
#define GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_zslices_zonly
#include "gridsort_front_tproc.h"

/* with ghosts and not periodic */
#define GRIDSORT_FRONT_TPROC_GHOST
#undef GRIDSORT_FRONT_TPROC_PERIODIC

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_tricl
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_bounds
#include "gridsort_front_tproc.h"

/*#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#define GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_zonly
#include "gridsort_front_tproc.h"*/

/* no ghosts and periodic */
#undef GRIDSORT_FRONT_TPROC_GHOST
#define GRIDSORT_FRONT_TPROC_PERIODIC

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_periodic
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_periodic_tricl
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_periodic_bounds
#include "gridsort_front_tproc.h"

/*#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#define GRIDSORT_FRONT_TPROC_ZSLICES
#define GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_periodic_zslices_zonly
#include "gridsort_front_tproc.h"*/

/* with ghosts and periodic */
#define GRIDSORT_FRONT_TPROC_GHOST
#define GRIDSORT_FRONT_TPROC_PERIODIC

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_periodic
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_periodic_tricl
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_periodic_bounds
#include "gridsort_front_tproc.h"

/*#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#define GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_periodic_zslices
#include "gridsort_front_tproc.h"*/


#define VSIZE(_v_)  fcs_sqrt((_v_)[0] * (_v_)[0] + (_v_)[1] * (_v_)[1] + (_v_)[2] * (_v_)[2])

static fcs_float get_ghost_factor(fcs_float *v0, fcs_float *v1, fcs_float *v2, fcs_float ghost_range)
{
  fcs_float f, n[3];


  n[0] = v1[1] * v2[2] - v1[2] * v2[1];
  n[1] = v1[2] * v2[0] - v1[0] * v2[2];
  n[2] = v1[0] * v2[1] - v1[1] * v2[0];

  f = VSIZE(n) * ghost_range / (n[0] * v0[0] + n[1] * v0[1] + n[2] * v0[2]);

  if (f < 0) f *= -1.0;

  return f;
}


static void invert_3x3(fcs_float *v0, fcs_float *v1, fcs_float *v2, fcs_float *iv)
{
  fcs_float det;


  det = v0[0] * v1[1] * v2[2] + v1[0] * v2[1] * v0[2] + v2[0] * v0[1] * v1[2] - v2[0] * v1[1] * v0[2] - v1[0] * v0[1] * v2[2] - v0[0] * v2[1] * v1[2];

  iv[0] = (v1[1] * v2[2] - v2[1] * v1[2]) / det;
  iv[1] = (v2[1] * v0[2] - v0[1] * v2[2]) / det;
  iv[2] = (v0[1] * v1[2] - v1[1] * v0[2]) / det;

  iv[3] = (v2[0] * v1[2] - v1[0] * v2[2]) / det;
  iv[4] = (v0[0] * v2[2] - v2[0] * v0[2]) / det;
  iv[5] = (v1[0] * v0[2] - v0[0] * v1[2]) / det;

  iv[6] = (v1[0] * v2[1] - v2[0] * v1[1]) / det;
  iv[7] = (v2[0] * v0[1] - v0[0] * v2[1]) / det;
  iv[8] = (v0[0] * v1[1] - v1[0] * v0[1]) / det;
}


#ifdef PRINT_FORWARD_SORTED
static void print_particles(fcs_int n, fcs_float *xyz, int size, int rank, MPI_Comm comm)
{
  const int root = 0;
  fcs_int max_n, i, j;
  fcs_float *in_xyz;
  MPI_Status status;
  int in_count;


  MPI_Reduce(&n, &max_n, 1, FCS_MPI_INT, MPI_MAX, root, comm);

  in_xyz = malloc(max_n * 3 * sizeof(fcs_float));

  if (rank == root)
  {
    for (i = 0; i < size; ++i)
    {
      if (i == root) MPI_Sendrecv(xyz, 3 * n, FCS_MPI_FLOAT, root, 0, in_xyz, 3 * max_n, FCS_MPI_FLOAT, root, 0, comm, &status);
      else MPI_Recv(in_xyz, 3 * max_n, FCS_MPI_FLOAT, i, 0, comm, &status);
      MPI_Get_count(&status, FCS_MPI_FLOAT, &in_count);
      
      in_count /= 3;
      
      for (j = 0; j < in_count; ++j) printf("%" FCS_LMOD_INT "d  %" FCS_LMOD_FLOAT "f  %" FCS_LMOD_FLOAT "f  %" FCS_LMOD_FLOAT "f\n", i, in_xyz[j * 3 + 0], in_xyz[j * 3 + 1], in_xyz[j * 3 + 2]);
    }

  } else
  {
    MPI_Send(xyz, 3 * n, FCS_MPI_FLOAT, root, 0, comm);
  }
  
  free(in_xyz);
}
#endif


void fcs_gridsort_create(fcs_gridsort_t *gs)
{
  gs->d.box_base[0] = gs->d.box_base[1] = gs->d.box_base[2] = 0;
  gs->d.box_a[0] = gs->d.box_a[1] = gs->d.box_a[2] = 0;
  gs->d.box_b[0] = gs->d.box_b[1] = gs->d.box_b[2] = 0;
  gs->d.box_c[0] = gs->d.box_c[1] = gs->d.box_c[2] = 0;
  gs->d.periodicity[0] = gs->d.periodicity[1] = gs->d.periodicity[2] = -1;

  gs->d.lower_bounds[0] = gs->d.lower_bounds[1] = gs->d.lower_bounds[2] = -1;
  gs->d.upper_bounds[0] = gs->d.upper_bounds[1] = gs->d.upper_bounds[2] = -1;
  
  gs->d.bounds = NULL;

  gs->local_nzslices = gs->ghost_nzslices = gs->max_ghost_nzslices = 0;

  gs->minalloc = 0;
  gs->overalloc = 0;

  gs->noriginal_particles = gs->max_noriginal_particles = 0;
  gs->original_positions = NULL;
  gs->original_charges = NULL;

  gs->nsorted_particles = gs->max_nsorted_particles = 0;
  gs->sorted_positions = NULL;
  gs->sorted_charges = NULL;
  gs->sorted_indices = NULL;

  gs->nsorted_real_particles = gs->nsorted_ghost_particles = 0;

  gs->nresort_particles = -1;

  gs->max_particle_move = -1;
  gs->nprocs = -1;
  gs->procs = NULL;

  gs->cache = NULL;
}


#ifdef GRIDSORT_PROCLIST
static void release_proclist(fcs_int *nprocs, int **procs);
#endif

void fcs_gridsort_destroy(fcs_gridsort_t *gs)
{
#ifdef GRIDSORT_PROCLIST
  release_proclist(&gs->nprocs, &gs->procs);
#endif
}


void fcs_gridsort_set_system(fcs_gridsort_t *gs, fcs_float *box_base, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_int *periodicity)
{
  fcs_int i;
  
  
  for (i = 0; i < 3; ++i)
  {
    gs->d.box_base[i] = box_base[i];
    gs->d.box_a[i] = box_a[i];
    gs->d.box_b[i] = box_b[i];
    gs->d.box_c[i] = box_c[i];

    if (periodicity) gs->d.periodicity[i] = periodicity[i];
  }
}


void fcs_gridsort_set_bounds(fcs_gridsort_t *gs, fcs_float *lower_bounds, fcs_float *upper_bounds)
{
  fcs_int i;
  
  
  for (i = 0; i < 3; ++i)
  {
    gs->d.lower_bounds[i] = lower_bounds[i];
    gs->d.upper_bounds[i] = upper_bounds[i];
  }
}


void fcs_gridsort_set_zslices(fcs_gridsort_t *gs, fcs_int local_nzslices, fcs_int ghost_nzslices)
{
  gs->local_nzslices = local_nzslices;
  gs->ghost_nzslices = ghost_nzslices;
}


void fcs_gridsort_set_minalloc(fcs_gridsort_t *gs, fcs_int minalloc)
{
  gs->minalloc = minalloc;
}


void fcs_gridsort_set_overalloc(fcs_gridsort_t *gs, fcs_float overalloc)
{
  gs->overalloc = overalloc;
}


void fcs_gridsort_set_particles(fcs_gridsort_t *gs, fcs_int nparticles, fcs_int max_nparticles, fcs_float *positions, fcs_float *charges)
{
  gs->noriginal_particles = nparticles;
  gs->max_noriginal_particles = max_nparticles;
  gs->original_positions = positions;
  gs->original_charges = charges;
}


void fcs_gridsort_set_max_particle_move(fcs_gridsort_t *gs, fcs_float max_particle_move)
{
  gs->max_particle_move = max_particle_move;
}


void fcs_gridsort_get_sorted_particles(fcs_gridsort_t *gs, fcs_int *nparticles, fcs_int *max_nparticles, fcs_float **positions, fcs_float **charges, fcs_gridsort_index_t **indices)
{
  if (nparticles) *nparticles = gs->nsorted_particles;
  if (max_nparticles) *max_nparticles = gs->max_nsorted_particles;
  if (positions) *positions = gs->sorted_positions;
  if (charges) *charges = gs->sorted_charges;
  if (indices) *indices = gs->sorted_indices;
}


void fcs_gridsort_get_real_particles(fcs_gridsort_t *gs, fcs_int *nparticles, fcs_float **positions, fcs_float **charges, fcs_gridsort_index_t **indices)
{
  if (nparticles) *nparticles = gs->nsorted_real_particles;
  if (positions) *positions = gs->sorted_positions;
  if (charges) *charges = gs->sorted_charges;
  if (indices) *indices = gs->sorted_indices;
}


void fcs_gridsort_get_ghost_particles(fcs_gridsort_t *gs, fcs_int *nparticles, fcs_float **positions, fcs_float **charges, fcs_gridsort_index_t **indices)
{
  if (nparticles) *nparticles = gs->nsorted_ghost_particles;
  if (gs->nsorted_ghost_particles > 0)
  {
    if (positions) *positions = gs->sorted_positions + (3 * gs->nsorted_real_particles);
    if (charges) *charges = gs->sorted_charges + gs->nsorted_real_particles;
    if (indices) *indices = gs->sorted_indices + gs->nsorted_real_particles;

  } else
  {
    if (positions) *positions = NULL;
    if (charges) *charges = NULL;
    if (indices) *indices = NULL;
  }
}


void fcs_gridsort_set_cache(fcs_gridsort_t *gs, fcs_gridsort_cache_t *cache)
{
  gs->cache = cache;
}


static void create_cache(fcs_gridsort_cache_t *cache)
{
  *cache = malloc(sizeof(**cache));
}


static void release_bounds(fcs_float **bounds);

void fcs_gridsort_release_cache(fcs_gridsort_cache_t *cache)
{
  if (*cache == FCS_GRIDSORT_CACHE_NULL) return;

  release_bounds(&(*cache)->bounds);

  if (*cache) free(*cache);
  (*cache) = FCS_GRIDSORT_CACHE_NULL;
}


static void update_cache(fcs_gridsort_t *gs)
{
  if (gs->cache == NULL) return;

  if (*gs->cache == FCS_GRIDSORT_CACHE_NULL) create_cache(gs->cache);

  memcpy(*gs->cache, &gs->d, sizeof(gs->d));

  gs->d.bounds = NULL;
}


#define float3_is_equal(_x_, _y_)  (z_fp_is_equal((_x_)[0], (_y_)[0]) && z_fp_is_equal((_x_)[1], (_y_)[1]) && z_fp_is_equal((_x_)[2], (_y_)[2]))
#define int3_is_equal(_x_, _y_)    ((_x_)[0] == (_y_)[0] && (_x_)[1] == (_y_)[1] && (_x_)[2] == (_y_)[2])

static fcs_float *get_cache_bounds(fcs_gridsort_t *gs, int size, int rank, MPI_Comm comm)
{
  fcs_int local_invalid, global_invalid;

  if (gs->cache == NULL || *gs->cache == FCS_GRIDSORT_CACHE_NULL) return NULL;

  local_invalid = (float3_is_equal((*gs->cache)->box_base, gs->d.box_base) &&
                   float3_is_equal((*gs->cache)->box_a, gs->d.box_a) && float3_is_equal((*gs->cache)->box_b, gs->d.box_b) && float3_is_equal((*gs->cache)->box_c, gs->d.box_c) &&
                   int3_is_equal((*gs->cache)->periodicity, gs->d.periodicity) &&
                   float3_is_equal((*gs->cache)->lower_bounds, gs->d.lower_bounds) && float3_is_equal((*gs->cache)->upper_bounds, gs->d.upper_bounds))?0:1;

  MPI_Allreduce(&local_invalid, &global_invalid, 1, FCS_MPI_INT, MPI_SUM, comm);

  if (global_invalid > 0) release_bounds(&(*gs->cache)->bounds);
  
  return (*gs->cache)->bounds;
}

#undef float3_is_equal
#undef int3_is_equal


static fcs_float *setup_bounds(fcs_float *local_bounds, fcs_float *min_bounds, fcs_float *max_bounds, int *cart_dims, int *cart_coords, int size, int rank, MPI_Comm comm)
{
  int root_rank, color, key;
  int root_coords[] = { 0, 0, 0 };

  MPI_Comm line_comm;
  int line_comm_size, line_comm_rank, line_comm_root, iamroot;

  fcs_int dim, i, j;
  fcs_float *bounds, *gather_bounds = NULL;


  bounds = malloc((cart_dims[0] + cart_dims[1] + cart_dims[2] + 3) * sizeof(fcs_float));

  MPI_Cart_rank(comm, root_coords, &root_rank);

  if (cart_coords[1] == root_coords[1] && cart_coords[2] == root_coords[2])
  {
    color = 0;
    key = cart_coords[0];
    dim = 0;

  } else if (cart_coords[2] == root_coords[2] && cart_coords[0] == root_coords[0])
  {
    color = 0;
    key = cart_dims[0] + cart_coords[1];
    dim = 1;

  } else if (cart_coords[0] == root_coords[0] && cart_coords[1] == root_coords[1])
  {
    color = 0;
    key = cart_dims[0] + cart_dims[1] + cart_coords[2];
    dim = 2;

  } else
  {
    color = MPI_UNDEFINED;
    key = 0;
    dim = -1;
  }

  MPI_Comm_split(comm, color, key, &line_comm);

  if (line_comm != MPI_COMM_NULL)
  {
    MPI_Comm_size(line_comm, &line_comm_size);
    MPI_Comm_rank(line_comm, &line_comm_rank);

    DEBUG_CMD(
      printf(DEBUG_PRINT_PREFIX "%d: line comm: %d of %d, dim: %" FCS_LMOD_INT "d, bound: %" FCS_LMOD_FLOAT "f\n", rank, line_comm_rank, line_comm_size, dim, local_bounds[dim]);
    );

    iamroot = (cart_coords[0] == root_coords[0] && cart_coords[1] == root_coords[1] && cart_coords[2] == root_coords[2])?line_comm_rank:0;

    MPI_Allreduce(&iamroot, &line_comm_root, 1, MPI_INT, MPI_SUM, line_comm);

    if (line_comm_rank == line_comm_root) gather_bounds = malloc(line_comm_size * sizeof(fcs_float));

    MPI_Gather(&local_bounds[dim], 1, FCS_MPI_FLOAT, gather_bounds, 1, FCS_MPI_FLOAT, line_comm_root, line_comm);

    if (line_comm_rank == line_comm_root)
    {
      j = 0;
      bounds[j] = min_bounds[0]; ++j;
      bounds[j] = local_bounds[0]; ++j;
      for (i = 0; i < cart_dims[0] - 1; ++i, ++j) bounds[j] = z_max(gather_bounds[j - 1], bounds[j - 1]);

      bounds[j] = min_bounds[1]; ++j;
      bounds[j] = local_bounds[1]; ++j;
      for (i = 0; i < cart_dims[1] - 1; ++i, ++j) bounds[j] = z_max(gather_bounds[j - 3], bounds[j - 1]);

      bounds[j] = min_bounds[2]; ++j;
      bounds[j] = local_bounds[2]; ++j;
      for (i = 0; i < cart_dims[2] - 1; ++i, ++j) bounds[j] = z_max(gather_bounds[j - 5], bounds[j - 1]);
    }

    if (line_comm_rank == line_comm_root) free(gather_bounds);
    
    MPI_Comm_free(&line_comm);
  }

  MPI_Bcast(bounds, cart_dims[0] + cart_dims[1] + cart_dims[2] + 3, FCS_MPI_FLOAT, root_rank, comm);

  DEBUG_CMD(
    if (rank == 0)
    {
      printf(DEBUG_PRINT_PREFIX "cart_dims: %dx%dx%d\n", cart_dims[0], cart_dims[1], cart_dims[2]);
      printf(DEBUG_PRINT_PREFIX " 0-bounds:");
      for (i = 0; i <= cart_dims[0]; ++i) printf("  %" FCS_LMOD_FLOAT "f", bounds[i]);
      printf("\n");
      printf(DEBUG_PRINT_PREFIX " 1-bounds:");
      for (i = 0; i <= cart_dims[1]; ++i) printf("  %" FCS_LMOD_FLOAT "f", bounds[i + cart_dims[0] + 1]);
      printf("\n");
      printf(DEBUG_PRINT_PREFIX " 2-bounds:");
      for (i = 0; i <= cart_dims[2]; ++i) printf("  %" FCS_LMOD_FLOAT "f", bounds[i + cart_dims[0] + 1 + cart_dims[1] + 1]);
      printf("\n");
    }
  );

  return bounds;
}


static void release_bounds(fcs_float **bounds)
{
  if (*bounds) free(*bounds);
  (*bounds) = NULL;
}


#define BOUNDS_XYZ2COORDS(_c_, _xyz_, _base_, _b_, _cd_, _p_, _bs_) Z_MOP( \
  if (periodicity) bounds_xyz2coords_ghost_periodic((_c_), (_xyz_), (_base_), (_b_), (_cd_), (_p_), (_bs_)); \
  else bounds_xyz2coords_ghost((_c_), (_xyz_), (_base_), (_b_), (_cd_)); \
)

#define BOUNDS_X(_b_, _x_)  ((_b_)[(_x_)])
#define BOUNDS_Y(_b_, _x_)  ((_b_)[cart_dims[0] + 1 + (_x_)])
#define BOUNDS_Z(_b_, _x_)  ((_b_)[cart_dims[0] + cart_dims[1] + 2 + (_x_)])


static void setup_max_nparts(fcs_int *max_nparts, fcs_float *ghost_f, fcs_float zslices_ghost_f, int *cart_dims, int *cart_coords, fcs_int *periodicity, fcs_float *bounds, fcs_float *box_size, MPI_Comm comm)
{
  int origin[3], low[3], high[3];
  fcs_float xyz[3], base[3];
  fcs_int local_nparts[3];


  if (bounds)
  {
    xyz[0] = BOUNDS_X(bounds, cart_coords[0]);
    xyz[1] = BOUNDS_Y(bounds, cart_coords[1]);
    xyz[2] = BOUNDS_Z(bounds, cart_coords[2]);

    base[0] = base[1] = base[2] = 0;
    BOUNDS_XYZ2COORDS(origin, xyz, base, bounds, cart_dims, periodicity, box_size);

    base[0] = 2 * ghost_f[0] * box_size[0];
    base[1] = 2 * ghost_f[1] * box_size[1];
    base[2] = 2 * ghost_f[2] * box_size[2];

    BOUNDS_XYZ2COORDS(low, xyz, base, bounds, cart_dims, periodicity, box_size);

    base[0] = -2 * ghost_f[0] * box_size[0];
    base[1] = -2 * ghost_f[1] * box_size[1];
    base[2] = -2 * ghost_f[2] * box_size[2];
    BOUNDS_XYZ2COORDS(high, xyz, base, bounds, cart_dims, periodicity, box_size);

    local_nparts[0] = z_max(origin[0] - low[0], high[0] - origin[0]) + 1;
    local_nparts[1] = z_max(origin[1] - low[1], high[1] - origin[1]) + 1;
    local_nparts[2] = z_max(origin[2] - low[2], high[2] - origin[2]) + 1;

    MPI_Allreduce(local_nparts, max_nparts, 3, FCS_MPI_INT, MPI_MAX, comm);

  } else
  {
    max_nparts[0] = (fcs_int) fcs_ceil(2.0 * ghost_f[0] * cart_dims[0]) + 1;
    max_nparts[1] = (fcs_int) fcs_ceil(2.0 * ghost_f[1] * cart_dims[1]) + 1;
    max_nparts[2] = (fcs_int) fcs_ceil(2.0 * z_max(ghost_f[2], zslices_ghost_f) * cart_dims[2]) + 1;
  }

  max_nparts[3] = max_nparts[0] * max_nparts[1] * max_nparts[2];

/*  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);
  printf("%d: origin: %d,%d,%d - low: %d,%d,%d - high: %d,%d,%d\n", comm_rank, origin[0], origin[1], origin[2], low[0], low[1], low[2], high[0], high[1], high[2]);
  printf("%d: local_nparts: %d,%d,%d - max_nparts: %d,%d,%d\n", comm_rank, (int) local_nparts[0], (int) local_nparts[1], (int) local_nparts[2], (int) max_nparts[0], (int) max_nparts[1], (int) max_nparts[2]);*/
}


#if defined(GRIDSORT_FRONT_TPROC_RANK_CACHE) || defined(GRIDSORT_PROCLIST)
static void low_high_coords(int *low_coords, int *high_coords, fcs_float *ghost_f, fcs_float zslices_ghost_f, fcs_float *move_f, int *cart_dims, int *cart_coords, fcs_int *periodicity, fcs_float *bounds, fcs_float *box_size)
{
  fcs_float low_xyz[3], high_xyz[3], base[3];
  int ghost_move_coords[3];


  if (bounds)
  {
    if (move_f)
    {
      low_xyz[0] = BOUNDS_X(bounds, cart_coords[0]);
      low_xyz[1] = BOUNDS_Y(bounds, cart_coords[1]);
      low_xyz[2] = BOUNDS_Z(bounds, cart_coords[2]);

      high_xyz[0] = BOUNDS_X(bounds, cart_coords[0] + 1);
      high_xyz[1] = BOUNDS_Y(bounds, cart_coords[1] + 1);
      high_xyz[2] = BOUNDS_Z(bounds, cart_coords[2] + 1);

      base[0] = (      ghost_f[0]                   + move_f[0]) * box_size[0];
      base[1] = (      ghost_f[1]                   + move_f[1]) * box_size[1];
      base[2] = (z_max(ghost_f[2], zslices_ghost_f) + move_f[2]) * box_size[2];

    } else
    {
      low_xyz[0] = BOUNDS_X(bounds, 0);
      low_xyz[1] = BOUNDS_Y(bounds, 0);
      low_xyz[2] = BOUNDS_Z(bounds, 0);

      high_xyz[0] = BOUNDS_X(bounds, cart_dims[0]);
      high_xyz[1] = BOUNDS_Y(bounds, cart_dims[1]);
      high_xyz[2] = BOUNDS_Z(bounds, cart_dims[2]);

      base[0] =       ghost_f[0]                   * box_size[0];
      base[1] =       ghost_f[1]                   * box_size[1];
      base[2] = z_max(ghost_f[2], zslices_ghost_f) * box_size[2];
    }

    BOUNDS_XYZ2COORDS(low_coords, low_xyz, base, bounds, cart_dims, periodicity, box_size);

    base[0] *= -1;
    base[1] *= -1;
    base[2] *= -1;

    BOUNDS_XYZ2COORDS(high_coords, high_xyz, base, bounds, cart_dims, periodicity, box_size);

  } else
  {
    if (move_f)
    {
      ghost_move_coords[0] = (int) fcs_ceil((ghost_f[0] + move_f[0]) * cart_dims[0]);
      ghost_move_coords[1] = (int) fcs_ceil((ghost_f[1] + move_f[1]) * cart_dims[1]);
      ghost_move_coords[2] = (int) fcs_ceil((z_max(ghost_f[2], zslices_ghost_f) + move_f[2]) * cart_dims[2]);

      low_coords[0] = cart_coords[0] - ghost_move_coords[0];
      low_coords[1] = cart_coords[1] - ghost_move_coords[1];
      low_coords[2] = cart_coords[2] - ghost_move_coords[2];

      high_coords[0] = cart_coords[0] + ghost_move_coords[0];
      high_coords[1] = cart_coords[1] + ghost_move_coords[1];
      high_coords[2] = cart_coords[2] + ghost_move_coords[2];

    } else
    {
      ghost_move_coords[0] = (int) fcs_ceil(ghost_f[0] * cart_dims[0]);
      ghost_move_coords[1] = (int) fcs_ceil(ghost_f[1] * cart_dims[1]);
      ghost_move_coords[2] = (int) fcs_ceil(z_max(ghost_f[2], zslices_ghost_f) * cart_dims[2]);

      low_coords[0] = 0 - ghost_move_coords[0];
      low_coords[1] = 0 - ghost_move_coords[1];
      low_coords[2] = 0 - ghost_move_coords[2];

      high_coords[0] = cart_dims[0] - 1 + ghost_move_coords[0];
      high_coords[1] = cart_dims[1] - 1 + ghost_move_coords[1];
      high_coords[2] = cart_dims[2] - 1 + ghost_move_coords[2];
    }
  }

/*  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  printf("%d: coords: %d,%d,%d - %d,%d,%d\n", comm_rank, low_coords[0], low_coords[1], low_coords[2], high_coords[0], high_coords[1], high_coords[2]);*/
}
#endif


#ifdef GRIDSORT_FRONT_TPROC_RANK_CACHE
static int *setup_rank_cache(fcs_int *rank_cache_offset, fcs_int *rank_cache_sizes, fcs_int *rank_cache_size, fcs_float *ghost_f, fcs_float zslices_ghost_f, fcs_float *move_f, int *cart_dims, int *cart_coords, fcs_int *periodicity, fcs_float *bounds, fcs_float *box_size)
{
  fcs_int i;
  int *rank_cache;

  int low_coords[3], high_coords[3];


  low_high_coords(low_coords, high_coords, ghost_f, zslices_ghost_f, move_f, cart_dims, cart_coords, periodicity, bounds, box_size);

  if (!periodicity[0]) { low_coords[0] = z_max(0, low_coords[0]); high_coords[0] = z_min(high_coords[0], cart_dims[0] - 1); }
  if (!periodicity[1]) { low_coords[1] = z_max(0, low_coords[1]); high_coords[1] = z_min(high_coords[1], cart_dims[1] - 1); }
  if (!periodicity[2]) { low_coords[2] = z_max(0, low_coords[2]); high_coords[2] = z_min(high_coords[2], cart_dims[2] - 1); }

  rank_cache_sizes[0] = high_coords[0] - low_coords[0] + 1;
  rank_cache_sizes[1] = high_coords[1] - low_coords[1] + 1;
  rank_cache_sizes[2] = high_coords[2] - low_coords[2] + 1;

  *rank_cache_offset = (0 - low_coords[0]) + rank_cache_sizes[0] * ((0 - low_coords[1]) + rank_cache_sizes[1] * (0 - low_coords[2]));

  *rank_cache_size = rank_cache_sizes[0] * rank_cache_sizes[1] * rank_cache_sizes[2];

  rank_cache = malloc(*rank_cache_size * sizeof(int));

  for (i = 0; i < *rank_cache_size; ++i) rank_cache[i] = MPI_UNDEFINED;

  return rank_cache;
}
#endif


#ifdef GRIDSORT_PROCLIST
static int *setup_proclist(fcs_int *nprocs, fcs_int max_nprocs, fcs_float *ghost_f, fcs_float zslices_ghost_f, fcs_float *move_f, int *cart_dims, int *cart_coords, fcs_int *periodicity, fcs_float *bounds, fcs_float *box_size, MPI_Comm comm)
{
  int *procs, low_coords[3], high_coords[3], offsets[3], nums[3], x[3], coords[3];
  fcs_int local_skip, global_skip;


  low_high_coords(low_coords, high_coords, ghost_f, zslices_ghost_f, move_f, cart_dims, cart_coords, periodicity, bounds, box_size);

#define MAP_COORD(_c_, _d_)  ((((_c_) % (_d_)) + (_d_)) % (_d_))

  if (!periodicity[0]) { offsets[0] = z_max(0, low_coords[0]); nums[0] = z_min(high_coords[0], cart_dims[0] - 1) - offsets[0] + 1; }
  else { offsets[0] = MAP_COORD(low_coords[0],  cart_dims[0]); nums[0] = z_min(high_coords[0] - low_coords[0] + 1, cart_dims[0]); }

  if (!periodicity[1]) { offsets[1] = z_max(0, low_coords[1]); nums[1] = z_min(high_coords[1], cart_dims[1] - 1) - offsets[1] + 1; }
  else { offsets[1] = MAP_COORD(low_coords[1],  cart_dims[1]); nums[1] = z_min(high_coords[1] - low_coords[1] + 1, cart_dims[1]); }

  if (!periodicity[2]) { offsets[2] = z_max(0, low_coords[2]); nums[2] = z_min(high_coords[2], cart_dims[2] - 1) - offsets[2] + 1; }
  else { offsets[2] = MAP_COORD(low_coords[2],  cart_dims[2]); nums[2] = z_min(high_coords[2] - low_coords[2] + 1, cart_dims[2]); }

#undef MAP_COORD

/*  printf("%d,%d,%d - %d,%d,%d -> %d,%d,%d - %d,%d,%d\n",
    low_coords[0], low_coords[1], low_coords[2], high_coords[0], high_coords[1], high_coords[2],
    offsets[0], offsets[1], offsets[2], nums[0], nums[1], nums[2]
  );*/

  local_skip = (nums[0] * nums[1] * nums[2] > max_nprocs)?1:0;
  MPI_Allreduce(&local_skip, &global_skip, 1, FCS_MPI_INT, MPI_SUM, comm);
  if (global_skip > 0)
  {
    *nprocs = nums[0] * nums[1] * nums[2];
    return NULL;
  }

  procs = malloc(nums[0] * nums[1] * nums[2] * sizeof(int));

  *nprocs = 0;

  for (x[0] = 0; x[0] < nums[0]; ++x[0])
  for (x[1] = 0; x[1] < nums[1]; ++x[1])
  for (x[2] = 0; x[2] < nums[2]; ++x[2])
  {
    coords[0] = offsets[0] + x[0];
    coords[1] = offsets[1] + x[1];
    coords[2] = offsets[2] + x[2];

    MPI_Cart_rank(comm, coords, &procs[*nprocs]);
/*    printf("%d\n", procs[*nprocs]);*/
    ++(*nprocs);
  }

  return procs;
}


static void release_proclist(fcs_int *nprocs, int **procs)
{
  *nprocs = -1;
  if (*procs) free(*procs);
  *procs = NULL;
}


#ifdef GRIDSORT_PROCLIST_VERIFY
static void verify_proclist(fcs_int nprocs, int *procs, int size, int rank, MPI_Comm comm)
{
  fcs_int i, nrecvs, failed, global_failed;
  int *sendtags, *recvtags;
  
  
  sendtags = malloc(size * sizeof(int));
  recvtags = malloc(size * sizeof(int));

  for (i = 0; i < size; ++i) sendtags[i] = 0;
  for (i = 0; i < nprocs; ++i) sendtags[procs[i]] = 1;
  
  MPI_Alltoall(sendtags, 1, MPI_INT, recvtags, 1, MPI_INT, comm);
  
  nrecvs = 0;
  for (i = 0; i < size; ++i) nrecvs += recvtags[i];

  failed = (nrecvs != nprocs);
  for (i = 0; i < nprocs; ++i) failed |= (recvtags[procs[i]] == 0);

  DEBUG_CMD(
    printf(DEBUG_PRINT_PREFIX "%d: verify proclist: %s (%d vs. %d)\n", rank, (failed)?"FAILED":"OK", (int) nrecvs, (int) nprocs);
  );

  free(sendtags);
  free(recvtags);

  failed = (failed)?1:0;
  
  MPI_Allreduce(&failed, &global_failed, 1, FCS_MPI_INT, MPI_SUM, comm);

  if (global_failed)
  {
    printf("%d: nprocs: %d\n", rank, (int) nprocs);
    for (i = 0; i < nprocs; ++i)
      printf("%d: %d: %d\n", rank, (int) i, procs[i]);

    printf("%d: sendtags / recvtags\n", rank);
    for (i = 0; i < size; ++i)
      printf("%d: %d: %d / %d\n", rank, (int) i, sendtags[i], recvtags[i]);

    MPI_Abort(MPI_COMM_WORLD, 28);
  }
}
#endif

#endif /* GRIDSORT_PROCLIST */


fcs_int fcs_gridsort_sort_forward(fcs_gridsort_t *gs, fcs_float ghost_range, MPI_Comm comm)
{
  int comm_size, comm_rank;

  fcs_gridsort_index_t *original_indices, index_rank;

  fcs_int i;

  fcs_forw_elements_t sin0, sout0;

  fcs_int periodicity[3];
  int cart_dims[3], cart_periods[3], cart_coords[3], topo_status;

  fcs_float iv[9];
  fcs_float grid_data[GRID_DATA_LAST];

  fcs_int max_nparts[4];

#ifdef GRIDSORT_FRONT_TPROC_RANK_CACHE
  int *rank_cache;
  fcs_int rank_cache_offset, rank_cache_sizes[3], rank_cache_size;
#endif

  void *grid_tproc_data[] = { &comm, grid_data, cart_dims, cart_coords, periodicity, NULL, NULL, &max_nparts[3], NULL, NULL };

  fcs_forw_tproc_t tproc;
  fcs_forw_tproc_f *tproc_func = NULL;
  fcs_forw_tprocs_mod_f *tprocs_mod_func = NULL;

  fcs_int with_ghost, with_periodic, with_triclinic, with_bounds, with_zslices;
  fcs_float ghost_f[3], zslices_ghost_range, zslices_ghost_f;

#if defined(GRIDSORT_FRONT_TPROC_RANK_CACHE) || defined(GRIDSORT_PROCLIST)
  fcs_float max_particle_move, move_f[3];
#endif

  fcs_float *min_bounds, max_bounds[3], box_size[3];

#ifdef ALLTOALLV_PACKED
  fcs_int local_packed, global_packed, original_packed;
#endif

  fcs_forw_slint_t old_minalloc;
  double old_overalloc;
  
#ifdef DO_TIMING
  double t[2] = { 0, 0 };
#endif


  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  gs->nsorted_particles = gs->max_nsorted_particles = 0;
  gs->sorted_positions = NULL;
  gs->sorted_charges = NULL;
  gs->sorted_indices = NULL;

  MPI_Topo_test(comm, &topo_status);

  if (topo_status != MPI_CART)
  {
    fprintf(stderr, "ERROR: no Cartesian communicator available for gridsort!\n");
    return -1;
  }

  MPI_Cart_get(comm, 3, cart_dims, cart_periods, cart_coords);

  periodicity[0] = (gs->d.periodicity[0])?cart_periods[0]:0;
  periodicity[1] = (gs->d.periodicity[1])?cart_periods[1]:0;
  periodicity[2] = (gs->d.periodicity[2])?cart_periods[2]:0;

  with_ghost = (ghost_range > 0);
  with_periodic = (periodicity[0] || periodicity[1] || periodicity[2]);
  with_triclinic = z_is_triclinic(gs->d.box_a, gs->d.box_b, gs->d.box_c);
  with_bounds = (gs->d.lower_bounds[0] >= 0 && gs->d.lower_bounds[1] >= 0 && gs->d.lower_bounds[2] >= 0 && gs->d.upper_bounds[0] >= 0 && gs->d.upper_bounds[1] >= 0 && gs->d.upper_bounds[2] >= 0);
  with_zslices = (gs->local_nzslices > 0);

  INFO_CMD(
    if (comm_rank == 0)
    {
      printf(INFO_PRINT_PREFIX "gridsort forward settings:\n");
      printf(INFO_PRINT_PREFIX " ghost: %s\n", with_ghost?"yes":"no");
      printf(INFO_PRINT_PREFIX " periodic: %s\n", with_periodic?"yes":"no");
      printf(INFO_PRINT_PREFIX " triclinic: %s\n", with_triclinic?"yes":"no");
      printf(INFO_PRINT_PREFIX " bounds: %s\n", with_bounds?"yes":"no");
      printf(INFO_PRINT_PREFIX " zslices: %s\n", with_zslices?"yes":"no");
      printf(INFO_PRINT_PREFIX " cartesian grid:\n");
      printf(INFO_PRINT_PREFIX "  dims: %dx%dx%d\n", cart_dims[0], cart_dims[1], cart_dims[2]);
      printf(INFO_PRINT_PREFIX "  periods: %dx%dx%d\n", cart_periods[0], cart_periods[1], cart_periods[2]);
    }
  );

  if (with_bounds && with_zslices)
  {
    fprintf(stderr, "ERROR: either lower and upper bounds OR zslices can be used for gridsort\n");
    return -1;
  }

  if (with_bounds && with_triclinic)
  {
    fprintf(stderr, "ERROR: box shape ([%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f] x [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f]"
      " x [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f]) not supported when using lower and upper bounds for gridsort\n",
      gs->d.box_a[0], gs->d.box_a[1], gs->d.box_a[2], gs->d.box_b[0], gs->d.box_b[1], gs->d.box_b[2], gs->d.box_c[0], gs->d.box_c[1], gs->d.box_c[2]);
    return -1;
  }

  if (with_zslices)
  {
    if (with_triclinic)
    {
      fprintf(stderr, "ERROR: box shape ([%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f] x [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f]"
        " x [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f]) not supported when using zslices for gridsort\n",
        gs->d.box_a[0], gs->d.box_a[1], gs->d.box_a[2], gs->d.box_b[0], gs->d.box_b[1], gs->d.box_b[2], gs->d.box_c[0], gs->d.box_c[1], gs->d.box_c[2]);
      return -1;
    }
    
    if (with_ghost)
    {
      fprintf(stderr, "WARNING: non-zero ghost_range (%" FCS_LMOD_FLOAT "f) is not supported when using zslices for gridsort\n", ghost_range);
      ghost_range = 0.0;
      with_ghost = 0;
    }
  }

  ghost_f[0] = get_ghost_factor(gs->d.box_a, gs->d.box_b, gs->d.box_c, ghost_range);
  ghost_f[1] = get_ghost_factor(gs->d.box_b, gs->d.box_c, gs->d.box_a, ghost_range);
  ghost_f[2] = get_ghost_factor(gs->d.box_c, gs->d.box_a, gs->d.box_b, ghost_range);

  INFO_CMD(
    if (comm_rank == 0)
      printf(INFO_PRINT_PREFIX "ghost_range: %" FCS_LMOD_FLOAT "f, ghost_f: %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f\n",
        ghost_range, ghost_f[0], ghost_f[1], ghost_f[2]);
  );

#if defined(GRIDSORT_FRONT_TPROC_RANK_CACHE) || defined(GRIDSORT_PROCLIST)
  MPI_Allreduce(&gs->max_particle_move, &max_particle_move, 1, FCS_MPI_FLOAT, MPI_MAX, comm);

  INFO_CMD(
    if (comm_rank == 0)
      printf(INFO_PRINT_PREFIX "max_particle_move: %" FCS_LMOD_FLOAT "f\n", max_particle_move);
  );

  move_f[0] = get_ghost_factor(gs->d.box_a, gs->d.box_b, gs->d.box_c, z_max(0, max_particle_move));
  move_f[1] = get_ghost_factor(gs->d.box_b, gs->d.box_c, gs->d.box_a, z_max(0, max_particle_move));
  move_f[2] = get_ghost_factor(gs->d.box_c, gs->d.box_a, gs->d.box_b, z_max(0, max_particle_move));

  INFO_CMD(
    if (comm_rank == 0)
      printf(INFO_PRINT_PREFIX "max_particle_move: %" FCS_LMOD_FLOAT "f, move_f: %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f\n",
        max_particle_move, move_f[0], move_f[1], move_f[2]);
  );
#endif

  if (with_zslices)
  {
    zslices_ghost_range = gs->ghost_nzslices * gs->d.box_c[2] / (cart_dims[2] * gs->local_nzslices);

    zslices_ghost_f = get_ghost_factor(gs->d.box_c, gs->d.box_a, gs->d.box_b, zslices_ghost_range);

    gs->max_ghost_nzslices = (fcs_int) fcs_ceil(z_max(ghost_range, zslices_ghost_range) / (gs->d.box_c[2] / (cart_dims[2] * gs->local_nzslices)));

    INFO_CMD(
      if (comm_rank == 0)
        printf(INFO_PRINT_PREFIX "zslices_ghost_range: %" FCS_LMOD_FLOAT "f, zslices_ghost_f: %" FCS_LMOD_FLOAT "f, max_ghost_nzslices: %" FCS_LMOD_INT "d\n",
          zslices_ghost_range, zslices_ghost_f, gs->max_ghost_nzslices);
    );

  } else zslices_ghost_f = 0;

  grid_data[GRID_DATA_BASE + 0] = gs->d.box_base[0];
  grid_data[GRID_DATA_BASE + 1] = gs->d.box_base[1];
  grid_data[GRID_DATA_BASE + 2] = gs->d.box_base[2];
  grid_data[GRID_DATA_LOW + 0] = gs->d.box_base[0] + gs->d.box_a[0] * ghost_f[0] + gs->d.box_b[0] * ghost_f[1] + gs->d.box_c[0] * ghost_f[2];
  grid_data[GRID_DATA_LOW + 1] = gs->d.box_base[1] + gs->d.box_a[1] * ghost_f[0] + gs->d.box_b[1] * ghost_f[1] + gs->d.box_c[1] * ghost_f[2];
  grid_data[GRID_DATA_LOW + 2] = gs->d.box_base[2] + gs->d.box_a[2] * ghost_f[0] + gs->d.box_b[2] * ghost_f[1] + gs->d.box_c[2] * ghost_f[2];
  grid_data[GRID_DATA_HIGH + 0] = gs->d.box_base[0] - gs->d.box_a[0] * ghost_f[0] - gs->d.box_b[0] * ghost_f[1] - gs->d.box_c[0] * ghost_f[2];
  grid_data[GRID_DATA_HIGH + 1] = gs->d.box_base[1] - gs->d.box_a[1] * ghost_f[0] - gs->d.box_b[1] * ghost_f[1] - gs->d.box_c[1] * ghost_f[2];
  grid_data[GRID_DATA_HIGH + 2] = gs->d.box_base[2] - gs->d.box_a[2] * ghost_f[0] - gs->d.box_b[2] * ghost_f[1] - gs->d.box_c[2] * ghost_f[2];

  if (with_zslices)
  {
    grid_data[GRID_DATA_ZSLICES_LOW + 0] = 0;
    grid_data[GRID_DATA_ZSLICES_LOW + 1] = 0;
    grid_data[GRID_DATA_ZSLICES_LOW + 2] = gs->d.box_base[2] + gs->d.box_c[2] * zslices_ghost_f;
    grid_data[GRID_DATA_ZSLICES_HIGH + 0] = 0;
    grid_data[GRID_DATA_ZSLICES_HIGH + 1] = 0;
    grid_data[GRID_DATA_ZSLICES_HIGH + 2] = gs->d.box_base[2] - gs->d.box_c[2] * zslices_ghost_f;
  }

  if (with_bounds)
  {
    min_bounds = gs->d.box_base;
    max_bounds[0] = gs->d.box_base[0] + gs->d.box_a[0];
    max_bounds[1] = gs->d.box_base[1] + gs->d.box_b[1];
    max_bounds[2] = gs->d.box_base[2] + gs->d.box_c[2];

    box_size[0] = gs->d.box_a[0];
    box_size[1] = gs->d.box_b[1];
    box_size[2] = gs->d.box_c[2];

    DEBUG_CMD(
      printf(DEBUG_PRINT_PREFIX "%d: my bounds: %" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f - %" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f\n",
        comm_rank, gs->d.lower_bounds[0], gs->d.lower_bounds[1], gs->d.lower_bounds[2], gs->d.upper_bounds[0], gs->d.upper_bounds[1], gs->d.upper_bounds[2]);
    );

    gs->d.bounds = get_cache_bounds(gs, comm_size, comm_rank, comm);

    DEBUG_CMD(
      printf(DEBUG_PRINT_PREFIX "%d: cached bounds: %p\n", comm_rank, gs->d.bounds);
    );

    if (gs->d.bounds == NULL) gs->d.bounds = setup_bounds(gs->d.upper_bounds, min_bounds, max_bounds, cart_dims, cart_coords, comm_size, comm_rank, comm);

    for (i = 0; i <= cart_dims[0]; ++i) gs->d.bounds[i] -= gs->d.box_base[0];
    for (i = 0; i <= cart_dims[1]; ++i) gs->d.bounds[i + cart_dims[0] + 1] -= gs->d.box_base[1];
    for (i = 0; i <= cart_dims[2]; ++i) gs->d.bounds[i + cart_dims[0] + cart_dims[1] + 2] -= gs->d.box_base[2];

    grid_tproc_data[5] = gs->d.bounds;
    grid_tproc_data[6] = box_size;

    gs->sub_box_base[0] = gs->d.box_base[0] + BOUNDS_X(gs->d.bounds, cart_coords[0]);
    gs->sub_box_base[1] = gs->d.box_base[1] + BOUNDS_Y(gs->d.bounds, cart_coords[1]);
    gs->sub_box_base[2] = gs->d.box_base[2] + BOUNDS_Z(gs->d.bounds, cart_coords[2]);

    gs->sub_box_a[0] = gs->d.box_base[0] + BOUNDS_X(gs->d.bounds, cart_coords[0] + 1) - BOUNDS_X(gs->d.bounds, cart_coords[0]);
    gs->sub_box_a[1] = 0;
    gs->sub_box_a[2] = 0;

    gs->sub_box_b[0] = 0;
    gs->sub_box_b[1] = gs->d.box_base[1] + BOUNDS_Y(gs->d.bounds, cart_coords[1] + 1) - BOUNDS_Y(gs->d.bounds, cart_coords[1]);
    gs->sub_box_b[2] = 0;

    gs->sub_box_c[0] = 0;
    gs->sub_box_c[1] = 0;
    gs->sub_box_c[2] = gs->d.box_base[2] + BOUNDS_Z(gs->d.bounds, cart_coords[2] + 1) - BOUNDS_Z(gs->d.bounds, cart_coords[2]);

  } else
  {
    invert_3x3(gs->d.box_a, gs->d.box_b, gs->d.box_c, iv);
  
    grid_data[GRID_DATA_A + 0] = cart_dims[0] * iv[0];
    grid_data[GRID_DATA_A + 1] = cart_dims[0] * iv[3];
    grid_data[GRID_DATA_A + 2] = cart_dims[0] * iv[6];
    grid_data[GRID_DATA_B + 0] = cart_dims[1] * iv[1];
    grid_data[GRID_DATA_B + 1] = cart_dims[1] * iv[4];
    grid_data[GRID_DATA_B + 2] = cart_dims[1] * iv[7];
    grid_data[GRID_DATA_C + 0] = cart_dims[2] * iv[2];
    grid_data[GRID_DATA_C + 1] = cart_dims[2] * iv[5];
    grid_data[GRID_DATA_C + 2] = cart_dims[2] * iv[8];

    DEBUG_CMD(
      if (comm_rank == 0)
      {
        printf(DEBUG_PRINT_PREFIX "A: %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f\n", grid_data[GRID_DATA_A + 0], grid_data[GRID_DATA_A + 1], grid_data[GRID_DATA_A + 2]);
        printf(DEBUG_PRINT_PREFIX "B: %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f\n", grid_data[GRID_DATA_B + 0], grid_data[GRID_DATA_B + 1], grid_data[GRID_DATA_B + 2]);
        printf(DEBUG_PRINT_PREFIX "C: %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f\n", grid_data[GRID_DATA_C + 0], grid_data[GRID_DATA_C + 1], grid_data[GRID_DATA_C + 2]);
      }
    );

    gs->sub_box_a[0] = (gs->d.box_a[0] / cart_dims[0]);
    gs->sub_box_a[1] = (gs->d.box_a[1] / cart_dims[0]);
    gs->sub_box_a[2] = (gs->d.box_a[2] / cart_dims[0]);

    gs->sub_box_b[0] = (gs->d.box_b[0] / cart_dims[1]);
    gs->sub_box_b[1] = (gs->d.box_b[1] / cart_dims[1]);
    gs->sub_box_b[2] = (gs->d.box_b[2] / cart_dims[1]);

    gs->sub_box_c[0] = (gs->d.box_c[0] / cart_dims[2]);
    gs->sub_box_c[1] = (gs->d.box_c[1] / cart_dims[2]);
    gs->sub_box_c[2] = (gs->d.box_c[2] / cart_dims[2]);

    gs->sub_box_base[0] = gs->d.box_base[0] + (cart_coords[0] * gs->sub_box_a[0]) + (cart_coords[1] * gs->sub_box_b[0]) + (cart_coords[2] * gs->sub_box_c[0]);
    gs->sub_box_base[1] = gs->d.box_base[1] + (cart_coords[0] * gs->sub_box_a[1]) + (cart_coords[1] * gs->sub_box_b[1]) + (cart_coords[2] * gs->sub_box_c[1]);
    gs->sub_box_base[2] = gs->d.box_base[2] + (cart_coords[0] * gs->sub_box_a[2]) + (cart_coords[1] * gs->sub_box_b[2]) + (cart_coords[2] * gs->sub_box_c[2]);
  }

  index_rank = GRIDSORT_INDEX_VAL_PROC(comm_rank);

  fcs_forw_SL_DEFCON(mpi.rank) = comm_rank;

  fcs_forw_mpi_datatypes_init();

  original_indices = malloc(gs->noriginal_particles * sizeof(fcs_gridsort_index_t));

  for (i = 0; i < gs->noriginal_particles; ++i) original_indices[i] = index_rank + i;

  fcs_forw_elem_set_size(&sin0, gs->noriginal_particles);
  fcs_forw_elem_set_max_size(&sin0, gs->noriginal_particles);
  fcs_forw_elem_set_keys(&sin0, original_indices);
  fcs_forw_elem_set_data(&sin0, gs->original_positions, gs->original_charges);

  fcs_forw_elem_set_size(&sout0, 0);
  fcs_forw_elem_set_max_size(&sout0, 0);
  fcs_forw_elem_set_keys(&sout0, NULL);
  fcs_forw_elem_set_data(&sout0, NULL, NULL);

  setup_max_nparts(max_nparts, ghost_f, zslices_ghost_f, cart_dims, cart_coords, periodicity, gs->d.bounds, box_size, comm);

  INFO_CMD(
    if (comm_rank == 0)
      printf(INFO_PRINT_PREFIX "max_nparts = %" FCS_LMOD_INT "d * %" FCS_LMOD_INT "d * %" FCS_LMOD_INT "d = %" FCS_LMOD_INT "d\n", max_nparts[0], max_nparts[1], max_nparts[2], max_nparts[3]);
  );

#ifdef GRIDSORT_FRONT_TPROC_RANK_CACHE
  rank_cache = setup_rank_cache(&rank_cache_offset, rank_cache_sizes, &rank_cache_size, ghost_f, zslices_ghost_f, (gs->max_particle_move >= 0)?move_f:NULL, cart_dims, cart_coords, periodicity, gs->d.bounds, box_size);

  grid_tproc_data[8] = rank_cache + rank_cache_offset;
  grid_tproc_data[9] = rank_cache_sizes;
#endif

  if (with_bounds)
  {
    if (with_ghost)
    {
      if (with_periodic) tprocs_mod_func = gridsort_front_tproc_ghost_periodic_bounds;
      else tprocs_mod_func = gridsort_front_tproc_ghost_bounds;

    } else
    {
      if (with_periodic) tproc_func = gridsort_front_tproc_periodic_bounds;
      else tproc_func = gridsort_front_tproc_bounds;
    }

  } else if (with_zslices)
  {
/* ghost_range AND zslices are currently not supported */
/*    if (with_ghost)
    {
      if (ghost_range <= zslices_ghost_range)
      {
        if (with_periodic) tproc_func = gridsort_front_tproc_ghost_periodic_zslices;
        else tproc_func = gridsort_front_tproc_zslices_zonly;

      } else
      {
        if (with_periodic) tproc_func = gridsort_front_tproc_ghost_periodic;
        else tproc_func = gridsort_front_tproc_ghost_zonly;
      }
    } else
*/
    {
/* zslices (without ghost_range) do not require any periodicty */
/*      if (with_periodic) tproc_func = gridsort_front_tproc_periodic_zslices_zonly;
      else*/ tprocs_mod_func = gridsort_front_tproc_zslices_zonly;
    }
  } else
  {
    if (with_triclinic)
    {
      if (with_ghost)
      {
        if (with_periodic) tprocs_mod_func = gridsort_front_tproc_ghost_periodic_tricl;
        else tprocs_mod_func = gridsort_front_tproc_ghost_tricl;

      } else
      {
        if (with_periodic) tproc_func = gridsort_front_tproc_periodic_tricl;
        else tproc_func = gridsort_front_tproc_tricl;
      }

    } else
    {
      if (with_ghost)
      {
        if (with_periodic) tprocs_mod_func = gridsort_front_tproc_ghost_periodic;
        else tprocs_mod_func = gridsort_front_tproc_ghost;

      } else
      {
        if (with_periodic) tproc_func = gridsort_front_tproc_periodic;
        else tproc_func = gridsort_front_tproc;
      }
    }
  }

  if (tproc_func) fcs_forw_tproc_create_tproc(&tproc, tproc_func, fcs_forw_TPROC_RESET_NULL, fcs_forw_TPROC_EXDEF_NULL);
  else if (tprocs_mod_func) fcs_forw_tproc_create_tprocs_mod(&tproc, tprocs_mod_func, fcs_forw_TPROC_RESET_NULL, fcs_forw_TPROC_EXDEF_NULL);
  else
  {
    fprintf(stderr, "ERROR: no tproc function selected in gridsort, something went seriously wrong!\n");
    return -1;
  }

#ifdef GRIDSORT_PROCLIST
  release_proclist(&gs->nprocs, &gs->procs);

  if (max_particle_move >= 0)
  {
    gs->procs = setup_proclist(&gs->nprocs, comm_size / 2, ghost_f, zslices_ghost_f, move_f, cart_dims, cart_coords, periodicity, gs->d.bounds, box_size, comm);

#ifdef GRIDSORT_PROCLIST_VERIFY
    if (gs->procs) verify_proclist(gs->nprocs, gs->procs, comm_size, comm_rank, comm);
#endif

    DEBUG_CMD(
      printf(DEBUG_PRINT_PREFIX "%d: nprocs: %" fcs_forw_slint_fmt "\n", comm_rank, gs->nprocs);
    );

    INFO_CMD(
    if (comm_rank == 0)
      printf(INFO_PRINT_PREFIX "%d: proclist: nprocs: %" FCS_LMOD_INT "d, procs: %p\n", comm_rank, gs->nprocs, gs->procs);
    );

# ifdef GRIDSORT_FRONT_PROCLIST
    if (gs->procs) fcs_forw_tproc_set_proclists(&tproc, gs->nprocs, gs->procs, gs->nprocs, gs->procs, comm_size, comm_rank, comm);
# endif
  }
#endif

  DEBUG_CMD(
    printf(DEBUG_PRINT_PREFIX "%d: keys in: %" fcs_forw_slint_fmt "\n", comm_rank, sin0.size);
  );

#if 0
  for (i = 0; i < sin0.size; ++i)
  {
    printf("%d: %f,%f,%f\n", comm_rank, sin0.data0[3 * i + 0], sin0.data0[3 * i + 1], sin0.data0[3 * i + 2]);
  }
#endif

#ifdef ALLTOALLV_PACKED
  local_packed = ALLTOALLV_PACKED(comm_size, sin0.size);
  MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
  original_packed = fcs_forw_SL_DEFCON(meas.packed); fcs_forw_SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

  old_minalloc = fcs_forw_SL_DEFCON(meas.minalloc);
  fcs_forw_SL_DEFCON(meas.minalloc) = gs->minalloc;
  old_overalloc = fcs_forw_SL_DEFCON(meas.overalloc);
  fcs_forw_SL_DEFCON(meas.overalloc) = gs->overalloc;

#ifndef ALLTOALL_SPECIFIC_IN_PLACE
  TIMING_SYNC(comm); TIMING_START(t[1]);
  fcs_forw_mpi_elements_alltoall_specific(&sin0, &sout0, NULL, tproc, grid_tproc_data, comm_size, comm_rank, comm);
  TIMING_SYNC(comm); TIMING_STOP(t[1]);
#else
  fcs_forw_elements_alloc(&sout0, 10 * sin0.max_size, SLCM_ALL);
  fcs_forw_elements_ncopy(&sin0, &sout0, sin0.size);
  sout0.size = sin0.size;
  TIMING_SYNC(comm); TIMING_START(t[1]);
  fcs_forw_mpi_elements_alltoall_specific(&sout0, NULL, NULL, tproc, grid_tproc_data, comm_size, comm_rank, comm);
  TIMING_SYNC(comm); TIMING_STOP(t[1]);
#endif

#ifdef ALLTOALLV_PACKED
  fcs_forw_SL_DEFCON(meas.packed) = original_packed;
#endif

  fcs_forw_SL_DEFCON(meas.minalloc) = old_minalloc;
  fcs_forw_SL_DEFCON(meas.overalloc) = old_overalloc;

  DEBUG_CMD(
    printf(DEBUG_PRINT_PREFIX "%d: keys out: %" fcs_forw_slint_fmt "\n", comm_rank, sout0.size);
  );

#if 0
  for (i = 0; i < sout0.size; ++i)
  {
/*    if ((sout0.keys[i] >= 0) &&
        (sout0.data0[3 * i + 0] < gs->d.lower_bounds[0] || sout0.data0[3 * i + 0] > gs->d.upper_bounds[0] ||
         sout0.data0[3 * i + 1] < gs->d.lower_bounds[1] || sout0.data0[3 * i + 1] > gs->d.upper_bounds[1] ||
         sout0.data0[3 * i + 2] < gs->d.lower_bounds[2] || sout0.data0[3 * i + 2] > gs->d.upper_bounds[2]))
      printf("%d: %f,%f,%f not in bounds %f,%f,%f - %f,%f,%f\n", comm_rank, sout0.data0[3 * i + 0], sout0.data0[3 * i + 1], sout0.data0[3 * i + 2], gs->d.lower_bounds[0], gs->d.lower_bounds[1], gs->d.lower_bounds[2], gs->d.upper_bounds[0], gs->d.upper_bounds[1], gs->d.upper_bounds[2]);*/
    printf("%d: %f,%f,%f\n", comm_rank, sout0.data0[3 * i + 0], sout0.data0[3 * i + 1], sout0.data0[3 * i + 2]);
  }
#endif

  fcs_forw_tproc_free(&tproc);

  free(original_indices);

#ifdef GRIDSORT_FRONT_TPROC_RANK_CACHE
  free(rank_cache);
#endif

  update_cache(gs);

  release_bounds(&gs->d.bounds);

  gs->nsorted_particles = sout0.size;
  gs->max_nsorted_particles = sout0.max_size;
  gs->sorted_indices = sout0.keys;
  gs->sorted_positions = sout0.data0;
  gs->sorted_charges = sout0.data1;

  gs->nsorted_real_particles = sout0.size;

#ifdef PRINT_FORWARD_SORTED
  print_particles(gs->nsorted_particles, gs->sorted_positions, comm_size, comm_rank, comm);
#endif

  fcs_forw_mpi_datatypes_release();
  
  TIMING_SYNC(comm); TIMING_STOP(t[0]);

  TIMING_CMD(
    if (comm_rank == 0)
      printf(TIMING_PRINT_PREFIX "fcs_gridsort_sort_forward: %f  %f\n", t[0], t[1]);
  );

  return 0;
}

  
#undef BOUNDS_X
#undef BOUNDS_Y
#undef BOUNDS_Z

#undef BOUNDS_XYZ2COORDS


static void fcs_forw_sendrecv(fcs_forw_elements_t *sb, int scount, int sdispl, int dst, fcs_forw_elements_t *rb, int rcount, int rdispl, int src, MPI_Comm comm, int *received)
{
  MPI_Status status;

  MPI_Sendrecv(&sb->keys[sdispl], scount, fcs_forw_sl_key_type_mpi, dst, 0, &rb->keys[rdispl], rcount, fcs_forw_sl_key_type_mpi, src, 0, comm, &status);

  MPI_Get_count(&status, fcs_forw_sl_key_type_mpi, received);

  MPI_Sendrecv(&sb->data0[3 * sdispl], 3 * scount, fcs_forw_sl_data0_type_mpi, dst, 0, &rb->data0[3 * rdispl], 3 * rcount, fcs_forw_sl_data0_type_mpi, src, 0, comm, &status);
  MPI_Sendrecv(&sb->data1[sdispl], scount, fcs_forw_sl_data1_type_mpi, dst, 0, &rb->data1[rdispl], rcount, fcs_forw_sl_data1_type_mpi, src, 0, comm, &status);

}


static void get_neighbors(int *neighbors, int *periodic, int size, int rank, MPI_Comm comm)
{
  int dims[3], periods[3], coords[3];

  MPI_Cart_shift(comm, 0, 1, &neighbors[0], &neighbors[1]);
  MPI_Cart_shift(comm, 1, 1, &neighbors[2], &neighbors[3]);
  MPI_Cart_shift(comm, 2, 1, &neighbors[4], &neighbors[5]);

  MPI_Cart_get(comm, 3, dims, periods, coords);

  periodic[0] = (periods[0] && coords[0] == 0);
  periodic[1] = (periods[0] && coords[0] == dims[0] - 1);
  periodic[2] = (periods[1] && coords[1] == 0);
  periodic[3] = (periods[1] && coords[1] == dims[1] - 1);
  periodic[4] = (periods[2] && coords[2] == 0);
  periodic[5] = (periods[2] && coords[2] == dims[2] - 1);

/*  printf("%d: coords: %d,%d,%d\n", rank, coords[0], coords[1], coords[2]);
  printf("%d: shift: %d,%d  %d,%d  %d,%d\n", rank, neighbors[0], neighbors[1], neighbors[2], neighbors[3], neighbors[4], neighbors[5]);
  printf("%d: periodic: %d,%d  %d,%d  %d,%d\n", rank, periodic[0], periodic[1], periodic[2], periodic[3], periodic[4], periodic[5]);*/
}


static void set_ghosts(fcs_forw_elements_t *s, int count, int displ)
{
  fcs_int i;

  for (i = displ; i < displ + count; ++i)
  {
    if (s->keys[i] >= 0) s->keys[i] = GRIDSORT_GHOST_BASE;
  }
}


static void set_periodics(fcs_forw_elements_t *s, int count, int displ, int d)
{
  fcs_int i;

  for (i = displ; i < displ + count; ++i)
  {
    if (s->keys[i] >= 0) s->keys[i] = GRIDSORT_GHOST_BASE;
    s->keys[i] |= GRIDSORT_PERIODIC_SET(1, d);
  }
}


static int gridsort_split_0b(fcs_forw_elements_t *s, fcs_forw_slint_t x, void *data)
{
  double *bounds = data;
  if (s->data0[3 * x + 0] <= bounds[0])
  {
    if (s->data0[3 * x + 0] >= bounds[1]) return 2;
    return 1;
  }
  if (s->data0[3 * x + 0] >= bounds[1]) return 3;
  return 0;
}

static int gridsort_split_1b(fcs_forw_elements_t *s, fcs_forw_slint_t x, void *data)
{
  double *bounds = data;
  if (s->data0[3 * x + 1] <= bounds[2])
  {
    if (s->data0[3 * x + 1] >= bounds[3]) return 2;
    return 1;
  }
  if (s->data0[3 * x + 1] >= bounds[3]) return 3;
  return 0;
}

static int gridsort_split_2b(fcs_forw_elements_t *s, fcs_forw_slint_t x, void *data)
{
  double *bounds = data;
  if (s->data0[3 * x + 2] <= bounds[4])
  {
    if (s->data0[3 * x + 2] >= bounds[5]) return 2;
    return 1;
  }
  if (s->data0[3 * x + 2] >= bounds[5]) return 3;
  return 0;
}


fcs_int fcs_gridsort_create_ghosts(fcs_gridsort_t *gs, fcs_float ghost_range, MPI_Comm comm)
{
  fcs_int d;
  int comm_size, comm_rank;

  int cart_dims[3], cart_periods[3], cart_coords[3], topo_status;

  fcs_int with_ghost, with_triclinic, with_bounds;

  int counts[4], displs[4], neighbors[6], periodic[6], scounts[2], rcounts[2], received;
  
  fcs_float bounds[6];

  fcs_forw_split_generic_t sg_gridsort_split[] = { fcs_forw_SPLIT_GENERIC_INIT_TPROC(gridsort_split_0b), fcs_forw_SPLIT_GENERIC_INIT_TPROC(gridsort_split_1b), fcs_forw_SPLIT_GENERIC_INIT_TPROC(gridsort_split_2b) };
  
  fcs_forw_elements_t s, sx;

  MPI_Status status;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  MPI_Topo_test(comm, &topo_status);

  if (topo_status != MPI_CART)
  {
    fprintf(stderr, "ERROR: no Cartesian communicator available for gridsort!\n");
    return -1;
  }

  MPI_Cart_get(comm, 3, cart_dims, cart_periods, cart_coords);

  with_ghost = (ghost_range > 0);
  with_triclinic = z_is_triclinic(gs->d.box_a, gs->d.box_b, gs->d.box_c);
  with_bounds = (gs->d.lower_bounds[0] >= 0 && gs->d.lower_bounds[1] >= 0 && gs->d.lower_bounds[2] >= 0 && gs->d.upper_bounds[0] >= 0 && gs->d.upper_bounds[1] >= 0 && gs->d.upper_bounds[2] >= 0);

  INFO_CMD(
    if (comm_rank == 0)
    {
      printf(INFO_PRINT_PREFIX "gridsort create ghosts settings:\n");
      printf(INFO_PRINT_PREFIX " ghost: %s\n", with_ghost?"yes":"no");
      printf(INFO_PRINT_PREFIX " triclinic: %s\n", with_triclinic?"yes":"no");
      printf(INFO_PRINT_PREFIX " bounds: %s\n", with_bounds?"yes":"no");
      printf(INFO_PRINT_PREFIX " cartesian grid:\n");
      printf(INFO_PRINT_PREFIX "  dims: %dx%dx%d\n", cart_dims[0], cart_dims[1], cart_dims[2]);
      printf(INFO_PRINT_PREFIX "  periods: %dx%dx%d\n", cart_periods[0], cart_periods[1], cart_periods[2]);
    }
  );
  
  if (!with_ghost) return 0;

  if (with_triclinic)
  {
    fprintf(stderr, "ERROR: box shape ([%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f] x [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f]"
      " x [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f]) not yet supported for fcs_gridsort_create_ghosts\n",
      gs->d.box_a[0], gs->d.box_a[1], gs->d.box_a[2], gs->d.box_b[0], gs->d.box_b[1], gs->d.box_b[2], gs->d.box_c[0], gs->d.box_c[1], gs->d.box_c[2]);
    return -1;
  }

  get_neighbors(neighbors, periodic, comm_size, comm_rank, comm);

  fcs_forw_elem_set_size(&s, gs->nsorted_particles);
  fcs_forw_elem_set_max_size(&s, gs->max_nsorted_particles);
  fcs_forw_elem_set_keys(&s, gs->sorted_indices);
  fcs_forw_elem_set_data(&s, gs->sorted_positions, gs->sorted_charges);

  fcs_forw_elements_alloc(&sx, 1, SLCM_ALL);

  if (with_bounds)
  {
    bounds[0] = gs->d.lower_bounds[0] + ghost_range;
    bounds[1] = gs->d.upper_bounds[0] - ghost_range;
    bounds[2] = gs->d.lower_bounds[1] + ghost_range;
    bounds[3] = gs->d.upper_bounds[1] - ghost_range;
    bounds[4] = gs->d.lower_bounds[2] + ghost_range;
    bounds[5] = gs->d.upper_bounds[2] - ghost_range;

  } else
  {
    bounds[0] = gs->d.box_base[0] + (gs->d.box_a[0] *  cart_coords[0]      / cart_dims[0]) + ghost_range;
    bounds[1] = gs->d.box_base[0] + (gs->d.box_a[0] * (cart_coords[0] + 1) / cart_dims[0]) - ghost_range;
    bounds[2] = gs->d.box_base[1] + (gs->d.box_b[1] *  cart_coords[1]      / cart_dims[1]) + ghost_range;
    bounds[3] = gs->d.box_base[1] + (gs->d.box_b[1] * (cart_coords[1] + 1) / cart_dims[1]) - ghost_range;
    bounds[4] = gs->d.box_base[2] + (gs->d.box_c[2] *  cart_coords[2]      / cart_dims[2]) + ghost_range;
    bounds[5] = gs->d.box_base[2] + (gs->d.box_c[2] * (cart_coords[2] + 1) / cart_dims[2]) - ghost_range;
  }

  for (d = 0; d < 3; ++d)
  {
    fcs_forw_split_generic_count_ip(&s, &sg_gridsort_split[d], bounds, counts, 4);

    fcs_forw_counts2displs(4, counts, displs);
    fcs_forw_split_generic_rearrange_ip(&s, &sx, &sg_gridsort_split[d], bounds, counts, displs, 4);

    fcs_forw_counts2displs(4, counts, displs);

/*    printf("%d: counts: %d  %d  %d  %d\n", comm_rank, counts[0], counts[1], counts[2], counts[3]);
    printf("%d: displs: %d  %d  %d  %d\n", comm_rank, displs[0], displs[1], displs[2], displs[3]);*/

    scounts[0] = counts[1] + counts[2];
    scounts[1] = counts[2] + counts[3];
    rcounts[0] = rcounts[1] = 0;

    MPI_Sendrecv(&scounts[0], 1, MPI_INT, neighbors[2 * d + 0], 0, &rcounts[0], 1, MPI_INT, neighbors[2 * d + 1], 0, comm, &status);
    MPI_Sendrecv(&scounts[1], 1, MPI_INT, neighbors[2 * d + 1], 0, &rcounts[1], 1, MPI_INT, neighbors[2 * d + 0], 0, comm, &status);

/*    printf("%d: scounts: %d  %d\n", comm_rank, scounts[0], scounts[1]);
    printf("%d: rcounts: %d  %d\n", comm_rank, rcounts[0], rcounts[1]);*/

    fcs_forw_elements_realloc(&s, fcs_forw_elem_get_size(&s) + rcounts[0] + rcounts[1], SLCM_ALL);

    fcs_forw_sendrecv(&s, scounts[0], displs[1], neighbors[2 * d + 0], &s, rcounts[0], s.size, neighbors[2 * d + 1], comm, &received);
    if (periodic[2 * d + 1]) set_periodics(&s, rcounts[0], s.size, 2 * d + 0); else set_ghosts(&s, rcounts[0], s.size);
    s.size += rcounts[0];

    fcs_forw_sendrecv(&s, scounts[1], displs[2], neighbors[2 * d + 1], &s, rcounts[1], s.size, neighbors[2 * d + 0], comm, &received);
    if (periodic[2 * d + 0]) set_periodics(&s, rcounts[1], s.size, 2 * d + 1); else set_ghosts(&s, rcounts[1], s.size);
    s.size += rcounts[1];
  }

  fcs_forw_elements_free(&sx);

  gs->nsorted_particles = s.size;
  gs->max_nsorted_particles = s.max_size;
  gs->sorted_indices = s.keys;
  gs->sorted_positions = s.data0;
  gs->sorted_charges = s.data1;

  gs->nsorted_real_particles = s.size;

  return 0;
}


static void separate_ghosts(fcs_int nparticles, fcs_gridsort_index_t *indices, fcs_float *positions, fcs_float *charges, fcs_int *real_nparticles, fcs_int *ghost_nparticles)
{
  fcs_int l, h;
  fcs_gridsort_index_t ti;
  fcs_float tf;


  l = 0;
  h = nparticles - 1;

  while (1)
  {
    while (l < h)
    if (GRIDSORT_IS_GHOST(indices[l])) break; else ++l;

    while (l < h)
    if (GRIDSORT_IS_GHOST(indices[h])) --h; else break;

    if (l >= h) break;

    z_swap(indices[l], indices[h], ti);
    
    z_swap(positions[3 * l + 0], positions[3 * h + 0], tf);
    z_swap(positions[3 * l + 1], positions[3 * h + 1], tf);
    z_swap(positions[3 * l + 2], positions[3 * h + 2], tf);

    z_swap(charges[l], charges[h], tf);

    ++l;
    --h;
  }

  if (l < nparticles) *real_nparticles = l + ((GRIDSORT_IS_GHOST(indices[l]))?0:1);
  else *real_nparticles = 0;
  *ghost_nparticles = nparticles - *real_nparticles;
}


void fcs_gridsort_separate_ghosts(fcs_gridsort_t *gs, fcs_int *real_nparticles, fcs_int *ghost_nparticles)
{
  separate_ghosts(gs->nsorted_particles, gs->sorted_indices, gs->sorted_positions, gs->sorted_charges, &gs->nsorted_real_particles, &gs->nsorted_ghost_particles);

  if (real_nparticles) *real_nparticles = gs->nsorted_real_particles;
  if (ghost_nparticles) *ghost_nparticles = gs->nsorted_ghost_particles;
}


#define ZSLICES_HEAD \
  fcs_int _i, slice, next, pos, end; \
  fcs_float tf; \
  fcs_gridsort_index_t ti;

#define ZSLICES_COUNT(_n_, _ps_, _nz_, _zc_) \
  for (_i = 0; _i < (_nz_); ++_i) (_zc_)[_i] = 0; \
  for (_i = 0; _i < (_n_); ++_i) { \
    slice = PARTICLE_TO_ZSLICE((_ps_), _i); \
    ++(_zc_)[slice]; \
  }
  
#define ZSLICES_REARRANGE(_n_, _is_, _ps_, _cs_, _nz_, _zc_, _zd_) \
  for (end = 0, _i = 0; _i < (_nz_); ++_i) { \
    pos = (_zd_)[_i]; end = pos + (_zc_)[_i]; \
    while (pos < end) { \
      slice = PARTICLE_TO_ZSLICE((_ps_), pos); \
      while (slice != _i) { \
        next = PARTICLE_TO_ZSLICE((_ps_), (_zd_)[slice]); \
        if (next != slice) { \
          z_swap((_is_)[pos], (_is_)[(_zd_)[slice]], ti); \
          z_swap((_ps_)[3 * pos + 0], (_ps_)[3 * (_zd_)[slice] + 0], tf); \
          z_swap((_ps_)[3 * pos + 1], (_ps_)[3 * (_zd_)[slice] + 1], tf); \
          z_swap((_ps_)[3 * pos + 2], (_ps_)[3 * (_zd_)[slice] + 2], tf); \
          z_swap((_cs_)[pos], (_cs_)[(_zd_)[slice]], tf); \
        } \
        --(_zc_)[slice]; ++(_zd_)[slice]; slice = next; \
      } \
      ++pos; \
    } \
  }


static void separate_zslices(fcs_gridsort_t *gs, fcs_int nparticles, fcs_gridsort_index_t *indices, fcs_float *positions, fcs_float *charges, fcs_int *zslices_low_ghost_nparticles, fcs_int *zslices_real_nparticles, fcs_int *zslices_high_ghost_nparticles)
{
  ZSLICES_HEAD

/*#define SORT_OUT_XY*/
  
#ifdef SORT_OUT_XY
  const fcs_int nz = gs->local_nzslices + (2 * gs->max_ghost_nzslices) + 1;
#else
  const fcs_int nz = gs->local_nzslices + (2 * gs->max_ghost_nzslices);
#endif

  fcs_int counts[nz], displs[nz];
  fcs_float o, f;
  fcs_int i, displ;
#ifdef SORT_OUT_XY
  fcs_float h[2]
#endif

  o = gs->sub_box_base[2];
  f = gs->local_nzslices / gs->sub_box_c[2];

#ifdef SORT_OUT_XY
  h[0] = gs->sub_box_base[0] + gs->sub_box_a[0];
  h[1] = gs->sub_box_base[1] + gs->sub_box_b[1];
# define PARTICLE_TO_ZSLICE(_ps_, _x_) \
  (((_ps_)[3 * (_x_) + 0] < gs->sub_box_base[0] || h[0] <= (_ps_)[3 * (_x_) + 0] || (_ps_)[3 * (_x_) + 1] < gs->sub_box_base[1] || h[1] <= (_ps_)[3 * (_x_) + 1])?(nz - 1):((fcs_int) (((_ps_)[3 * (_x_) + 2] - o) * f + gs->max_ghost_nzslices)))
#else
# define PARTICLE_TO_ZSLICE(_ps_, _x_) \
  ((fcs_int) (((_ps_)[3 * (_x_) + 2] - o) * f + gs->max_ghost_nzslices))
#endif
  ZSLICES_COUNT(nparticles, positions, nz, counts);

  DEBUG_CMD(
    printf(DEBUG_PRINT_PREFIX "counts =");
    for (i = 0; i < nz; ++i) printf(" %" FCS_LMOD_INT "d ", counts[i]);
    printf("\n");
  );

  /* make rearranged displacements (1. low ghost zslices, 2. real zslices, 3. high ghost zslices, 4. low ghost zslices, 5. above high ghost zslices, 6. everything outside x-y-range) */
  displ = 0;
#define SET_COUNT_DISPL(_i_)  Z_MOP(displs[_i_] = displ; displ += counts[_i_];)
  for (i = 0; i < gs->ghost_nzslices; ++i) SET_COUNT_DISPL(gs->max_ghost_nzslices - gs->ghost_nzslices + i);                       /* 1. */
  for (i = 0; i < gs->ghost_nzslices; ++i) SET_COUNT_DISPL(gs->max_ghost_nzslices + gs->local_nzslices + i);                       /* 2. */
  for (i = 0; i < gs->local_nzslices; ++i) SET_COUNT_DISPL(gs->max_ghost_nzslices + i);                                            /* 3. */
  for (i = gs->ghost_nzslices; i < gs->max_ghost_nzslices; ++i) SET_COUNT_DISPL(gs->max_ghost_nzslices - gs->ghost_nzslices + i);  /* 4. */
  for (i = gs->ghost_nzslices; i < gs->max_ghost_nzslices; ++i) SET_COUNT_DISPL(gs->max_ghost_nzslices + gs->local_nzslices + i);  /* 5. */
# ifdef SORT_OUT_XY
  SET_COUNT_DISPL(nz - 1);                                                                                                         /* 6. */
# endif
#undef SET_COUNT_DISPL

  DEBUG_CMD(
    printf(DEBUG_PRINT_PREFIX "displs =");
    for (i = 0; i < nz; ++i) printf(" %" FCS_LMOD_INT "d ", displs[i]);
    printf("\n");

    printf(DEBUG_PRINT_PREFIX "displ = %" FCS_LMOD_INT "d vs. nparticles = %" FCS_LMOD_INT "d\n", displ, nparticles);
  );

  if (zslices_low_ghost_nparticles)
    for (i = 0; i < gs->ghost_nzslices; ++i) { zslices_low_ghost_nparticles[i] = counts[gs->max_ghost_nzslices - gs->ghost_nzslices + i]; }

  if (zslices_real_nparticles)
    for (i = 0; i < gs->local_nzslices; ++i) zslices_real_nparticles[i] = counts[gs->max_ghost_nzslices + i];

  if (zslices_high_ghost_nparticles)
    for (i = 0; i < gs->ghost_nzslices; ++i) { zslices_high_ghost_nparticles[i] = counts[gs->max_ghost_nzslices + gs->local_nzslices + i]; }

  ZSLICES_REARRANGE(nparticles, indices, positions, charges, nz, counts, displs);

#undef PARTICLE_TO_ZSLICE

#if 0
  if (zslices_ghost_nparticles)
  {
    end = 0;
    pos = 0;
    for (i = 0; i < gs->ghost_nzslices; ++i)
    {
      end += zslices_ghost_nparticles[i];
      for (; pos < end; ++pos) printf("%" FCS_LMOD_INT "d  %" FCS_LMOD_FLOAT "f  %" FCS_LMOD_FLOAT "f  %" FCS_LMOD_FLOAT "f\n",
        i, positions[3 * pos + 0], positions[3 * pos + 1], positions[3 * pos + 2]);
    }
  }
#endif
}


void fcs_gridsort_separate_zslices(fcs_gridsort_t *gs, fcs_int *zslices_low_ghost_nparticles, fcs_int *zslices_real_nparticles, fcs_int *zslices_high_ghost_nparticles)
{
  separate_zslices(gs, gs->nsorted_real_particles, gs->sorted_indices, gs->sorted_positions, gs->sorted_charges, zslices_low_ghost_nparticles, zslices_real_nparticles, zslices_high_ghost_nparticles);
  
  if (gs->nsorted_ghost_particles > 0)
  {
    separate_zslices(gs, gs->nsorted_ghost_particles, gs->sorted_indices + gs->nsorted_real_particles, gs->sorted_positions + (3 * gs->nsorted_real_particles), gs->sorted_charges + gs->nsorted_real_particles,
      zslices_low_ghost_nparticles, NULL, zslices_high_ghost_nparticles);
  }
}


static int gridsort_random_tproc(fcs_forw_elements_t *s, fcs_forw_slint_t x, void *data)
{
  int comm_size = ((int *) data)[0];

  return rand() % comm_size;
}


static void gridsort_random_reset(void *data)
{
  unsigned int seed = ((int *) data)[1];

  srand(seed);
}


fcs_int fcs_gridsort_sort_random(fcs_gridsort_t *gs, MPI_Comm comm)
{
  int comm_size, comm_rank;

  fcs_gridsort_index_t *original_indices, index_rank;

  fcs_int i;

  fcs_forw_elements_t sin0, sout0;

  fcs_forw_tproc_t tproc;
  
  int random_data[2];

  fcs_forw_slint_t old_minalloc;
  double old_overalloc;
  
#ifdef DO_TIMING
  double t[2] = { 0, 0 };
#endif
  

  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  gs->nsorted_particles = gs->max_nsorted_particles = 0;
  gs->sorted_positions = NULL;
  gs->sorted_charges = NULL;
  gs->sorted_indices = NULL;

  index_rank = GRIDSORT_INDEX_VAL_PROC(comm_rank);

  fcs_forw_SL_DEFCON(mpi.rank) = comm_rank;

  fcs_forw_mpi_datatypes_init();

  original_indices = malloc(gs->noriginal_particles * sizeof(fcs_gridsort_index_t));

  for (i = 0; i < gs->noriginal_particles; ++i) original_indices[i] = index_rank + i;

  fcs_forw_elem_set_size(&sin0, gs->noriginal_particles);
  fcs_forw_elem_set_max_size(&sin0, gs->noriginal_particles);
  fcs_forw_elem_set_keys(&sin0, original_indices);
  fcs_forw_elem_set_data(&sin0, gs->original_positions, gs->original_charges);

  fcs_forw_elem_set_size(&sout0, 0);
  fcs_forw_elem_set_max_size(&sout0, 0);
  fcs_forw_elem_set_keys(&sout0, NULL);
  fcs_forw_elem_set_data(&sout0, NULL, NULL);

  random_data[0] = comm_size;
  random_data[1] = (comm_rank + 1) * 2501;

  fcs_forw_tproc_create_tproc(&tproc, gridsort_random_tproc, gridsort_random_reset, fcs_forw_TPROC_EXDEF_NULL);

  old_minalloc = fcs_forw_SL_DEFCON(meas.minalloc);
  fcs_forw_SL_DEFCON(meas.minalloc) = gs->minalloc;
  old_overalloc = fcs_forw_SL_DEFCON(meas.overalloc);
  fcs_forw_SL_DEFCON(meas.overalloc) = gs->overalloc;

  TIMING_SYNC(comm); TIMING_START(t[1]);
  fcs_forw_mpi_elements_alltoall_specific(&sin0, &sout0, NULL, tproc, random_data, comm_size, comm_rank, comm);
  TIMING_SYNC(comm); TIMING_STOP(t[1]);

  fcs_forw_SL_DEFCON(meas.minalloc) = old_minalloc;
  fcs_forw_SL_DEFCON(meas.overalloc) = old_overalloc;

  fcs_forw_tproc_free(&tproc);

  free(original_indices);

  gs->nsorted_particles = sout0.size;
  gs->max_nsorted_particles = sout0.max_size;
  gs->sorted_indices = sout0.keys;
  gs->sorted_positions = sout0.data0;
  gs->sorted_charges = sout0.data1;

  gs->nsorted_real_particles = sout0.size;

  fcs_forw_mpi_datatypes_release();

  TIMING_SYNC(comm); TIMING_STOP(t[0]);

  TIMING_CMD(if (comm_rank == 0) printf(TIMING_PRINT_PREFIX "fcs_gridsort_sort_random: %f  %f\n", t[0], t[1]););
  
  return 0;
}


#define V_ADD_ASSIGN(_v0_, _x_, _v1_)  Z_MOP((_v0_)[0] += (_x_) * (_v1_)[0]; (_v0_)[1] += (_x_) * (_v1_)[1]; (_v0_)[2] += (_x_) * (_v1_)[2];)
#define V_SUB_ASSIGN(_v0_, _x_, _v1_)  Z_MOP((_v0_)[0] -= (_x_) * (_v1_)[0]; (_v0_)[1] -= (_x_) * (_v1_)[1]; (_v0_)[2] -= (_x_) * (_v1_)[2];)

void fcs_gridsort_unfold_periodic_particles(fcs_int nparticles, fcs_gridsort_index_t *indices, fcs_float *positions, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c)
{
  fcs_int with_triclinic, i, x;


  with_triclinic = z_is_triclinic(box_a, box_b, box_c);

  if (with_triclinic)
  {
    for (i = 0; i < nparticles; ++i)
    {
      if (GRIDSORT_IS_GHOST(indices[i]))
      {
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 0))) V_ADD_ASSIGN(&positions[3 * i], x, box_a);
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 1))) V_SUB_ASSIGN(&positions[3 * i], x, box_a);
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 2))) V_ADD_ASSIGN(&positions[3 * i], x, box_b);
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 3))) V_SUB_ASSIGN(&positions[3 * i], x, box_b);
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 4))) V_ADD_ASSIGN(&positions[3 * i], x, box_c);
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 5))) V_SUB_ASSIGN(&positions[3 * i], x, box_c);
        
        indices[i] = GRIDSORT_GHOST_BASE;
      }
    }
  } else
  {
    for (i = 0; i < nparticles; ++i)
    {
      if (GRIDSORT_IS_GHOST(indices[i]))
      {
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 0))) positions[3 * i + 0] += x * box_a[0];
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 1))) positions[3 * i + 0] -= x * box_a[0];
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 2))) positions[3 * i + 1] += x * box_b[1];
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 3))) positions[3 * i + 1] -= x * box_b[1];
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 4))) positions[3 * i + 2] += x * box_c[2];
        if ((x = GRIDSORT_PERIODIC_GET(indices[i], 5))) positions[3 * i + 2] -= x * box_c[2];
        
        indices[i] = GRIDSORT_GHOST_BASE;
      }
    }
  }
}


static int gridsort_fcs_back_fp_tproc(fcs_back_fp_elements_t *s, fcs_back_fp_slint_t x, void *data)
{
  if (!GRIDSORT_INDEX_IS_VALID(s->keys[x])) return MPI_PROC_NULL;

  return GRIDSORT_INDEX_GET_PROC(s->keys[x]);
}


static int gridsort_fcs_back_f__tproc(fcs_back_f__elements_t *s, fcs_back_f__slint_t x, void *data)
{
  if (!GRIDSORT_INDEX_IS_VALID(s->keys[x])) return MPI_PROC_NULL;

  return GRIDSORT_INDEX_GET_PROC(s->keys[x]);
}


static int gridsort_fcs_back__p_tproc(fcs_back__p_elements_t *s, fcs_back__p_slint_t x, void *data)
{
  if (!GRIDSORT_INDEX_IS_VALID(s->keys[x])) return MPI_PROC_NULL;

  return GRIDSORT_INDEX_GET_PROC(s->keys[x]);
}


fcs_int fcs_gridsort_sort_backward(fcs_gridsort_t *gs,
                                   fcs_float *sorted_field, fcs_float *sorted_potentials,
                                   fcs_float *original_field, fcs_float *original_potentials,
                                   fcs_int set_values,
                                   MPI_Comm comm)
{
  int comm_size, comm_rank;

  fcs_int i, j, type;

  fcs_back_fp_elements_t sin0, sout0;
  fcs_back_f__elements_t sin1, sout1;
  fcs_back__p_elements_t sin2, sout2;
  
  fcs_back_fp_tproc_t tproc0;
  fcs_back_f__tproc_t tproc1;
  fcs_back__p_tproc_t tproc2;

  const fcs_gridsort_index_t index_mask = 0x00000000FFFFFFFFLL;

#ifdef ALLTOALLV_PACKED
  fcs_int local_packed, global_packed, original_packed;
#endif

#ifdef DO_TIMING
  double t[2] = { 0, 0 };
#endif
  

  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if (original_field && original_potentials) type = 0;
  else if (original_field && !original_potentials) type = 1;
  else type = 2;

  switch (type)
  {
    case 0:
      fcs_back_fp_SL_DEFCON(mpi.rank) = comm_rank;

      fcs_back_fp_mpi_datatypes_init();

      fcs_back_fp_elem_set_size(&sin0, gs->nsorted_real_particles);
      fcs_back_fp_elem_set_max_size(&sin0, gs->nsorted_real_particles);
      fcs_back_fp_elem_set_keys(&sin0, gs->sorted_indices);
      fcs_back_fp_elem_set_data(&sin0, sorted_field, sorted_potentials);

      fcs_back_fp_elem_set_size(&sout0, 0);
      fcs_back_fp_elem_set_max_size(&sout0, 0);
      fcs_back_fp_elem_set_keys(&sout0, NULL);
      fcs_back_fp_elem_set_data(&sout0, NULL, NULL);

      fcs_back_fp_tproc_create_tproc(&tproc0, gridsort_fcs_back_fp_tproc, fcs_back_fp_TPROC_RESET_NULL, fcs_back_fp_TPROC_EXDEF_NULL);

#ifdef GRIDSORT_BACK_PROCLIST
      if (gs->procs) fcs_back_fp_tproc_set_proclists(&tproc0, gs->nprocs, gs->procs, gs->nprocs, gs->procs, comm_size, comm_rank, comm);
#endif

#ifdef ALLTOALLV_PACKED
      local_packed = ALLTOALLV_PACKED(comm_size, sin0.size);
      MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
      original_packed = fcs_back_fp_SL_DEFCON(meas.packed); fcs_back_fp_SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

      TIMING_SYNC(comm); TIMING_START(t[1]);
      fcs_back_fp_mpi_elements_alltoall_specific(&sin0, &sout0, NULL, tproc0, NULL, comm_size, comm_rank, comm);
      TIMING_SYNC(comm); TIMING_STOP(t[1]);

#ifdef ALLTOALLV_PACKED
      fcs_back_fp_SL_DEFCON(meas.packed) = original_packed;
#endif

      fcs_back_fp_tproc_free(&tproc0);

      if (gs->noriginal_particles != sout0.size)
        fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" fcs_back_fp_slint_fmt "!\n", comm_rank, gs->noriginal_particles, sout0.size);

      if (set_values)
      {
        for (i = 0; i < sout0.size; ++i)
        {
          j = sout0.keys[i] & index_mask;

          original_field[3 * j + 0] = sout0.data0[3 * i + 0];
          original_field[3 * j + 1] = sout0.data0[3 * i + 1];
          original_field[3 * j + 2] = sout0.data0[3 * i + 2];

          original_potentials[j] = sout0.data1[i];
        }
      } else
      {
        for (i = 0; i < sout0.size; ++i)
        {
          j = sout0.keys[i] & index_mask;

          original_field[3 * j + 0] += sout0.data0[3 * i + 0];
          original_field[3 * j + 1] += sout0.data0[3 * i + 1];
          original_field[3 * j + 2] += sout0.data0[3 * i + 2];

          original_potentials[j] += sout0.data1[i];
        }
      }

      fcs_back_fp_elements_free(&sout0);

      fcs_back_fp_mpi_datatypes_release();

      break;

    case 1:
      fcs_back_f__SL_DEFCON(mpi.rank) = comm_rank;

      fcs_back_f__mpi_datatypes_init();

      fcs_back_f__elem_set_size(&sin1, gs->nsorted_real_particles);
      fcs_back_f__elem_set_max_size(&sin1, gs->nsorted_real_particles);
      fcs_back_f__elem_set_keys(&sin1, gs->sorted_indices);
      fcs_back_f__elem_set_data(&sin1, sorted_field);

      fcs_back_f__elem_set_size(&sout1, 0);
      fcs_back_f__elem_set_max_size(&sout1, 0);
      fcs_back_f__elem_set_keys(&sout1, NULL);
      fcs_back_f__elem_set_data(&sout1, NULL);

      fcs_back_f__tproc_create_tproc(&tproc1, gridsort_fcs_back_f__tproc, fcs_back_f__TPROC_RESET_NULL, fcs_back_f__TPROC_EXDEF_NULL);

#ifdef GRIDSORT_BACK_PROCLIST
      if (gs->procs) fcs_back_f__tproc_set_proclists(&tproc1, gs->nprocs, gs->procs, gs->nprocs, gs->procs, comm_size, comm_rank, comm);
#endif

#ifdef ALLTOALLV_PACKED
      local_packed = ALLTOALLV_PACKED(comm_size, sin1.size);
      MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
      original_packed = fcs_back_f__SL_DEFCON(meas.packed); fcs_back_f__SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

      TIMING_SYNC(comm); TIMING_START(t[1]);
      fcs_back_f__mpi_elements_alltoall_specific(&sin1, &sout1, NULL, tproc1, NULL, comm_size, comm_rank, comm);
      TIMING_SYNC(comm); TIMING_STOP(t[1]);

#ifdef ALLTOALLV_PACKED
      fcs_back_f__SL_DEFCON(meas.packed) = original_packed;
#endif

      fcs_back_f__tproc_free(&tproc1);

      if (gs->noriginal_particles != sout1.size)
        fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" fcs_back_f__slint_fmt "!\n", comm_rank, gs->noriginal_particles, sout1.size);

      if (set_values)
      {
        for (i = 0; i < sout1.size; ++i)
        {
          j = sout1.keys[i] & index_mask;
     
          original_field[3 * j + 0] = sout1.data0[3 * i + 0];
          original_field[3 * j + 1] = sout1.data0[3 * i + 1];
          original_field[3 * j + 2] = sout1.data0[3 * i + 2];
        }

      } else
      {
        for (i = 0; i < sout1.size; ++i)
        {
          j = sout1.keys[i] & index_mask;
     
          original_field[3 * j + 0] += sout1.data0[3 * i + 0];
          original_field[3 * j + 1] += sout1.data0[3 * i + 1];
          original_field[3 * j + 2] += sout1.data0[3 * i + 2];
        }
      }

      fcs_back_f__elements_free(&sout1);

      fcs_back_f__mpi_datatypes_release();

      break;

    case 2:
      fcs_back__p_SL_DEFCON(mpi.rank) = comm_rank;

      fcs_back__p_mpi_datatypes_init();

      fcs_back__p_elem_set_size(&sin2, gs->nsorted_real_particles);
      fcs_back__p_elem_set_max_size(&sin2, gs->nsorted_real_particles);
      fcs_back__p_elem_set_keys(&sin2, gs->sorted_indices);
      fcs_back__p_elem_set_data(&sin2, sorted_potentials);

      fcs_back__p_elem_set_size(&sout2, 0);
      fcs_back__p_elem_set_max_size(&sout2, 0);
      fcs_back__p_elem_set_keys(&sout2, NULL);
      fcs_back__p_elem_set_data(&sout2, NULL);

      fcs_back__p_tproc_create_tproc(&tproc2, gridsort_fcs_back__p_tproc, fcs_back__p_TPROC_RESET_NULL, fcs_back__p_TPROC_EXDEF_NULL);

#ifdef GRIDSORT_BACK_PROCLIST
      if (gs->procs) fcs_back__p_tproc_set_proclists(&tproc2, gs->nprocs, gs->procs, gs->nprocs, gs->procs, comm_size, comm_rank, comm);
#endif

#ifdef ALLTOALLV_PACKED
      local_packed = ALLTOALLV_PACKED(comm_size, sin2.size);
      MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
      original_packed = fcs_back__p_SL_DEFCON(meas.packed); fcs_back__p_SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

      TIMING_SYNC(comm); TIMING_START(t[1]);
      fcs_back__p_mpi_elements_alltoall_specific(&sin2, &sout2, NULL, tproc2, NULL, comm_size, comm_rank, comm);
      TIMING_SYNC(comm); TIMING_STOP(t[1]);

#ifdef ALLTOALLV_PACKED
      fcs_back__p_SL_DEFCON(meas.packed) = original_packed;
#endif

      fcs_back__p_tproc_free(&tproc2);

      if (gs->noriginal_particles != sout2.size)
        fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" fcs_back__p_slint_fmt "!\n", comm_rank, gs->noriginal_particles, sout2.size);

      if (set_values)
      {
        for (i = 0; i < sout2.size; ++i)
        {
          j = sout2.keys[i] & index_mask;
     
          original_potentials[j] = sout2.data1[i];
        }
      } else
      {
        for (i = 0; i < sout2.size; ++i)
        {
          j = sout2.keys[i] & index_mask;

          original_potentials[j] += sout2.data1[i];
        }
      }

      fcs_back__p_mpi_datatypes_release();

      break;
  }

  TIMING_SYNC(comm); TIMING_STOP(t[0]);

  TIMING_CMD(
    if (comm_rank == 0)
      printf(TIMING_PRINT_PREFIX "fcs_gridsort_sort_backward: %f  %f\n", t[0], t[1]);
  );

  return 0;
}


fcs_int fcs_gridsort_prepare_resort(fcs_gridsort_t *gs,
                                    fcs_float *sorted_field, fcs_float *sorted_potentials,
                                    fcs_float *original_field, fcs_float *original_potentials,
                                    MPI_Comm comm)
{
  int comm_size, comm_rank;
  fcs_int i, j, nresort_particles;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  
  if (sorted_field == NULL && sorted_potentials == NULL)
  {
    gs->nresort_particles = gs->nsorted_real_particles;
    return 1;
  }

  for (j = 0, i = 0; i < gs->nsorted_real_particles; ++i)
  if (gs->sorted_indices[i] >= 0) ++j;

  nresort_particles = j;

  i = (nresort_particles > gs->max_noriginal_particles)?1:0;
  MPI_Allreduce(&i, &j, 1, FCS_MPI_INT, MPI_SUM, comm);

  DEBUG_CMD(
    if (i) printf(DEBUG_PRINT_PREFIX "%d: Resort disabled, because max. array size %" FCS_LMOD_INT "d too small to store %" FCS_LMOD_INT "d particles\n", comm_rank, gs->max_noriginal_particles, nresort_particles);
  );

  if (j > 0)
  {
    INFO_CMD(
      if (comm_rank == 0) printf(INFO_PRINT_PREFIX "Resort disabled, because max. array sizes on %" FCS_LMOD_INT "d of %d process(es) to small!\n", j, comm_size);
    );
    return 0;
  }

  j = 0;
  for (i = 0; i < gs->nsorted_real_particles; ++i)
  {
/*    printf("sorted_indices[%" FCS_LMOD_INT "d]: " GRIDSORT_INDEX_STR "\n", i, GRIDSORT_INDEX_PARAM(gs->sorted_indices[i]));*/

    if (gs->sorted_indices[i] < 0) continue;

    gs->sorted_indices[j] = gs->sorted_indices[i];

    gs->original_positions[3 * j + 0] = gs->sorted_positions[3 * i + 0];
    gs->original_positions[3 * j + 1] = gs->sorted_positions[3 * i + 1];
    gs->original_positions[3 * j + 2] = gs->sorted_positions[3 * i + 2];

    gs->original_charges[j] = gs->sorted_charges[i];

    original_field[3 * j + 0] = sorted_field[3 * i + 0];
    original_field[3 * j + 1] = sorted_field[3 * i + 1];
    original_field[3 * j + 2] = sorted_field[3 * i + 2];

    original_potentials[j] = sorted_potentials[i];

    ++j;
  }

  gs->nresort_particles = j;

  return 1;
}


void fcs_gridsort_free(fcs_gridsort_t *gs)
{
  fcs_forw_elements_t s0;


  fcs_forw_elem_set_keys(&s0, gs->sorted_indices);
  fcs_forw_elem_set_data(&s0, gs->sorted_positions, gs->sorted_charges);

  fcs_forw_elements_free(&s0);

  gs->nsorted_particles = gs->max_nsorted_particles = 0;
  gs->sorted_indices = NULL;
  gs->sorted_positions = NULL;
  gs->sorted_charges = NULL;
}
