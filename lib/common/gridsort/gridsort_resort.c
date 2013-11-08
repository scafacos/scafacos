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

#include "sl_back_f_.h"
#include "sl_back__p.h"
#include "sl_back_x.h"

#ifdef HAVE_ZMPI_ATASP_H
# include "zmpi_atasp.h"
#endif

#include "z_tools.h"
#include "common.h"
#include "gridsort.h"
#include "gridsort_resort.h"


/*#define RESORT_13FLOATS*/


#ifdef RESORT_13FLOATS

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

#endif /* RESORT_13FLOATS */


static int gridsort_fcs_back_x_tproc(fcs_back_x_elements_t *s, fcs_back_x_slint_t x, void *data)
{
  if (!GRIDSORT_INDEX_IS_VALID(s->keys[x])) return MPI_PROC_NULL;

  return GRIDSORT_INDEX_GET_PROC(s->keys[x]);
}


void fcs_gridsort_resort_create(fcs_gridsort_resort_t *gridsort_resort, fcs_gridsort_t *gs, MPI_Comm comm)
{
  int comm_size, comm_rank;

  fcs_int i;

  fcs_back_x_elements_t sin, sout;

  fcs_back_x_tproc_t tproc;

  fcs_resort_index_t *resort_indices;
  
#ifdef ALLTOALLV_PACKED
  fcs_int local_packed, global_packed, original_packed;
#endif

#ifdef DO_TIMING
  double t[4] = { 0, 0, 0, 0 };
#endif


  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  fcs_resort_create(gridsort_resort);

  fcs_resort_set_original_particles(*gridsort_resort, gs->noriginal_particles);
  fcs_resort_set_sorted_particles(*gridsort_resort, gs->nresort_particles);

  if (gs->nresort_particles < 0)
  {
    fcs_resort_set_sorted_particles(*gridsort_resort, gs->noriginal_particles);
    return;
  }

  TIMING_SYNC(comm); TIMING_START(t[1]);

  resort_indices = fcs_resort_indices_alloc(gs->nresort_particles);

  fcs_resort_indices_init(gs->nresort_particles, resort_indices, comm_rank);

  TIMING_SYNC(comm); TIMING_STOP(t[1]);

/*  printf("nresort_particles = %" FCS_LMOD_INT "d\n", gs->nresort_particles);
  for (i = 0; i < gs->nresort_particles; ++i)
  {
    printf(" %" FCS_LMOD_INT "d: " GRIDSORT_INDEX_STR "  " FCS_RESORT_INDEX_STR "\n",
      i, GRIDSORT_INDEX_PARAM(gs->sorted_indices[i]), FCS_RESORT_INDEX_PARAM(resort_indices[i]));
  }*/
  
  fcs_back_x_SL_DEFCON(mpi.rank) = comm_rank;

  fcs_back_x_mpi_datatypes_init();

  fcs_back_x_elem_set_size(&sin, gs->nresort_particles);
  fcs_back_x_elem_set_max_size(&sin, gs->nresort_particles);
  fcs_back_x_elem_set_keys(&sin, gs->sorted_indices);
  fcs_back_x_elem_set_data(&sin, resort_indices);

  fcs_back_x_elem_set_size(&sout, 0);
  fcs_back_x_elem_set_max_size(&sout, 0);
  fcs_back_x_elem_set_keys(&sout, NULL);
  fcs_back_x_elem_set_data(&sout, NULL);

  fcs_back_x_tproc_create_tproc(&tproc, gridsort_fcs_back_x_tproc, fcs_back_x_TPROC_RESET_NULL, fcs_back_x_TPROC_EXDEF_NULL);

#ifdef GRIDSORT_RESORT_PROCLIST
  if (gs->procs) fcs_back_x_tproc_set_proclists(&tproc, gs->nprocs, gs->procs, gs->nprocs, gs->procs, comm_size, comm_rank, comm);
#endif

#ifdef ALLTOALLV_PACKED
  local_packed = ALLTOALLV_PACKED(comm_size, sin.size);
  MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
  original_packed = fcs_back_x_SL_DEFCON(meas.packed); fcs_back_x_SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

  TIMING_SYNC(comm); TIMING_START(t[2]);

  fcs_back_x_mpi_elements_alltoall_specific(&sin, &sout, NULL, tproc, NULL, comm_size, comm_rank, comm);

  TIMING_SYNC(comm); TIMING_STOP(t[2]);

#ifdef ALLTOALLV_PACKED
  fcs_back_x_SL_DEFCON(meas.packed) = original_packed;
#endif

  fcs_back_x_tproc_free(&tproc);

  fcs_resort_indices_free(resort_indices);

  if (gs->noriginal_particles != sout.size)
    fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" fcs_back_x_slint_fmt "!\n", comm_rank, gs->noriginal_particles, sout.size);

  fcs_back_x_mpi_datatypes_release();

/*  printf("noriginal_particles = %" fcs_back_x_slint_fmt "\n", sout.size);
  for (i = 0; i < sout.size; ++i)
  {
    printf(" %" FCS_LMOD_INT "d: " GRIDSORT_INDEX_STR "  " FCS_RESORT_INDEX_STR "\n",
      i, GRIDSORT_INDEX_PARAM(sout.keys[i]), FCS_RESORT_INDEX_PARAM(sout.data0[i]));
  }*/

  fcs_resort_alloc_indices(*gridsort_resort);
  
  resort_indices = fcs_resort_get_indices(*gridsort_resort);

  TIMING_SYNC(comm); TIMING_START(t[3]);

  for (i = 0; i < gs->noriginal_particles; ++i) resort_indices[GRIDSORT_INDEX_GET_POS(sout.keys[i])] = sout.data0[i];

  TIMING_SYNC(comm); TIMING_STOP(t[3]);

  fcs_back_x_elements_free(&sout);

#ifdef GRIDSORT_RESORT_PROCLIST
  if (gs->procs) fcs_resort_set_proclists(*gridsort_resort, gs->nprocs, gs->procs);
#endif

  TIMING_SYNC(comm); TIMING_STOP(t[0]);

  TIMING_CMD(
    if (comm_rank == 0)
      printf(TIMING_PRINT_PREFIX "fcs_gridsort_resort_create: %f  %f  %f  %f\n", t[0], t[1], t[2], t[3]);
  );
}


void fcs_gridsort_resort_destroy(fcs_gridsort_resort_t *gridsort_resort)
{
  fcs_resort_destroy(gridsort_resort);
}


void fcs_gridsort_resort_print(fcs_gridsort_resort_t gridsort_resort, MPI_Comm comm)
{
  printf("RESORT:\n");
  fcs_resort_print(gridsort_resort);
}


fcs_int fcs_gridsort_resort_is_available(fcs_gridsort_resort_t gridsort_resort)
{
  return fcs_resort_is_available(gridsort_resort);
}


fcs_int fcs_gridsort_resort_get_original_particles(fcs_gridsort_resort_t gridsort_resort)
{
  return fcs_resort_get_original_particles(gridsort_resort);
}


fcs_int fcs_gridsort_resort_get_sorted_particles(fcs_gridsort_resort_t gridsort_resort)
{
  return fcs_resort_get_sorted_particles(gridsort_resort);
}


#ifdef RESORT_13FLOATS

static void copy_floats(fcs_int n, fcs_float *src, fcs_float *dst, fcs_int x, fcs_gridsort_index_t *indices)
{
  fcs_int i, j, k;


  if (indices)
  {
    for (i = 0; i < n; ++i)
    {
      k = GRIDSORT_INDEX_GET_POS(indices[i]);

      for (j = 0; j < x; ++j) dst[k * x + j] = src[i * x + j];
    }

  } else
  {
    memcpy(dst, src, n * x * sizeof(fcs_float));
  }
}


static void resort_1float(fcs_gridsort_resort_t gridsort_resort, fcs_float *src, fcs_float *dst, MPI_Comm comm)
{
  int comm_size, comm_rank;
  fcs_back__p_elements_t sin, sout;
  fcs_back__p_tproc_t tproc;

#ifdef ALLTOALLV_PACKED
  fcs_int local_packed, global_packed, original_packed;
#endif

#ifdef DO_TIMING
  double t[3] = { 0, 0, 0 };
#endif


  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  fcs_back__p_SL_DEFCON(mpi.rank) = comm_rank;

  fcs_back__p_mpi_datatypes_init();

  fcs_back__p_elem_set_size(&sin, fcs_resort_get_original_particles(gridsort_resort));
  fcs_back__p_elem_set_max_size(&sin, fcs_resort_get_original_particles(gridsort_resort));
  fcs_back__p_elem_set_keys(&sin, fcs_resort_get_indices(gridsort_resort));
  fcs_back__p_elem_set_data(&sin, src);

  fcs_back__p_elem_set_size(&sout, 0);
  fcs_back__p_elem_set_max_size(&sout, 0);
  fcs_back__p_elem_set_keys(&sout, NULL);
  fcs_back__p_elem_set_data(&sout, NULL);

  fcs_back__p_tproc_create_tproc(&tproc, gridsort_fcs_back__p_tproc, fcs_back__p_TPROC_RESET_NULL, fcs_back__p_TPROC_EXDEF_NULL);

#ifdef GRIDSORT_RESORT_PROCLIST
  if (gridsort_resort->nprocs >= 0) fcs_back__p_tproc_set_proclists(&tproc, gridsort_resort->nprocs, gridsort_resort->procs, gridsort_resort->nprocs, gridsort_resort->procs, comm_size, comm_rank, comm);
#endif

#ifdef ALLTOALLV_PACKED
  local_packed = ALLTOALLV_PACKED(comm_size, sin.size);
  MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
  original_packed = fcs_back__p_SL_DEFCON(meas.packed); fcs_back__p_SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

  TIMING_SYNC(comm); TIMING_START(t[1]);

  fcs_back__p_mpi_elements_alltoall_specific(&sin, &sout, NULL, tproc, NULL, comm_size, comm_rank, comm);

  TIMING_SYNC(comm); TIMING_STOP(t[1]);

#ifdef ALLTOALLV_PACKED
  fcs_back__p_SL_DEFCON(meas.packed) = original_packed;
#endif

  fcs_back__p_tproc_free(&tproc);

  if (gridsort_resort->nsorted_particles != sout.size)
    fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" fcs_back__p_slint_fmt "!\n", comm_rank, gridsort_resort->nsorted_particles, sout.size);

  TIMING_SYNC(comm); TIMING_START(t[2]);

  copy_floats(sout.size, sout.data1, (dst)?dst:src, 1, sout.keys);

  TIMING_SYNC(comm); TIMING_STOP(t[2]);

  fcs_back__p_elements_free(&sout);

  fcs_back__p_mpi_datatypes_release();

  TIMING_SYNC(comm); TIMING_STOP(t[0]);

  TIMING_CMD(
    if (comm_rank == 0)
      printf(TIMING_PRINT_PREFIX "resort_1float: %f  %f  %f\n", t[0], t[1], t[2]);
  );
}


static void resort_3floats(fcs_gridsort_resort_t gridsort_resort, fcs_float *src, fcs_float *dst, MPI_Comm comm)
{
  int comm_size, comm_rank;
  fcs_back_f__elements_t sin, sout;
  fcs_back_f__tproc_t tproc;

#ifdef ALLTOALLV_PACKED
  fcs_int local_packed, global_packed, original_packed;
#endif

#ifdef DO_TIMING
  double t[3] = { 0, 0, 0 };
#endif


  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  fcs_back_f__SL_DEFCON(mpi.rank) = comm_rank;

  fcs_back_f__mpi_datatypes_init();

  fcs_back__p_elem_set_size(&sin, fcs_resort_get_original_particles(gridsort_resort));
  fcs_back__p_elem_set_max_size(&sin, fcs_resort_get_original_particles(gridsort_resort));
  fcs_back__p_elem_set_keys(&sin, fcs_resort_get_indices(gridsort_resort));
  fcs_back_f__elem_set_data(&sin, src);

  fcs_back_f__elem_set_size(&sout, 0);
  fcs_back_f__elem_set_max_size(&sout, 0);
  fcs_back_f__elem_set_keys(&sout, NULL);
  fcs_back_f__elem_set_data(&sout, NULL);

  fcs_back_f__tproc_create_tproc(&tproc, gridsort_fcs_back_f__tproc, fcs_back_f__TPROC_RESET_NULL, fcs_back_f__TPROC_EXDEF_NULL);

#ifdef GRIDSORT_RESORT_PROCLIST
  if (gridsort_resort->nprocs >= 0) fcs_back_f__tproc_set_proclists(&tproc, gridsort_resort->nprocs, gridsort_resort->procs, gridsort_resort->nprocs, gridsort_resort->procs, comm_size, comm_rank, comm);
#endif

#ifdef ALLTOALLV_PACKED
  local_packed = ALLTOALLV_PACKED(comm_size, sin.size);
  MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
  original_packed = fcs_back_f__SL_DEFCON(meas.packed); fcs_back_f__SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

  TIMING_SYNC(comm); TIMING_START(t[1]);

  fcs_back_f__mpi_elements_alltoall_specific(&sin, &sout, NULL, tproc, NULL, comm_size, comm_rank, comm);

  TIMING_SYNC(comm); TIMING_STOP(t[1]);

#ifdef ALLTOALLV_PACKED
  fcs_back_f__SL_DEFCON(meas.packed) = original_packed;
#endif

  fcs_back_f__tproc_free(&tproc);

  if (gridsort_resort->nsorted_particles != sout.size)
    fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" fcs_back_f__slint_fmt "!\n", comm_rank, gridsort_resort->nsorted_particles, sout.size);

  TIMING_SYNC(comm); TIMING_START(t[2]);

  copy_floats(sout.size, sout.data0, (dst)?dst:src, 3, sout.keys);

  TIMING_SYNC(comm); TIMING_STOP(t[2]);

  fcs_back_f__elements_free(&sout);

  fcs_back_f__mpi_datatypes_release();

  TIMING_SYNC(comm); TIMING_STOP(t[0]);

  TIMING_CMD(
    if (comm_rank == 0)
      printf(TIMING_PRINT_PREFIX "resort_3float: %f  %f  %f\n", t[0], t[1], t[2]);
  );
}

#endif /* RESORT_13FLOATS */


void fcs_gridsort_resort_ints(fcs_gridsort_resort_t gridsort_resort, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  fcs_resort_resort_ints(gridsort_resort, src, dst, n, comm);
}


void fcs_gridsort_resort_floats(fcs_gridsort_resort_t gridsort_resort, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
#ifdef RESORT_13FLOATS
  if (fcs_resort_get_indices(gridsort_resort) != NULL)
  switch (n)
  {
    case 1:
      resort_1float(gridsort_resort, src, dst, comm);
      return;
    case 3:
      resort_3floats(gridsort_resort, src, dst, comm);
      return;
  }
#endif

  fcs_resort_resort_floats(gridsort_resort, src, dst, n, comm);
}


void fcs_gridsort_resort_bytes(fcs_gridsort_resort_t gridsort_resort, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  fcs_resort_resort_bytes(gridsort_resort, src, dst, n, comm);
}
