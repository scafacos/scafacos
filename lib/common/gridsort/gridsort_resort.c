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


static int gridsort_back_f__tproc(back_f__elements_t *s, back_f__slint_t x, void *data)
{
  if (!GRIDSORT_INDEX_IS_VALID(s->keys[x])) return MPI_PROC_NULL;

  return GRIDSORT_INDEX_GET_PROC(s->keys[x]);
}


static int gridsort_back__p_tproc(back__p_elements_t *s, back__p_slint_t x, void *data)
{
  if (!GRIDSORT_INDEX_IS_VALID(s->keys[x])) return MPI_PROC_NULL;

  return GRIDSORT_INDEX_GET_PROC(s->keys[x]);
}


static int gridsort_back_x_tproc(back_x_elements_t *s, back_x_slint_t x, void *data)
{
  if (!GRIDSORT_INDEX_IS_VALID(s->keys[x])) return MPI_PROC_NULL;

  return GRIDSORT_INDEX_GET_PROC(s->keys[x]);
}


void fcs_gridsort_resort_create(fcs_gridsort_resort_t *gsr, fcs_gridsort_t *gs, MPI_Comm comm)
{
  int comm_size, comm_rank;

  fcs_int i;

  back_x_elements_t sin, sout;

  back_x_tproc_t tproc;

  fcs_gridsort_index_t *sorted_indices2, index_rank;
  
#ifdef ALLTOALLV_PACKED
  fcs_int local_packed, global_packed, original_packed;
#endif


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  *gsr = malloc(sizeof(**gsr));
  
  if (gs->nresort_particles < 0)
  {
    (*gsr)->noriginal_particles = gs->noriginal_particles;
    (*gsr)->nsorted_particles = gs->nresort_particles;
    (*gsr)->indices = NULL;

    return;
  }

  sorted_indices2 = malloc(gs->nresort_particles * sizeof(fcs_gridsort_index_t));

  index_rank = GRIDSORT_INDEX_VAL_PROC(comm_rank);

  for (i = 0; i < gs->nresort_particles; ++i) sorted_indices2[i] = index_rank + i;

/*  printf("nresort_particles = %" FCS_LMOD_INT "d\n", gs->nresort_particles);
  for (i = 0; i < gs->nresort_particles; ++i)
  {
    printf(" %" FCS_LMOD_INT "d: " GRIDSORT_INDEX_STR "  " GRIDSORT_INDEX_STR "\n", i, GRIDSORT_INDEX_PARAM(gs->sorted_indices[i]), GRIDSORT_INDEX_PARAM(sorted_indices2[i]));
  }*/
  
  back_x_SL_DEFCON(mpi.rank) = comm_rank;

  back_x_mpi_datatypes_init();

  back_x_elem_set_size(&sin, gs->nresort_particles);
  back_x_elem_set_max_size(&sin, gs->nresort_particles);
  back_x_elem_set_keys(&sin, gs->sorted_indices);
  back_x_elem_set_data(&sin, sorted_indices2);

  back_x_elem_set_size(&sout, 0);
  back_x_elem_set_max_size(&sout, 0);
  back_x_elem_set_keys(&sout, NULL);
  back_x_elem_set_data(&sout, NULL);

  back_x_tproc_create_tproc(&tproc, gridsort_back_x_tproc, back_x_TPROC_RESET_NULL, back_x_TPROC_EXDEF_NULL);

#ifdef ALLTOALLV_PACKED
  local_packed = ALLTOALLV_PACKED(comm_size, sin.size);
  MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
  original_packed = back_x_SL_DEFCON(meas.packed); back_x_SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

  TIMING_SYNC(comm); TIMING_START(t[1]);
  back_x_mpi_elements_alltoall_specific(&sin, &sout, NULL, tproc, NULL, comm_size, comm_rank, comm);
  TIMING_SYNC(comm); TIMING_STOP(t[1]);

#ifdef ALLTOALLV_PACKED
  back_x_SL_DEFCON(meas.packed) = original_packed;
#endif

  back_x_tproc_free(&tproc);

  if (gs->noriginal_particles != sout.size)
    fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" back_x_slint_fmt "!\n", comm_rank, gs->noriginal_particles, sout.size);

  back_x_mpi_datatypes_release();

/*  printf("noriginal_particles = %" back_x_slint_fmt "\n", sout.size);
  for (i = 0; i < sout.size; ++i)
  {
    printf(" %" FCS_LMOD_INT "d: " GRIDSORT_INDEX_STR "  " GRIDSORT_INDEX_STR "\n",
      i, GRIDSORT_INDEX_PARAM(sout.keys[i]), GRIDSORT_INDEX_PARAM(sout.data0[i]));
  }*/

  (*gsr)->noriginal_particles = gs->noriginal_particles;
  (*gsr)->nsorted_particles = gs->nresort_particles;
  (*gsr)->indices = malloc(gs->noriginal_particles * sizeof(fcs_gridsort_index_t));

  for (i = 0; i < gs->noriginal_particles; ++i) (*gsr)->indices[GRIDSORT_INDEX_GET_POS(sout.keys[i])] = sout.data0[i];
  
  back_x_elements_free(&sout);
  
  free(sorted_indices2);
}


void fcs_gridsort_resort_destroy(fcs_gridsort_resort_t *gsr)
{
  back_x_elements_t sout;


  if (gsr == FCS_GRIDSORT_RESORT_NULL) return;

  back_x_elem_set_size(&sout, 0);
  back_x_elem_set_max_size(&sout, 0);
  back_x_elem_set_keys(&sout, NULL);
  back_x_elem_set_data(&sout, (*gsr)->indices);

  back_x_elements_free(&sout);

  free(*gsr);

  *gsr = FCS_GRIDSORT_RESORT_NULL;
}


void fcs_gridsort_resort_print(fcs_gridsort_resort_t gsr, MPI_Comm comm)
{
  fcs_int i;


  printf("RESORT:\n");

  if (gsr->indices)
    for (i = 0; i < gsr->noriginal_particles; ++i)
      printf("%" FCS_LMOD_INT "d: " GRIDSORT_INDEX_STR "\n", i,  GRIDSORT_INDEX_PARAM(gsr->indices[i]));
  else
    printf("-> ID\n");
}


fcs_int fcs_gridsort_resort_is_available(fcs_gridsort_resort_t gsr)
{
  return (gsr != FCS_GRIDSORT_RESORT_NULL && gsr->indices != NULL);
}


fcs_int fcs_gridsort_resort_get_original_particles(fcs_gridsort_resort_t gsr)
{
  return gsr->noriginal_particles;
}


fcs_int fcs_gridsort_resort_get_sorted_particles(fcs_gridsort_resort_t gsr)
{
  if (gsr->nsorted_particles < 0) return gsr->noriginal_particles;

  return gsr->nsorted_particles;
}


static void copy_ints(fcs_int n, fcs_int *src, fcs_int *dst, fcs_int x, fcs_gridsort_index_t *indices)
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
    memcpy(dst, src, n * x * sizeof(fcs_int));
  }
}


static void pack_ints(fcs_int n, fcs_int *src, fcs_int x, fcs_gridsort_index_t *indices, char *dst)
{
  fcs_int i, j;
  size_t offset;
  

  offset = 0;
  for (i = 0; i < n; ++i)
  {
    *((fcs_gridsort_index_t *) (dst + offset)) = indices[i];
    offset += sizeof(fcs_gridsort_index_t);
    
    for (j = 0; j < x; ++j)
    {
      *((fcs_int *) (dst + offset)) = src[i * x + j];
      offset += sizeof(fcs_int);
    }
  }
}


static void unpack_ints(fcs_int n, char *src, fcs_int *dst, fcs_int x)
{
  fcs_int i, j, k;
  size_t offset;
  

  offset = 0;
  for (i = 0; i < n; ++i)
  {
    k = *((fcs_gridsort_index_t *) (src + offset));
    offset += sizeof(fcs_gridsort_index_t);
    
    for (j = 0; j < x; ++j)
    {
      dst[k * x + j] = *((fcs_int *) (src + offset));
      offset += sizeof(fcs_int);
    }
  }
}


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


static void pack_floats(fcs_int n, fcs_float *src, fcs_int x, fcs_gridsort_index_t *indices, char *dst)
{
  fcs_int i, j;
  size_t offset;
  

  offset = 0;
  for (i = 0; i < n; ++i)
  {
    *((fcs_gridsort_index_t *) (dst + offset)) = indices[i];
    offset += sizeof(fcs_gridsort_index_t);
    
    for (j = 0; j < x; ++j)
    {
      *((fcs_float *) (dst + offset)) = src[i * x + j];
      offset += sizeof(fcs_float);
    }
  }
}


static void unpack_floats(fcs_int n, char *src, fcs_float *dst, fcs_int x)
{
  fcs_int i, j, k;
  size_t offset;
  

  offset = 0;
  for (i = 0; i < n; ++i)
  {
    k = *((fcs_gridsort_index_t *) (src + offset));
    offset += sizeof(fcs_gridsort_index_t);
    
    for (j = 0; j < x; ++j)
    {
      dst[k * x + j] = *((fcs_float *) (src + offset));
      offset += sizeof(fcs_float);
    }
  }
}


static void copy_bytes(fcs_int n, char *src, char *dst, char x, fcs_gridsort_index_t *indices)
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


static void pack_bytes(fcs_int n, char *src, fcs_int x, fcs_gridsort_index_t *indices, char *dst)
{
  fcs_int i, j;
  size_t offset;
  

  offset = 0;
  for (i = 0; i < n; ++i)
  {
    *((fcs_gridsort_index_t *) (dst + offset)) = indices[i];
    offset += sizeof(fcs_gridsort_index_t);
    
    for (j = 0; j < x; ++j)
    {
      *((char *) (dst + offset)) = src[i * x + j];
      offset += sizeof(char);
    }
  }
}


static void unpack_bytes(fcs_int n, char *src, char *dst, fcs_int x)
{
  fcs_int i, j, k;
  size_t offset;
  

  offset = 0;
  for (i = 0; i < n; ++i)
  {
    k = *((fcs_gridsort_index_t *) (src + offset));
    offset += sizeof(fcs_gridsort_index_t);
    
    for (j = 0; j < x; ++j)
    {
      dst[k * x + j] = *((char *) (src + offset));
      offset += sizeof(char);
    }
  }
}


static int resort_tproc(void *b, ZMPI_Count x, void *tproc_data)
{
  size_t type_size = *((fcs_int *) tproc_data);
  
  fcs_gridsort_index_t idx = *((fcs_gridsort_index_t *) (((char *) b) + (x * type_size)));

  return GRIDSORT_INDEX_GET_PROC(idx);
}


static void resort_ints(fcs_gridsort_resort_t gsr, fcs_int *src, fcs_int *dst, fcs_int x, int size, int rank, MPI_Comm comm)
{
  void *send, *recv;
  size_t type_size = sizeof(fcs_gridsort_index_t) + x * sizeof(fcs_int);

  MPI_Datatype type;
  int lengths[2] = { 1, x };
  MPI_Aint displs[2] = { 0, sizeof(fcs_gridsort_index_t) };
  MPI_Datatype types[2] = { FCS_MPI_GRIDSORT_INDEX, FCS_MPI_INT };

  ZMPI_Tproc tproc;
  int received;


  send = malloc(gsr->noriginal_particles * type_size);
  recv = malloc(gsr->nsorted_particles * type_size);

  pack_ints(gsr->noriginal_particles, src, x, gsr->indices, send);

  MPI_Type_create_struct(2, lengths, displs, types, &type);
  MPI_Type_commit(&type);

  ZMPI_Tproc_create_tproc(&tproc, resort_tproc, ZMPI_TPROC_RESET_NULL, ZMPI_TPROC_EXDEF_NULL);

  ZMPI_Alltoall_specific(send, gsr->noriginal_particles, type, recv, gsr->nsorted_particles, type, tproc, &type_size, &received, comm);

  ZMPI_Tproc_free(&tproc);

  MPI_Type_free(&type);

  unpack_ints(gsr->nsorted_particles, recv, (dst)?dst:src, x);
  
  free(send);
  free(recv);
}


static void resort_floats(fcs_gridsort_resort_t gsr, fcs_float *src, fcs_float *dst, fcs_int x, int size, int rank, MPI_Comm comm)
{
  void *send, *recv;
  size_t type_size = sizeof(fcs_gridsort_index_t) + x * sizeof(fcs_float);

  MPI_Datatype type;
  int lengths[2] = { 1, x };
  MPI_Aint displs[2] = { 0, sizeof(fcs_gridsort_index_t) };
  MPI_Datatype types[2] = { FCS_MPI_GRIDSORT_INDEX, FCS_MPI_FLOAT };

  ZMPI_Tproc tproc;
  int received;


  send = malloc(gsr->noriginal_particles * type_size);
  recv = malloc(gsr->nsorted_particles * type_size);

  pack_floats(gsr->noriginal_particles, src, x, gsr->indices, send);

  MPI_Type_create_struct(2, lengths, displs, types, &type);
  MPI_Type_commit(&type);

  ZMPI_Tproc_create_tproc(&tproc, resort_tproc, ZMPI_TPROC_RESET_NULL, ZMPI_TPROC_EXDEF_NULL);

  ZMPI_Alltoall_specific(send, gsr->noriginal_particles, type, recv, gsr->nsorted_particles, type, tproc, &type_size, &received, comm);

  ZMPI_Tproc_free(&tproc);

  MPI_Type_free(&type);

  unpack_floats(gsr->nsorted_particles, recv, (dst)?dst:src, x);
  
  free(send);
  free(recv);
}


static void resort_bytes(fcs_gridsort_resort_t gsr, void *src, void *dst, fcs_int x, int size, int rank, MPI_Comm comm)
{
  void *send, *recv;
  size_t type_size = sizeof(fcs_gridsort_index_t) + x * sizeof(char);

  MPI_Datatype type;
  int lengths[2] = { 1, x };
  MPI_Aint displs[2] = { 0, sizeof(fcs_gridsort_index_t) };
  MPI_Datatype types[2] = { FCS_MPI_GRIDSORT_INDEX, MPI_BYTE };

  ZMPI_Tproc tproc;
  int received;


  send = malloc(gsr->noriginal_particles * type_size);
  recv = malloc(gsr->nsorted_particles * type_size);

  pack_bytes(gsr->noriginal_particles, src, x, gsr->indices, send);

  MPI_Type_create_struct(2, lengths, displs, types, &type);
  MPI_Type_commit(&type);

  ZMPI_Tproc_create_tproc(&tproc, resort_tproc, ZMPI_TPROC_RESET_NULL, ZMPI_TPROC_EXDEF_NULL);

  ZMPI_Alltoall_specific(send, gsr->noriginal_particles, type, recv, gsr->nsorted_particles, type, tproc, &type_size, &received, comm);

  ZMPI_Tproc_free(&tproc);

  MPI_Type_free(&type);

  unpack_bytes(gsr->nsorted_particles, recv, (dst)?dst:src, x);
  
  free(send);
  free(recv);
}


static void resort_1float(fcs_gridsort_resort_t gsr, fcs_float *src, fcs_float *dst, int size, int rank, MPI_Comm comm)
{
  back__p_elements_t sin, sout;

  back__p_tproc_t tproc;

#ifdef ALLTOALLV_PACKED
  fcs_int local_packed, global_packed, original_packed;
#endif


  back__p_SL_DEFCON(mpi.rank) = rank;

  back__p_mpi_datatypes_init();

  back__p_elem_set_size(&sin, gsr->noriginal_particles);
  back__p_elem_set_max_size(&sin, gsr->noriginal_particles);
  back__p_elem_set_keys(&sin, gsr->indices);
  back__p_elem_set_data(&sin, src);

  back__p_elem_set_size(&sout, 0);
  back__p_elem_set_max_size(&sout, 0);
  back__p_elem_set_keys(&sout, NULL);
  back__p_elem_set_data(&sout, NULL);

  back__p_tproc_create_tproc(&tproc, gridsort_back__p_tproc, back__p_TPROC_RESET_NULL, back__p_TPROC_EXDEF_NULL);

#ifdef ALLTOALLV_PACKED
  local_packed = ALLTOALLV_PACKED(size, sin.size);
  MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
  original_packed = back__p_SL_DEFCON(meas.packed); back__p_SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

  TIMING_SYNC(comm); TIMING_START(t[1]);
  back__p_mpi_elements_alltoall_specific(&sin, &sout, NULL, tproc, NULL, size, rank, comm);
  TIMING_SYNC(comm); TIMING_STOP(t[1]);

#ifdef ALLTOALLV_PACKED
  back__p_SL_DEFCON(meas.packed) = original_packed;
#endif

  back__p_tproc_free(&tproc);

  if (gsr->nsorted_particles != sout.size)
    fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" back__p_slint_fmt "!\n", rank, gsr->nsorted_particles, sout.size);

  copy_floats(sout.size, sout.data1, (dst)?dst:src, 1, sout.keys);

  back__p_elements_free(&sout);

  back__p_mpi_datatypes_release();
}


static void resort_3floats(fcs_gridsort_resort_t gsr, fcs_float *src, fcs_float *dst, int size, int rank, MPI_Comm comm)
{
  back_f__elements_t sin, sout;

  back_f__tproc_t tproc;

#ifdef ALLTOALLV_PACKED
  fcs_int local_packed, global_packed, original_packed;
#endif


  back_f__SL_DEFCON(mpi.rank) = rank;

  back_f__mpi_datatypes_init();

  back_f__elem_set_size(&sin, gsr->noriginal_particles);
  back_f__elem_set_max_size(&sin, gsr->noriginal_particles);
  back_f__elem_set_keys(&sin, gsr->indices);
  back_f__elem_set_data(&sin, src);

  back_f__elem_set_size(&sout, 0);
  back_f__elem_set_max_size(&sout, 0);
  back_f__elem_set_keys(&sout, NULL);
  back_f__elem_set_data(&sout, NULL);

  back_f__tproc_create_tproc(&tproc, gridsort_back_f__tproc, back_f__TPROC_RESET_NULL, back_f__TPROC_EXDEF_NULL);

#ifdef ALLTOALLV_PACKED
  local_packed = ALLTOALLV_PACKED(size, sin.size);
  MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
  original_packed = back_f__SL_DEFCON(meas.packed); back_f__SL_DEFCON(meas.packed) = (global_packed > 0);
#endif

  TIMING_SYNC(comm); TIMING_START(t[1]);
  back_f__mpi_elements_alltoall_specific(&sin, &sout, NULL, tproc, NULL, size, rank, comm);
  TIMING_SYNC(comm); TIMING_STOP(t[1]);

#ifdef ALLTOALLV_PACKED
  back_f__SL_DEFCON(meas.packed) = original_packed;
#endif

  back_f__tproc_free(&tproc);

  if (gsr->nsorted_particles != sout.size)
    fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" back_f__slint_fmt "!\n", rank, gsr->nsorted_particles, sout.size);

  copy_floats(sout.size, sout.data0, (dst)?dst:src, 3, sout.keys);

  back_f__elements_free(&sout);

  back_f__mpi_datatypes_release();
}


void fcs_gridsort_resort_ints(fcs_gridsort_resort_t gsr, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  int comm_size, comm_rank;


  if (gsr->indices == NULL)
  {
    if (dst) copy_ints(gsr->noriginal_particles, dst, src, n, NULL);
    
    return;
  }

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  
  resort_ints(gsr, src, dst, n, comm_size, comm_rank, comm);
}


void fcs_gridsort_resort_floats(fcs_gridsort_resort_t gsr, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  int comm_size, comm_rank;


  if (gsr->indices == NULL)
  {
    if (dst) copy_floats(gsr->noriginal_particles, dst, src, n, NULL);
    
    return;
  }

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  
  switch (n)
  {
    case 1:
      resort_1float(gsr, src, dst, comm_size, comm_rank, comm);
      break;
    case 3:
      resort_3floats(gsr, src, dst, comm_size, comm_rank, comm);
      break;
    default:
      resort_floats(gsr, src, dst, n, comm_size, comm_rank, comm);
      break;
  }
}


void fcs_gridsort_resort_bytes(fcs_gridsort_resort_t gsr, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  int comm_size, comm_rank;


  if (gsr->indices == NULL)
  {
    if (dst) copy_bytes(gsr->noriginal_particles, dst, src, n, NULL);
    
    return;
  }

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  
  resort_bytes(gsr, src, dst, n, comm_size, comm_rank, comm);
}
