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

#ifdef HAVE_ZMPI_ATASP_H
# include "zmpi_atasp.h"
#endif

#include "resort.h"


#if defined(FCS_ENABLE_DEBUG_RESORT)
# define DO_DEBUG
# define DEBUG_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define DEBUG_CMD(_cmd_)  Z_NOP()
#endif
#define DEBUG_PRINT_PREFIX  "RESORT_DEBUG: "

#if defined(FCS_ENABLE_INFO_RESORT)
# define DO_INFO
# define INFO_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define INFO_CMD(_cmd_)  Z_NOP()
#endif
#define INFO_PRINT_PREFIX  "RESORT_INFO: "

#if defined(FCS_ENABLE_TIMING_RESORT)
# define DO_TIMING
# define TIMING_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define TIMING_CMD(_cmd_)  Z_NOP()
#endif
#define TIMING_PRINT_PREFIX  "RESORT_TIMING: "


#define RESORT_PROCLIST


fcs_resort_index_t *fcs_resort_indices_alloc(fcs_int nindices)
{
  return malloc(nindices * sizeof(fcs_resort_index_t));
}


void fcs_resort_indices_free(fcs_resort_index_t *indices)
{
  if (indices) free(indices);
}


void fcs_resort_indices_init(fcs_int nindices, fcs_resort_index_t *indices, int rank)
{
  fcs_int i;
  fcs_resort_index_t index_rank;


  index_rank = FCS_RESORT_INDEX_VAL_PROC(rank);

  for (i = 0; i < nindices; ++i) indices[i] = index_rank + i;
}


void fcs_resort_indices_print(fcs_int nindices, fcs_resort_index_t *indices)
{
  fcs_int i;


  for (i = 0; i < nindices; ++i)
  {
    printf(" %" FCS_LMOD_INT "d: " FCS_RESORT_INDEX_STR "\n", i, FCS_RESORT_INDEX_PARAM(indices[i]));
  }
}


void fcs_resort_create(fcs_resort_t *resort)
{
  *resort = malloc(sizeof(**resort));

  (*resort)->noriginal_particles = -1;
  (*resort)->nsorted_particles = -1;
  (*resort)->indices = NULL;

  (*resort)->nprocs = -1;
  (*resort)->procs = NULL;
}


void fcs_resort_destroy(fcs_resort_t *resort)
{
  if (*resort == FCS_RESORT_NULL) return;

  fcs_resort_free_indices(*resort);

  fcs_resort_set_proclists(*resort, -1, NULL);

  free(*resort);

  *resort = FCS_RESORT_NULL;
}


void fcs_resort_print(fcs_resort_t resort)
{
  fcs_int i;


  printf("fcs_resort: %p -> [%" FCS_LMOD_INT "d, %" FCS_LMOD_INT "d, %p]\n", resort, resort->noriginal_particles, resort->nsorted_particles, resort->indices);

  if (resort->indices)
    for (i = 0; i < resort->noriginal_particles; ++i)
      printf("%" FCS_LMOD_INT "d: " FCS_RESORT_INDEX_STR "\n", i, FCS_RESORT_INDEX_PARAM(resort->indices[i]));
  else
    printf("-> ID\n");
}


fcs_int fcs_resort_is_available(fcs_resort_t resort)
{
  return (resort != FCS_RESORT_NULL && resort->indices != NULL);
}


void fcs_resort_set_original_particles(fcs_resort_t resort, fcs_int original_particles)
{
  resort->noriginal_particles = original_particles;
}


fcs_int fcs_resort_get_original_particles(fcs_resort_t resort)
{
  return resort->noriginal_particles;
}


void fcs_resort_set_sorted_particles(fcs_resort_t resort, fcs_int sorted_particles)
{
  resort->nsorted_particles = sorted_particles;
}


fcs_int fcs_resort_get_sorted_particles(fcs_resort_t resort)
{
  if (resort->nsorted_particles < 0) return resort->noriginal_particles;

  return resort->nsorted_particles;
}


void fcs_resort_alloc_indices(fcs_resort_t resort)
{
  if (resort->noriginal_particles < 0) return;

  resort->indices = fcs_resort_indices_alloc(resort->noriginal_particles);
}


void fcs_resort_free_indices(fcs_resort_t resort)
{
  fcs_resort_indices_free(resort->indices);

  resort->indices = NULL;
}


fcs_resort_index_t *fcs_resort_get_indices(fcs_resort_t resort)
{
  return resort->indices;
}


void fcs_resort_set_proclists(fcs_resort_t resort, fcs_int nprocs, int *procs)
{
  fcs_int i;


  resort->nprocs = nprocs;

  if (resort->procs) free(resort->procs);
  resort->procs = NULL;

  if (nprocs < 0) return;

  resort->procs = malloc(nprocs * sizeof(int));

  for (i = 0; i < nprocs; ++i) resort->procs[i] = procs[i];
}


static void copy_ints(fcs_int n, fcs_int *src, fcs_int *dst, fcs_int x, fcs_resort_index_t *indices)
{
  fcs_int i, j, k;


  if (indices)
  {
    for (i = 0; i < n; ++i)
    {
      k = FCS_RESORT_INDEX_GET_POS(indices[i]);

      for (j = 0; j < x; ++j) dst[k * x + j] = src[i * x + j];
    }

  } else
  {
    memcpy(dst, src, n * x * sizeof(fcs_int));
  }
}


static void pack_ints(fcs_int n, fcs_int *src, fcs_int x, fcs_resort_index_t *indices, char *dst)
{
  fcs_int i, j;
  size_t offset;
  

  offset = 0;
  for (i = 0; i < n; ++i)
  {
    *((fcs_resort_index_t *) (dst + offset)) = indices[i];
    offset += sizeof(fcs_resort_index_t);
    
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
    k = FCS_RESORT_INDEX_GET_POS(*((fcs_resort_index_t *) (src + offset)));
    offset += sizeof(fcs_resort_index_t);
    
    for (j = 0; j < x; ++j)
    {
      dst[k * x + j] = *((fcs_int *) (src + offset));
      offset += sizeof(fcs_int);
    }
  }
}


static void copy_floats(fcs_int n, fcs_float *src, fcs_float *dst, fcs_int x, fcs_resort_index_t *indices)
{
  fcs_int i, j, k;


  if (indices)
  {
    for (i = 0; i < n; ++i)
    {
      k = FCS_RESORT_INDEX_GET_POS(indices[i]);

      for (j = 0; j < x; ++j) dst[k * x + j] = src[i * x + j];
    }

  } else
  {
    memcpy(dst, src, n * x * sizeof(fcs_float));
  }
}


static void pack_floats(fcs_int n, fcs_float *src, fcs_int x, fcs_resort_index_t *indices, char *dst)
{
  fcs_int i, j;
  size_t offset;
  

  offset = 0;
  for (i = 0; i < n; ++i)
  {
    *((fcs_resort_index_t *) (dst + offset)) = indices[i];
    offset += sizeof(fcs_resort_index_t);
    
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
    k = FCS_RESORT_INDEX_GET_POS(*((fcs_resort_index_t *) (src + offset)));
    offset += sizeof(fcs_resort_index_t);
    
    for (j = 0; j < x; ++j)
    {
      dst[k * x + j] = *((fcs_float *) (src + offset));
      offset += sizeof(fcs_float);
    }
  }
}


static void copy_bytes(fcs_int n, char *src, char *dst, char x, fcs_resort_index_t *indices)
{
  fcs_int i, j, k;


  if (indices)
  {
    for (i = 0; i < n; ++i)
    {
      k = FCS_RESORT_INDEX_GET_POS(indices[i]);

      for (j = 0; j < x; ++j) dst[k * x + j] = src[i * x + j];
    }

  } else
  {
    memcpy(dst, src, n * x * sizeof(fcs_float));
  }
}


static void pack_bytes(fcs_int n, char *src, fcs_int x, fcs_resort_index_t *indices, char *dst)
{
  fcs_int i, j;
  size_t offset;
  

  offset = 0;
  for (i = 0; i < n; ++i)
  {
    *((fcs_resort_index_t *) (dst + offset)) = indices[i];
    offset += sizeof(fcs_resort_index_t);
    
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
    k = FCS_RESORT_INDEX_GET_POS(*((fcs_resort_index_t *) (src + offset)));
    offset += sizeof(fcs_resort_index_t);
    
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
  
  fcs_resort_index_t idx = *((fcs_resort_index_t *) (((char *) b) + (x * type_size)));

  return FCS_RESORT_INDEX_GET_PROC(idx);
}


static void resort_ints(fcs_resort_t resort, fcs_int *src, fcs_int *dst, fcs_int x, MPI_Comm comm)
{
  void *send, *recv;
  size_t type_size = sizeof(fcs_resort_index_t) + x * sizeof(fcs_int);

  MPI_Datatype type;
  int lengths[2] = { 1, x };
  MPI_Aint displs[2] = { 0, sizeof(fcs_resort_index_t) };
  MPI_Datatype types[2] = { FCS_MPI_RESORT_INDEX, FCS_MPI_INT };

  ZMPI_Tproc tproc;
  int received;


  send = malloc(resort->noriginal_particles * type_size);
  recv = malloc(resort->nsorted_particles * type_size);

  pack_ints(resort->noriginal_particles, src, x, resort->indices, send);

  MPI_Type_create_struct(2, lengths, displs, types, &type);
  MPI_Type_commit(&type);

  ZMPI_Tproc_create_tproc(&tproc, resort_tproc, ZMPI_TPROC_RESET_NULL, ZMPI_TPROC_EXDEF_NULL);

#ifdef RESORT_PROCLIST
  if (resort->nprocs >= 0) ZMPI_Tproc_set_proclists(&tproc, resort->nprocs, resort->procs, resort->nprocs, resort->procs, comm);
#endif

  ZMPI_Alltoall_specific(send, resort->noriginal_particles, type, recv, resort->nsorted_particles, type, tproc, &type_size, &received, comm);

  ZMPI_Tproc_free(&tproc);

  MPI_Type_free(&type);

  unpack_ints(resort->nsorted_particles, recv, (dst)?dst:src, x);
  
  free(send);
  free(recv);
}


static void resort_floats(fcs_resort_t resort, fcs_float *src, fcs_float *dst, fcs_int x, MPI_Comm comm)
{
  void *send, *recv;
  size_t type_size = sizeof(fcs_resort_index_t) + x * sizeof(fcs_float);

  MPI_Datatype type;
  int lengths[2] = { 1, x };
  MPI_Aint displs[2] = { 0, sizeof(fcs_resort_index_t) };
  MPI_Datatype types[2] = { FCS_MPI_RESORT_INDEX, FCS_MPI_FLOAT };

  ZMPI_Tproc tproc;
  int received;


  send = malloc(resort->noriginal_particles * type_size);
  recv = malloc(resort->nsorted_particles * type_size);

  pack_floats(resort->noriginal_particles, src, x, resort->indices, send);

  MPI_Type_create_struct(2, lengths, displs, types, &type);
  MPI_Type_commit(&type);

  ZMPI_Tproc_create_tproc(&tproc, resort_tproc, ZMPI_TPROC_RESET_NULL, ZMPI_TPROC_EXDEF_NULL);

#ifdef RESORT_PROCLIST
  if (resort->nprocs >= 0) ZMPI_Tproc_set_proclists(&tproc, resort->nprocs, resort->procs, resort->nprocs, resort->procs, comm);
#endif

  ZMPI_Alltoall_specific(send, resort->noriginal_particles, type, recv, resort->nsorted_particles, type, tproc, &type_size, &received, comm);

  ZMPI_Tproc_free(&tproc);

  MPI_Type_free(&type);

  unpack_floats(resort->nsorted_particles, recv, (dst)?dst:src, x);
  
  free(send);
  free(recv);
}


static void resort_bytes(fcs_resort_t resort, void *src, void *dst, fcs_int x, MPI_Comm comm)
{
  void *send, *recv;
  size_t type_size = sizeof(fcs_resort_index_t) + x * sizeof(char);

  MPI_Datatype type;
  int lengths[2] = { 1, x };
  MPI_Aint displs[2] = { 0, sizeof(fcs_resort_index_t) };
  MPI_Datatype types[2] = { FCS_MPI_RESORT_INDEX, MPI_BYTE };

  ZMPI_Tproc tproc;
  int received;


  send = malloc(resort->noriginal_particles * type_size);
  recv = malloc(resort->nsorted_particles * type_size);

  pack_bytes(resort->noriginal_particles, src, x, resort->indices, send);

  MPI_Type_create_struct(2, lengths, displs, types, &type);
  MPI_Type_commit(&type);

  ZMPI_Tproc_create_tproc(&tproc, resort_tproc, ZMPI_TPROC_RESET_NULL, ZMPI_TPROC_EXDEF_NULL);

#ifdef RESORT_PROCLIST
  if (resort->nprocs >= 0) ZMPI_Tproc_set_proclists(&tproc, resort->nprocs, resort->procs, resort->nprocs, resort->procs, comm);
#endif

  ZMPI_Alltoall_specific(send, resort->noriginal_particles, type, recv, resort->nsorted_particles, type, tproc, &type_size, &received, comm);

  ZMPI_Tproc_free(&tproc);

  MPI_Type_free(&type);

  unpack_bytes(resort->nsorted_particles, recv, (dst)?dst:src, x);
  
  free(send);
  free(recv);
}


void fcs_resort_resort_ints(fcs_resort_t resort, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  if (resort->indices == NULL)
  {
    if (dst) copy_ints(resort->noriginal_particles, src, dst, n, NULL);
    
    return;
  }
  
  resort_ints(resort, src, dst, n, comm);
}


void fcs_resort_resort_floats(fcs_resort_t resort, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  if (resort->indices == NULL)
  {
    if (dst) copy_floats(resort->noriginal_particles, src, dst, n, NULL);
    
    return;
  }

  resort_floats(resort, src, dst, n, comm);
}


void fcs_resort_resort_bytes(fcs_resort_t resort, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  if (resort->indices == NULL)
  {
    if (dst) copy_bytes(resort->noriginal_particles, src, dst, n, NULL);
    
    return;
  }
  
  resort_bytes(resort, src, dst, n, comm);
}
