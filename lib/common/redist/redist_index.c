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

#include "../gridsort/include/sl_back_fp.h"
#include "../gridsort/include/sl_back_f_.h"
#include "../gridsort/include/sl_back__p.h"

#include "z_tools.h"
#include "redist_common.h"
#include "redist_index.h"


fcs_redist_index_t *fcs_redist_indices_alloc(fcs_int nindices)
{
  return malloc(nindices * sizeof(fcs_redist_index_t));
}


void fcs_redist_indices_free(fcs_redist_index_t *indices)
{
  if (indices) free(indices);
}


void fcs_redist_indices_print(fcs_int nindices, fcs_redist_index_t *indices)
{
  fcs_int i;


  for (i = 0; i < nindices; ++i)
  {
    printf("fcs_redist_index: %" FCS_LMOD_INT "d: " FCS_REDIST_INDEX_STR "\n", i, FCS_REDIST_INDEX_PARAM(indices[i]));
  }
}


void fcs_redist_indices_init(fcs_int nindices, fcs_redist_index_t *indices, int rank)
{
  fcs_int i;
  fcs_redist_index_t index_rank;


  index_rank = FCS_REDIST_INDEX_VAL_PROC(rank);

  for (i = 0; i < nindices; ++i) indices[i] = index_rank + i;
}


#define GRIDSORT_INDEX_IS_VALID(_x_)  FCS_REDIST_INDEX_IS_VALID(_x_)
#define GRIDSORT_INDEX_GET_PROC(_x_)  FCS_REDIST_INDEX_GET_PROC(_x_)
#define GRIDSORT_INDEX_GET_POS(_x_)   FCS_REDIST_INDEX_GET_POS(_x_)

#define fcs_gridsort_index_t  fcs_redist_index_t
#include "../gridsort/gridsort_back_tproc_sl.h"


void fcs_redist_indices_sort_back_results(fcs_int nindices, fcs_redist_index_t *indices, fcs_float *field, fcs_float *potentials, fcs_int nback_results, fcs_float *back_field, fcs_float *back_potentials, MPI_Comm comm)
{
  int comm_size, comm_rank;

  fcs_int i, j, type;

  fcs_back_fp_elements_t sin0, sout0;
  fcs_back_f__elements_t sin1, sout1;
  fcs_back__p_elements_t sin2, sout2;

  fcs_back_fp_tproc_t tproc0;
  fcs_back_f__tproc_t tproc1;
  fcs_back__p_tproc_t tproc2;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if (back_field && back_potentials) type = 0;
  else if (back_field && !back_potentials) type = 1;
  else if (!back_field && back_potentials) type = 2;
  else type = 3;

  switch (type)
  {
    case 0:
      fcs_back_fp_SL_DEFCON(mpi.rank) = comm_rank;

      fcs_back_fp_mpi_datatypes_init();

      fcs_back_fp_elem_set_size(&sin0, nindices);
      fcs_back_fp_elem_set_max_size(&sin0, nindices);
      fcs_back_fp_elem_set_keys(&sin0, indices);
      fcs_back_fp_elem_set_data(&sin0, field, potentials);

      fcs_back_fp_elem_null(&sout0);

      fcs_back_fp_tproc_create_tproc(&tproc0, gridsort_back_fp_tproc, fcs_back_fp_TPROC_RESET_NULL, fcs_back_fp_TPROC_EXDEF_NULL);

      fcs_back_fp_mpi_elements_alltoall_specific(&sin0, &sout0, NULL, tproc0, NULL, comm_size, comm_rank, comm);

      fcs_back_fp_tproc_free(&tproc0);

      if (nback_results != sout0.size)
        fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got %" fcs_back_fp_slint_fmt "!\n", comm_rank, nback_results, sout0.size);

      for (i = 0; i < sout0.size; ++i)
      {
        j = FCS_REDIST_INDEX_GET_POS(sout0.keys[i]);

        back_field[3 * j + 0] = sout0.data0[3 * i + 0];
        back_field[3 * j + 1] = sout0.data0[3 * i + 1];
        back_field[3 * j + 2] = sout0.data0[3 * i + 2];

        back_potentials[j] = sout0.data1[i];
      }

      fcs_back_fp_elements_free(&sout0);

      fcs_back_fp_mpi_datatypes_release();

      break;

    case 1:
      fcs_back_f__SL_DEFCON(mpi.rank) = comm_rank;

      fcs_back_f__mpi_datatypes_init();

      fcs_back_f__elem_set_size(&sin1, nindices);
      fcs_back_f__elem_set_max_size(&sin1, nindices);
      fcs_back_f__elem_set_keys(&sin1, indices);
      fcs_back_f__elem_set_data(&sin1, field);

      fcs_back_f__elem_null(&sout1);

      fcs_back_f__tproc_create_tproc(&tproc1, gridsort_back_f__tproc, fcs_back_f__TPROC_RESET_NULL, fcs_back_f__TPROC_EXDEF_NULL);

      fcs_back_f__mpi_elements_alltoall_specific(&sin1, &sout1, NULL, tproc1, NULL, comm_size, comm_rank, comm);

      fcs_back_f__tproc_free(&tproc1);

      if (nback_results != sout1.size)
        fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got %" fcs_back_f__slint_fmt "!\n", comm_rank, nback_results, sout1.size);

      for (i = 0; i < sout1.size; ++i)
      {
        j = FCS_REDIST_INDEX_GET_POS(sout0.keys[i]);

        back_field[3 * j + 0] = sout1.data0[3 * i + 0];
        back_field[3 * j + 1] = sout1.data0[3 * i + 1];
        back_field[3 * j + 2] = sout1.data0[3 * i + 2];
      }

      fcs_back_f__elements_free(&sout1);

      fcs_back_f__mpi_datatypes_release();

      break;

    case 2:
      fcs_back__p_SL_DEFCON(mpi.rank) = comm_rank;

      fcs_back__p_mpi_datatypes_init();

      fcs_back__p_elem_set_size(&sin2, nindices);
      fcs_back__p_elem_set_max_size(&sin2, nindices);
      fcs_back__p_elem_set_keys(&sin2, indices);
      fcs_back__p_elem_set_data(&sin2, potentials);

      fcs_back__p_elem_null(&sout2);

      fcs_back__p_tproc_create_tproc(&tproc2, gridsort_back__p_tproc, fcs_back__p_TPROC_RESET_NULL, fcs_back__p_TPROC_EXDEF_NULL);

      fcs_back__p_mpi_elements_alltoall_specific(&sin2, &sout2, NULL, tproc2, NULL, comm_size, comm_rank, comm);

      fcs_back__p_tproc_free(&tproc2);

      if (nback_results != sout2.size)
        fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got only %" fcs_back__p_slint_fmt "!\n", comm_rank, nback_results, sout2.size);

      for (i = 0; i < sout2.size; ++i)
      {
        j = FCS_REDIST_INDEX_GET_POS(sout0.keys[i]);

        back_potentials[j] = sout2.data1[i];
      }

      fcs_back__p_elements_free(&sout2);

      fcs_back__p_mpi_datatypes_release();

      break;
  }
}
