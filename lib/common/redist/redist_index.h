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


#ifndef __REDIST_INDEX_H__
#define __REDIST_INDEX_H__


#ifdef __cplusplus
extern "C" {
#endif


#include <mpi.h>


typedef long long fcs_redist_index_t;
#define FCS_MPI_REDIST_INDEX  MPI_LONG_LONG


#define FCS_REDIST_INDEX_IS_VALID(_x_)       ((_x_) >= 0)
#define FCS_REDIST_INDEX_VAL_PROC(_proc_)    (((fcs_redist_index_t) (_proc_)) << 32)
#define FCS_REDIST_INDEX_VAL_POS(_pos_)      ((fcs_redist_index_t) (_pos_))
#define FCS_REDIST_INDEX_VAL(_proc_, _pos_)  (FCS_REDIST_INDEX_VAL_PROC(_proc_) + FCS_REDIST_INDEX_VAL_POS(_pos_))
#define FCS_REDIST_INDEX_GET_PROC(_x_)       ((_x_) >> 32)
#define FCS_REDIST_INDEX_GET_POS(_x_)        ((_x_) & 0x00000000FFFFFFFFLL)

#define FCS_REDIST_INDEX_STR         "(%lld,%lld)"
#define FCS_REDIST_INDEX_PARAM(_x_)  (FCS_REDIST_INDEX_IS_VALID(_x_)?FCS_REDIST_INDEX_GET_PROC(_x_):-1), (FCS_REDIST_INDEX_IS_VALID(_x_)?FCS_REDIST_INDEX_GET_POS(_x_):-1)


fcs_redist_index_t *fcs_redist_indices_alloc(fcs_int nindices);
void fcs_redist_indices_free(fcs_redist_index_t *indices);
void fcs_redist_indices_print(fcs_int nindices, fcs_redist_index_t *indices);
void fcs_redist_indices_init(fcs_int nindices, fcs_redist_index_t *indices, int rank);
void fcs_redist_indices_sort_back_results(fcs_int nindices, fcs_redist_index_t *indices, fcs_float *field, fcs_float *potentials, fcs_int nback_results, fcs_float *back_field, fcs_float *back_potentials, MPI_Comm comm);


#ifdef __cplusplus
}
#endif


#endif /* __REDIST_INDEX_H__ */
