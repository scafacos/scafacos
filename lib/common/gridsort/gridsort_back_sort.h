/*
  Copyright (C) 2011, 2012, 2013, 2014, 2015 Michael Hofmann
  
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

/* GRIDSORT_BACK_SORT_NAME */


#define CONCAT_(_a_, _b_)  _a_##_b_
#define CONCAT(_a_, _b_)   CONCAT_(_a_, _b_)

#define FCS_BACK(_x_)    CONCAT(GRIDSORT_BACK_SL_PREFIX, _x_)


static void CONCAT(GRIDSORT_PREFIX, CONCAT(GRIDSORT_BACK_PREFIX, GRIDSORT_BACK_SORT_NAME))(
  fcs_int nsorted_particles, fcs_gridsort_index_t *sorted_indices, fcs_float *sorted_field, fcs_float *sorted_potentials,
  fcs_int noriginal_particles, fcs_gridsort_index_t *original_indices, fcs_float *original_field, fcs_float *original_potentials,
  fcs_int nprocs, int *procs, double *t, int comm_size, int comm_rank, MPI_Comm comm)
{
  FCS_BACK(elements_t) sin, sout;

#if !SORT_BACKWARD_LOCAL_INPLACE
  FCS_BACK(elements_t) sx;
#endif

  FCS_BACK(tproc_t) tproc;

  FCS_BACK(permute_generic_t) tloc = FCS_BACK(PERMUTE_GENERIC_INIT_TLOC)(CONCAT(GRIDSORT_PREFIX, CONCAT(GRIDSORT_BACK_PREFIX, GRIDSORT_BACK_TLOC_NAME)));

#ifdef ALLTOALLV_PACKED
  fcs_int local_packed, global_packed, original_packed;
#endif


  FCS_BACK(SL_DEFCON)(mpi.rank) = comm_rank;

  FCS_BACK(mpi_datatypes_init)();

  FCS_BACK(elem_set_size)(&sin, nsorted_particles);
  FCS_BACK(elem_set_max_size)(&sin, nsorted_particles);
  FCS_BACK(elem_set_keys)(&sin, sorted_indices);
  FCS_BACK(elem_set_data)(&sin, sorted_field, sorted_potentials);

#if SORT_BACKWARD_LOCAL_INPLACE
  original_indices = malloc(noriginal_particles * sizeof(fcs_gridsort_index_t));

  FCS_BACK(elem_set_size)(&sout, noriginal_particles);
  FCS_BACK(elem_set_max_size)(&sout, noriginal_particles);
  FCS_BACK(elem_set_keys)(&sout, original_indices);
  FCS_BACK(elem_set_data)(&sout, original_field, original_potentials);
#else
  FCS_BACK(elem_null)(&sout);
#endif

  FCS_BACK(tproc_create_tproc)(&tproc, CONCAT(GRIDSORT_PREFIX, CONCAT(GRIDSORT_BACK_PREFIX, GRIDSORT_BACK_TPROC_NAME)), FCS_BACK(TPROC_RESET_NULL), FCS_BACK(TPROC_EXDEF_NULL));

#ifdef GRIDSORT_BACK_PROCLIST
  if (procs) FCS_BACK(tproc_set_proclists)(tproc, nprocs, procs, nprocs, procs, comm_size, comm_rank, comm);
#endif

#ifdef ALLTOALLV_PACKED
  local_packed = ALLTOALLV_PACKED(comm_size, FCS_BACK(elem_get_size(&sin)));
  MPI_Allreduce(&local_packed, &global_packed, 1, FCS_MPI_INT, MPI_SUM, comm);
  original_packed = FCS_BACK(SL_DEFCON)(meas.packed); FCS_BACK(SL_DEFCON)(meas.packed) = (global_packed > 0);
#endif

  TIMING_SYNC(comm); TIMING_START(*t);
  FCS_BACK(mpi_elements_alltoall_specific)(&sin, &sout, NULL, tproc, NULL, comm_size, comm_rank, comm);
  TIMING_SYNC(comm); TIMING_STOP(*t);

#ifdef ALLTOALLV_PACKED
  FCS_BACK(SL_DEFCON)(meas.packed) = original_packed;
#endif

  FCS_BACK(tproc_free)(&tproc);

  if (noriginal_particles != FCS_BACK(elem_get_size(&sout)))
    fprintf(stderr, "%d: error: wanted %" FCS_LMOD_INT "d particles, but got %" FCS_BACK(slint_fmt) "!\n", comm_rank, noriginal_particles, FCS_BACK(elem_get_size(&sout)));

#if SORT_BACKWARD_LOCAL_INPLACE
  FCS_BACK(permute_generic_ip)(&sout, NULL, &tloc, NULL);

  free(indices);
#else
  FCS_BACK(elem_set_size)(&sx, FCS_BACK(elem_get_size)(&sout));
  FCS_BACK(elem_set_max_size)(&sx, FCS_BACK(elem_get_max_size)(&sout));
  FCS_BACK(elem_set_keys)(&sx, FCS_BACK(elem_get_keys)(&sout));
  FCS_BACK(elem_set_data)(&sx, original_field, original_potentials);

  FCS_BACK(permute_generic_db)(&sout, &sx, &tloc, NULL);

  FCS_BACK(elements_free)(&sout);
#endif

  FCS_BACK(mpi_datatypes_release)();
}


#undef FCS_BACK

#undef CONCAT_
#undef CONCAT

#undef GRIDSORT_BACK_SORT_NAME
