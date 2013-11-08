/*
 *  Copyright (C) 2011, 2012, 2013 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS/FMM.
 *  
 *  ScaFaCoS/FMM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS/FMM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "sl_back_idx.h"

#ifdef HAVE_ZMPI_ATASP_H
# include "zmpi_atasp.h"
#endif

#include "config_fmm_sort.h"
#include "rename_fmm_sort.h"

#include "common.h"

#include "mpi_fmm_resort.h"


void fmm_resort_create(fcs_fmm_resort_t *fmm_resort, fcs_int nparticles, MPI_Comm comm)
{
  *fmm_resort = malloc(sizeof(**fmm_resort));

  fcs_resort_create(&(*fmm_resort)->resort);

  fcs_resort_set_original_particles((*fmm_resort)->resort, nparticles);

  (*fmm_resort)->comm = comm;
}


void fmm_resort_destroy(fcs_fmm_resort_t *fmm_resort)
{
  if (*fmm_resort == FCS_FMM_RESORT_NULL) return;

  fcs_resort_destroy(&(*fmm_resort)->resort);
  
  free(*fmm_resort);

  *fmm_resort = FCS_RESORT_NULL;
}


static void mpi_fmm_sort_idx(pint_t ntotal, pint_t nin, pint_t nout, back_idx_slkey_t *addr, fcs_resort_index_t *indices, double *t, int size, int rank, MPI_Comm comm)
{
  pint_t i, base_addr;
  back_idx_slkey_t lh_addrs[2];

  back_idx_elements_t sin, sout;

  fcs_resort_index_t *indices2, index_rank;


  TIMING_SYNC(comm); TIMING_START(t[0]);

  indices2 = malloc(nin * sizeof(fcs_resort_index_t));

  index_rank = FCS_RESORT_INDEX_VAL_PROC(rank);

  for (i = 0; i < nin; ++i) indices2[i] = index_rank + i;

/*  printf("nin = %" PARAM_INTEGER_FMT "\n", nin);
  for (i = 0; i < nin; ++i)
  {
    printf(" %" PARAM_INTEGER_FMT ": %" INTEGER_FMT "  " FCS_RESORT_INDEX_STR "\n", i, addr[i], FCS_RESORT_INDEX_PARAM(indices2[i]));
  }*/
  
  back_idx_elem_set_size(&sin, nin);
  back_idx_elem_set_max_size(&sin, nin);
  back_idx_elem_set_keys(&sin, addr);
  back_idx_elem_set_data(&sin, indices2);

  back_idx_elements_alloc(&sout, nout, SLCM_KEYS|SLCM_DATA);

  back_idx_mpi_datatypes_init();

  base_addr = 0;
  MPI_Exscan(&nout, &base_addr, 1, PARAM_INTEGER_MPI, MPI_SUM, comm);

  lh_addrs[0] = base_addr;
  lh_addrs[1] = base_addr + nout - 1;

  TIMING_SYNC(comm); TIMING_START(t[1]);
  back_idx_mpi_sort_back(&sin, &sout, NULL, lh_addrs, ntotal, size, rank, comm);
  TIMING_SYNC(comm); TIMING_STOP(t[1]);

  if (nout != back_idx_elem_get_size(&sout))
    fprintf(stderr, "%d: error: wanted %" PARAM_INTEGER_FMT " particles, but got only %" back_idx_slint_fmt "!\n", rank, nout, back_idx_elem_get_size(&sout));

  back_idx_mpi_datatypes_release();

  for (i = 0; i < nout; ++i) indices[i] = sout.data0[i];

  back_idx_elements_free(&sout);
  
  free(indices2);

  TIMING_SYNC(comm); TIMING_STOP(t[0]);
}


void mpi_fmm_resort_init(fcs_fmm_resort_t *fmm_resort, pint_t *ntotal, pint_t *nlocal, back_idx_slkey_t *addr)
{
  MPI_Comm comm;
  int size, rank;
#ifdef DO_MPI_INIT
  int flag;
#endif

#ifdef DO_TIMING
  double t[4];
#endif


#ifdef DO_MPI_INIT
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(NULL, NULL);
#endif

  comm = (*fmm_resort)->comm;

  if (comm == MPI_COMM_NULL) comm = MPI_COMM_WORLD;

  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size((*fmm_resort)->comm, &size);
  MPI_Comm_rank((*fmm_resort)->comm, &rank);

#define I_AM_MASTER  (rank == 0)

/*  printf("fmm_resort_init: %p (%" FCS_LMOD_INT "d,%" FCS_LMOD_INT "d,%p), %" PARAM_INTEGER_FMT ", %" PARAM_INTEGER_FMT ", %p\n",
    *fmm_resort, (*fmm_resort)->noriginal_particles, (*fmm_resort)->nsorted_particles, (*fmm_resort)->indices, *ntotal, *nlocal, addr);*/

  fcs_resort_set_sorted_particles((*fmm_resort)->resort, *nlocal);

  fcs_resort_alloc_indices((*fmm_resort)->resort);

  TIMING_SYNC(comm); TIMING_START(t[1]);

  mpi_fmm_sort_idx(*ntotal, *nlocal, fcs_resort_get_original_particles((*fmm_resort)->resort), addr, fcs_resort_get_indices((*fmm_resort)->resort),
#ifdef DO_TIMING
    &t[2],
#else
    NULL,
#endif
    size, rank, (*fmm_resort)->comm);

  TIMING_SYNC(comm); TIMING_STOP(t[1]);

  TIMING_SYNC(comm); TIMING_STOP(t[0]);
  
  TIMING_CMD(
    if (I_AM_MASTER) printf(TIMING_PRINT_PREFIX "mpi_fmm_resort_init: %f  %f  %f  %f\n", t[0], t[1], t[2], t[3]);
  );
#undef I_AM_MASTER
}

void mpi_fmm_resort_init_(fcs_fmm_resort_t *fmm_resort, pint_t *ntotal, pint_t *nlocal, back_idx_slkey_t *addr)
{
  mpi_fmm_resort_init(fmm_resort, ntotal, nlocal, addr);
}
