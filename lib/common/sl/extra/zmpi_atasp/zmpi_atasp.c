/*
 *  Copyright (C) 2011, 2012, 2013 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS.
 *  
 *  ScaFaCoS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  

 *  
 *  SL - Sorting Library, michael <dot> hofmann <at> informatik <dot> tu-chemnitz <dot> de
 */


#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "spec_core.h"
#include "zmpi_atasp.h"


struct _ZMPI_Tproc
{
  spec_tproc_t spec_tproc;

};


static void _ZMPI_Tproc_create(ZMPI_Tproc *tproc)
{
  *tproc = z_alloc(1, sizeof(*tproc));
}


static void _ZMPI_Tproc_free(ZMPI_Tproc *tproc)
{
  z_free(*tproc);
  
  *tproc = NULL;
}


int ZMPI_Tproc_create_tproc(ZMPI_Tproc *tproc, ZMPI_TPROC_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, ZMPI_Tproc_exdef exdef) /* zmpi_func ZMPI_Tproc_create_tproc */
{
  if (exdef != ZMPI_TPROC_EXDEF_NULL && exdef->type != 1) return 1;

  _ZMPI_Tproc_create(tproc);

  spec_tproc_create(&(*tproc)->spec_tproc, tfn, NULL, NULL, NULL);

  spec_tproc_set_reset(&(*tproc)->spec_tproc, rfn);

  if (exdef != ZMPI_TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tproc(&(*tproc)->spec_tproc, exdef->tproc_count_db, exdef->tproc_count_ip, exdef->tproc_rearrange_db, exdef->tproc_rearrange_ip);
  
  return MPI_SUCCESS;
}


int ZMPI_Tproc_create_tproc_mod(ZMPI_Tproc *tproc, ZMPI_TPROC_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, ZMPI_Tproc_exdef exdef) /* zmpi_func ZMPI_Tproc_create_tproc_mod */
{
  if (exdef != ZMPI_TPROC_EXDEF_NULL && exdef->type != 2) return 1;

  _ZMPI_Tproc_create(tproc);

  spec_tproc_create(&(*tproc)->spec_tproc, NULL, tfn, NULL, NULL);

  spec_tproc_set_reset(&(*tproc)->spec_tproc, rfn);

  if (exdef != ZMPI_TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tproc_mod(&(*tproc)->spec_tproc, exdef->tproc_mod_count_db, exdef->tproc_mod_count_ip, exdef->tproc_mod_rearrange_db, exdef->tproc_mod_rearrange_ip);
  
  return MPI_SUCCESS;
}


int ZMPI_Tproc_create_tprocs(ZMPI_Tproc *tproc, ZMPI_TPROCS_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, ZMPI_Tproc_exdef exdef) /* zmpi_func ZMPI_Tproc_create_tprocs */
{
  if (exdef != ZMPI_TPROC_EXDEF_NULL && exdef->type != 3) return 1;

  _ZMPI_Tproc_create(tproc);
  
  spec_tproc_create(&(*tproc)->spec_tproc, NULL, NULL, tfn, NULL);

  spec_tproc_set_reset(&(*tproc)->spec_tproc, rfn);

  if (exdef != ZMPI_TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tprocs(&(*tproc)->spec_tproc, exdef->tprocs_count_db, exdef->tprocs_count_ip, exdef->tprocs_rearrange_db, exdef->tprocs_rearrange_ip);
  
  return MPI_SUCCESS;
}


int ZMPI_Tproc_create_tprocs_mod(ZMPI_Tproc *tproc, ZMPI_TPROCS_MOD_FN *tfn, ZMPI_TPROC_RESET_FN *rfn, ZMPI_Tproc_exdef exdef) /* zmpi_func ZMPI_Tproc_create_tprocs_mod */
{
  if (exdef != ZMPI_TPROC_EXDEF_NULL && exdef->type != 4) return 1;

  _ZMPI_Tproc_create(tproc);
  
  spec_tproc_create(&(*tproc)->spec_tproc, NULL, NULL, NULL, tfn);

  spec_tproc_set_reset(&(*tproc)->spec_tproc, rfn);
  
  if (exdef != ZMPI_TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tprocs_mod(&(*tproc)->spec_tproc, exdef->tprocs_mod_count_db, exdef->tprocs_mod_count_ip, exdef->tprocs_mod_rearrange_db, exdef->tprocs_mod_rearrange_ip);
  
  return MPI_SUCCESS;
}


int ZMPI_Tproc_free(ZMPI_Tproc *tproc) /* zmpi_func ZMPI_Tproc_free */
{
  _ZMPI_Tproc_free(tproc);
  
  return MPI_SUCCESS;
}


int ZMPI_Tproc_set_neighbors(ZMPI_Tproc *tproc, int nneighbors, int *neighbors, MPI_Comm comm) /* zmpi_func ZMPI_Tproc_set_neighbors */
{
  int comm_size, comm_rank;

#ifdef SPEC_PROCLIST
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  spec_tproc_set_proclist(&(*tproc)->spec_tproc, nneighbors, neighbors, -1, NULL, comm_size, comm_rank, comm);
  
  return MPI_SUCCESS;
#else
  return 1;
#endif
}


int ZMPI_Alltoall_specific(void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype, ZMPI_Tproc tproc, void *tproc_data, int *received, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_specific */
{
  int exit_code;
  int comm_size, comm_rank;

  spec_elem_t sb, rb, *b, *xb;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if (sbuf == MPI_IN_PLACE && rbuf == MPI_IN_PLACE)
  {
#ifdef ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE
    fprintf(ZMPI_ALLTOALL_SPECIFIC_ERROR_FILE, "%d: MPI_Alltoall_specific_alltoallv_spec: error: either sbuf or rbuf can be MPI_IN_PLACE!\n", comm_rank);
#endif
    exit_code = 1;
    goto exit;
  }

  if (sbuf != MPI_IN_PLACE)
  {
    sb.buf = sbuf;
    sb.count = sb.max_count = scount;
    sb.mpi_type = stype;
    zmpil_create(&sb.zmpil_type, stype);
    
    b = &sb;
  }

  if (rbuf != MPI_IN_PLACE)
  {
    rb.buf = rbuf;
    rb.count = rb.max_count = rcount;
    rb.mpi_type = rtype;
    zmpil_create(&rb.zmpil_type, rtype);

    b = &rb;
  }

  xb = NULL;

  if (sbuf == MPI_IN_PLACE || rbuf == MPI_IN_PLACE)
  {
    b->count = scount;
    b->max_count = rcount;

    exit_code = spec_alltoallv_ip(b, xb, tproc->spec_tproc, tproc_data, comm_size, comm_rank, comm);
    *received = b->count;

  } else {

    exit_code = spec_alltoallv_db(&sb, &rb, xb, tproc->spec_tproc, tproc_data, comm_size, comm_rank, comm);
    *received = rb.count;
  }

  zmpil_destroy(&sb.zmpil_type);

  zmpil_destroy(&rb.zmpil_type);

exit:

  return exit_code;  
}
