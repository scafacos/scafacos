/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015 Michael Hofmann
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
#include "spec_common.h"


#ifndef SPEC_ALLTOALLV_TRACE_IF
# define SPEC_ALLTOALLV_TRACE_IF  (rank == -1)
#endif


/*#define SPEC_PRINT*/


spint_t spec_alltoallv_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_alltoallv_db */
{
  spint_t exit_code = SPEC_EXIT_SUCCESS;

  spint_t i;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;

  int *scounts, *sdispls, *rcounts, *rdispls;

  spec_elem_t _xb;

  spint_t stotal, rtotal;

  SPEC_DECLARE_TPROC_REARRANGE_DB
  SPEC_DECLARE_TPROC_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_REARRANGE_DB
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB

#ifdef Z_PACK_TIMING
  double t[5] = { 0, 0, 0, 0, 0 };
#endif


  Z_TRACE_IF(SPEC_ALLTOALLV_TRACE_IF, "spec_alltoallv_db");

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  /* check supported tproc functions */
  if (!spec_check_tproc_support(tproc, 1, 1, 1, 1, rank, "spec_alltoallv_db"))
  {
    spec_elem_set_n(rb, 0);
    exit_code = SPEC_EXIT_FAILED;
    goto exit;
  }

  scounts = z_alloc(4 * size, sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

  /* setup tproc buffers */
  spec_tproc_setup(tproc, sb, &procs, &mods);

#ifdef SPEC_PRINT
  printf("input\n");
  spec_print(tproc, tproc_data, sb);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  spec_make_counts(tproc, tproc_data, sb, 0, size, scounts, procs);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

#ifdef SPEC_PRINT
  printf("after count\n");
  spec_print(tproc, tproc_data, sb);
#endif

  /* redistribute local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[2]);

  spec_redistribute_counts(scounts, rcounts,
#ifdef SPEC_PROCLISTS
    tproc->nsend_procs, tproc->send_procs, tproc->nrecv_procs, tproc->recv_procs,
#endif
    size, rank, comm);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[2]);

  sdispls[0] = rdispls[0] = 0;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  }

  stotal = sdispls[size - 1] + scounts[size - 1];
  rtotal = rdispls[size - 1] + rcounts[size - 1];

  /* check size of receive buffer */
  if (!spec_check_buffer_size(rb, rtotal, 1, comm, rank, "spec_alltoallv_db", "receive"))
  {
    spec_elem_set_n(rb, rtotal);
    exit_code = SPEC_EXIT_FAILED;
    goto free_and_exit;
  }

  if (xb == NULL || spec_elem_get_nmax(xb) < stotal)
  {
    spec_elem_copy_type(sb, &_xb);
    spec_elem_alloc_tmp(&_xb, stotal);

    xb = &_xb;
  }

  spec_elem_copy_type(sb, xb);

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* local rearrange */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[3]);

  if (tproc->tproc)
  {
    if (tproc->tproc_ext.rearrange_db) tproc->tproc_ext.rearrange_db(tproc_data, sb, xb, sdispls);
    else SPEC_DO_TPROC_REARRANGE_DB(tproc->tproc, tproc_data, sb, xb, sdispls);

  } else if (tproc->tproc_mod)
  {
    if (tproc->tproc_mod_ext.rearrange_db) tproc->tproc_mod_ext.rearrange_db(tproc_data, sb, xb, sdispls, mods);
    else SPEC_DO_TPROC_MOD_REARRANGE_DB(tproc->tproc_mod, tproc_data, sb, xb, sdispls, mods);

  } else if (tproc->tprocs)
  {
    if (tproc->tprocs_ext.rearrange_db) tproc->tprocs_ext.rearrange_db(tproc_data, sb, xb, sdispls, procs);
    else SPEC_DO_TPROCS_REARRANGE_DB(tproc->tprocs, tproc_data, sb, xb, sdispls, procs);

  } else if (tproc->tprocs_mod)
  {
    if (tproc->tprocs_mod_ext.rearrange_db) tproc->tprocs_mod_ext.rearrange_db(tproc_data, sb, xb, sdispls, procs, mods);
    else SPEC_DO_TPROCS_MOD_REARRANGE_DB(tproc->tprocs_mod, tproc_data, sb, xb, sdispls, procs, mods);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

/*#ifdef SPEC_PRINT
  printf("after rearrange\n");
  spec_print(tproc, tproc_data, xb);
#endif*/

  /* redistribute with alltoallv */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[4]);

  sdispls[0] = 0;
  for (i = 1; i < size; ++i) sdispls[i] = sdispls[i - 1] + scounts[i - 1];

#ifdef SPEC_PROCLISTS
# ifdef spec_elem_alltoallv_proclists_db
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    spec_elem_alltoallv_proclists_db(xb, scounts, sdispls, tproc->nsend_procs, tproc->send_procs, rb, rcounts, rdispls, tproc->nrecv_procs, tproc->recv_procs, size, rank, comm);

  } else
# endif
#endif
  {
    spec_elem_alltoallv_db(xb, scounts, sdispls, rb, rcounts, rdispls, size, rank, comm);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  spec_elem_set_n(rb, rtotal);

  if (xb == &_xb) spec_elem_free_tmp(&_xb);

free_and_exit:

  /* free tproc buffers */
  spec_tproc_release(&procs, &mods);

  z_free(scounts);

exit:
  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[0]);

  Z_TIMING_PRINT(0, __func__, sizeof(t) / sizeof(t[0]), t, rank);

#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
  if (spec_timing)
    for (i = 0; i < sizeof(t) / sizeof(t[0]); ++i) spec_timing[i] = t[i];
#endif

  return exit_code;
}


spint_t spec_alltoallv_ip(spec_elem_t *b, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_alltoallv_ip */
{
  spint_t exit_code = SPEC_EXIT_SUCCESS;

  spint_t i;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;

  int *scounts, *sdispls, *rcounts, *rdispls;

  spec_elem_t _xb, *tb;

  spint_t stotal, rtotal;

  spint_t local_alltoallv_inplace, global_alltoallv_inplace;

  SPEC_DECLARE_TPROC_REARRANGE_IP
  SPEC_DECLARE_TPROC_REARRANGE_DB
  SPEC_DECLARE_TPROC_MOD_REARRANGE_IP
  SPEC_DECLARE_TPROC_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_REARRANGE_DB
  SPEC_DECLARE_TPROCS_REARRANGE_IP
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP

#ifdef Z_PACK_TIMING
  double t[5] = { 0, 0, 0, 0, 0 };
#endif


  Z_TRACE_IF(SPEC_ALLTOALLV_TRACE_IF, "spec_alltoallv_ip");

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  /* check supported tproc functions */
  if (!spec_check_tproc_support(tproc, 1, 1, 1, 1, rank, "spec_alltoallv_ip"))
  {
    spec_elem_set_n(b, 0);
    exit_code = SPEC_EXIT_FAILED;
    goto exit;
  }

  scounts = z_alloc(4 * size, sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

  /* setup tproc buffers */
  spec_tproc_setup(tproc, b, &procs, &mods);

#ifdef SPEC_PRINT
  printf("input\n");
  spec_print(tproc, tproc_data, b);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  spec_make_counts(tproc, tproc_data, b, 1, size, scounts, procs);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

#ifdef SPEC_PRINT
  printf("after count\n");
  spec_print(tproc, tproc_data, b);
#endif

  /* redistribute local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[2]);

  spec_redistribute_counts(scounts, rcounts,
#ifdef SPEC_PROCLISTS
    tproc->nsend_procs, tproc->send_procs, tproc->nrecv_procs, tproc->recv_procs,
#endif
    size, rank, comm);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[2]);

  sdispls[0] = rdispls[0] = 0;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  }

  stotal = sdispls[size - 1] + scounts[size - 1];
  rtotal = rdispls[size - 1] + rcounts[size - 1];

  /* check size of buffer */
  if (!spec_check_buffer_size(b, z_max(stotal, rtotal), 0, comm, rank, "spec_alltoallv_ip", ""))
  {
    spec_elem_set_n(b, z_max(stotal, rtotal));
    exit_code = SPEC_EXIT_FAILED;
    goto free_and_exit;
  }

  if (xb == NULL)
  {
    spec_elem_copy_type(b, &_xb);
    spec_elem_alloc_tmp(&_xb, 1);

    xb = &_xb;
  }

  tb = NULL;
  if (spec_elem_get_nmax(xb) >= stotal) tb = xb;

  spec_elem_copy_type(b, xb);

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* local rearrange */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[3]);

  if (tproc->tproc)
  {
    if (tb)
    {
      if (tproc->tproc_ext.rearrange_db) tproc->tproc_ext.rearrange_db(tproc_data, b, tb, sdispls);
      else SPEC_DO_TPROC_REARRANGE_DB(tproc->tproc, tproc_data, b, tb, sdispls);

    } else
    {
      if (tproc->tproc_ext.rearrange_ip) tproc->tproc_ext.rearrange_ip(tproc_data, b, xb, sdispls, scounts, size);
      else SPEC_DO_TPROC_REARRANGE_IP(tproc->tproc, tproc_data, b, xb, sdispls, scounts, size);
    }

  } else if (tproc->tproc_mod)
  {
    if (tb)
    {
      if (tproc->tproc_mod_ext.rearrange_db) tproc->tproc_mod_ext.rearrange_db(tproc_data, b, tb, sdispls, mods);
      else SPEC_DO_TPROC_MOD_REARRANGE_DB(tproc->tproc_mod, tproc_data, b, tb, sdispls, mods);

    } else
    {
      if (tproc->tproc_mod_ext.rearrange_ip) tproc->tproc_mod_ext.rearrange_ip(tproc_data, b, xb, sdispls, scounts, size, mods);
      else SPEC_DO_TPROC_MOD_REARRANGE_IP(tproc->tproc_mod, tproc_data, b, xb, sdispls, scounts, size, mods);
    }

  } else if (tproc->tprocs)
  {
    if (tb)
    {
      if (tproc->tprocs_ext.rearrange_db) tproc->tprocs_ext.rearrange_db(tproc_data, b, tb, sdispls, procs);
      else SPEC_DO_TPROCS_REARRANGE_DB(tproc->tprocs, tproc_data, b, tb, sdispls, procs);

    } else
    {
      if (tproc->tprocs_ext.rearrange_ip) tproc->tprocs_ext.rearrange_ip(tproc_data, b, xb, sdispls, scounts, size, procs);
      else SPEC_DO_TPROCS_REARRANGE_IP(tproc->tprocs, tproc_data, b, xb, sdispls, scounts, size, procs);
    }

  } else if (tproc->tprocs_mod)
  {
    if (tb)
    {
      if (tproc->tprocs_mod_ext.rearrange_db) tproc->tprocs_mod_ext.rearrange_db(tproc_data, b, tb, sdispls, procs, mods);
      else SPEC_DO_TPROCS_MOD_REARRANGE_DB(tproc->tprocs_mod, tproc_data, b, tb, sdispls, procs, mods);

    } else
    {
      if (tproc->tprocs_mod_ext.rearrange_ip) tproc->tprocs_mod_ext.rearrange_ip(tproc_data, b, xb, sdispls, scounts, size, procs, mods);
      else SPEC_DO_TPROCS_MOD_REARRANGE_IP(tproc->tprocs_mod, tproc_data, b, xb, sdispls, scounts, size, procs, mods);
    }
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

#ifdef SPEC_PRINT
  printf("after rearrange\n");
  spec_print(tproc, tproc_data, b);
#endif

  /* redistribute with alltoallv */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[4]);

  sdispls[0] = 0;
  for (i = 1; i < size; ++i) sdispls[i] = sdispls[i - 1] + scounts[i - 1];

  if (tb) local_alltoallv_inplace = 0;
  else
  {
    if (spec_elem_get_nmax(xb) >= rtotal) local_alltoallv_inplace = 0;
    else local_alltoallv_inplace = 1;
  }

  MPI_Allreduce(&local_alltoallv_inplace, &global_alltoallv_inplace, 1, MPI_SPINT, MPI_SUM, comm);

  /* no process wants inplace? */
  if (global_alltoallv_inplace == 0)
  {
    /* double-buffered alltoallv  */
    if (tb)
    {
#ifdef SPEC_PROCLISTS
# ifdef spec_elem_alltoallv_proclists_db
      if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
      {
        spec_elem_alltoallv_proclists_db(tb, scounts, sdispls, tproc->nsend_procs, tproc->send_procs, b, rcounts, rdispls, tproc->nrecv_procs, tproc->recv_procs, size, rank, comm);

      } else
# endif
#endif
      {
        spec_elem_alltoallv_db(tb, scounts, sdispls, b, rcounts, rdispls, size, rank, comm);
      }
    } else
    {
#ifdef SPEC_PROCLISTS
# ifdef spec_elem_alltoallv_proclists_db
      if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
      {
        spec_elem_alltoallv_proclists_db(b, scounts, sdispls, tproc->nsend_procs, tproc->send_procs, xb, rcounts, rdispls, tproc->nrecv_procs, tproc->recv_procs, size, rank, comm);

      } else
# endif
#endif
      {
        spec_elem_alltoallv_db(b, scounts, sdispls, xb, rcounts, rdispls, size, rank, comm);
      }
      spec_elem_ncopy_at(xb, 0, b, 0, rtotal);
    }

  } else
  {
    if (tb) spec_elem_ncopy_at(tb, 0, b, 0, stotal);
#ifdef SPEC_PROCLISTS
# ifdef spec_elem_alltoallv_proclists_ip
    if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
    {
      spec_elem_alltoallv_proclists_ip(b, xb, scounts, sdispls, tproc->nsend_procs, tproc->send_procs, rcounts, rdispls, tproc->nrecv_procs, tproc->recv_procs, size, rank, comm);

    } else
# endif
#endif
    {
#ifdef spec_elem_alltoallv_ip
      spec_elem_alltoallv_ip(b, xb, scounts, sdispls, rcounts, rdispls, size, rank, comm);
#else
#ifdef SPEC_ERROR_FILE
      fprintf(SPEC_ERROR_FILE, "%d: spec_alltoallv_ip: error: required in-place support for alltoallv not available", rank);
#endif
      exit_code = SPEC_EXIT_FAILED;
      goto free_and_exit;
#endif
    }
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  spec_elem_set_n(b, rtotal);

  if (xb == &_xb) spec_elem_free_tmp(&_xb);

free_and_exit:

  /* free tproc buffers */
  spec_tproc_release(&procs, &mods);

  z_free(scounts);

exit:
  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[0]);

  Z_TIMING_PRINT(0, __func__, sizeof(t) / sizeof(t[0]), t, rank);

#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
  if (spec_timing)
    for (i = 0; i < sizeof(t) / sizeof(t[0]); ++i) spec_timing[i] = t[i];
#endif

  return exit_code;
}



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#ifdef HAVE_ZMPI_TOOLS_H
# include "zmpi_tools.h"
#endif

#include "spec_core.h"
#include "spec_common.h"


#ifndef SPEC_COMMON_TRACE_IF
# define SPEC_COMMON_TRACE_IF  (rank == -1)
#endif


#ifdef SPEC_ERROR_FILE
static const char *spec_tproc_name(spint_t id)
{
  switch (id)
  {
    case 0: return "none";
    case 1: return "tproc";
    case 2: return "tproc_mod";
    case 3: return "tprocs";
    case 4: return "tprocs_mod";
  }

  return "unsupported";
}
#endif


static spint_t spec_tproc_supported(spec_tproc_t tproc, spint_t have_tproc, spint_t have_tproc_mod, spint_t have_tprocs, spint_t have_tprocs_mod, spint_t *id)
{
  if (tproc != SPEC_TPROC_NULL)
  {
    if (tproc->tproc)
    {
      *id = 1;
      return have_tproc;

    } else if (tproc->tproc_mod)
    {
      *id = 2;
      return have_tproc_mod;

    } else if (tproc->tprocs)
    {
      *id = 3;
      return have_tprocs;

    } else if (tproc->tprocs_mod)
    {
      *id = 4;
      return have_tprocs_mod;
    }
  }

  *id = 0;
  return 0;
}


spint_t spec_check_tproc_support(spec_tproc_t tproc, spint_t have_tproc, spint_t have_tproc_mod, spint_t have_tprocs, spint_t have_tprocs_mod, int rank, const char *name) /* sp_func spec_check_tproc_support */
{
  spint_t id;

  if (!spec_tproc_supported(tproc, have_tproc, have_tproc_mod, have_tprocs, have_tprocs_mod, &id))
  {
#ifdef SPEC_ERROR_FILE
    fprintf(SPEC_ERROR_FILE, "%d: spec_alltoallv_db: error: target process function type '%s' is not supported\n", rank, spec_tproc_name(id));
#endif
    return 0;
  }

  return 1;
}


spint_t spec_check_buffer_size(spec_elem_t *b, spint_t min_size, spint_t allocatable, MPI_Comm comm, int rank, const char *name, const char *buf_name) /* sp_func spec_check_buffer_size */
{
  spint_t local_buffer_exit, global_buffer_exit;


  local_buffer_exit = (min_size > spec_elem_get_nmax(b))?1:0;

#ifdef spec_elem_alloc_rbuf
  if (local_buffer_exit > 0 && spec_elem_alloc_rbuf(b) && allocatable)
  {
    spec_elem_free_buf(b);
    spec_elem_alloc_buf(b, min_size);
    local_buffer_exit = 0;
  }
#endif

#ifdef SPEC_GLOBAL_EXIT_ON_ERROR
  MPI_Allreduce(&local_buffer_exit, &global_buffer_exit, 1, MPI_SPINT, MPI_SUM, comm);
#else
  global_buffer_exit = local_buffer_exit;
#endif

  if (global_buffer_exit > 0)
  {
#ifdef SPEC_ERROR_FILE
    fprintf(SPEC_ERROR_FILE,
# ifdef SPEC_GLOBAL_EXIT_ON_ERROR
      "%d: %s: error: one %s%sbuffer is too small (local: %" spint_fmt " vs. %" spint_fmt")\n",
# else
      "%d: %s: error: local %s%sbuffer too small (%" spint_fmt " vs. %" spint_fmt")\n",
# endif
      rank, name, buf_name, ((strlen(buf_name) > 0)?" ":""), (spint_t) spec_elem_get_nmax(b), min_size);
#endif
    return 0;
  }

  return 1;
}


void spec_tproc_setup(spec_tproc_t tproc, spec_elem_t *b, spec_proc_t **procs, spec_elem_t **mods) /* sp_func spec_tproc_setup */
{
  spint_t tprocs_max;
  spint_t alloc_mods;


  tprocs_max = 0;
  alloc_mods = 0;

  if (tproc->tproc)
  {
    /* nothing */

  } else if (tproc->tproc_mod)
  {
    alloc_mods = 1;

  } else if (tproc->tprocs)
  {
    tprocs_max = tproc->max_tprocs;

  } else if (tproc->tprocs_mod)
  {
    tprocs_max = tproc->max_tprocs;
    alloc_mods = 1;
  }

  if (procs)
  {
    if (tprocs_max) *procs = z_alloc(tprocs_max, sizeof(spec_proc_t));
    else *procs = NULL;
  }

  if (mods)
  {
    if (alloc_mods)
    {
      *mods = z_alloc(1, sizeof(spec_elem_t));

      spec_elem_copy_type(b, *mods);
      spec_elem_alloc_buf(*mods, 2 * z_max(1, tprocs_max));

    } else *mods = NULL;
  }
}


void spec_tproc_release(spec_proc_t **procs, spec_elem_t **mods) /* sp_func spec_tproc_release */
{
  if (procs && *procs) z_free(*procs);

  if (mods && *mods)
  {
    spec_elem_free_buf(*mods);
    z_free(*mods);
  }
}


spint_t spec_make_counts(spec_tproc_t tproc, spec_tproc_data_t tproc_data, spec_elem_t *b, int ip, int size, int *counts, spec_proc_t *procs) /* sp_func spec_make_counts */
{
  spint_t i;

  SPEC_DECLARE_TPROC_COUNT_DB
  SPEC_DECLARE_TPROC_MOD_COUNT_DB
  SPEC_DECLARE_TPROCS_COUNT_DB
  SPEC_DECLARE_TPROCS_MOD_COUNT_DB

  SPEC_DECLARE_TPROC_COUNT_IP
  SPEC_DECLARE_TPROC_MOD_COUNT_IP
  SPEC_DECLARE_TPROCS_COUNT_IP
  SPEC_DECLARE_TPROCS_MOD_COUNT_IP

  for (i = 0; i < size; ++i) counts[i] = 0;

  if (!ip)
  {
    if (tproc->tproc)
    {
      if (tproc->tproc_ext.count_db) tproc->tproc_ext.count_db(tproc_data, b, counts);
      else SPEC_DO_TPROC_COUNT_DB(tproc->tproc, tproc_data, b, counts);

    } else if (tproc->tproc_mod)
    {
      if (tproc->tproc_mod_ext.count_db) tproc->tproc_mod_ext.count_db(tproc_data, b, counts);
      else SPEC_DO_TPROC_MOD_COUNT_DB(tproc->tproc_mod, tproc_data, b, counts);

    } else if (tproc->tprocs)
    {
      if (tproc->tprocs_ext.count_db) tproc->tprocs_ext.count_db(tproc_data, b, counts, procs);
      else SPEC_DO_TPROCS_COUNT_DB(tproc->tprocs, tproc_data, b, counts, procs);

    } else if (tproc->tprocs_mod)
    {
      if (tproc->tprocs_mod_ext.count_db) tproc->tprocs_mod_ext.count_db(tproc_data, b, counts, procs);
      else SPEC_DO_TPROCS_MOD_COUNT_DB(tproc->tprocs_mod, tproc_data, b, counts, procs);
    }

  } else {

    if (tproc->tproc)
    {
      if (tproc->tproc_ext.count_ip) tproc->tproc_ext.count_ip(tproc_data, b, counts);
      else SPEC_DO_TPROC_COUNT_IP(tproc->tproc, tproc_data, b, counts);

    } else if (tproc->tproc_mod)
    {
      if (tproc->tproc_mod_ext.count_ip) tproc->tproc_mod_ext.count_ip(tproc_data, b, counts);
      else SPEC_DO_TPROC_MOD_COUNT_IP(tproc->tproc_mod, tproc_data, b, counts);

    } else if (tproc->tprocs)
    {
      if (tproc->tprocs_ext.count_ip) tproc->tprocs_ext.count_ip(tproc_data, b, counts, procs);
      else SPEC_DO_TPROCS_COUNT_IP(tproc->tprocs, tproc_data, b, counts, procs);

    } else if (tproc->tprocs_mod)
    {
      if (tproc->tprocs_mod_ext.count_ip) tproc->tprocs_mod_ext.count_ip(tproc_data, b, counts, procs);
      else SPEC_DO_TPROCS_MOD_COUNT_IP(tproc->tprocs_mod, tproc_data, b, counts, procs);
    }
  }

  return 0;
}


/* sp_var spec_redistribute_counts_type spec_redistribute_counts_proclists_type */
spint_t spec_redistribute_counts_type = SPEC_REDISTRIBUTE_COUNTS_DEFAULT;
spint_t spec_redistribute_counts_proclists_type = SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_DEFAULT;


spint_t spec_redistribute_counts(int *scounts, int *rcounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm) /* sp_func spec_redistribute_counts */
{
  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "spec_redistribute_counts");

#ifdef HAVE_ZMPI_ALLTOALL_INT

#ifdef SPEC_PROCLISTS
  if ((nsend_procs >= 0 || nrecv_procs >= 0) && spec_redistribute_counts_proclists_type != SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_IGNORE)
  {
    memset(rcounts, 0, size * sizeof(int));

    switch (spec_redistribute_counts_proclists_type)
    {
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_DEFAULT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_DEFAULT -> isendirecv");
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ISENDIRECV:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ISENDIRECV");
        ZMPI_Alltoall_int_proclists_isendirecv(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ALLTOALLV:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ALLTOALLV");
        ZMPI_Alltoall_int_proclists_alltoallv(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT");
        ZMPI_Alltoall_int_proclists_put(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_ALLOC");
        ZMPI_Alltoall_int_proclists_put_alloc(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES");
        ZMPI_Alltoall_int_proclists_put_2phases(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES_ALLOC");
        ZMPI_Alltoall_int_proclists_put_2phases_alloc(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
    }

  } else
#endif
  {
    switch (spec_redistribute_counts_type)
    {
      case SPEC_REDISTRIBUTE_COUNTS_DEFAULT:
#ifdef SPEC_REDISTRIBUTE_COUNTS_2STEP_THRESHOLD
        if (size >= SPEC_REDISTRIBUTE_COUNTS_2STEP_THRESHOLD)
        {
          Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_DEFAULT -> 2step");
          ZMPI_Alltoall_int_2step(scounts, rcounts, comm);

        } else
#endif
        {
          Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_DEFAULT -> alltoall");
          ZMPI_Alltoall_int_alltoall(scounts, rcounts, comm);
        }
        break;
      case SPEC_REDISTRIBUTE_COUNTS_ALLTOALL:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_ALLTOALL");
        ZMPI_Alltoall_int_alltoall(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_2STEP:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_2STEP");
        ZMPI_Alltoall_int_2step(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT");
        ZMPI_Alltoall_int_put(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_ALLOC");
        ZMPI_Alltoall_int_put_alloc(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES");
        ZMPI_Alltoall_int_put_2phases(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES_ALLOC");
        ZMPI_Alltoall_int_put_2phases_alloc(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES");
        ZMPI_Alltoall_int_put_3phases(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES_ALLOC");
        ZMPI_Alltoall_int_put_3phases_alloc(scounts, rcounts, comm);
        break;
    }
  }

#else /* HAVE_ZMPI_ALLTOALL_INT */

  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "MPI_Alltoall");
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);

#endif /* HAVE_ZMPI_ALLTOALL_INT */

  return 0;
}


/* sp_var spec_reduce_scatter_counts_type spec_reduce_scatter_counts_proclists_type */
spint_t spec_reduce_scatter_counts_type = SPEC_REDUCE_SCATTER_COUNTS_DEFAULT;
spint_t spec_reduce_scatter_counts_proclists_type = SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_DEFAULT;


spint_t spec_reduce_scatter_counts(int *scounts, int *rcounts, int ncounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm) /* sp_func spec_reduce_scatter_counts */
{
  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "spec_reduce_scatter_counts");

#ifdef HAVE_ZMPI_REDUCE_SCATTER_BLOCK_INTSUM

#ifdef SPEC_PROCLISTS
  if ((nsend_procs >= 0 || nrecv_procs >= 0) && spec_reduce_scatter_counts_proclists_type != SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_IGNORE)
  {
    memset(rcounts, 0, ncounts * sizeof(int));

    switch (spec_reduce_scatter_counts_proclists_type)
    {
      case SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_DEFAULT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_DEFAULT -> isendirecv");
      case SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ISENDIRECV:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ISENDIRECV");
        ZMPI_Reduce_scatter_block_intsum_proclists_isendirecv(scounts, nsend_procs, send_procs, rcounts, ncounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ALLTOALLV:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ALLTOALLV");
        ZMPI_Reduce_scatter_block_intsum_proclists_alltoallv(scounts, nsend_procs, send_procs, rcounts, ncounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ACCUMULATE:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ACCUMULATE");
        ZMPI_Reduce_scatter_block_intsum_proclists_accumulate(scounts, nsend_procs, send_procs, rcounts, ncounts, nrecv_procs, recv_procs, comm);
        break;
    }

  } else
#endif
  {
    switch (spec_reduce_scatter_counts_type)
    {
      case SPEC_REDUCE_SCATTER_COUNTS_DEFAULT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_DEFAULT -> redscat");
      case SPEC_REDUCE_SCATTER_COUNTS_REDSCAT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_REDSCAT");
        ZMPI_Reduce_scatter_block(scounts, rcounts, ncounts, MPI_INT, MPI_SUM, comm);
        break;
      case SPEC_REDUCE_SCATTER_COUNTS_ACCUMULATE:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_ACCUMULATE");
        ZMPI_Reduce_scatter_block_intsum_accumulate(scounts, rcounts, ncounts, comm);
        break;
    }
  }

#elif defined(HAVE_ZMPI_REDUCE_SCATTER_BLOCK)

  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "ZMPI_Reduce_scatter_block");
  ZMPI_Reduce_scatter_block(scounts, rcounts, ncounts, MPI_INT, MPI_SUM, comm);

#else

#ifdef SPEC_ERROR_FILE
  fprintf(SPEC_ERROR_FILE, "%d: spec_reduce_scatter_counts: error: no implementation available\n", rank);
#endif

#endif

  return 0;
}


/* sp_var spec_prefix_counts_type spec_prefix_counts_proclists_type */
spint_t spec_prefix_counts_type = SPEC_PREFIX_COUNTS_DEFAULT;
spint_t spec_prefix_counts_proclists_type = SPEC_PREFIX_COUNTS_PROCLISTS_DEFAULT;


spint_t spec_prefix_counts(int *scounts, int *rcounts, int ncounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm) /* sp_func spec_prefix_counts */
{
  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "spec_prefix_counts");

#ifdef SPEC_PROCLISTS
  if ((nsend_procs >= 0 || nrecv_procs >= 0) && spec_prefix_counts_proclists_type != SPEC_PREFIX_COUNTS_PROCLISTS_IGNORE)
  {
    memset(rcounts, 0, ncounts * sizeof(int));

    switch (spec_prefix_counts_proclists_type)
    {
    }

  } else
#endif
  {
    switch (spec_prefix_counts_type)
    {
      case SPEC_PREFIX_COUNTS_SCAN:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_PREFIX_COUNTS_SCAN");
        MPI_Scan(scounts, rcounts, ncounts, MPI_INT, MPI_SUM, comm);
        break;
    }
  }

  return 0;
}



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "spec_core.h"
#include "spec_common.h"


#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
double *spec_timing = NULL; /* sp_var spec_timing */
#endif


/*#define DUMMY_TEST*/

#ifdef DUMMY_TEST

static int dummy_spec_tproc(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data)
{
  return SPEC_PROC_NONE;
}

SPEC_DEFINE_TPROC(dummy, dummy_spec_tproc);

static int dummy_spec_tproc_mod(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, spec_elem_buf_t mod)
{
  return SPEC_PROC_NONE;
}

SPEC_DEFINE_TPROC_MOD(dummy, dummy_spec_tproc_mod);

static int dummy_spec_tprocs(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, int *procs)
{
  return 0;
}

SPEC_DEFINE_TPROCS(dummy, dummy_spec_tprocs);

static int dummy_spec_tprocs_mod(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, int *procs, spec_elem_buf_t mod)
{
  return 0;
}

SPEC_DEFINE_TPROCS_MOD(dummy, dummy_spec_tprocs_mod);

#endif


static void spec_tproc_unset_tproc(spec_tproc_t tproc)
{
  tproc->tproc = NULL;
}


static void spec_tproc_unset_ext_tproc(spec_tproc_t tproc)
{
  spec_tproc_ext_t tproc_ext = SPEC_EXT_PARAM_TPROC_NULL;

  tproc->tproc_ext = tproc_ext;
}


static void spec_tproc_unset_tproc_mod(spec_tproc_t tproc)
{
  tproc->tproc_mod = NULL;
}


static void spec_tproc_unset_ext_tproc_mod(spec_tproc_t tproc)
{
  spec_tproc_mod_ext_t tproc_mod_ext = SPEC_EXT_PARAM_TPROC_MOD_NULL;

  tproc->tproc_mod_ext = tproc_mod_ext;
}


static void spec_tproc_unset_tprocs(spec_tproc_t tproc)
{
  tproc->tprocs = NULL;
}


static void spec_tproc_unset_ext_tprocs(spec_tproc_t tproc)
{
  spec_tprocs_ext_t tprocs_ext = SPEC_EXT_PARAM_TPROCS_NULL;

  tproc->tprocs_ext = tprocs_ext;
}


static void spec_tproc_unset_tprocs_mod(spec_tproc_t tproc)
{
  tproc->tprocs_mod = NULL;
}


static void spec_tproc_unset_ext_tprocs_mod(spec_tproc_t tproc)
{
  spec_tprocs_mod_ext_t tprocs_mod_ext = SPEC_EXT_PARAM_TPROCS_MOD_NULL;

  tproc->tprocs_mod_ext = tprocs_mod_ext;
}


spint_t spec_tproc_create(spec_tproc_t *tproc, spec_tproc_f *func, spec_tproc_mod_f *func_mod, spec_tprocs_f *func_s, spec_tprocs_mod_f *func_s_mod, spint_t max_tprocs) /* sp_func spec_tproc_create */
{
  *tproc = z_alloc(1, sizeof(struct _spec_tproc_t));

  (*tproc)->max_tprocs = max_tprocs;

  (*tproc)->tproc = func;
  spec_tproc_unset_ext_tproc(*tproc);

  (*tproc)->tproc_mod = func_mod;
  spec_tproc_unset_ext_tproc_mod(*tproc);

  (*tproc)->tprocs = func_s;
  spec_tproc_unset_ext_tprocs(*tproc);

  (*tproc)->tprocs_mod = func_s_mod;
  spec_tproc_unset_ext_tprocs_mod(*tproc);

  (*tproc)->reset = NULL;

#ifdef SPEC_PROCLISTS
  (*tproc)->nsend_procs = -1;
  (*tproc)->send_procs = NULL;
  (*tproc)->nrecv_procs = -1;
  (*tproc)->recv_procs = NULL;
#endif

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_destroy(spec_tproc_t *tproc) /* sp_func spec_tproc_destroy */
{
#ifdef SPEC_PROCLISTS
  spec_tproc_set_proclists(*tproc, -1, NULL, -1, NULL, 0, -1, MPI_COMM_NULL);
#endif

  z_free(*tproc);

  *tproc = NULL;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_duplicate(spec_tproc_t tproc, spec_tproc_t *newtproc) /* sp_func spec_tproc_duplicate */
{
  *newtproc = z_alloc(1, sizeof(struct _spec_tproc_t));

  memcpy(*newtproc, tproc, sizeof(struct _spec_tproc_t));

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tproc(spec_tproc_t tproc, spec_tproc_f *func) /* sp_func spec_tproc_set_tproc */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  tproc->tproc = func;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tproc(spec_tproc_t tproc, const spec_tproc_ext_t *tproc_ext) /* sp_func spec_tproc_set_ext_tproc */
{
  tproc->tproc_ext = *tproc_ext;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tproc_mod(spec_tproc_t tproc, spec_tproc_mod_f *func_mod) /* sp_func spec_tproc_set_tproc_mod */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  tproc->tproc_mod = func_mod;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tproc_mod(spec_tproc_t tproc, const spec_tproc_mod_ext_t *tproc_mod_ext) /* sp_func spec_tproc_set_ext_tproc_mod */
{
  tproc->tproc_mod_ext = *tproc_mod_ext;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tprocs(spec_tproc_t tproc, spec_tprocs_f *func_s, spint_t max_tprocs) /* sp_func spec_tproc_set_tprocs */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  tproc->max_tprocs = max_tprocs;

  tproc->tprocs = func_s;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tprocs(spec_tproc_t tproc, const spec_tprocs_ext_t *tprocs_ext) /* sp_func spec_tproc_set_ext_tprocs */
{
  tproc->tprocs_ext = *tprocs_ext;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tprocs_mod(spec_tproc_t tproc, spec_tprocs_mod_f *func_s_mod, spint_t max_tprocs) /* sp_func spec_tproc_set_tprocs_mod */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  tproc->max_tprocs = max_tprocs;

  tproc->tprocs_mod = func_s_mod;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tprocs_mod(spec_tproc_t tproc, const spec_tprocs_mod_ext_t *tprocs_mod_ext) /* sp_func spec_tproc_set_ext_tprocs_mod */
{
  tproc->tprocs_mod_ext = *tprocs_mod_ext;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_reset(spec_tproc_t tproc, spec_tproc_reset_f *reset) /* sp_func spec_tproc_set_reset */
{
  tproc->reset = reset;

  return SPEC_EXIT_SUCCESS;
}


#ifdef SPEC_PROCLISTS

void spec_make_recv_proclist(spint_t nsend_procs, sproc_t *send_procs, spint_t *nrecv_procs, sproc_t **recv_procs, int size, int rank, MPI_Comm comm) /* sp_func spec_make_recv_proclist */
{
  spint_t i, j;
  int *s, *r;


  s = z_alloc(2 * size, sizeof(int));
  r = s + 1 * size;

  for (i = 0; i < size; ++i) s[i] = 0;
  for (i = 0; i < nsend_procs; ++i) s[send_procs[i]] = 1;

  MPI_Alltoall(s, 1, MPI_INT, r, 1, MPI_INT, comm);

  for (j = 0, i = 0; i < size; ++i) j += r[i];

  *nrecv_procs = j;
  *recv_procs = z_alloc(*nrecv_procs, sizeof(sproc_t));

  for (j = 0, i = 0; i < size; ++i)
  if (r[i])
  {
    (*recv_procs)[j] = i;
    ++j;
  }

  z_free(s);
}


spint_t spec_tproc_set_proclists(spec_tproc_t tproc, spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs, int size, int rank, MPI_Comm comm) /* sp_func spec_tproc_set_proclists */
{
  spint_t i;

  if (tproc->send_procs) z_free(tproc->send_procs);
  if (tproc->recv_procs) z_free(tproc->recv_procs);

  tproc->nsend_procs = -1;
  tproc->send_procs = NULL;
  tproc->nrecv_procs = -1;
  tproc->recv_procs = NULL;

  if (nsend_procs >= 0)
  {
    tproc->nsend_procs = nsend_procs;
    tproc->send_procs = z_alloc(nsend_procs, sizeof(sproc_t));

    for (i = 0; i < nsend_procs; ++i) tproc->send_procs[i] = send_procs[i];
  }

  if (nrecv_procs >= 0)
  {
    tproc->nrecv_procs = nrecv_procs;
    tproc->recv_procs = z_alloc(nrecv_procs, sizeof(sproc_t));

    for (i = 0; i < nrecv_procs; ++i) tproc->recv_procs[i] = recv_procs[i];

  } else if (nsend_procs >= 0)
  {
    if (comm == MPI_COMM_NULL) return 1;

    spec_make_recv_proclist(nsend_procs, send_procs, &tproc->nrecv_procs, &tproc->recv_procs, size, rank, comm);
  }

/*  printf("%d: send_procs (%" spint_fmt ") = ", rank, tproc->nsend_procs);
  for (i = 0; i < tproc->nsend_procs; ++i) printf("  %" sproc_fmt, tproc->send_procs[i]);
  printf("\n");

  printf("%d: recv_procs (%" spint_fmt ") = ", rank, tproc->nrecv_procs);
  for (i = 0; i < tproc->nrecv_procs; ++i) printf("  %" sproc_fmt, tproc->recv_procs[i]);
  printf("\n");*/

  return SPEC_EXIT_SUCCESS;
}

#endif


spint_t spec_print(spec_tproc_t tproc, spec_tproc_data_t tproc_data, spec_elem_t *b) /* sp_func spec_print */
{
  spint_t i, j, n;
  spec_proc_t p;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;


  /* setup tproc buffers */
  spec_tproc_setup(tproc, b, &procs, &mods);

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);


  if (tproc->tproc)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      p = tproc->tproc(spec_elem_get_buf(b), i, tproc_data);
      printf("%" spint_fmt ": %" spec_proc_fmt "\n", i, p);
    }

  } else if (tproc->tproc_mod)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      p = tproc->tproc_mod(spec_elem_get_buf(b), i, tproc_data, spec_elem_get_buf(mods));
      printf("%" spint_fmt ": %" spec_proc_fmt "\n", i, p);
    }

  } else if (tproc->tprocs)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      tproc->tprocs(spec_elem_get_buf(b), i, tproc_data, &n, procs);
      printf("%" spint_fmt ":", i);
      for (j = 0; j < n; ++j) printf(" %" spec_proc_fmt, procs[j]);
      printf("\n");
    }

  } else if (tproc->tprocs_mod)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      tproc->tprocs_mod(spec_elem_get_buf(b), i, tproc_data, &n, procs, spec_elem_get_buf(mods));
      printf("%" spint_fmt ":", i);
      for (j = 0; j < n; ++j) printf(" %" spec_proc_fmt, procs[j]);
      printf("\n");
    }
  }

  /* free tproc buffers */
  spec_tproc_release(&procs, &mods);

  return 0;
}



#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "spec_core.h"
#include "spec_common.h"


#ifndef SPEC_SENDRECV_TRACE_IF
# define SPEC_SENDRECV_TRACE_IF  (rank == -1)
#endif


/*#define SPEC_PRINT*/


/* sp_var spec_sendrecv_aux spec_sendrecv_aux_size spec_sendrecv_send_requests spec_sendrecv_receive_requests */
void *spec_sendrecv_aux = NULL;
spint_t spec_sendrecv_aux_size = -16*1024;
spint_t spec_sendrecv_send_requests = 10;
spint_t spec_sendrecv_receive_requests = 10;


#define SPEC_SENDRECV_SEND_ALL_AT_ONCE  0


spint_t spec_sendrecv_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_sendrecv_db */
{
  spint_t exit_code = SPEC_EXIT_SUCCESS;

  spint_t i, j, nprocs;
  spec_proc_t p;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;

  int *scounts, rcounts[3], ii;

  spint_t scount, sdispl, rdispl, rdispl_old, rdisplpart;
  spint_t stotal, rtotal, sdone, rdone, nfullrecvs, npartrecvs;

  const spint_t sreqs_max = spec_sendrecv_send_requests;
  const spint_t rreqs_max = spec_sendrecv_receive_requests;

/*  printf("%d: requests: %d / %d\n", rank, (int) sreqs_max, (int) rreqs_max);*/

  MPI_Request reqs[(rreqs_max + sreqs_max) * spec_elem_x];
  MPI_Status stats[(rreqs_max + sreqs_max) * spec_elem_x];
  int inds[(rreqs_max + sreqs_max) * spec_elem_x], completed_nreqs;
  spint_t completed_nsends, completed_nrecvs;
  spint_t reqs_i[(rreqs_max + sreqs_max) * 2], reqs_tmp;
  spint_t rreqs_nfree, rreqs_free[rreqs_max], sreqs_nfree, sreqs_free[sreqs_max];
  const spint_t rreqs_base = 0;
  const spint_t sreqs_base = rreqs_max;

  const spint_t aux_size_min = 1;
  const spint_t aux_size_max = -1;
  spint_t aux_size, aux_n, *aux_bases, *aux_displs, *aux_queue, aux_queue_size, aux_queue_first, aux_queue_next, aux_done, aux_tmp, aux_alloc;
  spec_elem_t aux;

  SPEC_DECLARE_TPROC_SENDRECV_DB

#define _MAX_ROUNDS 10

#ifdef MAX_ROUNDS
  spint_t max_rounds = MAX_ROUNDS;
#endif

#ifdef Z_PACK_TIMING
  double tt, t[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
#endif


  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  /* check supported tproc functions */
  if (!spec_check_tproc_support(tproc, 1, 1, 1, 1, rank, "spec_sendrecv_buffer_db"))
  {
    spec_elem_set_n(rb, 0);
    exit_code = SPEC_EXIT_FAILED;
    goto exit;
  }

  scounts = z_alloc(size * 3, sizeof(int));

  /* setup tproc buffers */
  spec_tproc_setup(tproc, sb, &procs, &mods);

#ifdef SPEC_PRINT
  printf("input\n");
  spec_print(&tproc, tproc_data, sb);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  spec_make_counts(tproc, tproc_data, sb, 0, size, scounts, procs);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

  /* make total, full, and partial receives */

#ifdef SPEC_PROCLISTS
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    aux_n = tproc->nsend_procs;

  } else
#endif
  {
    aux_n = size - 1;
  }

  if (spec_sendrecv_aux_size > 0) aux_size = spec_elem_sizefor(sb, spec_sendrecv_aux_size) / z_max(1, aux_n);
  else aux_size = spec_elem_sizefor(sb, -spec_sendrecv_aux_size);

  if (aux_size_min > 0) aux_size = z_max(aux_size, aux_size_min);
  if (aux_size_max > 0) aux_size = z_min(aux_size, aux_size_max);

  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "aux_n: %" spint_fmt ", aux_size: %" spint_fmt, aux_n, aux_size);

  scount = spec_elem_get_n(sb);

  rdisplpart = scounts[rank];

  stotal = 0;
  for (i = size - 1; i >= 0; --i)
  {
    stotal += scounts[i];
    scounts[3 * i + 0] = scounts[i];
    scounts[3 * i + 1] = scounts[i] / aux_size;
    scounts[3 * i + 2] = (scounts[i] % aux_size > 0)?1:0;

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "scount[%" spint_fmt "]: %d,%d,%d", i, scounts[3 * i + 0], scounts[3 * i + 1], scounts[3 * i + 2]);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[2]);

  spec_reduce_scatter_counts(scounts, rcounts, 3,
#ifdef SPEC_PROCLISTS
    tproc->nsend_procs, tproc->send_procs, tproc->nrecv_procs, tproc->recv_procs,
#endif
    size, rank, comm);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[2]);

  rtotal = rcounts[0];
  nfullrecvs = rcounts[1] - scounts[3 * rank + 1];
  npartrecvs = rcounts[2] - scounts[3 * rank + 2];

  rdisplpart += nfullrecvs * aux_size;

  /* check size of receive buffer */
  if (!spec_check_buffer_size(rb, rtotal, 1, comm, rank, "spec_sendrecv_buffer_db", "receive"))
  {
    spec_elem_set_n(rb, rtotal);
    exit_code = SPEC_EXIT_FAILED;
    goto free_and_exit;
  }

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "rtotal = %" spint_fmt ", nfullrecvs = %" spint_fmt ", npartrecvs = %" spint_fmt, rtotal, nfullrecvs, npartrecvs);

  /* redistribute */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[3]);

#define TAG_FULL  0
#define TAG_PART  spec_elem_x

#define AUX_BASE(_p_)         aux_bases[_p_]
#define AUX_DISPL(_p_)        aux_displs[_p_]
#define AUX_DISPL_BEGIN(_p_)  AUX_BASE(_p_)
#define AUX_DISPL_END(_p_)    AUX_BASE((_p_) + 1)
#define AUX_DISPL_SET(_p_, _d_)  aux_displs[_p_] = (_d_)
#define AUX_DISPL_INC(_p_)    ++aux_displs[_p_]
#define AUX_SIZE(_p_)         AUX_DISPL(_p_) - AUX_DISPL_BEGIN(_p_)
#define AUX_ENQUEUE(_p_)      Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "enqueue aux send of %" spint_fmt " to process %" spec_proc_fmt, AUX_SIZE(_p_), (_p_)); \
  aux_queue[aux_queue_next] = (_p_); ++aux_queue_next; aux_queue_next %= aux_queue_size;)
#define AUX_DEQUEUE()         (aux_tmp = aux_queue[aux_queue_first], ++aux_queue_first, aux_queue_first %= aux_queue_size, aux_tmp)
#define AUX_QUEUED()          ((aux_queue_next + aux_queue_size - aux_queue_first) % aux_queue_size)

#define REQS                          reqs
#define REQS_N                        (rreqs_max + sreqs_max) * spec_elem_x
#define REQS_PTR(_r_, _x_)            &reqs[(_r_) * spec_elem_x + (_x_)]

#define REQS_COMPLETE(_r_)            reqs_i[2 * (_r_) + 1]
#define REQS_COMPLETE_RESET(_r_)      reqs_i[2 * (_r_) + 1] = spec_elem_x
#define REQS_COMPLETE_ONE(_r_)        --reqs_i[2 * (_r_) + 1]
#define REQS_COMPLETED(_r_)           (REQS_COMPLETE(_r_) == 0)
#define REQS_COMPLETED_FIRST(_r_)     (REQS_COMPLETE(_r_) == spec_elem_x - 1)

#define REQS_NCUR()                   (sreqs_max - sreqs_nfree + rreqs_max - rreqs_nfree)

#define REQS_SEND_DST_SET(_r_, _p_)   reqs_i[2 * (_r_) + 0] = (_p_)
#define REQS_SEND_DST(_r_)            reqs_i[2 * (_r_) + 0]
#define REQS_SEND_SIZE(_r_)           AUX_SIZE(REQS_SEND_DST(_r_))
#define REQS_SEND_NFREE()             sreqs_nfree
#define REQS_SEND_FREE(_r_)           Z_MOP(sreqs_free[sreqs_nfree] = (_r_); ++sreqs_nfree; reqs_i[2 * (_r_) + 0] = reqs_i[2 * (_r_) + 1] = -2501;)

#define REQS_RECV(_r_)                ((_r_) < rreqs_max)
#define REQS_RECV_SIZE(_r_)           reqs_i[2 * (_r_) + 0]
#define REQS_RECV_SIZE_SET(_r_, _s_)  reqs_i[2 * (_r_) + 0] = (_s_)
#define REQS_RECV_NFREE()             rreqs_nfree
#define REQS_RECV_FREE(_r_)           Z_MOP(rreqs_free[rreqs_nfree] = (_r_); ++rreqs_nfree; reqs_i[2 * (_r_) + 0] = reqs_i[2 * (_r_) + 1] = -2501;)
#define REQS_RECV_FULL_FIRST(_r_)     (REQS_COMPLETE(_r_) == 1)

#define SEND_AUX_FIRST(_p_)  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "send_aux_first: process: %" spec_proc_fmt ", size: %" spint_fmt ", displ: %" spint_fmt, (_p_), AUX_SIZE(_p_), AUX_BASE(_p_)); \
  --sreqs_nfree; reqs_tmp = sreqs_free[sreqs_nfree]; \
  spec_elem_isend_first(&aux, AUX_BASE(_p_), AUX_SIZE(_p_), _p_, (AUX_SIZE(_p_) == aux_size)?TAG_FULL:TAG_PART, REQS_PTR(reqs_tmp, 0), size, rank, comm); \
  REQS_COMPLETE_RESET(reqs_tmp); \
  REQS_SEND_DST_SET(reqs_tmp, _p_); \
)

#define SEND_AUX_NEXT(_p_, _r_)  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "send_aux_next: process: %" spec_proc_fmt ", req: %" spint_fmt ", size: %" spint_fmt ", displ: %" spint_fmt, (_p_), (_r_), AUX_SIZE(_p_), AUX_BASE(_p_)); \
  spec_elem_isend_next(&aux, AUX_BASE(_p_), AUX_SIZE(_p_), _p_, (AUX_SIZE(_p_) == aux_size)?TAG_FULL:TAG_PART, REQS_PTR(_r_, 0), size, rank, comm); \
)

#define RECV_FULL_FIRST()  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv_full_first: displ: %" spint_fmt, rdispl); \
  --rreqs_nfree; reqs_tmp = rreqs_free[rreqs_nfree]; \
  spec_elem_irecv_first(rb, rdispl, aux_size, MPI_ANY_SOURCE, TAG_FULL, REQS_PTR(reqs_tmp, 0), size, rank, comm); \
  REQS_COMPLETE_RESET(reqs_tmp); \
  REQS_RECV_SIZE_SET(reqs_tmp, rdispl); \
  rdispl += aux_size; \
)

#define RECV_FULL_NEXT(_p_, _n_, _r_)  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv_full_next: process: %" spec_proc_fmt ", req: %" spint_fmt ", size: %" spint_fmt ", displ: %" spint_fmt, (_p_), (_r_), (_n_), REQS_RECV_SIZE(_r_)); \
  spec_elem_irecv_next(rb, REQS_RECV_SIZE(_r_), _n_, _p_, TAG_FULL, REQS_PTR(_r_, 0), size, rank, comm); \
  REQS_RECV_SIZE_SET(_r_, _n_); \
)

#define RECV_PART_FIRST()  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv_part_first: displ: %" spint_fmt, rdisplpart); \
  --rreqs_nfree; reqs_tmp = rreqs_free[rreqs_nfree]; \
  spec_elem_irecv_first(rb, rdisplpart, aux_size, MPI_ANY_SOURCE, TAG_PART, REQS_PTR(reqs_tmp, 0), size, rank, comm); \
  REQS_COMPLETE_RESET(reqs_tmp); \
)

#define RECV_PART_NEXT(_p_, _n_, _r_)  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv_part_next: process: %" spec_proc_fmt ", req: %" spint_fmt ", size: %" spint_fmt ", displ: %" spint_fmt, (_p_), (_r_), (_n_), rdisplpart); \
  spec_elem_irecv_next(rb, rdisplpart, _n_, _p_, TAG_PART, REQS_PTR(_r_, 0), size, rank, comm); \
  REQS_RECV_SIZE_SET(_r_, _n_); \
  rdisplpart += _n_; \
)

  p = SPEC_PROC_NONE;
  nprocs = -1;

  sdone = rdone = 0;
  sdispl = rdispl = 0;

  rreqs_nfree = rreqs_max;
  for (i = 0; i < rreqs_max; ++i)
  {
    rreqs_free[i] = rreqs_base + i;
    for (j = 0; j < spec_elem_x; ++j) reqs[rreqs_free[i] * spec_elem_x + j] = MPI_REQUEST_NULL;
    reqs_i[2 * rreqs_free[i] + 0] = reqs_i[2 * rreqs_free[i] + 1] = -2501;
  }
  sreqs_nfree = sreqs_max;
  for (i = 0; i < sreqs_max; ++i)
  {
    sreqs_free[i] = sreqs_base + i;
    for (j = 0; j < spec_elem_x; ++j) reqs[sreqs_free[i] * spec_elem_x + j] = MPI_REQUEST_NULL;
    reqs_i[2 * sreqs_free[i] + 0] = reqs_i[2 * sreqs_free[i] + 1] = -2501;
  }

  aux_queue_size = aux_n + 1;  /* aux queue is a circular buffer implemented with only two indices 'next' (to write) and 'first' (to read), thus the queue has to be one slot larger than max. number of required slots to prevent a full queue */
  aux_bases = z_alloc(2 * size + 1 + aux_queue_size, sizeof(spint_t));
  aux_displs = aux_bases + size + 1;

  for (i = 0; i < size; ++i) aux_bases[i] = 0;
#ifdef SPEC_PROCLISTS
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    for (i = 0; i < tproc->nsend_procs; ++i) aux_bases[tproc->send_procs[i]] = aux_size;

  } else
#endif
  {
    for (i = 0; i < size; ++i) aux_bases[i] = aux_size;
  }
  aux_bases[rank] = 0;
  j = 0;
  for (i = 0; i < size; ++i)
  {
    aux_displs[i] = j;
    j += aux_bases[i];
    aux_bases[i] = aux_displs[i];
  }
  aux_bases[size] = j;

  aux_queue = aux_displs + size;
  aux_queue_first = aux_queue_next = 0;
  aux_done = 0;

  spec_elem_unset(&aux);

  spec_elem_copy_type(sb, &aux);

  spec_elem_set_nmax(&aux, 0);
#ifdef spec_elem_alloc_tmp_from_block
  if (spec_sendrecv_aux)
  {
    i = (spec_sendrecv_aux_size > 0)?spec_sendrecv_aux_size:(-spec_sendrecv_aux_size * aux_n);
    spec_elem_alloc_tmp_from_block(&aux, spec_sendrecv_aux, i);
  }
#endif

  aux_alloc = 0;
  if (spec_elem_get_nmax(&aux) < aux_size * aux_n)
  {
    aux_alloc = 1;
    spec_elem_alloc_tmp(&aux, aux_n * aux_size);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[4]);

  Z_TIMING_DECL(const int ottmps = 5;);

  while (1)
  {
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "send: %" spint_fmt " < %" spint_fmt ", recv: %" spint_fmt " < %" spint_fmt "", sdone, stotal, rdone, rtotal);

    if (!(sdone < stotal || rdone < rtotal)) break;

#ifdef MAX_ROUNDS
    if (max_rounds-- <= 0) break;
#endif

    /* if there are no more data elements, then skip processing loop */
    if (sdispl >= scount) goto do_send;

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "loop over data elements");

    Z_TIMING_START(tt);

    rdispl_old = rdispl;

    if (tproc->tproc)
    {
#if 0
      if (tproc->tproc_ext.sendrecv_db) p = tproc->tproc_ext.sendrecv_db(tproc_data, sb, rb, scount, &sdispl, &rdispl, &aux, aux_displs, aux_size, aux_queue, &aux_queue_next, aux_queue_size, rank, p);
      else SPEC_DO_TPROC_SENDRECV_DB(tproc->tproc, tproc_data, sb, rb, scount, &sdispl, &rdispl, &aux, aux_displs, aux_size, aux_queue, &aux_queue_next, aux_queue_size, rank, p);
#else
      while (sdispl < scount)
      {
        if (p == SPEC_PROC_NONE) p = tproc->tproc(spec_elem_get_buf(sb), sdispl, tproc_data);

        if (p != SPEC_PROC_NONE)
        {
          if (p == rank)
          {
            spec_elem_copy_at(sb, sdispl, rb, rdispl);
            ++rdispl;

          } else
          {
            if (AUX_DISPL(p) >= AUX_DISPL_END(p)) break;

            spec_elem_copy_at(sb, sdispl, &aux, AUX_DISPL(p));

            AUX_DISPL_INC(p);

            if (AUX_DISPL(p) >= AUX_DISPL_END(p)) AUX_ENQUEUE(p);
          }
        }

        p = SPEC_PROC_NONE;
        ++sdispl;
      }
#endif

    } else if (tproc->tproc_mod)
    {
      while (sdispl < scount)
      {
        if (p == SPEC_PROC_NONE) p = tproc->tproc_mod(spec_elem_get_buf(sb), sdispl, tproc_data, spec_elem_get_buf(mods));

        if (p != SPEC_PROC_NONE)
        {
          if (p == rank)
          {
            spec_elem_copy_at(mods, 0, rb, rdispl);
            ++rdispl;

          } else
          {
            if (AUX_DISPL(p) >= AUX_DISPL_END(p)) break;

            spec_elem_copy_at(mods, 0, &aux, AUX_DISPL(p));

            AUX_DISPL_INC(p);

            if (AUX_DISPL(p) >= AUX_DISPL_END(p)) AUX_ENQUEUE(p);
          }
        }

        p = SPEC_PROC_NONE;
        ++sdispl;
      }

    } else if (tproc->tprocs)
    {
      while (sdispl < scount)
      {
        if (nprocs < 0) { tproc->tprocs(spec_elem_get_buf(sb), sdispl, tproc_data, &nprocs, procs); --nprocs; }

        while (nprocs >= 0)
        {
          if (procs[nprocs] == rank)
          {
            spec_elem_copy_at(sb, sdispl, rb, rdispl);
            ++rdispl;

          } else
          {
            if (AUX_DISPL(procs[nprocs]) >= AUX_DISPL_END(procs[nprocs])) break;

            spec_elem_copy_at(sb, sdispl, &aux, AUX_DISPL(procs[nprocs]));

            AUX_DISPL_INC(procs[nprocs]);

            if (AUX_DISPL(procs[nprocs]) >= AUX_DISPL_END(procs[nprocs])) AUX_ENQUEUE(procs[nprocs]);
          }

          --nprocs;
        }

        if (nprocs >= 0) break;

        ++sdispl;
      }

    } else if (tproc->tprocs_mod)
    {
      while (sdispl < scount)
      {
        if (nprocs < 0) { tproc->tprocs_mod(spec_elem_get_buf(sb), sdispl, tproc_data, &nprocs, procs, spec_elem_get_buf(mods)); --nprocs; }

        while (nprocs >= 0)
        {
          if (procs[nprocs] == rank)
          {
            spec_elem_copy_at(mods, nprocs, rb, rdispl);
            ++rdispl;

          } else
          {
            if (AUX_DISPL(procs[nprocs]) >= AUX_DISPL_END(procs[nprocs])) break;

            spec_elem_copy_at(mods, nprocs, &aux, AUX_DISPL(procs[nprocs]));

            AUX_DISPL_INC(procs[nprocs]);

            if (AUX_DISPL(procs[nprocs]) >= AUX_DISPL_END(procs[nprocs])) AUX_ENQUEUE(procs[nprocs]);
          }

          --nprocs;
        }

        if (nprocs >= 0) break;

        ++sdispl;
      }
    }

    sdone += rdispl - rdispl_old;
    rdone += rdispl - rdispl_old;

    Z_TIMING_STOP_ADD(tt, t[ottmps + 0]);

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "scount: %" spint_fmt ", sdispl: %" spint_fmt ", rdispl: %" spint_fmt ", rdisplpart: %" spint_fmt, scount, sdispl, rdispl, rdisplpart);
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "stotal: %" spint_fmt ", sdone: %" spint_fmt ", rtotal: %" spint_fmt ", rdone: %" spint_fmt, stotal, sdone, rtotal, rdone);
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "aux_queue: first: %" spint_fmt ", next: %" spint_fmt, aux_queue_first, aux_queue_next);

    Z_TIMING_START(tt);

    if (!aux_done && sdispl >= scount)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "no more local data elements, enqueueing all remaining aux sends");

      for (i = 0; i < size; ++i)
      {
        j = (rank + i) % size;
        if (AUX_SIZE(j) > 0 && AUX_SIZE(j) < aux_size) AUX_ENQUEUE((spec_proc_t) j);
      }

      aux_done = 1;

      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "aux_queue: first: %" spint_fmt ", next: %" spint_fmt, aux_queue_first, aux_queue_next);
    }

    Z_TIMING_STOP_ADD(tt, t[ottmps + 1]);

do_send:
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "A: sreqs_nfree: %" spint_fmt ", aux_queued: %" spint_fmt "", REQS_SEND_NFREE(), AUX_QUEUED());

    Z_TIMING_START(tt);

    /* initiate queued sends */
    while (REQS_SEND_NFREE() > 0 && AUX_QUEUED() > 0)
    {
      i = AUX_DEQUEUE();
      SEND_AUX_FIRST((spec_proc_t) i);
#if SPEC_SENDRECV_SEND_ALL_AT_ONCE
      SEND_AUX_NEXT((spec_proc_t) i, reqs_tmp);
#endif
    }

    Z_TIMING_STOP_ADD(tt, t[ottmps + 2]);

do_recv:
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "B: sreqs_nfree: %" spint_fmt ", rreqs_nfree: %" spint_fmt "", REQS_SEND_NFREE(), REQS_RECV_NFREE());

    Z_TIMING_START(tt);

    /* initiate partial recv */
    if (REQS_RECV_NFREE() > 0 && npartrecvs > 0)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "initiate part first recv");
      RECV_PART_FIRST();
      --npartrecvs;
      npartrecvs *= -1; /* negative value signals that a part recv is active */
    }

    /* initiate remaining full recvs */
    while (REQS_RECV_NFREE() > 0 && nfullrecvs > 0)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "initiate full first recv");
      RECV_FULL_FIRST();
      --nfullrecvs;
    }

    Z_TIMING_STOP_ADD(tt, t[ottmps + 3]);

/*    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "reqs:");
    for (i = 0; i < rreqs_max + sreqs_max; ++i) Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "  (%" spint_fmt " / %" spint_fmt ")", reqs_i[2 * i + 0], reqs_i[2 * i + 1]);*/

do_check_requests:
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "C: sreqs_nfree: %" spint_fmt ", rreqs_nfree: %" spint_fmt "", REQS_SEND_NFREE(), REQS_RECV_NFREE());

    if (REQS_NCUR() == 0)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "no requests, continue loop");
      continue;
    }

    Z_TIMING_START(tt);

#if 0
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "waitany:");
    MPI_Waitany(REQS_N, REQS, &inds[0], &stats[0]);
    completed_nreqs = (inds[0] == MPI_UNDEFINED)?0:1;
#else
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "waitsome:");
    MPI_Waitsome(REQS_N, REQS, &completed_nreqs, inds, stats);
#endif

    Z_TIMING_STOP_ADD(tt, t[ottmps + 4]);
    Z_TIMING_CMD(++t[ottmps + 5];);

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "completed requests: %d", completed_nreqs);

    Z_TIMING_START(tt);

    completed_nsends = completed_nrecvs = 0;

    for (j = 0; j < completed_nreqs; ++j)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "check request #%" spint_fmt ": index: %d", j, inds[j]);

      if (inds[j] == MPI_UNDEFINED) continue;

      i = inds[j] / spec_elem_x;

      REQS_COMPLETE_ONE(i);

      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "%s-request %" spint_fmt ": (%" spint_fmt " / %" spint_fmt ")", (REQS_RECV(i))?"recv":"send", i, reqs_i[2 * i + 0], reqs_i[2 * i + 1]);

      if (REQS_RECV(i)) /* recv-req done */
      {
        Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv-request %" spint_fmt ": tag: %d, first: %s, completed: %s", i, stats[j].MPI_TAG, REQS_COMPLETED_FIRST(i)?"yes":"no", REQS_COMPLETED(i)?"yes":"no");

        /* if first component of receive, then receive next */
        if (REQS_COMPLETED_FIRST(i))
        {
          spec_elem_get_recv_count(rb, &stats[j], &ii);
          if (stats[j].MPI_TAG == TAG_FULL) RECV_FULL_NEXT(stats[j].MPI_SOURCE, (spint_t) ii, i);
          else RECV_PART_NEXT(stats[j].MPI_SOURCE, (spint_t) ii, i);
        }

        /* nothing was complete, thus continue with next completed request */
        if (!REQS_COMPLETED(i)) continue;

        /* count completed receives */
        ++completed_nrecvs;

        /* register recv */
        rdone += REQS_RECV_SIZE(i);

        /* free recv-req */
        REQS_RECV_FREE(i);

        /* set part recv inactive */
        if (stats[j].MPI_TAG >= TAG_PART) npartrecvs *= -1;

      } else /* send-req done */
      {
        Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "send-request %" spint_fmt " first: %s, completed: %s", i, REQS_COMPLETED_FIRST(i)?"yes":"no", REQS_COMPLETED(i)?"yes":"no");

#if !(SPEC_SENDRECV_SEND_ALL_AT_ONCE)
        /* if first component of send, then send next */
        if (REQS_COMPLETED_FIRST(i)) SEND_AUX_NEXT((spec_proc_t) REQS_SEND_DST(i), i);
#endif

        /* nothing was complete, thus continue with next completed request */
        if (!REQS_COMPLETED(i)) continue;

        /* count completed sends */
        ++completed_nsends;

        /* register send */
        sdone += REQS_SEND_SIZE(i);

        /* reset aux */
        AUX_DISPL_SET(REQS_SEND_DST(i), AUX_DISPL_BEGIN(REQS_SEND_DST(i)));

        /* free send-req */
        REQS_SEND_FREE(i);

        Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "sreqs_nfree: %" spint_fmt "", sreqs_nfree);
      }
    }

    Z_TIMING_STOP_ADD(tt, t[ottmps + 6]);

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "completed sends: %" spint_fmt ", completed receives: %" spint_fmt "", completed_nsends, completed_nrecvs);

    /* if nothing was completed, then wait for next requests */
    if (completed_nsends == 0 && completed_nrecvs == 0) goto do_check_requests;

    /* if only receives were completed, then start new receives */
    if (completed_nsends == 0 && completed_nrecvs > 0) goto do_recv;
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  if (aux_alloc) spec_elem_free_tmp(&aux);

  z_free(aux_bases);

#if 1
  Z_TIMING_CMD(
    const int nttmps = (sizeof(t) / sizeof(double)) - ottmps;
    double ttmps[nttmps];
    MPI_Allreduce(&t[ottmps], ttmps, nttmps, MPI_DOUBLE, MPI_MAX, comm);
    for (i = 0; i < nttmps; ++i) t[ottmps + i] = ttmps[i];
  );
#endif

  spec_elem_set_n(rb, rtotal);

free_and_exit:

  /* free tproc buffers */
  spec_tproc_release(&procs, &mods);

  z_free(scounts);

exit:
  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[0]);

  Z_TIMING_PRINT(0, __func__, sizeof(t) / sizeof(t[0]), t, rank);

#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
  if (spec_timing)
    for (i = 0; i < sizeof(t) / sizeof(t[0]); ++i) spec_timing[i] = t[i];
#endif

  return exit_code;
}

#undef TAG_FULL
#undef TAG_PART

#undef AUX_BASE
#undef AUX_DISPL
#undef AUX_DISPL_BEGIN
#undef AUX_DISPL_END
#undef AUX_DISPL_SET
#undef AUX_DISPL_INC
#undef AUX_SIZE
#undef AUX_ENQUEUE
#undef AUX_DEQUEUE
#undef AUX_QUEUED

#undef REQS
#undef REQS_N
#undef REQS_PTR

#undef REQS_COMPLETE
#undef REQS_COMPLETE_RESET
#undef REQS_COMPLETE_ONE
#undef REQS_COMPLETED

#undef REQS_NCUR

#undef REQS_SEND_DST_SET
#undef REQS_SEND_DST
#undef REQS_SEND_SIZE
#undef REQS_SEND_NFREE
#undef REQS_SEND_FREE

#undef REQS_RECV
#undef REQS_RECV_SIZE
#undef REQS_RECV_SIZE_SET
#undef REQS_RECV_NFREE
#undef REQS_RECV_FREE
#undef REQS_RECV_FULL_FIRST

#undef SEND_AUX
#undef RECV_FULL_FIRST
#undef RECV_FULL_NEXT
#undef RECV_PART_FIRST
#undef RECV_PART_NEXT
