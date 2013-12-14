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

#ifdef HAVE_ZMPI_TOOLS_H
# include "zmpi_tools.h"
#endif

#include "spec_core.h"
#include "spec_common.h"


#ifndef SPEC_ALLTOALLVG_TRACE_IF
# define SPEC_ALLTOALLVG_TRACE_IF  (rank == -1)
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


  Z_TRACE_IF(SPEC_ALLTOALLVG_TRACE_IF, "spec_alltoallv_db");

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
  spec_print(&tproc, tproc_data, sb);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  spec_make_counts(sb, tproc, tproc_data, 0, size, scounts, procs);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

#ifdef SPEC_PRINT
  printf("after count\n");
  spec_print(&tproc, tproc_data, sb);
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
    spec_elem_alloc_buf(&_xb, stotal);

    xb = &_xb;
  }

  spec_elem_copy_type(sb, xb);

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* local rearrange */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[3]);

  if (tproc->tproc)
  {
    if (tproc->tproc_rearrange_db) tproc->tproc_rearrange_db(sb, xb, tproc_data, sdispls);
    else SPEC_DO_TPROC_REARRANGE_DB(tproc->tproc, tproc_data, sb, xb, sdispls);

  } else if (tproc->tproc_mod)
  {
    if (tproc->tproc_mod_rearrange_db) tproc->tproc_mod_rearrange_db(sb, xb, tproc_data, sdispls, mods);
    else SPEC_DO_TPROC_MOD_REARRANGE_DB(tproc->tproc_mod, tproc_data, sb, xb, sdispls, mods);

  } else if (tproc->tprocs)
  {
    if (tproc->tprocs_rearrange_db) tproc->tprocs_rearrange_db(sb, xb, tproc_data, sdispls, procs);
    else SPEC_DO_TPROCS_REARRANGE_DB(tproc->tprocs, tproc_data, sb, xb, sdispls, procs);

  } else if (tproc->tprocs_mod)
  {
    if (tproc->tprocs_mod_rearrange_db) tproc->tprocs_mod_rearrange_db(sb, xb, tproc_data, sdispls, procs, mods);
    else SPEC_DO_TPROCS_MOD_REARRANGE_DB(tproc->tprocs_mod, tproc_data, sb, xb, sdispls, procs, mods);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

#ifdef SPEC_PRINT
  printf("after rearrange\n");
  spec_print(&tproc, tproc_data, xb);
#endif

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

  if (xb == &_xb) spec_elem_free_buf(&_xb);

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


  Z_TRACE_IF(SPEC_ALLTOALLVG_TRACE_IF, "spec_alltoallv_ip");

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
  spec_print(&tproc, tproc_data, b);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  spec_make_counts(b, tproc, tproc_data, 1, size, scounts, procs);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

#ifdef SPEC_PRINT
  printf("after count\n");
  spec_print(&tproc, tproc_data, &b);
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
    spec_elem_alloc_buf(&_xb, 1);

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
      if (tproc->tproc_rearrange_db) tproc->tproc_rearrange_db(b, tb, tproc_data, sdispls);
      else SPEC_DO_TPROC_REARRANGE_DB(tproc->tproc, tproc_data, b, tb, sdispls);

    } else
    {
      if (tproc->tproc_rearrange_ip) tproc->tproc_rearrange_ip(b, xb, tproc_data, sdispls, scounts, size);
      else SPEC_DO_TPROC_REARRANGE_IP(tproc->tproc, tproc_data, b, xb, sdispls, scounts, size);
    }

  } else if (tproc->tproc_mod)
  {
    if (tb)
    {
      if (tproc->tproc_mod_rearrange_db) tproc->tproc_mod_rearrange_db(b, tb, tproc_data, sdispls, mods);
      else SPEC_DO_TPROC_MOD_REARRANGE_DB(tproc->tproc_mod, tproc_data, b, tb, sdispls, mods);

    } else
    {
      if (tproc->tproc_mod_rearrange_ip) tproc->tproc_mod_rearrange_ip(b, xb, tproc_data, sdispls, scounts, size, mods);
      else SPEC_DO_TPROC_MOD_REARRANGE_IP(tproc->tproc_mod, tproc_data, b, xb, sdispls, scounts, size, mods);
    }

  } else if (tproc->tprocs)
  {
    if (tb)
    {
      if (tproc->tprocs_rearrange_db) tproc->tprocs_rearrange_db(b, tb, tproc_data, sdispls, procs);
      else SPEC_DO_TPROCS_REARRANGE_DB(tproc->tprocs, tproc_data, b, tb, sdispls, procs);

    } else
    {
      if (tproc->tprocs_rearrange_ip) tproc->tprocs_rearrange_ip(b, xb, tproc_data, sdispls, scounts, size, procs);
      else SPEC_DO_TPROCS_REARRANGE_IP(tproc->tprocs, tproc_data, b, xb, sdispls, scounts, size, procs);
    }

  } else if (tproc->tprocs_mod)
  {
    if (tb)
    {
      if (tproc->tprocs_mod_rearrange_db) tproc->tprocs_mod_rearrange_db(b, tb, tproc_data, sdispls, procs, mods);
      else SPEC_DO_TPROCS_MOD_REARRANGE_DB(tproc->tprocs_mod, tproc_data, b, tb, sdispls, procs, mods);

    } else
    {
      if (tproc->tprocs_mod_rearrange_ip) tproc->tprocs_mod_rearrange_ip(b, xb, tproc_data, sdispls, scounts, size, procs, mods);
      else SPEC_DO_TPROCS_MOD_REARRANGE_IP(tproc->tprocs_mod, tproc_data, b, xb, sdispls, scounts, size, procs, mods);
    }
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

#ifdef SPEC_PRINT
  printf("after rearrange\n");
  spec_print(&tproc, tproc_data, &b);
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

  if (xb == &_xb) spec_elem_free_buf(&_xb);

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


#ifndef SPEC_COMMONG_TRACE_IF
# define SPEC_COMMONG_TRACE_IF  (rank == -1)
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


spint_t spec_make_counts(spec_elem_t *b, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int ip, int size, int *counts, spec_proc_t *procs) /* sp_func spec_make_counts */
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
      if (tproc->tproc_count_db) tproc->tproc_count_db(b, tproc_data, counts);
      else SPEC_DO_TPROC_COUNT_DB(tproc->tproc, tproc_data, b, counts);

    } else if (tproc->tproc_mod)
    {
      if (tproc->tproc_mod_count_db) tproc->tproc_mod_count_db(b, tproc_data, counts);
      else SPEC_DO_TPROC_MOD_COUNT_DB(tproc->tproc_mod, tproc_data, b, counts);

    } else if (tproc->tprocs)
    {
      if (tproc->tprocs_count_db) tproc->tprocs_count_db(b, tproc_data, counts, procs);
      else SPEC_DO_TPROCS_COUNT_DB(tproc->tprocs, tproc_data, b, counts, procs);

    } else if (tproc->tprocs_mod)
    {
      if (tproc->tprocs_mod_count_db) tproc->tprocs_mod_count_db(b, tproc_data, counts, procs);
      else SPEC_DO_TPROCS_MOD_COUNT_DB(tproc->tprocs_mod, tproc_data, b, counts, procs);
    }

  } else {

    if (tproc->tproc)
    {
      if (tproc->tproc_count_ip) tproc->tproc_count_ip(b, tproc_data, counts);
      else SPEC_DO_TPROC_COUNT_IP(tproc->tproc, tproc_data, b, counts);

    } else if (tproc->tproc_mod)
    {
      if (tproc->tproc_mod_count_ip) tproc->tproc_mod_count_ip(b, tproc_data, counts);
      else SPEC_DO_TPROC_MOD_COUNT_IP(tproc->tproc_mod, tproc_data, b, counts);

    } else if (tproc->tprocs)
    {
      if (tproc->tprocs_count_ip) tproc->tprocs_count_ip(b, tproc_data, counts, procs);
      else SPEC_DO_TPROCS_COUNT_IP(tproc->tprocs, tproc_data, b, counts, procs);

    } else if (tproc->tprocs_mod)
    {
      if (tproc->tprocs_mod_count_ip) tproc->tprocs_mod_count_ip(b, tproc_data, counts, procs);
      else SPEC_DO_TPROCS_MOD_COUNT_IP(tproc->tprocs_mod, tproc_data, b, counts, procs);
    }
  }

  return 0;
}


/* sp_var spec_redistribute_counts_type spec_redistribute_counts_proclists_type */
spint_t spec_redistribute_counts_type = SPEC_REDISTRIBUTE_COUNTS_DEFAULT;
spint_t spec_redistribute_counts_proclists_type = SPEC_REDISTRIBUTE_COUNTS_ALLTOALL;


spint_t spec_redistribute_counts(int *scounts, int *rcounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm) /* sp_func spec_redistribute_counts */
{
  Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "spec_redistribute_counts");

#ifdef HAVE_ZMPI_ALLTOALL_INT

#ifdef SPEC_PROCLISTS
  if (nsend_procs >= 0 || nrecv_procs >= 0)
  {
    memset(rcounts, 0, size * sizeof(int));

    switch (spec_redistribute_counts_proclists_type)
    {
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ISENDIRECV:
        Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ISENDIRECV");
        ZMPI_Alltoall_int_c2c_proclists_isendirecv(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ALLTOALLV:
        Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ALLTOALLV");
        ZMPI_Alltoall_int_c2c_proclists_alltoallv(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT:
        Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT");
        ZMPI_Alltoall_int_c2c_proclists_put(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_ALLOC:
        Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_ALLOC");
        ZMPI_Alltoall_int_c2c_proclists_put_alloc(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES:
        Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES");
        ZMPI_Alltoall_int_c2c_proclists_put_2phases(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES_ALLOC:
        Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES_ALLOC");
        ZMPI_Alltoall_int_c2c_proclists_put_2phases_alloc(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
    }
  } else
#endif
  {
#ifdef SPEC_MPI_ALLTOALL_2STEP_THRESHOLD
    if (size >= SPEC_MPI_ALLTOALL_2STEP_THRESHOLD)
    {
      Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_MPI_ALLTOALL_2STEP_THRESHOLD");
      ZMPI_Alltoall_int_c2c_2step(scounts, rcounts, comm);

    } else
#endif
    {
      switch (spec_redistribute_counts_type)
      {
        case SPEC_REDISTRIBUTE_COUNTS_ALLTOALL:
          Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_ALLTOALL");
          ZMPI_Alltoall_int_c2c_alltoall(scounts, rcounts, comm);
          break;
        case SPEC_REDISTRIBUTE_COUNTS_2STEP:
          Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_2STEP");
          ZMPI_Alltoall_int_c2c_2step(scounts, rcounts, comm);
          break;
        case SPEC_REDISTRIBUTE_COUNTS_PUT:
          Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT");
          ZMPI_Alltoall_int_c2c_put(scounts, rcounts, comm);
          break;
        case SPEC_REDISTRIBUTE_COUNTS_PUT_ALLOC:
          Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_ALLOC");
          ZMPI_Alltoall_int_c2c_put_alloc(scounts, rcounts, comm);
          break;
        case SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES:
          Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES");
          ZMPI_Alltoall_int_c2c_put_2phases(scounts, rcounts, comm);
          break;
        case SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES_ALLOC:
          Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES_ALLOC");
          ZMPI_Alltoall_int_c2c_put_2phases_alloc(scounts, rcounts, comm);
          break;
        case SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES:
          Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES");
          ZMPI_Alltoall_int_c2c_put_3phases(scounts, rcounts, comm);
          break;
        case SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES_ALLOC:
          Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES_ALLOC");
          ZMPI_Alltoall_int_c2c_put_3phases_alloc(scounts, rcounts, comm);
          break;
      }
    }
  }

#else /* HAVE_ZMPI_ALLTOALL_INT */

  Z_TRACE_IF(SPEC_COMMONG_TRACE_IF, "MPI_Alltoall");
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);

#endif /* HAVE_ZMPI_ALLTOALL_INT */

  return 0;
}


/* sp_var spec_redistribute_displs_type */
spint_t spec_redistribute_displs_type = SPEC_REDISTRIBUTE_DISPLS_DEFAULT;


spint_t spec_redistribute_displs(int *sdispls, int stotal, int *rdispls,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm) /* sp_func spec_redistribute_displs */
{
  int i;

  MPI_Win win;


  if (spec_redistribute_displs_type == SPEC_REDISTRIBUTE_DISPLS_COUNTS)
  {
    spec_redistribute_counts(sdispls, rdispls,
#ifdef SPEC_PROCLISTS
      nsend_procs, send_procs, nrecv_procs, recv_procs,
#endif
      size, rank, comm);

  } else if (spec_redistribute_displs_type == SPEC_REDISTRIBUTE_DISPLS_PUT)
  {
    MPI_Win_create(rdispls, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &win);
    MPI_Win_fence(0, win);

#ifdef SPEC_PROCLISTS
    if (nsend_procs >= 0 || nrecv_procs >= 0)
    {
      for (i = 0; i < nsend_procs; ++i)
        if (sdispls[send_procs[i]] != ((send_procs[i] + 1 < size)?sdispls[send_procs[i] + 1]:stotal))
          MPI_Put(&sdispls[send_procs[i]], 1, MPI_INT, send_procs[i], rank, 1, MPI_INT, win);

    } else
#endif
    {
      for (i = 0; i < size - 1; ++i)
        if (sdispls[i] != sdispls[i + 1]) MPI_Put(&sdispls[i], 1, MPI_INT, i, rank, 1, MPI_INT, win);
      if (sdispls[size - 1] != stotal) MPI_Put(&sdispls[size - 1], 1, MPI_INT, i, rank, 1, MPI_INT, win);
    }

    MPI_Win_fence(0, win);
    MPI_Win_free(&win);
  }

  return 0;
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


static void spec_tproc_unset_tproc(spec_tproc_t *tproc)
{
  (*tproc)->tproc = NULL;
}


static void spec_tproc_unset_ext_tproc(spec_tproc_t *tproc)
{
  (*tproc)->tproc_count_db = NULL;
  (*tproc)->tproc_count_ip = NULL;
  (*tproc)->tproc_rearrange_db = NULL;
  (*tproc)->tproc_rearrange_ip = NULL;
  (*tproc)->tproc_indices_db = NULL;
}


static void spec_tproc_unset_tproc_mod(spec_tproc_t *tproc)
{
  (*tproc)->tproc_mod = NULL;
}


static void spec_tproc_unset_ext_tproc_mod(spec_tproc_t *tproc)
{
  (*tproc)->tproc_mod_count_db = NULL;
  (*tproc)->tproc_mod_count_ip = NULL;
  (*tproc)->tproc_mod_rearrange_db = NULL;
  (*tproc)->tproc_mod_rearrange_ip = NULL;
}


static void spec_tproc_unset_tprocs(spec_tproc_t *tproc)
{
  (*tproc)->tprocs = NULL;
}


static void spec_tproc_unset_ext_tprocs(spec_tproc_t *tproc)
{
  (*tproc)->tprocs_count_db = NULL;
  (*tproc)->tprocs_count_ip = NULL;
  (*tproc)->tprocs_rearrange_db = NULL;
  (*tproc)->tprocs_rearrange_ip = NULL;
  (*tproc)->tprocs_indices_db = NULL;
}


static void spec_tproc_unset_tprocs_mod(spec_tproc_t *tproc)
{
  (*tproc)->tprocs_mod = NULL;
}


static void spec_tproc_unset_ext_tprocs_mod(spec_tproc_t *tproc)
{
  (*tproc)->tprocs_mod_count_db = NULL;
  (*tproc)->tprocs_mod_count_ip = NULL;
  (*tproc)->tprocs_mod_rearrange_db = NULL;
  (*tproc)->tprocs_mod_rearrange_ip = NULL;
}


spint_t spec_tproc_create(spec_tproc_t *tproc, spec_tproc_f *func, spec_tproc_mod_f *func_mod, spec_tprocs_f *func_s, spec_tprocs_mod_f *func_s_mod, spint_t max_tprocs) /* sp_func spec_tproc_create */
{
  *tproc = z_alloc(1, sizeof(struct _spec_tproc_t));

  (*tproc)->max_tprocs = max_tprocs;

  (*tproc)->tproc = func;
  spec_tproc_unset_ext_tproc(tproc);

  (*tproc)->tproc_mod = func_mod;
  spec_tproc_unset_ext_tproc_mod(tproc);

  (*tproc)->tprocs = func_s;
  spec_tproc_unset_ext_tprocs(tproc);

  (*tproc)->tprocs_mod = func_s_mod;
  spec_tproc_unset_ext_tprocs_mod(tproc);

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
  spec_tproc_set_proclists(tproc, -1, NULL, -1, NULL, 0, -1, MPI_COMM_NULL);
#endif

  z_free(*tproc);

  *tproc = NULL;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_duplicate(spec_tproc_t *tproc, spec_tproc_t *newtproc) /* sp_func spec_tproc_duplicate */
{
  *newtproc = z_alloc(1, sizeof(struct _spec_tproc_t));

  memcpy(*newtproc, *tproc, sizeof(struct _spec_tproc_t));

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tproc(spec_tproc_t *tproc, spec_tproc_f *func) /* sp_func spec_tproc_set_tproc */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  (*tproc)->tproc = func;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tproc(spec_tproc_t *tproc, spec_tproc_count_f *func_count_db, spec_tproc_count_f *func_count_ip, spec_tproc_rearrange_db_f *func_rearrange_db, spec_tproc_rearrange_ip_f *func_rearrange_ip, spec_tproc_indices_db_f *func_indices_db) /* sp_func spec_tproc_set_ext_tproc */
{
  (*tproc)->tproc_count_db = func_count_db;
  (*tproc)->tproc_count_ip = func_count_ip;
  (*tproc)->tproc_rearrange_db = func_rearrange_db;
  (*tproc)->tproc_rearrange_ip = func_rearrange_ip;
  (*tproc)->tproc_indices_db = func_indices_db;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tproc_mod(spec_tproc_t *tproc, spec_tproc_mod_f *func_mod) /* sp_func spec_tproc_set_tproc_mod */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  (*tproc)->tproc_mod = func_mod;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tproc_mod(spec_tproc_t *tproc, spec_tproc_mod_count_f *func_count_db, spec_tproc_mod_count_f *func_count_ip, spec_tproc_mod_rearrange_db_f *func_rearrange_db, spec_tproc_mod_rearrange_ip_f *func_rearrange_ip) /* sp_func spec_tproc_set_ext_tproc_mod */
{
  (*tproc)->tproc_mod_count_db = func_count_db;
  (*tproc)->tproc_mod_count_ip = func_count_ip;
  (*tproc)->tproc_mod_rearrange_db = func_rearrange_db;
  (*tproc)->tproc_mod_rearrange_ip = func_rearrange_ip;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tprocs(spec_tproc_t *tproc, spec_tprocs_f *func_s, spint_t max_tprocs) /* sp_func spec_tproc_set_tprocs */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  (*tproc)->max_tprocs = max_tprocs;

  (*tproc)->tprocs = func_s;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tprocs(spec_tproc_t *tproc, spec_tprocs_count_f *func_count_db, spec_tprocs_count_f *func_count_ip, spec_tprocs_rearrange_db_f *func_rearrange_db, spec_tprocs_rearrange_ip_f *func_rearrange_ip, spec_tprocs_indices_db_f *func_indices_db) /* sp_func spec_tproc_set_ext_tprocs */
{
  (*tproc)->tprocs_count_db = func_count_db;
  (*tproc)->tprocs_count_ip = func_count_ip;
  (*tproc)->tprocs_rearrange_db = func_rearrange_db;
  (*tproc)->tprocs_rearrange_ip = func_rearrange_ip;
  (*tproc)->tprocs_indices_db = func_indices_db;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tprocs_mod(spec_tproc_t *tproc, spec_tprocs_mod_f *func_s_mod, spint_t max_tprocs) /* sp_func spec_tproc_set_tprocs_mod */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  (*tproc)->max_tprocs = max_tprocs;

  (*tproc)->tprocs_mod = func_s_mod;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tprocs_mod(spec_tproc_t *tproc, spec_tprocs_mod_count_f *func_count_db, spec_tprocs_mod_count_f *func_count_ip, spec_tprocs_mod_rearrange_db_f *func_rearrange_db, spec_tprocs_mod_rearrange_ip_f *func_rearrange_ip) /* sp_func spec_tproc_set_ext_tprocs_mod */
{
  (*tproc)->tprocs_mod_count_db = func_count_db;
  (*tproc)->tprocs_mod_count_ip = func_count_ip;
  (*tproc)->tprocs_mod_rearrange_db = func_rearrange_db;
  (*tproc)->tprocs_mod_rearrange_ip = func_rearrange_ip;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_reset(spec_tproc_t *tproc, spec_tproc_reset_f *reset) /* sp_func spec_tproc_set_reset */
{
  (*tproc)->reset = reset;
  
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


spint_t spec_tproc_set_proclists(spec_tproc_t *tproc, spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs, int size, int rank, MPI_Comm comm) /* sp_func spec_tproc_set_proclists */
{
  spint_t i;

  if ((*tproc)->send_procs) z_free((*tproc)->send_procs);
  if ((*tproc)->recv_procs) z_free((*tproc)->recv_procs);

  (*tproc)->nsend_procs = -1;
  (*tproc)->send_procs = NULL;
  (*tproc)->nrecv_procs = -1;
  (*tproc)->recv_procs = NULL;

  if (nsend_procs >= 0)
  {
    (*tproc)->nsend_procs = nsend_procs;
    (*tproc)->send_procs = z_alloc(nsend_procs, sizeof(sproc_t));

    for (i = 0; i < nsend_procs; ++i) (*tproc)->send_procs[i] = send_procs[i];
  }

  if (nrecv_procs >= 0)
  {
    (*tproc)->nrecv_procs = nrecv_procs;
    (*tproc)->recv_procs = z_alloc(nrecv_procs, sizeof(sproc_t));

    for (i = 0; i < nrecv_procs; ++i) (*tproc)->recv_procs[i] = recv_procs[i];

  } else if (nsend_procs >= 0)
  {
    if (comm == MPI_COMM_NULL) return 1;

    spec_make_recv_proclist(nsend_procs, send_procs, &(*tproc)->nrecv_procs, &(*tproc)->recv_procs, size, rank, comm);
  }

/*  printf("%d: send_procs (%" spint_fmt ") = ", rank, (*tproc)->nsend_procs);
  for (i = 0; i < (*tproc)->nsend_procs; ++i) printf("  %" sproc_fmt, (*tproc)->send_procs[i]);
  printf("\n");

  printf("%d: recv_procs (%" spint_fmt ") = ", rank, (*tproc)->nrecv_procs);
  for (i = 0; i < (*tproc)->nrecv_procs; ++i) printf("  %" sproc_fmt, (*tproc)->recv_procs[i]);
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
      tproc->tprocs_mod(spec_elem_get_buf(b), i, tproc_data, &n, procs, mods);
      printf("%" spint_fmt ":", i);
      for (j = 0; j < n; ++j) printf(" %" spec_proc_fmt, procs[j]);
      printf("\n");
    }
  }
  
  /* free tproc buffers */
  spec_tproc_release(&procs, &mods);

  return 0;
}
