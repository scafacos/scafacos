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


/*#define SPEC_PRINT*/


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


spint_t spec_alltoallv_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_alltoallv_db */
{
  spint_t exit_code = SPEC_EXIT_SUCCESS;

  spint_t i;

  spec_proc_t *procs = NULL;
  int *scounts, *sdispls, *rcounts, *rdispls;

  spint_t tprocs_max;
  spec_elem_t _sd, *sd, _xb;

  spint_t stotal, rtotal;

  spint_t local_buffer_exit, global_buffer_exit;

  SPEC_DECLARE_TPROC_COUNT_DB
  SPEC_DECLARE_TPROC_MOD_COUNT_DB
  SPEC_DECLARE_TPROCS_COUNT_DB
  SPEC_DECLARE_TPROCS_MOD_COUNT_DB
  SPEC_DECLARE_TPROC_REARRANGE_DB
  SPEC_DECLARE_TPROC_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_REARRANGE_DB
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB

#ifdef SPEC_PROCLIST
  int *scounts2, *rcounts2;
#endif

#ifdef Z_PACK_TIMING
  double t[5] = { 0, 0, 0, 0, 0 };
#endif


  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  scounts = z_alloc(4 * size, sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

  /* determine max. duplicates */

  tprocs_max = 0;

  if (tproc->tprocs) tprocs_max = tproc->tprocs(spec_elem_get_buf(sb), -1, tproc_data, NULL);
  else if (tproc->tprocs_mod) tprocs_max = tproc->tprocs_mod(spec_elem_get_buf(sb), -1, tproc_data, NULL, NULL);

  sd = &_sd;

  if (tprocs_max < 0)
  {
    tprocs_max *= -1;
    sd = NULL;
  }

  if (tprocs_max) procs = z_alloca(tprocs_max, sizeof(spec_proc_t));

  if (tproc->tproc_mod || tproc->tprocs_mod)
  {
    spec_elem_alloc_buf(&_sd, z_max(1, tprocs_max));
    spec_elem_copy_type(sb, &_sd);
  }

#ifdef SPEC_PRINT
  printf("input\n");
  spec_print(&tproc, tproc_data, sb);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  for (i = 0; i < size; ++i) scounts[i] = 0;

  if (tproc->tproc)
  {
    if (tproc->tproc_count_db) tproc->tproc_count_db(sb, tproc_data, scounts);
    else SPEC_DO_TPROC_COUNT_DB(tproc->tproc, tproc_data, sb, scounts);

  } else if (tproc->tproc_mod)
  {
    if (tproc->tproc_mod_count_db) tproc->tproc_mod_count_db(sb, tproc_data, scounts);
    else SPEC_DO_TPROC_MOD_COUNT_DB(tproc->tproc_mod, tproc_data, sb, scounts);

  } else if (tproc->tprocs)
  {
    if (tproc->tprocs_count_db) tproc->tprocs_count_db(sb, tproc_data, scounts, procs);
    else SPEC_DO_TPROCS_COUNT_DB(tproc->tprocs, tproc_data, sb, scounts, procs);

  } else if (tproc->tprocs_mod)
  {
    if (tproc->tprocs_mod_count_db) tproc->tprocs_mod_count_db(sb, tproc_data, scounts, procs);
    else SPEC_DO_TPROCS_MOD_COUNT_DB(tproc->tprocs_mod, tproc_data, sb, scounts, procs);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

#ifdef SPEC_PRINT
  printf("after count\n");
  spec_print(&tproc, tproc_data, sb);
#endif

  /* redistribute local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[2]);

#ifdef SPEC_PROCLIST
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    scounts2 = z_alloc(2 * size, sizeof(int));
    rcounts2 = scounts2 + 1 * size;

    for (i = 0; i < size; ++i)
    {
      scounts2[i] = rcounts2[i] = 0;
      sdispls[i] = rdispls[i] = i;
      rcounts[i] = 0;
    }

    for (i = 0; i < tproc->nsend_procs; ++i) scounts2[tproc->send_procs[i]] = 1;
    for (i = 0; i < tproc->nrecv_procs; ++i) rcounts2[tproc->recv_procs[i]] = 1;

/*    printf("%d: scounts2 = ", rank);
    for (i = 0; i < size; ++i) printf("  %d (%d)", scounts2[i], scounts[i]);
    printf("\n");*/

#ifdef HAVE_ZMPI_ALLTOALLV_PROCLISTS
    ZMPI_Alltoallv_proclists(scounts, scounts2, sdispls, MPI_INT, tproc->nsend_procs, tproc->send_procs, rcounts, rcounts2, rdispls, MPI_INT, tproc->nrecv_procs, tproc->recv_procs, comm);
#else
    MPI_Alltoallv(scounts, scounts2, sdispls, MPI_INT, rcounts, rcounts2, rdispls, MPI_INT, comm);
#endif

/*    printf("%d: rcounts2 = ", rank);
    for (i = 0; i < size; ++i) printf("  %d (%d)", rcounts2[i], rcounts[i]);
    printf("\n");*/
    
    z_free(scounts2);

  } else
#endif
  {
#ifdef HAVE_ZMPI_ALLTOALL_2STEP
# ifdef SPEC_MPI_ALLTOALL_2STEP_THRESHOLD
    if (size >= SPEC_MPI_ALLTOALL_2STEP_THRESHOLD)
      ZMPI_Alltoall_2step_int(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
    else
# endif
#endif
      MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[2]);

  sdispls[0] = rdispls[0] = 0;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  }

  stotal = sdispls[size - 1] + scounts[size - 1];
  rtotal = rdispls[size - 1] + rcounts[size - 1];

  local_buffer_exit = (rtotal > spec_elem_get_nmax(rb))?1:0;

#ifdef spec_elem_alloc_rbuf
  if (local_buffer_exit > 0 && spec_elem_alloc_rbuf(rb))
  {
    spec_elem_alloc_buf(rb, rtotal);
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
      "%d: spec_alltoallv_db: error: one receive buffer is too small (local: %" spint_fmt " vs. %" spint_fmt")\n",
# else
      "%d: spec_alltoallv_db: error: local receive buffer too small (%" spint_fmt " vs. %" spint_fmt")\n",
# endif
      rank, (spint_t) spec_elem_get_nmax(rb), rtotal);
#endif
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
    if (tproc->tproc_mod_rearrange_db) tproc->tproc_mod_rearrange_db(sb, xb, tproc_data, sdispls, sd);
    else SPEC_DO_TPROC_MOD_REARRANGE_DB(tproc->tproc_mod, tproc_data, sb, xb, sdispls, sd);

  } else if (tproc->tprocs)
  {
    if (tproc->tprocs_rearrange_db) tproc->tprocs_rearrange_db(sb, xb, tproc_data, sdispls, procs);
    else SPEC_DO_TPROCS_REARRANGE_DB(tproc->tprocs, tproc_data, sb, xb, sdispls, procs);

  } else if (tproc->tprocs_mod)
  {
    if (tproc->tprocs_mod_rearrange_db) tproc->tprocs_mod_rearrange_db(sb, xb, tproc_data, sdispls, procs, sd);
    else SPEC_DO_TPROCS_MOD_REARRANGE_DB(tproc->tprocs_mod, tproc_data, sb, xb, sdispls, procs, sd);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

#ifdef SPEC_PRINT
  printf("after rearrange\n");
  spec_print(&tproc, tproc_data, xb);
#endif

  /* alltoallv */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[4]);

  sdispls[0] = 0;
  for (i = 1; i < size; ++i) sdispls[i] = sdispls[i - 1] + scounts[i - 1];

#ifdef SPEC_PROCLIST
# ifdef spec_elem_alltoallv_proclists_db
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    spec_elem_alltoallv_proclists_db(xb, scounts, sdispls, tproc->nsend_procs, tproc->send_procs, rb, rcounts, rdispls, tproc->nrecv_procs, tproc->recv_procs, size, rank, comm);
  }
  else
# endif
#endif
  {
    spec_elem_alltoallv_db(xb, scounts, sdispls, rb, rcounts, rdispls, size, rank, comm);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  spec_elem_set_n(rb, rtotal);

  if (xb == &_xb) spec_elem_free_buf(&_xb);

free_and_exit:

  if (tprocs_max) z_freea(procs);

  if (tproc->tproc_mod || tproc->tprocs_mod) spec_elem_free_buf(&_sd);

  z_free(scounts);

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
  int *scounts, *sdispls, *rcounts, *rdispls;

  spint_t tprocs_max;
  spec_elem_t _sd, *sd, _xb, *tb;

  spint_t stotal, rtotal;

  spint_t local_buffer_exit, global_buffer_exit;
  spint_t local_alltoallv_inplace, global_alltoallv_inplace;

  SPEC_DECLARE_TPROC_COUNT_IP
  SPEC_DECLARE_TPROC_REARRANGE_IP
  SPEC_DECLARE_TPROC_REARRANGE_DB
  SPEC_DECLARE_TPROC_MOD_COUNT_IP
  SPEC_DECLARE_TPROC_MOD_REARRANGE_IP
  SPEC_DECLARE_TPROC_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_COUNT_IP
  SPEC_DECLARE_TPROCS_REARRANGE_DB
  SPEC_DECLARE_TPROCS_REARRANGE_IP
  SPEC_DECLARE_TPROCS_MOD_COUNT_IP
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP

#ifdef SPEC_PROCLIST
  int *scounts2, *rcounts2;
#endif

#ifdef Z_PACK_TIMING
  double t[5] = { 0, 0, 0, 0, 0 };
#endif


  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  scounts = z_alloc(4 * size, sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

  /* determine max. duplicates */

  tprocs_max = 0;

  if (tproc->tprocs) tprocs_max = tproc->tprocs(spec_elem_get_buf(b), -1, tproc_data, NULL);
  else if (tproc->tprocs_mod) tprocs_max = tproc->tprocs_mod(spec_elem_get_buf(b), -1, tproc_data, NULL, NULL);

  sd = &_sd;

  if (tprocs_max < 0)
  {
    tprocs_max *= -1;
    sd = NULL;
  }

  if (tprocs_max) procs = z_alloca(tprocs_max, sizeof(spec_proc_t));

  if (tproc->tproc_mod || tproc->tprocs_mod)
  {
    spec_elem_alloc_buf(&_sd, z_max(1, tprocs_max));
    spec_elem_copy_type(b, &_sd);
  }

#ifdef SPEC_PRINT
  printf("input\n");
  spec_print(&tproc, tproc_data, b);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  for (i = 0; i < size; ++i) scounts[i] = 0;

  if (tproc->tproc)
  {
    if (tproc->tproc_count_ip) tproc->tproc_count_ip(b, tproc_data, scounts);
    else SPEC_DO_TPROC_COUNT_IP(tproc->tproc, tproc_data, b, scounts);

  } else if (tproc->tproc_mod)
  {
    if (tproc->tproc_mod_count_ip) tproc->tproc_mod_count_ip(b, tproc_data, scounts);
    else SPEC_DO_TPROC_MOD_COUNT_IP(tproc->tproc_mod, tproc_data, b, scounts);

  } else if (tproc->tprocs)
  {
    if (tproc->tprocs_count_ip) tproc->tprocs_count_ip(b, tproc_data, scounts, procs);
    else SPEC_DO_TPROCS_COUNT_IP(tproc->tprocs, tproc_data, b, scounts, procs);

  } else if (tproc->tprocs_mod)
  {
    if (tproc->tprocs_mod_count_ip) tproc->tprocs_mod_count_ip(b, tproc_data, scounts, procs);
    else SPEC_DO_TPROCS_MOD_COUNT_IP(tproc->tprocs_mod, tproc_data, b, scounts, procs);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

#ifdef SPEC_PRINT
  printf("after count\n");
  spec_print(&tproc, tproc_data, &b);
#endif

  /* redistribute local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[2]);

#ifdef SPEC_PROCLIST
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    scounts2 = z_alloc(2 * size, sizeof(int));
    rcounts2 = scounts2 + 1 * size;

    for (i = 0; i < size; ++i)
    {
      scounts2[i] = rcounts2[i] = 0;
      sdispls[i] = rdispls[i] = i;
    }

    for (i = 0; i < tproc->nsend_procs; ++i) scounts2[tproc->send_procs[i]] = 1;
    for (i = 0; i < tproc->nrecv_procs; ++i) rcounts2[tproc->recv_procs[i]] = 1;

#ifdef HAVE_ZMPI_ALLTOALLV_PROCLISTS
    ZMPI_Alltoallv_proclists(scounts, scounts2, sdispls, MPI_INT, tproc->nsend_procs, tproc->send_procs, rcounts, rcounts2, rdispls, MPI_INT, tproc->nrecv_procs, tproc->recv_procs, comm);
#else
    MPI_Alltoallv(scounts, scounts2, sdispls, MPI_INT, rcounts, rcounts2, rdispls, MPI_INT, comm);
#endif

    z_free(scounts2);

  } else
#endif
  {
#ifdef HAVE_ZMPI_ALLTOALL_2STEP
# ifdef SPEC_MPI_ALLTOALL_2STEP_THRESHOLD
    if (size >= SPEC_MPI_ALLTOALL_2STEP_THRESHOLD)
      ZMPI_Alltoall_2step_int(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
    else
# endif
#endif
      MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[2]);

  sdispls[0] = rdispls[0] = 0;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  }

  stotal = sdispls[size - 1] + scounts[size - 1];
  rtotal = rdispls[size - 1] + rcounts[size - 1];

  local_buffer_exit = (z_max(stotal, rtotal) > spec_elem_get_nmax(b))?1:0;

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
      "%d: spec_alltoallv_ip: error: one buffer is too small (local: %" spint_fmt " vs. %" spint_fmt")\n",
# else
      "%d: spec_alltoallv_ip: error: local buffer too small (%" spint_fmt " vs. %" spint_fmt")\n",
# endif
      rank, z_max(stotal, rtotal), (spint_t) spec_elem_get_nmax(b));
#endif
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
      if (tproc->tproc_mod_rearrange_db) tproc->tproc_mod_rearrange_db(b, tb, tproc_data, sdispls, sd);
      else SPEC_DO_TPROC_MOD_REARRANGE_DB(tproc->tproc_mod, tproc_data, b, tb, sdispls, sd);

    } else
    {
      if (tproc->tproc_mod_rearrange_ip) tproc->tproc_mod_rearrange_ip(b, xb, tproc_data, sdispls, scounts, size, sd);
      else SPEC_DO_TPROC_MOD_REARRANGE_IP(tproc->tproc_mod, tproc_data, b, xb, sdispls, scounts, size, sd);
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
      if (tproc->tprocs_mod_rearrange_db) tproc->tprocs_mod_rearrange_db(b, tb, tproc_data, sdispls, procs, sd);
      else SPEC_DO_TPROCS_MOD_REARRANGE_DB(tproc->tprocs_mod, tproc_data, b, tb, sdispls, procs, sd);

    } else
    {
      if (tproc->tprocs_mod_rearrange_ip) tproc->tprocs_mod_rearrange_ip(b, xb, tproc_data, sdispls, scounts, size, procs, sd);
      else SPEC_DO_TPROCS_MOD_REARRANGE_IP(tproc->tprocs_mod, tproc_data, b, xb, sdispls, scounts, size, procs, sd);
    }
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

#ifdef SPEC_PRINT
  printf("after rearrange\n");
  spec_print(&tproc, tproc_data, &b);
#endif

  /* alltoallv */
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
    if (tb) spec_elem_alltoallv_db(tb, scounts, sdispls, b, rcounts, rdispls, size, rank, comm);
    else
    {
      spec_elem_alltoallv_db(b, scounts, sdispls, xb, rcounts, rdispls, size, rank, comm);
      spec_elem_ncopy_at(xb, 0, b, 0, rtotal);
    }

  } else
  {
#ifdef spec_elem_alltoallv_ip

    /* in-place alltoallv  */
    if (tb) spec_elem_ncopy_at(tb, 0, b, 0, stotal);
    spec_elem_alltoallv_ip(b, xb, scounts, sdispls, rcounts, rdispls, size, rank, comm);

#else

#ifdef SPEC_ERROR_FILE
    fprintf(SPEC_ERROR_FILE, "%d: spec_alltoallv_ip: error: required in-place MPI_Alltoallv support not available, use the ZMPI All-to-all In-place library!", rank);
#endif
    exit_code = SPEC_EXIT_FAILED;
    goto free_and_exit;

#endif
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  spec_elem_set_n(b, rtotal);

  if (xb == &_xb) spec_elem_free_buf(&_xb);

free_and_exit:

  if (tprocs_max) z_freea(procs);

  if (tproc->tproc_mod || tproc->tprocs_mod) spec_elem_free_buf(&_sd);

  z_free(scounts);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[0]);

  Z_TIMING_PRINT(0, __func__, sizeof(t) / sizeof(t[0]), t, rank);

#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
  if (spec_timing)
    for (i = 0; i < sizeof(t) / sizeof(t[0]); ++i) spec_timing[i] = t[i];
#endif

  return exit_code;
}


#undef PROCLIST_IFELSE



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "spec_core.h"


#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
double *spec_timing = NULL; /* sp_var spec_timing */
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


spint_t spec_tproc_create(spec_tproc_t *tproc, spec_tproc_f *func, spec_tproc_mod_f *func_mod, spec_tprocs_f *func_s, spec_tprocs_mod_f *func_s_mod) /* sp_func spec_tproc_create */
{
  *tproc = z_alloc(1, sizeof(struct _spec_tproc_t));

  (*tproc)->tproc = func;
  spec_tproc_unset_ext_tproc(tproc);

  (*tproc)->tproc_mod = func_mod;
  spec_tproc_unset_ext_tproc_mod(tproc);

  (*tproc)->tprocs = func_s;
  spec_tproc_unset_ext_tprocs(tproc);

  (*tproc)->tprocs_mod = func_s_mod;
  spec_tproc_unset_ext_tprocs_mod(tproc);

  (*tproc)->reset = NULL;
  
#ifdef SPEC_PROCLIST
  (*tproc)->nsend_procs = -1;
  (*tproc)->send_procs = NULL;
  (*tproc)->nrecv_procs = -1;
  (*tproc)->recv_procs = NULL;
#endif

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_destroy(spec_tproc_t *tproc) /* sp_func spec_tproc_destroy */
{
#ifdef SPEC_PROCLIST
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


spint_t spec_tproc_set_ext_tproc(spec_tproc_t *tproc, spec_tproc_count_f *func_count_db, spec_tproc_count_f *func_count_ip, spec_tproc_rearrange_db_f *func_rearrange_db, spec_tproc_rearrange_ip_f *func_rearrange_ip) /* sp_func spec_tproc_set_ext_tproc */
{
  (*tproc)->tproc_count_db = func_count_db;
  (*tproc)->tproc_count_ip = func_count_ip;
  (*tproc)->tproc_rearrange_db = func_rearrange_db;
  (*tproc)->tproc_rearrange_ip = func_rearrange_ip;

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


spint_t spec_tproc_set_tprocs(spec_tproc_t *tproc, spec_tprocs_f *func_s) /* sp_func spec_tproc_set_tprocs */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  (*tproc)->tprocs = func_s;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tprocs(spec_tproc_t *tproc, spec_tprocs_count_f *func_count_db, spec_tprocs_count_f *func_count_ip, spec_tprocs_rearrange_db_f *func_rearrange_db, spec_tprocs_rearrange_ip_f *func_rearrange_ip) /* sp_func spec_tproc_set_ext_tprocs */
{
  (*tproc)->tprocs_count_db = func_count_db;
  (*tproc)->tprocs_count_ip = func_count_ip;
  (*tproc)->tprocs_rearrange_db = func_rearrange_db;
  (*tproc)->tprocs_rearrange_ip = func_rearrange_ip;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tprocs_mod(spec_tproc_t *tproc, spec_tprocs_mod_f *func_s_mod) /* sp_func spec_tproc_set_tprocs_mod */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

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


#ifdef SPEC_PROCLIST

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


spint_t spec_print(spec_tproc_t *tproc, spec_tproc_data_t tproc_data, spec_elem_t *b) /* sp_func spec_print */
{
  spint_t i, j, n;
  spec_proc_t p, *procs = NULL;

  spint_t tprocs_max, tprocs_mod;
  spec_elem_t sd;
  
  
  tprocs_max = 0;
  tprocs_mod = 1;

  if ((*tproc)->tprocs) tprocs_max = (*tproc)->tprocs(spec_elem_get_buf(b), -1, tproc_data, NULL);
  else if ((*tproc)->tprocs_mod) tprocs_max = (*tproc)->tprocs_mod(spec_elem_get_buf(b), -1, tproc_data, NULL, NULL);

  if (tprocs_max < 0)
  {
    tprocs_max *= -1;
    tprocs_mod = 0;
  }

  if (tprocs_max)
  {
    procs = z_alloca(tprocs_max, sizeof(spec_proc_t));

    if ((*tproc)->tprocs_mod) spec_elem_alloc_buf(&sd, tprocs_max);
  }

  /* reset tproc if necessary */
  if ((*tproc)->reset) (*tproc)->reset(tproc_data);


  if ((*tproc)->tproc)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      p = (*tproc)->tproc(spec_elem_get_buf(b), i, tproc_data);
      printf("%" spint_fmt ": %" spec_proc_fmt "\n", i, p);
    }

  } else if ((*tproc)->tprocs)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      n = (*tproc)->tprocs(spec_elem_get_buf(b), i, tproc_data, procs);
      printf("%" spint_fmt ":", i);
      for (j = 0; j < n; ++j) printf(" %" spec_proc_fmt, procs[j]);
      printf("\n");
    }

  } else if ((*tproc)->tprocs_mod)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      n = (*tproc)->tprocs_mod(spec_elem_get_buf(b), i, tproc_data, procs, (tprocs_mod)?&sd:NULL);
      printf("%" spint_fmt ":", i);
      for (j = 0; j < n; ++j) printf(" %" spec_proc_fmt, procs[j]);
      printf("\n");
    }
  }
  
  if (tprocs_max)
  {
    z_freea(procs);

    if ((*tproc)->tprocs_mod) spec_elem_free_buf(&sd);
  }

  return 0;
}
