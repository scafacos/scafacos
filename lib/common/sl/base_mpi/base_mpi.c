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


/* sl_macro MB_REDUCEBCAST_THRESHOLD */
/* sl_macro MB_TRACE_IF */


#include "sl_common.h"


#define REDUCEBCAST_ROOT  0


#if !defined(MB_REDUCEBCAST_THRESHOLD) && defined(GLOBAL_REDUCEBCAST_THRESHOLD)
# define MB_REDUCEBCAST_THRESHOLD  GLOBAL_REDUCEBCAST_THRESHOLD
#endif

#ifndef MB_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MB_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MB_TRACE_IF  (SL_PROC_RANK == -1)
# endif
#endif


slint_t mpi_binning_create(global_bins_t *gb, slint_t max_nbins, slint_t max_nbinnings, elements_t *s, slint_t nelements, slint_t docounts, slint_t doweights, binning_t *bm, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_create */
{
  gb->bm = bm;

  binning_create(&gb->lb, max_nbins, max_nbinnings, s, nelements, docounts, doweights, bm);

  gb->bcws = gb->lb.bcws;

#if defined(elem_weight) && defined(sl_weight_intequiv)
  gb->cws = z_alloc(gb->lb.max_nbinnings * gb->lb.cw_factor * gb->bm->max_nbins, sizeof(slweight_t));
  gb->prefix_cws = z_alloc(gb->lb.max_nbinnings * gb->lb.cw_factor, sizeof(slweight_t));
#else
  gb->cs = z_alloc(gb->lb.max_nbinnings * 1 * gb->bm->max_nbins, sizeof(slint_t));
  gb->prefix_cs = z_alloc(gb->lb.max_nbinnings * 1, sizeof(slint_t));
# ifdef elem_weight
  if (gb->lb.doweights)
  {
    gb->ws = z_alloc(gb->lb.max_nbinnings * 1 * gb->bm->max_nbins, sizeof(slweight_t));
    gb->prefix_ws = z_alloc(gb->lb.max_nbinnings * 1, sizeof(slweight_t));

  } else gb->ws = gb->prefix_ws = NULL;
# endif
#endif

  return 0;
}


slint_t mpi_binning_destroy(global_bins_t *gb, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_destroy */
{
  binning_destroy(&gb->lb);

#if defined(elem_weight) && defined(sl_weight_intequiv)
  z_free(gb->cws);
#else
  z_free(gb->cs);
# ifdef elem_weight
  if (gb->lb.doweights)
    z_free(gb->ws);
# endif
#endif

  return 0;
}


slint_t mpi_binning_pre(global_bins_t *gb, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_pre */
{
  binning_pre(&gb->lb);

  return 0;
}


slint_t mpi_binning_exec_reset(global_bins_t *gb, slint_t do_bins, slint_t do_prefixes, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_exec_reset */
{
  binning_exec_reset(&gb->lb, do_bins, do_prefixes);

  return 0;
}


slint_t mpi_binning_exec_local(global_bins_t *gb, slint_t b, slint_t do_bins, slint_t do_prefixes, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_exec_local */
{
  return binning_exec(&gb->lb, b, do_bins, do_prefixes);
}


slint_t mpi_binning_exec_global(global_bins_t *gb, slint_t do_bins, slint_t do_prefixes, slint_t root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_exec_global */
{
  slint_t b, k;


  if (do_bins)
  {
#if defined(elem_weight) && defined(sl_weight_intequiv)
    Z_TRACE_IF(MB_TRACE_IF, "%sreducing %" slint_fmt " ints/weights", ((root < 0)?"all-":""), (gb->lb.nbinnings * gb->lb.cw_factor * gb->bm->nbins));
#else
    Z_TRACE_IF(MB_TRACE_IF, "%sreducing %" slint_fmt " ints", ((root < 0)?"all-":""), (gb->lb.nbinnings * 1 * gb->bm->nbins));
# ifdef elem_weight
    if (gb->lb.doweights)
      Z_TRACE_IF(MB_TRACE_IF, "%sreducing %" slint_fmt " weights", ((root < 0)?"all-":""), (gb->lb.nbinnings * 1 * gb->bm->nbins));
# endif
#endif

    if (root < 0)
    {
#if defined(elem_weight) && defined(sl_weight_intequiv)
      sl_MPI_Allreduce(gb->lb.cws, gb->cws, gb->lb.nbinnings * gb->lb.cw_factor * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, comm, size, rank);
#else
      sl_MPI_Allreduce(gb->lb.cs, gb->cs, gb->lb.nbinnings * 1 * gb->bm->nbins, int_mpi_datatype, MPI_SUM, comm, size, rank);
# ifdef elem_weight
      if (gb->lb.doweights)
        sl_MPI_Allreduce(gb->lb.ws, gb->ws, gb->lb.nbinnings * 1 * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, comm, size, rank);
# endif
#endif

    } else
    {
#if defined(elem_weight) && defined(sl_weight_intequiv)
      MPI_Reduce(gb->lb.cws, gb->cws, gb->lb.nbinnings * gb->lb.cw_factor * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, root, comm);
#else
      MPI_Reduce(gb->lb.cs, gb->cs, gb->lb.nbinnings * 1 * gb->bm->nbins, int_mpi_datatype, MPI_SUM, root, comm);
# ifdef elem_weight
      if (gb->lb.doweights)
        MPI_Reduce(gb->lb.ws, gb->ws, gb->lb.nbinnings * 1 * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, root, comm);
# endif
#endif
    }

    if (root < 0 || root == rank)
    {
      for (b = 0; b < gb->lb.nbinnings; ++b)
      {
        Z_TRACE_ARRAY_IF(MB_TRACE_IF, k, gb->bm->nbins, " %" slcount_fmt, gb_counts(gb, b, 0)[k], "%" slint_fmt ": counts =", b);
#ifdef elem_weight
        if (gb->lb.doweights)
          Z_TRACE_ARRAY_IF(MB_TRACE_IF, k, gb->bm->nbins, " %" slweight_fmt, gb_weights(gb, b, 0)[k], "%" slint_fmt ": weights =", b);
#endif
      }
    }
  }

  if (do_prefixes)
  {
#if defined(elem_weight) && defined(sl_weight_intequiv)
    Z_TRACE_IF(MB_TRACE_IF, "exscan with %" slint_fmt " ints/weights", (gb->lb.nbinnings * gb->lb.cw_factor));
#else
    if (gb->lb.docounts)
      Z_TRACE_IF(MB_TRACE_IF, "exscan with %" slint_fmt " ints", (gb->lb.nbinnings * 1));
# ifdef elem_weight
    if (gb->lb.doweights)
      Z_TRACE_IF(MB_TRACE_IF, "exscan with %" slint_fmt " weights", (gb->lb.nbinnings * 1));
# endif
#endif

    if (gb->lb.docounts)
      Z_TRACE_ARRAY_IF(MB_TRACE_IF, b, gb->lb.nbinnings, " %" slcount_fmt, lb_prefix_count(&gb->lb, b)[0], "local prefix_counts =");
#ifdef elem_weight
    if (gb->lb.doweights)
      Z_TRACE_ARRAY_IF(MB_TRACE_IF, b, gb->lb.nbinnings, " %" slweight_fmt, lb_prefix_weight(&gb->lb, b)[0], "local prefix_weights =");
#endif

#if defined(elem_weight) && defined(sl_weight_intequiv)
    MPI_Exscan(gb->lb.prefix_cws, gb->prefix_cws, gb->lb.nbinnings * gb->lb.cw_factor, weight_mpi_datatype, MPI_SUM, comm);
    if (rank == 0) for (b = 0; b < gb->lb.nbinnings * gb->lb.cw_factor; ++b) gb->prefix_cws[b] = 0;
#else
    if (gb->lb.docounts)
    {
      MPI_Exscan(gb->lb.prefix_cs, gb->prefix_cs, gb->lb.nbinnings * 1, int_mpi_datatype, MPI_SUM, comm);
      if (rank == 0) for (b = 0; b < gb->lb.nbinnings * 1; ++b) gb->prefix_cs[b] = 0;
    }
# ifdef elem_weight
    if (gb->lb.doweights)
    {
      MPI_Exscan(gb->lb.prefix_ws, gb->prefix_ws, gb->lb.nbinnings * 1, weight_mpi_datatype, MPI_SUM, comm);
      if (rank == 0) for (b = 0; b < gb->lb.nbinnings * 1; ++b) gb->prefix_ws[b] = 0;
    }
# endif
#endif

    if (gb->lb.docounts)
      Z_TRACE_ARRAY_IF(MB_TRACE_IF, b, gb->lb.nbinnings, " %" slcount_fmt, gb_prefix_count(gb, b)[0], "global prefix_counts =");
#ifdef elem_weight
    if (gb->lb.doweights)
      Z_TRACE_ARRAY_IF(MB_TRACE_IF, b, gb->lb.nbinnings, " %" slweight_fmt, gb_prefix_weight(gb, b)[0], "global prefix_weights =");
#endif
  }

  return 0;
}


slint_t mpi_binning_refine(global_bins_t *gb, slint_t b, slint_t k, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_refine */
{
  return binning_refine(&gb->lb, b, k, sp, s);
}


slint_t mpi_binning_hit(global_bins_t *gb, slint_t b, slint_t k, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_hit */
{
  binning_hit(&gb->lb, b, k, sp, s);
  
  return 0;
}


slint_t mpi_binning_finalize(global_bins_t *gb, slint_t b, slint_t dc, slweight_t dw, slint_t lc_min, slint_t lc_max, slcount_t *lcs, slweight_t *lws, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_finalize */
{
  binning_finalize(&gb->lb, b, dc, dw, lc_min, lc_max, lcs, lws, sp, s);
  
  return 0;
}


slint_t mpi_binning_post(global_bins_t *gb, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_post */
{
  binning_post(&gb->lb);

  return 0;
}


#undef REDUCEBCAST_ROOT



#include "sl_common.h"

#ifdef HAVE_SPI_KERNEL_INTERFACE_H
# include <spi/kernel_interface.h>
#endif
#ifdef HAVE_COMMON_BGP_PERSONALITY_H
# include <common/bgp_personality.h>
#endif
#ifdef HAVE_COMMON_BGP_PERSONALITY_INLINES_H
# include <common/bgp_personality_inlines.h>
#endif

/* sl_ifdef SL_USE_MPI sl_global sl_context CONTEXT_BEGIN mpi */
const MPI_Datatype default_mpi_int_datatype = MPI_DATATYPE_NULL;      /* sl_global sl_context sl_var default_mpi_int_datatype */
const MPI_Datatype default_mpi_key_datatype = MPI_DATATYPE_NULL;      /* sl_global sl_context sl_var default_mpi_key_datatype */
const MPI_Datatype default_mpi_pkey_datatype = MPI_DATATYPE_NULL;     /* sl_global sl_context sl_var default_mpi_pkey_datatype */
const MPI_Datatype default_mpi_pwkey_datatype = MPI_DATATYPE_NULL;    /* sl_global sl_context sl_var default_mpi_pwkey_datatype */
const MPI_Datatype default_mpi_index_datatype = MPI_DATATYPE_NULL;    /* sl_global sl_context sl_var default_mpi_index_datatype */
const MPI_Datatype default_mpi_weight_datatype = MPI_DATATYPE_NULL;   /* sl_global sl_context sl_var default_mpi_weight_datatype */
const MPI_Datatype default_mpi_data_datatype[SL_DATA_NMAX + 1] =      /* sl_global sl_context sl_var default_mpi_data_datatype */
{                                                                     /* sl_context */
#define xelem_call_data      MPI_DATATYPE_NULL,                       /* sl_context */
#define xelem_call_data_not  MPI_DATATYPE_NULL,                       /* sl_context */
#include "sl_xelem_call.h"                                            /* sl_context */
  MPI_DATATYPE_NULL                                                   /* sl_context */
};                                                                    /* sl_context */

const int default_mpi_rank = -2;  /* sl_global sl_context sl_var default_mpi_rank */
/* sl_endif sl_global sl_context CONTEXT_END mpi */


slint_t mpi_datatypes_init() /* sl_proto, sl_func mpi_datatypes_init */
{
  slpwkey_t pwk;
  int pwk_blenghts[2] = { key_pure_size_mpi, 1 };
  MPI_Datatype pwk_types[2] = { key_pure_type_mpi, sl_weight_type_mpi };
  MPI_Aint pwk_displs[2], base;

  /* intern integer */
  if (sl_int_size_mpi > 1)
  {
    MPI_Type_contiguous(sl_int_size_mpi, sl_int_type_mpi, &int_mpi_datatype);
    MPI_Type_commit(&int_mpi_datatype);

  } else int_mpi_datatype = sl_int_type_mpi;

#define xelem_call \
  if (xelem_size_mpi > 1) \
  { \
    MPI_Type_contiguous(xelem_size_mpi, xelem_type_mpi, &xelem_mpi_datatype); \
    MPI_Type_commit(&xelem_mpi_datatype); \
\
  } else xelem_mpi_datatype = xelem_type_mpi;
#include "sl_xelem_call.h"

  /* pure key */
  if (key_pure_size_mpi > 1)
  {
    MPI_Type_contiguous(key_pure_size_mpi, key_pure_type_mpi, &pkey_mpi_datatype);
    MPI_Type_commit(&pkey_mpi_datatype);

  } else pkey_mpi_datatype = key_pure_type_mpi;

  /* pure weighted key */
  MPI_Get_address(&pwk, &base);
  MPI_Get_address(&pwk.pkey, &pwk_displs[0]); pwk_displs[0] -= base;
  MPI_Get_address(&pwk.weight, &pwk_displs[1]); pwk_displs[1] -= base;

  MPI_Type_create_struct(2, pwk_blenghts, pwk_displs, pwk_types, &pwkey_mpi_datatype);
  MPI_Type_commit(&pwkey_mpi_datatype);

  /* weight (intern) */
  if (sl_weight_size_mpi > 1)
  {
    MPI_Type_contiguous(sl_weight_size_mpi, sl_weight_type_mpi, &weight_mpi_datatype);
    MPI_Type_commit(&weight_mpi_datatype);

  } else weight_mpi_datatype = sl_weight_type_mpi;

  return 0;
}


slint_t mpi_datatypes_release() /* sl_proto, sl_func mpi_datatypes_release */
{
  if (int_mpi_datatype != sl_int_type_mpi) MPI_Type_free(&int_mpi_datatype);

#define xelem_call \
  if (xelem_mpi_datatype != xelem_type_mpi) MPI_Type_free(&xelem_mpi_datatype);
#include "sl_xelem_call.h"

  if (pkey_mpi_datatype != key_pure_type_mpi) MPI_Type_free(&pkey_mpi_datatype);

  MPI_Type_free(&pwkey_mpi_datatype);

  if (weight_mpi_datatype != sl_weight_type_mpi) MPI_Type_free(&weight_mpi_datatype);

  return 0;
}


slint_t mpi_get_grid_properties(slint_t ndims, slint_t *dims, slint_t *pos, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_get_grid_properties */
{
  const slint_t max_ndims = 4;

  slint_t i, _dims[max_ndims], _pos[max_ndims], _rdims[max_ndims], _rpos[max_ndims];

  slint_t remap[max_ndims];

#ifdef HAVE__BGP_PERSONALITY_T

  _BGP_Personality_t personality;

  Kernel_GetPersonality(&personality, sizeof(personality));

  _dims[0] = personality.Network_Config.Xnodes;
  _dims[1] = personality.Network_Config.Ynodes;
  _dims[2] = personality.Network_Config.Znodes;

  _pos[0] = personality.Network_Config.Xcoord;
  _pos[1] = personality.Network_Config.Ycoord;
  _pos[2] = personality.Network_Config.Zcoord;

  _pos[3] = Kernel_PhysicalProcessorID();

  switch (personality.Kernel_Config.ProcessConfig)
  {
    case _BGP_PERS_PROCESSCONFIG_SMP:
      _dims[3] = 1;
      break;
    case _BGP_PERS_PROCESSCONFIG_VNM:
      _dims[3] = 4;
      break;
    case _BGP_PERS_PROCESSCONFIG_2x2:
      _dims[3] = 2;
      break;
    default:
      _dims[3] = 1;
      break;
  }

#else

/*#define STATIC*/

# ifndef STATIC

  int mpi_dims[4] = { 0, 0, 0, 0 };

  MPI_Dims_create(size, ndims, mpi_dims);

  for (i = 0; i < ndims; ++i) _dims[i] = mpi_dims[ndims - i - 1];
  for (i = ndims; i < 4; ++i) _dims[i] = 1;

# else
  switch (size)
  {
    case 4:   _dims[0] = 2;    _dims[1] = 2; _dims[2] = 1; _dims[3] = 1; break;
    case 8:   _dims[0] = 2;    _dims[1] = 2; _dims[2] = 2; _dims[3] = 1; break;
/*    case 16:  _dims[0] = 2;    _dims[1] = 2; _dims[2] = 2; _dims[3] = 2; break;*/
    case 16:  _dims[0] = 4;    _dims[1] = 2; _dims[2] = 2; _dims[3] = 1; break;
    case 32:  _dims[0] = 4;    _dims[1] = 4; _dims[2] = 2; _dims[3] = 1; break;
/*    case 32:  _dims[0] = 8;    _dims[1] = 4; _dims[2] = 1; _dims[3] = 1; break;*/
    case 64:  _dims[0] = 4;    _dims[1] = 4; _dims[2] = 4; _dims[3] = 1; break;
/*    case 64:  _dims[0] = 8;    _dims[1] = 8; _dims[2] = 1; _dims[3] = 1; break;*/
    case 128: _dims[0] = 8;    _dims[1] = 4; _dims[2] = 4; _dims[3] = 1; break;
    case 512: _dims[0] = 8;    _dims[1] = 8; _dims[2] = 8; _dims[3] = 1; break;
    default:  _dims[0] = size; _dims[1] = 1; _dims[2] = 1; _dims[3] = 1; break;
  }
# endif

#undef STATIC

  _pos[3] = (rank / (1))                                  % _dims[3];
  _pos[2] = (rank / (1 * _dims[3]))                       % _dims[2];
  _pos[1] = (rank / (1 * _dims[3] * _dims[2]))            % _dims[1];
  _pos[0] = (rank / (1 * _dims[3] * _dims[2] * _dims[1])) % _dims[0];

#endif

  for (i = 0; i < max_ndims; ++i) remap[i] = i;

  /* remap dimensions */
  for (i = 0; i < max_ndims; ++i)
  {
    _rdims[i] = _dims[remap[i]];
    _rpos[i] = _pos[remap[i]];
  }

  /* reduce dimensions (lowest significant first) */
  for (i = max_ndims - 1; i >= ndims; --i)
  {
    _rpos[i - 1] = _rpos[i - 1] * _rdims[i] + _rpos[i];
    _rdims[i - 1] *= _rdims[i];

    _rpos[i] = 0;
    _rdims[i] = 1;
  }

  for (i = 0; i < max_ndims; ++i)
  {
    dims[i] = _rdims[i];
    pos[i] = _rpos[i];
  }

/*  printf("%d: <%" slint_fmt ",%" slint_fmt ",%" slint_fmt ",%" slint_fmt "> of <%" slint_fmt "x%" slint_fmt "x%" slint_fmt"x%" slint_fmt ">\n", rank, pos[0], pos[1], pos[2], pos[3], dims[0], dims[1], dims[2], dims[3]);*/

  return 0;
}


slint_t mpi_subgroups_create(slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_subgroups_create */
{
  slint_t i, nrtsize;


  nrtsize = (slint_t) pow((double) size, 1.0 / (double) nsubgroups);

  for (i = 0; i < nsubgroups; ++i)
  {
    if (i == 0) MPI_Comm_dup(comm, &sub_comms[i]);
    else MPI_Comm_split(sub_comms[i - 1], sub_ranks[i - 1] * nrtsize / sub_sizes[i - 1], sub_ranks[i - 1], &sub_comms[i]);

    MPI_Comm_size(sub_comms[i], &sub_sizes[i]);
    MPI_Comm_rank(sub_comms[i], &sub_ranks[i]);
  }

  return 0;
}


slint_t mpi_subgroups_delete(slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_subgroups_delete */
{
  slint_t i;

  for (i = 0; i < nsubgroups; ++i) MPI_Comm_free(&sub_comms[i]);

  return 0;
}


/* sl_macro MC_REDUCEBCAST_THRESHOLD */
/* sl_macro MC_REDUCEBCAST_ROOT */

#if !defined(MC_REDUCEBCAST_THRESHOLD) && defined(GLOBAL_REDUCEBCAST_THRESHOLD)
# define MC_REDUCEBCAST_THRESHOLD  GLOBAL_REDUCEBCAST_THRESHOLD
#endif

#ifndef MC_REDUCEBCAST_ROOT
# define MC_REDUCEBCAST_ROOT  0
#endif

int sl_MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank) /* sl_proto, sl_func sl_MPI_Allreduce */
{
#ifdef MC_REDUCEBCAST_THRESHOLD
  if (size >= MC_REDUCEBCAST_THRESHOLD)
  {
    MPI_Reduce(sendbuf, recvbuf, count, datatype, op, MC_REDUCEBCAST_ROOT, comm);
    MPI_Bcast(recvbuf, count, datatype, MC_REDUCEBCAST_ROOT, comm);
    return 0;
  }
#endif
  return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}


/* sl_macro MC_ALLTOALL_INT_2STEP_THRESHOLD */

#if !defined(MC_ALLTOALL_INT_2STEP_THRESHOLD) && defined(GLOBAL_ALLTOALL_INT_2STEP_THRESHOLD)
# define MC_ALLTOALL_INT_2STEP_THRESHOLD  GLOBAL_ALLTOALL_INT_2STEP_THRESHOLD
#endif

#ifdef HAVE_ZMPI_TOOLS_H
# include "zmpi_tools.h"
#endif

int sl_MPI_Alltoall_int(void *sendbuf, int sendcount, void *recvbuf, int recvcount, MPI_Comm comm, int size, int rank) /* sl_proto, sl_func sl_MPI_Alltoall_int */
{
#if defined(HAVE_ZMPI_ALLTOALL_2STEP) && defined(MC_ALLTOALL_INT_2STEP_THRESHOLD)
  if (size >= MC_ALLTOALL_INT_2STEP_THRESHOLD)
  {
    return ZMPI_Alltoall_2step_int(sendbuf, sendcount, MPI_INT, recvbuf, recvcount, MPI_INT, comm);
  }
#endif
  return MPI_Alltoall(sendbuf, sendcount, MPI_INT, recvbuf, recvcount, MPI_INT, comm);
}



/* sl_macro ME_TRACE_IF */

#include "sl_common.h"


/* sl_ifdef SL_USE_MPI sl_context CONTEXT_BEGIN me */
const void *default_me_sendrecv_replace_mem = NULL;          /* sl_global sl_context sl_var default_me_sendrecv_replace_mem */
const slint_t default_me_sendrecv_replace_memsize = 0;       /* sl_global sl_context sl_var default_me_sendrecv_replace_memsize */
const slint_t default_me_sendrecv_replace_mpi_maxsize = -1;  /* sl_global sl_context sl_var default_me_sendrecv_replace_mpi_maxsize */
/* sl_endif sl_context CONTEXT_END me */

#ifndef ME_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define ME_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define ME_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_elements_keys_init_from_file(elements_t *s, char *filename, slint from, slint to, slint const_bytes_per_line, slint root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_keys_init_from_file */
{
  slint_t i, m, n, p, q;

  elements_t e;
  
  int send_counts[size], send_dipls[size];

  n = to - from + 1;

  if (root < 0 || size == 1)
  {
    m = (n / size) + ((n % size) > rank);
    p = ((int) (n / size)) * rank + z_min(rank, n % size);

    q = elements_keys_init_from_file(s, 1, filename, p, p + m - 1, const_bytes_per_line);

  } else
  {
    if (rank == root)
    {
      elements_keys_init_from_file(&e, 0, filename, from, to, const_bytes_per_line);

      for (i = 0; i < size; i++)
      {
        send_counts[i] = (n / size) + ((n % size) > i);
        send_counts[i] *= key_size_mpi;

        send_dipls[i] = ((int) (n / size)) * i + z_min(i, n % size);
        send_dipls[i] *= key_size_mpi;
      }
    }

    m = (n / size) + ((n % size) > rank);

    elements_alloc(s, m, SLCM_ALL);

    MPI_Scatterv(e.keys, send_counts, send_dipls, key_type_mpi, s->keys, m, key_type_mpi, root, comm);
    
    q = m;
    
    if (rank == root) elements_free(&e);
  }

  return q;
}


slint mpi_elements_validate_order(elements_t *s, slint n, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_validate_order */
{
  slint local_result = 0, global_result;
  slkey_pure_t pure_keys[2];
  MPI_Status status;

  if (size < 0) MPI_Comm_size(comm, &size);
  if (rank < 0) MPI_Comm_rank(comm, &rank);

  if (size > 1)
  {
    /* send lowest key to the left neighbor */
    pure_keys[0] = key_purify(s[0].keys[0]);

    if (rank == 0) MPI_Recv(&pure_keys[1], 1, pkey_mpi_datatype, rank + 1, 0, comm, &status);
    else if (rank + 1 == size) MPI_Send(&pure_keys[0], 1, pkey_mpi_datatype, rank - 1, 0, comm);
    else MPI_Sendrecv(&pure_keys[0], 1, pkey_mpi_datatype, rank - 1, 0, &pure_keys[1], 1, pkey_mpi_datatype, rank + 1, 0, comm, &status);

    if (rank + 1 < size) local_result += (key_pure_cmp_gt(key_purify(s[n - 1].keys[s[n - 1].size - 1]), pure_keys[1]) != 0);

    /* send highest key to the right neighbor */
    pure_keys[0] = key_purify(s[n - 1].keys[s[n - 1].size - 1]);

    if (rank == 0) MPI_Send(&pure_keys[0], 1, pkey_mpi_datatype, rank + 1, 0, comm);
    else if (rank + 1 == size) MPI_Recv(&pure_keys[1], 1, pkey_mpi_datatype, rank - 1, 0, comm, &status);
    else MPI_Sendrecv(&pure_keys[0], 1, pkey_mpi_datatype, rank + 1, 0, &pure_keys[1], 1, pkey_mpi_datatype, rank - 1, 0, comm, &status);

    if (rank > 0) local_result += (key_pure_cmp_lt(key_purify(s[0].keys[0]), pure_keys[1]) != 0);

    /* reduce the local results of the validation between neighbor-processes */
    MPI_Allreduce(&local_result, &global_result, 1, int_mpi_datatype, MPI_MAX, comm);
    /* exit if the validation fails */
    if (global_result > 0) return 1;
  }

  /* start the process-internal validation */
  local_result = (elements_validate_order(s, n) != 0);
  /* reduce the local results */
  MPI_Allreduce(&local_result, &global_result, 1, int_mpi_datatype, MPI_MAX, comm);

  return global_result;
}


slint_t mpi_linear_exchange_pure_keys(slkey_pure_t *in, slkey_pure_t *out, int size, int rank, MPI_Comm comm)  /* sl_proto, sl_func mpi_linear_exchange_pure_keys */
{
  MPI_Status status;

  out[0] = in[0];
  out[1] = in[1];

  /* exchange to the left */
  MPI_Sendrecv(&in[0], 1, pkey_mpi_datatype, (rank - 1 >= 0)?(rank - 1):MPI_PROC_NULL, 0,
               &out[1], 1, pkey_mpi_datatype, (rank + 1 < size)?(rank + 1):MPI_PROC_NULL, 0,
               comm, &status);

  /* exchange to the right */
  MPI_Sendrecv(&in[1], 1, pkey_mpi_datatype, (rank + 1 < size)?(rank + 1):MPI_PROC_NULL, 0,
               &out[0], 1, pkey_mpi_datatype, (rank - 1 >= 0)?(rank - 1):MPI_PROC_NULL, 0,
               comm, &status);

  if (rank - 1 < 0) out[0] = in[0];
  if (rank + 1 >= size) out[1] = in[1];

  return 0;
}


slint_t mpi_elements_check_order(elements_t *s, slint_t nelements, slint_t *orders, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_check_order */
{
  slint_t _orders[2], my_orders[2];
  slkey_pure_t my_keys[2], other_keys[2];

  if (orders == NULL) orders = _orders;

  /* check local order */
  my_orders[0] = (elements_validate_order(s, nelements) == 0);

  /* check global order (independent of local order) */
  elements_get_minmax_keys(s, nelements, my_keys);
  mpi_linear_exchange_pure_keys(my_keys, other_keys, size, rank, comm);
  my_orders[1] = (key_pure_cmp_le(other_keys[0], my_keys[0]) && key_pure_cmp_le(my_keys[1], other_keys[1]));

  MPI_Allreduce(my_orders, orders, 2, int_mpi_datatype, MPI_LAND, comm);

  return (orders[0] && orders[1]);
}


slint_t mpi_check_global_order(slkey_pure_t local_min, slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_check_global_order */
{
  slint_t local_order, global_order;
  slkey_pure_t my_keys[2], other_keys[2];


  my_keys[0] = local_min;
  my_keys[1] = local_max;
  
  mpi_linear_exchange_pure_keys(my_keys, other_keys, size, rank, comm);

  local_order = (key_pure_cmp_le(other_keys[0], my_keys[0]) && key_pure_cmp_le(my_keys[1], other_keys[1]));

  global_order = -1;
  
  if (root < 0) MPI_Allreduce(&local_order, &global_order, 1, int_mpi_datatype, MPI_LAND, comm);
  else MPI_Reduce(&local_order, &global_order, 1, int_mpi_datatype, MPI_LAND, root, comm);
  
  return global_order;
}


slint_t mpi_elements_digest_sum(elements_t *s, slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_digest_sum */
{
  unsigned int lsum;

#ifdef ELEMENTS_DIGEST_SUM_ADD
# define MPI_SUM_OP  MPI_SUM
#else
# define MPI_SUM_OP  MPI_BXOR
#endif

  elements_digest_sum(s, nelements, components, &lsum);

  MPI_Allreduce(&lsum, sum, 1, MPI_UNSIGNED, MPI_SUM_OP, comm);

  return 0;

#undef MPI_SUM_OP
}


/* deprecated */
unsigned int mpi_elements_crc32(elements_t *s, slint_t n, slint_t keys, slint_t data, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_crc32 */
{
  unsigned int crc32;

  mpi_elements_digest_sum(s, n, (keys?(SLCM_KEYS):0)|(data?(SLCM_DATA):0), &crc32, size, rank, comm);

  return crc32;
}


slint_t mpi_elements_digest_hash(elements_t *s, slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_digest_hash */
{
  slint_t hash_size = elements_digest_hash(NULL, 0, 0, NULL);
#ifndef HAVE_MPI_IN_PLACE
  void *hash_tmp;
#endif

  if (!hash) return hash_size;

#ifdef HAVE_MPI_IN_PLACE

  elements_digest_hash(s, nelements, components, hash);

  MPI_Allreduce(MPI_IN_PLACE, hash, hash_size, MPI_BYTE, MPI_BXOR, comm);

#else

  hash_tmp = z_alloca(hash_size, 1);

  elements_digest_hash(s, nelements, components, hash_tmp);

  MPI_Allreduce(hash_tmp, hash, hash_size, MPI_BYTE, MPI_BXOR, comm);

  z_freea(hash_tmp);

#endif

  return 0;
}


slint_t mpi_elements_get_counts(elements_t *s, slint_t *clocal, slint_t *cglobal, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_get_counts */
{
  slint_t lc, gc;


  if (!clocal) clocal = &lc;
  if (!cglobal) cglobal = &gc;
 
  *clocal = s->size;

  if (root < 0) MPI_Allreduce(clocal, cglobal, 1, int_mpi_datatype, MPI_SUM, comm);
  else MPI_Reduce(clocal, cglobal, 1, int_mpi_datatype, MPI_SUM, root, comm);

  return *cglobal;
}


slweight_t mpi_elements_get_weights(elements_t *s, slweight_t *wlocal, slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_get_weights */
{
  slweight_t lw, gw;


  if (!wlocal) wlocal = &lw;
  if (!wglobal) wglobal = &gw;

#ifdef elem_weight
  *wlocal = elements_get_weight(s);

  if (root < 0) MPI_Allreduce(wlocal, wglobal, 1, weight_mpi_datatype, MPI_SUM, comm);
  else MPI_Reduce(wlocal, wglobal, 1, weight_mpi_datatype, MPI_SUM, root, comm);
#else
  *wlocal = *wglobal = 0.0;
#endif
  return *wglobal;
}


slint_t mpi_elements_get_counts_and_weights(elements_t *s, slint_t nelements, slint_t *counts, slweight_t *weights, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_get_counts_and_weights */
{
  slint_t i;
  slweight_t lcw[2], gcw[2];


  lcw[0] = lcw[1] = gcw[0] = gcw[1] = 0.0;

  for (i = 0; i < nelements; ++i)
  {
    lcw[0] += s[i].size;
#ifdef elem_weight
    lcw[1] += elements_get_weight(&s[i]);
#endif
  }

  if (root < 0)
  {
    MPI_Allreduce(lcw, gcw, elem_weight_ifelse(2, 1), weight_mpi_datatype, MPI_SUM, comm);

  } else
  {
    MPI_Reduce(lcw, gcw, elem_weight_ifelse(2, 1), weight_mpi_datatype, MPI_SUM, root, comm);
  }
  
  counts[0] = lcw[0];
  counts[1] = gcw[0];
#ifdef elem_weight
  weights[0] = lcw[1];
  weights[1] = gcw[1];
#endif
    
  return 0;
}


slint_t mpi_elements_sendrecv_replace(elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm)  /* sl_proto, sl_func mpi_elements_sendrecv_replace */
{
  MPI_Status status;
  elements_t xs;
  int offset, current, maxc;

  slint_t replace, min_byte;
  slint_t smaxsize, rmaxsize;


  min_byte = key_byte;
#define xelem_call \
  if (min_byte < xelem_byte) min_byte = xelem_byte;
#include "sl_xelem_call.h"

  Z_TRACE_IF(ME_TRACE_IF, "me_sendrecv_replace_memsize: %" slint_fmt " vs. min_byte: %" slint_fmt "", SL_DEFCON(me.sendrecv_replace_memsize), min_byte);

  if (SL_DEFCON(me.sendrecv_replace_mem) == NULL || SL_DEFCON(me.sendrecv_replace_memsize) < min_byte)
  {
    replace = 1;
    rmaxsize = SL_DEFCON(me.sendrecv_replace_mpi_maxsize);

  } else
  {
    replace = 0;
    rmaxsize = SL_DEFCON(me.sendrecv_replace_memsize);
  }

  MPI_Sendrecv(&rmaxsize, 1, int_mpi_datatype, dest, sendtag, &smaxsize, 1, int_mpi_datatype, source, recvtag, comm, &status);

  Z_ASSERT(smaxsize != 0 && rmaxsize != 0);

  /* we make sure that total MPI message sizes fit into "int" (halfed), this "workaround" was (at least) necessary for JUROPA */

#define xelem_call \
  maxc = INT_MAX / 2 / xelem_byte; \
  if (smaxsize > 0) maxc = z_min(maxc, smaxsize / xelem_byte); \
  if (rmaxsize > 0) maxc = z_min(maxc, rmaxsize / xelem_byte); \
  Z_TRACE_IF(ME_TRACE_IF, xelem_name_str ": count: %d, maxc: %d", count, maxc); \
\
  offset = 0; \
  xelem_buf(&xs) = SL_DEFCON(me.sendrecv_replace_mem); \
\
  while (offset < count) \
  { \
    current = z_min(count - offset, maxc); \
    Z_TRACE_IF(ME_TRACE_IF, "current: %d", current); \
    if (replace) \
    { \
      Z_TRACE_IF(ME_TRACE_IF, "sendrecv_replace: %d @ %d to %d / from %d", current, offset, dest, source); \
      MPI_Sendrecv_replace(xelem_buf_at(s, offset), current, xelem_mpi_datatype, dest, sendtag, source, recvtag, comm, &status); \
    } else \
    { \
      Z_TRACE_IF(ME_TRACE_IF, "sendrecv: %d @ %d to %d from %d", current, offset, dest, source); \
      MPI_Sendrecv(xelem_buf_at(s, offset), current, xelem_mpi_datatype, dest, sendtag, xelem_buf(&xs), current, xelem_mpi_datatype, source, recvtag, comm, &status); \
      Z_TRACE_IF(ME_TRACE_IF, "ncopy: %d from 0 to %d", current, offset); \
      xelem_ncopy_at(&xs, 0, s, offset, current); \
    } \
    offset += current; \
  }
#include "sl_xelem_call.h"

  return 0;
}



/* sl_macro MEAS_TRACE_IF */

#include "sl_common.h"

#include "spec_core.h"


/* sl_ifdef SL_USE_MPI sl_context CONTEXT_BEGIN meas */
const double default_meas_t[2] = { 0 };     /* sl_global sl_context sl_var default_meas_t */
const slint_t default_meas_max_nprocs = 8;  /* sl_global sl_context sl_var default_meas_max_nprocs */
const slint_t default_meas_packed = 0;      /* sl_global sl_context sl_var default_meas_packed */
const slint_t default_meas_minalloc = 0;    /* sl_global sl_context sl_var default_meas_minalloc */
const double default_meas_overalloc = 0;    /* sl_global sl_context sl_var default_meas_overalloc */
/* sl_endif sl_context CONTEXT_END meas */


#ifndef MEAS_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MEAS_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MEAS_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


struct _tproc_t
{
  spec_tproc_t spec_tproc;
};


static void _tproc_create(tproc_t *tproc)
{
  *tproc = z_alloc(1, sizeof(*tproc));
}


static void _tproc_free(tproc_t *tproc)
{
  z_free(*tproc);
  
  *tproc = NULL;
}


slint_t tproc_create_tproc(tproc_t *tproc, tproc_f *tfn, tproc_reset_f *rfn, tproc_exdef exdef) /* sl_proto, sl_func tproc_create_tproc */
{
  if (exdef != TPROC_EXDEF_NULL && exdef->type != 1) return 1;

  _tproc_create(tproc);

  spec_tproc_create(&(*tproc)->spec_tproc, tfn, NULL, NULL, NULL);

  spec_tproc_set_reset(&(*tproc)->spec_tproc, rfn);

  if (exdef != TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tproc(&(*tproc)->spec_tproc, exdef->tproc_count_db, exdef->tproc_count_ip, exdef->tproc_rearrange_db, exdef->tproc_rearrange_ip);

  return 0;
}


slint_t tproc_create_tproc_mod(tproc_t *tproc, tproc_mod_f *tfn, tproc_reset_f *rfn, tproc_exdef exdef) /* sl_proto, sl_func tproc_create_tproc_mod */
{
  if (exdef != TPROC_EXDEF_NULL && exdef->type != 2) return 1;

  _tproc_create(tproc);

  spec_tproc_create(&(*tproc)->spec_tproc, NULL, tfn, NULL, NULL);

  spec_tproc_set_reset(&(*tproc)->spec_tproc, rfn);

  if (exdef != TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tproc_mod(&(*tproc)->spec_tproc, exdef->tproc_mod_count_db, exdef->tproc_mod_count_ip, exdef->tproc_mod_rearrange_db, exdef->tproc_mod_rearrange_ip);

  return 0;
}


slint_t tproc_create_tprocs(tproc_t *tproc, tprocs_f *tfn, tproc_reset_f *rfn, tproc_exdef exdef) /* sl_proto, sl_func tproc_create_tprocs */
{
  if (exdef != TPROC_EXDEF_NULL && exdef->type != 3) return 1;

  _tproc_create(tproc);
  
  spec_tproc_create(&(*tproc)->spec_tproc, NULL, NULL, tfn, NULL);

  spec_tproc_set_reset(&(*tproc)->spec_tproc, rfn);

  if (exdef != TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tprocs(&(*tproc)->spec_tproc, exdef->tprocs_count_db, exdef->tprocs_count_ip, exdef->tprocs_rearrange_db, exdef->tprocs_rearrange_ip);
  
  return 0;
}


slint_t tproc_create_tprocs_mod(tproc_t *tproc, tprocs_mod_f *tfn, tproc_reset_f *rfn, tproc_exdef exdef) /* sl_proto, sl_func tproc_create_tprocs_mod */
{
  if (exdef != TPROC_EXDEF_NULL && exdef->type != 4) return 1;

  _tproc_create(tproc);
  
  spec_tproc_create(&(*tproc)->spec_tproc, NULL, NULL, NULL, tfn);

  spec_tproc_set_reset(&(*tproc)->spec_tproc, rfn);
  
  if (exdef != TPROC_EXDEF_NULL)
    spec_tproc_set_ext_tprocs_mod(&(*tproc)->spec_tproc, exdef->tprocs_mod_count_db, exdef->tprocs_mod_count_ip, exdef->tprocs_mod_rearrange_db, exdef->tprocs_mod_rearrange_ip);
  
  return 0;
}


slint_t tproc_free(tproc_t *tproc) /* sl_proto, sl_func tproc_free */
{
  spec_tproc_destroy(&(*tproc)->spec_tproc);

  _tproc_free(tproc);
  
  return 0;
}


slint_t tproc_set_proclists(tproc_t *tproc, slint_t nsend_procs, int *send_procs, slint_t nrecv_procs, int *recv_procs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func tproc_set_proclists */
{
  spec_tproc_set_proclists(&(*tproc)->spec_tproc, nsend_procs, send_procs, nrecv_procs, recv_procs, size, rank, comm);

  return 0;
}


slint_t tproc_verify(tproc_t tproc, void *data, elements_t *s, int proc) /* sl_proto, sl_func tproc_verify */
{
  /* FIXME */

/*  slint_t i, j, n;
  int procs[SL_DEFCON(meas.max_nprocs)];


  for (i = 0; i < s->size; ++i)
  {
    n = tproc->tprocs_mod(s, i, data, procs, NULL);
    for (j = 0; j < n; ++j)
    if (procs[j] == proc) n = 0;

    if (n > 0) return 0;
  }*/

  return 1;
}


slint_t mpi_elements_alltoall_specific_alloc_size(slint_t n) /* sl_func mpi_elements_alltoall_specific_alloc_size */
{
  if (SL_DEFCON(meas.overalloc) < 0) n = (slint_t) (n * (1.0 - SL_DEFCON(meas.overalloc)));
  else n += SL_DEFCON(meas.overalloc);

  return z_max(n, SL_DEFCON(meas.minalloc));
}


slint_t mpi_elements_alltoall_specific(elements_t *sin, elements_t *sout, elements_t *xs, tproc_t tproc, void *data, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoall_specific */
{
  slint_t r, original_mea_packed;

  original_mea_packed = SL_DEFCON(mea.packed); SL_DEFCON(mea.packed) = SL_DEFCON(meas.packed);

  if (sin == NULL || sout == NULL)
    r = spec_alltoallv_ip((sin)?sin:sout, xs, tproc->spec_tproc, data, size, rank, comm);
  else
    r = spec_alltoallv_db(sin, sout, xs, tproc->spec_tproc, data, size, rank, comm);

  SL_DEFCON(mea.packed) = original_mea_packed;

  return r;
}



/* sl_macro MEA_TRACE_IF */
/* sl_macro MEA_INPLACE_ALLTOALLV */

#include "sl_common.h"

#ifdef HAVE_ALLTOALLV_H
# include "alltoallv.h"
#endif

#include "dash_core.h"
#include "dash_sched_a2av.h"
#include "dash_sched_a2av_aux.h"
#include "dash_exec_sl.h"


/*#define MEA_INPLACE_ALLTOALLV*/

/* sl_ifdef SL_USE_MPI sl_context CONTEXT_BEGIN mea */
const slint_t default_mea_packed = 0;     /* sl_global sl_context sl_var default_mea_packed */
const slint_t default_mea_db_packed = 0;  /* sl_global sl_context sl_var default_mea_db_packed */
const slint_t default_mea_ip_packed = 0;  /* sl_global sl_context sl_var default_mea_ip_packed */
/* sl_endif sl_context CONTEXT_END mea */

#ifndef MEA_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MEA_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MEA_TRACE_IF  (SL_PROC_RANK == -1)
# endif
#endif


static slint_t _mpi_elements_alltoallv_db_packed(elements_t *sbuf, int *scounts, int *sdispls, elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm, slint_t *lbs, slint_t *extents)
{
  slint_t i, _lbs[2] = { 0, 0 }, _extents[2] = { -1, -1 };
  slint_t local_exit, global_exit;
  packed_elements_t spackbuf, rpackbuf;
  MPI_Datatype packed_type;


  if (lbs == NULL) lbs = _lbs;
  if (extents == NULL) extents = _extents;

  if (extents[0] < 0) get_displcounts_extent(size, sdispls, scounts, &lbs[0], &extents[0]);
  if (extents[1] < 0) get_displcounts_extent(size, rdispls, rcounts, &lbs[1], &extents[1]);

  if (elements_alloc_packed(&spackbuf, extents[0]) != 0 || elements_alloc_packed(&rpackbuf, extents[1]) != 0)
  {
    elements_free_packed(&spackbuf);
    elements_free_packed(&rpackbuf);
    local_exit = 1;
  
  } else local_exit = 0;

  MPI_Allreduce(&local_exit, &global_exit, 1, int_mpi_datatype, MPI_SUM, comm);

  Z_TRACE_IF(MEA_TRACE_IF, "local_exit: %" slint_fmt ", global_exit: %" slint_fmt, local_exit, global_exit);

  if (global_exit) return -1;

  pelem_add(&spackbuf, -lbs[0]);
  pelem_add(&rpackbuf, -lbs[1]);

  mpi_elements_packed_datatype_create(&packed_type, 0);

  for (i = 0; i < size; ++i) elem_npack_at(sbuf, sdispls[i], &spackbuf, sdispls[i], scounts[i]);

  MPI_Alltoallv(spackbuf.elements, scounts, sdispls, packed_type, rpackbuf.elements, rcounts, rdispls, packed_type, comm);

  for (i = 0; i < size; ++i) pelem_nunpack_at(&rpackbuf, rdispls[i], rbuf, rdispls[i], rcounts[i]);

  mpi_elements_packed_datatype_destroy(&packed_type);

  pelem_add(&spackbuf, lbs[0]);
  pelem_add(&rpackbuf, lbs[1]);

  elements_free_packed(&spackbuf);
  elements_free_packed(&rpackbuf);

  return 0;
}


slint_t mpi_elements_alltoallv_db_packed(elements_t *sbuf, int *scounts, int *sdispls, elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_db_packed */
{
  return _mpi_elements_alltoallv_db_packed(sbuf, scounts, sdispls, rbuf, rcounts, rdispls, size, rank, comm, NULL, NULL);
}


slint_t mpi_elements_alltoallv_db(elements_t *sbuf, int *scounts, int *sdispls, elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_db */
{
  slint_t lbs[2] = { 0, 0 }, extents[2] = { -1, -1 };
  slint_t r = -1;

  Z_TRACE_IF(MEA_TRACE_IF, "packed: %" slint_fmt " or %" slint_fmt, SL_DEFCON(mea.packed), SL_DEFCON(mea.db_packed));

  if (SL_DEFCON(mea.packed) || SL_DEFCON(mea.db_packed))
    r = _mpi_elements_alltoallv_db_packed(sbuf, scounts, sdispls, rbuf, rcounts, rdispls, size, rank, comm, lbs, extents);

  if (r == 0) return 0;

#define xelem_call \
  MPI_Alltoallv(xelem_buf(sbuf), scounts, sdispls, xelem_mpi_datatype, xelem_buf(rbuf), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  return 0;
}


static slint_t _mpi_elements_alltoallv_ip_packed(elements_t *s, elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm, slint_t *lbs, slint_t *extents)
{
  slint_t i, _lbs[2] = { 0, 0 }, _extents[2] = { -1, -1 };
  slint_t local_exit, global_exit;
  void *block;
  slint_t blocksize;
  packed_elements_t packbuf, spackbuf, rpackbuf;
  MPI_Datatype packed_type;


  if (lbs == NULL) lbs = _lbs;
  if (extents == NULL) extents = _extents;

  if (extents[0] < 0) get_displcounts_extent(size, sdispls, scounts, &lbs[0], &extents[0]);
  if (extents[1] < 0) get_displcounts_extent(size, rdispls, rcounts, &lbs[1], &extents[1]);

  Z_TRACE_IF(MEA_TRACE_IF, "send: lb: %" slint_fmt ", extent: %" slint_fmt, lbs[0], extents[0]);
  Z_TRACE_IF(MEA_TRACE_IF, "receive: lb: %" slint_fmt ", extent: %" slint_fmt, lbs[1], extents[1]);

  if (sx)
  {
    elements_alloc_block(sx, &block, &blocksize, 8, -1);

    elements_alloc_packed_from_block(&packbuf, block, blocksize, 8, -1);

    Z_TRACE_IF(MEA_TRACE_IF, "alloc block size: %" slint_fmt ", alloc packed elements: %" slint_fmt, blocksize, packbuf.size);

    local_exit = (packbuf.size < extents[0] + extents[1])?1:0;

  } else local_exit = 1;

  MPI_Allreduce(&local_exit, &global_exit, 1, int_mpi_datatype, MPI_SUM, comm);

  Z_TRACE_IF(MEA_TRACE_IF, "local_exit: %" slint_fmt ", global_exit: %" slint_fmt, local_exit, global_exit);

  if (global_exit) return -1;

  pelem_assign(&packbuf, &spackbuf);
  pelem_assign_at(&packbuf, extents[0], &rpackbuf);

  mpi_elements_packed_datatype_create(&packed_type, 0);

  pelem_add(&spackbuf, -lbs[0]);
  pelem_add(&rpackbuf, -lbs[1]);

  for (i = 0; i < size; ++i) elem_npack_at(s, sdispls[i], &spackbuf, sdispls[i], scounts[i]);

  MPI_Alltoallv(spackbuf.elements, scounts, sdispls, packed_type, rpackbuf.elements, rcounts, rdispls, packed_type, comm);

  for (i = 0; i < size; ++i) pelem_nunpack_at(&rpackbuf, rdispls[i], s, rdispls[i], rcounts[i]);

  mpi_elements_packed_datatype_destroy(&packed_type);

  return 0;
}


slint_t mpi_elements_alltoallv_ip_packed(elements_t *s, elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_ip_packed */
{
  return _mpi_elements_alltoallv_ip_packed(s, sx, scounts, sdispls, rcounts, rdispls, size, rank, comm, NULL, NULL);
}


static slint_t _mpi_elements_alltoallv_ip_double(elements_t *s, elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm, slint_t *lbs, slint_t *extents)
{
  slint_t i, _lbs[2] = { 0, 0 }, _extents[2] = { -1, -1 };
  slint_t local_exit, global_exit;
  elements_t _sx, _bx;
  void *block;
  slint_t blocksize;


  if (lbs == NULL) lbs = _lbs;
  if (extents == NULL) extents = _extents;

  if (extents[1] < 0) get_displcounts_extent(size, rdispls, rcounts, &lbs[1], &extents[1]);

  Z_TRACE_IF(MEA_TRACE_IF, "receive: lb: %" slint_fmt ", extent: %" slint_fmt, lbs[1], extents[1]);

  if (sx)
  {
    elements_alloc_block(sx, &block, &blocksize, 8, -1);

    Z_TRACE_IF(MEA_TRACE_IF, "alloc block size: %" slint_fmt, blocksize);

    if (blocksize >= extents[1] * elem_get_max_byte())
    {
      elem_set_block(&_bx, block);
      elem_set_block_size(&_bx, blocksize);
      sx = &_bx;
    }
      
    if (elem_get_block(sx))
    {
      Z_TRACE_IF(MEA_TRACE_IF, "block mem: size: %" slint_fmt ", max byte per component: %" slint_fmt, elem_get_block_size(sx), elem_get_max_byte());

      local_exit = (elem_get_block_size(sx) < extents[1] * elem_get_max_byte())?1:0;

    } else
    {
      Z_TRACE_IF(MEA_TRACE_IF, "elements mem: size: %" slint_fmt, elem_get_max_size(sx));

      local_exit = (elem_get_max_size(sx) < extents[1])?1:0;
    }

  } else local_exit = 1;


  MPI_Allreduce(&local_exit, &global_exit, 1, int_mpi_datatype, MPI_SUM, comm);

  Z_TRACE_IF(MEA_TRACE_IF, "local_exit: %" slint_fmt ", global_exit: %" slint_fmt, local_exit, global_exit);

  if (global_exit) return -1;
  
  if (elem_get_block(sx))
  {
#define xelem_call \
  xelem_set(&_sx, elem_get_block(sx)); \
  xelem_add(&_sx, -lbs[1]); \
  MPI_Alltoallv(xelem_buf(s), scounts, sdispls, xelem_mpi_datatype, xelem_get(&_sx), rcounts, rdispls, xelem_mpi_datatype, comm); \
  for (i = 0; i < size; ++i) xelem_ncopy_at(&_sx, rdispls[i], s, rdispls[i], rcounts[i]);
#include "sl_xelem_call.h"

  } else
  {
    elem_assign_at(sx, -lbs[1], &_sx);

#define xelem_call \
  MPI_Alltoallv(xelem_buf(s), scounts, sdispls, xelem_mpi_datatype, xelem_buf(&_sx), rcounts, rdispls, xelem_mpi_datatype, comm); \
  for (i = 0; i < size; ++i) xelem_ncopy_at(&_sx, rdispls[i], s, rdispls[i], rcounts[i]);
#include "sl_xelem_call.h"
  }
  
  return 0;
}


slint_t mpi_elements_alltoallv_ip_double(elements_t *s, elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_ip_double */
{
  return _mpi_elements_alltoallv_ip_double(s, sx, scounts, sdispls, rcounts, rdispls, size, rank, comm, NULL, NULL);
}


slint_t mpi_elements_alltoallv_ip_mpi(elements_t *s, elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_ip_mpi */
{
#if defined(HAVE_ALLTOALLV_H)

#define xelem_call \
  alltoallv_inplace_ds_aux = xelem_buf(sx); \
  alltoallv_inplace_ds_aux_size = sx->size * xelem_byte; \
  alltoallv_inplace_ds(MPI_IN_PLACE, scounts, sdispls, xelem_mpi_datatype, s, rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  alltoallv_inplace_ds_aux = NULL;
  alltoallv_inplace_ds_aux_size = 0;

  return 0;

#else

  return -1;

#endif
}


slint_t mpi_elements_alltoallv_ip_dash(elements_t *s, elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_ip_dash */
{
  ds_t ds;
  ds_sched_t sched;
  ds_sched_a2av_aux_t sched_a2av_aux;
  ds_exec_t exec;

  slint_t aux_size = 1000;  
  slint_t aux_blocks = -8;

  elements_t _sx, _bx;


  if (sx == NULL)
  {
    sx = &_sx;
    elements_alloc(sx, aux_size, SLCM_ALL);

  } else if (elem_get_block(sx))
  {
    sx = &_bx;
    elements_alloc_from_block(&_bx, elem_get_block(sx), elem_get_block_size(sx), 8, -1, SLCM_ALL);
  }

  ds_exec_sl_create(&exec);

  dsint_t s_addr_id = ds_exec_sl_add_address(&exec, s);
  dsint_t sx_addr_id = ds_exec_sl_add_address(&exec, sx);

  ds_sched_a2av_create(&sched);

  dsint_t s_buf_id = ds_sched_add_buffer(&sched, s_addr_id);
  ds_sched_set_send(&sched, s_buf_id, size, scounts, sdispls, NULL, 0);
  ds_sched_set_recv(&sched, s_buf_id, size, rcounts, rdispls, NULL, 0);

  dsint_t sx_buf_id = ds_sched_add_buffer(&sched, sx_addr_id);
  ds_sched_set_aux(&sched, sx_buf_id, sx->size);

  if (aux_blocks < 0) aux_blocks = size / -aux_blocks;
  aux_blocks = z_minmax(1, aux_blocks, size);

  ds_aux_static_create(&sched_a2av_aux, sx->size, aux_blocks);

  ds_sched_a2av_set_aux(&sched, &sched_a2av_aux);

  ds_create(&ds, &sched, &exec);
  
  ds.comm = comm;
  ds.comm_size = size;
  ds.comm_rank = rank;

  ds_run(&ds);

  ds_destroy(&ds);

  ds_aux_static_destroy(&sched_a2av_aux);

  ds_sched_a2av_destroy(&sched);

  ds_exec_sl_destroy(&exec);

  if (sx == &_sx) elements_free(sx);

  return 0;
}


slint_t mpi_elements_alltoallv_ip(elements_t *s, elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_ip */
{
  slint_t lbs[2] = { 0, 0 }, extents[2] = { -1, -1 };
  slint_t r = -1;

  Z_TRACE_IF(MEA_TRACE_IF, "packed: %" slint_fmt " or %" slint_fmt, SL_DEFCON(mea.packed), SL_DEFCON(mea.ip_packed));

  if (SL_DEFCON(mea.packed) || SL_DEFCON(mea.ip_packed))
    r = _mpi_elements_alltoallv_ip_packed(s, sx, scounts, sdispls, rcounts, rdispls, size, rank, comm, lbs, extents);

  if (r == 0) return 0;

  r = _mpi_elements_alltoallv_ip_double(s, sx, scounts, sdispls, rcounts, rdispls, size, rank, comm, lbs, extents);

  if (r == 0) return 0;

#if defined(MEA_INPLACE_ALLTOALLV) && defined(HAVE_ALLTOALLV_H)

  r = mpi_elements_alltoallv_ip_mpi(s, sx, scounts, sdispls, rcounts, rdispls, size, rank, comm);

#else

  r = mpi_elements_alltoallv_ip_dash(s, sx, scounts, sdispls, rcounts, rdispls, size, rank, comm);
  
#endif

  return r;
}



/* sl_macro MEAP_TRACE_IF */

#include "sl_common.h"

#include "zmpi_tools.h"


#ifndef MEAP_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MEAP_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MEAP_TRACE_IF  (SL_PROC_RANK == -1)
# endif
#endif


#ifdef HAVE_ZMPI_ALLTOALLV_PROCLISTS

slint_t mpi_elements_alltoallv_proclists_db(elements_t *sbuf, int *scounts, int *sdispls, int nsendprocs, int *sendprocs, elements_t *rbuf, int *rcounts, int *rdispls, int nrecvprocs, int *recvprocs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_proclists_db */
{
#define xelem_call \
  ZMPI_Alltoallv_proclists(xelem_buf(sbuf), scounts, sdispls, xelem_mpi_datatype, nsendprocs, sendprocs, xelem_buf(rbuf), rcounts, rdispls, xelem_mpi_datatype, nrecvprocs, recvprocs, comm);
#include "sl_xelem_call.h"

  return 0;
}

#endif



#include "sl_common.h"


slint_t mpi_elements_packed_datatype_create(MPI_Datatype *pdt, slint_t structured) /* sl_proto, sl_func mpi_elements_packed_datatype_create */
{
  packed_element_t pe;
  
  int blengths[elem_n];
  MPI_Aint displs[elem_n];
  MPI_Datatype types[elem_n];

  int i = 0;
  
  MPI_Aint base;

  if (structured)
  {
    MPI_Get_address(&pe, &base);

#define xelem_index_not
#define xelem_call \
  blengths[i] = xelem_size_mpi; \
  MPI_Get_address(&pe.xelem_name_packed, &displs[i]); displs[i] -= base; \
  types[i] = xelem_type_mpi; \
  i++;
#include "sl_xelem_call.h"

    MPI_Type_create_struct(1 + data_n, blengths, displs, types, pdt);

  } else MPI_Type_contiguous(pelem_byte, MPI_BYTE, pdt);

  MPI_Type_commit(pdt);

  return 0;
}


slint_t mpi_elements_packed_datatype_destroy(MPI_Datatype *pdt) /* sl_proto, sl_func mpi_elements_packed_datatype_destroy */
{
  MPI_Type_free(pdt);

  return 0;
}


/* - basierend auf "bitonic merge-exchange", gefunden in [taxonomy] mit Verweis auf [Alekseyev 1969] und [Knuth 1973]
   - binäre Suche durch [Tridgell,Brent] Übergang zur Verwendung von "bisection" beschrieben in [Zhou,Tridgell]
*/

#include "sl_common.h"


/*#define CHECK_EQUAL*/

#ifndef MFE_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MFE_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MFE_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_find_exact_equal(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *ex_start, slint_t *ex_size, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_find_exact_equal */
{
  const int tag = 1;

  elements_t interval;

  slint_t low_cmp_i, high_cmp_i, *my_cmp_i;
#ifdef CHECK_EQUAL
  slint_t *other_cmp_i;
#endif
  slkey_pure_t low_cmp_key, high_cmp_key, *my_cmp_key, *other_cmp_key;
  
  slint_t local_ex_start, local_ex_size;

  MPI_Status status;


  if (ex_start == NULL) ex_start = &local_ex_start;
  if (ex_size == NULL) ex_size = &local_ex_size;

  if (rank != high_rank)
  {
    my_cmp_i = &low_cmp_i;
#ifdef CHECK_EQUAL
    other_cmp_i = &high_cmp_i;
#endif

    my_cmp_key = &low_cmp_key;
    other_cmp_key = &high_cmp_key;

  } else
  {
    my_cmp_i = &high_cmp_i;
#ifdef CHECK_EQUAL
    other_cmp_i = &low_cmp_i;
#endif

    my_cmp_key = &high_cmp_key;
    other_cmp_key = &low_cmp_key;
  }

#ifdef CHECK_EQUAL
  *my_cmp_i = s->size;
  MPI_Sendrecv(my_cmp_i, 1, int_mpi_datatype, other_rank, tag, other_cmp_i, 1, int_mpi_datatype, other_rank, tag, comm, &status);
  if (*other_cmp_i != s->size)
  {
    fprintf(stderr, "mpi_find_exact_equal: error: element lists not equal (%" sl_int_type_fmt " vs. %" sl_int_type_fmt ")\n", *my_cmp_i, *other_cmp_i);
    return -1;
  }
#endif

  elem_assign(s, &interval);

  while (interval.size > 0)
  {
    low_cmp_i = interval.size / 2;
    high_cmp_i = interval.size - low_cmp_i - 1;
    
    *my_cmp_key = key_purify(*elem_key_at(&interval, *my_cmp_i));

    MPI_Sendrecv(my_cmp_key, key_pure_size_mpi, key_pure_type_mpi, other_rank, tag, other_cmp_key, key_pure_size_mpi, key_pure_type_mpi, other_rank, tag, comm, &status);

    if (key_pure_cmp_le(low_cmp_key, high_cmp_key))
    {
      if (rank != high_rank) elem_add(&interval, low_cmp_i + 1);
      
      interval.size = high_cmp_i;

    } else
    {
      if (rank == high_rank) elem_add(&interval, high_cmp_i + 1);
      
      interval.size = low_cmp_i;
    }
  }

  if (rank != high_rank)
  {
    *ex_start = interval.keys - s->keys;
    *ex_size = s->size - *ex_start;

  } else
  {
    *ex_start = 0;
    *ex_size = interval.keys - s->keys;
  }

  return *ex_size;
}


slint_t mpi_find_exact(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, slint_t *ex_start, slint_t *ex_sizes, slint_t *nx_move, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_find_exact */
{
  const int tag = 1;

  slint_t low_sizes[2], high_sizes[2], *my_sizes, *other_sizes;
  slint_t low_diff, high_diff;
  slint_t low_int[2], high_int[2], int_size;
  slint_t local_ex_start, local_ex_sizes[2], local_nx_move;
  
  elements_t interval;

  MPI_Status status;


  if (ex_start == NULL) ex_start = &local_ex_start;
  if (ex_sizes == NULL) ex_sizes = local_ex_sizes;
  if (nx_move == NULL) nx_move = &local_nx_move;

  if (rank != high_rank)
  {
    my_sizes = low_sizes;
    other_sizes = high_sizes;

  } else
  {
    my_sizes = high_sizes;
    other_sizes = low_sizes;
  }
  
  my_sizes[0] = s->size;
  if (dst_size) my_sizes[1] = *dst_size; else my_sizes[1] = s->size;

  MPI_Sendrecv(my_sizes, 2, int_mpi_datatype, other_rank, tag, other_sizes, 2, int_mpi_datatype, other_rank, tag, comm, &status);

  low_diff = low_sizes[1] - low_sizes[0];
  high_diff = high_sizes[1] - high_sizes[0];

  /* correct diffs */  
  if (low_diff < -low_sizes[0])
  {
    if (high_diff < -high_sizes[0])
    {
      low_diff = 0;
      high_diff = 0;

    } else
    {
      low_diff = -z_min(high_diff, low_sizes[0]);
      high_diff = -low_diff;
    }

  } else
  {
    if (high_diff < -high_sizes[0] || low_diff + high_diff >= 0)
    {
      high_diff = -z_min(low_diff, high_sizes[0]);
      low_diff = -high_diff;

    } else if (low_diff + high_diff < 0)
    {
      fprintf(stderr, "mpi_find_exact2: error: destination sizes too small (%" sl_int_type_fmt " + %" sl_int_type_fmt " < 0!)\n", low_diff, high_diff);
      return -1;
    }
  }

  low_sizes[1] = low_diff + low_sizes[0];
  high_sizes[1] = high_diff + high_sizes[0];

  /* FIXME: APP_ASSERT(low_diff == -high_diff) */
  
  low_int[0] = 0;
  low_int[1] = low_sizes[0] - 1;
  
  high_int[0] = 0;
  high_int[1] = high_sizes[0] - 1;
  
  if (low_diff > 0) high_int[0] += low_diff;
  if (high_diff > 0) low_int[1] -= high_diff;
  
  int_size = z_min(low_int[1] - low_int[0] + 1, high_int[1] - high_int[0] + 1);
  
  low_int[0] = low_int[1] - int_size + 1;
  high_int[1] = high_int[0] + int_size - 1;
  
  if (rank != high_rank) elem_assign_at(s, low_int[0], &interval);
  else elem_assign_at(s, high_int[0], &interval);
  
  interval.size = int_size;
  
  mpi_find_exact_equal(&interval, other_rank, high_rank, ex_start, ex_sizes, size, rank, comm);

  low_int[0] = low_int[1] - *ex_sizes + 1;
  high_int[1] = high_int[0] + *ex_sizes - 1;

  if (low_diff > 0) high_int[0] -= low_diff;
  if (high_diff > 0) low_int[1] += high_diff;

  if (rank != high_rank)
  {
    if (dst_size) *dst_size = low_sizes[1];

    *ex_start = low_int[0];

    ex_sizes[0] = low_int[1] - low_int[0] + 1;
    ex_sizes[1] = high_int[1] - high_int[0] + 1;

    *nx_move = 0;

  } else
  {
    if (dst_size) *dst_size = high_sizes[1];

    *ex_start = high_int[0];

    ex_sizes[0] = high_int[1] - high_int[0] + 1;
    ex_sizes[1] = low_int[1] - low_int[0] + 1;

    *nx_move = high_diff;
  }
  
  return 0;
}


#undef CHECK_EQUAL



/* sl_macro MM2_ELEMENTS_SENDRECV_REPLACE */
/* sl_macro MM2_TRACE_IF */
/* sl_macro MM2_PRINT_TIMINGS */

#define MM2_ELEMENTS_SENDRECV_REPLACE

#include "sl_common.h"


#ifdef SLDEBUG
# define CHECK_ORDER
#endif

#ifndef MM2_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MM2_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MM2_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_merge2(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, merge2x_f m2, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_merge2 */
{
  const int tag = 1;

  slint_t ex_start, ex_sizes[2], nx_move, ex_size;
  elements_t s0, s1;

  MPI_Status status;

#ifdef CHECK_ORDER
  slint_t check_order;
#endif


  Z_TRACE_IF(MM2_TRACE_IF, "starting mpi_merge2");

  /* sl_tid rti_tid_mpi_merge2 */

  rti_treset(rti_tid_mpi_merge2_find);       /* sl_tid */
  rti_treset(rti_tid_mpi_merge2_moveright);  /* sl_tid */
  rti_treset(rti_tid_mpi_merge2_exchange);   /* sl_tid */
  rti_treset(rti_tid_mpi_merge2_moveleft);   /* sl_tid */
  rti_treset(rti_tid_mpi_merge2_local);     /* sl_tid */

  rti_tclear(rti_tid_mpi_merge2);

  if (other_rank < 0 || other_rank >= size) return -1;

  if (rank == other_rank) return 0;

  rti_tstart(rti_tid_mpi_merge2);

#ifdef CHECK_ORDER
  check_order = elements_validate_order(s, 1);
  if (check_order) Z_ERROR("input order failed at %" slint_fmt "", check_order);
#endif

  Z_TRACE_IF(MM2_TRACE_IF, "find_exact: s->size = %" slint_fmt ", other_rank / high_rank = %" slint_fmt " / %" slint_fmt, s->size, other_rank, high_rank);

  rti_tstart(rti_tid_mpi_merge2_find);
  mpi_find_exact(s, other_rank, high_rank, dst_size, &ex_start, ex_sizes, &nx_move, size, rank, comm);
  rti_tstop(rti_tid_mpi_merge2_find);

  Z_TRACE_IF(MM2_TRACE_IF, "find_exact: ex_start = %" slint_fmt ", ex_sizes = { %" slint_fmt ", %" slint_fmt " }, nx_move = %" slint_fmt, ex_start, ex_sizes[0], ex_sizes[1], nx_move);

  /* move the nx-block to the right (before exchange) */
  rti_tstart(rti_tid_mpi_merge2_moveright);

  if (nx_move > 0 && s->size - ex_sizes[0] > 0)
  {
    Z_TRACE_IF(MM2_TRACE_IF, "moving right %" slint_fmt "", nx_move);

    if (rank != high_rank) elem_nmove_at(s, 0, s, nx_move, s->size - ex_sizes[0]);
    else elem_nmove_at(s, ex_sizes[0], s, ex_sizes[0] + nx_move, s->size - ex_sizes[0]);
  }

  rti_tstop(rti_tid_mpi_merge2_moveright);

  /* exchange elements */
  rti_tstart(rti_tid_mpi_merge2_exchange);

  elem_assign_at(s, ex_start, &s0);
  ex_size = z_min(ex_sizes[0], ex_sizes[1]);

  if (ex_size > 0)
  {
    Z_TRACE_IF(MM2_TRACE_IF, "exchanging %" slint_fmt " elements at %" slint_fmt "", ex_size, ex_start);

#ifdef MM2_ELEMENTS_SENDRECV_REPLACE
    mpi_elements_sendrecv_replace(&s0, ex_size, other_rank, tag, other_rank, tag, size, rank, comm);
#else
#define xelem_call \
    MPI_Sendrecv_replace(xelem_buf(&s0), ex_size, xelem_mpi_datatype, other_rank, tag, other_rank, tag, comm, &status);
#include "sl_xelem_call.h"
#endif
  }

  elem_add(&s0, ex_size);

  if (ex_size < ex_sizes[0])
  {
    ex_size = ex_sizes[0] - ex_size;
    
    Z_TRACE_IF(MM2_TRACE_IF, "sending %" slint_fmt " at %" slint_fmt "", ex_size, (slint_t) (s0.keys - s->keys));

#define xelem_call \
    MPI_Send(xelem_buf(&s0), ex_size, xelem_mpi_datatype, other_rank, tag, comm);
#include "sl_xelem_call.h"

  } else if (ex_size < ex_sizes[1])
  {
    ex_size = ex_sizes[1] - ex_size;

    Z_TRACE_IF(MM2_TRACE_IF, "receiving %" slint_fmt " at %" slint_fmt "", ex_size, (slint_t) (s0.keys - s->keys));

#define xelem_call \
    MPI_Recv(xelem_buf(&s0), ex_size, xelem_mpi_datatype, other_rank, tag, comm, &status);
#include "sl_xelem_call.h"
  }

  rti_tstop(rti_tid_mpi_merge2_exchange);

  /* move the nx-block to the left (after exchange) */
  rti_tstart(rti_tid_mpi_merge2_moveleft);

  if (nx_move < 0 && s->size - ex_sizes[0] > 0)
  {
    Z_TRACE_IF(MM2_TRACE_IF, "moving left %" slint_fmt "", nx_move);

    if (rank != high_rank) elem_nmove_at(s, 0, s, nx_move, s->size - ex_sizes[0]);
    else elem_nmove_at(s, ex_sizes[0], s, ex_sizes[0] + nx_move, s->size - ex_sizes[0]);
  }

  rti_tstop(rti_tid_mpi_merge2_moveleft);

  /* prepare the local merge2 */
  if (rank != high_rank)
  {
    elem_assign_at(s, 0, &s0);
    s0.size = s->size - ex_sizes[0];
    
    elem_assign_at(s, s0.size, &s1);
    s1.size = ex_sizes[1];

  } else
  {
    elem_assign_at(s, 0, &s0);
    s0.size = ex_sizes[1];
    
    elem_assign_at(s, s0.size, &s1);
    s1.size = s->size - ex_sizes[0];
  }

#ifdef CHECK_ORDER
  check_order = elements_validate_order(&s0, 1);
  if (check_order) Z_ERROR("intermediate lower order failed at %" slint_fmt "", check_order);
  check_order = elements_validate_order(&s1, 1);
  if (check_order) Z_ERROR("intermediate higher order failed at %" slint_fmt "", check_order);
#endif

  s->size = s0.size + s1.size;

  /* local merge */
  rti_tstart(rti_tid_mpi_merge2_local);

  if (s0.size > 0 && s1.size > 0 && m2 != NULL)
  {
    Z_TRACE_IF(MM2_TRACE_IF, "local merge2 %" slint_fmt " with %" slint_fmt "", s0.size, s1.size);

    m2(&s0, &s1, xs);
  }

  rti_tstop(rti_tid_mpi_merge2_local);

#ifdef CHECK_ORDER
  check_order = elements_validate_order(s, 1);
  if (check_order) Z_ERROR("output order failed at %" slint_fmt "", check_order);
#endif

  rti_tstop(rti_tid_mpi_merge2);

#if defined(MM2_PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (MM2_PRINT_TIMINGS)
  {
    printf("%d: mpi_merge2: %f\n", rank, rti_tlast(rti_tid_mpi_merge2));
    printf("%d: mpi_merge2: find: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_find));
    printf("%d: mpi_merge2: move-right: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_moveright));
    printf("%d: mpi_merge2: exchange: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_exchange));
    printf("%d: mpi_merge2: move-left: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_moveleft));
    printf("%d: mpi_merge2: local: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_local));
  }
#endif

  return 0;
}


#undef CHECK_ORDER



/* sl_macro MMK_TRACE_IF */
/* sl_macro MMK_TRACE_FAIL */
/* sl_macro MMK_EQUAL_SYNC */
/* sl_macro MMK_SORTED_SYNC */
/* sl_macro MMK_SYNC */
/* sl_macro MMK_PRINT_TIMINGS */


#include "sl_common.h"


#ifdef SLDEBUG
# define CHECK_ORDER
#endif

/*#define CHECK_ORDER_BREAK*/


#ifndef MMK_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MMK_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MMK_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_mergek_equal(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_mergek_equal */
{
  slint_t stage, other_rank, up, high_rank;

#ifdef CHECK_ORDER
  slint_t local_order, global_order;
#endif

#ifdef MMK_TRACE_FAIL
  elements_t sin;
  char trace_file[128];
#endif


  Z_TRACE_IF(MMK_TRACE_IF, "starting mpi_mergek_equal");

  /* sl_tid rti_tid_mpi_mergek_equal */

  rti_treset(rti_tid_mpi_mergek_equal_while);         /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_equal_while_merge2);  /* sl_tid */

  rti_tstart(rti_tid_mpi_mergek_equal);

  if (size < 0) MPI_Comm_size(comm, &size);
  if (rank < 0) MPI_Comm_rank(comm, &rank);

#ifdef MMK_EQUAL_SYNC
    MPI_Barrier(comm);
#endif

  stage = 0;

  rti_tstart(rti_tid_mpi_mergek_equal_while);

  if (size > 1)
  while ((other_rank = (sn)(size, rank, stage, snd, &up)) >= 0)
  {
    Z_TRACE_IF(MMK_TRACE_IF, "stage: %" slint_fmt ", %d %s %" slint_fmt, stage, rank, ((rank == z_min(rank, other_rank))?(up?"<-":"->"):(up?"->":"<-")), other_rank);

#ifdef CHECK_ORDER
    local_order = elements_validate_order(s, 1);
    if (local_order) Z_ERROR("input order failed at %" slint_fmt "", local_order);

    local_order = (local_order)?1:0;
    MPI_Allreduce(&local_order, &global_order, 1, int_mpi_datatype, MPI_SUM, comm);
    if (global_order)
    {
# ifdef CHECK_ORDER_BREAK
      Z_ERROR("break: input order failed on %" slint_fmt " process(es)", global_order);
      break;
# endif
    }
#endif

#ifdef MMK_TRACE_FAIL
    elements_alloc(&sin, s->size, SLCM_ALL);
    elem_ncopy(s, &sin, s->size);
#endif

    rti_tstart(rti_tid_mpi_mergek_equal_while_merge2);
    high_rank = (up)?(z_min(rank, other_rank)):(z_max(rank, other_rank));
    mpi_merge2(s, other_rank, high_rank, NULL, m2x, xs, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_equal_while_merge2);

#ifdef CHECK_ORDER
    local_order = elements_validate_order(s, 1);
    if (local_order) Z_ERROR("output order failed at %" slint_fmt "", local_order);

    local_order = (local_order)?1:0;
    MPI_Allreduce(&local_order, &global_order, 1, int_mpi_datatype, MPI_SUM, comm);
    if (global_order)
    {
#ifdef MMK_TRACE_FAIL
      sprintf(trace_file, "mergek_trace_rank%d_stage%" slint_fmt "_other%" slint_fmt "_high%" slint_fmt ".in", rank, stage, other_rank, high_rank);
      elements_save_keys_to_file(&sin, trace_file);
      sprintf(trace_file, "mergek_trace_rank%d_stage%" slint_fmt "_other%" slint_fmt "_high%" slint_fmt ".out", rank, stage, other_rank, high_rank);
      elements_save_keys_to_file(s, trace_file);

# ifdef CHECK_ORDER_BREAK
      elements_free(&sin);
# endif
#endif

# ifdef CHECK_ORDER_BREAK
      Z_ERROR("break: output order failed on %" slint_fmt " process(es)", global_order);
      break;
# endif
    }
#endif

#ifdef MMK_TRACE_FAIL
    elements_free(&sin);
#endif

#ifdef MMK_EQUAL_SYNC
    MPI_Barrier(comm);
#endif

    ++stage;
  }

  rti_tstop(rti_tid_mpi_mergek_equal_while);

  rti_tstop(rti_tid_mpi_mergek_equal);

#if defined(MMK_PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (MMK_PRINT_TIMINGS)
  {
    printf("%d: mpi_mergek_equal: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_equal));
    printf("%d: mpi_mergek_equal: while: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_equal_while));
    printf("%d: mpi_mergek_equal:  merge2: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_equal_while_merge2));
    printf("%d: mpi_mergek_equal: stages: %" slint_fmt "\n", rank, stage);
  }
#endif

  return stage;
}


slint_t mpi_mergek_sorted(elements_t *s, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_mergek_sorted */
{
  slint_t stage, done;


  Z_TRACE_IF(MMK_TRACE_IF, "starting mpi_mergek_sorted");

  /* sl_tid rti_tid_mpi_mergek_sorted */

  rti_treset(rti_tid_mpi_mergek_sorted_while);          /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_sorted_while_check);    /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_sorted_while_oddeven);  /* sl_tid */

  rti_tstart(rti_tid_mpi_mergek_sorted);

  if (size < 0) MPI_Comm_size(comm, &size);
  if (rank < 0) MPI_Comm_rank(comm, &rank);

#ifdef MMK_SORTED_SYNC
  MPI_Barrier(comm);
#endif

  rti_tstart(rti_tid_mpi_mergek_sorted_while);

  stage = 0;

  /* use alternating odd/even rounds (simliar to odd-even-tranpose) until sorted */
  if (size > 1)
  while (1)
  {
    rti_tstart(rti_tid_mpi_mergek_sorted_while_check);
    done = mpi_check_global_order(key_purify(s->keys[0]), key_purify(s->keys[s->size - 1]), -1, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_sorted_while_check);

    if (done) break;

    Z_TRACE_IF(MMK_TRACE_IF, "stage: %" slint_fmt " order failed, do %s", stage, ((stage % 2)?"odd":"even"));

    rti_tstart(rti_tid_mpi_mergek_sorted_while_oddeven);
    mpi_mergek_equal(s, (stage % 2)?sn_odd:sn_even, NULL, m2x, xs, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_sorted_while_oddeven);

#ifdef MMK_SORTED_SYNC
    MPI_Barrier(comm);
#endif

    ++stage;
  }

  rti_tstop(rti_tid_mpi_mergek_sorted_while);

  rti_tstop(rti_tid_mpi_mergek_sorted);

#if defined(MMK_PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (MMK_PRINT_TIMINGS)
  {
    printf("%d: mpi_mergek_sorted: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_sorted));
    printf("%d: mpi_mergek_sorted: while: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_sorted_while));
    printf("%d: mpi_mergek_sorted:  check: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_sorted_while_check));
    printf("%d: mpi_mergek_sorted:  oddeven: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_sorted_while_oddeven));
    printf("%d: mpi_mergek_sorted: stages: %" slint_fmt "\n", rank, stage);
  }
#endif

  return stage;
}


slint_t mpi_mergek_sorted2(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_mergek_sorted2 */
{
  slint_t stage, other_rank, up, high_rank, done;


  Z_TRACE_IF(MMK_TRACE_IF, "starting mpi_mergek_sorted2");

  /* sl_tid rti_tid_mpi_mergek_sorted2 */

  rti_treset(rti_tid_mpi_mergek_sorted2_while);          /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_sorted2_while_check);    /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_sorted2_while_oddeven);  /* sl_tid */

  rti_tstart(rti_tid_mpi_mergek_sorted2);

  if (size < 0) MPI_Comm_size(comm, &size);
  if (rank < 0) MPI_Comm_rank(comm, &rank);

#ifdef MMK_SORTED_SYNC
  MPI_Barrier(comm);
#endif

  rti_tstart(rti_tid_mpi_mergek_sorted2_while);

  stage = 0;

  if (size > 1)
  while ((other_rank = (sn)(size, rank, stage, snd, &up)) >= 0)
  {
    rti_tstart(rti_tid_mpi_mergek_sorted2_while_check);
    done = mpi_check_global_order(key_purify(s->keys[0]), key_purify(s->keys[s->size - 1]), -1, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_sorted2_while_check);

    if (done) break;

    Z_TRACE_IF(MMK_TRACE_IF, "stage: %" slint_fmt " order failed, do %s", stage, ((stage % 2)?"odd":"even"));

    rti_tstart(rti_tid_mpi_mergek_sorted2_while_oddeven);
    high_rank = (up)?(z_min(rank, other_rank)):(z_max(rank, other_rank));
    mpi_merge2(s, other_rank, high_rank, NULL, m2x, xs, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_sorted2_while_oddeven);

#ifdef MMK_SORTED_SYNC
    MPI_Barrier(comm);
#endif

    ++stage;
  }

  rti_tstop(rti_tid_mpi_mergek_sorted2_while);

  rti_tstop(rti_tid_mpi_mergek_sorted2);

#if defined(MMK_PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (MMK_PRINT_TIMINGS)
  {
    printf("%d: mpi_mergek_sorted2: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_sorted2));
    printf("%d: mpi_mergek_sorted2: while: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_sorted2_while));
    printf("%d: mpi_mergek_sorted2:  check: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_sorted2_while_check));
    printf("%d: mpi_mergek_sorted2:  oddeven: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_sorted2_while_oddeven));
    printf("%d: mpi_mergek_sorted2: stages: %" slint_fmt "\n", rank, stage);
  }
#endif

  return stage;
}


slint_t mpi_mergek(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_mergek */
{
  slint_t stage, done;

  /* FIXME: use something more efficient that works for all kind of (non-equal) distributions? */


  Z_TRACE_IF(MMK_TRACE_IF, "starting mpi_mergek");

  /* sl_tid rti_tid_mpi_mergek */

  rti_treset(rti_tid_mpi_mergek_equalike);       /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_while);          /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_while_check);    /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_while_oddeven);  /* sl_tid */

  rti_tstart(rti_tid_mpi_mergek);

  /* assume equal distribution */
  rti_tstart(rti_tid_mpi_mergek_equalike);
  stage = mpi_mergek_equal(s, sn, snd, m2x, xs, size, rank, comm);
  rti_tstop(rti_tid_mpi_mergek_equalike);

/*  return finalize;*/

#ifdef MMK_SYNC
  MPI_Barrier(comm);
#endif

  rti_tstart(rti_tid_mpi_mergek_while);

  /* use alternating odd/even rounds (simliar to odd-even-tranpose) until sorted */
  if (size > 1)
  while (1)
  {
    rti_tstart(rti_tid_mpi_mergek_while_check);
    done = mpi_check_global_order(key_purify(s->keys[0]), key_purify(s->keys[s->size - 1]), -1, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_while_check);

    if (done) break;

    Z_TRACE_IF(MMK_TRACE_IF, "stage: %" slint_fmt " order failed, do %s", stage, ((stage % 2)?"odd":"even"));

    rti_tstart(rti_tid_mpi_mergek_while_oddeven);
    mpi_mergek_equal(s, (stage % 2)?sn_odd:sn_even, NULL, m2x, xs, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_while_oddeven);

#ifdef MMK_SYNC
    MPI_Barrier(comm);
#endif

    ++stage;
  }

  rti_tstop(rti_tid_mpi_mergek_while);

  rti_tstop(rti_tid_mpi_mergek);

#if defined(MMK_PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (MMK_PRINT_TIMINGS)
  {
    printf("%d: mpi_mergek: %f\n", rank, rti_tlast(rti_tid_mpi_mergek));
    printf("%d: mpi_mergek: equal: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_equalike));
    printf("%d: mpi_mergek: while: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_while));
    printf("%d: mpi_mergek:  check: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_while_check));
    printf("%d: mpi_mergek:  oddeven: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_while_oddeven));
    printf("%d: mpi_mergek: stages: %" slint_fmt "\n", rank, stage);
  }
#endif

  return stage;
}


slint_t mpi_mergek_equal2(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms) /* sl_proto, sl_func mpi_mergek_equal2 */
{
  slint_t stage = 0, other_ranks[2], ups[2], high_ranks[2];

  rti_tstart(rti_tid_mpi_mergek);

  if ((comms[0] != MPI_COMM_NULL && sizes[0] > 1) || (comms[1] != MPI_COMM_NULL && sizes[1] > 1))
  while (1)
  {
    other_ranks[0] = sn(sizes[0], ranks[0], stage, snd, &ups[0]);
    other_ranks[1] = sn(sizes[1], ranks[1], stage, snd, &ups[1]);
    
/*    printf("%d: doing merge2 %d vs. %d and %d vs. %d\n", sl_mpi_rank, ranks[0], other_ranks[0], ranks[1], other_ranks[1]);*/

    if ((other_ranks[0] < 0) && (other_ranks[1] < 0)) break;

    if (other_ranks[0] >= 0)
    {
      high_ranks[0] = (ups[0])?(z_min(ranks[0], other_ranks[0])):(z_max(ranks[0], other_ranks[0]));
      Z_TRACE_IF(MMK_TRACE_IF, "#0 merge2 with %" slint_fmt ", high: %" slint_fmt ", size = %d, rank = %d", other_ranks[0], high_ranks[0], sizes[0], ranks[0]);
      mpi_merge2(&s[0], other_ranks[0], high_ranks[0], NULL, m2x, xs, sizes[0], ranks[0], comms[0]);
    }

    if (other_ranks[1] >= 0)
    {
      high_ranks[1] = (ups[1])?(z_min(ranks[1], other_ranks[1])):(z_max(ranks[1], other_ranks[1]));
      Z_TRACE_IF(MMK_TRACE_IF, "#1 merge2 with %" slint_fmt ", high: %" slint_fmt ", size = %d, rank = %d", other_ranks[1], high_ranks[1], sizes[1], ranks[1]);
      mpi_merge2(&s[1], other_ranks[1], high_ranks[1], NULL, m2x, xs, sizes[1], ranks[1], comms[1]);
    }

    rti_tstop(rti_tid_mpi_mergek_merge2);

    ++stage;
  }

  rti_tstop(rti_tid_mpi_mergek);

  return stage;
}


#undef CHECK_ORDER
#undef CHECK_ORDER_BREAK



/* sl_macro MPEG_TRACE_IF */


#include "sl_common.h"


/* config */
/*#define PRINT_SCOUNTS_RCOUNTS
#define PRINT_STATS*/
/*#define PRINT_TIMINGS  0*/


#ifndef MPEG_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MPEG_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MPEG_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_partition_exact_generic(elements_t *s, partcond_t *pcond, binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_exact_generic */
{
  splitter_t sp;
#ifdef PRINT_SCOUNTS_RCOUNTS
  slint_t i, j;
#endif
#ifdef PRINT_STATS
  int sdispls[size];
#endif


  rti_treset(rti_tid_mpi_partition_exact_generic);          /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_generic_select);   /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_generic_rcounts);  /* sl_tid */


  rti_tstart(rti_tid_mpi_partition_exact_generic);

  rti_tstart(rti_tid_mpi_partition_exact_generic_select);

  sp.displs = scounts;
  mpi_select_exact_generic_grouped(s, 1, pcond, comm, MPI_COMM_NULL, bm, &sp, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_exact_generic_select);

  rti_tstart(rti_tid_mpi_partition_exact_generic_rcounts);

  /* create scounts from sdispls */
  displs2counts(size, scounts, NULL, s->size);

  /* create rcounts if necessary */
  if (rcounts) mpi_xcounts2ycounts_all2all(scounts, rcounts, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_exact_generic_rcounts);

  rti_tstop(rti_tid_mpi_partition_exact_generic);

#ifdef PRINT_SCOUNTS_RCOUNTS
  printf("%d: scounts:", rank);
  for (i = 0, j = 0; i < size; ++i) { printf(" %d ", scounts[i]); j += scounts[i]; }
  printf(" = %" sl_int_type_fmt "\n", j);
  printf("%d: rcounts:", rank);
  if (rcounts) for (i = 0, j = 0; i < size; ++i) { printf(" %d ", rcounts[i]); j += rcounts[i]; }
  printf("\n");
#endif

#ifdef PRINT_STATS
  counts2displs(size, scounts, sdispls);
  mpi_select_stats(s, size, sdispls, size, rank, comm);
#endif

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == PRINT_TIMINGS)
  {
    printf("%d: mpi_partition_exact_generic: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_generic));
    printf("%d: mpi_partition_exact_generic: select: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_generic_select));
    printf("%d: mpi_partition_exact_generic: rcounts: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_generic_rcounts));
  }
#endif

  return 0;
}


#undef PRINT_SCOUNTS_RCOUNTS
#undef PRINT_STATS
#undef PRINT_TIMINGS



/* sl_macro MPER_TRACE_IF */


#include "sl_common.h"


#ifndef MPER_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MPER_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MPER_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_partition_exact_radix(elements_t *s, partcond_t *pcond, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_exact_radix */
{
  binning_t bm;


  binning_radix_create(&bm, rhigh, rlow, rwidth, sorted|SL_SORTED_OUT);

  mpi_partition_exact_generic(s, pcond, &bm, scounts, rcounts, size, rank, comm);

  binning_radix_destroy(&bm);

  return 0;
}



/* sl_macro MPERG_TRACE_IF */


#include "sl_common.h"


/*#define PRINT_TIMINGS*/

#ifndef MPERG_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MPERG_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MPERG_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


#ifdef SL_INDEX

slint_t mpi_partition_exact_radix_ngroups(elements_t *s, partcond_t *pcond, slint_t ngroups, MPI_Comm *group_comms, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_exact_radix_ngroups */
{
  slint_t g, i;

  slint_t lgcounts[2], coffset, coffset_add;
#ifdef elem_weight
  slweight_t lgweights[2];
#endif

  slint_t sx_max_size;
  elements_t _sx;
  elements_t *es, ed0, ed1;
  
  int group_sizes[ngroups];
  int group_ranks[ngroups];
  
  MPI_Comm master_comms[ngroups];
  int master_sizes[ngroups];
  int master_ranks[ngroups];
  
  partcond_p group_pconds[ngroups];
  partcond_t merged_pconds[ngroups];

  int group_sdispls[size], group_scounts[size];
  int current_scounts[size], current_sdispls[size], current_rcounts[size], current_rdispls[size];

  int *_scounts, *_rcounts;

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif


  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups);                  /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_pconds);           /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_idxin);            /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_up);               /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down);             /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down_select);      /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down_alltoall);    /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down_x2suby);      /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down_merge);       /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_idxout);           /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_idxout_loop);      /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_idxout_alltoall);  /* sl_tid */


  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups);

  if (scounts == NULL) _scounts = z_alloca(size, sizeof(int)); else _scounts = scounts;
  if (rcounts == NULL) _rcounts = z_alloca(size, sizeof(int)); else _rcounts = rcounts;

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_pconds);

#ifdef elem_weight
  doweights = ((pcond->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

  /* make absolute pconds */
#ifdef elem_weight
  if (doweights) mpi_elements_get_counts_and_weights(s, 1, lgcounts, lgweights, -1, size, rank, comm);
  else
#endif
    mpi_elements_get_counts(s, &lgcounts[0], &lgcounts[1], -1, size, rank, comm);

  init_partconds(1, pcond, size, lgcounts[1], elem_weight_ifelse(doweights?lgweights[1]:0, 0));

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_pconds);

  sx_max_size = lgcounts[1];
  
  if (pcond->pcm & SLPC_COUNTS_LH) sx_max_size = z_min(sx_max_size, pcond->count_high - pcond->count_low);
  if (pcond->pcm & SLPC_COUNTS_MM) sx_max_size = z_minmax(pcond->count_min, sx_max_size, pcond->count_max);

  if (sx == NULL || sx->size < 2 * sx_max_size)
  {
    sx = &_sx;

    /* allocate scratch elements (need keys, indices and weights) */
    elements_alloc(sx, 2 * sx_max_size, SLCM_KEYS|SLCM_INDICES|SLCM_WEIGHTS);
  }

  /* use scratch for intermediate elements */  
  elem_assign(sx, &ed0);
  elem_assign_at(&ed0, sx_max_size, &ed1);

  es = s;

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_idxin);

  /* assign source ranks */
  for (i = 0; i < es->size; ++i) es->indices[i] = rank;

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_idxin);


  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_up);

  /* collect and create merged pconds (bottom-up) */
  for (g = ngroups - 1; g >= 0; --g)
  {
    MPI_Comm_size(group_comms[g], &group_sizes[g]);
    MPI_Comm_rank(group_comms[g], &group_ranks[g]);
    
    if (g < ngroups - 1)
    {
      MPI_Comm_split(group_comms[g], (group_ranks[g + 1] == 0)?0:MPI_UNDEFINED, group_ranks[g], &master_comms[g]);
      if (master_comms[g] != MPI_COMM_NULL)
      {
        MPI_Comm_size(master_comms[g], &master_sizes[g]);
        MPI_Comm_rank(master_comms[g], &master_ranks[g]);

      } else master_ranks[g] = -1;

      MPI_Bcast(&master_sizes[g], 1, MPI_INT, 0, group_comms[g + 1]);

    } else
    {
      master_comms[g] = group_comms[g];
      master_sizes[g] = group_sizes[g];
      master_ranks[g] = group_ranks[g];
    }
    
    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": group: %d of %d, master: %d of %d", g, group_ranks[g], group_sizes[g], master_ranks[g], master_sizes[g]);

    group_pconds[g] = z_alloca(master_sizes[g], sizeof(partcond_t));

    if (master_comms[g] != MPI_COMM_NULL) mpi_allgather_partconds((g < ngroups - 1)?&merged_pconds[g + 1]:pcond, group_pconds[g], master_sizes[g], master_ranks[g], master_comms[g]);

    if (g < ngroups - 1) mpi_bcast_partconds(master_sizes[g], group_pconds[g], 0, group_sizes[g + 1], group_ranks[g + 1], group_comms[g + 1]);

    merge_partconds(group_pconds[g], master_sizes[g], &merged_pconds[g]);
    
    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": merged_pconds: %" slint_fmt "", g, merged_pconds[g].pcm);
    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": merged_pconds: count: min/max: %f/%f - low/high: %f/%f", g, merged_pconds[g].count_min, merged_pconds[g].count_max, merged_pconds[g].count_low, merged_pconds[g].count_high);
#ifdef elem_weight
    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": merged_pconds: weight: min/max: %f/%f - low/high: %f/%f", g, merged_pconds[g].weight_min, merged_pconds[g].weight_max, merged_pconds[g].weight_low, merged_pconds[g].weight_high);
#endif
  }

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_up);

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down);

  coffset = 0;

  /* do selects (top-down) */
  for (g = 0; g < ngroups; ++g)
  {
    if (g > 0)
    {
      if (group_pconds[g][0].pcm & SLPC_COUNTS_LH)
      {
        coffset_add = 0;
        MPI_Exscan(&es->size, &coffset_add, 1, int_mpi_datatype, MPI_SUM, group_comms[g - 1]);
        MPI_Bcast(&coffset_add, 1, int_mpi_datatype, 0, group_comms[g]);
    
        Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": count offset = %" slint_fmt " + %" slint_fmt " = %" slint_fmt, g, coffset, coffset_add, coffset + coffset_add);
  
        coffset += coffset_add;

        for (i = 0; i < master_sizes[g]; ++i)
        {
          Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": correct group pcond %" slint_fmt ": count: min/max: %f/%f -> %f/%f",
            g, i, group_pconds[g][i].count_low, group_pconds[g][i].count_high, group_pconds[g][i].count_low - coffset, group_pconds[g][i].count_high - coffset);

          group_pconds[g][i].count_low -= coffset;
          group_pconds[g][i].count_high -= coffset;
        }
      }

#ifdef elem_weight
      if (group_pconds[g][0].pcm & SLPC_WEIGHTS_LH)
      {
        Z_ERROR("correction of weighted low/high group pconds not implemented!");
      }
#endif
    }
    
    rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down_select);
    mpi_select_exact_radix(es, 1, master_sizes[g], group_pconds[g], rhigh, rlow, rwidth, SL_SORTED_IN, group_sdispls, group_sizes[g], group_ranks[g], group_comms[g]);
    rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down_select);

    Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, master_sizes[g], "%d  ", group_sdispls[i], "%" slint_fmt ": group_sdispls = ", g);

    rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down_x2suby);

    if (g < ngroups - 1)
    {
      /* create scounts from sdispls */
      displs2counts(master_sizes[g], group_sdispls, group_scounts, es->size);

      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, master_sizes[g], "%d  ", group_scounts[i], "%" slint_fmt ": group_scounts = ", g);

      mpi_xcounts2ycounts_grouped(group_scounts, master_sizes[g], current_rcounts, group_comms[g + 1], master_comms[g], group_sizes[g], group_ranks[g], group_comms[g]);
    
      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, group_sizes[g], "%d  ", current_rcounts[i], "%" slint_fmt ": current_rcounts = ", g);

/*#define SPARSE_SCOUNTS_FROM_RCOUNTS*/

      /* make scounts from rcounts */
#ifdef SPARSE_SCOUNTS_FROM_RCOUNTS
      mpi_xcounts2ycounts_sparse(current_rcounts, current_scounts, es->size, group_sizes[g], group_ranks[g], group_comms[g]);
#else
      mpi_xcounts2ycounts_all2all(current_rcounts, current_scounts, group_sizes[g], group_ranks[g], group_comms[g]);
#endif

#undef SPARSE_SCOUNTS_FROM_RCOUNTS

      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, group_sizes[g], "%d  ", current_scounts[i], "%" slint_fmt ": current_scounts = ", g);

    } else
    {
      /* create scounts from sdispls */
      displs2counts(master_sizes[g], group_sdispls, current_scounts, es->size);

      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, group_sizes[g], "%d  ", current_scounts[i], "%" slint_fmt ": current_scounts = ", g);

      /* make rcounts from scounts */
      mpi_xcounts2ycounts_all2all(current_scounts, current_rcounts, group_sizes[g], group_ranks[g], group_comms[g]);

      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, group_sizes[g], "%d  ", current_rcounts[i], "%" slint_fmt ": current_rcounts = ", g);
    }

    rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down_x2suby);

    counts2displs(group_sizes[g], current_scounts, current_sdispls);
    counts2displs(group_sizes[g], current_rcounts, current_rdispls);

    Z_ASSERT(current_sdispls[group_sizes[g] - 1] + current_scounts[group_sizes[g] - 1] <= pcond->count_max);
    Z_ASSERT(current_rdispls[group_sizes[g] - 1] + current_rcounts[group_sizes[g] - 1] <= pcond->count_max);

    rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down_alltoall);
#define xelem_data_not
#define xelem_data_weight
#define xelem_call \
    MPI_Alltoallv(xelem_buf(es), current_scounts, current_sdispls, xelem_mpi_datatype, xelem_buf(&ed0), current_rcounts, current_rdispls, xelem_mpi_datatype, group_comms[g]);
#include "sl_xelem_call.h"
    rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down_alltoall);

    ed0.size = current_rdispls[group_sizes[g] - 1] + current_rcounts[group_sizes[g] - 1];

    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": received: %" slint_fmt, g, ed0.size);

    rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down_merge);

    /* merge received elements (but not in the last round) */
    if (g < ngroups - 1)
    {
      mergep_heap_int(&ed0, &ed1, group_sizes[g], current_rdispls, NULL);

      ed1.size = ed0.size;

      es = &ed1;
    
    } else es = &ed0;

    rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down_merge);

    z_freea(group_pconds[g]);
  }

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down);

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_idxout);

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_idxout_loop);

  /* create rcounts from source ranks (indexes) */
  for (i = 0; i < size; ++i) _rcounts[i] = 0;
  for (i = 0; i < es->size; ++i) ++_rcounts[es->indices[i]];

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, size, "%d  ", _rcounts[i], "rcounts = ");

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_idxout_loop);

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_idxout_alltoall);

  if (scounts)
  {
    /* make scounts from rcounts */
    mpi_xcounts2ycounts_all2all(_rcounts, _scounts, size, rank, comm);

    Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, size, "%d  ", _scounts[i], "scounts = ");
  }

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_idxout_alltoall);

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_idxout);

  if (sx == &_sx) elements_free(sx);

  if (scounts == NULL) z_freea(_scounts);
  if (rcounts == NULL) z_freea(_rcounts);

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups);

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == 0)
  {
    printf("%d: mpi_partition_exact_radix_ngroups: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups));
    printf("%d: mpi_partition_exact_radix_ngroups: pconds: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_pconds));
    printf("%d: mpi_partition_exact_radix_ngroups: idxin: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_idxin));
    printf("%d: mpi_partition_exact_radix_ngroups: up: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_up));
    printf("%d: mpi_partition_exact_radix_ngroups: down: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down));
    printf("%d: mpi_partition_exact_radix_ngroups:  select: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down_select));
    printf("%d: mpi_partition_exact_radix_ngroups:  alltoall: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down_alltoall));
    printf("%d: mpi_partition_exact_radix_ngroups:  x2suby: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down_x2suby));
    printf("%d: mpi_partition_exact_radix_ngroups:  merge: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down_merge));
    printf("%d: mpi_partition_exact_radix_ngroups: idxout: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_idxout));
    printf("%d: mpi_partition_exact_radix_ngroups:  loop: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_idxout_loop));
    printf("%d: mpi_partition_exact_radix_ngroups:  alltoall: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_idxout_alltoall));
  }
#endif

  return 0;

#undef doweights
}

#endif /* SL_INDEX */


slint_t mpi_partition_exact_radix_2groups(elements_t *s, partcond_t *pcond, MPI_Comm group_comm, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_exact_radix_2groups */
{
  slint_t i, j;

  slint_t lgcounts[2], coffset;
#ifdef elem_weight
  slweight_t lgweights[2];
#endif

  int group_size, group_rank;
  MPI_Comm master_comm;
  int master_size, master_rank;

  slint_t nparts;
  int group_sdispls[size], group_scounts[size];

  partcond_t *group_pconds, group_pcond;

  int sdispls[size], rdispls[size];
  slint_t rcounts_total, nsubelements;

  elements_t _sx;

  slint_t sub_elements_sources[size], sub_elements_sizes[size], sub_elements_size;
  elements_t sub_elements[size];

  int *sub_sdispls;

  int *_scounts, *_rcounts;

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif


  rti_treset(rti_tid_mpi_partition_exact_radix_2groups);            /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_pconds);     /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_select1st);  /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_x2suby);     /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_alltoall);   /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_select2nd);  /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_subx2y);     /* sl_tid */


  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups);

  if (scounts == NULL) _scounts = z_alloca(size, sizeof(int)); else _scounts = scounts;
  if (rcounts == NULL) _rcounts = z_alloca(size, sizeof(int)); else _rcounts = rcounts;

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_pconds);

#ifdef elem_weight
  doweights = ((pcond->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

  /* make absolute pconds */
#ifdef elem_weight
  if (doweights) mpi_elements_get_counts_and_weights(s, 1, lgcounts, lgweights, -1, size, rank, comm);
  else
#endif
    mpi_elements_get_counts(s, &lgcounts[0], &lgcounts[1], -1, size, rank, comm);

  init_partconds(1, pcond, size, lgcounts[1], elem_weight_ifelse(doweights?lgweights[1]:0, 0));

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_pconds);

  /* init group */
  if (group_comm != MPI_COMM_NULL)
  {
    MPI_Comm_size(group_comm, &group_size);
    MPI_Comm_rank(group_comm, &group_rank);

  } else
  {
    group_size = 1;
    group_rank = 0;
  }
  
  /* init master */
  MPI_Comm_split(comm, (group_rank == 0)?0:MPI_UNDEFINED, rank, &master_comm);
  if (master_comm != MPI_COMM_NULL)
  {
    MPI_Comm_size(master_comm, &master_size);
    MPI_Comm_rank(master_comm, &master_rank);

  } else master_rank = -1;

  /* distribute num. of masters */
  if (group_comm != MPI_COMM_NULL) MPI_Bcast(&master_size, 1, MPI_INT, 0, group_comm);

  Z_TRACE_IF(MPERG_TRACE_IF, "%d: group: %d of %d, master: %d of %d", rank, group_rank, group_size, master_rank, master_size);

  nparts = master_size;

  /* create group partcond */
  group_pconds = z_alloca(group_size, sizeof(partcond_t));

  mpi_gather_partconds_grouped(pcond, group_comm, MPI_COMM_NULL, group_pconds, NULL, size, rank, comm);

  merge_partconds(group_pconds, group_size, &group_pcond);

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_select1st);

  /* perform 1st grouped select */  
  mpi_select_exact_radix_grouped(s, 1, &group_pcond, master_comm, group_comm, rhigh, rlow, rwidth, SL_SORTED_IN, group_sdispls, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_select1st);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, nparts, "%d  ", group_sdispls[i], "group_sdispls = ");

  /* create scounts from sdispls */
  displs2counts(nparts, group_sdispls, group_scounts, s->size);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, nparts, "%d  ", group_scounts[i], "group_scounts = ");

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_x2suby);

  mpi_xcounts2ycounts_grouped(group_scounts, nparts, _rcounts, group_comm, master_comm, size, rank, comm);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, size, "%d  ", _rcounts[i], "rcounts = ");

#define SPARSE_SCOUNTS_FROM_RCOUNTS

  /* make scounts from rcounts */
#ifdef SPARSE_SCOUNTS_FROM_RCOUNTS
  mpi_xcounts2ycounts_sparse(_rcounts, _scounts, s->size, size, rank, comm);
#else
  mpi_xcounts2ycounts_all2all(_rcounts, _scounts, size, rank, comm);
#endif

#undef SPARSE_SCOUNTS_FROM_RCOUNTS

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, size, "%d  ", _scounts[i], "scounts = ");

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_x2suby);

  /* create displs from counts */
  sdispls[0] = rdispls[0] = 0;
  nsubelements = (_rcounts[0] != 0);;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + _scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + _rcounts[i - 1];

    /* determine number of sub lists to receive */
    nsubelements += (_rcounts[i] != 0);
  }
  rcounts_total = rdispls[size - 1] + _rcounts[size - 1];

  Z_TRACE_IF(MPERG_TRACE_IF, "rcounts_total = %" slint_fmt, rcounts_total);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, size, "%d  ", sdispls[i], "sdispls = ");
  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, size, "%d  ", rdispls[i], "rdispls = ");

  if (sx == NULL || sx->size < rcounts_total)
  {
    sx = &_sx;

    /* allocate elements for 1st redistribution (need keys and weights) */
    elements_alloc(sx, rcounts_total, SLCM_KEYS|SLCM_WEIGHTS);
  }

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_alltoall);

  /* 1st redistribution */
#define xelem_index_not
#define xelem_data_not
#define xelem_data_weight
#define xelem_call \
  MPI_Alltoallv(xelem_buf(s), _scounts, sdispls, xelem_mpi_datatype, xelem_buf(sx), _rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_alltoall);

  /* create sub lists */
  j = 0;
  sub_elements_size = 0;
  for (i = 0; i < size; ++i)
  if (_rcounts[i] != 0)
  {
    sub_elements_sources[j] = i;
    sub_elements_sizes[j] = _rcounts[i];

    sub_elements_size += sub_elements_sizes[j];

    elem_assign_at(sx, rdispls[i], &sub_elements[j]);
    sub_elements[j].size = _rcounts[i];

    ++j;
  }

  sub_sdispls = z_alloca(group_size * nsubelements, sizeof(int));

  if (group_pconds[0].pcm & SLPC_COUNTS_LH)
  {
    coffset = 0;
    MPI_Exscan(&sub_elements_size, &coffset, 1, int_mpi_datatype, MPI_SUM, comm);
    MPI_Bcast(&coffset, 1, int_mpi_datatype, 0, group_comm);

    Z_TRACE_IF(MPERG_TRACE_IF, "count offset = %" slint_fmt, coffset);

    for (i = 0; i < group_size; ++i)
    {
      Z_TRACE_IF(MPERG_TRACE_IF, "correct group pcond %" slint_fmt ": count: min/max: %f/%f -> %f/%f",
        i, group_pconds[i].count_low, group_pconds[i].count_high, group_pconds[i].count_low - coffset, group_pconds[i].count_high - coffset);

      group_pconds[i].count_low -= coffset;
      group_pconds[i].count_high -= coffset;
    }
  }

#ifdef elem_weight
  if (group_pconds[0].pcm & SLPC_WEIGHTS_LH)
  {
    Z_ERROR("correction of weighted low/high group pconds not implemented!");
  }
#endif

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_select2nd);

  /* perform 2nd select */
  mpi_select_exact_radix(sub_elements, nsubelements, group_size, group_pconds, rhigh, rlow, rwidth, 0, sub_sdispls, group_size, group_rank, group_comm);

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_select2nd);

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_subx2y);

  mpi_subxdispls2ycounts(nsubelements, sub_sdispls, sub_elements_sources, sub_elements_sizes, group_comm, group_size, _rcounts, size, rank, comm);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, size, "%d  ", _rcounts[i], "rcounts = ");

  if (scounts)
  {
    /* make scounts from rcounts */
    mpi_xcounts2ycounts_all2all(_rcounts, _scounts, size, rank, comm);

    Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, i, size, "%d  ", _scounts[i], "scounts = ");
  }

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_subx2y);

  z_freea(sub_sdispls);

  if (sx == &_sx) elements_free(sx);

  z_freea(group_pconds);

  if (master_comm != MPI_COMM_NULL) MPI_Comm_free(&master_comm);

  if (scounts == NULL) z_freea(_scounts);
  if (rcounts == NULL) z_freea(_rcounts);

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups);

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == 0)
  {
    printf("%d: mpi_partition_exact_radix_2groups: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups));
    printf("%d: mpi_partition_exact_radix_2groups: pconds: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_pconds));
    printf("%d: mpi_partition_exact_radix_2groups: select1st: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_select1st));
    printf("%d: mpi_partition_exact_radix_2groups: x2suby: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_x2suby));
    printf("%d: mpi_partition_exact_radix_2groups: alltoall: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_alltoall));
    printf("%d: mpi_partition_exact_radix_2groups: select2nd: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_select2nd));
    printf("%d: mpi_partition_exact_radix_2groups: subx2y: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_subx2y));
  }
#endif

  return 0;

#undef doweights
}


#undef PRINT_TIMINGS



/* sl_macro MPS_TRACE_IF */


#include "sl_common.h"


/* config */
/*#define PRINT_SCOUNTS_RCOUNTS*/
/*#define PRINT_STATS*/
/*#define PRINT_TIMINGS  0*/


#ifndef MPS_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MPS_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MPS_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_partition_sample_regular(elements_t *s, partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_sample_regular */
{
  partcond_t pconds[size];

#ifdef PRINT_SCOUNTS_RCOUNTS
  slint_t i, j;
#endif
#ifdef PRINT_STATS
  int sdispls[size];
#endif


  rti_treset(rti_tid_mpi_partition_sample);          /* sl_tid */
  rti_treset(rti_tid_mpi_partition_sample_select);   /* sl_tid */
  rti_treset(rti_tid_mpi_partition_sample_rcounts);  /* sl_tid */


  rti_tstart(rti_tid_mpi_partition_sample);

  rti_tstart(rti_tid_mpi_partition_sample_select);

  splitter_t sp;
  sp.displs = scounts;

  mpi_allgather_partconds(pcond, pconds, size, rank, comm);

  mpi_select_sample_regular(s, size, pconds, size - 1, &sp, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_sample_select);
  
  rti_tstart(rti_tid_mpi_partition_sample_rcounts);

  /* create scounts from sdispls */
  displs2counts(size, scounts, NULL, s->size);

  /* create rcounts if necessary */
  if (rcounts) mpi_xcounts2ycounts_all2all(scounts, rcounts, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_sample_rcounts);

  rti_tstop(rti_tid_mpi_partition_sample);

#ifdef PRINT_SCOUNTS_RCOUNTS
  printf("%d: scounts:", rank);
  for (i = 0, j = 0; i < size; ++i) { printf(" %d ", scounts[i]); j += scounts[i]; }
  printf(" = %" sl_int_type_fmt "\n", j);
  printf("%d: rcounts:", rank);
  if (rcounts)
  {
    for (i = 0, j = 0; i < size; ++i) { printf(" %d ", rcounts[i]); j += rcounts[i]; }
    printf(" = %" sl_int_type_fmt "\n", j);
  }
#endif

#ifdef PRINT_STATS
  counts2displs(size, scounts, sdispls);
  mpi_select_stats(s, size, sdispls, size, rank, comm);
#endif

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == PRINT_TIMINGS)
  {
    printf("%d: mpi_partition_sample: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_sample));
    printf("%d: mpi_partition_sample: select: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_sample_select));
    printf("%d: mpi_partition_sample: rcounts: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_sample_rcounts));
  }
#endif

  return 0;
}


#undef PRINT_SCOUNTS_RCOUNTS
#undef PRINT_STATS
#undef PRINT_TIMINGS



#include "sl_common.h"


/*#define LOCAL_NCOPY*/
#define NONBLOCKING

#define MPI_REBALANCE_ALLTOALLV_TAG  1


slint_t mpi_rebalance(elements_t *s0, elements_t *s1, slint_t stable, slint_t *dst_size, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_rebalance */
{
  int exit_code = 0;
  int i;
  slint_t s, t, u, v;
  slint_t *sizes, local_sizes[2];
  int *sendcounts, *senddispls, *recvcounts, *recvdispls;
  slint_t send_rank, send_size, recv_rank, recv_size;

  
  sizes = z_alloc(size * 2, sizeof(slint_t));

  local_sizes[0] = s0->size;
  if (dst_size) local_sizes[1] = *dst_size; else local_sizes[1] = s0->size;

  MPI_Allgather(local_sizes, 2, int_mpi_datatype, sizes, 2, int_mpi_datatype, comm);

  s = t = 0;
  for (i = 0; i < size; ++i)
  {
    s += sizes[i * 2];
    if (sizes[i * 2 + 1] < 0) ++t; else s -= sizes[i * 2 + 1];
  }

  /* correct oversized 'sizes' */
  if (s < 0)
  {
    for (i = size - 1; i >= 0; --i)
    if (sizes[i * 2 + 1] > 0)
    {
      v = z_min(-s, sizes[i * 2 + 1]);
      sizes[i * 2 + 1] -= v;
      s += v;
    }
  }

  /* FIXME: APP_ASSERT(s >= 0); */

  /* correct autosize 'sizes' */
  if (t > 0)
  {
    u = s % t;
    for (i = 0; i < size; ++i)
    if (sizes[i * 2 + 1] < 0)
    {
      sizes[i * 2 + 1] = s / t;
      if (u > 0)
      {
        ++sizes[i * 2 + 1];
        --u;
      }
    }

  } else if (s > 0)
  {
    fprintf(stderr, "%d: mpi_rebalance: error: destination sizes too small (%" sl_int_type_fmt " > 0)\n", rank, s);
    exit_code = -1;
    goto free_and_exit;
  }

/*  printf("%d here: sizes = [", rank);
  for (i = 0; i < size; ++i) printf(" %" sl_int_type_fmt ":%" sl_int_type_fmt " ", sizes[i * 2], sizes[i * 2 + 1]);
  printf("]\n");*/

  sendcounts = z_alloc(4 * size, sizeof(int));
  senddispls = sendcounts + 1 * size;
  recvcounts = sendcounts + 2 * size;
  recvdispls = sendcounts + 3 * size;

  memset(sendcounts, 0, size * sizeof(int));
  memset(recvcounts, 0, size * sizeof(int));

  send_rank = 0;
  send_size = -1;
  recv_rank = 0;
  recv_size = -1;

  while (send_rank < size && recv_rank < size)
  {
    if (send_size < 0)
    {
      send_size = sizes[send_rank * 2];
      if (!stable) send_size -= sizes[send_rank * 2 + 1];
      if (send_size < 0) send_size = 0;
    }

    if (recv_size < 0)
    {
      recv_size = sizes[recv_rank * 2 + 1];
      if (!stable) recv_size -= sizes[recv_rank * 2];
      if (recv_size < 0) recv_size = 0;
    }

    s = z_min(send_size, recv_size);

    if (send_rank == rank) sendcounts[recv_rank] += s;
    if (recv_rank == rank) recvcounts[send_rank] += s;

    send_size -= s;
    recv_size -= s;

    if (send_size <= 0)
    {
      ++send_rank;
      send_size = -1;
    }

    if (recv_size <= 0)
    {
      ++recv_rank;
      recv_size = -1;
    }
  }
    
  /* FIXME: APP_ASSERT(send_rank >= size && recv_rank >= size); */

  if (stable) senddispls[0] = recvdispls[0] = 0;
  else
  {
    senddispls[0] = sizes[rank * 2 + 1];
    recvdispls[0] = sizes[rank * 2];
  }
  for (i = 1; i < size; ++i)
  {
    senddispls[i] = senddispls[i - 1] + sendcounts[i - 1];
    recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
  }
  if (!stable)
  {
    sendcounts[rank] = recvcounts[rank] = z_min(sizes[rank * 2], sizes[rank * 2 + 1]);
    senddispls[rank] = recvdispls[rank] = 0;
  }

/*  printf("%d here: sendcounts = [", rank);
  for (i = 0; i < size; ++i) printf(" %d ", sendcounts[i]);
  printf("]\n");

  printf("%d here: senddispls = [", rank);
  for (i = 0; i < size; ++i) printf(" %d ", senddispls[i]);
  printf("]\n");

  printf("%d here: recvcounts = [", rank);
  for (i = 0; i < size; ++i) printf(" %d ", recvcounts[i]);
  printf("]\n");

  printf("%d here: recvdispls = [", rank);
  for (i = 0; i < size; ++i) printf(" %d ", recvdispls[i]);
  printf("]\n");*/

  
  if (s1) /* out-of-place */
  {
#ifdef LOCAL_NCOPY
    if (!stable)
    {
      elem_ncopy_at(s0, senddispls[rank], s1, recvdispls[rank], sendcounts[rank]);
      sendcounts[rank] = recvcounts[rank] = 0;
    }
#endif

  } else /* in-place */
  {
    /* unstable case: local elements don't need to move  */
    if (!stable) sendcounts[rank] = recvcounts[rank] = 0;
  }
  
  mpi_rebalance_alltoallv(s0, sendcounts, senddispls, s1, recvcounts, recvdispls, size, rank, comm);


  if (dst_size) *dst_size = sizes[rank * 2 + 1];

  if (s1) s1->size = sizes[rank * 2 + 1];
  else s0->size = sizes[rank * 2 + 1];

  z_free(sendcounts);

free_and_exit:
  z_free(sizes);

  return exit_code;
}


slint_t mpi_rebalance_alltoallv(elements_t *sbuf, int *scounts, int *sdispls, elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_rebalance_alltoallv */
{
  slint_t i;

#ifdef NONBLOCKING
  slint_t j, nreqs;
  const slint_t request_per_element = 1 + data_n;
  MPI_Request reqs[(size - 1) * request_per_element];
  MPI_Status stats[(size - 1) * request_per_element];
#else
  MPI_Status status;
#endif
  
  if (sbuf != NULL) /* out-of-place */
  {
#ifdef LOCAL_NCOPY
    if (scounts[rank] > 0) elem_ncopy_at(sbuf, sdispls[rank], rbuf, rdispls[rank], scounts[rank]);
    scounts[rank] = rcounts[rank] = 0;
#endif

#define xelem_call \
    MPI_Alltoallv(xelem_buf(sbuf), scounts, sdispls, xelem_mpi_datatype, xelem_buf(rbuf), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  } else /* in-place */
  {
    sbuf = rbuf;

#ifdef NONBLOCKING
    /* send */
    nreqs = 0;
    for (i = 1; i < size; ++i)
    {
      j = (rank + i) % size;
      if (scounts[j] > 0)
      {
#define xelem_call \
        MPI_Isend(xelem_buf_at(sbuf, sdispls[j]), scounts[j], xelem_mpi_datatype, j, MPI_REBALANCE_ALLTOALLV_TAG, comm, &reqs[nreqs]); ++nreqs;
#include "sl_xelem_call.h"
      }
    }
    MPI_Waitall(nreqs, reqs, stats);
    
    /* local move */
    if (scounts[rank] > 0) elem_nmove_at(sbuf, sdispls[rank], sbuf, rdispls[rank], scounts[rank]);
    
    /* receive */
    nreqs = 0;
    for (i = 1; i < size; ++i)
    {
      j = (rank - i + size) % size;
      if (rcounts[j] > 0)
      {
#define xelem_call \
        MPI_Irecv(xelem_buf_at(rbuf, rdispls[j]), rcounts[j], xelem_mpi_datatype, j, MPI_REBALANCE_ALLTOALLV_TAG, comm, &reqs[nreqs]); ++nreqs;
#include "sl_xelem_call.h"
      }
    }
    MPI_Waitall(nreqs, reqs, stats);

#else /* NONBLOCKING */

    /* send right */
    for (i = rank + 1; i < size; ++i)
    if (scounts[i] > 0)
    {
#define xelem_call \
      MPI_Send(xelem_buf_at(sbuf, sdispls[i]), scounts[i], xelem_mpi_datatype, i, MPI_REBALANCE_ALLTOALLV_TAG, comm);
#include "sl_xelem_call.h"
    }

    /* move right (local) */
    if (scounts[rank] > 0 && sdispls[rank] < rdispls[rank]) elem_nmove_at(sbuf, sdispls[rank], sbuf, rdispls[rank], scounts[rank]);
    
    /* receive left */
    for (i = rank - 1; i >= 0; --i)
    if (rcounts[i] > 0)
    {
#define xelem_call \
      MPI_Recv(xelem_buf_at(sbuf, rdispls[i]), rcounts[i], xelem_mpi_datatype, i, MPI_REBALANCE_ALLTOALLV_TAG, comm, &status);
#include "sl_xelem_call.h"
    }

    /* send left */
    for (i = rank - 1; i >= 0; --i)
    if (scounts[i] > 0)
    {
#define xelem_call \
      MPI_Send(xelem_buf_at(sbuf, sdispls[i]), scounts[i], xelem_mpi_datatype, i, MPI_REBALANCE_ALLTOALLV_TAG, comm);
#include "sl_xelem_call.h"
    }

    /* move left (local) */
    if (scounts[rank] > 0 && sdispls[rank] > rdispls[rank]) elem_nmove_at(sbuf, sdispls[rank], sbuf, rdispls[rank], scounts[rank]);

    /* receive right */
    for (i = rank + 1; i < size; ++i)
    if (rcounts[i] > 0)
    {
#define xelem_call \
      MPI_Recv(xelem_buf_at(sbuf, rdispls[i]), rcounts[i], xelem_mpi_datatype, i, MPI_REBALANCE_ALLTOALLV_TAG, comm, &status);
#include "sl_xelem_call.h"
    }

#endif /* NONBLOCKING */
  }

  return 0;
}


#undef LOCAL_NCOPY
#undef NONBLOCKING
#undef MPI_REBALANCE_ALLTOALLV_TAG



/* sl_macro MSC_TRACE_IF */


#include "sl_common.h"


#ifndef MSC_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSC_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSC_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


void mpi_partcond_set_even(partcond_t *pcond, slint_t pcm, slint_t ntotal, double nimba, double wtotal, double wimba, int size, int rank) /* sl_proto, sl_func mpi_partcond_set_even */
{
  slint_t npproc = ntotal / size;
#ifdef elem_weight
  double wpproc = wtotal / size;
#endif  


  pcond->pcm = 0;

  if (pcm & SLPC_COUNTS_LH)
  {
    pcond->pcm |= SLPC_COUNTS_LH;
    pcond->count_low = npproc * (rank + 0) - (npproc * 0.5 * nimba);
    pcond->count_low = z_max(pcond->count_low, 0);
    
    pcond->count_high = npproc * (rank + 1) + (npproc * 0.5 * nimba);
    pcond->count_high = z_min(pcond->count_high, ntotal);
  }

  if (pcm & SLPC_COUNTS_MM)
  {
    pcond->pcm |= SLPC_COUNTS_MM;
    pcond->count_min = floor(npproc * (1.0 - nimba));
    pcond->count_max =  ceil(npproc * (1.0 + nimba));
  }

#ifdef elem_weight
  if (pcm & SLPC_WEIGHTS_LH)
  {
    pcond->pcm |= SLPC_WEIGHTS_LH;
    pcond->weight_low =  wpproc * (rank + 0) - (wpproc * 0.5 * wimba);
    pcond->weight_low = z_max(pcond->weight_low, 0);
    
    pcond->weight_high = wpproc * (rank + 1) + (wpproc * 0.5 * wimba);
    pcond->weight_high = z_min(pcond->weight_high, wtotal);
  }

  if (pcm & SLPC_WEIGHTS_MM)
  {
    pcond->pcm |= SLPC_WEIGHTS_MM;
    pcond->weight_min = floor(wpproc * (1.0 - wimba));
    pcond->weight_max =  ceil(wpproc * (1.0 + wimba));
  }
#endif
}


slint_t init_partconds(slint_t npconds, partcond_t *pconds, slint_t nparts, slint_t total_count, slweight_t total_weight) /* sl_proto, sl_func init_partconds */
{
  slint_t i;


  for (i = 0; i < npconds; ++i)
  {
    /* set default values and determine local (count/weight) limits */
    if (!(pconds[i].pcm & SLPC_COUNTS_MM)) { pconds[i].count_min = 0.0; pconds[i].count_max = (double) total_count; }
    if (!(pconds[i].pcm & SLPC_COUNTS_LH)) { pconds[i].count_low = 0.0; pconds[i].count_high = (double) total_count; }

    if (pconds[i].count_min < 0.0) pconds[i].count_min *= -total_count / nparts;
    if (pconds[i].count_max < 0.0) pconds[i].count_max *= -total_count / nparts;
    if (pconds[i].count_low < 0.0) pconds[i].count_low *= -total_count;
    if (pconds[i].count_high < 0.0) pconds[i].count_high *= -total_count;

/*    if (sl_mpi_rank == 0) printf("before: %e/%e - %e/%e\n", pconds[i].count_min, pconds[i].count_max, pconds[i].count_low, pconds[i].count_high);*/

    pconds[i].count_min = ceil(pconds[i].count_min);
    pconds[i].count_max = floor(pconds[i].count_max);
    pconds[i].count_low = ceil(pconds[i].count_low);
    pconds[i].count_high = floor(pconds[i].count_high);

/*    if (sl_mpi_rank == 0) printf("before: %e/%e - %e/%e\n", pconds[i].count_min, pconds[i].count_max, pconds[i].count_low, pconds[i].count_high);*/

#ifdef elem_weight
    if (!(pconds[i].pcm & SLPC_WEIGHTS_MM)) { pconds[i].weight_min = 0.0; pconds[i].weight_max = total_weight; }
    if (!(pconds[i].pcm & SLPC_WEIGHTS_LH)) { pconds[i].weight_low = 0.0; pconds[i].weight_high = total_weight; }

    if (pconds[i].weight_min < 0.0) pconds[i].weight_min *= -total_weight / nparts;
    if (pconds[i].weight_max < 0.0) pconds[i].weight_max *= -total_weight / nparts;
    if (pconds[i].weight_low < 0.0) pconds[i].weight_low *= -total_weight;
    if (pconds[i].weight_high < 0.0) pconds[i].weight_high *= -total_weight;

/*    if (sl_mpi_rank == 0) printf("before: %e/%e - %e/%e\n", pconds[i].weight_min, pconds[i].weight_max, pconds[i].weight_low, pconds[i].weight_high);*/
#endif
  }

  return 0;
}


slint_t init_partconds_intern(slint_t npconds, partcond_intern_t *pci, partcond_t *pc, slint_t nparts, slint_t total_count, slweight_t total_weight) /* sl_proto, sl_func init_partconds_intern */
{
  slint_t i;

  double avg_count = ((double) total_count) / nparts;
  slint_t sum_count[2] = { 0, 0 };
  slint_t tot_count[2] = { total_count, total_count };
#ifdef elem_weight
  double avg_weight = ((double) total_weight) / nparts;
  slweight_t sum_weight[2] = { 0.0, 0.0 };
  slweight_t tot_weight[2] = { total_weight, total_weight };
#endif

  for (i = 0; i < npconds; ++i)
  {
    pci[i].pcm = pc[i].pcm;

    Z_TRACE_IF(MSC_TRACE_IF, "IN partcond %" slint_fmt ": %" slint_fmt, i, pc[i].pcm);

    Z_TRACE_IF(MSC_TRACE_IF, "IN partcond count %" slint_fmt ": %f / %f / %f / %f", i, pc[i].count_min, pc[i].count_max, pc[i].count_low, pc[i].count_high);

    if (pci[i].pcm & SLPC_COUNTS_MM)
    {
      pci[i].count_min = (slint_t) z_round((pc[i].count_min < 0.0)?(-pc[i].count_min * avg_count):pc[i].count_min);
      pci[i].count_max = (slint_t) z_round((pc[i].count_max < 0.0)?(-pc[i].count_max * avg_count):pc[i].count_max);

      /* check min/max consistency */
      if (pci[i].count_min > pci[i].count_max) pci[i].count_min = pci[i].count_max = (pci[i].count_min + pci[i].count_max) / 2;

      sum_count[0] += pci[i].count_min;
      sum_count[1] += pci[i].count_max;

    } else { pci[i].count_min = 0; pci[i].count_max = total_count; }

    if (pci[i].pcm & SLPC_COUNTS_LH)
    {
      pci[i].count_low  = (slint_t) z_round((pc[i].count_low  < 0.0)?(-pc[i].count_low  * total_count):pc[i].count_low);
      pci[i].count_high = (slint_t) z_round((pc[i].count_high < 0.0)?(-pc[i].count_high * total_count):pc[i].count_high);
      
      /* FIXME: low/high consistency not checked */

    } else { pci[i].count_low = 0; pci[i].count_high = total_count; }

    Z_TRACE_IF(MSC_TRACE_IF, "OUT partcond count %" slint_fmt ": %" slint_fmt " / %" slint_fmt " / %" slint_fmt " / %" slint_fmt, i, pci[i].count_min, pci[i].count_max, pci[i].count_low, pci[i].count_high);

#ifdef elem_weight
    Z_TRACE_IF(MSC_TRACE_IF, "IN partcond weight %" slint_fmt ": %f / %f / %f / %f", i, pc[i].weight_min, pc[i].weight_max, pc[i].weight_low, pc[i].weight_high);

    if (pci[i].pcm & SLPC_WEIGHTS_MM)
    {
      /* round only if weight is an integral type */
      if (0.5 != (slweight_t) 0.5)
      {
        pci[i].weight_min = (slweight_t) z_round((pc[i].weight_min < 0.0)?(-pc[i].weight_min * avg_weight):pc[i].weight_min);
        pci[i].weight_max = (slweight_t) z_round((pc[i].weight_max < 0.0)?(-pc[i].weight_max * avg_weight):pc[i].weight_max);

      } else
      {
        pci[i].weight_min = (pc[i].weight_min < 0.0)?(-pc[i].weight_min * avg_weight):pc[i].weight_min;
        pci[i].weight_max = (pc[i].weight_max < 0.0)?(-pc[i].weight_max * avg_weight):pc[i].weight_max;
      }
      
      /* check min/max consistency */
      if (pci[i].weight_min > pci[i].weight_max) pci[i].weight_min = pci[i].weight_max = (pci[i].weight_min + pci[i].weight_max) / 2.0;

      sum_weight[0] += pci[i].weight_min;
      sum_weight[1] += pci[i].weight_max;

    } else { pci[i].weight_min = 0.0; pci[i].weight_max = total_weight; }

    if (pci[i].pcm & SLPC_WEIGHTS_LH)
    {
      /* round only if weight is an integral type */
      if (0.5 != (slweight_t) 0.5)
      {
        pci[i].weight_low  = (slweight_t) z_round((pc[i].weight_low  < 0.0)?(-pc[i].weight_low  * total_weight):pc[i].weight_low);
        pci[i].weight_high = (slweight_t) z_round((pc[i].weight_high < 0.0)?(-pc[i].weight_high * total_weight):pc[i].weight_high);

      } else
      {
        pci[i].weight_low  = (pc[i].weight_low  < 0.0)?(-pc[i].weight_low  * total_weight):pc[i].weight_low;
        pci[i].weight_high = (pc[i].weight_high < 0.0)?(-pc[i].weight_high * total_weight):pc[i].weight_high;
      }

      /* FIXME: low/high consistency not checked */

    } else { pci[i].weight_low = 0; pci[i].weight_high = total_weight; }

    Z_TRACE_IF(MSC_TRACE_IF, "OUT partcond weight %" slint_fmt ": %" slweight_fmt " / %" slweight_fmt " / %" slweight_fmt " / %" slweight_fmt, i, pci[i].weight_min, pci[i].weight_max, pci[i].weight_low, pci[i].weight_high);
#endif
  }

  if (!(pci[0].pcm & SLPC_COUNTS_MM)) total_count = 0;

  Z_TRACE_IF(MSC_TRACE_IF, "total_count = %" slint_fmt ", sum_count = %" slint_fmt " / %" slint_fmt, total_count, sum_count[0], sum_count[1]);

#ifdef elem_weight
  if (!(pci[0].pcm & SLPC_WEIGHTS_MM)) total_weight = 0.0;

  Z_TRACE_IF(MSC_TRACE_IF, "total_weight = %" slweight_fmt ", sum_weight = %" slweight_fmt " / %" slweight_fmt, total_weight, sum_weight[0], sum_weight[1]);
#endif

  if (sum_count[0] <= total_count) sum_count[0] = -1;
  if (sum_count[1] >= total_count) sum_count[1] = -1;
#ifdef elem_weight
  if (sum_weight[0] <= total_weight) sum_weight[0] = -1;
  if (sum_weight[1] >= total_weight) sum_weight[1] = -1;
#endif

  if (sum_count[0] > 0 || sum_count[1] > 0
#ifdef elem_weight
   || sum_weight[0] > 0.0 || sum_weight[1] > 0.0
#endif
   )
  for (i = 0; i < npconds; ++i)
  {
    Z_TRACE_IF(MSC_TRACE_IF, "%" slint_fmt ": in: count_min/max: %" slint_fmt " / %" slint_fmt, i, pci[i].count_min, pci[i].count_max);

    if (sum_count[0] > 0)
    {
      sum_count[0] -= pci[i].count_min;
      pci[i].count_min = (tot_count[0] * pci[i].count_min) / (sum_count[0] + pci[i].count_min);
      tot_count[0] -= pci[i].count_min;
    }
    if (sum_count[1] > 0)
    {
      sum_count[1] -= pci[i].count_max;
      pci[i].count_max = (tot_count[1] * pci[i].count_max) / (sum_count[1] + pci[i].count_max);
      tot_count[1] -= pci[i].count_max;
    }

    Z_TRACE_IF(MSC_TRACE_IF, "%" slint_fmt ": out: count_min/max: %" slint_fmt " / %" slint_fmt, i, pci[i].count_min, pci[i].count_max);

#ifdef elem_weight
    if (sum_weight[0] > 0)
    {
      sum_weight[0] -= pci[i].weight_min;
      pci[i].weight_min = (tot_weight[0] * pci[i].weight_min) / (sum_weight[0] + pci[i].weight_min);
      tot_weight[0] -= pci[i].weight_min;
    }
    if (sum_weight[1] > 0)
    {
      sum_weight[1] -= pci[i].weight_max;
      pci[i].weight_max = (tot_weight[1] * pci[i].weight_max) / (sum_weight[1] + pci[i].weight_max);
      tot_weight[1] -= pci[i].weight_max;
    }
#endif
  }

  return 0;
}


slint_t merge_partconds(partcond_t *pconds_in, slint_t npconds_in, partcond_t *pcond_out)  /* sl_proto, sl_func merge_partconds */
{
  slint_t i;

  pcond_out->pcm = pconds_in[0].pcm;

  pcond_out->count_min = 0;
  pcond_out->count_max = 0;
  pcond_out->count_low = pconds_in[0].count_low;
  pcond_out->count_high = pconds_in[npconds_in - 1].count_high;
  for (i = 0; i < npconds_in; ++i)
  {
    pcond_out->count_min += pconds_in[i].count_min;
    pcond_out->count_max += pconds_in[i].count_max;
  }

#ifdef elem_weight  
  pcond_out->weight_min = 0;
  pcond_out->weight_max = 0;
  pcond_out->weight_low = pconds_in[0].weight_low;
  pcond_out->weight_high = pconds_in[npconds_in - 1].weight_high;
  for (i = 0; i < npconds_in; ++i)
  {
    pcond_out->weight_min += pconds_in[i].weight_min;
    pcond_out->weight_max += pconds_in[i].weight_max;
  }
#endif

  return 0;
}


slint_t mpi_gather_partconds_grouped(partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, partcond_t *pconds_out, slint_t *npconds_out, int size, int rank, MPI_Comm comm)  /* sl_proto, sl_func mpi_gather_partconds_grouped */
{
  int _npconds_out = 1;

  if (npconds_out) _npconds_out = *npconds_out;

  if (_npconds_out < 0)
  {
    if (pcond_in_comm != MPI_COMM_NULL) MPI_Comm_size(pcond_in_comm, &_npconds_out);
    if (pconds_out_comm != MPI_COMM_NULL) MPI_Bcast(&_npconds_out, 1, MPI_INT, 0, pconds_out_comm);
  }
  
  if (pcond_in_comm != MPI_COMM_NULL)
  {
    if (pconds_out)
    {
      if (pcond_in) MPI_Allgather(pcond_in, sizeof(partcond_t), MPI_BYTE, pconds_out, sizeof(partcond_t), MPI_BYTE, pcond_in_comm);
    }

  }

  if (pconds_out_comm != MPI_COMM_NULL && pconds_out) MPI_Bcast(pconds_out, _npconds_out * sizeof(partcond_t), MPI_BYTE, 0, pconds_out_comm);

  if (npconds_out) *npconds_out = _npconds_out;

  return 0;
}


slint_t mpi_gather_partconds(partcond_t *pcond_in, partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_gather_partconds */
{
/*  if (comm != MPI_COMM_NULL)*/
  MPI_Gather(pcond_in, sizeof(partcond_t), MPI_BYTE, pconds_out, sizeof(partcond_t), MPI_BYTE, root, comm);

  return 0;
}


slint_t mpi_allgather_partconds(partcond_t *pcond_in, partcond_t *pconds_out, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_allgather_partconds */
{
/*  if (comm != MPI_COMM_NULL)*/
  MPI_Allgather(pcond_in, sizeof(partcond_t), MPI_BYTE, pconds_out, sizeof(partcond_t), MPI_BYTE, comm);

  return 0;
}


slint_t mpi_bcast_partconds(slint_t npconds, partcond_t *pconds, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_bcast_partconds */
{
/*  if (comm != MPI_COMM_NULL)*/
  MPI_Bcast(pconds, npconds * sizeof(partcond_t), MPI_BYTE, root, comm);

  return 0;
}


slint_t mpi_post_check_partconds(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_post_check_partconds */
{
  slint_t i, j;
  int ssdispls[nparts], sscounts[nparts], rrcounts[nparts], ss;

  for (i = 0; i < nparts; ++i)
  {
    ssdispls[i] = 0;
    for (j = 0; j < nelements; ++j) ssdispls[i] += sdispls[i * nelements + j];
  }

  ss = 0;
  for (j = 0; j < nelements; ++j) ss += s[j].size;
  displs2counts(nparts, ssdispls, sscounts, ss);

  MPI_Reduce(sscounts, rrcounts, nparts,  MPI_INT, MPI_SUM, 0, comm);

  j = -1;
  if (rank == 0)
  for (i = 0; i < nparts; ++i)
  {
    Z_TRACE_IF(MSC_TRACE_IF, "%" slint_fmt " verifying %d against [%f  %f]", i, rrcounts[i], pconds[i].count_min, pconds[i].count_max);
  
    if (rrcounts[i] < pconds[i].count_min || rrcounts[i] > pconds[i].count_max) j = i;
  }

  MPI_Bcast(&j, 1, int_mpi_datatype, 0, comm);
  
  return j;
}


slint_t mpi_post_check_partconds_intern(elements_t *s, slint_t nelements, slint_t nparts, partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_post_check_partconds_intern */
{
  slint_t i, j;
  int ssdispls[nparts], sscounts[nparts], rrcounts[nparts], ss;

  for (i = 0; i < nparts; ++i)
  {
    ssdispls[i] = 0;
    for (j = 0; j < nelements; ++j) ssdispls[i] += sdispls[i * nelements + j];
  }

  ss = 0;
  for (j = 0; j < nelements; ++j) ss += s[j].size;
  displs2counts(nparts, ssdispls, sscounts, ss);

  MPI_Reduce(sscounts, rrcounts, nparts,  MPI_INT, MPI_SUM, 0, comm);

  j = -1;
  if (rank == 0)
  for (i = 0; i < nparts; ++i)
  {
    Z_TRACE_IF(MSC_TRACE_IF, "%" slint_fmt " verifying %d against [%" slint_fmt "  %" slint_fmt "]", i, rrcounts[i], pci[i].count_min, pci[i].count_max);
  
    if (rrcounts[i] < pci[i].count_min || rrcounts[i] > pci[i].count_max) j = i;
  }

  MPI_Bcast(&j, 1, int_mpi_datatype, 0, comm);
  
  return j;
}


slint_t mpi_select_stats(elements_t *s, slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_stats */
{
  slint_t i;
#ifdef elem_weight
  slint_t j;
#endif

  double v, vmin, vmax;

  slint_t partial_counts[nparts + 1];
#ifndef HAVE_MPI_IN_PLACE
  slint_t partial_counts_tmp[nparts + 1];
#endif

#ifdef elem_weight
  slweight_t partial_weights[nparts + 1];
# ifndef HAVE_MPI_IN_PLACE
  slweight_t partial_weights_tmp[nparts + 1];
# endif
#endif


  partial_counts[nparts] = 0;
#ifdef elem_weight
  partial_weights[nparts] = 0.0;
#endif

  for (i = 0; i < nparts; ++i)
  {
    partial_counts[i] = ((i < nparts - 1)?sdispls[i + 1]:s->size) - sdispls[i];
    partial_counts[nparts] += partial_counts[i];

#ifdef elem_weight
    partial_weights[i] = 0.0;
    for (j = sdispls[i]; j < ((i < nparts - 1)?sdispls[i + 1]:s->size); ++j) partial_weights[i] += elem_weight(s, j);
    partial_weights[nparts] += partial_weights[i];
#endif
  }

#ifdef HAVE_MPI_IN_PLACE
  /* recvbuf requires workaround for an in-place/aliased-buffer-check-bug in mpich2 (fixed with rev 5518) */
  MPI_Reduce((rank == 0)?MPI_IN_PLACE:partial_counts, (rank == 0)?partial_counts:NULL, nparts + 1, int_mpi_datatype, MPI_SUM, 0, comm);
# ifdef elem_weight
  MPI_Reduce((rank == 0)?MPI_IN_PLACE:partial_weights, (rank == 0)?partial_weights:NULL, nparts + 1, weight_mpi_datatype, MPI_SUM, 0, comm);
# endif
#else
  MPI_Reduce(partial_counts, partial_counts_tmp, nparts + 1, int_mpi_datatype, MPI_SUM, 0, comm);
# define partial_counts   partial_counts_tmp
# ifdef elem_weight
  MPI_Reduce(partial_weights, partial_weights_tmp, nparts + 1, weight_mpi_datatype, MPI_SUM, 0, comm);
# define partial_weights  partial_weights_tmp
# endif
#endif

  if (rank == 0)
  {
    printf("%d: count total: %" sl_int_type_fmt "\n", rank, partial_counts[nparts]);
    v = 0.0;
    vmin = 1.0;
    vmax = 0.0;
    for (i = 0; i < nparts; ++i)
    {
/*      printf("%d: %" sl_int_type_fmt " %" sl_int_type_fmt " / %f - %" sl_int_type_fmt " / %f\n", rank, i, partial_counts[i], (double) partial_counts[i] / partial_counts[nparts], (partial_counts[nparts] / nparts) - partial_counts[i], fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts])));*/
      v += fabs(((double) partial_counts[nparts] / nparts) - partial_counts[i]);
      if (fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts])) < vmin) vmin = fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts]));
      if (fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts])) > vmax) vmax = fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts]));
    }
    printf("%d: count min/max: %f / %f\n", rank, vmin, vmax);
    printf("%d: count average: %f - %f / %f\n", rank, (double) partial_counts[nparts] / nparts, v / nparts, v / partial_counts[nparts]);

#ifdef elem_weight
    printf("%d: weight total: %" slweight_fmt "\n", rank, partial_weights[nparts]);
    v = 0.0;
    vmin = 1.0;
    vmax = 0.0;
    for (i = 0; i < nparts; ++i)
    {
/*      printf("%d: %" sl_int_type_fmt " %f / %f - %f / %f\n", rank, i, partial_weights[i], partial_weights[i] / partial_weights[nparts], (partial_weights[nparts] / nparts) - partial_weights[i], fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts])));*/
      v += fabs((partial_weights[nparts] / nparts) - partial_weights[i]);
      if (fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts])) < vmin) vmin = fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts]));
      if (fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts])) > vmax) vmax = fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts]));
    }
    printf("%d: weight min/max: %f / %f\n", rank, vmin, vmax);
    printf("%d: weight average: %f - %f / %f\n", rank, (double) (partial_weights[nparts] / nparts), v / nparts, (double) (v / partial_weights[nparts]));
#endif
  }
  
  return 0;

#undef partial_counts
#undef partial_weights
}



/* sl_macro_global MSEG_ROOT */
/* sl_macro_global MSEG_BORDER_UPDATE_REDUCTION */
/* sl_macro_global MSEG_DISABLE_BEST_CHOICE */
/* sl_macro_global MSEG_DISABLE_MINMAX */
/* sl_macro_global MSEG_ENABLE_OPTIMZED_LOWHIGH */
/* sl_macro_global MSEG_FORWARD_ONLY */
/* sl_macro_global MSEG_INFO */
/* sl_macro_global MSEG_TRACE_IF */


#include "sl_common.h"


/* config */
/*#define SYNC_ON_INIT
#define SYNC_ON_EXIT*/

/*#define PRINT_SDISPLS*/
/*#define PRINT_STATS*/
/*#define PRINT_TIMINGS  0*/

/*#define VERIFY*/

typedef struct _border_info_t {
  slint_t crange[2], ccurrent[2];
#ifdef elem_weight
  slweight_t wrange[2], wcurrent[2];
#endif

} border_info_t;

#define LO  0
#define HI  1

/* sl_ifdef SL_USE_MPI sl_context CONTEXT_BEGIN mseg */
#ifdef MSEG_ROOT
const int default_mseg_root = -1;  /* sl_ifdef MSEG_ROOT sl_endif sl_global sl_context sl_var default_mseg_root */
#endif

#ifdef MSEG_BORDER_UPDATE_REDUCTION
const double default_mseg_border_update_count_reduction = 0.0;  /* sl_ifdef MSEG_BORDER_UPDATE_REDUCTION sl_global sl_context sl_var default_mseg_border_update_count_reduction */
# ifdef elem_weight
const double default_mseg_border_update_weight_reduction = 0.0;  /* sl_endif sl_global sl_context sl_var default_mseg_border_update_weight_reduction */
# endif
#endif

#ifdef MSEG_FORWARD_ONLY
const slint_t default_mseg_forward_only = 0;  /* sl_ifdef MSEG_FORWARD_ONLY sl_endif sl_global sl_context sl_var default_mseg_forward_only */
#endif

#ifdef MSEG_INFO
const slint_t default_mseg_info_rounds = 0;             /* sl_ifdef MSEG_INFO sl_global sl_context sl_var default_mseg_info_rounds */
const slint_t *default_mseg_info_finish_rounds = NULL;  /* sl_global sl_context sl_var default_mseg_info_finish_rounds */
const double default_mseg_info_finish_rounds_avg = 0;   /* sl_global sl_context sl_var default_mseg_info_finish_rounds_avg */
# ifdef elem_weight
const slweight_t default_mseg_info_total_weights = 0;   /* sl_endif sl_global sl_context sl_var default_mseg_info_total_weights */
# endif
#endif

const slint_t default_mseg_binnings = -1;  /* sl_global sl_context sl_var default_mseg_binnings */

const slint_t default_mseg_finalize_mode = SL_MSEG_FM_EXACT;  /* sl_global sl_context sl_var default_mseg_finalize_mode */
/* sl_endif sl_context CONTEXT_END mseg */

#ifndef MSEG_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSEG_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSEG_TRACE_IF  (SL_PROC_RANK == -1)
# endif
#endif

#define MSEG_ASSERT_IF  (rank == 0)


static void border_init(slint_t docounts, slint_t doweights, border_info_t *bi, slint_t current, slint_t tc, slweight_t tw)
{
  if (docounts)
  {
    bi[0].crange[0] = 0;
    bi[0].crange[1] = tc;

    Z_TRACE_IF(MSEG_TRACE_IF, "count range: %" slint_fmt " - %" slint_fmt "",
      bi[0].crange[0], bi[0].crange[1]);

    bi[0].ccurrent[LO] = (current == 0 || current != 1)?0:tc;
    bi[0].ccurrent[HI] = (current == 1 || current != 0)?tc:0;

    Z_TRACE_IF(MSEG_TRACE_IF, "count[low/high]: %" slint_fmt " / %" slint_fmt,
      bi[0].ccurrent[LO], bi[0].ccurrent[HI]);
  }

#ifdef elem_weight
  if (doweights)
  {
    bi[0].wrange[0] = 0;
    bi[0].wrange[1] = tw;

    Z_TRACE_IF(MSEG_TRACE_IF, "weight range: %" slweight_fmt " - %" slweight_fmt,
      bi[0].wrange[0], bi[0].wrange[1]);

    bi[0].wcurrent[LO] = (current == 0 || current != 1)?0:tw;
    bi[0].wcurrent[HI] = (current == 1 || current != 0)?tw:0;

    Z_TRACE_IF(MSEG_TRACE_IF, "weight[low/high]: %" slweight_fmt " / %" slweight_fmt,
      bi[0].wcurrent[LO], bi[0].wcurrent[HI]);

  }
#endif
}


static void border_update(slint_t docounts, slint_t doweights, border_info_t *bi, partcond_intern_t *pc, slint_t dir, slint_t reduction)
{
  slint_t ccurrent_new[2] = { -1, -1 };
#ifdef elem_weight
  slweight_t wcurrent_new[2] = { -1, -1 };
#endif

#ifdef MSEG_BORDER_UPDATE_REDUCTION
  slint_t count_reduction;
# ifdef elem_weight
  slweight_t weight_reduction;
# endif
#endif


#if 1
  /* if counts are used and the current interval is already as small as possible, then skip the update */
  if (docounts && bi[0].ccurrent[LO] == bi[0].ccurrent[HI])
  {
    Z_TRACE_IF(MSEG_TRACE_IF, "skip border_update");
    return;
  }
#endif

  /* init from range */
  if (docounts)
  {
    ccurrent_new[LO] = bi[0].crange[0];
    ccurrent_new[HI] = bi[0].crange[1];
  }
#ifdef elem_weight
  if (doweights)
  {
    wcurrent_new[LO] = bi[0].wrange[0];
    wcurrent_new[HI] = bi[0].wrange[1];
  }
#endif

  /* forward */
  if (dir > 0)
  {
    if (docounts)
    {
      /* init from min/max */
      if (pc[0].pcm & SLPC_COUNTS_MM)
      {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
        if (reduction)
        {
          count_reduction = z_round((bi[-1].ccurrent[HI] - bi[-1].ccurrent[LO]) * 0.5 * SL_DEFCON(mseg.border_update_count_reduction));
          Z_ASSERT(count_reduction >= 0);

          Z_TRACE_IF(MSEG_TRACE_IF, "forward: count_reduction: %" slint_fmt, count_reduction);

          ccurrent_new[LO] = z_min(bi[-1].ccurrent[LO] + count_reduction + pc[0].count_min, bi[0].ccurrent[HI]);
          ccurrent_new[HI] = z_max(bi[-1].ccurrent[HI] - count_reduction + pc[0].count_max, bi[0].ccurrent[LO]);

          Z_TRACE_IF(MSEG_TRACE_IF, "forward: count[low/high]: min(%" slint_fmt " + %" slint_fmt " + %" slint_fmt ", %" slint_fmt "), max(%" slint_fmt " - % "slint_fmt " + %" slint_fmt ", %" slint_fmt ")",
            bi[-1].ccurrent[LO], count_reduction, pc[0].count_min, bi[0].ccurrent[HI], bi[-1].ccurrent[HI], count_reduction, pc[0].count_max, bi[0].ccurrent[LO]);

        } else
#endif
        {
          ccurrent_new[LO] = bi[-1].ccurrent[LO] + pc[0].count_min;
          ccurrent_new[HI] = bi[-1].ccurrent[HI] + pc[0].count_max;

          Z_TRACE_IF(MSEG_TRACE_IF, "forward: count[low/high]: %" slint_fmt " + %" slint_fmt ", %" slint_fmt " + %" slint_fmt,
            bi[-1].ccurrent[LO], pc[0].count_min, bi[-1].ccurrent[HI], pc[0].count_max);
        }
      }
    }

#ifdef elem_weight
    if (doweights)
    {
      /* init from min/max */
      if (pc[0].pcm & SLPC_WEIGHTS_MM)
      {
# ifdef MSEG_BORDER_UPDATE_REDUCTION
        if (reduction)
        {
          weight_reduction = (bi[-1].wcurrent[HI] - bi[-1].wcurrent[LO]) * 0.5 * SL_DEFCON(mseg.border_update_weight_reduction);
          Z_ASSERT(weight_reduction >= 0.0);

          Z_TRACE_IF(MSEG_TRACE_IF, "forward: weight_reduction: %" slweight_fmt, weight_reduction);
  
          wcurrent_new[LO] = z_min(bi[-1].wcurrent[LO] + weight_reduction + pc[0].weight_min, bi[0].wcurrent[HI]);
          wcurrent_new[HI] = z_max(bi[-1].wcurrent[HI] - weight_reduction + pc[0].weight_max, bi[0].wcurrent[LO]);

          Z_TRACE_IF(MSEG_TRACE_IF, "weight_reduction: weight[low/high]: min(%" slweight_fmt " + %" slweight_fmt " + %" slweight_fmt ", %" slweight_fmt "), max(%" slweight_fmt " - %" slweight_fmt " + %" slweight_fmt ", %" slweight_fmt ")",
            bi[-1].wcurrent[LO], weight_reduction, pc[0].weight_min, bi[0].wcurrent[HI], bi[-1].wcurrent[HI], weight_reduction, pc[0].weight_max, bi[0].wcurrent[LO]);

        } else
# endif
        {
          wcurrent_new[LO] = bi[-1].wcurrent[LO] + pc[0].weight_min;
          wcurrent_new[HI] = bi[-1].wcurrent[HI] + pc[0].weight_max;

          Z_TRACE_IF(MSEG_TRACE_IF, "forward: weight[low/high]: %" slweight_fmt " + %" slweight_fmt ", %" slweight_fmt " + %" slweight_fmt,
            bi[-1].wcurrent[LO], pc[0].weight_min, bi[-1].wcurrent[HI], pc[0].weight_max);
        }
      }

    } else
#endif
    {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
      if (reduction) Z_TRACE_IF(MSEG_TRACE_IF, "");
#endif
      Z_TRACE_IF(MSEG_TRACE_IF, "");
    }

  } else /* backward */
  {
    if (docounts)
    {
      /* init from min/max */
      if (pc[0].pcm & SLPC_COUNTS_MM)
      {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
        if (reduction)
        {
          count_reduction = z_round((bi[1].ccurrent[HI] - bi[1].ccurrent[LO]) * 0.5 * SL_DEFCON(mseg.border_update_count_reduction));
          Z_ASSERT(count_reduction >= 0);

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: count_reduction: %" slint_fmt, count_reduction);

          ccurrent_new[LO] = z_min(bi[1].ccurrent[LO] + count_reduction - pc[1].count_max, bi[0].ccurrent[HI]);
          ccurrent_new[HI] = z_max(bi[1].ccurrent[HI] - count_reduction - pc[1].count_min, bi[0].ccurrent[LO]);

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: count[low/high]: min(%" slint_fmt " + %" slint_fmt " - %" slint_fmt ", %" slint_fmt "), max(%" slint_fmt " - %" slint_fmt " + %" slint_fmt ", %" slint_fmt ")",
            bi[1].ccurrent[LO], count_reduction, pc[1].count_max, bi[0].ccurrent[HI], bi[1].ccurrent[HI], count_reduction, pc[1].count_min, bi[0].ccurrent[LO]);

        } else
#endif
        {
          ccurrent_new[LO] = bi[1].ccurrent[LO] - pc[1].count_max;
          ccurrent_new[HI] = bi[1].ccurrent[HI] - pc[1].count_min;

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: count[low/high]: %" slint_fmt " - %" slint_fmt ", %" slint_fmt " - %" slint_fmt "",
            bi[1].ccurrent[LO], pc[1].count_max, bi[1].ccurrent[HI], pc[1].count_min);
        }
      }
    }

#ifdef elem_weight
    if (doweights)
    {
      /* init from min/max */
      if (pc[0].pcm & SLPC_WEIGHTS_MM)
      {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
        if (reduction)
        {
          weight_reduction = (bi[1].wcurrent[1] - bi[1].wcurrent[LO]) * 0.5 * SL_DEFCON(mseg.border_update_weight_reduction);
          Z_ASSERT(weight_reduction >= 0.0);

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: weight_reduction: %" slweight_fmt, weight_reduction);

          wcurrent_new[LO] = z_min(bi[1].wcurrent[LO] + weight_reduction - pc[1].weight_max, bi[0].wcurrent[HI]);
          wcurrent_new[HI] = z_max(bi[1].wcurrent[HI] - weight_reduction - pc[1].weight_min, bi[0].wcurrent[LO]);

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: weight[low/high]: min(%" slweight_fmt " + %" slweight_fmt " - %" slweight_fmt ", %" slweight_fmt "), max(%" slweight_fmt " - %" slweight_fmt " - %" slweight_fmt ", %" slweight_fmt ")",
            bi[1].wcurrent[LO], weight_reduction, pc[1].weight_max, bi[0].wcurrent[HI], bi[1].wcurrent[HI], weight_reduction, pc[1].weight_min, bi[0].wcurrent[LO]);

        } else
#endif
        {
          wcurrent_new[LO] = bi[1].wcurrent[LO] - pc[1].weight_max;
          wcurrent_new[HI] = bi[1].wcurrent[HI] - pc[1].weight_min;

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: weight[low/high]: %" slweight_fmt " - %" slweight_fmt ", %" slweight_fmt " - %" slweight_fmt,
            bi[1].wcurrent[LO], pc[1].weight_max, bi[1].wcurrent[HI], pc[1].weight_min);
        }
      }

    } else
#endif
    {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
      if (reduction) Z_TRACE_IF(MSEG_TRACE_IF, "");
#endif
      Z_TRACE_IF(MSEG_TRACE_IF, "");
    }
  }

  /* check against low/high */
  if (docounts)
  {
    if (pc[0].pcm & SLPC_COUNTS_LH)
    {
      if (ccurrent_new[LO] < pc[1].count_low)  ccurrent_new[LO] = pc[1].count_low;
      if (ccurrent_new[HI] > pc[0].count_high) ccurrent_new[HI] = pc[0].count_high;
    }

    /* fit to range (NOT REQUIRED!?) */
/*    ccurrent_new[LO] = z_minmax(bi[0].crange[0], ccurrent_new[LO], bi[0].crange[1]);
    ccurrent_new[HI] = z_minmax(bi[0].crange[0], ccurrent_new[HI], bi[0].crange[1]);*/

    bi[0].ccurrent[LO] = z_max(bi[0].ccurrent[LO], ccurrent_new[LO]);
    bi[0].ccurrent[HI] = z_min(bi[0].ccurrent[HI], ccurrent_new[HI]);

    Z_TRACE_IF(MSEG_TRACE_IF, "current count[low/high]: %" slint_fmt " / %" slint_fmt, bi[0].ccurrent[LO], bi[0].ccurrent[HI]);
  }

#ifdef elem_weight
  if (doweights)
  {
    if (pc[0].pcm & SLPC_WEIGHTS_LH)
    {
      if (wcurrent_new[LO] < pc[1].weight_low)  wcurrent_new[LO] = pc[1].weight_low;
      if (wcurrent_new[HI] > pc[0].weight_high) wcurrent_new[HI] = pc[0].weight_high;
    }

    /* fit to range (NOT REQUIRED!?) */
/*    wcurrent_new[LO] = z_minmax(bi[0].wrange[0], wcurrent_new[LO], bi[0].wrange[1]);
    wcurrent_new[HI] = z_minmax(bi[0].wrange[0], wcurrent_new[HI], bi[0].wrange[1]);*/

    bi[0].wcurrent[LO] = z_max(bi[0].wcurrent[LO], wcurrent_new[LO]);
    bi[0].wcurrent[HI] = z_min(bi[0].wcurrent[HI], wcurrent_new[HI]);

    Z_TRACE_IF(MSEG_TRACE_IF, "current weight[low/high]: %" slweight_fmt " / %" slweight_fmt, bi[0].wcurrent[LO], bi[0].wcurrent[HI]);
  }
#endif
}


#ifdef MSEG_ENABLE_OPTIMZED_LOWHIGH

# define border_update_update(_dc_, _dw_, _bi_, _pc_, _dir_, _red_)        Z_NOP()
# ifdef elem_weight
# define border_change_change(_dc_, _dw_, _bi_, _gcs_, _gc_, _gws_, _gw_)  Z_MOP( \
  if (_dc_) { (_bi_)[0].crange[0] += (_gcs_); (_bi_)[0].crange[1] = (_bi_)[0].crange[0] + (_gc_); } \
  if (_dw_) { (_bi_)[0].wrange[0] += (_gws_); (_bi_)[0].wrange[1] = (_bi_)[0].wrange[0] + (_gw_); })
# else
# define border_change_change(_dc_, _dw_, _bi_, _gcs_, _gc_, _gws_, _gw_)  Z_MOP( \
  if (_dc_) { (_bi_)[0].crange[0] += (_gcs_); (_bi_)[0].crange[1] = (_bi_)[0].crange[0] + (_gc_); })
# endif

#else /* MSEG_ENABLE_OPTIMZED_LOWHIGH */

static void border_change(slint_t docounts, slint_t doweights, border_info_t *bi, slint_t gcs, slint_t gc, slweight_t gws, slweight_t gw)
{
  if (docounts)
  {
    Z_TRACE_IF(MSEG_TRACE_IF, "change: gcs = %" slint_fmt ", gc = %" slint_fmt "", gcs, gc);

    bi[0].crange[0] += gcs;
    bi[0].crange[1] = bi[0].crange[0] + gc;

    Z_TRACE_IF(MSEG_TRACE_IF, "counts_range: %" slint_fmt "  %" slint_fmt "", bi[0].crange[0], bi[0].crange[1]);

    bi[0].ccurrent[LO] = z_minmax(bi[0].crange[0], bi[0].ccurrent[LO], bi[0].crange[1]);
    bi[0].ccurrent[HI] = z_minmax(bi[0].crange[0], bi[0].ccurrent[HI], bi[0].crange[1]);

    Z_TRACE_IF(MSEG_TRACE_IF, "count[low/high]: %" slint_fmt " / %" slint_fmt,
      bi[0].ccurrent[LO], bi[0].ccurrent[HI]);
  }

#ifdef elem_weight
  if (doweights)
  {
    Z_TRACE_IF(MSEG_TRACE_IF, "change: gws = %" slweight_fmt ", gc = %" slweight_fmt, gws, gw);

    bi[0].wrange[0] += gws;
    bi[0].wrange[1] = bi[0].wrange[0] + gw;

    Z_TRACE_IF(MSEG_TRACE_IF, "weights_range: %" slweight_fmt "  %" slweight_fmt, bi[0].wrange[0], bi[0].wrange[1]);

    bi[0].wcurrent[LO] = z_minmax(bi[0].wrange[0], bi[0].wcurrent[LO], bi[0].wrange[1]);
    bi[0].wcurrent[HI] = z_minmax(bi[0].wrange[0], bi[0].wcurrent[HI], bi[0].wrange[1]);

    Z_TRACE_IF(MSEG_TRACE_IF, "weight[low/high]: %" slweight_fmt " / %" slweight_fmt,
      bi[0].wcurrent[LO], bi[0].wcurrent[HI]);

  } else
#endif
  { Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); }

#if 0
  Z_TRACE_IF(MSEG_TRACE_IF, "range diff 0: %" slint_fmt "-%" slint_fmt " | %" slint_fmt "-%" slint_fmt,
    bi[0].crange[0] - bi[-1].crange[1], bi[0].crange[0] - bi[-1].crange[0],
    bi[1].crange[0] - bi[ 0].crange[0], bi[1].crange[1] - bi[ 0].crange[0]);
  Z_TRACE_IF(MSEG_TRACE_IF, "range diff 1: %" slint_fmt "-%" slint_fmt " | %" slint_fmt "-%" slint_fmt,
    bi[0].crange[1] - bi[-1].crange[1], bi[0].crange[1] - bi[-1].crange[0],
    bi[1].crange[0] - bi[ 0].crange[1], bi[1].crange[1] - bi[ 0].crange[1]);
#endif
}

# define border_update_update(_dc_, _dw_, _bi_, _pc_, _dir_, _red_)        border_update(_dc_, _dw_, _bi_, _pc_, _dir_, _red_)
# define border_change_change(_dc_, _dw_, _bi_, _gcs_, _gc_, _gws_, _gw_)  border_change(_dc_, _dw_, _bi_, _gcs_, _gc_, _gws_, _gw_)

#endif /* MSEG_ENABLE_OPTIMZED_LOWHIGH */


slint_t mpi_select_exact_generic_bulk(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_generic_bulk */
{
  const slint_t max_nborders = nparts - 1;
  slint_t border_lo, border_hi, nborders_removed;
  slint_t borders[max_nborders], border_bins[max_nborders];

  border_info_t border_infos_[1 + max_nborders + 1], *border_infos = border_infos_ + 1;

  slint_t total_counts;
#ifdef elem_weight
  slweight_t total_weights;
#endif

  slint_t pcm;
  partcond_intern_t pci[nparts];

#if defined(elem_weight) && defined(sl_weight_intequiv)
  slweight_t current_cw[4];
# define current_clo  current_cw[0]
# define current_chi  current_cw[1]
# define current_wlo  current_cw[2]
# define current_whi  current_cw[3]
#else
  slint_t current_c[2];
# define current_clo  current_c[0]
# define current_chi  current_c[1]
# ifdef elem_weight
  slweight_t current_w[2];
#  define current_wlo  current_w[0]
#  define current_whi  current_w[1]
# endif
#endif

  slint_t round, direction, refine, finalize, do_binning_bins, do_binning_prefixes;

  slint_t nothing, lc_min = -1, lc_max = -1;
#ifdef elem_weight
  slint_t lcw2gcw;
#endif

#if defined(elem_weight) && defined(sl_weight_intequiv)
  slweight_t final_lcws[2], final_gcws[2];
# define final_lcs  final_lcws[0]
# define final_gcs  final_gcws[0]
# define final_lws  final_lcws[1]
# define final_gws  final_gcws[1]
#else
  slint_t final_lcs, final_gcs;
# ifdef elem_weight
  slweight_t final_lws, final_gws;
# endif
#endif

  slint_t gc = 0, gcs = 0;
  slint_t final_mc = -1, final_dc = -1;
#ifdef elem_weight
  slweight_t gw = 0, gws = 0;
  slweight_t final_mw = -1, final_dw = -1;
#endif

  slint_t i, j, k, ix;

  slint_t docounts;
#ifdef elem_weight
  slint_t doweights;
# ifdef sl_weight_intequiv
  slint_t weight_factor;
# endif
#else
# define doweights  0
#endif

#ifdef VERIFY
  slint_t v;
#endif

  global_bins_t gb;


  Z_TRACE_IF(MSEG_TRACE_IF, "starting mpi_select_exact_generic");

  /* sl_tid rti_tid_mpi_select_exact_generic rti_tid_mpi_select_exact_generic_sync_init rti_tid_mpi_select_exact_generic_sync_exit */

  rti_treset(rti_tid_mpi_select_exact_generic_while);                    /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check);              /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_bins);         /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_bins_local);   /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_bins_global);  /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_round1);       /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_pre);          /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_part);         /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_part_root);    /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_final);        /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_final_root);   /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_post);         /* sl_tid */

  rti_tstart(rti_tid_mpi_select_exact_generic);

  rti_tstart(rti_tid_mpi_select_exact_generic_sync_init);
#ifdef SYNC_ON_INIT
  MPI_Barrier(comm);
#endif
  rti_tstop(rti_tid_mpi_select_exact_generic_sync_init);

#ifdef VERIFY
  v = elements_validate_order(s, 1);
  
  Z_TRACE_IF(MSEG_TRACE_IF, "elements order: %s (%" slint_fmt ")", (v > 0)?"FAILED":"SUCCESS", v);
#endif

  pcm = pconds->pcm;
#ifdef MSEG_DISABLE_MINMAX
  pcm &= ~(SLPC_COUNTS_MM|SLPC_WEIGHTS_MM);
#endif

  docounts = ((pcm & (SLPC_COUNTS_MM|SLPC_COUNTS_LH)) != 0);
#ifdef elem_weight
  doweights = ((pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
# ifdef sl_weight_intequiv
  weight_factor = 1 + (doweights != 0);
# endif
#endif

  mpi_binning_create(&gb, max_nborders, SL_DEFCON(mseg.binnings), s, nelements, docounts, doweights, bm, size, rank, comm);

  /* init parts */
  border_lo = 0;
  border_hi = max_nborders - 1;
  for (i = border_lo; i <= border_hi; ++i)
  {
    borders[i] = i;
    border_bins[i] = 0;
  }

  /* reset splitter */
  sp->n = nparts * nelements;
  splitter_reset(sp);

#ifdef MSEG_INFO
  SL_DEFCON(mseg.info_finish_rounds_avg) = 0;
#endif

  rti_tstart(rti_tid_mpi_select_exact_generic_while);

  direction = 1;

  round = 0;
  while (border_lo <= border_hi && (docounts || doweights))
  {
    ++round;

    Z_TRACE_IF(MSEG_TRACE_IF, "ROUND: %" slint_fmt ", %s, %" slint_fmt " border(s)", round, (direction > 0)?"forward":"backward", border_hi - border_lo + 1);

    nborders_removed = 0;

    mpi_binning_pre(&gb, size, rank, comm);

    Z_TRACE_IF(MSEG_TRACE_IF, "ROUND: %" slint_fmt ", bm_nbins: %" slint_fmt, round, gb.bm->nbins);

    finalize = (gb.bm->nbins <= 1);

    rti_tstart(rti_tid_mpi_select_exact_generic_while_check);

    i = (direction > 0)?border_lo:border_hi;
    while ((direction > 0)?(i <= border_hi):(i >= border_lo))
    {
      Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ": PART: %" slint_fmt ",%" slint_fmt ": %s", round, i, borders[i], ((direction > 0)?"forward":"backward"));

      ix = i;

      rti_tstart(rti_tid_mpi_select_exact_generic_while_check_bins);

      do_binning_bins = (!finalize || (finalize && round == 1 && i == border_lo));
      do_binning_prefixes = (finalize && SL_DEFCON(mseg.finalize_mode) == SL_MSEG_FM_EXACT);

      if (do_binning_bins || do_binning_prefixes)
      {
        mpi_binning_exec_reset(&gb, do_binning_bins, do_binning_prefixes, size, rank, comm);

        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_bins_local);

        while ((direction > 0)?(ix <= border_hi):(ix >= border_lo))
        {
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": bin %" slint_fmt, ix, borders[ix], border_bins[ix]);

          if (mpi_binning_exec_local(&gb, border_bins[ix], do_binning_bins, do_binning_prefixes, size, rank, comm) < 0)
          {
            Z_TRACE_IF(MSEG_TRACE_IF, "break");
            break;
          }

          ix += direction;
        }

        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_bins_local);

        Z_TRACE_IF(MSEG_TRACE_IF, "global %d (i = %" slint_fmt ", ix = %" slint_fmt ")", abs(ix - i), i, ix);

        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_bins_global);

        mpi_binning_exec_global(&gb, do_binning_bins, do_binning_prefixes,
#ifdef MSEG_ROOT
          SL_DEFCON(mseg.root),
#else
          -1,
#endif
          size, rank, comm);

        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_bins_global);

      } else ix += direction;

      rti_tstop(rti_tid_mpi_select_exact_generic_while_check_bins);

      if (round == 1 && i == border_lo)
      {
        /* do initialization */
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_round1);

#ifdef MSEG_ROOT
        if (SL_DEFCON(mseg.root) < 0 || SL_DEFCON(mseg.root) == rank)
#endif
        {
          total_counts = 0;
#ifdef elem_weight
          total_weights = 0;
#endif
          for (j = 0; j < gb.bm->nbins; ++j)
          {
            if (docounts) total_counts += *gb_counts(&gb, border_bins[i], j);
#ifdef elem_weight
            if (doweights) total_weights += *gb_weights(&gb, border_bins[i], j);
#endif
          }

          if (docounts) Z_TRACE_IF(MSEG_TRACE_IF, "total_counts = %" slint_fmt, total_counts);
#ifdef elem_weight
          if (doweights) Z_TRACE_IF(MSEG_TRACE_IF, "total_weights = %" slweight_fmt , total_weights);
#endif

#ifdef MSEG_INFO
# ifdef elem_weight
          if (doweights) SL_DEFCON(mseg.info_total_weights) = total_weights;
# endif
#endif

          init_partconds_intern(nparts, pci, pconds, nparts, total_counts, elem_weight_ifelse(total_weights, 0));

          /* init lowest and highest part (sentinels) */
          Z_TRACE_IF(MSEG_TRACE_IF, "init lowest border:");
          border_init(docounts, doweights, &border_infos[border_lo - 1], 0, total_counts, elem_weight_ifelse(total_weights, 0));
          Z_TRACE_IF(MSEG_TRACE_IF, "init highest border:");
          border_init(docounts, doweights, &border_infos[border_hi + 1], 1, total_counts, elem_weight_ifelse(total_weights, 0));

#ifdef MSEG_BORDER_UPDATE_REDUCTION
          /* init+update forwards */
          for (j = border_lo; j <= border_hi; ++j)
          {
            Z_TRACE_IF(MSEG_TRACE_IF, "init border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_init(docounts, doweights, &border_infos[borders[j]], -1, total_counts, elem_weight_ifelse(total_weights, 0));

            Z_TRACE_IF(MSEG_TRACE_IF, "update update %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_update(docounts, doweights, &border_infos[borders[j]], &pci[borders[j]], 1, 0);
          }
#endif

          /* [init+]update backwards */
          for (j = border_hi; j >= border_lo; --j)
          {
#ifndef MSEG_BORDER_UPDATE_REDUCTION
            Z_TRACE_IF(MSEG_TRACE_IF, "init border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_init(docounts, doweights, &border_infos[borders[j]], -1, total_counts, elem_weight_ifelse(total_weights, 0));
#endif
            Z_TRACE_IF(MSEG_TRACE_IF, "update border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_update(docounts, doweights, &border_infos[borders[j]], &pci[borders[j]], -1, 1);
          }
        }
        
        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_round1);
      }

do_partitioning:

      rti_tstart(rti_tid_mpi_select_exact_generic_while_check_pre);

#ifdef MSEG_ROOT
      if (SL_DEFCON(mseg.root) < 0 || SL_DEFCON(mseg.root) == rank)
#endif
      {
        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": check", i, borders[i]);

        if (docounts)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": crange: %" slint_fmt " - %" slint_fmt, i, borders[i], border_infos[borders[i]].crange[0], border_infos[borders[i]].crange[1]);
#ifdef elem_weight
        if (doweights)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": wrange: %" slweight_fmt " - %" slweight_fmt, i, borders[i], border_infos[borders[i]].wrange[0], border_infos[borders[i]].wrange[1]);
#endif

        border_update_update(docounts, doweights, &border_infos[borders[i]], &pci[borders[i]], direction, 1);

        if (docounts)
        {
          current_clo = border_infos[borders[i]].ccurrent[LO] - border_infos[borders[i]].crange[0];
          current_chi = border_infos[borders[i]].ccurrent[HI] - border_infos[borders[i]].crange[0];

          Z_ASSERT_IF(MSEG_ASSERT_IF, current_clo <= current_chi);
          Z_ASSERT_IF(MSEG_ASSERT_IF, 0 <= current_clo);
        
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": currents count: %" slcount_fmt " - %" slcount_fmt " (range: %" slcount_fmt ")",
            i, borders[i], current_clo, current_chi, current_chi - current_clo);

        } else { current_clo = -1; current_chi = 1; }

#ifdef elem_weight
        if (doweights)
        {
          current_wlo = border_infos[borders[i]].wcurrent[LO] - border_infos[borders[i]].wrange[0];
          current_whi = border_infos[borders[i]].wcurrent[HI] - border_infos[borders[i]].wrange[0];

          Z_ASSERT_IF(MSEG_ASSERT_IF, current_wlo <= current_whi);

          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": currents weight: %" slweight_fmt " - %" slweight_fmt " (range: %" slweight_fmt ")",
            i, borders[i], current_wlo, current_whi, current_whi - current_wlo);

        } else { current_wlo = -1; current_whi = 1; }
#endif
      }

      rti_tstop(rti_tid_mpi_select_exact_generic_while_check_pre);

      refine = 0;

      if (!finalize)
      {
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_part);

#ifdef MSEG_ROOT
        if (SL_DEFCON(mseg.root) < 0 || SL_DEFCON(mseg.root) == rank)
#endif
        {
          gcs = 0;
#ifdef elem_weight
          gws = 0;
#endif

          k = 0;
          
          while (1)
          {
            /* check for HIT */

            /* HIT if max count already skipped */
            if (k == 0 && current_chi < 0) break;

            /* if between min/max counts */
            if (current_clo <= 0 && current_chi >= 0)
            {
#ifdef elem_weight
              if (doweights)
              {
                Z_TRACE_IF(MSEG_TRACE_IF, "go to next: %d && %d", (current_chi > 0), (current_wlo > 0));

                /* go to next if max count not reached AND min weight not reached */
                if (current_chi > 0 && current_wlo > 0) goto donthit;
              }
#endif

#ifndef MSEG_DISABLE_BEST_CHOICE
              /* look ahead for a better stop */
              if (k < gb.bm->nbins && (!docounts || current_chi - *gb_counts(&gb, border_bins[i], k) >= 0))
              {
#ifdef elem_weight
                if (doweights)
                {
                  /* continue if weights will improve */
                  if (z_abs(current_wlo + current_whi) > z_abs(current_wlo + current_whi - 2 * *gb_weights(&gb, border_bins[i], k))) goto donthit;

                } else
#endif
                {
                  /* continue if counts will improve */
                  if (z_abs(current_clo + current_chi) > z_abs(current_clo + current_chi - 2 * *gb_counts(&gb, border_bins[i], k))) goto donthit;
                }
              }
#endif

              /* HIT if there is no better stop */
              break;
            }

donthit:

            /* HIT in the worst case */
            if (k >= gb.bm->nbins) break;

            /* skip k-th bin */
            
            if (docounts)
            {
              gc = *gb_counts(&gb, border_bins[i], k);

              current_clo -= gc;
              current_chi -= gc;

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": k = %" slint_fmt ", currents count: %" slcount_fmt " - %" slcount_fmt ", gc = %" slint_fmt ", gcs = %" slint_fmt,
                i, borders[i], k, current_clo, current_chi, gc, gcs);

            } else gc = 0;

#ifdef elem_weight
            if (doweights)
            {
              gw = *gb_weights(&gb, border_bins[i], k);

              current_wlo -= gw;
              current_whi -= gw;

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": k = %" slint_fmt ", currents weight: %" slweight_fmt " - %" slweight_fmt ", gw = %" slweight_fmt ", gws = %" slweight_fmt,
                i, borders[i], k, current_wlo, current_whi, gw, gws);

            } else gw = 0;
#endif

            /* check for REFINE */

            /* stop and refine if max count is skipped OR (min count AND max weight is skipped) */
            if (current_chi < 0
#ifdef elem_weight
              /* '(!docounts || current_clo < 0)' is omitted, because if !docounts then current_clo = -1 */
              /* 'doweights &&' is omitted, because if !doweights then current_whi = 1 */
              || (current_clo < 0 && current_whi < 0)
#endif
              )
            {
              /* stop for REFINE if we do not know the counts, or */
              /* if counts are known and there are more than one elements to refine (otherwise a HIT follows next) */
              if (!docounts || gc > 1)
              {
                refine = 1;
                break;
              }
            }

            gcs += gc;
            gc = 0;

#ifdef elem_weight
            gws += gw;
            gw = 0;
#endif
            
            ++k;
          }

          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": %s k = %" slint_fmt, i, borders[i], (refine)?"REFINE":"HIT", k);

          Z_ASSERT_IF(MSEG_ASSERT_IF, ((refine && k < gb.bm->nbins) || (!refine && k <= gb.bm->nbins)));

          k = (k + 1) * ((refine)?-1:1);
        }

#ifdef MSEG_ROOT
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_part_root);
        if (SL_DEFCON(mseg.root) >= 0) MPI_Bcast(&k, 1, int_mpi_datatype, SL_DEFCON(mseg.root), comm);
        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_part_root);
#endif

        refine = (k < 0);
        if (k < 0) k = -k;
        --k;

        /* refine or hit */
        if (refine)
        {
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": refine bin: %" slint_fmt " @ k = %" slint_fmt, i, borders[i], border_bins[i], k);

          border_bins[i] = mpi_binning_refine(&gb, border_bins[i], k, sp, borders[i] + 1, size, rank, comm);
          
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": new bin: %" slint_fmt, i, borders[i], border_bins[i]);

        } else
        {
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": hit bin: %" slint_fmt " @ k = %" slint_fmt, i, borders[i], border_bins[i], k);

          mpi_binning_hit(&gb, border_bins[i], k, sp, borders[i] + 1, size, rank, comm);
        }

        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_part);

      } else
      {
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_final);

#ifdef MSEG_ROOT
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_final_root);
        if (SL_DEFCON(mseg.root) >= 0)
        {
#if defined(elem_weight) && defined(sl_weight_intequiv)
          MPI_Bcast(current_cw, 4, weight_mpi_datatype, SL_DEFCON(mseg.root), comm);
#else
          if (docounts)
            MPI_Bcast(current_c, 2, int_mpi_datatype, SL_DEFCON(mseg.root), comm);
# ifdef elem_weight
          if (doweights)
            MPI_Bcast(current_w, 2, weight_mpi_datatype, SL_DEFCON(mseg.root), comm);
# endif
#endif
        }
        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_final_root);
#endif

        switch (SL_DEFCON(mseg.finalize_mode))
        {
          case SL_MSEG_FM_ALLORNOTHING:
            Z_TRACE_IF(MSEG_TRACE_IF, "finalize mode: all or nothing");
#ifdef elem_weight
            if (doweights)
            {
              nothing = (current_wlo < ((border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0]) - current_whi));
              Z_TRACE_IF(MSEG_TRACE_IF, "weight: %" slweight_fmt " vs. %" slweight_fmt " -> %s", current_wlo, (border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0]) - current_whi, ((nothing)?"NOTHING":"ALL"));

            } else
#endif
            {
              nothing = (current_clo < ((border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0]) - current_chi));
              Z_TRACE_IF(MSEG_TRACE_IF, "count: %" slint_fmt " vs. %" slint_fmt " -> %s", (slint_t) current_clo, (slint_t) ((border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0]) - current_chi), ((nothing)?"NOTHING":"ALL"));
            }

            if (nothing)
            {
              final_mc = final_dc = 0;
#ifdef elem_weight
              final_mw = final_dw = 0;
#endif
              lc_min = lc_max = 0;

            } else
            {
#ifdef elem_weight
              if (doweights)
              {
                final_mw = final_dw = border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0];

              } else
#endif
              {
                final_mc = final_dc = border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0];
              }
              lc_min = lc_max = border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0];
            }

#ifdef elem_weight
            if (doweights)
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": weights: final_mw = %" slweight_fmt ", final_dw = %" slweight_fmt, i, borders[i], final_mw, final_dw);
            else
#endif
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": counts: final_mc = %" slint_fmt ", final_dc = %" slint_fmt, i, borders[i], final_mc, final_dc);

            if (docounts)
              final_gcs = (nothing)?0:(border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0]);
#ifdef elem_weight
            if (doweights)
              final_gws = (nothing)?0:(border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0]);

            lcw2gcw = 0;
#endif
            break;
          
          case SL_MSEG_FM_MIDDLE:
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": finalize mode: middle (CURRENTLY MISSING!)", i, borders[i]);
            break;
          
          case SL_MSEG_FM_EXACT:
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": finalize mode: exact", i, borders[i]);

#ifdef elem_weight
            if (doweights)
            {
              /* middle of min/max weight */
              final_mw = (current_wlo + current_whi) / 2.0;

              /* min. part of weight to contribute */
              final_dw = z_max(0, final_mw - gb_prefix_weight(&gb, border_bins[i])[0]);

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": weights: prefix_weight: %" slweight_fmt ", final_mw = %" slweight_fmt ", final_dw = %" slweight_fmt,
                i, borders[i], gb_prefix_weight(&gb, border_bins[i])[0], final_mw, final_dw);

            } else
#endif
            {
              /* middle of min/max count */
              final_mc = (current_clo + current_chi) / 2;

              /* min. part of count to contribute */
              final_dc = z_max(0, final_mc - gb_prefix_count(&gb, border_bins[i])[0]);

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": counts: prefix_count: %" slcount_fmt ", final_mc = %" slint_fmt ", final_dc = %" slint_fmt,
                i, borders[i], gb_prefix_count(&gb, border_bins[i])[0], final_mc, final_dc);
            }

            lc_min = current_clo - gb_prefix_count(&gb, border_bins[i])[0];
            lc_max = current_chi - gb_prefix_count(&gb, border_bins[i])[0];

#ifdef elem_weight
            lcw2gcw = 1;
#endif
            break;

          default:
            Z_ERROR("finalize mode %" slint_fmt " not supported!", SL_DEFCON(mseg.finalize_mode));
            break;
        }

        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": lc_min = %" slint_fmt ", lc_max = %" slint_fmt, i, borders[i], lc_min, lc_max);

        mpi_binning_finalize(&gb, border_bins[i], elem_weight_ifelse(0, final_dc), elem_weight_ifelse(final_dw, 0), lc_min, lc_max, &final_lcs, elem_weight_ifelse(&final_lws, NULL), sp, borders[i] + 1, size, rank, comm);

        if (docounts)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": lcs_final = %" slcount_fmt, i, borders[i], final_lcs);
#ifdef elem_weight
        if (doweights)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": lws_final = %" slweight_fmt, i, borders[i], final_lws);
#endif

        gcs = gc = 0;
#ifdef elem_weight
        gws = gw = 0;
#endif

        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": next border: %" slint_fmt " <= %" slint_fmt " + %" slint_fmt " <= %" slint_fmt,
          i, borders[i], border_lo, i, direction, border_hi);

        /* if the next open border is really the _next_ border */
        if (border_lo <= i + direction && i + direction <= border_hi && borders[i + direction] == borders[i] + direction)
        {
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": next border: %" slint_fmt " == %" slint_fmt " + %" slint_fmt,
            i, borders[i], borders[i + direction], borders[i], direction);

#ifdef elem_weight
          if (doweights)
          {
            if (lcw2gcw)
            {
              /* need to determine the exact global counts/weights from the local counts/weights */
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": get gcw from lcw", i, borders[i]);
              
# ifdef MSEG_ROOT
              if (SL_DEFCON(mseg.root) >= 0)
              {
#  ifdef sl_weight_intequiv
                MPI_Reduce(final_lcws, final_gcws, weight_factor, weight_mpi_datatype, MPI_SUM, SL_DEFCON(mseg.root), comm);
#  else
                MPI_Reduce(&final_lcs, &final_gcs, 1, int_mpi_datatype, MPI_SUM, SL_DEFCON(mseg.root), comm);
                MPI_Reduce(&final_lws, &final_gws, 1, weight_mpi_datatype, MPI_SUM, SL_DEFCON(mseg.root), comm);
#  endif
              } else
# endif
              {
# ifdef sl_weight_intequiv
                sl_MPI_Allreduce(final_lcws, final_gcws, weight_factor, weight_mpi_datatype, MPI_SUM, comm, size, rank);
# else
                sl_MPI_Allreduce(&final_lcs, &final_gcs, 1, int_mpi_datatype, MPI_SUM, comm, size, rank);
                sl_MPI_Allreduce(&final_lws, &final_gws, 1, weight_mpi_datatype, MPI_SUM, comm, size, rank);
# endif
              }
            }

          } else
#endif
          {
            /* global counts is just what we selected above */
            final_gcs = final_mc;
          }

#ifdef MSEG_ROOT
          if (SL_DEFCON(mseg.root) < 0 || SL_DEFCON(mseg.root) == rank)
#endif
          {
            gcs = final_gcs;
#ifdef elem_weight
            gws = final_gws;
#endif
          }

          if (docounts)
          {
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": gcs = %" slint_fmt, i, borders[i], gcs);
/*            Z_ASSERT_IF(MSEG_ASSERT_IF, current_clo <= gcs && gcs <= current_chi);*/ /* FIXME: only if exact */
          }
#ifdef elem_weight
          if (doweights)
          {
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": gws = %" slweight_fmt, i, borders[i], gws);
/*            Z_ASSERT_IF(MSEG_ASSERT_IF, current_wlo <= gws && gws <= current_whi);*/  /* FIXME: only if exact */
          }
#endif
        }

        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_final);
      }

      rti_tstart(rti_tid_mpi_select_exact_generic_while_check_post);

#ifdef MSEG_ROOT
      if (SL_DEFCON(mseg.root) < 0 || SL_DEFCON(mseg.root) == rank)
#endif
      {
        border_change_change(docounts, doweights, &border_infos[borders[i]], gcs, gc, elem_weight_ifelse(gws, 0), elem_weight_ifelse(gw, 0));
      }
      
      Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": %s", i, borders[i], (refine)?"REFINE":"REMOVE");

      if (refine)
      {
        borders[i - nborders_removed * direction] = borders[i];
        border_bins[i - nborders_removed * direction] = border_bins[i];

      } else
      {
        ++nborders_removed;

#ifdef MSEG_INFO
        if (SL_DEFCON(mseg.info_finish_rounds)) SL_DEFCON(mseg.info_finish_rounds)[borders[i]] = round;
        SL_DEFCON(mseg.info_finish_rounds_avg) += round;
#endif
      }

      rti_tstop(rti_tid_mpi_select_exact_generic_while_check_post);

      i += direction;
      
      Z_TRACE_IF(MSEG_TRACE_IF, "do partitioning: %" slint_fmt " vs. %" slint_fmt, i, ix);
      
      if (i != ix) goto do_partitioning;
    }
    rti_tstop(rti_tid_mpi_select_exact_generic_while_check);

    mpi_binning_post(&gb, size, rank, comm);

    /* restrict the parts */
    if (direction > 0) border_hi -= nborders_removed;
    else border_lo += nborders_removed;

    Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ": remove %" slint_fmt", lo: %" slint_fmt ", hi: %" slint_fmt "", round, nborders_removed, border_lo, border_hi);

#ifdef MSEG_FORWARD_ONLY
    if (SL_DEFCON(mseg.forward_only))
    {
      /* do not change direction, but if there are min-max bounds, then perform a separate backward pass to update the (remaining) borders */
      if (pci->pcm & (SLPC_COUNTS_MM|SLPC_WEIGHTS_MM))
      {
        for (i = border_hi; i >= border_lo; --i) border_update_update(docounts, doweights, &border_infos[borders[i]], &pci[borders[i]], direction, -1);
      }

    } else
#endif
    {
      /* change direction */
      direction *= -1;
    }
  }

#ifdef MSEG_INFO  
  SL_DEFCON(mseg.info_rounds) = round;
  if (size > 1) SL_DEFCON(mseg.info_finish_rounds_avg) /= (size - 1); else SL_DEFCON(mseg.info_finish_rounds_avg) = 0.0;
#endif

  rti_tstop(rti_tid_mpi_select_exact_generic_while);

  mpi_binning_destroy(&gb, size, rank, comm);

  rti_tstop(rti_tid_mpi_select_exact_generic);

  rti_tstart(rti_tid_mpi_select_exact_generic_sync_exit);
#ifdef SYNC_ON_EXIT
  MPI_Barrier(comm);
#endif
  rti_tstop(rti_tid_mpi_select_exact_generic_sync_exit);

#ifdef VERIFY
  if (size == 1) v = -1;
  else v = mpi_post_check_partconds_intern(s, nelements, nparts, pci, sp->displs, size, rank, comm);
  
  Z_ASSERT_IF(MSEG_ASSERT_IF, v < 0);
  
  Z_NOTICE_IF(rank == 0, "post_check_partconds: %s (%" slint_fmt ")", (v >= 0)?"FAILED":"SUCCESS", v);
#endif

#ifdef PRINT_SDISPLS
  printf("%d: sdispls:", rank);
  for (i = 0; i < nparts; ++i) printf(" %d ", sp->displs[i]);
  printf("\n");
#endif

#ifdef PRINT_STATS
  mpi_select_stats(s, nparts, sp->displs, size, rank, comm);
#endif

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == PRINT_TIMINGS)
  {
    printf("%d: mpi_select_exact_generic: %f\n", rank, rti_tlast(rti_tid_mpi_select_exact_generic));
    printf("%d: mpi_select_exact_generic: sync init: %f\n", rank, rti_tlast(rti_tid_mpi_select_exact_generic_sync_init));
    printf("%d: mpi_select_exact_generic: while: %f\n", rank, rti_tlast(rti_tid_mpi_select_exact_generic_while));
    printf("%d: mpi_select_exact_generic:  check: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check));
    printf("%d: mpi_select_exact_generic:   bins: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_bins));
    printf("%d: mpi_select_exact_generic:    local: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_bins_local));
    printf("%d: mpi_select_exact_generic:    global: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_bins_global));
    printf("%d: mpi_select_exact_generic:   round1: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_round1));
    printf("%d: mpi_select_exact_generic:   pre: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_pre));
    printf("%d: mpi_select_exact_generic:   part: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_part));
    printf("%d: mpi_select_exact_generic:    root: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_part_root));
    printf("%d: mpi_select_exact_generic:   final: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_final));
    printf("%d: mpi_select_exact_generic:    root: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_final_root));
    printf("%d: mpi_select_exact_generic:   post: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_post));
    printf("%d: mpi_select_exact_generic: sync exit: %f\n", rank, rti_tlast(rti_tid_mpi_select_exact_generic_sync_exit));
    printf("%d: mpi_select_exact_generic: rounds: %" slint_fmt "\n", rank, round);
  }
#endif

  return 0;

#undef doweights
}


slint_t mpi_select_exact_generic_grouped(elements_t *s, slint_t nelements, partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_generic_grouped */
{
  slint_t npconds = -1;
  partcond_t *pconds;


  mpi_gather_partconds_grouped(pcond, pcond_comm, group_comm, NULL, &npconds, size, rank, comm);

  pconds = z_alloca(npconds, sizeof(partcond_t));
  
  mpi_gather_partconds_grouped(pcond, pcond_comm, group_comm, pconds, &npconds, size, rank, comm);

  mpi_select_exact_generic_bulk(s, nelements, npconds, pconds, bm, sp, size, rank, comm);
  
  z_freea(pconds);
  
  return 0;
}


#if 0
slint_t mpi_select_exact_generic(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_generic */
{
  const slint_t max_nborders = nparts - 1;
  slint_t border_lo, border_hi, nborders_removed;
/*  slint_t borders[max_nborders], border_bins[max_nborders];

  border_info_t border_infos_[1 + max_nborders + 1], *border_infos = border_infos_ + 1, border_info_old;*/

  slint_t round, direction;

  slint_t i;
/*  slint_t i, j, k, l;*/

  slint_t curr_l, curr_g, curr_p;
  slint_t next_l, next_g, next_p;
  
  slint_t binning_at_once = 1;
  

  border_lo = 0;
  border_hi = max_nborders - 1;

  direction = 1;

  round = 0;
  while (border_lo <= border_hi)
  {
    ++round;

    nborders_removed = 0;

    i = (direction > 0)?border_lo:border_hi;

    next_l = i;
    next_g = -1;
    next_p = -1;
    
    while ((direction > 0)?(i <= border_hi):(i >= border_lo))
    {
      curr_l = next_l;
      curr_g = next_g;
      curr_p = next_p;


      if (border_lo <= curr_g && curr_g <= border_hi)
      {
        /* init global binning at curr_g */
      
        next_p = curr_g;
      }

      if (border_lo <= curr_p && curr_p <= border_hi)
      {
        /* wait global binning at curr_p */
        
      
        /* partitioning at curr_p */
        
        i += binning_at_once * direction;
      }
      
      if (border_lo <= curr_l && curr_l <= border_hi)
      {
        /* local binning at curr_l */

        next_l += binning_at_once * direction;
        next_g = curr_l;
      }
    }

    /* restrict the parts */
    if (direction > 0) border_hi -= nborders_removed;
    else border_lo += nborders_removed;

    /* change direction */
    direction *= -1;
  }
  
  return 0;
}
#endif


#undef SYNC_ON_INIT
#undef SYNC_ON_EXIT
#undef PRINT_SDISPLS
#undef PRINT_STATS
#undef PRINT_TIMINGS
#undef VERIFY
#undef LO
#undef HI



#include "sl_common.h"


slint_t mpi_select_exact_radix(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_radix */
{
  binning_t bm;
  splitter_t sp;

  binning_radix_create(&bm, rhigh, rlow, rwidth, sorted|SL_SORTED_IN);

  sp.displs = sdispls;
  mpi_select_exact_generic_bulk(s, nelements, nparts, pconds, &bm, &sp, size, rank, comm);

  binning_radix_destroy(&bm);

  return 0;
}


slint_t mpi_select_exact_radix_grouped(elements_t *s, slint_t nelements, partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_radix_grouped */
{
  slint_t npconds = -1;
  partcond_t *pconds;


  mpi_gather_partconds_grouped(pcond, pcond_comm, group_comm, NULL, &npconds, size, rank, comm);

  pconds = z_alloca(npconds, sizeof(partcond_t));

  mpi_gather_partconds_grouped(pcond, pcond_comm, group_comm, pconds, &npconds, size, rank, comm);

  mpi_select_exact_radix(s, nelements, npconds, pconds, rhigh, rlow, rwidth, sorted, sdispls, size, rank, comm);

  z_freea(pconds);

  return 0;
}



/* sl_macro MSS_ROOT */
/* sl_macro MSS_TRACE_IF */


#include "sl_common.h"


/*#define sl_pivot_equal(_n_, _i_, _r_)  ((slint_t) (((((double) (_i_) + 1) * (_n_)) / ((_r_) + 1)) + 0.5))*/
#define sl_pivot_equal(_n_, _i_, _r_)  (((_i_) * (_n_)) / ((_r_) + 1))

/*#define MSS_ROOT*/

/* sl_ifdef SL_USE_MPI sl_context CONTEXT_BEGIN mss */
#ifdef MSS_ROOT
const int default_mss_root = -1;  /* sl_ifdef MSS_ROOT sl_endif sl_global sl_context sl_var default_mss_root */
#endif
/* sl_endif sl_context CONTEXT_END mss */

#ifndef MSS_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSS_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSS_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_select_sample_regular(elements_t *s, slint_t nparts, partcond_t *pconds, slint_t nsamples, splitter_t *sp, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_sample_regular */
{
  slint_t i, j;
#ifdef elem_weight
  slweight_t w, wi, wold;
#endif

  const slint_t nslocal = nsamples;
#ifdef MSS_ROOT
  const slint_t nsglobal = nslocal * size;
#endif

  const slint_t nsplitter = nparts - 1;
  slint_t splitter_skip = 0;

#ifdef elem_weight
  slpwkey_t lwskeys[nslocal];
# ifdef MSS_ROOT
  slpwkey_t gwskeys[nsglobal];
# endif
#endif
  slpkey_t lskeys[nslocal];
#ifdef MSS_ROOT
  slpkey_t gskeys[nsglobal];
#endif

  slpkey_t spkeys[nsplitter + 1];
#ifdef elem_weight
  int nlspkeys;
  int rcounts[size], rdispls[size];
  slpkey_t lspkeys[nsplitter];
#endif
  slpkey_t lspkey;

  slint_t lgcounts[2];
#ifdef elem_weight
  slweight_t lgweights[2];
#endif

  partcond_intern_t pci[nparts];

  elements_t gs, e;

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif


  sp->displs[0] = 0;

  if (nparts < 2) return 0;

#ifdef elem_weight
  doweights = ((pconds->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

#ifdef elem_weight
  if (doweights) mpi_elements_get_counts_and_weights(s, 1, lgcounts, lgweights, -1, size, rank, comm);
  else
#endif
    mpi_elements_get_counts(s, &lgcounts[0], &lgcounts[1], -1, size, rank, comm);

  init_partconds_intern(nparts, pci, pconds, nparts, lgcounts[1], elem_weight_ifelse(doweights?lgweights[1]:0, 0));

  Z_TRACE_IF(MSS_TRACE_IF, "counts: %" slint_fmt " / %" slint_fmt, lgcounts[0], lgcounts[1]);
#ifdef elem_weight
  if (doweights)
    Z_TRACE_IF(MSS_TRACE_IF, "weights: %" slweight_fmt " / %" slweight_fmt "", lgweights[0], lgweights[1]);
#endif

  Z_TRACE_IF(MSS_TRACE_IF, "local samples: %" slint_fmt, nslocal);

  /* select local samples */
#ifdef elem_weight
  if (doweights)
  {
    j = 0;
    w = wold = 0;

    for (i = 0; i < nslocal; ++i)
    {
      wi = (i + 1) * lgweights[0] / (nslocal + 1);

      while (w < wi && j < lgcounts[0])
      {
        w += elem_weight(s, j);
        ++j;
      }
    
      if (j < lgcounts[0]) lwskeys[i].pkey = *key_get_pure(elem_key_at(s, j));
      else lwskeys[i].pkey = *key_get_pure(elem_key_at(s, j - 1)) + 1;

      lwskeys[i].weight = w - wold;
    
      wold = w;

      Z_TRACE_IF(MSS_TRACE_IF, "local sample %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, j, lwskeys[i].pkey, lwskeys[i].weight);
    }

  } else
#endif
    for (i = 0; i < nslocal; ++i)
  {
    j = sl_pivot_equal(s->size, i + 1, nslocal);
  
    lskeys[i] = *key_get_pure(elem_key_at(s, j));

    Z_TRACE_IF(MSS_TRACE_IF, "local sample %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt, i, j, lskeys[i]);
  }

#ifdef MSS_ROOT
  /* with root-process (p^2 scaling!!!) */
  if (SL_DEFCON(mss.root) >= 0)
  {
    /* gather local samples at root */
#ifdef elem_weight
    if (doweights) MPI_Gather(lwskeys, nslocal, pwkey_mpi_datatype, gwskeys, nslocal, pwkey_mpi_datatype, SL_DEFCON(mss.root), comm);
    else
#endif
      MPI_Gather(lskeys, nslocal, pkey_mpi_datatype, gskeys, nslocal, pkey_mpi_datatype, SL_DEFCON(mss.root), comm);

    if (rank == SL_DEFCON(mss.root))
    {
      Z_TRACE_IF(MSS_TRACE_IF, "global samples: %" slint_fmt, nsglobal);

      /* prepare global samples */
      elements_alloc(&gs, nsglobal, SLCM_ALL);

      gs.size = nsglobal;

      Z_TRACE_IF(MSS_TRACE_IF, "unsorted global samples");

#ifdef elem_weight
      if (doweights)
      {
        for (i = 0; i < nsglobal; ++i)
        {
          key_set_pure(elem_key_at(&gs, i), gwskeys[i].pkey);
          elem_weight_set(&gs, i, gwskeys[i].weight);

          Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, *key_get_pure(elem_key_at(&gs, i)), elem_weight(&gs, i));
        }

      } else
#endif
        for (i = 0; i < nsglobal; ++i)
      {
        key_set_pure(elem_key_at(&gs, i), gskeys[i]);

        Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt, i, *key_get_pure(elem_key_at(&gs, i)));
      }

      /* sort global samples */
      sort_radix(&gs, NULL, -1, -1, -1);

      Z_TRACE_IF(MSS_TRACE_IF, "sorted global samples");

#ifdef elem_weight
      if (doweights)
      {
        for (i = 0; i < nsglobal; ++i)
        {
          Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, *key_get_pure(elem_key_at(&gs, i)), elem_weight(&gs, i));
        }

      } else
#endif
        for (i = 0; i < nsglobal; ++i)
      {
        Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt, i, *key_get_pure(elem_key_at(&gs, i)));
      }

      Z_TRACE_IF(MSS_TRACE_IF, "splitters: %" slint_fmt, nsplitter);

      /* select splitters from global samples */
#ifdef elem_weight
      if (doweights)
      {
        j = 0;
        w = 0;

        for (i = 0; i < nsplitter; ++i)
        {
          wi = (i + 1) * lgweights[1] / (nsplitter + 1) - (lgweights[1] / (nsplitter + 1) / 2);

          while (w < wi && j < nsglobal)
          {
            w += elem_weight(&gs, j);
            ++j;
          }

          if (j > 0 && (wi - w + elem_weight(&gs, j - 1)) <= (w - wi))
          {
            w -= elem_weight(&gs, j - 1);
            --j;
          }
          spkeys[i] = *key_get_pure(elem_key_at(&gs, j - 1));

          Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt " (target: %" slweight_fmt ")", i, j - 1, spkeys[i], w, wi);
        }
    
      } else
#endif
        for (i = 0; i < nsplitter; ++i)
      {
        j = sl_pivot_equal(gs.size, i + 1, nsplitter);
    
        spkeys[i] = *key_get_pure(elem_key_at(&gs, j));

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt, i, j, spkeys[i]);
      }

      elements_free(&gs);
    }

    /* broadcast splitters */
    MPI_Bcast(&spkeys, nsplitter, pkey_mpi_datatype, SL_DEFCON(mss.root), comm);

  } else
#endif
  {
    /* without root-process */
    elements_t xs;

    elements_alloc(&gs, nslocal, SLCM_ALL);
    elements_alloc(&xs, nslocal / 2 + 1, SLCM_ALL);

    gs.size = nslocal;

    Z_TRACE_IF(MSS_TRACE_IF, "unsorted global samples");

#ifdef elem_weight
    if (doweights)
    {
      for (i = 0; i < gs.size; ++i)
      {
        key_set_pure(elem_key_at(&gs, i), lwskeys[i].pkey);
        elem_weight_set(&gs, i, lwskeys[i].weight);

        Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, *key_get_pure(elem_key_at(&gs, i)), elem_weight(&gs, i));
      }

    } else
#endif
      for (i = 0; i < gs.size; ++i)
    {
      key_set_pure(elem_key_at(&gs, i), lskeys[i]);

      Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt, i, *key_get_pure(elem_key_at(&gs, i)));
    }

    /* parallel sort samples */
    mpi_mergek(&gs, sn_batcher, NULL, merge2_basic_auto_01_x, &xs, size, rank, comm);

    Z_TRACE_IF(MSS_TRACE_IF, "sorted global samples");

#ifdef elem_weight
    if (doweights)
    {
      for (i = 0; i < gs.size; ++i)
      {
        Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, *key_get_pure(elem_key_at(&gs, i)), elem_weight(&gs, i));
      }

    } else
#endif
      for (i = 0; i < gs.size; ++i)
    {
      Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt, i, *key_get_pure(elem_key_at(&gs, i)));
    }

#ifdef elem_weight
    if (doweights)
    {
      nlspkeys = 0;
    
      slweight_t my_sum, my_first, next_first;
      slweight_t prev_slast, next_sfirst;
      MPI_Status status;

      my_sum = 0;
      for (i = 0; i < gs.size; ++i) my_sum += elem_weight(&gs, i);
    
      prev_slast = 0;
      MPI_Exscan(&my_sum, &prev_slast, 1, weight_mpi_datatype, MPI_SUM, comm);
    
      my_first = elem_weight(&gs, 0);
      next_first = 0;
      MPI_Sendrecv(&my_first, 1, weight_mpi_datatype, (rank > 0)?(rank - 1):MPI_PROC_NULL, 0, &next_first, 1, weight_mpi_datatype, (rank + 1 < size)?(rank + 1):MPI_PROC_NULL, 0, comm, &status);

      next_sfirst = prev_slast + my_sum + next_first;

      Z_TRACE_IF(MSS_TRACE_IF, "weights: my_sum: %" slweight_fmt ", my_first: %" slweight_fmt ", next_first: %" slweight_fmt ", prev_slast: %" slweight_fmt ", next_sfirst: %" slweight_fmt, my_sum, my_first, next_first, prev_slast, next_sfirst);

      j = 0;
      w = prev_slast;

      for (i = 0; i < nsplitter; ++i)
      {
        wi = (i + 1) * lgweights[1] / (nsplitter + 1) - (lgweights[1] / (nsplitter + 1) / 2);

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt ": wi: %" slweight_fmt, i, wi);

        if (wi > next_sfirst || wi < prev_slast)
        {
          Z_TRACE_IF(MSS_TRACE_IF, " -> continue");
          continue;
        }

        while (w < wi && j < gs.size)
        {
          w += elem_weight(&gs, j);
          ++j;
        }

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt ": A: j = %" slint_fmt ", w = %" slweight_fmt ", wi: %" slweight_fmt, i, j, w, wi);

        /* step back? */
        if (j > 0 && (w >= wi) && (wi - (w - elem_weight(&gs, j - 1))) <= (w - wi))
        {
          w -= elem_weight(&gs, j - 1);
          --j;
        }

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt ": B: j = %" slint_fmt ", w = %" slweight_fmt ", wi: %" slweight_fmt ", next_sfirst: %" slweight_fmt, i, j, w, wi, next_sfirst);

        if (w < wi)
        {
          if (next_sfirst < wi || (next_sfirst >= wi && (next_sfirst - wi) <= (wi - w)))
          {
            w = next_sfirst;
            ++j;
          }
        }

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt ": C: j = %" slint_fmt ", w = %" slweight_fmt ", wi: %" slweight_fmt, i, j, w, wi);

        if (j > 0 && j - 1 < gs.size)
        {
          lspkeys[nlspkeys] = *key_get_pure(elem_key_at(&gs, j - 1));

          Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt " (target: %" slweight_fmt ")", i, j - 1, lspkeys[nlspkeys], w, wi);

          ++nlspkeys;

        } else Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt " @ %" slint_fmt ": skip", i, j - 1);
      }
      
      Z_TRACE_IF(MSS_TRACE_IF, "local splitters: %d", nlspkeys);

      MPI_Allgather(&nlspkeys, 1, MPI_INT, rcounts, 1, MPI_INT, comm);

      rdispls[0] = 0;
      for (i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + rcounts[i - 1];

      MPI_Allgatherv(lspkeys, nlspkeys, pkey_mpi_datatype, spkeys, rcounts, rdispls, pkey_mpi_datatype, comm);

    } else
#endif
    {
      lspkey = *key_get_pure(elem_key_at(&gs, 0));

      Z_TRACE_IF(MSS_TRACE_IF, "local splitter: %" key_pure_type_fmt, lspkey);
  
      MPI_Allgather(&lspkey, 1, pkey_mpi_datatype, spkeys, 1, pkey_mpi_datatype, comm);
  
      splitter_skip = 1;
    }

    elements_free(&gs);
    elements_free(&xs);
  }

  /* determine splitting positions from splitters with binary search */
  for (i = 0; i < nsplitter; ++i)
  {
    elem_assign_at(s, sp->displs[i], &e);
    e.size = s->size - sp->displs[i];

    sp->displs[i + 1] = sp->displs[i] + sl_search_binary_lt(&e, spkeys[splitter_skip + i]);

    Z_TRACE_IF(MSS_TRACE_IF, "displs %" slint_fmt ": %d (< %" key_pure_type_fmt ")", i + 1, sp->displs[i + 1], spkeys[splitter_skip + i]);

/*    w = 0;
    for (j = 0; j < sp->displs[i + 1] - sp->displs[i]; ++j) w += elem_weight(&e, j);

    Z_TRACE_IF(MSS_TRACE_IF, " -> count: %d, weight: %" slweight_fmt, sp->displs[i + 1] - sp->displs[i], w);*/
  }
  
  return 0;

#undef doweights
}


#undef sl_pivot_equal



/* sl_macro MSM_TRACE_IF */
/* sl_macro MSM_VERIFY */
/* sl_macro MSM_VERIFY_OUTPUT */

#include "sl_common.h"

/*#define MSM_VERIFY
#define MSM_VERIFY_OUTPUT*/

/* sl_ifdef SL_USE_MPI sl_context CONTEXT_BEGIN msm */
const double default_msm_t[4] = { 0 };  /* sl_global sl_context sl_var default_msm_t */

const slint_t default_msm_sync = 0;  /* sl_global sl_context sl_var default_msm_sync */
/* sl_endif sl_context CONTEXT_END msm */


#ifndef MSM_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSM_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSM_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


#ifdef SLDEBUG
# define CHECK_ORDER
#endif


slint_t mpi_sort_merge(elements_t *s0, elements_t *s1, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_merge */
{
  merge2x_f m2x = merge2_memory_adaptive;
  sortnet_f sn = sn_batcher;

#ifdef CHECK_ORDER
  slint_t stages;
  slint_t rorders[2];
#endif

#ifdef key_integer
  sort_radix(s0, xs, -1, -1, -1);
#else
  sort_quick(s0, xs);
#endif

#ifdef CHECK_ORDER
  stages =
#endif
    mpi_mergek(s0, sn, NULL, m2x, xs, size, rank, comm);

#ifdef CHECK_ORDER
  mpi_elements_check_order(s0, 1, rorders, size, rank, comm);

  if (rank == 0) printf("%d: %s (%" slint_fmt ", %" slint_fmt ") with %" slint_fmt " stages\n", rank, (rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1], stages);
#endif

  return 0;
}


slint_t mpi_sort_merge2(elements_t *s0, elements_t *s1, elements_t *xs, slint_t merge_type, slint_t sort_type, double *times, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_merge2 */
{
  double _times[2];

  merge2x_f m2x = merge2_memory_adaptive;
  sortnet_f sn = sn_batcher;

#ifdef CHECK_ORDER
  slint_t stages;
  slint_t rorders[2];
#endif


  if (times == NULL) times = _times;

  MPI_Barrier(comm);

  times[0] = z_time_get_s();
  sort_radix(s0, xs, -1, -1, -1);

  MPI_Barrier(comm);
  times[0] = z_time_get_s() - times[0];

  switch (merge_type)
  {
    case 0:
      m2x = merge2_basic_auto_01_x;
      break;
    case 1:
      m2x = merge2_compo_tridgell;
      break;
    case 2:
      m2x = merge2_compo_hula;
      break;
    default:
      m2x = merge2_memory_adaptive;
      break;
  }

  switch (sort_type)
  {
    case 0:
      sn = sn_batcher;
      break;
    case 1:
      sn = sn_bitonic;
      break;
    case 2:
      sn = sn_odd_even_trans;
      break;
  }

  times[1] = z_time_get_s();

#ifdef CHECK_ORDER
  stages =
#endif
    mpi_mergek(s0, sn, NULL, m2x, xs, size, rank, comm);

  MPI_Barrier(comm);
  times[1] = z_time_get_s() - times[1];

#ifdef CHECK_ORDER
  mpi_elements_check_order(s0, 1, rorders, size, rank, comm);

  if (rank == 0) printf("%d: %s (%" slint_fmt ", %" slint_fmt ") with %" slint_fmt " stages\n", rank, (rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1], stages);
#endif

  return 0;
}


#ifdef key_integer

slint_t mpi_sort_merge_radix(elements_t *s0, elements_t *s1, elements_t *xs, slint_t merge_type, slint_t sort_type, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_merge_radix */
{
  merge2x_f m2x = merge2_memory_adaptive;
  sortnet_f sn = sn_batcher;

#ifdef MSM_VERIFY
  slint_t rorders[2];
#endif


  if (SL_DEFCON(msm.sync)) MPI_Barrier(comm);
  SL_DEFCON(msm.t)[0] = z_time_get_s();
  sort_radix(s0, xs, -1, -1, -1);
  SL_DEFCON(msm.t)[0] = z_time_get_s() - SL_DEFCON(msm.t)[0];

  switch (merge_type)
  {
    case 0:
      m2x = merge2_basic_straight_01_x;
      break;
    case 1:
      m2x = merge2_compo_tridgell;
      break;
    case 2:
      m2x = merge2_compo_hula;
      break;
    default:
      m2x = merge2_memory_adaptive;
      break;
  }

  switch (sort_type)
  {
    case 0:
      sn = sn_batcher;
      break;
    case 1:
      sn = sn_bitonic;
      break;
    case 2:
      sn = sn_odd_even_trans;
      break;
  }

  if (SL_DEFCON(msm.sync)) MPI_Barrier(comm);
  SL_DEFCON(msm.t)[1] = z_time_get_s();

  mpi_mergek(s0, sn, NULL, m2x, xs, size, rank, comm);

  if (SL_DEFCON(msm.sync)) MPI_Barrier(comm);
  SL_DEFCON(msm.t)[1] = z_time_get_s() - SL_DEFCON(msm.t)[1];

#ifdef MSM_VERIFY
  mpi_elements_check_order(s0, 1, rorders, size, rank, comm);

# ifndef MSM_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    Z_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;
}

#endif


#undef CHECK_ORDER



/* sl_macro MSP_TRACE_IF */
/* sl_macro MSP_VERIFY */
/* sl_macro MSP_VERIFY_OUTPUT */


#include "sl_common.h"

/*#define MSP_VERIFY
#define MSP_VERIFY_OUTPUT*/


/* sl_ifdef SL_USE_MPI sl_context CONTEXT_BEGIN msp */
const double default_msp_t[4] = { 0 };  /* sl_global sl_context sl_var default_msp_t */

const slint_t default_msp_sync = 0;  /* sl_global sl_context sl_var default_msp_sync */

const partcond_t *default_msp_r_pc = NULL;  /* sl_global sl_context sl_var default_msp_r_pc */
/* sl_endif sl_context CONTEXT_END msp */


#ifndef MSP_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSP_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSP_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_sort_partition(elements_t *s0, elements_t *s1, elements_t *xs, slint_t part_type, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition */
{
  const slint_t paex_rhigh = -1;
  const slint_t paex_rlow = -1;
  const slint_t paex_rwidth = 3;

  partcond_t pc;
  
  double imba = 0.01;
  
  int scounts[size], sdispls[size], rcounts[size], rdispls[size];

  const slint_t max_nsubs = 4;
  slint_t nsubs;
  MPI_Comm sub_comms[max_nsubs];
  int sub_sizes[max_nsubs], sub_ranks[max_nsubs];


  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[0] = z_time_get_s();
#ifdef key_integer
  sort_radix(s0, xs, -1, -1, -1);
#else
  sort_quick(s0, xs);
#endif
  SL_DEFCON(msp.t)[0] = z_time_get_s() - SL_DEFCON(msp.t)[0];

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[1] = z_time_get_s();

  pc.pcm = SLPC_COUNTS_MM;
  pc.count_min = s0->size * (1.0 - imba);
  pc.count_max = s0->size * (1.0 + imba);

  if (part_type == 1)
  {
#ifdef key_integer
    mpi_partition_exact_radix(s0, &pc, paex_rhigh, paex_rlow, paex_rwidth, SL_SORTED_IN, scounts, rcounts, size, rank, comm);
#else
#endif
  
  } else if (part_type == 2)
  {
    nsubs = 2;

    mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
#ifdef key_integer
    mpi_partition_exact_radix_2groups(s0, &pc, sub_comms[1], NULL, paex_rhigh, paex_rlow, paex_rwidth, scounts, rcounts, size, rank, comm);
#else
#endif
    mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);
  
  }
#ifdef SL_INDEX
    else if (part_type < 0)
  {
    nsubs = (-part_type <= max_nsubs)?-part_type:max_nsubs;

    mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
#ifdef key_integer
    mpi_partition_exact_radix_ngroups(s0, &pc, nsubs, sub_comms, NULL, paex_rhigh, paex_rlow, paex_rwidth, scounts, rcounts, size, rank, comm);
#else
#endif
    mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);
  }
#endif

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[1] = z_time_get_s() - SL_DEFCON(msp.t)[1];


  SL_DEFCON(msp.t)[2] = z_time_get_s();
  
  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

#define xelem_index_not
#define xelem_call \
    MPI_Alltoallv(xelem_buf(s0), scounts, sdispls, xelem_mpi_datatype, xelem_buf(s1), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"
  
  SL_DEFCON(msp.t)[2] = z_time_get_s() - SL_DEFCON(msp.t)[2];

  s1->size = rdispls[size - 1] + rcounts[size - 1];

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[3] = z_time_get_s();

#ifdef key_integer
  sort_radix(s1, xs, -1, -1, -1);
#else
  sort_quick(s1, xs);
#endif

  SL_DEFCON(msp.t)[3] = z_time_get_s() - SL_DEFCON(msp.t)[3];

  return 0;
}


#ifdef key_integer

slint_t mpi_sort_partition_radix(elements_t *s0, elements_t *s1, elements_t *xs, slint_t part_type, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition_radix */
{
  partcond_t pc;
  
  double imba = 0.01;
  
  int scounts[size], sdispls[size], rcounts[size], rdispls[size];

  const slint_t max_nsubs = 4;
  slint_t nsubs;
  MPI_Comm sub_comms[max_nsubs];
  int sub_sizes[max_nsubs], sub_ranks[max_nsubs];

#ifdef MSP_VERIFY
  slint_t rorders[2];
#endif


  /* local sort */
  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[0] = z_time_get_s();
  sort_radix(s0, xs, rhigh, rlow, -1);
  SL_DEFCON(msp.t)[0] = z_time_get_s() - SL_DEFCON(msp.t)[0];

  /* partitioning */
  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[1] = z_time_get_s();

  if (SL_DEFCON(msp.r_pc)) pc = *SL_DEFCON(msp.r_pc);
  else
  {
    pc.pcm = SLPC_COUNTS_MM;
    pc.count_min = -(1.0 - imba);
    pc.count_max = -(1.0 + imba);
  }

#ifdef MSEG_BORDER_UPDATE_REDUCTION
  SL_DEFCON(mseg.border_update_count_reduction) = 0.25;
#endif

  if (part_type == 1)
  {
    mpi_partition_exact_radix(s0, &pc, rhigh, rlow, rwidth, SL_SORTED_IN, scounts, rcounts, size, rank, comm);
  
  } else if (part_type == 2)
  {
    nsubs = 2;

    mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
    mpi_partition_exact_radix_2groups(s0, &pc, sub_comms[1], NULL, rhigh, rlow, rwidth, scounts, rcounts, size, rank, comm);
    mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);
  
  }
#ifdef SL_INDEX
    else if (part_type < 0)
  {
    nsubs = (-part_type <= max_nsubs)?-part_type:max_nsubs;

    mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
    mpi_partition_exact_radix_ngroups(s0, &pc, nsubs, sub_comms, NULL, rhigh, rlow, rwidth, scounts, rcounts, size, rank, comm);
    mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);
  }
#endif

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[1] = z_time_get_s() - SL_DEFCON(msp.t)[1];

  /* all-to-all */
  SL_DEFCON(msp.t)[2] = z_time_get_s();
  
  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

#define xelem_index_not
#define xelem_call \
    MPI_Alltoallv(xelem_buf(s0), scounts, sdispls, xelem_mpi_datatype, xelem_buf(s1), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  /* local merge */
  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[2] = z_time_get_s() - SL_DEFCON(msp.t)[2];

  s0->size = s1->size = rdispls[size - 1] + rcounts[size - 1];

#define HOW_TO_MERGEP  3

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[3] = z_time_get_s();
#if HOW_TO_MERGEP == 0
  sort_radix(s1, s0, rhigh, rlow, -1);
  elem_ncopy(s1, s0, s1->size);
#elif HOW_TO_MERGEP == 1
  mergep_heap_int(s1, s0, size, rdispls, rcounts);
#elif HOW_TO_MERGEP == 2
  mergep_2way_ip_int(s1, s0, size, rdispls, merge2_basic_straight_01_x);
  elem_ncopy(s1, s0, s1->size);
#else
  mergep_2way_ip_int_rec(s1, s0, size, rdispls, merge2_basic_straight_01_x);
  elem_ncopy(s1, s0, s1->size);
#endif
  SL_DEFCON(msp.t)[3] = z_time_get_s() - SL_DEFCON(msp.t)[3];

#undef HOW_TO_MERGEP

#ifdef MSP_VERIFY
  mpi_elements_check_order(s0, 1, rorders, size, rank, comm);

# ifndef MSP_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    SL_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;
}


slint_t mpi_sort_partition_exact_radix(elements_t *s, elements_t *sx, partcond_t *pcond, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition_exact_radix */
{
  elements_t _sx;

  int scounts[size], sdispls[size], rcounts[size], rdispls[size];

#ifdef MSP_VERIFY
  slint_t rorders[2];
#endif


  if (sx == NULL || sx->size < s->size)
  {
    elements_alloc(&_sx, (slint_t) pcond->count_max, SLCM_ALL);
    sx = &_sx;
  }

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[0] = z_time_get_s();
  sort_radix(s, sx, rhigh, rlow, -1);
  SL_DEFCON(msp.t)[0] = z_time_get_s() - SL_DEFCON(msp.t)[0];

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[1] = z_time_get_s();
  mpi_partition_exact_radix(s, pcond, rhigh, rlow, rwidth, SL_SORTED_IN, scounts, rcounts, size, rank, comm);
  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[1] = z_time_get_s() - SL_DEFCON(msp.t)[1];

  SL_DEFCON(msp.t)[2] = z_time_get_s();
  
  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

#define xelem_call \
    MPI_Alltoallv(xelem_buf(s), scounts, sdispls, xelem_mpi_datatype, xelem_buf(sx), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"
  
  SL_DEFCON(msp.t)[2] = z_time_get_s() - SL_DEFCON(msp.t)[2];

  s->size = sx->size = rdispls[size - 1] + rcounts[size - 1];

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[3] = z_time_get_s();
  mergep_heap_int(sx, s, size, rdispls, rcounts);
  SL_DEFCON(msp.t)[3] = z_time_get_s() - SL_DEFCON(msp.t)[3];

  if (sx == &_sx) elements_free(sx);

#ifdef MSP_VERIFY
  mpi_elements_check_order(s, 1, rorders, size, rank, comm);

# ifndef MSP_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    SL_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;
}

#if 0
#ifdef SL_INDEX

slint_t mpi_sort_partition_exact_radix_ngroups(elements_t *s, elements_t *sx, partcond_t *pcond, slint_t ngroups, MPI_Comm *group_comms, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition_exact_radix_ngroups */
{
  slint_t g, i;

  slint_t lgcounts[2], coffset, coffset_add;;
#ifdef elem_weight
  slweight_t lgweights[2];
#endif

  elements_t _sx;
  
  int group_sizes[ngroups];
  int group_ranks[ngroups];
  
  MPI_Comm master_comms[ngroups];
  int master_sizes[ngroups];
  int master_ranks[ngroups];
  
  partcond_p group_pconds[ngroups];
  partcond_t merged_pconds[ngroups];

  int group_sdispls[size], group_scounts[size];
  int current_scounts[size], current_sdispls[size], current_rcounts[size], current_rdispls[size];

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif

#ifdef MSP_VERIFY
  slint_t rorders[2];
#endif

  double _t;


  SL_DEFCON(msp.t)[0] = SL_DEFCON(msp.t)[1] = SL_DEFCON(msp.t)[2] = SL_DEFCON(msp.t)[3] = 0;

  if (sx == NULL || sx->size < s->size)
  {
    elements_alloc(&_sx, pcond->count_max, SLCM_ALL);
    sx = &_sx;
  }

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  _t = z_time_get_s();
  sort_radix(s, sx, rhigh, rlow, -1);
  SL_DEFCON(msp.t)[0] += z_time_get_s() - _t;

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  _t = z_time_get_s();

#ifdef elem_weight
  doweights = ((pcond->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

  /* make absolute pconds */
#ifdef elem_weight
  if (doweights) mpi_elements_get_counts_and_weights(s, 1, lgcounts, lgweights, -1, size, rank, comm);
  else
#endif
    mpi_elements_get_counts(s, &lgcounts[0], &lgcounts[1], -1, size, rank, comm);

  init_partconds(1, pcond, size, lgcounts[1], elem_weight_ifelse(doweights?lgweights[1]:0, 0));

  /* collect and create merged pconds (bottom-up) */
  for (g = ngroups - 1; g >= 0; --g)
  {
    MPI_Comm_size(group_comms[g], &group_sizes[g]);
    MPI_Comm_rank(group_comms[g], &group_ranks[g]);
    
    if (g < ngroups - 1)
    {
      MPI_Comm_split(group_comms[g], (group_ranks[g + 1] == 0)?0:MPI_UNDEFINED, group_ranks[g], &master_comms[g]);
      if (master_comms[g] != MPI_COMM_NULL)
      {
        MPI_Comm_size(master_comms[g], &master_sizes[g]);
        MPI_Comm_rank(master_comms[g], &master_ranks[g]);

      } else master_ranks[g] = -1;

      MPI_Bcast(&master_sizes[g], 1, MPI_INT, 0, group_comms[g + 1]);

    } else
    {
      master_comms[g] = group_comms[g];
      master_sizes[g] = group_sizes[g];
      master_ranks[g] = group_ranks[g];
    }
    
    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": group: %d of %d, master: %d of %d", g, group_ranks[g], group_sizes[g], master_ranks[g], master_sizes[g]);

    group_pconds[g] = z_alloca(master_sizes[g], sizeof(partcond_t));

    if (master_comms[g] != MPI_COMM_NULL) mpi_allgather_partconds((g < ngroups - 1)?&merged_pconds[g + 1]:pcond, group_pconds[g], master_sizes[g], master_ranks[g], master_comms[g]);

    if (g < ngroups - 1) mpi_bcast_partconds(master_sizes[g], group_pconds[g], 0, group_sizes[g + 1], group_ranks[g + 1], group_comms[g + 1]);

    merge_partconds(group_pconds[g], master_sizes[g], &merged_pconds[g]);
    
    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": merged_pconds: %" slint_fmt "", g, merged_pconds[g].pcm);
    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": merged_pconds: count: min/max: %f/%f - low/high: %f/%f", g, merged_pconds[g].count_min, merged_pconds[g].count_max, merged_pconds[g].count_low, merged_pconds[g].count_high);
#ifdef elem_weight
    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": merged_pconds: weight: min/max: %f/%f - low/high: %f/%f", g, merged_pconds[g].weight_min, merged_pconds[g].weight_max, merged_pconds[g].weight_low, merged_pconds[g].weight_high);
#endif
  }

  /* do selects (top-down) */
  for (g = 0; g < ngroups; ++g)
  {
    if (g > 0)
    {
      if (group_pconds[g][0].pcm & SLPC_COUNTS_LH)
      {
        coffset_add = 0;
        MPI_Exscan(&s->size, &coffset_add, 1, int_mpi_datatype, MPI_SUM, group_comms[g - 1]);
        MPI_Bcast(&coffset_add, 1, int_mpi_datatype, 0, group_comms[g]);
    
        Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": count offset = %" slint_fmt " + %" slint_fmt " = %" slint_fmt, g, coffset, coffset_add, coffset + coffset_add);
  
        coffset += coffset_add;

        for (i = 0; i < master_sizes[g]; ++i)
        {
          Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": correct group pcond %" slint_fmt ": count: min/max: %f/%f -> %f/%f",
            g, i, group_pconds[g][i].count_low, group_pconds[g][i].count_high, group_pconds[g][i].count_low - coffset, group_pconds[g][i].count_high - coffset);

          group_pconds[g][i].count_low -= coffset;
          group_pconds[g][i].count_high -= coffset;
        }
      }

#ifdef elem_weight
      if (group_pconds[g][0].pcm & SLPC_WEIGHTS_LH)
      {
        Z_ERROR("correction of weighted low/high group pconds not implemented!");
      }
#endif
    }

    mpi_select_exact_radix(s, 1, master_sizes[g], group_pconds[g], rhigh, rlow, rwidth, SL_SORTED_IN, group_sdispls, group_sizes[g], group_ranks[g], group_comms[g]);

/*    Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": group_sdispls = ", "%d  ", i, master_sizes[g], group_sdispls, g);*/

    if (g < ngroups - 1)
    {
      /* create scounts from sdispls */
      displs2counts(master_sizes[g], group_sdispls, group_scounts, s->size);

/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": group_scounts = ", "%d  ", i, master_sizes[g], group_scounts, g);*/

      mpi_xcounts2ycounts_grouped(group_scounts, master_sizes[g], current_rcounts, group_comms[g + 1], master_comms[g], group_sizes[g], group_ranks[g], group_comms[g]);
    
/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": current_rcounts = ", "%d  ", i, group_sizes[g], current_rcounts, g);*/

/*#define SPARSE_SCOUNTS_FROM_RCOUNTS*/

      /* make scounts from rcounts */
#ifdef SPARSE_SCOUNTS_FROM_RCOUNTS
      mpi_xcounts2ycounts_sparse(current_rcounts, current_scounts, es->size, group_sizes[g], group_ranks[g], group_comms[g]);
#else
      mpi_xcounts2ycounts_all2all(current_rcounts, current_scounts, group_sizes[g], group_ranks[g], group_comms[g]);
#endif

#undef SPARSE_SCOUNTS_FROM_RCOUNTS

/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": current_scounts = ", "%d  ", i, group_sizes[g], current_scounts, g);*/

    } else
    {
      /* create scounts from sdispls */
      displs2counts(master_sizes[g], group_sdispls, current_scounts, s->size);

/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": current_scounts = ", "%d  ", i, group_sizes[g], current_scounts, g);*/

      /* make rcounts from scounts */
      mpi_xcounts2ycounts_all2all(current_scounts, current_rcounts, group_sizes[g], group_ranks[g], group_comms[g]);

/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": current_rcounts = ", "%d  ", i, group_sizes[g], current_rcounts, g);*/
    }

    counts2displs(group_sizes[g], current_scounts, current_sdispls);
    counts2displs(group_sizes[g], current_rcounts, current_rdispls);

    Z_ASSERT(current_sdispls[group_sizes[g] - 1] + current_scounts[group_sizes[g] - 1] <= pcond->count_max);
    Z_ASSERT(current_rdispls[group_sizes[g] - 1] + current_rcounts[group_sizes[g] - 1] <= pcond->count_max);

    SL_DEFCON(msp.t)[1] += z_time_get_s() - _t;

    _t = z_time_get_s();

#define xelem_call \
    MPI_Alltoallv(xelem_buf(s), current_scounts, current_sdispls, xelem_mpi_datatype, xelem_buf(sx), current_rcounts, current_rdispls, xelem_mpi_datatype, group_comms[g]);
#include "sl_xelem_call.h"

    sx->size = current_rdispls[group_sizes[g] - 1] + current_rcounts[group_sizes[g] - 1];

    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": received: %" slint_fmt, g, sx->size);

    SL_DEFCON(msp.t)[2] += z_time_get_s() - _t;

    _t = z_time_get_s();

    /* merge received elements (but not in the last round) */
    mergep_heap_int(sx, s, group_sizes[g], current_rdispls, NULL);

    s->size = sx->size;

    SL_DEFCON(msp.t)[3] += z_time_get_s() - _t;

    _t = z_time_get_s();

    z_freea(group_pconds[g]);
  }

  SL_DEFCON(msp.t)[1] += z_time_get_s() - _t;

  if (sx == &_sx) elements_free(sx);

#ifdef MSP_VERIFY
  mpi_elements_check_order(s, 1, rorders, size, rank, comm);

# ifndef MSP_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    SL_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;

#undef doweights
}

#endif


slint_t mpi_sort_partition_exact_radix_2groups(elements_t *s, elements_t *sx, partcond_t *pcond, MPI_Comm group_comm, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition_exact_radix_2groups */
{
  slint_t i, j;

  slint_t lgcounts[2], coffset;
#ifdef elem_weight
  slweight_t lgweights[2];
#endif

  int group_size, group_rank;
  MPI_Comm master_comm;
  int master_size, master_rank;

  slint_t nparts;
  int group_sdispls[size], group_scounts[size];

  partcond_t *group_pconds, group_pcond;

  int scounts[size], rcounts[size];
  int sdispls[size], rdispls[size];
  slint_t rcounts_total, nsubelements;

  elements_t _sx;

  slint_t sub_elements_sources[size], sub_elements_sizes[size], sub_elements_size;
  elements_t sub_elements[size];

  int *sub_sdispls;

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif

#ifdef MSP_VERIFY
  slint_t rorders[2];
#endif


  Z_TRACE_IF(MSP_TRACE_IF, "s: %" slint_fmt " (max %" slint_fmt "), sx: %" slint_fmt " (max: %" slint_fmt ")", s->size, s->max_size, (sx)?sx->size:0, (sx)?sx->max_size:0);

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[0] = z_time_get_s();
  sort_radix(s, sx, rhigh, rlow, -1);
  SL_DEFCON(msp.t)[0] = z_time_get_s() - SL_DEFCON(msp.t)[0];

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[1] = z_time_get_s();

#ifdef elem_weight
  doweights = ((pcond->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

  /* make absolute part.-conds */
#ifdef elem_weight
  if (doweights) mpi_elements_get_counts_and_weights(s, 1, lgcounts, lgweights, -1, size, rank, comm);
  else
#endif
    mpi_elements_get_counts(s, &lgcounts[0], &lgcounts[1], -1, size, rank, comm);

  init_partconds(1, pcond, size, lgcounts[1], elem_weight_ifelse(doweights?lgweights[1]:0, 0));

  /* init group */
  if (group_comm != MPI_COMM_NULL)
  {
    MPI_Comm_size(group_comm, &group_size);
    MPI_Comm_rank(group_comm, &group_rank);

  } else
  {
    group_size = 1;
    group_rank = 0;
  }
  
  /* init master */
  MPI_Comm_split(comm, (group_rank == 0)?0:MPI_UNDEFINED, rank, &master_comm);
  if (master_comm != MPI_COMM_NULL)
  {
    MPI_Comm_size(master_comm, &master_size);
    MPI_Comm_rank(master_comm, &master_rank);

  } else master_rank = -1;

  /* distribute num. of masters */
  if (group_comm != MPI_COMM_NULL) MPI_Bcast(&master_size, 1, MPI_INT, 0, group_comm);

  Z_TRACE_IF(MSP_TRACE_IF, "%d: group: %d of %d, master: %d of %d", rank, group_rank, group_size, master_rank, master_size);

  nparts = master_size;

  /* create group partcond */
  group_pconds = z_alloca(group_size, sizeof(partcond_t));

  mpi_gather_partconds_grouped(pcond, group_comm, MPI_COMM_NULL, group_pconds, NULL, size, rank, comm);

  merge_partconds(group_pconds, group_size, &group_pcond);

  Z_TRACE_IF(MSP_TRACE_IF, "1st: select grouped");

  /* perform 1st grouped select */  
  mpi_select_exact_radix_grouped(s, 1, &group_pcond, master_comm, group_comm, rhigh, rlow, rwidth, SL_SORTED_IN, group_sdispls, size, rank, comm);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "group_sdispls = ", "%d  ", i, nparts, group_sdispls);*/

  /* create scounts from sdispls */
  displs2counts(nparts, group_sdispls, group_scounts, s->size);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "group_scounts = ", "%d  ", i, nparts, group_scounts);*/

  mpi_xcounts2ycounts_grouped(group_scounts, nparts, rcounts, group_comm, master_comm, size, rank, comm);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "rcounts = ", "%d  ", i, size, rcounts);*/

#define SPARSE_SCOUNTS_FROM_RCOUNTS

  /* make scounts from rcounts */
#ifdef SPARSE_SCOUNTS_FROM_RCOUNTS
  mpi_xcounts2ycounts_sparse(rcounts, scounts, s->size, size, rank, comm);
#else
  mpi_xcounts2ycounts_all2all(rcounts, scounts, size, rank, comm);
#endif

#undef SPARSE_SCOUNTS_FROM_RCOUNTS

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "scounts = ", "%d  ", i, size, scounts);*/

  /* create displs from counts */
  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "sdispls = ", "%d  ", i, size, sdispls);
  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "rdispls = ", "%d  ", i, size, rdispls);*/
  
  /* determine number of sub lists to receive */
  nsubelements = 0;
  for (i = 0; i < size; ++i) nsubelements += (rcounts[i] != 0);

  Z_TRACE_IF(MSP_TRACE_IF, "nsubelements = %" slint_fmt, nsubelements);

  rcounts_total = rdispls[size - 1] + rcounts[size - 1];

  Z_TRACE_IF(MSP_TRACE_IF, "rcounts_total = %" slint_fmt, rcounts_total);

  if (sx == NULL || sx->size < rcounts_total)
  {
    sx = &_sx;

    /* allocate elements for 1st redistribution (need keys and weights) */
    elements_alloc(sx, rcounts_total, SLCM_KEYS|SLCM_WEIGHTS);
  }

  /* 1st redistribution */
#define xelem_call \
  MPI_Alltoallv(xelem_buf(s), scounts, sdispls, xelem_mpi_datatype, xelem_buf(sx), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  /* create sub lists */
  j = 0;
  sub_elements_size = 0;
  for (i = 0; i < size; ++i)
  if (rcounts[i] != 0)
  {
    sub_elements_sources[j] = i;
    sub_elements_sizes[j] = rcounts[i];

    sub_elements_size += sub_elements_sizes[j];

    elem_assign_at(sx, rdispls[i], &sub_elements[j]);
    sub_elements[j].size = rcounts[i];
    
    Z_TRACE_IF(MSP_TRACE_IF, "sub-elements: %" slint_fmt ": source: %" slint_fmt ", size: %" slint_fmt ", displs: %d", j, sub_elements_sources[j], sub_elements_sizes[j], rdispls[i]);

    ++j;
  }

  sub_sdispls = z_alloca(group_size * nsubelements, sizeof(int));

  if (group_pconds[0].pcm & SLPC_COUNTS_LH)
  {
    coffset = 0;
    MPI_Exscan(&sub_elements_size, &coffset, 1, int_mpi_datatype, MPI_SUM, comm);
    MPI_Bcast(&coffset, 1, int_mpi_datatype, 0, group_comm);

    Z_TRACE_IF(MSP_TRACE_IF, "count offset = %" slint_fmt, coffset);

    for (i = 0; i < group_size; ++i)
    {
      Z_TRACE_IF(MSP_TRACE_IF, "correct group pcond %" slint_fmt ": count: min/max: %f/%f -> %f/%f",
        i, group_pconds[i].count_low, group_pconds[i].count_high, group_pconds[i].count_low - coffset, group_pconds[i].count_high - coffset);

      group_pconds[i].count_low -= coffset;
      group_pconds[i].count_high -= coffset;
    }
  }

#ifdef elem_weight
  if (group_pconds[0].pcm & SLPC_WEIGHTS_LH)
  {
    Z_ERROR("correction of weighted low/high group pconds not implemented!");
  }
#endif

  Z_TRACE_IF(MSP_TRACE_IF, "2nd: select with %" slint_fmt " sub-elements", nsubelements);

  /* perform 2nd select */
  mpi_select_exact_radix(sub_elements, nsubelements, group_size, group_pconds, rhigh, rlow, rwidth, SL_SORTED_IN, sub_sdispls, group_size, group_rank, group_comm);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "sub_displs = ", "%d  ", i, nsubelements * group_size, sub_sdispls);*/

  mpi_subxdispls2ycounts(nsubelements, sub_sdispls, sub_elements_sources, sub_elements_sizes, group_comm, group_size, rcounts, size, rank, comm);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "rcounts = ", "%d  ", i, size, rcounts);*/

  z_freea(sub_sdispls);
  z_freea(group_pconds);

  if (master_comm != MPI_COMM_NULL) MPI_Comm_free(&master_comm);

  mpi_xcounts2ycounts_all2all(rcounts, scounts, size, rank, comm);

  SL_DEFCON(msp.t)[1] = z_time_get_s() - SL_DEFCON(msp.t)[1];

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "scounts = ", "%d  ", i, size, scounts);*/

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[2] = z_time_get_s();

  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "rdispls = ", "%d  ", i, size, rdispls);
  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "sdispls = ", "%d  ", i, size, sdispls);*/

#define xelem_call \
  MPI_Alltoallv(xelem_buf(s), scounts, sdispls, xelem_mpi_datatype, xelem_buf(sx), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  SL_DEFCON(msp.t)[2] = z_time_get_s() - SL_DEFCON(msp.t)[2];

  s->size = sx->size = rdispls[size - 1] + rcounts[size - 1];

  if (SL_DEFCON(msp.sync)) MPI_Barrier(comm);
  SL_DEFCON(msp.t)[3] = z_time_get_s();
  mergep_heap_int(sx, s, size, rdispls, rcounts);
  SL_DEFCON(msp.t)[3] = z_time_get_s() - SL_DEFCON(msp.t)[3];

  if (sx == &_sx) elements_free(sx);

#ifdef MSP_VERIFY
  mpi_elements_check_order(s, 1, rorders, size, rank, comm);

# ifndef MSP_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    SL_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;

#undef doweights
}

#endif

#endif



/* sl_macro MSS_TRACE_IF */
/* sl_macro MSS_VERIFY */

#include "sl_common.h"

/*#define MSS_VERIFY*/

/* sl_ifdef SL_USE_MPI sl_context CONTEXT_BEGIN mssp */
const double default_mssp_i_t[3] = { 0 };  /* sl_global sl_context sl_var default_mssp_i_t */
const double default_mssp_p_t[3] = { 0 };  /* sl_global sl_context sl_var default_mssp_p_t */
const double default_mssp_b_t[3] = { 0 };  /* sl_global sl_context sl_var default_mssp_b_t */

const slint_t default_mssp_sync = 0;    /* sl_global sl_context sl_var default_mssp_sync */
const slint_t default_mssp_i_sync = 0;  /* sl_global sl_context sl_var default_mssp_i_sync */
const slint_t default_mssp_p_sync = 0;  /* sl_global sl_context sl_var default_mssp_p_sync */
const slint_t default_mssp_b_sync = 0;  /* sl_global sl_context sl_var default_mssp_b_sync */

const slint_t default_mssp_back_packed = 0;  /* sl_global sl_context sl_var default_mssp_back_packed */
/* sl_endif sl_context CONTEXT_END mssp */


#ifndef MSS_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSS_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSS_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


static inline slint_t key2proc(slpkey_t k, slint_t nmm, slpkey_t *mm)
{
  slint_t l = 0;
  slint_t h = nmm - 1;
  slint_t m = (l + h) / 2;

  while (l < h)
  {
    if (sl_key_pure_cmp_le(k, mm[2 * m + 1])) h = m;
    else l = m + 1;

    m = (l + h) / 2;
  }

  return m;
}


#ifdef key_integer

slint_t mpi_sort_insert_radix(elements_t *s0, elements_t *s1, elements_t *xs, slpkey_t *mmkeys, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_insert_radix */
{
  slint_t i, j;

  slpkey_t all_mmkeys[2 * size];
  
  int scounts[size], sdispls[size], rcounts[size], rdispls[size];


  if (SL_DEFCON(mssp.sync) || SL_DEFCON(mssp.i_sync)) MPI_Barrier(comm);

  /* local rearrange */
  SL_DEFCON(mssp.i_t)[0] = z_time_get_s();

  MPI_Allgather(mmkeys, 2, pkey_mpi_datatype, all_mmkeys, 2, pkey_mpi_datatype, comm);
  
  for (i = 0; i < size; ++i) scounts[i] = 0;

  for (i = 0; i < s0->size; ++i)
  {
    j = key2proc(key_purify(s0->keys[i]), size, all_mmkeys);
    ++scounts[j];
  }

  sdispls[0] = 0;
  for (i = 1; i < size; ++i) sdispls[i] = sdispls[i - 1] + scounts[i - 1];
  
  for (i = 0; i < s0->size; ++i)
  {
    j = key2proc(key_purify(s0->keys[i]), size, all_mmkeys);
    elem_copy_at(s0, i, s1, sdispls[j]);
    ++sdispls[j];
  }

  if (SL_DEFCON(mssp.sync) || SL_DEFCON(mssp.i_sync)) MPI_Barrier(comm);
  SL_DEFCON(mssp.i_t)[0] = z_time_get_s() - SL_DEFCON(mssp.i_t)[0];


  /* all-to-all */
  SL_DEFCON(mssp.i_t)[1] = z_time_get_s();

  sl_MPI_Alltoall_int(scounts, 1, rcounts, 1, comm, size, rank);

  sdispls[0] = rdispls[0] = 0;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  }

  s0->size = s1->size = rdispls[size - 1] + rcounts[size - 1];

  mpi_elements_alltoallv_db(s1, scounts, sdispls, s0, rcounts, rdispls, size, rank, comm);

  if (SL_DEFCON(mssp.sync) || SL_DEFCON(mssp.i_sync)) MPI_Barrier(comm);
  SL_DEFCON(mssp.i_t)[1] = z_time_get_s() - SL_DEFCON(mssp.i_t)[1];


  /* local sort */
  SL_DEFCON(mssp.i_t)[2] = z_time_get_s();

  sort_radix(s0, s1, rhigh, rlow, -1);

  if (SL_DEFCON(mssp.sync) || SL_DEFCON(mssp.i_sync)) MPI_Barrier(comm);
  SL_DEFCON(mssp.i_t)[2] = z_time_get_s() - SL_DEFCON(mssp.i_t)[2];

  return 0;
}


slint_t mpi_sort_presorted_radix(elements_t *s0, elements_t *s1, elements_t *xs, slint_t merge_type, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_presorted_radix */
{
  merge2x_f m2x;


  if (SL_DEFCON(mssp.sync) || SL_DEFCON(mssp.p_sync)) MPI_Barrier(comm);

  /* local sort */
  SL_DEFCON(mssp.p_t)[0] = z_time_get_s();

  sort_radix(s0, xs, rhigh, rlow, -1);

  if (SL_DEFCON(mssp.sync) || SL_DEFCON(mssp.p_sync)) MPI_Barrier(comm);
  SL_DEFCON(mssp.p_t)[0] = z_time_get_s() - SL_DEFCON(mssp.p_t)[0];


  /* merge sorted */
  SL_DEFCON(mssp.p_t)[1] = z_time_get_s();

  switch (merge_type)
  {
    case 0:
      m2x = merge2_basic_straight_01_x;
      break;
    case 1:
      m2x = merge2_compo_tridgell;
      break;
    case 2:
      m2x = merge2_compo_hula;
      break;
    default:
      m2x = merge2_memory_adaptive;
      break;
  }

  mpi_mergek_sorted(s0, m2x, xs, size, rank, comm);

  if (SL_DEFCON(mssp.sync) || SL_DEFCON(mssp.p_sync)) MPI_Barrier(comm);
  SL_DEFCON(mssp.p_t)[1] = z_time_get_s() - SL_DEFCON(mssp.p_t)[1];

  return 0;
}

#endif


typedef struct _tproc_sort_back_data
{
  slint_t nmm;
  slpkey_t *mm;

} tproc_sort_back_data;


static int tproc_sort_back(elements_t *s, slint_t x, void *data) /* sl_func tproc_sort_back */
{
  tproc_sort_back_data *d = data;
  int l, h, m;


  if (key_purify(s->keys[x]) < 0) return MPI_PROC_NULL;

  l = 0;
  h = d->nmm - 1;
  m = (l + h) / 2;

  while (l < h)
  {
    if (key_pure_cmp_le(key_purify(s->keys[x]), d->mm[2 * m + 1])) h = m;
    else l = m + 1;

    m = (l + h) / 2;
  }

  return m;
}


#if 0
static int tproc_mod_sort_back(elements_t *s, slint_t x, void *data, elements_t *mod) /* sl_func tproc_mod_sort_back */
{
  tproc_sort_back_data *d = data;
  int l, h, m;


  if (key_purify(s->keys[x]) < 0) return MPI_PROC_NULL;

  l = 0;
  h = d->nmm - 1;
  m = (l + h) / 2;

  while (l < h)
  {
    if (key_pure_cmp_le(key_purify(s->keys[x]), d->mm[2 * m + 1])) h = m;
    else l = m + 1;

    m = (l + h) / 2;
  }
  
  elem_copy_at(s, x, mod, 0);

  return m;
}


static slint_t tprocs_sort_back(elements_t *s, slint_t x, void *data, int *procs) /* sl_func tprocs_sort_back */
{
  tproc_sort_back_data *d = data;
  int l, h, m;


  if (x < 0) return 1;

  if (key_purify(s->keys[x]) < 0) return 0;

  l = 0;
  h = d->nmm - 1;
  m = (l + h) / 2;

  while (l < h)
  {
    if (key_pure_cmp_le(key_purify(s->keys[x]), d->mm[2 * m + 1])) h = m;
    else l = m + 1;

    m = (l + h) / 2;
  }

  *procs = m;

  return 1;
}


static slint_t tprocs_mod_sort_back(elements_t *s, slint_t x, void *data, int *procs, elements_t *mod) /* sl_func tprocs_mod_sort_back */
{
  tproc_sort_back_data *d = data;
  int l, h, m;


  if (x < 0) return -1;

  if (key_purify(s->keys[x]) < 0) return 0;

  l = 0;
  h = d->nmm - 1;
  m = (l + h) / 2;

  while (l < h)
  {
    if (key_pure_cmp_le(key_purify(s->keys[x]), d->mm[2 * m + 1])) h = m;
    else l = m + 1;

    m = (l + h) / 2;
  }

  *procs = m;

  elem_copy_at(s, x, mod, 0);

  return 1;
}
#endif


/* sl_func sort_back_tproc_count_db sort_back_tproc_count_ip sort_back_tproc_rearrange_db sort_back_tproc_rearrange_ip */
/* sl_var _sort_back sort_back */
TPROC_EXDEF_DEFINE_TPROC(sort_back, tproc_sort_back, static)


slint_t mpi_sort_back(elements_t *sin, elements_t *sout, elements_t *sx, slpkey_t *lh, slint_t ntotal, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_back */
{
  slint_t i, j, original_meas_packed;

  slpkey_t lhs[2 * size];

  tproc_t tproc;
  tproc_sort_back_data data;
  
  elements_t *perm_src, *perm_dst, *perm_tmp, _sx, _bx;


  if (SL_DEFCON(mssp.sync) || SL_DEFCON(mssp.b_sync)) MPI_Barrier(comm);
  SL_DEFCON(mssp.b_t)[0] = z_time_get_s();

  MPI_Allgather(lh, 2, pkey_mpi_datatype, lhs, 2, pkey_mpi_datatype, comm);

  SL_DEFCON(mssp.b_t)[0] = z_time_get_s() - SL_DEFCON(mssp.b_t)[0];

/*  elements_print_keys(sin);*/

  if (sx != NULL && elem_get_block(sx))
  {
    elements_alloc_from_block(&_bx, elem_get_block(sx), elem_get_block_size(sx), 8, -1, SLCM_ALL);

    Z_TRACE_IF(MSS_TRACE_IF, "block mem: size: %" slint_fmt ", max byte per component: %" slint_fmt ", #elements: %" slint_fmt, elem_get_block_size(sx), elem_get_max_byte(), _bx.size);
    
    if (_bx.size > 0) sx = &_bx;
  }

  /* all-to-all specific */
  SL_DEFCON(mssp.b_t)[1] = z_time_get_s();

  data.nmm = size;
  data.mm = lhs;

  tproc_create_tproc(&tproc, tproc_sort_back, TPROC_RESET_NULL, sort_back);

  original_meas_packed = SL_DEFCON(meas.packed); SL_DEFCON(meas.packed) = SL_DEFCON(mssp.back_packed);

  mpi_elements_alltoall_specific(sin, sout, sx, tproc, &data, size, rank, comm);

  SL_DEFCON(meas.packed) = original_meas_packed;

  tproc_free(&tproc);

  SL_DEFCON(mssp.b_t)[1] = z_time_get_s() - SL_DEFCON(mssp.b_t)[1];

/*  elements_print_keys(sout);*/

  /* local permute back */
  SL_DEFCON(mssp.b_t)[2] = z_time_get_s();

  if (sout != NULL) perm_src = sout;
  else perm_src = sin;

  if (sx != NULL && sx->max_size >= perm_src->size) perm_dst = sx;
  else perm_dst = perm_src;

  if (perm_dst == perm_src && sx == NULL)
  {
    elements_alloc(&_sx, 1, SLCM_ALL);
    perm_tmp = &_sx;

  } else perm_tmp = sx;

/*  printf("%d: perm: %p, %p, %p\n", rank, perm_src, perm_dst, perm_tmp);*/

  if (perm_src != perm_dst)
  {
    for (i = 0; i < perm_src->size; ++i)
    {
      j = key_purify(perm_src->keys[i]) - lh[0];
      elem_copy_at(perm_src, i, perm_dst, j);
    }
    perm_dst->size = perm_src->size;

  } else
  {
    for (i = 0; i < perm_src->size; ++i)
    {
      j = key_purify(perm_src->keys[i]) - lh[0];
      while (j != i)
      {
        elem_xchange_at(perm_src, i, perm_src, j, perm_tmp);
        j = key_purify(perm_src->keys[i]) - lh[0];
      }
    }
  }

  if (sout != NULL)
  {
    if (perm_dst != sout) elem_ncopy(perm_dst, sout, perm_dst->size);

  } else
  {
    if (perm_dst != sin) elem_ncopy(perm_dst, sin, perm_dst->size);
  }

  if (perm_tmp == &_sx) elements_free(&_sx);

  if (SL_DEFCON(mssp.sync) || SL_DEFCON(mssp.b_sync)) MPI_Barrier(comm);
  SL_DEFCON(mssp.b_t)[2] = z_time_get_s() - SL_DEFCON(mssp.b_t)[2];

  return 0;
}



/* sl_macro MX2Y_TRACE_IF */

#include "sl_common.h"


#ifndef MX2Y_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MX2Y_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MX2Y_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_xcounts2ycounts_all2all(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_xcounts2ycounts_all2all */
{
  MPI_Alltoall(xcounts, 1, MPI_INT, ycounts, 1, MPI_INT, comm);
  
  return 0;
}


slint_t mpi_xcounts2ycounts_sparse(int *xcounts, int *ycounts, slint_t ytotal, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_xcounts2ycounts_sparse */
{
#define MPI_XCOUNTS2YCOUNTS_SPARSE_TAG  1

  slint_t i;

  MPI_Request req;
  MPI_Status status;

  int ycount;


  for (i = 0; i < size; ++i) ycounts[i] = 0;

  for (i = 0; i < size; ++i)
  if (xcounts[i] != 0)
  {
    MPI_Isend(&xcounts[i], 1, MPI_INT, i, MPI_XCOUNTS2YCOUNTS_SPARSE_TAG, comm, &req);
    MPI_Request_free(&req);
  }
  
  while (ytotal > 0)
  {
    MPI_Recv(&ycount, 1, MPI_INT, MPI_ANY_SOURCE, MPI_XCOUNTS2YCOUNTS_SPARSE_TAG, comm, &status);
    
    ytotal -= ycount;
    ycounts[status.MPI_SOURCE] = ycount;
  }
  
  return 0;

#undef MPI_XCOUNTS2YCOUNTS_SPARSE_TAG
}


slint_t mpi_xcounts2ycounts_grouped(int *xcounts, slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_xcounts2ycounts_grouped */
{
  slint_t i;
  int group_size, group_rank;

  int *all_group_xcounts = NULL;
  MPI_Datatype rdt, rdtr;

  int group_sizes[nxcounts];

  int all_xcounts[size];

  int mscounts[nxcounts], msdispls[nxcounts];
  int mrcounts[nxcounts], mrdispls[nxcounts];

  slint_t dist_total, dist_from, dist_current;

  int my_dist_plan_size = 0;
  slint_t my_dist_plan[size * 2];
  slint_t dist_plan_size, dist_plan[size * 2];
  
  MPI_Request req;
  MPI_Status status;

  
  if (group_comm != MPI_COMM_NULL)
  {
    MPI_Comm_size(group_comm, &group_size);
    MPI_Comm_rank(group_comm, &group_rank);

  } else
  {
    group_size = 1;
    group_rank = 0;
  }

  if (group_rank == 0)
  {
    /* init space and receive type to gather xcounts (transposed) at group-root */
    all_group_xcounts = z_alloca(nxcounts * group_size, sizeof(int));

    MPI_Type_vector(nxcounts, 1, group_size, MPI_INT, &rdt);
    MPI_Type_create_resized(rdt, 0, sizeof(int), &rdtr);
    MPI_Type_commit(&rdtr);

  } else rdtr = MPI_DATATYPE_NULL; /* workaround for an OpenMPI bug (recvtype on non-root process is used in gather routine "intra_binomial") */

  /* gather xcounts (transposed) at group-root */
  MPI_Gather(xcounts, nxcounts, MPI_INT, all_group_xcounts, 1, rdtr, 0, group_comm);

  if (group_rank == 0)
  {
    MPI_Type_free(&rdt);
    MPI_Type_free(&rdtr);

/*    printf("%d: all_group_scounts = ", rank);
    for (i = 0; i < nparts * group_size; ++i) printf("%d  ", all_group_scounts[i]);
    printf("\n");*/

    /* distribute group sizes */
    MPI_Allgather(&group_size, 1, MPI_INT, group_sizes, 1, MPI_INT, master_comm);

    /* re-distribute xcounts among the masters */
    for (i = 0; i < nxcounts; ++i)
    {
      mscounts[i] = group_size;
      msdispls[i] = (i == 0)?0:(msdispls[i - 1] + mscounts[i - 1]);
      mrcounts[i] = group_sizes[i];
      mrdispls[i] = (i == 0)?0:(mrdispls[i - 1] + mrcounts[i - 1]);
    }

    MPI_Alltoallv(all_group_xcounts, mscounts, msdispls, MPI_INT, all_xcounts, mrcounts, mrdispls, MPI_INT, master_comm);

    z_freea(all_group_xcounts);

/*    printf("%d: all_scounts = ", rank);
    for (i = 0; i < size; ++i) printf("%d  ", all_scounts[i]);
    printf("\n");*/
  }

  /* distribute xcounts to the group members */
  if (group_rank == 0)
  {
    /* determine total number of elements to distribute in this group */
    dist_total = 0;
    for (i = 0; i < size; ++i) dist_total += all_xcounts[i];

    /* determine elements for every group member */
    dist_from = 0;
    for (i = 0; i < group_size; ++i)
    {
      /* distribute equally */
      dist_current = (slint_t) ceil((double) dist_total / (group_size - i));

      dist_total -= dist_current;

/*      printf("%d: %" sl_int_type_fmt " receives %" sl_int_type_fmt "\n", rank, i, dist_current);*/

      /* determine how the current elements are distributed (make a distribution plan) */
      dist_plan_size = 0;
      while (dist_current > 0 && dist_from < size)
      {
        if (all_xcounts[dist_from] > 0)
        {
          dist_plan[dist_plan_size + 0] = z_min(dist_current, all_xcounts[dist_from]);
          dist_plan[dist_plan_size + 1] = dist_from;

          all_xcounts[dist_from] -= dist_plan[dist_plan_size + 0];
          dist_current -= dist_plan[dist_plan_size + 0];
          
          dist_plan_size += 2;
          
        } else
        {
          ++dist_from;
          continue;
        }
      }

/*      printf("%d: plan for %" sl_int_type_fmt ": ", rank, i);
      for (j = 0; j < dist_plan_size; j += 2) printf("(%" sl_int_type_fmt ",%" sl_int_type_fmt ")  ", dist_plan[j], dist_plan[j + 1]);
      printf("\n");*/

      /* send (or copy) a distribution plan */
      if (i == group_rank)
      {
        memcpy(my_dist_plan, dist_plan, sizeof(slint_t) * dist_plan_size);
        my_dist_plan_size = dist_plan_size;

      } else
      {
        MPI_Isend(dist_plan, dist_plan_size, int_mpi_datatype, i, 0, group_comm, &req);
        MPI_Request_free(&req);
      }

/*      printf("%d: my plan: ", rank);
      for (j = 0; j < my_dist_plan_size; j += 2) printf("(%" sl_int_type_fmt ",%" sl_int_type_fmt ")  ", my_dist_plan[j], my_dist_plan[j + 1]);
      printf("\n");*/
    }

  } else
  {
    /* receive the distribution plan */
    MPI_Recv(my_dist_plan, size * 2, int_mpi_datatype, 0, 0, group_comm, &status);
    MPI_Get_count(&status, int_mpi_datatype, &my_dist_plan_size);
  }

  /* make ycounts from the distribution plan */
  for (i = 0; i < size; ++i) ycounts[i] = 0;

  for (i = 0; i < my_dist_plan_size; i += 2) ycounts[my_dist_plan[i + 1]] = my_dist_plan[i + 0];

  return 0;
}


slint_t mpi_subxdispls2ycounts(slint_t nsubs, int *sub_xdispls, slint_t *sub_sources, slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_subxdispls2ycounts */
{
  slint_t i, j, k;
  slint_t all_nsubs[sub_size];

  int your_sub_plans[sub_size * nsubs * 2], my_sub_plan[2 * 2 * size];
  int sub_plan_scounts[sub_size], sub_plan_sdispls[sub_size], sub_plan_rcounts[sub_size], sub_plan_rdispls[sub_size];


  Z_TRACE_IF(MX2Y_TRACE_IF, "nsubs = %" slint_fmt "", nsubs);
  Z_TRACE_ARRAY_IF(MX2Y_TRACE_IF, i, nsubs * sub_size, "%d  ", sub_xdispls[i], "sub_xdispls = ");
  Z_TRACE_ARRAY_IF(MX2Y_TRACE_IF, i, nsubs, "%" slint_fmt "  ", sub_sources[i], "sub_sources = ");
  Z_TRACE_ARRAY_IF(MX2Y_TRACE_IF, i, nsubs, "%" slint_fmt "  ", sub_sizes[i], "sub_sizes = ");

  MPI_Allgather(&nsubs, 1, int_mpi_datatype, all_nsubs, 1, int_mpi_datatype, sub_comm);
  
  for (i = 0; i < sub_size; ++i)
  {
    for (j = 0; j < nsubs; ++j)
    {
      your_sub_plans[(i * nsubs + j) * 2 + 0] = sub_sources[j];
      your_sub_plans[(i * nsubs + j) * 2 + 1] = ((i + 1 < sub_size)?sub_xdispls[(i + 1) * nsubs + j]:sub_sizes[j]) - sub_xdispls[i * nsubs + j];
    }
    sub_plan_scounts[i] = nsubs * 2;
    sub_plan_sdispls[i] = i * nsubs * 2;
    
    sub_plan_rcounts[i] = all_nsubs[i] * 2;
    sub_plan_rdispls[i] = (i == 0)?0:(sub_plan_rdispls[i - 1] + sub_plan_rcounts[i - 1]);
  }
  
  MPI_Alltoallv(your_sub_plans, sub_plan_scounts, sub_plan_sdispls, MPI_INT, my_sub_plan, sub_plan_rcounts, sub_plan_rdispls, MPI_INT, sub_comm);

  for (i = 0; i < size; ++i) ycounts[i] = 0;
  
  k = 0;
  for (i = 0; i < sub_size; i++)
  for (j = 0; j < all_nsubs[i]; ++j)
  {
    ycounts[my_sub_plan[k + 0]] += my_sub_plan[k + 1];
    k += 2;
  }
  
  return 0;
}
