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
#include <mpi.h>

#include "config_fmm_sort.h"
#include "rename_fmm_sort.h"

/*#define DO_VERBOSE
#define DO_VALIDATE
#define DO_CHECKSUM
#define DO_TIMING*/
#define DO_TIMING_SYNC

#include "common.h"

#include "mpi_fmm_sort.h"


/*#define FMM_SORT_RADIX_1BIT*/

/*#define SORT_BACK_ONE2NCHARGES*/

/*#define SL_OUTPUT_TO_FILE*/

#define DO_MPI_INIT

#define MPI_SENDRECV_REPLACE_MAX_SIZE  10000000

#define ALLTOALLV_PACKED(_p_, _n_)  ((_p_) >= 1024 && (_n_) <= 200)

#define MIN_AUXMEM_ALLOC  1*1024*1024

#ifdef ALLTOALLV_PACKED
# define ALLTOALLV_PACKED_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define ALLTOALLV_PACKED_CMD(_cmd_)  Z_NOP()
#endif

#ifdef FMM_SORT_RADIX_1BIT
# define fmm_sort_radix(_prefix_, _s_, _sx_, _h_, _l_, _w_)  _prefix_##sort_radix_1bit(_s_, _sx_, _h_, _l_)
#else
# define fmm_sort_radix(_prefix_, _s_, _sx_, _h_, _l_, _w_)  _prefix_##sort_radix(_s_, _sx_, _h_, _l_, _w_)
#endif

#ifdef WITH_FCOMM
# define FCOMM_IFELSE(_if_, _else_) _if_
#else
# define FCOMM_IFELSE(_if_, _else_) _else_
#endif


int mpi_fmm_sort_front_part = 0;
int mpi_fmm_sort_back_part = 0;

int mpi_fmm_sort_front_merge_presorted = 0;


static INTEGER_C mpi_log2_floor(INTEGER_C v)
{
  INTEGER_C x = 0;

  v >>= 1;

  while (v)
  {
    x++;
    v >>= 1;
  }

  return x;
}


#define AUXMEM_NBLOCKS  3

static void auxmem_init(void *mem0, void *mem1, pint_t *mem_sizes, void **auxmem_blocks, slint_t *auxmem_sizes, slint_t *auxmem_max)
{
  if (mem_sizes && mem_sizes[0] > 0 && mem0) { auxmem_blocks[0] = mem0; auxmem_sizes[0] = mem_sizes[0]; }
  else { auxmem_blocks[0] = NULL; auxmem_sizes[0] = 0; }

  if (mem_sizes && mem_sizes[1] > 0 && mem1) { auxmem_blocks[1] = mem1; auxmem_sizes[1] = mem_sizes[1]; }
  else { auxmem_blocks[1] = NULL; auxmem_sizes[1] = 0; }

  *auxmem_max = (auxmem_sizes[0] < auxmem_sizes[1])?1:0;

#ifdef MIN_AUXMEM_ALLOC
  if (auxmem_sizes[0] < MIN_AUXMEM_ALLOC && auxmem_sizes[1] < MIN_AUXMEM_ALLOC)
  {
    auxmem_blocks[2] = malloc(MIN_AUXMEM_ALLOC);
    auxmem_sizes[2] = MIN_AUXMEM_ALLOC;
    *auxmem_max = 2;
  } else
#endif
  { auxmem_blocks[2] = NULL; auxmem_sizes[2] = 0; }
}


static void auxmem_release(void **auxmem_blocks, slint_t *auxmem_sizes)
{
  if (auxmem_blocks[2]) free(auxmem_blocks[2]);
}


#ifndef NO_SL_FRONT

static void mpi_fmm_sort_front_merge_body(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type,
 FINT8_TYPE_C *fcomm)
{
  front_slint_t front_type, i, highest;

#ifndef NOT_sl_front_xqsa0
  front_xqsa0_elements_t s0, *sx0, smem0;
  front_xqsa0_merge2x_f m20 = front_xqsa0_merge2_compo_hula;  /* front_xqsa0_merge2_basic_straight_01_x front_xqsa0_merge2_compo_tridgell front_xqsa0_merge2_compo_hula */
#endif
#ifndef NOT_sl_front_xqsaI
  front_xqsaI_elements_t s1, *sx1, smem1;
  front_xqsaI_merge2x_f m21 = front_xqsaI_merge2_compo_hula;  /* front_xqsaI_merge2_basic_straight_01_x front_xqsaI_merge2_compo_tridgell front_xqsaI_merge2_compo_hula */
#endif
#ifndef NOT_sl_front_xqsaX
  front_xqsaX_elements_t s2, *sx2, smem2;
  front_xqsaX_merge2x_f m22 = front_xqsaX_merge2_compo_hula;  /* front_xqsaX_merge2_basic_straight_01_x front_xqsaX_merge2_compo_tridgell front_xqsaX_merge2_compo_hula */
#endif
#ifndef NOT_sl_front_xq_a0
  front_xq_a0_elements_t s3, *sx3, smem3;
  front_xq_a0_merge2x_f m23 = front_xq_a0_merge2_compo_hula;  /* front_xq_a0_merge2_basic_straight_01_x front_xq_a0_merge2_compo_tridgell front_xq_a0_merge2_compo_hula */
#endif
#ifndef NOT_sl_front_xq_aI
  front_xq_aI_elements_t s4, *sx4, smem4;
  front_xq_aI_merge2x_f m24 = front_xq_aI_merge2_compo_hula;  /* front_xq_aI_merge2_basic_straight_01_x front_xq_aI_merge2_compo_tridgell front_xq_aI_merge2_compo_hula */
#endif
#ifndef NOT_sl_front_xq_aX
  front_xq_aX_elements_t s5, *sx5, smem5;
  front_xq_aX_merge2x_f m25 = front_xq_aX_merge2_compo_hula;  /* front_xq_aX_merge2_basic_straight_01_x front_xq_aX_merge2_compo_tridgell front_xq_aX_merge2_compo_hula */
#endif

  void *auxmem_blocks[AUXMEM_NBLOCKS];
  slint_t auxmem_sizes[AUXMEM_NBLOCKS], auxmem_max;

#ifdef DO_TIMING
  double t[3];
  front_slint_t presorted_merge_rounds = 0;
#endif

#ifdef VALIDATE
  front_slint_t o, l;
# ifdef CHECKSUM
  unsigned int crc32_in, crc32_out;
# endif
#endif

  MPI_Comm comm;
  int size, rank, world_rank;
#ifdef DO_MPI_INIT
  int flag;
#endif

#ifdef SL_OUTPUT_TO_FILE
  char output_file_str[32];
  FILE *sl_debug_fstream;
#endif


#ifdef DO_MPI_INIT
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(NULL, NULL);
#endif

  comm = FCOMM_IFELSE(((fcomm)?MPI_Comm_f2c(*fcomm):MPI_COMM_NULL), MPI_COMM_NULL);

  if (comm == MPI_COMM_NULL) comm = MPI_COMM_WORLD;

  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

#define I_AM_MASTER  ((comm == MPI_COMM_SELF)?(world_rank == 0):(rank == 0))

  if (type && *type >= 0)
  {
    switch (*type)
    {
      case 0: front_type = 1; break;
      case 1: front_type = 0; break;
    }

  } else
  {
    if (scr != NULL) front_type = 0;
    else front_type = 1;
  }
  
  front_type *= 3;

  if (addr_desc == NULL)
  {
    INFO_CMD(
      if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: addr_desc = %p\n", addr_desc);
    );
    front_type += 1;

  } else
  {
    INFO_CMD(
      if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: addr_desc = %" PARAM_INTEGER_FMT ", %" PARAM_INTEGER_FMT "\n", addr_desc[0], addr_desc[1]);
    );
    if (addr_desc[0] == 0) front_type += 0;
    else if (addr_desc[0] == sizeof(INTEGER_C) && addr_desc[1] <= 0) front_type += 1;
    else front_type += 2;
  }

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: type = %" front_slint_fmt "\n", front_type);
  );

  if (I_AM_MASTER)
  {
    switch (front_type)
    {
#ifdef NOT_sl_front_xqsa0
      case 0:
        fprintf(stderr, "error: sl_front_xqsa0 required!!!\n");
        break;
#endif
#ifdef NOT_sl_front_xqsaI
      case 1:
        fprintf(stderr, "error: sl_front_xqsaI required!!!\n");
        break;
#endif
#ifdef NOT_sl_front_xqsaX
      case 2:
        fprintf(stderr, "error: sl_front_xqsaX required!!!\n");
        break;
#endif
#ifdef NOT_sl_front_xq_a0
      case 3:
        fprintf(stderr, "error: sl_front_xq_a0 required!!!\n");
        break;
#endif
#ifdef NOT_sl_front_xq_aI
      case 4:
        fprintf(stderr, "error: sl_front_xq_aI required!!!\n");
        break;
#endif
#ifdef NOT_sl_front_xq_aX
      case 5:
        fprintf(stderr, "error: sl_front_xq_aX required!!!\n");
        break;
#endif
    }
  }

#ifdef SL_OUTPUT_TO_FILE
  sprintf(output_file_str, "sort_front.debug.%d", rank);
  front_xqsa0(front_xqsa0_sl_debug_fstream = )
  front_xqsaI(front_xqsaI_sl_debug_fstream = )
  front_xqsaX(front_xqsaX_sl_debug_fstream = )
  front_xq_a0(front_xq_a0_sl_debug_fstream = )
  front_xq_aI(front_xq_aI_sl_debug_fstream = )
  front_xq_aX(front_xq_aX_sl_debug_fstream = ) sl_debug_fstream = fopen(output_file_str, "w");
#endif

  front_xqsa0(SL_FMM_VAR(front_xqsa0_SL_DEFCON(mpi.rank)) = )
  front_xqsaI(SL_FMM_VAR(front_xqsaI_SL_DEFCON(mpi.rank)) = )
  front_xqsaX(SL_FMM_VAR(front_xqsaX_SL_DEFCON(mpi.rank)) = )
  front_xq_a0(SL_FMM_VAR(front_xq_a0_SL_DEFCON(mpi.rank)) = )
  front_xq_aI(SL_FMM_VAR(front_xq_aI_SL_DEFCON(mpi.rank)) = )
  front_xq_aX(SL_FMM_VAR(front_xq_aX_SL_DEFCON(mpi.rank)) = ) rank;

  SL_FMM_CONFIG_VAR(fmm_front_aX) = addr_desc[0];

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: mem0: %p, mem1: %p, mem_sizes: %p (%" PARAM_INTEGER_FMT ", %" PARAM_INTEGER_FMT ")\n", mem0, mem1, mem_sizes, (mem_sizes?mem_sizes[0]:-1), (mem_sizes?mem_sizes[1]:-1));
  );

  auxmem_init(mem0, mem1, mem_sizes, auxmem_blocks, auxmem_sizes, &auxmem_max);

  if (auxmem_sizes[auxmem_max] > 0)
  {
    switch (front_type)
    {
#ifndef NOT_sl_front_xqsa0
      case 0:
        front_xqsa0_elements_alloc_from_blocks(&smem0, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx0 = &smem0;
        break;
#endif
#ifndef NOT_sl_front_xqsaI
      case 1:
        front_xqsaI_elements_alloc_from_blocks(&smem1, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx1 = &smem1;
        break;
#endif
#ifndef NOT_sl_front_xqsaX
      case 2:
        front_xqsaX_elements_alloc_from_blocks(&smem2, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx2 = &smem2;
        break;
#endif
#ifndef NOT_sl_front_xq_a0
      case 3:
        front_xq_a0_elements_alloc_from_blocks(&smem3, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx3 = &smem3;
        break;
#endif
#ifndef NOT_sl_front_xq_aI
      case 4:
        front_xq_aI_elements_alloc_from_blocks(&smem4, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx4 = &smem4;
        break;
#endif
#ifndef NOT_sl_front_xq_aX
      case 5:
        front_xq_aX_elements_alloc_from_blocks(&smem5, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx5 = &smem5;
        break;
#endif
    }

  } else
  {
    front_xqsa0(sx0 = NULL;)
    front_xqsaI(sx1 = NULL;)
    front_xqsaX(sx2 = NULL;)
    front_xq_a0(sx3 = NULL;)
    front_xq_aI(sx4 = NULL;)
    front_xq_aX(sx5 = NULL;)
  }

  INFO_CMD(
    if (I_AM_MASTER)
    {
      printf(INFO_PRINT_PREFIX "front: sizeof(front_slint_t) = %d byte - #elements: %" PARAM_INTEGER_FMT " - mem0: %p / %" PARAM_INTEGER_FMT " bytes, mem1: %p / %" PARAM_INTEGER_FMT " bytes\n", (int) sizeof(front_(slint_t)), *n, mem0, (mem_sizes)?mem_sizes[0]:0, mem1, (mem_sizes)?mem_sizes[1]:0);
      printf(INFO_PRINT_PREFIX "front: ibox:      sizeof(front_slkey_t)   = %d byte\n", (int) sizeof(front_(slkey_t)));
      printf(INFO_PRINT_PREFIX "front:  xyz: %d * sizeof(front_sldata0_t) = %d * %d byte\n", (int) front_(sl_data0_size_c), (int) front_(sl_data0_size_c), (int) sizeof(front_(sldata0_t)));
      printf(INFO_PRINT_PREFIX "front:    q: %d * sizeof(front_sldata1_t) = %d * %d byte\n", (int) front_(sl_data1_size_c), (int) front_(sl_data1_size_c), (int) sizeof(front_(sldata1_t)));
      printf(INFO_PRINT_PREFIX "front:  scr: %d * sizeof(front_sldata2_t) = %d * %d byte\n", (int) front_(sl_data2_size_c), (int) front_(sl_data2_size_c), (int) sizeof(front_(sldata2_t)));
/*      printf(INFO_PRINT_PREFIX "front: addr:\n");
      switch (front_type)
      {
        case 0:
        case 3:
          printf(INFO_PRINT_PREFIX "front: addr: %d * sizeof(front_sldata3_t) = %d * %d byte\n", (int) front_xqsa0_sl_data3_size_c, (int) front_xqsa0_sl_data3_size_c, (int) sizeof(front_xqsa0_sldata3_t));
          break;
        case 1:
        case 4:
          printf(INFO_PRINT_PREFIX "front: addr: %d * sizeof(front_sldata3_t) = %d * %d byte\n", (int) front_xqsaI_sl_data3_size_c, (int) front_xqsaI_sl_data3_size_c, (int) sizeof(front_xqsaI_sldata3_t));
          break;
        case 2:
        case 5:
          printf(INFO_PRINT_PREFIX "front: addr: %d * sizeof(front_sldata3_t) = %d * %d byte\n", (int) front_xqsaX_sl_data3_size_c, (int) front_xqsaX_sl_data3_size_c, (int) sizeof(front_xqsaX_sldata3_t));
         break;
      }
      printf(INFO_PRINT_PREFIX "front: presorted: %d\n", mpi_fmm_sort_front_merge_presorted);*/
    }
  );

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0:
      front_xqsa0_elem_set_size(&s0, *n);
      front_xqsa0_elem_set_max_size(&s0, *n);
      front_xqsa0_elem_set_keys(&s0, ibox);
      front_xqsa0_elem_set_data(&s0, xyz, q, scr);
      break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1:
      front_xqsaI_elem_set_size(&s1, *n);
      front_xqsaI_elem_set_max_size(&s1, *n);
      front_xqsaI_elem_set_keys(&s1, ibox);
      front_xqsaI_elem_set_data(&s1, xyz, q, scr, addr);
      break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2:
      front_xqsaX_elem_set_size(&s2, *n);
      front_xqsaX_elem_set_max_size(&s2, *n);
      front_xqsaX_elem_set_keys(&s2, ibox);
      front_xqsaX_elem_set_data(&s2, xyz, q, scr, addr);
      break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3:
      front_xq_a0_elem_set_size(&s3, *n);
      front_xq_a0_elem_set_max_size(&s3, *n);
      front_xq_a0_elem_set_keys(&s3, ibox);
      front_xq_a0_elem_set_data(&s3, xyz, q);
      break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4:
      front_xq_aI_elem_set_size(&s4, *n);
      front_xq_aI_elem_set_max_size(&s4, *n);
      front_xq_aI_elem_set_keys(&s4, ibox);
      front_xq_aI_elem_set_data(&s4, xyz, q, addr);
      break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5:
      front_xq_aX_elem_set_size(&s5, *n);
      front_xq_aX_elem_set_max_size(&s5, *n);
      front_xq_aX_elem_set_keys(&s5, ibox);
      front_xq_aX_elem_set_data(&s5, xyz, q, addr);
      break;
#endif
  }

#if defined(VALIDATE) && defined(CHECKSUM)
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0: crc32_in = front_xqsa0_mpi_elements_crc32(&s0, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1: crc32_in = front_xqsaI_mpi_elements_crc32(&s1, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2: crc32_in = front_xqsaX_mpi_elements_crc32(&s2, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3: crc32_in = front_xq_a0_mpi_elements_crc32(&s3, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4: crc32_in = front_xq_aI_mpi_elements_crc32(&s4, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5: crc32_in = front_xq_aX_mpi_elements_crc32(&s5, 1, 1, 1, size, rank, comm); break;
#endif
  }
#endif

  if (*subx != 0) for (i = 0; i < *n; i++) ibox[i] -= *subx;

  if (depth != NULL)
  {
    if (*depth > 0) highest = 3 * *depth - 1; else highest = -1;
  }

  SL_FMM_CONFIG_VAR(fmm_front_key_mask) = ~((front_(slkey_t)) 0);
  if (addr_desc != NULL && addr_desc[1] > 0) SL_FMM_CONFIG_VAR(fmm_front_key_mask) = ~(SL_FMM_CONFIG_VAR(fmm_front_key_mask) << (sizeof(front_(slkey_t)) * 8 - addr_desc[1]));

  INFO_CMD(
    if (I_AM_MASTER)
    {
      if (depth == NULL) printf(INFO_PRINT_PREFIX "front: starting local radix-sort (lowest 3 bit)\n");
      else printf(INFO_PRINT_PREFIX "front: starting local radix-sort (depth = %" PARAM_INTEGER_FMT " -> sorting bits [0..%" front_slint_fmt "])\n", *depth, highest);
    }
  );

  TIMING_SYNC(comm); TIMING_START(t[1]);

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0:
      if (depth == NULL) front_xqsa0_sort_radix_iter(&s0, sx0, 1, 2, 0, -1);
      else fmm_sort_radix(front_xqsa0_, &s0, sx0, highest, -1, -1);
      break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1:
      if (depth == NULL) front_xqsaI_sort_radix_iter(&s1, sx1, 1, 2, 0, -1);
      else fmm_sort_radix(front_xqsaI_, &s1, sx1, highest, -1, -1);
      break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2:
      if (depth == NULL) front_xqsaX_sort_radix_iter(&s2, sx2, 1, 2, 0, -1);
      else fmm_sort_radix(front_xqsaX_, &s2, sx2, highest, -1, -1);
      break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3:
      if (depth == NULL) front_xq_a0_sort_radix_iter(&s3, sx3, 1, 2, 0, -1);
      else fmm_sort_radix(front_xq_a0_, &s3, sx3, highest, -1, -1);
      break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4:
      if (depth == NULL) front_xq_aI_sort_radix_iter(&s4, sx4, 1, 2, 0, -1);
      else fmm_sort_radix(front_xq_aI_, &s4, sx4, highest, -1, -1);
      break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5:
      if (depth == NULL) front_xq_aX_sort_radix_iter(&s5, sx5, 1, 2, 0, -1);
      else fmm_sort_radix(front_xq_aX_, &s5, sx5, highest, -1, -1);
      break;
#endif
  }

  TIMING_SYNC(comm); TIMING_STOP(t[1]);

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: local radix-sort done\n");
  );

#ifdef VALIDATE
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0: l = front_xqsa0_elements_validate_order(&s0, 1); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1: l = front_xqsaI_elements_validate_order(&s1, 1); break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2: l = front_xqsaX_elements_validate_order(&s2, 1); break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3: l = front_xq_a0_elements_validate_order(&s3, 1); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4: l = front_xq_aI_elements_validate_order(&s4, 1); break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5: l = front_xq_aX_elements_validate_order(&s5, 1); break;
#endif
    default: l = 1;
  }
#endif

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0: front_xqsa0_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1: front_xqsaI_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2: front_xqsaX_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3: front_xq_a0_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4: front_xq_aI_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5: front_xq_aX_mpi_datatypes_init(); break;
#endif
  }

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: starting parallel sort (num. proc. = %d)\n", size);
  );

#ifdef MPI_SENDRECV_REPLACE_MAX_SIZE
  front_xqsa0(SL_FMM_VAR(front_xqsa0_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = )
  front_xqsaI(SL_FMM_VAR(front_xqsaI_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = )
  front_xqsaX(SL_FMM_VAR(front_xqsaX_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = )
  front_xq_a0(SL_FMM_VAR(front_xq_a0_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = )
  front_xq_aI(SL_FMM_VAR(front_xq_aI_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = )
  front_xq_aX(SL_FMM_VAR(front_xq_aX_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = ) MPI_SENDRECV_REPLACE_MAX_SIZE;
#endif

  TIMING_SYNC(comm); TIMING_START(t[2]);

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(front_xqsa0_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(front_xqsa0_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m20 = front_xqsa0_merge2_memory_adaptive;
      }
      if (mpi_fmm_sort_front_merge_presorted)
      {
        TIMING_DECL(presorted_merge_rounds =) front_xqsa0_mpi_mergek_sorted2(&s0, front_xqsa0_sn_batcher, NULL, m20, sx0, size, rank, comm);
/*        TIMING_DECL(presorted_merge_rounds =) front_xqsa0_mpi_mergek_sorted2(&s0, front_xqsa0_sn_batcher, NULL, front_xqsa0_merge2_basic_straight_01_x, NULL, size, rank, comm);*/
      }
      else front_xqsa0_mpi_mergek(&s0, front_xqsa0_sn_batcher, NULL, m20, sx0, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(front_xqsaI_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(front_xqsaI_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m21 = front_xqsaI_merge2_memory_adaptive;
      }
      if (mpi_fmm_sort_front_merge_presorted)
      {
        TIMING_DECL(presorted_merge_rounds =) front_xqsaI_mpi_mergek_sorted2(&s1, front_xqsaI_sn_batcher, NULL, m21, sx1, size, rank, comm);
/*        TIMING_DECL(presorted_merge_rounds =) front_xqsaI_mpi_mergek_sorted2(&s1, front_xqsaI_sn_batcher, NULL, front_xqsaI_merge2_basic_straight_01_x, NULL, size, rank, comm);*/
      }
      else front_xqsaI_mpi_mergek(&s1, front_xqsaI_sn_batcher, NULL, m21, sx1, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(front_xqsaX_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(front_xqsaX_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m22 = front_xqsaX_merge2_memory_adaptive;
      }
      if (mpi_fmm_sort_front_merge_presorted)
      {
        TIMING_DECL(presorted_merge_rounds =) front_xqsaX_mpi_mergek_sorted2(&s2, front_xqsaX_sn_batcher, NULL, m22, sx2, size, rank, comm);
/*        TIMING_DECL(presorted_merge_rounds =) front_xqsaX_mpi_mergek_sorted2(&s2, front_xqsaX_sn_batcher, NULL, front_xqsaX_merge2_basic_straight_01_x, NULL, size, rank, comm);*/
      }
      else front_xqsaX_mpi_mergek(&s2, front_xqsaX_sn_batcher, NULL, m22, sx2, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(front_xq_a0_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(front_xq_a0_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m23 = front_xq_a0_merge2_memory_adaptive;
      }
      if (mpi_fmm_sort_front_merge_presorted)
      {
        TIMING_DECL(presorted_merge_rounds =) front_xq_a0_mpi_mergek_sorted2(&s3, front_xq_a0_sn_batcher, NULL, m23, sx3, size, rank, comm);
/*        TIMING_DECL(presorted_merge_rounds =) front_xq_a0_mpi_mergek_sorted2(&s3, front_xq_a0_sn_batcher, NULL, front_xq_a0_merge2_basic_straight_01_x, NULL, size, rank, comm);*/
      }
      else front_xq_a0_mpi_mergek(&s3, front_xq_a0_sn_batcher, NULL, m23, sx3, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(front_xq_aI_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(front_xq_aI_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m24 = front_xq_aI_merge2_memory_adaptive;
      }
      if (mpi_fmm_sort_front_merge_presorted)
      {
        TIMING_DECL(presorted_merge_rounds =) front_xq_aI_mpi_mergek_sorted2(&s4, front_xq_aI_sn_batcher, NULL, m24, sx4, size, rank, comm);
/*        TIMING_DECL(presorted_merge_rounds =) front_xq_aI_mpi_mergek_sorted2(&s4, front_xq_aI_sn_batcher, NULL, front_xq_aI_merge2_basic_straight_01_x, NULL, size, rank, comm);*/
      }
      else front_xq_aI_mpi_mergek(&s4, front_xq_aI_sn_batcher, NULL, m24, sx4, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(front_xq_aX_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(front_xq_aX_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m25 = front_xq_aX_merge2_memory_adaptive;
      }
      if (mpi_fmm_sort_front_merge_presorted)
      {
        TIMING_DECL(presorted_merge_rounds =) front_xq_aX_mpi_mergek_sorted2(&s5, front_xq_aX_sn_batcher, NULL, m25, sx5, size, rank, comm);
/*        TIMING_DECL(presorted_merge_rounds =) front_xq_aX_mpi_mergek_sorted2(&s5, front_xq_aX_sn_batcher, NULL, front_xq_aX_merge2_basic_straight_01_x, NULL, size, rank, comm);*/
      }
      else front_xq_aX_mpi_mergek(&s5, front_xq_aX_sn_batcher, NULL, m25, sx5, size, rank, comm);
      break;
#endif
  }

  TIMING_SYNC(comm); TIMING_STOP(t[2]);

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: parallel sort done\n");
  );

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0: front_xqsa0_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1: front_xqsaI_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2: front_xqsaX_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3: front_xq_a0_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4: front_xq_aI_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5: front_xq_aX_mpi_datatypes_release(); break;
#endif
  }

  if (*subx != 0) for (i = 0; i < *n; i++) ibox[i] += *subx;

  auxmem_release(auxmem_blocks, auxmem_sizes);

#ifdef VALIDATE
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0: o = front_xqsa0_mpi_elements_validate_order(&s0, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1: o = front_xqsaI_mpi_elements_validate_order(&s1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2: o = front_xqsaX_mpi_elements_validate_order(&s2, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3: o = front_xq_a0_mpi_elements_validate_order(&s3, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4: o = front_xq_aI_mpi_elements_validate_order(&s4, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5: o = front_xq_aX_mpi_elements_validate_order(&s5, 1, size, rank, comm); break;
#endif
    default: o = 1;
  }
  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: sorting_order: %s - local_order: %s\n", (!o)?"success":"FAILED", (!l)?"success":"failed");
  );

# ifdef CHECKSUM
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0: crc32_out = front_xqsa0_mpi_elements_crc32(&s0, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1: crc32_out = front_xqsaI_mpi_elements_crc32(&s1, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2: crc32_out = front_xqsaX_mpi_elements_crc32(&s2, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3: crc32_out = front_xq_a0_mpi_elements_crc32(&s3, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4: crc32_out = front_xq_aI_mpi_elements_crc32(&s4, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5: crc32_out = front_xq_aX_mpi_elements_crc32(&s5, 1, 1, 1, size, rank, comm); break;
#endif
  }
  INFO_CMD(
    if (I_AM_MASTER) printf("front: checksum: %s (%X vs. %X)\n", (crc32_in == crc32_out)?"success":"FAILED", crc32_in, crc32_out);
  );
# endif
#endif

#ifdef SL_OUTPUT_TO_FILE
  fclose(sl_debug_fstream);
  front_xqsa0(front_xqsa0_sl_debug_fstream = )
  front_xqsaI(front_xqsaI_sl_debug_fstream = )
  front_xqsaX(front_xqsaX_sl_debug_fstream = )
  front_xq_a0(front_xq_a0_sl_debug_fstream = )
  front_xq_aI(front_xq_aI_sl_debug_fstream = )
  front_xq_aX(front_xq_aX_sl_debug_fstream = ) sl_debug_fstream = NULL;
#endif

  TIMING_SYNC(comm); TIMING_STOP(t[0]);

  TIMING_CMD(
    if (I_AM_MASTER) printf(TIMING_PRINT_PREFIX "mpi_fmm_sort_front_merge_body: %f  %f  %f  %" front_slint_fmt "\n", t[0], t[1], t[2], presorted_merge_rounds);
  );
#undef I_AM_MASTER
}


#ifdef WITH_SORT_FRONT_LOAD

static void mpi_fmm_sort_front_part_body(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type,
 front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax,
 pint_t rebalance_only,
 FINT8_TYPE_C *fcomm)
{
  front_slint_t front_type, i, highest, nlocal;

#ifndef NOT_sl_front_xqsaIl
  front_xqsaIl_elements_t s0, *sx0, smem0, *bx0, bmem0;
  front_xqsaIl_partcond_t pc0;
#endif
#ifndef NOT_sl_front_xq_aIl
  front_xq_aIl_elements_t s1, *sx1, smem1, *bx1, bmem1;
  front_xq_aIl_partcond_t pc1;
#endif
#ifndef NOT_sl_front_xqsaI
  front_xqsaI_elements_t s2, *sx2, smem2, *bx2, bmem2;
  front_xqsaI_partcond_t pc2;
#endif
#ifndef NOT_sl_front_xq_aI
  front_xq_aI_elements_t s3, *sx3, smem3, *bx3, bmem3;
  front_xq_aI_partcond_t pc3;
#endif

  void *auxmem_blocks[AUXMEM_NBLOCKS];
  slint_t auxmem_sizes[AUXMEM_NBLOCKS], auxmem_max;

  front_xqsaIl_slweight_t local_weights, total_weights;

#ifdef ALLTOALLV_PACKED
  int local_packed, global_packed, original_packed;
#endif

#ifdef DO_TIMING
  double t[6];
#endif

#ifdef VALIDATE
  front_slint_t o, l;
# ifdef CHECKSUM
  unsigned int crc32_in, crc32_out;
# endif
#endif

  MPI_Comm comm;
  int size, rank, world_rank;
#ifdef DO_MPI_INIT
  int flag;
#endif

  REAL_C _imba = 0;
  pint_t _nmin = *n;
  pint_t _nmax = *n;

  int *scounts, *sdispls, *rcounts, *rdispls;

  
#ifdef DO_MPI_INIT
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(NULL, NULL);
#endif

  comm = FCOMM_IFELSE(((fcomm)?MPI_Comm_f2c(*fcomm):MPI_COMM_NULL), MPI_COMM_NULL);

  if (comm == MPI_COMM_NULL) comm = MPI_COMM_WORLD;

  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

#define I_AM_MASTER  ((comm == MPI_COMM_SELF)?(world_rank == 0):(rank == 0))

  if (type && *type >= 2)
  {
    switch (*type - 2)
    {
      case 0: front_type = 1; break;
      case 1: front_type = 0; break;
    }

  } else
  {
    if (scr != NULL) front_type = 0;
    else front_type = 1;
  }

  if (load == NULL) front_type += 2;

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: type = %" front_slint_fmt "\n", front_type);
  );

  if (I_AM_MASTER)
  {
    switch (front_type)
    {
#ifdef NOT_sl_front_xqsaIl
      case 0:
        fprintf(stderr, "error: sl_front_xqsaIl required!!!\n");
        break;
#endif
#ifdef NOT_sl_front_xq_aIl
      case 1:
        fprintf(stderr, "error: sl_front_xq_aIl required!!!\n");
        break;
#endif
    }
  }
  
  if (imba == NULL) imba = &_imba;
  if (nmin == NULL) nmin = &_nmin;
  if (nmax == NULL) nmax = &_nmax;

  front_xqsaI(SL_FMM_VAR(front_xqsaIl_SL_DEFCON(mpi.rank)) = )
  front_xq_aI(SL_FMM_VAR(front_xq_aIl_SL_DEFCON(mpi.rank)) = )
  front_xqsaI(SL_FMM_VAR(front_xqsaI_SL_DEFCON(mpi.rank)) = )
  front_xq_aI(SL_FMM_VAR(front_xq_aI_SL_DEFCON(mpi.rank)) = ) rank;

  auxmem_init(mem0, mem1, mem_sizes, auxmem_blocks, auxmem_sizes, &auxmem_max);

  if (auxmem_sizes[auxmem_max] > 0)
  {
    switch (front_type)
    {
#ifndef NOT_sl_front_xqsaIl
      case 0:
        front_xqsaIl_elements_alloc_from_blocks(&smem0, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx0 = &smem0;
        front_xqsaIl_elem_set_block(&bmem0, auxmem_blocks[auxmem_max]);
        front_xqsaIl_elem_set_block_size(&bmem0, auxmem_sizes[auxmem_max]);
        bx0 = &bmem0;
        break;
#endif
#ifndef NOT_sl_front_xq_aIl
      case 1:
        front_xq_aIl_elements_alloc_from_blocks(&smem1, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx1 = &smem1;
        front_xq_aIl_elem_set_block(&bmem1, auxmem_blocks[auxmem_max]);
        front_xq_aIl_elem_set_block_size(&bmem1, auxmem_sizes[auxmem_max]);
        bx1 = &bmem1;
        break;
#endif
#ifndef NOT_sl_front_xqsaI
      case 2:
        front_xqsaI_elements_alloc_from_blocks(&smem2, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx2 = &smem2;
        front_xqsaI_elem_set_block(&bmem2, auxmem_blocks[auxmem_max]);
        front_xqsaI_elem_set_block_size(&bmem2, auxmem_sizes[auxmem_max]);
        bx2 = &bmem2;
        break;
#endif
#ifndef NOT_sl_front_xq_aI
      case 3:
        front_xq_aI_elements_alloc_from_blocks(&smem3, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx3 = &smem3;
        front_xq_aI_elem_set_block(&bmem3, auxmem_blocks[auxmem_max]);
        front_xq_aI_elem_set_block_size(&bmem3, auxmem_sizes[auxmem_max]);
        bx3 = &bmem3;
        break;
#endif
    }

  } else
  {
    front_xqsaIl(sx0 = bx0 = NULL;)
    front_xq_aIl(sx1 = bx1 = NULL;)
    front_xqsaI(sx2 = bx2 = NULL;)
    front_xq_aI(sx3 = bx3 = NULL;)
  }


  INFO_CMD(
    if (I_AM_MASTER)
    {
      printf(INFO_PRINT_PREFIX "front: sizeof(front_slint_t) = %d byte - #elements: %" PARAM_INTEGER_FMT " - mem0: %p / %" PARAM_INTEGER_FMT " bytes, mem1: %p / %" PARAM_INTEGER_FMT " bytes\n", (int) sizeof(front_(slint_t)), *n, mem0, (mem_sizes)?mem_sizes[0]:0, mem1, (mem_sizes)?mem_sizes[1]:0);
      printf(INFO_PRINT_PREFIX "front: ibox:      sizeof(front_slkey_t)   = %d byte\n", (int) sizeof(front_(slkey_t)));
      printf(INFO_PRINT_PREFIX "front:  xyz: %d * sizeof(front_sldata0_t) = %d * %d byte\n", (int) front_(sl_data0_size_c), (int) front_(sl_data0_size_c), (int) sizeof(front_(sldata0_t)));
      printf(INFO_PRINT_PREFIX "front:    q: %d * sizeof(front_sldata1_t) = %d * %d byte\n", (int) front_(sl_data1_size_c), (int) front_(sl_data1_size_c), (int) sizeof(front_(sldata1_t)));
      printf(INFO_PRINT_PREFIX "front:  scr: %d * sizeof(front_sldata2_t) = %d * %d byte\n", (int) front_(sl_data2_size_c), (int) front_(sl_data2_size_c), (int) sizeof(front_(sldata2_t)));
      printf(INFO_PRINT_PREFIX "front: addr:\n");
      printf(INFO_PRINT_PREFIX "front: load: %d * sizeof(front_sldata4_t) = %d * %d byte\n", (int) front_xqsaIl_sl_data4_size_c, (int) front_xqsaIl_sl_data4_size_c, (int) sizeof(front_xqsaIl_sldata4_t));
      printf(INFO_PRINT_PREFIX "front: imba: %f, nmin: %" PARAM_INTEGER_FMT ", nmax: %" PARAM_INTEGER_FMT "\n", *imba, *nmin, *nmax);
      printf(INFO_PRINT_PREFIX "front: rebalance_only: %" PARAM_INTEGER_FMT "\n", rebalance_only);
    }
  );

  nlocal = *n;

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0:
      front_xqsaIl_elem_set_size(&s0, nlocal);
      front_xqsaIl_elem_set_max_size(&s0, *nmax);
      front_xqsaIl_elem_set_keys(&s0, ibox);
      front_xqsaIl_elem_set_data(&s0, xyz, q, scr, addr, load);
      break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1:
      front_xq_aIl_elem_set_size(&s1, nlocal);
      front_xq_aIl_elem_set_max_size(&s1, *nmax);
      front_xq_aIl_elem_set_keys(&s1, ibox);
      front_xq_aIl_elem_set_data(&s1, xyz, q, addr, load);
      break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2:
      front_xqsaI_elem_set_size(&s2, nlocal);
      front_xqsaI_elem_set_max_size(&s2, *nmax);
      front_xqsaI_elem_set_keys(&s2, ibox);
      front_xqsaI_elem_set_data(&s2, xyz, q, scr, addr);
      break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3:
      front_xq_aI_elem_set_size(&s3, nlocal);
      front_xq_aI_elem_set_max_size(&s3, *nmax);
      front_xq_aI_elem_set_keys(&s3, ibox);
      front_xq_aI_elem_set_data(&s3, xyz, q, addr);
      break;
#endif
  }

#if defined(VALIDATE) && defined(CHECKSUM)
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0: crc32_in = front_xqsaIl_mpi_elements_crc32(&s0, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1: crc32_in = front_xq_aIl_mpi_elements_crc32(&s1, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2: crc32_in = front_xqsaI_mpi_elements_crc32(&s2, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3: crc32_in = front_xq_aI_mpi_elements_crc32(&s3, 1, 1, 1, size, rank, comm); break;
#endif
  }
#endif

  if (*subx != 0) for (i = 0; i < nlocal; i++) ibox[i] -= *subx;

  if (depth != NULL)
  {
    if (*depth > 0) highest = 3 * *depth - 1; else highest = -1;
  }

  INFO_CMD(
    if (I_AM_MASTER)
    {
      if (depth == NULL) printf(INFO_PRINT_PREFIX "front: starting local radix-sort (lowest 3 bit)\n");
      else printf(INFO_PRINT_PREFIX "front: starting local radix-sort (depth = %" PARAM_INTEGER_FMT " -> sorting bits [0..%" front_slint_fmt "])\n", *depth, highest);
    }
  );

  TIMING_SYNC(comm); TIMING_START(t[1]);

  if (!rebalance_only)
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0:
      if (depth == NULL) front_xqsaIl_sort_radix_iter(&s0, sx0, 1, 2, 0, -1);
      else fmm_sort_radix(front_xqsaIl_, &s0, sx0, highest, -1, -1);
      break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1:
      if (depth == NULL) front_xq_aIl_sort_radix_iter(&s1, sx1, 1, 2, 0, -1);
      else fmm_sort_radix(front_xq_aIl_, &s1, sx1, highest, -1, -1);
      break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2:
      if (depth == NULL) front_xqsaI_sort_radix_iter(&s2, sx2, 1, 2, 0, -1);
      else fmm_sort_radix(front_xqsaI_, &s2, sx2, highest, -1, -1);
      break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3:
      if (depth == NULL) front_xq_aI_sort_radix_iter(&s3, sx3, 1, 2, 0, -1);
      else fmm_sort_radix(front_xq_aI_, &s3, sx3, highest, -1, -1);
      break;
#endif
  }

  TIMING_SYNC(comm); TIMING_STOP(t[1]);

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: local radix-sort done\n");
  );

#ifdef VALIDATE
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0: l = front_xqsaIl_elements_validate_order(&s0, 1); break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1: l = front_xq_aIl_elements_validate_order(&s1, 1); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2: l = front_xqsaI_elements_validate_order(&s2, 1); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3: l = front_xq_aI_elements_validate_order(&s3, 1); break;
#endif
    default: l = 1;
  }
#endif

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0: front_xqsaIl_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1: front_xq_aIl_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2: front_xqsaI_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3: front_xq_aI_mpi_datatypes_init(); break;
#endif
  }

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: starting parallel sort (num. proc. = %d)\n", size);
  );

  scounts = malloc(size * 4 * sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

#define REDUCTION  0.25

  TIMING_SYNC(comm); TIMING_START(t[2]);

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0:
      pc0.pcm = SLPC_COUNTS_MM|SLPC_WEIGHTS_MM;
      pc0.count_min = *nmin;
      pc0.count_max = *nmax;
      pc0.weight_min = -(1.0 - *imba);
      pc0.weight_max = -(1.0 + *imba);
      SL_FMM_VAR(front_xqsaIl_SL_DEFCON(mseg.border_update_weight_reduction)) = REDUCTION;
      front_xqsaIl_mpi_partition_exact_radix(&s0, &pc0, highest, -1, 3, SL_SORTED_IN, scounts, NULL, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1:
      pc1.pcm = SLPC_COUNTS_MM|SLPC_WEIGHTS_MM;
      pc1.count_min = *nmin;
      pc1.count_max = *nmax;
      pc1.weight_min = -(1.0 - *imba);
      pc1.weight_max = -(1.0 + *imba);
      SL_FMM_VAR(front_xq_aIl_SL_DEFCON(mseg.border_update_weight_reduction)) = REDUCTION;
      front_xq_aIl_mpi_partition_exact_radix(&s1, &pc1, highest, -1, 3, SL_SORTED_IN, scounts, NULL, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2:
      pc2.pcm = SLPC_COUNTS_MM;
      pc2.count_min = *nmin;
      pc2.count_max = *nmax;
      front_xqsaI_mpi_partition_exact_radix(&s2, &pc2, highest, -1, 3, SL_SORTED_IN, scounts, NULL, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3:
      pc3.pcm = SLPC_COUNTS_MM;
      pc3.count_min = *nmin;
      pc3.count_max = *nmax;
      front_xq_aI_mpi_partition_exact_radix(&s3, &pc3, highest, -1, 3, SL_SORTED_IN, scounts, NULL, size, rank, comm);
      break;
#endif
  }

  TIMING_SYNC(comm); TIMING_STOP(t[2]);

  TIMING_SYNC(comm); TIMING_START(t[3]);

  front_(sl_MPI_Alltoall_int)(scounts, 1, rcounts, 1, comm, size, rank);

  TIMING_SYNC(comm); TIMING_STOP(t[3]);

  ALLTOALLV_PACKED_CMD(
    local_packed = ALLTOALLV_PACKED(size, nlocal);
    MPI_Allreduce(&local_packed, &global_packed, 1, MPI_INT, MPI_SUM, comm);
  );

  front_xqsaIl_counts2displs(size, scounts, sdispls);
  front_xqsaIl_counts2displs(size, rcounts, rdispls);

  nlocal = rdispls[size - 1] + rcounts[size - 1];

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(front_xqsaIl_SL_DEFCON(mea.packed)); SL_FMM_VAR(front_xqsaIl_SL_DEFCON(mea.packed)) = (global_packed > 0);
      );
      TIMING_SYNC(comm); TIMING_START(t[4]);
      front_xqsaIl_mpi_elements_alltoallv_ip(&s0, bx0, scounts, sdispls, rcounts, rdispls, size, rank, comm);
      TIMING_SYNC(comm); TIMING_STOP(t[4]);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(front_xqsaIl_SL_DEFCON(mea.packed)) = original_packed;
      );
      front_xqsaIl_elem_set_size(&s0, nlocal);
      TIMING_SYNC(comm); TIMING_START(t[5]);
      if (!rebalance_only) front_xqsaIl_mergep_2way_ip_int(&s0, sx0, size, rdispls, front_xqsaIl_merge2_memory_adaptive);
      TIMING_SYNC(comm); TIMING_STOP(t[5]);
      break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(front_xq_aIl_SL_DEFCON(mea.packed)); SL_FMM_VAR(front_xq_aIl_SL_DEFCON(mea.packed)) = (global_packed > 0);
      );
      TIMING_SYNC(comm); TIMING_START(t[4]);
      front_xq_aIl_mpi_elements_alltoallv_ip(&s1, bx1, scounts, sdispls, rcounts, rdispls, size, rank, comm);
      TIMING_SYNC(comm); TIMING_STOP(t[4]);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(front_xq_aIl_SL_DEFCON(mea.packed)) = original_packed;
      );
      front_xq_aIl_elem_set_size(&s1, nlocal);
      TIMING_SYNC(comm); TIMING_START(t[5]);
      if (!rebalance_only) front_xq_aIl_mergep_2way_ip_int(&s1, sx1, size, rdispls, front_xq_aIl_merge2_memory_adaptive);
      TIMING_SYNC(comm); TIMING_STOP(t[5]);
      break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(front_xqsaI_SL_DEFCON(mea.packed)); SL_FMM_VAR(front_xqsaI_SL_DEFCON(mea.packed)) = (global_packed > 0);
      );
      TIMING_SYNC(comm); TIMING_START(t[4]);
      front_xqsaI_mpi_elements_alltoallv_ip(&s2, bx2, scounts, sdispls, rcounts, rdispls, size, rank, comm);
      TIMING_SYNC(comm); TIMING_STOP(t[4]);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(front_xqsaI_SL_DEFCON(mea.packed)) = original_packed;
      );
      front_xqsaI_elem_set_size(&s2, nlocal);
      TIMING_SYNC(comm); TIMING_START(t[5]);
      if (!rebalance_only) front_xqsaI_mergep_2way_ip_int(&s2, sx2, size, rdispls, front_xqsaI_merge2_memory_adaptive);
      TIMING_SYNC(comm); TIMING_STOP(t[5]);
      break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(front_xq_aI_SL_DEFCON(mea.packed)); SL_FMM_VAR(front_xq_aI_SL_DEFCON(mea.packed)) = (global_packed > 0);
      );
      TIMING_SYNC(comm); TIMING_START(t[4]);
      front_xq_aI_mpi_elements_alltoallv_ip(&s3, bx3, scounts, sdispls, rcounts, rdispls, size, rank, comm);
      TIMING_SYNC(comm); TIMING_STOP(t[4]);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(front_xq_aI_SL_DEFCON(mea.packed)) = original_packed;
      );
      front_xq_aI_elem_set_size(&s3, nlocal);
      TIMING_SYNC(comm); TIMING_START(t[5]);
      if (!rebalance_only) front_xq_aI_mergep_2way_ip_int(&s3, sx3, size, rdispls, front_xq_aI_merge2_memory_adaptive);
      TIMING_SYNC(comm); TIMING_STOP(t[5]);
      break;
#endif
  }

  free(scounts);

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: parallel sort done\n");
  );

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0: front_xqsaIl_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1: front_xq_aIl_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2: front_xqsaI_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3: front_xq_aI_mpi_datatypes_release(); break;
#endif
  }

  if (*subx != 0) for (i = 0; i < nlocal; i++) ibox[i] += *subx;

  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0: front_xqsaIl_mpi_elements_get_weights(&s0, &local_weights, &total_weights, -1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1: front_xq_aIl_mpi_elements_get_weights(&s1, &local_weights, &total_weights, -1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2: local_weights = 1.0; total_weights = size; break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3: local_weights = 1.0; total_weights = size; break;
#endif
  }

  *imba = (local_weights / (total_weights / size)) - 1.0;
  *nmin = nlocal;

  DEBUG_CMD(
    printf(DEBUG_PRINT_PREFIX "front: %d: imba = %f (local_weights = %f, total_weights = %f), nmin = %" PARAM_INTEGER_FMT "\n", rank, *imba, local_weights, total_weights, *nmin);
  );

  auxmem_release(auxmem_blocks, auxmem_sizes);

#ifdef VALIDATE
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0: o = front_xqsaIl_mpi_elements_validate_order(&s0, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1: o = front_xq_aIl_mpi_elements_validate_order(&s1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2: o = front_xqsaI_mpi_elements_validate_order(&s2, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3: o = front_xq_aI_mpi_elements_validate_order(&s3, 1, size, rank, comm); break;
#endif
    default: o = 1;
  }
  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: sorting_order: %s - local_order: %s\n", (!o)?"success":"FAILED", (!l)?"success":"failed");
  );

# ifdef CHECKSUM
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsaIl
    case 0: crc32_out = front_xqsaIl_mpi_elements_crc32(&s0, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aIl
    case 1: crc32_out = front_xq_aIl_mpi_elements_crc32(&s1, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 2: crc32_out = front_xqsaI_mpi_elements_crc32(&s2, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 3: crc32_out = front_xq_aI_mpi_elements_crc32(&s3, 1, 1, 1, size, rank, comm); break;
#endif
  }
  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "front: checksum: %s (%X vs. %X)\n", (crc32_in == crc32_out)?"success":"FAILED", crc32_in, crc32_out);
  );
# endif
#endif

  TIMING_SYNC(comm); TIMING_STOP(t[0]);

  TIMING_CMD(
    if (I_AM_MASTER) printf(TIMING_PRINT_PREFIX "mpi_fmm_sort_front_part_body: %f  %f  %f  %f  %f  %f\n", t[0], t[1], t[2], t[3], t[4], t[5]);
  );
#undef I_AM_MASTER
}

#endif /* WITH_SORT_FRONT_LOAD */


void mpi_fmm_sort_front_mem(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
#ifdef WITH_SORT_FRONT_LOAD
  if (type && *type >= 2)
    mpi_fmm_sort_front_part_body(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, load, imba, nmin, nmax, 0, FCOMM_IFELSE(fcomm, NULL));
  else if (mpi_fmm_sort_front_part)
    mpi_fmm_sort_front_part_body(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, NULL, NULL, NULL, NULL, 0, FCOMM_IFELSE(fcomm, NULL));
  else
#endif
    mpi_fmm_sort_front_merge_body(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, FCOMM_IFELSE(fcomm, NULL));
}

void mpi_fmm_sort_front_mem_(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  mpi_fmm_sort_front_mem(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type
#ifdef WITH_SORT_FRONT_LOAD
    , load, imba, nmin, nmax
#endif
#ifdef WITH_FCOMM
    , fcomm
#endif
  );
}


void mpi_fmm_sort_front(
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
#ifdef WITH_SORT_FRONT_LOAD
  if (type && *type >= 2)
    mpi_fmm_sort_front_part_body(NULL, NULL, NULL, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, load, imba, nmin, nmax, 0, FCOMM_IFELSE(fcomm, NULL));
  else if (mpi_fmm_sort_front_part)
    mpi_fmm_sort_front_part_body(NULL, NULL, NULL, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, NULL, NULL, NULL, NULL, 0, FCOMM_IFELSE(fcomm, NULL));
  else
#endif
    mpi_fmm_sort_front_merge_body(NULL, NULL, NULL, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, FCOMM_IFELSE(fcomm, NULL));
}

void mpi_fmm_sort_front_(
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  mpi_fmm_sort_front(depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type
#ifdef WITH_SORT_FRONT_LOAD
    , load, imba, nmin, nmax
#endif
#ifdef WITH_FCOMM
    , fcomm
#endif
  );
}


void mpi_fmm_sort_front_3bit_mem(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
#ifdef WITH_SORT_FRONT_LOAD
  if (type && *type >= 2)
    mpi_fmm_sort_front_part_body(mem0, mem1, mem_sizes, NULL, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, load, imba, nmin, nmax, 0, FCOMM_IFELSE(fcomm, NULL));
  else if (mpi_fmm_sort_front_part)
    mpi_fmm_sort_front_part_body(mem0, mem1, mem_sizes, NULL, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, NULL, NULL, NULL, NULL, 0, FCOMM_IFELSE(fcomm, NULL));
  else
#endif
    mpi_fmm_sort_front_merge_body(mem0, mem1, mem_sizes, NULL, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, FCOMM_IFELSE(fcomm, NULL));
}

void mpi_fmm_sort_front_3bit_mem_(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  mpi_fmm_sort_front_3bit_mem(mem0, mem1, mem_sizes, subx, n, ibox, xyz, q, addr_desc, addr, scr, type
#ifdef WITH_SORT_FRONT_LOAD
    , load, imba, nmin, nmax
#endif
#ifdef WITH_FCOMM
    , fcomm
#endif
  );
}


void mpi_fmm_sort_front_3bit(
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
#ifdef WITH_SORT_FRONT_LOAD
  if (type && *type >= 2)
    mpi_fmm_sort_front_part_body(NULL, NULL, NULL, NULL, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, load, imba, nmin, nmax, 0, FCOMM_IFELSE(fcomm, NULL));
  else if (mpi_fmm_sort_front_part)
    mpi_fmm_sort_front_part_body(NULL, NULL, NULL, NULL, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, NULL, NULL, NULL, NULL, 0, FCOMM_IFELSE(fcomm, NULL));
  else
#endif
    mpi_fmm_sort_front_merge_body(NULL, NULL, NULL, NULL, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, FCOMM_IFELSE(fcomm, NULL));
}

void mpi_fmm_sort_front_3bit_(
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  mpi_fmm_sort_front_3bit(subx, n, ibox, xyz, q, addr_desc, addr, scr, type
#ifdef WITH_SORT_FRONT_LOAD
   , load, imba, nmin, nmax
#endif
#ifdef WITH_FCOMM
   , fcomm
#endif
  );
}


void mpi_fmm_sort_front_rebalance_mem(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
#ifdef WITH_SORT_FRONT_LOAD
  if (type && *type >= 2)
    mpi_fmm_sort_front_part_body(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, load, imba, nmin, nmax, 1, FCOMM_IFELSE(fcomm, NULL));
  else if (mpi_fmm_sort_front_part)
    mpi_fmm_sort_front_part_body(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, NULL, NULL, NULL, NULL, 1, FCOMM_IFELSE(fcomm, NULL));
  else
#endif
    mpi_fmm_sort_front_merge_body(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, FCOMM_IFELSE(fcomm, NULL));
}

void mpi_fmm_sort_front_rebalance_mem_(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  mpi_fmm_sort_front_rebalance_mem(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type
#ifdef WITH_SORT_FRONT_LOAD
   , load, imba, nmin, nmax
#endif
#ifdef WITH_FCOMM
   , fcomm
#endif
  );
}


void mpi_fmm_sort_front_rebalance(
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
#ifdef WITH_SORT_FRONT_LOAD
  if (type && *type >= 2)
    mpi_fmm_sort_front_part_body(NULL, NULL, NULL, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, load, imba, nmin, nmax, 1, FCOMM_IFELSE(fcomm, NULL));
  else if (mpi_fmm_sort_front_part)
    mpi_fmm_sort_front_part_body(NULL, NULL, NULL, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, NULL, NULL, NULL, NULL, 1, FCOMM_IFELSE(fcomm, NULL));
  else
#endif
    mpi_fmm_sort_front_merge_body(NULL, NULL, NULL, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type, FCOMM_IFELSE(fcomm, NULL));
}

void mpi_fmm_sort_front_rebalance_(
 pint_t *depth,
 pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type
#ifdef WITH_SORT_FRONT_LOAD
 , front_xqsaIl_sldata4_t *load, REAL_C *imba, pint_t *nmin, pint_t *nmax
#endif
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  mpi_fmm_sort_front_rebalance(depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type
#ifdef WITH_SORT_FRONT_LOAD
   , load, imba, nmin, nmax
#endif
#ifdef WITH_FCOMM
   , fcomm
#endif
  );
}

#endif /* NO_SL_FRONT */


#ifndef NO_SL_BACK

static void mpi_fmm_sort_back_merge_body(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *ntotal, pint_t *n, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, pint_t *type,
 FINT8_TYPE_C *fcomm)
{
  back_slint_t back_type, highest;

#ifndef NOT_sl_back_qxpg
  back_qxpg_elements_t s0, *sx0, smem0;
  back_qxpg_merge2x_f m20 = back_qxpg_merge2_compo_hula;  /* back_qxpg_merge2_basic_straight_01_x back_qxpg_merge2_compo_tridgell back_qxpg_merge2_compo_hula */
#endif
#ifndef NOT_sl_back_qx_g
  back_qx_g_elements_t s1, *sx1, smem1;
  back_qx_g_merge2x_f m21 = back_qx_g_merge2_compo_hula;  /* back_qx_g_merge2_basic_straight_01_x back_qx_g_merge2_compo_tridgell back_qx_g_merge2_compo_hula */
#endif
#ifndef NOT_sl_back_q_pg
  back_q_pg_elements_t s2, *sx2, smem2;
  back_q_pg_merge2x_f m22 = back_q_pg_merge2_compo_hula;  /* back_q_pg_merge2_basic_straight_01_x back_q_pg_merge2_compo_tridgell back_q_pg_merge2_compo_hula */
#endif
#ifndef NOT_sl_back_q__g
  back_q__g_elements_t s3, *sx3, smem3;
  back_q__g_merge2x_f m23 = back_q__g_merge2_compo_hula;  /* back_q__g_merge2_basic_straight_01_x back_q__g_merge2_compo_tridgell back_q__g_merge2_compo_hula */
#endif

  void *auxmem_blocks[AUXMEM_NBLOCKS];
  slint_t auxmem_sizes[AUXMEM_NBLOCKS], auxmem_max;

#ifdef DO_TIMING
  double t[3];
#endif

#ifdef VALIDATE
  back_slint_t o, l;
# ifdef CHECKSUM
  unsigned int crc32_in, crc32_out;
# endif
#endif

  MPI_Comm comm;
  int size, rank, world_rank;
#ifdef DO_MPI_INIT
  int flag;
#endif

#ifdef SL_OUTPUT_TO_FILE
  char output_file_str[32];
  FILE *sl_debug_fstream;
#endif


#ifdef DO_MPI_INIT
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(NULL, NULL);
#endif

  comm = FCOMM_IFELSE(((fcomm)?MPI_Comm_f2c(*fcomm):MPI_COMM_NULL), MPI_COMM_NULL);

  if (comm == MPI_COMM_NULL) comm = MPI_COMM_WORLD;

  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

#define I_AM_MASTER  ((comm == MPI_COMM_SELF)?(world_rank == 0):(rank == 0))

  if (type && *type >= 0) back_type = *type;
  else
  {
    if (xyz != NULL && pot != NULL) back_type = 0;
    else if (xyz != NULL && pot == NULL) back_type = 1;
    else if (xyz == NULL && pot != NULL) back_type = 2;
    else if (xyz == NULL && pot == NULL) back_type = 3;
    else back_type = -1;
  }

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: type = %" back_slint_fmt "\n", back_type);
  );

#ifdef SL_OUTPUT_TO_FILE
  sprintf(output_file_str, "sort_back.debug.%d", rank);
  back_qxpg(back_qxpg_sl_debug_fstream = )
  back_qx_g(back_qx_g_sl_debug_fstream = )
  back_q_pg(back_q_pg_sl_debug_fstream = )
  back_q__g(back_q__g_sl_debug_fstream = ) sl_debug_fstream = fopen(output_file_str, "w");
#endif

  back_qxpg(SL_FMM_VAR(back_qxpg_SL_DEFCON(mpi.rank)) = )
  back_qx_g(SL_FMM_VAR(back_qx_g_SL_DEFCON(mpi.rank)) = )
  back_q_pg(SL_FMM_VAR(back_q_pg_SL_DEFCON(mpi.rank)) = )
  back_q__g(SL_FMM_VAR(back_q__g_SL_DEFCON(mpi.rank)) = ) rank;
  
  auxmem_init(mem0, mem1, mem_sizes, auxmem_blocks, auxmem_sizes, &auxmem_max);

  if (auxmem_sizes[auxmem_max] > 0)
  {
    switch (back_type)
    {
#ifndef NOT_sl_back_qxpg
      case 0:
        back_qxpg_elements_alloc_from_blocks(&smem0, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx0 = &smem0;
        break;
#endif
#ifndef NOT_sl_back_qx_g
      case 1:
        back_qx_g_elements_alloc_from_blocks(&smem1, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx1 = &smem1;
        break;
#endif
#ifndef NOT_sl_back_q_pg
      case 2:
        back_q_pg_elements_alloc_from_blocks(&smem2, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx2 = &smem2;
        break;
#endif
#ifndef NOT_sl_back_q__g
      case 3:
        back_q__g_elements_alloc_from_blocks(&smem3, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx3 = &smem3;
        break;
#endif
    }

  } else
  {
    back_qxpg(sx0 = NULL;)
    back_qx_g(sx1 = NULL;)
    back_q_pg(sx2 = NULL;)
    back_q__g(sx3 = NULL;)
  }

  INFO_CMD(
    if (I_AM_MASTER)
    {
      printf(INFO_PRINT_PREFIX "back: sizeof(back_slint_t) = %d bytes - #elements: %" PARAM_INTEGER_FMT " - mem0: %p / %" PARAM_INTEGER_FMT " bytes, mem1: %p / %" PARAM_INTEGER_FMT " bytes\n", (int) sizeof(back_(slint_t)), *n, mem0, (mem_sizes)?mem_sizes[0]:0, mem1, (mem_sizes)?mem_sizes[1]:0);
      printf(INFO_PRINT_PREFIX "back: addr:      sizeof(back_slkey_t)   = %d byte\n", (int) sizeof(back_(slkey_t)));
      printf(INFO_PRINT_PREFIX "back:    q: %d * sizeof(back_sldata0_t) = %d * %d byte\n", (int) back_(sl_data0_size_c), (int) back_(sl_data0_size_c), (int) sizeof(back_(sldata0_t)));
      printf(INFO_PRINT_PREFIX "back:  xyz: %d * sizeof(back_sldata1_t) = %d * %d byte\n", (int) back_(sl_data1_size_c), (int) back_(sl_data1_size_c), (int) sizeof(back_(sldata1_t)));
      printf(INFO_PRINT_PREFIX "back:  pot: %d * sizeof(back_sldata2_t) = %d * %d byte\n", (int) back_(sl_data2_size_c), (int) back_(sl_data2_size_c), (int) sizeof(back_(sldata2_t)));
      printf(INFO_PRINT_PREFIX "back: grad: %d * sizeof(back_sldata3_t) = %d * %d byte\n", (int) back_(sl_data3_size_c), (int) back_(sl_data3_size_c), (int) sizeof(back_(sldata3_t)));
    }
  );

  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0:
      back_qxpg_elem_set_size(&s0, *n);
      back_qxpg_elem_set_max_size(&s0, *n);
      back_qxpg_elem_set_keys(&s0, addr);
      back_qxpg_elem_set_data(&s0, q, xyz, pot, grad);
      break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1:
      back_qx_g_elem_set_size(&s1, *n);
      back_qx_g_elem_set_max_size(&s1, *n);
      back_qx_g_elem_set_keys(&s1, addr);
      back_qx_g_elem_set_data(&s1, q, xyz, grad);
      break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2:
      back_q_pg_elem_set_size(&s2, *n);
      back_q_pg_elem_set_max_size(&s2, *n);
      back_q_pg_elem_set_keys(&s2, addr);
      back_q_pg_elem_set_data(&s2, q, pot, grad);
      break;
#endif
#ifndef NOT_sl_back_q__g
    case 3:
      back_q__g_elem_set_size(&s3, *n);
      back_q__g_elem_set_max_size(&s3, *n);
      back_q__g_elem_set_keys(&s3, addr);
      back_q__g_elem_set_data(&s3, q, grad);
      break;
#endif
  }

#if defined(VALIDATE) && defined(CHECKSUM)
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: crc32_in = back_qxpg_mpi_elements_crc32(&s0, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: crc32_in = back_qx_g_mpi_elements_crc32(&s1, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: crc32_in = back_q_pg_mpi_elements_crc32(&s2, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: crc32_in = back_q__g_mpi_elements_crc32(&s3, 1, 1, 1, size, rank, comm); break;
#endif
  }
#endif

#ifdef SORT_BACK_ONE2NCHARGES
  highest = mpi_log2_floor(*ntotal);
#else
  highest = mpi_log2_floor(*ntotal - 1);
#endif

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: starting local radix-sort (total ncharges = %" PARAM_INTEGER_FMT " -> sorting bits [0..%" back_slint_fmt "])\n", *ntotal, highest);
  );

  TIMING_SYNC(comm); TIMING_START(t[1]);

  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: fmm_sort_radix(back_qxpg_, &s0, sx0, highest, -1, -1); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: fmm_sort_radix(back_qx_g_, &s1, sx1, highest, -1, -1); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: fmm_sort_radix(back_q_pg_, &s2, sx2, highest, -1, -1); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: fmm_sort_radix(back_q__g_, &s3, sx3, highest, -1, -1); break;
#endif
  }

  TIMING_SYNC(comm); TIMING_STOP(t[1]);

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: local radix-sort done\n");
  );

#ifdef VALIDATE
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: l = back_qxpg_elements_validate_order(&s0, 1); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: l = back_qx_g_elements_validate_order(&s1, 1); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: l = back_q_pg_elements_validate_order(&s2, 1); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: l = back_q__g_elements_validate_order(&s3, 1); break;
#endif
    default: l = 1;
  }
#endif

  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: back_qxpg_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: back_qx_g_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: back_q_pg_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: back_q__g_mpi_datatypes_init(); break;
#endif
  }

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: starting parallel sort (num. proc. = %d)\n", size);
  );

#ifdef MPI_SENDRECV_REPLACE_MAX_SIZE
  back_qxpg(SL_FMM_VAR(back_qxpg_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = )
  back_qx_g(SL_FMM_VAR(back_qx_g_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = )
  back_q_pg(SL_FMM_VAR(back_q_pg_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = )
  back_q__g(SL_FMM_VAR(back_q__g_SL_DEFCON(me.sendrecv_replace_mpi_maxsize)) = ) MPI_SENDRECV_REPLACE_MAX_SIZE;
#endif

  TIMING_SYNC(comm); TIMING_START(t[2]);

  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(back_qxpg_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(back_qxpg_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m20 = back_qxpg_merge2_memory_adaptive;
      }
      back_qxpg_mpi_mergek(&s0, back_qxpg_sn_batcher, NULL, m20, sx0, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(back_qx_g_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(back_qx_g_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m21 = back_qx_g_merge2_memory_adaptive;
      }
      back_qx_g_mpi_mergek(&s1, back_qx_g_sn_batcher, NULL, m21, sx1, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(back_q_pg_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(back_q_pg_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m22 = back_q_pg_merge2_memory_adaptive;
      }
      back_q_pg_mpi_mergek(&s2, back_q_pg_sn_batcher, NULL, m22, sx2, size, rank, comm);
      break;
#endif
#ifndef NOT_sl_back_q__g
    case 3:
      if (auxmem_sizes[auxmem_max] > 0)
      {
        SL_FMM_VAR(back_q__g_SL_DEFCON(me.sendrecv_replace_mem)) = auxmem_blocks[auxmem_max];
        SL_FMM_VAR(back_q__g_SL_DEFCON(me.sendrecv_replace_memsize)) = auxmem_sizes[auxmem_max];

        m23 = back_q__g_merge2_memory_adaptive;
      }
      back_q__g_mpi_mergek(&s3, back_q__g_sn_batcher, NULL, m23, sx3, size, rank, comm);
      break;
#endif
  }

  TIMING_SYNC(comm); TIMING_STOP(t[2]);

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: parallel sort done\n");
  );

  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: back_qxpg_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: back_qx_g_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: back_q_pg_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: back_q__g_mpi_datatypes_release(); break;
#endif
  }

  auxmem_release(auxmem_blocks, auxmem_sizes);

#ifdef VALIDATE
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: o = back_qxpg_mpi_elements_validate_order(&s0, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: o = back_qx_g_mpi_elements_validate_order(&s1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: o = back_q_pg_mpi_elements_validate_order(&s2, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: o = back_q__g_mpi_elements_validate_order(&s3, 1, size, rank, comm); break;
#endif
    default: o = 1;
  }
  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: sorting_order: %s - local_order: %s\n", (!o)?"success":"FAILED", (!l)?"success":"failed");
  );

# ifdef CHECKSUM
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: crc32_out = back_qxpg_mpi_elements_crc32(&s0, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: crc32_out = back_qx_g_mpi_elements_crc32(&s1, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: crc32_out = back_q_pg_mpi_elements_crc32(&s2, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: crc32_out = back_q__g_mpi_elements_crc32(&s3, 1, 1, 1, size, rank, comm); break;
#endif
  }
  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: checksum: %s (%X vs. %X)\n", (crc32_in == crc32_out)?"success":"FAILED", crc32_in, crc32_out);
  );
# endif
#endif

#ifdef SL_OUTPUT_TO_FILE
  fclose(sl_debug_fstream);
  back_qxpg(back_qxpg_sl_debug_fstream = )
  back_qx_g(back_qx_g_sl_debug_fstream = )
  back_q_pg(back_q_pg_sl_debug_fstream = )
  back_q__g(back_q__g_sl_debug_fstream = ) sl_debug_fstream = NULL;
#endif

  TIMING_SYNC(comm); TIMING_STOP(t[0]);

  TIMING_CMD(
    if (I_AM_MASTER) printf(TIMING_PRINT_PREFIX "mpi_fmm_sort_back_merge_body: %f  %f  %f\n", t[0], t[1], t[2]);
  );
#undef I_AM_MASTER
}


static void mpi_fmm_sort_back_part_body(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *ntotal, pint_t *nin, pint_t *nout, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, back_(sldata4_t) *load, pint_t *type,
 FINT8_TYPE_C *fcomm)
{
  pint_t base_addr;
  back_(slkey_t) lh_addrs[2];

  back_slint_t back_type, nlocal, max_nlocal, sxalloc;

#ifndef NOT_sl_back_qxpg
  back_qxpg_elements_t s0, *sx0, smem0, *bx0, bmem0;
#endif
#ifndef NOT_sl_back_qx_g
  back_qx_g_elements_t s1, *sx1, smem1, *bx1, bmem1;
#endif
#ifndef NOT_sl_back_q_pg
  back_q_pg_elements_t s2, *sx2, smem2, *bx2, bmem2;
#endif
#ifndef NOT_sl_back_q__g
  back_q__g_elements_t s3, *sx3, smem3, *bx3, bmem3;
#endif
#ifndef NOT_sl_back_qxpgl
  back_qxpgl_elements_t s4, *sx4, smem4, *bx4, bmem4;
#endif
#ifndef NOT_sl_back_qx_gl
  back_qx_gl_elements_t s5, *sx5, smem5, *bx5, bmem5;
#endif
#ifndef NOT_sl_back_q_pgl
  back_q_pgl_elements_t s6, *sx6, smem6, *bx6, bmem6;
#endif
#ifndef NOT_sl_back_q__gl
  back_q__gl_elements_t s7, *sx7, smem7, *bx7, bmem7;
#endif

  void *auxmem_blocks[AUXMEM_NBLOCKS];
  slint_t auxmem_sizes[AUXMEM_NBLOCKS], auxmem_max;

#ifdef ALLTOALLV_PACKED
  int local_packed, global_packed, original_packed;
#endif

#ifdef DO_TIMING
  double t[2], *mssp_b_t;
#endif

#ifdef VALIDATE
  back_slint_t o;
# ifdef CHECKSUM
  unsigned int crc32_in, crc32_out;
# endif
#endif

  MPI_Comm comm;
  int size, rank, world_rank;
#ifdef DO_MPI_INIT
  int flag;
#endif

#ifdef SL_OUTPUT_TO_FILE
  char output_file_str[32];
  FILE *sl_debug_fstream;
#endif


#ifdef DO_MPI_INIT
  MPI_Initialized(&flag);
  if (!flag) MPI_Init(NULL, NULL);
#endif

  comm = FCOMM_IFELSE(((fcomm)?MPI_Comm_f2c(*fcomm):MPI_COMM_NULL), MPI_COMM_NULL);

  if (comm == MPI_COMM_NULL) comm = MPI_COMM_WORLD;

  TIMING_SYNC(comm); TIMING_START(t[0]);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

#define I_AM_MASTER  ((comm == MPI_COMM_SELF)?(world_rank == 0):(rank == 0))

  if (type && *type >= 0) back_type = *type;
  else
  {
    if (xyz != NULL && pot != NULL) back_type = 0;
    else if (xyz != NULL && pot == NULL) back_type = 1;
    else if (xyz == NULL && pot != NULL) back_type = 2;
    else if (xyz == NULL && pot == NULL) back_type = 3;
    else back_type = -1;

    if (load != NULL && back_type >= 0) back_type += 4;
  }

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: type = %" back_slint_fmt "\n", back_type);
  );

  if (*nout < 0) *nout = *nin;

#ifdef SL_OUTPUT_TO_FILE
  sprintf(output_file_str, "sort_back.debug.%d", rank);
  back_qxpg(back_qxpg_sl_debug_fstream = )
  back_qx_g(back_qx_g_sl_debug_fstream = )
  back_q_pg(back_q_pg_sl_debug_fstream = )
  back_q__g(back_q__g_sl_debug_fstream = )
  back_qxpgl(back_qxpgl_sl_debug_fstream = )
  back_qx_gl(back_qx_gl_sl_debug_fstream = )
  back_q_pgl(back_q_pgl_sl_debug_fstream = )
  back_q__gl(back_q__gl_sl_debug_fstream = ) sl_debug_fstream = fopen(output_file_str, "w");
#endif

  back_qxpg(SL_FMM_VAR(back_qxpg_SL_DEFCON(mpi.rank)) = )
  back_qx_g(SL_FMM_VAR(back_qx_g_SL_DEFCON(mpi.rank)) = )
  back_q_pg(SL_FMM_VAR(back_q_pg_SL_DEFCON(mpi.rank)) = )
  back_q__g(SL_FMM_VAR(back_q__g_SL_DEFCON(mpi.rank)) = )
  back_qxpgl(SL_FMM_VAR(back_qxpgl_SL_DEFCON(mpi.rank)) = )
  back_qx_gl(SL_FMM_VAR(back_qx_gl_SL_DEFCON(mpi.rank)) = )
  back_q_pgl(SL_FMM_VAR(back_q_pgl_SL_DEFCON(mpi.rank)) = )
  back_q__gl(SL_FMM_VAR(back_q__gl_SL_DEFCON(mpi.rank)) = ) rank;

  sxalloc = 0;
  
  auxmem_init(mem0, mem1, mem_sizes, auxmem_blocks, auxmem_sizes, &auxmem_max);

  if (auxmem_sizes[auxmem_max] > 0)
  {
    switch (back_type)
    {
#ifndef NOT_sl_back_qxpg
      case 0:
        back_qxpg_elements_alloc_from_blocks(&smem0, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx0 = &smem0;
        back_qxpg_elem_set_block(&bmem0, auxmem_blocks[auxmem_max]);
        back_qxpg_elem_set_block_size(&bmem0, auxmem_sizes[auxmem_max]);
        bx0 = &bmem0;
        break;
#endif
#ifndef NOT_sl_back_qx_g
      case 1:
        back_qx_g_elements_alloc_from_blocks(&smem1, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx1 = &smem1;
        back_qx_g_elem_set_block(&bmem1, auxmem_blocks[auxmem_max]);
        back_qx_g_elem_set_block_size(&bmem1, auxmem_sizes[auxmem_max]);
        bx1 = &bmem1;
        break;
#endif
#ifndef NOT_sl_back_q_pg
      case 2:
        back_q_pg_elements_alloc_from_blocks(&smem2, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx2 = &smem2;
        back_q_pg_elem_set_block(&bmem2, auxmem_blocks[auxmem_max]);
        back_q_pg_elem_set_block_size(&bmem2, auxmem_sizes[auxmem_max]);
        bx2 = &bmem2;
        break;
#endif
#ifndef NOT_sl_back_q__g
      case 3:
        back_q__g_elements_alloc_from_blocks(&smem3, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx3 = &smem3;
        back_q__g_elem_set_block(&bmem3, auxmem_blocks[auxmem_max]);
        back_q__g_elem_set_block_size(&bmem3, auxmem_sizes[auxmem_max]);
        bx3 = &bmem3;
        break;
#endif
#ifndef NOT_sl_back_qxpgl
      case 4:
        back_qxpgl_elements_alloc_from_blocks(&smem4, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx4 = &smem4;
        back_qxpgl_elem_set_block(&bmem4, auxmem_blocks[auxmem_max]);
        back_qxpgl_elem_set_block_size(&bmem4, auxmem_sizes[auxmem_max]);
        bx4 = &bmem4;
        break;
#endif
#ifndef NOT_sl_back_qx_gl
      case 5:
        back_qx_gl_elements_alloc_from_blocks(&smem5, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx5 = &smem5;
        back_qx_gl_elem_set_block(&bmem5, auxmem_blocks[auxmem_max]);
        back_qx_gl_elem_set_block_size(&bmem5, auxmem_sizes[auxmem_max]);
        bx5 = &bmem5;
        break;
#endif
#ifndef NOT_sl_back_q_pgl
      case 6:
        back_q_pgl_elements_alloc_from_blocks(&smem6, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx6 = &smem6;
        back_q_pgl_elem_set_block(&bmem6, auxmem_blocks[auxmem_max]);
        back_q_pgl_elem_set_block_size(&bmem6, auxmem_sizes[auxmem_max]);
        bx6 = &bmem6;
        break;
#endif
#ifndef NOT_sl_back_q__gl
      case 7:
        back_q__gl_elements_alloc_from_blocks(&smem7, AUXMEM_NBLOCKS, auxmem_blocks, auxmem_sizes, 8, -1, SLCM_KEYS|SLCM_DATA);
        sx7 = &smem7;
        back_q__gl_elem_set_block(&bmem7, auxmem_blocks[auxmem_max]);
        back_q__gl_elem_set_block_size(&bmem7, auxmem_sizes[auxmem_max]);
        bx7 = &bmem7;
        break;
#endif
    }

  } else
  {
    back_qxpg(sx0 = bx0 = NULL;)
    back_qx_g(sx1 = bx1 = NULL;)
    back_q_pg(sx2 = bx2 = NULL;)
    back_q__g(sx3 = bx3 = NULL;)
    back_qxpgl(sx4 = bx4 = NULL;)
    back_qx_gl(sx5 = bx5 = NULL;)
    back_q_pgl(sx6 = bx6 = NULL;)
    back_q__gl(sx7 = bx7 = NULL;)

    sxalloc = 1;
  
    switch (back_type)
    {
#ifndef NOT_sl_back_qxpg
      case 0:
        back_qxpg_elements_alloc(&smem0, z_max(*nin, *nout), SLCM_KEYS|SLCM_DATA);
        sx0 = &smem0;
        break;
#endif
#ifndef NOT_sl_back_qx_g
      case 1:
        back_qx_g_elements_alloc(&smem1, z_max(*nin, *nout), SLCM_KEYS|SLCM_DATA);
        sx1 = &smem1;
        break;
#endif
#ifndef NOT_sl_back_q_pg
      case 2:
        back_q_pg_elements_alloc(&smem2, z_max(*nin, *nout), SLCM_KEYS|SLCM_DATA);
        sx2 = &smem2;
        break;
#endif
#ifndef NOT_sl_back_q__g
      case 3:
        back_q__g_elements_alloc(&smem3, z_max(*nin, *nout), SLCM_KEYS|SLCM_DATA);
        sx3 = &smem3;
        break;
#endif
#ifndef NOT_sl_back_qxpgl
      case 4:
        back_qxpgl_elements_alloc(&smem4, z_max(*nin, *nout), SLCM_KEYS|SLCM_DATA);
        sx4 = &smem4;
        break;
#endif
#ifndef NOT_sl_back_qx_gl
      case 5:
        back_qx_gl_elements_alloc(&smem5, z_max(*nin, *nout), SLCM_KEYS|SLCM_DATA);
        sx5 = &smem5;
        break;
#endif
#ifndef NOT_sl_back_q_pgl
      case 6:
        back_q_pgl_elements_alloc(&smem6, z_max(*nin, *nout), SLCM_KEYS|SLCM_DATA);
        sx6 = &smem6;
        break;
#endif
#ifndef NOT_sl_back_q__gl
      case 7:
        back_q__gl_elements_alloc(&smem7, z_max(*nin, *nout), SLCM_KEYS|SLCM_DATA);
        sx7 = &smem7;
        break;
#endif
    }
  }

  INFO_CMD(
    if (I_AM_MASTER)
    {
      printf(INFO_PRINT_PREFIX "back: sizeof(back_slint_t) = %d bytes - #elements in: %" PARAM_INTEGER_FMT " - #elements out: %" PARAM_INTEGER_FMT " - mem0: %p / %" PARAM_INTEGER_FMT " bytes, mem1: %p / %" PARAM_INTEGER_FMT " bytes\n", (int) sizeof(back_(slint_t)), *nin, *nout, mem0, (mem_sizes)?mem_sizes[0]:0, mem1, (mem_sizes)?mem_sizes[1]:0);
      printf(INFO_PRINT_PREFIX "back: addr:      sizeof(back_slkey_t)   = %d byte\n", (int) sizeof(back_(slkey_t)));
      printf(INFO_PRINT_PREFIX "back:    q: %d * sizeof(back_sldata0_t) = %d * %d byte\n", (int) back_(sl_data0_size_c), (int) back_(sl_data0_size_c), (int) sizeof(back_(sldata0_t)));
      printf(INFO_PRINT_PREFIX "back:  xyz: %d * sizeof(back_sldata1_t) = %d * %d byte\n", (int) back_(sl_data1_size_c), (int) back_(sl_data1_size_c), (int) sizeof(back_(sldata1_t)));
      printf(INFO_PRINT_PREFIX "back:  pot: %d * sizeof(back_sldata2_t) = %d * %d byte\n", (int) back_(sl_data2_size_c), (int) back_(sl_data2_size_c), (int) sizeof(back_(sldata2_t)));
      printf(INFO_PRINT_PREFIX "back: grad: %d * sizeof(back_sldata3_t) = %d * %d byte\n", (int) back_(sl_data3_size_c), (int) back_(sl_data3_size_c), (int) sizeof(back_(sldata3_t)));
      printf(INFO_PRINT_PREFIX "back: load: %d * sizeof(back_sldata4_t) = %d * %d byte\n", (int) back_(sl_data4_size_c), (int) back_(sl_data4_size_c), (int) sizeof(back_(sldata4_t)));
    }
  );

  nlocal = *nin;
  max_nlocal = z_max(*nin, *nout);

  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0:
      back_qxpg_elem_set_size(&s0, nlocal);
      back_qxpg_elem_set_max_size(&s0, max_nlocal);
      back_qxpg_elem_set_keys(&s0, addr);
      back_qxpg_elem_set_data(&s0, q, xyz, pot, grad);
      break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1:
      back_qx_g_elem_set_size(&s1, nlocal);
      back_qx_g_elem_set_max_size(&s1, max_nlocal);
      back_qx_g_elem_set_keys(&s1, addr);
      back_qx_g_elem_set_data(&s1, q, xyz, grad);
      break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2:
      back_q_pg_elem_set_size(&s2, nlocal);
      back_q_pg_elem_set_max_size(&s2, max_nlocal);
      back_q_pg_elem_set_keys(&s2, addr);
      back_q_pg_elem_set_data(&s2, q, pot, grad);
      break;
#endif
#ifndef NOT_sl_back_q__g
    case 3:
      back_q__g_elem_set_size(&s3, nlocal);
      back_q__g_elem_set_max_size(&s3, max_nlocal);
      back_q__g_elem_set_keys(&s3, addr);
      back_q__g_elem_set_data(&s3, q, grad);
      break;
#endif
#ifndef NOT_sl_back_qxpgl
    case 4:
      back_qxpgl_elem_set_size(&s4, nlocal);
      back_qxpgl_elem_set_max_size(&s4, max_nlocal);
      back_qxpgl_elem_set_keys(&s4, addr);
      back_qxpgl_elem_set_data(&s4, q, xyz, pot, grad, load);
      break;
#endif
#ifndef NOT_sl_back_qx_gl
    case 5:
      back_qx_gl_elem_set_size(&s5, nlocal);
      back_qx_gl_elem_set_max_size(&s5, max_nlocal);
      back_qx_gl_elem_set_keys(&s5, addr);
      back_qx_gl_elem_set_data(&s5, q, xyz, grad, load);
      break;
#endif
#ifndef NOT_sl_back_q_pgl
    case 6:
      back_q_pgl_elem_set_size(&s6, nlocal);
      back_q_pgl_elem_set_max_size(&s6, max_nlocal);
      back_q_pgl_elem_set_keys(&s6, addr);
      back_q_pgl_elem_set_data(&s6, q, pot, grad, load);
      break;
#endif
#ifndef NOT_sl_back_q__gl
    case 7:
      back_q__gl_elem_set_size(&s7, nlocal);
      back_q__gl_elem_set_max_size(&s7, max_nlocal);
      back_q__gl_elem_set_keys(&s7, addr);
      back_q__gl_elem_set_data(&s7, q, grad, load);
      break;
#endif
  }

#if defined(VALIDATE) && defined(CHECKSUM)
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: crc32_in = back_qxpg_mpi_elements_crc32(&s0, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: crc32_in = back_qx_g_mpi_elements_crc32(&s1, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: crc32_in = back_q_pg_mpi_elements_crc32(&s2, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: crc32_in = back_q__g_mpi_elements_crc32(&s3, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qxpgl
    case 4: crc32_in = back_qxpgl_mpi_elements_crc32(&s4, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qx_gl
    case 5: crc32_in = back_qx_gl_mpi_elements_crc32(&s5, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q_pgl
    case 6: crc32_in = back_q_pgl_mpi_elements_crc32(&s6, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q__gl
    case 7: crc32_in = back_q__gl_mpi_elements_crc32(&s7, 1, 1, 1, size, rank, comm); break;
#endif
  }
#endif

  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: back_qxpg_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: back_qx_g_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: back_q_pg_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: back_q__g_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_qxpgl
    case 4: back_qxpgl_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_qx_gl
    case 5: back_qx_gl_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_q_pgl
    case 6: back_q_pgl_mpi_datatypes_init(); break;
#endif
#ifndef NOT_sl_back_q__gl
    case 7: back_q__gl_mpi_datatypes_init(); break;
#endif
  }

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: starting parallel back sort (num. proc. = %d)\n", size);
  );

  TIMING_SYNC(comm); TIMING_START(t[1]);

  ALLTOALLV_PACKED_CMD(
    local_packed = ALLTOALLV_PACKED(size, nlocal);
    MPI_Allreduce(&local_packed, &global_packed, 1, MPI_INT, MPI_SUM, comm);
  );

  base_addr = 0;
  MPI_Exscan(nout, &base_addr, 1, PARAM_INTEGER_MPI, MPI_SUM, comm);

  lh_addrs[0] = base_addr;
  lh_addrs[1] = base_addr + *nout - 1;

  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(back_qxpg_SL_DEFCON(mssp.back_packed)); SL_FMM_VAR(back_qxpg_SL_DEFCON(mssp.back_packed)) = (global_packed > 0);
      );
      back_qxpg_mpi_sort_back(&s0, NULL, bx0, lh_addrs, -1, size, rank, comm);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(back_qxpg_SL_DEFCON(mssp.back_packed)) = original_packed;
      );
      break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(back_qx_g_SL_DEFCON(mssp.back_packed)); SL_FMM_VAR(back_qx_g_SL_DEFCON(mssp.back_packed)) = (global_packed > 0);
      );
      back_qx_g_mpi_sort_back(&s1, NULL, bx1, lh_addrs, -1, size, rank, comm);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(back_qx_g_SL_DEFCON(mssp.back_packed)) = original_packed;
      );
      break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(back_q_pg_SL_DEFCON(mssp.back_packed)); SL_FMM_VAR(back_q_pg_SL_DEFCON(mssp.back_packed)) = (global_packed > 0);
      );
      back_q_pg_mpi_sort_back(&s2, NULL, bx2, lh_addrs, -1, size, rank, comm);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(back_q_pg_SL_DEFCON(mssp.back_packed)) = original_packed;
      );
      break;
#endif
#ifndef NOT_sl_back_q__g
    case 3:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(back_q__g_SL_DEFCON(mssp.back_packed)); SL_FMM_VAR(back_q__g_SL_DEFCON(mssp.back_packed)) = (global_packed > 0);
      );
      back_q__g_mpi_sort_back(&s3, NULL, bx3, lh_addrs, -1, size, rank, comm);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(back_q__g_SL_DEFCON(mssp.back_packed)) = original_packed;
      );
      break;
#endif
#ifndef NOT_sl_back_qxpgl
    case 4:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(back_qxpgl_SL_DEFCON(mssp.back_packed)); SL_FMM_VAR(back_qxpgl_SL_DEFCON(mssp.back_packed)) = (global_packed > 0);
      );
      back_qxpgl_mpi_sort_back(&s4, NULL, bx4, lh_addrs, -1, size, rank, comm);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(back_qxpgl_SL_DEFCON(mssp.back_packed)) = original_packed;
      );
      break;
#endif
#ifndef NOT_sl_back_qx_gl
    case 5:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(back_qx_gl_SL_DEFCON(mssp.back_packed)); SL_FMM_VAR(back_qx_gl_SL_DEFCON(mssp.back_packed)) = (global_packed > 0);
      );
      back_qx_gl_mpi_sort_back(&s5, NULL, bx5, lh_addrs, -1, size, rank, comm);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(back_qx_gl_SL_DEFCON(mssp.back_packed)) = original_packed;
      );
      break;
#endif
#ifndef NOT_sl_back_q_pgl
    case 6:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(back_q_pgl_SL_DEFCON(mssp.back_packed)); SL_FMM_VAR(back_q_pgl_SL_DEFCON(mssp.back_packed)) = (global_packed > 0);
      );
      back_q_pgl_mpi_sort_back(&s6, NULL, bx6, lh_addrs, -1, size, rank, comm);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(back_q_pgl_SL_DEFCON(mssp.back_packed)) = original_packed;
      );
      break;
#endif
#ifndef NOT_sl_back_q__gl
    case 7:
      ALLTOALLV_PACKED_CMD(
        original_packed = SL_FMM_VAR(back_q__gl_SL_DEFCON(mssp.back_packed)); SL_FMM_VAR(back_q__gl_SL_DEFCON(mssp.back_packed)) = (global_packed > 0);
      );
      back_q__gl_mpi_sort_back(&s7, NULL, bx7, lh_addrs, -1, size, rank, comm);
      ALLTOALLV_PACKED_CMD(
        SL_FMM_VAR(back_q__gl_SL_DEFCON(mssp.back_packed)) = original_packed;
      );
      break;
#endif
  }

  TIMING_SYNC(comm); TIMING_STOP(t[1]);

  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: parallel sort back done\n");
  );

  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: back_qxpg_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: back_qx_g_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: back_q_pg_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: back_q__g_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_qxpgl
    case 4: back_qxpgl_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_qx_gl
    case 5: back_qx_gl_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_q_pgl
    case 6: back_q_pgl_mpi_datatypes_release(); break;
#endif
#ifndef NOT_sl_back_q__gl
    case 7: back_q__gl_mpi_datatypes_release(); break;
#endif
  }

  if (sxalloc)
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: back_qxpg_elements_free(sx0); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: back_qx_g_elements_free(sx1); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: back_q_pg_elements_free(sx2); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: back_q__g_elements_free(sx3); break;
#endif
#ifndef NOT_sl_back_qxpgl
    case 4: back_qxpgl_elements_free(sx4); break;
#endif
#ifndef NOT_sl_back_qx_gl
    case 5: back_qx_gl_elements_free(sx5); break;
#endif
#ifndef NOT_sl_back_q_pgl
    case 6: back_q_pgl_elements_free(sx6); break;
#endif
#ifndef NOT_sl_back_q__gl
    case 7: back_q__gl_elements_free(sx7); break;
#endif
  }

  auxmem_release(auxmem_blocks, auxmem_sizes);

#ifdef VALIDATE
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: o = back_qxpg_mpi_elements_validate_order(&s0, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: o = back_qx_g_mpi_elements_validate_order(&s1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: o = back_q_pg_mpi_elements_validate_order(&s2, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: o = back_q__g_mpi_elements_validate_order(&s3, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qxpgl
    case 4: o = back_qxpgl_mpi_elements_validate_order(&s4, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qx_gl
    case 5: o = back_qx_gl_mpi_elements_validate_order(&s5, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q_pgl
    case 6: o = back_q_pgl_mpi_elements_validate_order(&s6, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q__gl
    case 7: o = back_q__gl_mpi_elements_validate_order(&s7, 1, size, rank, comm); break;
#endif
    default: o = 1;
  }
  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: sorting_order: %s\n", (!o)?"success":"FAILED");
  );

# ifdef CHECKSUM
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: crc32_out = back_qxpg_mpi_elements_crc32(&s0, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: crc32_out = back_qx_g_mpi_elements_crc32(&s1, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: crc32_out = back_q_pg_mpi_elements_crc32(&s2, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: crc32_out = back_q__g_mpi_elements_crc32(&s3, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qxpgl
    case 4: crc32_out = back_qxpgl_mpi_elements_crc32(&s4, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_qx_gl
    case 5: crc32_out = back_qx_gl_mpi_elements_crc32(&s5, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q_pgl
    case 6: crc32_out = back_q_pgl_mpi_elements_crc32(&s6, 1, 1, 1, size, rank, comm); break;
#endif
#ifndef NOT_sl_back_q__gl
    case 7: crc32_out = back_q__gl_mpi_elements_crc32(&s7, 1, 1, 1, size, rank, comm); break;
#endif
  }
  INFO_CMD(
    if (I_AM_MASTER) printf(INFO_PRINT_PREFIX "back: checksum: %s (%X vs. %X)\n", (crc32_in == crc32_out)?"success":"FAILED", crc32_in, crc32_out);
  );
# endif
#endif

#ifdef SL_OUTPUT_TO_FILE
  fclose(sl_debug_fstream);
  back_qxpg(back_qxpg_sl_debug_fstream = )
  back_qx_g(back_qx_g_sl_debug_fstream = )
  back_q_pg(back_q_pg_sl_debug_fstream = )
  back_q__g(back_q__g_sl_debug_fstream = )
  back_qxpgl(back_qxpgl_sl_debug_fstream = )
  back_qx_gl(back_qx_gl_sl_debug_fstream = )
  back_q_pgl(back_q_pgl_sl_debug_fstream = )
  back_q__gl(back_q__gl_sl_debug_fstream = ) sl_debug_fstream = NULL;
#endif

  TIMING_SYNC(comm); TIMING_STOP(t[0]);

#ifdef DO_TIMING
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: mssp_b_t = SL_FMM_VAR(back_qxpg_SL_DEFCON(mssp.b_t)); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: mssp_b_t = SL_FMM_VAR(back_qx_g_SL_DEFCON(mssp.b_t)); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: mssp_b_t = SL_FMM_VAR(back_q_pg_SL_DEFCON(mssp.b_t)); break;;
#endif
#ifndef NOT_sl_back_q__g
    case 3: mssp_b_t = SL_FMM_VAR(back_q__g_SL_DEFCON(mssp.b_t)); break;
#endif
#ifndef NOT_sl_back_qxpgl
    case 4: mssp_b_t = SL_FMM_VAR(back_qxpgl_SL_DEFCON(mssp.b_t)); break;
#endif
#ifndef NOT_sl_back_qx_gl
    case 5: mssp_b_t = SL_FMM_VAR(back_qx_gl_SL_DEFCON(mssp.b_t)); break;
#endif
#ifndef NOT_sl_back_q_pgl
    case 6: mssp_b_t = SL_FMM_VAR(back_q_pgl_SL_DEFCON(mssp.b_t)); break;;
#endif
#ifndef NOT_sl_back_q__gl
    case 7: mssp_b_t = SL_FMM_VAR(back_q__gl_SL_DEFCON(mssp.b_t)); break;
#endif
  }
#endif

  TIMING_CMD(
    if (I_AM_MASTER) printf(TIMING_PRINT_PREFIX "mpi_fmm_sort_back_part_body: %f  %f  %f  %f  %f\n", t[0], t[1], mssp_b_t[0], mssp_b_t[1], mssp_b_t[2]);
  );
#undef I_AM_MASTER
}


void mpi_fmm_sort_back_mem(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *ntotal, pint_t *nin, pint_t *nout, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, back_(sldata4_t) *load, pint_t *type
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  if (*nout >= 0 || mpi_fmm_sort_back_part)
    mpi_fmm_sort_back_part_body(mem0, mem1, mem_sizes, ntotal, nin, nout, addr, q, xyz, pot, grad, load, type, FCOMM_IFELSE(fcomm, NULL));
  else
    mpi_fmm_sort_back_merge_body(mem0, mem1, mem_sizes, ntotal, nin, addr, q, xyz, pot, grad, type, FCOMM_IFELSE(fcomm, NULL));
}

void mpi_fmm_sort_back_mem_(
 void *mem0, void *mem1, pint_t *mem_sizes,
 pint_t *ntotal, pint_t *nin, pint_t *nout, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, back_(sldata4_t) *load, pint_t *type
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  mpi_fmm_sort_back_mem(mem0, mem1, mem_sizes, ntotal, nin, nout, addr, q, xyz, pot, grad, load, type
#ifdef WITH_FCOMM
   , fcomm
#endif
  );
}


void mpi_fmm_sort_back(
 pint_t *ntotal, pint_t *nin, pint_t *nout, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, back_(sldata4_t) *load, pint_t *type
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  if (*nout >= 0 || mpi_fmm_sort_back_part)
    mpi_fmm_sort_back_part_body(NULL, NULL, NULL, ntotal, nin, nout, addr, q, xyz, pot, grad, load, type, FCOMM_IFELSE(fcomm, NULL));
  else
    mpi_fmm_sort_back_merge_body(NULL, NULL, NULL, ntotal, nin, addr, q, xyz, pot, grad, type, FCOMM_IFELSE(fcomm, NULL));
}

void mpi_fmm_sort_back_(
 pint_t *ntotal, pint_t *nin, pint_t *nout, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, back_(sldata4_t) *load, pint_t *type
#ifdef WITH_FCOMM
 , FINT8_TYPE_C *fcomm
#endif
)
{
  mpi_fmm_sort_back(ntotal, nin, nout, addr, q, xyz, pot, grad, load, type
#ifdef WITH_FCOMM
   , fcomm
#endif
  );
}

#endif /* NO_SL_BACK */
