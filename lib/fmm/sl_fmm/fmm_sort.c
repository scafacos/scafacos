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

#include "config_fmm_sort.h"
#include "rename_fmm_sort.h"

/*#define DO_VERBOSE
#define DO_VALIDATE
#define DO_CHECKSUM
#define DO_TIMING*/

#include "common.h"

#include "fmm_sort.h"


/*#define FMM_SORT_RADIX_1BIT*/

/*#define SORT_BACK_ONE2NCHARGES*/

/*#define SL_OUTPUT_TO_FILE*/

int SL_FMM_CONFIG_VAR(fmm_front_aX);
INTEGER_C SL_FMM_CONFIG_VAR(fmm_front_key_mask);


#ifdef FMM_SORT_RADIX_1BIT
# define fmm_sort_radix(_prefix_, _s_, _sx_, _h_, _l_, _w_)  _prefix_##sort_radix_1bit(_s_, _sx_, _h_, _l_)
#else
# define fmm_sort_radix(_prefix_, _s_, _sx_, _h_, _l_, _w_)  _prefix_##sort_radix(_s_, _sx_, _h_, _l_, _w_)
#endif


static INTEGER_C log2_floor(INTEGER_C v)
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


#ifndef NO_SL_FRONT

static void fmm_sort_front_body(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *depth, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type)
{
  typedef front_(slint_t) front_slint_t;

  front_slint_t front_type, i, highest;

#ifndef NOT_sl_front_xqsa0
  front_xqsa0_elements_t s0, *sx0, smem0;
#endif
#ifndef NOT_sl_front_xqsaI
  front_xqsaI_elements_t s1, *sx1, smem1;
#endif
#ifndef NOT_sl_front_xqsaX
  front_xqsaX_elements_t s2, *sx2, smem2;
#endif
#ifndef NOT_sl_front_xq_a0
  front_xq_a0_elements_t s3, *sx3, smem3;
#endif
#ifndef NOT_sl_front_xq_aI
  front_xq_aI_elements_t s4, *sx4, smem4;
#endif
#ifndef NOT_sl_front_xq_aX
  front_xq_aX_elements_t s5, *sx5, smem5;
#endif

  void *mems[] = { mem0, mem1 };
  front_slint_t mem_sizes_sl[2];

#ifdef VALIDATE
  front_slint_t o;
# ifdef CHECKSUM
  unsigned int crc32_in, crc32_out;
# endif
#endif

#ifdef SL_OUTPUT_TO_FILE
  char output_file_str[32];
  FILE *sl_debug_fstream;
#endif


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
    INFO_CMD(printf(INFO_PRINT_PREFIX "front: addr_desc = %p\n", addr_desc););
    front_type += 1;

  } else
  {
    INFO_CMD(printf(INFO_PRINT_PREFIX "front: addr_desc = %" PARAM_INTEGER_FMT ", %" PARAM_INTEGER_FMT "\n", addr_desc[0], addr_desc[1]););
    if (addr_desc[0] == 0) front_type += 0;
    else if (addr_desc[0] == sizeof(INTEGER_C) && addr_desc[1] <= 0) front_type += 1;
    else front_type += 2;
  }

  INFO_CMD(printf(INFO_PRINT_PREFIX "front: type = %" front_slint_fmt "\n", front_type););

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

#ifdef SL_OUTPUT_TO_FILE
  sprintf(output_file_str, "sort_front.debug");
  front_xqsa0(front_xqsa0_sl_debug_fstream = )
  front_xqsaI(front_xqsaI_sl_debug_fstream = )
  front_xqsaX(front_xqsaX_sl_debug_fstream = )
  front_xq_a0(front_xq_a0_sl_debug_fstream = )
  front_xq_aI(front_xq_aI_sl_debug_fstream = )
  front_xq_aX(front_xq_aX_sl_debug_fstream = ) sl_debug_fstream = fopen(output_file_str, "w");
#endif

  SL_FMM_CONFIG_VAR(fmm_front_aX) = addr_desc[0];

  if (mem_sizes && ((mems[0] && mem_sizes[0] > 0) || (mems[1] && mem_sizes[1])))
  {
    mem_sizes_sl[0] = mem_sizes[0];
    mem_sizes_sl[1] = mem_sizes[1];
    switch (front_type)
    {
#ifndef NOT_sl_front_xqsa0
      case 0:
        front_xqsa0_elements_alloc_from_blocks(&smem0, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
        sx0 = &smem0;
        break;
#endif
#ifndef NOT_sl_front_xqsaI
      case 1:
        front_xqsaI_elements_alloc_from_blocks(&smem1, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
        sx1 = &smem1;
        break;
#endif
#ifndef NOT_sl_front_xqsaX
      case 2:
        front_xqsaX_elements_alloc_from_blocks(&smem2, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
        sx2 = &smem2;
        break;
#endif
#ifndef NOT_sl_front_xq_a0
      case 3:
        front_xq_a0_elements_alloc_from_blocks(&smem3, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
        sx3 = &smem3;
        break;
#endif
#ifndef NOT_sl_front_xq_aI
      case 4:
        front_xq_aI_elements_alloc_from_blocks(&smem4, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
        sx4 = &smem4;
        break;
#endif
#ifndef NOT_sl_front_xq_aX
      case 5:
        front_xq_aX_elements_alloc_from_blocks(&smem5, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
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
    printf(INFO_PRINT_PREFIX "front: type = %" front_slint_fmt "\n", front_type);
    printf(INFO_PRINT_PREFIX "front: sizeof(front_slint_t) = %d byte - #elements: %" PARAM_INTEGER_FMT " - mem0: %p / %" PARAM_INTEGER_FMT " bytes, mem1: %p / %" PARAM_INTEGER_FMT " bytes\n", (int) sizeof(front_(slint_t)), *n, mem0, (mem_sizes)?mem_sizes[0]:0, mem1, (mem_sizes)?mem_sizes[1]:0);
    printf(INFO_PRINT_PREFIX "front: ibox:      sizeof(front_slkey_t)   = %d byte\n", (int) sizeof(front_(slkey_t)));
    printf(INFO_PRINT_PREFIX "front:  xyz: %d * sizeof(front_sldata0_t) = %d * %d byte\n", (int) front_(sl_data0_size_c), (int) front_(sl_data0_size_c), (int) sizeof(front_(sldata0_t)));
    printf(INFO_PRINT_PREFIX "front:    q: %d * sizeof(front_sldata1_t) = %d * %d byte\n", (int) front_(sl_data1_size_c), (int) front_(sl_data1_size_c), (int) sizeof(front_(sldata1_t)));
    printf(INFO_PRINT_PREFIX "front:  scr: %d * sizeof(front_sldata2_t) = %d * %d byte\n", (int) front_(sl_data2_size_c), (int) front_(sl_data2_size_c), (int) sizeof(front_(sldata2_t)));
/*    switch (front_type)
    {
      case 0:
      case 3:
        printf(INFO_PRINT_PREFIX "front: addr: %d * sizeof(front_xqsa0_sldata3_t) = %d * %d byte\n", (int) front_xqsa0_sl_data3_size_c, (int) front_xqsa0_sl_data3_size_c, (int) sizeof(front_xqsa0_sldata3_t));
        break;
      case 1:
      case 4:
        printf(INFO_PRINT_PREFIX "front: addr: %d * sizeof(front_xqsaI_sldata3_t) = %d * %d byte\n", (int) front_xqsaI_sl_data3_size_c, (int) front_xqsaI_sl_data3_size_c, (int) sizeof(front_xqsaI_sldata3_t));
        break;
      case 2:
      case 5:
        printf(INFO_PRINT_PREFIX "front: addr: %d * sizeof(front_xqsaX_sldata3_t) = %d * %d byte\n", (int) front_xqsaX_sl_data3_size_c, (int) front_xqsaX_sl_data3_size_c, (int) sizeof(front_xqsaX_sldata3_t));
        break;
    }*/
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
    case 0: crc32_in = front_xqsa0_elements_crc32(&s0, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1: crc32_in = front_xqsaI_elements_crc32(&s1, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2: crc32_in = front_xqsaX_elements_crc32(&s2, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3: crc32_in = front_xq_a0_elements_crc32(&s3, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4: crc32_in = front_xq_aI_elements_crc32(&s4, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5: crc32_in = front_xq_aX_elements_crc32(&s5, 1, 1, 1); break;
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
    if (depth == NULL) printf(INFO_PRINT_PREFIX "front: starting local radix-sort (lowest 3 bit)\n");
    else printf(INFO_PRINT_PREFIX "front: starting local radix-sort (depth = %" PARAM_INTEGER_FMT " -> sorting bits [0..%" front_slint_fmt "])\n", *depth, highest);
  );

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

  INFO_CMD(printf(INFO_PRINT_PREFIX "front: local radix-sort done\n"););

  if (*subx != 0) for (i = 0; i < *n; i++) ibox[i] += *subx;

#ifdef VALIDATE
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0: o = front_xqsa0_elements_validate_order(&s0, 1); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1: o = front_xqsaI_elements_validate_order(&s1, 1); break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2: o = front_xqsaX_elements_validate_order(&s2, 1); break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3: o = front_xq_a0_elements_validate_order(&s3, 1); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4: o = front_xq_aI_elements_validate_order(&s4, 1); break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5: o = front_xq_aX_elements_validate_order(&s5, 1); break;
#endif
    default: o = 1;
  }
  INFO_CMD(printf(INFO_PRINT_PREFIX "front: sorting_order: %s\n", (!o)?"success":"FAILED"););

# ifdef CHECKSUM
  switch (front_type)
  {
#ifndef NOT_sl_front_xqsa0
    case 0: crc32_out = front_xqsa0_elements_crc32(&s0, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xqsaI
    case 1: crc32_out = front_xqsaI_elements_crc32(&s1, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xqsaX
    case 2: crc32_out = front_xqsaX_elements_crc32(&s2, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xq_a0
    case 3: crc32_out = front_xq_a0_elements_crc32(&s3, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xq_aI
    case 4: crc32_out = front_xq_aI_elements_crc32(&s4, 1, 1, 1); break;
#endif
#ifndef NOT_sl_front_xq_aX
    case 5: crc32_out = front_xq_aX_elements_crc32(&s5, 1, 1, 1); break;
#endif
  }
  INFO_CMD(printf(INFO_PRINT_PREFIX "front: checksum: %s (%X vs. %X)\n", (crc32_in == crc32_out)?"success":"FAILED", crc32_in, crc32_out););
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
}


void fmm_sort_front_mem(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *depth, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type)
{
  fmm_sort_front_body(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type);
}

void fmm_sort_front_mem_(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *depth, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type)
{
  fmm_sort_front_mem(mem0, mem1, mem_sizes, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type);
}


void fmm_sort_front(pint_t *depth, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type)
{
  fmm_sort_front_body(NULL, NULL, NULL, depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type);
}

void fmm_sort_front_(pint_t *depth, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type)
{
  fmm_sort_front(depth, subx, n, ibox, xyz, q, addr_desc, addr, scr, type);
}


void fmm_sort_front_3bit_mem(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type)
{
  fmm_sort_front_body(mem0, mem1, mem_sizes, NULL, subx, n, ibox, xyz, q, addr_desc, addr, scr, type);
}

void fmm_sort_front_3bit_mem_(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type)
{
  fmm_sort_front_3bit_mem(mem0, mem1, mem_sizes, subx, n, ibox, xyz, q, addr_desc, addr, scr, type);
}


void fmm_sort_front_3bit(pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type)
{
  fmm_sort_front_body(NULL, NULL, NULL, NULL, subx, n, ibox, xyz, q, addr_desc, addr, scr, type);
}

void fmm_sort_front_3bit_(pint_t *subx, pint_t *n, front_(slkey_t) *ibox, front_(sldata0_t) *xyz, front_(sldata1_t) *q, pint_t *addr_desc, void *addr, front_(sldata2_t) *scr, pint_t *type)
{
  fmm_sort_front_3bit(subx, n, ibox, xyz, q, addr_desc, addr, scr, type);
}

#endif /* NO_SL_FRONT */


#ifndef NO_SL_BACK

void fmm_sort_back_mem(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *n, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, pint_t *type)
{
  typedef back_(slint_t) back_slint_t;

  back_slint_t back_type, highest;

#ifndef NOT_sl_back_qxpg
  back_qxpg_elements_t s0, *sx0, smem0;
#endif
#ifndef NOT_sl_back_qx_g
  back_qx_g_elements_t s1, *sx1, smem1;
#endif
#ifndef NOT_sl_back_q_pg
  back_q_pg_elements_t s2, *sx2, smem2;
#endif
#ifndef NOT_sl_back_q__g
  back_q__g_elements_t s3, *sx3, smem3;
#endif

  void *mems[] = { mem0, mem1 };
  back_slint_t mem_sizes_sl[2];

#ifdef VALIDATE
  back_slint_t o;
# ifdef CHECKSUM
  unsigned int crc32_in, crc32_out;
# endif
#endif

#ifdef SL_OUTPUT_TO_FILE
  char output_file_str[32];
  FILE *sl_debug_fstream;
#endif


  if (type && *type >= 0) back_type = *type;
  else
  {
    if (xyz != NULL && pot != NULL) back_type = 0;
    else if (xyz != NULL && pot == NULL) back_type = 1;
    else if (xyz == NULL && pot != NULL) back_type = 2;
    else if (xyz == NULL && pot == NULL) back_type = 3;
    else back_type = -1;
  }

  INFO_CMD(printf(INFO_PRINT_PREFIX "back: type = %" back_slint_fmt "\n", back_type););

#ifdef SL_OUTPUT_TO_FILE
  sprintf(output_file_str, "sort_back.debug");
  back_qxpg(back_qxpg_sl_debug_fstream = )
  back_qx_g(back_qx_g_sl_debug_fstream = )
  back_q_pg(back_q_pg_sl_debug_fstream = )
  back_q__g(back_q__g_sl_debug_fstream = ) sl_debug_fstream = fopen(output_file_str, "w");
#endif

  if (mem_sizes && ((mems[0] && mem_sizes[0] > 0) || (mems[1] && mem_sizes[1])))
  {
    mem_sizes_sl[0] = mem_sizes[0];
    mem_sizes_sl[1] = mem_sizes[1];
    switch (back_type)
    {
#ifndef NOT_sl_back_qxpg
      case 0:
        back_qxpg_elements_alloc_from_blocks(&smem0, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
        sx0 = &smem0;
        break;
#endif
#ifndef NOT_sl_back_qx_g
      case 1:
        back_qx_g_elements_alloc_from_blocks(&smem1, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
        sx1 = &smem1;
        break;
#endif
#ifndef NOT_sl_back_q_pg
      case 2:
        back_q_pg_elements_alloc_from_blocks(&smem2, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
        sx2 = &smem2;
        break;
#endif
#ifndef NOT_sl_back_q__g
      case 3:
        back_q__g_elements_alloc_from_blocks(&smem3, 2, mems, mem_sizes_sl, 8, *n, SLCM_KEYS|SLCM_DATA);
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
    printf(INFO_PRINT_PREFIX "back: sizeof(back_slint_t) = %d bytes - #elements: %" PARAM_INTEGER_FMT " - mem0: %p / %" PARAM_INTEGER_FMT " bytes, mem1: %p / %" PARAM_INTEGER_FMT " bytes\n", (int) sizeof(back_slint_t), *n, mem0, (mem_sizes)?mem_sizes[0]:0, mem1, (mem_sizes)?mem_sizes[1]:0);
    printf(INFO_PRINT_PREFIX "back: addr:      sizeof(back_slkey_t)   = %d byte\n", (int) sizeof(back_(slkey_t)));
    printf(INFO_PRINT_PREFIX "back:    q: %d * sizeof(back_sldata0_t) = %d * %d byte\n", (int) back_(sl_data0_size_c), (int) back_(sl_data0_size_c), (int) sizeof(back_(sldata0_t)));
    printf(INFO_PRINT_PREFIX "back:  xyz: %d * sizeof(back_sldata1_t) = %d * %d byte\n", (int) back_(sl_data1_size_c), (int) back_(sl_data1_size_c), (int) sizeof(back_(sldata1_t)));
    printf(INFO_PRINT_PREFIX "back:  pot: %d * sizeof(back_sldata2_t) = %d * %d byte\n", (int) back_(sl_data2_size_c), (int) back_(sl_data2_size_c), (int) sizeof(back_(sldata2_t)));
    printf(INFO_PRINT_PREFIX "back: grad: %d * sizeof(back_sldata3_t) = %d * %d byte\n", (int) back_(sl_data3_size_c), (int) back_(sl_data3_size_c), (int) sizeof(back_(sldata3_t)));
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
    case 0: crc32_in = back_qxpg_elements_crc32(&s0, 1, 1, 1); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: crc32_in = back_qx_g_elements_crc32(&s1, 1, 1, 1); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: crc32_in = back_q_pg_elements_crc32(&s2, 1, 1, 1); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: crc32_in = back_q__g_elements_crc32(&s3, 1, 1, 1); break;
#endif
  }
#endif

#ifdef SORT_BACK_ONE2NCHARGES
  highest = log2_floor(*n);
#else
  highest = log2_floor(*n - 1);
#endif

  INFO_CMD(printf(INFO_PRINT_PREFIX "back: starting local radix-sort (ncharges = %" PARAM_INTEGER_FMT " -> sorting bits [0..%" back_slint_fmt "])\n", *n, highest););

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

  INFO_CMD(printf(INFO_PRINT_PREFIX "back: local radix-sort done\n"););

#ifdef VALIDATE
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: o = back_qxpg_elements_validate_order(&s0, 1); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: o = back_qx_g_elements_validate_order(&s1, 1); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: o = back_q_pg_elements_validate_order(&s2, 1); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: o = back_q__g_elements_validate_order(&s3, 1); break;
#endif
    default: o = 1;
  }
  INFO_CMD(printf(INFO_PRINT_PREFIX "back: sorting_order: %s\n", (!o)?"success":"FAILED"););

# ifdef CHECKSUM
  switch (back_type)
  {
#ifndef NOT_sl_back_qxpg
    case 0: crc32_out = back_qxpg_elements_crc32(&s0, 1, 1, 1); break;
#endif
#ifndef NOT_sl_back_qx_g
    case 1: crc32_out = back_qx_g_elements_crc32(&s1, 1, 1, 1); break;
#endif
#ifndef NOT_sl_back_q_pg
    case 2: crc32_out = back_q_pg_elements_crc32(&s2, 1, 1, 1); break;
#endif
#ifndef NOT_sl_back_q__g
    case 3: crc32_out = back_q__g_elements_crc32(&s3, 1, 1, 1); break;
#endif
  }
  INFO_CMD(printf(INFO_PRINT_PREFIX "back: checksum: %s (%X vs. %X)\n", (crc32_in == crc32_out)?"success":"FAILED", crc32_in, crc32_out););
# endif
#endif

#ifdef SL_OUTPUT_TO_FILE
  fclose(sl_debug_fstream);
  back_qxpg(back_qxpg_sl_debug_fstream = )
  back_qx_g(back_qx_g_sl_debug_fstream = )
  back_q_pg(back_q_pg_sl_debug_fstream = )
  back_q__g(back_q__g_sl_debug_fstream = ) sl_debug_fstream = NULL;
#endif
}

void fmm_sort_back_mem_(void *mem0, void *mem1, pint_t *mem_sizes, pint_t *n, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, pint_t *type)
{
  fmm_sort_back_mem(mem0, mem1, mem_sizes, n, addr, q, xyz, pot, grad, type);
}


void fmm_sort_back(pint_t *n, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, pint_t *type)
{
  fmm_sort_back_mem(NULL, NULL, NULL, n, addr, q, xyz, pot, grad, type);
}

void fmm_sort_back_(pint_t *n, back_(slkey_t) *addr, back_(sldata0_t) *q, back_(sldata1_t) *xyz, back_(sldata2_t) *pot, back_(sldata3_t) *grad, pint_t *type)
{
  fmm_sort_back(n, addr, q, xyz, pot, grad, type);
}

#endif /* NO_SL_BACK */


#if 0

void slcheck_fortran2c_types(double *f2c_fint, double *f2c_key, double *f2c_real)
{
  int error = 0;

  if (*f2c_fint != sizeof(blue_slint_t)) { fprintf(stderr, "WARNING: fortran integer = %d vs. slint_t = %d\nTry to change the FINT_DEFAULT flag or adjust file 'fortran2c_types.h'!\n", (int) *f2c_fint, (int) sizeof(blue_slint_t)); ++error; }
  if (*f2c_key != sizeof(blue_slkey_t)) { fprintf(stderr, "WARNING: fortran keys = %d vs. slkey_t = %d\nTry to change the INTEGER32 flag or adjust file 'fortran2c_types.h'!\n", (int) *f2c_key, (int) sizeof(blue_slkey_t)); ++error; }
  if (*f2c_real != sizeof(blue_sldata0_t)) { fprintf(stderr, "WARNING: fortran real = %d vs. sldata0_t = %d\nTry to change the REAL32 flag or adjust file 'fortran2c_types.h'!\n", (int) *f2c_real, (int) sizeof(blue_sldata0_t)); ++error; }

  if (error) fprintf(stderr, "WARNING: There seems to be a problem between Fortran and C data types.\n");
}

void slcheck_fortran2c_types_(double *f2c_fint, double *f2c_key, double *f2c_real)
{
  slcheck_fortran2c_types(f2c_fint, f2c_key, f2c_real);
}

#endif
