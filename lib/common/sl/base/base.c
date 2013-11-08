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


/* sl_macro B_TRACE_IF */


#include "sl_common.h"


#ifndef B_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define B_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define B_TRACE_IF  (SL_PROC_RANK == -1)
# endif
#endif


slint_t binning_create(local_bins_t *lb, slint_t max_nbins, slint_t max_nbinnings, elements_t *s, slint_t nelements, slint_t docounts, slint_t doweights, binning_t *bm) /* sl_proto, sl_func binning_create */
{
  slint_t j;

  bm->docounts = docounts;
#ifdef elem_weight
  bm->doweights = doweights;
#endif
  lb->bm = bm;

  lb->nbins = 0;
  lb->max_nbins = max_nbins;

  lb->nbinnings = 0;
  lb->max_nbinnings = (max_nbinnings > 0)?z_min(max_nbinnings,max_nbins):max_nbins;

  lb->nelements = nelements;

  lb->docounts = docounts;
#ifdef elem_weight
  lb->doweights = doweights;
#endif

  lb->bins0 = z_alloc(lb->max_nbins * lb->nelements, sizeof(bin_t));
  lb->bins1 = z_alloc(lb->max_nbins * lb->nelements, sizeof(bin_t));

  lb->bcws = z_alloc(lb->max_nbins, sizeof(slint_t));

#if defined(elem_weight) && defined(sl_weight_intequiv)
  lb->cw_factor = (lb->docounts != 0) + (lb->doweights != 0);
  lb->w_index = (lb->docounts != 0);
  lb->cws = z_alloc(lb->max_nbinnings * lb->cw_factor * lb->bm->max_nbins, sizeof(slweight_t));
  lb->bin_cw_factor = 1 + (lb->doweights != 0);
  lb->bin_cws = z_alloc(lb->max_nbinnings * lb->bin_cw_factor * lb->nelements * bm->max_nbins, sizeof(slweight_t));
  lb->prefix_cws = z_alloc(lb->max_nbinnings * lb->cw_factor, sizeof(slweight_t));
#else
  lb->cs = z_alloc(lb->max_nbinnings * 1 * lb->bm->max_nbins, sizeof(slint_t));
  lb->bin_cs = z_alloc(lb->max_nbinnings * 1 * lb->nelements * bm->max_nbins, sizeof(slint_t));
  lb->prefix_cs = z_alloc(lb->max_nbinnings * 1, sizeof(slint_t));
# ifdef elem_weight
  if (lb->doweights)
  {
    lb->ws = z_alloc(lb->max_nbinnings * 1 * lb->bm->max_nbins, sizeof(slweight_t));
    lb->bin_ws = z_alloc(lb->max_nbinnings * 1 * lb->nelements * bm->max_nbins, sizeof(slweight_t));
    lb->prefix_ws = z_alloc(lb->max_nbinnings * 1, sizeof(slweight_t));

  } else lb->ws = lb->bin_ws = lb->prefix_ws = NULL;
# endif
#endif

  lb->last_exec_b = -1;

  lb->bins = lb->bins0;
  lb->bins_new = lb->bins1;
  
  if (lb->max_nbins > 0)
  {
    lb->nbins = 1;
    for (j = 0; j < lb->nelements; ++j) elem_assign(&s[j], &lb->bins[0 * lb->nelements + j].s);
  }

  return 0;
}


slint_t binning_destroy(local_bins_t *lb) /* sl_proto, sl_func binning_destroy */
{
  z_free(lb->bins0);
  z_free(lb->bins1);
  
  z_free(lb->bcws);

#if defined(elem_weight) && defined(sl_weight_intequiv)
  z_free(lb->cws);
  z_free(lb->bin_cws);
  z_free(lb->prefix_cws);
#else
  z_free(lb->cs);
  z_free(lb->bin_cs);
  z_free(lb->prefix_cs);
# ifdef elem_weight
  if (lb->doweights)
  {
    z_free(lb->ws);
    z_free(lb->bin_ws);
    z_free(lb->prefix_ws);
  }
# endif
#endif

  return 0;
}


slint_t binning_pre(local_bins_t *lb) /* sl_proto, sl_func binning_pre */
{
  lb->bm->pre(lb->bm);

  lb->nbins_new = 0;
  lb->last_new_b = lb->last_new_k = -1;

  return 0;
}


slint_t binning_exec_reset(local_bins_t *lb, slint_t do_bins, slint_t do_prefixes) /* sl_proto, sl_func binning_exec_reset */
{
  lb->nbinnings = 0;

  return 0;
}


slint_t binning_exec(local_bins_t *lb, slint_t b, slint_t do_bins, slint_t do_prefixes) /* sl_proto, sl_func binning_exec */
{
  slint_t j;
  slkey_pure_t k;

  slcount_t *counts, *bin_counts, *prefix_count;
#ifdef elem_weight
  slweight_t *weights, *bin_weights, *prefix_weight;
#endif


  Z_TRACE_IF(B_TRACE_IF, "do_bins: %" slint_fmt ", do_prefixes: %" slint_fmt ", b: %" slint_fmt ", last_exec_b: %" slint_fmt ", nbinnings: %" slint_fmt " of %" slint_fmt,
    do_bins, do_prefixes, b, lb->last_exec_b, lb->nbinnings, lb->max_nbinnings);

  if (lb->last_exec_b == b) return 0;

  if (lb->nbinnings >= lb->max_nbinnings) return -1;
  
  lb->bcws[b] = lb->nbinnings;
  
  ++lb->nbinnings;
  
  lb->last_exec_b = b;

  if (do_bins)
  {
    counts = lb_counts(lb, b, 0);
    bin_counts = lb_bin_counts(lb, b, 0, 0);
#ifdef elem_weight
    weights = lb_weights(lb, b, 0);
    bin_weights = lb_bin_weights(lb, b, 0, 0);
#endif

#ifdef elem_weight
    if (lb->doweights)
    {
      for (k = 0; k < lb->bm->nbins; ++k) counts[k] = weights[k] = 0.0;
  
      for (j = 0; j < lb->nelements; ++j)
      for (k = 0; k < lb->bm->nbins; ++k) bin_counts[j * lb->bm->nbins + k] = bin_weights[j * lb->bm->nbins + k] = 0.0;

      /* for every list of elements */
      for (j = 0; j < lb->nelements; ++j)
      {
        Z_TRACE_IF(B_TRACE_IF, "bin %" slint_fmt ",%" slint_fmt ": size = %" slint_fmt, b, j, lb->bins[b * lb->nelements + j].s.size);
    
        lb->bm->exec(lb->bm, &lb->bins[b * lb->nelements + j], bin_counts, elem_weight_ifelse(bin_weights, NULL));

        lb->bins[b * lb->nelements + j].weight = 0;

        for (k = 0; k < lb->bm->nbins; ++k)
        {
          if (lb->docounts) counts[k] += bin_counts[k];

          weights[k] += bin_weights[k];
          lb->bins[b * lb->nelements + j].weight += bin_weights[k];
        }

        Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slcount_fmt, bin_counts[k], "%" slint_fmt ",%" slint_fmt ": bin_counts =", b, j);
        Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slweight_fmt, bin_weights[k], "%" slint_fmt ",%" slint_fmt ": bin_weights =", b, j);

        bin_counts += lb->bm->nbins;
        bin_weights += lb->bm->nbins;
      }

    } else
#endif
    {
      for (k = 0; k < lb->bm->nbins; ++k) counts[k] = 0;
  
      for (j = 0; j < lb->nelements; ++j)
      for (k = 0; k < lb->bm->nbins; ++k) bin_counts[j * lb->bm->nbins + k] = 0.0;

      /* for every list of elements */
      for (j = 0; j < lb->nelements; ++j)
      {
        Z_TRACE_IF(B_TRACE_IF, "bin %" slint_fmt ",%" slint_fmt ": size = %" slint_fmt, b, j, lb->bins[b * lb->nelements + j].s.size);
    
        lb->bm->exec(lb->bm, &lb->bins[b * lb->nelements + j], bin_counts, elem_weight_ifelse(bin_weights, NULL));

        for (k = 0; k < lb->bm->nbins; ++k)
        {
          counts[k] += bin_counts[k];
        }
    
        Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slcount_fmt, bin_counts[k], "%" slint_fmt ",%" slint_fmt ": bin_counts =", b, j);

        bin_counts += lb->bm->nbins;
      }
    }

    if (lb->docounts)
      Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slcount_fmt, counts[k], "%" slint_fmt ": counts =", b);
#ifdef elem_weight
    if (lb->doweights)
      Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slweight_fmt, weights[k], "%" slint_fmt ": weights =", b);
#endif
  }
  
  if (do_prefixes)
  {
    prefix_count = lb_prefix_count(lb, b);
#ifdef elem_weight
    prefix_weight = lb_prefix_weight(lb, b);
#endif

#ifdef elem_weight
    if (lb->doweights)
    {
      *prefix_count = 0;
      *prefix_weight = 0;

      for (j = 0; j < lb->nelements; ++j)
      {
        *prefix_count += lb_bin_count(lb, b, j);
        *prefix_weight += lb_bin_weight(lb, b, j);
      }
    } else
#endif
    {
      *prefix_count = 0;

      for (j = 0; j < lb->nelements; ++j)
      {
        *prefix_count += lb_bin_count(lb, b, j);
      }
    }

    if (lb->docounts)
      Z_TRACE_IF(B_TRACE_IF, "%" slint_fmt ": prefix_count = %" slcount_fmt, b, *prefix_count);
#ifdef elem_weight
    if (lb->doweights)
      Z_TRACE_IF(B_TRACE_IF, "%" slint_fmt ": prefix_weight = %" slweight_fmt, b, *prefix_weight);
#endif
  }

  return 0;
}


slint_t binning_refine(local_bins_t *lb, slint_t b, slint_t k, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_refine */
{
  slint_t j;
  bin_t *new_bin = NULL;

  if (lb->last_new_b != b || lb->last_new_k != k)
  {
    /* update last_new_... */
    lb->last_new_b = b;
    lb->last_new_k = k;
    
    new_bin = &lb->bins_new[lb->nbins_new * lb->nelements];
    
    ++lb->nbins_new;

  } else Z_TRACE_IF(B_TRACE_IF, "no new bin, with b = %" slint_fmt " and k = %" slint_fmt, b, k);

  /* create new bin */
  for (j = 0; j < lb->nelements; ++j)
  {
    lb->bm->refine(lb->bm, &lb->bins[b * lb->nelements + j], k, lb_bin_counts(lb, b, j, 0), lb_bin_weights(lb, b, j, 0), sp, s * lb->nelements + j, new_bin);

    if (new_bin)
    {
      Z_TRACE_IF(B_TRACE_IF, "new bin %td count: %" slint_fmt, new_bin - lb->bins_new, new_bin->s.size);
#ifdef elem_weight
      if (lb->doweights)
        Z_TRACE_IF(B_TRACE_IF, "new bin %td count: %" slweight_fmt, new_bin - lb->bins_new, new_bin->weight);
#endif
      ++new_bin;

    }
  }

/*  Z_TRACE_IF(B_TRACE_IF, "b: %" slint_fmt ", k: %" slint_fmt ", returning %" slint_fmt, b, k, lb->nbins_new - 1);*/

  return lb->nbins_new - 1;
}


slint_t binning_hit(local_bins_t *lb, slint_t b, slint_t k, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_hit */
{
  slint_t j;

  for (j = 0; j < lb->nelements; ++j)
  {
    lb->bm->hit(lb->bm, &lb->bins[b * lb->nelements + j], k, lb_bin_counts(lb, b, j, 0), sp, s * lb->nelements + j);
  }

  return 0;
}


slint_t binning_finalize(local_bins_t *lb, slint_t b, slint_t dc, slweight_t dw, slint_t lc_min, slint_t lc_max, slcount_t *lcs, slweight_t *lws, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_finalize */
{
  slint_t j, dc_left;
#ifdef elem_weight
  slweight_t dw_left;
#else
# define dw_left  0
#endif


  Z_TRACE_IF(B_TRACE_IF, "b: %" slint_fmt ", dc: %" slint_fmt ", dw: %" slweight_fmt ", lc_min: %" slint_fmt ", lc_max: %" slint_fmt, b, dc, dw, lc_min, lc_max);

  dc_left = dc;
  *lcs = 0.0;
#ifdef elem_weight
  dw_left = dw;
  *lws = 0.0;
#endif

  for (j = 0; j < lb->nelements; ++j)
  {
#ifdef elem_weight
    if (lb->doweights)
      dw_left = dw - *lws;
    else
#endif
      dc_left = dc - *lcs;

    if (lb->bm->finalize(lb->bm, &lb->bins[b * lb->nelements + j], dc_left, dw_left, lc_min - *lcs, lc_max - *lcs, lcs, lws, sp, s * lb->nelements + j)) break;
  }

  return 0;
}


slint_t binning_post(local_bins_t *lb) /* sl_proto, sl_func binning_post */
{
  lb->bm->post(lb->bm);
  
  lb->last_exec_b = -1;

  lb->nbins = lb->nbins_new;

  if (lb->bins == lb->bins0)
  {
    lb->bins = lb->bins1;
    lb->bins_new = lb->bins0;

  } else
  {
    lb->bins = lb->bins0;
    lb->bins_new = lb->bins1;
  }

  return 0;
}



/* sl_macro BR_TRACE_IF */


#include "sl_common.h"


#ifndef BR_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define BR_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define BR_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t binning_radix_create(binning_t *bm, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted) /* sl_proto, sl_func binning_radix_create */
{
  bm->pre = binning_radix_pre;
  bm->exec = binning_radix_exec;
  bm->refine = binning_radix_refine;
  bm->hit = binning_radix_hit;
  bm->finalize = binning_radix_finalize;
  bm->post = binning_radix_post;

  bm->sorted = sorted;

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;
  if (rwidth < 0) rwidth = sort_radix_width_default;

  bm->nbins = 0;
  bm->max_nbins = z_powof2_typed(rwidth, slkey_pure_t);

  bm->bd.radix.rhigh = rhigh;
  bm->bd.radix.rlow = rlow;
  bm->bd.radix.rwidth = rwidth;

  elements_alloc(&bm->bd.radix.sx, 1, SLCM_ALL);

  return 0;
}


slint_t binning_radix_destroy(binning_t *bm) /* sl_proto, sl_func binning_radix_destroy */
{
  elements_free(&bm->bd.radix.sx);

  return 0;
}


slint_t binning_radix_pre(binning_t *bm) /* sl_proto, sl_func binning_radix_pre */
{
  bm->bd.radix.rcurrent = z_min(bm->bd.radix.rwidth, bm->bd.radix.rhigh - bm->bd.radix.rlow + 1);
  bm->bd.radix.rhigh -= (bm->bd.radix.rcurrent > 0)?bm->bd.radix.rcurrent - 1:bm->bd.radix.rhigh;

  bm->nbins = (bm->bd.radix.rcurrent > 0)?z_powof2(bm->bd.radix.rcurrent):1;
  bm->bd.radix.bit_mask = bm->nbins - 1;
  
  return 0;
}


slint_t binning_radix_exec(binning_t *bm, bin_t *bin, slcount_t *counts, slweight_t *weights) /* sl_proto, sl_func binning_radix_exec */
{
  elements_t xi, end;
  slkey_pure_t k;
  slint_t i, *c, special;


  elem_assign_at(&bin->s, bin->s.size, &end);

  if (bm->nbins > 1)
  {
#ifdef elem_weight
    if (bm->doweights)
    {
      /* counts and weights in every class */
      for (elem_assign(&bin->s, &xi); xi.keys < end.keys; elem_inc(&xi))
      {
        k = key_radix_key2class(key_purify(*xi.keys), bm->bd.radix.rhigh, bm->bd.radix.bit_mask);
        counts[k] += 1;
        weights[k] += elem_weight(&xi, 0);
      }

    } else
#endif
    {
      /* counts in every class */
      if (bm->sorted & SL_SORTED_IN)
      {
        special = ((bm->bd.radix.bit_mask << bm->bd.radix.rhigh) < 0);
      
        Z_TRACE_IF(BR_TRACE_IF, "bitmask = %" key_pure_type_fmt " -> %s", (bm->bd.radix.bit_mask << bm->bd.radix.rhigh), special?"special signedness handling":"normal");

/*        if (BR_TRACE_IF) elements_print_keys(&bin->s);*/

        if (special)
        {
          special = sl_search_binary_sign_switch(&bin->s);
        
          Z_TRACE_IF(BR_TRACE_IF, "sign switch @ %" slint_fmt, special);
          
          /* make bit mask 1 bit smaller (and erase the sign bit) */
          bm->bd.radix.bit_mask = (bm->bd.radix.bit_mask >> 1) & ~(~((slkey_pure_t) 0) << (sizeof(slkey_pure_t)*8-1));

          elem_assign(&bin->s, &xi);
          xi.size = special;
          
          for (i = 0; i < bm->bd.radix.bit_mask; ++i)
          {
            k = ((slkey_pure_t) i + 1) << bm->bd.radix.rhigh;
            counts[i] = sl_search_binary_lt_bmask(&xi, k, bm->bd.radix.bit_mask << bm->bd.radix.rhigh);
            Z_TRACE_IF(BR_TRACE_IF, "i = %" slint_fmt " of %" key_pure_type_fmt " << %" slint_fmt ", k = %" key_pure_type_fmt ", searching in %" slint_fmt " -> %" slcount_fmt, i, bm->bd.radix.bit_mask, bm->bd.radix.rhigh, k, xi.size, counts[i]);
            elem_add(&xi, (slint_t) counts[i]);
            xi.size -= counts[i];
          }
          counts[i] = xi.size;
          
          elem_assign_at(&bin->s, special, &xi);
          xi.size = bin->s.size - special;

          for (i = 0; i < bm->bd.radix.bit_mask; ++i)
          {
            k = ((slkey_pure_t) i + 1) << bm->bd.radix.rhigh;
            counts[i + bm->bd.radix.bit_mask + 1] = sl_search_binary_lt_bmask(&xi, k, bm->bd.radix.bit_mask << bm->bd.radix.rhigh);
            Z_TRACE_IF(BR_TRACE_IF, "i = %" slint_fmt " of %" key_pure_type_fmt " << %" slint_fmt ", k = %" key_pure_type_fmt ", searching in %" slint_fmt " -> %" slcount_fmt, i, bm->bd.radix.bit_mask, bm->bd.radix.rhigh, k, xi.size, counts[i + bm->bd.radix.bit_mask + 1]);
            elem_add(&xi, (slint_t) counts[i + bm->bd.radix.bit_mask + 1]);
            xi.size -= counts[i + bm->bd.radix.bit_mask + 1];
          }
          counts[i + bm->bd.radix.bit_mask + 1] = xi.size;

        } else
        {
          elem_assign(&bin->s, &xi);
          for (i = 0; i < bm->bd.radix.bit_mask; ++i)
          {
            k = ((slkey_pure_t) i + 1) << bm->bd.radix.rhigh;
            counts[i] = sl_search_binary_lt_bmask(&xi, k, bm->bd.radix.bit_mask << bm->bd.radix.rhigh);
            Z_TRACE_IF(BR_TRACE_IF, "i = %" slint_fmt " of %" key_pure_type_fmt " << %" slint_fmt ", k = %" key_pure_type_fmt ", searching in %" slint_fmt " -> %" slcount_fmt, i, bm->bd.radix.bit_mask, bm->bd.radix.rhigh, k, xi.size, counts[i]);
            elem_add(&xi, (slint_t) counts[i]);
            xi.size -= counts[i];
          }
          counts[i] = xi.size;
        }
        
      } else
      {
        for (elem_assign(&bin->s, &xi); xi.keys < end.keys; elem_inc(&xi))
        {
          ++counts[key_radix_key2class(key_purify(*xi.keys), bm->bd.radix.rhigh, bm->bd.radix.bit_mask)];
        }
      }
    }

    if (!(bm->sorted & SL_SORTED_IN) && (bm->sorted & SL_SORTED_OUT))
    {
      c = z_alloca(bm->bd.radix.bit_mask + 1, sizeof(slint_t));
      for (i = 0; i < bm->bd.radix.bit_mask + 1; ++i) c[i] = counts[i];
      splitx_radix(&bin->s, &bm->bd.radix.sx, bm->bd.radix.bit_mask + 1, bm->bd.radix.rhigh, c);
      z_freea(c);
    }

  } else
  {
    /* total counts */
    counts[0] += bin->s.size;

#ifdef elem_weight
    /* total weights */
    if (bm->doweights)
    {
      for (elem_assign(&bin->s, &xi); xi.keys < end.keys; elem_inc(&xi)) weights[0] += elem_weight(&xi, 0);
    }
#endif
  }

  Z_TRACE_ARRAY_IF(BR_TRACE_IF, k, bm->nbins, " %" slcount_fmt, counts[k], "counts of %" slint_fmt " @ %p:", bin->s.size, bin->s.keys);
#ifdef elem_weight
  if (bm->doweights)
    Z_TRACE_ARRAY_IF(BR_TRACE_IF, k, bm->nbins, " %" slweight_fmt, weights[k], "weights of %" slint_fmt " @ %p:", bin->s.size, bin->s.keys);
#endif

  return 0;
}


slint_t binning_radix_refine(binning_t *bm, bin_t *bin, slint_t k, slcount_t *counts, slweight_t *weights, splitter_t *sp, slint_t s, bin_t *new_bin) /* sl_proto, sl_func binning_radix_refine */
{
  slint_t l, lcs;

  lcs = 0;
  for (l = 0; l < k; ++l) lcs += counts[l];

  if (new_bin)
  {
    elem_assign_at(&bin->s, lcs, &new_bin->s);
    new_bin->s.size = counts[k];
  
    Z_TRACE_IF(BR_TRACE_IF, "new bin count: %" slint_fmt, new_bin->s.size);

#ifdef elem_weight
    if (bm->doweights)
    {
      new_bin->weight = weights[k];

      Z_TRACE_IF(BR_TRACE_IF, "new bin weight: %" slweight_fmt, new_bin->weight);
    }
#endif
  }

  sp->displs[s] += lcs;

  Z_TRACE_IF(BR_TRACE_IF, "displs[%" slint_fmt "] += %" slint_fmt " = %d", s, lcs, sp->displs[s]);

  return 0;
}


slint_t binning_radix_hit(binning_t *bm, bin_t *bin, slint_t k, slcount_t *counts, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_radix_hit */
{
  slint_t l;

  for (l = 0; l < k; ++l) sp->displs[s] += counts[l];

  Z_TRACE_IF(BR_TRACE_IF, "displs[%" slint_fmt "] += ... = %d", s, sp->displs[s]);

  return 0;
}


slint_t binning_radix_finalize(binning_t *bm, bin_t *bin, slint_t dc, slweight_t dw, slint_t lc_min, slint_t lc_max, slcount_t *lcs, slweight_t *lws, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_radix_finalize */
{
  slint_t lc, r;
#ifdef elem_weight
  elements_t xi, end;
  slweight_t lw, lw_old;
#endif


  Z_TRACE_IF(BR_TRACE_IF, "bin size: %" slint_fmt ", dc = %" slint_fmt ", lc: %" slint_fmt " - %" slint_fmt ", *lcs = %" slcount_fmt, bin->s.size, dc, lc_min, lc_max, *lcs);
#ifdef elem_weight
  if (bm->doweights)
    Z_TRACE_IF(BR_TRACE_IF, "bin weight: %" slweight_fmt ", dw = %" slweight_fmt ", lc: %" slint_fmt " - %" slint_fmt ", *lws = %" slweight_fmt, bin->weight, dw, lc_min, lc_max, *lws);
#endif

  r = 0;

#ifdef elem_weight
  if (bm->doweights)
  {
    lc = 0;
    lw = 0.0;

    if (bin->s.size <= lc_min || (dw >= bin->weight && bin->s.size <= lc_max))
    {
      lc = bin->s.size;
      lw = bin->weight;

    } else
    {
      if (0 < lc_max)
      {
        elem_assign_at(&bin->s, bin->s.size, &end);

        lw = dw;

        for (elem_assign(&bin->s, &xi); xi.keys < end.keys; elem_inc(&xi))
        {
          ++lc;
          lw_old = lw;
          lw -= elem_weight(&xi, 0);
        
          if (lc <= lc_min) continue;

          if (lc > lc_max || (lw < 0.0 && fabs(lw_old) < fabs(lw)))
          {
            lw = lw_old;
            --lc;
            break;
          }
        }
      
        lw = dw - lw;
      }

      r = 1;
    }

  } else
#endif
  {
    lc = z_min(dc, bin->s.size);
    
    r = (lc >= dc);
  }

  *lcs += lc;
  Z_TRACE_IF(BR_TRACE_IF, "*lcs = %" slcount_fmt " + %" slint_fmt " = %" slcount_fmt, (slcount_t) (*lcs - lc), lc, *lcs);
#ifdef elem_weight
  if (bm->doweights)
  {
    *lws += lw;
    Z_TRACE_IF(BR_TRACE_IF, "*lws = %" slweight_fmt " + %" slweight_fmt " = %" slweight_fmt, *lws - lw, lw, *lws);
  }
#endif

  sp->displs[s] += lc;

  Z_TRACE_IF(BR_TRACE_IF, "displs[%" slint_fmt "] += %" slint_fmt " = %d", s, lc, sp->displs[s]);

  return r;
}


slint_t binning_radix_post(binning_t *bm) /* sl_proto, sl_func binning_radix_post */
{
  --bm->bd.radix.rhigh;

  return 0;
}



/* sl_macro E_TRACE_IF */

#include "sl_common.h"


#ifndef E_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define E_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define E_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t elements_alloc(elements_t *s, slint_t nelements, slcint_t components) /* sl_proto, sl_func elements_alloc */
{
  slint_t failed = 0;


  if (s == NULL) return -1;

  elem_null(s);

  if (nelements == 0) return 0;

  /* TODO: set size to 0 and only max_size to nelements (watch out for places where the old behavior is expected!) */
  s->size = s->max_size = nelements;

#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  xelem_buf(s) = z_alloc(nelements, xelem_size_c * sizeof(xelem_type_c)); \
  if (xelem_buf(s) == NULL) failed = 1; \
}
#include "sl_xelem_call.h"

  if (failed)
  {
    elements_free(s);
    return -1;
  }
  
  return 0;
}


slint_t elements_free(elements_t *s) /* sl_proto, sl_func elements_free */
{
  if (s == NULL) return -1;

#define xelem_call          z_free(xelem_buf(s));
#include "sl_xelem_call.h"

  elem_null(s);

  return 0;
}


slint_t elements_realloc(elements_t *s, slint_t nelements, slcint_t components) /* sl_proto, sl_func elements_realloc */
{
  slint_t failed = 0;
  slint_t old_size, old_max_size;


  if (s == NULL) return -1;

  if (nelements == s->max_size) return 0;

  old_size = s->size;
  old_max_size = s->max_size;

  s->size = z_min(s->size, nelements);
  s->max_size = nelements;

#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  xelem_buf(s) = z_realloc(xelem_buf(s), nelements, xelem_size_c * sizeof(xelem_type_c)); \
  if (xelem_buf(s) == NULL) failed = 1; \
}
#include "sl_xelem_call.h"

  if (failed)
  {
    /* try to restore */
    s->size = old_size;
    s->max_size = old_max_size;

#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  xelem_buf(s) = z_realloc(xelem_buf(s), nelements, xelem_size_c * sizeof(xelem_type_c)); \
}
#include "sl_xelem_call.h"
    return -1;
  }

  return 0;
}


slint_t elements_alloca(elements_t *s, slint_t nelements, slcint_t components) /* sl_proto, sl_func elements_alloca */
{
  slint_t failed = 0;


  if (s == NULL) return -1;

  elem_null(s);

  if (nelements == 0) return 0;

  s->size = s->max_size = nelements;

#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  xelem_buf(s) = z_alloca(nelements, xelem_size_c * sizeof(xelem_type_c)); \
  if (xelem_buf(s) == NULL) failed = 1; \
}
#include "sl_xelem_call.h"

  if (failed)
  {
    elements_freea(s);
    return -1;
  }
  
  return 0;
}


slint_t elements_freea(elements_t *s) /* sl_proto, sl_func elements_freea */
{
  if (s == NULL) return -1;

#define xelem_call          z_freea(xelem_buf(s));
#include "sl_xelem_call.h"

  elem_null(s);

  return 0;
}


/* FIXME: alignment! */
slint_t elements_alloc_from_blocks(elements_t *s, slint_t nblocks, void **blocks, slint_t *blocksizes, slint_t alignment, slint_t nmax, slcint_t components) /* sl_proto, sl_func elements_alloc_from_blocks */
{
  slint_t csizes[2 + data_nmax], ctargets[2 + data_nmax];
  slint_t cbytes[nblocks], nelements[nblocks];
  slint_t max_nelements, max_i, cur_min;
  void *_blocks[nblocks];

  slint_t i, j, n;


  if (s == NULL) return -1;

  elem_null(s);

  for (i = 0; i < nblocks; ++i)
  {
    cbytes[i] = 0;
    _blocks[i] = blocks[i];
    nelements[i] = -1;

    Z_TRACE_IF(E_TRACE_IF, "block %" slint_fmt ": %" slint_fmt " @ %p", i, blocksizes[i], _blocks[i]);
  }

  n = 1;
  i = 0;
#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  csizes[i] = xelem_byte; \
  ctargets[i] = 0; \
  cbytes[0] += csizes[i]; \
  n *= nblocks; \
  ++i; \
}
#include "sl_xelem_call.h"

  Z_TRACE_IF(E_TRACE_IF, "n: %" slint_fmt, n);

  if (cbytes[0] > 0) nelements[0] = (_blocks[0])?(blocksizes[0] / cbytes[0]):0;
  else nelements[0] = -1;

  max_nelements = nelements[0];
  max_i = 0;

  Z_TRACE_IF(E_TRACE_IF, "nelements[0]: %" slint_fmt ", max nelements: %" slint_fmt, nelements[0], max_nelements);

  /* find best combination (scales with nblocks^ncomponents) */
  for (i = 1; i < n; ++i)
  {
    if (nmax >= 0 && max_nelements >= nmax) break;

    j = -1;
    do
    {
      ++j;

      cbytes[ctargets[j]] -= csizes[j];
      
      if (cbytes[ctargets[j]] > 0) nelements[ctargets[j]] = (_blocks[ctargets[j]])?(blocksizes[ctargets[j]] / cbytes[ctargets[j]]):0;
      else nelements[ctargets[j]] = -1;

      ctargets[j] = (ctargets[j] + 1) % nblocks;

      cbytes[ctargets[j]] += csizes[j];

      if (cbytes[ctargets[j]] > 0) nelements[ctargets[j]] = (_blocks[ctargets[j]])?(blocksizes[ctargets[j]] / cbytes[ctargets[j]]):0;
      else nelements[ctargets[j]] = -1;

    } while (ctargets[j] == 0);

    cur_min = -1;
    for (j = 0; j < nblocks; ++j) if (cur_min < 0) cur_min = nelements[j]; else if (nelements[j] >= 0) cur_min = z_min(cur_min, nelements[j]);

    if (cur_min > max_nelements)
    {
      max_nelements = cur_min;
      max_i = i;
    }

    Z_TRACE_IF(E_TRACE_IF, "%" slint_fmt ": cur: %" slint_fmt ", max: %" slint_fmt "(%" slint_fmt ")", i, cur_min, max_nelements, max_i);
  }

  Z_ASSERT(max_nelements >= 0);

  Z_TRACE_IF(E_TRACE_IF, "max i: %" slint_fmt ", max nelements: %" slint_fmt, max_i, max_nelements);

  if (nmax >= 0 && max_nelements > nmax) max_nelements = nmax;

  if (max_nelements == 0) return 0;

  s->size = s->max_size = max_nelements;

#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  xelem_buf(s) = (xelem_type_c *) _blocks[max_i % nblocks]; \
  _blocks[max_i % nblocks] = xelem_buf_at(s, max_nelements); \
  max_i /= nblocks; \
}
#include "sl_xelem_call.h"

  return 0;
}


slint_t elements_alloc_from_block(elements_t *s, void *block, slint_t blocksize, slint_t alignment, slint_t nmax, slcint_t components) /* sl_proto, sl_func elements_alloc_from_block */
{
  return elements_alloc_from_blocks(s, 1, &block, &blocksize, alignment, nmax, components);
}


slint_t elements_alloc_block(elements_t *s, void **block, slint_t *blocksize, slint_t alignment, slint_t maxblocksize) /* sl_proto, sl_func elements_alloc_block */
{
  char *next = NULL;


  if (elem_get_block(s))
  {
    *block = elem_get_block(s);
    *blocksize = elem_get_block_size(s);

    return 0;
  }

  *block = NULL;
  *blocksize = 0;

#define xelem_call  if (xelem_buf(s)) \
{ \
  if (next == (char *) xelem_buf(s)) *blocksize += s->max_size * xelem_byte; \
  else if (*blocksize < s->max_size * xelem_byte) \
  { \
    *block = xelem_buf(s); \
    *blocksize = s->max_size * xelem_byte; \
    next = *block; \
  } \
  next += s->max_size * xelem_byte; \
  if (maxblocksize >= 0 && *blocksize >= maxblocksize) goto exit; \
}
#include "sl_xelem_call.h"

exit:

  return 0;
}


slint_t elements_copy(elements_t *s, elements_t *d) /* sl_proto, sl_func elements_copy */
{
  elem_copy(s, d);
  
  return 0;
}


slint_t elements_copy_at(elements_t *s, slint_t sat, elements_t *d, slint_t dat) /* sl_proto, sl_func elements_copy_at */
{
  elem_copy_at(s, sat, d, dat);
  
  return 0;
}


slint_t elements_ncopy(elements_t *s, elements_t *d, slint_t n) /* sl_proto, sl_func elements_ncopy */
{
  elem_ncopy(s, d, n);
  
  return 0;
}


slint_t elements_nmove(elements_t *s, elements_t *d, slint_t n) /* sl_proto, sl_func elements_nmove */
{
  elem_nmove(s, d, n);
  
  return 0;
}


slint_t elements_printf(elements_t *s, const char *prefix) /* sl_proto, sl_func elements_printf */
{
  if (s == NULL) return -1;

  printf("%s: [%" sl_int_type_fmt ", %" sl_int_type_fmt
#define xelem_call          ", %p"
#include "sl_xelem_call.h"
    "]\n", prefix, s->size, s->max_size
#define xelem_call          , xelem_buf(s)
#include "sl_xelem_call.h"
    );

  return 0;
}


slint_t elements_extract(elements_t *src, slint_t nelements, elements_t *dst0, elements_t *dst1) /* sl_proto, sl_func elements_extract */
{
  elements_t s;

  if (src == NULL) return -1;

  s = *src;

  if (dst0 != NULL)
  {
    elem_assign(&s, dst0);
    dst0->size = z_min(s.size, nelements);

    if (dst0->size <= 0) elem_null(dst0);
  }

  if (dst1 != NULL)
  {
    elem_assign_at(&s, nelements, dst1);
    dst1->size = z_max(s.size - nelements, 0);

    if (dst1->size <= 0) elem_null(dst1);
  }

  return 0;
}


slint_t elements_touch(elements_t *s) /* sl_proto, sl_func elements_touch */
{
  elements_t _s, end, t;

  elements_alloc(&t, 1, SLCM_ALL);

  elem_assign_at(s, s->size, &end);

  for (elem_assign(s, &_s); _s.keys < end.keys; elem_inc(&_s)) elem_copy(&_s, &t);

  elements_free(&t);

  return 0;
}


slint_t elements_digest_sum(elements_t *s, slint_t nelements, slcint_t components, unsigned int *sum) /* sl_proto, sl_func elements_digest_sum */
{
  slint_t i, j;
  unsigned int ssum = 0;

#ifdef ELEMENTS_DIGEST_SUM_ADD
# define SUM_OP(_a_,_b_)  _a_ + _b_
#else
# define SUM_OP(_a_,_b_)  _a_ ^ _b_
#endif

  *sum = 0;
  for (j = 0; j < nelements; j++)
  {
    for (i = 0; i < s[j].size; i++)
    {
#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  z_digest_sum_buffer((const void *) xelem_buf_at(&s[j], i), xelem_byte, (void *) &ssum); \
  *sum = SUM_OP(*sum, ssum); \
}
#include "sl_xelem_call.h"
    }
  }

#undef SUM_OP

  return 0;
}


/* deprecated */
unsigned int elements_crc32(elements_t *s, slint nelements, slint_t keys, slint_t data) /* sl_proto, sl_func elements_crc32 */
{
  unsigned int crc32;
  
  elements_digest_sum(s, nelements, (keys?(SLCM_KEYS):0)|(data?(SLCM_DATA):0), &crc32);

  return crc32;
}


slint_t elements_digest_hash(elements_t *s, slint_t nelements, slcint_t components, void *hash) /* sl_proto, sl_func elements_digest_hash */
{
  slint_t i, j;

  if (!hash) return z_digest_hash_read(NULL, NULL);

  void *hdl;

  z_digest_hash_open(&hdl);

  for (j = 0; j < nelements; j++)
  {
    for (i = 0; i < s[j].size; i++)
    {
#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  z_digest_hash_write(hdl, (const void *) xelem_buf_at(&s[j], i), xelem_byte); \
}
#include "sl_xelem_call.h"
    }
  }

  z_digest_hash_read(hdl, hash);

  z_digest_hash_close(hdl);

  return 0;
}


slint_t elements_random_exchange(elements_t *s, slint_t rounds, elements_t *xs) /* sl_proto, sl_func elements_random_exchange */
{
  slint_t i, j, k = 0;
  elements_t txs;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  j = 0;
  elem_copy(s, xs);

  for (i = 0; i < rounds; i++)
  {
    k = z_rand() % s->size;
    elem_copy_at(s, k, s, j);
    j = k;
  }

  elem_copy_at(xs, 0, s, k);

  if (xs == &txs) elements_free(xs);

  return 0;
}


slint_t elements_keys_init_seed(unsigned long s) /* sl_proto, sl_func elements_keys_init_seed */
{
#ifdef key_val_srand
  key_val_srand(s);
#endif

  z_urandom_seed(s);
  z_nrandom_seed(s);

  return 0;
}

#define KEY_SET_VARIABLES() \
  slint_t ks_j, ks_m; \
  slkey_pure_t *ks_pk; \
  key_set_f ks_func; \
  key_set_data_t ks_data

#define KEY_SET_SET_INIT(_d_)       Z_NOP()
#define KEY_SET_SET(_k_, _d_)       Z_MOP(key_set_pure((_k_), *((slkey_pure_t *) (_d_)));)

#define KEY_SET_SET_FUNC_INIT(_d_)  Z_MOP(ks_func = (key_set_f) ((void **) (_d_))[0]; ks_data = (key_set_data_t) ((void **) (_d_))[1];)
#define KEY_SET_SET_FUNC(_k_, _d_)  Z_MOP(ks_func(key_get_pure(_k_), ks_data);)

#define KEY_SET_RAND_INIT(_d_)      Z_MOP( \
  if (_d_) { \
    if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand_minmax is missing!"); \
    ks_pk = (slkey_pure_t *) (_d_); \
  } else { \
    if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand is missing!"); \
  })
#define KEY_SET_RAND(_k_, _d_)      Z_MOP(if (_d_) key_set_pure((_k_), key_val_rand_minmax(ks_pk[0], ks_pk[1])); else key_set_pure((_k_), key_val_rand());)

#define KEY_SET_RAND_QUAD_INIT(_d_)  Z_MOP( \
  if (_d_) { \
    if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand_minmax is missing!"); \
    ks_pk = (slkey_pure_t *) (_d_); \
  } else { \
    if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand is missing!"); \
  })
#define KEY_SET_RAND_QUAD(_k_, _d_)  Z_MOP(if (_d_) key_set_pure((_k_), key_val_rand_minmax(ks_pk[0], ks_pk[1])); else key_set_pure((_k_), key_val_rand()); key_set_pure((_k_), (key_purify(*(_k_)) * key_purify(*(_k_))));)

#ifdef key_integer
#define KEY_SET_RAND_AND_INIT(_d_)  Z_MOP( \
  ks_m = *((slint_t *) (_d_)); \
  if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand is missing!");)
#define KEY_SET_RAND_AND(_k_, _d_)  Z_MOP(key_set_pure((_k_), key_val_rand()); for (ks_j = 0; ks_j < ks_m; ++ks_j) key_set_pure((_k_), (key_purify(*(_k_)) & key_val_rand()));)
#endif

#define KEY_SET_URAND_INIT(_d_)  Z_MOP(ks_pk = (slkey_pure_t *) (_d_);)
#define KEY_SET_URAND(_k_, _d_)  Z_MOP(key_set_pure((_k_), (slkey_pure_t) (ks_pk[0] + z_urandom() * (ks_pk[1] - ks_pk[0]) + 0.5));)

#define KEY_SET_NRAND_INIT(_d_)  Z_MOP(ks_pk = (slkey_pure_t *) (_d_);)
#define KEY_SET_NRAND(_k_, _d_)  Z_MOP(key_set_pure((_k_), (slkey_pure_t) (ks_pk[2] + z_nrandom() * ks_pk[3])); key_set_pure((_k_), z_minmax(ks_pk[0], key_purify(*(_k_)), ks_pk[1]));)


slint_t elements_keys_init(elements_t *s, keys_init_type_t t, keys_init_data_t d) /* sl_proto, sl_func elements_keys_init */
{
  slint_t i;

  KEY_SET_VARIABLES();

  if (s == NULL) return -1;

  switch (t)
  {
    case SL_EKIT_SET:
      KEY_SET_SET_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_SET(&s->keys[i], d);
      break;
    case SL_EKIT_SET_FUNC:
      KEY_SET_SET_FUNC_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_SET_FUNC(&s->keys[i], d);
      break;
    case SL_EKIT_RAND:
      KEY_SET_RAND_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_RAND(&s->keys[i], d);
      break;
    case SL_EKIT_RAND_QUAD:
      KEY_SET_RAND_QUAD_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_RAND_QUAD(&s->keys[i], d);
      break;
#ifdef key_integer
    case SL_EKIT_RAND_AND:
      KEY_SET_RAND_AND_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_RAND_AND(&s->keys[i], d);
      break;
#endif
    case SL_EKIT_URAND:
      KEY_SET_URAND_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_URAND(&s->keys[i], d);
      break;
    case SL_EKIT_NRAND:
      KEY_SET_NRAND_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_NRAND(&s->keys[i], d);
      break;
  }

  return 0;
}


slint_t elements_keys_init_randomized(elements_t *s, slint_t nkeys, keys_init_type_t t, keys_init_data_t d) /* sl_proto, sl_func elements_keys_init_randomized */
{
  slint_t i, j;

  KEY_SET_VARIABLES();

  if (s == NULL) return -1;

  i = 0;
  switch (t)
  {
    case SL_EKIT_SET:
      KEY_SET_SET_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_SET(&s->keys[i], d); }
      break;
    case SL_EKIT_SET_FUNC:
      KEY_SET_SET_FUNC_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_SET_FUNC(&s->keys[i], d); }
      break;
    case SL_EKIT_RAND:
      KEY_SET_RAND_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_RAND(&s->keys[i], d); }
      break;
    case SL_EKIT_RAND_QUAD:
      KEY_SET_RAND_QUAD_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_RAND_QUAD(&s->keys[i], d); }
      break;
#ifdef key_integer
    case SL_EKIT_RAND_AND:
      KEY_SET_RAND_AND_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_RAND_AND(&s->keys[i], d); }
      break;
#endif
    case SL_EKIT_URAND:
      KEY_SET_URAND_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_URAND(&s->keys[i], d); }
      break;
    case SL_EKIT_NRAND:
      KEY_SET_NRAND_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_NRAND(&s->keys[i], d); }
      break;
  }

  return 0;
}


#undef KEY_SET_VARIABLES
#undef KEY_SET_SET_INIT
#undef KEY_SET_SET
#undef KEY_SET_SET_FUNC_INIT
#undef KEY_SET_SET_FUNC
#undef KEY_SET_RAND_INIT
#undef KEY_SET_RAND
#undef KEY_SET_RAND_QUAD_INIT
#undef KEY_SET_RAND_QUAD
#undef KEY_SET_RAND_AND_INIT
#undef KEY_SET_RAND_AND
#undef KEY_SET_URAND_INIT
#undef KEY_SET_URAND
#undef KEY_SET_NRAND_INIT
#undef KEY_SET_NRAND


#define LINE_LENGTH  1024

slint_t elements_keys_init_from_file(elements_t *s, slint_t data, char *filename, slint_t from, slint_t to, slint_t const_bytes_per_line) /* sl_proto, sl_func elements_keys_init_from_file */
{
  FILE *inputfile;
  char buffer[LINE_LENGTH];
  slint_t i = 0, line = 0;
  slint_t bytes_per_line;
  slkey_pure_t inkey;

  elements_alloc(s, to - from + 1, SLCM_ALL|((data)?0:(~SLCM_DATA)));

/*  printf("opening '%s'\n", filename);*/

  inputfile = fopen(filename, "r");

/*  printf("inputfile = %p\n", inputfile);*/

  if (!inputfile) { return -1; }

#ifdef key_integer
  if (const_bytes_per_line)
  {
    fgets(buffer, LINE_LENGTH, inputfile);
    bytes_per_line = ftell(inputfile);
    rewind(inputfile);

    fseek(inputfile, from * bytes_per_line, SEEK_SET);

    line = from;
    
  } else while (line < from)
  {
    line++;
    if (!fgets(buffer, LINE_LENGTH, inputfile)) break;
  }

  while((i < s->size) && (line <= to))
  {
    if (!fgets(buffer, LINE_LENGTH, inputfile)) break;

#ifdef key_pure_type_fmt
    sscanf(buffer, "%" key_pure_type_fmt, &inkey);
    key_set_pure(&s->keys[i], inkey);
#endif

/*    printf("line: %d - input: '%s'", line, buffer);
    printf("i sscanf'd %d: %ld\n", r, s->keys[i]);*/

    ++line;

    ++i;
  }
#endif


  fclose(inputfile);

  return i;
}

#undef LINE_LENGTH


slint_t elements_keys_save_to_file(elements_t *s, char *filename) /* sl_proto, sl_func elements_keys_save_to_file */
{
  FILE *outputfile;
  slint_t i = 0;

  printf("opening '%s'\n", filename);

  outputfile = fopen(filename, "w");

  printf("inputfile = %p\n", outputfile);

  if (!outputfile) { return -1; }

  for (i = 0; i < s->size; ++i)
#ifdef key_pure_type_fmt
    fprintf(outputfile, "%" key_pure_type_fmt "\n", key_purify(s->keys[i]));
#else
    ;
#endif

  fclose(outputfile);

  return 0;
}


#define evo_body \
  slint_t i, j, k = 0, l = -1; \
  if (s == NULL) return -1; \
  for (j = 0; j < n; j++) \
  { \
    k++; \
    if ((l >= 0) && (s[j].size > 0)) \
    if (key_pure_cmp_lt(km(s[j].keys[0]), km(s[l].keys[s[l].size - 1]))) return j; \
    for (i = 1; i < s[j].size; i++, k++) \
    if (key_pure_cmp_lt(km(s[j].keys[i]), km(s[j].keys[i - 1]))) return k; \
    if (s[j].size > 0) l = j; \
  }

slint_t elements_validate_order(elements_t *s, slint_t n) /* sl_proto, sl_func elements_validate_order */
{

#define km(k) key_purify(k)

  evo_body

#undef km

  return 0;
}


slint_t elements_validate_order_bmask(elements_t *s, slint_t n, slkey_pure_t bmask) /* sl_proto, sl_func elements_validate_order_bmask */
{

#define  km(k) (key_purify(k) & bmask)

  evo_body

#undef km

  return 0;
}


slint_t elements_validate_order_weight(elements_t *s, slint_t n, slkey_pure_t weight) /* sl_proto, sl_func elements_validate_order_weight */
{

#define  km(k) (key_purify(k) / weight)

  evo_body

#undef km

  return 0;
}

#undef evo_body


#if defined(HAVE_GMP_H)
# include <gmp.h>
#endif


slint_t elements_keys_stats(elements_t *s, slkey_pure_t *stats) /* sl_proto, sl_func elements_keys_stats */
{
  slint_t i;
  slkey_pure_t k, kmin, kmax;


  if (s->size <= 0) return 0;

#ifdef HAVE_GMP_H
# ifdef key_integer
#  define stats_t                                    mpz_t
#  define stats_init(_x_)                            mpz_init(_x_)
#  define stats_free(_x_)                            mpz_clear(_x_)
#  define stats_set(_x_, _v_)                        Z_MOP(if (sizeof(_v_) <= sizeof(long)) \
                                                           { if (key_integer_unsigned) mpz_set_ui(_x_, (unsigned long) _v_); else mpz_set_si(_x_, (long) _v_); } \
                                                           else \
                                                           { if (key_integer_unsigned) z_gmp_mpz_set_ull(_x_, (unsigned long long) _v_); else z_gmp_mpz_set_sll(_x_, (long long) _v_); })
#  define stats_get(_x_, _type_)                     (_type_) ((sizeof(_type_) <= sizeof(long))? \
                                                               (key_integer_unsigned?mpz_get_ui(_x_):mpz_get_si(_x_)): \
                                                               (key_integer_unsigned?z_gmp_mpz_get_ull(_x_):z_gmp_mpz_get_sll(_x_)))
#  define stats_add(_x_, _v_)                        mpz_add(_x_, _x_, _v_);
#  define stats_addsqr(_x_, _v_, _t_)                Z_MOP(mpz_mul(_t_, _v_, _v_); mpz_add(_x_, _x_, _t_);)
#  define stats_avg(_sum_, _n_, _t_, _type_)         (_type_) ((mpz_tdiv_q_ui(_t_, _sum_, (unsigned long) _n_) < (_n_ / 2))? \
                                                               (stats_get(_t_, _type_)): \
                                                               (stats_get(_t_, _type_) + 1))
#  define stats_std(_sum_, _sqr_, _n_, _t_, _type_)  (_type_) (mpz_mul(_t_, _sum_, _sum_), \
                                                               mpz_tdiv_q_ui(_t_, _t_, (unsigned long) _n_), \
                                                               mpz_sub(_t_, _sqr_, _t_), \
                                                               mpz_tdiv_q_ui(_t_, _t_, (unsigned long) (_n_ - 1)), \
                                                               mpz_sqrt(_t_, _t_), stats_get(_t_, _type_))
#  define stats_print(_x_)                           gmp_printf(#_x_ ": %Zd\n", _x_)
# endif
#else
# define stats_t                                     long double
# define stats_init(_x_)                             Z_NOP()
# define stats_free(_x_)                             Z_NOP()
# define stats_set(_x_, _v_)                         _x_ = (stats_t) _v_
# define stats_get(_x_, _type_)                      (_type_) _x_
# define stats_add(_x_, _v_)                         _x_ += (stats_t) _v_
# define stats_addsqr(_x_, _v_, _t_)                 _t_ = (stats_t) _v_ * (stats_t) _v_, _x_ += _t_
# define stats_avg(_sum_, _n_, _t_, _type_)          (_type_) (_sum_ / _n_ + 0.5)
# define stats_std(_sum_, _sqr_, _n_, _t_, _type_)   (_type_) sqrtl((_sqr_ - ((_sum_ * _sum_) / (stats_t) _n_)) / (stats_t) (_n_ - 1))
# define stats_print(_x_)                            printf(#_x_ ": %Lf\n", _x_)
#endif

  stats_t x, xsum, xsqr, xt;

/*  elements_print_keys(s);*/

  kmin = kmax = key_purify(s->keys[0]);

  stats_init(x);
  stats_init(xsum);
  stats_init(xsqr);
  stats_init(xt);

  stats_set(x, 0);
  stats_set(xsum, 0);
  stats_set(xsqr, 0);
  stats_set(xt, 0);

  for (i = 0; i < s->size; i++)
  {
    k = key_purify(s->keys[i]);

/*    printf("key: %" sl_key_type_fmt "\n", k);*/

    if (k < kmin) kmin = k;
    if (k > kmax) kmax = k;

    stats_set(x, k);
    stats_add(xsum, x);
    stats_addsqr(xsqr, x, xt);

/*    stats_print(x);
    stats_print(xsum);
    stats_print(xsqr);*/
  }

  stats[SL_EKS_MIN] = kmin;
  stats[SL_EKS_MAX] = kmax;
  stats[SL_EKS_SUM] = stats_get(xsum, slkey_pure_t);
  stats[SL_EKS_AVG] = stats_avg(xsum, s->size, xt, slkey_pure_t);

/*  stats_print(ksum);
  stats_print(ksqr);*/

  if (s->size > 1) stats[SL_EKS_STD] = stats_std(xsum, xsqr, s->size, xt, slkey_pure_t);
  else stats[SL_EKS_STD] = 0;

  stats_free(x);
  stats_free(xsum);
  stats_free(xsqr);
  stats_free(xt);

  return 0;
}

#undef stats_t
#undef stats_init
#undef stats_free
#undef stats_set
#undef stats_get
#undef stats_add
#undef stats_addsqr
#undef stats_avg
#undef stats_std
#undef stats_print


slint_t elements_keys_stats_print(elements_t *s) /* sl_proto, sl_func elements_keys_stats_print */
{
  slkey_pure_t stats[SL_EKS_SIZE];


  elements_keys_stats(s, stats);

  printf("min: %" key_pure_type_fmt "\n", stats[SL_EKS_MIN]);
  printf("max: %" key_pure_type_fmt "\n", stats[SL_EKS_MAX]);
  printf("sum: %" key_pure_type_fmt "\n", stats[SL_EKS_SUM]);
  printf("avg: %" key_pure_type_fmt "\n", stats[SL_EKS_AVG]);
  printf("std: %" key_pure_type_fmt "\n", stats[SL_EKS_STD]);

  return 0;
}


slint_t elements_print_keys(elements_t *s) /* sl_proto, sl_func elements_print_keys */
{
  slint_t i;

  if (s == NULL) return -1;

  for (i = 0; i < s->size; i++)
  {
    printf(" [%3" sl_int_type_fmt "] @ %p = ", i, &s->keys[i]);
/*    key_printf(s->keys[i]);*/
#ifdef key_pure_type_fmt
    printf("%" key_pure_type_fmt "\n", key_purify(s->keys[i]));
#else
    printf("\n");
#endif
  }

  return 0;
}


slint_t elements_print_all(elements_t *s) /* sl_proto, sl_func elements_print_all */
{
  slint_t i;

  if (s == NULL) return -1;

  for (i = 0; i < s->size; i++)
  {
    printf(" [%3" slint_fmt "]", i);
#ifdef SL_INDEX
    printf(", idx @ %p = %" sl_index_type_fmt, &s->indices[i], s->indices[i]);
#endif
#ifdef key_pure_type_fmt
    printf(", key @ %p = %" key_pure_type_fmt, &s->keys[i], key_purify(s->keys[i]));
#endif
#ifdef elem_weight
    printf(", weight = %" slweight_fmt, elem_weight(s, i));
#endif
    printf("\n");
  }

  return 0;
}


slweight_t elements_get_weight(elements_t *s) /* sl_proto, sl_func elements_get_weight */
{
#ifdef elem_weight
  slint_t i;
  slweight_t w = 0.0;

  for (i = 0; i < s->size; ++i) w += elem_weight(s, i);

  return w;
#else
  return 0.0;
#endif
}


slint_t elements_get_minmax_keys(elements_t *s, slint_t nelements, slkey_pure_t *minmaxkeys) /* sl_proto, sl_func elements_get_minmax_keys */
{
  slint_t i, j;

  if (s == NULL || nelements < 1) return -1;
  
  minmaxkeys[0] = minmaxkeys[1] = key_purify(s[0].keys[0]);
  
  for (j = 0; j < nelements; ++j)
  for (i = 0; i < s[j].size; ++i)
  {
    if (key_pure_cmp_lt(key_purify(s[j].keys[i]), minmaxkeys[0])) minmaxkeys[0] = key_purify(s[j].keys[i]);
    if (key_pure_cmp_gt(key_purify(s[j].keys[i]), minmaxkeys[1])) minmaxkeys[1] = key_purify(s[j].keys[i]);
  }

  return 0;
}



#include "sl_common.h"


slint_t elements_alloc_packed(packed_elements_t *s, slint_t nelements) /* sl_proto, sl_func elements_alloc_packed */
{
  if (s == NULL) return -1;

  s->size = s->max_size = 0;
  s->elements = NULL;

  if (nelements == 0) return 0;

  s->size = s->max_size = nelements;
  
  s->elements = z_alloc(nelements, pelem_byte);

  if (s->elements != NULL) return 0;

  elements_free_packed(s);

  return -1;
}


slint_t elements_free_packed(packed_elements_t *s) /* sl_proto, sl_func elements_free_packed */
{
  if (s == NULL) return -1;

  z_free(s->elements);

  s->size = s->max_size = 0;

  s->elements = NULL;

  return 0;
}


slint_t elements_alloc_packed_from_block(packed_elements_t *s, void *block, slint_t blocksize, slint_t alignment, slint_t nmax) /* sl_proto, sl_func elements_alloc_packed_from_block */
{
  if (s == NULL) return -1;

  s->size = s->max_size = blocksize / pelem_byte;

  if (nmax >= 0 && s->size > nmax) s->size = s->max_size = nmax;

  s->elements = block;

  return 0;
}


slint_t elements_pack_indexed(elements_t *s, packed_elements_t *d, slindex_t *rindx, slindex_t *windx) /* sl_proto, sl_func elements_pack_indexed */
{
  slint_t i;

  if (s->size > d->size) return -1;

#define DO_INDEXED \
  for (i = 0; i < s->size; ++i) \
  { \
    key_copy_at(s->keys, READ_INDEX(i), &d->elements[WRITE_INDEX(i)].key, 0); \
    data_copy_at(s, READ_INDEX(i), &d->elements[WRITE_INDEX(i)], 0); \
  }

  if (rindx == NULL)
  {
    if (windx == NULL)
    {
#define READ_INDEX(_i_)   _i_
#define WRITE_INDEX(_i_)  _i_
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX

    } else
    {
#define READ_INDEX(_i_)   _i_
#define WRITE_INDEX(_i_)  windx[_i_]
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX
    }

  } else
  {
    if (windx == NULL)
    {
#define READ_INDEX(_i_)   rindx[_i_]
#define WRITE_INDEX(_i_)  _i_
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX

    } else
    {
#define READ_INDEX(_i_)   rindx[_i_]
#define WRITE_INDEX(_i_)  windx[_i_]
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX
    }
  }

#undef DO_INDEXED

  return 0;
}


slint_t elements_pack(elements_t *s, packed_elements_t *d) /* sl_proto, sl_func elements_pack */
{
  return elements_pack_indexed(s, d, NULL, NULL);
}


slint_t elements_pack_at(elements_t *s, slint_t sat, packed_elements_t *d, slint_t dat) /* sl_proto, sl_func elements_pack_at */
{
  slint_t i;

  if (s->size > d->size) return -1;

  for (i = 0; i < s->size; ++i)
  {
    key_copy_at(s->keys, sat + i, &d->elements[dat + i].key, 0);
    data_copy_at(s, sat +i, &d->elements[dat + i], 0);
  }

  return 0;
}


slint_t elements_unpack_indexed(packed_elements_t *s, elements_t *d, slindex_t *rindx, slindex_t *windx) /* sl_proto, sl_func elements_unpack_indexed */
{
  slint_t i;

  if (s->size > d->size) return -1;

#define DO_INDEXED \
  for (i = 0; i < s->size; ++i) \
  { \
    key_copy_at(&s->elements[READ_INDEX(i)].key, 0, d->keys, WRITE_INDEX(i)); \
    data_copy_at(&s->elements[READ_INDEX(i)], 0, d, WRITE_INDEX(i)); \
  }

  if (rindx == NULL)
  {
    if (windx == NULL)
    {
#define READ_INDEX(_i_)   _i_
#define WRITE_INDEX(_i_)  _i_
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX

    } else
    {
#define READ_INDEX(_i_)   _i_
#define WRITE_INDEX(_i_)  windx[_i_]
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX
    }

  } else
  {
    if (windx == NULL)
    {
#define READ_INDEX(_i_)   rindx[_i_]
#define WRITE_INDEX(_i_)  _i_
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX

    } else
    {
#define READ_INDEX(_i_)   rindx[_i_]
#define WRITE_INDEX(_i_)  windx[_i_]
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX
    }
  }

#undef DO_INDEXED

  return 0;
}


slint_t elements_unpack(packed_elements_t *s, elements_t *d) /* sl_proto, sl_func elements_unpack */
{
  return elements_unpack_indexed(s, d, NULL, NULL);
}


slint_t elements_unpack_at(packed_elements_t *s, slint_t sat, elements_t *d, slint_t dat) /* sl_proto, sl_func elements_unpack_at */
{
  slint_t i;

  if (s->size > d->size) return -1;

  for (i = 0; i < s->size; ++i)
  {
    key_copy_at(&s->elements[sat + i].key, 0, d->keys, dat + i);
    data_copy_at(&s->elements[sat + i], 0, d, dat + i);
  }
  
  return 0;
}


slint_t elements_unpack_keys(packed_elements_t *s, slkey_t *k) /* sl_proto, sl_func elements_unpack_keys */
{
  slint_t i;

  for (i = 0; i < s->size; ++i) key_copy_at(&s->elements[i].key, 0, k, i);
  
  return 0;
}



#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "z_pack.h"
#include "local_generic_heap.h"


#define LGH_TRACE_IF  (z_mpi_rank == -1)


void lgh_create(lgh_t *lgh, lghint_t size) /* lgh_func lgh_create */
{
  Z_TRACE_IF(LGH_TRACE_IF, "size: %" lghint_fmt, size);

  lgh->nallocs = 0;

  lgh->segments = z_alloc(1, sizeof(lgh_segment_t));
  lgh->segments->next = NULL;
  lgh->segments->offset = 0;
  lgh->segments->size = size;
}


void lgh_destroy(lgh_t *lgh) /* lgh_func lgh_destroy */
{
  lgh_segment_t *seg;


  Z_ASSERT(lgh->nallocs == 0);
  Z_ASSERT(lgh->segments != NULL);
  Z_ASSERT(lgh->segments != NULL && lgh->segments->next == NULL);

  Z_TRACE_IF(LGH_TRACE_IF, "nallocs: %" lghint_fmt, lgh->nallocs);

  lgh->nallocs = 0;
  while (lgh->segments)
  {
    seg = lgh->segments->next;
    z_free(lgh->segments);
    lgh->segments = seg;
  }
}


lgh_segment_t *lgh_alloc(lgh_t *lgh, lghint_t size) /* lgh_func lgh_alloc */
{
  return lgh_alloc_minmax(lgh, size, size);
}


lgh_segment_t *lgh_alloc_minmax(lgh_t *lgh, lghint_t min, lghint_t max) /* lgh_func lgh_alloc_minmax */
{
  lgh_segment_t *curr, *prev, *best, *best_prev, *seg;
  lghint_t best_size;


  if (min > max) return NULL;

  Z_TRACE_IF(LGH_TRACE_IF, "searching for minmax: %" lghint_fmt "-%" lghint_fmt, min, max);

  best = best_prev = NULL;
  best_size = min - 1;

  prev = NULL;
  curr = lgh->segments;
  while (curr != NULL && best_size != max)
  {
    if (best_size < max)
    {
      if (curr->size > best_size)
      {
        best = curr;
        best_prev = prev;
        best_size = curr->size;
      }

    } else
    {
      if (max <= curr->size && curr->size < best_size)
      {
        best = curr;
        best_prev = prev;
        best_size = curr->size;
      }
    }

    prev = curr;
    curr = curr->next;

/*    if (!local_heap_alloc_search) continue;*/
  }

  max = z_min(max, best_size);
  curr = best;
  prev = best_prev;

  if (curr == NULL) return NULL;

  if (curr->size <= max)
  {
    seg = curr;

    if (prev == NULL) lgh->segments = curr->next;
    else prev->next = curr->next;

  } else
  {
    seg = z_alloc(1, sizeof(lgh_segment_t));
    seg->next = NULL;
    seg->offset = curr->offset;
    seg->size = max;

    curr->offset += max;
    curr->size -= max;
  }

  ++lgh->nallocs;

  Z_TRACE_IF(LGH_TRACE_IF, "allocated segment %p: [%" lghint_fmt ", %" lghint_fmt "]", seg, seg->offset, seg->size);

  Z_TRACE_IF(LGH_TRACE_IF, "nallocs OUT: %" lghint_fmt, lgh->nallocs);

  return seg;
}


void lgh_free(lgh_t *lgh, lgh_segment_t *seg) /* lgh_func lgh_free */
{
  lgh_segment_t *prev, *curr;


  if (seg == NULL) return;

  Z_TRACE_IF(LGH_TRACE_IF, "nallocs IN: %" lghint_fmt, lgh->nallocs);

  Z_TRACE_IF(LGH_TRACE_IF, "free seg: %p: [%" lghint_fmt ", %" lghint_fmt"]", seg, seg->offset, seg->size);

  prev = NULL;
  curr = lgh->segments;

/*  if (local_heap_free_search)*/
  while (curr != NULL && curr->offset < seg->offset)
  {
    prev = curr;
    curr = curr->next;
  }

  if (curr != NULL && seg->offset + seg->size == curr->offset)
  {
    curr->offset -= seg->size;
    curr->size += seg->size;

    z_free(seg);

    seg = curr;

  } else seg->next = curr;

  if (prev == NULL) lgh->segments = seg;
  else
  {
    if (prev->offset + prev->size == seg->offset)
    {
      prev->size += seg->size;
      prev->next = seg->next;

      z_free(seg);

    } else prev->next = seg;
  }

  --lgh->nallocs;

  Z_TRACE_IF(LGH_TRACE_IF, "nallocs OUT: %" lghint_fmt, lgh->nallocs);
}


#undef LGH_TRACE_IF



#include "sl_common.h"


slint merge2_basic_auto_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_auto_01_x */
{
  if (z_min(s0->size, s1->size) <= slint_sqrt(s0->size + s1->size)) return merge2_basic_sbin_01_x(s0, s1, sx);

  return merge2_basic_straight_01_x(s0, s1, sx);
}


slint merge2_basic_01_x(elements_t *s0, elements_t *s1, elements_t *sx, m2x_func _x0_1, m2x_func _0x_1) /* sl_proto, sl_func merge2_basic_01_x */
{
  elements_t s, e, d, x;

  /* if one list is empty, there is nothing to do */
  if (s0->size == 0 || s1->size == 0) return 0;

  if (sx == NULL)
  {
    elements_alloc(&x, z_min(s0->size, s1->size), SLCM_ALL);
    sx = &x;

  } else if (sx->size < z_min(s0->size, s1->size)) return -1;

  if (s0->size <= s1->size)
  {
    /* initialize */
    elem_assign(s0, &s); elem_assign_at(s0, s0->size, &e);
    elem_assign(sx, &d);

    /* skip already well-placed elements of s0 */
    while (s.keys != e.keys)
    if (key_cmp_le(*s.keys, *e.keys)) elem_inc(&s); else break;

    /* evacuate wrong-placed elements of s0 */
    d.size = s.size = e.keys - s.keys;
    elem_ncopy(&s, &d, d.size);

    /* call merge2 subroutine */
    (_x0_1)(s1, &d, &s);

  } else
  {
    /* initialize */
    elem_assign_at(s1, s1->size - 1, &s); elem_assign_at(s1, -1, &e);
    elem_assign(sx, &d);

    /* skip already well-placed elements of s1 */
    while (s.keys != e.keys)
    if (key_cmp_gt(*s.keys, *e.keys)) elem_dec(&s); else break;

    /* evacuate wrong-placed elements of s1 */
    d.size = e.size = s.keys - e.keys;
    elem_inc(&e);
    elem_ncopy(s1, &d, d.size);

    /* call merge2 subroutine */
    (_0x_1)(s0, &d, &e);
  }

  if (sx == &x) elements_free(&x);

  return 0;
}


slint merge2_basic_01_X(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t, m2X_func _X0_1, m2X_func _0X_1) /* sl_proto, sl_func merge2_basic_01_X */
{
  elements_t s, e, d, x;

  /* if one list is empty, there is nothing to do */
  if (s0->size == 0 || s1->size == 0) return 0;

  if (sx->size < z_min(s0->size, s1->size)) return -1;

  if (t == NULL)
  {
    elements_alloc(&x, 1, SLCM_ALL);
    t = &x;
  }

  if (s0->size <= s1->size)
  {
    /* initialize */
    elem_assign(s0, &s); elem_assign_at(s0, s0->size, &e);
    elem_assign(sx, &d);

    /* skip already well-placed elements of s0 */
    while (s.keys != e.keys)
    if (key_cmp_le(*s.keys, *e.keys)) elem_inc(&s); else break;

    /* evacuate wrong-placed elements of s0 */
    d.size = s.size = e.keys - s.keys;
    elem_nxchange_ro0(&s, &d, d.size, t);

    /* call merge2 subroutine */
    (_X0_1)(s1, &d, &s, t);

  } else
  {
    /* initialize */
    elem_assign_at(s1, s1->size - 1, &s); elem_assign_at(s1, -1, &e);
    elem_assign(sx, &d);

    /* skip already well-placed elements of s1 */
    while (s.keys != e.keys)
    if (key_cmp_gt(*s.keys, *e.keys)) elem_dec(&s); else break;

    /* evacuate wrong-placed elements of s1 */
    d.size = e.size = s.keys - e.keys;
    elem_inc(&e);
    elem_nxchange_ro0(s1, &d, d.size, t);

    /* call merge2 subroutine */
    (_0X_1)(s0, &d, &e, t);
  }

  if (t == &x) elements_free(&x);

  return 0;
}


#define the_merge2_01(s0, s1, sx)  merge2_basic_sseq_01(s0, s1, sx)


/*slint merge2_simplify_s0(elements_t *s0, elements_t *s1, elements_t *sx, slint s1elements)
{
  return 0;
}*/


slint merge2_simplify_s1(elements_t *s0, elements_t *s1, elements_t *sx, slint s1elements) /* sl_proto, sl_func merge2_simplify_s1 */
{
  slint m, m0, m1;

  elements_t _s0, _s1, x;


  s1elements = z_min(s1elements, s1->size);

/*  printf("simplifying %d elements from s1\n", s1elements);*/

  if (s1elements == 0) return 0;

  /** find the s1elements highest elements of s0 and s1 **/

  m = s1elements;
  m0 = m1 = 0;

  elem_assign_at(s0, s0->size - 1, &_s0);
  elem_assign_at(s1, s1->size - 1, &_s1);

  while (m-- > 0 && _s0.keys >= s0->keys)
  if (key_cmp_ge(*_s0.keys, *_s1.keys))
  {
    elem_dec(&_s0);
    ++m0;
  } else
  {
    elem_dec(&_s1);
  }

  elem_inc(&_s0);

  m1 = s1elements - m0;
  elem_assign_at(s1, s1->size - m1, &_s1);

/*  printf("highest %d elements, %d from s0, %d from s1\n", s1elements, m0, m1);*/

  /** bring the s1elements highest elements to the end of s1 **/

  elem_assign_at(s1, s1->size - s1elements, &x);

  elem_nxchange(&_s0, &x, m0, sx);

  /* merge the highest elements */
  x.size = m0;
  _s1.size = m1;

/*  printf("merging x(%d) @ %p & _s1(%d) @ %p\n", x.size, x.keys, _s1.size, _s1.keys);*/

  the_merge2_01(&x, &_s1, sx);

/*  elements_print_keys(s1);*/

  /* merge the highest elements of s0 */
  elem_assign(s0, &x); x.size -= m0;

  _s0.size = m0;

/*  printf("merging x(%d) @ %p & _s0(%d) @ %p\n", x.size, x.keys, _s0.size, _s0.keys);*/

  the_merge2_01(&x, &_s0, sx);

/*  elements_print_keys(s0);*/

  s1->size -= s1elements;

  return 0;
}


#undef the_merge2_01


slint merge2_memory_adaptive(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_memory_adaptive */
{
  slint v;

  v = merge2_basic_auto_01_x(s0, s1, sx);

  if (v >= 0) return v;

  v = merge2_compo_tridgell(s0, s1, sx);

  if (v >= 0) return v;

  v = merge2_compo_hula(s0, s1, sx);
  
  if (v >= 0) return v;

  return merge2_compo_hula(s0, s1, NULL);
}


/*
   "Practical In-Place Merging" (adapted)
   Bing-Chao Huang, Michael A. Langston
   Communications of the ACM, March 1988, Volume 31, Number 3
*/


#include "sl_common.h"


#ifdef key_integer
# define sort_again(s, sx)    sort_radix(s, sx, -1, -1, -1)
#else
# define sort_again(s, sx)    sort_quick(s, sx)
#endif

#define the_merge2_01(s0, s1, sx)       merge2_basic_sbin_01(s0, s1, sx)
#define the_merge2_0X_1(s0, s1, sx, t)  merge2_basic_straight_0X_1(s0, s1, sx, t)
#define the_merge2_X0_1(s0, s1, sx, t)  merge2_basic_straight_X0_1(s0, s1, sx, t)


static slint_t merge2_compo_hula_(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_func merge2_compo_hula_ */
{
/*  slint_t i, j, k;*/

/*  slint_t n0, n1, t0, t1, e0, e1;*/

/*  elements_t _src0, _src1, _sx, undersized;*/
  elements_t wa;

  slint_t s, m, m0, m1, s0mod, u;

  elements_t _s0, _s1;

  elements_t sa, sb, sc, sd, sf, sg, x;

  elements_t first, second, last, current_head, current_tail, lowest_head, lowest_tail;

  s = (slint_t) slint_floor(slint_sqrt(s0->size + s1->size));

  if (s0->size < s || s1->size < s)
  {
    /* do an easier merge with bisection and rotates */
    the_merge2_01(s0, s1, sx);

/*    printf("easier merge2 done\n");*/

    return 0;
  }

  elem_assign(s0, &_s0);
  elem_assign(s1, &_s1);

  /* bring the size of s1 to a mulitple of s */
  merge2_simplify_s1(&_s0, &_s1, sx, s1->size % s);

/*  printf("_s0:\n"); elements_printf(&_s0); elements_print_keys(&_s0);
  printf("_s1:\n"); elements_printf(&_s1); elements_print_keys(&_s1);*/

  /* now: every list contains at least s elements and the size of s1 is a multiple of s */

  /** detect the s biggest elements **/
  elem_assign_at(&_s0, _s0.size - 1, &sa);
  elem_assign_at(&_s1, _s1.size - 1, &sb);

  m0 = m1 = 0;
  m = s;
  while (m-- > 0)
  if (key_cmp_ge(*sa.keys, *sb.keys))
  {
    elem_dec(&sa);
    ++m0;
  } else
  {
    elem_dec(&sb);
    ++m1;
  }

/*  printf("highest %d elements: %d from s0, %d from s1\n", s, m0, m1);*/

  /* now A & B are the workarea, C & D are the biggest elements */
  elem_inc(&sa); sa.size = m0;
  elem_inc(&sb); sb.size = m1;
  elem_assign_at(&_s0, _s0.size - s, &sc); sc.size = m1;
  elem_assign_at(&_s1, _s1.size - s, &sd); sd.size = m0;

/*  printf("sd:"); elements_printf(&sd);
  printf("sc:"); elements_printf(&sc);
  printf("sb:"); elements_printf(&sb);*/

  /* merge C and D forming block E at the end of s1 */
  the_merge2_0X_1(&sd, &sc, &sb, sx);

/*  printf("_s0:\n"); elements_printf(&_s0); elements_print_keys(&_s0);
  printf("_s1:\n"); elements_printf(&_s1); elements_print_keys(&_s1);*/

  /* after the merge, the workarea is in one piece at the end of s0 */
  elem_assign(&sc, &wa); wa.size = s;

/*  printf("wa:\n"); elements_printf(&wa); elements_print_keys(&wa);*/

  /* handle the undersized block at the beginning of s0 */
  s0mod = _s0.size % s;
  if (s0mod > 0)
  {
    elem_assign(&_s0, &sf); sf.size = s0mod; /* block F */
    elem_assign(&_s1, &sg); sg.size = s; /* block G */
    elem_assign_at(&wa, s - s0mod, &x); x.size = s0mod; /* ending part of the workarea */

/*    printf("sg:"); elements_printf(&sg);
    printf("sf:"); elements_printf(&sf);
    printf("x:"); elements_printf(&x);*/

    /* merge F & G forming H & I */
    the_merge2_X0_1(&sg, &sf, &x, sx);

/*    printf("_s0:\n"); elements_printf(&_s0); elements_print_keys(&_s0);
    printf("_s1:\n"); elements_printf(&_s1); elements_print_keys(&_s1);*/

    /* exchanging H (replaced F) and the part of the workarea at the front of s0 */
    elem_nxchange(&sf, &x, s0mod, sx);

/*    printf("_s0:\n"); elements_printf(&_s0); elements_print_keys(&_s0);
    printf("_s1:\n"); elements_printf(&_s1); elements_print_keys(&_s1);*/

    /* skip the finished undersized block at the front of s0 */
    elem_add(&_s0, s0mod); _s0.size -= s0mod;
  }

  /* bring the workarea to the front (behind an already finished undersized block) */
  elem_nxchange_ro0(&_s0, &wa, s, sx);
  elem_assign(&_s0, &wa); wa.size = s;

  /** now: all undersized blocks are handled, the sizes of _s0 and _s1 are multiples of s **/

/*  printf("merging elements of size %d @ %p and %d @ %p\n", _s0.size, _s0.keys, _s1.size, _s1.keys);
  printf("_s0:\n"); elements_printf(&_s0); elements_print_keys(&_s0);
  printf("_s1:\n"); elements_printf(&_s1); elements_print_keys(&_s1);*/


  /** sort the blocks **/

  /* the first block starts behind the workarea */
  elem_assign_at(&_s0, s, &first);
  elem_assign_at(&_s1, _s1.size - s, &last);

  /* while more the one block remains */
  while (first.keys < last.keys)
  {
    /* take the last block as lowest */
    elem_assign(&last, &lowest_head);
    elem_assign_at(&lowest_head, s - 1, &lowest_tail);

    /* start search at the first block */
    elem_assign(&first, &current_head);
    elem_assign_at(&current_head, s - 1, &current_tail);

/*    printf("starting, first @ %p\n", first.keys);
    printf("starting, lowest @ %p\n", lowest_head.keys);
    printf("starting, current @ [%p,%p]\n", current_head.keys, current_tail.keys);
    printf("starting, last @ %p\n", last.keys);*/

    while (current_tail.keys < last.keys)
    {
/*      printf("lowest: [%d,%d]\n", *lowest_head.keys, *lowest_tail.keys);
      printf("current: [%d,%d]\n", *current_head.keys, *current_tail.keys);*/

      if (key_cmp_lt(*current_tail.keys, *lowest_tail.keys) || (key_cmp_eq(*current_tail.keys, *lowest_tail.keys) && key_cmp_lt(*current_head.keys, *lowest_head.keys)))
      {
        elem_assign(&current_head, &lowest_head);
        elem_assign(&current_tail, &lowest_tail);

/*        printf("taken!\n");*/
      }

      elem_add(&current_head, s);
      elem_add(&current_tail, s);
    }

/*    printf("first %d @ %p\n", *first.keys, first.keys);
    printf("lowest %d @ %p\n", *lowest_head.keys, lowest_head.keys);*/

    /* bring the lowest block to the front (if not already there) */
    if (first.keys != lowest_head.keys) elem_nxchange(&first, &lowest_head, s, sx);

/*    printf("exchange done\n");*/

    /* continue sort on the remaining blocks */
    elem_add(&first, s);
  }

/*  printf("after blocksort\n");
  printf("_s0:\n"); elements_printf(&_s0); elements_print_keys(&_s0);
  printf("_s1:\n"); elements_printf(&_s1); elements_print_keys(&_s1);
  printf("wa:\n"); elements_printf(&wa); elements_print_keys(&wa);*/


  /** merge the blocks **/

  elem_assign_at(&_s0, s, &first); first.size = s;
  elem_assign_at(&_s1, _s1.size, &last);


  while (first.keys < last.keys)
  {
    elem_assign(&first, &lowest_head);
    elem_assign_at(&first, first.size - 1, &lowest_tail);

    elem_assign_at(&first, first.size, &second); second.size = s;

    while (second.keys < last.keys)
    if (key_cmp_le(*lowest_tail.keys, *second.keys))
    {
      first.size += s;
      elem_add(&lowest_tail, s);
      elem_add(&second, s);

    } else break;

    /* there is no more second block to merge with */
    if (second.keys >= last.keys) break;

/*    printf("merging blocks:\n");
    printf("first:\n"); elements_printf(&first); elements_print_keys(&first);
    printf("second:\n"); elements_printf(&second); elements_print_keys(&second);
    printf("wa:\n"); elements_printf(&wa); elements_print_keys(&wa);*/

    /* blockmerge */
    u = merge2_basic_straight_X0_1u(&first, &second, &wa, sx);

    elem_assign_at(&second, s - u, &first); first.size = u;
    elem_assign_at(&first, -s, &wa); wa.size = s;

/*    printf("after merging blocks:\n");
    printf("first:\n"); elements_printf(&first); elements_print_keys(&first);
    printf("wa:\n"); elements_printf(&wa); elements_print_keys(&wa);*/
  }

  /* bring the workarea to the end */
  elem_nxchange_ro0(&first, &wa, first.size, sx);
  elem_add(&wa, first.size);

/*  printf("starting final sort:\n");
  printf("wa:\n"); elements_printf(&wa); elements_print_keys(&wa);*/

  /* sort the workarea */
  sort_again(&wa, sx);

  return 0;
}


slint_t merge2_compo_hula(elements_t *s0, elements_t *s1, elements_t *xs) /* sl_proto, sl_func merge2_compo_hula */
{
  elements_t txs;

  if ((s0 == NULL) || (s1 == NULL)) return -1;

  /* if one location is empty, we are finished */
  if ((s0->size == 0) || (s1->size == 0)) return 0;

  /* calloc mode? */
  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;

    /* this is really "in-place", need exact one additional element of eXtraspace */
    elements_alloc(xs, 1, SLCM_ALL);
  }

  merge2_compo_hula_(s0, s1, xs);

  /* was in calloc mode? */
  if (xs == &txs) elements_free(xs);

  return 0;
}


#undef the_merge2_01
#undef the_merge2_0X_1
#undef the_merge2_X0_1



#include "sl_common.h"


#define the_ncopy      elem_ncopy
#define the_nmove      elem_nmove
#define the_rotate     elem_rotate


#define the_search_lt  sl_search_sequential_lt
#define the_search_le  sl_search_sequential_le
#define the_search_gt  sl_search_sequential_gt
#define the_search_ge  sl_search_sequential_ge


slint_t merge2_basic_sseq_x0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sseq_x0_1 */
{
  slint_t n;
  elements_t src0, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size, &src1e);
  elem_assign(sx, &dst);

  /* process until one of the srcs is empty */
  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_le(&src0, key_purify(*src1.keys));

    if (n > 0)
    {
      the_nmove(&src0, &dst, n);
      elem_add(&src0, n);
      elem_add(&dst, n);
      src0.size -= n;
    }

    elem_copy(&src1, &dst);
    elem_inc(&src1);
    elem_inc(&dst);
  }

  /* copy the remaining elements of s1 to dst */
  src1.size = src1e.keys - src1.keys;
  if (src1.size > 0) elem_ncopy(&src1, &dst, src1.size);

  return 0;
}


slint_t merge2_basic_sseq_0x_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sseq_0x_1 */
{
  slint_t n;
  elements_t src0, src0e, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign_at(s1, s1->size - 1, &src1); elem_assign_at(s1, -1, &src1e);
  elem_assign_at(sx, sx->size, &dst);

  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_gt(&src0, key_purify(*src1.keys));

    if (n > 0)
    {
      elem_sub(&dst, n);
      elem_sub(&src0e, n);
      the_nmove(&src0e, &dst, n);
      src0.size -= n;
    }

    elem_dec(&dst);
    elem_copy(&src1, &dst);
    elem_dec(&src1);
  }

  /* copy the remaining elements of s1 to the front */
  src1.size = src1.keys - src1e.keys;
  if (src1.size > 0) elem_ncopy(&src1, s0, src1.size);

  return 0;
}


slint_t merge2_basic_sseq_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sseq_01_x */
{
  return merge2_basic_01_x(s0, s1, sx, merge2_basic_sseq_x0_1, merge2_basic_sseq_0x_1);
}


slint_t merge2_basic_sseq_01(elements_t *s0, elements_t *s1, elements_t *t) /* sl_proto, sl_func merge2_basic_sseq_01 */
{
  slint_t k;
  elements_t x, _s0, _s1, last;

  if (t == NULL)
  {
    elements_alloc(&x, 1, SLCM_ALL);
    t = &x;
  }

  elem_assign(s0, &_s0);
  elem_assign(s1, &_s1);

  elem_assign_at(s1, s1->size - 1, &last);

  while (_s0.size > 0 && _s1.size > 0)
  if (_s0.size <= _s1.size)
  {
    k = the_search_lt(&_s1, key_purify(*_s0.keys));

    the_rotate(&_s0, _s0.size, k, t);

    elem_add(&_s0, k + 1);
    _s0.size -= 1;
    elem_add(&_s1, k);
    _s1.size -= k;

  } else
  {
    k = the_search_gt(&_s0, key_purify(*last.keys));

    elem_sub(&_s1, k);

    the_rotate(&_s1, k, _s1.size, t);

    elem_sub(&last, k + 1);
    _s0.size -= k;
    _s1.size -= 1;
  }

  if (t == &x) elements_free(&x);

  return 0;
}


#undef the_search_lt
#undef the_search_le
#undef the_search_gt
#undef the_search_ge


#define the_search_lt  sl_search_binary_lt
#define the_search_le  sl_search_binary_le
#define the_search_gt  sl_search_binary_gt
#define the_search_ge  sl_search_binary_ge


slint_t merge2_basic_sbin_x0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sbin_x0_1 */
{
  slint_t n;
  elements_t src0, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size, &src1e);
  elem_assign(sx, &dst);

  /* process until one of the srcs is empty */
  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_le(&src0, key_purify(*src1.keys));

    if (n > 0)
    {
      the_nmove(&src0, &dst, n);
      elem_add(&src0, n);
      elem_add(&dst, n);
      src0.size -= n;
    }

    elem_copy(&src1, &dst);
    elem_inc(&src1);
    elem_inc(&dst);
  }

  /* copy the remaining elements of s1 to dst */
  src1.size = src1e.keys - src1.keys;
  if (src1.size > 0) the_ncopy(&src1, &dst, src1.size);

  return 0;
}


slint_t merge2_basic_sbin_0x_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sbin_0x_1 */
{
  slint_t n;
  elements_t src0, src0e, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign_at(s1, s1->size - 1, &src1); elem_assign_at(s1, -1, &src1e);
  elem_assign_at(sx, sx->size, &dst);

  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_gt(&src0, key_purify(*src1.keys));

    if (n > 0)
    {
      elem_sub(&dst, n);
      elem_sub(&src0e, n);
      the_nmove(&src0e, &dst, n);
      src0.size -= n;
    }

    elem_dec(&dst);
    elem_copy(&src1, &dst);
    elem_dec(&src1);
  }

  /* copy the remaining elements of s1 to the front */
  src1.size = src1.keys - src1e.keys;
  elem_sub(&src1, src1.size - 1);
  if (src1.size > 0) the_ncopy(&src1, s0, src1.size);

  return 0;
}


slint_t merge2_basic_sbin_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sbin_01_x */
{
  return merge2_basic_01_x(s0, s1, sx, merge2_basic_sbin_x0_1, merge2_basic_sbin_0x_1);
}


slint_t merge2_basic_sbin_01(elements_t *s0, elements_t *s1, elements_t *t) /* sl_proto, sl_func merge2_basic_sbin_01 */
{
  slint_t k;
  elements_t x, _s0, _s1, last;

  if (t == NULL)
  {
    elements_alloc(&x, 1, SLCM_ALL);
    t = &x;
  }

  elem_assign(s0, &_s0);
  elem_assign(s1, &_s1);

  elem_assign_at(s1, s1->size - 1, &last);

  while (_s0.size > 0 && _s1.size > 0)
  if (_s0.size <= _s1.size)
  {
    k = the_search_lt(&_s1, key_purify(*_s0.keys));

    the_rotate(&_s0, _s0.size, k, t);

    elem_add(&_s0, k + 1);
    _s0.size -= 1;
    elem_add(&_s1, k);
    _s1.size -= k;

  } else
  {
    k = the_search_gt(&_s0, key_purify(*last.keys));

    elem_sub(&_s1, k);

    the_rotate(&_s1, k, _s1.size, t);

    elem_sub(&last, k + 1);
    _s0.size -= k;
    _s1.size -= 1;
  }

  if (t == &x) elements_free(&x);

  return 0;
}


#undef the_search_lt
#undef the_search_le
#undef the_search_gt
#undef the_search_ge


/* using "hybrid-search" for binary-merge from [Hwang, Lin] (Knuth v3) */

#define the_search_lt  sl_search_hybrid_lt
#define the_search_le  sl_search_hybrid_le
#define the_search_gt  sl_search_hybrid_gt
#define the_search_ge  sl_search_hybrid_ge


slint_t merge2_basic_shyb_x0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_shyb_x0_1 */
{
  slint_t n, l;
  elements_t src0, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size, &src1e);
  elem_assign(sx, &dst);

  l = z_powof2((slint_t) (log((double) s0->size / (double) s1->size) / log(2.0)));

  /* process until one of the srcs is empty */
  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_le(&src0, key_purify(*src1.keys), l);

    if (n > 0)
    {
      the_nmove(&src0, &dst, n);
      elem_add(&src0, n);
      elem_add(&dst, n);
      src0.size -= n;
    }

    elem_copy(&src1, &dst);
    elem_inc(&src1);
    elem_inc(&dst);
  }

  /* copy the remaining elements of s1 to dst */
  src1.size = src1e.keys - src1.keys;
  if (src1.size > 0) elem_ncopy(&src1, &dst, src1.size);

  return 0;
}


slint_t merge2_basic_shyb_0x_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_shyb_0x_1 */
{
  slint_t n, l;
  elements_t src0, src0e, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign_at(s1, s1->size - 1, &src1); elem_assign_at(s1, -1, &src1e);
  elem_assign_at(sx, sx->size, &dst);

  l = z_powof2((slint_t) (log((double) s0->size / (double) s1->size) / log(2.0)));

  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_gt(&src0, key_purify(*src1.keys), l);

    if (n > 0)
    {
      elem_sub(&dst, n);
      elem_sub(&src0e, n);
      the_nmove(&src0e, &dst, n);
      src0.size -= n;
    }

    elem_dec(&dst);
    elem_copy(&src1, &dst);
    elem_dec(&src1);
  }

  /* copy the remaining elements of s1 to the front */
  src1.size = src1.keys - src1e.keys;
  if (src1.size > 0) elem_ncopy(&src1, s0, src1.size);

  return 0;
}


slint_t merge2_basic_shyb_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_shyb_01_x */
{
  return merge2_basic_01_x(s0, s1, sx, merge2_basic_shyb_x0_1, merge2_basic_shyb_0x_1);
}


slint_t merge2_basic_shyb_01(elements_t *s0, elements_t *s1, elements_t *t) /* sl_proto, sl_func merge2_basic_shyb_01 */
{
  slint_t k, l;
  elements_t x, _s0, _s1, last;

  if (t == NULL)
  {
    elements_alloc(&x, 1, SLCM_ALL);
    t = &x;
  }

  elem_assign(s0, &_s0);
  elem_assign(s1, &_s1);

  elem_assign_at(s1, s1->size - 1, &last);

  while (_s0.size > 0 && _s1.size > 0)
  if (_s0.size <= _s1.size)
  {
    l = z_powof2((slint_t) (log((double) _s1.size / (double) _s0.size) / log(2.0)));

    k = the_search_lt(&_s1, key_purify(*_s0.keys), l);

    the_rotate(&_s0, _s0.size, k, t);

    elem_add(&_s0, k + 1);
    _s0.size -= 1;
    elem_add(&_s1, k);
    _s1.size -= k;

  } else
  {
    l = z_powof2((slint_t) (log((double) _s0.size / (double) _s1.size) / log(2.0)));

    k = the_search_gt(&_s0, key_purify(*last.keys), l);

    elem_sub(&_s1, k);

    the_rotate(&_s1, k, _s1.size, t);

    elem_sub(&last, k + 1);
    _s0.size -= k;
    _s1.size -= 1;
  }

  if (t == &x) elements_free(&x);

  return 0;
}


#undef the_search_lt
#undef the_search_le
#undef the_search_gt
#undef the_search_ge


#undef the_ncopy
#undef the_nmove
#undef the_rotate



#include "sl_common.h"


/* front2back, overriding the content of the extraspace */
slint_t merge2_basic_straight_x0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_straight_x0_1 */
{
  elements_t src0, src0e, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size, &src1e);
  elem_assign(sx, &dst);

  /* process until one of the srcs is empty */
  while (src0.keys != src0e.keys && src1.keys != src1e.keys)
  {
    if (key_cmp_le(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_inc(&src0);
    } else
    {
      elem_copy(&src1, &dst);
      elem_inc(&src1);
    }
    elem_inc(&dst);
  }

  /* copy the remaining elements of the 1st src */
  src1.size = src1e.keys - src1.keys;
  if (src1.size > 0) elem_ncopy(&src1, &dst, src1.size);

  return 0;
}


/* back2front, overriding the content of the extraspace */
slint_t merge2_basic_straight_0x_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_straight_0x_1 */
{
  elements_t src0, src0e, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign_at(s0, s0->size - 1, &src0); elem_assign_at(s0, -1, &src0e);
  elem_assign_at(s1, s1->size - 1, &src1); elem_assign_at(s1, -1, &src1e);
  elem_assign_at(sx, s1->size - 1, &dst);

  /* process until one of the srcs is empty */
  while (src0.keys != src0e.keys && src1.keys != src1e.keys)
  {
    if (key_cmp_gt(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_dec(&src0);
    } else
    {
      elem_copy(&src1, &dst);
      elem_dec(&src1);
    }
    elem_dec(&dst);
  }

  /* copy the remaining elements of the 1st src */
  src1.size = src1.keys - src1e.keys;
  if (src1.size > 0) elem_ncopy(s1, s0, src1.size);

  return 0;
}


slint_t merge2_basic_straight_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_straight_01_x */
{
  return merge2_basic_01_x(s0, s1, sx, merge2_basic_straight_x0_1, merge2_basic_straight_0x_1);
}


slint_t merge2_basic_straight_x_0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_straight_x_0_1 */
{
  elements_t dst, dste, src0e, src1e;

  elem_assign_at(s0, s0->size, &src0e);
  elem_assign_at(s1, s1->size, &src1e);
  elem_assign(sx, &dst); elem_assign_at(sx, sx->size, &dste);

  while (s0->keys != src0e.keys && s1->keys != src1e.keys && dst.keys != dste.keys)
  {
    if (key_cmp_le(*s0->keys, *s1->keys))
    {
      elem_copy(s0, &dst);
      elem_inc(s0);

    } else
    {
      elem_copy(s1, &dst);
      elem_inc(s1);
    }
    elem_inc(&dst);
  }

  s0->size = src0e.keys - s0->keys;
  s1->size = src1e.keys - s1->keys;

  return dst.keys - sx->keys;
}


slint_t merge2_basic_straight_X0_1(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t) /* sl_proto, sl_func merge2_basic_straight_X0_1 */
{
  elements_t src0, src0e, src1, src1l, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size - 1, &src1l);
  elem_assign(sx, &dst);

  /* create the hole */
  elem_copy(&dst, t);

  while (src0.keys != src0e.keys && src1.keys != src1l.keys)
  {
    if (key_cmp_le(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_inc(&dst);
      elem_copy(&dst, &src0);
      elem_inc(&src0);

    } else
    {
      elem_copy(&src1, &dst);
      elem_inc(&dst);
      elem_copy(&dst, &src1);
      elem_inc(&src1);
    }
  }

  /* now: either src0 is empty or src1 contains just one element */

  while (src0.keys != src0e.keys)
  if (key_cmp_le(*src0.keys, *src1.keys))
  {
    /* remove one element from src0 */
    elem_copy(&src0, &dst);
    elem_inc(&dst);
    elem_inc(&src0);

  } else break; /* the last of src1 has to move */

  /* close the hole at dst */
  elem_copy(t, &dst);

  /* exchange the remaining elements at src1 with dst */
  src1.size = src1l.keys - src1.keys + 1;
  elem_nxchange_ro0(&src1, &dst, src1.size, t);

  return 0;
}


slint_t merge2_basic_straight_0X_1(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t) /* sl_proto, sl_func merge2_basic_straight_0X_1 */
{
  elements_t src0, src1, src0e, src1l, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign_at(s0, s0->size - 1, &src0); elem_assign_at(s0, -1, &src0e);
  elem_assign_at(s1, s1->size - 1, &src1); elem_assign(s1, &src1l);
  elem_assign_at(sx, s1->size - 1, &dst);

  /* create the hole */
  elem_copy(&dst, t);

  while (src0.keys != src0e.keys && src1.keys != src1l.keys)
  {
    if (key_cmp_gt(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_dec(&dst);
      elem_copy(&dst, &src0);
      elem_dec(&src0);

    } else
    {
      elem_copy(&src1, &dst);
      elem_dec(&dst);
      elem_copy(&dst, &src1);
      elem_dec(&src1);
    }
  }

  /* now: either src0 is empty or src1 contains just one element */

  while (src0.keys != src0e.keys)
  if (key_cmp_gt(*src0.keys, *src1.keys))
  {
    /* remove one (maybe the last) from src0 */
    elem_copy(&src0, &dst);
    elem_dec(&dst);
    elem_dec(&src0);

  } else break; /* the last of src1 has to move */

  /* close the hole at dst */
  elem_copy(t, &dst);

  /* xchange the remaining elements at s1 with s0 */
  src1.size = src1.keys - src1l.keys + 1;
  elem_inc(&src0);
  elem_nxchange_ro0(s1, &src0, src1.size, t);

  return 0;
}


slint_t merge2_basic_straight_01_X(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t) /* sl_proto, sl_func merge2_basic_straight_01_X */
{
  return merge2_basic_01_X(s0, s1, sx, t, merge2_basic_straight_X0_1, merge2_basic_straight_0X_1);
}


slint_t merge2_basic_straight_X0_1u(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t) /* sl_proto, sl_func merge2_basic_straight_X0_1u */
{
  elements_t src0, src0e, src1, src1l, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size - 1, &src1l);
  elem_assign(sx, &dst);

  /* create the hole */
  elem_copy(&dst, t);

  while (src0.keys != src0e.keys && src1.keys != src1l.keys)
  {
    if (key_cmp_le(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_inc(&dst);
      elem_copy(&dst, &src0);
      elem_inc(&src0);

    } else
    {
      elem_copy(&src1, &dst);
      elem_inc(&dst);
      elem_copy(&dst, &src1);
      elem_inc(&src1);
    }
  }

  /* now: either src0 is empty or src1 contains just one element */

  while (src0.keys != src0e.keys)
  if (key_cmp_le(*src0.keys, *src1.keys))
  {
    /* remove one element from src0 */
    elem_copy(&src0, &dst);
    elem_inc(&dst);
    elem_inc(&src0);

  } else
  {
    /* the last of src1 has to move */
    elem_copy(&src1, &dst);

    /* let the hole (dst) point to the last element of s1 */
    elem_assign(&src1, &dst);

    /* let src1 point behind src1l for computing the right number of unfinished elements of s1 (0!) */
    elem_inc(&src1);

    break;
  }

  /* close the hole at dst */
  elem_copy(t, &dst);

  /* return the number of unfinished elements of s1 */
  return src1l.keys - src1.keys + 1;
}


/* slightly optimized (performance and complexity of implementation?) merge2 from [Tridgell] */


#include "sl_common.h"


#define the_merge_x_0_1  merge2_basic_straight_x_0_1


typedef struct _block
{
  slint_t rank, size, quan;

  elements_t s;

} block;


slint_t merge2_compo_tridgell(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_compo_tridgell */
{
  slint_t blocksize = (slint_t) slint_sqrt(s0->size + s1->size);

  slint_t s0blocks = (s0->size + blocksize - 1) / blocksize;
  slint_t s1blocks = (s1->size + blocksize - 1) / blocksize;

  slint_t s0mod = s0->size % blocksize;
  slint_t s1offset = (blocksize - s0mod) % blocksize;

  slint_t s0rank = 0;
  slint_t s1rank = s0blocks;

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  block *blocks = NULL;
#else
  block blocks[s0blocks + s1blocks + 3];
#endif

  elements_t src0, src0b, src1, src1b, dst, _sx;

  slint_t nblocks, cblock, the_empty_blocks[3], empty_blocks;

  slint_t i, n, t, s;


/*  printf("jetzt gehts los!\n");
  printf("blocksize = %d\n", blocksize);
  printf("s0size = %d, s0blocks = %d, s0mod = %d\n", s0->size, s0blocks, s0mod);
  printf("s1size = %d, s1blocks = %d, s1offset = %d\n", s1->size, s1blocks, s1offset);*/

  /* if one list is empty, there is nothing to do */
  if (s0->size == 0 || s1->size == 0) return 0;

  if (sx == NULL)
  {
    sx = &_sx;
    elements_alloc(sx, 3 * blocksize, SLCM_ALL);
  }

  /* if the supplied eXtra-space is less than 3*blocksize */
  if (sx->size < 3 * blocksize) return -1;

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  blocks = z_alloc(s0blocks + s1blocks + 3, sizeof(block));
#endif

  elem_assign(s0, &src0);
  elem_assign(s1, &src1);

/*  printf("src0.size = %d - src1.size = %d\n", src0.size, src1.size);*/

  /* skip the first 's1offset' elements in s1 when building free blocks */
  elem_assign_at(s0, 0, &src0b);
  elem_assign_at(s1, s1offset, &src1b);

  /* first free block can be undersized */
  blocks[0].rank = -3;
  elem_assign_at(sx, 0 * blocksize, &blocks[0].s);
  blocks[0].s.size = blocksize;
  blocks[0].size = 0;

  blocks[1].rank = -2;
  elem_assign_at(sx, 1 * blocksize, &blocks[1].s);
  blocks[1].s.size = blocksize;
  blocks[1].size = 0;

  blocks[2].rank = -1;
  elem_assign_at(sx, 2 * blocksize, &blocks[2].s);
  blocks[2].s.size = blocksize;
  blocks[2].size = 0;

  cblock = -1;
  nblocks = 3;

  /* merge form s0 and s1 to the destination blocks */
  while (src0.size > 0 && src1.size > 0)
  {
    /* take new free block */
    cblock++;

/*    printf("src0.size = %d - src1.size = %d\n", src0.size, src1.size);
    printf("filling %dth block @ %p\n", cblock, blocks[cblock].s.keys);*/

    /* merge to the new free block */
    blocks[cblock].size = the_merge_x_0_1(&src0, &src1, &blocks[cblock].s);

/*    printf("one merge done, |src0| = %d, |src1| = %d\n", src0.size, src1.size);
    printf("free: %d & %d\n", src0.keys - src0b.keys, src1.keys - src1b.keys);*/

    /* append a new free and complete block from s0 */
    if (src0.keys - src0b.keys >= blocksize)
    {
/*      printf("adding free #%d from s0 @ %p\n", nblocks, src0b.keys);*/

      elem_assign(&src0b, &blocks[nblocks].s);
      blocks[nblocks].s.size = blocksize;
      blocks[nblocks].rank = s0rank++;
      blocks[nblocks].size = 0;

      elem_add(&src0b, blocksize);

/*      printf("new from s0 may start @ %p\n", src0b.keys);*/

      nblocks++;
    }

    /* append a new free and complete block from s1 */
    if (src1.keys - src1b.keys >= blocksize)
    {
/*      printf("adding free #%d from s1 @ %p\n", nblocks, src1b.keys);*/

      elem_assign(&src1b, &blocks[nblocks].s);
      blocks[nblocks].s.size = blocksize;
      blocks[nblocks].rank = s1rank++;
      blocks[nblocks].size = 0;

      elem_add(&src1b, blocksize);

      nblocks++;
    }
  }

/*  for (i = 0; i < nblocks; ++i)
  {
    printf("block %d @ %p: rank = %d, content: %d of %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size);
    printf("block %d: content = %d\n", i, blocks[i].size);
    elements_print_keys(&blocks[i].s);
  }*/

/*  printf("one is empty, |src0| = %d, |src1| = %d\n", src0.size, src1.size);*/

/*  printf("elements left: s0 -> %d, s1 -> %d\n", src0.size, src1.size);
  for (i = 0; i < nblocks; ++i) printf("block %d @ %p: rank = %d, content: %d of %d, first = %d, last = %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size, blocks[i].s.keys[0], blocks[i].s.keys[blocks[i].size - 1]);*/

  /* prepare the last half-filled block */
  elem_assign_at(&blocks[cblock].s, blocks[cblock].size, &dst);
  dst.size = blocks[cblock].s.size - blocks[cblock].size;

/*  printf("cblock %d (%d,%d) needs %d\n", cblock, blocks[cblock].s.size, blocks[cblock].size, dst.size);*/

  /* copy all the remaining elements from s0 to the destination blocks */
  while (src0.size > 0)
  {
    n = z_min(src0.size, dst.size);

/*    printf("filling %dth block @ %p with %d elements\n", cblock, blocks[cblock].s.keys, n);*/

    /* fill up with elements from s0 */
    elem_ncopy(&src0, &dst, n);
    elem_add(&src0, n);
    src0.size -= n;
    blocks[cblock].size += n;

/*    printf("copy done, |src0| = %d, |src1| = %d\n", src0.size, src1.size);
    printf("free: %d & %d\n", src0.keys - src0b.keys, src1.keys - src1b.keys);*/

    /* append a new free and full-sized block from s0 */
    if (src0.keys - src0b.keys >= blocksize)
    {
/*      printf("adding free #%d from s0 @ %p\n", nblocks, src0b.keys);*/

      elem_assign(&src0b, &blocks[nblocks].s);
      blocks[nblocks].s.size = blocksize;
      blocks[nblocks].rank = s0rank++;
      blocks[nblocks].size = 0;

      elem_add(&src0b, blocksize);

/*      printf("setting src0b to %p\n", src0b.keys);*/

      nblocks++;
    }

    /* take a new free block */
    cblock++;

    elem_assign(&blocks[cblock].s, &dst);
  }

  /* add an undersized free block from the end of s1, if possible */
  if (src1.keys - src1b.keys > 0)
  {
/*    printf("adding free #%d from s1 @ %p size %d\n", nblocks, src1b.keys, src1.keys - src1b.keys);*/

    elem_assign(&src1b, &blocks[nblocks].s);
    blocks[nblocks].s.size = src1.keys - src1b.keys;
    blocks[nblocks].rank = s1rank++;
    blocks[nblocks].size = 0;

    nblocks++;
  }

  empty_blocks = 0;
/*  printf("elements left: s0 -> %d, s1 -> %d\n", src0.size, src1.size);*/

/*  printf("elements left: s0 -> %d, s1 -> %d\n", src0.size, src1.size);
  for (i = 0; i < nblocks; ++i) printf("block %d @ %p: rank = %d, content: %d of %d, first = %d, last = %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size, blocks[i].s.keys[0], blocks[i].s.keys[blocks[i].size - 1]);*/

  /* if there is a (broken) block at the end of s0 and the beginning of s1 */
  if (s1offset > 0)
  {
/*    printf("handle broken block #%d of size %d, break into %d:%d, from %p to %p\n", s0blocks - 1, blocks[s0blocks - 1].size, s0mod, z_min(s1offset, src1.keys - s1->keys), blocks[s0blocks - 1].s.keys, src0b.keys);*/

    /* bring this bock in the right position */
    elem_ncopy(&blocks[s0blocks - 1].s, &src0b, s0mod);
    elem_assign_at(&blocks[s0blocks - 1].s, s0mod, &dst);
    elem_ncopy(&dst, s1, z_min(s1offset, src1.keys - s1->keys));

    /* mark the source of this block as empty */
    blocks[s0blocks - 1].size = 0;
  }

/*  printf("elements left: s0 -> %d, s1 -> %d\n", src0.size, src1.size);
  for (i = 0; i < nblocks; ++i) printf("block %d @ %p: rank = %d, content: %d of %d, first = %d, last = %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size, blocks[i].s.keys[0], blocks[i].s.keys[blocks[i].size - 1]);*/

  /* notice the possibly empty last 3 blocks */
  if (nblocks > 0)
  if (blocks[nblocks - 1].size == 0) the_empty_blocks[empty_blocks++] = nblocks - 1;
  if (nblocks > 1)
  if (blocks[nblocks - 2].size == 0) the_empty_blocks[empty_blocks++] = nblocks - 2;
  if (nblocks > 2)
  if (blocks[nblocks - 3].size == 0) the_empty_blocks[empty_blocks++] = nblocks - 3;

  /* notice the empty block (created by handling a broken block at the end of s0 and the beginning of s1),
     BUT only when he wasn't one of the former 3 empty block (and therefore already mentioned) */
  if (blocks[s0blocks - 1].size == 0 && (s0blocks - 1 < nblocks - 3)) the_empty_blocks[empty_blocks++] = s0blocks - 1;

/*  printf("empty_blocks at the end: %d -> %d, %d, %d\n", empty_blocks, the_empty_blocks[0], the_empty_blocks[1], the_empty_blocks[2]);

  for (i = 0; i < nblocks; ++i)
  {
    printf("block %d @ %p: rank = %d, content: %d of %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size);
    printf("block %d: content = %d\n", i, blocks[i].size);
    elements_print_keys(&blocks[i].s);
  }*/

  /* bring all blocks to their right positions */

  /* 1st, start cycling at empty blocks, this (max 3) cycles ending in the extraspace (sx) */
  for (i = 0; i < empty_blocks; ++i)
  {
    t = the_empty_blocks[i];
    s = blocks[t].rank;

/*    printf("empty cycle %d starting at %d, size = %d\n", i, the_empty_blocks[i], blocks[the_empty_blocks[i]].size);*/

    while (s >= 0)
    {
/*      printf("copying %d elements from block #%d to #%d\n", blocks[t].s.size, s, t);*/

      elem_ncopy(&blocks[s].s, &blocks[t].s, blocks[t].s.size);
      blocks[t].size = -1; /* mark the target block as 'done' */
      blocks[s].size = 0; /* mark the source block as 'empty' */

/*      printf("source: %d (%d), target: %d (%d), size = %d, first = %d, last = %d\n", s, blocks[s].rank, t, blocks[t].rank, blocks[t].s.size, blocks[t].s.keys[0], blocks[t].s.keys[blocks[t].s.size - 1]);*/

      t = s;
      s = blocks[t].rank;
    }
  }

/*  for (i = 0; i < nblocks; ++i) printf("block %d @ %p: rank = %d, content: %d of %d, first = %d, last = %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size, blocks[i].s.keys[0], blocks[i].s.keys[blocks[i].size - 1]);*/

  /* now: the 3 blocks residing in the extraspace (sx) are at there correct positions -> extraspace (sx) is free for use */

  /* 2nd, start cycling at full blocks */
  for (i = nblocks - 1; i > 2; --i)
  if (blocks[i].size > 0) /* skip already well-placed blocks (marked as 'done' -> size == -1) */
  {
/*    printf("cycle starting at %d, size = %d\n", i, blocks[i].size);*/

    /* copy the content of the current block to an empty location, if necessary */
    elem_ncopy(&blocks[i].s, sx, blocks[i].size);

    t = i;
    s = blocks[t].rank;

    /* cycle until we stuck at the current block ('i' is where we started) */
    while (s != i)
    {
/*      printf("copying %d elements from block #%d to #%d\n", blocks[t].s.size, s, t);*/

      elem_ncopy(&blocks[s].s, &blocks[t].s, blocks[t].s.size);
      blocks[t].size = -1; /* mark the target block as 'done' */
      blocks[s].size = 0; /* mark the source block as 'empty' */

/*      printf("source: %d (%d), target: %d (%d), size = %d, first = %d, last = %d\n", s, blocks[s].rank, t, blocks[t].rank, blocks[t].s.size, blocks[t].s.keys[0], blocks[t].s.keys[blocks[t].s.size - 1]);*/

      t = s;
      s = blocks[t].rank;
    }

    /* finished this cycle by copying from the empty location */
    elem_ncopy(sx, &blocks[t].s, blocks[t].s.size);
    blocks[t].size = -1;
  }
  
#ifdef NO_VARIABLE_LENGTH_ARRAYS
  z_free(blocks);
#endif

  if (sx == &_sx) elements_free(sx);

  return 0;
}


#undef the_merge_x_0_1



/* sl_macro MP2W_TRACE_IF */


#include "sl_common.h"


#ifndef MP2W_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MP2W_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MP2W_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mergep_2way_ip_int(elements_t *s, elements_t *sx, slint_t p, int *displs, merge2x_f m2x) /* sl_proto, sl_func mergep_2way_ip_int */
{
  slint_t i, step, counts[p];
  elements_t s0, s1;


  Z_TRACE_IF(MP2W_TRACE_IF, "merging %" slint_fmt " sub-sequences", p);

  for (i = 0; i < p - 1; ++i) counts[i] = displs[i + 1] - displs[i];
  counts[p - 1] = s->size - displs[p - 1];

  Z_TRACE_ARRAY_IF(MP2W_TRACE_IF, i, p, " %d", displs[i], "displs =");
  Z_TRACE_ARRAY_IF(MP2W_TRACE_IF, i, p, " %" slint_fmt, counts[i], "counts =");
  
  for (step = 1; step < p; step *= 2)
  {
    for (i = 0; i < p; i += 2 * step)
    {
      if (i + step < p)
      {
        elem_assign_at(s, displs[i], &s0);
        s0.size = counts[i];
        elem_assign_at(s, displs[i + step], &s1);
        s1.size = counts[i + step];
        
        Z_TRACE_IF(MP2W_TRACE_IF, "merging %" slint_fmt " and %" slint_fmt, s0.size, s1.size);

        if (s0.size > 0 && s1.size > 0) m2x(&s0, &s1, sx);
        
        counts[i] = s0.size + s1.size;
      }
    }
  }
  
  return 0;
}


static inline void rec(elements_t *s, elements_t *sx, slint_t p, int offset, int *displs, merge2x_f m2x)
{
  slint_t p0, p1;
  elements_t s0, s1;

  p0 = p / 2;
  p1 = p - p0;

  elem_assign(s, &s0); s0.size = displs[p0] - offset;
  if (p0 > 1) rec(&s0, sx, p0, offset, displs, m2x);

  elem_assign_at(s, displs[p0] - offset, &s1); s1.size = s->size - (displs[p0] - offset);
  if (p1 > 1) rec(&s1, sx, p1, displs[p0], displs + p0, m2x);

  m2x(&s0, &s1, sx);
}


slint_t mergep_2way_ip_int_rec(elements_t *s, elements_t *sx, slint_t p, int *displs, merge2x_f m2x) /* sl_proto, sl_func mergep_2way_ip_int_rec */
{
  if (p > 1) rec(s, sx, p, 0, displs, m2x);
  
  return 0;
}



#include "sl_common.h"


/* - heap gleich im dst-Feld -> heap umdrehen (Min bewegt sich aufs Ende zu)?
*/


#define HEAPIFY_INTRO() \
  hkx = heap_keys[j]; \
  hsx = heap_sources[j];

#define HEAPIFY_EXTRO() \
  heap_keys[j] = hkx; \
  heap_sources[j] = hsx;

#define HEAPIFY_UP() \
  while (j > 0) \
  { \
    k = (j - 1) / 2; \
    if (key_pure_cmp_le(heap_keys[k], hkx)) break; \
    heap_keys[j] = heap_keys[k]; \
    heap_sources[j] = heap_sources[k]; \
    j = k; \
  }

#define HEAPIFY_DOWN() \
  while ((k = 2 * j + 1) < heap_size) \
  { \
    if (k + 1 < heap_size) \
    if (key_pure_cmp_gt(heap_keys[k], heap_keys[k + 1])) ++k; \
    if (key_pure_cmp_gt(heap_keys[k], hkx)) break; \
    heap_keys[j] = heap_keys[k]; \
    heap_sources[j] = heap_sources[k]; \
    j = k; \
  }


slint_t mergep_heap_int(elements_t *s, elements_t *d, slint_t p, int *displs, int *counts) /* sl_proto, sl_func mergep_heap_int */
{
  slkey_pure_t heap_keys[p], hkx;
  slint_t heap_sources[p], hsx;
  slint_t heap_size = 0;

  slint_t i, j, k;
  
  int local_displs[p], local_counts[p];
  
  elements_t dst;


  memcpy(local_displs, displs, p * sizeof(int));

  if (counts == NULL)
  {
    for (i = 0; i < p - 1; ++i) local_counts[i] = local_displs[i + 1] - local_displs[i];
    local_counts[p - 1] = s->size - local_displs[p - 1];

  } else memcpy(local_counts, counts, p * sizeof(int));

  for (i = 0; i < p; ++i)
  if (local_counts[i] > 0)
  {
    heap_keys[heap_size] = key_purify(*elem_key_at(s, local_displs[i]));
    heap_sources[heap_size] = i;

    j = heap_size;
    ++heap_size;

    HEAPIFY_INTRO();
    HEAPIFY_UP();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  elem_assign(d, &dst);

  while (heap_size > 0)
  {
    i = heap_sources[0];

    /* copy min element */
    elem_copy_at(s, local_displs[i], &dst, 0);

    elem_inc(&dst);

    --local_counts[i];
    ++local_displs[i];
    
    if (local_counts[i] > 0) heap_keys[0] = key_purify(*elem_key_at(s, local_displs[i]));
    else 
    {
      heap_keys[0] = heap_keys[heap_size - 1];
      heap_sources[0] = heap_sources[heap_size - 1];
      --heap_size;
    }

    j = 0;
    HEAPIFY_INTRO();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  return 0;
}


slint_t mergep_heap_int_idx(elements_t *s, elements_t *d, slint_t p, int *displs, int *counts) /* sl_proto, sl_func mergep_heap_int_idx */
{
  slkey_pure_t heap_keys[p], hkx;
  slint_t heap_sources[p], hsx;
  slint_t heap_size = 0;

  slint_t i, j, k;
  
  int local_displs[p], local_counts[p];
  
  elements_t dst;


  memcpy(local_displs, displs, p * sizeof(int));

  if (counts == NULL)
  {
    for (i = 0; i < p - 1; ++i) local_counts[i] = local_displs[i + 1] - local_displs[i];
    local_counts[p - 1] = s->size - local_displs[p - 1];

  } else memcpy(local_counts, counts, p * sizeof(int));

  for (i = 0; i < p; ++i)
  if (local_counts[i] > 0)
  {
    heap_keys[heap_size] = key_purify(*elem_key_at(s, local_displs[i]));
    heap_sources[heap_size] = i;

    j = heap_size;
    ++heap_size;

    HEAPIFY_INTRO();
    HEAPIFY_UP();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  elem_assign(d, &dst);

  while (heap_size > 0)
  {
    i = heap_sources[0];

    /* copy min element */
    elem_copy_at(s, local_displs[i], &dst, 0);

#ifdef SL_INDEX
    *dst.indices = local_displs[i];
#endif

    elem_inc(&dst);

    --local_counts[i];
    ++local_displs[i];
    
    if (local_counts[i] > 0) heap_keys[0] = key_purify(*elem_key_at(s, local_displs[i]));
    else 
    {
      heap_keys[0] = heap_keys[heap_size - 1];
      heap_sources[0] = heap_sources[heap_size - 1];
      --heap_size;
    }

    j = 0;
    HEAPIFY_INTRO();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  return 0;
}


slint_t mergep_heap_idx(elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts) /* sl_proto, sl_func mergep_heap_idx */
{
  slkey_pure_t heap_keys[p], hkx;
  slint_t heap_sources[p], hsx;
  slint_t heap_size = 0;

  slint_t i, j, k;
  
  slindex_t local_displs[p], local_counts[p];
  
  elements_t dst;


  memcpy(local_displs, displs, p * sizeof(slindex_t));

  if (counts == NULL)
  {
    for (i = 0; i < p - 1; ++i) local_counts[i] = local_displs[i + 1] - local_displs[i];
    local_counts[p - 1] = s->size - local_displs[p - 1];

  } else memcpy(local_counts, counts, p * sizeof(slindex_t));

  for (i = 0; i < p; ++i)
  if (local_counts[i] > 0)
  {
    heap_keys[heap_size] = key_purify(*elem_key_at(s, local_displs[i]));
    heap_sources[heap_size] = i;

    j = heap_size;
    ++heap_size;

    HEAPIFY_INTRO();
    HEAPIFY_UP();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  elem_assign(d, &dst);

  while (heap_size > 0)
  {
    i = heap_sources[0];

    /* copy min element */
    elem_copy_at(s, local_displs[i], &dst, 0);
    
#ifdef SL_INDEX
    *dst.indices = local_displs[i];
#endif

    elem_inc(&dst);

    --local_counts[i];
    ++local_displs[i];
    
    if (local_counts[i] > 0) heap_keys[0] = key_purify(*elem_key_at(s, local_displs[i]));
    else 
    {
      heap_keys[0] = heap_keys[heap_size - 1];
      heap_sources[0] = heap_sources[heap_size - 1];
      --heap_size;
    }

    j = 0;
    HEAPIFY_INTRO();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  return 0;
}


slint_t mergep_heap_unpack_idx(packed_elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts) /* sl_proto, sl_func mergep_heap_unpack_idx */
{
  slkey_pure_t heap_keys[p], hkx;
  slint_t heap_sources[p], hsx;
  slint_t heap_size = 0;

  slint_t i, j, k;
  
  slindex_t local_displs[p], local_counts[p];
  
  elements_t dst;


  memcpy(local_displs, displs, p * sizeof(slindex_t));

  if (counts == NULL)
  {
    for (i = 0; i < p - 1; ++i) local_counts[i] = local_displs[i + 1] - local_displs[i];
    local_counts[p - 1] = s->size - local_displs[p - 1];

  } else memcpy(local_counts, counts, p * sizeof(slindex_t));

  for (i = 0; i < p; ++i)
  if (local_counts[i] > 0)
  {
    heap_keys[heap_size] = key_purify(*pelem_key_at(s, local_displs[i]));
    heap_sources[heap_size] = i;

    j = heap_size;
    ++heap_size;

    HEAPIFY_INTRO();
    HEAPIFY_UP();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  elem_assign(d, &dst);

  while (heap_size > 0)
  {
    i = heap_sources[0];

    /* copy min element */
    pelem_unpack_at(s, local_displs[i], &dst, 0);
    
#ifdef SL_INDEX 
    *dst.indices = local_displs[i];
#endif

    elem_inc(&dst);

    --local_counts[i];
    ++local_displs[i];
    
    if (local_counts[i] > 0) heap_keys[0] = key_purify(*pelem_key_at(s, local_displs[i]));
    else 
    {
      heap_keys[0] = heap_keys[heap_size - 1];
      heap_sources[0] = heap_sources[heap_size - 1];
      --heap_size;
    }

    j = 0;
    HEAPIFY_INTRO();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  return 0;
}


#ifdef SL_INDEX 
slint_t mergep_heap_unpack_idxonly(packed_elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts) /* sl_proto, sl_func mergep_heap_unpack_idxonly */
{
  slkey_pure_t heap_keys[p], hkx;
  slint_t heap_sources[p], hsx;
  slint_t heap_size = 0;

  slint_t i, j, k;
  
  slindex_t local_displs[p], local_counts[p];
  
  elements_t dst;


  memcpy(local_displs, displs, p * sizeof(slindex_t));

  if (counts == NULL)
  {
    for (i = 0; i < p - 1; ++i) local_counts[i] = local_displs[i + 1] - local_displs[i];
    local_counts[p - 1] = s->size - local_displs[p - 1];

  } else memcpy(local_counts, counts, p * sizeof(slindex_t));

  for (i = 0; i < p; ++i)
  if (local_counts[i] > 0)
  {
    heap_keys[heap_size] = key_purify(*pelem_key_at(s, local_displs[i]));
    heap_sources[heap_size] = i;

    j = heap_size;
    ++heap_size;

    HEAPIFY_INTRO();
    HEAPIFY_UP();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  elem_assign(d, &dst);

  while (heap_size > 0)
  {
    i = heap_sources[0];

    /* copy min element */
/*    pelem_unpack_at(s, local_displs[i], &dst, 0);*/
    
    *dst.indices = local_displs[i];

/*    elem_inc(&dst);*/
    index_inc(dst.indices);

    --local_counts[i];
    ++local_displs[i];
    
    if (local_counts[i] > 0) heap_keys[0] = key_purify(*pelem_key_at(s, local_displs[i]));
    else 
    {
      heap_keys[0] = heap_keys[heap_size - 1];
      heap_sources[0] = heap_sources[heap_size - 1];
      --heap_size;
    }

    j = 0;
    HEAPIFY_INTRO();
    HEAPIFY_DOWN();
    HEAPIFY_EXTRO();
  }

  return 0;
}
#endif


#undef HEAPIFY_INTRO
#undef HEAPIFY_EXTRO
#undef HEAPIFY_UP
#undef HEAPIFY_DOWN



#include "sl_common.h"


/*#define DUMMY_TEST*/

#ifdef DUMMY_TEST

slint_t dummy_tloc(elements_t *b, slint_t x, void *tloc_data)
{
  return 0;
}

permute_generic_t pg_dummy_tloc = PERMUTE_GENERIC_INIT_TLOC(dummy_tloc);
PERMUTE_GENERIC_DEFINE_TLOC(dummy_tloc)
permute_generic_t pg_dummy_ext_tloc_ext = PERMUTE_GENERIC_INIT_EXT_TLOC(dummy_tloc);

slint_t dummy_tloc_mod(elements_t *b, slint_t x, void *tloc_data, elements_t *mod)
{
  return 0;
}

permute_generic_t pg_dummy_tloc_mod = PERMUTE_GENERIC_INIT_TLOC_MOD(dummy_tloc_mod);
PERMUTE_GENERIC_DEFINE_TLOC_MOD(dummy_tloc_mod)
permute_generic_t pg_dummy_ext_tloc_mod = PERMUTE_GENERIC_INIT_EXT_TLOC_MOD(dummy_tloc_mod);

#endif


slint_t permute_generic_db(elements_t *s, elements_t *d, permute_generic_t *pg, void *pg_data) /* sl_proto, sl_func permute_generic_db */
{
  SPEC_DECLARE_TLOC_REARRANGE_DB
  SPEC_DECLARE_TLOC_MOD_REARRANGE_DB
  
  elements_t mod, *ib = NULL;
  
  
  if (pg->type == 2) { elements_alloc(&mod, 1, SLCM_ALL); ib = &mod; }

  switch (pg->type)
  {
    case 1:
      if (pg->tloc_rearrange_db) pg->tloc_rearrange_db(s, d, pg_data);
      else SPEC_DO_TLOC_REARRANGE_DB(pg->tloc, pg_data, s, d);
      break;
    case 2:
      if (pg->tloc_mod_rearrange_db) pg->tloc_mod_rearrange_db(s, d, pg_data, ib);
      else SPEC_DO_TLOC_MOD_REARRANGE_DB(pg->tloc_mod, pg_data, s, d, ib);
      break;
  }

  if (pg->type == 2) elements_free(&mod);
  
  return 0;
}


slint_t permute_generic_ip(elements_t *s, elements_t *x, permute_generic_t *pg, void *pg_data) /* sl_proto, sl_func permute_generic_ip */
{
  SPEC_DECLARE_TLOC_REARRANGE_IP
  SPEC_DECLARE_TLOC_MOD_REARRANGE_IP

  elements_t mod, *ib = NULL;
  
  
  if (pg->type == 2) { elements_alloc(&mod, 1, SLCM_ALL); ib = &mod; }

  switch (pg->type)
  {
    case 1:
      if (pg->tloc_rearrange_ip) pg->tloc_rearrange_ip(s, x, pg_data);
      else SPEC_DO_TLOC_REARRANGE_IP(pg->tloc, pg_data, s, x);
      break;
    case 2:
      if (pg->tloc_mod_rearrange_ip) pg->tloc_mod_rearrange_ip(s, x, pg_data, ib);
      else SPEC_DO_TLOC_MOD_REARRANGE_IP(pg->tloc_mod, pg_data, s, x, ib);
      break;
  }

  if (pg->type == 2) elements_free(&mod);

  return 0;
}



#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include "z_pack.h"
#include "prx.h"


void prx_seed(prxint_t seed) /* prx_func prx_seed */
{
  z_srand(seed);
}


void prx_permutation(prxint_t *permutation, prxint_t n, prx_type_t type) /* prx_func prx_permutation */
{
  prxint_t i, j, t;


  switch (type)
  {
    case PRX_FISHER_YATES_SHUFFLE:
      for (i = 0; i < n; ++i) permutation[i] = i;
      for (i = n - 1; i > 0; --i)
      {
        j = (prxint_t) z_rand_minmax(0, i);
        t = permutation[i]; permutation[i] = permutation[j]; permutation[j] = t;
      }
      break;

    default:
      fprintf(stderr, "ERROR: unknown prx type '%d'\n", (int) type);
  }
}


struct _prx_enumerate_t
{
  prx_type_t type;
  prxint_t n;
  prxint_t *permutation;
};


void prx_enumerate_create(prx_enumerate_t *enumerate, prxint_t n, prx_type_t type) /* prx_func prx_enumerate_create */
{
  *enumerate = z_alloc(1, sizeof(**enumerate));

  (*enumerate)->type = type;
  (*enumerate)->n = n;
  (*enumerate)->permutation = NULL;

  switch (type)
  {
    case PRX_FISHER_YATES_SHUFFLE:
      (*enumerate)->permutation = z_alloc(n, sizeof(prxint_t));
      prx_permutation((*enumerate)->permutation, n, type);
      break;

    default:
      fprintf(stderr, "ERROR: unknown prx type '%d'\n", (int) type);
      prx_enumerate_destroy(enumerate);
      return;
  }

/*  prx_enumerate_print(*enumerate);*/
}


void prx_enumerate_destroy(prx_enumerate_t *enumerate) /* prx_func prx_enumerate_destroy */
{
  if (enumerate == PRX_ENUMERATE_NULL) return;

  if ((*enumerate)->permutation) z_free((*enumerate)->permutation);

  z_free(*enumerate);

  *enumerate = PRX_ENUMERATE_NULL;
}


void prx_enumerate_print(prx_enumerate_t enumerate) /* prx_func prx_enumerate_print */
{
  prxint_t i;


  printf("enumerate %p:\n", enumerate);
  for (i = 0; i < enumerate->n; ++i) printf("  %" prxint_fmt " -> %" prxint_fmt "\n", i, prx_enumerate(enumerate, i));
}


prxint_t prx_enumerate(prx_enumerate_t enumerate, prxint_t i) /* prx_func prx_enumerate */
{
  switch (enumerate->type)
  {
    case PRX_FISHER_YATES_SHUFFLE:
      return enumerate->permutation[i];

    default:
      fprintf(stderr, "ERROR: unknown prx type '%d'\n", (int) enumerate->type);
  }

  return -1;
}



#include "sl_common.h"


/*
  ..._lt versions return:
   - max number of elements less than k
   - index i with s->keys[i-1] < k <= s->keys[i]
  ..._le versions return:
   - max number of elements less than or equal k
   - index i with s->keys[i-1] <= k < s->keys[i]
*/


/* sequential search routines, max comparisons: O(n) */

slint sl_search_sequential_lt(elements_t *s, slpkey_t k) /* sl_proto, sl_func sl_search_sequential_lt */
{
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign(s->keys, _s);
  key_assign_at(s->keys, s->size, _e);

  while (_s != _e)
  if (key_pure_cmp_lt(key_purify(*_s), k)) key_inc(_s); else break;

  return s->size - (_e - _s);
}


slint sl_search_sequential_le(elements_t *s, slpkey_t k) /* sl_proto, sl_func sl_search_sequential_le */
{
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign(s->keys, _s);
  key_assign_at(s->keys, s->size, _e);

  while (_s != _e)
  if (key_pure_cmp_le(key_purify(*_s), k)) key_inc(_s); else break;

  return s->size - (_e - _s);
}


slint sl_search_sequential_gt(elements_t *s, slpkey_t k) /* sl_proto, sl_func sl_search_sequential_gt */
{
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, s->size - 1, _s);
  key_assign_at(s->keys, -1, _e);

  while (_s != _e)
  if (key_pure_cmp_gt(key_purify(*_s), k)) key_dec(_s); else break;

  return s->size - (_s - _e);
}


slint sl_search_sequential_ge(elements_t *s, slpkey_t k) /* sl_proto, sl_func sl_search_sequential_ge */
{
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, s->size - 1, _s);
  key_assign_at(s->keys, -1, _e);

  while (_s != _e)
  if (key_pure_cmp_ge(key_purify(*_s), k)) key_dec(_s); else break;

  return s->size - (_s - _e);
}


slint sl_search_p_sequential_lt(elements_t *s, slpkey_t *k) /* sl_proto, sl_func sl_search_p_sequential_lt */
{
  return sl_search_sequential_lt(s, *k);
}


slint sl_search_p_sequential_le(elements_t *s, slpkey_t *k) /* sl_proto, sl_func sl_search_p_sequential_le */
{
  return sl_search_sequential_le(s, *k);
}


slint sl_search_p_sequential_gt(elements_t *s, slpkey_t *k) /* sl_proto, sl_func sl_search_p_sequential_gt */
{
  return sl_search_sequential_gt(s, *k);
}


slint sl_search_p_sequential_ge(elements_t *s, slpkey_t *k) /* sl_proto, sl_func sl_search_p_sequential_ge */
{
  return sl_search_sequential_ge(s, *k);
}


/* binary search routines, max comparisons: O(\log_2 n) */

slint sl_search_binary_lt(elements_t *s, slpkey_t k) /* sl_proto, sl_func sl_search_binary_lt */
{
  slint le, ri, mi;

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_le(k, key_purify(s->keys[mi]))) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


slint sl_search_binary_le(elements_t *s, slpkey_t k) /* sl_proto, sl_func sl_search_binary_le */
{
  slint le, ri, mi;

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_lt(k, key_purify(s->keys[mi]))) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


slint sl_search_binary_gt(elements_t *s, slpkey_t k) /* sl_proto, sl_func sl_search_binary_gt */
{
  return s->size - sl_search_binary_le(s, k);
}


slint sl_search_binary_ge(elements_t *s, slpkey_t k) /* sl_proto, sl_func sl_search_binary_ge */
{
  return s->size - sl_search_binary_lt(s, k);
}


slint sl_search_p_binary_lt(elements_t *s, slpkey_t *k) /* sl_proto, sl_func sl_search_p_binary_lt */
{
  return sl_search_binary_lt(s, *k);
}


slint sl_search_p_binary_le(elements_t *s, slpkey_t *k) /* sl_proto, sl_func sl_search_p_binary_le */
{
  return sl_search_binary_le(s, *k);
}


slint sl_search_p_binary_gt(elements_t *s, slpkey_t *k) /* sl_proto, sl_func sl_search_p_binary_gt */
{
  return sl_search_binary_gt(s, *k);
}


slint sl_search_p_binary_ge(elements_t *s, slpkey_t *k) /* sl_proto, sl_func sl_search_p_binary_ge */
{
  return sl_search_binary_ge(s, *k);
}


slint_t sl_search_binary_lt_bmask(elements_t *s, slpkey_t k, slpkey_t bmask) /* sl_proto, sl_func sl_search_binary_lt_bmask */
{
  slint le, ri, mi;

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_le(k, key_purify(s->keys[mi]) & bmask)) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


slint_t sl_search_binary_le_bmask(elements_t *s, slpkey_t k, slpkey_t bmask) /* sl_proto, sl_func sl_search_binary_le_bmask */
{
  slint le, ri, mi;

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_lt(k, key_purify(s->keys[mi]) & bmask)) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


slint_t sl_search_binary_sign_switch(elements_t *s) /* sl_proto, sl_func sl_search_binary_sign_switch */
{
  slint le, ri, mi;

  if (s->size <= 0) return 0;
  
  if (key_pure_cmp_ge(key_purify(s->keys[0]), 0) && key_pure_cmp_ge(key_purify(s->keys[s->size - 1]), 0)) return s->size;

  if (key_pure_cmp_lt(key_purify(s->keys[0]), 0) && key_pure_cmp_lt(key_purify(s->keys[s->size - 1]), 0)) return 0;
  
  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_ge(key_purify(s->keys[mi]), 0)) le = mi + 1;
    else ri = mi - 1;
  }
  
  return le;
}


/* hybrid search routines combining sequential steps of size t with a following binary search */

slint sl_search_hybrid_lt(elements_t *s, slpkey_t k, slint t) /* sl_proto, sl_func sl_search_hybrid_lt */
{
  slint n;
  elements_t x;
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, t - 1, _s);
  key_assign_at(s->keys, s->size, _e);

  while (_s < _e)
  if (key_pure_cmp_lt(key_purify(*_s), k)) key_add(_s, t); else break;

  n = (_s - s->keys) - (t - 1);

  elem_assign_at(s, n, &x);
  x.size = (z_min(_s, _e) - s->keys) - n;

  return n + sl_search_binary_lt(&x, k);
}


slint sl_search_hybrid_le(elements_t *s, slpkey_t k, slint t) /* sl_proto, sl_func sl_search_hybrid_le */
{
  slint n;
  elements_t x;
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, t - 1, _s);
  key_assign_at(s->keys, s->size, _e);

  while (_s < _e)
  if (key_pure_cmp_le(key_purify(*_s), k)) key_add(_s, t); else break;

  n = (_s - s->keys) - (t - 1);

  elem_assign_at(s, n, &x);
  x.size = (z_min(_s, _e) - s->keys) - n;

  return n + sl_search_binary_le(&x, k);
}


slint sl_search_hybrid_gt(elements_t *s, slpkey_t k, slint t) /* sl_proto, sl_func sl_search_hybrid_gt */
{
  slint n;
  elements_t x;
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, s->size - t, _s);
  key_assign_at(s->keys, -1, _e);

  while (_s > _e)
  if (key_pure_cmp_gt(key_purify(*_s), k)) key_sub(_s, t); else break;

  n = (z_max(_e, _s) - s->keys) + 1;

  elem_assign_at(s, n, &x);
  x.size = ((_s - s->keys) + t) - n;

  return s->size - (n + sl_search_binary_le(&x, k));
}


slint sl_search_hybrid_ge(elements_t *s, slpkey_t k, slint t) /* sl_proto, sl_func sl_search_hybrid_ge */
{
  slint n;
  elements_t x;
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, s->size - t, _s);
  key_assign_at(s->keys, -1, _e);

  while (_s > _e)
  if (key_pure_cmp_ge(key_purify(*_s), k)) key_sub(_s, t); else break;

  n = (z_max(_e, _s) - s->keys) + 1;

  elem_assign_at(s, n, &x);
  x.size = ((_s - s->keys) + t) - n;

  return s->size - (n + sl_search_binary_lt(&x, k));
}


slint sl_search_p_hybrid_lt(elements_t *s, slpkey_t *k, slint t) /* sl_proto, sl_func sl_search_p_hybrid_lt */
{
  return sl_search_hybrid_lt(s, *k, t);
}


slint sl_search_p_hybrid_le(elements_t *s, slpkey_t *k, slint t) /* sl_proto, sl_func sl_search_p_hybrid_le */
{
  return sl_search_hybrid_le(s, *k, t);
}


slint sl_search_p_hybrid_gt(elements_t *s, slpkey_t *k, slint t) /* sl_proto, sl_func sl_search_p_hybrid_gt */
{
  return sl_search_hybrid_gt(s, *k, t);
}


slint sl_search_p_hybrid_ge(elements_t *s, slpkey_t *k, slint t) /* sl_proto, sl_func sl_search_p_hybrid_ge */
{
  return sl_search_hybrid_ge(s, *k, t);
}



#include "sl_common.h"

sl_context_t sl_default_context =  /* sl_global sl_var sl_default_context */
{
#include "sl_context_init.h"
};

/* sl_context CONTEXT_BEGIN sl */
const int default_sl_dummy_rank = -2;  /* sl_global sl_context sl_var default_sl_dummy_rank */
/* sl_context CONTEXT_END sl */

const rti_t default_rti = { };  /* sl_ifdef SL_USE_RTI sl_endif sl_global sl_context sl_var default_rti */


slint ilog2c(slint x) /* sl_proto, sl_func ilog2c */
{
  slint l = 0;
  while (x /= 2) ++l;
  return l;
}


slint ilog2f(slint x) /* sl_proto, sl_func ilog2f */
{
  slint l = 0;
  if (x <= 1) return 0;
  --x;
  while (x /= 2) ++l;
  return l + 1;
}


slint print_bits(slint v) /* sl_proto, sl_func print_bits */
{
  slint i;

  for (i = sizeof(v) * 8 - 1; i >= 0; i--) printf("%d", (v & (1L << i)) != 0);

  return 0;
}


slint pivot_random(elements_t *s) /* sl_proto, sl_func pivot_random */
{
  return z_rand() % s->size;
}


slint_t counts2displs(slint_t n, int *counts, int *displs) /* sl_proto, sl_func counts2displs */
{
  slint_t i, j;

  if (displs)
  {
    displs[0] = 0;
    for (i = 1; i < n; ++i) displs[i] = displs[i - 1] + counts[i - 1];

  } else
  {
    j = 0;
    for (i = 0; i < n; ++i) counts[i] = (j += counts[i]) - counts[i];
  }
  
  return 0;
}


slint_t displs2counts(slint_t n, int *displs, int *counts, slint_t total_counts) /* sl_proto, sl_func displs2counts */
{
  slint_t i;

  if (counts)
  {
    for (i = 0; i < n - 1; ++i) counts[i] = displs[i + 1] - displs[i];
    counts[n - 1] = total_counts - displs[n - 1];

  } else
  {
    for (i = 0; i < n - 1; ++i) displs[i] = displs[i + 1] - displs[i];
    displs[n - 1] = total_counts - displs[n - 1];
  }
  
  return 0;
}


void get_displcounts_extent(slint_t n, int *displs, int *counts, slint_t *lb, slint_t *extent) /* sl_proto, sl_func get_displcounts_extent */
{
  slint_t i, _lb, _ub;

  if (n <= 0) return;

  _lb = displs[0];
  _ub = displs[0] + counts[0];

  for (i = 1; i < n; ++i)
  {
    if (displs[i] < _lb) _lb = displs[i];
    if (displs[i] + counts[i] > _ub) _ub = displs[i] + counts[i];
  }
  
  if (lb) *lb = _lb;
  if (extent) *extent = _ub - _lb;
}


/* calculates the intersection of [a0_start,a0_end] with [a1_start,a1_end]
   - storing the intersection in [ia_start,ia_end]
   - returning 'ia_end - ia_start + 1'
     x > 0: intersection consists of x elements
     x = 0: no intersection, areas touching
     x < 0: no intersection, x elements inbetween */
/*slint intersect_areas(slint a0_start, slint a0_end, slint a1_start, slint a1_end, slint *ia_start, slint *ia_end)
{
  slint temp_ia_start, temp_ia_end;

  if (ia_start == NULL) ia_start = &temp_ia_start;
  if (ia_end == NULL) ia_end = &temp_ia_end;

  *ia_start = z_max(a0_start, a1_start);
  *ia_end = z_min(a0_end, a1_end);

  return *ia_end - *ia_start + 1;
}*/

/* calculates the remaining area(s) when discarding [a1_start,a1_end] from [a0_start,a0_end]
   - storing the remaining area(s) in [a2_start,a2_end] and [a3_start,a3_end]
   - 'a2_end < a1_start && a1_end < a3_start' is always true
   - returning the added size of [a2_start,a2_end] and [a3_start,a3_end] */
/*slint subtract_areas(slint a0_start, slint a0_end, slint a1_start, slint a1_end, slint *a2_start, slint *a2_end, slint *a3_start, slint *a3_end)
{
  slint temp_a2_start, temp_a2_end, temp_a3_start, temp_a3_end;

  if (a2_start == NULL) a2_start = &temp_a2_start;
  if (a2_end == NULL) a2_end = &temp_a2_end;
  if (a3_start == NULL) a3_start = &temp_a3_start;
  if (a3_end == NULL) a3_end = &temp_a3_end;

  *a2_start = a0_start;
  *a2_end = z_min(a1_start - 1, a0_end);

  *a3_start = z_max(a1_end + 1, a0_start);
  *a3_end = a0_end;

  return z_max(0, *a2_end - *a2_start + 1) + z_max(0, *a3_end - *a3_start + 1);
}*/



#include "sl_common.h"

#include <stdarg.h>


#define ELEM_REVERSE  0  /* 0 = elem_reverse_aio, 1 = elem_reverse_obo */
#define ELEM_ROTATE   0  /* 0 = elem_rotate_3rev_aio, 1 = elem_rotate_3rev_obo, 2 = elem_rotate_cycles_aio, 3 = elem_rotate_cycles_obo */


void elem_set_data(elements_t *e, ...) /* sl_proto, sl_func elem_set_data */
{
  va_list ap;

  va_start(ap, e);

#define xelem_call_data  xelem_buf(e) = va_arg(ap, xelem_sltype_t *);
#include "sl_xelem_call.h"

  va_end(ap);
}


slint_t elem_get_max_byte() /* sl_proto, sl_func elem_get_max_byte */
{
  slint_t x = 0;

#define xelem_call  if (xelem_byte > x) x = xelem_byte;
#include "sl_xelem_call.h"
  
  return x;
}


#if ELEM_REVERSE == 0 || ELEM_ROTATE == 0
static slint_t elem_reverse_aio(elements_t *e, elements_t *t) /* sl_func elem_reverse_aio */
{
  elements_t front, back, end;

  elem_assign(e, &front);
  elem_assign_at(e, e->size - 1, &back);
  elem_assign_at(e, e->size / 2, &end);

  while (front.keys < end.keys)
  {
    elem_xchange(&front, &back, t);
    elem_inc(&front);
    elem_dec(&back);
  }

  return 0;
}
#endif


#if ELEM_REVERSE == 1 || ELEM_ROTATE == 1
static slint_t elem_reverse_obo(elements_t *e, elements_t *t) /* sl_func elem_reverse_obo */
{
  elements_t front, back, end;

  elem_assign(e, &front);
  elem_assign_at(e, e->size - 1, &back);
  elem_assign_at(e, e->size / 2, &end);

#define xelem_call \
  while (xelem_buf(&front) > xelem_buf(&end)) \
  { \
    xelem_xchange(&front, &back, t); \
    xelem_inc(&front); \
    xelem_dec(&back); \
  }
#include "sl_xelem_call.h"

  return 0;
}
#endif


slint_t elem_reverse(elements_t *e, elements_t *t) /* sl_proto, sl_func elem_reverse */
{
#if ELEM_REVERSE == 0
  return elem_reverse_aio(e, t);
#elif ELEM_REVERSE == 1
  return elem_reverse_obo(e, t);
#else
# error ELEM_REVERSE unknown
#endif
}


slint_t elem_nxchange_at(elements_t *e0, slint_t at0, elements_t *e1, slint_t at1, slint_t n, elements_t *t) /* sl_proto, sl_func elem_nxchange_at */
{
  elements_t _e0, _e1, end;

  elem_assign_at(e0, at0, &_e0);
  elem_assign_at(e1, at1, &_e1);
  elem_assign_at(&_e0, n, &end);

  while (_e0.keys < end.keys)
  {
    elem_xchange(&_e0, &_e1, t);
    elem_inc(&_e0);
    elem_inc(&_e1);
  }

  return 0;
}


slint_t elem_nxchange(elements_t *e0, elements_t *e1, slint_t n, elements_t *t) /* sl_proto, sl_func elem_nxchange */
{
  return elem_nxchange_at(e0, 0, e1, 0, n, t);
}


slint_t elem_nxchange_ro0(elements_t *e0, elements_t *e1, slint_t n, elements_t *t) /* sl_proto, sl_func elem_nxchange_ro0 */
{
  elements_t _e0, _e1, end;

  elem_assign(e0, &_e0);
  elem_assign(e1, &_e1);
  elem_assign_at(e1, n, &end);

  elem_copy(&_e1, t); /* create the hole */
  elem_copy(&_e0, &_e1);
  elem_inc(&_e1);

  while (_e1.keys < end.keys)
  {
    elem_copy(&_e1, &_e0);
    elem_inc(&_e0);
    elem_copy(&_e0, &_e1);
    elem_inc(&_e1);
  }

  elem_copy(t, &_e0); /* close the hole */

  return 0;
}


#if ELEM_ROTATE == 0
static slint_t elem_rotate_3rev_aio(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_func elem_rotate_3rev_aio */
{
  elements_t _e;

  if (m == 0 || n == 0) return 0;

  /* reverse 2nd part */
  elem_assign_at(e, m, &_e);
  _e.size = n;
  elem_reverse_aio(&_e, t);

  /* reverse 1st part */
  elem_assign(e, &_e);
  _e.size = m;
  elem_reverse_aio(&_e, t);

  /* reverse all */
  _e.size = m + n;
  elem_reverse_aio(&_e, t);

  return 0;
}
#endif


#if ELEM_ROTATE == 1
static slint_t elem_rotate_3rev_obo(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_func elem_rotate_3rev_obo */
{
  elements_t _e;

  if (m == 0 || n == 0) return 0;

  /* reverse 2nd part */
  elem_assign_at(e, m, &_e);
  _e.size = n;
  elem_reverse_obo(&_e, t);

  /* reverse 1st part */
  elem_assign(e, &_e);
  _e.size = m;
  elem_reverse_obo(&_e, t);

  /* reverse all */
  _e.size = m + n;
  elem_reverse_obo(&_e, t);

  return 0;
}
#endif


#if ELEM_ROTATE == 2
static slint_t elem_rotate_cycles_aio(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_func elem_rotate_cycles_aio */
{
  slint_t k;
  elements_t start, half, current, next;

  if (m == 0 || n == 0) return 0;

  elem_assign(e, &start);
  elem_assign_at(e, n, &half);

  k = m + n;
  while (k > 0)
  {
    elem_copy(&start, t);
    elem_assign(&start, &current);

    while (1)
    {
      if (current.keys < half.keys) elem_assign_at(&current, m, &next);
      else elem_assign_at(&current, -n, &next);

      if (next.keys == start.keys) break;

      elem_copy(&next, &current);
      elem_assign(&next, &current);
      --k;
    }

    elem_copy(t, &current);
    --k;

    elem_inc(&start);
  }

  return 0;
}
#endif


#if ELEM_ROTATE == 3
static slint_t elem_rotate_cycles_obo(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_func elem_rotate_cycles_obo */
{
  slint_t k;
  elements_t start, half, current, next;

  if (m == 0 || n == 0) return 0;

  elem_assign(e, &start);
  elem_assign_at(e, n, &half);

#define xelem_call \
  k = m + n; \
  while (k > 0) \
  { \
    xelem_copy(&start, t); \
    xelem_assign(&start, &current); \
\
    while (1) \
    { \
      if (xelem_buf(&current) < xelem_buf(&half)) xelem_assign_at(&current, m, &next); \
      else xelem_assign_at(&current, -n, &next); \
\
      if (xelem_buf(&next) == xelem_buf(&start)) break; \
\
      xelem_copy(&next, &current); \
      xelem_assign(&next, &current); \
      --k; \
    } \
\
    xelem_copy(t, &current); \
    --k; \
\
    xelem_inc(&start); \
  }
#include "sl_xelem_call.h"

  return 0;
}
#endif


slint_t elem_rotate(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_proto, sl_func elem_rotate */
{
#if ELEM_ROTATE == 0
  return elem_rotate_3rev_aio(e, m, n, t);
#elif ELEM_ROTATE == 1
  return elem_rotate_3rev_obo(e, m, n, t);
#elif ELEM_ROTATE == 2
  return elem_rotate_cycles_aio(e, m, n, t);
#elif ELEM_ROTATE == 3
  return elem_rotate_cycles_obo(e, m, n, t);
#else
# error ELEM_ROTATE unknown
#endif
}


/* retain order of the 1st part, do it back2front (for the sake of simplicity) */
slint_t elem_rotate_ro0(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_proto, sl_func elem_rotate_ro0 */
{
  elements_t e0, e1;

  elem_assign_at(e, m, &e0);
  elem_assign_at(e, m + n, &e1);

  elem_copy(&e1, t);
  elem_copy(&e0, &e1);

  while (e0.keys > e->keys)
  {
    elem_dec(&e1);
    elem_copy(&e1, &e0);
    elem_dec(&e0);
    elem_copy(&e0, &e1);
  }

  elem_copy(t, &e0);

  return 0;
}


/* retain order of the 2nd part, do it front2back (for the sake of simplicity) */
slint_t elem_rotate_ro1(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_proto, sl_func elem_rotate_ro1 */
{
  elements_t e0, e1, end;

  elem_assign(e, &e0);
  elem_assign_at(e, m, &e1);
  elem_assign_at(e, m + n, &end);

  elem_copy(&e0, t);
  elem_copy(&e1, &e0);

  while (e1.keys < end.keys)
  {
    elem_inc(&e0);
    elem_copy(&e0, &e1);
    elem_inc(&e1);
    elem_copy(&e1, &e0);
  }

  elem_copy(t, &e1);

  return 0;
}


#undef ELEM_REVERSE
#undef ELEM_ROTATE



#include "sl_common.h"


static void make_counts(elements_t *s, slint_t ncounts, slint_t *counts)
{
  slint_t i;
  key_type_c *k;

  memset(counts, 0, ncounts * sizeof(slint_t));

  for (i = 0, k = s->keys; i < s->size; ++i, ++k) ++counts[key_purify(*k)];
}


static void make_displs(slint_t ncounts, slint_t *counts, slint_t *displs)
{
  slint_t i;

  displs[0] = 0;
  for (i = 1; i < ncounts; ++i) displs[i] = displs[i - 1] + counts[i - 1];
}


static void make_counts2displs(slint_t ncounts, slint_t *countsdispls)
{
  slint_t i, s, t;

  s = 0;
  for (i = 0; i < ncounts; ++i)
  {
    t = countsdispls[i];
    countsdispls[i] = s;
    s += t;
  }
}


slint_t sort_counting_use_displs(elements_t *s, elements_t *d, slint_t ndispls, slint_t *displs)  /* sl_proto, sl_func sort_counting_use_displs */
{
  slint_t i;

  for (i = 0; i < s->size; ++i)
  {
    elem_copy_at(s, i, d, displs[key_purify(s->keys[i])]);
    ++displs[key_purify(s->keys[i])];
  }

  d->size = s->size;
  
  return 0;
}


slint_t sort_counting_use_counts(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts)  /* sl_proto, sl_func sort_counting_use_counts */
{
  slint_t r;
  slint_t *displs = NULL;

  if (counts == NULL)
  {
    displs = z_alloc(ncounts, sizeof(slint_t));
    make_counts(s, ncounts, displs);
    make_counts2displs(ncounts, displs);

  } else
  {
    displs = counts;

    if (ncounts < 0)
    {
      ncounts *= -1;
      displs += ncounts;
      make_displs(ncounts, counts, displs);

    } else make_counts2displs(ncounts, displs);
  }
  
  r = sort_counting_use_displs(s, d, ncounts, displs);

  if (counts == NULL) z_free(displs);

  return r;
}


slint_t sort_counting_get_counts(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts)  /* sl_proto, sl_func sort_counting_get_counts */
{
  slint_t r;
  slint_t *displs = NULL;
  
  if (counts == NULL)
  {
    displs = z_alloc(ncounts, sizeof(slint_t));
    make_counts(s, ncounts, displs);
    make_counts2displs(ncounts, displs);

  } else
  {
    if (ncounts < 0)
    {
      ncounts *= -1;
      displs = counts + ncounts;

    } else displs = z_alloc(ncounts, sizeof(slint_t));

    make_counts(s, ncounts, counts);
    make_displs(ncounts, counts, displs);
  }

  r = sort_counting_use_displs(s, d, ncounts, displs);
  
  if (counts == NULL || displs != counts + ncounts) z_free(displs);
  
  return r;
}


slint_t sort_counting(elements_t *s, elements_t *d, slint_t ncounts)  /* sl_proto, sl_func sort_counting */
{
  return sort_counting_use_counts(s, d, ncounts, NULL);
}



#include "sl_common.h"


#if 0
static void hs_heapify_i(elements_t *s, slint i, slint size, elements_t *xs)
{
  slint l, r, largest;

  elem_copy_at(s, i, xs, 0);

  while ((l = 2 * i + 1) < size)
  {
    largest = l;
    r = 2 * i + 2;

    if (r < size)
    if (key_cmp_gt(s->keys[r], s->keys[l])) largest = r;

    if (key_cmp_gt(*xs->keys, s->keys[largest])) break;

    elem_copy_at(s, largest, s, i);

    i = largest;
  }

  elem_copy_at(xs, 0, s, i);
}


static slint sort_heap_i(elements_t *s, elements_t *xs) /* sl_func sort_heap_i */
{
  slint i;
  elements_t txs;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size < 2) return 0;

  /* build the heap */
  for (i = s->size / 2 - 1; i >= 0; --i) hs_heapify_i(s, i, s->size, xs);

  /* extract the maxima */
  for (i = s->size - 1; i > 0; --i)
  {
    elem_xchange_at(s, 0, s, i, xs);
    hs_heapify_i(s, 0, i, xs);
  }

  if (xs == &txs) elements_free(xs);

  return 0;
}
#endif


static slint sort_heap_p(elements_t *s, elements_t *xs) /* sl_func sort_heap_p */
{
  elements_t txs, last, si, sj, l;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size < 2) return 0;

  elem_assign_at(s, s->size, &last);

  /* build the heap */
  for (elem_assign_at(s, s->size / 2 - 1, &si); si.keys >= s->keys; elem_dec(&si))
  {
    elem_assign(&si, &sj);

    elem_copy(&sj, xs);

    while (1)
    {
      elem_assign_at(s, (sj.keys - s->keys) * 2 + 1, &l); /* l = i * 2 + 1; */

      if (l.keys > last.keys) break;

      if (l.keys < last.keys)
      if (key_cmp_gt(l.keys[1], *l.keys)) elem_inc(&l);

      if (key_cmp_gt(*xs->keys, *l.keys)) break;

      elem_copy(&l, &sj);

      elem_assign(&l, &sj);
    }

    elem_copy(xs, &sj);
  }

  /* extract the maxima */
  for (elem_assign_at(s, s->size - 1, &si); si.keys > s->keys;)
  {
    elem_copy(&si, xs);
    elem_copy(s, &si);

    /* si becomes the last element in the heap */
    elem_dec(&si);

    elem_assign(s, &sj);

    while (1)
    {
      elem_assign_at(s, (sj.keys - s->keys) * 2 + 1, &l); /* l = i * 2 + 1; */

      if (l.keys > si.keys) break;

      if (l.keys < si.keys)
      if (key_cmp_gt(l.keys[1], *l.keys)) elem_inc(&l);

      if (key_cmp_gt(*xs->keys, *l.keys)) break;

      elem_copy(&l, &sj);

      elem_assign(&l, &sj);
    }

    elem_copy(xs, &sj);
  }

  if (xs == &txs) elements_free(xs);

  return 0;
}


slint sort_heap(elements_t *s, elements_t *xs) /* sl_proto, sl_func sort_heap */
{
  return sort_heap_p(s, xs);
}



#include "sl_common.h"


#ifdef key_integer

slint_t sort_insert_bmask_kernel(elements_t *s, elements_t *sx, slkey_pure_t bmask) /* sl_proto, sl_func sort_insert_bmask_kernel */
{
  slint_t i, j;

  for (i = 1; i < s->size; i++)
  {
    if (key_pure_cmp_lt(key_purify(s->keys[i]) & bmask, key_purify(s->keys[i - 1]) & bmask))
    {
      j = i - 1;
      elem_copy_at(s, i, sx, 0);

      do
      {
        elem_copy_at(s, j, s, j + 1);
        if (--j < 0) break;

      } while (key_pure_cmp_lt(key_purify(*sx->keys) & bmask, key_purify(s->keys[j]) & bmask));

      elem_copy_at(sx, 0, s, j + 1);
    }
  }

  return 0;
}

#endif


static slint_t sort_insert_kernel(elements_t *s, elements_t *sx) /* sl_func sort_insert_kernel */
{
  slint_t i, j;

  for (i = 1; i < s->size; i++)
  {
    if (key_cmp_lt(s->keys[i], s->keys[i - 1]))
    {
      j = i - 1;
      elem_copy_at(s, i, sx, 0);

      do
      {
        elem_copy_at(s, j, s, j + 1);
        if (--j < 0) break;

      } while (key_cmp_lt(*sx->keys, s->keys[j]));

      elem_copy_at(sx, 0, s, j + 1);
    }
  }

  return 0;
}


slint_t sort_insert(elements_t *s, elements_t *sx) /* sl_proto, sl_func sort_insert */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_insert);

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);
  }

  sort_insert_kernel(s, sx);

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_insert);

  return 0;
}



#include "sl_common.h"


static slint_t sort_permute_forward_(elements_t *s, elements_t *sx, slint_t *perm, slint_t offset) /* sl_func sort_permute_forward_ */
{
  elements_t src, end;

  slint_t j, i, *ia, *ja;

  elem_assign(s, &src);
  elem_assign_at(s, s->size, &end);

  i = offset;
  ia = perm;
  while (src.keys != end.keys)
  {
    while (*ia != i)
    {
      j = *ia - offset;

      elem_xchange_at(&src, 0, s, j, sx);

      ja = perm + j;
      *ia = *ja;
      *ja = j + offset;
    }

    ia++;
    i++;
    elem_inc(&src);
  }

  return 0;
}


static slint_t sort_permute_forward_masked(elements_t *s, elements_t *sx, slint_t *perm, slint_t offset, slint_t bit_mask) /* sl_func sort_permute_forward_masked */
{
  elements_t src, end;

  slint_t inv_bit_mask = bit_mask ^ -1L;
  slint_t j, i, *ia, *ja;

  elem_assign(s, &src);
  elem_assign_at(s, s->size, &end);

  i = 0;
  ia = perm;
  while (src.keys != end.keys)
  {
    if ((*ia & bit_mask) == 0)
    {
      j = *ia - offset;

      while (i != j)
      {
        elem_xchange_at(&src, 0, s, j, sx);

        ja = perm + j;
        j = *ja - offset;
        *ja |= bit_mask;
      }
    }

    *ia &= inv_bit_mask;

    ia++;
    i++;
    elem_inc(&src);
  }

  return 0;
}


slint_t sort_permute_forward(elements_t *s, elements_t *sx, slint_t *perm, slint_t offset, slint_t mask_bit) /* sl_proto, sl_func sort_permute_forward */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_permute_forward);

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);
  }

  if (mask_bit < 0) sort_permute_forward_(s, sx, perm, offset);
  else sort_permute_forward_masked(s, sx, perm, offset, z_powof2_typed(z_min(mask_bit, (slint_t) (sizeof(slint_t) * 8) - 1), slint_t));

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_permute_forward);

  return 0;
}


static slint_t sort_permute_backward_(elements_t *s, elements_t *sx, slint_t *perm, slint_t offset) /* sl_func sort_permute_backward_ */
{
  elements_t src, end, e0, e1, *from, *to, *t;

  slint_t i, j, k, *ia, *ja;

  elem_assign(s, &src);
  elem_assign_at(s, s->size, &end);

  from = &e0;
  to = &e1;

  i = offset;
  ia = perm;
  while (src.keys != end.keys)
  {
    if (*ia != i)
    {
      elem_copy(&src, sx);
      elem_assign(&src, to);

      ja = ia;
      j = i;

      while (i != (k = *ja))
      {
        elem_assign_at(s, k - offset, from);
        elem_copy(from, to);

        t = to;
        to = from;
        from = t;

        *ja = j;
        ja = perm + (j = k) - offset;
      }

      elem_copy(sx, to);
      *ja = j;
    }

    ia++;
    i++;
    elem_inc(&src);
  }

  return 0;
}


static slint_t sort_permute_backward_masked(elements_t *s, elements_t *sx, slint_t *perm, slint_t offset, slint_t bit_mask) /* sl_func sort_permute_backward_masked */
{
  elements_t src, end, e0, e1, *from, *to, *t;

  slint_t inv_bit_mask = bit_mask ^ -1L;

  slint_t i, k, *ia, *ja;

  elem_assign(s, &src);
  elem_assign_at(s, s->size, &end);

  from = &e0;
  to = &e1;

  i = offset;
  ia = perm;
  while (src.keys != end.keys)
  {
    if ((*ia & bit_mask) == 0)
    {
      elem_copy(&src, sx);
      elem_assign(&src, to);

      ja = ia;

      while (i != (k = *ja))
      {
        elem_assign_at(s, k - offset, from);
        elem_copy(from, to);

        t = to;
        to = from;
        from = t;

        *ja |= bit_mask;
        ja = perm + k - offset;
      }

      elem_copy(sx, to);

      *ja |= bit_mask;
    }

    *ia &= inv_bit_mask;

    ia++;
    i++;
    elem_inc(&src);
  }

  return 0;
}


slint_t sort_permute_backward(elements_t *s, elements_t *sx, slint_t *perm, slint_t offset, slint_t mask_bit) /* sl_proto, sl_func sort_permute_backward */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_permute_backward);

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);
  }

  if (mask_bit < 0) sort_permute_backward_(s, sx, perm, offset);
  else sort_permute_backward_masked(s, sx, perm, offset, z_powof2_typed(z_min(mask_bit, (slint_t) (sizeof(slint_t) * 8) - 1), slint_t));

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_permute_backward);

  return 0;
}


/* quicksort from "Introduction to Algorithms"
   - iterative version with own stack requires pre-allocated stack with max. depth -> bad
   - half-recursive version (qs_halfrec) using "tail recursion" (pg. 162)
   - reducing max. stack depth to \Theta(lg n) proposed by problem 7-4c (pg. 162)
   - hybrid version using own stack of fixed size and recursive calls if the stack exceeds
*/


#include "sl_common.h"


#define SORT_QUICK  1  /* 0 = sort_quick_i, 1 = sort_quick_p, 2 = sort_quick_h */


#if SORT_QUICK == 0
static void qs_halfrec_i(elements_t *s, slint l, slint h, elements_t *xs)
{
  slint q, i, j;

  elements_t x;

  while (l < h)
  {
    /* partitioning the array */
    elem_assign_at(s, h, &x);

    for (i = j = l; j < h; j++)
    if (key_cmp_le(s->keys[j], *x.keys))
    {
      /* exchange i and j */
      elem_xchange_at(s, i, s, j, xs);
      i++;
    }
    /* exchange i and h */
    elem_xchange_at(s, i, s, h, xs);

    q = i;

    /* call recursive with the smaller part to reduce max. stack depth */
    if (q - l < h - q)
    {
      q--;
      if (l < q) qs_halfrec_i(s, l, q, xs);
      l = q + 2;

    } else
    {
      q++;
      if (q < h) qs_halfrec_i(s, q, h, xs);
      h = q - 2;
    }
  }
}


static slint sort_quick_i(elements_t *s, elements_t *xs) /* sl_func sort_quick_i */
{
  elements_t txs;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size < 2) return 0;

  qs_halfrec_i(s, 0, s->size - 1, xs);

  if (xs == &txs) elements_free(xs);

  return 0;
}
#endif


#if SORT_QUICK == 1
static void qs_halfrec_p(elements_t xl, elements_t xh, elements_t *xs)
{
  elements_t xi, xj;

  while (xl.keys < xh.keys)
  {
    /* partitioning the array */
    elem_assign(&xl, &xi);
    elem_assign(&xl, &xj);

    for (; xj.keys < xh.keys;)
    {
      if (key_cmp_le(*xj.keys, *xh.keys))
      {
        /* exchange i and j */
        elem_xchange(&xi, &xj, xs);
        elem_inc(&xi);
      }
      elem_inc(&xj);
    }

    /* exchange i and h */
    elem_xchange(&xi, &xh, xs);

    /* call recursive with the smaller part to reduce max. stack depth */
    if (xi.keys - xl.keys < xh.keys - xi.keys)
    {
      elem_dec(&xi);
      if (xl.keys < xi.keys) qs_halfrec_p(xl, xi, xs);
      elem_assign_at(&xi, 2, &xl);

    } else
    {
      elem_inc(&xi);
      if (xi.keys < xh.keys) qs_halfrec_p(xi, xh, xs);
      elem_assign_at(&xi, -2, &xh);
    }
  }
}


static slint sort_quick_p(elements_t *s, elements_t *xs) /* sl_func sort_quick_p */
{
  elements_t xh, txs;

  if (s == NULL) return -1;

  rti_tstart(rti_tid_sort_quick);

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size > 1)
  {
    elem_assign_at(s, s->size - 1, &xh);
    qs_halfrec_p(*s, xh, xs);
  }

  if (xs == &txs) elements_free(xs);

  rti_tstop(rti_tid_sort_quick);

  return 0;
}
#endif


#if SORT_QUICK == 2
static void qs_hybrid(elements_t xl, elements_t xh, elements_t *xs)
{
#define stack_size 32

  elements_t xi, xj;
  struct { elements_t xl, xh; } stack[stack_size], *sp = stack, *lsp = stack + stack_size;

#define push(l, h)  ((sp < lsp)?(sp->xl = l, sp->xh = h, ++sp, 1):0)
#define pop(l, h)   ((sp > stack)?(--sp, l = sp->xl, h = sp->xh, 1):0)

  while (1)
  {
    /* if the current part is to small */
    if (xl.keys >= xh.keys)
    if (!pop(xl, xh)) break; /* pop a new part from the stack, or break if empty */

    /* partitioning the array */
    elem_assign(&xl, &xi);
    elem_assign(&xl, &xj);

    for (; xj.keys < xh.keys;)
    {
      if (key_cmp_le(*xj.keys, *xh.keys))
      {
        /* exchange i and j */
        elem_xchange(&xi, &xj, xs);
        elem_inc(&xi);
      }
      elem_inc(&xj);
    }

    /* put the pivot in the middle */
    elem_xchange(&xi, &xh, xs);

    /* push smaller part on the stack (call recursive if full) to reduce max. stack depth */
    if (xi.keys - xl.keys < xh.keys - xi.keys)
    {
      elem_dec(&xi);

      if (xl.keys < xi.keys)
      if (!push(xl, xi)) qs_hybrid(xl, xi, xs); /* call if push fails */

      elem_assign_at(&xi, 2, &xl);

    } else
    {
      elem_inc(&xi);

      if (xi.keys < xh.keys)
      if (!push(xi, xh)) qs_hybrid(xi, xh, xs); /* call if push fails */

      elem_assign_at(&xi, -2, &xh);
    }
  }

#undef push
#undef pop
#undef stack_size
}


static slint sort_quick_h(elements_t *s, elements_t *xs) /* sl_func sort_quick_h */
{
  elements_t xh, txs;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size < 2) return 0;

  elem_assign_at(s, s->size - 1, &xh);
  qs_hybrid(*s, xh, xs);

  if (xs == &txs) elements_free(xs);

  return 0;
}
#endif


slint sort_quick(elements_t *s, elements_t *xs) /* sl_proto, sl_func sort_quick */
{
#if SORT_QUICK == 0
  return sort_quick_i(s, xs);
#elif SORT_QUICK == 1
  return sort_quick_p(s, xs);
#elif SORT_QUICK == 2
  return sort_quick_h(s, xs);
#else
# error SORT_QUICK unknown
#endif
}


#undef SORT_QUICK



/* sl_macro SR_IP_INSERTSORT */
#define SR_IP_INSERTSORT

/* sl_macro SR_DB_INSERTSORT */
#define SR_DB_INSERTSORT

/* sl_macro SR_MA_INSERTSORT */
#define SR_MA_INSERTSORT


#include "sl_common.h"

/* sl_context CONTEXT_BEGIN sr */
const slint_t default_sr_ip_threshold = sort_radix_ip_threshold;  /* sl_global sl_context sl_var default_sr_ip_threshold */
const slint_t default_sr_db_threshold = sort_radix_db_threshold;  /* sl_global sl_context sl_var default_sr_db_threshold */
const slint_t default_sr_ma_threshold = sort_radix_db_threshold;  /* sl_global sl_context sl_var default_sr_ma_threshold */
/* sl_context CONTEXT_END sr */


#ifdef key_integer


#define max_nclasses(_width_, _type_) (z_powof2_typed(_width_, _type_))


static void rs_rec_ip(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth)
{
  slkey_pure_t bit_mask, nclasses;

  slint_t i, j, k, current_width, c[max_nclasses(sort_radix_width_max, slkey_pure_t)];
  elements_t xi, end, parts[max_nclasses(sort_radix_width_max, slkey_pure_t)];

  elem_assign_at(s, s->size, &end);

  current_width = z_min(rwidth, rhigh - rlow + 1);
  rhigh -= current_width - 1;

  nclasses = z_powof2_typed(current_width, slkey_pure_t);
  bit_mask = nclasses - 1;


  /* zero all counter */
  for (i = 0; i < nclasses; i++) c[i] = 0;

  /* count the number of elements in every class */
  for (elem_assign(s, &xi); xi.keys < end.keys; elem_inc(&xi)) ++c[key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask)];

  /* compute the target of every class */
  elem_assign(s, &parts[0]);
  for (i = 1; i < nclasses; i++) elem_assign_at(&parts[i - 1], c[i - 1], &parts[i]);;

  /* split the elements */
  elem_assign(s, &end);
  for (i = 0; i < nclasses; i++)
  {
    elem_add(&end, c[i]);

    elem_assign(&parts[i], &xi);

    while (xi.keys < end.keys)
    {
      j = key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

      while (j != i)
      {
        k = key_radix_key2class(key_purify(*parts[j].keys), rhigh, bit_mask);

        if (k != j) elem_xchange(&xi, &parts[j], sx);

        elem_inc(&parts[j]);

        j = k;
      }

      elem_inc(&xi);
    }
  }

  --rhigh;

  if (rhigh >= rlow)
  {
#ifdef SR_IP_INSERTSORT
    bit_mask = 0;
    if (rhigh - rlow + 1 <= key_radix_high) bit_mask = z_powof2_typed(rhigh - rlow + 1, slkey_pure_t);
    bit_mask = (bit_mask - 1) << rlow;
#endif

    elem_assign(s, &xi);
    for (i = 0; i < nclasses; i++)
    {
      xi.size = c[i];

#if defined(SR_IP_INSERTSORT)
      if (xi.size > SL_DEFCON(sr.ip_threshold)) rs_rec_ip(&xi, sx, rhigh, rlow, rwidth);
      else
      {
        if (xi.size > 1)
# ifdef SR_IP_INSERTSORT
          sort_insert_bmask_kernel(&xi, sx, bit_mask);
# endif
      }
#else
      if (xi.size > 1) rs_rec_ip(&xi, sx, rhigh, rlow, rwidth);
#endif

      elem_add(&xi, c[i]);
    }
  }
}


slint_t sort_radix_ip(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth) /* sl_proto, sl_func sort_radix_ip */
{
  elements_t _sx;


  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_radix);

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);

  } else if (sx->size < 1) return -1;

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;
  if (rwidth <= 0) rwidth = sort_radix_ip_width_default;

  rs_rec_ip(s, sx, rhigh, rlow, z_min(rwidth, sort_radix_width_max));

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_radix);

  return 0;
}


static void rs_rec_db(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t switchdb)
{
  slkey_pure_t bit_mask, nclasses;

  slint_t i, j, current_width, c[max_nclasses(sort_radix_width_max, slkey_pure_t)];
  elements_t xi, xj, end, parts[max_nclasses(sort_radix_width_max, slkey_pure_t)];

  elem_assign_at(s, s->size, &end);

  current_width = z_min(rwidth, rhigh - rlow + 1);
  rhigh -= current_width - 1;

  nclasses = z_powof2_typed(current_width, slkey_pure_t);
  bit_mask = nclasses - 1;


  /* zero all counter */
  for (i = 0; i < nclasses; i++) c[i] = 0;

  /* count the number of elements in every class */
  for (elem_assign(s, &xi); xi.keys < end.keys; elem_inc(&xi)) ++c[key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask)];

  /* compute the target of every class */
  elem_assign(sx, &parts[0]);
  for (i = 1; i < nclasses; i++) elem_assign_at(&parts[i - 1], c[i - 1], &parts[i]);

  /* split the elements */
  elem_assign(s, &xi);
  elem_assign_at(s, s->size, &end);
  while (xi.keys < end.keys)
  {
    j = key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

    elem_copy(&xi, &parts[j]);

    elem_inc(&xi);
    elem_inc(&parts[j]);
  }

  --rhigh;

  if (rhigh >= rlow)
  {
#ifdef SR_DB_INSERTSORT
    bit_mask = 0;
    if (rhigh - rlow + 1 <= key_radix_high) bit_mask = z_powof2_typed(rhigh - rlow + 1, slkey_pure_t);
    bit_mask = (bit_mask - 1) << rlow;
#endif

    elem_assign(s, &xi);
    elem_assign(sx, &xj);
    for (i = 0; i < nclasses; i++)
    {
      xi.size = xj.size = c[i];

#ifdef SR_DB_INSERTSORT
      if (c[i] > SL_DEFCON(sr.db_threshold)) rs_rec_db(&xj, &xi, rhigh, rlow, rwidth, (!switchdb));
      else
      {
        if (c[i] > 1) sort_insert_bmask_kernel(&xj, &xi, bit_mask);
        if (switchdb) elem_ncopy(&xj, &xi, c[i]);
      }

      elem_add(&xi, c[i]);
      elem_add(&xj, c[i]);
#else
      if (c[i] > 1) rs_rec_db(&xj, &xi, rhigh, rlow, rwidth, (!switchdb));
#endif
    }

  } else elem_ncopy(sx, s, s->size);
}


slint_t sort_radix_db(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth) /* sl_proto, sl_func sort_radix_db */
{
  elements_t _sx;


  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_radix);

  if (sx == NULL)
  {
    sx = &_sx;
    elements_alloc(sx, s->size, SLCM_ALL);

  } else if (sx->size < s->size) return -1;

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;
  if (rwidth <= 0) rwidth = sort_radix_db_width_default;

  rs_rec_db(s, sx, rhigh, rlow, z_min(rwidth, sort_radix_width_max), 1);

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_radix);

  return 0;
}


static void rs_rec_ma_db(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t switchdb)
{
  slkey_pure_t bit_mask, nclasses;

  slint_t i, j, current_width, c[max_nclasses(sort_radix_width_max, slkey_pure_t)];
  elements_t xi, xj, end, parts[max_nclasses(sort_radix_width_max, slkey_pure_t)];

  elem_assign_at(s, s->size, &end);

  current_width = z_min(rwidth, rhigh - rlow + 1);
  rhigh -= current_width - 1;

  nclasses = z_powof2_typed(current_width, slkey_pure_t);
  bit_mask = nclasses - 1;


  /* zero all counter */
  for (i = 0; i < nclasses; i++) c[i] = 0;

  /* count the number of elements in every class */
  for (elem_assign(s, &xi); xi.keys < end.keys; elem_inc(&xi)) ++c[key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask)];

  /* compute the target of every class */
  elem_assign(sx, &parts[0]);
  for (i = 1; i < nclasses; i++) elem_assign_at(&parts[i - 1], c[i - 1], &parts[i]);

  /* split the elements */
  elem_assign(s, &xi);
  elem_assign_at(s, s->size, &end);
  while (xi.keys < end.keys)
  {
    j = key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

    elem_copy(&xi, &parts[j]);

    elem_inc(&xi);
    elem_inc(&parts[j]);
  }

  --rhigh;

  if (rhigh >= rlow)
  {
#ifdef SR_MA_INSERTSORT
    bit_mask = 0;
    if (rhigh - rlow + 1 <= key_radix_high) bit_mask = z_powof2_typed(rhigh - rlow + 1, slkey_pure_t);
    bit_mask = (bit_mask - 1) << rlow;
#endif

    elem_assign(s, &xi);
    elem_assign(sx, &xj);
    for (i = 0; i < nclasses; i++)
    {
      xi.size = xj.size = c[i];

#ifdef SR_MA_INSERTSORT
      if (c[i] > SL_DEFCON(sr.ma_threshold)) rs_rec_ma_db(&xj, &xi, rhigh, rlow, rwidth, (!switchdb));
      else
      {
        if (c[i] > 1) sort_insert_bmask_kernel(&xj, &xi, bit_mask);
        if (switchdb) elem_ncopy(&xj, &xi, c[i]);
      }

      elem_add(&xi, c[i]);
      elem_add(&xj, c[i]);
#else
      if (c[i] > 1) rs_rec_ma_db(&xj, &xi, rhigh, rlow, rwidth, (!switchdb));
#endif
    }

  } else elem_ncopy(sx, s, s->size);
}


static void rs_rec_ma(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth)
{
  slkey_pure_t bit_mask, nclasses;

  slint_t i, j, k, current_width, c[max_nclasses(sort_radix_width_max, slkey_pure_t)];
  elements_t xi, end, parts[max_nclasses(sort_radix_width_max, slkey_pure_t)];

  elem_assign_at(s, s->size, &end);

  current_width = z_min(rwidth, rhigh - rlow + 1);
  rhigh -= current_width - 1;

  nclasses = z_powof2_typed(current_width, slkey_pure_t);
  bit_mask = nclasses - 1;


  /* zero all counter */
  for (i = 0; i < nclasses; i++) c[i] = 0;

  /* count the number of elements in every class */
  for (elem_assign(s, &xi); xi.keys < end.keys; elem_inc(&xi)) ++c[key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask)];

  /* compute the target of every class */
  elem_assign(s, &parts[0]);
  for (i = 1; i < nclasses; i++) elem_assign_at(&parts[i - 1], c[i - 1], &parts[i]);;

  /* split the elements */
  elem_assign(s, &end);
  for (i = 0; i < nclasses; i++)
  {
    elem_add(&end, c[i]);

    elem_assign(&parts[i], &xi);

    while (xi.keys < end.keys)
    {
      j = key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

      while (j != i)
      {
        k = key_radix_key2class(key_purify(*parts[j].keys), rhigh, bit_mask);

        if (k != j) elem_xchange(&xi, &parts[j], sx);

        elem_inc(&parts[j]);

        j = k;
      }

      elem_inc(&xi);
    }
  }

  --rhigh;

  if (rhigh >= rlow)
  {
#ifdef SR_MA_INSERTSORT
    bit_mask = 0;
    if (rhigh - rlow + 1 <= key_radix_high) bit_mask = z_powof2_typed(rhigh - rlow + 1, slkey_pure_t);
    bit_mask = (bit_mask - 1) << rlow;
#endif

    elem_assign(s, &xi);
    for (i = 0; i < nclasses; i++)
    {
      xi.size = c[i];

#ifdef SR_MA_INSERTSORT
      if (xi.size > SL_DEFCON(sr.ma_threshold))
#else
      if (xi.size > 1)
#endif
      {
        if (xi.size > sx->size) rs_rec_ma(&xi, sx, rhigh, rlow, rwidth);
        else rs_rec_ma_db(&xi, sx, rhigh, rlow, rwidth, 1);
      }
#ifdef SR_MA_INSERTSORT
        else
      {
        if (xi.size > 1) sort_insert_bmask_kernel(&xi, sx, bit_mask);
      }
#endif

      elem_add(&xi, c[i]);
    }
  }
}


slint_t sort_radix_ma(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth) /* sl_proto, sl_func sort_radix_ma */
{
  elements_t _sx;


  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_radix);

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);

  } else if (sx->size < 1) return -1;

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;
  if (rwidth <= 0) rwidth = sort_radix_width_default;

  rs_rec_ma(s, sx, rhigh, rlow, z_min(rwidth, sort_radix_width_max));

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_radix);

  return 0;
}


slint_t sort_radix(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth) /* sl_proto, sl_func sort_radix */
{
  if (sx && sx->size >= s->size) return sort_radix_db(s, sx, rhigh, rlow, rwidth);
  else return sort_radix_ip(s, sx, rhigh, rlow, rwidth);
}


#undef max_nclasses


#endif /* key_integer */


/* simple MSDF 1-bit radix-sort */

#include "sl_common.h"


#ifdef key_integer


static slint_t split2_b_1brs(elements_t *s, elements_t *sx, slkey_pure_t bmask)
{
  elements_t xl, xh;

  elem_assign(s, &xl);
  elem_assign_at(s, s->size - 1, &xh);

  while (1)
  {
    while (xl.keys < xh.keys)
    if (key_purify(*xl.keys) & bmask) break; else elem_inc(&xl);

    while (xl.keys < xh.keys)
    if (key_purify(*xh.keys) & bmask) elem_dec(&xh); else break;

    if (xl.keys >= xh.keys) break;

    elem_xchange(&xl, &xh, sx);
    elem_inc(&xl);
    elem_dec(&xh);
  }

  return xl.keys - s->keys + ((key_purify(*xl.keys) & bmask) == 0);
}


slint_t sort_radix_1bit_kernel(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow) /* sl_proto, sl_func sort_radix_1bit_kernel */
{
  slkey_pure_t bmask;

  elements_t xl, xh;

  slint_t n0, n1;

  elem_assign(s, &xl);

  while (xl.size > 1)
  {
    bmask = z_powof2_typed(rhigh, slkey_pure_t);

    n0 = split2_b_1brs(&xl, sx, bmask);
    n1 = xl.size - n0;

    if (rhigh <= rlow) break;

    rhigh--;

    xl.size = n0;
    
    if (n0 <= n1)
    {
      sort_radix_1bit_kernel(&xl, sx, rhigh, rlow);

      elem_add(&xl, n0);
      xl.size = n1;

    } else
    {
      elem_assign_at(&xl, n0, &xh);
      xh.size = n1;

      sort_radix_1bit_kernel(&xh, sx, rhigh, rlow);
    }
  }

  return 0;
}


slint sort_radix_1bit(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow) /* sl_proto, sl_func sort_radix_1bit */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);
  }

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;

  sort_radix_1bit_kernel(s, sx, rhigh, rlow);

  if (sx == &_sx) elements_free(sx);

  return 0;
}


#endif /* key_integer */



/* sl_macro SRI_INSERTSORT */
#define SRI_INSERTSORT


#include "sl_common.h"

/* sl_context CONTEXT_BEGIN sri */
const slint_t default_sri_threshold = sort_radix_iter_threshold;  /* sl_global sl_context sl_var default_sri_threshold */
/* sl_context CONTEXT_END sri */


#ifdef key_integer


#define max_nclasses (z_powof2_typed(sort_radix_width_max, slkey_pure_t))


/* insert sort */
static void rs_iter_insertsort(elements_t *s, elements_t *sx, slint_t presorted, slint_t rhigh, slint_t rlow)
{
  slint_t i, j;

  slkey_pure_t class_mask = 0, bmask = 0;


  if (presorted && (rhigh < key_radix_high)) class_mask = ~(z_powof2_typed(rhigh + 1, slkey_pure_t) - 1);

  if (rhigh - rlow + 1 <= key_radix_high) bmask = z_powof2_typed(rhigh - rlow + 1, slkey_pure_t);

  bmask = class_mask | (bmask - 1) << rlow;

  for (i = 1; i < s->size; i++)
  {
    if (key_pure_cmp_lt(key_purify(s->keys[i]) & bmask, key_purify(s->keys[i - 1]) & bmask))
    {
      j = i - 1;
      elem_copy_at(s, i, sx, 0);

      do
      {
        elem_copy_at(s, j, s, j + 1);
        if (--j < 0) break;

      } while (key_pure_cmp_lt(key_purify(*sx->keys) & bmask, key_purify(s->keys[j]) & bmask));

      elem_copy_at(sx, 0, s, j + 1);
    }
  }
}


static void rs_iter(elements_t *s, elements_t *sx, slint_t presorted, slint_t rhigh, slint_t rlow, slint_t rwidth)
{
  slkey_pure_t class_mask, bit_mask, current_class, nclasses;

  slint_t i, j, k, current_width, c[max_nclasses];
  elements_t xi, xj, end, parts[max_nclasses];
  
  slint_t rhigh_old = rhigh;


  if (presorted && (rhigh < key_radix_high)) class_mask = ~(z_powof2_typed(rhigh + 1, slkey_pure_t) - 1);
  else class_mask = 0;

  elem_assign_at(s, s->size, &end);

  while (rhigh >= rlow)
  {
    current_width = z_min(rwidth, rhigh - rlow + 1);
    rhigh -= current_width - 1;

    nclasses = z_powof2_typed(current_width, slkey_pure_t);
    bit_mask = nclasses - 1;

    elem_assign(s, &xi);

    while (xi.keys < end.keys)
    {
      elem_assign(&xi, &xj);

      current_class = class_mask & key_purify(*xi.keys);

      /* zero all counters */
      for (i = 0; i < nclasses; ++i) c[i] = 0;

      /* counting the number of sub-class elements in the current class */
      do
      {
        i = key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

        ++c[i];

        elem_inc(&xi);

        /* stop if we are at the end */
        if (xi.keys >= end.keys) break;

      } while (current_class == (class_mask & key_purify(*xi.keys)));

      /* if the current class has enough elements, perform split */
      if (xi.keys - xj.keys > 
#ifdef SRI_INSERTSORT
        SL_DEFCON(sri.threshold)
#else
        1
#endif
        )
      {
        /* compute the target of every sub-class */
        elem_assign(&xj, &parts[0]);
        for (i = 1; i < nclasses; i++) elem_assign_at(&parts[i - 1], c[i - 1], &parts[i]);;

        for (i = 0; i < nclasses; i++)
        {
          elem_add(&xj, c[i]);
          elem_assign(&parts[i], &xi);

          while (xi.keys < xj.keys)
          {
            j = key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

            while (j != i)
            {
              k = key_radix_key2class(key_purify(*parts[j].keys), rhigh, bit_mask);

              if (k != j) elem_xchange(&xi, &parts[j], sx);

              elem_inc(&parts[j]);

              j = k;
            }

            elem_inc(&xi);
          }
        }
      }
#ifdef SRI_INSERTSORT
        else
      {
        xj.size = xi.keys - xj.keys;

        rs_iter_insertsort(&xj, sx, presorted, rhigh_old, rlow);
      }
#endif
    }

    class_mask |= (bit_mask << rhigh);

    --rhigh;
  }
}


slint_t sort_radix_iter(elements_t *s, elements_t *sx, slint_t presorted, slint_t rhigh, slint_t rlow, slint_t rwidth) /* sl_proto, sl_func sort_radix_iter */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_radix_iter);

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);
  }

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;
  if (rwidth <= 0) rwidth = sort_radix_iter_width_default;

  rs_iter(s, sx, presorted, rhigh, rlow, z_min(rwidth, sort_radix_width_max));

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_radix_iter);

  return 0;
}


#undef max_nclasses


#endif /* key_integer */


/* a sorting_network function returns the counterpart of processor 'rank' (0..size-1) in the compare/merge-exchange operation in the 'stage' (0..)
   of the described sorting network consisting of 'size' number of participants, 'snp' may be pointing to additional preferences of the sorting-network
   - returning  < 0: there is no such 'stage' in the network
   - returning >= 0: as the rank of the counterpart, == rank or >= size means that the processor of the counterpart doesn't exist */


#include "sl_common.h"


/* simple hypercube, low -> high */
slint sn_hypercube_lh(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_hypercube_lh */
{
  slint stages = ilog2f(size);

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return rank ^ z_powof2(stage);
}

/* simple hypercube, high -> low */
slint sn_hypercube_hl(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_hypercube_hl */
{
  slint stages = ilog2f(size);

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return rank ^ z_powof2(stages - 1 - stage);
}


/* Knuth, pg. 241, 5.3.4 ex. 37 / Akl, pg. 41, 3.2 */
slint sn_odd_even_trans(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_odd_even_trans */
{
  slint stages = size;

  if (stages <= 2) --stages;

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return z_max(0, z_min(size - 1, ((stage % 2 == rank % 2)?rank + 1:rank - 1)));
}


/* single odd-stage */
slint sn_odd(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_odd */
{
  slint stages = 1;

  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return z_max(0, z_min(size - 1, ((0 == rank % 2)?rank + 1:rank - 1)));
}


/* single even-stage */
slint sn_even(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_even */
{
  slint stages = 1;

  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return z_max(0, z_min(size - 1, ((1 == rank % 2)?rank + 1:rank - 1)));
}


/* by Batcher / Knuth, pg. 225, 5.3.4 / Akl, pg. 23, 2.3 */
slint sn_batcher(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_batcher */
{
  slint stages, _p, p, q, r, d;

  _p = ilog2f(size);
  stages = (_p * (_p + 1)) / 2;

  /* if there are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  for (_p = p = z_powof2(_p - 1); p > 0; p /= 2)
  {
    q = _p;
    r = 0;
    d = p;

    do
    {
      if (stage <= 0)
      {
        if ((rank & p) == r) return rank + d;
        if ((rank - d >= 0) && (((rank - d) & p) == r)) return rank - d;
        return size;

      } else --stage;

      d = q - p;
      q /= 2;
      r = p;
    } while (d != 0); /* the original condition was: q != p before 'd = q - p' */
  }

  return -1;
}


/* by Batcher / Knuth, pg. 232, 5.3.4 / Akl, pg. 29, 2.4 */
slint sn_bitonic(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_bitonic */
{
  slint stages, p, i = 0;

  p = ilog2f(size);
  stages = (p * (p + 1)) / 2;

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  for (i = 0; (i < p) && (stage > i); stage -= ++i);

  if (up != NULL) *up = (z_powof2(i + 1) & rank);

  return rank ^ z_powof2(i - stage);
}


/* connected */
slint sn_connected(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_connected */
{
  sn_func *snfs = (sn_func *) snp;
  slint i = 0, stages = 0, cp = -1;

  while (snfs[i] != NULL) stages += (snfs[i++])(size, rank, -1, NULL, NULL);

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  i = 0;
  while ((snfs[i] != NULL) && (cp <= -1) && (stage >= 0))
  {
    cp = (snfs[i])(size, rank, stage, NULL, up);
    stage -= (snfs[i])(size, rank, -1, NULL, NULL);
    ++i;
  }

  return cp;
}



#include "sl_common.h"


/*#define DUMMY_TEST*/

#ifdef DUMMY_TEST

void dummy_reset(void *tproc_data)
{
}

int dummy_tproc(elements_t *b, slint_t x, void *tproc_data)
{
  return 0;
}

split_generic_t sg_dummy_tproc = SPLIT_GENERIC_INIT_TPROC(dummy_tproc, dummy_reset);
SPLIT_GENERIC_DEFINE_TPROC(dummy_tproc)
split_generic_t sg_ext_dummy_tproc = SPLIT_GENERIC_INIT_EXT_TPROC(dummy_tproc, dummy_reset);

int dummy_tproc_mod(elements_t *b, slint_t x, void *tproc_data, elements_t *mod)
{
  return 0;
}

split_generic_t sg_dummy_tproc_mod = SPLIT_GENERIC_INIT_TPROC_MOD(dummy_tproc_mod, dummy_reset);
SPLIT_GENERIC_DEFINE_TPROC_MOD(dummy_tproc_mod)
split_generic_t sg_ext_dummy_tproc_mod = SPLIT_GENERIC_INIT_EXT_TPROC_MOD(dummy_tproc_mod, dummy_reset);

int dummy_tprocs(elements_t *b, slint_t x, void *tproc_data, int *procs)
{
  return 0;
}

split_generic_t sg_dummy_tprocs = SPLIT_GENERIC_INIT_TPROCS(dummy_tprocs, dummy_reset);
SPLIT_GENERIC_DEFINE_TPROCS(dummy_tprocs)
split_generic_t sg_ext_dummy_tprocs = SPLIT_GENERIC_INIT_EXT_TPROCS(dummy_tprocs, dummy_reset);

int dummy_tprocs_mod(elements_t *b, slint_t x, void *tproc_data, int *procs, elements_t *mod)
{
  return 0;
}

split_generic_t sg_dummy_tprocs_mod = SPLIT_GENERIC_INIT_TPROCS_MOD(dummy_tprocs_mod, dummy_reset);
SPLIT_GENERIC_DEFINE_TPROCS_MOD(dummy_tprocs_mod)
split_generic_t sg_ext_dummy_tprocs_mod = SPLIT_GENERIC_INIT_EXT_TPROCS_MOD(dummy_tprocs_mod, dummy_reset);

#endif


static void _split_generic_count_db(elements_t *s, split_generic_t *sg, void *sg_data, int *counts, int *procs)
{
  SPEC_DECLARE_TPROC_COUNT_DB
  SPEC_DECLARE_TPROC_MOD_COUNT_DB
  SPEC_DECLARE_TPROCS_COUNT_DB
  SPEC_DECLARE_TPROCS_MOD_COUNT_DB


  switch (sg->type)
  {
    case 1:
      if (sg->tproc_count_db) sg->tproc_count_db(s, sg_data, counts);
      else SPEC_DO_TPROC_COUNT_DB(sg->tproc, sg_data, s, counts);
      break;
    case 2:
      if (sg->tproc_mod_count_db) sg->tproc_mod_count_db(s, sg_data, counts);
      else SPEC_DO_TPROC_MOD_COUNT_DB(sg->tproc_mod, sg_data, s, counts);
      break;
    case 3:
      if (sg->tprocs_count_db) sg->tprocs_count_db(s, sg_data, counts, procs);
      else SPEC_DO_TPROCS_COUNT_DB(sg->tprocs, sg_data, s, counts, procs);
      break;
    case 4:
      if (sg->tprocs_mod_count_db) sg->tprocs_mod_count_db(s, sg_data, counts, procs);
      else SPEC_DO_TPROCS_MOD_COUNT_DB(sg->tprocs_mod, sg_data, s, counts, procs);
      break;
  }
}


static void _split_generic_rearrange_db(elements_t *s, elements_t *d, split_generic_t *sg, void *sg_data, int *displs, int *procs, elements_t *m)
{
  SPEC_DECLARE_TPROC_REARRANGE_DB
  SPEC_DECLARE_TPROC_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_REARRANGE_DB
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB


  switch (sg->type)
  {
    case 1:
      if (sg->tproc_rearrange_db) sg->tproc_rearrange_db(s, d, sg_data, displs);
      else SPEC_DO_TPROC_REARRANGE_DB(sg->tproc, sg_data, s, d, displs);
      break;
    case 2:
      if (sg->tproc_mod_rearrange_db) sg->tproc_mod_rearrange_db(s, d, sg_data, displs, m);
      else SPEC_DO_TPROC_MOD_REARRANGE_DB(sg->tproc_mod, sg_data, s, d, displs, m);
      break;
    case 3:
      if (sg->tprocs_rearrange_db) sg->tprocs_rearrange_db(s, d, sg_data, displs, procs);
      else SPEC_DO_TPROCS_REARRANGE_DB(sg->tprocs, sg_data, s, d, displs, procs);
      break;
    case 4:
      if (sg->tprocs_mod_rearrange_db) sg->tprocs_mod_rearrange_db(s, d, sg_data, displs, procs, m);
      else SPEC_DO_TPROCS_MOD_REARRANGE_DB(sg->tprocs_mod, sg_data, s, d, displs, procs, m);
      break;
  }
}


static void _split_generic_count_ip(elements_t *s, split_generic_t *sg, void *sg_data, int *counts, int *procs)
{
  SPEC_DECLARE_TPROC_COUNT_IP
  SPEC_DECLARE_TPROC_MOD_COUNT_IP
  SPEC_DECLARE_TPROCS_COUNT_IP
  SPEC_DECLARE_TPROCS_MOD_COUNT_IP


  switch (sg->type)
  {
    case 1:
      if (sg->tproc_count_ip) sg->tproc_count_ip(s, sg_data, counts);
      else SPEC_DO_TPROC_COUNT_IP(sg->tproc, sg_data, s, counts);
      break;
    case 2:
      if (sg->tproc_mod_count_ip) sg->tproc_mod_count_ip(s, sg_data, counts);
      else SPEC_DO_TPROC_MOD_COUNT_IP(sg->tproc_mod, sg_data, s, counts);
      break;
    case 3:
      if (sg->tprocs_count_ip) sg->tprocs_count_ip(s, sg_data, counts, procs);
      else SPEC_DO_TPROCS_COUNT_IP(sg->tprocs, sg_data, s, counts, procs);
      break;
    case 4:
      if (sg->tprocs_mod_count_ip) sg->tprocs_mod_count_ip(s, sg_data, counts, procs);
      else SPEC_DO_TPROCS_MOD_COUNT_IP(sg->tprocs_mod, sg_data, s, counts, procs);
      break;
  }
}


static void _split_generic_rearrange_ip(elements_t *s, elements_t *x, split_generic_t *sg, void *sg_data, int *counts, int *displs, slint_t n, int *procs, elements_t *m)
{
  SPEC_DECLARE_TPROC_REARRANGE_IP
  SPEC_DECLARE_TPROC_MOD_REARRANGE_IP
  SPEC_DECLARE_TPROCS_REARRANGE_IP
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP


  switch (sg->type)
  {
    case 1:
      if (sg->tproc_rearrange_ip) sg->tproc_rearrange_ip(s, x, sg_data, displs, counts, n);
      else SPEC_DO_TPROC_REARRANGE_IP(sg->tproc, sg_data, s, x, displs, counts, n);
      break;
    case 2:
      if (sg->tproc_mod_rearrange_ip) sg->tproc_mod_rearrange_ip(s, x, sg_data, displs, counts, n, m);
      else SPEC_DO_TPROC_MOD_REARRANGE_IP(sg->tproc_mod, sg_data, s, x, displs, counts, n, m);
      break;
    case 3:
      if (sg->tprocs_rearrange_ip) sg->tprocs_rearrange_ip(s, x, sg_data, displs, counts, n, procs);
      else SPEC_DO_TPROCS_REARRANGE_IP(sg->tprocs, sg_data, s, x, displs, counts, n, procs);
      break;
    case 4:
      if (sg->tprocs_mod_rearrange_ip) sg->tprocs_mod_rearrange_ip(s, x, sg_data, displs, counts, n, procs, m);
      else SPEC_DO_TPROCS_MOD_REARRANGE_IP(sg->tprocs_mod, sg_data, s, x, displs, counts, n, procs, m);
      break;
  }
}


slint_t split_generic_db(elements_t *s, elements_t *d, split_generic_t *sg, void *sg_data, slint_t n) /* sl_proto, sl_func split_generic_db */
{
  slint_t i;

  int *procs = NULL;
  int *countsdispls;

  slint_t tprocs_max;
  elements_t _m, *m;
  
  
  countsdispls = z_alloc(n, sizeof(int));

  tprocs_max = 0;

  if (sg->type == 3) tprocs_max = sg->tprocs(s, -1, sg_data, NULL);
  else if (sg->type == 4) tprocs_max = sg->tprocs_mod(s, -1, sg_data, NULL, NULL);

  m = &_m;

  if (tprocs_max < 0)
  {
    tprocs_max *= -1;
    m = NULL;
  }

  if (tprocs_max) procs = z_alloca(tprocs_max, sizeof(int));

  if (sg->type == 2 || sg->type == 4) elements_alloc(&_m, z_max(1, tprocs_max), SLCM_ALL);

  if (sg->reset) sg->reset(sg_data);

  for (i = 0; i < n; ++i) countsdispls[i] = 0;

  _split_generic_count_db(s, sg, sg_data, countsdispls, procs);

  counts2displs(n, countsdispls, NULL);

  _split_generic_rearrange_db(s, d, sg, sg_data, countsdispls, procs, m);

  if (tprocs_max) z_freea(procs);

  if (sg->type == 2 || sg->type == 4) elements_free(&_m);

  z_free(countsdispls);

  return 0;
}


slint_t split_generic_ip(elements_t *s, elements_t *d, split_generic_t *sg, void *sg_data, slint_t n) /* sl_proto, sl_func split_generic_ip */
{
  slint_t i;

  int *procs = NULL;
  int *counts, *displs;

  slint_t tprocs_max;
  elements_t _m, *m;
  
  
  counts = z_alloc(2 * n, sizeof(int));
  displs = counts + n;

  tprocs_max = 0;

  if (sg->type == 3) tprocs_max = sg->tprocs(s, -1, sg_data, NULL);
  else if (sg->type == 4) tprocs_max = sg->tprocs_mod(s, -1, sg_data, NULL, NULL);

  m = &_m;

  if (tprocs_max < 0)
  {
    tprocs_max *= -1;
    m = NULL;
  }

  if (tprocs_max) procs = z_alloca(tprocs_max, sizeof(int));

  if (sg->type == 2 || sg->type == 4) elements_alloc(&_m, z_max(1, tprocs_max), SLCM_ALL);

  if (sg->reset) sg->reset(sg_data);

  for (i = 0; i < n; ++i) counts[i] = 0;

  _split_generic_count_ip(s, sg, sg_data, counts, procs);

  counts2displs(n, counts, displs);

  _split_generic_rearrange_ip(s, d, sg, sg_data, counts, displs, n, procs, m);

  if (tprocs_max) z_freea(procs);

  if (sg->type == 2 || sg->type == 4) elements_free(&_m);

  z_free(counts);

  return 0;
}


slint_t split_generic_count_db(elements_t *s, split_generic_t *sg, void *sg_data, int *counts, slint_t n) /* sl_proto, sl_func split_generic_count_db */
{
  slint_t i;

  int *procs = NULL;

  slint_t tprocs_max;
  
  
  tprocs_max = 0;

  if (sg->type == 3) tprocs_max = sg->tprocs(s, -1, sg_data, NULL);
  else if (sg->type == 4) tprocs_max = sg->tprocs_mod(s, -1, sg_data, NULL, NULL);

  if (tprocs_max < 0) tprocs_max *= -1;

  if (tprocs_max) procs = z_alloca(tprocs_max, sizeof(int));

  if (sg->reset) sg->reset(sg_data);

  for (i = 0; i < n; ++i) counts[i] = 0;

  _split_generic_count_db(s, sg, sg_data, counts, procs);

  if (tprocs_max) z_freea(procs);

  return 0;
}


slint_t split_generic_count_ip(elements_t *s, split_generic_t *sg, void *sg_data, int *counts, slint_t n) /* sl_proto, sl_func split_generic_count_ip */
{
  slint_t i;

  int *procs = NULL;

  slint_t tprocs_max;
  
  
  tprocs_max = 0;

  if (sg->type == 3) tprocs_max = sg->tprocs(s, -1, sg_data, NULL);
  else if (sg->type == 4) tprocs_max = sg->tprocs_mod(s, -1, sg_data, NULL, NULL);

  if (tprocs_max < 0) tprocs_max *= -1;

  if (tprocs_max) procs = z_alloca(tprocs_max, sizeof(int));

  if (sg->reset) sg->reset(sg_data);

  for (i = 0; i < n; ++i) counts[i] = 0;

  _split_generic_count_ip(s, sg, sg_data, counts, procs);

  if (tprocs_max) z_freea(procs);

  return 0;
}


slint_t split_generic_rearrange_db(elements_t *s, elements_t *d, split_generic_t *sg, void *sg_data, int *counts, slint_t n) /* sl_proto, sl_func split_generic_rearrange_db */
{
  int *procs = NULL;

  slint_t tprocs_max;
  elements_t _m, *m;
  
  
  tprocs_max = 0;

  if (sg->type == 3) tprocs_max = sg->tprocs(s, -1, sg_data, NULL);
  else if (sg->type == 4) tprocs_max = sg->tprocs_mod(s, -1, sg_data, NULL, NULL);

  m = &_m;

  if (tprocs_max < 0)
  {
    tprocs_max *= -1;
    m = NULL;
  }

  if (tprocs_max) procs = z_alloca(tprocs_max, sizeof(int));

  if (sg->type == 2 || sg->type == 4) elements_alloc(&_m, z_max(1, tprocs_max), SLCM_ALL);

  if (sg->reset) sg->reset(sg_data);

  counts2displs(n, counts, NULL);

  _split_generic_rearrange_db(s, d, sg, sg_data, counts, procs, m);

  if (tprocs_max) z_freea(procs);

  if (sg->type == 2 || sg->type == 4) elements_free(&_m);

  return 0;
}


slint_t split_generic_rearrange_ip(elements_t *s, elements_t *d, split_generic_t *sg, void *sg_data, int *counts, int *displs, slint_t n) /* sl_proto, sl_func split_generic_rearrange_ip */
{
  int *procs = NULL;

  slint_t tprocs_max;
  elements_t _m, *m;
  
  
  tprocs_max = 0;

  if (sg->type == 3) tprocs_max = sg->tprocs(s, -1, sg_data, NULL);
  else if (sg->type == 4) tprocs_max = sg->tprocs_mod(s, -1, sg_data, NULL, NULL);

  m = &_m;

  if (tprocs_max < 0)
  {
    tprocs_max *= -1;
    m = NULL;
  }

  if (tprocs_max) procs = z_alloca(tprocs_max, sizeof(int));

  if (sg->type == 2 || sg->type == 4) elements_alloc(&_m, z_max(1, tprocs_max), SLCM_ALL);

  if (sg->reset) sg->reset(sg_data);

  _split_generic_rearrange_ip(s, d, sg, sg_data, counts, displs, n, procs, m);

  if (tprocs_max) z_freea(procs);

  if (sg->type == 2 || sg->type == 4) elements_free(&_m);

  return 0;
}



#include "sl_common.h"


slint_t splitter_reset(splitter_t *sp) /* sl_proto, sl_func splitter_reset */
{
  slint_t i;


  for (i = 0; i < sp->n; ++i) sp->displs[i] = 0;

  return 0;
}



#include "sl_common.h"


slint_t splitx_radix(elements_t *s, elements_t *sx, slint_t nclasses, slint_t shl, slint_t *counts) /* sl_proto, sl_func splitx_radix */
{
#define max_nclasses  nclasses

  slkey_pure_t bit_mask = nclasses - 1;

  slint_t i, j, k, c[max_nclasses];
  elements_t xi, end, parts[max_nclasses];

  elem_assign_at(s, s->size, &end);
  
  
  if (counts == NULL)
  {
    counts = c;
  
    for (i = 0; i < nclasses; i++) counts[i] = 0;

    for (elem_assign(s, &xi); xi.keys < end.keys; elem_inc(&xi)) ++counts[key_radix_key2class(key_purify(*xi.keys), shl, bit_mask)];
  }

  elem_assign(s, &parts[0]);
  for (i = 1; i < nclasses; i++) elem_assign_at(&parts[i - 1], counts[i - 1], &parts[i]);;

  elem_assign(s, &end);
  for (i = 0; i < nclasses; i++)
  {
    elem_add(&end, counts[i]);

    elem_assign(&parts[i], &xi);

    while (xi.keys < end.keys)
    {
      j = key_radix_key2class(key_purify(*xi.keys), shl, bit_mask);

      while (j != i)
      {
        k = key_radix_key2class(key_purify(*parts[j].keys), shl, bit_mask);

        if (k != j) elem_xchange(&xi, &parts[j], sx);

        elem_inc(&parts[j]);

        j = k;
      }

      elem_inc(&xi);
    }
  }

  return 0;
}

#undef max_nclasses


slint split2_lt_ge(elements_t *s, slkey_pure_t *k, elements_t *t) /* sl_proto, sl_func split2_lt_ge */
{
  elements_t low, high;

  elem_assign(s, &low);
  elem_assign_at(s, s->size - 1, &high);

  while (1)
  {
    while (low.keys < high.keys)
    if (key_pure_cmp_lt(key_purify(*low.keys), *k)) elem_inc(&low); else break;

    while (low.keys < high.keys)
    if (key_pure_cmp_ge(key_purify(*high.keys), *k)) elem_dec(&high); else break;

    if (low.keys >= high.keys) break;

    elem_copy(&low, t);
    elem_copy(&high, &low);
    elem_copy(t, &high);
    elem_inc(&low);
    elem_dec(&high);
  }

  return (low.keys - s->keys) + (key_pure_cmp_lt(key_purify(*low.keys), *k));
}


slint split2_le_gt(elements_t *s, slkey_pure_t *k, elements_t *t) /* sl_proto, sl_func split2_le_gt */
{
  elements_t low, high;

  elem_assign(s, &low);
  elem_assign_at(s, s->size - 1, &high);

  while (1)
  {
    while (low.keys < high.keys)
    if (key_pure_cmp_le(key_purify(*low.keys), *k)) elem_inc(&low); else break;

    while (low.keys < high.keys)
    if (key_pure_cmp_gt(key_purify(*high.keys), *k)) elem_dec(&high); else break;

    if (low.keys >= high.keys) break;

    elem_copy(&low, t);
    elem_copy(&high, &low);
    elem_copy(t, &high);
    elem_inc(&low);
    elem_dec(&high);
  }

  return (low.keys - s->keys) + (key_pure_cmp_le(key_purify(*low.keys), *k));
}


slint split3_lt_eq_gt(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle) /* sl_proto, sl_func split3_lt_eq_gt */
{
  elements_t low, high, middle;

/*  slint i;*/

  *nlt = *nle = 0;

  if (s->size <= 0) return 0;

  elem_assign(s, &low);
  elem_assign(s, &middle);
  elem_assign_at(s, s->size, &high);

  while (low.keys < high.keys)
  if (key_pure_cmp_eq(key_purify(*low.keys), *k))
  {
    elem_xchange(&middle, &low, t);
    elem_inc(&middle);
    elem_inc(&low);

  } else if (key_pure_cmp_lt(key_purify(*low.keys), *k))
  {
    elem_inc(&low);

  } else
  {
    elem_dec(&high);
    elem_xchange(&high, &low, t);
  }

  *nle = low.keys - s->keys;
  *nlt = *nle - (middle.keys - s->keys);

  return 0;
}


slint split3_lt_eq_gt_old(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle) /* sl_proto, sl_func split3_lt_eq_gt_old */
{
  elements_t low, high, middle;

/*  slint i;*/

  *nlt = *nle = 0;

  if (s->size <= 0) return 0;

  elem_assign(s, &low);
  elem_assign(s, &middle);
  elem_assign_at(s, s->size - 1, &high);

  while (1)
  {
    while (low.keys < high.keys)
    if (key_pure_cmp_lt(key_purify(*low.keys), *k)) elem_inc(&low);
    else if (key_pure_cmp_eq(key_purify(*low.keys), *k))
    {
      elem_copy(&middle, t);
      elem_copy(&low, &middle);
      elem_copy(t, &low);
      elem_inc(&low);
      elem_inc(&middle);
    } else break;

    while (low.keys < high.keys)
    if (key_pure_cmp_gt(key_purify(*high.keys), *k)) elem_dec(&high); else break;

    if (low.keys >= high.keys) break;

    if (key_pure_cmp_eq(key_purify(*high.keys), *k))
    {
      elem_copy(&middle, t);
      elem_copy(&high, &middle);
      elem_inc(&middle);

    } else elem_copy(&high, t);

    if (low.keys >= middle.keys)
    {
      elem_copy(&low, &high);
      elem_copy(t, &low);

    } else elem_copy(t, &high);

    elem_inc(&low);
    elem_dec(&high);
  }

  *nle = (low.keys - s->keys) + (key_pure_cmp_le(key_purify(*low.keys), *k));
  *nlt = *nle - (middle.keys - s->keys) - (key_pure_cmp_eq(key_purify(*middle.keys), *k));

  return 0;
}


slint split2_b(elements_t *s, elements_t *sx, slkey_pure_t bmask) /* sl_proto, sl_func split2_b */
{
  elements_t xl, xh;

  elem_assign(s, &xl);
  elem_assign_at(s, s->size - 1, &xh);

  while (1)
  {
    while (xl.keys < xh.keys)
    if (key_purify(*xl.keys) & bmask) break; else elem_inc(&xl);

    while (xl.keys < xh.keys)
    if (key_purify(*xh.keys) & bmask) elem_dec(&xh); else break;

    if (xl.keys >= xh.keys) break;

    elem_xchange(&xl, &xh, sx);
    elem_inc(&xl);
    elem_dec(&xh);
  }

  return xl.keys - s->keys + ((key_purify(*xl.keys) & bmask) == 0);
}


slint splitk_k2c_af(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data) /* sl_proto, sl_func splitk_k2c_af */
{
  slint i;
  elements_t xi, end;

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  elements_t *parts = NULL;
#else
  elements_t parts[k];
#endif

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  parts = z_alloca(k, sizeof(elements_t));
#endif

  /* compute the target of every class */
  elem_assign_at(s, c[0], &parts[0]);
  for (i = 1; i < k; i++) elem_assign_at(&parts[i - 1], c[i], &parts[i]);

  /* permute the keys home */
  elem_assign_at(s, s->size, &end);
  for (elem_assign(s, &xi); xi.keys < end.keys; elem_add(&xi, c[i]))
  {
    while (1)
    {
      i = (k2c)(xi.keys, 0, k2c_data);

      elem_dec(&parts[i]);

      if (xi.keys >= parts[i].keys) break;

      elem_xchange(&parts[i], &xi, sx);
    }
  }

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  z_freea(parts);
#endif

  return 0;
}


slint splitk_k2c(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data) /* sl_proto, sl_func splitk_k2c */
{
  slint i, j, l;
  elements_t xi, end;

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  elements_t *parts = NULL;
#else
  elements_t parts[k];
#endif

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  parts = z_alloca(k, sizeof(elements_t));
#endif

  /* compute the target of every class */
  elem_assign(s, &parts[0]);
  for (i = 1; i < k; i++) elem_assign_at(&parts[i - 1], c[i - 1], &parts[i]);;

  elem_assign(s, &end);
  for (i = 0; i < k; i++)
  {
    elem_add(&end, c[i]);

    elem_assign(&parts[i], &xi);

    while (xi.keys < end.keys)
    {
      j = (k2c)(xi.keys, 0, k2c_data);

      while (j != i)
      {
        l = (k2c)(parts[j].keys, 0, k2c_data);

        if (l != j) elem_xchange(&xi, &parts[j], sx);

        elem_inc(&parts[j]);

        j = l;
      }

      elem_inc(&xi);
    }
  }

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  z_freea(parts);
#endif

  return 0;
}


slint splitk_k2c_count(elements_t *s, slint k, slint *c, k2c_func k2c, void *k2c_data) /* sl_proto, sl_func splitk_k2c_count */
{
  slint i;

  elements_t xi, end;

  elem_assign(s, &xi);
  elem_assign_at(s, s->size, &end);

  for (i = 0; i < k; i++) c[i] = 0;

  while (xi.keys < end.keys)
  {
    i = (k2c)(xi.keys, 0, k2c_data);

    c[i]++;

    elem_inc(&xi);
  }

  return 0;
}
