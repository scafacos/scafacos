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


#ifndef __SL_PROTOS_H__
#define __SL_PROTOS_H__


/* base/base.c */
slint_t SL_PROTO(binning_create)(local_bins_t *lb, slint_t max_nbins, slint_t max_nbinnings, elements_t *s, slint_t nelements, slint_t docounts, slint_t doweights, binning_t *bm);
slint_t SL_PROTO(binning_destroy)(local_bins_t *lb);
slint_t SL_PROTO(binning_pre)(local_bins_t *lb);
slint_t SL_PROTO(binning_exec_reset)(local_bins_t *lb, slint_t do_bins, slint_t do_prefixes);
slint_t SL_PROTO(binning_exec)(local_bins_t *lb, slint_t b, slint_t do_bins, slint_t do_prefixes);
slint_t SL_PROTO(binning_refine)(local_bins_t *lb, slint_t b, slint_t k, splitter_t *sp, slint_t s);
slint_t SL_PROTO(binning_hit)(local_bins_t *lb, slint_t b, slint_t k, splitter_t *sp, slint_t s);
slint_t SL_PROTO(binning_finalize)(local_bins_t *lb, slint_t b, slint_t dc, slweight_t dw, slint_t lc_min, slint_t lc_max, slcount_t *lcs, slweight_t *lws, splitter_t *sp, slint_t s);
slint_t SL_PROTO(binning_post)(local_bins_t *lb);
slint_t SL_PROTO(binning_radix_create)(binning_t *bm, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted);
slint_t SL_PROTO(binning_radix_destroy)(binning_t *bm);
slint_t SL_PROTO(binning_radix_pre)(binning_t *bm);
slint_t SL_PROTO(binning_radix_exec)(binning_t *bm, bin_t *bin, slcount_t *counts, slweight_t *weights);
slint_t SL_PROTO(binning_radix_refine)(binning_t *bm, bin_t *bin, slint_t k, slcount_t *counts, slweight_t *weights, splitter_t *sp, slint_t s, bin_t *new_bin);
slint_t SL_PROTO(binning_radix_hit)(binning_t *bm, bin_t *bin, slint_t k, slcount_t *counts, splitter_t *sp, slint_t s);
slint_t SL_PROTO(binning_radix_finalize)(binning_t *bm, bin_t *bin, slint_t dc, slweight_t dw, slint_t lc_min, slint_t lc_max, slcount_t *lcs, slweight_t *lws, splitter_t *sp, slint_t s);
slint_t SL_PROTO(binning_radix_post)(binning_t *bm);
slint_t SL_PROTO(elements_alloc)(elements_t *s, slint_t nelements, slcint_t components);
slint_t SL_PROTO(elements_free)(elements_t *s);
slint_t SL_PROTO(elements_realloc)(elements_t *s, slint_t nelements, slcint_t components);
slint_t SL_PROTO(elements_alloca)(elements_t *s, slint_t nelements, slcint_t components);
slint_t SL_PROTO(elements_freea)(elements_t *s);
slint_t SL_PROTO(elements_alloc_from_blocks)(elements_t *s, slint_t nblocks, void **blocks, slint_t *blocksizes, slint_t alignment, slint_t nmax, slcint_t components);
slint_t SL_PROTO(elements_alloc_from_block)(elements_t *s, void *block, slint_t blocksize, slint_t alignment, slint_t nmax, slcint_t components);
slint_t SL_PROTO(elements_block_alloc)(elements_t *s, slint_t nelements, slcint_t components);
slint_t SL_PROTO(elements_block_free)(elements_t *s);
slint_t SL_PROTO(elements_alloc_block)(elements_t *s, void **block, slint_t *blocksize, slint_t alignment, slint_t maxblocksize);
slint_t SL_PROTO(elements_copy)(elements_t *s, elements_t *d);
slint_t SL_PROTO(elements_copy_at)(elements_t *s, slint_t sat, elements_t *d, slint_t dat);
slint_t SL_PROTO(elements_ncopy)(elements_t *s, elements_t *d, slint_t n);
slint_t SL_PROTO(elements_nmove)(elements_t *s, elements_t *d, slint_t n);
slint_t SL_PROTO(elements_printf)(elements_t *s, const char *prefix);
slint_t SL_PROTO(elements_extract)(elements_t *src, slint_t nelements, elements_t *dst0, elements_t *dst1);
slint_t SL_PROTO(elements_touch)(elements_t *s);
slint_t SL_PROTO(elements_digest_sum)(elements_t *s, slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(elements_crc32)(elements_t *s, slint nelements, slint_t keys, slint_t data);
slint_t SL_PROTO(elements_digest_hash)(elements_t *s, slint_t nelements, slcint_t components, void *hash);
slint_t SL_PROTO(elements_random_exchange)(elements_t *s, slint_t rounds, elements_t *xs);
slint_t SL_PROTO(elements_keys_init_seed)(unsigned long s);
slint_t SL_PROTO(elements_keys_init)(elements_t *s, keys_init_type_t t, keys_init_data_t d);
slint_t SL_PROTO(elements_keys_init_randomized)(elements_t *s, slint_t nkeys, keys_init_type_t t, keys_init_data_t d);
slint_t SL_PROTO(elements_keys_init_from_file)(elements_t *s, slint_t data, char *filename, slint_t from, slint_t to, slint_t const_bytes_per_line);
slint_t SL_PROTO(elements_keys_save_to_file)(elements_t *s, char *filename);
slint_t SL_PROTO(elements_validate_order)(elements_t *s, slint_t n);
slint_t SL_PROTO(elements_validate_order_bmask)(elements_t *s, slint_t n, slkey_pure_t bmask);
slint_t SL_PROTO(elements_validate_order_weight)(elements_t *s, slint_t n, slkey_pure_t weight);
slint_t SL_PROTO(elements_keys_stats)(elements_t *s, slkey_pure_t *stats);
slint_t SL_PROTO(elements_keys_stats_print)(elements_t *s);
slint_t SL_PROTO(elements_print_keys)(elements_t *s);
slint_t SL_PROTO(elements_print_all)(elements_t *s);
slweight_t SL_PROTO(elements_get_weight)(elements_t *s);
slint_t SL_PROTO(elements_get_minmax_keys)(elements_t *s, slint_t nelements, slkey_pure_t *minmaxkeys);
slint_t SL_PROTO(elements_alloc_packed)(packed_elements_t *s, slint_t nelements);
slint_t SL_PROTO(elements_free_packed)(packed_elements_t *s);
slint_t SL_PROTO(elements_alloc_packed_from_block)(packed_elements_t *s, void *block, slint_t blocksize, slint_t alignment, slint_t nmax);
slint_t SL_PROTO(elements_pack_indexed)(elements_t *s, packed_elements_t *d, slindex_t *rindx, slindex_t *windx);
slint_t SL_PROTO(elements_pack)(elements_t *s, packed_elements_t *d);
slint_t SL_PROTO(elements_pack_at)(elements_t *s, slint_t sat, packed_elements_t *d, slint_t dat);
slint_t SL_PROTO(elements_unpack_indexed)(packed_elements_t *s, elements_t *d, slindex_t *rindx, slindex_t *windx);
slint_t SL_PROTO(elements_unpack)(packed_elements_t *s, elements_t *d);
slint_t SL_PROTO(elements_unpack_at)(packed_elements_t *s, slint_t sat, elements_t *d, slint_t dat);
slint_t SL_PROTO(elements_unpack_keys)(packed_elements_t *s, slkey_t *k);
slint SL_PROTO(merge2_basic_auto_01_x)(elements_t *s0, elements_t *s1, elements_t *sx);
slint SL_PROTO(merge2_basic_01_x)(elements_t *s0, elements_t *s1, elements_t *sx, m2x_func _x0_1, m2x_func _0x_1);
slint SL_PROTO(merge2_basic_01_X)(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t, m2X_func _X0_1, m2X_func _0X_1);
slint SL_PROTO(merge2_simplify_s1)(elements_t *s0, elements_t *s1, elements_t *sx, slint s1elements);
slint SL_PROTO(merge2_memory_adaptive)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_compo_hula)(elements_t *s0, elements_t *s1, elements_t *xs);
slint_t SL_PROTO(merge2_basic_sseq_x0_1)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_sseq_0x_1)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_sseq_01_x)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_sseq_01)(elements_t *s0, elements_t *s1, elements_t *t);
slint_t SL_PROTO(merge2_basic_sbin_x0_1)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_sbin_0x_1)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_sbin_01_x)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_sbin_01)(elements_t *s0, elements_t *s1, elements_t *t);
slint_t SL_PROTO(merge2_basic_shyb_x0_1)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_shyb_0x_1)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_shyb_01_x)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_shyb_01)(elements_t *s0, elements_t *s1, elements_t *t);
slint_t SL_PROTO(merge2_basic_straight_x0_1)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_straight_0x_1)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_straight_01_x)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_straight_x_0_1)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(merge2_basic_straight_X0_1)(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint_t SL_PROTO(merge2_basic_straight_0X_1)(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint_t SL_PROTO(merge2_basic_straight_01_X)(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint_t SL_PROTO(merge2_basic_straight_X0_1u)(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint_t SL_PROTO(merge2_compo_tridgell)(elements_t *s0, elements_t *s1, elements_t *sx);
slint_t SL_PROTO(mergep_2way_ip_int)(elements_t *s, elements_t *sx, slint_t p, int *displs, merge2x_f m2x);
slint_t SL_PROTO(mergep_2way_ip_int_rec)(elements_t *s, elements_t *sx, slint_t p, int *displs, merge2x_f m2x);
slint_t SL_PROTO(mergep_heap_int)(elements_t *s, elements_t *d, slint_t p, int *displs, int *counts);
slint_t SL_PROTO(mergep_heap_int_idx)(elements_t *s, elements_t *d, slint_t p, int *displs, int *counts);
slint_t SL_PROTO(mergep_heap_idx)(elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts);
slint_t SL_PROTO(mergep_heap_unpack_idx)(packed_elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts);
slint_t SL_PROTO(mergep_heap_unpack_idxonly)(packed_elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts);
slint_t SL_PROTO(permute_generic_db)(elements_t *s, elements_t *d, permute_generic_t *pg, void *pg_data);
slint_t SL_PROTO(permute_generic_ip)(elements_t *s, elements_t *x, permute_generic_t *pg, void *pg_data);
slint SL_PROTO(sl_search_sequential_lt)(elements_t *s, slpkey_t k);
slint SL_PROTO(sl_search_sequential_le)(elements_t *s, slpkey_t k);
slint SL_PROTO(sl_search_sequential_gt)(elements_t *s, slpkey_t k);
slint SL_PROTO(sl_search_sequential_ge)(elements_t *s, slpkey_t k);
slint SL_PROTO(sl_search_p_sequential_lt)(elements_t *s, slpkey_t *k);
slint SL_PROTO(sl_search_p_sequential_le)(elements_t *s, slpkey_t *k);
slint SL_PROTO(sl_search_p_sequential_gt)(elements_t *s, slpkey_t *k);
slint SL_PROTO(sl_search_p_sequential_ge)(elements_t *s, slpkey_t *k);
slint SL_PROTO(sl_search_binary_lt)(elements_t *s, slpkey_t k);
slint SL_PROTO(sl_search_binary_le)(elements_t *s, slpkey_t k);
slint SL_PROTO(sl_search_binary_gt)(elements_t *s, slpkey_t k);
slint SL_PROTO(sl_search_binary_ge)(elements_t *s, slpkey_t k);
slint SL_PROTO(sl_search_p_binary_lt)(elements_t *s, slpkey_t *k);
slint SL_PROTO(sl_search_p_binary_le)(elements_t *s, slpkey_t *k);
slint SL_PROTO(sl_search_p_binary_gt)(elements_t *s, slpkey_t *k);
slint SL_PROTO(sl_search_p_binary_ge)(elements_t *s, slpkey_t *k);
slint_t SL_PROTO(sl_search_binary_lt_bmask)(elements_t *s, slpkey_t k, slpkey_t bmask);
slint_t SL_PROTO(sl_search_binary_le_bmask)(elements_t *s, slpkey_t k, slpkey_t bmask);
slint_t SL_PROTO(sl_search_binary_sign_switch)(elements_t *s);
slint SL_PROTO(sl_search_hybrid_lt)(elements_t *s, slpkey_t k, slint t);
slint SL_PROTO(sl_search_hybrid_le)(elements_t *s, slpkey_t k, slint t);
slint SL_PROTO(sl_search_hybrid_gt)(elements_t *s, slpkey_t k, slint t);
slint SL_PROTO(sl_search_hybrid_ge)(elements_t *s, slpkey_t k, slint t);
slint SL_PROTO(sl_search_p_hybrid_lt)(elements_t *s, slpkey_t *k, slint t);
slint SL_PROTO(sl_search_p_hybrid_le)(elements_t *s, slpkey_t *k, slint t);
slint SL_PROTO(sl_search_p_hybrid_gt)(elements_t *s, slpkey_t *k, slint t);
slint SL_PROTO(sl_search_p_hybrid_ge)(elements_t *s, slpkey_t *k, slint t);
slint SL_PROTO(ilog2c)(slint x);
slint SL_PROTO(ilog2f)(slint x);
slint SL_PROTO(print_bits)(slint v);
slint SL_PROTO(pivot_random)(elements_t *s);
slint_t SL_PROTO(counts2displs)(slint_t n, int *counts, int *displs);
slint_t SL_PROTO(displs2counts)(slint_t n, int *displs, int *counts, slint_t total_counts);
void SL_PROTO(get_displcounts_extent)(slint_t n, int *displs, int *counts, slint_t *lb, slint_t *extent);
void SL_PROTO(elem_set_data)(elements_t *e, ...);
slint_t SL_PROTO(elem_get_max_byte)();
void SL_PROTO(elem_set)(elements_t *s, elements_t *d);
void SL_PROTO(elem_set_at)(elements_t *s, slint_t sat, elements_t *d);
slint_t SL_PROTO(elem_reverse)(elements_t *e, elements_t *t);
slint_t SL_PROTO(elem_nxchange_at)(elements_t *e0, slint_t at0, elements_t *e1, slint_t at1, slint_t n, elements_t *t);
slint_t SL_PROTO(elem_nxchange)(elements_t *e0, elements_t *e1, slint_t n, elements_t *t);
slint_t SL_PROTO(elem_nxchange_ro0)(elements_t *e0, elements_t *e1, slint_t n, elements_t *t);
slint_t SL_PROTO(elem_rotate)(elements_t *e, slint_t m, slint_t n, elements_t *t);
slint_t SL_PROTO(elem_rotate_ro0)(elements_t *e, slint_t m, slint_t n, elements_t *t);
slint_t SL_PROTO(elem_rotate_ro1)(elements_t *e, slint_t m, slint_t n, elements_t *t);
slint_t SL_PROTO(sort_counting_use_displs)(elements_t *s, elements_t *d, slint_t ndispls, slint_t *displs);
slint_t SL_PROTO(sort_counting_use_counts)(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts);
slint_t SL_PROTO(sort_counting_get_counts)(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts);
slint_t SL_PROTO(sort_counting)(elements_t *s, elements_t *d, slint_t ncounts);
slint SL_PROTO(sort_heap)(elements_t *s, elements_t *xs);
slint_t SL_PROTO(sort_insert_bmask_kernel)(elements_t *s, elements_t *sx, slkey_pure_t bmask);
slint_t SL_PROTO(sort_insert_kernel)(elements_t *s, elements_t *sx);
slint_t SL_PROTO(sort_insert)(elements_t *s, elements_t *sx);
slint_t SL_PROTO(sort_permute_forward)(elements_t *s, elements_t *sx, slint_t *perm, slint_t offset, slint_t mask_bit);
slint_t SL_PROTO(sort_permute_backward)(elements_t *s, elements_t *sx, slint_t *perm, slint_t offset, slint_t mask_bit);
slint SL_PROTO(sort_quick)(elements_t *s, elements_t *xs);
slint_t SL_PROTO(sort_radix_ip)(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth);
slint_t SL_PROTO(sort_radix_db)(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth);
slint_t SL_PROTO(sort_radix_ma)(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth);
slint_t SL_PROTO(sort_radix)(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth);
slint_t SL_PROTO(sort_radix_1bit_kernel)(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow);
slint SL_PROTO(sort_radix_1bit)(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow);
slint_t SL_PROTO(sort_radix_iter)(elements_t *s, elements_t *sx, slint_t presorted, slint_t rhigh, slint_t rlow, slint_t rwidth);
slint SL_PROTO(sn_hypercube_lh)(slint size, slint rank, slint stage, void *snp, slint *up);
slint SL_PROTO(sn_hypercube_hl)(slint size, slint rank, slint stage, void *snp, slint *up);
slint SL_PROTO(sn_odd_even_trans)(slint size, slint rank, slint stage, void *snp, slint *up);
slint SL_PROTO(sn_odd)(slint size, slint rank, slint stage, void *snp, slint *up);
slint SL_PROTO(sn_even)(slint size, slint rank, slint stage, void *snp, slint *up);
slint SL_PROTO(sn_batcher)(slint size, slint rank, slint stage, void *snp, slint *up);
slint SL_PROTO(sn_bitonic)(slint size, slint rank, slint stage, void *snp, slint *up);
slint SL_PROTO(sn_connected)(slint size, slint rank, slint stage, void *snp, slint *up);
slint_t SL_PROTO(split_generic_db)(elements_t *s, elements_t *d, split_generic_t *sg, void *sg_data, slint_t n);
slint_t SL_PROTO(split_generic_ip)(elements_t *s, elements_t *d, split_generic_t *sg, void *sg_data, slint_t n);
slint_t SL_PROTO(split_generic_count_db)(elements_t *s, split_generic_t *sg, void *sg_data, int *counts, slint_t n);
slint_t SL_PROTO(split_generic_count_ip)(elements_t *s, split_generic_t *sg, void *sg_data, int *counts, slint_t n);
slint_t SL_PROTO(split_generic_rearrange_db)(elements_t *s, elements_t *d, split_generic_t *sg, void *sg_data, int *counts, slint_t n);
slint_t SL_PROTO(split_generic_rearrange_ip)(elements_t *s, elements_t *d, split_generic_t *sg, void *sg_data, int *counts, int *displs, slint_t n);
slint_t SL_PROTO(splitter_reset)(splitter_t *sp);
slint_t SL_PROTO(splitx_radix)(elements_t *s, elements_t *sx, slint_t nclasses, slint_t shl, slint_t *counts);
slint SL_PROTO(split2_lt_ge)(elements_t *s, slkey_pure_t *k, elements_t *t);
slint SL_PROTO(split2_le_gt)(elements_t *s, slkey_pure_t *k, elements_t *t);
slint SL_PROTO(split3_lt_eq_gt)(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle);
slint SL_PROTO(split3_lt_eq_gt_old)(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle);
slint SL_PROTO(split2_b)(elements_t *s, elements_t *sx, slkey_pure_t bmask);
slint SL_PROTO(splitk_k2c_af)(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data);
slint SL_PROTO(splitk_k2c)(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data);
slint SL_PROTO(splitk_k2c_count)(elements_t *s, slint k, slint *c, k2c_func k2c, void *k2c_data);


#endif /* __SL_PROTOS_H__ */
