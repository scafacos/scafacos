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


#ifndef __RENAME_FMM_SORT_H__
#define __RENAME_FMM_SORT_H__


#define SL_FMM_CONCAT(_a_, _b_)   SL_FMM_CONCAT_(_a_, _b_)
#define SL_FMM_CONCAT_(_a_, _b_)  _a_##_b_

#ifdef SL_FMM_PREFIX
# define SL_FMM_FUNC(_f_)   SL_FMM_CONCAT(SL_FMM_PREFIX, _f_)
# define SL_FMM_VAR(_v_)    SL_FMM_CONCAT(SL_FMM_PREFIX, _v_)
# define SL_FMM_TYPE(_t_)   SL_FMM_CONCAT(SL_FMM_PREFIX, _t_)
# define SL_FMM_MACRO(_m_)  SL_FMM_CONCAT(SL_FMM_PREFIX, _m_)
#else
# define SL_FMM_FUNC(_f_)   _f_
# define SL_FMM_VAR(_v_)    _v_
# define SL_FMM_TYPE(_t_)   _t_
# define SL_FMM_MACRO(_m_)  _m_
#endif


/* front_xqsa0 */
#define front_xqsa0_slint_t    SL_FMM_TYPE(front_xqsa0_slint_t)
#define front_xqsa0_slint_fmt  SL_FMM_MACRO(front_xqsa0_slint_fmt)

#define front_xqsa0_slkey_t     SL_FMM_TYPE(front_xqsa0_slkey_t)
#define front_xqsa0_sldata0_t   SL_FMM_TYPE(front_xqsa0_sldata0_t)
#define front_xqsa0_sldata1_t   SL_FMM_TYPE(front_xqsa0_sldata1_t)
#define front_xqsa0_sldata2_t   SL_FMM_TYPE(front_xqsa0_sldata2_t)
#define front_xqsa0_sldata3_t   SL_FMM_TYPE(front_xqsa0_sldata3_t)
#define front_xqsa0_sldata4_t   SL_FMM_TYPE(front_xqsa0_sldata4_t)
#define front_xqsa0_slweight_t  SL_FMM_TYPE(front_xqsa0_slweight_t)

#define front_xqsa0_elements_t  SL_FMM_TYPE(front_xqsa0_elements_t)
#define front_xqsa0_merge2x_f   SL_FMM_TYPE(front_xqsa0_merge2x_f)
#define front_xqsa0_partcond_t  SL_FMM_TYPE(front_xqsa0_partcond_t)

#define front_xqsa0_elem_set_size        SL_FMM_FUNC(front_xqsa0_elem_set_size)
#define front_xqsa0_elem_set_max_size    SL_FMM_FUNC(front_xqsa0_elem_set_max_size)
#define front_xqsa0_elem_set_keys        SL_FMM_FUNC(front_xqsa0_elem_set_keys)
#define front_xqsa0_elem_set_data        SL_FMM_FUNC(front_xqsa0_elem_set_data)
#define front_xqsa0_elem_set_block       SL_FMM_FUNC(front_xqsa0_elem_set_block)
#define front_xqsa0_elem_set_block_size  SL_FMM_FUNC(front_xqsa0_elem_set_block_size)

#define front_xqsa0_elements_alloc_from_blocks  SL_FMM_FUNC(front_xqsa0_elements_alloc_from_blocks)
#define front_xqsa0_sort_radix_iter             SL_FMM_FUNC(front_xqsa0_sort_radix_iter)
#define front_xqsa0_sort_radix                  SL_FMM_FUNC(front_xqsa0_sort_radix)
#define front_xqsa0_merge2_compo_hula           SL_FMM_FUNC(front_xqsa0_merge2_compo_hula)
#define front_xqsa0_merge2_memory_adaptive      SL_FMM_FUNC(front_xqsa0_merge2_memory_adaptive)
#define front_xqsa0_sn_batcher                  SL_FMM_FUNC(front_xqsa0_sn_batcher)
#define front_xqsa0_mpi_datatypes_init          SL_FMM_FUNC(front_xqsa0_mpi_datatypes_init)
#define front_xqsa0_mpi_datatypes_release       SL_FMM_FUNC(front_xqsa0_mpi_datatypes_release)
#define front_xqsa0_mpi_mergek_sorted2          SL_FMM_FUNC(front_xqsa0_mpi_mergek_sorted2)
#define front_xqsa0_mpi_mergek                  SL_FMM_FUNC(front_xqsa0_mpi_mergek)
#define front_xqsa0_mpi_partition_exact_radix   SL_FMM_FUNC(front_xqsa0_mpi_partition_exact_radix)
#define front_xqsa0_sl_MPI_Alltoall_int         SL_FMM_FUNC(front_xqsa0_sl_MPI_Alltoall_int)
#define front_xqsa0_counts2displs               SL_FMM_FUNC(front_xqsa0_counts2displs)
#define front_xqsa0_mpi_elements_alltoallv_ip   SL_FMM_FUNC(front_xqsa0_mpi_elements_alltoallv_ip)
#define front_xqsa0_mpi_elements_get_weights    SL_FMM_FUNC(front_xqsa0_mpi_elements_get_weights)
#define front_xqsa0_mergep_2way_ip_int          SL_FMM_FUNC(front_xqsa0_mergep_2way_ip_int)


/* front_xqsaI */
#define front_xqsaI_slint_t    SL_FMM_TYPE(front_xqsaI_slint_t)
#define front_xqsaI_slint_fmt  SL_FMM_MACRO(front_xqsaI_slint_fmt)

#define front_xqsaI_slkey_t     SL_FMM_TYPE(front_xqsaI_slkey_t)
#define front_xqsaI_sldataI_t   SL_FMM_TYPE(front_xqsaI_sldataI_t)
#define front_xqsaI_sldata1_t   SL_FMM_TYPE(front_xqsaI_sldata1_t)
#define front_xqsaI_sldata2_t   SL_FMM_TYPE(front_xqsaI_sldata2_t)
#define front_xqsaI_sldata3_t   SL_FMM_TYPE(front_xqsaI_sldata3_t)
#define front_xqsaI_sldata4_t   SL_FMM_TYPE(front_xqsaI_sldata4_t)
#define front_xqsaI_slweight_t  SL_FMM_TYPE(front_xqsaI_slweight_t)

#define front_xqsaI_elements_t  SL_FMM_TYPE(front_xqsaI_elements_t)
#define front_xqsaI_merge2x_f   SL_FMM_TYPE(front_xqsaI_merge2x_f)
#define front_xqsaI_partcond_t  SL_FMM_TYPE(front_xqsaI_partcond_t)

#define front_xqsaI_elem_set_size        SL_FMM_FUNC(front_xqsaI_elem_set_size)
#define front_xqsaI_elem_set_max_size    SL_FMM_FUNC(front_xqsaI_elem_set_max_size)
#define front_xqsaI_elem_set_keys        SL_FMM_FUNC(front_xqsaI_elem_set_keys)
#define front_xqsaI_elem_set_data        SL_FMM_FUNC(front_xqsaI_elem_set_data)
#define front_xqsaI_elem_set_block       SL_FMM_FUNC(front_xqsaI_elem_set_block)
#define front_xqsaI_elem_set_block_size  SL_FMM_FUNC(front_xqsaI_elem_set_block_size)

#define front_xqsaI_elements_alloc_from_blocks  SL_FMM_FUNC(front_xqsaI_elements_alloc_from_blocks)
#define front_xqsaI_sort_radix_iter             SL_FMM_FUNC(front_xqsaI_sort_radix_iter)
#define front_xqsaI_sort_radix                  SL_FMM_FUNC(front_xqsaI_sort_radix)
#define front_xqsaI_merge2_compo_hula           SL_FMM_FUNC(front_xqsaI_merge2_compo_hula)
#define front_xqsaI_merge2_memory_adaptive      SL_FMM_FUNC(front_xqsaI_merge2_memory_adaptive)
#define front_xqsaI_sn_batcher                  SL_FMM_FUNC(front_xqsaI_sn_batcher)
#define front_xqsaI_mpi_datatypes_init          SL_FMM_FUNC(front_xqsaI_mpi_datatypes_init)
#define front_xqsaI_mpi_datatypes_release       SL_FMM_FUNC(front_xqsaI_mpi_datatypes_release)
#define front_xqsaI_mpi_mergek_sorted2          SL_FMM_FUNC(front_xqsaI_mpi_mergek_sorted2)
#define front_xqsaI_mpi_mergek                  SL_FMM_FUNC(front_xqsaI_mpi_mergek)
#define front_xqsaI_mpi_partition_exact_radix   SL_FMM_FUNC(front_xqsaI_mpi_partition_exact_radix)
#define front_xqsaI_sl_MPI_Alltoall_int         SL_FMM_FUNC(front_xqsaI_sl_MPI_Alltoall_int)
#define front_xqsaI_counts2displs               SL_FMM_FUNC(front_xqsaI_counts2displs)
#define front_xqsaI_mpi_elements_alltoallv_ip   SL_FMM_FUNC(front_xqsaI_mpi_elements_alltoallv_ip)
#define front_xqsaI_mpi_elements_get_weights    SL_FMM_FUNC(front_xqsaI_mpi_elements_get_weights)
#define front_xqsaI_mergep_2way_ip_int          SL_FMM_FUNC(front_xqsaI_mergep_2way_ip_int)


/* front_xqsaX */
#define front_xqsaX_slint_t    SL_FMM_TYPE(front_xqsaX_slint_t)
#define front_xqsaX_slint_fmt  SL_FMM_MACRO(front_xqsaX_slint_fmt)

#define front_xqsaX_slkey_t     SL_FMM_TYPE(front_xqsaX_slkey_t)
#define front_xqsaX_sldataX_t   SL_FMM_TYPE(front_xqsaX_sldataX_t)
#define front_xqsaX_sldata1_t   SL_FMM_TYPE(front_xqsaX_sldata1_t)
#define front_xqsaX_sldata2_t   SL_FMM_TYPE(front_xqsaX_sldata2_t)
#define front_xqsaX_sldata3_t   SL_FMM_TYPE(front_xqsaX_sldata3_t)
#define front_xqsaX_sldata4_t   SL_FMM_TYPE(front_xqsaX_sldata4_t)
#define front_xqsaX_slweight_t  SL_FMM_TYPE(front_xqsaX_slweight_t)

#define front_xqsaX_elements_t  SL_FMM_TYPE(front_xqsaX_elements_t)
#define front_xqsaX_merge2x_f   SL_FMM_TYPE(front_xqsaX_merge2x_f)
#define front_xqsaX_partcond_t  SL_FMM_TYPE(front_xqsaX_partcond_t)

#define front_xqsaX_elem_set_size        SL_FMM_FUNC(front_xqsaX_elem_set_size)
#define front_xqsaX_elem_set_max_size    SL_FMM_FUNC(front_xqsaX_elem_set_max_size)
#define front_xqsaX_elem_set_keys        SL_FMM_FUNC(front_xqsaX_elem_set_keys)
#define front_xqsaX_elem_set_data        SL_FMM_FUNC(front_xqsaX_elem_set_data)
#define front_xqsaX_elem_set_block       SL_FMM_FUNC(front_xqsaX_elem_set_block)
#define front_xqsaX_elem_set_block_size  SL_FMM_FUNC(front_xqsaX_elem_set_block_size)

#define front_xqsaX_elements_alloc_from_blocks  SL_FMM_FUNC(front_xqsaX_elements_alloc_from_blocks)
#define front_xqsaX_sort_radix_iter             SL_FMM_FUNC(front_xqsaX_sort_radix_iter)
#define front_xqsaX_sort_radix                  SL_FMM_FUNC(front_xqsaX_sort_radix)
#define front_xqsaX_merge2_compo_hula           SL_FMM_FUNC(front_xqsaX_merge2_compo_hula)
#define front_xqsaX_merge2_memory_adaptive      SL_FMM_FUNC(front_xqsaX_merge2_memory_adaptive)
#define front_xqsaX_sn_batcher                  SL_FMM_FUNC(front_xqsaX_sn_batcher)
#define front_xqsaX_mpi_datatypes_init          SL_FMM_FUNC(front_xqsaX_mpi_datatypes_init)
#define front_xqsaX_mpi_datatypes_release       SL_FMM_FUNC(front_xqsaX_mpi_datatypes_release)
#define front_xqsaX_mpi_mergek_sorted2          SL_FMM_FUNC(front_xqsaX_mpi_mergek_sorted2)
#define front_xqsaX_mpi_mergek                  SL_FMM_FUNC(front_xqsaX_mpi_mergek)
#define front_xqsaX_mpi_partition_exact_radix   SL_FMM_FUNC(front_xqsaX_mpi_partition_exact_radix)
#define front_xqsaX_sl_MPI_Alltoall_int         SL_FMM_FUNC(front_xqsaX_sl_MPI_Alltoall_int)
#define front_xqsaX_counts2displs               SL_FMM_FUNC(front_xqsaX_counts2displs)
#define front_xqsaX_mpi_elements_alltoallv_ip   SL_FMM_FUNC(front_xqsaX_mpi_elements_alltoallv_ip)
#define front_xqsaX_mpi_elements_get_weights    SL_FMM_FUNC(front_xqsaX_mpi_elements_get_weights)
#define front_xqsaX_mergep_2way_ip_int          SL_FMM_FUNC(front_xqsaX_mergep_2way_ip_int)


/* front_xqsaIl */
#define front_xqsaIl_slint_t    SL_FMM_TYPE(front_xqsaIl_slint_t)
#define front_xqsaIl_slint_fmt  SL_FMM_MACRO(front_xqsaIl_slint_fmt)

#define front_xqsaIl_slkey_t     SL_FMM_TYPE(front_xqsaIl_slkey_t)
#define front_xqsaIl_sldataIl_t   SL_FMM_TYPE(front_xqsaIl_sldataIl_t)
#define front_xqsaIl_sldata1_t   SL_FMM_TYPE(front_xqsaIl_sldata1_t)
#define front_xqsaIl_sldata2_t   SL_FMM_TYPE(front_xqsaIl_sldata2_t)
#define front_xqsaIl_sldata3_t   SL_FMM_TYPE(front_xqsaIl_sldata3_t)
#define front_xqsaIl_sldata4_t   SL_FMM_TYPE(front_xqsaIl_sldata4_t)
#define front_xqsaIl_slweight_t  SL_FMM_TYPE(front_xqsaIl_slweight_t)

#define front_xqsaIl_elements_t  SL_FMM_TYPE(front_xqsaIl_elements_t)
#define front_xqsaIl_merge2x_f   SL_FMM_TYPE(front_xqsaIl_merge2x_f)
#define front_xqsaIl_partcond_t  SL_FMM_TYPE(front_xqsaIl_partcond_t)

#define front_xqsaIl_elem_set_size        SL_FMM_FUNC(front_xqsaIl_elem_set_size)
#define front_xqsaIl_elem_set_max_size    SL_FMM_FUNC(front_xqsaIl_elem_set_max_size)
#define front_xqsaIl_elem_set_keys        SL_FMM_FUNC(front_xqsaIl_elem_set_keys)
#define front_xqsaIl_elem_set_data        SL_FMM_FUNC(front_xqsaIl_elem_set_data)
#define front_xqsaIl_elem_set_block       SL_FMM_FUNC(front_xqsaIl_elem_set_block)
#define front_xqsaIl_elem_set_block_size  SL_FMM_FUNC(front_xqsaIl_elem_set_block_size)

#define front_xqsaIl_elements_alloc_from_blocks  SL_FMM_FUNC(front_xqsaIl_elements_alloc_from_blocks)
#define front_xqsaIl_sort_radix_iter             SL_FMM_FUNC(front_xqsaIl_sort_radix_iter)
#define front_xqsaIl_sort_radix                  SL_FMM_FUNC(front_xqsaIl_sort_radix)
#define front_xqsaIl_merge2_compo_hula           SL_FMM_FUNC(front_xqsaIl_merge2_compo_hula)
#define front_xqsaIl_merge2_memory_adaptive      SL_FMM_FUNC(front_xqsaIl_merge2_memory_adaptive)
#define front_xqsaIl_sn_batcher                  SL_FMM_FUNC(front_xqsaIl_sn_batcher)
#define front_xqsaIl_mpi_datatypes_init          SL_FMM_FUNC(front_xqsaIl_mpi_datatypes_init)
#define front_xqsaIl_mpi_datatypes_release       SL_FMM_FUNC(front_xqsaIl_mpi_datatypes_release)
#define front_xqsaIl_mpi_mergek_sorted2          SL_FMM_FUNC(front_xqsaIl_mpi_mergek_sorted2)
#define front_xqsaIl_mpi_mergek                  SL_FMM_FUNC(front_xqsaIl_mpi_mergek)
#define front_xqsaIl_mpi_partition_exact_radix   SL_FMM_FUNC(front_xqsaIl_mpi_partition_exact_radix)
#define front_xqsaIl_sl_MPI_Alltoall_int         SL_FMM_FUNC(front_xqsaIl_sl_MPI_Alltoall_int)
#define front_xqsaIl_counts2displs               SL_FMM_FUNC(front_xqsaIl_counts2displs)
#define front_xqsaIl_mpi_elements_alltoallv_ip   SL_FMM_FUNC(front_xqsaIl_mpi_elements_alltoallv_ip)
#define front_xqsaIl_mpi_elements_get_weights    SL_FMM_FUNC(front_xqsaIl_mpi_elements_get_weights)
#define front_xqsaIl_mergep_2way_ip_int          SL_FMM_FUNC(front_xqsaIl_mergep_2way_ip_int)


/* front_xq_a0 */
#define front_xq_a0_slint_t    SL_FMM_TYPE(front_xq_a0_slint_t)
#define front_xq_a0_slint_fmt  SL_FMM_MACRO(front_xq_a0_slint_fmt)

#define front_xq_a0_slkey_t     SL_FMM_TYPE(front_xq_a0_slkey_t)
#define front_xq_a0_sldata0_t   SL_FMM_TYPE(front_xq_a0_sldata0_t)
#define front_xq_a0_sldata1_t   SL_FMM_TYPE(front_xq_a0_sldata1_t)
#define front_xq_a0_sldata2_t   SL_FMM_TYPE(front_xq_a0_sldata2_t)
#define front_xq_a0_sldata3_t   SL_FMM_TYPE(front_xq_a0_sldata3_t)
#define front_xq_a0_sldata4_t   SL_FMM_TYPE(front_xq_a0_sldata4_t)
#define front_xq_a0_slweight_t  SL_FMM_TYPE(front_xq_a0_slweight_t)

#define front_xq_a0_elements_t  SL_FMM_TYPE(front_xq_a0_elements_t)
#define front_xq_a0_merge2x_f   SL_FMM_TYPE(front_xq_a0_merge2x_f)
#define front_xq_a0_partcond_t  SL_FMM_TYPE(front_xq_a0_partcond_t)

#define front_xq_a0_elem_set_size        SL_FMM_FUNC(front_xq_a0_elem_set_size)
#define front_xq_a0_elem_set_max_size    SL_FMM_FUNC(front_xq_a0_elem_set_max_size)
#define front_xq_a0_elem_set_keys        SL_FMM_FUNC(front_xq_a0_elem_set_keys)
#define front_xq_a0_elem_set_data        SL_FMM_FUNC(front_xq_a0_elem_set_data)
#define front_xq_a0_elem_set_block       SL_FMM_FUNC(front_xq_a0_elem_set_block)
#define front_xq_a0_elem_set_block_size  SL_FMM_FUNC(front_xq_a0_elem_set_block_size)

#define front_xq_a0_elements_alloc_from_blocks  SL_FMM_FUNC(front_xq_a0_elements_alloc_from_blocks)
#define front_xq_a0_sort_radix_iter             SL_FMM_FUNC(front_xq_a0_sort_radix_iter)
#define front_xq_a0_sort_radix                  SL_FMM_FUNC(front_xq_a0_sort_radix)
#define front_xq_a0_merge2_compo_hula           SL_FMM_FUNC(front_xq_a0_merge2_compo_hula)
#define front_xq_a0_merge2_memory_adaptive      SL_FMM_FUNC(front_xq_a0_merge2_memory_adaptive)
#define front_xq_a0_sn_batcher                  SL_FMM_FUNC(front_xq_a0_sn_batcher)
#define front_xq_a0_mpi_datatypes_init          SL_FMM_FUNC(front_xq_a0_mpi_datatypes_init)
#define front_xq_a0_mpi_datatypes_release       SL_FMM_FUNC(front_xq_a0_mpi_datatypes_release)
#define front_xq_a0_mpi_mergek_sorted2          SL_FMM_FUNC(front_xq_a0_mpi_mergek_sorted2)
#define front_xq_a0_mpi_mergek                  SL_FMM_FUNC(front_xq_a0_mpi_mergek)
#define front_xq_a0_mpi_partition_exact_radix   SL_FMM_FUNC(front_xq_a0_mpi_partition_exact_radix)
#define front_xq_a0_sl_MPI_Alltoall_int         SL_FMM_FUNC(front_xq_a0_sl_MPI_Alltoall_int)
#define front_xq_a0_counts2displs               SL_FMM_FUNC(front_xq_a0_counts2displs)
#define front_xq_a0_mpi_elements_alltoallv_ip   SL_FMM_FUNC(front_xq_a0_mpi_elements_alltoallv_ip)
#define front_xq_a0_mpi_elements_get_weights    SL_FMM_FUNC(front_xq_a0_mpi_elements_get_weights)
#define front_xq_a0_mergep_2way_ip_int          SL_FMM_FUNC(front_xq_a0_mergep_2way_ip_int)


/* front_xq_aI */
#define front_xq_aI_slint_t    SL_FMM_TYPE(front_xq_aI_slint_t)
#define front_xq_aI_slint_fmt  SL_FMM_MACRO(front_xq_aI_slint_fmt)

#define front_xq_aI_slkey_t     SL_FMM_TYPE(front_xq_aI_slkey_t)
#define front_xq_aI_sldataI_t   SL_FMM_TYPE(front_xq_aI_sldataI_t)
#define front_xq_aI_sldata1_t   SL_FMM_TYPE(front_xq_aI_sldata1_t)
#define front_xq_aI_sldata2_t   SL_FMM_TYPE(front_xq_aI_sldata2_t)
#define front_xq_aI_sldata3_t   SL_FMM_TYPE(front_xq_aI_sldata3_t)
#define front_xq_aI_sldata4_t   SL_FMM_TYPE(front_xq_aI_sldata4_t)
#define front_xq_aI_slweight_t  SL_FMM_TYPE(front_xq_aI_slweight_t)

#define front_xq_aI_elements_t  SL_FMM_TYPE(front_xq_aI_elements_t)
#define front_xq_aI_merge2x_f   SL_FMM_TYPE(front_xq_aI_merge2x_f)
#define front_xq_aI_partcond_t  SL_FMM_TYPE(front_xq_aI_partcond_t)

#define front_xq_aI_elem_set_size        SL_FMM_FUNC(front_xq_aI_elem_set_size)
#define front_xq_aI_elem_set_max_size    SL_FMM_FUNC(front_xq_aI_elem_set_max_size)
#define front_xq_aI_elem_set_keys        SL_FMM_FUNC(front_xq_aI_elem_set_keys)
#define front_xq_aI_elem_set_data        SL_FMM_FUNC(front_xq_aI_elem_set_data)
#define front_xq_aI_elem_set_block       SL_FMM_FUNC(front_xq_aI_elem_set_block)
#define front_xq_aI_elem_set_block_size  SL_FMM_FUNC(front_xq_aI_elem_set_block_size)

#define front_xq_aI_elements_alloc_from_blocks  SL_FMM_FUNC(front_xq_aI_elements_alloc_from_blocks)
#define front_xq_aI_sort_radix_iter             SL_FMM_FUNC(front_xq_aI_sort_radix_iter)
#define front_xq_aI_sort_radix                  SL_FMM_FUNC(front_xq_aI_sort_radix)
#define front_xq_aI_merge2_compo_hula           SL_FMM_FUNC(front_xq_aI_merge2_compo_hula)
#define front_xq_aI_merge2_memory_adaptive      SL_FMM_FUNC(front_xq_aI_merge2_memory_adaptive)
#define front_xq_aI_sn_batcher                  SL_FMM_FUNC(front_xq_aI_sn_batcher)
#define front_xq_aI_mpi_datatypes_init          SL_FMM_FUNC(front_xq_aI_mpi_datatypes_init)
#define front_xq_aI_mpi_datatypes_release       SL_FMM_FUNC(front_xq_aI_mpi_datatypes_release)
#define front_xq_aI_mpi_mergek_sorted2          SL_FMM_FUNC(front_xq_aI_mpi_mergek_sorted2)
#define front_xq_aI_mpi_mergek                  SL_FMM_FUNC(front_xq_aI_mpi_mergek)
#define front_xq_aI_mpi_partition_exact_radix   SL_FMM_FUNC(front_xq_aI_mpi_partition_exact_radix)
#define front_xq_aI_sl_MPI_Alltoall_int         SL_FMM_FUNC(front_xq_aI_sl_MPI_Alltoall_int)
#define front_xq_aI_counts2displs               SL_FMM_FUNC(front_xq_aI_counts2displs)
#define front_xq_aI_mpi_elements_alltoallv_ip   SL_FMM_FUNC(front_xq_aI_mpi_elements_alltoallv_ip)
#define front_xq_aI_mpi_elements_get_weights    SL_FMM_FUNC(front_xq_aI_mpi_elements_get_weights)
#define front_xq_aI_mergep_2way_ip_int          SL_FMM_FUNC(front_xq_aI_mergep_2way_ip_int)


/* front_xq_aX */
#define front_xq_aX_slint_t    SL_FMM_TYPE(front_xq_aX_slint_t)
#define front_xq_aX_slint_fmt  SL_FMM_MACRO(front_xq_aX_slint_fmt)

#define front_xq_aX_slkey_t     SL_FMM_TYPE(front_xq_aX_slkey_t)
#define front_xq_aX_sldataX_t   SL_FMM_TYPE(front_xq_aX_sldataX_t)
#define front_xq_aX_sldata1_t   SL_FMM_TYPE(front_xq_aX_sldata1_t)
#define front_xq_aX_sldata2_t   SL_FMM_TYPE(front_xq_aX_sldata2_t)
#define front_xq_aX_sldata3_t   SL_FMM_TYPE(front_xq_aX_sldata3_t)
#define front_xq_aX_sldata4_t   SL_FMM_TYPE(front_xq_aX_sldata4_t)
#define front_xq_aX_slweight_t  SL_FMM_TYPE(front_xq_aX_slweight_t)

#define front_xq_aX_elements_t  SL_FMM_TYPE(front_xq_aX_elements_t)
#define front_xq_aX_merge2x_f   SL_FMM_TYPE(front_xq_aX_merge2x_f)
#define front_xq_aX_partcond_t  SL_FMM_TYPE(front_xq_aX_partcond_t)

#define front_xq_aX_elem_set_size        SL_FMM_FUNC(front_xq_aX_elem_set_size)
#define front_xq_aX_elem_set_max_size    SL_FMM_FUNC(front_xq_aX_elem_set_max_size)
#define front_xq_aX_elem_set_keys        SL_FMM_FUNC(front_xq_aX_elem_set_keys)
#define front_xq_aX_elem_set_data        SL_FMM_FUNC(front_xq_aX_elem_set_data)
#define front_xq_aX_elem_set_block       SL_FMM_FUNC(front_xq_aX_elem_set_block)
#define front_xq_aX_elem_set_block_size  SL_FMM_FUNC(front_xq_aX_elem_set_block_size)

#define front_xq_aX_elements_alloc_from_blocks  SL_FMM_FUNC(front_xq_aX_elements_alloc_from_blocks)
#define front_xq_aX_sort_radix_iter             SL_FMM_FUNC(front_xq_aX_sort_radix_iter)
#define front_xq_aX_sort_radix                  SL_FMM_FUNC(front_xq_aX_sort_radix)
#define front_xq_aX_merge2_compo_hula           SL_FMM_FUNC(front_xq_aX_merge2_compo_hula)
#define front_xq_aX_merge2_memory_adaptive      SL_FMM_FUNC(front_xq_aX_merge2_memory_adaptive)
#define front_xq_aX_sn_batcher                  SL_FMM_FUNC(front_xq_aX_sn_batcher)
#define front_xq_aX_mpi_datatypes_init          SL_FMM_FUNC(front_xq_aX_mpi_datatypes_init)
#define front_xq_aX_mpi_datatypes_release       SL_FMM_FUNC(front_xq_aX_mpi_datatypes_release)
#define front_xq_aX_mpi_mergek_sorted2          SL_FMM_FUNC(front_xq_aX_mpi_mergek_sorted2)
#define front_xq_aX_mpi_mergek                  SL_FMM_FUNC(front_xq_aX_mpi_mergek)
#define front_xq_aX_mpi_partition_exact_radix   SL_FMM_FUNC(front_xq_aX_mpi_partition_exact_radix)
#define front_xq_aX_sl_MPI_Alltoall_int         SL_FMM_FUNC(front_xq_aX_sl_MPI_Alltoall_int)
#define front_xq_aX_counts2displs               SL_FMM_FUNC(front_xq_aX_counts2displs)
#define front_xq_aX_mpi_elements_alltoallv_ip   SL_FMM_FUNC(front_xq_aX_mpi_elements_alltoallv_ip)
#define front_xq_aX_mpi_elements_get_weights    SL_FMM_FUNC(front_xq_aX_mpi_elements_get_weights)
#define front_xq_aX_mergep_2way_ip_int          SL_FMM_FUNC(front_xq_aX_mergep_2way_ip_int)


/* front_xq_aIl */
#define front_xq_aIl_slint_t    SL_FMM_TYPE(front_xq_aIl_slint_t)
#define front_xq_aIl_slint_fmt  SL_FMM_MACRO(front_xq_aIl_slint_fmt)

#define front_xq_aIl_slkey_t     SL_FMM_TYPE(front_xq_aIl_slkey_t)
#define front_xq_aIl_sldataIl_t   SL_FMM_TYPE(front_xq_aIl_sldataIl_t)
#define front_xq_aIl_sldata1_t   SL_FMM_TYPE(front_xq_aIl_sldata1_t)
#define front_xq_aIl_sldata2_t   SL_FMM_TYPE(front_xq_aIl_sldata2_t)
#define front_xq_aIl_sldata3_t   SL_FMM_TYPE(front_xq_aIl_sldata3_t)
#define front_xq_aIl_sldata4_t   SL_FMM_TYPE(front_xq_aIl_sldata4_t)
#define front_xq_aIl_slweight_t  SL_FMM_TYPE(front_xq_aIl_slweight_t)

#define front_xq_aIl_elements_t  SL_FMM_TYPE(front_xq_aIl_elements_t)
#define front_xq_aIl_merge2x_f   SL_FMM_TYPE(front_xq_aIl_merge2x_f)
#define front_xq_aIl_partcond_t  SL_FMM_TYPE(front_xq_aIl_partcond_t)

#define front_xq_aIl_elem_set_size        SL_FMM_FUNC(front_xq_aIl_elem_set_size)
#define front_xq_aIl_elem_set_max_size    SL_FMM_FUNC(front_xq_aIl_elem_set_max_size)
#define front_xq_aIl_elem_set_keys        SL_FMM_FUNC(front_xq_aIl_elem_set_keys)
#define front_xq_aIl_elem_set_data        SL_FMM_FUNC(front_xq_aIl_elem_set_data)
#define front_xq_aIl_elem_set_block       SL_FMM_FUNC(front_xq_aIl_elem_set_block)
#define front_xq_aIl_elem_set_block_size  SL_FMM_FUNC(front_xq_aIl_elem_set_block_size)

#define front_xq_aIl_elements_alloc_from_blocks  SL_FMM_FUNC(front_xq_aIl_elements_alloc_from_blocks)
#define front_xq_aIl_sort_radix_iter             SL_FMM_FUNC(front_xq_aIl_sort_radix_iter)
#define front_xq_aIl_sort_radix                  SL_FMM_FUNC(front_xq_aIl_sort_radix)
#define front_xq_aIl_merge2_compo_hula           SL_FMM_FUNC(front_xq_aIl_merge2_compo_hula)
#define front_xq_aIl_merge2_memory_adaptive      SL_FMM_FUNC(front_xq_aIl_merge2_memory_adaptive)
#define front_xq_aIl_sn_batcher                  SL_FMM_FUNC(front_xq_aIl_sn_batcher)
#define front_xq_aIl_mpi_datatypes_init          SL_FMM_FUNC(front_xq_aIl_mpi_datatypes_init)
#define front_xq_aIl_mpi_datatypes_release       SL_FMM_FUNC(front_xq_aIl_mpi_datatypes_release)
#define front_xq_aIl_mpi_mergek_sorted2          SL_FMM_FUNC(front_xq_aIl_mpi_mergek_sorted2)
#define front_xq_aIl_mpi_mergek                  SL_FMM_FUNC(front_xq_aIl_mpi_mergek)
#define front_xq_aIl_mpi_partition_exact_radix   SL_FMM_FUNC(front_xq_aIl_mpi_partition_exact_radix)
#define front_xq_aIl_sl_MPI_Alltoall_int         SL_FMM_FUNC(front_xq_aIl_sl_MPI_Alltoall_int)
#define front_xq_aIl_counts2displs               SL_FMM_FUNC(front_xq_aIl_counts2displs)
#define front_xq_aIl_mpi_elements_alltoallv_ip   SL_FMM_FUNC(front_xq_aIl_mpi_elements_alltoallv_ip)
#define front_xq_aIl_mpi_elements_get_weights    SL_FMM_FUNC(front_xq_aIl_mpi_elements_get_weights)
#define front_xq_aIl_mergep_2way_ip_int          SL_FMM_FUNC(front_xq_aIl_mergep_2way_ip_int)


/* back_qxpg */
#define back_qxpg_slint_t      SL_FMM_TYPE(back_qxpg_slint_t)
#define back_qxpg_slint_fmt    SL_FMM_MACRO(back_qxpg_slint_fmt)

#define back_qxpg_slkey_t    SL_FMM_TYPE(back_qxpg_slkey_t)
#define back_qxpg_sldata0_t  SL_FMM_TYPE(back_qxpg_sldata0_t)
#define back_qxpg_sldata1_t  SL_FMM_TYPE(back_qxpg_sldata1_t)
#define back_qxpg_sldata2_t  SL_FMM_TYPE(back_qxpg_sldata2_t)
#define back_qxpg_sldata3_t  SL_FMM_TYPE(back_qxpg_sldata3_t)
#define back_qxpg_sldata4_t  SL_FMM_TYPE(back_qxpg_sldata4_t)

#define back_qxpg_elements_t  SL_FMM_FUNC(back_qxpg_elements_t)
#define back_qxpg_merge2x_f   SL_FMM_FUNC(back_qxpg_merge2x_f)

#define back_qxpg_elem_set_size        SL_FMM_FUNC(back_qxpg_elem_set_size)
#define back_qxpg_elem_get_size        SL_FMM_FUNC(back_qxpg_elem_get_size)
#define back_qxpg_elem_set_max_size    SL_FMM_FUNC(back_qxpg_elem_set_max_size)
#define back_qxpg_elem_set_keys        SL_FMM_FUNC(back_qxpg_elem_set_keys)
#define back_qxpg_elem_set_data        SL_FMM_FUNC(back_qxpg_elem_set_data)
#define back_qxpg_elem_set_block       SL_FMM_FUNC(back_qxpg_elem_set_block)
#define back_qxpg_elem_set_block_size  SL_FMM_FUNC(back_qxpg_elem_set_block_size)

#define back_qxpg_elements_alloc_from_blocks  SL_FMM_FUNC(back_qxpg_elements_alloc_from_blocks)
#define back_qxpg_sort_radix                  SL_FMM_FUNC(back_qxpg_sort_radix)
#define back_qxpg_merge2_compo_hula           SL_FMM_FUNC(back_qxpg_merge2_compo_hula)
#define back_qxpg_merge2_memory_adaptive      SL_FMM_FUNC(back_qxpg_merge2_memory_adaptive)
#define back_qxpg_sn_batcher                  SL_FMM_FUNC(back_qxpg_sn_batcher)
#define back_qxpg_mpi_mergek                  SL_FMM_FUNC(back_qxpg_mpi_mergek)
#define back_qxpg_mpi_datatypes_init          SL_FMM_FUNC(back_qxpg_mpi_datatypes_init)
#define back_qxpg_mpi_datatypes_release       SL_FMM_FUNC(back_qxpg_mpi_datatypes_release)
#define back_qxpg_elements_alloc              SL_FMM_FUNC(back_qxpg_elements_alloc)
#define back_qxpg_elements_free               SL_FMM_FUNC(back_qxpg_elements_free)
#define back_qxpg_mpi_sort_back               SL_FMM_FUNC(back_qxpg_mpi_sort_back)


/* back_qx_g */
#define back_qx_g_slint_t      SL_FMM_TYPE(back_qx_g_slint_t)
#define back_qx_g_slint_fmt    SL_FMM_MACRO(back_qx_g_slint_fmt)

#define back_qx_g_slkey_t    SL_FMM_TYPE(back_qx_g_slkey_t)
#define back_qx_g_sldata0_t  SL_FMM_TYPE(back_qx_g_sldata0_t)
#define back_qx_g_sldata1_t  SL_FMM_TYPE(back_qx_g_sldata1_t)
#define back_qx_g_sldata2_t  SL_FMM_TYPE(back_qx_g_sldata2_t)
#define back_qx_g_sldata3_t  SL_FMM_TYPE(back_qx_g_sldata3_t)
#define back_qx_g_sldata4_t  SL_FMM_TYPE(back_qx_g_sldata4_t)

#define back_qx_g_elements_t  SL_FMM_FUNC(back_qx_g_elements_t)
#define back_qx_g_merge2x_f   SL_FMM_FUNC(back_qx_g_merge2x_f)

#define back_qx_g_elem_set_size        SL_FMM_FUNC(back_qx_g_elem_set_size)
#define back_qx_g_elem_get_size        SL_FMM_FUNC(back_qx_g_elem_get_size)
#define back_qx_g_elem_set_max_size    SL_FMM_FUNC(back_qx_g_elem_set_max_size)
#define back_qx_g_elem_set_keys        SL_FMM_FUNC(back_qx_g_elem_set_keys)
#define back_qx_g_elem_set_data        SL_FMM_FUNC(back_qx_g_elem_set_data)
#define back_qx_g_elem_set_block       SL_FMM_FUNC(back_qx_g_elem_set_block)
#define back_qx_g_elem_set_block_size  SL_FMM_FUNC(back_qx_g_elem_set_block_size)

#define back_qx_g_elements_alloc_from_blocks  SL_FMM_FUNC(back_qx_g_elements_alloc_from_blocks)
#define back_qx_g_sort_radix                  SL_FMM_FUNC(back_qx_g_sort_radix)
#define back_qx_g_merge2_compo_hula           SL_FMM_FUNC(back_qx_g_merge2_compo_hula)
#define back_qx_g_merge2_memory_adaptive      SL_FMM_FUNC(back_qx_g_merge2_memory_adaptive)
#define back_qx_g_sn_batcher                  SL_FMM_FUNC(back_qx_g_sn_batcher)
#define back_qx_g_mpi_mergek                  SL_FMM_FUNC(back_qx_g_mpi_mergek)
#define back_qx_g_mpi_datatypes_init          SL_FMM_FUNC(back_qx_g_mpi_datatypes_init)
#define back_qx_g_mpi_datatypes_release       SL_FMM_FUNC(back_qx_g_mpi_datatypes_release)
#define back_qx_g_elements_alloc              SL_FMM_FUNC(back_qx_g_elements_alloc)
#define back_qx_g_elements_free               SL_FMM_FUNC(back_qx_g_elements_free)
#define back_qx_g_mpi_sort_back               SL_FMM_FUNC(back_qx_g_mpi_sort_back)


/* back_q_pg */
#define back_q_pg_slint_t      SL_FMM_TYPE(back_q_pg_slint_t)
#define back_q_pg_slint_fmt    SL_FMM_MACRO(back_q_pg_slint_fmt)

#define back_q_pg_slkey_t    SL_FMM_TYPE(back_q_pg_slkey_t)
#define back_q_pg_sldata0_t  SL_FMM_TYPE(back_q_pg_sldata0_t)
#define back_q_pg_sldata1_t  SL_FMM_TYPE(back_q_pg_sldata1_t)
#define back_q_pg_sldata2_t  SL_FMM_TYPE(back_q_pg_sldata2_t)
#define back_q_pg_sldata3_t  SL_FMM_TYPE(back_q_pg_sldata3_t)
#define back_q_pg_sldata4_t  SL_FMM_TYPE(back_q_pg_sldata4_t)

#define back_q_pg_elements_t  SL_FMM_FUNC(back_q_pg_elements_t)
#define back_q_pg_merge2x_f   SL_FMM_FUNC(back_q_pg_merge2x_f)

#define back_q_pg_elem_set_size        SL_FMM_FUNC(back_q_pg_elem_set_size)
#define back_q_pg_elem_get_size        SL_FMM_FUNC(back_q_pg_elem_get_size)
#define back_q_pg_elem_set_max_size    SL_FMM_FUNC(back_q_pg_elem_set_max_size)
#define back_q_pg_elem_set_keys        SL_FMM_FUNC(back_q_pg_elem_set_keys)
#define back_q_pg_elem_set_data        SL_FMM_FUNC(back_q_pg_elem_set_data)
#define back_q_pg_elem_set_block       SL_FMM_FUNC(back_q_pg_elem_set_block)
#define back_q_pg_elem_set_block_size  SL_FMM_FUNC(back_q_pg_elem_set_block_size)

#define back_q_pg_elements_alloc_from_blocks  SL_FMM_FUNC(back_q_pg_elements_alloc_from_blocks)
#define back_q_pg_sort_radix                  SL_FMM_FUNC(back_q_pg_sort_radix)
#define back_q_pg_merge2_compo_hula           SL_FMM_FUNC(back_q_pg_merge2_compo_hula)
#define back_q_pg_merge2_memory_adaptive      SL_FMM_FUNC(back_q_pg_merge2_memory_adaptive)
#define back_q_pg_sn_batcher                  SL_FMM_FUNC(back_q_pg_sn_batcher)
#define back_q_pg_mpi_mergek                  SL_FMM_FUNC(back_q_pg_mpi_mergek)
#define back_q_pg_mpi_datatypes_init          SL_FMM_FUNC(back_q_pg_mpi_datatypes_init)
#define back_q_pg_mpi_datatypes_release       SL_FMM_FUNC(back_q_pg_mpi_datatypes_release)
#define back_q_pg_elements_alloc              SL_FMM_FUNC(back_q_pg_elements_alloc)
#define back_q_pg_elements_free               SL_FMM_FUNC(back_q_pg_elements_free)
#define back_q_pg_mpi_sort_back               SL_FMM_FUNC(back_q_pg_mpi_sort_back)


/* back_q__g */
#define back_q__g_slint_t      SL_FMM_TYPE(back_q__g_slint_t)
#define back_q__g_slint_fmt    SL_FMM_MACRO(back_q__g_slint_fmt)

#define back_q__g_slkey_t    SL_FMM_TYPE(back_q__g_slkey_t)
#define back_q__g_sldata0_t  SL_FMM_TYPE(back_q__g_sldata0_t)
#define back_q__g_sldata1_t  SL_FMM_TYPE(back_q__g_sldata1_t)
#define back_q__g_sldata2_t  SL_FMM_TYPE(back_q__g_sldata2_t)
#define back_q__g_sldata3_t  SL_FMM_TYPE(back_q__g_sldata3_t)
#define back_q__g_sldata4_t  SL_FMM_TYPE(back_q__g_sldata4_t)

#define back_q__g_elements_t  SL_FMM_FUNC(back_q__g_elements_t)
#define back_q__g_merge2x_f   SL_FMM_FUNC(back_q__g_merge2x_f)

#define back_q__g_elem_set_size        SL_FMM_FUNC(back_q__g_elem_set_size)
#define back_q__g_elem_get_size        SL_FMM_FUNC(back_q__g_elem_get_size)
#define back_q__g_elem_set_max_size    SL_FMM_FUNC(back_q__g_elem_set_max_size)
#define back_q__g_elem_set_keys        SL_FMM_FUNC(back_q__g_elem_set_keys)
#define back_q__g_elem_set_data        SL_FMM_FUNC(back_q__g_elem_set_data)
#define back_q__g_elem_set_block       SL_FMM_FUNC(back_q__g_elem_set_block)
#define back_q__g_elem_set_block_size  SL_FMM_FUNC(back_q__g_elem_set_block_size)

#define back_q__g_elements_alloc_from_blocks  SL_FMM_FUNC(back_q__g_elements_alloc_from_blocks)
#define back_q__g_sort_radix                  SL_FMM_FUNC(back_q__g_sort_radix)
#define back_q__g_merge2_compo_hula           SL_FMM_FUNC(back_q__g_merge2_compo_hula)
#define back_q__g_merge2_memory_adaptive      SL_FMM_FUNC(back_q__g_merge2_memory_adaptive)
#define back_q__g_sn_batcher                  SL_FMM_FUNC(back_q__g_sn_batcher)
#define back_q__g_mpi_mergek                  SL_FMM_FUNC(back_q__g_mpi_mergek)
#define back_q__g_mpi_datatypes_init          SL_FMM_FUNC(back_q__g_mpi_datatypes_init)
#define back_q__g_mpi_datatypes_release       SL_FMM_FUNC(back_q__g_mpi_datatypes_release)
#define back_q__g_elements_alloc              SL_FMM_FUNC(back_q__g_elements_alloc)
#define back_q__g_elements_free               SL_FMM_FUNC(back_q__g_elements_free)
#define back_q__g_mpi_sort_back               SL_FMM_FUNC(back_q__g_mpi_sort_back)


/* back_qxpgl */
#define back_qxpgl_slint_t      SL_FMM_TYPE(back_qxpgl_slint_t)
#define back_qxpgl_slint_fmt    SL_FMM_MACRO(back_qxpgl_slint_fmt)

#define back_qxpgl_slkey_t    SL_FMM_TYPE(back_qxpgl_slkey_t)
#define back_qxpgl_sldata0_t  SL_FMM_TYPE(back_qxpgl_sldata0_t)
#define back_qxpgl_sldata1_t  SL_FMM_TYPE(back_qxpgl_sldata1_t)
#define back_qxpgl_sldata2_t  SL_FMM_TYPE(back_qxpgl_sldata2_t)
#define back_qxpgl_sldata3_t  SL_FMM_TYPE(back_qxpgl_sldata3_t)
#define back_qxpgl_sldata4_t  SL_FMM_TYPE(back_qxpgl_sldata4_t)

#define back_qxpgl_elements_t  SL_FMM_FUNC(back_qxpgl_elements_t)
#define back_qxpgl_merge2x_f   SL_FMM_FUNC(back_qxpgl_merge2x_f)

#define back_qxpgl_elem_set_size        SL_FMM_FUNC(back_qxpgl_elem_set_size)
#define back_qxpgl_elem_get_size        SL_FMM_FUNC(back_qxpgl_elem_get_size)
#define back_qxpgl_elem_set_max_size    SL_FMM_FUNC(back_qxpgl_elem_set_max_size)
#define back_qxpgl_elem_set_keys        SL_FMM_FUNC(back_qxpgl_elem_set_keys)
#define back_qxpgl_elem_set_data        SL_FMM_FUNC(back_qxpgl_elem_set_data)
#define back_qxpgl_elem_set_block       SL_FMM_FUNC(back_qxpgl_elem_set_block)
#define back_qxpgl_elem_set_block_size  SL_FMM_FUNC(back_qxpgl_elem_set_block_size)

#define back_qxpgl_elements_alloc_from_blocks  SL_FMM_FUNC(back_qxpgl_elements_alloc_from_blocks)
#define back_qxpgl_sort_radix                  SL_FMM_FUNC(back_qxpgl_sort_radix)
#define back_qxpgl_merge2_compo_hula           SL_FMM_FUNC(back_qxpgl_merge2_compo_hula)
#define back_qxpgl_merge2_memory_adaptive      SL_FMM_FUNC(back_qxpgl_merge2_memory_adaptive)
#define back_qxpgl_sn_batcher                  SL_FMM_FUNC(back_qxpgl_sn_batcher)
#define back_qxpgl_mpi_mergek                  SL_FMM_FUNC(back_qxpgl_mpi_mergek)
#define back_qxpgl_mpi_datatypes_init          SL_FMM_FUNC(back_qxpgl_mpi_datatypes_init)
#define back_qxpgl_mpi_datatypes_release       SL_FMM_FUNC(back_qxpgl_mpi_datatypes_release)
#define back_qxpgl_elements_alloc              SL_FMM_FUNC(back_qxpgl_elements_alloc)
#define back_qxpgl_elements_free               SL_FMM_FUNC(back_qxpgl_elements_free)
#define back_qxpgl_mpi_sort_back               SL_FMM_FUNC(back_qxpgl_mpi_sort_back)


/* back_qx_gl */
#define back_qx_gl_slint_t      SL_FMM_TYPE(back_qx_gl_slint_t)
#define back_qx_gl_slint_fmt    SL_FMM_MACRO(back_qx_gl_slint_fmt)

#define back_qx_gl_slkey_t    SL_FMM_TYPE(back_qx_gl_slkey_t)
#define back_qx_gl_sldata0_t  SL_FMM_TYPE(back_qx_gl_sldata0_t)
#define back_qx_gl_sldata1_t  SL_FMM_TYPE(back_qx_gl_sldata1_t)
#define back_qx_gl_sldata2_t  SL_FMM_TYPE(back_qx_gl_sldata2_t)
#define back_qx_gl_sldata3_t  SL_FMM_TYPE(back_qx_gl_sldata3_t)
#define back_qx_gl_sldata4_t  SL_FMM_TYPE(back_qx_gl_sldata4_t)

#define back_qx_gl_elements_t  SL_FMM_FUNC(back_qx_gl_elements_t)
#define back_qx_gl_merge2x_f   SL_FMM_FUNC(back_qx_gl_merge2x_f)

#define back_qx_gl_elem_set_size        SL_FMM_FUNC(back_qx_gl_elem_set_size)
#define back_qx_gl_elem_get_size        SL_FMM_FUNC(back_qx_gl_elem_get_size)
#define back_qx_gl_elem_set_max_size    SL_FMM_FUNC(back_qx_gl_elem_set_max_size)
#define back_qx_gl_elem_set_keys        SL_FMM_FUNC(back_qx_gl_elem_set_keys)
#define back_qx_gl_elem_set_data        SL_FMM_FUNC(back_qx_gl_elem_set_data)
#define back_qx_gl_elem_set_block       SL_FMM_FUNC(back_qx_gl_elem_set_block)
#define back_qx_gl_elem_set_block_size  SL_FMM_FUNC(back_qx_gl_elem_set_block_size)

#define back_qx_gl_elements_alloc_from_blocks  SL_FMM_FUNC(back_qx_gl_elements_alloc_from_blocks)
#define back_qx_gl_sort_radix                  SL_FMM_FUNC(back_qx_gl_sort_radix)
#define back_qx_gl_merge2_compo_hula           SL_FMM_FUNC(back_qx_gl_merge2_compo_hula)
#define back_qx_gl_merge2_memory_adaptive      SL_FMM_FUNC(back_qx_gl_merge2_memory_adaptive)
#define back_qx_gl_sn_batcher                  SL_FMM_FUNC(back_qx_gl_sn_batcher)
#define back_qx_gl_mpi_mergek                  SL_FMM_FUNC(back_qx_gl_mpi_mergek)
#define back_qx_gl_mpi_datatypes_init          SL_FMM_FUNC(back_qx_gl_mpi_datatypes_init)
#define back_qx_gl_mpi_datatypes_release       SL_FMM_FUNC(back_qx_gl_mpi_datatypes_release)
#define back_qx_gl_elements_alloc              SL_FMM_FUNC(back_qx_gl_elements_alloc)
#define back_qx_gl_elements_free               SL_FMM_FUNC(back_qx_gl_elements_free)
#define back_qx_gl_mpi_sort_back               SL_FMM_FUNC(back_qx_gl_mpi_sort_back)


/* back_q_pgl */
#define back_q_pgl_slint_t      SL_FMM_TYPE(back_q_pgl_slint_t)
#define back_q_pgl_slint_fmt    SL_FMM_MACRO(back_q_pgl_slint_fmt)

#define back_q_pgl_slkey_t    SL_FMM_TYPE(back_q_pgl_slkey_t)
#define back_q_pgl_sldata0_t  SL_FMM_TYPE(back_q_pgl_sldata0_t)
#define back_q_pgl_sldata1_t  SL_FMM_TYPE(back_q_pgl_sldata1_t)
#define back_q_pgl_sldata2_t  SL_FMM_TYPE(back_q_pgl_sldata2_t)
#define back_q_pgl_sldata3_t  SL_FMM_TYPE(back_q_pgl_sldata3_t)
#define back_q_pgl_sldata4_t  SL_FMM_TYPE(back_q_pgl_sldata4_t)

#define back_q_pgl_elements_t  SL_FMM_FUNC(back_q_pgl_elements_t)
#define back_q_pgl_merge2x_f   SL_FMM_FUNC(back_q_pgl_merge2x_f)

#define back_q_pgl_elem_set_size        SL_FMM_FUNC(back_q_pgl_elem_set_size)
#define back_q_pgl_elem_get_size        SL_FMM_FUNC(back_q_pgl_elem_get_size)
#define back_q_pgl_elem_set_max_size    SL_FMM_FUNC(back_q_pgl_elem_set_max_size)
#define back_q_pgl_elem_set_keys        SL_FMM_FUNC(back_q_pgl_elem_set_keys)
#define back_q_pgl_elem_set_data        SL_FMM_FUNC(back_q_pgl_elem_set_data)
#define back_q_pgl_elem_set_block       SL_FMM_FUNC(back_q_pgl_elem_set_block)
#define back_q_pgl_elem_set_block_size  SL_FMM_FUNC(back_q_pgl_elem_set_block_size)

#define back_q_pgl_elements_alloc_from_blocks  SL_FMM_FUNC(back_q_pgl_elements_alloc_from_blocks)
#define back_q_pgl_sort_radix                  SL_FMM_FUNC(back_q_pgl_sort_radix)
#define back_q_pgl_merge2_compo_hula           SL_FMM_FUNC(back_q_pgl_merge2_compo_hula)
#define back_q_pgl_merge2_memory_adaptive      SL_FMM_FUNC(back_q_pgl_merge2_memory_adaptive)
#define back_q_pgl_sn_batcher                  SL_FMM_FUNC(back_q_pgl_sn_batcher)
#define back_q_pgl_mpi_mergek                  SL_FMM_FUNC(back_q_pgl_mpi_mergek)
#define back_q_pgl_mpi_datatypes_init          SL_FMM_FUNC(back_q_pgl_mpi_datatypes_init)
#define back_q_pgl_mpi_datatypes_release       SL_FMM_FUNC(back_q_pgl_mpi_datatypes_release)
#define back_q_pgl_elements_alloc              SL_FMM_FUNC(back_q_pgl_elements_alloc)
#define back_q_pgl_elements_free               SL_FMM_FUNC(back_q_pgl_elements_free)
#define back_q_pgl_mpi_sort_back               SL_FMM_FUNC(back_q_pgl_mpi_sort_back)


/* back_q__gl */
#define back_q__gl_slint_t      SL_FMM_TYPE(back_q__gl_slint_t)
#define back_q__gl_slint_fmt    SL_FMM_MACRO(back_q__gl_slint_fmt)

#define back_q__gl_slkey_t    SL_FMM_TYPE(back_q__gl_slkey_t)
#define back_q__gl_sldata0_t  SL_FMM_TYPE(back_q__gl_sldata0_t)
#define back_q__gl_sldata1_t  SL_FMM_TYPE(back_q__gl_sldata1_t)
#define back_q__gl_sldata2_t  SL_FMM_TYPE(back_q__gl_sldata2_t)
#define back_q__gl_sldata3_t  SL_FMM_TYPE(back_q__gl_sldata3_t)
#define back_q__gl_sldata4_t  SL_FMM_TYPE(back_q__gl_sldata4_t)

#define back_q__gl_elements_t  SL_FMM_FUNC(back_q__gl_elements_t)
#define back_q__gl_merge2x_f   SL_FMM_FUNC(back_q__gl_merge2x_f)

#define back_q__gl_elem_set_size        SL_FMM_FUNC(back_q__gl_elem_set_size)
#define back_q__gl_elem_get_size        SL_FMM_FUNC(back_q__gl_elem_get_size)
#define back_q__gl_elem_set_max_size    SL_FMM_FUNC(back_q__gl_elem_set_max_size)
#define back_q__gl_elem_set_keys        SL_FMM_FUNC(back_q__gl_elem_set_keys)
#define back_q__gl_elem_set_data        SL_FMM_FUNC(back_q__gl_elem_set_data)
#define back_q__gl_elem_set_block       SL_FMM_FUNC(back_q__gl_elem_set_block)
#define back_q__gl_elem_set_block_size  SL_FMM_FUNC(back_q__gl_elem_set_block_size)

#define back_q__gl_elements_alloc_from_blocks  SL_FMM_FUNC(back_q__gl_elements_alloc_from_blocks)
#define back_q__gl_sort_radix                  SL_FMM_FUNC(back_q__gl_sort_radix)
#define back_q__gl_merge2_compo_hula           SL_FMM_FUNC(back_q__gl_merge2_compo_hula)
#define back_q__gl_merge2_memory_adaptive      SL_FMM_FUNC(back_q__gl_merge2_memory_adaptive)
#define back_q__gl_sn_batcher                  SL_FMM_FUNC(back_q__gl_sn_batcher)
#define back_q__gl_mpi_mergek                  SL_FMM_FUNC(back_q__gl_mpi_mergek)
#define back_q__gl_mpi_datatypes_init          SL_FMM_FUNC(back_q__gl_mpi_datatypes_init)
#define back_q__gl_mpi_datatypes_release       SL_FMM_FUNC(back_q__gl_mpi_datatypes_release)
#define back_q__gl_elements_alloc              SL_FMM_FUNC(back_q__gl_elements_alloc)
#define back_q__gl_elements_free               SL_FMM_FUNC(back_q__gl_elements_free)
#define back_q__gl_mpi_sort_back               SL_FMM_FUNC(back_q__gl_mpi_sort_back)


/* back_idx */
#define back_idx_slint_t      SL_FMM_TYPE(back_idx_slint_t)
#define back_idx_slint_fmt    SL_FMM_MACRO(back_idx_slint_fmt)

#define back_idx_slkey_t    SL_FMM_TYPE(back_idx_slkey_t)
#define back_idx_sldata0_t  SL_FMM_TYPE(back_idx_sldata0_t)
#define back_idx_sldata1_t  SL_FMM_TYPE(back_idx_sldata1_t)
#define back_idx_sldata2_t  SL_FMM_TYPE(back_idx_sldata2_t)
#define back_idx_sldata3_t  SL_FMM_TYPE(back_idx_sldata3_t)
#define back_idx_sldata4_t  SL_FMM_TYPE(back_idx_sldata4_t)

#define back_idx_elements_t  SL_FMM_FUNC(back_idx_elements_t)
#define back_idx_merge2x_f   SL_FMM_FUNC(back_idx_merge2x_f)

#define back_idx_elem_set_size        SL_FMM_FUNC(back_idx_elem_set_size)
#define back_idx_elem_get_size        SL_FMM_FUNC(back_idx_elem_get_size)
#define back_idx_elem_set_max_size    SL_FMM_FUNC(back_idx_elem_set_max_size)
#define back_idx_elem_set_keys        SL_FMM_FUNC(back_idx_elem_set_keys)
#define back_idx_elem_set_data        SL_FMM_FUNC(back_idx_elem_set_data)
#define back_idx_elem_set_block       SL_FMM_FUNC(back_idx_elem_set_block)
#define back_idx_elem_set_block_size  SL_FMM_FUNC(back_idx_elem_set_block_size)

#define back_idx_elements_alloc_from_blocks  SL_FMM_FUNC(back_idx_elements_alloc_from_blocks)
#define back_idx_sort_radix                  SL_FMM_FUNC(back_idx_sort_radix)
#define back_idx_merge2_compo_hula           SL_FMM_FUNC(back_idx_merge2_compo_hula)
#define back_idx_merge2_memory_adaptive      SL_FMM_FUNC(back_idx_merge2_memory_adaptive)
#define back_idx_sn_batcher                  SL_FMM_FUNC(back_idx_sn_batcher)
#define back_idx_mpi_mergek                  SL_FMM_FUNC(back_idx_mpi_mergek)
#define back_idx_mpi_datatypes_init          SL_FMM_FUNC(back_idx_mpi_datatypes_init)
#define back_idx_mpi_datatypes_release       SL_FMM_FUNC(back_idx_mpi_datatypes_release)
#define back_idx_elements_alloc              SL_FMM_FUNC(back_idx_elements_alloc)
#define back_idx_elements_free               SL_FMM_FUNC(back_idx_elements_free)
#define back_idx_mpi_sort_back               SL_FMM_FUNC(back_idx_mpi_sort_back)


/* sl_fmm */
/*#define fmm_sort_front                     SL_FMM_FUNC(fmm_sort_front)
#define fmm_sort_front_                    SL_FMM_FUNC(fmm_sort_front_)
#define fmm_sort_front_mem                 SL_FMM_FUNC(fmm_sort_front_mem)
#define fmm_sort_front_mem_                SL_FMM_FUNC(fmm_sort_front_mem_)
#define fmm_sort_front_3bit                SL_FMM_FUNC(fmm_sort_front_3bit)
#define fmm_sort_front_3bit_               SL_FMM_FUNC(fmm_sort_front_3bit_)
#define fmm_sort_front_3bit_mem            SL_FMM_FUNC(fmm_sort_front_3bit_mem)
#define fmm_sort_front_3bit_mem_           SL_FMM_FUNC(fmm_sort_front_3bit_mem_)
#define fmm_sort_back                      SL_FMM_FUNC(fmm_sort_back)
#define fmm_sort_back_                     SL_FMM_FUNC(fmm_sort_back_)
#define fmm_sort_back_mem                  SL_FMM_FUNC(fmm_sort_back_mem)
#define fmm_sort_back_mem_                 SL_FMM_FUNC(fmm_sort_back_mem_)
#define mpi_fmm_sort_front                 SL_FMM_FUNC(mpi_fmm_sort_front)
#define mpi_fmm_sort_front_                SL_FMM_FUNC(mpi_fmm_sort_front_)
#define mpi_fmm_sort_front_mem             SL_FMM_FUNC(mpi_fmm_sort_front_mem)
#define mpi_fmm_sort_front_mem_            SL_FMM_FUNC(mpi_fmm_sort_front_mem_)
#define mpi_fmm_sort_front_3bit            SL_FMM_FUNC(mpi_fmm_sort_front_3bit)
#define mpi_fmm_sort_front_3bit_           SL_FMM_FUNC(mpi_fmm_sort_front_3bit_)
#define mpi_fmm_sort_front_3bit_mem        SL_FMM_FUNC(mpi_fmm_sort_front_3bit_mem)
#define mpi_fmm_sort_front_3bit_mem_       SL_FMM_FUNC(mpi_fmm_sort_front_3bit_mem_)
#define mpi_fmm_sort_front_rebalance       SL_FMM_FUNC(mpi_fmm_sort_front_rebalance)
#define mpi_fmm_sort_front_rebalance_      SL_FMM_FUNC(mpi_fmm_sort_front_rebalance_)
#define mpi_fmm_sort_front_rebalance_mem   SL_FMM_FUNC(mpi_fmm_sort_front_rebalance_mem)
#define mpi_fmm_sort_front_rebalance_mem_  SL_FMM_FUNC(mpi_fmm_sort_front_rebalance_mem_)
#define mpi_fmm_sort_back                  SL_FMM_FUNC(mpi_fmm_sort_back)
#define mpi_fmm_sort_back_                 SL_FMM_FUNC(mpi_fmm_sort_back_)
#define mpi_fmm_sort_back_mem              SL_FMM_FUNC(mpi_fmm_sort_back_mem)
#define mpi_fmm_sort_back_mem_             SL_FMM_FUNC(mpi_fmm_sort_back_mem_)
#define mpi_fmm_resort_init                SL_FMM_FUNC(mpi_fmm_resort_init)
#define mpi_fmm_resort_init_               SL_FMM_FUNC(mpi_fmm_resort_init_)*/
#define fmm_resort_create                  SL_FMM_FUNC(fmm_resort_create)
#define fmm_resort_destroy                 SL_FMM_FUNC(fmm_resort_destroy)

#define mpi_fmm_sort_front_part             SL_FMM_VAR(mpi_fmm_sort_front_part)
#define mpi_fmm_sort_back_part              SL_FMM_VAR(mpi_fmm_sort_back_part)
#define mpi_fmm_sort_front_merge_presorted  SL_FMM_VAR(mpi_fmm_sort_front_merge_presorted)


#endif /* __RENAME_FMM_SORT_H__ */
