
#ifndef __RENAME_PEPC_SORT_H__
#define __RENAME_PEPC_SORT_H__


#define SL_PEPC_CONCAT(_a_, _b_)   SL_PEPC_CONCAT_(_a_, _b_)
#define SL_PEPC_CONCAT_(_a_, _b_)  _a_##_b_

#ifdef SL_PEPC_PREFIX
# define SL_PEPC_FUNC(_f_)   SL_PEPC_CONCAT(SL_PEPC_PREFIX, _f_)
# define SL_PEPC_VAR(_v_)    SL_PEPC_CONCAT(SL_PEPC_PREFIX, _v_)
# define SL_PEPC_TYPE(_t_)   SL_PEPC_CONCAT(SL_PEPC_PREFIX, _t_)
# define SL_PEPC_MACRO(_m_)  SL_PEPC_CONCAT(SL_PEPC_PREFIX, _m_)
#else
# define SL_PEPC_FUNC(_f_)   _f_
# define SL_PEPC_VAR(_v_)    _v_
# define SL_PEPC_TYPE(_t_)   _t_
# define SL_PEPC_MACRO(_m_)  _m_
#endif


/* pepckeys */
#define pepckeys_slint_t    SL_PEPC_TYPE(pepckeys_slint_t)
#define pepckeys_slint_fmt  SL_PEPC_MACRO(pepckeys_slint_fmt)

#define pepckeys_slkey_t    SL_PEPC_TYPE(pepckeys_slkey_t)
#define pepckeys_slindex_t  SL_PEPC_TYPE(pepckeys_slindex_t)
#define pepckeys_sldata0_t  SL_PEPC_TYPE(pepckeys_sldata0_t)

#define pepckeys_elements_t  SL_PEPC_TYPE(pepckeys_elements_t)
#define pepckeys_partcond_t  SL_PEPC_TYPE(pepckeys_partcond_t)

#define pepckeys_elem_set_size      SL_PEPC_FUNC(pepckeys_elem_set_size)
#define pepckeys_elem_set_max_size  SL_PEPC_FUNC(pepckeys_elem_set_max_size)
#define pepckeys_elem_set_keys      SL_PEPC_FUNC(pepckeys_elem_set_keys)
#define pepckeys_elem_set_indices   SL_PEPC_FUNC(pepckeys_elem_set_indices)
#define pepckeys_elem_set_data      SL_PEPC_FUNC(pepckeys_elem_set_data)

#define pepckeys_mpi_get_grid_properties            SL_PEPC_FUNC(pepckeys_mpi_get_grid_properties)
#define pepckeys_mpi_datatypes_init                 SL_PEPC_FUNC(pepckeys_mpi_datatypes_init)
#define pepckeys_sort_radix                         SL_PEPC_FUNC(pepckeys_sort_radix)
#define pepckeys_elements_validate_order            SL_PEPC_FUNC(pepckeys_elements_validate_order)
#define pepckeys_mpi_partition_exact_radix          SL_PEPC_FUNC(pepckeys_mpi_partition_exact_radix)
#define pepckeys_mpi_subgroups_create               SL_PEPC_FUNC(pepckeys_mpi_subgroups_create)
#define pepckeys_mpi_subgroups_delete               SL_PEPC_FUNC(pepckeys_mpi_subgroups_delete)
#define pepckeys_mpi_partition_exact_radix_2groups  SL_PEPC_FUNC(pepckeys_mpi_partition_exact_radix_2groups)
#define pepckeys_mpi_partition_exact_radix_ngroups  SL_PEPC_FUNC(pepckeys_mpi_partition_exact_radix_ngroups)
#define pepckeys_mpi_partition_sample_regular       SL_PEPC_FUNC(pepckeys_mpi_partition_sample_regular)
#define pepckeys_counts2displs                      SL_PEPC_FUNC(pepckeys_counts2displs)
#define pepckeys_mergep_heap_idx                    SL_PEPC_FUNC(pepckeys_mergep_heap_idx)
#define pepckeys_mpi_datatypes_release              SL_PEPC_FUNC(pepckeys_mpi_datatypes_release)
#define pepckeys_mpi_elements_validate_order        SL_PEPC_FUNC(pepckeys_mpi_elements_validate_order)
#define pepckeys_elements_ncopy                     SL_PEPC_FUNC(pepckeys_elements_ncopy)


/* pepcparts */
#define pepcparts_slint_t    SL_PEPC_TYPE(pepcparts_slint_t)
#define pepcparts_slint_fmt  SL_PEPC_MACRO(pepcparts_slint_fmt)

#define pepcparts_slkey_t          SL_PEPC_TYPE(pepcparts_slkey_t)
#define pepcparts_sl_key_type_fmt  SL_PEPC_MACRO(pepcparts_sl_key_type_fmt)

#define pepcparts_sldata0_t   SL_PEPC_TYPE(pepcparts_sldata0_t)
#define pepcparts_sldata1_t   SL_PEPC_TYPE(pepcparts_sldata1_t)
#define pepcparts_sldata2_t   SL_PEPC_TYPE(pepcparts_sldata2_t)
#define pepcparts_sldata3_t   SL_PEPC_TYPE(pepcparts_sldata3_t)
#define pepcparts_sldata4_t   SL_PEPC_TYPE(pepcparts_sldata4_t)
#define pepcparts_sldata5_t   SL_PEPC_TYPE(pepcparts_sldata5_t)
#define pepcparts_sldata6_t   SL_PEPC_TYPE(pepcparts_sldata6_t)
#define pepcparts_sldata7_t   SL_PEPC_TYPE(pepcparts_sldata7_t)
#define pepcparts_sldata8_t   SL_PEPC_TYPE(pepcparts_sldata8_t)
#define pepcparts_sldata9_t   SL_PEPC_TYPE(pepcparts_sldata9_t)
#define pepcparts_sldata10_t  SL_PEPC_TYPE(pepcparts_sldata10_t)
#define pepcparts_sldata11_t  SL_PEPC_TYPE(pepcparts_sldata11_t)
#define pepcparts_sldata12_t  SL_PEPC_TYPE(pepcparts_sldata12_t)
#define pepcparts_slindex_t   SL_PEPC_TYPE(pepcparts_slindex_t)

#define pepcparts_elements_t         SL_PEPC_TYPE(pepcparts_elements_t)
#define pepcparts_packed_elements_t  SL_PEPC_TYPE(pepcparts_packed_elements_t)
#define pepcparts_partcond_t         SL_PEPC_TYPE(pepcparts_partcond_t)

#define pepcparts_elem_set_size      SL_PEPC_FUNC(pepcparts_elem_set_size)
#define pepcparts_elem_set_max_size  SL_PEPC_FUNC(pepcparts_elem_set_max_size)
#define pepcparts_elem_set_keys      SL_PEPC_FUNC(pepcparts_elem_set_keys)
#define pepcparts_elem_set_data      SL_PEPC_FUNC(pepcparts_elem_set_data)
#define pepcparts_elem_set_indices   SL_PEPC_FUNC(pepcparts_elem_set_indices)

#define pepcparts_pelem_set_size      SL_PEPC_FUNC(pepcparts_pelem_set_size)
#define pepcparts_pelem_set_max_size  SL_PEPC_FUNC(pepcparts_pelem_set_max_size)
#define pepcparts_pelem_set_elements  SL_PEPC_FUNC(pepcparts_pelem_set_elements)

#define pepcparts_mpi_get_grid_properties               SL_PEPC_FUNC(pepcparts_mpi_get_grid_properties)
#define pepcparts_mpi_datatypes_init                    SL_PEPC_FUNC(pepcparts_mpi_datatypes_init)
#define pepcparts_elements_pack_indexed                 SL_PEPC_FUNC(pepcparts_elements_pack_indexed)
#define pepcparts_mpi_elements_packed_datatype_create   SL_PEPC_FUNC(pepcparts_mpi_elements_packed_datatype_create)
#define pepcparts_mpi_elements_packed_datatype_destroy  SL_PEPC_FUNC(pepcparts_mpi_elements_packed_datatype_destroy)
#define pepcparts_mergep_heap_unpack_idxonly            SL_PEPC_FUNC(pepcparts_mergep_heap_unpack_idxonly)
#define pepcparts_elements_unpack_indexed               SL_PEPC_FUNC(pepcparts_elements_unpack_indexed)
#define pepcparts_elements_unpack_keys                  SL_PEPC_FUNC(pepcparts_elements_unpack_keys)
#define pepcparts_mergep_heap_unpack_idx                SL_PEPC_FUNC(pepcparts_mergep_heap_unpack_idx)
#define pepcparts_mpi_elements_validate_order           SL_PEPC_FUNC(pepcparts_mpi_elements_validate_order)
#define pepcparts_mpi_datatypes_release                 SL_PEPC_FUNC(pepcparts_mpi_datatypes_release)


/* sl_pepc */
#define sl_pepc_check_fortran2c_types   SL_PEPC_FUNC(sl_pepc_check_fortran2c_types)
#define sl_pepc_check_fortran2c_types_  SL_PEPC_FUNC(sl_pepc_check_fortran2c_types_)
#define sl_pepc_border_stats            SL_PEPC_FUNC(sl_pepc_border_stats)
#define sl_pepc_receive_stats           SL_PEPC_FUNC(sl_pepc_receive_stats)
#define sl_pepc_sort_keys               SL_PEPC_FUNC(sl_pepc_sort_keys)
#define sl_pepc_sort_keys_              SL_PEPC_FUNC(sl_pepc_sort_keys_)
#define sl_pepc_sort_parts              SL_PEPC_FUNC(sl_pepc_sort_parts)
#define sl_pepc_sort_parts_             SL_PEPC_FUNC(sl_pepc_sort_parts_)


#endif /* __RENAME_PEPC_SORT_H__ */
