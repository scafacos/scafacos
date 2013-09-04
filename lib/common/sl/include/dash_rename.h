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


#ifndef __DASH_RENAME_H__
#define __DASH_RENAME_H__


#define DS_CONCAT(_a_, _b_)           DS_CONCAT_(_a_, _b_)
#define DS_CONCAT_(_a_, _b_)          _a_##_b_

#define DS_CONCONCAT(_a_, _b_, _c_)   DS_CONCONCAT_(_a_, _b_, _c_)
#define DS_CONCONCAT_(_a_, _b_, _c_)  _a_##_b_##_c_

#ifdef DS_PREFIX
# define DS_VAR(_v_)   DS_CONCAT(DS_PREFIX, _v_)
# define DS_FUNC(_f_)  DS_CONCAT(DS_PREFIX, _f_)
#else
# define DS_VAR(_v_)   _v_
# define DS_FUNC(_f_)  _f_
#endif


/* dash_aux_heap.c */
#define ds_aux_heap_create  DS_FUNC(ds_aux_heap_create)
#define ds_aux_heap_destroy  DS_FUNC(ds_aux_heap_destroy)

/* dash_aux_static.c */
#define ds_aux_static_create  DS_FUNC(ds_aux_static_create)
#define ds_aux_static_destroy  DS_FUNC(ds_aux_static_destroy)

/* dash_common.c */
#define ds_sort_dsints  DS_FUNC(ds_sort_dsints)

/* dash_core.c */
#define ds_core_run_sync  DS_VAR(ds_core_run_sync)
#define ds_times  DS_VAR(ds_times)
#define ds_create  DS_FUNC(ds_create)
#define ds_destroy  DS_FUNC(ds_destroy)
#define ds_run  DS_FUNC(ds_run)

/* dash_exec.c */
#define ds_exec_create  DS_FUNC(ds_exec_create)
#define ds_exec_destroy  DS_FUNC(ds_exec_destroy)
#define ds_exec_add_address  DS_FUNC(ds_exec_add_address)
#define ds_exec_make  DS_FUNC(ds_exec_make)

/* dash_exec_mpi.c */
#define ds_exec_mpi_create  DS_FUNC(ds_exec_mpi_create)
#define ds_exec_mpi_destroy  DS_FUNC(ds_exec_mpi_destroy)
#define ds_exec_mpi_add_address  DS_FUNC(ds_exec_mpi_add_address)
#define ds_exec_mpi_add_type  DS_FUNC(ds_exec_mpi_add_type)
#define ds_exec_mpi_sizefor  DS_FUNC(ds_exec_mpi_sizefor)
#define ds_exec_mpi_extent  DS_FUNC(ds_exec_mpi_extent)

/* dash_sched_a2a_sym.c */
#define ds_sched_a2a_sym_create  DS_FUNC(ds_sched_a2a_sym_create)
#define ds_sched_a2a_sym_destroy  DS_FUNC(ds_sched_a2a_sym_destroy)

/* dash_sched_a2av_aux.c */
#define ds_aux_create  DS_FUNC(ds_aux_create)
#define ds_aux_destroy  DS_FUNC(ds_aux_destroy)
#define ds_aux_pre_acquire  DS_FUNC(ds_aux_pre_acquire)
#define ds_aux_acquire  DS_FUNC(ds_aux_acquire)
#define ds_aux_post_acquire  DS_FUNC(ds_aux_post_acquire)
#define ds_aux_get_count  DS_FUNC(ds_aux_get_count)
#define ds_aux_get_displ  DS_FUNC(ds_aux_get_displ)
#define ds_aux_vacate  DS_FUNC(ds_aux_vacate)
#define ds_aux_accept_recv  DS_FUNC(ds_aux_accept_recv)

/* dash_sched_a2av.c */
#define ds_sched_a2av_aux_blocks  DS_VAR(ds_sched_a2av_aux_blocks)
#define ds_sched_a2av_create  DS_FUNC(ds_sched_a2av_create)
#define ds_sched_a2av_destroy  DS_FUNC(ds_sched_a2av_destroy)
#define ds_sched_a2av_set_aux  DS_FUNC(ds_sched_a2av_set_aux)
#define ds_sched_a2av_skip_sym  DS_FUNC(ds_sched_a2av_skip_sym)

/* dash_sched_a2av_sym.c */
#define ds_sched_a2av_sym_create  DS_FUNC(ds_sched_a2av_sym_create)
#define ds_sched_a2av_sym_destroy  DS_FUNC(ds_sched_a2av_sym_destroy)

/* dash_sched.c */
#define ds_sched_create  DS_FUNC(ds_sched_create)
#define ds_sched_destroy  DS_FUNC(ds_sched_destroy)
#define ds_sched_add_buffer  DS_FUNC(ds_sched_add_buffer)
#define ds_sched_set_send  DS_FUNC(ds_sched_set_send)
#define ds_sched_set_recv  DS_FUNC(ds_sched_set_recv)
#define ds_sched_set_aux  DS_FUNC(ds_sched_set_aux)


#endif /* __DASH_RENAME_H__ */
