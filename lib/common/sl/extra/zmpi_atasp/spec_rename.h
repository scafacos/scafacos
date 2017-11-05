/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Michael Hofmann
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


#ifndef __SPEC_RENAME_H__
#define __SPEC_RENAME_H__


#define SP_CONCAT(_a_, _b_)           SP_CONCAT_(_a_, _b_)
#define SP_CONCAT_(_a_, _b_)          _a_##_b_

#define SP_CONCONCAT(_a_, _b_, _c_)   SP_CONCONCAT_(_a_, _b_, _c_)
#define SP_CONCONCAT_(_a_, _b_, _c_)  _a_##_b_##_c_

#ifdef SP_PREFIX
# define SP_VAR(_v_)   SP_CONCAT(SP_PREFIX, _v_)
# define SP_FUNC(_f_)  SP_CONCAT(SP_PREFIX, _f_)
#else
# define SP_VAR(_v_)   _v_
# define SP_FUNC(_f_)  _f_
#endif


/* spec_alltoallv.c */
#define spec_alltoallv_db  SP_FUNC(spec_alltoallv_db)
#define spec_alltoallv_ip  SP_FUNC(spec_alltoallv_ip)

/* spec_alltoallw.c */
#define spec_alltoallw_db  SP_FUNC(spec_alltoallw_db)

/* spec_common.c */
#define spec_check_tproc_support  SP_FUNC(spec_check_tproc_support)
#define spec_check_buffer_size  SP_FUNC(spec_check_buffer_size)
#define spec_tproc_setup  SP_FUNC(spec_tproc_setup)
#define spec_tproc_release  SP_FUNC(spec_tproc_release)
#define spec_make_counts  SP_FUNC(spec_make_counts)
#define spec_redistribute_counts_type  SP_VAR(spec_redistribute_counts_type)
#define spec_redistribute_counts_proclists_type  SP_VAR(spec_redistribute_counts_proclists_type)
#define spec_redistribute_counts  SP_FUNC(spec_redistribute_counts)
#define spec_reduce_scatter_counts_type  SP_VAR(spec_reduce_scatter_counts_type)
#define spec_reduce_scatter_counts_proclists_type  SP_VAR(spec_reduce_scatter_counts_proclists_type)
#define spec_reduce_scatter_counts  SP_FUNC(spec_reduce_scatter_counts)
#define spec_prefix_counts_type  SP_VAR(spec_prefix_counts_type)
#define spec_prefix_counts_proclists_type  SP_VAR(spec_prefix_counts_proclists_type)
#define spec_prefix_counts  SP_FUNC(spec_prefix_counts)

/* spec_core.c */
#define spec_timing  SP_VAR(spec_timing)
#define spec_tproc_create  SP_FUNC(spec_tproc_create)
#define spec_tproc_destroy  SP_FUNC(spec_tproc_destroy)
#define spec_tproc_duplicate  SP_FUNC(spec_tproc_duplicate)
#define spec_tproc_set_tproc  SP_FUNC(spec_tproc_set_tproc)
#define spec_tproc_set_ext_tproc  SP_FUNC(spec_tproc_set_ext_tproc)
#define spec_tproc_set_tproc_mod  SP_FUNC(spec_tproc_set_tproc_mod)
#define spec_tproc_set_ext_tproc_mod  SP_FUNC(spec_tproc_set_ext_tproc_mod)
#define spec_tproc_set_tprocs  SP_FUNC(spec_tproc_set_tprocs)
#define spec_tproc_set_ext_tprocs  SP_FUNC(spec_tproc_set_ext_tprocs)
#define spec_tproc_set_tprocs_mod  SP_FUNC(spec_tproc_set_tprocs_mod)
#define spec_tproc_set_ext_tprocs_mod  SP_FUNC(spec_tproc_set_ext_tprocs_mod)
#define spec_tproc_set_reset  SP_FUNC(spec_tproc_set_reset)
#define spec_make_recv_proclist  SP_FUNC(spec_make_recv_proclist)
#define spec_tproc_set_proclists  SP_FUNC(spec_tproc_set_proclists)
#define spec_print  SP_FUNC(spec_print)

/* spec_put.c */
#define spec_put_db  SP_FUNC(spec_put_db)
#define spec_put_2phases_db  SP_FUNC(spec_put_2phases_db)

/* spec_sendrecv.c */
#define spec_sendrecv_aux  SP_VAR(spec_sendrecv_aux)
#define spec_sendrecv_aux_size  SP_VAR(spec_sendrecv_aux_size)
#define spec_sendrecv_send_requests  SP_VAR(spec_sendrecv_send_requests)
#define spec_sendrecv_receive_requests  SP_VAR(spec_sendrecv_receive_requests)
#define spec_sendrecv_db  SP_FUNC(spec_sendrecv_db)


#endif /* __SPEC_RENAME_H__ */
