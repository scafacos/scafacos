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


#ifndef __ZMPI_TOOLS_RENAME_H__
#define __ZMPI_TOOLS_RENAME_H__


#ifndef ZMPI_RENAME

#define ZMPI_RENAME

#define ZMPI_CONCAT(_a_, _b_)           ZMPI_CONCAT_(_a_, _b_)
#define ZMPI_CONCAT_(_a_, _b_)          _a_##_b_

#define ZMPI_CONCONCAT(_a_, _b_, _c_)   ZMPI_CONCONCAT_(_a_, _b_, _c_)
#define ZMPI_CONCONCAT_(_a_, _b_, _c_)  _a_##_b_##_c_

#ifdef ZMPI_PREFIX
# define ZMPI_VAR(_v_)   ZMPI_CONCAT(ZMPI_PREFIX, _v_)
# define ZMPI_FUNC(_f_)  ZMPI_CONCAT(ZMPI_PREFIX, _f_)
#else
# define ZMPI_VAR(_v_)   _v_
# define ZMPI_FUNC(_f_)  _f_
#endif

#endif /* ZMPI_RENAME */


/* zmpi_alltoall_2step.c */
#define ZMPI_Alltoall_2step_int  ZMPI_FUNC(ZMPI_Alltoall_2step_int)

/* zmpi_alltoall_int.c */
#define ZMPI_Alltoall_int_alltoall  ZMPI_FUNC(ZMPI_Alltoall_int_alltoall)
#define ZMPI_Alltoall_int_2step  ZMPI_FUNC(ZMPI_Alltoall_int_2step)
#define ZMPI_Alltoall_int_put  ZMPI_FUNC(ZMPI_Alltoall_int_put)
#define ZMPI_Alltoall_int_put_alloc  ZMPI_FUNC(ZMPI_Alltoall_int_put_alloc)
#define ZMPI_Alltoall_int_put_2phases  ZMPI_FUNC(ZMPI_Alltoall_int_put_2phases)
#define ZMPI_Alltoall_int_put_2phases_alloc  ZMPI_FUNC(ZMPI_Alltoall_int_put_2phases_alloc)
#define ZMPI_Alltoall_int_put_3phases  ZMPI_FUNC(ZMPI_Alltoall_int_put_3phases)
#define ZMPI_Alltoall_int_put_3phases_alloc  ZMPI_FUNC(ZMPI_Alltoall_int_put_3phases_alloc)
#define ZMPI_Alltoall_int_proclists_isendirecv  ZMPI_FUNC(ZMPI_Alltoall_int_proclists_isendirecv)
#define ZMPI_Alltoall_int_proclists_alltoallv  ZMPI_FUNC(ZMPI_Alltoall_int_proclists_alltoallv)
#define ZMPI_Alltoall_int_proclists_put  ZMPI_FUNC(ZMPI_Alltoall_int_proclists_put)
#define ZMPI_Alltoall_int_proclists_put_alloc  ZMPI_FUNC(ZMPI_Alltoall_int_proclists_put_alloc)
#define ZMPI_Alltoall_int_proclists_put_2phases  ZMPI_FUNC(ZMPI_Alltoall_int_proclists_put_2phases)
#define ZMPI_Alltoall_int_proclists_put_2phases_alloc  ZMPI_FUNC(ZMPI_Alltoall_int_proclists_put_2phases_alloc)

/* zmpi_alltoallx_proclists.c */
#define ZMPI_Alltoall_proclists_isendirecv  ZMPI_FUNC(ZMPI_Alltoall_proclists_isendirecv)
#define ZMPI_Alltoall_proclists  ZMPI_FUNC(ZMPI_Alltoall_proclists)
#define ZMPI_Alltoallv_proclists_isendirecv  ZMPI_FUNC(ZMPI_Alltoallv_proclists_isendirecv)
#define ZMPI_Alltoallv_proclists  ZMPI_FUNC(ZMPI_Alltoallv_proclists)
#define ZMPI_Alltoallw_proclists_isendirecv  ZMPI_FUNC(ZMPI_Alltoallw_proclists_isendirecv)
#define ZMPI_Alltoallw_proclists  ZMPI_FUNC(ZMPI_Alltoallw_proclists)

/* zmpi_reduce_scatter_block.c */
#define ZMPI_Reduce_scatter_block  ZMPI_FUNC(ZMPI_Reduce_scatter_block)

/* zmpi_reduce_scatter_block_intsum.c */
#define ZMPI_Reduce_scatter_block_intsum_accumulate  ZMPI_FUNC(ZMPI_Reduce_scatter_block_intsum_accumulate)
#define ZMPI_Reduce_scatter_block_intsum_proclists_isendirecv  ZMPI_FUNC(ZMPI_Reduce_scatter_block_intsum_proclists_isendirecv)
#define ZMPI_Reduce_scatter_block_intsum_proclists_alltoallv  ZMPI_FUNC(ZMPI_Reduce_scatter_block_intsum_proclists_alltoallv)
#define ZMPI_Reduce_scatter_block_intsum_proclists_accumulate  ZMPI_FUNC(ZMPI_Reduce_scatter_block_intsum_proclists_accumulate)


#endif /* __ZMPI_TOOLS_RENAME_H__ */
