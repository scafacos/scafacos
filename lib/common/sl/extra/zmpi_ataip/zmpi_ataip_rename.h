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


#ifndef __ZMPI_ATAIP_RENAME_H__
#define __ZMPI_ATAIP_RENAME_H__


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


/* zmpi_ataip.c */
#define ZMPI_Alltoallv_inplace_sym_type  ZMPI_VAR(ZMPI_Alltoallv_inplace_sym_type)
#define ZMPI_Alltoallv_inplace_aux_type  ZMPI_VAR(ZMPI_Alltoallv_inplace_aux_type)
#define ZMPI_Alltoallv_inplace_aux  ZMPI_VAR(ZMPI_Alltoallv_inplace_aux)
#define ZMPI_Alltoallv_inplace_aux_size  ZMPI_VAR(ZMPI_Alltoallv_inplace_aux_size)
#define ZMPI_Alltoallv_inplace_aux_static_blocks  ZMPI_VAR(ZMPI_Alltoallv_inplace_aux_static_blocks)
#define ZMPI_Alltoallv_inplace_sync_on_init  ZMPI_VAR(ZMPI_Alltoallv_inplace_sync_on_init)
#define ZMPI_Alltoallv_inplace_sync_on_exit  ZMPI_VAR(ZMPI_Alltoallv_inplace_sync_on_exit)
#define ZMPI_Alltoallv_inplace_sync_run  ZMPI_VAR(ZMPI_Alltoallv_inplace_sync_run)
#define ZMPI_Alltoall_inplace_symmetric_sym_type  ZMPI_VAR(ZMPI_Alltoall_inplace_symmetric_sym_type)
#define ZMPI_Alltoallv_inplace_symmetric_sym_type  ZMPI_VAR(ZMPI_Alltoallv_inplace_symmetric_sym_type)
#define ZMPI_Alltoallv_inplace_symmetric_aux  ZMPI_VAR(ZMPI_Alltoallv_inplace_symmetric_aux)
#define ZMPI_Alltoallv_inplace_symmetric_aux_size  ZMPI_VAR(ZMPI_Alltoallv_inplace_symmetric_aux_size)
#define ZMPI_Alltoallv_inplace_times  ZMPI_VAR(ZMPI_Alltoallv_inplace_times)
#define ZMPI_Alltoallv_inplace  ZMPI_FUNC(ZMPI_Alltoallv_inplace)
#define ZMPI_Neighbor_alltoallv_inplace  ZMPI_FUNC(ZMPI_Neighbor_alltoallv_inplace)
#define ZMPI_Alltoall_inplace_symmetric  ZMPI_FUNC(ZMPI_Alltoall_inplace_symmetric)
#define ZMPI_Neighbor_alltoall_inplace_symmetric  ZMPI_FUNC(ZMPI_Neighbor_alltoall_inplace_symmetric)
#define ZMPI_Alltoallv_inplace_symmetric  ZMPI_FUNC(ZMPI_Alltoallv_inplace_symmetric)
#define ZMPI_Neighbor_alltoallv_inplace_symmetric  ZMPI_FUNC(ZMPI_Neighbor_alltoallv_inplace_symmetric)


#endif /* __ZMPI_ATAIP_RENAME_H__ */
