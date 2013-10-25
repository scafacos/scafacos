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


#ifndef __SL_CONFIG_H__
#define __SL_CONFIG_H__


#include "config_fmm_sort.h"


/* standard (SL) integer data type */
#define sl_int_type_c             SL_INTEGER_C
#define sl_int_type_mpi           SL_INTEGER_MPI
#define sl_int_size_mpi           1
#define sl_int_type_fmt           SL_INTEGER_FMT


/* key section */
#define sl_key_type_c             INTEGER_C
#define sl_key_type_mpi           INTEGER_MPI
#define sl_key_size_mpi           1

#define sl_key_integer
#define sl_key_type_fmt           INTEGER_FMT

#define sl_key_purify(_k_)        (_k_ & SL_FMM_CONFIG_VAR(fmm_front_key_mask))


/* data section */
#define SL_DATA0                  /* xyz */
#define sl_data0_type_c           REAL_C
#define sl_data0_size_c           3
#define sl_data0_type_mpi         REAL_MPI
#define sl_data0_size_mpi         3

#define SL_DATA1                  /* q */
#define sl_data1_type_c           REAL_C
#define sl_data1_size_c           1
#define sl_data1_type_mpi         REAL_MPI
#define sl_data1_size_mpi         1

#undef SL_DATA2                   /* scr */
#define sl_data2_type_c           INTEGER_C
#define sl_data2_size_c           1
#define sl_data2_type_mpi         INTEGER_MPI
#define sl_data2_size_mpi         1

#undef SL_DATA3                   /* addr */
#define sl_data3_type_c           INTEGER_C
#define sl_data3_size_c           1
#define sl_data3_type_mpi         INTEGER_MPI
#define sl_data3_size_mpi         1

#endif /* __SL_CONFIG_H__ */
