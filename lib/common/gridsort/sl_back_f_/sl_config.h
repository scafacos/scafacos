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


#define sl_key_type_c              long long
#define sl_key_type_mpi            MPI_LONG_LONG
#define sl_key_size_mpi            1
#define sl_key_type_fmt            "lld"

#define sl_key_integer

#define SL_DATA0
#define sl_data0_type_c            fcs_float
#define sl_data0_size_c            3
#define sl_data0_type_mpi          FCS_MPI_FLOAT
#define sl_data0_size_mpi          3

/*#define SL_DATA1*/
#define sl_data1_type_c            fcs_float
#define sl_data1_size_c            1
#define sl_data1_type_mpi          FCS_MPI_FLOAT
#define sl_data1_size_mpi          1

#endif /* __SL_CONFIG_H__ */
