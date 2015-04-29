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


#ifndef __SL_DEPRECATED_H__
#define __SL_DEPRECATED_H__


/* sl_macro MPI_BINNING_REDUCEBCAST_THRESHOLD */
#ifdef MB_REDUCEBCAST_THRESHOLD
# define MPI_BINNING_REDUCEBCAST_THRESHOLD  MB_REDUCEBCAST_THRESHOLD
#endif


#define key_mpi_datatype     sl_default_context.mpi.key_datatype
#define int_mpi_datatype     sl_default_context.mpi.int_datatype
#define pkey_mpi_datatype    sl_default_context.mpi.pkey_datatype
#define pwkey_mpi_datatype   sl_default_context.mpi.pwkey_datatype
#define index_mpi_datatype   sl_default_context.mpi.index_datatype
#define weight_mpi_datatype  sl_default_context.mpi.weight_datatype
#define data_mpi_datatype    sl_default_context.mpi.data_datatype

#undef sl_mpi_rank
#define sl_mpi_rank   SL_PROC_RANK

#endif /* __SL_DEPRECATED_H__ */
